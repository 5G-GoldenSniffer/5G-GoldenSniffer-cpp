#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include "NRscrambling.h"
#include "pbch.h"
#include "ofdm.h"

namespace {

const uint32_t PBCH_A = 32;          // payload with the timing bits, 38.212 7.1.1
const uint32_t PBCH_K = PBCH_A + 24; // with the CRC
const uint32_t PBCH_E = 864;         // 432 QPSK REs
const uint32_t PBCH_M = PBCH_A - 3;  // scrambled payload bits for L_max 4 or 8 (38.212 7.1.2)
const int PBCH_list_size = 8;
const int PBCH_hypotheses = 2;       // DM-RS hypotheses decoded, from the best match down

// 38.212 Table 7.1.1-1
const int G[32] = {16, 23, 18, 17, 8, 30, 10, 6, 24, 7, 0, 5, 3, 2, 1, 4,
                   9, 11, 12, 13, 14, 15, 19, 20, 21, 22, 25, 26, 27, 28, 29, 31};

// position in the interleaved payload a' of every bit of the payload abar (38.212 7.1.1): abar[0] is the
// BCCH-BCH-Message choice, abar[1..23] the MIB (SFN MSBs in abar[1..6]), abar[24..27] the 4 LSBs of the
// SFN, abar[28] the half frame bit, abar[29..31] the MSB of k_SSB and two reserved bits (L_max 4 or 8)
struct PayloadMap {
    int pos[PBCH_A];
    PayloadMap()
    {
        int j_SFN = 0, j_HRF = 10, j_SSB = 11, j_other = 14;
        for(uint32_t i = 0; i < PBCH_A; i++){
            if((i >= 1 && i <= 6) || (i >= 24 && i <= 27))
                pos[i] = G[j_SFN++];
            else if(i == 28)
                pos[i] = G[j_HRF];
            else if(i >= 29)
                pos[i] = G[j_SSB++];
            else
                pos[i] = G[j_other++];
        }
    }
};
const PayloadMap payload_map;

// payload bits not scrambled (38.212 7.1.2): half frame bit, 3rd and 2nd LSBs of the SFN
bool unscrambled(int i)
{
    return i == G[10] || i == G[7] || i == G[8];
}

unsigned read_bits(const uint8_t *b, int n)
{
    unsigned v = 0;
    for(int i = 0; i < n; i++)
        v = (v<<1) | b[i];
    return v;
}

void write_bits(uint8_t *b, unsigned v, int n)
{
    for(int i = 0; i < n; i++)
        b[i] = (v>>(n-1-i)) & 1;
}

// REs of the PBCH DM-RS and data in the SS block, as l*SSB_carriers+k in mapping order (38.211 Table 7.4.3.1-1:
// first k, then l), and the DM-RS sequence (38.211 7.4.1.4.1) for i_SSB_bar
void PBCH_REs(uint32_t PCI, std::vector<uint16_t> &dmrs, std::vector<uint16_t> &data)
{
    const uint32_t v = PCI%4;
    dmrs.clear();
    data.clear();
    for(uint32_t l = 1; l < SSB_symbols; l++)
        for(uint32_t k = 0; k < SSB_carriers; k++){
            if(l == 2 && k >= 48 && k < 192) // SSS
                continue;
            (k%4 == v ? dmrs : data).push_back(l*SSB_carriers + k);
        }
}

void PBCH_DMRS(std::vector<cfloat> &r, uint32_t PCI, int i_bar)
{
    const size_t N_DMRS = 144;
    std::vector<uint8_t> c(2*N_DMRS);
    scrambling(c.data(), 2*N_DMRS, PBCH_DMRS_c_init(PCI, i_bar));
    r.resize(N_DMRS);
    for(size_t m = 0; m < N_DMRS; m++)
        r[m] = cfloat(1.0f-2.0f*c[2*m], 1.0f-2.0f*c[2*m+1])*(float)M_SQRT1_2;
}

// runs of DM-RS spaced by 4 subcarriers in the same symbol: [begin,end) in dmrs
std::vector<std::pair<size_t,size_t>> DMRS_runs(const std::vector<uint16_t> &dmrs)
{
    std::vector<std::pair<size_t,size_t>> runs;
    size_t m = 0;
    while(m < dmrs.size()){
        size_t e = m+1;
        while(e < dmrs.size() && dmrs[e] == dmrs[e-1]+4 && dmrs[e]/SSB_carriers == dmrs[m]/SSB_carriers)
            e++;
        runs.push_back({m, e});
        m = e;
    }
    return runs;
}

// how well the DM-RS hypothesis r matches: correlation between adjacent DM-RS after removing r, robust
// to the phase slope of a timing error
float DMRS_match(const cfloat *ssb, const std::vector<uint16_t> &dmrs, const std::vector<cfloat> &r,
                 const std::vector<std::pair<size_t,size_t>> &runs)
{
    cfloat acc(0.0f, 0.0f);
    float norm = 0.0f;
    for(const auto &run : runs)
        for(size_t m = run.first; m+1 < run.second; m++){
            cfloat z0 = ssb[dmrs[m]]*std::conj(r[m]);
            cfloat z1 = ssb[dmrs[m+1]]*std::conj(r[m+1]);
            acc += z1*std::conj(z0);
            norm += std::abs(z0)*std::abs(z1);
        }
    return norm > 0.0f ? std::abs(acc)/norm : 0.0f;
}

// channel at every subcarrier of the PBCH: least squares at the DM-RS, averaged over 3 adjacent DM-RS and
// linearly interpolated within each run of DM-RS (held constant beyond its ends)
void PBCH_channel(const cfloat *ssb, const std::vector<uint16_t> &dmrs, const std::vector<cfloat> &r,
                  const std::vector<std::pair<size_t,size_t>> &runs, uint32_t PCI, cfloat *H)
{
    const int v = PCI%4;
    std::vector<cfloat> h, hs;
    for(const auto &run : runs){
        const size_t np = run.second - run.first;
        h.resize(np);
        hs.resize(np);
        for(size_t i = 0; i < np; i++)
            h[i] = ssb[dmrs[run.first+i]]*std::conj(r[run.first+i]); // |r| == 1
        smooth3(h.data(), hs.data(), np);
        const int base = dmrs[run.first] - dmrs[run.first]%SSB_carriers; // l*SSB_carriers
        const int k_first = dmrs[run.first]%SSB_carriers;
        const int k_last = dmrs[run.second-1]%SSB_carriers;
        // the run covers the PBCH subcarriers from k_first-v to k_last+3-v
        for(int k = k_first-v; k <= k_last+3-v; k++){
            cfloat Hk;
            if(k <= k_first)
                Hk = hs[0];
            else if(k >= k_last)
                Hk = hs[np-1];
            else{
                int i = (k-k_first)/4;
                float a = ((k-k_first)%4)/4.0f;
                Hk = hs[i]*(1.0f-a) + hs[i+1]*a;
            }
            H[base+k] = Hk;
        }
    }
}

} // namespace

uint32_t PBCH_DMRS_c_init(uint32_t PCI, uint32_t i_bar)
{
    return ((i_bar+1)*(PCI/4+1) << 11) + ((i_bar+1) << 6) + PCI%4;
}

int SSB_first_symbol(int mu, int SSB_index)
{
    // Case A: {2,8}+14n, n=0,1; Case C: {2,8}+14n, n=0..3
    return 2 + 6*(SSB_index%2) + 14*(SSB_index/2);
}

PBCHDecoder::PBCHDecoder(int mu) : L_max(mu == 0 ? 4 : 8), pc(PBCH_K, PBCH_E, 9, true) {}

const std::vector<uint8_t> &PBCHDecoder::scrambling_PCI(uint32_t PCI)
{
    // 38.211 7.3.3.1 (v*M_bit offsets, v < 8) and 38.212 7.1.2 (v*M offsets, v < 4), both with c_init = PCI
    if(PCI != c_PCI){
        c.resize(8*PBCH_E);
        scrambling(c.data(), c.size(), PCI);
        c_PCI = PCI;
    }
    return c;
}

bool PBCHDecoder::decode(const cfloat *ssb, uint32_t PCI, MIB &mib)
{
    std::vector<uint16_t> dmrs, data;
    PBCH_REs(PCI, dmrs, data);
    const auto runs = DMRS_runs(dmrs);

    // DM-RS hypotheses i_SSB_bar (SS block index, plus the half frame for L_max=4), best match first
    std::vector<std::vector<cfloat>> r(8);
    std::vector<std::pair<float,int>> match(8);
    for(int i_bar = 0; i_bar < 8; i_bar++){
        PBCH_DMRS(r[i_bar], PCI, i_bar);
        match[i_bar] = {DMRS_match(ssb, dmrs, r[i_bar], runs), i_bar};
    }
    std::sort(match.begin(), match.end(), [](const auto &a, const auto &b){ return a.first > b.first; });

    const std::vector<uint8_t> &c = scrambling_PCI(PCI);
    std::vector<cfloat> H(SSB_symbols*SSB_carriers);
    std::vector<float> llr(PBCH_E);
    std::vector<std::vector<uint8_t>> cand;
    for(int hyp = 0; hyp < PBCH_hypotheses; hyp++){
        const int i_bar = match[hyp].second;
        PBCH_channel(ssb, dmrs, r[i_bar], runs, PCI, H.data());
        // matched filter: the LLR scale doesn't matter to the min-sum decoder
        for(size_t j = 0; j < data.size(); j++){
            cfloat z = ssb[data[j]]*std::conj(H[data[j]]);
            llr[2*j] = z.real();
            llr[2*j+1] = z.imag();
        }
        // descrambling, 38.211 7.3.3.1
        const uint32_t v = i_bar & (L_max-1);
        for(uint32_t i = 0; i < PBCH_E; i++)
            if(c[i + v*PBCH_E])
                llr[i] = -llr[i];

        pc.decode_list(llr.data(), PBCH_list_size, cand);
        for(const auto &a2 : cand){
            if(crc24c(a2.data(), PBCH_A, 0) != read_bits(a2.data()+PBCH_A, 24))
                continue;
            // payload descrambling (38.212 7.1.2): v from the 3rd and 2nd LSBs of the SFN, which aren't scrambled
            const uint32_t v_SFN = 2*a2[G[7]] + a2[G[8]];
            uint8_t a1[PBCH_A];
            for(uint32_t i = 0, j = 0; i < PBCH_A; i++)
                a1[i] = unscrambled(i) ? a2[i] : a2[i] ^ c[(j++) + v_SFN*PBCH_M];
            uint8_t abar[PBCH_A];
            for(uint32_t i = 0; i < PBCH_A; i++)
                abar[i] = a1[payload_map.pos[i]];
            if(abar[0] != 0) // messageClassExtension
                continue;
            const int HRF = abar[28];
            if(L_max == 4 && HRF != (i_bar>>2)) // the half frame is also in the DM-RS
                continue;
            mib.SFN = (read_bits(abar+1, 6)<<4) | read_bits(abar+24, 4);
            mib.HRF = HRF;
            mib.SSB_index = i_bar & (L_max-1);
            mib.scs_common = abar[7];
            mib.k_SSB = (abar[29]<<4) | read_bits(abar+8, 4);
            mib.dmrs_typeA_pos = abar[12] ? 3 : 2;
            mib.coreset0 = read_bits(abar+13, 4);
            mib.ss0 = read_bits(abar+17, 4);
            mib.cell_barred = abar[21] == 0;      // ENUMERATED {barred, notBarred}
            mib.intra_freq_resel = abar[22] == 0; // ENUMERATED {allowed, notAllowed}
            return true;
        }
    }
    return false;
}

void PBCHDecoder::encode(cfloat *ssb, const MIB &mib, uint32_t PCI)
{
    uint8_t abar[PBCH_A] = {0};
    write_bits(abar+1, mib.SFN>>4, 6);
    abar[7] = mib.scs_common;
    write_bits(abar+8, mib.k_SSB & 0x0f, 4);
    abar[12] = mib.dmrs_typeA_pos == 3;
    write_bits(abar+13, mib.coreset0, 4);
    write_bits(abar+17, mib.ss0, 4);
    abar[21] = !mib.cell_barred;
    abar[22] = !mib.intra_freq_resel;
    write_bits(abar+24, mib.SFN & 0x0f, 4);
    abar[28] = mib.HRF;
    abar[29] = (mib.k_SSB>>4) & 1;

    const std::vector<uint8_t> &c = scrambling_PCI(PCI);
    uint8_t a1[PBCH_A], a2[PBCH_K];
    for(uint32_t i = 0; i < PBCH_A; i++)
        a1[payload_map.pos[i]] = abar[i];
    const uint32_t v_SFN = 2*((mib.SFN>>2) & 1) + ((mib.SFN>>1) & 1);
    for(uint32_t i = 0, j = 0; i < PBCH_A; i++)
        a2[i] = unscrambled(i) ? a1[i] : a1[i] ^ c[(j++) + v_SFN*PBCH_M];
    write_bits(a2+PBCH_A, crc24c(a2, PBCH_A, 0), 24);

    std::vector<uint8_t> e(PBCH_E);
    pc.encode(a2, e.data());
    const int i_bar = L_max == 4 ? mib.SSB_index + 4*mib.HRF : mib.SSB_index;
    const uint32_t v = i_bar & (L_max-1);
    for(uint32_t i = 0; i < PBCH_E; i++)
        e[i] ^= c[i + v*PBCH_E];

    std::vector<uint16_t> dmrs, data;
    PBCH_REs(PCI, dmrs, data);
    std::vector<cfloat> r;
    PBCH_DMRS(r, PCI, i_bar);
    for(size_t m = 0; m < dmrs.size(); m++)
        ssb[dmrs[m]] = r[m];
    for(size_t j = 0; j < data.size(); j++)
        ssb[data[j]] = cfloat(1.0f-2.0f*e[2*j], 1.0f-2.0f*e[2*j+1])*(float)M_SQRT1_2;
}

bool PBCHDecoder::selftest(bool verbose)
{
    std::mt19937 rng(12345);
    std::normal_distribution<float> gauss(0.0f, 1.0f);
    std::uniform_real_distribution<float> unif(0.0f, 1.0f);
    const int trials = 20;
    const float SNR_dB = 0.0f;
    const float sigma = std::pow(10.0f, -SNR_dB/20.0f)*(float)M_SQRT1_2; // per real dimension
    std::vector<cfloat> ssb(SSB_symbols*SSB_carriers);
    bool is_ok = true;
    for(int t = 0; t <= trials; t++){
        const bool noise_only = t == trials; // must not decode
        MIB tx;
        const uint32_t PCI = rng()%1008;
        tx.SFN = rng()%1024;
        tx.HRF = rng()%2;
        tx.SSB_index = rng()%L_max;
        tx.scs_common = rng()%2;
        tx.k_SSB = rng()%32;
        tx.dmrs_typeA_pos = 2 + rng()%2;
        tx.coreset0 = rng()%16;
        tx.ss0 = rng()%16;
        tx.cell_barred = rng()%2;
        tx.intra_freq_resel = rng()%2;
        std::fill(ssb.begin(), ssb.end(), cfloat(0.0f, 0.0f));
        if(!noise_only)
            encode(ssb.data(), tx, PCI);
        // channel: gain, phase and the phase slope of up to one sample of timing error at 128 samples per symbol
        const float gain = 0.5f + unif(rng), phi = 2.0f*(float)M_PI*unif(rng), tau = 2.0f*unif(rng)-1.0f;
        for(size_t l = 0; l < SSB_symbols; l++)
            for(size_t k = 0; k < SSB_carriers; k++){
                cfloat h = std::polar(gain, phi - 2.0f*(float)M_PI*tau*(float)k/128.0f);
                ssb[l*SSB_carriers+k] = ssb[l*SSB_carriers+k]*h + gain*sigma*cfloat(gauss(rng), gauss(rng));
            }
        MIB rx;
        bool decoded = decode(ssb.data(), PCI, rx);
        bool ok;
        if(noise_only)
            ok = !decoded;
        else
            ok = decoded && rx.SFN == tx.SFN && rx.HRF == tx.HRF && rx.SSB_index == tx.SSB_index
              && rx.scs_common == tx.scs_common && rx.k_SSB == tx.k_SSB && rx.dmrs_typeA_pos == tx.dmrs_typeA_pos
              && rx.coreset0 == tx.coreset0 && rx.ss0 == tx.ss0 && rx.cell_barred == tx.cell_barred
              && rx.intra_freq_resel == tx.intra_freq_resel;
        if(verbose)
            std::cout << "PBCH selftest #" << t << (noise_only ? " (noise only)" : "") << ": " << (ok ? "ok" : "FAILED") << std::endl;
        is_ok = is_ok && ok;
    }
    return is_ok;
}

bool CORESET0_config(const MIB &mib, int mu, int k120, size_t N, CORESET0 &cs)
{
    // 38.213 Tables 13-1 to 13-4: multiplexing pattern, N_RB, symbols, offset [RB]; N_RB == 0: reserved
    struct Entry { int8_t pattern, N_RB, duration, offset; };
    static const Entry T13_1[16] = { // {SS/PBCH block, PDCCH} SCS {15, 15} kHz
        {1,24,2,0}, {1,24,2,2}, {1,24,2,4}, {1,24,3,0}, {1,24,3,2}, {1,24,3,4}, {1,48,1,12}, {1,48,1,16},
        {1,48,2,12}, {1,48,2,16}, {1,48,3,12}, {1,48,3,16}, {1,96,1,38}, {1,96,2,38}, {1,96,3,38}, {0,0,0,0}};
    static const Entry T13_2[16] = { // {15, 30} kHz
        {1,24,2,5}, {1,24,2,6}, {1,24,2,7}, {1,24,2,8}, {1,24,3,5}, {1,24,3,6}, {1,24,3,7}, {1,24,3,8},
        {1,48,1,18}, {1,48,1,20}, {1,48,2,18}, {1,48,2,20}, {1,48,3,18}, {1,48,3,20}, {0,0,0,0}, {0,0,0,0}};
    static const Entry T13_3[16] = { // {30, 15} kHz
        {1,48,1,2}, {1,48,1,6}, {1,48,2,2}, {1,48,2,6}, {1,48,3,2}, {1,48,3,6}, {1,96,1,28}, {1,96,2,28},
        {1,96,3,28}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
    static const Entry T13_4[16] = { // {30, 30} kHz
        {1,24,2,0}, {1,24,2,1}, {1,24,2,2}, {1,24,2,3}, {1,24,2,4}, {1,24,3,0}, {1,24,3,1}, {1,24,3,2},
        {1,24,3,3}, {1,24,3,4}, {1,48,1,12}, {1,48,1,14}, {1,48,1,16}, {1,48,2,12}, {1,48,2,14}, {1,48,2,16}};
    const Entry *T = mu == 0 ? (mib.scs_common ? T13_2 : T13_1) : (mib.scs_common ? T13_4 : T13_3);
    const Entry &e = T[mib.coreset0 & 0x0f];

    cs = CORESET0{};
    cs.valid = mib.k_SSB <= 23 && e.N_RB > 0; // FR1: k_SSB > 23 means no CORESET#0 (38.213 13)
    if(!cs.valid)
        return false;
    cs.pattern = e.pattern;
    cs.N_RB = e.N_RB;
    cs.duration = e.duration;
    cs.offset_RB = e.offset;

    // k_SSB is in 15 kHz units in FR1; the common RBs have the CORESET subcarrier spacing
    cs.same_scs = mib.scs_common == mu;
    const double g0 = (double)N/2 + k120 - 120;        // subcarrier 0 of the SS block in the grid
    const double k_sc = mib.k_SSB*15.0/(15<<mu);        // in SS block subcarriers
    cs.CRB_grid_ofs = std::fmod(std::fmod(k_sc - g0, 12.0) + 12.0, 12.0);
    cs.k_start = g0 - std::fmod(k_sc, 12.0) - 12.0*cs.offset_RB;

    // 38.213 Table 13-11 (FR1, multiplexing pattern 1): O, M, first symbol (-1: 0 for even SS block
    // indices, the CORESET duration for odd ones)
    struct SS0 { int8_t O; float M; int8_t first_symbol; };
    static const SS0 T13_11[16] = {
        {0,1,0}, {0,0.5f,-1}, {2,1,0}, {2,0.5f,-1}, {5,1,0}, {5,0.5f,-1}, {7,1,0}, {7,0.5f,-1},
        {0,2,0}, {5,2,0}, {0,1,1}, {0,1,2}, {2,1,1}, {2,1,2}, {5,1,1}, {5,1,2}};
    const SS0 &s = T13_11[mib.ss0 & 0x0f];
    const int mu_c = mib.scs_common;
    const int slots = 10<<mu_c;
    const int slot = (s.O<<mu_c) + (int)std::floor(mib.SSB_index*s.M);
    cs.n0 = slot%slots;
    cs.frame_parity = (slot/slots)%2;
    cs.first_symbol = s.first_symbol >= 0 ? s.first_symbol : (mib.SSB_index%2 ? cs.duration : 0);
    return true;
}

std::string MIB_to_string(const MIB &mib, int mu, int k120, size_t N)
{
    std::ostringstream out;
    out << "MIB: SFN " << mib.SFN << ", SS block #" << mib.SSB_index << " in half frame " << mib.HRF
        << ", SCS common " << (15<<mib.scs_common) << " kHz, k_SSB " << mib.k_SSB
        << ", DM-RS type A position " << mib.dmrs_typeA_pos << ", cell " << (mib.cell_barred ? "barred" : "not barred")
        << ", pdcch-ConfigSIB1 " << mib.coreset0 << "/" << mib.ss0 << "\n";
    CORESET0 cs;
    if(!CORESET0_config(mib, mu, k120, N, cs)){
        out << "CORESET#0: none" << (mib.k_SSB > 23 ? " (k_SSB > 23)" : " (reserved configuration)");
        return out.str();
    }
    out << "CORESET#0: " << cs.N_RB << " RB x " << cs.duration << " symbol" << (cs.duration > 1 ? "s" : "")
        << ", " << cs.offset_RB << " RB below the SS block";
    if(cs.same_scs)
        out << ", subcarriers " << cs.k_start << ".." << cs.k_start + 12*cs.N_RB - 1 << " of the grid";
    else
        out << ", at a different SCS than the SS block (not in the grid)";
    out << "; Type0-PDCCH in slots " << cs.n0 << "," << cs.n0+1 << " of " << (cs.frame_parity ? "odd" : "even")
        << " frames from symbol " << cs.first_symbol;
    if(cs.same_scs && cs.CRB_grid_ofs != 0.0)
        out << "; the grid starts " << cs.CRB_grid_ofs << " subcarriers into a common RB";
    return out.str();
}
