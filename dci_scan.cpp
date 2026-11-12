#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <thread>
#include <mutex>
#include <optional>
#include <atomic>
#include <numeric>
#include <pthread.h>
#include <fftw3.h>
#include "nr_common.h"
#include "gf2solve.h"
#include "NRscrambling.h"
#include "DCIfactory.h"
#include "pipeline.h"
#include "ofdm.h"
#include "dci_scan.h"
#include "util.h"

// ---------- DCI size statistics ----------
// number of decoded DCIs per message length A (index A-A_min), shared by all workers
static std::atomic<uint32_t> DCI_size_hits[NA_max];
// when a counter reaches this value all counters are halved, so the ordering tracks changes in the cell
static const uint32_t DCI_size_hits_decay = 1024;

static void count_DCI_size(uint32_t A)
{
    // only the thread that brings the counter exactly to the threshold performs the decay
    if(DCI_size_hits[A-A_min].fetch_add(1, std::memory_order_relaxed)+1 != DCI_size_hits_decay)
        return;
    for(uint32_t i = 0; i < NA_max; i++){
        uint32_t v = DCI_size_hits[i].load(std::memory_order_relaxed);
        while(!DCI_size_hits[i].compare_exchange_weak(v, v/2, std::memory_order_relaxed));
    }
}

// DCI sizes tried by the blind search (--dci-sizes); empty: all sizes A_min..A_max
std::vector<uint32_t> DCI_sizes;

// fill A_order with the DCI sizes to try (DCI_sizes, or all of A_min..A_max), most frequently decoded first (ties: shorter first)
static void sort_DCI_sizes(std::vector<uint32_t> &A_order)
{
    // snapshot the counters so the comparator sees consistent values while other workers update them
    uint32_t hits[NA_max];
    for(uint32_t i = 0; i < NA_max; i++)
        hits[i] = DCI_size_hits[i].load(std::memory_order_relaxed);
    if(DCI_sizes.empty()){
        A_order.resize(NA_max);
        std::iota(A_order.begin(), A_order.end(), A_min);
    }
    else
        A_order = DCI_sizes;
    std::stable_sort(A_order.begin(), A_order.end(),
                     [&hits](uint32_t a, uint32_t b){ return hits[a-A_min] > hits[b-A_min]; });
}

// ---------- Blind PDCCH search ----------
// maximum CORESET duration tested, as max_duration in GoldenSniffer.m (TS 38.211 allows up to 3)
static const size_t PDCCH_max_duration = 2;
// amplitude correlation between the first symbols of a slot above which they belong to the same CORESET (power_detect_REGs.m)
static const float CORESET_duration_corr_thres = 0.975f;
// DM-RS PDCCH: 38.211 7.4.1.3.2 sc=1+4*k
static const int PDCCH_sc_spacing = 4;
static const int PDCCH_sc_offset = 1;
// DM-RS carriers needed to estimate c_init from a single symbol: 2*18 = 36 >= 32 bits, i.e. 1 CCE in 1 symbol
static const int PDCCH_min_N_DMRS = 18;
// with fewer DM-RS, the scrambling ID is found from all the CORESET symbols (PDCCH_DMRS_solve_underdet.m) if every symbol
// has at least this many DM-RS: 9 DM-RS = 18 bits leave 2^12 solutions per symbol (duration 2, AL1)
static const int PDCCH_underdet_min_N_DMRS = 9;
// guard on the solutions enumerated per symbol and rotation, if the equations are not independent
static const int PDCCH_underdet_max_free = 14;
// subcarriers of a CCE (6 REGs) in a CORESET with duration 1
static const int Nsc_CCE_1sym = 6*12;
// CORESET#0: a CCE is used if the phase rotation between its adjacent DM-RS (after removing the expected
// sequence) is consistent; noise or other DM-RS give about 1/sqrt(N_DMRS-1) of it
static const float CORESET0_DMRS_coherence = 0.7f;

// CORESET duration from the amplitude correlation of the first symbols of the slot starting at sym (power_detect_REGs.m)
static size_t CORESET_duration(const cfloat *Y, size_t N, size_t sym, const std::vector<bool> &is_Y_valid)
{
    size_t duration = 1;
    while(duration < PDCCH_max_duration && is_Y_valid[sym+duration]){
        // correlation between the amplitudes of the next symbol and the sum of those already in the CORESET
        double s_y = 0, s_s = 0, y_y = 0;
        for(size_t k = 0; k < N; k++){
            float s = 0;
            for(size_t l = 0; l < duration; l++)
                s += std::abs(Y[(sym+l)*N + k]);
            float y = std::abs(Y[(sym+duration)*N + k]);
            s_y += s*y;
            s_s += s*s;
            y_y += y*y;
        }
        if(s_s == 0 || y_y == 0 || s_y < CORESET_duration_corr_thres*std::sqrt(s_s*y_y))
            break;
        duration++;
    }
    return duration;
}

// c_init of the PDCCH DM-RS in symbol sym_l of the frame: (2^17*(sym_l+1)*(2*NID+1)+2*NID) mod 2^31, 38.211 7.4.1.3.1
static uint32_t PDCCH_DMRS_c_init(size_t sym_l, uint32_t NID)
{
    return ((((uint64_t)(sym_l+1)*(2*NID+1))<<17) + 2*NID) & 0x7FFFFFFF;
}

// QPSK rotation of the N_DMRS DM-RS of the unwrapped REs C from their known bits c_ref, -1 if no rotation
// matches at least 3/4 of the bits (a wrong rotation agrees on about half of them)
static int DMRS_rotation(std::vector<uint8_t> &c_est, cfloat *C, int N_DMRS, const uint8_t *c_ref)
{
    int best_rot = -1, best_agree = -1;
    for(int rot = 0; rot < 4; rot++){
        DMRS_demod(c_est, rot, C, N_DMRS, PDCCH_sc_spacing, PDCCH_sc_offset);
        int agree = 0;
        for(int i = 0; i < 2*N_DMRS; i++)
            agree += c_est[i] == c_ref[i];
        if(agree > best_agree){
            best_agree = agree;
            best_rot = rot;
        }
    }
    return 4*best_agree >= 3*2*N_DMRS ? best_rot : -1;
}

// size of the DCI formats 1_0 and 0_0 in a common search space of CORESET#0, whose size is the
// bandwidth used for the frequency domain assignment (38.212 7.3.1.0, 7.3.1.2.1)
static uint32_t fallback_DCI_size(int N_RB)
{
    uint32_t N_FDA = 0;
    while((1u<<N_FDA) < (uint32_t)(N_RB*(N_RB+1)/2))
        N_FDA++;
    return 28 + N_FDA;
}

// blind search of the PDCCH candidates in a busy run of a CORESET (blind_search_for_DCIs.m)
class PDCCHSearch {
public:
    PDCCHSearch(const Job &job, BoundedQueue<Result> &out_q, const cfloat *Y, size_t N, cfloat *C_est,
                std::vector<uint8_t> &c_est, std::vector<uint8_t> &c_ref, const std::vector<uint8_t> &x1,
                const std::vector<uint32_t> &x2_mask, std::vector<uint8_t> &x2_est, std::vector<uint8_t> &b_tilde, std::vector<uint8_t> &a_est,
                const std::vector<uint32_t> &A_order, DCIfactory &DCIcodec)
        : job(job), out_q(out_q), Y(Y), N(N), C_est(C_est), c_est(c_est), c_ref(c_ref), x1(x1), x2_mask(x2_mask), x2_est(x2_est),
          b_tilde(b_tilde), a_est(a_est), A_order(A_order), DCIcodec(DCIcodec) {}

    // search the RB-aligned run of subcarriers [k_ini,k_fin) of the CORESET with the given duration starting at symbol sym;
    // the subcarriers of the DCIs found are set in sc_used
    void run(size_t sym_, size_t duration_, int k_ini, int k_fin, std::vector<uint8_t> &sc_used_)
    {
        sc_used = &sc_used_;
        sym = sym_;
        duration = duration_;
        run_ini = k_ini;
        underdet_cache.clear();
        // phase unwrap the QPSK REs of the run in every symbol of the CORESET
        for(size_t l = 0; l < duration; l++)
            unwrap_QPSK(C_est + l*N, Y, sym+l, N, k_ini, k_fin, 1); // 1==k_step
        // DM-RS scrambling ID of the whole run, from the first symbol where it is consistent
        run_NID_ok = false;
        for(size_t l = 0; l < duration && !run_NID_ok; l++)
            run_NID_ok = find_NID(run_NID, k_ini, (k_fin-k_ini)/PDCCH_sc_spacing, l);
        search(k_ini, k_fin);
    }

    // CORESET#0 starting at symbol sym0 (38.211 7.3.2.2): interleaved CCE-to-REG mapping with REG bundles of
    // 6 REGs, R = 2 and n_shift = PCI, so every CCE is a block of 6/duration RBs scattered over the CORESET,
    // and DM-RS scrambled with the PCI and indexed from the first subcarrier of the CORESET (7.4.1.3.2).
    // The DM-RS are known, so every busy CCE is equalized with its own DM-RS based channel estimate (no
    // phase unwrap), then the candidates of every AL are put together from the CCEs. The subcarriers of the
    // DCIs found are set in sc_used.
    void search_CORESET0(size_t sym0, const CORESET0Grid &cs, float busy_power, std::vector<uint8_t> &sc_used)
    {
        const int d = cs.duration;
        const int Nsc_CCE = Nsc_CCE_1sym/d;
        const int N_DMRS = Nsc_CCE/PDCCH_sc_spacing; // per CCE and symbol
        const int E_sym = 6*N_DMRS;                  // data bits per CCE and symbol
        const int N_CCE = cs.N_RB*d/6;
        const int C = N_CCE/2;
        sym = sym0;
        duration = d;

        // REG bundle (CCE block in frequency) of every CCE: f(x) = (r*C + c + n_shift) mod (N_CCE), x = c*R + r
        std::vector<int> bundle(N_CCE);
        for(int x = 0; x < N_CCE; x++)
            bundle[x] = ((x%2)*C + x/2 + (int)cs.PCI) % N_CCE;

        // busy CCEs, from the power summed over the CORESET symbols
        std::vector<uint8_t> ok(N_CCE, 0);
        for(int j = 0; j < N_CCE; j++){
            const int k0 = cs.k_start + Nsc_CCE*bundle[j];
            float P = 0;
            for(int l = 0; l < d; l++)
                for(int k = k0; k < k0+Nsc_CCE; k++)
                    P += std::norm(Y[(sym0+l)*N + k]);
            ok[j] = P > Nsc_CCE*d*busy_power;
        }
        if(std::find(ok.begin(), ok.end(), 1) == ok.end())
            return;

        // data bits of every busy CCE whose DM-RS are the expected ones in every symbol
        CCE_bits.resize((size_t)N_CCE*d*E_sym);
        CCE_buf.resize(2*N_DMRS);
        for(int l = 0; l < d; l++){
            scrambling(c_ref.data(), 2*3*cs.N_RB, PDCCH_DMRS_c_init(sym0+l, cs.PCI));
            const cfloat *Yl = Y + (sym0+l)*N;
            for(int j = 0; j < N_CCE; j++){
                if(!ok[j])
                    continue;
                const int k0 = cs.k_start + Nsc_CCE*bundle[j];
                const uint8_t *c = c_ref.data() + 2*(Nsc_CCE*bundle[j]/PDCCH_sc_spacing);
                // least squares channel at the DM-RS: the expected ones have consistent phases across the CCE
                cfloat *h = CCE_buf.data(), *hs = h + N_DMRS;
                cfloat acc(0.0f, 0.0f);
                float norm = 0.0f;
                for(int m = 0; m < N_DMRS; m++){
                    h[m] = Yl[k0 + PDCCH_sc_spacing*m + PDCCH_sc_offset]*cfloat(1.0f-2.0f*c[2*m], -(1.0f-2.0f*c[2*m+1]));
                    if(m > 0){
                        acc += h[m]*std::conj(h[m-1]);
                        norm += std::abs(h[m])*std::abs(h[m-1]);
                    }
                }
                if(norm == 0.0f || std::abs(acc) < CORESET0_DMRS_coherence*norm){
                    ok[j] = 0;
                    continue;
                }
                // averaged over 3 DM-RS, linearly interpolated at the data REs (k0+4m, +2, +3), held at the edges
                smooth3(h, hs, N_DMRS);
                uint8_t *b = &CCE_bits[((size_t)j*d + l)*E_sym];
                for(int m = 0; m < N_DMRS; m++)
                    for(int i : {-1, 1, 2}){ // data REs around the DM-RS of index m
                        cfloat H;
                        if(i < 0)
                            H = m > 0 ? 0.75f*hs[m] + 0.25f*hs[m-1] : hs[m];
                        else if(m+1 < N_DMRS)
                            H = (1.0f-0.25f*i)*hs[m] + 0.25f*i*hs[m+1];
                        else
                            H = hs[m];
                        cfloat z = Yl[k0 + PDCCH_sc_spacing*m + PDCCH_sc_offset + i]*std::conj(H);
                        *b++ = z.real() < 0;
                        *b++ = z.imag() < 0;
                    }
            }
        }

        // DCI sizes: the fallback formats of the common search spaces first
        const uint32_t A0 = fallback_DCI_size(cs.N_RB);
        std::vector<uint32_t> sizes(1, A0);
        for(uint32_t A : A_order)
            if(A != A0)
                sizes.push_back(A);

        // candidates of CCEs n..n+AL-1 (n multiple of AL, 38.213 10.1), from the highest AL down
        std::vector<uint8_t> used(N_CCE, 0);
        std::vector<int> cces;
        for(int log2AL = log2AL_max; log2AL >= (int)log2AL_min; log2AL--){
            const int AL = 1<<log2AL;
            const uint32_t E = 108u<<log2AL;
            for(int n = 0; n + AL <= N_CCE; n += AL){
                bool usable = true;
                for(int j = n; j < n+AL; j++)
                    usable = usable && ok[j] && !used[j];
                if(!usable)
                    continue;
                // the data REs are mapped frequency first, then time (38.211 7.3.2.5)
                cces.assign(AL, 0);
                std::iota(cces.begin(), cces.end(), n);
                std::sort(cces.begin(), cces.end(), [&](int a, int b){ return bundle[a] < bundle[b]; });
                size_t pos = 0;
                for(int l = 0; l < d; l++)
                    for(int j : cces){
                        std::copy_n(&CCE_bits[((size_t)j*d + l)*E_sym], E_sym, b_tilde.begin() + pos);
                        pos += E_sym;
                    }
                for(uint32_t A : sizes){
                    if(A > A_max[log2AL] || A+16+DCI_min_redundancy > E)
                        continue;
                    uint32_t RNTI_CRC_est;
                    if(!DCIcodec.decode(RNTI_CRC_est, a_est, b_tilde, log2AL, A, cs.PCI))
                        continue;
                    push_DCI(cs.k_start + Nsc_CCE*bundle[cces.front()], cs.k_start + Nsc_CCE*(bundle[cces.back()]+1),
                             log2AL, cs.PCI, RNTI_CRC_est, A, true);
                    for(int j : cces){
                        used[j] = 1;
                        const int k0 = cs.k_start + Nsc_CCE*bundle[j];
                        std::fill(sc_used.begin() + k0, sc_used.begin() + k0 + Nsc_CCE, 1);
                    }
                    break;
                }
            }
        }
    }

private:
    const Job &job;
    BoundedQueue<Result> &out_q;
    const cfloat *Y;
    const size_t N;
    // C_est + l*N holds the unwrapped symbol sym+l of the run, starting from subcarrier run_ini
    cfloat *C_est;
    std::vector<uint8_t> &c_est, &c_ref;
    const std::vector<uint8_t> &x1;
    const std::vector<uint32_t> &x2_mask;
    std::vector<uint8_t> &x2_est, &b_tilde, &a_est;
    const std::vector<uint32_t> &A_order;
    DCIfactory &DCIcodec;

    size_t sym = 0, duration = 1;
    int run_ini = 0;
    // CORESET#0: unwrapped REs of a CCE in a symbol, data bits of every CCE and symbol
    std::vector<cfloat> CCE_buf;
    std::vector<uint8_t> CCE_bits;
    // subcarriers of the DCIs found in the slot (set by run)
    std::vector<uint8_t> *sc_used = nullptr;
    bool run_NID_ok = false;
    uint32_t run_NID = 0;
    // find_NID_underdet results (scrambling ID, -1 if none) by k0*32+N_DMRS, valid while C_est is unchanged: the
    // same short candidates come back in every pair and triple of the interleaved search and after every DCI found
    std::unordered_map<uint32_t, int32_t> underdet_cache;

    // unwrapped symbol l of the CORESET, from the RB-aligned subcarrier k0
    cfloat *C_at(int k0, size_t l) { return C_est + l*N + (k0 - run_ini); }

    // scrambling ID of the N_DMRS DM-RS starting at k0 in symbol l, trying the 4 QPSK rotations
    // TODO: the scrambling ID (here) and the DCI (DCIfactory::decode) need error-free hard decisions, which limits the
    // search at low SNR; fine for the intended line of sight, high SNR reception (soft decisions or bit flipping would help)
    bool find_NID(uint32_t &NID, int k0, int N_DMRS, size_t l)
    {
        if(N_DMRS < PDCCH_min_N_DMRS)
            return false;
        // the DM-RS index is counted from the start of the grid: Y[k0+sc_offset] contains r[k0/sc_spacing]
        size_t n_ofs = k0/PDCCH_sc_spacing;
        for(int rot = 0; rot < 4; rot++){
            DMRS_demod(c_est, rot, C_at(k0, l), N_DMRS, PDCCH_sc_spacing, PDCCH_sc_offset);
            uint32_t c_init = 0;
            if(c_init_estimate(c_est, 2*N_DMRS, n_ofs, c_init, x1, x2_est)){
                NID = (c_init>>1) & 0xFFFF;
                if(c_init == PDCCH_DMRS_c_init(sym+l, NID))
                    return true;
            }
        }
        return false;
    }

    // QPSK rotation of symbol l from its known DM-RS, -1 if no rotation matches at least 3/4 of the bits
    int match_rot(int k0, int N_DMRS, size_t l, uint32_t NID)
    {
        size_t n_ofs = k0/PDCCH_sc_spacing;
        scrambling(c_ref.data(), 2*(n_ofs+N_DMRS), PDCCH_DMRS_c_init(sym+l, NID));
        return DMRS_rotation(c_est, C_at(k0, l), N_DMRS, c_ref.data() + 2*n_ofs);
    }

    // scrambling IDs of the N_DMRS DM-RS starting at k0 in symbol l, when their 2*N_DMRS bits are fewer than the
    // 31 of c_init: all the solutions of the underdetermined system x2 = c_est ^ x1 are enumerated, and only those
    // with the structure of a PDCCH DM-RS c_init are kept, for all the 4 QPSK rotations (PDCCH_DMRS_solve_underdet.m)
    void solve_NID_underdet(std::vector<uint32_t> &NIDs, int k0, int N_DMRS, size_t l)
    {
        const size_t n_ofs = k0/PDCCH_sc_spacing;
        const int M = 2*N_DMRS;
        uint32_t row[2*PDCCH_min_N_DMRS];
        uint8_t rhs[2*PDCCH_min_N_DMRS];
        int pivot_bit[2*PDCCH_min_N_DMRS];
        NIDs.clear();
        for(int rot = 0; rot < 4; rot++){
            DMRS_demod(c_est, rot, C_at(k0, l), N_DMRS, PDCCH_sc_spacing, PDCCH_sc_offset);
            // c_init is even: only the bits 1..30 are unknown
            for(int i = 0; i < M; i++){
                row[i] = x2_mask[2*n_ofs+i] & ~1U;
                rhs[i] = c_est[i] ^ x1[2*n_ofs+i];
            }
            // reduced row echelon form
            int rank = 0;
            uint32_t pivots = 0;
            for(int b = 30; b >= 1 && rank < M; b--){
                int p = rank;
                while(p < M && !((row[p]>>b) & 1))
                    p++;
                if(p == M)
                    continue;
                std::swap(row[p], row[rank]);
                std::swap(rhs[p], rhs[rank]);
                for(int i = 0; i < M; i++)
                    if(i != rank && ((row[i]>>b) & 1)){
                        row[i] ^= row[rank];
                        rhs[i] ^= rhs[rank];
                    }
                pivot_bit[rank++] = b;
                pivots |= 1U << b;
            }
            bool consistent = true;
            for(int i = rank; i < M; i++) // the remaining rows are all zero
                if(rhs[i])
                    consistent = false;
            if(!consistent || 30-rank > PDCCH_underdet_max_free)
                continue;
            // particular solution with the free bits at zero, and a null space vector for every free bit
            uint32_t sol = 0;
            for(int r = 0; r < rank; r++)
                if(rhs[r])
                    sol |= 1U << pivot_bit[r];
            uint32_t z[30];
            int nz = 0;
            for(int b = 1; b <= 30; b++){
                if((pivots>>b) & 1)
                    continue;
                uint32_t v = 1U << b;
                for(int r = 0; r < rank; r++)
                    if((row[r]>>b) & 1)
                        v |= 1U << pivot_bit[r];
                z[nz++] = v;
            }
            // Gray code enumeration of the 2^nz solutions
            for(uint32_t i = 0; ; ){
                uint32_t NID = (sol>>1) & 0xFFFF;
                if(sol == PDCCH_DMRS_c_init(sym+l, NID))
                    NIDs.push_back(NID);
                if(++i == (1U << nz))
                    break;
                sol ^= z[__builtin_ctz(i)];
            }
        }
        std::sort(NIDs.begin(), NIDs.end());
        NIDs.erase(std::unique(NIDs.begin(), NIDs.end()), NIDs.end());
    }

    // the only scrambling ID consistent with the short DM-RS starting at k0 in all the symbols of the CORESET
    bool find_NID_underdet(uint32_t &NID, int k0, int N_DMRS)
    {
        if(N_DMRS < PDCCH_underdet_min_N_DMRS)
            return false;
        auto [it, is_new] = underdet_cache.try_emplace(k0*32 + N_DMRS, -1);
        if(!is_new){
            NID = it->second;
            return it->second >= 0;
        }
        std::vector<uint32_t> common, NIDs_l, tmp;
        for(size_t l = 0; l < duration; l++){
            solve_NID_underdet(NIDs_l, k0, N_DMRS, l);
            if(l == 0)
                common.swap(NIDs_l);
            else{
                tmp.clear();
                std::set_intersection(common.begin(), common.end(), NIDs_l.begin(), NIDs_l.end(), std::back_inserter(tmp));
                common.swap(tmp);
            }
            if(common.empty())
                return false;
        }
        if(common.size() != 1)
            return false;
        NID = it->second = common[0];
        return true;
    }

    // report the DCI in a_est, of the CORESET starting at symbol sym, spanning the subcarriers [k_ini,k_fin)
    void push_DCI(int k_ini, int k_fin, uint32_t log2AL, uint32_t NID, uint32_t RNTI, uint32_t A, bool CORESET0 = false)
    {
        count_DCI_size(A);
        Result r;
        r.kind = RES_DCI;
        r.LFN = job.LFN;
        r.rx_ofs = job.rx_ofs;
        r.slot = sym/14;
        r.k_ini = k_ini;
        r.k_fin = k_fin;
        r.AL = 1<<log2AL;
        r.duration = duration;
        r.DCI_NID = NID;
        r.RNTI = RNTI;
        r.DCI_len = A;
        r.CORESET0 = CORESET0;
        memcpy(r.DCI, a_est.data(), A);
        out_q.push(std::move(r));
    }

    // try and decode the interleaved candidate with AL 2^log2AL made of the RB-aligned parts [k0,k0+len) (in increasing
    // frequency): the scrambling ID from the first part where it is found, every part and symbol must have its DM-RS
    bool try_interleaved(const std::vector<std::pair<int,int>> &parts, uint32_t log2AL)
    {
        uint32_t NID = 0;
        bool NID_ok = false;
        for(size_t i = 0; i < parts.size() && !NID_ok; i++){
            const int N_DMRS = parts[i].second/PDCCH_sc_spacing;
            for(size_t l = 0; l < duration && !NID_ok; l++)
                NID_ok = find_NID(NID, parts[i].first, N_DMRS, l);
            if(!NID_ok && duration > 1 && N_DMRS < PDCCH_min_N_DMRS)
                NID_ok = find_NID_underdet(NID, parts[i].first, N_DMRS);
        }
        if(!NID_ok)
            return false;
        // the data REs are mapped frequency first (over all the parts), then time (38.211 7.3.2.5)
        size_t pos = 0;
        for(size_t l = 0; l < duration; l++)
            for(const auto &p : parts){
                const int N_DMRS = p.second/PDCCH_sc_spacing;
                int rot = match_rot(p.first, N_DMRS, l, NID);
                if(rot < 0)
                    return false;
                PDCCH_demod(b_tilde.data() + pos, rot, C_at(p.first, l), N_DMRS, PDCCH_sc_spacing, PDCCH_sc_offset);
                pos += 6*N_DMRS;
            }
        const uint32_t E = 108u<<log2AL;
        for(uint32_t A : A_order){
            if(A > A_max[log2AL] || A+16+DCI_min_redundancy > E)
                continue;
            uint32_t RNTI_CRC_est;
            if(DCIcodec.decode(RNTI_CRC_est, a_est, b_tilde, log2AL, A, NID)){
                push_DCI(parts.front().first, parts.back().first + parts.back().second, log2AL, NID, RNTI_CRC_est, A);
                for(const auto &p : parts)
                    std::fill(sc_used->begin() + p.first, sc_used->begin() + p.first + p.second, 1);
                return true;
            }
        }
        return false;
    }

public:
    // interleaved candidates (REG bundles of 6 REGs, so every CCE is a block of 72/duration subcarriers) made of the RB-aligned
    // runs [k_ini,k_fin) left by the non-interleaved search, sorted by frequency (enum_candidates.m): pairs of runs of the same
    // length, and triples of runs of the same length or with one as long as the other two. After a DCI is found, its
    // subcarriers are removed from the runs and the candidates are enumerated again (update_candidates.m).
    void search_interleaved(size_t sym_, size_t duration_, std::vector<std::pair<int,int>> runs, std::vector<uint8_t> &sc_used_)
    {
        sc_used = &sc_used_;
        sym = sym_;
        duration = duration_;
        run_ini = 0; // C_est holds the unwrapped REs at their subcarrier index
        run_NID_ok = false;
        underdet_cache.clear();
        const int Nsc_CCE = Nsc_CCE_1sym/duration;
        for(const auto &r : runs)
            for(size_t l = 0; l < duration; l++)
                unwrap_QPSK(C_est + l*N + r.first, Y, sym+l, N, r.first, r.second, 1);
        std::vector<std::pair<int,int>> parts;
        bool found = true;
        while(found && runs.size() >= 2){
            found = false;
            std::vector<int> len(runs.size());
            int total = 0;
            for(size_t i = 0; i < runs.size(); i++)
                total += len[i] = runs[i].second - runs[i].first;
            for(int log2AL = log2AL_max; log2AL >= 1 && !found; log2AL--){
                const int Nsc_cand = Nsc_CCE<<log2AL;
                if(Nsc_cand > total)
                    continue;
                // pairs: each candidate takes half of its subcarriers from the same position in both runs
                for(size_t i = 0; i < runs.size() && !found; i++)
                    for(size_t j = i+1; j < runs.size() && !found; j++){
                        if(len[i] != len[j] || (2*len[i]) % Nsc_cand)
                            continue;
                        const int half = Nsc_cand/2;
                        for(int c = 0; c < len[i]/half && !found; c++){
                            parts = {{runs[i].first + c*half, half}, {runs[j].first + c*half, half}};
                            found = try_interleaved(parts, log2AL);
                        }
                    }
                // triples making a single candidate
                for(size_t i = 0; i < runs.size() && !found; i++)
                    for(size_t j = i+1; j < runs.size() && !found; j++)
                        for(size_t k = j+1; k < runs.size() && !found; k++){
                            const int li = len[i], lj = len[j], lk = len[k];
                            if(li+lj+lk != Nsc_cand || !((li == lj && li == lk) || li == lj+lk || lj == li+lk || lk == li+lj))
                                continue;
                            parts = {{runs[i].first, li}, {runs[j].first, lj}, {runs[k].first, lk}};
                            found = try_interleaved(parts, log2AL);
                        }
            }
            if(!found)
                break;
            // remove the parts of the DCI found from the runs
            std::vector<std::pair<int,int>> left;
            for(const auto &r : runs){
                std::vector<std::pair<int,int>> pieces = {r};
                for(const auto &p : parts){
                    std::vector<std::pair<int,int>> next;
                    for(const auto &q : pieces){
                        const int a = std::max(q.first, p.first), b = std::min(q.second, p.first + p.second);
                        if(a >= b){
                            next.push_back(q);
                            continue;
                        }
                        if(q.first < a)
                            next.push_back({q.first, a});
                        if(b < q.second)
                            next.push_back({b, q.second});
                    }
                    pieces.swap(next);
                }
                for(const auto &q : pieces)
                    if(q.second - q.first >= Nsc_CCE)
                        left.push_back(q);
            }
            runs.swap(left);
        }
    }

private:
    // decode the non-interleaved candidate with AL 2^log2AL starting at subcarrier k0 with DM-RS scrambling ID NID
    bool decode_candidate(int k0, uint32_t log2AL, uint32_t NID)
    {
        const int Nsc_cand = (Nsc_CCE_1sym<<log2AL)/duration; // subcarriers per symbol
        const int N_DMRS = Nsc_cand/PDCCH_sc_spacing;
        const uint32_t E = 108u<<log2AL;
        const uint32_t E_sym = E/duration; // == 6*N_DMRS

        // every symbol has its own rotation; the data REs are mapped frequency first, then time (38.211 7.3.2.5)
        for(size_t l = 0; l < duration; l++){
            int rot = match_rot(k0, N_DMRS, l, NID);
            if(rot < 0)
                return false;
            PDCCH_demod(b_tilde.data() + l*E_sym, rot, C_at(k0, l), N_DMRS, PDCCH_sc_spacing, PDCCH_sc_offset);
        }

        // DCI sizes from the most frequent to the least frequent
        // TODO: at AL1, a DCI of size A is occasionally (~1e-3) also a solution for a nearby size, with a wrong RNTI
        for(uint32_t A : A_order){
            if(A > A_max[log2AL] || A+16+DCI_min_redundancy > E)
                continue;
            uint32_t RNTI_CRC_est;
            if(DCIcodec.decode(RNTI_CRC_est, a_est, b_tilde, log2AL, A, NID)){
                push_DCI(k0, k0 + Nsc_cand, log2AL, NID, RNTI_CRC_est, A);
                std::fill(sc_used->begin() + k0, sc_used->begin() + k0 + Nsc_cand, 1);
                return true;
            }
        }
        return false;
    }

    // try and decode the non-interleaved candidate with AL 2^log2AL starting at subcarrier k0
    bool try_candidate(int k0, uint32_t log2AL)
    {
        const int N_DMRS = ((Nsc_CCE_1sym<<log2AL)/duration)/PDCCH_sc_spacing;

        // scrambling ID hypotheses, each tried once: the one of the whole run, the one of the candidate alone
        // (adjacent DCIs can use different scrambling IDs, and a longer run is more exposed to bit errors),
        // and for short candidates the one from all the CORESET symbols
        uint32_t tried[3];
        int ntried = 0;
        auto try_NID = [&](uint32_t NID){
            for(int i = 0; i < ntried; i++)
                if(tried[i] == NID)
                    return false;
            tried[ntried++] = NID;
            return decode_candidate(k0, log2AL, NID);
        };
        if(run_NID_ok && try_NID(run_NID))
            return true;
        uint32_t NID;
        for(size_t l = 0; l < duration; l++)
            if(find_NID(NID, k0, N_DMRS, l))
                return try_NID(NID);
        if(duration > 1 && N_DMRS < PDCCH_min_N_DMRS && find_NID_underdet(NID, k0, N_DMRS))
            return try_NID(NID);
        return false;
    }

    // candidates in [k_ini,k_fin) from the highest AL down, aligned to the start and to the end of the range;
    // a decoded DCI is removed and the parts on its two sides are searched again (enum_candidates.m, update_candidates.m)
    void search(int k_ini, int k_fin)
    {
        const int Nsc_CCE = Nsc_CCE_1sym/duration;
        const int N_CCE = (k_fin-k_ini)/Nsc_CCE;
        if(N_CCE == 0)
            return;
        uint32_t log2AL = log2AL_max;
        while(log2AL > log2AL_min && (1<<log2AL) > N_CCE)
            log2AL--;
        for(;; log2AL--){
            const int Nsc_cand = Nsc_CCE<<log2AL;
            const int Ncand = (k_fin-k_ini)/Nsc_cand;
            for(int from_end = 0; from_end < 2; from_end++){
                if(from_end && (k_fin-k_ini)%Nsc_cand == 0) // same candidates as from the start
                    break;
                for(int i = 0; i < Ncand; i++){
                    int k0 = from_end ? k_fin-(i+1)*Nsc_cand : k_ini+i*Nsc_cand;
                    if(try_candidate(k0, log2AL)){
                        search(k_ini, k0);
                        search(k0+Nsc_cand, k_fin);
                        return;
                    }
                }
            }
            if(log2AL == log2AL_min)
                break;
        }
    }
};

// c_init of the NZP-CSI-RS in symbol sym_l of the frame: (2^10*(sym_l+1)*(2*NID+1)+NID) mod 2^31, 38.211 7.4.1.5.2
static uint32_t CSI_RS_c_init(size_t sym_l, uint32_t NID)
{
    return ((((uint64_t)(sym_l+1)*(2*NID+1))<<10) + NID) & 0x7FFFFFFF;
}

// ---------- Processing ----------
// report a NZP-CSI-RS found in symbol sym (of the frame) starting at subcarrier k0
static void push_csi(BoundedQueue<Result>& out_q, const Job& job, ResultKind kind, size_t sym, size_t k0, uint32_t NID)
{
    Result r;
    r.kind = kind;
    r.LFN = job.LFN;
    r.rx_ofs = job.rx_ofs;
    r.slot = sym/14;
    r.symbol = sym%14;
    r.k_ini = k0;
    r.DCI_NID = NID;
    r.DCI_len = 0;
    out_q.push(std::move(r));
}

void process(const Job& job, BoundedQueue<Result>& out_q, cfloat *C_est, std::vector<uint8_t> &c_est, std::vector<uint8_t> &c_ref, std::vector<uint8_t> &x1, std::vector<uint32_t> &x2_mask, std::vector<uint8_t> &x2_est, std::vector<uint8_t> &b_tilde, std::vector<uint8_t> &a_est, std::vector<uint32_t> &A_order, DCIfactory &DCIcodec) {

    const size_t NY = job.Y->size();
    const cfloat* Y = job.Y->data();
    const std::vector<cfloat> &Y_CP_corr = job.Y_CP_corr;
    const int mu = job.mu;
    const size_t Nsymb_frame = 140<<mu;
    const std::size_t N = NY/Nsymb_frame;
    const std::vector<bool> &is_Y_valid = job.is_Y_valid;
    const float noise_est = job.noise_est;
    const size_t Nsymb_slot = 14;
    const size_t Nsc_RB = 12;
    const float RE_busy_thres = 3.3;//std::log(1e3); // power ratio wrt noise power
    const float RE_free_thres = 1.8;//RE_busy_thres-1;
    size_t sym, sym_pdcch;

    // DCI sizes to try, most frequent first; refreshed once per frame
    sort_DCI_sizes(A_order);

    sym_pdcch = 0;
    for(sym = 0; sym < Nsymb_frame; sym++){
        if(sym == sym_pdcch + Nsymb_slot)
            sym_pdcch += Nsymb_slot;
        
        if(sym != sym_pdcch
        && std::abs(Y_CP_corr[sym_pdcch]) < CP_corr_thres // empty PDCCH
        && std::abs(Y_CP_corr[sym]) > CP_corr_thres){     // possible NZP-CSI-RS

            //std::cout<<"LFN: "<<LFN<<" CP_corr["<<sym<<"]: "<<Y_CP_corr[sym]<<std::endl;
            
            // FFT
            if(is_Y_valid[sym]){
            
                // test for NZP-CSI-RS density 3 - TS 38.211 Table 7.4.1.5.3-1 row 1 (noCDM)
                size_t density = 3;
                for(size_t k0 = 0; k0 < Nsc_RB/density; k0++){ // one comb offset per CSI-RS subcarrier spacing
                    uint32_t c_init = 0;
                    bool is_ok = csi_check_noCDM(c_init,0,k0,N,density,C_est,Y,sym,N,c_est,x1,x2_est);
                    if(is_ok){
                        // c_init must be the one of an NZP-CSI-RS in this symbol
                        uint32_t NID_est = c_init&0x03FF;
                        if(c_init == CSI_RS_c_init(sym, NID_est))
                            push_csi(out_q, job, RES_TRS, sym, k0, NID_est);
                    }
                }
                
                // test for NZP-CSI-RS density 1  - TS 38.211 Table 7.4.1.5.3-1 row 2 (noCDM)
                density = 1; // or 0 for 0.5
                for(size_t k0 = 0; k0 < Nsc_RB; k0++){
                    uint32_t c_init = 0;
                    bool is_ok = csi_check_noCDM(c_init,0,k0,N,density,C_est,Y,sym,N,c_est,x1,x2_est);
                    if(is_ok){
                        // c_init must be the one of an NZP-CSI-RS in this symbol
                        uint32_t NID_est = c_init&0x03FF;
                        if(c_init == CSI_RS_c_init(sym, NID_est))
                            push_csi(out_q, job, RES_CSI, sym, k0, NID_est);
                    }
                }
                
                // test for NZP-CSI-RS density 1 - TS 38.211 Table 7.4.1.5.3-1 rows 3,4,5,6,7,9,11,13,16 (CDM2: fd2)
                // TODO (RES_CSI_FD2): the MATLAB algorithm, a channel estimate per RB for each relative phase of the
                // two CDM REs, is in the git history of csi_check_CDM2 (removed: it always returned false)
                
                // test for NZP-CSI-RS density 1 - TS 38.211 Table 7.4.1.5.3-1 rows 8,10,12,14,17 (CDM4: fd2,td2)
                // TODO
                
                // test for NZP-CSI-RS density 1 - TS 38.211 Table 7.4.1.5.3-1 rows 15,18 (CDM8: fd2,td4)
                // TODO
            }
        }
    }
    
    // the frames tracked before the first SS block have no noise estimate: every subcarrier would be busy and
    // the whole band a single run, so their PDCCH is not searched (the NZP-CSI-RS search above doesn't need it)
    if(!(noise_est > 0.0f))
        return;

    // blind PDCCH search in the CORESET at the beginning of every slot
    PDCCHSearch pdcch(job, out_q, Y, N, C_est, c_est, c_ref, x1, x2_mask, x2_est, b_tilde, a_est, A_order, DCIcodec);
    const CORESET0Grid &cs0 = job.cs0;
    // TODO: sc_used_type0 is never read nor cleared (the Type0-PDCCH CORESET#0 doesn't start at symbol 0, so it can't
    // overlap the run search below)
    std::vector<uint8_t> sc_used(N), sc_used_type0(N);
    // a CORESET starting at symbol s of the frame is searched if its symbols are there and s is not empty
    auto CORESET_ok = [&](size_t s, size_t d){
        if(s+d > Nsymb_frame || std::abs(Y_CP_corr[s]) <= CP_corr_thres)
            return false;
        for(size_t l = 0; l < d; l++)
            if(!is_Y_valid[s+l])
                return false;
        return true;
    };
    for(sym = 0; sym < Nsymb_frame; sym += Nsymb_slot){
        // CORESET#0 (known from the MIB): from symbol 0 for all the common search spaces, and in the Type0-PDCCH
        // slots from the first symbol given by searchSpaceZero; its DCIs are left out of the search below
        std::fill(sc_used.begin(), sc_used.end(), 0);
        if(cs0.valid){
            if(CORESET_ok(sym, cs0.duration))
                pdcch.search_CORESET0(sym, cs0, RE_busy_thres*noise_est, sc_used);
            const int slot = sym/Nsymb_slot;
            if(cs0.first_symbol > 0 && job.LFN%2 == cs0.frame_parity && (slot == cs0.n0 || slot == cs0.n0+1)
            && cs0.first_symbol + cs0.duration <= (int)Nsymb_slot && CORESET_ok(sym+cs0.first_symbol, cs0.duration))
                pdcch.search_CORESET0(sym+cs0.first_symbol, cs0, RE_busy_thres*noise_est, sc_used_type0);
        }

        if(std::abs(Y_CP_corr[sym]) <= CP_corr_thres || !is_Y_valid[sym])
            continue;

        // the CORESET duration estimated from the amplitudes (power_detect_REGs.m) fails when the PDCCH is a small part of
        // the band, and PDSCH in the next symbol can fake it: both durations are searched, the estimated one first, each
        // on the subcarriers left by the DCIs already found
        const size_t d_est = CORESET_duration(Y, N, sym, is_Y_valid);
        // the 6-RB groups of a CORESET start at a CRB multiple of 6 (38.211 7.3.2.2), so its CCEs of 6/duration RBs
        // are aligned to the grid, assuming it starts at CRB 0 as for the DM-RS index (power_detect_REGs.m)
        // TODO: rb-Offset-r16 is not aligned to this grid
        // runs of busy subcarriers: from the power of the REs summed over the CORESET duration, or, for a duration not
        // given by the estimate, busy in every symbol (a longer CORESET only where all its symbols are occupied)
        auto find_runs = [&](size_t duration, bool every_symbol){
            const int Nsc_CCE = Nsc_CCE_1sym/duration;
            std::vector<std::pair<int,int>> runs;
            bool RE_busy = false;
            int k_ini = 0;
            for(size_t k = 0; k <= N; k++){
                bool busy = false, free = true; // k == N closes a run reaching the band edge
                if(k < N && !sc_used[k]){
                    float P = 0, P_min = INFINITY;
                    for(size_t l = 0; l < duration; l++){
                        const float P_l = std::norm(Y[(sym+l)*N + k]);
                        P += P_l;
                        P_min = std::min(P_min, P_l);
                    }
                    busy = every_symbol ? P_min > RE_busy_thres*noise_est : P > duration*RE_busy_thres*noise_est;
                    free = every_symbol ? P_min < RE_free_thres*noise_est : P < duration*RE_free_thres*noise_est;
                }
                if(!RE_busy && busy){
                    k_ini = k;
                    RE_busy = true;
                }
                else if(RE_busy && free){
                    RE_busy = false;
                    // PDCCHs are made of whole CCEs: round the run edges to the closest CCE boundary
                    int run_ini = (k_ini + Nsc_CCE/2)/Nsc_CCE*Nsc_CCE;
                    int run_fin = ((int)k + Nsc_CCE/2)/Nsc_CCE*Nsc_CCE;
                    if(run_fin - run_ini >= Nsc_CCE)
                        runs.push_back({run_ini, run_fin});
                }
            }
            return runs;
        };
        // non-interleaved candidates in every run, then interleaved ones over the runs left
        auto search_duration = [&](size_t duration, bool every_symbol){
            for(const auto &r : find_runs(duration, every_symbol))
                pdcch.run(sym, duration, r.first, r.second, sc_used);
            pdcch.search_interleaved(sym, duration, find_runs(duration, every_symbol), sc_used);
        };
        search_duration(d_est, false);
        for(size_t d = 1; d <= PDCCH_max_duration; d++)
            if(d != d_est && CORESET_ok(sym, d))
                search_duration(d, d > d_est);
    }
}

// ---------- Worker ----------
std::mutex fftw_plan_mutex;
int DCIdecoder_n_workers = std::thread::hardware_concurrency()/2;

void worker(BoundedQueue<Job>& in_q, BoundedQueue<Result>& out_q, size_t N_FFT, const DCIfactory &DCIcodec_shared) {

    // --worker-cpus (default P-cores 8,10,12,14 + E-cores 16:31; --rx-cpus are left to streaming and FFT)
    pin_self(worker_cpus);

    // shares the matrices, with its own decoding buffers
    DCIfactory DCIcodec(DCIcodec_shared);

    // phase unwrap / scrambling ID guessing, one symbol per CORESET symbol
    fftwf_complex *C_est_fftw = fftwf_alloc_complex(PDCCH_max_duration*N_FFT);
    cfloat *C_est = reinterpret_cast<cfloat*>(C_est_fftw);
    std::vector<uint8_t> c_est(2*N_FFT);
    // known DM-RS sequence
    std::vector<uint8_t> c_ref(2*N_FFT);
    // the x1 sequence is precomputed over this many points to handle a worst case
    std::vector<uint8_t> x1(2*N_FFT);
    scrambling_x1(x1.data(), 2*N_FFT);
    // and x2 as a linear function of c_init, for the short DM-RS
    std::vector<uint32_t> x2_mask(2*N_FFT);
    scrambling_x2_masks(x2_mask.data(), 2*N_FFT);
    std::vector<uint8_t> x2_est(2*N_FFT+1631);
    
    // DCI / RNTI decoding buffers
    // a single candidate, up to the 108*16 bits of AL16
    std::vector<uint8_t> b_tilde(108*16);
    std::vector<uint8_t> a_est(140);
    std::vector<uint32_t> A_order;

    // loop
    while (true) {
        std::optional<Job> job_opt = in_q.pop();
        if (!job_opt) break;
        process(*job_opt, out_q, C_est, c_est, c_ref, x1, x2_mask, x2_est, b_tilde, a_est, A_order, DCIcodec); // plus the thread's variables
        Result done;
        done.kind = RES_FRAME_DONE;
        done.rx_ofs = job_opt->rx_ofs;
        done.LFN = job_opt->LFN;
        job_opt.reset();   // the buffer goes back to the pool right away
        out_q.push(std::move(done));
    }
}

// ---------- FFT worker ----------
// the symbols process() reads, so the others are not demodulated: in every slot whose first symbol has CP
// correlation, the CORESET symbols from symbol 0 (both durations searched, and CORESET#0) and the Type0-PDCCH
// CORESET#0 of the slot; in a slot with an empty first symbol, the symbols with CP correlation (NZP-CSI-RS)
static void symbols_to_demodulate(const FFTJob &fj, size_t Nsymb_frame, std::vector<uint8_t> &need)
{
    need.assign(Nsymb_frame, 0);
    const CORESET0Grid &cs0 = fj.cs0;
    const bool pdcch = fj.noise_est > 0.0f; // no PDCCH search before the first noise estimate
    auto busy = [&](size_t s){ return s < Nsymb_frame && std::abs(fj.Y_CP_corr[s]) > CP_corr_thres; };
    auto mark = [&](size_t s, size_t d){
        for(size_t l = 0; l < d && s+l < Nsymb_frame; l++)
            need[s+l] = 1;
    };
    for(size_t sym = 0; sym < Nsymb_frame; sym += Nsymb_slot){
        if(!busy(sym)){
            for(size_t l = 1; l < Nsymb_slot; l++)
                if(busy(sym+l))
                    need[sym+l] = 1;
        }
        else if(pdcch)
            mark(sym, std::max(PDCCH_max_duration, cs0.valid ? (size_t)cs0.duration : (size_t)0));
        // the Type0-PDCCH CORESET#0 is searched whether or not symbol 0 is busy
        const int slot = sym/Nsymb_slot;
        if(pdcch && cs0.valid && cs0.first_symbol > 0 && fj.LFN%2 == cs0.frame_parity && (slot == cs0.n0 || slot == cs0.n0+1)
        && busy(sym + cs0.first_symbol))
            mark(sym + cs0.first_symbol, cs0.duration);
    }
}

void fft_worker(BoundedQueue<FFTJob>& fft_q, BoundedQueue<Job>& in_q, BufferPool& pool,
                size_t N_FFT, size_t Nsymb_frame, size_t CP_slack, size_t N, int grid_ofs,
                double Ts, const std::vector<float>& phase_comp, const std::vector<size_t>& sym_ofs)
{
    // --rx-cpus (default P-cores 0,2,4,6), shared with the RX thread
    pin_self(rx_cpus);

    const size_t NY = Nsymb_frame * N_FFT;
    fftwf_complex* Y_local_fftw = fftwf_alloc_complex(NY);
    cfloat*        Y_local      = reinterpret_cast<cfloat*>(Y_local_fftw);
    fftwf_complex* tmpbuf_fftw  = fftwf_alloc_complex(N_FFT);
    cfloat*        tmpbuf       = reinterpret_cast<cfloat*>(tmpbuf_fftw);

    std::vector<fftwf_plan> plans(Nsymb_frame);
    {
        std::lock_guard<std::mutex> lk(fftw_plan_mutex);
        for (size_t n = 0; n < Nsymb_frame; n++)
            plans[n] = fftwf_plan_dft_1d(N_FFT,
                           Y_local_fftw + n * N_FFT, tmpbuf_fftw,
                           FFTW_FORWARD, FFTW_MEASURE);
    }

    std::vector<uint8_t> need;
    while (true) {
        std::optional<FFTJob> job_opt = fft_q.pop();
        if (!job_opt) break;
        FFTJob& fj = *job_opt;

        // only the symbols process() reads; those left invalid aren't copied to the job either
        symbols_to_demodulate(fj, Nsymb_frame, need);
        for (size_t sym = 0; sym < Nsymb_frame; sym++)
            if (fj.is_Y_valid[sym] || need[sym])
                std::copy(fj.Y.begin() + sym * N_FFT, fj.Y.begin() + (sym + 1) * N_FFT, Y_local + sym * N_FFT);

        for (size_t sym = 0; sym < Nsymb_frame; sym++) {
            if (!fj.is_Y_valid[sym] && need[sym]) {
                ofdm_fft(sym, Y_local, phase_comp[sym], fj.Cfo_eF,
                         sym_ofs[sym], Ts, N_FFT, plans[sym], tmpbuf, CP_slack);
                fj.is_Y_valid[sym] = true;
            }
        }

        BufferPool::Handle buf = pool.acquire();
        if (!buf) break;
        buf->set_size(N * Nsymb_frame);
        cfloat* p = buf->data();
        for (size_t sym = 0; sym < Nsymb_frame; sym++)
            if (fj.is_Y_valid[sym])
                for (size_t i = 0; i < N; i++)
                    p[sym * N + i] = Y_local[sym * N_FFT + grid_ofs + i];

        Job job(fj.LFN, fj.rx_ofs, std::move(buf), fj.Y_CP_corr, fj.mu, fj.is_Y_valid, fj.noise_est);
        job.cs0 = fj.cs0;
        in_q.push(std::move(job));
    }

    {
        std::lock_guard<std::mutex> lk(fftw_plan_mutex);
        for (size_t n = 0; n < Nsymb_frame; n++)
            fftwf_destroy_plan(plans[n]);
    }
    fftwf_free(tmpbuf_fftw);
    fftwf_free(Y_local_fftw);
}

bool csi_check_noCDM(uint32_t &c_init, size_t grid_ofs, size_t k0, size_t N, size_t density, cfloat *C_est, const cfloat *Y, size_t sym,
               size_t N_FFT, std::vector<uint8_t> &c_est, const std::vector<uint8_t> &x1, std::vector<uint8_t> &x2_est){
    int k_ini = grid_ofs + k0;
    int k_fin = k_ini + N;
    const size_t Nsc_RB = 12;
    int sc_spacing = density ? Nsc_RB/density : 2*Nsc_RB;
    int Nsc_CSI = N/sc_spacing;
    unwrap_QPSK(C_est,Y,sym,N_FFT,k_ini,k_fin,sc_spacing);
    int c_est_len = 2*Nsc_CSI;
    bool is_ok = false;
    for(int rot = 0; rot < 4; rot++){
        DMRS_demod(c_est,rot,C_est,Nsc_CSI,1,0); // 0:1:Nsc_CSI
        is_ok = c_init_estimate(c_est,c_est_len,0,c_init,x1,x2_est);
        if(is_ok)
            break;
    }
    return is_ok;
}
