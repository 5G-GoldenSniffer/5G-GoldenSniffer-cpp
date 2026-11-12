#ifndef PBCH_H_
#define PBCH_H_

#include <cstdint>
#include <string>
#include <vector>
#include "nr_common.h"
#include "polar.h"

// subcarriers and symbols of an SS/PBCH block
const size_t SSB_carriers = 240;
const size_t SSB_symbols = 4;

// content of a PBCH: MIB (38.331) plus the timing bits of the PBCH payload (38.212 7.1.1)
struct MIB {
    uint32_t SFN;            // system frame number of the frame carrying the SS block
    int HRF;                 // half frame of the SS block
    int SSB_index;           // SS block index (from the PBCH DM-RS)
    int scs_common;          // subCarrierSpacingCommon in FR1: 0 (15 kHz) or 1 (30 kHz)
    int k_SSB;               // ssb-SubcarrierOffset plus its MSB from the PBCH payload (FR1), 15 kHz units
    int dmrs_typeA_pos;      // dmrs-TypeA-Position: 2 or 3
    int coreset0;            // pdcch-ConfigSIB1.controlResourceSetZero
    int ss0;                 // pdcch-ConfigSIB1.searchSpaceZero
    bool cell_barred;
    bool intra_freq_resel;   // intraFreqReselection allowed
};

// first symbol of SS block SSB_index in its half frame (38.213 4.1, Case A for mu=0, Case C for mu=1)
int SSB_first_symbol(int mu, int SSB_index);
// c_init of the PBCH DM-RS of cell PCI for i_SSB_bar (38.211 7.4.1.4.1)
uint32_t PBCH_DMRS_c_init(uint32_t PCI, uint32_t i_bar);

// CORESET#0 and Type0-PDCCH CSS from pdcch-ConfigSIB1 (38.213 13, FR1, bands with minimum channel bandwidth
// 5 or 10 MHz), for a receiver grid of N subcarriers at the SS block subcarrier spacing whose center is k120
// subcarriers below subcarrier 120 of the SS block
struct CORESET0 {
    bool valid;              // false if k_SSB > 23 (no CORESET#0) or the table entry is reserved
    int pattern;             // SS/PBCH block and CORESET multiplexing pattern
    int N_RB;                // CORESET size [RB]
    int duration;            // CORESET symbols
    int offset_RB;           // from the first RB of the CORESET to the common RB overlapping the first RB of the SS block
    bool same_scs;           // the CORESET has the SS block subcarrier spacing, so it is in the receiver grid
    double k_start;          // first subcarrier of the CORESET in the receiver grid (same_scs only)
    double CRB_grid_ofs;     // first subcarrier of the grid modulo 12 with respect to the common RBs (0: aligned)
    // Type0-PDCCH monitoring occasions, multiplexing pattern 1 (38.213 Table 13-11): slots n0 and n0+1
    // of the frames with SFN mod 2 == frame_parity, CORESET starting at first_symbol
    int n0;
    int frame_parity;
    int first_symbol;
};
bool CORESET0_config(const MIB &mib, int mu, int k120, size_t N, CORESET0 &cs);
// one line on the MIB and one on CORESET#0, for the console
std::string MIB_to_string(const MIB &mib, int mu, int k120, size_t N);

// PBCH decoding (38.211 7.3.3, 7.4.1.4; 38.212 7.1): DM-RS based channel estimate, CRC-aided list
// decoding of the BCH; the SS block index is found from the DM-RS
class PBCHDecoder {
public:
    explicit PBCHDecoder(int mu); // L_max 4 for mu=0 (Case A below 3 GHz), 8 for mu=1 (Case C)

    // ssb[l*SSB_carriers+k]: subcarrier k of symbol l (0: PSS) of the SS block of cell PCI
    bool decode(const cfloat *ssb, uint32_t PCI, MIB &mib);
    // random MIBs through a noisy channel
    bool selftest(bool verbose);
    // PBCH and its DM-RS of a MIB (only SFN, HRF, SSB_index, k_SSB and the MIB fields are used)
    void encode(cfloat *ssb, const MIB &mib, uint32_t PCI);

private:
    int L_max;
    PolarCode pc; // K = 56, E = 864
    uint32_t c_PCI = 0xffffffff;
    std::vector<uint8_t> c; // scrambling sequence with c_init = PCI, for all the v
    const std::vector<uint8_t> &scrambling_PCI(uint32_t PCI);
};

#endif
