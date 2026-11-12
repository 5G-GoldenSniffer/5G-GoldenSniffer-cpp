#ifndef _DCIFACTORY_H
#define _DCIFACTORY_H 1

#include <cstdint>
#include <memory>
#include <vector>
#include "gf2solve.h"
#include "polar.h"

const uint32_t log2AL_min = 0;
const uint32_t log2AL_max = 4;
const uint32_t N_log2AL = log2AL_max-log2AL_min+1;
const uint32_t A_min = 20;
const uint32_t A_max[N_log2AL] = {84, 140, 140, 140, 140};
const uint32_t NA_max = A_max[N_log2AL-1]-A_min+1;
// minimum number of redundant bits E-(A+16) for a (AL, A) hypothesis: decode() accepts any consistent
// solution, so a hypothesis with r redundant bits is a false positive with probability ~2^-r
const uint32_t DCI_min_redundancy = 16;

// decoding matrices of every (AL, A) hypothesis: built by DCIfactory::init, then read-only and shared by its copies
struct DCItables {
    mzd_t *G_A   [N_log2AL][NA_max];
    mzp_t *P_A   [N_log2AL][NA_max];
    mzd_t *G_B   [N_log2AL][NA_max];
    mzp_t *P_B   [N_log2AL][NA_max];
    mzv_t *b0    [N_log2AL][NA_max];
};

class DCIfactory {
public:
    // build the DCI encoding/decoding matrices; the copies of this object share them, each with its own buffers
    void init(void);
    // runs a random input self-test
    bool selftest(bool verbose);
    // dci encode the message a(Mbit) into b_tilde, with parameters rnti_crc, rnti_scr, N_ID, with initializations crc_init(=0xffffff) and x1_init(=0x1)
    void encode(std::vector<uint8_t> &b_tilde, const std::vector<uint8_t> &a, uint32_t rnti_crc, uint32_t rnti_scr, uint32_t N_ID, uint32_t Mbit, uint32_t crc_init, uint32_t x1_init);
    // try and decode b_tilde as a DCI with AL 2^log2AL, message length A, DMRS_scrID N_ID
    bool decode(uint32_t &RNTI_CRC, std::vector<uint8_t> &a, const std::vector<uint8_t> &b_tilde, uint32_t log2AL, uint32_t A, uint32_t N_ID);

private:
    std::shared_ptr<const DCItables> tables;

    // DCI encoding buffers (the polar code itself is in polar.cpp)
    // maximum size of a DCI: 140 bits message + 24 bits CRC
    static constexpr uint32_t DCIEnc_Kmax = 164;
    uint8_t DCIEnc_c[DCIEnc_Kmax];

    static constexpr uint32_t DCIEnc_Emax = 1728; // AL16(CCE) * 6(RB per CCE) * 12(sc per RB) * 3/4 (data subcarriers) * 2 (bits per subcarrier)
    uint8_t DCIEnc_f[DCIEnc_Emax];
    uint8_t DCIEnc_c_DCI[DCIEnc_Emax];
    uint8_t DCIEnc_x1[DCIEnc_Emax];
    // polar code of the last encode, reused while K and E don't change
    std::shared_ptr<const PolarCode> polar;

    // DCI decoding: the x2 scrambling sequence of the last N_ID of every AL (b_NID_key = N_ID+1, 0 if none), and the
    // right-hand side b
    uint32_t b_NID_key[N_log2AL] = {};
    uint8_t b_NID[N_log2AL][DCIEnc_Emax];
    uint8_t b[DCIEnc_Emax];
};

#endif
