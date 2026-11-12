#ifndef POLAR_H_
#define POLAR_H_

#include <cstdint>
#include <vector>

// CRC24C of 38.212 5.1 over a[0..A-1], with the shift register initialized to init
// (0xffffff is the same as the 24 ones prepended to a DCI, 38.212 7.3.2; 0 for the BCH)
uint32_t crc24c(const uint8_t *a, uint32_t A, uint32_t init);

// polar code of 38.212 5.3.1 without parity check bits (DCI, BCH), with the rate matching of
// 5.4.1.1-5.4.1.2 and no coded bits interleaving (I_BIL = 0)
class PolarCode {
public:
    // K bits (CRC included) into E bits, N = 2^n with n <= n_max; I_IL: input bits interleaving (5.3.1.1)
    PolarCode(uint32_t K, uint32_t E, uint32_t n_max, bool I_IL);

    uint32_t K, E, N;

    // c[K] -> e[E]
    void encode(const uint8_t *c, uint8_t *e) const;
    // successive cancellation list decoding of the LLRs llr_e[E] (> 0: bit 0), min-sum approximation:
    // up to L candidate c[K], most likely first (the caller checks the CRC)
    void decode_list(const float *llr_e, int L, std::vector<std::vector<uint8_t>> &c_list) const;

private:
    uint32_t n;
    std::vector<uint16_t> PI;    // input interleaver: c'[k] = c[PI[k]]
    std::vector<uint16_t> J;     // sub-block interleaver: y[n] = d[J[n]]
    std::vector<uint16_t> info;  // positions of u carrying c'[0..K-1], increasing
    std::vector<uint8_t> frozen; // N flags
    enum { REPETITION, PUNCTURING, SHORTENING } rm;
};

#endif
