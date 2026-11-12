#ifndef _NRSCRAMBLING_H
#define _NRSCRAMBLING_H 1

#include <cstddef>
#include <cstdint>
#include <vector>

// scrambling sequence according to 38.211 5.2.1
void scrambling(uint8_t *c, size_t M_PN, uint32_t c_init);
// only the sequence x1
void scrambling_x1(uint8_t *c, size_t M_PN);
// x2[1600+n] (n < M_PN) as a linear function of c_init: bit i of mask[n] is set if c_init bit i contributes
void scrambling_x2_masks(uint32_t *mask, size_t M_PN);
// estimate c_init from c_est (using the cached x1, the offset n_ofs and x2_est as temporary storage)
bool c_init_estimate(const std::vector<uint8_t> &c_est, size_t c_est_len, size_t n_ofs, uint32_t &c_init, const std::vector<uint8_t> &x1, std::vector<uint8_t> &x2_est);

#endif
