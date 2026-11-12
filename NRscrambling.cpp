#include <cstdint>
#include <vector>
#include "NRscrambling.h"

void scrambling_x1(uint8_t *c, size_t M_PN)
{
    uint32_t x1 = 0x5e485840;
    for(size_t n = 0; n < M_PN; n++) {
        c[n] = x1 & 0x01;
        uint8_t fb_x1 = ((x1 >> 3) ^ x1) & 0x01;
        x1 = (x1 >> 1) | (fb_x1 << 30);
    }
}

void scrambling(uint8_t *c, size_t M_PN, uint32_t c_init)
{
    uint32_t x1 = 0x5e485840;
    uint32_t x2 = 0;
    static const uint32_t x2_NC[31]={
        0x0099110e,0x0132221c,0x02644438,0x04c88870,0x099110e0,0x132221c0,0x26444380,0x4c888700,
        0x19110e0f,0x32221c1e,0x6444383c,0x48887077,0x1110e0e1,0x2221c1c2,0x44438384,0x08870707,
        0x110e0e0e,0x221c1c1c,0x44383838,0x0870707f,0x10e0e0fe,0x21c1c1fc,0x438383f8,0x070707ff,
        0x0e0e0ffe,0x1c1c1ffc,0x38383ff8,0x70707ff0,0x60e0ffef,0x41c1ffd1,0x0383ffad};
    for(int i = 0; i < 31; i++) {
        if(__builtin_parity(x2_NC[i] & c_init))
            x2 |= (1U << i);
    }
    for(size_t n = 0; n < M_PN; n++) {
        c[n] = (x1 ^ x2) & 0x01;
        uint8_t fb_x1 = ((x1 >> 3) ^ x1) & 0x01;
        uint8_t fb_x2 = __builtin_parity(x2 & 0x0F);
        x1 = (x1 >> 1) | (fb_x1 << 30);
        x2 = (x2 >> 1) | (fb_x2 << 30);
    }
}

void scrambling_x2_masks(uint32_t *mask, size_t M_PN)
{
    const size_t Nc = 1600;
    std::vector<uint32_t> m(Nc+M_PN);
    // the initial state of x2 is c_init: x2[i] = c_init bit i
    for(size_t i = 0; i < 31 && i < Nc+M_PN; i++)
        m[i] = 1U << i;
    // x2(n+31) = x2(n+3)+x2(n+2)+x2(n+1)+x2(n)
    for(size_t n = 31; n < Nc+M_PN; n++)
        m[n] = m[n-28]^m[n-29]^m[n-30]^m[n-31];
    for(size_t n = 0; n < M_PN; n++)
        mask[n] = m[Nc+n];
}

bool c_init_estimate(const std::vector<uint8_t> &c_est, size_t c_est_len, size_t n_ofs, uint32_t &c_init, const std::vector<uint8_t> &x1, std::vector<uint8_t> &x2_est)
{
    size_t x2_ofs = 1600+2*n_ofs;
    for(size_t i = 0; i < c_est_len; i++)
        x2_est[x2_ofs+i] = c_est[i] ^ x1[2*n_ofs+i];

    // does it satisfy the x2 recursion?
    bool is_ok;
    if(c_est_len < 32)
        is_ok = false;
    else{
        is_ok = true;
        for(int i = c_est_len-32; i>=0; i--){
            if(x2_est[x2_ofs+i+31]^x2_est[x2_ofs+i+3]^x2_est[x2_ofs+i+2]^x2_est[x2_ofs+i+1]^x2_est[x2_ofs+i]){
                is_ok = false;
                break;
            }
        }
    }
    if(is_ok){
        // recover c_init
        for(int i = x2_ofs-1; i >= 0; i--)
            x2_est[i] = x2_est[i+31]^x2_est[i+3]^x2_est[i+2]^x2_est[i+1];
        
        c_init = 0;
        for(int i = 0; i < 31; i++)
            c_init |= x2_est[i]<<i;
    }
    return is_ok;
}
