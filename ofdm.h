#ifndef OFDM_H_
#define OFDM_H_

#include <cstddef>
#include <cstdint>
#include <vector>
#include <fftw3.h>
#include "nr_common.h"

// FFT with Cfo compensation, CP correction and discard
void ofdm_fft(size_t sym, cfloat *Y, float phase_comp, float Cfo_eF, size_t sym_ofs, double Ts, size_t N_FFT, fftwf_plan Y_symb_to_freq, cfloat *Y_tmpbuf, size_t CP_slack);
// QPSK phase unwrap
void unwrap_QPSK(cfloat *C_est, const cfloat *Y, size_t sym, size_t N_FFT, int k_ini, int k_fin, int k_step);
// demodulate a DM-RS sequence with an hypothesis on the phase rotation exp(1i*rot*pi/2)
void DMRS_demod(std::vector<uint8_t> &c_est, int rot, cfloat *C_est, int N_DMRS, int sc_spacing, int sc_offset);
// channel estimates at n DM-RS averaged over 3 adjacent ones (2 at the ends): hs[m] = mean(h[m-1..m+1])
inline void smooth3(const cfloat *h, cfloat *hs, size_t n)
{
    for(size_t m = 0; m < n; m++){
        const size_t m0 = m > 0 ? m-1 : 0, m1 = m+1 < n ? m+1 : n-1;
        cfloat s(0.0f, 0.0f);
        for(size_t i = m0; i <= m1; i++)
            s += h[i];
        hs[m] = s/(float)(m1-m0+1);
    }
}
// demodulate the PDCCH for the PDCCH with the DM-RS above
void PDCCH_demod(uint8_t *b_tilde, int rot, cfloat *C_est, int N_DMRS, int sc_spacing, int sc_offset);

#endif
