#ifndef DCI_SCAN_H_
#define DCI_SCAN_H_

#include <cstdint>
#include <cstddef>
#include <vector>
#include <mutex>
#include "nr_common.h"
#include "pipeline.h"

// tests if there is a NZP-CSI-RS (row 1 and 2) with density "density" starting at carrier k0 in symbol sym, with a lot of temporary storage variables
bool csi_check_noCDM(uint32_t &c_init, size_t grid_ofs, size_t k0, size_t N, size_t density, cfloat *C_est, const cfloat *Y, size_t sym, size_t N_FFT, std::vector<uint8_t> &c_est, const std::vector<uint8_t> &x1, std::vector<uint8_t> &x2_est);

// DCI sizes tried by the blind search, sorted (--dci-sizes); empty: all sizes A_min..A_max
extern std::vector<uint32_t> DCI_sizes;

// ---------- Worker ----------
extern std::mutex fftw_plan_mutex;
extern int DCIdecoder_n_workers;

class DCIfactory;
// DCIcodec: initialized and self-tested, each worker decodes with a copy (shared matrices, own buffers)
void worker(BoundedQueue<Job>& in_q, BoundedQueue<Result>& out_q, size_t N_FFT, const DCIfactory &DCIcodec);

void fft_worker(BoundedQueue<FFTJob>& fft_q, BoundedQueue<Job>& in_q, BufferPool& pool,
                size_t N_FFT, size_t Nsymb_frame, size_t CP_slack, size_t N, int grid_ofs,
                double Ts, const std::vector<float>& phase_comp, const std::vector<size_t>& sym_ofs);

#endif
