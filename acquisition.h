#ifndef ACQUISITION_H_
#define ACQUISITION_H_

#include <vector>
#include <mutex>
#include <condition_variable>
#include "nr_common.h"

// -------------------------------------- ACQUISITION -----------------------------------------

// output of SSB_search
typedef struct {
    bool valid;
    double f_SSB;
    long long ofs;
    size_t NID2; // PSS
    size_t NID1; // SSS
    float SIR_est_dB;
    float Cfo_est;
    size_t symSSB;
} SSB_info_t;

const size_t N_PSS = 3;
const size_t N_SSS = 336;
const size_t SS_carriers = 128; // with a null first subcarrier, so 55:182 instead of 56:182
const size_t PSS_len = 128+9; // including the cyclic prefix
extern std::vector<std::vector<cfloat>> PSS_freq;
extern std::vector<std::vector<cfloat>> PSS_time;
extern std::vector<float> PSS_time_energy;
extern std::vector<std::vector<std::vector<cfloat>>> SSS_freq;

// SSB
void SS_init(void);
// ACQ_window_size: search window in samples after the decimation by decim (sample rate / (128*SCS))
// n_threads search threads: the caller of SSB_search and n_threads-1 helpers (joined by ACQ_corr_done)
void ACQ_corr_init(size_t ACQ_window_size, size_t decim, int n_threads);
void ACQ_corr_done(void);
const float PSS_threshold = 0.25;
const float SSS_threshold = 0.3;
const float H_est_SS_corr_thres = 0.6;
SSB_info_t SSB_search(const std::vector<cfloat> *r, long long r_ofs,
                double freq_rx, double bandwidth, double samp_rate, int mu, bool benchmark);

// Persistent worker thread
// All the following variables are protected by worker_mutex.
extern std::mutex worker_mutex;
extern std::condition_variable worker_cv;
extern std::condition_variable result_cv;

// Worker thread state
enum WorkerState { IDLE, PROCESSING, RESULT_READY };
extern WorkerState worker_state;
extern bool worker_should_exit;

// Pointer to the buffer to process (protected by worker_mutex)
extern std::vector<cfloat>* buffer_to_process;
extern long long buffer_to_process_rx_ofs;

// Parameters for SSB_search (protected by worker_mutex)
extern double worker_freq_rx;
extern double worker_bandwidth;
extern double worker_samp_rate;
extern int worker_mu;

// Result (protected by worker_mutex)
extern SSB_info_t worker_result;

// worker thread loop
void SSB_worker_thread();

#endif
