#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <atomic>
#include <string>
#include <algorithm>
#include <pthread.h>
#include <fftw3.h>
#include "nr_common.h"
#include "LP_order_32.h"
#include "NRscrambling.h"
#include "gscn.h"
#include "util.h"
#include "acquisition.h"
#include "pbch.h"

std::vector<std::vector<cfloat>> PSS_freq(N_PSS,std::vector<cfloat>(SS_carriers));
std::vector<std::vector<cfloat>> PSS_time(N_PSS,std::vector<cfloat>(PSS_len));
std::vector<float> PSS_time_energy(N_PSS);
std::vector<std::vector<std::vector<cfloat>>> SSS_freq(N_PSS,
    std::vector<std::vector<cfloat>>(N_SSS,std::vector<cfloat>(SS_carriers)));

// definitions of the variables shared with main (declared extern in acquisition.h)
std::mutex worker_mutex;
std::condition_variable worker_cv;
std::condition_variable result_cv;
WorkerState worker_state = IDLE;
bool worker_should_exit = false;
std::vector<cfloat>* buffer_to_process = nullptr;
long long buffer_to_process_rx_ofs = 0;
double worker_freq_rx;
double worker_bandwidth;
double worker_samp_rate;
int worker_mu;
SSB_info_t worker_result;

// worker thread loop
void SSB_worker_thread() {
    while (true) {
        std::unique_lock<std::mutex> lock(worker_mutex);
        
        // Wait for work to do or an exit request
        worker_cv.wait(lock, [] { 
            return buffer_to_process != nullptr || worker_should_exit; 
        });
        
        if (worker_should_exit) {
            break;
        }
        
        // Take the parameters while holding the lock
        auto* buf = buffer_to_process;
        auto rx_ofs = buffer_to_process_rx_ofs;
        auto freq = worker_freq_rx;
        auto bw = worker_bandwidth;
        auto srate = worker_samp_rate;
        auto mu_val = worker_mu;
        
        // Signal that processing is in progress and reset for the next cycle
        worker_state = PROCESSING;
        buffer_to_process = nullptr;
        // Release the lock to process
        lock.unlock();
        // Notify that buffer_to_process has been taken (in case main is waiting)
        result_cv.notify_one();
        
        SSB_info_t result = SSB_search(buf, rx_ofs, freq, bw, srate, mu_val, false);

        // Reacquire the lock to store the result
        lock.lock();
        worker_result = result;        
        worker_state = RESULT_READY;
        lock.unlock();
        
        // Notify that the result is ready
        result_cv.notify_one();
    }
}

// the first sufficient power of 2, possibly reduced to 3/4
size_t ACQ_corrFFTsize;

// FFT of the three PSS over ACQ_corrFFTsize points
fftwf_complex *ACQ_corr_PSS_fftw[N_PSS];
cfloat *ACQ_corr_PSS[N_PSS];

// FFTW plans shared by the search threads: each thread executes them on its own buffers with
// fftwf_execute_dft (thread-safe, unlike the planner); all the buffers come from fftwf_alloc_complex,
// so they have the alignment of the buffers used for planning
static fftwf_plan ACQ_corr_z_to_freq;  // in-place FFT of z, ACQ_corrFFTsize points
static fftwf_plan ACQ_corr_xc_to_time; // in-place IFFT of the product with a PSS, ACQ_corrFFTsize points
static fftwf_plan SS_buf_to_freq;      // in-place FFT of one PSS/SSS/PBCH symbol, SS_carriers points
// fine timing of SSB_search (only run by the thread calling it): one symbol at the full sample rate.
// Planned once here: the FFTW planner is not thread-safe and the other threads plan while searches run
static fftwf_complex *ACQ_fine_buf_fftw;
static fftwf_plan ACQ_fine_to_freq;    // in-place FFT, SS_carriers*decim points

// scratch buffers of one search thread
struct SSBSearchCtx {
    std::vector<cfloat> r_rot; // r_window moved to the GSCN
    std::vector<cfloat> z;     // r_rot filtered and decimated
    // z is copied into this buffer, and the plan ACQ_corr_z_to_freq is executed for the FFT
    fftwf_complex *corr_z_fftw;
    cfloat *corr_z;
    // in these three buffers the products of ACQ_corr_PSS[PSS]
    // and corr_z are computed, and the plan ACQ_corr_xc_to_time is executed for the IFFT
    fftwf_complex *corr_z_xc_PSS_fftw[N_PSS];
    cfloat *corr_z_xc_PSS[N_PSS];
    // in this buffer we compute the FFTs of the symbols containing PSS/SSS/center of the PBCH
    fftwf_complex *SS_buf_fftw;
    cfloat *SS_buf;
};
static std::vector<SSBSearchCtx> ACQ_ctx; // [0]: the thread calling SSB_search, [1..]: the helper threads

// best SS block found at one GSCN
struct SSBCandidate {
    float SIR_est_dB = -1000.0;
    float Cfo_est = 0.0;
    size_t NID2 = -1;
    size_t NID1 = -1;
    size_t symSSB = -1;
    long long ofs = -1;
    long long coarse_ofs = -1; // timing of the PSS in r_raw
};

// one SSB_search call, shared by the search threads (each takes the next GSCN from next_GSCN)
struct SSBJob {
    const cfloat *r_raw;
    size_t r_size;
    long long r_ofs;
    double freq_rx;
    int mu;
    bool benchmark;
    float Ts, SCS, SSB_BW;
    size_t decim, N, L;
    std::vector<float> LP_rev; // decimating filter (LP_order_32.h) reversed, L*decim taps: one dot product per output
    size_t N_FFT, N_CP, CP_slack_corr, CP_slack_FFT;
    size_t iSSB_max;
    size_t symSSB_ofs[8];
    size_t symSSB_val[8];
    std::vector<int> GSCN_range;
    std::atomic<size_t> next_GSCN{0};
    std::vector<SSBCandidate> best; // per GSCN
    std::vector<std::string> logs;  // per GSCN, printed in GSCN order
};

// search thread pool: SSB_search publishes a job, then the helpers and the caller take its GSCNs
static std::vector<std::thread> ACQ_helpers;
static std::mutex ACQ_pool_mutex;
static std::condition_variable ACQ_pool_cv;      // new job or exit request
static std::condition_variable ACQ_pool_done_cv; // all the helpers are done with the job
static SSBJob *ACQ_pool_job = nullptr;
static uint64_t ACQ_pool_gen = 0;   // jobs published so far
static size_t ACQ_pool_running = 0; // helpers still working on the current job
static bool ACQ_pool_exit = false;

static void search_GSCN(SSBSearchCtx& ctx, SSBJob& J, size_t iG);

static void run_job(SSBSearchCtx& ctx, SSBJob& J)
{
    for(size_t iG; (iG = J.next_GSCN.fetch_add(1)) < J.GSCN_range.size(); )
        search_GSCN(ctx, J, iG);
}

static void SSB_helper_thread(size_t i_ctx)
{
    // --worker-cpus, those of the DCI workers, idle during ACQUIRE
    // (--rx-cpus are left to streaming and to the SSB worker)
    pin_self(worker_cpus);

    uint64_t gen = 0;
    while(true){
        SSBJob *job;
        {
            std::unique_lock<std::mutex> lock(ACQ_pool_mutex);
            ACQ_pool_cv.wait(lock, [&gen]{ return ACQ_pool_exit || ACQ_pool_gen != gen; });
            if(ACQ_pool_exit)
                return;
            gen = ACQ_pool_gen;
            job = ACQ_pool_job;
        }
        run_job(ACQ_ctx[i_ctx], *job);
        {
            std::lock_guard<std::mutex> lock(ACQ_pool_mutex);
            if(--ACQ_pool_running == 0)
                ACQ_pool_done_cv.notify_one();
        }
    }
}

void SS_init(void)
{
    const size_t SS_period = 127;
    
    // m-sequence
    std::vector<uint8_t> x(SS_period);
    x[0] = 0; x[1] = 1; x[2] = 1; x[3] = 0; x[4] = 1; x[5] = 1; x[6] = 1;
    for(size_t i = 7; i < SS_period; i++)
        x[i] = x[i-3] ^ x[i-7];
    
    size_t fft_size = SS_carriers;
    fftwf_complex *buf = fftwf_alloc_complex(fft_size);
    fftwf_plan plan = fftwf_plan_dft_1d(fft_size, buf, buf, FFTW_BACKWARD, FFTW_ESTIMATE);
    for(size_t NID2 = 0; NID2 < N_PSS; NID2++){
        PSS_freq[NID2][0] = cfloat(0.0,0.0); // -64
        for(size_t n = 0; n < SS_period; n++) // -63:63
            PSS_freq[NID2][1+n] = cfloat(1-2*x[(n+NID2*43) % SS_period],0);
        
        // n=[0,1:127]  => [-64,-63:63]
        for(size_t n = 0; n < 128; n++){
            buf[n][0] = PSS_freq[NID2][(n+64)&0x7F].real();
            buf[n][1] = PSS_freq[NID2][(n+64)&0x7F].imag();
        }
        fftwf_execute(plan);
        
        float temp = 0.0;
        for(size_t n = 0; n < PSS_len; n++){
            PSS_time[NID2][n] = cfloat(buf[(n-9)&0x7F][0],buf[(n-9)&0x7F][1]);
            temp += std::norm(PSS_time[NID2][n]);
        }
        PSS_time_energy[NID2] = temp;
    }
    fftwf_destroy_plan(plan);
    fftwf_free(buf);
    
    std::vector<uint8_t> x0(SS_period), x1(SS_period);
    x0[0] = 1; x1[0] = 1;
    x0[1] = 0; x1[1] = 0;
    x0[2] = 0; x1[2] = 0;
    x0[3] = 0; x1[3] = 0;
    x0[4] = 0; x1[4] = 0;
    x0[5] = 0; x1[5] = 0;
    x0[6] = 0; x1[6] = 0;
    for(size_t i = 7; i < SS_period; i++){
        x0[i] = x0[i-3] ^ x0[i-7];
        x1[i] = x1[i-6] ^ x1[i-7];
    }
    for(size_t NID2 = 0; NID2 < N_PSS; NID2++)
    for(size_t NID1 = 0; NID1 < N_SSS; NID1++){
        SSS_freq[NID2][NID1][0] = cfloat(0.0,0.0);
        for(size_t n = 0; n < SS_period; n++)
            SSS_freq[NID2][NID1][1+n] = cfloat(
                (1-2*x0[(n+(15*(NID1/112)+5*NID2))%SS_period])*(1-2*x1[(n+(NID1%112))%SS_period]),0);
    }
}

void ACQ_corr_init(size_t ACQ_window_size, size_t decim, int n_threads)
{
    size_t required_size = ACQ_window_size + LP_order-1 + PSS_len-1;
    ACQ_corrFFTsize = 1024;
    while(ACQ_corrFFTsize < required_size)
        ACQ_corrFFTsize *= 2;    
    if(ACQ_corrFFTsize*3/4 >= required_size)
        ACQ_corrFFTsize = ACQ_corrFFTsize/4*3;
    std::cout << "correlation FFT size: " << ACQ_corrFFTsize << ", SSB search threads: " << n_threads << std::endl;
    
    fftwf_complex *PSS_with_zero_padding = fftwf_alloc_complex(ACQ_corrFFTsize);
    for(size_t NID2 = 0; NID2 < N_PSS; NID2++){
        // time-reversed and conjugate
        for(size_t n = 0; n < PSS_len; n++){
            PSS_with_zero_padding[n][0] =  PSS_time[NID2][PSS_len-1-n].real();
            PSS_with_zero_padding[n][1] = -PSS_time[NID2][PSS_len-1-n].imag();
        }
        for(size_t n = PSS_len; n < ACQ_corrFFTsize; n++){
            PSS_with_zero_padding[n][0] = 0.0;
            PSS_with_zero_padding[n][1] = 0.0;
        }
        ACQ_corr_PSS_fftw[NID2] = fftwf_alloc_complex(ACQ_corrFFTsize);
        ACQ_corr_PSS[NID2] = reinterpret_cast<cfloat*>(ACQ_corr_PSS_fftw[NID2]);
        fftwf_plan plan = fftwf_plan_dft_1d(ACQ_corrFFTsize,
            PSS_with_zero_padding, ACQ_corr_PSS_fftw[NID2], FFTW_FORWARD, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);
    }
    fftwf_free(PSS_with_zero_padding);
    
    ACQ_ctx.resize(n_threads);
    for(auto& ctx : ACQ_ctx){
        ctx.corr_z_fftw = fftwf_alloc_complex(ACQ_corrFFTsize);
        ctx.corr_z = reinterpret_cast<cfloat*>(ctx.corr_z_fftw);
        for(size_t NID2 = 0; NID2 < N_PSS; NID2++){
            ctx.corr_z_xc_PSS_fftw[NID2] = fftwf_alloc_complex(ACQ_corrFFTsize);
            ctx.corr_z_xc_PSS[NID2] = reinterpret_cast<cfloat*>(ctx.corr_z_xc_PSS_fftw[NID2]);
        }
        ctx.SS_buf_fftw = fftwf_alloc_complex(SS_carriers);
        ctx.SS_buf = reinterpret_cast<cfloat*>(ctx.SS_buf_fftw);
    }
    ACQ_corr_z_to_freq = fftwf_plan_dft_1d(ACQ_corrFFTsize,
         ACQ_ctx[0].corr_z_fftw, ACQ_ctx[0].corr_z_fftw, FFTW_FORWARD, FFTW_MEASURE);
    ACQ_corr_xc_to_time = fftwf_plan_dft_1d(ACQ_corrFFTsize,
        ACQ_ctx[0].corr_z_xc_PSS_fftw[0], ACQ_ctx[0].corr_z_xc_PSS_fftw[0], FFTW_BACKWARD, FFTW_MEASURE);
    SS_buf_to_freq = fftwf_plan_dft_1d(SS_carriers, ACQ_ctx[0].SS_buf_fftw, ACQ_ctx[0].SS_buf_fftw, FFTW_FORWARD, FFTW_MEASURE);
    ACQ_fine_buf_fftw = fftwf_alloc_complex(SS_carriers*decim);
    ACQ_fine_to_freq = fftwf_plan_dft_1d(SS_carriers*decim, ACQ_fine_buf_fftw, ACQ_fine_buf_fftw, FFTW_FORWARD, FFTW_ESTIMATE);

    for(int i = 1; i < n_threads; i++)
        ACQ_helpers.emplace_back(SSB_helper_thread, (size_t)i);
}

void ACQ_corr_done(void)
{
    {
        std::lock_guard<std::mutex> lock(ACQ_pool_mutex);
        ACQ_pool_exit = true;
    }
    ACQ_pool_cv.notify_all();
    for(auto& t : ACQ_helpers)
        t.join();
    ACQ_helpers.clear();

    fftwf_destroy_plan(ACQ_fine_to_freq);
    fftwf_free(ACQ_fine_buf_fftw);
    fftwf_destroy_plan(SS_buf_to_freq);
    fftwf_destroy_plan(ACQ_corr_xc_to_time);
    fftwf_destroy_plan(ACQ_corr_z_to_freq);
    for(auto& ctx : ACQ_ctx){
        fftwf_free(ctx.SS_buf_fftw);
        for(size_t NID2 = 0; NID2 < N_PSS; NID2++)
            fftwf_free(ctx.corr_z_xc_PSS_fftw[NID2]);
        fftwf_free(ctx.corr_z_fftw);
    }
    ACQ_ctx.clear();
    for(size_t NID2 = 0; NID2 < N_PSS; NID2++)
        fftwf_free(ACQ_corr_PSS_fftw[NID2]);
}

// searches the SS block at GSCN J.GSCN_range[iG]: result in J.best[iG], log in J.logs[iG]
static void search_GSCN(SSBSearchCtx& ctx, SSBJob& J, size_t iG)
{
    const cfloat *r_raw = J.r_raw;
    const long long r_ofs = J.r_ofs;
    const double freq_rx = J.freq_rx;
    const int mu = J.mu;
    const bool benchmark = J.benchmark;
    const float Ts = J.Ts, SCS = J.SCS, SSB_BW = J.SSB_BW;
    const size_t decim = J.decim, N = J.N, L = J.L;
    const size_t N_FFT = J.N_FFT, N_CP = J.N_CP, CP_slack_corr = J.CP_slack_corr, CP_slack_FFT = J.CP_slack_FFT;
    const size_t iSSB_max = J.iSSB_max;
    const size_t *symSSB_ofs = J.symSSB_ofs;
    const size_t *symSSB_val = J.symSSB_val;
    const int GSCN = J.GSCN_range[iG];
    SSBCandidate& best = J.best[iG];
    std::ostringstream log;

    std::vector<cfloat>& r_rot = ctx.r_rot;
    r_rot.resize(J.r_size);
    cfloat *ACQ_corr_z = ctx.corr_z;
    cfloat **ACQ_corr_z_xc_PSS = ctx.corr_z_xc_PSS;
    cfloat *SS_buf = ctx.SS_buf;
    bool PSS_found_at_GSCN;
    double f_SSB = GSCN_to_freq(GSCN) - SCS; // the center of 56:182 is 119, not 120
    PSS_found_at_GSCN = false;
    
    // TIME-DOMAIN PROCESSING
    //
    // multiply r[i] by the complex exponential that
    // brings the SSB to baseband (subcarriers -63 to 63
    // of a 128-point FFT)
    cfloat wi = std::exp(cfloat(0.0f,2*M_PI*(freq_rx-f_SSB)*Ts));
    cfloat w(1.0,0.0);
    for(size_t i = 0; i < N*decim; i++){
        r_rot[i] = r_raw[i] * w;
        w *= wi;
    }
    
    // decimating low-pass filter (polyphase with decim phases)
    //
    std::vector<cfloat>& z = ctx.z;
    z.assign(N+L-1, cfloat(0.0,0.0));
    if(decim == 1)
        z = r_rot;
    else{
        // z[m] = sum_k LP[k]*r_rot[m*decim+decim-1-k] over the samples in r_rot, i.e. a dot product of the
        // reversed filter with K = L*decim consecutive samples starting at (m+1)*decim-K
        const long K = (long)(L*decim), R = (long)(N*decim);
        const float *h = J.LP_rev.data();
        const float *r = reinterpret_cast<const float*>(r_rot.data());
        for(size_t m = 0; m < N+L-1; m++){
            const long base = (long)(m+1)*(long)decim - K;
            const long t0 = std::max(0L, -base), t1 = std::min(K, R - base);
            float re = 0.0f, im = 0.0f;
            #pragma omp simd reduction(+:re,im)
            for(long t = t0; t < t1; t++){
                re += h[t]*r[2*(base+t)];
                im += h[t]*r[2*(base+t)+1];
            }
            z[m] = cfloat(re, im);
        }
    }
    //std::cout << "N: " << N << std::endl;
    //save_fc32_debug(z, "z");

    // correlate z with the three PSS_time (128 samples):
    // - ACQ_corrFFTsize-point FFT of z, with zero padding
    // - multiplication by the conjugate FFTs of PSS_time[0:2] (precomputed)
    // - ACQ_corrFFTsize-point IFFT, normalization by the energy of the
    //   sequence (sliding window)
    // - search for the maximum of the normalized correlation:
    //     store GSCN, NID2, and delay of the correlation peak (in the decimated signal)

    for(size_t n = 0; n < N+L-1; n++)
        ACQ_corr_z[n] = z[n];
    for(size_t n = N+L-1; n < ACQ_corrFFTsize; n++)
        ACQ_corr_z[n] = cfloat(0.0,0.0);
    fftwf_execute_dft(ACQ_corr_z_to_freq, ctx.corr_z_fftw, ctx.corr_z_fftw);
    for(size_t NID2 = 0; NID2 < N_PSS; NID2++){
    
        for(size_t n = 0; n < ACQ_corrFFTsize; n++)
            ACQ_corr_z_xc_PSS[NID2][n] = ACQ_corr_z[n]*ACQ_corr_PSS[NID2][n];
        fftwf_execute_dft(ACQ_corr_xc_to_time, ctx.corr_z_xc_PSS_fftw[NID2], ctx.corr_z_xc_PSS_fftw[NID2]);
        
        int max_pos = -1;
        float max_val = 0.0; // not normalized!
        float z_energy_max_val = 0.0;
        float z_energy = 0.0;
        for(size_t n = 0; n < N+L-1+PSS_len-1-4*(N_FFT+N_CP); n++){
            float val = std::norm(ACQ_corr_z_xc_PSS[NID2][n]);
            z_energy += std::norm(z[n]) * (ACQ_corrFFTsize*ACQ_corrFFTsize) * PSS_time_energy[NID2];
            if(n >= PSS_len)
                z_energy -= std::norm(z[n-PSS_len]) * (ACQ_corrFFTsize*ACQ_corrFFTsize) * PSS_time_energy[NID2];
            // n >= PSS_len-1: the whole PSS with its CP is in the window, so max_pos >= N_CP and the
            // symbols read below start inside z
            if(n >= PSS_len-1 && val/z_energy > PSS_threshold && val>max_val){
                max_val = val;
                z_energy_max_val = z_energy;
                max_pos = n - (PSS_len-9 - 1); // points to the start of the PSS symbol, after the CP
            }
        }

        // this threshold is set so that it is exceeded even when we receive several overlapping synchronous PSS
        if((z_energy_max_val>0.0 && max_val/z_energy_max_val > PSS_threshold) || benchmark){
            if(benchmark){
                max_pos = N_CP;
            }
            if(!PSS_found_at_GSCN){
                log << "GSCN #"<<GSCN <<"("<<(f_SSB+SCS)/1e6<<"MHz):" << std::endl;
                PSS_found_at_GSCN = true;
            }
            
            // FREQUENCY DOMAIN PROCESSING
            // - CFO estimate based on the cyclic prefix
            // - compensate the estimated CFO on PSS and SSS
            // - channel estimate with the PSS: H_est_PSS computing the 128-point FFT at
            //   the correlation peak
            // - SSS identification by searching the maximum-likelihood SSS sequence
            //   given H_est_PSS
            // - fine CFO estimate comparing the phase of H_est_SSS with that of H_est_PSS
            // - extraction of the PBCH DM-RS (visible subcarriers), application of the
            //   phase-based method to find iSSB and the frame start offset

            // estimate the CFO from the autocorrelation of the PBCH cyclic prefix (4 symbols)
            float Cfo_est_num = 0.0, Cfo_est_den = 0.0;
            for(size_t sym = 0; sym < 4; sym++){
                size_t ofs_sym = max_pos-N_CP+sym*(N_FFT+N_CP);
                cfloat p_num = cfloat(0.0,0.0);
                float p_den1 = 0.0;
                float p_den2 = 0.0;
                for(size_t n = CP_slack_corr; n < N_CP-CP_slack_corr; n++){
                    p_num += z[ofs_sym+n+N_FFT] * std::conj(z[ofs_sym+n]);
                    p_den1 += std::norm(z[ofs_sym+n]);
                    p_den2 += std::norm(z[ofs_sym+N_FFT+n]);
                }
                float CP_corr_sym = std::abs(p_num)/sqrtf(p_den1*p_den2);
                float Cfo_est_sym = std::arg(p_num)/(2.*M_PI*N_FFT/(float)SSB_BW);
                Cfo_est_num += CP_corr_sym*Cfo_est_sym;
                Cfo_est_den += CP_corr_sym;
            }
            float Cfo_est = Cfo_est_num/Cfo_est_den;

            log<<"  PSS #"<<NID2<<" max_corr: "<<max_val/z_energy_max_val<< "@"<<(r_ofs/decim+max_pos-N_CP)%(19200<<mu)
               <<" Cfo:"<<Cfo_est<<"Hz"<<std::endl;

            // symbol s of the SS block (0: PSS) with the CFO compensated, as subcarriers -64..63 of the SS block.
            // Only the relative phases across the subcarriers of a symbol matter below (every use is an absolute
            // value or a product of two symbols of the same block), so the mixer starts from phase 0 in every
            // symbol and the common phase of the CP_slack_FFT early start isn't removed
            const cfloat wi = std::exp(cfloat(0.0f,-2.0*M_PI*Cfo_est/(float)SSB_BW));
            auto SSB_symbol = [&](size_t s, std::vector<cfloat> &Yk){
                const size_t ofs = max_pos + s*(N_FFT+N_CP) - CP_slack_FFT;
                cfloat w(1.0f, 0.0f);
                for(size_t n = 0; n < N_FFT; n++){
                    SS_buf[n] = w*z[ofs+n];
                    w *= wi;
                }
                fftwf_execute_dft(SS_buf_to_freq, ctx.SS_buf_fftw, ctx.SS_buf_fftw);
                Yk.resize(N_FFT);
                for(size_t k = 0; k < N_FFT; k++)
                    Yk[k] = SS_buf[(k+64)&0x7F];
            };

            // H_est_PSS: fft(z[max_pos+(0:N_FFT-1)])./PSS_freq[NID2]
            std::vector<cfloat> H_est_PSS;
            SSB_symbol(0, H_est_PSS);
            float H_est_PSS_energy = 0.0;
            for(size_t k = 0; k < N_FFT; k++){ // the first sample contains noise
                H_est_PSS[k] *= PSS_freq[NID2][k];
                H_est_PSS_energy += std::norm(H_est_PSS[k]);
            }

            std::vector<cfloat> Y_SSS;
            SSB_symbol(2, Y_SSS);
            float Y_SSS_energy = 0.0;
            for(size_t k = 0; k < N_FFT; k++)
                Y_SSS_energy += std::norm(Y_SSS[k]);
            
            // compute the correlation between Y_SSS./SSS_freq[NID2][NID1] and H_est_PSS, for each NID1
            std::vector<float> SSS_corr(N_SSS);
            for (size_t NID1 = 0; NID1 < N_SSS; NID1++){
                cfloat tmp(0.0,0.0);
                for(size_t k = 0; k < N_FFT; k++)
                    tmp += Y_SSS[k]*SSS_freq[NID2][NID1][k] * std::conj(H_est_PSS[k]);
                SSS_corr[NID1] = std::abs(tmp)/sqrtf(H_est_PSS_energy*Y_SSS_energy);
            }
            //log << "SSS_corr=[";
            //for(size_t n = 0; n < N_SSS; n++)
            //    log << SSS_corr[n] << ((n<N_SSS-1)?",":"];");
            //log << std::endl;
            
            float NID1_max_corr_val = 0.0;
            size_t NID1_max_corr_pos = -1;
            for (size_t NID1 = 0; NID1 < N_SSS; NID1++){
                if(SSS_corr[NID1] > NID1_max_corr_val){
                    NID1_max_corr_pos = NID1;
                    NID1_max_corr_val = SSS_corr[NID1];
                }
            }
            if(NID1_max_corr_val > SSS_threshold || benchmark){
                if(benchmark){
                    NID1_max_corr_pos = 0;
                }
                log << "    SSS";
                for(size_t NID1 = 0; NID1 < N_SSS; NID1++){
                    if(SSS_corr[NID1] > SSS_threshold)
                        log << " #" << NID1<<"("<<SSS_corr[NID1]<<")";
                    if(NID1 == NID1_max_corr_pos)
                        log << "*";
                }
                log << std::endl;
                
                size_t NID1 = NID1_max_corr_pos;
                size_t cellid = 3*NID1+NID2;
                
                // H_est_SSS
                std::vector<cfloat> H_est_SSS(N_FFT);
                float H_est_SSS_energy = 0.0;
                cfloat H_est_SS_corr(0.0,0.0);
                for(size_t k = 0; k < N_FFT; k++){
                    H_est_SSS[k] = Y_SSS[k] * SSS_freq[NID2][NID1][k];
                    H_est_SS_corr += H_est_SSS[k] * std::conj(H_est_PSS[k]);
                    H_est_SSS_energy += std::norm(H_est_SSS[k]);
                }
                H_est_SS_corr /= sqrtf(H_est_PSS_energy*H_est_SSS_energy);
                float SIR_est_dB = 10.0*log10(std::abs(H_est_SS_corr)/(1.0-std::abs(H_est_SS_corr)));
                
                if(std::abs(H_est_SS_corr) > H_est_SS_corr_thres || benchmark){
                    //save_fc32_debug(H_est_PSS,"H_est_PSS");
                    //save_fc32_debug(H_est_SSS,"H_est_SSS");
                    
                    // the PBCH DM-RS subcarriers are those with index%4 == cellid%4
                    // the SS subcarriers are 56:182
                    // we see -64 to 63 relative to subcarrier 119, i.e. 55:182
                    // m_min = ceil((55-v)/4) ----> v+4*m_min>=55
                    // m_max = floor((182-v)/4) --> v+4*m_max<=182
                    size_t v = cellid%4;
                    size_t m_min = ceil((55.0-v)/4.0);
                    size_t m_max = floor((182.0-v)/4.0);
                    // positions of the DM-RS subcarriers in the FFT after fftshift: v+4*m-55
                    
                    // DM-RS sequences (only the segments of 32 that fall between subcarriers 55 and 182)
                    std::vector<uint8_t> c(120);
                    std::vector<std::vector<cfloat>> rSSB(8,std::vector<cfloat>(32));
                    
                    for (size_t iSSB = 0; iSSB < iSSB_max; iSSB++){
                        uint32_t c_init_SSB = PBCH_DMRS_c_init(cellid, iSSB);
                        scrambling(c.data(),120,c_init_SSB);
                        for(size_t m = m_min; m <= m_max; m++)
                            rSSB[iSSB][m-m_min] = cfloat(1.0-2.0*c[2*m],1.0-2.0*c[2*m+1]);
                    }
                    
                    // the DM-RS of the first PBCH symbol (that of the third one was weighted 0 and is not computed)
                    std::vector<cfloat> Y_PBCH1;
                    SSB_symbol(1, Y_PBCH1);

                    // compute the correlation between Y_PBCH1[v+4*(m_min:m_max)] and rSSB[iSSB].*H_est_SSS[v+4*(m_min:m_max)]
                    // is H_est_SSS needed? interference is not being taken into account!
                    float Y_PBCH1_DMRS_energy = 0.0;
                    float H_est_SSS_DMRS_energy = 0.0;
                    for(size_t m = m_min; m <= m_max; m++){
                        Y_PBCH1_DMRS_energy += std::norm(Y_PBCH1[v+4*m-55]);
                        H_est_SSS_DMRS_energy += 2*std::norm(H_est_SSS[v+4*m-55]); // |rSSB|^2==2!
                    }
                    
                    std::vector<float> iSSB_corr(iSSB_max);
                    for(size_t iSSB = 0; iSSB < iSSB_max; iSSB++){
                        cfloat tmp(0.0,0.0);
                        for(size_t m = m_min; m <= m_max; m++)
                            tmp += Y_PBCH1[v+4*m-55] * std::conj(H_est_SSS[v+4*m-55]*rSSB[iSSB][m-m_min]);
                        iSSB_corr[iSSB] = std::abs(tmp)/sqrtf(Y_PBCH1_DMRS_energy * H_est_SSS_DMRS_energy);
                    }
                    float iSSB_max_corr_val = -1.0;
                    size_t iSSB_max_corr_pos = -1;
                    for(size_t iSSB = 0; iSSB < iSSB_max; iSSB++){
                        if(iSSB_corr[iSSB] > iSSB_max_corr_val){
                            iSSB_max_corr_pos = iSSB;
                            iSSB_max_corr_val = iSSB_corr[iSSB];
                        }
                    }
                    
                    size_t iSSB = iSSB_max_corr_pos;
                    size_t rel_ofs = (r_ofs/decim + max_pos -N_CP - symSSB_ofs[iSSB]) % (19200<<mu);
                    long long ofs = r_ofs + decim * (max_pos - N_CP - symSSB_ofs[iSSB] - L/2);
                    log << "      cellid: "<<cellid<< " SSB SIR: "<<SIR_est_dB<< "dB"<<
                           " iSSB: "<<iSSB << "("<<iSSB_max_corr_val<<")@ofs: "<<rel_ofs<<std::endl;
                    
                    //std::cout << "iSSB_corr: [";
                    //for(size_t i = 0; i < iSSB_max; i++){
                    //    std::cout << iSSB_corr[i];
                    //    if (i == iSSB_max_corr_pos)
                    //        std::cout << "*";
                    //    if (i == iSSB_max-1)
                    //        std::cout << "]" << std::endl;
                    //    else
                    //        std::cout << " ";
                    //}
                    
                    if(SIR_est_dB > best.SIR_est_dB){
                        best.coarse_ofs = (long long)decim * ((long long)max_pos - (long long)(L/2) - (long long)N_CP);
                        best.ofs = ofs;
                        best.NID2 = NID2;
                        best.NID1 = NID1;
                        best.SIR_est_dB = SIR_est_dB;
                        best.Cfo_est = Cfo_est;
                        best.symSSB = symSSB_val[iSSB];
                    }
                }
            }
        }
    }
    J.logs[iG] = log.str();
}

SSB_info_t SSB_search(const std::vector<cfloat> *r_window, long long r_ofs,
                double freq_rx, double bandwidth, double samp_rate, int mu, bool benchmark)
{
    const cfloat *r_raw = r_window->data();
    const float Ts = 1.0/samp_rate;
    float SCS = 15000<<mu;
    float SSB_BW = SCS*128;
    size_t decim = samp_rate/SSB_BW;
    //std::cout << "decim_rate: " << decim << std::endl;
    size_t N = r_window->size()/decim; // input length N*decim
    const float *LP;
    switch(decim){
        case   2: LP =   LP2; break; case   3: LP =   LP3; break;
        case   4: LP =   LP4; break; case   6: LP =   LP6; break;
        case   8: LP =   LP8; break; case  12: LP =  LP12; break;
        case  16: LP =  LP16; break; case  24: LP =  LP24; break;
        case  32: LP =  LP32; break; case  48: LP =  LP48; break;
        case  64: LP =  LP64; break; case  96: LP =  LP96; break;
        case 128: LP = LP128; break; default:  LP =  NULL; break;
    }
    size_t L = (decim==1) ? 1 : LP_order;

    const size_t N_FFT = SS_carriers;
    const size_t N_CP = (N_FFT*9)/128;
    const size_t N_CP_07 = (N_FFT*(9+(1<<mu)))/128;
    const size_t CP_slack_corr = N_CP/9;
    const size_t CP_slack_FFT = 4*N_CP/9;
    //std::cout << "N_FFT: "<<N_FFT << " N_CP: "<<N_CP << " N_CP_07: "<<N_CP_07 << " CP_slack: "<<CP_slack << std::endl;
    size_t iSSB_max = -1;
    std::vector<size_t> symSSB_ofs(8);
    std::vector<size_t> symSSB_val(8);
    switch(mu){
    case 0: // one N_CP_07 every 7 symbols
        // mu=0 n28 Case A: {2,8}+14*n
        iSSB_max = 8;
        symSSB_ofs[0] =     N_FFT+N_CP_07 +   (N_FFT+N_CP);symSSB_val[0]= 2; //  2 Nhf=0
        symSSB_ofs[1] =  2*(N_FFT+N_CP_07)+ 6*(N_FFT+N_CP);symSSB_val[1]= 8; //  8
        symSSB_ofs[2] =  3*(N_FFT+N_CP_07)+13*(N_FFT+N_CP);symSSB_val[2]=16; // 16
        symSSB_ofs[3] =  4*(N_FFT+N_CP_07)+18*(N_FFT+N_CP);symSSB_val[3]=22; // 22
        symSSB_ofs[4] = 11*(N_FFT+N_CP_07)+61*(N_FFT+N_CP);symSSB_val[4]=72; //  2 Nhf=1
        symSSB_ofs[5] = 12*(N_FFT+N_CP_07)+66*(N_FFT+N_CP);symSSB_val[5]=78; //  8
        symSSB_ofs[6] = 13*(N_FFT+N_CP_07)+73*(N_FFT+N_CP);symSSB_val[6]=86; // 16
        symSSB_ofs[7] = 14*(N_FFT+N_CP_07)+78*(N_FFT+N_CP);symSSB_val[7]=92; // 22
        break;
    case 1: // one N_CP_07 every 14 symbols
        // mu=1 n77 Case C: {2,8}+14*n
        iSSB_max = 8;
        symSSB_ofs[0] =   (N_FFT+N_CP_07)+   (N_FFT+N_CP);symSSB_val[0]= 2; //  2
        symSSB_ofs[1] =   (N_FFT+N_CP_07)+ 7*(N_FFT+N_CP);symSSB_val[1]= 8; //  8
        symSSB_ofs[2] = 2*(N_FFT+N_CP_07)+14*(N_FFT+N_CP);symSSB_val[2]=16; // 16
        symSSB_ofs[3] = 2*(N_FFT+N_CP_07)+20*(N_FFT+N_CP);symSSB_val[3]=22; // 22
        symSSB_ofs[4] = 3*(N_FFT+N_CP_07)+27*(N_FFT+N_CP);symSSB_val[4]=30; // 30
        symSSB_ofs[5] = 3*(N_FFT+N_CP_07)+33*(N_FFT+N_CP);symSSB_val[5]=36; // 36
        symSSB_ofs[6] = 4*(N_FFT+N_CP_07)+40*(N_FFT+N_CP);symSSB_val[6]=44; // 44
        symSSB_ofs[7] = 4*(N_FFT+N_CP_07)+46*(N_FFT+N_CP);symSSB_val[7]=50; // 50
        break;
    }

    // the GSCNs are searched in parallel (search_GSCN) by this thread and by the helper threads
    SSBJob J;
    J.r_raw = r_raw;
    J.r_size = r_window->size();
    J.r_ofs = r_ofs;
    J.freq_rx = freq_rx;
    J.mu = mu;
    J.benchmark = benchmark;
    J.Ts = Ts; J.SCS = SCS; J.SSB_BW = SSB_BW;
    J.decim = decim; J.N = N; J.L = L;
    if(LP)
        J.LP_rev.assign(std::reverse_iterator<const float*>(LP + L*decim), std::reverse_iterator<const float*>(LP));
    J.N_FFT = N_FFT; J.N_CP = N_CP; J.CP_slack_corr = CP_slack_corr; J.CP_slack_FFT = CP_slack_FFT;
    J.iSSB_max = iSSB_max;
    std::copy(symSSB_ofs.begin(), symSSB_ofs.end(), J.symSSB_ofs);
    std::copy(symSSB_val.begin(), symSSB_val.end(), J.symSSB_val);
    freq_range_to_GSCN(freq_rx, bandwidth, mu, J.GSCN_range);
    J.best.resize(J.GSCN_range.size());
    J.logs.resize(J.GSCN_range.size());
    {
        std::lock_guard<std::mutex> lock(ACQ_pool_mutex);
        ACQ_pool_job = &J;
        ACQ_pool_gen++;
        ACQ_pool_running = ACQ_helpers.size();
    }
    ACQ_pool_cv.notify_all();
    run_job(ACQ_ctx[0], J);
    {
        std::unique_lock<std::mutex> lock(ACQ_pool_mutex);
        ACQ_pool_done_cv.wait(lock, []{ return ACQ_pool_running == 0; });
        ACQ_pool_job = nullptr;
    }

    // best SS block over all the GSCNs (ties: the first GSCN, as in a sequential search)
    float best_SIR_est_dB = -1000.0;
    float best_Cfo_est = 0.0;
    size_t best_GSCN = -1;
    size_t best_NID2 = -1;
    size_t best_NID1 = -1;
    size_t best_symSSB = -1;
    long long best_ofs = -1;
    long long coarse_ofs = -1; // timing of the best PSS in r_raw, -1 if no SS block was found
    std::string log;
    for(size_t iG = 0; iG < J.GSCN_range.size(); iG++){
        log += J.logs[iG];
        const SSBCandidate& c = J.best[iG];
        if(c.SIR_est_dB > best_SIR_est_dB){
            coarse_ofs = c.coarse_ofs;
            best_GSCN = J.GSCN_range[iG];
            best_ofs = c.ofs;
            best_NID2 = c.NID2;
            best_NID1 = c.NID1;
            best_SIR_est_dB = c.SIR_est_dB;
            best_Cfo_est = c.Cfo_est;
            best_symSSB = c.symSSB;
        }
    }
    
    // fine timing over +-decim samples; skipped when no SS block was found or the FFT windows do not fit in r_window
    const long long M_FFT = N_FFT*decim;
    const long long fine_first = coarse_ofs - (long long)decim - (long long)(CP_slack_FFT*decim);
    const long long fine_last = coarse_ofs + (long long)decim - (long long)(CP_slack_FFT*decim) + M_FFT; // exclusive
    if(!benchmark && best_SIR_est_dB > -1000.0 && fine_first >= 0 && fine_last <= (long long)r_window->size()){
        float max_SNR_dB_est = -1000.0;
        int fine_ofs_est = 0;
        cfloat *tmpbuf = reinterpret_cast<cfloat*>(ACQ_fine_buf_fftw);
        cfloat wi = std::exp(cfloat(0.0,-2.0*M_PI*best_Cfo_est*Ts));
        for(int fine_ofs = -((int)decim); fine_ofs <= (int)decim; fine_ofs++){
            cfloat w(1.0,0.0);
            for(size_t i = 0; i < N_FFT*decim; i++){
                tmpbuf[i] = r_raw[coarse_ofs + fine_ofs + i - CP_slack_FFT*decim]*w;
                w *= wi;
            }
            fftwf_execute(ACQ_fine_to_freq);
            const long long TRACK_SSB_k120 = (long long)round((GSCN_to_freq(best_GSCN)-freq_rx)/SCS);
            // FFT bin of SSB subcarrier n (0..239), wrapped into [0,M_FFT) also for negative frequencies
            auto bin = [&](long long n){ return (size_t)(((TRACK_SSB_k120-120+n)%M_FFT + M_FFT)%M_FFT); };
            float noise_est_num = 0.0;
            for(long long n = 0; n < 56; n++)
                noise_est_num += std::norm(tmpbuf[bin(n)]);
            for(long long n = 183; n < 240; n++)
                noise_est_num += std::norm(tmpbuf[bin(n)]);
            float power_est_num = 0.0;
            for(long long n = 57; n < 183; n++)
                power_est_num += std::norm(tmpbuf[bin(n)]);
            float SNR_dB_est = 10.0*log10((power_est_num/127.0)/(noise_est_num/113.0));
            //std::cout<<"freq_ofs: "<<freq_ofs<<" fine_ofs: "<<fine_ofs<<" SNR_dB_est:"<<SNR_dB_est<<std::endl;
            if(SNR_dB_est > max_SNR_dB_est){
                max_SNR_dB_est = SNR_dB_est;
                fine_ofs_est = fine_ofs;
            }
        }
        best_ofs += fine_ofs_est;
    }
    
    SSB_info_t SSB_info;
    if(best_SIR_est_dB > -1000.0 && !benchmark){
        if(verbose)
            std::cout << std::endl << log;
        SSB_info.valid = true;
        SSB_info.f_SSB = GSCN_to_freq(best_GSCN);
        SSB_info.ofs = best_ofs;
        SSB_info.NID2 = best_NID2;
        SSB_info.NID1 = best_NID1;
        SSB_info.SIR_est_dB = best_SIR_est_dB;
        SSB_info.Cfo_est = best_Cfo_est;
        SSB_info.symSSB = best_symSSB;
    }else{
        SSB_info.valid = false;
    }
    return SSB_info;
}
