#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread.hpp>
#include <boost/program_options.hpp>
#include <csignal>
#include <atomic>
#include <thread>
#include <mutex>
#include <vector>
#include <optional>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <limits>
#include <chrono>
#include <algorithm>
#include <pthread.h>
#include <fftw3.h>
#include "nr_common.h"
#include "LP_order_32.h"
#include "util.h"
#include "gscn.h"
#include "ofdm.h"
#include "pipeline.h"
#include "acquisition.h"
#include "dci_scan.h"
#include "rnti_log.h"
#include "summary.h"
#ifdef WITH_WEBUI
#include "webui.h"
#endif
#include "fallback_dci.h"
#include "pbch.h"
#include "DCIfactory.h"

// command line options
namespace po = boost::program_options;

// right now, the interface is dumb and you exit with CTRL-C (SIGINT)
static std::atomic<bool> stop_signal_called(false);
void sig_int_handler(int) { stop_signal_called = true; }

size_t BW_tab[13]=      {   5, 10, 15,  20,  25,  30,  40,  50,  60,  70,  80,  90, 100}; // MHz
size_t N_RB_tab[2][13]= {{ 25, 52, 79, 106, 133, 160, 216, 270,   0,   0,   0,   0,   0}, // 15kHz
                         { 11, 24, 38,  51,  65,  78, 106, 133, 162, 189, 217, 245, 273}};// 30kHz

// --------------------------------------------------------------------------------------------

int UHD_SAFE_MAIN(int argc, char *argv[]) {
    uhd::set_thread_priority_safe(0,true); // real-time prio 0 (0-99)

    const char* wisdom_file = ".fftw_wisdom";
    if (fftwf_import_wisdom_from_filename(wisdom_file))
        std::cout << "Loaded FFTW wisdom from " << wisdom_file << std::endl;
    else
        std::cout << "Couldn't load FFTW wisdom file " << wisdom_file << std::endl;

    std::string device_args(""); // serial=...
    double bandwidth = 30e6; // sampling rate
    double freq_rx = 1940e6; // initial receive frequency
    double SIR_threshold_dB = -1.0; // minimum SSB SIR for tracking (affected by multiple superposed PSS/SSS)
    double AGC_target = 0.15; // reasonable default with commercial networks
    int AGC_max_rounds = 3; // 3 rounds seems more than enough
    double AGC_Tskip = 0.16; // skip 160ms with a high DC transient
    double AGC_round_delay = 0.32; // 160ms observation interval
    int ACQ_max_fails = 10; // consecutive SSB searches without a cell before the AGC runs again
    std::string file_name(""); // IQ samples file (replaces the USRP)
    std::string file_format(""); // sc16 or fc32
    double file_rate = 0.0; // sample rate of the file (0: same rate used with a USRP)
    bool rnti_log = false; // per-RNTI DCI logs in logs/
    bool estimate_throughput = false; // DL throughput from the DCI resource allocations
    int fallback_BWP = 0; // BWP size [RB] of the fallback DCIs (formats 1_0 and 0_0)
    std::string webui_dest(""); // [host:]port of the dashboard bridge (webui/server.py)
    int fft_n_threads = 3; // OFDM demodulation threads (fft_worker)
    int acq_n_threads = 4; // SSB search threads (GSCNs searched in parallel)
    int mu_opt = -1; // numerology (-1: from the frequency)
    std::string dci_sizes_opt(""); // comma-separated DCI sizes to try (empty: all)
    // default for a 16 P-core (HT) + 16 E-core CPU: P-cores 0..15 (even: one per core, odd: HT siblings,
    // left unused), E-cores 16..31; streaming and FFT on 0,2,4,6, DCI processing on the others
    std::string rx_cpus_opt("0-6:2"); // RX thread, SSB worker, fft_workers
    std::string worker_cpus_opt("8-14:2,16-31"); // DCI workers, SSB search helpers
    // global DCIdecoder_n_workers = std::thread::hardware_concurrency() - 8; // how many matrix Y workers
    
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&device_args)->default_value(""), "USRP device args")
        ("file", po::value<std::string>(&file_name), "read IQ samples from this file instead of a USRP")
        ("format", po::value<std::string>(&file_format), "IQ file format: sc16 (interleaved int16) or fc32 (interleaved float) [default: file extension]")
        ("rate", po::value<double>(&file_rate), "IQ file sample rate in Hz [default: the rate chosen for a USRP]")
        ("bandwidth", po::value<double>(&bandwidth)->default_value(bandwidth), "Bandwidth in Hz")
        ("freq-rx", po::value<double>(&freq_rx)->default_value(freq_rx), "RX frequency in Hz (center frequency of the file with --file)")
        ("mu", po::value<int>(&mu_opt), "numerology of the SS block and of the received grid: 0 (15 kHz, SS block Case A) or 1 (30 kHz, Case C, e.g. n41) [default: 1 above 3 GHz, 0 below]")
        ("min_SSB_SIR", po::value<double>(&SIR_threshold_dB)->default_value(SIR_threshold_dB), "minimum SSB SIR for tracking [dB]")
        ("AGC_target", po::value<double>(&AGC_target)->default_value(AGC_target), "AGC target [<1]")
        ("AGC_rounds", po::value<int>(&AGC_max_rounds)->default_value(AGC_max_rounds), "AGC rounds [integer]")
        ("AGC_Tskip", po::value<double>(&AGC_Tskip)->default_value(AGC_Tskip), "AGC deaf time [s]")
        ("AGC_period", po::value<double>(&AGC_round_delay)->default_value(AGC_round_delay), "AGC period [s]")
        ("ACQ_max_fails", po::value<int>(&ACQ_max_fails)->default_value(ACQ_max_fails), "consecutive failed SSB searches before the AGC runs again [integer]")
        ("workers", po::value<int>(&DCIdecoder_n_workers)->default_value(DCIdecoder_n_workers), "processing workers [integer]")
        ("fft-threads", po::value<int>(&fft_n_threads)->default_value(fft_n_threads), "OFDM demodulation (FFT) threads [integer]")
        ("acq-threads", po::value<int>(&acq_n_threads)->default_value(acq_n_threads), "SSB search threads, each searching a part of the GSCNs [integer]")
        ("rx-cpus", po::value<std::string>(&rx_cpus_opt)->default_value(rx_cpus_opt), "CPUs of the RX thread, SSB search worker and FFT threads: list as taskset -c (e.g. 0,2,4,6 or 0-6:2) or \"all\"")
        ("worker-cpus", po::value<std::string>(&worker_cpus_opt)->default_value(worker_cpus_opt), "CPUs of the DCI workers and SSB search helpers: list as taskset -c (e.g. 8-14:2,16-31) or \"all\"")
        ("verbose", po::bool_switch(&verbose), "print every DCI, CSI-RS and tracked frame instead of a summary per second")
        ("rnti-log", po::bool_switch(&rnti_log), "write one DCI log per RNTI in logs-<date>-<time>/<RNTI>.log (closed after 10 s without DCIs)")
        ("estimate-throughput", po::bool_switch(&estimate_throughput), "estimate the DL/UL throughput in the summary (TBS of each DCI with --verbose); needs --fallback")
        ("fallback", po::value<int>(&fallback_BWP), "decode fallback DCIs (formats 1_0 and 0_0) for a BWP of this size [RB]")
        ("dci-sizes", po::value<std::string>(&dci_sizes_opt), "comma-separated DCI sizes [bits] to try in the blind search, e.g. 39,41,50 [default: all]")
        ("webui", po::value<std::string>(&webui_dest), "send the per-second summary as JSON over UDP to [host:]port (dashboard: webui/server.py)"
#ifndef WITH_WEBUI
                                                        "; not in this build: set WEBUI = 1 in config.mk"
#endif
        )
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return EXIT_SUCCESS;
    }

    if(DCIdecoder_n_workers < 1){
        std::cout << "--workers: at least one worker is needed" << std::endl;
        return EXIT_FAILURE;
    }
#ifndef WITH_WEBUI
    if(vm.count("webui")){
        std::cout << "--webui: built without the web UI, set WEBUI = 1 in config.mk and rebuild" << std::endl;
        return EXIT_FAILURE;
    }
#endif
    if(fft_n_threads < 1){
        std::cout << "--fft-threads: at least one thread is needed" << std::endl;
        return EXIT_FAILURE;
    }
    if(acq_n_threads < 1){
        std::cout << "--acq-threads: at least one thread is needed" << std::endl;
        return EXIT_FAILURE;
    }

    {
        std::string err;
        if(!parse_cpu_list(rx_cpus_opt, rx_cpus, err)){
            std::cout << "--rx-cpus: " << err << std::endl;
            return EXIT_FAILURE;
        }
        if(!parse_cpu_list(worker_cpus_opt, worker_cpus, err)){
            std::cout << "--worker-cpus: " << err << std::endl;
            return EXIT_FAILURE;
        }
        cpu_set_t both;
        CPU_AND(&both, &rx_cpus, &worker_cpus);
        if(CPU_COUNT(&both))
            std::cout << "Warning: --rx-cpus and --worker-cpus share CPUs " << cpu_list_string(both)
                      << ": the DCI workers will compete with the real-time RX thread" << std::endl;
        // the threads spawned from here on (codec init, RX thread, SSB worker) inherit rx_cpus;
        // the fft_workers, DCI workers and SSB helpers pin themselves
        pin_self(rx_cpus);
        if(verbose)
            std::cout << "RX CPUs: " << cpu_list_string(rx_cpus) << ", worker CPUs: " << cpu_list_string(worker_cpus) << std::endl;
    }

    if(vm.count("mu") && mu_opt != 0 && mu_opt != 1){
        std::cout << "--mu: only 0 (15 kHz) and 1 (30 kHz) are supported" << std::endl;
        return EXIT_FAILURE;
    }
    if(mu_opt == 0 && freq_rx > 3e9){
        // Case A above 3 GHz has 8 SS blocks per half frame, while mu=0 is handled with 4 (L_max, SSB_search)
        std::cout << "--mu 0: 15 kHz SS blocks above 3 GHz (Case A with 8 SS blocks) are not supported" << std::endl;
        return EXIT_FAILURE;
    }

    if(vm.count("fallback") && (fallback_BWP < 1 || fallback_BWP > 275)){
        std::cout << "--fallback: the BWP size must be in 1..275 RBs" << std::endl;
        return EXIT_FAILURE;
    }
    if(estimate_throughput && !vm.count("fallback")){
        std::cout << "--estimate-throughput needs the DCI format to decode: use --fallback BWPsize" << std::endl;
        return EXIT_FAILURE;
    }
    const int throughput_BWP = estimate_throughput ? fallback_BWP : 0;

    if(vm.count("dci-sizes")){
        std::istringstream ss(dci_sizes_opt);
        std::string tok;
        while(std::getline(ss, tok, ',')){
            size_t end = 0;
            long A = -1;
            try{ A = std::stol(tok, &end); } catch(const std::exception&){ end = 0; }
            if(end == 0 || tok.find_first_not_of(" \t", end) != std::string::npos || A < (long)A_min || A > (long)A_max[N_log2AL-1]){
                std::cout << "--dci-sizes: \"" << tok << "\" is not a DCI size in " << A_min << ".." << A_max[N_log2AL-1] << std::endl;
                return EXIT_FAILURE;
            }
            DCI_sizes.push_back((uint32_t)A);
        }
        std::sort(DCI_sizes.begin(), DCI_sizes.end());
        DCI_sizes.erase(std::unique(DCI_sizes.begin(), DCI_sizes.end()), DCI_sizes.end());
        if(DCI_sizes.empty()){
            std::cout << "--dci-sizes: no DCI size given" << std::endl;
            return EXIT_FAILURE;
        }
    }

    const bool from_file = vm.count("file");
    if(from_file){
        for(const char* opt : {"args", "AGC_target", "AGC_rounds", "AGC_Tskip", "AGC_period", "ACQ_max_fails"})
            if(!vm[opt].defaulted()){
                std::cout << "--" << opt << " is a USRP option and can't be used with --file" << std::endl;
                return EXIT_FAILURE;
            }
        if(file_format.empty()){
            size_t dot = file_name.rfind('.');
            if(dot != std::string::npos)
                file_format = file_name.substr(dot+1);
        }
        if(file_format != "sc16" && file_format != "fc32"){
            std::cout << "unknown IQ file format \"" << file_format << "\", use --format sc16 or --format fc32" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if(vm.count("format") || vm.count("rate")){
        std::cout << "--format and --rate can only be used with --file" << std::endl;
        return EXIT_FAILURE;
    }
    if(ACQ_max_fails < 1){
        std::cout << "--ACQ_max_fails must be at least 1" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream iq_file;
    if(from_file){
        iq_file.open(file_name, std::ios::binary);
        if(!iq_file){
            std::cout << "can't open " << file_name << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << "Reading " << file_format << " IQ samples from " << file_name << std::endl;
    }

    uhd::usrp::multi_usrp::sptr usrp;
    uhd::rx_streamer::sptr rx_stream;
    std::string model;
    double mcr = 0.0;
    if(!from_file){
    std::cout << "\nCreating USRP device..." << std::endl;
    usrp = uhd::usrp::multi_usrp::make(device_args);
    std::cout << "Using Device: " << usrp->get_pp_string() << std::endl;

    model = usrp->get_mboard_name();
    uhd::stream_args_t stream_args;
    if(model == "x410" || model == "n310" || model == "x310"){
        stream_args = uhd::stream_args_t("fc32", "sc16");
        mcr = usrp->get_master_clock_rate();
    }
    else if(model == "B200mini" || model == "B210"){
        stream_args = uhd::stream_args_t("fc32", "sc12");
        mcr = 0.0; // any rate is possible on the b200mini/B210
    }
    else
        std::cout <<"[WARNING] add this model to the stream_args setting routine: "<<model << std::endl << std::flush;
    //stream_args.args["spp"] = "1200";
    stream_args.channels = {0};
    rx_stream = usrp->get_rx_stream(stream_args);
    } // !from_file

    // supported:
    // - Case A (FR1 fc<3GHz SCS=15kHz mu=0), L_max = 4
    // - Case C (FR1 SCS=30kHz mu=1), L_max = 8: above 3 GHz, or below with --mu 1 (unpaired bands above 1.88 GHz, e.g. n41)
    // TODO: Case B (30 kHz, e.g. n5, n66) and the Case C bands with L_max = 4 (paired below 3 GHz) are not supported
    const int mu = mu_opt >= 0 ? mu_opt : (freq_rx>3e9 ? 1 : 0);
    
    int BW_id = -1;
    for(int i = 0; i < 13; i++)
        if(BW_tab[i] == std::round(bandwidth/1e6)){
            BW_id=i;
            break;
        }
    if(BW_id==-1 || N_RB_tab[mu][BW_id]==0){
        std::cout << "unsupported (bandwidth,frequency range) combination."<<std::endl;
        exit(EXIT_FAILURE);
    }
        
    const size_t N_RB = N_RB_tab[mu][BW_id];
    const size_t Nsc_RB = 12;
    const size_t N = N_RB*Nsc_RB;
    
    // NR constants
    const double Tf = 10e-3;
    const size_t SCS = 15000<<mu;
    const size_t samp_rate_unit = 128*SCS;
    size_t ACQ_decim;
    double samp_rate;
    if(from_file && file_rate > 0.0){
        // the file sample rate is given: it must be an integer multiple of 128*SCS covering the bandwidth
        ACQ_decim = std::round(file_rate/samp_rate_unit);
        samp_rate = samp_rate_unit*ACQ_decim;
        bool decim_ok = (ACQ_decim == 1);
        for(size_t i = 0; i < sizeof(LP_decim_rates)/sizeof(LP_decim_rates[0]); i++)
            if(LP_decim_rates[i] == ACQ_decim)
                decim_ok = true;
        if(std::abs(file_rate-samp_rate) > 1.0 || !decim_ok || samp_rate < bandwidth){
            std::cout << "unsupported file sample rate: it must be 128*SCS*decim >= bandwidth, with decim in LP_order_32.h" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(bandwidth <= samp_rate_unit){
        ACQ_decim = 1;
        samp_rate = samp_rate_unit;
    }
    else{
        bool ACQ_decim_found = false;
        for(size_t i = 0; i < sizeof(LP_decim_rates)/sizeof(LP_decim_rates[0]); i++)
            if(samp_rate_unit*LP_decim_rates[i] >= bandwidth){
                ACQ_decim = LP_decim_rates[i];
                samp_rate = samp_rate_unit*ACQ_decim;
                if (mcr == 0.0 || std::abs(1.0-samp_rate/mcr*std::round(mcr/samp_rate))<1e-6){
                    ACQ_decim_found = true;
                    break;
                }
            }
        if(!ACQ_decim_found){
            std::cout << "unsupported bandwidth (or master clock rate), increase the maximum decim rate in LP_order_32.h"<<std::endl;
            exit(EXIT_FAILURE);
        }
    }
    if(from_file)
        std::cout << "IQ file sample rate: " << samp_rate/1e6 << " MHz" << std::endl;
    if(model == "B200mini" || model == "B210"){
        mcr = std::floor(61.44e6/samp_rate)*samp_rate;
        usrp->set_master_clock_rate(mcr);
    }

    const double Ts = 1.0/samp_rate;
    const size_t samples_per_frame = samp_rate_unit/100*ACQ_decim;
    const size_t N_FFT = ACQ_decim*128;
    const size_t N_CP = ACQ_decim*9;
    const size_t N_CP_07 = ACQ_decim*(9+(1<<mu));
    const size_t CP_slack = ACQ_decim; // worked well in nr_demo
    const size_t Nsymb_frame = 140<<mu;

    // the receive frequency should be an integer number of subcarriers away from the
    // frequency of the GSCNs in the band. otherwise, print a warning
    // and set the receive frequency closest to the requested one that
    // satisfies the constraint
    {
        std::vector<int> GSCN_range;
        //std::cout<<"bandwidth:"<<bandwidth<<std::endl;
        freq_range_to_GSCN(freq_rx, bandwidth, mu, GSCN_range);
        if(GSCN_range.size() == 0){
            std::cout<<"[ERROR] No known GSCN's in the receive bandwidth, enable scan mode or change receive parameters."<<std::endl;
            return EXIT_FAILURE;
        }
        else{
            double f120 = GSCN_to_freq(GSCN_range[0]);
            double fSCS = (double)SCS;
            int k120 = round((f120-freq_rx)/fSCS);
            if(std::abs(freq_rx-(f120-k120*fSCS)) > 0.01*fSCS && from_file){
                // a recording can't be retuned: the residual offset is left to the CFO estimate
                std::cout<<"[WARNING] The file center frequency is incompatible with the GSCN frequency raster"<<std::endl;
            }
            else if(std::abs(freq_rx-(f120-k120*fSCS)) > 0.01*fSCS){
                std::cout<<"[WARNING] The selected receive frequency is incompatible with the GSCN frequency raster"<<std::endl;
                std::cout<<"          The closest GSCN-compatible frequency is "<<(f120-k120*fSCS)/1e6<<"MHz, overriding freq-rx."<<std::endl;
                freq_rx = f120-k120*fSCS;
            }
        }
    }
    // phase_comp depends on the final freq_rx
    std::vector<float> phase_comp(Nsymb_frame);
    std::vector<size_t> sym_ofs(Nsymb_frame);
    size_t tmp = 0;
    // std::cout << "N_CP:"<<N_CP<<" N_CP_07:"<<N_CP_07<<std::endl;
    for(size_t n = 0; n < Nsymb_frame; n++){
        tmp += (n%(7<<mu)) ? N_CP : N_CP_07;
        sym_ofs[n] = tmp;
        phase_comp[n] = 2.0*M_PI*std::remainder(freq_rx*Ts*tmp,1.0); // in double: freq_rx*Ts*tmp is up to ~1e8 cycles
        tmp += N_FFT;
    }
    //assert(tmp==samp_rate/100);

    // PBCH decoder (MIB, SFN, SS block index), checked before any thread is started
    PBCHDecoder pbch(mu);
    if(!pbch.selftest(false)){
        std::cout << "PBCH decoder selftest failed" << std::endl;
        return EXIT_FAILURE;
    }

    const size_t NY_MAX = samp_rate/100;
    const size_t ALIGN_BYTES = 64;
    const size_t POOL_INITIAL = 4;
    const size_t POOL_MAX = 16;
    //const size_t MAX_OUT_QUEUE = 4096;
    BufferPool pool(NY_MAX, ALIGN_BYTES, POOL_INITIAL, POOL_MAX);
    BoundedQueue<Job>    in_q;                // unbounded: the pool provides backpressure
    BoundedQueue<Result> out_q;
    // bounded: every FFTJob holds a whole frame of samples. With --file the RX thread waits for room,
    // live it drops the frame (RxStatus::frames_dropped) rather than stall the stream
    const size_t FFT_Q_MAX = 16;
    BoundedQueue<FFTJob> fft_q(FFT_Q_MAX);

    // DCI codec matrices, built once and shared by the workers; a worker that can't decode would silently drop its frames
    std::cout << "Initializing DCI decoding workers..." << std::flush;
    DCIfactory DCIcodec;
    {
        // the init threads inherit the affinity of this thread: give them the worker CPUs too, idle for now
        cpu_set_t init_cpus;
        CPU_OR(&init_cpus, &rx_cpus, &worker_cpus);
        pin_self(init_cpus);
        DCIcodec.init();
        pin_self(rx_cpus);
    }
    if(!DCIcodec.selftest(false)){
        std::cout << " DCI codec selftest failed" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << " done"<<std::endl;
    
    int grid_ofs = N_FFT/2-N/2; // offset of the grid start in the FFT

    fftwf_complex *Y_tmpbuf_fftw = fftwf_alloc_complex(N_FFT);
    cfloat *Y_tmpbuf = reinterpret_cast<cfloat*>(Y_tmpbuf_fftw);
    const size_t NY = Nsymb_frame*N_FFT;
    fftwf_complex *Y_fftw = fftwf_alloc_complex(NY);
    std::vector<fftwf_plan> Y_symb_to_freq(Nsymb_frame);
    for(size_t n = 0; n < Nsymb_frame; n++)
        Y_symb_to_freq[n] = fftwf_plan_dft_1d(N_FFT, Y_fftw+n*N_FFT, Y_tmpbuf_fftw, FFTW_FORWARD, FFTW_MEASURE);
    std::vector<bool> is_Y_valid(Nsymb_frame);
    cfloat *Y = reinterpret_cast<cfloat*>(Y_fftw);
    fftwf_complex *Y_CP_fftw = fftwf_alloc_complex(N_CP_07);
    cfloat *Y_CP = reinterpret_cast<cfloat*>(Y_CP_fftw);
    size_t Y_CP_ofs = 0;
    std::vector<cfloat> Y_CP_corr(Nsymb_frame);
    float noise_est = 0.0;
    
    std::signal(SIGINT, &sig_int_handler);

    // absolute offset of the first sample in rx_buf
    long long rx_ofs = 0;

    // initial RX gain and range of possible values
    uhd::gain_range_t rx_gain_range;
    double gain_rx = 0.0;

    std::vector<cfloat> rx_buf;
    std::vector<int16_t> file_buf; // sc16 file samples, before the conversion to fc32
    if(from_file){
        rx_buf.resize(16384);
        if(file_format == "sc16")
            file_buf.resize(2*rx_buf.size());
    }
    else{
    std::cout << "Setting RX Rate: " << samp_rate/1e6 << " MHz";
    usrp->set_rx_rate(samp_rate);
    double actual_samp_rate = usrp->get_rx_rate();
    std::cout << " (set to: " << actual_samp_rate/1e6 << " MHz)" << std::endl;
    if(fabs(actual_samp_rate-samp_rate)>1.0){
        std::cout << "The required sample rate is not possible with the default master clock rate." << std::endl;
        return EXIT_FAILURE;
    }

    if(model == "B200mini" || model == "B210"){
        std::cout << "Setting RX Bandwidth: " << bandwidth/1e6 << " MHz";
        usrp->set_rx_bandwidth(bandwidth);
        double actual_rx_bandwidth = usrp->get_rx_bandwidth();
        std::cout << " (set to: " << actual_rx_bandwidth/1e6 << " MHz)" << std::endl;
    }

    size_t max_samps_per_packet = rx_stream->get_max_num_samps();
    std::cout << "Max samples per packet: " << max_samps_per_packet << std::endl;
    rx_buf.resize(max_samps_per_packet);

    // after a tuning request, the state must always go to AGC
    std::cout << "Setting RX Freq: " << freq_rx/1e6 << " MHz";
    uhd::tune_request_t tune_request_rx(freq_rx);
    tune_request_rx.args = uhd::device_addr_t("mode_n=integer");
    usrp->set_rx_freq(tune_request_rx);
    while (!usrp->get_rx_sensor("lo_locked",0).to_bool()){
        struct timespec ts;
        ts.tv_sec = 0;
        ts.tv_nsec = 10000000; // 10ms
        nanosleep(&ts, NULL);
        std::cout<<".";
    }
    std::cout << " (set to: " << usrp->get_rx_freq()/1e6 << " MHz)" << std::endl;

    rx_gain_range = usrp->get_rx_gain_range();
    gain_rx = rx_gain_range.step()*round(rx_gain_range.stop()/(4*rx_gain_range.step()));
    usrp->set_rx_gain(gain_rx);
    } // !from_file

    // everything that can fail is done before the first thread is started: returning (or throwing) with
    // joinable threads would abort the program
    std::optional<RNTILogger> rnti_logger; // throws if logs/ can't be made
    if(rnti_log){
        try{
            rnti_logger.emplace("logs", samp_rate, mu, 10.0); // close files after 10 s of stream time without DCIs
        }
        catch(const std::exception& e){
            std::cout << e.what() << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << "RNTI logs in " << rnti_logger->dir() << "/" << std::endl;
    }
    WebUI* webui_ptr = nullptr; // the summary also goes to the web UI if set
#ifdef WITH_WEBUI
    std::optional<WebUI> webui;
    if(!webui_dest.empty()){
        std::ostringstream cfg; // static run configuration, repeated in every datagram
        std::string source = from_file ? file_name : model;
        for(char& ch : source) // keep the JSON string trivially valid
            if(ch == '"' || ch == '\\' || (unsigned char)ch < 0x20)
                ch = '_';
        cfg << "{\"source\":\"" << source << "\",\"file\":" << (from_file ? "true" : "false")
            << ",\"freq_rx_hz\":" << std::fixed << std::setprecision(0) << freq_rx
            << ",\"bandwidth_hz\":" << bandwidth << ",\"samp_rate_hz\":" << samp_rate
            << ",\"mu\":" << mu << ",\"scs_hz\":" << SCS << ",\"n_rb\":" << N_RB
            << ",\"fallback_bwp\":" << (vm.count("fallback") ? std::to_string(fallback_BWP) : "null")
            << ",\"throughput\":" << (throughput_BWP > 0 ? "true" : "false") << "}";
        try{
            webui.emplace(webui_dest, cfg.str());
        }
        catch(const std::exception& e){
            std::cout << e.what() << std::endl;
            return EXIT_FAILURE;
        }
        std::cout << "Web UI updates to UDP " << webui_dest << std::endl;
        webui_ptr = &*webui;
    }
#endif

    // DCI workers, each decoding with its own copy of DCIcodec
    std::vector<std::thread> workers;
    workers.reserve(DCIdecoder_n_workers);
    for (int i = 0; i < DCIdecoder_n_workers; ++i)
        workers.emplace_back(worker, std::ref(in_q), std::ref(out_q), N_FFT, std::cref(DCIcodec));

    // receiver states
    enum { AGC=0, ACQUIRE, TRACK } state;
    
    // AGC parameters
    float AGC_max_real = 0, AGC_min_real = 0, AGC_max_imag = 0, AGC_min_imag = 0;

    // AGC state variables
    state = AGC;
    bool AGC_complete = false;
    int AGC_rounds = AGC_max_rounds;
    long long rx_ofs_last_AGC_round = rx_ofs;

    // initialize the PSS/SSS sequences
    SS_init();

    // ACQuisition parameters
    double T_SSB = 2*Tf;                 // SSB periodicity (minimum 20ms)
    size_t ACQ_bufsize = std::ceil(1.5*T_SSB*samp_rate); // a bit more than the periodicity to avoid borders
    ACQ_corr_init(ACQ_bufsize/ACQ_decim, ACQ_decim, acq_n_threads);

    std::chrono::duration<double, std::milli> toc_tic;
    {
        auto tic = std::chrono::high_resolution_clock::now();
        std::vector<cfloat> tmp_buf(ACQ_bufsize);
        SSB_search(&tmp_buf,0LL,freq_rx,bandwidth,samp_rate, mu, true);
        auto toc = std::chrono::high_resolution_clock::now();
        toc_tic = toc-tic;
    }
    std::cout << "toc: " << toc_tic.count() << "ms" << std::endl;
    
    size_t ACQ_period = 1.25*toc_tic.count()/1000*samp_rate;
    if(ACQ_period < ACQ_bufsize)
        ACQ_period = ACQ_bufsize;
    size_t ACQ_samples_to_skip = ACQ_period - ACQ_bufsize;
    if(from_file) // no real-time constraint: RX waits for each SSB search instead of skipping samples
        ACQ_samples_to_skip = 0;
    std::cout << "ACQ period: "<<ACQ_period << " process window: "<<ACQ_bufsize << " (skip: "<<ACQ_samples_to_skip<<")"<<std::endl;

    // Initialize global parameters for the worker
    worker_freq_rx = freq_rx;
    worker_bandwidth = bandwidth;
    worker_samp_rate = samp_rate;
    worker_mu = mu;
    // Create the worker thread BEFORE starting the streaming
    std::thread ACQ_thread{ SSB_worker_thread };

    // ACQuisition buffer and state
    std::vector<cfloat> ACQ_buf1;
    ACQ_buf1.reserve(ACQ_bufsize);
    std::vector<cfloat> ACQ_buf2;
    ACQ_buf2.reserve(ACQ_bufsize);

    std::vector<cfloat>* ACQ_buf = &ACQ_buf1;
    bool ACQ_buf_is_buf1 = true;
    long long ACQ_buf_rx_ofs;
    long long last_ACQ_buf_rx_ofs = -1;
    long long next_ACQ_buf_rx_ofs = -1;
    
    size_t ACQ_samples_in_buf;
    size_t ACQ_skipped_samples;
    bool ACQ_is_skipping;

    int ACQ_fails = 0; // consecutive SSB searches without a cell

    // --file: the SFNs whose PDCCH was searched, as runs of consecutive frames shipped with a noise estimate
    // and a known SFN; printed at the end, so an evaluation can count the reference DCIs over the same frames
    struct SFNSegment { int first, last, frames; };
    std::vector<SFNSegment> SFN_segments;
    bool SFN_segment_open = false;    // the next frame may extend SFN_segments.back()
    int frames_before_SFN_known = 0;  // searched, but shipped with the local frame counter
    std::vector<cfloat>* ACQ_buf_in_search = nullptr; // last buffer handed to the SSB worker
    bool ACQ_drop_result = false; // the search in flight started before ACQ_start: its result is stale

    // (re)start the SSB search from an empty buffer
    auto ACQ_start = [&]() {
        ACQ_samples_in_buf = 0;
        ACQ_skipped_samples = -1;
        ACQ_buf_rx_ofs = -1;
        ACQ_is_skipping = false;
        last_ACQ_buf_rx_ofs = -1;
        {
            // the samples searched so far precede the restart (other gain, or an rx_ofs before an overflow):
            // a buffer not picked up yet is withdrawn, a result is dropped, now or when the search in flight ends
            std::lock_guard<std::mutex> lock(worker_mutex);
            buffer_to_process = nullptr;
            if(worker_state == RESULT_READY)
                worker_state = IDLE;
            ACQ_drop_result = (worker_state == PROCESSING);
        }
        // the buffer of the search in flight is left alone until it ends
        ACQ_buf_is_buf1 = !(ACQ_drop_result && ACQ_buf_in_search == &ACQ_buf1);
        ACQ_buf = ACQ_buf_is_buf1 ? &ACQ_buf1 : &ACQ_buf2;
        ACQ_buf->clear();
        next_ACQ_buf_rx_ofs = -1;
        ACQ_fails = 0;
        SFN_segment_open = false;
        state = ACQUIRE;
    };
    // run the AGC again, then ACQUIRE (after a tracking loss or a failed acquisition: the signal level
    // may have changed abruptly); with --file there is no gain to control, so ACQUIRE right away
    auto AGC_restart = [&]() {
        if(from_file){
            ACQ_start();
            return;
        }
        AGC_complete = false;
        AGC_rounds = AGC_max_rounds;
        AGC_max_real = AGC_min_real = AGC_max_imag = AGC_min_imag = 0;
        rx_ofs_last_AGC_round = rx_ofs;
        SFN_segment_open = false;
        state = AGC;
    };

    if(from_file){ // no gain to control, start acquiring right away
        AGC_complete = true;
        ACQ_start();
    }

    // TRACKING
    long long last_frame_rx_ofs = 0; // set when leaving the ACQ state
    long long next_frame_rx_ofs = 0; // expected start of the next frame
    size_t Y_ofs; // offset of the next sample to insert in the buffer
    size_t Y_cur_symb; // index of the next symbol to fill (sym_ofs[Y_cur_symb]-CP_slack+(0:N_FFT-1))
    bool is_in_symbol; // whether the next sample from the buffer goes into Y or not
    float Cfo_eD,Cfo_eF2,Cfo_eF,Cfo_rounds_since_last_update;
    float Sfo_eD,Sfo_eF2,Sfo_eF,Sfo_rounds_since_last_update;
    // the expected SS block is in the frame if both the correlation and the phase coherence between
    // H_est_SSS and H_est_PSS exceed this (noise alone gives about 1/sqrt(127) = 0.09 for both)
    const float SS_presence_thres = 0.3f;
    float Kpf=0.1;     // 0.2667
    float Kif=Kpf*Kpf/4.;
    float Cfo_eD_num=0.0,Cfo_eD_den=0.0;

    // used for tracking    
    uint16_t LFN;
    size_t TRACK_SSB_sym = -1;
    int32_t TRACK_SSB_k120 = 0;
    int32_t TRACK_NID2 = -1;
    int32_t TRACK_NID1 = -1;
    std::vector<cfloat> H_est_PSS(SS_carriers);
    std::vector<cfloat> H_est_SSS(SS_carriers);
    fftwf_complex *TRACK_H_ref_fftw = fftwf_alloc_complex(SS_carriers);
    cfloat *TRACK_H_ref = reinterpret_cast<cfloat*>(TRACK_H_ref_fftw);
    bool TRACK_H_ref_is_valid = false;
    fftwf_complex *TRACK_H_tmp_fftw = fftwf_alloc_complex(4*SS_carriers);
    cfloat *TRACK_H_tmp = reinterpret_cast<cfloat*>(TRACK_H_tmp_fftw);
    fftwf_plan TRACK_H_tmp_to_time = fftwf_plan_dft_1d(4*SS_carriers, TRACK_H_tmp_fftw, TRACK_H_tmp_fftw, FFTW_BACKWARD, FFTW_MEASURE);
    int32_t TRACK_samples_to_skip = 0;
    int32_t TRACK_samples_to_skip_cumsum = 0;
    // PBCH: once the MIB is known, LFN is the SFN and the PBCH is checked again every PBCH_check_period SS blocks
    const int PBCH_check_period = 100;
    bool TRACK_MIB_valid = false;
    MIB TRACK_MIB{};
    int TRACK_PBCH_countdown = 0;
    CORESET0Grid TRACK_cs0; // CORESET#0 for the DCI workers, once the MIB is known
    std::vector<cfloat> SSB_grid(SSB_symbols*SSB_carriers);
    // the PBCH puts the SS block in another symbol of the frame (SS block index, half frame): the frames
    // restart from the frame of that SS block, TRACK_realign_delta samples after the current one
    bool TRACK_realign = false;
    long long TRACK_realign_delta = 0;
    uint32_t TRACK_realign_SFN = 0;
    size_t TRACK_realign_sym = 0;
    long long TRACK_rx_ofs_of_first_sample = 0;

    
    // frames are demodulated in parallel, one per thread: downstream order doesn't matter
    // (the DCI workers already complete out of order, Summary bins results by rx_ofs)
    std::vector<std::thread> fft_thr;
    fft_thr.reserve(fft_n_threads);
    for (int i = 0; i < fft_n_threads; ++i)
        fft_thr.emplace_back(fft_worker,
            std::ref(fft_q), std::ref(in_q), std::ref(pool),
            N_FFT, Nsymb_frame, CP_slack, N, grid_ofs, Ts,
            std::cref(phase_comp), std::cref(sym_ofs));

    // wall-clock time of rx_ofs 0, for the RNTI logs; written before the first Result is pushed,
    // read by the main thread after popping one (out_q's mutex orders the two)
    std::chrono::system_clock::time_point stream_start;

    // receiver statistics for the per-second summary, sent to out_q as a RES_TICK every second
    // of stream time (rx_ofs/samp_rate); frames are counted in the second of the rx_ofs they are
    // shipped with, which is also the rx_ofs of their RES_FRAME_DONE
    const long long samples_per_sec = std::llround(samp_rate);
    long long tick_sec = 0;
    RxStatus rx_st{};
    double snr_lin_sum = 0.0;
    auto send_tick = [&](bool partial) {
        Result r;
        r.kind = RES_TICK;
        r.rx_ofs = tick_sec*samples_per_sec;
        rx_st.sec = tick_sec;
        rx_st.partial = partial;
        rx_st.state = state;
        rx_st.PCI = (state == TRACK) ? 3*TRACK_NID1 + TRACK_NID2 : -1;
        rx_st.f_SSB = freq_rx + (double)TRACK_SSB_k120*SCS;
        rx_st.snr_dB = rx_st.ssb_found ? 10.0*log10(snr_lin_sum/rx_st.ssb_found) : 0.0f;
        rx_st.cfo_Hz = (state == TRACK) ? Cfo_eF : 0.0f;
        rx_st.sfo_ppm = (state == TRACK && rx_ofs > TRACK_rx_ofs_of_first_sample)
                      ? (double)TRACK_samples_to_skip_cumsum/(double)(rx_ofs-TRACK_rx_ofs_of_first_sample)*1e6 : 0.0;
        rx_st.gain_dB = from_file ? NAN : gain_rx;
        rx_st.SFN = (state == TRACK && TRACK_MIB_valid) ? LFN : -1;
        r.rx = rx_st;
        out_q.push(std::move(r));
        rx_st = RxStatus{};
        snr_lin_sum = 0.0;
        tick_sec++;
    };

    // ---------------------------------------------------------------------------------------- RX thread
    // The RX thread owns the UHD streamer and the AGC/ACQUIRE/TRACK state machine. When it exits
    // (CTRL-C or RX timeout) it tears down the downstream pipeline and closes out_q, which ends
    // the result-printing loop in the main thread.
    auto rx_loop = [&]() {
        // start streaming
        uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
        stream_start = std::chrono::system_clock::now();
        if(!from_file){
            std::cout << "Starting streaming..." << std::endl;
            stream_cmd.stream_now = true;
            rx_stream->issue_stream_cmd(stream_cmd);
        }
    
        uhd::rx_metadata_t rx_md;

        std::cout << "Streaming (Press Ctrl+C to stop)..." << std::endl;

        while (!stop_signal_called) {

            // RX samples
            size_t rx_buf_len;
            if(from_file){
                if(file_format == "sc16"){
                    iq_file.read(reinterpret_cast<char*>(file_buf.data()), file_buf.size()*sizeof(int16_t));
                    rx_buf_len = iq_file.gcount()/(2*sizeof(int16_t));
                    for(size_t n = 0; n < rx_buf_len; n++) // same scaling as the UHD sc16->fc32 converter
                        rx_buf[n] = cfloat(file_buf[2*n], file_buf[2*n+1])*(1.0f/32768.0f);
                }
                else{
                    iq_file.read(reinterpret_cast<char*>(rx_buf.data()), rx_buf.size()*sizeof(cfloat));
                    rx_buf_len = iq_file.gcount()/sizeof(cfloat);
                }
                if(rx_buf_len == 0){
                    std::cout << "End of file" << std::endl;
                    break;
                }
            }
            else
                rx_buf_len = rx_stream->recv(&rx_buf.front(), rx_buf.size(), rx_md, 0.5);
            if (rx_md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
                std::cerr << "RX Timeout" << std::endl;
                break;
            }
            if (rx_md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
                if(verbose)
                    std::cerr << "O" << std::flush;
                rx_st.overflows++;
                if(state == TRACK)
                    rx_st.track_losses++;
                if(AGC_complete)
                    ACQ_start();
                else{
                    rx_ofs_last_AGC_round = rx_ofs;
                    state = AGC;
                }
                continue;
            }
            switch(state){
            case AGC:
        
                // ignore AGC_Tskip seconds to avoid DC transients
                if(rx_ofs > rx_ofs_last_AGC_round + AGC_Tskip*samp_rate){
                    for (size_t n = 0; n < rx_buf_len; n++) {
                        const cfloat sample = rx_buf[n];
                        if (sample.real() > AGC_max_real) AGC_max_real = sample.real();
                        if (sample.real() < AGC_min_real) AGC_min_real = sample.real();
                        if (sample.imag() > AGC_max_imag) AGC_max_imag = sample.imag();
                        if (sample.imag() < AGC_min_imag) AGC_min_imag = sample.imag();
                    }
                }
                // adjust RX gain every AGC_round_delay seconds
                if (rx_ofs > rx_ofs_last_AGC_round + AGC_round_delay*samp_rate){
                    //std::cout << "real: [" << AGC_min_real << "," << AGC_max_real
                    //          << "] imag: [" << AGC_min_imag << "," << AGC_max_imag << "]" << std::endl << std::flush;

                    float max_ampl = std::max(std::max(AGC_max_real,-AGC_min_real),std::max(AGC_max_imag,-AGC_min_imag));
                    float delta_rx_gain = rx_gain_range.step()*round(20.0*log10(AGC_target/max_ampl)/rx_gain_range.step());
                    gain_rx = std::min(rx_gain_range.stop(),std::max(gain_rx+delta_rx_gain,rx_gain_range.start()));
                    usrp->set_rx_gain(gain_rx);
                    //std::cout << "Set RX Gain (delta_rx_gain: " << delta_rx_gain << "): " << usrp->get_rx_gain() << " dB" << std::endl;
                    // the next round measures the peak at the new gain
                    AGC_max_real = AGC_min_real = AGC_max_imag = AGC_min_imag = 0;
                
                    rx_ofs_last_AGC_round += AGC_round_delay*samp_rate;
                    if(--AGC_rounds==0){
                        std::cout << "[AGC] RX Gain: " << usrp->get_rx_gain() << " dB" << std::endl;
                        AGC_complete = true;
                        ACQ_start(); // drops a search result taken at another gain
                    }
                }
                break; // state: AGC

            case ACQUIRE:

                // Check if SSB_search worker has a result ready
                {
                    std::unique_lock<std::mutex> lock(worker_mutex);
                    if (worker_state == RESULT_READY) {
                        SSB_info_t SSB_info = worker_result;
                        worker_state = IDLE;
                        lock.unlock();
                        const bool stale = ACQ_drop_result; // search started before the last ACQ_start
                        ACQ_drop_result = false;
                    
                        if(!stale && SSB_info.valid && SSB_info.SIR_est_dB > SIR_threshold_dB){
                            // next_frame_rx_ofs >= rx_ofs + rx_buf_len
                            last_frame_rx_ofs = SSB_info.ofs;
                            size_t k_next = 1 + (rx_ofs + rx_buf_len - last_frame_rx_ofs)/samples_per_frame;
                            next_frame_rx_ofs = last_frame_rx_ofs + k_next * samples_per_frame;
                            if(verbose)
                                std::cout <<"last frame rx_ofs: "<<last_frame_rx_ofs<<" rx_ofs: "<<rx_ofs <<" next_frame_rx_ofs: "<<next_frame_rx_ofs<<std::endl;
                            Y_ofs = 0;
                            Y_cur_symb = 0;
                            is_in_symbol = false;
                            Y_CP_ofs = 0;
                            Cfo_eD = 0.0; Cfo_eF2 = SSB_info.Cfo_est; Cfo_eF = Kpf * Cfo_eD + Cfo_eF2;
                            Sfo_eD = 0.0; Sfo_eF2 = 0;                Sfo_eF = Kpf * Sfo_eD + Sfo_eF2;
                            Sfo_rounds_since_last_update = 1.0f;
                            Cfo_rounds_since_last_update = 1.0f;
                            Cfo_eD_num = 0.0f;Cfo_eD_den = 0.0f;
                            TRACK_SSB_k120 = (int)round((SSB_info.f_SSB-freq_rx)/SCS);
                            TRACK_SSB_sym = SSB_info.symSSB;
                            TRACK_NID2 = SSB_info.NID2;
                            TRACK_NID1 = SSB_info.NID1;
                            TRACK_samples_to_skip = 0;
                            TRACK_samples_to_skip_cumsum = 0;
                            TRACK_rx_ofs_of_first_sample = rx_ofs;
                            TRACK_MIB_valid = false;
                            TRACK_cs0 = CORESET0Grid{};
                            TRACK_PBCH_countdown = 0;
                            TRACK_realign = false;
                            TRACK_H_ref_is_valid = false; // the reference belongs to the previous cell/timing
                            noise_est = 0.0;
                            if(verbose)
                                std::cout<<"k120:"<<TRACK_SSB_k120<<" sym:"<<TRACK_SSB_sym
                                         <<" NID2:"<<TRACK_NID2<<" NID1:"<<TRACK_NID1<<std::endl;
                            else
                                std::cout<<"Cell found: PCI "<<3*TRACK_NID1+TRACK_NID2<<", SSB at "<<std::fixed<<std::setprecision(3)
                                         <<SSB_info.f_SSB/1e6<<" MHz, SIR "<<std::setprecision(1)<<SSB_info.SIR_est_dB<<" dB"
                                         <<std::defaultfloat<<std::setprecision(6)<<std::endl;
                            rx_st.acquisitions++;
                            LFN = 0;
                            state = TRACK;
                            break;
                        }
                        if(!stale && !from_file && ++ACQ_fails >= ACQ_max_fails){
                            std::cout << "Acquisition failed (" << ACQ_fails << " SSB searches without a cell), adjusting the RX gain again" << std::endl;
                            rx_st.acq_failures++;
                            AGC_restart();
                            break;
                        }
                    }
                }

                if(!ACQ_is_skipping){
                    // add samples to the buffer up to ACQ_bufsize
                    if(ACQ_samples_in_buf == 0){ // && ACQ_buf_rx_ofs==(long long)-1){
                        //std::cout<<"NEW ACQ_buf_rx_ofs: "<<rx_ofs<<"("<<ACQ_buf_rx_ofs<<")"<<std::endl;
                        ACQ_buf_rx_ofs = rx_ofs;
                    }
                
                    size_t samples_to_copy = std::min(rx_buf_len, ACQ_bufsize - ACQ_samples_in_buf);
                    ACQ_buf->insert(ACQ_buf->end(), rx_buf.begin(), rx_buf.begin()+samples_to_copy);
                    ACQ_samples_in_buf += samples_to_copy;
                
                    if(ACQ_samples_in_buf == ACQ_bufsize){
                    
                        // debug
                        if(last_ACQ_buf_rx_ofs==(long long)-1){
                            if(verbose)
                                std::cout <<"rx_ofs init: "<< ACQ_buf_rx_ofs <<std::endl;
                            last_ACQ_buf_rx_ofs = ACQ_buf_rx_ofs;
                        }
                        else if(ACQ_buf_rx_ofs>last_ACQ_buf_rx_ofs){
                            if(verbose){
                                std::cout <<"rx_ofs incr: "<< ACQ_buf_rx_ofs - last_ACQ_buf_rx_ofs;
                                if(ACQ_buf_rx_ofs!=next_ACQ_buf_rx_ofs)
                                    std::cout << " [OFF BY "<< ACQ_buf_rx_ofs - next_ACQ_buf_rx_ofs<<" SAMPLES]"<<std::endl;
                                else
                                    std::cout <<std::endl;
                            }
                            last_ACQ_buf_rx_ofs = next_ACQ_buf_rx_ofs;
                        }
                        next_ACQ_buf_rx_ofs = last_ACQ_buf_rx_ofs + ACQ_bufsize+ACQ_samples_to_skip;
                    
                        // Hand the buffer over to the worker under lock.
                        // If the worker has not taken the previous buffer yet, wait.
                        {
                            std::unique_lock<std::mutex> lock(worker_mutex);
                            // The worker sets buffer_to_process = nullptr when it takes it.
                            // If it is still != nullptr, or if the previous search is still running, the worker
                            // is late: the buffer it is reading is the one that will be filled after the swap
                            if(buffer_to_process != nullptr || worker_state == PROCESSING){
                                std::cerr << "[WARNING] worker too slow, waiting..." << std::endl;
                                result_cv.wait(lock, [] { return buffer_to_process == nullptr && worker_state != PROCESSING; });
                            }
                            if(worker_state == RESULT_READY && ACQ_drop_result){
                                // the stale search ended in this very iteration: nobody will read its result
                                worker_state = IDLE;
                                ACQ_drop_result = false;
                            }
                            buffer_to_process = ACQ_buf;
                            ACQ_buf_in_search = ACQ_buf;
                            buffer_to_process_rx_ofs = ACQ_buf_rx_ofs;
                            //assert(ACQ_buf->size() == ACQ_bufsize);
                            worker_cv.notify_one();
                            // reading a file is faster than real time: wait for the result, so that the
                            // buffer under search is never overwritten and no result is lost
                            if(from_file)
                                result_cv.wait(lock, [] { return worker_state == RESULT_READY; });
                        }
                    
                        // Swap the buffers
                        ACQ_buf = ACQ_buf_is_buf1 ? &ACQ_buf2 : &ACQ_buf1;
                        ACQ_buf_is_buf1 = !ACQ_buf_is_buf1;
                        ACQ_buf->clear();
                        ACQ_samples_in_buf = 0;
                    
                        if(ACQ_samples_to_skip == 0){
                            if(samples_to_copy < rx_buf_len){
                                // the extra samples go into the other buffer
                                size_t remaining_samples = rx_buf_len - samples_to_copy;
                                ACQ_buf->insert(ACQ_buf->end(), rx_buf.begin()+samples_to_copy, rx_buf.begin()+rx_buf_len);
                                ACQ_samples_in_buf += remaining_samples;
                            }
                            //std::cout<<"+ACQ_buf_rx_ofs: "<<rx_ofs+samples_to_copy<<"("<<ACQ_buf_rx_ofs
                            //         <<") delta: "<<rx_ofs+samples_to_copy-ACQ_buf_rx_ofs<<std::endl;
                            ACQ_buf_rx_ofs = rx_ofs + samples_to_copy;
                        }
                        else{ // ACQ_samples_to_skip>0
                            size_t samples_to_skip = std::min(ACQ_samples_to_skip, rx_buf_len - samples_to_copy);
                            if(samples_to_skip < ACQ_samples_to_skip){
                                ACQ_is_skipping = true;
                                ACQ_skipped_samples = samples_to_skip;
                            }
                            else{
                                size_t remaining_samples2 = rx_buf_len - samples_to_copy - samples_to_skip;
                                ACQ_buf->insert(ACQ_buf->end(), rx_buf.begin()+samples_to_copy+samples_to_skip, rx_buf.begin()+rx_buf_len);
                                ACQ_samples_in_buf += remaining_samples2;
                                //std::cout<<"-ACQ_buf_rx_ofs: "<<rx_ofs+samples_to_copy+samples_to_skip<<"("<<ACQ_buf_rx_ofs
                                //         <<") delta: "<<rx_ofs + samples_to_copy + samples_to_skip-ACQ_buf_rx_ofs<<std::endl;
                                ACQ_buf_rx_ofs = rx_ofs + samples_to_copy + samples_to_skip; // OR HERE???
                            }
                        }
                    }
                }
                else{ // ACQ_is_skipping
                    size_t samples_to_skip2 = std::min(rx_buf_len, ACQ_samples_to_skip - ACQ_skipped_samples);
                    ACQ_skipped_samples += samples_to_skip2;
                    if(ACQ_skipped_samples == ACQ_samples_to_skip){
                        ACQ_is_skipping = false;
                        ACQ_buf->insert(ACQ_buf->end(), rx_buf.begin()+samples_to_skip2, rx_buf.begin()+rx_buf_len);
                        ACQ_samples_in_buf += rx_buf_len - samples_to_skip2;
                        //std::cout<<"/ACQ_buf_rx_ofs: "<<rx_ofs+samples_to_skip2<<"("<<ACQ_buf_rx_ofs
                        //         <<") delta: "<<rx_ofs + samples_to_skip2-ACQ_buf_rx_ofs
                        //         <<"(rx_buf_len: "<<rx_buf_len<<",samples_to_skip: "<<samples_to_skip2
                        //         <<",ACQ_samps_in_buf:"<<ACQ_samples_in_buf <<")"<<std::endl;
                        ACQ_buf_rx_ofs = rx_ofs + samples_to_skip2;
                    }
                }
                break; // state: ACQUIRE
            
            case TRACK:
                // TRACKING
                // - given GSCN, NID1, NID2, iSSB, we can track the frequency
                //   by evaluating the FFTs of only the symbols containing PSS and SSS:
                //   - the phase error between their channel estimates gives the CFO
                //   - comparing the phase difference of the SS with the one found in the
                //     acquisition phase gives the timing error
            
                // next_frame_rx_ofs is the index of the first sample of the first cyclic prefix of the frame
                if(rx_ofs + (long long)rx_buf_len > next_frame_rx_ofs){
                    long int rx_buf_ofs; // rx_buf_ofs is the offset within rx_buf of the next sample to use
                    long int ofs_in_frame; // ofs_in_frame is the offset within the frame
                    if(rx_ofs < next_frame_rx_ofs){
                        rx_buf_ofs = next_frame_rx_ofs - rx_ofs;
                        ofs_in_frame = 0;
                    }
                    else{
                        rx_buf_ofs = 0;
                        ofs_in_frame = rx_ofs - next_frame_rx_ofs;
                    }
                    //std::cout << "rx_buf_ofs: "<<rx_buf_ofs<<" rx_ofs: "<<rx_ofs<<" ofs_in_frame:"<<ofs_in_frame<<" sym_ofs[Y_cur_symb]:"<<sym_ofs[Y_cur_symb]<<std::endl;
                    while(rx_buf_ofs < (long long)rx_buf_len){
                        if(!is_in_symbol){ // N_CP

                            // if there are samples to skip to keep the frame alignment and we are in the first cyclic prefix:
                            // skip TRACK_samples_to_skip samples in the receive buffer
                            if(Y_cur_symb==0 && (TRACK_samples_to_skip>0 && rx_buf_ofs+TRACK_samples_to_skip < (long long)rx_buf_len)){
                                // std::cout << "+++ TRACK_sts:"<<TRACK_samples_to_skip << std::endl;
                                TRACK_samples_to_skip = 0;
                            }
                            // update ofs_in_frame and Y_CP_ofs as if -TRACK_samples_to_skip samples had already been received
                            if(Y_cur_symb==0 && (TRACK_samples_to_skip<0 && rx_buf_ofs==0)){
                                // std::cout << "--- TRACK_sts:"<< TRACK_samples_to_skip << std::endl;
                                // std::cout << "ofs_in_frame: " << ofs_in_frame << " TRACK_samples_to_skip: " << TRACK_samples_to_skip << std::endl;
                                ofs_in_frame -= TRACK_samples_to_skip;
                                Y_CP_ofs -= TRACK_samples_to_skip;
                                TRACK_samples_to_skip = 0;
                            }
                        
                            // Compute how many samples to skip (CP)
                            size_t symbol_start = sym_ofs[Y_cur_symb] - CP_slack;
                            size_t samples_to_skip3 = std::min(symbol_start - ofs_in_frame, rx_buf_len - rx_buf_ofs);
                            for(size_t j = 0; j < samples_to_skip3; j++)
                                Y_CP[Y_CP_ofs + j] = rx_buf[rx_buf_ofs+j];
                            Y_CP_ofs += samples_to_skip3; // index in Y_CP
                            rx_buf_ofs += samples_to_skip3; // rx_buf_ofs: index in rx_buf
                            ofs_in_frame += samples_to_skip3;
                            if(ofs_in_frame == (long int)symbol_start)
                                is_in_symbol = true;
                        }
                        else{ // is_in_symbol N_FFT
                        
                            // Compute how many samples to copy into the current symbol
                            size_t symbol_end = sym_ofs[Y_cur_symb] - CP_slack + N_FFT;
                            size_t samples_to_copy2 = std::min(symbol_end - ofs_in_frame, rx_buf_len - rx_buf_ofs);
                            std::copy(&rx_buf[rx_buf_ofs],&rx_buf[rx_buf_ofs+samples_to_copy2],Y+Y_ofs);
                            //for(size_t j = 0; j < samples_to_copy2; j++)
                            //    Y[Y_ofs+j] = rx_buf[rx_buf_ofs+j];
                            Y_ofs += samples_to_copy2;
                            rx_buf_ofs += samples_to_copy2;
                            ofs_in_frame += samples_to_copy2;
                        
                            if(ofs_in_frame == (long int)symbol_end){
                                // compute the correlation between Y_CP[Y_CP_ofs-N_CP+2*CP_slack:Y_CP_ofs-1] and Y[Y_ofs-N_CP+2*CP_slack:Y_ofs-1]
                                cfloat corr_num(0.0,0.0);
                                float corr_den1 = 0.0;
                                float corr_den2 = 0.0;
                                for(size_t n = 0; n < N_CP-2*CP_slack; n++){
                                    corr_num += Y[Y_ofs-1-n] * std::conj(Y_CP[Y_CP_ofs-1-n]);
                                    corr_den1 += std::norm(Y_CP[Y_CP_ofs-1-n]);
                                    corr_den2 += std::norm(Y[Y_ofs-1-n]);
                                }
                                Y_CP_corr[Y_cur_symb] = corr_num/sqrtf(corr_den1*corr_den2);
                                if(std::abs(Y_CP_corr[Y_cur_symb]) > CP_corr_thres){
                                    Cfo_eD_num += std::arg(Y_CP_corr[Y_cur_symb])/(2.0*M_PI*float(N_FFT)*Ts);
                                    Cfo_eD_den += 1.0f;
                                }
                                //float corr_val = std::abs(corr_num)/sqrtf(corr_den1*corr_den2);
                                //if(corr_val>0.9)
                                //    std::cout<<"CP corr["<<Y_cur_symb<<"]="<<std::abs(corr_num)<<"/sqrtf("<< corr_den1<<"*"<<corr_den2<<")="<<corr_val<<std::endl;
                                Y_CP_ofs = 0;
                                is_in_symbol = false;
                                Y_cur_symb++;
                                //if(Y_cur_symb==Nsymb_frame && Y_ofs!=NY)
                                //    std::cout << "Y_ofs:"<<Y_ofs<<" NY:"<<NY<<std::endl;
                            }
                        
                            // if the frame is full:
                            if(Y_ofs == NY){
                            
                                // compute the FFT of the symbols containing PSS and SSS
                                // - multiply them by w[n], which compensates Cfo_est and scales the amplitudes (floats with int16 amplitudes)
                                // - call fftw_execute(Y_symb_to_freq) and compensate the phase of the term depending on CP_slack + phase_comp
                                // - compute H_est_PSS and H_est_SSS and:
                                //    - the phase difference between the estimates is due to the residual CFO
                                //    - the first average H_est is stored as the reference
                                //    - the following ones are equalized (ZF?) with the initial estimate and the timing
                                //      error is estimated from the maximum with parabolic interpolation of the IFFT
                                size_t sym = TRACK_SSB_sym; // PSS
                                if(!is_Y_valid[sym]){
                                    ofdm_fft(sym,Y,phase_comp[sym],Cfo_eF,sym_ofs[sym],Ts,N_FFT,Y_symb_to_freq[sym],Y_tmpbuf,CP_slack);
                                    is_Y_valid[sym] = true;
                                }
                                // SSS
                                if(!is_Y_valid[sym+2]){
                                    ofdm_fft(sym+2,Y,phase_comp[sym+2],Cfo_eF,sym_ofs[sym+2],Ts,N_FFT,Y_symb_to_freq[sym+2],Y_tmpbuf,CP_slack);
                                    is_Y_valid[sym+2] = true;
                                }
                            
                                // PSS_freq[NID2][n]
                                // SSS_freq[NID2][NID1][n]
                                cfloat H_est_SS_corr(0.0,0.0);
                                cfloat H_est_SS_coh(0.0,0.0); // the same with every subcarrier weighted 1
                                float H_est_PSS_energy = 0.0;
                                float H_est_SSS_energy = 0.0;
                                for(size_t n = 0; n < SS_carriers; n++){
                                    // in the specs the PSS sequence spans -64 to 62 relative to k120
                                    // but the PSS sequence was built with a leading zero, so to align them we start from -65
                                    H_est_PSS[n] = Y[sym*N_FFT    +((TRACK_SSB_k120-65+n+N_FFT/2)%N_FFT)]*PSS_freq[TRACK_NID2][n];
                                    H_est_PSS_energy += std::norm(H_est_PSS[n]);
                                    H_est_SSS[n] = Y[(sym+2)*N_FFT+((TRACK_SSB_k120-65+n+N_FFT/2)%N_FFT)]*SSS_freq[TRACK_NID2][TRACK_NID1][n];
                                    H_est_SSS_energy += std::norm(H_est_SSS[n]);
                                    const cfloat z = H_est_SSS[n] * std::conj(H_est_PSS[n]);
                                    H_est_SS_corr += z;
                                    if(std::abs(z) > 0.0f)
                                        H_est_SS_coh += z/std::abs(z);
                                }
                                // with the SS block, H_est_SSS/H_est_PSS has the same phase on all its 127 subcarriers
                                // (the CFO rotation between the two symbols). The correlation alone is scale invariant:
                                // a few strong subcarriers (DC/LO spur, narrowband interferer) push it to 1 in frames
                                // without the SS block, while they move the phase coherence only by 1/127 each
                                const float H_est_SS_corr_val = std::abs(H_est_SS_corr)/sqrtf(H_est_PSS_energy*H_est_SSS_energy);
                                const float H_est_SS_coh_val = std::abs(H_est_SS_coh)/127.0f;
                                //std::cout<<"H_est_SS_corr: "<<std::abs(H_est_SS_corr)/sqrtf(H_est_PSS_energy*H_est_SSS_energy)
                                //         <<" CP corr PSS: "<<std::abs(Y_CP_corr[sym])<<" SSS: "<<std::abs(Y_CP_corr[sym+2]) << std::endl;
                            
                                if (H_est_SS_corr_val > SS_presence_thres && H_est_SS_coh_val > SS_presence_thres){
                                    // the frame contains the expected SS block
                                
                                    // update the noise variance estimate
                                    float noise_est_num = 0.0;
                                    float noise_est_den = 0.0;
                                    // FFT bin of subcarrier n (0..239) of the SS block, wrapped into [0,N_FFT) also
                                    // when the block starts below the lowest bin (signed: N_FFT may not be a power of 2)
                                    auto SSB_bin = [&](long n){
                                        const long M = (long)N_FFT;
                                        return (size_t)(((M/2 + TRACK_SSB_k120 - 120 + n)%M + M)%M);
                                    };
                                    for(size_t n = 0; n < 56; n++)
                                        noise_est_num += std::norm(Y[sym*N_FFT+SSB_bin(n)]);
                                    for(size_t n = 183; n < 240; n++)
                                        noise_est_num += std::norm(Y[sym*N_FFT+SSB_bin(n)]);
                                    noise_est_den += 113.0;
                                
                                    for(size_t n = 48; n < 56; n++)
                                        noise_est_num += std::norm(Y[(sym+2)*N_FFT+SSB_bin(n)]);
                                    for(size_t n = 183; n < 192; n++)
                                        noise_est_num += std::norm(Y[(sym+2)*N_FFT+SSB_bin(n)]);
                                    noise_est_den += 17.0;
                                    float noise_est_new = noise_est_num/noise_est_den;
                                
                                    if(noise_est/noise_est_new + noise_est_new/noise_est > 3.0)
                                        noise_est = noise_est_new;
                                    else
                                        noise_est = 0.95*noise_est + 0.05*noise_est_new;
                                
                                    rx_st.ssb_found++;
                                    snr_lin_sum += sqrtf(H_est_PSS_energy*H_est_SSS_energy)/(127.0*noise_est);
                                    if(verbose)
                                        std::cout<<"LFN: "<<std::setw(4)<<LFN<<" H_est_SS_corr: "<<H_est_SS_corr_val<<" coh: "<<H_est_SS_coh_val
                                        //         <<" PSS CP corr: "<<std::abs(Y_CP_corr[sym])<<" SSS CP corr: "<<std::abs(Y_CP_corr[sym+2]) 
                                                 << " SNR: "<< 10.0*log10(sqrtf(H_est_PSS_energy*H_est_SSS_energy)/(127.0*noise_est));

                                    //save_fc32_debug(H_est_PSS,"H_est");
                                    //save_fc32_debug(H_est_SSS,"H_est");

                                    cfloat alpha = H_est_SS_corr/std::abs(H_est_SS_corr); // to combine the two estimates constructively
                                    //std::cout << " angle(H_est_SS_corr): "<<std::arg(H_est_SS_corr);
                                    //float Cfo_eD_SS = std::arg(H_est_SS_corr)/(4.0*M_PI*float(N_FFT+N_CP)*Ts); // CHECK!!!!
                                    //std::cout << " Cfo_eD(SS): "<<Cfo_eD_SS;
                                    if(Cfo_eD_den > 0.0){
                                        //Cfo_eD = Cfo_eD_SS;
                                        Cfo_eD = Cfo_eD_num/Cfo_eD_den - Cfo_eF;
                                        Cfo_eD_num = 0.0f; Cfo_eD_den = 0.0f;
                                        Cfo_eF2 += Kif*Cfo_rounds_since_last_update*Cfo_eD;
                                        Cfo_eF = Kpf*Cfo_rounds_since_last_update*Cfo_eD + Cfo_eF2;
                                        if(verbose)
                                            std::cout << " Cfo_eF: "<<Cfo_eF;
                                        Cfo_rounds_since_last_update = 1.0f; // rounds since the last measurement
                                    }
                                    if(!TRACK_H_ref_is_valid){
                                        for(size_t n = 0; n < SS_carriers; n++)
                                            TRACK_H_ref[n] = alpha*H_est_PSS[n]+H_est_SSS[n];
                                        //save_fc32_debug(TRACK_H_ref,SS_carriers,"TRACK_H_ref");
                                        TRACK_H_ref_is_valid=true;
                                        if(verbose)
                                            std::cout << std::endl;
                                    }
                                    else{ // TRACK_H_ref holds a channel estimate based on the SS
                                        size_t OSF = 4;
                                        for(size_t n = 0; n < SS_carriers/2; n++) // -64:-1
                                            TRACK_H_tmp[OSF*SS_carriers-SS_carriers/2+n] = std::conj(TRACK_H_ref[n])*(alpha*H_est_PSS[n]+H_est_SSS[n]);
                                        for(size_t n = SS_carriers/2; n < OSF*SS_carriers-SS_carriers/2; n++)
                                            TRACK_H_tmp[n] = cfloat(0.0,0.0);
                                        for(size_t n = SS_carriers/2; n < SS_carriers; n++)
                                            TRACK_H_tmp[n-SS_carriers/2] = std::conj(TRACK_H_ref[n])*(alpha*H_est_PSS[n]+H_est_SSS[n]);
                                        //save_fc32_debug(TRACK_H_tmp,4*SS_carriers,"Hf");
                                        fftwf_execute(TRACK_H_tmp_to_time);
                                        //save_fc32_debug(TRACK_H_tmp,4*SS_carriers,"Ht");
                                        float max_val = 0.0;
                                        int max_pos = -1;
                                        for(size_t n = 0; n < OSF*SS_carriers; n++)
                                            if(std::abs(TRACK_H_tmp[n]) > max_val){
                                                max_val = std::abs(TRACK_H_tmp[n]);
                                                max_pos = (n+OSF*SS_carriers/2)%(OSF*SS_carriers)-OSF*SS_carriers/2;
                                            }
                                        float yp = std::abs(TRACK_H_tmp[(max_pos+1)%(OSF*SS_carriers)]);
                                        float y0 = max_val;
                                        float ym = std::abs(TRACK_H_tmp[(max_pos-1)%(OSF*SS_carriers)]);
                                        float time_ofs_est = (yp-ym)/2.0/(2.0*y0-ym-yp); // Mengali
                                        if(time_ofs_est > 0.5)
                                            time_ofs_est = 0.5;
                                        else if(time_ofs_est < -0.5)
                                            time_ofs_est = -0.5;
                                        Sfo_eD = (max_pos+time_ofs_est)/(float)OSF; // OSF==4
                                        Sfo_eF2 += Kif*Sfo_rounds_since_last_update*Sfo_eD;
                                        Sfo_eF = Kpf*Sfo_rounds_since_last_update*Sfo_eD + Sfo_eF2;
                                        Sfo_rounds_since_last_update = 1.0f; // rounds since the last measurement
                                        if(verbose)
                                            std::cout<<" Sfo_eF: "<<Sfo_eF; //<<std::endl;
                                        TRACK_samples_to_skip = std::max(-2,std::min((int)std::round(Sfo_eF*N_FFT/SS_carriers),2));
                                    
                                        TRACK_samples_to_skip_cumsum += TRACK_samples_to_skip;
                                        if(verbose){
                                            std::cout<<" sts: "<<TRACK_samples_to_skip;
                                            if(TRACK_samples_to_skip > 2 || TRACK_samples_to_skip < -2)
                                                std::cout << "*****************";
                                            std::cout<<std::endl;
                                        }
                                    }

                                    // PBCH: MIB and SFN until decoded, then a check every PBCH_check_period SS blocks
                                    if(TRACK_PBCH_countdown > 0)
                                        TRACK_PBCH_countdown--;
                                    else{
                                        const long k_lo = (long)N_FFT/2 + TRACK_SSB_k120 - 120; // subcarrier 0 of the SS block
                                        if(k_lo < 0 || k_lo + (long)SSB_carriers > (long)N_FFT){
                                            std::cout << "The SS block is not entirely in the received band: PBCH not decoded" << std::endl;
                                            TRACK_PBCH_countdown = std::numeric_limits<int>::max();
                                        }
                                        else{
                                            for(size_t l = 1; l < SSB_symbols; l += 2)
                                                if(!is_Y_valid[sym+l]){
                                                    ofdm_fft(sym+l,Y,phase_comp[sym+l],Cfo_eF,sym_ofs[sym+l],Ts,N_FFT,Y_symb_to_freq[sym+l],Y_tmpbuf,CP_slack);
                                                    is_Y_valid[sym+l] = true;
                                                }
                                            for(size_t l = 0; l < SSB_symbols; l++)
                                                std::copy(Y + (sym+l)*N_FFT + k_lo, Y + (sym+l)*N_FFT + k_lo + SSB_carriers,
                                                          SSB_grid.begin() + l*SSB_carriers);
                                            MIB mib;
                                            if(pbch.decode(SSB_grid.data(), 3*TRACK_NID1 + TRACK_NID2, mib)){
                                                rx_st.pbch_ok++;
                                                const bool changed = !TRACK_MIB_valid || mib.SSB_index != TRACK_MIB.SSB_index
                                                    || mib.HRF != TRACK_MIB.HRF || mib.scs_common != TRACK_MIB.scs_common
                                                    || mib.k_SSB != TRACK_MIB.k_SSB || mib.dmrs_typeA_pos != TRACK_MIB.dmrs_typeA_pos
                                                    || mib.coreset0 != TRACK_MIB.coreset0 || mib.ss0 != TRACK_MIB.ss0
                                                    || mib.cell_barred != TRACK_MIB.cell_barred;
                                                if(changed || verbose)
                                                    std::cout << (verbose ? "\n" : "") << MIB_to_string(mib, mu, TRACK_SSB_k120, N) << std::endl;
                                                if(changed){
                                                    // CORESET#0 is searched if it is entirely in the grid of the received frames
                                                    CORESET0 cs;
                                                    TRACK_cs0 = CORESET0Grid{};
                                                    if(CORESET0_config(mib, mu, TRACK_SSB_k120, N, cs)){
                                                        if(cs.same_scs && cs.k_start == std::floor(cs.k_start) && cs.k_start >= 0
                                                        && cs.k_start + 12*cs.N_RB <= (double)N){
                                                            TRACK_cs0.valid = true;
                                                            TRACK_cs0.k_start = (int)cs.k_start;
                                                            TRACK_cs0.N_RB = cs.N_RB;
                                                            TRACK_cs0.duration = cs.duration;
                                                            TRACK_cs0.PCI = 3*TRACK_NID1 + TRACK_NID2;
                                                            TRACK_cs0.n0 = cs.n0;
                                                            TRACK_cs0.frame_parity = cs.frame_parity;
                                                            TRACK_cs0.first_symbol = cs.first_symbol;
                                                        }
                                                        else
                                                            std::cout << "CORESET#0 is not in the grid of the received frames: not searched" << std::endl;
                                                    }
                                                }
                                                const size_t sym_SSB = mib.HRF*(70<<mu) + SSB_first_symbol(mu, mib.SSB_index);
                                                if(sym_SSB == TRACK_SSB_sym){
                                                    if(TRACK_MIB_valid && LFN != mib.SFN)
                                                        std::cout << "SFN " << mib.SFN << " instead of " << LFN << std::endl;
                                                    LFN = mib.SFN;
                                                }
                                                else{
                                                    TRACK_realign = true;
                                                    TRACK_realign_delta = (long long)sym_ofs[TRACK_SSB_sym] - (long long)sym_ofs[sym_SSB];
                                                    TRACK_realign_SFN = mib.SFN;
                                                    TRACK_realign_sym = sym_SSB;
                                                }
                                                TRACK_MIB = mib;
                                                TRACK_MIB_valid = true;
                                                TRACK_PBCH_countdown = PBCH_check_period;
                                            }
                                            else
                                                rx_st.pbch_fail++;
                                        }
                                    }
                                } // SS block present
                                else{
                                    rx_st.ssb_missed++;
                                    if(verbose)
                                        std::cout << "LFN: "<<std::setw(4)<<LFN<<" H_est_SS_corr: "<<H_est_SS_corr_val<<" coh: "<<H_est_SS_coh_val<<std::endl;
                                    Cfo_rounds_since_last_update++;
                                    Sfo_rounds_since_last_update++;
                                }
                            
                                // ------------------------------------------------------------------------------ pass Y to fft_worker thread
                                FFTJob ffj;
                                ffj.LFN        = LFN;
                                ffj.rx_ofs     = rx_ofs;
                                ffj.Y.assign(Y, Y + NY);
                                ffj.is_Y_valid = is_Y_valid;
                                ffj.Y_CP_corr  = Y_CP_corr;
                                ffj.mu         = mu;
                                ffj.noise_est  = noise_est;
                                ffj.Cfo_eF     = Cfo_eF;
                                ffj.cs0        = TRACK_cs0;
                                if(from_file ? fft_q.push(std::move(ffj)) : fft_q.try_push(std::move(ffj))){
                                    rx_st.frames_shipped++;
                                    // LFN is the SFN once the MIB is decoded, except in the frame that needs a
                                    // realignment (shipped with the old counter); no PDCCH search without noise_est
                                    const bool SFN_known = TRACK_MIB_valid && !TRACK_realign;
                                    if(noise_est > 0.0f && SFN_known){
                                        if(SFN_segment_open && LFN == (SFN_segments.back().last + 1)%1024){
                                            SFN_segments.back().last = LFN;
                                            SFN_segments.back().frames++;
                                        }
                                        else
                                            SFN_segments.push_back({LFN, LFN, 1});
                                        SFN_segment_open = true;
                                    }
                                    else{
                                        if(noise_est > 0.0f)
                                            frames_before_SFN_known++;
                                        SFN_segment_open = false;
                                    }
                                }
                                else{
                                    SFN_segment_open = false;
                                    rx_st.frames_dropped++;
                                    if(verbose)
                                        std::cerr << "D" << std::flush;
                                }

                                LFN++;
                                if(LFN == 1024){
                                    if(verbose)
                                        std::cout << "Estimated Sfo: "<<((double)TRACK_samples_to_skip_cumsum/(double)(rx_ofs-TRACK_rx_ofs_of_first_sample))*1e6 << "ppm Estimated Cfo: "<<-Cfo_eF/freq_rx*1e6<<"ppm"<<std::endl;
                                    LFN = 0;
                                }
                            
                                // prepare for the next frame
                                Y_cur_symb = 0;
                                Y_ofs = 0;
                                for(sym = 0; sym < Nsymb_frame; sym++)
                                    is_Y_valid[sym] = false;
                                ofs_in_frame -= samples_per_frame;
                                const long long frame_rx_ofs = next_frame_rx_ofs; // start of the frame just shipped
                                next_frame_rx_ofs += samples_per_frame + TRACK_samples_to_skip;
                            
                                if(Cfo_rounds_since_last_update > 32.0){
                                    std::cout << "Tracking lost (no SS block for 32 frames), "
                                              << (from_file ? "searching again" : "adjusting the RX gain and searching again") << std::endl;
                                    rx_st.track_losses++;
                                    AGC_restart();
                                }
                                else if(TRACK_realign){
                                    // restart the frames from the frame of the decoded SS block, as at the end of
                                    // ACQUIRE: the samples up to the start of the next of those frames are dropped
                                    TRACK_realign = false;
                                    const long long F = frame_rx_ofs + TRACK_realign_delta;
                                    const long long buf_end = rx_ofs + (long long)rx_buf_len;
                                    const long long k_next = F >= buf_end ? 0 : 1 + (buf_end - F)/(long long)samples_per_frame;
                                    next_frame_rx_ofs = F + k_next*(long long)samples_per_frame;
                                    std::cout << "Frame timing corrected: SS block #" << TRACK_MIB.SSB_index << " of half frame " << TRACK_MIB.HRF
                                              << " is in symbol " << TRACK_realign_sym << " of the frame, not " << TRACK_SSB_sym << std::endl;
                                    LFN = (TRACK_realign_SFN + k_next)%1024;
                                    TRACK_SSB_sym = TRACK_realign_sym;
                                    TRACK_samples_to_skip = 0;
                                    Y_CP_ofs = 0;
                                    rx_buf_ofs = rx_buf_len;
                                }
                                else{
                                    if(TRACK_samples_to_skip > 0){
                                        if(rx_buf_ofs+TRACK_samples_to_skip < (long)rx_buf_len){
                                            rx_buf_ofs += TRACK_samples_to_skip;
                                            TRACK_samples_to_skip = 0;
                                        }
                                        else{
                                            TRACK_samples_to_skip -= rx_buf_len-rx_buf_ofs;
                                            rx_buf_ofs = rx_buf_len;
                                        }
                                    }
                                    else if(TRACK_samples_to_skip < 0){
                                        if(rx_buf_ofs+TRACK_samples_to_skip >= 0){
                                            rx_buf_ofs += TRACK_samples_to_skip;
                                            TRACK_samples_to_skip = 0;
                                        }
                                        else{
                                            TRACK_samples_to_skip += rx_buf_ofs;
                                            rx_buf_ofs = 0;
                                        }
                                    }
                                }
                            } // Y_ofs==NY
                        } // is_in_symbol
                    } // rx_buf_ofs < rx_buf_len
                } // rx_ofs + rx_buf_len > next_frame_rx_ofs
                break; // state: TRACK
            } // switch(state)

            // update rx_ofs for the next loop
            rx_ofs += rx_buf_len;

            // close the seconds of stream time that are over
            while(tick_sec < rx_ofs/samples_per_sec)
                send_tick(false);
        }
        send_tick(true);
    
        if(!from_file){
            std::cout << "Stopping streaming..." << std::endl;
            stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
            rx_stream->issue_stream_cmd(stream_cmd);
        }

        /* tell the ACQUISITION worker to exit */ {
            std::lock_guard<std::mutex> lock(worker_mutex);
            worker_should_exit = true;
            worker_cv.notify_one();
        }
        if(ACQ_thread.joinable())
            ACQ_thread.join();

        fft_q.close();
        for (std::size_t i = 0; i < fft_thr.size(); ++i) fft_thr[i].join();

        // stop the PROCESSING workers
        in_q.close();
        pool.close();
        for (std::size_t i = 0; i < workers.size(); ++i) workers[i].join();
        out_q.close();
    };

    // the RX thread inherits the real-time priority and the --rx-cpus affinity of main
    std::thread rx_thr(rx_loop);

    // the main thread only prints results: drop it back to normal scheduling
    {
        sched_param sp{};
        sp.sched_priority = 0;
        pthread_setschedparam(pthread_self(), SCHED_OTHER, &sp);
    }

    // --------------- print results until the RX thread shuts the pipeline down
    // with --verbose the summary only feeds the web UI
    Summary summary(samp_rate, mu, stream_start, throughput_BWP, !verbose, webui_ptr);
    while (true) {
        std::optional<Result> res_opt = out_q.pop();
        if (!res_opt) break;
        const Result& r = *res_opt;
        if(rnti_logger && r.kind == RES_DCI)
            rnti_logger->log(r, stream_start);
        if(!verbose || webui_ptr)
            summary.add(r);
        if(!verbose)
            continue;
        std::ostringstream line; // single write, to avoid interleaving with the RX thread output
        if(r.kind == RES_DCI){
            line<<"  #"<<(int)r.slot<<"["<<(int)r.k_ini<<":"<<(int)r.k_fin-1<<"] AL"<<(int)r.AL<<" d"<<(int)r.duration<<(r.CORESET0 ? " C0" : "")<<" NID_est: 0x"
                <<std::hex<<r.DCI_NID<<" RNTI(A:"<<std::dec<<(int)r.DCI_len<<"): 0x"<<std::hex<<r.RNTI
                <<" DCI: "<<std::dec;
            for(size_t j = 0; j < r.DCI_len; j++)
                line << (char)('0'+r.DCI[j]);
            Grant g;
            if(throughput_BWP > 0 && decode_fallback(r, throughput_BWP, g)){
                line<<" RB "<<g.RB_start<<"+"<<g.L_RB<<" MCS "<<g.MCS<<" RV "<<g.RV<<" TBS ";
                if(g.MCS <= 28)
                    line<<std::lround(g.TBS);
                else
                    line<<"? (retx)";
            }
            line << '\n';
        }
        else if(r.kind == RES_TRS || r.kind == RES_CSI || r.kind == RES_CSI_FD2){
            const char* name = r.kind == RES_TRS ? "TRS" : r.kind == RES_CSI ? "CSI" : "CSIfd2";
            line<<"  "<<name<<"["<<r.slot*14+r.symbol<<"] k0:"<<r.k_ini<<" NID_est: "<<r.DCI_NID<<'\n';
        }
        std::cout << line.str() << std::flush;
    }
    summary.flush();
    rnti_logger.reset(); // close all the RNTI logs

//shutdown:

    rx_thr.join();
    ACQ_corr_done();

    // --file: the frames whose PDCCH was searched, by SFN (first..last wraps through 1023 -> 0 if last < first)
    if(from_file){
        for(const auto& seg : SFN_segments)
            std::cout << "SFN range: " << seg.first << ".." << seg.last << " (" << seg.frames << " frames)" << std::endl;
        if(SFN_segments.empty())
            std::cout << "SFN range: none (no frame searched with a known SFN)" << std::endl;
        if(frames_before_SFN_known)
            std::cout << "SFN range: " << frames_before_SFN_known
                      << " frames searched before the SFN was known (their DCIs have local frame numbers)" << std::endl;
    }
    
    fftwf_free(Y_CP_fftw);
    fftwf_free(Y_tmpbuf_fftw);
    fftwf_destroy_plan(TRACK_H_tmp_to_time);
    fftwf_free(TRACK_H_tmp_fftw);
    fftwf_free(TRACK_H_ref_fftw);
    
    for(size_t n = 0; n < Nsymb_frame; n++)
        fftwf_destroy_plan(Y_symb_to_freq[n]);
    fftwf_free(Y_fftw);
    
    // Save FFTW wisdom for next run
    if (fftwf_export_wisdom_to_filename(wisdom_file))
        std::cout << "Saved FFTW wisdom to " << wisdom_file << std::endl;
    else
        std::cerr << "Warning: Failed to save FFTW wisdom" << std::endl;
    
    std::cout << " Done!" << std::endl;
    return EXIT_SUCCESS;
}
