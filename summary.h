#ifndef SUMMARY_H_
#define SUMMARY_H_

#include <chrono>
#include <cstdint>
#include <map>
#include <set>
#include <vector>
#include "pipeline.h"

class WebUI;

// Per-second console summary (default output; --verbose prints the raw results instead).
// Results are binned by stream time (rx_ofs/samp_rate), so a file replay is summarized like a
// live capture. Second s is printed once the RX thread has closed it (RES_TICK) and every frame
// it shipped in that second has been processed (RES_FRAME_DONE), so late DCIs from slow workers
// still land in the right second. Each second goes to the console and/or, as one JSON datagram,
// to the --webui dashboard bridge.
class Summary {
public:
    // stream_start: wall-clock time of rx_ofs 0, read only when a second is printed
    // fallback_BWP > 0: estimate the DL/UL throughput from the fallback DCIs of a BWP of that many RBs
    // console: print the summary to stdout; webui: also send it as JSON (nullptr: don't)
    Summary(double samp_rate, int mu, const std::chrono::system_clock::time_point& stream_start, int fallback_BWP = 0,
            bool console = true, WebUI* webui = nullptr);

    void add(const Result& r);
    void flush(); // print every pending second (end of run)

private:
    struct RNTIStats {
        unsigned DL = 0; // DCI first bit 1
        unsigned UL = 0; // DCI first bit 0
        std::set<unsigned> NIDs; // PDCCH DM-RS scrambling IDs
        double DL_bits = 0, UL_bits = 0; // estimated TBS sums
    };
    struct Second {
        bool have_tick = false;
        RxStatus rx{};
        unsigned frames_done = 0;
        unsigned DL = 0, UL = 0;
        double DL_bits = 0, UL_bits = 0; // estimated TBS sums
        std::set<std::pair<unsigned, unsigned>> PDCCH_slots; // (LFN, slot)
        std::map<uint16_t, RNTIStats> RNTIs;
    };

    // all the RNTIs of a second, busiest first
    static std::vector<std::pair<uint16_t, const RNTIStats*>> busiest_first(const Second& s);
    void report(long long sec, const Second& s) const;
    void print(long long sec, const Second& s) const;
    void send_json(long long sec, const Second& s) const; // WITH_WEBUI builds only

    long long sps_;     // samples per second
    int slots_per_frame_;
    int fallback_BWP_;
    bool console_;
    WebUI* webui_;
    const std::chrono::system_clock::time_point& stream_start_;
    std::map<long long, Second> pending_;
    long long next_sec_ = 0; // seconds before this one were already printed
};

#endif
