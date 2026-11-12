#ifndef RNTI_LOG_H_
#define RNTI_LOG_H_

#include <chrono>
#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "pipeline.h"

// Per-RNTI DCI log (--rnti-log): one file <dir>/<RNTI hex>.log per RNTI, one line per DCI
// "{date-time} {LFN:04}.{slot:02} {RNTI} {NID} {AL} {duration} {DCI bits}", where date-time is the
// local time of the DCI (stream start + rx_ofs/samp_rate + slot*slot duration) as YYYYmmdd-HHMMSS.uuuuuu.
// A file is closed when its RNTI hasn't been seen for timeout_samples of stream time (rx_ofs,
// so a file replay behaves like a live capture) and reopened in append mode if the RNTI comes
// back; files are truncated the first time their RNTI shows up in a run.
class RNTILogger {
public:
    // dir_prefix gets the current local time appended, e.g. "logs" -> "logs-20261008-095230"
    RNTILogger(const std::string& dir_prefix, double samp_rate, int mu, double timeout_s);
    ~RNTILogger() { close_all(); }

    const std::string& dir() const { return dir_; }
    // stream_start: wall-clock time of the sample with rx_ofs 0
    void log(const Result& r, std::chrono::system_clock::time_point stream_start);
    void close_all() { files_.clear(); }

private:
    struct Entry {
        std::ofstream f;
        long long last_rx_ofs;
    };
    void close_idle();

    std::string dir_;
    double samp_rate_;
    double slot_s_; // slot duration [s]
    long long timeout_;
    long long newest_rx_ofs_ = 0;
    long long last_sweep_rx_ofs_ = 0;
    std::unordered_map<uint16_t, Entry> files_; // open files
    std::unordered_set<uint16_t> seen_;          // RNTIs already logged in this run
};

#endif
