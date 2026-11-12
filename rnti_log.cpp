#include <cstdio>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include "rnti_log.h"

// local time as YYYYmmdd-HHMMSS, plus .uuuuuu if with_us
static std::string format_time(std::chrono::system_clock::time_point t, bool with_us)
{
    using namespace std::chrono;
    const auto us = duration_cast<microseconds>(t.time_since_epoch()).count();
    std::time_t secs = us/1000000;
    std::tm tm;
    localtime_r(&secs, &tm);
    char buf[32];
    size_t n = std::strftime(buf, sizeof(buf), "%Y%m%d-%H%M%S", &tm);
    if(with_us)
        snprintf(buf+n, sizeof(buf)-n, ".%06lld", (long long)(us%1000000));
    return buf;
}

RNTILogger::RNTILogger(const std::string& dir_prefix, double samp_rate, int mu, double timeout_s)
    : dir_(dir_prefix + "-" + format_time(std::chrono::system_clock::now(), false)),
      samp_rate_(samp_rate), slot_s_(1e-3/(1<<mu)), timeout_((long long)(timeout_s*samp_rate))
{
    std::error_code ec;
    std::filesystem::create_directories(dir_, ec);
    if(ec)
        throw std::runtime_error("can't create RNTI log directory " + dir_ + ": " + ec.message());
}

void RNTILogger::log(const Result& r, std::chrono::system_clock::time_point stream_start)
{
    // results arrive slightly out of order from the workers: track the newest one
    if(r.rx_ofs > newest_rx_ofs_)
        newest_rx_ofs_ = r.rx_ofs;

    auto it = files_.find(r.RNTI);
    if(it == files_.end()){
        char name[16];
        snprintf(name, sizeof(name), "%04x.log", r.RNTI);
        const std::string path = dir_ + "/" + name;
        const bool first = seen_.insert(r.RNTI).second;
        it = files_.emplace(r.RNTI, Entry{std::ofstream(path, first ? std::ios::trunc : std::ios::app), r.rx_ofs}).first;
        if(!it->second.f)
            std::cerr << "[WARNING] can't open RNTI log " << path << std::endl;
    }
    Entry& e = it->second;
    if(r.rx_ofs > e.last_rx_ofs)
        e.last_rx_ofs = r.rx_ofs;

    const auto t = stream_start + std::chrono::duration_cast<std::chrono::system_clock::duration>(
                                      std::chrono::duration<double>(r.rx_ofs/samp_rate_ + r.slot*slot_s_));
    char head[64];
    snprintf(head, sizeof(head), " %04u.%02u 0x%04x 0x%04x %u %u ",
             (unsigned)r.LFN, (unsigned)r.slot, (unsigned)r.RNTI, (unsigned)r.DCI_NID, (unsigned)r.AL, (unsigned)r.duration);
    std::string line = format_time(t, true) + head;
    for(size_t j = 0; j < r.DCI_len; j++)
        line += (char)('0'+r.DCI[j]);
    line += '\n';
    e.f << line;

    // sweep idle files at most once per timeout/10 of stream time
    if(newest_rx_ofs_ - last_sweep_rx_ofs_ >= timeout_/10){
        last_sweep_rx_ofs_ = newest_rx_ofs_;
        close_idle();
    }
}

void RNTILogger::close_idle()
{
    for(auto it = files_.begin(); it != files_.end(); ){
        if(newest_rx_ofs_ - it->second.last_rx_ofs > timeout_)
            it = files_.erase(it); // ofstream destructor flushes and closes
        else
            ++it;
    }
}
