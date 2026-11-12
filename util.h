#ifndef UTIL_H_
#define UTIL_H_

#include <string>
#include <thread>
#include <vector>
#include <sched.h>
#include <sys/stat.h>
#include "nr_common.h"

// --verbose: per-frame, per-DCI and SSB search debug output instead of the per-second summary
extern bool verbose;

// useful for debug
inline bool file_exists(const std::string& filename)
{
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

void save_fc32_debug(const cfloat* data, size_t length, const std::string& base_name = "debug");
void save_fc32_debug(const std::vector<cfloat>& data, const std::string& base_name = "debug");

// CPU affinity of the threads (--rx-cpus, --worker-cpus):
// RX thread, SSB worker and fft_workers on rx_cpus; DCI workers and SSB search helpers on worker_cpus
extern cpu_set_t rx_cpus;
extern cpu_set_t worker_cpus;

// parses a CPU list as taskset -c does ("0,2,4,6", "8-14:2,16-31") or "all" (the CPUs the process
// may use at startup) into cpus; CPUs the process may not use are dropped with a warning.
// Returns false, with a message in err, on a syntax error or when no usable CPU is left.
bool parse_cpu_list(const std::string& list, cpu_set_t& cpus, std::string& err);
std::string cpu_list_string(const cpu_set_t& cpus);
// pins the calling thread to cpus
void pin_self(const cpu_set_t& cpus);

#endif
