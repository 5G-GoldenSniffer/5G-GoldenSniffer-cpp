#include <cerrno>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <netdb.h>
#include <sys/socket.h>
#include <unistd.h>
#include "webui.h"

WebUI::WebUI(const std::string& dest, const std::string& cfg_json) : dest_(dest), cfg_(cfg_json)
{
    std::string host = "127.0.0.1", port = dest;
    size_t colon = dest.rfind(':');
    if(colon != std::string::npos){
        if(colon > 0)
            host = dest.substr(0, colon);
        port = dest.substr(colon+1);
    }
    if(host.size() > 2 && host.front() == '[' && host.back() == ']') // [IPv6]:port
        host = host.substr(1, host.size()-2);
    if(port.empty() || port.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("--webui: expected [host:]port, got \"" + dest + "\"");

    addrinfo hints{}, *res = nullptr;
    hints.ai_family = AF_UNSPEC;
    hints.ai_socktype = SOCK_DGRAM;
    int err = getaddrinfo(host.c_str(), port.c_str(), &hints, &res);
    if(err)
        throw std::runtime_error("--webui " + dest + ": " + gai_strerror(err));
    for(addrinfo* ai = res; ai; ai = ai->ai_next){
        fd_ = socket(ai->ai_family, ai->ai_socktype | SOCK_NONBLOCK, ai->ai_protocol);
        if(fd_ < 0)
            continue;
        if(connect(fd_, ai->ai_addr, ai->ai_addrlen) == 0)
            break;
        close(fd_);
        fd_ = -1;
    }
    freeaddrinfo(res);
    if(fd_ < 0)
        throw std::runtime_error("--webui " + dest + ": can't open a UDP socket: " + strerror(errno));
}

WebUI::~WebUI()
{
    if(fd_ >= 0)
        close(fd_);
}

void WebUI::send(const std::string& json)
{
    // errors (ECONNREFUSED while the bridge isn't running, EAGAIN) are ignored: the next second retries
    (void)::send(fd_, json.data(), json.size(), MSG_DONTWAIT | MSG_NOSIGNAL);
}

std::string WebUI::cpu_json()
{
    // /proc/stat: "cpu[N] user nice system idle iowait irq softirq steal ..." in USER_HZ ticks
    std::vector<CPUTimes> now;
    std::ifstream stat("/proc/stat");
    std::string line;
    while(std::getline(stat, line) && line.compare(0, 3, "cpu") == 0){
        std::istringstream in(line.substr(line.find(' ')));
        unsigned long long v, idle = 0, total = 0;
        for(int i = 0; i < 8 && in >> v; i++){
            total += v;
            if(i == 3 || i == 4) // idle, iowait
                idle += v;
        }
        now.push_back({total - idle, total});
    }
    // /proc/self/stat: utime and stime are fields 14 and 15, after the parenthesized command name
    unsigned long long proc = 0;
    std::ifstream self("/proc/self/stat");
    if(std::getline(self, line)){
        std::istringstream in(line.substr(line.rfind(')') + 2));
        std::string f;
        unsigned long long ut = 0, st = 0;
        for(int i = 3; i <= 15 && in >> f; i++){
            if(i == 14) ut = std::stoull(f);
            if(i == 15) st = std::stoull(f);
        }
        proc = ut + st;
    }

    std::string out = "null";
    if(now.size() > 1 && now.size() == cpu_prev_.size() && now[0].total > cpu_prev_[0].total){
        auto pct = [&](size_t i){
            const double dt = now[i].total - cpu_prev_[i].total;
            return dt > 0 ? 100.0*(now[i].busy - cpu_prev_[i].busy)/dt : 0.0;
        };
        char buf[32];
        std::ostringstream j;
        snprintf(buf, sizeof(buf), "%.1f", pct(0));
        j << "{\"system\":" << buf;
        snprintf(buf, sizeof(buf), "%.1f", 100.0*(proc - proc_prev_)/(now[0].total - cpu_prev_[0].total));
        j << ",\"process\":" << buf << ",\"cores\":[";
        for(size_t i = 1; i < now.size(); i++){
            snprintf(buf, sizeof(buf), "%s%.0f", i > 1 ? "," : "", pct(i));
            j << buf;
        }
        j << "]}";
        out = j.str();
    }
    cpu_prev_ = std::move(now);
    proc_prev_ = proc;
    return out;
}
