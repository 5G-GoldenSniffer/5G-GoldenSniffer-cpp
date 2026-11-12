#ifndef WEBUI_H_
#define WEBUI_H_

#include <string>
#include <vector>

// --webui [host:]port: one JSON datagram per summary second, sent over UDP to the dashboard
// bridge (webui/server.py). Sending never blocks the summary: a missing listener or a full
// socket buffer just drops the datagram.
class WebUI {
public:
    static constexpr size_t max_datagram = 60000; // below the 65507-byte UDP payload limit

    // dest: "host:port", ":port" or "port" (host defaults to 127.0.0.1); throws on a bad address
    // cfg_json: static run configuration, a JSON object sent in every datagram as "cfg"
    WebUI(const std::string& dest, const std::string& cfg_json);
    ~WebUI();
    WebUI(const WebUI&) = delete;
    WebUI& operator=(const WebUI&) = delete;

    const std::string& dest() const { return dest_; }
    const std::string& cfg() const { return cfg_; }
    void send(const std::string& json);
    // CPU load since the previous call, as a JSON object (null on the first call): "system" and
    // "process" (this program) in % of all the cores, "cores" the per-core load in %
    std::string cpu_json();

private:
    int fd_ = -1;
    std::string dest_, cfg_;
    struct CPUTimes { unsigned long long busy = 0, total = 0; };
    std::vector<CPUTimes> cpu_prev_; // [0] all the cores, [1+n] core n
    unsigned long long proc_prev_ = 0;
};

#endif
