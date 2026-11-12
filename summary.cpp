#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "summary.h"
#include "fallback_dci.h"
#ifdef WITH_WEBUI
#include "webui.h"
#endif

static const char* const state_name[3] = {"AGC", "ACQUIRE", "TRACK"};

Summary::Summary(double samp_rate, int mu, const std::chrono::system_clock::time_point& stream_start, int fallback_BWP,
                 bool console, WebUI* webui)
    : sps_(std::llround(samp_rate)), slots_per_frame_(10<<mu), fallback_BWP_(fallback_BWP), console_(console), webui_(webui),
      stream_start_(stream_start) {}

void Summary::add(const Result& r)
{
    long long sec = (r.kind == RES_TICK) ? r.rx.sec : r.rx_ofs/sps_;
    if(sec < next_sec_) // can't happen if the frame accounting is right: count it in the next second
        sec = next_sec_;
    Second& s = pending_[sec];

    switch(r.kind){
    case RES_TICK:
        s.have_tick = true;
        s.rx = r.rx;
        break;
    case RES_FRAME_DONE:
        s.frames_done++;
        break;
    case RES_DCI: {
        s.PDCCH_slots.insert({r.LFN, r.slot});
        RNTIStats& st = s.RNTIs[r.RNTI];
        if(DCI_is_DL(r)){
            st.DL++;
            s.DL++;
        }
        else{
            st.UL++;
            s.UL++;
        }
        st.NIDs.insert(r.DCI_NID);
        Grant g;
        if(fallback_BWP_ > 0 && decode_fallback(r, fallback_BWP_, g) && g.RV == 0){ // retransmissions (RV != 0) not counted
            (g.DL ? st.DL_bits : st.UL_bits) += g.TBS;
            (g.DL ? s.DL_bits : s.UL_bits) += g.TBS;
        }
        break;
    }
    default: // NZP-CSI-RS/TRS are only printed with --verbose
        break;
    }

    // print the seconds that are complete, in order
    while(!pending_.empty()){
        auto it = pending_.begin();
        if(it->first != next_sec_ || !it->second.have_tick || it->second.frames_done < it->second.rx.frames_shipped)
            break;
        report(it->first, it->second);
        pending_.erase(it);
        next_sec_++;
    }
}

void Summary::flush()
{
    for(auto& p : pending_)
        report(p.first, p.second);
    pending_.clear();
}

std::vector<std::pair<uint16_t, const Summary::RNTIStats*>> Summary::busiest_first(const Second& s)
{
    std::vector<std::pair<uint16_t, const RNTIStats*>> RNTIs;
    for(const auto& p : s.RNTIs)
        RNTIs.push_back({p.first, &p.second});
    std::sort(RNTIs.begin(), RNTIs.end(), [](const auto& a, const auto& b){
        const unsigned na = a.second->DL + a.second->UL, nb = b.second->DL + b.second->UL;
        return na != nb ? na > nb : a.first < b.first; });
    return RNTIs;
}

void Summary::report(long long sec, const Second& s) const
{
    if(s.rx.partial && s.frames_done == 0) // nothing to report in the last fraction of a second
        return;
    if(console_)
        print(sec, s);
#ifdef WITH_WEBUI
    if(webui_)
        send_json(sec, s);
#endif
}

void Summary::print(long long sec, const Second& s) const
{
    std::ostringstream out;
    char buf[256];

    // local time of the start of the second
    std::time_t t = std::chrono::system_clock::to_time_t(stream_start_ + std::chrono::seconds(sec));
    std::tm tm;
    localtime_r(&t, &tm);
    std::strftime(buf, sizeof(buf), "%H:%M:%S", &tm);
    out << buf;
    snprintf(buf, sizeof(buf), " +%llds%s ", sec, s.rx.partial ? " (end)" : "");
    out << buf;

    if(!s.have_tick){
        out << " (no receiver status)\n";
    }
    else if(s.rx.state != 2){
        out << " " << state_name[s.rx.state < 3 ? s.rx.state : 1];
        if(s.rx.state == 0)
            out << "  adjusting RX gain";
        else
            out << "  searching for a cell";
        if(!std::isnan(s.rx.gain_dB)){
            snprintf(buf, sizeof(buf), "  gain %.1f dB", s.rx.gain_dB);
            out << buf;
        }
        if(s.rx.track_losses)
            out << "  [tracking lost]";
        if(s.rx.acq_failures)
            out << "  [acquisition failed]";
        if(s.rx.overflows)
            out << "  [" << s.rx.overflows << " overflows]";
        out << '\n';
    }
    else{
        snprintf(buf, sizeof(buf), " TRACK  PCI %d @ %.3f MHz  SSB in %u/%u frames",
                 s.rx.PCI, s.rx.f_SSB/1e6, s.rx.ssb_found, s.rx.ssb_found + s.rx.ssb_missed);
        out << buf;
        if(s.rx.SFN >= 0){
            snprintf(buf, sizeof(buf), "  SFN %d", s.rx.SFN);
            out << buf;
        }
        else
            out << "  SFN ?";
        if(s.rx.pbch_fail){
            snprintf(buf, sizeof(buf), "  [PBCH %u/%u]", s.rx.pbch_ok, s.rx.pbch_ok + s.rx.pbch_fail);
            out << buf;
        }
        if(s.rx.ssb_found){
            snprintf(buf, sizeof(buf), "  SNR %.1f dB", s.rx.snr_dB);
            out << buf;
        }
        snprintf(buf, sizeof(buf), "  CFO %+.0f Hz  SFO %+.2f ppm", s.rx.cfo_Hz, s.rx.sfo_ppm);
        out << buf;
        if(!std::isnan(s.rx.gain_dB)){
            snprintf(buf, sizeof(buf), "  gain %.1f dB", s.rx.gain_dB);
            out << buf;
        }
        if(s.rx.acquisitions)
            out << "  [acquired]";
        if(s.rx.track_losses)
            out << "  [" << s.rx.track_losses << "x tracking lost]";
        if(s.rx.overflows)
            out << "  [" << s.rx.overflows << " overflows]";
        if(s.rx.frames_dropped)
            out << "  [" << s.rx.frames_dropped << " frames dropped]";
        out << '\n';
    }

    if(s.frames_done){
        snprintf(buf, sizeof(buf), "    %u frames, PDCCH in %zu/%u slots, %u DCIs (DL %u, UL %u) from %zu RNTIs\n",
                 s.frames_done, s.PDCCH_slots.size(), s.frames_done*slots_per_frame_, s.DL+s.UL, s.DL, s.UL, s.RNTIs.size());
        out << buf;
    }
    // bits over the frames actually processed (10 ms each), so partial seconds aren't underestimated
    const double T = s.frames_done*10e-3;
    if(fallback_BWP_ > 0 && T > 0){
        snprintf(buf, sizeof(buf), "    throughput DL ~%.2f Mbit/s, UL ~%.2f Mbit/s (fallback DCIs)\n", s.DL_bits/T/1e6, s.UL_bits/T/1e6);
        out << buf;
    }

    // only the busiest RNTIs, so the console stays readable (the JSON for the web UI is not limited here)
    const size_t max_RNTIs = 41;
    const auto RNTIs = busiest_first(s);
    for(size_t i = 0; i < RNTIs.size() && i < max_RNTIs; i++){
        const auto& p = RNTIs[i];
        const uint16_t RNTI = p.first;
        const RNTIStats& st = *p.second;
        const unsigned count = st.DL + st.UL;
        const char* name = RNTI == 0xffff ? "SI-RNTI" : RNTI == 0xfffe ? "P-RNTI" : "";
        snprintf(buf, sizeof(buf), "    RNTI 0x%04x %-7s scrID", RNTI, name);
        out << buf;
        bool first = true;
        for(unsigned nid : st.NIDs){
            snprintf(buf, sizeof(buf), "%s0x%04x", first ? " " : ",", nid);
            out << buf;
            first = false;
        }
        snprintf(buf, sizeof(buf), "  %5u DCI%s  DL %5u  UL %5u", count, count == 1 ? " " : "s", st.DL, st.UL);
        out << buf;
        if(fallback_BWP_ > 0 && T > 0 && (st.DL_bits > 0 || st.UL_bits > 0)){
            snprintf(buf, sizeof(buf), "  DL ~%.2f UL ~%.2f Mbit/s", st.DL_bits/T/1e6, st.UL_bits/T/1e6);
            out << buf;
        }
        out << '\n';
    }
    if(RNTIs.size() > max_RNTIs)
        out << "    ... " << RNTIs.size() - max_RNTIs << " more RNTIs\n";

    std::cout << out.str() << std::flush;
}

#ifdef WITH_WEBUI
// JSON number, null for NaN/inf (not valid JSON)
static std::string jnum(double x, const char* fmt = "%.6g")
{
    if(!std::isfinite(x))
        return "null";
    char buf[32];
    snprintf(buf, sizeof(buf), fmt, x);
    return buf;
}

// One datagram per second, schema "nr-scanner/1" (documented in webui/README.md)
void Summary::send_json(long long sec, const Second& s) const
{
    const long long t_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        (stream_start_ + std::chrono::seconds(sec)).time_since_epoch()).count();
    const double T = s.frames_done*10e-3;
    const bool tput = fallback_BWP_ > 0 && T > 0;
    std::ostringstream j;
    j << "{\"schema\":\"nr-scanner/1\",\"t\":" << t_ms << ",\"sec\":" << sec
      << ",\"partial\":" << (s.rx.partial ? "true" : "false")
      << ",\"cfg\":" << webui_->cfg() << ",\"cpu\":" << webui_->cpu_json();
    if(s.have_tick){
        const RxStatus& rx = s.rx;
        j << ",\"rx\":{\"state\":\"" << state_name[rx.state < 3 ? rx.state : 1] << "\""
          << ",\"pci\":" << (rx.state == 2 ? std::to_string(rx.PCI) : "null")
          << ",\"sfn\":" << (rx.state == 2 && rx.SFN >= 0 ? std::to_string(rx.SFN) : "null")
          << ",\"f_ssb_hz\":" << (rx.state == 2 ? jnum(rx.f_SSB, "%.1f") : "null")
          << ",\"ssb_found\":" << rx.ssb_found << ",\"ssb_missed\":" << rx.ssb_missed
          << ",\"snr_db\":" << (rx.state == 2 && rx.ssb_found ? jnum(rx.snr_dB, "%.2f") : "null")
          << ",\"cfo_hz\":" << (rx.state == 2 ? jnum(rx.cfo_Hz, "%.1f") : "null")
          << ",\"sfo_ppm\":" << (rx.state == 2 ? jnum(rx.sfo_ppm, "%.4f") : "null")
          << ",\"gain_db\":" << jnum(rx.gain_dB, "%.1f")
          << ",\"frames_shipped\":" << rx.frames_shipped << ",\"frames_dropped\":" << rx.frames_dropped
          << ",\"overflows\":" << rx.overflows << ",\"acquisitions\":" << rx.acquisitions
          << ",\"track_losses\":" << rx.track_losses << ",\"acq_failures\":" << rx.acq_failures << "}";
    }
    else
        j << ",\"rx\":null";
    j << ",\"frames\":" << s.frames_done << ",\"slots\":" << s.frames_done*slots_per_frame_
      << ",\"pdcch_slots\":" << s.PDCCH_slots.size()
      << ",\"dci\":{\"dl\":" << s.DL << ",\"ul\":" << s.UL << "}"
      << ",\"throughput\":";
    if(tput)
        j << "{\"dl_bps\":" << jnum(s.DL_bits/T, "%.0f") << ",\"ul_bps\":" << jnum(s.UL_bits/T, "%.0f") << "}";
    else
        j << "null";
    j << ",\"rnti_count\":" << s.RNTIs.size() << ",\"rntis\":[";

    // busiest first, as many as fit in one datagram (the rest only counts in rnti_count)
    std::string out = j.str(), entry;
    const std::string tail_fmt = "],\"rntis_truncated\":false}";
    bool truncated = false, first = true;
    for(const auto& p : busiest_first(s)){
        const RNTIStats& st = *p.second;
        std::ostringstream e;
        e << (first ? "" : ",") << "{\"rnti\":" << p.first << ",\"nids\":[";
        bool first_nid = true;
        for(unsigned nid : st.NIDs){
            e << (first_nid ? "" : ",") << nid;
            first_nid = false;
        }
        e << "],\"dl\":" << st.DL << ",\"ul\":" << st.UL;
        if(tput)
            e << ",\"dl_bps\":" << jnum(st.DL_bits/T, "%.0f") << ",\"ul_bps\":" << jnum(st.UL_bits/T, "%.0f");
        e << "}";
        entry = e.str();
        if(out.size() + entry.size() + tail_fmt.size() + 1 > WebUI::max_datagram){
            truncated = true;
            break;
        }
        out += entry;
        first = false;
    }
    out += truncated ? "],\"rntis_truncated\":true}" : tail_fmt;
    webui_->send(out);
}
#endif
