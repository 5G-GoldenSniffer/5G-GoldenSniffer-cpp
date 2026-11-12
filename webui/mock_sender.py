#!/usr/bin/env python3
"""Synthetic nr-scanner/1 datagrams, to work on the dashboard without a USRP.

    python3 webui/mock_sender.py [--to 127.0.0.1:9870] [--backfill 600]
"""
import argparse, json, math, random, socket, time

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--to", default="127.0.0.1:9870")
ap.add_argument("--backfill", type=int, default=0, help="send this many past seconds at once first")
a = ap.parse_args()
host, _, port = a.to.rpartition(":")
sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
dest = (host or "127.0.0.1", int(port))

cfg = {"source": "n310", "file": False, "freq_rx_hz": 3650000000, "bandwidth_hz": 30000000, "samp_rate_hz": 46080000,
       "mu": 1, "scs_hz": 30000, "n_rb": 78, "fallback_bwp": 51, "throughput": True}
ues = {}
start = time.time() - a.backfill
state, cfo = "AGC", 400.0

def second(sec):
    global state, cfo
    rx = {"state": state, "pci": None, "f_ssb_hz": None, "ssb_found": 0, "ssb_missed": 0, "snr_db": None, "cfo_hz": None,
          "sfo_ppm": None, "gain_db": 42.0, "frames_shipped": 0, "frames_dropped": 0, "overflows": 0, "acquisitions": 0, "track_losses": 0, "acq_failures": 0}
    msg = {"schema": "nr-scanner/1", "t": int((start + sec) * 1000), "sec": sec, "partial": False, "cfg": cfg, "rx": rx,
           "frames": 0, "slots": 0, "pdcch_slots": 0, "dci": {"dl": 0, "ul": 0}, "throughput": None, "rnti_count": 0, "rntis": []}
    cores = [min(100, max(0, round(random.gauss(85 if i in (0, 2) else 35 if i < 16 else 25, 6)))) for i in range(32)]
    msg["cpu"] = {"system": round(sum(cores) / 32, 1), "process": round(sum(cores) / 32 - 3, 1), "cores": cores}
    if state == "AGC":
        if sec >= 1: state = "ACQUIRE"
        return msg
    if state == "ACQUIRE":
        if random.random() < 0.6:
            state = "TRACK"; rx["acquisitions"] = 1
        else:
            return msg
    rx["state"] = "TRACK"
    if random.random() < 0.004:  # lost during the second: reported in ACQUIRE
        state = rx["state"] = "ACQUIRE"; rx["track_losses"] = 1
    load = 0.5 + 0.45 * math.sin(sec / 90)
    frames = 100
    found = random.randint(48, 50)
    cfo += random.gauss(0, 3)
    rx.update(pci=271, f_ssb_hz=3649440000.0, ssb_found=found, ssb_missed=frames - found,
              snr_db=round(18 + 3 * math.sin(sec / 40) + random.gauss(0, 0.6), 2), cfo_hz=round(cfo, 1),
              sfo_ppm=round(-0.42 + random.gauss(0, 0.02), 4), frames_shipped=frames,
              overflows=1 if random.random() < 0.01 else 0)
    # UEs come and go
    for r in list(ues):
        if random.random() < 0.02: del ues[r]
    while len(ues) < int(3 + 25 * load):
        ues[random.randint(0x4601, 0x4fff)] = (random.choice([1, 2, 3, 4]), random.choice([0x03ef, 0x0123]))
    rntis, dl_t, ul_t, dl_b, ul_b = [], 0, 0, 0.0, 0.0
    for r, (w, nid) in ues.items():
        dl = int(random.expovariate(1 / (w * 20 * load))); ul = int(random.expovariate(1 / (w * 12 * load)))
        if dl + ul == 0: continue
        db, ub = dl * random.uniform(800, 9000), ul * random.uniform(300, 4000)
        rntis.append({"rnti": r, "nids": [nid], "dl": dl, "ul": ul, "dl_bps": round(db), "ul_bps": round(ub)})
        dl_t += dl; ul_t += ul; dl_b += db; ul_b += ub
    si = random.randint(8, 12)
    rntis.append({"rnti": 0xffff, "nids": [271], "dl": si, "ul": 0, "dl_bps": si * 2000, "ul_bps": 0})
    dl_t += si; dl_b += si * 2000
    rntis.sort(key=lambda x: (-(x["dl"] + x["ul"]), x["rnti"]))
    msg.update(frames=frames, slots=frames * 20, pdcch_slots=min(frames * 20, int((dl_t + ul_t) * 0.8)),
               dci={"dl": dl_t, "ul": ul_t}, throughput={"dl_bps": round(dl_b), "ul_bps": round(ul_b)},
               rnti_count=len(rntis), rntis=rntis)
    return msg

sec = 0
for sec in range(a.backfill):
    sock.sendto(json.dumps(second(sec)).encode(), dest)
    time.sleep(0.002)
sec = a.backfill
while True:
    sock.sendto(json.dumps(second(sec)).encode(), dest)
    sec += 1
    time.sleep(max(0, start + sec - time.time()))
