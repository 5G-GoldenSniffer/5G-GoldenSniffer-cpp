# Web dashboard

```
./main ... --webui [host:]port  --UDP JSON, 1/s-->  webui/server.py  --HTTP + SSE-->  browser (index.html)
```

`--webui` is in the program only when built with `WEBUI = 1` in `config.mk` (or `make WEBUI=1`).

Browsers can't receive UDP, so `server.py` (Python 3 stdlib only) receives the datagrams, keeps the
last hour, and serves the dashboard plus a Server-Sent Events stream.

```bash
python3 webui/server.py                         # UDP 127.0.0.1:9870 -> http://127.0.0.1:8080/
./main $(cat device-args.example) --freq-rx 3650e6 --webui 9870
python3 webui/server.py --http 0.0.0.0:8080     # expose to the LAN (no authentication)
python3 webui/mock_sender.py --backfill 600     # synthetic data, to work on the UI without a radio
```

`--webui` works with or without `--verbose` (with it the console gets the raw lines and the
summary only goes to the web UI). The scanner never blocks on the socket: if the bridge isn't
running the datagrams are dropped. `index.html` is self-contained (no CDN), so it works offline.

Other endpoints: `/api/latest` (last datagram), `/api/history` (buffered datagrams, JSON array).

## Datagram schema `nr-scanner/1`

One JSON object per second of stream time, the same data as the console summary
(`Summary::send_json` in `summary.cpp`). `null` marks a value that doesn't apply (not tracking,
no SSB in that second, no gain with `--file`, throughput disabled).

| field | meaning |
|---|---|
| `schema` | `"nr-scanner/1"` |
| `t` | wall-clock start of the second, Unix ms (stream start + `sec`; synthetic with `--file`) |
| `sec` | second index since the stream started |
| `partial` | last, incomplete second of the run |
| `cfg` | run configuration, same in every datagram: `source` (USRP model or file path), `file`, `freq_rx_hz`, `bandwidth_hz`, `samp_rate_hz`, `mu`, `scs_hz`, `n_rb`, `fallback_bwp` (RB or null), `throughput` (bool) |
| `cpu` | CPU load over the interval since the previous datagram (null in the first): `system` and `process` (the scanner) in % of all the cores, `cores` per-core load in % (`WebUI::cpu_json`, from `/proc/stat`) |
| `rx` | receiver status (null if the RX tick is missing): `state` (`AGC`/`ACQUIRE`/`TRACK`, at the end of the second), `pci`, `sfn` (from the MIB, null until decoded), `f_ssb_hz`, `ssb_found`/`ssb_missed` (frames), `snr_db`, `cfo_hz`, `sfo_ppm`, `gain_db`, `frames_shipped`, `frames_dropped` (live only: frames skipped because the workers were behind), `overflows`, `acquisitions`, `track_losses`, `acq_failures` (counts in this second; both a tracking loss and a failed acquisition send the receiver back to `AGC`) |
| `frames` | frames processed in this second (10 ms each) |
| `slots`, `pdcch_slots` | slots processed / slots with at least one DCI |
| `dci` | `{dl, ul}` DCI counts (first DCI bit 1 = DL; SI-RNTI and P-RNTI DCIs are always DL) |
| `throughput` | `{dl_bps, ul_bps}` over the processed frames, or null without `--estimate-throughput` |
| `rnti_count` | RNTIs seen in this second |
| `rntis` | busiest first: `{rnti, nids: [DM-RS scrambling IDs], dl, ul, dl_bps?, ul_bps?}` (`*_bps` only with throughput) |
| `rntis_truncated` | the list was cut to keep the datagram under 60000 bytes; `rnti_count` is still exact |

A breaking change bumps the schema version; new fields may be added to `nr-scanner/1`.
