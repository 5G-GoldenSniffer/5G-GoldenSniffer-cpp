#!/usr/bin/env python3
"""Dashboard bridge for the NR scanner (stdlib only).

Listens for the per-second JSON datagrams sent by `./main --webui [host:]port` and serves the
dashboard (index.html) plus a Server-Sent Events stream of the datagrams. Browsers can't read
UDP, so this is the piece in between.

    python3 webui/server.py                       # UDP 127.0.0.1:9870, HTTP 127.0.0.1:8080
    ./main ... --webui 9870
    open http://127.0.0.1:8080/

    python3 webui/server.py --http 0.0.0.0:8080   # reachable from other machines (no auth!)

Endpoints: /  (dashboard), /events (SSE: one "history" event with the buffered seconds, then
one "message" per datagram), /api/latest, /api/history (JSON), /a.png /b.png /c.png (default
images of the three image slots, if present next to index.html).
"""
import argparse
import collections
import json
import os
import queue
import socket
import sys
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_IMAGES = ("/a.png", "/b.png", "/c.png")  # optional, next to index.html: default images of the three slots


def hostport(s, default_host):
    host, _, port = s.rpartition(":")
    return (host.strip("[]") or default_host, int(port))


class Hub:
    """History ring buffer + fan-out to the SSE clients."""

    def __init__(self, maxlen):
        self.lock = threading.Lock()
        self.history = collections.deque(maxlen=maxlen)
        self.clients = set()

    def publish(self, line):
        with self.lock:
            self.history.append(line)
            for q in list(self.clients):
                try:
                    q.put_nowait(line)
                except queue.Full:  # stalled client: drop it, it reconnects and gets the history
                    self.clients.discard(q)
                    # the queue is full and this thread holds the lock: make room for the end marker
                    # instead of waiting for a client that isn't reading
                    try:
                        q.get_nowait()
                    except queue.Empty:
                        pass
                    try:
                        q.put_nowait(None)
                    except queue.Full:
                        pass

    def subscribe(self):
        q = queue.Queue(maxsize=300)
        with self.lock:
            self.clients.add(q)
            return q, list(self.history)

    def unsubscribe(self, q):
        with self.lock:
            self.clients.discard(q)


def udp_loop(sock, hub, verbose):
    while True:
        data, addr = sock.recvfrom(65535)
        try:
            msg = json.loads(data)
            if not str(msg.get("schema", "")).startswith("nr-scanner/"):
                raise ValueError("unknown schema")
        except (ValueError, AttributeError) as e:
            print(f"dropped datagram from {addr[0]}:{addr[1]}: {e}", file=sys.stderr)
            continue
        if verbose:
            print(f"+{msg.get('sec')}s {len(data)} bytes", file=sys.stderr)
        hub.publish(json.dumps(msg, separators=(",", ":")))


def make_handler(hub):
    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def log_message(self, fmt, *args):  # quiet
            pass

        def send_body(self, code, ctype, body):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            path = self.path.split("?", 1)[0]
            if path in ("/", "/index.html"):
                with open(os.path.join(HERE, "index.html"), "rb") as f:
                    self.send_body(200, "text/html; charset=utf-8", f.read())
            elif path == "/api/latest":
                with hub.lock:
                    last = hub.history[-1] if hub.history else "null"
                self.send_body(200, "application/json", last.encode())
            elif path == "/api/history":
                with hub.lock:
                    body = "[" + ",".join(hub.history) + "]"
                self.send_body(200, "application/json", body.encode())
            elif path == "/events":
                self.events()
            elif path in DEFAULT_IMAGES and os.path.isfile(os.path.join(HERE, path[1:])):
                with open(os.path.join(HERE, path[1:]), "rb") as f:
                    self.send_body(200, "image/png", f.read())
            else:
                self.send_body(404, "text/plain", b"not found\n")

        def events(self):
            q, history = hub.subscribe()
            try:
                self.send_response(200)
                self.send_header("Content-Type", "text/event-stream")
                self.send_header("Cache-Control", "no-store")
                self.send_header("Connection", "keep-alive")
                self.end_headers()
                self.wfile.write(("retry: 2000\nevent: history\ndata: [" + ",".join(history) + "]\n\n").encode())
                self.wfile.flush()
                while True:
                    try:
                        line = q.get(timeout=15)
                    except queue.Empty:
                        self.wfile.write(b": keepalive\n\n")
                    else:
                        if line is None:
                            break
                        self.wfile.write(("data: " + line + "\n\n").encode())
                    self.wfile.flush()
            except (BrokenPipeError, ConnectionResetError):
                pass
            finally:
                hub.unsubscribe(q)
                self.close_connection = True

    return Handler


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--udp", default="127.0.0.1:9870", help="[host:]port to receive the datagrams on (default %(default)s)")
    ap.add_argument("--http", default="127.0.0.1:8080", help="[host:]port to serve the dashboard on (default %(default)s)")
    ap.add_argument("--history", type=int, default=3600, help="seconds kept for newly opened pages (default %(default)s)")
    ap.add_argument("-v", "--verbose", action="store_true", help="log every datagram")
    a = ap.parse_args()

    hub = Hub(a.history)
    udp_addr = hostport(a.udp, "127.0.0.1")
    sock = socket.socket(socket.AF_INET6 if ":" in udp_addr[0] else socket.AF_INET, socket.SOCK_DGRAM)
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_RCVBUF, 1 << 20)
    sock.bind(udp_addr)
    threading.Thread(target=udp_loop, args=(sock, hub, a.verbose), daemon=True).start()

    http_addr = hostport(a.http, "127.0.0.1")
    ThreadingHTTPServer.daemon_threads = True
    if ":" in http_addr[0]:
        ThreadingHTTPServer.address_family = socket.AF_INET6
    srv = ThreadingHTTPServer(http_addr, make_handler(hub))
    print(f"UDP {udp_addr[0]}:{udp_addr[1]} -> http://{http_addr[0]}:{http_addr[1]}/", file=sys.stderr)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
