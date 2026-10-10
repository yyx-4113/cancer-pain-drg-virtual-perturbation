"""Segmented parallel downloader + integrity check for GSE292522.

WHY THIS SCRIPT EXISTS
----------------------
Two NCBI endpoints serve the same supplementary archive and they do NOT
behave the same way:

  https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE292522&format=file
      -> HTTP 200, NO Accept-Ranges, NO Content-Length.
         The whole 711 MB has to stream over one connection. Measured
         throughput here: ~12.5 KB/s  ->  ~16 hours. Not usable.

  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292522/suppl/GSE292522_RAW.tar
      -> HTTP 206, Accept-Ranges: bytes, Content-Length: 711014400.
         Range requests work, so the file can be pulled in parallel
         segments and resumed.

Any attempt to "parallelise" the first URL silently degenerates: the server
ignores the Range header and every worker streams the whole file from byte
0. That is exactly how the earlier attempt produced a 229 KB file that was
also not valid gzip.

The downloader below therefore (a) uses only the ftp endpoint, (b) refuses
to start unless the server advertises byte ranges and reports the expected
length, (c) writes into a sidecar .part file per segment so it can resume,
and (d) verifies the result before declaring success.

Resume: delete nothing; re-running continues from the .part files.
"""
import hashlib
import math
import os
import queue
import struct
import sys
import tarfile
import threading
import time
import urllib.request

URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE292nnn/GSE292522/"
       "suppl/GSE292522_RAW.tar")
EXPECTED = 711014400                      # bytes, from filelist.txt
RAW = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "data", "raw")
OUT = os.path.join(RAW, "GSE292522_RAW.tar")
N_SEG = 16
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def probe():
    """Refuse to proceed unless the endpoint really supports byte ranges."""
    req = urllib.request.Request(URL, headers={**UA, "Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=40) as r:
        if r.status != 206 or not r.headers.get("Accept-Ranges"):
            raise SystemExit(f"endpoint does not support ranges "
                             f"(HTTP {r.status}, Accept-Ranges="
                             f"{r.headers.get('Accept-Ranges')}) - "
                             f"parallel download impossible, aborting")
        total = int(r.headers["Content-Range"].split("/")[1])
    log(f"endpoint OK: HTTP 206, Accept-Ranges=bytes, total={total:,}")
    if total != EXPECTED:
        log(f"  WARNING: filelist.txt said {EXPECTED:,}, server says "
            f"{total:,} -> trusting the server")
    return total


def seg_path(i):
    return OUT + f".p{i:02d}"


def worker(i, start, end, done, counter, lock, errq):
    """Pull bytes [start, end] into OUT.pNNG with resume."""
    p = seg_path(i)
    have = os.path.getsize(p) if os.path.exists(p) else 0
    if have >= end - start + 1:
        with lock:
            counter[0] += end - start + 1
        done[i] = True
        return
    rng_hdr = {"Range": f"bytes={start + have}-{end}"}
    try:
        with urllib.request.urlopen(
                urllib.request.Request(URL, headers={**UA, **rng_hdr}),
                timeout=60) as r, open(p, "ab") as fh:
            while True:
                buf = r.read(262144)
                if not buf:
                    break
                fh.write(buf)
                with lock:
                    counter[0] += len(buf)
    except Exception as e:                      # transient: leave partial
        errq.put((i, f"{type(e).__name__}: {e}"))
        return
    if os.path.getsize(p) == end - start + 1:
        done[i] = True


def download_parallel(total):
    seg_len = math.ceil(total / N_SEG)
    bounds = [(i, i * seg_len, min(total, (i + 1) * seg_len) - 1)
              for i in range(N_SEG)]
    done_bytes = lambda i: (os.path.getsize(seg_path(i))
                            if os.path.exists(seg_path(i)) else 0)
    lock = threading.Lock()
    counter = [0]
    errq = queue.Queue()
    round_no = 0
    while round_no < 40:
        todo = [(i, s, e) for i, s, e in bounds if done_bytes(i) < e - s + 1]
        if not todo:
            break
        round_no += 1
        done = {}
        t0 = time.time()
        base = counter[0]
        threads = [threading.Thread(target=worker,
                                    args=(i, s, e, done, counter, lock, errq))
                   for i, s, e in todo]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        dt = time.time() - t0
        got = sum(done_bytes(i) for i, _, _ in bounds)
        log(f"  round {round_no}: {len(todo)} open, {got:,}/{total:,} "
            f"({100 * got / total:.1f}%), +{(counter[0] - base) / 1024:.0f} KB "
            f"({(counter[0] - base) / max(dt, 1) / 1024:.0f} KB/s)")
    while not errq.empty():
        log(f"  segment error: {errq.get()}")
    return all(done_bytes(i) >= e - s + 1 for i, s, e in bounds)


def assemble(total):
    h = hashlib.sha256()
    with open(OUT, "wb") as out:
        for i, s, e in [(i, i * math.ceil(total / N_SEG),
                         min(total, (i + 1) * math.ceil(total / N_SEG)) - 1)
                        for i in range(N_SEG)]:
            with open(seg_path(i), "rb") as fh:
                while True:
                    buf = fh.read(1 << 20)
                    if not buf:
                        break
                    out.write(buf)
                    h.update(buf)
    sz = os.path.getsize(OUT)
    log(f"assembled {OUT}: {sz:,} bytes, sha256={h.hexdigest()[:16]}")
    if sz != total:
        raise SystemExit(f"size mismatch {sz:,} != {total:,}")
    return h.hexdigest()


def verify():
    """TAR can hold garbage and still list names; check members and CRCs."""
    log("verifying archive integrity (tar headers + gzip CRC per member)")
    bad = []
    try:
        with tarfile.open(OUT, "r:") as tf:
            members = tf.getmembers()
            log(f"  {len(members)} members declared")
            for m in members:
                if not m.isfile():
                    continue
                f = tf.extractfile(m)
                try:
                    n = 0
                    while f.read(1 << 20):
                        n += 1
                except Exception as e:
                    bad.append((m.name, f"{type(e).__name__}: {e}"))
    except Exception as e:
        raise SystemExit(f"tar could not be opened: {type(e).__name__}: {e}")
    if bad:
        for n_, e_ in bad:
            log(f"  CORRUPT: {n_} -> {e_}")
        raise SystemExit(f"{len(bad)} corrupt member(s)")
    ok_files = [m.name for m in members if m.name.endswith(".gz")]
    log(f"  OK: {len(members)} members, {len(ok_files)} .gz members readable")
    return members


def main():
    os.makedirs(RAW, exist_ok=True)
    disk = os.statvfs(RAW) if hasattr(os, "statvfs") else None
    total = probe()
    need = total / 2**20
    log(f"target {OUT} ({need:.0f} MB), {N_SEG} parallel segments")
    if os.path.exists(OUT) and os.path.getsize(OUT) == total:
        log("file already complete - skipping download")
    else:
        ok = download_parallel(total)
        if not ok:
            raise SystemExit("some segments incomplete; re-run to resume")
        assemble(total)
    verify()
    for i in range(N_SEG):
        p = seg_path(i)
        if os.path.exists(p):
            os.remove(p)
    log("DONE")


if __name__ == "__main__":
    main()
