#!/usr/bin/env python3
"""Does simply NOT hard-trimming short leading G-runs fix the overtrim?

The claim under test: fastprep hard-trims leading G-runs of length 1..=n
(head_dropout — the shipped default was 4 when this was written; it is now 0,
commit 3d76a31) unconditionally. On fastsim ground truth 64.4% of those r=1
cases are pure biology, so the hard trim overtrims them.

This runs fastsim from a profile, then fastprep in three configurations:
  A  head_dropout 4     (pinned; the then-shipped default)
  B  head dropout off   (--head-dropout 0)
  C  boundary-anchored + head dropout off
and scores exact / overtrim / undertrim against the GT.
"""
import argparse, csv, gzip, os, struct, subprocess, sys, tempfile, collections

FS = "tools/fastsim/target/release/fastsim"
FP = "tools/fastprep/target/release/fastprep"


def gt_rows(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out


def read_fastq(path):
    out, seen = [], collections.Counter()
    if not os.path.exists(path):
        return out
    pid = mate = None
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                f = line[1:].split()[0].split(":")
                x, y = int(f[5]), int(f[6])
                pid = x + 10000 * y
                mate = int(line.split()[1][0])
            elif i % 4 == 1:
                out.append((pid, mate, line.rstrip("\n")))
    return out


def parse_pxmt(path):
    out = {}
    try:
        buf = open(path, "rb").read()
    except FileNotFoundError:
        return out
    if buf[:4] != b"PXMT":
        return out
    off, n = 5, len(buf)
    while off + 2 <= n:
        (idlen,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + idlen].decode("utf-8", "replace").split()[0]; off += idlen
        if off + 22 > n:
            break
        mate = buf[off] + 1; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out[(rid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def core_to_pid(paths):
    out = {}
    for p in paths:
        if not os.path.exists(p):
            continue
        with gzip.open(p, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    core = line[1:].split()[0]
                    f = core.split(":")
                    pid = int(f[5]) + 10000 * int(f[6])
                    out[(core, int(line.split()[1][0]))] = pid
    return out


def score(pxmt, coremap, gt, tol):
    c = collections.Counter(exact=0, overtrim=0, undertrim=0,
                            dropped=0, unmatched=0)
    for (rid, mate), (status, t5, t3) in pxmt.items():
        pid = coremap.get((rid, mate))
        if pid is None or pid not in gt:
            c["unmatched"] += 1
            continue
        h, t = (gt[pid][0], gt[pid][1]) if mate == 1 else (gt[pid][2], gt[pid][3])
        if status != 0:
            c["dropped"] += 1
            continue
        if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
            c["exact"] += 1
        elif t5 > h or t3 > t:
            c["overtrim"] += 1
        else:
            c["undertrim"] += 1
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--pairs", type=int, default=100000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    import json
    prof = json.load(open(a.profile))
    rl = int(prof.get("read_len") or 150)
    ins = prof.get("insert") or {}
    lens = sorted(int(k) for k in ins) if ins else []
    imin = max(1, min(lens[0] if lens else 1, rl))

    work = tempfile.mkdtemp(prefix="overtrim_probe_")
    out = f"{work}/sim"
    subprocess.run([FS, "--pairs", str(a.pairs), "--profile", a.profile,
                    "--out", out, "--seed", "42", "--read-len", str(rl),
                    "--insert-min", str(imin), "--drop-rate", "0", "--force"],
                   check=True, capture_output=True)
    gt = gt_rows(f"{out}.gt.tsv")
    sim1, sim2 = f"{out}_R1.fastq.gz", f"{out}_R2.fastq.gz"
    coremap = core_to_pid([sim1, sim2])
    print(f"profile={os.path.basename(a.profile)} pairs={a.pairs} gt={len(gt)}")

    # NOTE: --head-dropout 0 is now accepted (dfae9b7) and disables the rule.
    # Variant A pins --head-dropout 4 explicitly rather than relying on "the
    # default": 4 was the shipped default when this was written, but the
    # shipped default is now 0 (commit 3d76a31), so an unpinned A would
    # silently collapse into B.
    variants = [
        ("A head-dropout 4", ["--head-dropout", "4"]),
        ("B head-dropout 0 (rule disabled)", ["--head-dropout", "0"]),
        ("C head-dropout 1 (min)", ["--head-dropout", "1"]),
        ("D head-dropout 2", ["--head-dropout", "2"]),
        ("E hd0 + boundary-anchored",
         ["--head-dropout", "0", "--boundary-anchored"]),
        ("F hd0 + boundary-anchored + interior off",
         ["--head-dropout", "0", "--boundary-anchored", "--no-interior"]),
    ]
    rows = []
    for name, extra in variants:
        d = f"{work}/v{len(rows)}"
        os.makedirs(d, exist_ok=True)
        px = f"{d}/c0.pxmt"
        subprocess.run([FP, "-1", sim1, "-2", sim2, "-o", f"{d}/r1.fq.gz",
                        "-p", f"{d}/r2.fq.gz", "--metadata", px, "--no-dedup",
                        "--no-mean-q", "-q"] + extra,
                       check=True, capture_output=True)
        p = parse_pxmt(px)
        for tol in (0, 2):
            c = score(p, coremap, gt, tol)
            tot = c["exact"] + c["overtrim"] + c["undertrim"]
            rows.append((name, tol, c, tot))
            print(f"  {name:38s} tol{tol}: exact {c['exact']/tot*100:6.2f}% "
                  f"over {c['overtrim']/tot*100:6.2f}% "
                  f"under {c['undertrim']/tot*100:6.2f}%  "
                  f"(n={tot} dropped={c['dropped']} unmatch={c['unmatched']})")

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["variant", "tol", "exact", "overtrim", "undertrim",
                    "dropped", "unmatched", "scored"])
        for name, tol, c, tot in rows:
            w.writerow([name, tol, c["exact"], c["overtrim"], c["undertrim"],
                        c["dropped"], c["unmatched"], tot])
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
