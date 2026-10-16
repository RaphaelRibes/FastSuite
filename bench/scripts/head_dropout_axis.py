#!/usr/bin/env python3
"""Does the extended head_dropout axis change the per-dataset winners?

The old sweep could never express head_dropout 0 (effective_head_dropout
floored the value at 1) and never emitted --head-dropout at all, so every
winning config ran at whatever the CLI default was then — 4 (the shipped
value at the time; it is now 0, commit 3d76a31). This re-scores each
dataset's WINNING config from winners.tsv with ONLY head_dropout overridden
(0..4) and scores the PXMT trim decisions against the sim ground truth per
mate at tolerance 0 and 2 — the exact / overtrim / undertrim convention of
overtrim_real.py and the paper's tables. All three sim variants (v0/v1/v2)
are scored where they exist.

One TSV row per (dataset, sim variant, head_dropout, tolerance):
  domain acc variant head_dropout tol exact overtrim undertrim dropped scored
  cfg_min_run cfg_config_id read_len

`head_dropout` is the value that actually RAN after the min_run clamp — the
same invariant as sweep.py: the recorded column is what the argv contained,
never what the config authored. min_len is passed as the winner's recorded
min_len_effective (the exact arg the winning run used); read_len is recorded
per row so a cross-variant read-length mismatch would be visible.

Each invocation writes its OWN output file (header included); shards never
share an append. Merge with a header-once concat.

Usage (on the desktop, where the corpus lives):
  python3 head_dropout_axis.py --out /tmp/hd_axis/hd.0.tsv --stride 8 --offset 0
"""
import argparse, collections, csv, gzip, os, struct, subprocess, sys

W = "/mnt/Fast/fastsuite-scratch/work"
R = "/mnt/Fast/fastsuite-scratch/results"
FP = "/home/raphael/FastSuite/tools/fastprep/target/release/fastprep"

FIELDS = ["domain", "acc", "variant", "head_dropout", "tol", "exact",
          "overtrim", "undertrim", "dropped", "scored", "cfg_min_run",
          "cfg_config_id", "read_len"]


def gt_rows(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
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
        rid = buf[off:off + idlen].decode("utf-8", "replace").split()[0]
        off += idlen
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


def score(pxmt, coremap, gt, tol):
    c = collections.Counter(exact=0, overtrim=0, undertrim=0, dropped=0)
    for (rid, mate), (status, t5, t3) in pxmt.items():
        pid = coremap.get((rid, mate))
        if pid is None or pid not in gt:
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


def eff_hd(hd, min_run):
    """Mirror sweep.py's effective_head_dropout after the floor-0 fix.

    fastprep refuses head_dropout >= min_run; 0 (rule off) is a real value and
    must survive. The recorded column is this post-clamp value.
    """
    return max(0, min(int(hd), int(min_run) - 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--winners", default=f"{R}/winners.tsv")
    ap.add_argument("--fastprep", default=FP)
    ap.add_argument("--variants", nargs="+", default=["v0", "v1", "v2"])
    ap.add_argument("--hd", nargs="+", type=int, default=[0, 1, 2, 3, 4])
    ap.add_argument("--limit", type=int, default=0)
    # shard the corpus across processes: each takes every Nth dataset
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--offset", type=int, default=0)
    a = ap.parse_args()

    winners = list(csv.DictReader(open(a.winners), delimiter="\t"))
    if a.stride > 1:
        winners = winners[a.offset::a.stride]
    if a.limit:
        winners = winners[:a.limit]
    print(f"datasets: {len(winners)}", flush=True)

    with open(a.out, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(FIELDS)

    for i, w in enumerate(winners, 1):
        dom, acc = w["domain"], w["run_accession"]
        mr = int(w["min_run"])
        for var in a.variants:
            d0 = f"{W}/{dom}/{acc}/{var}"
            gtp = f"{d0}/sim.{var}.gt.tsv"
            s1, s2 = (f"{d0}/sim.{var}_R1.fastq.gz", f"{d0}/sim.{var}_R2.fastq.gz")
            if not (os.path.exists(gtp) and os.path.exists(s1)
                    and os.path.exists(s2)):
                print(f"  [{i}] {acc} {var}: missing sim/gt; skip", flush=True)
                continue
            gt = gt_rows(gtp)
            coremap = core_to_pid([s1, s2])
            import gzip as _gz
            read_len = 0
            with _gz.open(s1, "rt") as fh:
                for line in fh:
                    if not line.startswith("@"):
                        read_len = len(line.rstrip("\n"))
                        break

            # mirror sweep.py build_argv exactly, with the winner's own
            # min_len_effective and the head_dropout override
            base = ["--force", "-q", "--gz-level", "6",
                    "--min-run", str(mr),
                    "--mercy", str(w["mercy"]),
                    "--phred-sat", str(w["phred_sat"]),
                    "--gate-shape", str(w["gate_shape"]),
                    "--min-len", str(w["min_len_effective"]),
                    "--mean-q", str(w["mean_q"]),
                    "--adapter-error-rate", str(w["adapter_err"]),
                    "--min-adapter-overlap", str(w["min_overlap"])]
            if w.get("r2_reverse") in ("True", "true", "1"):
                base += ["--r2-scan-reverse"]
            if w.get("interior") in ("True", "true", "1"):
                base += ["--interior"]

            for hd in a.hd:
                hd_eff = eff_hd(hd, mr)
                argv = base + ["--head-dropout", str(hd_eff)]
                d = (f"/tmp/hd_axis/{a.offset}_{dom}_{acc}_{var}_hd{hd}")
                os.makedirs(d, exist_ok=True)
                px = f"{d}/c.pxmt"
                r = subprocess.run([a.fastprep, "-1", s1, "-2", s2,
                                    "-o", f"{d}/r1.fq.gz", "-p", f"{d}/r2.fq.gz",
                                    "--metadata", px] + argv,
                                   capture_output=True)
                if r.returncode != 0:
                    print(f"  [{i}] {acc} {var} hd={hd}: rc={r.returncode} "
                          f"{r.stderr.decode()[:160]}", flush=True)
                    continue
                p = parse_pxmt(px)
                for tol in (0, 2):
                    c = score(p, coremap, gt, tol)
                    tot = c["exact"] + c["overtrim"] + c["undertrim"]
                    with open(a.out, "a", newline="") as fh:
                        csv.writer(fh, delimiter="\t",
                                   lineterminator="\n").writerow(
                            [dom, acc, var, hd_eff, tol, c["exact"],
                             c["overtrim"], c["undertrim"], c["dropped"],
                             tot, mr, w.get("config_id", ""), read_len])
                subprocess.run(["rm", "-rf", d], capture_output=True)
        if i % 10 == 0:
            print(f"  [{i}/{len(winners)}] {acc}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
