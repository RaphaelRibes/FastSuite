#!/usr/bin/env python3
"""Does the head_dropout / boundary-anchored change a real-corpus number?

Runs the Phase D per-dataset winning config on the SIMULATED reads (whose
ground truth records the exact artefact window per mate) and scores
exact / overtrim / undertrim, varying only the setting under test.

Inputs on the desktop:
  winners.tsv                      per-dataset winning config from the sweep
  work/<domain>/<acc>/v0/sim.v0_*  the simulated pair + its ground truth

Variants (all inherit the dataset's winning config for every other knob):
  baseline          the shipped default (head_dropout 4)
  hd1               --head-dropout 1
  hd0               --head-dropout 0
  ba                --boundary-anchored
  hd1+ba            both
"""
import argparse, collections, csv, gzip, os, struct, subprocess, sys

W = "/mnt/Fast/fastsuite-scratch/work"
R = "/mnt/Fast/fastsuite-scratch/results"
FP = "/home/raphael/FastSuite/tools/fastprep/target/release/fastprep"


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


def _minlen(v):
    """--min-len takes an integer N (bases) or a fraction written `FL`
    (e.g. `0.9L`); a bare `0.9` is rejected by the parser."""
    f = float(v)
    return f"{v}L" if f < 1.0 else str(int(f))


def eff_hd(cfg):
    """Mirror sweep.py's effective_head_dropout.

    NOTE: the sweep floors this at 1, so its grid could never express
    head_dropout 0 even if a config authored it. That is part of why the
    default went unexamined.
    """
    hd = int(cfg.get("head_dropout") or 4)
    return max(1, min(hd, int(cfg["min_run"]) - 1))


CONF = {
    "min_run": ("--min-run", str), "mercy": ("--mercy", str),
    "phred_sat": ("--phred-sat", str), "gate_shape": ("--gate-shape", str),
    "min_len": ("--min-len", _minlen), "mean_q": ("--mean-q", str),
    "adapter_err": ("--adapter-error-rate", str),
    "min_overlap": ("--min-adapter-overlap", str),
}

# Each variant REPLACES the baseline's trailing `--head-dropout N` pair.
VARIANTS = [
    ("baseline", None),
    ("hd1", 1),
    ("hd0", 0),
    ("ba", None, ["--boundary-anchored"]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--variants", nargs="+", default=[v[0] for v in VARIANTS])
    # shard the corpus across processes: each takes every Nth dataset
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--offset", type=int, default=0)
    a = ap.parse_args()

    winners = list(csv.DictReader(open(f"{R}/winners.tsv"), delimiter="\t"))
    if a.stride > 1:
        winners = winners[a.offset::a.stride]
    if a.limit:
        winners = winners[:a.limit]
    print(f"datasets: {len(winners)}", flush=True)

    if a.offset == 0 and not os.path.exists(a.out):
        with open(a.out, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(
                ["domain", "acc", "variant", "tol", "exact", "overtrim",
                 "undertrim", "dropped", "scored", "cfg_min_run",
                 "cfg_head_dropout"])
    elif a.offset != 0:
        # wait for the header shard to create the file before appending
        import time
        for _ in range(600):
            if os.path.exists(a.out):
                break
            time.sleep(1)

    for i, w in enumerate(winners, 1):
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{W}/{dom}/{acc}/v0"
        gtp = f"{v0}/sim.v0.gt.tsv"
        s1, s2 = f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz"
        if not (os.path.exists(gtp) and os.path.exists(s1)):
            print(f"  [{i}] {acc}: missing sim/gt; skip", flush=True)
            continue
        gt = gt_rows(gtp)
        coremap = core_to_pid([s1, s2])

        # mirror sweep.py build_argv exactly, minus the head-dropout flag
        # (each variant supplies its own), so the baseline reproduces the
        # published fastprep-winner number rather than a relaxed variant.
        base = ["--force", "-q", "--gz-level", "6"]
        for key, (flag, conv) in CONF.items():
            if w.get(key) not in (None, ""):
                base += [flag, conv(w[key])]
        if w.get("r2_reverse") in ("True", "true", "1"):
            base += ["--r2-scan-reverse"]
        if w.get("interior") in ("True", "true", "1"):
            base += ["--interior"]
        base += ["--head-dropout", str(eff_hd(w))]

        for v in VARIANTS:
            name, hd, *rest = v
            if name not in a.variants:
                continue
            extra = list(rest[0]) if rest else []
            # rebuild from the pristine sweep argv each time so variants never
            # inherit each other's head-dropout
            argv = list(base)
            if hd is not None:
                argv = argv[:-2] + ["--head-dropout", str(hd)]
            d = f"/tmp/overtrim_real/{a.offset}_{dom}_{acc}_{name.replace('+', '_')}"
            os.makedirs(d, exist_ok=True)
            px = f"{d}/c.pxmt"
            r = subprocess.run([FP, "-1", s1, "-2", s2, "-o", f"{d}/r1.fq.gz",
                                "-p", f"{d}/r2.fq.gz", "--metadata", px]
                               + argv + extra, capture_output=True)
            if r.returncode != 0:
                print(f"  [{i}] {acc} {name}: rc={r.returncode} "
                      f"{r.stderr.decode()[:160]}", flush=True)
                continue
            p = parse_pxmt(px)
            for tol in (0, 2):
                c = score(p, coremap, gt, tol)
                tot = c["exact"] + c["overtrim"] + c["undertrim"]
                with open(a.out, "a", newline="") as fh:
                    csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(
                        [dom, acc, name, tol, c["exact"], c["overtrim"],
                         c["undertrim"], c["dropped"], tot,
                         w.get("min_run", ""), w.get("head_dropout", "")])
        if i % 10 == 0:
            print(f"  [{i}/{len(winners)}] {acc}", flush=True)
        for name in a.variants:
            subprocess.run(["rm", "-rf",
                            f"/tmp/overtrim_real/{a.offset}_{dom}_{acc}_{name.replace('+', '_')}"],
                           capture_output=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
