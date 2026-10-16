#!/usr/bin/env python3
"""Per-mate breakdown of the leading-G hard trim's cost.

The profiles carry head artefact on ONE mate only (mid: R2 has 18.5% h=1,
R1 has none), while leading biological G's appear on both. An unconditional
rule therefore overtrims R1 and correctly trims R2. If that asymmetry is
general, a PER-MATE rule beats any global setting.
"""
import csv, os, subprocess, sys, tempfile, json, collections, struct, gzip
sys.path.insert(0, "bench/scripts")
import overtrim_head_dropout as H

FS = "tools/fastsim/target/release/fastsim"
FP = "tools/fastprep/target/release/fastprep"


def score_by_mate(pxmt, coremap, gt, tol):
    c = collections.defaultdict(lambda: collections.Counter(
        exact=0, overtrim=0, undertrim=0, dropped=0))
    for (rid, mate), (status, t5, t3) in pxmt.items():
        pid = coremap.get((rid, mate))
        if pid is None or pid not in gt:
            continue
        h, t = (gt[pid][0], gt[pid][1]) if mate == 1 else (gt[pid][2], gt[pid][3])
        d = c[mate]
        if status != 0:
            d["dropped"] += 1
            continue
        if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
            d["exact"] += 1
        elif t5 > h or t3 > t:
            d["overtrim"] += 1
        else:
            d["undertrim"] += 1
    return c


def run(prof_path, label, pairs=30000):
    prof = json.load(open(prof_path))
    rl = int(prof.get("read_len") or 150)
    ins = prof.get("insert") or {}
    lens = sorted(int(k) for k in ins) if ins else []
    imin = max(1, min(lens[0] if lens else 1, rl))
    work = tempfile.mkdtemp(prefix="pm_")
    out = f"{work}/sim"
    subprocess.run([FS, "--pairs", str(pairs), "--profile", prof_path,
                    "--out", out, "--seed", "42", "--read-len", str(rl),
                    "--insert-min", str(imin), "--drop-rate", "0", "--force"],
                   check=True, capture_output=True)
    gt = H.gt_rows(f"{out}.gt.tsv")
    sim1, sim2 = f"{out}_R1.fastq.gz", f"{out}_R2.fastq.gz"
    coremap = H.core_to_pid([sim1, sim2])
    print(f"\n=== {label}  (head artefact R1={prof['head'].get('R1')} "
          f"R2={prof['head'].get('R2')})")
    for vi, (name, extra) in enumerate([
            ("head-dropout 4 (default)", []),
            ("head-dropout 1", ["--head-dropout", "1"]),
            ("head-dropout 2", ["--head-dropout", "2"]),
            ("head-dropout 0 (off)", ["--head-dropout", "0"])]):
        d = f"{work}/v{vi}"
        os.makedirs(d, exist_ok=True)
        px = f"{d}/c0.pxmt"
        r = subprocess.run([FP, "-1", sim1, "-2", sim2, "-o", f"{d}/r1.fq.gz",
                            "-p", f"{d}/r2.fq.gz", "--metadata", px, "--no-dedup",
                            "--no-mean-q", "-q"] + extra,
                           capture_output=True)
        if r.returncode != 0:
            print(f"  !! {name} rc={r.returncode}: "
                  f"{r.stderr.decode()[:300]}")
            continue
        c = score_by_mate(H.parse_pxmt(px), coremap, gt, 0)
        for mate in (1, 2):
            d2 = c[mate]
            tot = d2["exact"] + d2["overtrim"] + d2["undertrim"]
            if not tot:
                continue
            print(f"  {name:24s} R{mate}: exact {d2['exact']/tot*100:6.2f}% "
                  f"over {d2['overtrim']/tot*100:6.2f}% "
                  f"under {d2['undertrim']/tot*100:6.2f}%  (n={tot})")


if __name__ == "__main__":
    run("bench/fixtures/overtrim/mid.json", "mid")
    run("bench/fixtures/overtrim/low.json", "low")
