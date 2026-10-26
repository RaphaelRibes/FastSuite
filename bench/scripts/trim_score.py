#!/usr/bin/env python3
"""Score fastprep's trim against the simulator's ground truth (no fastp needed).

Committed 2026-10-26 with the post-fix re-run (bench/results/fix2026/).
TRIM-ACCURACY.md cites `bench/scripts/overtrim_vs_fastp.py` as the producer of
the four-configuration table; that file has never been in this repository, so
this scorer was written to re-derive the fastprep rows from the committed GT.
It scores one configuration at a time (the binary's defaults or a wrapper).

For each dataset: read the GT (pair_id, r1_h, r1_t, r2_h, r2_t, r1_adapter,
r2_adapter), run fastprep on the sim reads, then join output reads back by
header and infer the observed cut from the emitted sequence:

    out_seq == in_seq[t5 : L - t3]   ->   t5 = offset, t3 = L - offset - len(out)

A read whose observed (t5,t3) equals the drawn (h+adapter, t) within tol is
exact; more trimmed is overtrim, less is undertrim. A dropped read is
`dropped_or_missing`.
"""
import argparse, csv, gzip, os, subprocess, sys, tempfile
from collections import Counter

def read_fastq(path):
    """yield (name_without_mate_suffix, seq) for a (possibly gzipped) FASTQ."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            yield h[1:].rstrip("\n"), s

def key(name):
    # strip the mate suffix fastprep may add ("@id/1") and take the first token
    n = name.split()[0]
    return n.rsplit("/", 1)[0] if n.endswith(("/1", "/2")) else n

def load_gt(path):
    """GT rows in file order; the simulator writes pair_id 0..N-1 in that order."""
    gt = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gt.append({k: int(row[k]) for k in
                       ("r1_h", "r1_t", "r2_h", "r2_t", "r1_adapter", "r2_adapter")})
    return gt

def find_offset(seq, out):
    """Return t5 if out is a contiguous substring of seq (prefix-anchored search)."""
    i = seq.find(out)
    return i if i >= 0 else None

def score(dataset, fastprep, tol):
    work, domain, acc = dataset
    v0 = os.path.join(work, domain, acc, "v0")
    gt_path = os.path.join(v0, "sim.v0.gt.tsv")
    r1 = os.path.join(v0, "sim.v0_R1.fastq.gz")
    r2 = os.path.join(v0, "sim.v0_R2.fastq.gz")
    if not (os.path.exists(gt_path) and os.path.exists(r1)):
        return None
    gt = load_gt(gt_path)
    # The GT is row-ordered, the FASTQ is Illumina-named: join by ORDER, using
    # each output header to find its input index (fastprep preserves headers and
    # input order for the reads it keeps).
    r1rec = list(read_fastq(r1))
    r2rec = list(read_fastq(r2))
    ins = {1: [s for _, s in r1rec], 2: [s for _, s in r2rec]}
    idx = {"1": {key(n): i for i, (n, _) in enumerate(r1rec)},
           "2": {key(n): i for i, (n, _) in enumerate(r2rec)}}
    counts = Counter()
    with tempfile.TemporaryDirectory() as d:
        o1, o2 = os.path.join(d, "o1.fq.gz"), os.path.join(d, "o2.fq.gz")
        rep = os.path.join(d, "r.json")
        cp = subprocess.run([fastprep, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
                             "--json", rep, "-q", "--force"],
                            capture_output=True, text=True)
        if cp.returncode != 0 and not os.path.exists(o1):
            return None
        outs = {}
        for mate, p in (("1", o1), ("2", o2)):
            if os.path.exists(p):
                for n, s in read_fastq(p):
                    k = idx[mate].get(key(n))
                    if k is not None:
                        outs[(k, mate)] = s
    for i, g in enumerate(gt):
        for m, hcol, tcol, acol in ((1, "r1_h", "r1_t", "r1_adapter"),
                                    (2, "r2_h", "r2_t", "r2_adapter")):
            seq = ins[m][i] if i < len(ins[m]) else None
            if seq is None:
                continue
            mate = str(m)
            want5 = g[hcol] + g[acol]
            want3 = g[tcol]
            out = outs.get((i, mate))
            if out is None:
                counts["dropped_or_missing"] += 1
                continue
            if len(out) > len(seq):
                counts["unaligned"] += 1
                continue
            off = find_offset(seq, out)
            if off is None:
                counts["unaligned"] += 1
                continue
            t3 = len(seq) - off - len(out)
            if abs(off - want5) <= tol and abs(t3 - want3) <= tol:
                counts["exact"] += 1
            elif off > want5 or t3 > want3:
                counts["overtrim"] += 1
            else:
                counts["undertrim"] += 1
    return counts

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--fastprep", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=int, default=0)
    args = ap.parse_args()

    datasets = []
    for domain in sorted(os.listdir(args.work)):
        dd = os.path.join(args.work, domain)
        if not os.path.isdir(dd):
            continue
        for acc in sorted(os.listdir(dd)):
            datasets.append((args.work, domain, acc))

    total = Counter()
    rows = []
    for i, ds in enumerate(datasets, 1):
        c = score(ds, args.fastprep, args.tol)
        if c is None:
            print(f"  {ds[2]}: no GT / failed", flush=True)
            continue
        total.update(c)
        rows.append((ds[1], ds[2], c["exact"], c["overtrim"], c["undertrim"],
                     c["unaligned"], c["dropped_or_missing"]))
        print(f"  [{i}/{len(datasets)}] {ds[1]}/{ds[2]}: exact {c['exact']} "
              f"over {c['overtrim']} under {c['undertrim']} dropped {c['dropped_or_missing']}",
              flush=True)

    n = sum(total.values())
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["domain", "run_accession", "exact", "overtrim", "undertrim",
                    "unaligned", "dropped_or_missing"])
        for r in rows:
            w.writerow(r)
    print(f"\nTOTAL mates {n}")
    for k in ("exact", "overtrim", "undertrim", "unaligned", "dropped_or_missing"):
        print(f"  {k:18s} {total[k]:8d}  {100.0*total[k]/n if n else 0:.2f}%")
    scored = n - total["dropped_or_missing"] - total["unaligned"]
    print(f"  exact/scored       {100.0*total['exact']/scored if scored else 0:.2f}%  (denominator {scored})")

if __name__ == "__main__":
    main()
