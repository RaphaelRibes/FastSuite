#!/usr/bin/env python3
"""Two-axis A/B on the adapter axis: undertrim vs overtrim on one corpus.

Builds a single-end corpus with a realistic insert distribution straddling the
read length, so a known share of reads carries 1..33 bases of adapter
readthrough, and scores fastprep's output per read against the exact truth:

    exact   cut == t
    under   cut < t          (residual adapter left in the output read)
    over    cut > t          (real insert removed)

Both axes come from the same run, which is the point: a knob that buys
undertrim without its overtrim price quoted is not a measurement.

Usage
-----
    python3 bench/scripts/undertrim_ab.py [--reads 200000] [--out DIR]
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from undertrim_probe import TRUSEQ_R1, dna, parse_pxmt, write_fastq  # noqa: E402

L = 150


def build_corpus(rng, n):
    """-> records, {id: (t, cls)} with cls in clean/plain/g/n/sub"""
    recs, truth = [], {}
    for k in range(n):
        rid = f"Q{k}"
        # 70% of fragments are longer than the read (no readthrough at all);
        # 30% fall short and carry readthrough, uniformly spread over 1..33.
        if rng.random() < 0.70:
            t, cls = 0, "clean"
        else:
            t = rng.randint(1, 33)
            cls = "plain"
        body = dna(rng, L - t)
        tail = list(TRUSEQ_R1[:t])
        if t and cls == "plain":
            r = rng.random()
            if r < 0.10:                      # NovaSeq no-signal: last bases G
                g = rng.randint(3, 8)
                for j in range(max(0, len(tail) - g), len(tail)):
                    tail[j] = "G"
                cls = "g"
            elif r < 0.15:                    # a couple of ambiguous calls
                for j in rng.sample(range(len(tail)), min(2, len(tail))):
                    tail[j] = "N"
                cls = "n"
            elif r < 0.18:                    # substitution errors in the adapter
                for _ in range(rng.randint(3, 4)):
                    j = rng.randrange(len(tail))
                    tail[j] = rng.choice([b for b in "ACGT" if b != tail[j]])
                cls = "sub"
        recs.append((rid, body + "".join(tail)))
        truth[rid] = (t, cls)
    return recs, truth


CONFIGS = [
    ("shipped (mo=14, eps=0.2)", []),
    ("mo=12", ["--min-adapter-overlap", "12"]),
    ("mo=10", ["--min-adapter-overlap", "10"]),
    ("mo=8", ["--min-adapter-overlap", "8"]),
    ("mo=10, eps=0.3", ["--min-adapter-overlap", "10", "--adapter-error-rate", "0.3"]),
    ("mo=14, eps=0.3", ["--adapter-error-rate", "0.3"]),
    ("mo=6", ["--min-adapter-overlap", "6"]),
]


def run(fp, d, extra, tag):
    o = os.path.join(d, f"o_{tag}.fq")
    px = os.path.join(d, f"o_{tag}.pxmt")
    for p in (o, px):
        if os.path.exists(p):
            os.remove(p)
    cmd = [fp, "-1", os.path.join(d, "in.fq"), "-o", o, "--metadata", px,
           "--adapter-preset", "illumina", "--no-dedup", "--no-mean-q",
           "--min-len", "0", "-j", "1", "-q", *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"fastprep failed: {r.stderr[-1200:]}")
    return parse_pxmt(px)


def score(px, truth):
    out = {}
    for rid, (t, cls) in truth.items():
        st, t5a, t3a, t5p, t3p, _ = px[(rid, 0)]
        cut = t3a + t3p
        if st != 0:
            b = "dropped"
        elif cut == t:
            b = "exact"
        elif cut < t:
            b = "under"
        else:
            b = "over"
        out.setdefault(cls, {}).setdefault(b, 0)
        out[cls][b] += 1
        out.setdefault("ALL", {}).setdefault(b, 0)
        out["ALL"][b] += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--reads", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=31)
    ap.add_argument("--out", default="bench/results/undertrim-ab.json")
    a = ap.parse_args()
    fp = os.path.abspath(a.fastprep)
    rng = random.Random(a.seed)
    d = tempfile.mkdtemp(prefix="utab_")
    recs, truth = build_corpus(rng, a.reads)
    write_fastq(os.path.join(d, "in.fq"), recs)
    n_clean = sum(1 for t, c in truth.values() if c == "clean")
    print(f"corpus: {a.reads} reads, {n_clean} adapter-free, "
          f"{a.reads - n_clean} with readthrough "
          f"({100*(a.reads-n_clean)/a.reads:.1f}%)\n")
    print(f"{'config':26s} {'exact':>7} {'under':>7} {'over':>6} | "
          f"{'und/carrier':>11} {'false/cmp':>10} | g under | n under | sub under")
    report = {"reads": a.reads, "clean": n_clean, "configs": {}}
    for name, extra in CONFIGS:
        px = run(fp, d, extra, name.split()[0].replace("=", "").replace(",", ""))
        s = score(px, truth)
        carriers = a.reads - n_clean
        und = s["ALL"].get("under", 0)
        fp_clean = s["clean"].get("over", 0) + s["clean"].get("under", 0)
        row = dict(name=name, extra=extra, buckets=s,
                   undertrim_rate=und / carriers,
                   clean_false_trim_rate=fp_clean / n_clean)
        report["configs"][name] = row
        print(f"{name:26s} {s['ALL'].get('exact',0):7d} {und:7d} "
              f"{s['ALL'].get('over',0):6d} | {100*und/carriers:10.2f}% "
              f"{100*fp_clean/n_clean:9.3f}% | "
              f"{s.get('g',{}).get('under',0):7d} | {s.get('n',{}).get('under',0):7d} | "
              f"{s.get('sub',{}).get('under',0):7d}")
    with open(a.out, "w") as f:
        json.dump(report, f, indent=1)
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
