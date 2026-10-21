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
from undertrim_probe import (  # noqa: E402
    TRUSEQ_R1, TRUSEQ_R2, dna, parse_pxmt, write_fastq,
)

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]

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
    # UT-3: the corroborated short tail. Factor 50-100 is the measured sweet
    # spot on this corpus; factor 5 is recorded as the arm that FAILS the
    # false-trim gate, because that contrast is the whole argument for the
    # corroboration.
    ("UT-3 st=6, enr=50", ["--adapter-short-tail", "6", "--adapter-min-enrichment", "50"]),
    ("UT-3 st=6, enr=100", ["--adapter-short-tail", "6", "--adapter-min-enrichment", "100"]),
    ("UT-3 st=6, enr=5 (fails the gate)", ["--adapter-short-tail", "6", "--adapter-min-enrichment", "5"]),
]


# --------------------------------------------------------------------------
# Paired mode.
#
# The single-end corpus measures the ordinary matcher. It cannot price the
# PE-overlap path, because that path needs two mates to pin a boundary. The
# paired corpus therefore builds mate pairs with a KNOWN insert I, and scores
# BOTH mates against the truth.
#
# Geometry (read length L, insert I, readthrough t = L - I):
#
#   R1 = frag[0..I] + adapter1[:t]              (adapter1 = TruSeq/R1)
#   R2 = revcomp(frag[0..I]) + adapter2[:t]     (adapter2 = TruSeq/R2)
#
# Then revcomp(R2) = revcomp(adapter2[:t]) + frag[0..I], whose suffix of
# length I is exactly frag[0..I] -- R1's prefix of length I. The detector
# pins the boundary at I, so t3_1 = L - I = t.
#
# R2's tail is TruSeq/R2 (a set member) rather than a reverse complement, so
# the ordinary matcher CAN trim it: the two paths are then directly
# comparable, and any difference is the PE path's doing.
# --------------------------------------------------------------------------

PAIRED_CONFIGS = [
    ("shipped (mo=14, eps=0.2)", []),
    ("mo=12", ["--min-adapter-overlap", "12"]),
    ("mo=10", ["--min-adapter-overlap", "10"]),
    ("mo=8", ["--min-adapter-overlap", "8"]),
]


def build_corpus_paired(rng, n):
    """-> (pairs, {id: (t, cls)}); pairs is [(id, r1_seq, r2_seq)]."""
    pairs, truth = [], {}
    for k in range(n):
        rid = f"P{k}"
        if rng.random() < 0.70:
            # No readthrough: the fragment is at least as long as the read.
            # The insert may exceed L (the mates do not overlap at all) or
            # equal it (they touch exactly), so the clean population covers
            # both "no overlap" and "full overlap, zero readthrough".
            i = L + rng.randint(0, 8)
            frag = dna(rng, i)
            r1 = frag[:L]
            r2 = revcomp(frag[i - L:i])
            t, cls = 0, "clean"
        else:
            t = rng.randint(1, 33)
            cls = "plain"
            frag = dna(rng, L - t)
            tail1 = list(TRUSEQ_R1[:t])
            tail2 = list(TRUSEQ_R2[:t])
            if t and cls == "plain":
                r = rng.random()
                if r < 0.10:                      # NovaSeq no-signal: last bases G
                    g = rng.randint(3, 8)
                    for j in range(max(0, len(tail1) - g), len(tail1)):
                        tail1[j] = "G"
                    for j in range(max(0, len(tail2) - g), len(tail2)):
                        tail2[j] = "G"
                    cls = "g"
                elif r < 0.15:                    # a couple of ambiguous calls
                    for j in rng.sample(range(len(tail1)), min(2, len(tail1))):
                        tail1[j] = "N"
                    for j in rng.sample(range(len(tail2)), min(2, len(tail2))):
                        tail2[j] = "N"
                    cls = "n"
                elif r < 0.18:                    # substitution errors in the adapter
                    for _ in range(rng.randint(3, 4)):
                        j = rng.randrange(len(tail1))
                        tail1[j] = rng.choice([b for b in "ACGT" if b != tail1[j]])
                    for _ in range(rng.randint(3, 4)):
                        j = rng.randrange(len(tail2))
                        tail2[j] = rng.choice([b for b in "ACGT" if b != tail2[j]])
                    cls = "sub"
            r1 = frag + "".join(tail1)
            r2 = revcomp(frag) + "".join(tail2)
        pairs.append((rid, r1, r2))
        truth[rid] = (t, cls)
    return pairs, truth


def run_paired(fp, d, extra, tag, pe_mode, evidence=None):
    """Run fastprep on the pair; return the PXMT dict (both mates)."""
    o1 = os.path.join(d, f"o1_{tag}.fq")
    o2 = os.path.join(d, f"o2_{tag}.fq")
    px = os.path.join(d, f"o_{tag}.pxmt")
    for p in (o1, o2, px):
        if os.path.exists(p):
            os.remove(p)
    cmd = [fp, "-1", os.path.join(d, "r1.fq"), "-2", os.path.join(d, "r2.fq"),
           "-o", o1, "-p", o2, "--metadata", px,
           "--adapter-preset", "illumina", "--pe-overlap", pe_mode,
           "--no-dedup", "--no-mean-q", "--min-len", "0", "-j", "1", "-q", *extra]
    if evidence:
        cmd += ["--adapter-evidence", evidence]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"fastprep failed ({tag}): {r.stderr[-1200:]}")
    return parse_pxmt(px)


def score_paired(px, truth):
    """Bucket both mates against the truth.

    Keys: "ALL", "R1", "R2" and the per-class buckets, each holding
    exact / under / over / dropped counts.
    """
    out = {}

    def bump(scope, cls, bucket):
        out.setdefault(scope, {}).setdefault(cls, {}).setdefault(bucket, 0)
        out[scope][cls][bucket] += 1

    for rid, (t, cls) in truth.items():
        for mate, scope in ((0, "R1"), (1, "R2")):
            st, t5a, t3a, t5p, t3p, _ = px[(rid, mate)]
            cut = t5a + t3a + t5p + t3p
            if st != 0:
                b = "dropped"
            elif cut == t:
                b = "exact"
            elif cut < t:
                b = "under"
            else:
                b = "over"
            bump(scope, cls, b)
            bump(scope, "ALL", b)
            bump("ALL", cls, b)
            bump("ALL", "ALL", b)
    return out


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
    ap.add_argument("--paired", action="store_true",
                    help="build mate pairs and score both mates (PE path)")
    ap.add_argument("--pe-mode", default="always", choices=["never", "auto", "always"],
                    help="--pe-overlap mode for the paired corpus (default always)")
    ap.add_argument("--evidence", default=None, choices=["floor", "boundary"],
                    help="--adapter-evidence policy (default: the tool's own default)")
    a = ap.parse_args()
    if a.paired:
        return main_paired(a)
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


def write_pairs(d, pairs):
    with open(os.path.join(d, "r1.fq"), "w") as f1, \
         open(os.path.join(d, "r2.fq"), "w") as f2:
        for rid, s1, s2 in pairs:
            f1.write(f"@{rid}\n{s1}\n+\n{'I' * len(s1)}\n")
            f2.write(f"@{rid}\n{s2}\n+\n{'I' * len(s2)}\n")


def main_paired(a):
    """The paired two-axis A/B: both mates scored against a known insert."""
    fp = os.path.abspath(a.fastprep)
    rng = random.Random(a.seed)
    d = tempfile.mkdtemp(prefix="utabp_")
    pairs, truth = build_corpus_paired(rng, a.reads)
    write_pairs(d, pairs)
    n_clean = sum(1 for t, c in truth.values() if c == "clean")
    carriers = a.reads - n_clean
    n_mates = 2 * a.reads
    print(f"paired corpus: {a.reads} pairs, {n_clean} without readthrough, "
          f"{carriers} with readthrough ({100*carriers/a.reads:.1f}%), "
          f"--pe-overlap {a.pe_mode}, "
          f"--adapter-evidence {a.evidence or 'default'}\n")
    print(f"{'config':26s} {'exact':>7} {'under':>7} {'over':>6} | "
          f"{'und/carrier':>11} {'false/cmp':>10} | {'R1 und':>7} {'R2 und':>7}")

    report: dict = {"mode": "paired", "pe_mode": a.pe_mode,
                    "evidence": a.evidence, "reads": a.reads,
                    "clean_pairs": n_clean, "carriers": carriers,
                    "n_mates": n_mates, "configs": {}}
    for name, extra in PAIRED_CONFIGS:
        tag = name.split()[0].replace("=", "").replace(",", "")
        px = run_paired(fp, d, extra, tag, a.pe_mode, a.evidence)
        s = score_paired(px, truth)
        allm = s["ALL"]["ALL"]        # scope "ALL", class "ALL"
        und = allm.get("under", 0)
        # The clean axis: a mate with NO readthrough that got cut anyway.
        fp_clean = (s["ALL"]["clean"].get("over", 0)
                    + s["ALL"]["clean"].get("under", 0))
        n_clean_mates = 2 * n_clean
        undertrim_rate = und / (2 * carriers)
        clean_false_rate = fp_clean / n_clean_mates
        row = dict(name=name, extra=extra, buckets=s,
                   undertrim_rate=undertrim_rate,
                   clean_false_trim_rate=clean_false_rate)
        report["configs"][name] = row
        print(f"{name:26s} {allm.get('exact',0):7d} {und:7d} "
              f"{allm.get('over',0):6d} | {100*undertrim_rate:10.2f}% "
              f"{100*clean_false_rate:9.3f}% | "
              f"{s['R1']['ALL'].get('under',0):7d} {s['R2']['ALL'].get('under',0):7d}")
    with open(a.out, "w") as f:
        json.dump(report, f, indent=1)
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
