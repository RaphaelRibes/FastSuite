#!/usr/bin/env python3
"""prefix_features.py — can any per-read feature beat the blind 5' G-run rule?

The campaign priced two policies and found nothing in between: trimming every
leading 1-4 bp G-run (`--head-dropout 4`, 7:1 score +0.6749) or trimming none
(`--head-dropout 0`, +0.1680 at the anchor but +1.3180 at 1:1). Every quality
and margin flag the tool already ships is inert on this corpus, so the only way
to do better is a signal that separates an injected dark-cycle run from real
biology *per read*.

This script measures that separability directly on the ambiguous set: mates whose
read starts with 1-4 G bases. For each such mate it records what a trimmer can
see (run length, the first base after the run, GC of the next 10 bases, mean Q of
the run and of the next 10 bases, read length, the mate's own 5' run) and the
ground-truth label (was a 5' run injected on this mate).

Gain of a rule is in mate-equivalents against always trimming the run:

    gain = kept_biology - w_under * kept_artefact

Keeping a clean mate saves one overtrim unit; keeping an artefact costs w_under
units instead of the one unit trimming would have cost. A rule is worth shipping
only if it is positive at the anchor, not at some r chosen after the fact.

    prefix_features.py --stride 5 --out /tmp/prefix.tsv
"""
import argparse
import csv
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import overtrim_vs_fastp as ovf  # noqa: E402

FIELDS = ["dataset", "mate", "obs_len", "inj_len", "next_base", "gc10",
          "q_run", "q_next", "mate_run", "read_len"]


def leading_run(seq):
    n = 0
    for b in seq:
        if b in "gG":
            n += 1
        else:
            break
    return n


def read_fastq(path):
    """(read id, seq, qual) for one file, upper-cased sequence."""
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip("\n").upper()
            fh.readline()
            qual = fh.readline().rstrip("\n")
            yield h[1:].split()[0], seq, qual


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=5,
                    help="take every Nth winner: a separability probe does not "
                         "need the whole corpus, but the rule it suggests does")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit(f"{a.out} exists; move it aside")

    winners = list(csv.DictReader(open(f"{ovf.R}/winners.tsv"), delimiter="\t"))[::a.stride]
    print(f"datasets: {len(winners)}", flush=True)
    total = 0
    done = 0
    with open(a.out, "w") as out:
        out.write("\t".join(FIELDS) + "\n")
        for w in winners:
            dom, acc = w["domain"], w["run_accession"]
            v0 = f"{ovf.W}/{dom}/{acc}/v0"
            gtp = f"{v0}/sim.v0.gt.tsv"
            p1, p2 = f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz"
            if not (os.path.exists(gtp) and os.path.exists(p1) and os.path.exists(p2)):
                continue
            gt = ovf.gt_rows(gtp)
            coremap = ovf.core_to_pid([p1, p2])

            # pass 1: the 5' run of every mate, so pairs can be compared
            runs = {}
            for mate, path in ((1, p1), (2, p2)):
                for rid, seq, _ in read_fastq(path):
                    runs[(rid, mate)] = leading_run(seq)

            # pass 2: features for the ambiguous mates only
            n = 0
            for mate, path in ((1, p1), (2, p2)):
                for rid, seq, qual in read_fastq(path):
                    obs = runs[(rid, mate)]
                    if not 1 <= obs <= 4:
                        continue
                    pid = coremap.get((rid, mate))
                    if pid is None or pid not in gt:
                        continue
                    g = gt[pid]
                    inj = g[0] if mate == 1 else g[2]
                    nxt = seq[obs] if obs < len(seq) else "N"
                    win = seq[obs:obs + 10]
                    gc10 = (win.count("G") + win.count("C")) / len(win) if win else 0.0
                    qs = [ord(c) - 33 for c in qual]
                    q_run = sum(qs[:obs]) / obs if obs else 0.0
                    qn = qs[obs:obs + 10]
                    q_next = sum(qn) / len(qn) if qn else 0.0
                    other = runs.get((rid, 3 - mate), 0)
                    out.write(f"{dom}/{acc}\t{mate}\t{obs}\t{inj}\t{nxt}\t{gc10:.3f}\t"
                              f"{q_run:.1f}\t{q_next:.1f}\t{other}\t{len(seq)}\n")
                    n += 1
            total += n
            done += 1
            print(f"  [{done}/{len(winners)}] {dom}/{acc}: {n} ambiguous mates", flush=True)

    if total == 0:
        sys.exit("no ambiguous mates found: check the corpus path and winners.tsv")
    print(f"wrote {a.out}: {total} ambiguous mates over {done} datasets")


if __name__ == "__main__":
    main()
