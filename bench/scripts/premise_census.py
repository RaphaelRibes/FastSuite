#!/usr/bin/env python3
"""Premise test for the two opt-in flags, over the WHOLE corpus, reading the
raw reads and the ground truth only (no trimmer involved).

For every mate, find the terminal homopolymer at each physical end and ask the
ground truth whether a correct trimmer should move that end at all
(`gt_trim > 0`: yes, there is an injected artefact there; `gt_trim == 0`: no).
That is exactly the question a run-length-only rule has to answer, so the table
is the ceiling for any threshold at that end -- per BASE and per run length.

Also censuses terminal-N runs: how many mates carry one, and whether the GT
declares the boundary moved (`gt_trim > 0`) or not (`gt_trim == 0`, i.e.
trimming those bases is an OVERTRIM).
"""
import collections, csv, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import overtrim_vs_fastp as H
import trim_report as T

R, W = H.R, H.W
CAP = 15                       # run lengths > CAP are pooled as f">{CAP}"

tab = collections.defaultdict(lambda: collections.Counter())
n_mates = 0
alln = 0
allsame = 0


def end_run(seq, end):
    b = seq[0] if end == "5" else seq[-1]
    n = 0
    it = seq if end == "5" else reversed(seq)
    for c in it:
        if c == b:
            n += 1
        else:
            break
    return b, n


winners = list(csv.DictReader(open(f"{R}/winners.tsv"), delimiter="\t"))
for i, w in enumerate(winners, 1):
    v0 = f"{W}/{w['domain']}/{w['run_accession']}/v0"
    gtp, s1, s2 = (f"{v0}/sim.v0.gt.tsv", f"{v0}/sim.v0_R1.fastq.gz",
                   f"{v0}/sim.v0_R2.fastq.gz")
    if not (os.path.exists(gtp) and os.path.exists(s1)):
        print(f"SKIP missing {w['domain']}/{w['run_accession']}", file=sys.stderr)
        continue
    gt = T.gt_full(gtp)
    for path, mate in ((s1, 1), (s2, 2)):
        for pid, m, seq in H.read_fastq_keyed(path):
            g = gt.get(pid)
            if g is None:
                continue
            h, t = (g[0], g[1]) if m == 1 else (g[2], g[3])
            n_mates += 1
            if seq and seq.count(seq[0]) == len(seq):
                allsame += 1
                if seq[0] == "N":
                    alln += 1
            for end, gt_trim in (("5", h), ("3", t)):
                b, r = end_run(seq, end)
                if r == 0:
                    continue
                key = (end, b, r if r <= CAP else CAP + 1)
                tab[key]["n"] += 1
                if gt_trim > 0:
                    tab[key]["gt_moved"] += 1
    if i % 20 == 0:
        print(f"  [{i}/{len(winners)}]", file=sys.stderr, flush=True)

print(f"mates scanned: {n_mates}   all-same-base reads: {allsame}   "
      f"all-N reads: {alln}")
print("\n=== 3' terminal homopolymer: P(GT moves that end | base, run length) ===")
print(f"{'base':>5} {'run':>5} {'mates':>10} {'gt_moved':>10} {'P(moved)':>9}")
for r in list(range(3, CAP + 1)) + [CAP + 1]:
    for b in "ACGTN":
        s = tab[("3", b, r)]
        if not s["n"]:
            continue
        rr = f">{CAP}" if r == CAP + 1 else str(r)
        print(f"{b:>5} {rr:>5} {s['n']:10d} {s['gt_moved']:10d} "
              f"{s['gt_moved'] / s['n']:9.4f}")
print("\n=== 5' terminal homopolymer: P(GT moves that end | base, run length) ===")
print(f"{'base':>5} {'run':>5} {'mates':>10} {'gt_moved':>10} {'P(moved)':>9}")
for r in list(range(3, CAP + 1)) + [CAP + 1]:
    for b in "ACGTN":
        s = tab[("5", b, r)]
        if not s["n"]:
            continue
        rr = f">{CAP}" if r == CAP + 1 else str(r)
        print(f"{b:>5} {rr:>5} {s['n']:10d} {s['gt_moved']:10d} "
              f"{s['gt_moved'] / s['n']:9.4f}")

print("\n=== terminal-N summary (what --trim-terminal-n would cut) ===")
for end in ("5", "3"):
    for r in list(range(1, 8)) + [CAP, CAP + 1]:
        s = tab[(end, "N", r)]
        if not s["n"]:
            continue
        rr = f">{CAP}" if r == CAP + 1 else str(r)
        print(f"  {end}' N-run len {rr:>3}: mates {s['n']:9d}  "
              f"GT says the boundary MOVED (trim correct) {s['gt_moved']:8d} "
              f"({s['gt_moved'] / s['n'] * 100:5.2f}%)  "
              f"GT says KEEP (trim = overtrim) {s['n'] - s['gt_moved']:9d}")
