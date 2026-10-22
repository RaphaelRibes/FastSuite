#!/usr/bin/env python3
"""prefix_policy.py — is the 5' G-run question a per-LIBRARY question?

prefix_features.py showed the decision is undecidable per read. But the injected
fraction is not the same in every library: on the 21-library probe it runs from
0% to 63%, and at the 7:1 anchor a library is better off keeping its 1-4 bp
leading runs unless at least one in eight of them is injected. A single global
--head-dropout applies one answer to all of them.

This script asks whether the answer is estimable from the reads alone, which is
what a tool would need. Under biology the rate of reads beginning with 1-4 G
bases should follow the library's own G frequency (p - p^5 for a run of 1-4);
any excess over that is the injected share. So:

    est_injected = observed_rate_1_4 - (p_G - p_G^5)

with p_G measured over the same reads. The rule a tool could implement is then
"trim this library's short 5' runs iff est_injected > w_over / (w_over + w_under)
= 1/8", i.e. the cost model's own break-even, applied per library instead of per
read.

    prefix_policy.py --out /tmp/policy.tsv
"""
import argparse
import collections
import csv
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import overtrim_vs_fastp as ovf  # noqa: E402

W_OVER, W_UNDER = 1.0, 7.0
THRESH = W_OVER / (W_OVER + W_UNDER)


def leading_run(seq):
    n = 0
    for b in seq:
        if b in "gG":
            n += 1
        else:
            break
    return n


def scan(path, acc):
    mates = 0
    run_hist = collections.Counter()
    g = total = 0
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip("\n").upper()
            fh.readline()
            fh.readline()
            mates += 1
            run_hist[min(leading_run(seq), 8)] += 1
            g += seq.count("G")
            total += len(seq)
    return mates, run_hist, (g / total if total else 0.25)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stride", type=int, default=5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit(f"{a.out} exists; move it aside")

    winners = list(csv.DictReader(open(f"{ovf.R}/winners.tsv"), delimiter="\t"))[::a.stride]
    rows = []
    for w in winners:
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{ovf.W}/{dom}/{acc}/v0"
        p1 = f"{v0}/sim.v0_R1.fastq.gz"
        if not os.path.exists(p1):
            continue
        mates, hist, pg = scan(p1, acc)
        amb = sum(hist[k] for k in (1, 2, 3, 4))
        obs = amb / mates if mates else 0.0
        under_bio = pg - pg ** 5          # expected 1-4 bp leading rate under biology
        est = max(0.0, min(1.0, obs - under_bio))
        row = {"dataset": f"{dom}/{acc}", "mates": mates, "amb": amb,
               "obs_rate": obs, "p_g": pg, "exp_bio": under_bio,
               "est_injected": est}
        # The whole leading-run histogram, so a self-calibrated estimator (p
        # fitted from the tail the injection barely touches) can be compared to
        # the base-composition one without another pass over the corpus.
        for k in range(0, 9):
            row[f"run_{k}"] = hist.get(k, 0)
        rows.append(row)
        print(f"  {dom}/{acc}: obs {100*obs:5.2f}%  pG {pg:.3f}  bio {100*under_bio:5.2f}%"
              f"  est {100*est:5.2f}%", flush=True)

    if not rows:
        sys.exit("no datasets scanned")
    with open(a.out, "w") as fh:
        fh.write("\t".join(rows[0].keys()) + "\n")
        for r in rows:
            fh.write("\t".join(f"{v:.6f}" if isinstance(v, float) else str(v)
                               for v in r.values()) + "\n")
    print(f"wrote {a.out} ({len(rows)} libraries); break-even injected share is {THRESH:.3f}")


if __name__ == "__main__":
    main()
