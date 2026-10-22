#!/usr/bin/env python3
"""P(the G-run is an injected artefact | the observed run has length r).

This is the number that decides whether undertrim is recoverable at all.

A trimmer sees only the READ. The simulator's GT additionally records `r1_h` /
`r1_t` — the 5'/3' G-run length it INJECTED. So:

  observed run of length r, and the read's first r bases are all G
  -> the run is artefact if the injected length matches r
  -> it is biology if the injected length is 0 (the biology itself starts with G)
  -> an injected run of length k < r means both are present

On REAL data the trimmer cannot see `r1_h`. So the posterior below is the best
any trimmer could do from the sequence alone, and the Bayes-optimal rule is to
trim when P > 0.5. That is the ceiling for this bucket.

Usage: trim_posterior.py [--limit N] [--end 5|3]
"""
import argparse, collections, gzip, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import overtrim_vs_fastp as H


def leading_g(seq):
    n = 0
    for b in seq:
        if b in "gG": n += 1
        else: break
    return n


def trailing_g(seq):
    n = 0
    for b in reversed(seq):
        if b in "gG": n += 1
        else: break
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--end", choices=["5", "3"], default="5")
    ap.add_argument("--out", default=None,
                    help="write the per-run table as TSV, so the paper's table "
                         "can be generated instead of typed")
    a = ap.parse_args()

    winners = list(H.csv.DictReader(open(f"{H.R}/winners.tsv"), delimiter="\t"))
    if a.limit: winners = winners[:a.limit]

    # observed run length -> [n_artefact, n_biology]
    tab = collections.defaultdict(lambda: [0, 0])

    for i, w in enumerate(winners, 1):
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{H.W}/{dom}/{acc}/v0"
        gtp = f"{v0}/sim.v0.gt.tsv"
        s1, s2 = f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz"
        if not (os.path.exists(gtp) and os.path.exists(s1)): continue
        gt = {}
        with open(gtp) as fh:
            for r in H.csv.DictReader(fh, delimiter="\t"):
                gt[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                         int(r["r2_h"]), int(r["r2_t"]))
        # read the sim's own headers to recover pair_id and mate
        for path, mate in ((s1, 1), (s2, 2)):
            hdr = pid = None
            with gzip.open(path, "rt") as fh:
                for j, line in enumerate(fh):
                    k = j % 4
                    if k == 0:
                        hdr = line[1:].strip()
                        parts = hdr.split()[0].split(":")
                        x, y = int(parts[5]), int(parts[6])
                        pid = x + 10000 * y
                    elif k == 1:
                        g = gt.get(pid)
                        if g is None: continue
                        inj = g[0] if mate == 1 else g[2]
                        if a.end == "3":
                            inj = g[1] if mate == 1 else g[3]
                            obs = trailing_g(line.strip())
                        else:
                            obs = leading_g(line.strip())
                        if obs == 0: continue
                        tab[obs][0 if inj > 0 else 1] += 1

    print(f"\n  {a.end}' observed G-run length vs whether an artefact was injected")
    print(f"  {'r':>4} {'artefact':>10} {'biology':>10} {'total':>11} {'P(artefact|r)':>14}  Bayes action")
    tot_a = tot_b = 0
    for r in sorted(tab):
        if r > 12: continue
        na, nb = tab[r]
        t = na + nb; tot_a += na; tot_b += nb
        p = na / t if t else 0
        act = "TRIM" if p > 0.5 else "keep"
        bar = "#" * int(p * 30)
        print(f"  {r:>4} {na:>10} {nb:>10} {t:>11} {p:>13.3f}   {act:4s} {bar}")
    t = tot_a + tot_b
    print(f"\n  overall P(artefact) = {tot_a/t:.3f}   ({tot_a} artefact / {tot_b} biology of {t})")
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("r\tartefact\tbiology\ttotal\tp_artefact\n")
            for r in sorted(tab):
                if r > 12:
                    continue
                na, nb = tab[r]
                tt = na + nb
                fh.write(f"{r}\t{na}\t{nb}\t{tt}\t{na/tt:.4f}\n" if tt else "")
        print(f"  wrote {a.out}")
    reach = sum(tab[r][0] for r in tab if tab[r][0] + tab[r][1] and
                tab[r][0] / (tab[r][0] + tab[r][1]) > 0.5)
    print(f"  artefacts a Bayes-optimal rule of r alone could recover: {reach} "
          f"({reach/tot_a*100:.1f}% of injected runs seen)")


if __name__ == "__main__":
    main()
