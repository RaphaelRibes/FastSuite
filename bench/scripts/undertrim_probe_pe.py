#!/usr/bin/env python3
"""Paired-end half of the undertrim probe: which refusal cases cost a trim.

Constructs mate pairs with a known insert I < L, so each mate carries
L - I bases of its own adapter readthrough, then runs fastprep with
--pe-overlap auto (the shipped default whenever an explicit adapter set is
configured) against --pe-overlap never. The difference is what the geometry
bought.

Usage: python3 bench/scripts/undertrim_probe_pe.py [--n 1000]
"""
import argparse
import json
import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from undertrim_probe import TRUSEQ_R1, TRUSEQ_R2, dna, parse_pxmt, run_fastprep  # noqa: E402

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def write_pair(d1, d2, pairs):
    with open(os.path.join(d1, "r1.fq"), "w") as f1, open(os.path.join(d1, "r2.fq"), "w") as f2:
        for rid, s1, s2 in pairs:
            f1.write(f"@{rid}\n{s1}\n+\n{'I'*len(s1)}\n")
            f2.write(f"@{rid}\n{s2}\n+\n{'I'*len(s2)}\n")


def run_paired(fp, d, mode, extra=()):
    o1 = os.path.join(d, f"o1_{mode}.fq")
    o2 = os.path.join(d, f"o2_{mode}.fq")
    p1 = os.path.join(d, f"r1_{mode}.pxmt")
    for p in (o1, o2, p1):
        if os.path.exists(p):
            os.remove(p)
    cmd = [fp, "-1", os.path.join(d, "r1.fq"), "-2", os.path.join(d, "r2.fq"),
           "-o", o1, "-p", o2, "--metadata", p1,
           "--adapter-preset", "illumina", "--pe-overlap", mode,
           "--no-dedup", "--no-mean-q", "--min-len", "0", "-j", "1", "-q", *extra]
    r = __import__("subprocess").run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"fastprep failed ({mode}): {r.stderr[-1500:]}")
    px = parse_pxmt(p1)
    return ({k[0]: v for k, v in px.items() if k[1] == 0},
            {k[0]: v for k, v in px.items() if k[1] == 1})


def case(name, fp, rng, n, L, insert, g_contam=0, n_contam=0, adapter1=TRUSEQ_R1, adapter2=TRUSEQ_R2):
    pairs, truth = [], {}
    for k in range(n):
        frag = dna(rng, insert)
        t = L - insert
        a1 = list(adapter1[:t] if t <= len(adapter1) else adapter1 + "A" * (t - len(adapter1)))
        a2 = list(adapter2[:t] if t <= len(adapter2) else adapter2 + "T" * (t - len(adapter2)))
        for _ in range(n_contam):
            if a1:
                a1[rng.randrange(len(a1))] = "N"
        s1 = list(frag + "".join(a1))
        s2 = list(revcomp(frag) + "".join(a2))
        if g_contam:
            for j in range(len(s1) - g_contam, len(s1)):
                s1[j] = "G"
            for j in range(len(s2) - g_contam, len(s2)):
                s2[j] = "G"
        rid = f"{name}_{k}"
        pairs.append((rid, "".join(s1), "".join(s2)))
        truth[rid] = t
    d = tempfile.mkdtemp(prefix="utpe_")
    write_pair(d, d, pairs)
    out = {}
    for mode in ("auto", "never", "always"):
        p1, p2 = run_paired(fp, d, mode)
        for mate, px in (("r1", p1), ("r2", p2)):
            ex = none = part = 0
            for rid, t in truth.items():
                st, t5a, t3a, t5p, t3p, _ = px[rid]
                cut = t3a + t3p
                if st != 0:
                    continue
                if cut == t:
                    ex += 1
                elif cut == 0:
                    none += 1
                else:
                    part += 1
            out[(mode, mate)] = dict(exact=ex, none=none, partial=part)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--out", default="bench/results/undertrim-probe-pe.json")
    a = ap.parse_args()
    fp = os.path.abspath(a.fastprep)
    rng = random.Random(11)
    L = 150
    cases = [
        ("tail30_clean", dict(insert=120)),
        ("tail20_clean", dict(insert=130)),
        ("tail10_short", dict(insert=140)),
        ("tail6_short", dict(insert=144)),
        ("tail3_short", dict(insert=147)),
        ("tail20_g6", dict(insert=130, g_contam=6)),
        ("tail20_g4", dict(insert=130, g_contam=4)),
        ("tail20_n2", dict(insert=130, n_contam=2)),
    ]
    report = {"cases": {}}
    print(f"{'case':14s} {'mode':6s} {'mate':4s} {'exact':>6} {'none':>6} {'partial':>7}")
    for name, kw in cases:
        res = case(name, fp, rng, a.n, L, **kw)
        report["cases"][name] = {f"{m}/{mate}": v for (m, mate), v in res.items()}
        for (mode, mate), v in sorted(res.items()):
            print(f"{name:14s} {mode:6s} {mate:4s} {v['exact']:6d} {v['none']:6d} {v['partial']:7d}")
        print()
    with open(a.out, "w") as f:
        json.dump(report, f, indent=1)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
