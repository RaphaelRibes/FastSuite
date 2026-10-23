#!/usr/bin/env python3
"""Measure how the shipped --pe-correct gate acts on a flat high-Q corpus.

The posterior rule abstains when the two mates' qualities are equal (the odds
are exactly even), so on a library whose overlap disagreements all carry high,
similar quality it should rewrite little or nothing. This probe quantifies
that: simulate a corpus from the real-profile mandrill fastsim profile (whose
per-base quality is uniformly high), run the shipped fastprep binary with
--pe-correct at several gates, and record bases_corrected.

The corpus is regenerated deterministically (fastsim, fixed seed), so every
number in the committed TSV is re-derivable from this script alone:

    /usr/bin/python3 bench/scripts/pe_correct_gate_probe.py \
        --fastsim tools/fastsim/target/release/fastsim \
        --fastprep tools/fastprep/target/release/fastprep \
        --profile tools/fastsim/profiles/mandrill_E358.json \
        --pairs 50000 --seed 11 \
        --out bench/results/pe-correct-gate-probe.tsv

The TSV is the committed evidence; the paper quotes the TSV, never this run.
"""
import argparse
import csv
import json
import os
import subprocess
import tempfile


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    assert r.returncode == 0, f"failed: {' '.join(cmd)}\n{r.stderr[-500:]}"
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastsim", required=True)
    ap.add_argument("--fastprep", required=True)
    ap.add_argument("--profile", required=True)
    ap.add_argument("--pairs", type=int, default=50000)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--gates", default="0.0,0.5,2.0")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    gates = [float(g) for g in args.gates.split(",")]

    version = run([args.fastprep, "--version"]).stdout.strip().splitlines()[-1]
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "corpus")
        run([args.fastsim, "--pairs", str(args.pairs), "--profile", args.profile,
             "--seed", str(args.seed), "-o", prefix])
        r1, r2 = f"{prefix}_R1.fastq.gz", f"{prefix}_R2.fastq.gz"
        for f in (r1, r2):
            assert os.path.getsize(f) > 0, f"corpus half missing: {f}"

        rows = []
        for g in gates:
            o1, o2 = os.path.join(tmp, f"g{g}_R1.fq"), os.path.join(tmp, f"g{g}_R2.fq")
            rep = os.path.join(tmp, f"g{g}.json")
            run([args.fastprep, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
                 "-j", "1", "--force", "-q", "--no-dedup",
                 f"--pe-correct={g}", "--json", rep])
            m = json.load(open(rep))["metrics"]
            rows.append({
                "gate": g,
                "pairs": args.pairs,
                "seed": args.seed,
                "profile": os.path.basename(args.profile),
                "fastprep": version,
                "bases_corrected": int(m.get("bases_corrected", 0)),
            })

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out}")
    for r in rows:
        print(f"  gate {r['gate']}: {r['bases_corrected']} bases corrected")


if __name__ == "__main__":
    main()
