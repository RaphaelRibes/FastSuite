#!/usr/bin/env python3
"""Cell B -- synthetic PE-overlap cohort: build, then score with/without.

Run this ON THE DESKTOP (the corpus lives at /mnt/Fast/fastsuite-scratch):

    python3 bench/work/cellB/build_and_run_cellB.py

What it does
------------
1. Regenerates the cohort deterministically:

       fastsim --pairs 400000 --profile <repo>/tools/fastsim/profiles/cellB-pe-overlap.json \
                --out <cohort>/cohort/work/synthetic/cellB-pe-overlap/v0/sim.v0 \
                --seed 4242 --force --no-gzip --drop-rate 0 --dup-rate 0

   `--seed 4242` is fixed, so the GT/FASTQ pair is byte-identical every time.
   `--no-gzip` writes plain FASTQ: `trim_report.py` reads the corpus FASTQs with
   `gzip.open(..., "rt")`, which on this system rejects a plain file, so the
   cohort is written plain to keep the reader uniform.
   `--drop-rate 0 --dup-rate 0` remove the two classes that have nothing to do
   with the PE question (a forced-too-short fragment and a verbatim replay) --
   the profile already carries `dup_rate: 0.0`.

2. Places the cohort into the layout `bench/scripts/trim_report.py` expects --
   it reads `<work>/<domain>/<acc>/v0/{sim.v0.gt.tsv,sim.v0_R1.fastq.gz,
   sim.v0_R2.fastq.gz}` and a `<results>/winners.tsv`:

       <cohort>/cohort/work/<domain>/<acc>/v0/sim.v0_*
       <cohort>/cohort/results/winners.tsv   one row: the published fastprep
                                             winner config for that dataset (so
                                             the baselines reproduce the real
                                             harness)

3. Runs the trim report TWICE, with and without the mode:

       FASTPREP_VARIANTS='[["fp-pe",["--pe-overlap"]]]' \\
       python3 trim_report.py --mode decomp --work <cohort>/work \\
           --results <cohort>/results --out <out>/cellB-with-pe.tsv

       FASTPREP_VARIANTS='' \\
       python3 trim_report.py --mode decomp --work <cohort>/work \\
           --results <cohort>/results --out <out>/cellB-without-pe.tsv

   The report REFUSES to overwrite an existing TSV, so each run gets its own
   --out path. The second run's `fastprep-hd0` row is the with/without
   comparison's baseline; `fp-pe` is the mode's row in the first run.

4. Prints the headline numbers: per tool and tolerance, the exact/overtrim/
   undertrim_recoverable/undertrim_irreducible/dropped decomposition, the
   per-mate split, and the readthrough / in_window slices.

This script writes nothing into the repo tree except the cohort itself (under
--cohort-root) and the report TSVs (under --out-root).
"""
from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

# --- the cohort's shape ----------------------------------------------------
PAIRS = 400_000
SEED = 4242
# The dataset key inside the corpus layout. A synthetic cohort is one dataset;
# it is named after what it is, not after an accession, so the domain column
# says `synthetic` (the trim report only uses `domain` for the per-domain
# retention table, which this run does not write).
DOMAIN = "synthetic"
ACCESSION = "cellB-pe-overlap"
# The published fastprep winner config for this cohort. It is a single dataset,
# so one row is enough; the values are the suite's shipped default winner
# (min_run 12, mercy 3, phred_sat 20, gate_shape 1.0, min_len 0.9L, mean_q 25,
# adapter_err 0.2, min_overlap 14) -- the same config the corpus sweep picked on
# 100 of 101 datasets, and the one the report's published baselines reproduce.
WINNER = {
    "domain": DOMAIN,
    "run_accession": ACCESSION,
    "config_id": "cellB",
    "min_run": "12",
    "mercy": "3",
    "phred_sat": "20",
    "gate_shape": "1.0",
    "min_len": "0.9",
    "mean_q": "25.0",
    "adapter_err": "0.2",
    "min_overlap": "14",
    "r2_reverse": "0",
    "interior": "0",
    "min_len_effective": "0.9L",
}
WINNER_COLS = list(WINNER)

# --- what to run -----------------------------------------------------------
VARIANTS_JSON = '[["fp-pe",["--pe-overlap"]]]'


def sh(cmd: list[str], **kw) -> subprocess.CompletedProcess:
    print("+", " ".join(cmd), flush=True)
    return subprocess.run(cmd, **kw)


def build_cohort(fastsim: str, profile: str, v0: str) -> None:
    os.makedirs(v0, exist_ok=True)
    prefix = os.path.join(v0, "sim.v0")
    sh([fastsim, "--pairs", str(PAIRS), "--profile", profile,
        "--out", prefix, "--seed", str(SEED), "--force",
        "--drop-rate", "0", "--dup-rate", "0"])
    for suffix in ("_R1.fastq.gz", "_R2.fastq.gz", ".gt.tsv"):
        p = prefix + suffix
        if not os.path.exists(p):
            raise SystemExit(f"fastsim did not produce {p}")
        print(f"  {p}  {os.path.getsize(p)} bytes")


def write_winners(results: str) -> None:
    os.makedirs(results, exist_ok=True)
    path = os.path.join(results, "winners.tsv")
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=WINNER_COLS, delimiter="\t",
                            lineterminator="\n")
        wr.writeheader()
        wr.writerow(WINNER)
    print(f"wrote {path}")


def run_report(tag: str, trim_report: str, work: str, results: str, out_root: str,
               scratch: str, variants: str) -> str:
    """Run trim_report.py once; return the TSV path it wrote."""
    os.makedirs(out_root, exist_ok=True)
    os.makedirs(scratch, exist_ok=True)
    out = os.path.join(out_root, f"cellB-{tag}.tsv")
    env = dict(os.environ)
    env["FASTSCRATCH"] = scratch
    if variants:
        env["FASTPREP_VARIANTS"] = variants
    else:
        env.pop("FASTPREP_VARIANTS", None)
    r = sh([sys.executable, trim_report, "--mode", "decomp",
            "--work", work, "--results", results,
            "--out", out, "--tol", "0"], env=env)
    if r.returncode != 0:
        raise SystemExit(f"trim_report ({tag}) failed rc={r.returncode}")
    return out


def show(path: str) -> None:
    """Print the main table and the two side tables of one report."""
    print(f"\n===== {path} =====")
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            print("\t".join(row))
    for suffix, label in (("-permate.tsv", "per mate"),
                          ("-pe.tsv", "PE slice")):
        side = path.rsplit(".", 1)[0] + suffix
        if os.path.exists(side):
            print(f"--- {label}: {side} ---")
            with open(side) as fh:
                for row in csv.reader(fh, delimiter="\t"):
                    print("\t".join(row))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fastsim", default=os.environ.get(
        "FASTSIM", os.path.expanduser("~/FastSuite/tools/fastsim/target/release/fastsim")))
    ap.add_argument("--profile", default=os.path.join(
        REPO, "tools", "fastsim", "profiles", "cellB-pe-overlap.json"))
    ap.add_argument("--trim-report", default=os.path.join(
        REPO, "bench", "scripts", "trim_report.py"))
    ap.add_argument("--cohort-root", default=os.environ.get(
        "CELLB_ROOT",
        "/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap/cellB"))
    ap.add_argument("--out-root", default=os.environ.get(
        "CELLB_OUT", "/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap/results"))
    ap.add_argument("--scratch-root", default=os.environ.get(
        "CELLB_SCRATCH", "/mnt/Fast/fastsuite-scratch/bench-work/20261019-pe-overlap/trim-cellB"))
    ap.add_argument("--skip-sim", action="store_true",
                    help="reuse an existing cohort instead of re-simulating")
    ap.add_argument("--show-only", action="store_true",
                    help="print the already-written TSVs and exit")
    a = ap.parse_args()

    cohort = os.path.join(a.cohort_root, "cohort")
    work = os.path.join(cohort, "work")
    results = os.path.join(cohort, "results")
    v0 = os.path.join(work, DOMAIN, ACCESSION, "v0")

    if a.show_only:
        for tag in ("with-pe", "without-pe"):
            show(os.path.join(a.out_root, f"cellB-{tag}.tsv"))
        return

    for p in (a.fastsim, a.profile, a.trim_report):
        if not os.path.exists(p):
            raise SystemExit(f"missing {p}")

    if not a.skip_sim:
        build_cohort(a.fastsim, a.profile, v0)
    write_winners(results)

    # The report's `build_posterior` walks the SAME work root, so the posterior
    # it measures is this cohort's own (a synthetic cohort's posterior is not the
    # corpus's -- that is the point of a dedicated cohort).
    with_pe = run_report("with-pe", a.trim_report, work, results, a.out_root,
                         os.path.join(a.scratch_root, "with-pe"), VARIANTS_JSON)
    without_pe = run_report("without-pe", a.trim_report, work, results, a.out_root,
                            os.path.join(a.scratch_root, "without-pe"), "")

    show(with_pe)
    show(without_pe)
    print(f"\nwrote {with_pe}")
    print(f"wrote {without_pe}")


if __name__ == "__main__":
    main()
