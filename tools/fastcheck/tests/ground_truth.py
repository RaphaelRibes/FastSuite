#!/usr/bin/env python3
"""Ground-truth recovery test for fastcheck (plan Task 13, PLAN-QC §7.3).

Generates paired reads from a synthetic profile with `fastsim`, runs `fastcheck`
on them, and asserts every recovered number against the simulator's *exact*
ground truth (`sim.gt.tsv`) within a stated tolerance. This is the difference
between "the numbers look plausible" and "the numbers are right": the ground
truth is the simulator's own per-pair record, not another tool's output.

No real data is involved (docs/DECISIONS.md D-11): the profile
(`tests/gt_profile.json`) is hand-written round numbers describing a mechanism,
so the expected recovery is computable by hand.

The ten checks and their tolerances are the plan's table
(`docs/plans/2026-10-12_fastcheck-exactness.md`, "Tolerances the ground-truth
test asserts"):

  #  quantity                          ground truth              tolerance
  1  per-position mean Q               30.0 (flat pos_q)         <= 0.05
  2  whole-read GC                     0.5 (profile gc)          +/- 0.03
  3  adapter rate (pairs)              gt.tsv adapter columns    +/- 1% of pairs
  4  3' G-run histogram, bin 5         gt.tsv r*_t == 5          +/- 4%
  5  3' G-run histogram, bin 10        gt.tsv r*_t == 10         +/- 4%
  6  5' head histogram, bin 1          convolution model         +/- 1.5%
  7  5' head histogram, bin 2          convolution model         +/- 1.5%
  8  5' head histogram, bin 3          convolution model         +/- 1.5%
  9  insert mean                       gt.tsv insert, in the     +/- 2%
                                       overlap-measurable set
  10 duplication level                 0 <= level <= dup_rate    lower bound

`mean_read_mean_q` is deliberately NOT asserted: it bins the floor of each
read's mean, biasing it low by ~0.5 (measured 29.53 for a true 30.0). The
per-position mean from the exact histogram is the unbiased quantity.

The 5' head model is the convolution the boundary rule implies: a read whose
injected head is k reports a leading G-run of exactly k, and a read with no
injected head (h=0) reports whatever biological run happens to start the read,
so observed(k) = P(h=k) + P(h=0) * P(random run = k), with P(run=k) =
pG^k * (1 - pG) and pG = gc/2.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# fastcheck's overlap anchor (src/insert.rs ANCHOR = 12). The documented
# measurement ceiling is `insert <= 2L - ANCHOR`.
ANCHOR = 12

GT_TSV_COLUMNS = [
    "pair_id", "insert", "r1_h", "r1_t", "r2_h", "r2_t",
    "r1_adapter", "r2_adapter", "full_flag", "is_dup", "r1_n", "r2_n",
]


def repo_root() -> Path:
    # tools/fastcheck/tests/ground_truth.py -> repo root is four parents up.
    return Path(__file__).resolve().parents[3]


def resolve_bin(explicit, env_key, default: Path) -> Path:
    cand = explicit or os.environ.get(env_key) or str(default)
    p = Path(cand)
    if not p.is_file():
        sys.exit(
            f"ground_truth: {p} does not exist.\n"
            f"  Build it first:  pixi run build   (or pass --{env_key.lower()})\n"
        )
    return p


def run(cmd, label):
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        sys.exit(f"ground_truth: {label} failed (exit {proc.returncode}): {' '.join(map(str, cmd))}")
    return proc


def histogram_to_pmf(h: dict) -> dict:
    total = sum(h.values())
    if total <= 0:
        sys.exit("ground_truth: profile carries an empty histogram")
    return {int(k): v / total for k, v in h.items()}


def load_gt(path: Path):
    """Return the per-pair rows of sim.gt.tsv as dicts of ints."""
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    if header != GT_TSV_COLUMNS:
        sys.exit(f"ground_truth: unexpected gt.tsv header: {header}")
    return [{c: int(v) for c, v in zip(header, row.split("\t"))} for row in lines[1:]]


def main() -> int:
    ap = argparse.ArgumentParser(description="fastcheck ground-truth recovery test")
    ap.add_argument("--pairs", type=int, default=20000, help="read pairs to simulate")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--insert-mean", type=float, default=200.0)
    ap.add_argument("--insert-sd", type=float, default=40.0)
    ap.add_argument("--dup-rate", type=float, default=0.02)
    ap.add_argument("--profile", type=Path, default=None,
                    help="synthetic profile (default: tests/gt_profile.json)")
    ap.add_argument("--workdir", type=Path, default=None,
                    help="where to write the generated reads (default: a temp dir)")
    ap.add_argument("--fastsim", default=None)
    ap.add_argument("--fastcheck", default=None)
    ap.add_argument("--keep", action="store_true", help="keep the work directory")
    args = ap.parse_args()

    root = repo_root()
    profile = args.profile or (Path(__file__).resolve().parent / "gt_profile.json")
    fastsim = resolve_bin(args.fastsim, "FASTSIM", root / "tools/fastsim/target/release/fastsim")
    fastcheck = resolve_bin(args.fastcheck, "FASTCHECK", root / "tools/fastcheck/target/release/fastcheck")

    prof = json.loads(profile.read_text())
    read_len = int(prof["read_len"])
    gc = float(prof["gc"])
    dup_rate = args.dup_rate
    adapter_seq = prof["adapter"]["seq"]

    work = args.workdir or Path(tempfile.mkdtemp(prefix="fastcheck-gt-"))
    work.mkdir(parents=True, exist_ok=True)
    prefix = work / "sim"

    # 1. Generate reads with exact ground truth.
    run(
        [
            str(fastsim), "--pairs", str(args.pairs), "--profile", str(profile),
            "--insert-mean", str(args.insert_mean), "--insert-sd", str(args.insert_sd),
            "--gc", str(gc), "--seed", str(args.seed), "--no-gzip",
            "--dup-rate", str(dup_rate), "-o", str(prefix),
        ],
        "fastsim",
    )
    # 2. QC them. --insert-hist needs the adapter to reach the readthrough regime.
    fc_json = work / "fc.json"
    pos_tsv = work / "pos.tsv"
    run(
        [
            str(fastcheck), str(prefix) + "_R1.fastq", str(prefix) + "_R2.fastq",
            "--json", str(fc_json), "--positions-tsv", str(pos_tsv),
            "--insert-hist", "--insert-pairs", str(args.pairs),
            "--adapter-seq", adapter_seq, "-q",
        ],
        "fastcheck",
    )

    gt = load_gt(Path(str(prefix) + ".gt.tsv"))
    report = json.loads(fc_json.read_text())
    metrics = report["metrics"]
    polyx = report["polyx"]
    insert_hist = report.get("insert_hist")

    pairs = len(gt)
    reads = 2 * pairs  # both mates are QC'd together

    # Ground-truth quantities.
    gt_adapter_pairs = sum(1 for r in gt if r["r1_adapter"] > 0 or r["r2_adapter"] > 0)
    gt_tail = {}
    for k in (5, 10):
        gt_tail[k] = sum((r["r1_t"] == k) + (r["r2_t"] == k) for r in gt)
    # Insert mean over the OVERLAP-MEASURABLE WINDOW: read_len < insert <=
    # 2*read_len - ANCHOR (the documented measurement ceiling, ANCHOR = 12).
    # The window is the support BOTH sides can express: the tool reports a
    # histogram of integer inserts and therefore cannot reproduce the plan's
    # original per-pair ceiling (b1 + b2 - ANCHOR, where b = L - h - t varies
    # read by read). Comparing the tool's whole bin range against that per-pair
    # set mixes two supports and reads as a ~2.2% recovery error that is really
    # a basis difference, so both numbers are computed and both are printed.
    ins_lo, ins_hi = read_len, 2 * read_len - ANCHOR
    gt_window = [r["insert"] for r in gt if ins_lo < r["insert"] <= ins_hi]
    gt_insert_mean = sum(gt_window) / len(gt_window) if gt_window else 0.0
    gt_meas_perpair = [
        r["insert"]
        for r in gt
        if r["insert"] > read_len
        and r["insert"] <= (read_len - r["r1_h"] - r["r1_t"])
        + (read_len - r["r2_h"] - r["r2_t"])
        - ANCHOR
    ]
    gt_perpair_mean = (sum(gt_meas_perpair) / len(gt_meas_perpair)
                       if gt_meas_perpair else 0.0)
    fc_bins = {int(k): v for k, v in (insert_hist or {}).get("bins", {}).items()}
    fc_window = {i: c for i, c in fc_bins.items() if ins_lo < i <= ins_hi}
    fc_insert_mean = (sum(i * c for i, c in fc_window.items()) / sum(fc_window.values())
                      if fc_window else 0.0)

    # 5' head convolution model, from the profile's head histogram.
    head_pmf = histogram_to_pmf(prof["head"]["R1"])
    pg = gc / 2.0
    p_head0 = head_pmf.get(0, 0.0)

    def head_model(k: int) -> float:
        p_run = pg ** k * (1.0 - pg)  # P(biological leading run == k)
        return head_pmf.get(k, 0.0) + p_head0 * p_run

    # fastcheck's per-position mean Q comes from the positions TSV (the JSON
    # report carries no per-position quality).
    tsv_lines = pos_tsv.read_text().splitlines()
    tsv_header = tsv_lines[0].split("\t")
    mq_col = tsv_header.index("mean_q")
    pos_mean_q = [float(line.split("\t")[mq_col]) for line in tsv_lines[1:]]

    # Ground-truth per-position mean Q: the profile's flat curve.
    pos_q_r1 = prof["pos_q"]["R1"]

    def p3g(k):
        return polyx["3p"]["G"].get(str(k), 0)

    def p5g(k):
        return polyx["5p"]["G"].get(str(k), 0)

    checks = []  # (name, ground-truth string, fastcheck string, deviation string, tolerance, passed)

    def add(name, gt_v, got_v, dev_str, tol, ok):
        checks.append((name, str(gt_v), str(got_v), dev_str, tol, ok))

    # 1. per-position mean Q (max absolute deviation, in Phred, over the positions).
    dev = max(abs(v - g) for v, g in zip(pos_mean_q, pos_q_r1))
    add("per-position mean Q (max dev)", "30.0", f"{max(pos_mean_q):.4f} max",
        f"{dev:.4f}", "<= 0.05", dev <= 0.05)

    # 2. whole-read GC.
    got_gc = metrics["gc_fraction"]
    add("whole-read GC", gc, f"{got_gc:.4f}", f"{abs(got_gc - gc):.4f}",
        "+/- 0.03", abs(got_gc - gc) <= 0.03)

    # 3. adapter rate (pairs).
    got_ad = insert_hist["measured_from_adapter"]
    add("adapter rate (pairs)", gt_adapter_pairs, got_ad,
        f"{abs(got_ad - gt_adapter_pairs)}",
        f"+/- 1% of {pairs}", abs(got_ad - gt_adapter_pairs) <= 0.01 * pairs)

    # 4-5. 3' G-run histogram, bins 5 and 10 (per read, both mates).
    for k in (5, 10):
        g, o = gt_tail[k], p3g(k)
        rel = abs(o - g) / g if g else math.inf
        add(f"3' G-run histogram, bin {k}", g, o, f"{rel * 100:.3f}%",
            "+/- 4%", g > 0 and rel <= 0.04)

    # 6-8. 5' head histogram, bins 1..3, against the convolution model.
    for k in (1, 2, 3):
        model = head_model(k) * reads
        o = p5g(k)
        rel = abs(o - model) / model if model else math.inf
        add(f"5' head histogram, bin {k} (convolution)", f"{model:.1f}", o,
            f"{rel * 100:.3f}%", "+/- 1.5%", model > 0 and rel <= 0.015)

    # 9. insert mean: both sides restricted to the SAME overlap window (above).
    rel = abs(fc_insert_mean - gt_insert_mean) / gt_insert_mean if gt_insert_mean else math.inf
    add("insert mean (overlap window)", f"{gt_insert_mean:.3f}", f"{fc_insert_mean:.3f}",
        f"{rel * 100:.3f}%", "+/- 2%", gt_insert_mean > 0 and rel <= 0.02)
    print(
        f"  note: mixed-support comparison for reference — tool's full bin range "
        f"{metrics['insert_mean']:.3f} vs the per-pair measurable set "
        f"{gt_perpair_mean:.3f} = "
        f"{abs(metrics['insert_mean'] - gt_perpair_mean) / gt_perpair_mean * 100:.3f}%"
        f" (a basis difference: the tool's bins include inserts below the read "
        f"length, which the per-pair set excludes by construction)"
    )

    # 10. duplication level is a documented LOWER BOUND on the PCR duplicate rate.
    got_dup = metrics["duplication_level"]
    add("duplication level (lower bound)", f"0..{dup_rate}", f"{got_dup:.4f}",
        "-", f"0 <= level <= {dup_rate}", 0.0 <= got_dup <= dup_rate)

    # Report.
    width = max(len(c[0]) for c in checks)
    print(f"ground-truth recovery: {pairs} pairs, read_len {read_len}, seed {args.seed}")
    print(f"{'check'.ljust(width)}  {'ground truth':>20}  {'fastcheck':>12}  {'dev':>9}  tolerance")
    failures = []
    for name, gt_v, got_v, dev_s, tol, ok in checks:
        mark = "ok " if ok else "FAIL"
        print(f"{name.ljust(width)}  {gt_v:>20}  {got_v:>12}  {dev_s:>9}  {tol}  [{mark}]")
        if not ok:
            failures.append((name, gt_v, got_v, dev_s, tol))

    if failures:
        print(f"\nGROUND TRUTH FAILED: {len(failures)} of {len(checks)} checks failed:",
              file=sys.stderr)
        for name, gt_v, got_v, dev, tol in failures:
            print(f"  {name}: ground truth {gt_v}, fastcheck {got_v} (tolerance {tol})", file=sys.stderr)
        if not args.keep and args.workdir is None:
            shutil.rmtree(work, ignore_errors=True)
        return 1

    print(f"\nGROUND TRUTH OK: all {len(checks)} checks within tolerance")
    if not args.keep and args.workdir is None:
        shutil.rmtree(work, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
