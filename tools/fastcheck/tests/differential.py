#!/usr/bin/env python3
"""Differential test: fastcheck vs the independent Python reference.

Runs `fastcheck --json` and `reference_qc.py` on the same input and compares
every metric, the length histogram, the per-read mean-Q histogram and the
full per-position table. Exits non-zero with a printed diff on any mismatch.

Usage: differential.py FILE [--max-pos N]
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
FASTCHECK = HERE.parent / "target" / "release" / "fastcheck"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--max-pos", type=int, default=1000)
    ap.add_argument("--fastcheck", default=str(FASTCHECK))
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as td:
        tsv = Path(td) / "pos.tsv"
        js = Path(td) / "fastcheck.json"
        subprocess.run(
            [args.fastcheck, args.path, "--json", str(js), "--positions-tsv", str(tsv),
             "--max-read-len", str(args.max_pos), "-q"],
            check=True,
        )
        mine = json.loads(js.read_text())
        got_rows = [l.split("\t") for l in tsv.read_text().splitlines()[1:]]

    ref = json.loads(
        subprocess.run(
            [sys.executable, str(HERE / "reference_qc.py"), args.path, "--max-pos", str(args.max_pos)],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    )

    failures = []

    # 1. Scalar metrics, exact comparison (no tolerance: both are exact integer
    #    work, so a difference is a bug, not noise).
    #    `partial` is excluded because the reference is never budgeted.
    mine_m = mine["metrics"]
    for key, ref_val in ref["metrics"].items():
        if key == "partial":
            continue
        got = mine_m.get(key)
        if got is None:
            failures.append(f"metric {key}: missing from fastcheck")
        elif isinstance(ref_val, float):
            if abs(float(got) - ref_val) > 1e-9:
                failures.append(f"metric {key}: fastcheck={got} reference={ref_val}")
        elif got != ref_val:
            failures.append(f"metric {key}: fastcheck={got} reference={ref_val}")

    # 2. Histograms, key by key.
    for block in ("read_length_hist", "read_mean_q_hist"):
        a, b = mine[block], ref[block]
        if a != b:
            keys = sorted(set(a) | set(b), key=lambda k: int(k))
            for k in keys:
                if a.get(k) != b.get(k):
                    failures.append(f"{block}[{k}]: fastcheck={a.get(k)} reference={b.get(k)}")

    # 3. Per-position table, every column.
    ref_rows = ref["positions"]
    if len(got_rows) != len(ref_rows):
        failures.append(f"positions: fastcheck has {len(got_rows)} rows, reference {len(ref_rows)}")
    else:
        for gr, rr in zip(got_rows, ref_rows):
            pos, a, c, g, t, n, depth, mean_q, q10, q50, q90 = gr
            for name, got, want in (
                ("a", int(a), rr["a"]), ("c", int(c), rr["c"]), ("g", int(g), rr["g"]),
                ("t", int(t), rr["t"]), ("n", int(n), rr["n"]), ("depth", int(depth), rr["depth"]),
            ):
                if got != want:
                    failures.append(f"pos {pos} {name}: fastcheck={got} reference={want}")
            if abs(float(mean_q) - rr["mean_q"]) > 1e-3:
                failures.append(f"pos {pos} mean_q: fastcheck={mean_q} reference={rr['mean_q']:.6f}")
            for name, got in (("q10", q10), ("q50", q50), ("q90", q90)):
                if abs(float(got) - rr[name]) > 1e-9:
                    failures.append(f"pos {pos} {name}: fastcheck={got} reference={rr[name]}")

    if failures:
        print(f"DIFFERENTIAL FAILED ({len(failures)} differences) on {args.path}")
        for f in failures[:40]:
            print("  " + f)
        return 1

    print(
        f"differential OK on {args.path}: {mine_m['reads']} reads, "
        f"{len(ref['positions'])} positions, every metric identical to the reference"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
