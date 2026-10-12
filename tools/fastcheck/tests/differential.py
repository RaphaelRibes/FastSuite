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
        got_rows = [l.split("\t") for l in tsv.read_text().splitlines()]

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

    # 3. Per-position table, every column, read BY NAME from the header so adding a
    #    column (g_frac, polyx3_frac, ...) cannot break the comparison.
    ref_rows = ref["positions"]
    header = got_rows[0]
    col = {name: i for i, name in enumerate(header)}
    need = ["pos", "a", "c", "g", "t", "n", "depth", "mean_q", "q10", "q50", "q90"]
    missing = [n for n in need if n not in col]
    if missing:
        failures.append(f"positions TSV missing columns: {missing}")
    if len(got_rows) - 1 != len(ref_rows):
        failures.append(
            f"positions: fastcheck has {len(got_rows) - 1} rows, reference {len(ref_rows)}"
        )
    elif not missing:
        for gr, rr in zip(got_rows[1:], ref_rows):
            g = lambda name: gr[col[name]]
            for name in ("a", "c", "g", "t", "n", "depth"):
                got_i, want_i = int(g(name)), rr[name]
                if got_i != want_i:
                    failures.append(f"pos {g('pos')} {name}: fastcheck={got_i} reference={want_i}")
            if abs(float(g("mean_q")) - rr["mean_q"]) > 1e-3:
                failures.append(
                    f"pos {g('pos')} mean_q: fastcheck={g('mean_q')} reference={rr['mean_q']:.6f}"
                )
            for name in ("q10", "q50", "q90"):
                if abs(float(g(name)) - rr[name]) > 1e-9:
                    failures.append(f"pos {g('pos')} {name}: fastcheck={g(name)} reference={rr[name]}")

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
