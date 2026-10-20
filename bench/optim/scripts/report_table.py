"""Emit the optimisation report's before/after table from a campaign result file.

Reads the aggregated `results.tsv` (one row per cell, `kind` = baseline | new |
comparator) and prints a Typst `table(...)` body for the six-tool summary the
report carries. One row per tool: its largest measured size, single thread, at
the tool's default compression, before and after.

The point of generating it rather than typing it is that the numbers cannot drift
from the file, and a missing cell shows up as a visible gap instead of a silently
remembered figure.

Usage: python3 report_table.py /path/to/results.tsv [--size 5000000] [--threads 1]
"""

import argparse
import sys


def load(path):
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def pick(rows, tool, size, threads):
    out = {}
    for r in rows:
        if r["tool"] != tool or r["kind"] not in ("baseline", "new"):
            continue
        if int(r["reads"]) != size or int(r["threads"]) != threads:
            continue
        if r.get("gz_level", "") not in ("", "6"):
            continue
        out[r["kind"]] = r
    return out.get("baseline"), out.get("new")


def fmt_wall(r):
    return f"{float(r['wall_s']):.2f}"


def fmt_rss(r):
    return f"{int(r['peak_rss_kb']) / 1024:.0f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--size", type=int, default=5_000_000)
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    rows = load(args.results)
    tools = ["fastcheck", "fastprep", "fastsim", "fastmaps", "fasthelper", "fastdedup"]

    print("  table(")
    print("    columns: (auto, auto, auto, auto, 1fr),")
    print("    table.header([*Outil*], [*Avant, s*], [*Après, s*], "
          "[*RSS avant / après, Mo*], [*Rapport de contrôle*]),")
    missing = []
    for tool in tools:
        b, n = pick(rows, tool, args.size, args.threads)
        if not b or not n:
            missing.append(tool)
            print(f"    [{tool}], [---], [---], [---], [cellule absente],")
            continue
        ratio = float(b["cpu_over_wall"]), float(n["cpu_over_wall"])
        check = (f"débit {float(b['wall_s']) / float(n['wall_s']):.2f}× ; "
                 f"cpu/mur {ratio[0]:.2f} → {ratio[1]:.2f}")
        print(f"    [{tool}], [{fmt_wall(b)}], [{fmt_wall(n)}], "
              f"[{fmt_rss(b)} / {fmt_rss(n)}], [{check}],")
    print(f"  ),")
    if missing:
        print(f"// MISSING CELLS: {', '.join(missing)}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
