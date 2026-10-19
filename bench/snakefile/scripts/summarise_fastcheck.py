#!/usr/bin/env python3
"""Summarise a fastcheck QC report into one row of the benchmark table.

Declared as an input of the fastcheck rule so editing it forces a rerun rather
than silently reusing a stale result (a script referenced only from `shell:` is
invisible to Snakemake).

Usage:
  summarise_fastcheck.py REPORT.json [--out TABLE.tsv]

Reads the `fastsuite.report.v1` JSON that `fastcheck --json` writes and emits a
tab-separated row: dataset id, the measurements the papers quote, and the file
size the run was charged with.
"""
import argparse
import json
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("report")
    ap.add_argument("--out", default=None)
    ap.add_argument("--ds", default=None, help="dataset id (domain/accession)")
    args = ap.parse_args()

    with open(args.report) as fh:
        rep = json.load(fh)

    # The report is a flat dict of counts; the exact keys depend on the
    # modules that ran, so emit whatever is present rather than assuming.
    def g(*names):
        for n in names:
            if n in rep:
                return rep[n]
        return ""

    row = {
        "dataset": args.ds or "",
        "reads": g("reads", "read_count", "total_reads"),
        "bases": g("bases", "base_count"),
        "mean_len": g("mean_len", "mean_length"),
        "q20": g("q20", "q20_pct"),
        "q30": g("q30", "q30_pct"),
        "gc": g("gc_pct", "gc"),
        "n_pct": g("n_pct", "n_fraction"),
        "adapter_pct": g("adapter_pct", "adapter"),
        "report_version": rep.get("version", rep.get("schema", "")),
    }
    cols = list(row)
    out = args.out
    write_header = out is not None and not os.path.exists(out)
    if out is None:
        print("\t".join(cols))
        print("\t".join(str(row[c]) for c in cols))
        return 0
    with open(out, "a") as fh:
        if write_header:
            fh.write("\t".join(cols) + "\n")
        fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
