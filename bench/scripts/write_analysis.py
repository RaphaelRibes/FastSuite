#!/usr/bin/env python3
"""Write the per-config retention + delta table from a decomposition TSV.

Reads bench/results/phase-d/trim-report-flags-<date>.tsv and emits a NEW
timestamped TSV (never overwrites). Every number is computed here from the
decomposition, so nothing is transcribed by hand.
"""
import csv, sys

src = sys.argv[1]
out = sys.argv[2]
rows = list(csv.DictReader(open(src), delimiter="\t"))
per = {(r["tool"], int(r["tol"])): r for r in rows}

ORDER = ["fastprep-hd0", "fp-3p6", "fp-n", "fp-3p6-n", "fastprep-hd4"]
LABEL = {"fastprep-hd0": "a_baseline_hd0",
         "fp-3p6": "b_plus_min-run-3p_6",
         "fp-n": "c_plus_trim-terminal-n",
         "fp-3p6-n": "d_plus_both",
         "fastprep-hd4": "control_hd4"}
INT = ("exact", "overtrim", "undertrim_recoverable", "undertrim_irreducible",
       "dropped", "unaligned", "scored")

hdr = ["config", "tol"] + list(INT) + [
    "total", "retention_pct", "exact_rate_on_kept_pct",
    "d_exact", "d_overtrim", "d_undertrim_total", "d_dropped",
    "exact_rate_delta_pp_same_kept_basis"]
outrows = []
for tol in (0, 2):
    base = {k: int(per[("fastprep-hd0", tol)][k]) for k in INT}
    for tag in ORDER:
        r = per[(tag, tol)]
        s = {k: int(r[k]) for k in INT}
        tot = s["scored"] + s["dropped"]
        und = s["undertrim_recoverable"] + s["undertrim_irreducible"]
        b_und = base["undertrim_recoverable"] + base["undertrim_irreducible"]
        outrows.append([
            LABEL[tag], tol] + [s[k] for k in INT] + [
            tot, f"{100 * s['scored'] / tot:.4f}",
            f"{100 * s['exact'] / s['scored']:.4f}",
            s["exact"] - base["exact"], s["overtrim"] - base["overtrim"],
            und - b_und, s["dropped"] - base["dropped"],
            f"{100 * (s['exact'] - base['exact']) / base['scored']:+.4f}"])

with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(hdr)
    w.writerows(outrows)
print(f"wrote {out} ({len(outrows)} rows)")
print(f"source: {src}")
