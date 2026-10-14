#!/usr/bin/env python3
"""C1 aggregate — sweep rows -> pilot.tsv (winner per dataset).

Every sweep row is one (dataset, variant, config). A winner is chosen per
dataset on the MEAN score across the 3 variants, because a config that only
wins on one seed has not demonstrated it trims this data well — it got lucky.
The variants are three independent draws from the same fitted profile, so the
mean across them is the honest estimator.

The shipped default (config_id 0) is reported alongside every winner so the
memo can state, per dataset, whether the sweep beat fastprep's own default at
all — the number the pilot is really about.

Outputs:
  results/pilot.tsv    one row per dataset: the winner, the default, the gap
  results/winners.tsv  just the winners, the format fastp_compare.py consumes

Usage:
  aggregate.py --rows results/sweep_rows.tsv --outdir results/
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import sys
from collections import defaultdict

WINNER_FIELDS = ["domain", "run_accession", "winner_config_id", "winner_score",
                 "winner_keep_recall", "winner_keep_precision", "winner_over_trim",
                 "winner_false_keep", "winner_trim_mae",
                 "default_config_id", "default_score", "default_keep_recall",
                 "default_keep_precision", "default_over_trim", "default_false_keep",
                 "n_variants", "n_configs", "beats_default", "gap_vs_default",
                 "score_spread_across_variants"]

FLOAT_COLS = {"winner_score", "winner_keep_recall", "winner_keep_precision",
              "winner_over_trim", "winner_false_keep", "winner_trim_mae",
              "default_score", "default_keep_recall", "default_keep_precision",
              "default_over_trim", "default_false_keep", "gap_vs_default",
              "score_spread_across_variants"}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # (domain, acc, config_id) -> list of per-variant metric dicts
    by_cfg: dict[tuple[str, str, int], list[dict]] = defaultdict(list)
    with open(args.rows, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["domain"], r["run_accession"], int(r["config_id"]))
            by_cfg[key].append(r)

    # group configs by dataset
    by_ds: dict[tuple[str, str], dict[int, list[dict]]] = defaultdict(dict)
    for (domain, acc, cid), rows in by_cfg.items():
        by_ds[(domain, acc)][cid] = rows

    def mean(rows, col):
        vals = []
        for r in rows:
            try:
                vals.append(float(r[col]))
            except (KeyError, ValueError):
                pass
        return statistics.fmean(vals) if vals else None

    def config_summary(cid: int, rows: list[dict]) -> dict:
        scores = [float(r["score"]) for r in rows if r.get("score") not in (None, "")]
        sc = mean(rows, "score")
        # the effective min_len is per-config, not a mean: all variants of one
        # config share the same read_len, so the first row's value is the value.
        ml_eff = next((r.get("min_len_effective", "") for r in rows if r.get("min_len_effective")), "")
        return dict(
            config_id=cid,
            score=sc if sc is not None else float("-inf"),
            keep_recall=mean(rows, "keep_recall"),
            keep_precision=mean(rows, "keep_precision"),
            over_trim=mean(rows, "over_trim_rate"),
            false_keep=mean(rows, "false_keep_rate"),
            trim_mae=mean(rows, "trim_mae"),
            spread=(max(scores) - min(scores)) if len(scores) > 1 else 0.0,
            min_len_effective=ml_eff,
            params={k: rows[0].get(k, "") for k in
                    ("min_run", "mercy", "phred_sat", "gate_shape", "min_len",
                     "mean_q", "adapter_err", "min_overlap", "r2_reverse", "interior")},
        )

    pilot_rows = []
    for (domain, acc), cfgs in sorted(by_ds.items()):
        summaries = [config_summary(cid, rows) for cid, rows in sorted(cfgs.items())]
        summaries = [s for s in summaries if s["score"] is not None]
        if not summaries:
            continue
        winner = max(summaries, key=lambda s: s["score"])
        default = next((s for s in summaries if s["config_id"] == 0), summaries[0])
        n_variants = len({r["variant"] for rows in cfgs.values() for r in rows})

        row = {
            "domain": domain, "run_accession": acc,
            "winner_config_id": winner["config_id"],
            "winner_score": winner["score"],
            "winner_keep_recall": winner["keep_recall"],
            "winner_keep_precision": winner["keep_precision"],
            "winner_over_trim": winner["over_trim"],
            "winner_false_keep": winner["false_keep"],
            "winner_trim_mae": winner["trim_mae"],
            "default_config_id": default["config_id"],
            "default_score": default["score"],
            "default_keep_recall": default["keep_recall"],
            "default_keep_precision": default["keep_precision"],
            "default_over_trim": default["over_trim"],
            "default_false_keep": default["false_keep"],
            "n_variants": n_variants,
            "n_configs": len(summaries),
            "beats_default": int(winner["score"] > default["score"]),
            "gap_vs_default": winner["score"] - default["score"],
            "score_spread_across_variants": winner["spread"],
        }
        row.update({f"winner_{k}": v for k, v in winner["params"].items()})
        row["winner_min_len_effective"] = winner.get("min_len_effective", "")
        pilot_rows.append(row)

    # pilot.tsv: the summary columns first, the winning argv's params after
    pilot_cols = WINNER_FIELDS + [f"winner_{k}" for k in
                                  ("min_run", "mercy", "phred_sat", "gate_shape", "min_len",
                                   "mean_q", "adapter_err", "min_overlap", "r2_reverse", "interior")]
    pilot_path = os.path.join(args.outdir, "pilot.tsv")
    with open(pilot_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pilot_cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in pilot_rows:
            out = dict(r)
            for c in FLOAT_COLS:
                if isinstance(out.get(c), float):
                    out[c] = f"{out[c]:.4f}"
            w.writerow(out)

    # winners.tsv: exactly the columns fastp_compare.py consumes.
    # min_len_effective is carried through so C3 replays the winner with the
    # argv the sweep actually ran, not a re-derivation of it.
    winners_path = os.path.join(args.outdir, "winners.tsv")
    win_cols = ["domain", "run_accession", "config_id"] + \
               ["min_run", "mercy", "phred_sat", "gate_shape", "min_len",
                "mean_q", "adapter_err", "min_overlap", "r2_reverse", "interior",
                "min_len_effective"]
    with open(winners_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=win_cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in pilot_rows:
            out = {"domain": r["domain"], "run_accession": r["run_accession"],
                   "config_id": r["winner_config_id"]}
            out.update({k: r.get(f"winner_{k}", "") for k in win_cols[3:]})
            # the effective value comes from the sweep row of the winning config
            out["min_len_effective"] = r.get("winner_min_len_effective", "")
            w.writerow(out)

    n_beat = sum(r["beats_default"] for r in pilot_rows)
    print(f"aggregate: {len(pilot_rows)} datasets -> {pilot_path}")
    print(f"aggregate: winner beat the shipped default on {n_beat}/{len(pilot_rows)} datasets")
    print(f"aggregate: winners -> {winners_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
