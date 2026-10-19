#!/usr/bin/env python3
"""E1 — the label table: one row per (dataset, config), QC-side features, the
outcome as the label.

The training matrix this builds has a hard asymmetry baked into its shape:

  FEATURES  what fastcheck measures on the INPUT (the QC report's `metrics`),
            plus the config's own knobs — nothing that requires running the
            trimmer to know.
  LABEL     how well that config actually did, scored against fastsim's ground
            truth by sweep.py. This side needs the GT; the features side does
            not.

The whole point of the split is that at predict time only the left side
exists. A model that could see the right side before predicting has learned
nothing about generalisation — that is what bench/scripts/leak_test.py checks,
and it is why this script refuses to pass any GT-derived column through.

Two label conventions, both from the sweep rows (never a re-derivation):

  --label winner   the winner flag: 1 for the config the C1 aggregate picked
                   for this dataset (the pilot.tsv/winners.tsv winner, resolved
                   on the MEAN score across variants, the same rule
                   aggregate.py uses), 0 otherwise. Phase D's results feed the
                   exact same shape.
  --label quantile 1 for the top-Q fraction of configs by score within the
                   dataset (e.g. --label quantile --quantile 0.1), 0 otherwise.
                   The ranking view: what a per-dataset ranker should learn.

Inputs (all TSV, read at run time — never committed, never fabricated):
  --rows     sweep_rows.tsv      one row per (dataset, variant, config)
  --pilot    pilot.tsv           one row per dataset (the winner's config id)
  --qc-dir   DIR                 the fastcheck QC report per dataset, named
                                 <domain>/<run_accession>.json — the report of
                                 the REAL input the profile was fitted from
                                 (e.g. real.report.json). Phase D's equivalent.
  --gt       DIR (optional)      GT files <domain>/<run_accession>/sim.v*.gt.tsv
                                 — only used to carry the GT regime
                                 (full/dup rate) as a STUDY-LEVEL context
                                 feature, never as a per-row label. Omit and
                                 the label table is built without them.

Outputs:
  --out      labels.tsv          the training matrix
  --features features.txt        the feature column list (the contract for the
                                 model: these columns, in this order)

Usage:
  labels.py --rows results/sweep_rows.tsv --pilot results/pilot.tsv \\
            --qc-dir qc/ --out results/labels.tsv --features results/features.txt

The QC features are exactly the scalars fastcheck reports in `metrics`
(QC_FEATURES below), which is what `fastprep --auto-params <qc.json>` would
have in hand at predict time. Anything else would be a field the predictor
cannot see.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import defaultdict

# The QC fields a predictor can actually read at predict time: fastcheck's
# `metrics` block on the REAL input. These are the ONLY non-config columns
# allowed in the feature matrix. A field is here if and only if it appears in
# the fastcheck report of the dataset before trimming happens.
QC_FEATURES = [
    "read_len",          # mean_len (the profile's read length)
    "gc_fraction",
    "mean_read_mean_q",
    "reads_q20_frac",
    "reads_q30_frac",
    "polyx5_frac_gt_1",
    "polyx3_frac_gt_5",
    "head_g_frac",
    "head_dropout",
    "duplication_level",
    "kmer_top_frac",
    "bases_non_acgt_frac",   # derived: bases_non_acgt / bases, bounded [0,1]
]

# The config knobs: the axis the model chooses over. Every one is an input the
# user could pass to fastprep; none is measured on the output.
CONFIG_FEATURES = [
    "min_run", "mercy", "phred_sat", "gate_shape", "min_len", "mean_q",
    "adapter_err", "min_overlap", "r2_reverse", "interior",
]

# The outcome. score is the sweep's single winner number (sweep.py docstring);
# the rest decompose it so the ranker can see WHERE a config lost.
OUTCOME_FEATURES = [
    "score", "keep_precision", "keep_recall", "trim_mae", "over_trim_rate",
    "under_trim_rate", "false_keep_rate", "drop_accuracy",
]

# The label: the target. Exactly one of these is 1 per dataset row.
LABEL_COLS = {
    "winner": "is_winner",
    "quantile": "in_top",
}

FLOAT_COLS = QC_FEATURES + CONFIG_FEATURES + OUTCOME_FEATURES


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def to_float(v, default=0.0):
    try:
        f = float(v)
        return default if f != f else f  # NaN -> default
    except (TypeError, ValueError):
        return default


def load_qc_features(qc_dir: str) -> dict[tuple[str, str], dict[str, float]]:
    """fastcheck `metrics` per dataset, keyed by (domain, run_accession).

    One report per dataset: the report of the REAL reads the fitted profile
    came from. A dataset whose report is missing contributes no QC features —
    the row is DROPPED, because a label row without the predict-time features
    is exactly the unlearnable case the leak test bans (a model trained on
    rows whose features it could never see is guessing).
    """
    out: dict[tuple[str, str], dict[str, float]] = {}
    if not os.path.isdir(qc_dir):
        return out
    for domain in sorted(os.listdir(qc_dir)):
        dpath = os.path.join(qc_dir, domain)
        if not os.path.isdir(dpath):
            continue
        # Two layouts, both seen in the wild:
        #   <domain>/<acc>.json                 (a flat QC dir)
        #   <domain>/<acc>/real.report.json     (the campaign's work/ layout —
        #                                      the report of the REAL reads the
        #                                      fitted profile came from)
        entries = []
        for name in sorted(os.listdir(dpath)):
            sub = os.path.join(dpath, name)
            if os.path.isdir(sub):
                rep = os.path.join(sub, "real.report.json")
                if os.path.exists(rep):
                    entries.append((name, rep))
            elif name.endswith(".json"):
                entries.append((name[: -len(".json")], sub))
        for acc, path in entries:
            try:
                with open(path) as fh:
                    doc = json.load(fh)
            except (OSError, json.JSONDecodeError):
                continue
            m = doc.get("metrics", {})
            if not m:
                continue
            bases = to_float(m.get("bases"), 0.0)
            non_acgt = to_float(m.get("bases_non_acgt"), 0.0)
            feat = {
                "read_len": to_float(m.get("mean_len")),
                "gc_fraction": to_float(m.get("gc_fraction")),
                "mean_read_mean_q": to_float(m.get("mean_read_mean_q")),
                "reads_q20_frac": to_float(m.get("reads_q20_frac")),
                "reads_q30_frac": to_float(m.get("reads_q30_frac")),
                "polyx5_frac_gt_1": to_float(m.get("polyx5_frac_gt_1")),
                "polyx3_frac_gt_5": to_float(m.get("polyx3_frac_gt_5")),
                "head_g_frac": to_float(m.get("head_g_frac")),
                "head_dropout": to_float(m.get("head_dropout")),
                "duplication_level": to_float(m.get("duplication_level")),
                "kmer_top_frac": to_float(m.get("kmer_top_frac")),
                "bases_non_acgt_frac": (non_acgt / bases) if bases > 0 else 0.0,
            }
            out[(domain, acc)] = feat
    return out


def load_gt_regime(gt_dir: str) -> dict[tuple[str, str], dict[str, float]]:
    """Study-level GT context: the artefact regime the GT defines.

    Used ONLY as two study-level columns (gt_full_frac, gt_dup_frac), never as
    a per-row feature a row-level split could leak through: they are constant
    within a study, and the split is by study, so a test study's rows carry
    values no train study contributed. They describe the DATA, not the trim —
    at predict time a fastcheck report of the real input has no GT, so these
    two columns are study context, deliberately excluded from FEATURES below.
    """
    out: dict[tuple[str, str], dict[str, float]] = {}
    if not os.path.isdir(gt_dir):
        return out
    for domain in sorted(os.listdir(gt_dir)):
        dpath = os.path.join(gt_dir, domain)
        if not os.path.isdir(dpath):
            continue
        for acc in sorted(os.listdir(dpath)):
            apath = os.path.join(dpath, acc)
            if not os.path.isdir(apath):
                continue
            n = 0
            full = 0
            dup = 0
            for name in sorted(os.listdir(apath)):
                if not name.endswith(".gt.tsv"):
                    continue
                with open(os.path.join(apath, name), newline="") as fh:
                    rd = csv.reader(fh, delimiter="\t")
                    header = next(rd, None)
                    if not header or header[0] != "pair_id":
                        continue
                    for r in rd:
                        if not r:
                            continue
                        n += 1
                        full += int(r[8])
                        dup += int(r[9])
            if n:
                out[(domain, acc)] = {"gt_full_frac": full / n,
                                      "gt_dup_frac": dup / n}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", required=True, help="sweep_rows.tsv")
    ap.add_argument("--pilot", required=True, help="pilot.tsv (winners)")
    ap.add_argument("--qc-dir", required=True, help="per-dataset fastcheck reports")
    ap.add_argument("--out", required=True, help="labels.tsv")
    ap.add_argument("--features", help="write the feature column list here")
    ap.add_argument("--gt", help="GT dir (optional study-level context)")
    ap.add_argument("--label", choices=["winner", "quantile"], default="winner")
    ap.add_argument("--quantile", type=float, default=0.1,
                    help="top fraction within a dataset for --label quantile")
    args = ap.parse_args()

    rows = read_tsv(args.rows)
    pilot = read_tsv(args.pilot)
    qc = load_qc_features(args.qc_dir)
    gt = load_gt_regime(args.gt) if args.gt else {}

    # The winner per dataset: the config the aggregate picked (mean score
    # across variants), read from pilot.tsv so the label is the campaign's own
    # verdict, never a re-derivation that could drift from it.
    winners: dict[tuple[str, str], int] = {}
    for r in pilot:
        try:
            winners[(r["domain"], r["run_accession"])] = int(r["winner_config_id"])
        except (KeyError, ValueError):
            continue

    # group sweep rows per dataset, keep only datasets with a QC report.
    by_ds: dict[tuple[str, str], list[dict]] = defaultdict(list)
    skipped_no_qc: set[tuple[str, str]] = set()
    for r in rows:
        key = (r["domain"], r["run_accession"])
        if key not in qc:
            skipped_no_qc.add(key)
            continue
        by_ds[key].append(r)

    label_col = LABEL_COLS[args.label]
    feature_cols = QC_FEATURES + CONFIG_FEATURES
    out_cols = ["domain", "run_accession", "config_id"] + feature_cols + \
               OUTCOME_FEATURES + [label_col] + \
               (["gt_full_frac", "gt_dup_frac"] if gt else [])

    n_rows = 0
    n_datasets = 0
    missing_winner: set[tuple[str, str]] = set()
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for (domain, acc), ds_rows in sorted(by_ds.items()):
            qf = qc[(domain, acc)]
            n_datasets += 1
            win = winners.get((domain, acc))
            if args.label == "winner" and win is None:
                missing_winner.add((domain, acc))
            # per-dataset score list for the quantile label
            scores = [to_float(r.get("score"), float("-inf")) for r in ds_rows]
            if scores:
                k = max(1, int(len(scores) * args.quantile))
                cutoff = sorted(scores, reverse=True)[min(k, len(scores)) - 1]
            else:
                cutoff = float("-inf")
            for r in ds_rows:
                out = {"domain": domain, "run_accession": acc,
                       "config_id": r.get("config_id", "")}
                # QC features come from the report (constant per dataset);
                # config features come from the sweep row.
                out.update({c: f"{qf[c]:.6g}" for c in QC_FEATURES})
                out.update({c: f"{to_float(r.get(c)):.6g}" for c in CONFIG_FEATURES})
                out.update({c: f"{to_float(r.get(c)):.6g}" for c in OUTCOME_FEATURES})
                score = to_float(r.get("score"), float("-inf"))
                if args.label == "winner":
                    out[label_col] = 1 if win is not None and \
                        r.get("config_id") != "" and \
                        int(r["config_id"]) == win else 0
                else:
                    out[label_col] = 1 if score >= cutoff else 0
                if gt and (domain, acc) in gt:
                    out.update({c: f"{gt[(domain, acc)][c]:.6g}"
                                for c in ("gt_full_frac", "gt_dup_frac")})
                w.writerow(out)
                n_rows += 1

    if args.features:
        os.makedirs(os.path.dirname(os.path.abspath(args.features)), exist_ok=True)
        with open(args.features, "w") as fh:
            fh.write("\n".join(feature_cols) + "\n")

    print(f"labels: {n_rows} rows over {n_datasets} datasets -> {args.out}")
    if skipped_no_qc:
        print(f"labels: {len(skipped_no_qc)} datasets dropped (no QC report): "
              f"{sorted(skipped_no_qc)[:5]}{'...' if len(skipped_no_qc) > 5 else ''}",
              file=sys.stderr)
    if missing_winner:
        print(f"labels: WARNING {len(missing_winner)} datasets have no winner in "
              f"pilot.tsv: {sorted(missing_winner)[:5]}", file=sys.stderr)
    if args.features:
        print(f"labels: {len(feature_cols)} feature columns -> {args.features}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
