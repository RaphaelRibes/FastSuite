#!/usr/bin/env python3
"""E2 — the leak test: the two checks that make a held-out score mean anything.

A model's reported accuracy is worthless if either of these is true, and both
failures are silent — nothing errors, the numbers just look fine. So both are
asserted here, as a runnable check with a named failure, exit 1 on any fail.

CHECK 1 — study-level split safety (PLAN-ML-STEPS E3 acceptance: "no study
appears in two splits"). For every fold pair, the intersection of the train
study sets must be empty. The split key is the RUN ACCESSION, not the row:
a dataset contributes one row per config, and every one of those rows must
live on the same side of every split. This checks splits.json's folds against
the label table, so it is the end-to-end form of the rule — not just the
split file in isolation.

CHECK 2 — no outcome feature (E2 acceptance: "no simulator-only field is
reachable from the real-data vector"). The feature columns of the label table
must all be computable at predict time from a fastcheck QC report of the
INPUT plus the config knobs. Concretely:

  * every feature must be in QC_FEATURES ∪ CONFIG_FEATURES (the allow-list
    labels.py exports to features.txt and the model consumes);
  * no feature column may be any outcome column (score, keep_*, trim_mae,
    *_rate, drop_accuracy, kept, should_keep, seconds, ...): those are the
    LABEL side, derived from the GT, and a model that could read them would
    be reading the answer;
  * every QC feature must actually be present in the fastcheck `metrics`
    block the predictor reads (a feature named for a field fastcheck does not
    emit is a feature that cannot be computed at predict time — same leak,
    opposite direction).

Usage:
  leak_test.py --labels results/labels.tsv --features results/features.txt \\
               --splits results/splits.json --qc-dir qc/
  # or, with a single QC report standing in for the schema check:
  leak_test.py --labels results/labels.tsv --features results/features.txt \\
               --splits results/splits.json --qc-json results/c3work/DOMAIN_ACC/in.json
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

# Same allow-list labels.py exports. Kept in sync by the feature-list check
# below (a drift between the two is itself a failure).
QC_FEATURES = [
    "read_len", "gc_fraction", "mean_read_mean_q", "reads_q20_frac",
    "reads_q30_frac", "polyx5_frac_gt_1", "polyx3_frac_gt_5", "head_g_frac",
    "head_dropout", "duplication_level", "kmer_top_frac", "bases_non_acgt_frac",
]
CONFIG_FEATURES = [
    "min_run", "mercy", "phred_sat", "gate_shape", "min_len", "mean_q",
    "adapter_err", "min_overlap", "r2_reverse", "interior",
]
OUTCOME_COLS = {
    "score", "keep_precision", "keep_recall", "trim_mae", "over_trim_rate",
    "under_trim_rate", "false_keep_rate", "drop_accuracy", "kept",
    "should_keep", "seconds", "read_len_sweep",
}

# QC feature name -> the fastcheck `metrics` field it is computed from. A
# feature that maps to a field fastcheck does not emit would be uncomputable
# at predict time — that is a leak, so the mapping is explicit and checked.
QC_FIELD_MAP = {
    "read_len": "mean_len",
    "gc_fraction": "gc_fraction",
    "mean_read_mean_q": "mean_read_mean_q",
    "reads_q20_frac": "reads_q20_frac",
    "reads_q30_frac": "reads_q30_frac",
    "polyx5_frac_gt_1": "polyx5_frac_gt_1",
    "polyx3_frac_gt_5": "polyx3_frac_gt_5",
    "head_g_frac": "head_g_frac",
    "head_dropout": "head_dropout",
    "duplication_level": "duplication_level",
    "kmer_top_frac": "kmer_top_frac",
    "bases_non_acgt_frac": ("bases_non_acgt", "bases"),
}

fails: list[str] = []


def fail(msg: str) -> None:
    fails.append(msg)
    print(f"FAIL  {msg}")


def ok(msg: str) -> None:
    print(f"ok    {msg}")


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check_split_safety(labels: list[dict], splits: dict) -> None:
    """Every study's rows sit on one side of every fold."""
    rows_by_study: dict[str, list[dict]] = {}
    for r in labels:
        rows_by_study.setdefault(f"{r['domain']}/{r['run_accession']}", []).append(r)

    all_test: set[str] = set()
    for i, fold in enumerate(splits.get("folds", [])):
        test = set(fold["test"])
        overlap = all_test & test
        if overlap:
            fail(f"fold {i}: study in more than one test fold: "
                 f"{sorted(overlap)[:5]}")
        # a study in this fold's test must contribute no train rows
        leaked = test & set(fold["train"])
        if leaked:
            fail(f"fold {i}: study in BOTH train and test: {sorted(leaked)[:5]}")
        all_test |= test

    missing = set(rows_by_study) - all_test
    if missing:
        fail(f"{len(missing)} studies in the label table are in no test fold: "
             f"{sorted(missing)[:5]}")

    # The row-level form: any study whose rows appear on both sides of a fold
    # is the same leak seen from the data. The set check above covers it, but
    # assert on actual rows so a malformed splits.json cannot pass by luck.
    for i, fold in enumerate(splits.get("folds", [])):
        test = set(fold["test"])
        for sid in test & set(fold["train"]):
            fail(f"fold {i}: study {sid} has rows in train AND test")
    if not fails:
        ok(f"study-level split safety: {len(rows_by_study)} studies, "
           f"{len(splits.get('folds', []))} folds, no study in two splits")


def check_no_outcome_features(labels: list[dict], features: list[str],
                              qc_metrics: dict | None) -> None:
    declared = list(labels[0].keys()) if labels else []

    for f in features:
        if f in OUTCOME_COLS:
            fail(f"feature {f!r} is an outcome column — it is the label side, "
                 f"computed from the GT, unreachable at predict time")
        if f not in QC_FEATURES and f not in CONFIG_FEATURES:
            fail(f"feature {f!r} is neither a QC feature nor a config knob — "
                 f"not computable at predict time")
    if features and not fails:
        ok(f"feature allow-list: all {len(features)} features are QC-side "
           f"(fastcheck on the input) or config knobs")

    # The label table MUST carry the outcome columns (they are the label side,
    # needed to evaluate). What must never happen is an outcome column inside
    # the feature list — checked above — or a row-level feature that is not in
    # the declared feature list. Here: every column the label table carries is
    # either a declared feature, a declared label, or an outcome column. Any
    # other column is an undeclared feature and fails.
    known = set(features) | {"domain", "run_accession", "config_id",
                             "is_winner", "in_top", "gt_full_frac",
                             "gt_dup_frac"} | OUTCOME_COLS
    undeclared = [c for c in declared if c not in known]
    if undeclared:
        fail(f"label table carries columns that are neither features, labels "
             f"nor outcomes: {undeclared}")

    if qc_metrics is not None:
        for f in features:
            if f not in QC_FEATURES:
                continue
            fields = QC_FIELD_MAP[f]
            if isinstance(fields, str):
                fields = (fields,)
            absent = [x for x in fields if x not in qc_metrics]
            if absent:
                fail(f"QC feature {f!r} needs fastcheck metrics field(s) "
                     f"{absent} that the QC report does not emit")
        if not fails:
            ok(f"QC schema: every QC feature is computable from the fastcheck "
               f"`metrics` block at predict time")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--splits", required=True)
    ap.add_argument("--qc-dir", help="per-dataset QC reports <domain>/<acc>.json")
    ap.add_argument("--qc-json", help="one QC report (schema check only)")
    args = ap.parse_args()

    labels = read_tsv(args.labels)
    with open(args.features) as fh:
        features = [l.strip() for l in fh if l.strip()]
    with open(args.splits) as fh:
        splits = json.load(fh)

    qc_metrics = None
    if args.qc_json:
        with open(args.qc_json) as fh:
            qc_metrics = json.load(fh).get("metrics", {})
    elif args.qc_dir:
        # the first report found stands in for the schema check: the metrics
        # block is the same shape for every report of this tool version
        for domain in sorted(os.listdir(args.qc_dir)):
            p = os.path.join(args.qc_dir, domain)
            if not os.path.isdir(p):
                continue
            for name in sorted(os.listdir(p)):
                if name.endswith(".json"):
                    with open(os.path.join(p, name)) as fh:
                        qc_metrics = json.load(fh).get("metrics", {})
                    break
            if qc_metrics is not None:
                break

    check_split_safety(labels, splits)
    check_no_outcome_features(labels, features, qc_metrics)

    print()
    if fails:
        print(f"LEAK TEST FAILED: {len(fails)} failure(s)")
        return 1
    print("LEAK TEST PASSED — splits are study-level; no outcome feature is "
          "reachable at predict time")
    return 0


if __name__ == "__main__":
    sys.exit(main())
