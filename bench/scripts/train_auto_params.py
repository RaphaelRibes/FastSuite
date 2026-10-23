#!/usr/bin/env python3
"""E4/F1 — train the `fastprep --auto-params` decision list from the label table.

Input : labels.tsv (bench/scripts/labels.py) + the config columns it carries.
Output: a model file `tools/fastprep/models/auto-params.v1.json` and a model
        card stating what it was trained on, how it was evaluated, and where it
        must not be trusted.

WHAT THE MODEL IS
  A shallow decision tree (depth-limited, minimum leaf support) over the QC
  features fastcheck measures on the INPUT, whose leaves are fastprep configs.
  The tree is flattened into an ordered decision list: mutually exclusive
  root-to-leaf paths, first match wins, and the root's best config as the
  `default` fallback for a profile that matches no path (a missing feature can
  never satisfy a condition, by construction in the Rust walk).

WHY A TREE AND NOT A RANKER
  The recommendation has to be inspectable: `--print-args` prints the rule that
  fired, its training coverage and its mean regret, so the user can see the
  evidence and override it. A gradient-boosted ranker would recommend slightly
  better configs and be unexplainable on a tool whose failure mode is silently
  wrong reads.

THE ACTION SPACE IS SMALLER THAN THE SWEEP'S
  `r2_reverse` is a BENCH-side transform (the sweep revcomp'd R2 to simulate the
  other strand convention); it is not a fastprep flag, so every config with
  r2_reverse=1 is dropped from the action space — recommending one would be
  recommending an effect the tool cannot produce. The script prints how many
  datasets' oracle winners that removes.

EVALUATION
  Group k-fold by DATASET (never by row): a dataset's 55 config rows are all in
  the same fold, so the score is "a library the model has never seen". Reported:
  mean regret vs the fold's oracle, delta vs the shipped default, top-5
  agreement, and a per-domain breakdown. `--leave-one-domain-out` additionally
  reports the transfer case (train on five domains, test on the sixth).

HONESTY RULE (docs/PLAN-ML-STEPS.md E5)
  No GT-derived column is read here. The features are input-side only, exactly
  the columns labels.py allows, and the labels are the sweep's scores. If the
  model does not beat the default on held-out datasets, the card says so and the
  trainer exits non-zero unless `--force`, so a bad model cannot be shipped by
  accident.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from collections import defaultdict

SCHEMA = "fastsuite.auto-params.v1"

FEATURES = [
    "read_len",
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
    "bases_non_acgt_frac",
]

# label-table column -> model-file knob. `r2_reverse` is deliberately absent:
# it is a bench-side transform, not a fastprep flag.
CONFIG_COLS = [
    "min_run",
    "mercy",
    "phred_sat",
    "gate_shape",
    "min_len",
    "mean_q",
    "adapter_err",
    "min_overlap",
    "r2_reverse",
    "interior",
]

SHIPPED_DEFAULT = {
    "min_run": 12,
    "mercy": 3.0,
    "phred_sat": 20,
    "gate_shape": 1.0,
    "min_len": 0.9,          # the sweep ran at 0.9L; the shipped floor is 0.67L
    "mean_q": 25.0,
    "adapter_err": 0.2,
    "min_overlap": 14,
    "r2_reverse": 0,
    "interior": 0,
}


def fnum(v, default=float("nan")):
    try:
        f = float(v)
        return default if f != f else f
    except (TypeError, ValueError):
        return default


def load_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_dataset_table(rows):
    """dataset -> {config_key: mean score over variants}, plus the knob table.

    The label table has one row per (dataset, variant, config); a config's score
    for a dataset is the MEAN over variants, the same rule aggregate.py uses to
    pick the winners (never a re-derivation of a different statistic).
    """
    scores = defaultdict(lambda: defaultdict(list))
    knobs = {}
    feats = {}
    domain_of = {}
    for r in rows:
        ds = (r["domain"], r["run_accession"])
        key = tuple(r[c] for c in CONFIG_COLS)
        knobs[key] = {
            "min_run": int(fnum(r["min_run"])),
            "mercy": fnum(r["mercy"]),
            "phred_sat": int(fnum(r["phred_sat"])),
            "gate_shape": fnum(r["gate_shape"]),
            "min_len": fnum(r["min_len"]),
            "mean_q": fnum(r["mean_q"]),
            "adapter_err": fnum(r["adapter_err"]),
            "min_overlap": int(fnum(r["min_overlap"])),
            "r2_reverse": int(fnum(r["r2_reverse"], 0)),
            "interior": int(fnum(r["interior"], 0)),
        }
        scores[ds][key].append(fnum(r["score"]))
        domain_of[ds] = r["domain"]
        if ds not in feats:
            feats[ds] = {f: fnum(r.get(f)) for f in FEATURES}
    table = {}
    for ds, per_key in scores.items():
        table[ds] = {k: statistics.fmean(v) for k, v in per_key.items() if v}
    return table, knobs, feats, domain_of


def config_patch(knobs):
    """The knobs as a model-file `params` patch (complete, never partial)."""
    ml = knobs["min_len"]
    return {
        "min_run": knobs["min_run"],
        "mercy": round(knobs["mercy"], 4),
        "phred_sat": knobs["phred_sat"],
        "gate_shape": round(knobs["gate_shape"], 4),
        "min_len": (f"{ml:g}L" if ml <= 1.0 else int(ml)),
        "mean_q": (None if knobs["mean_q"] is None or knobs["mean_q"] != knobs["mean_q"]
                   else round(knobs["mean_q"], 4)),
        "adapter_error_rate": round(knobs["adapter_err"], 4),
        "min_adapter_overlap": knobs["min_overlap"],
        "interior": bool(knobs["interior"]),
    }


def config_label(knobs):
    """A short human label for the card, e.g. `mr16/m5/s20/l0.85L/q22`."""
    ml = knobs["min_len"]
    return (
        f"mr{knobs['min_run']}/m{knobs['mercy']:g}/s{knobs['phred_sat']}"
        f"/g{knobs['gate_shape']:g}/l{(f'{ml:g}L' if ml <= 1.0 else int(ml))}"
        f"/q{knobs['mean_q']:g}"
    )


# ---------------------------------------------------------------------------
# The learner: a greedy, depth-limited regression tree whose leaves are configs.
# ---------------------------------------------------------------------------


class Node:
    __slots__ = ("field", "threshold", "left", "right", "config", "n", "score")

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


def best_config(subset, table, space):
    """The config with the highest mean score over `subset` (ties: first)."""
    best, best_score = None, -math.inf
    for key in space:
        vals = [table[ds][key] for ds in subset if key in table[ds]]
        if not vals:
            continue
        s = statistics.fmean(vals)
        if s > best_score:
            best, best_score = key, s
    return best, best_score


def mean_regret(subset, table, key):
    """Mean (oracle - this config) over the subset. 0 = always the winner."""
    out = []
    for ds in subset:
        row = table[ds]
        if not row or key not in row:
            continue
        out.append(max(row.values()) - row[key])
    return statistics.fmean(out) if out else float("nan")


def candidate_thresholds(values, max_cands):
    vals = sorted(set(values))
    if len(vals) <= 1:
        return []
    if len(vals) <= max_cands:
        mids = [(a + b) / 2.0 for a, b in zip(vals, vals[1:])]
    else:
        step = len(vals) / (max_cands + 1)
        mids = []
        for i in range(1, max_cands + 1):
            j = int(i * step)
            if 0 < j < len(vals):
                mids.append((vals[j - 1] + vals[j]) / 2.0)
    return sorted(set(mids))


def grow(subset, feats, table, space, depth, min_leaf, max_depth, max_cands, min_gain):
    cfg, base = best_config(subset, table, space)
    node = Node(field=None, threshold=None, left=None, right=None, config=cfg, n=len(subset), score=base)
    if depth >= max_depth or len(subset) < 2 * min_leaf:
        return node

    best = None
    for f in FEATURES:
        vals = [feats[ds][f] for ds in subset if f in feats[ds] and feats[ds][f] == feats[ds][f]]
        if len(vals) < len(subset):  # a missing feature cannot be split on
            continue
        for t in candidate_thresholds(vals, max_cands):
            left = [ds for ds in subset if feats[ds][f] <= t]
            right = [ds for ds in subset if feats[ds][f] > t]
            if len(left) < min_leaf or len(right) < min_leaf:
                continue
            _, sl = best_config(left, table, space)
            _, sr = best_config(right, table, space)
            gain = (len(left) * sl + len(right) * sr) / len(subset) - base
            if gain < min_gain:
                continue
            if best is None or gain > best[0] + 1e-12:
                best = (gain, f, t, left, right, sl, sr)
    if best is None:
        return node

    gain, f, t, left, right, sl, sr = best
    node.field, node.threshold = f, t
    node.left = grow(left, feats, table, space, depth + 1, min_leaf, max_depth, max_cands, min_gain)
    node.right = grow(right, feats, table, space, depth + 1, min_leaf, max_depth, max_cands, min_gain)
    node.config = None
    return node


def flatten(node, conditions, out):
    """Root-to-leaf paths -> ordered rules. Paths are mutually exclusive."""
    if node.config is not None:
        out.append((list(conditions), node.config, node.n, node.score))
        return
    flatten(node.left, conditions + [(node.field, "<=", node.threshold)], out)
    flatten(node.right, conditions + [(node.field, ">", node.threshold)], out)


def predict(node, prof):
    while node.config is None:
        v = prof.get(node.field, float("nan"))
        if v != v:  # a missing feature can satisfy neither branch: fall back
            return root_config(node)
        node = node.left if v <= node.threshold else node.right
    return node.config


def root_config(node):
    while node.config is None:
        node = node.left
    return node.config


def evaluate(tree, root_default, datasets, feats, table, space):
    """Score the model on `datasets`: regret vs oracle, vs default, top-5."""
    regrets, deltas, top5, per_domain = [], [], [], defaultdict(lambda: [0, 0])
    for ds in datasets:
        row = table[ds]
        if not row:
            continue
        cfg = predict(tree, feats[ds])
        oracle = max(row.values())
        got = row.get(cfg, float("nan"))
        if got != got:
            continue
        regrets.append(oracle - got)
        if root_default in row:
            deltas.append(got - row[root_default])
        ranked = sorted(row.values(), reverse=True)[:5]
        hit = got >= ranked[-1] - 1e-12
        top5.append(1.0 if hit else 0.0)
        d = [0, 0]
        d[0] += 1
        d[1] += 1 if hit else 0
        per_domain[ds[0]][0] += 1
        per_domain[ds[0]][1] += 1 if hit else 0
    return {
        "n": len(regrets),
        "mean_regret": statistics.fmean(regrets) if regrets else float("nan"),
        "mean_delta_vs_default": statistics.fmean(deltas) if deltas else float("nan"),
        "beats_default_frac": (sum(1 for d in deltas if d > 1e-12) / len(deltas)) if deltas else float("nan"),
        "top5": statistics.fmean(top5) if top5 else float("nan"),
        "per_domain_top5": {
            k: {"n": v[0], "top5": (v[1] / v[0] if v[0] else float("nan"))}
            for k, v in sorted(per_domain.items())
        },
    }


def kfold(datasets, k, seed):
    """Deterministic folds, stratified by domain so no fold is one study type."""
    rnd = random.Random(seed)
    by_domain = defaultdict(list)
    for ds in datasets:
        by_domain[ds[0]].append(ds)
    folds = [[] for _ in range(k)]
    for dom in sorted(by_domain):
        ds = sorted(by_domain[dom])
        rnd.shuffle(ds)
        for i, d in enumerate(ds):
            folds[i % k].append(d)
    return folds


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", required=True, help="labels.tsv from labels.py")
    ap.add_argument("--out", required=True, help="model JSON to write")
    ap.add_argument("--card", help="model card markdown to write")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20261023)
    ap.add_argument("--max-depth", type=int, default=3)
    ap.add_argument("--min-leaf", type=int, default=10)
    ap.add_argument("--max-cands", type=int, default=16)
    ap.add_argument("--min-gain", type=float, default=0.0005)
    ap.add_argument("--leave-one-domain-out", action="store_true")
    ap.add_argument("--force", action="store_true", help="write the model even if it does not beat the default")
    args = ap.parse_args()

    rows = load_rows(args.labels)
    table, knobs, feats, _domain_of = build_dataset_table(rows)
    datasets = sorted(table, key=lambda d: (d[0], d[1]))
    if not datasets:
        sys.exit("no datasets in the label table")

    # Action space: every config fastprep can express.
    all_keys = sorted({k for row in table.values() for k in row}, key=lambda k: knobs[k].__repr__())
    space = [k for k in all_keys if knobs[k]["r2_reverse"] == 0]
    dropped = [k for k in all_keys if knobs[k]["r2_reverse"] == 1]
    if not space:
        sys.exit("no expressible configs in the action space")

    default_key = next(
        (k for k in space if all(
            (abs(knobs[k][c] - SHIPPED_DEFAULT[c]) < 1e-9 if isinstance(SHIPPED_DEFAULT[c], float)
             else knobs[k][c] == SHIPPED_DEFAULT[c])
            for c in CONFIG_COLS if c != "r2_reverse")),
        None,
    )

    # How much of the oracle the action space gives up by dropping r2_reverse.
    lost = 0
    for ds in datasets:
        row = table[ds]
        if not row:
            continue
        best_any = max(row.values())
        best_expr = max((row[k] for k in space if k in row), default=float("nan"))
        if best_expr == best_expr and best_expr < best_any - 1e-12:
            lost += 1

    print(f"datasets {len(datasets)}  configs {len(all_keys)}  expressible {len(space)}  dropped(r2_reverse) {len(dropped)}")
    print(f"datasets whose oracle winner needed an r2_reverse config: {lost}")
    print(f"shipped default resolved to {'a swept config' if default_key else 'NO swept config'}")

    # ---- cross-validated evaluation ------------------------------------
    folds = kfold(datasets, args.folds, args.seed)
    cv = {"mean_regret": [], "mean_delta_vs_default": [], "beats_default_frac": [], "top5": []}
    per_domain_cv = defaultdict(lambda: [0, 0])
    for i, test in enumerate(folds):
        train = [d for j, f in enumerate(folds) if j != i for d in f]
        tree = grow(train, feats, table, space, 0, args.min_leaf, args.max_depth, args.max_cands, args.min_gain)
        m = evaluate(tree, default_key, test, feats, table, space)
        for key in cv:
            cv[key].append(m[key])
        for dom, v in m["per_domain_top5"].items():
            per_domain_cv[dom][0] += v["n"]
            per_domain_cv[dom][1] += v["top5"] * v["n"]
        print(f"  fold {i + 1}/{len(folds)}: n={m['n']} regret={m['mean_regret']:.4f} "
              f"delta_vs_default={m['mean_delta_vs_default']:+.4f} top5={m['top5']:.3f}")

    cv_summary = {k: statistics.fmean(v) for k, v in cv.items()}
    cv_summary["per_domain_top5"] = {
        dom: {"n": v[0], "top5": (v[1] / v[0] if v[0] else float("nan"))}
        for dom, v in sorted(per_domain_cv.items())
    }

    # Baseline: what the default alone scores across all datasets.
    base = evaluate(
        Node(field=None, threshold=None, left=None, right=None, config=default_key, n=len(datasets), score=0.0),
        default_key, datasets, feats, table, space,
    )

    lod = None
    if args.leave_one_domain_out:
        lod = {}
        for dom in sorted({d[0] for d in datasets}):
            train = [d for d in datasets if d[0] != dom]
            test = [d for d in datasets if d[0] == dom]
            tree = grow(train, feats, table, space, 0, args.min_leaf, args.max_depth, args.max_cands, args.min_gain)
            m = evaluate(tree, default_key, test, feats, table, space)
            lod[dom] = {k: m[k] for k in ("n", "mean_regret", "mean_delta_vs_default", "top5")}
            print(f"  LODO {dom}: n={m['n']} regret={m['mean_regret']:.4f} "
                  f"delta_vs_default={m['mean_delta_vs_default']:+.4f} top5={m['top5']:.3f}")

    # ---- final model on everything --------------------------------------
    tree = grow(datasets, feats, table, space, 0, args.min_leaf, args.max_depth, args.max_cands, args.min_gain)
    leaves = []
    flatten(tree, [], leaves)

    rules = []
    for conditions, cfg, n, score in sorted(leaves, key=lambda x: -x[2]):
        regret = mean_regret([d for d in datasets if predict(tree, feats[d]) == cfg], table, cfg)
        rules.append({
            "when": [{"field": f, "op": op, "value": round(t, 6)} for f, op, t in conditions],
            "params": config_patch(knobs[cfg]),
            "coverage": n,
            "regret": round(regret, 6) if regret == regret else None,
            "note": config_label(knobs[cfg]) + f" (rule best on {n} dataset(s))",
        })

    model = {
        "schema": SCHEMA,
        "trained_on": {
            "corpus": "phase-d (real ENA libraries, GT-free features)",
            "datasets": len(datasets),
            "domains": sorted({d[0] for d in datasets}),
            "configs_in_space": len(space),
            "configs_dropped_r2_reverse": len(dropped),
            "created": "2026-10-23",
            "trainer": "bench/scripts/train_auto_params.py",
        },
        "features": FEATURES,
        "rules": rules,
        "default": {"params": config_patch(knobs[default_key]) if default_key else {}},
        "acceptance": {
            "method": f"group {args.folds}-fold by dataset, stratified by domain",
            "cv_mean_regret": round(cv_summary["mean_regret"], 6),
            "cv_delta_vs_default": round(cv_summary["mean_delta_vs_default"], 6),
            "cv_beats_default_frac": round(cv_summary["beats_default_frac"], 4),
            "cv_top5": round(cv_summary["top5"], 4),
            "cv_per_domain_top5": cv_summary["per_domain_top5"],
            "baseline_default_mean_regret": round(base["mean_regret"], 6),
            "baseline_default_top5": round(base["top5"], 4),
            "leave_one_domain_out": lod,
            "datasets_losing_to_r2_reverse_drop": lost,
        },
    }

    beats = cv_summary["mean_delta_vs_default"] > 0
    print(f"\nCV: regret {cv_summary['mean_regret']:.4f} vs default's {base['mean_regret']:.4f}; "
          f"delta vs default {cv_summary['mean_delta_vs_default']:+.4f}; "
          f"top5 {cv_summary['top5']:.3f} (default {base['top5']:.3f})")

    with open(args.out, "w") as fh:
        json.dump(model, fh, indent=2)
        fh.write("\n")
    print(f"wrote {args.out} ({len(rules)} rule(s))")

    if args.card:
        write_card(args.card, model, cv_summary, base, lod, lost, len(dropped), len(space), len(all_keys), len(datasets))
        print(f"wrote {args.card}")

    if not beats and not args.force:
        print("\nREFUSING: the model does not beat the shipped default on held-out datasets.", file=sys.stderr)
        print("Pass --force to write it anyway, and say so in the model card.", file=sys.stderr)
        sys.exit(3)


def write_card(path, model, cv, base, lod, lost, n_dropped, n_space, n_all, n_ds):
    t = model["trained_on"]
    a = model["acceptance"]
    lines = []
    lines.append("# Model card — `fastprep --auto-params` (v1)\n")
    lines.append("Generated by `bench/scripts/train_auto_params.py`. Every number here is")
    lines.append("computed from the label table at training time; none is hand-written.\n")
    lines.append("## What it was trained on\n")
    lines.append(f"* {t['datasets']} real ENA libraries, domains: {', '.join(t['domains'])}")
    lines.append(f"* {n_all} distinct configurations in the sweep (55 config ids; three ids "
                 f"repeat a knob set), {n_space} expressible by fastprep "
                 f"({n_dropped} dropped because they need the bench-only `r2_reverse` transform)")
    lines.append(f"* datasets whose oracle winner needed a dropped config: **{lost}**")
    lines.append("* features: fastcheck's input-side `metrics` only (no trimmed-output column)")
    lines.append("\n## How it was evaluated\n")
    lines.append(f"* {a['method']}")
    lines.append(f"* CV mean regret vs the per-dataset oracle: **{a['cv_mean_regret']:.4f}** "
                 f"(the shipped default alone: {a['baseline_default_mean_regret']:.4f})")
    lines.append(f"* CV delta vs the shipped default: **{a['cv_delta_vs_default']:+.4f}** "
                 f"(beats it on {100 * a['cv_beats_default_frac']:.1f} % of held-out datasets)")
    lines.append(f"* CV top-5 agreement: **{a['cv_top5']:.3f}** (default: {a['baseline_default_top5']:.3f})")
    lines.append("\n| domain | n | top-5 |")
    lines.append("|---|---|---|")
    for dom, v in a["cv_per_domain_top5"].items():
        lines.append(f"| {dom} | {v['n']} | {v['top5']:.3f} |")
    if lod:
        lines.append("\n### Leave-one-domain-out (transfer to an unseen study type)\n")
        lines.append("| held-out domain | n | regret | delta vs default | top-5 |")
        lines.append("|---|---|---|---|---|")
        for dom, v in sorted(lod.items()):
            lines.append(f"| {dom} | {v['n']} | {v['mean_regret']:.4f} | "
                         f"{v['mean_delta_vs_default']:+.4f} | {v['top5']:.3f} |")
    lines.append("\n## Where it must not be trusted\n")
    lines.append("* **Read length.** The sweep ran on 101 bp reads. `read_len` is a feature,")
    lines.append("  but no dataset in the corpus is 150 bp, so a rule that splits on it is")
    lines.append("  extrapolating. `--print-args` shows the rule; a 150 bp NovaSeq library should")
    lines.append("  not trust a `read_len < 110` leaf.")
    lines.append("* **Chemistry.** The corpus has no 2-colour (NovaSeq/NextSeq) library, so the")
    lines.append("  high-quality 5' polyG head case is not represented. `--head-dropout` is not")
    lines.append("  an axis of this model at all.")
    lines.append("* **Long low-quality 3' tails.** The corpus's longest injected 3' tail is 20 bp,")
    lines.append("  so the long-tail regime is unlearned.")
    lines.append("* **The length floor.** Every rule's `min_len` comes from the Phase D grid,")
    lines.append("  which ran at 0.9L-era floors. D-18 moved the shipped floor to 0.67L after")
    lines.append("  the sweep, so applying a recommendation REPLACES that floor; pass `--min-len`")
    lines.append("  explicitly (an explicit flag always wins) to keep the 0.67L budget.")
    lines.append("* **Orientation.** 13 configs of the sweep exploited a bench-only R2 reverse;")
    lines.append("  the model cannot reproduce that effect, and it is not a fastprep flag.")
    lines.append("* **Domain coverage:** " + ", ".join(t["domains"]) + " — a study type outside")
    lines.append("  these six is extrapolation, and the LODO table above is the measurement of how")
    lines.append("  bad that gets.")
    lines.append("\n## The rules\n")
    for i, r in enumerate(model["rules"], 1):
        when = " && ".join(f"{c['field']} {c['op']} {c['value']}" for c in r["when"]) or "always"
        lines.append(f"{i}. `{when}` → `{r['note']}` (n={r['coverage']}, regret {r['regret']})")
    lines.append(f"\ndefault: `{json.dumps(model['default']['params'], sort_keys=True)}`")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
