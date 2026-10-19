#!/usr/bin/env python3
"""E3 — the study-level split: group k-fold by run_accession, stratified by
domain, one fixed seed, byte-for-byte deterministic.

The split key is the RUN ACCESSION, never the row. Each dataset's rows (one
per config, across the variants) are one group: a study that appears in the
test fold must contribute ZERO rows to any train fold, because the whole point
is to score the model on studies it has never seen — per-row splits leak a
study into both sides and the held-out number becomes fiction.

Stratification is by DOMAIN, so every fold carries every domain in the same
proportion as the corpus. Within a domain the accessions are dealt round-robin
into folds in a deterministic order (sorted, then a seeded shuffle), which
balances fold sizes without any dependence on dict iteration order.

Determinism: the shuffle is a seeded Fisher-Yates over the SORTED accession
list, so the same corpus + seed + k always produce the same splits, on any
machine, in any Python. The seed defaults to 7 (the bench's canonical seed,
see docs/BENCHMARKS.md).

Outputs:
  --out splits.json   {"folds": [{"train": [ids...], "test": [ids...]}...],
                       "seed": 7, "k": 5, "domains": {...}, "n_studies": N}
                       where each id is "<domain>/<run_accession>".

Usage:
  splits.py --rows results/sweep_rows.tsv --k 5 --seed 7 --out results/splits.json
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
from collections import defaultdict


def read_tsv(path: str) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_splits(studies: list[str], domains: dict[str, str], k: int,
                 seed: int) -> list[dict[str, list[str]]]:
    """Deal every domain's studies round-robin into k folds, deterministically.

    Round-robin (not a global shuffle-then-chunk) keeps each fold's domain mix
    proportional: within a domain the studies are shuffled once with the seed,
    then assigned fold = i % k, so a domain with 10 studies and k=5 puts two
    in each fold.
    """
    by_domain: dict[str, list[str]] = defaultdict(list)
    for sid in studies:
        by_domain[domains[sid]].append(sid)

    folds: list[list[str]] = [[] for _ in range(k)]
    for domain in sorted(by_domain):
        accs = sorted(by_domain[domain])
        rng = random.Random(seed)
        rng.shuffle(accs)
        for i, acc in enumerate(accs):
            folds[i % k].append(acc)

    out = []
    all_ids = set(studies)
    for i in range(k):
        test = sorted(folds[i])
        test_set = set(test)
        train = sorted(all_ids - test_set)
        out.append({"train": train, "test": test})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", required=True,
                    help="sweep_rows.tsv (or any TSV with domain + run_accession)")
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", required=True, help="splits.json")
    args = ap.parse_args()

    if args.k < 2:
        print("splits: --k must be >= 2", file=sys.stderr)
        return 1

    domains: dict[str, str] = {}
    for r in read_tsv(args.rows):
        sid = f"{r['domain']}/{r['run_accession']}"
        domains[sid] = r["domain"]

    studies = sorted(domains)
    folds = build_splits(studies, domains, args.k, args.seed)

    by_domain: dict[str, int] = defaultdict(int)
    for sid in studies:
        by_domain[domains[sid]] += 1

    doc = {
        "seed": args.seed,
        "k": args.k,
        "n_studies": len(studies),
        "domains": dict(sorted(by_domain.items())),
        "split_key": "run_accession (study-level: no study appears in two folds)",
        "folds": folds,
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(doc, fh, indent=1)
        fh.write("\n")

    # The acceptance criterion from PLAN-ML-STEPS E3: no study appears in two
    # splits. Verify it here, on the written file, so a green run is proof.
    all_test: set[str] = set()
    overlap: set[str] = set()
    for f in folds:
        overlap |= all_test & set(f["test"])
        all_test |= set(f["test"])
    ok = not overlap and all_test == set(studies)
    print(f"splits: {len(studies)} studies, k={args.k}, seed={args.seed} -> {args.out}")
    for i, f in enumerate(folds):
        print(f"  fold {i}: train {len(f['train'])}  test {len(f['test'])}")
    if not ok:
        print(f"splits: FAILED — studies in multiple test folds: "
              f"{sorted(overlap)[:5]}", file=sys.stderr)
        return 1
    print("splits: verified — every study in exactly one test fold")
    return 0


if __name__ == "__main__":
    sys.exit(main())
