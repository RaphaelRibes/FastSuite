#!/usr/bin/env python3
"""Weighted trim scoring at a pre-registered over:under cost ratio.

`trim_accuracy.py` counts mates as exact / overtrim / undertrim and leaves
DROPPED reads in an unweighted `dropped_or_missing` bucket, so "drop everything"
and "trim everything" both look harmless there. This script turns the same
per-read decisions into ONE number under the cost model that was fixed in
`bench/weights/trim-weights-v1.json` *before* the sweep ran, so a rule cannot be
judged on the axis it happens to improve.

THE COST MODEL (v1; read from the weights file, never hardcoded here)
    w_over   = 1     trimming real biology away, or dropping a read that should
                     have been kept, costs 1 unit
    w_under  = r     leaving artefact behind, or keeping a read that should have
                     been dropped, costs r units; r is swept over
                     ladder_r = [1, 3, 5, 7, 10] and 7 is the pre-registered
                     anchor ("we prefer less data than wrong data")
    unaligned        has its own column and is NEVER weighted into the score

PER-MATE CLASSES (each mate gets exactly one; the raw counts are weight-free)
    exact          |t5-h| <= tol and |t3-t| <= tol            (trim_accuracy.classify)
    overtrim       t5 > h or t3 > t                     real biology removed
    undertrim      t5 < h or t3 < t                     artefact left behind
    drop_correct   tool dropped AND the GT says drop: the mate is entirely
                   artefact (r1_full/r2_full, or the pair's full_flag) OR its
                   ideal keep window [h, min(I, L-t)) is shorter than the
                   min-len floor (--min-len-frac, default 0.9L — fastprep's
                   shipped --min-len 0.9L)
    drop_wrong     tool dropped a mate whose window is >= the floor and which
                   is not entirely artefact
    keep_wrong     tool kept a mate the GT says is entirely artefact
    unaligned      a kept mate whose cuts cannot describe a window of the read
                   (t5 + t3 > recorded read length), or t5/t3 unknown. It is
                   reported on its own axis and its PAIR is removed from both
                   the numerator and the denominator, so it can never move the
                   score in either direction.

THE SCORE (per dataset and pooled; every row stamps weights_version)
    mixed_trim = a pair whose mates are split over/under. Counted once per pair
                 (its two mates are moved out of the over/under counts) at
                 w_under + w_over.
    numerator  = exact + drop_correct
                 - w_over  * (overtrim_rem  + drop_wrong + mixed_trim)
                 - w_under * (undertrim_rem + keep_wrong + mixed_trim)
    score      = numerator / pairs_scored

    overtrim_rem = overtrim  - mixed_trim      (see `mixed_trim` above)
    undertrim_rem = undertrim - mixed_trim

    The divisor is `pairs` as the model specifies, while the numerator is a sum
    over MATES, so a perfect scorer reads 2.0 (1.0 per mate). Only the absolute
    scale is affected: every ratio, every ladder comparison and r* are invariant
    to it. `pairs_scored` is `pairs` minus the pairs holding an unaligned mate
    (equal to `pairs` on PXMT input, where unaligned cannot arise); both are
    emitted so the exclusion is auditable rather than silent.

BREAK-EVEN
    Priced per rule against a baseline on the SAME corpus:
        over_units  = overtrim  + drop_wrong       (weight-1 side, mixed mates
        under_units = undertrim + keep_wrong        already inside over/undertrim)
        r* = -delta_over / delta_under   (in reads)
    A variant is net-positive at r iff r > r* when it trades overtrim in for
    undertrim (delta_under < 0 < delta_over). When both deltas share a sign the
    rule dominates / is dominated at EVERY r > 0 and r* does not exist; those
    rows say `dominates` / `dominated` and leave r_star empty instead of
    printing a meaningless negative ratio. delta_score_at_anchor is the honest
    companion number: r* prices only the over/under trade, so a variant that
    also loses correct decisions can sit below the baseline at the anchor even
    when the trade itself is favourable.

SECONDARY, NEVER IN THE SCORE
    SECONDARY_bases_over  insert bases destroyed, summed over mates:
                          max(0, t5-h) + max(0, window_end - (L-t3))
    SECONDARY_bases_under artefact bases left: max(0, h-t5) + max(0, (L-t3) - window_end)
    Kept as their own clearly-prefixed columns; the primary score is per-read
    because a 3 bp slip and a 60 bp over-trim are not the same event.

INPUTS (both per-read; at least one is required)
    --pxmt PXMT --gt GT   fastprep metadata (PXMT v3, one record per read:
                          u16 id_len, id bytes, u8 mate (0-based), u8 status
                          (0 = kept), u32 t5_adapter, t3_adapter, t5_poly,
                          t3_poly, u32 read_len, u8 n_runs, n_runs*(u32,u32)).
                          Layout mirrors trim_accuracy.parse_pxmt, re-implemented
                          here because the score also needs read_len and the
                          window validity check. Cuts are t5a+t5p / t3a+t3p.
    --cuts TSV --gt GT    one row per decision: read_id, mate, status, t5, t3
                          [, read_len]. t5/t3 = `NA` means unaligned. This is the
                          fastp-style path (trim_accuracy infers its cuts by
                          aligning the trimmed read), and it is what makes the
                          unaligned axis testable without a corpus.

    A dataset is identified by its GT (so two --label values scored against the
    same GT are two rules on one corpus, which is what break-even compares).
    An aggregate-only mode over results/trim-accuracy.tsv is deliberately NOT
    offered: that table has one undifferentiated dropped_or_missing column, so
    drop_correct vs drop_wrong cannot be recovered from it and the score would
    be a guess.

Usage
    score_weighted.py --label fastprep-a --pxmt c12.pxmt --gt sim.v0.gt.tsv \\
                      --label fastprep-b --pxmt c31.pxmt --gt sim.v0.gt.tsv \\
                      --baseline fastprep-a --variant fastprep-b --out score.tsv
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import struct
import sys

WEIGHTS_DEFAULT = "bench/weights/trim-weights-v1.json"
MAGIC = b"PXMT"
CLASSES = ("exact", "overtrim", "undertrim", "drop_correct", "drop_wrong",
           "keep_wrong", "unaligned")

# row_type, weights_version, scope, dataset, label, baseline, variant, r,
# w_over, w_under, anchor_r, tol, min_len_frac, min_len_floor, read_len,
# pairs, pairs_scored, unaligned_pairs, mates_scored, mates_unscored,
# dup_pairs_skipped, exact, overtrim, undertrim, drop_correct, drop_wrong,
# keep_wrong, unaligned, mixed_trim, over_units, under_units,
# weighted_numerator, score, delta_over, delta_under, r_star,
# break_even_status, variant_better_at_anchor, delta_score_at_anchor,
# SECONDARY_bases_over, SECONDARY_bases_under
ROW_FIELDS = [
    "row_type", "weights_version", "scope", "dataset", "label", "baseline",
    "variant", "r", "w_over", "w_under", "anchor_r", "tol", "min_len_frac",
    "min_len_floor", "read_len", "pairs", "pairs_scored", "unaligned_pairs",
    "mates_scored", "mates_unscored", "dup_pairs_skipped",
    "clean_mates", "clean_false_trims", "clean_false_trim_rate",
] + list(CLASSES) + [
    "mixed_trim", "over_units", "under_units", "weighted_numerator", "score",
    "delta_over", "delta_under", "r_star", "break_even_status",
    "variant_better_at_anchor", "delta_score_at_anchor",
    "SECONDARY_bases_over", "SECONDARY_bases_under",
]


# ---------------------------------------------------------------- weights
def resolve_weights(path):
    """An explicit path wins; the default is resolved next to this script.

    `bench/weights/...` is the documented default, but the script must also work
    from any cwd, so the repo-relative form is tried as a fallback.
    """
    if os.path.exists(path):
        return path
    here = os.path.dirname(os.path.abspath(__file__))
    alt = os.path.join(here, os.pardir, os.pardir, path)
    alt = os.path.normpath(alt)
    if os.path.exists(alt):
        return alt
    raise SystemExit(f"weights file not found: {path} (also tried {alt})")


def load_weights(path):
    with open(path) as fh:
        w = json.load(fh)
    for key in ("w_over", "ladder_r"):
        if key not in w:
            raise SystemExit(f"{path}: missing {key!r}")
    # the pre-registered ladder must include the anchor, or "the anchor" is a
    # number nobody ever scored
    anchor = float(w.get("anchor_r", 7.0))
    rungs = [float(r) for r in w["ladder_r"]]
    if anchor not in rungs:
        rungs = sorted(rungs + [anchor])
    return {
        "path": path,
        "version": str(w.get("weights_version", w.get("version", "unknown"))),
        "w_over": float(w["w_over"]),
        "anchor_r": anchor,
        "ladder_r": rungs,
        "raw": w,
    }


# ---------------------------------------------------------------- readers
def read_key(rid, mate):
    """PXMT read id (Illumina coordinates) + 1-based mate -> (pair_id, mate).

    fastsim writes `@LH00916:<run>:<fc>:<lane>:<tile>:<x>:<y> <mate>:N:0:<idx>`
    with x = pair_id % 10000 and y = pair_id // 10000, so the pair id is exact
    with no occurrence counting and no reliance on file order — the same key
    trim_accuracy.header_key() uses.
    """
    fields = rid.split()[0].split(":")
    if len(fields) < 7:
        raise SystemExit(f"unexpected read id in PXMT: {rid!r}")
    return (int(fields[5]) + 10_000 * int(fields[6]), mate)


def parse_pxmt(path):
    """-> {(pair_id, mate): (status, t5, t3, read_len)}; status 0 = kept.

    Header: b"PXMT" + one version byte, then records of
      u16 id_len | id | u8 mate (0-BASED) | u8 status |
      u32 t5_adapter | u32 t3_adapter | u32 t5_poly | u32 t3_poly |
      u32 read_len | u8 n_runs | n_runs * (u32 start, u32 len)
    The per-run detail is irrelevant to a per-read score, so it is skipped by
    the byte count rather than parsed.
    """
    try:
        with open(path, "rb") as fh:
            buf = fh.read()
    except FileNotFoundError:
        raise SystemExit(f"PXMT not found: {path}")
    if buf[:4] != MAGIC:
        raise SystemExit(f"{path}: not PXMT (magic {buf[:4]!r})")
    out = {}
    off, n = 5, len(buf)
    while off + 2 <= n:
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len].decode("utf-8", "replace"); off += id_len
        if off + 22 > n:
            break
        # PXMT stores the mate 0-based (Mate::R1 => 0u8); normalise to 1-based
        # to match the GT's r1_/r2_ columns and the FASTQ suffix.
        mate = buf[off] + 1; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        (rlen,) = struct.unpack_from("<I", buf, off); off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out[read_key(rid, mate)] = (status, t5a + t5p, t3a + t3p, rlen)
    return out


def read_cuts(path):
    """-> {(pair_id, mate): (status, t5, t3, read_len)} from a decisions TSV.

    Columns: read_id, mate (1-based), status (0 = kept), t5, t3 [, read_len].
    `NA` / `*` / empty in t5 or t3 marks an unaligned decision (the fastp path:
    trim_accuracy's infer_cuts returns None when the trimmed read is not a pure
    substring of the original).
    """
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mate = int(r["mate"])
            pid = read_key(r["read_id"], mate)[0]
            status = int(r.get("status") or 0)
            t5 = _opt_int(r.get("t5"))
            t3 = _opt_int(r.get("t3"))
            rlen = _opt_int(r.get("read_len")) or 0
            out[(pid, mate)] = (status, t5, t3, rlen)
    return out


def _opt_int(v):
    if v is None:
        return None
    v = v.strip()
    if v == "" or v.upper() in ("NA", "*", "NAN", "NONE"):
        return None
    return int(v)


def read_gt(path):
    """-> {pair_id: {"insert": I, 1: (h, t, full), 2: (h, t, full), "is_dup": bool}}

    Schema (tools/fastsim/src/truth.rs): pair_id, insert, r1_h, r1_t, r2_h,
    r2_t, r1_adapter, r2_adapter, full_flag, is_dup [, r1_n, r2_n, r1_full,
    r2_full]. `full_flag` is PAIR-level (both mates entirely artefact); when the
    per-mate r1_full/r2_full columns are present they win, because a pair can
    have one all-artefact mate and one with biology and only the per-mate flags
    say which.
    """
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid = int(r["pair_id"])
            pair_full = int(r.get("full_flag", 0) or 0) == 1
            per = {}
            for mate, pre in ((1, "r1"), (2, "r2")):
                full = pair_full
                flag = r.get(f"{pre}_full")
                if flag not in (None, ""):
                    full = int(flag) == 1
                per[mate] = (int(r[f"{pre}_h"]), int(r[f"{pre}_t"]), full)
            out[pid] = {
                "insert": int(r["insert"]),
                1: per[1],
                2: per[2],
                "is_dup": int(r.get("is_dup", 0) or 0) == 1,
            }
    return out


# ---------------------------------------------------------------- classify
def classify_cuts(t5, t3, h, t, tol):
    """Byte-for-byte the semantics of trim_accuracy.classify()."""
    if t5 is None or t3 is None:
        return "unaligned"
    if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
        return "exact"
    if t5 > h or t3 > t:
        return "overtrim"
    return "undertrim"


def classify_mate(dec, row, mate, read_len, min_len_floor, tol):
    """-> (class, bases_over, bases_under) for one mate.

    Drop correctness needs the GT's full flag AND the min-len floor, which is
    exactly the pair of facts trim_accuracy's `dropped_or_missing` bucket throws
    away.
    """
    status, t5, t3, rlen = dec
    h, t, full = row[mate]
    rl = rlen or read_len
    if rl <= 0:
        raise SystemExit("read length unknown: pass --read-len (or a read_len "
                         "column in --cuts)")
    # the correct kept window is [h, min(I, L-t)); its length decides whether a
    # drop was the right call
    window_end = min(row["insert"], rl - t)
    window = max(0, window_end - h)

    if status != 0:
        if full or window < min_len_floor:
            return "drop_correct", 0, 0
        return "drop_wrong", 0, 0
    if full:
        # kept a mate the GT calls entirely artefact: no window to measure, and
        # scoring its boundaries would price the same mistake twice
        return "keep_wrong", 0, 0
    if t5 is None or t3 is None:
        return "unaligned", 0, 0
    if t5 + t3 > rl:
        # a kept read keeps L-t5-t3 bases; negative means the recorded cuts do
        # not describe this read at all (a length mismatch), so there is no
        # window to compare — unaligned, not an over-trim
        return "unaligned", 0, 0
    cls = classify_cuts(t5, t3, h, t, tol)
    if window_end <= h:
        return cls, 0, 0
    over = max(0, t5 - h) + max(0, window_end - (rl - t3))
    under = max(0, h - t5) + max(0, (rl - t3) - window_end)
    return cls, over, under


# ---------------------------------------------------------------- scoring
def buckets():
    c = collections.Counter({k: 0 for k in CLASSES})
    c["mixed_trim"] = 0
    return c


def score_dataset(decisions, gt, read_len, min_len_floor, tol, skip_dup=False):
    """-> per-dataset counts + the derived axes.

    Keys present in `decisions` but not in `gt` (or vice versa) are reported as
    `mates_unscored` rather than silently dropped or guessed: a missing key is a
    join failure, not a trim decision.
    """
    c = buckets()
    bases_over = bases_under = 0
    mates_unscored = 0
    dup_pairs_skipped = 0
    # the pre-registered false-trim ceiling is measured on CLEAN reads: mates
    # whose GT window is the whole read, i.e. nothing should have been cut
    clean_mates = clean_false_trims = 0
    by_pid = collections.defaultdict(dict)
    for (pid, mate), dec in decisions.items():
        by_pid[pid][mate] = dec

    pairs = 0
    unaligned_pairs = 0
    mixed = 0
    for pid in sorted(by_pid):
        row = gt.get(pid)
        if row is None:
            mates_unscored += len(by_pid[pid])
            continue
        if skip_dup and row["is_dup"]:
            # a PCR duplicate is a dedup decision, not a trim decision; excluded
            # from the trim score rather than scored as a wrong drop
            dup_pairs_skipped += 1
            continue
        pairs += 1
        classes = {}
        for mate in (1, 2):
            dec = by_pid[pid].get(mate)
            if dec is None:
                mates_unscored += 1
                continue
            cls, bo, bu = classify_mate(dec, row, mate, read_len,
                                        min_len_floor, tol)
            c[cls] += 1
            classes[mate] = cls
            bases_over += bo
            bases_under += bu
            h_gt, t_gt, full_gt = row[mate]
            if h_gt == 0 and t_gt == 0 and not full_gt:
                clean_mates += 1
                st, t5d, t3d, _ = dec
                # an unaligned classification is a join/length failure, not a
                # trim decision, so it never enters the false-trim ceiling
                if cls != "unaligned" and st == 0 and (t5d or 0) + (t3d or 0) > 0:
                    clean_false_trims += 1
        if "unaligned" in classes.values():
            unaligned_pairs += 1
        elif len(classes) == 2 and set(classes.values()) == {"overtrim", "undertrim"}:
            mixed += 1
    c["mixed_trim"] = mixed
    return {
        "counts": c,
        "pairs": pairs,
        "pairs_scored": pairs - unaligned_pairs,
        "unaligned_pairs": unaligned_pairs,
        "mates_scored": sum(c[k] for k in CLASSES),
        "mates_unscored": mates_unscored,
        "dup_pairs_skipped": dup_pairs_skipped,
        "bases_over": bases_over,
        "bases_under": bases_under,
        "clean_mates": clean_mates,
        "clean_false_trims": clean_false_trims,
    }


def over_units(counts):
    """Weight-1 side of the loss, in reads. Mixed mates already sit inside
    overtrim/undertrim, so they are not added again."""
    return counts["overtrim"] + counts["drop_wrong"]


def under_units(counts):
    return counts["undertrim"] + counts["keep_wrong"]


def weighted_numerator(counts, w_under, w_over):
    """exact + drop_correct - w_over*over - w_under*under.

    The mixed-pair term of the model, (w_under + w_over) * mixed_trim, is folded
    in by moving those mates out of the raw over/under counts first; the two
    forms are algebraically the same number, and the raw counts stay as
    trim_accuracy reports them.
    """
    m = counts["mixed_trim"]
    over_rem = counts["overtrim"] - m
    under_rem = counts["undertrim"] - m
    return (counts["exact"] + counts["drop_correct"]
            - w_over * (over_rem + counts["drop_wrong"] + m)
            - w_under * (under_rem + counts["keep_wrong"] + m))


def score_at(counts, r, w_over, pairs_scored):
    num = weighted_numerator(counts, r, w_over)
    return num, (num / pairs_scored if pairs_scored else 0.0)


def break_even(d_over, d_under, anchor):
    """-> (r_star|None, status, variant_better_at_anchor).

    score_v - score_b = (delta_positive - d_over - r*d_under) / pairs, so the
    variant is net-positive iff r*d_under + d_over < 0. Solving that against the
    anchor keeps the verdict and r* reported together rather than letting a
    reader infer the direction.
    """
    if d_under == 0:
        if d_over < 0:
            return None, "dominates", True      # better on the over side, tied under
        if d_over > 0:
            return None, "dominated", False
        return None, "identical", False
    r_star = -d_over / d_under
    if d_over < 0 and d_under < 0:
        return None, "dominates", True          # better on both axes
    if d_over > 0 and d_under > 0:
        return None, "dominated", False
    if d_under < 0:
        # trades overtrim in for undertrim: wins only above r*
        return r_star, "above", anchor > r_star
    return r_star, "below", anchor < r_star     # wins only below r*


# ---------------------------------------------------------------- output
def fmt(v):
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        return repr(v)
    return str(v)


def emit_rows(w, args_cfg, measurements, out_fh):
    """measurements: {(dataset, label): score_dataset() result}."""
    w_over = w["w_over"]
    rows = []

    def base_row(**kw):
        row = {k: "" for k in ROW_FIELDS}
        row.update({
            "weights_version": w["version"],
            "w_over": w_over,
            "anchor_r": w["anchor_r"],
            "tol": args_cfg["tol"],
            "min_len_frac": args_cfg["min_len_frac"],
            "min_len_floor": args_cfg["min_len_floor"],
            "read_len": args_cfg["read_len"],
        })
        row.update(kw)
        return row

    def score_rows(scope, dataset, label, m):
        for r in w["ladder_r"]:
            num, sc = score_at(m["counts"], r, w_over, m["pairs_scored"])
            rows.append(base_row(
                row_type="score", scope=scope, dataset=dataset, label=label,
                r=r, w_under=r,
                pairs=m["pairs"], pairs_scored=m["pairs_scored"],
                unaligned_pairs=m["unaligned_pairs"],
                mates_scored=m["mates_scored"],
                mates_unscored=m["mates_unscored"],
                dup_pairs_skipped=m["dup_pairs_skipped"],
                clean_mates=m["clean_mates"],
                clean_false_trims=m["clean_false_trims"],
                clean_false_trim_rate=(m["clean_false_trims"] / m["clean_mates"]
                                       if m["clean_mates"] else ""),
                mixed_trim=m["counts"]["mixed_trim"],
                over_units=over_units(m["counts"]),
                under_units=under_units(m["counts"]),
                weighted_numerator=num, score=sc,
                SECONDARY_bases_over=m["bases_over"],
                SECONDARY_bases_under=m["bases_under"],
                **{k: m["counts"][k] for k in CLASSES},
            ))

    for (ds, label) in sorted(measurements):
        score_rows("dataset", ds, label, measurements[(ds, label)])

    # pooled per label over every dataset it was measured on (counts summed, then
    # scored — not a mean of per-dataset scores, which would weight a 200k-read
    # corpus the same as a 1k one)
    labels = sorted({lab for (_ds, lab) in measurements})
    pooled = {}
    for label in labels:
        agg = buckets()
        acc = {"counts": agg, "pairs": 0, "pairs_scored": 0, "unaligned_pairs": 0,
               "mates_scored": 0, "mates_unscored": 0, "dup_pairs_skipped": 0,
               "bases_over": 0, "bases_under": 0, "clean_mates": 0,
               "clean_false_trims": 0}
        for (ds, lab), m in sorted(measurements.items()):
            if lab != label:
                continue
            for k in CLASSES + ("mixed_trim",):
                agg[k] += m["counts"][k]
            for k in ("pairs", "pairs_scored", "unaligned_pairs", "mates_scored",
                      "mates_unscored", "dup_pairs_skipped", "bases_over",
                      "bases_under", "clean_mates", "clean_false_trims"):
                acc[k] += m[k]
        pooled[label] = acc
        score_rows("AGGREGATE", "*", label, acc)

    # ---- break-even, per common dataset and pooled over them
    base_label = args_cfg["baseline"]
    for var in args_cfg["variants"]:
        common = sorted(ds for (ds, lab) in measurements if lab == base_label
                        and (ds, var) in measurements)
        for ds in common:
            b, v = measurements[(ds, base_label)], measurements[(ds, var)]
            rows.append(breakeven_row(base_row, ds, base_label, var, b, v, w))
        if common:
            b, v = pooled[base_label], pooled[var]
            rows.append(breakeven_row(base_row, "*", base_label, var, b, v, w))

    writer = csv.DictWriter(out_fh, fieldnames=ROW_FIELDS, delimiter="\t",
                            lineterminator="\n", extrasaction="raise")
    writer.writeheader()
    for row in rows:
        writer.writerow({k: fmt(row.get(k)) for k in ROW_FIELDS})
    return rows


def breakeven_row(base_row, dataset, base_label, var, b, v, w):
    d_over = over_units(v["counts"]) - over_units(b["counts"])
    d_under = under_units(v["counts"]) - under_units(b["counts"])
    r_star, status, better = break_even(d_over, d_under, w["anchor_r"])
    _, sb = score_at(b["counts"], w["anchor_r"], w["w_over"], b["pairs_scored"])
    _, sv = score_at(v["counts"], w["anchor_r"], w["w_over"], v["pairs_scored"])
    return base_row(
        row_type="break_even", scope="dataset" if dataset != "*" else "AGGREGATE",
        dataset=dataset, baseline=base_label,
        variant=var, delta_over=d_over, delta_under=d_under, r_star=r_star,
        break_even_status=status, variant_better_at_anchor=better,
        delta_score_at_anchor=sv - sb,
    )


# ---------------------------------------------------------------- main
def build_inputs(a):
    kinds = [("pxmt", p) for p in a.pxmt] + [("cuts", p) for p in a.cuts]
    if not kinds:
        raise SystemExit("no input: pass --pxmt PATH and/or --cuts PATH "
                         "(each with a matching --gt)")
    n = len(kinds)
    if len(a.gt) != n:
        raise SystemExit(f"{n} input(s) but {len(a.gt)} --gt value(s); a GT is "
                         f"required per input (the score needs full_flag and "
                         f"insert, which only the GT carries)")
    if a.label and len(a.label) != n:
        raise SystemExit(f"{n} input(s) but {len(a.label)} --label value(s)")
    if a.name and len(a.name) != n:
        raise SystemExit(f"{n} input(s) but {len(a.name)} --name value(s)")
    out = []
    for i, (kind, path) in enumerate(kinds):
        gt = a.gt[i]
        label = a.label[i] if a.label else os.path.splitext(
            os.path.basename(path))[0]
        name = a.name[i] if a.name else os.path.basename(
            os.path.dirname(os.path.abspath(gt))) or "dataset"
        out.append({"kind": kind, "path": path, "gt": gt, "label": label,
                    "dataset": name})
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Weighted trim score at r in the pre-registered ladder, "
                    "plus the break-even r* against a baseline rule.")
    ap.add_argument("--weights", default=WEIGHTS_DEFAULT,
                    help=f"cost model (default {WEIGHTS_DEFAULT})")
    ap.add_argument("--pxmt", action="append", default=[],
                    help="PXMT metadata for one (dataset, label); repeatable")
    ap.add_argument("--cuts", action="append", default=[],
                    help="decisions TSV for one (dataset, label); repeatable")
    ap.add_argument("--gt", action="append", default=[],
                    help="sim.vN.gt.tsv, one per input, same order")
    ap.add_argument("--label", action="append", default=[],
                    help="rule name per input (default: input file stem)")
    ap.add_argument("--name", action="append", default=[],
                    help="dataset name per input (default: GT's parent dir)")
    ap.add_argument("--baseline", default=None, help="label to price variants against")
    ap.add_argument("--variant", action="append", default=[],
                    help="label to compare with the baseline; repeatable")
    ap.add_argument("--min-len-frac", type=float, default=0.9,
                    help="min-len floor as a fraction of the read length "
                         "(fastprep ships 0.9L; default 0.9)")
    ap.add_argument("--tol", type=int, default=0,
                    help="boundary tolerance in bp for 'exact' (0 strict, 2 forgiving)")
    ap.add_argument("--read-len", type=int, default=0,
                    help="read length L; 0 = infer from PXMT read_len")
    ap.add_argument("--skip-dup-pairs", action="store_true",
                    help="exclude is_dup=1 pairs (dedup, not trim) from the score")
    ap.add_argument("--out", default="-", help="output TSV (default stdout)")
    a = ap.parse_args(argv)

    w = load_weights(resolve_weights(a.weights))
    inputs = build_inputs(a)
    if a.baseline and a.baseline not in {i["label"] for i in inputs}:
        raise SystemExit(f"--baseline {a.baseline!r} is not one of the labels: "
                         f"{sorted({i['label'] for i in inputs})}")
    for v in a.variant:
        if v not in {i["label"] for i in inputs}:
            raise SystemExit(f"--variant {v!r} is not one of the labels: "
                             f"{sorted({i['label'] for i in inputs})}")

    parsed = []
    for spec in inputs:
        dec = parse_pxmt(spec["path"]) if spec["kind"] == "pxmt" else read_cuts(spec["path"])
        parsed.append((spec, dec, read_gt(spec["gt"])))

    read_len = a.read_len
    if not read_len:
        lengths = collections.Counter()
        for _spec, dec, _gt in parsed:
            for (_k, (_s, _t5, _t3, rl)) in dec.items():
                if rl:
                    lengths[rl] += 1
        read_len = lengths.most_common(1)[0][0] if lengths else 0
    if not read_len:
        raise SystemExit("read length unknown: pass --read-len (the min-len "
                         "floor is a fraction of it)")
    floor = int(a.min_len_frac * read_len)

    measurements = {}
    for spec, dec, gt in parsed:
        key = (spec["dataset"], spec["label"])
        if key in measurements:
            raise SystemExit(f"duplicate measurement for dataset "
                             f"{spec['dataset']!r} label {spec['label']!r}")
        measurements[key] = score_dataset(dec, gt, read_len, floor, a.tol,
                                          a.skip_dup_pairs)

    cfg = {"tol": a.tol, "min_len_frac": a.min_len_frac, "min_len_floor": floor,
           "read_len": read_len, "baseline": a.baseline, "variants": a.variant}
    if a.out in ("", "-"):
        emit_rows(w, cfg, measurements, sys.stdout)
    else:
        d = os.path.dirname(os.path.abspath(a.out))
        os.makedirs(d, exist_ok=True)
        with open(a.out, "w", newline="") as fh:
            emit_rows(w, cfg, measurements, fh)
        print(f"wrote {a.out} (weights_version {w['version']}, "
              f"r ladder {w['ladder_r']}, anchor {w['anchor_r']}, "
              f"min_len_floor {floor} = {a.min_len_frac:g}L of {read_len})",
              file=sys.stderr)


if __name__ == "__main__":
    main()
