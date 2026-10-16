#!/usr/bin/env python3
"""Can paired-end structure break the artefact/biology G-run ambiguity?

fastprep trims the whole terminal G-run it observes. When the biology next to
the run is itself G's, those biological G's are trimmed too (docs/
OVERTRIM-DIAGNOSIS.md): a run of h artefact G's + k biological G's is
indistinguishable from h+k artefact G's using base identity alone.

This script tests whether the MATE breaks the tie:

  (a) mate overlap — a run position the other mate corroborates as G
      (complement C in the mate's orientation) is biology, not artefact;
      an artefact is a per-mate read-end phenomenon and is not corroborated
      at the same fragment position by the mate.
  (b) insert arithmetic — where the artefact run begins at the fragment
      boundary, h/t are exactly the run beyond the insert: without adapter
      readthrough h + t = L - I per mate, so knowing I plus one clean run
      boundary pins the other end.

Three analyses, one code path:

  probe   synthetic paired reads with KNOWN truth and forced ambiguous
          boundaries (biology edge G's merged into the artefact run), in
          five pad-free insert cohorts (adapter pads between biology and
          tail are excluded BY DESIGN: the mate-side evidence of a head
          artefact slot can land in the mate's adapter pad, whose first
          base complements to G and falsely corroborates — a pitfall the
          sim analysis measures separately):
            head_full    I = 150-h < L: merged head run (h artefact G's +
                         k biology G's), the whole fragment is in both
                         mates -> the rule applies at full strength
            head_partial L < I = 220 < 2L: same ambiguity, fragment edges
                         outside the mate's reach -> the rule abstains
            head_none    I = 380 >= 2L: no mate overlap at all
            tail_abuts   I = 138 = L-t: the 3' artefact begins exactly at
                         the fragment end (fastsim's fragment-end regime)
            tail_under   I = 220: the fragment continues under the 3'
                         artefact slots (fastsim's collapse regime)
          Runs fastprep, scores its cut against the truth, then applies the
          overlap-aware add-back rule and the insert rule and re-scores.

  sim     the four fastsim datasets in scratch (low/mid/high/mandril),
          scored exactly like bench/scripts/overtrim_sensitivity.py
          (PXMT adapter+polyX totals vs GT h/t). Rule = the mate-overlap
          add-back only (the insert rule needs the artefact to begin at the
          fragment boundary, which fastsim's model does not guarantee).

  bench   throughput of fastprep on 100k synthetic pairs and the measured
          cost of the rule per pair.

The add-back rule (the primary rule):

  For each read end, walk the observed terminal G-run from its INTERIOR edge
  outward and count how many run bases the mate corroborates as G
  (mate base == C at the corresponding fragment position). Those are
  biological G's fastprep would otherwise trim; add them back:

      corrected_cut = fastprep_cut - min(add_back, fastprep_cut)

  The walk stops at the first non-corroborated base. The rule never trims
  MORE than fastprep did, so it can only move overtrim -> exact (or ->
  undertrim where the corroboration is a coincidence; measured below).

Mate correspondence: R1[j] corresponds to rc(R2)[j - d]. The script finds d
by exhaustive shift search with mismatch counting (the same identity
fastcheck/src/insert.rs uses to measure the insert), so no oracle
information is used for the correspondence itself.

Insert measurement: `measure_insert` is a faithful re-implementation of
fastcheck's `--insert-hist` per-pair identity  I = b1 + b2 - overlap  with
b = L - (terminal homopolymer runs >= 3) and overlap mismatch tolerance 3.
fastcheck therefore already measures I on the pairs where the mates
overlap; no new machinery is required.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import os
import random
import struct
import subprocess
import sys
import time

import numpy as np

FASTPREP = "/home/raphael/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"
SIM_ROOT = "/home/raphael/.hermes/cache/scratch/overtrim"
WORK = "/home/raphael/.hermes/cache/scratch/overtrim_paired"
ADAPTER = b"CTGTCTCTTATACACATCT"
COMP_TAB = bytes.maketrans(b"ACGTN", b"TGCAN")

# base codes: A C G T N
CODE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    CODE[_b] = _i
COMP_CODE = np.array([3, 2, 1, 0, 4], dtype=np.uint8)  # A<->T, C<->G, N


# ---------------------------------------------------------------- PXMT parse
def parse_pxmt(path: str) -> dict[tuple[int, int], tuple]:
    """-> {(pair_id, mate): (status, t5a, t3a, t5p, t3p)}; mate 0-based."""
    out: dict[tuple[int, int], tuple] = {}
    buf = open(path, "rb").read()
    if buf[:4] != b"PXMT":
        raise ValueError(f"{path}: not PXMT")
    off, n = 5, len(buf)
    while off < n:
        (id_len,) = struct.unpack_from("<H", buf, off)
        off += 2
        rid = buf[off:off + id_len]
        off += id_len
        mate = buf[off]
        off += 1
        status = buf[off]
        off += 1
        t5a, t3a, t5p, t3p, _rl = struct.unpack_from("<IIIII", buf, off)
        off += 20
        n_runs = buf[off]
        off += 1
        off += 8 * n_runs
        try:
            parts = rid.split(b" ", 1)[0].decode("ascii").lstrip("@").split(":")
            x, y = int(parts[-2]), int(parts[-1])
            pid = y * 10_000 + x
        except (ValueError, UnicodeDecodeError, IndexError):
            pid = -1
        out[(pid, mate)] = (status, t5a, t3a, t5p, t3p)
    return out


def classify(t5: int, t3: int, h: int, t: int, tol: int) -> str:
    if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
        return "exact"
    if t5 > h or t3 > t:
        return "overtrim"
    return "undertrim"


# ------------------------------------------------------- mate correspondence
def encode(seq: bytes) -> np.ndarray:
    return CODE[np.frombuffer(seq, dtype=np.uint8)]


def find_delta(c1: np.ndarray, c2: np.ndarray, min_span: int = 24,
               min_rate: float = 0.75) -> int | None:
    """Find the shift d such that mate1[j] pairs with rc(mate2)[j - d].

    Exhaustive shift search counting mate-consistent positions (A<->T,
    C<->G), the same identity fastcheck uses. N on either side is ignored.
    Returns None when no candidate is confident (e.g. no overlap).
    """
    n = len(c1)
    rc2 = COMP_CODE[c2][::-1]
    best_d, best_score, best_span = None, -1, 0
    for d in range(-n + 1, n):
        lo = max(0, d)
        hi = min(n, d + n)
        if hi - lo < min_span:
            continue
        a = c1[lo:hi]
        b = rc2[lo - d:hi - d]
        valid = (a != 4) & (b != 4)
        vspan = int(valid.sum())
        if vspan < min_span:
            continue
        score = int((a[valid] == b[valid]).sum())
        if score > best_score:
            best_d, best_score, best_span = d, score, vspan
    if best_d is None or best_span == 0:
        return None
    if best_score / best_span < min_rate:
        return None
    return best_d


# ------------------------------------------- fastcheck's insert measurement
def _terminal_run(seq: bytes, from_start: bool) -> int:
    first = seq[0] if from_start else seq[-1]
    if first not in b"ACGT":
        return 0
    i, n = 0, len(seq)
    while i < n and (seq[i] if from_start else seq[n - 1 - i]) == first:
        i += 1
    return i


def measure_insert(s1: bytes, s2: bytes, anchor: int = 12,
                   max_mm: int = 3) -> int | None:
    """fastcheck --insert-hist per-pair identity: I = b1 + b2 - overlap.

    Terminal homopolymer runs >= 3 are excluded from the biology lengths
    (MIN_HOMOPOLYMER in fastcheck/src/insert.rs); the overlap matches R1's
    biology suffix against revcomp(R2)'s biology prefix with up to 3
    mismatches. Returns None when no overlap >= anchor verifies (the
    censored-above regime).
    """
    l1, l2 = len(s1), len(s2)
    tr1 = _terminal_run(s1, False)
    tr1 = tr1 if tr1 >= 3 else 0
    lr1 = _terminal_run(s1, True)
    lr1 = lr1 if lr1 >= 3 else 0
    tr2 = _terminal_run(s2, False)
    tr2 = tr2 if tr2 >= 3 else 0
    lr2 = _terminal_run(s2, True)
    lr2 = lr2 if lr2 >= 3 else 0
    b1, b2 = l1 - lr1 - tr1, l2 - lr2 - tr2
    if b1 <= 0 or b2 <= 0:
        return None
    rc2 = s2[::-1].translate(COMP_TAB)
    for o in range(min(b1, b2), anchor - 1, -1):
        seg1 = s1[l1 - tr1 - o:l1 - tr1]
        seg2 = rc2[tr2:tr2 + o]
        if len(seg1) != o or len(seg2) != o:
            continue
        mm = 0
        ok = True
        for x, y in zip(seg1, seg2):
            if x != y and x != 78 and y != 78:  # 'N'
                mm += 1
                if mm > max_mm:
                    ok = False
                    break
        if ok:
            return b1 + b2 - o
    return None


# ------------------------------------------------------- the add-back rule
def terminal_runs(seq: bytes) -> tuple[int, int]:
    n = len(seq)
    r5 = 0
    while r5 < n and seq[r5] == 71:  # 'G'
        r5 += 1
    r3 = 0
    while r3 < n and seq[n - 1 - r3] == 71:
        r3 += 1
    return r5, r3


def addback_counts(c1: np.ndarray, c2: np.ndarray, delta: int) -> dict:
    """Per-mate, per-end corroboration counts (the add-back rule).

    For each terminal G-run, walk from the run's interior edge outward and
    count run bases whose corresponding mate base is C (the fragment base
    is G -> biology). Stop at the first base without that evidence.
    """
    n = len(c1)
    rc1 = COMP_CODE[c1][::-1]
    rc2 = COMP_CODE[c2][::-1]
    out = {}
    for which, (own, mate_rc) in enumerate(((c1, rc2), (c2, rc1))):
        n_ = n
        r5 = 0
        while r5 < n_ and own[r5] == 2:
            r5 += 1
        r3 = 0
        while r3 < n_ and own[n_ - 1 - r3] == 2:
            r3 += 1
        # 5' run [0, r5): walk DOWN from interior edge r5-1 toward 0
        a5 = 0
        for j in range(r5 - 1, -1, -1):
            k = j - delta
            if 0 <= k < n and mate_rc[k] == 2:
                a5 += 1
            else:
                break
        # 3' run [n-r3, n): walk UP from interior edge n-r3 toward n-1
        a3 = 0
        for j in range(n - r3, n):
            k = j - delta
            if 0 <= k < n and mate_rc[k] == 2:
                a3 += 1
            else:
                break
        out[which] = (a5, a3)
    return out


# ---------------------------------------------------------- probe generator
def build_pair(rng: random.Random, pid: int, cohort: str, L: int = 150):
    """One synthetic pair with an explicit fragment and forced boundaries.

    Layout follows fastsim's build_read: [h Gs][biology][adapter?][t Gs].
    Each biology edge adjacent to an artefact run is FORCED to be a G- or
    C-run of length 0..2 with a non-run guard base, so the merge (the
    ambiguous case) happens or not by construction, exactly the case
    /tmp/probe_overtrim.py established: the error appears only when the
    boundary base is itself a G.
    """
    I = {"head_full": None, "head_partial": 220, "head_none": 380,
         "tail_abuts": 138, "tail_under": 220}[cohort]
    if cohort == "head_full":
        # merged head run, full mate coverage: symmetric h so that
        # I = L - h leaves no adapter pad and both mates see the fragment.
        h = rng.choice([0, 1, 2])
        h1 = h2 = h
        t1 = t2 = 0
        I = L - h
    elif cohort in ("head_partial", "head_none"):
        h1, h2 = rng.choice([0, 1, 2]), rng.choice([0, 1, 2])
        t1 = t2 = 0
    else:
        # tail cohorts: the 3' artefact sits at the fragment boundary
        # (tail_abuts, I = L - t) or past the read window (tail_under).
        h1 = h2 = 0
        t1 = t2 = 12
    b1, b2 = min(I, L - h1 - t1), min(I, L - h2 - t2)
    assert L - h1 - t1 - b1 == 0 and L - h2 - t2 - b2 == 0, "probe is pad-free"

    frag = bytearray(rng.choice(b"ACGT") for _ in range(I))

    def force(start: int, base: int, k: int, guard: int) -> None:
        base = int(base)
        for q in range(max(0, start), min(start + k, I)):
            frag[q] = base
        if 0 <= guard < I:
            frag[guard] = int(rng.choice([b for b in b"ACGT" if b != base]))

    # fragment 5' edge [0, nA): G merges with R1's head artefact; C merges
    # with R2's tail artefact (complement) where R2 reaches position 0.
    nA = rng.choice([0, 1, 2])
    force(0, rng.choice([71, 67]), nA, nA)                  # guard at nA
    # fragment 3' edge [I-nB, I): C merges with R2's head artefact; G
    # merges with R1's tail artefact where R1 reaches I-1.
    nB = rng.choice([0, 1, 2])
    force(I - nB, rng.choice([71, 67]), nB, I - nB - 1)     # guard interior
    # interior edges for the collapse cohorts (visible biology edge inside
    # the fragment): G merges with R1's tail; C merges with R2's tail.
    if b1 < I:
        nC1 = rng.choice([0, 1, 2])
        force(b1 - nC1, rng.choice([71, 67]), nC1, b1 - nC1 - 1)
    if b2 < I:
        nC2 = rng.choice([0, 1, 2])
        force(I - b2, rng.choice([71, 67]), nC2, I - b2 + nC2)

    bio1 = bytes(frag[:b1])
    bio2 = bytes(frag[I - b2:]).translate(COMP_TAB)[::-1]

    def mate_read(h: int, t: int, bio: bytes) -> bytes:
        window = L - h - t
        b = min(len(bio), window)
        pad = window - b
        seq = bytearray(b"G" * h)
        seq += bio[:b]
        seq += (ADAPTER * ((pad // len(ADAPTER)) + 1))[:pad]
        seq += b"G" * t
        assert len(seq) == L
        return bytes(seq)

    r1 = mate_read(h1, t1, bio1)
    r2 = mate_read(h2, t2, bio2)
    truth: dict = dict(pid=pid, cohort=cohort, I=I, h1=h1, t1=t1, h2=h2, t2=t2,
                 a1=L - h1 - t1 - b1, a2=L - h2 - t2 - b2)
    # per-mate artefact geometry: 'adapter' (adapter between biology and
    # tail), 'abuts' (the fragment ends where the artefact begins), or
    # 'under' (the fragment continues under the artefact slots).
    for m, (b, a) in ((1, (b1, L - h1 - t1 - b1)), (2, (b2, L - h2 - t2 - b2))):
        truth[f"geom{m}"] = ("adapter" if a > 0 else
                             "abuts" if I <= b else "under")
    # observed run lengths (= artefact + merged biological G's)
    r51, r31 = terminal_runs(r1)
    r52, r32 = terminal_runs(r2)
    truth["runs"] = (r51, r31, r52, r32)
    return r1, r2, truth


def write_fastq_pair(path1: str, path2: str, pairs) -> None:
    with open(path1, "wb") as f1, open(path2, "wb") as f2:
        for pid, r1, r2 in pairs:
            x, y = pid % 10_000, pid // 10_000
            hdr = f"@LH00916:1:SIMFLOW1:1:1101:{x}:{y}"
            q = b"I" * len(r1)
            f1.write(f"{hdr} 1:N:0:ATCGATCG\n".encode() + r1 + b"\n+\n" + q + b"\n")
            f2.write(f"{hdr} 2:N:0:ATCGATCG\n".encode() + r2 + b"\n+\n" + q + b"\n")


def run_fastprep(in1: str, in2: str, work: str, tag: str,
                 extra: list[str]) -> dict:
    out1 = f"{work}/{tag}_R1.fastq.gz"
    out2 = f"{work}/{tag}_R2.fastq.gz"
    pxmt = f"{work}/{tag}.pxmt"
    argv = [FASTPREP, "-1", in1, "-2", in2, "-o", out1, "-p", out2,
            "--metadata", pxmt, "--force", "-q", "-j", "4"] + extra
    t0 = time.time()
    p = subprocess.run(argv, capture_output=True, text=True)
    secs = time.time() - t0
    if p.returncode != 0:
        raise SystemExit(f"fastprep rc={p.returncode}: {p.stderr[-400:]}")
    return dict(pxmt=parse_pxmt(pxmt), secs=secs)


# ---------------------------------------------------------------- scoring
TOLS = (0, 2)


def new_res() -> dict:
    return {"before": {t: collections.Counter() for t in TOLS},
            "after": {t: collections.Counter() for t in TOLS},
            "fired": 0, "fixed": 0, "broken": 0, "fixed_sites": collections.Counter()}


def rate_block(counter: collections.Counter) -> str:
    n = max(1, counter["exact"] + counter["overtrim"] + counter["undertrim"])
    return (f"exact {counter['exact'] / n * 100:.1f}%  "
            f"overtrim {counter['overtrim'] / n * 100:.1f}%  "
            f"undertrim {counter['undertrim'] / n * 100:.1f}%  "
            f"(scored {n}, dropped {counter['dropped']})")


def insert_rule_cuts(r1: bytes, r2: bytes, i_est: int | None):
    """Rule (b): insert arithmetic in the fragment-end geometry.

    Without adapter readthrough h + t = L - I per mate. This ALONE does not
    split the run: every h in [max(0, budget-r3), min(r5, budget)] is
    self-consistent with some fragment. It resolves the boundary only when
    one end carries NO G-run at all (then that end's artefact is 0 and the
    other end is pinned), or when there is no artefact at all. Returns
    ((c5, c3), (c5, c3)) or None when the rule does not apply.
    """
    if i_est is None:
        return None
    out = []
    for seq in (r1, r2):
        L = len(seq)
        r5, r3 = terminal_runs(seq)
        budget = L - i_est
        if budget < 0 or budget > r5 + r3:
            return None  # not the fragment+artefact geometry
        if budget == 0:
            out.append((0, 0))
        elif r3 == 0 and budget <= r5:
            out.append((budget, 0))
        elif r5 == 0 and budget <= r3:
            out.append((0, budget))
        else:
            return None
    return tuple(out)


# ---------------------------------------------------------------- probe cmd
def cmd_probe(args) -> dict:
    os.makedirs(WORK, exist_ok=True)
    rng = random.Random(args.seed)
    cohorts = ["head_full", "head_partial", "head_none", "tail_abuts",
               "tail_under"]
    truth_rows: list[dict] = []
    pairs = []
    for cohort in cohorts:
        for _ in range(args.pairs_per_cohort):
            pid = len(pairs)
            r1, r2, truth = build_pair(rng, pid, cohort)
            truth["_r1"], truth["_r2"] = r1, r2
            truth_rows.append(truth)
            pairs.append((pid, r1, r2))
    in1, in2 = f"{WORK}/probe_R1.fastq", f"{WORK}/probe_R2.fastq"
    write_fastq_pair(in1, in2, pairs)
    run = run_fastprep(in1, in2, WORK, "probe",
                       ["--no-dedup", "--min-len", "80"])
    px = run["pxmt"]

    for trow in truth_rows:
        c1, c2 = encode(trow["_r1"]), encode(trow["_r2"])
        d = find_delta(c1, c2)
        trow["_ab"] = addback_counts(c1, c2, d) if d is not None else None
        trow["_delta_ok"] = d is not None
        i_est = measure_insert(trow["_r1"], trow["_r2"])
        trow["_i_est"] = i_est
        trow["_ins"] = insert_rule_cuts(trow["_r1"], trow["_r2"], i_est)

    def score(rows, ab, ins, geom=None) -> dict:
        """ab/ins select which corrections to apply; geom optionally filters
        mates by their artefact geometry tag."""
        r = new_res()
        for trow in rows:
            for mate in (0, 1):
                if geom is not None and trow[f"geom{mate + 1}"] != geom:
                    continue
                h, t = trow[f"h{mate + 1}"], trow[f"t{mate + 1}"]
                exp3 = t  # trim_accuracy semantics: the correct cut is
                # exactly h at 5' and t at 3'; an untrimmed adapter is
                # neither required nor penalised (matches the sim scorer)
                s = px.get((trow["pid"], mate))
                if s is None or s[0] != 0:
                    for tol in TOLS:
                        r["before"][tol]["dropped"] += 1
                        r["after"][tol]["dropped"] += 1
                    continue
                _, t5a, t3a, t5p, t3p = s
                t5, t3 = t5a + t5p, t3a + t3p
                c5, c3 = t5, t3
                if ab and trow["_ab"] is not None:
                    a5, a3 = trow["_ab"][mate]
                    if a5 or a3:
                        r["fired"] += 1
                    c5 = max(0, t5 - min(a5, t5))
                    c3 = max(0, t3 - min(a3, t3))
                elif ins and trow["_ins"] is not None:
                    c5, c3 = trow["_ins"][mate]
                    if (c5, c3) != (t5, t3):
                        r["fired"] += 1
                for tol in TOLS:
                    cb = classify(t5, t3, h, exp3, tol)
                    ca = classify(c5, c3, h, exp3, tol)
                    r["before"][tol][cb] += 1
                    r["after"][tol][ca] += 1
                    if tol == 0:
                        if cb == "overtrim" and ca == "exact":
                            r["fixed"] += 1
                            r["fixed_sites"]["5p" if t5 > h else "3p"] += 1
                        if cb != "undertrim" and ca == "undertrim":
                            r["broken"] += 1
        return r

    out = {"secs_fastprep": run["secs"], "n_pairs": len(truth_rows)}
    print(f"[probe] {len(truth_rows)} pairs, fastprep {run['secs']:.2f}s")
    cohorts = ["head_full", "head_partial", "head_none", "tail_abuts",
               "tail_under"]
    for cohort in cohorts:
        rows = [t for t in truth_rows if t["cohort"] == cohort]
        n = len(rows)
        meas = sum(1 for t in rows if t["_i_est"] is not None)
        exact_i = sum(1 for t in rows if t["_i_est"] == t["I"])
        within2 = sum(1 for t in rows
                      if t["_i_est"] is not None and abs(t["_i_est"] - t["I"]) <= 2)
        b_none = score(rows, ab=False, ins=False)
        b_ab = score(rows, ab=True, ins=False)
        b_ins = score(rows, ab=False, ins=True)
        print(f"\n== cohort {cohort}  (I={rows[0]['I']}, n={n})")
        print(f"  correspondence found: {sum(1 for t in rows if t['_delta_ok'])}/{n}"
              f"   fastcheck-style I measured: {meas}/{n}, exact {exact_i}/{n},"
              f" within±2 {within2}/{n}")
        print(f"  fastprep alone            tol0 {rate_block(b_none['before'][0])}")
        print(f"  {'':24s}  tol2 {rate_block(b_none['before'][2])}")
        print(f"  + mate-overlap add-back   tol0 {rate_block(b_ab['after'][0])}")
        print(f"  {'':24s}  tol2 {rate_block(b_ab['after'][2])}")
        print(f"    fired on {b_ab['fired']} mates, fixed {b_ab['fixed']} overtrims "
              f"({dict(b_ab['fixed_sites'])}), broke {b_ab['broken']} into undertrim")
        print(f"  + insert arithmetic       tol0 {rate_block(b_ins['after'][0])}")
        print(f"  {'':24s}  tol2 {rate_block(b_ins['after'][2])}")
        print(f"    fired on {b_ins['fired']} mates, fixed {b_ins['fixed']} overtrims, "
              f"broke {b_ins['broken']} into undertrim")
        out[cohort] = {
            "before_tol0": rate_block(b_none["before"][0]),
            "before_tol2": rate_block(b_none["before"][2]),
            "addback_tol0": rate_block(b_ab["after"][0]),
            "addback_tol2": rate_block(b_ab["after"][2]),
            "addback_fired": b_ab["fired"], "addback_fixed": b_ab["fixed"],
            "addback_broken": b_ab["broken"],
            "insert_tol0": rate_block(b_ins["after"][0]),
            "insert_tol2": rate_block(b_ins["after"][2]),
            "insert_fired": b_ins["fired"], "insert_fixed": b_ins["fixed"],
            "insert_broken": b_ins["broken"],
            "i_measured": meas, "i_exact": exact_i, "i_within2": within2,
        }

    # pooled per-mate artefact geometry (the physical regime each mate is in)
    print("\n== pooled by artefact geometry (per mate)")
    out["by_geometry"] = {}
    for geom in ("abuts", "adapter", "under"):
        b_none = score(truth_rows, ab=False, ins=False, geom=geom)
        b_ab = score(truth_rows, ab=True, ins=False, geom=geom)
        b_ins = score(truth_rows, ab=False, ins=True, geom=geom)
        n_mates = (b_none["before"][0]["exact"] + b_none["before"][0]["overtrim"]
                   + b_none["before"][0]["undertrim"])
        print(f"  {geom:8s} n={n_mates}")
        print(f"    fastprep alone        tol0 {rate_block(b_none['before'][0])}")
        print(f"    + mate-overlap        tol0 {rate_block(b_ab['after'][0])}"
              f"   (fixed {b_ab['fixed']}, broke {b_ab['broken']})")
        print(f"    + insert arithmetic   tol0 {rate_block(b_ins['after'][0])}"
              f"   (fixed {b_ins['fixed']}, broke {b_ins['broken']})")
        out["by_geometry"][geom] = {
            "n": n_mates,
            "before_tol0": rate_block(b_none["before"][0]),
            "addback_tol0": rate_block(b_ab["after"][0]),
            "addback_fixed": b_ab["fixed"], "addback_broken": b_ab["broken"],
            "insert_tol0": rate_block(b_ins["after"][0]),
            "insert_fixed": b_ins["fixed"], "insert_broken": b_ins["broken"],
        }
    return out


# ---------------------------------------------------------------- sim data
def read_fastq_seqs(path: str) -> dict[int, bytes]:
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rb") as fh:
        i, key = 0, None
        for line in fh:
            if i % 4 == 0:
                c = line[1:].split()[0].decode()
                x, y = int(c.split(":")[-2]), int(c.split(":")[-1])
                key = y * 10_000 + x
            elif i % 4 == 1:
                out[key] = line.rstrip()
            i += 1
    return out


def read_gt(path: str) -> dict[int, dict]:
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = dict(
                h1=int(r["r1_h"]), t1=int(r["r1_t"]),
                h2=int(r["r2_h"]), t2=int(r["r2_t"]),
                I=int(r["insert"]),
                a1=int(r["r1_adapter"]), a2=int(r["r2_adapter"]))
    return out


def cmd_sim(args) -> dict:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    pooled = new_res()
    for regime in args.regimes:
        work = os.path.join(args.sim_root, regime)
        g = read_gt(f"{work}/sim.v0.gt.tsv")
        seqs1 = read_fastq_seqs(f"{work}/sim.v0_R1.fastq.gz")
        seqs2 = read_fastq_seqs(f"{work}/sim.v0_R2.fastq.gz")
        run = run_fastprep(f"{work}/sim.v0_R1.fastq.gz",
                           f"{work}/sim.v0_R2.fastq.gz",
                           WORK, f"sim_{regime}", ["--no-dedup"])
        px = run["pxmt"]
        r = new_res()
        n_corr = 0
        ov_events_i = collections.Counter()   # overtrim events by insert class
        pad_split = collections.Counter()     # rule effects, pad vs no-pad mates
        for pid, gt in g.items():
            c1 = c2 = d = ab = None
            s1, s2 = seqs1.get(pid), seqs2.get(pid)
            if s1 is not None and s2 is not None:
                c1, c2 = encode(s1), encode(s2)
                d = find_delta(c1, c2)
                if d is not None:
                    n_corr += 1
                    ab = addback_counts(c1, c2, d)
            for mate in (0, 1):
                h, t = (gt["h1"], gt["t1"]) if mate == 0 else (gt["h2"], gt["t2"])
                s = px.get((pid, mate))
                if s is None or s[0] != 0:
                    for tol in TOLS:
                        r["before"][tol]["dropped"] += 1
                        r["after"][tol]["dropped"] += 1
                    continue
                _, t5a, t3a, t5p, t3p = s
                t5, t3 = t5a + t5p, t3a + t3p
                c5, c3 = t5, t3
                if ab is not None:
                    a5, a3 = ab[mate]
                    if a5 or a3:
                        r["fired"] += 1
                    c5 = max(0, t5 - min(a5, t5))
                    c3 = max(0, t3 - min(a3, t3))
                for tol in TOLS:
                    cb = classify(t5, t3, h, t, tol)
                    ca = classify(c5, c3, h, t, tol)
                    r["before"][tol][cb] += 1
                    r["after"][tol][ca] += 1
                    if tol == 0:
                        pad = (gt["a1"] if mate == 0 else gt["a2"]) > 0
                        key = "pad" if pad else "nopad"
                        if cb == "overtrim" and ca == "exact":
                            r["fixed"] += 1
                            r["fixed_sites"]["5p" if t5 > h else "3p"] += 1
                            pad_split[f"fixed_{key}"] += 1
                        if cb != "undertrim" and ca == "undertrim":
                            r["broken"] += 1
                            pad_split[f"broken_{key}"] += 1
                        if cb == "overtrim":
                            cls = ("I<L" if gt["I"] < 150 else
                                   "I<2L" if gt["I"] < 300 else "I>=2L")
                            ov_events_i[cls] += 1
        print(f"\n== sim {regime} ({len(g)} pairs, fastprep {run['secs']:.2f}s)")
        print(f"  correspondence found on {n_corr}/{len(g)} pairs")
        print(f"  fastprep alone          tol0 {rate_block(r['before'][0])}")
        print(f"  {'':22s}  tol2 {rate_block(r['before'][2])}")
        print(f"  fastprep + add-back     tol0 {rate_block(r['after'][0])}")
        print(f"  {'':22s}  tol2 {rate_block(r['after'][2])}")
        print(f"  rule fired on {r['fired']} mates; fixed {r['fixed']} overtrims "
              f"({dict(r['fixed_sites'])}); broke {r['broken']} into undertrim")
        print(f"  pad split: {dict(pad_split)}")
        print(f"  overtrim events by insert class: {dict(ov_events_i)}")
        out[regime] = {
            "before_tol0": rate_block(r["before"][0]),
            "before_tol2": rate_block(r["before"][2]),
            "after_tol0": rate_block(r["after"][0]),
            "after_tol2": rate_block(r["after"][2]),
            "fired": r["fired"], "fixed": r["fixed"], "broken": r["broken"],
            "corr": n_corr, "pairs": len(g),
            "overtrim_by_insert": dict(ov_events_i),
            "pad_split": dict(pad_split),
        }
        for tol in TOLS:
            pooled["before"][tol].update(r["before"][tol])
            pooled["after"][tol].update(r["after"][tol])
        pooled["fired"] += r["fired"]
        pooled["fixed"] += r["fixed"]
        pooled["broken"] += r["broken"]
        pooled["fixed_sites"].update(r["fixed_sites"])
    print(f"\n== pooled over {args.regimes}")
    print(f"  fastprep alone       tol0 {rate_block(pooled['before'][0])}")
    print(f"  {'':18s}  tol2 {rate_block(pooled['before'][2])}")
    print(f"  fastprep + add-back  tol0 {rate_block(pooled['after'][0])}")
    print(f"  {'':18s}  tol2 {rate_block(pooled['after'][2])}")
    print(f"  fired {pooled['fired']}, fixed {pooled['fixed']} "
          f"({dict(pooled['fixed_sites'])}), broke {pooled['broken']}")
    out["pooled"] = {
        "before_tol0": rate_block(pooled["before"][0]),
        "before_tol2": rate_block(pooled["before"][2]),
        "after_tol0": rate_block(pooled["after"][0]),
        "after_tol2": rate_block(pooled["after"][2]),
        "fired": pooled["fired"], "fixed": pooled["fixed"],
        "broken": pooled["broken"],
    }
    return out


# ---------------------------------------------------------------- bench
def cmd_bench(args) -> dict:
    os.makedirs(WORK, exist_ok=True)
    rng = random.Random(args.seed)
    cohorts = ["head_full", "head_partial", "head_none", "tail_abuts",
               "tail_under"]
    pairs = []
    for i in range(args.pairs):
        r1, r2, _ = build_pair(rng, i, cohorts[i % len(cohorts)])
        pairs.append((i, r1, r2))
    in1, in2 = f"{WORK}/bench_R1.fastq", f"{WORK}/bench_R2.fastq"
    write_fastq_pair(in1, in2, pairs)
    print(f"[bench] {len(pairs)} synthetic pairs written")

    run = run_fastprep(in1, in2, WORK, "bench", ["--no-dedup", "--min-len", "120"])
    pairs_per_s = len(pairs) / run["secs"]
    print(f"  fastprep paired path: {run['secs']:.2f}s "
          f"({pairs_per_s:,.0f} pairs/s wall, 4 threads)")

    t0 = time.time()
    subprocess.run([FASTPREP, "-1", in1, "-o", f"{WORK}/bench_se.fastq.gz",
                    "--metadata", f"{WORK}/bench_se.pxmt", "--force", "-q",
                    "-j", "4", "--no-dedup", "--min-len", "120"],
                   capture_output=True, text=True)
    se_secs = time.time() - t0
    print(f"  fastprep single-mate path: {se_secs:.2f}s "
          f"({len(pairs) / se_secs:,.0f} reads/s wall)")

    seqs = [(p[1], p[2]) for p in pairs]
    t0 = time.time()
    n_ok = 0
    for r1, r2 in seqs:
        c1, c2 = encode(r1), encode(r2)
        d = find_delta(c1, c2)
        if d is not None:
            addback_counts(c1, c2, d)
            n_ok += 1
    rule_secs = time.time() - t0
    per_pair_us = rule_secs / len(seqs) * 1e6
    print(f"  add-back rule (Python, 1 thread): {rule_secs:.2f}s "
          f"({per_pair_us:.0f} us/pair, correspondence found on {n_ok} pairs)")
    return dict(pairs=len(pairs), fastprep_secs=round(run["secs"], 2),
                fastprep_pairs_per_s=round(pairs_per_s),
                single_secs=round(se_secs, 2),
                rule_secs=round(rule_secs, 2),
                rule_us_per_pair=round(per_pair_us))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Paired-end resolution of the artefact/biology G-run ambiguity")
    ap.add_argument("mode", choices=["probe", "sim", "bench", "all"],
                    nargs="?", default="all")
    ap.add_argument("--seed", type=int, default=20261016)
    ap.add_argument("--pairs-per-cohort", type=int, default=600)
    ap.add_argument("--pairs", type=int, default=100_000)
    ap.add_argument("--sim-root", default=SIM_ROOT)
    ap.add_argument("--regimes", nargs="+",
                    default=["low", "mid", "high", "mandril"])
    ap.add_argument("--out", default=None, help="write the JSON summary here")
    args = ap.parse_args()

    summary = {}
    if args.mode in ("probe", "all"):
        summary["probe"] = cmd_probe(args)
    if args.mode in ("sim", "all"):
        summary["sim"] = cmd_sim(args)
    if args.mode in ("bench", "all"):
        summary["bench"] = cmd_bench(args)
    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(summary, fh, indent=1)
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
