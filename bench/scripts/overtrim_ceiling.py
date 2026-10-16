#!/usr/bin/env python3
"""Bayes-optimal ceiling for run-length trimming decisions (fastprep overtrim).

QUESTION. fastprep trims a poly-G artefact by cutting the ENTIRE observed run,
i.e. its decision rule is f(r) = r where r is the observed run length at the
read end.  Measured overtrim (strict tolerance 0, vs fastsim ground truth) is
29.5%.  Is ANY rule of the observed run length alone able to fix that, or is
f(r) already optimal given base identity alone?

WHAT THIS SCRIPT COMPUTES (per regime, per end, and mate-joint):

  1. The empirical conditional P(true artefact length | observed run length r)
     from fastsim's exact ground truth (sim.v0.gt.tsv: pair_id insert r1_h
     r1_t r2_h r2_t r1_adapter r2_adapter full_flag is_dup).
  2. Every candidate rule f(r) -> cut — constant, f(r)=r, f(r)=r-c,
     f(r)=min(r,k), f(r)=round(a*r), posterior-thresholded mode — scored at
     tolerance 0 and 2 (exact / overtrim / undertrim, mean signed error).
  3. The Bayes-optimal rules (posterior mode at tol 0; best covered-mass cut at
     tol 2; posterior mean) and the CEILING  sum_r max_cut P(exact | r) P(r),
     which is the best strict exact rate ANY function of r can reach.
  4. The class flow (exact/over/under -> exact/over/under) between the current
     rule and the optimum, which answers whether overtrim can drop without an
     equal undertrim increase.

TWO ARMS.

  fastsim arm  — reads simulated by tools/fastsim from the synthetic profiles
                 in bench/fixtures/overtrim/ (see profiles_overtrim.py),
                 scored against the simulator's GT.  Authoritative arm.
  merge arm    — an ANALYTIC sensitivity model where the biological G-run is
                 allowed to continue the artefact run (r = h + b, b ~ Geom),
                 which fastsim's boundary rule forbids (artefact.rs: when the
                 artefact is > 0 the base next to it is forced non-G, so a
                 simulated run is exactly the artefact unless the artefact is
                 zero).  This is the mechanism /tmp/probe_overtrim.py
                 demonstrated on real fastprep behaviour.  It is NOT fastsim
                 ground truth; it is here to bound what the ceiling looks like
                 when merging is allowed.

An optional fastprep reference run scores the real binary on the same data and
prints its empirical cut policy E[cut | r], so f(r)=r can be checked against
what the tool actually does.

Deliverable writeup: docs/OVERTRIM-CEILING.md.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import struct
import subprocess
import sys
from collections import Counter, defaultdict
from typing import Any

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FASTSIM = os.path.join(ROOT, "tools", "fastsim", "target", "release", "fastsim")
FASTPREP = os.path.join(ROOT, "tools", "fastprep", "target", "release", "fastprep")
FIXTURES = os.path.join(ROOT, "bench", "fixtures", "overtrim")
PROFILES_PY = os.path.join(ROOT, "bench", "scripts", "profiles_overtrim.py")
WORK = os.path.join(os.path.expanduser("~"), ".hermes", "cache", "scratch", "overtrim_ceiling")

RMAX = 25          # largest observed run length bucket kept (rest -> "25+")
CMAX = 25          # largest candidate cut considered by the Bayes rules
TOLS = (0, 2)

# fastprep reference config: the shipped default (bench/scripts/sweep.py base).
FP_BASE = ["--min-run", "12", "--mercy", "3", "--phred-sat", "20",
           "--gate-shape", "1.0", "--min-len", "0.9L", "--mean-q", "25.0",
           "--adapter-error-rate", "0.2", "--min-adapter-overlap", "14",
           "--gz-level", "6", "--head-dropout", "4"]


# ----------------------------------------------------------------- fixtures
def ensure_fixtures() -> None:
    if os.path.exists(os.path.join(FIXTURES, "mid.json")):
        return
    subprocess.run([sys.executable, PROFILES_PY], check=True, cwd=ROOT)


# ----------------------------------------------------------------- simulate
def simulate(regime: str, pairs: int, seed: int) -> str:
    outdir = os.path.join(WORK, regime)
    os.makedirs(outdir, exist_ok=True)
    prefix = os.path.join(outdir, "sim.v0")
    gt = prefix + ".gt.tsv"
    if os.path.exists(gt):
        # reuse only if the cached dataset matches the request
        with open(gt) as fh:
            n_rows = sum(1 for _ in fh) - 1
        if n_rows == pairs:
            return outdir
    profile = os.path.join(FIXTURES, f"{regime}.json")
    argv = [FASTSIM, "--pairs", str(pairs), "--profile", profile,
            "-o", prefix, "--seed", str(seed), "-f", "-q"]
    p = subprocess.run(argv, capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f"fastsim {regime} failed: {p.stderr[-400:]}")
    return outdir


def load_gt(path: str) -> dict[int, tuple[int, int, int, int]]:
    out: dict[int, tuple[int, int, int, int]] = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        want = {"pair_id", "r1_h", "r1_t", "r2_h", "r2_t"}
        if not want.issubset(rd.fieldnames or []):
            raise SystemExit(f"unexpected GT columns: {rd.fieldnames}")
        for r in rd:
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out


# ----------------------------------------------------------------- scanning
def run_len(seq: bytes, end: int) -> int:
    """Observed G-run at the read end (end=5 leading, end=3 trailing)."""
    if end == 5:
        i = 0
        while i < len(seq) and seq[i] == 71:  # 'G'
            i += 1
        return i
    i = len(seq) - 1
    n = 0
    while i >= 0 and seq[i] == 71:
        i -= 1
        n += 1
    return n


def scan_mates(path: str, mate: int):
    """Yield (pair_id, r5, r3) for every record of one FASTQ(.gz)."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rb") as fh:
        while True:
            name = fh.readline()
            if not name:
                return
            seq = fh.readline().strip()
            fh.readline()
            fh.readline()
            bare = name[1:].split(None, 1)[0].decode("ascii")
            parts = bare.split(":")
            x, y = int(parts[-2]), int(parts[-1])
            yield y * 10_000 + x, run_len(seq, 5), run_len(seq, 3)


# ----------------------------------------------------------------- tables
class EndTable:
    """Weighted counts of (observed r, true artefact length) for one end."""

    def __init__(self, name: str):
        self.name = name
        self.joint: dict[int, Any] = defaultdict(Counter)  # r -> {h: w}

    def add(self, r: int, h: float, w: float = 1.0) -> None:
        self.joint[min(r, RMAX)][min(h, RMAX)] += w

    def post(self) -> dict[int, dict[int, float]]:
        """P(h | r) -> {r: {h: p}}."""
        out = {}
        for r, c in self.joint.items():
            tot = sum(c.values())
            out[r] = {h: w / tot for h, w in c.items()} if tot else {}
        return out

    def total(self) -> float:
        return sum(sum(c.values()) for c in self.joint.values())

    def p_r(self) -> dict[int, float]:
        tot = self.total()
        return {r: sum(c.values()) / tot for r, c in self.joint.items()}


class MateTable:
    """Weighted counts of ((r5, r3), (h, t)) over mates."""

    def __init__(self):
        self.cells: dict[tuple[int, int], Any] = defaultdict(Counter)

    def add(self, r5: int, r3: int, h: int, t: int, w: float = 1.0) -> None:
        self.cells[(min(r5, RMAX), min(r3, RMAX))][
            (min(h, RMAX), min(t, RMAX))] += w

    def total(self) -> float:
        return sum(sum(c.values()) for c in self.cells.values())


# ----------------------------------------------------------------- rules
def mode_cut(post: dict[int, float]) -> int:
    """argmax_h P(h|r); ties -> smallest h (the conservative cut)."""
    return min(post, key=lambda h: (-post[h], h))


def make_rules(tab: EndTable, tol: int) -> list[tuple[str, object]]:
    """Candidate f(r) -> cut.  Bayes rules derived from this end's table."""
    post = tab.post()
    rules: list[tuple[str, object]] = []

    rules.append(("cut=r (current)", lambda r: r))
    rules.append(("cut=0 (never)", lambda r: 0))
    for c in (1, 2, 3):
        rules.append((f"cut=max(0,r-{c})", lambda r, c=c: max(0, r - c)))
    for k in (1, 2, 3, 4, 6, 8, 12):
        rules.append((f"cut=min(r,{k})", lambda r, k=k: min(r, k)))
    for a in (0.25, 0.5, 0.75, 1.25, 1.5, 2.0):
        rules.append((f"cut=round({a}*r)",
                      lambda r, a=a: int(math.floor(a * r + 0.5))))
    for c in (1, 2, 3, 6):
        rules.append((f"cut=min({c},r) if r>0", lambda r, c=c: min(c, r) if r > 0 else 0))

    # Bayes-optimal cut at this tolerance: maximise P(|h - cut| <= tol | r).
    def bayes_tol(r: int, post=post, tol=tol) -> int:
        p = post.get(r, {})
        if not p:
            return 0
        best_c, best_p = 0, -1.0
        for c in range(0, CMAX + 1):
            mass = sum(w for h, w in p.items() if abs(h - c) <= tol)
            if mass > best_p + 1e-15:
                best_c, best_p = c, mass
        return best_c

    rules.append((f"bayes-cut@tol{tol}", bayes_tol))

    if tol == 0:
        rules.append(("bayes-mode", lambda r, post=post: mode_cut(post.get(r, {0: 1.0}))))
        rules.append(("bayes-mean",
                      lambda r, post=post: int(math.floor(
                          sum(h * w for h, w in post.get(r, {0: 1.0}).items()) + 0.5))))
        for th in (0.25, 0.5, 0.75, 0.9, 0.99):
            def thresh(r: int, post=post, th=th) -> int:
                p = post.get(r, {})
                if not p:
                    return 0
                c = mode_cut(p)
                return c if p[c] >= th else 0
            rules.append((f"mode if P>={th} else 0", thresh))
    return rules


# ----------------------------------------------------------------- scoring
def classify(cut5: int, cut3: int, h: int, t: int, tol: int) -> str:
    """Same convention as bench/scripts/overtrim_sensitivity.py."""
    if abs(cut5 - h) <= tol and abs(cut3 - t) <= tol:
        return "exact"
    if cut5 > h or cut3 > t:
        return "overtrim"
    return "undertrim"


def score_end(tab: EndTable, rule, tol: int) -> dict:
    c = Counter()
    err = 0.0
    for r, hc in tab.joint.items():
        cut = rule(r)
        for h, w in hc.items():
            d = cut - h
            if abs(d) <= tol:
                c["exact"] += w
            elif d > 0:
                c["overtrim"] += w
            else:
                c["undertrim"] += w
            err += d * w
    n = tab.total()
    return dict(exact=c["exact"] / n, over=c["overtrim"] / n,
                under=c["undertrim"] / n, err=err / n)


def score_mate(mt: MateTable, rule5, rule3, tol: int) -> dict:
    c = Counter()
    for (r5, r3), htc in mt.cells.items():
        c5, c3 = rule5(r5), rule3(r3)
        for (h, t), w in htc.items():
            c[classify(c5, c3, h, t, tol)] += w
    n = mt.total()
    return dict(exact=c["exact"] / n, over=c["overtrim"] / n,
                under=c["undertrim"] / n)


def ceiling_end(tab: EndTable, tol: int) -> float:
    """sum_r max_cut P(|h-cut| <= tol | r) P(r): the per-end ceiling."""
    tot = tab.total()
    acc = 0.0
    for r, hc in tab.joint.items():
        if tol == 0:
            acc += max(hc.values())
        else:
            best = 0.0
            for c in range(0, CMAX + 1):
                mass = sum(w for h, w in hc.items() if abs(h - c) <= tol)
                best = max(best, mass)
            acc += best
    return acc / tot


def ceiling_mate_joint(mt: MateTable, tol: int) -> float:
    """Ceiling for a rule reading BOTH ends' runs: per (r5,r3) cell pick the
    cut pair maximising joint exact mass.  Upper bound for end-composed rules."""
    tot = mt.total()
    acc = 0.0
    for cell, htc in mt.cells.items():
        if tol == 0:
            acc += max(htc.values())
            continue
        hs = sorted({h for h, _ in htc})
        ts = sorted({t for _, t in htc})
        h_lo, h_hi = hs[0], hs[-1]
        t_lo, t_hi = ts[0], ts[-1]
        best = 0.0
        for c5 in range(max(0, h_lo - tol), h_hi + tol + 1):
            for c3 in range(max(0, t_lo - tol), t_hi + tol + 1):
                mass = sum(w for (h, t), w in htc.items()
                           if abs(h - c5) <= tol and abs(t - c3) <= tol)
                if mass > best:
                    best = mass
        acc += best
    return acc / tot


def class_flow(tab: EndTable, rule_a, rule_b, tol: int) -> dict:
    """Transition counts between the classes induced by two rules (one end)."""
    flow: dict[tuple[str, str], float] = defaultdict(float)
    for r, hc in tab.joint.items():
        ca, cb = rule_a(r), rule_b(r)
        for h, w in hc.items():
            da, db = ca - h, cb - h

            def cls(d):
                if abs(d) <= tol:
                    return "exact"
                return "overtrim" if d > 0 else "undertrim"
            flow[(cls(da), cls(db))] += w
    return {f"{a}->{b}": w / tab.total() for (a, b), w in flow.items()}


# ----------------------------------------------------------------- arms
def build_fastsim_arm(outdir: str, kept: set[tuple[int, int]] | None = None
                      ) -> tuple[EndTable, EndTable, MateTable, int, dict]:
    gt = load_gt(os.path.join(outdir, "sim.v0.gt.tsv"))
    t5, t3, mt = EndTable("5'"), EndTable("3'"), MateTable()
    # per-mate splits: R1 and R2 carry DIFFERENT artefact priors (mid: R1 head
    # is all-zero, R2 head has the 1 bp dropout), so a rule that may also know
    # which mate it reads has a different posterior and a different ceiling.
    split = {"5p1": EndTable("5'/R1"), "5p2": EndTable("5'/R2"),
             "3p1": EndTable("3'/R1"), "3p2": EndTable("3'/R2")}
    seen: set[tuple[int, int]] = set()
    for mate, f in ((1, "sim.v0_R1.fastq.gz"), (2, "sim.v0_R2.fastq.gz")):
        path = os.path.join(outdir, f)
        if not os.path.exists(path):
            continue
        for pid, r5, r3 in scan_mates(path, mate):
            if pid not in gt:
                continue
            if kept is not None and (pid, mate - 1) not in kept:
                continue
            h1, t1, h2, t2 = gt[pid]
            h, t = (h1, t1) if mate == 1 else (h2, t2)
            t5.add(r5, h)
            t3.add(r3, t)
            split[f"5p{mate}"].add(r5, h)
            split[f"3p{mate}"].add(r3, t)
            mt.add(r5, r3, h, t)
            seen.add((pid, mate - 1))
    return t5, t3, mt, len(seen), split


def build_merge_arm(profile_path: str, n: float = 1.0) -> tuple[EndTable, EndTable, MateTable]:
    """ANALYTIC merge-permissive model: r = h + b, b ~ Geom(p_G), h/t drawn
    from the profile's head/tail marginals.  R2 priors (R1 head is all-zero in
    the mid profile).  Not fastsim GT — a sensitivity bound only."""
    with open(profile_path) as fh:
        p = json.load(fh)
    p_g = p.get("gc", 0.47) / 2.0

    def geom(k: int) -> float:
        return (1 - p_g) * (p_g ** k)

    def prior(hist: dict) -> dict[int, float]:
        tot = sum(hist.values())
        return {int(k): v / tot for k, v in hist.items()}

    h_pri = prior(p["head"]["R2"])
    t_pri = prior(p["tail"]["R2"])
    t5, t3, mt = EndTable("5'"), EndTable("3'"), MateTable()
    for r in range(0, RMAX + 1):
        for h, ph in h_pri.items():
            if 0 <= r - h <= RMAX and r >= h:
                w = ph * geom(r - h) * n
                if w > 0:
                    t5.add(r, h, w)
        for t_, pt in t_pri.items():
            if 0 <= r - t_ <= RMAX and r >= t_:
                w = pt * geom(r - t_) * n
                if w > 0:
                    t3.add(r, t_, w)
    # mate cells = independent product of the two ends
    p5 = t5.p_r()
    p3 = t3.p_r()
    post5, post3 = t5.post(), t3.post()
    for r5, q5 in p5.items():
        for r3, q3 in p3.items():
            for h, a in post5[r5].items():
                for t_, b in post3[r3].items():
                    mt.add(r5, r3, h, t_, q5 * q3 * a * b * n)
    return t5, t3, mt


# ----------------------------------------------------------------- fastprep ref
def parse_pxmt(path: str) -> dict[tuple[int, int], tuple[int, int, int]]:
    out: dict[tuple[int, int], tuple[int, int, int]] = {}
    buf = open(path, "rb").read()
    if buf[:4] != b"PXMT":
        return out
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
        nr = buf[off]
        off += 1
        off += 8 * nr
        parts = rid.split(b" ", 1)[0].decode("ascii").lstrip("@").split(":")
        pid = int(parts[-1]) * 10_000 + int(parts[-2])
        out[(pid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def fastprep_reference(outdir: str, gt: dict) -> dict | None:
    if not os.path.exists(FASTPREP):
        return None
    pxmt = os.path.join(outdir, "fastprep_default.pxmt")
    for f in (pxmt, os.path.join(outdir, "z1.fastq.gz"),
              os.path.join(outdir, "z2.fastq.gz"), os.path.join(outdir, "z.json")):
        if os.path.exists(f):
            os.unlink(f)
    argv = [FASTPREP,
            "-1", os.path.join(outdir, "sim.v0_R1.fastq.gz"),
            "-2", os.path.join(outdir, "sim.v0_R2.fastq.gz"),
            "-o", os.path.join(outdir, "z1.fastq.gz"),
            "-p", os.path.join(outdir, "z2.fastq.gz"),
            "--metadata", pxmt, "--json", os.path.join(outdir, "z.json"),
            "--force", "-q", "-j", "4"] + FP_BASE
    p = subprocess.run(argv, capture_output=True, text=True)
    if p.returncode != 0:
        print(f"  fastprep reference failed: {p.stderr[-200:]}", file=sys.stderr)
        return None
    px = parse_pxmt(pxmt)
    # empirical policy: observed run r -> cut actually applied
    pol5: dict[int, Any] = defaultdict(Counter)
    pol3: dict[int, Any] = defaultdict(Counter)
    cls: dict[int, Any] = {0: Counter(), 2: Counter()}
    kept_ids: set[tuple[int, int]] = set()
    for mate, f in ((1, "sim.v0_R1.fastq.gz"), (2, "sim.v0_R2.fastq.gz")):
        for pid, r5, r3 in scan_mates(os.path.join(outdir, f), mate):
            s = px.get((pid, mate - 1))
            if s is None or s[0] != 0:
                continue
            kept_ids.add((pid, mate - 1))
            _, c5, c3 = s
            pol5[min(r5, RMAX)][min(c5, RMAX)] += 1
            pol3[min(r3, RMAX)][min(c3, RMAX)] += 1
            h1, t1, h2, t2 = gt.get(pid, (0, 0, 0, 0))
            h, t = (h1, t1) if mate == 1 else (h2, t2)
            for tol in TOLS:
                cls[tol][classify(c5, c3, h, t, tol)] += 1
    out = {"kept": len(kept_ids)}
    for tol in TOLS:
        n = max(1, sum(cls[tol].values()))
        out[f"tol{tol}"] = {k: cls[tol][k] / n for k in ("exact", "overtrim", "undertrim")}
    out["policy5"] = {r: dict(c) for r, c in sorted(pol5.items())}
    out["policy3"] = {r: dict(c) for r, c in sorted(pol3.items())}
    return out


# ----------------------------------------------------------------- report
def fmt_pct(x: float) -> str:
    return f"{100.0 * x:6.2f}%"


def print_cond(tab: EndTable, label: str) -> None:
    post = tab.post()
    p_r = tab.p_r()
    hs = sorted({h for c in tab.joint.values() for h in c})
    cols = [h for h in hs if h <= 8] + ([">8"] if any(h > 8 for h in hs) else [])
    print(f"\n  P(true artefact | observed run r) — {label}")
    print(f"  {'r':>3s} {'P(r)':>8s} " + " ".join(f"{str(c):>6s}" for c in cols) +
          f"   {'argmax':>6s} {'E[h|r]':>7s}")
    for r in range(0, RMAX + 1):
        if r not in post:
            continue
        p = post[r]
        cells = []
        for c in cols:
            if c == ">8":
                cells.append(sum(w for h, w in p.items() if h > 8))
            else:
                cells.append(p.get(c, 0.0))
        am = mode_cut(p)
        em = sum(h * w for h, w in p.items())
        tag = f"{r}" if r < RMAX else f"{r}+"
        print(f"  {tag:>3s} {p_r[r]:8.5f} " +
              " ".join(f"{v:6.3f}" for v in cells) +
              f"   {am:6d} {em:7.3f}")


def print_rules(rows: list[tuple[str, dict, dict]], title: str) -> None:
    print(f"\n  {title}")
    print(f"  {'rule':26s} {'ex@0':>7s} {'ov@0':>7s} {'un@0':>7s} | "
          f"{'ex@2':>7s} {'ov@2':>7s} {'un@2':>7s} | {'err@0':>7s}")
    for name, r0, r2 in rows:
        err = f"{r0['err']:+7.3f}" if "err" in r0 else "      -"
        print(f"  {name:26s} {fmt_pct(r0['exact'])} {fmt_pct(r0['over'])} "
              f"{fmt_pct(r0['under'])} | {fmt_pct(r2['exact'])} {fmt_pct(r2['over'])} "
              f"{fmt_pct(r2['under'])} | {err}")


def analyse_arm(name: str, t5: EndTable, t3: EndTable, mt: MateTable,
                results: dict, fp: dict | None = None,
                split: dict | None = None) -> None:
    print(f"\n===== arm: {name} =====  mates scored: {mt.total():.0f}")
    print_cond(t5, "5' end, h = true head artefact")
    print_cond(t3, "3' end, t = true tail artefact")

    per_end = {}
    for lab, tab in (("5p", t5), ("3p", t3)):
        rows = []
        for tol in TOLS:
            for rname, rule in make_rules(tab, tol):
                rows.append((tol, rname, rule))
        seen = set()
        out_rows = []
        for tol, rname, rule in rows:
            if rname in seen:
                continue
            seen.add(rname)
            r0 = score_end(tab, rule, 0)
            r2 = score_end(tab, rule, 2)
            out_rows.append((rname, r0, r2))
        print_rules(out_rows, f"{lab} end rules (all mates)")
        per_end[lab] = {n: {"tol0": a, "tol2": b} for n, a, b in out_rows}

    # mate-joint: same rule family applied to both ends
    joint_rows = []
    rule_names = [n for n, _ in make_rules(t5, 0)]
    for rname in rule_names:
        r5 = dict((n, f) for n, f in make_rules(t5, 0)).get(rname) or \
             dict((n, f) for n, f in make_rules(t5, 2)).get(rname)
        r3 = dict((n, f) for n, f in make_rules(t3, 0)).get(rname) or \
             dict((n, f) for n, f in make_rules(t3, 2)).get(rname)
        r0 = score_mate(mt, r5, r3, 0)
        r2 = score_mate(mt, r5, r3, 2)
        joint_rows.append((rname, r0, r2))
    print_rules(joint_rows, "mate-joint rules (same f on both ends)")

    ceil5_0, ceil5_2 = ceiling_end(t5, 0), ceiling_end(t5, 2)
    ceil3_0, ceil3_2 = ceiling_end(t3, 0), ceiling_end(t3, 2)
    cj0, cj2 = ceiling_mate_joint(mt, 0), ceiling_mate_joint(mt, 2)
    print(f"\n  ceilings (any function of the observed run alone):")
    print(f"    5' tol0 {fmt_pct(ceil5_0)}   3' tol0 {fmt_pct(ceil3_0)}")
    print(f"    5' tol2 {fmt_pct(ceil5_2)}   3' tol2 {fmt_pct(ceil3_2)}")
    print(f"    mate-joint tol0 {fmt_pct(cj0)}  tol2 {fmt_pct(cj2)}")
    split_ceil = {}
    if split:
        c51, c52 = ceiling_end(split["5p1"], 0), ceiling_end(split["5p2"], 0)
        c31, c32 = ceiling_end(split["3p1"], 0), ceiling_end(split["3p2"], 0)
        w51 = split["5p1"].total() or 1.0
        w52 = split["5p2"].total() or 1.0
        m5 = (c51 * w51 + c52 * w52) / (w51 + w52)
        w31 = split["3p1"].total() or 1.0
        w32 = split["3p2"].total() or 1.0
        m3 = (c31 * w31 + c32 * w32) / (w31 + w32)
        print(f"    mate-aware (rule may know R1 vs R2): "
              f"5' tol0 {fmt_pct(m5)} (R1 {fmt_pct(c51)}, R2 {fmt_pct(c52)})   "
              f"3' tol0 {fmt_pct(m3)} (R1 {fmt_pct(c31)}, R2 {fmt_pct(c32)})")
        split_ceil = {"end5_R1_tol0": c51, "end5_R2_tol0": c52,
                      "end3_R1_tol0": c31, "end3_R2_tol0": c32,
                      "end5_mateaware_tol0": m5, "end3_mateaware_tol0": m3}

    cur = dict((n, (a, b)) for n, a, b in joint_rows)["cut=r (current)"]
    bay = dict((n, (a, b)) for n, a, b in joint_rows).get("bayes-mode", cur)
    print(f"\n  current f(r)=r mate-joint: tol0 exact {fmt_pct(cur[0]['exact'])} "
          f"over {fmt_pct(cur[0]['over'])} under {fmt_pct(cur[0]['under'])}")
    print(f"  bayes-mode    mate-joint: tol0 exact {fmt_pct(bay[0]['exact'])} "
          f"over {fmt_pct(bay[0]['over'])} under {fmt_pct(bay[0]['under'])}")

    flow5 = class_flow(t5, lambda r: r,
                       dict((n, f) for n, f in make_rules(t5, 0))["bayes-mode"], 0)
    flow3 = class_flow(t3, lambda r: r,
                       dict((n, f) for n, f in make_rules(t3, 0))["bayes-mode"], 0)
    print(f"\n  class flow f(r)=r -> bayes-mode at tol 0:")
    print(f"    5': " + "  ".join(f"{k} {fmt_pct(v)}" for k, v in sorted(flow5.items())))
    print(f"    3': " + "  ".join(f"{k} {fmt_pct(v)}" for k, v in sorted(flow3.items())))

    if fp:
        print(f"\n  fastprep default config on this dataset (kept {fp['kept']} mates):")
        for tol in TOLS:
            d = fp[f"tol{tol}"]
            print(f"    tol{tol}: exact {fmt_pct(d['exact'])} over {fmt_pct(d['overtrim'])} "
                  f"under {fmt_pct(d['undertrim'])}")
        for lab in ("policy5", "policy3"):
            pol = fp[lab]
            modal = {r: max(c.items(), key=lambda kv: kv[1])[0] for r, c in pol.items()
                     if sum(c.values()) >= 20}
            print(f"    {lab} modal cut by observed run (n>=20): {modal}")

    results[name] = {
        "mates": mt.total(),
        "ceilings": {"end5_tol0": ceil5_0, "end5_tol2": ceil5_2,
                     "end3_tol0": ceil3_0, "end3_tol2": ceil3_2,
                     "joint_tol0": cj0, "joint_tol2": cj2, **split_ceil},
        "per_end_rules": per_end,
        "joint_rules": {n: {"tol0": a, "tol2": b} for n, a, b in joint_rows},
        "flow5": flow5, "flow3": flow3,
        "fastprep": fp,
    }


# ----------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--regimes", nargs="+", default=["mid"])
    ap.add_argument("--merge", action="store_true", default=True,
                    help="also run the analytic merge-permissive arm (default on)")
    ap.add_argument("--no-merge", dest="merge", action="store_false")
    ap.add_argument("--fastprep", action="store_true", default=True)
    ap.add_argument("--no-fastprep", dest="fastprep", action="store_false")
    args = ap.parse_args()

    ensure_fixtures()
    os.makedirs(WORK, exist_ok=True)
    results: dict[str, Any] = {"pairs": args.pairs, "seed": args.seed}

    for regime in args.regimes:
        outdir = simulate(regime, args.pairs, args.seed)
        gt = load_gt(os.path.join(outdir, "sim.v0.gt.tsv"))
        t5, t3, mt, n, split = build_fastsim_arm(outdir)
        print(f"[{regime}] fastsim: {len(gt)} pairs in GT, {n} mate observations scored")
        fp = fastprep_reference(outdir, gt) if args.fastprep else None
        analyse_arm(f"fastsim/{regime}", t5, t3, mt, results, fp, split)

        if args.merge:
            m5, m3, mm = build_merge_arm(os.path.join(FIXTURES, f"{regime}.json"),
                                         n=max(1.0, mt.total()))
            analyse_arm(f"merge/{regime}", m5, m3, mm, results)

    out_path = os.path.join(WORK, "ceiling.json")
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=1, default=str)
    print(f"\nwrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
