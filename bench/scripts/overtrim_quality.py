#!/usr/bin/env python3
"""Can base quality discriminate an artefact G-run from biology at a G-run
boundary well enough to reduce fastprep's strict-tolerance overtrim?

Reproduces the full analysis in docs/OVERTRIM-QUALITY.md. Four parts:

  model   Audit fastsim's quality model at artefact runs: what quality does
          the simulator actually give run positions (first head base,
          positions 2..h, the whole tail run) versus biology?

  probe   Controlled probe grid that survives fastprep's length budget
          (the earlier /tmp/probe_q_split.py failed because its reads were
          budget-dropped before any boundary decision). Short runs (3-6 bp),
          long reads, generous --min-len so the budget never fires, and PXMT
          status bytes checked per read. Varies run length, artefact
          quality, first-biology-base quality and interior quality, at both
          ends. Measures (a) whether fastprep's cut moves with quality at
          all, (b) whether the cut can ever land inside the observed run.

  rule    Offline evaluation of a quality-boundary rule (trim only the
          low-quality segment of the end G-run) against fastsim ground
          truth, in three quality worlds: as simulated (weak signal),
          artefact-quality rewritten low (2-colour "no signal" world, the
          world where quality SHOULD help), and a 3'-decay world (biology
          quality decays toward the read end: the rule's false-positive
          stress test). Reports exact/overtrim/undertrim at tol 0 and 2
          versus fastprep's current cuts, plus the false-positive cost.

  cost    Throughput before/after adding the rule to fastprep itself
          (--q-boundary, default off) on a synthetic 100k-pair file.

Usage:
  python3 bench/scripts/overtrim_quality.py --part all
  python3 bench/scripts/overtrim_quality.py --part probe --verbose

Everything runs locally against the in-repo fastsim/fastprep binaries.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import json
import os
import random
import statistics
import struct
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FASTPREP = os.environ.get(
    "FASTPREP", os.path.join(REPO, "tools/fastprep/target/release/fastprep"))
FASTSIM = os.environ.get(
    "FASTSIM", os.path.join(REPO, "tools/fastsim/target/release/fastsim"))
FIXTURES = os.path.join(REPO, "bench/fixtures/overtrim")
WORK = os.environ.get(
    "OVERTRIM_Q_WORK", "/home/raphael/.hermes/cache/scratch/overtrim-quality")

# The Phase-D reference numbers (desktop, 101 fitted real-library profiles,
# bench/results/phase-d/trim-accuracy.tsv, fastprep-winner, scored per-mate):
PHASE_D = {0: dict(exact=66.5, over=29.5, under=4.0),
           2: dict(exact=98.4, over=1.3, under=0.3)}

# fastprep config mirroring bench/scripts/overtrim_sensitivity.py BASE.
BASE_CFG = ["--min-run", "12", "--mercy", "3", "--phred-sat", "20",
            "--gate-shape", "1.0", "--min-len", "0.9L", "--mean-q", "25.0",
            "--head-dropout", "1"]


# --------------------------------------------------------------- primitives
def write_fastq(path, reads):
    """reads: iterable of (name, seq, quals) with raw Phred ints."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "wb") as fh:
        for name, seq, quals in reads:
            fh.write(name.encode() + b"\n" + seq.encode() + b"\n+\n")
            fh.write(bytes(q + 33 for q in quals) + b"\n")


def parse_pxmt(path):
    """-> {read_id_core: (status, t5, t3, t5a, t3a, t5p, t3p)} (mate-agnostic;
    probe read names are unique per read)."""
    out = {}
    if not os.path.exists(path):
        return out
    buf = open(path, "rb").read()
    if buf[:4] != b"PXMT":
        return out
    off = 5
    while off + 2 <= len(buf):
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len].decode("utf-8", "replace"); off += id_len
        if off + 22 > len(buf):
            break
        mate = buf[off] + 1; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        (_rl,) = struct.unpack_from("<I", buf, off); off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        key = rid.split()[0]
        out[key] = (status, t5a + t5p, t3a + t3p, t5a, t3a, t5p, t3p)
    return out


def parse_pxmt_paired(path):
    """-> {(pair_id, mate): (status, t5, t3, t5p, t3p)} keyed by fastsim
    coordinate pair_id = x + 10000*y (see bench/scripts/trim_accuracy.py)."""
    out = {}
    if not os.path.exists(path):
        return out
    buf = open(path, "rb").read()
    if buf[:4] != b"PXMT":
        return out
    off = 5
    while off + 2 <= len(buf):
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len]; off += id_len
        if off + 22 > len(buf):
            break
        mate = buf[off] + 1; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        (_rl,) = struct.unpack_from("<I", buf, off); off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        try:
            fields = rid.split(b" ", 1)[0].decode("ascii").lstrip("@").split(":")
            x, y = int(fields[-2]), int(fields[-1])
            pid = y * 10_000 + x
        except (UnicodeDecodeError, IndexError, ValueError):
            continue
        out[(pid, mate)] = (status, t5a + t5p, t3a + t3p, t5p, t3p)
    return out


DROP_REASONS = {0: "kept", 1: "drop:hard", 2: "drop:poly5", 3: "drop:poly3",
                4: "drop:poly_int", 5: "drop:meanq", 6: "drop:adapter",
                7: "drop:duplicate"}


def run_fastprep(reads_path, out_prefix, extra, paired=False):
    """Run fastprep; return parsed PXMT. Dedup/mean-q(when asked)/auto-phred
    off; adapter off unless the caller passes it."""
    r1 = reads_path if not paired else reads_path + "_R1.fastq.gz"
    argv = [FASTPREP, "-1", r1]
    if paired:
        argv += ["-2", reads_path + "_R2.fastq.gz",
                 "-p", out_prefix + "_R2.fastq.gz"]
    argv += ["-o", out_prefix + "_R1.fastq.gz", "--metadata", out_prefix + ".pxmt",
             "--force", "-q", "--no-dedup", "--no-auto-phred", "--no-adapter",
             "-j", "2"]
    argv += extra
    rc = subprocess.run(argv, capture_output=True, text=True)
    if rc.returncode != 0:
        raise SystemExit(f"fastprep failed rc={rc.returncode}: {rc.stderr[-400:]}")
    return (parse_pxmt_paired(out_prefix + ".pxmt") if paired
            else parse_pxmt(out_prefix + ".pxmt"))


def classify(t5, t3, h, t, tol):
    if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
        return "exact"
    if t5 > h or t3 > t:
        return "overtrim"
    return "undertrim"


def score_rows(rows, tol):
    """rows: iterable of (t5, t3, h, t). -> rates over scored rows."""
    c = collections.Counter()
    over_bp = under_bp = 0
    for t5, t3, h, t in rows:
        c[classify(t5, t3, h, t, tol)] += 1
        if t5 > h or t3 > t:
            over_bp += max(t5 - h, 0) + max(t3 - t, 0)
        if t5 < h or t3 < t:
            under_bp += max(h - t5, 0) + max(t - t3, 0)
    n = sum(c.values()) or 1
    return dict(exact=c["exact"] / n, over=c["overtrim"] / n,
                under=c["undertrim"] / n, n=n,
                over_bp=over_bp / n, under_bp=under_bp / n)


def fmt(r):
    return (f"exact {r['exact']*100:5.1f}%  over {r['over']*100:5.1f}%  "
            f"under {r['under']*100:5.1f}%  (n={r['n']}, "
            f"over {r['over_bp']:.2f} bp/read, under {r['under_bp']:.2f} bp/read)")


# ------------------------------------------------- part 1: fastsim model
def simulate(profile_name, pairs, seed, out_prefix):
    prof = os.path.join(FIXTURES, profile_name + ".json")
    argv = [FASTSIM, "--pairs", str(pairs), "--profile", prof,
            "--out", out_prefix, "--seed", str(seed), "--no-adapter",
            "-f", "-q"]
    rc = subprocess.run(argv, capture_output=True, text=True)
    if rc.returncode != 0:
        raise SystemExit(f"fastsim failed: {rc.stderr[-400:]}")
    return out_prefix


def load_gt(path):
    out = {}
    with open(path, newline="") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            pid = int(r["pair_id"])
            out[pid] = (int(r["r1_h"]), int(r["r1_t"]),
                        int(r["r2_h"]), int(r["r2_t"]),
                        int(r.get("full_flag", 0)))
    return out


def read_fastq_index(path):
    """-> {(pair_id, mate): (seq, quals)} for fastsim headers."""
    out = {}
    with gzip.open(path, "rt") as fh:
        lines = fh
        while True:
            h = lines.readline()
            if not h:
                break
            seq = lines.readline().rstrip("\n")
            lines.readline()
            qual = lines.readline().rstrip("\n")
            fields = h[1:].split()[0].split(":")
            mate = int(h.split()[1][0])
            pid = int(fields[-1]) * 10_000 + int(fields[-2])
            out[(pid, mate)] = (seq, [ord(c) - 33 for c in qual])
    return out


def part_model(work, verbose):
    print("=== PART 1: what quality does fastsim give artefact runs? ===")
    prefix = os.path.join(work, "model", "sim")
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
    if not os.path.exists(prefix + ".gt.tsv"):
        simulate("mid", 20000, 7, prefix)
    gt = load_gt(prefix + ".gt.tsv")
    reads = {}
    for m in (1, 2):
        reads.update(read_fastq_index(f"{prefix}_R{m}.fastq.gz"))

    # positions of interest per mate: head run [0,h), tail run [L-t, L),
    # biology = [h, L-t). Report mean Q per run position and at the two
    # boundary bases. Split tails by length so the boundary base quality is
    # measured where the overtrim question actually lives.
    acc = collections.defaultdict(list)
    for (pid, mate), (seq, q) in reads.items():
        if gt[pid][4] == 1:
            continue
        h, t = gt[pid][0:2] if mate == 1 else gt[pid][2:4]
        L = len(seq)
        if h > 0:
            acc["head_pos0"].append(q[0])
        if h > 1:
            acc["head_pos1..h-1"].extend(q[1:h])
        if t > 0:
            acc["tail_run"].extend(q[L - t:])
            acc["tail_first_base"].append(q[L - t])
            acc["tail_last_base"].append(q[L - 1])
        if h < L - t:
            acc["biology"].extend(q[h:L - t])
        if t > 0 and L - t - 1 >= h:
            acc["bio_base_before_tail"].append(q[L - t - 1])
        if h > 0 and h < L - t:
            acc["bio_base_after_head"].append(q[h])
        if t >= 2:
            acc["tail_pos_L-t+1.."].extend(q[L - t + 1:])
    print(f"  dataset: {sum(1 for _ in gt)} pairs, adapter-free, "
          f"profile mid.json (tail_q_collapse = 3.0, pos_q = 37 flat)")
    for k in ("head_pos0", "head_pos1..h-1", "tail_run", "tail_pos_L-t+1..",
              "tail_first_base", "tail_last_base", "bio_base_after_head",
              "bio_base_before_tail", "biology"):
        v = acc.get(k, [])
        if v:
            print(f"  {k:24s} n={len(v):7d}  mean Q={statistics.mean(v):5.2f}  "
                  f"sd={statistics.pstdev(v):4.2f}  min={min(v)} max={max(v)}")
    # Verdict lines the doc quotes.
    if acc.get("head_pos0") and acc.get("head_pos1..h-1"):
        d = statistics.mean(acc["head_pos0"]) - statistics.mean(acc["head_pos1..h-1"])
        print(f"  -> first head base minus deeper head bases: {d:+.2f} Phred")
    if acc.get("tail_run") and acc.get("biology"):
        d = statistics.mean(acc["tail_run"]) - statistics.mean(acc["biology"])
        print(f"  -> tail-run bases minus biology: {d:+.2f} Phred")

    # The mid profile's head CDF only has mass at h in {0,1}, so positions
    # 2..h never occur there. Author a head-probe profile (same quality
    # model, head CDF at h = 2..6) and audit those positions directly.
    hp = os.path.join(work, "model", "headprobe.json")
    prof = json.load(open(os.path.join(FIXTURES, "mid.json")))
    prof["head"] = {"R1": {"2": 400, "3": 400, "4": 400, "5": 400, "6": 400},
                    "R2": {"2": 400, "3": 400, "4": 400, "5": 400, "6": 400}}
    prof["tail"] = {"R1": {"0": 10000}, "R2": {"0": 10000}}
    prof["joint_tail"] = {"0,0": 10000}
    json.dump(prof, open(hp, "w"))
    prefix2 = os.path.join(work, "model", "heads")
    if not os.path.exists(prefix2 + ".gt.tsv"):
        rc = subprocess.run([FASTSIM, "--pairs", "10000", "--profile", hp,
                             "--out", prefix2, "--seed", "8", "--no-adapter",
                             "-f", "-q"], capture_output=True, text=True)
        if rc.returncode != 0:
            raise SystemExit(f"fastsim(headprobe) failed: {rc.stderr[-300:]}")
    gt2 = load_gt(prefix2 + ".gt.tsv")
    reads2 = {}
    for m in (1, 2):
        reads2.update(read_fastq_index(f"{prefix2}_R{m}.fastq.gz"))
    a0, adeep, abio = [], [], []
    for (pid, mate), (seq, q) in reads2.items():
        h = gt2[pid][0] if mate == 1 else gt2[pid][2]
        if h >= 2:
            a0.append(q[0])
            adeep.extend(q[1:h])
            abio.extend(q[h:h + 3])
    print(f"  head-probe profile (h = 2..6), no deep_head_q in profile: "
          f"pos0 n={len(a0)} mean Q={statistics.mean(a0):.2f}; "
          f"positions 2..h n={len(adeep)} mean Q={statistics.mean(adeep):.2f}; "
          f"first biology bases mean Q={statistics.mean(abio):.2f}")

    # Does quality decay through the run? Measure tail-run quality by offset
    # from the run start (boundary side) and from the read end.
    by_off = collections.defaultdict(list)
    for (pid, mate), (seq, q) in reads.items():
        if gt[pid][4] == 1:
            continue
        t = gt[pid][1] if mate == 1 else gt[pid][3]
        L = len(seq)
        for i in range(t):
            off = min(i, 3)          # offset from run start (boundary side)
            by_off[off].append(q[L - t + i])
    print("  tail-run quality by offset from the run's boundary base: "
          + "  ".join(f"off{o}: {statistics.mean(v):.2f}"
                      for o, v in sorted(by_off.items())))

    # Boundary-pair separation: is the last artefact base distinguishable
    # from the first biology base at the same boundary?
    sep_t = sep_h = sep_h1 = 0
    n_t = n_h = n_h1 = 0
    for (pid, mate), (seq, q) in reads.items():
        if gt[pid][4] == 1:
            continue
        h, t = gt[pid][0:2] if mate == 1 else gt[pid][2:4]
        L = len(seq)
        if t > 0 and L - t - 1 >= h:
            n_t += 1
            if q[L - t] < q[L - t - 1]:
                sep_t += 1
        if h >= 2:
            n_h += 1
            if q[h - 1] < q[h]:
                sep_h += 1
        if h == 1 and L > 1:
            n_h1 += 1
            if q[0] < q[1]:
                sep_h1 += 1
    print(f"  P(last artefact base < first biology base) -- the boundary "
          f"step the rule needs:")
    print(f"    tail boundary (t>0):        {sep_t}/{n_t} = "
          f"{sep_t/max(1,n_t)*100:.1f}%")
    print(f"    head boundary (h>=2):       {sep_h}/{n_h} = "
          f"{sep_h/max(1,n_h)*100:.1f}%")
    print(f"    head boundary (h=1 only):   {sep_h1}/{n_h1} = "
          f"{sep_h1/max(1,n_h1)*100:.1f}%  <- the dominant 5' case; no signal")

    # Same, but with the measured deep-head quality present (as the real
    # fitted profiles carry): positions 2..h must then collapse while
    # position 0 stays normal.
    hp2 = os.path.join(work, "model", "headprobe_deep.json")
    prof["deep_head_q"] = {"R1": 25.0, "R2": 25.0}
    prof["run_q_mean"] = {"R1": {"head": 25.0, "tail": 25.0},
                          "R2": {"head": 25.0, "tail": 25.0}}
    json.dump(prof, open(hp2, "w"))
    prefix3 = os.path.join(work, "model", "heads_deep")
    if not os.path.exists(prefix3 + ".gt.tsv"):
        rc = subprocess.run([FASTSIM, "--pairs", "10000", "--profile", hp2,
                             "--out", prefix3, "--seed", "8", "--no-adapter",
                             "-f", "-q"], capture_output=True, text=True)
        if rc.returncode != 0:
            raise SystemExit(f"fastsim(headprobe_deep) failed: {rc.stderr[-300:]}")
    gt3 = load_gt(prefix3 + ".gt.tsv")
    reads3 = {}
    for m in (1, 2):
        reads3.update(read_fastq_index(f"{prefix3}_R{m}.fastq.gz"))
    b0, bdeep, bbio = [], [], []
    for (pid, mate), (seq, q) in reads3.items():
        h = gt3[pid][0] if mate == 1 else gt3[pid][2]
        if h >= 2:
            b0.append(q[0])
            bdeep.extend(q[1:h])
            bbio.extend(q[h:h + 3])
    print(f"  head-probe profile (h = 2..6), deep_head_q=25: "
          f"pos0 n={len(b0)} mean Q={statistics.mean(b0):.2f}; "
          f"positions 2..h n={len(bdeep)} mean Q={statistics.mean(bdeep):.2f}; "
          f"first biology bases mean Q={statistics.mean(bbio):.2f}")
    print(f"  -> the model's head contract: position 0 of the run is ALWAYS "
          f"normal quality; positions 2..h get deep_head_q when the profile "
          f"carries it (run_q_mean.head / deep_head_q) and mean(pos_q) -- i.e. "
          f"NO collapse -- when it does not (the bench fixtures do not).")
    return acc


# ------------------------------------------------- part 2: probe grid
def part_probe(work, verbose):
    print()
    print("=== PART 2: controlled probe grid (budget-safe) ===")
    rnd = random.Random(20261016)
    NO_G = "ACT"

    # grid: geometry x run x k(bio G's at boundary) x q_art x q_bio x q_int
    runs = (3, 4, 5, 6)
    ks = (0, 1, 2)
    q_arts = (2, 10, 20, 30, 40)
    q_bios = (2, 10, 20, 30, 40)
    q_ints = (15, 35)
    body_len = 120  # 3-6 + 2 + 120 = 128..128bp reads; min-len 30 -> budget >= 98

    specs = []
    for geo in ("head", "tail"):
        for rlen in runs:
            for k in ks:
                for qa in q_arts:
                    for qb in q_bios:
                        for qi in q_ints:
                            specs.append((geo, rlen, k, qa, qb, qi))
    reps = 4
    reads = []
    truth = {}
    for idx, (geo, rlen, k, qa, qb, qi) in enumerate(specs):
        for rep in range(reps):
            body = "".join(rnd.choice(NO_G) for _ in range(body_len))
            bio_g = "G" * k
            run = "G" * rlen
            if geo == "head":
                seq = run + bio_g + body
                quals = [qa] * rlen + [qb] * k + [qi] * body_len
                h, t = rlen, 0
            else:
                seq = body + bio_g + run
                quals = [qi] * body_len + [qb] * k + [qa] * rlen
                h, t = 0, rlen
            name = (f"@P{idx}_{rep}:{geo}:{rlen}:{k}:{qa}:{qb}:{qi}"
                    f" 1:N:0:ATCGATCG")
            rid = name[1:].split()[0]
            truth[rid] = dict(geo=geo, rlen=rlen, k=k, qa=qa, qb=qb, qi=qi,
                              h=h, t=t, R=(rlen + k))
            reads.append((name, seq, quals))
    os.makedirs(os.path.join(work, "probe"), exist_ok=True)
    fq = os.path.join(work, "probe", "grid.fastq")
    write_fastq(fq, reads)
    print(f"  {len(specs)} configs x {reps} reps = {len(reads)} reads; "
          f"runs {min(runs)}-{max(runs)} bp; reads 128 bp; --min-len 30 "
          f"(budget 98 bp, never fires)")

    for mercy in ("1", "3"):
        px = run_fastprep(fq, os.path.join(work, f"probe/g{mercy}"),
                          ["--min-run", "3", "--head-dropout", "1",
                           "--min-len", "30", "--mercy", mercy,
                           "--no-mean-q"])
        status = collections.Counter(v[0] for v in px.values())
        print(f"  --mercy {mercy}: PXMT status bytes: "
              + ", ".join(f"{DROP_REASONS.get(k, k)}={v}"
                          for k, v in sorted(status.items()))
              + f"  (of {len(truth)} reads)")
        missing = len(truth) - len(px)
        if missing:
            print(f"  WARNING: {missing} reads missing from PXMT")

        # per-read cut analysis, scored only on kept reads
        by = collections.defaultdict(list)   # (geo,qa,qb) -> list[err_vs_run]
        overs = collections.Counter()
        split_seen = 0
        n_scored = 0
        run_only = 0
        for rid, gt_ in truth.items():
            v = px.get(rid)
            if v is None or v[0] != 0:
                continue
            n_scored += 1
            status_, t5, t3 = v[0], v[1], v[2]
            cut = t5 if gt_["geo"] == "head" else t3
            err_run = cut - gt_["R"]      # >0: overshot past observed run
            # boundary error vs truth: for head, truth h = rlen; biological
            # G's at the boundary are GT biology -> overtrim when cut > rlen.
            err_truth = cut - gt_["rlen"]
            by[(gt_["geo"], gt_["qa"], gt_["qb"])].append(err_run)
            overs[err_truth] += 1
            if cut < gt_["R"]:
                split_seen += 1
            if cut == gt_["R"]:
                run_only += 1
        print(f"  scored {n_scored} kept reads (status==0); "
              f"cut inside observed run (would be a quality/base split): "
              f"{split_seen}; cut == observed run exactly: {run_only}")
        hist = ", ".join(f"{k:+d}:{v}" for k, v in sorted(overs.items()))
        print(f"  cut - truth_run_length histogram: {hist}")

        # does the cut move with quality? mean overshoot (cut - observed run)
        # grouped by artefact quality and boundary-base quality.
        print(f"  mean (cut - observed_run) by [geo, q_art, q_bio_first]:")
        table = collections.defaultdict(list)
        for (geo, qa, qb), vals in by.items():
            table[(geo, qa, qb)].extend(vals)
        for geo in ("head", "tail"):
            print(f"    {geo}: " + "  ".join(
                f"qa={qa}:[" + ",".join(
                    f"{statistics.mean(table[(geo,qa,qb)]):+.2f}"
                    for qb in q_bios) + "]"
                for qa in q_arts) + "   (rows qa, cols qb)")
        # cleanest single-number answer: variance of the cut explained by q
        head_cuts = [v for (geo, qa, qb), vals in by.items() if geo == "head"
                     for v in vals]
        if head_cuts:
            print(f"    head cuts: mean {statistics.mean(head_cuts):+.2f}, "
                  f"sd {statistics.pstdev(head_cuts):.2f} vs observed run")
    return truth, px


# ------------------------------------------------- part 3: rule evaluation
def quality_run_cut(seq, q, q_cut, end):
    """The quality-boundary rule: length of the maximal end G-run segment
    whose bases all have q < q_cut. `end` is 'left' or 'right'."""
    n = len(seq)
    if end == "left":
        i = 0
        while i < n and seq[i] == "G" and q[i] < q_cut:
            i += 1
        return i
    i = n - 1
    c = 0
    while i >= 0 and seq[i] == "G" and q[i] < q_cut:
        c += 1
        i -= 1
    return c


def part_rule(work, verbose):
    print()
    print("=== PART 3: quality-boundary rule vs fastprep (fastsim GT) ===")
    datasets = []
    for name in ("low", "mid", "high"):
        prefix = os.path.join(work, "rule", name)
        os.makedirs(prefix, exist_ok=True)
        sim = os.path.join(prefix, "sim")
        if not os.path.exists(sim + ".gt.tsv"):
            simulate(name, 20000, 11, sim)
        datasets.append((name, sim))

    q_cuts = (2, 5, 10, 15, 20, 25, 30, 33, 35, 36, 37, 38, 40)
    results = {}
    for name, sim in datasets:
        gt = load_gt(sim + ".gt.tsv")
        reads = {}
        for m in (1, 2):
            reads.update(read_fastq_index(f"{sim}_R{m}.fastq.gz"))

        # ---- quality worlds -------------------------------------------
        # as-is:     FASTSIM'S OWN MODEL output, unmodified. The headline:
        #            any separation here is fastsim's, not the probe's.
        # strong-Q2: probe construction -- artefact run bases forced to Q2
        #            (2-colour "no signal" world: the signal SHOULD exist).
        # decay:     probe construction -- biology decays 37->22 toward the
        #            3' end, artefact = local quality - 12.
        # overlap:   probe construction stress -- biology decays 37->25 with
        #            deterministic +-3 jitter, artefact = local - 8 with the
        #            same jitter: the two quality distributions OVERLAP.
        def world_as_is(seq, q, h, t, jit):
            return q

        def world_strong(seq, q, h, t, jit):
            qq = list(q)
            for i in range(h):
                qq[i] = 2
            for i in range(len(seq) - t, len(seq)):
                qq[i] = 2
            return qq

        def world_decay(seq, q, h, t, jit):
            L = len(seq)
            qq = [max(0, 37 - 15 * i // max(1, L - 1)) for i in range(L)]
            for i in list(range(h)) + list(range(L - t, L)):
                qq[i] = max(0, qq[i] - 12)
            return qq

        def world_overlap(seq, q, h, t, jit):
            L = len(seq)
            qq = [max(0, 37 - 12 * i // max(1, L - 1) + jit[i]) for i in range(L)]
            for i in list(range(h)) + list(range(L - t, L)):
                qq[i] = max(0, qq[i] - 8 + jit[i])
            return qq

        worlds = (("as-is(fastsim-model)", world_as_is),
                  ("strong-Q2(probe)", world_strong),
                  ("decay(probe)", world_decay),
                  ("overlap(probe)", world_overlap))

        for wname, wfun in worlds:
            # rewrite qualities per world and run fastprep ON the world's own
            # reads, so the baseline is measured in the same world as the rule.
            wq = {}
            for (pid, mate), (seq, q) in reads.items():
                h, t = gt[pid][0:2] if mate == 1 else gt[pid][2:4]
                jit = [((i * 2654435761 + pid * 40503 + mate * 97) % 7) - 3
                       for i in range(len(seq))]
                wq[(pid, mate)] = (seq, wfun(seq, q, h, t, jit))
            wdir = os.path.join(work, "rule", f"{name}_{wname.split('(')[0]}")
            os.makedirs(wdir, exist_ok=True)
            for m in (1, 2):
                rq = [(f"@LH00916:1:SIMFLOW1:1:1101:{pid%10000}:{pid//10000} "
                       f"{m}:N:0:ATCG", seq, q)
                      for (pid, mate), (seq, q) in sorted(wq.items())
                      if mate == m]
                write_fastq(os.path.join(wdir, f"w_R{m}.fastq.gz"), rq)
            px = run_fastprep(os.path.join(wdir, "w"), os.path.join(wdir, "o"),
                              BASE_CFG, paired=True)

            base_rows, per_read = [], []
            dropped = 0
            for (pid, mate), (seq, q) in wq.items():
                if gt[pid][4] == 1:      # full-artefact pairs: truth is DROP
                    continue
                h, t = gt[pid][0:2] if mate == 1 else gt[pid][2:4]
                v = px.get((pid, mate))
                if v is None or v[0] != 0:
                    dropped += 1
                    continue
                _, t5, t3, t5p, t3p = v
                base_rows.append((t5, t3, h, t))
                per_read.append((seq, q, t5p, t3p, h, t))
            print(f"\n  --- dataset {name} / world {wname}: "
                  f"{len(base_rows)} mates scored, {dropped} dropped")
            print(f"      fastprep now (tol 0): {fmt(score_rows(base_rows, 0))}")
            print(f"      fastprep now (tol 2): {fmt(score_rows(base_rows, 2))}")
            results[(name, wname, "base", "na")] = base_rows

            # rule cuts are per (read, q_cut); reuse across modes
            cuts = {}
            for i, (seq, q, t5p, t3p, h, t) in enumerate(per_read):
                for c in q_cuts:
                    cuts[(i, c)] = (quality_run_cut(seq, q, c, "left"),
                                    quality_run_cut(seq, q, c, "right"))
            for mode in ("replace", "refine"):
                for c in q_cuts:
                    rows = []
                    for i, (seq, q, t5p, t3p, h, t) in enumerate(per_read):
                        hq, tq = cuts[(i, c)]
                        if mode == "refine":
                            n5 = hq if hq > 0 else t5p
                            n3 = tq if tq > 0 else t3p
                        else:
                            n5, n3 = hq, tq
                        rows.append((n5, n3, h, t))
                    results[(name, wname, mode, c)] = rows

            show_sweep = wname.startswith("as-is") or wname.startswith("overlap")
            for mode in ("replace", "refine"):
                if mode == "refine" and not show_sweep:
                    continue
                if show_sweep:
                    print(f"      {mode} sweep (tol 0 exact/over/under | "
                          f"tol 2 exact):")
                for c in q_cuts:
                    r0 = score_rows(results[(name, wname, mode, c)], 0)
                    r2 = score_rows(results[(name, wname, mode, c)], 2)
                    if show_sweep:
                        print(f"        q_cut={c:2d}: "
                              f"{r0['exact']*100:5.1f}/{r0['over']*100:5.1f}/"
                              f"{r0['under']*100:5.1f} | "
                              f"{r2['exact']*100:5.1f}")

            best = max(q_cuts, key=lambda c: score_rows(
                results[(name, wname, "replace", c)], 0)["exact"])
            r0 = score_rows(results[(name, wname, "replace", best)], 0)
            r2 = score_rows(results[(name, wname, "replace", best)], 2)
            print(f"      rule best replace q_cut={best} (tol 0): {fmt(r0)}")
            print(f"      rule best replace q_cut={best} (tol 2): {fmt(r2)}")
            results[(name, wname, "best", best)] = None

            # False-positive cost, both directions, on mates fastprep had
            # EXACT at tol 0: (1) artefact left behind (undertrim), (2)
            # biology cut away where fastprep was right (overtrim).
            base_exact = [i for i, (t5, t3, h, t) in enumerate(base_rows)
                          if classify(t5, t3, h, t, 0) == "exact"]
            u_reads = o_reads = 0
            u_bp = o_bp = 0
            for i in base_exact:
                seq, q, t5p, t3p, h, t = per_read[i]
                hq, tq = cuts[(i, best)]
                cl = classify(hq, tq, h, t, 0)
                if cl == "undertrim":
                    u_reads += 1
                    u_bp += max(h - hq, 0) + max(t - tq, 0)
                elif cl == "overtrim":
                    o_reads += 1
                    o_bp += max(hq - h, 0) + max(tq - t, 0)
            if base_exact:
                print(f"      false-positive cost at q_cut={best} on "
                      f"{len(base_exact)} baseline-exact mates:")
                print(f"        -> undertrim (artefact left): {u_reads} "
                      f"({u_reads/len(base_exact)*100:.1f}%), "
                      f"{u_bp/max(1,u_reads):.2f} bp each")
                print(f"        -> overtrim (biology cut):    {o_reads} "
                      f"({o_reads/len(base_exact)*100:.1f}%), "
                      f"{o_bp/max(1,o_reads):.2f} bp each")
    return results


# ------------------------------------------------- part 4: cost
def part_cost(work, verbose):
    print()
    print("=== PART 4: throughput before/after the rule (100k pairs) ===")
    d = os.path.join(work, "cost")
    os.makedirs(d, exist_ok=True)
    fq1, fq2 = os.path.join(d, "big_R1.fastq.gz"), os.path.join(d, "big_R2.fastq.gz")
    if not (os.path.exists(fq1) and os.path.exists(fq2)):
        rnd = random.Random(99)
        NO_G = "ACGT"
        for path, mate in ((fq1, 1), (fq2, 2)):
            reads = []
            for i in range(100_000):
                body = [rnd.choice(NO_G) for _ in range(150)]
                quals = [37] * 150
                # sprinkle end runs at low quality: ~25% leading G (1-2 bp),
                # ~10% trailing G (3-12 bp). Low quality makes the rule walk
                # the whole run, which is its maximal-cost path.
                if rnd.random() < 0.25:
                    for j in range(rnd.randint(1, 2)):
                        body[j] = "G"
                        quals[j] = 10
                if rnd.random() < 0.10:
                    for j in range(150 - rnd.randint(3, 12), 150):
                        body[j] = "G"
                        quals[j] = 10
                seq = "".join(body)
                reads.append((f"@C{i}:{i%10000}:{i//10000} {mate}:N:0:ATCG",
                              seq, quals))
            write_fastq(path, reads)

    def time_run(tag, extra, reps=3):
        best = None
        for rep in range(reps):
            out = os.path.join(d, f"t_{tag}_{rep}")
            t0 = time.time()
            run_fastprep(os.path.join(d, "big"), out, extra, paired=True)
            dt = time.time() - t0
            best = dt if best is None else min(best, dt)
        rate = 200_000 / best
        print(f"  {tag:24s} best {best:.2f}s  ->  {rate:,.0f} mates/s")
        return rate

    base = time_run("fastprep (rule off)", BASE_CFG)
    has_flag = "--q-boundary" in subprocess.run(
        [FASTPREP, "--help"], capture_output=True, text=True).stdout
    if has_flag:
        on = time_run("fastprep (--q-boundary 25)",
                      BASE_CFG + ["--q-boundary", "25"])
        print(f"  throughput delta: {(on/base - 1)*100:+.1f}%")
    else:
        print("  binary lacks --q-boundary; rule cost measured analytically")
    return base


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="all",
                    choices=["model", "probe", "rule", "cost", "all"])
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    if args.part in ("model", "all"):
        part_model(WORK, args.verbose)
    if args.part in ("probe", "all"):
        part_probe(WORK, args.verbose)
    if args.part in ("rule", "all"):
        part_rule(WORK, args.verbose)
    if args.part in ("cost", "all"):
        part_cost(WORK, args.verbose)
    return 0


if __name__ == "__main__":
    sys.exit(main())
