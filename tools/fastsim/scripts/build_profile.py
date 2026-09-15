#!/usr/bin/env python3
"""Build a simulation profile JSON from polyg-profiler output.

Aggregates every stats/E*.json (one per pair) into a single profile:
- head/tail G-run length histograms (R1/R2 separate), exact counts
- joint R1×R2 tail histogram (capped at 64 like the profiler)
- per-position mean quality curves (R1/R2)
- per-read mean-Q histograms (for per-read baseline variance)
- tail quality collapse (overall mean Q minus tail-run mean Q)
- full-length spike rate + quality (from tail_len_hist at read_len)
- GC fraction

Usage: build_profile.py STATS_DIR OUT_JSON [--read-len 150]
"""

import argparse
import glob
import json
import os
import sys

JOINT_CAP = 64  # must match profile::JOINT_CAP


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("stats_dir", help="directory of stats/E*.json")
    ap.add_argument("out", help="output profile JSON path")
    ap.add_argument("--read-len", type=int, default=150)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.stats_dir, "*.json")))
    if not files:
        print(f"no JSONs found in {args.stats_dir}", file=sys.stderr)
        return 1
    print(f"aggregating {len(files)} files")

    # Accumulators (counts per mate).
    head = {"R1": {}, "R2": {}}
    tail = {"R1": {}, "R2": {}}
    joint = {}
    qsum = {"R1": [0] * args.read_len, "R2": [0] * args.read_len}
    n_reads = {"R1": 0, "R2": 0}
    read_meanq = {"R1": {}, "R2": {}}
    # run quality totals
    run_q = {"R1": {"head_sum": 0, "head_n": 0, "tail_sum": 0, "tail_n": 0},
             "R2": {"head_sum": 0, "head_n": 0, "tail_sum": 0, "tail_n": 0}}
    # per-position per-base counts [A,C,G,T,N]
    comp = {"R1": [[0] * 5 for _ in range(args.read_len)],
            "R2": [[0] * 5 for _ in range(args.read_len)]}

    def add_hist(dst, src):
        for k, v in src.items():
            dst[k] = dst.get(k, 0) + v

    for fp in files:
        with open(fp) as f:
            v = json.load(f)
        L = v.get("read_len", args.read_len)
        if L != args.read_len:
            print(f"warning: {fp} read_len {L} != {args.read_len}, skipping")
            continue
        for mate in ("R1", "R2"):
            s = v[f"{mate.lower()}_stats"]
            add_hist(head[mate], s["head_len_hist"])
            add_hist(tail[mate], s["tail_len_hist"])
            add_hist(read_meanq[mate], s["read_meanq_hist"])
            n_reads[mate] += s["n_reads"]
            for i, q in enumerate(s["qsum"][:args.read_len]):
                qsum[mate][i] += q
            for i, c in enumerate(s["comp"][:args.read_len]):
                for b in range(5):
                    comp[mate][i][b] += c[b]
            rq = s["run_q"]
            for k in ("head_sum", "head_n", "tail_sum", "tail_n"):
                run_q[mate][k] += rq[k]
        # joint (already capped by profiler; enforce cap defensively)
        for k, cnt in v["joint_tail"].items():
            a, b = k.split(",")
            a, b = min(int(a), JOINT_CAP), min(int(b), JOINT_CAP)
            key = f"{a},{b}"
            joint[key] = joint.get(key, 0) + cnt

    # ---- derive profile fields ----
    # Measured per-position means (unconditional: include tail-run collapse).
    pos_q_uncond = {}
    for mate in ("R1", "R2"):
        n = n_reads[mate]
        pos_q_uncond[mate] = [round(q / n, 3) for q in qsum[mate]] if n else [0] * args.read_len

    # Run-quality means (head vs tail) from run_q totals.
    # NOTE: run_q.head includes the FIRST base of each head run, which is
    # normal quality (a 1bp head is just a biological G). The simulator
    # models that first base with the out-of-run curve and collapses only
    # DEEPER head bases, so it needs the deep head mean: exclude the first
    # base of every head run from the average.
    run_q_mean = {}
    deep_head_q = {}
    for mate in ("R1", "R2"):
        head_n = run_q[mate]["head_n"]
        head_sum = run_q[mate]["head_sum"]
        all_head = round(head_sum / head_n, 3) if head_n else None
        run_q_mean[mate] = {
            "head": all_head,
            "tail": round(run_q[mate]["tail_sum"] / run_q[mate]["tail_n"], 3)
            if run_q[mate]["tail_n"]
            else None,
        }
        # deep head mean: remove one first base per NON-full head read (the
        # first base is normal Q) AND exclude full-length reads entirely —
        # those 150 all-G bases are handled by the full-length path with
        # full_length_q, and including them here would bias the deep mean
        # (they are not "first base + deep bases").
        h_hist = head[mate]
        full_count = h_hist.get(str(args.read_len), 0)
        non_full_heads = sum(v for k, v in h_hist.items() if 1 <= int(k) < args.read_len)
        n_head_reads = non_full_heads + full_count
        head_bases = sum(int(k) * v for k, v in h_hist.items())
        # deep bases = all non-full head bases minus one first base per head
        deep_n = head_bases - full_count * args.read_len - non_full_heads
        if deep_n > 0 and all_head is not None:
            # full-length quality proxy = tail-run mean (as full_q below)
            full_q_m = (
                run_q[mate]["tail_sum"] / run_q[mate]["tail_n"]
                if run_q[mate]["tail_n"]
                else 20.0
            )
            head_sum_nf = head_sum - full_count * args.read_len * full_q_m
            first_q = pos_q_uncond[mate][0] if n_reads[mate] else all_head
            deep_head_q[mate] = round((head_sum_nf - non_full_heads * first_q) / deep_n, 3)
        else:
            deep_head_q[mate] = all_head

    # Deconvolve pos_q into the OUT-OF-RUN position curve.
    # measured[i] = (1-p_i-q_i)*out[i] + p_i*tail_q + q_i*head_q
    #   p_i = frac of reads whose TAIL covers position i (from tail CDF)
    #   q_i = frac of reads whose HEAD covers position i AND is NOT the first
    #         head base. The FIRST base of a head run is normal quality
    #         (measured: R1 heads are 1bp and Q 38.3 == outside; R2 1bp heads
    #         are likewise normal) — only DEEPER head bases collapse. So
    #         position 0 of a head is out-of-run, positions 1..h-1 collapse.
    # => out[i] = (measured[i] - p_i*tail_q - q_i*head_q) / (1 - p_i - q_i)
    pos_q = {}
    for mate in ("R1", "R2"):
        n = n_reads[mate]
        if n == 0:
            pos_q[mate] = [0] * args.read_len
            continue
        head_q = deep_head_q[mate]
        tail_q = run_q_mean[mate]["tail"]
        t_hist = tail[mate]
        t_total = sum(t_hist.values())
        h_hist = head[mate]
        h_total = sum(h_hist.values())
        out = []
        for i in range(args.read_len):
            # tail covers position i (0-based) iff t >= read_len - i
            p_i = sum(v for k, v in t_hist.items() if int(k) >= args.read_len - i) / t_total if t_total else 0.0
            # DEEP head coverage at position i: position 0 of any head is the
            # FIRST base (normal quality, excluded); positions 1..h-1 are
            # deep (collapsed). So q_i = P(h >= i+1) for i>=1, 0 at i=0.
            if i == 0:
                q_i = 0.0
            else:
                q_i = sum(v for k, v in h_hist.items() if int(k) >= i + 1) / h_total if h_total else 0.0
            meas = pos_q_uncond[mate][i]
            denom = 1 - p_i - q_i
            if denom <= 1e-6:
                out.append(meas)  # degenerate: everything covered here
            else:
                num = meas - p_i * (tail_q or 0.0) - q_i * (head_q or 0.0)
                out.append(round(num / denom, 3))
        pos_q[mate] = out

    collapse = {}
    for mate in ("R1", "R2"):
        n = n_reads[mate]
        overall = sum(qsum[mate]) / (n * args.read_len) if n else 0.0
        tq = run_q[mate]["tail_sum"] / run_q[mate]["tail_n"] if run_q[mate]["tail_n"] else overall
        collapse[mate] = round(max(0.0, overall - tq), 3)

    # full-length spike: reads whose tail covers the whole read (t == read_len)
    full_rate = {}
    full_q = {}
    for mate in ("R1", "R2"):
        t = tail[mate]
        n = n_reads[mate]
        full_count = t.get(str(args.read_len), 0)
        full_rate[mate] = full_count / n if n else 0.0
        # quality proxy: tail-run mean. Cross-checked against the tail
        # aggregate (full-length = ~9% of R2 tail bases): the measured tail
        # mean 37.19 is consistent with Q_full ~36.5 + non-full tails ~37.26.
        full_q[mate] = round(run_q[mate]["tail_sum"] / run_q[mate]["tail_n"], 3) if run_q[mate]["tail_n"] else 20

    # GC fraction (G+C over all positions, R1+R2 pooled)
    gc_num = gc_den = 0
    for mate in ("R1", "R2"):
        for pos in range(args.read_len):
            c = comp[mate][pos]
            gc_num += c[1] + c[2]  # C + G
            gc_den += sum(c)
    gc = round(gc_num / gc_den, 4) if gc_den else 0.5

    # spike rate = R2's (R1 has none); keep per-mate + a combined rate for the
    # full_length_rate field (decide() samples one rate for the pair).
    combined_full_rate = max(full_rate["R1"], full_rate["R2"])

    profile = {
        "read_len": args.read_len,
        "gc": gc,
        "geometry": "tail3p_both",
        "head": head,
        "tail": tail,
        "joint_tail": joint,
        "tail_q_collapse": collapse,
        "pos_q": pos_q,
        "pos_q_uncond": pos_q_uncond,
        "run_q_mean": run_q_mean,
        "deep_head_q": deep_head_q,
        "read_meanq": read_meanq,
        "full_length_rate": combined_full_rate,
        "full_length_q": {"R1": full_q["R1"], "R2": full_q["R2"]},
        "adapter": {"seq": "CTGTCTCTTATACACATCT", "mode": "readthrough"},
    }

    with open(args.out, "w") as f:
        json.dump(profile, f, indent=1, sort_keys=True)
    print(f"wrote {args.out}")
    print(f"  R1 reads: {n_reads['R1']:,}  R2 reads: {n_reads['R2']:,}")
    print(f"  GC: {gc:.3f}")
    print(f"  collapse: R1 {collapse['R1']}  R2 {collapse['R2']}")
    print(f"  full-length rate: R1 {full_rate['R1']:.5f}  R2 {full_rate['R2']:.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
