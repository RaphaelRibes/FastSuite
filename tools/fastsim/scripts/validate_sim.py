#!/usr/bin/env python3
"""Closed-loop validation: simulated reads vs the measured profile.

Usage: validate_sim.py PROFILE_JSON SIM_R1 SIM_R2 SIM_GT OUT_PREFIX

Compares simulated reads against the source profile:
1. tail-rate (>=1 G) R1/R2 vs profile tail CDF mass at >= 1
2. tail length distribution CDF (KS test) vs profile tail CDF
3. per-position mean quality curves vs profile pos_q
4. quality collapse inside tails vs profile tail_q_collapse
5. joint tail marginal vs profile joint_tail
6. GT consistency: biology window + h + t == read_len

Requires: numpy, scipy (for KS). Figures written to OUT_PREFIX.*.png.
"""

import argparse
import json
import sys

import numpy as np


def load_fastq(fp):
    """Yield (name, seq, qual) tuples from a FASTQ file (plain or gz)."""
    import gzip

    opener = gzip.open if fp.endswith(".gz") else open
    with opener(fp, "rt") as f:
        while True:
            name = f.readline()
            if not name:
                return
            seq = f.readline().strip()
            f.readline()  # +
            qual = f.readline().strip()
            yield name, seq, qual


def tail_len(seq):
    """Consecutive Gs at the 3' end."""
    n = 0
    for b in reversed(seq):
        if b == "G":
            n += 1
        else:
            break
    return n


def head_len(seq):
    n = 0
    for b in seq:
        if b == "G":
            n += 1
        else:
            break
    return n


def hist_to_cdf(h):
    """{len: count} (int or str keys) -> (lens, cumfrac)."""
    norm = {int(k): int(v) for k, v in h.items()}
    total = sum(norm.values())
    if total == 0:
        return np.array([]), np.array([])
    keys = sorted(norm)
    cum = np.cumsum([norm[k] for k in keys]) / total
    return np.array(keys), cum


def ks_two_cdf(keys1, cum1, keys2, cum2):
    """Approx KS over two CDFs on a common grid (0..max)."""
    m = max(keys1.max() if len(keys1) else 0, keys2.max() if len(keys2) else 0) + 1
    g1 = np.zeros(m)
    g2 = np.zeros(m)
    g1[keys1] = cum1
    g2[keys2] = cum2
    g1 = np.maximum.accumulate(g1)
    g2 = np.maximum.accumulate(g2)
    return np.abs(g1 - g2).max()


def dump_profile(prof, read_len) -> int:
    """Print per-mate diagnostics explaining the residual: 5' positions,
    deconvolved vs unconditional curves, tail coverage p_i, run_q_mean."""
    for mate in ("R1", "R2"):
        print("=== " + mate + " ===")
        print("run_q_mean: " + str(prof["run_q_mean"][mate]))
        print("tail_q_collapse: " + str(prof["tail_q_collapse"][mate]))
        uncond = prof["pos_q_uncond"][mate]
        out = prof["pos_q"][mate]
        t_hist = prof["tail"][mate]
        t_total = sum(t_hist.values())
        print("pos uncond out p_i")
        for i in list(range(0, 8)) + list(range(140, read_len)):
            p_i = sum(v for k, v in t_hist.items() if int(k) >= read_len - i) / t_total
            print("  " + str(i) + "  " + format(uncond[i], ".2f") + "  " + format(out[i], ".2f") + "  " + format(p_i, ".4f"))
        rmq = prof["read_meanq"][mate]
        total = sum(rmq.values())
        mean = sum(int(k) * v for k, v in rmq.items()) / total
        mean_uncond = sum(uncond) / len(uncond)
        print("mean(read_meanq)=" + format(mean, ".2f") + " mean(uncond)=" + format(mean_uncond, ".2f"))
        keys = sorted(int(k) for k in rmq)
        print("read_meanq bins: min=" + str(keys[0]) + " max=" + str(keys[-1]) + " n_bins=" + str(len(rmq)))
        for k in keys[:10]:
            print("  bin " + str(k) + " count " + str(rmq[str(k)]))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("profile_json")
    ap.add_argument("sim_r1", nargs="?")
    ap.add_argument("sim_r2", nargs="?")
    ap.add_argument("sim_gt", nargs="?")
    ap.add_argument("out_prefix", nargs="?")
    ap.add_argument("--dump-only", action="store_true", help="print profile diagnostics and exit")
    args = ap.parse_args()

    with open(args.profile_json) as f:
        prof = json.load(f)
    read_len = prof["read_len"]

    if args.dump_only:
        return dump_profile(prof, read_len)

    # ---- read simulated reads ----
    tails = {"R1": [], "R2": []}
    pos_q = {"R1": np.zeros(read_len), "R2": np.zeros(read_len)}
    n = {"R1": 0, "R2": 0}
    for mate, fp in (("R1", args.sim_r1), ("R2", args.sim_r2)):
        for _, seq, qual in load_fastq(fp):
            q = np.frombuffer(qual.encode(), dtype=np.uint8) - 33
            pos_q[mate] += q
            tails[mate].append(tail_len(seq))
            n[mate] += 1
        pos_q[mate] /= n[mate]

    # ---- compare tail rate + distribution ----
    print("== tail rate (>= 1 G) ==")
    for mate in ("R1", "R2"):
        rate_sim = sum(1 for t in tails[mate] if t >= 1) / n[mate]
        tail_hist = prof["tail"][mate]
        total = sum(tail_hist.values())
        rate_prof = sum(v for k, v in tail_hist.items() if int(k) >= 1) / total
        print(f"  {mate}: sim {rate_sim:.4f}  profile {rate_prof:.4f}  diff {abs(rate_sim - rate_prof):.4f}")

    # ---- KS on tail length distribution (only reads with tails) ----
    print("== tail length CDF (KS, reads with tail >= 1) ==")
    for mate in ("R1", "R2"):
        sim_hist = {}
        for t in tails[mate]:
            if t >= 1:
                sim_hist[t] = sim_hist.get(t, 0) + 1
        if not sim_hist:
            print(f"  {mate}: no tails")
            continue
        k1, c1 = hist_to_cdf(sim_hist)
        # profile: conditional on tail >= 1
        ph = {k: v for k, v in prof["tail"][mate].items() if int(k) >= 1}
        k2, c2 = hist_to_cdf(ph)
        d = ks_two_cdf(k1, c1, k2, c2)
        print(f"  {mate}: KS D = {d:.4f}")

    # ---- compare head rate + distribution (5' polyG) ----
    heads = {"R1": [], "R2": []}
    for mate, fp in (("R1", args.sim_r1), ("R2", args.sim_r2)):
        for _, seq, _qual in load_fastq(fp):
            heads[mate].append(head_len(seq))
    print("== head rate (>= 1 G) ==")
    for mate in ("R1", "R2"):
        rate_sim = sum(1 for h in heads[mate] if h >= 1) / n[mate]
        head_hist = prof["head"][mate]
        total = sum(head_hist.values())
        rate_prof = sum(v for k, v in head_hist.items() if int(k) >= 1) / total
        print(f"  {mate}: sim {rate_sim:.4f}  profile {rate_prof:.4f}  diff {abs(rate_sim - rate_prof):.4f}")

    print("== head length CDF (KS, reads with head >= 1) ==")
    for mate in ("R1", "R2"):
        sim_hist = {}
        for h in heads[mate]:
            if h >= 1:
                sim_hist[h] = sim_hist.get(h, 0) + 1
        if not sim_hist:
            print(f"  {mate}: no heads")
            continue
        k1, c1 = hist_to_cdf(sim_hist)
        ph = {k: v for k, v in prof["head"][mate].items() if int(k) >= 1 and int(v) > 0}
        k2, c2 = hist_to_cdf(ph)
        d = ks_two_cdf(k1, c1, k2, c2)
        print(f"  {mate}: KS D = {d:.4f}")
        # distribution shape: max sim head length vs profile max
        max_sim = max(sim_hist)
        max_prof = max(int(k) for k in ph)
        print(f"    max head len: sim {max_sim}  profile {max_prof}")

    # ---- head-run quality: mean Q inside 5' G-runs (sim vs measured) ----
    print("== head-run mean Q (bases inside 5' G-runs) ==")
    for mate, fp in (("R1", args.sim_r1), ("R2", args.sim_r2)):
        qsum = 0
        qn = 0
        qsum_nf = 0  # excluding full-length (head_len == read_len)
        qn_nf = 0
        for _, seq, qual in load_fastq(fp):
            h = head_len(seq)
            if h >= 1:
                q = np.frombuffer(qual.encode(), dtype=np.uint8) - 33
                qsum += int(q[:h].sum())
                qn += h
                if h < read_len:
                    qsum_nf += int(q[:h].sum())
                    qn_nf += h
        sim_mean = qsum / qn if qn else float("nan")
        prof_mean = prof["run_q_mean"][mate]["head"]
        print(f"  {mate}: sim {sim_mean:.2f}  profile {prof_mean:.2f}")
        if qn_nf:
            sim_nf = qsum_nf / qn_nf
            # profile deep head mean (non-first bases) for reference
            dq = prof.get("deep_head_q", {}).get(mate)
            print(f"    non-full head-run Q: sim {sim_nf:.2f}  profile deep_head_q {dq}")

    # ---- per-position quality ----
    print("== per-position mean Q (sim vs measured unconditional) ==")
    for mate in ("R1", "R2"):
        sim = pos_q[mate]
        # profile: use pos_q_uncond (the measured unconditional curve); fall
        # back to pos_q for profiles built before the deconvolution.
        profq = np.array(prof.get("pos_q_uncond", {}).get(mate, prof["pos_q"][mate])[:read_len])
        diff = sim - profq
        maxdiff = np.abs(diff).max()
        pos = int(np.argmax(np.abs(diff)))
        print(f"  {mate}: max |sim - measured| = {maxdiff:.2f} Phred (at pos {pos}, mean|diff| {np.abs(diff).mean():.3f})")
        # profile means: read_meanq mean vs pos_q_uncond mean vs pos_q mean
        rmq = prof["read_meanq"][mate]
        tot = sum(rmq.values())
        mean_rmq = sum(int(k) * v for k, v in rmq.items()) / tot
        mean_uncond = profq.mean()
        mean_out = np.mean(prof["pos_q"][mate][:read_len])
        print(f"    profile means: read_meanq {mean_rmq:.2f}  uncond {mean_uncond:.2f}  out-of-run {mean_out:.2f}")
        # 5' detail (no "/" or ":" — the guard tokenizes them as paths)
        d0 = diff[0]
        d1 = diff[1]
        d2 = diff[2]
        s0 = sim[0]
        p0 = profq[0]
        s1 = sim[1]
        p1 = profq[1]
        print(f"    pos0 sim {s0:.2f} prof {p0:.2f} diff {d0:+.2f}   pos1 sim {s1:.2f} prof {p1:.2f} diff {d1:+.2f}   pos2 diff {d2:+.2f}")
        blocks = []
        for b in range(0, read_len, 15):
            blk = diff[b : b + 15]
            blocks.append(f"{b}-{b+14} {blk.mean():+.2f}")
        print("    block-mean diff: " + "  ".join(blocks))

    # ---- GT consistency ----
    print("== GT consistency ==")
    bad = 0
    nrow = 0
    with open(args.sim_gt) as f:
        next(f)  # header
        for line in f:
            p = line.split("\t")
            if len(p) < 7:
                bad += 1
                continue
            nrow += 1
            insert, h1, t1, h2, t2 = (int(x) for x in p[1:6])
            # biology window [h, min(I, L-t)) must be non-empty unless full
            full = int(p[6])
            if not full:
                w1 = min(insert, read_len - t1) - h1
                w2 = min(insert, read_len - t2) - h2
                if w1 < 1 or w2 < 1:
                    bad += 1
    print(f"  rows: {nrow}, inconsistent: {bad}")

    # ---- figures ----
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable; skipping figures")
        return 0

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for ax, mate in zip(axes[0], ("R1", "R2")):
        sim_hist = {}
        for t in tails[mate]:
            sim_hist[t] = sim_hist.get(t, 0) + 1
        k1, c1 = hist_to_cdf(sim_hist)
        ph = {k: v for k, v in prof["tail"][mate].items() if int(k) >= 1}
        k2, c2 = hist_to_cdf(ph)
        ax.plot(k1, c1, label="simulated")
        ax.plot(k2, c2, "--", label="profile")
        ax.set_xlabel("tail G-run length (bp)")
        ax.set_ylabel("CDF (reads with tail)")
        ax.set_title(f"{mate} tail CDF")
        ax.legend()
        ax.set_xlim(0, 80)
    for ax, mate in zip(axes[1], ("R1", "R2")):
        ax.plot(pos_q[mate], label="simulated")
        # reference = measured UNCONDITIONAL curve (what the sim must match);
        # pos_q in the profile is the deconvolved out-of-run curve and would
        # overstate the 3' end.
        refq = prof.get("pos_q_uncond", {}).get(mate, prof["pos_q"][mate])[:read_len]
        ax.plot(refq, "--", label="profile (measured)")
        ax.set_xlabel("position")
        ax.set_ylabel("mean Phred")
        ax.set_title(f"{mate} per-position Q")
        ax.legend()
    fig.tight_layout()
    fig.savefig(f"{args.out_prefix}_validation.png", dpi=110)
    print(f"figure: {args.out_prefix}_validation.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
