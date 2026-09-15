#!/usr/bin/env python3
"""Diagnose per-position Q diff: where is the max? Tail/head run fractions."""
import json

sim = json.load(open("/scratch/users/ribesr/novaseq-sim/validation/sim_profiled.json"))
prof = json.load(open("/scratch/users/ribesr/novaseq-sim/profiles/mandril_E.json"))

for mate in ("R1", "R2"):
    s = sim[f"{mate.lower()}_stats"]
    n = s["n_reads"]
    sim_q = [q / n for q in s["qsum"]]
    prof_q = prof["pos_q"][mate]
    diffs = [abs(sim_q[i] - prof_q[i]) for i in range(150)]
    mx = max(diffs)
    pos = diffs.index(mx)
    print(f"{mate}: max diff {mx:.2f} at position {pos} (3' end is 149)")
    print(f"   sim_q[145..149] = {[round(sim_q[i],1) for i in range(145,150)]}")
    print(f"   prof_q[145..149] = {[round(prof_q[i],1) for i in range(145,150)]}")
    # run fractions from the tail CDF
    t = prof["tail"][mate]
    tot = sum(t.values())
    cdf = sorted((int(k), v / tot) for k, v in t.items())
    # P(tail >= 151-i) at position i (0-based): tail covers [L-t, L)
    run_frac = []
    for i in range(150):
        # position i in tail iff t >= 150-i
        f = sum(v for k, v in cdf if k >= 150 - i)
        run_frac.append(f)
    print(f"   run_frac[145..149] = {[round(run_frac[i],3) for i in range(145,150)]}")
    # expected unconditional if sim used out-of-run curve + collapse:
    # out-of-run ~38.2 flat; expected_uncond[i] = out - run_frac[i]*collapse
    print()
