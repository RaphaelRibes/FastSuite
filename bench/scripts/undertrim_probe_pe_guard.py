#!/usr/bin/env python3
"""The decisive UT-1 guard A/B: does the PE detector's own floor discard the
boundary it just pinned?

WHY THIS PROBE EXISTS
---------------------
`detect_pe_overlap` pins the insert boundary from mate geometry alone. The
boundary is then handed to `readthrough_is_adapter`, which asks the general
matcher whether the tail past that boundary is adapter. The general matcher
applies `min_overlap_se` (14 by default) to the tail it is given -- so a
10-base tail is PROVABLY unmatchable there, and the detector's perfect
boundary is thrown away.

The probe makes that claim falsifiable:

  R1 = fragment + TruSeqR1[:t]     (real adapter readthrough)
  R2 = revcomp(fragment) + a RANDOM t-base tail (never adapter)

R2 carries NO adapter at all. Any exact cut on R2 can therefore only come
from the PE-overlap path -- R2 is a pure probe of detector + guard, with the
ordinary matcher unable to contribute.

Arms (t, --min-adapter-overlap):
    A  t=20, floor 14  -> both mates cut exactly   (the guard passes)
    B  t=10, floor 14  -> 0/200                    (the floor blocks the guard)
    C  t=10, floor 10  -> both mates cut exactly   (same geometry, lower floor)

Arm B vs arm C is the whole experiment: identical geometry, one flag apart.

Usage
-----
    python3 bench/scripts/undertrim_probe_pe_guard.py [--n 200] [--out JSON]
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from undertrim_probe import TRUSEQ_R1, dna, parse_pxmt  # noqa: E402

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def write_pair(d, pairs):
    """Write R1/R2 FASTQ files. `pairs` is [(id, r1_seq, r2_seq)]."""
    with open(os.path.join(d, "r1.fq"), "w") as f1, \
         open(os.path.join(d, "r2.fq"), "w") as f2:
        for rid, s1, s2 in pairs:
            f1.write(f"@{rid}\n{s1}\n+\n{'I' * len(s1)}\n")
            f2.write(f"@{rid}\n{s2}\n+\n{'I' * len(s2)}\n")


def run_paired(fp, d, floor, mode="always"):
    """Run fastprep on the pair; return {(id, mate): (t5_adapter, t3_adapter)}."""
    o1 = os.path.join(d, f"o1_{floor}.fq")
    o2 = os.path.join(d, f"o2_{floor}.fq")
    px = os.path.join(d, f"o_{floor}.pxmt")
    for p in (o1, o2, px):
        if os.path.exists(p):
            os.remove(p)
    cmd = [fp, "-1", os.path.join(d, "r1.fq"), "-2", os.path.join(d, "r2.fq"),
           "-o", o1, "-p", o2, "--metadata", px,
           "--adapter-preset", "illumina", "--pe-overlap", mode,
           "--no-dedup", "--no-mean-q", "--min-len", "0", "-j", "1", "-q",
           "--min-adapter-overlap", str(floor)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"fastprep failed (floor {floor}): {r.stderr[-1500:]}")
    out = parse_pxmt(px)
    return {(k[0], k[1]): (v[1], v[2]) for k, v in out.items()}


def build_arm(rng, n, L, t):
    """Build n mate pairs with insert I = L - t.

    R1 = fragment + TruSeqR1[:t]        (the readthrough is REAL adapter)
    R2 = revcomp(fragment) + RANDOM tail (R2 carries NO adapter at all)

    The geometry makes the two mates share exactly the fragment:
    `revcomp(R2) = revcomp(random_tail) + fragment`, so R1's prefix of
    length I equals the suffix of `revcomp(R2)` of length I -- the detector
    pins the boundary at I = L - t, and t3_1 = L - I = t.

    Why R2's tail is random: `detect_pe_overlap_trim` checks the guard on
    R1 (the mate whose boundary the detector pinned) and then applies the
    boundary to BOTH mates. So R2's own tail never has to be adapter. R2
    carries none, which makes an exact cut on R2 impossible for the ordinary
    matcher to produce: it is evidence of the PE path and nothing else.
    """
    pairs, truth = [], {}
    for k in range(n):
        rid = f"G{t}_{k}"
        frag = dna(rng, L - t)
        r1 = frag + TRUSEQ_R1[:t]
        r2 = revcomp(frag) + dna(rng, t)
        pairs.append((rid, r1, r2))
        truth[rid] = t
    return pairs, truth


def score(cuts, truth, mate):
    """exact = cut == t; none = cut == 0; under = 0 < cut < t; over = cut > t."""
    n = len(truth)
    exact = none = under = over = 0
    for rid, t in truth.items():
        t5, t3 = cuts.get((rid, mate), (0, 0))
        cut = t5 + t3
        if cut == t:
            exact += 1
        elif cut == 0:
            none += 1
        elif cut < t:
            under += 1
        else:
            over += 1
    return dict(n=n, exact=exact, no_trim=none, partial=under, past_truth=over)


def arm(fp, rng, n, L, t, floor, label):
    d = tempfile.mkdtemp(prefix="utg_")
    pairs, truth = build_arm(rng, n, L, t)
    write_pair(d, pairs)
    cuts = run_paired(fp, d, floor)
    r1 = score(cuts, truth, 0)
    r2 = score(cuts, truth, 1)
    return dict(arm=label, t=t, floor=floor, r1=r1, r2=r2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=31)
    ap.add_argument("--L", type=int, default=150)
    ap.add_argument("--out", default="bench/results/undertrim-probe-pe-guard.json")
    a = ap.parse_args()
    fp = os.path.abspath(a.fastprep)
    rng = random.Random(a.seed)

    print(f"== UT-1 guard probe: n={a.n} per arm, L={a.L}")
    print(f"{'arm':>4} {'t':>3} {'floor':>5} {'R1 exact':>9} {'R2 exact':>9} "
          f"{'R2 none':>8} {'R2 partial':>10}")

    arms = [
        ("A", 20, 14, "A_t20_mo14_PE_fires_when_guard_passes"),
        ("B", 10, 14, "B_t10_mo14_guard_floor_blocks"),
        ("C", 10, 10, "C_t10_mo10_floor_lowered_guard_passes"),
    ]
    report: dict = {
        "purpose": (
            "Decisive A/B: does detect_pe_overlap pin short-tail boundaries, "
            "or does readthrough_is_adapter's min_overlap_se floor block "
            "them? R2 carries a RANDOM tail, so any exact R2 cut is "
            "PE-path-only evidence."
        ),
        "config": {
            "read_len": a.L,
            "n_per_arm": a.n,
            "adapter": "TRUSEQ_R1",
            "pe_mode": "always",
            "set": "--adapter-preset illumina",
            "seed": a.seed,
        },
        "fastprep": fp,
        "arms": {},
    }
    for label, t, floor, name in arms:
        r = arm(fp, rng, a.n, a.L, t, floor, label)
        r["name"] = name
        report["arms"][name] = r
        print(f"{label:>4} {t:3d} {floor:5d} {r['r1']['exact']:9d} "
              f"{r['r2']['exact']:9d} {r['r2']['no_trim']:8d} "
              f"{r['r2']['partial']:10d}")

    # The verdict, stated so a reader does not have to infer it.
    A = report["arms"]["A_t20_mo14_PE_fires_when_guard_passes"]
    B = report["arms"]["B_t10_mo14_guard_floor_blocks"]
    C = report["arms"]["C_t10_mo10_floor_lowered_guard_passes"]
    guard_blocks = B["r2"]["exact"] == 0 and C["r2"]["exact"] > 0
    report["verdict"] = {
        "arm_A_both_mates_cut": A["r1"]["exact"] == a.n and A["r2"]["exact"] == a.n,
        "arm_B_nothing_cut": B["r2"]["exact"] == 0,
        "arm_C_both_mates_cut": C["r1"]["exact"] == a.n and C["r2"]["exact"] == a.n,
        "the_floor_is_the_blocker": guard_blocks,
    }
    print("\nverdict:")
    for k, v in report["verdict"].items():
        print(f"  {k}: {v}")

    with open(a.out, "w") as f:
        json.dump(report, f, indent=1)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
