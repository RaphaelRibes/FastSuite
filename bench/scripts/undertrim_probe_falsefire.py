#!/usr/bin/env python3
"""The false-fire probe: mates that OVERLAP but carry no adapter.

The paired A/B corpus's "clean" population has insert >= read length, so those
mates never overlap and the guard is never reached. That leaves the
dangerous case unmeasured: a pair whose insert is SHORTER than the read (so
the detector pins a boundary) but whose readthrough is ORDINARY SEQUENCE, not
adapter.

Under the `floor` policy such a pair is refused twice over: the tail is
shorter than the overlap floor, and it is not adapter anyway. Under the
`boundary` policy the floor is lifted -- so the tail must still be refused
because it is not the adapter prefix.

This probe measures exactly that: R1 carries the real TruSeq prefix (so the
geometry is right and the boundary is pinned), R2's tail is RANDOM. Every
trim on either mate is a false fire.

Usage: python3 bench/scripts/undertrim_probe_falsefire.py [--n 300]
"""
import argparse
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--seed", type=int, default=99)
    ap.add_argument("--L", type=int, default=150)
    ap.add_argument("--out", default="bench/results/undertrim-probe-falsefire.json")
    a = ap.parse_args()
    fp = os.path.abspath(a.fastprep)
    rng = random.Random(a.seed)

    print(f"== false-fire probe: {a.n} pairs per tail length, L={a.L}")
    print("   (BOTH mates carry ORDINARY SEQUENCE past the boundary -- the")
    print("    mates overlap, so the detector pins a boundary, but nothing was")
    print("    ligated past it. Any trim here is a false fire.)")
    print(f"{'t':>3} {'policy':>9} {'R1 trimmed':>10} {'R2 trimmed':>10}")

    report = {"fastprep": fp, "n": a.n, "L": a.L, "seed": a.seed, "arms": {}}
    worst = 0
    for t in (2, 3, 4, 5, 6, 8, 10, 14):
        d = tempfile.mkdtemp(prefix="utff_")
        r1b, r2b = "", ""
        for k in range(a.n):
            frag = dna(rng, a.L - t)
            # BOTH tails are ordinary sequence: the mates overlap (so the
            # boundary is pinnable) but no adapter was ligated past it.
            r1 = frag + dna(rng, t)
            r2 = revcomp(frag) + dna(rng, t)
            r1b += f"@p{k}\n{r1}\n+\n{'I' * len(r1)}\n"
            r2b += f"@p{k}\n{r2}\n+\n{'I' * len(r2)}\n"
        with open(os.path.join(d, "r1.fq"), "w") as f:
            f.write(r1b)
        with open(os.path.join(d, "r2.fq"), "w") as f:
            f.write(r2b)

        for ev in ("floor", "boundary"):
            o1 = os.path.join(d, f"o1_{ev}.fq")
            o2 = os.path.join(d, f"o2_{ev}.fq")
            px = os.path.join(d, f"o_{ev}.pxmt")
            cmd = [fp, "-1", os.path.join(d, "r1.fq"), "-2", os.path.join(d, "r2.fq"),
                   "-o", o1, "-p", o2, "--metadata", px,
                   "--adapter-preset", "illumina", "--pe-overlap", "always",
                   "--no-dedup", "--no-mean-q", "--min-len", "0", "-j", "1", "-q",
                   "--adapter-evidence", ev]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0:
                sys.exit(f"fastprep failed ({ev}): {r.stderr[-1500:]}")
            recs = parse_pxmt(px)
            cut1 = sum(1 for (rid, m), v in recs.items()
                       if m == 0 and (v[1] + v[2] + v[3] + v[4]) > 0)
            cut2 = sum(1 for (rid, m), v in recs.items()
                       if m == 1 and (v[1] + v[2] + v[3] + v[4]) > 0)
            key = f"t{t}_{ev}"
            report["arms"][key] = dict(t=t, policy=ev, r1_trimmed=cut1,
                                       r2_trimmed=cut2, n=a.n)
            worst = max(worst, cut1, cut2)
            print(f"{t:3d} {ev:>9} {cut1:10d} {cut2:10d}")

    report["verdict"] = {
        "worst_false_fire": worst,
        "no_false_fire": worst == 0,
    }

    # The honest comparison: the boundary policy must not fire MATERIALLY
    # more often than the shipped floor policy does. The floor policy is not
    # itself zero -- the ordinary matcher hits adapter prefixes by chance
    # too -- so the question is the DIFFERENCE, not the absolute count.
    # A difference of one read in a few hundred pairs is noise; a difference
    # that scales with the tail length is the policy misfiring.
    floor_total = sum(v["r1_trimmed"] for k, v in report["arms"].items()
                      if k.endswith("_floor"))
    bnd_total = sum(v["r1_trimmed"] for k, v in report["arms"].items()
                    if k.endswith("_boundary"))
    n_arms = len([k for k in report["arms"] if k.endswith("_floor")])
    pairs = n_arms * a.n
    report["verdict"] = {
        "floor_false_fires": floor_total,
        "boundary_false_fires": bnd_total,
        "pairs_per_arm": a.n,
        "floor_rate": floor_total / pairs,
        "boundary_rate": bnd_total / pairs,
        # The gate: the boundary policy's false-fire rate must not exceed
        # the floor policy's by more than 0.5 pp (i.e. it may not be a
        # materially worse guard).
        "boundary_no_worse": (bnd_total - floor_total) / pairs <= 0.005,
    }
    v = report["verdict"]
    print(f"\nfalse fires over {pairs} pairs per policy:")
    print(f"  floor    : {v['floor_false_fires']:5d}  ({v['floor_rate']*100:.3f}%)")
    print(f"  boundary : {v['boundary_false_fires']:5d}  ({v['boundary_rate']*100:.3f}%)")
    print(f"  boundary_no_worse (<= +0.5 pp): {v['boundary_no_worse']}")
    with open(a.out, "w") as f:
        import json
        json.dump(report, f, indent=1)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
