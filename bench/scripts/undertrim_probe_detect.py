#!/usr/bin/env python3
"""Third undertrim probe: the discovery path (what happens with no explicit
--adapter-preset), the 5' adapter case, and readthrough that runs past the
adapter into the downstream primer.

Usage: python3 bench/scripts/undertrim_probe_detect.py [--n 4000]
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from undertrim_probe import TRUSEQ_R1, dna, parse_pxmt, write_fastq  # noqa: E402

CUSTOM = "AACGTTGCAACGTTGCAACGT"          # 21 bp, in no catalogue
PRIMER = "AATGATACGGCGACCACCGAGATCTACAC"   # P5 flowcell oligo, 29 bp


def run(fp, d, name, extra):
    out1 = os.path.join(d, f"{name}.fq")
    px = os.path.join(d, f"{name}.pxmt")
    cmd = [fp, "-1", os.path.join(d, "in.fq"), "-o", out1, "--metadata", px,
           "--no-dedup", "--no-mean-q", "--min-len", "0", "-j", "1"]
    if extra is not None:
        cmd += extra
    else:
        cmd += ["-q"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"fastprep failed ({name}): {r.stderr[-1200:]}")
    return parse_pxmt(px), r.stderr


def score(px, ids, truth):
    ex = none = part = 0
    for i in ids:
        st, t5a, t3a, t5p, t3p, _ = px[(i, 0)]
        cut = t3a + t3p
        truth_i = truth[i]
        if cut == truth_i:
            ex += 1
        elif cut == 0:
            none += 1
        else:
            part += 1
    return dict(exact=ex, none=none, partial=part)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--n", type=int, default=4000)
    ap.add_argument("--out", default="bench/results/undertrim-probe-detect.json")
    a = ap.parse_args()
    fp = os.path.abspath(a.fastprep)
    rng = random.Random(23)
    rep = {}
    d = tempfile.mkdtemp(prefix="utd_")

    # ---- 1. off-catalogue adapter, auto-detection vs explicit FASTA ---------
    recs, ids, truth = [], [], {}
    for k in range(a.n):
        seq = dna(rng, 127) + CUSTOM
        rid = f"X{k}"
        ids.append(rid); truth[rid] = len(CUSTOM)
        recs.append((rid, seq))
    write_fastq(os.path.join(d, "in.fq"), recs)
    px, err = run(fp, d, "auto_offcat", None)
    rep["off_catalogue_auto"] = score(px, ids, truth)
    rep["off_catalogue_auto"]["stderr_head"] = err.strip().splitlines()[:2]
    fa = os.path.join(d, "custom.fa")
    open(fa, "w").write(f">custom\n{CUSTOM}\n")
    px, _ = run(fp, d, "fasta_offcat", ["--adapters", fa, "-q"])
    rep["off_catalogue_fasta"] = score(px, ids, truth)
    print("1. off-catalogue adapter (21 bp, 100% of reads carry it)")
    print("   auto-detect  :", rep["off_catalogue_auto"])
    print("   --adapters   :", rep["off_catalogue_fasta"])

    # ---- 2. discovery threshold: a mixed population -------------------------
    print("\n2. mixed population, f = share of reads with a full 33 bp readthrough")
    rep["detect_threshold"] = {}
    for f in (0.002, 0.005, 0.01, 0.02, 0.05):
        recs, ids, truth, carrier = [], [], {}, set()
        for k in range(a.n):
            rid = f"M{k}"
            if rng.random() < f:
                seq = dna(rng, 117) + TRUSEQ_R1
                truth[rid] = len(TRUSEQ_R1); carrier.add(rid)
            else:
                seq = dna(rng, 150); truth[rid] = 0
            ids.append(rid); recs.append((rid, seq))
        write_fastq(os.path.join(d, "in.fq"), recs)
        px, err = run(fp, d, f"mix{f}", None)
        car = [r for r in ids if r in carrier]
        non = [r for r in ids if r not in carrier]
        detected = "NO adapter trimming will be applied" not in err
        rep["detect_threshold"][f] = dict(
            detected=detected, carriers_trimmed=score(px, car, truth)["exact"],
            carriers=len(car), clean_trimmed=sum(
                1 for r in non
                if px[(r, 0)][2] + px[(r, 0)][4] > 0))
        r_ = rep["detect_threshold"][f]
        print(f"   f={f*100:4.1f}%  detected={str(detected):5s} "
              f"carriers exact {r_['carriers_trimmed']}/{r_['carriers']}  "
              f"clean falsely trimmed {r_['clean_trimmed']}")

    # ---- 3. readthrough past the adapter into the primer --------------------
    recs, ids, truth = [], [], {}
    for k in range(1000):
        rid = f"P{k}"
        seq = dna(rng, 150 - len(TRUSEQ_R1) - len(PRIMER)) + TRUSEQ_R1 + PRIMER
        ids.append(rid); truth[rid] = len(TRUSEQ_R1) + len(PRIMER)
        recs.append((rid, seq))
    write_fastq(os.path.join(d, "in.fq"), recs)
    px, _ = run(fp, d, "primer", ["--adapter-preset", "illumina", "-q"])
    rep["adapter_plus_primer"] = score(px, ids, truth)
    print("\n3. 33 bp adapter + primer readthrough:", rep["adapter_plus_primer"],
          f"(truth cut = {len(TRUSEQ_R1)+len(PRIMER)})")

    # ---- 4. 5' adapter (chimeric / small-RNA style) -------------------------
    recs, ids, truth = [], [], {}
    for k in range(1000):
        rid = f"F{k}"
        seq = TRUSEQ_R1 + dna(rng, 150 - 33)
        ids.append(rid); truth[rid] = 0
        recs.append((rid, seq))
    write_fastq(os.path.join(d, "in.fq"), recs)
    px, _ = run(fp, d, "fiveprime", ["--adapter-preset", "illumina", "-q"])
    t5 = sum(1 for r in ids if px[(r, 0)][1] > 0)
    t3 = sum(1 for r in ids if px[(r, 0)][2] > 0)
    rep["five_prime_adapter"] = dict(reads=len(ids), t5_trimmed=t5, t3_trimmed=t3,
        note="a read whose 5' end IS the adapter matches at p=0, so t3 = read_len and "
             "the record is emitted with a zero-length sequence")
    print("\n4. read starting with the adapter (5' contamination): "
          f"t5>0 on {t5}/{len(ids)}, t3>0 on {t3}/{len(ids)}")

    with open(a.out, "w") as f:
        json.dump(rep, f, indent=1)
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
