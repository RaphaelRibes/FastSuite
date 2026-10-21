#!/usr/bin/env python3
"""Undertrim probe: what the shipped adapter matcher actually removes.

Pure black-box measurement through the fastprep binary and its PXMT metadata
(the same per-read trim record bench/scripts/overtrim_vs_fastp.py reads), so
nothing here re-implements the matcher.

Four arms
---------
A  tail ladder       read = biology(150-t) + TruSeqR1[0:t]      t = 1..40
B  error arm         same, with e substitutions inside the t adapter bases
C  contamination     same, with the last g bases replaced by G (NovaSeq
                     no-signal) or n bases replaced by N inside the adapter
D  clean-read arm    adapter-free random reads, swept over
                     --min-adapter-overlap to price the false-trim rate of a
                     lower floor (the overtrim side of any floor change)

Usage
    python3 bench/scripts/undertrim_probe.py [--fastprep PATH] [--n 2000]
Writes a JSON summary next to the stdout report.
"""
import argparse
import json
import os
import random
import struct
import subprocess
import sys
import tempfile

TRUSEQ_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRUSEQ_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


def dna(rng, n, gc=0.47):
    return "".join(
        rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(n)
    )


def write_fastq(path, records):
    with open(path, "w") as f:
        for rid, seq in records:
            f.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")


def parse_pxmt(path):
    """{(read_id, mate): (status, t5_adapter, t3_adapter, t5_poly, t3_poly, read_len)}

    Paired runs write both mates into ONE file, and a pair shares its read id,
    so the key is (id, mate). Single-end records carry mate 0.
    """
    out = {}
    with open(path, "rb") as f:
        blob = f.read()
    assert blob[:4] == b"PXMT", "not a PXMT file"
    i = 5
    while i < len(blob):
        (idlen,) = struct.unpack_from("<H", blob, i)
        i += 2
        rid = blob[i:i + idlen].decode()
        i += idlen
        mate, status = blob[i], blob[i + 1]
        i += 2
        t5a, t3a, t5p, t3p, rlen = struct.unpack_from("<IIIII", blob, i)
        i += 20
        nruns = blob[i]
        i += 1 + 8 * nruns
        out[(rid, mate)] = (status, t5a, t3a, t5p, t3p, rlen)
    return out


def run_fastprep(fp, seq_path, extra=()):
    d = os.path.dirname(seq_path)
    out1 = os.path.join(d, "out.fq")
    px = os.path.join(d, "out.pxmt")
    for p in (out1, px):
        if os.path.exists(p):
            os.remove(p)
    cmd = [fp, "-1", seq_path, "-o", out1, "--metadata", px,
           "--adapter-preset", "illumina", "--no-dedup", "--no-mean-q",
           "--min-len", "0", "-j", "1", "-q", *extra]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"fastprep failed: {r.returncode}\n{r.stderr[-2000:]}")
    return parse_pxmt(px), out1


def summarize(px, ids, t, label, extra=()):
    n = len(ids)
    exact = none = over = under = dropped = 0
    for rid in ids:
        st, t5a, t3a, t5p, t3p, _ = px[(rid, 0)]
        if st != 0:
            dropped += 1
            continue
        cut = t3a + t3p
        if cut == t:
            exact += 1
        elif cut == 0:
            none += 1
        elif cut < t:
            under += 1
        else:
            over += 1
    return dict(arm=label, tail=t, n=n, exact_tail=exact, no_trim=none,
                partial=under, past_truth=over, dropped=dropped)


def arm_tail(fp, rng, n, t, errs=0, g_contam=0, n_contam=0, adapter=TRUSEQ_R1,
             L=150, tag="A", extra=()):
    recs, ids, meta = [], [], {}
    for k in range(n):
        t_eff = min(t, L)
        bio_len = L - t_eff
        bio = dna(rng, bio_len)
        tail = list(adapter[:t_eff])
        # substitutions inside the adapter portion
        for _ in range(errs):
            if not tail:
                break
            j = rng.randrange(len(tail))
            tail[j] = rng.choice([b for b in "ACGT" if b != tail[j]])
        if n_contam and tail:
            for j in rng.sample(range(len(tail)), min(n_contam, len(tail))):
                tail[j] = "N"
        seq = list(bio + "".join(tail))
        if g_contam:
            for j in range(len(seq) - g_contam, len(seq)):
                seq[j] = "G"
        rid = f"{tag}{t}_{errs}_{g_contam}_{n_contam}_{k}"
        ids.append(rid)
        recs.append((rid, "".join(seq)))
        meta[rid] = t_eff
    d = tempfile.mkdtemp(prefix="ut_")
    sq = os.path.join(d, "in.fq")
    write_fastq(sq, recs)
    px, _ = run_fastprep(fp, sq, extra)
    rows = []
    for t_eff in sorted(set(meta.values())):
        sub = [r for r in ids if meta[r] == t_eff]
        r = summarize(px, sub, t_eff, tag)
        if extra:
            r["extra"] = list(extra)
        rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default="bench/results/undertrim-probe.json")
    a = ap.parse_args()
    fp = os.path.abspath(a.fastprep)
    rng = random.Random(a.seed)
    report = {"fastprep": fp, "n_per_cell": a.n, "arms": {}}

    # ---- Arm A: tail ladder -------------------------------------------------
    print("== Arm A: tail ladder (perfect adapter, no errors)")
    print(f"{'t':>3} {'exact':>7} {'none':>7} {'partial':>8} {'past':>6}")
    A = []
    for t in list(range(1, 21)) + [25, 30, 33, 40]:
        rs = arm_tail(fp, rng, a.n, t, tag="A")
        r = rs[0]
        A.append(r)
        print(f"{t:3d} {r['exact_tail']:7d} {r['no_trim']:7d} {r['partial']:8d} {r['past_truth']:6d}")
    report["arms"]["A_tail_ladder"] = A

    # ---- Arm A2: the same ladder at a lower floor ---------------------------
    print("\n== Arm A2: ladder re-run at --min-adapter-overlap 10 and 8")
    A2 = []
    print(f"{'mo':>3} " + " ".join(f"t={t:<2}" for t in range(6, 15)))
    for mo in (10, 8):
        cells, rows_ = [], []
        for t in range(6, 15):
            r = arm_tail(fp, rng, a.n, t, extra=["--min-adapter-overlap", str(mo)], tag=f"A2{mo}")[0]
            r["min_overlap"] = mo
            cells.append(r["exact_tail"] / r["n"])
            A2.append(r)
        print(f"{mo:3d} " + " ".join(f"{c*100:4.0f}" for c in cells))
    report["arms"]["A2_floor"] = A2

    # ---- Arm B: substitution errors inside the adapter ----------------------
    print("\n== Arm B: e substitutions inside the adapter tail")
    B = []
    print(f"{'t':>3} " + " ".join(f"e={e:<2}" for e in range(0, 6)))
    for t in (14, 20, 33):
        cells = []
        for e in range(0, 6):
            r = arm_tail(fp, rng, a.n, t, errs=e, tag="B")[0]
            cells.append(r["exact_tail"] / r["n"])
            B.append(r)
        print(f"{t:3d} " + " ".join(f"{c*100:5.1f}" for c in cells))
    report["arms"]["B_errors"] = B

    # ---- Arm C: contamination ----------------------------------------------
    print("\n== Arm C: contamination of a 20-base adapter tail")
    C = []
    for g in (0, 2, 4, 6, 8, 10, 12):
        r = arm_tail(fp, rng, a.n, 20, g_contam=g, tag="Cg")[0]
        C.append(r)
        print(f"  G-run {g:2d} -> exact {r['exact_tail']:5d}  none {r['no_trim']:5d}  partial {r['partial']:5d}")
    for nN in (1, 2, 3):
        r = arm_tail(fp, rng, a.n, 20, n_contam=nN, tag="Cn")[0]
        C.append(r)
        print(f"  N x{nN}   -> exact {r['exact_tail']:5d}  none {r['no_trim']:5d}  partial {r['partial']:5d}")
    # adapter-only read (insert 0)
    r = arm_tail(fp, rng, 200, 150, tag="C0")
    C.append(r)
    print(f"  adapter-only (insert 0, t=150) -> {r}")
    report["arms"]["C_contam"] = C

    # ---- Arm D: clean reads (false-trim price of the floor) -----------------
    print("\n== Arm D: adapter-free reads, false-trim rate vs --min-adapter-overlap")
    D = []
    recs, ids = [], []
    for k in range(max(a.n, 20000)):
        rid = f"D{k}"
        ids.append(rid)
        recs.append((rid, dna(rng, 150)))
    d = tempfile.mkdtemp(prefix="utc_")
    sq = os.path.join(d, "clean.fq")
    write_fastq(sq, recs)
    for mo in (20, 14, 12, 10, 8, 6, 5, 4):
        px, _ = run_fastprep(fp, sq, extra=["--min-adapter-overlap", str(mo)])
        trimmed = sum(1 for r in ids if px[(r, 0)][0] == 0 and (px[(r, 0)][2] + px[(r, 0)][4]) > 0)
        row = dict(min_overlap=mo, n=len(ids), false_trim=trimmed,
                   rate=trimmed / len(ids))
        D.append(row)
        print(f"  min_overlap {mo:2d}: {trimmed:6d}/{len(ids)} = {100*row['rate']:.3f}%"
              f" (1 in {len(ids)/max(trimmed,1):.0f})")
    report["arms"]["D_clean"] = D

    with open(a.out, "w") as f:
        json.dump(report, f, indent=1)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
