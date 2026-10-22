#!/usr/bin/env python3
"""Decompose FastPrep's undertrim at tolerance 0, per mate.

Table 7 of the paper states that the 10.3% undertrim splits into 5' runs of 1-4,
3' runs of 1-5 and 3' runs of 6 or more, with only the last being recoverable.
Those three numbers were typed from an ad-hoc analysis and no committed artifact
carried them, so this script recomputes them from the same corpus, binary and
config as the four-tool table (the `fastprep-hd0` arm at tolerance 0).

For each mate the ground truth wants trimmed and FastPrep left, the mate is
classified by the SIDE that undershot and by the length of the GT trim on that
side. A mate that undershoots on both sides is attributed to the 5' side, which
is where 99% of the overtrim also lives; the count of such mates is reported so
the arbitrariness is visible rather than buried.

Output: one TSV row per (side, gt-length) cell plus a `pooled` summary, so the
shares in the paper can be re-derived with awk rather than copied.
"""
import collections
import csv
import importlib.util
import os
import pathlib
import subprocess
import sys

HARNESS = pathlib.Path(__file__).resolve().parent / "overtrim_vs_fastp.py"
spec = importlib.util.spec_from_file_location("ovf", HARNESS)
ovf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ovf)

BUCKETS = [("5p", 1, 4, "5' runs of 1-4"),
           ("3p", 1, 5, "3' runs of 1-5"),
           ("3p", 6, 10 ** 9, "3' runs of 6 or more")]


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "undertrim-decomp.tsv"
    if os.path.exists(out):
        sys.exit(f"{out} exists; move it aside rather than mixing two runs")

    winners = list(csv.DictReader(open(f"{ovf.R}/winners.tsv"), delimiter="\t"))
    if not winners:
        sys.exit(f"no winners under {ovf.R}")
    print(f"datasets: {len(winners)}", flush=True)

    cell = collections.Counter()          # (side, gtlen) -> mates
    total = collections.Counter()         # every scored mate, for the pp delta
    under = 0
    both_sides = 0
    per_ds = []
    ran = 0
    for w in winners:
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{ovf.W}/{dom}/{acc}/v0"
        gtp, s1, s2 = (f"{v0}/sim.v0.gt.tsv", f"{v0}/sim.v0_R1.fastq.gz",
                       f"{v0}/sim.v0_R2.fastq.gz")
        if not (os.path.exists(gtp) and os.path.exists(s1)):
            continue
        gt = ovf.gt_rows(gtp)
        coremap = ovf.core_to_pid([s1, s2])
        d = os.path.join(os.environ.get("FASTSCRATCH_OVF", "/tmp/ovf"), f"ud_{dom}_{acc}")
        os.makedirs(d, exist_ok=True)
        base = ["--force", "-q", "--gz-level", "6"]
        for key, (flag, conv) in ovf.CONF.items():
            if w.get(key) not in (None, ""):
                base += [flag, conv(w[key])]
        if w.get("min_len") not in (None, ""):
            base += ["--min-len", ovf._minlen(w)]
        if w.get("r2_reverse") in ("True", "true", "1"):
            base += ["--r2-scan-reverse"]
        if w.get("interior") in ("True", "true", "1"):
            base += ["--interior"]
        # hd0: the rule that does not trim a leading G, i.e. the arm whose
        # undertrim is the 10.3% Table 7 decomposes.
        px = f"{d}/hd0.pxmt"
        r = subprocess.run([ovf.FP, "-1", s1, "-2", s2, "-o", f"{d}/R1.fq.gz",
                            "-p", f"{d}/R2.fq.gz", "--metadata", px,
                            "--head-dropout", "0"] + base,
                           capture_output=True)
        if r.returncode != 0:
            print(f"  {acc}: rc={r.returncode} {r.stderr[:120]!r}", flush=True)
            continue
        ran += 1
        p = ovf.parse_pxmt(px)
        n_u = 0
        for (rid, mate), (status, t5, t3) in p.items():
            pid = coremap.get((rid, mate))
            if pid is None or pid not in gt:
                continue
            h, t = (gt[pid][0], gt[pid][1]) if mate == 1 else (gt[pid][2], gt[pid][3])
            if status != 0:
                continue
            total["scored"] += 1
            u5, u3 = t5 < h, t3 < t
            if not (u5 or u3):
                continue
            under += 1
            n_u += 1
            if u5 and u3:
                both_sides += 1
            if u5:
                cell[("5p", h)] += 1
            else:
                cell[("3p", t)] += 1
        per_ds.append((dom, acc, n_u))
        subprocess.run(["rm", "-rf", d], capture_output=True)

    if ran == 0 or total["scored"] == 0:
        sys.exit(f"scored nothing: ran {ran} datasets, {total['scored']} mates")

    n = total["scored"]
    rows = []
    for side, lo, hi, label in BUCKETS:
        m = sum(v for (s, g), v in cell.items() if s == side and lo <= g <= hi)
        rows.append((label, side, lo, hi, m, 100 * m / n, 100 * m / under))
    print(f"scored {n:,} mates over {ran} datasets; undertrim {under:,} "
          f"({100 * under / n:.2f}%), of which both-sides {both_sides:,}", flush=True)

    with open(out, "w") as fh:
        fh.write("component\tside\tlo\thi\tmates\tpct_of_mates\tpct_of_undertrim\n")
        for label, side, lo, hi, m, pm, pu in rows:
            fh.write(f"{label}\t{side}\t{lo}\t{hi}\t{m}\t{pm:.4f}\t{pu:.2f}\n")
        fh.write(f"pooled\tside\t0\t0\t{under}\t{100 * under / n:.4f}\t100.00\n")
    for label, side, lo, hi, m, pm, pu in rows:
        print(f"  {label:<22} {m:>10,} mates  {pu:6.2f}% of undertrim  {pm:7.3f}% of mates")
    # What recovering the 3' tail entirely would buy in exact rate.
    tail = [r for r in rows if r[0].startswith("3' runs of 6")][0]
    print(f"  recovering all -> +{tail[5]:.3f} pp exact (the paper's 'about 0.09 pp')")
    print(f"wrote {out}")
    if both_sides:
        print(f"note: {both_sides} mates undershot on both sides and are counted 5'")


if __name__ == "__main__":
    main()
