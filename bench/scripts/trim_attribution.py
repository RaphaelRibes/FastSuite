#!/usr/bin/env python3
"""Attribute the undertrim bucket to its cause.

The trim comparison scores each read into exact / overtrim / undertrim /
unaligned (`overtrim_vs_fastp.py::bucket`). Undertrim means the trimmer kept
sequence the ground truth says is artefact: `t5 < h` or `t3 < t`.

This asks WHICH artefact. The GT carries, per mate, the adapter length
(`r{1,2}_adapter`) and the ambiguous-base count (`r{1,2}_n`), so a failing read
can be attributed to an untrimmed adapter, an untrimmed N-run, both, or
neither (a non-zero biology boundary with no adapter and no N — i.e. a polyG
tail the scan missed).

Usage:
  trim_attribution.py [--limit N] [--config-id ID] [--head-dropout V] [--tol 0]

Runs fastprep exactly as the harness does, so the totals reproduce the TSV.
"""
import argparse, collections, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import overtrim_vs_fastp as H

FP = os.environ.get("FASTPREP",
                    "/home/raphael/FastSuite/tools/fastprep/target/release/fastprep")


def gt_full(path):
    """-> {pair_id: (r1_h, r1_t, r2_h, r2_t, r1_adapter, r2_adapter,
                    r1_n, r2_n)}

    The harness's own gt_rows() keeps only the four boundaries; the undertrim
    attribution needs the adapter length and the ambiguous-base count too.
    """
    out = {}
    with open(path) as fh:
        for r in H.csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (
                int(r["r1_h"]), int(r["r1_t"]), int(r["r2_h"]), int(r["r2_t"]),
                int(r["r1_adapter"]), int(r["r2_adapter"]),
                int(r["r1_n"]), int(r["r2_n"]),
            )
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--config-id", default=None,
                    help="use this config row from the sweep instead of the winner")
    ap.add_argument("--head-dropout", default="0")
    ap.add_argument("--tol", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    winners = list(H.csv.DictReader(open(f"{H.R}/winners.tsv"), delimiter="\t"))
    if a.config_id is not None:
        rows = list(H.csv.DictReader(open(f"{H.R}/sweep_rows.tsv"), delimiter="\t"))
        pick = [r for r in rows if r["config_id"] == a.config_id]
        if not pick:
            sys.exit(f"config {a.config_id} not in sweep_rows.tsv")
        cfg = pick[0]
        winners = [dict(w, **{k: cfg[k] for k in cfg if k in
                              ('min_run','mercy','phred_sat','gate_shape','min_len',
                               'mean_q','adapter_err','min_overlap','r2_reverse',
                               'interior','head_dropout')}) for w in winners]
    if a.limit:
        winners = winners[:a.limit]

    cause = collections.Counter()
    direction = collections.Counter()
    per_domain = collections.defaultdict(collections.Counter)
    glen = collections.defaultdict(collections.Counter)   # G-run length -> count
    tot = collections.Counter()

    for i, w in enumerate(winners, 1):
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{H.W}/{dom}/{acc}/v0"
        gtp, s1, s2 = (f"{v0}/sim.v0.gt.tsv", f"{v0}/sim.v0_R1.fastq.gz",
                       f"{v0}/sim.v0_R2.fastq.gz")
        if not (os.path.exists(gtp) and os.path.exists(s1)):
            continue
        gt = gt_full(gtp)
        coremap = H.core_to_pid([s1, s2])
        d = f"/tmp/attr/{dom}_{acc}"
        os.makedirs(d, exist_ok=True)

        base = ["--force", "-q", "--gz-level", "6"]
        for key, (flag, conv) in H.CONF.items():
            if w.get(key) not in (None, ""):
                base += [flag, conv(w[key])]
        if w.get("min_len") not in (None, ""):
            base += ["--min-len", H._minlen(w)]
        if str(w.get("r2_reverse")) in ("True", "true", "1"):
            base += ["--r2-scan-reverse"]
        if str(w.get("interior")) in ("True", "true", "1"):
            base += ["--interior"]
        base += ["--head-dropout", str(a.head_dropout)]

        px = f"{d}/p.pxmt"
        r = subprocess.run([FP, "-1", s1, "-2", s2, "-o", f"{d}/o1.fq.gz",
                            "-p", f"{d}/o2.fq.gz", "--metadata", px] + base,
                           capture_output=True)
        if r.returncode != 0:
            print(f"  [{i}] {acc}: rc={r.returncode}", flush=True)
            continue
        p = H.parse_pxmt(px)

        for (rid, mate), (status, t5, t3) in p.items():
            pid = coremap.get((rid, mate))
            if pid is None or pid not in gt:
                continue
            g = gt[pid]
            h, t = (g[0], g[1]) if mate == 1 else (g[2], g[3])
            ad = g[4] if mate == 1 else g[5]      # adapter length, this mate
            n = g[6] if mate == 1 else g[7]       # ambiguous bases, this mate

            if status != 0:
                tot["dropped"] += 1
                continue
            if abs(t5 - h) <= a.tol and abs(t3 - t) <= a.tol:
                tot["exact"] += 1
                continue
            if t5 > h or t3 > t:
                tot["overtrim"] += 1
                continue
            # ---- undertrim ----
            tot["undertrim"] += 1
            five = t5 < h
            three = t3 < t
            if five and three: direction["both ends"] += 1
            elif five: direction["5' only"] += 1
            else: direction["3' only"] += 1

            has_ad = ad > 0
            has_n = n > 0
            if has_ad and has_n: k = "adapter + N"
            elif has_ad: k = "adapter"
            elif has_n: k = "N-run"
            else: k = "neither (polyG/boundary?)"
            cause[k] += 1
            per_domain[dom][k] += 1

            # For the G-run cases the length IS the recoverability question: the
            # Bayes analysis found P(artefact | run length) is a coin flip at
            # r=1 and approaches certainty by r>=4, so the length distribution
            # of the undertrim says how much of it is reachable at all.
            if not has_ad and not has_n:
                if five and not three: glen["5' G-run length"][h] += 1
                elif three and not five: glen["3' G-run length"][t] += 1
                elif five and three:
                    glen["5' G-run length"][h] += 1
                    glen["3' G-run length"][t] += 1

    n_under = tot["undertrim"]
    print(f"\n  scored {sum(tot.values())} mates: "
          f"exact {tot['exact']}, overtrim {tot['overtrim']}, "
          f"undertrim {n_under}, dropped {tot['dropped']}")
    if n_under:
        print(f"\n  UNDERTRIM by cause (of {n_under}):")
        for k, v in cause.most_common():
            print(f"    {k:26s} {v:8d}  {v/n_under*100:5.1f}%")
        print(f"\n  UNDERTRIM by end:")
        for k, v in direction.most_common():
            print(f"    {k:26s} {v:8d}  {v/n_under*100:5.1f}%")
        for end in ("5' G-run length", "3' G-run length"):
            c = glen.get(end)
            if not c: continue
            tt = sum(c.values())
            print(f"\n  {end} among undertrim (n={tt}):")
            for L in sorted(c, key=lambda x: int(x))[:8]:
                share = c[L]/tt*100
                bar = "#" * int(share/2)
                print(f"    r={int(L):>3}  {c[L]:8d}  {share:5.1f}%  {bar}")
            big = sum(v for L, v in c.items() if int(L) >= 4)
            print(f"    (r>=4, where P(artefact) approaches 1: {big} = {big/tt*100:.1f}%)")
        print(f"\n  UNDERTRIM by domain (share of that domain's undertrim):")
        for dom in sorted(per_domain):
            c = per_domain[dom]
            tt = sum(c.values())
            top = ", ".join(f"{k} {v/tt*100:.0f}%" for k, v in c.most_common(3))
            print(f"    {dom:18s} n={tt:8d}   {top}")


if __name__ == "__main__":
    main()
