#!/usr/bin/env python3
"""auto_hd_sweep.py — price a PER-LIBRARY head-dropout policy end to end.

The 5' G-run decision is undecidable per read (prefix_features) but the injected
share is a library property (prefix_policy): it runs from 0% to 63% across the
probe, and at 7:1 a library should keep its 1-4 bp leading runs unless at least
one in eight of them is injected. One global `--head-dropout` applies one answer
to all of them, which is why the corpus-level score is a compromise.

This script measures the policy a tool could actually ship: for each library it
prescans the first reads of R1, estimates the injected share as the excess of the
1-4 bp leading-G rate over what the library's own G frequency predicts
(p - p^5), and then runs the REAL binary with `--head-dropout 4` when the excess
clears 1/8 and `--head-dropout 0` when it does not. Scoring is `exactness_sweep`'s
own (`score_weighted` at the 7:1 anchor, drops counted), so the arms are directly
comparable to @tbl-trim-acc's campaign.

The estimator never sees the ground truth: it reads only the FASTQ.

    auto_hd_sweep.py --raw /tmp/auto.jsonl --jobs 6
"""
import argparse
import csv
import gzip
import json
import os
import sys
import tempfile
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exactness_sweep as es   # noqa: E402
import overtrim_vs_fastp as ovf  # noqa: E402
import score_weighted as sw    # noqa: E402

W_OVER, W_UNDER = 1.0, 7.0
THRESH = W_OVER / (W_OVER + W_UNDER)      # 1/8: break-even injected share
PRESCAN = 100_000


def leading_run(seq):
    n = 0
    for b in seq:
        if b in "gG":
            n += 1
        else:
            break
    return n


def prescan(path, max_reads=PRESCAN):
    """p_G over the reads, and the rate of reads starting with 1-4 G bases."""
    mates = 0
    amb = 0
    g = total = 0
    with gzip.open(path, "rt") as fh:
        while mates < max_reads:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().rstrip("\n").upper()
            fh.readline()
            fh.readline()
            mates += 1
            if 1 <= leading_run(seq) <= 4:
                amb += 1
            g += seq.count("G")
            total += len(seq)
    if mates == 0 or total == 0:
        raise RuntimeError(f"prescan read nothing from {path}")
    pg = g / total
    obs = amb / mates
    exp_bio = max(0.0, min(1.0, pg - pg ** 5))
    est_all = max(0.0, min(1.0, obs - exp_bio))
    # The break-even (1 in 8) is about the mates the rule acts on, i.e. the
    # ambiguous class, not about all mates. est_all is a share of ALL mates, so
    # dividing by the ambiguous rate converts it back to the class the decision
    # is about. Without this the estimator reads 5.7% for a library whose
    # ambiguous class is 29% injected and keeps a library it should trim.
    est_cls = min(1.0, est_all / obs) if obs > 0 else 0.0
    return {"mates": mates, "p_g": pg, "obs_1_4": obs, "exp_bio": exp_bio,
            "est_injected_all": est_all, "est_injected_cls": est_cls,
            "trim": est_cls > THRESH}


def run_one(task):
    (fp, fp_md5, ds, r1, r2, gt, cfg, workroot, min_frac, anchor) = task
    d = tempfile.mkdtemp(dir=workroot, prefix="auto-")
    try:
        census = prescan(r1)
        hd = 4 if census["trim"] else 0
        cfg = list(cfg)
        if "--head-dropout" in cfg:
            i = cfg.index("--head-dropout")
            del cfg[i:i + 2]
        cfg += ["--head-dropout", str(hd)]
        if anchor and "--boundary-anchored" not in cfg:
            cfg += ["--boundary-anchored"]
        o1, o2 = os.path.join(d, "o1.fq"), os.path.join(d, "o2.fq")
        px = os.path.join(d, "cuts.pxmt")
        argv = [fp, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
                "--metadata", px, "-j", "1"] + cfg
        t0 = time.time()
        r = subprocess.run(argv, capture_output=True, text=True)
        elapsed = time.time() - t0
        if r.returncode != 0:
            return {"ds": ds, "error": (r.stderr or r.stdout)[-400:], "census": census}
        dec = sw.parse_pxmt(px)
        grow = sw.read_gt(gt)
        rl = max((x[3] or 0) for x in dec.values())
        floor = int(min_frac * rl)
        return {"ds": ds, "hd": hd, "census": census, "read_len": rl,
                "elapsed_s": round(elapsed, 3), "fp_md5": fp_md5,
                "tols": {"0": sw.score_dataset(dec, grow, rl, floor, 0)}}
    except Exception as exc:
        return {"ds": ds, "error": f"{type(exc).__name__}: {exc}"}
    finally:
        shutil.rmtree(d, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--work", default="/tmp/auto-hd")
    ap.add_argument("--anchor", action="store_true",
                    help="add --boundary-anchored, so the policy is comparable "
                         "to the shipped default rather than to the unanchored "
                         "arms of the campaign")
    a = ap.parse_args()
    if os.path.exists(a.raw):
        sys.exit(f"{a.raw} exists; move it aside")

    winners = list(csv.DictReader(open(f"{ovf.R}/winners.tsv"), delimiter="\t"))
    tasks = []
    for w in winners:
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{ovf.W}/{dom}/{acc}/v0"
        r1, r2, gt = f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz", f"{v0}/sim.v0.gt.tsv"
        if not all(os.path.exists(p) for p in (r1, r2, gt)):
            continue
        tasks.append((ovf.FP, es.md5_of(ovf.FP), f"{dom}/{acc}", r1, r2, gt,
                      es.winner_flags(w), a.work, 0.9, a.anchor))
    if not tasks:
        sys.exit("no datasets resolved")
    os.makedirs(a.work, exist_ok=True)
    print(f"{len(tasks)} datasets", flush=True)

    recs = []
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        for i, rec in enumerate(ex.map(run_one, tasks), 1):
            recs.append(rec)
            if "error" in rec:
                print(f"  FAIL {rec['ds']}: {rec['error'][:160]}", flush=True)
            if i % 10 == 0:
                print(f"  {i}/{len(tasks)}", flush=True)
    with open(a.raw, "w") as fh:
        for rec in recs:
            fh.write(json.dumps(rec) + "\n")

    ok = [r for r in recs if "error" not in r]
    if not ok:
        sys.exit("every dataset failed")
    tot = {"exact": 0, "overtrim": 0, "undertrim": 0, "drop_correct": 0,
           "drop_wrong": 0, "keep_wrong": 0}
    mates = 0
    trimmed_libs = sum(1 for r in ok if r["hd"] == 4)
    for r in ok:
        tot.update({k: tot[k] + r["tols"]["0"]["counts"][k] for k in tot})
        mates += r["tols"]["0"]["mates_scored"]
    score = (tot["exact"] - tot["overtrim"] - 7 * tot["undertrim"]
             + tot["drop_correct"] - tot["drop_wrong"]) / mates
    print(f"\ndatasets {len(ok)}  mates {mates:,}  libraries trimmed {trimmed_libs}"
          f"  kept {len(ok) - trimmed_libs}")
    print(f"  exact {100*tot['exact']/mates:.2f}%  over {100*tot['overtrim']/mates:.2f}%"
          f"  under {100*tot['undertrim']/mates:.2f}%  score@7 {score:+.4f}")
    print(f"wrote {a.raw}")


if __name__ == "__main__":
    main()
