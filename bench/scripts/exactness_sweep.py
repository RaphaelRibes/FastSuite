#!/usr/bin/env python3
"""exactness_sweep.py — the Tier-2/3 arm sweep over the GT corpus.

For every simulated dataset and every arm, run fastprep ONCE, read its trim
decisions back from PXMT (v3, unchanged) and score them against the simulator
ground truth with `score_weighted`'s model (the pre-registered 7:1 weights,
ladder plus per-rule break-even). Rows are appended to a JSONL raw file as they
complete, so an interrupted run resumes and `--score-only` rebuilds the TSV from
what finished; nothing is silently dropped.

Arms (all default-off in the tool; this script is what arms them):

    base  shipped default — the reference every arm is priced against
    ba    --boundary-anchored                 Tier 2
    m1    --overtrim-margin 1                  Tier 3a, N=1
    m2    --overtrim-margin 2                  Tier 3a, N=2
    stl   --short-tail-longest                 Tier 3b
    pe    --pe-ambiguous-smallest              Tier 3c
    all   --boundary-anchored --overtrim-margin 2 --short-tail-longest
          --pe-ambiguous-smallest

Every raw row records the machine, the fastprep binary md5 and the exact argv:
a bench row that cannot say which build produced it is not evidence.

    exactness_sweep.py --limit 3 --jobs 2      # smoke
    exactness_sweep.py --jobs 10               # the campaign
    exactness_sweep.py --score-only            # rebuild the TSV
"""

import argparse
import glob
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import score_weighted as sw  # noqa: E402

ARMS = {
    "base": [],
    "ba": ["--boundary-anchored"],
    "m1": ["--overtrim-margin", "1"],
    "m2": ["--overtrim-margin", "2"],
    "stl": ["--short-tail-longest"],
    "pe": ["--pe-ambiguous-smallest"],
    "all": ["--boundary-anchored", "--overtrim-margin", "2",
            "--short-tail-longest", "--pe-ambiguous-smallest"],
}


def discover(corpus):
    """-> [{"ds", "r1", "r2", "gt"}] for every sim.vN ground truth in `corpus`."""
    out = []
    for gt in sorted(glob.glob(os.path.join(corpus, "**", "*.gt.tsv"),
                               recursive=True)):
        d = os.path.dirname(gt)
        base = os.path.basename(gt)[: -len(".gt.tsv")]
        r1 = os.path.join(d, base + "_R1.fastq.gz")
        r2 = os.path.join(d, base + "_R2.fastq.gz")
        if not (os.path.exists(r1) and os.path.exists(r2)):
            continue
        out.append({
            "ds": os.path.relpath(d, corpus).replace(os.sep, "/"),
            "r1": r1, "r2": r2, "gt": gt,
        })
    return out


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_one(task):
    fp, fp_md5, ds, r1, r2, gt, arm, tols, workroot, min_frac = task
    d = tempfile.mkdtemp(dir=workroot, prefix=arm + "-")
    try:
        o1 = os.path.join(d, "o1.fq")
        o2 = os.path.join(d, "o2.fq")
        px = os.path.join(d, "cuts.pxmt")
        argv = [fp, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
                "--metadata", px, "-j", "1", "--quiet"] + ARMS[arm]
        t0 = time.time()
        r = subprocess.run(argv, capture_output=True, text=True)
        elapsed = time.time() - t0
        if r.returncode != 0:
            return {"ds": ds, "arm": arm, "fp_md5": fp_md5,
                    "error": (r.stderr or r.stdout)[-500:]}
        dec = sw.parse_pxmt(px)
        row = sw.read_gt(gt)
        rl = max((v[3] or 0) for v in dec.values())
        if rl <= 0:
            return {"ds": ds, "arm": arm, "fp_md5": fp_md5,
                    "error": "PXMT carries no read_len"}
        floor = int(min_frac * rl)
        tols_out = {}
        for tol in tols:
            tols_out[str(tol)] = sw.score_dataset(dec, row, rl, floor, tol)
        return {"ds": ds, "arm": arm, "read_len": rl, "min_len_floor": floor,
                "elapsed_s": round(elapsed, 3), "host": socket.gethostname(),
                "fp_md5": fp_md5, "argv": argv[1:], "tols": tols_out}
    except Exception as exc:  # a broken dataset must not kill the sweep
        return {"ds": ds, "arm": arm, "fp_md5": fp_md5,
                "error": f"{type(exc).__name__}: {exc}"}
    finally:
        shutil.rmtree(d, ignore_errors=True)


def run(args, datasets):
    done = set()
    if os.path.exists(args.raw):
        with open(args.raw) as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if "error" not in rec:
                    done.add((rec["ds"], rec["arm"]))
        print(f"resume: {len(done)} (dataset, arm) results already on disk")

    fp = os.path.abspath(args.fastprep)
    fp_md5 = md5_of(fp)
    os.makedirs(os.path.dirname(os.path.abspath(args.raw)), exist_ok=True)
    workroot = args.work
    os.makedirs(workroot, exist_ok=True)
    tols = [int(t) for t in args.tol.split(",")]

    tasks = [(fp, fp_md5, d["ds"], d["r1"], d["r2"], d["gt"], arm, tols,
              workroot, args.min_len_frac)
             for d in datasets for arm in args.arms
             if (d["ds"], arm) not in done]
    print(f"{len(datasets)} datasets x {len(args.arms)} arms = "
          f"{len(datasets) * len(args.arms)} runs; {len(tasks)} to do now")
    if not tasks:
        return

    t0 = time.time()
    errors = 0
    with open(args.raw, "a") as fh, \
            ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futs = [pool.submit(run_one, t) for t in tasks]
        for i, fut in enumerate(as_completed(futs), 1):
            rec = fut.result()
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if "error" in rec:
                errors += 1
                print(f"  FAIL {rec['ds']} {rec['arm']}: {rec['error'][:200]}")
            if i % 25 == 0 or i == len(tasks):
                rate = (time.time() - t0) / i
                print(f"  {i}/{len(tasks)} runs, {rate:.2f} s/run, "
                      f"eta {(len(tasks) - i) * rate / 60:.1f} min", flush=True)
    print(f"done: {len(tasks)} runs, {errors} failed, "
          f"{(time.time() - t0) / 60:.1f} min")


def score(args):
    recs = []
    with open(args.raw) as fh:
        for line in fh:
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                continue
            if "error" not in rec:
                recs.append(rec)
    if not recs:
        raise SystemExit(f"no successful rows in {args.raw}")
    hosts = sorted({r["host"] for r in recs})
    md5s = sorted({r["fp_md5"] for r in recs})
    print(f"rows: {len(recs)}  hosts: {hosts}  fastprep md5: {md5s}")
    labels = sorted({r["arm"] for r in recs})
    rl = Counter(r["read_len"] for r in recs).most_common(1)[0][0]
    floor = Counter(r["min_len_floor"] for r in recs).most_common(1)[0][0]
    print(f"read_len {rl}, min_len floor {floor}, arms {labels}")
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)

    with open(args.weights) as fh:
        w = sw.load_weights(args.weights)

    summary = []
    for tol in [int(t) for t in args.tol.split(",")]:
        meas = {}
        for r in recs:
            m = r["tols"].get(str(tol))
            if m:
                meas[(r["ds"], r["arm"])] = m
        cfg = {"tol": tol, "min_len_frac": args.min_len_frac,
               "min_len_floor": floor, "read_len": rl,
               "baseline": args.baseline,
               "variants": [l for l in labels if l != args.baseline]}
        path = args.out.replace(".tsv", f".tol{tol}.tsv")
        with open(path, "w", newline="") as fh:
            sw.emit_rows(w, cfg, meas, fh)
        print(f"wrote {path}  ({len(meas)} dataset x arm measurements)")

        # ---- the acceptance table, printed: per arm, pooled counts
        for lab in labels:
            # mixed_trim is a PAIR counter folded into the numerator, so it is
            # summed for the score but never added to the mate denominator
            agg = {k: 0 for k in sw.CLASSES + ("mixed_trim",)}
            extra = {k: 0 for k in ("pairs_scored", "clean_mates",
                                    "clean_false_trims", "bases_over",
                                    "bases_under")}
            m5 = m3 = 0
            for (ds, al), m in meas.items():
                if al != lab:
                    continue
                for k in sw.CLASSES + ("mixed_trim",):
                    agg[k] += m["counts"][k]
                for k in extra:
                    extra[k] += m[k]
            mates = sum(agg[k] for k in sw.CLASSES) or 1
            num, sc = sw.score_at(agg, float(w["anchor_r"]), w["w_over"],
                                  extra["pairs_scored"] or 1)
            summary.append((tol, lab, 100.0 * agg["exact"] / mates,
                            100.0 * agg["overtrim"] / mates,
                            100.0 * agg["undertrim"] / mates,
                            (extra["clean_false_trims"] / extra["clean_mates"]
                             if extra["clean_mates"] else 0.0),
                            sc, num))
    print("\n== pooled acceptance table (per mate, % of mates scored) ==")
    print(f"{'tol':>3} {'arm':>5} {'exact':>7} {'over':>7} {'under':>7} "
          f"{'cleanFT':>8} {'score@7':>10}")
    base_score = {}
    for tol, lab, ex, ov, un, cft, sc, num in summary:
        if lab == args.baseline:
            base_score[tol] = sc
        d = sc - base_score.get(tol, sc)
        print(f"{tol:>3} {lab:>5} {ex:>6.2f}% {ov:>6.2f}% {un:>6.2f}% "
              f"{100 * cft:>7.3f}% {sc:>10.4f}  ({d:+.4f} vs base)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", default="/mnt/Fast/fastsuite-scratch/work-gt2026")
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--work", default="/mnt/Fast/fastsuite-scratch/work-exactness")
    ap.add_argument("--raw", default="bench/results/exactness/sweep-raw.jsonl")
    ap.add_argument("--out", default="bench/results/exactness/sweep.tsv")
    ap.add_argument("--weights", default=os.path.join(
        HERE, os.pardir, "weights", "trim-weights-v1.json"))
    ap.add_argument("--arms", default="base,ba,m1,m2,stl,pe,all")
    ap.add_argument("--tol", default="0,2")
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--min-len-frac", type=float, default=0.9)
    ap.add_argument("--baseline", default="base")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--score-only", action="store_true")
    args = ap.parse_args(argv)
    args.arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    unknown = [a for a in args.arms if a not in ARMS]
    if unknown:
        raise SystemExit(f"unknown arm(s) {unknown}; known: {sorted(ARMS)}")

    datasets = discover(args.corpus)
    if args.limit:
        datasets = datasets[: args.limit]
    if not datasets:
        raise SystemExit(f"no datasets under {args.corpus}")
    print(f"corpus {args.corpus}: {len(datasets)} datasets")
    if not args.score_only:
        run(args, datasets)
    score(args)


if __name__ == "__main__":
    main()
