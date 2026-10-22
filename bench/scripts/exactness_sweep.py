#!/usr/bin/env python3
"""exactness_sweep.py — the Tier-2/3 arm sweep over the GT corpus.

For every simulated dataset and every arm, run fastprep ONCE, read its trim
decisions back from PXMT (v3, unchanged) and score them against the simulator
ground truth with `score_weighted`'s pre-registered model (7:1, ladder plus
per-rule break-even). Rows are appended to a JSONL raw file as they complete, so
an interrupted run resumes and `--score-only` rebuilds the TSV from what
finished; nothing is silently dropped.

**Every arm runs on that dataset's OWN swept winner config** (the phase-D
`winners.tsv` row), assembled by the same helpers the authoritative harness uses
(`overtrim_vs_fastp.CONF` / `_minlen` / `eff_hd`). That is what makes an arm's
delta a property of the arm rather than of a config I made up: `base` is the
winner plus nothing.

Arms (all default-off in the tool; this script is what arms them):

    base  the swept winner, no extra flag — the reference every arm is priced
          against
    ba    --boundary-anchored                     Tier 2
    m1    --overtrim-margin 1                      Tier 3a, N=1
    m2    --overtrim-margin 2                      Tier 3a, N=2
    poly  --boundary-anchored --overtrim-margin 2  the polyX package
    st6   --adapter-short-tail 6                   short-tail path ON, shortest
                                                   wins = today's semantics on
                                                   that path (reference for stl)
    stl   --adapter-short-tail 6 --short-tail-longest   Tier 3b
    peov  --adapter-preset illumina --pe-overlap always  reference for pe
    pe    --adapter-preset illumina --pe-overlap always
          --pe-ambiguous-smallest                  Tier 3c

No single `all` arm spans the axes: stl and pe each need a prerequisite flag, so
a combined arm would price the prerequisite rather than the rule. Each Tier-3b/3c
arm is priced against its own reference (st6 / peov), which isolates the rule.

Every raw row records the machine, the fastprep binary md5 and the exact argv: a
bench row that cannot say which build produced it is not evidence.

    exactness_sweep.py --limit 3 --jobs 2      # smoke
    exactness_sweep.py --jobs 10 --variants v0 # the campaign
    exactness_sweep.py --score-only            # rebuild the TSV
"""

import argparse
import csv
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
import score_weighted as sw          # noqa: E402
import overtrim_vs_fastp as ovf      # noqa: E402  (the authoritative flag set)

ARMS = {
    "base": [],
    "ba": ["--boundary-anchored"],
    "m1": ["--overtrim-margin", "1"],
    "m2": ["--overtrim-margin", "2"],
    "poly": ["--boundary-anchored", "--overtrim-margin", "2"],
    "st6": ["--adapter-short-tail", "6"],
    "stl": ["--adapter-short-tail", "6", "--short-tail-longest"],
    "peov": ["--adapter-preset", "illumina", "--pe-overlap", "always"],
    "pe": ["--adapter-preset", "illumina", "--pe-overlap", "always",
           "--pe-ambiguous-smallest"],
}

# arm -> the arm it is priced against when it enables a prerequisite flag
REFERENCE = {"stl": "st6", "pe": "peov"}


def winner_flags(row):
    """The swept winner's own flags, exactly as the authoritative harness builds
    them (`overtrim_vs_fastp.py`): reproduce the config, never re-guess it."""
    base = ["--force", "-q", "--gz-level", "6"]
    for key, (flag, conv) in ovf.CONF.items():
        val = row.get(key)
        if val not in (None, ""):
            base += [flag, conv(val)]
    if row.get("min_len") not in (None, ""):
        base += ["--min-len", ovf._minlen(row)]
    if row.get("r2_reverse") in ("True", "true", "1"):
        base += ["--r2-scan-reverse"]
    if row.get("interior") in ("True", "true", "1"):
        base += ["--interior"]
    base += ["--head-dropout", str(ovf.eff_hd(row))]
    return base


def load_winners(path):
    with open(path) as fh:
        return {(r["domain"], r["run_accession"]): r
                for r in csv.DictReader(fh, delimiter="\t")}


def discover(corpus, variants):
    """-> [{"ds", "dom", "acc", "v", "r1", "r2", "gt"}] for every sim ground truth."""
    out = []
    for gt in sorted(glob.glob(os.path.join(corpus, "**", "*.gt.tsv"),
                               recursive=True)):
        d = os.path.dirname(gt)
        stem = os.path.basename(gt)[: -len(".gt.tsv")]      # sim.v0
        v = stem.split(".")[-1]                             # v0
        if variants and v not in variants:
            continue
        r1 = os.path.join(d, stem + "_R1.fastq.gz")
        r2 = os.path.join(d, stem + "_R2.fastq.gz")
        if not (os.path.exists(r1) and os.path.exists(r2)):
            continue
        rel = os.path.relpath(d, corpus).replace(os.sep, "/")
        parts = rel.split("/")
        out.append({"ds": rel, "dom": parts[0], "acc": parts[1], "v": v,
                    "r1": r1, "r2": r2, "gt": gt})
    return out


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_one(task):
    (fp, fp_md5, ds, r1, r2, gt, arm, cfg, tols, workroot, min_frac) = task
    d = tempfile.mkdtemp(dir=workroot, prefix=arm + "-")
    try:
        o1 = os.path.join(d, "o1.fq")
        o2 = os.path.join(d, "o2.fq")
        px = os.path.join(d, "cuts.pxmt")
        argv = ([fp, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
                 "--metadata", px, "-j", "1"] + cfg + ARMS[arm])
        t0 = time.time()
        r = subprocess.run(argv, capture_output=True, text=True)
        elapsed = time.time() - t0
        if r.returncode != 0:
            return {"ds": ds, "arm": arm, "fp_md5": fp_md5,
                    "error": (r.stderr or r.stdout)[-500:]}
        dec = sw.parse_pxmt(px)
        grow = sw.read_gt(gt)
        rl = max((x[3] or 0) for x in dec.values())
        if rl <= 0:
            return {"ds": ds, "arm": arm, "fp_md5": fp_md5,
                    "error": "PXMT carries no read_len"}
        floor = int(min_frac * rl)
        tols_out = {str(tol): sw.score_dataset(dec, grow, rl, floor, tol)
                    for tol in tols}
        return {"ds": ds, "arm": arm, "read_len": rl, "min_len_floor": floor,
                "elapsed_s": round(elapsed, 3), "host": socket.gethostname(),
                "fp_md5": fp_md5, "argv": argv[1:], "tols": tols_out}
    except Exception as exc:      # one broken dataset must not kill the sweep
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

    winners = load_winners(args.winners)
    fp = os.path.abspath(args.fastprep)
    fp_md5 = md5_of(fp)
    os.makedirs(os.path.dirname(os.path.abspath(args.raw)), exist_ok=True)
    os.makedirs(args.work, exist_ok=True)
    tols = [int(t) for t in args.tol.split(",")]

    tasks, missing = [], []
    for d in datasets:
        w = winners.get((d["dom"], d["acc"]))
        if w is None:
            missing.append(f"{d['dom']}/{d['acc']}")
            continue
        cfg = winner_flags(w)
        for arm in args.arms:
            if (d["ds"], arm) in done:
                continue
            tasks.append((fp, fp_md5, d["ds"], d["r1"], d["r2"], d["gt"],
                          arm, cfg, tols, args.work, args.min_len_frac))
    if missing:
        print(f"WARNING: {len(set(missing))} accession(s) have no winner row: "
              f"{sorted(set(missing))[:5]}")
    print(f"{len(datasets)} datasets x {len(args.arms)} arms; "
          f"{len(tasks)} runs to do now")
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
                print(f"  FAIL {rec['ds']} {rec['arm']}: {rec['error'][:300]}")
            if i % 50 == 0 or i == len(tasks):
                rate = (time.time() - t0) / i
                print(f"  {i}/{len(tasks)} runs, {rate:.2f} s/run, "
                      f"eta {(len(tasks) - i) * rate / 60:.1f} min", flush=True)
    print(f"done: {len(tasks)} runs, {errors} failed, "
          f"{(time.time() - t0) / 60:.1f} min")
    if errors:
        raise SystemExit(f"{errors} run(s) failed — fix before scoring")


def pooled(meas, label):
    agg = {k: 0 for k in sw.CLASSES + ("mixed_trim",)}
    extra = {k: 0 for k in ("pairs_scored", "clean_mates", "clean_false_trims",
                            "bases_over", "bases_under")}
    for (ds, al), m in meas.items():
        if al != label:
            continue
        for k in sw.CLASSES + ("mixed_trim",):
            agg[k] += m["counts"][k]
        for k in extra:
            extra[k] += m[k]
    return agg, extra


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
    labels = [a for a in ARMS if a in {r["arm"] for r in recs}]
    rl = Counter(r["read_len"] for r in recs).most_common(1)[0][0]
    floor = Counter(r["min_len_floor"] for r in recs).most_common(1)[0][0]
    print(f"read_len {rl}, min_len floor {floor}, arms {labels}")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    w = sw.load_weights(args.weights)

    for tol in [int(t) for t in args.tol.split(",")]:
        meas = {}
        for r in recs:
            m = r["tols"].get(str(tol))
            if m:
                meas[(r["ds"], r["arm"])] = m
        cfg = {"tol": tol, "min_len_frac": args.min_len_frac,
               "min_len_floor": floor, "read_len": rl, "baseline": "base",
               "variants": [l for l in labels if l != "base"]}
        path = args.out.replace(".tsv", f".tol{tol}.tsv")
        with open(path, "w", newline="") as fh:
            sw.emit_rows(w, cfg, meas, fh)
        print(f"\nwrote {path}  ({len(meas)} dataset x arm measurements)")

        # ---- acceptance table: pooled per arm, plus each Tier-3 arm against
        # its own reference so a prerequisite flag is never priced as the rule
        print(f"\n== tol {tol}: pooled per mate ==")
        print(f"{'arm':>5} {'n_ds':>5} {'exact':>7} {'over':>7} {'under':>7} "
              f"{'drop_w':>7} {'cleanFT':>8} {'score@7':>10} {'vs base':>9} "
              f"{'vs ref':>9}")
        base_score = None
        for lab in labels:
            agg, extra = pooled(meas, lab)
            mates = sum(agg[k] for k in sw.CLASSES) or 1
            num, sc = sw.score_at(agg, float(w["anchor_r"]), w["w_over"],
                                  extra["pairs_scored"] or 1)
            if lab == "base":
                base_score = sc
            ref = REFERENCE.get(lab, "base")
            ragg, rextra = pooled(meas, ref)
            _, rsc = sw.score_at(ragg, float(w["anchor_r"]), w["w_over"],
                                 rextra["pairs_scored"] or 1)
            n_ds = len({ds for (ds, a) in meas if a == lab})
            cft = (extra["clean_false_trims"] / extra["clean_mates"]
                   if extra["clean_mates"] else 0.0)
            print(f"{lab:>5} {n_ds:>5} "
                  f"{100 * agg['exact'] / mates:>6.2f}% "
                  f"{100 * agg['overtrim'] / mates:>6.2f}% "
                  f"{100 * agg['undertrim'] / mates:>6.2f}% "
                  f"{100 * agg['drop_wrong'] / mates:>6.2f}% "
                  f"{100 * cft:>7.3f}% {sc:>10.4f} "
                  f"{sc - (base_score if base_score is not None else sc):>+9.4f} "
                  f"{sc - rsc:>+9.4f}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--corpus", default="/mnt/Fast/fastsuite-scratch/work-gt2026")
    ap.add_argument("--winners",
                    default="bench/results/phase-d/winners.tsv",
                    help="per-dataset swept winner config (phase D)")
    ap.add_argument("--fastprep", default="tools/fastprep/target/release/fastprep")
    ap.add_argument("--work", default="/mnt/Fast/fastsuite-scratch/work-exactness")
    ap.add_argument("--raw", default="bench/results/exactness/sweep-raw.jsonl")
    ap.add_argument("--out", default="bench/results/exactness/sweep.tsv")
    ap.add_argument("--weights", default=os.path.join(
        HERE, os.pardir, "weights", "trim-weights-v1.json"))
    ap.add_argument("--arms", default="base,ba,m1,m2,poly,st6,stl,peov,pe")
    ap.add_argument("--tol", default="0,2")
    ap.add_argument("--variants", default="v0",
                    help="comma list of sim variants (default v0: the corpus "
                         "the published baseline used); 'all' for every one")
    ap.add_argument("--jobs", type=int, default=10)
    ap.add_argument("--min-len-frac", type=float, default=0.9)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--score-only", action="store_true")
    args = ap.parse_args(argv)
    args.arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    unknown = [a for a in args.arms if a not in ARMS]
    if unknown:
        raise SystemExit(f"unknown arm(s) {unknown}; known: {sorted(ARMS)}")
    variants = [] if args.variants == "all" else args.variants.split(",")

    datasets = discover(args.corpus, variants)
    if args.limit:
        datasets = datasets[: args.limit]
    if not datasets:
        raise SystemExit(f"no datasets under {args.corpus} for {args.variants}")
    print(f"corpus {args.corpus}: {len(datasets)} datasets "
          f"({args.variants or 'all variants'})")
    if not args.score_only:
        run(args, datasets)
    score(args)


if __name__ == "__main__":
    main()
