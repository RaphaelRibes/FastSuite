#!/usr/bin/env python3
"""C3 — the winner vs fastp on the REAL reads (PLAN-ML-LOCAL C3).

Both tools run on the real 100k-pair unit and both outputs are fastcheck-ed,
so the comparison is on fastcheck's exact numbers rather than on either tool's
self-reported statistics. Three rows per dataset:

    fastprep-winner   the C1 winner's argv, verbatim
    fastp-default     fastp with no flags beyond -I/-O
    fastp-best        fastp with a best-effort tuning pass (--cut_front
                      --cut_tail --qualified_quality_phred 20 --length_required
                      30 --low_complexity_filter), the settings its own
                      documentation recommends for NovaSeq-like polyG data

Quality columns come from the fastcheck report of each output (GC, mean Q,
polyX5 fraction, reads surviving). Speed and memory are wall time and peak RSS
of the tool process itself.

Note on fairness: fastp's defaults are genuinely good and this comparison is
not a rigged race — fastp-best is included precisely so that "fastp loses" can
never mean "fastp was misconfigured". If the winner only beats fastp-default
and loses to fastp-best, the memo says so.

Usage:
  fastp_compare.py --fastcheck PATH --fastprep PATH --fastp PATH \\
                   --units ROOT --winners winners.tsv --row fastp-comparison.tsv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import resource
import subprocess
import sys
import time

ROW_FIELDS = ["domain", "run_accession", "tool", "variant", "argv",
              "reads_in", "reads_out", "mean_q", "gc_fraction",
              "polyx5_frac_gt_1", "polyx3_frac_gt_5", "duplication_level",
              "seconds", "peak_rss_mb", "status", "note"]

# fastp's best-effort pass: the settings its documentation recommends for
# NovaSeq-style libraries with polyG at the 5' end.
FASTP_BEST = ["--cut_front", "--cut_front_window_size", "5",
              "--cut_front_mean_quality", "20",
              "--cut_tail", "--cut_tail_window_size", "5",
              "--cut_tail_mean_quality", "20",
              "--qualified_quality_phred", "20",
              "--length_required", "30",
              "--low_complexity_filter"]


def peak_rss_mb() -> float:
    """Peak RSS of this process AND its reaped children, in MB."""
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    return ru.ru_maxrss / 1024.0


def run_timed(argv: list[str], log: str) -> tuple[int, float, str]:
    t0 = time.time()
    with open(log, "w") as fh:
        p = subprocess.run(argv, stdout=fh, stderr=subprocess.STDOUT)
    return p.returncode, time.time() - t0, log


def fastcheck_report(fastcheck: str, r1: str, r2: str | None, out_json: str,
                     threads: int) -> dict | None:
    argv = [fastcheck, "--json", out_json, "-q", "-j", str(threads), "--"]
    argv += [r1, r2] if r2 else [r1]
    p = subprocess.run(argv, capture_output=True)
    if p.returncode != 0 or not os.path.exists(out_json):
        return None
    with open(out_json) as fh:
        return json.load(fh)


def reads_out_of(report: dict | None) -> int | None:
    if not report:
        return None
    v = report.get("metrics", {}).get("reads")
    return int(v) if isinstance(v, (int, float)) else None


def load_winners(path: str) -> dict[str, dict]:
    """winners.tsv -> {domain/acc: row} (the C1 aggregate's winner rows)."""
    out: dict[str, dict] = {}
    if not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[f"{r['domain']}/{r['run_accession']}"] = r
    return out


def winner_argv(w: dict) -> list[str]:
    """Rebuild fastprep's argv from the winning config's columns."""
    a = ["--min-run", str(w["min_run"]), "--mercy", str(w["mercy"]),
         "--phred-sat", str(w["phred_sat"]), "--gate-shape", str(w["gate_shape"]),
         "--mean-q", str(w["mean_q"]),
         "--adapter-error-rate", str(w["adapter_err"]),
         "--min-adapter-overlap", str(w["min_overlap"])]
    ml = w["min_len"]
    try:
        f = float(ml)
        a += ["--min-len", str(int(f)) if f.is_integer() else f"{f}L"]
    except ValueError:
        a += ["--min-len", str(ml)]
    if str(w.get("r2_reverse", "0")) == "1":
        a.append("--r2-scan-reverse")
    if str(w.get("interior", "0")) == "1":
        a.append("--interior")
    return a


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastcheck", required=True)
    ap.add_argument("--fastprep", required=True)
    ap.add_argument("--fastp", required=True)
    ap.add_argument("--units", required=True)
    ap.add_argument("--winners", required=True)
    ap.add_argument("--row", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--work", default=None, help="scratch dir (default: beside --row)")
    args = ap.parse_args()

    work = args.work or os.path.join(os.path.dirname(args.row) or ".", "c3work")
    os.makedirs(work, exist_ok=True)
    os.makedirs(os.path.dirname(args.row) or ".", exist_ok=True)
    if not os.path.exists(args.row):
        with open(args.row, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(ROW_FIELDS)

    winners = load_winners(args.winners)
    if not winners:
        print("fastp_compare: no winners loaded — run aggregate.py first", file=sys.stderr)
        return 1

    def emit(domain, acc, tool, variant, argv, reads_in, report, secs, status, note=""):
        m = report.get("metrics", {}) if report else {}
        row = [domain, acc, tool, variant, " ".join(argv),
               reads_in, reads_out_of(report),
               m.get("mean_read_mean_q", ""), m.get("gc_fraction", ""),
               m.get("polyx5_frac_gt_1", ""), m.get("polyx3_frac_gt_5", ""),
               m.get("duplication_level", ""),
               f"{secs:.2f}", f"{peak_rss_mb():.1f}", status, note]
        with open(args.row, "a", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(row)

    n_units = 0
    for domain in sorted(os.listdir(args.units)):
        dpath = os.path.join(args.units, domain)
        if not os.path.isdir(dpath):
            continue
        for acc in sorted(os.listdir(dpath)):
            adir = os.path.join(dpath, acc)
            r1 = os.path.join(adir, f"{acc}.R1.fastq.gz")
            r2 = os.path.join(adir, f"{acc}.R2.fastq.gz")
            if not (os.path.exists(r1) and os.path.exists(r2)):
                continue
            key = f"{domain}/{acc}"
            if key not in winners:
                continue
            w = winners[key]
            tag = f"{domain}_{acc}"
            d = os.path.join(work, tag)
            os.makedirs(d, exist_ok=True)

            base = fastcheck_report(args.fastcheck, r1, r2, os.path.join(d, "in.json"), args.threads)
            reads_in = reads_out_of(base)

            # ---- 1. the winner ----
            o1 = os.path.join(d, "win_R1.fastq.gz")
            o2 = os.path.join(d, "win_R2.fastq.gz")
            argv = [args.fastprep, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
                    "-j", str(args.threads), "--force", "-q"] + winner_argv(w)
            rc, secs, log = run_timed(argv, os.path.join(d, "win.log"))
            rep = fastcheck_report(args.fastcheck, o1, o2, os.path.join(d, "win.json"), args.threads) \
                if rc == 0 else None
            emit(domain, acc, "fastprep-winner", w.get("config_id", ""), argv,
                 reads_in, rep, secs, "ok" if rc == 0 else "fail",
                 "" if rc == 0 else f"rc={rc}")

            # ---- 2. fastp defaults ----
            o1 = os.path.join(d, "fp_R1.fastq.gz")
            o2 = os.path.join(d, "fp_R2.fastq.gz")
            argv = [args.fastp, "-i", r1, "-I", r2, "-o", o1, "-O", o2,
                    "-w", str(args.threads)]
            rc, secs, log = run_timed(argv, os.path.join(d, "fp.log"))
            rep = fastcheck_report(args.fastcheck, o1, o2, os.path.join(d, "fp.json"), args.threads) \
                if rc == 0 else None
            emit(domain, acc, "fastp-default", "", argv, reads_in, rep, secs,
                 "ok" if rc == 0 else "fail", "" if rc == 0 else f"rc={rc}")

            # ---- 3. fastp best-effort ----
            o1 = os.path.join(d, "fpb_R1.fastq.gz")
            o2 = os.path.join(d, "fpb_R2.fastq.gz")
            argv = [args.fastp, "-i", r1, "-I", r2, "-o", o1, "-O", o2,
                    "-w", str(args.threads)] + FASTP_BEST
            rc, secs, log = run_timed(argv, os.path.join(d, "fpb.log"))
            rep = fastcheck_report(args.fastcheck, o1, o2, os.path.join(d, "fpb.json"), args.threads) \
                if rc == 0 else None
            emit(domain, acc, "fastp-best", "", argv, reads_in, rep, secs,
                 "ok" if rc == 0 else "fail", "" if rc == 0 else f"rc={rc}")

            for f in os.listdir(d):
                if f.endswith(".fastq.gz"):
                    try:
                        os.unlink(os.path.join(d, f))
                    except OSError:
                        pass
            n_units += 1
            print(f"fastp_compare: {tag} done ({n_units})")

    print(f"fastp_compare: {n_units} datasets -> {args.row}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
