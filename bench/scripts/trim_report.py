#!/usr/bin/env python3
"""Decompose the undertrim bucket into recoverable vs irreducible, and report
per-dataset / per-domain read RETENTION.

Why this exists
---------------
`docs/TRIM-ACCURACY.md` quotes undertrim as a single percentage (10.3% for
fastprep-hd-off at tolerance 0). Read plainly that says "the trimmer fails one
time in ten", which overstates the real deficiency by ~100x. Most of that bucket
is a G-run the simulator injected but a real trimmer cannot distinguish from
biology: `bench/scripts/trim_posterior.py` shows P(artefact | observed run
length) sits at 0.27-0.38 for short 5' runs, so trimming them is wrong ~73% of
the time. Only the tail where P(artefact | r) crosses 0.5 is genuinely
recoverable by a run-length-only rule.

This script splits the bucket into exactly two classes and prints, per tool:

    exact / overtrim / undertrim-recoverable / undertrim-irreducible / dropped

Definition of the split (the whole point of the script)
-------------------------------------------------------
  * RECOVERABLE  = a failing (undertrim) mate whose ground-truth artefact is a
    clean G-run -- no adapter, no ambiguous base, per the GT's per-mate
    `r{1,2}_adapter` / `r{1,2}_n` columns -- AND whose observed run length r at
    at least one failing end has P(artefact | r) > 0.5 in the measured
    posterior. A run-length-only rule COULD have caught it; today's default
    `min_run` (12) does not.
  * IRREDUCIBLE  = every other failing mate. That is the short/mid G-runs where
    P <= 0.5 (trimming would be wrong more often than right), plus the adapter
    and N-run cases the posterior does not model. Keeping these bases is the
    Bayes-optimal decision, and the ground truth can only label them artefact
    because a simulator injected them.

The posterior is `trim_posterior.py`'s, reused rather than reinvented: this
module imports its `leading_g` / `trailing_g` and replicates its tabulation
verbatim -- observed run length read off the RAW read, labelled "artefact" iff
the GT's injected run length (`r1_h`/`r1_t`) is > 0 at that end, tabulated per
end (5' and 3'). `trim_posterior.py` has no importable function returning the
table (its `main()` does everything inline), so the *maths* is replicated and
the *primitives* imported; the header of that script states this posterior is
the ceiling for any run-length-only rule.

Second output: RETENTION
------------------------
Dropped mates are scored by nobody, so two tools with different drop rates are
measured on different kept-sets. This script reads the same TSV that produced
the headline table and reports, per dataset and per domain, kept vs total
(`kept = scored`, `total = scored + dropped`, `retention = kept / total`) for
fastprep-hd0 against fastp-default -- i.e. how much of fastprep's exact-rate
lead sits on a smaller kept-set. Retention needs only the TSV, so it works even
when the ground-truth corpus is absent; the decomposition needs the corpus and
is skipped (with a clear message) when it is missing.

Outputs (never overwrite an existing .tsv unless --force)
---------------------------------------------------------
  <results>/trim-report-<YYYYMMDD>.tsv     the per-tool decomposition
  <results>/trim-retention-<YYYYMMDD>.tsv  per-domain + per-dataset retention

Third output: the PE-overlap slice
----------------------------------
Phase 1 of the PE-overlap plan (`docs/plans/…-pe-overlap-with-without.md`)
needs the SAME counters cut three extra ways, all of them additive: every
column the published table has keeps its name and its value.

  per mate   exact_R1 / exact_R2, overtrim_R1 / overtrim_R2,
             undertrim_rec_R1 / _R2, undertrim_irr_R1 / _R2 -- so the hypothesis
             "the effect lives on a specific mate" can be read off the table
             instead of re-running anything. The R1+R2 sums MUST equal the
             pooled column (asserted by `bench/scripts/trim_report_check.py`).
  by GT      readthrough (r1_adapter>0 or r2_adapter>0) and in_window
             (L <= insert <= 2L-12) -- the two conditions the mode keys on,
             from the gt's per-pair columns, per MATE (a pair is readthrough if
             either mate carried adapter).
  pe_trims   pairs where the tool's observed (5',3') cuts differ from the
             `fastprep-hd0` baseline for the same (pair, mate). With no
             --pe-overlap flag in fastprep yet this is ~0 BY CONSTRUCTION for
             the published configs; the column exists so that when the flag
             lands (Task 2.4: FASTPREP_VARIANTS='[["fp-pe",["--pe-overlap"]]]')
             it has somewhere to show up.

Usage
-----
  trim_report.py                       # retention + decomposition, full corpus
  trim_report.py --mode retention      # TSV only, no corpus needed
  trim_report.py --limit 2             # scope / time the run on N datasets
  trim_report.py --work W --results R  # override corpus / results roots
"""
import argparse
import collections
import csv
import datetime
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import overtrim_vs_fastp as H  # the harness: paths, keying, PXMT, scoring

try:
    # reuse the posterior's own run-length primitives (its tabulation is
    # replicated below because main() exposes no function returning the table).
    from trim_posterior import leading_g, trailing_g
    _POSTERIOR_SRC = "imported trim_posterior.leading_g/trailing_g"
except Exception:  # pragma: no cover - fallback keeps the script runnable
    def leading_g(seq):
        n = 0
        for b in seq:
            if b in "gG":
                n += 1
            else:
                break
        return n

    def trailing_g(seq):
        n = 0
        for b in reversed(seq):
            if b in "gG":
                n += 1
            else:
                break
        return n
    _POSTERIOR_SRC = "replicated trim_posterior maths (import failed)"

TOOLS = ("fastprep-hd0", "fastp-default", "fastp-tuned", "fastprep-hd4")

# The PE-overlap slice is built on a pair's ground-truth insert, so it needs
# the two gt columns below. read_len is optional (older gt files lack it) and
# falls back to --read-len; `insert` is not -- a gt without it cannot place a
# pair in the window, and that slice is then skipped for that file with a
# message instead of crashing (the decomposition itself never needs it).
GT_COLS = ("insert", "r1_adapter", "r2_adapter")
DEFAULT_READ_LEN = 150
DEFAULT_ANCHOR = 12


# --------------------------------------------------------------------------
# MEASUREMENT-ONLY extra fastprep configs (opt-in flags under test)
# --------------------------------------------------------------------------
# FASTPREP_VARIANTS is a JSON list of [tag, [args…]] entries. Each entry runs
# fastprep once more per dataset with the winner's own config, `--head-dropout
# 0` (the published baseline) and these extra args appended. UNSET/empty (the
# default) registers NOTHING: the decomposition then runs the same four
# configs, in the same order, with the same argv, as the published table.
#   FASTPREP_VARIANTS='[["fp-3p6",["--min-run-3p","6"]], \
#                       ["fp-n",["--trim-terminal-n"]]]' python3 trim_report.py …
VARIANTS = []


def load_variants():
    raw = os.environ.get("FASTPREP_VARIANTS", "").strip()
    if not raw:
        return []
    import json
    out = []
    for tag, args in json.loads(raw):
        out.append((str(tag), [str(a) for a in args]))
    return out


def gt_full(path):
    """-> {pair_id: (r1_h, r1_t, r2_h, r2_t, r1_ad, r2_ad, r1_n, r2_n)}

    The harness's gt_rows() keeps only the four boundaries; the split needs the
    adapter length and the ambiguous-base count too (as trim_attribution.py).
    """
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (
                int(r["r1_h"]), int(r["r1_t"]), int(r["r2_h"]), int(r["r2_t"]),
                int(r["r1_adapter"]), int(r["r2_adapter"]),
                int(r["r1_n"]), int(r["r2_n"]),
            )
    return out


def gt_meta(path):
    """-> (read_len, header) for a gt.tsv, or (None, header) if it has no read_len.

    The PE slice needs a per-pair read length; fastsim writes one into some gt
    files and not others. Absent the column the caller falls back to --read-len
    (150) rather than failing, since the decomposition never needs it.
    """
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if "read_len" not in header:
        return None, header
    rl = None
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                rl = int(f[header.index("read_len")])
                break
            except (ValueError, IndexError):
                continue
    return rl, header


def gt_flags(path, read_len):
    """-> (flags, skipped) for the PE slice of one dataset.

    flags is {pair_id: (readthrough, in_window)} -- a PAIR-level answer to two
    pair-level questions, applied to BOTH mates when counting (a pair is
    "readthrough" if either mate carried adapter, so the pair's mates are all
    readthrough mates for slicing purposes).

    skipped is a short reason string when the gt cannot answer: a missing
    `insert` column (older gt files) means a pair cannot be placed in the
    window, and an all-empty insert column means the same. In both cases the
    caller drops the in_window slice for this file and says so -- it must not
    crash the run, and it must not silently count the file as "no window
    pairs" (which would bias the slice towards the datasets that do have the
    column).
    """
    rl, header = gt_meta(path)
    if "insert" not in header:
        return {}, "no `insert` column"
    if rl is None:
        rl = read_len
    ceiling = 2 * rl - DEFAULT_ANCHOR
    flags, seen_insert, n = {}, 0, 0
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        fields = rd.fieldnames or []
        if "r1_adapter" not in fields or "r2_adapter" not in fields:
            return {}, "no adapter columns"
        for r in rd:
            n += 1
            try:
                ins = int(r["insert"])
                rt = int(r["r1_adapter"]) > 0 or int(r["r2_adapter"]) > 0
            except (ValueError, TypeError):
                continue
            seen_insert += 1
            flags[int(r["pair_id"])] = (rt, rl <= ins <= ceiling)
    if seen_insert == 0:
        return {}, "empty `insert` column"
    return flags, ""


# --------------------------------------------------------------------------
# Phase A -- the measured posterior P(artefact | observed run length)
# --------------------------------------------------------------------------
def build_posterior(winners, work, limit):
    """Tabulate observed G-run length vs whether the sim injected a G-run.

    Verbatim trim_posterior.py logic: obs is read off the RAW read; a read is
    labelled artefact iff the GT's injected length is > 0 at that end.
    """
    if limit:
        winners = winners[:limit]
    tab5 = collections.defaultdict(lambda: [0, 0])
    tab3 = collections.defaultdict(lambda: [0, 0])
    used = 0
    for w in winners:
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{work}/{dom}/{acc}/v0"
        gtp, s1, s2 = (f"{v0}/sim.v0.gt.tsv", f"{v0}/sim.v0_R1.fastq.gz",
                       f"{v0}/sim.v0_R2.fastq.gz")
        if not (os.path.exists(gtp) and os.path.exists(s1)):
            continue
        used += 1
        gt = gt_full(gtp)
        for path, mate in ((s1, 1), (s2, 2)):
            for pid, m, seq in H.read_fastq_keyed(path):
                g = gt.get(pid)
                if g is None:
                    continue
                inj5 = g[0] if mate == 1 else g[2]
                inj3 = g[1] if mate == 1 else g[3]
                o5 = leading_g(seq)
                if o5 > 0:
                    tab5[o5][0 if inj5 > 0 else 1] += 1
                o3 = trailing_g(seq)
                if o3 > 0:
                    tab3[o3][0 if inj3 > 0 else 1] += 1
    p5 = {r: na / (na + nb) for r, (na, nb) in tab5.items()}
    p3 = {r: na / (na + nb) for r, (na, nb) in tab3.items()}
    return p5, p3, tab5, tab3, used


# --------------------------------------------------------------------------
# Phase B -- per-tool decomposition
# --------------------------------------------------------------------------
def classify(t5, t3, h, t, ad, n, obs5, obs3, p5, p3, tol):
    """-> 'exact' | 'overtrim' | 'rec' | 'irr' | None (unknown/unaligned)."""
    if t5 is None:
        return "unaligned"
    if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
        return "exact"
    if t5 > h or t3 > t:
        return "overtrim"
    # ---- undertrim: split it ----
    five, three = t5 < h, t3 < t
    if ad > 0 or n > 0:
        return "irr"                      # adapter / N-run: not the posterior's model
    if five and p5.get(obs5, 0.0) > 0.5:
        return "rec"
    if three and p3.get(obs3, 0.0) > 0.5:
        return "rec"
    return "irr"


def decomp(winners, work, results, p5, p3, limit, run_fastp=True,
           read_len=DEFAULT_READ_LEN):
    """Run every tool on every dataset; return per (tool, tol) counters.

    Returns (c, cdom, cmate, cpe). `c`/`cdom` are the published decomposition;
    `cmate` (task 1.1) adds exact_R1/_R2, overtrim_R1/_R2, undertrim_rec_R1/_R2
    and undertrim_irr_R1/_R2 alongside them; `cpe` (tasks 1.2 + 1.3) adds
    pe_trims (cut decisions that differ from the fastprep-hd0 baseline) and the
    readthrough / in_window slices. None of the three changes what `c` counts:
    each is written into a fresh Counter keyed by its own prefix.
    """
    if limit:
        winners = winners[:limit]
    c = collections.defaultdict(collections.Counter)      # (tool,tol) -> Counter
    # MEASUREMENT-ONLY: per-domain counters, so a config's retention can be
    # read per domain (the kept-set question is concentrated: virome/rna_seq).
    # Written to $FASTPREP_DOMAIN_TSV when set; inert otherwise.
    cdom = collections.defaultdict(collections.Counter)   # (tool,tol,domain)
    cmate = collections.defaultdict(collections.Counter)  # (tool,tol) -> mate-split
    cpe = collections.defaultdict(collections.Counter)    # (tool,tol) -> PE slice
    # (pair, mate) -> (t5, t3) from the no-PE baseline, per dataset: the
    # reference pe_trims compares every config against. fastp has no baseline
    # row of its own, so its pe_trims is left at 0 and reported as "n/a" --
    # see the note in the report header.
    base_cuts = {}
    # (domain, acc) -> reason the in_window slice was skipped (task 1.2)
    pe_skipped = {}
    for i, w in enumerate(winners, 1):
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{work}/{dom}/{acc}/v0"
        gtp, s1, s2 = (f"{v0}/sim.v0.gt.tsv", f"{v0}/sim.v0_R1.fastq.gz",
                       f"{v0}/sim.v0_R2.fastq.gz")
        if not (os.path.exists(gtp) and os.path.exists(s1)):
            continue
        gt = gt_full(gtp)
        raw = {(pid, m): seq for pid, m, seq in
               (H.read_fastq_keyed(s1) + H.read_fastq_keyed(s2))}
        coremap = H.core_to_pid([s1, s2])
        gflags, gskip = gt_flags(gtp, read_len)
        if gskip:
            pe_skipped[(dom, acc)] = gskip
        d = f"/tmp/trim_report/{dom}_{acc}"
        os.makedirs(d, exist_ok=True)

        base = ["--force", "-q", "--gz-level", "6"]
        for key, (flag, conv) in H.CONF.items():
            if w.get(key) not in (None, ""):
                base += [flag, conv(w[key])]
        has_minlen = w.get("min_len") not in (None, "")
        if has_minlen:
            base += ["--min-len", H._minlen(w)]
        if str(w.get("r2_reverse")) in ("True", "true", "1"):
            base += ["--r2-scan-reverse"]
        if str(w.get("interior")) in ("True", "true", "1"):
            base += ["--interior"]

        # ---- fastprep-hd0 and fastprep-hd4 (PXMT gives exact cuts) ----
        # The two published configs keep their exact argv; VARIANTS add opt-in
        # flag sets on top of the newer of them (head-dropout 0).
        fp_tools = [("fastprep-hd0", ["--head-dropout", "0"]),
                    ("fastprep-hd4", ["--head-dropout", str(H.eff_hd(w))])]
        fp_tools += [(tag, ["--head-dropout", "0"] + args) for tag, args in VARIANTS]
        for tag, hd_args in fp_tools:
            argv = list(base) + list(hd_args)
            # A VARIANT may carry its own `--min-len` to OVERRIDE the winner's
            # (a min-len LADDER is expressed this way). clap rejects a duplicate
            # `--min-len` (rc=2), and the append hook cannot express an override,
            # so the winner's `--min-len <val>` pair is dropped for that variant
            # only. Inert when the variant has no `--min-len`: argv == base +
            # hd_args, so the published baselines keep their exact argv.
            if has_minlen and "--min-len" in hd_args:
                k = argv.index("--min-len")
                del argv[k:k + 2]
            px = f"{d}/{tag}.pxmt"
            r = subprocess.run([H.FP, "-1", s1, "-2", s2, "-o", f"{d}/{tag}_1.fq.gz",
                                "-p", f"{d}/{tag}_2.fq.gz", "--metadata", px]
                               + argv,
                               capture_output=True)
            if r.returncode != 0:
                print(f"  [{i}] {acc} {tag}: rc={r.returncode}", file=sys.stderr,
                      flush=True)
                continue
            p = H.parse_pxmt(px)
            cuts = {}
            for (rid, mate), (status, t5, t3) in p.items():
                pid = coremap.get((rid, mate))
                if pid is None or pid not in gt:
                    continue
                cuts[(pid, mate)] = (status, t5, t3)
            if tag == "fastprep-hd0":
                base_cuts = cuts
            for (pid, mate), (status, t5, t3) in cuts.items():
                g = gt[pid]
                h, t = (g[0], g[1]) if mate == 1 else (g[2], g[3])
                ad = g[4] if mate == 1 else g[5]
                n = g[6] if mate == 1 else g[7]
                obs5 = leading_g(raw.get((pid, mate), ""))
                obs3 = trailing_g(raw.get((pid, mate), ""))
                rt, win = gflags.get(pid, (None, None))
                for tol in (0, 2):
                    k = "dropped" if status != 0 else classify(
                        t5, t3, h, t, ad, n, obs5, obs3, p5, p3, tol)
                    c[(tag, tol)][k] += 1
                    cdom[(tag, tol, dom)][k] += 1
                    if mate == 1:
                        cmate[(tag, tol)][k + "_R1"] += 1
                    else:
                        cmate[(tag, tol)][k + "_R2"] += 1
                    # pe_trims: the cut decision this config made vs the
                    # baseline's, on the same (pair, mate). Both configs are
                    # in the same dataset, so base_cuts is populated by the time
                    # any variant is scored (fp_tools runs hd0 first).
                    if base_cuts:
                        b = base_cuts.get((pid, mate))
                        if b is not None and (status, t5, t3) != b:
                            cpe[(tag, tol)]["pe_trims"] += 1
                    if rt is not None:
                        cpe[(tag, tol)][f"{k}_readthrough_{int(rt)}"] += 1
                    if win is not None:
                        cpe[(tag, tol)][f"{k}_in_window_{int(win)}"] += 1

        # ---- fastp-default and fastp-tuned (cuts recovered by alignment) ----
        for tag, extra in ((("fastp-default", []), ("fastp-tuned", H.FASTP_TUNED))
                           if run_fastp else ()):
            o1, o2 = f"{d}/{tag}_1.fastq.gz", f"{d}/{tag}_2.fastq.gz"
            r = subprocess.run([H.FASTP, "-i", s1, "-I", s2, "-o", o1, "-O", o2,
                                "-w", "4", "-z", "6"] + extra, capture_output=True)
            if r.returncode != 0:
                print(f"  [{i}] {acc} {tag}: rc={r.returncode}", file=sys.stderr,
                      flush=True)
                continue
            tmap = {(pid, m): seq for pid, m, seq in
                    (H.read_fastq_keyed(o1) + H.read_fastq_keyed(o2))}
            for (pid, mate), seq in raw.items():
                if pid not in gt:
                    continue
                g = gt[pid]
                h, t = (g[0], g[1]) if mate == 1 else (g[2], g[3])
                ad = g[4] if mate == 1 else g[5]
                n = g[6] if mate == 1 else g[7]
                rt, win = gflags.get(pid, (None, None))
                tr = tmap.get((pid, mate))
                for tol in (0, 2):
                    if tr is None:
                        k = "dropped"
                    else:
                        cuts = H.infer_cuts(seq, tr)
                        if cuts is None:
                            k = "unaligned"
                        else:
                            k = classify(cuts[0], cuts[1], h, t, ad, n,
                                         leading_g(seq), trailing_g(seq), p5, p3, tol)
                    c[(tag, tol)][k] += 1
                    cdom[(tag, tol, dom)][k] += 1
                    if mate == 1:
                        cmate[(tag, tol)][k + "_R1"] += 1
                    else:
                        cmate[(tag, tol)][k + "_R2"] += 1
                    if rt is not None:
                        cpe[(tag, tol)][f"{k}_readthrough_{int(rt)}"] += 1
                    if win is not None:
                        cpe[(tag, tol)][f"{k}_in_window_{int(win)}"] += 1
        # measurement-only: every config's PXMT/FASTQ for this dataset has been
        # scored above, so free it — a 102-dataset × 6-config run otherwise
        # fills /tmp (overtrim_vs_fastp.py does the same). Numbers unaffected.
        subprocess.run(["rm", "-rf", d], capture_output=True)
        if i % 10 == 0:
            print(f"  [{i}/{len(winners)}] {acc}", file=sys.stderr, flush=True)
    return c, cdom, cmate, cpe, pe_skipped


# --------------------------------------------------------------------------
# Retention -- from the TSV, no corpus required
# --------------------------------------------------------------------------
def retention(tsv):
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    per_ds = collections.defaultdict(lambda: collections.defaultdict(dict))
    for r in rows:
        if int(r["tol"]) != 0:
            continue
        kept, dropped = int(r["scored"]), int(r["dropped"])
        per_ds[r["domain"]][r["acc"]][r["tool"]] = (kept, dropped, kept + dropped)
    return per_ds


def cmd_retention(tsv, out):
    per_ds = retention(tsv)
    dom = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    lines = []
    lines.append("=== RETENTION (mates kept / total), tolerance 0 ===")
    lines.append("total = kept + dropped; a tool's exact rate is earned on its own "
                 "kept-set")
    lines.append(f"{'domain':18} {'fastprep-hd0 retained':>22} "
                 f"{'fastp-default retained':>23} {'kept-set gap (pp)':>17}")
    hdr = ["domain", "dataset", "tool", "kept", "dropped", "total",
           "retention_pct"]
    tsv_rows = []
    for d in sorted(per_ds):
        for ds in per_ds[d]:
            for tool, (k, dr, tot) in per_ds[d][ds].items():
                dom[d][tool][0] += k
                dom[d][tool][1] += dr
                tsv_rows.append([d, ds, tool, k, dr, tot,
                                 f"{100 * k / tot:.2f}"])
        k0, d0 = dom[d]["fastprep-hd0"]
        k1, d1 = dom[d]["fastp-default"]
        r0 = 100 * k0 / (k0 + d0)
        r1 = 100 * k1 / (k1 + d1)
        lines.append(f"{d:18} {r0:14.1f}% ({k0}) {r1:18.1f}% ({k1}) "
                     f"{r0 - r1:17.1f}")
    K0 = sum(dom[d]["fastprep-hd0"][0] for d in dom)
    D0 = sum(dom[d]["fastprep-hd0"][1] for d in dom)
    K1 = sum(dom[d]["fastp-default"][0] for d in dom)
    D1 = sum(dom[d]["fastp-default"][1] for d in dom)
    lines.append(f"{'ALL':18} {100 * K0 / (K0 + D0):14.1f}% ({K0}) "
                 f"{100 * K1 / (K1 + D1):18.1f}% ({K1}) "
                 f"{100 * K0 / (K0 + D0) - 100 * K1 / (K1 + D1):17.1f}")
    # widest per-dataset asymmetries, to name the cases the caveat bites on
    gaps = []
    for d in per_ds:
        for ds in per_ds[d]:
            a = per_ds[d][ds].get("fastprep-hd0")
            b = per_ds[d][ds].get("fastp-default")
            if not a or not b:
                continue
            ga = 100 * a[0] / a[2]
            gb = 100 * b[0] / b[2]
            gaps.append((ga - gb, d, ds, ga, gb, a[2]))
    gaps.sort()
    lines.append("")
    lines.append("widest per-dataset kept-set asymmetries "
                 "(fastprep-hd0 retained - fastp-default retained):")
    for gap, d, ds, ga, gb, tot in gaps[:8]:
        lines.append(f"  {d:18} {ds:12} {ga:5.1f}% vs {gb:5.1f}%  "
                     f"gap {gap:6.1f} pp  (n={tot})")
    text = "\n".join(lines)
    print(text)
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(hdr)
        wr.writerows(tsv_rows)
    print(f"\nwrote {out} ({len(tsv_rows)} rows)")


def cmd_decomp(winners, work, results, limit, tol_show, out, run_fastp=True):
    print(f"building posterior ({_POSTERIOR_SRC}) ...", file=sys.stderr, flush=True)
    p5, p3, tab5, tab3, used = build_posterior(winners, work, limit)
    print(f"  posterior from {used} datasets; "
          f"5' crossing at r={next((r for r in sorted(p5) if p5[r] > 0.5), None)}, "
          f"3' crossing at r={next((r for r in sorted(p3) if p3[r] > 0.5), None)}",
          file=sys.stderr, flush=True)
    print("running tools ...", file=sys.stderr, flush=True)
    c, cdom = decomp(winners, work, results, p5, p3, limit,
                     run_fastp=run_fastp)

    lines = []
    lines.append("=== UNDERTRIM DECOMPOSITION (per tool) ===")
    lines.append(f"posterior: P(artefact | observed run length); {_POSTERIOR_SRC}")
    lines.append("recoverable = clean G-run (no adapter, no N) whose observed run "
                 "length has P>0.5;")
    lines.append("irreducible = every other failing mate (P<=0.5, or adapter/N-run)")
    hdr = ["tool", "tol", "exact", "overtrim", "undertrim_recoverable",
           "undertrim_irreducible", "dropped", "unaligned", "scored"]
    tsv_rows = []
    for tol in tol_show:
        lines.append("")
        lines.append(f"-- tolerance {tol} " + "-" * 60)
        lines.append(f"{'tool':17} {'exact':>12} {'overtrim':>12} "
                     f"{'und-rec':>12} {'und-irr':>12} {'dropped':>12}")
        for tag in (["fastprep-hd0", "fastprep-hd4"] + [t for t, _ in VARIANTS]
                    + ["fastp-default", "fastp-tuned"]):
            s = c[(tag, tol)]
            tot = s["exact"] + s["overtrim"] + s["rec"] + s["irr"] + s["unaligned"]
            if tot == 0:
                continue
            lines.append(
                f"{tag:17} {s['exact']:12d} {s['overtrim']:12d} "
                f"{s['rec']:12d} {s['irr']:12d} {s['dropped']:12d}")
            tsv_rows.append([tag, tol, s["exact"], s["overtrim"], s["rec"],
                             s["irr"], s["dropped"], s["unaligned"], tot])
    text = "\n".join(lines)
    print("\n" + text)
    with open(out, "w", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(hdr)
        wr.writerows(tsv_rows)
    print(f"\nwrote {out} ({len(tsv_rows)} rows)")

    # MEASUREMENT-ONLY per-domain retention (kept = scored), when requested.
    dtsv = os.environ.get("FASTPREP_DOMAIN_TSV", "").strip()
    if dtsv:
        drows = []
        for (tag, tol, dom), s in cdom.items():
            kept = (s["exact"] + s["overtrim"] + s["rec"] + s["irr"]
                    + s["unaligned"])
            dr = s["dropped"]
            if kept + dr == 0:
                continue
            drows.append([dom, tag, tol, kept, dr, kept + dr,
                          f"{100 * kept / (kept + dr):.2f}"])
        drows.sort(key=lambda r: (r[0], r[2], r[1]))
        with open(dtsv, "w", newline="") as fh:
            wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
            wr.writerow(["domain", "tool", "tol", "kept", "dropped", "total",
                         "retention_pct"])
            wr.writerows(drows)
        print(f"wrote {dtsv} ({len(drows)} rows)")
    return p5, p3


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=H.W, help="corpus root (…/work)")
    ap.add_argument("--results", default=H.R, help="results root (…/results)")
    ap.add_argument("--tsv", default=None,
                    help="overtrim TSV for retention (default <results>/overtrim-vs-fastp.tsv)")
    ap.add_argument("--mode", choices=["both", "retention", "decomp"], default="both")
    ap.add_argument("--limit", type=int, default=0,
                    help="scope the decomposition to the first N datasets (timing)")
    ap.add_argument("--no-fastp", action="store_true",
                    help="skip the two fastp configs (time); decomposition then "
                         "covers fastprep configs only")
    ap.add_argument("--tol", type=int, action="append", default=None,
                    help="tolerance to print (repeatable; default 0 and 2)")
    ap.add_argument("--out", default=None, help="decomposition TSV (default dated)")
    ap.add_argument("--out-retention", default=None, help="retention TSV (default dated)")
    ap.add_argument("--force", action="store_true", help="overwrite existing outputs")
    a = ap.parse_args()

    global VARIANTS
    VARIANTS = load_variants()
    if VARIANTS:
        print(f"extra fastprep configs (FASTPREP_VARIANTS): "
              f"{[t for t, _ in VARIANTS]}", file=sys.stderr)

    date = datetime.date.today().strftime("%Y%m%d")
    tsv = a.tsv or f"{a.results}/overtrim-vs-fastp.tsv"
    out_dec = a.out or f"{a.results}/trim-report-{date}.tsv"
    out_ret = a.out_retention or f"{a.results}/trim-retention-{date}.tsv"
    tols = a.tol if a.tol else [0, 2]

    for p in (out_dec, out_ret):
        if os.path.exists(p) and not a.force:
            sys.exit(f"trim_report: {p} exists (this writes a fresh report per "
                     f"run); pass a new path or --force.")

    if a.mode in ("both", "retention"):
        if not os.path.exists(tsv):
            sys.exit(f"trim_report: retention TSV not found: {tsv}")
        cmd_retention(tsv, out_ret)

    if a.mode in ("both", "decomp"):
        wpath = f"{a.results}/winners.tsv"
        if not os.path.exists(wpath):
            sys.exit(f"trim_report: winners.tsv not found: {wpath}")
        winners = list(csv.DictReader(open(wpath), delimiter="\t"))
        if not os.path.exists(f"{a.work}"):
            sys.exit(f"trim_report: corpus root not found: {a.work} "
                     f"(decomposition needs the ground-truth corpus; "
                     f"use --mode retention without it)")
        print("")
        cmd_decomp(winners, a.work, a.results, a.limit, tols, out_dec,
                   run_fastp=not a.no_fastp)


if __name__ == "__main__":
    main()
