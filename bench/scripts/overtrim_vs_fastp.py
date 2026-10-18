#!/usr/bin/env python3
"""fastprep (hd=4 / hd=0) vs fastp (default / tuned) on the SAME 101 real
libraries, same ground truth, same harness.

Why this exists: the published trim-accuracy table scored fastprep from the
sweep's stored PXMT and fastp from its own run, and its fastprep baseline
(66.5%) does not reproduce in a fresh harness (75.35%). Any claim of the form
"fastp beats fastprep" has to be measured with one harness or it is comparing
two methodologies, not two tools.

Scoring:
  fastprep  exact per-mate trim positions, read back from its PXMT metadata.
  fastp     reports no trim positions, so the trimmed read is aligned back
            onto its original: a pure trim is a substring, so the offset and
            the tail length are recovered exactly. Reads fastp dropped are
            counted separately and excluded from the rates.

Read keying is identical for both: the pair id comes from the sim header
(x + 10000*y), verified against the GT, so neither tool is advantaged by the
pairing.
"""
import argparse, collections, csv, gzip, os, shlex, struct, subprocess, sys

W = "/mnt/Fast/fastsuite-scratch/work"
R = "/mnt/Fast/fastsuite-scratch/results"
FP = "/home/raphael/FastSuite/tools/fastprep/target/release/fastprep"
FASTP = "/home/raphael/.local/bin/fastp"

# MEASUREMENT-ONLY hook: extra fastprep args, shell-quoted, appended after the
# winner's own config and after `--head-dropout`. UNSET (the default) is a
# no-op: every invocation below is then byte-for-byte the published one, so
# the four baseline configs stay reproducible to the mate. Set it to score an
# opt-in flag on the same corpus, same ground truth, same read keying:
#   FASTPREP_EXTRA_ARGS="--min-run-3p 6 --trim-terminal-n" python3 ...
FP_EXTRA = shlex.split(os.environ.get("FASTPREP_EXTRA_ARGS", ""))

FASTP_TUNED = ["--cut_front", "--cut_front_window_size", "5",
               "--cut_front_mean_quality", "20", "--cut_tail",
               "--cut_tail_window_size", "5", "--cut_tail_mean_quality", "20",
               "--qualified_quality_phred", "20", "--length_required", "30",
               "--low_complexity_filter"]

COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def gt_rows(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out


def read_fastq_keyed(path):
    """-> [(pair_id, mate, seq)] keyed by the sim header's pair id."""
    out = []
    if not os.path.exists(path):
        return out
    pid = mate = None
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                parts = line[1:].split()
                f = parts[0].split(":")
                pid = int(f[5]) + 10000 * int(f[6])
                mate = int(parts[1][0])
            elif i % 4 == 1 and pid is not None:
                out.append((pid, mate, line.rstrip("\n")))
    return out


def core_to_pid(paths):
    out = {}
    for p in paths:
        if not os.path.exists(p):
            continue
        with gzip.open(p, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    core = line[1:].split()[0]
                    f = core.split(":")
                    out[(core, int(line.split()[1][0]))] = (
                        int(f[5]) + 10000 * int(f[6]))
    return out


def parse_pxmt(path):
    out = {}
    try:
        buf = open(path, "rb").read()
    except FileNotFoundError:
        return out
    if buf[:4] != b"PXMT":
        return out
    off, n = 5, len(buf)
    while off + 2 <= n:
        (idlen,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + idlen].decode("utf-8", "replace").split()[0]
        off += idlen
        if off + 22 > n:
            break
        mate = buf[off] + 1; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out[(rid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def infer_cuts(orig, trimmed):
    if trimmed == "":
        return (len(orig), 0)
    i = orig.find(trimmed)
    if i >= 0:
        return (i, len(orig) - i - len(trimmed))
    rc = "".join(COMP.get(c, "N") for c in reversed(trimmed))
    j = orig.find(rc)
    if j >= 0:
        return (j, len(orig) - j - len(rc))
    return None


def bucket(c, t5, t3, h, t, tol):
    if t5 is None:
        c["unaligned"] += 1
    elif abs(t5 - h) <= tol and abs(t3 - t) <= tol:
        c["exact"] += 1
    elif t5 > h or t3 > t:
        c["overtrim"] += 1
    else:
        c["undertrim"] += 1


def newb():
    return collections.Counter(exact=0, overtrim=0, undertrim=0,
                               unaligned=0, dropped=0)


def score_fastprep(pxmt, coremap, gt, tol):
    c = newb()
    for (rid, mate), (status, t5, t3) in pxmt.items():
        pid = coremap.get((rid, mate))
        if pid is None or pid not in gt:
            continue
        h, t = (gt[pid][0], gt[pid][1]) if mate == 1 else (gt[pid][2], gt[pid][3])
        if status != 0:
            c["dropped"] += 1
            continue
        bucket(c, t5, t3, h, t, tol)
    return c


def score_fastp(orig, trimmed, gt, tol):
    tmap = {(p, m): s for p, m, s in trimmed}
    c = newb()
    for pid, mate, seq in orig:
        if pid not in gt:
            continue
        h, t = (gt[pid][0], gt[pid][1]) if mate == 1 else (gt[pid][2], gt[pid][3])
        tr = tmap.get((pid, mate))
        if tr is None:
            c["dropped"] += 1
            continue
        cuts = infer_cuts(seq, tr)
        if cuts is None:
            c["unaligned"] += 1
            continue
        bucket(c, cuts[0], cuts[1], h, t, tol)
    return c


def eff_hd(cfg):
    hd = int(cfg.get("head_dropout") or 4)
    return max(1, min(hd, int(cfg["min_run"]) - 1))


def _minlen(w):
    """The exact `--min-len` value the sweep ran with, from winners.tsv.

    sweep.py does NOT pass the grid's raw min_len. An integer grid value is a
    base count authored against a 150 bp reference; `scale_min_len()` converts
    it to a FRACTION of the actual read length, clamped to <= 1.0, and the
    resolved string is recorded in the winner's `min_len_effective` column.
    Passing the raw integer as a literal base count (the old `_minlen`) runs a
    different configuration whenever the clamp bites — SRR19649201 as-swept
    `--min-len 1` vs raw `--min-len 50` (75.9% vs 97.7% exact), DRR000432
    (78.0% vs 99.7%). Read the recorded column so the harness reproduces the
    swept run.
    """
    eff = w.get("min_len_effective")
    if eff not in (None, ""):
        return eff
    raise SystemExit(f"winner {w.get('run_accession')}: no min_len_effective column")


CONF = {
    "min_run": ("--min-run", str), "mercy": ("--mercy", str),
    "phred_sat": ("--phred-sat", str), "gate_shape": ("--gate-shape", str),
    "mean_q": ("--mean-q", str),
    "adapter_err": ("--adapter-error-rate", str),
    "min_overlap": ("--min-adapter-overlap", str),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing --out instead of refusing")
    a = ap.parse_args()

    # Refuse to append — checked FIRST, before any work, so a re-run against an
    # existing output costs nothing. `emit` below opens --out in mode "a" so
    # that sharded runs (--stride/--offset) can share one file. The consequence,
    # when the file already exists from a previous run, is that this run's rows
    # are SILENTLY MIXED with that run's. That is not hypothetical: it produced
    # a file in which every key appeared twice — once from the stale
    # pre-regeneration corpus and once from the new one — and aggregating it
    # returned numbers that were a blend of two different corpora. Same class of
    # silent-stale-data bug as the chain.sh artefact reuse.
    if a.offset == 0 and os.path.exists(a.out) and not a.force:
        sys.exit(
            f"overtrim_vs_fastp: {a.out} already exists and this script APPENDS "
            f"to its output (sharded runs share one file), so re-running would "
            f"mix the old rows with the new ones. Delete it, pass a new --out, "
            f"or pass --force to overwrite."
        )
    if a.offset == 0 and a.force and os.path.exists(a.out):
        os.remove(a.out)

    winners = list(csv.DictReader(open(f"{R}/winners.tsv"), delimiter="\t"))
    winners = winners[a.offset::a.stride]
    print(f"datasets: {len(winners)}", flush=True)

    # Write the header (the append-refusal happens at the top of main(), before
    # any work, so that a re-run against an existing path costs nothing).
    if a.offset == 0:
        with open(a.out, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(
                ["domain", "acc", "tool", "tol", "exact", "overtrim",
                 "undertrim", "unaligned", "dropped", "scored"])
    else:
        import time
        for _ in range(600):
            if os.path.exists(a.out):
                break
            time.sleep(1)

    def emit(dom, acc, tool, c):
        for tol in (0, 2):
            sub = c[tol]
            tot = sub["exact"] + sub["overtrim"] + sub["undertrim"] + sub["unaligned"]
            with open(a.out, "a", newline="") as fh:
                csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(
                    [dom, acc, tool, tol, sub["exact"], sub["overtrim"],
                     sub["undertrim"], sub["unaligned"], sub["dropped"], tot])

    for i, w in enumerate(winners, 1):
        dom, acc = w["domain"], w["run_accession"]
        v0 = f"{W}/{dom}/{acc}/v0"
        gtp, s1, s2 = (f"{v0}/sim.v0.gt.tsv", f"{v0}/sim.v0_R1.fastq.gz",
                       f"{v0}/sim.v0_R2.fastq.gz")
        if not (os.path.exists(gtp) and os.path.exists(s1)):
            continue
        gt = gt_rows(gtp)
        coremap = core_to_pid([s1, s2])
        d = os.path.join(os.environ.get("FASTSCRATCH_OVF", "/tmp/ovf"), f"{a.offset}_{dom}_{acc}")
        os.makedirs(d, exist_ok=True)

        base = ["--force", "-q", "--gz-level", "6"]
        for key, (flag, conv) in CONF.items():
            if w.get(key) not in (None, ""):
                base += [flag, conv(w[key])]
        # min-len is resolved by sweep.py's scale_min_len and recorded as
        # min_len_effective; reproduce the swept run from that column, not the
        # raw grid integer.
        if w.get("min_len") not in (None, ""):
            base += ["--min-len", _minlen(w)]
        if w.get("r2_reverse") in ("True", "true", "1"):
            base += ["--r2-scan-reverse"]
        if w.get("interior") in ("True", "true", "1"):
            base += ["--interior"]

        # ---- fastprep: the then-shipped value (head_dropout 4), then 0
        for tag, hdv in (("fastprep-hd4", None), ("fastprep-hd0", 0)):
            argv = list(base)
            argv += ["--head-dropout", str(eff_hd(w) if hdv is None else hdv)]
            px = f"{d}/{tag}.pxmt"
            r = subprocess.run([FP, "-1", s1, "-2", s2, "-o", f"{d}/{tag}_R1.fq.gz",
                                "-p", f"{d}/{tag}_R2.fq.gz", "--metadata", px]
                               + argv + FP_EXTRA, capture_output=True)
            if r.returncode != 0:
                print(f"  [{i}] {acc} {tag}: rc={r.returncode}", flush=True)
                continue
            p = parse_pxmt(px)
            emit(dom, acc, tag, {t: score_fastprep(p, coremap, gt, t)
                                 for t in (0, 2)})

        # ---- fastp: shipped default and tuned, scored by alignment
        orig = read_fastq_keyed(s1) + read_fastq_keyed(s2)
        for tag, extra in (("fastp-default", []), ("fastp-tuned", FASTP_TUNED)):
            o1, o2 = f"{d}/{tag}_R1.fastq.gz", f"{d}/{tag}_R2.fastq.gz"
            r = subprocess.run([FASTP, "-i", s1, "-I", s2, "-o", o1, "-O", o2,
                                "-w", "4", "-z", "6"] + extra,
                               capture_output=True)
            if r.returncode != 0:
                print(f"  [{i}] {acc} {tag}: rc={r.returncode}", flush=True)
                continue
            trimmed = read_fastq_keyed(o1) + read_fastq_keyed(o2)
            emit(dom, acc, tag, {t: score_fastp(orig, trimmed, gt, t)
                                 for t in (0, 2)})
        subprocess.run(["rm", "-rf", d], capture_output=True)
        if i % 10 == 0:
            print(f"  [{i}/{len(winners)}] {acc}", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
