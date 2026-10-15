#!/usr/bin/env python3
"""C1 — the fastprep config sweep, scored exactly against fastsim's ground truth.

One simulated dataset (R1/R2 + gt.tsv) is trimmed by every config in the grid;
each config's trim decisions are read back from fastprep's `--metadata` (PXMT,
one record per read, documented in tools/fastprep/src/metadata.rs) and compared
to the truth the GT defines:

    keep window per mate  = [h, min(I, L-t))
    full_flag=1           -> the correct answer is to DROP the pair
    is_dup=1              -> the correct answer is to drop it too (a PCR
                             duplicate re-sequences an already-seen fragment)

so accuracy is absolute, not relative to another tool. Scores per config:

  keep_precision   kept pairs whose keep window is (near-)correct / kept pairs
  keep_recall      correctly-kept pairs / all pairs that should be kept
  trim_mae         mean |kept_left_boundary - h| over kept R1 reads (bp)
  over_trim_rate   kept pairs that lost >= 25 % of their biology window
  under_trim_rate  kept pairs that left >= 10 bp of artefact at either end
  drop_accuracy    1 - |kept - should_keep| / should_keep  (full+dup class)
  score            the single number the winner is chosen on:
                       keep_recall * keep_precision - 2*over_trim_rate
                   (recall and precision both matter; over-trimming destroys
                   real biology and is penalised twice as hard as it is common)

The grid: 55 authored configs x 3 variants (PLAN-ML-STEPS C1 said "60"; the
authored grid is 55, and the count is reported against the authored length so
a short grid is a visible mismatch, never a silent "55/60 failed"). fastprep's shipped
default (small-trim-rest-max) is known to be the swept winner from the earlier
dedup bench, so the grid explores AROUND it on the axes that matter for the
artefact regime the profile describes: polyX aggressiveness (min-run, mercy,
phred-sat, gate-shape), the length floor (min-len), the quality gate
(mean-q), and adapter handling (error rate, overlap).

Usage:
  sweep.py --acc ACC --domain D --variant V --fastprep PATH \\
           --r1 sim.v0_R1.fastq.gz --r2 sim.v0_R2.fastq.gz --gt sim.v0.gt.tsv \\
           --work DIR --row ROW.tsv --threads 4 --configs 60
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import struct
import subprocess
import sys
import time

# ---------------------------------------------------------------- the grid
# (min_run, mercy, phred_sat, gate_shape, min_len, mean_q, adapter_err, min_overlap, r2_reverse, interior)
# min_len is a FRACTION of the read length for the fractions and a base count
# for the integers. The integers are scaled by scale_min_len() below, because
# the pilot's units are 110 bp: a hardcoded 30 means "keep 27% of a short read"
# on 110 bp data and "20%" on 150 bp data, which would make the two read
# regimes non-comparable.
DEFAULTS = dict(min_run=12, mercy=3, phred_sat=20, gate_shape=1.0,
                min_len=0.9, mean_q=25.0, adapter_err=0.2, min_overlap=14,
                r2_reverse=False, interior=False)

# The shipped default is config #0 and is always evaluated first.
GRID = [
    DEFAULTS,                                        # 0  shipped default
    # --- polyX aggressiveness axis ---
    dict(min_run=8),  dict(min_run=10), dict(min_run=14), dict(min_run=16), dict(min_run=20),
    dict(mercy=1),    dict(mercy=2),   dict(mercy=5),    dict(mercy=8),
    dict(phred_sat=15), dict(phred_sat=25), dict(phred_sat=30), dict(phred_sat=35),
    dict(gate_shape=0.5), dict(gate_shape=1.5), dict(gate_shape=2.0),
    # --- length floor axis ---
    dict(min_len=0.8), dict(min_len=0.85), dict(min_len=0.95), dict(min_len=30), dict(min_len=50),
    # --- quality gate axis ---
    dict(mean_q=20.0), dict(mean_q=22.0), dict(mean_q=28.0), dict(mean_q=30.0), dict(mean_q=15.0),
    # --- adapter axis ---
    dict(adapter_err=0.1), dict(adapter_err=0.3),
    dict(min_overlap=8),  dict(min_overlap=10), dict(min_overlap=20),
    # --- scan direction (Nextera/NovaSeq put polyG at R2's 5' end) ---
    dict(r2_reverse=True), dict(r2_reverse=True, min_run=10),
    dict(r2_reverse=True, mean_q=20.0),
    # --- interior scan (off in the shipped default) ---
    dict(interior=True), dict(interior=True, min_run=10),
    # --- head dropout (NovaSeq first-cycle 5' G) ---
    dict(min_run=10, head_dropout=4), dict(min_run=8, head_dropout=6),
    # --- combinations that a per-axis sweep cannot reach ---
    dict(min_run=10, mercy=2, phred_sat=25),
    dict(min_run=10, mercy=2, phred_sat=25, mean_q=22.0),
    dict(min_run=14, mercy=5, phred_sat=25, min_len=0.85),
    dict(min_run=8,  mercy=5, phred_sat=30, r2_reverse=True),
    dict(min_run=10, mean_q=22.0, r2_reverse=True, min_len=0.85),
    dict(min_run=12, mercy=2, phred_sat=20, r2_reverse=True, min_len=0.85),
    dict(min_run=16, mercy=8, phred_sat=30, min_len=30),
    dict(min_run=10, mercy=2, phred_sat=25, r2_reverse=True, interior=True),
    dict(min_run=8,  mercy=3, phred_sat=20, min_len=50, mean_q=20.0),
    dict(min_run=12, mercy=3, phred_sat=25, r2_reverse=True, min_overlap=10),
    dict(min_run=10, mercy=5, phred_sat=25, r2_reverse=True, mean_q=22.0),
    dict(min_run=14, mercy=2, phred_sat=20, min_len=0.95, r2_reverse=True),
    dict(min_run=8,  mercy=8, phred_sat=15, r2_reverse=True, interior=True),
    dict(min_run=10, mercy=2, phred_sat=25, min_len=30, r2_reverse=True),
    dict(min_run=12, mercy=5, phred_sat=30, r2_reverse=True, min_overlap=8),
    dict(min_run=10, mercy=3, phred_sat=20, min_len=0.9, r2_reverse=True, head_dropout=4),
]

ROW_FIELDS = ["domain", "run_accession", "variant", "config_id", "min_run", "mercy",
              "phred_sat", "gate_shape", "min_len", "mean_q", "adapter_err",
              "min_overlap", "r2_reverse", "interior", "keep_precision", "keep_recall",
              "trim_mae", "over_trim_rate", "under_trim_rate", "false_keep_rate",
              "drop_accuracy", "kept", "should_keep", "score", "seconds",
              "read_len", "min_len_effective"]

# PXMT record layout (little-endian, v3).
MAGIC = b"PXMT"
STATUS_KEPT = 0


# ---------------------------------------------------------------- PXMT parse
def parse_pxmt(path: str) -> dict[tuple[int, int], tuple[int, int, int]]:
    """-> {(pair_id, mate): (status, t5_total, t3_total)}.

    t5/t3 fold adapter and polyX contributions together, which is what the
    keep window needs: the truth says where biology starts and ends, and both
    trim sources move that boundary.

    fastsim's read IDs are Illumina coordinates
    `@LH00916:<run>:<fc>:<lane>:<tile>:<x>:<y> <mate>:N:0:<idx>` with
    x = pair_id % 10000 and y = pair_id // 10000 (tools/fastsim/src/fastq.rs),
    so pair_id is recovered exactly from the header — that is the join key
    against gt.tsv, whose first column is pair_id.
    """
    out: dict[tuple[int, int], tuple[int, int, int]] = {}
    with open(path, "rb") as fh:
        buf = fh.read()
    if buf[:4] != MAGIC:
        raise ValueError(f"{path}: not PXMT (magic {buf[:4]!r})")
    off = 5
    n = len(buf)
    while off < n:
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len]; off += id_len
        mate = buf[off]; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p, _rlen = struct.unpack_from("<IIIII", buf, off); off += 20
        n_runs = buf[off]; off += 1
        off += n_runs * 8
        out[coord_to_pid(rid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def coord_to_pid(rid: bytes, mate: int) -> tuple[int, int]:
    """Illumina coordinate header -> (pair_id, mate).

    `@LH00916:RUN:FC:LANE:TILE:x:y M:N:0:IDX` with x = pid % 10000,
    y = pid // 10000. Returns (-1, mate) when the header is not that shape,
    which the scorer then treats as unmatched.
    """
    try:
        head = rid.split(b" ", 1)[0].decode("ascii")
        parts = head.lstrip("@").split(":")
        if len(parts) < 7:
            return (-1, mate)
        x, y = int(parts[-2]), int(parts[-1])
        return (y * 10_000 + x, mate)
    except (ValueError, UnicodeDecodeError):
        return (-1, mate)


# ---------------------------------------------------------------- GT parse
def parse_gt(path: str) -> list[tuple]:
    """-> list of (pid, insert, r1_h, r1_t, r2_h, r2_t, full_flag, is_dup)."""
    rows = []
    with open(path, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        header = next(rd)
        # PAIRED: pair_id insert r1_h r1_t r2_h r2_t r1_adapter r2_adapter full_flag is_dup
        i_ins, i_h1, i_t1, i_h2, i_t2, i_full, i_dup = 1, 2, 3, 4, 5, 8, 9
        if header[0] != "pair_id":
            raise ValueError(f"{path}: unexpected GT header {header}")
        for r in rd:
            if not r:
                continue
            rows.append((int(r[0]), int(r[i_ins]), int(r[i_h1]), int(r[i_t1]),
                         int(r[i_h2]), int(r[i_t2]), int(r[i_full]), int(r[i_dup])))
    return rows


# ---------------------------------------------------------------- scoring
def score_config(meta: dict, gt: list[tuple], read_len: int) -> dict:
    kept = 0
    should_keep = 0
    false_keep = 0
    correct_window = 0
    trim_err_sum = 0
    trim_err_n = 0
    over_trim = 0
    under_trim = 0

    for (pid, ins, h1, t1, h2, t2, full, dup) in gt:
        should = (full == 0 and dup == 0)
        if should:
            should_keep += 1

        s1 = meta.get((pid, 0))
        s2 = meta.get((pid, 1))
        if s1 is None or s2 is None:
            continue  # a read fastprep never recorded: treat as dropped
        # A pair survives only if BOTH mates are kept (fastprep PE semantics:
        # a dropped mate drops the pair unless singletons are requested, and
        # this sweep requests none).
        if s1[0] != STATUS_KEPT or s2[0] != STATUS_KEPT:
            continue
        if not should:
            # A kept pair the truth says to drop (full artefact or PCR dup).
            # Counted as a false keep: it is scored through drop_accuracy and
            # excluded from keep/window metrics, so recall is never inflated by
            # keeping what should go.
            false_keep += 1
            continue
        kept += 1

        # Observed keep window on each mate.
        lo1, hi1 = s1[1], read_len - s1[2]
        lo2, hi2 = s2[1], read_len - s2[2]
        # Truth: [h, min(I, L-t))
        k1a, k1b = h1, min(ins, read_len - t1)
        k2a, k2b = h2, min(ins, read_len - t2)

        # trim MAE on the 5' boundary of R1 (where the polyG artefact lives)
        trim_err_sum += abs(lo1 - k1a)
        trim_err_n += 1

        # A near-correct window: boundaries within 2 bp on both mates.
        if (abs(lo1 - k1a) <= 2 and abs(hi1 - k1b) <= 2 and
                abs(lo2 - k2a) <= 2 and abs(hi2 - k2b) <= 2):
            correct_window += 1

        bio = max(k1b - k1a, 1)
        obs = max(hi1 - lo1, 0)
        if obs < 0.75 * bio:
            over_trim += 1
        # under-trim: >=10 bp of artefact left at either end of either mate
        if (lo1 - k1a >= 10 or k1b - hi1 >= 10 or
                lo2 - k2a >= 10 or k2b - hi2 >= 10):
            under_trim += 1

    keep_precision = (correct_window / kept) if kept else 0.0
    keep_recall = (kept / should_keep) if should_keep else 0.0
    trim_mae = (trim_err_sum / trim_err_n) if trim_err_n else 0.0
    over_trim_rate = (over_trim / kept) if kept else 0.0
    under_trim_rate = (under_trim / kept) if kept else 0.0
    # false keeps are relative to the whole GT, not to kept: keeping a pair
    # that should be dropped is an error about the population, not about the
    # kept set.
    n_gt = len(gt)
    false_keep_rate = (false_keep / n_gt) if n_gt else 0.0
    # drop accuracy: how close the kept count is to the true keep count.
    # The deviation counts both ways — a pair dropped that should be kept and a
    # pair kept that should be dropped — so a trimmer cannot buy accuracy by
    # dropping everything.
    if should_keep:
        drop_accuracy = max(0.0, 1.0 - (abs(kept - should_keep) + false_keep) / should_keep)
    else:
        drop_accuracy = 1.0
    score = keep_recall * keep_precision - 2.0 * over_trim_rate - 0.5 * false_keep_rate
    return dict(keep_precision=keep_precision, keep_recall=keep_recall, trim_mae=trim_mae,
                over_trim_rate=over_trim_rate, under_trim_rate=under_trim_rate,
                drop_accuracy=drop_accuracy, false_keep_rate=false_keep_rate,
                kept=kept, should_keep=should_keep, score=score)


GZ_LEVEL = 6  # every tool compresses its output at the same level

# ---------------------------------------------------------------- fastprep
def min_len_arg(v: float) -> str:
    """fastprep takes `N` (bases) or `0.6L` (fraction of the read length).

    A fraction below 1.0 means "keep if >= this fraction survives"; the shipped
    default is 0.9L. Integers are passed as bases.
    """
    if isinstance(v, int) or float(v).is_integer():
        return str(int(v))
    return f"{v}L"


def scale_min_len(cfg: dict, read_len: int) -> float:
    """Resolve a config's min_len against the ACTUAL read length.

    A fraction (0.9) is passed through as `0.9L` and fastprep resolves it. An
    integer in the grid is a base count authored against a 150 bp reference, so
    it is converted to the equivalent fraction of read_len — otherwise 30 means
    "keep 20% of the read" at 150 bp but "keep 27%" at 110 bp, and the two read
    regimes would not be comparable.
    """
    ml = cfg["min_len"]
    if isinstance(ml, int) or (isinstance(ml, float) and float(ml).is_integer()):
        if not read_len:
            return float(ml)
        # The scaled fraction must stay in fastprep's accepted range (0, 1].
        # A 36 bp dataset scales the grid's 50 into 1.3889L, which fastprep
        # correctly refuses — and the refusal aborted the whole sweep on that
        # dataset (53/55 scored, 2 configs lost). Clamp to 1.0 ("keep every
        # base"): a min-len floor above the read length is a no-op floor.
        return round(min(1.0, float(ml) / read_len), 4)
    return ml


def build_argv(cfg: dict, work: str, cid: int, read_len: int = 150) -> list[str]:
    a = ["-1", os.path.join(work, "in_R1.fastq.gz"),
         "-2", os.path.join(work, "in_R2.fastq.gz"),
         "-o", os.path.join(work, f"c{cid}_R1.fastq.gz"),
         "-p", os.path.join(work, f"c{cid}_R2.fastq.gz"),
         "--metadata", os.path.join(work, f"c{cid}.pxmt"),
         "--json", os.path.join(work, f"c{cid}.json"),
         "--force", "-q",
         "--min-run", str(cfg["min_run"]),
         "--mercy", str(cfg["mercy"]),
         "--phred-sat", str(cfg["phred_sat"]),
         "--gate-shape", str(cfg["gate_shape"]),
         "--min-len", min_len_arg(scale_min_len(cfg, read_len)),
         "--mean-q", str(cfg["mean_q"]),
         "--adapter-error-rate", str(cfg["adapter_err"]),
         "--min-adapter-overlap", str(cfg["min_overlap"]),
         "--gz-level", str(GZ_LEVEL)]  # equal compression with every other tool
    if cfg.get("r2_reverse"):
        a.append("--r2-scan-reverse")
    if cfg.get("interior"):
        a.append("--interior")
    return a


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc", required=True)
    ap.add_argument("--domain", required=True)
    ap.add_argument("--variant", type=int, required=True)
    ap.add_argument("--fastprep", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--row", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--configs", type=int, default=len(GRID))
    args = ap.parse_args()

    work = os.path.join(args.work, f"sweep.v{args.variant}")
    os.makedirs(work, exist_ok=True)

    # The simulated reads are the sweep's input; link them once.
    for src, dst in ((args.r1, "in_R1.fastq.gz"), (args.r2, "in_R2.fastq.gz")):
        d = os.path.join(work, dst)
        if not os.path.exists(d) or os.path.realpath(d) != os.path.realpath(src):
            try:
                os.unlink(d)
            except FileNotFoundError:
                pass
            os.link(src, d)

    gt = parse_gt(args.gt)
    # read_len comes from the GT's own artefact clamp (the profile's read length).
    read_len = max(max(r[2] + r[3], r[4] + r[5]) for r in gt) if gt else 150
    # The GT clamps artefact lengths to read_len, so recover it from the sim FASTQ.
    import gzip
    with gzip.open(os.path.join(work, "in_R1.fastq.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            read_len = len(line.rstrip("\n"))
            break

    # Sweep is the ONE script here whose append is correct, and it must stay an
    # append: run_campaign.sh launches one sweep process per (dataset, variant)
    # chain CONCURRENTLY against a single shared rows file, so each chain owns a
    # disjoint slice of the table and writes one line per config with O_APPEND.
    # The header is written only by whoever creates the file (first chain in),
    # which is safe because writerow of a single short line is atomic at this
    # size under O_APPEND.
    #
    # The distinct failure this avoids is the one in closed_loop.py/fastp_compare.py:
    # those are SEQUENTIAL generators over the same units, so a re-run must
    # overwrite. Sweep is a parallel partitioner, so a re-run must append. Same
    # file mode, opposite contracts — the caller decides, and neither script
    # should guess.
    new_file = not os.path.exists(args.row)
    if new_file:
        with open(args.row, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(ROW_FIELDS)

    n_done = 0
    for cid, delta in enumerate(GRID[: args.configs]):
        cfg = dict(DEFAULTS)
        cfg.update(delta)
        rowfile = os.path.join(work, f"c{cid}.row")
        if os.path.exists(rowfile):  # resumable
            with open(rowfile) as fh:
                print(fh.read().rstrip("\n"))
            n_done += 1
            continue

        argv = [args.fastprep] + build_argv(cfg, work, cid, read_len) + ["-j", str(args.threads)]
        t0 = time.time()
        p = subprocess.run(argv, capture_output=True, text=True)
        secs = time.time() - t0
        if p.returncode != 0:
            print(f"sweep: config {cid} failed rc={p.returncode}: {p.stderr[-300:]}",
                  file=sys.stderr)
            continue

        meta = parse_pxmt(os.path.join(work, f"c{cid}.pxmt"))
        m = score_config(meta, gt, read_len)

        vals = [args.domain, args.acc, args.variant, cid,
                cfg["min_run"], cfg["mercy"], cfg["phred_sat"], cfg["gate_shape"],
                cfg["min_len"], cfg["mean_q"], cfg["adapter_err"], cfg["min_overlap"],
                int(cfg["r2_reverse"]), int(cfg["interior"]),
                f"{m['keep_precision']:.4f}", f"{m['keep_recall']:.4f}",
                f"{m['trim_mae']:.2f}", f"{m['over_trim_rate']:.4f}",
                f"{m['under_trim_rate']:.4f}", f"{m['false_keep_rate']:.4f}",
                f"{m['drop_accuracy']:.4f}",
                m["kept"], m["should_keep"], f"{m['score']:.4f}", f"{secs:.2f}",
                read_len, min_len_arg(scale_min_len(cfg, read_len))]
        line = "\t".join(str(v) for v in vals)
        with open(rowfile, "w") as fh:
            fh.write(line + "\n")
        # append to the shared rows file too — this is the campaign's real
        # input (aggregate.py reads it); without the append every chain ran
        # and the aggregate saw only the header. O_APPEND makes concurrent
        # chains safe: each row is one write, the chains never share a cursor.
        with open(args.row, "a") as fh:
            fh.write(line + "\n")
        print(line)
        n_done += 1

        # The per-config outputs are bulk; keep only the metadata + json.
        for f in (f"c{cid}_R1.fastq.gz", f"c{cid}_R2.fastq.gz"):
            try:
                os.unlink(os.path.join(work, f))
            except FileNotFoundError:
                pass

    print(f"sweep: {args.acc} v{args.variant}: {n_done}/{args.configs} configs scored "
          f"(grid has {len(GRID)})", file=sys.stderr)
    return 0 if n_done == args.configs else 1


if __name__ == "__main__":
    sys.exit(main())
