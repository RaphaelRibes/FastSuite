#!/usr/bin/env python3
"""Benchmark overtrim / undertrim / exact against simulator ground truth.

The simulator knows the exact composition of every read it emits, so the
correct trim is known exactly rather than estimated. For a read

    [h Gs][biology][adapter?][t Gs]      length L

the correct trim is exactly h from the 5' end and t from the 3' end. A read
with h = t = 0 must be left untouched. Whether a trimmer removed the whole
adapter and the whole polyG is therefore a binary, checkable fact.

Per mate, with observed cuts t5 (front) and t3 (back):

    exact     |t5-h| <= tol and |t3-t| <= tol
    overtrim  t5 > h or t3 > t        real biology trimmed away
    undertrim t5 < h or t3 < t        artefact left behind

Scored at tolerance 0 (strict) and tolerance 2 (forgiving); both reported.

TOOLS (gzip level 6, 4 threads each):
    fastprep-winner   the per-dataset winning config from the Phase D sweep,
                      scored from its PXMT metadata, which records the exact
                      trim applied per read
    fastp-default     fastp with no trimming flags
    fastp-tuned       fastp with cut_front/cut_tail/quality/length/complexity

fastp reports no trim positions, so its cuts are recovered by aligning the
trimmed read back onto the original. A pure trim is a substring of the
original, so the alignment is exact.

READ IDENTITY. fastsim emits illumina headers of the form

    @LH00916:1:SIMFLOW1:1:1101:<tile>:<cycle> <mate>:N:0:<tag>

field 5 is the constant 1101, field 6 is the tile (0..9) and field 7 is the
cycle (0..9999), so (tile, cycle, mate) is unique within a mate file and is
used as the read key throughout. Two earlier revisions keyed reads wrongly
-- one as tile+10000*occurrence, one from the wrong colon fields -- and each
produced a phantom "unaligned" bucket by comparing a read against a
different read's ground truth. The key is asserted unique before scoring.

CAVEAT carried into any paper: this measures agreement with the simulator's
artefact model, which is the same model the closed loop tests. It is exact
ground truth for trimming, not wet-lab truth.
"""
import csv, gzip, os, struct, subprocess, collections

FASTP = os.environ.get("FASTP", "/home/raphael/.local/bin/fastp")
WORK = os.environ.get("TRIM_WORK", "/mnt/Fast/fastsuite-scratch/work")
RESULTS = os.environ.get("TRIM_RESULTS", "/mnt/Fast/fastsuite-scratch/results")
OUTROW = f"{RESULTS}/trim-accuracy.tsv"

FASTP_TUNED = ["--cut_front", "--cut_front_window_size", "5",
               "--cut_front_mean_quality", "20", "--cut_tail",
               "--cut_tail_window_size", "5", "--cut_tail_mean_quality", "20",
               "--qualified_quality_phred", "20", "--length_required", "30",
               "--low_complexity_filter"]

TOLS = (0, 2)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


# ---------------------------------------------------------------- primitives
def header_key(line):
    """(tile, cycle, mate) from a fastsim header line."""
    core = line[1:].split()
    fields = core[0].split(":")
    if len(fields) < 7:
        raise SystemExit(f"unexpected header: {line.rstrip()}")
    tile = int(fields[5])
    cyc = int(fields[6])
    mate = int(core[1][0])
    return (tile, cyc, mate)


def read_fastq(path):
    """-> [(key, sequence)] in file order; key = (tile, cycle, mate)."""
    out = []
    if not os.path.exists(path):
        return out
    seen = collections.Counter()
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                key = header_key(line)
                seen[key] += 1
            elif i % 4 == 1:
                out.append((key, line.rstrip("\n")))
    dup = [k for k, v in seen.items() if v > 1]
    if dup:
        raise SystemExit(f"{path}: {len(dup)} duplicate keys, first {dup[0]}")
    return out


def gt_map(v0):
    """-> {(tile, cycle, mate): (h, t)} keyed from the FASTQ header.

    The GT file lists one row per emitted read in emission order, and the
    FASTQ files are in that same order, so they are replayed in lockstep.
    """
    rows = []
    with open(f"{v0}/sim.v0.gt.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append((int(r["r1_h"]), int(r["r1_t"]),
                         int(r["r2_h"]), int(r["r2_t"])))
    out = {}
    for path in (f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz"):
        reads = read_fastq(path)
        if len(reads) != len(rows):
            raise SystemExit(f"{path}: {len(reads)} reads vs {len(rows)} GT rows")
        for (key, _), (h1, t1, h2, t2) in zip(reads, rows):
            tile, cyc, mate = key
            out[key] = (h1, t1) if mate == 1 else (h2, t2)
    return out


def parse_pxmt(path):
    """-> {(read_id_core, mate): (status, t5, t3)}; status 0 = kept."""
    out = {}
    try:
        buf = open(path, "rb").read()
    except FileNotFoundError:
        return out
    if buf[:4] != b"PXMT":
        return out
    off, n = 5, len(buf)
    while off + 2 <= n:
        (id_len,) = struct.unpack_from("<H", buf, off); off += 2
        rid = buf[off:off + id_len].decode("utf-8", "replace").split()[0]; off += id_len
        if off + 22 > n:
            break
        # PXMT encodes the mate 0-based (Mate::R1 => 0u8); the FASTQ suffix is
        # 1-based, so normalise to 1-based here to match header_key.
        mate = buf[off] + 1; off += 1
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        (_rl,) = struct.unpack_from("<I", buf, off); off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out[(rid, mate)] = (status, t5a + t5p, t3a + t3p)
    return out


def pxmt_core_map(paths):
    """-> {(read_id_core, mate): (tile, cycle, mate)} for PXMT key lookup."""
    out = {}
    for path in paths:
        if not os.path.exists(path):
            continue
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    core = line[1:].split()
                    out[(core[0], mate_of(line))] = header_key(line)
    return out


def mate_of(line):
    core = line[1:].split()
    return int(core[1][0])


def infer_cuts(orig, trimmed):
    """(t5, t3) for a pure trim, or None when the read is not a substring."""
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


def classify(t5, t3, h, t, tol):
    if t5 is None:
        return "unaligned"
    if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
        return "exact"
    if t5 > h or t3 > t:
        return "overtrim"
    return "undertrim"


def run(argv, log):
    with open(log, "w") as fh:
        return subprocess.run(argv, stdout=fh, stderr=subprocess.STDOUT).returncode


def buckets():
    return collections.Counter(exact=0, overtrim=0, undertrim=0,
                               unaligned=0, dropped_or_missing=0)


# ---------------------------------------------------------------- scorers
def score_inferred(orig, trimmed, gt, tol):
    """fastp: recover its cuts by aligning the trimmed read onto the original."""
    tmap = {k: s for k, s in trimmed}
    c = buckets()
    for key, seq in orig:
        if key not in gt:
            continue
        h, t = gt[key]
        tr = tmap.get(key)
        if tr is None:
            c["dropped_or_missing"] += 1
            continue
        cuts = infer_cuts(seq, tr)
        if cuts is None:
            c["unaligned"] += 1
            continue
        t5, t3 = cuts
        c[classify(t5, t3, h, t, tol)] += 1
    return c


def score_pxmt(pxmt, coremap, gt, tol):
    """fastprep: exact trim positions from its own metadata."""
    c = buckets()
    for (rid, mate), (status, t5, t3) in pxmt.items():
        key = coremap.get((rid, mate))
        if key is None or key not in gt:
            continue
        h, t = gt[key]
        if status != 0:
            c["dropped_or_missing"] += 1
            continue
        c[classify(t5, t3, h, t, tol)] += 1
    return c


def main():
    datasets = []
    with open(f"{RESULTS}/winners.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            datasets.append((r["domain"], r["run_accession"], r["config_id"]))
    print(f"trim-accuracy: {len(datasets)} datasets", flush=True)

    # always start from a clean file: a previous partial run may be on disk
    if os.path.exists(OUTROW):
        os.remove(OUTROW)
    with open(OUTROW, "w", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(
            ["domain", "run_accession", "tool", "tol", "exact", "overtrim",
             "undertrim", "unaligned", "dropped_or_missing", "mates_scored"])

    def emit(domain, acc, tool, tol, c):
        with open(OUTROW, "a", newline="") as fh:
            csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(
                [domain, acc, tool, tol, c["exact"], c["overtrim"],
                 c["undertrim"], c["unaligned"], c["dropped_or_missing"],
                 sum(c.values())])

    for i, (domain, acc, cid) in enumerate(datasets, 1):
        v0 = f"{WORK}/{domain}/{acc}/v0"
        if not os.path.exists(f"{v0}/sim.v0.gt.tsv"):
            print(f"  {acc}: no GT", flush=True); continue
        gt = gt_map(v0)
        sim1, sim2 = f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz"
        orig = read_fastq(sim1) + read_fastq(sim2)
        d = f"{RESULTS}/trimwork/{domain}_{acc}"
        os.makedirs(d, exist_ok=True)

        # ---- fastprep: the sweep already ran every config, so the winner's
        # metadata is on disk at sweep.v0/c<cid>.pxmt.
        px = parse_pxmt(f"{v0}/sweep.v0/c{cid}.pxmt")
        if px:
            coremap = pxmt_core_map([sim1, sim2])
            for tol in TOLS:
                emit(domain, acc, "fastprep-winner", tol,
                     score_pxmt(px, coremap, gt, tol))
        else:
            print(f"  {acc}: no PXMT for config {cid}", flush=True)

        # ---- fastp: run on the same simulated input
        for tool, extra in (("fastp-default", []), ("fastp-tuned", FASTP_TUNED)):
            o1, o2 = f"{d}/{tool}_R1.fastq.gz", f"{d}/{tool}_R2.fastq.gz"
            if not os.path.exists(o1):
                rc = run([FASTP, "-i", sim1, "-I", sim2, "-o", o1, "-O", o2,
                          "-w", "4", "-z", "6"] + extra, f"{d}/{tool}.log")
                if rc != 0:
                    print(f"  {acc} {tool}: rc={rc}", flush=True); continue
            trimmed = read_fastq(o1) + read_fastq(o2)
            for tol in TOLS:
                emit(domain, acc, tool, tol, score_inferred(orig, trimmed, gt, tol))
        print(f"  [{i}/{len(datasets)}] {acc}: done", flush=True)

    print(f"wrote {OUTROW}")


if __name__ == "__main__":
    main()
