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

READ IDENTITY -- POSITIONAL, AND WHY. fastsim's headers are

    @LH00916:1:SIMFLOW1:1:1101:<tile>:<cycle> <mate>:N:0:<tag>

where <tile> takes only 10 values (0..9) and each repeats 10000 times per
mate file. No header field is a stable read identifier: (tile, cycle, mate)
collides 10000 times per file. Reads are therefore paired POSITIONALLY --
both fastp and fastprep preserve input order and only drop reads, so the
i-th surviving trimmed read belongs to the i-th surviving original. The GT
file is in the same emission order and is replayed positionally.

Three earlier revisions keyed reads by header fields or by
tile+10000*occurrence; each mismatched reads against the wrong ground truth
and produced a phantom "unaligned" bucket of ~5% (fastp) that was pure
scorer artefact. The stored trimwork/ fastp outputs from those runs were
also found to be internally inconsistent (most trimmed reads were not
substrings of any original), so this scorer re-runs fastp from scratch and
scores only freshly generated output.

Pairing is a single forward walk: for each trimmed read, advance the
original pointer until an original contains it; originals skipped on the way
are reads the tool dropped. A trimmed read that matches no remaining
original is counted unaligned and never advances the pointer, which keeps
the walk from desynchronising.

CAVEAT carried into any paper: this measures agreement with the simulator's
artefact model, which is the same model the closed loop tests. It is exact
ground truth for trimming, not wet-lab truth.
"""
import csv, gzip, os, struct, subprocess, collections

FASTP = os.environ.get("FASTP", "/home/raphael/.local/bin/fastp")
WORK = os.environ.get("TRIM_WORK", "/mnt/Fast/fastsuite-scratch/work")
RESULTS = os.environ.get("TRIM_RESULTS", "/mnt/Fast/fastsuite-scratch/results")
OUTROW = f"{RESULTS}/trim-accuracy.tsv"
TRIMWORK = f"{RESULTS}/trimwork2"

FASTP_TUNED = ["--cut_front", "--cut_front_window_size", "5",
               "--cut_front_mean_quality", "20", "--cut_tail",
               "--cut_tail_window_size", "5", "--cut_tail_mean_quality", "20",
               "--qualified_quality_phred", "20", "--length_required", "30",
               "--low_complexity_filter"]

TOLS = (0, 2)
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


# ---------------------------------------------------------------- primitives
def read_fastq(path):
    """-> [sequence] in file order."""
    out = []
    if not os.path.exists(path):
        return out
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out.append(line.rstrip("\n"))
    return out


def read_ids(path):
    """-> [(read_id_core, mate)] in file order, for PXMT pairing."""
    out = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                core = line[1:].split()
                out.append((core[0], int(core[1][0])))
    return out


def gt_rows(v0):
    """-> [(h1,t1,h2,t2)] in emission order."""
    out = []
    with open(f"{v0}/sim.v0.gt.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.append((int(r["r1_h"]), int(r["r1_t"]),
                        int(r["r2_h"]), int(r["r2_t"])))
    return out


def parse_pxmt(path):
    """-> [(rid, mate, status, t5, t3)] in PXMT (input) order."""
    out = []
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
        mate = buf[off] + 1; off += 1   # PXMT mate is 0-based
        status = buf[off]; off += 1
        t5a, t3a, t5p, t3p = struct.unpack_from("<IIII", buf, off); off += 16
        (_rl,) = struct.unpack_from("<I", buf, off); off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out.append((rid, mate, status, t5a + t5p, t3a + t3p))
    return out


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
def score_positional(orig, trimmed, gt, tol):
    """Pair by a single forward walk; see module docstring."""
    c = buckets()
    oi = 0
    matched = 0
    for t in trimmed:
        while oi < len(orig) and t not in orig[oi]:
            oi += 1
        if oi >= len(orig):
            c["unaligned"] += 1
            continue
        cuts = infer_cuts(orig[oi], t)
        h, tt = gt[oi]
        c[classify(cuts[0], cuts[1], h, tt, tol)] += 1
        matched += 1
        oi += 1
    c["dropped_or_missing"] = len(orig) - matched
    return c


def score_pxmt(pxmt, ids, gt, tol):
    """fastprep: exact trim positions, paired positionally per mate."""
    c = buckets()
    by_mate = collections.defaultdict(list)
    for rid, mate, status, t5, t3 in pxmt:
        by_mate[mate].append((status, t5, t3))
    for mate in (1, 2):
        entries = by_mate.get(mate, [])
        file_ids = [i for i in ids if i[1] == mate]
        if len(entries) != len(file_ids):
            raise SystemExit(f"PXMT mate{mate}: {len(entries)} entries vs "
                             f"{len(file_ids)} reads")
        for (status, t5, t3), idx in zip(entries, range(len(file_ids))):
            h, t = gt[idx][0:2] if mate == 1 else gt[idx][2:4]
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

    os.makedirs(TRIMWORK, exist_ok=True)
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
        rows = gt_rows(v0)
        sim1, sim2 = f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz"
        ids = read_ids(sim1) + read_ids(sim2)
        d = f"{TRIMWORK}/{domain}_{acc}"
        os.makedirs(d, exist_ok=True)

        # ---- fastprep: winner config metadata from the sweep
        px = parse_pxmt(f"{v0}/sweep.v0/c{cid}.pxmt")
        if px:
            for tol in TOLS:
                emit(domain, acc, "fastprep-winner", tol,
                     score_pxmt(px, ids, rows, tol))
        else:
            print(f"  {acc}: no PXMT for config {cid}", flush=True)

        # ---- fastp: fresh run, never reuse stored output
        for tool, extra in (("fastp-default", []), ("fastp-tuned", FASTP_TUNED)):
            o1, o2 = f"{d}/{tool}_R1.fastq.gz", f"{d}/{tool}_R2.fastq.gz"
            for f in (o1, o2):
                if os.path.exists(f):
                    os.remove(f)
            rc = run([FASTP, "-i", sim1, "-I", sim2, "-o", o1, "-O", o2,
                      "-w", "4", "-z", "6"] + extra, f"{d}/{tool}.log")
            if rc != 0:
                print(f"  {acc} {tool}: rc={rc}", flush=True); continue
            for tol in TOLS:
                c = buckets()
                for path_o, path_t, mate in ((sim1, o1, 1), (sim2, o2, 2)):
                    orig = read_fastq(path_o)
                    trimmed = read_fastq(path_t)
                    gt = [r[0:2] if mate == 1 else r[2:4] for r in rows]
                    c.update(score_positional(orig, trimmed, gt, tol))
                emit(domain, acc, tool, tol, c)
        print(f"  [{i}/{len(datasets)}] {acc}: done", flush=True)

    print(f"wrote {OUTROW}")


if __name__ == "__main__":
    main()
