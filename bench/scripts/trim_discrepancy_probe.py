#!/usr/bin/env python3
"""Reconstruct the published vs fresh fastprep trim-accuracy rows side by side.

Why this exists: trim-accuracy.tsv (published) scores fastprep from the Phase D
sweep's STORED PXMT while overtrim-vs-fastp.tsv (fresh) re-runs fastprep on the
current sim. The two disagree (66.50% vs 75.59% exact at tol 0) while fastp
reproduces to the digit in both. This probe, per dataset, reproduces both rows
from their sources and then diffs the two PXMTs per read, to localise which
input changed between the two measurements.

Per dataset it reports:
  A. stored sweep PXMT (sweep.v0/c<cid>.pxmt) scored vs CURRENT GT
     -> must equal the published trim-accuracy.tsv row
  B. fresh run of the CURRENT binary with the fresh harness's argv
     -> must equal the overtrim-vs-fastp.tsv fastprep-hd4 row
  C. fresh run with the SWEEP's argv semantics (scale_min_len-scaled --min-len)
     -> what the published row would be if re-measured consistently
  per-read diffs between A's PXMT and B's/C's PXMT, GT marginals, and a
  content check (are the bases a stale cut removed even present in the
  current read?).

Usage: trim_discrepancy_probe.py D [D ...]   with D = domain/acc
"""
import collections, csv, gzip, os, struct, subprocess, sys

W = "/mnt/Fast/fastsuite-scratch/work"
R = "/mnt/Fast/fastsuite-scratch/results"
FP = "/home/raphael/FastSuite/tools/fastprep/target/release/fastprep"
TOL0, TOL2 = 0, 2


# ------------------------------------------------------------------ readers
def gt_rows(path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["pair_id"])] = (int(r["r1_h"]), int(r["r1_t"]),
                                      int(r["r2_h"]), int(r["r2_t"]))
    return out


def read_fastq(path):
    """-> {(pair_id, mate): (core, seq)}; header key = x + 10000*y."""
    out = {}
    core = pid = mate = None
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                parts = line[1:].split()
                f = parts[0].split(":")
                pid = int(f[5]) + 10000 * int(f[6])
                mate = int(parts[1][0])
                core = parts[0]
            elif i % 4 == 1:
                out[(pid, mate)] = (core, line.rstrip("\n"))
    return out


def parse_pxmt(path):
    """-> {(read_id_core, mate1based): (status, t5, t3, read_len)}."""
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
        (rl,) = struct.unpack_from("<I", buf, off); off += 4
        nrun = buf[off]; off += 1
        off += 8 * nrun
        out[(rid, mate)] = (status, t5a + t5p, t3a + t3p, rl)
    return out


# ------------------------------------------------------------------ scoring
def score(pxmt, coremap, gt, tol):
    c = collections.Counter(exact=0, overtrim=0, undertrim=0,
                            unaligned=0, dropped=0, unpaired=0)
    for (rid, mate), (status, t5, t3, _rl) in pxmt.items():
        key = coremap.get((rid, mate))
        if key is None or key[0] not in gt:
            c["unpaired"] += 1
            continue
        row = gt[key[0]]
        h, t = (row[0], row[1]) if mate == 1 else (row[2], row[3])
        if status != 0:
            c["dropped"] += 1
            continue
        if abs(t5 - h) <= tol and abs(t3 - t) <= tol:
            c["exact"] += 1
        elif t5 > h or t3 > t:
            c["overtrim"] += 1
        else:
            c["undertrim"] += 1
    return c


def bucket_str(c, tol):
    scored = c["exact"] + c["overtrim"] + c["undertrim"]
    return (f"tol{tol} exact={c['exact']} over={c['overtrim']} "
            f"under={c['undertrim']} dropped={c['dropped']} scored={scored}")


# -------------------------------------------------------------- argv builds
def _minlen_fresh(v):
    f = float(v)
    return f"{v}L" if f < 1.0 else str(int(f))


def scale_min_len_sweep(v, read_len):
    """sweep.py scale_min_len + min_len_arg: integers become a fraction of the
    read length, clamped to <= 1.0; min_len_arg emits `N` for integers and
    `0.5556L` for fractions."""
    ml = float(v)
    if ml.is_integer():
        if not read_len:
            return str(int(ml))
        scaled = round(min(1.0, ml / read_len), 4)
        return str(int(scaled)) if float(scaled).is_integer() else f"{scaled}L"
    return f"{v}L"


CONF = {
    "min_run": ("--min-run",), "mercy": ("--mercy",),
    "phred_sat": ("--phred-sat",), "gate_shape": ("--gate-shape",),
    "mean_q": ("--mean-q",), "adapter_err": ("--adapter-error-rate",),
    "min_overlap": ("--min-adapter-overlap",),
}


def base_flags(w):
    a = ["--force", "-q", "--gz-level", "6"]
    for key, (flag,) in CONF.items():
        if w.get(key) not in (None, ""):
            a += [flag, str(w[key])]
    if w.get("r2_reverse") in ("True", "true", "1"):
        a.append("--r2-scan-reverse")
    if w.get("interior") in ("True", "true", "1"):
        a.append("--interior")
    return a


def argv_fresh_harness(w, s1, s2, outdir, tag, hd):
    """Exactly what overtrim_vs_fastp.py emits."""
    a = base_flags(w)
    if w.get("min_len") not in (None, ""):
        a += ["--min-len", _minlen_fresh(w["min_len"])]
    a += ["--head-dropout", str(hd)]
    return [FP, "-1", s1, "-2", s2, "-o", f"{outdir}/{tag}_R1.fq.gz",
            "-p", f"{outdir}/{tag}_R2.fq.gz", "--metadata", f"{outdir}/{tag}.pxmt"] + a


def argv_sweep(w, s1, s2, outdir, tag, read_len):
    """What sweep.py's build_argv emitted at campaign time (no --head-dropout:
    the campaign's build_argv did not emit the flag, fastprep's CLI default
    applied). --min-len is scale_min_len'd."""
    a = base_flags(w)
    if w.get("min_len") not in (None, ""):
        a += ["--min-len", scale_min_len_sweep(w["min_len"], read_len)]
    return [FP, "-1", s1, "-2", s2, "-o", f"{outdir}/{tag}_R1.fq.gz",
            "-p", f"{outdir}/{tag}_R2.fq.gz", "--metadata", f"{outdir}/{tag}.pxmt"] + a


# ------------------------------------------------------------------- probe
def probe(dom_acc, outdir="/tmp/trimprobe"):
    dom, acc = dom_acc.split("/")
    v0 = f"{W}/{dom}/{acc}/v0"
    gtp, s1, s2 = f"{v0}/sim.v0.gt.tsv", f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz"
    gt = gt_rows(gtp)
    fq = read_fastq(s1)
    fq.update(read_fastq(s2))
    coremap = {(core, m): k for k, (core, _s) in fq.items() for m in [k[1]]}

    w = None
    with open(f"{R}/winners.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["domain"] == dom and r["run_accession"] == acc:
                w = r
    assert w, f"no winner row for {dom_acc}"
    cid = w["config_id"]

    # read length: sweep.py uses the FIRST read of in_R1 (fastsim emits a fixed
    # read length per dataset)
    read_len = len(next(iter(fq.values()))[1])

    print(f"== {dom_acc}  winner config {cid}  read_len={read_len}  "
          f"min_len={w.get('min_len')} min_len_effective={w.get('min_len_effective')}")
    print(f"   mtime sim R1: {os.path.getmtime(s1):.0f}  "
          f"stored pxmt: {os.path.getmtime(f'{v0}/sweep.v0/c{cid}.pxmt'):.0f}")

    # GT marginals
    marg = collections.Counter()
    for pid, (h1, t1, h2, t2) in gt.items():
        for h, t in ((h1, t1), (h2, t2)):
            marg["mates"] += 1
            if h == 0 and t == 0:
                marg["00"] += 1
            if h > 0:
                marg["h>0"] += 1
            if t > 0:
                marg["t>0"] += 1
    print(f"   GT: mates={marg['mates']} (h,t)=(0,0): {marg['00']} "
          f"({100*marg['00']/marg['mates']:.1f}%)  h>0: {marg['h>0']}  t>0: {marg['t>0']}")

    # A: stored sweep PXMT vs current GT
    stale = parse_pxmt(f"{v0}/sweep.v0/c{cid}.pxmt")
    print(f"   A stored pxmt records: {len(stale)}")
    for tol in (TOL0, TOL2):
        print("   A " + bucket_str(score(stale, coremap, gt, tol), tol))

    os.makedirs(f"{outdir}/{dom}_{acc}", exist_ok=True)
    d = f"{outdir}/{dom}_{acc}"
    hd = max(1, min(4, int(w["min_run"]) - 1))  # overtrim_vs_fastp.eff_hd, head_dropout absent

    # B: fresh harness argv
    r = subprocess.run(argv_fresh_harness(w, s1, s2, d, "fresh", hd),
                       capture_output=True)
    print(f"   B rc={r.returncode}  min-len flag: "
          f"{argv_fresh_harness(w, s1, s2, d, 'fresh', hd)}".split("--head-dropout")[0][-30:])
    fresh = parse_pxmt(f"{d}/fresh.pxmt")

    # C: sweep argv semantics
    r = subprocess.run(argv_sweep(w, s1, s2, d, "sweepargv", read_len),
                       capture_output=True)
    sweeppx = parse_pxmt(f"{d}/sweepargv.pxmt")

    for tag, px in (("A-stored", stale), ("B-fresh", fresh), ("C-sweepargv", sweeppx)):
        for tol in (TOL0, TOL2):
            c = score(px, coremap, gt, tol)
            print(f"   {tag} " + bucket_str(c, tol))
            with open("/tmp/trimprobe/summary.tsv", "a", newline="") as fh:
                csv.writer(fh, delimiter="\t", lineterminator="\n").writerow(
                    [dom, acc, cid, w.get("min_len"), w.get("min_len_effective"),
                     tag, tol, c["exact"], c["overtrim"], c["undertrim"],
                     c["dropped"],
                     c["exact"] + c["overtrim"] + c["undertrim"]])

    # per-read diffs: stale vs fresh
    keys = set(stale) | set(fresh)
    same = diff_cut = diff_status = only_stale = only_fresh = 0
    deltas = collections.Counter()
    content_ok = content_tot = 0
    for k in keys:
        a_, b_ = stale.get(k), fresh.get(k)
        if a_ is None:
            only_fresh += 1
            continue
        if b_ is None:
            only_stale += 1
            continue
        if (a_[0], a_[1], a_[2]) == (b_[0], b_[1], b_[2]):
            same += 1
            continue
        if a_[0] != b_[0]:
            diff_status += 1
        else:
            diff_cut += 1
            if a_[0] == 0 and b_[0] == 0:
                deltas[(a_[1] - b_[1], a_[2] - b_[2])] += 1
        # content check: did the stale cut remove bases that exist in the
        # current read? t3 cut > 0 -> the current read's last t3 bases should
        # be G if the cut is content-consistent.
        if a_[0] == 0 and a_[2] > 0:
            key = coremap.get(k)
            seq = fq.get(key, (None, ""))[1]
            if seq:
                content_tot += 1
                if seq[len(seq) - a_[2]:].count("G") == a_[2]:
                    content_ok += 1
    print(f"   diff stale-vs-fresh: same={same} cut-diff={diff_cut} "
          f"status-diff={diff_status} only-stale={only_stale} only-fresh={only_fresh}")
    top = deltas.most_common(6)
    print(f"   (dt5,dt3) stale-fresh top: {top}")

    # concrete examples where the cuts differ (generation-shift signature)
    shown = 0
    for k in keys:
        a_, b_ = stale.get(k), fresh.get(k)
        if a_ and b_ and (a_[0], a_[1], a_[2]) != (b_[0], b_[1], b_[2]) and shown < 5:
            key = coremap.get(k)
            row = gt.get(key[0])
            hh, tt = ((row[0], row[1]) if key[1] == 1 else (row[2], row[3])) if row else (None, None)
            seq = fq.get(key, (None, ""))[1]
            g5 = len(seq) - len(seq.lstrip("G"))
            print(f"   ex {key} GT=({hh},{tt}) stale=(st{a_[0]},{a_[1]},{a_[2]}) "
                  f"fresh=(st{b_[0]},{b_[1]},{b_[2]}) leading_G_now={g5}")
            shown += 1

    # is the stale cut consistent with the CURRENT read's own head run?
    for name, px in (("stale", stale), ("fresh", fresh)):
        ok = tot = 0
        for k, (st, t5, t3, _rl) in px.items():
            key = coremap.get(k)
            seq = fq.get(key, (None, ""))[1]
            if st == 0 and t5 > 0 and seq:
                tot += 1
                if seq[:t5].count("G") == t5:
                    ok += 1
        print(f"   {name} t5>0 cuts removing only Gs in the CURRENT read: {ok}/{tot}")
    if content_tot:
        print(f"   stale t3>0 cuts whose removed bases are all G in the CURRENT "
              f"read: {content_ok}/{content_tot}")

    # read-length consistency: PXMT record length vs current read length
    rl_mismatch = 0
    for k, (_s, _t5, _t3, rl) in stale.items():
        key = coremap.get(k)
        if key in fq and rl != len(fq[key][1]):
            rl_mismatch += 1
    print(f"   stale PXMT read_len != current read length: {rl_mismatch}/{len(stale)}")


def main():
    ds = sys.argv[1:] or ["bacterial_isolate/DRR014283"]
    for d in ds:
        probe(d)


if __name__ == "__main__":
    main()
