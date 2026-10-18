#!/usr/bin/env python3
"""Mechanism diagnostic: what does each opt-in flag change, per mate, on ONE
dataset? Runs fastprep three times (baseline / --min-run-3p 6 /
--trim-terminal-n / both) on the same dataset and joins each PXMT cut to the
sim ground truth, so the transition can be attributed rather than guessed.

Usage: diag_flags.py <domain> <acc>
"""
import collections, csv, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import overtrim_vs_fastp as H
import trim_report as T

dom = sys.argv[1]
acc = sys.argv[2]
R = H.R
w = next(x for x in csv.DictReader(open(f"{R}/winners.tsv"), delimiter="\t")
         if x["domain"] == dom and x["run_accession"] == acc)
v0 = f"{H.W}/{dom}/{acc}/v0"
gtp, s1, s2 = (f"{v0}/sim.v0.gt.tsv", f"{v0}/sim.v0_R1.fastq.gz",
               f"{v0}/sim.v0_R2.fastq.gz")
gt = T.gt_full(gtp)
coremap = H.core_to_pid([s1, s2])
raw = {(pid, m): seq for pid, m, seq in
       (H.read_fastq_keyed(s1) + H.read_fastq_keyed(s2))}

base = ["--force", "-q", "--gz-level", "6"]
for key, (flag, conv) in H.CONF.items():
    if w.get(key) not in (None, ""):
        base += [flag, conv(w[key])]
if w.get("min_len") not in (None, ""):
    base += ["--min-len", H._minlen(w)]
if str(w.get("r2_reverse")) in ("True", "true", "1"):
    base += ["--r2-scan-reverse"]
if str(w.get("interior")) in ("True", "true", "1"):
    base += ["--interior"]

VARIANTS = [("base", []), ("3p6", ["--min-run-3p", "6"]),
            ("tn", ["--trim-terminal-n"]),
            ("both", ["--min-run-3p", "6", "--trim-terminal-n"])]

d = f"/tmp/fpdiag/{dom}_{acc}"
os.makedirs(d, exist_ok=True)
res = {}
for tag, extra in VARIANTS:
    px = f"{d}/{tag}.pxmt"
    r = subprocess.run([H.FP, "-1", s1, "-2", s2, "-o", f"{d}/{tag}_1.fq.gz",
                        "-p", f"{d}/{tag}_2.fq.gz", "--metadata", px]
                       + base + ["--head-dropout", "0"] + extra,
                       capture_output=True)
    if r.returncode != 0:
        sys.exit(f"{tag}: rc={r.returncode} {r.stderr[-400:]}")
    res[tag] = H.parse_pxmt(px)
    n = len(res[tag])
    print(f"{tag:6} mates_in_pxmt={n}", file=sys.stderr)

TOL = 0
for vtag in ("3p6", "tn", "both"):
    print(f"\n################ variant {vtag} vs baseline "
          f"({dom}/{acc}, tol={TOL}) ################")
    trans = collections.Counter()
    changed3 = collections.Counter()
    changed5 = collections.Counter()
    mag = collections.Counter()
    n_run_at_end = collections.Counter()
    for key in set(res["base"]) & set(res[vtag]):
        rid, mate = key
        pid = coremap.get(key)
        if pid is None or pid not in gt:
            continue
        g = gt[pid]
        h, t = (g[0], g[1]) if mate == 1 else (g[2], g[3])
        ad = g[4] if mate == 1 else g[5]
        nn = g[6] if mate == 1 else g[7]
        b_stat, b5, b3 = res["base"][key]
        v_stat, v5, v3 = res[vtag][key]
        cb = H.newb(); H.bucket(cb, b5, b3, h, t, TOL)
        cv = H.newb(); H.bucket(cv, v5, v3, h, t, TOL)
        bc = [k for k, n in cb.items() if n]
        vc = [k for k, n in cv.items() if n]
        trans[(bc[0], vc[0])] += 1
        if (v5, v3) == (b5, b3):
            continue
        seq = raw.get((pid, mate), "")
        mag[(v5 - b5) + (v3 - b3)] += 1
        if v3 != b3:
            changed3[(v3 > b3, b3 == t, v3 == t, nn > 0,
                      T.trailing_g(seq), seq[-1] if seq else "?")] += 1
        if v5 != b5:
            changed5[(v5 > b5, b5 == h, v5 == h, nn > 0,
                      T.leading_g(seq), seq[0] if seq else "?")] += 1
    print("transition matrix (baseline class -> variant class):")
    for (a, b), n in sorted(trans.items(), key=lambda kv: -kv[1]):
        print(f"   {a:10} -> {b:10} {n:9d}")
    print(f"mates whose (t5,t3) moved: {sum(mag.values())}")
    print("distribution of the trim delta (t5+t3 change, base=0):")
    for k, n in sorted(mag.items())[:12]:
        print(f"   delta {k:+3d} : {n:9d}")
    print("3' moves by (increased?, base already at GT, now exact, GT n>0, "
          "raw trailing-G run, terminal base):")
    for k, n in sorted(changed3.items(), key=lambda kv: -kv[1])[:12]:
        print(f"   incr={str(k[0]):5} base_at_gt={str(k[1]):5} "
              f"now_exact={str(k[2]):5} gt_n={str(k[3]):5} "
              f"raw_trailG_run={k[4]:3} base={k[5]} : {n:9d}")
    print("5' moves by (increased?, base already at GT, now exact, GT n>0, "
          "raw leading-G run, terminal base):")
    for k, n in sorted(changed5.items(), key=lambda kv: -kv[1])[:12]:
        print(f"   incr={str(k[0]):5} base_at_gt={str(k[1]):5} "
              f"now_exact={str(k[2]):5} gt_n={str(k[3]):5} "
              f"raw_leadG_run={k[4]:3} base={k[5]} : {n:9d}")

# ---- raw N content of the dataset, and how the GT accounts for it ----
print(f"\n################ raw N content ({dom}/{acc}) ################")
nstat = collections.Counter()
leadn = collections.Counter()
trailn = collections.Counter()
for key, seq in raw.items():
    pid, mate = key
    if pid not in gt:
        continue
    g = gt[pid]
    nn = g[6] if mate == 1 else g[7]
    ad = g[4] if mate == 1 else g[5]
    h, t = (g[0], g[1]) if mate == 1 else (g[2], g[3])
    l5, l3 = T.leading_g(seq), T.trailing_g(seq)
    cnt_n = seq.count("N")
    nstat[("hasN", cnt_n > 0), ("gt_n>0", nn > 0), ("gt_ad>0", ad > 0), ("gt_t>0", t > 0)] += 1
    if l5 or l3:
        leadn[(l5, l3, nn, h, t)] += 1
for k, n in sorted(nstat.items(), key=lambda kv: -kv[1])[:20]:
    print(f"   {k} : {n}")
print("   (l5,l3, gt_n, gt_h, gt_t) for reads with a leading/trailing G run:")
for k, n in sorted(leadn.items(), key=lambda kv: -kv[1])[:15]:
    print(f"   {k} : {n}")

# ---- how many mates have a raw TERMINAL N run, and what the GT wants there ----
print("\nterminal-N census (raw read starts/ends with N):")
tn = collections.Counter()
for key, seq in raw.items():
    pid, mate = key
    if pid not in gt:
        continue
    g = gt[pid]
    h, t = (g[0], g[1]) if mate == 1 else (g[2], g[3])
    nn = g[6] if mate == 1 else g[7]
    ad = g[4] if mate == 1 else g[5]
    lead = 0
    while lead < len(seq) and seq[lead] == "N":
        lead += 1
    trail = 0
    while trail < len(seq) and seq[-1 - trail] == "N":
        trail += 1
    if lead or trail:
        tn[(("lead", lead) if lead else ("trail", trail), "gt_h" if lead else "gt_t",
            h if lead else t, "gt_n", nn, "gt_ad", ad)] += 1
for k, n in sorted(tn.items(), key=lambda kv: -kv[1])[:25]:
    print(f"   {k} : {n}")
