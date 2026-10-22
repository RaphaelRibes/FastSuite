#!/usr/bin/env python3
"""Prove the 0.4.1 default is the measured configuration, on a dataset where the
knob actually bites.

Runs the same dataset/winnner config three ways -- no flag, `--head-dropout 4`,
`--head-dropout 0` -- and compares the PXMT streams byte for byte. The claim
under test is that a bare invocation IS the 7:1 arm, so its PXMT must equal the
explicit `4` run and differ from the `0` run. A dataset where the two values
produce identical output proves nothing, so this asserts the discriminating
case first.
"""
import csv, hashlib, os, subprocess, sys

FP = "/home/raphael/FastSuite/tools/fastprep/target/release/fastprep"
W = "/mnt/Fast/fastsuite-scratch/work-gt2026"
R = "/mnt/Fast/fastsuite-scratch/results"
DOM, ACC = "bacterial_isolate", "DRR014197"
OUT = "/tmp/hdflip"

CONF = {
    "min_run": ("--min-run", str), "mercy": ("--mercy", str),
    "phred_sat": ("--phred-sat", str), "gate_shape": ("--gate-shape", str),
    "mean_q": ("--mean-q", str), "adapter_err": ("--adapter-error-rate", str),
    "min_overlap": ("--min-adapter-overlap", str),
}


def main():
    w = next(r for r in csv.DictReader(open(f"{R}/winners.tsv"), delimiter="\t")
             if r["run_accession"] == ACC and r["domain"] == DOM)
    base = ["--force", "-q", "--gz-level", "6"]
    for key, (flag, conv) in CONF.items():
        if w.get(key) not in (None, ""):
            base += [flag, conv(w[key])]
    # sweep.py records the RESOLVED floor (a fraction like "0.85L"), not the
    # grid's raw integer: passing the raw value runs a different configuration
    eff = w.get("min_len_effective") or w.get("min_len")
    if not eff:
        sys.exit(f"winner {ACC}: no min_len_effective column")
    base += ["--min-len", eff]
    if w.get("r2_reverse") in ("True", "true", "1"):
        base += ["--r2-scan-reverse"]
    if w.get("interior") in ("True", "true", "1"):
        base += ["--interior"]

    v0 = f"{W}/{DOM}/{ACC}/v0"
    s1, s2 = f"{v0}/sim.v0_R1.fastq.gz", f"{v0}/sim.v0_R2.fastq.gz"
    os.makedirs(OUT, exist_ok=True)
    res = {}
    for tag, extra in (("bare", []), ("hd4", ["--head-dropout", "4"]),
                       ("hd0", ["--head-dropout", "0"])):
        px = f"{OUT}/{tag}.pxmt"
        argv = [FP, "-1", s1, "-2", s2, "-o", f"{OUT}/{tag}_R1.fq.gz",
                "-p", f"{OUT}/{tag}_R2.fq.gz", "--metadata", px] + base + extra
        r = subprocess.run(argv, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"{tag}: rc={r.returncode} {r.stderr[-400:]}")
        res[tag] = hashlib.md5(open(px, "rb").read()).hexdigest()
    print(f"dataset {DOM}/{ACC}  config_id {w['config_id']}  "
          f"read_len {open(f'{v0}/sim.v0.gt.tsv').readline().split(chr(9))[3] if True else '?'}")
    for k, v in res.items():
        print(f"  {k:5} pxmt md5 {v}")
    assert res["bare"] == res["hd4"], "a bare invocation is NOT the 4 configuration"
    assert res["bare"] != res["hd0"], "4 and 0 produce identical output here: pick another dataset"
    print("OK: bare == --head-dropout 4 != --head-dropout 0")


if __name__ == "__main__":
    main()
