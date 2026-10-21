#!/usr/bin/env python3
"""One driver for the two-axis undertrim report.

WHY THIS EXISTS
---------------
A trimming change is only a measurement when it is quoted on BOTH axes:
residual adapter (undertrim) and false trims on clean reads (overtrim). The
two-axis discipline lives in `undertrim_ab.py`, but every fix so far has been
measured by remembering to run several probes and compare the numbers by
eye. This driver runs them all and prints one table with a PASS/FAIL against
declared thresholds, so "the fix works" cannot be claimed on one axis.

It is the gate for Phase 1 of the round-3 plan: no default flips until this
script passes.

Usage
-----
    python3 bench/scripts/undertrim_report.py [--fastprep PATH] [--quick]
    python3 bench/scripts/undertrim_report.py --paired-only

Exit code 0 when every threshold holds, 1 otherwise.
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# --------------------------------------------------------------------------
# Thresholds. These are the CONTRACT; change them only with a reason in the
# commit message, and never to make a failing binary pass.
# --------------------------------------------------------------------------
RESIDUAL_SE_MAX = 0.25     # residual / carriers, single-end
RESIDUAL_PE_MAX = 0.25     # residual / carriers, paired
FALSE_CLEAN_MAX = 0.0002   # false trims / clean reads
PE_EXACT_MIN = 0.90       # the PE probe's exact rate at a short tail
GUARD_B_MUST_BE = 0        # arm B (floor 14, t=10) must cut nothing

FASTPREP = "tools/fastprep/target/release/fastprep"


def run_script(script, args, out):
    """Run one probe script; return its parsed JSON (or {} on failure)."""
    cmd = [sys.executable, os.path.join(HERE, script)] + args
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  !! {script} failed rc={r.returncode}")
        print(r.stderr[-1500:])
        return {}
    try:
        with open(out) as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        print(f"  !! cannot read {out}: {e}")
        return {}


def section(title):
    print(f"\n{'=' * 72}\n{title}\n{'=' * 72}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastprep", default=FASTPREP)
    ap.add_argument("--quick", action="store_true",
                    help="smaller corpora (a smoke test, not a measurement)")
    ap.add_argument("--paired-only", action="store_true")
    ap.add_argument("--out", default="bench/results/undertrim-report.json")
    ap.add_argument("--save-baseline", action="store_true",
                    help="write this run as the baseline the next run compares against")
    a = ap.parse_args()

    fp = os.path.abspath(a.fastprep)
    if not os.path.exists(fp):
        print(f"fastprep binary not found: {fp}\n"
              f"build it first: pixi run build")
        return 1
    d = tempfile.mkdtemp(prefix="utrep_")
    n_se = 20000 if a.quick else 200000
    n_pe = 100 if a.quick else 500
    checks = []
    report: dict = {"fastprep": fp, "quick": a.quick, "checks": {}}

    def record(name, value, ok, detail):
        checks.append((name, value, ok, detail))
        report["checks"][name] = dict(value=value, ok=bool(ok), detail=detail)

    # ---------------------------------------------------------------- guard
    if not a.paired_only:
        section("1. The PE guard (is the boundary still discarded by the floor?)")
        g = run_script("undertrim_probe_pe_guard.py",
                       ["--fastprep", fp, "--n", "200" if a.quick else "200",
                        "--out", os.path.join(d, "guard.json")],
                       os.path.join(d, "guard.json"))
        v = g.get("verdict", {})
        print(f"  arm A (t=20, floor 14) both mates cut : {v.get('arm_A_both_mates_cut')}")
        print(f"  arm B (t=10, floor 14) nothing cut     : {v.get('arm_B_nothing_cut')}")
        print(f"  arm C (t=10, floor 10) both mates cut : {v.get('arm_C_both_mates_cut')}")
        # Arm B must cut nothing: it is the floor blocking a pinned boundary.
        record("guard.arm_B_nothing_cut", v.get("arm_B_nothing_cut"),
               v.get("arm_B_nothing_cut") is True,
               "arm B (t=10, floor 14) must cut nothing")
        # The mechanism UT-1 addresses must still be present.
        record("guard.the_floor_is_the_blocker", v.get("the_floor_is_the_blocker"),
               v.get("the_floor_is_the_blocker") is True,
               "the mechanism UT-1 addresses is still present")

    # ------------------------------------------------------------- PE probe
    if not a.paired_only:
        section("2. The PE probe (tail10/tail6 exact rate)")
        pe = run_script("undertrim_probe_pe.py",
                        ["--fastprep", fp, "--n", str(n_pe),
                         "--out", os.path.join(d, "pe.json")],
                        os.path.join(d, "pe.json"))
        cases = pe.get("cases", {})
        # schema: cases -> <case name> -> "<pe-mode>/<mate>" -> {exact, none, partial}
        for name, row in cases.items():
            if not ("tail10" in name or "tail6" in name):
                continue
            # "always/r1" is the PE path on the mate that carries the adapter;
            # "never/r1" is the ordinary matcher alone.
            cell = row.get("always/r1")
            if not cell:
                continue
            ex = cell.get("exact", 0)
            # The PE probe's cells carry no "n"; every case is built with the
            # same --n, so that is the denominator.
            rate = ex / n_pe
            print(f"  {name:24s} always/r1 exact {ex}/{n_pe} = {rate*100:.1f}%")
            record(f"pe.{name}.exact_rate", rate, rate >= PE_EXACT_MIN,
                   f"short-tail exact rate must be >= {PE_EXACT_MIN}")

    # ------------------------------------------------------------ SE two-axis
    if not a.paired_only:
        section("3. The single-end two-axis A/B")
        se = run_script("undertrim_ab.py",
                        ["--fastprep", fp, "--reads", str(n_se),
                         "--out", os.path.join(d, "se.json")],
                        os.path.join(d, "se.json"))
        shipped = se.get("configs", {}).get("shipped (mo=14, eps=0.2)", {})
        if shipped:
            ur = shipped.get("undertrim_rate")
            fr = shipped.get("clean_false_trim_rate")
            print(f"  residual / carrier : {ur*100:.2f}%  (max {RESIDUAL_SE_MAX*100:.0f}%)")
            print(f"  false trims / clean : {fr*100:.3f}%  (max {FALSE_CLEAN_MAX*100:.3f}%)")
            record("se.residual_rate", ur, ur is not None and ur <= RESIDUAL_SE_MAX,
                   f"residual/carrier <= {RESIDUAL_SE_MAX}")
            record("se.clean_false_trim_rate", fr,
                   fr is not None and fr <= FALSE_CLEAN_MAX,
                   f"false trims/clean <= {FALSE_CLEAN_MAX}")

    # ----------------------------------------------------------- PE two-axis
    section("4. The paired two-axis A/B")
    peab = run_script("undertrim_ab.py",
                      ["--fastprep", fp, "--paired", "--reads", str(n_se),
                       "--out", os.path.join(d, "peab.json")],
                      os.path.join(d, "peab.json"))
    shipped = peab.get("configs", {}).get("shipped (mo=14, eps=0.2)", {})
    if shipped:
        ur = shipped.get("undertrim_rate")
        fr = shipped.get("clean_false_trim_rate")
        print(f"  residual / carrier : {ur*100:.2f}%  (max {RESIDUAL_PE_MAX*100:.0f}%)")
        print(f"  false trims / clean : {fr*100:.3f}%  (max {FALSE_CLEAN_MAX*100:.3f}%)")
        record("pe.residual_rate", ur, ur is not None and ur <= RESIDUAL_PE_MAX,
               f"paired residual/carrier <= {RESIDUAL_PE_MAX}")
        record("pe.clean_false_trim_rate", fr,
               fr is not None and fr <= FALSE_CLEAN_MAX,
               f"paired false trims/clean <= {FALSE_CLEAN_MAX}")

    # ---------------------------------------------------------------- verdict
    section("VERDICT")
    failed = 0
    for name, value, ok, detail in checks:
        status = "PASS" if ok else "FAIL"
        if not ok:
            failed += 1
        shown = value if not isinstance(value, float) else f"{value*100:.3f}%"
        print(f"  [{status}] {name:36s} {str(shown):>12s}   {detail}")
    print(f"\n{len(checks) - failed}/{len(checks)} checks passed")

    # ---- baseline comparison ------------------------------------------------
    # A fix is judged against the "before" numbers, not only against the
    # thresholds: a change that improves the residual but worsens the clean
    # axis has still lost, and only the two columns side by side show that.
    base_path = os.path.join(HERE, "..", "results", "undertrim-baseline.json")
    baseline = {}
    if os.path.exists(base_path):
        try:
            with open(base_path) as f:
                baseline = json.load(f)
        except (OSError, json.JSONDecodeError):
            baseline = {}
    if baseline:
        section("AGAINST THE BASELINE")
        print(f"  {'metric':32s} {'before':>12s} {'now':>12s} {'delta':>12s}")
        for name, value, _ok, _detail in checks:
            b = baseline.get("checks", {}).get(name, {}).get("value")
            if b is None or not isinstance(value, (int, float)) \
                    or not isinstance(b, (int, float)):
                continue
            delta = value - b
            # Direction is per metric: a residual rate is better when LOWER,
            # an exact rate is better when HIGHER.
            higher_is_better = "exact_rate" in name
            if abs(delta) < 1e-12:
                arrow = ""
            elif (delta > 0) == higher_is_better:
                arrow = "  better"
            else:
                arrow = "  WORSE"
            print(f"  {name:32s} {b:12.6f} {value:12.6f} {delta:+12.6f}{arrow}")
    else:
        print(f"\n(no baseline at {base_path} -- run with --save-baseline to create one)")

    if a.save_baseline:
        out = os.path.join(HERE, "..", "results", "undertrim-baseline.json")
        with open(out, "w") as f:
            json.dump(report, f, indent=1)
        print("baseline written to", out)

    with open(a.out, "w") as f:
        json.dump(report, f, indent=1)
    print("wrote", a.out)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
