#!/usr/bin/env python3
"""Self-test for the weighted trim score: hand-built decisions, hand-arithmetic.

The scorer's correctness is not "the number looks plausible". A synthetic PXMT
is written byte by byte from the v3 layout (with a GT that says what the right
answer is), every class is planted at least once, and the weighted score is
asserted against the arithmetic done here in the comments. Offline: no corpus,
no binaries, no network.

Fixtures
    fx9     9 pairs covering exact / overtrim / undertrim / mixed_trim /
            drop_correct (twice: once via the GT full flag, once via the
            min-len floor) / drop_wrong / keep_wrong / unaligned
    ds0     4 pairs measured by two rules (base, var) whose deltas are set so
            that r* is a non-trivial 2.0 rather than the 1.0 every
            mate-conserving swap produces
    fxtol   1 pair that is overtrim at --tol 0 and exact at --tol 2
"""
import csv
import os
import struct
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_weighted as sw
import trim_accuracy as ta  # reference semantics for classify()

L = 150
FLOOR = int(0.9 * L)  # 135, fastprep's shipped --min-len 0.9L

GT_HEADER = ("pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tr1_adapter\tr2_adapter"
             "\tfull_flag\tis_dup\tr1_full\tr2_full\n")

# pair_id insert r1_h r1_t r2_h r2_t ad ad full dup r1f r2f
FX9_GT = [
    (0, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),      # exact both mates
    (1, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),      # overtrim both mates
    (2, 300, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0),      # undertrim both mates
    (3, 300, 4, 6, 4, 6, 0, 0, 0, 0, 0, 0),      # mixed: mate1 over, mate2 under
    (4, 300, 100, 100, 100, 100, 0, 0, 1, 0, 1, 1),   # full -> must drop
    (5, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),       # window 80 < floor -> must drop
    (6, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),      # window 150 -> must keep
    (7, 300, 100, 100, 100, 100, 0, 0, 1, 0, 1, 1),   # full but kept
    (8, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),      # unaligned both mates
]

# (pid, mate, status, t5, t3) — mate is 1-based here, converted on write
FX9_DEC = [
    (0, 1, 0, 0, 0), (0, 2, 0, 0, 0),            # exact / exact
    (1, 1, 0, 5, 0), (1, 2, 0, 5, 0),            # 5 > h=0 -> overtrim
    (2, 1, 0, 0, 3), (2, 2, 0, 0, 3),            # 0 < h=5 -> undertrim
    (3, 1, 0, 9, 6),                             # 9 > h=4 -> overtrim
    (3, 2, 0, 4, 0),                             # 0 < t=6 -> undertrim  => mixed
    (4, 1, 1, 0, 0), (4, 2, 1, 0, 0),            # dropped, full      -> drop_correct
    (5, 1, 1, 0, 0), (5, 2, 1, 0, 0),            # dropped, 80 < 135  -> drop_correct
    (6, 1, 1, 0, 0), (6, 2, 1, 0, 0),            # dropped, 150 >= 135 -> drop_wrong
    (7, 1, 0, 0, 0), (7, 2, 0, 0, 0),            # kept, full -> keep_wrong
    (8, 1, 0, 100, 100), (8, 2, 0, 100, 100),    # 200 > L -> unaligned
]

DS0_GT = [
    (0, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),      # b0 exact in base
    (1, 300, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0),      # b1 undertrim in base
    (2, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),      # b2 overtrim in both
    (3, 300, 100, 100, 100, 100, 0, 0, 1, 0, 1, 1),   # b3 drop_correct in both
]
DS0_BASE = [
    (0, 1, 0, 0, 0), (0, 2, 0, 0, 0),            # exact, exact
    (1, 1, 0, 0, 0), (1, 2, 0, 0, 0),            # undertrim, undertrim
    (2, 1, 0, 5, 0), (2, 2, 0, 5, 0),            # overtrim, overtrim
    (3, 1, 1, 0, 0), (3, 2, 1, 0, 0),            # drop_correct, drop_correct
]
DS0_VAR = [
    (0, 1, 1, 0, 0), (0, 2, 1, 0, 0),            # dropped, window 150 -> drop_wrong x2
    (1, 1, 0, 8, 8), (1, 2, 0, 8, 8),            # 8 > h=5 -> overtrim x2 (was undertrim)
    (2, 1, 0, 5, 0), (2, 2, 0, 5, 0),            # overtrim unchanged
    (3, 1, 1, 0, 0), (3, 2, 1, 0, 0),            # drop_correct unchanged
]

FXTOL_GT = [(0, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)]
FXTOL_DEC = [(0, 1, 0, 2, 0), (0, 2, 0, 0, 0)]   # off by exactly 2 bp


# ---------------------------------------------------------------- builders
def header(pid, mate):
    return f"@LH00916:1:SIMFLOW1:1:1101:{pid % 10000}:{pid // 10000} {mate}:N:0:ATCGATCG".encode()


def rec(pid, mate, status, t5=0, t3=0, rlen=L, runs=()):
    """One PXMT v3 record: u16 id_len | id | u8 mate (0-based) | u8 status |
    u32 t5_adapter | u32 t3_adapter | u32 t5_poly | u32 t3_poly | u32 read_len |
    u8 n_runs | n_runs*(u32, u32)."""
    h = header(pid, mate)
    out = struct.pack("<H", len(h)) + h + bytes([mate - 1, status])
    out += struct.pack("<IIIII", t5, 0, 0, t3, rlen)
    out += bytes([len(runs)]) + b"".join(struct.pack("<II", a, b) for a, b in runs)
    return out


def write_pxmt(path, decs, runs_by_pid=None):
    buf = b"PXMT\x03"
    for pid, mate, status, t5, t3 in decs:
        runs = (runs_by_pid or {}).get(pid, ())
        buf += rec(pid, mate, status, t5, t3, runs=runs)
    with open(path, "wb") as fh:
        fh.write(buf)
    return path


def write_gt(path, rows):
    with open(path, "w", newline="") as fh:
        fh.write(GT_HEADER)
        for r in rows:
            fh.write("\t".join(str(v) for v in r) + "\n")
    return path


fails = []


def check(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol if isinstance(want, float) else got == want
    print(f"{'ok  ' if ok else 'FAIL'} {name}: got {got} want {want}")
    if not ok:
        fails.append(name)


# ---------------------------------------------------------------- unit checks
# classify() must agree with the scorer it replaces over the whole small grid
_bad = 0
for tol in (0, 2):
    for t5 in (0, 2, 3, 5, 9):
        for t3 in (0, 2, 3, 5, 6):
            for h in (0, 4, 5):
                for t in (0, 5, 6):
                    if sw.classify_cuts(t5, t3, h, t, tol) != ta.classify(t5, t3, h, t, tol):
                        _bad += 1
check("classify_cuts == trim_accuracy.classify over 900 cells", _bad, 0)
check("classify_cuts(None,..) is unaligned", sw.classify_cuts(None, 0, 0, 0, 0), "unaligned")

# break-even sign cases
check("break_even(+4,-2) r*", sw.break_even(4, -2, 7.0)[0], 2.0)
check("break_even(+4,-2) verdict", sw.break_even(4, -2, 7.0)[1], "above")
check("break_even(+4,-2) better at 7", sw.break_even(4, -2, 7.0)[2], True)
check("break_even(+4,-2) at r=1 worse", sw.break_even(4, -2, 1.0)[2], False)
check("break_even(-2,+4) r*", sw.break_even(-2, 4, 7.0)[0], 0.5)
check("break_even(-2,+4) verdict", sw.break_even(-2, 4, 7.0)[1], "below")
check("break_even(both negative) dominates", sw.break_even(-3, -1, 7.0)[1], "dominates")
check("break_even(both negative) r* undefined", sw.break_even(-3, -1, 7.0)[0], None)
check("break_even(both positive) dominated", sw.break_even(3, 1, 7.0)[1], "dominated")
check("break_even(over only, 0 under) dominates", sw.break_even(-5, 0, 7.0)[1], "dominates")
check("break_even(nothing changed)", sw.break_even(0, 0, 7.0)[1], "identical")

# ---------------------------------------------------------------- CLI fixture
tmp = tempfile.mkdtemp(prefix="scoreweighted-")
d_fx9, d_ds0, d_tol = (os.path.join(tmp, n) for n in ("fx9", "ds0", "fxtol"))
for d in (d_fx9, d_ds0, d_tol):
    os.makedirs(d)

gt_fx9 = write_gt(os.path.join(d_fx9, "sim.v0.gt.tsv"), FX9_GT)
gt_ds0 = write_gt(os.path.join(d_ds0, "sim.v0.gt.tsv"), DS0_GT)
gt_tol = write_gt(os.path.join(d_tol, "sim.v0.gt.tsv"), FXTOL_GT)

px_fx9 = write_pxmt(os.path.join(d_fx9, "a.pxmt"), FX9_DEC,
                    runs_by_pid={1: ((0, 12), (40, 9))})   # run payload must be skipped
px_base = write_pxmt(os.path.join(d_ds0, "base.pxmt"), DS0_BASE)
px_var = write_pxmt(os.path.join(d_ds0, "var.pxmt"), DS0_VAR)
px_tol = write_pxmt(os.path.join(d_tol, "t.pxmt"), FXTOL_DEC)

out1 = os.path.join(tmp, "score.tsv")
repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
cmd = [sys.executable, os.path.join(repo, "bench", "scripts", "score_weighted.py"),
       "--label", "one", "--pxmt", px_fx9, "--gt", gt_fx9,
       "--label", "base", "--pxmt", px_base, "--gt", gt_ds0,
       "--label", "var", "--pxmt", px_var, "--gt", gt_ds0,
       "--baseline", "base", "--variant", "var",
       "--out", out1]
res = subprocess.run(cmd, capture_output=True, text=True, cwd=repo)
check("CLI exit 0", res.returncode, 0)
if res.returncode != 0:
    print(res.stdout, res.stderr)

with open(out1, newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))


def pick(scope, dataset, label, r=None, row_type="score"):
    for row in rows:
        if (row["row_type"] == row_type and row["scope"] == scope
                and row["dataset"] == dataset and row["label"] == label
                and (r is None or float(row["r"]) == r)):
            return row
    raise AssertionError(f"no row {scope}/{dataset}/{label}/{r}/{row_type}")


fx = pick("dataset", "fx9", "one")
check("fx9 pairs", int(fx["pairs"]), 9)
check("fx9 pairs_scored (unaligned pair removed)", int(fx["pairs_scored"]), 8)
check("fx9 unaligned pair not scored", int(fx["unaligned_pairs"]), 1)
check("fx9 mates_scored", int(fx["mates_scored"]), 18)
check("fx9 exact", int(fx["exact"]), 2)
check("fx9 overtrim (raw, incl. the mixed mate)", int(fx["overtrim"]), 3)
check("fx9 undertrim (raw, incl. the mixed mate)", int(fx["undertrim"]), 3)
check("fx9 mixed_trim pairs", int(fx["mixed_trim"]), 1)
check("fx9 drop_correct (full flag + floor)", int(fx["drop_correct"]), 4)
check("fx9 drop_wrong", int(fx["drop_wrong"]), 2)
check("fx9 keep_wrong", int(fx["keep_wrong"]), 2)
check("fx9 unaligned", int(fx["unaligned"]), 2)
check("fx9 over_units (weight-1 side)", int(fx["over_units"]), 5)
check("fx9 under_units (r side)", int(fx["under_units"]), 5)
check("fx9 min_len_floor", int(fx["min_len_floor"]), FLOOR)
check("fx9 read_len inferred from PXMT", int(fx["read_len"]), L)
check("fx9 weights_version stamped", fx["weights_version"], "v1")
check("fx9 all rows stamp the version", {r["weights_version"] for r in rows}, {"v1"})

# numerator(r) = exact + drop_correct - 1*(overtrim+drop_wrong) - r*(undertrim+keep_wrong)
#             = 2 + 4 - 1*5 - r*5 = 1 - 5r        (mixed mates already inside the raw counts)
WANT = {1.0: -4.0, 3.0: -14.0, 5.0: -24.0, 7.0: -34.0, 10.0: -49.0}
for r, num in WANT.items():
    row = pick("dataset", "fx9", "one", r)
    check(f"fx9 numerator at r={r:g}", float(row["weighted_numerator"]), num)
    check(f"fx9 score at r={r:g}", float(row["score"]), num / 8.0)
check("fx9 score at the anchor r=7", float(pick("dataset", "fx9", "one", 7.0)["score"]), -4.25)
check("fx9 score at r=7 exact float", float(pick("dataset", "fx9", "one", 7.0)["score"]), -34 / 8)
check("fx9 ladder is the pre-registered one",
      sorted(float(r["r"]) for r in rows if r["row_type"] == "score"
             and r["dataset"] == "fx9"), [1.0, 3.0, 5.0, 7.0, 10.0])
check("fx9 secondary bases destroyed > 0", int(fx["SECONDARY_bases_over"]) > 0, True)
check("fx9 secondary artefact bases left > 0", int(fx["SECONDARY_bases_under"]) > 0, True)

# ---- the false-trim ceiling is measured on CLEAN reads only: GT window = the
# whole read. fx9 has pids 0,1,5,6,8 clean (10 mates); only pid1 is cut (5 bp
# on each mate), so 2/10. A dropped clean read is a drop error, not a cut. pid8's
# impossible 200 bp cut is unaligned and stays out of the ceiling.
check("fx9 clean mates", int(fx["clean_mates"]), 10)
check("fx9 clean-read false trims (pid1's two mates)", int(fx["clean_false_trims"]), 2)
check("fx9 clean-read false-trim rate", float(fx["clean_false_trim_rate"]), 0.2)

# unaligned must not be folded in: pid8's two unaligned mates are inside
# mates_scored (18) but in no weighted term, and its pair leaves the denominator
# (pairs 9 -> pairs_scored 8) instead of being scored as exact or as an error
check("fx9 unaligned never in the numerator",
      float(pick("dataset", "fx9", "one", 7.0)["weighted_numerator"]), -34.0)
check("fx9 18 mates scored but 8 pairs in the denominator",
      (int(fx["mates_scored"]), int(fx["pairs_scored"])), (18, 8))

# ---- the aggregate pools counts, and stamps the version too
agg = pick("AGGREGATE", "*", "one", 7.0)
check("AGGREGATE fx9 score == dataset score (single corpus)",
      float(agg["score"]), -4.25)
check("AGGREGATE fx9 pairs", int(agg["pairs"]), 9)

# ---- break-even: base vs var on ds0
# base: exact 2, undertrim 2, overtrim 2, drop_correct 2
#       over_units 2, under_units 2, positive 4
# var:  overtrim 4, drop_wrong 2, drop_correct 2
#       over_units 6, under_units 0, positive 2
# delta_over = +4, delta_under = -2  =>  r* = 4/2 = 2.0, variant wins for r > 2
be = [r for r in rows if r["row_type"] == "break_even"
      and r["dataset"] == "ds0" and r["variant"] == "var"][0]
check("be baseline label", be["baseline"], "base")
check("be delta_over", int(be["delta_over"]), 4)
check("be delta_under", int(be["delta_under"]), -2)
check("be r*", float(be["r_star"]), 2.0)
check("be verdict above r*", be["break_even_status"], "above")
check("be variant better at the anchor", int(be["variant_better_at_anchor"]), 1)
# base score at 7: (2 + 2 - 1*2 - 7*2)/4 = -12/4 = -3.0
# var  score at 7: (0 + 2 - 1*6 - 7*0)/4 =  -4/4 = -1.0 ; delta = +2.0
check("be base score at r=7", float(pick("dataset", "ds0", "base", 7.0)["score"]), -3.0)
check("be var score at r=7", float(pick("dataset", "ds0", "var", 7.0)["score"]), -1.0)
check("be delta_score_at_anchor", float(be["delta_score_at_anchor"]), 2.0)
check("be aggregate break-even present",
      [r["dataset"] for r in rows if r["row_type"] == "break_even"], ["ds0", "*"])

# ---- --tol plumbing: the same fixture scores differently at 0 and 2
outs = {}
for tol in (0, 2):
    p = os.path.join(tmp, f"tol{tol}.tsv")
    r2 = subprocess.run([sys.executable,
                         os.path.join(repo, "bench", "scripts", "score_weighted.py"),
                         "--label", "t", "--pxmt", px_tol, "--gt", gt_tol,
                         "--tol", str(tol), "--out", p],
                        capture_output=True, text=True, cwd=repo)
    check(f"tol {tol} exit 0", r2.returncode, 0)
    with open(p, newline="") as fh:
        outs[tol] = list(csv.DictReader(fh, delimiter="\t"))
def pick_row(rows_, r):
    return [x for x in rows_ if float(x["r"]) == r][0]


check("tol 0: off-by-2 is overtrim", outs[0][0]["overtrim"], "1")
check("tol 0: exact", outs[0][0]["exact"], "1")
check("tol 2: off-by-2 is exact", outs[2][0]["exact"], "2")
check("tol 2: nothing overtrimmed", outs[2][0]["overtrim"], "0")
# one pair: tol0 numerator at r=7 = 1 - 1*1 = 0 ; tol2 = 2 - 0 = 2
check("tol 0 score at r=7", float(pick_row(outs[0], 7.0)["score"]), 0.0)
check("tol 2 score at r=7", float(pick_row(outs[2], 7.0)["score"]), 2.0)

# ---- the --cuts path: the same decisions as a TSV, `NA` for the unaligned pair
cuts_path = os.path.join(tmp, "fx9.cuts.tsv")
with open(cuts_path, "w", newline="") as fh:
    wtr = csv.writer(fh, delimiter="\t", lineterminator="\n")
    wtr.writerow(["read_id", "mate", "status", "t5", "t3", "read_len"])
    for pid, mate, status, t5, t3 in FX9_DEC:
        rid = header(pid, mate).decode().split()[0]
        if pid == 8:
            wtr.writerow([rid, mate, status, "NA", "NA", L])
        else:
            wtr.writerow([rid, mate, status, t5, t3, L])
p_cuts = os.path.join(tmp, "cuts.tsv")
r3 = subprocess.run([sys.executable,
                     os.path.join(repo, "bench", "scripts", "score_weighted.py"),
                     "--label", "onecuts", "--cuts", cuts_path, "--gt", gt_fx9,
                     "--out", p_cuts],
                    capture_output=True, text=True, cwd=repo)
check("cuts mode exit 0", r3.returncode, 0)
if r3.returncode != 0:
    print(r3.stdout, r3.stderr)
with open(p_cuts, newline="") as fh:
    cut_rows = list(csv.DictReader(fh, delimiter="\t"))
ck = cut_rows[0]
check("cuts mode reproduces every class count",
      [ck[k] for k in sw.CLASSES + ("mixed_trim",)],
      [fx[k] for k in sw.CLASSES + ("mixed_trim",)])
check("cuts mode NA is unaligned, not folded in",
      (int(ck["unaligned"]), int(ck["pairs_scored"])), (2, 8))
check("cuts mode score at r=7",
      float([x for x in cut_rows if float(x["r"]) == 7.0][0]["score"]), -4.25)

# ---- the weights file is the source of truth, not a constant in the script
w = sw.load_weights(os.path.join(repo, "bench", "weights", "trim-weights-v1.json"))
check("weights version", w["version"], "v1")
check("weights anchor", w["anchor_r"], 7.0)
check("weights ladder", w["ladder_r"], [1.0, 3.0, 5.0, 7.0, 10.0])
check("weights w_over", w["w_over"], 1.0)
check("default --weights resolves from the repo",
      os.path.basename(sw.resolve_weights(sw.WEIGHTS_DEFAULT)),
      "trim-weights-v1.json")

print()
if fails:
    print(f"WEIGHTED SCORE SELF-TEST FAILED: {len(fails)}")
    for f in fails:
        print(f"  - {f}")
    sys.exit(1)
print("WEIGHTED SCORE SELF-TEST PASSED")
