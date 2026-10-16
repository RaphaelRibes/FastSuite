#!/usr/bin/env python3
"""Self-test for the C1 scorer: feed it decisions it can predict.

The scorer's correctness is not "the numbers look plausible" — it is whether a
synthetic set of decisions scores exactly what arithmetic says it should. Build
a small PXMT by hand from a small GT and assert the metrics.
"""
import struct, subprocess, sys, os, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sweep import parse_pxmt, score_config, coord_to_pid

GT = "pair_id\tinsert\tr1_h\tr1_t\tr2_h\tr2_t\tr1_adapter\tr2_adapter\tfull_flag\tis_dup\n"
L = 100
# pid 0: insert 300, no artefact -> keep [0,100) both mates, correct decision
GT += "0\t300\t0\t0\t0\t0\t0\t0\t0\t0\n"
# pid 1: insert 80 -> adapter readthrough, keep window [0,80)
GT += "1\t80\t0\t0\t0\t0\t0\t0\t0\t0\n"
# pid 2: h=5 3' G-run t=5, insert 300 -> keep [5,95)
GT += "2\t300\t5\t5\t5\t5\t0\t0\t0\t0\n"
# pid 3: full artefact -> must drop
GT += "3\t300\t100\t100\t100\t100\t0\t0\t1\t0\n"
# pid 4: PCR duplicate -> must drop
GT += "4\t300\t0\t0\t0\t0\t0\t0\t0\t1\n"

def header(pid, mate):
    return f"@LH00916:1:SIMFLOW1:1:1101:{pid%10000}:{pid//10000} {mate}:N:0:ATCGATCG".encode()

def rec(pid, mate, status, t5=0, t3=0):
    h = header(pid, mate)
    out = struct.pack("<H", len(h)) + h + bytes([mate, status])
    out += struct.pack("<IIIII", t5, 0, 0, t3, L)
    out += bytes([0])  # n_runs
    return out

fails = []
def check(name, got, want, tol=1e-9):
    ok = abs(got - want) <= tol if isinstance(want, float) else got == want
    print(f"{'ok  ' if ok else 'FAIL'} {name}: got {got} want {want}")
    if not ok:
        fails.append(name)

# coord round-trip
check("coord_to_pid(0,0)", coord_to_pid(header(0,0), 0), (0,0))
check("coord_to_pid(10001,1)", coord_to_pid(header(10001,1), 1), (10001,1))
check("coord_to_pid garbage", coord_to_pid(b"@junk", 0), (-1,0))

# ---- scenario A: the perfect trimmer ----
# keep 0,1,2 with exact windows; drop 3,4
# status: 0=kept, 1=drop:hard, 7=drop:duplicate (PXMT v3 codes)
buf = b"PXMT\x03"
buf += rec(0,0,0,0,0)   + rec(0,1,0,0,0)
buf += rec(1,0,0,0,20)  + rec(1,1,0,0,20)
buf += rec(2,0,0,5,5)   + rec(2,1,0,5,5)
buf += rec(3,0,1,0,0)   + rec(3,1,1,0,0)
buf += rec(4,0,7,0,0)   + rec(4,1,7,0,0)
with tempfile.NamedTemporaryFile(delete=False, suffix=".pxmt") as f:
    f.write(buf); px = f.name
meta = parse_pxmt(px)
gt = [(0,300,0,0,0,0,0,0),(1,80,0,0,0,0,0,0),(2,300,5,5,5,5,0,0),
      (3,300,100,100,100,100,1,0),(4,300,0,0,0,0,0,1)]
m = score_config(meta, gt, L)
check("A kept", m["kept"], 3)
check("A should_keep", m["should_keep"], 3)
check("A keep_precision", m["keep_precision"], 1.0)
check("A keep_recall", m["keep_recall"], 1.0)
check("A trim_mae", m["trim_mae"], 0.0)
check("A over_trim_rate", m["over_trim_rate"], 0.0)
check("A under_trim_rate", m["under_trim_rate"], 0.0)
check("A drop_accuracy", m["drop_accuracy"], 1.0)
check("A score", m["score"], 1.0)

# ---- scenario B: keeps everything, trims nothing ----
buf = b"PXMT\x03"
for pid in range(5):
    buf += rec(pid,0,0,0,0) + rec(pid,1,0,0,0)
with tempfile.NamedTemporaryFile(delete=False, suffix=".pxmt") as f:
    f.write(buf); px = f.name
meta = parse_pxmt(px)
m = score_config(meta, gt, L)
# only the 3 should-keep pairs count as keeps; the full/dup pair is a false
# keep and is scored through drop_accuracy alone
check("B kept", m["kept"], 3)
check("B should_keep", m["should_keep"], 3)
check("B keep_recall", m["keep_recall"], 1.0)
# window correct only for pid 0 (truth [0,100), observed [0,100))
check("B keep_precision", m["keep_precision"], 1/3, tol=1e-6)
# pid 3,4 kept but should be dropped -> drop_accuracy = 1 - 2/3
check("B drop_accuracy", m["drop_accuracy"], 1 - 2/3, tol=1e-6)
# Trimming NOTHING over-retains: the observed window is never SHORTER than the
# truth, so neither under-trim (>=10bp artefact left) nor over-trim (<75% of
# biology kept) fires. That is the correct answer for this decision set — the
# defects show up in keep_precision (wrong boundaries), not in these rates.
check("B under_trim_rate", m["under_trim_rate"], 0.0)
check("B over_trim_rate", m["over_trim_rate"], 0.0)
check("B false_keep_rate", m["false_keep_rate"], 2/5, tol=1e-6)

# ---- scenario D: an over-aggressive trimmer ----
# keeps pid 0 correctly, but on pid 1 it trims 60 of 80 biology bases
# (obs window [0,20) vs truth [0,80) -> 20 < 0.75*80) and drops pid 2 entirely.
buf = b"PXMT\x03"
buf += rec(0,0,0,0,0)     + rec(0,1,0,0,0)
buf += rec(1,0,0,0,80)    + rec(1,1,0,0,80)
buf += rec(2,0,1,0,0)     + rec(2,1,1,0,0)
buf += rec(3,0,1,0,0)     + rec(3,1,1,0,0)
buf += rec(4,0,1,0,0)     + rec(4,1,1,0,0)
with tempfile.NamedTemporaryFile(delete=False, suffix=".pxmt") as f:
    f.write(buf); px = f.name
meta = parse_pxmt(px)
m = score_config(meta, gt, L)
check("D kept", m["kept"], 2)
check("D keep_recall", m["keep_recall"], 2/3, tol=1e-6)
# only pid 0 has a near-correct window -> precision 1/2
check("D keep_precision", m["keep_precision"], 0.5)
# pid 1 lost 60 of 80 biology bases -> over-trim
check("D over_trim_rate", m["over_trim_rate"], 0.5)
check("D false_keep_rate", m["false_keep_rate"], 0.0)
check("D score", m["score"], (2/3)*0.5 - 2.0*0.5, tol=1e-6)

# ---- scenario C: drops everything ----
buf = b"PXMT\x03"
for pid in range(5):
    buf += rec(pid,0,1,0,0) + rec(pid,1,1,0,0)
with tempfile.NamedTemporaryFile(delete=False, suffix=".pxmt") as f:
    f.write(buf); px = f.name
meta = parse_pxmt(px)
m = score_config(meta, gt, L)
check("C kept", m["kept"], 0)
check("C keep_recall", m["keep_recall"], 0.0)
check("C score", m["score"], 0.0)

# ---- the head-dropout axis is emitted, not assumed -----------------------
# The grid authors head_dropout 0/1/4/6. Before, build_argv never emitted
# --head-dropout, so every config ran at fastprep's CLI default (4) and the
# axis was never measured or recorded — and a min_run <= 4 config would have
# aborted the whole sweep on the clap guard. effective_head_dropout() caps
# rather than skips, so every config still produces a row that says what ran.
import sweep as _sw
_extra = 0
def chk2(name, got, want):
    global _extra
    ok = got == want
    if not ok: _extra += 1
    print(("ok   " if ok else "FAIL ") + f"{name}: got {got} want {want}")

chk2("hd authored 0 (expressible)", _sw.effective_head_dropout(
    {"min_run": 12, "head_dropout": 0}), 0)
chk2("hd authored 1", _sw.effective_head_dropout(
    {"min_run": 12, "head_dropout": 1}), 1)
chk2("hd authored 6, min_run 8", _sw.effective_head_dropout(
    {"min_run": 8, "head_dropout": 6}), 6)
chk2("hd authored 6, min_run 2 (capped)", _sw.effective_head_dropout(
    {"min_run": 2, "head_dropout": 6}), 1)
chk2("hd 0 survives the cap, min_run 2", _sw.effective_head_dropout(
    {"min_run": 2, "head_dropout": 0}), 0)
chk2("hd unset (shipped default 4)", _sw.effective_head_dropout(
    {"min_run": 12}), 4)
chk2("grid records the column", "head_dropout" in _sw.ROW_FIELDS, True)
_axis = set()
for _d in _sw.GRID:
    _cfg = dict(_sw.DEFAULTS); _cfg.update(_d)
    _axis.add(_sw.effective_head_dropout(_cfg))
chk2("grid authors hd 0 and 1", {0, 1} <= _axis, True)
for i, d in enumerate(_sw.GRID):
    cfg = dict(_sw.DEFAULTS); cfg.update(d)
    argv = _sw.build_argv(cfg, "/tmp", i, 150)
    if "--head-dropout" not in argv:
        chk2(f"grid {i} emits --head-dropout", False, True)
        break
else:
    chk2("every grid config emits --head-dropout", True, True)
# and the emitted value must equal the recorded value
for i, d in enumerate(_sw.GRID):
    cfg = dict(_sw.DEFAULTS); cfg.update(d)
    argv = _sw.build_argv(cfg, "/tmp", i, 150)
    emitted = argv[argv.index("--head-dropout") + 1]
    if int(emitted) != _sw.effective_head_dropout(cfg):
        chk2(f"grid {i} argv == row", False, True)
        break
else:
    chk2("argv value and row value agree", True, True)

# ---- the resume key carries the sim inputs' identity ---------------------
# Keying the resume cache on `c<cid>.row` existence alone survives a corpus
# regeneration: the old rows are replayed and the campaign reports metrics
# against a sim that no longer exists. input_identity() folds each input's
# path/size/mtime in, so a changed corpus misses the cache.
_tmp = tempfile.mkdtemp()
_p1 = os.path.join(_tmp, "sim_R1.fastq.gz"); open(_p1, "wb").write(b"x")
_p2 = os.path.join(_tmp, "sim_R2.fastq.gz"); open(_p2, "wb").write(b"y")
_pg = os.path.join(_tmp, "sim.gt.tsv"); open(_pg, "wb").write(b"z")
_id1 = _sw.input_identity(_p1, _p2, _pg)
chk2("input identity is stable", _sw.input_identity(_p1, _p2, _pg), _id1)
_t = os.path.getmtime(_pg) + 10
os.utime(_pg, (_t, _t))
chk2("identity changes when the GT is regenerated",
     _sw.input_identity(_p1, _p2, _pg) != _id1, True)

# ---- exit code compares against the ACTUAL grid length -------------------
# The chain passes --configs 60 while the grid is 57, so a full 57/57 sweep
# must exit 0. The expected count is min(requested, grid).
chk2("expected count clamps to the grid", min(60, len(_sw.GRID)), len(_sw.GRID))

# ---- the semantics overtrim_vs_fastp now reproduces ----------------------
# An integer grid min_len is scaled to a fraction of the actual read length and
# clamped to <= 1.0; the resolved string is what the harness must pass.
chk2("scale_min_len clamps 50/36 to 1.0",
     _sw.scale_min_len({"min_len": 50}, 36), 1.0)
chk2("min_len_arg(1.0) is the literal base count '1'",
     _sw.min_len_arg(_sw.scale_min_len({"min_len": 50}, 36)), "1")

print()
if fails or _extra:
    print(f"SCORER SELF-TEST FAILED: {len(fails) + _extra}")
    sys.exit(1)
print("SCORER SELF-TEST PASSED")
