#!/usr/bin/env python3
"""Headroom census for PE-overlap consensus correction, with per-base truth.

Why this exists
---------------
The trim-accuracy harness scores CUT POSITIONS against fastsim ground truth.
The quality axis this probe addresses is different: BASE ERRORS left in the
output. fastp's `--correction` is its only error-reducing feature (a rigid
rule on mismatches inside the PE overlap); fastprep has NO base correction at
all today. Before implementing a quality-weighted posterior consensus, measure
the ceiling and its price:

  * how much of the injected error mass sits in the mate-overlap region
    (where two independent observations of the same fragment base exist),
  * how much of it fastp --correction actually fixes against real truth, and
    how often it CORRUPTS a correct base (the second axis -- never one side),
  * how much a posterior rule fixes/misfixes on the same disagreements, as a
    function of its log-odds gate (threshold sweep, not one operating point).

Truth is exact: the corpus is BUILT here (the undertrim_ab.py pattern), so
every injected substitution is known per (pair, mate, position), and the
error model matches fastsim's (`p = 10^(-q/10)`, uniform over the other 3
bases), so posterior numbers are self-consistent with the simulator.

Coordinates
-----------
Mate 2 carries the reverse complement of the fragment. ALL comparisons in
this probe happen in FRAGMENT orientation: the mate-2 base at read position
`ins-1-fpos` is complemented before it meets mate 1's base or the truth
base. The first version of this file compared read-orientation bases and
reported ~100% of shared positions as disagreements; verify_truth() now
asserts mate agreement on the pre-error reads so that class of bug cannot
return silently.

Scoring rules baked in (references/measurement-pitfalls.md):
  * BOTH axes always: residual error bases AND correct bases corrupted.
  * Rates print beside the counts they came from and beside retention
    (reads out per arm) -- never at unequal retention silently.
  * Mate agreement is not truth: two mates agreeing on the same wrong base
    are counted from truth, not from agreement (agree_both_wrong).
  * Corrections never change retention (they rewrite bases in place), so the
    correction arms are directly comparable to their own no-correction arm.
  * The output TSV is overwritten each run; the run asserts its own census
    counts and an identity (no-op) arm before printing anything (a probe
    that cannot fail is not a probe).

Usage
-----
    python3 bench/scripts/consensus_headroom.py [--pairs 20000] [--out TSV]
        [--fastp /path/to/fastp] [--fastprep /path/to/fastprep]
"""
import argparse
import math
import os
import random
import subprocess
import tempfile

COMP = str.maketrans("ACGTN", "TGCAN")
CBASE = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}   # single-base complement
TRUSEQ_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"   # TruSeq R1 readthrough prefix
TRUSEQ_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"   # TruSeq R2 readthrough prefix

L = 150

# NovaSeq-ish binned quality with a 3' decay over the last 40 bases.
Q_BINS = [(2, 0.004), (12, 0.016), (23, 0.10), (37, 0.88)]

# Posterior gates to sweep: |log10 odds(true=b1 vs true=b2)| >= thr to act.
POSTERIOR_THRS = (0.5, 1.0, 2.0, 3.0, 5.0, 10.0)


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]

def dna(n: int, rng: random.Random) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))

def tile(s: str, n: int) -> str:
    return (s * ((n // len(s)) + 1))[:n]

def draw_q(pos: int, rng: random.Random) -> int:
    bins = Q_BINS
    if pos > L - 40:
        shift = (pos - (L - 40)) / 40.0
        bins = [(2, 0.004 + 0.02 * shift), (12, 0.016 + 0.10 * shift),
                (23, 0.10 + 0.15 * shift), (37, 0.88 - 0.27 * shift)]
    r = rng.random()
    acc = 0.0
    for q, w in bins:
        acc += w
        if r <= acc:
            return q
    return bins[-1][0]

def q_to_p(q: int) -> float:
    return 10.0 ** (-q / 10.0)

def build_corpus(n_pairs: int, seed: int, tmpdir: str, err_scale: float = 1.0):
    """PE corpus with exact per-base truth and a realistic overlap mix.

    Insert regimes (read length 150):
      short   I in [80, 150)   full overlap + adapter readthrough   25%
      medium  I in [150, 299)  partial overlap                      45%
      long    I in [300, 420]  no overlap (consensus cannot fire)   30%
    Every emitted read is exactly L bases: the readthrough tail tiles the
    TruSeq prefix as needed. No G-runs, no Ns: this probe isolates BASE
    ERRORS (the trim axes have their own harnesses).

    err_scale multiplies the injection probability while the QUALITY STRING
    stays unchanged -- a deliberate Q-miscalibration model (the claim `p_err =
    10^(-q/10)` becomes false). Both the fastp rule and the posterior rule
    still read the same Q strings, so the sweep measures how each degrades
    when the Q values lie: the oracle-trap stress test.
    """
    rng = random.Random(seed)
    r1_path = os.path.join(tmpdir, "in_R1.fastq")
    r2_path = os.path.join(tmpdir, "in_R2.fastq")
    truth = {}
    with open(r1_path, "w") as f1, open(r2_path, "w") as f2:
        for pid in range(n_pairs):
            u = rng.random()
            if u < 0.25:
                ins = rng.randrange(80, 150)
            elif u < 0.70:
                ins = rng.randrange(150, 299)
            else:
                ins = rng.randrange(300, 421)
            frag = dna(ins, rng)
            if ins >= L:
                pre1 = frag[:L]
                pre2 = revcomp(frag[ins - L:ins])
            else:
                pre1 = frag + tile(TRUSEQ_R1, L - ins)
                pre2 = revcomp(frag) + tile(TRUSEQ_R2, L - ins)
            assert len(pre1) == L and len(pre2) == L, "read construction broke"
            mates = []
            err = []
            for s in (pre1, pre2):
                q = [draw_q(i, rng) for i in range(L)]
                e = []
                b = list(s)
                for i in range(L):
                    if rng.random() < min(1.0, err_scale * q_to_p(q[i])):
                        true_base = b[i]
                        b[i] = rng.choice([x for x in "ACGT" if x != true_base])
                        e.append((i, true_base, b[i]))
                mates.append(("".join(b), q))
                err.append(e)
            truth[pid] = dict(frag=frag, ins=ins, pre=[pre1, pre2],
                              seq=[mates[0][0], mates[1][0]],
                              qual=[mates[0][1], mates[1][1]], err=err)
            for mate in (0, 1):
                fh = f1 if mate == 0 else f2
                fh.write(f"@pair{pid}/{mate + 1}\n{mates[mate][0]}\n+\n"
                         f"{''.join(chr(x + 33) for x in mates[mate][1])}\n")
    return r1_path, r2_path, truth


def verify_truth(truth):
    """Self-check invariants -- run before anything is measured.

    (a) On the PRE-ERROR reads the two mates agree at 100% of shared
        fragment positions (after complementing mate 2). This is the check
        that catches orientation/offset bugs in every consumer below.
    (b) The emitted reads differ from the pre-error reads at exactly the
        recorded error positions -- error records and sequences cannot drift.
    """
    for pid, t in truth.items():
        ins = t["ins"]
        lo, hi = max(0, ins - L), min(ins, L)
        for fpos in range(lo, hi):
            j1, j2 = fpos, ins - 1 - fpos
            assert t["pre"][0][j1] == t["frag"][fpos], (pid, fpos, "R1 pre mismatch")
            assert CBASE[t["pre"][1][j2]] == t["frag"][fpos], \
                (pid, fpos, "R2 pre not reverse-complement of truth")
        for mate in (0, 1):
            observed = {i for i in range(L) if t["seq"][mate][i] != t["pre"][mate][i]}
            recorded = {p for p, _, _ in t["err"][mate]}
            assert observed == recorded, (pid, mate, "error records vs sequences drift")


def read_fastq(path):
    """Strict 4-line FASTQ parse. Blank lines are tolerated ONLY between
    records; a malformed record raises with its line number. (A parser that
    silently skips odd lines desyncs and reports plausible-looking rates over
    garbage -- the failure mode this replaces produced an IndexError on a
    one-character line read as a header.) Mate identity is recovered by the
    caller from which file the record came from: tools rewrite `/1` suffixes."""
    with open(path, "rt") as f:
        lines = f.read().split("\n")
    if lines and lines[-1] == "":
        lines.pop()          # trailing newline artifact only; not every blank
    out = {}
    i = 0
    while i < len(lines):
        if lines[i].strip() == "":
            i += 1
            continue
        if not lines[i].startswith("@"):
            raise ValueError(f"{path}:{i + 1}: expected '@' header, got {lines[i]!r}")
        if i + 3 >= len(lines):
            raise ValueError(f"{path}:{i + 1}: truncated record (no qualities)")
        name = lines[i][1:].split()
        if not name:
            raise ValueError(f"{path}:{i + 1}: empty read name in {lines[i]!r}")
        if len(lines[i + 1]) != len(lines[i + 3]):
            raise ValueError(f"{path}:{i + 1}: seq/qual length mismatch")
        out[name[0]] = (lines[i + 1], lines[i + 3])
        i += 4
    return out


def best_offset(out_seq, ref):
    """Locate an output read on the pre-error read: end-only trimming makes
    the output a window of `ref` up to base corrections. Wrong offsets die
    after a few characters (~75% mismatch), so the scan stays cheap."""
    n = len(out_seq)
    limit = max(3, n // 12)          # corrections are rare; anything above is unmappable
    best = (limit + 1, -1)
    for off in range(0, len(ref) - n + 1):
        mm = 0
        for i in range(n):
            if out_seq[i] != ref[off + i]:
                mm += 1
                if mm > limit:
                    break
        if mm < best[0]:
            best = (mm, off)
            if mm == 0:
                break
    return best


def rule_census(truth, min_ovl=30):
    """Truth-known price of correction rules on the overlap region, in BASE
    units over both mates, all in FRAGMENT orientation.

    before_err      : error bases in the shared region, no correction
                      (must equal errors_in_shared -- asserted in main)
    fastp rule      : correct only when one mate is >=Q30 and the other <=Q14;
                      the high-quality base wins
    posterior rule  : log10 posterior odds between the two observed bases
                      ((1-p1)*p2/3 vs (1-p2)*p1/3; the /3 cancels), acting
                      only past a gate |odds| >= thr -- swept, not fixed
    fixed           : error bases the rule REMOVES
    corrupted       : correct input bases the rule DESTROYS (the second axis)
    after_err       : error bases left after the rule
    Agreements where both mates carry the same wrong base are left alone by
    every rule (no evidence to act on) and counted from truth.
    """
    c: dict = dict(pairs_overlap=0, pairs_no_overlap=0,
             shared_positions=0, shared_bases=0,
             injected_errors=0, errors_in_shared=0,
             agree=0, agree_both_wrong=0, disagree=0,
             dis_one_sided=0, dis_both_wrong=0, before_err=0,
             fastp=dict(fired=0, fixed=0, corrupted=0, after_err=0),
             post={thr: dict(fired=0, fixed=0, corrupted=0, after_err=0)
                   for thr in POSTERIOR_THRS})

    # every rule starts from the no-correction baseline: each position's
    # `before` errors are pre-added to every rule's after_err, and apply()
    # records only the DELTA at positions it fires on. Agreement positions
    # (including agree_both_wrong) therefore stay in every ledger -- a rule
    # that cannot act there must still carry those errors.
    all_rules = [c["fastp"]] + [c["post"][thr] for thr in POSTERIOR_THRS]

    def apply(s: dict, before: int, pick: str, true_b: str):
        after = 2 * (pick != true_b)     # both mates end up carrying `pick`
        s["fired"] += 1
        s["after_err"] += after - before
        if after < before:
            s["fixed"] += before - after
        elif after > before:
            s["corrupted"] += after - before

    for pid, t in truth.items():
        ins = t["ins"]
        lo, hi = max(0, ins - L), min(ins, L)
        n_inj = sum(len(e) for e in t["err"])
        c["injected_errors"] += n_inj
        if hi - lo < min_ovl:
            c["pairs_no_overlap"] += 1
            continue
        c["pairs_overlap"] += 1
        err_pos = [{p for p, _, _ in e} for e in t["err"]]
        for mate in (0, 1):
            for p, _, _ in t["err"][mate]:
                rp = p if mate == 0 else ins - 1 - p
                if lo <= rp < hi:
                    c["errors_in_shared"] += 1
        for fpos in range(lo, hi):
            j1, j2 = fpos, ins - 1 - fpos
            if not (0 <= j1 < L and 0 <= j2 < L):
                continue
            c["shared_positions"] += 1
            c["shared_bases"] += 2
            b1, q1 = t["seq"][0][j1], t["qual"][0][j1]
            b2 = CBASE[t["seq"][1][j2]]          # fragment orientation
            q2 = t["qual"][1][j2]
            true_b = t["frag"][fpos]
            before = (b1 != true_b) + (b2 != true_b)
            c["before_err"] += before
            for s in all_rules:              # no-correction baseline for this position
                s["after_err"] += before
            if b1 == b2:
                c["agree"] += 1
                if b1 != true_b:
                    c["agree_both_wrong"] += 1
                continue
            c["disagree"] += 1
            one_sided = (j1 in err_pos[0]) != (j2 in err_pos[1])
            c["dis_one_sided" if one_sided else "dis_both_wrong"] += 1
            # ---- fastp's documented rigid rule
            if q1 >= 30 and q2 <= 14:
                pick = b1
            elif q2 >= 30 and q1 <= 14:
                pick = b2
            else:
                pick = None
            if pick is not None:
                apply(c["fastp"], before, pick, true_b)
            # ---- posterior rule, swept over log-odds gates
            p1, p2 = q_to_p(q1), q_to_p(q2)
            lodds = math.log10((1 - p1) * p2) - math.log10((1 - p2) * p1)
            for thr in POSTERIOR_THRS:
                if lodds >= thr:
                    pick2 = b1
                elif lodds <= -thr:
                    pick2 = b2
                else:
                    continue
                apply(c["post"][thr], before, pick2, true_b)
    return c


def score_tool(truth, out1, out2):
    """Per-base error rate of the OUTPUT against the pre-error read (the exact
    molecule the instrument read), plus corruption: output wrong where the
    input base was already correct, plus the tool's own edits (output base
    differs from the observed INPUT base) and the Q it wrote at edited bases.

    Mate identity comes from which file a record was read from (tools rewrite
    `/1` suffixes), never from the name. Fragment positions (read position <
    insert length) are counted separately from adapter readthrough so the
    error-rate denominator cannot be inflated by leftover tail bases. Values
    are (seq, qual) pairs or bare seq strings; Q stats cover fragment bases
    only and stay empty when qual is absent."""
    s = dict(bases=0, frag_bases=0, errs=0, frag_errs=0, corrupted=0,
             mapped=0, unmapped=0, edits=0, frag_edits=0,
             q_out_sum=0, q_out_n=0, q30_out=0,
             q_edit_in_sum=0, q_edit_out_sum=0)
    for mate, out in ((0, out1), (1, out2)):
        for name, val in out.items():
            seq, qual = val if isinstance(val, tuple) else (val, None)
            tok = name.split()[0]
            pid = int(tok.replace("pair", "").split("/")[0])
            t = truth[pid]
            ref = t["pre"][mate]
            obs = t["seq"][mate]
            mm, off = best_offset(seq, ref)
            if off < 0:
                s["unmapped"] += 1
                continue
            s["mapped"] += 1
            s["bases"] += len(seq)
            err_pos = {p for p, _, _ in t["err"][mate]}
            for i, ch in enumerate(seq):
                j = off + i                      # position in the input read
                on_frag = j < t["ins"]           # adapter readthrough is not fragment
                if on_frag:
                    s["frag_bases"] += 1
                    if qual is not None:
                        qo = ord(qual[i]) - 33
                        s["q_out_sum"] += qo
                        s["q_out_n"] += 1
                        if qo >= 30:
                            s["q30_out"] += 1
                if ch != obs[j]:                 # the tool rewrote this base
                    s["edits"] += 1
                    if on_frag:
                        s["frag_edits"] += 1
                        if qual is not None:
                            s["q_edit_in_sum"] += t["qual"][mate][j]
                            s["q_edit_out_sum"] += ord(qual[i]) - 33
                if ch == ref[j]:
                    continue
                s["errs"] += 1
                if on_frag:
                    s["frag_errs"] += 1
                if j not in err_pos:
                    s["corrupted"] += 1          # the tool broke a correct base
    return s


def div(a, b):
    return a / b if b else 0.0


def posterior_correct(truth, thr, min_ovl=30):
    """Prototype correction arm: at overlap disagreements, pick the base with
    the higher two-hypothesis posterior and rewrite ONLY the losing mate's base
    and Q (winner untouched). Acts when |log10 odds| >= thr.

    Odds = P(true=b1 | obs) / P(true=b2 | obs) = (1-p1)*p2 / ((1-p2)*p1)
    (wrong-base-uniform /3 factors cancel). The rewritten base's Q is the
    posterior Q of the consensus, capped at 40: evidence-based rescore, not
    inflation. Pairs below min_ovl shared bases are left alone -- same filter
    the census applies, so this arm must reproduce its numbers exactly.

    Overlap geometry is truth-known -- this arm isolates the CORRECTION rule;
    overlap detection quality is a separate lever. Asserted against the census
    in main: edits == fired, corrupted == census corrupted, and the frag error
    drop == fixed - corrupted."""
    out1, out2 = {}, {}
    for pid, t in truth.items():
        ins = t["ins"]
        lo, hi = max(0, ins - L), min(ins, L)
        if hi - lo < min_ovl:
            # emit untouched so retention and scoring stay comparable
            out1[f"pair{pid}/1"] = (t["seq"][0],
                                    "".join(chr(x + 33) for x in t["qual"][0]))
            out2[f"pair{pid}/2"] = (t["seq"][1],
                                    "".join(chr(x + 33) for x in t["qual"][1]))
            continue
        s1, s2 = list(t["seq"][0]), list(t["seq"][1])
        q1, q2 = list(t["qual"][0]), list(t["qual"][1])
        for fpos in range(lo, hi):
            j1, j2 = fpos, ins - 1 - fpos
            if not (0 <= j1 < L and 0 <= j2 < L):
                continue
            b1, b2 = s1[j1], CBASE[s2[j2]]      # fragment orientation
            if b1 == b2:
                continue
            p1, p2 = q_to_p(q1[j1]), q_to_p(q2[j2])
            lodds = math.log10((1 - p1) * p2) - math.log10((1 - p2) * p1)
            if abs(lodds) < thr:
                continue
            qn = min(40, int(10 * math.log10(1 + 10 ** abs(lodds))))
            if lodds > 0:                       # mate 1 wins, mate 2 loses
                s2[j2] = CBASE[b1]
                q2[j2] = qn
            else:                               # mate 2 wins, mate 1 loses
                s1[j1] = b2
                q1[j1] = qn
        out1[f"pair{pid}/1"] = ("".join(s1), "".join(chr(x + 33) for x in q1))
        out2[f"pair{pid}/2"] = ("".join(s2), "".join(chr(x + 33) for x in q2))
    return out1, out2




def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261022)
    ap.add_argument("--err-scale", type=float, default=1.0,
                    help="multiply injection probability, Q strings unchanged "
                         "(Q-miscalibration stress; 1.0 = calibrated)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--fastp", default="/tmp/consensus-probe/fastp")
    ap.add_argument("--fastprep", default=os.path.expanduser(
        "~/hermes_scratch/FastSuite/tools/fastprep/target/release/fastprep"))
    args = ap.parse_args()

    out_tsv = args.out or os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "bench", "results", "consensus-headroom.tsv")
    tmp = tempfile.mkdtemp(prefix="consensus-headroom-")
    r1, r2, truth = build_corpus(args.pairs, args.seed, tmp, args.err_scale)
    assert truth, "corpus build produced no pairs"
    verify_truth(truth)

    cen = rule_census(truth)
    assert cen["injected_errors"] > 0, "no errors injected -- probe is vacuous"
    assert cen["pairs_overlap"] > 0, "no overlapping pairs -- geometry broken"
    assert cen["disagree"] > 0, "no mate disagreements -- errors missing"
    assert cen["shared_bases"] > 0, "no shared bases -- overlap mapping broken"
    # ledger invariants: the census must reconcile with the injection metadata
    assert cen["before_err"] == cen["errors_in_shared"], \
        (f"census counted {cen['before_err']} error bases but metadata says "
         f"{cen['errors_in_shared']} -- coordinate bug is back")
    assert cen["before_err"] / max(1, cen["shared_bases"]) < 0.02, \
        "implausible error rate on shared bases -- comparison is broken"
    assert cen["dis_one_sided"] + cen["dis_both_wrong"] == cen["disagree"]
    for label, s in [("fastp-rule", cen["fastp"])] + \
                    [(f"posterior@{t}", cen["post"][t]) for t in POSTERIOR_THRS]:
        assert s["after_err"] == cen["before_err"] - s["fixed"] + s["corrupted"], \
            (label, s, cen["before_err"])
        assert s["fixed"] == 0 or s["fired"] > 0

    # identity arm: raw input reads must score exactly as injected, no corruption
    ident = score_tool(truth,
                       {f"pair{p}/1": t["seq"][0] for p, t in truth.items()},
                       {f"pair{p}/2": t["seq"][1] for p, t in truth.items()})
    assert ident["errs"] == cen["injected_errors"], \
        (f"identity arm found {ident['errs']} errors, injected "
         f"{cen['injected_errors']} -- scorer disagrees with truth")
    assert ident["corrupted"] == 0, "identity arm corrupts without touching a base"
    assert ident["unmapped"] == 0 and ident["bases"] == 2 * args.pairs * L

    rows = []

    def run_fastp(tag, extra):
        o1, o2 = os.path.join(tmp, f"{tag}_R1.fastq"), os.path.join(tmp, f"{tag}_R2.fastq")
        cmd = [args.fastp, "-i", r1, "-I", r2, "-o", o1, "-O", o2, "-w", "2",
               "--json", os.path.join(tmp, f"{tag}.json"),
               "--html", os.path.join(tmp, f"{tag}.html")] + extra
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 0, f"fastp {tag} failed: {r.stderr[-400:]}"
        return read_fastq(o1), read_fastq(o2)

    def run_fastprep(tag, extra):
        o1, o2 = os.path.join(tmp, f"{tag}_R1.fastq"), os.path.join(tmp, f"{tag}_R2.fastq")
        cmd = [args.fastprep, "-1", r1, "-2", r2, "-o", o1, "-p", o2,
               "-j", "2", "--force", "-q", "--no-dedup"] + extra
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 0, f"fastprep {tag} failed: {r.stderr[-400:]}"
        return read_fastq(o1), read_fastq(o2)

    arms = [
        ("fastp-default", lambda: run_fastp("fpdef", [])),
        ("fastp-correct", lambda: run_fastp("fpcor", ["--correction"])),
        # filters relaxed so retention matches fastp's; trimming stays on
        ("fastprep-default", lambda: run_fastprep(
            "fpadef", ["--no-mean-q", "--min-len", "1", "--head-dropout", "0"])),
        # The shipped feature: identical run + the opt-in flag. Trim decisions
        # run first and correction only rewrites the doubly-retained overlap,
        # so retention and trim boundaries must match fastprep-default exactly
        # (asserted below on reads_out/bases_out).
        ("fastprep-correct@2.0", lambda: run_fastprep(
            "fpac20", ["--no-mean-q", "--min-len", "1", "--head-dropout", "0",
                       "--pe-correct=2.0"])),
        ("fastprep-correct@0.5", lambda: run_fastprep(
            "fpac05", ["--no-mean-q", "--min-len", "1", "--head-dropout", "0",
                       "--pe-correct=0.5"])),
        # prototype consensus: posterior gate + evidence-based Q rescore on the
        # losing base only. No trimming, so retention is trivially matched and
        # frag_err_per_kb is the column comparable against the tool arms.
        ("posterior@2.0", lambda: posterior_correct(truth, 2.0)),
        ("posterior@0.5", lambda: posterior_correct(truth, 0.5)),
    ]
    for tag, fn in arms:
        o1, o2 = fn()
        sc = score_tool(truth, o1, o2)
        n_out = len(o1) + len(o2)
        rows.append(dict(
            arm=tag, reads_in=2 * args.pairs, reads_out=n_out,
            dropped=2 * args.pairs - n_out,
            bases_out=sc["bases"], frag_bases=sc["frag_bases"],
            err_out=sc["errs"], frag_errs=sc["frag_errs"],
            frag_err_per_kb=round(1000.0 * sc["frag_errs"]
                                  / max(1, sc["frag_bases"]), 4),
            corrupted=sc["corrupted"], edits=sc["edits"],
            mean_q_frag=round(div(sc["q_out_sum"], sc["q_out_n"]), 3),
            q30_frac=round(div(sc["q30_out"], sc["q_out_n"]), 5),
            edit_q_in=round(div(sc["q_edit_in_sum"], sc["frag_edits"]), 2),
            edit_q_out=round(div(sc["q_edit_out_sum"], sc["frag_edits"]), 2),
            mapped=sc["mapped"], unmapped=sc["unmapped"],
        ))
        assert sc["mapped"] > 0, f"{tag}: no output read mapped back to truth"
        if tag.startswith("posterior@"):
            cr = cen["post"][float(tag.split("@")[1])]
            # the end-to-end arm must reproduce the rule census exactly:
            # one rewrite per fired position, same corruption count, and the
            # frag error drop must equal fixed - corrupted in base units.
            assert sc["edits"] == cr["fired"], (tag, sc["edits"], cr["fired"])
            assert sc["corrupted"] == cr["corrupted"], (tag, sc["corrupted"], cr)
            assert ident["frag_errs"] - sc["frag_errs"] \
                == cr["fixed"] - cr["corrupted"], (tag, cr)

    assert all(r["bases_out"] > 0 for r in rows), "an arm produced no bases"
    # Correction rewrites bases in place: it must never move a trim boundary.
    # Identical retention AND identical emitted base counts between the
    # default and corrected arms is the structural proof.
    by_arm = {r["arm"]: r for r in rows}
    for tgt in ("fastprep-correct@2.0", "fastprep-correct@0.5"):
        if tgt in by_arm:
            assert by_arm[tgt]["reads_out"] == by_arm["fastprep-default"]["reads_out"], \
                f"{tgt}: retention differs from fastprep-default"
            assert by_arm[tgt]["bases_out"] == by_arm["fastprep-default"]["bases_out"], \
                f"{tgt}: emitted base count differs (trim boundary moved)"

    tot = cen["injected_errors"]
    print(f"# pairs={args.pairs} injected_errors={tot} "
          f"overlap_pairs={cen['pairs_overlap']}/{args.pairs} "
          f"shared_positions={cen['shared_positions']} "
          f"shared_mate_bases={cen['shared_bases']}")
    print(f"# error mass inside >=30bp overlap: {cen['errors_in_shared']}/{tot} "
          f"= {100.0*cen['errors_in_shared']/tot:.1f}%  (consensus ceiling)")
    print(f"# disagreements={cen['disagree']} "
          f"(one-sided={cen['dis_one_sided']}, both-wrong={cen['dis_both_wrong']}) "
          f"agree_both_wrong={cen['agree_both_wrong']}")
    print(f"# RULE CENSUS on shared bases (before={cen['before_err']} error bases "
          f"over {cen['shared_bases']} mate-bases):")

    def row(label, s):
        prec = 100.0 * s["fixed"] / max(1, s["fixed"] + s["corrupted"])
        return (f"#   {label:<16} fired={s['fired']:<6} fixed={s['fixed']:<6} "
                f"corrupted={s['corrupted']:<6} after_err={s['after_err']:<6} "
                f"precision={prec:.1f}%")
    print(row("fastp-rule", cen["fastp"]))
    for thr in POSTERIOR_THRS:
        print(row(f"posterior@{thr}", cen["post"][thr]))

    # ---- ceiling sweep: how much error mass becomes reachable as the
    # minimum overlap required to call a correctable overlap drops (fastp
    # requires 30). Pure injection metadata; no tool runs involved.
    tot_inj = cen["injected_errors"]

    def shared_err(min_ovl):
        pairs = shared = err_in = 0
        for t in truth.values():
            ins = t["ins"]
            lo, hi = max(0, ins - L), min(ins, L)
            if hi - lo < min_ovl:
                continue
            pairs += 1
            shared += hi - lo
            for mate in (0, 1):
                for p, _, _ in t["err"][mate]:
                    rp = p if mate == 0 else ins - 1 - p
                    if lo <= rp < hi:
                        err_in += 1
        return pairs, shared, err_in

    sweep = [(m,) + shared_err(m) for m in (15, 20, 25, 30, 40, 60)]
    print("# CEILING SWEEP (min overlap -> correctable error mass):")
    print("#   min_ovl\tpairs\tshared_positions\terrors_in_shared\tpct_of_injected")
    for m, pairs, shared, err_in in sweep:
        print(f"#   {m}\t{pairs}\t{shared}\t{err_in}\t{100.0*err_in/max(1,tot_inj):.2f}")
    hist = []
    for lo_e, hi_e in ((1, 20), (20, 30), (30, 60), (60, 120), (120, 151)):
        pairs = err_in = 0
        for t in truth.values():
            ins = t["ins"]
            lo, hi = max(0, ins - L), min(ins, L)
            ovl = hi - lo
            if not (lo_e <= ovl < hi_e):
                continue
            pairs += 1
            for mate in (0, 1):
                for p, _, _ in t["err"][mate]:
                    rp = p if mate == 0 else ins - 1 - p
                    if lo <= rp < hi:
                        err_in += 1
        hist.append((f"{lo_e}-{hi_e - 1}", pairs, err_in))
    print("# ERROR MASS BY OVERLAP LENGTH:")
    for label, pairs, err_in in hist:
        print(f"#   ovl {label}bp\tpairs={pairs}\terrors_in_shared={err_in}")

    hdr = ["arm", "reads_in", "reads_out", "dropped", "bases_out", "frag_bases",
           "err_out", "frag_errs", "frag_err_per_kb", "corrupted", "edits",
           "mean_q_frag", "q30_frac", "edit_q_in", "edit_q_out",
           "mapped", "unmapped"]
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in hdr) + "\n")
    rules_path = out_tsv.replace(".tsv", "-rules.tsv")
    with open(rules_path, "w") as f:
        f.write("rule\tfired\tfixed\tcorrupted\tafter_err\tprecision\n")
        for label, s in [("fastp-rule", cen["fastp"])] + \
                        [(f"posterior@{t}", cen["post"][t]) for t in POSTERIOR_THRS]:
            prec = 100.0 * s["fixed"] / max(1, s["fixed"] + s["corrupted"])
            f.write(f"{label}\t{s['fired']}\t{s['fixed']}\t{s['corrupted']}\t"
                    f"{s['after_err']}\t{prec:.2f}\n")
    sweep_path = out_tsv.replace(".tsv", "-sweep.tsv")
    with open(sweep_path, "w") as f:
        f.write("min_ovl\tpairs\tshared_positions\terrors_in_shared\tpct_of_injected\n")
        for m, pairs, shared, err_in in sweep:
            f.write(f"{m}\t{pairs}\t{shared}\t{err_in}\t{100.0*err_in/max(1,tot_inj):.2f}\n")
    hist_path = out_tsv.replace(".tsv", "-ovl-hist.tsv")
    with open(hist_path, "w") as f:
        f.write("overlap_bp\tpairs\terrors_in_shared\n")
        for label, pairs, err_in in hist:
            f.write(f"{label}\t{pairs}\t{err_in}\n")
    print(f"# wrote {out_tsv}, {rules_path}, {sweep_path}, {hist_path}")
    print("\t".join(hdr))
    for r in rows:
        print("\t".join(str(r[k]) for k in hdr))
    # cross-check: arm-level net fix vs census emulation of fastp's rule
    by = {r["arm"]: r for r in rows}
    if "fastp-default" in by and "fastp-correct" in by:
        print(f"# cross-check: fastp --correction net fix at arm level = "
              f"{by['fastp-default']['err_out'] - by['fastp-correct']['err_out']} "
              f"error bases; census fastp-rule emulation fixed "
              f"{cen['fastp']['fixed']} (corrupted {cen['fastp']['corrupted']})")


if __name__ == "__main__":
    main()
