#!/usr/bin/env python3
"""Generate the FastSuite pathology-corpus fixtures (bench/fixtures/*.fastq).

Every fixture here is small, deterministic (fixed seed) and hand-checkable. Each
one exists to reproduce — at the level of a committed input, not a prose bug
report — one of the five defects the Phase D campaign found only by measuring
real sequencing data. See bench/fixtures/README.md for the mapping and for the
rationale.

Run:  python3 bench/fixtures/make_fixtures.py [OUTDIR]
Default OUTDIR is the directory this script lives in.

The generator is committed alongside the fixtures so the corpus is reproducible
and reviewable rather than opaque binary blobs. The fixtures it writes are the
bytes the integration test (bench/tests/pathology_test.sh) consumes.
"""
import os
import random
import sys

# Read length shared by most fixtures. Small enough to keep files tiny, long
# enough that a 12-mer window has room (defect 3) and that a boundary base has a
# stable position (defect 1).
L = 80
# Read length for the k-mer saturation fixture: 150 bp gives 139 windows per
# read, matching the profile the defect was measured on (defect 3).
KL = 150

ADAPTER = "CTGTCTCTTATACACATCT"


def write_fastq(path, records):
    """records: iterable of (name, seq, qual). Plain (uncompressed) FASTQ."""
    with open(path, "w") as f:
        for name, seq, qual in records:
            f.write("@%s\n%s\n+\n%s\n" % (name, seq, qual))


def rand_dna(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


# ---------------------------------------------------------------------------
# Defect 1 — fastsim 3' boundary rule.
#
# The guard was `t < read_len`, true at t == 0, so with NO tail artefact
# (tail = {0: N}) EVERY read's last base was forced non-G. Real libraries are
# ~23.5% G at the last position; the simulator produced 0.42%. The fix guards
# the rule with `t > 0`.
#
# This fixture is the PROFILE that drives fastsim into that regime: a tail
# histogram that is entirely 0 (no tail artefact) with a GC-rich genome, so the
# last biology base is very often G. The regression test simulates reads from
# this profile and asserts the last-position G fraction stays near the genome GC
# (~25%), not pinned to ~0%.
def make_notail_profile(path):
    import json

    prof = {
        "read_len": L,
        "gc": 0.70,  # GC-rich so ~35% of last bases are G before the bug
        "geometry": "tail3p_both",
        "head": {"R1": {"0": 100}, "R2": {"0": 100}},
        "tail": {"R1": {"0": 100}, "R2": {"0": 100}},  # NO tail artefact
        "joint_tail": {"0,0": 100},
        "tail_q_collapse": {"R1": 0.0, "R2": 0.0},
        "pos_q": {"R1": [30] * L, "R2": [30] * L},
        "read_meanq": {"R1": {"30": 100}, "R2": {"30": 100}},
        "full_length_rate": 0.0,
        "full_length_q": {"R1": 20, "R2": 20},
        "adapter": {"seq": ADAPTER, "mode": "none"},
    }
    with open(path, "w") as f:
        json.dump(prof, f, indent=1)
    return prof


# ---------------------------------------------------------------------------
# Defect 2 — fastcheck profile dup_rate used accs.first() (first mate only)
# while the report uses the merged R1+R2 fingerprint table. On DRR000466 the two
# mates carried different duplicate rates (R1 0.0086, R2 0.0182, merged 0.0135),
# so the profile disagreed with the report by construction.
#
# The fix computes dup_rate over the MERGED accumulator. This fixture is a
# paired input whose two mates have very different duplicate rates, so a
# first-mate-only profile and a merged profile give visibly different numbers.
def make_paired_dup(outdir):
    rng = random.Random(11)
    r1, r2 = [], []
    # 300 unique pairs: R1 and R2 each carry a distinct sequence (shared within a
    # pair only). R2's duplicate class stays near 0 for this block.
    for _ in range(300):
        a = rand_dna(rng, L)
        b = rand_dna(rng, L)
        r1.append(a)
        r2.append(b)
    # 300 more pairs whose R2 is the SAME sequence repeated -> R2 duplicate rate
    # ~0.5 while R1 stays diverse. Merged rate lands ~0.5 too, but a first-mate
    # (R1) reading would be ~0.0 — the exact disagreement the defect produced.
    hot = rand_dna(rng, L)
    for _ in range(300):
        r1.append(rand_dna(rng, L))
        r2.append(hot)
    write_fastq(os.path.join(outdir, "paired_dup_R1.fastq"),
                [("p%d" % i, s, "I" * L) for i, s in enumerate(r1)])
    write_fastq(os.path.join(outdir, "paired_dup_R2.fastq"),
                [("p%d" % i, s, "I" * L) for i, s in enumerate(r2)])
    return len(r1)


# ---------------------------------------------------------------------------
# Defect 3 — fastcheck k-mer counter ran an O(n) singleton eviction on EVERY
# insert once the map hit its 1M cap. On real data the 12-mer space saturates
# and every insert paid a 1M-entry rehash (40 min of CPU on DRR509062, frozen
# RSS). The fix evicts lazily (amortized O(1)).
#
# NOTE (verified, and stated in the README): this pathology is NOT reproducible
# at the fixture level. A synthetic stream that saturates the 1M cap finishes in
# well under a second even with the pre-fix every-insert eviction, because the
# synthetic key distribution never keeps a large count>1 population resident at
# the cap the way DRR509062's real 12-mer spectrum does. The fixture below is
# still committed: it saturates the 12-mer space (a necessary condition of the
# defect) and the test asserts the post-fix invariant — the scan completes fast
# and emits a k-mer spectrum. It is a guard against the amortization regressing
# into a per-insert scan, not a reproduction of the 40-minute hang.
def make_kmer_saturation(path):
    rng = random.Random(5)
    with open(path, "w") as f:
        for i in range(8000):
            # Fully random 150 bp reads: 139 overlapping 12-mer windows each,
            # ~1.076M DISTINCT 12-mers total — past the 1M cap, so the counter
            # saturates exactly as the defect requires.
            seq = rand_dna(rng, KL)
            f.write("@k%d\n%s\n+\n%s\n" % (i, seq, "I" * KL))


# ---------------------------------------------------------------------------
# Defect 4 — Phred+64 libraries decoded as Phred+33. Real mean_q of 62-68 is
# physically impossible; two campaign datasets (DRR014039, DRR006305) carried
# quality bytes 66-105. Decoding them as +33 gave a mean Q near 52 instead of
# the plausible ~20-40 band.
#
# fastcheck does NOT guess +64 from the byte census (a census entirely >= 64 is
# arithmetically consistent with BOTH offsets, so guessing would be a lie); it
# reports the ambiguity as an assumption and `--phred-offset 64` settles it. The
# fixture is a FASTQ whose quality bytes sit in the 66-105 range. The test runs
# fastcheck twice: once on auto (asserts the +33/+64 ambiguity is *declared*,
# certain == 0) and once with --phred-offset 64 (asserts the mean Q lands in the
# plausible 20-40 band and the forced offset is reported consistent).
def make_phred64(path):
    rng = random.Random(23)
    with open(path, "w") as f:
        for i in range(2000):
            seq = rand_dna(rng, L)
            # Real Q values 2..41, encoded as Phred+64 -> bytes 66..105.
            qual = "".join(chr(64 + rng.randint(2, 41)) for _ in range(L))
            f.write("@r%d\n%s\n+\n%s\n" % (i, seq, qual))


# ---------------------------------------------------------------------------
# Defect 5 — N-heavy padding dominated the duplicate count. Real libraries carry
# N-runs that are the MOST duplicated sequences (NNNN... x276 on one dataset),
# while the simulator emits zero Ns (its only Ns are adapter-None padding). The
# fingerprint counts every sequence including N-runs, so a corpus full of
# identical N-runs reports them as the top duplicated sequence.
#
# The fix's regression value is to lock in that N-runs ARE counted as duplicates
# (not filtered, not special-cased away): the duplication_level must reflect
# them and the most-repeated sequence must be the N-run. If a future change
# silently dropped non-ACGT sequences from the fingerprint, this fixture fails.
def make_nruns(path):
    rng = random.Random(7)
    records = []
    # 400 identical 80 bp N-runs: the dominant duplicated sequence.
    for i in range(400):
        records.append(("nrun%d" % i, "N" * L, "I" * L))
    # 600 distinct ACGT reads so the level is a real mixture, not 1.0.
    for i in range(600):
        records.append(("u%d" % i, rand_dna(rng, L), "I" * L))
    write_fastq(path, records)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    os.makedirs(outdir, exist_ok=True)

    make_notail_profile(os.path.join(outdir, "notail_profile.json"))
    n_pairs = make_paired_dup(outdir)
    make_kmer_saturation(os.path.join(outdir, "kmer_saturation.fastq"))
    make_phred64(os.path.join(outdir, "phred64.fastq"))
    make_nruns(os.path.join(outdir, "nruns.fastq"))

    total = 0
    for name in sorted(os.listdir(outdir)):
        p = os.path.join(outdir, name)
        if os.path.isfile(p) and name != "make_fixtures.py":
            sz = os.path.getsize(p)
            total += sz
            print("  %-28s %8d B" % (name, sz))
    print("  %-28s %8d B" % ("TOTAL", total))
    print("wrote fixtures to %s (%d pairs in paired_dup)" % (outdir, n_pairs))


if __name__ == "__main__":
    main()
