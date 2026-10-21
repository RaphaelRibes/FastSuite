#!/usr/bin/env python3
"""Task 2.5: is fastdedup's KEEP SET stable while its ORDER is not?

The plan asserts it is: "the KEEP SET is identical across modes (same
sequences, same counts; sorted-sequence MD5 f07fb28e..., 204 675 distinct
each), but the output order differs". Before documenting "order is not part of
the contract", verify both halves — a claim that the set is stable is exactly
the sort of thing that stops being true.

Runs the same input at several thread counts and compares:
  * raw file md5            (order-sensitive)
  * sorted-records md5      (order-insensitive)
  * record count and the set of sequence digests
"""
import hashlib
import gzip
import os
import shutil
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FD = os.path.join(REPO, "tools", "fastdedup", "target", "release", "fastdedup")
D = "/tmp/dedup-order"


def make_corpus(n_pairs=40000, seed=5):
    """Paired FASTQ with a controlled duplicate rate."""
    shutil.rmtree(D, ignore_errors=True)
    os.makedirs(D, exist_ok=True)
    src1 = os.path.join(D, "in_R1.fq.gz")
    src2 = os.path.join(D, "in_R2.fq.gz")
    if os.path.exists(src1) and os.path.exists(src2):
        return src1, src2
    x = seed
    with gzip.open(src1, "wt") as f1, gzip.open(src2, "wt") as f2:
        for i in range(n_pairs):
            # Every 5th pair is a duplicate of the previous one, so the dedup
            # path has real work and the keep set is a strict subset.
            idx = i - (1 if i % 5 == 4 else 0)
            for fh, mate in ((f1, 1), (f2, 2)):
                chars = []
                for _ in range(60):
                    x ^= (x << 13) & 0xFFFFFFFFFFFFFFFF
                    x ^= x >> 7
                    x ^= (x << 17) & 0xFFFFFFFFFFFFFFFF
                    chars.append("ACGT"[(x >> 5) & 3])
                # Deterministic by (idx, mate) so duplicate pairs match exactly.
                h = hashlib.blake2b(f"{idx}:{mate}".encode(), digest_size=8).digest()
                seq = "".join(chars[:20]) + "".join(
                    "ACGT"[h[k % 8] >> (k % 5) & 3] for k in range(40))
                # The SAME id on both mates: fastdedup reads a differing id as
                # a desynchronised pair, not as R1/R2 of one read.
                fh.write(f"@r{idx}\n{seq}\n+\n{'I' * 60}\n")
    return src1, src2


def read_records(path):
    """(id, seq) pairs, for both plain and gzip."""
    op = gzip.open if path.endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            out.append((h.rstrip("\n"), s))
    return out


def digest_rows(rows):
    h = hashlib.sha256()
    for r in rows:
        h.update(("\t".join(r) + "\n").encode())
    return h.hexdigest()


def file_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    n_pairs = int(sys.argv[1]) if len(sys.argv) > 1 else 40000
    if not os.path.exists(FD):
        raise SystemExit(f"build fastdedup first: {FD}")
    r1, r2 = make_corpus(n_pairs)
    runs = {}
    for j in (1, 2, 4, 8):
        out = os.path.join(D, f"j{j}")
        shutil.rmtree(out, ignore_errors=True)
        os.makedirs(out, exist_ok=True)
        o1 = os.path.join(out, "o1.fq")
        o2 = os.path.join(out, "o2.fq")
        r = subprocess.run(
            [FD, "-1", r1, "-2", r2, "-o", o1, "-p", o2, "-j", str(j),
             "--force", "-q"], capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr[-600:], file=sys.stderr)
            raise SystemExit(f"j{j} failed")
        rows = read_records(o1)
        runs[j] = dict(rows=rows, n=len(rows), sorted_digest=digest_rows(sorted(rows)),
                       raw_digest=digest_rows(rows), file_md5=file_md5(o1))
        print(f"j{j}: {len(rows):,} records  raw {runs[j]['raw_digest'][:16]}  "
              f"sorted {runs[j]['sorted_digest'][:16]}  file {runs[j]['file_md5'][:16]}")

    print()
    sorted_set = {v["sorted_digest"] for v in runs.values()}
    raw_set = {v["raw_digest"] for v in runs.values()}
    md5_set = {v["file_md5"] for v in runs.values()}
    counts = {v["n"] for v in runs.values()}

    print(f"distinct record COUNTS      : {counts}")
    print(f"distinct SORTED digests     : {len(sorted_set)}")
    print(f"distinct RAW-ORDER digests  : {len(raw_set)}")
    print(f"distinct file md5s          : {len(md5_set)}")
    print()
    if len(sorted_set) == 1:
        print("The KEEP SET is stable across -j.")
    else:
        print("THE KEEP SET ITSELF CHANGES across -j. 'Order is not part of the")
        print("contract' would NOT be enough to make these runs equivalent.")
    if len(raw_set) == 1 and len(md5_set) == 1:
        print("Output order was ALSO stable on this corpus -- the order difference")
        print("the plan recorded did not reproduce here. That does not license a")
        print("byte-identity gate: a property that holds on one corpus is not a")
        print("contract. See the scale note in the evidence doc.")
    else:
        print(f"Output order DIFFERS across -j ({len(raw_set)} distinct orderings),")
        print("reproducing the plan's finding: only the SET is contractual.")


if __name__ == "__main__":
    main()
