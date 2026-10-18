#!/usr/bin/env python3
"""Phase 0 of the PE-overlap plan: where does the corpus actually live, insert-wise?

fastprep's paired-end overlap mode can only act when the mates overlap, i.e. when
the insert is short enough for readthrough. fastsim's insert is measurable only
below 2L-12, and a pair above that is censored. So the mode's window is:

    L <= I <= 2L - 12        (overlap detectable AND measurable)
    I <  L                   (readthrough for sure; insert < read length)
    I > 2L - 12              (no shared anchor: nothing can be measured)

This script buckets every pair in every committed gt.tsv by that window and
reports the three bucket fractions per domain, plus how often a mate actually
carried adapter readthrough (r1_adapter>0 or r2_adapter>0) -- the condition the
mode is supposed to exploit.

Usage:  python3 insert_window_survey.py --root /mnt/Fast/fastsuite-scratch/work \
            --anchor 12 --out insert_window.tsv
"""
import argparse
import collections
import pathlib


def read_pairs(path: pathlib.Path):
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        need = ("insert", "r1_adapter", "r2_adapter", "read_len")
        missing = [c for c in need if c not in idx]
        if missing:
            # older gt files carry no read_len column; infer it is not needed
            # because we take L from the CLI argument.
            missing = [c for c in missing if c != "read_len"]
            if missing:
                raise SystemExit(f"{path}: missing columns {missing}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            try:
                yield int(f[idx["insert"]]), int(f[idx["r1_adapter"]]), int(f[idx["r2_adapter"]])
            except ValueError:
                continue


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="directory containing */<acc>/v0/*.gt.tsv")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--anchor", type=int, default=12, help="the 12-base anchor (fastcheck insert.rs)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    L = args.read_len
    ceiling = 2 * L - args.anchor
    root = pathlib.Path(args.root)

    files = sorted(root.glob("*/*/v0/*.gt.tsv"))
    if not files:
        raise SystemExit(f"no gt.tsv under {root}/*/*/v0/")

    # domain -> counter
    domains: dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    all_files = sorted(files)
    if not all_files:
        raise SystemExit(f"no gt.tsv under {root}/*/*/v0/")

    for path in all_files:
        domain = path.relative_to(root).parts[0]
        n = below = window = censored = readthrough = 0
        for insert, a1, a2 in read_pairs(path):
            n += 1
            if insert < L:
                below += 1
            elif insert <= ceiling:
                window += 1
            else:
                censored += 1
            if a1 > 0 or a2 > 0:
                readthrough += 1
        c = domains[domain]
        c["files"] += 1
        c["pairs"] += n
        c["below_L"] += below
        c["window"] += window
        c["censored"] += censored
        c["readthrough"] += readthrough

    rows = []
    for domain in sorted(domains):
        c = domains[domain]
        n = c["pairs"]
        rows.append(
            (domain, c["files"], n, c["below_L"], c["window"], c["censored"], c["readthrough"])
        )

    out = args.out
    lines = [
        "# Phase 0 survey: insert distribution per domain",
        f"# read_len={L}  anchor={args.anchor}  ceiling=2L-{args.anchor}={ceiling}",
        "# buckets: below_L (I<L), window (L<=I<=ceiling), censored (I>ceiling)",
        "domain\tfiles\tpairs\tbelow_L\twindow\tcensored\treadthrough\t"
        "pct_below_L\tpct_window\tpct_censored\tpct_readthrough",
    ]
    for domain, files, n, below, window, cens, rt in rows:
        pc = lambda x: (100.0 * x / n) if n else 0.0
        lines.append(
            f"{domain}\t{files}\t{n}\t{below}\t{window}\t{cens}\t{rt}\t"
            f"{pc(below):.2f}\t{pc(window):.2f}\t{pc(cens):.2f}\t{pc(rt):.2f}"
        )
    text = "\n".join(lines) + "\n"
    if out:
        pathlib.Path(out).write_text(text)
        print(f"wrote {out}")
    print(text)
    print(f"# files scanned: {len(all_files)}   ceiling: {ceiling} bp")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
