"""Materialise a fasthelper aggregation tree from ONE real report.

fasthelper's cost is a function of the number of sample directories and of the
bytes in their reports, not of sequencing data, so the tree is synthesised from
one genuine fastcheck v1 report (`qc.json` in every sample dir, with its
`sample` field rewritten) plus one genuine fastmaps grid + v1 wrapper on every
tenth sample, which exercises both the ordinary-report path and the grid path.

This is labelled synthetic in the report; it is not a claim about any library.
"""

import argparse
import json
import os
import shutil
import subprocess


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"{cmd[0]} failed: {r.stderr[-400:]}")
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--samples", type=int, required=True)
    ap.add_argument("--fastcheck", required=True)
    ap.add_argument("--fastmaps-bin", required=True)
    ap.add_argument("--reads-dir", required=True)
    args = ap.parse_args()

    out = args.out
    seed_dir = os.path.join(out, "_seed")
    os.makedirs(seed_dir, exist_ok=True)

    r1 = os.path.join(args.reads_dir, sorted(f for f in os.listdir(args.reads_dir) if f.endswith("_R1.fastq.gz"))[0])
    r2 = r1.replace("_R1.fastq.gz", "_R2.fastq.gz")

    report = os.path.join(seed_dir, "qc.json")
    if not os.path.exists(report):
        run([args.fastcheck, r1, r2, "--json", report, "-q", "-j", "2"])

    grid = os.path.join(seed_dir, "heatmap.json")
    wrapper = os.path.join(seed_dir, "heatmap.report.json")
    if not os.path.exists(grid):
        run([args.fastmaps_bin, r1, r2, "-o", seed_dir,
             "--json", grid, "--force", "-q"])

    with open(report) as fh:
        template = fh.read()

    for i in range(args.samples):
        name = f"lib{i:05d}"
        sample_dir = os.path.join(out, name)
        os.makedirs(sample_dir, exist_ok=True)
        # Rewrite only the sample field: the metrics stay realistic and the
        # aggregation still has to parse a full v1 document per sample.
        body = template.replace('"sample": "', f'"sample": "{name}', 1)
        with open(os.path.join(sample_dir, "qc.json"), "w") as fh:
            fh.write(body)
        if i % 10 == 0:
            shutil.copyfile(grid, os.path.join(sample_dir, "heatmap.json"))
            shutil.copyfile(wrapper, os.path.join(sample_dir, "heatmap.report.json"))

    # Remove the seed dir: it would otherwise be walked as a 100001th sample.
    shutil.rmtree(seed_dir, ignore_errors=True)
    print(f"helper_tree: {args.samples} samples under {out}")


if __name__ == "__main__":
    main()
