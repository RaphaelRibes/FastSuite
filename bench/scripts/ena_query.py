#!/usr/bin/env python3
"""Build the FastSuite corpus manifest from the ENA portal API (PLAN-ML-LOCAL B1).

Queries the ENA portal search API, stratified by domain, keeping one row per
read run with fastq_ftp + fastq_md5 + read_count (plus fastq_bytes for a
sanity bound in B2). Output: manifests/real.tsv

Deterministic and re-runnable: domains are queried in a fixed order, rows are
sorted by accession within each domain, and an identical query set produces a
byte-identical manifest (ENA returns stable ordering for a fixed sort key).

Usage:
  python3 ena_query.py                     # build manifests/real.tsv (100/domain)
  python3 ena_query.py --stats             # per-domain sanity stats, no query
  python3 ena_query.py --per-domain 50 -o manifests/small.tsv
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

API = "https://www.ebi.ac.uk/ena/portal/api/search"
RESULT = "read_run"
FIELDS = "run_accession,fastq_ftp,fastq_md5,fastq_bytes,read_count,library_strategy,sample_title"
FORMAT = "tsv"
TIMEOUT = 60
MAX_RETRIES = 4

# Fixed domain order (deterministic output).
DOMAINS = [
    # generalised: broad bacterial WGS pool (no human taxon filter, no niche
    # filter) — the "anything goes" reference population.
    ("generalised", 'library_strategy="WGS"'),
    # metagenomics: WGS from a metagenomic library source.
    ("metagenomics", 'library_strategy="WGS" AND library_source="METAGENOMIC"'),
    # human gut: WGS of samples titled gut/stool/feces/fecal.
    ("human_gut", 'library_strategy="WGS" AND (sample_title="*gut*" OR sample_title="*stool*" OR sample_title="*feces*" OR sample_title="*fecal*")'),
    # virome: WGS of samples titled virome/viral (host-read depletion style).
    ("virome", 'library_strategy="WGS" AND (sample_title="*virome*" OR sample_title="*viral*")'),
    # bacterial isolate: WGS from a genomic library source (isolates).
    ("bacterial_isolate", 'library_strategy="WGS" AND library_source="GENOMIC"'),
    # rna_seq: transcriptomics.
    ("rna_seq", 'library_strategy="RNA-Seq"'),
]

# Deterministic, stable ordering + a hard limit well above the per-domain
# target so post-filters (dropping rows without ftp/md5/read_count) still
# leave >= target rows.
OVERSAMPLE = 3

HEADER = ["domain", "run_accession", "fastq_ftp", "fastq_md5", "fastq_bytes", "read_count", "library_strategy", "sample_title"]


def query_ena(query: str, limit: int) -> list[list[str]]:
    """One ENA portal search call -> rows (list of field values), retries on 5xx/network."""
    params = {
        "result": RESULT,
        "query": query,
        "fields": FIELDS,
        "format": FORMAT,
        "limit": str(limit),
        "sort": "run_accession",
    }
    url = API + "?" + urllib.parse.urlencode(params)
    last_err: Exception | None = None
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "fastsuite-ena-query/1.0"})
            with urllib.request.urlopen(req, timeout=TIMEOUT) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            rows = [line.split("\t") for line in text.splitlines() if line.strip()]
            # First row is the header the API echoes back.
            if rows and rows[0][0] == "run_accession":
                rows = rows[1:]
            return rows
        except Exception as exc:  # noqa: BLE001 - retry any transport/server error
            last_err = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"ENA query failed after {MAX_RETRIES} attempts: {last_err}")


def complete_row(row: list[str]) -> bool:
    """A usable row needs a run accession, an FTP path, an md5 and a read count."""
    if len(row) < len(FIELDS.split(",")):
        return False
    acc, ftp, md5, _bytes, reads = row[:5]
    return bool(acc and ftp and md5 and reads) and all(p.startswith("ftp.") for p in ftp.split(";") if p)


def build_manifest(per_domain: int, out_path: str) -> dict[str, int]:
    seen: dict[str, str] = {}  # accession -> domain (global dedup across domains)
    rows_out: list[dict[str, str]] = []
    counts: dict[str, int] = {}

    for domain, query in DOMAINS:
        raw = query_ena(query, limit=per_domain * OVERSAMPLE)
        kept = 0
        for row in raw:
            if not complete_row(row):
                continue
            acc = row[0]
            if acc in seen:  # global dedup: a run keeps its first domain only
                continue
            seen[acc] = domain
            rows_out.append(dict(zip(HEADER, [domain] + row)))
            kept += 1
            if kept >= per_domain:
                break
        counts[domain] = kept
        print(f"  {domain:18s} {kept:3d}/{per_domain} rows", file=sys.stderr)

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp = out_path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER, delimiter="\t")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    os.replace(tmp, out_path)
    return counts


def load_manifest(path: str) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def stats(path: str) -> None:
    """Per-domain sanity: count, byte spread, read-count spread, library strategies."""
    rows = load_manifest(path)
    by_domain: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        by_domain[r["domain"]].append(r)

    print(f"manifest: {path}")
    print(f"total rows: {len(rows)}  accessions: {len({r['run_accession'] for r in rows})}")
    print()
    print(f"{'domain':18s} {'n':>4s} {'bytes p10/med/p90':>22s} {'reads p10/med/p90':>22s} strategies")
    for domain, _ in DOMAINS:
        rs = by_domain.get(domain, [])
        if not rs:
            print(f"{domain:18s} {'0':>4s}  (none)")
            continue
        sizes = sorted(int(r["fastq_bytes"]) for r in rs if r["fastq_bytes"].isdigit())
        counts = sorted(int(r["read_count"]) for r in rs if r["read_count"].isdigit())
        strategies = sorted({r["library_strategy"] for r in rs if r["library_strategy"]})

        def pct(xs, p):
            return xs[min(len(xs) - 1, int(p * len(xs)))] if xs else 0

        bs = f"{pct(sizes, .1)/1e6:.0f}/{pct(sizes, .5)/1e6:.0f}/{pct(sizes, .9)/1e6:.0f}M"
        cs = f"{pct(counts, .1)/1e6:.1f}/{pct(counts, .5)/1e6:.1f}/{pct(counts, .9)/1e6:.1f}M"
        print(f"{domain:18s} {len(rs):4d} {bs:>22s} {cs:>22s} {','.join(strategies)}")


def verify(path: str, per_domain: int) -> int:
    """Hard verification: floors, uniqueness, completeness. Exit code 0 = pass."""
    rows = load_manifest(path)
    errs: list[str] = []
    accs = [r["run_accession"] for r in rows]
    dupes = {a for a in accs if accs.count(a) > 1}
    if dupes:
        errs.append(f"duplicate accessions: {sorted(dupes)[:10]}")
    for r in rows:
        missing = [f for f in ("fastq_ftp", "fastq_md5", "read_count") if not r.get(f)]
        if missing:
            errs.append(f"{r['run_accession']}: missing {','.join(missing)}")
    by_domain: dict[str, int] = defaultdict(int)
    for r in rows:
        by_domain[r["domain"]] += 1
    for domain, _ in DOMAINS:
        n = by_domain.get(domain, 0)
        if n < per_domain:
            errs.append(f"domain {domain}: {n} rows < floor {per_domain}")

    print(f"verify: {len(rows)} rows, {len(accs)} unique accessions")
    for domain, _ in DOMAINS:
        print(f"  {domain:18s} {by_domain.get(domain, 0):4d}")
    if errs:
        print(f"\nFAIL ({len(errs)} problems):")
        for e in errs[:20]:
            print(f"  - {e}")
        return 1
    print("\nOK: floors met, no duplicates, every row complete")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-domain", type=int, default=100, help="datasets per domain (default 100)")
    ap.add_argument("-o", "--out", default=os.path.join("manifests", "real.tsv"), help="output TSV path")
    ap.add_argument("--stats", action="store_true", help="print per-domain stats for an existing manifest and exit")
    ap.add_argument("--verify", action="store_true", help="verify an existing manifest and exit")
    args = ap.parse_args()

    if args.stats:
        stats(args.out)
        return 0
    if args.verify:
        return verify(args.out, args.per_domain)

    print(f"building manifest -> {args.out} ({args.per_domain}/domain)", file=sys.stderr)
    build_manifest(args.per_domain, args.out)
    print(f"wrote {args.out}", file=sys.stderr)
    print()
    return verify(args.out, args.per_domain)


if __name__ == "__main__":
    sys.exit(main())
