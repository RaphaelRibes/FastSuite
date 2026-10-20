"""Turn the bench's results.tsv into the before/after table each paper carries.

The papers' existing tables are rendered from their own committed CSVs by their
own `bench-table.typ`; this round's numbers are a different question (same tool,
same machine, same reads, two binary revisions), so they get their own table
rather than being spliced into a comparator table where the readings would no
longer be comparable.

Usage: python3 papers.py <results.tsv> [tool]
Prints Typst table markup on stdout; with no tool argument, prints one block per
tool.
"""

import sys
from collections import defaultdict

COLUMNS = ("tool", "kind", "reads", "threads", "gz_level", "wall_s", "cpu_s",
           "cpu_over_wall", "peak_rss_kb", "out_bytes", "files", "rc")


def load(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if tuple(header[:len(COLUMNS)]) != COLUMNS:
            raise SystemExit(f"unexpected header: {header}")
        rows = []
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(zip(COLUMNS, parts)))
        return rows


def fmt_reads(n):
    n = int(n)
    if n >= 1_000_000:
        return f"{n / 1_000_000:g} M"
    return f"{n / 1000:g} k"


def table_for(rows, tool):
    cells = defaultdict(dict)
    for r in rows:
        if r["tool"] != tool:
            continue
        key = (r["reads"], r["threads"], r["gz_level"])
        cells[key][r["kind"]] = r
    if not cells:
        return None
    lines = []
    lines.append("  table(")
    lines.append("    columns: (auto, auto, auto, auto, auto, auto),")
    lines.append("    table.header([*reads*], [*cores*], [*level*], "
                 "[*before, s*], [*after, s*], [*peak RSS before / after, MB*]),")
    for (reads, threads, gz) in sorted(cells, key=lambda k: (int(k[0]), int(k[1]), int(k[2] or 0))):
        pair = cells[(reads, threads, gz)]
        if "baseline" not in pair or "new" not in pair:
            continue
        b, n = pair["baseline"], pair["new"]
        rss_b = int(b["peak_rss_kb"]) / 1024
        rss_n = int(n["peak_rss_kb"]) / 1024
        lines.append(
            f"    [{fmt_reads(reads)}], [{threads}], [{gz or 'default'}], "
            f"[{float(b['wall_s']):.2f}], [{float(n['wall_s']):.2f}], "
            f"[{rss_b:.0f} / {rss_n:.0f}],"
        )
    lines.append("  ),")
    return "\n".join(lines)


def main():
    rows = load(sys.argv[1])
    tools = [sys.argv[2]] if len(sys.argv) > 2 else sorted({r["tool"] for r in rows})
    for tool in tools:
        block = table_for(rows, tool)
        if block:
            print(f"// {tool}")
            print(block)
            print()


if __name__ == "__main__":
    main()
