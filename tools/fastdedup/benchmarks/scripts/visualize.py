#!/usr/bin/env python3
"""
Visualize FastDedup benchmark results from benchmark.tsv.

Generates PNG graphs for:
  1. Wall time scaling (threads × variant)
  2. Peak memory scaling (threads × variant)
  3. I/O operations by mode
  4. SE vs PE comparison
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# Dark theme
plt.style.use("dark_background")
COLORS = {"normal": "#00bfff", "abundance": "#ff8c00", "abundance_tempdir": "#32cd32"}
LINE_STYLES = {"SE": "-", "PE": "--"}
LABEL_MAP = {
    "100000": "100K", "1000000": "1M", "10000000": "10M",
    "100K": "100K", "1M": "1M", "10M": "10M",
}


def load_data(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    return df


def plot_wall_time(df: pd.DataFrame, output_dir: str):
    """Wall time vs threads, faceted by dataset × io_mode."""

    datasets = sorted(df["dataset"].unique())
    io_modes = df["io"].unique()
    variants = df["variant"].unique()
    threads = sorted(df["threads"].unique())

    fig, axes = plt.subplots(
        len(datasets), len(io_modes),
        figsize=(5 * len(io_modes), 4 * len(datasets)),
        sharex="col", sharey="row",
    )
    if len(datasets) == 1:
        axes = axes.reshape(1, -1)
    if len(io_modes) == 1:
        axes = axes.reshape(-1, 1)

    for row, ds in enumerate(datasets):
        for col, io in enumerate(io_modes):
            ax = axes[row][col]
            sub = df[(df["dataset"] == ds) & (df["io"] == io)]

            for mode in ["SE", "PE"]:
                sub_mode = sub[sub["mode"] == mode]
                for var in variants:
                    sub_var = sub_mode[sub_mode["variant"] == var]
                    if sub_var.empty:
                        continue
                    sub_var = sub_var.sort_values("threads")
                    ax.plot(
                        sub_var["threads"].astype(int),
                        sub_var["wall_sec"],
                        marker="o", markersize=4,
                        linestyle=LINE_STYLES[mode],
                        color=COLORS[var],
                        label=f"{var} {mode}" if mode == "SE" and row == 0 and col == 0 else "",
                    )

            ax.set_title(f"{LABEL_MAP.get(ds, ds)} — {io}")
            ax.set_xlabel("Threads")
            ax.set_ylabel("Wall time (s)")
            ax.grid(alpha=0.15)
            ax.set_xticks(threads)

    # Single legend for the whole figure
    handles = []
    for var in variants:
        for mode in ["SE", "PE"]:
            h = plt.Line2D([], [], color=COLORS[var], linestyle=LINE_STYLES[mode],
                           label=f"{var} ({mode})")
            handles.append(h)
    fig.legend(handles=handles, loc="upper right", fontsize=9)

    fig.suptitle("Wall Time Scaling", fontsize=14, y=1.02)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "wall_time.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_memory(df: pd.DataFrame, output_dir: str):
    """Peak RSS vs threads, faceted by dataset × io_mode."""

    datasets = sorted(df["dataset"].unique())
    io_modes = df["io"].unique()
    variants = df["variant"].unique()
    threads = sorted(df["threads"].unique())

    fig, axes = plt.subplots(
        len(datasets), len(io_modes),
        figsize=(5 * len(io_modes), 4 * len(datasets)),
        sharex="col", sharey="row",
    )
    if len(datasets) == 1:
        axes = axes.reshape(1, -1)
    if len(io_modes) == 1:
        axes = axes.reshape(-1, 1)

    for row, ds in enumerate(datasets):
        for col, io in enumerate(io_modes):
            ax = axes[row][col]
            sub = df[(df["dataset"] == ds) & (df["io"] == io)]

            for mode in ["SE", "PE"]:
                sub_mode = sub[sub["mode"] == mode]
                for var in variants:
                    sub_var = sub_mode[sub_mode["variant"] == var]
                    if sub_var.empty:
                        continue
                    sub_var = sub_var.sort_values("threads")
                    ax.plot(
                        sub_var["threads"].astype(int),
                        sub_var["max_rss_kb"] / 1024,  # Convert to MB
                        marker="s", markersize=4,
                        linestyle=LINE_STYLES[mode],
                        color=COLORS[var],
                        label=f"{var} {mode}" if mode == "SE" and row == 0 and col == 0 else "",
                    )

            ax.set_title(f"{LABEL_MAP.get(ds, ds)} — {io}")
            ax.set_xlabel("Threads")
            ax.set_ylabel("Peak RSS (MB)")
            ax.grid(alpha=0.15)
            ax.set_xticks(threads)

    handles = []
    for var in variants:
        for mode in ["SE", "PE"]:
            h = plt.Line2D([], [], color=COLORS[var], linestyle=LINE_STYLES[mode],
                           label=f"{var} ({mode})")
            handles.append(h)
    fig.legend(handles=handles, loc="upper right", fontsize=9)

    fig.suptitle("Peak Memory Usage (RSS)", fontsize=14, y=1.02)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "memory.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_io(df: pd.DataFrame, output_dir: str):
    """I/O blocks vs threads, one plot per dataset, grouped by io_mode."""

    datasets = sorted(df["dataset"].unique())

    fig, axes = plt.subplots(
        1, len(datasets),
        figsize=(6 * len(datasets), 5),
        sharey=True,
    )
    if len(datasets) == 1:
        axes = [axes]

    io_colors = {"u2u": "#00bfff", "c2u": "#ff8c00", "c2c": "#32cd32"}

    for idx, ds in enumerate(datasets):
        ax = axes[idx]
        sub = df[(df["dataset"] == ds) & (df["variant"] == "normal") & (df["mode"] == "SE")]

        # Group by io_mode, sum I/O blocks
        io_groups = sub.groupby("io")[["io_in_blocks", "io_out_blocks"]].sum()
        x = range(len(io_groups))
        width = 0.35

        ax.bar([i - width / 2 for i in x], io_groups["io_in_blocks"] / 1e6,
               width, label="Read blocks (M)", color=io_colors.values(), alpha=0.7)
        ax.bar([i + width / 2 for i in x], io_groups["io_out_blocks"] / 1e6,
               width, label="Write blocks (M)", color=io_colors.values(), alpha=0.3)

        ax.set_title(LABEL_MAP.get(ds, ds))
        ax.set_xticks(list(x))
        ax.set_xticklabels(io_groups.index)
        ax.set_ylabel("I/O blocks (M)")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.15)

    fig.suptitle("Filesystem I/O by Input/Output Mode", fontsize=14)
    fig.tight_layout()
    out_path = os.path.join(output_dir, "io.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


def plot_variant_comparison(df: pd.DataFrame, output_dir: str):
    """Direct comparison of variants for a representative config."""

    # Pick 100K u2u SE as a representative config
    sub = df[
        (df["dataset"] == df["dataset"].unique()[0])
        & (df["io"] == "u2u")
        & (df["mode"] == "SE")
    ]

    if sub.empty:
        print("  Skipping variant comparison — no data for representative config")
        return

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Wall time
    for var in sub["variant"].unique():
        s = sub[sub["variant"] == var].sort_values("threads")
        axes[0].plot(s["threads"].astype(int), s["wall_sec"],
                     marker="o", color=COLORS[var], label=var)
    axes[0].set_title("Wall Time (100K u2u SE)")
    axes[0].set_xlabel("Threads")
    axes[0].set_ylabel("Seconds")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.15)

    # Memory
    for var in sub["variant"].unique():
        s = sub[sub["variant"] == var].sort_values("threads")
        axes[1].plot(s["threads"].astype(int), s["max_rss_kb"] / 1024,
                     marker="s", color=COLORS[var], label=var)
    axes[1].set_title("Peak RSS (100K u2u SE)")
    axes[1].set_xlabel("Threads")
    axes[1].set_ylabel("MB")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.15)

    # I/O
    for var in sub["variant"].unique():
        s = sub[sub["variant"] == var].sort_values("threads")
        axes[2].plot(s["threads"].astype(int), (s["io_in_blocks"] + s["io_out_blocks"]) / 1e6,
                     marker="^", color=COLORS[var], label=var)
    axes[2].set_title("Total I/O Blocks (100K u2u SE)")
    axes[2].set_xlabel("Threads")
    axes[2].set_ylabel("Total blocks (M)")
    axes[2].legend(fontsize=8)
    axes[2].grid(alpha=0.15)

    fig.tight_layout()
    out_path = os.path.join(output_dir, "variant_comparison.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate benchmark visualization graphs")
    parser.add_argument("--tsv", default=None,
                        help="Path to benchmark.tsv (default: benchmarks/benchmark.tsv)")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory for PNG files (default: benchmarks/)")
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)

    tsv_path = args.tsv or os.path.join(base_dir, "benchmark.tsv")
    output_dir = args.output_dir or base_dir

    if not os.path.exists(tsv_path):
        print(f"❌ {tsv_path} not found — run the Snakemake pipeline first")
        sys.exit(1)

    print(f"Loading {tsv_path} ...")
    df = load_data(tsv_path)
    print(f"  {len(df)} rows loaded")
    print(f"  Datasets: {sorted(df['dataset'].unique())}")
    print(f"  I/O modes: {list(df['io'].unique())}")
    print(f"  Modes: {list(df['mode'].unique())}")
    print(f"  Variants: {list(df['variant'].unique())}")

    os.makedirs(output_dir, exist_ok=True)
    print(f"\nGenerating graphs in {output_dir} ...")

    plot_wall_time(df, output_dir)
    plot_memory(df, output_dir)
    plot_io(df, output_dir)
    plot_variant_comparison(df, output_dir)

    print(f"\n✅ Graphs saved to {output_dir}")


if __name__ == "__main__":
    main()
