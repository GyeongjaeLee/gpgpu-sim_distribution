#!/usr/bin/env python3
"""
Plot per-direction link stats from AccelSim results CSV files as 3D bar charts.

Reads *_directions.csv (or *_nm_directions.csv with --near) and the
corresponding *_summary.csv files produced by parse_accelsim_results.py.

Supply one or more directions CSV files directly as positional arguments.
The script finds the matching *_summary.csv automatically (same prefix).

Layout (within one figure per structure):
  rows    = bandwidth × benchmark  (bw outer, benchmark inner)
  columns = routing scheme

Each cell is one 3D bar chart: Src axis × Dst axis, bar height = usage % or avg_sat %.
Separate figures are produced for each structure.

──────────────────────────────────────────────────────────────────────────────
Usage examples
──────────────────────────────────────────────────────────────────────────────

# 1. Single benchmark, selected structures/BWs/routings, save as PNG
python plot_accelsim_link_stats.py \
    bfs_BG_H3_bas+mina_directions.csv \
    --structure B100_Global \
    --bandwidth B200+HBM3e \
    --routing baseline near_min_adaptive near_min_random \
    --near-min-k 2 --near-min-p 1.0 \
    --metric both \
    -o bfs_link.png

# 2. Multiple benchmarks from separate CSVs, specific routing
python plot_accelsim_link_stats.py \
    bfs_BG_H3_bas_directions.csv \
    gemm_BG_H3_bas_directions.csv \
    --routing fixed_min \
    --metric util \
    -o multi_link_util.png

# 3. Near-min direction breakdown (use _nm_directions.csv)
python plot_accelsim_link_stats.py \
    bfs_RU_H4_nm_directions.csv \
    --near \
    --metric sat \
    -o bfs_nearmin_sat.png

# 4. All configs from CSV, both metrics, interactive display
python plot_accelsim_link_stats.py \
    results_directions.csv \
    --all-configs
"""

import argparse
import csv
import os
import re
from collections import defaultdict
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

from experiments_loader import load_experiments

# ── Experiment tables ─────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
STRUCTURES, BANDWIDTHS_DICT = load_experiments(os.path.join(_HERE, "experiments.csv"))
BANDWIDTHS = list(BANDWIDTHS_DICT.keys())

ROUTING_CHOICES = [
    "baseline", "min_oblivious", "min_adaptive",
    "near_min_adaptive", "near_min_random", "fixed_min", 
    "ugal", "valiant",
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def routing_to_key(routing: str, near_min_k: Optional[int] = None, near_min_p: Optional[float] = None) -> str:
    """Convert routing name (+ optional K and P values) to a directory/CSV key."""
    if routing == "near_min_adaptive":
        k = near_min_k if near_min_k is not None else 2
        p = near_min_p if near_min_p is not None else 1.0
        return f"{routing}_nmk{k}_nmp{p:.1f}"
    elif routing == "near_min_random":
        k = near_min_k if near_min_k is not None else 2
        return f"{routing}_nmk{k}"
    return routing


def canonical_router_name(name: str) -> str:
    m = re.match(r'^(Xbar\d+|MC\d+)', name)
    return m.group(1) if m else name


def short_bench(bench: str) -> str:
    """rodinia-3.1:bfs-rodinia-3.1 → bfs"""
    app = bench.split(":")[-1] if ":" in bench else bench
    return app.split("-")[0]


def short_routing(rk: str) -> str:
    """Convert long routing keys into compact chart titles."""
    abbrev = {
        "baseline":      "bas",
        "min_oblivious": "mino",
        "min_adaptive":  "mina",
        "fixed_min":     "fixm",
        "ugal":          "ug",
        "valiant":       "val",
    }
    if rk in abbrev:
        return abbrev[rk]
        
    m_adp = re.match(r"^near_min_adaptive_nmk(\d+)_nmp([\d.]+)$", rk)
    if m_adp:
        # e.g., near_min_adaptive_nmk2_nmp1.0 -> nma_k2_p10
        return f"nma_k{m_adp.group(1)}_p{m_adp.group(2).replace('.', '')}"
        
    m_rnd = re.match(r"^near_min_random_nmk(\d+)$", rk)
    if m_rnd:
        # e.g., near_min_random_nmk2 -> nmr_k2
        return f"nmr_k{m_rnd.group(1)}"
        
    return rk[:5]


# ── CSV loading ───────────────────────────────────────────────────────────────

def load_directions(path: str) -> dict:
    """
    Returns dict: (benchmark, struct, bw, routing) → list of
        {'src', 'dst', 'count', 'pct', 'avg_sat'} dicts.
    """
    data: dict = defaultdict(list)
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            key = (row["benchmark"], row["struct"], row["bw"], row["routing"])
            data[key].append({
                "src":     row["src"],
                "dst":     row["dst"],
                "count":   int(row["count"]),
                "pct":     float(row["pct"]),
                "avg_sat": float(row["avg_sat"]),
            })
    return dict(data)


def load_nm_directions(path: str) -> dict:
    """
    Returns dict: (benchmark, struct, bw, routing) → list of
        {'src', 'dst', 'link_type', 'count', 'pct', 'avg_sat'} dicts.
    """
    data: dict = defaultdict(list)
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            key = (row["benchmark"], row["struct"], row["bw"], row["routing"])
            data[key].append({
                "src":       row["src"],
                "dst":       row["dst"],
                "link_type": row["link_type"],
                "count":     int(row["count"]),
                "pct":       float(row["pct"]),
                "avg_sat":   float(row["avg_sat"]),
            })
    return dict(data)


def load_summary(path: str) -> dict:
    """
    Returns dict: (benchmark, struct, bw, routing) → summary row dict.
    """
    data: dict = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            key = (row["benchmark"], row["struct"], row["bw"], row["routing"])
            data[key] = row
    return data


# ── Router list ───────────────────────────────────────────────────────────────

def build_router_list(dir_entries_list: list) -> list:
    names = set()
    for entries in dir_entries_list:
        for e in entries:
            names.add(canonical_router_name(e["src"]))
            names.add(canonical_router_name(e["dst"]))
    xbars = sorted([n for n in names if n.startswith("Xbar")],
                   key=lambda x: int(re.search(r"\d+", x).group()))
    mcs   = sorted([n for n in names if n.startswith("MC")],
                   key=lambda x: int(re.search(r"\d+", x).group()))
    return xbars + mcs


# ── Annotation text ───────────────────────────────────────────────────────────

def make_annotation(summary_row: Optional[dict], mode: str) -> str:
    if summary_row is None:
        return ""
    lines = []

    if mode == "util":
        parts = []
        for lt, short in [("XBAR_XBAR", "X-X"), ("XBAR_MC", "X-M"), ("MC_MC", "M-M")]:
            pct_key = f"{lt}_pct"
            if pct_key in summary_row and summary_row[pct_key]:
                try:
                    v = float(summary_row[pct_key])
                    if v > 0:
                        parts.append(f"{short}={v:.1f}%")
                except ValueError:
                    pass
        if parts:
            lines.append("Util: " + " ".join(parts))
    else:  # sat
        parts = []
        for lt, short in [("XBAR_XBAR", "X-X"), ("XBAR_MC", "X-M"), ("MC_MC", "M-M")]:
            sat_key = f"{lt}_sat"
            if sat_key in summary_row and summary_row[sat_key]:
                try:
                    v = float(summary_row[sat_key])
                    parts.append(f"{short}={v:.1f}%")
                except ValueError:
                    pass
        if parts:
            lines.append("Sat: " + " ".join(parts))

    # L2 hit rate + cycles
    perf = []
    try:
        lr = summary_row.get("l2_hit_rate", "")
        if lr:
            perf.append(f"L2hit={float(lr)*100:.1f}%")
    except ValueError:
        pass
    try:
        cyc = summary_row.get("tot_cycle", "")
        if cyc:
            perf.append(f"Cyc={int(cyc):,}")
    except ValueError:
        pass
    if perf:
        lines.append(" ".join(perf))

    # Near-min ratio (present for near_min_adaptive/random)
    try:
        nm = summary_row.get("nm_ratio", "")
        if nm:
            lines.append(f"NM-ratio={float(nm):.1f}%")
    except ValueError:
        pass

    return "\n".join(lines)


# ── 3D bar chart ──────────────────────────────────────────────────────────────

def plot_3d_bars(ax, dir_entries, router_list, value_key, vmax, cmap, norm):
    n = len(router_list)
    dst_list    = list(reversed(router_list))
    src_to_idx  = {name: i for i, name in enumerate(router_list)}
    dst_to_idx  = {name: i for i, name in enumerate(dst_list)}

    grid = np.full((n, n), np.nan)
    for e in dir_entries:
        src = canonical_router_name(e["src"])
        dst = canonical_router_name(e["dst"])
        if src in src_to_idx and dst in dst_to_idx:
            grid[src_to_idx[src], dst_to_idx[dst]] = e[value_key]

    for si in range(n):
        for di in range(n):
            val = grid[si, di]
            if np.isnan(val):
                continue
            ax.bar3d(si, di, 0, 0.6, 0.6, val,
                     color=cmap(norm(val)), alpha=0.88,
                     edgecolor="k", linewidth=0.2)

    ax.set_xticks(np.arange(n) + 0.3)
    ax.set_xticklabels(router_list, rotation=35, ha="right",
                       fontsize=10, fontweight="bold")
    ax.set_yticks(np.arange(n) + 0.3)
    ax.set_yticklabels(dst_list, rotation=-20, ha="left",
                       fontsize=10, fontweight="bold")
    ax.set_xlabel("Src", fontsize=12, fontweight="bold", labelpad=-2)
    ax.set_ylabel("Dst", fontsize=12, fontweight="bold", labelpad=-2)
    ax.set_zlabel("%",   fontsize=12, fontweight="bold", labelpad=0)
    ax.set_zlim(0, vmax * 1.1 if vmax > 0 else 1)
    ax.tick_params(axis="z", labelsize=9)
    ax.tick_params(axis="x", pad=-7)
    ax.tick_params(axis="y", pad=-7)
    ax.view_init(elev=28, azim=-55)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("input", nargs="+", metavar="CSV",
                    help="One or more *_directions.csv (or *_nm_directions.csv) files")
    ap.add_argument("--benchmark", nargs="+", metavar="BENCH",
                    help="Benchmark specs (default: all found in CSV)")
    ap.add_argument("--structure", nargs="+", choices=list(STRUCTURES),
                    metavar="STRUCT",
                    help="Structures to plot (default: all found in CSV)")
    ap.add_argument("--bandwidth", nargs="+", choices=BANDWIDTHS,
                    metavar="BW",
                    help="Bandwidths to plot (default: all found in CSV)")
    ap.add_argument("--routing", nargs="+", choices=ROUTING_CHOICES,
                    metavar="ROUTING", default=["baseline"])
    ap.add_argument("--near-min-k", nargs="+", type=int, metavar="K",
                    default=[2], help="near_min_k values (default: 2)")
    ap.add_argument("--near-min-p", nargs="+", type=float, metavar="P",
                    default=[1.0], help="near_min_penalty p values (default: 1.0)")
    ap.add_argument("--all-configs", action="store_true",
                    help="Use all structures and bandwidths from experiments.csv")
    ap.add_argument("--metric", choices=["both", "util", "sat"], default="both")
    ap.add_argument("--near", action="store_true",
                    help="Plot near-min direction breakdown (_nm_directions.csv). "
                         "Auto-detected if input filename contains _nm_directions.")
    ap.add_argument("--output", "-o", default=None,
                    help="Save figure to file. With --metric=both, "
                         "_util/_sat suffixes are appended.")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args()

    # ── Auto-detect near mode from filename ───────────────────────────────
    if not args.near:
        args.near = any("_nm_directions" in p for p in args.input)

    # ── Load all provided CSV files ───────────────────────────────────────
    merged_dirs: dict = defaultdict(list)
    for csv_path in args.input:
        if not os.path.isfile(csv_path):
            print(f"[WARN] Not found: {csv_path}")
            continue
        loaded = load_nm_directions(csv_path) if args.near else load_directions(csv_path)
        for k, entries in loaded.items():
            merged_dirs[k].extend(entries)

    all_directions = dict(merged_dirs)
    if not all_directions:
        print("No data in provided CSV files.")
        return

    # ── Auto-find matching summary CSVs ───────────────────────────────────
    all_summary: dict = {}
    for csv_path in args.input:
        for suffix in ("_nm_directions.csv", "_directions.csv"):
            if csv_path.endswith(suffix):
                sum_path = csv_path[: -len(suffix)] + "_summary.csv"
                if os.path.isfile(sum_path):
                    all_summary.update(load_summary(sum_path))
                break

    # ── Determine filter sets ─────────────────────────────────────────────
    # Expand routing keys to match exactly how they are formatted in CSV
    routing_keys: list = []
    for r in args.routing:
        if r == "near_min_adaptive":
            for k in args.near_min_k:
                for p in args.near_min_p:
                    routing_keys.append(routing_to_key(r, near_min_k=k, near_min_p=p))
        elif r == "near_min_random":
            for k in args.near_min_k:
                routing_keys.append(routing_to_key(r, near_min_k=k))
        else:
            routing_keys.append(r)
            
    routing_keys = list(dict.fromkeys(routing_keys)) # Remove duplicates

    # Infer available values from CSV
    csv_benchmarks = sorted({k[0] for k in all_directions})
    csv_structs    = list(dict.fromkeys(k[1] for k in all_directions))  # preserve CSV order
    csv_bws        = list(dict.fromkeys(k[2] for k in all_directions))
    csv_routings   = list(dict.fromkeys(k[3] for k in all_directions))

    if args.all_configs:
        filter_structs = csv_structs
        filter_bws     = csv_bws
    else:
        filter_structs = args.structure if args.structure else csv_structs
        filter_bws     = args.bandwidth if args.bandwidth else csv_bws

    filter_benchmarks = args.benchmark if args.benchmark else csv_benchmarks
    filter_routings   = routing_keys   if routing_keys   else csv_routings

    # ── Verify at least some data matches filters ─────────────────────────
    visible_entries = [
        v for k, v in all_directions.items()
        if k[0] in filter_benchmarks
        and k[1] in filter_structs
        and k[2] in filter_bws
        and k[3] in filter_routings
    ]
    if not visible_entries:
        print("No matching data after applying filters.")
        return

    # ── Plot loop: one figure per structure ───────────────────────────────
    for struct in filter_structs:
        # Rows = bw × benchmark (bw outer, benchmark inner — preserve CSV order)
        row_keys = []
        for bw in filter_bws:
            for bench in filter_benchmarks:
                # Check if any routing has data for this (bench, struct, bw)
                if any((bench, struct, bw, rk) in all_directions for rk in filter_routings):
                    row_keys.append((bw, bench))

        # Cols = routing
        col_keys = filter_routings  # already in user-specified order

        if not row_keys:
            continue

        # Router list scoped to this structure only
        struct_entries = [
            v for k, v in all_directions.items()
            if k[1] == struct
            and k[0] in filter_benchmarks
            and k[2] in filter_bws
            and k[3] in filter_routings
        ]
        router_list = build_router_list(struct_entries)

        n_rows = len(row_keys)
        n_cols = len(col_keys)

        # Determine metrics to plot
        metric_list = []
        if args.metric in ("both", "util"):
            metric_list.append(("pct",     "Link Usage (%)"))
        if args.metric in ("both", "sat"):
            metric_list.append(("avg_sat", "Avg Saturation (%)"))

        for value_key, metric_label in metric_list:
            if args.near:
                metric_label = metric_label.replace(
                    "Link Usage", "Near-Min Direction Usage").replace(
                    "Avg Saturation", "Near-Min Direction Avg Sat")

            # Global z-range for this (struct, metric)
            vmax = 1.0
            for (bw, bench), rk in (
                    (rk_tuple, rk)
                    for rk_tuple in row_keys
                    for rk in col_keys):
                entries = all_directions.get((bench, struct, bw, rk), [])
                for e in entries:
                    if e[value_key] > vmax:
                        vmax = e[value_key]

            cmap = matplotlib.colormaps["coolwarm"]
            norm = colors.Normalize(vmin=0, vmax=vmax)

            cell_w = max(4.5, 3.5 + 0.25 * len(router_list))
            cell_h = max(4.5, 3.5 + 0.25 * len(router_list))
            fig_w  = cell_w * n_cols + 1.8
            fig_h  = cell_h * n_rows + 1.5

            fig = plt.figure(figsize=(fig_w, fig_h))
            fig.suptitle(
                f"{metric_label}  —  {struct}",
                fontsize=26, fontweight="bold", y=0.998)

            for ri, (bw, bench) in enumerate(row_keys):
                for ci, rk in enumerate(col_keys):
                    idx = ri * n_cols + ci + 1
                    ax  = fig.add_subplot(n_rows, n_cols, idx, projection="3d")

                    entries = all_directions.get((bench, struct, bw, rk))
                    if not entries:
                        ax.set_visible(False)
                        continue

                    plot_3d_bars(ax, entries, router_list,
                                 value_key, vmax, cmap, norm)

                    # Cell title: col header on row 0, row label always
                    row_label = f"{bw}\n{short_bench(bench)}"
                    col_label = short_routing(rk)
                    title = f"{row_label}\n[{col_label}]" if ri > 0 else f"{col_label}\n{row_label}"
                    ax.set_title(title, fontsize=18, fontweight="bold", pad=10)

                    # Annotation box
                    ann_mode  = "util" if value_key == "pct" else "sat"
                    sum_row   = all_summary.get((bench, struct, bw, rk))
                    ann_text  = make_annotation(sum_row, ann_mode)
                    if ann_text:
                        ax.text2D(0.5, 0.94, ann_text,
                                  transform=ax.transAxes,
                                  fontsize=14, ha="center", va="top",
                                  family="monospace",
                                  bbox=dict(boxstyle="round,pad=0.4",
                                            facecolor="lightyellow",
                                            edgecolor="gray", alpha=0.9))

            # Colorbar
            sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            cbar_ax = fig.add_axes([0.95, 0.15, 0.015, 0.7])
            cbar    = fig.colorbar(sm, cax=cbar_ax)
            cbar.set_label(metric_label, fontsize=10)
            cbar.ax.tick_params(labelsize=8)

            plt.subplots_adjust(left=0.02, right=0.93, top=0.93,
                                bottom=0.02, wspace=0.08, hspace=0.30)

            # Save or show
            if args.output:
                base, ext = os.path.splitext(args.output)
                if not ext:
                    ext = ".png"
                suffix = ""
                if args.metric == "both":
                    suffix += "_util" if value_key == "pct" else "_sat"
                if len(filter_structs) > 1:
                    suffix += f"_{struct}"
                outpath = f"{base}{suffix}{ext}"
                fig.savefig(outpath, dpi=args.dpi, bbox_inches="tight")
                print(f"Saved: {outpath}")
                plt.close(fig)
            else:
                plt.show()
                plt.close(fig)


if __name__ == "__main__":
    main()
