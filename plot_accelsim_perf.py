#!/usr/bin/env python3
"""
Compare normalized performance and L2 BW across benchmarks from AccelSim summary CSVs.

Reads one or more *_summary.csv files produced by parse_accelsim_results.py.
Saves two PNG files automatically: {prefix}_perf.png and {prefix}_bw.png.

Reference combination for normalization (worst expected):
  - struct   : first in experiments.csv order among selected (top = worst)
  - bandwidth: first in experiments.csv order among selected (top = worst)
  - routing  : 'baseline' if selected, else first selected routing

Normalized values (both higher = better, reference = 1.0):
  - Performance = ref_lat / this_lat   (speedup; ref_lat = ref cycles / clock_freq)
  - L2 BW       = this_bw / ref_bw

Clock frequencies:
  H100            → 1980 MHz
  B100_* (B200)   → 1965 MHz
  Rubin_Ultra     → 2380 MHz

X-axis layout:
  Benchmark → (struct+bw group) → routing bars

Usage examples:

  # Compare two benchmarks across structures and routings
  python plot_accelsim_perf.py bfs_BG_H3_bas.csv gemm_BG_H3_bas.csv \\
      --structure B100_Global H100 \\
      --bandwidth B200+HBM3e \\
      --routing baseline min_adaptive near_min_adaptive \\
      --near-min-p 0.0 1.0 \\
      -o comparison

  # Single benchmark, all configs, baseline vs near-min
  python plot_accelsim_perf.py mem_bw_summary.csv \\
      --structure B100_Local B100_Global B100_Core_Rotate \\
      --bandwidth B200+HBM3e Shoreline_1x Shoreline_2x \\
      --routing baseline near_min_adaptive \\
      --near-min-p 0.0 \\
      -o mem_bw_sweep

  # Omit --structure/--bandwidth to include all rows found in the CSV
  python plot_accelsim_perf.py bfs_BG_H3_bas.csv \\
      --routing baseline near_min_adaptive \\
      --near-min-p 0.0
"""

import argparse
import csv
import os
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt

import re

from experiments_loader import (load_experiments,
                                 STRUCT_ABBREV, BW_ABBREV,
                                 bench_abbrev)

# ── Constants ──────────────────────────────────────────────────────────────────

_HERE = os.path.dirname(os.path.abspath(__file__))
STRUCTURES, BANDWIDTHS = load_experiments(os.path.join(_HERE, "experiments.csv"))

ROUTING_CHOICES = [
    "baseline", "min_oblivious", "min_adaptive",
    "near_min_adaptive", "ugal", "valiant",
]

CLOCK_FREQ_MHZ: dict = {
    "H100":             1980,
    "B100_Local":       1965,
    "B100_Global":      1965,
    "B100_Core_Rotate": 1965,
    "Rubin_Ultra":      2380,
}

Y_LIM = (0.8, 1.8)


# ── Helpers ────────────────────────────────────────────────────────────────────

def routing_to_key(routing: str, near_min_p: Optional[float] = None) -> str:
    if routing == "near_min_adaptive":
        p = near_min_p if near_min_p is not None else 1.0
        return f"near_min_p{p:.1f}"
    return routing


def sb_label(struct: str, bw: str) -> str:
    return f"{STRUCT_ABBREV.get(struct, struct[:3])}/{BW_ABBREV.get(bw, bw[:3])}"


def clock_mhz(struct: str) -> float:
    return float(CLOCK_FREQ_MHZ.get(struct, 1965))


def pct_annotation(val: float) -> str:
    """e.g. 1.0 → 'ref', 1.15 → '+15.0%', 0.92 → '-8.0%'"""
    pct = (val - 1.0) * 100.0
    if abs(pct) < 0.05:
        return "ref"
    return f"{pct:+.1f}%"


def rk_full_name(rk: str) -> str:
    """Full human-readable routing name for legend."""
    names = {
        "baseline":      "Baseline",
        "min_oblivious": "Min Oblivious",
        "min_adaptive":  "Min Adaptive",
        "ugal":          "UGAL",
        "valiant":       "Valiant",
    }
    if rk in names:
        return names[rk]
    m = re.match(r"^near_min_p([\d.]+)$", rk)
    if m:
        return f"Near-Min Adaptive (p={m.group(1)})"
    return rk


# ── Data loading ───────────────────────────────────────────────────────────────

def load_csvs(paths: list,
              filter_structs:  Optional[set],
              filter_bws:      Optional[set],
              filter_routings: Optional[set]) -> list:
    rows = []
    for path in paths:
        if not os.path.isfile(path):
            print(f"[WARN] Not found: {path}")
            continue
        with open(path, newline='') as f:
            for row in csv.DictReader(f):
                struct  = row.get('struct', '')
                bw      = row.get('bw', '')
                routing = row.get('routing', '')
                if filter_structs  and struct  not in filter_structs:  continue
                if filter_bws      and bw      not in filter_bws:      continue
                if filter_routings and routing not in filter_routings:  continue
                try:
                    tot_cycle = int(row['tot_cycle'])
                    l2_bw     = float(row['l2_bw_total'])
                except (ValueError, KeyError):
                    continue
                rows.append({
                    'benchmark': row['benchmark'],
                    'struct':    struct,
                    'bw':        bw,
                    'routing':   routing,
                    'tot_cycle': tot_cycle,
                    'l2_bw':     l2_bw,
                })
    return rows


# ── Normalization ──────────────────────────────────────────────────────────────

def find_reference(structs_bws: list, active_rks: list) -> tuple:
    """
    Reference = (worst struct+bw, worst routing).
    Worst struct+bw: first in experiments.csv order (top row = worst performance).
    Worst routing  : 'baseline' if present, else first in active_rks.
    """
    struct_rank = {s: i for i, s in enumerate(STRUCTURES)}
    bw_rank     = {b: i for i, b in enumerate(BANDWIDTHS)}
    sorted_sb   = sorted(structs_bws,
                         key=lambda sb: (struct_rank.get(sb[0], 999),
                                         bw_rank.get(sb[1], 999)))
    worst_sb = sorted_sb[0]
    worst_rk = 'baseline' if 'baseline' in active_rks else active_rks[0]
    return worst_sb, worst_rk


def build_normalized(rows: list, benches: list,
                     structs_bws: list, active_rks: list) -> dict:
    """
    Returns norm[(bench, sb, rk)] = (speedup, bw_ratio).
      speedup  = ref_lat / this_lat   (higher is better, ref = 1.0)
      bw_ratio = this_bw / ref_bw     (higher is better, ref = 1.0)
    """
    # Raw lookup
    raw: dict = {}
    for r in rows:
        key = (r['benchmark'], (r['struct'], r['bw']), r['routing'])
        lat = r['tot_cycle'] / (clock_mhz(r['struct']) * 1e6)
        raw[key] = (lat, r['l2_bw'])

    worst_sb, worst_rk = find_reference(structs_bws, active_rks)
    print(f"Reference: struct+bw={sb_label(*worst_sb)}  routing={worst_rk}")

    norm: dict = {}
    for bench in benches:
        ref_key = (bench, worst_sb, worst_rk)
        if ref_key in raw:
            ref_lat, ref_bw = raw[ref_key]
        else:
            # Fallback: use worst observed values for this benchmark
            bench_vals = [(lat, bw) for (b, _, __), (lat, bw) in raw.items() if b == bench]
            if not bench_vals:
                ref_lat, ref_bw = 1.0, 1.0
            else:
                ref_lat = max(v[0] for v in bench_vals)
                ref_bw  = min((v[1] for v in bench_vals if v[1] > 0), default=1.0)
            print(f"  [WARN] Reference ({sb_label(*worst_sb)}, {worst_rk}) not found "
                  f"for {bench_abbrev(bench)} — using worst observed.")

        for sb in structs_bws:
            for rk in active_rks:
                key = (bench, sb, rk)
                if key not in raw:
                    continue
                lat, bw = raw[key]
                norm[key] = (
                    ref_lat / max(lat, 1e-30),
                    bw      / max(ref_bw, 1e-30),
                )
    return norm


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('input', nargs='+', metavar='CSV',
                    help='*_summary.csv files from parse_accelsim_results.py')
    ap.add_argument('--structure', nargs='+', choices=list(STRUCTURES), metavar='STRUCT')
    ap.add_argument('--bandwidth', nargs='+', choices=list(BANDWIDTHS), metavar='BW')
    ap.add_argument('--routing',   nargs='+', choices=ROUTING_CHOICES,  metavar='ROUTING')
    ap.add_argument('--near-min-p', nargs='+', type=float, default=None, metavar='P',
                    help='near_min_adaptive p values (default: all found in CSV)')
    ap.add_argument('--bench-order', nargs='+', metavar='BENCH',
                    help='Explicit benchmark display order (full spec)')
    ap.add_argument('--output', '-o', default=None)
    ap.add_argument('--dpi',     type=int,   default=150)
    ap.add_argument('--figsize', nargs=2, type=float, metavar=('W', 'H'))
    args = ap.parse_args()

    # First pass: build rk_list for non-near-min routings
    rk_list: list = []
    near_min_requested = False
    if args.routing:
        for r in args.routing:
            if r == 'near_min_adaptive':
                near_min_requested = True
                if args.near_min_p:
                    for p in args.near_min_p:
                        rk_list.append(routing_to_key(r, p))
            else:
                rk_list.append(r)

    # Load (without near_min filter if p not specified — resolve from CSV)
    routing_filter = set(rk_list) if (rk_list and not near_min_requested) else None
    rows = load_csvs(
        args.input,
        filter_structs  = set(args.structure) if args.structure else None,
        filter_bws      = set(args.bandwidth) if args.bandwidth else None,
        filter_routings = routing_filter,
    )
    if not rows:
        print("No data after filtering.")
        return

    # Add all near_min_p* found in CSV if requested without explicit p values
    if near_min_requested and not args.near_min_p:
        csv_nm = list(dict.fromkeys(
            r['routing'] for r in rows if r['routing'].startswith('near_min_p')))
        rk_list.extend(csv_nm)

    # Post-filter rows to match rk_list (handles the near_min case above)
    if rk_list:
        keep = set(rk_list)
        rows = [r for r in rows if r['routing'] in keep]
    if not rows:
        print("No data after filtering.")
        return

    # Ordered lists
    seen_benches = list(dict.fromkeys(r['benchmark'] for r in rows))
    if args.bench_order:
        benches  = [b for b in args.bench_order if b in seen_benches]
        benches += [b for b in seen_benches if b not in benches]
    else:
        benches = sorted(seen_benches)

    structs_bws = list(dict.fromkeys((r['struct'], r['bw']) for r in rows))

    csv_rks    = list(dict.fromkeys(r['routing'] for r in rows))
    active_rks = [rk for rk in (rk_list or csv_rks) if rk in csv_rks]
    if not active_rks:
        print("No matching routing data.")
        return

    # Normalize
    norm = build_normalized(rows, benches, structs_bws, active_rks)

    # ── Layout ────────────────────────────────────────────────────────────
    n_r  = len(active_rks)
    n_sb = len(structs_bws)
    n_b  = len(benches)

    bar_w   = 0.25
    sb_gap  = 0.30
    bch_gap = 1.0

    sb_w  = n_r * bar_w
    bch_w = n_sb * sb_w + (n_sb - 1) * sb_gap

    bar_xs  = {}
    sb_cxs  = {}
    bch_cxs = []
    sep_xs  = []

    x0 = 0.0
    for bi in range(n_b):
        bch_cxs.append(x0 + bch_w / 2)
        for si in range(n_sb):
            sb_x0 = x0 + si * (sb_w + sb_gap)
            sb_cxs[(bi, si)] = sb_x0 + sb_w / 2
            for ri in range(n_r):
                bar_xs[(bi, si, ri)] = sb_x0 + ri * bar_w + bar_w / 2
        if bi < n_b - 1:
            sep_xs.append(x0 + bch_w + bch_gap / 2)
        x0 += bch_w + bch_gap

    x_max = x0 - bch_gap

    cmap      = matplotlib.colormaps['tab10']
    rk_color  = {rk: cmap(i / max(n_r, 1)) for i, rk in enumerate(active_rks)}

    # ── Figure setup ──────────────────────────────────────────────────────
    # Width scales with number of bars; height is fixed
    auto_w = max(14.0, n_b * (bch_w + bch_gap) * 0.95 + 3.0)
    figw   = args.figsize[0] if args.figsize else auto_w
    figh   = args.figsize[1] if args.figsize else 10.0

    # Auto output prefix: first input CSV stem
    if args.output:
        out_base = os.path.splitext(args.output)[0]
    else:
        out_base = os.path.splitext(os.path.basename(args.input[0]))[0]

    plot_specs = [
        (0, "Performance",  "Normalized Performance\n(higher = better)",
         f"{out_base}_perf.png"),
        (1, "L2 BW",        "Normalized L2 Total BW\n(higher = better)",
         f"{out_base}_bw.png"),
    ]

    for vi, fig_title, ylabel, out_path in plot_specs:
        fig, ax = plt.subplots(1, 1, figsize=(figw, figh))
        legend_added: set = set()

        for bi, bench in enumerate(benches):
            for si, sb in enumerate(structs_bws):
                for ri, rk in enumerate(active_rks):
                    val = norm.get((bench, sb, rk))
                    if val is None:
                        continue
                    x = bar_xs[(bi, si, ri)]
                    h = val[vi]
                    lbl = rk_full_name(rk) if rk not in legend_added else '_'
                    ax.bar(x, h, width=bar_w, color=rk_color[rk],
                           edgecolor='white', linewidth=0.4, label=lbl)
                    legend_added.add(rk)

                    # Annotation above bar
                    ann_y = min(h, Y_LIM[1]) + 0.01
                    ax.text(x, ann_y, pct_annotation(h),
                            ha='center', va='bottom',
                            fontsize=13, fontweight='bold',
                            rotation=90, clip_on=False)

        # Reference line
        ax.axhline(1.0, color='#333', linestyle='--', linewidth=1.2,
                   alpha=0.7, zorder=0)

        # Benchmark separators
        for sx in sep_xs:
            ax.axvline(sx, color='gray', linestyle=':', linewidth=1.0, alpha=0.5)

        # Struct+bw sub-labels (first level)
        for bi in range(n_b):
            for si in range(n_sb):
                ax.text(sb_cxs[(bi, si)], -0.04,
                        sb_label(*structs_bws[si]),
                        ha='center', va='top',
                        fontsize=15, fontweight='bold', rotation=30,
                        transform=ax.get_xaxis_transform(), clip_on=False)

        # Benchmark labels (second level)
        for bi, bench in enumerate(benches):
            ax.text(bch_cxs[bi], -0.16,
                    bench_abbrev(bench),
                    ha='center', va='top',
                    fontsize=20, fontweight='bold',
                    transform=ax.get_xaxis_transform(), clip_on=False)

        ax.set_xticks([])
        ax.set_xlim(-bch_gap / 2, x_max + bch_gap / 2)
        ax.set_ylim(*Y_LIM)
        ax.set_ylabel(ylabel, fontsize=20, fontweight='bold', labelpad=10)
        ax.set_title(fig_title, fontsize=24, fontweight='bold', pad=12)
        ax.tick_params(axis='y', labelsize=13, width=1.5)
        for tick in ax.get_yticklabels():
            tick.set_fontweight('bold')
        ax.legend(title='Routing', fontsize=10, title_fontsize=11,
                  loc='upper right', framealpha=0.85)
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        ax.spines[['top', 'right']].set_visible(False)
        ax.spines[['left', 'bottom']].set_linewidth(1.5)

        plt.subplots_adjust(bottom=0.22)
        fig.savefig(out_path, dpi=args.dpi, bbox_inches='tight')
        print(f"Saved: {out_path}")
        plt.close(fig)


if __name__ == '__main__':
    main()
