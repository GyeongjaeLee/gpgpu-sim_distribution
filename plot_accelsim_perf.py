#!/usr/bin/env python3
"""
Compare normalized performance and L2 BW across benchmarks.
Normalization is ALWAYS performed relative to the 'baseline' routing 
within each specific benchmark and struct+bw combination.
(Benchmark names are now derived from the input CSV filenames.)
"""

import argparse
import csv
import os
import re
from typing import Optional

import matplotlib
import matplotlib.pyplot as plt

from experiments_loader import (load_experiments,
                                 STRUCT_ABBREV, BW_ABBREV,
                                 bench_abbrev)

# ── Constants ──────────────────────────────────────────────────────────────────

_HERE = os.path.dirname(os.path.abspath(__file__))
STRUCTURES, BANDWIDTHS = load_experiments(os.path.join(_HERE, "experiments.csv"))

ROUTING_CHOICES = [
    "baseline", "min_oblivious", "min_adaptive",
    "near_min_adaptive", "near_min_random", "fixed_min", 
    "ugal", "valiant",
]

CLOCK_FREQ_MHZ: dict = {
    "H100":             1980,
    "B100_Local":       1965,
    "B100_Global":      1965,
    "B100_Core_Rotate": 1965,
    "Rubin_Ultra":      2380,
}

Y_LIM_NORM = (0.8, 1.8)

# ── Helpers ────────────────────────────────────────────────────────────────────

def routing_to_key(routing: str, near_min_k: Optional[int] = None, near_min_p: Optional[float] = None) -> str:
    if routing == "near_min_adaptive":
        k = near_min_k if near_min_k is not None else 2
        p = near_min_p if near_min_p is not None else 1.0
        return f"{routing}_nmk{k}_nmp{p:.1f}"
    elif routing == "near_min_random":
        k = near_min_k if near_min_k is not None else 2
        return f"{routing}_nmk{k}"
    return routing

def sb_label(struct: str, bw: str) -> str:
    return f"{STRUCT_ABBREV.get(struct, struct[:3])}/{BW_ABBREV.get(bw, bw[:3])}"

def clock_mhz(struct: str) -> float:
    return float(CLOCK_FREQ_MHZ.get(struct, 1965))

def pct_annotation(val: float) -> str:
    pct = (val - 1.0) * 100.0
    if abs(pct) < 0.05: return "ref"
    return f"{pct:+.1f}%"

def rk_full_name(rk: str) -> str:
    names = {
        "baseline":      "Baseline",
        "min_oblivious": "Min Oblivious",
        "min_adaptive":  "Min Adaptive",
        "fixed_min":     "Fixed Min (Fabric Pref)",
        "ugal":          "UGAL",
        "valiant":       "Valiant",
    }
    if rk in names: return names[rk]
    m_adp = re.match(r"^near_min_adaptive_nmk(\d+)_nmp([\d.]+)$", rk)
    if m_adp: return f"Near-Min Adp (k={m_adp.group(1)}, p={m_adp.group(2)})"
    m_rnd = re.match(r"^near_min_random_nmk(\d+)$", rk)
    if m_rnd: return f"Near-Min Rnd (k={m_rnd.group(1)})"
    return rk

# ── Data loading ───────────────────────────────────────────────────────────────

def load_csvs(paths: list, f_structs: Optional[set], f_bws: Optional[set], f_routings: Optional[set]) -> list:
    rows = []
    for path in paths:
        if not os.path.isfile(path): continue
        
        file_bench_name = os.path.splitext(os.path.basename(path))[0]
        
        with open(path, newline='') as f:
            for row in csv.DictReader(f):
                s, b, r = row.get('struct', ''), row.get('bw', ''), row.get('routing', '')
                if f_structs and s not in f_structs: continue
                if f_bws and b not in f_bws: continue
                if f_routings and r not in f_routings: continue
                try:
                    rows.append({
                        'benchmark': file_bench_name, 
                        'struct': s, 'bw': b, 'routing': r,
                        'tot_cycle': int(row['tot_cycle']), 'l2_bw': float(row['l2_bw_total']),
                    })
                except (ValueError, KeyError): continue
    return rows

# ── Normalization & Data Build ─────────────────────────────────────────────────

def build_plot_data(rows: list, benches: list, structs_bws: list, active_rks: list, no_normalize: bool) -> dict:
    raw: dict = {}
    for r in rows:
        key = (r['benchmark'], (r['struct'], r['bw']), r['routing'])
        lat = r['tot_cycle']  # 🔥 변경됨: 시간 변환 제거, 원본 사이클 그대로 사용
        raw[key] = (lat, r['l2_bw'])

    data_dict: dict = {}
    for bench in benches:
        for sb in structs_bws:
            if not no_normalize:
                ref_key = (bench, sb, 'baseline')
                if ref_key in raw:
                    ref_lat, ref_bw = raw[ref_key]
                else:
                    group_vals = [v for k, v in raw.items() if k[0] == bench and k[1] == sb]
                    if not group_vals: continue
                    ref_lat = max(v[0] for v in group_vals)
                    ref_bw  = min((v[1] for v in group_vals if v[1] > 0), default=1.0)
                    print(f"  [INFO] 'baseline' not found for {bench}/{sb_label(*sb)}. Using group worst as ref.")

            for rk in active_rks:
                key = (bench, sb, rk)
                if key not in raw: continue
                lat, bw = raw[key]
                
                if no_normalize:
                    # 그대로 반환 (사이클, BW)
                    data_dict[key] = (lat, bw)
                else:
                    # 정규화 반환 (Speedup, Normalized BW)
                    # 동일한 클럭 기준이므로 (ref_cycle / cycle)은 그대로 정확한 Speedup이 됩니다.
                    data_dict[key] = (
                        ref_lat / max(lat, 1e-30),
                        bw / max(ref_bw, 1e-30)
                    )
    return data_dict

# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input', nargs='+')
    ap.add_argument('--structure', nargs='+', choices=list(STRUCTURES))
    ap.add_argument('--bandwidth', nargs='+', choices=list(BANDWIDTHS))
    ap.add_argument('--routing',   nargs='+', choices=ROUTING_CHOICES)
    ap.add_argument('--near-min-k', nargs='+', type=int, default=[2])
    ap.add_argument('--near-min-p', nargs='+', type=float, default=[1.0])
    ap.add_argument('--benchmark', nargs='+')
    ap.add_argument('--no-normalize', action='store_true', help='Plot raw values instead of normalizing to baseline')
    ap.add_argument('--output', '-o', default=None)
    ap.add_argument('--dpi',     type=int,   default=150)
    args = ap.parse_args()

    rk_list = []
    if args.routing:
        for r in args.routing:
            if r == 'near_min_adaptive':
                for k in args.near_min_k:
                    for p in args.near_min_p: rk_list.append(routing_to_key(r, k, p))
            elif r == 'near_min_random':
                for k in args.near_min_k: rk_list.append(routing_to_key(r, k))
            else: rk_list.append(r)
    
    load_rk_filter = set(rk_list) | {'baseline'} if rk_list else None

    rows = load_csvs(args.input, set(args.structure) if args.structure else None,
                     set(args.bandwidth) if args.bandwidth else None, load_rk_filter)
    
    if not rows: return

    seen_benches = list(dict.fromkeys(r['benchmark'] for r in rows))
    benches = [b for b in (args.benchmark or []) if b in seen_benches]
    benches += [b for b in seen_benches if b not in benches]

    structs_bws = list(dict.fromkeys((r['struct'], r['bw']) for r in rows))
    csv_rks = list(dict.fromkeys(r['routing'] for r in rows))
    active_rks = [rk for rk in (rk_list or csv_rks) if rk in csv_rks]

    plot_data = build_plot_data(rows, benches, structs_bws, active_rks, args.no_normalize)

    # ── Layout ────────────────────────────────────────────────────────────
    n_r, n_sb, n_b = len(active_rks), len(structs_bws), len(benches)
    bar_w, sb_gap, bch_gap = 0.22, 0.25, 1.0
    sb_w = n_r * bar_w
    bch_w = n_sb * sb_w + (n_sb - 1) * sb_gap

    bar_xs, sb_cxs, bch_cxs, sep_xs = {}, {}, [], []
    x0 = 0.0
    for bi in range(n_b):
        bch_cxs.append(x0 + bch_w / 2)
        for si in range(n_sb):
            sb_x0 = x0 + si * (sb_w + sb_gap)
            sb_cxs[(bi, si)] = sb_x0 + sb_w / 2
            for ri in range(n_r): bar_xs[(bi, si, ri)] = sb_x0 + ri * bar_w + bar_w / 2
        if bi < n_b - 1: sep_xs.append(x0 + bch_w + bch_gap / 2)
        x0 += bch_w + bch_gap

    rk_color = {rk: plt.cm.tab10(i % 10) for i, rk in enumerate(active_rks)}
    figw = max(14.0, n_b * (bch_w + bch_gap) * 0.95 + 3.0)
    out_base = args.output or os.path.splitext(os.path.basename(args.input[0]))[0]
    
    out_prefix = f"{out_base}_raw" if args.no_normalize else out_base

    # 🔥 Plot 라벨 수정: Execution Time -> Total Cycles
    plot_specs = [
        (0, "Cycles" if args.no_normalize else "Performance", 
            "Total Cycles" if args.no_normalize else "Speedup (vs. Baseline)", 
            f"{out_prefix}_perf.png"),
        (1, "L2 BW", 
            "L2 Total BW (Raw)" if args.no_normalize else "Normalized L2 BW (vs. Baseline)", 
            f"{out_prefix}_bw.png")
    ]

    for vi, base_title, ylabel, out_path in plot_specs:
        fig, ax = plt.subplots(1, 1, figsize=(figw, 9))
        legend_added = set()
        
        # Y축 자동 스케일링 (정규화가 꺼졌을 경우)
        if args.no_normalize:
            all_vals = [v[vi] for v in plot_data.values()]
            max_val = max(all_vals) if all_vals else 1.0
            current_ylim = (0, max_val * 1.15)
            title = f"Raw {base_title} " + ("(Lower is Better)" if vi == 0 else "(Higher is Better)")
        else:
            current_ylim = Y_LIM_NORM
            title = f"Normalized {base_title}"

        for bi, bench in enumerate(benches):
            for si, sb in enumerate(structs_bws):
                for ri, rk in enumerate(active_rks):
                    val = plot_data.get((bench, sb, rk))
                    if val is None: continue
                    
                    x, h = bar_xs[(bi, si, ri)], val[vi]
                    ax.bar(x, h, width=bar_w, color=rk_color[rk], edgecolor='white', linewidth=0.4, 
                           label=rk_full_name(rk) if rk not in legend_added else '_')
                    legend_added.add(rk)
                    
                    # 🔥 어노테이션 포맷 변경 (사이클은 정수형으로 출력)
                    if args.no_normalize:
                        if vi == 0:
                            ann_text = f"{int(h)}"  # 사이클은 소수점 없이 정수로
                        else:
                            ann_text = f"{h:.2f}"
                        ann_y = h + (current_ylim[1] * 0.01)
                    else:
                        ann_text = pct_annotation(h)
                        ann_y = min(h, current_ylim[1]) + 0.01

                    ax.text(x, ann_y, ann_text, ha='center', va='bottom', 
                            fontsize=10, fontweight='bold', rotation=90, clip_on=False)

        if not args.no_normalize:
            ax.axhline(1.0, color='#333', linestyle='--', linewidth=1.2, alpha=0.7)
            
        for sx in sep_xs: ax.axvline(sx, color='gray', linestyle=':', linewidth=1.0, alpha=0.4)
        
        for bi, bench in enumerate(benches):
            ax.text(bch_cxs[bi], -0.18 * current_ylim[1] if args.no_normalize else -0.18, 
                    bench_abbrev(bench), ha='center', va='top', 
                    fontsize=16, fontweight='bold', transform=ax.get_xaxis_transform())

        ax.set_xticks([])
        ax.set_xlim(-bch_gap/2, x0 - bch_gap + bch_gap/2)
        ax.set_ylim(*current_ylim)
        ax.set_ylabel(ylabel, fontsize=18, fontweight='bold')
        ax.set_title(title, fontsize=22, fontweight='bold', pad=20)
        ax.legend(title='Routing', loc='upper right', framealpha=0.9)
        ax.grid(axis='y', linestyle=':', alpha=0.4)
        plt.subplots_adjust(bottom=0.25)
        fig.savefig(out_path, dpi=args.dpi, bbox_inches='tight')
        plt.close(fig)

if __name__ == '__main__':
    main()