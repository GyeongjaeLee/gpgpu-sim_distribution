#!/usr/bin/env python3
"""
Parse AccelSim simulation logs for HBMNet experiments and write CSV files.

For each (structure × bandwidth × routing) combination:
  - Reads ../../util/job_launching/apps/define-all-apps.yml to find active args
  - Discovers the log file under
      ../../sim_run_12.8/{app}/{sanitized_args}/{config}-SASS/{app}*.o{N}
    where {sanitized_args} matches AccelSim's get_argfoldername() output.
    Falls back to arg* glob if YAML is unavailable. Picks the largest N suffix.
  - Parses per-kernel network and cache stats, then aggregates across all kernels.
  - Global metrics (gpu_tot_sim_cycle, L2_BW_total, latency, inj_rate) are taken 
    from the LAST occurrence in the log (cumulative totals).

Outputs:
  {prefix}_summary.csv       – one row per (struct, bw, routing)
  {prefix}_directions.csv    – one row per (struct, bw, routing, src, dst)
  {prefix}_nm_directions.csv – near-min direction breakdown
"""

import argparse
import csv
import glob
import hashlib
import os
import re
from typing import Optional

from experiments_loader import (load_experiments,
                                bench_abbrev, route_abbrev, make_csv_prefix)

# ── Paths ─────────────────────────────────────────────────────────────────────
_HERE            = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR         = os.path.normpath(os.path.join(_HERE, "..", ".."))
SIM_RUN_BASE     = os.path.join(ROOT_DIR, "sim_run_12.8")
DEFINE_APPS_YAML = os.path.join(ROOT_DIR, "util", "job_launching", "apps", "define-all-apps.yml")

STRUCTURES, BANDWIDTHS = load_experiments(os.path.join(_HERE, "experiments.csv"))

ROUTING_CHOICES = [
    "baseline", "min_oblivious", "min_adaptive",
    "near_min_adaptive", "near_min_random", "fixed_min", 
    "ugal", "valiant",
]

LINK_TYPES     = ["XBAR_XBAR", "XBAR_HBM", "XBAR_MC", "MC_HBM", "MC_MC"]
SAT_LINK_TYPES = ["XBAR_XBAR", "XBAR_MC", "MC_MC"]


# ── Helpers ───────────────────────────────────────────────────────────────────

def routing_to_key(routing: str, near_min_k: Optional[int] = None, near_min_p: Optional[float] = None) -> str:
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

def _parse_active_args_fallback(yaml_path: str, group_name: str, app_name: str) -> list:
    args_list = []
    in_group, in_app, app_col = False, False, None
    with open(yaml_path) as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            stripped = line.lstrip()
            if not stripped or stripped.startswith("#"): continue
            col = len(line) - len(stripped)

            if not in_group:
                if col == 0 and stripped == group_name + ":": in_group = True
                continue
            if col == 0: break

            if not in_app:
                m = re.match(r'^(\s*)-\s+' + re.escape(app_name) + r'\s*:\s*$', line)
                if m:
                    in_app, app_col = True, col
                continue
            if col <= app_col and re.match(r'^\s*-\s+\S', line): break
            m = re.match(r'^\s+-\s+args:\s+(.+)$', line)
            if m: args_list.append(m.group(1).strip())
    return args_list

def get_active_args(benchmark: str, yaml_path: str = DEFINE_APPS_YAML) -> list:
    if not os.path.isfile(yaml_path): return []
    parts = benchmark.split(":", 1)
    if len(parts) != 2: return []
    group_name, app_name = parts

    try:
        import yaml
        with open(yaml_path) as f: content = yaml.safe_load(f)
        group = content.get(group_name)
        if group:
            for exec_entry in group.get("execs", []):
                if not isinstance(exec_entry, dict): continue
                app_data = exec_entry.get(app_name)
                if app_data is None: continue
                return [str(cfg["args"]).strip() for cfg in app_data if isinstance(cfg, dict) and "args" in cfg]
        return []
    except ImportError:
        pass
    return _parse_active_args_fallback(yaml_path, group_name, app_name)

def sanitize_args_for_path(args: str) -> str:
    if not args or not args.strip(): return "NO_ARGS"
    s = args.strip()
    foldername = re.sub(r"[^a-z^A-Z^0-9]", "_", s)
    if len(s) > 256: foldername = "hashed_args_" + hashlib.md5(s.encode()).hexdigest()
    return foldername

def find_log_file(app_name: str, config_name: str, sim_run_base: str, active_args: Optional[list] = None) -> Optional[str]:
    candidates = []
    if active_args:
        for args in active_args:
            args_dir = sanitize_args_for_path(args)
            pattern  = os.path.join(sim_run_base, app_name, args_dir, config_name, f"{app_name}*.o*")
            for path in glob.glob(pattern):
                m = re.search(r'\.o(\d+)$', path)
                if m: candidates.append((int(m.group(1)), path))

    if not candidates:
        pattern = os.path.join(sim_run_base, app_name, "*", config_name, f"{app_name}*.o*")
        for path in glob.glob(pattern):
            m = re.search(r'\.o(\d+)$', path)
            if m: candidates.append((int(m.group(1)), path))

    return max(candidates, key=lambda x: x[0])[1] if candidates else None

# ── Log parser ────────────────────────────────────────────────────────────────

def parse_log(path: str) -> Optional[dict]:
    with open(path) as f: text = f.read()

    # Global metrics
    all_cycles = re.findall(r'^gpu_tot_sim_cycle\s*=\s*(\d+)', text, re.MULTILINE)
    tot_cycle = int(all_cycles[-1]) if all_cycles else None

    all_bw = re.findall(r'^L2_BW_total\s*=\s*([\d.]+)', text, re.MULTILINE)
    l2_bw_total = float(all_bw[-1]) if all_bw else None

    # L2
    total_l2_acc  = sum(int(m.group(1)) for m in re.finditer(r'^L2_total_cache_accesses\s*=\s*(\d+)', text, re.MULTILINE))
    total_l2_miss = sum(int(m.group(1)) for m in re.finditer(r'^L2_total_cache_misses\s*=\s*(\d+)', text, re.MULTILINE))
    l2_hit_rate = (1 - total_l2_miss / total_l2_acc) if total_l2_acc > 0 else None

    # Link stats
    link_starts = [m.start() for m in re.finditer(r'=== Link Type Utilization ===', text)]
    if not link_starts: return None

    num_kernels = len(link_starts)
    last_link_start = link_starts[-1]
    m_end = re.search(r'Miss link avg saturation:.*', text[last_link_start:])
    last_link_end = last_link_start + m_end.end() if m_end else len(text)
    link_block = text[last_link_start:last_link_end]

    link_util: dict = {}
    for lm in re.finditer(r'^\s+(XBAR_XBAR|XBAR_HBM|XBAR_MC|MC_HBM|MC_MC):\s*(\d+)\s*\(([\d.]+)%\)', link_block, re.MULTILINE):
        link_util[lm.group(1)] = {'count': int(lm.group(2)), 'pct': float(lm.group(3))}

    m2 = re.search(r'Total link traversals:\s*(\d+)', link_block)
    total_traversals = int(m2.group(1)) if m2 else 0

    dir_entries = []
    for dm in re.finditer(r'^\s{4}(\S+(?:\([^)]*\))?)\s*->\s*(\S+(?:\([^)]*\))?):\s*count=(\d+)\s*\(([\d.]+)%\)\s*avg_sat=([\d.]+)%', link_block, re.MULTILINE):
        dir_entries.append({
            'src': canonical_router_name(dm.group(1)), 'dst': canonical_router_name(dm.group(2)),
            'count': int(dm.group(3)), 'pct': float(dm.group(4)), 'avg_sat': float(dm.group(5)),
        })

    m2 = re.search(r'Miss link avg saturation:\s*XBAR_XBAR=([\d.]+)%\s*XBAR_MC=([\d.]+)%\s*MC_MC=([\d.]+)%', link_block)
    agg_sat = {
        'XBAR_XBAR': float(m2.group(1)) if m2 else 0.0,
        'XBAR_MC':   float(m2.group(2)) if m2 else 0.0,
        'MC_MC':     float(m2.group(3)) if m2 else 0.0,
    }

    # ── Traffic Class / Latency / Injection Stats (LAST occurrence) ─────────
    pkt_lat_avg, pkt_lat_max = None, None
    inj_rate_avg, inj_rate_max = None, None

    tc0_starts = [m.start() for m in re.finditer(r'====== Traffic class 0 ======', text)]
    if tc0_starts:
        last_tc0 = text[tc0_starts[-1]:]
        
        # Regex captures average, and optionally looks ahead for the maximum value if present
        m_lat = re.search(r'Packet latency average = ([\d.]+)(?:[^\n]*\n[^\n]*minimum[^\n]*\n[^\n]*maximum = ([\d.]+))?', last_tc0)
        if m_lat:
            pkt_lat_avg = float(m_lat.group(1))
            if m_lat.group(2): pkt_lat_max = float(m_lat.group(2))

        m_inj = re.search(r'Injected packet rate average = ([\d.]+)(?:[^\n]*\n[^\n]*minimum[^\n]*\n[^\n]*maximum = ([\d.]+))?', last_tc0)
        if m_inj:
            inj_rate_avg = float(m_inj.group(1))
            if m_inj.group(2): inj_rate_max = float(m_inj.group(2))

    # Near-min
    nm_starts = [m.start() for m in re.finditer(r'=== Near-Min Adaptive Decisions ===', text)]
    near_min: Optional[dict] = None
    nm_dir_entries: list = []
    if nm_starts:
        last_nm_start = nm_starts[-1]
        hdr_end = text.index('\n', last_nm_start) + 1
        next_sec = re.search(r'^===', text[hdr_end:], re.MULTILINE)
        nm_end   = hdr_end + next_sec.start() if next_sec else len(text)
        nm_block = text[last_nm_start:nm_end]

        nm_m = re.search(r'Near-min decisions:\s*total=(\d+)\s+min=(\d+)\s+non-min=(\d+)\s+near-min ratio=([\d.]+)%', nm_block)
        nm_path = re.search(r'Near-min path usage:\s*(\d+)\s*/\s*(\d+)\s*\(([\d.]+)%\)', nm_block)
        if nm_m:
            non_min_total = int(nm_m.group(3))
            near_min = {
                'total': int(nm_m.group(1)), 'min': int(nm_m.group(2)), 'non_min': non_min_total,
                'ratio': float(nm_m.group(4)),
                'path_used': int(nm_path.group(1)) if nm_path else 0, 'path_total': int(nm_path.group(2)) if nm_path else 0,
                'path_pct': float(nm_path.group(3)) if nm_path else 0.0,
            }
            for dm in re.finditer(r'^\s+(\S+)\s+->\s+(\S+(?:\([^)]*\))?)\s+\[(\w+)\]:\s*(\d+)\s+avg_sat=([\d.]+)%', nm_block, re.MULTILINE):
                count = int(dm.group(4))
                nm_dir_entries.append({
                    'src': canonical_router_name(dm.group(1)), 'dst': canonical_router_name(dm.group(2)),
                    'link_type': dm.group(3), 'count': count,
                    'pct': count / non_min_total * 100 if non_min_total > 0 else 0.0, 'avg_sat': float(dm.group(5)),
                })

    return {
        'num_kernels':      num_kernels,
        'tot_cycle':        tot_cycle,
        'l2_bw_total':      l2_bw_total,
        'l2_accesses':      total_l2_acc,
        'l2_misses':        total_l2_miss,
        'l2_hit_rate':      l2_hit_rate,
        'total_traversals': total_traversals,
        'pkt_lat_avg':      pkt_lat_avg,
        'pkt_lat_max':      pkt_lat_max,
        'inj_rate_avg':     inj_rate_avg,
        'inj_rate_max':     inj_rate_max,
        'link_util':        link_util,
        'agg_sat':          agg_sat,
        'dir_entries':      dir_entries,
        'near_min':         near_min,
        'nm_dir_entries':   nm_dir_entries,
    }


# ── CSV helpers ───────────────────────────────────────────────────────────────

SUMMARY_FIELDS = (
    ['benchmark', 'struct', 'bw', 'routing', 'num_kernels',
     'tot_cycle', 'l2_bw_total',
     'l2_accesses', 'l2_misses', 'l2_hit_rate',
     'total_traversals', 'pkt_lat_avg', 'pkt_lat_max', 'inj_rate_avg', 'inj_rate_max']
    + [f'{lt}_count' for lt in LINK_TYPES]
    + [f'{lt}_pct'   for lt in LINK_TYPES]
    + [f'{lt}_sat'   for lt in SAT_LINK_TYPES]
    + ['nm_total', 'nm_min', 'nm_non_min', 'nm_ratio',
       'nm_path_used', 'nm_path_total', 'nm_path_pct']
)

DIR_FIELDS = ['benchmark', 'struct', 'bw', 'routing', 'src', 'dst', 'count', 'pct', 'avg_sat']
NM_DIR_FIELDS = ['benchmark', 'struct', 'bw', 'routing', 'src', 'dst', 'link_type', 'count', 'pct', 'avg_sat']

def _summary_row(benchmark, struct, bw, rk, data):
    row = {
        'benchmark':        benchmark,
        'struct':           struct,
        'bw':               bw,
        'routing':          rk,
        'num_kernels':      data['num_kernels'],
        'tot_cycle':        data['tot_cycle'],
        'l2_bw_total':      data['l2_bw_total'],
        'l2_accesses':      data['l2_accesses'],
        'l2_misses':        data['l2_misses'],
        'l2_hit_rate':      f"{data['l2_hit_rate']:.6f}" if data['l2_hit_rate'] is not None else '',
        'total_traversals': data['total_traversals'],
        'pkt_lat_avg':      f"{data['pkt_lat_avg']:.4f}" if data['pkt_lat_avg'] is not None else '',
        'pkt_lat_max':      f"{data['pkt_lat_max']:.4f}" if data['pkt_lat_max'] is not None else '',
        'inj_rate_avg':     f"{data['inj_rate_avg']:.6f}" if data['inj_rate_avg'] is not None else '',
        'inj_rate_max':     f"{data['inj_rate_max']:.6f}" if data['inj_rate_max'] is not None else '',
    }
    for lt in LINK_TYPES:
        lu = data['link_util'].get(lt, {'count': 0, 'pct': 0.0})
        row[f'{lt}_count'] = lu['count']
        row[f'{lt}_pct']   = f"{lu['pct']:.4f}"
    for lt in SAT_LINK_TYPES:
        row[f'{lt}_sat'] = f"{data['agg_sat'].get(lt, 0.0):.4f}"
    nm = data['near_min']
    if nm:
        row.update({
            'nm_total':      nm['total'], 'nm_min': nm['min'], 'nm_non_min': nm['non_min'],
            'nm_ratio':      f"{nm['ratio']:.4f}",
            'nm_path_used':  nm['path_used'], 'nm_path_total': nm['path_total'],
            'nm_path_pct':   f"{nm['path_pct']:.4f}",
        })
    return row

# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--benchmark', required=True, help='Benchmark spec, e.g. rodinia-3.1:bfs-rodinia-3.1')
    ap.add_argument('--structure', nargs='+', choices=list(STRUCTURES), metavar='STRUCT', required=True)
    ap.add_argument('--bandwidth', nargs='+', choices=list(BANDWIDTHS), metavar='BW', required=True)
    ap.add_argument('--routing', nargs='+', choices=ROUTING_CHOICES, metavar='ROUTING', default=['baseline'])
    ap.add_argument('--near-min-k', nargs='+', type=int, metavar='K', default=[2])
    ap.add_argument('--near-min-p', nargs='+', type=float, metavar='P', default=[1.0])
    ap.add_argument('--all-configs', action='store_true')
    ap.add_argument('--output', '-o', default=None)
    ap.add_argument('--sim-run-dir', default=None)
    args = ap.parse_args()

    sim_run_base = os.path.abspath(args.sim_run_dir) if args.sim_run_dir else SIM_RUN_BASE

    if args.all_configs: structures, bandwidths = list(STRUCTURES), list(BANDWIDTHS)
    else: structures, bandwidths = args.structure, args.bandwidth

    routing_keys: list = []
    for r in args.routing:
        if r == 'near_min_adaptive':
            for k in args.near_min_k:
                for p in args.near_min_p: routing_keys.append(routing_to_key(r, near_min_k=k, near_min_p=p))
        elif r == 'near_min_random':
            for k in args.near_min_k: routing_keys.append(routing_to_key(r, near_min_k=k))
        else: routing_keys.append(r)
    routing_keys = list(dict.fromkeys(routing_keys))

    bench = args.benchmark
    if args.output is None: args.output = make_csv_prefix(bench, structures, bandwidths, routing_keys)

    combos = [(s, b, rk) for s in structures for b in bandwidths for rk in routing_keys]
    app_name = bench.split(':')[-1] if ':' in bench else bench
    active_args = get_active_args(bench)

    summary_rows, dir_rows, nm_dir_rows, missing = [], [], [], []

    for struct, bw, rk in combos:
        config_name = f"{struct}_{bw}_{rk}-SASS"
        label       = f"{struct} / {bw} / {rk}"
        log_path = find_log_file(app_name, config_name, sim_run_base, active_args)
        if log_path is None:
            missing.append(label)
            continue

        data = parse_log(log_path)
        if data is None:
            missing.append(label)
            continue

        summary_rows.append(_summary_row(bench, struct, bw, rk, data))

        for e in data['dir_entries']:
            dir_rows.append({'benchmark': bench, 'struct': struct, 'bw': bw, 'routing': rk, **e})
        for e in data['nm_dir_entries']:
            nm_dir_rows.append({'benchmark': bench, 'struct': struct, 'bw': bw, 'routing': rk, **e})

    # ── Write CSV files ───────────────────────────────────────────────────
    with open(args.output + '_summary.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS, extrasaction='ignore')
        w.writeheader(); w.writerows(summary_rows)
    with open(args.output + '_directions.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=DIR_FIELDS, extrasaction='ignore')
        w.writeheader(); w.writerows(dir_rows)
    with open(args.output + '_nm_directions.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=NM_DIR_FIELDS, extrasaction='ignore')
        w.writeheader(); w.writerows(nm_dir_rows)

    if missing:
        print(f"\n[WARN] Missing logs ({len(missing)}):")
        for m in missing: print(f"  {m}")

if __name__ == '__main__':
    main()