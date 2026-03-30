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
  - Global metrics (gpu_tot_sim_cycle, L2_BW_total) are taken from the LAST
    occurrence in the log (cumulative totals).

Aggregation rules:
  L2 hit rate   : sum(L2_total_cache_accesses/misses) across kernels
  Link util     : sum counts, recompute pct
  Link avg_sat  : weighted average (weight = per-kernel link-type traversal count)
  Per-direction : same weighted scheme
  Near-min      : summed across kernels (near_min_adaptive/random only)

Outputs (prefix auto-generated from options if --output is omitted):
  {prefix}_summary.csv       – one row per (struct, bw, routing)
  {prefix}_directions.csv    – one row per (struct, bw, routing, src, dst)
  {prefix}_nm_directions.csv – near-min direction breakdown (near_min_adaptive/random only)

Usage examples:

  # Single structure/bandwidth/routing, auto-generated prefix
  python parse_accelsim_results.py \
      --benchmark rodinia-3.1:bfs-rodinia-3.1 \
      --structure B100_Global \
      --bandwidth B200+HBM3e \
      --routing baseline min_adaptive near_min_adaptive near_min_random \
      --near-min-k 2 \
      --near-min-p 0.0 1.0

  # Multiple structures and bandwidths with explicit output prefix
  python parse_accelsim_results.py \
      --benchmark polybench:polybench-gemm \
      --structure B100_Local B100_Global H100 \
      --bandwidth B200+HBM3e Shoreline_1x \
      --routing baseline near_min_adaptive \
      --near-min-k 2 --near-min-p 1.0 \
      --output gemm_results

  # All configs from experiments.csv
  python parse_accelsim_results.py \
      --benchmark GPU_Microbenchmark:mem_bw \
      --all-configs \
      --routing baseline min_adaptive fixed_min
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
SAT_LINK_TYPES = ["XBAR_XBAR", "XBAR_MC", "MC_MC"]  # only these have avg_sat data


# ── Helpers ───────────────────────────────────────────────────────────────────

def routing_to_key(routing: str, near_min_k: Optional[int] = None, near_min_p: Optional[float] = None) -> str:
    """Convert routing name (+ optional K and P values) to a directory key."""
    if routing == "near_min_adaptive":
        k = near_min_k if near_min_k is not None else 2
        p = near_min_p if near_min_p is not None else 1.0
        return f"{routing}_nmk{k}_nmp{p:.1f}"
    elif routing == "near_min_random":
        k = near_min_k if near_min_k is not None else 2
        return f"{routing}_nmk{k}"
    return routing


def canonical_router_name(name: str) -> str:
    """Strip coordinate suffix: MC0(c0r0) → MC0"""
    m = re.match(r'^(Xbar\d+|MC\d+)', name)
    return m.group(1) if m else name


def _parse_active_args_fallback(yaml_path: str, group_name: str, app_name: str) -> list:
    """
    Pure-Python fallback parser for define-all-apps.yml.

    Reads the file line-by-line and extracts non-commented `args:` values
    belonging to the requested group / app.  Works without pyyaml.

    Relevant YAML shape (indentation signals hierarchy):
        {group_name}:           ← group header (col 0)
            execs:              ← col 4
                - {app_name}:   ← col 8  (list item with app key)
                    - args:  VAL        ← col 12 (active)
                    # - args:  VAL      ← commented out (skipped)
    """
    args_list = []
    in_group = False
    in_app   = False
    app_col  = None   # indentation column where the app "- {app_name}:" appears

    with open(yaml_path) as f:
        for raw_line in f:
            line    = raw_line.rstrip("\n")
            stripped = line.lstrip()

            # Blank lines: keep state
            if not stripped:
                continue

            # Commented lines: skip entirely
            if stripped.startswith("#"):
                continue

            col = len(line) - len(stripped)  # indentation level

            # ── group detection ──────────────────────────────────────────
            if not in_group:
                if col == 0 and stripped == group_name + ":":
                    in_group = True
                continue

            # Inside the group — a new col-0 key means we left the group
            if col == 0:
                break

            # ── app detection ────────────────────────────────────────────
            if not in_app:
                # Match:  - {app_name}:
                m = re.match(r'^(\s*)-\s+' + re.escape(app_name) + r'\s*:\s*$', line)
                if m:
                    in_app  = True
                    app_col = col
                continue

            # Inside the app block — detect when we leave it
            # A new exec list item at the same indentation level → left app block
            if col <= app_col and re.match(r'^\s*-\s+\S', line):
                break

            # Match active args lines:  - args:  VALUE
            m = re.match(r'^\s+-\s+args:\s+(.+)$', line)
            if m:
                args_list.append(m.group(1).strip())

    return args_list


def get_active_args(benchmark: str, yaml_path: str = DEFINE_APPS_YAML) -> list:
    """
    Return the list of active (non-commented) args strings for the given
    benchmark spec ("group:app_name") from define-all-apps.yml.

    Example:
        "rodinia-3.1:bfs-rodinia-3.1" → ["./data/graph1MW_6.txt"]

    Tries pyyaml first; falls back to the built-in line parser.
    Returns [] if the YAML file is missing or the benchmark is not found.
    """
    if not os.path.isfile(yaml_path):
        return []

    parts = benchmark.split(":", 1)
    if len(parts) != 2:
        return []
    group_name, app_name = parts

    # ── Try pyyaml (preferred) ────────────────────────────────────────────
    try:
        import yaml  # type: ignore[import]
        with open(yaml_path) as f:
            content = yaml.safe_load(f)
        group = content.get(group_name)
        if group:
            for exec_entry in group.get("execs", []):
                if not isinstance(exec_entry, dict):
                    continue
                app_data = exec_entry.get(app_name)
                if app_data is None:
                    continue
                return [
                    str(cfg["args"]).strip()
                    for cfg in app_data
                    if isinstance(cfg, dict) and "args" in cfg
                ]
        return []
    except ImportError:
        pass  # fall through to built-in parser

    # ── Built-in fallback ─────────────────────────────────────────────────
    return _parse_active_args_fallback(yaml_path, group_name, app_name)


def sanitize_args_for_path(args: str) -> str:
    """
    Convert an args string to the AccelSim run-directory name.

    Exactly mirrors AccelSim's common.get_argfoldername():
      re.sub(r"[^a-z^A-Z^0-9]", "_", args.strip())
    with MD5 fallback for args longer than 256 characters.

    Example: "./data/graph1MW_6.txt" → "__data_graph1MW_6_txt"
    """
    if not args or not args.strip():
        return "NO_ARGS"
    s = args.strip()
    foldername = re.sub(r"[^a-z^A-Z^0-9]", "_", s)
    if len(s) > 256:
        foldername = "hashed_args_" + hashlib.md5(s.encode()).hexdigest()
    return foldername


def find_log_file(app_name: str, config_name: str, sim_run_base: str,
                  active_args: Optional[list] = None) -> Optional[str]:
    """
    Search for a simulation log file and return the path with the largest .oN suffix.

    AccelSim places output at:
        {sim_run_base}/{app_name}/{sanitized_args}/{config_name}/{app_name}*.o{N}
    where sanitized_args = sanitize_args_for_path(args)
    e.g. "./data/graph1MW_6.txt" → "__data_graph1MW_6_txt"

    Falls back to a broad wildcard glob if nothing found via active_args.
    """
    candidates = []

    if active_args:
        for args in active_args:
            args_dir = sanitize_args_for_path(args)
            pattern  = os.path.join(sim_run_base, app_name,
                                    args_dir, config_name, f"{app_name}*.o*")
            for path in glob.glob(pattern):
                m = re.search(r'\.o(\d+)$', path)
                if m:
                    candidates.append((int(m.group(1)), path))

    # Fallback: any single-level directory under app_name
    if not candidates:
        pattern = os.path.join(sim_run_base, app_name, "*",
                               config_name, f"{app_name}*.o*")
        for path in glob.glob(pattern):
            m = re.search(r'\.o(\d+)$', path)
            if m:
                candidates.append((int(m.group(1)), path))

    if not candidates:
        return None
    return max(candidates, key=lambda x: x[0])[1]


# ── Log parser ────────────────────────────────────────────────────────────────

def parse_log(path: str) -> Optional[dict]:
    """
    Parse a single AccelSim log and return metrics.

    gpu_tot_sim_cycle, L2_BW_total, link utilization, per-direction stats,
    avg_sat, and near-min counters are all cumulative globals — only the LAST
    printed occurrence of each is used (final totals).

    L2 cache stats (accesses / misses) are printed per-kernel incrementally,
    so they are summed across all kernel outputs.
    """
    with open(path) as f:
        text = f.read()

    # ── Global metrics: last occurrence only (cumulative) ────────────────
    all_cycles = re.findall(r'^gpu_tot_sim_cycle\s*=\s*(\d+)', text, re.MULTILINE)
    tot_cycle = int(all_cycles[-1]) if all_cycles else None

    all_bw = re.findall(r'^L2_BW_total\s*=\s*([\d.]+)', text, re.MULTILINE)
    l2_bw_total = float(all_bw[-1]) if all_bw else None

    # ── L2: sum all per-kernel occurrences ────────────────────────────────
    total_l2_acc  = sum(int(m.group(1)) for m in re.finditer(
        r'^L2_total_cache_accesses\s*=\s*(\d+)', text, re.MULTILINE))
    total_l2_miss = sum(int(m.group(1)) for m in re.finditer(
        r'^L2_total_cache_misses\s*=\s*(\d+)', text, re.MULTILINE))
    l2_hit_rate = (1 - total_l2_miss / total_l2_acc) if total_l2_acc > 0 else None

    # ── Link stats: LAST block only (cumulative totals) ───────────────────
    link_starts = [m.start() for m in re.finditer(r'=== Link Type Utilization ===', text)]
    if not link_starts:
        return None

    num_kernels     = len(link_starts)
    last_link_start = link_starts[-1]

    # Find end of last link block at "Miss link avg saturation: ..."
    m_end = re.search(r'Miss link avg saturation:.*', text[last_link_start:])
    last_link_end = last_link_start + m_end.end() if m_end else len(text)
    link_block = text[last_link_start:last_link_end]

    # ── Parse link-block ──────────────────────────────────────────────────
    link_util: dict = {}
    for lm in re.finditer(
            r'^\s+(XBAR_XBAR|XBAR_HBM|XBAR_MC|MC_HBM|MC_MC):\s*(\d+)\s*\(([\d.]+)%\)',
            link_block, re.MULTILINE):
        link_util[lm.group(1)] = {
            'count': int(lm.group(2)), 'pct': float(lm.group(3))}

    m2 = re.search(r'Total link traversals:\s*(\d+)', link_block)
    total_traversals = int(m2.group(1)) if m2 else 0

    dir_entries = []
    for dm in re.finditer(
            r'^\s{4}(\S+(?:\([^)]*\))?)\s*->\s*(\S+(?:\([^)]*\))?):\s*'
            r'count=(\d+)\s*\(([\d.]+)%\)\s*avg_sat=([\d.]+)%',
            link_block, re.MULTILINE):
        dir_entries.append({
            'src':     canonical_router_name(dm.group(1)),
            'dst':     canonical_router_name(dm.group(2)),
            'count':   int(dm.group(3)),
            'pct':     float(dm.group(4)),
            'avg_sat': float(dm.group(5)),
        })

    m2 = re.search(
        r'Miss link avg saturation:\s*'
        r'XBAR_XBAR=([\d.]+)%\s*XBAR_MC=([\d.]+)%\s*MC_MC=([\d.]+)%',
        link_block)
    agg_sat = {
        'XBAR_XBAR': float(m2.group(1)) if m2 else 0.0,
        'XBAR_MC':   float(m2.group(2)) if m2 else 0.0,
        'MC_MC':     float(m2.group(3)) if m2 else 0.0,
    }

    # ── Near-min: find last === Near-Min Adaptive Decisions === block ────
    nm_starts = [m.start() for m in re.finditer(
        r'=== Near-Min Adaptive Decisions ===', text)]
    near_min: Optional[dict] = None
    nm_dir_entries: list = []

    if nm_starts:
        last_nm_start = nm_starts[-1]
        # Block ends at the next === section header
        hdr_end = text.index('\n', last_nm_start) + 1
        next_sec = re.search(r'^===', text[hdr_end:], re.MULTILINE)
        nm_end   = hdr_end + next_sec.start() if next_sec else len(text)
        nm_block = text[last_nm_start:nm_end]

        nm_m = re.search(
            r'Near-min decisions:\s*total=(\d+)\s+min=(\d+)\s+non-min=(\d+)'
            r'\s+near-min ratio=([\d.]+)%',
            nm_block)
        nm_path = re.search(
            r'Near-min path usage:\s*(\d+)\s*/\s*(\d+)\s*\(([\d.]+)%\)',
            nm_block)

        if nm_m:
            non_min_total = int(nm_m.group(3))
            near_min = {
                'total':      int(nm_m.group(1)),
                'min':        int(nm_m.group(2)),
                'non_min':    non_min_total,
                'ratio':      float(nm_m.group(4)),
                'path_used':  int(nm_path.group(1))   if nm_path else 0,
                'path_total': int(nm_path.group(2))   if nm_path else 0,
                'path_pct':   float(nm_path.group(3)) if nm_path else 0.0,
            }

            # Direction breakdown:  Xbar1 -> MC2(c0r2) [XBAR_MC]: 2942 avg_sat=4.38962%
            for dm in re.finditer(
                    r'^\s+(\S+)\s+->\s+(\S+(?:\([^)]*\))?)\s+\[(\w+)\]:\s*'
                    r'(\d+)\s+avg_sat=([\d.]+)%',
                    nm_block, re.MULTILINE):
                count = int(dm.group(4))
                nm_dir_entries.append({
                    'src':       canonical_router_name(dm.group(1)),
                    'dst':       canonical_router_name(dm.group(2)),
                    'link_type': dm.group(3),
                    'count':     count,
                    'pct':       count / non_min_total * 100 if non_min_total > 0 else 0.0,
                    'avg_sat':   float(dm.group(5)),
                })

    return {
        'num_kernels':      num_kernels,
        'tot_cycle':        tot_cycle,
        'l2_bw_total':      l2_bw_total,
        'l2_accesses':      total_l2_acc,
        'l2_misses':        total_l2_miss,
        'l2_hit_rate':      l2_hit_rate,
        'total_traversals': total_traversals,
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
     'total_traversals']
    + [f'{lt}_count' for lt in LINK_TYPES]
    + [f'{lt}_pct'   for lt in LINK_TYPES]
    + [f'{lt}_sat'   for lt in SAT_LINK_TYPES]
    + ['nm_total', 'nm_min', 'nm_non_min', 'nm_ratio',
       'nm_path_used', 'nm_path_total', 'nm_path_pct']
)

DIR_FIELDS = ['benchmark', 'struct', 'bw', 'routing',
              'src', 'dst', 'count', 'pct', 'avg_sat']

NM_DIR_FIELDS = ['benchmark', 'struct', 'bw', 'routing',
                 'src', 'dst', 'link_type', 'count', 'pct', 'avg_sat']


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
            'nm_total':      nm['total'],
            'nm_min':        nm['min'],
            'nm_non_min':    nm['non_min'],
            'nm_ratio':      f"{nm['ratio']:.4f}",
            'nm_path_used':  nm['path_used'],
            'nm_path_total': nm['path_total'],
            'nm_path_pct':   f"{nm['path_pct']:.4f}",
        })
    return row


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument('--benchmark', required=True,
                    help='Benchmark spec, e.g. rodinia-3.1:bfs-rodinia-3.1')
    ap.add_argument('--structure', nargs='+', choices=list(STRUCTURES),
                    metavar='STRUCT', required=True)
    ap.add_argument('--bandwidth', nargs='+', choices=list(BANDWIDTHS),
                    metavar='BW', required=True)
    ap.add_argument('--routing', nargs='+', choices=ROUTING_CHOICES,
                    metavar='ROUTING', default=['baseline'])
    ap.add_argument('--near-min-k', nargs='+', type=int, metavar='K',
                    default=[2], help='Near-min routing budget (default: 2)')
    ap.add_argument('--near-min-p', nargs='+', type=float, metavar='P',
                    default=[1.0], help='Near-min routing penalty multiplier (default: 1.0)')
    ap.add_argument('--all-configs', action='store_true',
                    help='Use all structures and bandwidths from experiments.csv')
    ap.add_argument('--output', '-o', default=None,
                    help='Output CSV prefix (default: auto-generated from options, '
                         'e.g. bfs_BG_H3_bas)')
    ap.add_argument('--sim-run-dir', default=None,
                    help='Override sim_run_12.8 base path')
    args = ap.parse_args()

    sim_run_base = os.path.abspath(args.sim_run_dir) if args.sim_run_dir else SIM_RUN_BASE

    if args.all_configs:
        structures = list(STRUCTURES)
        bandwidths = list(BANDWIDTHS)
    else:
        structures = args.structure
        bandwidths = args.bandwidth

    # Expand routing keys dynamically based on combinations
    routing_keys: list = []
    for r in args.routing:
        if r == 'near_min_adaptive':
            for k in args.near_min_k:
                for p in args.near_min_p:
                    routing_keys.append(routing_to_key(r, near_min_k=k, near_min_p=p))
        elif r == 'near_min_random':
            for k in args.near_min_k:
                routing_keys.append(routing_to_key(r, near_min_k=k))
        else:
            routing_keys.append(r)
            
    # Remove duplicates
    routing_keys = list(dict.fromkeys(routing_keys))

    bench = args.benchmark
    if args.output is None:
        args.output = make_csv_prefix(bench, structures, bandwidths, routing_keys)

    combos = [
        (s, b, rk)
        for s in structures
        for b in bandwidths
        for rk in routing_keys
    ]

    app_name = bench.split(':')[-1] if ':' in bench else bench

    # Resolve active args from define-all-apps.yml
    active_args = get_active_args(bench)
    if active_args:
        print(f"Active args: {active_args}  (from define-all-apps.yml)")
    else:
        print(f"[WARN] No active args found for '{bench}' — using arg* glob fallback")

    print(f"Benchmark  : {bench}  (app={app_name})")
    print(f"sim_run    : {sim_run_base}")
    print(f"Combos     : {len(combos)}\n")

    summary_rows:  list = []
    dir_rows:      list = []
    nm_dir_rows:   list = []
    missing:       list = []

    for struct, bw, rk in combos:
        config_name = f"{struct}_{bw}_{rk}-SASS"
        label       = f"{struct} / {bw} / {rk}"

        log_path = find_log_file(app_name, config_name, sim_run_base, active_args)
        if log_path is None:
            arg_hint = sanitize_args_for_path(active_args[0]) if active_args else "*"
            print(f"[MISS] {label}  →  no log in {app_name}/{arg_hint}/{config_name}/")
            missing.append(label)
            continue

        print(f"[PARSE] {label}")
        print(f"        {log_path}")

        data = parse_log(log_path)
        if data is None:
            print(f"  [WARN] No network stats blocks found — skipping.")
            missing.append(label)
            continue

        print(f"  kernels={data['num_kernels']}  "
              f"cycles={data['tot_cycle']}  "
              f"L2_hit={data['l2_hit_rate']:.3f}" if data['l2_hit_rate'] is not None
              else f"  kernels={data['num_kernels']}  cycles={data['tot_cycle']}  L2_hit=N/A")

        summary_rows.append(_summary_row(bench, struct, bw, rk, data))

        for e in data['dir_entries']:
            dir_rows.append({
                'benchmark': bench,
                'struct':    struct,
                'bw':        bw,
                'routing':   rk,
                'src':       e['src'],
                'dst':       e['dst'],
                'count':     e['count'],
                'pct':       f"{e['pct']:.4f}",
                'avg_sat':   f"{e['avg_sat']:.4f}",
            })

        for e in data['nm_dir_entries']:
            nm_dir_rows.append({
                'benchmark': bench,
                'struct':    struct,
                'bw':        bw,
                'routing':   rk,
                'src':       e['src'],
                'dst':       e['dst'],
                'link_type': e['link_type'],
                'count':     e['count'],
                'pct':       f"{e['pct']:.4f}",
                'avg_sat':   f"{e['avg_sat']:.4f}",
            })

    # ── Write CSV files ───────────────────────────────────────────────────
    summary_path = args.output + '_summary.csv'
    dir_path     = args.output + '_directions.csv'

    with open(summary_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS, extrasaction='ignore')
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\nWrote: {summary_path}  ({len(summary_rows)} rows)")

    with open(dir_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=DIR_FIELDS, extrasaction='ignore')
        w.writeheader()
        w.writerows(dir_rows)
    print(f"Wrote: {dir_path}  ({len(dir_rows)} rows)")

    nm_dir_path = args.output + '_nm_directions.csv'
    with open(nm_dir_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=NM_DIR_FIELDS, extrasaction='ignore')
        w.writeheader()
        w.writerows(nm_dir_rows)
    print(f"Wrote: {nm_dir_path}  ({len(nm_dir_rows)} rows)")

    if missing:
        print(f"\n[WARN] Missing logs ({len(missing)}):")
        for m in missing:
            print(f"  {m}")


if __name__ == '__main__':
    main()
