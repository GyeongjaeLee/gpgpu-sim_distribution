"""
Load experiment parameters from experiments.csv.

CSV format (tab-separated, two sections):

  Structure  Core  HBM  SM_per_Xbar
  B100_Local  1    4    74
  ...

  Bandwidth         GPU-to-GPU  GPU-to-HBM  HBM-to-HBM  TSV   L2_per_HBM
  B200+HBM3e        5           0.5         0.5          0.5   32
  ...

Rules:
  - Section header rows start with "Structure" or "Bandwidth".
  - Blank rows are ignored.
  - Structure:  Core = num_xbars,  HBM = total HBM stacks K.
                hbm_per_side = K // (num_xbars * 2).
                SM_per_Xbar  = SMs per Xbar (concentration).
  - Bandwidth:  L2_per_HBM  = L2 slices per HBM stack.
  - Adding a new row automatically exposes it as a CLI --structure /
    --bandwidth choice in all run/plot/gen scripts.
"""

import csv
import os
import re as _re
from collections import OrderedDict


# ── Abbreviation tables (shared by all run/parse/plot scripts) ────────────────

STRUCT_ABBREV = {
    "B100_Local":        "BL",
    "H100":              "H1",
    "B100_Global":       "BG",
    "B100_Core_Rotate":  "BCR",
    "Rubin_Ultra":       "RU",
}

BW_ABBREV = {
    "B200+HBM3e":       "H3",
    "Rubin_Ultra+HBM4": "H4",
    "Shoreline_1x":     "S1",
    "Shoreline_1.5x":   "S15",
    "Shoreline_2x":     "S2",
    "G2G_1.5x":         "G15",
    "G2G_1x":           "G1",
    "G2G_0.8x":         "G08",
    "G2G_0.5x":         "G05",
}

ROUTING_ABBREV = {
    "baseline":      "bas",
    "min_oblivious": "mino",
    "min_adaptive":  "mina",
    "ugal":          "ug",
    "valiant":       "val",
}


def bench_abbrev(benchmark: str) -> str:
    """
    Extract the distinctive part of a benchmark spec by removing
    tokens shared with the group name.

    "rodinia-3.1:bfs-rodinia-3.1"                 → "bfs"
    "polybench:polybench-gemm"                    → "gemm"
    "Deepbench_nvidia_tencore:gemm_bench-tencore" → "gemm_bench"
    "GPU_Microbenchmark:mem_bw"                   → "mem_bw"
    """
    if ':' not in benchmark:
        return benchmark
    group, app = benchmark.split(':', 1)
    group_tokens = set(_re.split(r'[-_.]', group.lower()))
    app_parts = _re.split(r'[-_.]', app)
    kept = [p for p in app_parts
            if p.lower() not in group_tokens and not _re.match(r'^\d+$', p)]
    return '_'.join(kept) if kept else app


def route_abbrev(rk: str) -> str:
    """Abbreviate a routing key (including near_min_pX.X variants)."""
    if rk in ROUTING_ABBREV:
        return ROUTING_ABBREV[rk]
    m = _re.match(r"near_min_p([\d.]+)", rk)
    if m:
        return "nm" + m.group(1).replace(".", "")
    return rk[:5]


def make_csv_prefix(benchmark: str, structures: list,
                    bandwidths: list, routing_keys: list) -> str:
    """
    Build an auto-generated CSV output prefix from option abbreviations.

    "rodinia-3.1:bfs-rodinia-3.1", ["B100_Global"], ["B200+HBM3e"], ["baseline"]
    → "bfs_BG_H3_bas"
    """
    ba  = bench_abbrev(benchmark)
    sa  = '+'.join(STRUCT_ABBREV.get(s, s[:3]) for s in structures)
    bwa = '+'.join(BW_ABBREV.get(b, b[:3])     for b in bandwidths)
    ra  = '+'.join(route_abbrev(rk)             for rk in routing_keys)
    return f"{ba}_{sa}_{bwa}_{ra}"


def load_experiments(csv_path: str = None):
    """
    Parse experiments.csv and return (STRUCTURES, BANDWIDTHS).

    STRUCTURES : OrderedDict[str, {
                    "num_xbars":    int,
                    "hbm_per_side": int,
                    "sm_per_xbar":  int,   # SMs per Xbar router
                 }]
    BANDWIDTHS : OrderedDict[str, {
                    "gpu_gpu":    float,
                    "gpu_hbm":    float,
                    "hbm_hbm":   float,
                    "tsv":        float,
                    "l2_per_hbm": int,    # L2 slices per HBM stack
                 }]
    """
    if csv_path is None:
        csv_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "experiments.csv")

    structures: OrderedDict = OrderedDict()
    bandwidths: OrderedDict = OrderedDict()
    section = None

    with open(csv_path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for raw_row in reader:
            row = [c.strip() for c in raw_row]
            while row and not row[-1]:
                row.pop()
            if not row:
                continue

            first = row[0]
            if first == "Structure":
                section = "structure"
                continue
            if first == "Bandwidth":
                section = "bandwidth"
                continue

            if section == "structure" and len(row) >= 3:
                try:
                    num_xbars   = int(row[1])
                    total_hbm   = int(row[2])
                    sm_per_xbar = int(row[3]) if len(row) >= 4 else 74
                except ValueError:
                    continue
                structures[first] = {
                    "num_xbars":    num_xbars,
                    "hbm_per_side": total_hbm // (num_xbars * 2),
                    "sm_per_xbar":  sm_per_xbar,
                }

            elif section == "bandwidth" and len(row) >= 5:
                try:
                    bandwidths[first] = {
                        "gpu_gpu":    float(row[1]),
                        "gpu_hbm":    float(row[2]),
                        "hbm_hbm":    float(row[3]),
                        "tsv":        float(row[4]),
                        "l2_per_hbm": int(row[5]) if len(row) >= 6 else 32,
                    }
                except ValueError:
                    continue

    if not structures:
        raise ValueError(f"No structures parsed from {csv_path}")
    if not bandwidths:
        raise ValueError(f"No bandwidths parsed from {csv_path}")

    return structures, bandwidths
