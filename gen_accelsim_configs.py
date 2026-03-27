#!/usr/bin/env python3
"""
Generate AccelSim (GPGPU-Sim) config pairs for HBMNet experiments.

For each (structure × bandwidth × routing) combination, writes:

  <out_dir>/<structure>/<bandwidth>/<routing_key>/
      gpgpusim.config              – L2-slice-scaled copy of base gpgpusim.config
      config_blackwell_islip.icnt  – topology / bandwidth / routing patched .icnt

Scaling rule
------------
  num_l2_slices = num_xbars * hbm_per_side * 2 * L2_PER_HBM   (L2_PER_HBM = 32)
  gpgpu_n_mem   = num_l2_slices / SUBPART_PER_MEM               (SUBPART_PER_MEM = 2)

  gpgpu_n_sub_partition_per_mchannel stays 2 (must be power-of-2).
  H100 (96 mem channels) may need manual tuning — see note in code.

Bandwidth columns (match experiments.csv)
-----------------------------------------
  gpu_gpu  → xbar_xbar_bandwidth
  gpu_hbm  → xbar_hbm_bandwidth = xbar_mc_bandwidth  (SM↔HBM hit / SM↔MC miss entry)
  hbm_hbm  → mc_mc_bandwidth                         (MC↔MC fabric)
  tsv      → mc_hbm_bandwidth                        (TSV: MC↔HBM stack, separate from fabric!)

  All *_bandwidth values = column_value × PORTS_PER_UNIT (default 7).

TSV note
--------
  When tsv > ~3 you may want hbm_internal_speedup > internal_speedup so the
  HBM router crossbar is not the bottleneck instead of the TSV links.

Routing keys
------------
  baseline         routing_function = baseline
  min_oblivious    routing_function = hybrid, hybrid_routing = min_oblivious
  min_adaptive     routing_function = hybrid, hybrid_routing = min_adaptive
  near_min_p<P>    routing_function = hybrid, hybrid_routing = near_min_adaptive,
                   near_min_penalty = P
  ugal             routing_function = hybrid, hybrid_routing = ugal
  valiant          routing_function = hybrid, hybrid_routing = valiant

Examples
--------
  # All structures × two bandwidths × three routings
  python3 gen_accelsim_configs.py \\
      --structure B100_Global Rubin_Ultra \\
      --bandwidth B100+HBM3e B100+HBM4e \\
      --routing baseline min_adaptive near_min_adaptive \\
      --near-min-p 1.0

  # Every combination defined in the tables
  python3 gen_accelsim_configs.py --all --routing baseline min_adaptive

  # Dry-run to preview paths
  python3 gen_accelsim_configs.py --all --routing baseline --dry-run
"""

import argparse
import os
import re
import sys
from collections import OrderedDict
from typing import Optional

# ── Constants ─────────────────────────────────────────────────────────────────

PORTS_PER_UNIT   = 7    # parallel ports per bandwidth-unit
L2_PER_HBM       = 32   # L2 slices per HBM stack (kept fixed across structures)
SUBPART_PER_MEM  = 2    # gpgpu_n_sub_partition_per_mchannel (must be power-of-2)

# ── Topology table ────────────────────────────────────────────────────────────
# K (total HBM stacks) = num_xbars * hbm_per_side * 2
# num_l2_slices = K * L2_PER_HBM
# gpgpu_n_mem   = num_l2_slices // SUBPART_PER_MEM

STRUCTURES = OrderedDict([
    #  key                  xbars  hbm/side  gpgpu_n_mem  num_l2_slices
    ("B100_Local",       {"num_xbars": 1, "hbm_per_side": 2, "n_mem":  64, "num_l2_slices": 128}),
    ("H100",             {"num_xbars": 1, "hbm_per_side": 3, "n_mem":  96, "num_l2_slices": 192}),
    ("B100_Global",      {"num_xbars": 2, "hbm_per_side": 2, "n_mem": 128, "num_l2_slices": 256}),
    ("B100_Core_Rotate", {"num_xbars": 2, "hbm_per_side": 3, "n_mem": 192, "num_l2_slices": 384}),
    ("Rubin_Ultra",      {"num_xbars": 4, "hbm_per_side": 2, "n_mem": 256, "num_l2_slices": 512}),
])

# ── Bandwidth table ───────────────────────────────────────────────────────────
# Matches experiments.csv columns (Structure | Core | HBM | … | Bandwidth | GPU-GPU | GPU-HBM | HBM-HBM | TSV)

BANDWIDTHS = OrderedDict([
    # key               gpu_gpu  gpu_hbm  hbm_hbm  tsv
    ("B100+HBM3e",    {"gpu_gpu": 10,  "gpu_hbm": 1.0,  "hbm_hbm": 1.0,  "tsv": 1.0}),
    ("B100+HBM4",     {"gpu_gpu": 10,  "gpu_hbm": 2.0,  "hbm_hbm": 2.0,  "tsv": 2.0}),   # HBM3e × 2
    ("B100+HBM4e",    {"gpu_gpu": 10,  "gpu_hbm": 4.0,  "hbm_hbm": 4.0,  "tsv": 4.0}),
    ("Shoreline_1x",  {"gpu_gpu": 10,  "gpu_hbm": 3.33, "hbm_hbm": 3.89, "tsv": 3.33}),
    ("Shoreline_1.5x",{"gpu_gpu": 10,  "gpu_hbm": 3.33, "hbm_hbm": 3.89, "tsv": 5.0}),
    ("Shoreline_2x",  {"gpu_gpu": 10,  "gpu_hbm": 3.33, "hbm_hbm": 3.89, "tsv": 6.67}),
    ("G2G_1.5x",      {"gpu_gpu": 15,  "gpu_hbm": 3.33, "hbm_hbm": 3.89, "tsv": 5.0}),
    ("G2G_1x",        {"gpu_gpu": 10,  "gpu_hbm": 3.33, "hbm_hbm": 3.89, "tsv": 5.0}),
    ("G2G_0.8x",      {"gpu_gpu": 8,   "gpu_hbm": 3.33, "hbm_hbm": 3.89, "tsv": 5.0}),
    ("G2G_0.5x",      {"gpu_gpu": 5,   "gpu_hbm": 3.33, "hbm_hbm": 3.89, "tsv": 5.0}),
])

# ── Routing ───────────────────────────────────────────────────────────────────

ROUTING_CHOICES = [
    "baseline", "min_oblivious", "min_adaptive",
    "near_min_adaptive", "ugal", "valiant",
]


def routing_to_key(routing: str, near_min_p: Optional[float] = None) -> str:
    if routing == "near_min_adaptive":
        p = near_min_p if near_min_p is not None else 1.0
        return f"near_min_p{p:.1f}"
    return routing


def routing_key_to_overrides(key: str) -> dict:
    if key == "baseline":
        return {"routing_function": "baseline"}
    ov: dict = {"routing_function": "hybrid"}
    m = re.match(r"near_min_p([\d.]+)$", key)
    if m:
        ov["hybrid_routing"]   = "near_min_adaptive"
        ov["near_min_penalty"] = m.group(1)
        return ov
    ov["hybrid_routing"] = key
    return ov


# ── Config patching ───────────────────────────────────────────────────────────

def _patch_file(src: str, dst: str, overrides: "dict[str, str]") -> None:
    """Rewrite src → dst replacing 'key = ...' lines per overrides, appending unknowns."""
    with open(src) as f:
        lines = f.readlines()
    written: set[str] = set()
    result = []
    for line in lines:
        s = line.strip()
        if s.startswith("//") or not s:
            result.append(line)
            continue
        m = re.match(r"^(\w+)\s*=", s)
        if m and m.group(1) in overrides:
            result.append(f"{m.group(1)} = {overrides[m.group(1)]};\n")
            written.add(m.group(1))
            continue
        result.append(line)
    for k, v in overrides.items():
        if k not in written:
            result.append(f"{k} = {v};\n")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as f:
        f.writelines(result)


def _patch_gpgpusim(src: str, dst: str, n_mem: int) -> None:
    """Rewrite src → dst updating -gpgpu_n_mem."""
    with open(src) as f:
        lines = f.readlines()
    result = []
    for line in lines:
        if re.match(r"^-gpgpu_n_mem\s+\d+", line.strip()):
            result.append(f"-gpgpu_n_mem {n_mem}\n")
        else:
            result.append(line)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as f:
        f.writelines(result)


# ── Override builders ─────────────────────────────────────────────────────────

def _icnt_overrides(struct: str, bw: str, routing_key: str) -> "dict[str, str]":
    s  = STRUCTURES[struct]
    b  = BANDWIDTHS[bw]
    p  = PORTS_PER_UNIT
    ov = {
        # topology
        "num_xbars":           str(s["num_xbars"]),
        "hbm_per_side":        str(s["hbm_per_side"]),
        "num_l2_slices":       str(s["num_l2_slices"]),
        # bandwidth  (TSV = mc_hbm, separate from MC fabric = mc_mc)
        "xbar_xbar_bandwidth": str(round(b["gpu_gpu"] * p)),
        "xbar_hbm_bandwidth":  str(round(b["gpu_hbm"] * p)),
        "xbar_mc_bandwidth":   str(round(b["gpu_hbm"] * p)),
        "mc_hbm_bandwidth":    str(round(b["tsv"]     * p)),   # TSV bottleneck
        "mc_mc_bandwidth":     str(round(b["hbm_hbm"] * p)),
    }
    ov.update(routing_key_to_overrides(routing_key))
    return ov


# ── Paths ─────────────────────────────────────────────────────────────────────

_HERE         = os.path.dirname(os.path.abspath(__file__))
_BASE_DIR     = os.path.join(_HERE, "configs", "tested-cfgs", "SM100_B200_fabric")
BASE_ICNT     = os.path.join(_BASE_DIR, "config_blackwell_islip.icnt")
BASE_GPGPUSIM = os.path.join(_BASE_DIR, "gpgpusim.config")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "--structure", nargs="+", choices=list(STRUCTURES), metavar="STRUCT",
        help="Structure(s). Choices: " + ", ".join(STRUCTURES),
    )
    ap.add_argument(
        "--bandwidth", nargs="+", choices=list(BANDWIDTHS), metavar="BW",
        help="Bandwidth class(es). Choices: " + ", ".join(BANDWIDTHS),
    )
    ap.add_argument(
        "--routing", nargs="+", choices=ROUTING_CHOICES, metavar="ROUTING",
        default=["baseline"],
        help="Routing function(s). Default: baseline",
    )
    ap.add_argument(
        "--near-min-p", nargs="+", type=float, metavar="P", default=[1.0],
        help="near_min_penalty value(s) for near_min_adaptive. Default: 1.0",
    )
    ap.add_argument(
        "--all", action="store_true",
        help="Use all structures and bandwidths defined in the tables",
    )
    ap.add_argument(
        "--out-dir",
        default=os.path.join(_HERE, "configs", "experiments"),
        metavar="DIR",
        help="Output base directory (default: configs/experiments/)",
    )
    ap.add_argument(
        "--dry-run", action="store_true",
        help="Print planned output paths without writing files",
    )
    args = ap.parse_args()

    if args.all:
        structures = list(STRUCTURES)
        bandwidths = list(BANDWIDTHS)
    else:
        if not args.structure or not args.bandwidth:
            ap.error("Specify --structure and --bandwidth, or use --all")
        structures = args.structure
        bandwidths = args.bandwidth

    # Expand near_min_adaptive into per-P keys
    routing_keys: list[str] = []
    for r in args.routing:
        if r == "near_min_adaptive":
            for p in args.near_min_p:
                routing_keys.append(routing_to_key(r, p))
        else:
            routing_keys.append(r)

    combos = [
        (s, b, rk)
        for s  in structures
        for b  in bandwidths
        for rk in routing_keys
    ]

    print(f"Generating {len(combos)} config pair(s) → {args.out_dir}/")
    print()

    for struct, bw, rk in combos:
        s     = STRUCTURES[struct]
        n_mem = s["n_mem"]
        n_l2  = s["num_l2_slices"]

        dir_name  = f"SM100_{struct}_{bw}_{rk}"
        out_dir   = os.path.join(args.out_dir, dir_name)
        icnt_dst  = os.path.join(out_dir, "config_blackwell_islip.icnt")
        gpgpu_dst = os.path.join(out_dir, "gpgpusim.config")

        label = f"  {dir_name:<60}  n_mem={n_mem:3d}  n_l2={n_l2:3d}"

        if args.dry_run:
            print(f"[DRY] {label}")
            continue

        ov = _icnt_overrides(struct, bw, rk)
        _patch_file(BASE_ICNT, icnt_dst, ov)
        _patch_gpgpusim(BASE_GPGPUSIM, gpgpu_dst, n_mem)
        print(f"[OK]  {label}")

    if not args.dry_run:
        print(f"\nDone. {len(combos)} pair(s) written.")


if __name__ == "__main__":
    main()
