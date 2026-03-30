#!/usr/bin/env python3
"""
Generate AccelSim (GPGPU-Sim) config pairs for HBMNet experiments.

For each (structure × bandwidth × routing) combination, writes:

  1) HW Configs (gpgpusim.config, .icnt):
     ./configs/tested-cfgs/SM100_<...>/
  
  2) Trace Configs (*.config):
     ../configs/tested-cfgs/SM100_<...>/
     
  3) YAML definitions:
     Inserts into ../../util/job_launching/configs/define-standard-cfgs.yml
     exactly below the 'B200:' block.
"""

import argparse
import os
import re
import shutil
from typing import Optional

from experiments_loader import load_experiments

# ── Constants ─────────────────────────────────────────────────────────────────

PORTS_PER_UNIT     = 7
N_SUB_PARTITION    = 2   # -gpgpu_n_sub_partition_per_mchannel (fixed)

# ── Topology / bandwidth tables loaded from experiments.csv ───────────────────
# All values derived from per-structure / per-bandwidth primitives:
#   K            = num_xbars * hbm_per_side * 2
#   n_clusters   = num_xbars * sm_per_xbar          (from structure)
#   n_mem        = l2_per_hbm * K / N_SUB_PARTITION (from bandwidth)
#   gpgpu_n_clusters = n_clusters
#   gpgpu_n_mem      = n_mem

_HERE = os.path.dirname(os.path.abspath(__file__))
STRUCTURES, BANDWIDTHS = load_experiments(os.path.join(_HERE, "experiments.csv"))

# ── Routing ───────────────────────────────────────────────────────────────────
ROUTING_CHOICES = [
    "baseline", "min_oblivious", "min_adaptive",
    "near_min_adaptive", "near_min_random", "fixed_min", 
    "ugal", "valiant",
]

def routing_to_key(routing: str, near_min_k: Optional[int] = None, near_min_p: Optional[float] = None) -> str:
    """Generate a unique directory key based on routing and its parameters."""
    if routing == "near_min_adaptive":
        k = near_min_k if near_min_k is not None else 2
        p = near_min_p if near_min_p is not None else 1.0
        return f"{routing}_nmk{k}_nmp{p:.1f}"
    elif routing == "near_min_random":
        k = near_min_k if near_min_k is not None else 2
        return f"{routing}_nmk{k}"
    return routing

def routing_key_to_overrides(key: str) -> dict:
    """Convert the routing key back into config overrides."""
    # Base overrides for all hybrid routing cases
    ov: dict = {"routing_function": "hybrid", "is_fabric": "1"}
    
    if key == "baseline":
        ov["is_fabric"] = "0"
        ov["hybrid_routing"] = "baseline"
        return ov
        
    m_adp = re.match(r"^near_min_adaptive_nmk(\d+)_nmp([\d.]+)$", key)
    if m_adp:
        ov["hybrid_routing"]   = "near_min_adaptive"
        ov["near_min_k"]       = m_adp.group(1)
        ov["near_min_penalty"] = m_adp.group(2)
        return ov

    m_rnd = re.match(r"^near_min_random_nmk(\d+)$", key)
    if m_rnd:
        ov["hybrid_routing"]   = "near_min_random"
        ov["near_min_k"]       = m_rnd.group(1)
        return ov
        
    # All other routings (min_adaptive, fixed_min, ugal, etc.)
    ov["hybrid_routing"] = key
    return ov

# ── Config patching ───────────────────────────────────────────────────────────
def _patch_file(src: str, dst: str, overrides: "dict[str, str]") -> None:
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

def _patch_gpgpusim(src: str, dst: str, n_mem: int, n_clusters: int) -> None:
    """Rewrite src → dst updating -gpgpu_n_mem and -gpgpu_n_clusters."""
    with open(src) as f:
        lines = f.readlines()
    result = []
    for line in lines:
        if re.match(r"^-gpgpu_n_mem\s+\d+", line.strip()):
            result.append(f"-gpgpu_n_mem {n_mem}\n")
        elif re.match(r"^-gpgpu_n_clusters\s+\d+", line.strip()):
            result.append(f"-gpgpu_n_clusters {n_clusters}\n")
        else:
            result.append(line)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as f:
        f.writelines(result)

def _icnt_overrides(struct: str, bw: str, routing_key: str) -> "dict[str, str]":
    s_info = STRUCTURES[struct]
    b      = BANDWIDTHS[bw]
    p      = PORTS_PER_UNIT
    ov = {
        "num_xbars":           str(s_info["num_xbars"]),
        "hbm_per_side":        str(s_info["hbm_per_side"]),
        "sm_per_xbar":         str(s_info["sm_per_xbar"]),
        "l2_per_hbm":          str(b["l2_per_hbm"]),
        "xbar_xbar_bandwidth": str(round(b["gpu_gpu"] * p)),
        "xbar_hbm_bandwidth":  str(round(b["gpu_hbm"] * p)),
        "xbar_mc_bandwidth":   str(round(b["gpu_hbm"] * p)),
        "mc_hbm_bandwidth":    str(round(b["tsv"]     * p)),
        "mc_mc_bandwidth":     str(round(b["hbm_hbm"] * p)),
    }
    ov.update(routing_key_to_overrides(routing_key))
    return ov

# ── Paths ─────────────────────────────────────────────────────────────────────

# 1. HW Configs
HW_BASE_DIR     = "./configs/tested-cfgs/SM100_B200_fabric"
BASE_ICNT       = os.path.join(HW_BASE_DIR, "config_blackwell_islip.icnt")
BASE_GPGPUSIM   = os.path.join(HW_BASE_DIR, "gpgpusim.config")
HW_OUT_BASE_DIR = "./configs/tested-cfgs"

# 2. Trace Configs
TRACE_SRC_DIR      = "../configs/tested-cfgs/SM100_B200"
TRACE_OUT_BASE_DIR = "../configs/tested-cfgs"

# 3. YAML
YAML_DEF_FILE = "../../util/job_launching/configs/define-standard-cfgs.yml"

# ── CLI ───────────────────────────────────────────────────────────────────────
def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--structure", nargs="+", choices=list(STRUCTURES), metavar="STRUCT")
    ap.add_argument("--bandwidth", nargs="+", choices=list(BANDWIDTHS), metavar="BW")
    ap.add_argument("--routing", nargs="+", choices=ROUTING_CHOICES, metavar="ROUTING", default=["baseline"])
    ap.add_argument("--near-min-k", nargs="+", type=int, metavar="K", default=[2])
    ap.add_argument("--near-min-p", nargs="+", type=float, metavar="P", default=[1.0])
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if args.all:
        structures = list(STRUCTURES)
        bandwidths = list(BANDWIDTHS)
    else:
        if not args.structure or not args.bandwidth:
            ap.error("Specify --structure and --bandwidth, or use --all")
        structures = args.structure
        bandwidths = args.bandwidth

    routing_keys: list[str] = []
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

    # Remove any duplicates (especially from near_min_random expanding multiple 'p' values)
    routing_keys = list(dict.fromkeys(routing_keys))

    combos = [(s, b, rk) for s in structures for b in bandwidths for rk in routing_keys]

    print(f"Generating {len(combos)} config pair(s)...\n")
    
    yaml_payload_lines = []

    for struct, bw, rk in combos:
        s_info = STRUCTURES[struct]
        b_info = BANDWIDTHS[bw]
        
        dir_name = f"SM100_{struct}_{bw}_{rk}"
        yaml_key = f"{struct}_{bw}_{rk}" 
        
        K = s_info["num_xbars"] * s_info["hbm_per_side"] * 2
        n_clusters = s_info["num_xbars"] * s_info["sm_per_xbar"]
        n_mem = b_info["l2_per_hbm"] * K // N_SUB_PARTITION

        hw_out_dir = os.path.join(HW_OUT_BASE_DIR, dir_name)
        icnt_dst   = os.path.join(hw_out_dir, "config_blackwell_islip.icnt")
        gpgpu_dst  = os.path.join(hw_out_dir, "gpgpusim.config")

        trace_out_dir = os.path.join(TRACE_OUT_BASE_DIR, dir_name)

        yaml_path = f"$GPGPUSIM_ROOT/configs/tested-cfgs/{dir_name}/gpgpusim.config"
        
        yaml_payload_lines.append((yaml_key, yaml_path))

        if args.dry_run:
            continue
        
        os.makedirs(hw_out_dir, exist_ok=True)
        os.makedirs(trace_out_dir, exist_ok=True)

        if os.path.exists(TRACE_SRC_DIR):
            for item in os.listdir(TRACE_SRC_DIR):
                src_file = os.path.join(TRACE_SRC_DIR, item)
                dst_file = os.path.join(trace_out_dir, item)
                if os.path.isfile(src_file):
                    shutil.copy2(src_file, dst_file)
        else:
            print(f"[ERROR] Trace source dir not found: {TRACE_SRC_DIR}")

        if os.path.exists(BASE_ICNT) and os.path.exists(BASE_GPGPUSIM):
            ov = _icnt_overrides(struct, bw, rk)
            _patch_file(BASE_ICNT, icnt_dst, ov)
            _patch_gpgpusim(BASE_GPGPUSIM, gpgpu_dst, n_mem, n_clusters)
        else:
            print(f"[ERROR] HW Base configs not found in {HW_BASE_DIR}")

        print(f"[OK] {dir_name} (Clusters: {n_clusters})")
        print(f"  ├─ HW: {hw_out_dir}")
        print(f"  └─ Trace: {trace_out_dir}")
        print("-" * 60)

    if args.dry_run:
        print("[DRY RUN] Finished without writing files.")
        return

    if os.path.exists(YAML_DEF_FILE):
        print(f"Updating YAML: {YAML_DEF_FILE}")
        with open(YAML_DEF_FILE, "r") as f:
            yaml_lines = f.readlines()
            
        existing_content = "".join(yaml_lines)
        
        insert_str = ""
        appended_count = 0
        for key_name, path in yaml_payload_lines:
            if f"{key_name}:" not in existing_content:
                insert_str += f"{key_name}:\n"
                insert_str += f"    base_file: \"{path}\"\n\n"
                appended_count += 1
                
        if insert_str:
            insert_idx = -1
            for i, line in enumerate(yaml_lines):
                if line.strip() == "B200:":
                    insert_idx = i + 2 
                    break
            
            if insert_idx != -1:
                yaml_lines.insert(insert_idx, insert_str)
                with open(YAML_DEF_FILE, "w") as f:
                    f.writelines(yaml_lines)
                print(f"SUCCESS: Inserted {appended_count} new entries right below 'B200:' block.")
            else:
                print("WARNING: Could not find 'B200:' in yaml. Appending to bottom instead.")
                with open(YAML_DEF_FILE, "a") as f:
                    f.write("\n" + insert_str)
        else:
            print("INFO: All YAML entries already exist. No changes made.")
    else:
        print(f"[ERROR] YAML file not found: {YAML_DEF_FILE}")

if __name__ == "__main__":
    main()