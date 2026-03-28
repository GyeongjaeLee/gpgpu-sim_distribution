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
    "near_min_adaptive", "ugal", "valiant",
]

def routing_to_key(routing: str, near_min_p: Optional[float] = None) -> str:
    if routing == "near_min_adaptive":
        p = near_min_p if near_min_p is not None else 1.0
        return f"near_min_p{p:.1f}"
    return routing

def routing_key_to_overrides(key: str) -> dict:
    if key == "baseline":
        return {"routing_function": "baseline", "is_fabric": "0"}
    ov: dict = {"routing_function": "hybrid", "is_fabric": "1"}
    m = re.match(r"near_min_p([\d.]+)$", key)
    if m:
        ov["hybrid_routing"]   = "near_min_adaptive"
        ov["near_min_penalty"] = m.group(1)
        return ov
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

# 1. HW Configs (생성 위치: 현재 경로 내 ./configs/...)
HW_BASE_DIR     = "./configs/tested-cfgs/SM100_B200_fabric"
BASE_ICNT       = os.path.join(HW_BASE_DIR, "config_blackwell_islip.icnt")
BASE_GPGPUSIM   = os.path.join(HW_BASE_DIR, "gpgpusim.config")
HW_OUT_BASE_DIR = "./configs/tested-cfgs"

# 2. Trace Configs (명시적 상대 경로: ../configs/...)
TRACE_SRC_DIR      = "../configs/tested-cfgs/SM100_B200"
TRACE_OUT_BASE_DIR = "../configs/tested-cfgs"

# 3. YAML 정의 파일 (명시적 상대 경로)
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
            for p in args.near_min_p:
                routing_keys.append(routing_to_key(r, p))
        else:
            routing_keys.append(r)

    combos = [(s, b, rk) for s in structures for b in bandwidths for rk in routing_keys]

    print(f"Generating {len(combos)} config pair(s)...\n")
    
    yaml_payload_lines = []

    for struct, bw, rk in combos:
        s_info = STRUCTURES[struct]
        
        # 실제 디렉토리명
        dir_name = f"SM100_{struct}_{bw}_{rk}"
        # YAML에 들어갈 키 이름 (B200: 처럼 앞부분 생략)
        yaml_key = f"{struct}_{bw}_{rk}" 
        
        # --- 새로 추가된 부분: 클러스터(SM) 수 동적 계산 ---
        # Structure의 num_xbars 값에 Xbar 당 SM 수를 곱합니다.
        b_info = BANDWIDTHS[bw]
        K = s_info["num_xbars"] * s_info["hbm_per_side"] * 2
        n_clusters = s_info["num_xbars"] * s_info["sm_per_xbar"]
        n_mem = b_info["l2_per_hbm"] * K // N_SUB_PARTITION

        # 1. HW Config 아웃풋 경로
        hw_out_dir = os.path.join(HW_OUT_BASE_DIR, dir_name)
        icnt_dst   = os.path.join(hw_out_dir, "config_blackwell_islip.icnt")
        gpgpu_dst  = os.path.join(hw_out_dir, "gpgpusim.config")

        # 2. Trace Config 아웃풋 경로
        trace_out_dir = os.path.join(TRACE_OUT_BASE_DIR, dir_name)

        # 3. YAML에 기록될 포맷 ($GPGPUSIM_ROOT 기준)
        yaml_path = f"$GPGPUSIM_ROOT/configs/tested-cfgs/{dir_name}/gpgpusim.config"
        
        # 페이로드 임시 저장 (중복 방지는 나중에 처리)
        yaml_payload_lines.append((yaml_key, yaml_path))

        if args.dry_run:
            continue
        
        # --- 디렉토리 생성 ---
        os.makedirs(hw_out_dir, exist_ok=True)
        os.makedirs(trace_out_dir, exist_ok=True)

        # --- Trace Config 복사 ---
        if os.path.exists(TRACE_SRC_DIR):
            for item in os.listdir(TRACE_SRC_DIR):
                src_file = os.path.join(TRACE_SRC_DIR, item)
                dst_file = os.path.join(trace_out_dir, item)
                if os.path.isfile(src_file):
                    shutil.copy2(src_file, dst_file)
        else:
            print(f"[ERROR] Trace source dir not found: {TRACE_SRC_DIR}")

        # --- HW Config 패치 및 생성 ---
        if os.path.exists(BASE_ICNT) and os.path.exists(BASE_GPGPUSIM):
            ov = _icnt_overrides(struct, bw, rk)
            _patch_file(BASE_ICNT, icnt_dst, ov)
            # n_clusters 인자를 추가로 넘겨주어 gpgpusim.config를 수정합니다.
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

    # --- YAML 파일 특정 위치(B200 바로 아래)에 삽입 로직 ---
    if os.path.exists(YAML_DEF_FILE):
        print(f"Updating YAML: {YAML_DEF_FILE}")
        with open(YAML_DEF_FILE, "r") as f:
            yaml_lines = f.readlines()
            
        existing_content = "".join(yaml_lines)
        
        # 삽입할 내용 문자열 생성 (중복 배제)
        insert_str = ""
        appended_count = 0
        for key_name, path in yaml_payload_lines:
            if f"{key_name}:" not in existing_content:
                insert_str += f"{key_name}:\n"
                insert_str += f"    base_file: \"{path}\"\n\n"
                appended_count += 1
                
        if insert_str:
            # B200: 위치 찾기
            insert_idx = -1
            for i, line in enumerate(yaml_lines):
                if line.strip() == "B200:":
                    # 그 다음 줄이 base_file: ... 일 테니, 그 다음다음 줄(i+2)에 넣기 위함
                    insert_idx = i + 2 
                    break
            
            if insert_idx != -1:
                # 찾은 위치 바로 아래에 삽입
                yaml_lines.insert(insert_idx, insert_str)
                with open(YAML_DEF_FILE, "w") as f:
                    f.writelines(yaml_lines)
                print(f"SUCCESS: Inserted {appended_count} new entries right below 'B200:' block.")
            else:
                # 혹시라도 B200: 을 못 찾으면 맨 밑에 추가
                print("WARNING: Could not find 'B200:' in yaml. Appending to bottom instead.")
                with open(YAML_DEF_FILE, "a") as f:
                    f.write("\n" + insert_str)
        else:
            print("INFO: All YAML entries already exist. No changes made.")
    else:
        print(f"[ERROR] YAML file not found: {YAML_DEF_FILE}")

if __name__ == "__main__":
    main()