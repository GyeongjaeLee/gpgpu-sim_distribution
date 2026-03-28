#!/usr/bin/env python3
"""
Run AccelSim trace-based simulations for HBMNet experiments.

For each (structure × bandwidth × routing) combination:
  1. Reads the hit rate for the given benchmark from hitrate.csv
  2. Patches baseline_ratio in the generated .icnt config
     (config must already exist — run gen_accelsim_configs.py first)
  3. Runs run_simulations.py sequentially (one at a time)
  4. Retries once on failure

Job names are auto-abbreviated, e.g.:
  B100_Global + B200+HBM3e + near_min_p0.0 + bfs  →  "BGH3nm0_bfs"

Usage examples:

  # Single structure/bandwidth, multiple routings
  python run_accelsim.py \\
      --benchmark rodinia-3.1:bfs-rodinia-3.1 \\
      --structure B100_Global \\
      --bandwidth B200+HBM3e \\
      --routing baseline min_adaptive near_min_adaptive \\
      --near-min-p 0.0 1.0

  # Multiple structures and bandwidths
  python run_accelsim.py \\
      --benchmark polybench:polybench-gemm \\
      --structure B100_Local B100_Global H100 \\
      --bandwidth B200+HBM3e Shoreline_1x \\
      --routing baseline near_min_adaptive \\
      --near-min-p 0.0

  # All configs from experiments.csv (dry-run to preview)
  python run_accelsim.py \\
      --benchmark GPU_Microbenchmark:mem_bw \\
      --all-configs \\
      --routing baseline \\
      --dry-run
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from typing import Optional

from experiments_loader import (load_experiments,
                                STRUCT_ABBREV, BW_ABBREV, ROUTING_ABBREV,
                                bench_abbrev, route_abbrev)

# ── Paths ─────────────────────────────────────────────────────────────────────
_HERE          = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR       = os.path.normpath(os.path.join(_HERE, "..", ".."))  # cwd for run_simulations.py
HITRATE_CSV    = os.path.join(_HERE, "hitrate.csv")
HW_OUT_BASE    = os.path.join(_HERE, "configs", "tested-cfgs")
RUN_SIM_SCRIPT = os.path.join(ROOT_DIR, "util", "job_launching", "run_simulations.py")
TRACE_BASE     = os.path.join(ROOT_DIR, "hw_run", "traces", "device-0", "12.8")

# ── Experiment tables ─────────────────────────────────────────────────────────
STRUCTURES, BANDWIDTHS = load_experiments(os.path.join(_HERE, "experiments.csv"))

ROUTING_CHOICES = [
    "baseline", "min_oblivious", "min_adaptive",
    "near_min_adaptive", "ugal", "valiant",
]

# ── Helpers ───────────────────────────────────────────────────────────────────

def routing_to_key(routing: str, near_min_p: Optional[float] = None) -> str:
    if routing == "near_min_adaptive":
        p = near_min_p if near_min_p is not None else 1.0
        return f"near_min_p{p:.1f}"
    return routing


def job_name(struct: str, bw: str, rk: str, bench: str) -> str:
    """Build a short job name for -N: {struct_abbrev}{bw_abbrev}{route_abbrev}_{bench_abbrev}."""
    sa  = STRUCT_ABBREV.get(struct, struct[:3])
    ba  = BW_ABBREV.get(bw, bw[:3])
    ra  = route_abbrev(rk)
    app = bench_abbrev(bench)
    return f"{sa}{ba}{ra}_{app}"


def load_hitrates(path: str) -> dict:
    """Parse hitrate.csv → {benchmark: float}. Missing/empty → None."""
    rates: dict = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            bench = row["Benchmark"].strip()
            val   = row["Hit_rate"].strip()
            rates[bench] = float(val) if val else None
    return rates


def patch_baseline_ratio(icnt_path: str, ratio: float) -> None:
    """Overwrite baseline_ratio in the given .icnt file."""
    with open(icnt_path) as f:
        lines = f.readlines()
    result = []
    patched = False
    for line in lines:
        if re.match(r"^\s*baseline_ratio\s*=", line):
            result.append(f"baseline_ratio = {ratio};\n")
            patched = True
        else:
            result.append(line)
    if not patched:
        result.append(f"baseline_ratio = {ratio};\n")
    with open(icnt_path, "w") as f:
        f.writelines(result)


def run_once(cmd: list) -> int:
    """Run cmd from ROOT_DIR (../../) and return its exit code."""
    print(f"\n[RUN] (cwd={ROOT_DIR})")
    print(f"      {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=ROOT_DIR)
    return result.returncode


def run_with_retry(cmd: list, retries: int = 1) -> bool:
    """Run cmd; retry up to `retries` times on failure. Returns True if success."""
    for attempt in range(retries + 1):
        rc = run_once(cmd)
        if rc == 0:
            return True
        if attempt < retries:
            print(f"[WARN] Command failed (exit {rc}). Retrying ({attempt+1}/{retries})...")
        else:
            print(f"[ERROR] Command failed after {retries+1} attempt(s). Skipping.")
    return False


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--benchmark", required=True,
                    help="Benchmark spec, e.g. rodinia-3.1:bfs-rodinia-3.1")
    ap.add_argument("--structure", nargs="+", choices=list(STRUCTURES),
                    metavar="STRUCT", required=True)
    ap.add_argument("--bandwidth", nargs="+", choices=list(BANDWIDTHS),
                    metavar="BW", required=True)
    ap.add_argument("--routing", nargs="+", choices=ROUTING_CHOICES,
                    metavar="ROUTING", default=["baseline"])
    ap.add_argument("--near-min-p", nargs="+", type=float, metavar="P",
                    default=[1.0])
    ap.add_argument("--all-configs", action="store_true",
                    help="Use all structures and bandwidths from experiments.csv")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print actions without running or patching files")
    args = ap.parse_args()

    if args.all_configs:
        structures = list(STRUCTURES)
        bandwidths = list(BANDWIDTHS)
    else:
        structures = args.structure
        bandwidths = args.bandwidth

    # Build routing key list
    routing_keys: list = []
    for r in args.routing:
        if r == "near_min_adaptive":
            for p in args.near_min_p:
                routing_keys.append(routing_to_key(r, p))
        else:
            routing_keys.append(r)

    # Load hit rates
    hitrates = load_hitrates(HITRATE_CSV)
    bench    = args.benchmark
    hit_rate = hitrates.get(bench)
    if hit_rate is None:
        print(f"[WARN] No hit rate found for '{bench}' in {HITRATE_CSV}. "
              f"baseline_ratio will not be modified.")

    # Build combo list
    combos = [
        (s, b, rk)
        for s in structures
        for b in bandwidths
        for rk in routing_keys
    ]
    print(f"Benchmark    : {bench}")
    print(f"Hit rate     : {hit_rate}")
    print(f"Combinations : {len(combos)}\n")

    failures = []

    for struct, bw, rk in combos:
        config_name = f"{struct}_{bw}_{rk}"
        dir_name = f"SM100_{config_name}"
        icnt_path = os.path.join(HW_OUT_BASE, dir_name, "config_blackwell_islip.icnt")

        print(f"─── {dir_name}")

        # Check that generated config exists
        if not os.path.isfile(icnt_path):
            print(f"[ERROR] Config not found: {icnt_path}")
            print("  → Run gen_accelsim_configs.py first. Skipping.")
            failures.append((struct, bw, rk, "config missing"))
            continue

        # Patch baseline_ratio
        if hit_rate is not None:
            if args.dry_run:
                print(f"  [DRY] patch baseline_ratio = {hit_rate} → {icnt_path}")
            else:
                patch_baseline_ratio(icnt_path, hit_rate)
                print(f"  Patched baseline_ratio = {hit_rate}")

        # Build run_simulations.py command
        cfg_name = f"{config_name}-SASS"
        jname    = job_name(struct, bw, rk, bench)
        cmd = [
            sys.executable, RUN_SIM_SCRIPT,
            "-B", bench,
            "-C", cfg_name,
            "-T", TRACE_BASE,
            "-N", jname,
        ]

        if args.dry_run:
            print(f"  [DRY] {' '.join(cmd)}")
            continue

        ok = run_with_retry(cmd, retries=1)
        if not ok:
            failures.append((struct, bw, rk, "run_simulations failed"))

    # Summary
    print(f"\n{'='*60}")
    print(f"Done. {len(combos) - len(failures)}/{len(combos)} succeeded.")
    if failures:
        print("Failed combos:")
        for struct, bw, rk, reason in failures:
            print(f"  SM100_{struct}_{bw}_{rk}: {reason}")


if __name__ == "__main__":
    main()
