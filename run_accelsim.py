#!/usr/bin/env python3
"""
Run AccelSim trace-based simulations for HBMNet experiments.

For each (structure × bandwidth × routing) combination:
  1. Reads the hit rate for the given benchmark from hitrate.csv
  2. Patches baseline_ratio in the generated .icnt config
  3. Runs run_simulations.py sequentially (one at a time)
  4. Retries once on failure

Usage example:
  python run_accelsim.py \\
      --benchmark rodinia-3.1:bfs-rodinia-3.1 \\
      --structure B100_Global \\
      --bandwidth B200+HBM3e \\
      --routing baseline min_adaptive
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from typing import Optional

from experiments_loader import load_experiments

# ── Paths ─────────────────────────────────────────────────────────────────────
_HERE          = os.path.dirname(os.path.abspath(__file__))
HITRATE_CSV    = os.path.join(_HERE, "hitrate.csv")
HW_OUT_BASE    = os.path.join(_HERE, "configs", "tested-cfgs")
RUN_SIM_SCRIPT = os.path.join(_HERE, "..", "..", "util", "job_launching", "run_simulations.py")
TRACE_BASE     = os.path.join(_HERE, "hw_run", "traces", "device-0", "12.8")

# ── Experiment tables ─────────────────────────────────────────────────────────
STRUCTURES, BANDWIDTHS = load_experiments(os.path.join(_HERE, "experiments.csv"))

ROUTING_CHOICES = [
    "baseline", "min_oblivious", "min_adaptive",
    "near_min_adaptive", "ugal", "valiant",
]

# ── Abbreviation tables (for -N job name) ─────────────────────────────────────
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


# ── Helpers ───────────────────────────────────────────────────────────────────

def routing_to_key(routing: str, near_min_p: Optional[float] = None) -> str:
    if routing == "near_min_adaptive":
        p = near_min_p if near_min_p is not None else 1.0
        return f"near_min_p{p:.1f}"
    return routing


def route_abbrev(rk: str) -> str:
    if rk in ROUTING_ABBREV:
        return ROUTING_ABBREV[rk]
    m = re.match(r"near_min_p([\d.]+)", rk)
    if m:
        return "nm" + m.group(1).replace(".", "")
    return rk[:4]


def job_name(struct: str, bw: str, rk: str, bench: str) -> str:
    """Build a short job name for -N: {struct_abbrev}{bw_abbrev}{route_abbrev}_{app}."""
    sa = STRUCT_ABBREV.get(struct, struct[:3])
    ba = BW_ABBREV.get(bw, bw[:3])
    ra = route_abbrev(rk)
    # bench is like "rodinia-3.1:bfs-rodinia-3.1" — use part after ':'
    app = bench.split(":")[-1] if ":" in bench else bench
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
    """Run cmd and return its exit code."""
    print(f"\n[RUN] {' '.join(cmd)}")
    result = subprocess.run(cmd)
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
        dir_name = f"SM100_{struct}_{bw}_{rk}"
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
        cfg_name = f"{dir_name}-SASS"
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
