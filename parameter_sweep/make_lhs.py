"""
batch_driver.py — LHS + structured diagnostic sweep for DFM Darcy-flow data.

Usage
-----
    python3.11 batch_driver.py [sweep_config.yaml]

Workflow
--------
1. Read sweep_config.yaml (default if no arg given).
2. Generate LHS block in (k_m, k_ratio, p32) space; derive k_f = k_m * k_ratio.
3. Append structured diagnostic block (factorial k_ratio × p32 at fixed k_m).
4. Write combined lhs_params.csv with columns:
       domain_index, k_f, k_m, k_ratio, p32, sample_type
5. Create working directories.
6. Run dfm_driver.py in parallel over all domain indices.
"""

import os
import sys
import shutil
import subprocess
import time
import multiprocessing as mp
import itertools

import numpy as np
import pandas as pd
import yaml
from scipy.stats.qmc import LatinHypercube

mp.set_start_method("fork")


# ── Config ────────────────────────────────────────────────────────────────────

def load_config(path="sweep_config.yaml"):
    with open(path) as f:
        return yaml.safe_load(f)


# ── Table generation ──────────────────────────────────────────────────────────

def generate_lhs_block(cfg):
    """Draw LHS samples in (k_m, k_ratio, p32) space.

    Returns a DataFrame with columns: k_f, k_m, k_ratio, p32, sample_type.
    """
    sweep  = cfg["sweep"]
    params = cfg["parameters"]
    n      = sweep["num_of_experiments"]
    seed   = sweep.get("seed", None)

    param_names = ["k_m", "k_ratio", "p32"]
    sampler     = LatinHypercube(d=len(param_names), seed=seed)
    unit        = sampler.random(n=n)   # (n, 3) in [0, 1)

    rows = {}
    for col, name in enumerate(param_names):
        lo = params[name]["log10_min"]
        hi = params[name]["log10_max"]
        rows[name] = 10.0 ** (lo + (hi - lo) * unit[:, col])

    df = pd.DataFrame(rows)
    df["k_f"]         = df["k_m"] * df["k_ratio"]
    df["sample_type"] = "lhs"
    return df[["k_f", "k_m", "k_ratio", "p32", "sample_type"]]


def generate_diagnostic_block(cfg):
    """Build a factorial grid over k_ratio × p32 at fixed k_m.

    Returns a DataFrame with the same columns as the LHS block.
    """
    diag = cfg.get("diagnostic", {})
    if not diag.get("enabled", False):
        return pd.DataFrame(columns=["k_f", "k_m", "k_ratio", "p32", "sample_type"])

    k_m_fixed     = float(diag["k_m_fixed"])
    k_ratio_levels = [float(v) for v in diag["k_ratio_levels"]]
    p32_levels     = [float(v) for v in diag["p32_levels"]]

    records = []
    for k_ratio, p32 in itertools.product(k_ratio_levels, p32_levels):
        records.append({
            "k_f":         k_m_fixed * k_ratio,
            "k_m":         k_m_fixed,
            "k_ratio":     k_ratio,
            "p32":         p32,
            "sample_type": "diagnostic",
        })
    return pd.DataFrame(records)


def build_and_save_param_table(cfg, out_csv="lhs_params.csv"):
    """Combine LHS and diagnostic blocks, assign domain indices, write CSV."""
    lhs   = generate_lhs_block(cfg)
    diag  = generate_diagnostic_block(cfg)

    combined = pd.concat([lhs, diag], ignore_index=True)
    combined.insert(0, "domain_index", np.arange(1, len(combined) + 1))
    combined.to_csv(out_csv, index=False)

    n_lhs  = len(lhs)
    n_diag = len(diag)
    print(f"--> Parameter table: {n_lhs} LHS  +  {n_diag} diagnostic  =  {len(combined)} total")
    print(f"--> Written to {out_csv}")
    return combined


# ── Per-sample runner ─────────────────────────────────────────────────────────

def run_dfnworks(args):
    sample_index, cfg = args
    home    = os.getcwd()
    cleanup = cfg.get("cleanup", {})
    start   = time.time()

    stdout_file = f"x{sample_index:04d}.out"
    cmd = f"python3.11 dfm_driver.py {sample_index} > {stdout_file} 2>&1"
    print(f">> [{sample_index:04d}] launching")

    try:
        subprocess.call(cmd, shell=True)
        os.chdir(home)

        if cleanup.get("remove_jobdir", True):
            jobdir = f"pressure_x{sample_index:02d}"
            if os.path.isdir(jobdir):
                shutil.rmtree(jobdir)

        if cleanup.get("remove_logfile", True):
            logfile = f"pressure_x{sample_index:02d}.log"
            if os.path.exists(logfile):
                os.remove(logfile)

        if cleanup.get("remove_stdout", True) and os.path.exists(stdout_file):
            os.remove(stdout_file)

        print(f"   [{sample_index:04d}] done in {time.time() - start:.1f}s")

    except Exception as e:
        print(f"   [{sample_index:04d}] FAILED after {time.time() - start:.1f}s — {e}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    config_path = sys.argv[1] if len(sys.argv) > 1 else "sweep_config.yaml"
    cfg = load_config(config_path)

    sweep    = cfg["sweep"]
    num_jobs = sweep["num_jobs"]

    print(f"--> Config: {config_path}")

    for d in cfg.get("directories", []):
        os.makedirs(d, exist_ok=True)

    table    = build_and_save_param_table(cfg)
    n_total  = len(table)
#    print(f"--> Running {n_total} samples on {num_jobs} workers\n")
#
#    args = [(int(idx), cfg) for idx in table["domain_index"]]
#    pool = mp.Pool(num_jobs)
#    pool.map(run_dfnworks, args, chunksize=1)
#    pool.close()
#    pool.join()
#    pool.terminate()
#
    print(f"\n--> Sweep complete ({n_total} samples).")
