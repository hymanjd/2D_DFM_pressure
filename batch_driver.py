import itertools
import shutil
import os
import subprocess
import pandas as pd
import numpy as np
from scipy.stats.qmc import LatinHypercube
import random
import multiprocessing as mp
import time 

mp.set_start_method("fork")

def run_dfnworks(sample_index):
    """
    Run the DFN→Pflotran workflow for one sample, timing the run.
    """
    home = os.getcwd()
    start_time = time.time()
    
    cmd = f"python3.11 dfm_driver.py {sample_index} > x{sample_index:02d}.out"
    print(f">> {cmd}")
    
    try:
        subprocess.call(cmd, shell=True)
        os.chdir(home)
        
        jobname = f"pressure_x{sample_index:02d}"
        shutil.rmtree(jobname)
        # Clean up
        os.remove(f'pressure_x{sample_index:02d}.log')
        try:
            os.remove(f'x{sample_index:02d}.out')
        except FileNotFoundError:
            pass
        
        elapsed = time.time() - start_time
        print(f" index: {sample_index} done in {elapsed:.1f}s")
    except Exception as e:
        elapsed = time.time() - start_time
        print(f" index: {sample_index} failed after {elapsed:.1f}s")
        print(e)

 
## MAIN ##

########### LHS SAMPLING ###########
# Num of samples should fill the space appropriately.
num_of_experiments = 2
num_jobs = 2

print(f"-> running {num_of_experiments} samples")

dirs = ['mas_files', 'dfn_info', 'pflotran_files', 'h5_files']
for d in dirs:
    try:
        os.makedirs(d, exist_ok=True)
    except Exception as e:
        print(e)
        continue

job_id = list(range(1, num_of_experiments + 1))
pool = mp.Pool(num_jobs)
pool.map(run_dfnworks, job_id, chunksize=1)
pool.close()
pool.join()
pool.terminate()
