import os 
import pdb
import sys 
import errno
import argparse
import subprocess

import os.path as osp

from config import jobs

def makedir(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

# sanity check
if len(jobs) > 500:
    print("Cannot submit more than 500 jobs simultaneously. Please re-run with specific arguments.")
    sys.exit()

# submit jobs
for job in jobs:

    # create log file
    makedir(osp.dirname(f"./log/{job['Name']}.out"))
    try:
        os.remove(f"./log/{job['Name']}.out")
    except Exception as e:
        pass

    # create job queue command
    cmd = ["sbatch", "-J", job['Name']]

    # time specification
    cmd += ["--time 0:15:0"]

    # task specification
    cmd += ["--ntasks 1"]

    # core specification
    cmd += ["--cpus-per-task 2"]

    # memory specification per cpu
    cmd += ["--mem-per-cpu 4GB"]

    # call job  
    cmd += ["./2_1__1000G_LD_expansion.sh", str(job["chrom"]), job["target_dir"], 
        job["bimfile_path"], job["rsids_dir"]]

    #submit job
    print(' '.join(cmd))
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
    print(result.stdout)


print("Script finished")
