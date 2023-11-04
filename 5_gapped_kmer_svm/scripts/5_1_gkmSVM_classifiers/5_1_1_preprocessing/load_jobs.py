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
    cmd += ["--time 3:0:0"]

    # cpu specification
    cmd += ["--cpus-per-task 8"]

    # memory specification per cpu
    cmd += ["--mem-per-cpu 8GB"]

    # call job script 
    cmd += ["./init_job.sh", job["cluster"], job["OCR_dir"], 
        job["candidate_neg_path"], job["target_dir"], 
        job["genome_path"], job["fold_split_file"], 
        job["num_gc_bins"]]

    #submit job
    print(' '.join(cmd))
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
    print(result.stdout)


print("Script finished")
