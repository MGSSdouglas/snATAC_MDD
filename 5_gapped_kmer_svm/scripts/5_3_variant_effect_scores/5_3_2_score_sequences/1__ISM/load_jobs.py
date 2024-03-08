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
    cmd += ["--time 0:40:0"]

    # task specification
    cmd += ["--ntasks 1"]

    # core specification
    cmd += ["--cpus-per-task 16"]

    # memory specification per cpu
    cmd += ["--mem-per-cpu 8GB"]

    # call job script 
    cmd += ["./1_1__gkmpredict.sh"]
    cmd += [str(x) for x in [job["cluster"], job["fold_id"], job["model_path"], \
        job["model_scripts_dir"], job["cdsnp_dir"], job["cdsnp_score_dir"]]]
    
    # submit job
    print(' '.join(cmd))

    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
    print(result.stdout)

print("Script complete.")