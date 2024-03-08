import os 
import errno
import subprocess
import sys 
import argparse
import os.path as osp
import pdb

from config import job_params

def makedir(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

# create jobs
jobs = []
for param in job_params:
    jobs.append({
        "Name": f"{param['cell_type']}_{param['fold_id']}",
        "cell_type": param["cell_type"],
        "fold_id": param["fold_id"],
        "model_base_path": param["model_base_path"],
        "gkmsvm_base_path": param["gkmsvm_base_path"],
        "deltasvm_base_path": param["deltasvm_output_base_path"],
        "gkmpredict_base_path": param["gkmpredict_base_path"],
        "snp_sequence_base_path": param["snp_sequence_base_path"],
        "deltasvm_script_path": param["deltasvm_script_path"],
        "broad_cluster": param["broad_cluster"]
    })

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
    cmd += ["--time 0:20:0"]

    # task specification
    cmd += ["--ntasks 1"]

    # core specification
    cmd += ["--cpus-per-task 2"]

    # memory specification per cpu
    cmd += ["--mem-per-cpu 6GB"]

    # call job script 
    cmd += ["./generate_deltaSVM_scores.sh"]
    cmd += [str(x) for x in [job["cell_type"], job["fold_id"], job["model_base_path"], \
        job["deltasvm_base_path"], job["gkmpredict_base_path"], \
        job["snp_sequence_base_path"], job["gkmsvm_base_path"], \
        job["deltasvm_script_path"], job["broad_cluster"]]]

    # submit job
    print(' '.join(cmd))
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
    print(result.stdout)

print("Script complete.")