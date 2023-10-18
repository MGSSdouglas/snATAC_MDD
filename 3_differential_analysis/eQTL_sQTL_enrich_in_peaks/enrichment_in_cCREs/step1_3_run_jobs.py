import os 
import errno
import subprocess
import sys 
import argparse
import os.path as osp
import pdb
from step14_1_job_input_preparation import job_params


# ###### DEBUGGING ######
# time python step14_analysis_QTL_peakset_P2G.py --peakset-name ct_M_0.05 --qtl-name eQTL --out-dir ./step14.output
# real    8m35.095s
# user    7m17.320s
# sys     1m10.920s
# specs:
# mem=${mem:-"4GB"}
# cpu=${cpu:-12}
# salloc --time="$time" --ntasks=1 --cpus-per-task="$cpu" --account=ctb-liyue --mem-per-cpu="$mem"



# config
verbose = True

# #####################
# ##### DEBUGGING #####
# #####################
# job_params = job_params[0:2]


if len(job_params) > 500:
    print("Cannot submit more than 500 jobs simultaneously. Please re-run with specific arguments.")
    sys.exit()

# submit jobs
for job in job_params:
    
    cmd = ["sbatch", "-J", job['Name']] # SLURM command
    cmd += ["--time 0:15:0"] # time specification
    cmd += ["--cpus-per-task 8"] # cpu specification
    cmd += ["--mem-per-cpu 8GB"] # memory specification per cpu

    # call job script 
    job_script_name = "./step14_2_job.sh"
    cmd += [job_script_name, job["peakset_name"], job["qtl_name"], job["out_dir"]]
    
    if verbose:
        print("Command: ", " ".join(cmd))

    # submit job
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)

    if verbose:
        print("Result:", result.stdout)
    
    

print("Script complete.")
