#!/bin/bash

# Job script
echo "Job started at: `date`"
echo "Job ID: $SLURM_JOBID"

source "${GKMSVM_ENV_DIR}/p38/bin/activate"

python 1__preprocessing.py \
    --cluster $1 --ocr-dir $2 \
    --candidate-negative-file $3 \
    --target-dir $4 --genome-path $5 \
    --fold-split-file $6 --num-gc-bins $7 \

echo "Script finished"
