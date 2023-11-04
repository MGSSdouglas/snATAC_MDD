#!/bin/bash

# Job script
echo "Job started at: `date`"
echo "Job ID: $SLURM_JOBID"

cluster=$1
fold_id=$2
data_dir=$3
predict_dir=$4

gkmtrain="${5}/gkmtrain"
gkmpredict="${5}/gkmpredict"

models_dir=$6

./1__train.sh $cluster $fold_id $data_dir \
    $gkmtrain $models_dir

./2__eval.sh $cluster $fold_id $data_dir $predict_dir \
    $gkmpredict $models_dir

echo "Script finished"
