#!/bin/bash

cluster=$1
fold_id=$2
data_dir="$3/$cluster/$fold_id"
gkmtrain=$4
model_save_path_prefix="$5/$cluster/$fold_id/$cluster.$fold"

train_pos="$data_dir/CV.train.positives.fa"
train_neg="$data_dir/CV.train.negatives.fa"

# default hyperparameters for gkmsvm
t="4"
l="11"
k="7"
d="3"
M="50"
H="50"
c="1"
e="0.001"

# train model and save the trained model
$gkmtrain -t $t -l $l -k $k -d $d -M $M -H $H -c $c -e $e -T 16 \
    $train_pos $train_neg $model_save_path_prefix

echo "Training finished"
