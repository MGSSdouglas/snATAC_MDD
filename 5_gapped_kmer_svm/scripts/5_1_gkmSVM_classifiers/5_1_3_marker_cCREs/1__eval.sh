#!/bin/bash

cluster=$1
fold_id=$2
data_dir="$3/$cluster/$fold_id"
predict_dir="$4/$cluster/$fold_id"
gkmpredict=$5
model_path="$6/$cluster/$fold_id/$cluster.$fold_id.model.txt"

test_pos="$data_dir/marker.test.pos.fa"
test_neg="$data_dir/marker.test.neg.fa"

test_pos_pred="$predict_dir/marker.test.pos.pred.txt"
test_neg_pred="$predict_dir/marker.test.neg.pred.txt"

$gkmpredict $test_pos $model_load_path $test_pos_pred
$gkmpredict $test_neg $model_load_path $test_neg_pred

echo "Eval finished"