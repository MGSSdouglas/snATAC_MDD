#!/bin/bash

echo "chr,start,end,concordance" > imputation_metrics_v2.csv
echo "chr,start,end" > empty_intervals_reference_v2.csv
echo "chr,start,end" > empty_intervals_study_v2.csv

for log in ./IMPUTE2_logs/*.log; do
	chr=$(echo $log | egrep -o 'chr[0-9]{1,2}')
	start_pos=$(echo $log | egrep -o '_[0-9]*_' | egrep -o '[0-9]*')
	end_pos=$(echo $log | egrep -o '_[0-9]*\.' | egrep -o '[0-9]*')
	if grep -q "There are no SNPs in the imputation interval" $log; then
		echo "$chr,$start_pos,$end_pos" >> empty_intervals_reference_v2.csv	
	elif grep -q "ERROR: There are no type 2 SNPs after applying the command-line settings for this run" $log; then
		echo "$chr,$start_pos,$end_pos" >> empty_intervals_study_v2.csv
	else
		concord_line=$(grep "\[ >= 0.0\]" $log)
		concordance=${concord_line##* }
		echo "$chr,$start_pos,$end_pos,$concordance" >> imputation_metrics_v2.csv
	fi
done
