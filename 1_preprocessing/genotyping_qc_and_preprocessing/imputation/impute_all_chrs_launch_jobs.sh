#!/bin/bash

mkdir IMPUTE2_results

interval_size=3000000

#iterate over chromosomes in file containing chromosome sizes
while IFS=$'\t' read -r chr size; do
	mkdir IMPUTE2_results/$chr
	floor_num_intervals=$(($size / $interval_size - 1))
	for interval in $(seq 0 $floor_num_intervals); do 
		begin_pos=$(($interval_size*$interval))
		if [[ $interval -eq $floor_num_intervals ]]
		then
			end_pos=$size
		else
			end_pos=$(($begin_pos + $interval_size))
		fi
		echo "$chr $begin_pos $end_pos"
		bash impute2_job.sh $chr $begin_pos $end_pos
		sleep 3s
	done
done < hg19_chr_sizes.txt
