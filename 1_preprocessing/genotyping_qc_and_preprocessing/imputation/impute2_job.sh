#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH --account=rrg-gturecki
#SBATCH --time=00:30:00
#SBATCH --mem=10G
#SBATCH -o impute2_job_"$1"_"$2"_"$3".log
module load StdEnv/2020
module load impute2/2.3.2

mkdir IMPUTE2_results/$1/imputed_$2_$3

impute2 -use_prephased_g \
	-known_haps_g ../SHAPEIT_phasing/shapeit_phased_data/$1/phased_$1.haps \
	-h ../SHAPEIT_phasing/1000GP_Phase3/1000GP_Phase3_$1.hap.gz \
	-sample_g ../SHAPEIT_phasing/shapeit_phased_data/$1/phased_$1.sample \
	-l ../SHAPEIT_phasing/1000GP_Phase3/1000GP_Phase3_$1.legend.gz  \
	-m ../SHAPEIT_phasing/1000GP_Phase3/genetic_map_$1_combined_b37.txt \
	-int $2 $3 \
	-Ne 20000 \
	-seed 1 \
	-o IMPUTE2_results/$1/imputed_$2_$3/imputed_$2_$3  \
exit 0
EOT
