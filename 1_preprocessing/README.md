## Main
This directory contains the scripts developed for demultiplexing male and female combined libraries and genotyping external data

## Schematic
<img width="1184" alt="Screenshot 2023-10-20 at 3 19 56 PM" src="https://github.com/MGSSdouglas/snATAC_MDD/assets/60046859/6bbd9fef-c1d7-4055-9eea-dbc6f1a112ef">

Schematic overview of the snATAC-seq study design in 84 subjects from nuclei extraction and multiplexing to downstream data analyses

## We created snATAC-seq objects from 10x fragment files 
[scripts](https://github.com/MGSSdouglas/snATAC_MDD/tree/main/1_preprocessing/snATAC_preprocessing)

## To demultiplex snATAC-seq male and female cells 
[scripts](https://github.com/MGSSdouglas/snATAC_MDD/tree/main/1_preprocessing/snATAC_preprocessing/07_demultiplex_barcodes_by_variants.sh)

## To preprocess external genotypes from genotyping arrays
[scripts](https://github.com/MGSSdouglas/snATAC_MDD/tree/main/1_preprocessing/genotyping_qc_and_preprocessing)

## To match demultiplexed subjects with external genotypes
[scripts](https://github.com/MGSSdouglas/snATAC_MDD/tree/main/1_preprocessing/genotyping_qc_and_preprocessing/snATAC_genotype_comparison)

## To plot correlations between snATAC-seq and genotypes
[scripts](https://github.com/MGSSdouglas/snATAC_MDD/blob/main/1_preprocessing/snATAC_preprocessing/08_plot_genotype_correlations.R)

## To perform snATAC-seq quality control
[scripts](https://github.com/MGSSdouglas/snATAC_MDD/tree/main/1_preprocessing/quality_control)
