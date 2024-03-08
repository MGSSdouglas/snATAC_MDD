import pdb
import os 
import sys 
import pickle
import os.path as osp

import pandas as pd 
import numpy as np 
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import statsmodels as sm
import pyranges as pr

from Bio import SeqIO
from statsmodels.graphics.gofplots import qqplot_2samples
from tqdm import tqdm
from copy import deepcopy
from pprint import pprint


"""
Configuration
"""

# temporary filepaths
temp_dir = "./temp"
temp_data_dir = f"{temp_dir}/data"
os.makedirs(temp_data_dir, exist_ok=True)

# output filepaths
output_dir = "./output"
output_data_dir = f"{output_dir}/data"
output_plot_dir = f"{output_dir}/plots"
os.makedirs(output_data_dir, exist_ok=True)

null_sequences_count = 10

def aggregate_v1(x):
    """This function is used to aggregate rows of pandas dataframe
    """
    if len(np.unique(x)) != 1:
        return ",".join([str(y) for y in x])
    else:
        return [str(y) for y in x][0]

def aggregate_v2(x):
    """This function is used to aggregate rows of pandas dataframe
    """
    return ",".join([str(y) for y in x])

#############################################
## dbSNP155 (UCSC) and UCSC toolset preparation
#############################################

dbSNP155_path = "/home/dcakma3/scratch/UCSC_genome_browser_tracks/dbSnp155.bb"
'''dbSNP155 bigBed information
table bigDbSnp
"Variant summary data extracted from dbSNP, 2019 and later"
    (
    string   chrom;            "Reference sequence chromosome or scaffold"
    uint     chromStart;       "Start position in chrom"
    uint     chromEnd;         "End position in chrom"
    string   name;             "dbSNP Reference SNP (rs) identifier"
    lstring  ref;              "Reference allele; usually major allele, but may be minor allele"
    int      altCount;         "Number of alternate alleles (usually 1)"
    lstring[altCount]  alts;   "Alternate (non-reference) alleles; may include major allele"
    uint     shiftBases;       "Bases by which an indel placement could be shifted left or right (display shows thin line over uncertain region, shifts minimal representation right)"
    int      freqSourceCount;  "Number of projects reporting frequencies in current dbSNP build"
    double[freqSourceCount] minorAlleleFreq;  "Minor allele frequency, i.e. second highest allele frequency, from each frequency source; NaN if no data from project"
    lstring[freqSourceCount] majorAllele;     "Allele most frequently observed by each source"
    lstring[freqSourceCount] minorAllele;     "Allele second most frequently observed by each source"
    uint     maxFuncImpact;    "Sequence Ontology (SO) ID number for greatest functional impact on gene; 0 if no SO terms are annotated"
    enum(snv, mnv, ins, del, delins, identity) class;  "Variation class/type"
    lstring  ucscNotes;        "Interesting or anomalous properties noted by UCSC"
    bigint _dataOffset;        "Offset into bigDbSnpDetails file for line with more info"
    int _dataLen;              "Length of line in bigDbSnpDetails"
    )
'''

dbSNP155Details_path = "/home/dcakma3/scratch/UCSC_genome_browser_tracks/dbSnp155Details.tab.gz"
'''dbSNP155Details table information
table dbSnpDetails
"dbSNP annotations that are too lengthy to include in bigDbSnp file; for variant details page"
    (
    string    name;            "rs# ID of variant"
    int freqSourceCount;       "Number of frequency sources"
    string[freqSourceCount] alleleCounts;  "Array of each source's |-sep list of allele:count"
    int[freqSourceCount] alleleTotals; "Array of each source's total number of chromosomes sampled; may be > sum of observed counts and differ across variants."
    int soTermCount;           "Number of distinct SO terms annotated on variant"
    int[soTermCount] soTerms;  "SO term numbers annotated on RefSeq transcripts"
    int clinVarCount;          "Number of ClinVar accessions associated with variant"
    string[clinVarCount] clinVarAccs;  "ClinVar accessions associated with variant"
    string[clinVarCount] clinVarSigs;  "ClinVar significance for each accession"
    int submitterCount;                "Number of organizations/projects that reported variant"
    string[submitterCount] submitters; "dbSNP 'handles' of submitting orgs/projects"
    int pubMedIdCount;         "Number of PubMed-indexed publications associated with variant"
    int [pubMedIdCount]pubMedIds;      "PMIDs of associated publications"
    )
'''

# # load tab delimited dbSNPDetails file
# dbSNP155Details = pd.read_csv(dbSNP155Details_path, sep="\t", compression="gzip")


'''for dbSNP155Details and dbSNP155.bb file
For columns that contain lists of allele frequency data, the order of projects providing the data listed is as follows:

1000Genomes
dbGaP_PopFreq
TOPMED
KOREAN
SGDP_PRJ
Qatari
NorthernSweden
Siberian
TWINSUK
TOMMO
ALSPAC
GENOME_DK
GnomAD
GoNL
Estonian
Vietnamese
Korea1K
HapMap
PRJEB36033
HGDP_Stanford
Daghestan
PAGE_STUDY
Chileans
MGP
PRJEB37584
GoESP
ExAC
GnomAD_exomes
FINRISK
PharmGKB
PRJEB37766
'''

bigBedNamedItems_path = "/home/dcakma3/scratch/ucsc_cmdline_utils/bigBedNamedItems"
'''Run as follows
bigBedNamedItems -nameFile dbSnp155.bb myIds.txt dbSnp155.myIds.bed
'''

# load pipeline data structure
with open(osp.join(temp_data_dir, "merged_data.v4.pkl"), "rb") as f:
    pipeline_data = pickle.load(f)

# extract unique rsids from pipeline data structures
rsids = []
for cell_type in pipeline_data:
    snp_table = pipeline_data[cell_type]["original_snp_tables"]
    rsids += snp_table["Name"].tolist()
rsids = np.unique(rsids)
# save them to a textfile
rsid_save_path = f"{temp_data_dir}/merged_data.v4.unique.rsid"
with open(rsid_save_path, "w") as f:
    for rsid in rsids:
        f.write(f"{rsid}\n")
        
## run bigBedNamedItems on statistically significant SNP rsids
cmd = []
cmd.append(bigBedNamedItems_path)
cmd.append("-nameFile")
cmd.append(dbSNP155_path)
cmd.append(rsid_save_path)
cmd.append(f"{rsid_save_path}.dbSNP155.bed")
print(" ".join(cmd))


# result = subprocess.run(cmd, shell=True, capture_output=True) 
# TODO: above line does not work. run the printed command manually in the terminal

# dbSNP155 results 
dbsnp_results = pd.read_csv(f"{rsid_save_path}.dbSNP155.bed", sep="\t", header=None)
dbsnp_results.columns = ["chrom", "chromStart", "chromEnd", "name", "ref", "altCount", "alts", "shiftBases", "freqSourceCount", "minorAlleleFreq", "majorAllele", "minorAllele", "maxFuncImpact", "class", "ucscNotes", "_dataOffset", "_dataLen"]
# only retain chr1, chr2 to chr22, chrX, chrY
dbsnp_results = dbsnp_results[dbsnp_results.chrom.isin(
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
)]
# add "dbSNP155__" prefix to column names
dbsnp_results.columns = ["dbSNP155__" + col for col in dbsnp_results.columns]
# rename "dbSNP155__name" to "Name"
dbsnp_results = dbsnp_results.rename(columns={"dbSNP155__name": "Name"})
# set index to Name
dbsnp_results.set_index("Name", inplace=True)

# join pipeline and dbsnp_results table
for cell_type in pipeline_data:
    snp_table = pipeline_data[cell_type]["original_snp_tables"].set_index("Name")
    snp_table = snp_table.join(dbsnp_results, how="left").reset_index()
    pipeline_data[cell_type]["original_snp_tables"] = snp_table

# save the resulting pipeline_data structure 
with open(f"{temp_data_dir}/merged_data.v5.pkl", "wb") as f:
    pickle.dump(pipeline_data, f)


print("Script finished")
