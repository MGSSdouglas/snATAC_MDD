import pdb 
import os
import pickle
import subprocess

import os.path as osp
import pandas as pd 
import numpy as np 

from tqdm import tqdm


###### configuration 


# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# UCSC version dbSNP155 paths
ucsc_dbsnp1155 = f"{GKMSVM_RAW_DATA_DIR}/dbsnp_v155/dbSnp155.bb"
bigBedNamedItems = f"{GKMSVM_BIN_DIR}/ucsc_utils/bigBedNamedItems"

# configure paths
data_id = "MDD_GWAS_finemapping"
raw_data_dir = osp.join(GKMSVM_RAW_DATA_DIR, data_id, sep = "/")
processed_data_dir = osp.join(GKMSVM_PREPARED_DATA_DIR, data_id, sep = "/")
os.makedirs(processed_data_dir, exist_ok=True)

chromosomes = list(range(1, 23))

###### Parse finemapped snps per effect per chromosome

finemap_base_filepath = f"{raw_data_dir}/MDD"
rsid_mapping_basepath = f"{raw_data_dir}/MDD/ukb"
finemap_ALL_snps_filepath = f"{processed_data_dir}/finemap.all.hg19_information.tsv"

finemapped_snp_table = {
    "Chromosome": [],
    "Start": [],
    "End": [],
    "A1": [],
    "A2": [],
    "rsid": [],
    "Effect": [],
    "pip": [],
    "beta": []
}

num_snps_input_to_finemapping = 0
for chrom in tqdm(chromosomes, desc="Chromosomes"):

    # read genome location to rsid mapping file
    mapping_filepath = os.path.join(rsid_mapping_basepath, f"{chrom}.rsid")
    mapping_df = pd.read_table(mapping_filepath, header=None)
    mapping_dct = dict(zip(mapping_df[0], mapping_df[1]))

    num_snps_input_to_finemapping += len(mapping_dct)

    # read finemap table
    finemap_filepath = os.path.join(finemap_base_filepath, f"MDD_{chrom}.cs")
    table = pd.read_table(finemap_filepath)

    # parse finemap table
    pip, beta, chromosome, start, end, a1, a2, effect, rsid  = [], [], [], [], [], [], [], [], []
    for idx, row in tqdm(table.iterrows(), desc="Effects"):

        snps = row["cs"].strip("][").split(', ')
        snps = [x.strip("\'") for x in snps]

        rsid += [mapping_dct[x] for x in snps]

        for snp in snps:
            temp = snp.split(".")
            chromosome += [f"chr{temp[0]}"]
            start += [int(temp[1])]
            end += [int(temp[1])+1]
            a1 += [temp[2]]
            a2 += [temp[3]]
            effect += [f"{temp[0]}.{idx}"]

        pip += [float(x) for x in row["pip"].strip('][').split(', ')]
        beta += [float(x) for x in row["beta"].strip('][').split(', ')]
        

    # add snps to final table 
    finemapped_snp_table["Chromosome"] += chromosome
    finemapped_snp_table["Start"] += start
    finemapped_snp_table["End"] += end
    finemapped_snp_table["rsid"] += rsid
    finemapped_snp_table["A1"] += a1
    finemapped_snp_table["A2"] += a2
    finemapped_snp_table["Effect"] += effect
    finemapped_snp_table["pip"] += pip
    finemapped_snp_table["beta"] += beta


# current genomic coordinates are 1-based. Convert them to 0-based representation
finemapped_snp_table = pd.DataFrame(finemapped_snp_table)
finemapped_snp_table["Start"] = finemapped_snp_table["Start"] - 1
finemapped_snp_table["End"] = finemapped_snp_table["End"] - 1

# save rsids as txt file 
rsid_txt_filepath = osp.join(processed_data_dir, "finemap.all.rsids.txt")
with open(rsid_txt_filepath, "w") as f:
    f.write("\n".join(finemapped_snp_table["rsid"].tolist()))


###### Add hg38 locs to SNPs by querying dbSNP155

dbsnp155_outfilepath = osp.join(processed_data_dir, "snps.hg38_locs.bed")

# form bigBedNamedItems command to query hg38 locations txt files of snps 
cmd = [bigBedNamedItems, "-nameFile", ucsc_dbsnp1155, rsid_txt_filepath, dbsnp155_outfilepath]
result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
print(result.stdout)

# load hg38 locations
hg38_locs = pd.read_csv(dbsnp155_outfilepath, sep="\t", header=None)

# remove problematic chromosomes
hg38_locs = hg38_locs[hg38_locs.iloc[:,0].isin([f"chr{i}" for i in range(1, 22+1)])]

# following information is needed down below

'''
reference: https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/lib/bigDbSnp.as

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

dbSNP155bb_cols = ["chrom", "chromStart", "chromEnd", "name", "ref", "altCount", \
    "alts", "shiftBases", "freqSourceCount", "minorAlleleFreq", "majorAllele", "minorAllele", \
    "maxFuncImpact", "class", "ucscNotes", "_dataOffset", "_dataLen"]

'''
reference: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=dbSnp155Composite

For columns that contain lists of allele frequency data, the 
order of projects providing the data listed is as follows:
'''

snp_freq_source_order = ["1000Genomes", "dbGaP_PopFreq", "TOPMED", "KOREAN", "SGDP_PRJ", \
    "Qatari", "NorthernSweden", "Siberian", "TWINSUK", "TOMMO", "ALSPAC", "GENOME_DK", \
    "GnomAD", "GoNL", "Estonian", "Vietnamese", "Korea1K", "HapMap", "PRJEB36033", \
    "HGDP_Stanford", "Daghestan", "PAGE_STUDY", "Chileans", "MGP", "PRJEB37584", "GoESP", \
    "ExAC", "GnomAD_exomes", "FINRISK", "PharmGKB", "PRJEB37766"]


# add ucsc dbsnp155 table columns to the current table
hg38_locs.columns = dbSNP155bb_cols

# join by rsid
finemapped_snp_table["name"] = finemapped_snp_table["rsid"]
finemapped_snp_table = finemapped_snp_table.merge(hg38_locs, how="left", on="name")

# retain a subset of columns and save
finemapped_snp_table = finemapped_snp_table[['Chromosome', 'Start', 'End', 'A1', 'A2', 'rsid', 'Effect', 'pip',
       'beta', 'name', 'chrom', 'chromStart', 'chromEnd']]

# add hg19 coords
finemapped_snp_table["Chromosome_hg19"] = finemapped_snp_table["Chromosome"]
finemapped_snp_table["Start_hg19"] = finemapped_snp_table["Start"]
finemapped_snp_table["End_hg19"] = finemapped_snp_table["End"]

# add hg38 coords
finemapped_snp_table["Chromosome_hg38"] = finemapped_snp_table["chrom"]
finemapped_snp_table["Start_hg38"] = finemapped_snp_table["chromStart"]
finemapped_snp_table["End_hg38"] = finemapped_snp_table["chromEnd"]

# make Chromosome, Start and End based on hg38
finemapped_snp_table.drop(columns=["Chromosome", "Start", "End", "chrom", "chromStart", "chromEnd"], inplace=True)
finemapped_snp_table["Chromosome"] = finemapped_snp_table["Chromosome_hg38"]
finemapped_snp_table["Start"] = finemapped_snp_table["Start_hg38"]
finemapped_snp_table["End"] = finemapped_snp_table["End_hg38"]

# rename name:Name
finemapped_snp_table.rename(columns={"name":"Name"}, inplace=True)

# save table
finemapped_snp_table.to_csv(osp.join(processed_data_dir, \
    "snps.hg38_information.tsv"), index=False, sep="\t")


print("Script finished")
