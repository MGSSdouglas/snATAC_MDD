import pdb 
import os
import subprocess

import os.path as osp
import pandas as pd 
import numpy as np 


# ------- Configuration -------

# get environment variable values 
GKMSVM_WORKSPACE_DIR = os.environ.get("GKMSVM_WORKSPACE_DIR")
GKMSVM_ENV_DIR = os.environ.get("GKMSVM_ENV_DIR")
GKMSVM_SCRIPTS_DIR = os.environ.get("GKMSVM_SCRIPTS_DIR")
GKMSVM_RAW_DATA_DIR = os.environ.get("GKMSVM_RAW_DATA_DIR")
GKMSVM_PREPARED_DATA_DIR = os.environ.get("GKMSVM_PREPARED_DATA_DIR")
GKMSVM_MODEL_DIR = os.environ.get("GKMSVM_MODEL_DIR")
GKMSVM_TMP_DIR = os.environ.get("GKMSVM_TMP_DIR")
GKMSVM_BIN_DIR = os.environ.get("GKMSVM_BIN_DIR")

# create temporary path for current analysis
tmp_dir = osp.join(GKMSVM_TMP_DIR, "5_2_candidate_MDD_SNPs")
os.makedirs(tmp_dir, exist_ok=True)

# source filepaths
data_id = "snps/lead_MDD_GWAS_SNPs"
raw_data_dir = osp.join(GKMSVM_RAW_DATA_DIR, data_id)

# ucsc utils to query dbSNP locally 
data_id = "dbsnp155/dbSnp155.bb"
dbsnp155_filepath = osp.join(GKMSVM_RAW_DATA_DIR, data_id) 

script_id = "ucsc_utils/bigBedNamedItems"
bigBedNamedItems = osp.join(GKMSVM_BIN_DIR, script_id)

# dbSNP columns reference: https://github.com/ucscGenomeBrowser/kent/blob/master/src/hg/lib/bigDbSnp.as
dbsnp155_cols = ["chrom", "chromStart", "chromEnd", "name", "ref", "altCount", \
    "alts", "shiftBases", "freqSourceCount", "minorAlleleFreq", "majorAllele", "minorAllele", \
    "maxFuncImpact", "class", "ucscNotes", "_dataOffset", "_dataLen"]

# 1kG EUR plinkfiles
data_id = "1000G_EUR_Phase3_plinkfiles"
bimfile_dir = osp.join(GKMSVM_RAW_DATA_DIR, data_id)
bimfile_template = "1000G.EUR.QC.{}"

# target filepaths
data_id = "snps/lead_MDD_GWAS_SNPs"
prepared_data_dir = osp.join(GKMSVM_PREPARED_DATA_DIR, data_id)


# ------- lead SNPs from Als et al. MDD GWAS -------

# load snps
filepath = osp.join(raw_data_dir, "als_et_al.supptable2.xlsx")
snps = pd.read_excel(filepath)

# retain and rename relevant columns 
snps = snps[["index.SNP", "CHR", "P.x.index.SNP", "BP.x.index.SNP", "OR.x.index.SNP", "SE.index.SNP", "A1A2.index.SNP"]]
snps["A1"] =  [x.split("/")[0] for x in snps["A1A2.index.SNP"]]
snps["A2"] =  [x.split("/")[1] for x in snps["A1A2.index.SNP"]]
snps["CHR"] = snps.CHR.map(lambda x: f"chr{x}")
snps = snps.rename(columns={
    "index.SNP": "Name",
    "CHR": "Chromosome_hg19",
    "P.x.index.SNP": "pval",
    "BP.x.index.SNP": "End_hg19",
    "OR.x.index.SNP": "OR",
    "SE.index.SNP": "SE"
}).drop(columns=["A1A2.index.SNP"])
snps["Start_hg19"] = snps["End_hg19"] - 1

# save snp rsids as txt file 
rsids_filepath = osp.join(tmp_dir, "als_et_al.gws.rsids.txt")
with open(filepath, "w") as f:
    for rsid in snps["Name"].unique().tolist():
        f.write(f"{rsid}\n")

# form bigBedNamedItems command to query hg38 locations txt files of snps 
out_filepath = osp.join(tmp_dir, "als_et_al.gws.dbsnp.bed")
cmd = [bigBedNamedItems, "-nameFile", dbsnp155_filepath, rsids_filepath, out_filepath]
result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
print(result.stdout)

# read and prepare dbsnp out
hg38_locations = pd.read_csv(out_filepath,\
    sep="\t", header=None, names=dbsnp155_cols)
hg38_locations = hg38_locations[hg38_locations.iloc[:,0].isin([f"chr{i}" for i in range(1, 22+1)])]
hg38_locations = hg38_locations[["name", "chrom", "chromStart", "chromEnd"]]
hg38_locations = hg38_locations.rename(columns={
    "name": "Name",
    "chrom": "Chromosome",
    "chromStart": "Start",
    "chromEnd": "End"
})

# merge locations and snps
snps.set_index("Name", inplace=True)
hg38_locations.set_index("Name", inplace=True)
snps = snps.join(hg38_locations, how="left").reset_index()

# save snps
filepath = osp.join(prepared_data_dir, "als_et_al.gws.tsv")
snps.to_csv(filepath, sep="\t")


# ------- lead SNPs from Levey et al. MDD GWAS -------

# load snps
filepath = osp.join(raw_data_dir, "levey_et_al.supptable.tsv")
snps = pd.read_csv(filepath, sep="\t")
snps.rename(columns={"rsid": "Name"}, inplace=True)

# save snp rsids as txt file 
rsids_filepath = osp.join(tmp_dir, "levey_et_al.gws.rsids.txt")
with open(filepath, "w") as f:
    for rsid in snps["Name"].unique().tolist():
        f.write(f"{rsid}\n")

# form bigBedNamedItems command to query hg38 locations txt files of snps 
out_filepath = osp.join(tmp_dir, "levey_et_al.gws.dbsnp.bed")
cmd = [bigBedNamedItems, "-nameFile", dbsnp155_filepath, rsids_filepath, out_filepath]
result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
print(result.stdout)

# read and prepare dbsnp out
hg38_locations = pd.read_csv(out_filepath,\
    sep="\t", header=None, names=dbsnp155_cols)
hg38_locations = hg38_locations[hg38_locations.iloc[:,0].isin([f"chr{i}" for i in range(1, 22+1)])]
hg38_locations = hg38_locations[["name", "chrom", "chromStart", "chromEnd"]]
hg38_locations = hg38_locations.rename(columns={
    "name": "Name",
    "chrom": "Chromosome",
    "chromStart": "Start",
    "chromEnd": "End"
})

# merge locations and snps
snps.set_index("Name", inplace=True)
hg38_locations.set_index("Name", inplace=True)
snps = snps.join(hg38_locations, how="left").reset_index()

# group snps by chromosomes and query 1kG bimfiles
chrom_groups = snps.groupby("Chromosome")
bim_tables = []
for chrom, chrom_snps in chrom_groups:

    # read and prepare 1kG bimfile for cur chrom
    chrom = chrom[3:]
    filepath = bimfile_template.format(chrom)
    bim_table = bim_table.iloc[:, [0, 1, -2, -1]]
    bim_table.columns = ["Name", "A1", "A2"]

    # retain current snps in the bim and save
    bim_table = bim_table[bim_table["Name"].isin(chrom_snps["Name"].tolist())]
    bim_tables.append(bim_table)

# merge filtered bim file tables across chroms
bim_table = pd.concat(bim_tables, axis=0).reset_index(drop=True)

# merge with snp table 
snps.set_index("Name", inplace=True)
bim_table.set_index("Name", inplace=True)
snps = snps.join(bim_table, how="left").reset_index()

# save snps
filepath = osp.join(prepared_data_dir, "levey_et_al.gws.tsv")
snps.to_csv(filepath, sep="\t")

# save snps per chromosome 
for chrom, chrom_snps in snps.groupby("Chromosome"):
    filepath = osp.join(tmp_dir, f"levey_et_al.gws.{chrom[3:]}.rsids.txt")
    with open(filepath, "w") as f:
        for rsid in chrom_snps["Name"].unique().tolist():
            f.write(f"{rsid}\n")


# ------- lead SNPs from Howard et al. MDD GWAS -------

# load snps
filepath = osp.join(raw_data_dir, "howard_et_al.supptable1.tsv")
snps = pd.read_csv(filepath, sep="\t")
snps["Chromosome"] = "chr" + snps["Chromosome"].astype(str)
snps = snps[["Name", "Chromosome", "Start", "End"]].rename(columns={
    "Chromosome": "Chromosome_hg19",
    "Start": "Start_hg19",
    "End": "End_hg19"
})

# save snp rsids as txt file 
rsids_filepath = osp.join(tmp_dir, "howard_et_al.gws.rsids.txt")
with open(filepath, "w") as f:
    for rsid in snps["Name"].unique().tolist():
        f.write(f"{rsid}\n")

# form bigBedNamedItems command to query hg38 locations txt files of snps 
out_filepath = osp.join(tmp_dir, "howard_et_al.gws.dbsnp.bed")
cmd = [bigBedNamedItems, "-nameFile", dbsnp155_filepath, rsids_filepath, out_filepath]
result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
print(result.stdout)

# read and prepare dbsnp out
hg38_locations = pd.read_csv(out_filepath,\
    sep="\t", header=None, names=dbsnp155_cols)
hg38_locations = hg38_locations[hg38_locations.iloc[:,0].isin([f"chr{i}" for i in range(1, 22+1)])]
hg38_locations = hg38_locations[["name", "chrom", "chromStart", "chromEnd"]]
hg38_locations = hg38_locations.rename(columns={
    "name": "Name",
    "chrom": "Chromosome",
    "chromStart": "Start",
    "chromEnd": "End"
})

# merge locations and snps
snps.set_index("Name", inplace=True)
hg38_locations.set_index("Name", inplace=True)
snps = snps.join(hg38_locations, how="left").reset_index()

# group snps by chromosomes and query 1kG bimfiles
chrom_groups = snps.groupby("Chromosome")
bim_tables = []
for chrom, chrom_snps in chrom_groups:

    # read and prepare 1kG bimfile for cur chrom
    chrom = chrom[3:]
    filepath = bimfile_template.format(chrom)
    bim_table = bim_table.iloc[:, [0, 1, -2, -1]]
    bim_table.columns = ["Name", "A1", "A2"]

    # retain current snps in the bim and save
    bim_table = bim_table[bim_table["Name"].isin(chrom_snps["Name"].tolist())]
    bim_tables.append(bim_table)

# merge filtered bim file tables across chroms
bim_table = pd.concat(bim_tables, axis=0).reset_index(drop=True)

# merge with snp table 
snps.set_index("Name", inplace=True)
bim_table.set_index("Name", inplace=True)
snps = snps.join(bim_table, how="left").reset_index()

# save snps
filepath = osp.join(prepared_data_dir, "howard_et_al.gws.tsv")
snps.to_csv(filepath, sep="\t")

# save snps per chromosome 
for chrom, chrom_snps in snps.groupby("Chromosome"):
    filepath = osp.join(tmp_dir, f"howard_et_al.gws.{chrom[3:]}.rsids.txt")
    with open(filepath, "w") as f:
        for rsid in chrom_snps["Name"].unique().tolist():
            f.write(f"{rsid}\n")


print("Script finished")
