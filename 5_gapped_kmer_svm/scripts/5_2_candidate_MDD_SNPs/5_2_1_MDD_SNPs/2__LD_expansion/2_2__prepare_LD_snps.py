import os
import pdb
import pickle

import os.path as osp
import pandas as pd
import numpy as np

from Bio.SeqUtils import GC


"""
This script loads LD expanded SNPs and maps them to MDD GWAS lead SNPs
"""

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


'''
Load LD expanded SNPs
'''

chroms = list(range(1,22+1))
ld_expand_dir = osp.join(tmp_dir, "gws.ld_expansion")
tag_list_template = "chr{}.high_corr.tags.list"
retain_cols = ["SNP", "TAGS"]

lead2ld_tables = []
ld_snp_tables = []
for chrom in chroms:

    filename = filename_template.format(chrom)

    # read LD expansion table for the current chromosome
    filename = tag_list_template.format(chrom)
    filepath = osp.join(ld_expand_dir, filename)

    if not osp.exists(ld_expand_path):
        print(f"Warning: ld expansion for chr{chrom} not found")
        continue

    ld_snp_table = pd.read_table(filepath, sep="\s+")
    ld_snp_table = ld_snp_table[retain_cols]

    # splits tags and format the table such that 
    # each row corresponds to a (lead snp, ld snp) pair
    temp = {"lead_snps":[], "ld_expanded_snps":[]}
    for idx, row in ld_snp_table.iterrows():
        lead = row["SNP"]
        lds = row["TAGS"].strip().split("|")
        for ld in lds:
            temp["lead_snps"].append(lead)
            temp["ld_expanded_snps"].append(ld)
    ld_snp_table = pd.DataFrame(temp)

    lead2ld_tables.append(ld_snp_table)

    # read bim file of the chromosome
    bim_table = pd.read_csv(osp.join(bimfile_dir, bimfile_template.format(i)), sep="\t", header=None)
    bim_table = bim_table.iloc[:, [0, 1, -2, -1]]
    bim_table.columns = ["Chromosome", "Name", "A1", "A2"]

    # add chr prefix to Chromosome column
    bim_table["Chromosome"] = bim_table["Chromosome"].apply(lambda x: f"chr{x}")

    # retain only specified rsids
    bim_table = bim_table[bim_table["Name"].isin(
        ld_snp_table["ld_expanded_snps"].unique().tolist())].drop(columns=["Chromosome"])

    ld_snp_tables.append(bim_table)

# merge per chromosome tables
ld_snp_tables = pd.concat(ld_snp_tables, axis=0)
ld_snp_tables.reset_index(drop=True, inplace=True)
lead2ld_tables = pd.concat(lead2ld_tables, axis=0)

# convert lead to ld mapping to dictionary
lead2ld = {}
for idx, row in lead2ld_tables.iterrows():
    lead = row["lead_snps"]
    ld = row["ld_expanded_snps"]
    if lead not in lead2ld:
        lead2ld[lead] = []
    lead2ld[lead].append(ld)

'''
Add hg38 locs to ld expanded snps
'''
# save snp rsids as txt file 
rsids_filepath = osp.join(tmp_dir, "ld_expand.rsids.txt")
with open(filepath, "w") as f:
    for rsid in snps["Name"].unique().tolist():
        f.write(f"{rsid}\n")

# form bigBedNamedItems command to query hg38 locations txt files of snps 
out_filepath = osp.join(tmp_dir, "ld_expand.dbsnp.bed")
cmd = [bigBedNamedItems, "-nameFile", dbsnp155_filepath, rsids_filepath, out_filepath]
result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)
print(result.stdout)

# read and prepare dbsnp out
hg38_locations = pd.read_csv(out_filepath,\
    sep="\t", header=None, names=dbsnp155_cols)
hg38_locations = hg38_locations[hg38_locations.iloc[:,0].isin([f"chr{i}" for i in chroms])]
hg38_locations = hg38_locations[["name", "chrom", "chromStart", "chromEnd"]]
hg38_locations = hg38_locations.rename(columns={
    "name": "Name",
    "chrom": "Chromosome",
    "chromStart": "Start",
    "chromEnd": "End"
})

# merge with 1000G alleles
ld_snp_tables.set_index("Name", inplace=True)
hg38_locations.set_index("Name", inplace=True)
hg38_locations = ld_snp_tables.join(hg38_locations, how="left").reset_index()

'''
For each MDD GWAS, match lead snps and ld snps
'''

# howard et al
filepath = osp.join(prepared_data_dir, "howard_et_al.gws.tsv")
howard_snps = pd.read_csv(filepath, sep="\t")

cur_lead_snps = howard_snps["Name"].unique().tolist()
cur_ld_snps = [ld for lead in cur_lead_snps for ld in lead2ld[lead]]
cur_ld_snps = list(set(cur_ld_snps) - set(cur_lead_snps))

howard_ld_snps = hg38_locations[hg38_locations["Name"].isin(cur_ld_snps)]

# levey et al
filepath = osp.join(prepared_data_dir, "levey_et_al.gws.tsv")
levey_snps = pd.read_csv(filepath, sep="\t")

cur_lead_snps = levey_snps["Name"].unique().tolist()
cur_ld_snps = [ld for lead in cur_lead_snps for ld in lead2ld[lead]]
cur_ld_snps = list(set(cur_ld_snps) - set(cur_lead_snps))

levey_ld_snps = hg38_locations[hg38_locations["Name"].isin(cur_ld_snps)]

# als et al
filepath = osp.join(prepared_data_dir, "als_et_al.gws.tsv")
als_snps = pd.read_csv(filepath, sep="\t")

cur_lead_snps = als_snps["Name"].unique().tolist()
cur_ld_snps = [ld for lead in cur_lead_snps for ld in lead2ld[lead]]
cur_ld_snps = list(set(cur_ld_snps) - set(cur_lead_snps))

als_ld_snps = hg38_locations[hg38_locations["Name"].isin(cur_ld_snps)]

ld_expanded_snps = pd.concat([howard_ld_snps, levey_ld_snps, als_ld_snps], axis=0)


'''
LD friends (r2 > 0.8) for Als et al SNPs
'''


filepath = osp.join(raw_data_dir, "MDD_GWAS/als_et_al/als_et_al.supptable2.xlsx")
table = pd.read_excel(filepath)

# the loaded table provides LD friends (r2 > 0.6) for each index snp
ld_friends_06 = table["LD-friends(0.6).p0.001.index.SNP"].tolist()
ld_friends_06 = [x.split(",") for x in ld_friends_06]

# extract r2 and rsids of LD friends
ld_friends_06_r2 = []
ld_friends_06_rsid = []
ld_friends_06_signed_kb_dist = []
for ld_friend_list in ld_friends_06:
    ld_friends_06_r2.append([])
    ld_friends_06_rsid.append([])
    ld_friends_06_signed_kb_dist.append([])
    for ld_friend in ld_friend_list:
        # extract r2, rsid and signed distance in kb for the friend
        try:
            temp = ld_friend.split("(")
            ld_friends_06_rsid[-1].append(temp[0])
            temp = temp[1].split(")")[0]
            temp = temp.split("/")
            ld_friends_06_r2[-1].append(float(temp[0]))
            ld_friends_06_signed_kb_dist[-1].append(float(temp[1]))
        except:
            print("Invalid LD friend:", ld_friend)
            continue  


# subset LD friend snps having r2 > .8 
ld_friends_08_r2 = []
ld_friends_08_rsid = []
ld_friends_08_signed_kb_dist = []
for outer_idx in range(len(ld_friends_06)):
    ld_friends_08_r2.append([])
    ld_friends_08_rsid.append([])
    ld_friends_08_signed_kb_dist.append([])
    for r2, rsid, dist in zip(ld_friends_06_r2[outer_idx], ld_friends_06_rsid[outer_idx], ld_friends_06_signed_kb_dist[outer_idx]):
        if r2 > .8:
            ld_friends_08_r2[-1].append(r2)
            ld_friends_08_rsid[-1].append(rsid)
            ld_friends_08_signed_kb_dist[-1].append(dist)

# create a 2-level dictionary from the list of lists
# first level: dictionary with als_index_snp_rsid s as keys
# second level: dictionary with "r2", "rsid", "signed_kb_dist" as keys
ld_snps_dct = {}
ld_friend_rsids = []
for index_snp, friend_r2, friend_rsid, friend_dist in zip(als_index_snp_rsids, ld_friends_08_r2, ld_friends_08_rsid, ld_friends_08_signed_kb_dist):
    ld_snps_dct[index_snp] = {}
    ld_snps_dct[index_snp]["r2"] = friend_r2
    ld_snps_dct[index_snp]["rsid"] = friend_rsid
    ld_snps_dct[index_snp]["signed_kb_dist"] = friend_dist
    ld_friend_rsids += friend_rsid

# remove duplicate ld friends
ld_friend_rsids = list(set(ld_friend_rsids))

# remove index snps from ld friends
ld_friend_rsids = [x for x in ld_friend_rsids if x not in als_index_snp_rsids]

# save ld friends as a text file
rsids_filepath = f"{GKMSVM_PREPARED_DATA_DIR}/ld_friends.txt"
with open(rsids_filepath, "w") as f:
    for ld_friend in ld_friend_rsids:
        f.write(f"{ld_friend}\n")

os.makedirs(f"{GKMSVM_PREPARED_DATA_DIR}/ld_friend_chunks", exist_ok=True)

# split them to 1000 snp length chunks and save each to a txt file 
counter = 1
prev_i = 0
for i in range(1000, len(ld_friend_rsids), 1000):
    if prev_i + 1000 > len(ld_friend_rsids):
        i = len(ld_friend_rsids)
    temp = ld_friend_rsids[prev_i:i]
    prev_i = i
    with open(f"{GKMSVM_PREPARED_DATA_DIR}/ld_friend_chunks/ld_friends-{counter}.txt", "w") as f:
        for rsid in temp:
            f.write(f"{rsid}\n")  
    counter += 1

os.makedirs(f"{GKMSVM_PREPARED_DATA_DIR}/ld_friend_locs", exist_ok=True)

for i in tqdm(range(1, counter)):
    rsid_txt_filepath = f"{GKMSVM_PREPARED_DATA_DIR}/ld_friend_chunks/ld_friends-{i}.txt"
    dbsnp155_outfilepath = f"{GKMSVM_PREPARED_DATA_DIR}/ld_friend_locs/ld_friends-{i}.bed"

    # form bigBedNamedItems command to query hg38 locations txt files of snps 
    cmd = [bigBedNamedItems_filepath, "-nameFile", ucsc_dbsnp155_filepath, rsid_txt_filepath, dbsnp155_outfilepath]
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True)

# load hg38 locations of chunks and concatenate them together
hg38_locations = []
for i in tqdm(range(1, counter)):
    # load hg38 location of the snps in current chunk
    temp = pd.read_csv(f"{GKMSVM_PREPARED_DATA_DIR}/ld_friend_locs/ld_friends-{i}.bed",\
        sep="\t", header=None, names=dbSNP155bb_cols)
    # remove problematic chromosomes
    temp = temp[temp.iloc[:,0].isin([f"chr{i}" for i in range(1, 22+1)])]
    # append
    hg38_locations.append(temp)

hg38_locations = pd.concat(hg38_locations, axis=0)

# read and prepare dbsnp out
hg38_locations = pd.read_csv(out_filepath,\
    sep="\t", header=None, names=dbsnp155_cols)
hg38_locations = hg38_locations[hg38_locations.iloc[:,0].isin([f"chr{i}" for i in chroms])]
hg38_locations = hg38_locations[["name", "chrom", "chromStart", "chromEnd"]]
hg38_locations = hg38_locations.rename(columns={
    "name": "Name",
    "chrom": "Chromosome",
    "chromStart": "Start",
    "chromEnd": "End"
})

# get A1 and A2 alleles for each allele from 1kg .bim file
chrom_grouped = hg38_locations.groupby("chrom")
bim_tables = []
for chrom, chrom_df in chrom_grouped:
    i = chrom[3:]

    # read bim file of the chromosome
    bim_table = pd.read_csv(osp.join(ld_files_basepath, f"1000G.EUR.QC.{i}.bim"), sep="\t", header=None)
    bim_table = bim_table.iloc[:, [0, 1, -2, -1]]
    bim_table.columns = ["Chromosome", "Name", "A1", "A2"]

    # add chr prefix to Chromosome column
    bim_table["Chromosome"] = bim_table["Chromosome"].apply(lambda x: f"chr{x}")

    # retain only specified rsids
    bim_table = bim_table[bim_table["Name"].isin(chrom_df["name"].tolist())]
    bim_tables.append(bim_table)

# merge per chromosome tables
bim_table = pd.concat(bim_tables, axis=0)
bim_table.reset_index(drop=True, inplace=True)

# join location and allele information tables
hg38_locations.rename(columns={"name":"Name"}, inplace=True)

# join by rsid
ld_friend_table = hg38_locations.merge(bim_table, how="left", on="Name")

# merge ld friends with ld expanded snps
ld_expanded_snps.set_index("Name", inplace=True)
ld_friend_table.set_index("Name", inplace=True)
ld_expanded_snps = ld_expanded_snps.join(ld_friend_table, how="left").reset_index()

# save ld expanded snps
filepath = osp.join(GKMSVM_PREPARED_DATA_DIR, "snps/ld_expand.tsv")
snps.to_csv(filepath, sep="\t")

