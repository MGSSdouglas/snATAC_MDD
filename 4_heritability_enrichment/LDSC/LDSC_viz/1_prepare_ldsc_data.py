import pdb 
import os 
import sys
import matplotlib
matplotlib.use('agg')

import numpy as np 
import pandas as pd 
import os.path as osp


###### Configuration

data_dir = "./data/raw/LDSC_for_Doruk"

all_marker_data_dir = osp.join(data_dir, "Cts_LDSC_all_marker")
DAR_data_dir = osp.join(data_dir, "Cts_LDSC_DAR")

processed_data_dir = "./data/processed"
os.makedirs(processed_data_dir, exist_ok=True)


###### Load data



### Results on all and marker peaks

## get filenames
all_marker_ldsc_filenames = os.listdir(all_marker_data_dir)

# retain the ones having txt extension
all_marker_ldsc_filenames = [f for f in all_marker_ldsc_filenames if f.endswith(".txt")]

# split the ones containing broad and subclusters
subcluster_all_marker_ldsc_filenames = [f for f in all_marker_ldsc_filenames if "sub" in f]


### subcluster marker peak results 
subcluster_marker_ldsc_filenames = [f for f in subcluster_all_marker_ldsc_filenames if "marker" in f]


# extract results 
subcluster_marker_ldsc_results = {}
subcluster_marker_available_cell_types = []
for filename in subcluster_marker_ldsc_filenames:

    # construct the filepath
    filepath = osp.join(all_marker_data_dir, filename)

    # read the file
    table = pd.read_table(filepath, sep="\s+")

    # extract trait name 
    trait_name = filename.split("_")[0]

    # add trait name and peak type as column to the table
    table["Trait"] = trait_name
    table["PeakType"] = "marker"

    # rename some values in Name column 
    table["Name"] = table["Name"].apply(lambda x: cell_type_name_refactor_dct[x] if x in cell_type_name_refactor_dct else x)

    # add to dictionary
    subcluster_marker_ldsc_results[trait_name] = table

    # extract available cell types
    subcluster_marker_available_cell_types += table["Name"].unique().tolist()

# retain unique cell types
subcluster_marker_available_cell_types = list(set(subcluster_marker_available_cell_types))



### Sanity check: check if all cell types are available in all results 
for trait in subcluster_marker_ldsc_results:

    # get Name column of the table
    cell_types = subcluster_marker_ldsc_results[trait]["Name"].tolist()

    # check if all cell types are available
    try:
        assert all([c in cell_types for c in subcluster_marker_available_cell_types])
    except:
        print("[MARKER] Trait {} does not contain the following cell types {}".format(trait, [c for c in subcluster_marker_available_cell_types if c not in cell_types]))
        




### subcluster all peak results
subcluster_all_ldsc_filenames = [f for f in subcluster_all_marker_ldsc_filenames if "all" in f]

subcluster_all_ldsc_results = {}
subcluster_all_available_cell_types = []
for filename in subcluster_all_ldsc_filenames:
    
    # construct the filepath
    filepath = osp.join(all_marker_data_dir, filename)

    # read the file
    table = pd.read_table(filepath, sep="\s+")

    # extract trait name 
    trait_name = filename.split("_")[0]

    # add trait name and peak type as column to the table
    table["Trait"] = trait_name
    table["PeakType"] = "all"

    # rename some values in Name column 
    table["Name"] = table["Name"].apply(lambda x: cell_type_name_refactor_dct[x] if x in cell_type_name_refactor_dct else x)

    # add to dictionary
    subcluster_all_ldsc_results[trait_name] = table

    # extract available cell types
    subcluster_all_available_cell_types += table["Name"].unique().tolist()

# retain unique cell types
subcluster_all_available_cell_types = list(set(subcluster_all_available_cell_types))



### Sanity check: check if all cell types are available in all results
for trait in subcluster_all_ldsc_results:

    cell_types = subcluster_all_ldsc_results[trait]["Name"].tolist()

    # check if all cell types are available
    try:
        assert all([c in cell_types for c in subcluster_all_available_cell_types])
    except:
        print("[ALL] Trait {} does not contain the following cell types {}".format(trait, [c for c in subcluster_all_available_cell_types if c not in cell_types]))
        



### Results on DAR peaks

# get filenames
DAR_ldsc_filenames = os.listdir(DAR_data_dir)
# retain the ones having txt extension
DAR_ldsc_filenames = [f for f in DAR_ldsc_filenames if f.endswith(".txt")]

# DAR peak results
DAR_ldsc_results = {}
subcluster_DAR_available_cell_types = []
for filename in DAR_ldsc_filenames:

    # construct the filepath
    filepath = osp.join(DAR_data_dir, filename)

    # read the file
    table = pd.read_table(filepath, sep="\s+")

    # extract trait name 
    trait_name = filename.split("_")[0]

    # add trait name and peak type as column to the table
    table["Trait"] = trait_name
    table["PeakType"] = "DAR"

    # rename some values in Name column 
    table["Name"] = table["Name"].apply(lambda x: cell_type_name_refactor_dct[x] if x in cell_type_name_refactor_dct else x)

    # add to dictionary
    DAR_ldsc_results[trait_name] = table

    # extract available cell types
    subcluster_DAR_available_cell_types += table["Name"].unique().tolist()

# retain unique cell types
subcluster_DAR_available_cell_types = list(set(subcluster_DAR_available_cell_types))


### Sanity check: check if all cell types are available in all results
for trait in DAR_ldsc_results:

    cell_types = DAR_ldsc_results[trait]["Name"].tolist()

    # check if all cell types are available
    try:
        assert all([c in cell_types for c in subcluster_DAR_available_cell_types])
    except:
        print("[DAR] Trait {} does not contain the following cell types {}".format(trait, [c for c in subcluster_DAR_available_cell_types if c not in cell_types]))
        



###### Convert dictionaries to dataframes

# convert each dct to dataframe 
subcluster_marker_ldsc_results_df = pd.concat(subcluster_marker_ldsc_results.values(), ignore_index=True)
subcluster_all_ldsc_results_df = pd.concat(subcluster_all_ldsc_results.values(), ignore_index=True)
DAR_ldsc_results_df = pd.concat(DAR_ldsc_results.values(), ignore_index=True)

# concatenate all dataframes
all_ldsc_results_df = pd.concat([subcluster_marker_ldsc_results_df, subcluster_all_ldsc_results_df, DAR_ldsc_results_df], ignore_index=True)

# reset index 
all_ldsc_results_df.reset_index(inplace=True, drop=True)


###### Save the table
filepath = "all_ldsc_results.tsv"
filepath = osp.join(processed_data_dir, filepath)
all_ldsc_results_df.to_csv(filepath, sep="\t", index=False)


print("Script finished")
