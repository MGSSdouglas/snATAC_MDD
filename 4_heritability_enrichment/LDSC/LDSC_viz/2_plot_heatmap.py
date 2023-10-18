import pdb 
import os 
import sys
import matplotlib
matplotlib.use('agg')

import numpy as np 
import pandas as pd 
import os.path as osp
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.colors import colorConverter
from matplotlib.colors import to_hex
from matplotlib.colors import LogNorm
from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
from pprint import pprint
from config import *


###### Configuration

processed_data_dir = "./data/processed"

results_dir = "./results"
os.makedirs(results_dir, exist_ok=True)

# subcluster ordering for the heatmap
subcluster_ordering = ['Ast1', 'Ast2', 'Ast3', 'Ast4', 'End1', 'End2', 'ExN1', 'ExN2', 'ExN3', 'ExN4', 'ExN5', 'ExN6', 'ExN7', 'ExN8', 'ExN9', 'ExN10', 'ExN11', 'ExN12', 'InN3', 'InVIP', 'InLAMP5', 'InPV', 'InSST', 'Mic1', 'Mic2', 'Oli1', 'Oli2', 'Oli3', 'Oli4', 'Oli5', 'Oli6', 'Oli7', 'OPC1', 'OPC2', 'OPC3', 'OPC4']

# trait short name to full name mapping 
trait_short2full = {
    "ADHD": "ADHD",
    "SA": "Suicide Attempt",
    "BMD": "Bone Mineral Density",
    "BMI": "BMI",
    "NEU": "Neuroticism",
    "MDD": "MDD (Howard et al. - 2019)",
    "MDDALS": "MDD (Als et al. - 2023)",
    "HGT": "Height",
    "BIP": "Bipolar Disorder",
    "INS": "Insomnia",
    "SCZ": "Schizophrenia",
}

trait_short_name_ordering = ["MDDALS", "MDD", "INS", "NEU", "BIP", "SCZ", "ADHD", "SA", "BMD", "BMI", "HGT"]

###### Load data

# load all ldsc results 
filepath = osp.join(processed_data_dir, "all_ldsc_results.tsv")
combined_ldsc_results = pd.read_csv(filepath, sep="\t")

# retain rows whose Name is in subcluster_ordering
combined_ldsc_results = combined_ldsc_results[combined_ldsc_results["Name"].isin(subcluster_ordering)]

# rename Name col to old_names
combined_ldsc_results = combined_ldsc_results.rename(columns={"Name": "old_names"})

# map old names to new names
combined_ldsc_results["new_names"] = combined_ldsc_results["old_names"].map(old2new)

# split marker DAR and all
all_ldsc_results = combined_ldsc_results[combined_ldsc_results["PeakType"] == "all"]
marker_ldsc_results = combined_ldsc_results[combined_ldsc_results["PeakType"] == "marker"]
DAR_ldsc_results = combined_ldsc_results[combined_ldsc_results["PeakType"] == "DAR"]


###### Create heatmaps 

### all

# find number of cell types in the table 
num_cell_types = all_ldsc_results["old_names"].nunique()

# create a matrix of all the data
all_heatmap = np.zeros((num_cell_types, len(trait_short_name_ordering)), dtype=float)

# iterate over traits
cell_types_to_remove = []
for i, trait in enumerate(trait_short_name_ordering):

    index_shift_factor = 0

    # iterate over cell types
    for j, cell_type in enumerate(subcluster_ordering_new_names):

        row = all_ldsc_results[(all_ldsc_results["new_names"] == cell_type) & (all_ldsc_results["Trait"] == trait)]

        if len(row) == 0:
            print("[ALL] No row found for cell type {} and trait {}".format(cell_type, trait))
            index_shift_factor += 1
            cell_types_to_remove.append(cell_type)
            continue
        elif len(row) == 1:
            all_heatmap[j - index_shift_factor, i] = row["Coefficient_P_value"].values[0]
        else:
            print("[ALL] More than one row found for cell type {} and trait {}".format(cell_type, trait))
            pdb.set_trace()

# remove unavailable cell types from cell type names
new_cell_types = [x for x in subcluster_ordering_new_names if x not in cell_types_to_remove]

# map trait short names to full names
trait_full_name_ordering = [trait_short2full[x] for x in trait_short_name_ordering]

# create a dataframe from the matrix
all_heatmap_df = pd.DataFrame(all_heatmap, index=new_cell_types, columns=trait_full_name_ordering)

# normalize value in dataframe as -np.log10(p-value)
all_heatmap_df = -np.log10(all_heatmap_df)

### marker

# find number of cell types in the table
num_cell_types = marker_ldsc_results["old_names"].nunique()

# create a matrix of all the data
marker_heatmap = np.zeros((num_cell_types, len(trait_short_name_ordering)), dtype=float)

# iterate over traits
cell_types_to_remove = []
for i, trait in enumerate(trait_short_name_ordering):

    index_shift_factor = 0

    # iterate over cell types
    for j, cell_type in enumerate(subcluster_ordering_new_names):

        row = marker_ldsc_results[(marker_ldsc_results["new_names"] == cell_type) & (marker_ldsc_results["Trait"] == trait)]

        if len(row) == 0:
            print("[MARKER] No row found for cell type {} and trait {}".format(cell_type, trait))
            index_shift_factor += 1
            cell_types_to_remove.append(cell_type)
            continue
        elif len(row) == 1:
            marker_heatmap[j - index_shift_factor, i] = row["Coefficient_P_value"].values[0]
        else:
            print("[MARKER] More than one row found for cell type {} and trait {}".format(cell_type, trait))
            pdb.set_trace()

# remove unavailable cell types from cell type names
new_cell_types = [x for x in subcluster_ordering_new_names if x not in cell_types_to_remove]

# map trait short names to full names
trait_full_name_ordering = [trait_short2full[x] for x in trait_short_name_ordering]

# create a dataframe from the matrix
marker_heatmap_df = pd.DataFrame(marker_heatmap, index=new_cell_types, columns=trait_full_name_ordering)

# normalize value in dataframe as -np.log10(p-value)
marker_heatmap_df = -np.log10(marker_heatmap_df)


### DAR

# find number of cell types in the table
num_cell_types = DAR_ldsc_results["old_names"].nunique()

# create a matrix of all the data
DAR_heatmap = np.zeros((num_cell_types, len(trait_short_name_ordering)), dtype=float)

# iterate over traits
cell_types_to_remove = []
for i, trait in enumerate(trait_short_name_ordering):

    index_shift_factor = 0

    # iterate over cell types
    for j, cell_type in enumerate(subcluster_ordering_new_names):

        row = DAR_ldsc_results[(DAR_ldsc_results["new_names"] == cell_type) & (DAR_ldsc_results["Trait"] == trait)]

        if len(row) == 0:
            print("[DAR] No row found for cell type {} and trait {}".format(cell_type, trait))
            index_shift_factor += 1
            cell_types_to_remove.append(cell_type)
            continue
        elif len(row) == 1:
            DAR_heatmap[j - index_shift_factor, i] = row["Coefficient_P_value"].values[0]
        else:
            print("[DAR] More than one row found for cell type {} and trait {}".format(cell_type, trait))
            pdb.set_trace()

# remove unavailable cell types from cell type names
new_cell_types = [x for x in subcluster_ordering_new_names if x not in cell_types_to_remove]

# map trait short names to full names
trait_full_name_ordering = [trait_short2full[x] for x in trait_short_name_ordering]

# create a dataframe from the matrix
DAR_heatmap_df = pd.DataFrame(DAR_heatmap, index=new_cell_types, columns=trait_full_name_ordering)

# normalize value in dataframe as -np.log10(p-value)
DAR_heatmap_df = -np.log10(DAR_heatmap_df)


###### Harmonize heatmap columns based on all heatmap

# iterate over indices of all heatmap 
for i in range(all_heatmap_df.shape[0]):

    # get index
    cell_type = all_heatmap_df.index[i]

    ### marker 

    # check if the cell type is present among indices of marker heatmap 
    if cell_type not in marker_heatmap_df.index.tolist():

        # append a rowof -1s to marker heatmap
        marker_heatmap_df.loc[cell_type] = -1
    
    ### DAR

    # check if the cell type is present among indices of DAR heatmap
    if cell_type not in DAR_heatmap_df.index.tolist():

        # append a rowof -1s to DAR heatmap
        DAR_heatmap_df.loc[cell_type] = -1

# order the rows of marker heatmap based on all heatmap
marker_heatmap_df = marker_heatmap_df.loc[all_heatmap_df.index]
    
# order the rows of DAR heatmap based on all heatmap
DAR_heatmap_df = DAR_heatmap_df.loc[all_heatmap_df.index]

assert DAR_heatmap_df.shape == marker_heatmap_df.shape == all_heatmap_df.shape


###### Transpose heatmap data 

# transpose all heatmap
all_heatmap_df = all_heatmap_df.T

# transpose marker heatmap
marker_heatmap_df = marker_heatmap_df.T

# transpose DAR heatmap
DAR_heatmap_df = DAR_heatmap_df.T


###### Create masks for missing data in heatmaps

# create a mask for marker
marker_heatmap_mask = marker_heatmap_df == -1

# create a mask for DAR
DAR_heatmap_mask = DAR_heatmap_df == -1

###### Create annotations for heatmaps


### all
all_annots = np.chararray(all_heatmap_df.shape, itemsize=10, unicode=True)
all_annots[:] = ""
for i in range(all_heatmap_df.shape[0]):
    for j in range(all_heatmap_df.shape[1]):
        if all_heatmap_df.iloc[i, j] >= 3:
            all_annots[i, j] = "**"
        elif all_heatmap_df.iloc[i, j] >= 1.3:
            all_annots[i, j] = "*"

### marker 
marker_annots = np.chararray(marker_heatmap_df.shape, itemsize=10, unicode=True)
marker_annots[:] = ""
for i in range(marker_heatmap_df.shape[0]):
    for j in range(marker_heatmap_df.shape[1]):
        if not marker_heatmap_mask.iloc[i, j]:
            if marker_heatmap_df.iloc[i, j] >= 3:
                marker_annots[i, j] = "**"
            elif marker_heatmap_df.iloc[i, j] >= 1.3:
                marker_annots[i, j] = "*"

### DAR
DAR_annots = np.chararray(DAR_heatmap_df.shape, itemsize=10, unicode=True)
DAR_annots[:] = ""
for i in range(DAR_heatmap_df.shape[0]):
    for j in range(DAR_heatmap_df.shape[1]):
        if not DAR_heatmap_mask.iloc[i, j]:
            if DAR_heatmap_df.iloc[i, j] >= 3:
                DAR_annots[i, j] = "**"
            elif DAR_heatmap_df.iloc[i, j] >= 1.3:
                DAR_annots[i, j] = "*"


pdb.set_trace()



###### Plot heatmaps (coral)

# cmap = sns.color_palette("coolwarm", as_cmap=True)
cmap = sns.color_palette("light:coral", as_cmap=True)

color = '#9ca09e'
missing_cmap = ListedColormap([color] + plt.get_cmap('Blues')(range(9)).tolist())


### all

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

min_neglog10_pval = all_heatmap_df[all_heatmap_df>=0].min().min()
max_neglog10_pval = all_heatmap_df[all_heatmap_df>=0].max().max()

# plot heatmap 
sns.heatmap(
    data=all_heatmap_df,
    annot=all_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save 
filepath = osp.join(results_dir, "ldsc.heatmap.coral.beforeIPR.pdf")
plt.savefig(filepath, bbox_inches='tight')


### marker

min_neglog10_pval = marker_heatmap_df[marker_heatmap_df>=0].min().min()
max_neglog10_pval = marker_heatmap_df[marker_heatmap_df>=0].max().max()

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

# plot non-missing values as heatmap
sns.heatmap(
    data=marker_heatmap_df,
    annot=marker_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    mask=marker_heatmap_mask,
    ax=ax,
    rasterized=False
)

# plot missing values as dark gray_r 
sns.heatmap(
    data=marker_heatmap_df,
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=missing_cmap,
    cbar=False,
    mask=~marker_heatmap_mask,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save
filepath = osp.join(results_dir, "ldsc.heatmap.coral.marker.pdf")
plt.savefig(filepath, bbox_inches='tight')



### DAR

min_neglog10_pval = DAR_heatmap_df[DAR_heatmap_df>=0].min().min()
max_neglog10_pval = DAR_heatmap_df[DAR_heatmap_df>=0].max().max()

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

# plot non missing values as heatmap
sns.heatmap(
    data=DAR_heatmap_df,
    annot=DAR_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    mask=DAR_heatmap_mask,
    ax=ax,
    rasterized=False
)

# plot missing values as dark gray_r 
sns.heatmap(
    data=DAR_heatmap_df,
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=missing_cmap,
    cbar=False,
    mask=~DAR_heatmap_mask,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save
filepath = osp.join(results_dir, "ldsc.heatmap.coral.DAR.pdf")
plt.savefig(filepath, bbox_inches='tight')





###### Plot heatmaps (coolwarm)

cmap = sns.color_palette("coolwarm", as_cmap=True)
# cmap = sns.color_palette("light:coral", as_cmap=True)

color = '#9ca09e'
missing_cmap = ListedColormap([color] + plt.get_cmap('Blues')(range(9)).tolist())


### all

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

min_neglog10_pval = all_heatmap_df[all_heatmap_df>=0].min().min()
max_neglog10_pval = all_heatmap_df[all_heatmap_df>=0].max().max()

# plot heatmap 
sns.heatmap(
    data=all_heatmap_df,
    annot=all_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save 
filepath = osp.join(results_dir, "ldsc.heatmap.coolwarm.beforeIPR.pdf")
plt.savefig(filepath, bbox_inches='tight')


### marker

min_neglog10_pval = marker_heatmap_df[marker_heatmap_df>=0].min().min()
max_neglog10_pval = marker_heatmap_df[marker_heatmap_df>=0].max().max()

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

# plot non-missing values as heatmap
sns.heatmap(
    data=marker_heatmap_df,
    annot=marker_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    mask=marker_heatmap_mask,
    ax=ax,
    rasterized=False
)

# plot missing values as dark gray_r 
sns.heatmap(
    data=marker_heatmap_df,
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=missing_cmap,
    cbar=False,
    mask=~marker_heatmap_mask,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save
filepath = osp.join(results_dir, "ldsc.heatmap.coolwarm.marker.pdf")
plt.savefig(filepath, bbox_inches='tight')



### DAR

min_neglog10_pval = DAR_heatmap_df[DAR_heatmap_df>=0].min().min()
max_neglog10_pval = DAR_heatmap_df[DAR_heatmap_df>=0].max().max()

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

# plot non missing values as heatmap
sns.heatmap(
    data=DAR_heatmap_df,
    annot=DAR_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    mask=DAR_heatmap_mask,
    ax=ax,
    rasterized=False
)

# plot missing values as dark gray_r 
sns.heatmap(
    data=DAR_heatmap_df,
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=missing_cmap,
    cbar=False,
    mask=~DAR_heatmap_mask,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save
filepath = osp.join(results_dir, "ldsc.heatmap.coolwarm.DAR.pdf")
plt.savefig(filepath, bbox_inches='tight')



###### Plot heatmaps (custom coral with linear segmented colormap - v1)

# cmap = sns.color_palette("coolwarm", as_cmap=True)
cmap = sns.color_palette("light:coral", as_cmap=True)

final_color = "#f96445"
rgb = colorConverter.to_rgb(final_color)

color = '#9ca09e'
missing_cmap = ListedColormap([color] + plt.get_cmap('Blues')(range(9)).tolist())


### all

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

min_neglog10_pval = all_heatmap_df[all_heatmap_df>=0].min().min()
max_neglog10_pval = all_heatmap_df[all_heatmap_df>=0].max().max()

# sample 15 uniform values from [0, 1]
uniform_values = np.linspace(0, 1, 30)

# get 15 uniform values between [0, rgb[0]]
reds_ = np.linspace(0.95, rgb[0], 30)
greens_ = np.linspace(0.95, rgb[1], 30)
blues_ = np.linspace(0.95, rgb[2], 30)

# generate colormap
reds = [(x, reds_[i], reds_[i]) for i, x in enumerate(uniform_values)]
greens = [(x, greens_[i], greens_[i]) for i, x in enumerate(uniform_values)]
blues = [(x, blues_[i], blues_[i]) for i, x in enumerate(uniform_values)]

# color dct
cdict = {'red': reds, 'green': greens, 'blue': blues}

cmap = LinearSegmentedColormap('custom_cmap', cdict)

# plot heatmap 
sns.heatmap(
    data=all_heatmap_df,
    annot=all_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save 
filepath = osp.join(results_dir, "ldsc.heatmap.coral_custom_v1.beforeIPR.pdf")
plt.savefig(filepath, bbox_inches='tight')


### marker

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

min_neglog10_pval = marker_heatmap_df[marker_heatmap_df>=0].min().min()
max_neglog10_pval = marker_heatmap_df[marker_heatmap_df>=0].max().max()

# sample 15 uniform values from [0, 1]
uniform_values = np.linspace(0, 1, 30)

# get 15 uniform values between [0, rgb[0]]
reds_ = np.linspace(0.95, rgb[0], 30)
greens_ = np.linspace(0.95, rgb[1], 30)
blues_ = np.linspace(0.95, rgb[2], 30)

# generate colormap
reds = [(x, reds_[i], reds_[i]) for i, x in enumerate(uniform_values)]
greens = [(x, greens_[i], greens_[i]) for i, x in enumerate(uniform_values)]
blues = [(x, blues_[i], blues_[i]) for i, x in enumerate(uniform_values)]

# color dct
cdict = {'red': reds, 'green': greens, 'blue': blues}

cmap = LinearSegmentedColormap('custom_cmap', cdict)

# plot non-missing values as heatmap
sns.heatmap(
    data=marker_heatmap_df,
    annot=marker_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    mask=marker_heatmap_mask,
    ax=ax,
    rasterized=False
)

# plot missing values as dark gray_r 
sns.heatmap(
    data=marker_heatmap_df,
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=missing_cmap,
    cbar=False,
    mask=~marker_heatmap_mask,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save
filepath = osp.join(results_dir, "ldsc.heatmap.coral_custom_v1.marker.pdf")
plt.savefig(filepath, bbox_inches='tight')



### DAR

min_neglog10_pval = DAR_heatmap_df[DAR_heatmap_df>=0].min().min()
max_neglog10_pval = DAR_heatmap_df[DAR_heatmap_df>=0].max().max()

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

# sample 15 uniform values from [0, 1]
uniform_values = np.linspace(0, 1, 30)

# get 15 uniform values between [0, rgb[0]]
reds_ = np.linspace(0.95, rgb[0], 30)
greens_ = np.linspace(0.95, rgb[1], 30)
blues_ = np.linspace(0.95, rgb[2], 30)

# generate colormap
reds = [(x, reds_[i], reds_[i]) for i, x in enumerate(uniform_values)]
greens = [(x, greens_[i], greens_[i]) for i, x in enumerate(uniform_values)]
blues = [(x, blues_[i], blues_[i]) for i, x in enumerate(uniform_values)]

# color dct
cdict = {'red': reds, 'green': greens, 'blue': blues}

cmap = LinearSegmentedColormap('custom_cmap', cdict)

# plot non missing values as heatmap
sns.heatmap(
    data=DAR_heatmap_df,
    annot=DAR_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    mask=DAR_heatmap_mask,
    ax=ax,
    rasterized=False
)

# plot missing values as dark gray_r 
sns.heatmap(
    data=DAR_heatmap_df,
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=missing_cmap,
    cbar=False,
    mask=~DAR_heatmap_mask,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save
filepath = osp.join(results_dir, "ldsc.heatmap.coral_custom_v1.DAR.pdf")
plt.savefig(filepath, bbox_inches='tight')




###### Plot heatmaps (custom coral with linear segmented colormap - v2)

# cmap = sns.color_palette("coolwarm", as_cmap=True)
cmap = sns.color_palette("light:coral", as_cmap=True)

final_color = "#f96445"
rgb = colorConverter.to_rgb(final_color)

color = '#9ca09e'
missing_cmap = ListedColormap([color] + plt.get_cmap('Blues')(range(9)).tolist())


### all

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

min_neglog10_pval = all_heatmap_df[all_heatmap_df>=0].min().min()
max_neglog10_pval = all_heatmap_df[all_heatmap_df>=0].max().max()

# sample 15 uniform values from [0, 1]
uniform_values = np.linspace(0, 1, 30)

# define red green and blue break points 
red_break = 0.9582150101419878
green_break = 0.7768762677484786
blue_break = 0.7391480730223123

# get 3 uniform values between [0, break]
# get 27 uniform values between [break, rgb[0]]
reds_ = np.linspace(0.95, red_break, 3).tolist() + np.linspace(red_break, rgb[0], 27).tolist()
greens_ = np.linspace(0.95, green_break, 3).tolist() + np.linspace(green_break, rgb[1], 27).tolist()
blues_ = np.linspace(0.95, blue_break, 3).tolist() + np.linspace(blue_break, rgb[2], 27).tolist()

# convert colors to hex and print 

# hex_colors = [to_hex((x,y,z)) for x,y,z in zip(reds_, greens_, blues_)]
# pprint([(reds_[i], greens_[i], blues_[i], i, color) for i, color in enumerate(hex_colors)])
# pdb.set_trace()

# [(0.95, 0.95, 0.95, 0, '#f2f2f2'),
#  (0.9509127789046653, 0.930764029749831, 0.9265720081135902, 1, '#f2edec'),
#  (0.9518255578093305, 0.9115280594996619, 0.9031440162271804, 2, '#f3e8e6'),
#  (0.9527383367139959, 0.8922920892494929, 0.8797160243407708, 3, '#f3e4e0'),
#  (0.9536511156186612, 0.8730561189993238, 0.856288032454361, 4, '#f3dfda'),
#  (0.9545638945233266, 0.8538201487491548, 0.8328600405679513, 5, '#f3dad4'),
#  (0.9554766734279918, 0.8345841784989858, 0.8094320486815415, 6, '#f4d5ce'),
#  (0.9563894523326572, 0.8153482082488167, 0.7860040567951319, 7, '#f4d0c8'),
#  (0.9573022312373225, 0.7961122379986477, 0.7625760649087221, 8, '#f4cbc2'),
#  (0.9582150101419878, 0.7768762677484786, 0.7391480730223123, 9, '#f4c6bc'),
#  (0.9591277890466531, 0.7576402974983096, 0.7157200811359026, 10, '#f5c1b7'),
#  (0.9600405679513184, 0.7384043272481406, 0.6922920892494928, 11, '#f5bcb1'),
#  (0.9609533468559838, 0.7191683569979715, 0.668864097363083, 12, '#f5b7ab'),
#  (0.961866125760649, 0.6999323867478024, 0.6454361054766734, 13, '#f5b2a5'),
#  (0.9627789046653144, 0.6806964164976335, 0.6220081135902636, 14, '#f6ae9f'),
#  (0.9636916835699797, 0.6614604462474645, 0.598580121703854, 15, '#f6a999'),
#  (0.964604462474645, 0.6422244759972955, 0.5751521298174442, 16, '#f6a493'),
#  (0.9655172413793103, 0.6229885057471264, 0.5517241379310345, 17, '#f69f8d'),
#  (0.9664300202839756, 0.6037525354969573, 0.5282961460446247, 18, '#f69a87'),
#  (0.967342799188641, 0.5845165652467883, 0.504868154158215, 19, '#f79581'),
#  (0.9682555780933062, 0.5652805949966193, 0.48144016227180525, 20, '#f7907b'),
#  (0.9691683569979715, 0.5460446247464503, 0.4580121703853955, 21, '#f78b75'),
#  (0.9700811359026369, 0.5268086544962812, 0.4345841784989858, 22, '#f7866f'),
#  (0.9709939148073022, 0.5075726842461121, 0.411156186612576, 23, '#f88169'),
#  (0.9719066937119675, 0.48833671399594314, 0.38772819472616626, 24, '#f87d63'),
#  (0.9728194726166328, 0.46910074374577415, 0.3643002028397565, 25, '#f8785d'),
#  (0.9737322515212982, 0.44986477349560505, 0.34087221095334685, 26, '#f87357'),
#  (0.9746450304259635, 0.43062880324543606, 0.3174442190669371, 27, '#f96e51'),
#  (0.9755578093306287, 0.4113928329952671, 0.29401622718052733, 28, '#f9694b'),
#  (0.9764705882352941, 0.39215686274509803, 0.27058823529411763, 29, '#f96445')]


# generate colormap
reds = [(x, reds_[i], reds_[i]) for i, x in enumerate(uniform_values)]
greens = [(x, greens_[i], greens_[i]) for i, x in enumerate(uniform_values)]
blues = [(x, blues_[i], blues_[i]) for i, x in enumerate(uniform_values)]

# color dct
cdict = {'red': reds, 'green': greens, 'blue': blues}

cmap = LinearSegmentedColormap('custom_cmap', cdict)

# plot heatmap 
sns.heatmap(
    data=all_heatmap_df,
    annot=all_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save 
filepath = osp.join(results_dir, "ldsc.heatmap.coral_custom_v2.beforeIPR.pdf")
plt.savefig(filepath, bbox_inches='tight')


### marker

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

min_neglog10_pval = marker_heatmap_df[marker_heatmap_df>=0].min().min()
max_neglog10_pval = marker_heatmap_df[marker_heatmap_df>=0].max().max()

# sample 15 uniform values from [0, 1]
uniform_values = np.linspace(0, 1, 30)

# define red green and blue break points 
red_break = 0.9582150101419878
green_break = 0.7768762677484786
blue_break = 0.7391480730223123

# get 3 uniform values between [0, break]
# get 27 uniform values between [break, rgb[0]]
reds_ = np.linspace(0.95, red_break, 3).tolist() + np.linspace(red_break, rgb[0], 27).tolist()
greens_ = np.linspace(0.95, green_break, 3).tolist() + np.linspace(green_break, rgb[1], 27).tolist()
blues_ = np.linspace(0.95, blue_break, 3).tolist() + np.linspace(blue_break, rgb[2], 27).tolist()

# generate colormap
reds = [(x, reds_[i], reds_[i]) for i, x in enumerate(uniform_values)]
greens = [(x, greens_[i], greens_[i]) for i, x in enumerate(uniform_values)]
blues = [(x, blues_[i], blues_[i]) for i, x in enumerate(uniform_values)]

# color dct
cdict = {'red': reds, 'green': greens, 'blue': blues}

cmap = LinearSegmentedColormap('custom_cmap', cdict)

# plot non-missing values as heatmap
sns.heatmap(
    data=marker_heatmap_df,
    annot=marker_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    mask=marker_heatmap_mask,
    ax=ax,
    rasterized=False
)

# plot missing values as dark gray_r 
sns.heatmap(
    data=marker_heatmap_df,
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=missing_cmap,
    cbar=False,
    mask=~marker_heatmap_mask,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save
filepath = osp.join(results_dir, "ldsc.heatmap.coral_custom_v2.marker.pdf")
plt.savefig(filepath, bbox_inches='tight')



### DAR

min_neglog10_pval = DAR_heatmap_df[DAR_heatmap_df>=0].min().min()
max_neglog10_pval = DAR_heatmap_df[DAR_heatmap_df>=0].max().max()

# create a figure
fig, ax = plt.subplots(figsize=(10, 3))

# sample 15 uniform values from [0, 1]
uniform_values = np.linspace(0, 1, 30)

# define red green and blue break points 
red_break = 0.9582150101419878
green_break = 0.7768762677484786
blue_break = 0.7391480730223123

# get 3 uniform values between [0, break]
# get 27 uniform values between [break, rgb[0]]
reds_ = np.linspace(0.95, red_break, 3).tolist() + np.linspace(red_break, rgb[0], 27).tolist()
greens_ = np.linspace(0.95, green_break, 3).tolist() + np.linspace(green_break, rgb[1], 27).tolist()
blues_ = np.linspace(0.95, blue_break, 3).tolist() + np.linspace(blue_break, rgb[2], 27).tolist()

# generate colormap
reds = [(x, reds_[i], reds_[i]) for i, x in enumerate(uniform_values)]
greens = [(x, greens_[i], greens_[i]) for i, x in enumerate(uniform_values)]
blues = [(x, blues_[i], blues_[i]) for i, x in enumerate(uniform_values)]

# color dct
cdict = {'red': reds, 'green': greens, 'blue': blues}

cmap = LinearSegmentedColormap('custom_cmap', cdict)

# plot non missing values as heatmap
sns.heatmap(
    data=DAR_heatmap_df,
    annot=DAR_annots,
    fmt="s",
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=cmap,
    vmin=min_neglog10_pval, 
    vmax=max_neglog10_pval,
    cbar=True,
    mask=DAR_heatmap_mask,
    ax=ax,
    rasterized=False
)

# plot missing values as dark gray_r 
sns.heatmap(
    data=DAR_heatmap_df,
    square=True,
    linewidths=0.5,
    linecolor='black',
    cmap=missing_cmap,
    cbar=False,
    mask=~DAR_heatmap_mask,
    ax=ax,
    rasterized=False
)

ax.set_facecolor('white')

plt.tight_layout()

# save
filepath = osp.join(results_dir, "ldsc.heatmap.coral_custom_v2.DAR.pdf")
plt.savefig(filepath, bbox_inches='tight')



print("Script finished")
