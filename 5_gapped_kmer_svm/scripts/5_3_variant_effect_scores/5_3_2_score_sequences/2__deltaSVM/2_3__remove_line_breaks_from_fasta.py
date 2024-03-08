import pdb
import os.path as osp


# remove line breaks in each record of the fasta file
def remove_line_breaks_from_fasta(infile_path, outfile_path):
    with open(infile_path, "r") as f, open (outfile_path, "w") as output:
        lines = f.readlines()
        for i in range(len(lines)-1):
            # start line a new record
            if lines[i].startswith(">"):
                continue
            # end line of current record
            elif lines[i+1].startswith(">"):
                continue
            # intermediate lines of the current record
            else:
                lines[i] = lines[i].strip() # remove the trailing line break
        # write to the output file
        for line in lines:
            output.write(line)

seq_len = "201bp"

# analysis_type = "subcluster"
analysis_type = "broad"

# cell types to analyze
cell_types = ["Ast1", "Ast2", "Ast3", "Ast4", "End2", "ExN1", \
    "ExN1_L23", "ExN1_L24", "ExN1_L46", "ExN1_L56", "ExN2", "ExN2_L23", \
    "ExN2_L46", "ExN2_L56", "ExN3_L46", "ExN3_L56", "ExN4_L56", "In_LAMP5", \
    "InN3", "In_PV", "In_SST", "In_VIP", "Mic1", "Mic2", "End1", "Oli1", "Oli2", \
    "Oli3", "Oli4", "Oli5", "Oli6", "Oli7", "OPC1", "OPC2", "OPC3", "OPC4"]

broad2sub = {
    "Ast": ["Ast1", "Ast2", "Ast3", "Ast4"],
    "End": ["End2"],
    "ExN": ["ExN1", "ExN1_L23", "ExN1_L24", "ExN1_L46", "ExN1_L56", "ExN2", "ExN2_L23", "ExN2_L46", "ExN2_L56", "ExN3_L46", "ExN3_L56", "ExN4_L56"],
    "InN": ["In_LAMP5", "InN3", "In_PV", "In_SST", "In_VIP"],
    "Mic": ["Mic1", "Mic2", "End1"],
    "Oli": ["Oli1", "Oli2", "Oli3", "Oli4", "Oli5", "Oli6", "Oli7"],
    "OPC": ["OPC1", "OPC2", "OPC3", "OPC4"]
}

broad_clusters = []
for x in cell_types:
    if x in broad2sub["Ast"]:
        broad_clusters.append("Ast")
    elif x in broad2sub["End"]:
        broad_clusters.append("End")
    elif x in broad2sub["ExN"]:
        broad_clusters.append("ExN")
    elif x in broad2sub["InN"]:
        broad_clusters.append("InN")
    elif x in broad2sub["Mic"]:
        broad_clusters.append("Mic")
    elif x in broad2sub["Oli"]:
        broad_clusters.append("Oli")
    elif x in broad2sub["OPC"]:
        broad_clusters.append("OPC")
    else:
        raise ValueError("Cell type {} not found in broad2sub dictionary".format(x))

# SNP sequence path 
seq_base_path = f"/home/dcakma3/scratch/mdd-prepare_snps_in_peaks/{analysis_type}"

source_seq_filenames = []
for cell_type, broad_cluster in zip(cell_types, broad_clusters):
    source_seq_filenames.append(
        osp.join(seq_base_path, f"{broad_cluster}/{cell_type}/{seq_len}/snps.hg38_information.{seq_len}.a1"),
    )
    source_seq_filenames.append(
        osp.join(seq_base_path, f"{broad_cluster}/{cell_type}/{seq_len}/snps.hg38_information.{seq_len}.a2"),
    )
    source_seq_filenames.append(
        osp.join(seq_base_path, f"{broad_cluster}/{cell_type}/{seq_len}/snps.hg38_information.{seq_len}.a1.shuffled"),
    )
    source_seq_filenames.append(
        osp.join(seq_base_path, f"{broad_cluster}/{cell_type}/{seq_len}/snps.hg38_information.{seq_len}.a2.shuffled"),
    )

target_seq_filenames = [f"{x}.no_lb.fa" for x in source_seq_filenames]
source_seq_filenames = [f"{x}.fa" for x in source_seq_filenames]

for source, target in zip(source_seq_filenames, target_seq_filenames):
    remove_line_breaks_from_fasta(source, target)
