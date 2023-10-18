import pandas as pd 


coloring = {
    "Ast1":"#23B205",
    "Ast2":"#74C69D",
    "Ast3":"#2DF304",
    "Ast4":"#145F04",
    "End1":"#E3D2A7",
    "End2":"#F7E3B4",
    "ExN1":"#E7B321",
    "ExN1_L23":"#A6A3A1",
    "ExN1_L24":"#CCC7C4",
    "ExN1_L46":"#D0C75C",
    "ExN1_L56":"#FAB866",
    "ExN2":"#FAE605",
    "ExN2_L23":"#736F6C",
    "ExN2_L46":"#F7F7B5",
    "ExN2_L56":"#F29F3A",
    "ExN3_L46":"#E7D973",
    "ExN3_L56":"#D18221",
    "ExN4_L56":"#B36609",
    "In_LAMP5":"#53B1F3",
    "In_PV":"#4853A2",
    "In_SST":"#6E7EBD",
    "In_VIP":"#81B2F1",
    "InN3":"#C0D9FC",
    "Mic1":"#1792AA",
    "Mic2":"#83D6E6",
    "Mix1":"#D2D6CB",
    "Mix2":"#D7C7D5",
    "Oli1":"#C60E9D",
    "Oli2":"#B199C7",
    "Oli3":"#C3B2D2",
    "Oli4":"#7F56A4",
    "Oli5":"#E04EC0",
    "Oli6":"#E677CD",
    "Oli7":"#F4BBE7",
    "OPC1":"#F4978E",
    "OPC2":"#FFB4A2",
    "OPC3":"#F8651B",
    "OPC4":"#F08080"
}
subcluster_coloring = coloring



broad2sub = {
    "Ast": ["Ast1", "Ast2", "Ast3", "Ast4"],
    "End": ["End1", "End2"],
    "ExN": ["ExN1", "ExN2", "ExN3", "ExN4", "ExN5", "ExN6", "ExN7", "ExN8", "ExN9", "ExN10", "ExN11", "ExN12"],
    "InN": ["InLAMP5", "InN3", "InPV", "InSST", "InVIP"],
    "Mic": ["Mic1", "Mic2"],
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

sub2broad = {}
for k, v in broad2sub.items():
    for x in v:
        sub2broad[x] = k
        
