# Format scATAC-seq peaks for calculating LD scores


## Format marker peaks
```
for i in MarkerSubPeaks/*.tsv; do cut -f 1,3,4 $i | sed 1d | awk '{print $0"\t"$1"-"$2"-"$3}' >  hg38/*.bed ; done
```

## Format all peaks
```
for i in Ast End ExN InN Mic Oli OPC; do Rscript get_peaks.R BroadPeaks_without_constraint/$i\-reproduciblePeaks.gr.rds 250 hg38/$i\.bed; done
while read a; do Rscript get_peaks.R SubClusterPeaks_without_constraint/$a\-reproduciblePeaks.gr.rds 250 hg38/$a\.bed; done < sub.lst 
```

## Convert peak coordinates to hg19
```
for i in hg38/*.bed; do ~/utils/liftOver $i ~/utils/hg38ToHg19.over.chain.gz $(echo $i | sed "s/38/19/g") unmapped ; done
cat hg19/*_markerpeaks.tsv.bed | sort | uniq > hg19/sub_markerpeaks.tsv.bed
for i in Ast End ExN InN Mic Oli OPC; do cat hg19/$i.bed ; done > hg19/broad.bed
while read a; do cat hg19/$a.bed ; done < sub.lst > hg19/sub.bed
```
