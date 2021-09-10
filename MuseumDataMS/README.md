# Museum Lepidoptera samples

Aim: Pipeline to use museum Lep samples

Cross-species comparison of genomic data quality using lab protocol (designed by Carl) and Bioinformatic pipeline. 

Questions

1) What factors predict data quality (sample quality, sample preservation agent (e.g. dessication vs chemical), DNA concentrations.. what else?)


## Initial check of data qualty vs DNA extraction concentrations for E3 Aphantopus hyperantus

Worksheet and graphs on mac here: 

/Users/alexjvr/2021postdoc/MuseumDataQuality


```
E3 <- read.table("E3_data", header=T)

pdf("E3_Museumdata.figs.pdf")
ggplot(E3, aes(x=DNAExtraction_concentration.ng.ul., y=RawReadsTotal, color=Pop2))+ geom_point() + ggtitle("DNA conc vs total raw reads") + labs(x="DNA extraction (ng/ul)", y="Total raw reads")

ggplot(E3, aes(x=DNAExtraction_concentration.ng.ul., y=Pairs.that.were.too.short...20bp.min.length., color=Pop2))+ geom_point() + ggtitle("DNA conc vs nr short reads filtered") + labs(x="DNA extraction (ng/ul)", y="Nr Filtered Read Pairs (< 20bp)")

ggplot(E3, aes(x=DNAExtraction_concentration.ng.ul., y=mapped_merged.perc, color=Pop2))+ geom_point() +ggtitle("DNA conc vs % merged read pairs") + labs(x="DNA extraction (ng/ul)", y="% Merged Read Pairs (ATLAS)")

dev.off()
```

