# Museum Lepidoptera samples

Aim: Pipeline to use museum Lep samples

Cross-species comparison of genomic data quality using lab protocol (designed by Carl) and Bioinformatic pipeline. 

Questions

1) What factors predict data quality (sample quality, sample preservation agent (e.g. dessication vs chemical), DNA concentrations.. what else?)

2) Does the PMD profile change vs age and preservation method? 

3) Does PMD change He and nuc div? 

4) Does depth, sequence length or insert size affect estimates of He and nucleotide diversity? (within species) 

5) Can we efficiently compare between modern and museum data? Do we have to make any corrections? e.g. subsample to the same depth

expectation 1) Based on prelim data it looks like we have good museum data (nr of reads) but short reads. So less of the genome is covered compared to the modern data. 


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

![alt_txt][Fig1]

[Fig1]:https://user-images.githubusercontent.com/12142475/132857971-b267556e-9d9d-4f07-a182-cd2e09c93449.png

![alt_txt][Fig2]

[Fig2]:https://user-images.githubusercontent.com/12142475/132857979-6accd49b-8bf9-4b6a-a0c1-b4c9b34e123f.png

![alt_txt][Fig3]

[Fig3]:https://user-images.githubusercontent.com/12142475/132857982-463e57f2-1b61-4625-9df5-846a8c19303b.png





