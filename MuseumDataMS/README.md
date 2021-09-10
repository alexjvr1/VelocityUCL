# Museum Lepidoptera samples

Aim: Pipeline to use museum Lep samples

papers [here](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1755-0998.2011.03052.x?casa_token=YTs1POVOPK4AAAAA%3A7XFOKFfBpfpIGMPdhKVZlFNkla5EuUg40XKAgqg4FAVCpilKt1qG6QfobST675at_brqhhiOklqQ6w); [here](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12516?casa_token=yegfoQbNAncAAAAA%3AUxANb32iov7KiQJzBPu6edGendNDI6EfJ3hBl1-wPiz7hrQOaZf7zznOH9uXmDAnM2sq1gdP34HE0w): [here](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6594-0); [here](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.3065); [here](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13082?casa_token=SwEzkfoTKcoAAAAA%3ANY3ox93WiWfk8JwOGyYygbh1fViES9ReVP3c55BMrqQ9OoxjvZV3Bwle1AIu8jGTuEixneuCfc_sgA); [here](https://onlinelibrary.wiley.com/doi/full/10.1111/syen.12481): [here](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13245?casa_token=NNHXxQWaA7IAAAAA%3AbJtwuQJOVentokJAraPq9AQzYGSiUiT2GmZ3KzLdjv-zc9NW295hciVOWImnFZ7nqpzi6aROWZw-lA); [here](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13269); 

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





