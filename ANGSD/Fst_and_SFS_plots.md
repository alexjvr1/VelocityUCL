# Plot SFS

```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_Dec2021/03.1_SAF_0.8GenotypingRate

#MODC
MODC.sfs <- read.table("MODC.LR75.sfs", header=F)
library(reshape2)
colnames(MODC.sfs) <- seq(1:length(MODC.sfs))
MODC.sfs
MODC.sfs2 <- melt(MODC.sfs)
MODC.sfs2
library(ggplot2)
pdf("MODC.sfs.pdf")
ggplot(MODC.sfs2, aes(x=variable, y=value))+geom_point()
ggplot(MODC.sfs2[-1,], aes(x=variable, y=value))+geom_point()
dev.off()


#MODE
MODE.sfs <- read.table("MODE.LR75.sfs", header=F)
library(reshape2)
colnames(MODE.sfs) <- seq(1:length(MODE.sfs))
MODE.sfs
MODE.sfs2 <- melt(MODE.sfs)
MODE.sfs2
library(ggplot2)
pdf("MODE.sfs.pdf")
ggplot(MODE.sfs2, aes(x=variable, y=value))+geom_point()
ggplot(MODE.sfs2[-1,], aes(x=variable, y=value))+geom_point()
dev.off()

#MUS
MUS.sfs <- read.table("MUS.LR75.sfs", header=F)
library(reshape2)
colnames(MUS.sfs) <- seq(1:length(MUS.sfs))
MUS.sfs
MUS.sfs2 <- melt(MUS.sfs)
MUS.sfs2
library(ggplot2)
pdf("MUS.sfs.pdf")
ggplot(MUS.sfs2, aes(x=variable, y=value))+geom_point()
ggplot(MUS.sfs2[-1,], aes(x=variable, y=value))+geom_point()
dev.off()

```



# Calculate Fst

```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_Dec2021/SFS_0.8GenotypingRate

realSFS=/share/apps/genomics/angsd-0.935/bin/realSFS

#MODC-MODE
~/bin/angsd/misc/realSFS fst index MODC.saf.idx MODE.saf.idx -sfs MODC.MODE.folded.sfs -fstout MODC.MODE.fstout

$realSFS fst stats2 MODC.MODE.fstout.fst.idx -win 1000000 -step 250000 > slidingwindow.MODC.MODE_win1M.step0.25M.fst
$realSFS fst stats2 MODC.MODE.fstout.fst.idx -win 1000000 -step 1000000 > slidingwindow.MODC.MODE_win1M.step1M.fst


#MODC-MUS
~/bin/angsd/misc/realSFS fst index MODC.saf.idx MUS.saf.idx -sfs MODC.MUS.folded.sfs -fstout MODC.MUS.fstout

$realSFS fst stats2 MODC.MUS.fstout.fst.idx -win 1000000 -step 250000 > slidingwindow.MODC.MUS_win1M.step0.25M.fst
$realSFS fst stats2 MODC.MUS.fstout.fst.idx -win 1000000 -step 1000000 > slidingwindow.MODC.MUS_win1M.step1M.fst



#MODE-MUS
~/bin/angsd/misc/realSFS fst index MODE.saf.idx MUS.saf.idx -sfs MODE.MUS.folded.sfs -fstout MODE.MUS.fstout

$realSFS fst stats2 MODE.MUS.fstout.fst.idx -win 1000000 -step 250000 > slidingwindow.MODE.MUS_win1M.step0.25M.fst
$realSFS fst stats2 MODE.MUS.fstout.fst.idx -win 1000000 -step 1000000 > slidingwindow.MODE.MUS_win1M.step1M.fst


```


Global Fst
```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_Dec2021/SFS_0.8GenotypingRate

#MODC-MODE
$realSFS fst stats MODC.MODE.fstout.fst.idx 
	-> Assuming idxname:MODC.MODE.fstout.fst.idx
	-> Assuming .fst.gz file: MODC.MODE.fstout.fst.gz
	-> FST.Unweight[nObs:4171689]:0.009984 Fst.Weight:0.095347
0.009984	0.095347

#MODC-MUS


#MODE-MUS

```

Diversity Estimates 
```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_Dec2021/03.1_SAF_0.8GenotypingRate

##MODE
$realSFS MODE.LR761675.1.1.SEPT20.saf.idx -fold 1 > MODE.LR75.fold.sfs
$realSFS saf2theta MODE.LR761675.1.1.SEPT20.saf.idx -sfs MODE.LR75.fold.sfs -outname MODE.LR75.fold.thetas
 
##MODC
$realSFS MODC.LR761675.1.1.SEPT20.saf.idx -fold 1 > MODC.LR75.fold.sfs
$realSFS saf2theta MODC.LR761675.1.1.SEPT20.saf.idx -sfs MODC.LR75.fold.sfs -outname MODC.LR75.fold.thetas


```
