# Plot SFS

```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_Dec2021/03.1_SAF_0.8GenotypingRate

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
