# Recal tests

We've run into a problem with the recal module in ATLAS. Although it runs well on the modern samples, we are running into three main errors when using the museum samples: 

Two errors are returned when runs are stopped early: 

```
1. double free or corruption (!prev): 0x00000003f5ee7f00 ***

2. line 31: 251665 Segmentation fault 
```


A second issue is that although some runs complete, they produce unusable results: 

```
readGroup mate model quality position context

E3mus first qualFuncPosFuncContext 1.000000,0.000000 0.000000,0.000000 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000

E3mus second qualFuncPosFuncContext 1.149249,2.305794 0.423096,-0.676807 8.685975,-0.011506,0.030326,-3.649817,-14.060771,0.022166,-0.000305,-2.359712,4.599552,0.007064,0.025039,-3.001226,6.238543,0.020372,0.014092,-3.002097,-1.474109,0.001849,-0.044600,-0.315406



readGroup mate model quality position context

E3mus first qualFuncPosFuncContext 1.000000,0.000000 0.000000,0.000000 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000

E3mus second qualFuncPosFuncContext 1.000000,0.000000 0.000000,0.000000 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000
```
Vivian thought the runs might not be initialising properly. These outputs seem unrelated to the amount of data available as they are produced with pooled data as well. 



Based on discussions with Vivian and Daniel, we've decided to test pooling samples from a run. We need to test the following: 

1. Is the pooled recal comparable to the mean individually estimated recal within species? 

2. Do runs work if samples that ran to completion as individuals are pooled? (ie we can set up the pipeline to identify poolable samples rather than using trial and error)

3. How does variance in recal estimates within population compare between museum and modern populations of the same species? 

4. How does variance in recal compare between species? 


## 1. Is the pooled recal comparable to the mean individually estimated recal within species? 

I'm estimating pooled and mean recal from the modern populations for C3, D3 and E3. 

Individual runs are submitted as an array using modifications of the script: [04_ATLAS_LR.recal_ARRAY.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/04_ATLAS_LR.recal_ARRAY.sh)

Subsets of 5 and 10 individuals where I obtained usable recal estimates are pooled to estimate a pooled recal. I limited the number of individuals pooled in the modern dataset because using too much data either stops the run prematurely or the run takes a very long time to complete.

Script for pooled samples:[04b_ATLAS_recal.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/04b_ATLAS_recal.sh) 

When pooling samples we have to generate a pooled pmd file by concatenating all the pmd outputs together. Change the RGs if necessary. 

We also need to provide a text file containing all the RGs to be pooled for the analysis: [MergeTheseRGs.txt](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/MergeTheseRGs.txt)


### Results: 

Recal outputs a line of numbers for each RG and each read mate (first and second). To plot these we do the following: 

#1. Concatenate all the individual outputs together
```
cat PA*EM.txt > D3_concat_EM.txt
```

#2. Add the pooled estimate at the end of the concatenated file
```
#Check that this adds a header line + the recal scores for one sample. 
#Recal scores are printed out for each RG in the merged recal file, but they are the same for each RG. 
head -n 3 D3_MODC.merged_mergedReads_recalibrationEM.txt >> D3_concat_EM.txt

#Remove all header lines
sed -i '/^readGroup/d' D3_concat_EM.txt

#replace all commas with white space
sed -i 's:,: :g' D3_concat_EM.txt
```


#3. Download to the local computer

```
/Users/alexjvr/2022ResearchAssoc/Velocity/Velocity/ATLAS/recal.tests
```

#4. Plot in R
```
library(ggplot2)
library(tidyr)


## Import
D3.data <- read.table("D3.concat_and_pooled_MODC.5samples.txt", header=F)


## Check the file looks as expected
D3.data
      V1     V2                     V3       V4       V5       V6        V7
1 D3modc  first qualFuncPosFuncContext 1.205653 3.963455 0.987101 -0.563430
2 D3modc second qualFuncPosFuncContext 0.894571 1.017038 1.048818 -1.006856
3 D3modc  first qualFuncPosFuncContext 1.075654 2.965052 0.250877 -0.064447
4 D3modc second qualFuncPosFuncContext 0.963502 1.654850 1.454248 -1.261661
5 D3modc  first qualFuncPosFuncContext 1.198500 3.970225 0.789103 -0.446626
6 D3modc second qualFuncPosFuncContext 0.933076 1.407542 1.085903 -1.018628
7 pooled  first qualFuncPosFuncContext 1.233400 4.818856 0.566271 -0.336367
8 pooled second qualFuncPosFuncContext 1.035614 2.774342 1.150683 -1.000124
        V8        V9       V10      V11      V12       V13       V14      V15
1 1.546000 -1.256210 -0.465360 1.916629 2.162605 -1.876201 -0.520749 2.182440
2 1.554150 -1.208768 -0.949921 1.799370 2.158892 -1.843451 -1.014033 1.921691
3 1.366367 -2.035269 -0.802269 1.927776 2.249168 -1.678366 -0.774877 1.960602
4 1.631323 -1.109913 -0.788684 1.999599 2.352351 -1.916735 -0.812680 2.175060
5 1.512194 -1.362608 -0.682871 2.059782 2.273671 -1.706140 -0.411752 2.037303
6 1.459185 -1.263943 -0.764046 1.877952 2.280425 -1.622745 -0.954507 1.963918
7 1.778585 -1.074183 -0.500317 2.170472 2.343943 -1.551735 -0.460030 2.289950
8 1.843988 -0.966976 -0.729365 2.043402 2.393554 -1.417355 -0.800932 2.194779
       V16       V17       V18      V19      V20       V21       V22      V23
1 1.982016 -1.445051 -0.736458 1.757839 2.051968 -1.539481 -0.533273 1.896688
2 2.183993 -1.309434 -1.308306 1.615752 2.088104 -1.339623 -1.014109 1.774497
3 1.815404 -1.786892 -0.693061 1.690254 1.979666 -1.907492 -0.829189 1.838229
4 2.276607 -1.259518 -1.477554 1.758469 2.236132 -1.215711 -1.229573 1.771273
5 2.062245 -1.374439 -0.761343 1.734326 2.075187 -1.630160 -0.427580 1.995777
6 2.209132 -1.181627 -1.223465 1.747409 2.121742 -1.341757 -0.972142 1.869105
7 2.189389 -1.383429 -0.742299 1.984686 2.206470 -1.298542 -0.485930 2.171634
8 2.346280 -1.242370 -1.162647 1.906477 2.311779 -1.121378 -0.844369 2.056696
        V24        V25        V26        V27
1  1.554695  -2.351347  -1.331290   2.082627
2 -0.383266 -13.942028 -13.946869   1.529104
3  1.898253 -13.374007  -1.430130   0.653675
4  1.589262  -1.687660  -4.065619   2.015323
5  2.486214 -13.160077 -13.482991 -10.765245
6  1.565993  -2.052202 -13.510992   2.004051
7  1.987764  -1.538292  -1.778104   1.884101
8  1.198569  -2.419262  -1.937229   2.169192


## Add informative column headers
colnames(D3.data) <- c("readGroup", "mate", "model", "quality1", "quality2", "position1", "position2", "context1", "context2", "context3", "context4", "context5", "context6", "context7", "context8", "context9", "context10", "context11", "context12", "context13", "context14", "context15", "context16", "context17", "context18", "context19", "context20")


## Pivot to long format
D3.data_long <- gather(D3.data, variable, number, quality1:quality2, position1:position2, context1:context20, factor_key=T)

D3.data_long
    readGroup   mate                  model  variable     number
1      D3modc  first qualFuncPosFuncContext  quality1   1.205653
2      D3modc second qualFuncPosFuncContext  quality1   0.894571
3      D3modc  first qualFuncPosFuncContext  quality1   1.075654
4      D3modc second qualFuncPosFuncContext  quality1   0.963502
5      D3modc  first qualFuncPosFuncContext  quality1   1.198500
6      D3modc second qualFuncPosFuncContext  quality1   0.933076
7      pooled  first qualFuncPosFuncContext  quality1   1.233400

## Plot
## If needed subset to exclude the last three variables as they show high variance. 
pdf("D3_MODC.recal_variance.pdf")
ggplot(D3.data_long[1:1680,], aes(x=variable, y=number, shape=mate))+geom_boxplot()+geom_point(position=position_dodge(width=0.75), aes(colour=factor(readGroup)))
dev.off()
```


#### C3 modc
 
 
 
#### C3 mus



#### D3 modc

D3 pool based on the first 5 samples that worked. I tried to pool 10 samples but this run stopped with an error. 

Results: 

D3_MODC.recal_variance.pdf


![alt_txt][D3.modc.EM]

[D3.modc.EM]:https://user-images.githubusercontent.com/12142475/219612271-674402bf-ffb5-4d79-9572-5f052e6629c1.png




#### D3 mus





#### E3 modc 

E3 modc was run across several lanes. Thus we have the opportunity to compare variance between runs within a modern population. 

Results: 

E3_MODC_RG1_recal_variance.pdf

E3_MODC_RG2_recal_variance.pdf


##### RG1
![alt_txt][RG1_MODC]

[RG1_MODC]:https://user-images.githubusercontent.com/12142475/219434612-c4c4691d-b9c6-418e-9214-a157e39794e1.png


##### RG2
![alt_txt][RG2_MODC]

[RG2_MODC]: https://user-images.githubusercontent.com/12142475/219437583-a08e57dc-82c6-40d2-9119-74e6e56de567.png


#### E3 mus

Tried to pool various combinations of samples that worked individually but runs either stop without warning or yield 1.000.. results. 

Plotted the results of the samples that did worked individually so that we can get an idea of the variance in recal across samples: 

E3.mus.recal_variance.pdf		

E3.mus.recal_variance_excl_outliers.pdf


##### Full plot

![alt_txt][E3.mus.plot1]

[E3.mus.plot1]:https://user-images.githubusercontent.com/12142475/219616406-9ca03652-ff4a-41dd-959e-657a760a35bb.png


Oultiers found in second mate of these samples: 

AH-01-1900-17.RG1_recalibrationEM.txt

AH-01-1900-18.RG1_recalibrationEM.txt

AH-01-1900-26.RG1_recalibrationEM.txt

AH-01-1900-44.RG1_recalibrationEM.txt



##### Outliers removed

![alt_txt][E3.mus.plot2]

[E3.mus.plot2]:https://user-images.githubusercontent.com/12142475/219616368-96c9b90b-083a-4d84-9493-2f3c8a7bab4e.png
