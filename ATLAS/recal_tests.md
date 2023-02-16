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
