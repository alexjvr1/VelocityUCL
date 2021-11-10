# PE read merging

After plotting MapDamage profiles, I've found a strange distribution of G>A substitutions in the Ringlet data after merging PE reads. This pattern holds in the Brown Argus and Pararge aegeria museum data as well (run by Frances and Zhirou). 

![alt_txt][Fig1]

[Fig1]:https://user-images.githubusercontent.com/12142475/141117504-0e913f6f-fe53-4ccb-960f-d9c2220edf5b.png


Possible issues: 

1) Adapters haven't been filtered out properly

2) As are kept preferentially over Gs based on base quality during merging. Is there a bias in base quality? Are we more certain of As in general? 

3) Reads are merged before any error correction (e.g. filtering for low quality reads, or removing the first 10bp of the reads. 



