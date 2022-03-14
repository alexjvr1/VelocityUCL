# diploS/HIC

To identify signals of selection between populations we'll use the deep learning neural network in diploS/HIC. This was originally developed to differentiate

hard and soft selection in a single population. 

diploSHIC paper [here](https://academic.oup.com/g3journal/article/8/6/1959/6028059) 

[GitHub of diploS/HIC](https://github.com/kr-colab/diploSHIC)

[GitHub of demographic simulator discoal](https://github.com/kr-colab/discoal)

We'll use a modification of this method as implemented by [Kaut et al. 2020](https://www.nature.com/articles/s41586-020-2845-0?fbclid=IwAR3gk7HEGmPn5V0giczpdMycAgpu-Xttr8_cD550VPuW8tdQH6KSCYN_e60) on pairwise cichlid populations. 


# Set up work environment

We need three python3 packages for working with machine learning algorithms and genomic data: 
```
export PATH=/share/apps/python-3.6.4-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.6.4-shared/lib:$LD_LIBRARY_PATH

python3 -m pip install scikit-allel tensorflow keras --user
```

And install diploSHIC
```
git clone https://github.com/kern-lab/diploSHIC.git
cd diploSHIC
source /share/apps/source_files/python/python-3.9.5.source
python3 -m pip install scikit-allel tensorflow keras --user
python3 setup.py install --user

Then, run:

python3 diploSHIC.py -h

to verify the install went smoothly.
```


# Initial tests on a single population

We're testing the method on C3_Aricia_agestis MODC

We're assuming a constant population size over the last 150 years. Discoal can simulate demographic events and changes in population size as steps (not continuous). 


## Workflow

1. Test most likely demographic history with FastSimCoal

2. Simulate data in windows for Neutral, Hard selection and Soft selection

3. Generate feature vectors for all the simulations

4. Train and test the model

5. Generate feature vectors for the real data

6. Run the trained model on the real data. 



## Use diploSHIC

To run diploSHIC we need to install some python3 packages. The machine learning model is built in tensorflow using keras: 
```
export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH

python3 -m pip install scikit-allel tensorflow keras --user
```

```
export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH

diploSHIC=/SAN/ugi/LepGenomics/Software/diploSHIC/diploshic/diploSHIC

python3 $diploSHIC -h
```

