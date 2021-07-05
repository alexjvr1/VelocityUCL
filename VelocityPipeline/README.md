# Pipeline for Velocity project: From raw data to filtered variants. 

Use this pipeline by running scripts in the pipeline/ folder in the numbered order. These scripts generate submission scripts for 
each step that can be submitted to the server queue (i.e. qsub script.sh). 

We are using the University of Bristol BlueCrystal p3 server and the UCL pchuckle server

Inputs for each step should be submitted via the command line. 

Folders: 

1. pipeline #Scripts for generating submission scripts for each step. Inputs can be specified in the command line. 
2. wrapper #Scripts called by pipeline scripts. Wrapper scripts for generating Queue request and specifying inputs from command line for the tools to be called. 
3. tools #Scripts of tools or functions used in each step. These are called by the wrapper script

