# Data Transfer

Transfer data from UoB backups to UCL shared folder


1. Log onto UoB bluecrystalp3 (which has access to the projects folders)

2. Set up port forwarding to UCL server
```
ssh -l ajansen -L 3000:pchuckle.cs.ucl.ac.uk:22 ajansen@tails.cs.ucl.ac.uk
```

3. Open a new terminal and log onto UoB bluecp3

4. Run rsync in screen
```
screen

rsync /projects/Butterfly_genome_analysis/alex/C3_Aricia_agestis/00_raw_reads_museum/*gz -auve \
"ssh -p 3000" ajansen@localhost:/SAN/ugi/LepGenomics/C3_Aricia_agestis/00_raw_reads_museum/

#cntrl-A-D to log out of screen
```

Useful screen options
```
#lists active screens
screen -ls

#reconnect to screen
screen -r screenname

#kill a screen
screen -X -S screenname kill

```
