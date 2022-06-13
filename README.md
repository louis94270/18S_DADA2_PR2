# 18S DADA2 PR2

A metabarcoding pipeline using DADA2 and the PR2 reference database for the 18S marker 

## Required programs

1. R (https://www.r-project.org/)
2. DADA2 Bioconductor package (https://benjjneb.github.io/dada2/dada-installation.html)
3. Cutadapt (https://cutadapt.readthedocs.io/en/stable/installation.html)

## How to use the pipeline 

1. Download the R code on your local computer (or clone this directory)
2. Put all your demultiplexed sequences in a folder
3. Change file names to have the following pattern if not the case: Forward reads files end with "_R1_001.fastq", and Reverse reads files end with "_R2_001.fastq"
4. Change the different paths in the R code to match your computer set-up (line 11 and line 45)
5. Change the primer for the one you are using in the R code (line 30 and line 32)
6. You can start running the code ! 
7. Be carefull of looking at your quality profils before running line 100 to adapt your arguments, you can use the amplicon sequence tool I developped to look at your sequences and calculate the overlap to make sure you are not cutting to many base pairs 
8. Be carefull along the code some objects dimensions will need updating

The different steps are very well explained in these tutorials : 
1. https://astrobiomike.github.io/amplicon/dada2_workflow_ex
2. https://benjjneb.github.io/dada2/tutorial.html