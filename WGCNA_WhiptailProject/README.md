# Weighted Gene Co-expression Network Analysis (WGCNA) of the Checkered Whiptail Lizard, _Aspidoselis tesselata_

The code written in **WGNCA.R** was developed by Ally Swank to assess co-expression of genes across three different populations of _A. tesselata_ using two different capture methods. This work was done in collaboration with Dr. Randy Klabacka by Akila Abesinghe, Ally Swank, Camille Miceli, Gabriel Amorim de Albuquerque Silva, and Melissa Gathman.

## Directories, input, and output files:
1. **raw-data**: directory containing all of the input files to be run in the WGCNA.R script. 
    - **gene_count_matrix.csv**: a comma delimited output from the bionformatic processing of raw reads. Our gene count matrix was created using StringTie.
    - **PHENO_DATA.txt**: tab delimited data file that contains each sample ID in column 1 to match each sample ID in row 1 of the gene count matrix. This file will contain data for variables of interest. In this study, we are evaluating differential gene expression based on the location of collection and the type of capture method used.  
2. **results**: directory containing any modules of interest from the WGCNA that you wish to see a specific gene list for.
3. **outputs**: directory containing any output figures 

This pipeline could also be used to analyze co-expression in transcript read counts using the file **transcript_count_matrix.csv** in the **raw-data** directory.
