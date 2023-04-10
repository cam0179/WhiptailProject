# DESeq2 Analysis of the Checkered Whiptail Lizard, _Aspidoselis tesselata_

The code in the file **WhiptailCorePipeline.R** is the pipeline for differential expression analysis using RNA-seq data collected from three different populations of _A. tesselata_. This pipeline was created in collaboration with Dr. Randy Klabacka by Akila Abesinghe, Ally Swank, Camille Miceli, Gabriel Amorim de Albuquerque Silva, and Melissa Gathman. 

## Input files:
1. Gene count matrix that would be a comma delimited output from the bionformatic processing of raw reads. Our gene count matrix was created using StringTie. 
    - **gene_count_matrix.csv**
2. Phenotype tab delimited data file that contains each sample ID in column 1 to match each sample ID in row 1 of the gene count matrix. This file will contain data for variables of interest. In this study, we are evaluating differential gene expression based on the location of collection and the type of capture method used. 
    - **PHENO_DATA.txt**

This pipeline could also be used to analyze differences in transcript read counts using the file **transcript_count_matrix.csv**.
