#!/usr/bin/env bash

#SBATCH --job-name=HISAT2_idx
#SBATCH --output=logs/HISAT2-idx_%A.log
#SBATCH --cpus-per-task=20
#SBATCH --mem=100gb
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --qos=large
#SBATCH --constraint=dmc
#SBATCH --mail-user=gas0042@auburn.edu
#SBATCH --mail-type=ALL


source activate /home/aubclsb0305/.conda/envs/Hisat2

workDir=/scratch/WhiptailProject

baseName=${workDir}/Reference/P_raffonei/rPodRaf1
outIndex=${workDir}/Reference/P_raffonei/Hisat2_idx/rPodRaf1

mkdir -p ${workDir}/Reference/P_raffonei/Hisat2_idx/

extract_splice_sites.py ${baseName}.gtf > ${baseName}.ss
extract_exons.py ${baseName}.gtf > ${baseName}.exon

hisat2-build -p 20 --ss ${baseName}.ss --exon ${baseName}.exon -f ${baseName}.fna ${outIndex}