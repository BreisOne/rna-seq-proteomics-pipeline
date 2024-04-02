#!/bin/bash
#SBATCH -t  5:00:00 # execution time. Ex: 1 hour
#SBATCH -c 12 # number of cores. Ex:8
#SBATCH --mem 80GB 
#SBATCH --mail-user b.mascat@gmail.com
#SBATCH --mail-type ALL
# set up the software environment

#module load fastqc
module load gcc star
#module load samtools

#Create Genome Index

STAR --runMode genomeGenerate --genomeDir $LUSTRE/Reference_Genome_HS --genomeFastaFiles $LUSTRE/Reference_Genome_HS/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
