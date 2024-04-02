#!/bin/bash
#SBATCH -t  00:50:00 # execution time. Ex: 1 hour
#SBATCH -c 8 # number of cores. Ex:8
#SBATCH --mem 80GB
#SBATCH --mail-user b.mascat@gmail.com
#SBATCH --mail-type ALL

#Module load
module load gcc star

#Define Base variable

base=`basename $1 .fastq`

# Run STAR

STAR --runMode alignReads --genomeDir $LUSTRE/Reference_Genome_HS/genome-STAR-index --readFilesIn $LUSTRE/Run_2020-01-15/$1 --outFileNamePrefix $LUSTRE/rnaseq/results/STAR/2020-03-30_aligment/${base}_URT_ --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 10 --outFilterMismatchNmax 10 --alignIntronMin 21 --limitOutSJcollapsed 2000000 --outSAMunmapped Within

