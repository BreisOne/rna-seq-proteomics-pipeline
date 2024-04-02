#!/bin/bash
#SBATCH -t 03:30:00 # execution time. Ex: 1 hour
#SBATCH -c 12 # number of cores. Ex:8
#SBATCH --mem 80GB 
#SBATCH --mail-user b.mascat@gmail.com
#SBATCH --mail-type ALL
# set up the software environment

module load htseq

input=$LUSTRE/rnaseq/results/STAR/2020-03-30_aligment
#input=$LUSTRE/rnaseq/results/STAR/2019-12-20_aligment/Sin_reads_no_mapeados

gtf=$LUSTRE/annotation_genome_hs/Homo_sapiens.GRCh38.98.gtf
output_path=$LUSTRE/rnaseq/results/STAR/2020-03-30_aligment
#output_path=$LUSTRE/rnaseq/results/STAR/2019-12-20_aligment/Sin_reads_no_mapeados

echo "sample loaded $1"

htseq-count -m union -f bam -t exon -i gene_id  $input/${1}_Aligned.sortedByCoord.out.bam  $gtf > $output_path/${1}_countMatrix.txt




