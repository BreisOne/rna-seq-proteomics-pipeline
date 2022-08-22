#!/bin/bash
#SBATCH -t 0:30:00 # execution time. Ex:00:30:00
#SBATCH --mail-user b.mascat@gmail.com
#SBATCH --mail-type ALL
#SBATCH -n 8 # number of tasks. Ex:8
#SBATCH -N 1 # always specify this parameter

#Module load
module load gcc/6.4.0 R/3.6.0

#RScrypt
Rscript  4.MatrixCounts_merge.R
