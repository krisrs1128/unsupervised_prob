#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=lda
#################
#a file for job output, you can check job progress
#SBATCH --output=lda.out
#################
# a file for errors from the job
#SBATCH --error=lda.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=00:20:00
#################
#memory per node; default is 4000 MB per CPU
#SBATCH --mem=4000
#SBATCH --partition=normal,hns
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kriss1@stanford.edu
ml R/3.3.0
Rscript /scratch/users/kriss1/writing/research/unsupervised_prob/src/antibiotics-study/lda.R
