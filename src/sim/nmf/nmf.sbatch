#!/bin/bash
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling
#################
#set a job name
#SBATCH --job-name=nmf-sim
#################
#a file for job output, you can check job progress, append the job ID with %j to make it unique
#SBATCH --output=nmf-sim.%j.out
#################
# a file for errors from the job
#SBATCH --error=nmf-sim.%j.err
#################
#time you think you need; default is 2 hours
#format could be dd-hh:mm:ss, hh:mm:ss, mm:ss, or mm
#SBATCH --time=60:00
#################
#Quality of Service (QOS); think of it as job priority, there is also --qos=long for with a max job length of 7 days, qos normal is 48 hours.
# REMOVE "normal" and set to "long" if you want your job to run longer than 48 hours,
# NOTE- in the hns partition the default max run time is 7 days , so you wont need to include qos
#SBATCH --qos=normal
# We are submitting to the dev partition, there are several on sherlock: normal, gpu, owners, hns, bigmem (jobs requiring >64Gigs RAM)
# The more partitions you can submit to the less time you will wait, you can submit to multiple partitions with -p at once in comma separated format.
#SBATCH -p normal,hns
#################
#number of nodes you are requesting, the more you ask for the longer you wait
#SBATCH --nodes=1
#################
# --mem is memory per node; default is 4000 MB per CPU, remember to ask for enough mem to match your CPU request, since
# sherlock automatically allocates 8 Gigs of RAM/CPU, if you ask for 8 CPUs you will need 32 Gigs of RAM, so either
# leave --mem commented out or request >= to the RAM needed for your CPU request.
#SBATCH --mem=8000
#################
# Also, if you submit hundreds of jobs at once you will get hundreds of emails.
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
# Remember to change this to your email
#SBATCH --mail-user=kriss1@stanford.edu

ml git
ml gcc/4.9.1
ml R/3.3.0
ml python/3.3.2
ml unload intel/2015
ml load R/3.3.0

source /scratch/users/kriss1/writing/research/unsupervised_prob/.env
Rscript /scratch/users/kriss1/writing/research/unsupervised_prob/src/sim/nmf/nmf_expers.R
R -q -e 'Sys.sleep(3400)'