#!/bin/bash
#SBATCH --job-name=hisat # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL               # notifications 
#SBATCH --mail-user=splaisie@asu.edu # send-to address
#SBATCH -n 2
#SBATCH -p serial

cd /scratch/splaisie/XYAlign/RNA/
source activate SCCalign_v3 
snakemake -s rnaseq_data_processing_hisat2.snakefile --rerun-incomplete --latency-wait=60 -j 100 --cluster "sbatch -n 2 -p serial --mem=32G --mail-type=FAIL --mail-user=splaisie@asu.edu"
