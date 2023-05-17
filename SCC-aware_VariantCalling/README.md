# Overall Description

These Snakemake workflows are for sex chromosome complement informed (SCC-aware) whole-genome short-read sequence alignment and variant calling. For each sample, the genotype, identified using SCC-check pipeline (or reported sex), is used to determine which sex chromosome complement version of the reference genome is used for alignment/pseudoalignment. For females (no Y chromosome), we use the reference genome with Y chromosome hard-masked. For males (samples with a Y chromosome), we use the reference genome with the pseudoautosomal regions (PARs) hard masked on chromosome Y. A JSON file is used to hold all the file information in the previous directory (same config for SCC_check).

# System requirements

Minimum: 4 available threads; Recommended: >12 threads.

Minimum: ~16Gb RAM; Recommended: >24Gb RAM.

# Test data

For testing and development, we used whole genome sequencing data from the Precision V2 project from Genome in a Bottle project from NIST.
https://data.nist.gov/od/id/mds2-2336

Male (with Y chromosome) samples: 

HG002 (Ashkenazim son):
https://data.nist.gov/od/ds/ark:/88434/mds2-2336/input_fastqs/HG002.novaseq.pcr-free.35x.R1.fastq.gz

https://data.nist.gov/od/ds/ark:/88434/mds2-2336/input_fastqs/HG002.novaseq.pcr-free.35x.R2.fastq.gz

HG003 (Ashkenazim father): 
https://data.nist.gov/od/ds/ark:/88434/mds2-2336/input_fastqs/HG003.novaseq.pcr-free.35x.R1.fastq.gz

https://data.nist.gov/od/ds/ark:/88434/mds2-2336/input_fastqs/HG003.novaseq.pcr-free.35x.R2.fastq.gz

Female (without Y chromosome) samples: 

HG004 (Ashkenazim mother):
https://data.nist.gov/od/ds/ark:/88434/mds2-2336/input_fastqs/HG004.novaseq.pcr-free.35x.R1.fastq.gz

https://data.nist.gov/od/ds/ark:/88434/mds2-2336/input_fastqs/HG004.novaseq.pcr-free.35x.R2.fastq.gz

<TO-DO: need to add one more female>


# Running the workflow

Steps for running the pipeline: 
1) Created a config JSON for your DNA samples (see `custom_config`) 
2) Activate the conda environment (see main `SCC-alignment` page)
3) Open `.snakefile` files in a text editor and make sure that the name of your config JSON is set correctly
4) Samples with Y chromosome: (in `Y_samples` entry in DNA config JSON)
a) Test Snakemake pipeline: `snakemake -np -s Snakefile_XY_SNPs.snakefile`
b) Run Snakemake pipeline: `snakemake -s Snakefile_XY_SNPs.snakefile`

5) Samples without Y chromosome: (in `X_samples` entry in DNA config JSON)
a) Test Snakemake pipeline: `snakemake -np -s Snakefile_XX_SNPs.snakefile`
b) Run Snakemake pipeline: `snakemake -s Snakefile_XX_SNPs.snakefile`

Example of how to run on a high performance cluster using slurm workflow manager: 
```
#!/bin/bash
#SBATCH --job-name=testSC  # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=splaisie@asu.edu # send-to address
#SBATCH -t 4-00:00
#SBATCH -p serial
#SBATCH -n 2

source activate SCCalign_v3
snakemake -s Snakefile_XX_SNPs.snakefile -j 100 --rerun-incomplete --latency-wait=60 --cluster "sbatch -n 2 -p serial --mem=50G --mail-type=FAIL --mail-user=splaisie@asu.edu"
snakemake -s Snakefile_XY_SNPs.snakefile -j 100 --rerun-incomplete --latency-wait=60 --cluster "sbatch -n 2 -p serial --mem=50G --mail-type=FAIL --mail-user=splaisie@asu.edu"

```

6) Output vcfs will be populated in subdirectories
















# Citations 
If you use this code in your work, please cite its source and tools:

Plasier et al. 202x. SCC-aware analysis tutorial.

bwa: McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110

or 
 
Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191,

and 

SAMtools: Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352,

GATK: McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research, 20(9), 1297-1303.  http://www.genome.org/cgi/doi/10.1101/gr.107524.110.



