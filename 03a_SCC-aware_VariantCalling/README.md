# Overall Description

These Snakemake workflows are for sex chromosome complement informed (SCC-aware) whole-genome short-read sequence alignment and variant calling. For each sample, the genotype, identified using SCC-check pipeline (or reported sex), is used to determine which sex chromosome complement version of the reference genome is used for alignment/pseudoalignment. For females (no Y chromosome), we use the reference genome with Y chromosome hard-masked. For males (samples with a Y chromosome), we use the reference genome with the pseudoautosomal regions (PARs) hard masked on chromosome Y. A JSON file is used to hold all the file information in the previous directory (same config for SCC_check).

# Steps of the sex chromosome complement (SCC)-aware pipeline for variant calling
The rules in the SCC-aware variant calling pipeline SNPs_XX.snakefile and SNPs_XY.snakefile are as follows:
1. Samples with a Y chromosome must be specified in the `Y_samples` variable in the custom config passed in through the `configfile` parameter
2. Samples without a Y chromosome must be specified in the `X_samples` variable in the custom config passed in through the `configfile` parameter
3. SNPs_XY.snakefile will call variants on the `Y_samples` samples
4. SNPs_XX.snakefile will call variants on the `X_samples` samples
5. Each sample's reads will be aligned to the appropriate SCC reference genome
6. Read groups are added to the alignment files
7. `gatk HaplotypeCaller` will be used to call single nucleotide variants across all samples for all chromosomes
8. `gatk CombineGVCFs` will be used to combine the called variant files to show all variants across all samples in `Y_samples` or `X_samples`
9. `gatk GenotypeGVCFs` will be used to joint genotype all variants across all samples
10. `gatk SelectVariants` will be used to filter for biallelic single nucleotide polymorphisms (SNPs)

# Running the SCC-aware variant calling pipeline

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


# System recommendations

Minimum: 4 available threads; Recommended: >12 threads.

Minimum: ~16Gb RAM; Recommended: >24Gb RAM.

# Test data

## Whole Genome Sequencing data
For testing and development, you can use 40X Novaseq whole genome sequencing data from the Genome in a Bottle project from NIST. These samples are from families but having related individuals is not at all required for our pipeline.  We list these as test data because they are publicly available, high quality human samples with a known sex chromosome complement specifically intended for methods development.  

https://console.cloud.google.com/storage/browser/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x

### Male (with Y chromosome) samples: 


HG003 (Ashkenazim father):

https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG003.novaseq.pcr-free.40x.R1.fastq.gz?_ga=2.237371716.-57486372.1684265308

https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG003.novaseq.pcr-free.40x.R2.fastq.gz?_ga=2.237371716.-57486372.1684265308


HG006 (Chinese father):

https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG005.novaseq.pcr-free.40x.R1.fastq.gz?_ga=2.203778004.-57486372.1684265308

https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG005.novaseq.pcr-free.40x.R2.fastq.gz?_ga=2.266211275.-57486372.1684265308

### Female (without Y chromosome) samples: 


HG004 (Ashkenazim mother):

https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG004.novaseq.pcr-free.40x.R1.fastq.gz?_ga=2.266211275.-57486372.1684265308

https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG004.novaseq.pcr-free.40x.R2.fastq.gz?_ga=2.266211275.-57486372.1684265308


HG007 (Chinese mother):

https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG007.novaseq.pcr-free.40x.R1.fastq.gz?_ga=2.266211275.-57486372.1684265308

https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG007.novaseq.pcr-free.40x.R2.fastq.gz?_ga=2.229190616.-57486372.1684265308


# Citations 
If you use this code in your work, please cite its source and tools:

Plasier et al. 202x. SCC-aware analysis tutorial.
 
Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191,

and 

SAMtools: Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352,

GATK: McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research, 20(9), 1297-1303.  http://www.genome.org/cgi/doi/10.1101/gr.107524.110.



