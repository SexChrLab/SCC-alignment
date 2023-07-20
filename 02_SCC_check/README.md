# Overall Description

This Snakemake workflow is for identifying the sex chromosome complement of a set of samples using short-read sequencing data using comparative read depth. This information should be then be incorporated into the JSON config file (`DNA_samples.config.json`) for use in SCC-aware variant calling and/or SCC-aware gene quantification (`RNAseq_samples.config.json`).

# Techniques for estimating the sex chromosome complement of each sample

It is important to know if a sample possesses a Y chromosome or not to understand potential biases this can introduce into traditional read mapping-based genomic analyses. Here we present various methods you can use to gather evidence that you can to know the sex chromosome complement of your samples or confirm reported sex of the sample.  

We acknowledge that previous iterations of this read mapping depth approach were computationally intensive and were (somewhat justifiably) avoided. This pipeline attempts to bypass these limitations by subsampling the DNA sequencing data to 1X coverage prior to alignment to significantly reduce runtime. However, this subsampling restricts the inferential power of the analysis to simply ask, "Are Y chromosome reads present in the sequence data?". If additional information is needed to be inferred from the data (e.g. investigating X chromosome copy number (CN)), you should not subsample the data.

## Samples with DNA sequencing
In this module, we provide code that uses DNA sequencing data mapped to a Y PARs-masked reference genome to identify evidence of a Y chromosome in the sequence reads. This information can be used to validate reported sex of the sample metadata or as an independent investigation into the individuals' genotype.

If you have DNA sequencing data, we have provided code that uses the `indexcov` application to calculate relative read depth for chromosome Y in order to determine whether a Y chromosome is present in your sample.  We subsample the read files to ~1X coverage (assuming PE 150bp short reads) and quickly map these read to a Y PAR-masked reference.  If the sample possesses a Y chromosome the chrY copy number (CN) will be greater than ~0. This pipeline is built for speed to encourage its use, but is not accurate for determining chrX copy number. If chrX copy number is important for your use case, simply align all the reads (20-30X) and run indexcov following the programs documentation (https://github.com/brentp/goleft/tree/master/indexcov).  The subsequencing step is in the `subsample_reads` rule in `SCC-check.snakefile`, so to align all the reads simply remove this rule and set the input of the `map_reads` rule to be the original reads in the `reads/` directory as originally the input for `subsample_reads`.    

### Sex chromosome complement (SCC) check pipeline for DNA sequencing data

The rules in the SCC check pipeline for DNA sequencing data `SCC-check.snakefile` are as follows, for each sample: 
1) Subsample the reads using `seqtk` to create a representative set of reads that is faster to run
2) Map reads to the Y-PARs masked reference genome using `minimap2`
3) Index the resulting alignment (bam) files using `samtools index`
4) Calculate the relative read depth of chromosome Y using `goleft indexcov`

Once the pipeline has been run, your results will be in the directory you specified in your DNA custom config, in the "indexcov_dir" variable.  These results will include html files where you can view your samples according to their read depth on chrY as well as a ped file called `indexcov-indexcov.ped` which contains this information in table form.  We have provided a simple Python script to ask if the calculated chrY copy number is greater than 0.25, and if so, predict that this sample has a Y chromosome.  These predictions are output in a simple table that can be used to fill in the `Y_samples` (have Y chromosomes) and `X_samples` (do not have Y chromosomes) lists in the custom configs before doing SCC-aware variant calling or SCC-aware gene quantification (for samples where you have matching gene expression).  

Steps to run the SCC check pipeline for DNA sequencing data: 
1) Create a config JSON for your DNA samples (see `custom_config`)  
3) Activate the conda environment (see main `SCC-alignment` page)
4) Open `SCC-check.snakefile` in a text editor and make sure that the name and path of your config JSON is set correctly in the `configfile` variable
5) Test Snakemake pipeline: `snakemake -np -s SCC-check.snakefile`
6) Run Snakemake pipeline: `snakemake -s SCC-check.snakefile`

Example of how to run as a job submission on a high performance cluster using slurm workflow manager with SCCalign_v3 conda environment installed: 
```
#!/bin/bash
#SBATCH --job-name=sexcheck  # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=splaisie@asu.edu # send-to address
#SBATCH -t 4-00:00
#SBATCH -p serial
#SBATCH -n 2

source activate SCCalign_v3
snakemake -s SCC-check.snakefile -j 100 --rerun-incomplete --latency-wait=60 --cluster "sbatch -n 2 -p serial --mem=50G --mail-type=FAIL --mail-user=splaisie@asu.edu"
```

7) Copy Python script `inferred_SCC.py` into the directory containing `indexcov` results (including file named `indexcov-indexcov.ped`) and run to get a simple output file that summarizes which samples have a Y chromosome
``` 
python inferred_SCC.py
```


## Samples with RNA sequencing 

`NEEDS TO BE COMPLETED: add scripts that predict sex counting the alignments these genes using the Y PARs masked reference genome`
If you only have RNA sequencing data, we have used high expression of the XIST gene to indicate presence of at least 2 X chromosomes in the sample and high expression of Y chromosome genes whose expression is not limited to the testis to indicate the presence of a Y chromosome in the sample.  The genes we use include: DDX3Y, KDM5D, USP9Y, UTY, ZFY.   

# Citations 

If you use this script in your work, please cite the relevant softwares that make it work:

seqkit: https://github.com/lh3/seqtk

Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

SAMtools: Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

Indexcov: Pedersen, B. S., Collins, R. L., Talkowski, M. E., & Quinlan, A. R. (2017). Indexcov: fast coverage quality control for whole-genome sequencing. GigaScience, 6(11). https://doi.org/10.1093/gigascience/gix090
