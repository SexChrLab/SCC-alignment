`To do: Add threads and read orientation to config, softcode threads in snakefile(s).`

# Overall Description

These Snakemake workflows are for sex chromosome complement informed gene quantification for short read, paired end RNA sequencing data.  For each sample, the sex is used to determine which sex chromosome complement version of the reference genome is used for alignment/pseudoalignment.  For females (XX, no Y chromosome),  we use the reference genome with chromosome Y hard masked.  For males (XY, with a Y chromosome), we use the reference genome with the pseudoautosomal regions (PARs) hard masked on chromosome Y.  There are separate Snakemake pipelines for doing a full alignment (`hisat2`) with gene quantification (`featureCounts`) or instead doing pseudoalignment and quantification using SCC versions of the transcriptome (`salmon`) which is much faster.  

This page summarizes the steps of the workflows and how to use the output, but a detailed tutorial on how to perform the workflows with test data from Genome in a Bottle is given in the `tutorials` directory.

# Gene quantification pipeline with full alignment and gene counts

The Snakemake workflow `rnaseq_data_processing_hisat2.snakefile` does an alignment to sex chromsoome complement refernece genome using `hisat2` and quantifies gene expression using the genome annotation using `featureCounts`.  

## Steps for full alignment gene quantification 
1. `hisat2` is used to align the samples to the appropriate SCC reference genome (HISAT2 indexed files are available in the /references directory in the Docker)
2. `bamtools` is used to index the alignment output files
3. `featureCounts` is used to count reads that align to exon regions
4. `featureCounts` results output with both gene IDs and transcript IDs for use in further analysis

# Pseudoalignment and transcript quantification

The Snakemake workflow `rna_data_processing_salmon.snakefile` uses the `salmon` algorithm to quantify expression.  

## Steps for the pseudoalignment gene quantification 
1. `salmon quant` is used to map reads to SCC reference transcriptome (SCC reference transcriptome are provided in the /references directory in the Docker)
2. Quantified expression will be in the `quantified_rna_salmon` directory, with specific subdirectories for each sample
3. `quant.sf` file inside each result subdirectory will give the quantified expression for the transcripts in each sample

# Running the SCC-aware gene quantification workflows

Steps for running the pipeline: 

1) Create a config JSON for your RNA samples (see `custom_config`) 

2) Activate the conda environment (see main `SCC-alignment` page) and Docker environment if using

3) Open the `.snakefile` file in a text editor corresponding to which alignment procedure you would like to do and make sure that the name of your config JSON is set correctly in the `configfile` variable

4) To do full alignment with `hisat2` and quantify genes with `featureCounts`

a) Test Snakemake pipeline: `snakemake -np -s rnaseq_data_processing_hisat2_piped.snakefile`

b) Run Snakemake pipeline: `snakemake -s rnaseq_data_processing_hisat2_piped.snakefile`

c) Quantification results for each sample given in `feature_counts_rna_hisat` directory

5) To do pseudoalignment and gene quantification with Salmon:

a) Test Snakemake pipeline: `snakemake -np -s rnaseq_data_processing_salmon.snakefile`

b) Run Snakemake pipeline: `snakemake -s rnaseq_data_processing_salmon.snakefile`

c) Quantification results for each sample given in `quantified_rna_salmon` directory

Example of how to run on a high performance cluster using a conda environment, local copies of the SCC reference genomes, and slurm workflow manager are given in `.sbatch` files.

# Test data

For testing and development, you can use RNA sequencing data from the Genome in a Bottle project from NIST. These samples are from families but having related individuals is not at all required for our pipeline. We list these as test data because they are publicly available, high quality human samples with a known sex chromosome complement specifically intended for methods development.

Male sample: 

https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg004_gm24143.mrna.R1.fastq.gz

https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg004_gm24143.mrna.R2.fastq.gz

Female sample: 

https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg002_gm24385.mrna.R1.fastq.gz

https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg002_gm24385.mrna.R2.fastq.gz

# Manuals

Hisat2: https://daehwankimlab.github.io/hisat2/manual/

Salmon: https://salmon.readthedocs.io/en/latest/

# Citations

`Hisat2` RNA sequencing alignment algorithm: 

Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 2019 Aug;37(8):907-915. doi: 10.1038/s41587-019-0201-4. Epub 2019 Aug 2. PMID: 31375807; PMCID: PMC7605509.

`featureCounts` gene quantification algorithm from the `subread` package:

Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014 Apr 1;30(7):923-30. doi: 10.1093/bioinformatics/btt656. Epub 2013 Nov 13. PMID: 24227677.

`Salmon` RNA sequencing pseudoaligner and gene quantification algorithm:

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017 Apr;14(4):417-419. doi: 10.1038/nmeth.4197. Epub 2017 Mar 6. PMID: 28263959; PMCID: PMC5600148.
