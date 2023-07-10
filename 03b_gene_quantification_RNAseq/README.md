To do: Add threads and read orientation to config, softcode threads in snakefile(s).

# Overall Description

These Snakemake workflows are for sex chromosome complement informed gene quantification for short read, paired end RNA sequencing data.  For each sample, the sex is used to determine which sex chromosome complement version of the reference genome is used for alignment/pseudoalignment.  For females (XX, no Y chromosome),  we use the reference genome with chromosome Y hard masked.  For males (XY, with a Y chromosome), we use the reference genome with the pseudoautosomal regions (PARs) hard masked on chromosome Y.  There are separate Snakemake pipelines for doing a alignment (`hisat2`) with gene quantification (`featureCounts`) or doing pseudoalignment and quantification (`salmon`) which is much faster.  A JSON file is used to hold all the file information; a script that can be used to generate the sample specific information is provided (see `custom_config`).  

# Alignment and gene quantification

The Snakemake workflow `rnaseq_data_processing_hisat2.snakefile` does an alignment to sex chromsoome complement refernece genome using `hisat2` and quantifies gene expression using the genome annotation using `featureCounts`.  Interrim files produced during alignment and file processing are marked as temporary, so they will be generated and deleted when the step that needs them is completed so as not to fill up hard drive space unnecessarily.  The rules for running `featureCounts` are set up to output results quantifying gene expression on the exon regions and output with both the gene IDs and the transcript IDs.

# Pseudoalignment and transcript quantification

The Snakemake workflow `rna_data_processing_salmon.snakefile` uses the `salmon` algorithm to quantify expression.  Salmon indices are provided in `references/` for reference transcriptome with sequences from chromosome Y masked, so these are used according to the sex of the sample (as specified in the `X_samples` and `Y_samples` entries in the configuration JSON).  Quantified expression will be in the `quantified_rna_salmon` directory, with specific subdirectories for each sample.  The `quant.sf` file inside each result subdirectory will give the quantified expression for the transcripts in each sample.

# Test data

In the development of this tutorial we used RNA sequencing files collected by the Genome in a Bottle project with NIST. You can download these files for testing and development purposes. 

Male sample: 
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg004_gm24143.mrna.R1.fastq.gz
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg004_gm24143.mrna.R2.fastq.gz

Female sample: 
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg002_gm24385.mrna.R1.fastq.gz
https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg002_gm24385.mrna.R2.fastq.gz

# Running the workflow

Steps for running the pipeline: 

1) Created a config JSON for your RNA samples (see `custom_config`) 

2) Activate the conda environment (see main `SCC-alignment` page)

3) Open `.snakefile` file in a text editor corresponding to which alignment procedure you would like to do and make sure that the name of your config JSON is set correctly

4) To do full alignment with `hisat2` and quantify genes with `featureCounts`

a) Test Snakemake pipeline: `snakemake -np -s rnaseq_data_processing_hisat2_piped.snakefile`

b) Run Snakemake pipeline: `snakemake -s rnaseq_data_processing_hisat2_piped.snakefile`

c) Example of how to run on a high performance cluster using slurm workflow manager given in `run_preprocessing_hisat.sbatch`

d) Quantification results for each sample given in `feature_counts_rna` directory

5) To do pseudoalignment and gene quantification with Salmon:

a) Test Snakemake pipeline: `snakemake -np -s rnaseq_data_processing_salmon.snakefile`

b) Run Snakemake pipeline: `snakemake -s rnaseq_data_processing_salmon.snakefile`

c) Example of how to run on a high performance cluster using slurm workflow manager given in `run_preprocessing_salmon.sbatch`

d) Quantification results for each sample given in `quantified_rna` directory

# Citations

`Hisat2` RNA sequencing alignment algorithm: 

Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol. 2019 Aug;37(8):907-915. doi: 10.1038/s41587-019-0201-4. Epub 2019 Aug 2. PMID: 31375807; PMCID: PMC7605509.

`featureCounts` gene quantification algorithm from the `subread` package:

Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014 Apr 1;30(7):923-30. doi: 10.1093/bioinformatics/btt656. Epub 2013 Nov 13. PMID: 24227677.

`Salmon` RNA sequencing pseudoaligner and gene quantification algorithm:

Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. 2017 Apr;14(4):417-419. doi: 10.1038/nmeth.4197. Epub 2017 Mar 6. PMID: 28263959; PMCID: PMC5600148.
