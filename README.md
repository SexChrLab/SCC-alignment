# Sex chromosome complement (SCC) informed alignment pipelines

This repository of pipelines and code from Plaisier et al. 202x shows how to effectively account for sex chromosomes in common human genomic analyses: read alignment, variant calling, and gene expression quantification.

The description below is geared towards running these pipelines in a Linux environment and shows how to analyze short-read sequencing data from human samples with the CHM13v2.0 (telomere-to-telomere) or GRCh38 reference genomes masked to account the sex chromosome complement in the samples being sequenced and analyzed.  For those studying non-human species with a known sex determination system that uses sex chromosomes, we hope that the detailed explanations here and in the manuscript will help you to adapt these techniques using the published reference genome sequences; please open an issue or email the Sex Chromosomes Lab at Arizona State University if you need help (http://www.sexchrlab.org/).

# Reasons to make the switch to sex-chromosome complement aware genomics methods

Due to the unique biology of the sex chromosomes, accurate measurement of genetic variation and gene expression from the sex chromosomes can be difficult.  As a result, sex chromosomes are often excluded from genomic analyses, even though they are an essential part of the genetic variation within individuals and can have a significant impact on cell biology as it pertains to health and disease.  Immune disorders and cancer are known to have a sex bias in incidence and outcome. Sex linked differences in gene expression have been observed across tissues and consistently over the course of life. Omitting or incorrectly handling sex chromosomes in genomic analysis can potentially lead to false positives or false negatives from the sex chromosome genes in genomic studies.  While learning and accounting for sex chromosomes in your samples can require a few extra processing steps, the benefit is that you will get more accurate measurements of the variants and genes in your samples and better understand what is underlying the cell biology in the samples you are studying.

The pipelines and resources presented here aim to facilitate correct and accurate alignment of short read sequencing data to human reference genome adjusted to account for the sex chromosomes present in the sample.  For samples that do not contain a Y chromosome, the Y chromosome sequence will be masked so that no reads are incorrectly assigned to it.  For samples that contain a Y chromosome, the regions of the Y chromosome that have perfectly matching regions in the X chromosome (psuedoautosomal regions (PARs)) are masked on the Y chromosome so that reads in this ambiguous regions are not assigned randomly and variant calling on the sex chromosomes is adjusted to match the ploidy. 

# Knowing the sex chromosome complement in your samples
In order to do sex chromosome complement aware analysis, it is necessary to know the sex chromosome complement (SCC) present in each of the samples you are working with.  

When available, the reported sex from the individual from whom the samples being analyzed can be used to infer the sex chromosome complement.  For example, if you are studying tissue samples from a patient that is reported to be male, most males have Y chromosomes in the vast majority of cells in their body, so you could go forward and process those samples assuming they have a Y chromosome.  Similarly, for samples from individuals that are reported as female, you could go forward processing those samples as having no Y chromosome.  This approach makes assumption that do not always hold true; sex chromosomes have been shown to be lost in samples from individuals of old age and gains and losses have been reported in tumor samples, for example.   

To be more certain of the sex chromosome complement of your samples, you can use the provided sample sex chromosome complement checking module (`02_SCC_check`) to confirm presence or lack of a Y chromosome.  For DNA sequencing data, relative read depth on chromosome Y can be used to confirm the presence of a Y chromosome.  If you do not have DNA sequencing data and/or are instead working with RNA sequencing data and do not have or trust the reported sex,  `COMPLETE!!!!` 

# Overview of the SCC genomics analysis modules

1. Custom configuration files (`01_custom_config`)

In this module, we have provided scripts to create custom configuration files in JSON format.  You will need specific configs for DNA analysis and RNA analysis because different parameters are requested for specific SCC-informed analyses.

2. Sex chromosome complement check (`02_SCC_check`)

In this module, we provide code and guidelines to determine whether samples have a Y chromosome. For DNA samples, code is provided to use relative read depth on the Y chromosome. For RNA samples, we have `COMPLETE!!!`

3a. Sex chromosome complement-aware variant calling (`03a_SCC-aware_Variant Calling`)

In this module, we show how to use DNA sequencing data to determine variants across the genome. We use the inferred Sex Chromosome Complement (SCC) from module 2, or sex reported in sample metadata, to assign the appropriate SCC aware reference genome and downstream gentyping criteria for each sample. In short, we perform sex chromosome complement specific sequence alignment and determine genomic variants considering the biologically relevant ploidy levels across the genome using GATK.

3b. SCC-aware gene expression analyses (`03b_SCC_gene_quantification`) 

In this module, we use HISAT2 aligner with a gene quantification algorithm (`featureCounts`) or Salmon pseudoaligner to quantify the expression of genes in samples assayed with RNA sequencing.  

# Sex chromosome complement reference genomes

The basis of the sex chromosome complement aware genomics analysis presented in this repository is the use of versions of the human reference genome sequence adjusted to allow more accurate alignment of reads based on the sex chromosomes present in the sample.  Manuscripts describing sex chromosome complement alignment in detail are given in the Citations section below, but in summary we have created two versions of the human reference genome sequence that is used to align sequencing data to. 

1.  Chromosome Y masked reference genome sequence for samples that do not have a Y chromosome (such as XX female samples)
2.  Chromosome Y PARs masked reference for samples that do have a Y chromosome (such as XY male samples)

We have provided sex chromosome complement versions of the GRCh38 and CHM13v2 releases, but detailed instructions on how these reference genome sequences were made are provided in the `references` directory.  

A summary of how these were generated is as follows: 
1. Download the human genome reference sequence and transcriptome sequence
2. Hard mask the entire sequence of chromosome Y (switch all nucleotides to N) to align to samples without a Y chromosome
3. Determine the boundaries of the pseudoautosomal regions (PARs), regions of complete sequence identity (either annotated with the released genome sequence or determined by aligning the X and Y chromosome sequence to each other)
4. Hard mask the original reference genome and transcriptome sequence in the PAR regions on chromosome Y
5. Use index functions of alignment algorithms where needed with the Y-masked and Y PARs masked versions of the reference genome (such as `hisat2 build`)
6. Detailed instructions are given in the `references` directory

# Order of operations

This section describes the order in which our modules can be run to give you an idea of how to proceed given the type of data you have.

Detailed tutorials are provided in the `tutorials` directory.

In summary, first, you will need to set up a custom config that describes the necessary information about your experiment and data files.  Code and examples for this are given in `01_custom_config`  There will be specific fields required for the SCC check, SCC-aware variant calling, and gene quantification modules; the code provided will help you to generate a template for that json which you can manually fill in with the paths to your data and specific experimental details.  Most importantly, you will need to which samples have a Y chromosome and which do not. These lists are used to determine which sex chromosome complement reference genome to use and which ploidy to use for variant calling. You can use reported sex from sample metadata to make this list or assert the sex chromosome complement using the actual sequencing data using code in the `02_SCC_check` module.

Once custom config files are created, you will need to set up an environment to access the software needed to do the analysis and the sex chromosome complement reference genomes.  We have provided a Docker container and a conda environment that contains both of these to analyze human samples.  

From here, you can use DNA sequencing data to call variants using the `03a_SCC-aware_VariantCalling` module or use RNA sequencing data to quantify gene expression using the `03b_SCC_gene_quantification` module.  These modules are independent of each other, so they can be run in parallel if both DNA and RNA sequencing is available for your samples.

# Setting up your environment to run SCC aware genomics pipelines

## Clone repository

Start by getting a local copy of this Github repository which contains all the code for the modules by navigating to a directory you want to work in and issuing the clone command:

``` 
cd /path/to/local/directory/
git clone https://github.com/SexChrLab/SCC-alignment.git 
```

## Tools used to run workflow code
Our pipelines use Snakemake, a workflow management tool for Python, to streamline and parallelize processes. 

For more information on programming in Python: 
https://www.python.org/

For more information on Snakemake workflows: 
https://snakemake.readthedocs.io/en/stable/

If you use different workflow management systems or other scripting tools, the `shell` portion of the Snakemake rules can be used as a guide on how to issue commands with the appropriate parameters.

## Required packages/software and SCC reference genomes
We have assembled all of the required software needed to run our methods that can be loaded with conda package manager (exported into `SCCalign_v3.yml`). Here we describe how to run from a Docker container which also contains sex chromosome complement reference genomes as well as how to activate the conda environment on its own in a Linux environment.  

For more information on conda package manager including how to install it in the correct operating system: 
https://conda.io/projects/conda/en/latest/index.html

### Docker/Singularity container containing SCC reference genomes

Docker is a utility that allows you to use, store, and share a runtime environment with all the software installed properly. We have created a Docker image that contains SCC reference genomes and transcriptomes used for the analysis modules and conda installed so that required software can be easily loaded and used.  

For more information about Docker containers: 
https://docs.docker.com/

To obtain a copy of the Docker image to run our SCC aware analysis pipelines and have access SCC reference genomes we built, install Docker and use the following command from your terminal:
```
docker pull sbplaisier/omics:1.3
```

This Docker image is about 24 GB because of the reference genomes, so make sure you have enough space available.  

To make sure the image was pulled correctly, list the Docker images and make sure it is there and assigned an ID: 
```
docker images
```

You can familiarize yourself with this Docker container by running it interactively:
```
docker run -it <IMAGE ID>
```

We have installed the packages needed to run our pipelines in the base environment.  For example, you can test `snakemake` by looking at its documentation: 
```
snakemake --help
```

In the working environment, you can see that we are using a Linux environment and have included sex chromosome complement references from CHM13 version 2 (telomere-to-telomere) and HG38 (GRCh38) in the `/references` directory in a docker container: 

```
ls /references
```

To fully perform the analysis, you would attach a volume, meaning that you would make a local directory visible inside the Docker container: 
```
docker run -it -v /path/to/local/directory/:/data -t sbplaisier/omics:1.3
```

Once you attach a volume/bound directory, you can set the reference paths in your custom config to be inside `/reference` and the path to the reads to be inside `/data`.  

The environment that is activated when you first start the Docker container is the `gatk` environment.  To add the other packages needed to run the workflows we have provided, update the `gatk` environment with the environment yml provided in the repository
```
conda env update -f /data/SCCalign_v3.yml --prune
conda activate SCCalign_v3
```

If you are working in a high performance computing cluster or other environment that does not allow use of Docker but instead allows the use of Singularity containers, you can create a Singularity container from a Docker image and run the above steps within the Singularity container: 

```
singularity pull -F SCC_analysis.sif  docker://sbplaisier/omics:1.3
singularity shell -B /path/to/local/directory/:/data SCC_analysis.sif
```

For more information about Singularity (see 'Singularity and Docker' section): 
https://docs.sylabs.io/guides/3.5/user-guide/introduction.html

Now you should be ready to use `snakemake` to run the analysis workflows and have all the packages within the workflow ready to go!

### Installing the conda environment alone
If you would like to simply install the conda environment and generate your own SCC reference genomes, you can install conda to your local environment and then build a conda environment containing the necessary packages for SCC aware analysis using the environment yml file we have provided: 
```
conda env create -f SCC-alignment/SCCalign_v3.yml
```

Once the environment has been created, you can activate it and use all the software used by our workflow: 
```
conda activate SCCalign_v3
```

# Test Data (Genome in a Bottle)

If you would like to use publicly available data to set up and test our pipeline for sex chromosome complement aware genomics analysis techniques, we provide links to data from the Genome in a Bottle project at the National Institute for Standards and Technology (NIST).  You can read more about it here: 
https://www.nist.gov/programs-projects/genome-bottle

In summary, the Genome in a Bottle project aims to provide benchmarking data for the development of standards for genomic analysis. They are rigorously analyzing a specific set of human samples using many sequencing technologies to characterize genomic variants and gene expression differences.  They are making both the data and the results publicly available for the scientific community. 

The samples we will link to are derived from 2 family trios (son, mother, and father).  Our methods do *not* at all require samples to be from related individuals, but we do mark this description to help label samples as having Y chromosomes (sons and fathers) versus not having a Y chromosome (mothers).  Knowing or having the data to find out the sex chromosome complement genotype is essential to running our pipelines.  

We link to specific test files in the readme for `02_SCC_check`, `03a_SCC-aware_VariantCalling`, and `03b_SCC_gene_quantification`.

# Citations

Webster TH, Couse M, Grande BM, Karlins E, Phung TN, Richmond PA, Whitford W, Wilson MA. Identifying, understanding, and correcting technical artifacts on the sex chromosomes in next-generation sequencing data. Gigascience. 2019 Jul 1;8(7):giz074. doi: 10.1093/gigascience/giz074. PMID: 31289836; PMCID: PMC6615978.

Olney KC, Brotman SM, Andrews JP, Valverde-Vesling VA, Wilson MA. Reference genome and transcriptome informed by the sex chromosome complement of the sample increase ability to detect sex differences in gene expression from RNA-Seq data. Biol Sex Differ. 2020 Jul 21;11(1):42. doi: 10.1186/s13293-020-00312-9. PMID: 32693839; PMCID: PMC7374973.

Olney KC, Plaisier SB, Phung TN, Silasi M, Perley L, O'Bryan J, Ramirez L, Kliman HJ, Wilson MA. Sex differences in early and term placenta are conserved in adult tissues. Biol Sex Differ. 2022 Dec 22;13(1):74. doi: 10.1186/s13293-022-00470-y. PMID: 36550527; PMCID: PMC9773522.
