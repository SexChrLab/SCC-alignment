# Sex chromosome complement (SCC) informed alignment pipelines

This repository of pipelines and code from Plaisier et al. 202x shows how to effectively account for sex chromosomes in common human genomic analyses: read alignment and variant calling (SCC-aware_VariantCalling) and gene expression analyses (gene_quantification_RNAseq). 

# Reasons to make the switch to sex-chromsome complement aware genomics methods

Due to the unique biology of the sex chromosomes, accurate measurement of genetic variation and gene expression from the sex chromosomes can be difficult.  As a result, sex chromosomes are often excluded from genomic analyses, even though they are an essential part of the genetic variation within individuals and can have a significant impact on cell biology as it pertains to health and disease.  

The pipelines and resources presented here aim to facilite correct and accurate alignment of short read sequencing data to sex chromosome complement adjusted versions of the human reference genome.  For samples that do not contain a Y chromosome, the Y chromosome sequence will be masked so that no reads are incorrectly assigned to it.  For samples that contain a Y chromosome, the regions of the Y chromosome that have identical matching regions in the X chromosome (psuedoautosomal regions (PARs)) are masked so that reads in this ambiguous regions are not assigned randomly and variant calling on the sex chromosomes is adjusted to match the ploidy.

# Knowing the sex chromosome complement in your samples
In order to do sex chromosome complement aware analysis, it is necessary to know the sex chromosome complement (SCC) present in each of the samples you are working with.  

If you have DNA sequencing data (whole genome sequencing data for example), you can use the provided sample sex chromosome complement checking module (`SCC_check`) to confirm presence or lack of a Y chromosome.  If you do not have DNA sequencing data and are instead working with RNA sequencing data, you can use reported sex for your samples.  If you do not have reported sex, we provide suggestions on how to use gene expression to predict the sex chromosome complement of your samples. 

# Overview of the three SCC-informed analysis modules

1. Custom configuration files (`custom_config`)

In this module, we have provided scripts to create a custom configuration file in JSON format.  You will need specific configs for DNA analysis and RNA analysis.

2. Sex chromosome complement check (`SCC_check`)

In this module, we use whole genome resequencing (WGS) data mapped to a Y PARs-masked reference genome to identify evidence of a Y chromosome in the sequence reads. This information can be used to validate reported sex of the sample metadata or as an independent investigation into the individuals' genotype. We acknowledge that previous iterations of this read mapping depth approach were computationally intensive and were (somewhat justifiably) avoided. This pipeline attempt to bypass these limitations by subsampling the WGS data to 1X coverage prior to alignment to significantly reduce runtime. However, this subsampling restricts the inferential power of the analysis to simply ask, "Are Y chromosome reads present in the sequence data?". If additional information is needed to be inferred from the data (e.g. investigating X chromosome copy number (CN)), you should not use this workflow as it is not suited for this purpose.

3. Sex chromosome complement-aware variant calling (`SCC-aware_Variant Calling`)

In this module, we use whole genome sequencing (`WGS`) data to impute variants across the genome. We use the inferred Sex Chromosome Complement (SCC) from module 1, or sex reported in sample metadata, to assign the appropriate SCC aware reference genome and downstream gentyping criteria for each sample. In short, we map WGS reads using bwa/minimap2 and calculate variants considering the biologically relevant ploidy levels across the genome using GATK.

4. SCC-aware gene expression analyses (`gene_quantification_RNAseq`) 

In this module, we use HISAT2 aligner with a gene quantification algorithm (`featureCounts`) or Salmon pseudoaligner to quantify the expression of genes in samples assayed with RNA sequencing.  

# Test Data (Genome in a Bottle)

If you would like to use publicly available data to set up and test our pipeline for sex chromosome complement aware genomics analysis techniques, we provide links to data from the Genome in a Bottle project at the National Institute for Standards and Technology (NIST).  You can read more about it here: 
https://www.nist.gov/programs-projects/genome-bottle

In summary, the Genome in a Bottle project aims to provide benchmarking data for the development of standards for genomic analysis. They are rigorously analyzing a specific set of human samples using many sequencing technologies to characterize genomic variants and gene expression differences.  They are making both the data and the results publicly available for the scientific community.  

The samples we will link to are from 2 family trios (son, mother, and father).  Our methods do not at all require samples to be from related individuals, but we do mark this description to help label samples as having Y chromosomes (sons and fathers) versus not having a Y chromosome (mothers).  Knowing or having the data to find out the sex chromosome complement genotype is essential to running our pipelines.  


# Setting up your environment to run SCC aware genomics pipelines

## Clone repository

Start by get a local copy of this Github repository by navigating to a directory you want to work in and issuing the clone command:

``` 
cd /path/to/local/directory/
git clone https://github.com/SexChrLab/SCC-alignment.git 
```

## Workflows
Our pipelines use Snakemake, a workflow management tool for Python, to streamline and parallelize processes. 

For more information on programming in Python: 
https://www.python.org/

For more information on Snakemake workflows: 
https://snakemake.readthedocs.io/en/stable/

If you use different workflow management systems or other scripting tools, the `shell` portion of the Snakemake rules can be used as a guide on how to issue commands with the appropriate parameters.

## Required packages/software
We have assembled all of the required software needed to run our methods into a conda environment (`SCCalign_v3.yml`) which we have activated in a Docker container along with sex chromosome complement versions of the two latest releases of the human reference genome (HG38 and CHM13).  Here we describe how to activate the conda environment on its own in a Linux environment as well as how to run from a Docker container which contains sex chromosome complement reference genomes.  

For more information on conda package manager: 
https://conda.io/projects/conda/en/latest/index.html

If you would like to install/load conda to your local environment and build a conda environment containing the necessary packages for SCC aware analysis: 
```
conda env create -f /path/to/local/directory/SCC-alignment/SCCalign_v3.yml
```

## Docker container

Docker is a utility that allows you to use, store, and share a runtime environment with all the software installed properly.  

For more information about Docker containers: 
https://docs.docker.com/


To obtain a copy of the Docker image to run our SCC aware analysis pipelines, install or load Docker and use the following command from your terminal:
```
docker pull sbplaisier/omics:1.3
```

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

In the working environment, you can see that we are using a Linux environment and have included sex chromosome complement references from CHM13 version 2 (telomere-to-telomere) and HG38 (GRCh38) in the /references directory in a docker container: 

```
ls /references
```

To fully perform the analysis, you would attach a volume, meaning that you would make a local directory visible inside the Docker container: 
```
docker run -it -v /path/to/local/directory/:/data -t sbplaisier/omics:1.3
```

Once you attach a volume/bound directory, you can set the reference paths in your custom config to /reference.  Because the reference genomes are included inside the Docker, it is quite large to download (~25 GB).

The environment that is activated when you first start the Docker container is the `gatk` environment.  To add the other packages, update this with the environment yml provided in the repository
```
conda env update -f /data/SCCalign_v3.yml --prune
```

If you are working in a high performance computing cluster that does not allow use of Docker but instead allows the use of Singularity containers, you can create a Singularity container from a Docker image and run the above steps within the Singularity container 

```
singularity pull -F SCC_analysis.sif  docker://sbplaisier/omics:1.3
singularity singularity shell -B /path/to/local/directory/:/data SCC_analysis.sif
```

For more information about Singularity (see 'Singularity and Docker' section): 
https://docs.sylabs.io/guides/3.5/user-guide/introduction.html

Now you should be ready to use `snakemake` to run the analysis workflows and have all the packages within the workflow ready to go!

# Citations

Webster TH, Couse M, Grande BM, Karlins E, Phung TN, Richmond PA, Whitford W, Wilson MA. Identifying, understanding, and correcting technical artifacts on the sex chromosomes in next-generation sequencing data. Gigascience. 2019 Jul 1;8(7):giz074. doi: 10.1093/gigascience/giz074. PMID: 31289836; PMCID: PMC6615978.

Olney KC, Brotman SM, Andrews JP, Valverde-Vesling VA, Wilson MA. Reference genome and transcriptome informed by the sex chromosome complement of the sample increase ability to detect sex differences in gene expression from RNA-Seq data. Biol Sex Differ. 2020 Jul 21;11(1):42. doi: 10.1186/s13293-020-00312-9. PMID: 32693839; PMCID: PMC7374973.

Olney KC, Plaisier SB, Phung TN, Silasi M, Perley L, O'Bryan J, Ramirez L, Kliman HJ, Wilson MA. Sex differences in early and term placenta are conserved in adult tissues. Biol Sex Differ. 2022 Dec 22;13(1):74. doi: 10.1186/s13293-022-00470-y. PMID: 36550527; PMCID: PMC9773522.
