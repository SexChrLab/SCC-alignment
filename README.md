# SCC-alignment

This repository hosts a walkthrough, or tutorial, of pipelines and code from Plaisier et al. 202x, where we show how to effectively account for sex chromosomes in common human genomic analyses (sample sex chromosome complement checking (SCC_check), read alignment and variant calling (SCC-aware_VariantCalling), differential expression analyses (gene_quantification_RNAseq). This tutorial uses ```snakemake``` as a workflow manager to streamline and parallelize these processes. At a minimum, we recommend using the provided conda environment ```SCCalign_v3.yml``` to run these analyses. We also provide this working environment and associated data files in a docker computing enviroment (https://www.docker.com) and provide information to adapt this for singularity (https://sylabs.io/docs) for individuals that don't possess root priviledges on their respective computing ecosystem.

# Installing the conda evironment for analysis

Start by cloning the repository to your working directory:

``` 
git clone https://github.com/SexChrLab/SCC-alignment.git 
```

We suggest using ```mamba``` (https://github.com/mamba-org/mamba) to install the environment quickly, but can also be done by replacing ```mamba``` with ```conda``` in the folowing line:
```
mamba env create --name SCCalign_v3 --file=SCCalign/SCCalign_v3.yml 
```

# Installing the docker or singularity images

We have provided sex chromosome complement references from CHM13 version 2 (telomere-to-telomere) and HG38 (GRCh38) in the /references directory in a docker container: 

```
docker pull sbplaisier/omics:1.3
```

```
singularity pull -F SCC_analysis.sif  docker://sbplaisier/omics:1.3
```

Once you attach a volume/bound directory, you can set the paths in your config to /reference.  Because the reference genomes are included inside the Docker, it is quite large to download (~25 GB).

# Snakemake config file generation

We have provided scripts to create a custom configuration file in JSON format (see `custom_config`).  You will need specific configs for DNA analysis and RNA analysis.

# Overview of the three SCC-informed analysis modules

1. Sex chromosome complement check (`SCC_check`)

In this module, we use whole genome resequencing (WGS) data mapped to a Y PARs-masked reference genome to identify evidence of a Y chromosome in the sequence reads. This information can be used to validate reported sex of the sample metadata or as an independent investigation into the individuals' genotype. We acknowledge that previous iterations of this read mapping depth approach were computationally intensive and were (somewhat justifiably) avoided. This pipeline attempt to bypass these limitations by subsampling the WGS data to 1X coverage prior to alignment to significantly reduce runtime. However, this subsampling restricts the inferential power of the analysis to simply ask, "Are Y chromosome reads present in the sequence data?". If additional information is needed to be inferred from the data (e.g. investigating X chromosome copy number (CN)), you should not use this workflow as it is not suited for this purpose.

2. Sex chromosome complement-aware variant calling (`SCC-aware_Variant Calling`)

In this module, we use whole genome sequencing (`WGS`) data to impute variants across the genome. We use the inferred Sex Chromosome Complement (SCC) from module 1, or sex reported in sample metadata, to assign the appropriate reference genome and downstream gentyping criteria for each sample. In short, we map WGS reads using bwa/minimap2 and calculate variants considering the biologically relevant ploidy levels across the genome using GATK.

3. SCC-aware gene expression analyses (gene_quantification_RNAseq) 

In this module, we use HISAT2 aligner with a gene quantification algorithm (featureCounts) or Salmon pseudoaligner to quantify the expression of genes in samples assayed with RNA sequencing.  

# Additional resources

Our workflows are implemented in Snakemake, a workflow management tool for Python.  

For more information on programming in Python: 
https://www.python.org/

For more information on Snakemake workflows: 
https://snakemake.readthedocs.io/en/stable/

