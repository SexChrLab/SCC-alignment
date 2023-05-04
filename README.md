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

```
docker pull
```

```
singularity pull
```

# Snakemake config file generation

Add your working directory path to the JSON file, e.g. using sed (note the use of "\\" prior to symbols is necessary for literal interpretation)

```
sed -i 's/\~/\/scratch\/working_dir/g' SCC-analysis_info.json
python [updated] generate_json_config_addreadgroups.py 
cat SCC-analysis_info.json SCC-analysis_samples.json > SCC-analysis_config.json
```

To check that the config is formatted properly, we suggest using a JSON validator tool, such as https://jsonlint.com/ 

# Overview of the three SCC-informed analysis modules

1. Sex chromosome complement check (SCC_check)

In this module, we use whole genome resequencing (WGS) data mapped to a Y PARs-masked reference genome to identify evidence of a Y chromosome in the sequence reads. This information can be used to validate reported sex of the sample metadata or as an independent investigation into the individuals' genotype. We acknowledge that previous iterations of this read mapping depth approach were computationally intensive and were (somewhat justifiably) avoided. This pipeline attempt to bypass these limitations by subsampling the WGS data to 1X coverage prior to alignment to significantly reduce runtime. However, this subsampling restricts the inferential power of the analysis to simply ask, "Are Y chromosome reads present in the sequence data?". If additional information is needed to be inferred from the data (e.g. investigating X chromosome copy number (CN)), you should not use this workflow as it is not suited for this purpose.

2. Sex chromosome complement-aware variant calling (SCC-aware_Variant Calling)

In this module, we use whole genome resequencing (WGS) data to impute variants across the genome. We use the inferred Sex Chromosome Complement (SCC) from modeule 1, or sex reported in sample metadata, to assign the appropriate reference genome and downstream gentyping criteria for each sample. In short, we map WGS reads using bwa/minimap2 and calculate variants considering the biologically relevant ploidy levels across the genome using GATK.

3. SCC-aware gene expression analyses (gene_quantification_RNAseq) 

