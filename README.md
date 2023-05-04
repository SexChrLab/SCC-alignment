# SCC-alignment

This repository hosts a walkthrough, or tutorial, of pipelines and code from Plaisier et al. 202x, where we show how to effectively account for sex chromosomes in common human genomic analyses (sample sex chromosome complement checking (SCC_check), read alignment and variant calling (SCC-aware_VariantCalling), differential expression analyses (gene_quantification_RNAseq). This tutorial uses ```snakemake``` as a workflow manager to streamline and parallelize these processes. At a minimum, we recommend using the provided conda environment ```SCCalign_v3.yml``` to run these analyses. We also provide this working environment and associated data files in a docker computing enviroment (https://www.docker.com) and provide information to adapt this for singularity (https://sylabs.io/docs) for individuals that don't possess root priviledges on their respective computing ecosystem.

# Installing the conda evironment for analysis

We suggest using ```mamba``` (https://github.com/mamba-org/mamba) is install, but can also be done by replacing ```mamba``` with ```conda```
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
