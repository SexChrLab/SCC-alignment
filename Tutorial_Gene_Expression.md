# Tutorial for sex chromosome complement aware gene expression quantification

This page describes how to use the code provided to do sex chromosome complement aware gene expression quantification using the RNA sequencing data from Genome in a Bottle.  This tutorial will describe how to use two methods for gene expression quantification-- full alignment to sex chromosome complement reference genomes and feature counting as a means of quantification, and pseudoalignment and quantification to sex chromosome complement reference transcriptome.  The first is more comprehensive but takes more time and computational power, the second is faster and so is used quite frequently for gene expression analysis.

# Download RNA sequencing data for testing

Listed in `03b_gene_quantification_RNAseq` module, download RNA sequencing data for testing.

```

# navigate to your working directory
cd /path/working_directory/

# create a subdirectory called reads that contains the fastq files of your sequencing runs
mkdir reads

# download the test data fastqs
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg004_gm24143.mrna.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg004_gm24143.mrna.R2.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg002_gm24385.mrna.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/rna/illumina/mrna/hg002_gm24385.mrna.R2.fastq.gz
```

# Clone Repository to get required code

```
git clone https://github.com/SexChrLab/SCC-alignment.git
```

# Setting up environment using Docker to get required software

```
# pull Docker with reference genome sequences
docker pull sbplaisier/omics:1.3

# run Docker iteratively 
docker run -it -v /path/working_directory:/data -t sbplaisier/omics:1.3

conda env update -f /data/SCC-alignment-main/SCCalign_v3.yml --prune
conda activate SCCalign_v3
```

# Consider sex of the samples

From the information given by Genome In A Bottle, we see that samples HG002 is from a male and likely has a Y chromosome, while sample HG004 is from a female and likely has no Y chromosome.  We are going to use the `02_SCC_check` module for RNA sequencing data to verify whether the gene expression profile of these samples supports this description and then use the `03b_gene_quantification_RNAseq` module to measure gene expression aligning to the SCC-aware versions of the human reference genome. This tutorial will walk you through two methods of gene expression quantification.  

# Create custom config for RNA data

In order to provide the gene expression quantification tools with the information needed, we must create a custom configuration file. 

### Create sample file info table
In order to create a custom configuration file for RNA sequencing analysis, we need to create a sample file info table, let's call it `RNA_samples.csv`.

Open a text editor, enter the following lines, and save as `RNA_samples.csv` in your working directory (such `/data` if you mapped a drive and are working in a Docker):

```
HG004_AshkMother,hg004_gm24143.mrna.R1.fastq.gz,hg004_gm24143.mrna.R2.fastq.gz
HG002_AshkSon,hg002_gm24385.mrna.R1.fastq.gz,hg002_gm24385.mrna.R2.fastq.gz

```

### Get custom configuration file template
Next, create a template configuration file using `generate_custom_json_RNA.py`.
`FIX: make them input parametters!!!`

First open `01_custom_config\generate_custom_json_RNA.py`, check or change the following variables, and save:

```
# set this variable to the name of your sample csv
all_sample_ids = "RNA_samples.csv"

# set this variable to the desired output name
out_config_json = "RNA_samples.config.json"

# set this to be the directory where the fastq read are stored
input_directory = '/data/reads/'
```

Once this is done, you can run this script using `python`

```
python generate_custom_json_RNA.py
```

This will create a JSON file in the working directory that you will need to modify to add details specific to your system and experiment.

### Personalizing your custom configuration file

If you open the JSON file in a text editor, you will see entries in JSON format.  Many of these elements are predefined to provide coordinates and options for RNA analysis in the `02_SCC_check` module RNA workflow and the `03b_gene_quantification_RNAseq` workflows.  

Somewhere in the config JSON you will see entries for every sample listed in your sample table file.  This will contain the input directory you specified in the `fq_path` variable and the paired end read files you specified in the `fq1` and `fq2` variable and other information from the headers of the sequence files stored as read groups.  

### Indicate sex chromosome complement of your samples

A key element of the custom config is indicating which of your samples have a Y chromosome so that the appropriate sex chromosome complement refernece and ploidy are used.  The `ALL_samples` variable contains a list of all samples (first element in list of entries in your sample table).  The `X_samples` list is meant to contain whichever samples in `ALL_samples` do not have a Y chromosome and `Y_samples` is a list meant to contain the samples with a Y chromosome. You can leave these blank before the `02_SCC_check` module and then fill them in after the you have confirmed that which samples have a high relative read depth of chrY, or just go ahead and fill them using the reported sex  first and change them if you see a reason to. See module `02_SCC_check` for more info.

### Check other elements in the configuration file

The other variables you will need to verify or modify to proceed properly with the rest of analysis are:

1) Paths to your reference genome indices and genome annotation files (located inside /references inside the Docker container):

```
"-----------------Comment_Transcriptome_Index-----------------": "This section specifies the location of the index files for gene quantification tools",
    "HG38_Transcriptome_Index_HISAT_Path_female": "your/path/to/female/HG38/HISAT/index/",
    "HG38_Transcriptome_Index_HISAT_Path_male": "your/path/to/male/HG38/HISAT/index/",
    "HG38_Transcriptome_Index_SALMON_Path_female": "your/path/to/female/HG38/SALMON/index/",
    "HG38_Transcriptome_Index_SALMON_Path_male": "your/path/to/male/HG38/SALMON/index/",
    "CHM13_Transcriptome_Index_HISAT_Path_female": "your/path/to/female/CHM13/HISAT/index/",
    "CHM13_Transcriptome_Index_HISAT_Path_male": "your/path/to/male/CHM13/HISAT/index/",
    "CHM13_Transcriptome_Index_SALMON_Path_female": "your/path/to/female/CHM13/SALMON/index/",
    "CHM13_Transcriptome_Index_SALMON_Path_male": "your/path/to/male/CHM13/SALMON/index/",

    "-----------------Comment_Annotation-----------------": "This section specifies the location of the genome annotation file",
    "HG38_annotation_path": "your/path/to/HG38/annotation/gtf/gff",
    "CHM13_annotation_path": "your/path/to/CHM13/annotation/gtf/gff",
```

2) Stranded-ness of your sequencing files for Hisat2 if using (F = forward, R = reverse, options provided as or statement, keep only one, link to manual in module readme)

```
    "-----------------Comment_HISAT_Parameters-----------------": "This section specifies parameters for HISAT analysis of your data",
    "HISAT_strandedness": "F|R|RF|FR",
```

3) Set library type for Salmon if using (options provided, keep only one, link to manual in module readme)
```
    "-----------------Comment_SALMON_Parameters-----------------": "This section specifies parameters for SALMON analysis of your data",
    "SALMON_libtype": "I|O|M|S|U|F|R",
```

# Verify sex chromosome complement of samples

Now that your custom config file for RNA samples is created, we can use the `02_SCC_check` module RNA workflow to confirm or identify whether or not each sample has a Y chromosome (see readme for this module for more information).

First, navigate to the `02_SCC_check` directory, open `02_SCC_check/SCC-check.snakefile` in a text editor, and confirm that the name and path of your config JSON is set correctly in the `configfile` variable:
```
configfile: "RNA_samples.config.json"
```
The rest of the necessary information will be filled in from the config JSON.  

`FINISH!!!!: gene expression of sex chr genes`
`TO-DO!!!!: make the piped version the only one`

# Quantify gene expression using alignment and feature counts (Hisat2 and featureCounts)

Once you have indicated the sex chromosome complement of your samples and completed your RNA custom config, you are ready to quantify gene expression!  

Quantification using full alignment and feature counting is implemented in `03b_gene_quantification_RNAseq/rnaseq_data_processing_hisat2_piped.snakefile`

Before running this workflow, perform the following checks:
1) copy the custom RNA config file to the `03b_gene_quantification_RNAseq` directory: `cp 01_custom_config/RNA_samples.config.json 03b_gene_quantification_RNAseq/`
2) open `rnaseq_data_processing_hisat2_piped.snakefile` in a text editor
3) make sure the `configfile` variable is set to the name of your custom config JSON
4) make sure the `X_genome` and `Y_genome` variables are set to the SCC version of the reference genome according to your preference.  By default, these are set to the CHM13v2 (telomere-to-telomere) SCC genome reference sequences, but if you prefer to use SCC HG38 (GRCh38) references, uncomment the lines refering to the path to those set in the custom config and comment the lines indicating the CHM13v2 references like this:
```
# uncomment these if you want to run with the HG38 (GRCh38) human genome
X_genome = config["HG38_Transcriptome_Index_HISAT_Path_female"]
Y_genome = config["HG38_Transcriptome_Index_HISAT_Path_male"]
genome_annotation = config["HG38_annotation_path"]

# uncomment these if you want to run with the CHM13 (telomere-to-telomere) human genome
#X_genome = config["CHM13_Transcriptome_Index_HISAT_Path_female"]
#Y_genome = config["CHM13_Transcriptome_Index_HISAT_Path_male"]
#genome_annotation = config["CHM13_annotation_path"]
```

Once these checks and preferences have been indicated, you are ready to run the workflow. 

To make sure everything is setup correctly, do a dry run for both workflows and resolve any warnings or errors that come up: 
```
snakemake -np -s rnaseq_data_processing_hisat2_piped.snakefile
```

Once any errors are fixed, go ahead and run the Snakemake workflow for subsampling the reads and determining the relative read depth of chrY using the `indexcov` application: 
```
snakemake -s rnaseq_data_processing_hisat2_piped.snakefile --rerun-incomplete
```

Or if you are running on an HPC cluster with Slurm job submission manager, you may use something to this effect: 
```
#!/bin/bash
#SBATCH --job-name=SNPs  # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=user@domain.edu # send-to address
#SBATCH -t 1-00:00
#SBATCH -p serial
#SBATCH -q public
#SBATCH -n 2

source activate SCCalign_v3
cd /data/SCC-alignment-main/03b_gene_quantification_RNAseq
snakemake -s rnaseq_data_processing_hisat2_piped.snakefile -j 20 --rerun-incomplete --latency-wait=60

```

Once the workflow is completed, you will see quantified gene expression files in the `featureCounts/` directory.  See readme for `03b_gene_quantification_RNAseq` for more information.  

# Quantify gene expression using pseudoalignment (Salmon)

Quantification using full alignment and feature counting is implemented in `03b_gene_quantification_RNAseq/rnaseq_data_processing_salmon.snakefile`

Before running this workflow, perform the following checks:
1) copy the custom RNA config file to the `03b_gene_quantification_RNAseq` directory: `cp 01_custom_config/RNA_samples.config.json 03b_gene_quantification_RNAseq/`
2) open `rnaseq_data_processing_hisat2_piped.snakefile` in a text editor
3) make sure the `configfile` variable is set to the name of your custom config JSON
4) make sure the `X_genome` and `Y_genome` variables are set to the SCC version of the reference genome according to your preference.  By default, these are set to the CHM13v2 (telomere-to-telomere) SCC genome reference sequences, but if you prefer to use SCC HG38 (GRCh38) references, uncomment the lines refering to the path to those set in the custom config and comment the lines indicating the CHM13v2 references like this:
```
# uncomment these to use HG38 (GRCh38) SCC reference genomes
X_genome = config["HG38_Transcriptome_Index_SALMON_Path_female"]
Y_genome = config["HG38_Transcriptome_Index_SALMON_Path_male"]

# uncomment these to use CHM13 (telomere-to-telomere) SCC reference genomes
#X_genome = config["CHM13_Transcriptome_Index_SALMON_Path_female"]
#Y_genome = config["CHM13_Transcriptome_Index_SALMON_Path_male"]
```

Once these checks and preferences have been indicated, you are ready to run the workflow. 

To make sure everything is setup correctly, do a dry run for both workflows and resolve any warnings or errors that come up: 
```
snakemake -np -s rnaseq_data_processing_salmon.snakefile
```

Once any errors are fixed, go ahead and run the Snakemake workflow for subsampling the reads and determining the relative read depth of chrY using the `indexcov` application: 
```
snakemake -s rnaseq_data_processing_salmon.snakefile --rerun-incomplete
```

Or if you are running on an HPC cluster with Slurm job submission manager, you may use something to this effect: 
```
#!/bin/bash
#SBATCH --job-name=SNPs  # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=user@domain.edu # send-to address
#SBATCH -t 1-00:00
#SBATCH -p serial
#SBATCH -q public
#SBATCH -n 2

source activate SCCalign_v3
cd /data/SCC-alignment-main/03b_gene_quantification_RNAseq
snakemake -s rnaseq_data_processing_salmon.snakefile -j 20 --rerun-incomplete --latency-wait=60

```

Once the workflow is completed, you will see quantified gene expression files in the `quantified_rna_salmon/` directory.  See readme for `03b_gene_quantification_RNAseq` for more information.  


We hope that this tutorial helps you to see the various parts of the sex chromosome complement aware gene expression pipeline and allows you to make modifications for your own use.  If you have any issues, please contact the Sex Chromosome Lab for help (http://www.sexchrlab.org/).  
