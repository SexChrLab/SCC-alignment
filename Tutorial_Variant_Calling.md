# Tutorial

This page describes how to use the code provided to do sex-aware variant calling using the test data files from Genome in a Bottle using the telomere-to-telomere reference genome (CHM13v2).

# Download DNA sequencing data for testing

Listed in `03a_SCC-aware_VariantCalling` module, download DNA sequencing data for testing.

```

# navigate to your working directory
cd /path/working_directory/

# create a subdirectory called reads that contains the fastq files of your sequencing runs
mkdir reads

# download the test data fastqs
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG003.novaseq.pcr-free.40x.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG003.novaseq.pcr-free.40x.R2.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG005.novaseq.pcr-free.40x.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG005.novaseq.pcr-free.40x.R2.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG004.novaseq.pcr-free.40x.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG004.novaseq.pcr-free.40x.R2.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG007.novaseq.pcr-free.40x.R1.fastq.gz
wget https://storage.googleapis.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG007.novaseq.pcr-free.40x.R2.fastq.gz
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

From the information given by Genome In A Bottle, we see that samples HG003 and HG005 are from males and likely have a Y chromosome, while samples HG004 and HG007 are from females and likely have no Y chromosome.  We are going to use the `02_SCC_check` module to  verify that the sequencing data supports this description and then use the `03a_SCC-aware_VariantCalling` to call variants in these samples based on the SCC-aware versions of the human reference genome.

# Create custom config for DNA data

### Create sample file info table
In order to create a custom configuration file for DNA analysis, we need to create a sample file info table, let's call it `DNA_samples.csv`.

Open a text editor, enter the following lines, and save as `DNA_samples.csv` in your working directory (such `/data` if you mapped a drive and are working in a Docker):

```
HG003_AshkFather,HG003.novaseq.pcr-free.40x.R1.fastq.gz,HG003.novaseq.pcr-free.40x.R2.fastq.gz
HG005_ChinFather,HG005.novaseq.pcr-free.40x.R1.fastq.gz,HG005.novaseq.pcr-free.40x.R2.fastq.gz
HG004_AshkMother,HG004.novaseq.pcr-free.40x.R1.fastq.gz,HG004.novaseq.pcr-free.40x.R2.fastq.gz
HG007_ChinMother,HG007.novaseq.pcr-free.40x.R1.fastq.gz,HG007.novaseq.pcr-free.40x.R2.fastq.gz
```

### Get custom configuration file template
Next, create a template configuration file using `generate_custom_json_DNA.py`.
`FIX: make them input parametters!!!`

First open `01_custom_config\generate_custom_json_DNA.py`, change the following variables, and save:

```
# set this variable to the name of your sample csv
all_sample_ids = "DNA_samples.csv"

# set this variable to the desired output name
out_config_json = "DNA_samples.config.json"

# set this to be the directory where the fastq read are stored
input_directory = '/data/reads/'
```

Once this is done, you can run this script using `python`

```
python generate_custom_json_DNA.py
```

This will create a JSON file in the working directory that you will need to modify to add details specific to your system and experiment.

### Personalizing your custom configuration file

If you open the JSON file in a text editor, you will see entries in JSON format.  Many of these elements are predefined to provide coordinates and options for DNA analysis in the `02_SCC_check` module DNA workflow and the `03a_SCC-aware_VariantCalling` workflow.  

Somewhere in the config JSON you will see entries for every sample listed in your sample table file.  This will contain the input directory you specified in the `fq_path` variable and the paired end read files you specified in the `fq1` and `fq2` variable and other information from the headers of the sequence files stored as read groups.  

### Indicate sex chromosome complement of your samples

A key element of the custom config is indicating which of your samples have a Y chromosome so that the appropriate sex chromosome complement refernece and ploidy are used.  The `ALL_samples` variable contains a list of all samples (first element in list of entries in your sample table).  The `X_samples` list is meant to contain whichever samples in `ALL_samples` do not have a Y chromosome and `Y_samples` is a list meant to contain the samples with a Y chromosome. You can leave these blank before the `02_SCC_check` module and then fill them in after the you have confirmed that which samples have a high relative read depth of chrY, or just go ahead and fill them using the reported sex  first and change them if you see a reason to. See module `02_SCC_check` for more info.

### Check other elements in the configuration file

The other variables you will need to verify or modify to proceed properly with the rest of analysis are:

1) Threads: number of threads you would like to use for multiprocessing.
```
"-----------------Comment_RUN_INFO-----------------": "Run specifications",
    "threads": "6",
```

2) Paths to output directories for results of the `02_SCC_check` module
```
    "-----------------Comment_RESULTS_DIRECTORY-----------------": "Output directory for relative chrY read depth (indexcov) results.",
    "indexcov_dir": "/data/sex_check/indexcov",
```

3) Path to input directory containing your sequencing read files
```
    "-----------------Comment_FASTQ_DIRECTORIES-----------------": "Directory containing FASTQ files.",
    "reads": "/data/reads",
```

4) Path to output directories for your variant calling results
```
    "-----------------Comment_VCF_DIRECTORIES-----------------": "Output directories for called read vcf file(s).",
    "haplotyped_samples": "/data/SCC_snakemake/haplotyped_vcfs/",
    "genotyped_samples": "/data/SCC_snakemake/genotyped_vcfs/",
```

5) Paths to the SCC reference genomes: default values are set to paths in the Docker container
```
    "-----------------Comment_CLEANED_REFERENCE_GENOME-----------------": "Reference genome files.",
    "CHM13_X": "/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa",
    "CHM13_Y": "/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa",
    "GRCh38_X": "/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XX.fa",
    "GRCh38_Y": "/references/GENCODE/GRCh38.p12.genome.XY/GRCh38.p12.genome.XY.fa",
```

# Verify sex chromosome complement of samples

Now that your custom config file for DNA samples is created, we can use the `02_SCC_check` module DNA workflow to confirm or identify whether or not each sample has a Y chromosome (see readme for this module for more information).

First, navigate to the `02_SCC_check` directory, open `02_SCC_check/SCC-check.snakefile` in a text editor, and confirm that the name and path of your config JSON is set correctly in the `configfile` variable:
```
configfile: "DNA_samples.config.json"
```
The rest of the necessary information will be filled in from the config JSON.  

All of our workflows are written in Snakemake.  For more information on Snakemake, see the main readme for this repository

To make sure everything is setup correctly, do a dry run of the `SCC-check.snakefile` Snakemake workflow for subsampling the reads and determining the relative read depth of chrY using the `indexcov` application: 
```
snakemake -np -s 02_SCC_check/SCC-check.snakefile
```

Once any errors are fixed, go ahead and run the Snakemake workflow for subsampling the reads and determining the relative read depth of chrY using the `indexcov` application: 
```
snakemake -s 02_SCC_check/SCC-check.snakefile --rerun-incomplete
```

Or if you are running on an HPC cluster with Slurm job submission manager, you may use something to this effect: 
```
#!/bin/bash
#SBATCH --job-name=sexcheck  # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=ALL           # notifications for job done & fail
#SBATCH --mail-user=user@domain.edu # send-to address
#SBATCH -t 1-00:00
#SBATCH -p serial
#SBATCH -q public
#SBATCH -n 2

source activate SCCalign_v3
cd /data/SCC-alignment-main/02_SCC_check
snakemake -s SCC-check.snakefile -j 20 --rerun-incomplete --latency-wait=60

```

Upon successful completion of this workflow, you will see results in the `indexcov_dir` you specified in your custom config file.  This will include HTML files you can use to browse your results to see if the chrY read depth indicates the presence of chrY.  There should also be an output file called `indexcov-indexcov.ped`.  This file has specific information plotted in the HTML files.  We provided a simple script to make an easy look up table for your samples if you would like to use it.  

To do so, copy `02_SCC_check\inferred_SCC.py` into the directory containing your `indexcov` results (`indexcov-indexcov.ped`), and run with Python:
```
python inferred_SCC.py
```

```

# Perform variant calling

