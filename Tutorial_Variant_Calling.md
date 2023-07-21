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

# Create custom config for DNA data

# Verify sex chromosome complement of samples

# Perform variant calling

