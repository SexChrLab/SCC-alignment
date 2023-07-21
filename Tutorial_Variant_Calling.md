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
wget https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG003.novaseq.pcr-free.40x.R1.fastq.gz?_ga=2.237371716.-57486372.1684265308
wget https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG003.novaseq.pcr-free.40x.R2.fastq.gz?_ga=2.237371716.-57486372.1684265308
wget https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG005.novaseq.pcr-free.40x.R1.fastq.gz?_ga=2.203778004.-57486372.1684265308
wget https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG005.novaseq.pcr-free.40x.R2.fastq.gz?_ga=2.266211275.-57486372.1684265308

wget https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG004.novaseq.pcr-free.40x.R1.fastq.gz?_ga=2.266211275.-57486372.1684265308
wget https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG004.novaseq.pcr-free.40x.R2.fastq.gz?_ga=2.266211275.-57486372.1684265308
wget https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG007.novaseq.pcr-free.40x.R1.fastq.gz?_ga=2.266211275.-57486372.1684265308
wget https://storage.cloud.google.com/brain-genomics-public/research/sequencing/fastq/novaseq/wgs_pcr_free/40x/HG007.novaseq.pcr-free.40x.R2.fastq.gz?_ga=2.229190616.-57486372.1684265308
```

# Clone Repository to get required code

# Setting up environment using Docker to get required software

# Create custom config for DNA data

# Verify sex chromosome complement of samples

# Perform variant calling

