import json
from collections import defaultdict
import os
import gzip
import argparse

# Specifying inputs:
all_sample_ids = "DNA_samples.csv"
out_config_json = "DNA_samples.config.json"
input_directory = '/data/CEM/wilsonlab/projects/XYAlign_Tutorial/reads/'

'''
parser = argparse.ArgumentParser(description='Creates custom configuration JSON for SCC aware gene quantification')
parser.add_argument('-s','--sampleTable', help='Enter a csv containing columns for sample name, path to file directory, forward read file name, reverse read file name, and optionally has_Y', required=True)
parser.add_argument('-j','--outputJson', help='Enter path to output JSON', required=True)

args = vars(parser.parse_args())
'''

# create a dictionary to hold the info we will dump into a config json
data = {}

# add place holders for various information needed for variant calling
data["-----------------Comment_RUN_INFO-----------------"] =  "Run specifications"
data["threads"] =  "6"

data["-----------------Comment_RESULTS_DIRECTORY-----------------"] =  "Directories containing indexcov results."
data["indexcov_dir"] =  "/data/CEM/wilsonlab/projects/XYAlign_Tutorial/sex_check/indexcov"

data["-----------------Comment_SAMPLE_INFO_CONTROL-----------------"] =  "All control (known) samples to be analyzed."
data["samples_CONTROL"] =  ["HG002", "HG004"]

data["-----------------Comment_FASTQ_DIRECTORIES-----------------"] =  "Directory containing FASTQ files."
data["reads"] =  "/data/CEM/wilsonlab/projects/XYAlign_Tutorial/reads"

data["-----------------Comment_VCF_DIRECTORIES-----------------"] =  "Directories containing called read vcf file(s)."
data["haplotyped_samples"] =  "/scratch/splaisie/SCC-alignment-main_new/SCC-alignment-main/SCC_snakemake/haplotyped_vcfs/"
data["genotyped_samples"] =   "/scratch/splaisie/SCC-alignment-main_new/SCC-alignment-main/SCC_snakemake/genotyped_vcfs/"
	
data["-----------------Comment_CLEANED_REFERENCE_GENOME-----------------"] =  "Reference genome files."
data["CHM13_X"] =  "/data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa"
data["CHM13_Y"] =  "/data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa"
data["GRCh38_X"] =  "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XXonly/GRCh38.p12.genome.XX.fa"
data["GRCh38_Y"] =  "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome.XY/GRCh38.p12.genome.XY.fa"

data["-----------------Comment_CHROMOSOME_INFO-----------------"] =  "Chromosomes to run"
data["chrX"] =  ["chrX"]
data["chrY"] =  ["chrY"]
data["directory"] =  ["X", "Y"]
data["autosomes"] =  ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
data["XX_diploid"] =  ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
data["XY_diploid"] =  ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "PAR1", "PAR2"]
data["XY_haploid"] =  ["chrX_nonPAR", "chrY"]
data["all_chromosomes"] =  ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
data["all_regions"] =  ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX_PAR1", "chrX_PAR2", "chrX_nonPAR", "chrY"]

data["-----------------Comment_REFERENCE_COORDINATES-----------------"] =  "Reference genome feature coordinates."
data["GRCh38_PAR1"] =  ["chrX:10001-2781479"]
data["GRCh38_PAR2"] =  ["chrX:155701383-156030895"]
data["GRCh38_nonPAR"] =  ["chrX:2781480-155701382"]
data["CHM13_PAR1"] =  ["chrX:1-2394410"]
data["CHM13_PAR2"] =  ["chrX:153925834-154259566"]
data["CHM13_nonPAR"] =  ["chrX:2394411-153925833"]

data["-----------------Comment_diploid_FILTERING-----------------"] =  "DIPLOID: GATK genotype filtering options (XX samples)."
data["diploid_AN"] =  "4.0"
data["diploid_MQ"] =  "10.0"
data["diploid_QD"] =  "7.0"
data["diploid_DP1"] =  "10.0"
data["diploid_DP2"] =  "1000.0"

data["-----------------Comment_haploid_XY_FILTERING-----------------"] =  "HAPLOID: GATK genotype filtering options (XY samples)."
data["haploid_AN"] =  "2.0"
data["haploid_MQ"] =  "5.0"
data["haploid_QD"] =  "4.0"
data["haploid_DP1"] =  "5.0"
data["haploid_DP2"] =  "500.0"

data["-----------------Comment_Sample_Info-----------------"] = "This section lists the samples that are to be analyzed"

# create entries for samples that are to be analyzed

data['ALL_samples'] = []
data['X_samples'] = ["select samples with X chromosomes only from ALL_samples"]
data['Y_samples'] = ["select samples with Y chromosomes from ALL_samples"]

dna_samples_fastq_path = defaultdict(list)

with open(all_sample_ids, 'r') as f:
    for line in f:
        items = line.rstrip('\n').split(',')
        data['ALL_samples'].append(items[0])
        dna_samples_fastq_path[items[0]].append(items[1])
        dna_samples_fastq_path[items[0]].append(items[2])
        R1_file = items[1]
        R1_filesplit = R1_file.split("_")

for sample in data['ALL_samples']:
    read_group_info = {}
    fq_path = input_directory
    fq_1 = dna_samples_fastq_path[sample][0]
    fq_2 = dna_samples_fastq_path[sample][1]
    R1_file = fq_1 
    R1_filesplit = R1_file.split("_")

    # find pu
    with gzip.open(os.path.join(fq_path,fq_1), 'rt', encoding="utf8", errors='ignore') as f:
        first_line = f.readline()
        items = first_line.split(':')
        if (len(items) >= 3):
            pu = items[2] + '.' + items[3]
        else:
            pu = sample

    read_group_info[sample] = {
        'fq_path': fq_path,
        'fq1': fq_1,
        'fq2': fq_2,
        'ID': sample,
        'SM': sample,
        'LB': sample,
        'PU': pu,
        'PL': 'Illumina',
        'individual': sample,
        'tissue': 'NA'
    }

    data.update(read_group_info)

with open(out_config_json, 'w') as outfile:
    json.dump(data, outfile, indent=4)

