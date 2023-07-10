import json
from collections import defaultdict
import os
import gzip
import argparse

# Specifying inputs:
all_sample_ids = "RNA_samples.csv"
out_config_json = "RNAseq_samples.config.json"
input_directory = '/data/CEM/wilsonlab/downloaded_data/GIAB/IlluminaRNAseq/'

'''
parser = argparse.ArgumentParser(description='Creates custom configuration JSON for SCC aware gene quantification')
parser.add_argument('-s','--sampleTable', help='Enter a csv containing columns for sample name, path to file directory, forward read file name, reverse read file name, and optionally has_Y', required=True)
parser.add_argument('-j','--outputJson', help='Enter path to output JSON', required=True)

args = vars(parser.parse_args())
'''
# create a dictionary to hold the info we will dump into a config json
data = {}

# add place holders for various information needed for gene quantification tools
data["-----------------Comment_Transcriptome_Index-----------------"] = "This section specifies the location of the index files for gene quantification tools"
data["HG38_Transcriptome_Index_HISAT_Path_female"] = "your/path/to/female/HG38/HISAT/index/"
data["HG38_Transcriptome_Index_HISAT_Path_male"] = "your/path/to/male/HG38/HISAT/index/"
data["HG38_Transcriptome_Index_SALMON_Path_female"] = "your/path/to/female/HG38/SALMON/index/"
data["HG38_Transcriptome_Index_SALMON_Path_male"] = "your/path/to/male/HG38/SALMON/index/"
data["CHM13_Transcriptome_Index_HISAT_Path_female"] = "your/path/to/female/CHM13/HISAT/index/"
data["CHM13_Transcriptome_Index_HISAT_Path_male"] = "your/path/to/male/CHM13/HISAT/index/"
data["CHM13_Transcriptome_Index_SALMON_Path_female"] = "your/path/to/female/CHM13/SALMON/index/"
data["CHM13_Transcriptome_Index_SALMON_Path_male"] = "your/path/to/male/CHM13/SALMON/index/"

data["-----------------Comment_Annotation-----------------"] = "This section specifies the location of the genome annotation file"
data["HG38_annotation_path"] = "your/path/to/HG38/annotation/gtf/gff"
data["CHM13_annotation_path"] = "your/path/to/CHM13/annotation/gtf/gff"

data["-----------------Comment_HISAT_Parameters-----------------"] = "This section specifies parameters for HISAT analysis of your data"
data["HISAT_strandedness"] = "F|R|RF|FR"

data["-----------------Comment_SALMON_Parameters-----------------"] = "This section specifies parameters for SALMON analysis of your data"
data["SALMON_libtype"] = "I|O|M|S|U|F|R"

data["-----------------Comment_Sample_Info-----------------"] = "This section lists the samples that are to be analyzed"

# create entries for samples that are to be analyzed

data['ALL_samples'] = []
data['X_samples'] = ["select samples with X chromosomes only from ALL_samples"]
data['Y_samples'] = ["select samples with Y chromosomes from ALL_samples"]

rna_samples_fastq_path = defaultdict(list)

with open(all_sample_ids, 'r') as f:
    for line in f:
        items = line.rstrip('\n').split(',')
        data['ALL_samples'].append(items[0])
        rna_samples_fastq_path[items[0]].append(items[1])
        rna_samples_fastq_path[items[0]].append(items[2])
        R1_file = items[1]
        R1_filesplit = R1_file.split("_")

for sample in data['ALL_samples']:
    read_group_info = {}
    fq_path = input_directory
    fq_1 = rna_samples_fastq_path[sample][0]
    fq_2 = rna_samples_fastq_path[sample][1]
    R1_file = fq_1 
    R1_filesplit = R1_file.split("_")

    # find pu
    with gzip.open(os.path.join(fq_path,fq_1), 'rt', encoding="utf8", errors='ignore') as f:
        first_line = f.readline()
        items = first_line.split(':')
        pu = items[2] + '.' + items[3]

    read_group_info[sample] = {
        'fq_path': fq_path,
        'fq_1': fq_1,
        'fq_2': fq_2,
        'ID': sample,
        'SM': sample,
        'LB': sample,
        'PU': pu,
        'PL': 'Illumina'
    }

    data.update(read_group_info)

with open(out_config_json, 'w') as outfile:
    json.dump(data, outfile, indent=4)

