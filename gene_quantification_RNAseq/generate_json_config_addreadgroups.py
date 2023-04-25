import json
from collections import defaultdict
import os
import gzip

# Specifying inputs:
all_sample_ids = "samples.csv"
out_config_json = "RNAseq_samples.json"
input_directory = '/data/CEM/wilsonlab/downloaded_data/GIAB/IlluminaRNAseq/'

data = {}

data['all_rna_samples'] = []
dna_samples_fastq_path = defaultdict(list)

with open(all_sample_ids, 'r') as f:
    for line in f:
        items = line.rstrip('\n').split(',')
        data['all_rna_samples'].append(items[0])
        dna_samples_fastq_path[items[0]].append(items[1])
        dna_samples_fastq_path[items[0]].append(items[2])
        R1_file = items[1]
        R1_filesplit = R1_file.split("_")

for sample in data['all_rna_samples']:
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

