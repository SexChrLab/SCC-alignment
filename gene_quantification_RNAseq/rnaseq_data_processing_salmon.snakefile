import os

configfile: "RNAseq_samples.json"

# path to salmon indexes for sex chromosome complement genome

# uncomment these to use HG38 (GRCh38) SCC reference genomes
#X_genome = config["HG38_Transcriptome_Index_SALMON_Path_female"]
#Y_genome = config["HG38_Transcriptome_Index_SALMON_Path_male"]

# uncomment these to use CHM13 (telomere-to-telomere) SCC reference genomes
X_genome = config["CHM13_Transcriptome_Index_SALMON_Path_female"]
Y_genome = config["CHM13_Transcriptome_Index_SALMON_Path_male"]

rule all:
   input:
        # salmon psuedoalignment rules
        expand("quantified_rna_salmon/{sample_name}_salmon_quant_XY/", sample_name = config["Y_samples"]),
        expand("quantified_rna_salmon/{sample_name}_salmon_quant_XX/", sample_name = config["X_samples"])

rule salmon_paired_males:
    input:
        FASTQ1 = "fastq_files_rna/{sample_name}_R1.fastq.gz",
        FASTQ2 = "fastq_files_rna/{sample_name}_R2.fastq.gz"
    output:
        out_1 =  directory("quantified_rna_salmon/"+"{sample_name}_salmon_quant_XY/")
    params:
        SALMON_Index_male = Y_genome,
        libtype = "ISR", # LIBTYPE for paired reads
    shell:
        "salmon quant -i {params.SALMON_Index_male} -l {params.libtype} -1 {input.FASTQ1} -2 {input.FASTQ2} --validateMappings -o {output.out_1}"

rule salmon_paired_females:
    input:
        FASTQ1 = "fastq_files_rna/{sample_name}_R1.fastq.gz",
        FASTQ2 = "fastq_files_rna/{sample_name}_R2.fastq.gz"
    output:
        out_1 = directory("quantified_rna_salmon/{sample_name}_salmon_quant_XX/")
    params:
        SALMON_Index_female = X_genome,
        libtype = "ISR", # LIBTYPE for paired reads
    shell:
        "salmon quant -i {params.SALMON_Index_female} -l {params.libtype} -1 {input.FASTQ1} -2 {input.FASTQ2} --validateMappings -o {output.out_1}"
