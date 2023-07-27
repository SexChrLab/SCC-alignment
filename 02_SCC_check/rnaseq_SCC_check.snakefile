import os

configfile: "RNAseq_samples.json"

# path to salmon indexes for sex chromosome complement genome

# uncomment these to use HG38 (GRCh38) SCC reference genomes
#X_genome = config["HG38_Transcriptome_Index_SALMON_Path_female"]
Y_genome = "/data/CEM/wilsonlab/projects/XYAlign_Tutorial/PC_chrY_transcript_salmon_index" 

# uncomment these to use CHM13 (telomere-to-telomere) SCC reference genomes
#X_genome = config["CHM13_Transcriptome_Index_SALMON_Path_female"]
#Y_genome = config["CHM13_Transcriptome_Index_SALMON_Path_male"]

rule all:
   input:
        # use salmon to quickly look for chrY matches with the Y PARs masked genome
        expand("SCC_check_rna/{sample_name}_salmon_quant/", sample_name = config["ALL_samples"]),

rule salmon_paired_males:
    input:
        FASTQ1 = "fastq_files_rna/{sample_name}_R1.fastq.gz",
        FASTQ2 = "fastq_files_rna/{sample_name}_R2.fastq.gz"
    output:
        out_1 =  directory("SCC_check_rna/"+"{sample_name}_salmon_quant/")
    params:
        SALMON_Index_male = Y_genome,
        libtype = "ISR", # LIBTYPE for paired reads
    shell:
        "salmon quant -i {params.SALMON_Index_male} -l {params.libtype} -1 {input.FASTQ1} -2 {input.FASTQ2} --validateMappings -o {output.out_1}"
