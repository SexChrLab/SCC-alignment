import os

configfile: "RNAseq_samples.json"

# path to HISAT indexes for sex chromosome complement genome

# uncomment these if you want to run with the HG38 (GRCh38) human genome
#X_genome = config["HG38_Transcriptome_Index_HISAT_Path_female"]
#Y_genome = config["HG38_Transcriptome_Index_HISAT_Path_male"]
#genome_annotation = config["HG38_annotation_path"]

# uncomment these if you want to run with the CHM13 (telomere-to-telomere) human genome
X_genome = config["CHM13_Transcriptome_Index_HISAT_Path_female"]
Y_genome = config["CHM13_Transcriptome_Index_HISAT_Path_male"]
genome_annotation = config["CHM13_annotation_path"]

rule all:
   input:
        # hisat2 alignment rules
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam", sample_name = config["Y_samples"]),
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam", sample_name = config["X_samples"]),
        # index rules
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam.bai", sample_name = config["Y_samples"]),
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam.bai", sample_name = config["X_samples"]),
        # stats rules
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_XY.txt", sample_name = config["Y_samples"]),
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_XX.txt", sample_name = config["X_samples"]),
        # gene and transcript level quantification rules
        expand("feature_counts_rna_hisat/{sample_name}_HISAT_transcriptCounts_XY.txt", sample_name = config["Y_samples"]),
        expand("feature_counts_rna_hisat/{sample_name}_HISAT_transcriptCounts_XX.txt", sample_name = config["X_samples"]),
        expand("feature_counts_rna_hisat/{sample_name}_HISAT_geneCounts_XY.txt", sample_name = config["Y_samples"]),
        expand("feature_counts_rna_hisat/{sample_name}_HISAT_geneCounts_XX.txt", sample_name = config["X_samples"])

rule link_fastqs:
    input:
        original_R1 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq_1"],
        original_R2 = lambda wildcards: config[wildcards.sample_name]["fq_path"] + config[wildcards.sample_name]["fq_2"]
    output:
        R1_out = "fastq_files_rna/{sample_name}_R1.fastq.gz",
        R2_out = "fastq_files_rna/{sample_name}_R2.fastq.gz"
    shell:
        """
        ln -s {input.original_R1} {output.R1_out};
        ln -s {input.original_R2} {output.R2_out}
        """

rule HISAT_paired_males:
    input:
        FASTQ1 = "fastq_files_rna/{sample_name}_R1.fastq.gz",
        FASTQ2 = "fastq_files_rna/{sample_name}_R2.fastq.gz"
    output:
        "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam"
    params:
        HISAT_Index_male = Y_genome,
        threads = 4,
        id = lambda wildcards: config[wildcards.sample_name]["ID"],
        sm = lambda wildcards: config[wildcards.sample_name]["SM"],
        lb = lambda wildcards: config[wildcards.sample_name]["LB"],
        pu = lambda wildcards: config[wildcards.sample_name]["PU"],
        pl = lambda wildcards: config[wildcards.sample_name]["PL"]
    shell:
        """
        hisat2 --rna-strandness RF -p {params.threads} --dta --rg-id {params.id} --rg SM:{params.sm} --rg LB:{params.lb} --rg PU:{params.pu} --rg PL:{params.pl} -x {params.HISAT_Index_male} -1 {input.FASTQ1} -2 {input.FASTQ2} | samblaster | samtools sort -@ {params.threads} -O bam - -o {output};
        """
       
        
rule HISAT_paired_females:
    input:
        FASTQ1 = "fastq_files_rna/{sample_name}_R1.fastq.gz",
        FASTQ2 = "fastq_files_rna/{sample_name}_R2.fastq.gz"
    output:
        "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam"
    params:
        HISAT_Index_female = X_genome,
        threads = 4,
        id = lambda wildcards: config[wildcards.sample_name]["ID"],
        sm = lambda wildcards: config[wildcards.sample_name]["SM"],
        lb = lambda wildcards: config[wildcards.sample_name]["LB"],
        pu = lambda wildcards: config[wildcards.sample_name]["PU"],
        pl = lambda wildcards: config[wildcards.sample_name]["PL"]
    shell:
        """
        hisat2 --rna-strandness RF -p {params.threads} --dta --rg-id {params.id} --rg SM:{params.sm} --rg LB:{params.lb} --rg PU:{params.pu} --rg PL:{params.pl} -x {params.HISAT_Index_female} -1 {input.FASTQ1} -2 {input.FASTQ2} | samblaster | samtools sort -@ {params.threads} -O bam - -o {output};
        """
        
# males
rule index_bam_males:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam"
    output:
        BAI = "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam.bai"
    shell:
        "bamtools index -in {input.BAM}"

rule stats_bam_males:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam"
    output:
        stats = "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.txt"
    shell:
        "bamtools stats -in {input.BAM} > {output.stats}"

rule featureCounts_gene_males:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam"
    output:
        counts = "feature_counts_rna_hisat/{sample_name}_HISAT_geneCounts_XY.txt"
    params:
        GTF = genome_annotation,
        threads = 4
    shell:
        "featureCounts -T {params.threads} --primary -p -s 2 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

rule featureCounts_transcripts_males:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam"
    output:
        counts = "feature_counts_rna_hisat/{sample_name}_HISAT_transcriptCounts_XY.txt"
    params:
        GTF = genome_annotation,
        threads = 4
    shell:
        "featureCounts -T {params.threads} --primary -p -s 2 -t exon -g transcript_id -O -a {params.GTF} -o {output.counts} {input.BAM}"
        
# females
rule index_bam_females:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam"
    output:
        BAI = "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam.bai"
    shell:
        "bamtools index -in {input.BAM}"

rule stats_bam_females:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam"
    output:
        stats = "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.txt"
    shell:
        "bamtools stats -in {input.BAM} > {output.stats}"
 
rule featureCounts_gene_females:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam"
    output:
        counts = "feature_counts_rna_hisat/{sample_name}_HISAT_geneCounts_XX.txt"
    params:
        GTF = genome_annotation,
        threads = 4
    shell:
        "featureCounts -T {params.threads} --primary -p -s 2 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"
        
rule featureCounts_transcripts_females:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam"
    output:
        counts = "feature_counts_rna_hisat/{sample_name}_HISAT_transcriptCounts_XX.txt"
    params:
        GTF = genome_annotation,
        threads = 4
    shell:
        "featureCounts -T {params.threads} --primary -p -s 2 -t exon -g transcript_id -O -a {params.GTF} -o {output.counts} {input.BAM}"

