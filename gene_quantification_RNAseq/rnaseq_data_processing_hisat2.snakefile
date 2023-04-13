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

# hisat2 sometimes has problem locating the perl5 libraries
#   so we can set this as an environmental variable
#   should be set to the directory containing the lib/ folder in the perl5 installation
perl5lib_path = "PERL5LIB=/home/splaisie/mambaforge/pkgs/perl-5.32.1-2_h7f98852_perl5/" 

rule all:
   input:
        # hisat2 alignment rules
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_XY.sam", sample_name = config["Y_samples"]),
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_XX.sam", sample_name = config["X_samples"]),
        # mark PCR duplicates rules
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XY.bam.bai", sample_name = config["Y_samples"]),
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XX.bam.bai", sample_name = config["X_samples"]),
        # add read groups and alignment stats rules
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_stats_XY.txt",sample_name = config["Y_samples"]),
        expand("mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_stats_XX.txt",sample_name = config["X_samples"]),
        # gene and transcript level quantification rules
        expand("feature_counts_rna_hisat/{sample_name}_HISAT_transcriptCounts_XY.txt", sample_name = config["Y_samples"]),
        expand("feature_counts_rna_hisat/{sample_name}_HISAT_transcriptCounts_XX.txt", sample_name = config["X_samples"]),
        expand("feature_counts_rna_hisat/{sample_name}_HISAT_geneCounts_XY.txt", sample_name = config["Y_samples"]),
        expand("feature_counts_rna_hisat/{sample_name}_HISAT_geneCounts_XX.txt", sample_name = config["X_samples"])

rule HISAT_paired_males:
    input:
        FASTQ1 = "fastq_files_rna/{sample_name}_R1.fastq.gz",
        FASTQ2 = "fastq_files_rna/{sample_name}_R2.fastq.gz"
    output:
        out_1 =  temp("mapped_rna_hisat/"+"{sample_name}_HISAT_pair_XY.sam")
    params:
        HISAT_Index_male = Y_genome,
        perl_path = perl5lib_path
    shell:
        "{params.perl_path} hisat2 -q --rna-strandness RF -p 8 -x {params.HISAT_Index_male} -1 {input.FASTQ1} -2 {input.FASTQ2} -S {output.out_1}"
        
        
rule HISAT_paired_females:
    input:
        FASTQ1 = "fastq_files_rna/{sample_name}_R1.fastq.gz",
        FASTQ2 = "fastq_files_rna/{sample_name}_R2.fastq.gz"
    output:
        out_1 = temp("mapped_rna_hisat/{sample_name}_HISAT_pair_XX.sam")
    params:
        HISAT_Index_female = X_genome,
        perl_path = perl5lib_path
    shell:
        "{params.perl_path} hisat2 -q --rna-strandness RF -p 8 -x {params.HISAT_Index_female} -1 {input.FASTQ1} -2 {input.FASTQ2} -S {output.out_1}"
        
# males
rule samtools_view_males:
    input:
        SAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.sam"
    output:
        BAM = temp("mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam")
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule bam_sort_males:    
    input:
        IN_BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XY.bam"
    output:
        sort_BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_XY.bam"
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule MarkDups_males:
    input:
        sort_BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_XY.bam"
    output:
        BAM = temp("mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_XY.bam"),
        metrics = "stats_hisat/{sample_name}.XY.picard_mkdup_metrics.txt"
    shell:
        "picard -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule AddReadGrps_males:
    input:
        Read_BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_XY.bam"
    output:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XY.bam"
    params:
        id = lambda wildcards: config[wildcards.sample_name]["ID"],
        sm = lambda wildcards: config[wildcards.sample_name]["SM"],
        lb = lambda wildcards: config[wildcards.sample_name]["LB"],
        pu = lambda wildcards: config[wildcards.sample_name]["PU"],
        pl = lambda wildcards: config[wildcards.sample_name]["PL"]
    shell:
        "picard -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule index_bam_males:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XY.bam"
    output:
        BAI = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XY.bam.bai"
    shell:
        "bamtools index -in {input.BAM}"

rule stats_bam_males:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XY.bam"
    output:
        stats = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_stats_XY.txt"
    shell:
        "bamtools stats -in {input.BAM} > {output.stats}"

rule featureCounts_gene_males:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XY.bam"
    output:
        counts = "feature_counts_rna_hisat/{sample_name}_HISAT_geneCounts_XY.txt"
    params:
        GTF = genome_annotation,
    shell:
        "featureCounts -T 8 --primary -p -s 2 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"

rule featureCounts_transcripts_males:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XY.bam"
    output:
        counts = "feature_counts_rna_hisat/{sample_name}_HISAT_transcriptCounts_XY.txt"
    params:
        GTF = genome_annotation,
    shell:
        "featureCounts -T 8 --primary -p -s 2 -t exon -g transcript_id -O -a {params.GTF} -o {output.counts} {input.BAM}"
        
# females
rule samtools_view_females:
    input:
        SAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.sam"
    output:
        BAM = temp("mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam")
    shell:
        "samtools view -b {input.SAM} > {output.BAM}"

rule bam_sort_females:  
    input:
        IN_BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_XX.bam"
    output:
        sort_BAM = temp("mapped_rna_hisat/{sample_name}_HISAT_pair_sort_XX.bam")
    shell:
        "bamtools sort -in {input.IN_BAM} -out {output.sort_BAM}"

rule MarkDups_females:
    input:
        sort_BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_XX.bam"
    output:
        BAM = temp("mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_XX.bam"),
        metrics = "stats_hisat/{sample_name}.XY.picard_mkdup_metrics.txt"
    shell:
        "picard -Xmx14g MarkDuplicates I={input.sort_BAM} O={output.BAM} M={output.metrics} VALIDATION_STRINGENCY=LENIENT"

rule AddReadGrps_females:
    input:
        Read_BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_XX.bam"
    output:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XX.bam",
    params:
        id = lambda wildcards: config[wildcards.sample_name]["ID"],
        sm = lambda wildcards: config[wildcards.sample_name]["SM"],
        lb = lambda wildcards: config[wildcards.sample_name]["LB"],
        pu = lambda wildcards: config[wildcards.sample_name]["PU"],
        pl = lambda wildcards: config[wildcards.sample_name]["PL"]
    shell:
        "picard -Xmx14g AddOrReplaceReadGroups I={input.Read_BAM} O={output.BAM} RGID={params.id} RGPU={params.pu} RGSM={params.sm} RGPL={params.pl} RGLB={params.lb} VALIDATION_STRINGENCY=LENIENT"

rule index_bam_females:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XX.bam"
    output:
        BAI = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XX.bam.bai"
    shell:
        "bamtools index -in {input.BAM}"

rule stats_bam_females:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XX.bam"
    output:
        stats = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_stats_XX.txt"
    shell:
        "bamtools stats -in {input.BAM} > {output.stats}"
 
rule featureCounts_gene_females:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XX.bam"
    output:
        counts = "feature_counts_rna_hisat/{sample_name}_HISAT_geneCounts_XX.txt"
    params:
        GTF = genome_annotation
    shell:
        "featureCounts -T 8 --primary -p -s 2 -t exon -g gene_name -a {params.GTF} -o {output.counts} {input.BAM}"
        
rule featureCounts_transcripts_females:
    input:
        BAM = "mapped_rna_hisat/{sample_name}_HISAT_pair_sort_mkdup_rdgrp_XX.bam"
    output:
        counts = "feature_counts_rna_hisat/{sample_name}_HISAT_transcriptCounts_XX.txt"
    params:
        GTF = genome_annotation
    shell:
        "featureCounts -T 8 --primary -p -s 2 -t exon -g transcript_id -O -a {params.GTF} -o {output.counts} {input.BAM}"

