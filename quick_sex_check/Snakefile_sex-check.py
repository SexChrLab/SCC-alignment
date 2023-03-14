import os

#mamba create -n sex_check snakemake seqtk minimap2 samtools goleft pigz
#mamba env create --name sex_check -f sex_check.yml
#conda activate sex_check
#snakemake --use-conda --cores 4 -s Snakefile_sex-check.py -np

configfile: "configsex-check_AGAVE.json"

samples = config["samples_ALL"]
genome = config["genome"]
indexcov_dir = config["indexcov_dir"]

rule all:
	input:
#subsample_reads rule
		expand("subseq/{sample}_1x_R1.fq.gz", sample=samples),
		expand("subseq/{sample}_1x_R2.fq.gz", sample=samples),
#bwa map_reads rule
		expand("mapped/{sample}_1x.bam", sample=samples),
#index rule
		expand("mapped/{sample}_1x.bam.bai", sample=samples),
#indexcov rule
		"indexcov/index.html",

rule subsample_reads:
	input:
		fq1 = "reads/{sample}_R1.fastq.gz", 
		fq2 = "reads/{sample}_R2.fastq.gz", 
	output:
		fq1 = "subseq/{sample}_1x_R1.fq.gz", 
		fq2 = "subseq/{sample}_1x_R2.fq.gz", 
	params:
		threads = 4,
		seed = 100,
	shell:
		"""
		seqtk sample -s {params.seed} {input.fq1} 10000000 | pigz -c --best -p {params.threads} > {output.fq1};
		seqtk sample -s {params.seed} {input.fq2} 10000000 | pigz -c --best -p {params.threads} > {output.fq2};
		"""

rule map_reads:
	input:
		fq1 = "subseq/{sample}_1x_R1.fq.gz", 
		fq2 = "subseq/{sample}_1x_R2.fq.gz", 
	output:
		"mapped/{sample}_1x.bam",
	params:
		threads = 4,
		genome = genome,
	shell:
		"""
		minimap2 -ax sr {params.genome} {input.fq1} {input.fq2} -t {params.threads} | samtools sort -@ 4 -O bam - -o {output};
		"""

rule index_bam:
	input:
		"mapped/{sample}_1x.bam",
	output:
		"mapped/{sample}_1x.bam.bai",
	shell:
		"""
		samtools index {input};
		"""

rule indexcov:
	input:
		bam = expand("mapped/{sample}_1x.bam", sample=samples),
		bai = expand("mapped/{sample}_1x.bam.bai", sample=samples),
	params:
		indexcov_dir = indexcov_dir,
	output:
		"indexcov/index.html",
	shell:
		"""
		goleft indexcov --directory {params.indexcov_dir} {input.bam};
		"""
