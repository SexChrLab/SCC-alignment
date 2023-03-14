import os

#mamba create -n sex_check snakemake seqtk minimap2 samtools samblaster goleft pigz gatk4=4.2.6.1-0 picard=2.27.0-0
#mamba env create --name sex_check -f sex_check.yml
#mamba activate w
#snakemake --use-conda --cores 4 -s Snakefile_SNPs.py -np

#shorthand:
#X == no Y chromosome
#Y == Y present

configfile: "config_SNPs.json"

#sample sets
X_samples = config["X_samples"]
Y_samples = config["Y_samples"]
ALL_samples = config["ALL_samples"]
dir = config["directory"]

#genome sets
#X_genome = config["GRCh38_X"]
#Y_genome = config["GRCh38_Y"]
X_genome = config["CHM13_X"]
Y_genome = config["CHM13_Y"]

#chromosomes/intervals
all_chromosomes = config["all_chromosomes"]
diploid = config["diploid_chromosomes"]
autosomes = config["autosomes"]
#PAR1 = config["GRCh38_PAR1"]
#PAR2 = config["GRCh38_PAR2"]
#nonPAR = config["GRCh38_nonPAR"]
PAR1 = config["CHM13_PAR1"]
PAR2 = config["CHM13_PAR2"]
nonPAR = config["CHM13_nonPAR"]

rule all:
	input:
#minimap2_map_reads rule
#		expand("temp/X/{X}.bam", X=X_samples),
#		expand("temp/Y/{Y}.bam", Y=Y_samples),
#index_stat rules
		expand("mapped/X/{X}.bai", X=X_samples),
		expand("stats_X/{X}.bam.stats", X=X_samples),
		expand("mapped/Y/{Y}.bai", Y=Y_samples),
		expand("stats/Y/{Y}.bam.stats", Y=Y_samples),
#gatk_gvcfs rule
		expand("haplotyped_vcfs/X/{X}.{chr_n}.g.vcf.gz", X=X_samples, chr_n=diploid),
		expand("haplotyped_vcfs/Y/{Y}.{chr_n}.g.vcf.gz", Y=Y_samples, chr_n=all_chromosomes),
#gatk_combinegvcfs_auto rule
		expand("haplotyped_vcfs/temp/{X,Y}.{chr_n}.g.vcf.gz", X=X_samples, Y=Y_samples, chr_n=diploid),
#		expand("haplotyped_vcfs/temp/{Y}.{chr_n}.g.vcf.gz", Y=Y_samples, chr_n=all_chromosomes),
#		expand("haplotyped_vcfs/merged_vcfs/{chr_n}.gatk.combined.g.vcf.gz", chr_n=autosomes),
#gatk_genotypegvcf_auto rule
#		expand("genotyped_vcfs/{chr_n}.gatk.raw.vcf.gz", chr_n=autosomes),
#gatk_genotypegvcf_chrX rules
#		expand("genotyped_vcfs/{chr_n}.gatk.raw.vcf.gz", chr_n=autosomes),
#gatk_genotypegvcf_chrY rules

#hard_filter_auto rule
#		expand("genotyped_vcfs/{chr_n}.gatk.filtered.vcf.gz", chr_n=autosomes),
#		expand("genotyped_vcfs/{chr_n}.gatk.filtered.vcf.gz.tbi", chr_n=autosomes),
#subset_individuals_hard_filter_diploid rule
#		expand("genotyped_vcfs/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi", sample=X_samples, chr_n=diploid),
#gatk_selectheterozygous_diploid rule
#		expand("genotyped_vcfs/{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz", sample=X_samples, chr_n=diploid),

# ------------------
# Read mapping rules
# ------------------
rule minimap2_mapping_X:
	input:
		fq1 = "reads/{X}_R1.fastq.gz", 
		fq2 = "reads/{X}_R2.fastq.gz", 
	output:
		temp("temp/X/{X}.bam"),
	params:
		threads = 4,
		X_genome = X_genome,
	shell:
		"""
		minimap2 -ax sr {params.X_genome} {input.fq1} {input.fq2} -t {params.threads} | samblaster | samtools fixmate -@ {params.threads} - - | samtools sort -@ {params.threads} -O bam - -o {output} 2>/dev/null 
		"""

rule minimap2_mapping_Y:
	input:
		fq1 = "reads/{Y}_R1.fastq.gz", 
		fq2 = "reads/{Y}_R2.fastq.gz", 
	output:
		temp("temp/Y/{Y}.bam"),
	params:
		threads = 4,
		Y_genome = Y_genome,
	shell:
		"""
		minimap2 -ax sr {params.Y_genome} {input.fq1} {input.fq2} -t {params.threads} | samblaster | samtools fixmate -@ {params.threads} - - | samtools sort -@ {params.threads} -O bam - -o {output} 2>/dev/null 
		"""

# ------------------
# Add read groups to bam files
# ------------------
rule index_readgroups_X:
	input:
		"temp/X/{X}.bam",
	output:
		bam = "mapped/X/{X}.bam",
		bai = "mapped/X/{X}.bai",
	params:
		threads = 2,
		id = lambda wildcards: config[wildcards.X]["ID"],
		sm = lambda wildcards: config[wildcards.X]["SM"],
		lb = lambda wildcards: config[wildcards.X]["LB"],
		pu = lambda wildcards: config[wildcards.X]["PU"],
		pl = lambda wildcards: config[wildcards.X]["PL"],
	shell:
		"""
		picard AddOrReplaceReadGroups -I {input} -O {output.bam} -SORT_ORDER coordinate -RGID {params.id} -RGLB {params.lb} -RGPL {params.pl} -RGPU {params.pu} -RGSM {params.sm} -CREATE_INDEX True
		"""

rule index_readgroups_Y:
	input:
		"temp/Y/{Y}.bam",
	output:
		bam = "mapped/Y/{Y}.bam",
		bai = "mapped/Y/{Y}.bai",
	params:
		threads = 2,
		id = lambda wildcards: config[wildcards.Y]["ID"],
		sm = lambda wildcards: config[wildcards.Y]["SM"],
		lb = lambda wildcards: config[wildcards.Y]["LB"],
		pu = lambda wildcards: config[wildcards.Y]["PU"],
		pl = lambda wildcards: config[wildcards.Y]["PL"],
	shell:
		"""
		picard AddOrReplaceReadGroups -I {input} -O {output.bam} -SORT_ORDER coordinate -RGID {params.id} -RGLB {params.lb} -RGPL {params.pl} -RGPU {params.pu} -RGSM {params.sm} -CREATE_INDEX True
		"""

# ------------------
# Calculate read mapping stats for QC
# ------------------
rule stats_X:
	input:
		"mapped/X/{sample}.bam",
	output:
		"stats_X/{sample}.bam.stats",
	shell:
		"""
		samtools stats {input} | grep '^SN' | cut -f 2- > {output};
		"""

rule stats_Y:
	input:
		"mapped/Y/{sample}.bam",
	output:
		"stats/Y/{sample}.bam.stats",
	shell:
		"""
		samtools stats {input} | grep '^SN' | cut -f 2- > {output};
		"""

# ---------------------------------------
# Call haplotypes across all samples/chromosomes
# ---------------------------------------

rule gatk_gvcfs_X:
	input:
		bam = "mapped/X/{X}.bam",
		bai = "mapped/X/{X}.bai",
	output:
		"haplotyped_vcfs/X/{X}.{chr_n}.g.vcf.gz",
	params:
		X_genome = X_genome,
		chr_n = "{chr_n}",
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.X_genome} -I {input.bam} -L {params.chr_n} -ERC GVCF --output {output};
		"""

rule gatk_gvcfs_Y:
	input:
		bam = "mapped/Y/{Y}.bam",
		bai = "mapped/Y/{Y}.bai",
	output:
		"haplotyped_vcfs/Y/{Y}.{chr_n}.g.vcf.gz",
	params:
		Y_genome = Y_genome,
		chr_n = "{chr_n}",
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.Y_genome} -I {input.bam} -L {params.chr_n} -ERC GVCF --output {output};
		"""

## ---------------------------------------
## Compile all gvcfs
## ---------------------------------------
rule merge_XY_dirs:
	input:
		"haplotyped_vcfs/{dir}/{X,Y}.{chr_n}.g.vcf.gz",
	output:
		"haplotyped_vcfs/temp/{X,Y}.{chr_n}.g.vcf.gz",
	params:
		chr_n = "{chr_n}",
	threads:
		4
	shell:
		"""
		cp haplotyped_vcfs/X/{X,Y}.{chr_n}.g.vcf.gz haplotyped_vcfs/temp/{X,Y}.{chr_n}.g.vcf.gz;
		"""

## ---------------------------------------
## Processing of autosomes (chr1 to chr22)
## ---------------------------------------

#rule gatk_combinegvcfs_auto:
#	input:
#		gvcfs = expand("haplotyped_vcfs/{sample}.{{chr_n}}.g.vcf.gz", zip, sample=ALL_samples, chr_n=config["autosomes"]),
#	output:
#		"haplotyped_vcfs/merged_autosomes/{chr_n}.gatk.combined.g.vcf.gz",
#	params:
#		Y_genome = Y_genome,
#		chr_n = "{chr_n}",
#	threads:
#		1
#	run:
#		variant_files = []
#		for i in input.gvcfs:
#			variant_files.append("--variant " + i)
#		variant_files = " ".join(variant_files)
#		shell(
#			"""
#			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.Y_genome} {variant_files} --intervals {params.chr_n} -O {output}
#			""")

#rule gatk_genotypegvcf_auto:
#	input:
#		expand("haplotyped_vcfs/merged_vcfs/{chr_n}.gatk.combined.g.vcf.gz", chr_n=config["autosomes"]),
#	output:
#		"genotyped_vcfs/{chr_n}.gatk.raw.vcf.gz",
#	params:
#		genome = genome,
#	threads:
#		1
#	shell:
#		"""
#		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.genome} -V {input} -O {output}
#		"""

## ---------------------------------------
## Processing of chrX
## ---------------------------------------
#
#
#
#
#
#
#
#
#
## ---------------------------------------
## Processing of chrY
## ---------------------------------------
#
#rule gatk_combinegvcfs_chrY:
#	input:
#		gvcfs = expand("haplotyped_vcfs/{sample}.chrY.g.vcf.gz", zip, sample=Y_samples),
#	output:
#		"haplotyped_vcfs/merged_vcfs/chrY.gatk.combined.g.vcf.gz",
#	params:
#		genome = genome,
#		chr_n = "{chr_n}",
#	threads:
#		1
#	run:
#		variant_files = []
#		for i in input.gvcfs:
#			variant_files.append("--variant " + i)
#		variant_files = " ".join(variant_files)
#		shell(
#			"""
#			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.genome} {variant_files} --intervals {params.chr_n} -O {output}
#			""")
#
#rule gatk_genotypegvcf_chrY:
#	input:
#		expand("haplotyped_vcfs/merged_vcfs/{{chr_n}}.gatk.combinegvcf.g.vcf.gz", chr_n=config["autosomes"]),
#	output:
#		"genotyped_vcfs/{chr_n}.gatk.raw.vcf.gz",
#	params:
#		genome = genome,
#	threads:
#		1
#	shell:
#		"""
#		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.genome} -V {input} -O {output}
#		"""
#
#
#
#
## ------------------------
## Hard Filter on autosomes
## ------------------------
#
#rule hard_filter_diploid:
#	input:
#		"genotyped_vcfs/{chr_n}.gatk.called.raw.vcf.gz"
#	output:
#		vcf = "genotyped_vcfs/{chr_n}.gatk.called.hard.filter.vcf.gz",
#		idx = "genotyped_vcfs/{chr_n}.gatk.called.hard.filter.vcf.gz.tbi",
#	params:
#		genome = genome,
#	shell:
#		"""
#		gatk SelectVariants -R {params.genome} -V {input} -O {output} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -select "AN >= 4 && MQ > 40.0 && QD > 7.0 && DP >= 10.0 && DP <= 1000.0";
#		touch -c {output.idx};
#		"""
#
#rule subset_individuals_hard_filter_diploid:
#	input:
#		"genotyped_vcfs/{chr_n}.gatk.called.hard.filter.vcf.gz",
#	output:
#		vcf = "genotyped_vcfs/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz",
#		index = "genotyped_vcfs/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi",
#	params:
#		sample = lambda wildcards: config[wildcards.sample]["SM"],
#	shell:
#		"""
#		bcftools view -Oz -s {params.sample} {input} > {output.vcf};
#		tabix -p vcf {output.vcf};
#		"""

#After subsetting for each individual. In some individuals, the genotypes could be homozygous for the reference. This next rule is to remove these sites.
#
#rule gatk_selectheterozygous_diploid:
#	input:
#		"genotyped_vcfs/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz", 
#	output:
#		"diploid_chromosomes_vcfs/{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz",
#	params:
#		genome = genome,
#	shell:
#		"""
#		gatk SelectVariants -R {params.genome} -V {input} -O {output} -select "AC == 1";
#		"""
#