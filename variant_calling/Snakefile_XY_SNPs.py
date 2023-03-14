import os

#mamba create -n sex_check snakemake seqtk minimap2 samtools samblaster goleft pigz gatk4=4.2.6.1-0 picard=2.27.0-0
#mamba env create --name sex_check -f sex_check.yml
#mamba activate w
#snakemake --use-conda --cores 4 -s Snakefile_SNPs.py -np

#shorthand:
#X == no chrY present
#Y == chrY present

configfile: "config_SNPs.json"

#sample sets
X_samples = config["X_samples"]
Y_samples = config["Y_samples"]
ALL_samples = config["ALL_samples"]
#dir = config["directory"]

#chromosomes/intervals
all_chromosomes = config["all_chromosomes"]
diploid = config["XX_diploid"]
auto = config["autosomes"]
all_regions = config["all_regions"]

#Reference genome choices, only ONE may be used at a time
#GRCh38 reference genome (uncomment below if using)
#X_genome = config["GRCh38_X"]
#Y_genome = config["GRCh38_Y"]
#PAR1 = config["GRCh38_PAR1"]
#PAR2 = config["GRCh38_PAR2"]
#nonPAR = config["GRCh38_nonPAR"]

#CHM13_T2T reference (comment out below if NOT using)
X_genome = config["CHM13_X"]
Y_genome = config["CHM13_Y"]
PAR1 = config["CHM13_PAR1"]
PAR2 = config["CHM13_PAR2"]
nonPAR = config["CHM13_nonPAR"]

rule all:
	input:
########################Stage 1: Map reads and QC for all Y+ samples
#minimap2_map_reads rule
		expand("temp/Y/{Y}.bam", Y=Y_samples),
#index_stat rules
		expand("mapped/Y/{Y}.bai", Y=Y_samples),
		expand("stats/Y/{Y}.bam.stats", Y=Y_samples),

########################Stage 2: Call haplotypes for all genomic regions (ploidy-aware)
#HaplotypeCaller
#gatk_gvcfs_auto rule
		expand("haplotyped_vcfs/Y/{Y}.{chr}.g.vcf.gz", Y=Y_samples, chr=auto),
#gatk_gvcfs_PAR1 rule
		expand("haplotyped_vcfs/Y/chrX/{Y}.chrX_PAR1.g.vcf.gz", Y=Y_samples),
#gatk_gvcfs_PAR2 rule
		expand("haplotyped_vcfs/Y/chrX/{Y}.chrX_PAR2.g.vcf.gz", Y=Y_samples),
#gatk_gvcfs_chrX rule
		expand("haplotyped_vcfs/Y/chrX/{Y}.chrX_nonPAR.g.vcf.gz", Y=Y_samples),
#gatk_gvcfs_chrY rule
		expand("haplotyped_vcfs/Y/chrY_nonPAR/{Y}.chrY_nonPAR.g.vcf.gz", Y=Y_samples),

########################Stage 3: Merge gVCFs for all genomic regions
#CombineGVCFs
#gatk_combinegvcfs_auto rule
		expand("haplotyped_vcfs/Y/merged/{chr}.gatk.combined.gvcf.vcf.gz", chr=auto),
#gatk_combinegvcfs_PAR1
		expand("haplotyped_vcfs/Y/merged/chrX_PAR1.gatk.combined.gvcf.vcf.gz"),
#gatk_combinegvcfs_PAR2
		expand("haplotyped_vcfs/Y/merged/chrX_PAR2.gatk.combined.gvcf.vcf.gz"),
#gatk_combinegvcfs_chrX
		expand("haplotyped_vcfs/Y/merged/chrX_nonPAR.gatk.combined.gvcf.vcf.gz"),
#gatk_combinegvcfs_chrY
		expand("haplotyped_vcfs/Y/merged/chrY_nonPAR.gatk.combined.gvcf.vcf.gz"),

########################Stage 4: Compute joint genotypes across all samples (ploidy-aware)
#GenotypeGVCFs
#gatk_genotype_auto rule
		expand("genotyped_vcfs/Y/{chr_n}.gatk.genotyped.raw.vcf.gz", chr_n=auto),
#gatk_genotype_PAR1 rule
		"genotyped_vcfs/Y/chrX_PAR1.gatk.genotyped.raw.vcf.gz",
#gatk_genotype_PAR2 rule
		"genotyped_vcfs/Y/chrX_PAR2.gatk.genotyped.raw.vcf.gz",
#gatk_genotype_X rule
		"genotyped_vcfs/Y/chrX_nonPAR.gatk.genotyped.raw.vcf.gz",
#gatk_genotype_Y rule
		"genotyped_vcfs/Y/chrY.gatk.genotyped.raw.vcf.gz",

########################Stage 5: Merge all three chrX region vcf files
#MergeVcfs
#gatk_genotypegvcf_mergeX rule
		"genotyped_vcfs/Y/chrX.gatk.genotyped.raw.vcf.gz",

########################Stage 6: Downstream analysis-aware VCF filtering
#hard_filter_auto rule
#		expand("genotyped_vcfs/{chr_n}.gatk.filtered.vcf.gz", chr_n=diploid),
#subset_individuals_hard_filter_diploid rule
#		expand("genotyped_vcfs/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi", sample=Y_samples, chr_n=diploid),
#gatk_selectheterozygous_diploid rule
#		expand("diploid_vcfs/{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz", sample=Y_samples, chr_n=diploid),

# ------------------
# Read mapping rules
# ------------------
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

rule gatk_gvcfs_diploid_auto:
	input:
		bam = "mapped/Y/{Y}.bam",
		bai = "mapped/Y/{Y}.bai",
	output:
		"haplotyped_vcfs/Y/{Y}.{chr}.g.vcf.gz",
	params:
		Y_genome = Y_genome,
		ploidy = 2,
		chr = "{chr}",
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.Y_genome} -I {input.bam} -L {params.chr} -ERC GVCF --output {output};
		"""

rule gatk_gvcfs_diploid_PAR1:
	input:
		bam = "mapped/Y/{Y}.bam",
		bai = "mapped/Y/{Y}.bai",
	output:
		"haplotyped_vcfs/Y/chrX/{Y}.chrX_PAR1.g.vcf.gz",
	params:
		Y_genome = Y_genome,
		ploidy = 2,
		PAR1 = PAR1,
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.Y_genome} -I {input.bam} -L {params.PAR1} -ploidy {params.ploidy} -ERC GVCF --output {output};
		"""

rule gatk_gvcfs_diploid_PAR2:
	input:
		bam = "mapped/Y/{Y}.bam",
		bai = "mapped/Y/{Y}.bai",
	output:
		"haplotyped_vcfs/Y/chrX/{Y}.chrX_PAR2.g.vcf.gz",
	params:
		Y_genome = Y_genome,
		ploidy = 2,
		PAR2 = PAR2,
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.Y_genome} -I {input.bam} -L {params.PAR2} -ploidy {params.ploidy} -ERC GVCF --output {output};
		"""

rule gatk_gvcfs_haploid_chrX:
	input:
		bam = "mapped/Y/{Y}.bam",
		bai = "mapped/Y/{Y}.bai",
	output:
		"haplotyped_vcfs/Y/chrX/{Y}.chrX_nonPAR.g.vcf.gz",
	params:
		Y_genome = Y_genome,
		ploidy = 1,
		chrX = nonPAR,
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.Y_genome} -I {input.bam} -L {params.chrX} -ploidy {params.ploidy} -ERC GVCF --output {output};
		"""

rule gatk_gvcfs_haploid_chrY:
	input:
		bam = "mapped/Y/{Y}.bam",
		bai = "mapped/Y/{Y}.bai",
	output:
		"haplotyped_vcfs/Y/chrY_nonPAR/{Y}.chrY_nonPAR.g.vcf.gz",
	params:
		Y_genome = Y_genome,
		ploidy = 1,
		chrY = "chrY",
	threads:
		4
	shell:
		"""
		gatk HaplotypeCaller -R {params.Y_genome} -I {input.bam} -L {params.chrY} -ploidy {params.ploidy} -ERC GVCF --output {output};
		"""

## ---------------------------------------
## Combine all gvcfs
## ---------------------------------------

rule gatk_combinegvcfs_auto:
	input:
		expand("haplotyped_vcfs/Y/{Y}.{chr}.g.vcf.gz", Y=Y_samples, chr=auto),
	output:
		"haplotyped_vcfs/Y/merged/{chr}.gatk.combined.gvcf.vcf.gz",
	params:
		Y_genome = Y_genome,
		chr = "{chr}",
	threads:
		1
	run:
		variant_files = []
		for i in input:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.Y_genome} {variant_files} --intervals {params.chr} -O {output}
			""")

rule gatk_combinegvcfs_PAR1:
	input:
		expand("haplotyped_vcfs/Y/chrX/{Y}.chrX_PAR1.g.vcf.gz", Y=Y_samples),
	output:
		"haplotyped_vcfs/Y/merged/chrX_PAR1.gatk.combined.gvcf.vcf.gz",
	params:
		Y_genome = Y_genome,
		PAR1 = PAR1,
	threads:
		1
	shell:
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.Y_genome} -V {input} -O {output}
			"""

rule gatk_combinegvcfs_PAR2:
	input:
		expand("haplotyped_vcfs/Y/chrX/{Y}.chrX_PAR2.g.vcf.gz", Y=Y_samples),
	output:
		"haplotyped_vcfs/Y/merged/chrX_PAR2.gatk.combined.gvcf.vcf.gz",
	params:
		Y_genome = Y_genome,
		PAR2 = PAR2,
	threads:
		1
	shell:
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.Y_genome} -V {input} -O {output}
			"""

rule gatk_combinegvcfs_chrX:
	input:
		expand("haplotyped_vcfs/Y/chrX/{Y}.chrX_nonPAR.g.vcf.gz", Y=Y_samples),
	output:
		"haplotyped_vcfs/Y/merged/chrX_nonPAR.gatk.combined.gvcf.vcf.gz",
	params:
		Y_genome = Y_genome,
		nonPAR = nonPAR,
	threads:
		1
	shell:
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.Y_genome} -V {input} -O {output}
			"""

rule gatk_combinegvcfs_chrY:
	input:
		expand("haplotyped_vcfs/Y/chrY_nonPAR/{Y}.chrY_nonPAR.g.vcf.gz", Y=Y_samples),
	output:
		"haplotyped_vcfs/Y/merged/chrY_nonPAR.gatk.combined.gvcf.vcf.gz",
	params:
		Y_genome = Y_genome,
		chrY = "chrY",
	threads:
		1
	shell:
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.Y_genome} -V {input} -O {output}
			"""

## ---------------------------------------
## Call genotypes autosomes
## ---------------------------------------

rule gatk_genotypegvcf_auto:
	input:
		"haplotyped_vcfs/Y/merged/{chr_n}.gatk.combined.gvcf.vcf.gz",
	output:
		"genotyped_vcfs/Y/{chr_n}.gatk.genotyped.raw.vcf.gz",
	params:
		Y_genome = Y_genome,
		chr_n = "{chr_n}",
	threads:
		1
	shell:
		"""
		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.Y_genome} -L {params.chr_n} -V {input} -O {output}
		"""

## ---------------------------------------
## Call genotypes sex chromosomes, haploid X and Y, diploid PARs --> merge genotyped X vcfs
## ---------------------------------------

##gatk_genotype_PAR1 rule
rule gatk_genotypegvcf_PAR1:
	input:
		"haplotyped_vcfs/Y/merged/chrX_PAR1.gatk.combined.gvcf.vcf.gz",
	output:
		temp("genotyped_vcfs/Y/chrX_PAR1.gatk.genotyped.raw.vcf.gz"),
	params:
		Y_genome = Y_genome,
		ploidy = 2,
		PAR1 = PAR1,
	threads:
		1
	shell:
		"""
		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.Y_genome} -L {params.PAR1} -ploidy {params.ploidy} -V {input} -O {output}
		"""

##gatk_genotype_PAR2 rule
rule gatk_genotypegvcf_PAR2:
	input:
		"haplotyped_vcfs/Y/merged/chrX_PAR2.gatk.combined.gvcf.vcf.gz",
	output:
		temp("genotyped_vcfs/Y/chrX_PAR2.gatk.genotyped.raw.vcf.gz"),
	params:
		Y_genome = Y_genome,
		ploidy = 2,
		PAR2 = PAR2,
	threads:
		1
	shell:
		"""
		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.Y_genome} -L {params.PAR2} -ploidy {params.ploidy} -V {input} -O {output}
		"""

##gatk_genotype_X rule
rule gatk_genotypegvcf_XnonPAR:
	input:
		"haplotyped_vcfs/Y/merged/chrX_nonPAR.gatk.combined.gvcf.vcf.gz",
	output:
		temp("genotyped_vcfs/Y/chrX_nonPAR.gatk.genotyped.raw.vcf.gz"),
	params:
		Y_genome = Y_genome,
		nonPAR = nonPAR,
		ploidy = 1,
	threads:
		1
	shell:
		"""
		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.Y_genome} -L {params.nonPAR} -ploidy {params.ploidy} -V {input} -O {output}
		"""

#gatk_genotype_Y rule
rule gatk_genotypegvcf_YnonPAR:
	input:
		"haplotyped_vcfs/Y/merged/chrY_nonPAR.gatk.combined.gvcf.vcf.gz",
	output:
		"genotyped_vcfs/Y/chrY.gatk.genotyped.raw.vcf.gz",
	params:
		Y_genome = Y_genome,
		ploidy = 1,
		chrY = "chrY",
	threads:
		1
	shell:
		"""
		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.Y_genome} -L {params.chrY} -ploidy {params.ploidy} -V {input} -O {output}
		"""

#gatk_genotype_merge_X rule
rule gatk_genotypegvcf_mergeX:
	input:
		PAR1_vcf = "genotyped_vcfs/Y/chrX_PAR1.gatk.genotyped.raw.vcf.gz",
		nonPAR_vcf = "genotyped_vcfs/Y/chrX_nonPAR.gatk.genotyped.raw.vcf.gz",
		PAR2_vcf = "genotyped_vcfs/Y/chrX_PAR2.gatk.genotyped.raw.vcf.gz",
	output:
		"genotyped_vcfs/Y/chrX.gatk.genotyped.raw.vcf.gz",
	threads:
		2
	shell:
		"""
		picard MergeVcfs \
          I={input.PAR1_vcf} \
          I={input.nonPAR_vcf} \
		  I={input.PAR2_vcf} \
          O={output}
		"""

## ------------------------
## Hard Filter on diploid_chromosomes
## ------------------------

#rule hard_filter_diploid:
#	input:
#		"genotyped_vcfs/{chr_n}.gatk.genotyped.raw.vcf.gz",
#	output:
#		vcf = "genotyped_vcfs/{chr_n}.gatk.filtered.vcf.gz",
#		idx = "genotyped_vcfs/{chr_n}.gatk.filtered.vcf.gz.tbi",
#	params:
#		Y_genome = Y_genome,
#	shell:
#		"""
#		gatk SelectVariants -R {params.Y_genome} -V {input} -O {output} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -select "AN >= 4 && MQ > 40.0 && QD > 7.0 && DP >= 10.0 && DP <= 1000.0";
#		touch -c {output.idx};
#		"""
#
#rule subset_individuals_hard_filter_diploid:
#	input:
#		"genotyped_vcfs/{chr_n}.gatk.filtered.vcf.gz",
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
#
##After subsetting for each individual. In some individuals, the genotypes could be homozygous for the reference. This next rule is to remove these sites.
#
#rule gatk_selectheterozygous_diploid:
#	input:
#		"genotyped_vcfs/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz", 
#	output:
#		"diploid_vcfs/{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz",
#	params:
#		Y_genome = Y_genome,
#	shell:
#		"""
#		gatk SelectVariants -R {params.Y_genome} -V {input} -O {output} -select "AC == 1";
#		"""


#
#
#
#
#
#
#
#
#
#import os
#
#configfile: "config_SNPs.json"
#
#XY_samples = config["XY_samples"]
#genome = config["genome_XY"]
#all_chromosomes = config["all_chromosomes"]
#samples_ALL = config["samples_ALL"]
##autosomes = config["autosomes"]
#
#rule all:
#	input:
##minimap2_map_reads rule
#		expand("temp/{sample}.bam", sample=samples_ALL),
##index_stat rules
#		expand("mapped_reads/{sample}.bam", sample=XY_samples),
#		expand("mapped_reads/{sample}.bai", sample=XY_samples),
#		expand("mapped_reads/{sample}.bam.stats", sample=XY_samples),
##gatk_gvcfs rule
#		expand("haplotyped_vcfs/{sample}.{chr_n}.g.vcf.gz", sample=XY_samples, chr_n=all_chromosomes),
#
## ------------------
## Read mapping rules
## ------------------
#rule minimap2_map_reads:
#	input:
#		fq1 = "reads/{sample}_R1.fastq.gz", 
#		fq2 = "reads/{sample}_R2.fastq.gz", 
#	output:
#		bam = "temp/{sample}.bam",
#	params:
#		threads = 4,
#		genome = genome,
#	shell:
#		"""
#		minimap2 -ax sr {params.genome} {input.fq1} {input.fq2} -t {params.threads} | samblaster | samtools fixmate -@ {params.threads} - - | samtools sort -@ {params.threads} -O bam - -o {output.bam} 2>/dev/null 
#		"""
#
#rule index_group_reads:
#	input:
#		"temp/{sample}.bam",
#	output:
#		bam = "mapped_reads/{sample}.bam",
#		bai = "mapped_reads/{sample}.bai",
#	params:
#		threads = 2,
#		id = lambda wildcards: config[wildcards.sample]["ID"],
#		sm = lambda wildcards: config[wildcards.sample]["SM"],
#		lb = lambda wildcards: config[wildcards.sample]["LB"],
#		pu = lambda wildcards: config[wildcards.sample]["PU"],
#		pl = lambda wildcards: config[wildcards.sample]["PL"],
#	shell:
#		"""
#		picard AddOrReplaceReadGroups -I {input} -O {output.bam} -SORT_ORDER coordinate -RGID {params.id} -RGLB {params.lb} -RGPL {params.pl} -RGPU {params.pu} -RGSM {params.sm} -CREATE_INDEX True
#		"""
#
#rule stats_dna:
#	input:
#		"mapped_reads/{sample}.bam",
#	output:
#		"mapped_reads/{sample}.bam.stats",
#	shell:
#		"""
#		samtools stats {input} | grep '^SN' | cut -f 2- > {output};
#		"""
#
## -----------------
## Call variants dna
## -----------------
#
#rule gatk_gvcfs:
#	input:
#		bam = "mapped_reads/{sample}.bam",
#		bai = "mapped_reads/{sample}.bam.bai",
#	output:
#		"haplotyped_vcfs/{sample}.{chr_n}.g.vcf.gz",
#	params:
#		chr_n = config["all_chromosomes"],
#		genome = genome,
#	threads:
#		4
#	shell:
#		"""
#		gatk HaplotypeCaller -R {params.genome} -I {input.bam} --intervals {params.chr_n} --emit-ref-confidence GVCF -O {output};
#		"""
#
#
#
#rule gatk_gvcf_chrX_PARs:
#	input:
#		bam = "mapped_reads/{sample}.bam",
#		bai = "mapped_reads/{sample}.bam.bai",
#	output:
#		"called_reads/{sample}.chrX.g.vcf.gz",
#	params:
#		chrX = config["chrX"],
#		genome = genome,
#	threads:
#		4
#	shell:
#		"""
#		gatk HaplotypeCaller -R {params.genome} -I {input.bam} -L {params.X_chrom} --emit-ref-confidence GVCF -O {output};
#		"""
#
#rule gatk_gvcf_chrX_nonPAR:
#	input:
#		bam = "mapped_reads/{sample}.bam",
#		bai = "mapped_reads/{sample}.bam.bai",
#	output:
#		"called_reads/{sample}.chrX_nonPAR.g.vcf.gz",
#	params:
#		chrX = config["chrX"],
#		genome = genome,
#	threads:
#		4
#	shell:
#		"""
#		gatk HaplotypeCaller -R {params.genome} -I {input.bam} -L {params.X_chrom} --emit-ref-confidence GVCF -O {output};
#		"""
#
#rule gatk_gvcf_chrY:
#	input:
#		bam = "mapped_reads/{sample}.bam",
#		bai = "mapped_reads/{sample}.bam.bai",
#	output:
#		"called_reads/{sample}.chrX.g.vcf.gz",
#	params:
#		chrY = config["chrY"],
#		genome = genome,
#	threads:
#		4
#	shell:
#		"""
#		gatk HaplotypeCaller -R {params.genome} -I {input.bam} -L {params.chrY} --emit-ref-confidence GVCF -O {output};
#		"""
#
#
#
#--sample-ploidy 1
#
#