import os

#mamba env create --name SCCalign_v3 -f SCCalign_v3.yml
#mamba activate SCCalign_v3
#snakemake --use-conda --cores 4 -s Snakefile_XX_SNPs.py -np

#shorthand:
#X == no chrY present
#Y == chrY present

configfile: "SCC-analysis_config.json"

#threads
threads = config["threads"]

#sample sets
X_samples = config["X_samples"]
Y_samples = config["Y_samples"]
ALL_samples = config["ALL_samples"]
#dir = config["directory"]

#chromosomes/intervals
all_chromosomes = config["all_chromosomes"]
diploid = config["XX_diploid"]
autosomes = config["autosomes"]

#Filtering options
AN = config["diploid_AN"],
MQ = config["diploid_MQ"],
QD = config["diploid_QD"],
DP1 = config["diploid_DP1"],
DP2 = config["diploid_DP2"],

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
########################Stage 1: Map reads and QC for all Y- samples
#minimap2_map_reads rule
		expand("temp/X/{X}.bam", X=X_samples),
#index_stat rules
		expand("mapped/X/{X}.bai", X=X_samples),
		expand("stats/X/{X}.bam.stats", X=X_samples),

########################Stage 2: Call haplotypes for all genomic regions (ploidy-aware)
#HaplotypeCaller
#gatk_gvcfs rule
		expand("haplotyped_vcfs/X/{X}.{chr_n}.g.vcf.gz", X=X_samples, chr_n=diploid),

########################Stage 3: Merge gVCFs for all genomic regions
#CombineGVCFs
#gatk_combinegvcfs_auto rule
		expand("haplotyped_vcfs/X/{chr_n}.gatk.combined.gvcf.vcf.gz", chr_n=diploid),

########################Stage 4: Compute joint genotypes across all samples (ploidy-aware)
#GenotypeGVCFs
#gatk_genotypegvcfs rule
		expand("genotyped_vcfs/X/{chr_n}.gatk.genotyped.raw.vcf.gz", chr_n=diploid),

########################Stage 6: Downstream analysis-aware VCF filtering
#hard_filter rule
		expand("genotyped_vcfs/X/{chr_n}.gatk.filtered.vcf.gz", chr_n=diploid),
#subset_individuals_hard_filter_diploid rule
#		expand("genotyped_vcfs/X/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi", sample=X_samples, chr_n=diploid),
#gatk_selectheterozygous_diploid rule
#		expand("diploid_vcfs/X/{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz", sample=X_samples, chr_n=diploid),

# ------------------
# Read mapping rules
# ------------------
rule minimap2_mapping_X:
	input:
		fq1 = "reads/{X}_R1.fastq.gz", 
		fq2 = "reads/{X}_R2.fastq.gz", 
	output:
		"temp/X/{X}.bam",
	params:
		threads = threads,
		X_genome = X_genome,
	shell:
		"""
		minimap2 -ax sr {params.X_genome} {input.fq1} {input.fq2} -t {params.threads} | samblaster | samtools fixmate -@ {params.threads} - - | samtools sort -@ {params.threads} -O bam - -o {output} 2>/dev/null 
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
		threads = threads,
		id = lambda wildcards: config[wildcards.X]["ID"],
		sm = lambda wildcards: config[wildcards.X]["SM"],
		lb = lambda wildcards: config[wildcards.X]["LB"],
		pu = lambda wildcards: config[wildcards.X]["PU"],
		pl = lambda wildcards: config[wildcards.X]["PL"],
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
		"stats/X/{sample}.bam.stats",
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

## ---------------------------------------
## Combine all gvcfs
## ---------------------------------------

rule gatk_combinegvcfs_diploid:
	input:
		expand("haplotyped_vcfs/X/{X}.{chr_n}.g.vcf.gz", X=X_samples, chr_n=diploid),
	output:
		"haplotyped_vcfs/X/{chr_n}.gatk.combined.gvcf.vcf.gz",
	params:
		X_genome = X_genome,
		chr_n = "{chr_n}",
	threads:
		2
	run:
		variant_files = []
		for i in input:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.X_genome} {variant_files} --intervals {params.chr_n} -O {output}
			""")

## ---------------------------------------
## Genotype XX samples
## ---------------------------------------

rule gatk_genotypegvcf_diploid:
	input:
		"haplotyped_vcfs/X/{chr_n}.gatk.combined.gvcf.vcf.gz",
	output:
		"genotyped_vcfs/X/{chr_n}.gatk.genotyped.raw.vcf.gz",
	params:
		X_genome = X_genome,
	threads:
		2
	shell:
		"""
		gatk --java-options "-Xmx4g" GenotypeGVCFs -R {params.X_genome} -V {input} -O {output}
		"""

## ------------------------
## Hard Filter on diploid_chromosomes
## ------------------------

rule hard_filter_diploid:
	input:
		"genotyped_vcfs/X/{chr_n}.gatk.genotyped.raw.vcf.gz",
	output:
		vcf = "genotyped_vcfs/X/{chr_n}.gatk.filtered.vcf.gz",
		idx = "genotyped_vcfs/X/{chr_n}.gatk.filtered.vcf.gz.tbi",
	params:
		X_genome = X_genome,
		AN = AN,
		MQ = MQ,
		QD = QD,
		DP1 = DP1,
		DP2 = DP2,
	shell:
		"""
		gatk SelectVariants -R {params.X_genome} -V {input} -O {output.vcf} --select-type-to-include SNP --restrict-alleles-to BIALLELIC -select AN >= {params.AN} && MQ > {params.MQ} && QD > {params.QD} && DP >= {params.DP1} && DP <= {params.DP2}"
		touch -c {output.idx};
		"""

#rule subset_individuals_hard_filter_diploid:
#	input:
#		"genotyped_vcfs/X/{chr_n}.gatk.filtered.vcf.gz",
#	output:
#		vcf = "genotyped_vcfs/X/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz",
#		index = "genotyped_vcfs/X/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi",
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
#		"genotyped_vcfs/X/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz", 
#	output:
#		"diploid_vcfs/X/{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz",
#	params:
#		X_genome = X_genome,
#	shell:
#		"""
#		gatk SelectVariants -R {params.X_genome} -V {input} -O {output} -select "AC == 1";
#		"""
