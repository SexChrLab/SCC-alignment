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
autosomes = config["autosomes"]

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
#minimap2_map_reads rule
		expand("temp/X/{X}.bam", X=X_samples),
#index_stat rules
		expand("mapped/X/{X}.bai", X=X_samples),
		expand("stats/X/{X}.bam.stats", X=X_samples),
#gatk_gvcfs rule
		expand("haplotyped_vcfs/X/{X}.{chr_n}.g.vcf.gz", X=X_samples, chr_n=diploid),
#gatk_combinegvcfs_auto rule
		expand("haplotyped_vcfs/X/{chr_n}.gatk.combined.gvcf.vcf.gz", chr_n=diploid),
#gatk_genotypegvcfs rule
		expand("genotyped_vcfs/X/{chr_n}.gatk.genotyped.raw.vcf.gz", chr_n=diploid),
#hard_filter rule
		expand("genotyped_vcfs/X/{chr_n}.gatk.filtered.vcf.gz", chr_n=diploid),
#subset_individuals_hard_filter_diploid rule
		expand("genotyped_vcfs/X/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi", sample=X_samples, chr_n=diploid),
#gatk_selectheterozygous_diploid rule
		expand("diploid_vcfs/X/{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz", sample=X_samples, chr_n=diploid),

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
		1
	run:
		variant_files = []
		for i in input:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""
			gatk --java-options "-Xmx4g" CombineGVCFs -R {params.X_genome} {variant_files} --intervals {params.chr_n} -O {output}
			""")

rule gatk_genotypegvcf_diploid:
	input:
		"haplotyped_vcfs/X/{chr_n}.gatk.combined.gvcf.vcf.gz",
	output:
		"genotyped_vcfs/X/{chr_n}.gatk.genotyped.raw.vcf.gz",
	params:
		X_genome = X_genome,
	threads:
		1
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
	shell:
		"""
		gatk SelectVariants -R {params.X_genome} -V {input} -O {output} \
		#--select-type-to-include SNP --restrict-alleles-to BIALLELIC -select 
		"AN >= 4 && MQ > 10.0 && QD > 7.0 && DP >= 10.0 && DP <= 1000.0"
		touch -c {output.idx};
		"""

rule subset_individuals_hard_filter_diploid:
	input:
		"genotyped_vcfs/X/{chr_n}.gatk.filtered.vcf.gz",
	output:
		vcf = "genotyped_vcfs/X/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz",
		index = "genotyped_vcfs/X/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz.tbi",
	params:
		sample = lambda wildcards: config[wildcards.sample]["SM"],
	shell:
		"""
		bcftools view -Oz -s {params.sample} {input} > {output.vcf};
		tabix -p vcf {output.vcf};
		"""

#After subsetting for each individual. In some individuals, the genotypes could be homozygous for the reference. This next rule is to remove these sites.

rule gatk_selectheterozygous_diploid:
	input:
		"genotyped_vcfs/X/{chr_n}.gatk.called.hard.filter.{sample}.vcf.gz", 
	output:
		"diploid_vcfs/X/{chr_n}.gatk.called.hard.filter.het.{sample}.vcf.gz",
	params:
		X_genome = X_genome,
	shell:
		"""
		gatk SelectVariants -R {params.X_genome} -V {input} -O {output} -select "AC == 1";
		"""
