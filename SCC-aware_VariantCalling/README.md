# Overall Description

These Snakemake workflows are for sex chromosome complement informed (SCC-aware) whole-genome short-read sequence alignment and variant calling. For each sample, the genotype, identified using SCC-check pipeline (or reported sex), is used to determine which sex chromosome complement version of the reference genome is used for alignment/pseudoalignment. For females (no Y chromosome), we use the reference genome with Y chromosome hard-masked. For males (samples with a Y chromosome), we use the reference genome with the pseudoautosomal regions (PARs) hard masked on chromosome Y. A JSON file is used to hold all the file information in the previous directory (same config for SCC-check).

minimum requirements 13Gb RAM, recommended > 16Gb RAM
