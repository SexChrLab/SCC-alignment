# Overall Description

These Snakemake workflows are for sex chromosome complement informed (SCC-aware) whole-genome short-read sequence alignment and variant calling. For each sample, the genotype, identified using SCC-check pipeline (or reported sex), is used to determine which sex chromosome complement version of the reference genome is used for alignment/pseudoalignment. For females (no Y chromosome), we use the reference genome with Y chromosome hard-masked. For males (samples with a Y chromosome), we use the reference genome with the pseudoautosomal regions (PARs) hard masked on chromosome Y. A JSON file is used to hold all the file information in the previous directory (same config for SCC_check).

# System requirements

Minimum: 4 available threads; Recommended: >12 threads.

Minimum: ~16Gb RAM; Recommended: >24Gb RAM.

# Walkthrough
















# Citations 
If you use this code in your work, please cite its source and tools:

Plasier et al. 202x. SCC-aware analysis tutorial.

bwa: McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., & DePristo, M. A. (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. Genome Research, 20(9), 1297–1303. https://doi.org/10.1101/gr.107524.110

or 
 
Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191,

and 

SAMtools: Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352,

GATK: McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome research, 20(9), 1297-1303.  http://www.genome.org/cgi/doi/10.1101/gr.107524.110.



