# Overall Description

This Snakemake workflow is for identifying the sex chromosome complement of a set of samples using whole-genome short-read sequencing data using comparative read depth. This information should be then be incorporated into the JSON config file (SCC-analysis_config.json) for use in SCC-aware variant calling.

# Rationale

It is important to know if a sample possesses a Y chromosome or not to understand potential biases this can introduce into traditional read mapping-based genomic analyses. Here, were subsample the read files to ~1X coverage (assuming PE 150bp short reads) and quickly map these read to a Y PAR-masked reference. If the sample possesses a Y chromosome the chrY copy number (CN) will be greater than ~0. This pipeline is built for speed to encourage its use, but is not accurate for determining chrX CN. If chrX CN important for your use case, simply align all the reads (20-30X) and run indexcov following the programs documentation (https://github.com/brentp/goleft/tree/master/indexcov).

# Running the workflow

Add more here



# Citations 

If you use this script in your work, please cite the relevant softwares that make it work:

seqkit: https://github.com/lh3/seqtk

Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191

SAMtools: Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. https://doi.org/10.1093/bioinformatics/btp352

Indexcov: Pedersen, B. S., Collins, R. L., Talkowski, M. E., & Quinlan, A. R. (2017). Indexcov: fast coverage quality control for whole-genome sequencing. GigaScience, 6(11). https://doi.org/10.1093/gigascience/gix090
