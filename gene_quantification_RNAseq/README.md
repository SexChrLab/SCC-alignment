Pipe hisat w/ samtools and skip picard

# Overall Description

These Snakemake workflows are for sex chromosome complement informed gene quantification for short read, paired end RNA sequencing data.  For each sample, the sex is used to determine which sex chromosome complement version of the reference genome is used for alignment/pseudoalignment.  For females (XX, no Y chromosome),  we use the reference genome with chromosome Y hard masked.  For males (XY, with a Y chromosome), we use the reference genome with the pseudoautosomal regions (PARs) hard masked on chromosome Y.  There are separate Snakemake pipelines for doing a alignment (hisat2) with gene quantification (featureCounts) or doing pseudoalignment and quantification (salmon) which is much faster.  A JSON file is used to hold all the file information; a script that can be used to generate the sample specific information is provided.  

# Steps to running these workflows

Set up: 
1) Install conda (or use Docker)
2) Set up conda environment with provided environment yaml (`conda env create -f SCCalign_v3.yml`)

Analyze:
1) Generate configuration JSON used to specify the location of files and sex chromosome complement reference genomes
2) Decide what method you would like to use to quantify gene expression
a) full alignment with gene counts (reads mapping to each genes coordinates as a measure of expression) [long run time]
b) pseudoalignment and quantification with sex chromosome complement reference transcriptome [short run time]
3) Use Snakemake in conda environment `SCCalign_v3` to run the provided Snakemake workflow for your desired method (examples provided in `run_*.sbatch`)

# Generating a configuration JSON

The Python script `generate_json_config_addreadgroups.py` can be used to generate the sample information elements of a configuration JSON (see example RNAseq_samples.json).  To use this script, open it in an editor and change the variables under the comment `"Specifying inputs:"`.  The `all_sample_ids` variable should be set to a comma separated file that contains the sample ID, forward read fastq file, and reverse read fastq file (example provided in `samples.csv`).  If the files are in subfolders, include that in the path to the file.  The `out_config_json` variable should be set to the desired name of the output JSON to be generated.  The `input_directory` should be set to the directory containing the read sequence files (fastq/fastq.gz).  Thus, later in the workflow, the full path to the sequence files will be the `input_directory` variable contents concatenated to the path to the read_file.  Once these variables have been set, run the python script and check the output.  You should see a list of all the samples IDs in an entry called `ALL_samples` and then each sample will have a entry that contains the path information for the sequence files, ID, and read groups.  

Once that all looks good, you will need to additional information to the configuration json.  

1) add the paths to the reference genomes (template for this is given in `additional_info_config.json`, entries can be copied and editted according to where you are running the scripts)
2) `X_samples` and `Y_samples` lists can be filled using the results of the quick sex check or reported sex (samples can be copied and pasted from the `ALL_samples` list based on whether they were predicted to be XX into `X_samples` or predicted to be XY into `Y_samples`)


# Alignment and gene quantification

The Snakemake workflow `rnaseq_data_processing_hisat2.snakefile` does an alignment to sex chromsoome complement refernece genome using `hisat2` and quantifies gene expression using the genome annotation using `featureCounts`.  Interrim files produced during alignment and file processing are marked as temporary, so they will be generated and deleted when the step that needs them is completed so as not to fill up hard drive space unnecessarily.  The rules for running `featureCounts` are set up to output results quantifying gene expression on the exon regions and output with both the gene IDs and the transcript IDs.

# Pseudoalignment and transcript quantification

The Snakemake workflow `rna_data_processing_salmon.snakefile` uses the `salmon` algorithm to quantify expression.  Salmon indices are provided in `references/` for reference transcriptome with sequences from chromosome Y masked, so these are used according to the sex of the sample (as specified in the `X_samples` and `Y_samples` entries in the configuration JSON).  Quantified expression will be in the `quantified_rna_salmon` directory, with specific subdirectories for each sample.  The `quant.sf` file inside each result subdirectory will give the quantified expression for the transcripts in each sample.
