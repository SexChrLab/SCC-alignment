# Custom Config JSON

These scripts can be used to help you to generate custom config JSONs to hold the required information to run sex chromosome complement (SCC) methods on your sequencing data.  These scripts take a table of sample IDs with the sequencing files associated with those samples as input and produces a config json containing sample information and place holders for other variables used in the sex chromosome complement aware alignment protocols.  It is intended to have a separate config for DNA based analyses (`SCC_check` for DNA samples and `SCC-aware_VariantCalling`) and RNA based analysis (`gene_quantification_RNAseq`).

# Steps

Separately for each set of DNA or RNA samples:
1) Create sample info table csv
2) Run `generate_custom_json_*.py` script

# Sample info table file

As in the provide example files `RNA_samples.csv` and `DNA_samples.csv`, to run the custom config scripts, you must generate a comma separated file with no headers wherein each row contains: 

1) Sample ID:  a string used to identify a specific sample
2) Forward read fastq file: name of the file containing the forward reads for the sample
3) Reverse read fastq file: name of the file containing the reverse reads for the sample

You will specify the directory these fastq files are found in when you run the script to generate the custom config JSON (below).  When entering fields 2 and 3 in the lines, you can specify any subdirectory structure.  That is, for example, if all files are in `/data/rnaseq/files/` but the sequencing files are in subdirectories inside that main directory, you can enter `subdir1/sample1_R1.fastq.gz` and `subdir1/sample1_R2.fastq.gz`, you can enter the following row into your input table csv: 

`sample1,subdir1/sample1_R1.fastq.gz,subdir1/sample1_R1.fastq.gz`
 
And indicate `/data/rnaseq/files/` as the directory of the sequencing files in the next step (below).  This will tell the analysis pipeline to read your sample sequencing files from `/data/rnaseq/files/subdir1/sample1_R1.fastq.gz' and `/data/rnaseq/files/subdir1/sample1_R2.fastq.gz`.  

# DNA

For DNA analyses (SCC check and variant calling), create sample info table (comma separate values CSV) for your DNA sequencing samples.

The `generate_custom_json_DNA.py` Python script contains information needed to make a config JSON for SCC check with DNA sequencing data and for variant calling.

To run `generate_custom_json_DNA.py`:
1) open the script with a text editor
2) change the `all_sample_ids` variable on line 9 to include the name of your sample input table file csv
3) change the `inferred_SCC` variable on line 10 to the path to the output of SCC_check if you have already run it so that inferred SCC can be used to populate the samples with and without a Y chromosome for appropriate reference genome selection. Indicate an empty string if you have not done this step.
4) change the `out_config_json` variable on line 11 to the path to which your output JSON will be written.
5) change the `input_directory` variable to the directory containing your sequencing files (see above)
6) save changes to your script

Once these have been indicated, run `python generate_custom_json_DNA.py` on a terminal where python is installed (such as when the provided conda environment SCCalign_v3 is activated).

After this script is executed, open your output JSON in a text editor.  The output is in JSON format, which is essentially a text representation of a dictionary in Python, that will be used to indicate important information needed for SCC aware DNA analyses.  You should see something like the provided example file `DNA_samples.config.json`.  Some fields will be filled in using by the script reading information from the sequencing files, others will have to be filled in by hand to reflect the specifics of your own system.

To complete your DNA config JSON: 
1) Double check that the elements that have been filled look correct
2) Edit the paths to references point to your own paths for reference files
3) The field `Y_samples` should have samples that contain a Y chromosomes and `X_samples` should have samples that do not have a Y chromosome
a) It's best to copy from the `ALL_samples` field to avoid typos
b) These fields will be read when doing variant calling analysis so that alignment can be done with the appropriate SCC reference
c) If you are using this script to generate a config json before doing an SCC check, you can leave these blank and fill them in after the SCC check is complete
4) Fill in the path to your sequencing files

# RNA

For RNA sequencing analyses (alignment and gene quantification), create sample info table (comma separate values CSV) for your RNA sequencing samples.

The `generate_custom_json_RNA.py` script contains information needed to make a config JSON for RNA sequencing data for gene quantitation, including referneces for hisat2 or salmon.  

To run `generate_custom_json_RNA.py`:
1) open the script with a text editor
2) change the `all_sample_ids` variable on line 9 to include the name of your sample input table file csv
3) change the `out_config_json` variable on line 11 to the path to which your output JSON will be written.
4) change the `input_directory` variable to the directory containing your sequencing files (see above)
5) save changes to your script

Once these have been indicated, run the script on a terminal where python is installed (such as when the provided conda environment SCCalign_v3 is activated):
```
python generate_custom_json_RNA.py
```

After this script is executed, open your output JSON in a text editor.  You should see something like the provided file `RNA_samples.config.json`.  Some fields will be filled in using by the script reading information from the sequencing files, others will have to be filled in by hand to reflect the specifics of your own system.

To complete your RNA config JSON: 
1) Double check that the elements that have been filled look correct
2) Edit the paths to HISAT2 or Salmon indexed references to point to your own paths for reference files
3) The field `Y_samples` should have samples that contain a Y chromosomes and `X_samples` should have samples that do not have a Y chromosome
a) It's best to copy from the `ALL_samples` field to avoid typos
4) Fill in the path to your sequencing files
5) Fill in the options for HISAT2 or Salmon analysis (check the manuals for these programs for descriptions of these options (see below))

# Manuals

`HISAT2` manual: 
https://daehwankimlab.github.io/hisat2/manual/

`featureCounts` manual (`subread` package)
https://subread.sourceforge.net/

Salmon manual: 
https://salmon.readthedocs.io/en/latest/salmon.html
