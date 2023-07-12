# Instructions to create sex chromosome complement (SCC) reference genomes 

This section is a walk through how we generated sex chromosome complement (SCC) reference genomes.  This shows the steps for how we generated SCC versions of the latest human genome reference, CHM13v2 (telomere-to-telomere reference).  We used a similar process to generate SCC references from the HG38 (GRCh38) and transcriptome sequences and similar steps can be used to generate SCC versions of other reference genomes where the sex chromosomes have been chararacterized. 

# Output
The following are steps to making two versions of the human reference genome: 
1. A reference genome with pseudoautosomal regions 1 and 2 (PAR1 and PAR2) on the Y chromosome hard masked (replace with Ns) for alignment of sequencing reads in samples with a Y chromosome (such as XY males)
2. A reference genome where the entire Y chromosome is hard masked for alignment of sequencing reads in samples with no Y chromosome (such as XX females)

# Download reference genome sequence

The complete human reference genome sequence can be downloaded: 
https://www.ncbi.nlm.nih.gov/assembly/GCF_009914755.1/?shouldredirect=false
Select Genome FASTA from the Download Assembly menu.

# Identify PAR coordinates on X and Y
The coordinates of the PAR regions are provided with certain releases of the human reference genome sequence, but below we will show one method that can be used to figure out these coordinates. 

PAR coordinates will be identified by running `lastz` between the full sequences of chromosomes X and Y, retaining alignments with 100% sequence identity, then visually assessing the dotplot to see where the one to one alignments at each end of the chromosomes start and end.

```
# First extract chromosomes X and Y from reference genome
cd ../
samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna

grep ">" GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna
>CP068277.2 Homo sapiens isolate CHM13 chromosome 1
>CP068276.2 Homo sapiens isolate CHM13 chromosome 2
>CP068275.2 Homo sapiens isolate CHM13 chromosome 3
>CP068274.2 Homo sapiens isolate CHM13 chromosome 4
>CP068273.2 Homo sapiens isolate CHM13 chromosome 5
>CP068272.2 Homo sapiens isolate CHM13 chromosome 6
>CP068271.2 Homo sapiens isolate CHM13 chromosome 7
>CP068270.2 Homo sapiens isolate CHM13 chromosome 8
>CP068269.2 Homo sapiens isolate CHM13 chromosome 9
>CP068268.2 Homo sapiens isolate CHM13 chromosome 10
>CP068267.2 Homo sapiens isolate CHM13 chromosome 11
>CP068266.2 Homo sapiens isolate CHM13 chromosome 12
>CP068265.2 Homo sapiens isolate CHM13 chromosome 13
>CP068264.2 Homo sapiens isolate CHM13 chromosome 14
>CP068263.2 Homo sapiens isolate CHM13 chromosome 15
>CP068262.2 Homo sapiens isolate CHM13 chromosome 16
>CP068261.2 Homo sapiens isolate CHM13 chromosome 17
>CP068260.2 Homo sapiens isolate CHM13 chromosome 18
>CP068259.2 Homo sapiens isolate CHM13 chromosome 19
>CP068258.2 Homo sapiens isolate CHM13 chromosome 20
>CP068257.2 Homo sapiens isolate CHM13 chromosome 21
>CP068256.2 Homo sapiens isolate CHM13 chromosome 22
>CP068255.2 Homo sapiens isolate CHM13 chromosome X
>CP086569.2 Homo sapiens isolate NA24385 chromosome Y
>CP068254.1 Homo sapiens isolate CHM13 mitochondrion, complete genome

# CP068255.2 is chromosome X and CP086569.2 is chromosome Y

samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna CP068255.2 > T2T_CHM13_v2_SCC/T2T_chrX.fa
samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna CP086569.2 > T2T_CHM13_v2_SCC/T2T_chrY.fa

less GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna.fai
CP068277.2      248387328       52      80      81
CP068276.2      242696752       251492274       80      81
CP068275.2      201105948       497222788       80      81
CP068274.2      193574945       700842613       80      81
CP068273.2      182045439       896837297       80      81
CP068272.2      172126628       1081158356      80      81
CP068271.2      160567428       1255436619      80      81
CP068270.2      146259331       1418011192      80      81
CP068269.2      150617247       1566098817      80      81
CP068268.2      134758134       1718598833      80      81
CP068267.2      135127769       1855041497      80      81
CP068266.2      133324548       1991858417      80      81
CP068265.2      113566686       2126849575      80      81
CP068264.2      101161492       2241835898      80      81
CP068263.2      99753195        2344261962      80      81
CP068262.2      96330374        2445262125      80      81
CP068261.2      84276897        2542796682      80      81
CP068260.2      80542538        2628127094      80      81
CP068259.2      61707364        2709676467      80      81
CP068258.2      66210255        2772155227      80      81
CP068257.2      45090682        2839193164      80      81
CP068256.2      51324926        2884847533      80      81
CP068255.2      154259566       2936814073      80      81
CP086569.2      62460029        3093001938      80      81
CP068254.1      16569   3156242788      80      81


# Then run lastZ
lastz T2T_chrX.fa T2T_chrY.fa --identity=100 --exact=50 --ambiguous=iupac --notransition --nogapped --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_lastz_identity100_exact50_ambiguous_iupac_notransition_nogapped_step10.dotplot

```

We can determine the lengths of the full X and Y chromosomes.
```
# chrX
grep CP068255.2 ../GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna.fai
CP068255.2 154259566       2936814073      80      81

# chrY
grep CP086569.2 ../GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna.fai
CP086569.2 62460029        3093001938      80      81

```

From this we determine the coordinates of the two pseudoautosomal (PAR) regions on chromosomes X and Y:
PAR1:

- X chromosome: 1 - 2394365
- Y chromosome: 1 - 2458275

PAR2:

- X chromosome: 153928134 - 154259566
- Y chromosome: 62125118 - 62460029


# Mask Y chromosome and Y PARs
Run `bedtools maskfasta` to replace the PAR sequences on chromosome Y with Ns and the entire Y chromosome.

```
# Make bed files for the Y PAR coordinates and entire Y chromosome
vi T2T_chrY_PARs.txt
awk '{ print $1"\t"$2"\t"$3 }' T2T_chrY_PARs.txt > T2T_chrY_PARs.bed

vi T2T_chrY.txt
awk '{ print $1"\t"$2"\t"$3 }' T2T_chrY.txt > T2T_chrY.bed

bedtools maskfasta -fi T2T_chrY.fa -bed T2T_chrY_PARs.bed -fo T2T_chrY_YPARs_masked.fa

bedtools maskfasta -fi T2T_chrY.fa -bed T2T_chrY.bed -fo T2T_chrY_YHardMasked.fa

# Extract chr 1-22, M and X from T2T reference.
samtools faidx ../GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna CP068277.2 CP068276.2 CP068275.2 CP068274.2 CP068273.2 CP068272.2 CP068271.2 CP068270.2 CP068269.2 CP068268.2 CP068267.2 CP068266.2 CP068265.2 CP068264.2 CP068263.2 CP068262.2 CP068261.2 CP068260.2 CP068259.2 CP068258.2 CP068257.2 CP068256.2 CP068255.2 CP068254.1 > GCA_009914755.4_CHM13_T2T_v2.0_genomic_chr1-22_chrX_chrM.fa

# Merge Y hard masked fa and Y PARs masked fa, separately to make each of the SCC references.

# YHardMasked
cat GCA_009914755.4_CHM13_T2T_v2.0_genomic_chr1-22_chrX_chrM.fa T2T_chrY_YHardMasked.fasta > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_unsorted.fa

# YPARsMasked
cat GCA_009914755.4_CHM13_T2T_v2.0_genomic_chr1-22_chrX_chrM.fa T2T_chrY_YPARs_masked.fasta > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_unsorted.fa

# Reorder the chromosomes in each reference (chr 1-22, X, Y, M)

# YHardMasked
samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_unsorted.fa

samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_unsorted.fa CP068277.2 CP068276.2 CP068275.2 CP068274.2 CP068273.2 CP068272.2 CP068271.2 CP068270.2 CP068269.2 CP068268.2 CP068267.2 CP068266.2 CP068265.2 CP068264.2 CP068263.2 CP068262.2 CP068261.2 CP068260.2 CP068259.2 CP068258.2 CP068257.2 CP068256.2 CP068255.2 CP086569.2 CP068254.1 > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked.fa

# YPARsMasked
samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_unsorted.fa

samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_unsorted.fa CP068277.2 CP068276.2 CP068275.2 CP068274.2 CP068273.2 CP068272.2 CP068271.2 CP068270.2 CP068269.2 CP068268.2 CP068267.2 CP068266.2 CP068265.2 CP068264.2 CP068263.2 CP068262.2 CP068261.2 CP068260.2 CP068259.2 CP068258.2 CP068257.2 CP068256.2 CP068255.2 CP086569.2 CP068254.1 > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked.fa

# Rename the chromosomes in reference
sed 's/CP068277.2/chr1/g' GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked.fa  | sed 's/CP068276.2/chr2/g' | sed 's/CP068275.2/chr3/g' | sed 's/CP068274.2/chr4/g' | sed 's/CP068273.2/chr5/g' | sed 's/CP068272.2/chr6/g' | sed 's/CP068271.2/chr7/g' | sed 's/CP068270.2/chr8/g' | sed 's/CP068269.2/chr9/g' | sed 's/CP068268.2/chr10/g' | sed 's/CP068267.2/chr11/g' | sed 's/CP068266.2/chr12/g' | sed 's/CP068265.2/chr13/g' | sed 's/CP068264.2/chr14/g' | sed 's/CP068263.2/chr15/g' | sed 's/CP068262.2/chr16/g' | sed 's/CP068261.2/chr17/g' | sed 's/CP068260.2/chr18/g' | sed 's/CP068259.2/chr19/g' | sed 's/CP068258.2/chr20/g' | sed 's/CP068257.2/chr21/g' | sed 's/CP068256.2/chr22/g' | sed 's/CP068255.2/chrX/g' | sed 's/CP086569.2/chrY/g' | sed 's/CP068254.1/chrM/g' > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa

sed 's/CP068277.2/chr1/g' GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked.fa  | sed 's/CP068276.2/chr2/g' | sed 's/CP068275.2/chr3/g' | sed 's/CP068274.2/chr4/g' | sed 's/CP068273.2/chr5/g' | sed 's/CP068272.2/chr6/g' | sed 's/CP068271.2/chr7/g' | sed 's/CP068270.2/chr8/g' | sed 's/CP068269.2/chr9/g' | sed 's/CP068268.2/chr10/g' | sed 's/CP068267.2/chr11/g' | sed 's/CP068266.2/chr12/g' | sed 's/CP068265.2/chr13/g' | sed 's/CP068264.2/chr14/g' | sed 's/CP068263.2/chr15/g' | sed 's/CP068262.2/chr16/g' | sed 's/CP068261.2/chr17/g' | sed 's/CP068260.2/chr18/g' | sed 's/CP068259.2/chr19/g' | sed 's/CP068258.2/chr20/g' | sed 's/CP068257.2/chr21/g' | sed 's/CP068256.2/chr22/g' | sed 's/CP068255.2/chrX/g' | sed 's/CP086569.2/chrY/g' | sed 's/CP068254.1/chrM/g' > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_ChrNamesAdded.fa

grep ">" GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna
>CP068277.2 Homo sapiens isolate CHM13 chromosome 1
>CP068276.2 Homo sapiens isolate CHM13 chromosome 2
>CP068275.2 Homo sapiens isolate CHM13 chromosome 3
>CP068274.2 Homo sapiens isolate CHM13 chromosome 4
>CP068273.2 Homo sapiens isolate CHM13 chromosome 5
>CP068272.2 Homo sapiens isolate CHM13 chromosome 6
>CP068271.2 Homo sapiens isolate CHM13 chromosome 7
>CP068270.2 Homo sapiens isolate CHM13 chromosome 8
>CP068269.2 Homo sapiens isolate CHM13 chromosome 9
>CP068268.2 Homo sapiens isolate CHM13 chromosome 10
>CP068267.2 Homo sapiens isolate CHM13 chromosome 11
>CP068266.2 Homo sapiens isolate CHM13 chromosome 12
>CP068265.2 Homo sapiens isolate CHM13 chromosome 13
>CP068264.2 Homo sapiens isolate CHM13 chromosome 14
>CP068263.2 Homo sapiens isolate CHM13 chromosome 15
>CP068262.2 Homo sapiens isolate CHM13 chromosome 16
>CP068261.2 Homo sapiens isolate CHM13 chromosome 17
>CP068260.2 Homo sapiens isolate CHM13 chromosome 18
>CP068259.2 Homo sapiens isolate CHM13 chromosome 19
>CP068258.2 Homo sapiens isolate CHM13 chromosome 20
>CP068257.2 Homo sapiens isolate CHM13 chromosome 21
>CP068256.2 Homo sapiens isolate CHM13 chromosome 22
>CP068255.2 Homo sapiens isolate CHM13 chromosome X
>CP086569.2 Homo sapiens isolate NA24385 chromosome Y
>CP068254.1 Homo sapiens isolate CHM13 mitochondrion, complete genome

```

