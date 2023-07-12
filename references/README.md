# README.md


## Basic information
**Created by:** Angela Oill (amtarave@asu.edu)

**Created on:** 02-17-2022

**Last modified:** 03-02-2022

**Last modified by:** Elizabeth Borden (knodele@email.arizona.edu)

**Path to README.md:** `/data/CEM/shared/public_data/references/T2T_CHM13_v2/T2T_CHM13_v2_SCC/README.md`


## Background
This readme will serve to outline the steps to making sex chromosome complement versions of the telomere to telomere reference genome. For this, I will make two versions:
1. A reference genome with pseudoautosomal regions 1 and 2 (PAR1 and PAR2) on the Y chromosome hard masked (replace with Ns). This should be used when aligning XY individuals.
2. A reference genome where the entire Y chromosome is hard masked. This should be used when aligning XX individuals.

## Identify PAR coordinates on X and Y
PAR coordinates will be identified by running lastz between the telomere to telomere X and Y, retaining alignments with 100% sequence identity, then visually assessing the dotplot to see where the one to one alignments at each end of the chromosomes start and end.

Run lastz:
```
# First extract X and Y from reference genome
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

Lengths of the chromosomes:
```
# chrX
grep CP068255.2 ../GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna.fai
CP068255.2 154259566       2936814073      80      81

# chrY
grep CP086569.2 ../GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna.fai
CP086569.2 62460029        3093001938      80      81

```

**PAR1:**
- X chromosome: 1 - 2394365
- Y chromosome: 1 - 2458275

**PAR2:**
- X chromosome: 153928134 - 154259566
- Y chromosome: 62125118 - 62460029


## Mask Y chromosome and Y PARs
Run `bedtools maskfasta` to replace the PAR sequences on chromosome Y with Ns and the entire Y chromosome.

```
# Make bed files for the Y PAR coordinates and entire Y chromosome
vi T2T_chrY_PARs.txt
awk '{ print $1"\t"$2"\t"$3 }' T2T_chrY_PARs.txt > T2T_chrY_PARs.bed

vi T2T_chrY.txt
awk '{ print $1"\t"$2"\t"$3 }' T2T_chrY.txt > T2T_chrY.bed

bedtools maskfasta -fi T2T_chrY.fa -bed T2T_chrY_PARs.bed -fo T2T_chrY_YPARs_masked.fa

bedtools maskfasta -fi T2T_chrY.fa -bed T2T_chrY.bed -fo T2T_chrY_YHardMasked.fa
```

Extract chr 1-22, M and X from T2T reference.
```
samtools faidx ../GCA_009914755.4_CHM13_T2T_v2.0_genomic.fna CP068277.2 CP068276.2 CP068275.2 CP068274.2 CP068273.2 CP068272.2 CP068271.2 CP068270.2 CP068269.2 CP068268.2 CP068267.2 CP068266.2 CP068265.2 CP068264.2 CP068263.2 CP068262.2 CP068261.2 CP068260.2 CP068259.2 CP068258.2 CP068257.2 CP068256.2 CP068255.2 CP068254.1 > GCA_009914755.4_CHM13_T2T_v2.0_genomic_chr1-22_chrX_chrM.fa

```

Merge Y hard masked fa and Y PARs masked fa, separately to make each of the SCC references.
```
# YHardMasked
cat GCA_009914755.4_CHM13_T2T_v2.0_genomic_chr1-22_chrX_chrM.fa T2T_chrY_YHardMasked.fasta > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_unsorted.fa

# YPARsMasked
cat GCA_009914755.4_CHM13_T2T_v2.0_genomic_chr1-22_chrX_chrM.fa T2T_chrY_YPARs_masked.fasta > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_unsorted.fa

```

Reorder the chromosomes in each reference (chr 1-22, X, Y, M)
```
# YHardMasked
samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_unsorted.fa

samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_unsorted.fa CP068277.2 CP068276.2 CP068275.2 CP068274.2 CP068273.2 CP068272.2 CP068271.2 CP068270.2 CP068269.2 CP068268.2 CP068267.2 CP068266.2 CP068265.2 CP068264.2 CP068263.2 CP068262.2 CP068261.2 CP068260.2 CP068259.2 CP068258.2 CP068257.2 CP068256.2 CP068255.2 CP086569.2 CP068254.1 > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked.fa

# YPARsMasked
samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_unsorted.fa

samtools faidx GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked_unsorted.fa CP068277.2 CP068276.2 CP068275.2 CP068274.2 CP068273.2 CP068272.2 CP068271.2 CP068270.2 CP068269.2 CP068268.2 CP068267.2 CP068266.2 CP068265.2 CP068264.2 CP068263.2 CP068262.2 CP068261.2 CP068260.2 CP068259.2 CP068258.2 CP068257.2 CP068256.2 CP068255.2 CP086569.2 CP068254.1 > GCA_009914755.4_CHM13_T2T_v2.0_genomic_YPARsMasked.fa

```

Rename the chromosomes in reference
```
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

## Identify XTR coordinates on X and Y

```
lastz T2T_chrX.fa T2T_chrY.fa --identity=98 --exact=50 --ambiguous=iupac --notransition --nogapped --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_lastz_identity98_exact50_ambiguous_iupac_notransition_nogapped_step10.dotplot

lastz T2T_chrX.fa T2T_chrY.fa --identity=97 --exact=50 --ambiguous=iupac --notransition --nogapped --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_lastz_identity97_exact50_ambiguous_iupac_notransition_nogapped_step10.dotplot



lastz T2T_chrX.fa T2T_chrY.fa --exact=50 --ambiguous=iupac --notransition --nogapped --step=10 ‑‑format=general > T2T_chrX_chrY_lastz_exact50_ambiguous_iupac_notransition_nogapped_step10_general.txt

lastz T2T_chrX.fa T2T_chrY.fa --exact=50 --ambiguous=iupac --notransition --step=10 ‑‑format=general > T2T_chrX_chrY_lastz_exact50_ambiguous_iupac_notransition_step10_general.txt
lastz T2T_chrX.fa T2T_chrY.fa --exact=50 --ambiguous=iupac --step=10 ‑‑format=general > T2T_chrX_chrY_lastz_exact50_ambiguous_iupac_step10_general.txt
lastz T2T_chrX.fa T2T_chrY.fa --ambiguous=iupac ‑‑format=general > T2T_chrX_chrY_lastz_ambiguous_iupac_general.txt

lastz T2T_chrX.fa T2T_chrY.fa --identity=100 --exact=50 --ambiguous=iupac --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_lastz_identity100_exact50_ambiguous_iupac_step10.dotplot
lastz T2T_chrX.fa T2T_chrY.fa --identity=98 --exact=50 --ambiguous=iupac --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_identity98_lastz_exact50_ambiguous_iupac_step10.dotplot
lastz T2T_chrX.fa T2T_chrY.fa --identity=97 --exact=50 --ambiguous=iupac --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_identity97_lastz_exact50_ambiguous_iupac_step10.dotplot


lastz T2T_chrX.fa T2T_chrY.fa --exact=50 --ambiguous=iupac --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_lastz_exact50_ambiguous_iupac_step10.dotplot

lastz T2T_chrX.fa T2T_chrY.fa --ambiguous=iupac ‑‑format=general > T2T_chrX_chrY_lastz_ambiguous_iupac_general.txt

lastz T2T_chrX.fa T2T_chrY.fa --identity=100 --ambiguous=iupac ‑‑format=rdotplot > T2T_chrX_chrY_identity100_lastz_ambiguous_iupac.rdotplot
lastz T2T_chrX.fa T2T_chrY.fa --identity=100 --ambiguous=iupac ‑‑format=general > T2T_chrX_chrY_identity100_lastz_ambiguous_iupac_general.txt

lastz T2T_chrX.fa T2T_chrY.fa --identity=98 --ambiguous=iupac ‑‑format=rdotplot > T2T_chrX_chrY_identity98_lastz_ambiguous_iupac.rdotplot
lastz T2T_chrX.fa T2T_chrY.fa --identity=98 --ambiguous=iupac ‑‑format=general > T2T_chrX_chrY_identity98_lastz_ambiguous_iupac_general.txt


lastz T2T_chrX.fa T2T_chrY.fa --identity=100 ‑‑format=rdotplot > T2T_chrX_chrY_identity100_lastz.rdotplot
lastz T2T_chrX.fa T2T_chrY.fa --identity=100 ‑‑format=general > T2T_chrX_chrY_identity100_lastz_general.txt

#  ‑‑nogapped eliminates the computation of gapped alignments
# Using ‑‑notransition lowers seeding sensitivity and reduces runtime (by a factor of about 10 in this case)
# When removing these flags, we get alignments less than 100%

# So I dont think we should use these options (‑‑nogapped and ‑‑notransition)

# --exact= flag: Exact match extension (‑‑exact) simply extends the seed until a mismatch is found. If the resulting length is enough, the extended seed is kept as an HSP for further processing. Exact match extension is most useful when the target and query are expected to be very similar, e.g. when aligning short reads to a similar reference genome.
# --step= flag: A large step, e.g. step=100, could potentially speed up the inference process.

cd /data/CEM/shared/public_data/references/T2T_HG002/T2T_SCC
lastz ../NA24385_chrX.fasta ../NA24385_chrY.fasta --exact=50 --ambiguous=iupac --notransition --nogapped --step=10 ‑‑format=general > T2T_chrX_chrY_lastz_exact50_ambiguous_iupac_notransition_nogapped_step10_general.txt

# So the sequences have soft masking applied and I wonder if this is affecting the results
# make unmasked versions of X and Y
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' T2T_chrX.fa > T2T_chrX_all_upper.fa
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' T2T_chrY.fa > T2T_chrY_all_upper.fa

lastz T2T_chrX_all_upper.fa T2T_chrY_all_upper.fa --identity=100 --exact=50 --ambiguous=iupac --notransition --nogapped --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity100_exact50_ambiguous_iupac_notransition_nogapped_step10.dotplot
lastz T2T_chrX_all_upper.fa T2T_chrY_all_upper.fa --exact=50 --ambiguous=iupac --notransition --nogapped --step=10 ‑‑format=general > T2T_chrX_chrY_all_upper_lastz_exact50_ambiguous_iupac_notransition_nogapped_step10_general.txt
lastz T2T_chrX_all_upper.fa T2T_chrY_all_upper.fa --identity=98 --exact=50 --ambiguous=iupac --notransition --nogapped --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity98_exact50_ambiguous_iupac_notransition_nogapped_step10.dotplot
# same number of lines as with 100 percent

lastz T2T_chrX_all_upper.fa T2T_chrY_all_upper.fa --exact=50 --ambiguous=iupac --step=10 ‑‑format=general > T2T_chrX_chrY_all_upper_lastz_identity100_exact50_ambiguous_iupac_step10_general.txt

lastz T2T_chrX_all_upper.fa T2T_chrY_all_upper.fa --identity=100 --exact=50 --ambiguous=iupac --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity100_exact50_ambiguous_iupac_step10.dotplot

lastz T2T_chrX_all_upper.fa T2T_chrY_all_upper.fa --identity=98 --exact=50 --ambiguous=iupac --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity98_exact50_ambiguous_iupac_step10.dotplot

lastz T2T_chrX_all_upper.fa T2T_chrY_all_upper.fa --identity=97 --exact=50 --ambiguous=iupac --step=10 ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity97_exact50_ambiguous_iupac_step10.dotplot


lastz_32 ../T2T_chrX_all_upper.fa ../T2T_chrY_all_upper.fa --exact=50 --ambiguous=iupac --step=10 ‑‑notransition ‑‑format=general > lastz_32_T2T_chrX_chrY_all_upper_lastz_exact50_ambiguous_iupac_step10_notransition_general.txt
lastz_32 ../T2T_chrX_all_upper.fa ../T2T_chrY_all_upper.fa --exact=50 --ambiguous=iupac --step=10 ‑‑notransition ‑‑format=general --allocate:traceback=1.99G > lastz_32_T2T_chrX_chrY_all_upper_lastz_exact50_ambiguous_iupac_step10_notransition_general_traceback1.99G.txt


lastz_32 ../T2T_chrX.fa[unmask] ../T2T_chrY.fa[unmask] --step=10 --seed=match12 --notransition --exact=20 --match=1,5 --ambiguous=iupac --format=general > results_general.txt
FAILURE: in init_from_anchors(), structure size would exceed 2^32 (1972530576 + 128*246566321)
 consider raising scoring threshold (--hspthresh or --exact) or breaking your target sequence into smaller pieces
 lastz_32 ../T2T_chrX.fa[unmask] ../T2T_chrY.fa[unmask] --step=10 --seed=match12 --notransition --exact=50 --match=1,5 --ambiguous=iupac --format=general > results_general.txt


lastz_32 ../T2T_chrX.fa ../T2T_chrY.fa --step=10 --seed=match12 --notransition --exact=20 --match=1,5 --ambiguous=iupac --format=general > results_general_withmasking.txt
lastz_32 ../T2T_chrX.fa ../T2T_chrY.fa --step=10 --seed=match12 --notransition --exact=20 --match=1,5 --ambiguous=iupac --format=rdotplot --identity=100 > results_general_withmasking_identity100.dotplot
lastz_32 ../T2T_chrX.fa ../T2T_chrY.fa --step=10 --seed=match12 --notransition --exact=20 --match=1,5 --ambiguous=iupac --format=rdotplot --identity=98 > results_general_withmasking_identity98.dotplot
lastz_32 ../T2T_chrX.fa ../T2T_chrY.fa --step=10 --seed=match12 --notransition --exact=20 --match=1,5 --ambiguous=iupac --format=rdotplot --identity=97 > results_general_withmasking_identity97.dotplot


lastz ../T2T_chrX.fa[unmask] ../T2T_chrY.fa[unmask] --step=20 --seed=match12 --notransition --exact=20 --match=1,5 --ambiguous=iupac --format=rdotplot --identity=100 > results_step_20_seed_match12_notransition_exact20_match1_5_ambiguous_iupac_withmasking_identity100.dotplot

lastz ../T2T_chrX.fa[unmask] ../T2T_chrY.fa[unmask] --step=20 --seed=match12 --notransition --exact=50 --match=1,5 --ambiguous=iupac --format=rdotplot --identity=100 > results_step_20_seed_match12_notransition_exact50_match1_5_ambiguous_iupac_withmasking_identity100.dotplot

lastz ../T2T_chrX.fa[unmask] ../T2T_chrY.fa[unmask] --step=50 --seed=match12 --notransition --exact=50 --match=1,5 --ambiguous=iupac --format=rdotplot --identity=100 > results_step_50_seed_match12_notransition_exact50_match1_5_ambiguous_iupac_withmasking_identity100.dotplot

lastz_32 ../T2T_chrX.fa ../T2T_chrY.fa --step=10 --seed=match12 --exact=20 --match=1,5 --ambiguous=iupac --format=general > results_general_withmasking_02.txt


lastz ../T2T_chrX.fa[unmask] ../T2T_chrY.fa[unmask] --identity=98 --exact=50 --ambiguous=iupac --step=10 ‑‑notransition ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity98_exact50_ambiguous_iupac_step10_notransition.dotplot

lastz ../T2T_chrX.fa ../T2T_chrY.fa --identity=98 --exact=50 --ambiguous=iupac --step=10 ‑‑notransition ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity98_exact50_ambiguous_iupac_step10_notransition_softmasked.dotplot

lastz ../T2T_chrX.fa ../T2T_chrY.fa --identity=97 --exact=50 --ambiguous=iupac --step=10 ‑‑notransition ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity97_exact50_ambiguous_iupac_step10_notransition_softmasked.dotplot

lastz ../T2T_chrX.fa ../T2T_chrY.fa --identity=100 --exact=50 --ambiguous=iupac --step=10 ‑‑notransition ‑‑format=rdotplot > T2T_chrX_chrY_all_upper_lastz_identity100_exact50_ambiguous_iupac_step10_notransition_softmasked.dotplot

```

# Creation of Salmon indexes

Salmon version 1.8.0 used

Command: salmon index -t GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded.fa -i GCA_009914755.4_CHM13_T2T_v2.0_genomic_YHardMasked_ChrNamesAdded_salmon_index
