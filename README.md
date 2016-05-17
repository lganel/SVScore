# SVScore
SVScore is a VCF annotation tool which scores structural variants by predicted pathogenicity based on SNP-based CADD scores. For each variant, SVScore first defines important genomic intervals based on the variant type, breakend confidence intervals, and nearby gene/exon annotations. It then applies an operation to each interval to aggregate the CADD scores in that interval into an interval score. A score for a given operation defined as the maximum of all interval scores calculated using that operation.

## Usage
```
usage: ./svscore.pl [-dv] [-o op] [-t topnumber] [-g genefile] [-m geneannotationcolumn] [-p genestrandcolumn] [-e exonfile] [-c caddfile] -i vcf
    -i        Input VCF file. May be bgzip compressed (ending in .vcf.gz). Use "-i stdin" if using standard input
    -d        Debug mode, keeps intermediate and supporting files, displays progress
    -v        Verbose mode - show all calculated scores (left/right/span/ltrunc/rtrunc, as appropriate)
    -o        Comma-separated list of operations to perform on CADD score intervals (must be some combination of sum, max, mean, meanweighted, top\\d, and top\\dweighted - defaults to top10weighted)
    -g        Points to gene BED file (refGene.genes.b37.bed)
    -e        Points to exon BED file (refGene.exons.b37.bed)
    -m        Column number for gene name in gene BED file (4)
    -p        Column number for strand in gene BED file (5)
    -n        Column number for gene name in exon BED file (5 for refGene.exons.b37.bed, 4 otherwise)
    -c        Points to whole_genome_SNVs.tsv.gz (defaults to current directory)

    --help    Display this message
    --version Display version
```

## Output
SVScore outputs a VCF file with scores added to the INFO field of each variant. The VCF header is also updated to include those scores which are added. Each score field has the following format: SVSCORE\[op\](_[interval]), where [op] represents the operation used to calculate that score (see [Operations](#operations)) and [interval] represents the interval over which the score was calculated, which is one of left breakend, right breakend, span (for DEL/DUP), left truncation score (for INV/TRX variants which seem to truncate a gene on the left side, the interval is from the left breakend to the end of the gene), and right truncation score. Scores with no interval listed (such as SVSCOREMAX=) are the maximum over all intervals for that operation.

## Intervals
For each variant, scores are calculated over a number of intervals which varies by SV type. The intervals chosen for each SV type, are described in [Supported SV types and intervals](#supported-sv-types-and-intervals)
* LEFT - confidence interval around the left breakpoint
* RIGHT - confidence interval around the right breakpoint
* SPAN - from the most likely base in the left breakpoint confidence interval to the most likely base in the right breakpoint confidence interval
* LTRUNC - left truncation
* RTRUNC - right truncation

Truncation intervals are defined for each gene which seems to be truncated by a variant. The interval extends from the furthest upstream base of the furthest upstream breakend (LEFT for genes on the + strand, RIGHT for those on the - strand) to the end of the gene. Each truncation score is the maximum over all genes truncated by a variant.

## Supported SV types and intervals
|      | LEFT | RIGHT | SPAN | LTRUNC | RTRUNC | Notes
|:---:|:---:|:---:|:---:|:---:|:---:|:---
|DEL|X|X|X|X|X|
|DUP|X|X|X|||
|INV|X|X||X|X|
|BND|X|X||||
|TRX|X|X||X|X|Translocation
|INS|X|X||X|X|
|CNV|X|X|X|||
|MEI|X|X||||
To function correctly, SVScore requires that POS=END and CIPOS=CIEND for INS variants

LTRUNC and RTRUNC scores are only calculated when a breakend overlaps an exon or a breakend overlaps an intron which is not also touched by the opposite breakend.

## Operations
-o specifies the operation(s) used to calculate SVScores. These operations are applied to each interval of the SV (see [Supported SV types and intervals](#supported-sv-types-and-intervals)). This option takes an arbitrary-length, case insensitive, comma-separated list of operations from the following list:
* sum - reports the sum of each interval
* max - reports the maximum of each interval
* mean - reports the arithmetic mean of each interval
* meanweighted - reports the arithmetic mean of each interval with the left and right breakends weighted by the probability distribution in the PRPOS and PREND fields of the VCF INFO column 
* top[n] - reports the mean of the [n] largest scores in each interval. If the interval has less than [n] possible breakpoints, then the top[n] score of an interval is equal to the mean score
* top[n]weighted - reports the mean of the [n] largest scores in each interval with the left and right breakends weighted by the probability distribution in the PRPOS and PREND fields of the VCF INFO column

For weighted operations, if PRPOS is not found in the header, SVScore will calculate unweighted means with a warning. If PRPOS or PREND is missing from a variant but is present in the header, that variant will receive a score of -1 for all weighted operations

## Dependencies
* A Linux-like system with a Bash-like shell
* Perl v5.10.1 or later
* [vcfanno](https://www.github.com/brentp/vcfanno) v0.0.10
* [bedtools](https://www.github.com/arq5x/bedtools2)
* [svtools](https://github.com/hall-lab/svtools)
* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from v1.3 of [CADD](http://cadd.gs.washington.edu/download) 
* [tabix/bgzip](https://github.com/samtools/htslib)
* Your favorite hg19/GRCh37-based, tab-delimited, exon- and gene-describing BED files (optional). If not supplied, svscore.pl will automatically download RefSeq annotations (functionality courtesy of Colby Chiang).
  * SVScore expects custom gene annotation files to contain gene symbol/name in column 4 and strand information in column 5, though these can be changed with -m and -p
  * SVScore expects custom exon annotation files to contain gene symbol/name in column 4, though this can be changed with -n
The following must be in your path to use SVScore: svtools, vcfanno, bedtools, tabix
  
## Notes
If an input VCF file already has SVSCORE annotations in the INFO column, new annotations will overwrite old ones.

Input VCF files may be gzipped, but gzipped files must end with .gz. Uncompressed input files should not end with this suffix. Annotation files may be gzipped or unzipped. SVScore will zip/unzip files as necessary using bgzip and zcat.

For multiline variants, primary mate is considered the left breakend and the secondary mate is considered the right breakend.

If only one mate line of a multiline variant is present in the VCF file, left and right breakend scores are still calculated, as well as one truncation score if applicable (whether it is the left or right truncation score depends on whether the line describes a primary or secondary mate). There must be a CIEND interval in the INFO field for this to happen.

Variants with type DEL, DUP, or CNV which are over 1 Mb in length are automatically given a score of 100
