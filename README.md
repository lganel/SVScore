# SVScore
Prioritize structural variants based on annotations and C scores

# Usage
```
usage: ./svscore.pl [-ds] [-g genefile] [-m geneannotationcolumn] [-e exonfile] [-n exonannotationcolumn] [-c caddfile] vcf
    -d	      Debug (verbose) mode, keeps intermediate and supporting files
    -s	      Create/download supporting files and quit
    -c	      Points to whole_genome_SNVs.tsv.gz (defaults to current directory)
    -g	      Points to gene BED file (refGene.genes.b37.bed)
    -m	      Column number for annotation in gene BED file to be added to VCF (4)
    -e	      Points to exon BED file (refGene.exons.b37.bed)
    -n	      Column number for annotation in exon BED file to be added to VCF (4 if using -e, 5 otherwise)
    -w	      Weight CADD scores in breakends by probability distribution (requires PRPOS/PREND in INFO field)
    -o	      Specify operation to perform on CADD scores (must be sum, max, top[number], or all - defaults to all)

    --help    Display this message
    --version Display version
```

-o specifies the operation used to calculate SVScores. "sum" and "max" report the sum and maximum of each interval respectively, while "top[number]" reports the sum of the top [number] scores in each interval and "all" reports the maximum score, sum of all scores, and sum of the top 100 scores for each interval. This option is case insensitive. Under -o top[number] or -o all, if the interval is small than [number] positions (100 in the case of all), then the top[number] score is equal to the sum score.

# Dependencies
* A Linux-like system with a Bash-like shell

* Perl v5.8.7 or later

* [vcfanno](https://www.github.com/brentp/vcfanno) v0.0.7

* [bedtools](https://www.github.com/arq5x/bedtools2)

* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from v1.3 of [CADD](http://cadd.gs.washington.edu/download) 

* [tabix](https://github.com/samtools/htslib)

* Your favorite hg19/GRCh37-based, tab-delimited, exon- and gene-describing BED files (optional). If not supplied, svscore.pl will automatically download them (functionality courtesy of Colby Chiang).

  * If using your own exon file, you must use -n to specify which column contains the desired annotation (typically gene symbol or gene name).

  * SVScore expects custom gene annotation files to contain gene symbol/gene name in column 4 and strand information in column 5
  
# Notes
SVScore outputs a VCF file with the following scores added to the INFO field of each variant. The VCF header is also updated to include those scores which are added.
  * Under -o max or -o all
      1. SVSCOREMAX_LEFT for all variants - max C score within the left breakend
      2. SVSCOREMAX_RIGHT for all variants - max C score within the right breakend
      3. SVSCOREMAX_SPAN for DEL/DUP - max C score within the span
      4. SVSCOREMAX_LTRUNC for BND/INV variants whose left breakend seems to truncate a gene - max C score within gene downstream of beginning of left breakend 
      5. SVSCOREMAX_RTRUNC for BND/INV variants whose right breakend seems to truncate a gene - max C score within gene downstream of beginning of right breakend 
  * Under -o sum or -o all
      1. SVSCORESUM_LEFT for all variants - sum of C scores within the left breakend
      2. SVSCORESUM_RIGHT for all variants - sum of C scores within the right breakend
      3. SVSCORESUM_SPAN for DEL/DUP - sum of C scores within the span
      4. SVSCORESUM_LTRUNC for BND/INV variants whose left breakend seems to truncate a gene - sum of C scores within gene downstream of beginning of left breakend 
      5. SVSCORESUM_RTRUNC for BND/INV variants whose right breakend seems to truncate a gene - sum of C scores within gene downstream of beginning of right breakend 

SVScore creates a file of introns (unless a file called introns.bed already exists in the current directory) by subtracting the exon file from the gene file using bedtools. So, if there is already file called introns.bed in the current directory, rename it or SVScore will not work correctly.

Input VCF files may be gzipped, but gzipped files must end with .gz. Uncompressed input files should not end with this suffix.

For BND variants, primary mate is considered the left breakend and the secondary mate is considered the right breakend.

If only one mate line of a BND variant is present in the VCF file, left and right breakend scores are still calculated, as well as one truncation score if applicable (whether it is the left or right truncation score depends on whether the line describes a primary or secondary mate). There must be a CIEND interval in the INFO field for this to happen.

-w is ignored with a warning if PRPOS is not represented in the header
