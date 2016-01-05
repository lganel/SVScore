# SVScore
Prioritize structural variants based on CADD scores

# Usage
```
usage: ./svscore.pl [-dsvw] [-o op] [-g genefile] [-m geneannotationcolumn] [-p genestrandcolumn] [-e exonfile] [-c caddfile] vcf
    -d	      Debug mode, keeps intermediate and supporting files, displays progress
    -s	      Create/download supporting files and quit
    -v	      Verbose mode - show all calculated scores (left/right/span/ltrunc/rtrunc, as appropriate)
    -w	      Weight CADD scores in breakends by probability distribution (requires PRPOS/PREND in INFO field)
    -o	      Specify operation to perform on CADD scores (must be sum, max, top[number], or all - defaults to all)
    -g	      Points to gene BED file (refGene.genes.b37.bed)
    -m	      Column number for annotation in gene BED file to be added to VCF (4)
    -p	      Column number for strand in gene BED file (5)
    -e	      Points to exon BED file (refGene.exons.b37.bed)
    -c	      Points to whole_genome_SNVs.tsv.gz (defaults to current directory)

    --help    Display this message
    --version Display version
```

-o specifies the operation used to calculate SVScores. "sum" and "max" report the sum and maximum of each interval respectively, while "top" reports the sum of the top 100 scores in each interval (use -t to use a different number of scores) and "all" reports the maximum score, sum of all scores, and sum of the top 100 scores for each interval. This option is case insensitive. Under -o top[number] or -o all, if the interval is smaller than [number] positions, then the top[number] score is equal to the sum score.

# Dependencies
* A Linux-like system with a Bash-like shell

* Perl v5.10.1 or later

* [vcfanno](https://www.github.com/brentp/vcfanno) v0.0.8

* [bedtools](https://www.github.com/arq5x/bedtools2)

* [svtools](https://github.com/hall-lab/svtools)

* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from v1.3 of [CADD](http://cadd.gs.washington.edu/download) 

* [tabix/bgzip](https://github.com/samtools/htslib)

* Your favorite hg19/GRCh37-based, tab-delimited, exon- and gene-describing BED files (optional). If not supplied, svscore.pl will automatically download RefSeq annotations (functionality courtesy of Colby Chiang).

  * SVScore expects custom gene annotation files to contain gene symbol/gene name in column 4 and strand information in column 5 (though these can be changed with -m and -p)
  
# Notes
The following must be in your path to use SVScore: svtools/bin, vcfanno, bedtools, tabix

SVScore outputs a VCF file with scores added to the INFO field of each variant. The VCF header is also updated to include those scores which are added. Each score field has the following format: SVSCORE{$op}(_{$interval}). $op represents the operation used to calculate that score, which is one of max, sum, top (mean of top n scores in the interval), or mean. $interval represents the interval over which the score was calculated, which is one of left breakend, right breakend, span (for DEL/DUP), left truncation score (for INV/TRX variants which seem to truncate a gene on the left side, the interval is from the left breakend to the end of the gene), and right truncation score. Scores with no interval listed are the maximum over all intervals for that operation.

Input VCF files may be gzipped, but gzipped files must end with .gz. Uncompressed input files should not end with this suffix. Annotation files may be gzipped or unzipped. SVScore will zip/unzip files as necessary using bgzip and zcat.

For BND variants, primary mate is considered the left breakend and the secondary mate is considered the right breakend.

If only one mate line of a BND variant is present in the VCF file, left and right breakend scores are still calculated, as well as one truncation score if applicable (whether it is the left or right truncation score depends on whether the line describes a primary or secondary mate). There must be a CIEND interval in the INFO field for this to happen.

-w is ignored with a warning if PRPOS is not represented in the header

Variants with type DEL, DUP, MEI, or CNV which are over 1 Mb in length are automatically given a score of 100
