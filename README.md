# SVScore
SVScore is a VCF annotation tool which scores structural variants by predicted pathogenicity based on SNP-based CADD scores. For each variant, SVScore first defines important genomic intervals based on the variant type, breakend confidence intervals, and overlapping exon/intron annotations. It then applies an operation to each interval to aggregate the CADD scores in that interval into an interval score. A score for a given operation defined as the maximum of all interval scores calculated using that operation. SVScore is based on hg19/GRCh37.

## Usage
```
usage: ./svscore.pl [-dv] [-o op] [-e exonfile] [-f intronfile] [-c caddfile] -i vcf
    -i        Input VCF file. May be bgzip compressed (ending in .vcf.gz). Use \"-i stdin\" if using standard input
    -d        Debug mode, keeps intermediate and supporting files, displays progress
    -v        Verbose mode - show all calculated scores (left/right/span/ltrunc/rtrunc, as appropriate)
    -o        Comma-separated list of operations to perform on CADD score intervals (must be some combination of sum, max, mean, meanweighted, top\\d, and top\\dweighted - defaults to top10weighted)
    -e        Points to exon BED file (refGene.exons.bed)
    -f        Points to intron BED file (refGene.introns.bed)
    -c        Points to whole_genome_SNVs.tsv.gz (defaults to current directory)

    --help    Display this message
    --version Display version
```

## First Time Setup
After downloading SVScore, there are a few steps to follow before it is ready to use.
  1. Place the path to the SVScore installation directory in an environment variable called $SVSCOREDIR. Ideally, this should be done by appending the following line to the user's ~/.bashrc file or equivalent:
  ```
    export SVSCOREDIR=/path/to/installation/
  ```
  2. Test SVScore using `sh tests/test.sh path/to/whole_genome_SNVs.tsv.gz`
  3. Generate annotation files. For more on this, see [Annotation Files](#annotation-files).
  4. SVScore assumes the user's version of perl is installed in the default directory (`/usr/bin/perl`). If this is not the case, the first line of all .pl files should be changed to reflect the correct perl installation directory.

## Annotation Files
  * If planning to use the default (refGene) annotations, simply execute `./generateannotations.pl` to generate the annotation files required by SVScore. If planning to use a custom annotation track, `generateannotations.pl` can be used to generate custom annotation files, or the user can generate them manually (though this is not recommended).
  * If generating custom annotation files using `generateannotations.pl`, the user must supply an annotation track in which each line represents a transcript and contains the following columns: chromosome, transcript start position, transcript stop position, transcript strand, transcript name, exon start positions (comma-delimited), and exon stop positions (comma-delimited). Command line options must be used to specify each column number. To see usage instructions, execute `./generateannotations.pl --help`. `generateannotations.pl` will create two files in the current directory (or the directory specified using -o), named based on the prefix to the input file - [prefix].introns.bed and [prefix].exons.bed. These should be specified to SVScore using the -e and -f options.
  * **If generating custom annotation files for SVScore manually, users should ensure that each transcript has a unique name.** Annotation files should contain the following columns, in order:
    * Exon file:
      * 1 - Exon chromosome 
      * 2 - Exon start position
      * 3 - Exon stop position
      * 4 - Transcript name
      * 5 - Transcript start position
      * 6 - Transcript stop position
      * 7 - Transcript strand
    * Intron file:
      * 1 - Intron chromosome
      * 2 - Intron start position
      * 3 - Intron stop position
      * 4 - Transcript name
      * 5 - Intron number (arbitrary, but must be unique. Line number works well)

## Output
SVScore outputs a VCF file with scores added to the INFO field of each variant. The VCF header is also updated to include those scores which are added. Each score field has the following format: SVSCORE\[op\](_[interval]), where [op] represents the operation used to calculate that score (see [Operations](#operations)) and [interval] represents the interval over which the score was calculated, which is one of left breakend, right breakend, span (for DEL/DUP), left truncation score (for INV/DEL/INS variants which seem to truncate a transcript on the left side, the interval is from the most likely base of the left breakend to the end of the transcript), and right truncation score. Scores with no interval listed (such as SVSCOREMAX=) are the maximum over all intervals for that operation.

## Intervals
For each variant, scores are calculated over a number of intervals which varies by SV type. The intervals chosen for each SV type, are described in [Supported SV types and intervals](#supported-sv-types-and-intervals)
* LEFT - confidence interval around the left breakpoint
* RIGHT - confidence interval around the right breakpoint
* SPAN - from the most likely base in the left breakpoint confidence interval to the most likely base in the right breakpoint confidence interval
* LTRUNC - left truncation
* RTRUNC - right truncation

Truncation intervals are defined for each transcript which seems to be truncated by a variant. The interval extends from the most likely base of the furthest upstream breakend (LEFT for transcripts on the + strand, RIGHT for those on the - strand) to the end of the transcript. Each truncation score is the maximum over all transcripts truncated by a variant.

## Supported SV types and intervals
|      | LEFT | RIGHT | SPAN | LTRUNC | RTRUNC |
|:----:|:----:|:-----:|:----:|:------:|:------:|
| DEL  |  X   |   X   |  X   |    X   |    X   |
| DUP  |  X   |   X   |  X   |        |        |
| INV  |  X   |   X   |      |    X   |    X   |
| BND  |  X   |   X   |      |        |        |
| INS  |  X   |   X   |      |    X   |    X   |
| CNV  |  X   |   X   |  X   |        |        |
| MEI  |  X   |   X   |      |        |        |

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

For SPAN/LTRUNC/RTRUNC, these operations are applied to the scores of the bases in the interval. For LEFT/RIGHT intervals, the operations are applied to scores assigned to each possible breakpoint, which is calculated by taking the average of the 2 flanking bases (one on either side of the possible breakpoint)

## Dependencies
* A Linux-like system with a Bash-like shell, with sh, cd, sort (**With version sort, i.e. a -V option**), cat, zcat, rm, rmdir, diff, echo, chmod, grep, and awk
* Perl v5.10.1 or later
* [vcfanno](https://www.github.com/brentp/vcfanno) v0.0.11
* [svtools](https://github.com/hall-lab/svtools/releases/latest) v0.2.0 or later
* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from v1.3 of [CADD](http://cadd.gs.washington.edu/download) 
* [tabix/bgzip](https://github.com/samtools/htslib)
* (optional) An custom hg19/GRCh37-based gene/exon track. If one is not supplied, SVScore will download refGene annotations in `generateannotations.pl`. For more information, see [First Time Setup](#first-time-setup)

The following must be in your path to use SVScore: svtools, vcfanno, tabix
  
## Notes
If an input VCF file already has SVSCORE annotations in the INFO column, new annotations will overwrite old ones.

Input VCF files may be gzipped, but gzipped files must end with .gz. Uncompressed input files should not end with this suffix. Annotation files may be gzipped or unzipped. SVScore will zip/unzip files as necessary using bgzip and zcat.

For multiline variants, primary mate is considered the left breakend and the secondary mate is considered the right breakend.

If only one mate line of a multiline variant is present in the VCF file, left and right breakend scores are still calculated, as well as one truncation score if applicable (whether it is the left or right truncation score depends on whether the line describes a primary or secondary mate). There must be a CIEND interval in the INFO field for this to happen.

Variants with type DEL, DUP, or CNV which are over 1 Mb in length are automatically given a score of 100
