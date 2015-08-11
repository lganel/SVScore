# SVScore
Prioritize structural variants based on annotations and C scores

Still under active development

# Usage
```
usage: ./svscore.pl [-ds] [-g genefile] [-e exonfile] [-n exonannotationcolumn] [-c caddfile] vcf
    -d	      Debug (verbose) mode, keeps intermediate and supporting files
    -s	      Create/download supporting files and quit
    -c	      Points to whole_genome_SNVs.tsv.gz (defaults to current directory)
    -g	      Used to point to gene BED file (refGene.genes.b37.bed)
    -e	      Used to point to exon BED file (refGene.exons.b37.bed)
    -n	      Column number for annotation in exon BED file to be added to VCF (5)
    -w	      Weight CADD scores in breakends by probability distribution (requires PRPOS/PREND in INFO field)
    -o	      Specify operation to perform on CADD scores (must be sum, max, or both - defaults to both)

    --help    Display this message
    --version Display version
```

# Dependencies
* A Bash-like shell

* Perl v5.8.7 or later

* [vcfanno](https://www.github.com/brentp/vcfanno)

* [bedtools](https://www.github.com/arq5x/bedtools2)

* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from [CADD](http://cadd.gs.washington.edu/download) 

* [tabix](https://github.com/samtools/htslib)

* Your favorite hg19/GRCh37-based, tab-delimited, exon- and gene-describing BED files (optional). If not supplied, svscore.pl will automatically download them (functionality courtesy of Colby Chiang).

  * If using your own exon file, you must use -n to specify which column contains the desired annotation (typically gene symbol or gene name).

  * Files must be sorted and contain no header.
  
  * Gene annotation files must have at least five columns, representing, in order:
    1. Chromosome number (with no "chr" prefix)
    2. Start coordinate
    3. Stop coordinate
    4. Strand
    5. Gene symbol/name/other identifier


# Output
SVScore outputs a VCF file with the following scores added to the INFO field of each variant. The VCF header is also updated to include those scores which are added.
  * Under -o max or -o both
      1. SVSCOREMAX_LEFT for all variants - max C score within the left breakend
      2. SVSCOREMAX_RIGHT for all variants - max C score within the right breakend
      3. SVSCOREMAX_SPAN for DEL/DUP - max C score within the span
      4. SVSCOREMAX_LTRUNC for BND/INV variants whose left breakend seems to truncate a gene - max C score within gene downstream of beginning of left breakend 
      5. SVSCOREMAX_RTRUNC for BND/INV variants whose right breakend seems to truncate a gene - max C score within gene downstream of beginning of right breakend 
  * Under -o sum or -o both
      1. SVSCORESUM_LEFT for all variants - sum of C scores within the left breakend
      2. SVSCORESUM_RIGHT for all variants - sum of C scores within the right breakend
      3. SVSCORESUM_SPAN for DEL/DUP - sum of C scores within the span
      4. SVSCORESUM_LTRUNC for BND/INV variants whose left breakend seems to truncate a gene - sum of C scores within gene downstream of beginning of left breakend 
      5. SVSCORESUM_RTRUNC for BND/INV variants whose right breakend seems to truncate a gene - sum of C scores within gene downstream of beginning of right breakend 

# Notes
SVScore creates a file of introns (unless a file called introns.bed already exists in the current directory) by subtracting the exon file from the gene file using bedtools. So, if there is already file called introns.bed in the current directory, rename it or SVScore will not work correctly.

SVScore will clobber any files in the current directory with the name {prefix}header, where {prefix} is the input filename without the .vcf suffix

Input VCF files may be gzipped, but gzipped files must end with .gz. Uncompressed input files should not end with this suffix.

For BND variants, primary mate is considered the left breakend and the secondary mate is considered the right breakend.

If only one mate line of a BND variant is present in the VCF file, left and right breakend scores are still calculated, as well as one truncation score if applicable (whether it is the left or right truncation score depends on whether the line describes a primary or secondary mate). There must be a CIEND interval in the INFO field for this to happen.

-w is ignored with a warning if PRPOS is not represented in the header
