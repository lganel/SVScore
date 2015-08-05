# SVScore
Prioritize structural variants based on annotations and C scores

Still under active development

# Usage
```
usage: ./svscore.pl [-ds] [-g genefile] [-e exonfile] [-n exonannotationcolumn] [-C binarycaddfile] -c caddfile vcf
    -d	      Debug (verbose) mode, keeps intermediate and supporting files
    -s	      Create/download supporting files and quit
    -c	      Points to whole_genome_SNVs.tsv.gz
    -C	      Points to cadd_v1.2.idx (defaults to current directory if not provided; *.bin file must be in same directory as index)
    -g	      Used to point to gene BED file (refGene.genes.b37.bed)
    -e	      Used to point to exon BED file (refGene.exons.b37.bed)
    -n	      Column number for annotation in exon BED file to be added to VCF (5)

    --help    Display this message
    --version Display version
```

# Dependencies
* [vcfanno](https://www.github.com/brentp/vcfanno)
  * cadd_v1.3a.idx - click [here](https://s3.amazonaws.com/vcfanno/cadd_v1.3a.idx) to download
  * cadd_v1.3a.bin - click [here](https://s3.amazonaws.com/vcfanno/cadd_v1.3a.bin) to download

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


# Notes
SVScore creates a file of introns (unless a file called introns.bed already exists in the current directory) by subtracting the exon file from the gene file using bedtools. So, if there is already file called introns.bed in the current directory, rename it or SVScore will not work correctly.

SVScore will clobber any files in the current directory with the name {prefix}header, where {prefix} is the input filename without the .vcf suffix

Input VCF files may be gzipped, but gzipped files must end with .gz. Uncompressed input files should not end with this suffix.

For BND variants, primary mate is considered the left breakend and the secondary mate is considered the right breakend.

If only one mate line of a BND variant is present in the VCF file, left and right breakend scores are still calculated, as well as one truncation score if applicable (whether it is the left or right truncation score depends on whether the line describes a primary or secondary mate).

Currently, SVScore requires two files containing CADD scores - whole_genome_SNVs.tsv.gz (from the CADD link above) and cadd_v1.3a.bin (from the vcfanno link above). This may change in the future. Both files must be in the same directory as their respective index file, available for download via the same links. The files may have any name so long as these names are supplied using the -c and -C options and each file's name matches that of its index (up to the suffix).
