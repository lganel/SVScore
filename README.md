# SVScore
Prioritize structural variants based on annotations and C scores

Still under active development

# Changes in version 0.2
* Outputs sorted VCF to standard output, no longer creates *.scored.txt file

* Relaxed stringency on SVTYPE field

* File prefix now disregards path prefix to always work in current directory

* Removed natural sort requirement. SVScore now automatically sorts input VCF based on ID

* Returns a negative score for variants in regions with no C scores (e.g. GL contigs)

* Returns a score of 100 for deletions/duplications above 1 Mbp

# Usage
```
usage: ./svscore.pl [-ds] [-c caddfile] [-g genefile] [-e exonfile] [-n exonannotationcolumn] vcf
    -d        Debug (verbose) mode, keeps intermediate and supporting files
    -s        Create/download supporting files and quit
    -c        Used to point to whole_genome_SNVs.tsv.gz
    -g        Used to point to gene BED file
    -e        Used to point to exon BED file
    -n        Column number for annotation in exon BED file to be added to VCF

    --help    Display this message
    --version Display version
```

# Dependencies
* [vcfanno](https://www.github.com/brentp/vcfanno)

* [bedtools](https://www.github.com/arq5x/bedtools)

* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from [CADD](http://cadd.gs.washington.edu/download) 

* Your favorite tab-delimited, exon- and gene-describing BED files (optional). If not supplied, svscore.pl will automatically download them (functionality courtesy of Colby Chiang).

  * If using your own exon file, you must use -n to specify which column contains the desired annotation (typically gene symbol or gene name).

  * Files must be sorted and contain no header.
  
  * Gene annotation files must have at least five columns, representing, in order:

    1. Chromosome number (with no "chr" prefix)

    2. Start coordinate
    
    3. Stop coordinate
    
    4. Strand
    
    5. Gene symbol/name/other identifier


# Notes
SVScore creates a file of introns (unless a file called introns.bed already exists in the current directory) by subtracting the exon file from the gene file using bedtools

SVScore will clobber any files in the current directory with the name {prefix}header, where {prefix} is the input filename without the .vcf suffix

Input VCF files may be gzipped, but gzipped files must end with .gz. Uncompressed input files should not end with this suffix.
