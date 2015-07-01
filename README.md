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
usage: ./svscore.pl [-ds] [-c caddfile] [-g genefile] [-e exonfile] vcf
    -d        Debug (verbose) mode, keeps intermediate and supporting files
    -s        Create/download supporting files and quit
    -c        Used to point to whole_genome_SNVs.tsv.gz
    -g	      Used to point to gene annotation file
    -e	      Used to point to exon annotation file
    --help    Display this message
    --version Display version
```

# Dependencies
* [vcfanno](https://www.github.com/brentp/vcfanno)

* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from [CADD](http://cadd.gs.washington.edu/download) 

* Your favorite exon- and gene-describing BED files (optional). If you don't have any, svscore.pl will automatically download them (functionality courtesy of Colby Chiang). If using your own files, they must be sorted and contain the gene/exon name in column 5

* [bedtools](https://www.github.com/arq5x/bedtools)


# Notes
SVScore creates a file of introns (unless a file called introns.bed already exists in the current directory) by subtracting the exon file from the gene file using bedtools

SVScore will clobber any files in the current directory with the name {prefix}header, where {prefix} is the input filename without the .vcf suffix

Input VCF files may be gzipped, but gzipped files must end with .gz. Uncompressed input files should not end with this suffix.
