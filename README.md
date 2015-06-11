# SVScore
Prioritize structural variants based on annotations and C scores
Still under active development

# Usage
```
usage: ./svscore.pl [-dza] vcf
  -d	Debug (verbose) mode, keeps intermediate and supporting files
  -z	Indicates that vcf is gzipped
  -a	Indicates that vcf has already been annotated using vcfanno
  -s	Create/download supporting files and quit
```

# Dependencies
* [vcfanno](https://www.github.com/brentp/vcfanno)

* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from [CADD](http://cadd.gs.washington.edu/download) 

* Your favorite exon- and gene-describing BED files (optional). If you don't have any, svscore.pl will automatically download them (functionality courtesy of Colby Chiang)
