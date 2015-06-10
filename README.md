# SVScore
Prioritize structural variants based on annotations and C scores
Still under active development

# Usage
usage: ./svscore [-dza] vcf
  -d	Debug (verbose) mode
  -z	Indicates that vcf is gzipped
  -a	Indicates that vcf has already been annotated using vcfanno

# Dependencies
[vcfanno](https://www.github.com/brentp/vcfanno)
whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from [CADD](http://cadd.gs.washington.edu/download) 
