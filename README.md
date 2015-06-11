# SVScore
Prioritize structural variants based on annotations and C scores
Still under active development

# Usage
```
usage: ./svscore.pl [-dza] vcf
  -d	Debug (verbose) mode
  -z	Indicates that vcf is gzipped
  -a	Indicates that vcf has already been annotated using vcfanno
```

# Dependencies
* [vcfanno](https://www.github.com/brentp/vcfanno)

* whole_genome_SNVs.tsv.gz (and .tbi) - file of all possible hg19/GRCh37 SNVs and associated C scores from [CADD](http://cadd.gs.washington.edu/download) 

* Your favorite exon-describing BED file. If you don't have one, execute the following command:

```
curl  -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" | gzip -cdfq | awk '{gsub("^chr","",$3); n=int($9); split($10,start,",");split($11,end,","); for(i=1;i<=n;++i) {print $3,start[i],end[i],$2"."i,$13,$2; } }' OFS="\t"  | bedtools sort > refGene.exons.b37.bed
```
