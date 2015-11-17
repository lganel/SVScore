#!/bin/bash

if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1;
fi

grep -v "^#" NA12878.sv.svscore.vcf > NA12878.sv.svscore.noheader.vcf

../svscore.pl -dvc $1 NA12878.sv.vcf | grep -v "^#" > NA12878.sv.svscore.test.noheader.vcf

echo

OUTPUTLINES=`cat NA12878.sv.svscore.test.noheader.vcf | wc -l`
DIFFLINES=`diff NA12878.sv.svscore.test.noheader.vcf NA12878.sv.svscore.noheader.vcf | wc -l`

if [ $DIFFLINES = 0 ] && [ $OUTPUTLINES -gt 0 ]; then
  echo "TEST PASSED"
else
  echo "TEST FAILED"
fi

rm -rf svscoretmp/
rm -f refGene.exons.b37.bed refGene.genes.b37.bed* introns.refGene.genes.b37.bed.gz.refGene.exons.b37.bed.bed.gz* conf.toml NA12878.sv.svscore.noheader.vcf NA12878.sv.svscore.test.noheader.vcf
