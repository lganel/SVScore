#!/bin/bash

grep -v "^#" NA12878.sv.svscore.vcf > NA12878.sv.svscore.noheader.vcf

../svscore.pl -dc $1 NA12878.sv.vcf | grep -v "^#" > NA12878.sv.svscore.test.noheader.vcf

echo

OUTPUTLINES=`cat NA12878.sv.svscore.test.noheader.vcf | wc -l`
DIFFLINES=`diff NA12878.sv.svscore.test.noheader.vcf NA12878.sv.svscore.noheader.vcf | wc -l`

if [ $DIFFLINES = 0 ] && [ $OUTPUTLINES -gt 0 ]; then
  echo "TEST PASSED"
else
  echo "TEST FAILED"
fi

rm -rf svscoretmp/
rm -f refGene.exons.b37.bed refGene.genes.b37.bed introns.bed conf.toml NA12878.sv.svscore.noheader.vcf NA12878.sv.svscore.test.noheader.vcf
