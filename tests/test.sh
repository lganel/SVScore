#!/bin/bash

if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1;
fi

grep -v "^#" NA12878.sv.svscore.vcf > NA12878.sv.svscore.noheader.vcf

cd ..
./svscore.pl -dvc tests/$1 tests/NA12878.sv.vcf > tests/NA12878.sv.svscore.test.vcf
grep -v "^#" tests/NA12878.sv.svscore.test.vcf > tests/NA12878.sv.svscore.test.noheader.vcf
cd tests

echo

OUTPUTLINES=`cat NA12878.sv.svscore.test.noheader.vcf | wc -l`
DIFFLINES=`diff NA12878.sv.svscore.test.noheader.vcf NA12878.sv.svscore.noheader.vcf | wc -l`

if [ $DIFFLINES = 0 ] && [ $OUTPUTLINES -gt 0 ]; then
  echo "TEST PASSED"
else
  echo "TEST FAILED"
fi

rm -rf svscoretmp/
rm -f NA12878.sv.svscore.noheader.vcf NA12878.sv.svscore.test.noheader.vcf
