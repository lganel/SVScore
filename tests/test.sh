if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1
fi

echo -e "********Testing SVScore. Some warnings beginning with \"Warning: missing secondary multiline variant at ID:\" are normal********\n\n"

## Strip header from truth set
grep -v "^#" NA12878.sv.svscore.vcf > NA12878.sv.svscore.noheader.vcf
grep -v "^#" stresstest.svscore.vcf > stresstest.svscore.noheader.vcf

# Perform stress test
echo "********Stress Test********"
echo "Generating annotation files"
../generateannotations.pl -c 1 -a 2 -b 3 -t 4 -s 5 -e 6 -f 7 dummyannotations.bed
cd ..
./svscore.pl -f tests/dummyannotations.introns.bed -e tests/dummyannotations.exons.bed -o max,sum,top2,top2weighted,top3weighted,top4weighted,mean,meanweighted -dvc tests/dummyCADD.tsv.gz -i tests/stresstest.vcf | grep -v "^#" > tests/stresstest.svscore.test.vcf
cd tests
echo -e "\n\n********NA12878********"
echo "Generating annotation files"
../generateannotations.pl
cd ..
./svscore.pl -o max,sum,top5,top10,mean -dvc $1 -i tests/NA12878.sv.vcf | grep -v "^#" > tests/NA12878.sv.svscore.test.vcf
echo -e "\n\n"
cd tests

diff stresstest.svscore.test.vcf stresstest.svscore.noheader.vcf > stresstest.diff
diff NA12878.sv.svscore.test.vcf NA12878.sv.svscore.noheader.vcf > NA12878.diff

STRESSTESTOUTPUTLINES=`cat stresstest.svscore.test.vcf | wc -l`
STRESSTESTDIFFLINES=`cat stresstest.diff | wc -l`
NA12878OUTPUTLINES=`cat NA12878.sv.svscore.test.vcf | wc -l`
NA12878DIFFLINES=`cat NA12878.diff | wc -l`

if [ $STRESSTESTDIFFLINES = 0 ] && [ $STRESSTESTOUTPUTLINES -gt 0 ]; then
  echo "********STRESSTEST TEST PASSED********"
else
  echo "********STRESSTEST TEST FAILED********"
  if [ $STRESSTESTDIFFLINES -gt 0 ]; then
    echo "********STRESSTEST: DIFFLINES > 0********"
    cat stresstest.diff
  else
    echo "********STRESSTEST: NO OUTPUT********"
  fi
fi

if [ $NA12878DIFFLINES = 0 ] && [ $NA12878OUTPUTLINES -gt 0 ]; then
  echo "********NA12878 TEST PASSED********"
else
  echo "********NA12878 TEST FAILED********"
  if [ $NA12878DIFFLINES -gt 0 ]; then
    echo "********NA12878: DIFFLINES > 0********"
    cat NA12878.diff
  else
    echo "********NA12878: NO OUTPUT********"
  fi
fi

rm -rf ../svscoretmp/
rm -f NA12878.sv.svscore.noheader.vcf NA12878.sv.svscore.test.vcf NA12878.diff stresstest.svscore.noheader.vcf stresstest.svscore.test.vcf stresstest.diff dummyannotations.*.bed* refGene*
