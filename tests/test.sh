if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1
fi

echo -e "********Testing SVScore. Some warnings beginning with \"Warning: missing secondary multiline variant at ID:\" are normal********\n\n"

## Strip header from truth set (and strip vcfanno annotations from NA12878)
#grep -v "^#" tests/NA12878.sv.svscore.vcf | perl -ne 's/(left_|right_)?Intron[^;]*;//g; s/(left_|right_)?ExonTranscript=[^;]*;//g; print' > tests/NA12878.sv.svscore.noheader.vcf
#grep -v "^#" tests/stresstest.svscore.vcf > tests/stresstest.svscore.noheader.vcf

# Perform stress test
echo "********Stress Test********"
echo "Generating annotation files"
$SVSCOREDIR/generateannotations.pl -c 1 -a 2 -b 3 -t 4 -s 5 -e 6 -f 7  -o tests tests/dummyannotations.bed
$SVSCOREDIR/svscore.pl -f tests/dummyannotations.introns.bed -e tests/dummyannotations.exons.bed -o max,sum,top2,top2weighted,top3weighted,top4weighted,mean,meanweighted -dvc tests/dummyCADD.tsv.gz -i tests/stresstest.vcf | grep -v "^#" > tests/stresstest.svscore.test.vcf
echo -e "\n\n"

diff tests/stresstest.svscore.test.vcf tests/stresstest.svscore.vcf > tests/stresstest.diff

STRESSTESTOUTPUTLINES=`cat tests/stresstest.svscore.test.vcf | wc -l`
STRESSTESTDIFFLINES=`cat tests/stresstest.diff | wc -l`

if [ $STRESSTESTDIFFLINES = 0 ] && [ $STRESSTESTOUTPUTLINES -gt 0 ]; then
  echo "********STRESSTEST TEST PASSED********"
else
  echo "********STRESSTEST TEST FAILED********"
  if [ $STRESSTESTDIFFLINES -gt 0 ]; then
    echo "********STRESSTEST: DIFFLINES > 0********"
    cat tests/stresstest.diff
  else
    echo "********STRESSTEST: NO OUTPUT********"
  fi
fi

# NA12878
echo -e "********NA12878********"
echo "Generating annotation files"
$SVSCOREDIR/generateannotations.pl -o tests
$SVSCOREDIR/svscore.pl -o max,sum,top5,top10,mean -e tests/refGene.exons.bed -f tests/refGene.introns.bed -dvc $1 -i tests/NA12878.sv.vcf | grep -v "^#" | perl -ne 's/(left_|right_)?Intron[^;]*;//g; s/(left_|right_)?ExonTranscript=[^;]*;//g; print' > tests/NA12878.sv.svscore.test.vcf
echo -e "\n\n"

diff tests/NA12878.sv.svscore.test.vcf tests/NA12878.sv.svscore.vcf > tests/NA12878.diff

NA12878OUTPUTLINES=`cat tests/NA12878.sv.svscore.test.vcf | wc -l`
NA12878DIFFLINES=`cat tests/NA12878.diff | wc -l`

if [ $NA12878DIFFLINES = 0 ] && [ $NA12878OUTPUTLINES -gt 0 ]; then
  echo "********NA12878 TEST PASSED********"
else
  echo "********NA12878 TEST FAILED********"
  if [ $NA12878DIFFLINES -gt 0 ]; then
    echo "********NA12878: DIFFLINES > 0********"
    cat tests/NA12878.diff
  else
    echo "********NA12878: NO OUTPUT********"
  fi
fi

rm -rf svscoretmp/
rm -f tests/NA12878.sv.svscore.noheader.vcf tests/NA12878.sv.svscore.test.vcf tests/NA12878.diff tests/stresstest.svscore.noheader.vcf tests/stresstest.svscore.test.vcf tests/stresstest.diff tests/dummyannotations.*.bed* tests/refGene* conf.toml
