if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1
fi

echo -e "********Testing SVScore. Some warnings beginning with \"Warning: missing secondary multiline variant at ID:\" are normal********\n\n"

## Strip header from truth set (and strip vcfanno annotations from NA12878)
grep -v "^#" NA12878.sv.svscore.vcf | perl -ne 's/(left_|right_)?Intron[^;]*;//g; s/(left_|right_)?ExonTranscript=[^;]*;//g; print' >  NA12878.sv.svscore.noheader.vcf
grep -v "^#" stresstest.svscore.vcf > stresstest.svscore.noheader.vcf

# Perform stress test
echo "********Stress Test********"
echo "Generating annotation files"
../generateannotations.pl -c 1 -a 2 -b 3 -t 4 -s 5 -e 6 -f 7 dummyannotations.bed
../svscore.pl -f dummyannotations.introns.bed -e dummyannotations.exons.bed -h .. -o max,sum,top2,top2weighted,top3weighted,top4weighted,mean,meanweighted -dvc dummyCADD.tsv.gz -i stresstest.vcf | grep -v "^#" > stresstest.svscore.test.vcf
echo -e "\n\n"

diff stresstest.svscore.test.vcf stresstest.svscore.noheader.vcf > stresstest.diff

STRESSTESTOUTPUTLINES=`cat stresstest.svscore.test.vcf | wc -l`
STRESSTESTDIFFLINES=`cat stresstest.diff | wc -l`

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

# NA12878
echo -e "********NA12878********"
echo "Generating annotation files"
../generateannotations.pl
../svscore.pl -o max,sum,top5,top10,mean -e refGene.exons.bed -f refGene.introns.bed -h .. -dvc ../$1 -i NA12878.sv.vcf | grep -v "^#" | perl -ne 's/(left_|right_)?Intron[^;]*;//g; s/(left_|right_)?ExonTranscript=[^;]*;//g; print' > NA12878.sv.svscore.test.vcf
echo -e "\n\n"

diff NA12878.sv.svscore.test.vcf NA12878.sv.svscore.noheader.vcf > NA12878.diff

NA12878OUTPUTLINES=`cat NA12878.sv.svscore.test.vcf | wc -l`
NA12878DIFFLINES=`cat NA12878.diff | wc -l`

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

rm -rf /svscoretmp/
rm -f NA12878.sv.svscore.noheader.vcf NA12878.sv.svscore.test.vcf NA12878.diff stresstest.svscore.noheader.vcf stresstest.svscore.test.vcf stresstest.diff dummyannotations.*.bed* refGene* conf.toml
