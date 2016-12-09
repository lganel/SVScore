if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1
fi

echo "********Testing SVScore. Some warnings beginning with \"Warning: missing secondary multiline variant at ID:\" are normal********"
echo

# Perform stress test
echo "**************"
echo "*            *"
echo "* STRESSTEST *"
echo "*            *"
echo "**************"

echo "Generating annotation files"
$SVSCOREDIR/generateannotations.pl -c 1 -a 2 -b 3 -t 4 -s 5 -e 6 -f 7  -o $SVSCOREDIR/tests $SVSCOREDIR/tests/dummyannotations.bed
$SVSCOREDIR/svscore.pl -f $SVSCOREDIR/tests/dummyannotations.introns.bed -e $SVSCOREDIR/tests/dummyannotations.exons.bed -o max,sum,top2,top2weighted,top3weighted,top4weighted,mean,meanweighted -dvc $SVSCOREDIR/tests/dummyCADD.tsv.gz -i $SVSCOREDIR/tests/stresstest.vcf > $SVSCOREDIR/tests/stresstest.svscore.test.vcf
perl $SVSCOREDIR/tests/scorecomparison.pl STRESSTEST $SVSCOREDIR/tests/stresstest.svscore.vcf $SVSCOREDIR/tests/stresstest.svscore.test.vcf
echo
echo

# NA12878
echo "***********"
echo "*         *"
echo "* NA12878 *"
echo "*         *"
echo "***********"
echo "Generating annotation files"
$SVSCOREDIR/generateannotations.pl -o $SVSCOREDIR/tests
$SVSCOREDIR/svscore.pl -o max,sum,top5,top10,mean -e $SVSCOREDIR/tests/refGene.exons.bed -f $SVSCOREDIR/tests/refGene.introns.bed -dvc $1 -i $SVSCOREDIR/tests/NA12878.sv.vcf > $SVSCOREDIR/tests/NA12878.sv.svscore.test.vcf
perl $SVSCOREDIR/tests/scorecomparison.pl NA12878 $SVSCOREDIR/tests/NA12878.sv.svscore.vcf $SVSCOREDIR/tests/NA12878.sv.svscore.test.vcf

rm -rf svscoretmp/
rm -f $SVSCOREDIR/tests/NA12878.sv.svscore.test.vcf $SVSCOREDIR/tests/stresstest.svscore.test.vcf $SVSCOREDIR/tests/dummyannotations.*.bed* $SVSCOREDIR/tests/refGene* conf.toml
