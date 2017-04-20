if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1
fi

testdir=$(dirname "$(readlink -f "$0")")
svscoredir=$(readlink -f $testdir"/..")

echo "********Testing SVScore. Some warnings beginning with \"Warning: missing secondary multiline variant at ID:\" are normal********"
echo

# Perform stress test
echo "**************"
echo "*            *"
echo "* STRESSTEST *"
echo "*            *"
echo "**************"

echo "Generating annotation files"
$svscoredir/generateannotations.pl -c 1 -a 2 -b 3 -t 4 -s 5 -e 6 -f 7  -o $testdir $testdir/dummyannotations.bed
$svscoredir/svscore.pl -f $testdir/dummyannotations.introns.bed -e $testdir/dummyannotations.exons.bed -o max,sum,top2,top2weighted,top3weighted,top4weighted,mean,meanweighted -dvc $testdir/dummyCADD.tsv.gz -i $testdir/stresstest.vcf > $testdir/stresstest.svscore.test.vcf
perl $testdir/scorecomparison.pl STRESSTEST $testdir/stresstest.svscore.vcf $testdir/stresstest.svscore.test.vcf
echo
echo

# NA12878
echo "***********"
echo "*         *"
echo "* NA12878 *"
echo "*         *"
echo "***********"
echo "Generating annotation files"
$svscoredir/generateannotations.pl -o $testdir
$svscoredir/svscore.pl -o max,sum,top5,top10,mean -e $testdir/refGene.exons.bed -f $testdir/refGene.introns.bed -dvc $1 -i $testdir/NA12878.sv.vcf > $testdir/NA12878.sv.svscore.test.vcf
perl $testdir/scorecomparison.pl NA12878 $testdir/NA12878.sv.svscore.vcf $testdir/NA12878.sv.svscore.test.vcf

rm -rf svscoretmp/
rm -f $testdir/NA12878.sv.svscore.test.vcf $testdir/stresstest.svscore.test.vcf $testdir/dummyannotations.*.bed* $testdir/refGene* conf.toml
