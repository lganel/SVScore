if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1
fi

echo -e "********Testing SVScore. Some warnings beginning with \"Warning: missing secondary multiline variant at ID:\" are normal********\n\n"

# Perform stress test
echo "**************"
echo "*            *"
echo "* STRESSTEST *"
echo "*            *"
echo "**************"

echo "Generating annotation files"
$SVSCOREDIR/generateannotations.pl -c 1 -a 2 -b 3 -t 4 -s 5 -e 6 -f 7  -o tests tests/dummyannotations.bed
$SVSCOREDIR/svscore.pl -f tests/dummyannotations.introns.bed -e tests/dummyannotations.exons.bed -o max,sum,top2,top2weighted,top3weighted,top4weighted,mean,meanweighted -dvc tests/dummyCADD.tsv.gz -i tests/stresstest.vcf > tests/stresstest.svscore.test.vcf
perl tests/scorecomparison.pl STRESSTEST tests/stresstest.svscore.vcf tests/stresstest.svscore.test.vcf
echo
echo

# NA12878
echo "***********"
echo "*         *"
echo "* NA12878 *"
echo "*         *"
echo "***********"
echo "Generating annotation files"
$SVSCOREDIR/generateannotations.pl -o tests
$SVSCOREDIR/svscore.pl -o max,sum,top5,top10,mean -e tests/refGene.exons.bed -f tests/refGene.introns.bed -dvc $1 -i tests/NA12878.sv.vcf > tests/NA12878.sv.svscore.test.vcf
perl tests/scorecomparison.pl NA12878 tests/NA12878.sv.svscore.vcf tests/NA12878.sv.svscore.test.vcf

rm -rf svscoretmp/
rm -f tests/NA12878.sv.svscore.test.vcf tests/stresstest.svscore.test.vcf tests/dummyannotations.*.bed* tests/refGene* conf.toml
