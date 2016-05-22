if [ -z "$1" ]; then
  echo "Please provide the path to whole_genome_SNVs.tsv.gz as an argument"
  exit 1
fi
chmod +x generateannotations.pl svscore.pl tests/test.sh
cd tests
sh test.sh $1
cd ..
