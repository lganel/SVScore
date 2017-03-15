#!/usr/bin/perl -w

use strict;
use Getopt::Std;

sub main::HELP_MESSAGE();
$Getopt::Std::STANDARD_HELP_VERSION = 1;
$main::VERSION = '0.1';

my %options = ();
getopts('a:b:c:t:g:e:f:s:o:',\%options);

my $chromcolumn = (defined $options{'c'} ? $options{'c'} : 3);
my $startcolumn = (defined $options{'a'} ? $options{'a'} : 5);
my $stopcolumn = (defined $options{'b'} ? $options{'b'} : 6);
my $transcriptcolumn = (defined $options{'t'} ? $options{'t'} : 2);
my $exonstartcolumn = (defined $options{'e'} ? $options{'e'} : 10);
my $exonstopcolumn = (defined $options{'f'} ? $options{'f'} : 11);
my $strandcolumn = (defined $options{'s'} ? $options{'s'} : 4);
my $outputdir = (defined $options{'o'} ? $options{'o'} : "");
$outputdir .= "/" if $options{'o'} && $outputdir !~ /\/$/;
my ($inputfile,$prefix);
if (defined $ARGV[0]) {
  $inputfile = $ARGV[0];
  ($prefix) = ($ARGV[0] =~ /(?:.*\/)?(.*)\.[^.]*(?:\.gz)?/);
} else {
  $inputfile = 'curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz" | gzip -cdfq |';
  $prefix = "refGene";
}

my $intronfile = "$outputdir$prefix.introns.bed";
my $exonfile = "$outputdir$prefix.exons.bed";
open(IN, $inputfile) || die "Download/opening of refGene file failed";
open(INTRON,"> $intronfile") || die "Could not open $intronfile for writing";
open(EXON,"> $exonfile") || die "Could not open $exonfile for writing";

my %names = ();
my $intronnumber = 1;
while(my $line = <IN>) {
  next if $line =~ /^#/;
  my ($chrom,$start,$stop,$transcript,$exonstart,$exonstop,$strand) = (split(/\s+/,$line))[$chromcolumn-1,$startcolumn-1,$stopcolumn-1,$transcriptcolumn-1,$exonstartcolumn-1,$exonstopcolumn-1,$strandcolumn-1];
  if (exists $names{$transcript}) {
    my $i = 0;
    my $original = $transcript;
    do {
      $i++;
      $transcript = "$original.$i";
    } while (exists $names{$transcript});
  }
  $names{$transcript} = undef;
  $chrom =~ s/^chr//;
  my @exonstarts = split(/,/,$exonstart);
  my @exonstops = split(/,/,$exonstop);

  unless (@exonstarts==@exonstops) {
    warn "Transcript $transcript has an unequal number of exon starts and stops. Skipping $transcript";
    next;
  }

  if ($start < $exonstarts[0]) {
    print INTRON "$chrom\t$start\t$exonstarts[0]\t$transcript\t0\n" if $start < $exonstarts[0];
    $intronnumber++;
  }
  foreach my $i (0..$#exonstarts) {
    print EXON "$chrom\t$exonstarts[$i]\t$exonstops[$i]\t$transcript\t$start\t$stop\t$strand\n";
    unless ($i==$#exonstarts) {
      print INTRON "$chrom\t$exonstops[$i]\t$exonstarts[$i+1]\t$transcript\t$intronnumber\n";
      $intronnumber++;
    }
  }
  if ($stop > $exonstops[$#exonstops]) {
    print INTRON "$chrom\t$exonstops[$#exonstops]\t$stop\t$transcript\t" . $intronnumber . "\n";
    $intronnumber++;
  }
}

close EXON;
close INTRON;
close IN;

my $sortedintronfile = "$outputdir$prefix.introns.sort.bed";
my $sortedexonfile = "$outputdir$prefix.exons.sort.bed";

!system("sort -k1,1V -k2,2n -k3,3n $intronfile | uniq > $sortedintronfile;  mv -f $sortedintronfile $intronfile") || die "Intron file sorting failed";
!system("sort -k1,1V -k2,2n -k3,3n $exonfile | uniq > $sortedexonfile;  mv -f $sortedexonfile $exonfile") || die "Exon file sorting failed";

sub main::HELP_MESSAGE(){
  print STDERR "usage: ./generateannotations.pl [-o outputdir] [-c chromcolumn] [-a startcolumn] [-b stopcolumn] [-t transcriptcolumn] [-e exonstartcolumn] [-f exonstopcolumn] [-s strandcolumn] [file]
    -o	      Path to directory in which to save annotation files (current directory)
    -c	      Chromosome column (3)
    -a	      Transcript start column (5)
    -b	      Transcript stop column (6)
    -t	      Transcript name column (2)
    -e	      Exon start positions column (must be comma delimited) (10)
    -f	      Exon stop positions column (must be comma delimited) (11)
    -s	      Strand column (4)
    
    --help    Display this message
    --version Display version\n"
}

sub main::VERSION_MESSAGE() {
  print "generateannotations.pl version $main::VERSION\n";
}
