#!/usr/bin/perl -w


use strict;
use Getopt::Std;
use List::Util qw(max min);

$Getopt::Std::STANDARD_HELP_VERSION = 1; # Make --help and --version flags halt execution
$main::VERSION = '0.1';

my %options = ();
getopts('dzasc:',\%options);

my $compressed = defined $options{'z'};
my $debug = defined $options{'d'};
my $annotated = defined $options{'a'};
my $support = defined $options{'s'};
my $cadd = defined $options{'c'};

my $caddfile = ($cadd ? $options{'c'} : '/gscmnt/gc2719/halllab/src/gemini/data/whole_genome_SNVs.tsv.compressed.gz');

##TODO PRIORITY 2: Enable piping input through STDIN - use an option to specify input file rather than @ARGV
##TODO PRIORITY 1: Usage message

# Set up all necessary preprocessing to be taken care of before analysis can begin. This includes decompression (if necessary), annotation using vcfanno, and generation of intron/exon/gene files if necessary. May be a little slower than necessary in certain situations because some arguments are supplied by piping cat output rather than supplying filenames directly.
unless (-s 'refGene.exons.b37.bed') { # Generate exon file if necessary
  print STDERR "Generating exon file\n" if $debug;
  system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); n=int(\$9); split(\$10,start,\",\");split(\$11,end,\",\"); for(i=1;i<=n;++i) {print \$3,start[i],end[i],\$2\".\"i,\$13,\$2; } }' OFS=\"\t\" | sort -V -k 1,1 -k 2,2 -k 3,3 | uniq > refGene.exons.b37.bed");
}

unless (-s 'refGene.genes.b37.bed') { # Generate gene file if necessary
  print STDERR "Generating gene file\n" if $debug;
  system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); print \$3,\$5,\$6,\$4,\$13}' OFS=\"\\t\" | sort -V -k 1,1 -k 2,2 -k 3,3 | uniq > refGene.genes.b37.bed");
}

unless (-s 'introns.bed') { # Generate intron file if necessary - add column with unique intron ID equal to line number (assuming introns.bed has no header line) and sort
  print STDERR "Generating intron file\n" if $debug;
  system("bedtools subtract -a refGene.genes.b37.bed -b refGene.exons.b37.bed | sort -Vu -k 1,1 -k 2,2 -k 3,3 | awk '{print \$0 \"\t\" NR}' > introns.bed");
}

my $wroteconfig = 0;
# Write conf.toml file
unless ($annotated || -s "conf.toml") {
  print STDERR "Writing config file\n" if $debug;
  $wroteconfig = 1;
  open(CONFIG, "> conf.toml") || die "Could not open conf.toml: $!";
  print CONFIG "[[annotation]]\nfile=\"refGene.genes.b37.bed\"\nnames=[\"Gene\"]\ncolumns=[5]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"refGene.exons.b37.bed\"\nnames=[\"ExonGeneNames\"]\ncolumns=[5]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"introns.bed\"\nnames=[\"Intron\"]\ncolumns=[6]\nops=[\"uniq\"]";
  close CONFIG;
}

print STDERR "Preparing preprocessing command\n" if $debug;
die unless defined $ARGV[0];
my ($prefix) = ($ARGV[0] =~ /^(.*)\.vcf$/);
my $preprocessedfile;
if (!$compressed && $annotated) { # No need to copy the file if it's already preprocessed
  $preprocessedfile = $ARGV[0];
} else {
  $preprocessedfile = "$prefix.ann.vcf";
  # Create preprocessing command - annotation is done without normalization because REF and ALT nucleotides are not included in VCFs describing SVs
  my $preprocess = ($compressed ? "z": "") . "cat $ARGV[0]" . ($annotated ? "" : " | vcfanno -ends -natural-sort conf.toml -") . " > $preprocessedfile";

  # Execute preprocessing command
  print STDERR "Preprocessing command: $preprocess\n" if $debug;
  die "vcfanno failed: $!" if system("$preprocess");

  unlink "conf.toml" if $wroteconfig && !$debug;
}

die if $support;

print STDERR "Reading gene list\n" if $debug;
my %genes = (); # Symbol => (chrom, start, stop, strand)
open(GENES, "< refGene.genes.b37.bed") || die "Could not open refGene.genes.b37.bed: $!";
foreach my $geneline (<GENES>) { # Parse gene file, recording the chromosome, strand, 5'-most start coordinate, and 3'-most stop coordinate found 
  my ($genechrom, $genestart, $genestop, $genestrand, $genesymbol) = split(/\s+/,$geneline);
  if (defined $genes{$genesymbol}) { ## Assume chromosome and strand stay constant
    $genes{$genesymbol}->[1] = min($genes{$genesymbol}->[1], $genestart);
    $genes{$genesymbol}->[2] = max($genes{$genesymbol}->[2], $genestop);
  } else {
    $genes{$genesymbol} = [$genechrom, $genestart, $genestop, $genestrand];
  }
}

open(IN, "< $preprocessedfile") || die "Could not open $preprocessedfile: $!";

my %processedids = ();
my $outputfile = "$prefix.scored.txt" ;
open(OUT, "> $outputfile") || die "Could not open output file: $!";
print OUT "Span_score\tSpan_genenames\tSpan_exon\tLeft_score\tLeft_genenames\tLeft_exon\tRight_score\tRight_genenames\tRight_exon\n";

print STDERR "Entering loop\n" if $debug;
my $linenum = 1;
foreach my $vcfline (<IN>) {
  $linenum++ && next if ($vcfline =~ /^#/);
  print STDERR "$linenum," if $debug;

  # Parse line
  my ($leftchrom, $leftpos, $id, $alt, $info) = (split(/\s+/,$vcfline))[0..2, 4, 7];
  my ($mateid,$rightchrom,$rightpos,$rightstart,$rightstop,$leftstart,$leftstop);

  $info =~ /SVTYPE=(\w+);/;
  my $svtype = $1;
  my ($spanexongenenames,$spangenenames) = getfields($info,"ExonGeneNames","Gene");
  my ($leftexongenenames,$leftgenenames) = getfields($info,"left_ExonGeneNames","left_Gene");
  my @leftgenenames = split(/,/,$leftgenenames);
  my ($rightexongenenames,$rightgenenames) = getfields($info,"right_ExonGeneNames","right_Gene");
  my @rightgenenames = split(/,/,$rightgenenames);
  my $leftintrons = getfields($info, "left_Intron") if $svtype eq 'BND' || $svtype eq 'INV';
  my $rightintrons = getfields($info, "right_Intron") if $svtype eq 'BND' || $svtype eq 'INV';

  if ($svtype eq 'BND') {
    my ($mateid) = ($id =~ /(\d+)_(?:1|2)/);
    if (exists $processedids{$mateid}) { # Is this the second mate of this BND seen? If so, get annotations for first mate. Otherwise, store in %processedids for later when other mate is found
      ($rightexongenenames, $rightgenenames, $rightintrons) = @{$processedids{$mateid}}; # "Right side" comes from primary mate - assumes left and right treated equally later (currently true)
      delete $processedids{$mateid};
      ($rightchrom,$rightpos) = ($alt =~ /N?[\[\]](.*):(\d+)/); # Grab right (primary) breakend coordinates
    } else {
      $processedids{$mateid} = [$leftexongenenames,$leftgenenames, $leftintrons];
      $linenum++;
      next;
    }
  } else { ## DEL, DUP, INV, or INS
    $rightchrom = $leftchrom;
    $rightpos = getfields($info,"END");
  }

  # Calculate start/stop coordinates from POS and CI
  my ($cipos, $ciend) = getfields($info,"CIPOS","CIEND");
  my @cipos = split(/,/,$cipos);
  my @ciend = split(/,/,$ciend);
  $leftstart = $leftpos + $cipos[0] - 1;
  $leftstop = $leftpos + $cipos[1];
  $rightstart = $rightpos + $ciend[0] - 1;
  $rightstop = $rightpos + $ciend[1];

  my ($spanscore, $leftscore, $rightscore);

  # Calculate maximum C score depending on SV type
  if ($svtype eq "DEL" || $svtype eq "DUP") {
    if ($rightstop - $leftstart > 1000000) { ## Hack to avoid extracting huge regions from CADD file
      print OUT "$svtype too big at line $linenum: $leftstart-$rightstop";
      $linenum++;
      next;
    }
    $spanscore = maxcscore($caddfile, $leftchrom, $leftstart, $rightstop);
    $leftscore = maxcscore($caddfile, $leftchrom, $leftstart, $leftstop);
    $rightscore = maxcscore($caddfile, $rightchrom, $rightstart, $rightstop);
  } elsif ($svtype eq "INV" || $svtype eq "BND") {
    $leftscore = maxcscore($caddfile, $leftchrom, $leftstart, $leftstop);
    $rightscore = maxcscore($caddfile, $rightchrom, $rightstart, $rightstop);

    my %leftintrons = map {$_ => 1} (split(/,/,$leftintrons));
    my @rightintrons = split(/,/,$rightintrons);
    ## At worst, $leftintrons and $rightintrons are lists of introns. The only case in which the gene is not disrupted is if both lists are equal and nonempty, meaning that in every gene hit by this variant, both ends of the variant are confined to the same intron
    my $sameintrons = scalar (grep {$leftintrons{$_}} @rightintrons) == scalar @rightintrons && scalar @rightintrons > 0;
    if ((!$leftgenenames && !$rightgenenames) || $sameintrons) { # Either breakends don't hit genes or some genes are involved, but ends of variant are within same intron in each gene hit
      $spanscore = "${svtype}SameIntrons" if $sameintrons;
      $spanscore = "${svtype}NoGenes" unless $sameintrons;
    } else { # Consider variant to be truncating the gene(s)
      my @lefttruncationscores = ();
      foreach my $gene (split(/,/,$leftgenenames)) {
	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}}[1..3];	
	if ($genestrand eq '+') {
	  push @lefttruncationscores, maxcscore($caddfile, $leftchrom, max($genestart,$leftstart), $genestop); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
	} else { ## Minus strand
	  push @lefttruncationscores,maxcscore($caddfile, $leftchrom, $genestart, min($genestop,$leftstop)); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
	}
      }
      my @righttruncationscores = ();
      foreach my $gene (split(/,/,$rightgenenames)) {
	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}}[1..3];	
	if ($genestrand eq '+') {
	  push @righttruncationscores, maxcscore($caddfile, $rightchrom, max($genestart,$rightstart), $genestop); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
	} else { ## Minus strand
	  push @righttruncationscores,maxcscore($caddfile, $rightchrom, $genestart, min($genestop,$rightstop)); # Start from beginning of gene, stop at end of gene or breakend, whichever is further right (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
	}
      }
      my $lefttruncationscore = (@lefttruncationscores ? max(@lefttruncationscores) : "NoLeftGenes");

      my $righttruncationscore = (@righttruncationscores ? max(@righttruncationscores) : "NoRightGenes");
      $spanscore = "$lefttruncationscore,$righttruncationscore";
    }
    $spangenenames = "$svtype:IgnoredSpan";
    $spanexongenenames = "$svtype:IgnoredSpan";
  } elsif ($svtype eq "INS") { # leftscore is base before insertion, rightscore is base after insertion
    $leftscore = maxcscore($caddfile, $leftchrom, $leftstart-1, $leftstart-1);
    $rightscore = maxcscore($caddfile, $rightchrom, $rightstart+1, $rightstart+1);
    $spanscore = "INS";
    $spangenenames = "INS";
    $spanexongenenames = "INS";
  } else {
    die "Unrecognized SVTYPE $svtype at line $linenum of annotated VCF file\n";
  }

  # Multiplier for deletions and duplications which hit an exon, lower multiplier if one of these hits a gene but not an exon. Purposely not done for BND and INV
  if ($spanexongenenames && ($svtype eq "DEL" || $svtype eq "DUP")) { 
    $spanscore *= 1.5;
  } elsif($spangenenames && ($svtype eq "DEL" || $svtype eq "DUP")) {
    $spanscore *= 1.2;
  }

  # For all types except INS, multiply left and right scores respectively if exon/gene is hit
  if ($leftexongenenames && $svtype ne "INS") {
    $leftscore *= 1.5;
  } elsif ($leftgenenames && $svtype ne "INS") {
    $leftscore *= 1.2;
  }
  if ($rightexongenenames && $svtype ne "INS") {
    $rightscore *= 1.5;
  } elsif ($rightgenenames && $svtype ne "INS") {
    $rightscore *= 1.2;
  }

  # We don't like empty strings
  unless ($spangenenames) {
#    if ($svtype eq "BND" || $svtype eq "INV") {
#      $spangenenames = "$svtype:IgnoredSpan";
#    } else {
      $spangenenames = "NoGenesInSpan";
#    }
  }
  unless ($spanexongenenames) {
#    if ($svtype eq "BND" || $svtype eq "INV") {
#      $spangenenames = "$svtype:IgnoredSpan";
#    } else {
      $spanexongenenames = "NoExonsInSpan";
#    }
  }
  unless ($leftgenenames) {
    $leftgenenames = "NoGenesInLeftBreakend";
  }
  unless ($leftexongenenames) {
    $leftexongenenames = "NoExonsInLeftBreakend";
  }
  unless ($rightgenenames) {
    $rightgenenames = "NoGenesInRightBreakend";
  }
  unless ($rightexongenenames) {
    $rightexongenenames = "NoExonsInRightBreakend";
  } 
  print OUT "$spanscore\t$spangenenames\t$spanexongenenames\t$leftscore\t$leftgenenames\t$leftexongenenames\t$rightscore\t$rightgenenames\t$rightexongenenames\n";

  $linenum++;
}
unlink "$preprocessedfile" unless (!$compressed && $annotated) || $debug;

sub maxcscore { # Calculate maximum C score within a given region using CADD data
  my ($filename, $chrom, $start, $stop) = @_;
  my @scores = ();
  my $tabixoutput = `tabix $filename $chrom:$start-$stop`;
  my @tabixoutputlines = split(/\n/,$tabixoutput);
  foreach my $taboutline (@tabixoutputlines) {
    push @scores, split(/,/,(split(/\s+/,$taboutline))[4]);
  }
  
  return max(@scores);
}

sub getfields { # Parse info field of VCF line, getting fields specified in @_. $_[0] must be the info field itself. Returns list of field values
  my $info = shift @_;
  my @ans;
  foreach my $i (0..$#_) {
    my ($ann) = ($info =~ /(?:;|^)$_[$i]=(.*?)(?:;|$)/);
    push @ans, ($ann ? $ann : "");
  }
  if (@ans > 1) {
    return @ans;
  } else {
    return $ans[0];
  }
}

sub main::HELP_MESSAGE() {
  print "usage: ./svscore.pl [-dzas] [-c file] vcf
    -d	      Debug (verbose) mode, keeps intermediate and supporting files
    -z	      Indicates that vcf is gzipped
    -a	      Indicates that vcf has already been annotated using vcfanno
    -s	      Create/download supporting files and quit
    -c	      Used to point to whole_genome_SNVs.tsv.gz
    --help    Display this message
    --version Display version

    vcf must be sorted naturally (e.g. chromosome 1, chromosome 2,...,chromosome 9, chromosome 10,...)\n"
}

sub main::VERSION_MESSAGE() {
  print "SVScore version $main::VERSION\n";
}
