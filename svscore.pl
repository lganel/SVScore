#!/usr/bin/perl -w 
## Author: Liron Ganel
## Laboratory of Ira Hall, McDonnell Genome Institute
## Washington University in St. Louis
## Version 0.5

use strict;
use Getopt::Std;
use List::Util qw(max min sum);
use List::MoreUtils qw(pairwise uniq);
use Time::HiRes qw(gettimeofday);
use Math::Round qw(nearest);

sub main::HELP_MESSAGE(); # Declaration to make perl happy
$Getopt::Std::STANDARD_HELP_VERSION = 1; # Make --help and --version flags halt execution
$main::VERSION = '0.5';

my %possibleoperations = ("MAX", 1, "SUM", 1, "TOP", 1, "MEAN", 1, "MEANWEIGHTED", 1, "TOPWEIGHTED", 1); # Hash of supported operations
my %types = map {$_ => 1} ("DEL", "DUP", "INV", "BND", "TRX", "CNV", "MEI", "INS"); # Hash of supported svtypes
my %truncationtypes = ("INV","TRX","INS","DEL"); # All svtypes for which truncation scores are calculated
my %intervals = ("LEFT", 0, "RIGHT", 1, "SPAN", 2, "LTRUNC", 3, "RTRUNC", 4); # Hash of supported intervals

my %options = ();
getopts('dvc:g:e:m:o:p:i:',\%options);

# Parse command line options, set variables, check input parameters
my $debug = defined $options{'d'};
#my $support = defined $options{'s'};
#my $weight = defined $options{'w'};
my $verbose = defined $options{'v'};

&main::HELP_MESSAGE() && die unless defined $options{'i'};

my $caddfile = (defined $options{'c'} ? $options{'c'} : 'whole_genome_SNVs.tsv.gz');
die "Could not find $caddfile" unless -s $caddfile;
my $genefile = (defined $options{'g'} ? $options{'g'} : 'refGene.genes.b37.bed.gz');
my $geneanncolumn = (defined $options{'m'} && defined $options{'g'} ? $options{'m'} : 4);
$geneanncolumn = 4 && warn "Gene annotation column provided without nonstandard gene annotation file - defaulting to standard gene annotation column (4)" if defined $options{'m'} && !defined $options{'g'};
$geneanncolumn = 4 && warn "Gene annotation column must be greater than 3 - defaulting to standard gene annotation column (4)" if $geneanncolumn <= 3;
my $genestrandcolumn = (defined $options{'p'} ? $options{'p'} : 5);
$genestrandcolumn = 5 && warn "Gene strand column provided without nonstandard gene annotation file - defaulting to standard gene strand column (5)" if defined $options{'m'} && !defined $options{'g'};
$genestrandcolumn = 5 && warn "Gene annotation column must be greater than 3 - defaulting to standard gene strand column (5)" if $genestrandcolumn <= 3;
die "Gene annotation column cannot equal gene strand column" if $geneanncolumn==$genestrandcolumn;
my $exonfile = (defined $options{'e'} ? $options{'e'} : 'refGene.exons.b37.bed');
$options{'o'} =~ tr/[a-z]/[A-Z]/ if defined $options{'o'};
my $ops = (defined $options{'o'} ? $options{'o'} : 'TOP10WEIGHTED');
#my $topn = (defined $options{'t'} ? $options{'t'} : 100);
my @ops = uniq(split(/,/,$ops));
my %operations = (); # Hash of chosen operations
foreach my $i (0..$#ops) { # Populate %operations with chosen operations given by -o
  my $op = $ops[$i];
  $op =~ s/\d+//g; # Get rid of all numbers in $op for lookup in %possibleoperations
  unless (defined $possibleoperations{$op}) {
    warn "Unrecognized operation specified: $ops[$i]\n";
    &main::HELP_MESSAGE() && die;
  }
  $operations{$ops[$i]} = $i;
}
#my %operations = ($ops eq 'ALL' ? %possibleoperations : ($ops => 0)); # Hash of indices for each operation within lists in the values of %scores
#if (defined $operations{"TOP"}) { # Replace "TOP" with "TOP$topn" in %operations keys 
#  $operations{"TOP$topn"} = $operations{"TOP"};
#  delete $operations{"TOP"};
#}
my $inputfile = (defined $options{'i'} ? $options{'i'} : "");
my $compressed = ($inputfile =~ /\.gz$/);
my ($uncompressedgenefile, $compressedgenefile, $uncompressedexonfile, $alteredgenefile, $alteredexonfile, $headerfile, $sortedfile, $preprocessedfile, $bedpeout, $vcfout);

# Set up all necessary preprocessing to be taken care of before analysis can begin. This includes decompression, annotation using vcfanno, and generation of intron/exon/gene files, whichever are necessary
if ($exonfile eq 'refGene.exons.b37.bed' && !-s $exonfile) { # Generate exon file if necessary
  print STDERR "Generating exon file: $exonfile\n" if $debug;
  system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); n=int(\$9); split(\$10,start,\",\");split(\$11,end,\",\"); for(i=1;i<=n;++i) {print \$3,start[i],end[i],\$2\".\"i,\$13,\$2; } }' OFS=\"\t\" | sort -k 1,1V -k 2,2n | uniq > refGene.exons.b37.bed");
} elsif ($exonfile ne 'refGene.exons.b37.bed' && !-s $exonfile) {
  die "$exonfile not found or empty!";
}

my ($geneprefix) = ($genefile =~ /^(?:.*\/)?(.*)\.bed(?:\.gz)?$/);
my ($exonprefix) = ($exonfile =~ /^(?:.*\/)?(.*)\.bed(?:\.gz)?$/);
my $intronfile = "introns.$geneprefix.$exonprefix.bed.gz";
if ($genefile eq 'refGene.genes.b37.bed.gz' && !-s $genefile && !-s $intronfile) { # Generate gene file if necessary
  print STDERR "Generating gene file: $genefile\n" if $debug;
  system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); print \$3,\$5,\$6,\$13,\$4}' OFS=\"\\t\" | sort -k 1,1V -k 2,2n | uniq | bgzip -c > refGene.genes.b37.bed.gz");
} elsif ($genefile ne 'refGene.genes.b37.bed.gz' && !-s $genefile) {
  die "$genefile not found or empty!";
}

# Zip/unzip annotation files as necessary
if ($genefile =~ /\.gz/) {
  $compressedgenefile = $genefile;
  ($uncompressedgenefile) = ($genefile =~ /(.*)\.gz$/);
  unless (-s $uncompressedgenefile) {
    $alteredgenefile = 1;
    print STDERR "Unzipping $genefile\n" if $debug;
    if (system("zcat $genefile > $uncompressedgenefile")) {
      die "Unzipping $genefile failed: $!";
    }
  }
} else {
  $uncompressedgenefile = $genefile;
  $compressedgenefile = "$genefile.gz";
  unless(-s $compressedgenefile) {
    $alteredgenefile = 1;
    print STDERR "Zipping $genefile\n" if $debug;
    if (system("bgzip -c $genefile > $compressedgenefile")) {
      die "Compressing $genefile failed: $!";
    }
  }
}

if ($exonfile =~ /\.gz/) {
  ($uncompressedexonfile) = ($exonfile =~ /(.*)\.gz$/);
  unless (-s $uncompressedexonfile) {
    $alteredexonfile = 1;
    print STDERR "Unzipping $exonfile\n" if $debug;
    if (system("zcat $exonfile > $uncompressedexonfile")) {
      die "Unzipping $exonfile failed: $!";
    }
  }
} else {
  $uncompressedexonfile = $exonfile;
}

mkdir "svscoretmp" unless -d "svscoretmp";

unless (-s "$intronfile") { # Generate intron file if necessary - add column with unique intron ID equal to line number (assuming intron file has no header line) and sort
  print STDERR "Generating intron file\n" if $debug;
  system("bedtools subtract -a $uncompressedgenefile -b $uncompressedexonfile | sort -u -k1,1V -k2,2n | awk '{print \$0 \"\t\" NR}' | bgzip -c > $intronfile");
}

# Use tabix to index the annotation files
unless (-s "$intronfile.tbi") {
  print STDERR "Tabix indexing $intronfile\n" if $debug;
  if(system("tabix -p bed $intronfile")) {
    die "Tabix failed on $intronfile";
  }
}
unless (-s "$compressedgenefile.tbi") {
  print STDERR "Tabix indexing $compressedgenefile\n" if $debug;
  if(system("tabix -p bed $compressedgenefile")) {
    die "Tabix failed on $compressedgenefile";
  }
}

my $intronnumcolumn = `zcat $intronfile | head -n 1 | awk '{print NF}'`; # Figure out which column has the intron number
chomp $intronnumcolumn;

# Write conf.toml file
print STDERR "Writing config file\n" if $debug;
open(CONFIG, "> conf.toml") || die "Could not open conf.toml: $!";
print CONFIG "[[annotation]]\nfile=\"$compressedgenefile\"\nnames=[\"Gene\"]\ncolumns=[$geneanncolumn]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"$intronfile\"\nnames=[\"Intron\",\"IntronGene\"]\ncolumns=[$intronnumcolumn,4]\nops=[\"uniq\",\"uniq\"]\n";
close CONFIG;

# Create first preprocessing command - annotation is done without normalization because REF and ALT nucleotides are not included in VCFs describing SVs
my ($prefix,$tempfile);
my $time = gettimeofday();
if ($inputfile eq "stdin") {
  $tempfile = "svscoretmp/stdin$time.vcf";
  open(TEMP,">$tempfile") || die "Could not open $tempfile for writing; $!";
  print TEMP <STDIN>;
  close TEMP;
}

if ($inputfile) {
  print STDERR "Preparing preprocessing command\n" if $debug;
  if ($inputfile eq "stdin") {
    $prefix = "stdin";
    $inputfile = $tempfile;
  } else {
    ($prefix) = ($inputfile =~ /^(?:.*\/)?(.*)\.vcf(?:\.gz)?$/);
  }
  # Tag intermediate files with timestamp to avoid collisions
  $preprocessedfile = "svscoretmp/$prefix$time.preprocess.bedpe";
  my $sortedfile = "svscoretmp/$prefix$time.sort.vcf.gz";
  my $reorderout = "svscoretmp/$prefix$time.reorderheaderout.vcf";
  if ($compressed) {
    system("gunzip -c $inputfile > svscoretmp/$prefix$time.vcf") && die "Could not unzip $inputfile: $!";
    $inputfile = "svscoretmp/$prefix$time.vcf";
  }
  my $preprocess = "awk '\$0~\"^#\" {print \$0; next } { print \$0 | \"sort -k1,1V -k2,2n\" }' $inputfile | bgzip -c > $sortedfile; tabix -p vcf $sortedfile; vcfanno -ends conf.toml $sortedfile | perl reorderheader.pl stdin $inputfile > $reorderout"; # Sort, annotate, reorder header
  my $preprocess2 = "svtools vcftobedpe -i $reorderout > $preprocessedfile; rm -f $sortedfile $sortedfile.tbi $reorderout";
  print STDERR "Preprocessing commands:\n$preprocess\n$preprocess2\n" if $debug;
  if (system($preprocess) || system($preprocess2) || -z $preprocessedfile) {
    unless ($debug) {
      unlink $preprocessedfile;
      rmdir "svscoretmp";
    }
    die "Preprocessing failed: $!"
  }

  $bedpeout = "svscoretmp/$prefix$time.out.bedpe";
  $vcfout = "svscoretmp/$prefix$time.out.vcf";
  open(IN, "< $preprocessedfile") || die "Could not open $preprocessedfile: $!";
  open(OUT, "> $bedpeout") || die "Could not open output file: $!";

  # Update header
  open(HEADER, "grep \"^#\" $preprocessedfile |") || die "Error grabbing header: $!";
  my @newheader = ();
  my @oldheader = <HEADER>;
  close HEADER;

  if ($ops =~ /WEIGHTED/ && !(grep {/^##INFO=<ID=PRPOS,/} @oldheader)) { # If an op is weighted, but PRPOS is absent from the header, switch to unweighted operations with a warning
    warn "*****PRPOS not found in header - switching to unweighted operations*****\n";
    $ops =~ s/WEIGHTED//g;
    foreach my $op (keys %operations) {
      if ($op =~ /WEIGHTED/) {
	my $value = $operations{$op};
	delete $operations{$op};
	$op =~ s/WEIGHTED//;
	$operations{$op} = $value;
      }
    }
  }

  my $headerline;
  while(($headerline = (shift @oldheader)) !~ /^##INFO/) {
    push @newheader, $headerline;
  }
  unshift @oldheader, $headerline; # Return first info line to top of stack
  while(($headerline = (shift @oldheader)) =~ /^##INFO/) {
    push @newheader, $headerline;
  }
  unshift @oldheader, $headerline; # Return first format line to top of stack

  if (defined $operations{"MAX"}) {
    push @newheader, "##INFO=<ID=SVSCOREMAX,Number=1,Type=Float,Description=\"Maximum of SVSCORE_MAX fields of structural variant\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_LEFT,Number=1,Type=Float,Description=\"Maximum C score in left breakend of structural variant\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_RIGHT,Number=1,Type=Float,Description=\"Maximum C score in right breakend of structural variant\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_SPAN,Number=1,Type=Float,Description=\"Maximum C score in span of structural variant\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_LTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of left breakend to end of truncated gene\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_RTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of right breakend to end of truncated gene\">\n";
  }
  if (defined $operations{"SUM"}) {
    push @newheader, "##INFO=<ID=SVSCORESUM,Number=1,Type=Float,Description=\"Maximum of SVSCORE_SUM fields of structural variant\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_LEFT,Number=1,Type=Float,Description=\"Sum of C scores in left breakend of structural variant\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_RIGHT,Number=1,Type=Float,Description=\"Sum of C scores in right breakend of structural variant\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_SPAN,Number=1,Type=Float,Description=\"Sum of C scores in outer span of structural variant\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_LTRUNC,Number=1,Type=Float,Description=\"Sum of C scores from beginning of left breakend to end of truncated gene\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_RTRUNC,Number=1,Type=Float,Description=\"Sum of C scores from beginning of right breakend to end of truncated gene\">\n";
  }
  foreach my $op (keys %operations) {
    my $weighted = ($op =~ /WEIGHTED$/);
    if ($op =~ /^TOP\d+/) {
      my ($n) = ($op =~ /^TOP(\d+)/);
      push @newheader, "##INFO=<ID=SVSCORETOP$n,Number=1,Type=Float,Description=\"Maximum of SVSCORE_TOP$n fields of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${n}_LEFT,Number=1,Type=Float,Description=\"Mean of top $n C scores in left breakend of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${n}_RIGHT,Number=1,Type=Float,Description=\"Mean of top $n C scores in right breakend of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${n}_SPAN,Number=1,Type=Float,Description=\"Mean of top $n C scores in outer span of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${n}_LTRUNC,Number=1,Type=Float,Description=\"Mean of top $n C scores from beginning of left breakend to end of truncated gene" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${n}_RTRUNC,Number=1,Type=Float,Description=\"Mean of top $n C scores from beginning of right breakend to end of truncated gene" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
    }
   
    if ($op =~ /^MEAN/) {
      push @newheader, "##INFO=<ID=SVSCOREMEAN,Number=1,Type=Float,Description=\"Maximum of SVSCORE_MEAN fields of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_LEFT,Number=1,Type=Float,Description=\"Mean of C scores in left breakend of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_RIGHT,Number=1,Type=Float,Description=\"Mean of C scores in right breakend of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_SPAN,Number=1,Type=Float,Description=\"Mean of C scores in outer span of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_LTRUNC,Number=1,Type=Float,Description=\"Mean of C scores from beginning of left breakend to end of truncated gene" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_RTRUNC,Number=1,Type=Float,Description=\"Mean of C scores from beginning of right breakend to end of truncated gene" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
    }
  }
  push @newheader, @oldheader;
  foreach (uniq(@newheader)) {
    print OUT;
  }

}

print STDERR "Reading gene list\n" if $debug;
my %genes = (); # Symbol => (Chrom => (chrom, start, stop, strand)); Hash of hashes of arrays
open(GENES, "$uncompressedgenefile") || die "Could not open $genefile: $!";
foreach my $geneline (<GENES>) { # Parse gene file, recording the chromosome, strand, 5'-most start coordinate, and 3'-most stop coordinate found  (coordinates are 1-based)
  my ($genechrom, $genestart, $genestop, $genesymbol, $genestrand) = (split(/\s+/,$geneline))[0..2,$geneanncolumn-1,$genestrandcolumn-1];
  $genestart++; ## Move interval to 1-based (VCF) coordinates (don't correct $genestop because BED intervals are non-inclusive of end position, so corrections cancel out)
  if (defined $genes{$genesymbol}->{$genechrom}) { ## Assume strand stays constant
    $genes{$genesymbol}->{$genechrom}->[0] = min($genes{$genesymbol}->{$genechrom}->[0], $genestart);
    $genes{$genesymbol}->{$genechrom}->[1] = max($genes{$genesymbol}->{$genechrom}->[1], $genestop);
  } else {
    $genes{$genesymbol}->{$genechrom} = [$genestart, $genestop, $genestrand];
  }
}
close GENES;

print STDERR "Entering loop\n" if $debug;
while (my $inputline = <IN>) {
  if ($inputline =~ /^#/) {
    next;
  }
  # Parse line
  my @splitline = split(/\s+/,$inputline);
  my ($leftchrom, $leftstart, $leftstop, $rightchrom, $rightstart, $rightstop, $id, $svtype, $info_a, $info_b) = @splitline[0..6, 10, 12, 13];

  $svtype = substr($svtype, 0, 3);
  unless (exists $types{$svtype}) {
    warn "Unrecognized SVTYPE $svtype at line ",$.," of preprocessed VCF file\n";
    print OUT $inputline;
    next;
  }
  my ($probleft,$probright) = getfields($info_a,"PRPOS","PREND") if $ops =~ /WEIGHTED/;

  # Get vcfanno annotations
  my ($leftgenenames, $rightgenenames, $leftintrons, $rightintrons, $leftintrongenenames, $rightintrongenenames);
  if ($info_b eq ".") { # Single line variant in VCF
    ($leftgenenames,$rightgenenames) = getfields($info_a,"left_Gene","right_Gene");
    ($leftintrons,$rightintrons) = getfields($info_a,"left_Intron","right_Intron") if $truncationtypes{$svtype};
    ($leftintrongenenames,$rightintrongenenames) = getfields($info_a,"left_IntronGene","right_IntronGene");
  } else { # Multiline variant in VCF (possibly only one line of variant present)
    $leftgenenames = getfields($info_a,"Gene");
    $leftintrongenenames = getfields($info_a,"IntronGene");
    $leftintrons = getfields($info_a,"Intron") if $truncationtypes{$svtype};
    $rightintrons = getfields($info_b,"Intron") if $truncationtypes{$svtype};
    $rightgenenames = getfields($info_b,"Gene");
    $rightintrongenenames = getfields($info_b,"IntronGene");
  }
  my @leftgenenames = split(/\|/,$leftgenenames);
  my @rightgenenames = split(/\|/,$rightgenenames);
  
  my @probleft = split(/,/,$probleft) if $probleft;
  my @probright = split(/,/,$probright) if $probright;

  my %scores = (); # Interval => List of scores by op; e.g. (LEFT => (MAXLEFT, SUMLEFT, TOP100LEFT, MEANLEFT), RIGHT => (MAXRIGHT, SUMRIGHT, TOP100RIGHT, MEANRIGHT))

  #if ($svtype eq "INS") {
  #  $scores{"LEFT"} = cscoreop($caddfile, "", $ops, $leftchrom, $leftstart-10, $leftstart, "", $topn);
  #  $scores{"RIGHT"} = cscoreop($caddfile, "", $ops, $leftchrom, $rightstop, $rightstop+10, "", $topn);
  #} else {
  $scores{"LEFT"} = cscoreop($caddfile, $ops, $leftchrom, $leftstart, $leftstop, \@probleft);
  $scores{"RIGHT"} = cscoreop($caddfile, $ops, $rightchrom, $rightstart, $rightstop, \@probright);
  #}

  if ($svtype eq "DEL" || $svtype eq "DUP" || $svtype eq "CNV") {
    my ($pos,$end) = getfields($info_a,"POS","END");
    if ($rightstop - $leftstart > 1000000) {
      $scores{"SPAN"} = (100) x @ops;
    } else {
      $scores{"SPAN"} = cscoreop($caddfile, $ops, $leftchrom, $pos, $end, "");
    }
  }

  # Calculate truncation scores
  if (exists $truncationtypes{$svtype}) {
    my @leftintrons = split(/\|/,$leftintrons);
    my @leftintrongenenames = split(/\|/,$leftintrongenenames);
    my @rightintrons = split(/\|/,$rightintrons);
    my @rightintrongenenames = split(/\|/,$rightintrongenenames);
    my %leftintrons = map {$leftintrons[$_] => $leftintrongenenames[$_]} (0..$#leftintrons); # @leftintrons and @leftintrongenes should have the same number of elements if vcfanno is working as it should
    my %rightintrons = map {$rightintrons[$_] => $rightintrongenenames[$_]} (0..$#rightintrons);
    foreach my $intron (keys %leftintrons) { # Cancel out introns hit by both right and left breakends
      if (exists $rightintrons{$intron}) {
	delete $rightintrons{$intron};
	delete $leftintrons{$intron};
      }
    }
    @leftintrongenenames = uniq(values(%leftintrons));
    @rightintrongenenames = uniq(values(%rightintrons));
    if (@leftintrongenenames) {
      my $leftscore = truncationscore($leftchrom, $leftstart, $leftstop, \@leftintrongenenames, \%genes, $caddfile, $ops, \%operations);
      $scores{"LTRUNC"} = $leftscore if $leftscore;
    }
    if (@rightintrongenenames) {
      my $rightscore = truncationscore($rightchrom, $rightstart, $rightstop, \@rightintrongenenames, \%genes, $caddfile, $ops, \%operations);
      $scores{"RTRUNC"} = $rightscore if $rightscore;
    }
  }

  # This is an ugly loop which transposes %scores so that the keys are operations, not intervals
  my %scoresbyop = ();
  foreach my $interval (sort {$intervals{$a} <=> $intervals{$b}} keys %scores) { # LEFT, RIGHT, (SPAN, LTRUNC, RTRUNC)
    foreach my $op (@ops) { # MAX, SUM, TOP$topn, MEAN
      push @{$scoresbyop{$op}}, $scores{$interval}->[$operations{$op}];
    }
  }

# Calculate maxes and add to info, replacing existing SVSCORE fields if they exist
  my %maxscores = ();
  foreach my $op (keys %scoresbyop) {
    $info_a = replaceoraddfield($op, "", $info_a, \%operations, \%scoresbyop);
    $info_b = replaceoraddfield($op, "", $info_b, \%operations, \%scoresbyop);
    if ($verbose) {
      foreach my $interval (keys %scores) {
	$info_a = replaceoraddfield($op, $interval, $info_a, \%operations, \%scores);
	$info_b = replaceoraddfield($op, $interval, $info_b, \%operations, \%scores);
      }
    }
  }

  $splitline[12] = $info_a;
  $splitline[13] = $info_b;
  print OUT join("\t",@splitline) . "\n";

  print STDERR $.,", " if $debug;
}

close IN;
close OUT;

# Convert back to vcf, sort, and add header
system("svtools bedpetovcf -b $bedpeout > $vcfout");
system("grep \"^#\" $vcfout > $vcfout.header");
print `grep -v "^#" $vcfout | sort -k1,1V -k2,2n | cat $vcfout.header -`;
unlink "$vcfout.header";

# Clean up
unless ($debug) {
  unlink $preprocessedfile;
  unlink $vcfout;
  unlink $bedpeout;
  unlink $inputfile if $compressed || $prefix eq "stdin";
  
  if ($alteredgenefile) {
    if ($genefile =~ /\.gz$/) {
      unlink $uncompressedgenefile;
    } else {
      unlink $compressedgenefile;
    }
  }

  if ($alteredexonfile) {
    if ($exonfile =~ /\.gz$/) {
      unlink $uncompressedexonfile;
    }
  }

  # Delete svscoretmp if empty
  opendir(DIR,"svscoretmp");
  my @dir = readdir(DIR);
  my $dircount = @dir;
  foreach (@dir) {
    $dircount-- if ($_ eq '.' || $_ eq '..');
  }
  if ($dircount == 0) {
    if(system("rmdir svscoretmp")) {
      warn "Could not delete svscoretmp: $!";
    }
  }
  closedir(DIR);
}

sub cscoreop { # Apply operation(s) specified in $ops to C scores within a given region using CADD data. In VCF coordinates, $start is the base preceding the first possible breakpoint, and $stop is the base preceding the last possible breakpoint. In BED, $start is the first possible breakpoint, and $stop is the final possible breakpoint
  my ($filename, $ops, $chrom, $start, $stop, $prpos) = @_;
  my @ops = split(/,/,$ops);
#  pop @{$prpos} if $weight; ## TEMPORARY FIX - remove final element of @probdist to make length of @probdist equal to # of bases in interval
  my $weight = ($ops =~ /WEIGHTED/);
  my @prpos = @{$prpos} if $weight;
  if ($weight && !scalar @prpos) { ## Don't calculate scores if other scores are being weighted but this variant does not have a probability distribution
    my @res = (-1) x @ops;
    return \@res;
  }
  my %probdist = () if $weight; # Hash of pdf for each position 
  if ($weight) {
    foreach my $i ($start..$stop) { # Populate %probdist
      $probdist{$i} = $prpos[$i-$start];
    }
  }
  
  my (%bptscores,$res) = (); # %bptscores = {BEDcoordinate => Possiblebreakpointscore}
  my $stopinc = $stop+1; ## Increment end coordinate for tabix so that we capture the CADD score for the base following the interval to allow for calculation of a score for the final possible breakpoint
  my $tabixoutput = `tabix $filename $chrom:$start-$stopinc`;
  my @tabixoutputlines = split(/\n/,$tabixoutput);
#  print STDERR "$chrom:$start-$stop\n"; ## DEBUG
#  print STDERR Dumper(keys %probdist),"\n"; ## DEBUG
  
  my %basescores = (); # Hash from VCF position to list of scores for each position
  foreach my $line (@tabixoutputlines) { # Populate %basescores with tabix output
    my @split = split(/\s+/,$line);
    push @{$basescores{$split[1]}}, $split[5];
  }
  foreach my $pos (sort {$a <=> $b} keys %basescores) { # Replace values in %basescores with max score at each VCF position
#    push @basemaxscores, max(@{$basescores{$pos}});
    $basescores{$pos} = max(@{$basescores{$pos}});
    if (exists $basescores{$pos-1}) { # Calculate all scores for possible breakpoints by averaging the base scores of the two flanking bases of the possible breakpoint. Place scores in %bptscores
      $bptscores{$pos-1} = ($basescores{$pos-1} + $basescores{$pos}) / 2;
      delete $basescores{$pos-1}; # Save some memory - the basescores for all bases before $pos are now useless
    }
  }

#  print STDERR "\n$chrom:$start-$stop"; ## DEBUG
  unless (@tabixoutputlines && %bptscores) { # Short circuit if interval does not have enough base scores to calculate breakpoint scores (i.e. there are no 2 consecutive bases with scores in the interval)
    my @res = (-1) x @ops;
    return \@res; 
  }
#  if ($weight && scalar(keys %bptscores) != @probdist) { # If weighting with probdist, but C scores are not available for the entire interval, trim %probdist to contain only those BED positions which have possible breakpoint scores (i.e. base scores for two flanking bases are available)
##    my @newprobdist = ();
##    foreach my $pos (sort {$a <=> $b} keys %basescores) { # Get probability of positions in interval that have C scores, collect in @newprobdist
##      push @newprobdist,$probdist{$pos};
##    }
##    @probdist = @newprobdist; # Replace @probdist
#    for my $pos (keys %probdist) { # Get rid of probabilities for any possible breakpoint position without a score
#      delete $probdist{$pos} unless exists $bptscores{$pos};
#    }
#  }

#  print STDERR "Scores: ",join(", ",@basemaxscores),"\n"; ## DEBUG
#  print STDERR "PRPOS: ",join(", ",@probdist),"\n" if $weight; ## DEBUG

  my (@bptscores,@probdist,@weightedbptscores) = ();
  foreach my $pos (sort {$a <=> $b} keys %bptscores) { # Collapse %bptscores and %probdist into arrays, getting rid of positions in %probdist with no corresponding possible breakpoint scores
    push @bptscores, $bptscores{$pos};
    push @probdist, $probdist{$pos} if $weight;
  }

  if ($weight) { # Rescale probability distribution to add up to 1 and weight @bptscores
    my $sumprobs = sum(@probdist);
    unless ($sumprobs == 1) {
      foreach my $i (0..$#probdist) {
	$probdist[$i] = $probdist[$i] / $sumprobs;
      }
    }
    @weightedbptscores = pairwise {$a * $b} @bptscores, @probdist;
  }

#  print STDERR "\t\t\t@bptscores\n"; ## DEBUG

#    print STDERR "Probdist after normalization: ",join(',',@probdist),"\n"; ## DEBUG
#  print STDERR "Newscores: ",join(", ",@bptscores),"\n"; ## DEBUG
  foreach my $op (@ops) {
    my $scoresref = ($op =~ /WEIGHTED/ ? \@weightedbptscores : \@bptscores); # Use a reference to avoid copying arrays
    my $newscore;

    if ($op =~ /^TOP(\d+)/) {
      my $topn = min($1, scalar @{$scoresref});
      my @topn = (sort {$b <=> $a} @{$scoresref})[0..$topn-1];
      $scoresref = \@topn;
    }

    my @scores = @{$scoresref};
    if ($op eq "MAX") {
      $newscore = nearest(0.001,max(@scores));
    } elsif ($op eq "SUM" || ($op =~ /WEIGHTED$/ && $prpos)) { # Compute sum of @scores if $op is SUM, or if the op is weighted and we're on a breakend, in which case, @scores is @weightedbptscores (or a subset of it in the TOP case), so the sum of @scores is the weighted mean of @bptscores (or of the subset, if $op begins with TOP)
      $newscore = nearest(0.001,sum(@scores));
    } elsif ($op eq "MEAN" || $op =~ /^TOP\d+/ || ($op =~ /WEIGHTED$/ && !$prpos)) { # Compute mean of @scores if $op is MEAN or TOP, or if $op is MEANWEIGHTED or TOP\d+WEIGHTED and we're on the span/truncation score
      $newscore = nearest(0.001,sum(@scores)/scalar(@scores));
    } else {
      die "Error: Unrecognized operation: $op"
    }
    push @{$res}, $newscore;
  }
  return $res;
}

sub getfields { # Parse info field of VCF line, getting fields specified in @_. $_[0] must be the info field itself. Returns list of field values if more than one field is being requested; otherwise, returns a scalar value representing the requested field
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

sub truncationscore { # Calculate truncation score based on the coordinates of a given breakend and the names of the genes it hits
  my ($chrom, $start, $stop, $introngenesref, $genesref, $caddfile, $ops, $operationsref) = @_;
  my @ops = split (/,/,$ops);
  my @introngenes = @{$introngenesref};
  return "" unless @introngenes;
  my %operations = %{$operationsref};
  my %genes = %{$genesref};
  my %truncationscores;
  foreach my $gene (@introngenes) {
    my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$chrom}}[0..2];	
    my $cscoreopres;
    if ($genestrand eq '+') {
      $cscoreopres = cscoreop($caddfile, "", $ops, $chrom, max($genestart,$start),$genestop, ""); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
    } else {
      $cscoreopres = cscoreop($caddfile, "", $ops, $chrom, $genestart,min($genestop,$stop), ""); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but none of the supported operations are order-dependent)
    }
    foreach my $op (keys %operations) {
      push @{$truncationscores{$op}}, $cscoreopres->[$operations{$op}];
    }
  }
  
  my $res = [];
  foreach my $op (sort {$operations{$a} <=> $operations{$b}} keys %truncationscores) {
    push @{$res}, max(@{$truncationscores{$op}});
  }

  return $res;
}

sub replaceoraddfield {
  my ($op, $interval, $info, $operationsref, $scoresref) = @_; ## $scoresref points to either %scores or %scoresbyop, depending on whether $interval is defined
  my %scores = %{$scoresref};
  my %operations = %{$operationsref};
  my $field = "SVSCORE$op" . ($interval ? "_$interval" : "" );
  my $newscore = ($interval ? $scores{$interval}->[$operations{$op}] : max(@{$scores{$op}}));
  if ($info =~ /[\t;]$field=/) { ## Replace existing field
    $info =~ s/$field=[^\t;]*/$field=$newscore/;
  } elsif ($info ne "MISSING" && $info ne ".") { ## Info is present but $field needs to be added
    $info .= (";$field=$newscore");
  }
  return $info;
}

sub main::HELP_MESSAGE() {
  print STDERR "usage: ./svscore.pl [-dsvw] [-o op] [-t topnumber] [-g genefile] [-m geneannotationcolumn] [-p genestrandcolumn] [-e exonfile] [-c caddfile] -i vcf
    -i	      Input VCF file. May be bgzip compressed (ending in .vcf.gz). Use \"-i stdin\" if using standard input
    -d	      Debug mode, keeps intermediate and supporting files, displays progress
    -v	      Verbose mode - show all calculated scores (left/right/span/ltrunc/rtrunc, as appropriate)
    -o	      Comma-separated list of operations to perform on CADD score intervals (must be some combination of sum, max, mean, meanweighted, top\\d, or top\\dweighted - defaults to top10weighted)
    -g	      Points to gene BED file (refGene.genes.b37.bed)
    -m	      Column number for annotation in gene BED file to be added to VCF (4)
    -p	      Column number for strand in gene BED file (5)
    -e	      Points to exon BED file (refGene.exons.b37.bed)
    -c	      Points to whole_genome_SNVs.tsv.gz (defaults to current directory)

    --help    Display this message
    --version Display version\n"
}

sub main::VERSION_MESSAGE() {
  print "SVScore version $main::VERSION\n";
}
