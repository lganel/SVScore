#!/usr/bin/perl -w 
## Author: Liron Ganel
## Laboratory of Ira Hall, McDonnell Genome Institute
## Washington University in St. Louis
## Version 0.4

use strict;
use Getopt::Std;
use List::Util qw(max min sum);
use List::MoreUtils qw(pairwise uniq);
use Time::HiRes qw(gettimeofday);
use Math::Round qw(nearest);

sub main::HELP_MESSAGE(); # Declaration to make perl happy
$Getopt::Std::STANDARD_HELP_VERSION = 1; # Make --help and --version flags halt execution
$main::VERSION = '0.4';

my %possibleoperations = ("MAX", 0, "SUM", 1, "TOP", 2, "MEAN", 3); # Hash of supported operations
my %types = map {$_ => 1} ("DEL", "DUP", "INV", "BND", "TRX", "CNV", "MEI", "INS"); # Hash of supported svtypes
#my %onelinetruncationtypes = map {$_ => 1} ("INV"); # Hash of svtypes represented on one line for which truncation scores are calculated
#my %twolinetruncationtypes = map {$_ => 1} ("TRX", "INS"); # Hash of svtypes represented on two lines for which truncation scores are calculated
my %truncationtypes = ("INV","TRX","INS","DEL"); # All svtypes for which truncation scores are calculated
my %intervals = ("LEFT", 0, "RIGHT", 1, "SPAN", 2, "LTRUNC", 3, "RTRUNC", 4); # Hash of supported intervals

my %options = ();
getopts('dswvc:g:e:m:o:t:p:',\%options);

# Parse command line options, set variables, check input parameters
my $debug = defined $options{'d'};
my $support = defined $options{'s'};
my $weight = defined $options{'w'};
my $verbose = defined $options{'v'};

my $pipedinput = !-t STDIN; ## Tells whether input is coming to STDIN through a pipe
&main::HELP_MESSAGE() && die unless (defined $ARGV[0] || $pipedinput || $support);

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
my $ops = (defined $options{'o'} ? $options{'o'} : 'ALL');
my $topn = (defined $options{'t'} ? $options{'t'} : 100);
warn "Unrecognized operation specified: $ops\n" && die main::HELP_MESSAGE() unless ($ops eq "ALL" || defined $possibleoperations{$ops});
my %operations = ($ops eq 'ALL' ? %possibleoperations : ($ops => 0)); # Hash of indices for each operation within lists in the values of %scores
if ($ops eq 'ALL' || $ops eq 'TOP') {
  $operations{"TOP$topn"} = $operations{"TOP"};
  delete $operations{"TOP"};
}
my $compressed = ($ARGV[0] =~ /\.gz$/) if defined $ARGV[0];
my ($uncompressedgenefile, $compressedgenefile, $uncompressedexonfile, $alteredgenefile, $alteredexonfile, $headerfile, $sortedfile, $preprocessedfile, $bedpeout, $vcfout);

# Set up all necessary preprocessing to be taken care of before analysis can begin. This includes decompression, annotation using vcfanno, and generation of intron/exon/gene files, whichever are necessary. May be a little slower than necessary in certain situations because some arguments are supplied by piping cat output rather than supplying filenames directly.
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
my ($inputfile, $prefix, $time);
if ($pipedinput) {
  open(TEMP,">svscoretmp/pipedinput.tmp.vcf") || die "Could not open svscoretmp/pipedinput.tmp.vcf for writing; $!";
  print TEMP <STDIN>;
  close TEMP;
  $inputfile = "svscoretmp/pipedinput.tmp.vcf";
} else {
  $inputfile = $ARGV[0];
}
if (defined $inputfile) {
  print STDERR "Preparing preprocessing command\n" if $debug;
  ($prefix) = ($inputfile =~ /^(?:.*\/)?(.*)\.vcf(?:\.gz)?$/);

  # Tag intermediate files with timestamp to avoid collisions
  $time = gettimeofday();
  $preprocessedfile = "svscoretmp/$prefix$time.preprocess.bedpe";
  my $sortedfile = "svscoretmp/$prefix$time.sort.vcf.gz";
  my $reorderout = "svscoretmp/$prefix$time.reorderheaderout.vcf";
  if ($compressed) {
    $inputfile = "$prefix.vcf";
    system("gunzip -c $ARGV[0] > $prefix.vcf") && die "Could not unzip $ARGV[0]: $!";
  }
  my $preprocess = "awk '\$0~\"^#\" {print \$0; next } { print \$0 | \"sort -k1,1V -k2,2n\" }' $inputfile | bgzip -c > $sortedfile; tabix -p vcf $sortedfile; vcfanno -ends conf.toml $sortedfile | perl reorderheader.pl stdin $inputfile > $reorderout"; # Sort, annotate, reorder header
  my $preprocess2 = "vcftobedpe -i $reorderout > $preprocessedfile; rm -f $sortedfile $sortedfile.tbi $reorderout";
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
  open(OUT, "> $bedpeout") || die "Could not open output file: $!" unless $support;

  # Update header
  unless ($support) {
    open(HEADER, "grep \"^#\" $preprocessedfile |") || die "Error grabbing header: $!";
    my @newheader = ();
    my @oldheader = <HEADER>;
    close HEADER;
    die "Header indicates that SVScore has already been run on this file. Please remove these annotations and header lines to avoid confusion between old and new scores" if grep {/^##INFO=<ID=SVSCORE/} @oldheader;
    my $headerline;
    while(($headerline = (shift @oldheader)) !~ /^##INFO/) {
      push @newheader, $headerline;
    }
    unshift @oldheader, $headerline; # Return first info line to top of stack
    while(($headerline = (shift @oldheader)) =~ /^##INFO/) {
      push @newheader, $headerline;
    }
    unshift @oldheader, $headerline; # Return first format line to top of stack

    if ($ops eq 'MAX' || $ops eq 'ALL') {
      push @newheader, "##INFO=<ID=SVSCOREMAX_LEFT,Number=1,Type=Float,Description=\"Maximum C score in left breakend of structural variant" . ($weight ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMAX_RIGHT,Number=1,Type=Float,Description=\"Maximum C score in right breakend of structural variant" . ($weight ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMAX_SPAN,Number=1,Type=Float,Description=\"Maximum C score in outer span of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMAX_LTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of left breakend to end of truncated gene\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMAX_RTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of right breakend to end of truncated gene\">\n";
    }
    if ($ops eq 'SUM' || $ops eq 'ALL') {
      push @newheader, "##INFO=<ID=SVSCORESUM_LEFT,Number=1,Type=Float,Description=\"Sum of C scores in left breakend of structural variant" . ($weight ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORESUM_RIGHT,Number=1,Type=Float,Description=\"Sum of C scores in right breakend of structural variant" . ($weight ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORESUM_SPAN,Number=1,Type=Float,Description=\"Sum of C scores in outer span of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCORESUM_LTRUNC,Number=1,Type=Float,Description=\"Sum of C scores from beginning of left breakend to end of truncated gene\">\n";
      push @newheader, "##INFO=<ID=SVSCORESUM_RTRUNC,Number=1,Type=Float,Description=\"Sum of C scores from beginning of right breakend to end of truncated gene\">\n";
    }
    if ($ops eq 'TOP' || $ops eq 'ALL') {
      push @newheader, "##INFO=<ID=SVSCORETOP${topn}_LEFT,Number=1,Type=Float,Description=\"Mean of top $topn C scores in left breakend of structural variant" . ($weight ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${topn}_RIGHT,Number=1,Type=Float,Description=\"Mean of top $topn C scores in right breakend of structural variant" . ($weight ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${topn}_SPAN,Number=1,Type=Float,Description=\"Mean of top $topn C scores in outer span of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${topn}_LTRUNC,Number=1,Type=Float,Description=\"Mean of top $topn C scores from beginning of left breakend to end of truncated gene\">\n";
      push @newheader, "##INFO=<ID=SVSCORETOP${topn}_RTRUNC,Number=1,Type=Float,Description=\"Mean of top $topn C scores from beginning of right breakend to end of truncated gene\">\n";
    }
    if ($ops eq 'MEAN' || $ops eq 'ALL') {
      push @newheader, "##INFO=<ID=SVSCOREMEAN_LEFT,Number=1,Type=Float,Description=\"Mean of C scores in left breakend of structural variant" . ($weight ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_RIGHT,Number=1,Type=Float,Description=\"Mean of C scores in right breakend of structural variant" . ($weight ? ", weighted by probability distribution" : "") . "\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_SPAN,Number=1,Type=Float,Description=\"Mean of C scores in outer span of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_LTRUNC,Number=1,Type=Float,Description=\"Mean of C scores from beginning of left breakend to end of truncated gene\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMEAN_RTRUNC,Number=1,Type=Float,Description=\"Mean of C scores from beginning of right breakend to end of truncated gene\">\n";
    }
    push @newheader, @oldheader;
    foreach (@newheader) {
      print OUT;
    }
    if ($weight && !(grep {/^##INFO=<ID=PRPOS,/} @newheader)) {
      $weight = "";
      warn "PRPOS not found in header - SVScore will not weight CADD scores\n";
    }
  }

}

if ($support) {
  if ((defined $ARGV[0] || $pipedinput) && !$debug) {
    unlink $preprocessedfile;
  }
  die;
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
  $leftstart++; ## Make intervals inclusive in VCF (1-based)
  $rightstart++;

  $svtype = substr($svtype, 0, 3);
  unless (exists $types{$svtype}) {
    warn "Unrecognized SVTYPE $svtype at line ",$.," of preprocessed VCF file\n";
    print OUT $inputline;
    next;
  }
  my ($probleft,$probright) = getfields($info_a,"PRPOS","PREND") if $weight;

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
  my $localweight = $weight && $probleft;
  
  my @probleft = split(/,/,$probleft) if $probleft;
  my @probright = split(/,/,$probright) if $probright;

  my %scores = (); # Interval => List of scores by op; e.g. (LEFT => (MAXLEFT, SUMLEFT, TOP100LEFT, MEANLEFT), RIGHT => (MAXRIGHT, SUMRIGHT, TOP100RIGHT, MEANRIGHT))

  #if ($svtype eq "INS") {
  #  $scores{"LEFT"} = cscoreop($caddfile, "", $ops, $leftchrom, $leftstart-10, $leftstart, "", $topn);
  #  $scores{"RIGHT"} = cscoreop($caddfile, "", $ops, $leftchrom, $rightstop, $rightstop+10, "", $topn);
  #} else {
    $scores{"LEFT"} = cscoreop($caddfile, $localweight, $ops, $leftchrom, $leftstart, $leftstop, \@probleft, $topn);
    $scores{"RIGHT"} = cscoreop($caddfile, $localweight, $ops, $rightchrom, $rightstart, $rightstop, \@probright, $topn);
  #}

  if ($svtype eq "DEL" || $svtype eq "DUP" || $svtype eq "CNV") {
    my ($pos,$end) = getfields($info_a,"POS","END");
    if ($rightstop - $leftstart > 1000000) {
      $scores{"SPAN"} = ($ops eq "ALL" ? [100, 100, 100, 100] : [100]);
    } else {
      $scores{"SPAN"} = cscoreop($caddfile, "", $ops, $leftchrom, $pos, $end, "", $topn);
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
      my $leftscore = truncationscore($leftchrom, $leftstart, $leftstop, \@leftintrongenenames, \%genes, $caddfile, $ops, $topn, \%operations);
      $scores{"LTRUNC"} = $leftscore if $leftscore;
    }
    if (@rightintrongenenames) {
      my $rightscore = truncationscore($rightchrom, $rightstart, $rightstop, \@rightintrongenenames, \%genes, $caddfile, $ops, $topn, \%operations);
      $scores{"RTRUNC"} = $rightscore if $rightscore;
    }
  }

  # This is an ugly loop which transposes %scores so that the keys are operations, not intervals
  my %scoresbyop = ();
  foreach my $interval (sort {$intervals{$a} <=> $intervals{$b}} keys %scores) { # LEFT, RIGHT, (SPAN, LTRUNC, RTRUNC)
    foreach my $op (keys %operations) { # MAX, SUM, TOP$topn, MEAN
      push @{$scoresbyop{$op}}, $scores{$interval}->[$operations{$op}];
    }
  }

# Calculate maxes and add to info
  my %maxscores = ();
  foreach my $op (keys %scoresbyop) {
    $info_a .= (";SVSCORE${op}=" . max(@{$scoresbyop{$op}})) unless $info_a eq "MISSING";
    $info_b .= (";SVSCORE${op}=" . max(@{$scoresbyop{$op}})) unless $info_b eq "." || $info_b eq "MISSING";
    if ($verbose) {
      foreach my $interval (keys %scores) {
	$info_a .= ";SVSCORE${op}_$interval=$scores{$interval}->[$operations{$op}]" unless $info_a eq "MISSING";
	$info_b .= ";SVSCORE${op}_$interval=$scores{$interval}->[$operations{$op}]" unless ($info_b eq "." || $info_b eq "MISSING");
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
system("bedpetovcf -b $bedpeout > $vcfout");
system("grep \"^#\" $vcfout > $vcfout.header");
print `grep -v "^#" $vcfout | sort -k1,1V -k2,2n | cat $vcfout.header -`;
unlink "$vcfout.header";

# Clean up
unless ($debug) {
  unlink $preprocessedfile;
  unlink $vcfout;
  unlink $bedpeout;
  unlink $inputfile if $compressed || $pipedinput;
  
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

sub cscoreop { # Apply operation(s) specified in $ops to C scores within a given region using CADD data. [$start,$stop] represents an inclusive interval in VCF coordinates
  my ($filename, $weight, $ops, $chrom, $start, $stop, $probdist, $topn) = @_;
#  print STDERR join(",",$start..$stop),"\n"; ## DEBUG
  my @probdist = @{$probdist} if $weight;
  my (@scores,$res) = ();
  my $tabixoutput = `tabix $filename $chrom:$start-$stop`;
  my @tabixoutputlines = split(/\n/,$tabixoutput);
  return ($ops eq 'ALL' ? [-1, -1, -1, -1] : [-1]) unless @tabixoutputlines; # Short circuit if interval has no C scores
  
  my %allscores = (); # Hash from position to list of scores
  foreach my $line (@tabixoutputlines) { # Populate hash
    my @split = split(/\s+/,$line);
    push @{$allscores{$split[1]}}, $split[5];
  }
  foreach my $pos (sort {$a <=> $b}  keys %allscores) { # Populate @scores
    push @scores, max(@{$allscores{$pos}});
  }

#  print STDERR join(", ",@scores),"\n"; ## DEBUG
#  print STDERR join(", ",@probdist),"\n"; ## DEBUG
  @scores = pairwise {$a * $b}	@scores, @probdist if $weight;
#  print STDERR join(", ",@scores),"\n"; ## DEBUG
  if ($ops eq 'MAX') {
    $res = [max(@scores)];
  } elsif ($ops eq 'SUM') {
    $res = [nearest(0.001,sum(@scores))];
  } elsif ($ops=~/^TOP/) {
    $topn = min($topn, scalar @scores);
    my @topn = (sort {$b <=> $a} @scores)[0..$topn-1];
    $res = [nearest(0.001,sum(@topn) / scalar(@topn))];
  } elsif ($ops eq 'MEAN') {
    $res = [nearest(0.001,sum(@scores)/scalar(@scores))];
  } else {
    $topn = min($topn, scalar @scores);
    my @topn = (sort {$b <=> $a} @scores)[0..$topn-1];
    $res = [max(@scores), nearest(0.001,sum(@scores)), nearest(0.001,sum(@topn) / scalar(@topn)), nearest(0.001,sum(@scores)/scalar(@scores))];
  }
#  print STDERR "\n"; ## DEBUG
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
  my ($chrom, $start, $stop, $introngenesref, $genesref, $caddfile, $ops, $topn, $operationsref) = @_;
  my @introngenes = @{$introngenesref};
  return "" unless @introngenes;
  my %operations = %{$operationsref};
  my %genes = %{$genesref};
  my %truncationscores;
  foreach my $gene (@introngenes) {
    my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$chrom}}[0..2];	
    my $cscoreopres;
    if ($genestrand eq '+') {
      $cscoreopres = cscoreop($caddfile, "", $ops, $chrom, max($genestart,$start),$genestop, "", $topn); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
    } else {
      $cscoreopres = cscoreop($caddfile, "", $ops, $chrom, $genestart,min($genestop,$stop), "", $topn); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but none of the supported operations are order-dependent)
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

sub main::HELP_MESSAGE() {
  print STDERR "usage: ./svscore.pl [-dsvw] [-o op] [-t topnumber] [-g genefile] [-m geneannotationcolumn] [-p genestrandcolumn] [-e exonfile] [-c caddfile] vcf
    -d	      Debug mode, keeps intermediate and supporting files, displays progress
    -s	      Create/download supporting files and quit
    -v	      Verbose mode - show all calculated scores (left/right/span/ltrunc/rtrunc, as appropriate)
    -w	      Weight CADD scores in breakends by probability distribution (requires PRPOS/PREND in INFO field)
    -o	      Specify operation to perform on CADD scores (must be sum, max, top, or all - defaults to all)
    -t	      Number of bases for TOP to take mean of (under -o top or -o all) (100)
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
