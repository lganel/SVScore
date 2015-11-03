#!/usr/bin/perl -w 
## Author: Liron Ganel
## Laboratory of Ira Hall, McDonnell Genome Institute
## Washington University in St. Louis
## Version 0.4

use strict;
use Getopt::Std;
use List::Util qw(max min sum);
use List::MoreUtils qw(pairwise);
use Time::HiRes qw(gettimeofday);
use Math::Round qw(nearest);

sub main::HELP_MESSAGE(); # Declaration to make perl happy
$Getopt::Std::STANDARD_HELP_VERSION = 1; # Make --help and --version flags halt execution
$main::VERSION = '0.4';

my %possibleoperations = ("MAX", 0, "SUM", 1, "TOP", 2, "MEAN", 3); # Hash of supported operations
my %types = map {$_ => 1} ("DEL", "DUP", "INV", "BND", "TRX", "INS", "CNV", "MEI"); # Hash of supported svtypes

my %options = ();
getopts('dswvc:g:e:m:o:t:',\%options);

# Parse command line options, set variables
my $debug = defined $options{'d'};
my $support = defined $options{'s'};
my $weight = defined $options{'w'};
my $verbose = defined $options{'v'};

&main::HELP_MESSAGE() && die unless defined $ARGV[0] || $support;

my $caddfile = (defined $options{'c'} ? $options{'c'} : 'whole_genome_SNVs.tsv.gz');
die "Could not find $caddfile" unless -s $caddfile;
my $genefile = (defined $options{'g'} ? $options{'g'} : 'refGene.genes.b37.bed');
my $geneanncolumn = (defined $options{'m'} && defined $options{'g'} ? $options{'m'} : 4);
warn "Gene annotation column provided without nonstandard gene annotation file - defaulting to standard gene annotation column (4)" if defined $options{'m'} && !defined $options{'g'};
my $exonfile = (defined $options{'e'} ? $options{'e'} : 'refGene.exons.b37.bed');
#my $exonanncolumn = (defined $options{'n'} && defined $options{'e'} ? $options{'n'} : (defined $options{'e'} ? 4 : 5));
#warn "Exon annotation column provided without nonstandard exon annotation file - defaulting to standard exon annotation column (5)" if defined $options{'n'} && !defined $options{'e'};
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
my ($headerfile, $preprocessedfile, $bedpeout, $vcfout);

##TODO PRIORITY 2: Enable piping input through STDIN - use an option to specify input file rather than @ARGV

# Set up all necessary preprocessing to be taken care of before analysis can begin. This includes decompression, annotation using vcfanno, and generation of intron/exon/gene files, whichever are necessary. May be a little slower than necessary in certain situations because some arguments are supplied by piping cat output rather than supplying filenames directly.
if ($exonfile eq 'refGene.exons.b37.bed' && !-s $exonfile) { # Generate exon file if necessary
  print STDERR "Generating exon file\n" if $debug;
  system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); n=int(\$9); split(\$10,start,\",\");split(\$11,end,\",\"); for(i=1;i<=n;++i) {print \$3,start[i],end[i],\$2\".\"i,\$13,\$2; } }' OFS=\"\t\" | sort -k 1,1n -k 2,2n | uniq > refGene.exons.b37.bed");
} elsif ($exonfile ne 'refGene.exons.b37.bed' && !-s $exonfile) {
  die "$exonfile not found or empty!";
}

if ($genefile eq 'refGene.genes.b37.bed' && !-s $genefile) { # Generate gene file if necessary
  print STDERR "Generating gene file\n" if $debug;
  system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); print \$3,\$5,\$6,\$13,\$4}' OFS=\"\\t\" | sort -k 1,1n -k 2,2n | uniq > refGene.genes.b37.bed");
} elsif ($genefile ne 'refGene.genes.b37.bed' && !-s $genefile) {
  die "$genefile not found or empty!";
}

unless (-s 'introns.bed') { # Generate intron file if necessary - add column with unique intron ID equal to line number (assuming introns.bed has no header line) and sort
  print STDERR "Generating intron file\n" if $debug;
  system("bedtools subtract -a $genefile -b $exonfile | sort -u -k 1,1n -k 2,2n | awk '{print \$0 \"\t\" NR}' > introns.bed");
}

my $intronnumcolumn = `head -n 1 introns.bed | awk '{print NF}'`; # Figure out which column has the intron number
chomp $intronnumcolumn;

# Write conf.toml file
print STDERR "Writing config file\n" if $debug;
open(CONFIG, "> conf.toml") || die "Could not open conf.toml: $!";
print CONFIG "[[annotation]]\nfile=\"$genefile\"\nnames=[\"Gene\"]\ncolumns=[4]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"introns.bed\"\nnames=[\"Intron\"]\ncolumns=[$intronnumcolumn]\nops=[\"uniq\"]";
close CONFIG;

# Create first preprocessing command - annotation is done without normalization because REF and ALT nucleotides are not included in VCFs describing SVs
mkdir "svscoretmp" unless -d "svscoretmp";
my ($prefix, $time);
if (defined $ARGV[0]) {
  print STDERR "Preparing preprocessing command\n" if $debug;
  ($prefix) = ($ARGV[0] =~ /^(?:.*\/)?(.*)\.vcf(?:\.gz)?$/);

  # Tag intermediate files with timestamp to avoid collisions
  $time = gettimeofday();
  $preprocessedfile = "svscoretmp/$prefix.preprocess$time.bedpe";
  #$annfile = "svscoretmp/$prefix.ann$time.bedpe";
#  $headerfile = "svscoretmp/${prefix}header$time";
  my $preprocess = ($compressed ? "z": "") . "cat $ARGV[0] | awk '\$0~\"^#\" {print \$0; next } { print \$0 | \"sort -k1,1n -k2,2n\" }' | vcfanno -ends conf.toml - | vcftobedpe > $preprocessedfile"; #; grep '^#' $preprocessedfile > $headerfile"; # Sort, annotate, convert to BEDPE, grab header
  print STDERR "Preprocessing command: $preprocess\n" if $debug;
  if (system($preprocess)) {
#    unlink $annfile unless $debug;
    unless ($debug) {
#      unlink $headerfile;
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
#  open(HEADER, "> $headerfile") || die "Could not open $headerfile: $!";
    foreach (@newheader) {
      print OUT;
    }
    if ($weight && !(grep {/^##INFO=<ID=PRPOS,/} @newheader)) {
      $weight = "";
      warn "PRPOS not found in header - SVScore will not weight CADD scores\n";
    }
  }

  # Create and execute second preprocessing command
#  my $preprocess2 = "grep -v '^#' $annfile | sort -k 3,3 | cat $headerfile - > $preprocessedfile; rm -f $headerfile"; # Sort by ID, add new header, clean up
#  print STDERR "Preprocessing command 2: $preprocess2\n" if $debug;
#  if (system($preprocess2)) {
#    unlink $headerfile unless $debug;
#    unlink $annfile unless $debug;
#    unlink $preprocessedfile unless $debug;
#    die "Preprocessing2 failed: $!";
#  }
#  unlink "$annfile" || warn "Could not delete $annfile: $!" unless $debug;
}

if ($support) {
  unlink $preprocessedfile if (defined $ARGV[0] && !$debug);
  die;
}

print STDERR "Reading gene list\n" if $debug;
my %genes = (); # Symbol => (Chrom => (chrom, start, stop, strand)); Hash of hashes of arrays
open(GENES, "< $genefile") || die "Could not open $genefile: $!";
foreach my $geneline (<GENES>) { # Parse gene file, recording the chromosome, strand, 5'-most start coordinate, and 3'-most stop coordinate found 
  my ($genechrom, $genestart, $genestop, $genesymbol, $genestrand) = split(/\s+/,$geneline);
  if (defined $genes{$genesymbol}->{$genechrom}) { ## Assume strand stays constant
    $genes{$genesymbol}->{$genechrom}->[0] = min($genes{$genesymbol}->{$genechrom}->[0], $genestart);
    $genes{$genesymbol}->{$genechrom}->[1] = max($genes{$genesymbol}->{$genechrom}->[1], $genestop);
  } else {
    $genes{$genesymbol}->{$genechrom} = [$genestart, $genestop, $genestrand];
  }
}
close GENES;

#my $lastbndmateid = "";
print STDERR "Entering loop\n" if $debug;
while (my $inputline = <IN>) {
  if ($inputline =~ /^#/) {
#    print OUT $inputline;
    next;
  }
  # Parse line
  my @splitline = split(/\s+/,$inputline);
  my ($leftchrom, $leftstart, $leftstop, $rightchrom, $rightstart, $rightstop, $id, $svtype, $info_a, $info_b) = @splitline[0..6, 10, 12, 13];
#  my ($mateid,$rightchrom,$rightpos,$rightstart,$rightstop,$leftstart,$leftstop,@cipos,@ciend,$mateoutputline,@splitmateline,$mateinfo,$mateline);

  $svtype = substr($svtype, 0, 3);
  unless (exists $types{$svtype}) {
    warn "Unrecognized SVTYPE $svtype at line ",$.," of preprocessed VCF file\n";
    print $inputline;
    next;
  }
  my ($probleft,$probright) = getfields($info_a,"PRPOS","PREND") if $weight;
  my ($leftgenenames, $rightgenenames);
  if ($info_b eq ".") { 
    ($leftgenenames,$rightgenenames) = getfields($info_a,"left_Gene","right_Gene");
  } else {
    $leftgenenames = getfields($info_a,"Gene");
    $rightgenenames = getfields($info_b,"Gene");
  }
  my @leftgenenames = split(/\|/,$leftgenenames);
  my @rightgenenames = split(/\|/,$rightgenenames);
  my $localweight = $weight && $probleft;
  
  my ($leftintrons,$rightintrons) = getfields($info_a,"left_Intron","right_Intron") if $svtype eq 'INV';


#  if ($svtype eq 'BND') {
#    my ($mateid) = ($id =~ /^(\d+)_(?:1|2)/);
#    $mateid = $id unless $mateid;
#    my ($nextmateid) = ($i == $#inputlines ? ("") : ((split(/\s+/,$inputlines[$i+1]))[2] =~ /^(\d+)_(?:1|2)/)); # $nextmateid is set to the mateid of the next line, unless the current line is the final line in the file. In this case, $nextmateid is set to the empty string because we know it's not the first mate of a consecutive pair
#    my $firstmate = ($nextmateid && $mateid == $nextmateid); # First mate of a consecutive pair
#    my $secondmate = ($lastbndmateid && $mateid == $lastbndmateid); # Second mate of a consecutive pair
#    $singletonbnd = !($firstmate || $secondmate);
#    if ($secondmate || $singletonbnd) { # Is this the second mate of this BND seen or a singleton BND? If the second of a pair, get annotations for first mate. If a singleton, get rid of one set of annotations
#      $lastbndmateid = "";
#      if ($singletonbnd) { # Make sure CIEND is present - otherwise, skip line and print error
#	warn "Could not process variant $id because mate is absent and no CIEND field was found" && next unless $ciend;
#      } else { # Get mateinfo
#	my $mateline = $inputlines[$i-1];
#	@splitmateline = split(/\s+/,$mateline);
#	$mateinfo = $splitmateline[7];
#      }
#      if ($info =~ /SECONDARY/) { # Current line is secondary (right), so mate is primary (left)
#	($leftexongenenames, $leftgenenames, $leftintrons) = ($singletonbnd ? ("","","") : getfields($mateinfo,"ExonGeneNames","Gene","Intron")); # Get mate annotations
#	($ciend, $probright) = getfields($mateinfo, "CIPOS", "PRPOS") unless $singletonbnd; # Overwrite (possibly empty) ciend with mate's CIPOS and get mate's PRPOS if exists
#	($rightchrom, $rightpos) = ($leftchrom, $leftpos);
#	($leftchrom,$leftpos) = ($alt =~ /([\w.]+):(\d+)/);
#	($cipos, $ciend) = ($ciend, $cipos); # Switch $cipos and $ciend
#	($probleft, $probright) = ($probright, $probleft);
#      } else { # Current line is primary (left), so mate is secondary (right)
#	($rightexongenenames, $rightgenenames, $rightintrons) = ($singletonbnd ? ("","","") : getfields($mateinfo,"ExonGeneNames","Gene","Intron")); # Get mate annotations
#	($ciend, $probright) = getfields($mateinfo, "CIPOS", "PRPOS") unless $singletonbnd; # Overwrite (possibly empty) ciend with mate's CIPOS and get mate's PRPOS if exists
#	($rightchrom,$rightpos) = ($alt =~ /([\w.]+):(\d+)/); # Get right breakend coordinates if standard BND
#	unless ($rightchrom && $rightpos) { # Get right breakend coordinates if reclassified BND
#	  $rightchrom = $leftchrom;
#	  $rightpos = getfields($info, "END");
#	}
#      }
#      undef $mateline unless $singletonbnd;
#    } else { # Must be the first mate of a consecutive pair
#      $lastbndmateid = $mateid;
#      ($leftchrom, $leftpos, $id, $alt, $info,$mateid,$rightchrom,$rightpos,$rightstart,$rightstop,$leftstart,$leftstop,$svtype,$spanexongenenames,$spangenenames,$leftexongenenames,$leftgenenames,$rightexongenenames,$rightgenenames,@leftgenenames,@rightgenenames,$leftintrons,$rightintrons) = (); ## Get rid of old variables
#      next;
#    }
#  } else { # Reset $lastbndmateid if variant isn't a BND
#    $lastbndmateid = "";
#    $singletonbnd = 0;
#  }

  my @probleft = split(/,/,$probleft) if $probleft;
  my @probright = split(/,/,$probright) if $probright;

  # Calculate start/stop coordinates from POS and CI
#  unless ($svtype eq "BND") {
#    $rightchrom = $leftchrom;
#    $rightpos = getfields($info,"END");
#  }
#  @cipos = split(/,/,$cipos);
#  $leftstart = $leftpos + $cipos[0];
#  $leftstop = $leftpos + $cipos[1];
#  @ciend = split(/,/,$ciend);
#  $rightstart = $rightpos + $ciend[0];
#  $rightstop = $rightpos + $ciend[1];

  my %scores = (); # Interval => List of scores by op; e.g. (LEFT => (MAXLEFT, SUMLEFT, TOP100LEFT, MEANLEFT), RIGHT => (MAXRIGHT, SUMRIGHT, TOP100RIGHT, MEANRIGHT))


  $scores{"LEFT"} = cscoreop($caddfile, $localweight, $ops, $leftchrom, $leftstart, $leftstop, \@probleft, $topn);
  $scores{"RIGHT"} = cscoreop($caddfile, $localweight, $ops, $rightchrom, $rightstart, $rightstop, \@probright, $topn);

  if ($svtype eq "DEL" || $svtype eq "DUP" || $svtype eq "CNV" || $svtype eq "MEI") {
    if ($rightstop - $leftstart > 1000000) {
      $scores{"SPAN"} = ($ops eq "ALL" ? [100, 100, 100, 100] : [100]);
    } else {
      $scores{"SPAN"} = cscoreop($caddfile, "", $ops, $leftchrom, $leftstart, $rightstop, "", $topn);
    }
  }

  if ($svtype eq "INV" || $svtype eq "TRX") {
    my %leftintrons = map {$_ => 1} (split(/\|/,$leftintrons));
    my @rightintrons = split(/\|/,$rightintrons);
    ## At worst, $leftintrons and $rightintrons are lists of introns. The only case in which the gene is not disrupted is if both lists are equal and nonempty, meaning that in every gene hit by this variant, both ends of the variant are confined to the same intron
    my $sameintrons = scalar (grep {$leftintrons{$_}} @rightintrons) == scalar @rightintrons && scalar @rightintrons > 0;
    unless ($sameintrons) {
      my $leftscore = truncationscore($leftchrom, $leftstart, $leftstop, $leftgenenames, \%genes, $caddfile, $ops, $topn);
      $scores{"LTRUNC"} = $leftscore if $leftscore;
      my $rightscore = truncationscore($rightchrom, $rightstart, $rightstop, $rightgenenames, \%genes, $caddfile, $ops, $topn);
      $scores{"RTRUNC"} = $rightscore if $rightscore;
    }
  }

  # Calculate maximum C score depending on SV type
#  if ($svtype eq "DEL" || $svtype eq "DUP" || $svtype eq "CNV" || $svtype eq "MEI") {
#    my $spanscore;
#    if ($rightstop - $leftstart > 1000000) {
#      $spanscore = ($ops eq "ALL" ? [100, 100, 100] : 100);
#    } else {
#      $spanscore = cscoreop($caddfile, "", $ops, $leftchrom, $leftstart, $rightstop, "");
#    }
#    $leftscore = cscoreop($caddfile, $localweight, $ops, $leftchrom, $leftstart, $leftstop, \@probleft);
#    $rightscore = cscoreop($caddfile, $localweight, $ops, $rightchrom, $rightstart, $rightstop, \@probright);
#    $info .= ($ops eq "ALL" ? ";SVSCOREMAX_SPAN=$spanscore->[0];SVSCOREMAX_LEFT=$leftscore->[0];SVSCOREMAX_RIGHT=$rightscore->[0];SVSCORESUM_SPAN=$spanscore->[1];SVSCORESUM_LEFT=$leftscore->[1];SVSCORESUM_RIGHT=$rightscore->[1];SVSCORETOP${topn}_SPAN=$spanscore->[2];SVSCORETOP${topn}_LEFT=$leftscore->[2];SVSCORETOP${topn}_RIGHT=$rightscore->[2]" : ";SVSCORE${ops}_SPAN=$spanscore;SVSCORE${ops}_LEFT=$leftscore;SVSCORE${ops}_RIGHT=$rightscore");
#    undef $spanscore;
#  } elsif ($svtype eq "INV" || $svtype eq "BND") {
#    $leftscore = cscoreop($caddfile, $localweight, $ops, $leftchrom, $leftstart, $leftstop, \@probleft);
#    $rightscore = cscoreop($caddfile, $localweight, $ops, $rightchrom, $rightstart, $rightstop, \@probright);
#    my ($sameintrons,@lefttruncationscores,@lefttruncationscoressum,@lefttruncationscorestop,@righttruncationscores,@righttruncationscoressum,@righttruncationscorestop,$lefttruncationscore,$righttruncationscore,%leftintrons,@rightintrons) = ();
#    unless ($singletonbnd) {
#      %leftintrons = map {$_ => 1} (split(/\|/,$leftintrons));
#      @rightintrons = split(/\|/,$rightintrons);
#      ## At worst, $leftintrons and $rightintrons are lists of introns. The only case in which the gene is not disrupted is if both lists are equal and nonempty, meaning that in every gene hit by this variant, both ends of the variant are confined to the same intron
#      $sameintrons = scalar (grep {$leftintrons{$_}} @rightintrons) == scalar @rightintrons && scalar @rightintrons > 0;
#    }
#    if (($leftgenenames || $rightgenenames) && ($singletonbnd || !$sameintrons)) { # Gene is being truncated - left or right breakend hits a gene, and the breakends are not confined to the same introns (or, if the variant is a singleton BND, this latter condition is not necessary)
#      foreach my $gene (split(/\|/,$leftgenenames)) {
#	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$leftchrom}}[0..2];	
#	my $cscoreopres;
#	if ($genestrand eq '+') {
#	  $cscoreopres = cscoreop($caddfile, "", $ops, $leftchrom, max($genestart,$leftstart),$genestop, ""); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
#	} else {
#	  $cscoreopres = cscoreop($caddfile, "", $ops, $leftchrom, $genestart,min($genestop,$leftstop), ""); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
#	}
#	if ($ops eq 'ALL') {
#	  push @lefttruncationscores,$cscoreopres->[0] unless $cscoreopres->[0] == -1;
#	  push @lefttruncationscoressum,$cscoreopres->[1] unless $cscoreopres->[1] == -1;
#	  push @lefttruncationscorestop,$cscoreopres->[1] unless $cscoreopres->[2] == -1;
#	} else {
#	  push @lefttruncationscores,$cscoreopres unless $cscoreopres == -1;
#	}
#      }
#      $lefttruncationscore = ($ops eq 'ALL' ? [max(@lefttruncationscores), max(@lefttruncationscoressum), max(@lefttruncationscorestop)] : max(@lefttruncationscores)) if @lefttruncationscores;
#      foreach my $gene (split(/\|/,$rightgenenames)) {
#	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$rightchrom}}[0..2];	
#	my $cscoreopres;
#	if ($genestrand eq '+') {
#	  $cscoreopres = cscoreop($caddfile, "", $ops, $rightchrom, max($genestart,$rightstart),$genestop, ""); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
#	} else {
#	  $cscoreopres = cscoreop($caddfile, "", $ops, $rightchrom, $genestart,min($genestop,$rightstop), ""); # Start from beginning of gene, stop at end of gene or breakend, whichever is further right (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
#	}
#	if ($ops eq 'ALL') {
#	  push @righttruncationscores,$cscoreopres->[0] unless $cscoreopres->[0] == -1;
#	  push @righttruncationscoressum,$cscoreopres->[1] unless $cscoreopres->[1] == -1;
#	  push @righttruncationscorestop,$cscoreopres->[1] unless $cscoreopres->[2] == -1;
#	} else {
#	  push @righttruncationscores,$cscoreopres unless $cscoreopres == -1;
#	}
#      }
#      $righttruncationscore = ($ops eq 'ALL' ? [max(@righttruncationscores), max(@righttruncationscoressum), max(@righttruncationscorestop)] : max(@righttruncationscores)) if @righttruncationscores;
#    }
# }
#    ($sameintrons, %leftintrons,@rightintrons) = (); # Get rid of old variables

# Transpose %scores hash to make taking maxes easy
#  foreach my $x (@{$scores{$interval}}) { ## DEBUG
#    print "$x ";
#  }
  my %scoresbyop = ();
  foreach my $interval (keys %scores) { # LEFT, RIGHT, (SPAN, LTRUNC, RTRUNC)
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
	$info_b .= ";SVSCORE${op}_$interval=$scores{$interval}->[$operations{$op}]" unless $info_b eq "." || $info_b eq "MISSING";
      }
    }
  }


#  my $addtoinfo;
#    if ($ops eq "ALL") {
#      $addtoinfo = ";SVSCOREMAX_LEFT=$leftscore->[0];SVSCOREMAX_RIGHT=$rightscore->[0];SVSCORESUM_LEFT=$leftscore->[1];SVSCORESUM_RIGHT=$rightscore->[1];SVSCORETOP${topn}_LEFT=$leftscore->[2];SVSCORETOP${topn}_RIGHT=$rightscore->[2]" . (defined $lefttruncationscore ? ";SVSCOREMAX_LTRUNC=$lefttruncationscore->[0];SVSCORESUM_LTRUNC=$lefttruncationscore->[1];SVSCORETOP${topn}_LTRUNC=$lefttruncationscore->[2]" : "") . (defined $righttruncationscore ? ";SVSCOREMAX_RTRUNC=$righttruncationscore->[0];SVSCORESUM_RTRUNC=$righttruncationscore->[1];SVSCORETOP${topn}_RTRUNC=$righttruncationscore->[2]" : "");
#    } else {
#      $addtoinfo = ";SVSCORE${ops}_LEFT=$leftscore;SVSCORE${ops}_RIGHT=$rightscore" . (defined $lefttruncationscore ? ";SVSCORE${ops}_LTRUNC=$lefttruncationscore" : "") . (defined $righttruncationscore ? ";SVSCORE${ops}_RTRUNC=$righttruncationscore" : "");
#    }
#    $info .= $addtoinfo;
#    $mateinfo .= $addtoinfo if $svtype eq "BND" && !$singletonbnd;
#    (@lefttruncationscores,@righttruncationscores,@lefttruncationscoressum,@righttruncationscoressum,$lefttruncationscore,$righttruncationscore,$addtoinfo) = (); # Get rid of old variables
  

  # Multiplier for deletions and duplications which hit an exon, lower multiplier if one of these hits a gene but not an exon. Purposely not done for BND and INV. THESE WILL NEED TO BE PLACED IN ABOVE IF STATEMENT IN THE FUTURE
  #if ($spanexongenenames && ($svtype eq "DEL" || $svtype eq "DUP")) { 
  #  $spanscore *= 1.5;
  #} elsif($spangenenames && ($svtype eq "DEL" || $svtype eq "DUP")) {
  #  $spanscore *= 1.2;
  #}

  # For all types except INS, multiply left and right scores respectively if exon/gene is hit
  #if ($leftexongenenames && $svtype ne "INS") {
  #  $leftscore *= 1.5;
  #} elsif ($leftgenenames && $svtype ne "INS") {
  #  $leftscore *= 1.2;
  #}
  #if ($rightexongenenames && $svtype ne "INS") {
  #  $rightscore *= 1.5;
  #} elsif ($rightgenenames && $svtype ne "INS") {
  #  $rightscore *= 1.2;
  #}

#  my $outputline = "";
#  $mateoutputline = "" if $svtype eq "BND" && !$singletonbnd;
#  foreach my $i (0..$#splitline) { # Build output line
#    $outputline .= (($i == 7 ? $info : $splitline[$i]) . ($i < $#splitline ? "\t" : ""));
#    $mateoutputline .= (($i == 7 ? $mateinfo : $splitmateline[$i]) . ($i < $#splitline ? "\t" : "")) if $svtype eq "BND" && !$singletonbnd;
#  }
  $splitline[12] = $info_a;
  $splitline[13] = $info_b;
  print OUT join("\t",@splitline) . "\n";
#  push @outputlines, $mateoutputline if $svtype eq "BND" && !$singletonbnd;

#  ($leftchrom, $leftpos, $id, $alt, $info,$mateid,$rightchrom,$rightpos,$rightstart,$rightstop,$leftstart,$leftstop,$svtype,$spanexongenenames,$spangenenames,$leftexongenenames,$leftgenenames,$rightexongenenames,$rightgenenames,@leftgenenames,@rightgenenames,$leftintrons,$rightintrons,$rightchrom,$rightpos,$leftscore,$rightscore,$cipos,$ciend,@cipos,@ciend,$leftscore,$rightscore,$outputline,$mateoutputline,@splitmateline,$mateline,$mateinfo) = (); # Get rid of old variables

  print STDERR $.,", " if $debug;
}

# Extract chromosomes and starting positions for sorting 
#my @chroms;
#foreach my $line (@outputlines) {
#  push @chroms, (split(/\s+/,$line))[0];
#}
#my @starts;
#foreach my $line (@outputlines) {
#  push @starts, (split(/\s+/,$line))[1];
#}
#
# Sort and print
#foreach my $i (sort {$chroms[$a] cmp $chroms[$b] || $starts[$a] <=> $starts[$b]} (0..$#outputlines)) {
#  print "$outputlines[$i]\n";
#}

close OUT;
system("bedpetovcf -b $bedpeout > $vcfout");
system("grep \"^#\" $vcfout > $vcfout.header");
print `grep -v "^#" $vcfout | sort -k1,1V -k2,2n | cat $vcfout.header -`;
unlink "$vcfout.header";

unless ($debug) {
  unlink "$preprocessedfile";
  unlink "$vcfout";
  unlink "$bedpeout";
  rmdir "svscoretmp" || warn "Could not delete svscoretmp: $!";
}

sub cscoreop { # Apply operation(s) specified in $ops to C scores within a given region using CADD data
  my ($filename, $weight, $ops, $chrom, $start, $stop, $probdist, $topn) = @_;
  my @probdist = @{$probdist} if $weight;
  my (@scores,$res) = ();
#  if ($stop-$start>1000000) {## DEBUG
#    return ($ops eq 'ALL' ? [-1, -1, -1 -1] : [-1]) if ($stop-$start>1000000); # Short circuit if region is too big - this is usually caused by faulty annotations
#    print STDERR "Interval too big: $chrom:$start-$stop\n";
#  }
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

  @scores = pairwise {$a * $b}	@scores, @probdist if $weight;
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
  ($filename,$chrom,$start,$stop,$ops,$probdist,@probdist,@scores,$tabixoutput,@tabixoutputlines,$topn) = (); # Get rid of old variables
  return $res;
}

sub getfields { # Parse info field of VCF line, getting fields specified in @_. $_[0] must be the info field itself. Returns list of field values if more than one field is being requested; otherwise, returns a scalar value representing the requested field
  my $info = shift @_;
  my @ans;
  foreach my $i (0..$#_) {
    my ($ann) = ($info =~ /(?:;|^)$_[$i]=(.*?)(?:;|$)/);
    push @ans, ($ann ? $ann : "");
  }
  $info = undef; # Get rid of old variable
  if (@ans > 1) {
    return @ans;
  } else {
    return $ans[0];
  }
}

sub truncationscore { # Calculate truncation score based on the coordinates of a given breakend and the names of the genes it hits
  my ($chrom, $start, $stop, $genenames, $genesref, $caddfile, $ops, $topn) = @_;
  return "" unless $genenames;
  my %genes = %{$genesref};
  my (@truncationscores, @truncationscoressum, @truncationscorestop);
  foreach my $gene (split(/\|/,$genenames)) {
    my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$chrom}}[0..2];	
    my $cscoreopres;
    if ($genestrand eq '+') {
      $cscoreopres = cscoreop($caddfile, "", $ops, $chrom, max($genestart,$start),$genestop, "", $topn); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
    } else {
      $cscoreopres = cscoreop($caddfile, "", $ops, $chrom, $genestart,min($genestop,$stop), "", $topn); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but none of the supported operations are order-dependent)
    }
    if ($ops eq 'ALL') {
      push @truncationscores,$cscoreopres->[0] unless $cscoreopres->[0] == -1;
      push @truncationscoressum,$cscoreopres->[1] unless $cscoreopres->[1] == -1;
      push @truncationscorestop,$cscoreopres->[1] unless $cscoreopres->[2] == -1;
    } else {
      push @truncationscores,$cscoreopres unless $cscoreopres == -1;
    }
  }
  if (@truncationscores) {
    return ($ops eq 'ALL' ? [max(@truncationscores), max(@truncationscoressum), max(@truncationscorestop)] : [max(@truncationscores)]);
  } else {
    return "";
  }
}

sub main::HELP_MESSAGE() {
  print "usage: ./svscore.pl [-dsvw] [-o op] [-t topnumber] [-g genefile] [-m geneannotationcolumn] [-e exonfile] [-c caddfile] vcf
    -d	      Debug mode, keeps intermediate and supporting files, displays progress
    -s	      Create/download supporting files and quit
    -v	      Verbose mode - show all calculated scores (left/right/span/ltrunc/rtrunc, as appropriate)
    -w	      Weight CADD scores in breakends by probability distribution (requires PRPOS/PREND in INFO field)
    -o	      Specify operation to perform on CADD scores (must be sum, max, top, or all - defaults to all)
    -t	      Number of bases for TOP to take mean of (under -o top or -o all) (100)
    -g	      Points to gene BED file (refGene.genes.b37.bed)
    -m	      Column number for annotation in gene BED file to be added to VCF (4)
    -e	      Points to exon BED file (refGene.exons.b37.bed)
    -c	      Points to whole_genome_SNVs.tsv.gz (defaults to current directory)

    --help    Display this message
    --version Display version\n"
}

#    -n	      Column number for annotation in exon BED file to be added to VCF (4 if using -e, 5 otherwise)

sub main::VERSION_MESSAGE() {
  print "SVScore version $main::VERSION\n";
}
