#!/usr/bin/perl -w 
## Author: Liron Ganel
## Laboratory of Ira Hall, McDonnell Genome Institute
## Washington University in St. Louis
## Version 0.3

use strict;
use Getopt::Std;
use List::Util qw(max min sum);
use List::MoreUtils qw(pairwise);
use Time::HiRes qw(gettimeofday);

$Getopt::Std::STANDARD_HELP_VERSION = 1; # Make --help and --version flags halt execution
$main::VERSION = '0.3';

my %options = ();
getopts('dswc:g:e:n:o:',\%options);

my $debug = defined $options{'d'};
my $support = defined $options{'s'};
my $weight = defined $options{'w'};

&main::HELP_MESSAGE() && die unless defined $ARGV[0] || $support;

my $caddfile = (defined $options{'c'} ? $options{'c'} : 'whole_genome_SNVs.tsv.gz');
my $genefile = (defined $options{'g'} ? $options{'g'} : 'refGene.genes.b37.bed');
my $exonfile = (defined $options{'e'} ? $options{'e'} : 'refGene.exons.b37.bed');
my $exonanncolumn = (defined $options{'n'} && defined $exonfile ? $options{'n'} : 5);
warn "Exon annotation column provided without nonstandard exon annotation file - defaulting to standard exon annotation file" if defined $options{'n'} && !defined $options{'e'};
die "Nonstandard exon annotation file without column number - rerun with -n option" if !defined $options{'n'} && defined $options{'e'};
$options{'o'} =~ tr/[a-z]/[A-Z]/ if defined $options{'o'};
my $ops = (defined $options{'o'} ? $options{'o'} : 'BOTH');
die "Unrecognized operation specified: $ops" unless ($ops eq 'SUM' || $ops eq 'MAX' || $ops eq 'BOTH');
my $compressed = ($ARGV[0] =~ /\.gz$/) if defined $ARGV[0];
my ($annfile, $headerfile, $preprocessedfile);

##TODO PRIORITY 2: Enable piping input through STDIN - use an option to specify input file rather than @ARGV

# Set up all necessary preprocessing to be taken care of before analysis can begin. This includes decompression, annotation using vcfanno, and generation of intron/exon/gene files, whichever are necessary. May be a little slower than necessary in certain situations because some arguments are supplied by piping cat output rather than supplying filenames directly.
if ($exonfile eq 'refGene.exons.b37.bed' && !-s $exonfile) { # Generate exon file if necessary
  print STDERR "Generating exon file\n" if $debug;
  system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); n=int(\$9); split(\$10,start,\",\");split(\$11,end,\",\"); for(i=1;i<=n;++i) {print \$3,start[i],end[i],\$2\".\"i,\$13,\$2; } }' OFS=\"\t\" | sort -k 1,1 -k 2,2n | uniq > refGene.exons.b37.bed");
} elsif ($exonfile ne 'refGene.exons.b37.bed' && !-s $exonfile) {
  die "$exonfile not found or empty!";
}

if ($genefile eq 'refGene.genes.b37.bed' && !-s $genefile) { # Generate gene file if necessary
  print STDERR "Generating gene file\n" if $debug;
  system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); print \$3,\$5,\$6,\$4,\$13}' OFS=\"\\t\" | sort -k 1,1 -k 2,2n | uniq > refGene.genes.b37.bed");
} elsif ($genefile ne 'refGene.genes.b37.bed' && !-s $genefile) {
  die "$genefile not found or empty!";
}

unless (-s 'introns.bed') { # Generate intron file if necessary - add column with unique intron ID equal to line number (assuming introns.bed has no header line) and sort
  print STDERR "Generating intron file\n" if $debug;
  system("bedtools subtract -a $genefile -b $exonfile | sort -u -k 1,1 -k 2,2n | awk '{print \$0 \"\t\" NR}' > introns.bed");
}

my $intronnumcolumn = `head -n 1 introns.bed | awk '{print NF}'`;
chomp $intronnumcolumn;

# Write conf.toml file
print STDERR "Writing config file\n" if $debug;
open(CONFIG, "> conf.toml") || die "Could not open conf.toml: $!";
print CONFIG "[[annotation]]\nfile=\"$genefile\"\nnames=[\"Gene\"]\ncolumns=[5]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"$exonfile\"\nnames=[\"ExonGeneNames\"]\ncolumns=[$exonanncolumn]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"introns.bed\"\nnames=[\"Intron\"]\ncolumns=[$intronnumcolumn]\nops=[\"uniq\"]";
close CONFIG;

# Create first preprocessing command - annotation is done without normalization because REF and ALT nucleotides are not included in VCFs describing SVs
if (defined $ARGV[0]) {
  print STDERR "Preparing preprocessing command\n" if $debug;
  my ($prefix) = ($ARGV[0] =~ /^(?:.*\/)?(.*)\.vcf(?:\.gz)?$/);

  # Tag intermediate files with timestamp to avoid collisions
  my $time = gettimeofday();
  $preprocessedfile = "$prefix.preprocess$time.vcf";
  $annfile = "$prefix.ann$time.vcf";
  $headerfile = "${prefix}header$time";
  my $preprocess = ($compressed ? "z": "") . "cat $ARGV[0] | awk '\$0~\"^#\" {print \$0; next } { print \$0 | \"sort -k1,1 -k2,2n\" }' | vcfanno -ends conf.toml - > $annfile; grep '^#' $annfile > $headerfile"; # Sort, annotate, grab header
  print STDERR "Preprocessing command 1: $preprocess\n" if $debug;
  die "Preprocessing failed: $!" if system($preprocess);

  # Update header
  open(HEADER, "$headerfile") || die "Could not open $headerfile: $!";
  my @newheader = ();
  my @oldheader = <HEADER>;
  close HEADER;
  my $headerline;
  while(($headerline = (shift @oldheader)) !~ /^##INFO/) {
    push @newheader, $headerline;
  }
  unshift @oldheader, $headerline; # Return first info line to top of stack
  while(($headerline = (shift @oldheader)) =~ /^##INFO/) {
    push @newheader, $headerline;
  }
  unshift @oldheader, $headerline; # Return first format line to top of stack

  if ($ops eq 'MAX' || $ops eq 'BOTH') {
    push @newheader, "##INFO=<ID=SVSCOREMAX_LEFT,Number=1,Type=Float,Description=\"Maximum C score in left breakend of structural variant, weighted by probability distribution\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_RIGHT,Number=1,Type=Float,Description=\"Maximum C score in right breakend of structural variant, weighted by probability distribution\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_SPAN,Number=1,Type=Float,Description=\"Maximum C score in outer span of structural variant, weighted by probability distribution\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_LTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of left breakend to end of gene, weighted by probability distribution\">\n";
    push @newheader, "##INFO=<ID=SVSCOREMAX_RTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of right breakend to end of gene, weighted by probability distribution\">\n";
  }
  if ($ops eq 'SUM' || $ops eq 'BOTH') {
    push @newheader, "##INFO=<ID=SVSCORESUM_LEFT,Number=1,Type=Float,Description=\"Sum of C scores in left breakend of structural variant, weighted by probability distribution\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_RIGHT,Number=1,Type=Float,Description=\"Sum of C scores in right breakend of structural variant, weighted by probability distribution\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_SPAN,Number=1,Type=Float,Description=\"Sum of C scores in outer span of structural variant, weighted by probability distribution\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_LTRUNC,Number=1,Type=Float,Description=\"Sum of C scores from beginning of left breakend to end of gene, weighted by probability distribution\">\n";
    push @newheader, "##INFO=<ID=SVSCORESUM_RTRUNC,Number=1,Type=Float,Description=\"Sum of C scores from beginning of right breakend to end of gene, weighted by probability distribution\">\n";
  }
  push @newheader, @oldheader;
  open(HEADER, "> $headerfile") || die "Could not open $headerfile: $!";
  foreach (@newheader) {
    print HEADER;
  }
  if ($weight && !(grep {/^##INFO=<ID=PRPOS,/} @newheader)) {
    $weight = 0;
    warn "PRPOS not found in header - SVScore will not weight CADD scores\n";
  }

  # Create and execute second preprocessing command
  my $preprocess2 = "grep -v '^#' $annfile | sort -k 3,3 | cat $headerfile - > $preprocessedfile; rm -f $headerfile"; # Sort by ID, add new header, clean up
  print STDERR "Preprocessing command 2: $preprocess2\n" if $debug;
  die "Preprocessing2 failed: $!" if system("$preprocess2");
  unlink "$annfile" || warn "Could not delete $annfile: $!" unless $debug;
}

if ($support) {
  unlink $preprocessedfile if defined $ARGV[0];
  die;
}

print STDERR "Reading gene list\n" if $debug;
my %genes = (); # Symbol => (Chrom => (chrom, start, stop, strand)); Hash of hashes of arrays
open(GENES, "< $genefile") || die "Could not open $genefile: $!";
foreach my $geneline (<GENES>) { # Parse gene file, recording the chromosome, strand, 5'-most start coordinate, and 3'-most stop coordinate found 
  my ($genechrom, $genestart, $genestop, $genestrand, $genesymbol) = split(/\s+/,$geneline);
  if (defined $genes{$genesymbol}->{$genechrom}) { ## Assume strand stays constant
    $genes{$genesymbol}->{$genechrom}->[0] = min($genes{$genesymbol}->{$genechrom}->[0], $genestart);
    $genes{$genesymbol}->{$genechrom}->[1] = max($genes{$genesymbol}->{$genechrom}->[1], $genestop);
  } else {
    $genes{$genesymbol}->{$genechrom} = [$genestart, $genestop, $genestrand];
  }
}

open(IN, "< $preprocessedfile") || die "Could not open $preprocessedfile: $!";

my @outputlines;
my $lastbndmateid = "";
print STDERR "Reading input file\n" if $debug;
my @inputlines = <IN>;

print STDERR "Entering loop\n" if $debug;
foreach my $i (0..$#inputlines) {
  my $vcfline = $inputlines[$i];
  if ($vcfline =~ /^#/) {
    print $vcfline;
    next;
  }
  # Parse line
  my @splitline = split(/\s+/,$vcfline);
  my ($leftchrom, $leftpos, $id, $alt, $info) = @splitline[0..2, 4, 7];
  my ($mateid,$rightchrom,$rightpos,$rightstart,$rightstop,$leftstart,$leftstop,@cipos,@ciend,$mateoutputline,@splitmateline,$mateinfo,$mateline,$singletonbnd);

  my ($svtype) = ($info =~ /SVTYPE=(\w{3})/);
  my ($spanexongenenames,$spangenenames,$leftexongenenames,$leftgenenames,$rightexongenenames,$rightgenenames,$cipos,$ciend,$probleft,$probright) = getfields($info,"ExonGeneNames","Gene","left_ExonGeneNames","left_Gene","right_ExonGeneNames","right_Gene","CIPOS","CIEND","PRPOS","PREND");
  my @leftgenenames = split(/\|/,$leftgenenames);
  my @rightgenenames = split(/\|/,$rightgenenames);
  my $localweight = $weight && $probleft;
  
  my ($leftintrons,$rightintrons) = getfields($info,"left_Intron","right_Intron") if $svtype eq 'BND' || $svtype eq 'INV';


  if ($svtype eq 'BND') {
    my ($mateid) = ($id =~ /^(\d+)_(?:1|2)/);
    my ($nextmateid) = ($i == $#inputlines ? ("") : ((split(/\s+/,$inputlines[$i+1]))[2] =~ /^(\d+)_(?:1|2)/)); # $nextmateid is set to the mateid of the next line, unless the current line is the final line in the file. In this case, $nextmateid is set to the empty string because we know it's not the first mate of a consecutive pair
    my $firstmate = ($nextmateid && $mateid == $nextmateid); # First mate of a consecutive pair
    my $secondmate = ($lastbndmateid && $mateid == $lastbndmateid); # Second mate of a consecutive pair
    $singletonbnd = !($firstmate || $secondmate);
    if ($secondmate || $singletonbnd) { # Is this the second mate of this BND seen or a singleton BND? If the second of a pair, get annotations for first mate. If a singleton, get rid of one set of annotations
      $lastbndmateid = "";
      if ($singletonbnd) { # Make sure CIEND is present - otherwise, skip line and print error
	warn "Could not process variant $id because mate is absent and no CIEND field was found" && next unless $ciend;
      } else { # Get mateinfo
	my $mateline = $inputlines[$i-1];
	@splitmateline = split(/\s+/,$mateline);
	$mateinfo = $splitmateline[7];
      }
      if ($info =~ /SECONDARY/) { # Current line is secondary (right), so mate is primary (left)
	($leftexongenenames, $leftgenenames, $leftintrons) = ($singletonbnd ? ("","","") : getfields($mateinfo,"ExonGeneNames","Gene","Intron")); # Get mate annotations
	($ciend, $probright) = getfields($mateinfo, "CIPOS", "PRPOS") unless $singletonbnd; # Overwrite (possibly empty) ciend with mate's CIPOS and get mate's PRPOS if exists
	($rightchrom, $rightpos) = ($leftchrom, $leftpos);
	($leftchrom,$leftpos) = ($alt =~ /([\w.]+):(\d+)/);
	($cipos, $ciend) = ($ciend, $cipos); # Switch $cipos and $ciend
	($probleft, $probright) = ($probright, $probleft);
      } else { # Current line is primary (left), so mate is secondary (right)
	($rightexongenenames, $rightgenenames, $rightintrons) = ($singletonbnd ? ("","","") : getfields($mateinfo,"ExonGeneNames","Gene","Intron")); # Get mate annotations
	($ciend, $probright) = getfields($mateinfo, "CIPOS", "PRPOS") unless $singletonbnd; # Overwrite (possibly empty) ciend with mate's CIPOS and get mate's PRPOS if exists
	($rightchrom,$rightpos) = ($alt =~ /([\w.]+):(\d+)/);
      }
      undef $mateline unless $singletonbnd;
    } else { # Must be the first mate of a consecutive pair
      $lastbndmateid = $mateid;
      ($leftchrom, $leftpos, $id, $alt, $info,$mateid,$rightchrom,$rightpos,$rightstart,$rightstop,$leftstart,$leftstop,$svtype,$spanexongenenames,$spangenenames,$leftexongenenames,$leftgenenames,$rightexongenenames,$rightgenenames,@leftgenenames,@rightgenenames,$leftintrons,$rightintrons) = (); ## Get rid of old variables
      next;
    }
  } else { # Reset $lastbndmateid if variant isn't a BND
    $lastbndmateid = "";
    $singletonbnd = 0;
  }

#  print STDERR "$id\n"; ## DEBUG
  my @probleft = split(/,/,$probleft);
  my @probright = split(/,/,$probright) if $probright;

  # Calculate start/stop coordinates from POS and CI
  unless ($svtype eq "BND") {
    $rightchrom = $leftchrom;
    $rightpos = getfields($info,"END");
  }
  @cipos = split(/,/,$cipos);
  $leftstart = $leftpos + $cipos[0];
  $leftstop = $leftpos + $cipos[1];
  @ciend = split(/,/,$ciend);
  $rightstart = $rightpos + $ciend[0];
  $rightstop = $rightpos + $ciend[1];

  my ($leftscore,$rightscore);

  # Calculate maximum C score depending on SV type
  if ($svtype eq "DEL" || $svtype eq "DUP") {
#    print STDERR "Calculating spanscore:\n" if $debug; ## DEBUG
    my $spanscore;
    if ($rightstop - $leftstart > 1000000) {
      $spanscore = ($ops eq "BOTH" ? [100, 100] : 100);
    } else {
      $spanscore = cscoreop($caddfile, 0, $ops, $leftchrom, $leftstart, $rightstop, "");
    }
#    print STDERR "Calculating left score:\n" if $debug; ## DEBUG
    $leftscore = cscoreop($caddfile, $localweight, $ops, $leftchrom, $leftstart, $leftstop, \@probleft);
#    print STDERR "Calculating right score:\n" if $debug; ## DEBUG
    $rightscore = cscoreop($caddfile, $localweight, $ops, $rightchrom, $rightstart, $rightstop, \@probright);
    $info .= ($ops eq "BOTH" ? ";SVSCOREMAX_SPAN=$spanscore->[0];SVSCOREMAX_LEFT=$leftscore->[0];SVSCOREMAX_RIGHT=$rightscore->[0];SVSCORESUM_SPAN=$spanscore->[1];SVSCORESUM_LEFT=$leftscore->[1];SVSCORESUM_RIGHT=$rightscore->[1]" : ";SVSCORE${ops}_SPAN=$spanscore;SVSCORE${ops}_LEFT=$leftscore;SVSCORE${ops}_RIGHT=$rightscore");
    undef $spanscore;
  } elsif ($svtype eq "INV" || $svtype eq "BND") {
    $leftscore = cscoreop($caddfile, $localweight, $ops, $leftchrom, $leftstart, $leftstop, \@probleft);
    $rightscore = cscoreop($caddfile, $localweight, $ops, $rightchrom, $rightstart, $rightstop, \@probright);
    my ($sameintrons,@lefttruncationscores,@lefttruncationscoressum,@righttruncationscores,@righttruncationscoressum,$lefttruncationscore,$righttruncationscore,%leftintrons,@rightintrons) = ();
    unless ($singletonbnd) {
      %leftintrons = map {$_ => 1} (split(/\|/,$leftintrons));
      @rightintrons = split(/\|/,$rightintrons);
      ## At worst, $leftintrons and $rightintrons are lists of introns. The only case in which the gene is not disrupted is if both lists are equal and nonempty, meaning that in every gene hit by this variant, both ends of the variant are confined to the same intron
      $sameintrons = scalar (grep {$leftintrons{$_}} @rightintrons) == scalar @rightintrons && scalar @rightintrons > 0;
    }
    if (($leftgenenames || $rightgenenames) && ($singletonbnd || !$sameintrons)) { # Gene is being truncated - left or right breakend hits a gene, and the breakends are not confined to the same introns (or, if the variant is a singleton BND, this latter condition is not necessary)
      foreach my $gene (split(/\|/,$leftgenenames)) {
	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$leftchrom}}[0..2];	
	my $cscoreopres;
	if ($genestrand eq '+') {
	  $cscoreopres = cscoreop($caddfile, 0, $ops, $leftchrom, max($genestart,$leftstart),$genestop, ""); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
	} else {
	  $cscoreopres = cscoreop($caddfile, 0, $ops, $leftchrom, $genestart,min($genestop,$leftstop), ""); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
	}
	if ($ops eq 'BOTH') {
	  push @lefttruncationscores,$cscoreopres->[0] unless $cscoreopres->[0] == -1;
	  push @lefttruncationscoressum,$cscoreopres->[1] unless $cscoreopres->[1] == -1;
	} else {
	  push @lefttruncationscores,$cscoreopres unless $cscoreopres == -1;
	}
      }
      $lefttruncationscore = ($ops eq 'BOTH' ? [max(@lefttruncationscores), max(@lefttruncationscoressum)] : max(@lefttruncationscores)) if @lefttruncationscores;
      foreach my $gene (split(/\|/,$rightgenenames)) {
	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$rightchrom}}[0..2];	
	my $cscoreopres;
	if ($genestrand eq '+') {
	  $cscoreopres = cscoreop($caddfile, 0, $ops, $rightchrom, max($genestart,$rightstart),$genestop, ""); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
	} else {
	  $cscoreopres = cscoreop($caddfile, 0, $ops, $rightchrom, $genestart,min($genestop,$rightstop), ""); # Start from beginning of gene, stop at end of gene or breakend, whichever is further right (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
	}
	if ($ops eq 'BOTH') {
	  push @righttruncationscores,$cscoreopres->[0] unless $cscoreopres->[0] == -1;
	  push @righttruncationscoressum,$cscoreopres->[1] unless $cscoreopres->[1] == -1;
	} else {
	  push @righttruncationscores,$cscoreopres unless $cscoreopres == -1;
	}
      }
      $righttruncationscore = ($ops eq 'BOTH' ? [max(@righttruncationscores), max(@righttruncationscoressum)] : max(@righttruncationscores)) if @righttruncationscores;
    }
    ($sameintrons, %leftintrons,@rightintrons) = (); # Get rid of old variables
    my $addtoinfo;
    if ($ops eq "BOTH") {
      $addtoinfo = ";SVSCOREMAX_LEFT=$leftscore->[0];SVSCOREMAX_RIGHT=$rightscore->[0];SVSCORESUM_LEFT=$leftscore->[1];SVSCORESUM_RIGHT=$rightscore->[1]" . (defined $lefttruncationscore ? ";SVSCOREMAX_LTRUNC=$lefttruncationscore->[0];SVSCORESUM_LTRUNC=$lefttruncationscore->[1]" : "") . (defined $righttruncationscore ? ";SVSCOREMAX_RTRUNC=$righttruncationscore->[0];SVSCORESUM_RTRUNC=$righttruncationscore->[1]" : "");
    } else {
      $addtoinfo = ";SVSCORE${ops}_LEFT=$leftscore;SVSCORE${ops}_RIGHT=$rightscore" . (defined $lefttruncationscore ? ";SVSCORE${ops}_LTRUNC=$lefttruncationscore" : "") . (defined $righttruncationscore ? ";SVSCORE${ops}_RTRUNC=$righttruncationscore" : "");
    }
    $info .= $addtoinfo;
    $mateinfo .= $addtoinfo if $svtype eq "BND" && !$singletonbnd;
    (@lefttruncationscores,@righttruncationscores,@lefttruncationscoressum,@righttruncationscoressum,$lefttruncationscore,$righttruncationscore,$addtoinfo) = (); # Get rid of old variables
  } elsif ($svtype eq "INS") { # leftscore is base before insertion, rightscore is base after insertion
    $leftscore = cscoreop($caddfile, 0, $ops, $leftchrom, $leftstart-1, $leftstart-1, "");
    $rightscore = cscoreop($caddfile, 0, $ops, $rightchrom, $rightstart+1, $rightstart+1, "");
    $info .= ($ops eq "BOTH" ? ";SVSCOREMAX_LEFT=$leftscore->[0];SVSCORESUM_LEFT=$leftscore->[1];SVSCOREMAX_RIGHT=$rightscore->[0];SVSCORESUM_RIGHT=$rightscore->[1]" : ";SVSCORE${ops}_LEFT=$leftscore;SVSCORE${ops}_RIGHT=$rightscore");
  } else {
    die "Unrecognized SVTYPE $svtype at line ",$i+1," of annotated VCF file\n";
  }

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

  my $outputline = "";
  $mateoutputline = "" if $svtype eq "BND" && !$singletonbnd;
  foreach my $i (0..$#splitline) { # Build output line
    $outputline .= (($i == 7 ? $info : $splitline[$i]) . ($i < $#splitline ? "\t" : ""));
    $mateoutputline .= (($i == 7 ? $mateinfo : $splitmateline[$i]) . ($i < $#splitline ? "\t" : "")) if $svtype eq "BND" && !$singletonbnd;
  }
  push @outputlines, $outputline;
  push @outputlines, $mateoutputline if $svtype eq "BND" && !$singletonbnd;

  ($leftchrom, $leftpos, $id, $alt, $info,$mateid,$rightchrom,$rightpos,$rightstart,$rightstop,$leftstart,$leftstop,$svtype,$spanexongenenames,$spangenenames,$leftexongenenames,$leftgenenames,$rightexongenenames,$rightgenenames,@leftgenenames,@rightgenenames,$leftintrons,$rightintrons,$rightchrom,$rightpos,$leftscore,$rightscore,$cipos,$ciend,@cipos,@ciend,$leftscore,$rightscore,$outputline,$mateoutputline,@splitmateline,$mateline,$mateinfo) = (); # Get rid of old variables

  print STDERR $i+1, " " if $debug;
}

# Extract chromosomes and IDs for sorting 
my @chroms;
foreach my $line (@outputlines) {
  push @chroms, (split(/\s+/,$line))[0];
}
my @starts;
foreach my $line (@outputlines) {
  push @starts, (split(/\s+/,$line))[1];
}

# Sort and print
foreach my $i (sort {$chroms[$a] cmp $chroms[$b] || $starts[$a] <=> $starts[$b]} (0..$#outputlines)) {
  print "$outputlines[$i]\n";
}

unlink "$preprocessedfile" unless $debug;

sub cscoreop { # Apply operation specified in $ops to C scores within a given region using CADD data
  my ($filename, $weight, $ops, $chrom, $start, $stop, $probdist) = @_;
  my @probdist = @{$probdist} if $weight;
  my (@scores,$res) = ();
  my $tabixoutput = `tabix $filename $chrom:$start-$stop`;
#  warn "$chrom:$start-$stop\n" if $weight; ## DEBUG
  my @tabixoutputlines = split(/\n/,$tabixoutput);
  return ($ops eq 'BOTH' ? [-1, -1] : -1) unless (@tabixoutputlines == $stop-$start+1); # Short circuit if variant hits region with no C scores
  foreach my $taboutline (@tabixoutputlines) {
    push @scores, max(split(/,/,(split(/\s+/,$taboutline))[4]));
  }
#  warn "Scores: @scores\n\nProbdist: @probdist\n" if $weight; ## DEBUG
#  die "Scores: ",scalar @scores, "\t", "Probdist: ", scalar @probdist, "\tShould be: ",$stop-$start+1 if ($weight && scalar @scores != scalar @probdist); ## DEBUG
  @scores = pairwise {$a * $b}	@scores, @probdist if $weight;
  if ($ops eq 'MAX') {
    $res = max(@scores)
  } elsif ($ops eq 'SUM') {
    $res = sum(@scores)
  } else {
    $res = [max(@scores), sum(@scores)]
  }
  ($filename,$chrom,$start,$stop,$ops,$probdist,@probdist,@scores,$tabixoutput,@tabixoutputlines) = (); # Get rid of old variables
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

sub main::HELP_MESSAGE() {
  print "usage: ./svscore.pl [-ds] [-g genefile] [-e exonfile] [-n exonannotationcolumn] [-c caddfile] vcf
      -d        Debug (verbose) mode, keeps intermediate and supporting files
      -s        Create/download supporting files and quit
      -c        Points to whole_genome_SNVs.tsv.gz (defaults to current directory)
      -g        Used to point to gene BED file (refGene.genes.b37.bed)
      -e        Used to point to exon BED file (refGene.exons.b37.bed)
      -n        Column number for annotation in exon BED file to be added to VCF (5)
      -w        Weight CADD scores in breakends by probability distribution (requires PRPOS/PREND in INFO field)
      -o        Specify operation to perform on CADD scores (must be sum, max, or both - defaults to both)

      --help    Display this message
      --version Display version\n";
}

sub main::VERSION_MESSAGE() {
  print "SVScore version $main::VERSION\n";
}
