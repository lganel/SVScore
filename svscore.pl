#!/usr/bin/perl -w 
## Author: Liron Ganel
## Laboratory of Ira Hall, McDonnell Genome Institute
## Washington University in St. Louis
## Version 0.2

use strict;
use Getopt::Std;
use List::Util qw(max min);

$Getopt::Std::STANDARD_HELP_VERSION = 1; # Make --help and --version flags halt execution
$main::VERSION = '0.2';

my %options = ();
getopts('dsc:g:e:n:',\%options);

my $debug = defined $options{'d'};
my $support = defined $options{'s'};

&main::HELP_MESSAGE() && die unless defined $ARGV[0] || $support;

my $caddfile = (defined $options{'c'} ? $options{'c'} : '/gscmnt/gc2719/halllab/src/gemini/data/whole_genome_SNVs.tsv.compressed.gz');
my $genefile = (defined $options{'g'} ? $options{'g'} : 'refGene.genes.b37.bed');
my $exonfile = (defined $options{'e'} ? $options{'e'} : 'refGene.exons.b37.bed');
#my $geneanncolumn = (defined $options{'m'} && defined $genefile ? $options{'m'} : 5);
my $exonanncolumn = (defined $options{'n'} && defined $exonfile ? $options{'n'} : 5);
#warn "Gene annotation column provided without nonstandard gene annotation file - defaulting to standard gene annotation file" if defined $options{'m'} && !defined $options{'g'};
#die "Nonstandard gene annotation file without column number - rerun with -m option" if !defined $options{'m'} && defined $options{'g'};
warn "Exon annotation column provided without nonstandard exon annotation file - defaulting to standard exon annotation file" if defined $options{'n'} && !defined $options{'e'};
die "Nonstandard exon annotation file without column number - rerun with -n option" if !defined $options{'n'} && defined $options{'e'};
my $compressed = ($ARGV[0] =~ /\.gz$/) if defined $ARGV[0];
my $preprocessedfile;

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

# Create preprocessing command - annotation is done without normalization because REF and ALT nucleotides are not included in VCFs describing SVs
  if (defined $ARGV[0]) {
  print STDERR "Preparing preprocessing command\n" if $debug;
  my ($prefix) = ($ARGV[0] =~ /^(?:.*\/)?(.*)\.vcf(?:\.gz)?$/);
  $preprocessedfile = "$prefix.preprocess.vcf";
  die "Please rename ${prefix}header so SVScore does not overwrite it" if -s "${prefix}header";
  my $preprocess = ($compressed ? "z": "") . "cat $ARGV[0] | awk '\$0~\"^#\" {print \$0; next } { print \$0 | \"sort -k1,1 -k2,2n\" }' | vcfanno -ends conf.toml - > $prefix.ann.vcf; grep '^#' $prefix.ann.vcf > ${prefix}header"; # Sort, annotate, grab header
  print STDERR "Preprocessing command 1: $preprocess\n" if $debug;
  die "Preprocessing failed: $!" if system($preprocess);

  # Update header
  open(HEADER, "${prefix}header") || die "Could not open ${prefix}header: $!";
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

  push @newheader, "##INFO=<ID=SVSCORE_LEFT,Number=1,Type=Float,Description=\"Maximum C score in left breakend of structural variant\">\n";
  push @newheader, "##INFO=<ID=SVSCORE_RIGHT,Number=1,Type=Float,Description=\"Maximum C score in right breakend of structural variant\">\n";
  push @newheader, "##INFO=<ID=SVSCORE_SPAN,Number=1,Type=Float,Description=\"Maximum C score in outer span of structural variant\">\n";
  push @newheader, "##INFO=<ID=SVSCORE_LTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of left breakend to end of gene\">\n";
  push @newheader, "##INFO=<ID=SVSCORE_RTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of right breakend to end of gene\">\n";
  push @newheader, @oldheader;
  open(HEADER, "> ${prefix}header") || die "Could not open ${prefix}header: $!";
  foreach (@newheader) {
    print HEADER;
  }

  my $preprocess2 = "grep -v '^#' $prefix.ann.vcf | sort -k 3,3 | cat ${prefix}header - > $preprocessedfile; rm -f ${prefix}header"; # Sort by ID, add new header, clean up

# Execute preprocessing command
  print STDERR "Preprocessing command 2: $preprocess2\n" if $debug;
  die "Preprocessing2 failed: $!" if system("$preprocess2");
  unlink "$prefix.ann.vcf" || warn "Could not delete $prefix.ann.vcf: $!" unless $debug;
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

#my %processedids = ();
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
  my ($spanexongenenames,$spangenenames,$leftexongenenames,$leftgenenames,$rightexongenenames,$rightgenenames,$cipos,$ciend) = getfields($info,"ExonGeneNames","Gene","left_ExonGeneNames","left_Gene","right_ExonGeneNames","right_Gene","CIPOS","CIEND");
  my @leftgenenames = split(/\|/,$leftgenenames);
  my @rightgenenames = split(/\|/,$rightgenenames);
  my ($leftintrons,$rightintrons) = getfields($info,"left_Intron","right_Intron") if $svtype eq 'BND' || $svtype eq 'INV';


  if ($svtype eq 'BND') {
    my ($mateid) = ($id =~ /^(\d+)_(?:1|2)/);
    my ($nextmateid) = ($i == $#inputlines ? ("") : ((split(/\s+/,$inputlines[$i+1]))[2] =~ /^(\d+)_(?:1|2)/)); # $nextmateid is set to the mateid of the next line, unless the current line is the final line in the file. In this case, $nextmateid is set to the empty string because we know it's not the first mate of a consecutive pair
    my $firstmate = ($nextmateid && $mateid == $nextmateid); # First mate of a consecutive pair
    my $secondmate = ($lastbndmateid && $mateid == $lastbndmateid); # Second mate of a consecutive pair
    $singletonbnd = !($firstmate || $secondmate);
    if ($secondmate || $singletonbnd) { # Is this the second mate of this BND seen or a singleton BND? If the second of a pair, get annotations for first mate. If a singleton, get rid of one set of annotations
      $lastbndmateid = "";
      unless ($singletonbnd) {
	my $mateline = $inputlines[$i-1];
	@splitmateline = split(/\s+/,$mateline);
	$mateinfo = $splitmateline[7];
      }
      if ($info =~ /SECONDARY/) { # Current line is secondary (right), so mate is primary (left)
	($leftexongenenames, $leftgenenames, $leftintrons) = ($singletonbnd ? ("","","") : getfields($mateinfo,"ExonGeneNames","Gene","Intron"));
	($rightchrom, $rightpos) = ($leftchrom, $leftpos);
	($leftchrom,$leftpos) = ($alt =~ /([\w.]+):(\d+)/);
	($cipos, $ciend) = ($ciend, $cipos); # Switch $cipos and $ciend
      } else { # Current line is primary (left), so mate is secondary (right)
	($rightexongenenames, $rightgenenames, $rightintrons) = ($singletonbnd ? ("","","") : getfields($mateinfo,"ExonGeneNames","Gene","Intron"));
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

  # Calculate start/stop coordinates from POS and CI
  unless ($svtype eq "BND") {
    $rightchrom = $leftchrom;
    $rightpos = getfields($info,"END");
  }
  @cipos = split(/,/,$cipos);
  $leftstart = $leftpos + $cipos[0] - 1;
  $leftstop = $leftpos + $cipos[1];
  @ciend = split(/,/,$ciend);
  $rightstart = $rightpos + $ciend[0] - 1;
  $rightstop = $rightpos + $ciend[1];

  my ($leftscore, $rightscore);

  # Calculate maximum C score depending on SV type
  if ($svtype eq "DEL" || $svtype eq "DUP") {
    my $spanscore = ($rightstop - $leftstart > 1000000 ? 100 : maxcscore($caddfile, $leftchrom, $leftstart, $rightstop));
    $leftscore = maxcscore($caddfile, $leftchrom, $leftstart, $leftstop);
    $rightscore = maxcscore($caddfile, $rightchrom, $rightstart, $rightstop);
    $info .= ";SVSCORE_SPAN=$spanscore;SVSCORE_LEFT=$leftscore;SVSCORE_RIGHT=$rightscore";
    undef $spanscore;
  } elsif ($svtype eq "INV" || $svtype eq "BND") {
    $leftscore = maxcscore($caddfile, $leftchrom, $leftstart, $leftstop);
    $rightscore = maxcscore($caddfile, $rightchrom, $rightstart, $rightstop);
    my ($sameintrons,@lefttruncationscores,@righttruncationscores,$lefttruncationscore,$righttruncationscore,%leftintrons,@rightintrons) = ();
    unless ($singletonbnd) {
      %leftintrons = map {$_ => 1} (split(/\|/,$leftintrons));
      @rightintrons = split(/\|/,$rightintrons);
      ## At worst, $leftintrons and $rightintrons are lists of introns. The only case in which the gene is not disrupted is if both lists are equal and nonempty, meaning that in every gene hit by this variant, both ends of the variant are confined to the same intron
      $sameintrons = scalar (grep {$leftintrons{$_}} @rightintrons) == scalar @rightintrons && scalar @rightintrons > 0;
    }
    if (($leftgenenames || $rightgenenames) && ($singletonbnd || !$sameintrons)) { # Gene is being truncated - left or right breakend hits a gene, and the breakends are not confined to the same introns (or, if the variant is a singleton BND, this latter condition is not necessary)
      foreach my $gene (split(/\|/,$leftgenenames)) {
	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$leftchrom}}[0..2];	
	die "$id, $gene" unless defined $genestrand; ## DEBUG
	if ($genestrand eq '+') {
	  push @lefttruncationscores, maxcscore($caddfile, $leftchrom, max($genestart,$leftstart), $genestop); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
	} else { ## Minus strand
	  push @lefttruncationscores,maxcscore($caddfile, $leftchrom, $genestart, min($genestop,$leftstop)); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
	}
      }
      $lefttruncationscore = max(@lefttruncationscores) if @lefttruncationscores;
      foreach my $gene (split(/\|/,$rightgenenames)) {
	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}->{$rightchrom}}[0..2];	
	die "$id, $rightgenenames, $gene" unless defined $genestrand; ## DEBUG
	if ($genestrand eq '+') {
	  push @righttruncationscores, maxcscore($caddfile, $rightchrom, max($genestart,$rightstart), $genestop); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
	} else { ## Minus strand
	  push @righttruncationscores,maxcscore($caddfile, $rightchrom, $genestart, min($genestop,$rightstop)); # Start from beginning of gene, stop at end of gene or breakend, whichever is further right (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
	}
	($genestart,$genestop,$genestrand) = (); # Get rid of old variables
      }
      $righttruncationscore = max(@righttruncationscores) if @righttruncationscores;
    }
    ($sameintrons, %leftintrons,@rightintrons) = (); # Get rid of old variables
    $info .= ";SVSCORE_LEFT=$leftscore;SVSCORE_RIGHT=$rightscore" . (defined $lefttruncationscore ? ";SVSCORE_LTRUNC=$lefttruncationscore" : "") . (defined $righttruncationscore ? ";SVSCORE_RTRUNC=$righttruncationscore" : "");
    $mateinfo .= ";SVSCORE_LEFT=$leftscore;SVSCORE_RIGHT=$rightscore" . (defined $lefttruncationscore ? ";SVSCORE_LTRUNC=$lefttruncationscore" : "") . (defined $righttruncationscore ? ";SVSCORE_RTRUNC=$righttruncationscore" : "") if $svtype eq "BND" && !$singletonbnd;
    (@lefttruncationscores,@righttruncationscores,$lefttruncationscore,$righttruncationscore) = (); # Get rid of old variables
  } elsif ($svtype eq "INS") { # leftscore is base before insertion, rightscore is base after insertion
    $leftscore = maxcscore($caddfile, $leftchrom, $leftstart-1, $leftstart-1);
    $rightscore = maxcscore($caddfile, $rightchrom, $rightstart+1, $rightstart+1);
    $info .= ";SVSCORE_LEFT=$leftscore;SVSCORE_RIGHT=$rightscore";
  } else {
    die "Unrecognized SVTYPE $svtype at line " . $i+1 . " of annotated VCF file\n";
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

sub maxcscore { # Calculate maximum C score within a given region using CADD data
  my ($filename, $chrom, $start, $stop) = @_;
  my @scores = ();
  my $tabixoutput = `tabix $filename $chrom:$start-$stop`;
  my @tabixoutputlines = split(/\n/,$tabixoutput);
  foreach my $taboutline (@tabixoutputlines) {
    push @scores, split(/,/,(split(/\s+/,$taboutline))[4]);
  }
  
  my $max = max(@scores);
  $max = -1 unless defined $max;
  ($filename,$chrom,$start,$stop,@scores,$tabixoutput,@tabixoutputlines) = (); # Get rid of old variables
  return $max;
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
  print "usage: ./svscore.pl [-ds] [-c caddfile] [-g genefile] [-e exonfile] [-n exonannotationcolumn] vcf
    -d	      Debug (verbose) mode, keeps intermediate and supporting files
    -s	      Create/download supporting files and quit
    -c	      Used to point to whole_genome_SNVs.tsv.gz
    -g	      Used to point to gene BED file
    -e	      Used to point to exon BED file
    -n	      Column number for annotation in exon BED file to be added to VCF

    --help    Display this message
    --version Display version\n";
}

sub main::VERSION_MESSAGE() {
  print "SVScore version $main::VERSION\n";
}
