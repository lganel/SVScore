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

eval {
  my %possibleoperations = ("MAX", 1, "SUM", 1, "TOP", 1, "MEAN", 1, "MEANWEIGHTED", 1, "TOPWEIGHTED", 1); # Hash of supported operations
  my %types = map {$_ => 1} ("DEL", "DUP", "INV", "BND", "CNV", "MEI", "INS"); # Hash of supported svtypes
  my %truncationtypes = map {$_ => 1} ("INV","INS","DEL"); # All svtypes for which truncation scores are calculated
  my %intervals = ("LEFT", 0, "RIGHT", 1, "SPAN", 2, "LTRUNC", 3, "RTRUNC", 4); # Hash of supported intervals

  my %options = ();
  getopts('dvc:g:e:m:o:p:i:n:',\%options);

# Parse command line options, set variables, check input parameters
  my $debug = defined $options{'d'};
  my $verbose = defined $options{'v'};

  my $inputfile  = $options{'i'};
  &main::HELP_MESSAGE() && die unless defined $inputfile;
  die "Could not find input file $inputfile" unless (-s $inputfile || $inputfile eq "stdin");
  my $caddfile = (defined $options{'c'} ? $options{'c'} : 'whole_genome_SNVs.tsv.gz');
  die "Could not find CADD file $caddfile" unless -s $caddfile;
  my $genefile = (defined $options{'g'} ? $options{'g'} : 'refGene.genes.b37.bed.gz');
  my $geneanncolumn = (defined $options{'m'} ? $options{'m'} : 4);
  my $exonfile = (defined $options{'e'} ? $options{'e'} : 'refGene.exons.b37.bed');
  my $exonanncolumn = (defined $options{'n'} ? $options{'n'} : ($exonfile eq 'refGene.exons.b37.bed' ? 5 : 4));
  my $genestrandcolumn = (defined $options{'p'} ? $options{'p'} : 5);
  $options{'o'} =~ tr/[a-z]/[A-Z]/ if defined $options{'o'};
  my $ops = (defined $options{'o'} ? $options{'o'} : 'TOP10WEIGHTED');

  die "Gene annotation column cannot equal gene strand column" if $geneanncolumn==$genestrandcolumn;
  if ($geneanncolumn <= 3) {
    $geneanncolumn = 4;
    warn "Gene annotation column must be greater than 3 - defaulting to standard gene annotation column (4)";
  }
  if ($genestrandcolumn <= 3) {
    $genestrandcolumn = 5;
    warn "Gene annotation column must be greater than 3 - defaulting to standard gene strand column (5)";
  }
  if ($exonanncolumn <= 3) {
    $exonanncolumn = 5;
    warn "Exon annotation column must be greater than 3 - defaulting to standard exon annotation column (5)";
  }
  if (defined $options{'m'} && !defined $options{'g'}) {
    $geneanncolumn = 4;
    warn "Custom gene annotation column provided without custom gene file - defaulting to standard gene annotation column (4)";
  }
  if (defined $options{'p'} && !defined $options{'g'}) {
    $genestrandcolumn = 5;
    warn "Custom gene strand column provided without custom gene annotation file - defaulting to standard gene strand column (5)";
  }
  if (defined $options{'n'} && !defined $options{'e'}) {
    $exonanncolumn = 5;
    warn "Custom exon annotation column provided without custom exon file - defaulting to standard exon annotation column (5)";
  }

  my @ops = uniq(split(/,/,$ops));
  my %operations = (); # Hash of chosen operations
  foreach my $i (0..$#ops) { # Populate %operations with chosen operations given by -o
    my $op = $ops[$i];
    $op =~ s/\d+//g; # Get rid of all numbers in $op for lookup in %possibleoperations
    unless (defined $possibleoperations{$op} && ($ops[$i] !~ /TOP/ || $ops[$i] =~ /^TOP\d/)) { ## Make sure all operations are recognized (and that all TOP operations have a number)
      warn "Unrecognized operation specified: $ops[$i]\n";
      &main::HELP_MESSAGE() && die;
    }
    $operations{$ops[$i]} = $i;
  }

  my $compressed = ($inputfile =~ /\.gz$/);
  my ($uncompressedgenefile, $uncompressedexonfile, $compressedexonfile, $alteredgenefile, $alteredexonfile, $headerfile, $sortedfile, $preprocessedfile, $bedpeout, $vcfout);

# Set up all necessary preprocessing to be taken care of before analysis can begin. This includes decompression, annotation using vcfanno, and generation of intron/exon/gene files, whichever are necessary
  if ($exonfile eq 'refGene.exons.b37.bed' && !-s $exonfile) { # Generate exon file if necessary
    print STDERR "Generating exon file: $exonfile\n" if $debug;
    system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); n=int(\$9); split(\$10,start,\",\");split(\$11,end,\",\"); for(i=1;i<=n;++i) {print \$3,start[i],end[i],\$2\".\"i,\$13,\$2; } }' OFS=\"\t\" | sort -k 1,1V -k 2,2n | uniq > refGene.exons.b37.bed");
  } elsif ($exonfile ne 'refGene.exons.b37.bed' && !-s $exonfile) {
    die "Exon file $exonfile not found or empty!";
  }

  if ($genefile eq 'refGene.genes.b37.bed.gz' && !-s $genefile) { # Generate gene file if necessary
    print STDERR "Generating gene file: $genefile\n" if $debug;
    system("curl -s \"http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz\" | gzip -cdfq | awk '{gsub(\"^chr\",\"\",\$3); print \$3,\$5,\$6,\$13,\$4}' OFS=\"\\t\" | sort -k 1,1V -k 2,2n | uniq | bgzip -c > refGene.genes.b37.bed.gz");
  } elsif ($genefile ne 'refGene.genes.b37.bed.gz' && !-s $genefile) {
    die "Gene file $genefile not found or empty!";
  }
  my ($geneprefix) = ($genefile =~ /^(?:.*\/)?(.*)\.bed(?:\.gz)?$/);
  my ($exonprefix) = ($exonfile =~ /^(?:.*\/)?(.*)\.bed(?:\.gz)?$/);
  my $intronfile = "svscoretmp/introns.$geneprefix.$exonprefix.bed.gz";

  my @todelete = (); # List of files to delete at end of run

# Zip/unzip annotation files as necessary
  if ($genefile =~ /\.gz/) {
    ($uncompressedgenefile) = ($genefile =~ /(.*)\.gz$/);
    unless (-s $uncompressedgenefile) {
      push @todelete, $uncompressedgenefile unless $debug;
      print STDERR "Unzipping $genefile\n" if $debug;
      if (system("zcat $genefile > $uncompressedgenefile")) {
	die "Unzipping $genefile failed;;" . join(',',@todelete);
      }
    }
  } else {
    $uncompressedgenefile = $genefile;
  }

  if ($exonfile =~ /\.gz/) {
    ($uncompressedexonfile) = ($exonfile =~ /(.*)\.gz$/);
    unless (-s $uncompressedexonfile) {
      push @todelete, $uncompressedexonfile unless $debug;
      print STDERR "Unzipping $exonfile\n" if $debug;
      if (system("zcat $exonfile > $uncompressedexonfile")) {
	die "Unzipping $exonfile failed;;" . join(',',@todelete);
      }
    }
  } else {
    $uncompressedexonfile = $exonfile;
    $compressedexonfile = "$exonfile.gz";
    unless(-s $compressedexonfile) {
      push @todelete, $compressedexonfile unless $debug;
      print STDERR "Zipping $exonfile\n" if $debug;
      if (system("bgzip -c $exonfile > $compressedexonfile")) {
	die "Compressing $exonfile failed;;" . join(',',@todelete);
      }
    }
  }

  mkdir "svscoretmp" unless -d "svscoretmp";

  unless (-s "$intronfile") { # Generate intron file if necessary - add column with unique intron ID equal to line number (assuming intron file has no header line) and sort
    push @todelete, $intronfile unless $debug;
    print STDERR "Generating intron file\n" if $debug;
    system("bedtools subtract -a $uncompressedgenefile -b $uncompressedexonfile | sort -u -k1,1V -k2,2n | awk '{print \$0 \"\t\" NR}' | bgzip -c > $intronfile");
  }

# Use tabix to index the annotation files
  unless (-s "$intronfile.tbi") {
    push @todelete, "$intronfile.tbi" unless $debug;
    print STDERR "Tabix indexing $intronfile\n" if $debug;
    if(system("tabix -fp bed $intronfile")) {
      die "Tabix failed on $intronfile;;" . join(',',@todelete);
    }
  }
  unless (-s "$compressedexonfile.tbi") {
    push @todelete, "$compressedexonfile.tbi" unless $debug;
    print STDERR "Tabix indexing $compressedexonfile\n" if $debug;
    if(system("tabix -fp bed $compressedexonfile")) {
      die "Tabix failed on $compressedexonfile;;" . join(',',@todelete);
    }
  }

  my $intronnumcolumn = `zcat $intronfile | head -n 1 | awk '{print NF}'`; # Figure out which column has the intron number
  $intronnumcolumn = 6 unless $intronnumcolumn;
  chomp $intronnumcolumn;

# Write conf.toml file
  print STDERR "Writing config file\n" if $debug;
  open(CONFIG, "> conf.toml") || die "Could not open conf.toml: $!;;" . join(',',@todelete);
  push @todelete, "conf.toml" unless $debug;
  print CONFIG "[[annotation]]\nfile=\"$compressedexonfile\"\nnames=[\"ExonGene\"]\ncolumns=[$exonanncolumn]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"$intronfile\"\nnames=[\"Intron\",\"IntronGene\"]\ncolumns=[$intronnumcolumn,4]\nops=[\"concat\",\"concat\"]\n";
  close CONFIG;

# Create first preprocessing command - annotation is done without normalization because REF and ALT nucleotides are not included in VCFs describing SVs
  my ($prefix,$tempfile);
  my $time = gettimeofday();
  if ($inputfile eq "stdin") { # Write standard input to temp file if input file comes from STDIN
    $tempfile = "svscoretmp/stdin$time.vcf";
    open(TEMP,">$tempfile") || die "Could not open $tempfile for writing; $!;;" . join(',',@todelete);
    push @todelete, $tempfile unless $debug;
    print TEMP <STDIN>;
    close TEMP;
  }

  print STDERR "Preparing preprocessing command\n" if $debug;
  if ($inputfile eq "stdin") {
    $prefix = "stdin";
    $inputfile = $tempfile;
  } else {
    ($prefix) = ($inputfile =~ /^(?:.*\/)?(.*)\.vcf(?:\.gz)?$/);
  }

# Tag intermediate files with timestamp to avoid collisions
  $preprocessedfile = "svscoretmp/$prefix$time.preprocess.bedpe";
  $sortedfile = "svscoretmp/$prefix$time.sort.vcf.gz";
  if ($compressed) {
    if (system("gunzip -c $inputfile > svscoretmp/$prefix$time.vcf")) {
      die "Could not unzip $inputfile;;" . join(',',@todelete);
    } else {
      $inputfile = "svscoretmp/$prefix$time.vcf";
      push @todelete, $inputfile unless $debug;
    }
  }

# Make sure header is present in input file
  open(HEADERCHECK, "< $inputfile") || die "Could not open $inputfile: $!;;" . join(',',@todelete);
  my $firstline = <HEADERCHECK>;
  close HEADERCHECK;
  die "***Missing header on input file;;" . join(',',@todelete) unless($firstline =~ /^#/);

  my $preprocess = "awk '\$0~\"^#\" {print \$0; next } { print \$0 | \"sort -k1,1V -k2,2n\" }' $inputfile | bgzip -c > $sortedfile; vcfanno -ends conf.toml $sortedfile | perl reorderheader.pl stdin $inputfile | svtools vcftobedpe > $preprocessedfile"; # Sort, annotate, reorder header, convert to BEDPE
  push @todelete, $preprocessedfile, $sortedfile unless $debug;
  print STDERR "Preprocessing command:\n$preprocess\n" if $debug;
  die "Preprocessing failed;;" . join(',',@todelete) if (system($preprocess) || -z $preprocessedfile);

  $bedpeout = "svscoretmp/$prefix$time.out.bedpe";
  $vcfout = "svscoretmp/$prefix$time.out.vcf";
  open(IN, "< $preprocessedfile") || die "Could not open $preprocessedfile: $!;;" . join(',',@todelete);
  open(OUT, "> $bedpeout") || die "Could not open output file: $!;;" . join(',',@todelete);
  push @todelete, $bedpeout unless $debug;

# Update header
  open(HEADER, "grep \"^#\" $preprocessedfile |") || die "Error grabbing header: $!;;" . join(',',@todelete);
  my @newheader = ();
  my @oldheader = <HEADER>;
  close HEADER;

  if ($ops =~ /WEIGHTED/ && !(grep {/^##INFO=<ID=PRPOS,/} @oldheader)) { # If an op is weighted, but PRPOS is absent from the header, switch to unweighted operations with a warning
    warn "*****PRPOS not found in header - switching to unweighted operations*****\n";
    $ops =~ s/WEIGHTED//g;
    @ops = uniq(split(/,/,$ops));
    $ops = join(",",@ops);
    foreach my $op (keys %operations) {
      if ($op =~ /WEIGHTED/) {
	my $value = $operations{$op};
	delete $operations{$op};
	$op =~ s/WEIGHTED//;
	$operations{$op} = $value;
      }
    }
    # Close gaps in values(%operations) created by deleting entries
    my %revops = reverse %operations; # $index => $operation
    my $offset = 0; # Size of gap to close
    my %newops;
    foreach my $index (sort keys %revops) {
      $offset++ unless $index == 0 || exists $revops{$index-1}; # Increase offset if there is a gap
      $newops{$revops{$index}} = $index-$offset; # Close gap in %newops
    }
    %operations = %newops;
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
    push @newheader, "##INFO=<ID=SVSCOREMAX,Number=1,Type=Float,Description=\"Maximum of SVSCOREMAX fields of structural variant\">\n";
    if ($verbose) {
      push @newheader, "##INFO=<ID=SVSCOREMAX_LEFT,Number=1,Type=Float,Description=\"Maximum C score in left breakend of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMAX_RIGHT,Number=1,Type=Float,Description=\"Maximum C score in right breakend of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMAX_SPAN,Number=1,Type=Float,Description=\"Maximum C score in span of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMAX_LTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of left breakend to end of truncated gene\">\n";
      push @newheader, "##INFO=<ID=SVSCOREMAX_RTRUNC,Number=1,Type=Float,Description=\"Maximum C score from beginning of right breakend to end of truncated gene\">\n";
    }
  }
  if (defined $operations{"SUM"}) {
    push @newheader, "##INFO=<ID=SVSCORESUM,Number=1,Type=Float,Description=\"Maximum of SVSCORESUM fields of structural variant\">\n";
    if ($verbose) {
      push @newheader, "##INFO=<ID=SVSCORESUM_LEFT,Number=1,Type=Float,Description=\"Sum of C scores in left breakend of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCORESUM_RIGHT,Number=1,Type=Float,Description=\"Sum of C scores in right breakend of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCORESUM_SPAN,Number=1,Type=Float,Description=\"Sum of C scores in outer span of structural variant\">\n";
      push @newheader, "##INFO=<ID=SVSCORESUM_LTRUNC,Number=1,Type=Float,Description=\"Sum of C scores from beginning of left breakend to end of truncated gene\">\n";
      push @newheader, "##INFO=<ID=SVSCORESUM_RTRUNC,Number=1,Type=Float,Description=\"Sum of C scores from beginning of right breakend to end of truncated gene\">\n";
    }
  }
  foreach my $op (@ops) {
    my $weighted = ($op =~ /WEIGHTED$/);
    my ($n) = ($op =~ /^TOP(\d+)/);
    if ($op =~ /^TOP\d+/ || $op =~ /^MEAN/) {
      push @newheader, "##INFO=<ID=SVSCORE$op,Number=1,Type=Float,Description=\"Maximum of SVSCORE$op fields of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      if ($verbose) {
	push @newheader, "##INFO=<ID=SVSCORE${op}_LEFT,Number=1,Type=Float,Description=\"Mean of " . ($n ? "top $n " : "") . "C scores in left breakend of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
	push @newheader, "##INFO=<ID=SVSCORE${op}_RIGHT,Number=1,Type=Float,Description=\"Mean of " . ($n ? "top $n " : "") . "C scores in right breakend of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
	push @newheader, "##INFO=<ID=SVSCORE${op}_SPAN,Number=1,Type=Float,Description=\"Mean of " . ($n ? "top $n " : "") . "C scores in outer span of structural variant" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
	push @newheader, "##INFO=<ID=SVSCORE${op}_LTRUNC,Number=1,Type=Float,Description=\"Mean of " . ($n ? "top $n " : "") . "C scores from beginning of left breakend to end of truncated gene" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
	push @newheader, "##INFO=<ID=SVSCORE${op}_RTRUNC,Number=1,Type=Float,Description=\"Mean of " . ($n ? "top $n " : "") . "C scores from beginning of right breakend to end of truncated gene" . ($weighted ? ", weighted by probability distribution" : "") . "\">\n";
      }
    }
  }
  push @newheader, @oldheader;
  foreach (uniq(@newheader)) {
    print OUT;
  }

  print STDERR "Reading gene list\n" if $debug;
  my %genes = (); # Symbol => (Chrom => (chrom, start, stop, strand)); Hash of hashes of arrays
  open(GENES, "$uncompressedgenefile") || die "Could not open $genefile: $!;;" . join(',',@todelete);
  foreach my $geneline (<GENES>) { # Parse gene file, recording the chromosome, strand, and first and last bases (in VCF coordinates) of each gene
    my ($genechrom, $genestart, $genestop, $genesymbol, $genestrand) = (split(/\s+/,$geneline))[0..2,$geneanncolumn-1,$genestrandcolumn-1];
    $genestart++;
    if (defined $genes{$genesymbol}->{$genechrom}) { ## Assume strand stays constant
      $genes{$genesymbol}->{$genechrom}->[0] = min($genes{$genesymbol}->{$genechrom}->[0], $genestart);
      $genes{$genesymbol}->{$genechrom}->[1] = max($genes{$genesymbol}->{$genechrom}->[1], $genestop);
    } else {
      $genes{$genesymbol}->{$genechrom} = [$genestart, $genestop, $genestrand];
    }
  }
  close GENES;

# Perform sanity check to make sure annotation columns are specified correctly
  open(EXONS, "$uncompressedexonfile") || die "Could not open $exonfile: $!;;" . join(',',@todelete);

  my $count = 1;
  my $genesfound = 0;
  while(defined(my $line = <EXONS>) && $count <= 10) {
    $count++;
    my $line = <EXONS>;
    my $symbol = (split(/\s+/,$line))[$exonanncolumn-1];
    if (exists $genes{$symbol}) {
      $genesfound++;
      last;
    }
  }
  close EXONS;
  warn "************First 10 genes in exon file not found in genes file. Are your annotation columns (-m and -n) specified correctly?" unless $genesfound;

  print STDERR "Entering loop\n" if $debug;
  while (my $inputline = <IN>) {
    if ($inputline =~ /^#/) {
      next;
    }
    # Parse line
    my @splitline = split(/\s+/,$inputline);
    my ($leftchrom, $leftstart, $leftstop, $rightchrom, $rightstart, $rightstop, $id, $svtype, $info_a, $info_b) = @splitline[0..6, 10, 18, 19];
    my ($pos,$end) = getfields($info_a,"POS","END");

    $svtype = substr($svtype, 0, 3);
    unless (exists $types{$svtype}) {
      warn "Unrecognized SVTYPE $svtype at line ",$.," of preprocessed VCF file\n";
      print OUT $inputline;
      next;
    }
    my ($probleft,$probright);
    ($probleft,$probright) = getfields($info_a,"PRPOS","PREND") if $ops =~ /WEIGHTED/;
    $probleft = "" unless defined $probleft;
    $probright = "" unless defined $probright;

    my %scores = (); # Interval => List of scores by op; e.g. (LEFT => (MAXLEFT, SUMLEFT, TOP100LEFT, MEANLEFT), RIGHT => (MAXRIGHT, SUMRIGHT, TOP100RIGHT, MEANRIGHT))

    $scores{"LEFT"} = cscoreop($caddfile, $ops, $leftchrom, $leftstart, $leftstop, $probleft);
    $scores{"RIGHT"} = cscoreop($caddfile, $ops, $rightchrom, $rightstart, $rightstop, $probright);
    
    if ($svtype eq "DEL" || $svtype eq "DUP" || $svtype eq "CNV") {
      if ($rightstop - $leftstart > 1000000) {
	$scores{"SPAN"} = (100) x @ops;
      } else {
	$scores{"SPAN"} = cscoreop($caddfile, $ops, $leftchrom, $pos, $end, -1);
      }
    }

    # Calculate truncation scores
    if (exists $truncationtypes{$svtype}) {
      
      # Get vcfanno annotations
      my ($leftexongenenames, $rightexongenenames, $leftintrons, $rightintrons, $leftintrongenenames, $rightintrongenenames);
      if ($info_b eq ".") { # Single line variant in VCF
	($leftexongenenames,$rightexongenenames) = getfields($info_a,"left_ExonGene","right_ExonGene");
	($leftintrons,$rightintrons) = getfields($info_a,"left_Intron","right_Intron");
	($leftintrongenenames,$rightintrongenenames) = getfields($info_a,"left_IntronGene","right_IntronGene");
      } else { # Multiline variant in VCF (possibly only one line of variant present)
	$leftexongenenames = getfields($info_a,"ExonGene");
	$leftintrongenenames = getfields($info_a,"IntronGene");
	$leftintrons = getfields($info_a,"Intron");
	$rightintrons = getfields($info_b,"Intron");
	$rightexongenenames = getfields($info_b,"ExonGene");
	$rightintrongenenames = getfields($info_b,"IntronGene");
      }

      # We want to make sure we don't consider any variant which is contained within an intron. To do this, we want to cancel out any introns present in both $leftintrons and $rightintrons. However, because a breakend's confidence interval may extend into multiple introns, we can't just exclude any gene with an intron in both lists, or we would miss some truly gene-truncating variants. So, we start by considering truncated the genes whose exons are affected ($leftexongenenames and $rightexongenenames). We then count the number of instances of each gene in $leftintrongenenames and $rightintrongenenames. Next, we reduce the number of instances of each gene by 1 for each of its introns which occurs in both $leftintrons and $rightintrons, as this means this variant is confined to that intron, so that intron does not provide evidence of a truncation. Any genes left with at least one instance in %leftintrongenenames are then considered truncated, and are added back to %lefttruncatedgenes and %righttruncatedgenes. Genes contained within the span of a DEL are not considered truncating because they are already captured by the SPAN score
      # %lefttruncatedgenes and %righttruncatedgenes are {gene => 1}
      # %leftintrongenenames and %rightintrongenenames are {gene => # of instances in $leftgenenames or $rightgenenames}
      # %leftintrons and %rightintrons are {number of intron hit by the left or right breakend => name of gene}
      
      my %lefttruncatedgenes = map {$_ => 1} (split(/\|/,$leftexongenenames));
      my %righttruncatedgenes = map {$_ => 1} (split(/\|/,$rightexongenenames));

      my (%leftintrongenenames,%rightintrongenenames);
      my @leftintrons = split(/\|/,$leftintrons);
      my @leftintrongenenames = split(/\|/,$leftintrongenenames);
      my %leftintrons = map {$leftintrons[$_] => $leftintrongenenames[$_]} (0..$#leftintrons); # @leftintrons and @leftintrongenes should have the same number of elements if vcfanno is working as it should
      my @rightintrons = split(/\|/,$rightintrons);
      my @rightintrongenenames = split(/\|/,$rightintrongenenames);
      my %rightintrons = map {$rightintrons[$_] => $rightintrongenenames[$_]} (0..$#rightintrons);

      foreach (@leftintrongenenames) {
	$leftintrongenenames{$_}++; # Count instances of genes in @leftintrongenenames
      }
      foreach (@rightintrongenenames) {
	$rightintrongenenames{$_}++;
      }

      foreach my $intron (@leftintrons) { # Cancel out introns hit by both right and left breakends
	if (exists $rightintrons{$intron}) {
	  my $introngene = $leftintrons{$intron};
	  $leftintrongenenames{$introngene}--;
	  delete $leftintrongenenames{$introngene} unless $leftintrongenenames{$introngene}; # If there are no more instances of $introngene in %leftgenenames, delete the entry
	  $rightintrongenenames{$introngene}--;
	  delete $rightintrongenenames{$introngene} unless $rightintrongenenames{$introngene}; # If there are no more instances of $introngene in %rightgenenames, delete the entry
	}
      }

      foreach (keys %leftintrongenenames) { # Add remaining genes in %leftintrongenenames to %lefttruncatedgenes
	$lefttruncatedgenes{$_} = 1;
      }
      foreach (keys %rightintrongenenames) { # Add remaining genes in %rightintrongenenames to %righttruncatedgenes
	$righttruncatedgenes{$_} = 1;
      }

      my @lefttruncatedgenes = keys %lefttruncatedgenes;
      my @righttruncatedgenes = keys %righttruncatedgenes;
      $scores{"LTRUNC"} = truncationscore($leftchrom, $pos, \@lefttruncatedgenes, \%genes, $caddfile, $ops, \%operations) if @lefttruncatedgenes;
      $scores{"RTRUNC"} = truncationscore($rightchrom, $end, \@righttruncatedgenes, \%genes, $caddfile, $ops, \%operations) if @righttruncatedgenes;
    }

    # This is an ugly loop which transposes %scores so that the keys are operations, not intervals
    my %scoresbyop = ();
    foreach my $interval (sort {$intervals{$a} <=> $intervals{$b}} keys %scores) { # LEFT, RIGHT, (SPAN, LTRUNC, RTRUNC)
      foreach my $op (@ops) { # MAX, SUM, TOP\d, TOP\dWEIGHTED, MEAN, MEANWEIGHTED
	push @{$scoresbyop{$op}}, $scores{$interval}->[$operations{$op}];
      }
    }

# Calculate maxes and add to info, replacing existing SVSCORE fields if they exist
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

    $splitline[18] = $info_a;
    $splitline[19] = $info_b;
    print OUT join("\t",@splitline) . "\n";

    print STDERR $.,", " if $debug;
  }

  close IN;
  close OUT;

# Convert back to vcf, sort, and add header
  system("svtools bedpetovcf -b $bedpeout > $vcfout");
  system("grep \"^#\" $vcfout > $vcfout.header");
  print `grep -v "^#" $vcfout | sort -k1,1V -k2,2n | cat $vcfout.header -`;
  push @todelete, $vcfout,"$vcfout.header";

  die ";;" . join(',',@todelete); ## Signal to the rest of the program that the eval block terminated successfully
};

my ($error,$todelete) = split(/;;/,$@);
if ($todelete) {
  $todelete = (split(/ at /,$todelete))[0];
  my @todelete = split(/,/,$todelete);
  unlink @todelete;
}

# Delete svscoretmp if it's now empty
if (-d "svscoretmp") {
  opendir(DIR,"svscoretmp");
  my @dir = readdir(DIR);
  my $dircount = @dir;
  foreach (@dir) {
    $dircount-- if ($_ eq '.' || $_ eq '..');
  }
  if ($dircount == 0) {
    if(system("rmdir svscoretmp")) {
      warn "Could not delete svscoretmp";
    }
  }
  close DIR;
}

chomp $error;
print STDERR "$error\n" if $error; # Relay error from eval to STDERR

sub cscoreop { # Apply operation(s) specified in $ops to C scores within a given region using CADD data. For left/right breakends, $start is the first possible breakpoint, and $stop is the final possible breakpoint (in BED coordinates). For span/truncation scores, [$start,$stop] is an inclusive interval of all the affected bases (in VCF coordinates). cscoreop will return -1 for any weighted operation on a variant without a PRPOS field. Giving -1 as $prpos signifies that scores should not be weighted (i.e. the score is a span/truncation score), while giving "" as $prpos signifies that the variant has no PRPOS field
  my ($filename, $ops, $chrom, $start, $stop, $prpos) = @_;
  my @ops = uniq(split(/,/,$ops));
  my (@prpos,%probdist);
  my $spanortrunc = ($prpos eq -1);
  my $weight = ($ops =~ /WEIGHTED/ && $prpos && !$spanortrunc);
  if ($weight) {
    @prpos = split(/,/,$prpos);
    foreach my $i ($start..$stop) { # Populate %probdist
      $probdist{$i} = $prpos[$i-$start];
    }
  }
  
  my (%bptscores,$res) = (); # %bptscores = {BEDcoordinate => Possiblebreakpointscore}
  $stop++ unless $spanortrunc; ## Increment end coordinate for tabix so that we capture the CADD score for the base following the interval to allow for calculation of a score for the final possible breakpoint
  my $tabixoutput = `tabix $filename $chrom:$start-$stop`;
  my @tabixoutputlines = split(/\n/,$tabixoutput);
  
  my %basescores = (); # Hash from VCF position to list of scores for each position
  foreach my $line (@tabixoutputlines) { # Populate %basescores with tabix output
    my @split = split(/\s+/,$line);
    push @{$basescores{$split[1]}}, $split[5];
  }
  foreach my $pos (sort {$a <=> $b} keys %basescores) { # Replace values in %basescores with max score at each VCF position
    $basescores{$pos} = max(@{$basescores{$pos}});
    if (!$spanortrunc && exists $basescores{$pos-1}) { # Calculate all scores for possible breakpoints by averaging the base scores of the two flanking bases of the possible breakpoint. Place scores in %bptscores
      $bptscores{$pos-1} = ($basescores{$pos-1} + $basescores{$pos}) / 2;
      delete $basescores{$pos-1}; # Save some memory - the basescores for all bases before $pos are now useless
    }
  }

  my $scorehashref = ($spanortrunc ? \%basescores : \%bptscores); # Choose which scores hash to collapse based on whether the interval is a span/trunc or a left/right breakend
  my %scores = %{$scorehashref}; # Use reference to avoid copying hash directly. %scores is %basescores if the interval is span/trunc, or %bptscores if the interval is left/right

  unless (@tabixoutputlines && %scores) { # Short circuit if interval is a left/right breakend and there are not enough base scores to calculate breakpoint scores (i.e. there are no 2 consecutive bases with scores in the interval), or if the interval is a span or truncation and there are no base scores
    my @res = (-1) x @ops;
    return \@res; 
  }

  # @scores is keys %basescores (if the interval is span/trunc) or keys %bptscores (if the interval is left/right). Either way, @scores is ordered by position
  my (@scores,@probdist,@weightedbptscores) = ();
  foreach my $pos (sort {$a <=> $b} keys %scores) { # Collapse %bptscores and %probdist into arrays, getting rid of positions in %probdist with no corresponding possible breakpoint scores
    push @scores, $scores{$pos};
    push @probdist, $probdist{$pos} if $weight;
  }

  if ($weight) { # Rescale probability distribution to add up to 1 (to account for excluded bases with no C scores or faulty PRPOS annotation) and weight @scores
    my $normref = normalize(\@probdist);
    @probdist = @{$normref};
    @weightedbptscores = pairwise {$a * $b} @scores, @probdist;
  }

  foreach my $op (@ops) { ## Loop through all operations in @ops, appending the resulting score to @res
    my $weightedop = ($op =~ /WEIGHTED/ && !$spanortrunc);
    if ($weightedop && !$prpos) { ## Don't calculate scores for weighted ops if PRPOS is absent
      push @{$res}, -1;
      next;
    }

    my $opscorelistref = ($weightedop ? \@weightedbptscores : \@scores); # Use a reference to avoid copying arrays.

    if ($op =~ /^TOP(\d+)/) {
      my @optopscores = @{$opscorelistref}; # @optopscores is @weightedbptscores if $op is weighted and the interval is left/right, keys %bptscores if $op is not weighted and the interval is left/right, or keys %basescores if the interval is span/trunc
      my $topn = min($1, scalar @optopscores);
      my @topnindices = (sort {$optopscores[$b] <=> $optopscores[$a]} (0..$#optopscores))[0..$topn-1]; # Get indices of the greatest $topn positions in @optopscores. These represent the positions with the highest (product of probability and) breakpoint scores
      my @topnscores = @scores[@topnindices]; # Capture unweighted scores of $topn possible breakpoints in interval

      if ($weightedop) {
	my @topnprobdist = @probdist[@topnindices];
	# Rescale probability distribution to add up to 1 and weight @topnscores
	my $topnnormref = normalize(\@topnprobdist);
	@topnprobdist = @{$topnnormref};
	@topnscores = pairwise {$a * $b} @topnscores, @topnprobdist;
      }
      $opscorelistref = \@topnscores;
    }

    my $newscore;
    my @opscores = @{$opscorelistref}; # If $op is TOP, @opscores is @topnscores (which is weighted if $op is weighted). Otherwise, @opscores is @weightedbptscores if $op is weighted and the interval is left/right, keys %bptscores if $op is not weighted and the interval is left/right, or keys %basescores if the interval is span/trunc
    if ($op eq "MAX") {
      $newscore = nearest(0.001,max(@opscores));
    } elsif ($op eq "SUM" || $weightedop) { # Compute sum of @opscores if $op is SUM, or if the op is weighted and we're on a breakend, in which case, @opscores is @weightedbptscores (or a subset of it in the TOP case), so the sum of @opscores is the weighted mean of @scores (or of the subset, if $op begins with TOP)
      $newscore = nearest(0.001,sum(@opscores));
    } elsif ($op eq "MEAN" || $op =~ /^TOP\d+$/ || ($op =~ /WEIGHTED/ && $prpos eq -1)) { # Compute mean of @opscores if $op is MEAN or TOP, or if $op is MEANWEIGHTED or TOP\d+WEIGHTED and we're on the span/truncation score
      $newscore = nearest(0.001,sum(@opscores)/scalar(@opscores));
    } else {
      die "Error: Unrecognized operation: $op"
    }
    push @{$res}, $newscore;
  }
  return $res;
}

sub normalize { # Given an array reference, normalize the array so it sums to 1 and return a reference to the array
  my @ls = @{$_[0]};
  my $sum = sum(@ls);
  unless ($sum == 1) {
    foreach (0..$#ls) {
      $ls[$_] = $ls[$_] / $sum;
    }
  }
  return \@ls;
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
  my ($chrom, $pos, $introngenesref, $genesref, $caddfile, $ops, $operationsref) = @_;
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
      $cscoreopres = cscoreop($caddfile, $ops, $chrom, max($genestart,$pos),$genestop, -1); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
    } else {
      $cscoreopres = cscoreop($caddfile, $ops, $chrom, $genestart,min($genestop,$pos), -1); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but none of the supported operations are order-dependent)
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
  print STDERR "usage: ./svscore.pl [-dv] [-o op] [-t topnumber] [-g genefile] [-m geneannotationcolumn] [-p genestrandcolumn] [-e exonfile] [-c caddfile] -i vcf
    -i	      Input VCF file. May be bgzip compressed (ending in .vcf.gz). Use \"-i stdin\" if using standard input
    -d	      Debug mode, keeps intermediate and supporting files, displays progress
    -v	      Verbose mode - show all calculated scores (left/right/span/ltrunc/rtrunc, as appropriate)
    -o	      Comma-separated list of operations to perform on CADD score intervals (must be some combination of sum, max, mean, meanweighted, top\\d, and top\\dweighted - defaults to top10weighted)
    -g	      Points to gene BED file (refGene.genes.b37.bed)
    -e	      Points to exon BED file (refGene.exons.b37.bed)
    -m	      Column number for gene name in gene BED file (4)
    -p	      Column number for strand in gene BED file (5)
    -n	      Column number for gene name in exon BED file (5 for refGene.exons.b37.bed, 4 otherwise)
    -c	      Points to whole_genome_SNVs.tsv.gz (defaults to current directory)

    --help    Display this message
    --version Display version\n"
}

sub main::VERSION_MESSAGE() {
  print "SVScore version $main::VERSION\n";
}
