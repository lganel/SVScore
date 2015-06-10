#!/usr/bin/perl -w

# Use -d for debug mode
# Use -z to indicate that the given VCF file is compressed (ending in .vcf.gz)
# Use -a to indicate that the given VCF file is already annotated
# Requires vcfanno (http://www.github.com/brentp/vcfanno)
# Requires /gscmnt/gc2719/halllab/src/gemini/data/whole_genome_SNVs.tsv.compressed.gz and its *.tbi counterpart
# Requires refGene.exons.b37.bed and humangenes.bed (from UCSC Genome Browser knownGene table - chrom, txStart, txEnd, strand, kgXref.geneSymbol in that order). Neither should have a header line
# Creates file of introns using bedtools subtract on the above two files
# Annotation done without normalization because REF and ALT nucleotides are not included in VCFs describing SVs
# Assumes input vcf is sorted

use strict;
use Getopt::Std;
use List::Util qw(max min);

my %options = ();
getopts('dza',\%options);

my $compressed = defined $options{'z'};
my $debug = defined $options{'d'};
my $annotated = defined $options{'a'};

##TODO PRIORITY 2: Enable piping input through STDIN 
##TODO PRIORITY 2: Enable option for pointing to CADD file directory
##TODO PRIORITY 1: Work with VCFs instead of BEDPEs

# Set up all necessary preprocessing to be taken care of before analysis can begin. This includes decompression, annotation using vcfanno, and generation of intron file. May be a little slower than necessary in certain situations because some arguments are supplied by piping cat output rather than supplying filenames directly.
print "Preparing preprocessing command\n" if $debug;
$ARGV[0] =~ /^(.*)\.vcf$/;
my $prefix = $1;
my $preprocessedfile;
if (!$compressed && $annotated) { # No need to copy the file if it's already preprocessed
  $preprocessedfile = $ARGV[0];
} else {
  $preprocessedfile = "$prefix.ann.vcf";
  # Create preprocessing command
  my $preprocess = ($compressed ? "z": "") . "cat $ARGV[0]" . ($annotated ? "" : " | vcfanno -ends conf.toml -") . " > $preprocessedfile";

  my $wroteconfig = 0;
  # Write conf.toml file
  unless ($annotated || -s "conf.toml") {
    print "Writing config file\n" if $debug;
    $wroteconfig = 1;
    open(CONFIG, "> conf.toml") || die "Could not open conf.toml: $!";
    print CONFIG "[annotation]]\nfile=\"refGene.exons.b37.bed\"\nnames=[\"ExonGeneNames\"]\ncolumns=[5]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"humangenes.bed\"\nnames=[\"Gene\"]\ncolumns=[5]\nops=[\"uniq\"]\n\n[[annotation]]\nfile=\"introns.bed\"\nnames=[\"Intron\"]\ncolumns=[6]\nops=[\"uniq\"]";
    close CONFIG;
  }

  # Execute preprocessing command
  print "Preprocessing command: $preprocess\n" if $debug;
  die if system("preprocess");
}

unless (-s 'introns.bed') { # Generate introns file if necessary - add column with unique intron ID equal to line number (assuming introns.bed has no header line) and sort
  print "Generating introns file\n" if $debug;
  system("bedtools subtract -a humangenes.bed -b refGene.exons.b37.bed | bedtools sort -i - | uniq | awk '{print \$0 \"\t\" NR}' > introns.bed"); # Add column with unique intron ID equal to line number, assumes introns.bed has no header line
}

print "Reading gene list\n" if $debug;
my %genes = (); # Symbol => (chrom, start, stop, strand)
open(GENES, "< humangenes.bed") || die "Could not open humangenes.bed: $!";
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

print "Entering loop\n" if $debug;
my $linenum = 1;
foreach my $vcfline (<IN>) {
  print $linenum if $debug;
  $linenum++ && next if ($vcfline =~ /^#/);

  # Parse line
  my ($leftchrom, $leftpos, $id, $alt, $info) = (split(/\s+/,$vcfline))[0..2, 4, 7];

  # Calculate coordinates
  my ($mateid,$rightchrom,$rightpos,$rightstart,$rightstop,$leftstart,$leftstop);

  $info =~ /SVTYPE=(\w+);/;
  my $svtype = $1;
  my ($cipos, $ciend) = getfields($info,"CIPOS","CIEND");
  my @cipos = split(/,/,$cipos);
  my @ciend = split(/,/,$ciend);

  if ($svtype eq 'BND') {
    $id =~ /(\d+)_(?:1|2)/;
    my $mateid = $1;
    ($rightchrom,$rightpos) = ($alt =~ /(\w):(\d+)/);
  } else {
    $rightchrom = $leftchrom;
    $rightpos = getfields($info,"END");
  }

  $leftstart = $leftpos + $cipos[0] - 1;
  $leftstop = $leftpos + $cipos[1];
  $rightstart = $rightpos + $ciend[0] - 1;
  $rightstop = $rightpos + $ciend[1];

  $linenum++ && next if $svtype eq 'BND' && exists $processedids{$mateid}; ## Skip secondary BND lines
  $processedids{$mateid} = 1;
  my ($spanscore, $leftscore, $rightscore);
  # Check for overlapping exons
  my ($spanexongenenames,$spangenenames) = getfields($info,"ExonGeneNames","Gene");

  my ($leftexongenenames,$leftgenenames) = getfields($info,"left_ExonGeneNames","left_Gene");
  my @leftgenenames = split(/,/,$leftgenenames);

  my ($rightexongenenames,$rightgenenames) = getfields($info,"right_ExonGeneNames","right_Gene");
  my @rightgenenames = split(/,/,$rightgenenames);
   
  # Calculate maximum C score depending on SV type
  if ($svtype eq "DEL" || $svtype eq "DUP") {
    warn "Unequal left and right chromosome number in $svtype at line $linenum" unless $leftchrom eq $rightchrom;
    $spanscore = maxcscore($leftchrom, $leftstart, $rightstop);
    $leftscore = maxcscore($leftchrom, $leftstart, $leftstop);
    $rightscore = maxcscore($rightchrom, $rightstart, $rightstop);
    
  } elsif ($svtype eq "INV" || $svtype eq "BND") {
    warn "Unequal left and right chromosome number in $svtype at line $linenum" unless $leftchrom eq $rightchrom || $svtype eq "BND";
    $leftscore = maxcscore($leftchrom, $leftstart, $leftstop);
    $rightscore = maxcscore($rightchrom, $rightstart, $rightstop);

    my $leftintrons = getfields($info, "left_Intron");
    my $rightintrons = getfields($info, "right_Intron");
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
	  push @lefttruncationscores, maxcscore($leftchrom, max($genestart,$leftstart), $genestop); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
	} else { ## Minus strand
	  push @lefttruncationscores,maxcscore($leftchrom, $genestart, min($genestop,$leftstop)); # Start from beginning of gene, stop at end of gene or breakend, whichever is further left (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)
	}
      }
      my @righttruncationscores = ();
      foreach my $gene (split(/,/,$rightgenenames)) {
	my ($genestart,$genestop,$genestrand) = @{$genes{$gene}}[1..3];	
	print scalar @{$genes{$gene}} if $linenum == 18; ## DEBUG
	print @{$genes{$gene}}[4] if $linenum == 18; ## DEBUG
	if ($genestrand eq '+') {
	  push @righttruncationscores, maxcscore($rightchrom, max($genestart,$rightstart), $genestop); # Start from beginning of gene or breakend, whichever is further right, stop at end of gene
	} else { ## Minus strand
	  push @righttruncationscores,maxcscore($rightchrom, $genestart, min($genestop,$rightstop)); # Start from beginning of gene, stop at end of gene or breakend, whichever is further right (this is technically backwards, but it doesn't matter for the purposes of finding a maximum C score)

	}
      }
      my $lefttruncationscore = (@lefttruncationscores ? max(@lefttruncationscores) : "NoLeftGenes");

      print "\n" if $linenum == 18 && $righttruncationscores[0] eq ""; ## DEBUG
      my $righttruncationscore = (@righttruncationscores ? max(@righttruncationscores) : "NoRightGenes");
      $spanscore = "$lefttruncationscore,$righttruncationscore";
      $spangenenames = "$svtype:IgnoredSpan";
      $spanexongenenames = "$svtype:IgnoredSpan";
    }
  } elsif ($svtype eq "INS") { # leftscore is base before insertion, rightscore is base after insertion
    $leftscore = maxcscore($leftchrom, $leftstart-1, $leftstart-1);
    $rightscore = maxcscore($rightchrom, $rightstart+1, $rightstart+1);
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
    if ($svtype eq "BND" || $svtype eq "INV") {
      $spangenenames = "$svtype:IgnoredSpan";
    } else {
      $spangenenames = "NoGenesInSpan";
    }
  }
  unless ($spanexongenenames) {
    if ($svtype eq "BND" || $svtype eq "INV") {
      $spangenenames = "$svtype:IgnoredSpan";
    } else {
      $spanexongenenames = "NoExonsInSpan";
    }
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

unlink "$preprocessedfile" unless (!$compressed && $annotated);

sub maxcscore { # Calculate maximum C score within a given region using CADD data
  my ($chrom, $start, $stop) = @_;
  my @scores = ();
  my $tabixoutput = `tabix /gscmnt/gc2719/halllab/src/gemini/data/whole_genome_SNVs.tsv.compressed.gz $chrom:$start-$stop`;
  my @tabixoutputlines = split(/\n/,$tabixoutput);
  foreach my $taboutline (@tabixoutputlines) {
    push @scores, split(/,/,(split(/\s+/,$taboutline))[4]);
  }
  
  return max(@scores);
}

sub getfields { # Parse info field of VCF line, getting fields specified in @_. $_[0] must be the info field itself. Returns list of field values
  my $info = shift @_;
  my $pattern = "";
  foreach my $i (0..$#_) { # Build pattern
    $pattern .= (";" . $_[$i] . "=(.*?)");
  }
  $pattern .= "(?:;|\$)";
  my @ans = ($info =~ /$pattern/);
  @ans = ("") x @_ unless @ans;
  if (@ans > 1) {
    return @ans;
  } else {
    return $ans[0];
  }
}
