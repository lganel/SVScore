#!/usr/bin/perl -w

use strict;

my $usage = "reorderheader.pl version 0.0.1
Author: Liron Ganel

Reorders header lines of a VCF file based on the header of a reference file. Header lines not in the reference are placed at the end of the header (before the #CHROM line).

usage: ./reorderheader.pl file.vcf reference.vcf
	file.vcf: VCF whose header is to be reordered
	reference.vcf: VCF with desired order of header lines;
	
	substituting \"stdin\" for either file results in reading that file from standard input\n";


die $usage unless scalar @ARGV == 1 || scalar @ARGV == 2;
my ($file, $reffile) = @ARGV;
die "Standard input can only be used for one input file" if $reffile eq "stdin" && $file eq "stdin";
my ($reffilehandle, $filehandle);
if($reffile eq "stdin") {
  $reffilehandle=*STDIN;
} else {
  open($reffilehandle,"<$reffile") || die "Could not open $reffile: $!";
}

my (%reforder, $line) = ();
my $count = 0;
until(!defined($line = <$reffilehandle>) || $line =~ /^#CHROM/) {
  $reforder{$line} = $count;
  $count++;
}
close $reffilehandle;
die "Malformed VCF" unless defined $line;

if($file eq "stdin") {
  $filehandle=*STDIN;
} else {
  open($filehandle,"<$file") || die "Could not open $file: $!";
}
my @filelines = <$filehandle>;
my @filenoheader = grep {!/^#/} @filelines;
my @sortedhead = sort mysort (grep {/^#/} @filelines);

push @sortedhead, @filenoheader;
print foreach (@sortedhead);

sub mysort {
  if (defined($reforder{$a}) && defined($reforder{$b})) {
    return $reforder{$a} <=> $reforder{$b};
  } elsif ($a =~ /^#CHROM/) {
    return 1;
  } elsif ($b =~ /^#CHROM/) {
    return -1;
  } elsif (defined $reforder{$a}) {
    return -1;
  } elsif (defined $reforder{$b}) {
    return 1;
  } else {
    return $a cmp $b;
  }
}
