#!/usr/bin/perl -w

use strict;
use Getopt::Std;

# split.pl takes a VCF file and splits it into many files, each with approximately the given number of variant lines, ensuring that multiline variants all end up in the same file. The input VCF files can be bgzipped or uncompressed.

my $linenum = 1;
my $filenum = 1;

sub main::HELP_MESSAGE(); # Declaration to make perl happy
$Getopt::Std::STANDARD_HELP_VERSION = 1; # Make --help and --version flags halt execution
$main::VERSION = '0.1';

my %options = ();
getopts('n:o:',\%options);
my $linesperfile = (defined $options{'n'} ? $options{'n'} : 10000);
my $outputdir = (defined $options{'o'} ? $options{'o'} : ".");
mkdir $outputdir unless -d $outputdir;

die "Please provide a valid input file" unless defined $ARGV[0] && -s $ARGV[0];
my $infile = $ARGV[0];
my $catprefix = ($infile =~ /\.gz$/ ? "z" : "");
print STDERR "Grabbing header\n";
system("${catprefix}cat $infile | grep \"^#\" > $outputdir/header.tmp") && die "Could not extract header from $infile: $!"; # Extract header
print STDERR "Sorting file\n";
open(IN,"${catprefix}cat $infile | grep -v \"^#\" | sort -k3,3n |") || die "Could not open $infile: $!"; # Remove header, sort, and open

open(OUT,"| sort -k1,1V -k2,2n | cat tmp/header.tmp - | bgzip -c > $outputdir/split1.vcf.gz") || die "Cannot pipe to bgzip: $!"; 

my $lastid;
print STDERR "Starting file 1\n";

while(my $line = <IN>) {
  my $id = (split(/\s+/,$line))[2];
  if ($linenum >= $linesperfile) { # Need to split
    if ($id =~ /_/) { # Last line is multiline variant
      my ($event) = ($id =~ /^(.*)_/);
      my ($lastevent) = ($lastid =~ /^(.*)_/);
      while($lastevent ne $event) { # Keep printing until a new event is found
	print OUT $line;
	$lastevent = $event;
	$line = <IN>;
	$id = (split(/\s+/,$line))[2];
	if ($id =~ /_/) { # Get next line's event ID from its variant ID
	  ($event) = ($id =~ /^(.*)_/);
	} else {
	  $event = $id;
	}
      }
    }

    # Close old file and open a new one
    print OUT $line;
    close OUT;
    $filenum++;
    $linenum = 1;
    #my $filenumstring = ($filenum < 10 ? "0$filenum" : $filenum);
    open(OUT,"| sort -k1,1V -k2,2n | cat tmp/header.tmp - | bgzip -c > $outputdir/split$filenum.vcf.gz") || die "Cannot pipe to bgzip: $!"; 
    print STDERR "Starting file $filenum\n";
	
  } else { # Don't need to split
    print OUT $line;
    $linenum++;
  }
  $lastid = $id;
}
close OUT;

unlink("$outputdir/header.tmp") || die "Could not delete header.tmp: $!";
