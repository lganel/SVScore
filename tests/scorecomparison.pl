#!/usr/bin/perl -w

use strict;

my ($testname,$truthfile,$testfile) = @ARGV;

die "********$testname FAILED: $truthfile EMPTY OR NONEXISTENT - DID SVSCORE DOWNLOAD SUCCEED?********\n" unless -s $truthfile;
die "********$testname FAILED: $testfile EMPTY OR NONEXISTENT - IS SVSCORE PROVIDING OUTPUT?********\n" unless -s $testfile;

open(TRUTH,"< $truthfile") || die "Can't open $truthfile: $!";
open(TEST,"< $testfile") || die "Can't open $testfile: $!";

my $proximity = 0.001; # For two scores to be considered the same, the min score (between truth and test) must be within 0.1% of the max score
my $success = 1;
while (my $truthline = <TRUTH>) {
  my $testline = <TEST>;
  next if $truthline =~ /^#/;
  my @truthsplit = split(/\s+/,$truthline);
  my ($truthid, $truthinfo) = @truthsplit[2,7];
  my @truthscores = grep {/^SVSCORE/} split(/;/,$truthinfo); # Get SVScore annotations from truth file
  my @testsplit = split(/\s+/,$testline);
  my ($testid, $testinfo) = @testsplit[2,7];
  my @testscores = grep {/^SVSCORE/} split(/;/,$testinfo); # Get SVScore annotations from test file
  printmessage("$testname FAILED: IDS $truthid AND $testid DON'T MATCH AT LINE $. ($truthid in $truthfile and $testid in $testfile)",1) unless $truthid eq $testid;
  printmessage("$testname FAILED: ID $truthid LINE HAS DIFFERING NUMBERS OF FIELDS (",scalar @truthsplit," in $truthfile and ",scalar @testsplit," in $testfile)",1) unless @truthsplit == @testsplit;
  foreach my $i (0..6,8..$#truthsplit) {
    printmessage("$testname FAILED: ID $truthid LINE DIFFERS AT COLUMN ",$i+1," ($truthsplit[$i] in $truthfile and $testsplit[$i] in $testfile)",1) unless $truthsplit[$i] eq $testsplit[$i];
  }
  my %truthscores = createscorehash(\@truthscores); # Place scores in hash
  my %alltruthscores = %truthscores; # Copy for later use
  my %testscores = createscorehash(\@testscores);
  my %alltestscores = %testscores;
  my (@missingintruth, @missingintest, @different);
  foreach my $opinterval (keys %truthscores) {
    unless(defined $testscores{$opinterval}) { # Identify scores in truth file but missing in test file
      push @missingintest, $opinterval;
    } else {
      my ($max,$min) = ($truthscores{$opinterval} >= $testscores{$opinterval} ? ($truthscores{$opinterval}, $testscores{$opinterval}) : ($testscores{$opinterval}, $truthscores{$opinterval}));
      push @different, $opinterval unless $min == $max || ($min >= (1-$proximity) * $max); # Identify scores with more than $proximity percent difference between truth and test
      delete $testscores{$opinterval}; # Remove score from test hash
    }
    delete $truthscores{$opinterval}; # Remove score from truth hash
  }
  @missingintruth = keys %testscores; # Everything left in the test hash must be missing from the truth hash
  if (@missingintruth || @missingintest || @different) { # Print error messages on individual lines
    printmessage("$testname FAILED!") if $success; # If this is the first error, print failure message
    $success = 0;
    if (@missingintruth) {
      print "ID $truthid: ",join(', ',@missingintruth)," missing in $truthfile\n";
    }
    if (@missingintest) {
      print "ID $truthid: ",join(', ',@missingintest)," missing in $testfile\n";
    }
    if (@different) {
      foreach my $opinterval (@different) {
	print "ID $truthid: $opinterval = $alltruthscores{$opinterval} in $truthfile and $alltestscores{$opinterval} in $testfile\n";
      }
    }
    print "\n";
  }
}

printmessage("$testname PASSED") if $success;

sub createscorehash {
  my %scorehash = ();
  foreach my $scorestr (@{$_[0]}) {
    my @pair = split(/=/,$scorestr);
    die "Error in $scorestr" if @pair > 2;
    $scorehash{$pair[0]} = $pair[1];
  }
  return %scorehash;
}

sub printmessage {
  my ($message, $die) = @_;
  my $length = length($message);
  print "\n","*" x ($length+4),"\n*"," " x ($length+2),"*\n* $message *\n*"," " x ($length+2),"*\n", "*" x ($length+4),"\n";
  exit if $die;
}
