#! /usr/bin/perl -w
#
# File: filterLowMappability.pl
# Created: Tue Oct 15 15:01:08 2013
# $Id: $
#
# Copyright (C) 2013 by Per Unneberg
#
# Author: Per Unneberg
#
# Description:
#
# Read a gem mappability file and print sites below a given
# k-frequency (as defined in the GEM paper)

use strict;

if ($#ARGV!=0) {
  die "Usage: $0 [-k max k-frequency] mappabilityfile\n";
}
my $kfreq = 10;
my $fh;
my $seenchr=0;
my $scaffold = undef;
my $cref = {};
my $pos = 0;
open($fh, $ARGV[0]) || die "Can't open $ARGV[0]: $!\n";
open(OUT, "> sizes.txt") || die "Can't open sizes.txt: $!\n";
while (<$fh>) {
  chomp;
  if (/^~(scaffold_\d+)/ || /^~(chrM)/) {
    $scaffold = $1;
    if ($pos > 0) {
      print OUT $pos, "\n";
    }
    print  OUT $scaffold, "\t";
    $pos = 0;
    $seenchr=1;
  } else {
    if ($seenchr==0) {
      if (/^'(.*)'~\[.*-(\d+)\]/) {
	$cref->{$1} = $2;
      } else {
	next;	
      }
    } else {
      my @chars=split //, $_;
      foreach my $c (@chars) {
	$pos++;
	if ($c =~ /\x00/) {
	  $c = " ";
	  print STDERR "WARNING: $scaffold $pos: Converting NULL character to space\n";
	}
	if (!defined $c) {
	  print STDERR "WARNING: $scaffold $pos: $c not defined\n";
	  $c = " ";
	}
	if (!exists $cref->{$c}) {
	  print STDERR "WARNING: $scaffold $pos: '$c' doesn't exist in dictionary\n";
	  $c = " ";
	}
	if ($cref->{$c} < $kfreq && $cref->{$c} > 0) {
	  print $scaffold, "\t", $pos, "\n";
	}
      }
    }
  }
}
if ($pos > 0) {
  print OUT $pos, "\n";
}
close(OUT);
close($fh);
