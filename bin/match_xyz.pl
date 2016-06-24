#! /usr/bin/perl
use strict ;
use warnings ;

# Script to match a permuted (xyz_file2) and non-permuted xyz file (xyz_file1), 
# then print force_file2 in the atom order of the non-permuted xyz file.

die "Usage: match_xyz.pl xyz_file1 xyz_file2 force_file2\n" unless ( scalar(@ARGV) == 3 ) ;

open(IN1, "<$ARGV[0]") or die "Could not open $ARGV[0]\n" ;
open(IN2, "<$ARGV[1]") or die "Could not open $ARGV[1]\n" ;
open(FORCE, "<$ARGV[2]") or die "Could not open $ARGV[2]\n" ;

my $count = 0 ;
my @line_array ;
while ( <IN1> ) {
  $line_array[$count] = $_ ;
  $count++ ;
}

$count = 0 ;
my %line_hash2 ;
while ( <IN2> ) {
  $line_hash2{$_} = $count ;
  $count++ ;
}

foreach my $line ( @line_array ) {
  die "$line is not found in file 2\n" unless defined($line_hash2{$line}) ;
  $count = $line_hash2{$line} ;
#  printf("%03d: ",$count) ;
#  print $line ;
}


my @force_lines ;
while ( <FORCE> ) {
  push(@force_lines,$_ );
}

foreach my $j ( 2 .. scalar(@line_array) - 1 ) {
  my $line = $line_array[$j] ;
  $count = $line_hash2{$line} - 2 ;
  print "$count: $force_lines[3*$count]" ;
  print "$count: $force_lines[3*$count+1]" ;
  print "$count: $force_lines[3*$count+2]" ;
}
