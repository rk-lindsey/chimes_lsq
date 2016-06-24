#! /usr/bin/perl
# Sort a file according to a selected column in numerical order.
use strict ;
use warnings ;

my @lines ;
my $column = 8 ;

open ( IN, "<$ARGV[0]" ) or die "Could not open $ARGV[0]\n" ;

my $count = 0 ;
while ( <IN> ) {
  $lines[$count++] = $_ ;
}

my @sorted = sort by_column @lines ;

for my $j ( 0 .. $count - 1 ) {
  print $sorted[$j] ;
}

sub by_column
  {
    my @a_col = split(" ", $a) ;
    my @b_col = split(" ", $b) ;

    return $a_col[$column] <=> $b_col[$column] ;
  }

