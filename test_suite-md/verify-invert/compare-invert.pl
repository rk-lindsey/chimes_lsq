#! /usr/bin/perl
use strict ;
use warnings ;

# Compare inverted and non-inverted forces
open(my $file1, "<", $ARGV[0]) or die "Could not open $ARGV[0]" ;
open(my $file2, "<", $ARGV[1]) or die "Could not open $ARGV[1]" ;

my $count = 1 ;

for (1 .. 2) {
  $_ = <$file1> ; 
  $_ = <$file2> ;
}

while ( <$file1> ) {
  my $f1 = $_ ;
  my $f2 = <$file2> ;
  
  $f1 *= -1 ;

  if ( abs($f1 -$f2) > 1.0e-12 ) {
	 die "Bad comparison on line $count\n" ;
  }
  $count++ ;
}


