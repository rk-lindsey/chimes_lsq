#! /usr/bin/perl
use strict ;
use warnings ;

# Invert a force output file.
open(my $file1, "<", $ARGV[0]) or die "Could not open $ARGV[0]" ;

my $count = 1 ;

$_ = <$file1> ; 

while ( <$file1> ) {
  if ( /[0-9]+/ ) {
	 my $f1 = $_ ;
  
	 $f1 *= -1 ;
	 printf("%13.6e\n",$f1) ;
  }
}

print "\n\n" ;
