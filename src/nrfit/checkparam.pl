#! /usr/bin/perl
#
# Check optimization output.

use strict ;
use warnings ;

my $s = <> ;
my ($x) = split(" ", $s) ;

$s = <> ;
my ($y) = split(" ", $s) ;

if ( abs($x) > 1.0e-02 ) {
    print "Error: X is non-zero\n" ;
    exit(1) ;
} elsif ( abs($y) > 1.0e-02 ) {
    print "Error: Y is non-zero\n" ;
    exit(1) ;
}
print "Output is OK\n" ;
exit(0) ;

