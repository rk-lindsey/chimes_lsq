#! /usr/bin/perl
#  Print out potential values from a spline file.
#  Usage: spline_vals.pl < params.txt
#
use strict ;
use warnings ;
use POSIX qw/floor/ ;

my $rmin = <> ;
chomp $rmin ;

my $rmax = <> ;
chomp $rmax ;

my $sdelta = <> ;
chomp $sdelta ;

my $snum = floor(($rmax - $rmin) / $sdelta) ;

print "Number of parameters = $snum\n" ;

my (@dvdr, @pot) ;

for ( my $i = 0 ; $i < 3 ; $i++ ) {
    open(OUT, ">spline.$i.out\n") ;
    open(POT, ">splinepot.$i.out\n") ;
    foreach my $k  ( 0 .. $snum - 1 ) {
      $_ = <> ;
      my ( $idx, $val ) = split(" ") ;
      my $kidx = 2 * $k + 2 * $i * $snum ;
      if ( $idx != $kidx ) {
	die "Index error: $idx $kidx\n" ;
      }
      $_ = <> ;
      my $r = $rmin + $k * $sdelta ;
      $dvdr[$k] = $val ;
      printf OUT ("%13.6e %13.6e\n", $r, $val) ;
    }
    $pot[$snum] = 0.0 ;
    for ( my $k = $snum - 1 ; $k >= 0 ; $k-- ) {
      $pot[$k] = $pot[$k+1] - $dvdr[$k] * $sdelta ;
      my $r = $rmin + $k * $sdelta ;
      printf POT ("%13.6e %13.6e\n", $r, $pot[$k]) ;
    }
  }
