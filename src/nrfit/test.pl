#!/usr/bin/perl
open(OUT, ">error.out") ;
open(IN, "<params.in") ;
$_ = <IN> ;
$x = $_ ;
$_ = <IN> ;
$y = $_ ;
printf OUT ("%21.13e\n", $x * $x - 1.0 + 2.0 * $y * $y ) ;


