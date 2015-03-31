#! /usr/bin/perl
#
# Given parameters for the overcoordination term, fit other parameters using least squares.
#
use strict ;
use warnings ;

my $nframes = 49 ;

open(TEMPLATE, "<overcoord.template") or die "Could not open overcoord.template\n" ;
open(PARAMS, "<overcoord.in") or die "Could not open overcoord.in\n" ;

my $nparams = <PARAMS> ;
my @params ;

foreach my $i ( 0 .. $nparams - 1 ) {
  $params[$i] = <PARAMS> ;
}

open(OUT, ">splines_over.in") or die "Could not write to splines_over.in\n" ;
while ( <TEMPLATE> ) {
  print OUT ;
}
print OUT "over_params " ;
printf OUT ("%d ", $nparams) ;
foreach my $i ( 0 .. $nparams - 1 ) {
  printf OUT ("%21.16e ", $params[$i]) ;
}
print OUT "\n" ;

if ( system("splines_ls splines_over.in $nframes") ) {
  printf("Error: splines_ls failed\n") ;
}
if ( system("lsq.py A.txt b.txt > params.over.txt") ) {
  printf("Error: lsq.py failed\n") ;
}


