#! /usr/bin/perl
#
# Given parameters for the overcoordination term, fit other parameters using least squares.
# overcoord.template is a template file.
use strict ;
use warnings ;

# Set this to the number of MD frames to use in nonlinear least squares.
my $nframes = 49 ;

# Set this to the number of overcoordination parameters to fit nonlinearly.
my $nparams = 4 ;

open(TEMPLATE, "<overcoord.template") or die "Could not open overcoord.template\n" ;
open(PARAMS, "<overcoord.in") or die "Could not open overcoord.in\n" ;

my @params ;
my @sign = (1, -1, 1, -1) ;

foreach my $i ( 0 .. $nparams - 1 ) {
  $params[$i] = <PARAMS> ;
}

open(OUT, ">splines_over.in") or die "Could not write to splines_over.in\n" ;
while ( <TEMPLATE> ) {
  print OUT ;
}
print OUT "over_params " ;
printf OUT ("%d ", $nparams+1) ;

# The linear overcoordination parameter is given a fixed value of 50.0
# here.  The actual value is set by the linear least squares code.
print OUT "50.0 " ;
foreach my $i ( 0 .. $nparams - 1 ) {
  printf OUT ("%21.16e ", $params[$i]) ;
}
print OUT "\n" ;
open(ERROR, ">error.out") or die "Could not write to error.out\n" ;

# Don't allow parameters to change sign.
foreach my $i ( 0 .. $nparams - 1 ) {
  if ( $params[$i] * $sign[$i] < 0.0 ) {
    print ERROR "1.0e+06\n" ;
    exit(0) ;
  }
}

if ( system("splines_ls splines_over.in $nframes > splines_ls.out") ) {
  print ERROR "1.0e+06\n" ;
  exit(0) ;
}
if ( system("lsq.py A.txt b.txt params.header > params.over.txt") ) {
  print ERROR "1.0e+06\n" ;
  exit(0) ;
} else {
  open(PARAMS,"<params.over.txt") or die "Could not open params.over.txt\n" ;
  while ( <PARAMS> ) {
    if ( /RMS force error/ ) {
      my @f = split(" ") ;
      my $rmserr = $f[5] ;
      print ERROR $rmserr, "\n" ;
    }
  }
}

  


