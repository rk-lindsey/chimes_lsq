#! /usr/bin/perl
#
# Translate all coordinates in an xyz file by a given amount.
#
use strict ;
use warnings ;

my $dx = 10.0 ;
my $dy = 11.0 ;
my $dz = 120.0 ;

my $natoms = <> ;
print $natoms ;

$_ = <> ;
print $_ ;

for my $i ( 1 .. $natoms ) {
  $_ = <> ;
  my ( $ele, $x, $y, $z) = split(" ") ;
  $x += $dx ;
  $y += $dy ;
  $z += $dz ;
  printf("%s %21.14e %21.14e %21.14e\n", $ele, $x, $y, $z) ;
}

