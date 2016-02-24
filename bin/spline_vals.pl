#! /usr/bin/perl
#  Print out potential values from a spline file.
#  Usage: spline_vals.pl < params.txt
#
use strict ;
use warnings ;
use POSIX qw/floor/ ;

while ( <> ) {
  # Read initial comments.
  if ( !( $_ =~ /^\#/ ) )  {
    last ;
  }
  print
  ;
}

my ($token, $npair) = split(" ") ;

die "Error: npair expected on first non-comment line\n" 
  unless ($token eq "npair") ;

chomp $npair ;

if ( $npair <= 0 ) {
  die "A positive number of pairs is required\n" ;
}

my (@rmin, @rmax, @sdelta, @snum) ;
for my $i ( 0 .. $npair - 1 ) {
  $_ = <> ;
  ($rmin[$i], $rmax[$i], $sdelta[$i]) = split(" ") ;
  $snum[$i] = 2 + floor(($rmax[$i] - $rmin[$i]) / $sdelta[$i]) ;
  print "Rmin = $rmin[$i] Rmax = $rmax[$i] Delta = $sdelta[$i] N = $snum[$i]\n" ;
}

# Read past overcoordination parameters.
my $nover = <> ;
chomp $nover ;
print "Number of overcoordination parameters = $nover\n" ;
for my $i ( 1 .. $nover ) {
  $_ = <> ;
  print "Overcoordination parameter $i = $_" ;
}

my (@dvdr, @pot) ;

my $offset = 0 ;
for ( my $i = 0 ; $i < $npair ; $i++ ) {
    open(OUT, ">spline.$i.out\n") ;
    open(POT, ">splinepot.$i.out\n") ;
    foreach my $k  ( 0 .. $snum[$i] - 1 ) {
      $_ = <> || die "Could not read parameter\n" ;
      # Read forces from every other line.
      my ( $idx, $val ) = split(" ") ;
      my $kidx = 2 * $k + $offset ;
      # print "IDX $idx KIDX $kidx VAL $val\n" ;
      if ( $idx != $kidx ) {
	die "Index error: $idx $kidx\n" ;
      }
      $_ = <> || die "Could not read parameter\n" ;
      my $r = $rmin[$i] + $k * $sdelta[$i] ;
      $dvdr[$k] = $val ;
      printf OUT ("%13.6e %13.6e\n", $r, $val) ;
    }
    $offset += 2 * $snum[$i] ;
    $pot[$snum[$i]] = 0.0 ;
    for ( my $k = $snum[$i] - 1 ; $k >= 0 ; $k-- ) {
      $pot[$k] = $pot[$k+1] - $dvdr[$k] * $sdelta[$i] ;
      my $r = $rmin[$i] + $k * $sdelta[$i] ;
      printf POT ("%13.6e %13.6e\n", $r, $pot[$k]) ;
    }
  }
