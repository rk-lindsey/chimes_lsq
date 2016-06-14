#! /usr/bin/perl
# Print out potential values from a chebyshev params file.
#
# Usage: cheby_print.pl < params.cheby.txt
#
use strict ;
use warnings ;

my $inverse_r = 0 ;
my $morse     = 1 ;

my (@coef, @tn, @tnd) ;
my ($x, $force) ;

my $npair ;
my (@rmin, @rmax, @order, @lambda) ;

while ( <> ) {
  # Read to definition of pair type.
  if ( $_ =~ /npair/ ) {
    my @f = split(" ") ;
    $npair = $f[1] ;
    chomp $npair ;
    if ( $npair <= 0 ) {
      die "A positive number of pairs is required\n" ;
    }
    last ;
  }
}
$_ = <> ;
$_ = <> ;
if ( /pair_type chebyshev/ ) {
  print "Chebyshev pair type found\n" ;
} else {
  die "Error: chebyshev pair type was not found\n" ;
}

for my $i ( 0 .. $npair - 1 ) {
  $_ = <> ;
  ($rmin[$i], $rmax[$i], $lambda[$i],$order[$i]) = split(" ") ;
  print "Rmin = $rmin[$i] Rmax = $rmax[$i] Lambda = $lambda[$i] Order = $order[$i]\n" ;
}

while ( <> ) {
  if ( $_ =~ /least squares parameters/ )  {
    last ;
  }
  print ;
}

foreach my $ipr (0 .. $npair - 1) {
  open(OUT, ">cheby.$ipr.out\n") ;
  open(POT, ">chebypot.$ipr.out\n") ;
  foreach my $k (0 .. $order[$ipr]-1) {
    $_ = <> ;
    die "File ended prematurely\n" unless defined($_) ;
    my @f = split(" ") ;
    $coef[$k] = $f[1] ;
    printf("COEF $k $coef[$k]\n" ) ;
  }
  
  for ( my $i = 0 ; $i < 100 ; $i++ ) {
    my $r = $rmin[$ipr] + $i * ($rmax[$ipr]-$rmin[$ipr])/99 ;

    my ($ravg, $rdiff, $x) ;

    if ( $inverse_r > 0 ) {
      $ravg = 0.5 * ( 1.0 / $rmin[$ipr] + 1.0 / $rmax[$ipr] ) ;
      $rdiff = 0.5 * ( 1.0 / $rmin[$ipr] - 1.0 / $rmax[$ipr] ) ;
      $x = (1.0/$r - $ravg) / $rdiff ;
    } elsif ( $morse > 0 ) {
      my $xmax = exp(-$rmin[$ipr]/$lambda[$ipr]) ;
      my $xmin = exp(-$rmax[$ipr]/$lambda[$ipr]) ;
      $ravg = 0.5 * ( $xmax + $xmin ) ;
      $rdiff = 0.5 * ( $xmax - $xmin ) ;
      $x = ( exp(-$r/$lambda[$ipr]) - $ravg )/ $rdiff ;
    } else {
      $ravg = 0.5 * ( $rmin[$ipr] + $rmax[$ipr] ) ;
      $rdiff = 0.5 * ( $rmax[$ipr] - $rmin[$ipr] ) ;
      $x = ($r - $ravg) / $rdiff ;
    }
    if ( $x < -1.0000001 || $x > 1.000001 ) {
      die "Error: x = $x\n" ;
    }

    $tn[0] = 1.0 ;
    $tn[1] = $x ;
    $tnd[0] = 1.0 ;
    $tnd[1] = 2.0 * $x ;
    foreach my $k (2 .. $order[$ipr]) {
      $tn[$k] = 2.0 * $x * $tn[$k-1] - $tn[$k-2] ;
      $tnd[$k] = 2.0 * $x * $tnd[$k-1] - $tnd[$k-2] ;
    }
    for ( my $k = $order[$ipr] ; $k >= 1 ; $k-- ) {
      $tnd[$k] = $k * $tnd[$k-1] ;
    }
    $tnd[0] = 0 ;
#    for ( my $k= 0 ; $k <= $order ; $k++ ) {
#     print ("x = $x TND $k = $tnd[$k], TN[$k] = $tn[$k]\n") ;
#    }

    $force = 0.0 ;
    my $pot   = 0.0 ;
    if ( $morse > 0 ) {
      my $fcut0 = (1.0 - $r/$rmax[$ipr]) ;
      my $fcut = $fcut0 * $fcut0 * $fcut0 ;
      my $fcutprime = -3.0 * $fcut0 * $fcut0 / $rmax[$ipr] ;
      foreach my $k (0 .. $order[$ipr]-1) {

	die "Coef[$k] not defined\n" unless defined($coef[$k]) ;

	$force += $fcut * $coef[$k] * $tnd[$k+1] 
	  * (-exp(-$r/$lambda[$ipr]) / $lambda[$ipr]) / $rdiff 
	  + $fcutprime * $coef[$k] * $tn[$k+1] ;
	$pot += $fcut * $coef[$k] * $tn[$k+1] ;
      }
    } else {
      foreach my $k (0 .. $order[$ipr] - 1) {
	$force += ($rmax[$ipr]-$r) * $coef[$k] * $tn[$k] / ($r*$r) ;
      }
    }
    printf OUT ("%8.5f %8.5f\n", $r, $force) ;
    printf POT ("%8.5f %8.5f\n", $r, $pot) ;

  }
}

