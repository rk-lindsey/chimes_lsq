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

my $lambda = 1.25 ;

while ( <> ) {
  # Read initial comments.
  if ( !( $_ =~ /^\#/ ) )  {
    last ;
  }
  print
  ;
}
my $rmin = $_ ;
chomp $rmin ;

my $rmax = <> ;
chomp $rmax ;

my $order = <> ;
chomp $order ;

print "Rmin = $rmin\n" ;
print "Rmax = $rmax\n" ;
print "Order = $order\n" ;

# Read past overcoordination parameters.
my $nover = <> ;
chomp $nover ;
for my $i ( 1 .. $nover ) {
  $_ = <> ;
}

foreach my $atom (0 .. 2) {
  open(OUT, ">cheby.$atom.out\n") ;
  open(POT, ">chebypot.$atom.out\n") ;
  foreach my $k (0 .. $order-1) {
    $_ = <> ;
    die "File ended prematurely\n" unless defined($_) ;
    my @f = split(" ") ;
    $coef[$k] = $f[1] ;
#    printf("COEF $k $coef[$k]\n" ) ;
  }
  
  for ( my $i = 0 ; $i < 100 ; $i++ ) {
    my $r = $rmin + $i * ($rmax-$rmin)/99 ;

    my ($ravg, $rdiff, $x) ;

    if ( $inverse_r > 0 ) {
      $ravg = 0.5 * ( 1.0 / $rmin + 1.0 / $rmax ) ;
      $rdiff = 0.5 * ( 1.0 / $rmin - 1.0 / $rmax ) ;
      $x = (1.0/$r - $ravg) / $rdiff ;
    } elsif ( $morse > 0 ) {
      my $xmax = exp(-$rmin/$lambda) ;
      my $xmin = exp(-$rmax/$lambda) ;
      $ravg = 0.5 * ( $xmax + $xmin ) ;
      $rdiff = 0.5 * ( $xmax - $xmin ) ;
      $x = ( exp(-$r/$lambda) - $ravg )/ $rdiff ;
    } else {
      $ravg = 0.5 * ( $rmin + $rmax ) ;
      $rdiff = 0.5 * ( $rmax - $rmin ) ;
      $x = ($r - $ravg) / $rdiff ;
    }
    if ( $x < -1.0000001 || $x > 1.000001 ) {
      die "Error: x = $x\n" ;
    }

    $tn[0] = 1.0 ;
    $tn[1] = $x ;
    $tnd[0] = 1.0 ;
    $tnd[1] = 2.0 * $x ;
    foreach my $k (2 .. $order) {
      $tn[$k] = 2.0 * $x * $tn[$k-1] - $tn[$k-2] ;
      $tnd[$k] = 2.0 * $x * $tnd[$k-1] - $tnd[$k-2] ;
    }
    for ( my $k = $order ; $k >= 1 ; $k-- ) {
      $tnd[$k] = $k * $tnd[$k-1] ;
    }
    $tnd[0] = 0 ;
#    for ( my $k= 0 ; $k <= $order ; $k++ ) {
#     print ("x = $x TND $k = $tnd[$k], TN[$k] = $tn[$k]\n") ;
#    }

    $force = 0.0 ;
    my $pot   = 0.0 ;
    if ( $morse > 0 ) {
      my $fcut0 = (1.0 - $r/$rmax) ;
      my $fcut = $fcut0 * $fcut0 * $fcut0 ;
      my $fcutprime = -3.0 * $fcut0 * $fcut0 / $rmax ;
      foreach my $k (0 .. $order-1) {

	die "Coef[$k] not defined\n" unless defined($coef[$k]) ;

	$force += $fcut * $coef[$k] * $tnd[$k+1] 
	  * (-exp(-$r/$lambda) / $lambda) / $rdiff 
	  + $fcutprime * $coef[$k] * $tn[$k+1] ;
	$pot += $fcut * $coef[$k] * $tn[$k+1] ;
      }
    } else {
      foreach my $k (0 .. $order - 1) {
	$force += ($rmax-$r) * $coef[$k] * $tn[$k] / ($r*$r) ;
      }
    }
    printf OUT ("%8.5f %8.5f\n", $r, $force) ;
    printf POT ("%8.5f %8.5f\n", $r, $pot) ;

  }
}

