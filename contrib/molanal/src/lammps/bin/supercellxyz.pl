#! /usr/bin/perl
## Make a supercell from a lammps data file.  
## This program does not produce bond connectivity information 
## e.g. this will work with ReaxFF, LJ, or metals, but not DREIDING, etc
##     die "Usage:  supercelldata.pl nx ny nz data_file > output_file 
use strict ;
use warnings ;

my ($natom, $na, $nb, $nc, $copy,$ax, $bx, $cx) ;
my ($xlo, $xhi, $zlo, $zhi, $ylo, $yhi, $j,$count) ;

if ( scalar(@ARGV) >= 3 ) {
    $na = $ARGV[0] ;
    shift ;
    $nb = $ARGV[0] ;
    shift ;
    $nc = $ARGV[0] ;
    shift ;
} else {
    die "Usage:  supercelldata.pl nx ny nz data_file > output_file \n" ;
}

$copy = $na * $nb * $nc ;
$_ = <> ;
$_ =~ s/\#//g ;
print "#Supercell for $_\n" ;
$_ = <> ;
$_ = <> ;
($natom) = split(" ") ;
$natom *= $copy ;
print "$natom atoms\n\n" ;

$_ = <> ;
$_ = <> ;
($xlo, $xhi) = split(" ") ;
$_ = <> ;
$ax = $xhi - $xlo ;
$ax *= $na ;

($ylo, $yhi) = split(" ") ;
$bx = $yhi - $ylo ;
$bx *= $nb ;

$_ = <> ;
($zlo, $zhi) = split(" ") ;
$cx = $zhi - $zlo ;
$cx *= $nc ;

printf( "%11.4f %11.4f xlo xhi\n", $xlo, $xlo + $ax) ;
printf( "%11.4f %11.4f ylo yhi\n", $ylo, $ylo + $bx) ;
printf( "%11.4f %11.4f zlo zhi\n", $zlo, $zlo + $cx) ;

$_ = <> ;
foreach ( 0 .. 8 ) {
    $_ = <> ;
    print  $_ ;
}

$count = 0 ;
while ( <> ) {
    my ($k, $l, $index, $type, $xx, $yy, $zz,
	$x, $y, $z, $q) ;
    ($index, $type, $q, $x, $y, $z) = split(" ") ;
    for ( $j = 0 ; $j < $na ; $j++ ) {
	$xx = $x + $j * $ax ;
	for ( $k = 0 ; $k < $nb ; $k++ ) {
	    $yy = $y + $k * $bx ;
	    for ( $l = 0 ; $l < $nc ; $l++ ) {
		$zz = $z + $l * $cx ;
		$count++ ;
		printf  ("%d %d %11.4f %11.4f %11.4f %11.4f\n", $count,
			 $type, $q, $xx, $yy, $zz) ;
	    }
	}
    }
}

