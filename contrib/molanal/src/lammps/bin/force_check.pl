#! /usr/bin/perl
# Checks maximum forces for a LAMMPS fx fy fz dump file.
# Usage:  force_check.pl < dump_file

use strict ;
my ($q, $timestep,$natom,$indx,$type,$indx_max,$j, $force_max,
    $force_lim, $force_all) ;

$force_lim = 1.0e5 ;
$force_all = 0.0 ;
#printf("Step  Atom  Max |Q|\n") ;

while ( <> ) {
    if ( /ITEM: TIMESTEP/ ) {
	$timestep = <> ;
	chop $timestep ;
#	print "Checking timestep $timestep\n" ;
    }
    if ( /ITEM: NUMBER OF ATOMS/ ) {
	$natom = <> ;
#	print "Number of atoms = $natom\n" ;
    }
    if ( /ITEM: ATOMS id type fx fy fz/ ) {
	$force_max = 0 ;
	for ( $j = 0 ; $j < $natom ; $j++ ) {
	    my ($fx, $fy, $fz) ;
	    $_ = <> ;	
	    ($indx, $type, $fx, $fy, $fz) = split(" ") ;
	    $q = $fx * $fx + $fy * $fy + $fz * $fz ;
	    $q = sqrt($q) ;
	    if ( $q > $force_max ) {
		$force_max = abs($q) ;
		$indx_max = $indx ;
	    }
	}
#	print OUT ("%7d %7d %11.4f\n", $timestep, $indx_max, $force_max) ;
	if ( $force_max > $force_lim ) {
	    die "ERROR: FORCE WAS TOO LARGE: $force_max\n" ;
	} elsif ( $force_max > $force_all ) {
	    $force_all = $force_max ;
	}
    }
}
printf("Maximum force found = %11.4f (OK).\n", $force_all) ;


