#! /usr/bin/perl
# Checks maximum velocities for a LAMMPS vx vy vz dump file.
# Usage:  vel_check.pl < dump_file

use strict ;
my ($q, $timestep,$natom,$indx,$type,$indx_max,$j, $vel_max,
    $vel_lim, $vel_all) ;

$vel_lim = 10.0 ;
$vel_all = 0.0 ;
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
    if ( /ITEM: ATOMS id type vx vy vz/ ) {
	$vel_max = 0 ;
	for ( $j = 0 ; $j < $natom ; $j++ ) {
	    my ($fx, $fy, $fz) ;
	    $_ = <> ;	
	    ($indx, $type, $fx, $fy, $fz) = split(" ") ;
	    $q = $fx * $fx + $fy * $fy + $fz * $fz ;
	    $q = sqrt($q) ;
	    if ( $q > $vel_max ) {
		$vel_max = abs($q) ;
		$indx_max = $indx ;
	    }
	}
#	print OUT ("%7d %7d %11.4f\n", $timestep, $indx_max, $vel_max) ;
	if ( $vel_max > $vel_lim ) {
	    die "ERROR: SPEED WAS TOO LARGE: $vel_max\n" ;
	} elsif ( $vel_max > $vel_all ) {
	    $vel_all = $vel_max ;
	}
    }
}
printf("Maximum speed found = %11.4f (OK).\n", $vel_all) ;


