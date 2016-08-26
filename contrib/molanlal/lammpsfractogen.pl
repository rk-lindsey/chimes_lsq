#! /usr/bin/perl
# This program translates a LAMMPS fractional coordinate dump file 
# into a dftb gen file.
# The gen file can be read by molanal.new or gr_calc.
# Usage:  lammpsfractogen.pl <TRAJECTORY> <lammps.dat> <trajectory.gen>
# <lammps.dat> is a file containing the name of each element in order according to the
# lammps dump file in the first line.
#
# Larry Fried, 5/9/07
#
use strict ;

{
    my $natoms ;
    my ($xlow, $xhigh, $ylow, $yhigh, $zlow, $zhigh) ;
    my ($xbox, $ybox, $zbox) ;
    my ($j, $k, $type) ;
    my ($x, $y, $z) ;
    my $traj_file ;
    my $data_file ;
    my $gen_file ;
    my @element ;
    my $frame ;

    if ( $#ARGV != 2 ) {
	die 'Usage:  lammpstogen.pl <TRAJECTORY> <lammpsdat> <trajectory.gen>' ;
    }
    $traj_file = $ARGV[0] ;
    $data_file = $ARGV[1] ;
    $gen_file  = $ARGV[2] ;

    open DATA, "< $data_file" || die "Could not open $data_file\n" ;
    open(TRAJ, "< $traj_file") || die "Could not open $traj_file\n" ;
    open(GEN, "> $gen_file") || die "Could not open $gen_file\n" ;

    while ( <DATA> ) {
	($element[$j]) = split(" ") ;
	$j++ ;
    }

    while ( <TRAJ> ) {
	if ( /ITEM: NUMBER OF ATOMS/ ) {
	    $_ = <TRAJ> ;
	    $natoms = $_ ;
	    chop $natoms ;
	    print "Reading Frame = $frame\n" ;
	    $frame++ ;
	}
	if ( /ITEM: BOX BOUNDS/ ) {
	    $_ = <TRAJ> ;
	    ($xlow, $xhigh) = split(" ") ;
	    $_ = <TRAJ> ;
	    ($ylow, $yhigh) = split(" ") ;
	    $_ = <TRAJ> ;
	    ($zlow, $zhigh) = split(" ") ;
	    $xbox = $xhigh - $xlow ;
	    $ybox = $yhigh - $ylow ;
	    $zbox = $zhigh - $zlow ;
	}
	if ( /ITEM: ATOMS/ ) {
	    printf GEN ("%5d S\n", $natoms) ;
	    for ( $k = 0 ; $k <= $#element ; $k++ ) {
		print GEN '"', $element[$k], '" ' ;
	    }
	    print GEN "\n" ;
	    for ( $j = 0 ; $j < $natoms ; $j++ ) {
		$_ = <TRAJ> ;
		($k, $type, $x, $y, $z) = split(" ") ;
		$x *= $xbox ;
		$y *= $ybox ;
		$z *= $zbox ;

		$x += $xlow ;
		$y += $ylow ;
		$z += $zlow ;
		
		if ( $k > $natoms || $k <= 0 ) {
		    die "Inconsistent atom index found $k\n" ;
		}
		printf GEN ("%4d %4d %9.5f %9.5f %9.5f\n", $j+1, $type, $x, $y, $z) ;
	    }
	    printf GEN ("%9.5f %9.5f %9.5f\n", 0.0, 0.0, 0.0) ;
	    printf GEN ("%9.5f %9.5f %9.5f\n", $xbox, 0.0, 0.0) ;
	    printf GEN ("%9.5f %9.5f %9.5f\n", 0.0, $ybox, 0.0) ;
	    printf GEN ("%9.5f %9.5f %9.5f\n", 0.0, 0.0, $zbox) ;
	}
    }
}
