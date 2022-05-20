#! /usr/bin/perl
# Calculates the system dipole moment LAMMPS q xu yu zu dump file.
# Usage:  dipole.pl < dump_file > qmax_file

use strict ;

my ($qmax, $q, $timestep,$natom,$indx,$type,$indx_max,$j) ;
my ($xlo, $xhi, $ylo, $yhi, $zlo, $zhi, $xbox, $ybox, $zbox) ;
my ($x, $y, $z, $ix, $iy, $iz, $dipx, $dipy, $dipz,$dipole) ;

my $eps = 1.0e-10 ;

# 1 Debye = 0.208 1 e Angstrom 
my $debye = 0.208 ;

printf("Step  Dipx  Dipy  Dipz\n") ;

while ( <> ) {
    if ( /ITEM: TIMESTEP/ ) {
	$timestep = <> || die "Could not read timestep\n" ;;
	chop $timestep ;
#	print "Checking timestep $timestep\n" ;
    }
    if ( /ITEM: BOX BOUNDS/ ) {
	$_ = <> || die "Could not read box bounds\n" ;
	($xlo, $xhi) = split(" ") ;
	$_ = <> || die "Could not read box bounds\n" ;
	($ylo, $yhi) = split(" ") ;
	$_ = <> || die "Could not read box bounds\n" ;
	($zlo, $zhi) = split(" ") ;
	$xbox = $xhi - $xlo ;
	$ybox = $yhi - $ylo ;
	$zbox = $zhi - $zlo ;
    }
	
    if ( /ITEM: NUMBER OF ATOMS/ ) {
	$natom = <> || die "Could not read number of atoms\n" ;
#	print "Number of atoms = $natom\n" ;
    }
    if ( /ITEM: ATOMS id type q xu yu zu/ ) {
	$dipx = 0.0 ;
	$dipy = 0.0 ;
	$dipz = 0.0 ;
	for ( $j = 0 ; $j < $natom ; $j++ ) {
	    $_ = <> || die "Dump file ended prematurely\n" ;	
	    ($indx, $type, $q, $x, $y, $z) = split(" ") ;
	    $ix = &nint($x/$xbox) ;
	    $iy = &nint($y/$ybox) ;
	    $iz = &nint($z/$zbox) ;

	    $x -= $ix * $xbox ;
	    $y -= $iy * $ybox ;
	    $z -= $iz * $zbox ;
	    
	    if ( $x > 0.5 * $xbox + $eps ) {
		die "x is too big" ;
	    }
	    if ( $x < -0.5 * $xbox - $eps ) {
		die "x is too small" ;
	    }

	    if ( $y > 0.5 * $ybox + $eps ) {
		die "y is too big" ;
	    }
	    if ( $y < -0.5 * $ybox - $eps ) {
		die "y is too small" ;
	    }

	    if ( $z > 0.5 * $zbox + $eps ) {
		die "z is too big" ;
	    }
	    if ( $z < -0.5 * $zbox - $eps ) {
		die "z is too small" ;
	    }

	    $dipx += $q * $x ;
	    $dipy += $q * $y ;
	    $dipz += $q * $z ;
	}
#	$dipx /= $natom ;
#	$dipy /= $natom ;
#	$dipz /= $natom ;
	$dipole = sqrt($dipx * $dipx + $dipy * $dipy + $dipz * $dipz) ;
	printf("%7d %11.4f\n", $timestep, $dipole/$debye) ;
    }
}

sub nint()
# Nearest Integer Function.
{
    my $v ;
    $v = $_[0] ;

    if ( $v > 0.0 ) {
	return ( int($v + 0.5) ) ;
    } else {
	return ( int($v - 0.5) ) ;
    }
}
