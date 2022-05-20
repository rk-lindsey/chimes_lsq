#! /usr/bin/perl
#
# Program to wrap initial coordinates for OpenMX cartesian
# input into the first periodic simulation cell.  A supercell
# of the original system can be created by specifying sx, sy, and sz.
#
# Usage:  ./wrap_openmx.pl sx sy sz input.old > input.new
# 
# sx, sy, and sz are the number of supercells in the x, y, and z
# directions.  They are optional arguments.
# input.old is the original input file.
# input.new is the new input file with coordinates wrapped into the
# first (primitive) periodic simulation cell.

use strict ;
use warnings ;
use POSIX ;

my ($idx, $ele, $x, $y, $z, $nu, $nd, $xbox, $ybox, $zbox, $file) ;
my ($nx, $ny, $nz) ;
my ($sx, $sy, $sz) ;
my $eps = 1.0e-10 ;
sub wrap_box ;

if ( scalar(@ARGV) == 1 ) {
    $file = $ARGV[0] ;
    $sx = $sy = $sz = 1 ;
} elsif ( scalar(@ARGV) == 4 ) {
    ($file, $sx, $sy, $sz) = @ARGV ;
} else {
    die "Usage: ./wrap_openmx.pl input.old sx sy sz > input.new\n" ;
}


# Read input file.
open(IN, "<$file" ) || die "Could not open file $file\n" ;

$xbox = -1.0 ;
$ybox = -1.0 ;
$zbox = -1.0 ;

while ( <IN> ) {
    if ( /<Atoms\.UnitVectors/ && (!
	/<Atoms\.UnitVectors\.Unit/)  ) {
	# Scale up size of box for supercell if necessary.
	for ( my $j = 0 ; $j < 3 ; $j++ ) {
	    my ($xx, $yy, $zz) ;
	    $_ = <IN> || die "Input ended prematurely\n" ;
	    ($xx, $yy, $zz) = split(" ") ;
	    if ( abs($xx) > $eps ) {
		$xbox = $xx ;
		die "Error: non-orthorhombic box\n" 
		    if ( abs($yy) > $eps || abs($zz) > $eps ) ;
	    } elsif ( abs($yy) > $eps ) {
		$ybox = $yy ;
		die "Error: non-orthorhombic box\n" 
		    if ( abs($xx) > $eps || abs($zz) > $eps ) ;
	    } elsif ( abs($zz) > $eps ) {
		$zbox = $zz ;
		die "Error: non-orthorhombic box\n" 
		    if ( abs($xx) > $eps || abs($yy) > $eps ) ;
	    }
	}
	last ;
    }
}
if ( $xbox <= 0.0 || $ybox <= 0.0 || $xbox <= 0.0 ) {
    die "Error finding box: xbox, ybox, and zbox must be positive !\n" ;
}
close(IN) ;

open(IN, "<$file" ) || die "Could not open file\n" ;

my $count = 1 ;
while ( <IN> ) {
    if ( /<Atoms.SpeciesAndCoordinates/ ) {
	print ;
	while ( <IN> ) { 
	    if ( /Atoms.SpeciesAndCoordinates>/ ) {
		# End of species section.
		print ;
		last ;
	    } elsif ( /^\s*#/ ) {
		# Comment line.
		print ;
	    } elsif ( ! /\S/ ) {
		# White space line.
		print ;
	    } else {
		# A coordinate line.
		($idx, $ele, $x, $y, $z, $nu, $nd) =
		    split(" ") ;
		$x = wrap_box($x, $xbox) ;
		$y = wrap_box($y, $ybox) ;
		$z = wrap_box($z, $zbox) ;

		for ( my $lx = 0 ; $lx < $sx ; $lx++ ) {
		    my $xout = $x + $lx * $xbox ;
		    for ( my $ly = 0 ; $ly < $sy ; $ly++ ) {
			my $yout = $y + $ly * $ybox ;
			for ( my $lz = 0 ; $lz < $sz ; $lz++ ) {
			    my $zout = $z + $lz * $zbox ;
			    printf("%3d %2s %16.13f %16.13f %16.13f %3.1f %3.1f\n",
				   $count, $ele, $xout, $yout, $zout, 
				   $nu, $nd ) ;
			    $count++ ;
			}
		    }
		}
	    }
	}
    } elsif ( /<Atoms.UnitVectors/ ) {
	print ;
	# Scale up size of box for supercell if necessary.
	for ( my $j = 0 ; $j < 3 ; $j++ ) {
	    my ($xx, $yy, $zz) ;
	    $_ = <IN> || die "Input ended prematurely\n" ;
	    ($xx, $yy, $zz) = split(" ") ;
	    printf("%16.13f %16.13f %16.13f\n",
		   $sx * $xx, $sy * $yy, $sz * $zz) ;
	}
    } elsif ( /Atoms.Number\s/ ) {
	my ( $label, $natm ) ;
	($label, $natm) = split(" ") ;
	$natm *= $sx * $sy * $sz ;
	printf("%s %5d\n", $label, $natm) ;
    } else {
	# Pass through all other lines to the output file.
	print ;
    }
}

sub wrap_box()
# Wrap the given coordinate back into the box, 
# arguments: coordinate and box size.
# returns: The new coordinate.
{
    my ( $x, $nx, $box ) ;

    $x = $_[0] ;
    $box = $_[1] ;

    $nx = POSIX::floor($x/$box) ;
#    printf("X = $x, nx = $nx\n") ;

    $x -= $nx * $box ;
    if ( $x < 0.0 ) {
	die "Negative value found\n" ;
    } elsif ( $x > $box ) {
	die "Value bigger than box\n" ;
    }
    return($x) ;
}

    
    
