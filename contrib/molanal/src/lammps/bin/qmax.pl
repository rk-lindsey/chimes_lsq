#! /usr/bin/perl
# Calculates maximum charge for a LAMMPS q x y z dump file.
# Usage:  qmax.pl <type to find> < dump_file > qmax_file
# If the atom type to find is not specified, the maximum charge
# among all atom types is found.

use strict ;
use warnings ;

sub nint ;

my ($qmax, $q, $timestep,$natom,$indx,$type,$indx_max,$j,$tfind,
    $type_max,$qavg,$qcount,$xbox, $ybox, $zbox) ;
my $frame = 0 ;
my $dump_neighbor = 1 ;
my $rneighbor = 3.0 ;
my $frame_skip = 1000 ;

# How often dump files were written.
my $dump_freq = 20 ;

my @element_type = qw/C N O H/ ;

if ( scalar(@ARGV) > 0 ) {
    $tfind = $ARGV[0] ;
    shift @ARGV ;
    print "Looking for type $tfind\n" ;
    if ( ! ($tfind =~ /^[0-9]+/)  ) {
	die "Expecting an integer atom type\n" ;
    }
} else {
    $tfind = -1 ;
}
    
printf("Step  Atom  Type   Max |Q|   RMS Q\n") ;

while ( <> ) {
    if ( /ITEM: TIMESTEP/ ) {
	$timestep = <> ;
	chop $timestep ;
#	print "Checking timestep $timestep\n" ;
    }
    if ( /ITEM: BOX BOUNDS/ ) {
	my ($tmp,$lo, $hi) ;
	$_ = <> ;
	($lo, $hi) = split " " ;
	$xbox = $hi - $lo ;
	$_ = <> ;
	($lo, $hi) = split " " ;
	$ybox = $hi - $lo ;
	$_ = <> ;
	($lo, $hi) = split " " ;
	$zbox = $hi - $lo ;
    }
    if ( /ITEM: NUMBER OF ATOMS/ ) {
	$natom = <> ;
#	print "Number of atoms = $natom\n" ;
    }
    if ( /ITEM: ATOMS id type q/ ) {
	$qcount = 0 ;
	$qmax = 0 ;
	$qavg = 0 ;
	$type_max = 0 ;
	$frame++ ;
	my (@x, @y, @z,@qsave,@element,$xx,$yy,$zz) ;
	for ( $j = 0 ; $j < $natom ; $j++ ) {
	    $_ = <> ;	
	    ($indx, $type, $q,$xx,
	     $yy,$zz) = split(" ") ;
	    $x[$indx] = $xx ;
	    $y[$indx] = $yy ;
	    $z[$indx] = $zz ;
	    $qsave[$indx] = $q ;
	    $element[$indx] = $element_type[$type-1] ;
	    if ( $tfind >= 0 ) {
		if ( $type == $tfind ) {
		    if ( abs($q) > $qmax ) {
			$qmax = abs($q) ;
			$indx_max = $indx ;
			$type_max = $type ;
		    }
		    $qavg += $q * $q ;
		    $qcount++ ;
		}
	    } else {
		if ( abs($q) > $qmax ) {
		    $qmax = abs($q) ;
		    $indx_max = $indx ;
		    $type_max = $type ;
		    
		}
		$qavg += $q * $q ;
		$qcount++ ;
	    }
	}
	if ( $qcount > 0 ) {
	    $qavg /= $qcount ;
	    $qavg = sqrt($qavg) ;
	}
	printf("%7d %7d %4d %11.4f %11.4f\n", $timestep, $indx_max, 
	       $type_max, $qmax, $qavg) ;

	if ( $dump_neighbor && $frame % $frame_skip == 0 ) {
	    dump_xyz(\@x,\@y,\@z,\@qsave,\@element,$xbox,$ybox,$zbox,
		     $rneighbor,$frame,$dump_freq) ;
	}
    }
}



sub nint
# Nearest Integer Function.
{
    my $v ;
    $v = $_[0] ;

    die "Too many arguments to nint" if ( scalar (@$) > 1 ) ;

    if ( $v > 0.0 ) {
	return ( int($v + 0.5) ) ;
    } else {
	return ( int($v - 0.5) ) ;
    }
}

sub wrap_box
# Wrap atoms into the periodic simulation box.  Returns the wrapped
# coordinates.
{
    my ($x,$y,$z,$xbox,$ybox,$zbox) = @_ ;
    my ($dx, $dy, $dz) ;

    $dx = $xbox * nint($x/$xbox) ;
    $dy = $ybox * nint($y/$ybox) ;
    $dz = $zbox * nint($z/$zbox) ;

    $x -= $dx ;
    $y -= $dy ;
    $z -= $dz ;
    
    return($x,$y,$z) ;
}

sub dump_xyz
# Write an xyz file showing the neighborhood around the atom in the frame
# that has maximum charge.  All atoms within rneighbor of the atom
# with maximum charge will be output.  Atoms will be wrapped within
# period boundaries, which are assumed to be orthogonal.
{

    my ($x,$y,$z,$qsave,$element, $xbox,$ybox,$zbox,$rneighbor,
	$frame,$dump_freq) = @_ ;
    die "Incorrect args to dump_xyz" if scalar @_ != 11 ;
    my $atom_count = 0 ;
    my $lines ;
    my $filename = sprintf("qmax.neighbor.%07d.xyz", 
			   $frame*$dump_freq) ;
    open( my $frame_file, ">", $filename) 
	|| die "Could not open $filename\n" ;
    print STDERR "Opened file $filename\n" ;
    ($$x[$indx_max],$$y[$indx_max],$$z[$indx_max]) =
	wrap_box($$x[$indx_max],$$y[$indx_max],$$z[$indx_max],
		 $xbox, $ybox, $zbox) ;
    $lines = sprintf( "%2s %11.4e %11.4e %11.4e %11.4e\n",
		      $$element[$indx_max], 
		      $$x[$indx_max], 
		      $$y[$indx_max], 
		      $$z[$indx_max],
		      $$qsave[$indx_max]) ;
    for ( $j = 1 ; $j <= $natom ; $j++ ) {
	next if ( $j == $indx_max ) ;
	my $dx = $$x[$j] - $$x[$indx_max] ;
	my $dy = $$y[$j] - $$y[$indx_max] ;
	my $dz = $$z[$j] - $$z[$indx_max] ;

	my $dx2 = $xbox * nint($dx/$xbox) ;
	my $dy2 = $ybox * nint($dy/$ybox) ;
	my $dz2 = $zbox * nint($dz/$zbox) ;

	($dx,$dy,$dz) =
	    wrap_box($dx,$dy,$dz,$xbox,$ybox,$zbox);
	my $r2 = sqrt($dx*$dx + $dy*$dy + $dz*$dz) ;

	die "Element[$j] undefined\n" if ! defined($$element[$j]) ;
	if ( $r2 < $rneighbor ) {
	    $lines .= sprintf( "%2s %11.4e %11.4e %11.4e %11.4e\n",
			       $$element[$j], 
			       $$x[$j]-$dx2, 
			       $$y[$j]-$dy2,
			       $$z[$j]-$dz2,
			       $$qsave[$j]) ;
	    $atom_count++ ;
	    $r2 = norm3d($$x[$j]-$dx2-$$x[$indx_max], 
			 $$y[$j]-$dy2-$$y[$indx_max],
			 $$z[$j]-$dz2-$$z[$indx_max]) ;
	    if ( $r2 > $rneighbor + 1.0e-12 ) {
		die "Failed check on atom distance\n" ;
	    }
	}
    }
    print $frame_file "$atom_count\n\n" ;
    print $frame_file $lines ;
    close ($frame_file) ;
}

sub norm3d
{
    my ($x, $y, $z) = @_ ;

    return(sqrt($x*$x+$y*$y+$z*$z) ) ;
}


