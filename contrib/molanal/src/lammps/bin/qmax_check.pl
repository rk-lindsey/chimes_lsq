#! /usr/bin/perl
# Checks maximum charge for a LAMMPS q x y z dump file.
# Usage:  qmax_check.pl dump_file last_timestep 

use strict ;
my ($qmax, $q, $timestep,$natom,$indx,$type,$indx_max,$j,
    $qmax_lim, $qmax_all, $last_timestep) ;

$qmax_lim = 2.99 ;
$qmax_all = 0.0 ;
#printf("Step  Atom  Max |Q|\n") ;

$last_timestep = pop @ARGV ;

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
    if ( /ITEM: ATOMS id type q/ ) {
	$qmax = 0 ;
	for ( $j = 0 ; $j < $natom ; $j++ ) {
	    $_ = <> ;	
	    ($indx, $type, $q) = split(" ") ;
	    if ( abs($q) > $qmax ) {
		$qmax = abs($q) ;
		$indx_max = $indx ;
	    }
	}
#	print OUT ("%7d %7d %11.4f\n", $timestep, $indx_max, $qmax) ;
	if ( $qmax > $qmax_lim ) {
	    die "ERROR: CHARGE WAS TOO LARGE: $qmax\n" ;
	} elsif ( $qmax > $qmax_all ) {
	    $qmax_all = $qmax ;
	}
    }
}
if ( $timestep != $last_timestep ) {
    print "ERROR: RUN DID NOT TERMINATE AT EXPECTED STEP: $timestep\n" ;
}
else {
    print "Dump file terminated at step $timestep (OK)\n" ;
}
print "Maximum charge found = $qmax_all (OK).\n" ;


