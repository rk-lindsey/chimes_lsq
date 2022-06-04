#! /usr/bin/perl
# Calculates minimum distance between atoms in rdf file.
# Usage:  rmin.pl rdf_file > output

use strict ;

my ($j, $step, $nrdf, $r, $bin, $pop, $rmin,
    $file, $outname) ;


foreach $file ( @ARGV ) {
    print "Scanning $file\n" ;
    open (IN, "<$file" ) || die "Could not open $file\n" ;
    $outname = $file ;
    $outname =~ s/\.rdf// ;
    $outname = $outname . ".rmin" ;
    open (OUT, ">$outname") || die "Could not open $outname\n" ;
    print "Opened $outname\n" ;
    print OUT "  Step         Rmin\n" ;
    for ( $j = 0 ; $j < 3 ; $j++ ) {
	$_ = <IN> || die "Could not read header\n" ;
    }
    while ( <IN> ) {
	($step, $nrdf) = split(" ") ;
	$rmin = 0.0 ;
	for ( $j = 0 ; $j < $nrdf ; $j++ ) {
	    $_ = <IN> || die "File ended prematurely\n" ;
	    ($bin, $r, $pop) = split(" ") ;
	    if ( $rmin == 0.0 && $pop > 0 ) {
		$rmin = $r ;
	    }
	}
	printf OUT ("%9d %11.4f\n", $step, $rmin) ;
    }
}
