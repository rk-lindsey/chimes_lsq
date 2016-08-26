#! /usr/bin/perl
#
# This program splits the output of "joinmolecules.pl" into separate
# output files for plotting. Output files will be labeled by molecule
# and/or property, and will end with the suffix .txt.  It hasn't been
# tested yet for the output of findmolecules.pl.  In the meantime,
# you can run a single findmolecules.pl output through joinmolecules.pl,
# and then split it using splitmol.pl
#
# Usage: splitmol.pl <join file>
#
use strict ;
use warnings ;
my $conc_header = "Concentration history for" ;
my $conc_header_lines = 2 ;

my $poly_header = "Polymer history for" ;
my $poly_header_lines = 1 ;

my $molwt_header = "Molecular weight histogram" ;

my $poly_number = "Number of polymer molecules" ;

open(IN,"<$ARGV[0]") || die "Could not read $ARGV[0]\n" ;

while ( <IN> ) {
    if ( /$conc_header/ ) {
	my $name = $_ ;
	
	$name =~ s/$conc_header// ;
	$name =~ s/\s$// ;
	$name =~ s/^ // ;
	$name =~ s/\s/_/g ;
	$name =~ s/\(//g ;
	$name =~ s/\)//g ;

	for ( 1 .. $conc_header_lines ) {
	    $_ = <IN> ;
	}
	open (OUT, ">$name.hist.txt") || die "Could not open $name.hist.txt" ;
	print "Writing $name.hist.txt\n" ;

	print OUT "# $_" ;

	while ( <IN> ) {
	    if ( /^$/ ) {
		close OUT ;
		last ;
	    }
	    print OUT ;
	}
    } elsif ( /$poly_header/ ) {
	my $name = $_ ;
	
	$name =~ s/$poly_header// ;
	$name =~ s/\s$// ;
	$name =~ s/^ // ;
	$name =~ s/\s/_/g ;

	for ( 1 .. $poly_header_lines ) {
	    $_ = <IN> ;
	}
	open (OUT, ">poly.$name.hist.txt") || die "Could not open poly.$name.hist.txt" ;
	print "Writing poly.$name.hist.txt\n" ;

	print OUT "# $_" ;

	while ( <IN> ) {
	    if ( /^$/ ) {
		close OUT ;
		last ;
	    }
	    print OUT ;
	}
    } elsif ( /$poly_number/ ) {
	my $name = "poly.number.hist.txt" ;
	
	for ( 1 .. $poly_header_lines ) {
	    $_ = <IN> ;
	}
	open (OUT, "> $name") || die "Could not open $name\n" ;
	print "Writing $name\n" ;
	print OUT "# $_" ;

	while ( <IN> ) {
	    if ( /^$/ ) {
		close OUT ;
		last ;
	    }
	    print OUT ;
	}
    } elsif ( /$molwt_header/ ) {
	my $name = "molwt.histo.txt" ;
	for ( 1 .. $poly_header_lines ) {
	    $_ = <IN> ;
	}
	open (OUT, "> $name") || die "Could not open $name\n" ;
	print "Writing $name\n" ;
	print OUT "# $_" ;

	while ( <IN> ) {
	    if ( /^$/ ) {
		close OUT ;
		last ;
	    }
	    print OUT ;
	}
    }
}

