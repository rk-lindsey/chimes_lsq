#! /usr/bin/perl
#
# Grab specific atoms from a gen file, and print a new gen file
# with just those atoms to the standard output.
#
# Usage:  grabatom.pl <atom numbers> <gen file> > <new gen file>
#
use strict ;
use warnings ;

my $file  = pop ( @ARGV ) ;
my @atoms = @ARGV ;
my $natom2 = scalar(@atoms) ;
my $fh ;

local $, = " " ;
open ( $fh, "<$file" ) || die "Could not open $file\n" ;

while ( <$fh> ) {
    my $natom1 ;

    ($natom1) = split(" ") ;
    
    $_ =~ s/$natom1/$natom2/ ;
    print ;
    
    $_ = <$fh> ;
    print ;

    for ( my $j = 0 ; $j < $natom1 ; $j++ ) {
	$_ = <$fh> ;
	my @f = split(" ") ;
	my $idx = shift @f ;
	for ( my $k = 0 ; $k < $natom2 ; $k++ ) {
	    if ( $idx == $atoms[$k] + 1 ) {
		print $k+1, @f ;
		print "\n" ;
	    }
	}
    }
    for ( my $j = 0 ; $j < 4 ; $j++ ) {
	$_ = <$fh> ;
	print ;
    }
}

    
