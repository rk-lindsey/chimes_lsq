#!/usr/bin/perl
#
# This program translates an OPENMX md file into DFTB gen format.
# It assumes an orthorhombic box.
#
# Usage: omx_md_to_gen.pl file boxx boxy boxz

use strict ;
use warnings ;

my ($file, $lx, $ly, $lz, $natom, $time_step,
    %ele_idx, $i) ;

if ( scalar(@ARGV) != 4 ) {
    die "Usage: omx_md_to_gen.pl file boxx boxy boxz\n" ;
}
	   
open (INF, "<$ARGV[0]") || die "Can't open input dump file\n";
$file = $ARGV[0];
$file =~ s/.md//;
$file .= ".gen";
open (OUT, ">$file") || die "Can't open output $file\n";

($lx, $ly, $lz) = ($ARGV[1], $ARGV[2], $ARGV[3]) ;
if ( $lx <= 0.0 || $ly <= 0.0 || $lz <= 0.0 ) {
    die "Box dimensions must be positive\n" ;
}

my ( @ele, @x, @y, @z, @idx, @elements) ;

while ( $natom = <INF> ) {
    $time_step = <INF>;
    chomp($time_step);
    chomp($natom);
    for ($i = 0; $i < $natom; $i++) {
	$_ = <INF> ;
	($ele[$i], $x[$i], $y[$i], $z[$i]) = split(" ") ;
	if ( ! defined( $ele_idx{$ele[$i]} ) ) {
	    $ele_idx{$ele[$i]} = scalar(keys(%ele_idx))+1 ;
	    push(@elements, $ele[$i])
	}
	$idx[$i] = $ele_idx{$ele[$i]} ;
    }
    print (OUT "  $natom S # $time_step\n");
    for ( $i = 0 ; $i < scalar(@elements) ; $i++ ) {
	if ( $i < scalar(@elements)-1 ) {
	    print OUT "$elements[$i] " ;
	} else {
	    print OUT "$elements[$i]\n" ;
	}
    }
    for ($i = 0; $i < $natom; $i++) {
	printf(OUT "%4d %4d %12.6f %12.6f %12.6f\n",
	       $i+1, $idx[$i], $x[$i], $y[$i], $z[$i]) ;
    }
    printf (OUT "%12.6f %12.6f %12.6f\n",0,0,0);
    printf (OUT "%12.6f %12.6f %12.6f\n",$lx,0,0);
    printf (OUT "%12.6f %12.6f %12.6f\n",0,$ly,0);
    printf (OUT "%12.6f %12.6f %12.6f\n",0,0,$lz);
}

close (INF);
close (OUT);
