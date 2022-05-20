#! /usr/bin/perl
#
# Translate an XYZ file to a LAMMPS data file.
use strict ;
my ($natom,$x, $y, $z, $j, $ntype,$ntypes, $h, $atom) ;
my ($xbox, $ybox, $zbox) ;

$xbox = 30.0 ;
$ybox = 30.0 ;
$zbox = 30.0 ;

#my %type = (
#    "C" => 1,
#    "N" => 2,
#    "O" => 3,
#    "H" => 4
#) ;

# Modify this hash to set the LAMMPS atom type numbers.
my %type = (
    "O" => 1,
    "H" => 2
) ;

#Mass of atoms
my %mass = ( 
    "C" => 12,
    "O" => 15.999,
    "N" => 14,
    "H" => 1.008 ) ;

$natom = <> ;
chomp $natom ;

print "\# Lammps atomic data from xyz2lammps.pl\n\n" ;
print "$natom atoms\n\n" ;

print "0.0 $xbox xlo xhi\n" ;
print "0.0 $ybox ylo yhi\n" ;
print "0.0 $zbox zlo zhi\n\n" ;

$_ = <> ;

$ntypes = scalar keys(%type) ;
print "$ntypes atom types\n\n" ;
    
print "Masses\n\n" ;

foreach $h ( keys(%type) ) {
    printf("%d %12.5f\n", $type{$h}, $mass{$h}) ;
}
print "\n" ;

print "Atoms\n\n" ;

for ( $j = 1 ; $j <= $natom ; $j++ ) {
    $_ = <> || die "File ended prematurely\n" ;
    ($atom, $x, $y, $z) = split(" ") ;
    printf("%5d %3d   0.0 %12.5f %12.5f %12.5f\n", $j, $type{$atom}, $x, $y, $z) ;
}

