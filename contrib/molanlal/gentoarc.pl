#!/usr/bin/perl
## This is a program to translate DFTB gen files from multiple
## simulations into a single Biosym/MSI/Accelrys .arc format file.
## 
## Usage:  Then run
## gentoarc.pl <dir list>.  The output is placed in dftb.arc.
##
################################################################
##
## User defined parameters.
##
## Read in the number of atoms
use Math::Trig ;
    while ( <> ) {
($natoms,$pbcmode) = split(" ",$_) ;
## Read in the elements.
$_ = <> ;
@ele = split(" ",$_) ;
foreach $i (@ele) {
    $i =~ s/\"//g ;
}
$neletype = $#ele + 1 ;
print @ele, "\n", $neletype, "\n" ;
for ( $j = 0 ; $j < $natoms ; $j++ ) {
    $_ = <> ;
    ($idx, $eletype[$j], $x[$j], $y[$j], $z[$j]) = split(" ", $_) ;
    print $x[$j], " ", $y[$j], " ", $z[$j], "\n" ;
}
$_ = <> ;
$_ = <> ;
($ax, $ay, $az) = split(" ",$_) ;
print "$ax $ay $az\n" ;
$_ = <> ;
($bx, $by, $bz) = split(" ",$_) ;
print "$bx $by $bz\n" ;
$_ = <> || die ;
($cx, $cy, $cz) = split(" ",$_) ;
print "$cx $cy $cz\n" ;
$a = sqrt($ax * $ax + $ay * $ay + $az * $az) ;
$b = sqrt($bx * $bx + $by * $by + $bz * $bz) ;
$c = sqrt($cx * $cx + $cy * $cy + $cz * $cz) ;
print "a = $a, b = $b, c = $c\n" ;
if ( abs($bx) > 1.0e-06 || abs($cx) > 1.0e-06 || abs($cy) > 1.0e-06 ) {
    print $bx, " ", $cx, " ", $cy, "\n" ;
    die "Error: unit cell must be oriented so that c is in the z direction\n" ;
}
$adotb = $ax * $bx + $ay * $by + $az * $bz ;
$adotc = $ax * $cx + $ay * $cy + $az * $cz ;
$bdotc = $bx * $cx + $by * $cy + $bz * $cz ;

$pii = pi ;
print "pi = $pii\n" ;
$alpha = acos($bdotc/($b*$c)) * 180.0 / pi ;
print "Alpha = $alpha\n" ;
$beta = acos($adotc/($a*$c)) * 180.0 / pi ;
print "Beta = $beta\n" ;
$gamma = acos($adotb/($a*$b)) * 180.0 / pi ;
print "Gamma = $gamma\n" ;


##
## End of user defined parameters.
##
#################################################################
## Unit cell parameters.
$l = 0 ;
for ( $l = 0 ; $l <= $natoms ; $l++ ) {
    $j = $eletype[$l] - 1 ;
    $label[$l] = $ele[$j] . "$l" ;
    $label3[$l] = $ele[$j] ;
    $label2[$l] = $ele[$j] ;
    $label2[$l] =~ tr/A-Z/a-z/ ;
}
$date=`date` ;
$date =~ s/PST //;
open(OUT,">dftb.arc") ;
print OUT "!BIOSYM archive 3\nPBC=ON\n" ;
$alast = 0.0 ;
print OUT "\n!DATE $date" ;
printf OUT "PBC  %9.5f %9.5f %9.5f %8.4f  %8.4f  %8.4f (P1) \n",
    $a, $b, $c, $alpha, $beta, $gamma ;
for ( $k = 0 ; $k < $natoms ; $k++ ) {
    printf OUT "%-5s%15.9f%15.9f%15.9f%7s%6s  %6s %7.3f\n", $label[$k], 
    $x[$k], $y[$k], $z[$k], " XXX  ND", $label2[$k], $label3[$k], 0.0 ;
}
printf OUT "end\n" ;
}
print OUT "end\n" ;
    
