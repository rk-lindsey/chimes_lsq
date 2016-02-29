#!/usr/bin/perl

open (REF, "<force.txt")  || die "Need file: force.txt\n";
open (CLC, "<f_poly.dat") || die "Need file: f_poly.dat\n";
$df = 0;
while (<CLC>) {
  @line = split " ";
  $fc = $line[0];
  $line = <REF>;
  chomp($line);
  $fr = $line;
  $df += sqrt(($fc-$fr)**2);
  $dfc = sqrt(($fc-$fr)**2);
}
$df /= $.;
print "RMS error of force check of coefficients: $df\n"; 
close (REF);
close (CLC);
