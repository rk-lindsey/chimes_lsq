#!/usr/bin/perl
use Switch;

$cleanup = 0;
$CMD = "cp ../input.xyz .";
print "$CMD\n";
system($CMD);
$exec = "../../../bin/splines_md orig_3b_force.in";
print "Running reference calculation...\n";
print "$exec > ref.out";
system("$exec > ref.out");
print "\n";
system("mv 3b_results.dat 3b_ref.dat");
open (R3B, "<3b_ref.dat");
$line = <R3B>; # read in energy 
$rcount = 0;
@fref = ();
while (<R3B>) {
  chomp($_);
  $fref[$rcount] = $_;
  $rcount++;
}
$dr = 1e-5;
if (! -e "ref.xyz") {
  system("cp input.xyz ref.xyz");
}
open (REF, "<ref.xyz")   || die "Can't open ref.xyz\n";
@line = split " ", <REF>;
$natom = $line[0];
$second_line = <REF>;
for ($i = 0; $i < $natom; $i++) {
  @line = split " ", <REF>;
  $type[$i] = $line[0];
  $xc[$i] = $line[1];
  $yc[$i] = $line[2];
  $zc[$i] = $line[3];
  $vx[$i] = $line[4];
  $vy[$i] = $line[5];
  $vz[$i] = $line[6];
}
open (NFC, ">num_forces.dat") || die "Can't open num_forces.dat\n";
@fnum = ();
$ncount = 0;
for ($tag = 0; $tag < $natom; $tag++) {
  for ($c = 0; $c < 3; $c++) {
    switch ($c) {
      case(0) {
        $dx = $dr;
        $dy = 0;
        $dz = 0;
        $suffix = "x";
      }
      case(1) {
        $dx = 0;
        $dy = $dr;
        $dz = 0;
        $suffix = "y";
      }
      case(2) {
        $dx = 0;
        $dy = 0;
        $dz = $dr;
        $suffix = "z";
      }
    }
    $xc[$tag] += $dx;
    $yc[$tag] += $dy;
    $zc[$tag] += $dz;
    open (INP, ">input.xyz") || die "Can't open input.xyz\n";
    print (INP "$natom\n");
    print (INP "$second_line");
    for ($i = 0; $i < $natom; $i++) {
      print (INP "$type[$i] $xc[$i] $yc[$i] $zc[$i] $vx[$i] $vy[$i] $vz[$i]\n");
    }
    close (INP);
    $CMD = "$exec > p$tag\_$suffix.out";
    print "$CMD";
    system($CMD);
    print "\n";
    @line = split " ", `head -n 1 3b_results.dat`;
    $eP = $line[$#line];
    system("mv 3b_results.dat 3b_p$tag\_$suffix.dat");
    $xc[$tag] -= 2*$dx;
    $yc[$tag] -= 2*$dy;
    $zc[$tag] -= 2*$dz;
    open (INP, ">input.xyz") || die "Can't open input.xyz\n";
    print (INP "$natom\n");
    print (INP "$second_line");
    for ($i = 0; $i < $natom; $i++) {
      print (INP "$type[$i] $xc[$i] $yc[$i] $zc[$i] $vx[$i] $vy[$i] $vz[$i]\n");
    }
    $CMD = "$exec > n$tag\_$suffix.out";
    print "$CMD";
    system($CMD);
    print "\n";
    @line = split " ", `head -n 1 3b_results.dat`;
    $eN = $line[$#line];
    system("mv 3b_results.dat 3b_n$tag\_$suffix.dat");
    close (INP);
    $force = -($eP-$eN)/(2*$dr);
    print (NFC "$force\n");
    $fnum[$ncount] = $force;
    $ncount++;
  }
}
$rms = 0;
$len_ref = @fref;
$len_num = @fnum;
if ($len_ref != $len_num) {
  die "len_ref != $len_num: $len_ref $len_num\n";
}
for ($i = 0; $i < $len_ref; $i++) {
  $rms = ($fref[$i] - $fnum[$i])**2;
}
$rms /= $len_ref;
$rms = sqrt($rms);
print "RMS error in forces = $rms\n";
close (NFC);
if ($cleanup) {
  system("rm -f *.out 3b*.dat ref.xyz input.xyz");  
}
