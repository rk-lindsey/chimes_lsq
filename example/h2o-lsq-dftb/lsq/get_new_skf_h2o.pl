#!/usr/bin/perl
use constant autoang => 0.52917706;
use constant nfile  => 4;

open (PARM, "<$ARGV[0]") || die "Can't open parameter file\n";
for ($i = 0; $i < 5; $i++) {
  $line = <PARM>; # junk lines
}
@line = split " ", <PARM>;
$npair = $line[1];
for ($i = 0; $i < $npair; $i++) {
  @line = split " ", <PARM>;
  $rc[$i] = $line[1];
}
$line = <PARM>; # junk
@param = ();
for ($i = 0; $i < $npair; $i++) {
  for ($j = 0; $j < 8; $j++) {
    @line  = split " ", <PARM>;
    $param[$i][$j] = $line[$#line];
  }
}
close (PARM);

$file[0] = "O-O.skf";
$file[1] = "O-H.skf";
$file[2] = "H-O.skf";
$file[3] = "H-H.skf";
$mass[0] = 15.994;
$mass[1] = 1.0;
$mass[2] = 1.0;
$mass[3] = 1.008;
for ($i = 0; $i < nfile; $i++) {
  $rcut = $rc[$i]/autoang;
  $name = $file[$i];
  $data = $name;
  $data =~ s/.skf//;
  @atom = split "-", $data;
  if ($atom[0] !~ $atom[1]) {
    $hetero = 1;
    $erep_line = 2;
  } else {
    $hetero = 0;
    $erep_line = 3;
  }
  if (($atom[0] =~ /O/) && ($atom[1] =~ /O/)) {
    $pair_index = 0;
  } elsif ((($atom[0] =~ /H/) && ($atom[1] =~ /O/)) || (($atom[0] =~ /O/) && ($atom[1] =~ /H/))) {
    $pair_index = 1; 
  } elsif (($atom[0] =~ /H/) && ($atom[1] =~ /H/)) {
    $pair_index = 2;
  }
  system ("cp ~/slako/mio-1-1/$name .");
  open (INF, "<$name") || die "Can't open $name\n";
  open (OUT, ">tmp")   || die "Can't open tmp\n";
  while (<INF>) {
    if ($. == $erep_line) {
      print (OUT "$mass[$i] $param[$pair_index][0] $param[$pair_index][1] $param[$pair_index][2] $param[$pair_index][3] $param[$pair_index][4] $param[$pair_index][5] $param[$pair_index][6] $param[$pair_index][7] $rcut 3.5, 9*0.0,\n");
    } else {
      print (OUT "$_");
    }
  }
  system ("mv tmp $name");
  close (INF);
  close (OUT);
}
