#!/usr/bin/perl

system ("pwd");
open (INF, "<$ARGV[0]") || die "Need input file.\n";
open (OUT, ">input_dftb.xyzf") || die "Can't open input_dftb.xyzf.\n";
shift(@ARGV); # remove file name from ARGV
if ($ARGV[0] =~ /^\d+$/) {
  $nprint = $ARGV[0];
  shift(@ARGV);
} else {
  $nprint = 1;
}
print "nprint: $nprint\n";
$len = @ARGV;
if ($len == 0) {
  die "Need to input atom list\n";
}
print "$len elements in list.\n";
$atom_types = join " ", @ARGV;
print "atom types: $atom_types\n";

$frame = 1;
while (<INF>) {
  print "frame: $frame\n";
  $frame++;
  @line = split " ";
  $natom = $line[0];
  @line = split " ", <INF>;
  $la = $line[0];
  $lb = $line[1];
  $lc = $line[2];
  print (OUT "$natom\n");
  print (OUT "$la $lb $lc\n");
  open (TMP, ">tmp.gen") || die "Can't open tmp.gen.\n";
  print (TMP "$natom S\n");
  print (TMP "$atom_types\n");
  for ($i = 0; $i < $natom; $i++) {
    @line = split " ", <INF>; 
    $type[$i] = $line[0];
    $rx[$i] = $line[1];
    $ry[$i] = $line[2];
    $rz[$i] = $line[3];
    $fx_dft[$i] = $line[4];
    $fy_dft[$i] = $line[5];
    $fz_dft[$i] = $line[6];
    $index = $i+1;
    if ($type[$i] =~ /O/) {
      $itype = 1;
    } elsif ($type[$i] =~ /H/) {
      $itype = 2;
    } elsif ($type[$i] =~ /C/) {
      $itype = 3;
    } elsif ($type[$i] =~ /N/) {
      $itype = 4;
    } else {
      die "ERROR in atom types: $type[$i] not supported.\n";
    }
    print (TMP "$index $itype $rx[$i] $ry[$i] $rz[$i]\n");
  }
  print (TMP "0.0000 0.0000 0.0000\n");
  print (TMP "$la 0.0000 0.0000\n");
  print (TMP "0.0000 $lb 0.0000\n");
  print (TMP "0.0000 0.0000 $lc\n");
  close (TMP);
  $CMD = "/usr/gapps/polymers/dftb+_1.2_llnl/prg_dftb/_obj_x86_64-linux-ifort/dftb+ > tmp.out";
  print "$CMD\n";
  system($CMD);
  open (DET, "<detailed.out") || die "Can't open detailed.out\n";
  while ($data = <DET>) {
    if ($data =~ /Total Forces/) {
      for ($i = 0; $i < $natom; $i++) {
        @line = split " ", <DET>;
        $fx_dftb[$i] = $line[0];
        $fy_dftb[$i] = $line[1];
        $fz_dftb[$i] = $line[2];
      }
      last;
    }
  }
  close (DET);
  for ($i = 0; $i < $natom; $i++) {
    $dfx = $fx_dft[$i] - $fx_dftb[$i];
    $dfy = $fy_dft[$i] - $fy_dftb[$i];
    $dfz = $fz_dft[$i] - $fz_dftb[$i];
    print (OUT "$type[$i] $rx[$i] $ry[$i] $rz[$i] $dfx $dfy $dfz\n");
  }
}
close (INF);
close (OUT);
