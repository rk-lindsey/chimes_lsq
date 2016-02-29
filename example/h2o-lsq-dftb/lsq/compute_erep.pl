#!/usr/bin/perl
use constant nparm => 8;
use constant autoang => 0.5291772488820865;
use constant nlayers => 1;

sub nint {
  my $x = $_[0];
  my $n = int($x);
  if ( $x > 0 ) {
    if ( $x-$n > 0.5) {
      return $n+1;
    }
    else {
      return $n;
    }
  }
  else {
    if ( $n-$x > 0.5) {
      return $n-1;
    }
    else {
      return $n;
    }
  }
}

open (INF,  "<input.xyzf") || die "Can't open input.xyzf.\n";
open (PARM, "<$ARGV[0]")   || die "Need to input parameter file: $ARGV[0].\n"; 
open (OUT,  ">f_poly.dat") || die "Can't open f_poly.dat\n";
if ($ARGV[1] =~ /^\d+$/) {
  $nconf = $ARGV[1];
} else {
  die "Need positive integer for nconf.\n";
}
for ($i = 0; $i < 5; $i++) {
  $line = <PARM>; # junk lines
}
@line = split " ", <PARM>;
$npair = $line[$#line];
for ($i = 0; $i < $npair; $i++) {
  @line = split " ", <PARM>;
  $rc[$i] = $line[1];
}
$line = <PARM>; # junk line
for ($i = 0; $i < nparm*$npair; $i++) {
  @line = split " ", <PARM>;
  $parm[$i] = $line[1];
}
$conf = 0;
while (<INF>) {
  $conf++;
  @fx = ();
  @fy = ();
  @fz = ();
  if ($conf > $nconf) {
    last;
  }
  @line = split " ";
  $natom = $line[0];
  @line = split " ", <INF>;
  $la = $line[0];
  $lb = $line[1];
  $lc = $line[2];
  for ($i = 0; $i < $natom; $i++) {
    @line = split " ", <INF>;
    $type[$i] = $line[0];
    $rx[$i] = $line[1];
    $ry[$i] = $line[2];
    $rz[$i] = $line[3];
  }
  for ($i = 0; $i < $natom -1; $i++) {
    for ($j = $i+1; $j < $natom; $j++) {
      if (($type[$i] =~ /O/) && ($type[$j] =~ /O/)) {
        $pair_index = 0;
      } elsif ((($type[$i] =~ /H/) && ($type[$j] =~ /O/)) || (($type[$i] =~ /O/) && ($type[$j] =~ /H/))) {
        $pair_index = 1;
      } elsif (($type[$i] =~ /H/) && ($type[$j] =~ /H/)) {
        $pair_index = 2;
      } 
      $vstart = $pair_index*nparm;
      $dx = $rx[$i] - $rx[$j];
      $dx -= $la*nint($dx/$la);
      $dy = $ry[$i] - $ry[$j];
      $dy -= $lb*nint($dy/$lb);
      $dz = $rz[$i] - $rz[$j];
      $dz -= $lc*nint($dz/$lc);
      for ($n1=-1*nlayers;$n1<nlayers+1;$n1++) {
        for ($n2=-1*nlayers;$n2<nlayers+1;$n2++) {
          for ($n3=-1*nlayers;$n3<nlayers+1;$n3++) {
            $Rab[0]=$dx+$n1*$la;
            $Rab[1]=$dy+$n2*$lb;
            $Rab[2]=$dz+$n3*$lc;
            $rij = sqrt($Rab[0]**2 + $Rab[1]**2 + $Rab[2]**2);
            if ($rij <= $rc[$pair_index]) {
              $rijp = $rij/autoang; 
              $rcp = $rc[$pair_index]/autoang;
              $c2 = $parm[$vstart];
              $c3 = $parm[$vstart+1];
              $c4 = $parm[$vstart+2];
              $c5 = $parm[$vstart+3];
              $c6 = $parm[$vstart+4];
              $c7 = $parm[$vstart+5];
              $c8 = $parm[$vstart+6];
              $c9 = $parm[$vstart+7];
              $force =  2*$c2*($rcp-$rijp);
              $force += 3*$c3*($rcp-$rijp)**2;
              $force += 4*$c4*($rcp-$rijp)**3;
              $force += 5*$c5*($rcp-$rijp)**4;
              $force += 6*$c6*($rcp-$rijp)**5;
              $force += 7*$c7*($rcp-$rijp)**6;
              $force += 8*$c8*($rcp-$rijp)**7;
              $force += 9*$c9*($rcp-$rijp)**8;
              $fx[$i] += $force*$dx/$rij;
              $fx[$j] -= $force*$dx/$rij;
              $fy[$i] += $force*$dy/$rij;
              $fy[$j] -= $force*$dy/$rij;
              $fz[$i] += $force*$dz/$rij;
              $fz[$j] -= $force*$dz/$rij;
            }
          }
        }
      }
    }
  }
  for ($i = 0; $i < $natom; $i++) {
    print (OUT "$fx[$i]\n");
    print (OUT "$fy[$i]\n");
    print (OUT "$fz[$i]\n");
  }
}
close (INF);
close (OUT);
system("./check_erep.pl");  
