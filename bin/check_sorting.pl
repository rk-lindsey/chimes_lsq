#! /usr/bin/perl
# Checks the output of spline_md when the TESTING option in Cheby3B.C is turned on.
# Makes sure that all atom permutations are correct.
use strict ;
use warnings ;

my @ele ;
my @n ;
my @a ;
my @pair ;

my @ele_new ;
my @n_new ;
my @a_new ;
my @pair_new ;
while ( <> ) {
  if ( /Input:/ ) {
    $_ = <> || die "End of file" ; 
    chomp ;
    $_ =~ s/ele[1-3]//g ;
    my @f = split("=", $_) ;
    ($ele[0], $ele[1], $ele[2]) = ($f[1], $f[2], $f[3]) ;

    $_ = <> || die "End of file" ; 
    chomp ;
    $_ =~ s/n[1-3]+//g ;
    @f = split("=", $_) ;
    ($n[0], $n[1], $n[2]) = ($f[1], $f[2], $f[3]) ;

    $_ = <> || die "End of file" ; 
    chomp ;
    $_ =~ s/a[1-3]+//g ;
    @f = split("=", $_) ;
    ($a[0], $a[1], $a[2]) = ($f[1], $f[2], $f[3]) ;

    $_ = <> || die "End of file" ; 
    chomp ;
    $_ =~ s/ipair[1-3]+//g ;
    @f = split("=", $_) ;
    ($pair[0], $pair[1], $pair[2]) = ($f[1], $f[2], $f[3]) ;

    print("Input values:\n") ;
    print(" Elements : @ele\n") ;
    print(" Orders   : @n\n") ;
    print(" Atoms    : @a\n") ;
    print(" Pairs    : @pair\n") ;
  }
  if ( /Sorted/ ) {
    $_ = <> || die "End of file" ; 
    chomp ;
    $_ =~ s/ele[1-3]//g ;
    my @f = split("=", $_) ;
    ($ele_new[0], $ele_new[1], $ele_new[2]) = ($f[1], $f[2], $f[3]) ;

    foreach my $i ( 0 .. 1 ) {
      if ( $ele_new[$i] > $ele_new[$i+1] ) {
	die "Elements were not sorted\n" ;
      }
    }

    $_ = <> || die "End of file" ; 
    chomp ;
    $_ =~ s/n[1-3]+//g ;
    @f = split("=", $_) ;
    ($n_new[0], $n_new[1], $n_new[2]) = ($f[1], $f[2], $f[3]) ;

    $_ = <> || die "End of file" ; 
    chomp ;
    $_ =~ s/a[1-3]+//g ;
    @f = split("=", $_) ;
    ($a_new[0], $a_new[1], $a_new[2]) = ($f[1], $f[2], $f[3]) ;

    $_ = <> || die "End of file" ; 
    chomp ;
    $_ =~ s/ipair[1-3]+//g ;
    @f = split("=", $_) ;
    ($pair_new[0], $pair_new[1], $pair_new[2]) = ($f[1], $f[2], $f[3]) ;

    print("Sorted values\n") ;
    print(" Elements : @ele_new\n") ;
    print(" Orders   : @n_new\n") ;
    print(" Atoms    : @a_new\n") ;
    print(" Pairs    : @pair_new\n") ;

    my @matched ;
    foreach my $i ( 0 .. 2 ) {
      $matched[$i] = -1 ;
      foreach my $j ( 0 .. 2 ) {
	if ( $a_new[$i] == $a[$j] ) {
	  $matched[$i] = $j ;
	}
      }
      if ( $matched[$i] == -1 ) {
	die "Error: did not match atom index $i\n" ;
      }
      if ( $ele_new[$i] != $ele[$matched[$i]] ) {
	die "Error: elements did not correspond\n" ;
      }
    }

    if ( $pair_new[0] != pair_index($ele_new[0], $ele_new[1]) ) {
      die "Error: Pair 1 order did not correspond\n" ;
    }
    if ( $pair_new[1] != pair_index($ele_new[0], $ele_new[2]) ) {
      die "Error: Pair 2 order did not correspond\n" ;
    }
    if ( $pair_new[2] != pair_index($ele_new[1], $ele_new[2]) ) {
      die "Error: Pair 3 order did not correspond\n" ;
    }

    my $idx = $matched[0] + $matched[1] - 1 ;
    if ( $n_new[0] != $n[$idx] ) {
      print " Matched $matched[0] $matched[1]\n" ;
      print " N_new $n_new[0]\n" ;
      print " N $n[$idx] \n" ;
      die "1st Polynomial order did not correspond\n" ;
    }

    $idx = $matched[0] + $matched[2] - 1 ;
    if ( $n_new[1] != $n[$idx] ) {
      print " Matched $matched[0] $matched[2]\n" ;
      print " N_new $n_new[1]\n" ;
      print " N $n[$idx] \n" ;
      die "2nd Polynomial order did not correspond\n" ;
    }

    $idx = $matched[1] + $matched[2] - 1 ;
    if ( $n_new[2] != $n[$idx] ) {
      print " Matched $matched[1] $matched[2]\n" ;
      print " N_new $n_new[2]\n" ;
      print " N $n[$idx] \n" ;
      die "2nd Polynomial order did not correspond\n" ;
    }
    print "All atoms were matched\n" ;
  }
}

sub pair_index
{
  my ($ele0, $ele1) = @_ ;

  return $ele0 + $ele1 ;
}





