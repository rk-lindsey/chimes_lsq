#! /usr/bin/perl
#
# Gather OpenMX forces from output files.
# Create an "xyzf" file suitable for use in splines_ls
#
use strict ;
use warnings ;

if ( scalar(@ARGV) != 2 ) {
  die "Usage:  gather_omx_force.pl <basename> <xyzf file>\n";
}

my $basename = $ARGV[0] ;

my @files = glob("$basename.*.out") ;

open(OUT,">", $ARGV[1]) or die "Could not write to $ARGV[1]\n" ;

my ($ax, $by, $cz) ;

foreach my $f ( @files ) {
  open(OMX,"<",$f) or die "Could not open $f\n" ;
  while ( <OMX> ) {

    if( /<Atoms.UnitVectors/ ) {
      $_ = <OMX> ;
      my @f = split(" ") ;
      $ax = $f[0] ;

      $_ = <OMX> ;
      @f = split(" ") ;
      $by = $f[1] ;

      $_ = <OMX> ;
      @f = split(" ") ;
      $cz = $f[2] ;
      chomp $cz ;
    }

    if ( /<coordinates.forces/ ) {
      my $natom = <OMX> ;

      print OUT $natom ;
      die "Lattice parameters not found\n" 
	unless ( defined($ax) and defined($by) and defined($cz) ) ;

      print OUT "$ax $by $cz\n" ;
      for my $j ( 1 .. $natom ) {
	$_ = <OMX> ;
	die "Premature end of file" unless defined($_) ;
	my ($idx, $ele, $x, $y, $z, $fx, $fy, $fz) = split(" ") ;
	printf OUT ("%2s %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n",
		    $ele, $x, $y, $z, $fx, $fy, $fz);
      }
    }
  }
  undef $ax ;
  undef $by ;
  undef $cz ;
}

