#! /usr/bin/perl
#
# Translate a dftb format gen file into separate OpenMX DFT inputs.
# The spline_md code also creates a suitable gen file.
#
# One OpenMX input is created for each configuration in the gen file.
# These files can then be used to generate DFT forces from DFTB or force field 
# D trajectories.
#
# Usage:  gen_to_omx_inputs.pl <gen file> <openmx header> <openmx footer> <basename>
#
use strict ;
use warnings ;

die "Usage:  gen_to_omx_inputs.pl <gen file> <openmx header> <openmx footer> <basename>\n"
  unless ( scalar(@ARGV) == 4 ) ;

open(GEN,"<",$ARGV[0])    or die "Could not open $ARGV[0] for reading\n" ;
open(HEADER,"<",$ARGV[1]) or die "Could not open $ARGV[1] for reading\n" ;
open(FOOTER,"<",$ARGV[2]) or die "Could not open $ARGV[2] for reading\n" ;
my $basename = $ARGV[3] ;

my %valence = ( 
	        B => 3,
	        C => 4,
		N => 5,
		O => 6,
		F => 7,
		Ne => 8,
		H => 1,
	        Cl => 7
	      ) ;

my $header=join("", <HEADER>) ;
my $footer=join("", <FOOTER>) ;
close HEADER ;
close FOOTER ;

my $count = 1 ;
while ( 1 ) {
  
  my (@x, @y, @z, @type) ;

  $_ = <GEN> ;
  last unless defined($_) ;

  my @f = split(" ") ;
  my $natom = $f[0] ;

  $_ = <GEN> ;
  die "End of file reached prematurely\n" unless defined($_) ;

  my @elements = split(" ") ;
  for my $i ( 1 .. $natom ) {
    $_ = <GEN> ;
    die "End of file reached prematurely\n" unless defined($_) ;
    @f = split(" ") ;
    my $idx  = $f[0] ;
    $idx-- ;
    my $ele  = $f[1] ;
    $x[$idx] = $f[2] ;
    $y[$idx] = $f[3] ;
    $z[$idx] = $f[4] ;
    $type[$idx] = $elements[$ele-1] ;
    if ( ! defined($type[$idx]) ) {
      die "Type not defined\n" ;
    }
  }
  $_ = <GEN> ;
  my $avec = <GEN> ;
  my $bvec = <GEN> ;
  my $cvec = <GEN> ;
  chomp $avec ;
  chomp $bvec ;
  chomp $cvec ;
  
  die "End of file found prematurely\n" unless defined($cvec) ;

  my $filename = sprintf("%s.%03d.dat", $basename, $count) ;
  open(OMX,">",$filename) or die "Could not write to $filename\n" ;
  
  print OMX $header ;
  print OMX <<SEC1;
Atoms.Number $natom
Atoms.SpeciesAndCoordinates.Unit   Ang # Ang|AU
<Atoms.SpeciesAndCoordinates           # Unit=Ang.
SEC1

  for my $i ( 0 .. $natom - 1 ) {
    my $q = $valence{$type[$i]} ;
    die "Valence for $type[$i] not defined\n" unless defined($type[$i]) ;
    $q /= 2.0 ;
    printf OMX ("%3d %2s %21.13e %21.13e %21.13e %3.1f %3.1f\n",
	   $i+1, $type[$i], $x[$i], $y[$i], $z[$i],
	   $q, $q) ;
  }

print OMX<<SEC2;
Atoms.SpeciesAndCoordinates>
Atoms.UnitVectors.Unit Ang
<Atoms.UnitVectors
$avec
$bvec
$cvec  
Atoms.UnitVectors>
SEC2

  print OMX $footer ;
  close(OMX) ;
  $count++ ;
}

printf ("Created %d OpenMX inputs\n", $count-1) ;
