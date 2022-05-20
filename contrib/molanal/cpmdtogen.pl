#! /usr/bin/perl
# 
# A program to convert a CPMD trajectory to DFTB gen format.
# The gen format can then be read by other tools in this directory.
#
# Usage:  cpmdtogen.pl <TRAJECTORY> <cpmddat> <trajectory.gen>
# Modified by Nir Goldman, 1/23/07
#  allows for simulations such as MSSM where cell parameters may vary during the simulation
use Cwd;
use strict ;

{
    my ($ax, $ay, $az) ;
    my ($bx, $by, $bz) ;
    my ($cx, $cy, $cz) ;
    my @element ;
    my @count ;
    my ($j, $k, $l)  ;
    my $traj_file ;
    my $data_file ;
    my $gen_file ;
    my @atom_label ;
    my $natom ;
    my $step ;
    my ($x, $y, $z) ;
    my $cell_check;
    my $CELL;
    my $cell_step;
    my @line;
    my $prev_step;
    my $bohr = 0.5291772108 ;  # NIST VALUE


    if ( $#ARGV != 2 ) {
	die 'Usage:  cpmdtogen.pl <TRAJECTORY> <cpmddat> <trajectory.gen>' ;
    }
    $traj_file = $ARGV[0] ;
    $data_file = $ARGV[1] ;
    $gen_file  = $ARGV[2] ;

    open DATA, "< $data_file" || die "Could not open $data_file\n" ;
    open(TRAJ, "< $traj_file") || die "Could not open $traj_file\n" ;
    open(GEN, "> $gen_file") || die "Could not open $gen_file\n" ;

    $_ = <DATA> || die 'End of file in DATA file\n' ;
    ($ax, $ay, $az) = split(" ") ;
    $_ = <DATA> || die 'End of file in DATA file\n' ;
    ($bx, $by, $bz) = split(" ") ;
    $_ = <DATA> || die 'End of file in DATA file\n' ;
    ($cx, $cy, $cz) = split(" ") ;

    printf("Box vectors\n") ;
    printf("%7.3f %7.3f %7.3f\n", $ax, $ay, $az) ;
    printf("%7.3f %7.3f %7.3f\n", $bx, $by, $bz) ;
    printf("%7.3f %7.3f %7.3f\n", $cx, $cy, $cz) ;

    $CELL = "CELL";
    $cell_check = -e $CELL;
    if ($cell_check) {
      open (CELL, "<CELL") || die "DIE: can't open CELL\n";
    }

    $j = 0 ;
    while ( <DATA> ) {
	($element[$j], $count[$j]) = split(" ") ;
	$j++ ;
    }

    $l = 0 ; 
    for ( $j = 0 ; $j <= $#element ; $j++ ) {
	for ( $k = 0 ; $k < $count[$j] ; $k++ ) {
	    $atom_label[$l] = $j + 1  ;
	    $l++ ;
	}
    }
    $natom = $l ;
    
    printf("There are %d atoms\n", $natom) ;
    do {
	for ( $j = 0 ; $j < $natom ; $j++ ) {
	    if ( $_ = <TRAJ> ) {
		($step, $x, $y, $z) = split(" ") ;
                if ($. == 1) { # record first time step read in for checking purposes, below
                  $prev_step = $step;
                } 
                if ($. > $natom) { # ensure that $natom is set properly by checking value of time step in TRAJ file
                  if ($step != ($prev_step + 1)) {
                    print "Error in reading TRAJ file:\n";
                    print "step equals: $step\n";
                    printf ("     should equal: %d\n", $prev_step+1);
                    print "error could be in number of atoms in cpmd.dat file\n";
                    exit;
                  }
                }
		$x *= $bohr ;
		$y *= $bohr ;
		$z *= $bohr ;
	    } else {
		print "End of the trajectory file reached on step $step\n" ;
		exit(0) ;
	    }
	    if ( $j == 0 ) {
		printf GEN ("%5d S\n", $natom) ;
		for ( $k = 0 ; $k <= $#element ; $k++ ) {
		    print GEN '"', $element[$k], '" ' ;
		}
		print GEN "\n" ;
	    }
	    printf GEN ("%4d %4d %9.5f %9.5f %9.5f\n", $j+1, $atom_label[$j], $x, $y, $z) ;
	} 
        $prev_step = $step; # update record of previous time step
        if ($cell_check) {
          @line = split " ", <CELL>;
          $cell_step = $line[4];
          if ($cell_step != $step) {
            print "ERROR: cell_step != traj_step\n";
            print "  cell_step = $cell_step   traj_step = $step\n";
            exit;
          }
          $_ = <CELL> || die 'End of file in CELL file\n' ;
          ($ax, $ay, $az) = split(" ") ;
          $_ = <CELL> || die 'End of file in CELL file\n' ;
          ($bx, $by, $bz) = split(" ") ;
          $_ = <CELL> || die 'End of file in CELL file\n' ;
          ($cx, $cy, $cz) = split(" ") ;
          $ax *= $bohr;
          $ay *= $bohr;
          $az *= $bohr;
          $bx *= $bohr;
          $by *= $bohr;
          $bz *= $bohr;
          $cx *= $bohr;
          $cy *= $bohr;
          $cz *= $bohr;
        }
	printf GEN ("%9.5f %9.5f %9.5f\n", 0.0, 0.0, 0.0) ;
	printf GEN ("%9.5f %9.5f %9.5f\n", $ax, $ay, $az) ;
	printf GEN ("%9.5f %9.5f %9.5f\n", $bx, $by, $bz) ;
	printf GEN ("%9.5f %9.5f %9.5f\n", $cx, $cy, $cz) ;
    } while ( 1 ) ;
}
    
