#! /usr/bin/perl
# Run LAMMPS for diatomics wall calculation.
# Larry Fried 11/13/2012.
# 
# Usage:  diatomic.pl
# 
# The program will automatically calculate potential curves for the specified diatomics.
# Edit rname to control the force field files tested.
# Edit distances to control the distances tested.
# Edit molecules to change the molecules that are tested.

use strict ;
use warnings ;
use File::Copy ;

#####
##### USER-DEFINED PARAMETERS BELOW
#####

# Set to the desired LAMMPS executable.
my $lammps="/g/g17/fried/lammps/bin/lmp_chaos5-intel.wall" ;

# Choose interatomic distances.
my @distances = (0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65,
		 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 
		 1.15, 1.2, 1.25, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 
		 1.9, 2.0) ;

#my @distances = (1.0) ;

# Specify force field files and descriptive names.
my %rname = ( 
#              "ffield.reax.old" => "Old",
	      "ffield.reax.adri" => "New",
	      "ffield.reax.adri.wall" => "New Wall",
	      "ffield.reax.adri.wall.2" => "New Wall 2",
#              "ffield.reax.wall.3" => "Wall 3"
#	      "ffield.reax.wall.4" => "Wall 4"
) ;

## Define ordering of atom types for various ffield.reax files.
my %reax_type = ( 
    "C" => 1, 
    "H" => 2,
    "O" => 3,
    "N" => 4) ;

my %reax_type2 = ( 
    "C" => 1, 
    "H" => 2,
    "O" => 3,
    "N" => 7) ;


## Define which ordering is used in which ffield.reax file.
my %reax_type_hash = (
    "ffield.reax.old" => \%reax_type,
    "ffield.reax.adri" => \%reax_type2,
    "ffield.reax.adri.wall" => \%reax_type2,
    "ffield.reax.adri.wall.2" => \%reax_type2,
    "ffield.reax.wall.2" => \%reax_type,
    "ffield.reax.wall.3" => \%reax_type,
    "ffield.reax.wall.4" => \%reax_type 
) ;

# List of molecules to test.
my @molecules = ("CC", "HH", "OO", "NN", "CO", "CN", "CH", "NO", "NH", "OH") ;
#my @molecules = ("NH", "OH") ;
#my @molecules = ("NN") ;


####
#### END OF USER-DEFINED PARAMETERS
####

my @reax_files = keys %rname ;
my @reax_names = sort values %rname ;

#my @reax_files = ("ffield.reax.wall.2", "ffield.reax.old", "ffield.reax.adri", 
#		  "ffield.reax.wall.3") ;

#my @distances = (0.3) ;


my @energies ;

# Hash of references to energies.
# my %eref ; (NOT USED YET)

# Convert from Hartree to kcal/mol.
my $kcal = 627.5 ;

foreach my $atoms ( @molecules ) {
    my (@args, $atom1, $atom2, $id, $idu) ;

    &set_names($atoms, \@args, \$atom1, \$atom2, \$id, \$idu) ;

    printf "Testing molecule $id\n" ;

    foreach my $reax ( @reax_files ) {
	# Loop over force field files.

	&make_lammps_files($reax, \%reax_type_hash, \@args) ;

	@energies = &run_lammps(\@distances,\@args) ;
#	$eref{$reax} = \@energies ;
	my $ffield_name = $rname{$reax} ;
	print "Force field name = $ffield_name file = $reax\n" ;
	open(OUT,">pot.$id.${ffield_name}.out") || die "Could not open potential file\n" ;
	for ( my $j = 0 ; $j <= $#distances ; $j++ ) {
	    printf( OUT "%13.4e %13.4e\n", $distances[$j], $energies[$j]) ;
	}
    }
    &run_gnuplot($id,$idu) ;
}

sub run_gnuplot
## Use GNUPLOT to make a plot of the potential energy curves.
{
    my $id = $_[0] ;
    my $idu = $_[1] ;

    open(ABINIT,"<qcisdt/${id}\_abinitio.out" ) 
	|| die "Could not open ${id}\_abinitio.out\n" ;
    open(AOUT, ">abinit.out") 
	|| die "Could not open abinit.out\n" ;
    $_ = <ABINIT> ;
    $_ = <ABINIT> ;

    my ($r,$eoffset, $e) ;
    ($r, $eoffset) = split(" ") ;
    while ( <ABINIT> ) {
	($r, $e) = split(" ") ;
	if ( defined($r) && $r > 0.0 ) {
	    $e -= $eoffset ;
	    $e *= $kcal ;
	    printf (AOUT "%13.4e %13.4e\n", $r, $e) ;
	}
    }
    close (ABINIT) ;
    close (AOUT) ;
#    || die "Could not open ${id}\_abinitio.out\n" ;

    open(GNU, ">diatomic.gnu") ;
    print GNU <<EOF;
set term post color eps enhanced 20 linewidth 2
set encoding iso_8859_1
set size 0.7,1.2
set out "diatomic.$id.ps"
set xlabel "R({\305})"
set ylabel "E (kcal/mol)"
set title "$idu potential energy"
set yrange [:200]
set xrange [0.3:2.0]
set xtics 0.5
set mxtics 5
set ytics 100
set mytics 5
EOF

	print GNU "plot \"abinit.out\" u 1:2 w lp lw 3 ti \"Ab initio\" " ; 
    my $j = 1 ;
    foreach my $name ( @reax_names ) {
	print GNU ", \"pot.$id.$name.out\" u 1:2 w lp ti \"$name\"" ;
	$j++ ;
    }
    print GNU "\nquit\n" ;

    close(GNU) ;
    system("gnuplot < diatomic.gnu") == 0
	|| die "GNUplot returned an error\n" ;
}


sub run_lammps 
# Run lammps for the given distances, return an array of energies.
{
    my @energies ;

    my @distances = @{$_[0]} ;
    my @args      = @{$_[1]} ;

    foreach my $rx ( @distances ) {
	my $energy ;
	if ( scalar(@args) == 1 ) {
	    open(IN,"<diatomic.template") || die "Could not open xyz template\n" ;
	} else {
	    open(IN,"<diatomic2.template") || die "Could not open xyz template\n" ;
	}
	open(OUT,">diatomic.data") ;
	while ( <IN> ) {
	    if ( /<RX>/ ) {
		$_ =~ s/<RX>/$rx/ ;
	    }
	    print OUT ;
	}
	close(IN) ;
	close(OUT) ;
	system("srun -n 1 $lammps < in.diatomic > out.diatomic 2> err.diatomic") == 0
	    || die "Lammps failed\n" ;

	open(IN, "<out.diatomic") || die "Could not open out.diatomic\n" ;
	while ( <IN> ) {
	    if ( /PotEng/ ) {
		$_ = <IN> ;
		($energy) = split(" ") ;
#		print "ENERGY = $energy\n" ;
		push(@energies, $energy) ;
	    }
	}
	close(IN) ;
    }
    return @energies ;
}

sub set_names
# Given atoms in the molecule, set argument strings and id strings.
# Arguments 2-6 are passed by reference.

# Usage:     &set_names($atoms, \@args, \$atom1, \$atom2, \$id, \$idu) ;
#
{
    my $atoms = $_[0] ;
    my $args  = $_[1] ;
    my $atom1 = $_[2] ;
    my $atom2 = $_[3] ;
    my $id    = $_[4] ;
    my $idu   = $_[5] ;

    ($$atom1, $$atom2) = split(//, $atoms) ;
    if ( $$atom1 eq $$atom2 ) {
	$$args[0] = $$atom1 ;
    } else {
	$$args[0] = $$atom1 ;
	$$args[1] = $$atom2 ;
    }

    if ( scalar(@$args) == 1 ) {
	$$id = lc($$args[0]) . "2" ;
    } else {
	$$id = $$args[0] . $$args[1] ;
	$$id = lc($$id) ;
    }
    $$idu = uc($$id) ;
    $$idu =~ s/2/_{2}/ ;

}

sub make_lammps_files
# Give the atom types and a hash containing reaxff atom indices, write input files for LAMMPS.

{
    my $reax = $_[0] ;
    my $reax_type_hash = $_[1] ;
    my $args = $_[2] ;

    my @atom_type ;
    my $reax_order ;

    foreach my $type ( @$args ) {
	# Obtain atom order appropriate to the force field file.
	$reax_order = $$reax_type_hash{$reax} ;
	if ( defined($reax_order->{$type}) ) {
	    push(@atom_type, $reax_order->{$type}) ;
	} else {
	    die "Order for Atom type $type is not known\n" ;
	}
    }

    open (LMP,"<in.template") || die "Could not open Lammps input template\n" ;
    open (LMPOUT,">in.diatomic") || die "Could not open lammps input file\n" ;
    while ( <LMP> ) {
	if ( /<TYPES>/ ) {
	    chop ;
	    $_ =~ s/<TYPES>// ;
	    print LMPOUT ;
	    foreach my $type ( @atom_type ) {
		print LMPOUT (" $type") ;
	    }
	    print LMPOUT "\n" ;
	} else {
	    print LMPOUT ;
	}
    }
    copy($reax,"ffield.reax") ;
}
