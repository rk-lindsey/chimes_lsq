Tools for the analysis of reactive molecular simulations (molanal)

Larry Fried  11/10/06-06/13/13

These are molecular analysis (molanal) tools.  They are used as follows:

1.  molanal.new - 
  
Usage:  molanal.new <trajectory file> <charge file> > molanal.out

This file takes as input the trajectory and charge (CHR.DAT) files
created by dftb_llnl.  The original DFTB does not accumulate charge
data.  The trajectory format is very simple, so it should be easy to
write a converter from another program.  The program will loop over
each frame of the trajectory, and report which molecules are present
in the simulation.  An xyz output file (molanal.xyz) is created in
which the molecules are "re-constructed" from the raw trajectory.  The
charge file is optional.  If charges are given to the program, the
dipole moment wrt the center of mass and the net charge are reported
for each molecule.

As of 04/2012, molanal.new can be run in parallel with OPENMP threads.
Before execution of the program, enter the command "export
OMP_NUM_THREADS=4", or "setenv OMP_NUM_THREADS 4", depending on
whether you are running in the bash or tcsh shells, respectively.
Testing shows that threads can give a modest (factor of 2 or less)
speed-up on 8000 atoms.  Speed-ups may be better for larger
systems. If you are running in the batch system, you must request more
than 1 thread.  A simple way to do this is to run a 1 node job on a
production parallel system, such as atlas.llnl.gov.

As of 04/2012, molanal.new creates an output file called
"molanal.com".  This file contains the center of mass coordinate of
each molecule (or non-bonded atom) found in the simulation.
Coordinates are labeled with the name of the molecule, according to
the molanal.new naming convention.

molanal.new works for a non-orthorhombic unit cell.  Since the code
loops over a fixed number of periodic images, it could fail to give
correct results for a severely distorted cell.  For such problems, the
MAXWRAP parameter in the code should be increased. MAXWRAP is not used
for an orthorhombic unit cell.  For an orthorhombic unit cell, the
minimum image convention is used to find the nearest image of a given
atom.  For an orthorhombic unit cell, the x unit vector should be
specified first in the gen file, then the y unit vector, and finally
the z unit vector.

As of 04/2012, molanal.new has a test suite.  After compiling molanal.new,
move to the tests directory, and type "make all".  Tests of differing size
are compared to reference output.

As of 06/2013, molanal.new supports skipping steps for long MD simulations.
The file "skip.dat" must be present to generate the step skipping.  If
the file is not present, steps are not skipped.  Frames numbers reported
by molanal.new used to be shifted by the longest bond lifetime.  This
was to make the frame number start at 1.  The code has been changed so
that the frames are no longer shifted.  The frame number printed is the
same as the MD frame number.  Times reported by findmolecules.pl will
likewise be shifted from their former values.  This was done to make
the bookkeeping of time more simple and accurate.

bonds.dat:
Bond distances used to find molecules are taken from the file
bonds.dat.  As of 11/10/06, the bonds.dat file takes additional
input parameters.

First line: Minimum bond lifetime in MD frames (dump time intervals).
Two atoms will not be considered bonded unless they are continuously
within the cutoff distance for at least the given number of frames.  To
get the behavior of the previous code version, enter 0 here.  A
minimum bond lifetime of 10-30 fs seems reasonable.

New feature (NG, 11/07): In order to use different lifetimes for different 
bond types, set the minimum bond lifetime to be negative.  
Note: this feature requires the use of the file bond_times.dat (see below).

Second line: Number of copies for xyz file.  The specified number of
periodic copies of the simulation cell in each direction will be
printed to the xyz file.  This means that if you have natoms in the
cell, there will be natoms * (2*copies+1)^3 atoms output each frame to
the xyz file.  Adding copies helps to make the trajectory look better.

Third line:  The time between frames in the gen file in seconds.

Fourth line: If 1 is given, a directory called "molecules" will be
created containing a separate xyz file for each molecule found in the
simulation.  As of 06/18/2013, the molecules directory will not be
deleted if it already exists.  This allows multiple parallel
invocations of molanal.new to share the same molecules directory.
The directory should be deleted by hand before running molanal.new
in most situations.

The file naming is similar to that used in the molanal.out
output file.  If 0 is given, no per-molecule xyz files will be generated.
The xyz files will reflect the last frame containing a particular molecule.

Additional lines:  Bond distances. Here is the format:
H H 1.1
H O 1.2

etc.

Here is a sample bonds.dat file:
-- Begin bonds.dat
20
1
2.41e-15
1
H H 0.900
C H 1.400
N H 1.450
O H 1.350
C C 1.900
C N 1.800
C O 1.800
N N 1.750
N O 1.650
O O 1.700
-- End bonds.dat

bond_times.dat: (optional)
Use this file to set different bond cut times for different atom pairs (see
above comment on bonds.dat).  The format is similar to the last part of 
bonds.dat. The bond lifetimes are given in terms of the MD dump (frame) intervals:
H H 15
H O 41

etc.
In this example, 15 means that the H-H bond lifetime is equal to 15 dump intervals.
Given bonds.dat, this would correspond to 2.41e-15 * 15 seconds.

skip.dat:

This file controls the skipping of steps to analyze.  This can substantially speed
up the analysis of long MD simulations.  The stepping is controlled by two skip 
frequencies, a block length, and an offset.  

The first parameter is called "skip1" (step skip).
Every skip1 th step is a candidate for analysis.  The bond lifetimes in bonds.dat
are automatically changed to take account of the skip1.  If the skip1 is larger than
the bond lifetime, then an instantaneous bond criterion will be used.

The second parameter is called "skip2" (block skip).  Skip2 skips between blocks of MD steps.
The length of the block is specified as the third parameter.  For instance, if
the block length is 100 and skip2 is 5, then steps 1 to 100 will be analyzed.
Then the code will skip to analyze steps 501 to 600.  Next it analyzes 1001 to 1100, etc.

The fourth parameter is called the offset.  The offset controls which blocks to analyze.
For instance, in the above example, if the offset is 1, then steps 201 to 300 will be analyzed,
followed by 601 to 700, etc.

Here is an example skip.dat:
-----------------------
2  # step skip
10 # block skip
20 # block length
0  # block offset
------------------------
2.  gr_calc

Usage: gr_calc < gr_calc.in > gr_calc.out

This program calculates the radial distribution function for given element
pairs.  The program is written in fortran 90, so you will need an f90 compiler.
g95 or gfortran are freely available.

Here is an example gr_calc.in file:

-- Begin gr_calc.in
tatb_113cj5.out  ! Trajectory file.
4        ! input format type (4 = DFTB)
4 4      ! Desired atom types.
0.05     ! dr
2        ! number of super-cells added in each direction. (use 2 for TATB).
100.0    ! time blocking	
1        ! number of time blocks.
0.00025  ! dump time interval.
10       ! Diffusion time block.
-- End gr_calc.in

The first line is the trajectory file.  The second line is the input
format type.  ONLY DFTB (4) HAS BEEN TESTED.  The third line gives the
desired atom types.  So if the DFTB trajectory file lists

"C" "N" "O" "H"

Then 1 is "C", 2 is "N", 3 is "O" and 4 is "H". The example file
calculates the H-H radial distribution file.  The dr parameter is
the bin size for the radial distribution function.  The gr_calc program
works with non-orthorhombic unit cells.  For such cells, a super-cell
search is done to determine distances.  Convergence with respect to 
the number of super-cells should be checked.  1 is adequate for 
orthorhombic cells.  

The time blocking parameter gives the option of calculating G(r) for
multiple time slices of the trajectory.  Give a large number if you
just want one time slice.  The number of time blocks desired is then
given.  The time interval parameter is the time between configurations in
the trajectory file.  This is used in calculating the time blocks.

The diffusion time block is used in calculating the diffusion
coefficient and mean squared displacement.  The displacement of each
particle is placed into a bin with a width given (in steps) by the
diffusion time block.  Also, every nth time step is used for a new
initial condition in calculating the displacement, where n is also the
diffusion time block.  Diffusion constants are calculated by the
formula D = 1/6 d/dt < |r(t) - r(t0)| >^2.  The trajectory is divided
into 5 segments.  Each segment is used to calculate the time
derivative in the formula through finite differences.  This leads to
some scatter in the estimates, which is realistic given the difficulty
in calculating accurate diffusion constants.  The 1st estimate of the
diffusion constant starts from 0 time and thus includes the "inertial"
part of the system response.  It should not be used in reported
diffusion constants.

One reason for calculating the g(r) is to estimate a bond distance.
Two estimates are automatically printed out:  the second value where
g(r) = 1, and the first minimum of g(r) where g(r) < 1.  The second
was used by Goldman et al. in studying water.  The first was suggested
in a recent PRL on high temperature water.  Note that in some cases
one or both of the criteria can fail to give an answer.  In that case
there may not be a covalent bond between the species involved, or the
fraction of bonds may be very low.

3.  findmolecules.pl 

Usage: findmolecules.pl molanal.out > findmolecules.out

This program takes the output of molanal.new, and reports the 
average concentration, lifetime, dipole moment, and charge of each species.
The program has been modified as of 11/10/06 to calculate more quantities.
The code now also defines a "polymer" composite species.  This is useful
for simulations where a great number of polymeric species are formed.
All these species are now lumped into a single "polymer".  

There are a number of user-defined inputs inside the findmolecules.pl program.
As of 04/2012, these inputs may be specified in a separate file 
"findmolecules.cfg", which should be located in the working directory
of the computation.  Alternatively, the findmolecules.cfg file may
be omitted, and the may inputs specified by editing the beginning of the findmolecules.pl
program instead.

Here is the section of the program to edit:
##################################
#  USER-DEFINED VARIABLES
#
# The block_size is the number of frames to average concentrations over.
$block_size = 100 ;
# Molecules over this size will be counted as a "polymer".
$max_molecule_size = 10 ;
# Stop reading after this many frames. Or when end of file is found.
$stop_frame = 100000 ;
# Do not print out molecules with concentration less than this (in units of mol/cc).
$conc_floor = 1.0e-04 ;
# Minimum lifetime for a molecule (in ps).
$min_mol_life = 0.02 ;
# Number of iterations to use in finding reactions.
$react_iter_max = 10 ;
# Set to 1 in order to analyze rates, 0 otherwise.
$analyze_rates = 1 ;
# Reaction histories will be printed only for reactions with more than 
# this number of reaction events.
$reaction_count_floor = 10 ;
# If set to 1, calculate reverse reaction pairs.  This is somewhat slow right now.
$find_reverse_rxn = 0 ;
#  END USER-DEFINED VARIABLES
##################################

The block_size determines how many saved MD steps (frames) are used in
averaging output concentration histories.

Some tuning of max_molecule_size is desirable.  Any molecule with more
then max_molecule_size is counted as a "polymer".  This prevents an
explosion of complexity in simulations where very many polymeric species
are formed.  If max_molecule_size is set to a large number (greater than
the number of atoms in the simulation), no "polymer" species will be created
and all molecules will be tracked.

The trajectory file will only be read up to the given stop_frame.  This is
useful for short exploratory calculations while determining parameters.

Only molecules with concentration larger than conc_floor are printed out.

A bonded entity is counted as a "molecule" if is lasts longer than the
number of frames given by min_mol_life.  Otherwise, it is counted as a
"transition state".

If analyze_rates is set to 1, reaction rates will be found as well as
reaction histories.  If it is set to 0, no reaction rate analysis will be
performed.

findmolecules.pl prints the each frame it reads into findmolecules.log.
It is somewhat helpful to look at this file if the program takes a long
time to run.  The program run time is sensitive to the max_molecule_size
for some simulations.

findmolecules.pl has been updated as of 11/10/06 to print out 
reactions as well as species.  Reactions are identified based on
the loss or gain of molecules between two frames.  A new file
findmolecules.log, is created when the program is run.  This
file has a list of which reactions occur during a particular frame.

If find_reverse_rxn is set to 1, reverse reactions are found and the net
reactive flux for the forward and reverse reactions are determined.  The 
reactive flux for a reaction i is defined as natoms * |nforward - nreverse|,
where natoms is the number of atoms on one side of the reaction.  nforward
is the number of forward reactions, and nreverse is the number of reverse
reactions.

4.  molanal.pl  

Usage: molanal.pl <gen name> <chr name>. 

This program runs molanal.new on a set of numbered trajectory and
charge files.  Since dftb_llnl no longer creates numbered files, this
program is mostly obsolete.  THIS HAS NOT BEEN TESTED IN A LONG TIME.

5.  gentoarc.pl

Usage: gentoarc.pl < <gen file>

This program reads in dftb trajectory (gen) file, and creates a
corresponding MSI arc file.  This file can be read in by Materials
Studio or by Cerius2.


6. molanal

This is the "classic" version of the molecular analyzer.  It only
works for orthorhombic unit cells.  It does not analyze charges.  Bond
lengths are compiled in.  Currently it is faster than molanal.new.
THIS HAS NOT BEEN TESTED IN A LONG TIME.  molanal.new supports bond
lifetimes, which gives a much better bond definition.

7. cpmdtogen.pl

This program converts a CPMD trajectory to DFTB gen format. The 
trajectory must be contiguous (no overlapping restart data).
The gen format can then be read by other tools in this directory.

Usage:  cpmdtogen.pl <TRAJECTORY> <cpmddat> <trajectory.gen>

The TRAJECTORY file is the cpmd trajectory dump.  The cpmddat
file contains extra data used to create the gen file.  Here
is an example:

- begin cpmd.dat -
7.9273 0.0 0.0
0.0 7.9273 0.0
0.0 0.0 7.9273
F 54
H 54
- end cpmd.dat -

first three lines contain the unit cell vectors.  They do not have to be
orthogonal.  The next lines contain the element types and counts used
in the trajectory file.  In the example, there are 54 HF molecules.  The
first 54 atoms are fluorines and the next 54 atoms are hydrogens.

The file trajectory.gen is the DFTB gen-format file created as output.

8.  lammpsxytogen.pl

This program, written by Nir Goldman, converts a LAMMPS xyz dump file
into the DFTB gen format.  Please be careful not to use this program
on a LAMMPS fractional coordinate dump file.  The atom_types variable
in the lammpsxyztogen.pl program file must be set to be consistent
with the LAMMPS atom types.  The order of the atom types is
significant.

Usage:  lammpsxyztogen.pl <dump file>

The output is written with the same basename as the dump file, but
with the file extension ".gen" added.

9.  lammpsfractogen.pl

This program, written by Larry Fried, converts a LAMMPS fractional
dump file into the DFTB gen format.  This program has not been tested
recently.  Please be careful not to use this program on a LAMMPS xyz
coordinate dump file.  The atom types are specified in the lammps.dat
file.

Usage:  lammpsfractogen.pl <TRAJECTORY> <lammps.dat> <trajectory.gen>

<lammps.dat> is a file containing the name of each element in order
according to the # lammps dump file in the first line.

10.  openmx/wrap_openmx.pl

This program will the xyz coordinates of an OPENMX input file into the
primitive periodic cell.  OPENMX 3.6 does not wrap the initial
coordinates of atoms into the cell.  This can lead to some atoms being
left off the computational grid, with strange results.  Inputs should
be wrapped if the initial grid is specified with the Atoms.UnitVectors
command.  This script works only if the box is orthogonal, with the
three unit vectors pointing in the x, y, and z directions,
respectively.  A supercell of the original system can be created by
specifying sx, sy, and sz.

Usage:  ./wrap_openmx.pl sx sy sz input.old > input.new

sx, sy, and sz are the number of supercells in the x, y, and z
directions.  They are optional arguments.

11.  openmx/omx_md_to_gen.pl

This program will convert the "*.md" file produced by the openmx 3.6
program into a DFTB gen file.  The gen file can then be read by the 
molanal.new program.  This program is written for an orthorhombic
simulation cell.  The dimensions of the cell are given as command
line arguments.

Usage: omx_md_to_gen.pl file.md boxx boxy boxz

The output is written to "file.gen".  Boxx, boxy, and boxz are the
dimensions of the simulation cell.
     

12.  ReaxFF force fields for high density simulations.

Riad Manaa and Larry Fried have produced ReaxFF force fields for high 
density simulations.  See notes in reax/README.txt for information on
the force field files.

13.  LAMMPS MD analysis tools produced for TATB simulations.  Some of the
tools may be useful for more general purpose MD simulations as well.
See the lammps/bin/README.txt file.  There are also scripts to compare
LAMMPS potentials to quantum simulations.  See the lammps/diatomic/README.txt 
file.


13.  Parallel Analysis Driver (pad)

This program reads a list of commands, and executes the commands
across a number of MPI processes.

Usage: srun -n <nprocs> pad <infile>

infile is an input file with a list of shell commands to execute.
The first line of inline gives the number of commands per processor,
which is fixed.  Usually this is the number of sequential steps in
the molecular analysis.
  
The next lines give commands to be executed in parallel.  It is
assumed that each analysis program is MPI serial, but it may be
thread parallel.  If the analysis program is thread parallel, an
appropriate number of CPU's per mpi process must be allocated
in the msub script.

Here is an example input file:
---------------------------------
2
molanal.new Benzene.md.gen > molanal.1.out
findmolecules.pl molanal.1.out > findmolecules.1.out
molanal.new Benzene.md.gen > molanal.2.out
findmolecules.pl molanal.2.out > findmolecules.2.out
molanal.new Benzene.md.gen > molanal.3.out
findmolecules.pl molanal.3.out > findmolecules.3.out
---------------------------------

This could be executed in parallel by up to 3 mpi processes.  Each
process will execute 2 commands in sequence.

14. gencommand.pl

This program generates a list of commands which are executed in
parallel by the pad program.  gencommand.pl is a simple program
that is meant to be modified on a case by case basis, calling
various analysis tools.

Usage: gencommand.pl > commands.pad
       srun -n procs pad commands.pad



15. fixgen.pl

Usage: fixgen.pl <gen files>

In some cases dump files may overlap in time if the dump interval and
restart interval are different in the simulation.  This script fixes a
number of gen files so they start sequentially.  It relies on a time
comment in the dump file, which should say # time = xxxx, or #TIMESTEP
yyyy.  The script searches for the word time, case insensitive, and
extracts the next number.

A list of all files to be fixed needs to be given in sequential order
(e.g. dump.001, dump.002, etc.).  The first dump file (dump.001) is
truncated if need be so that it ends before the second (dump.002)
begins.  The original version of fixed dump files are saved with a
".old" prefix.  It is strongly recommended that you save backup copies
of the gen files before running this script to prevent possible data
loss.

The simplest way to use this program is to simply give it all the dump
files and let it fix them.  Unfortunately, this is very slow for a
large MD run with many gen files.  It is possible to fix the gen files
in parallel with the pad program (see above), but there are some
caveats.  fixgen.pl needs at least two file arguments.  The first one
will be the file to fix, and the second one is the next gen file to
read for comparison.  In order to avoid file conflicts, the files must
be fixed in two batches: first the odd-numbered files are used as
first arguments, while the even-numbered files are used as second
arguments.  In the second batch, the even-numbered files are used as first
arguments, while the odd-numbered files are used as second arguments.

Here is an example pad command file. This can be run with up to 5 processors
using pad.

_______________
1
/g/g17/fried/bin/fixgen.pl tatb.msst.9kpms.001.x.gen tatb.msst.9kpms.002.x.gen
/g/g17/fried/bin/fixgen.pl tatb.msst.9kpms.003.x.gen tatb.msst.9kpms.004.x.gen
/g/g17/fried/bin/fixgen.pl tatb.msst.9kpms.005.x.gen tatb.msst.9kpms.006.x.gen
/g/g17/fried/bin/fixgen.pl tatb.msst.9kpms.007.x.gen tatb.msst.9kpms.008.x.gen
/g/g17/fried/bin/fixgen.pl tatb.msst.9kpms.009.x.gen tatb.msst.9kpms.010.x.gen
_______________

16. splitmol.pl

This program splits the output of "joinmolecules.pl" into separate
output files for plotting. Output files will be labeled by molecule
and/or property, and will end with the suffix .txt.  It hasn't been
tested yet for the output of findmolecules.pl.  In the meantime,
you can run a single findmolecules.pl output through joinmolecules.pl,
and then split it using splitmol.pl

Usage: splitmol.pl <join file>

17. grabatom.pl

This program will extract specified atoms from a gen file, and write them
into a new gen file.  This is useful if you want to track only a subset of the
atoms in a simulation.  The atom numbers are used to identify which atoms in the 
gen file to extract.  There is no limit on the number of atoms that can be
extracted (although the command line may become very long).

Usage:  grabatom.pl <atom numbers> <gen file> > <new gen file>

