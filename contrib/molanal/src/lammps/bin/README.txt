Larry Fried:  12/10/2012.

This directory has LAMMPS executables and perl scripts to help in 2012
ReaxFF simulations of TATB using LAMMPS.  Coordinate dump files are
assumed to be of the type q xu yu zu.

Files are as follows.  Perl files should be read to check for usage
and required editing before using.

dipole.pl:  Calculates system dipole moment from a LAMMPS dump file.

force_check.pl: Checks forces in a LAMMPS force dump file for
reasonable values.

lmp_chaos5-intel.wall:  This is a version of LAMMPS compiled with Intel
compilers for the Chaos 5 LLNL operating system.  It implements custom
pair-wise wall parameters.

pdb2: This is a program by Richard Gee to translate lammps dump files
into PDB files that can be read with AtomEye.

qmax.pl :  Calculates maximum and RMS charges from a LAMMPS dump file.
Can find charges for specified atom types, or values calculated across
all atom types.

qmax_check.pl:  Checks the LAMMPS dump file to make sure that the
charges are reasonable.

restart2data: This is a program supplied with LAMMPS to translate a
restart file into a data input file.

rmin.pl:  Calculates the minimum distance between atoms as a function of time, 
based on a LAMMPS radial distribution function (rdf) file.

supercellxyz.pl: Given a LAMMPS initial coordinate data file, this
will generate a supercell with a specified number of copies in the x,
y, and z directions.

vel_check.pl:  Checks velocities in a LAMMPS velocity dump file for reasonable
values.  Assumes vx vy vz format.

xyz2lammps.pl: This will translate a JMOL/XMOL xyz file into a lammps
data file.  It is set to work with the elements H/C/N/O.

