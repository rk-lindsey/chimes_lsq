This is a collection of reaxff force field input files for use with
lammps.  Here are some notes on the files:

Recommended:  ffield.reax.adri.wall.4.  Must be run with lmp_chaos5-intel.wall
binary, which is a modified version of LAMMPS.


Files based on 2012 HCNO force field:

ffield.reax.adri: This is the "new" HCNO force field file prepared by
Adri Van Duin for LLNL in 1/2012 under a contract with LLNL.

ffield.reax.adri.wall: This is an initial attempt at using an
exponential inner wall to prevent close distances between atoms from
occurring during MD simulations.  This file should be run only with
reax/c in the 10/2012 or later versions of LAMMPS.

ffield.reax.adri.wall.2: LEF added separate pairwise inner wall
parameters to more closely match N-H and O-H diatomic potential energy
surfaces.  This file should be run only with the LLNL-modified version
of lammps lmp_chaos5-intel.wall.

ffield.reax.adri.wall.3: This version uses Njo's 1998
electronegativity equilibration method parameters.  Calculations with
this model are under-charged.  There is likely a unit/input convention
error in using Njo's parameters with ReaxFF.

ffield.reax.adri.wall.4: This version has an increased EEM hardness
parameter for H and N.  The EEM hardness parameters were set to keep
the maximum charge on any atom in a shock simulation slightly less
than 1.  This potential is currently recommended for shock work.

Files based on 2005 HCNO force field.

ffield.reax.old: This is the HCNO force field used in Van Duin's 2005
paper on RDX.  This can be used with either the fortran version of
ReaxFF or the C version.  We found crashes in MSST and NVT simulations
using this potential.  The reason was excessively large atomic charges.

ffield.reax.gamma: This is a version of ffield.reax.old with modified
charge screening parameters (gamma).  This ran MSST and NVT
simulations stably, but had charge values that were lower than
expected for simple molecules such as H2O.  Wall parameters were also added.

ffield.reax.wall: This is a version of ffield.reax.old with wall
parameters added.  The parameters for every atom are the same.

ffield.reax.wall.2: This is a version of ffield.reax.old with wall
parameters added.  The parameters for every atom are based on QCISD(T)
dissociation calculations.

ffield.reax.wall.3: The charge screening parameters (gamma) were reduced.

ffield.reax.wall.4: Njo's 1998 charge equilibration MKS parameters
were used in ffield.reax.wall.2.  The system charges appeared to be
too small.

ffield.reax.wall.5: Njo's 1998 charge equilibraiton Mulliken (Table
10) parameters were used in ffield.reax.wall.2.  The system charges
appeared to be too small.










