Force matching codes by Lucas Koziol and Larry Fried

Part A: Determining spline parameters.

1.) Compile splines_ls using make (modify Makefile in this directory).

2.) Code takes an input file called "input.xyzf". The format is:

Line 1: number of atoms.
Line 2: dimensions of box (Angs.)
Line 3 to natoms: Atomlabel xcoord ycoord zcoord xforce yforce zforce
Coords in Angs, Forces in atomic units (direct from OpenMX file).
Example:
//96
//8.6091309 8.6091309 8.6091309
//O 2.77746 1.36959 6.31692 0.001532305362 0.029788169404 -0.045470993055
//H 3.34223 2.07637 5.70246 -0.032494203931 -0.063952565017 0.041425795947
//H 4.36127 3.66403 7.07424 -0.042972843521 0.010603733297 0.000570193096
//...

3.) Run the code spline fitting code
    ./splines_ls 99 

99 here is the number of MD snapshots in the .xyzf file.
Once the code finishes (should take ~5 mins for 99 snapshots of
32-water cell), it will print two files "A.txt" and "b.txt". A.txt contains
derivatives of the force with respect to the linear parameters, while
b.txt contains the forces.  These are
used for input in a python script that does the SVD solution of Ax=b
and prints x. The x vector is a list of linear parameters.
There are two types of fitting that are done by splines_ls.
In spline fitting, grid points (two-body force values and 1st-derivatives 
on a grid) which is then read as input into the MD code.  In polynomial
fitting, polynomial coefficients for the two-body force are determined.
In should be noted that due to a historical accident, the parameter
values describe the negative of the pairwise force dV/dr.

4.) Run the SVD script, python/lsq.py A.txt b.txt params.header >
    params.txt params.header contains a header used in generating the
    spline or polynomial fit.  The params.header file is created by
    splines_ls.  For spline fitting, params.header contains the
    minimum spline distance, maximum spline distance, and spline step
    size in lines 1,2,3 respectively.  For chebyshev polynomial fitting,
    params.header contains the minimum spline distance, maximum spline
    distance, and number of polynomial coefficients in lines 1,2,3
    respectively.
    
    The lsq.py script also creates the file force.txt, which contains the 
    force on each atom as calculated with the best fit parameters.  This
    file can be used to test the consistency of the MD code with the reported
    forces.

5.) Delete the first two lines of params.txt, this just tells you the
RMS force error value and the maximum difference between fitted and DFT force.

Part B: calculating short-range part of the force.

6.) If you plot the forces (plot every 2nd point in params.txt), you will
see that they equal exactly 0 up until a certain (short) distance. 
Then they will be oscillating for a few 0.1 Angs, then they will become very 
well-defined. This is because the values of A are average over
interatomic distances sampled in the DFT-MD.

The spline force and potential can be conveniently plotted with the
code spline_vals.pl.  The negative of the force is placed in spline.[0-2].txt,
while the potential is placed in splinepot.[0-2].txt.  spline.0.txt is the
O-O interaction.  spline.1.txt is the O-H interaction, while spline.2.txt
is the H-H interactions.

The MD code splines_md supplements unoccupied values of the spline table with
a repulsive force to ensure that the MD simulation is stable.
Good energy conservation (4-5 digits) is achived with either the spline
table or with the polynomial force representation.

Part C: running Molecular dynamics.

7.) params.txt (with top two text lines removed) from output of
lsq.py should be in the working directory. Also the starting 
geometry and velocities should be included in input.xyz. The N-3, N-2,
and N-1 (N is the number of lines in params.txt) lines of params.txt are the squared charges and must be entered as doubles qoo, qoh, and qhh in functions.C.
This is now done automatically by the fitting routine lsq.py.  The Nth
line of params.txt is the magnitude of the three body overcoordination
interaction.

First line of input.xyz: number of atoms. 2nd line: 3 box 
dimensions in Angstroms.
3rd-number of atoms+3 line: label xcoord ycoord xvel yvel zvel.
Coordinates in Angs and velocities in Angs/MD-time unit 
(conversion from fs=0.0208025).

The MD program takes a configuration file splines_md.in.  This file
controls various aspects of the MD calculation.  Here is an example
file:
---------------------------------------------------
# Input for the spline_md code.
temperature 2000.0
deltat      0.125
nsteps      1000
nlayers     1
output_force false
read_force   false
spline_q     true
init_vel     false
chebyshev    false
rand_seed     12357
hoover_time   10
energy_freq  10
scale_freq   0
gen_freq 20
----------------------------------------------------

The temperature is used for Nose-Hoover NVT dynamics, velocity scaling
dynamics, or velocity initialization.  deltat is the MD time step in
fs.  nsteps is the number of MD steps.  nlayers is the number of
periodic replicas used in evaluating forces.  output_force is a flag
controlling whether forces are written to file.  In read_force is
true, the force.txt file created by lsq.py is read in and compared to
the current forces.  If spline_q is true, charges are taken from the
params.txt file.  Otherwise compiled-in defaults are used for charges.
If init_vel is true, velocities will be initialized.  Otherwise,
initial velocities are read from the input.xyz file.  If chebyshev is
true, chebyshev polynomials are used to evaluate the forces.
Otherwise, splines are used to evaluate the forces.  rand_seed is a
random number seed.  hoover_time is a timescale in fs for the coupling
between the Nose-Hoover thermostat and the system.  If hoover_time is
<= 0, no Hoover thermostat is used.  scale_freq controls how
frequently velocity is rescaled.  If scale_freq is 0, velocities are
not rescaled.  Velocity scaling and the Hoover thermostat should not
be used at the same time.  energy_freq controls how often the energy
and other thermodynamic quantities are output.  gen_freq controls how
often configurations are written into a DFTB gen file.

At the end of the run, positions and velocities are written to
output.xyz.  By copying output.xyz to input.xyz and running
splines_md, the trajectory may be continued further in time.

8.) compile splines_md using make

9.) Run the code by typing ./splines_md > splines_md.out.
The MD will print geometries and velocities at every step in
.xyz format, with MD info in the infoline. These snapshots can
be used as inputs for sequential run of code (put box 
dimensions in infoline instead of MD info). output.xyz is also
a movie of the trajectory that can be directly opened and 
visualized with molden.

11.) There is a fair amount of annotation in the program code
(splines_md.C and functions.C).

Work by Larry Fried: 2/24/2015

1.  Put a variable controlling the overcoordination subtraction in splines_ls.C.

2.  Changed lsq.py to output calculated forces.

Work by Larry Fried: 2/20/2015

The MD program in the mdcode directory was improved.  Related
functions in the main directory were not changed.  

1.  The Ewald sum calculation was greatly accelerated by
re-implementing the k space sum.

2.  The potential energy function for the force spline was
implemented.  MD calculations with the force spline conserve energy.

3.  Each term in the potential model was separated into a separate
routine.  Terms can be selected in ZCalc.

4.  Most terms in the potential were optimized to some degree.

5.  The behavior of "layers" for the force evaluation was changed.
The layers now always start from the minimum image displacement.  The
code should now work when nlayers = 0 (for a big system), but that hasn't
been tested.

6.  There are two test problems "input.small.xyz" and "input.big.xyz" in the 
mdcode directory.  You can select the test problem by copying these files to "input.xyz".

7.  Previous versions of the MD code do not conserve energy (try input.big.xyz).  The current version should be used for all calculations.




