Force matching codes by Lucas Koziol and Larry Fried

The supplied Makefiles will automatically compile the force fitting
program and run example calculations.  From the top-level directory,
simply type "make all".  Typing "make clean" will remove object files,
"make realclean" will restore the directory to only the original
source files.

Part A: Determining two-body parameters.

1.) Compile splines_ls using make in the src/lsq directory.

2.  Prepare input options for splines_ls in the file splines_ls.in.
Here is an example:
------------------------------------
nlayers 1
fit_coulomb true
fit_pover true
pair_type spline
subtract_coord false
smin 0.6 0.6 0.6
smax 6.0 6.0 6.0
sdelta 0.1 0.1 0.1
over_params 5 50.0 1.0165 -0.0657 5.0451 -3.6141
------------------------------------
For more examples see, example/h2o-lsq and example/lsq2.

nlayers is the number of x,y,z supercells used in evaluating
interactions.  Small unit cell should have nlayers >= 1.

fit_coulomb controls whether charges are determined by least
squares fitting.  Otherwise, charges are implicitly added to the
pair interaction but not otherwise considered.

fit_pover controls whether the ReaxFF linear overcoordination
parameter is found.

pair_type is the type of short-range pair interaction to use.
Choices are "spline" or "chebyshev".  "spline" computes the pair
interaction in terms of a cubic spline function.  "chebyshev" uses
polynomials of Morse functions to determine the interaction.

subtract_coord controls whether the ReaxFF overcoordination force
is subtracted from the forces before fitting.  This option should
not be used in conjunction with fit_pover.

smin is the minimum distance between atoms considered in the pair interaction.
     A parameter is given for each atom pair.

smax is the force cutoff distance used in the pair interaction.
     A parameter is given for each atom pair.

sdelta is the step size in angstroms between spline points.
     A parameter is given for each atom pair.

chebyshev_order is the number of terms in the Chebyshev force fitting.
It is not used in the example, because pair_type is set to spline.

3.) Code takes an input file called "input.xyzf". The format is:

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

An example input.xyzf file is located in example/h2o-lsq.

4.) Run the code spline fitting code
    ./splines_ls splines_md.in 99 

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

5.) Run the SVD script, python/lsq.py A.txt b.txt params.header >
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

The format of the params.txt file is as follows:

# Date  2016-06-14
# Number of columns =  34
# Number of rows    =  14114
# RMS force error =  13.9161502676
# max variable =  136.955400283
npair 3
3b_cheby false
pair_type chebyshev
 0.75000  6.00000  1.25000 10
 0.75000  6.00000  1.25000 10
 0.75000  6.00000  1.25000 10
coulomb true
fit_coulomb true
overcoord true
nover 5
  5.0000000000000e+01
  1.0165000000000e+00
 -6.5700000000000e-02
  5.0451000000000e+00
 -3.6141000000000e+00
fit_pover true
least squares parameters
0 62.2118312075
1 -23.0064319877
2 25.6115981525
3 -32.3001922969
4 2.13594798821
5 9.58482168479
6 36.259156255
7 -39.3514902706
8 -47.3118719056
9 -28.7304905567
10 -1.94218125182
11 28.490710697
12 26.7981298576
13 16.461114931
14 1.21081313796
15 4.56871613581
16 0.0147266482222
17 -0.808047108048
18 1.17821172921
19 0.875255025009
20 -10.7273078087
21 -17.4874987319
22 -35.0647688716
23 -36.6020682952
24 -36.9698808083
25 -27.8795218389
26 -20.534033962
27 -13.1578217039
28 -6.8693881773
29 -2.72199950168
30 136.955400283
31 -68.477702629
32 34.2388508942
33 25.2543658956

The params.txt file begins with comments.  The number of columns is the number
of fitting variables used in lsq.py programs.  The number of rows is the number
of force components fit to.  The number of columns should be much less than the number
of rows.  The RMS force error is given in units of kcal/mol-Ang.  The max variable
is the maximum fitting coefficient found.  If this value is very large, the
fit may be less reliable.

npair is the number of pair interactions.  The code has only been tested for water
(npair=3) at this time.  

3b_cheby controls the evaluation of 3-body chebyshev interactions.

The pair_type gives the pair interaction.  Possible values
are chebyshev or spline.  The following numbers are the minimum interaction distance,
maximum interaction distance, Morse distance scale, and maximum Chebyshev order.

coulomb controls the evaluation of coulomb interactions.  
If fit_coulomb is true, coulomb parameters are taken from the params.txt file.  Otherwise,
build-in charges are used.

If overcoord is true, a ReaxFF overcoordination interaction is calculated.
nover is the number of overcoordination parameters.  It should be 5 for ReaxFF.

fit_pover controls whether the specified first overcoordination parameter is used (false),
or whether the first overcoordination parameter is taken from the last listed parameter (true).

The least squares fitted parameters follow in order.  Pair parameters are given first,
then 3-body interaction parameters, then coulomb parameters, and finally the overcoordination parameter.

Part B: calculating short-range part of the force.

6.) If you plot the spline forces (plot every 2nd point in
params.txt), you will see that they equal nearly 0 up until a certain
(short) distance.  Then they will be oscillating for a few 0.1 Angs,
then they will become very well-defined. This is because the values of
A are average over interatomic distances sampled in the DFT-MD.

The spline force and potential can be conveniently plotted with the
code fmatch/bin/spline_vals.pl.  The negative of the force is placed
in spline.[0-2].txt, while the potential is placed in
splinepot.[0-2].txt.  spline.0.txt is the O-O interaction.
spline.1.txt is the O-H interaction, while spline.2.txt is the H-H
interactions.

The MD code splines_md supplements unoccupied values of the spline table with
a repulsive force to ensure that the MD simulation is stable.
Good energy conservation (4-5 digits) is achived with either the spline
table or with the polynomial force representation.

Part C: running Molecular dynamics.

7.) params.txt (with top two text lines removed) from output of lsq.py
should be in the working directory. 
All parameters specifying the force evaluation are placed in params.txt. This is 
different from earlier versions of the code, where params.txt had only least-squares
fitting parameters.

The filename of the force
parameter file can be used specified in splines_md.in. Also the
starting geometry and velocities should be included in input.xyz. The
N-3, N-2, and N-1 (N is the number of lines in params.txt) lines of
params.txt are the squared charges and must be entered as doubles qoo,
qoh, and qhh in functions.C.  This is now done automatically by the
fitting routine lsq.py.  The Nth line of params.txt is the magnitude
of the three body overcoordination interaction.

First line of input.xyz: number of atoms. 2nd line: 3 box 
dimensions in Angstroms.
3rd-number of atoms+3 line: label xcoord ycoord xvel yvel zvel.
Coordinates in Angs and velocities in Angs/MD-time unit 
(conversion from fs=0.0208025).

The MD program takes a configuration file given on the command line.  This file
controls various aspects of the MD calculation.  Here is an example
file:
---------------------------------------------------
# Input for the spline_md code.
temperature 2000.0
deltat      0.125
nsteps      100
nlayers     1
output_force false
read_force   false
init_vel     false
num_pressure false
params_file  params.cheby.txt
rand_seed     12357
hoover_time   10
energy_freq  10
scale_freq   0
gen_freq 20
----------------------------------------------------
The order of the options specified does not matter.  The code has
default values for all of the parameters.

The temperature is used for Nose-Hoover NVT dynamics, velocity scaling
dynamics, or velocity initialization.  

deltat is the MD time step in fs.  

nsteps is the number of MD steps.  

nlayers is the number of periodic replicas used in evaluating forces.  

output_force is a flag controlling whether forces are written to file.  

If read_force is true, the force.txt file created by lsq.py is read in
and compared to the current forces.  The input.xyzf file used in force
matching must be present.  It is read, and the first configuration is used
in the force comparison.

If spline_q is true, charges are taken from the params.txt file.
Otherwise compiled-in defaults are used for charges.

If init_vel is true, velocities will be initialized.  Otherwise,
initial velocities are read from the input.xyz file.  

pair_type specifies the type of short-range pair interaction.  Current
choices are chebyshev, spline, stillinger, or lennard-jones.

params_file specifies a file with force parameters.

rand_seed is a random number seed.  

hoover_time is a timescale in fs for the coupling between the
Nose-Hoover thermostat and the system.  If hoover_time is <= 0, no
Hoover thermostat is used.

scale_freq controls how frequently velocity is rescaled.  If
scale_freq is 0, velocities are not rescaled.  Velocity scaling and
the Hoover thermostat should not be used at the same time.

energy_freq controls how often the energy and other thermodynamic
quantities are output.

The trajectory is written to the file traj.gen, which is in DFTB gen format.
gen_freq controls how often configurations are written into a DFTB gen
file.

At the end of the run, positions and velocities are written to
output.xyz.  By copying output.xyz to input.xyz and running
splines_md, the trajectory may be continued further in time.



8.) compile splines_md using make

9.) Run the code by typing ./splines_md splines_md.in > splines_md.out.
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

6.  Previous versions of the MD code do not conserve energy.  The
current version should be used for all calculations.

7.  Support was put in for chebyshev polynomials of Morse functions.
Only the Chebyshev and spline potential options are well tested.


