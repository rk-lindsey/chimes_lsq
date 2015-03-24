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

Notes from Lucas Koziol:

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

3.) Run the code ./splines_ls 99 
99 here is the number of MD snapshots in the .xyzf file.

Once the code finishes (should take ~5 mins for 99 snapshots of 32-water cell),
it will print two files "A.txt" and "b.txt". These are used for input in 
a python script that does the SVD solution of Ax=b and prints x. The x vector
is a list of grid points (central-force values and 1st-derivatives on a grid)
which is then read as input into the MD code.

4.) Run the SVD script, python/lsq.py A.txt b.txt params.header > params.txt
    params.header contains a header used in generating the spline or polynomial
    fit.  The params.header file is created by splines_ls.  For spline fitting,
    params.header contains the minimum spline distance, maximum spline distance,
    and spline step size in lines 1,2,3 respectively.  For chebyshev fitting,
    params.header contains the minimum spline distance, maximum spline distance,
    and number of polynomial coefficients in lines 1,2,3 respectively.
    
    The lsq.py script also creates the file force.txt, which contains the force on
    each atom.

5.) Delete the first two lines of params.txt, this just tells you the 
X^2 value and the maximum difference between fitted and DFT force.


Part B: calculating short-range part of the force.

6.) If you plot the forces (plot every 2nd point in params.txt), you will
see that they equal exactly 0 up until a certain (short) distance. Then
they will be oscillating for a few 0.1 Angs, then they will become very 
well-defined. This is because the values of A are average over
interatomic distances sampled in the DFT-MD.

As a first step, you can bridge the spline part with an analytical 
Lennard-Jones with continuous values but not derivatives. This is already
coded; you will have to 
uncomment lines 363-365, 377-379, and 389-91. This will be enough to
make MD stable. The energy due to the splines is not included in this
version of the code. Integration of the spline force leads to overly
large energy fluctations (order of 0.01 units) compared to 10^-4
units with analytical potential/force. For production
you should fit the resulting spline function to some analytical 
form (see Koziol and Fried publication and Itzekov reference within).

Part C: running Molecular dynamics.

7.) params.txt (with top two text lines removed) from output of
lsq.py should be in the working directory. Also the starting 
geometry and velocities should be included in input.xyz. The last 
three lines of params.txt are the squared charges and must be
entered as doubles qoo, qoh, and qhh in functions.C 
First line of input.xyz: number of atoms. 2nd line: 3 box 
dimensions in Angstroms.
3rd-number of atoms+3 line: label xcoord ycoord xvel yvel zvel.
Coordinates in Angs and velocities in Angs/MD-time unit 
(conversion from fs=0.0208025).

8.) adjust temperature in first few lines of splines_md.C 

9.) compile splines_md using make

10.) Run the code by typing ./splines_md > output.xyz
The MD will print geometries and velocities at every step in
.xyz format, with MD info in the infoline. These snapshots can
be used as inputs for sequential run of code (put box 
dimensions in infoline instead of MD info). output.xyz is also
a movie of the trajectory that can be directly opened and 
visualized with molden.

11.) There is a fair amount of annotation in the program code
(splines_md.C and functions.C) including commented sections on how
to assign random initial velocities and remove velocity c.o.m.





