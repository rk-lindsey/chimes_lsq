NRFIT, a fitting program.
Larry Fried 10/14/2013

 nrfit is a driver for non-linear optimization of a general objective
 function.  Most of the algorithms are derived from Numerical Recipes
 in C.  The objective function is implemented through a system call.
 The system call can be a perl script, etc.  that executes a program
 and reports a value to be optimized.

 Usage:  nrfit -option option_file < fit.in >& fit.out
 


 Parameters that do not depend on the optimization algorithm are read
 from standard input (fit.in in the above example).

 Line 1 :  The name of the command that evaluates the objective function.
 Line 2 :  An input file which specifies parameters that the program in 
           line 1 will read.
 Line 3 :  An output file which specifies parameters that the program 
           in line 1 will write.  The output file should have a single 
           line, which is the value of the objective function.
 Line  4 : The number of parameters to optimize (N).
 Lines 5 : 4 + N : Initial values of the parameters.

 The options control which optimization methods are used.  A global
 optimization method may be followed by a local optimization method by
 specifying two options on the command line.

 The options are as follows:
  -p : Powell's local optimization method.  See Numerical Recipes.
  -d file : Davidson-Fletcher-Powell variable metric local optimization 
            method. See Numerical Recipes.
             File contents:
             Line 1:   Optimization function value tolerance.

  -b : The Downhill Simplex (amoeba) local optimization method.  
       See Numerical Recipes.
  -s file : Simulated annealing global optimization method.
             File contents:
             Line 1 :  Number of cycles
             Line 2 :  Initial temperature. Final temperature.
                       The initial temperature should be comparable
                       to the value of the function.  The final temperature
                       should be 2-3 orders of magnitude less than the
                       initial temperature.
             Line 3:3+N-1 :  Maximum MC step (1 line per parameter). 

  -a  file : Amoeba simulated annealing method.  See Numerical Recipes.
             File contents:  
             Line 1: Initial temperature.  Set to a value somewhat less
                     than a typical value of the objective function.
             Line 2: Number of blocks:  Set to 10 or greater.  The temperature
                       is reduced by a factor of 2 each block.
	       Line 3: Number of steps in a block.
  -r file : Record to record global optimization method.
            See G. Dueck, "New optimization hueristics: The great
            deluge algorithm and the record to record travel",
            J. Comp. Phys., 104, 86-92(1993).

             File contents:
             Line 1: Number of cycles to perform.  1 cycle is N trials,
                     where N is the number of variables to optimize.
             Line 2 : N + 1: Maximum MC step size.  1 line per variable.
             Line 2+N :  The deviation parameter.  This is an allowed
                         fractional tolerance from the best value 
                         obtained so far.


  -g file : Great deluge algorithm global optimization method.
            See G. Dueck, "New optimization hueristics: The great
            deluge algorithm and the record to record travel",
            J. Comp. Phys., 104, 86-92(1993).
            The algorithm as modified in "Comparison of the performance
            of modern heuristics for combinatorial optimization of
            real data", Computers Ops. Res., 20, 687-695(1995) is
            implemented.  The change in the water level is proportional
            to the decrease in the objective function.
            File contents:
             Line 1: Number of cycles to perform.  1 cycle is N trials,
                     where N is the number of variables to optimize.
             Line 2: Maximum MC step size.  1 line per variable.
             Line 2+N :  The water level change parameter.  Larger
                         values will converge more rapidly to 
                         a local minimum.  Smaller values will wander more.

Recommendations:  This is a re-working of a program used for many years to fit parameters
                  for the Cheetah thermochemical code.  The traditional use has been to
                  perform Powell optimization only, or to perform record to record optimization
                  followed by Powell optimization.  The simplex-type algorithms added recently
                  to the code are worth trying for problems with strongly coupled variables.
