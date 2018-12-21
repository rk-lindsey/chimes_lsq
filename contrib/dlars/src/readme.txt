dlars is a distributed memory, parallel implementation of the LARS
algorithm for LARS regression or LASSO regression.  It is designed to
optimize 1/2N || Ax - b ||^2_2 + lambda || x ||_1, where N is the
number of data entries (rows) in the A matrix.

Currently, storage of the A matrix is distributed across processes, and operations involving
the A matrix or the active subset of the A matrix are parallelized.  This allows very large
A matrices to be treated.  The A matrix size can be greater than the single-node memory
if more than one node is used.

The LARS algorithms solves linear equations for an active set of variables.  At the present
time this is not parallel, and the storage of the active set matrices is local.  The current
code should handle problems up to about 10,000 active (non-zero) variables, with no limit on
the number of non-active variables and the number of data entries.

The notation and implementation closely follows B. Efron, T. Hastie,
I. Johnstone, and R. Tibshirani, "Least Angle Regression", The Annals
of Statistics, 32, 407-499(2004).  The X and X_A matrices in the paper are
distributed, while other matrices are local.

Usage:
	srun -n <XX> dlars <A matrix> <b vector> <A dim> <options> > dlars.log
	
Inputs:
   A matrix:        Matrix of properties.
	b vector:        Vector of data values.
   A dim            File with dimensions of A.  The number of columns is given first, followed
	                 by the number of rows.
Outputs:
   Diagnostic standard output is given for each iteration.
	x.txt  :  The final set of fit parameters.
	Ax.txt :  The final estimated value of b (= Ax), given the parameters x.
	traj.txt: The 'trajectory' of x values for each iteration, along with values of
	          the RMS error and objective function.

Options:
   --algorithm=<alg>    The algorithm may be either LARS or LASSO.  LASSO always gives better answers than
	                     LARS for a given L1 norm, but requires more iterations.
   --iterations=<num>   Sets a limit on the number of iterations allowed.
	--lambda=<val>       Sets the value of lambda in the objective function.  By default, lambda = 0.0
	--max_norm=<val>     Set the maximum L1 norm of the solution.  This is based on the scaled variables.
	--normalize=<y or n> Specifies whether the A matrix and b vector are normalized prior to fitting.
	                     The default is to normalize.
	--split_files        If specified, split input files are read.  Instead of A.txt, A.0000.txt,
	                     A.0001.txt, etc. is read by each MPI process.  This can speed job execution
								for large A matrices.  The chimes_lsq code generates these files if the
								#SPLITFI# option is specified.
								A corresponding dim file is required for each A.xxxx.txt file, e.g. dim.xxxx.txt.
								The dim file gives:  The number of data columns, the starting row for the file,
								the ending row (inclusive) for the file, and the total number of rows in the A matrix.
								
   --help               Print a list of supported options.
