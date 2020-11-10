Orthant-Wise Limited-memory Quasi-Newton algorithm minimizes functions of the form

  f(w) = loss(w) + C |w|_1
  
where loss is an arbitrary differentiable convex loss function, and |w|_1 is the L1 norm of the weight (parameter) vector.  It is based on the L-BFGS Quasi-Newton algorithm, with modifications to deal with the fact that the L1 norm is not differentiable.  The algorithm is very fast, and capable of scaling efficiently to problems with millions of parameters.

The algorithm is described in detail and proven to converge to the optimal parameter vector in
  Galen Andrew and Jianfeng Gao, "Scalable training of L1-regularized log-linear models", in ICML 2007.
Please cite the paper if you use OWL-QN in published research.

This package contains a c++ class "OWLQN" which you can use to optimize the parameters of any L1-regularized convex differentiable loss you specify, and a standalone main method for training L1-regularized logistic regression or least-squares models.  A compiled win32 executable for the trainer is also provided.

For questions, suggestions, bug reports, etc. please contact Galen Andrew at galena@microsoft.com.


USING THE STANDALONE OWL-QN TRAINER FOR L1-REGULARIZED LOGISTIC REGRESSION OR LEAST-SQUARES MODELS

Logistic regression formulation finds weights w that minimize
  sum_i log_loss(w | x_i, y_i) + C |w|_1

Least-squares formulation finds weights w that minimize
  sum_i 0.5 |<w, x_i> - y_i|^2_2 + C |w|_1

Note that in both cases, there is no "intercept term" b.  In very high dimensional spaces, explicitly adding an intercept is usually unnecessary.  If you need an intercept, you can add a "constant feature" that always has the value 1, or you can subtract the mean from your output in the case of least-squares.

usage: [options] feature_file label_file regWeight output_file
  feature_file   input feature matrix in Matrix Market format (mxn real coordinate or array).
                   rows represent features for each instance
  label_file     input instance labels in Matrix Market format (mx1 real array).
                   rows contain single real value
                   for logistic regression problems, value must be 1 or -1"
  regWeight      coefficient of l1 regularizer
  output_file    output weight vector in Matrix Market format (1xm real array).
  
options:
  -ls            use least squares formulation (logistic regression is default
  -q             quiet.  Suppress all output
  -tol <value>   sets convergence tolerance (default is 1e-4)
  -m <value>     sets L-BFGS memory parameter (default is 10)
  -l2weight <value>
                 sets L2 regularization weight (default is 0)


USING THE OWLQN CLASS FOR ARBITRARY CONVEX DIFFERENTIABLE L1-REGULARIZED LOSS FUNCTIONS

The c++ class "OWLQN" can be used to optimize the parameters of any L1-regularized convex differentiable loss.  (In fact, you could use it to optimize the weights of any convex differentiable loss without the L1 regularization by simply setting the regularization weight to zero.  For example, you could use L2 regularization instead by incorporating the value and gradient of the L2 regularizer in your objective, just as is done in the included implementations of logistic regression and least-squares.)

To call the optimizer, you need to define your loss function by deriving from the abstract base class DifferentiableFunction. This class has a single method
	double Eval(const vector<double>& input, vector<double>& gradient) const.
It must compute and return the value of the loss function at the point defined by the vector "input", and also fill in the vector "gradient" with the value of the gradient at that point.  You do not need to compute the value of the L1 regularization term, only the unregularized loss.

Then you simply instantiate OWLQN, and call its method Minimize:
	void Minimize(const DifferentiableFunction& function, const vector<double>& initial, vector<double>& minimum, double l1weight = 1.0, double tol = 1e-4, int m = 10) const;
The arguments are:
	function:  the unregularized loss function to miminize
	initial:   the initial point at which to begin optimization
	minimum:   a vector in which OWLQN will return the minimizing vector
	l1weight:  the weight "C" of the L1 regularizer
	tol:       convergence tolerance
	m:         memory parameter of L-BFGS
	
Notes:
	The default convergence criterion is the average decrease in function value over the last ten iterations, relative to the current value: ((v[10] - v[0]) / 10) / abs(v[0]), where v[i] is the function value from i iterations ago.  It is possible to specify any convergence criterion you like by implementing the TerminationCriterion interface, and passing an instantiation of your criterion to the OWLQN constructor.
	The "m" parameter of L-BFGS specifies how many previous iterations to store for use in approximating the function's inverse Hessian.  The algorithm stores 2m vectors of length equal to the dimensionality of the problem.  Using a higher "m" will probably decrease the number of iterations required to reach convergence, but increases memory requirements significantly and (less significantly) time requirements.  m = 0 (which is not allowed by OWLQN) would reduce L-BFGS to gradient descent.  Even m=1 is better than gradient descent, but the best value may depend on the problem and your time/memory requirements.  For most problems, the default value of 10 is reasonable.  Please see Nocedal and Wright's "Numerical Optimization, Second Edition" (2006, Springer) for more information about L-BFGS.
	OWLQN outputs to stdout information at each iteration.  "new_value" is the function value at the current point.  "conv_crit" is the value of the convergence criterion.  "line_search" displays a dot for each extra function evaluation performed during the line search (if any).  To turn off output, set OWLQN's "quiet" parameter to true.
	If _DEBUG compiler flag is turned on, OWLQN will also output the directional derivative as computed in two different ways: with a finite-difference approximation and using the reported gradient of the loss.  If these are systematically different, it may indicate a bug in the computation of the gradient of the loss.


RELEASE NOTES:

June 2007: 
	Original release.
Oct 2007:
	Added support for user-specified convergence criteria, including the duality gap for L1-regularized logistic regression and least-squares in the standalone trainer.
	Fixed three minor bugs (that affected output and efficiency only, not functionality)
	Compiles now with g++ and MSVS.


FUTURE WORK

We will soon release an updated version with the following improvements:
	Support for non-convex loss functions (in which case the algorithm will be guaranteed to find a local optimum, not necessarily global)
	
Improvements that may occur in the longer term include:
	Classification/Regression mode for the standalone app (currently it trains models, but does not perform classification)
	Computation of approximate regularization path with warmstart technique