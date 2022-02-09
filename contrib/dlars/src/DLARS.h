
class DLARS
{
public:
	Matrix &X ;         // The matrix of data vs properties
	Matrix pre_con ;   // A pre-conditioning matrix for conjugate gradient.
	Vector y ;         // The vector of data
	Vector mu ;        // The predicted values of data
	Vector beta ;      // The scaled fitting coefficients
	Vector beta_A ;    // The fitting coefficients for the current active set.
	Vector c ;         // The correlation vector.
	IntVector A ;      // Indices of the active set of properties.
	IntVector A_last ;  // Indices of the last active set of properties.
	IntVector exclude ; // Indices of properties not to use because of degeneracy.
	int num_exclude ;   // The number of excluded properties
	IntVector sign ;   // Signs for each property in X_A
	Matrix G_A ;       // Active set Gram matrix.
	Matrix X_A ;       // Matrix formed by of X for active properties * signs
	Matrix chol ;       // Cholesky decomposition of G_A.
	Vector G_A_Inv_I ;  // G_A^-1 * I
	double A_A ;        // Normalization (I G_A^-1 I)^-1/2
	Vector u_A ;        // Unit vector for predicted values.  Step occurs along this director.
	Vector w_A ;        // Step direction vector for fitting coefficients
	double C_max ;      // Maximum correlation found on this iteration.
	double lambda ;     // Weighting for L1 norm in objective function.
	Vector a ;          // Update vector for correlation.
	double gamma ;      // LARS Step size ;
	double gamma_lasso ;  // Lasso constraint on LARS step size.
	double gamma_use ;  // The value of gamma to use on current step.  Based on gamma and gamma_lasso.
	int ndata ;       // Number of data items to fit to = X.dim1
	int nprops ;      // Number of properties to correlate = X.dim2
	int nactive ;     // The number of active properties to correlate <= nprops
	int remove_prop ;  // The index of the property to remove from the active set during a LASSO calculation.
	int add_prop ;  // The index of the property to add to the active set during a LASSO calculation.
	bool do_lasso ;   // If TRUE, do a lasso calculation.  If false, do a regular LARS calculation.
	bool solve_succeeded ; // If true, the solve of G_A succeeded.
	bool solve_con_grad ;  // If true, solve G_A by conjugate gradient instead of cholesky decomp.
	bool use_precondition  ; // If true, use preconditioning in conjugate gradient.
	bool distributed_solver ;  // If true, use a distributed G_A matrix and cholesky solver.
	double obj_func_val ;  // Latest value of the objective function.
	int iterations ;    // The number of solver iterations.
	ofstream trajfile ;  // Output file for the trajectory (solution history).

  DLARS(Matrix &Xin, Vector &yin, double lamin): X(Xin), y(Xin.dim1), mu(Xin.dim1),
		beta(Xin.dim2), c(Xin.dim2), A(0),  exclude(Xin.dim2, 0), sign(0), lambda(lamin)
		{
			do_lasso = false ;
			gamma_lasso = 1.0e20 ;
			gamma = 0.0 ;
			gamma_use = 0.0 ;
			ndata = X.dim1 ;
			nprops = X.dim2 ;
			nactive = 0 ;
			remove_prop = -1 ;
			add_prop = -1 ;
			num_exclude = 0 ;
			solve_succeeded = true ;
			solve_con_grad = false ;
			use_precondition = false ;
			
			iterations = 0 ;

			X_A.distribute(Xin) ;

			if ( RANK == 0 ) {
				trajfile.open("traj.txt") ;
				trajfile.precision(12) ;
				trajfile << scientific ;
			}
			for ( int j = 0 ; j < ndata ; j++ ) {
				y.set(j, yin.get(j)) ;
				mu.set(j, 0.0) ;
			}
			for ( int k = 0 ; k < nprops; k++ ) {
				beta.set(k, 0.0) ;
				c.set(k, 0.0) ;
			}
			y.shift = yin.shift ;
			for ( int k = 0 ; k < nprops ; k++ ) {
				X.scale[k] = Xin.scale[k] ;
				X.shift[k] = Xin.shift[k] ;
			}
		}

	int iteration() ;
	void build_G_A(Matrix &G_A_in, bool increment_G_A)	;
	bool build_G_A_here() ;
	void decrement_G_A(Matrix &G_A_in) ;
	void increment_G_A(Matrix &G_A_in)	;
	bool solve_G_A(bool use_incremental_updates) ;
	bool chol_backsub(Matrix &G_A_in, Matrix &chol_in) ;
	void increment_G_A_dist(Matrix &G_A_dist)	;
	void decrement_G_A_dist(Matrix &G_A_dist) ;	
	void predict() ;
	void predict_all() ;
	double sq_error() ;
	void objective_func() ;
	void correlation() ;
	void build_X_A() ;
	int restart(string filename) ;
	void broadcast_solution() ;
	bool solve_G_A_con_grad() ;
	void build_u_A() ;
	void reduce_active_set() ;
	void update_active_set() ;
	void update_step_gamma() ;
	void update_lasso_gamma() ;
	void update_beta() ;
	void print_unscaled(ostream &out)  ;
	
	void print_unshifted_mu(ostream &out)
	// Print the given prediction in unscaled units.
		{
			if ( RANK == 0 ) {
				//out << "Y constant offset = " << offset << endl ;
				for ( int j = 0 ; j < ndata ; j++ ) {
					out << mu.get(j) + y.shift << endl ;
				}
			}
		}
	
  void print_unshifted_mu(ostream &out, Vector &weights)
	// Print the given prediction in unscaled units.
	{
		if ( RANK == 0 ) {
			//out << "Y constant offset = " << offset << endl ;
			for ( int j = 0 ; j < ndata ; j++ ) {
				out << (mu.get(j) + y.shift)/weights.get(j) << endl ;
			}
		}
	}	

	void print_error(ostream &out)
	// Print the current fitting error and related parameters.
	{
		if ( RANK == 0 ) {
			out  << "L1 norm of solution: " << beta.l1norm() << " RMS Error: " << sqrt(sq_error() / ndata) << " Objective fn: " << obj_func_val << " Number of vars: " << A.dim << endl ;
		}
	}

	void print_restart()
		// Print the restart file
	{
		ofstream rst("restart.txt") ;

		if ( RANK == 0 && rst.is_open() ) {
			rst << scientific ;
			rst.precision(16) ;
			rst.width(24) ;
			rst << "Iteration " << iterations << endl ;
			print_error(rst) ;
			beta.print_sparse(rst) ;
			rst << "Exclude " << endl ;
			exclude.print_sparse(rst) ;
			rst << "Mu" << endl ;
			mu.print_sparse(rst) ;
		}
		rst.close() ;
	}



		
};
