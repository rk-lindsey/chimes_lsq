
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
	
	bool solve_G_A_con_grad()
	// Use the conjugate gradient method to find G_A_Inv_I and A_A
	// Does not use cholesky decomposition.
	// Returns true if the solution passes consistency tests.
		{

			const double eps_fail = 1.0e-04 ;  // Max allowed solution error.
			const double eps_con_grad = 1.0e-08 ;  // Conjugate gradient error tolerance.

			// Solve for G_A^-1 * unity
			Vector unity(nactive, 1.0) ;

			// if ( RANK == 0 ) {
			// 	// cout << "G_A_Inv_I guess = \n" ;
			// 	G_A_Inv_I.print(cout) ;
			// }
			bool use_ssor = false ;
			bool use_chol_precon = true ;

			if ( use_precondition ) {
				if ( use_ssor ) {
					// SSOR preconditioner

					pre_con.realloc(nactive, nactive) ;			
					Matrix K(nactive, nactive) ;
					double omega = 1.1 ;
					double omega_scale = sqrt(2.0-omega) ;
			
					for ( int i = 0 ; i < nactive ; i++ ) {
						for ( int j = 0 ; j <= i ; j++ ) {
							if ( i != j ) {
								double val = -omega_scale * sqrt(omega/G_A.get(i,i))
									* omega * G_A.get(i,j) / G_A.get(j,j) ;
								K.set(i,j,val) ;
							} else {
								double val = omega_scale * sqrt(omega/G_A.get(i,i))
									* (1.0-omega) ;
								K.set(i,i,val) ;
							}
						}
						for ( int j = i+1 ; j < nactive ; j++ ) {
							K.set(i,j,0.0) ;
						}
					}
					// pre_con = K^T * K
					for ( int i = 0 ; i < nactive ; i++ ) {
						for ( int j = 0 ; j < nactive ; j++ ) {
							double sum = 0.0 ;
							for ( int k = 0 ; k < nactive ; k++ ) {
								sum += K.get(k,i) * K.get(k,j) ;
							}
							// Approximate pre-conditioner.
							pre_con.set(i,j,sum) ;

							// Use diagonal matrix.
							//pre_con.set(i,j,0.0) ;
						}
						// Use diagonal matrix.
						// pre_con.set(i,i,1.0/G_A.get(i,i)) ;
					}

				} else if ( use_chol_precon ) {
					// Test: use the inverse of G_A to precondition.
					pre_con.resize(nactive, nactive) ;

					if ( nactive == A_last.dim + 1 && nactive > 2 ) {
						// Add 1 row to pre-conditioner.
						pre_con.set(nactive-1, nactive-1, 1.0) ;

						// DEBUG !
						//if ( RANK == 0 ) {
						//cout << "Updated preconditioner\n" ;
						//pre_con.print() ;
					    //}
						
					} else {
						// Full calculation of pre-conditioner.
						Matrix chol_precon(nactive, nactive) ;
						
						if ( ! G_A.cholesky(chol_precon) ) {
							if ( RANK == 0 ) cout << "Cholesky decomposition for pre-conditioning failed\n" ;
						}
					
						chol_precon.cholesky_invert(pre_con) ;
					}
				}
					
				if ( ! G_A.pre_con_grad(G_A_Inv_I, unity, pre_con, nactive+10, 10, eps_con_grad) ) {
					if ( RANK == 0 ) cout << "Pre-conditioned conjugate gradient failed\n" ;
					return false ;
				}
			} else { // ! use_precondition
			
				if ( ! G_A.con_grad(G_A_Inv_I, unity, nactive+10, 10, eps_con_grad) ) {
					if ( RANK == 0 ) cout << "Conjugate gradient failed\n" ;
					return false ; 
				} 
			}
			
			// if ( RANK == 0 ) {
			// 	cout << "G_A_Inv_I solution" << endl ;
			// 	G_A_Inv_I.print(cout) ;
			// }

			// Test to see if the solution worked.
			Vector test(nactive) ;
			G_A.dot(test, G_A_Inv_I) ;
			double errval = 0.0 ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				errval += fabs(test.get(j)-1.0) ;
				if ( fabs(test.get(j) - 1.0) > eps_fail ) {
					if ( RANK == 0 ) {
						cout << "Conjugate gradient solution test failed\n" ;
						cout << "Error = " << fabs(test.get(j) - 1.0) << endl ;
					}
					return false ;
				}
			}
			if ( nactive > 0 && RANK == 0 ) cout << "Cholesky error test = " << errval / nactive << endl ;
			
			
			A_A = 0.0 ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				A_A += G_A_Inv_I.get(j) ;
			}
			if ( A_A > 0.0 ) 
				A_A = 1.0 / sqrt(A_A) ;
			else {
				if ( RANK == 0 ) cout << "A_A Normalization failed" << endl ;
				return false ;
			}
			return true ;
		}

	void build_u_A()
		{
			const double eps_fail = 1.0e-04 ;
			w_A.realloc(nactive) ;
			u_A.realloc(ndata) ;
			a.realloc(nprops) ;

			G_A_Inv_I.scale(w_A, A_A) ;
#ifdef VERBOSE			
			cout << "w_A " << endl ;
			w_A.print() ;
#endif			
			
			X_A.dot(u_A,w_A) ;

#ifdef VERBOSE			
			cout << "U_A " << endl ;
			u_A.print() ;
#endif			
			
			double test = 0.0 ;
			for ( int j = 0 ; j < ndata ; j++ ) {
				test += u_A.get(j) * u_A.get(j) ;
			}
			test = sqrt(test) ;
			if ( fabs(test-1.0) > eps_fail ) {
				cout << "U_A norm test failed" << endl ;
				cout << "Norm = " << test << endl ;
				stop_run(1) ;
			}

			// Test X_A^T u__A = A_A * I
			Vector testv(nactive,0.0) ;
			X_A.dot_transpose(testv, u_A) ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				if ( fabs(testv.get(j) - A_A) > eps_fail ) {
					cout << "u_A test failed " << endl ;
					stop_run(1) ;
				}
			}
				
			X.dot_transpose(a, u_A) ;

#ifdef VERBOSE			
			cout << "a vector = " << endl ;
			a.print() ;
#endif			

		}
	
	void reduce_active_set() 
	// Reduce the active set of directions to those having maximum correlation.
	// See Eq. 3.6
		{
			// Undo the change in the active set.
			if ( RANK == 0 ) cout << "Will remove property " << A.get(remove_prop) << " from the active set" << endl ;
			A.remove(remove_prop) ;
			nactive = A.dim ;
			gamma_lasso = 1.0e20 ;
		}

	void update_active_set() 
	// Update the active set of directions to those having maximum correlation.
		{
			IntVector a_trial(nprops) ;
			//int count = 0 ;
			const double eps = 1.0e-6 ;

			// Save the last active set
			A_last.realloc(nactive) ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				A_last.set( j, A.get(j) ) ;
			}

			if ( do_lasso && gamma > gamma_lasso ) {
				reduce_active_set() ;
			} else if ( add_prop >= 0 ) {
				if ( RANK == 0 ) cout << "Adding property " << add_prop << " to the active set" << endl ;
				A.push(add_prop) ;
				nactive++ ;
			} else {
				// Either we are restarting or something strange has happened.
				// Search for the new active set.
				int count = A_last.dim ;
				a_trial = A_last ;

				for ( int j = 0 ; j < nprops ; j++ ) {
					if ( fabs( fabs(c.get(j)) - C_max ) < eps
							 && ! exclude.get(j) ) {
						int k ;
						// See if this index has occurred before.
						for ( k = 0 ; k < nactive ; k++ ) {
							if ( j == A_last.get(k) ) {
								break ;
							}
						}
						if ( k == nactive ) {
							a_trial.push(j) ;
							count++ ;
							if ( RANK == 0 ) cout << "Adding property " << j << " to the active set" << endl ;
							// Break to add only one property to the active set at a time.
							break ;
						}
					}
				}
				A.realloc(count) ;
				nactive = count ;
				for ( int j = 0 ; j < nactive ; j++ ) {
					A.set(j, a_trial.get(j)) ;
				}
			}
			if ( RANK == 0 ) cout << "New active set: " << endl ;
			A.print_all(cout) ;

#ifdef USE_MPI
			// Sync the active set to avoid possible divergence between processes.
			MPI_Bcast(&nactive, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
			MPI_Bcast(A.vec, nactive, MPI_INT, 0, MPI_COMM_WORLD) ;
#endif			
		}

	void update_step_gamma()
	// Update gamma, eq. 2.13.
		{
			double huge = 1.0e20 ;
			gamma = huge ;

			remove_prop = -1 ;
			add_prop = -1 ;
			
			if ( nactive < nprops ) {
				for ( int j = 0 ; j < nprops ; j++ ) {
					int k = 0 ;
					for ( k = 0 ; k < nactive ; k++ ) {
						if ( A.get(k) == j ) 
							break ;
					}
					if ( k != nactive ) continue ;
					double c1 = ( C_max - c.get(j) ) / (A_A - a.get(j) ) ;
					double c2 = ( C_max + c.get(j) ) / (A_A + a.get(j) ) ;

					if ( c1 > 0.0 && c1 < gamma ) {
						gamma = c1 ;
						add_prop = j ;
					}
					if ( c2 > 0.0 && c2 < gamma ) {
						gamma = c2 ;
						add_prop = j ;
					}
				}
			} else {
				// Active set = all variables.
				gamma = C_max / A_A ;
			}
			if ( RANK == 0 ) {
				cout << "Updated step gamma = " << gamma << endl ;
				if ( add_prop >= 0 )
					cout << "Gamma limited by property " << add_prop << endl ;
			}
			if ( do_lasso ) update_lasso_gamma() ;
		}

	void update_lasso_gamma()
	// Find the Lasso gamma step, which may be less than the LARS gamma step.
	// See Eq. 3.4 and 3.5
		{
			gamma_lasso = 1.0e20 ;
			const double eps = 1.0e-12 ;
			
			for ( int i = 0 ; i < nactive ; i++ ) {
				if ( fabs(w_A.get(i)) > 1.0e-40 ) {
					double gamma_i = -beta.get(A.get(i)) / ( sign.get(i) * w_A.get(i) ) ;
					if ( gamma_i > eps && gamma_i < gamma_lasso ) {
						gamma_lasso = gamma_i ;
						remove_prop = i ;
					}
				}
			}
			if ( RANK == 0 ) cout << "Lasso step gamma limit = " << gamma_lasso << endl ;
		}

	void update_beta()
	// Update the regression coefficients (beta)
		{
			if ( do_lasso && gamma > gamma_lasso ) {
				if ( RANK == 0 ) {
					cout << "LASSO is limiting gamma from " << gamma << " to " << gamma_lasso << endl ; 
					cout << "LASSO will set property " << A.get(remove_prop) << " to 0.0" << endl ;
				}

				//cout << "Current beta:" << endl ;
				//beta.print() ;
				//cout << "Current correlation: " << endl ;
				//c.print() ;

				gamma_use = gamma_lasso ;
			} else {
				gamma_use = gamma ;
			}
			for ( int j = 0 ; j < nactive ; j++ ) {
				int idx = A.get(j) ;
				double val = beta.get(idx) + w_A.get(j) * sign.get(j) * gamma_use ;
				beta.set(idx,val) ;
			}
			// If we are removing a property, set the value to exactly 0.
			// Check that the calculated value is close to 0.
			if ( do_lasso && gamma > gamma_lasso ) {
				if ( fabs(beta.get(A.get(remove_prop))) > 1.0e-08 ) {
					if ( RANK == 0 ) {
						cout << "Error: failed to set variable to zero when removing prop\n" ;
						stop_run(1) ;
					}
				}
				beta.set(A.get(remove_prop),0.0) ;
			}
			
#ifdef USE_MPI
			// Sync the beta values to avoid possible divergence between processes.
			MPI_Bcast(beta.vec, nprops, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
#endif
			
#ifdef VERBOSE			
			cout << "New beta: " ;
			beta.print(cout) ;

			cout << "Predicted mu: " << endl ;
			for ( int j = 0 ; j < ndata ; j++ ) {
				cout << mu.get(j) + gamma_use * u_A.get(j) << " " ;
			}
			cout << "\n" ;
#endif
			
		}

	void print_unscaled(ostream &out) 
	// Print the coefficients in unscaled units.
		{
			if ( RANK == 0 ) {
				double offset = y.shift ;
				Vector uns_beta(nprops,0.0) ;
				for ( int j = 0 ; j < nprops ; j++ ) {
					if ( X.scale[j] == 0.0 ) {
						out << "Error: scale factor = 0.0" << endl ;
						stop_run(1) ;
					}
					offset -= beta.get(j) * X.shift[j] / X.scale[j] ;
				}
				for ( int j = 0 ; j < nprops ; j++ ) {
					uns_beta.set(j, beta.get(j) / X.scale[j]) ;
				}
				if ( out.rdbuf() == cout.rdbuf() ) {
					uns_beta.print(cout) ;
				} else {
					for ( int j = 0 ; j < nprops ; j++ ) {
						out << uns_beta.get(j) << endl ;
					}
				}
			}
		}
	
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

	void broadcast_solution()
	// Broadcast results of solving G_A.
	{
#ifdef USE_MPI
		if ( RANK != 0 ) {
			G_A_Inv_I.realloc(nactive) ;
			A.realloc(nactive) ;
		}
		
		MPI_Bcast(&nactive, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
		MPI_Bcast(A.vec, nactive, MPI_INT, 0, MPI_COMM_WORLD) ;
		MPI_Bcast(G_A_Inv_I.vec, nactive, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
		MPI_Bcast(&A_A, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
#endif
	}


		
};
