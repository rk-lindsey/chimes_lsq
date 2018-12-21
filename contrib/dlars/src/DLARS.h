
class DLARS
{
public:
	Matrix X ;         // The matrix of data vs properties
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
	double obj_func_val ;  // Latest value of the objective function.
	
  DLARS(Matrix &Xin, Vector &yin, double lamin): X(Xin.dim1, Xin.dim2,true), y(Xin.dim1), mu(Xin.dim1),
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

			X.distribute() ;
			X_A.distribute() ;
			
			for ( int j = X.row_start ; j <= X.row_end ; j++ ) {
				for ( int k = 0 ; k < nprops; k++ ) {
					X.set(j, k, Xin.get(j,k) ) ;
				}
				for ( int j = 0 ; j < ndata ; j++ ) {
					y.set(j, yin.get(j)) ;
					mu.set(j, 0.0) ;
				}
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

	void predict() 
	// Calculated predicted values of y (mu hat, Eq. 1.2).  Update previous prediction
	// based on u_A and gamma_use.
		{
			if ( mu.size() != ndata ) {
				cout << "Error:  matrix dim mismatch" << endl ;
				exit(1) ;
			}

			if ( u_A.size() == 0 ) {
				// First iteration.
				for ( int j = 0 ; j < ndata ; j++ ) {
					mu.set(j,  0.0) ;
				}
			} else {
				for ( int j = 0 ; j < ndata ; j++ ) {
					mu.set(j,  mu.get(j) + gamma_use * u_A.get(j) ) ;
				}
			}
/**			
			for ( int j = 0 ; j < ndata ; j++ ) {
				double tmp = 0.0 ;
				for ( int k = 0 ; k < nprops ; k++ ) {
					tmp += X.get(j,k) * beta.get(k) ;
				}
				mu.set(j, tmp) ;
			}
**/
#ifdef VERBOSE			
			cout << "Mu = " << endl ;
			mu.print() ;
#endif
		}
	
	void predict_all() 
	// Calculated predicted values of y (mu hat, Eq. 1.2) with no updating based on u_A.
		{
			if ( mu.size() != ndata ) {
				cout << "Error:  matrix dim mismatch" << endl ;
				exit(1) ;
			}
			X.dot(mu, beta) ;

#ifdef VERBOSE			
			cout << "Mu = " << endl ;
			mu.print() ;
#endif
		}
	double sq_error() 
	// Squared error Eq. 1.3
		{
			double result = 0.0 ;
			for ( int j = 0 ; j < ndata ; j++ ) {
				result += (y.get(j)-mu.get(j)) * (y.get(j) - mu.get(j) ) ;
			}
			return(result) ;
		}

	void objective_func()
	// Calculate optimization objective function based on the requested
	// regularization parameter lambda.  This should be called after
	// predict_all() or predict().
	{
		obj_func_val = 0.5 * sq_error() / ndata + lambda * beta.l1norm() ;
	}
	void correlation() 
	// Calculate the correlation vector c, Eq. 2.1
		{
			C_max = -1.0 ;

			if ( gamma_use <= 0.0 ) {
				// First iteration.
				Vector ydiff(ndata,0.0) ;
				for ( int k = 0 ; k < ndata ; k++ ) {
					ydiff.set(k, y.get(k) - mu.get(k)) ;
				}

				X.dot_transpose(c, ydiff) ;
			} else {
				// c = c - gamma_use * a.
				c.add_mult(a, -gamma_use) ;
			}
					
			for ( int j = 0 ; j < nprops ; j++ ) {
				if ( fabs(c.get(j)) > C_max ) {
					// Only look for C_max if the coordinate has not been excluded.
					int i = 0 ;
					for ( ; i < exclude.dim ; i++ ) {
						if ( j == exclude.get(i) )
							break ;
					}
					if ( i == exclude.dim )
						C_max = fabs(c.get(j)) ;
				}
			}

			//cout << "New correlation: " << endl ;
			//c.print() ;
			//cout << "Max correlation:" << C_max << endl ;
		}
	void build_X_A() {
		// Calculate the sign and the X_A array.
		X_A.realloc(ndata, nactive) ;
		sign.realloc(nactive) ;

		// Calculate the sign of the correlations.
		for ( int j = 0 ; j < nactive ; j++ ) {
			if ( c.get( A.get(j) ) < 0 ) 
				sign.set( j, -1) ;
			else
				sign.set( j, 1) ;
		}

		// Calculate the X_A array.
		for ( int j = X_A.row_start ; j <= X_A.row_end ; j++ ) {
			for ( int k = 0 ; k < nactive ; k++ ) {
				double val = X.get( j, A.get(k) ) * sign.get(k) ;
				X_A.set(j,k, val) ;
			}
		}
	}

	void build_G_A()
	// Build the G_A matrix after X_A has been built.
	// Increment the matrix from previous versions if possible.
		{
			if ( A.dim == A_last.dim + 1 ) {
				// Use prior values to increment one row.
				increment_G_A() ;
			} else if ( A.dim == A_last.dim - 1 ) {
				// Use prior values to decrement one row.
				decrement_G_A() ;
			} else {
				// Unusual event: rebuild the array.
				G_A.realloc(nactive, nactive) ;
				for ( int j = 0 ; j < nactive ; j++ ) {
					for ( int k = 0 ; k <= j  ; k++ ) {
						
						double tmp = X_A.mult_T(j,k) ;
						G_A.set(j, k, tmp) ;
					}
				}
				for ( int j = 0 ; j < nactive ; j++ ) {
					for ( int k = j + 1 ; k < nactive  ; k++ ) {
						G_A.set(j, k, G_A.get(k,j) ) ;
					}
				}
			}
		}

	void increment_G_A()
	// Increment the G_A array by one extra column and one extra row.
		{
			int newc = 0 ;
			// Find the new index.
			if ( nactive != A.dim ) {
				cout << "Error: A dimension mismatch" << endl ;
				exit(1) ;
			}
			for ( ; newc < nactive ; newc++ ) {
				int k = 0 ;
				for ( ; k < A_last.dim ; k++ ) {
					if ( A_last.get(k) == A.get(newc) )
						break ;
				}
				if ( k == A_last.dim ) {
					break ;
				}
			}
			Matrix G_New(nactive, nactive) ;

			// Copy unchanged elements of the array.
			for ( int i = 0 ; i < newc ; i++ ) {
				for ( int j = 0 ; j <= i ; j++ ) {
					G_New.set(j,i, G_A.get(j,i)) ;
				}
			}
			// Copy shifted elements of the existing array.
			for ( int i = newc ; i < nactive - 1 ; i++ ) {
				for ( int j = 0 ; j < newc ; j++ ) {
					G_New.set(j,i+1, G_A.get(j,i)) ;
				}
				for ( int j = newc ; j <= i ; j++ ) {
					G_New.set(j+1,i+1, G_A.get(j,i)) ;
				}
			}
			// Calculate the new elements.
			for ( int k = 0 ; k < nactive ; k++ ) {
				double tmp = X_A.mult_T(newc, k) ;
				
				//for ( int l = 0 ; l < ndata ; l++ ) {
				//tmp += X_A.get(l, newc) * X_A.get(l, k) ;
				//}
				if ( newc <= k ) {
					G_New.set(newc, k, tmp) ;
				} else {
					G_New.set(k, newc, tmp) ;
				}
			}

			// Copy elements back into the reallocated array.
			G_A.realloc(nactive, nactive) ;
			for ( int i = 0 ; i < nactive ; i++ ) {
				for ( int j = 0 ; j <= i ; j++ ) {
					G_A.set(i,j, G_New.get(j,i)) ;
					G_A.set(j,i, G_New.get(j,i)) ;
				}
			}
		}

	
	void decrement_G_A()
	// Decrement the G_A array by one column and one row.
		{
			int delc = 0 ;
			// Find the new index.
			if ( nactive != A.dim ) {
				cout << "Error: A dimension mismatch" << endl ;
				exit(1) ;
			}
			for ( ; delc < A_last.dim ; delc++ ) {
				int k = 0 ;
				for ( ; k < A.dim ; k++ ) {
					if ( A_last.get(delc) == A.get(k) )
						break ;
				}
				if ( k == A.dim ) {
					break ;
				}
			}
			// delc is the index of the deleted column.
			if ( delc >= A_last.dim ) {
				cout << "Error: did not find deleted index" << endl ;
				exit(1) ;
			}

			Matrix G_New(nactive, nactive) ;

			// Copy unchanged elements of the array.
			for ( int i = 0 ; i < delc && i < nactive ; i++ ) {
				for ( int j = 0 ; j <= i ; j++ ) {
					G_New.set(j,i, G_A.get(j,i)) ;
				}
			}
			// Copy shifted elements of the existing array.
			for ( int i = delc + 1 ; i < nactive + 1 ; i++ ) {
				for ( int j = 0 ; j < delc ; j++ ) {
					if ( j <= i - 1 ) {
						G_New.set(j,i-1, G_A.get(j,i)) ;
					} else {
						G_New.set(i-1,j, G_A.get(j,i)) ;
					}
				}
				for ( int j = delc + 1 ; j < nactive + 1 && j <= i ; j++ ) {
					if ( j <= i ) {
						G_New.set(j-1,i-1, G_A.get(j,i)) ;
					} else {
						G_New.set(i-1,j-1, G_A.get(j,i)) ;
					}
				}
			}
			// Copy elements back into the reallocated array.
			G_A.realloc(nactive, nactive) ;
			for ( int i = 0 ; i < nactive ; i++ ) {
				for ( int j = 0 ; j <= i ; j++ ) {
					G_A.set(i,j, G_New.get(j,i)) ;
					G_A.set(j,i, G_New.get(j,i)) ;
				}
			}
		}

	bool solve_G_A()
	// Find G_A^-1 * I
		{
			G_A_Inv_I.realloc(nactive) ;
			const bool use_incremental_updates = true ;

			bool succeeded = false ;
			// If solve_succeeded == true, the last linear solve worked and
			// we can possibly update the cholesky decomposition.
			// Otherwise, the whole decomposition needs to be recalculated.
			if ( solve_succeeded && use_incremental_updates ) {
				if ( nactive == A_last.dim + 1 && nactive > 2 ) {
					Matrix chol0(chol) ;
					Vector G_row(nactive) ;
					for ( int j = 0 ; j < nactive ; j++ ) {
						G_row.set(j, G_A.get(nactive-1, j) ) ;
					}
					chol.realloc(nactive, nactive) ;
					succeeded = chol.cholesky_add_row(chol0, G_row) ;
					if ( succeeded ) {
						// Back-substitute using the updated cholesky matrix.
						succeeded = chol_backsub() ;
						if ( succeeded ) {
							solve_succeeded = true ;
							return true ;
						}
					} else {
						cout << "Failed to add a row to the Cholesky decomposition" << endl ;
					} 
				} else if ( nactive == A_last.dim - 1 && nactive > 2 ) {
					succeeded = chol.cholesky_remove_row(remove_prop) ;
					if ( succeeded ) {
						// Back-substitute using the updated cholesky matrix.
						succeeded = chol_backsub() ;
						if ( succeeded ) {
							solve_succeeded = true ;
							return true ;
						} else {
							cout << "Failed to remove a row from the Cholesky decomposition" << endl ;
						}
					}
				}
			}

			// Try non-incremental if incremental failed or not possible/requested.
			if ( ! succeeded ) {
				chol.realloc(nactive, nactive) ;
				if ( ! G_A.cholesky(chol) ) {
					cout << "Cholesky failed" << endl ;
					for ( int j = 0 ; j < nactive ; j++ ) {
						int k ;
						// See if this is a new index.
						for ( k = 0 ; k < A_last.dim ; k++ ) {
							if ( A.get(j) == A_last.get(k) ) {
								break ;
							}
						}
						if ( k == A_last.dim ) {
							// The index is new. Exclude it in the future.
							exclude.set( A.get(j), 1) ;
							++num_exclude ;
						}
					}
					solve_succeeded = false ;
					return false ;
				}
			}
			solve_succeeded = chol_backsub() ;
			return solve_succeeded ;
		}
	bool chol_backsub()
	// Perform back substitution on the cholesky matrix to find G_A_Inv_I and A_A
	// Returns true if the solution passes consistency tests.
		{

			//cout << "Cholesky " << endl ;
      //chol.print() ;
			
			// Solve for G_A^-1 * unity
			Vector unity(nactive, 1.0) ;
			chol.cholesky_sub(G_A_Inv_I, unity) ;

			//cout << "G_A_Inv_I " << endl ;
			//G_A_Inv_I.print() ;

			// Test to see if the solution worked.
			Vector test(nactive) ;
			G_A.dot(test, G_A_Inv_I) ;
			double errval = 0.0 ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				errval += fabs(test.get(j)-1.0) ;
				if ( fabs(test.get(j) - 1.0) > 1.0e-04 ) {
					cout << "Cholesky solution test failed\n" ;
					cout << "Error = " << fabs(test.get(j) - 1.0) << endl ;
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
				cout << "A_A Normalization failed" << endl ;
				return false ;
			}
			return true ;
		}
	
	void build_u_A()
		{
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
			if ( fabs(test-1.0) > 1.0e-04 ) {
				cout << "U_A norm test failed" << endl ;
				cout << "Norm = " << test << endl ;
				exit(1) ;
			}

			// Test X_A^T u__A = A_A * I
			Vector testv(nactive,0.0) ;
			X_A.dot_transpose(testv, u_A) ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				if ( fabs(testv.get(j) - A_A) > 1.0e-08 ) {
					cout << "u_A test failed " << endl ;
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
			const double eps = 1.0e-10 ;

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
				// Something strange has happened.
				// Search for the new active set.
				int count = 0 ;
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
							a_trial.set(count, j) ;
							if ( RANK == 0 ) cout << "Adding property " << j << " to the active set" << endl ;
							count++ ;
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
			A.print() ;

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

			for ( int i = 0 ; i < nactive ; i++ ) {
				if ( fabs(w_A.get(i)) > 1.0e-40 ) {
					double gamma_i = -beta.get(A.get(i)) / ( sign.get(i) * w_A.get(i) ) ;
					if ( gamma_i > 0.0 && gamma_i < gamma_lasso ) {
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
#ifdef USE_MPI
			// Sync the beta values to avoid possible divergence between processes.
			MPI_Bcast(beta.vec, nprops, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
#endif
			
#ifdef VERBOSE			
			cout << "New beta: " ;
			beta.print() ;

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
						exit(1) ;
					}
					offset -= beta.get(j) * X.shift[j] / X.scale[j] ;
				}
				for ( int j = 0 ; j < nprops ; j++ ) {
					uns_beta.set(j, beta.get(j) / X.scale[j]) ;
				}
				if ( out.rdbuf() == cout.rdbuf() ) {
					uns_beta.print() ;
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

	void print_error(ostream &out)
	{
		if ( RANK == 0 ) {
			out  << "L1 norm of solution: " << beta.l1norm() << " RMS Error: " << sqrt(sq_error() / ndata) << " Objective fn: " << obj_func_val << endl ;
		}
	}

	bool iteration()
	// Perform a single iteration of the LARS algorithm.
	// Return false when no more iterations can be performed.
		{

			if ( nactive >= nprops - num_exclude ) {
				// No more iterations are possible.
				return false ;
			}

			predict() ;
			correlation() ;

#ifdef VERBOSE
			if ( RANK == 0 ) {
				cout << "Pre-step beta: " << endl ;
				beta.print() ;
				cout << "Prediction: " << endl ;
				mu.print() ;
				cout << " Correlation: " << endl ;
				c.print() ;
				cout << "C_max: " << C_max << endl ;
			}
#endif			
			objective_func() ;

			print_error(cout) ;

			update_active_set() ;

			// build the X_A array.
			build_X_A() ;
			build_G_A() ;

			if ( ! solve_G_A() ) {
				remove_prop = -1 ;
				cout << "Iteration failed" << endl ;
				return false ;
			}
			
#ifdef VERBOSE
			if ( RANK == 0 ) {
				cout << "X_A" << endl ;
				X_A.print() ;
				cout << "G_A" << endl ;
				G_A.print() ;
				cout << "G_A_Inv " << endl ;
				G_A_Inv_I.print() ;
				cout << "A_A " << A_A << endl ;
			}
#endif			

			build_u_A() ;

			update_step_gamma() ;
			update_beta() ;

			if ( RANK == 0 ) {
				cout << "Beta: " << endl ;
				beta.print() ;
				cout << "Y constant offset = " << y.shift << endl ;
				cout << "Unscaled coefficients: " << endl ;
				print_unscaled(cout) ;
			}

			return true ;

		}
};
