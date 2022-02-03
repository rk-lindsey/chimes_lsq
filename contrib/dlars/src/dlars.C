/** Distributed LARS-LASSO algorithm 
		The notation and implementation closely follows
		B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, "Least Angle Regression",
    		The Annals of Statistics, 32, 407-499(2004).
    Larry Fried

		The X and X_A matrices are distributed.  Code works with MPI.

    12/2018
**/

// #define VERBOSE // Define for extra output

#include<math.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string.h>
#include<getopt.h>
#include <chrono>

#ifdef USE_MPI
#include <mpi.h>
#endif


int RANK ;
int NPROCS ;


using namespace std ;

#include "Vector.h"
#include "IntVector.h"
#include "Matrix.h"
#include "DLARS.h"

void display_usage(struct option* opt)
{
	cout << "Recognized options" << endl ;
	while ( opt->name != NULL ) {
		cout << "  --" << opt->name << " " ;
		if ( opt->has_arg == required_argument ) {
			cout << "ARG required" ;
		} else if ( opt->has_arg == no_argument ) {
			cout << "NO ARG" ;
		} else {
			cout << "ARG optional" ;
		}
		cout << endl ;
		++opt ;
	}
}
		
	
int main(int argc, char **argv)
{

#ifdef USE_MPI	
	MPI_Init     (&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
  MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
#else
	RANK = 0 ;
	NPROCS = 1 ;
#endif	

	if ( RANK == 0 ) {
		cout << "Distributed LARS algorithm" << endl ;
		cout << scientific ;
		cout.precision(12) ;
	}

	
	if ( argc < 4 ) {
		if ( RANK == 0 ) cout << "Not enough args " << endl ;
		stop_run(1) ;
	}
	string xname(argv[1]) ;
	string yname(argv[2]) ;
	string dname(argv[3]) ;

	static struct option long_options[] =
	{
		{"algorithm", required_argument, 0, 'a'},
		{"distributed-solver", required_argument, 0, 'd'},
		{"iterations", required_argument, 0, 'i'},
		{"lambda", required_argument, 0, 'l'},
		{"max_norm", required_argument, 0, 'm'},
		{"normalize", required_argument, 0, 'n'},
		{"con_grad", no_argument, 0, 'c'},
		{"precondition", no_argument, 0, 'p'},
		{"restart", required_argument, 0, 'r'},
		{"split_files", no_argument, 0, 's'},
		{"weights", required_argument, 0, 'w'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	} ;

	int option_index = 0 ;
	int opt_type ;

	// Default options.
	
	string algorithm("lasso") ;       // Algorithm to use: lasso or lars
	bool split_files = false ;        // Read input matrix from split files ?
	bool normalize=true ;             // Whether to normalize the X matrix.
	bool con_grad = false ;           // Whether to use congugate gradient algorithm to solve linear equations.

	bool use_precondition = false ;
	bool distributed_solver = false ;
	// Stopping criteria.  Default is to calculate all possible solutions.
	
	double max_beta_norm = 1.0e+50 ;  // Maximum L1 norm of solution 
	int max_iterations = 10000000 ;   // Maximum number of LARS iterations.
	double lambda = 0.0 ;             // L1 weighting factor.
	
	string weight_file("") ;
	string restart_file ;

	while (1) {
		// Colons in string indicate required arguments.
		opt_type = getopt_long(argc, argv, "a:d:i:l:m:n:cpr:sw:h", long_options, &option_index) ;
		if ( opt_type == -1 ) break ;
		switch ( opt_type ) {
		case 'a':
			algorithm = string(optarg) ;
			break ;
		case 'd':
			if ( optarg[0] == 'y' ) {
				distributed_solver = true ;
			} else if ( optarg[0] == 'n' ) {
				distributed_solver = false ;
			} else {
				cerr << "--distributed-solver arg should be y or n" ;
				stop_run(1) ;
			}
			break ;			
		case 'i':
			max_iterations = atoi(optarg) ;
			break ;
		case 'l':
		  lambda = atof(optarg) ;
			break ;			
		case 'm':
			max_beta_norm = atof(optarg) ;
			break ;
		case 'n':
			if ( optarg[0] == 'y' ) {
				normalize = true ;
			} else if ( optarg[0] == 'n' ) {
				normalize = false ;
			} else {
				cerr << "--normalize arg should be y or n" ;
				stop_run(1) ;
			}
			break ;
		case 'p':
			use_precondition = true ;
			break ;
		case 'c':
			con_grad = true ;
			break ;
		case 'r':
			restart_file = string(optarg) ;
			break ;
		case 's':
			split_files = true ;
			break ;
		case 'w':
			weight_file=string(optarg) ;
			break ;
		case 'h':
			// Help !
			display_usage(long_options) ;
			stop_run(0) ;
		default:
			if ( RANK == 0 ) cout << "Unrecognized option: " << (char) opt_type << endl ;
			abort() ;
		}
	}

	bool do_lasso = true ;
	if ( algorithm == "lars" ) {
		do_lasso = false ;
	} else if ( algorithm == "lasso" ) {
		do_lasso = true ;
	} else {
		if ( RANK == 0 ) cout << "Error: unrecognized algorithm: " << algorithm << endl ;
		stop_run(1) ;
	}
	
	if ( RANK == 0 ) {
		cout << " ...options read." << endl;
	}
		
	int nprops, ndata ;
	
	Matrix xmat ;
	if ( split_files ) {
		// Read the X matrix from multiple split files, as output by chimes_lsq
		
		if ( RANK == 0 ) {
			cout << " ...reading split xmat." << endl;
		}
		
		xmat.read_split_files(xname.c_str(), dname.c_str()) ;
		if ( RANK == 0 ) {
			cout << " finished." << endl;
		}		
		
		nprops= xmat.dim2 ;
		ndata = xmat.dim1 ;
	} else {
	
		if ( RANK == 0 ) {
			cout << " ...reading single xmat." << endl;
		}
			
		// Read the X matrix from a single file.
		ifstream xfile(xname) ;
		if ( ! xfile.is_open() ) {
			if ( RANK == 0 ) cout << "Could not open " << xname << endl ;
			stop_run(1) ;
		}
		ifstream dfile(dname) ;
		if ( ! dfile.is_open() ) {
			if ( RANK == 0 ) cout << "Error: could not open " << dname << endl ;
			stop_run(1) ;
		}
		dfile >> nprops >> ndata ;
		xmat.read(xfile, ndata, nprops, true) ;
	}
	
	if ( RANK == 0 ) {
		cout << " ...xmat read." << endl;
	}


	ifstream yfile(yname) ;
	if ( ! yfile.is_open() ) {
		if ( RANK == 0 ) cout << "Could not open " << yname << endl ;
		stop_run(1) ;
	}
	
	Vector yvec ;
	yvec.read(yfile, ndata) ;
	
	if ( RANK == 0 ) {
		cout << " ...yvec read." << endl;
	}	

	Vector weights ;

	if ( ! weight_file.empty() ) {
		
		ifstream weight_stream(weight_file) ;
		if ( ! weight_stream.is_open() ) {
			if ( RANK == 0 ) cout << "Could not open " << weight_file << endl ;
			stop_run(1) ;
		}
		weights.read(weight_stream, ndata) ;
		xmat.scale_rows(weights) ;
		yvec.scale(yvec,weights) ;
		
		if ( RANK == 0 ) {
			cout << " ...weights read." << endl;
		}			
	}
	
	if ( normalize ) {
		xmat.normalize() ;
		xmat.check_norm() ;
		
		if ( RANK == 0 ) {
			cout << " ...xmat normalized." << endl;
		}		
		yvec.normalize() ;
		yvec.check_norm() ;
		
		if ( RANK == 0 ) {
			cout << " ...yvec normalized." << endl;
		}			
	}

#if(0)	
	ofstream xtest("Xtest.txt") ;
	xmat.print(xtest) ;

	ofstream ytest("Ytest.txt") ;
	yvec.print(ytest) ;
#endif	
	
#ifdef VERBOSE	
	if ( RANK == 0 ) {
		cout << "Normalized X matrix" << endl ;
		xmat.print() ;

		cout << "Normalized Y vector" << endl ;
		for ( int j = 0 ; j < yvec.dim ; j++ ) {
			cout << yvec.get(j) << " " ;
		}
		cout << endl ;
	}
#endif	

	DLARS lars(xmat, yvec, lambda) ;
	lars.do_lasso = do_lasso ;
	if ( RANK == 0 ) {
		if ( do_lasso ) {
			cout << "Using the LASSO algorithm\n" ;
		} else {
			cout << "Using the LARS algorithm\n" ;
		}
	}

	lars.solve_con_grad = con_grad ;
	lars.use_precondition = use_precondition ;
	lars.distributed_solver = distributed_solver ;
	
	Vector last_beta(nprops) ; // Last good coefficients

	double last_obj_func = 1.0e50 ;

	const double eps = 1.0e-10 ;

	int j = 0 ;
	if ( ! restart_file.empty() ) {
		j = lars.restart(restart_file) ;
		last_obj_func = lars.obj_func_val ;
		last_beta = lars.beta ;
	}
	int last_status = 1 ;
	auto time1 = std::chrono::system_clock::now() ;
	
	for (  ; j + 1 <= max_iterations ; j++ ) {
		int status = lars.iteration() ;
		if ( status == 0 ) {
			if ( RANK == 0 ) cout << "Stopping: no more iterations possible" << endl ;
			break ;
		} else if ( status == -1 && last_status == 1 ) {
			if ( RANK == 0 ) cout << "Iteration failed: continuing" << endl ;
			continue ;
		} 

		if ( RANK == 0 ) cout << "Finished iteration " << j + 1 << endl ;

		// Check against exit conditions.
		if ( lars.beta.l1norm() > max_beta_norm ) {
			if ( RANK == 0 ) cout << "Stopping: Maximum L1 norm exceeded" << endl ;
			break ;
		}
		else if ( lars.obj_func_val > last_obj_func + eps ) {
			if ( RANK == 0 ) cout << "Stopping: Objective function increased" << endl ;
			break ;
		}
		last_beta = lars.beta ;
		last_obj_func = lars.obj_func_val ;
		last_status = status ;

		auto time2 = std::chrono::system_clock::now() ;
		std::chrono::duration<double> elapsed_seconds = time2 - time1 ;
		
		if ( RANK == 0 )
			cout << "Time for iteration " << j << " = " << elapsed_seconds.count() << " seconds " << endl ;

		time1 = std::chrono::system_clock::now() ;
	}

	if ( RANK == 0 ) 
	{
		cout.precision(12) ;
		cout << "Final values:" << endl ;
		cout << "Beta: " << endl ;
		last_beta.print(cout) ;
	}

	lars.beta = last_beta ;
	lars.predict_all() ;
	lars.correlation() ;
	
	if ( RANK == 0 ) cout << "Prediction: " << endl ;
	lars.mu.print(cout) ;

	if ( RANK == 0 ) cout << "Sq Error " << lars.sq_error() << endl ;

	
	if ( RANK == 0 ) 
	{
		// Print out the unscaled coefficients to X.txt.
		
		ofstream xtxt("x.txt") ;
		xtxt.precision(12) ;
		xtxt << scientific ;
		if ( ! xtxt.is_open() ) 
		{
			if ( RANK == 0 ) cout << "Error: could not open x.txt" << endl ;
			stop_run(1) ;
		}
		lars.print_unscaled(xtxt) ;

		// Print out the predicted values based on unscaled coefficients

		ofstream Axfile("Ax.txt") ;
		Axfile.precision(12) ;
		Axfile << scientific ;
		if ( ! Axfile.is_open() ) 
		{
			if ( RANK == 0 ) cout << "Error: could not open Ax.txt" << endl ;
			stop_run(1) ;
		}
		if ( ! weight_file.empty() )
			lars.print_unshifted_mu(Axfile, weights) ;
		else
			lars.print_unshifted_mu(Axfile) ;
	}
	
#ifdef VERBOSE
	if ( RANK == 0 ) cout << " Correlation: " << endl ;
	lars.c.print() ;
#endif	

#ifdef USE_MPI
	MPI_Finalize() ;
#endif		
}


void stop_run(int stat) 
{
#ifdef USE_MPI
	MPI_Abort(MPI_COMM_WORLD,stat) ;
#else
	exit(stat) ;
#endif
}

void DLARS::build_G_A()
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
		if ( build_G_A_here() ) {
			// Only store G_A on rank 0 for non-distributed solver.
			G_A.realloc(nactive, nactive) ;
		}
		for ( int j = 0 ; j < nactive ; j++ ) {
			Vector tmp(nactive) ;

			// All processes are used in X_A.mult_T because X_A is a distributed
			// matrix.
			X_A.mult_T_lower(j, tmp) ;
			if ( build_G_A_here() ) {
				for ( int k = 0 ; k <= j  ; k++ ) {
					G_A.set(j, k, tmp.get(k)) ;                                    
				}
			}
		}
		if ( build_G_A_here() ) {
			for ( int j = 0 ; j < nactive ; j++ ) {
				for ( int k = j + 1 ; k < nactive  ; k++ ) {
					G_A.set(j, k, G_A.get(k,j) ) ;
				}
			}
		}
	}
}

bool DLARS::build_G_A_here()
// Should the G_A matrix be created on this rank ?
{
	return ( distributed_solver || RANK == 0 ) ;
}


void DLARS::increment_G_A()
	// Increment the G_A array by one extra column and one extra row.
{
	if ( nactive != A.dim ) {
		cout << "Error: A dimension mismatch" << endl ;
		stop_run(1) ;
	}

	// Find the new index (newc).
	int newc = 0 ;
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
	Matrix G_New ;

	if ( build_G_A_here() ) {
		G_New.resize(nactive, nactive) ;

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
	}
	// Calculate the new elements.
	Vector tmp(nactive) ;
	X_A.mult_T(newc, tmp) ;
	for ( int k = 0 ; k < nactive ; k++ ) {
		//double tmp = X_A.mult_T(newc, k) ;
				
		//for ( int l = 0 ; l < ndata ; l++ ) {
		//tmp += X_A.get(l, newc) * X_A.get(l, k) ;
		//}
		if ( build_G_A_here() ) {
			if ( newc <= k ) {
				G_New.set(newc, k, tmp.get(k)) ;
			} else {
				G_New.set(k, newc, tmp.get(k)) ;
			}
		}
	}

	// Copy elements back into the reallocated array.
	if ( build_G_A_here() ) {
		G_A.realloc(nactive, nactive) ;
		for ( int i = 0 ; i < nactive ; i++ ) {
			for ( int j = 0 ; j <= i ; j++ ) {
				G_A.set(i,j, G_New.get(j,i)) ;
				G_A.set(j,i, G_New.get(j,i)) ;
			}
		}
	}
}

void DLARS::increment_G_A_dist(Matrix &G_A_dist)
// Increment the G_A_dist array by one extra column and one extra row.
// The array G_A_dist is assumed to be distributed.
{
	if ( nactive != A.dim ) {
		cout << "Error: A dimension mismatch" << endl ;
		stop_run(1) ;
	}

	if ( ! G_A_dist.distributed ) {
		cout << "Error in increment_G_A_dist: The matrix was not distributed\n" ;
		stop_run(1) ;
	}
	
	// Find the new index (newc).
	int newc = 0 ;
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
	Matrix G_New ;

	if ( build_G_A_here() ) {

		if ( G_A_dist.distributed ) {
			G_New.distribute(nactive) ;
		}
		G_New.realloc(nactive, nactive) ;
			

		Vector G_buf(nactive-1) ;

		// Copy unchanged elements of the array.
		for ( int j = 0 ; j < nactive - 1 ; j++ ) {
			if ( j == newc ) continue ;
			
			int j_new = ( j < newc ) ? j : j + 1 ;
			int rank_j_new = G_New.rank_from_row(j_new) ;
			int rank_j_dist = G_A_dist.rank_from_row(j) ;

			if ( rank_j_new == rank_j_dist ) {
				// Data lies on the same MPI rank.
				if ( G_New.row_start <= j_new && j_new <= G_New.row_end ) {
					for ( int i = 0 ; i < newc ; i++ ) {
						G_New.set(j_new,i, G_A_dist.get(j,i)) ;
					}
					for ( int i = newc ; i < nactive - 1 ; i++ ) {
						G_New.set(j_new,i+1, G_A_dist.get(j,i)) ;
					}					
				}
			} else {
				// Need to transfer data from another rank.

				if ( RANK == rank_j_dist ) {
					for ( int k = 0 ; k < nactive - 1 ; k++ ) {
						G_buf.set(k, G_A_dist.get(j,k) ) ;
					}
#ifdef USE_MPI
					MPI_Send(G_buf.vec, nactive-1, MPI_DOUBLE, rank_j_new, 0, MPI_COMM_WORLD) ;
#endif					
				} else if ( RANK == rank_j_new ) {
#ifdef USE_MPI					
					MPI_Recv(G_buf.vec, nactive-1, MPI_DOUBLE, rank_j_dist, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
#endif										
					for ( int i = 0 ; i < newc ; i++ ) {
						G_New.set(j_new,i, G_buf.get(i) ) ;
					}
					for ( int i = newc ; i < nactive-1 ; i++ ) {
						G_New.set(j_new,i+1, G_buf.get(i) ) ;
					}										
				}
			}
		}
	}
	
	// Calculate the new elements.
	Vector tmp(nactive) ;
	X_A.mult_T(newc, tmp) ;

	if ( build_G_A_here() ) {
		for ( int k = 0 ; k < nactive ; k++ ) {
		//double tmp = X_A.mult_T(newc, k) ;
				
		//for ( int l = 0 ; l < ndata ; l++ ) {
		//tmp += X_A.get(l, newc) * X_A.get(l, k) ;
		//}
			if ( G_New.row_start <= newc && newc <= G_New.row_end ) {
				G_New.set(newc, k, tmp.get(k)) ;
			}

			// Symmetric matrix.
			if ( G_New.row_start <= k && k <= G_New.row_end ) {
				G_New.set(k, newc, tmp.get(k)) ;
			}
		}

		// Copy elements back into the reallocated array.
		G_A_dist.realloc(nactive, nactive) ;
		for ( int j = G_A_dist.row_start ; j <= G_A_dist.row_end ; j++ ) {
			for ( int i = 0 ; i < nactive ; i++ ) {
				G_A_dist.set(j,i, G_New.get(j,i)) ;
			}
		}
	}
}


void DLARS::decrement_G_A()
	// Decrement the G_A array by one column and one row.
{
	if ( ! build_G_A_here() ) return ; 
			
	int delc = 0 ;
	// Find the new index.
	if ( nactive != A.dim ) {
		cout << "Error: A dimension mismatch" << endl ;
		stop_run(1) ;
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
		stop_run(1) ;
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


bool DLARS::solve_G_A(bool use_incremental_updates)
// Find G_A^-1 * I
{
			
	auto time1 = std::chrono::system_clock::now() ;

	G_A_Inv_I.realloc(nactive) ;

	Matrix G_A_dist ;
	Matrix chol_dist ;

	if ( distributed_solver ) {
		G_A_dist.distribute(nactive) ;
		chol_dist.distribute(nactive) ;
		G_A_dist.realloc(nactive, nactive) ;
		chol_dist.realloc(nactive, nactive) ;				
	}
	
		
	bool succeeded = false ;
	// If solve_succeeded == true, the last linear solve worked and
	// we can possibly update the cholesky decomposition.
	// Otherwise, the whole decomposition needs to be recalculated.
	if ( solve_con_grad ) {
		solve_succeeded = solve_G_A_con_grad() ;
		succeeded = solve_succeeded ;
		if ( ! solve_succeeded ) {
			if ( RANK == 0 ) {
				cout << "Conjugate gradient method failed. " << endl ;
				cout << "Trying cholesky instead \n" ;
			}
		} else {

#ifdef TIMING
			auto time2 = std::chrono::system_clock::now() ;
			std::chrono::duration<double> elapsed_seconds = time2 - time1 ;
			if ( RANK == 0 ) {
				cout << "Time solving equations = " << elapsed_seconds.count() << endl ;
			}
#endif
			return true ;
		}
	}
	if ( solve_succeeded && use_incremental_updates ) {
		if ( nactive == A_last.dim + 1 && nactive > 2 ) {
			Matrix chol0(chol) ;
			Vector G_row(nactive) ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				G_row.set(j, G_A.get(nactive-1, j) ) ;
			}
			chol.realloc(nactive, nactive) ;
			auto time1_add = std::chrono::system_clock::now() ;
			succeeded = chol.cholesky_add_row(chol0, G_row) ;
			auto time2_add = std::chrono::system_clock::now() ;
			std::chrono::duration<double> elapsed_seconds = time2_add - time1_add ;

#ifdef TIMING
			if ( RANK == 0 ) {
				cout << "Time adding cholesky row = " << elapsed_seconds.count() << endl ;
			}
#endif					
					
			if ( succeeded ) {
				// Back-substitute using the updated cholesky matrix.
				auto time1_back = std::chrono::system_clock::now() ;
				if ( distributed_solver ) {
					succeeded = chol_backsub(G_A_dist, chol_dist) ;
				} else if ( build_G_A_here() ) {
					succeeded = chol_backsub(G_A, chol) ;
				}
				if ( succeeded ) {
					solve_succeeded = true ;

					auto time2 = std::chrono::system_clock::now() ;
					std::chrono::duration<double> elapsed_seconds = time2 - time1 ;
					std::chrono::duration<double> backsub_seconds = time2 - time1_back ;

#ifdef TIMING							
					if ( RANK == 0 ) {
						cout << "Time back-substituting = " << backsub_seconds.count() << endl ;
						cout << "Time solving equations = " << elapsed_seconds.count() << endl ;
					}
#endif							
							
					return true ;
				}
			} else {
				if ( RANK == 0 ) {
					cout << "Failed to add a row to the Cholesky decomposition" << endl ;
					cout << "Will perform a non-incremental Cholesky decomposition" << endl ;
				}
						
			} 
		} else if ( nactive == A_last.dim - 1 && nactive > 2 ) {
			auto time1_back = std::chrono::system_clock::now() ;											
			succeeded = chol.cholesky_remove_row(remove_prop) ;
			if ( succeeded ) {
				// Back-substitute using the updated cholesky matrix.
				if ( distributed_solver ) {
					succeeded = chol_backsub(G_A_dist, chol_dist) ;
				} else if ( build_G_A_here() ) {
					succeeded = chol_backsub(G_A, chol) ;
				}
				if ( succeeded ) {
					solve_succeeded = true ;

					auto time2 = std::chrono::system_clock::now() ;
					std::chrono::duration<double> elapsed_seconds = time2 - time1 ;
					std::chrono::duration<double> backsub_seconds = time2 - time1_back ;														
#ifdef TIMING								
					if ( RANK == 0 ) {
						cout << "Time back-substituting = " << backsub_seconds.count() << endl ;							               				cout << "Time solving equations = " << elapsed_seconds.count() << endl ;
								
					}
#endif
							
					return true ;
				} else {
					if ( RANK == 0 ) {
						cout << "Failed to remove a row from the Cholesky decomposition" << endl ;
						cout << "Will perform a non-incremental decomposition" << endl ;
					}
				}
			}
		}
	}

	// Try non-incremental if incremental failed or not possible/requested.
	if ( ! succeeded ) {
		chol.realloc(nactive, nactive) ;
		auto time1_chol = std::chrono::system_clock::now() ;											

		if ( ! G_A.cholesky(chol) ) {
			if ( RANK == 0 ) cout << "Non-incremental Cholesky failed" << endl ;
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
				
		auto time2_chol = std::chrono::system_clock::now() ;
		std::chrono::duration<double> chol_seconds = time2_chol - time1_chol ;
#ifdef TIMING				
		if ( RANK == 0 ) {
			cout << "Time solving full cholesky = " << chol_seconds.count() << endl ;
		}
#endif

		if ( distributed_solver ) {
			auto time1_chol_dist = std::chrono::system_clock::now() ;				

			for ( int j = G_A_dist.row_start ; j <= G_A_dist.row_end ; j++ ) {
				for ( int k = 0 ; k < nactive ; k++ ) {
					G_A_dist.set(j,k, G_A.get(j,k)) ;
				}
			}
			G_A_dist.cholesky_distribute(chol_dist) ;
				
			auto time2_chol_dist = std::chrono::system_clock::now() ;								
			std::chrono::duration<double> chol_dist_seconds = time2_chol_dist - time1_chol_dist ;
			if ( RANK == 0 ) {
				cout << "Time solving distributed Cholesky = " << chol_dist_seconds.count() << endl ;
			}							
		}
	}
	if ( distributed_solver ) {
		succeeded = chol_backsub(G_A_dist, chol_dist) ;
	} else if ( build_G_A_here() ) {
		succeeded = chol_backsub(G_A, chol) ;
	}
	auto time2 = std::chrono::system_clock::now() ;
	std::chrono::duration<double> elapsed_seconds = time2 - time1 ;

	if ( RANK == 0 ) {
		cout << "Time solving equations = " << elapsed_seconds.count() << endl ;
	}			
	return solve_succeeded ;
}


bool DLARS::chol_backsub(Matrix &G_A_in, Matrix &chol_in)
	// Perform back substitution on the cholesky matrix to find G_A_Inv_I and A_A
	// Returns true if the solution passes consistency tests.
{

	// 			// DEBUG !!
	// 			if ( RANK == 0 ) cout << "Cholesky " << endl ;
	//       chol.print() ;

	// #ifdef USE_MPI			
	// 			MPI_Barrier(MPI_COMM_WORLD) ;
	// 			cout.flush() ;
	// #endif			
			
	// 			if ( RANK == 0 ) cout << "Cholesky dist" << endl ;
	// 			chol_in.print() ;

	// #ifdef USE_MPI			
	// 			MPI_Barrier(MPI_COMM_WORLD) ;
	// 			cout.flush() ;
	// #endif			
			
	const double eps_fail = 1.0e-04 ;  // Max allowed solution error.

	// Solve for G_A^-1 * unity
	Vector unity(nactive, 1.0) ;


	auto time1_sub = std::chrono::system_clock::now() ;								

	if ( distributed_solver ) {
		chol_in.cholesky_sub_distribute(G_A_Inv_I, unity) ;
	} else {
		chol_in.cholesky_sub(G_A_Inv_I, unity) ;
	}		
	auto time2_sub = std::chrono::system_clock::now() ;			
	std::chrono::duration<double> sub_seconds = time2_sub - time1_sub ;

	if ( RANK == 0 ) {
		if ( distributed_solver ) {
			cout << "Time substituting distributed Cholesky = " << sub_seconds.count() << endl ;
		}	else {
			cout << "Time substituting local Cholesky = " << sub_seconds.count() << endl ;
		}
	}			
	// if ( RANK == 0 ) cout << "Distributed G_A_Inv_I:\n" ;
	// G_A_Inv_I.print_all(cout) ;
			
	//			Vector G_A_Inv_I_cg(nactive, 0.0) ;
	//			G_A.con_grad(G_A_Inv_I_cg, unity, nactive, 3, 1.0e-08) ;
			
	//			cout << "G_A_Inv_I " << endl ;
	//			G_A_Inv_I.print(cout) ;

	//			cout << "G_A_Inv_I_cg " << endl ;
	//			G_A_Inv_I_cg.print(cout) ;
			
	// Test to see if the solution worked.
	Vector test(nactive) ;
	G_A_in.dot(test, G_A_Inv_I) ;
	double errval = 0.0 ;
	for ( int j = 0 ; j < nactive ; j++ ) {
		errval += fabs(test.get(j)-1.0) ;
		if ( fabs(test.get(j) - 1.0) > eps_fail ) {
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

