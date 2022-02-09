/** Distributed LARS-LASSO algorithm 
		The notation and implementation closely follows
		B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, "Least Angle Regression",
		The Annals of Statistics, 32, 407-499(2004).

		Larry Fried

		The X and X_A matrices are distributed among MPI processes.	 
		If the --distributed_solver argument is given, Cholesky decomposition and all 
    matrices are distributed.  This removes memory limitations in earlier versions
    of the program.

		2/2022
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
		cout << "	 --" << opt->name << " " ;
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
	MPI_Init		 (&argc, &argv);
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
		{"distributed_solver", required_argument, 0, 'd'},
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
	
	string algorithm("lasso") ;				// Algorithm to use: lasso or lars
	bool split_files = false ;				// Read input matrix from split files ?
	bool normalize=true ;							// Whether to normalize the X matrix.
	bool con_grad = false ;						// Whether to use congugate gradient algorithm to solve linear equations.

	bool use_precondition = false ;
	bool distributed_solver = false ;
	// Stopping criteria.	 Default is to calculate all possible solutions.
	
	double max_beta_norm = 1.0e+50 ;	// Maximum L1 norm of solution 
	int max_iterations = 10000000 ;		// Maximum number of LARS iterations.
	double lambda = 0.0 ;							// L1 weighting factor.
	
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
	
	for (	 ; j + 1 <= max_iterations ; j++ ) {
		int status = lars.iteration() ;
		if ( status == 0 ) {
			if ( RANK == 0 ) cout << "Stopping: no more iterations possible" << endl ;
			break ;
		} else if ( status == -1 && last_status == 1 ) {
			// if ( RANK == 0 ) cout << "Iteration failed: continuing" << endl ;
			if ( RANK == 0 ) cout << "Stopping: iteration failed" << endl ;			
			break ;
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
			cout << "Time for iteration " << j + 1 << " = " << elapsed_seconds.count() << " seconds " << endl ;

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


int DLARS::iteration()
	// Perform a single iteration of the LARS algorithm.
	// Return 0 when no more iterations can be performed.
	// Return 1 on success.
	// Return -1 on a failed iteration that could be recovered from.
{

	if ( nactive >= nprops - num_exclude ) {
		// No more iterations are possible.
		return 0 ;
	}

	iterations++ ;
	auto time1 = std::chrono::system_clock::now() ;			
	predict() ;
	auto time2 = std::chrono::system_clock::now() ;
	std::chrono::duration<double> elapsed_seconds = time2 - time1 ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time making prediction = " << elapsed_seconds.count() << endl ;
	}
#endif			
			
	objective_func() ;
	auto time3 = std::chrono::system_clock::now() ;
	elapsed_seconds = time3 - time2 ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time calculating objective function = " << elapsed_seconds.count() << endl ;
	}
#endif			

			
	if ( RANK == 0 ) {
		auto time_rst_1 = std::chrono::system_clock::now() ;

		trajfile << "Iteration " << iterations << endl ;
		print_error(trajfile) ;
		beta.print(trajfile) ;
		print_error(cout) ;

		if ( iterations % 10 == 0 ) {
			print_restart() ;
		}

		auto time_rst_2 = std::chrono::system_clock::now() ;		
		elapsed_seconds = time_rst_2 - time_rst_1 ;
#ifdef TIMING			
		if ( RANK == 0 ) {
			cout << "Time writing restart file = " << elapsed_seconds.count() << endl ;
		}
#endif					
	}

	auto time4a = std::chrono::system_clock::now() ;	
	correlation() ;
	auto time4b = std::chrono::system_clock::now() ;
	elapsed_seconds = time4b - time4a ;

#ifdef TIMING
	if ( RANK == 0 ) {
		cout << "Time calculating correlation = " << elapsed_seconds.count() << endl ;
	}
#endif
			
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
	update_active_set() ;
	auto time5 = std::chrono::system_clock::now() ;
	elapsed_seconds = time5 - time4b ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time updating active set = " << elapsed_seconds.count() << endl ;
	}
#endif			

	// build the X_A array.
	build_X_A() ;

	auto time6 = std::chrono::system_clock::now() ;
	elapsed_seconds = time6 - time5 ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time building X_A = " << elapsed_seconds.count() << endl ;
	}
#endif			

	build_G_A(G_A, true) ;

	auto time7 = std::chrono::system_clock::now() ;
	elapsed_seconds = time7 - time6 ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time building G_A = " << elapsed_seconds.count() << endl ;
	}
#endif			

	if ( build_G_A_here() && ! solve_G_A(true) ) {
		remove_prop = -1 ;
		add_prop = -1 ;
		cout << "Iteration failed" << endl ;
		return -1 ;
	}

	auto time8 = std::chrono::system_clock::now() ;
	elapsed_seconds = time8 - time7 ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time solving G_A = " << elapsed_seconds.count() << endl ;
	}
#endif			

	broadcast_solution() ;
			
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

	auto time9 = std::chrono::system_clock::now() ;
	elapsed_seconds = time9 - time8 ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time building u_A = " << elapsed_seconds.count() << endl ;
	}
#endif			

	update_step_gamma() ;

	auto time10 = std::chrono::system_clock::now() ;
	elapsed_seconds = time10 - time9 ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time updating gamma = " << elapsed_seconds.count() << endl ;
	}
#endif			
			
	update_beta() ;

	auto time11 = std::chrono::system_clock::now() ;
	elapsed_seconds = time11 - time10 ;

#ifdef TIMING			
	if ( RANK == 0 ) {
		cout << "Time updating beta = " << elapsed_seconds.count() << endl ;
	}
#endif			

	if ( RANK == 0 ) {
		cout << "Beta: " << endl ;
		beta.print(cout) ;
		//cout << "Y constant offset = " << y.shift << endl ;
		//cout << "Unscaled coefficients: " << endl ;
		//print_unscaled(cout) ;
	}

	if ( RANK == 0 ) {
		double total_mem = 0.0 ;
		int prec = cout.precision() ;
		cout << "Matrix memory usage on Rank 0:" << endl ;
		total_mem += X.print_memory("X") ;
		total_mem += X_A.print_memory("X_A") ;
		total_mem += G_A.print_memory("G_A") ;
		total_mem += chol.print_memory("chol") ;
		total_mem += pre_con.print_memory("pre_con") ;
		cout << "Total memory = " << std::fixed << std::setprecision(1)
				 << total_mem/1024.0 << " Gb " << endl ;
		cout.precision(prec) ;
		cout << std::scientific ;
	}
	return 1 ;

}

void DLARS::build_G_A(Matrix &G_A_in, bool increment)
	// Build the G_A matrix after X_A has been built.
	// Increment the matrix from previous versions if it is possible and the increment parameter is true.
{
	if ( distributed_solver && ! G_A_in.distributed ) {
		G_A_in.distribute(nactive) ;
	}	
	if ( increment && A.dim == A_last.dim + 1 ) {
		// Use prior values to increment one row.
		if ( distributed_solver ) {
			increment_G_A_dist(G_A_in) ;
		} else {
			increment_G_A(G_A_in) ;
		}
	} else if ( increment && A.dim == A_last.dim - 1 ) {
		// Use prior values to decrement one row.
		if ( distributed_solver ) {
			decrement_G_A_dist(G_A_in) ;
		} else {
			decrement_G_A(G_A_in) ;
		}
	} else {
		// Unusual event: rebuild the array.
		if ( build_G_A_here() ) {
			// Only store G_A on rank 0 for non-distributed solver.
			if ( distributed_solver ) {
				G_A_in.distribute(nactive) ;
			}
			G_A_in.realloc(nactive, nactive) ;
		}
		for ( int j = 0 ; j < nactive ; j++ ) {
			Vector tmp(nactive) ;

			// All processes are used in X_A.mult_T because X_A is a distributed
			// matrix.
			X_A.mult_T_lower(j, tmp) ;
			if ( build_G_A_here() ) {
				if ( G_A_in.row_start <= j && j <= G_A_in.row_end ) {
					for ( int k = 0 ; k <= j	; k++ ) {
						G_A_in.set(j, k, tmp.get(k)) ;																		
					}
				}
				// Transpose elements.
				for ( int k = 0 ; k < j	 ; k++ ) {
					if ( G_A_in.row_start <= k && k <= G_A_in.row_end ) {
						G_A_in.set(k, j, tmp.get(k)) ;
					}
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


void DLARS::increment_G_A(Matrix &G_A_in)
	// Increment the G_A_in array by one extra column and one extra row.
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
				G_New.set(j,i, G_A_in.get(j,i)) ;
			}
		}

		// Copy shifted elements of the existing array.
		for ( int i = newc ; i < nactive - 1 ; i++ ) {
			for ( int j = 0 ; j < newc ; j++ ) {
				G_New.set(j,i+1, G_A_in.get(j,i)) ;
			}
			for ( int j = newc ; j <= i ; j++ ) {
				G_New.set(j+1,i+1, G_A_in.get(j,i)) ;
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
		G_A_in.realloc(nactive, nactive) ;
		for ( int i = 0 ; i < nactive ; i++ ) {
			for ( int j = 0 ; j <= i ; j++ ) {
				G_A_in.set(i,j, G_New.get(j,i)) ;
				G_A_in.set(j,i, G_New.get(j,i)) ;
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


void DLARS::decrement_G_A(Matrix &G_A_in)
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
			G_New.set(j,i, G_A_in.get(j,i)) ;
		}
	}
	// Copy shifted elements of the existing array.
	for ( int i = delc + 1 ; i < nactive + 1 ; i++ ) {
		for ( int j = 0 ; j < delc ; j++ ) {
			if ( j <= i - 1 ) {
				G_New.set(j,i-1, G_A_in.get(j,i)) ;
			} else {
				G_New.set(i-1,j, G_A_in.get(j,i)) ;
			}
		}
		for ( int j = delc + 1 ; j < nactive + 1 && j <= i ; j++ ) {
			if ( j <= i ) {
				G_New.set(j-1,i-1, G_A_in.get(j,i)) ;
			} else {
				G_New.set(i-1,j-1, G_A_in.get(j,i)) ;
			}
		}
	}
	// Copy elements back into the reallocated array.
	G_A_in.realloc(nactive, nactive) ;
	for ( int i = 0 ; i < nactive ; i++ ) {
		for ( int j = 0 ; j <= i ; j++ ) {
			G_A_in.set(i,j, G_New.get(j,i)) ;
			G_A_in.set(j,i, G_New.get(j,i)) ;
		}
	}
}


void DLARS::decrement_G_A_dist(Matrix &G_A_dist)
// Decrement the G_A array by one column and one row.
// The G_A_in array is assumed to be distributed.
{
	if ( ! build_G_A_here() ) return ; 
			
	if ( nactive != A.dim ) {
		cout << "Error: A dimension mismatch" << endl ;
		stop_run(1) ;
	}
	if ( ! G_A_dist.distributed ) {
		cout << "Error in increment_G_A_dist: The matrix was not distributed\n" ;
		stop_run(1) ;
	}

	// Find the new index.
	int delc = 0 ;
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

	Matrix G_New ;
	if ( G_A_dist.distributed ) {
		G_New.distribute(nactive) ;
	}
	G_New.realloc(nactive, nactive) ;

	Vector G_buf(nactive+1) ;
		
	// Copy unchanged elements of the array.
	for ( int j = 0 ; j < nactive + 1 ; j++ ) {
		if ( j == delc ) continue ;

		int j_new = ( j < delc ) ? j : j - 1 ;
		int rank_j_new = G_New.rank_from_row(j_new) ;
		int rank_j_dist = G_A_dist.rank_from_row(j) ;

		if ( rank_j_new == rank_j_dist ) {
			// Data lies on the same MPI rank.
			if ( G_New.row_start <= j_new && j_new <= G_New.row_end ) {
				for ( int i = 0 ; i < delc ; i++ ) {
					G_New.set(j_new,i, G_A_dist.get(j,i)) ;
				}
				for ( int i = delc + 1 ; i < nactive + 1 ; i++ ) {
					G_New.set(j_new,i-1, G_A_dist.get(j,i)) ;
				}					
			}
		} else {
			// Need to transfer data from another rank.
			if ( RANK == rank_j_dist ) {
				for ( int k = 0 ; k < nactive + 1 ; k++ ) {
					G_buf.set(k, G_A_dist.get(j,k) ) ;
				}
#ifdef USE_MPI
				MPI_Send(G_buf.vec, nactive+1, MPI_DOUBLE, rank_j_new, 0, MPI_COMM_WORLD) ;
#endif					
			} else if ( RANK == rank_j_new ) {
#ifdef USE_MPI					
				MPI_Recv(G_buf.vec, nactive+1, MPI_DOUBLE, rank_j_dist, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
#endif										
				for ( int i = 0 ; i < delc ; i++ ) {
					G_New.set(j_new,i, G_buf.get(i) ) ;
				}
				for ( int i = delc + 1 ; i < nactive + 1 ; i++ ) {
					G_New.set(j_new,i-1, G_buf.get(i) ) ;
				}										
			}
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


bool DLARS::solve_G_A(bool use_incremental_updates)
// Find G_A^-1 * I, using either a local or distributed algorithm.
{
			
	auto time1 = std::chrono::system_clock::now() ;

	G_A_Inv_I.realloc(nactive) ;

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

			Vector G_row(nactive,0.0) ;
			int nactive1 = nactive - 1 ;
			
			if ( G_A.row_start <= nactive1 && nactive1 <= G_A.row_end ) {
				for ( int j = 0 ; j < nactive ; j++ ) {
					G_row.set(j, G_A.get(nactive1, j) ) ;
				}
			}
#ifdef USE_MPI
			if ( distributed_solver ) {
				int rank_nactive = G_A.rank_from_row(nactive1) ;
				MPI_Bcast(G_row.vec, nactive, MPI_DOUBLE, rank_nactive, MPI_COMM_WORLD) ;
			}
#endif

			if ( ! distributed_solver ) {
				Matrix chol0(chol) ;
				chol.realloc(nactive, nactive) ;
				auto time1_add = std::chrono::system_clock::now() ;

				succeeded = chol.cholesky_add_row(chol0, G_row) ;

				//cout << "Serial cholesky:\n" ;
				//chol.print() ;
			
				auto time2_add = std::chrono::system_clock::now() ;
				std::chrono::duration<double> elapsed_seconds = time2_add - time1_add ;
				
#ifdef TIMING
				if ( RANK == 0 ) {
					cout << "Time adding cholesky row = " << elapsed_seconds.count() << endl ;
				}
#endif
			} else { /* distributed_solver */
				if ( ! chol.distributed ) {
					chol.distribute(nactive) ;
				}
				Matrix chol0(chol) ;
				chol.realloc(nactive, nactive) ;
				auto time1_add = std::chrono::system_clock::now() ;

				succeeded = chol.cholesky_add_row_distribute(chol0, G_row) ;

				//cout << "Distributed cholesky:\n" ;
				//chol.print() ;

				auto time2_add = std::chrono::system_clock::now() ;
				std::chrono::duration<double> elapsed_seconds = time2_add - time1_add ;
				
#ifdef TIMING
				if ( RANK == 0 ) {
						cout << "Time adding distributed cholesky row = " << elapsed_seconds.count() << endl ;
				}
#endif
			}
			if ( succeeded ) {
				// Back-substitute using the updated cholesky matrix.
				auto time1_back = std::chrono::system_clock::now() ;
				if ( distributed_solver ) {
					succeeded = chol_backsub(G_A, chol) ;
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
			auto time1_rem = std::chrono::system_clock::now() ;				

			if ( distributed_solver ) {
				succeeded = chol.cholesky_remove_row_dist(remove_prop) ;
				auto time2_rem = std::chrono::system_clock::now() ;
				std::chrono::duration<double> rem_seconds = time2_rem - time1_rem ;
#ifdef TIMING								
					if ( RANK == 0 ) {
							cout << "Time removing a variable (distributed) = " << rem_seconds.count() << endl ;	
					}
#endif

			} else {
				succeeded = chol.cholesky_remove_row(remove_prop) ;

				auto time2_rem = std::chrono::system_clock::now() ;
				std::chrono::duration<double> rem_seconds = time2_rem - time1_rem ;
#ifdef TIMING								
					if ( RANK == 0 ) {
							cout << "Time removing a variable (serial) = " << rem_seconds.count() << endl ;	
					}
#endif					
			}

			auto time1_back = std::chrono::system_clock::now() ;
			if ( succeeded ) {
				// Back-substitute using the updated cholesky matrix.
				if ( build_G_A_here() ) {
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
		if ( ! distributed_solver ) {
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
			
		} else {
				
			auto time1_chol_dist = std::chrono::system_clock::now() ;				
			if ( ! chol.distributed ) {
				chol.distribute(nactive) ;
			}			
			chol.realloc(nactive, nactive) ;
			
			G_A.cholesky_distribute(chol) ;
				
			auto time2_chol_dist = std::chrono::system_clock::now() ;								
			std::chrono::duration<double> chol_dist_seconds = time2_chol_dist - time1_chol_dist ;
			if ( RANK == 0 ) {
				cout << "Time solving distributed Cholesky = " << chol_dist_seconds.count() << endl ;
			}							
		}
	}
	if ( distributed_solver ) {
		succeeded = chol_backsub(G_A, chol) ;
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

	//			// DEBUG !!
	//			if ( RANK == 0 ) cout << "Cholesky " << endl ;
	//			 chol.print() ;

	// #ifdef USE_MPI			
	//			MPI_Barrier(MPI_COMM_WORLD) ;
	//			cout.flush() ;
	// #endif			
			
	//			if ( RANK == 0 ) cout << "Cholesky dist" << endl ;
	//			chol_in.print() ;

	// #ifdef USE_MPI			
	//			MPI_Barrier(MPI_COMM_WORLD) ;
	//			cout.flush() ;
	// #endif			
			
	const double eps_fail = 1.0e-04 ;	 // Max allowed solution error.

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

#ifdef TIMING		
	if ( RANK == 0 ) {
		if ( distributed_solver ) {
			cout << "Time substituting distributed Cholesky = " << sub_seconds.count() << endl ;
		}	else {
			cout << "Time substituting local Cholesky = " << sub_seconds.count() << endl ;
		}
	}
#endif
	
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

void DLARS::predict() 
	// Calculated predicted values of y (mu hat, Eq. 1.2).	Update previous prediction
	// based on u_A and gamma_use.
{
	if ( mu.size() != ndata ) {
		cout << "Error:	 matrix dim mismatch" << endl ;
		stop_run(1) ;
	}

	if ( u_A.dim == 0 ) {
		// First iteration.
		X.dot(mu, beta) ;
	} else {
		for ( int j = 0 ; j < ndata ; j++ ) {
			mu.set(j,	 mu.get(j) + gamma_use * u_A.get(j) ) ;
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

void DLARS::predict_all() 
// Calculated predicted values of y (mu hat, Eq. 1.2) with no updating based on u_A.
{
	if ( mu.size() != ndata ) {
		cout << "Error:	 matrix dim mismatch" << endl ;
		stop_run(1) ;
	}
	X.dot(mu, beta) ;

#ifdef VERBOSE			
	cout << "Mu = " << endl ;
	mu.print() ;
#endif
}

double DLARS::sq_error()
// Squared error Eq. 1.3
{
	double result = 0.0 ;
	for ( int j = 0 ; j < ndata ; j++ ) {
		result += (y.get(j)-mu.get(j)) * (y.get(j) - mu.get(j) ) ;
	}
	return(result) ;
}

void DLARS::objective_func()
	// Calculate optimization objective function based on the requested
	// regularization parameter lambda.	 This should be called after
	// predict_all() or predict().
{
	obj_func_val = 0.5 * sq_error() / ndata + lambda * beta.l1norm() ;
}

void DLARS::correlation()		
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

void DLARS::build_X_A()
{		
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

int DLARS::restart(string filename)
// Restart from the given file.
{
	int iter ;
	ifstream inf(filename) ;
	if ( ! inf.good() ) {
		cout << "Could not open " << filename << " for restart" << endl ;
		stop_run(1) ;
	}
	while (1) {
		string s ;
		int iter_tmp ;
		// Get the iteration number.
		inf >> s >> iter_tmp ;
		//iter-- ;
		if ( inf.eof() || ! inf.good() ) break ;

		iter = iter_tmp ;
		// Get the objective function from the next line.
		for ( int j = 0 ; j < 15 ; j++ ) {
			inf >> s ;
			if ( j == 10 ) {
				obj_func_val = stod(s) ;
				//cout << "OBJ FUNC: " << obj_func_val << endl ;
			}
		}

		if ( inf.eof() || ! inf.good() ) {
			cout << "Could not read the objective function value from " << filename << endl ;
			stop_run(1) ;
		}
		A.clear() ;

		// Read all of the beta values.
		string line ;
		getline(inf,line) ;
		beta.read_sparse(inf) ;
		for ( int j = 0 ; j < nprops ; j++ ) {
			if ( fabs(beta.get(j)) > 0.0 ) {
				A.push(j) ;
			}
		}
		if ( inf.eof() || ! inf.good() ) {
			cout << "Could not read the parameter values from " << filename << endl ;
			stop_run(1) ;
		}

		getline(inf,line) ;
		if ( line.find("Exclude") != string::npos ) {
			exclude.read_sparse(inf) ;
		}

		if ( line.find("Mu") != string::npos ) {
			mu.read_sparse(inf) ;
		}
	}
		
	inf.close() ;
		
	nactive = A.dim ;
	iterations = iter - 1 ;

	bool con_grad_save = solve_con_grad ;
	solve_con_grad = false ;
	//predict_all() ;
	objective_func() ;

	if ( RANK == 0 ) cout << "Restart: calculating correlation\n" ;
	correlation() ;

	if ( RANK == 0 ) cout << "Restart: building X_A\n" ;
	build_X_A() ;

	if ( RANK == 0 ) cout << "Restart: building G_A\n" ;
	build_G_A(G_A, false) ;
	
	if ( build_G_A_here() ) {
		if ( RANK == 0 ) cout << "Restart: solving G_A\n" ;
		solve_G_A(false) ;
	}
	broadcast_solution() ;

	// Full calculation of pre-conditioner.
	if ( use_precondition ) {
		pre_con.resize(nactive, nactive) ;
		Matrix chol_precon(nactive, nactive) ;
						
		if ( ! G_A.cholesky(chol_precon) ) {
			cout << "Cholesky decomposition for pre-conditioning failed\n" ;
			stop_run(1) ;
		}
		chol_precon.cholesky_invert(pre_con) ;
	}
		
	solve_con_grad = con_grad_save ;
		
	return iter -1 ;
}

void DLARS::broadcast_solution()
	// Broadcast results of solving G_A.
{
#ifdef USE_MPI
	if ( ! distributed_solver ) {
		if ( RANK != 0 ) {
			G_A_Inv_I.realloc(nactive) ;
			A.realloc(nactive) ;
		}
		
		MPI_Bcast(&nactive, 1, MPI_INT, 0, MPI_COMM_WORLD) ;
		MPI_Bcast(A.vec, nactive, MPI_INT, 0, MPI_COMM_WORLD) ;
		MPI_Bcast(G_A_Inv_I.vec, nactive, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
		MPI_Bcast(&A_A, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) ;
	}
#endif
}

	
bool DLARS::solve_G_A_con_grad()
	// Use the conjugate gradient method to find G_A_Inv_I and A_A
	// Does not use cholesky decomposition.
	// Returns true if the solution passes consistency tests.
{

	const double eps_fail = 1.0e-04 ;	 // Max allowed solution error.
	const double eps_con_grad = 1.0e-08 ;	 // Conjugate gradient error tolerance.

	// Solve for G_A^-1 * unity
	Vector unity(nactive, 1.0) ;

	// if ( RANK == 0 ) {
	//	// cout << "G_A_Inv_I guess = \n" ;
	//	G_A_Inv_I.print(cout) ;
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
	//	cout << "G_A_Inv_I solution" << endl ;
	//	G_A_Inv_I.print(cout) ;
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

void DLARS::build_u_A()
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


void DLARS::reduce_active_set() 
	// Reduce the active set of directions to those having maximum correlation.
	// See Eq. 3.6
{
	// Undo the change in the active set.
	if ( RANK == 0 ) cout << "Will remove property " << A.get(remove_prop) << " from the active set" << endl ;
	A.remove(remove_prop) ;
	nactive = A.dim ;
	gamma_lasso = 1.0e20 ;
}

void DLARS::update_active_set() 
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


void DLARS::update_step_gamma()
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


void DLARS::update_lasso_gamma()
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


void DLARS::update_beta()
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

	
void DLARS::print_unscaled(ostream &out) 
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
