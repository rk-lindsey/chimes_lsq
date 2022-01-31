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
	
	// Stopping criteria.  Default is to calculate all possible solutions.
	
	double max_beta_norm = 1.0e+50 ;  // Maximum L1 norm of solution 
	int max_iterations = 10000000 ;   // Maximum number of LARS iterations.
	double lambda = 0.0 ;             // L1 weighting factor.
	
	string weight_file("") ;
	string restart_file ;

	while (1) {
		// Colons in string indicate required arguments.
		opt_type = getopt_long(argc, argv, "a:i:l:m:n:cpr:sw:h", long_options, &option_index) ;
		if ( opt_type == -1 ) break ;
		switch ( opt_type ) {
		case 'a':
			algorithm = string(optarg) ;
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

