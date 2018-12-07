/** Distributed LARS-LASSO algorithm 
		The notation and implementation closely follows
		B. Efron, T. Hastie, I. Johnstone, and R. Tibshirani, "Least Angle Regression",
    The Annals of Statistics, 32, 407-499(2004).
    Larry Fried

		Currently only works serially (non-distributed).
    11/2018
**/

// #define VERBOSE // Define for extra output

#include<math.h>
#include<iostream>
#include<fstream>
#include<string.h>
#include<getopt.h>

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
		cout.precision(6) ;
	}

	
	if ( argc < 4 ) {
		if ( RANK == 0 ) cout << "Not enough args " << endl ;
		exit(1) ;
	}
	string xname(argv[1]) ;
	string yname(argv[2]) ;
	string dname(argv[3]) ;

	static struct option long_options[] =
	{
		{"algorithm", required_argument, 0, 'a'},
		{"split_files", no_argument, 0, 's'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0}
	} ;

	int option_index = 0 ;
	int opt_type ;

	string algorithm("lasso") ;
	bool split_files = false ;
	
	while (1) {
		opt_type = getopt_long(argc, argv, "a:sh", long_options, &option_index) ;
		if ( opt_type == -1 ) break ;
		switch ( opt_type ) {
		case 'a':
			algorithm = string(optarg) ;
			break ;
		case 's':
			split_files = true ;
			break ;
		case 'h':
			display_usage(long_options) ;
			exit(0) ;
		default:
			abort() ;
		}
	}

	bool do_lasso ;
	if ( algorithm == "lars" ) {
		do_lasso = false ;
	} else if ( algorithm == "lasso" ) {
		do_lasso = true ;
	} else {
		if ( RANK == 0 ) cout << "Error: unrecognized algorithm: " << algorithm << endl ;
		exit(1) ;
	}
	
	int nprops, ndata ;
	double max_beta_norm = 1.0e+06 ;
	
	Matrix xmat ;
	if ( split_files ) {
		// Read the X matrix from multiple split files, as output by chimes_lsq
		xmat.read_split_files(xname.c_str(), dname.c_str()) ;
		nprops= xmat.dim2 ;
		ndata = xmat.dim1 ;
	} else {
		// Read the X matrix from a single file.
		ifstream xfile(xname) ;
		if ( ! xfile.is_open() ) {
			if ( RANK == 0 ) cout << "Could not open " << xname << endl ;
			exit(1) ;
		}
		ifstream dfile(dname) ;
		if ( ! dfile.is_open() ) {
			if ( RANK == 0 ) cout << "Error: could not open " << dname << endl ;
			exit(1) ;
		}
		dfile >> nprops >> ndata ;
		xmat.read(xfile, ndata, nprops, true) ;
	}
	xmat.normalize() ;
	xmat.check_norm() ;

	ifstream yfile(yname) ;
	if ( ! yfile.is_open() ) {
		if ( RANK == 0 ) cout << "Could not open " << yname << endl ;
		exit(1) ;
	}
	Vector yvec ;
	yvec.read(yfile, ndata) ;
	yvec.normalize() ;
	yvec.check_norm() ;

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

	DLARS lars(xmat, yvec) ;
	lars.do_lasso = do_lasso ;
	if ( RANK == 0 ) {
		if ( do_lasso ) {
			cout << "Using the LASSO algorithm\n" ;
		} else {
			cout << "Using the LARS algorithm\n" ;
		}
	}

	for ( int j = 0 ; ; j++ ) {
		if ( RANK == 0 ) cout << "Working on iteration " << j + 1 << endl ;
		if ( ! lars.iteration() ) 
			break ;
		if ( lars.beta.l1norm() > max_beta_norm * lars.nprops ) {
			break ;
		}
	}

	if ( RANK == 0 ) {
		cout.precision(16) ;		
		cout << "Final values:" << endl ;
		cout << "Beta: " << endl ;
		lars.beta.print() ;
	}

	lars.predict() ;
	lars.correlation() ;

	if ( RANK == 0 ) cout << "Prediction: " << endl ;
	lars.mu.print() ;

	if ( RANK == 0 ) cout << "Sq Error " << lars.sq_error() << endl ;

	// Print out the unscaled coefficients to X.txt.
	if ( RANK == 0 ) {
		ofstream xtxt("x.txt") ;
		if ( ! xtxt.is_open() ) {
			if ( RANK == 0 ) cout << "Error: could not open x.txt" << endl ;
			exit(1) ;
		}
		lars.print_unscaled(xtxt) ;

		ofstream Axfile("Ax.txt") ;
		if ( ! Axfile.is_open() ) {
			if ( RANK == 0 ) cout << "Error: could not open Ax.txt" << endl ;
			exit(1) ;
		}
		lars.print_unshifted_mu(Axfile) ;
	}
	
#ifdef VERBOSE
	if ( RANK == 0 ) cout << " Correlation: " << endl ;
	lars.c.print() ;
#endif	

}

