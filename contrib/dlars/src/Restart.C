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


extern int RANK ;
extern int NPROCS ;


using namespace std ;

#include "Vector.h"
#include "IntVector.h"
#include "Matrix.h"
#include "DLARS.h"
#include "Restart.h"

void Restart::read_vector(ifstream &inf, string name, Vector &vec)
{
	string line ;
	getline(inf,line) ;

#ifdef VERBOSE				
	if ( RANK == 0 ) cout << name << " line is: " << line << endl ;
#endif
	
	if ( line.find(name) != string::npos ) {
		vec.read_sparse(inf) ;
	} else {
		if ( RANK == 0 ) cout << "Could not find " << name << " in restart file\n" ;
		stop_run(1) ;
	}
}

void Restart::read_int_vector(ifstream &inf, string name, IntVector &vec)
{
	string line ;
	getline(inf,line) ;

#ifdef VERBOSE				
	if ( RANK == 0 ) cout << name << " line is: " << line << endl ;
#endif
	
	if ( line.find(name) != string::npos ) {
		vec.read_sparse(inf) ;
	} else {
		if ( RANK == 0 ) cout << "Could not find " << name << " in restart file\n" ;
		stop_run(1) ;
	}	
}

void Restart::read_matrix(ifstream &inf, string name, Matrix &mat, int dim1,
						 int dim2, bool distributed)
{
	string line ;
	
	getline(inf,line) ;

#ifdef VERBOSE				
	if ( RANK == 0 ) cout << "RANK " << RANK << name << " line is " << line << endl ;
#endif
	
	if ( line.find(name) != string::npos ) {
		if ( distributed ) {
			mat.read(inf, dim1, dim2, true, true) ;
#ifdef VERBOSE			
			if ( RANK == 0 ) cout << name << " = " << endl ;
			mat.print() ;
#endif			
		} else if ( RANK == 0 ) {
			mat.read(inf, dim1, dim2, false, false) ;
#ifdef VERBOSE						
			if ( RANK == 0 ) cout << name << " = " << endl ;
			mat.print() ;
#endif						
		}
	} else {
		cout << "Could not find " << name << " in restart file\n" ;
		cout << "Read the following line: " << line << endl ;
		stop_run(1) ;
	}

	// Read terminal newline.
	getline(inf,line) ;
}
