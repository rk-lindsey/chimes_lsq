#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "mpiexterns.h"
#include "leastSquares.h"

LeastSquaresProblem::LeastSquaresProblem(const char* matFilename, const char* bFilename) {
	ifstream matfile(matFilename);
	if (!matfile.good()) {
		cerr << "error opening matrix file " << matFilename << endl;
		exit(1);
	}

	string s;
	getline(matfile, s);
	if (!s.compare("%%MatrixMarket matrix array real general")) {
		skipEmptyAndComment(matfile, s);
		stringstream st(s);
		st >> m >> n;

		calc_storage() ;
		Amat.resize(mstore * n);

		for (size_t j=0; j<n; j++) {
			for (size_t i=0; i<m; i++) {
				double val;
				matfile >> val;
				if ( i >= mstart && i <= mend ) 
					A(i, j) = val;
			}
		}

		matfile.close();
	} else {
		matfile.close();
		cerr << "Unsupported matrix format \"" << s << "\" in " << matFilename << endl;
		exit(1);
	}

	ifstream bFile(bFilename);
	if (!bFile.good()) {
		cerr << "error opening y-value file " << bFilename << endl;
		exit(1);
	}
	getline(bFile, s);
	if (s.compare("%%MatrixMarket matrix array real general")) {
		bFile.close();
		cerr << "unsupported y-value file format \"" << s << "\" in " << bFilename << endl;
		exit(1);
	}

	skipEmptyAndComment(bFile, s);
	stringstream bst(s);
	size_t bNum, bCol;
	bst >> bNum >> bCol;
	if (bNum != m) {
		cerr << "number of y-values doesn't match number of instances in " << bFilename << endl;
		exit(1);
	} else if (bCol != 1) {
		cerr << "y-value matrix may not have more than one column" << endl;
		exit(1);
	}

	b.resize(m);
	for (size_t i=0; i<m; i++) {
		double val;
		bFile >> val;
		b[i] = val;
	}
	bFile.close();
}


double LeastSquaresObjective::Eval(const DblVec& input, DblVec& gradient)
{

	DblVec temp(problem.m,0.0);

	if (input.size() != problem.n) {
		cerr << "Error: input is not the correct size." << endl;
		exit(1);
	}

	for (size_t i= problem.mstart; i <= problem.mend; i++) {
		temp[i] = -problem.b[i];
	}

	double value = 0.0;
	for (size_t j=0; j<problem.n; j++) {
		if ( RANK == 0 ) {
			value += input[j] * input[j] * l2weight;
			gradient[j] = l2weight * input[j];
		} else {
			gradient[j] = 0.0 ;
		}
		for (size_t i= problem.mstart; i <= problem.mend; i++) {
			temp[i] += input[j] * problem.A(i,j);
		}
	}

	for (size_t i= problem.mstart ; i<= problem.mend; i++) {
		if (temp[i] == 0.0) continue;

		value += temp[i] * temp[i];
		for (size_t j=0; j<problem.n; j++) {
			gradient[j] += problem.A(i, j) * temp[i];
		}
	}

#ifdef USE_MPI
	DblVec send(problem.n+1,0.0);
	DblVec recv(problem.n+1,0.0);
	for ( size_t j = 0 ; j < problem.n ; j++ ) {
		send[j] = gradient[j] ;
	}
	send[problem.n] = value ;
	MPI_Allreduce(send.data(), recv.data(), problem.n+1,
								MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	for ( size_t j = 0 ; j < problem.n ; j++ ) {
		gradient[j] = recv[j] ;
	}
	value = recv[problem.n] ;
#endif

	// Scale down the gradient to make problem definition the
	// same as is used in SciKit-Learn.
	for ( size_t j = 0 ; j < problem.n ; j++ ) {
		gradient[j] /= (double) problem.m ;
	}
	return (0.5 * value / problem.m) + 1.0;
}

void LeastSquaresProblem::calc_storage()
{
	// Calculate which rows of the A matrix should be store by this
	// process.
	if ( RANK != NPROCS - 1 ) {
		mstore = m / NPROCS ;
		mstart = RANK * mstore ;
		mend   = (RANK+1) * mstore - 1 ;
	}
	else {
		mstore = m / NPROCS ;
		mstart = RANK * mstore ;
		mend = m - 1 ;
		mstore += m % NPROCS ;
	}
}
