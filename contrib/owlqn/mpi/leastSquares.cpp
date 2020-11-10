#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "mpiexterns.h"
#include "leastSquares.h"

LeastSquaresProblem::LeastSquaresProblem(const char* matFilename, const char* bFilename, bool split_files) {

	if ( ! split_files ) {
		// Read in matrix market format.
		ifstream matfile(matFilename);
		if (!matfile.good()) {
			cerr << "error opening matrix file " << matFilename << endl;
			exit_run(1);
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
			exit_run(1);
		}
		
		ifstream bFile(bFilename);
		if (!bFile.good()) {
			cerr << "error opening y-value file " << bFilename << endl;
			exit_run(1);
		}

		getline(bFile, s);
		if (s.compare("%%MatrixMarket matrix array real general")) {
			bFile.close();
			cerr << "unsupported y-value file format \"" << s << "\" in " << bFilename << endl;
			exit_run(1);
		}

		skipEmptyAndComment(bFile, s);
		stringstream bst(s);
		size_t bNum, bCol;
		bst >> bNum >> bCol;
		if (bNum != m) {
			cerr << "number of y-values doesn't match number of instances in " << bFilename << endl;
			exit_run(1);
		} else if (bCol != 1) {
			cerr << "y-value matrix may not have more than one column" << endl;
			exit_run(1);
		}

		b.resize(m);
		for (size_t i=0; i<m; i++) {
			double val;
			bFile >> val;
			b[i] = val;
		}
		bFile.close();
	}
	else {
		read_split_files(matFilename, bFilename) ;
	}
}

void LeastSquaresProblem::read_split_files(const char* matFilename, const char* bFilename)
// Read split file output from chimes_lsq.
{
	ifstream dim_file ;
	char name[80] ;

	// Find out how many files were written by chimes_lsq.
	int total_files = 0 ;
	ifstream test_file ;
	for ( int j = 0 ; j < NPROCS + 1 ; j++ ) {
		memset(name, 0, 80) ;
		sprintf(name, "dim.%04d.txt", j) ;
		test_file.open(name) ;
		if ( test_file.is_open() ) {
			++total_files ;
			test_file.close() ;
		} else {
			test_file.close() ;
			break ;
		}
	}
	if ( total_files > NPROCS ) {
		ErrorMsg("Not enough processes specified") ;
	}

	// Somewhat tricky logic to allow the number of processors to be greater than the number of files.
	// In that case, only some of the rows in the split A matrix are used.
	int proc_fac = NPROCS / total_files ;
	int rank_div = RANK / proc_fac ;
	int my_file = ( rank_div > total_files - 1 ) ? total_files - 1 : rank_div ;
	int my_offset = RANK % proc_fac ;

		
	cout << "RANK = " << RANK << " proc_fac = " << proc_fac << " my_file = " << my_file << " my_offset = " << my_offset << endl ;

	// Read matrix dimensions from the dim.*.txt file.
	sprintf(name, "dim.%04d.txt", my_file) ;
	dim_file.open(name) ;
	if ( ! dim_file.is_open() ) {
		cerr << "Could not open " + string(name) + "\n" ;
		exit_run(1) ;
	}
			
	int mdim, mdim2, ndim ;
		
	// Dimensions to use if NPROCS == total_files
	int mstart0, mend0, mstore0 ;
		
	dim_file >> n >> mstart0 >> mend0 >> m ;
	if ( ! dim_file.good() ) {
		cerr << "Error reading dim file\n" ;
		exit_run(1) ;
	}
	dim_file.close() ;

	if ( mstart0 <= mend0 ) {
		// mstore0 is the number of rows in the file.
		mstore0 = (mend0 - mstart0 + 1) ;
	} else {
		mstore0 = 0 ;
	}

	mstart = mstart0 + my_offset * (mstore0 / proc_fac) ;

	if ( rank_div < total_files ) {
		if ( my_offset == proc_fac - 1 )
			mend = mend0 ;
		else
			mend = mstart + (mstore0 / proc_fac) - 1 ;
		mstore = mend - mstart + 1 ;
	} else {
		mend = mstart - 1 ;
		mstore = 0 ;
	}

	cout << "RANK = " << RANK << " mstart = " << mstart << " mend = " << mend << " mstore = " << mstore << endl ;
		
	// Append the processor number to the A matrix name.
	char matFilename2[80] ;
	string str_filename(matFilename) ;
	std::size_t found = str_filename.find(".") ;
	if ( found == string::npos ) {
		cerr < "A matrix file name must end with a suffix" ;
		exit_run(1) ;
	}
		
	str_filename = str_filename.substr(0,found+1) ;
	sprintf(matFilename2, "%s%04d.txt", str_filename.data(), my_file) ;
	ifstream matfile(matFilename2);
	if (!matfile.good()) {
		cerr << "error opening matrix file " << matFilename << endl;
		exit_run(1);
	}

	Amat.resize(mstore * n);

	for (size_t i= mstart0 ; i <= mend0 ; i++) {
		for (size_t j=0; j<n; j++) {
			double val;
			matfile >> val;
			if ( i >= mstart && i <= mend ) 
				A(i, j) = val;
		}
	}
	if ( ! matfile.good() ) {
		cerr << "Error reading A matrix" ;
		exit_run(1) ;
	}
	matfile.close();

	// Open and read the b vector.
	ifstream bFile(bFilename);
	if (!bFile.good()) {
		cerr << "error opening y-value file " << bFilename << endl;
		exit_run(1);
	}
	b.resize(m);
	for (size_t i=0; i<m; i++) {
		double val;
		bFile >> val;
		b[i] = val;
		if ( ! bFile.good() ) {
			cerr << "Error reading b file" ;
			exit_run(1) ;
		}
	}
	bFile.close();
}


double LeastSquaresObjective::Eval(const DblVec& input, DblVec& gradient)
{

	DblVec temp(problem.m,0.0);

	if (input.size() != problem.n) {
		cerr << "Error: input is not the correct size." << endl;
		exit_run(1);
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



void LeastSquaresObjective::Eval_Ax(const DblVec& input, DblVec& Ax)
// Evaluate A * x
{

	DblVec temp(problem.m,0.0);

	if (input.size() != problem.n) {
		cerr << "Error: input is not the correct size." << endl;
		exit_run(1);
	}
	if ( Ax.size() != problem.m ) {
		cerr << "Error: output is not the correct size\n" ;
		exit_run(1) ;
	}
	
	double value = 0.0;
	for (size_t j=0; j<problem.n; j++) {
		for (size_t i= problem.mstart; i <= problem.mend; i++) {
			temp[i] += input[j] * problem.A(i,j);
		}
	}

#ifdef USE_MPI
	DblVec send(problem.m,0.0);
	DblVec recv(problem.m,0.0);
	for ( size_t j = 0 ; j < problem.m ; j++ ) {
		send[j] = temp[j] ;
	}
	MPI_Allreduce(send.data(), recv.data(), problem.m,
								MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
	for ( size_t j = 0 ; j < problem.m ; j++ ) {
		Ax[j] = recv[j] ;
	}
#else
	for ( size_t j = 0 ; j < problem.m ; j++ ) {
		Ax[j] = temp[j] ;
	}
#endif

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


void skipEmptyAndComment(std::ifstream& file, std::string& s)
{
	do {
		std::getline(file, s);
	} while (s.size() == 0 || s[0] == '%');
}
