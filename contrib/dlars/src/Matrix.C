#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string.h>
#include <getopt.h>
#include <chrono>

#ifdef USE_MPI
#include <mpi.h>
#endif

using namespace std ;

extern int RANK ;
extern int NPROCS ;


#include "Vector.h"
#include "IntVector.h"
#include "Matrix.h"

bool Matrix::cholesky(Matrix &chol)
// Calculate the cholesky decomposition of the current matrix, and store in chol.
// Uses the upper triangular variant, A = R^T * R.	R is calculated.
{
	double eps = 1.0e-10 ;
	if ( dim1 != dim2 ) {
		cout << "Error: Cholesky decomposition only works for square matrices " << endl ;
	}
	for ( int j = 0 ; j < dim1 ; j++ ) {
		for ( int k = 0 ; k < dim1 ; k++ ) {
			if ( fabs( get(j,k) - get(k,j) ) > 1.0e-10 ) {
				cout << "Error: Cholesky decomposition only works for symmetric matrices " << endl ;
			}
		}
	}
	for ( int j = 0 ; j < dim1 ; j++ ) {
		// Diagonal elements.
		double diag ;
		diag = get(j,j) ;
		for ( int k = 0 ; k < j ; k++ ) {
			diag -= chol.get(k,j) * chol.get(k,j) ;
		}
		if ( diag < eps * eps ) 
			return false ;
		diag = sqrt(diag) ;
		chol.set(j,j,diag) ;
		for ( int k = j+1 ; k < dim1 ; k++ ) {
			double offdiag = 0.0 ;
			offdiag = get(j,k) ;
			for ( int l = 0 ; l < j ; l++ ) {
				offdiag -= chol.get(l,j) * chol.get(l,k) ;
			}
			offdiag /= chol.get(j,j) ;
			chol.set(j,k,offdiag) ;
		}

	}
	for ( int j = 0 ; j < dim1 ; j++ ) {
		for ( int k = 0 ; k < j ; k++ ) {
			chol.set(j,k,0.0) ;
		}
	}
	return true ;
}

bool Matrix::cholesky_distribute(Matrix &chol)
// Calculate the cholesky decomposition of the current matrix, and store in chol.
// Uses the upper triangular variant, A = R^T * R.	R is calculated.
// This version assumes a distributed matrix.
{
	double eps = 1.0e-10 ;
	if ( dim1 != dim2 ) {
		cout << "Error: Cholesky decomposition only works for square matrices " << endl ;
		stop_run(1) ;
	}
	if ( ! distributed ) {
		cout << "Error: cholesky_distribute requires a distributed matrix" << endl ;
		stop_run(1) ;		
	}
	if ( ! chol.distributed ) {
		cout << "Error: The cholesky output matrix was not distributed" << endl ;
		stop_run(1) ;
	}	
	if ( dim1 != chol.dim1 ) {
		cout << "Error: dimension 1 mismatch in cholesky_distribute" << endl ;
		stop_run(1) ;
	}
	if ( dim2 != chol.dim2 ) {
		cout << "Error: dimension 1 mismatch in cholesky_distribute" << endl ;
		stop_run(1) ;
	}

	// Column-oriented algorithm:	 Find one column at a time.
	for ( int j = 0 ; j < dim1 ; j++ ) {
		int rank_j = rank_from_row(j) ;

		// Diagonal elements.
		double diag = 0.0 ;
		double diag0 =	( j <= row_end && j >= row_start ) ? get(j,j) : 0.0 ;

		for ( int k = row_start ; k < j && k <= row_end ; k++ ) {
			diag0 -= chol.get(k,j) * chol.get(k,j) ;
		}
#ifdef USE_MPI
		MPI_Reduce(&diag0, &diag, 1, MPI_DOUBLE, MPI_SUM, rank_j, MPI_COMM_WORLD) ;
		//		MPI_Allreduce(&diag0, &diag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;		
#else
		diag = diag0 ;
#endif
		if ( j <= row_end && j >= row_start ) {
			if ( diag < eps * eps ) {
				cout << "Distributed Cholesky failed: negative diagonal element\n" ;
				stop_run(1) ;
			}
			diag = sqrt(diag) ;			
			chol.set(j,j,diag) ;
		}
		
		for ( int k = j+1 ; k < dim1 ; k++ ) {
			double offdiag0 = ( j <= row_end && j >= row_start ) ? get(j,k) : 0.0 ;
			double offdiag = 0.0 ;
			for ( int l = row_start ; l < j && l <= row_end ; l++ ) {
				offdiag0 -= chol.get(l,j) * chol.get(l,k) ;
			}
#ifdef USE_MPI
			MPI_Reduce(&offdiag0, &offdiag, 1, MPI_DOUBLE, MPI_SUM, rank_j, MPI_COMM_WORLD) ;
			//			MPI_Allreduce(&offdiag0, &offdiag, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;			
#else
			offdiag = offdiag0 ;
#endif
			if ( j <= row_end && j >= row_start ) {
				offdiag /= chol.get(j,j) ;
				chol.set(j,k,offdiag) ;
			}
		}
	}
	for ( int j = row_start ; j <= row_end ; j++ ) {
		for ( int k = 0 ; k < j ; k++ ) {
			chol.set(j,k,0.0) ;
		}
	}
	return true ;
}


bool Matrix::cholesky_add_row(const Matrix &chol0, const Vector &newr)
	// Given cholesky decomposition of a matrix in chol0
	// find the cholesky decomposition of a new matrix formed by adding
	// the specified row (newr) to the original matrix (not given).
	// The new row is placed in the last index.
{
	if ( dim1 != dim2 ) {
		cout << "Error: Cholesky decomposition only works for square matrices " << endl ;
		stop_run(1) ;		
	}
	if ( newr.dim != dim1 || chol0.dim1 != dim1 - 1 || chol0.dim1 != chol0.dim2 ) {
		cout << "Dimension mismatch" ;
		stop_run(1) ;
	}
	Vector xtmp(dim1) ;
	double cdotc = 0.0 ;
	for ( int j = 0 ; j < dim1 - 1 ; j++ ) {
		double sum = newr.get(j) ;
#ifdef USE_OPENMP		
#pragma omp parallel for shared(chol0,xtmp,j) reduction(+:sum) default(none)
#endif						
		for ( int k = 0 ; k < j ; k++ ) {
			sum -= chol0.get(k,j) * xtmp.get(k) ;
		}
		xtmp.set(j, sum / chol0.get(j,j) );
		cdotc += xtmp.get(j) * xtmp.get(j) ;
	}
	double d = newr.get(dim1-1) - cdotc ;
	if ( d >= 0.0 ) {
		xtmp.set(dim1-1, sqrt(d) ) ;
	} else {
		cout << "Could not add row to Cholesky matrix " << endl ;
		return false ;
	}

#ifdef USE_OPENMP		
#pragma omp parallel for shared(chol0) default(none)
#endif					
	for ( int j = 0 ; j < dim1 - 1 ; j++ ) {
		for ( int k = 0 ; k <= j ; k++ ) {
			set(k,j, chol0.get(k,j)) ;
		}
		for ( int k = j + 1 ; k < dim1 ; k++ ) {
			set(k,j, 0.0) ;
		}
	}

#ifdef USE_OPENMP		
#pragma omp parallel for shared(xtmp) default(none)
#endif					
	for ( int j = 0 ; j < dim1 ; j++ ) {
		set(j, dim1-1, xtmp.get(j)) ;
	}
			
#ifdef USE_OPENMP		
#pragma omp parallel for default(none)
#endif					
	for ( int j = 0 ; j < dim1 - 1 ; j++ ) {
		set(dim1-1, j, 0.0) ;
	}
	return true ;
}
	
bool Matrix::cholesky_remove_row(int id ) 
	// Based on larscpp github package.	 Modified so that diagonal elements of the cholesky
	// matrix are always positive.
	// Remove a row from a Cholesky matrix.	 Update the dimension of the matrix
	// Use transpose matrix access functions because original code was written in terms of a lower
	// triangular matrix, whereas we use an upper triangular matrix.
{
	double a(0),b(0),c(0),s(0),tau(0);
	if ( dim1 != dim2 ) {
		cout << "Error: the cholesky matrix should be square" << endl ;
		stop_run(1) ;
	}
	int lth = dim1 - 1 ;

	for(int i=id; i < lth; ++i){

#ifdef USE_OPENMP		
#pragma omp parallel for shared(i) default(none)
#endif					
		for( int j=0; j<i; ++j) {
			setT(i,j, getT(i+1,j)) ;
		}

		a = getT(i+1,i);
		b = getT(i+1,i+1);
		if(b==0){
			setT(i,i,a) ;
			continue;
		}
		if( fabs(b) > fabs(a) ){
			tau = -a/b;
			s = 1/sqrt(1.0 + tau*tau);
			c = s*tau;
		} else {
			tau = -b/a;
			c = 1/sqrt(1.0 + tau*tau);
			s = c * tau;
		}
		if ( c*a - s*b > 0.0 ) {
			setT(i,i, c*a - s*b) ;
		} else {
			setT(i,i, s*b - c*a) ;
			s = -s ;
			c = -c ;
		}
		// L(i,i+1) = s*a + c*b;

#ifdef USE_OPENMP		
#pragma omp parallel for shared(i,lth,c,s) private(a,b) default(none)
#endif					
		for( int j=i+2; j<=lth; ++j){
			a = getT(j,i);
			b = getT(j, i+1);
			setT(j, i, c*a - s*b) ;
			setT(j, i+1, s*a + c*b) ;
		}
	}

#ifdef USE_OPENMP		
#pragma omp parallel for shared(lth) default(none)
#endif				
	for( int i=0; i<= lth; ++i)
		setT(lth, i, 0) ;

	// Shift the storage to accommodate a change of dimension.
	double *newmat = new double[dim1*dim2] ;



#ifdef USE_OPENMP
	//#pragma omp parallel for shared(lth,newmat,mat) default(none)
#pragma omp parallel for shared(lth,newmat) default(none)
#endif				
	for ( int i = 0 ; i < lth ; i++ ) {
		for ( int j = 0 ; j < lth ; j++ ) {
			newmat[i * (dim2-1) + j] = mat[i * dim2 + j] ;
		}
	}
	delete [] mat ;
	mat = newmat ;
		
	dim1-- ;
	dim2-- ;
	num_rows-- ;
	row_end-- ;
		
	return true ;
}
	
void Matrix::cholesky_sub(Vector &x, const Vector &b)
// Assuming that the current matrix holds a Cholesky decomposition
// in R^T * R form, calculate the solution vector x given the vector
// of input values b.
{
	if ( b.dim != x.dim || b.dim != dim1 || b.dim != dim2 ) {
		cout << "Dimension mismatch" ;
		stop_run(1) ;
	}

	double condition ;
	double diag_min = 1.0e50 ;
	double diag_max = 0.0 ;

	// Print an estimate of the condition number.
	for ( int j = 0 ; j < dim1 ; j++ ) {
		if ( get(j,j) > diag_max ) {
			diag_max = get(j,j) ;
		}
		if ( get(j,j) < diag_min ) {
			diag_min = get(j,j) ;
		}
		if ( get(j,j) < 0.0 ) {
			cout << "Warning: Negative diagonal of cholesky decomp found!\n" ;
		}
	}
	if ( diag_min > 0.0 ) {
		condition = diag_max * diag_max / ( diag_min * diag_min) ;
	} else {
		condition = 1.0e50 ;
	}

	int prec = cout.precision() ;
	cout << "Cholesky estimate of condition number = " <<
		std::scientific << std::setprecision(3) << condition << endl ;
	cout.precision(prec) ;
												
	Vector xtmp(dim1) ;
	Matrix& mat=*this ;
			
	for ( int j = 0 ; j < dim1 ; j++ ) {
		double sum = b.get(j) ;

#ifdef USE_OPENMP		
#pragma omp parallel for shared(xtmp,mat,j) reduction(+:sum) default(none)
#endif						
		for ( int k = 0 ; k < j ; k++ ) {
			sum -= mat.get(k,j) * xtmp.get(k) ;
		}

		xtmp.set(j, sum / mat.get(j,j) );
	}

	for ( int j = dim1 - 1 ; j >= 0 ; j-- ) {
		double sum = xtmp.get(j) ;

#ifdef USE_OPENMP		
#pragma omp parallel for shared(j,x) reduction(+:sum) default(none)
#endif						
		for ( int k = j + 1 ; k < dim1 ; k++ ) {
			sum -= get(j,k) * x.get(k) ;
		}

		x.set(j, sum / get(j,j) ) ;
	}
}

void Matrix::cholesky_sub_distribute(Vector &x, const Vector &b)
// Assuming that the current matrix holds a Cholesky decomposition
// in R^T * R form, calculate the solution vector x given the vector
// of input values b.
// The matrix is assumed to be distributed.
{
	if ( b.dim != x.dim || b.dim != dim1 || b.dim != dim2 ) {
		cout << "Dimension mismatch" ;
		stop_run(1) ;
	}
	if ( ! distributed ) {
		cout << "Error: cholesky_sub_distribute requires a distributed matrix" << endl ;
	}

	double condition ;
	double diag_min0 = 1.0e50 ;
	double diag_max0 = 0.0 ;
	double diag_min = 0.0 ;
	double diag_max = 0.0 ;
	// Print an estimate of the condition number.
	for ( int j = row_start ; j <= row_end ; j++ ) {
		if ( get(j,j) > diag_max0 ) {
			diag_max0 = get(j,j) ;
		}
		if ( get(j,j) < diag_min0 ) {
			diag_min0 = get(j,j) ;
		}
		if ( get(j,j) < 0.0 ) {
			cout << "Warning: Negative diagonal of cholesky decomp found!\n" ;
		}
	}
#ifdef USE_MPI
		MPI_Allreduce(&diag_min0, &diag_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD) ;
		MPI_Allreduce(&diag_max0, &diag_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD) ;		
#else
		diag_min = diag_min0 ;
		diag_max = diag_max0 ;
#endif
	
	if ( diag_min > 0.0 ) {
		condition = diag_max * diag_max / ( diag_min * diag_min) ;
	} else {
		condition = 1.0e50 ;
	}

	if ( RANK == 0 ) {
		int prec = cout.precision() ;
		cout << "Cholesky estimate of condition number = " <<
			std::scientific << std::setprecision(3) << condition << endl ;
		cout.precision(prec) ;
	}
												
	Vector xtmp(dim1) ;
	Matrix& mat=*this ;
			
	for ( int j = 0 ; j < dim1 ; j++ ) {
		double sum0 = ( row_start <= j && j <= row_end ) ? b.get(j) : 0.0 ;

		int end = min(j-1, row_end) ;
		
#ifdef USE_OPENMP		
#pragma omp parallel for shared(xtmp,mat,j,end) reduction(+:sum0) default(none)
#endif						
		for ( int k = row_start ; k <= end ; k++ ) {
			sum0 -= mat.get(k,j) * xtmp.get(k) ;
		}

		double sum = 0.0 ;
#ifdef USE_MPI
		MPI_Allreduce(&sum0, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else
		sum = sum0 ;
#endif

		if ( row_start <= j && j <= row_end ) {
			xtmp.set(j, sum / mat.get(j,j) );
		}
	}

	// This serializes with a row-distributed matrix.
	for ( int j = dim1 - 1 ; j >= 0 ; j-- ) {

		double sum = 0.0 ;
		if ( row_start <= j && j <= row_end ) {

			sum = xtmp.get(j) ;

#ifdef USE_OPENMP		
#pragma omp parallel for shared(j,x) reduction(+:sum) default(none)
#endif
			for ( int k = j + 1 ; k < dim1 ; k++ ) {
				sum -= get(j,k) * x.get(k) ;
			}
			sum /= get(j,j) ;
		}
		double sum1 = 0.0 ;
#ifdef USE_MPI
		int root ;
		root = rank_from_row(j) ;
		MPI_Allreduce(&sum, &sum1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else
		sum1 = sum ;
#endif
		x.set(j, sum1 ) ;
	}

}


bool Matrix::cholesky_add_row_distribute(const Matrix &chol0, const Vector &newr)
	// Given cholesky decomposition of a matrix in chol0
	// find the cholesky decomposition of a new matrix formed by adding
	// the specified row (newr) to the original matrix (not given).
	// The new row is placed in the last index.
{
	if ( dim1 != dim2 ) {
		cout << "Error: Cholesky decomposition only works for square matrices " << endl ;
		stop_run(1) ;
	}
	if ( ! distributed ) {
		cout << "Error: cholesky_add_row_distribute requires a distributed matrix" << endl ;
		stop_run(1) ;		
	}
	if ( ! chol0.distributed ) {
		cout << "Error: cholesky_add_row_distribute requires a distributed input matrix" << endl ;
		stop_run(1) ;		
	}	
	if ( newr.dim != dim1 || chol0.dim1 != dim1 - 1 || chol0.dim1 != chol0.dim2 ) {
		cout << "Dimension mismatch" ;
		stop_run(1) ;
	}
	
	Vector xtmp(dim1) ;
	Vector chol0_diag_tmp(dim1,0.0) ;
	Vector chol0_diag(dim1,0.0) ;	

	if ( dim1-1 > 0 ) {
		for ( int j = chol0.row_start ; j <= chol0.row_end ; j++ ) {
			chol0_diag_tmp.set(j, chol0.get(j,j)) ;
		}
#ifdef USE_MPI
		MPI_Allreduce(chol0_diag_tmp.vec, chol0_diag.vec, dim1-1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else
		for ( int j = 0 ; j < dim1-1 ; j++ ) {
			chol0_diag.set(j, chol0_diag_tmp.get(j)) ;
		}
#endif
	}
	
	double cdotc = 0.0 ;
	for ( int j = 0 ; j < dim1 - 1 ; j++ ) {
		double sum0 = ( row_start <= j && j <= row_end ) ? newr.get(j) : 0.0 ;
		double sum = 0.0 ;
		
		for ( int k = chol0.row_start ; k < j && k <= chol0.row_end ; k++ ) {
			sum0 -= chol0.get(k,j) * xtmp.get(k) ;
		}
#ifdef USE_MPI
		MPI_Allreduce(&sum0, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;		
#else
		sum = sum0 ;
#endif		
		xtmp.set(j, sum / chol0_diag.get(j) );
		cdotc += xtmp.get(j) * xtmp.get(j) ;
	}
	double d = newr.get(dim1-1) - cdotc ;
	if ( d >= 0.0 ) {
		xtmp.set(dim1-1, sqrt(d) ) ;
	} else {
		cout << "Could not add row to distributed Cholesky matrix " << endl ;
		return false ;
	}

	for ( int k = 0	 ; k < dim1 - 1 ; k++ ) {
		int chol0_rank_k = chol0.rank_from_row(k) ;
		int chol_rank_k	 = rank_from_row(k) ;
			
		if ( chol0_rank_k == chol_rank_k ) {
			// Data resides in the same rank. 
			if ( row_start <= k && k <= row_end ) {
				for ( int j = k ; j < dim1 - 1 ; j++ ) {
					set(k,j, chol0.get(k,j)) ;
				}
			}
		} else {
			Vector val(dim1-1,0.0) ;
			if ( RANK == chol0_rank_k ) {
				for ( int j = k ; j < dim1 - 1 ; j++ ) {					
					val.set(j,chol0.get(k,j)) ;
				}
#ifdef USE_MPI
				MPI_Send(val.vec, dim1-1, MPI_DOUBLE, chol_rank_k, 0, MPI_COMM_WORLD) ;
#endif
			} else if ( RANK == chol_rank_k ) {
#ifdef USE_MPI					
				MPI_Recv(val.vec, dim1-1, MPI_DOUBLE, chol0_rank_k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
				if ( row_start <= k && k <= row_end ) {
					for ( int j = k ; j < dim1 - 1 ; j++ ) {											
						set(k, j, val.get(j)) ;
					}
				}
#endif
			}
			//MPI_Barrier(MPI_COMM_WORLD) ;
		}
	}

	for ( int k = 1 ; k < dim1 ; k++ ) {	
		if ( row_start <= k && k <= row_end ) {
			for ( int j = k - 1 ; j >= 0 ; j-- ) {					
				set(k,j, 0.0) ;
			}
		}
	}

	for ( int j = row_start	 ; j <= row_end && j < dim1 ; j++ ) {
		set(j, dim1-1, xtmp.get(j)) ;
	}

	int dim11 = dim1-1 ;
	if ( row_start <= dim11 && dim11 <= row_end ) {
		for ( int j = 0 ; j < dim11 ; j++ ) {
			set(dim11, j, 0.0) ;
		}
	}
	
	return true ;
}

bool Matrix::cholesky_remove_row_dist(int id ) 
	// Based on larscpp github package.	 Modified so that diagonal elements of the cholesky
	// matrix are always positive.
	// Remove a row from a Cholesky matrix.	 Update the dimension of the matrix
	// This version assumes a distributed matrix.
{
	double a(0),b(0),c(0),s(0),tau(0);
	if ( dim1 != dim2 ) {
		cout << "Error: the cholesky matrix should be square" << endl ;
		stop_run(1) ;
	}
	if ( ! distributed ) {
		cout << "Error: cholesky_remove_dist assumes a distributed matrix\n" ;
		stop_run(1) ;
	}
	
	int lth = dim1 - 1 ;
	Vector atmp(dim1), btmp(dim1) ;

#ifdef USE_MPI	
	int rtag = 1, stag = 2 ;
#endif
	
	for(int i=id; i < lth; ++i) {

		for( int j= row_start ; j<i && j <= row_end ; ++j) {
			set(j,i, get(j,i+1)) ;
		}

		int rank_i = rank_from_row(i) ;
		int rank_i_1 = rank_from_row(i+1) ;

		if ( RANK == rank_i ) {
			a = get(i, i+1);
//			if ( rank_i != rank_i_1 ) {
// #ifdef USE_MPI				
//				MPI_Send(&a, 1, MPI_DOUBLE, rank_i_1, 1, MPI_COMM_WORLD) ;
// #endif				
//			}
//		} else if ( RANK == rank_i_1 ) {
//			if ( rank_i != rank_i_1 ) {
// #ifdef USE_MPI				
//				MPI_Recv(&a, 1, MPI_DOUBLE, rank_i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
// #endif				
//			}
		} else {
			a = 0.0 ;
		}
		// #ifdef USE_MPI		
		//		MPI_Bcast(&a, 1, MPI_DOUBLE, rank_i, MPI_COMM_WORLD) ;
		// #endif		

		if ( RANK == rank_i_1 ) {
			b = get(i+1,i+1);
//			if ( rank_i != rank_i_1 ) {
// #ifdef USE_MPI				
//				MPI_Send(&b, 1, MPI_DOUBLE, rank_i, 0, MPI_COMM_WORLD) ;
// #endif				
//			}
//		} else if ( RANK == rank_i ) {
//			if ( rank_i != rank_i_1 ) {
// #ifdef USE_MPI				
//				MPI_Recv(&b, 1, MPI_DOUBLE, rank_i_1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
// #endif				
//			}
		} else {
			b = 0.0 ;
		}
		// #ifdef USE_MPI		
		//		MPI_Bcast(&b, 1, MPI_DOUBLE, rank_i_1, MPI_COMM_WORLD) ;
		// #endif		

#ifdef USE_MPI
		if ( rank_i != rank_i_1 ) {
			// Exchange a and b between ranks.
			if ( RANK == rank_i ) {
				MPI_Sendrecv(&a, 1, MPI_DOUBLE, rank_i_1, stag,
										 &b, 1, MPI_DOUBLE, rank_i_1, rtag,
										 MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
			} else if ( RANK == rank_i_1 ) {
				MPI_Sendrecv(&b, 1, MPI_DOUBLE, rank_i, rtag,
										 &a, 1, MPI_DOUBLE, rank_i, stag,
										 MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
			}
		}
#endif		
		if( RANK == rank_i && b == 0.0 ) {
			set(i,i,a) ;
			continue;
		}

		if ( RANK == rank_i || RANK == rank_i_1 ) {
			if( fabs(b) > fabs(a) ){
				tau = -a/b;
				s = 1/sqrt(1.0 + tau*tau);
				c = s*tau;
			} else {
				tau = -b/a;
				c = 1/sqrt(1.0 + tau*tau);
				s = c * tau;
			}
			if ( c*a - s*b > 0.0 ) {
				if ( RANK == rank_i ) {
					set(i,i, c*a - s*b) ;
				}
			} else {
				if ( RANK == rank_i ) {
					set(i,i, s*b - c*a) ;
				}
				s = -s ;
				c = -c ;
			}
			// L(i,i+1) = s*a + c*b;

			for( int j=i+2; j<=lth; ++j){
				if ( RANK == rank_i ) {
					atmp.set(j, get(i, j)) ;
				}
				if ( RANK == rank_i_1 ) {
					btmp.set(j, get(i+1, j)) ;
				}
			}
			
			if ( RANK == rank_i && RANK == rank_i_1 ) {
				for( int j=i+2; j<=lth; ++j){					
					set(i, j, c*atmp.get(j) - s*btmp.get(j)) ;
					set(i+1, j, s*atmp.get(j) + c*btmp.get(j)) ;
				}
			} else if ( RANK == rank_i ) {
				// Sends and receives need to be matched.
#ifdef USE_MPI				
				MPI_Sendrecv(atmp.vec, dim1, MPI_DOUBLE, rank_i_1, stag,
										 btmp.vec, dim1, MPI_DOUBLE, rank_i_1, rtag,
										 MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
#endif				
				for( int j=i+2; j<=lth; ++j){					
					set(i, j, c*atmp.get(j) - s*btmp.get(j)) ;
				}					
			} else if ( RANK == rank_i_1 ) {
				// Sends and receives need to be matched.
#ifdef USE_MPI				
				MPI_Sendrecv(btmp.vec, dim1, MPI_DOUBLE, rank_i, rtag,
										 atmp.vec, dim1, MPI_DOUBLE, rank_i, stag,
										 MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
#endif				
				for( int j=i+2; j<=lth; ++j){										
					set(i+1, j, s*atmp.get(j) + c*btmp.get(j)) ;
				}
			}
		}
	}

	for( int i= row_start; i<= lth && i <= row_end ; ++i)
		set(i, lth, 0) ;

	Matrix newmat(lth, lth, true) ;

	for ( int i = 0 ; i < lth ; i++ ) {
		int rank_i_new = newmat.rank_from_row(i) ;
		int rank_i = rank_from_row(i) ;

		if ( rank_i == rank_i_new ) {
			if ( newmat.row_start <= i && i <= newmat.row_end ) {
				for ( int j = 0 ; j < lth ; j++ ) {
					newmat.set(i, j, get(i,j)) ;
				}
			}
		} else {
			Vector tmp(lth,0.0) ;
			if ( RANK == rank_i ) {
				for ( int l = 0 ; l < lth ; l++ ) {
					tmp.set( l, get(i, l) ) ;
				}
				MPI_Send(tmp.vec, lth, MPI_DOUBLE, rank_i_new, 0, MPI_COMM_WORLD) ;
			}
			else if ( RANK == rank_i_new ) {
				MPI_Recv(tmp.vec, lth, MPI_DOUBLE, rank_i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) ;
				for ( int l = 0 ; l < lth ; l++ ) {
					newmat.set( i, l, tmp.get(l) ) ;
				}
			}
		}
	}
	delete [] mat ;
	mat = newmat.mat ;
	newmat.mat = nullptr ;
		
	dim1			= newmat.dim1 ;
	dim2			= newmat.dim2 ;
	num_rows	= newmat.num_rows ;
	row_end		= newmat.row_end ;
	row_start = newmat.row_start ;
	
	return true ;
}
