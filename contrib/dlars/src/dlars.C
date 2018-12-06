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


#ifdef USE_MPI
#include <mpi.h>
#endif

int RANK ;
int NPROCS ;

using namespace std ;

class Vector {
public:
	double *vec ;
	int dim ;
	double shift ;

	Vector(int d1)
		{
			dim = d1 ;
			vec = new double[d1] ;
			shift = 0 ;
		}
	Vector(int d1, double val)
		{
			dim = d1 ;
			vec = new double[d1] ;
			shift = 0 ;
			for ( int j = 0 ; j < dim ; j++ ) {
				vec[j] = val ;
			}
		}
	Vector()
		{
			dim = 0 ;
			vec = NULL ;
			shift = 0 ;
		}
	~Vector()
		{
			
				delete[] vec ;
		}
	int size() const {
		return dim ;
	}
	void set(int i, double val) {
#ifdef DEBUG					
		if ( i >= dim ) {
			cout << "Vector set out of bounds" << endl ;
			exit(1) ;
		}
#endif					
		vec[i] = val ;
	}
	void realloc(int size) {
		// Reallocate the vector.
		
			delete [] vec ;
		vec = new double[size] ;
		dim = size ;
	}
	
	void read(ifstream &file, int dim0) 
		{
			dim = dim0 ;
			vec = new double[dim] ;
			for ( int i = 0 ; i < dim ; i++ ) {
				double val ;
				file >> val ;
				set(i, val) ;
			}
		}
	void normalize()
		{
			shift = 0 ;
			for ( int i = 0 ; i < dim ; i++ ) {
				shift += vec[i] ;
			}
			shift /= dim;
			for ( int i = 0 ; i < dim ; i++ ) {
				vec[i] -= shift ;
			}
		}
	void check_norm()
		{
			double test = 0.0 ;
			for ( int i = 0 ; i < dim ; i++ ) {
				test += vec[i] ;
			}
			if ( fabs(test) > 1.0e-06 ) {
				cout << "Vector was not normalized" << endl ;
				exit(1) ;
			}
		}
	double get(int idx) const 
		{
#ifdef DEBUG			
			if ( idx >= dim ) {
				cout <<  "Vector out of bounds" << endl ;
				exit(1) ;
			}
#endif			
			return vec[idx] ;
		}
	void add(int idx, double val) 
		{
			vec[idx] += val ;
		}
	void print() 
		{
			if ( RANK == 0 ) {
				cout << "[" << endl ;
				for ( int j = 0 ; j < dim ; j++ ) {
					cout << j << " " << vec[j] << endl ;
				}
				cout << "]" << endl ;
			}
		}
	void print(ostream &of) 
		{
			if ( RANK == 0 ) {
				for ( int j = 0 ; j < dim ; j++ ) {
					of << j << " " << vec[j] << endl ;
				}
			}
		}
	void scale(Vector &out, double val) 
	// Scale the vector by the given value, put result in Out.
		{
			for ( int j = 0 ; j < dim ; j++ ) {
				out.set(j, val * vec[j] ) ;
			}
		}
	double l1norm()
	// Returns L1 norm (sum of abs values).
		{
			double norm = 0 ;
			for ( int i = 0 ; i < dim ; i++ ) {
				norm += fabs(vec[i]) ;
			}
			return norm ;
		}
	void remove(int idx)
	// Remove the specified index from the vector.
		{
			if ( idx < 0 || idx >= dim ) {
				cout << "Error: bad index to remove from vector: " << idx << endl ;
			}
			for ( int i = idx ; i < dim - 1 ; i++ ) {
				vec[i] = vec[i+1] ;
			}
			--dim ;
		}
	
} ;

class IntVector {
public:
	int *vec ;
	int dim ;

	IntVector(int d1)
		{
			dim = d1 ;
			
				vec = new int[d1] ;
		}
	IntVector(int d1, int val)
		{
			dim = d1 ;
			vec = new int[d1] ;
			for ( int j = 0 ; j < dim ; j++ ) {
				vec[j] = val ;
			}
		}
	IntVector()
		{
			dim = 0 ;
			vec = NULL ;
		}
	~IntVector()
		{
			
				delete[] vec ;
		}
	void realloc(int size) {
		// Reallocate the vector.
		delete [] vec ;
		vec = new int[size] ;
		dim = size ;
	}

	int size() const {
		return dim ;
	}
	void set(int i, int val) {
#ifdef DEBUG					
		if ( i >= dim ) {
			cout << "IntVector set out of bounds" << endl ;
			exit(1) ;
		}
#endif					
		vec[i] = val ;
	}
	void read(ifstream &file, int dim0) 
		{
			dim = dim0 ;
			vec = new int[dim] ;
			for ( int i = 0 ; i < dim ; i++ ) {
				int val ;
				file >> val ;
				set(i, val) ;
			}
		}
	int get(int idx) const {
#ifdef DEBUG					
		if ( idx >= dim ) {
			cout << "IntVector index out of bounds" << endl ;
			exit(1) ;
		}
#endif					
		return vec[idx] ;
	}
	void add(int idx, int val) {
		vec[idx] += val ;
	}

	void print() 
		{
			if ( RANK == 0 ) {
				cout << "[" << endl ;
				for ( int j = 0 ; j < dim ; j++ ) {
					cout << j << " " << vec[j] << endl ;
				}
				cout << "]" << endl ;
			}
		}

		void print(ofstream &of) 
		{
			if ( RANK == 0 ) {
				for ( int j = 0 ; j < dim ; j++ ) {
					of << j << " " << vec[j] << endl ;
				}
			}
		}
	
	void remove(int idx)
	// Remove the specified index from the vector.
		{
			if ( idx < 0 || idx >= dim ) {
				cout << "Error: bad index to remove from vector: " << idx << endl ;
			}
			for ( int i = idx ; i < dim - 1 ; i++ ) {
				vec[i] = vec[i+1] ;
			}
			--dim ;
		}

} ;

class Matrix {
public:
	double *mat ;  // Elements of the matrix stored here.
	double *shift ; // Normalization shift of each column ;
	double *scale ; // Normalization scale factor for each column.
	int dim1, dim2 ; 			// dim1 is number of rows.  Dim2 is number of columns.
	bool distributed ;    // Is this matrix distributed over processes ?
	int row_start, row_end ;  // Starting and ending row for parallel calculations.
	int num_rows ;    // Number of rows stored on this process.
	
	Matrix(const Matrix &matin) {
		// Create a matrix that is a copy of another.
		dim1 = matin.dim1 ;
		dim2 = matin.dim2 ;
		distributed = matin.distributed ;
		row_start = matin.row_start ;
		row_end = matin.row_end ;
		num_rows = matin.num_rows ;
		
		mat = new double[num_rows * dim2] ;
		shift = new double[dim2] ;
		scale = new double[dim2] ;
		
		for ( int j = row_start ; j <= row_end ; j++ ) {
			for ( int k = 0 ; k < dim1 ; k++ ) {
				set(j,k, matin.get(j,k) ) ;
			}
		}
		for ( int j = 0 ; j < dim2 ; j++ ) {
			shift[j] = matin.shift[j] ;
			scale[j] = matin.scale[j] ;
		}
	}
	
	Matrix(int d1, int d2) 
		{
			dim1 = d1 ;
			dim2 = d2 ;
			mat = new double[d1 * d2] ;
			shift = new double[d2] ;
			scale = new double[d2] ;
			row_start = 0 ;
			row_end = dim1 - 1 ;
			num_rows = dim1 ;
			distributed = false ;
		}
		
	Matrix(int d1, int d2, bool will_dist)
	// Specify whether to distribute or not by will_dist argument
		{
			dim1 = d1 ;
			dim2 = d2 ;
			shift = new double[d2] ;
			scale = new double[d2] ;

			if ( ! will_dist ) {
				mat = new double[d1 * d2] ;
				row_start = 0 ;
				row_end = dim1 - 1 ;
				num_rows = dim1 ;
				distributed = false ;
			} else {
				distribute() ;
				mat = new double[num_rows * d2] ;
			}
		}
	
	Matrix()	
		{
			dim1 = 0 ;
			dim2 = 0 ;
			mat = NULL ;
			shift = NULL ;
			scale = NULL ;
			row_start = 0 ;
			row_end = -1 ;
			num_rows = 0 ;
			distributed = false ;
		}
	~Matrix() {
		delete [] mat ;
		delete [] shift ;
		delete [] scale ;
		dim1 = 0 ;
		dim2 = 0 ;
	}
	void read(std::ifstream &file, int dim01, int dim02)
	// Read a non-distributed matrix.
		{
			dim1 = dim01 ;
			dim2 = dim02 ;
			if ( mat != NULL ) {
				delete [] mat ;
			}
			mat = new double [dim1 * dim2] ;
			shift = new double[dim2] ;
			scale = new double[dim2] ;
			for ( int i = 0 ; i < dim1 ; i++ ) {
				for ( int j = 0 ; j < dim2 ; j++ ) {
					double val ;
					file >> val ;
					set(i, j, val) ;
				}
			}
		}
	void realloc(int d1, int d2) 
		{
			if ( dim1 > 0 && dim2 > 0 ) {
				delete [] mat ;
				delete [] scale ;
				delete [] shift ;
			}
			mat = new double[d1 * d2] ;
			scale = new double[d2] ;
			shift = new double[d2] ;
			dim1 = d1 ;
			dim2 = d2 ;
			if ( ! distributed ) {
				row_start = 0 ;
				row_end = dim1 - 1 ;
			} else {
				row_start = (dim1 * RANK) / NPROCS ;
				if ( RANK < NPROCS - 1 ) {
					row_end = (dim1 * (RANK+1)) / NPROCS - 1 ;
				} else {
					row_end = dim1 - 1 ;
				}
			}
			num_rows = row_end - row_start + 1 ;
		}
	inline double get(int i, int j) const
		{
#ifdef DEBUG						
			if ( i >= dim1 || j >= dim2 ) {
				cout << "Matrix out of bounds" << endl ;
				exit(1) ;
			}
			if ( distributed && (i > row_end || i < row_start) ) {
				cout << "Distributed matrix out of bounds " << endl ;
				exit(1) ;
			}
#endif						
			return(mat[(i-row_start) * dim2 + j]) ;
		}
	inline void set(int i, int j, double val) 
		{
#ifdef DEBUG						
			if ( i >= dim1 || j >= dim2 ) {
				cout << "Matrix set out of bounds" << endl ;
				exit(1) ;
			}
			if ( distributed && (i > row_end || i < row_start) ) {
				cout << "Distributed matrix set out of bounds " << endl ;
				exit(1) ;
			}
#endif						
			mat[(i-row_start) * dim2 + j] = val  ;
		}
	inline void setT(int i, int j, double val) {
		// Set a value with transposed indexing.
		set(j, i, val) ;
	}
	inline double getT(int i, int j) {
		// Get a value with transposed indexing.
		return get(j, i) ;
	}
	void distribute()
	// Settings to distribute a matrix across processes.
		{
						
			row_start = (dim1 * RANK) / NPROCS ;
			if ( RANK < NPROCS - 1 ) {
				row_end = (dim1 * (RANK+1)) / NPROCS - 1 ;
			} else {
				row_end = dim1 - 1 ;
			}
			distributed = true ;
			num_rows = row_end - row_start + 1 ;
		}
	double mult_T(int j, int k)
	// Calculate the j,k element of X^T * X, where X is the current matrix.
		{
			double tmp = 0.0 ;
			for ( int l = row_start ; l <= row_end ; l++ ) {
				tmp += get(l, j) * get(l, k) ;
			}
			double tmp2 ;
#ifdef USE_MPI
			MPI_Allreduce(&tmp, &tmp2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else
			tmp2 = tmp ;
#endif
			return tmp2 ;
		}
	void normalize()
	// Normalize matrix to according to Eq. 1.1
		{
			for ( int j = 0 ; j < dim2 ; j++ ) {
				shift[j] = 0.0 ;
				for ( int k = 0 ; k < dim1 ; k++ ) {
					shift[j] += get(k, j) ;
				}
				shift[j] /= dim1 ;
				for ( int k = 0 ; k < dim1 ; k++ ) {
					double val = get(k, j) - shift[j] ;
					set(k, j, val) ;
				}
			}
			for ( int j = 0 ; j < dim2 ; j++ ) {
				scale[j] = 0.0 ;
				for ( int k = 0 ; k < dim1 ; k++ ) {
					scale[j] += get(k, j) * get(k,j) ;
				}
				scale[j] = sqrt(scale[j]) ;
				for ( int k = 0 ; k < dim1 ; k++ ) {
					double val = get(k, j) / scale[j] ;
					set(k, j, val) ;
				}
			}
		}
	void check_norm() const
	// Normalize matrix to according to Eq. 1.1
		{
			for ( int j = 0 ; j < dim2 ; j++ ) {
				double test = 0.0 ;
				for ( int k = 0 ; k < dim1 ; k++ ) {
					test += get(k, j) ;
				}
				if ( fabs(test) > 1.0e-06 ) {
					cout << "Column " << j << " did not have 0 average" << endl ;
				}
			}
			for ( int j = 0 ; j < dim2 ; j++ ) {
				double test = 0.0 ;
				for ( int k = 0 ; k < dim1 ; k++ ) {
					test += get(k, j) * get(k,j) ;
				}
				if ( fabs(test-1.0) > 1.0e-06 ) {
					cout << "Column " << j << " did not have norm = 1 " << endl ;
				}
			}
		}
	bool cholesky(Matrix &chol)
	// Calculate the cholesky decomposition of the current matrix, and store in chol.
	// Uses the upper triangular variant, A = R^T * R.  R is calculated.
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

	bool cholesky_add_row(const Matrix &chol0, const Vector &newr)
	// Given cholesky decomposition of a matrix in chol0
	// find the cholesky decomposition of a new matrix formed by adding
	// the specified row (newr) to the original matrix (not given).
	// The new row is placed in the last index.
		{
			if ( dim1 != dim2 ) {
				cout << "Error: Cholesky decomposition only works for square matrices " << endl ;
			}
			if ( newr.dim != dim1 || chol0.dim1 != dim1 - 1 || chol0.dim1 != chol0.dim2 ) {
				cout << "Dimension mismatch" ;
				exit(1) ;
			}
			Vector xtmp(dim1) ;
			double cdotc = 0.0 ;
			for ( int j = 0 ; j < dim1 - 1 ; j++ ) {
				double sum = newr.get(j) ;
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

			for ( int j = 0 ; j < dim1 - 1 ; j++ ) {
				for ( int k = 0 ; k <= j ; k++ ) {
					set(k,j, chol0.get(k,j)) ;
				}
				for ( int k = j + 1 ; k < dim1 ; k++ ) {
					set(k,j, 0.0) ;
				}
			}
			for ( int j = 0 ; j < dim1 ; j++ ) {
				set(j, dim1-1, xtmp.get(j)) ;
			}
			for ( int j = 0 ; j < dim1 - 1 ; j++ ) {
				set(dim1-1, j, 0.0) ;
			}
			return true ;
		}
	
	bool cholesky_remove_row(int id ) {
		// Based on larscpp github package.  Modified so that diagonal elements of the cholesky
		// matrix are always positive.
		// Remove a row from a Cholesky matrix.  Update the dimension of the matrix
		// Use transpose matrix access functions because original code was written in terms of a lower
		// triangular matrix, whereas we use an upper triangular matrix.
		double a(0),b(0),c(0),s(0),tau(0);
		if ( dim1 != dim2 ) {
			cout << "Error: the cholesky matrix should be square" << endl ;
			exit(1) ;
		}
		int lth = dim1 - 1 ;
		for(int i=id; i < lth; ++i){
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
			for( int j=i+2; j<=lth; ++j){
				a = getT(j,i);
				b = getT(j, i+1);
				setT(j, i, c*a - s*b) ;
				setT(j, i+1, s*a + c*b) ;
			}
		}
		for( int i=0; i<= lth; ++i)
			setT(lth, i, 0) ;

		// Shift the storage to accommodate a change of dimension.
		for ( int i = 0 ; i < lth ; i++ ) {
			for ( int j = 0 ; j < lth ; j++ ) {
				mat[i * (dim2-1) + j] = mat[i * dim2 + j] ;
			}
		}
		dim1-- ;
		dim2-- ;
		num_rows-- ;
		row_end-- ;
		
		return true ;
	}
	
	void cholesky_sub(Vector &x, const Vector &b)
	// Assuming that the current matrix holds a Cholesky decomposition
	// in R^T * R form, calculate the solution vector x given the vector
	// of input values b.
		{
			if ( b.dim != x.dim || b.dim != dim1 || b.dim != dim2 ) {
				cout << "Dimension mismatch" ;
				exit(1) ;
			}
			Vector xtmp(dim1) ;
			for ( int j = 0 ; j < dim1 ; j++ ) {
				double sum = b.get(j) ;
				for ( int k = 0 ; k < j ; k++ ) {
					sum -= get(k,j) * xtmp.get(k) ;
				}
				xtmp.set(j, sum / get(j,j) );
			}
			for ( int j = dim1 - 1 ; j >= 0 ; j-- ) {
				double sum = xtmp.get(j) ;
				for ( int k = j + 1 ; k < dim1 ; k++ ) {
					sum -= get(j,k) * x.get(k) ;
				}
				x.set(j, sum / get(j,j) ) ;
			}
		}
			
	void print()
		{
			if ( RANK == 0 ) {
				cout << "[" << endl ;
				for ( int j = 0 ; j < dim1 ; j++ ) {
					cout << "[" << endl ;
					for ( int k = 0 ; k < dim2 ; k++ ) {
						cout << j << " " << k << " " << get(j,k) << endl ;
					}
					cout << "]" << endl ;
				}
				cout << "]" << endl ;
			}
		}

		void print(ofstream &of)
		{
			if ( RANK == 0 ) {
				for ( int j = 0 ; j < dim1 ; j++ ) {
					for ( int k = 0 ; k < dim2 ; k++ ) {
						of << get(j,k) << " " ;
					}
					of << endl ;
				}
			}
		}
	
	void dot(Vector &out, const Vector &in)
	// Find matrix * in = out
		{
			if ( out.dim != dim1 || in.dim != dim2 ) {
				cout << "Array dimension mismatch" << endl ;
				exit(1) ;
			}
			for ( int j = row_start ; j <= row_end ; j++ ) {
				double sum = 0.0 ;
				for ( int k = 0 ; k < dim2 ; k++ ) {
					sum += get(j,k) * in.get(k) ;
				}
				out.set(j, sum) ;
			}
#ifdef USE_MPI
			if ( distributed ) {
				Vector out2(dim1, 0.0) ;    // Temporary array to collect out.vec values.
				IntVector countv(NPROCS) ;  // The number of items to receive from each process.
				IntVector displs(NPROCS) ;  // Storage displacements in out2 for each process.

				int count = row_end - row_start + 1 ;

				MPI_Allgather(&count, 1, MPI_INT, countv.vec, 1, MPI_INT, MPI_COMM_WORLD) ;

				displs.set(0,0) ;
				for ( int j = 1 ; j < NPROCS ; j++ ) {
					displs.set(j, displs.get(j-1) + countv.get(j-1) ) ;
				}
				MPI_Allgatherv(&(out.vec[row_start]), countv.get(RANK),
											 MPI_DOUBLE, out2.vec, countv.vec, displs.vec, MPI_DOUBLE, MPI_COMM_WORLD) ;
				
				for ( int j = 0 ; j < dim1 ; j++ ) {
					out.set(j, out2.get(j) ) ;
				}
			}
#endif				
		}

	void dot_transpose(Vector &out, const Vector &in)
	// Find Transpose(matrix) * in = out
		{
			if ( out.dim != dim2 || in.dim != dim1 ) {
				cout << "Array dimension mismatch" << endl ;
				exit(1) ;
			}
			Vector sumv(dim2,0.0) ;			
			for ( int j = 0 ; j < dim2 ; j++ ) {
				for ( int k = row_start ; k <= row_end ; k++ ) {
					sumv.add(j, get(k,j) * in.get(k) ) ;
				}
			}
#ifdef USE_MPI
			if ( distributed ) {
				MPI_Allreduce(sumv.vec, out.vec, dim2,
											MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
			} else {
				for ( int j = 0 ; j < dim2 ; j++ ) {
					out.set(j, sumv.get(j)) ;
				}
			}
#else				
			for ( int j = 0 ; j < dim2 ; j++ ) {
				out.set(j, sumv.get(j)) ;
			}
#endif				
		}
} ;

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
	Vector a ;         
	double gamma ;      // LARS Step size ;
	double gamma_lasso ;  // Lasso constraint on LARS step size.
	double gamma_use ;  // The value of gamma to use on current step.  Based on gamma and gamma_lasso.
	int ndata ;       // Number of data items to fit to = X.dim1
	int nprops ;      // Number of properties to correlate = X.dim2
	int nactive ;     // The number of active properties to correlate <= nprops
	int remove_prop ;  // The index of the property to remove from the active set during a LASSO calculation.
	bool do_lasso ;   // If TRUE, do a lasso calculation.  If false, do a regular LARS calculation.
	bool solve_succeeded ; // If true, the solve of G_A succeeded.
	
	DLARS(Matrix &Xin, Vector &yin): X(Xin.dim1, Xin.dim2,true), y(Xin.dim1), mu(Xin.dim1),
																	 beta(Xin.dim2), c(Xin.dim2), A(0),  exclude(Xin.dim2, 0), sign(0) 
		{
			do_lasso = false ;
			gamma_lasso = 1.0e20 ;
			gamma = 0.0 ;
			gamma_use = 0.0 ;
			ndata = X.dim1 ;
			nprops = X.dim2 ;
			nactive = 0 ;
			remove_prop = -1 ;
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
	// Calculated predicted values of y (mu hat, Eq. 1.2)
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
	double sq_error() 
	// Squared error Eq. 1.3
		{
			double result = 0.0 ;
			for ( int j = 0 ; j < ndata ; j++ ) {
				result += (y.get(j)-mu.get(j)) * (y.get(j) - mu.get(j) ) ;
			}
			return(result) ;
		}
	void correlation() 
	// Calculate the correlation vector c, Eq. 2.1
		{
			C_max = -1.0 ;
			Vector ydiff(ndata,0.0) ;
			for ( int k = 0 ; k < ndata ; k++ ) {
				ydiff.set(k, y.get(k) - mu.get(k)) ;
			}

			X.dot_transpose(c, ydiff) ;
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

			// DEBUG !!
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
			const double eps = 1.0e-20 ;

			// DEBUG !!
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
			for ( int j = 0 ; j < nactive ; j++ ) {
				if ( fabs(test.get(j) - 1.0) > 1.0e-04 ) {
					cout << "Cholesky solution test failed\n" ;
					cout << "Error = " << fabs(test.get(j) - 1.0) << endl ;
					return false ;
				}
			}
			
			A_A = 0.0 ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				A_A += G_A_Inv_I.get(j) ;
			}
			if ( A_A > eps ) 
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
			} else {
				// Maintain the ordering of the last active set to allow cholesky updating.
				for ( int j = 0 ; j < nactive ; j++ ) {
					a_trial.set( j, A_last.get(j) ) ;
				}
				int count = nactive ;
				// Search for the new active set.
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
			int add_prop = -1 ;
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
					double cmin ;

					if ( c1 < 0.0 ) {
						cmin = ( c2 > 0.0 ) ? c2 : huge ;
					}
					else if ( c2 < 0.0 ) {
						cmin = ( c1 > 0.0 ) ? c1 : huge ;
					}
					else {
						cmin = ( c1 < c2 ) ? c1 : c2 ;
					}
					// DEBUG !!
					//cout << "Property " << j << " gamma = " << cmin << endl ;
					
					if ( cmin < gamma ) {
						gamma = cmin ;
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
			const double eps = 1.0e-10 ;

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

				// DEBUG !! 
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
				out << "Y constant offset = " << offset << endl ;
				out << "Unscaled coefficients: " << endl ;
				for ( int j = 0 ; j < nprops ; j++ ) {
					uns_beta.set(j, beta.get(j) / X.scale[j]) ;
				}
				if ( out == cout ) {
					uns_beta.print() ;
				} else {
					uns_beta.print(out) ;
				}
			}
		}

	
	void print_unshifted_mu(ostream &out)
	// Print the prediction in unscaled units.
		{
			if ( RANK == 0 ) {
				double offset = y.shift ;

				out << "Y constant offset = " << offset << endl ;
				for ( int j = 0 ; j < ndata ; j++ ) {
					out << mu.get(j) + y.shift << endl ;
				}
			}
		}

	bool iteration()
	// Perform a single iteration of the LARS algorithm.
	// Return false when no more iterations can be performed.
		{

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
			if ( RANK == 0 ) cout  << "L1 norm of solution: " << beta.l1norm() << " RMS Error: " << sqrt(sq_error() / ndata) << endl ;

			update_active_set() ;

			// build the X_A array.
			build_X_A() ;
			build_G_A() ;

			if ( ! solve_G_A() ) {
				remove_prop = -1 ;
				cout << "Iteration failed" << endl ;
				return true ;
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
				print_unscaled(cout) ;
			}

			if ( nactive < nprops - num_exclude ) {
				return true ;
			} else {
				return false ;
			}
		}
};
	
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

	if ( argc < 5 ) {
		if ( RANK == 0 ) cout << "Not enough args " << endl ;
		exit(1) ;
	}
	string xname(argv[1]) ;
	string yname(argv[2]) ;
	string dname(argv[3]) ;
	string algorithm(argv[4]) ;

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
	
	ifstream dfile(dname) ;
	if ( ! dfile.is_open() ) {
		if ( RANK == 0 ) cout << "Error: could not open " << dname << endl ;
		exit(1) ;
	}
	dfile >> nprops >> ndata ;
	
	ifstream xfile(xname) ;
	if ( ! xfile.is_open() ) {
		if ( RANK == 0 ) cout << "Could not open " << xname << endl ;
		exit(1) ;
	}
	Matrix xmat ;
	xmat.read(xfile, ndata, nprops) ;
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

#if(0)
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
#endif

