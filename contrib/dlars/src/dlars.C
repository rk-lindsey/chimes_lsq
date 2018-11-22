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
		if ( dim > 0 ) 
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
			cout << "[" << endl ;
			for ( int j = 0 ; j < dim ; j++ ) {
				cout << j << " " << vec[j] << endl ;
			}
			cout << "]" << endl ;
		}
	void scale(Vector &out, double val) 
	// Scale the vector by the given value, put result in Out.
		{
			for ( int j = 0 ; j < dim ; j++ ) {
				out.set(j, val * vec[j] ) ;
			}
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
		if ( dim > 0 ) 
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
			cout << "[" << endl ;
			for ( int j = 0 ; j < dim ; j++ ) {
				cout << j << " " << vec[j] << endl ;
			}
			cout << "]" << endl ;
		}

} ;

class Matrix {
public:
	double *mat ;  // Elements of the matrix stored here.
	double *shift ; // Normalization shift of each column ;
	double *scale ; // Normalization scale factor for each column.
	int dim1, dim2 ; 			// dim1 is number of rows.  Dim2 is number of columns.

	Matrix(int d1, int d2) 
		{
			dim1 = d1 ;
			dim2 = d2 ;
			mat = new double[d1 * d2] ;
			shift = new double[d2] ;
			scale = new double[d2] ;
		}
	Matrix()	
		{
			dim1 = 0 ;
			dim2 = 0 ;
			mat = NULL ;
			shift = NULL ;
			scale = NULL ;
		}
	~Matrix() {
		delete [] mat ;
		delete [] shift ;
		delete [] scale ;
		dim1 = 0 ;
		dim2 = 0 ;
		
	}
	void read(std::ifstream &file, int dim01, int dim02)
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
			}
			mat = new double[d1 * d2] ;
			dim1 = d1 ;
			dim2 = d2 ;
		}
	inline double get(int i, int j) 
		{
#ifdef DEBUG						
			if ( i >= dim1 || j >= dim2 ) {
				cout << "Matrix out of bounds" << endl ;
				exit(1) ;
			}
#endif						
			return(mat[i * dim2 + j]) ;
		}
	inline void set(int i, int j, double val) 
		{
#ifdef DEBUG						
			if ( i >= dim1 || j >= dim2 ) {
				cout << "Matrix set out of bounds" << endl ;
				exit(1) ;
			}
#endif						
			mat[i * dim2 + j] = val  ;
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
	void check_norm()
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
	void dot(Vector &out, const Vector &in)
	// Find matrix * in = out
		{
			if ( out.dim != dim1 || in.dim != dim2 ) {
				cout << "Array dimension mismatch" << endl ;
				exit(1) ;
			}
			for ( int j = 0 ; j < dim1 ; j++ ) {
				double sum = 0.0 ;
				for ( int k = 0 ; k < dim2 ; k++ ) {
					sum += get(j,k) * in.get(k) ;
				}
				out.set(j, sum) ;
			}
		}
	void dot_transpose(Vector &out, const Vector &in)
	// Find Transpose(matrix) * in = out
		{
			if ( out.dim != dim2 || in.dim != dim1 ) {
				cout << "Array dimension mismatch" << endl ;
				exit(1) ;
			}
			for ( int j = 0 ; j < dim2 ; j++ ) {
				double sum = 0.0 ;
				for ( int k = 0 ; k < dim1 ; k++ ) {
					sum += get(k,j) * in.get(k) ;
				}
				out.set(j, sum) ;
			}
		}
} ;

class DLARS
{
public:
	Matrix X ;         // The matrix of data vs properties
	Vector y ;         // The vector of data
	Vector mu ;        // The predicted values of data
	Vector beta ;      // The scaled fitting coefficients
	Vector c ;         // The correlation vector.
	IntVector A ;      // Indices of the active set of properties.
	IntVector A_last ;  // Indices of the last active set of properties.
	IntVector exclude ; // Indices of properties not to use because of degeneracy.
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
	double gamma ;  // Step size ;
	int ndata ;       // Number of data items to fit to = X.dim1
	int nprops ;      // Number of properties to correlate = X.dim2
	int nactive ;     // The number of active properties to correlate <= nprops
	
	DLARS(Matrix &Xin, Vector &yin): X(Xin.dim1, Xin.dim2), y(Xin.dim1), mu(Xin.dim1),
																	 beta(Xin.dim2), c(Xin.dim2), A(0),  exclude(Xin.dim2, 0), sign(0) 
		{
			ndata = X.dim1 ;
			nprops = X.dim2 ;
			nactive = 0 ;
			
			for ( int j = 0 ; j < ndata ; j++ ) {
				for ( int k = 0 ; k < nprops; k++ ) {
					X.set(j, k, Xin.get(j,k) ) ;
				}
				y.set(j, yin.get(j)) ;
				mu.set(j, 0.0) ;
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
			if ( mu.size() != ndata ||
				  beta.size() != nprops )  {
				cout << "Error:  matrix dim mismatch" << endl ;
				exit(1) ;
			}

			for ( int j = 0 ; j < ndata ; j++ ) {
				double tmp = 0.0 ;
				for ( int k = 0 ; k < nprops ; k++ ) {
					tmp += X.get(j,k) * beta.get(k) ;
				}
				mu.set(j, tmp) ;
			}
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
			for ( int j = 0 ; j < nprops ; j++ ) {
				c.set(j, 0.0) ;
				for ( int k = 0 ; k < ndata ; k++ ) {
					double val = X.get(k,j) * ( y.get(k) - mu.get(k) ) ;
					c.add(j, val) ;
				}
				if ( fabs(c.get(j)) > C_max ) {
					C_max = fabs(c.get(j)) ;
				}
			}
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
		for ( int j = 0 ; j < ndata ; j++ ) {
			for ( int k = 0 ; k < nactive ; k++ ) {
				double val = X.get( j, A.get(k) ) * sign.get(k) ;
				X_A.set(j,k, val) ;
			}
		}
	}

	void build_G_A()
	// Build the G_A matrix after X_A has been built.
		{
			G_A.realloc(nactive, nactive) ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				for ( int k = 0 ; k < nactive ; k++ ) {
					double tmp = 0.0 ;
					for ( int l = 0 ; l < ndata ; l++ ) {
						tmp += X_A.get(l, j) * X_A.get(l, k) ;
					}
					G_A.set(j, k, tmp) ;
				}
			}
		}
	bool solve_G_A()
	// Find G_A^-1 * I
		{
			const double eps = 1.0e-20 ;
			chol.realloc(nactive, nactive) ;
			G_A_Inv_I.realloc(nactive) ;

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
					}
				}
				return false ;
			}
			// Solve for G_A^-1 * unity
			Vector unity(nactive, 1.0) ;
			chol.cholesky_sub(G_A_Inv_I, unity) ;

			// Test to see if the solution worked.
			Vector test(nactive) ;
			G_A.dot(test, G_A_Inv_I) ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				if ( fabs(test.get(j) - 1.0) > 1.0e-06 ) {
					cout << "Cholesky solution test failed\n" ;
					exit(1) ;
				}
			}
			
			A_A = 0.0 ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				A_A += G_A_Inv_I.get(j) ;
			}
			if ( A_A > eps ) 
				A_A = 1.0 / sqrt(A_A) ;
			else {
				cout << "Normalization failed" << endl ;
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
	void update_active_set() 
	// Update the active set of directions to those having maximum correlation.
		{
			IntVector a_trial(nprops) ;
			int count = 0 ;
			const double eps = 1.0e-6 ;

			// Save the last active set
			A_last.realloc(nactive) ;
			for ( int j = 0 ; j < nactive ; j++ ) {
				A_last.set( j, A.get(j) ) ;
			}

			// Search for the new active set.
			// Some indices could be excluded if they are degenerate.
			for ( int j = 0 ; j < nprops ; j++ ) {
				if ( fabs( fabs(c.get(j)) - C_max ) < eps
						 && ! exclude.get(j) ) {
					a_trial.set(count, j) ;
					count++ ;
				}
			}
			A.realloc(count) ;
			nactive = count ;
			
			for ( int j = 0 ; j < nactive ; j++ ) {
				A.set(j, a_trial.get(j)) ;
			}
#ifdef VERBOSE			
			cout << "New active set: " << endl ;
			A.print() ;
#endif			
		}

	void update_step_gamma()
	// Update gamma, eq. 2.13.
		{
			double huge = 1.0e20 ;
			gamma = huge ;

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
					if ( cmin < gamma ) {
						gamma = cmin ;
					}
				}
			} else {
				// Active set = all variables.
				gamma = C_max / A_A ;
			}
			cout << "Updated step gamma = " << gamma << endl ;
		}
	void update_beta()
		{
			for ( int j = 0 ; j < nactive ; j++ ) {
				int idx = A.get(j) ;
				double val = beta.get(idx) + w_A.get(j) * sign.get(j) * gamma ;
				beta.set(idx,val) ;
			}
#ifdef VERBOSE			
			cout << "New beta: " ;
			beta.print() ;

			cout << "Predicted mu: " << endl ;
			for ( int j = 0 ; j < ndata ; j++ ) {
				cout << mu.get(j) + gamma * u_A.get(j) << " " ;
			}
			cout << "\n" ;
#endif			
		}

	void print_unscaled()
	// Print the coefficients in unscaled units.
		{
			double offset = y.shift ;
			Vector uns_beta(nprops,0.0) ;
			for ( int j = 0 ; j < nprops ; j++ ) {
				if ( X.scale[j] == 0.0 ) {
					cout << "Error: scale factor = 0.0" << endl ;
					exit(1) ;
				}
				offset -= beta.get(j) * X.shift[j] / X.scale[j] ;
			}
			cout << "Y constant offset = " << offset << endl ;
			cout << "Unscaled coefficients: " << endl ;
			for ( int j = 0 ; j < nprops ; j++ ) {
				uns_beta.set(j, beta.get(j) / X.scale[j]) ;
			}
			uns_beta.print() ;
		}
			
	void iteration()
	// Perform a single iteration of the LARS algorithm.
		{

			predict() ;
			correlation() ;

#ifdef VERBOSE
			cout << "Pre-step beta: " << endl ;
			beta.print() ;
			cout << "Prediction: " << endl ;
			mu.print() ;
			cout << " Correlation: " << endl ;
			c.print() ;
			cout << "C_max: " << C_max << endl ;
#endif			
			cout << "RMS Error " << sqrt(sq_error() / ndata) << endl ;

			update_active_set() ;

			// build the X_A array.
			build_X_A() ;
			build_G_A() ;
			if ( ! solve_G_A() ) {
				cout << "Iteration failed" << endl ;
				return ;
			}
			
#ifdef VERBOSE			
			cout << "X_A" << endl ;
			X_A.print() ;
			cout << "G_A" << endl ;
			G_A.print() ;
			cout << "G_A_Inv " << endl ;
			G_A_Inv_I.print() ;
			cout << "A_A " << A_A << endl ;
#endif			

			build_u_A() ;
			update_step_gamma() ;
			update_beta() ;

			cout << "Beta: " << endl ;
			beta.print() ;
			print_unscaled() ;
		}
};
	
int main(int argc, char **argv)
{
	cout << "Distributed LARS algorithm" << endl ;

	if ( argc < 4 ) {
		cout << "Not enough args " << endl ;
		exit(1) ;
	}
	string xname(argv[1]) ;
	string yname(argv[2]) ;
	string dname(argv[3]) ;

	int nprops, ndata ;
	ifstream dfile(dname) ;
	dfile >> nprops >> ndata ;
	
	ifstream xfile(xname) ;
	if ( ! xfile.is_open() ) {
		cout << "Could not open " << xname << endl ;
		exit(1) ;
	}
	Matrix xmat ;
	xmat.read(xfile, ndata, nprops) ;
	xmat.normalize() ;
	xmat.check_norm() ;

	ifstream yfile(yname) ;
	if ( ! yfile.is_open() ) {
		cout << "Could not open " << yname << endl ;
		exit(1) ;
	}
	Vector yvec ;
	yvec.read(yfile, ndata) ;
	yvec.normalize() ;
	yvec.check_norm() ;

#ifdef VERBOSE	
	cout << "Normalized X matrix" << endl ;
	xmat.print() ;

	cout << "Normalized Y vector" << endl ;
	for ( int j = 0 ; j < yvec.dim ; j++ ) {
		cout << yvec.get(j) << " " ;
	}
	cout << endl ;
#endif	

	DLARS lars(xmat, yvec) ;

	for ( int j = 0 ; j < xmat.dim2 ; j++ ) {
		cout << "Working on iteration " << j + 1 << endl ;
		lars.iteration() ;
	}

	cout << "Final values:" << endl ;
	cout << "Beta: " << endl ;
	lars.beta.print() ;

	lars.predict() ;
	lars.correlation() ;

	cout << "Prediction: " << endl ;
	lars.mu.print() ;

	cout << "Sq Error " << lars.sq_error() << endl ;

#ifdef VERBOSE
	cout << " Correlation: " << endl ;
	lars.c.print() ;
#endif	

}

