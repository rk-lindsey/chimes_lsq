
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
			for ( int k = 0 ; k < dim2 ; k++ ) {
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
			for ( int j = 0 ; j < dim2 ; j++ ) {
				shift[j] = 0.0 ;
				scale[j] = 1.0 ;
			}
		}
		
	Matrix(int d1, int d2, bool will_dist)
	// Specify whether to distribute or not by will_dist argument
		{
			dim1 = d1 ;
			dim2 = d2 ;
			shift = new double[d2] ;
			scale = new double[d2] ;
			for ( int j = 0 ; j < dim2 ; j++ ) {
				shift[j] = 0.0 ;
				scale[j] = 1.0 ;
			}
			
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
	void read(std::ifstream &file, int dim01, int dim02, bool is_distributed)
	// Read a non-split matrix from a single file.  The matrix is optionally distributed among processes.
		{
			dim1 = dim01 ;
			dim2 = dim02 ;
			if ( mat != NULL ) {
				delete [] mat ;
			}
			shift = new double[dim2] ;
			scale = new double[dim2] ;

			if ( ! is_distributed ) {
				distributed = false ;
				mat = new double [dim1 * dim2] ;
				for ( int i = 0 ; i < dim1 ; i++ ) {
					for ( int j = 0 ; j < dim2 ; j++ ) {
						double val ;
						file >> val ;
						set(i, j, val) ;
					}
				}
			} else {
				distribute() ;
				mat = new double [num_rows * dim2] ;
				for ( int i = 0 ; i < dim1 ; i++ ) {
					for ( int j = 0 ; j < dim2 ; j++ ) {
						double val ;
						file >> val ;
						if ( i >= row_start && i <= row_end ) {
							set(i, j, val) ;
						}
					}
				}
			}
			for ( int j = 0 ; j < dim2 ; j++ ) {
				shift[j] = 0.0 ;
				scale[j] = 1.0 ;
			}
		}
	
	void read_split_files(const char* matFilename, const char* dimFilename)
	// Read split file output from chimes_lsq.
	{
		ifstream dim_file ;
		char name[80] ;

		distributed = true ;
		
		string strDimFile(dimFilename) ;
		std::size_t found1 = strDimFile.find(".") ;
		if ( found1 == string::npos ) 
		{
			cerr << "A dimension file name must end with a suffix" ;
			stop_run(1) ;
		}
		string dim_ext  = strDimFile.substr(found1+1, string::npos) ;
		strDimFile = strDimFile.substr(0,found1) ;
		// Find out how many files were written by chimes_lsq.
		int total_files = 0 ;
		ifstream test_file ;
		for ( int j = 0 ; j < NPROCS + 1 ; j++ ) 
		{
			memset(name, 0, 80) ;
			sprintf(name, "%s.%04d.%s", strDimFile.c_str(), j, dim_ext.c_str()) ;
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
			if ( RANK == 0 ) {
				cout << "Not enough processes specified" << endl ;
				cout << "Counted more files than processes:" << endl;
				cout << "Total files: " << total_files << endl;
				cout << "Total procs: " << NPROCS << endl;
				stop_run(1) ;
			}
		}
		if ( total_files == 0 ) 
		{
			if ( RANK == 0 ) 
			{
				cout << "No dimension files were found" << endl ;
				stop_run(1) ;
			}
		}

		// Somewhat tricky logic to allow the number of processors to be greater than the number of files.
		// In that case, only some of the rows in the split A matrix are used.		
		int proc_fac = NPROCS / total_files ;
		int rank_div = RANK / proc_fac ;
		int my_file = ( rank_div > total_files - 1 ) ? total_files - 1 : rank_div ;
		int my_offset = RANK % proc_fac ;

#ifdef VERBOSE		
		cout << "RANK = " << RANK << " proc_fac = " << proc_fac << " my_file = " << my_file << " my_offset = " << my_offset << endl ;
#endif		

		// Read matrix dimensions from the dim.*.txt file.
		sprintf(name, "%s.%04d.%s", strDimFile.c_str(), my_file, dim_ext.c_str()) ;
		dim_file.open(name) ;
		if ( ! dim_file.is_open() ) {
			cerr << "Could not open " + string(name) + "\n" ;
			stop_run(1) ;
		}		
		
		// Dimensions to use if NPROCS == total_files
		int mstart0, mend0, mstore0 ;
		
		dim_file >> dim2 >> mstart0 >> mend0 >> dim1 ;
		if ( ! dim_file.good() ) {
			cerr << "Error reading dim file\n" ;
			stop_run(1) ;
		}
		dim_file.close() ;
		
		if ( mstart0 <= mend0 ) {
			// mstore0 is the number of rows in the file.
			mstore0 = (mend0 - mstart0 + 1) ;
		} else {
			mstore0 = 0 ;
		}

		row_start = mstart0 + my_offset * (mstore0 / proc_fac) ;
		
		if ( rank_div < total_files ) {
			if ( my_offset == proc_fac - 1 )
				row_end = mend0 ;
			else
				row_end = row_start + (mstore0 / proc_fac) - 1 ;
			num_rows= row_end - row_start + 1 ;
		} else {
			// Left-over process not used.
			row_start = dim1 + 1 ;
			row_end  = dim1 ;
			num_rows= 0 ;
		}

#ifdef VERBOSE		
		cout << "RANK = " << RANK << " row_start = " << row_start << " row_end = " << row_end << " num_rows = " << num_rows << endl ;
#endif		
		
		// Append the processor number to the A matrix name.
		char matFilename2[80] ;
		string str_filename(matFilename) ;
		std::size_t found = str_filename.find(".") ;
		if ( found == string::npos ) {
			cerr << "A matrix file name must end with a suffix" ;
			stop_run(1) ;
		}
		
		string mat_ext = str_filename.substr(found+1) ;
		str_filename = str_filename.substr(0,found+1) ;
		sprintf(matFilename2, "%s%04d.%s", str_filename.c_str(), my_file, mat_ext.c_str()) ;
		ifstream matfile(matFilename2);		
		if (!matfile.good()) {
			cerr << "error opening matrix file " << matFilename2 << endl;
			stop_run(1);
		}

		if ( mat != NULL ) {
			delete [] mat ;
		}
		mat = new double[num_rows * dim2];

		if ( shift != NULL ) {
			delete [] shift ;
		}
		shift = new double[dim2] ;

		if ( scale != NULL ) {
			delete [] scale ;
		}
		scale = new double[dim2] ;
		
		for ( int i = 0 ; i < dim2 ; i++ ) {
			shift[i] = 0.0 ;
			scale[i] = 1.0 ;
		}
		
		for (int i= mstart0 ; i <= mend0 ; i++) {
			for (int j=0; j< dim2 ; j++) {
				double val;
				matfile >> val;
				if ( i >= row_start && i <= row_end ) 
					set(i, j, val) ;
			}
		}
		if ( ! matfile.good() ) {
			cerr << "Error reading A matrix" ;
			stop_run(1) ;
		}
		matfile.close();
	}

	void realloc(int d1, int d2) 
		{
			if ( dim1 > 0 && dim2 > 0 ) {
				delete [] mat ;
				delete [] scale ;
				delete [] shift ;
			}
			scale = new double[d2] ;
			shift = new double[d2] ;
			if ( ! distributed ) {
				row_start = 0 ;
				row_end = d1 - 1 ;
				num_rows = row_end - row_start + 1 ;
			} else if ( d1 != dim1 ) {
				// Do not change distribution parameters unless dim1 changes.
				row_start = (d1 * RANK) / NPROCS ;
				if ( RANK < NPROCS - 1 ) {
					row_end = (d1 * (RANK+1)) / NPROCS - 1 ;
				} else {
					row_end = d1 - 1 ;
				}
				num_rows = row_end - row_start + 1 ;
			}
			dim1 = d1 ;
			dim2 = d2 ;
			mat = new double[num_rows * d2] ;
		}
	inline double get(int i, int j) const
		{
#ifdef DEBUG						
			if ( i >= dim1 || j >= dim2 ) {
				cout << "Matrix out of bounds" << endl ;
				stop_run(1) ;
			}
			if ( distributed && (i > row_end || i < row_start) ) {
				cout << "Distributed matrix out of bounds " << endl ;
				stop_run(1) ;
			}
#endif						
			return(mat[(i-row_start) * dim2 + j]) ;
		}
	inline void set(int i, int j, double val) 
		{
#ifdef DEBUG						
			if ( i >= dim1 || j >= dim2 ) {
				cout << "Matrix set out of bounds" << endl ;
				stop_run(1) ;
			}
			if ( distributed && (i > row_end || i < row_start) ) {
				cout << "Distributed matrix set out of bounds " << endl ;
				stop_run(1) ;
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
	void distribute(const Matrix &xin)
	// Settings to distribute a matrix across processes based on another matrix's distribution pattern.
		{
			dim1 = xin.dim1 ;
			row_start = xin.row_start ;
			row_end = xin.row_end ;
			distributed = xin.distributed ;
			num_rows = xin.num_rows ;
		}
	void scale_rows(const Vector& vals)
		// Multiply each row of the matrix by the values in vals.
	{
		if ( vals.dim != dim1 ) {
			cout << "Error in scale_rows: dimensions did not match" << endl ;
		}

#ifdef USE_OPENMP		
#pragma omp parallel for shared(vals) default(none)
#endif		
		for ( int j = row_start ; j <= row_end ; j++ ) {
			for ( int k = 0 ; k < dim2 ; k++ ) {
				set(j,k, get(j,k) * vals.get(j)) ;
			}
		}
		
	}
	double mult_T(int j, int k)
	// Calculate the j,k element of X^T * X, where X is the current matrix.
		{
			double tmp = 0.0 ;
#ifdef USE_OPENMP		
#pragma omp parallel for shared(k,j) reduction(+:tmp) default(none)
#endif					
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
			Vector tmp(dim2, 0.0) ;
			
			for ( int j = 0 ; j < dim2 ; j++ ) {
				double sum = 0.0 ;
#ifdef USE_OPENMP		
#pragma omp parallel for shared(j) reduction(+:sum) default(none)
#endif						
				for ( int k = row_start ; k <= row_end ; k++ ) {
					sum += get(k, j) ;
				}

				tmp.set(j,sum) ;
			}
#ifdef USE_MPI
			MPI_Allreduce(tmp.vec, shift, dim2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else			
			for ( int j = 0 ; j < dim2 ; j++ ) {
				shift[j] = tmp.get(j) ;
			}
#endif			
			
			for ( int j = 0 ; j < dim2 ; j++ ) {
				shift[j] /= dim1 ;
				for ( int k = row_start ; k <= row_end ; k++ ) {
					double val = get(k, j) - shift[j] ;
					set(k, j, val) ;
				}
			}
			for ( int j = 0 ; j < dim2 ; j++ ) {
				tmp.set(j, 0.0) ;
				for ( int k = row_start ; k <= row_end ; k++ ) {
					tmp.add(j, get(k, j) * get(k,j) );
				}
			}
#ifdef USE_MPI
			MPI_Allreduce(tmp.vec, scale, dim2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else
			for ( int j = 0 ; j < dim2 ; j++ ) {
				scale[j] = tmp.get(j) ;
			}
#endif
			
			for ( int j = 0 ; j < dim2 ; j++ ) {			
				scale[j] = sqrt(scale[j]) ;
				for ( int k = row_start ; k <= row_end ; k++ ) {
					double val = get(k, j) / scale[j] ;
					set(k, j, val) ;
				}
			}
		}
	void check_norm() const
	// Check normalization of matrix to according to Eq. 1.1
		{
			Vector test(dim2, 0.0) ;
			Vector test2(dim2, 0.0) ;
			
			for ( int j = 0 ; j < dim2 ; j++ ) {
				for ( int k = row_start ; k <= row_end ; k++ ) {
					test.add(j, get(k, j) ) ;
				}
			}
#ifdef USE_MPI
			MPI_Allreduce(test.vec, test2.vec, dim2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else
			for ( int j = 0 ; j < dim2 ; j++ ) {
				test2.set(j, test.get(j) ) ;
			}
#endif
			if ( RANK == 0 ) {
				for ( int j = 0 ; j < dim2 ; j++ ) {
					if ( fabs(test2.get(j)) > 1.0e-06 ) {
						cout << "Column " << j << " did not have 0 average" << endl ;
					}
				}
			}
			for ( int j = 0 ; j < dim2 ; j++ ) {
				test.set(j, 0.0) ;
				for ( int k = row_start ; k <= row_end ; k++ ) {
					test.add(j, get(k, j) * get(k,j) ) ;
				}
			}
#ifdef USE_MPI
			MPI_Allreduce(test.vec, test2.vec, dim2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else
			for ( int j = 0 ; j < dim2 ; j++ ) {
			    test2.set(j, test.get(j) ) ;
			}
#endif
			if ( RANK == 0 ) {
         	for ( int j = 0 ; j < dim2 ; j++ ) {
        			if ( fabs(test2.get(j)-1.0) > 1.0e-06 ) {
			          cout << "Column " << j << " did not have norm = 1 " << endl ;
    		       }
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
	
	bool cholesky_remove_row(int id ) {
		// Based on larscpp github package.  Modified so that diagonal elements of the cholesky
		// matrix are always positive.
		// Remove a row from a Cholesky matrix.  Update the dimension of the matrix
		// Use transpose matrix access functions because original code was written in terms of a lower
		// triangular matrix, whereas we use an upper triangular matrix.
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
	
	void cholesky_sub(Vector &x, const Vector &b)
	// Assuming that the current matrix holds a Cholesky decomposition
	// in R^T * R form, calculate the solution vector x given the vector
	// of input values b.
		{
			if ( b.dim != x.dim || b.dim != dim1 || b.dim != dim2 ) {
				cout << "Dimension mismatch" ;
				stop_run(1) ;
			}
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
			
	void print()
		{
			if ( RANK == 0 ) {
				cout << "[" << endl ;
				for ( int j = row_start ; j <= row_end ; j++ ) {
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
				stop_run(1) ;
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
				stop_run(1) ;
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
