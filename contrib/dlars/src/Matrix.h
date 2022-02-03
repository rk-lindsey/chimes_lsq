#ifdef USE_MKL
#include <mkl.h>
#endif
//#ifdef USE_BLAS
//#include <cblas.h>
//#endif

class Matrix {
public:
	double *mat ;  // Elements of the matrix stored here.
	double *shift ; // Normalization shift of each column ;
	double *scale ; // Normalization scale factor for each column.
	int dim1, dim2 ; 			// dim1 is number of rows.  Dim2 is number of columns.
	bool distributed ;    // Is this matrix distributed over processes ?
	int row_start, row_end ;  // Starting and ending row for parallel calculations.
	int num_rows ;    // Number of rows stored on this process.

	bool cholesky(Matrix &chol) ;
	bool cholesky_distribute(Matrix &chol) ;  
	bool cholesky_add_row(const Matrix &chol0, const Vector &newr)  ;
	void cholesky_sub(Vector &x, const Vector &b) ;
	bool cholesky_remove_row(int id ) ;
	void cholesky_sub_distribute(Vector &x, const Vector &b) ;
	int rank_from_row(int j) {
		for ( int k = 0 ; k < NPROCS ; k++ ) {
			int r_start = (dim1 * k) / NPROCS ;
			int r_end = 0 ;
			// Repeat assignment for each rank, see which works with row j.
			if ( k < NPROCS - 1 ) {
				r_end = (dim1 * (k+1)) / NPROCS - 1 ;
			} else {
				r_end = dim1 - 1 ;
			}
			if ( r_start <= j && j <= r_end ) {
				return k ;
			}
		}
		cout << "Error: Did not find a rank for row " << j << endl ;
		stop_run(1) ;
		return -1 ;
	}
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

		if ( RANK == 0 ) cout << "Assigning " << proc_fac << " MPI processes per file\n" ;
                  
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

		void resize(int d1, int d2)
		// Resize a matrix, preserving its values and setting any new values to 0.
		{
			if ( dim1 == 0 || dim2 == 0 ) {
				realloc(d1,d2) ;
				for ( int j = 0 ; j < d1 ; j++ ) {
					for ( int k = 0 ; k < d2 ; k++ ) {
						set(j,k,0.0) ;
					}
				}
				return ;
			}

			// New dimensions.
			int row_start1, row_end1, num_rows1 ;

			if ( ! distributed ) {
				row_start1 = 0 ;
				row_end1 = d1 - 1 ;
				num_rows1 = row_end1 - row_start1 + 1 ;
			} else if ( d1 != dim1 ) {
				// Do not change distribution parameters unless dim1 changes.
				row_start1 = (d1 * RANK) / NPROCS ;
				if ( RANK < NPROCS - 1 ) {
					row_end1 = (d1 * (RANK+1)) / NPROCS - 1 ;
				} else {
					row_end1 = d1 - 1 ;
				}
				num_rows1 = row_end1 - row_start1 + 1 ;
			} else {
				row_start1 = row_start ;
				row_end1 = row_end ;
				num_rows1 = num_rows ;
			}
			
			double *scale1 = new double[d2] ;
			double *shift1 = new double[d2] ;
			double *mat1   = new double[num_rows1 * d2] ;

			for ( int i = 0 ; i < d1 && i < dim1 ; i++ ) {
				for ( int j = 0 ; j < d2 && j < dim2 ; j++ ) {
					mat1[(i-row_start1) * d2 + j] = get(i,j) ;
				}
			}
			for ( int i = 0 ; i < d1 ; i++ ) {
				for ( int j = 0 ; j < d2 ; j++ ) {
					if ( i >= dim1 || j >= dim2 ) {
						mat1[(i-row_start1) * d2 + j] = 0.0 ;
					}
				}
			}
			dim1 = d1 ;
			dim2 = d2 ;
			row_start = row_start1 ;
			row_end = row_end1 ;
			num_rows = num_rows1 ;

			delete [] mat ;
			delete [] scale ;
			delete [] shift ;

			mat = mat1 ;
			scale = scale1 ;
			shift = shift1 ;
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
	// Settings to distribute a matrix across processes.  Assumes the number of rows have already been set in dim1.
	{
		distribute(dim1) ;
	}
	
	void distribute(int total_rows)
	// Settings to distribute a matrix across processes, given the total number of rows.
		{
			dim1 = total_rows ;
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


         void mult_T(int j, Vector& out)
	// Calculate the jth row of X^T * X, where X is the current matrix.
        // Places the result in out.
		{
                  Vector tmp(dim2, 0.0) ;
#ifdef USE_OPENMP		
#pragma omp parallel for shared(j,tmp) default(none)
#endif					
                  for ( int k = 0 ; k < dim2 ; k++ ) {
                    for ( int l = row_start ; l <= row_end ; l++ ) {
                      tmp.add(k, get(l, j) * get(l, k) );
                    }
                  }

#ifdef USE_MPI
                  MPI_Allreduce(tmp.vec, out.vec, dim2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else // USE_MPI
                  for ( int k = 0 ; k < dim2 ; k++ ) {
                    out.set(k, tmp.get(k)) ;
                  }
#endif // USE_MPI
		}

  
         void mult_T_lower(int j, Vector& out)
	// Calculate the jth row of X^T * X, where X is the current matrix.
        // Only the lower half k <= j are calculated in out[k].
        // Places the result in out.
		{
                  Vector tmp(j+1, 0.0) ;

#ifdef USE_OPENMP		
#pragma omp parallel for shared(j,tmp) default(none)
#endif	// USE_OPENMP				
                  for ( int k = 0 ; k <= j ; k++ ) {
                    for ( int l = row_start ; l <= row_end ; l++ ) {
                      tmp.add(k, get(l, k)  * get(l, j)  );
                    }
                  }
                  
#ifdef USE_MPI
                  MPI_Allreduce(tmp.vec, out.vec, j+1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else // USE_MPI
                  for ( int k = 0 ; k <= j ; k++ ) {
                    out.set(k, tmp.get(k)) ;
                  }
#endif // USE_MPI
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
			
	void print()
		{
			if ( distributed ) {
				if ( RANK == 0 ) cout << "[" << endl ;				
				for ( int proc = 0 ; proc < NPROCS ; proc++ ) {
#ifdef USE_MPI
					cout.flush() ;
					MPI_Barrier(MPI_COMM_WORLD) ;
#endif
					if ( RANK == proc ) {
						for ( int j = row_start ; j <= row_end ; j++ ) {
							cout << "[" << endl ;
							for ( int k = 0 ; k < dim2 ; k++ ) {
								cout << j << " " << k << " " << get(j,k) << endl ;
							}
							cout << "]" << endl ;
						}
					}
				}
				if ( RANK == 0 ) cout << "]" << endl ;								
			} else {
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
	
	void dot(Vector &out, const Vector &in) const
	// Find matrix * in = out
		{
			if ( out.dim != dim1 || in.dim != dim2 ) {
				cout << "Array dimension mismatch" << endl ;
				stop_run(1) ;
			}
			if ( ! distributed ) {
#ifdef USE_BLAS			
				cblas_dgemv(CblasRowMajor, CblasNoTrans, dim1, dim2, 1.0,
							mat, dim2, in.vec, 1, 0.0, out.vec, 1) ;
#else				
				for ( int j = 0 ; j < dim1 ; j++ ) {
					double sum = 0.0 ;
					for ( int k = 0 ; k < dim2 ; k++ ) {
						sum += get(j,k) * in.get(k) ;
					}
					out.set(j, sum) ;
				}
#endif			
			} else {
#ifdef USE_BLAS
				// Perform matrix-vector multiply on just the rows owned by this process.
				// Put results into offset indices of out.vec[].
				cblas_dgemv(CblasRowMajor, CblasNoTrans, num_rows, dim2, 1.0,
							mat, dim2, in.vec, 1, 0.0, out.vec + row_start, 1) ;
#else				
				for ( int j = row_start ; j <= row_end ; j++ ) {
					double sum = 0.0 ;
					for ( int k = 0 ; k < dim2 ; k++ ) {
						sum += get(j,k) * in.get(k) ;
					}
					out.set(j, sum) ;
				}
#endif // USE_BLAS
				
#ifdef USE_MPI
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

	void cholesky_invert(Matrix &out)
	// Invert a cholesky decomposed matrix, storing the result in out
	// This is done by solving a linear equation n times.
	{
		if ( out.dim2 != dim2 || out.dim1 != dim1 || dim1 != dim2 ) {
			cout << "Array dimension mismatch" << endl ;
			stop_run(1) ;
		}
		for ( int j = 0 ; j < dim1 ; j++ ) {
			Vector b(dim1, 0.0) ;
			Vector x(dim1, 0.0) ;
			b.set(j, 1.0) ;
			cholesky_sub(x, b) ;

			for ( int k = 0 ; k < dim1 ; k++ ) {
				out.set(k,j, x.get(k) ) ;
			}
		}
	}
		
	void dot_transpose(Vector &out, const Vector &in)
	// Find Transpose(matrix) * in = out
		{
			if ( out.dim != dim2 || in.dim != dim1 ) {
				cout << "Array dimension mismatch" << endl ;
				stop_run(1) ;
			}
			if ( ! distributed ) {
#ifdef USE_BLAS
				// Perform matrix-vector multiply on all rows.
				cblas_dgemv(CblasRowMajor, CblasTrans, dim1, dim2, 1.0,
							mat, dim2, in.vec, 1, 0.0, out.vec, 1) ;
#else								
				Vector sumv(dim2,0.0) ;			
				for ( int j = 0 ; j < dim2 ; j++ ) {
					for ( int k = 0 ; k < dim1 ; k++ ) {
						sumv.add(j, get(k,j) * in.get(k) ) ;
					}
				}
				for ( int j = 0 ; j < dim2 ; j++ ) {
					out.set(j, sumv.get(j)) ;
				}
#endif				
			} else {
				Vector sumv(dim2,0.0) ;			
#ifdef USE_BLAS
				// Perform matrix-vector multiply on only rows owned by this process.
				// The in.vec needs to be offset.
				cblas_dgemv(CblasRowMajor, CblasTrans, num_rows, dim2, 1.0,
							mat, dim2, in.vec + row_start, 1, 0.0, sumv.vec, 1) ;
#else				
				for ( int j = 0 ; j < dim2 ; j++ ) {
					for ( int k = row_start ; k <= row_end ; k++ ) {
						sumv.add(j, get(k,j) * in.get(k) ) ;
					}
				}
#endif // USE_BLAS
				
#ifdef USE_MPI
				MPI_Allreduce(sumv.vec, out.vec, dim2,
							  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD) ;
#else				
				for ( int j = 0 ; j < dim2 ; j++ ) {
					out.set(j, sumv.get(j)) ;
				}					
#endif // USE_MPI
			}
		}

	bool con_grad(Vector &x, const Vector &b, int max_iter, int max_restart,
								double tol)
	// Congugate gradient method to solver linear equations.
	// Solves Matrix * x = b.  The value given in x is used
	// as an initial guess.
	// max_iter is the maximum number of iterations to perform.
	// max_restart is the maximum number of restarts to perform.
	// tol is the tolerance.  The solution is returned in x.
		{
			int dim = x.dim ;
			
			Vector r2(dim), r1(dim) ;
			Vector x2(dim), x1(dim) ;
			Vector p2(dim), p1(dim) ;
			Vector tmp(dim), xstep(dim), tmp2(dim) ;

			int iter = 0 ;
			for ( int l = 0 ; l < max_restart ; l++ ) {
				double oldfunc = 1.0e50 ;

				dot(tmp, x) ;
				r1.assign_mult(b, tmp, -1.0) ;

				p1 = r1 ;
				x1 = x ;

				
				for ( int k = 0 ; k < max_iter ; k++ ) {

					++iter ;
					// Matrix dot vector.
					dot(tmp, p1) ;

					// Vector dot vector.
					double denom = p1.dot(tmp) ;
					if ( fabs(denom) < 1.0e-50 ) {
						if ( RANK == 0 ) cout << "con_grad failed: zero denominator\n" ;
						return false ;
					}
					double alpha = r1.dot(r1) / denom ;

					x2.assign_mult(x1, p1, alpha) ;
					r2.assign_mult(r1, tmp, -alpha) ;

					p1.scale(xstep,alpha) ;
					double err =  sqrt( r2.dot(r2) ) ;
					double stp = sqrt(xstep.dot(xstep)) ;

					dot(tmp2, x2) ;

					// func is the objective function for minimization, which should decrease.
					double func = 0.5 * tmp2.dot(x2) - x2.dot(b) ;
					
					if ( RANK == 0 ) { 
						cout << " Iter = " << iter << " Error = " << err << " Func = " << func << " Step = " << stp << endl ; 
					}
					
					if ( err < tol ) {
//						if ( RANK == 0 ) cout << "con_grad succeeded in " << k+1 << " iterations" << endl ;
						x = x2 ;
						return true ;
					}
					/* if ( func > oldfunc ) { */
					/* 	if ( RANK == 0 ) { */
					/* 		cout << "Failed to decrease conjugate gradient function\n" ; */
					/* 		cout << "Restarting\n" ; */
					/* 		break ; */
					/* 	} */
					/* } */

					// Fletcher-Reeves formula.
					double beta = r2.dot(r2) / r1.dot(r1) ;

					// Polak-Ribiere formula.
          // double beta = (r2.dot(r2) - r2.dot(r1)) / r1.dot(r1) ;

					p2.assign_mult(r2,p1,beta) ;

					// Reset for net iterations.
					x1 = x2 ;
					r1 = r2 ;
					p1 = p2 ;
					oldfunc = func ;
					
				}

				// Reset x for restart.
				// if ( RANK == 0 ) cout << "Restarting conjugate gradient " << endl ;
				x = x2 ;
				
			}
			if ( RANK == 0 ) {
				cout << "con_grad failed in " << iter << " iterations " << endl ;
			}
			return false ;
		} // End of con_grad


	bool pre_con_grad(Vector &x, const Vector &b, const Matrix &M_inv, int max_iter, int max_restart,
								double tol)
	// Preconditioned congugate gradient method to solver linear equations.
	// Solves Matrix * x = b.  The value given in x is used
	// as an initial guess.
	// The inverse of the pre-conditioner, M_inv, is explicitly specified, 
	// as opposed to specifying M.
	// max_iter is the maximum number of iterations to perform.
	// max_restart is the maximum number of restarts to perform.
	// tol is the tolerance.  The solution is returned in x.
		{
			int dim = x.dim ;
			
			Vector r2(dim), r1(dim) ;
			Vector x2(dim), x1(dim) ;
			Vector p2(dim), p1(dim) ;
			Vector tmp(dim), xstep(dim), tmp2(dim) ;
			Vector z2(dim), z1(dim) ;

			int iter = 0 ;
			for ( int l = 0 ; l < max_restart ; l++ ) {
				double oldfunc = 1.0e50 ;

				dot(tmp, x) ;
				r1.assign_mult(b, tmp, -1.0) ;

				x1 = x ;
				M_inv.dot(z1, r1) ;
				p1 = z1 ;
				
				for ( int k = 0 ; k < max_iter ; k++ ) {

					++iter ;
					// Matrix dot vector.
					dot(tmp, p1) ;

					// Vector dot vector.
					double denom = p1.dot(tmp) ;
					if ( fabs(denom) < 1.0e-50 ) {
						if ( RANK == 0 ) cout << "pre_con_grad failed: zero denominator\n" ;
						return false ;
					}
					double alpha = r1.dot(z1) / denom ;
					
					x2.assign_mult(x1, p1, alpha) ;
					r2.assign_mult(r1, tmp, -alpha) ;

					p1.scale(xstep,alpha) ;
					double err =  sqrt( r2.dot(r2) ) ;
					double stp = sqrt(xstep.dot(xstep)) ;

					dot(tmp2, x2) ;

					// func is the objective function for minimization, which should decrease.
					double func = 0.5 * tmp2.dot(x2) - x2.dot(b) ;
					
					if ( RANK == 0 ) { 
						cout << " Iter = " << iter << " Error = " << err << " Func = " << func << " Step = " << stp << endl ; 
					}
					
					if ( err < tol ) {
						if ( RANK == 0 ) cout << "pre_con_grad succeeded in " << k+1 << " iterations" << endl ;
						x = x2 ;
						return true ;
					}

					M_inv.dot(z2, r2) ;

					double beta = z2.dot(r2) / z1.dot(r1) ;

					p2.assign_mult(z2,p1,beta) ;

					// Reset for net iterations.
					x1 = x2 ;
					r1 = r2 ;
					p1 = p2 ;
					z1 = z2 ;
					oldfunc = func ;
					
				}

				// Reset x for restart.
				if ( RANK == 0 ) cout << "Restarting preconditioned conjugate gradient " << endl ;
				x = x2 ;
				
			}
			if ( RANK == 0 ) {
				cout << "pre_con_grad failed in " << iter << " iterations " << endl ;
			}
			return false ;
		} // End of pre_con_grad

	double memory()
	// Returns memory used on the current rank in MB.
	{
		return( (double) (row_end - row_start + 1) * dim2 * 8 / (1024.0 * 1024.0) ) ;
	}

        double print_memory(string name) 
        {
          double mem = memory() ;

          int prec = cout.precision() ;
          cout << "Memory in " << name << " = " << std::fixed << std::setprecision(1) << mem << " MB " << endl ;
          cout.precision(prec) ;
          cout << std::scientific ;
          
          return mem ;
        }
  
} ;
