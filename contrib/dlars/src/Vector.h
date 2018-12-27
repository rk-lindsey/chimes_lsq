
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

	Vector &operator=(const Vector& in) {
		if ( dim !=in.dim ) {
			delete [] vec ;
			vec = new double [in.dim] ;
			dim = in.dim ;
		}
		for ( int j = 0 ; j < dim ; j++ ) {
			vec[j] = in.vec[j] ;
		}
		return *this ;
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
	void scale(Vector &out, const Vector &vals) 
	// Multiply the vector by the given array of values, put result in Out.
		{
			for ( int j = 0 ; j < dim ; j++ ) {
				out.set(j, vals.get(j) * vec[j] ) ;
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

	void push(double val)
		// Add the value to the end of the vector.
	{
		double *newv = new double[dim+1] ;
		for ( int j = 0 ; j < dim ; j++ ) {
			newv[j] = vec[j] ;
		}
		newv[dim] = val ;
		delete [] vec ;
		vec = newv ;
		++dim ;
	}
	void add_mult(const Vector &in, double factor)
		// Set out = out + factor * in
	{
		if ( in.dim != dim ) {
			cout << "Error in add_mult: dim mismatch\n" ;
			exit(1) ;
		}
		for ( int k = 0 ; k < dim ; k++ ) {
			vec[k] += factor * in.get(k) ;
		}
	}
} ;
