
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

	IntVector &operator=(const IntVector& in) {
		if ( dim !=in.dim ) {
			delete [] vec ;
			vec = new int [in.dim] ;
			dim = in.dim ;
		}
		for ( int j = 0 ; j < dim ; j++ ) {
			vec[j] = in.vec[j] ;
		}
		return *this ;
	}
		
	int size() const {
		return dim ;
	}
	void set(int i, int val) {
#ifdef DEBUG					
		if ( i >= dim ) {
			cout << "IntVector set out of bounds" << endl ;
			stop_run(1) ;
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
		void read_sparse(ifstream &file)
		// Read a vector in sparse format.
		{
			string line ;

			clear() ;
			getline(file,line) ;

			size_t pos =  line.find('[') ;
			if ( pos == string::npos ) {
				cout << "Did not find '[' character in " + line << endl ;
				stop_run(1) ;
			}
			istringstream istr(line.substr(pos+1)) ;
			istr >> dim ;
			vec = new int[dim] ;
			for ( int i = 0 ; i < dim ; i++ ) {
				vec[i] = 0 ;
			}
			int idx, val ;
			for ( int i = 0 ; i < dim ; i++ ) {
				getline(file, line) ;
				if ( line.find(']') != string::npos ) break ;
				istringstream istr(line) ;
				istr >> idx >> val ;
				if ( idx < dim && idx >= 0 ) 
					set(idx, val) ;
				else {
					cout << "Error reading sparse vector " << endl ;
					cout << line ;
					stop_run(1) ;
				}
			}
		}
		
	int get(int idx) const {
#ifdef DEBUG					
		if ( idx >= dim ) {
			cout << "IntVector index out of bounds" << endl ;
			stop_run(1) ;
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
					if ( abs(vec[j]) > 0.0 ) 
						cout << j << " " << vec[j] << endl ;
				}
				cout << "]" << endl ;
			}
		}

		void print_sparse(ofstream &of) 
		{
			if ( RANK == 0 ) {
				of << "[ " << dim << endl ;
				for ( int j = 0 ; j < dim ; j++ ) {
					if ( abs(vec[j]) > 0.0 ) 
						of << j << " " << vec[j] << endl ;
				}
				of << "]" << endl ;
			}
		}

		void print_all(ostream &of)
		// Print all values.
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

	void push(int val)
		// Add the value to the end of the vector.
	{
		int *newv = new int[dim+1] ;
		for ( int j = 0 ; j < dim ; j++ ) {
			newv[j] = vec[j] ;
		}
		newv[dim] = val ;
		delete [] vec ;
		vec = newv ;
		++dim ;
	}
	void add_mult(const IntVector &in, int factor)
	// Set out = out + factor * in
	{
		if ( in.dim != dim ) {
			cout << "Error in add_mult: dim mismatch\n" ;
			stop_run(1) ;
		}
		for ( int k = 0 ; k < dim ; k++ ) {
			vec[k] += factor * in.get(k) ;
		}
	}
	void clear()
	// Clear all entries of the vector.
	{
		if ( dim > 0 ) 
			delete [] vec ;
		vec = NULL ;
		dim = 0 ;
	}
} ;
