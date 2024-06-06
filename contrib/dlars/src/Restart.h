
class Restart {
	// Helper functions for restarting.
public:
	static void read_vector(ifstream &inf, string name, Vector &vec) ;
	static void read_int_vector(ifstream &inf, string name, IntVector &vec) ;
	static void read_matrix(ifstream &inf, string name, Matrix &mat, int dim1,
							int dim2, bool distributed) ;
    template<typename T>
    static void read_scalar(ifstream &inf, string name, T &val) ;
} ;


template<typename T>
void Restart::read_scalar(ifstream &inf, string name, T &val)
{
	string line ;
	inf >> line ;

	if ( line.find(name) != string::npos ) {
		inf >> val ;
		if ( RANK == 0 ) cout << "Set " << name << " = " << val << " from restart file" << endl ;
	} else {
		if ( RANK == 0 ) cout << line ;
		cout << "Could not read " << name << " from restart file" << endl ;
		stop_run(1) ;
	}
}
