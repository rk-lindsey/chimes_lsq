class Cheby
{
public:
  JOB_CONTROL &CONTROLS ;
  FRAME &SYSTEM ;
  NEIGHBORS &NEIGHBOR_LIST ;
  vector<PAIRS> &FF_2BODY ;
  map<string,int> &PAIR_MAP ;
  vector<int> &INT_PAIR_MAP ;

  // Derivatives of the force with respect to the coefficients - 2 BODY
  void Deriv_2B(vector<vector <XYZ > > & FRAME_A_MATRIX) ;

  // Derivatives of the force with respect to the coefficients - 3 BODY
  void Deriv_3B(vector<vector <XYZ > > & FRAME_A_MATRIX, CLUSTER_LIST &TRIPS) ;

  // Derivatives of the force with respect to the coefficients - 4 BODY
  void Deriv_4B(vector<vector <XYZ > > & FRAME_A_MATRIX, int n_3b_cheby_terms, CLUSTER_LIST &QUADS) ;

  void Force_all(CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS) ;

Cheby(JOB_CONTROL &controls, FRAME &system, NEIGHBORS &neighbor_list, vector<PAIRS> &ff_2body, map<string,int> &pair_map, vector<int> &int_pair_map ) :
  CONTROLS{controls}, SYSTEM{system}, NEIGHBOR_LIST{neighbor_list}, FF_2BODY{ff_2body}, PAIR_MAP{pair_map}, INT_PAIR_MAP{int_pair_map} 
  {
	 // Constructor initializes references.
  } ;

  // Print 2-Body potential.
  void Print_2B(int ij, string PAIR_NAME, bool INCLUDE_FCUT, bool INCLUDE_CHARGES, bool INCLUDE_PENALTY, string FILE_TAG) ;

  // Print 3-Body potential.
  void Print_3B(CLUSTER_LIST &TRIPS, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, 
					 int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan)	;

} ;

