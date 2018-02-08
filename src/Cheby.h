void ZCalc_Cheby_Deriv   (JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

void ZCalc_3B_Cheby_Deriv(JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<TRIPLETS> & PAIR_TRIPLETS, 
								  vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, 
								  map<string,int> TRIAD_MAP, NEIGHBORS & NEIGHBOR_LIST);

void ZCalc_4B_Cheby_Deriv(JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, CLUSTER_LIST TRIPS,
								  CLUSTER_LIST QUADS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP ,
								  map<int,int>  INT_QUAD_MAP, NEIGHBORS & NEIGHBOR_LIST) ;

void ZCalc_Cheby_ALL(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, CLUSTER_LIST &TRIPS,
							CLUSTER_LIST &QUADS, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST) ;


//////////////////////////////////////////
//
//	FUNCTION HEADERS -- PRINTING OF POTENTIAL ENERGY SURFACE
//
//////////////////////////////////////////

void Print_Cheby             (vector<PAIR_FF> & FF_2BODY, int ij, string PAIR_NAME, bool INCLUDE_FCUT, bool INCLUDE_CHARGES, bool INCLUDE_PENALTY, string FILE_TAG);
void Print_3B_Cheby          (JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk);
void Print_3B_Cheby_Scan     (JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan);
void Print_Ternary_Cheby_Scan(JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan);
