#ifndef _Cheby_h
#define _Cheby_h

#ifndef CHECK_CHEBY_RANGE 
	#define CHECK_CHEBY_RANGE 1	// (true)
#endif


class Cheby
{
public:
  JOB_CONTROL &CONTROLS;
  FRAME &SYSTEM;
  NEIGHBORS &NEIGHBOR_LIST;
  vector<PAIRS> &FF_2BODY;
  vector<int> &INT_PAIR_MAP;

  double DERIV_CONST;

  // Derivatives of the force with respect to the coefficients - 2 BODY
  void Deriv_2B(A_MAT & FRAME_A_MATRIX);

  // Derivatives of the force with respect to the coefficients - 3 BODY
  void Deriv_3B(A_MAT & FRAME_A_MATRIX, CLUSTER_LIST &TRIPS);

  // Derivatives of the force with respect to the coefficients - 4 BODY
  void Deriv_4B(A_MAT & FRAME_A_MATRIX, int n_3b_cheby_terms, CLUSTER_LIST &QUADS);

  // Calculate the many-body Chebyshev force
  void Force_all(CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS);

  // Calculate the 3-body Chebyshev force
  void Force_3B(CLUSTER_LIST &TRIPS);

  // Calculate the 4-body Chebyshev force
  void Force_4B(CLUSTER_LIST &QUADS);

  // General function to Matches pair indices to the ij. ik, jk type pairs formed from the atoms in the cluster.
  void map_indices(CLUSTER & cluster, vector<string> & atom_type, vector<int> & pair_map);

  // General function to Matches pair indices to the ij. ik, jk type pairs formed from the atoms in the cluster.
  // Based on integer atom indices
  void map_indices_int(CLUSTER & cluster, vector<int> & atom_type_idx, vector<int> & pair_map);

  Cheby(JOB_CONTROL &controls, FRAME &system, NEIGHBORS &neighbor_list, vector<PAIRS> &ff_2body, vector<int> &int_pair_map ) : CONTROLS(controls), SYSTEM(system), NEIGHBOR_LIST(neighbor_list), FF_2BODY(ff_2body), INT_PAIR_MAP(int_pair_map)
  {
	 // Constructor initializes references.
	 DERIV_CONST = FF_2BODY[0].CHEBY_RANGE_HIGH - FF_2BODY[0].CHEBY_RANGE_LOW;	// Ranges should be the same for all types
	 DERIV_CONST /= 2.0; // i.e the width of the default cheby range
  };

  // Print 2-Body potential.
  void Print_2B(int ij, string PAIR_NAME, bool INCLUDE_FCUT, bool INCLUDE_CHARGES, bool INCLUDE_PENALTY, string FILE_TAG);

  // Print 3-Body potential.
  void Print_3B(CLUSTER_LIST &TRIPS, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan);

  // Return the transformation corresponding to the given string.
  // static member functions do not depend on a particular instance of the class.
  static Cheby_trans get_trans_type(string cheby_type);

  // Return a string corresponding to the Chebyshev transformation type.
  static string get_trans_string(Cheby_trans trans);

  static void set_cheby_params(double sminim, double smaxim, double lambda, Cheby_trans cheby_type, double &xminim, double &xmaxim, double &xdiff, double &xavg);


private:

  // Return a range-limited copy of x.
  inline double fix_val(double x);

  // Evaluate Chebyshev polynomials.
  void set_polys(int index, double *Tn, double *Tnd, const double rlen, double x_diff, double x_avg, 
							 double SNUM);

  // Does the cheby distance transformation with pre-calculated limits.										
  inline void transform(double rlen, double x_diff, double x_avg, double lambda, 
								Cheby_trans cheby_type, double & x, double &exprlen);

  // Get the index of an atom pair given the two atom numbers.
  inline int get_pair_index(int a1, int a2, const vector<int> &atomtype_idx, int natmtyp, 
									 const vector<int> &parent);

  // Set the chebyshev power for each atom pair in the triplet.
  inline void set_3b_powers(const TRIPLETS & FF_3BODY, const vector<int> &pair_index, int POWER_SET,
									 int & pow_ij, int & pow_ik, int & pow_jk ) ;
};

#endif // defined(_Cheby_h)
