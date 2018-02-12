
#ifndef _Quad_h  // Protects against double-inclusion.


struct PAIRS	// NEEDS UPDATING
{
	
	int    PAIRIDX;			// Index for the specific pair type (triplet type)
	string PRPR_NM;			// The order-specific atom pair. For example, OH, never HO, or something similiar "proper name"
	string PAIRTYP;			// Allowed values are CHEBYSHEV, DFTBPOLY, SPLINE, INVERSE_R, LJ, and STILLINGER
	string ATM1TYP;			// Atom chemistry (i.e. C, H, N, O...)
	string ATM2TYP;
	double ATM1CHG;			// Atom partial charge... used when charges are fixed
	double ATM2CHG;
	string CHRGSGN;			// Should the fitted charge on a given atom be negative or positive?
	double ATM1MAS;			// Atomic mass (i.e. 12.011 for C)
	double ATM2MAS;	
	double OLD_S_MINIM;		// Used for smoothing functions when 3B and 2B inner cutoffs differ. This should be the inner cutoffs used for the actual SVD
	double S_MINIM;			// Minimum allowed pair distance for fitting
	double S_MAXIM;			// Maximum allowed pair distance for fitting
	double S_DELTA;			// Fitting "grid" spacing (width)

	int N_CFG_CONTRIB;		// How many configurations actually contribute to fitting this pair??

	double CHEBY_RANGE_LOW;	//	When fitting to Cheby polynomials, pair distances are typically transformed to exist defined on range -1 to 1 (where Cheby poly's are defined), 
	double CHEBY_RANGE_HIGH;//  but under certain circumstances, it can be advantagous to only fit some sub range (i.e. -1 to 0). These variables define the transformation range
							//  and arr, by default, set to -1 and 1 for low and high, respectively.

	int    SNUM;			// Number of fitting parameters for pair ... WHY WOULD THIS BE DIFFERENT FOR DIFFERENT PAIR TYPES?***
	int    SNUM_3B_CHEBY;	// Number of fitting parameters for pair ... WHY WOULD YOU NEED BOTH SNUM AND THIS SPECIAL CHEBY ONE?***
	int    SNUM_4B_CHEBY;	// Number of fitting parameters for pair ... WHY WOULD YOU NEED BOTH SNUM AND THIS SPECIAL CHEBY ONE?***
	
	string CHEBY_TYPE;		// Are distances transformed into inverse-r type or morse-type distances?... or not at all? (default)
	double PENALTY_SCALE;	// For 2B Cheby potentials... "a" in vpenalty = a*(smin-penalty_dist-rlen)^3 ... default value is 1.0e8
	double PENALTY_DIST;	// For 2B Cheby potentials... "penalty_dist" in vpenalty = a*(smin-penalty_dist-rlen)^3 ... default value is 0.01
	double CUBIC_SCALE;		// Factor to multiply to the cubic penalty function, (1-rlen/smax)^3... default value is 1
	
	double LAMBDA;			// Morse lambda for CHEBYSHEV type pairs
	double MIN_FOUND_DIST;	// Minimum distance between pairs
	
	bool   USE_OVRPRMS;		// Should overbonding even be computed pair type
	string OVER_TO_ATM;		// Which atom is overbonding *to* being defined for... for example, overbonding to oxygen
	vector<double> OVRPRMS;	// [0] = P_OVERB; [1] = R_0_VAL; [2] = P_1_VAL; [3] = P_2_VAL; [4] = LAMBDA6
	
   vector<double> NBINS;				// Number of bins to use for ij, ik, and jk distances when building the 3B population histograms 

	FCUT FORCE_CUTOFF;	// "CUBIC" "COSINE" or "SIGMOID" currently supported

  // Only used in force field
	vector<double> 	PARAMS;
	vector<double>	POT_PARAMS;		// Used by splines to compute pressure by integrating spline eq's
	double			PAIR_CHRG;
	
PAIRS(): N_CFG_CONTRIB(0),OVRPRMS(5), NBINS(3) {}	// Just a constructor to allow the size of the OVRPRMS vector to be pre-specified
};


// We are not deriving PAIR_FF from PAIRS, because in that case vector<PAIR_FF> does not inheriting from vector<PAIRS>.
// This leads to awkwardness in functions that want vector<PAIR_FF> as an argument.  If PAIR_FF is simply a typedef,
// vector<PAIR_FF> and vector<PAIRS> can be used interchangably, which improves code reuse.  An alternative would
// be to use vector<PAIRS*>, but we are not using pointers much in the current code.
typedef PAIRS PAIR_FF ;

class CLUSTER 
// Structure for a generic cluster of interacting atoms.
// This is a generic class that specific clusters derive from.
{
public:
  int    INDX;
  int NATOMS ;
  int NPAIRS ;
  vector<string> ATOM_PAIRS; 	// Contains the different atom pair string names
  vector<string> ATOM_NAMES ; // Names of each atom in the quad

  FCUT FORCE_CUTOFF;		// "CUBIC" "COSINE" or "SIGMOID" currently supported
	
  int    N_CFG_CONTRIB;		// How many configurations actually contribute to fitting this cluster ??
	
  vector<double> MIN_FOUND;	// Testing two cases: 1. 
	
  vector<double> S_MAXIM ;	// A unique outer cutoff for 4B interactions... by default, is set to S_MAXIM
  vector<double> S_MINIM ;	// Similar for inner cutoff. This is useful when the 2-body 
  // is refit/extrapolated, thus has a s_min lower than the original fitted value
  // Values need to be specified for each contributing pair
  // [0] -> ij, [1] -> ik, [2] -> il, [3] -> jk, [4] -> jl, [5] -> kl 

  map<vector<int>,int> POP_HIST; // Population histogram that s used to set 3B behavior in unsampled regions
	 
  vector<double> NBINS ;				// Number of bins to use for ij, ik, and jk distances when building the population histograms 
  vector<double> BINWS; 				// Binwidths to use for ij, ik, and jk distances when building the population histograms 
	
  vector<double> 	PARAMS;

	
  int N_TRUE_ALLOWED_POWERS;	// How many UNIQUE sets of powers do we have?
  int N_ALLOWED_POWERS;		// How many total sets of powers do we have?
	
  vector<vector<int> > ALLOWED_POWERS;	// This will keep a list of the allowed polynomial powers for each coefficient
  map<vector<int>,int> ALLOWED_POWERS_MAP ;  // Use a map for searching efficiency of allowed polynomial powers.
  vector<vector<int> > UNIQUE_POWERS ;  // This is a list of unique polynomial powers for each coefficient
  vector<int> EQUIV_INDICES;	// For each set of allowed powers, what is the index of the first equivalent set? For example, for the set (OO, OH, OH), (1,0,1) and (1,1,0) are is equivalent
  vector<int>	PARAM_INDICES;	// For each of the set of allowed powers, what would be the index in the FF? for example, for a set of EQUIV_INDICES {0,0,2,3}, PARAM_INDICES would be {0, 0, 1, 2}

  // Build a list of atom names from the pair names and atom types.
  void atom_names_from_pairs(const vector<string> &atom_types) ;

  void build(int cheby_4b_order) ; // The the ALLOWED_POWERS, etc. for an interaction.
  void store_permutations(vector<int> &unsorted_powers) ; // Store all the permutations.
  void print()  ;  // Print the quad powers and element types.

  // Print special parameters to the header file.
  void print_special(ofstream &header, string QUAD_MAP_REVERSE, string output_mode) ;

  // Print the params file header for a cluster
  void print_header(ofstream &header) ;

  CLUSTER() {} 
  virtual ~CLUSTER() {} 

  // Sets up the histogram for TRIPLETS.  Returns true on success, false otherwise.
  bool init_histogram(vector<struct PAIRS> & pairs, map<string,int>& pair_map) ;

  // Increment the histogram with the given index vector.
  void increment_histogram(vector<int> &index) ;

  // Get the value of the histogram with the given index vector.  Return 0 if no entry is found.
  int get_histogram(vector<int> &index) ;

  // Return the maximum cutoff distance for the specified pair.
  double get_smaxim(PAIRS & FF_2BODY, string TYPE)  ;  

  // Return the minimum inner cutoff distance for the specified pair.
  double get_sminim(PAIRS & FF_2BODY, string TYPE) ;

CLUSTER(int natom, int npair): ATOM_PAIRS(npair), ATOM_NAMES(natom), MIN_FOUND(npair),
	 S_MAXIM(npair), S_MINIM(npair), NBINS(npair), BINWS(npair)
  {
	 NATOMS = natom ;
	 NPAIRS = npair ;
	 N_CFG_CONTRIB = 0 ;
	 N_TRUE_ALLOWED_POWERS = 0 ;
	 for ( int j = 0 ; j < npair ; j++ ) {
		MIN_FOUND[j] = -1 ;
	 }
	 for (int j=0; j < npair ; j++)
	 {
		S_MINIM[j] = -1;
		S_MAXIM[j] = -1;
	 }
	 FORCE_CUTOFF.TYPE = FCUT_TYPE::CUBIC;	
  }
private:
  // Recursively permute atom indices for each element type.
  void permute_atom_indices(int idx, vector<string> names, vector<int> &unsorted_powers, 
									 vector<int> perm, int unique_index, int equiv_index) ;

  // Interate over a power index in constructing the powers vector.
  void build_loop(int indx, int cheby_order, vector<int> powers) ;
};

typedef CLUSTER TRIPLETS ;
typedef CLUSTER QUADRUPLETS ;

class CLUSTER_LIST
// A group of clusters that represents an N-body interaction.
{
public:
  // The number of CLUSTERS in the list.
  int NCLUSTERS ;        

  // A vector containing all clusters.
  // Need to store cluster pointers so that polymorphism with different cluster types works.
  vector<CLUSTER> VEC ;

  // A map from pair types into the VEC index.
  map<string,int> MAP ;

  // A map from the VEC index into the pair types.
  map<int,string> MAP_REVERSE ;

  vector<int> INT_MAP ;

  vector<int> INT_MAP_REVERSE ;

  vector<string> EXCLUDE ;

  void build_pairs(vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP) ;

  int build_all(int cheby_order, vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP) ;

  void exclude() ;

  int make_id_int(vector<int>& index) ;

  void print_min_distances() ;

  static string tuplet_name(int natom, bool plural, bool caps) ;

  // Build the fast integer maps for the cluster list.  Used by the MD code.
  void build_int_maps(vector<string> ATOMTYPE, vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP) ;

  // Read the maps from the force field file (MD only).
  void read_maps(ifstream& paramfile, string line) ;

private:
  void build_pairs_loop(int index, vector<int> atom_index, 
								vector<string> ATOM_CHEMS, vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP, int &count) ;

  // Loop across atom types for each atom in the cluster, and build a corresponding loop.
  void build_int_maps_loop(int index, vector<int> atom_index, vector<string> ATOMTYPE,
									vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP) ;
} ;


// These are no longer inherited from CLUSTER so that CLUSTER_LIST works correctly.
typedef QUADRUPLETS QUAD_FF ;
typedef TRIPLETS TRIP_FF ;


#define _Quad_h
#endif // ifndef _Quad_h

