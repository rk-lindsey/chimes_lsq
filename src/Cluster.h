
#ifndef _Cluster_h  // Protects against double-inclusion.


#include "../imports/chimes_calculator/serial_interface/src/serial_chimes_interface.h"

class PAIRS	// NEEDS UPDATING
{
	public:
	  	int    PAIRIDX; 		      // Index for the specific pair type (triplet type)
	  	string PRPR_NM; 		      // The order-specific atom pair. For example, OH, never HO, or something similiar "proper name"
	  	string PAIRTYP; 		      // Only allowed values is CHEBYSHEV or LJ
	  	string ATM1TYP; 		      // Atom chemistry (i.e. C, H, N, O...)
	  	string ATM2TYP;
	  	int    ATM1TYPE_IDX;		      // The type index of atom 1
	  	int    ATM2TYPE_IDX;		      // The type index of atom 2
	  	double ATM1CHG; 		      // Atom partial charge... used when charges are fixed
	  	double ATM2CHG;
	  	string CHRGSGN; 		      // Should the fitted charge on a given atom be negative or positive?
	  	double ATM1MAS; 		      // Atomic mass (i.e. 12.011 for C)
	  	double ATM2MAS;       
	  	double S_MINIM; 		      // Minimum allowed pair distance for fitting
	  	double S_MAXIM; 		      // Maximum allowed pair distance for fitting
		
		double KILLLEN;			      // Pair distance below which to kill the simulation

	  	double X_MINIM; 		      // Minimum transformed Cheby
	  	double X_MAXIM; 		      // Maximum transformed Cheby
	  	double X_AVG;			      // Average of transformed Cheby limits.
	  	double X_DIFF;  		      // Difference between transformed cheby limits.

	  	int N_CFG_CONTRIB;		      // How many configurations actually contribute to fitting this pair??
		      
	  	double CHEBY_RANGE_LOW; 	      //      When fitting to Cheby polynomials, pair distances are typically transformed to exist defined on range -1 to 1 (where Cheby poly's are defined), 
	  	double CHEBY_RANGE_HIGH;	      //  but under certain circumstances, it can be advantagous to only fit some sub range (i.e. -1 to 0). These variables define the transformation range
	  					      //  and arr, by default, set to -1 and 1 for low and high, respectively.

	  	int    SNUM;			      // Number of fitting parameters for pair ... WHY WOULD THIS BE DIFFERENT FOR DIFFERENT PAIR TYPES?***
	  	int    SNUM_3B_CHEBY;		      // Number of fitting parameters for pair ... WHY WOULD YOU NEED BOTH SNUM AND THIS SPECIAL CHEBY ONE?***
	  	int    SNUM_4B_CHEBY;		      // Number of fitting parameters for pair ... WHY WOULD YOU NEED BOTH SNUM AND THIS SPECIAL CHEBY ONE?***
		      
	  	Cheby_trans CHEBY_TYPE; 	      // Are distances transformed into inverse-r type or morse-type distances?... or not at all? (default)
	  	double PENALTY_SCALE;		      // For 2B Cheby potentials... "a" in vpenalty = a*(smin-penalty_dist-rlen)^3 ... default value is 1.0e8
	  	double PENALTY_DIST;		      // For 2B Cheby potentials... "penalty_dist" in vpenalty = a*(smin-penalty_dist-rlen)^3 ... default value is 0.01
		      
	  	double LAMBDA;  		      // Morse lambda for CHEBYSHEV type pairs
	  	double MIN_FOUND_DIST;  	      // Minimum distance between pairs
		      	      
	  	FCUT FORCE_CUTOFF;		      // "CUBIC" "COSINE" or "SIGMOID" currently supported

	  	// Only used in force field
	  	vector<double>        PARAMS;
	  	vector<double>        POT_PARAMS;     // Used by splines to compute pressure by integrating spline eq's
	  	double  	      PAIR_CHRG;

		serial_chimes_interface chimes ;
		string serial_chimes_file ;
	 
	 PAIRS() : chimes(false)
	  	{
			KILLLEN = 0.0;
			N_CFG_CONTRIB = 0;
			SNUM = 0;
			SNUM_3B_CHEBY = 0;
			SNUM_4B_CHEBY = 0;
			CHEBY_RANGE_HIGH = 1.0;
			CHEBY_RANGE_LOW  = -1.0;

			// Set simple defaults.
			CHEBY_TYPE = Cheby_trans::NONE ;
			X_MINIM = 0.8 ;
			X_MAXIM = 5.0 ;
			X_AVG = (X_MINIM + X_MAXIM)/2.0 ;
			X_DIFF = (X_MAXIM-X_MINIM) ;
			S_MINIM = X_MINIM ;
			S_MAXIM = X_MAXIM ;
				 
	  	}     // Just a constructor to set some defaults

	  	// Set Chebyshev min/max vals.
	  	void set_cheby_vals();
};


// We are not deriving PAIR_FF from PAIRS, because in that case vector<PAIR_FF> does not inheriting from vector<PAIRS>.
// This leads to awkwardness in functions that want vector<PAIR_FF> as an argument.  If PAIR_FF is simply a typedef,
// vector<PAIR_FF> and vector<PAIRS> can be used interchangably, which improves code reuse.  An alternative would
// be to use vector<PAIRS*>, but we are not using pointers much in the current code.
typedef PAIRS PAIR_FF;

class CLUSTER 
// Structure for a generic cluster of interacting atoms.
// This is a generic class that specific clusters derive from.
{

public:
  int INDX;
  int NATOMS;
  int NPAIRS;


  vector<string> ATOM_PAIRS; 	// Contains the different atom pair string names
  vector<string> ATOM_NAMES; 	// Names of each atom in the cluster

  
  vector<int> ATOM_INDICES;  	// Index (element type) of each atom

  FCUT FORCE_CUTOFF;		// "CUBIC" "COSINE" or "SIGMOID" currently supported
	
  int N_CFG_CONTRIB;		// How many configurations actually contribute to fitting this cluster ??
	
  bool EXCLUDED;          	// If true, this is an excluded interaction.

  vector<double> MIN_FOUND;	// Testing two cases: 1. 
	
  vector<double> S_MAXIM;	// A unique outer cutoff for 4B interactions... by default, is set to S_MAXIM
  vector<double> S_MINIM;	// Similar for inner cutoff. This is useful when the 2-body 

  vector<double> X_MINIM;	// Minimum transformed Cheby
  vector<double> X_MAXIM;	// Maximum transformed Cheby
  vector<double> X_AVG;		// Average of transformed Cheby limits.
  vector<double> X_DIFF;	// Difference between transformed cheby limits.

  bool SPECIAL_S_MINIM;		// Determine whether s_maxim was set via user input.	  
  bool SPECIAL_S_MAXIM;		// Determine whether s_maxim was set via user input.
  				// Values need to be specified for each contributing pair
  				// [0] -> ij, [1] -> ik, [2] -> il, [3] -> jk, [4] -> jl, [5] -> kl 
	
  vector<double> PARAMS;
	
  int N_TRUE_ALLOWED_POWERS;	// How many UNIQUE sets of powers do we have?
  int N_ALLOWED_POWERS;		// How many total sets of powers do we have?
	
  vector<vector<int> > ALLOWED_POWERS;		// This will keep a list of the allowed polynomial powers for each coefficient
  map<vector<int>,int> ALLOWED_POWERS_MAP;	// Use a map for searching efficiency of allowed polynomial powers.
  vector<vector<int> > UNIQUE_POWERS;		// This is a list of unique polynomial powers for each coefficient
  vector<int>          EQUIV_INDICES;		// For each set of allowed powers, what is the index of the first equivalent set? For example, for the set (OO, OH, OH), (1,0,1) and (1,1,0) are is equivalent
  vector<int>	       PARAM_INDICES;		// For each of the set of allowed powers, what would be the index in the FF? for example, for a set of EQUIV_INDICES {0,0,2,3}, PARAM_INDICES would be {0, 0, 1, 2}

  vector<int> POWER_COUNT;  // This counts how many times a set of powers occurs due to permutations.

	////////////////////////
  	//    MEMBER FUNCTIONS
	////////////////////////

  void build(int cheby_4b_order); 				// The the ALLOWED_POWERS, etc. for an interaction.
  void store_permutations(vector<int> &unsorted_powers); 	// Store all the permutations.
  void print(bool md_mode) const;				// Print parameters for the cluster.  If md_mode is true, print the potential parameters.

  void print_special (ofstream &header, string QUAD_MAP_REVERSE, string output_mode);	// Print special parameters to the header file.
  void print_header  (ofstream &header);					// Print the params file header for a cluster
  void read_ff_params(ifstream &paramfile, const vector<string> &atomtype);	// Read the force field parameters for a cluster.

  CLUSTER() {} 
  virtual ~CLUSTER() {} 

  
  int match_atom_type_idx(string atm_typ);			// Return the atom type index of an atom in the cluster with the given name.
  inline double get_smaxim(string TYPE);			// Return the maximum cutoff distance for the specified pair.
  inline double get_sminim(string TYPE);			// Return the minimum inner cutoff distance for the specified pair.

  void set_default_smaxim(const vector<PAIRS> & FF_2BODY);	// Set defaults for outer polynomial ranges
  void set_default_sminim(const vector<PAIRS> & FF_2BODY);	// Set defaults for inner polynomial ranges

  void set_cheby_vals(vector<PAIRS> &FF_2BODY);			// Calculate Chebyshev xmin, xmax, xavg.

  void set_atom_indices(vector<string> ATOMTYPE, vector<int> ATOMTYPE_IDX);	// Set the ATOM_INDICES based on the given atomtype names and atomtype_idx indices.

CLUSTER(int natom, int npair): 
	ATOM_PAIRS(npair), 
	ATOM_NAMES(natom), 
	ATOM_INDICES(natom,-1), 
	MIN_FOUND(npair,1.0e10),
	S_MAXIM(npair,-1), 
	S_MINIM(npair,-1), 
    	X_MINIM(npair), 
	X_MAXIM(npair), 
	X_AVG(npair), 
	X_DIFF(npair)
  {
	 EXCLUDED = false;
	 SPECIAL_S_MAXIM = false;
	 SPECIAL_S_MINIM = false;
	 NATOMS = natom;
	 NPAIRS = npair;
	 N_CFG_CONTRIB = 0;
	 N_TRUE_ALLOWED_POWERS = 0;
	 FORCE_CUTOFF.TYPE = FCUT_TYPE::CUBIC;	
  }

  
  void map_indices_int(vector<int> & atom_type_idx, vector<int> & pair_map);	// Matches the allowed powers to the ij, ik, il... type pairs formed from the atoms  ai, aj, ak, al

private:
  // Recursively permute atom indices for each element type.
  void permute_atom_indices(int idx, vector<string> names, vector<int> &unsorted_powers, vector<int> perm, int unique_index, int equiv_index);

  // Interate over a power index in constructing the powers vector.
  void build_loop(int indx, int cheby_order, vector<int> powers);

};

typedef CLUSTER TRIPLETS;
typedef CLUSTER QUADRUPLETS;

class CLUSTER_LIST
// A group of clusters that represents an N-body interaction.
{
public:

  int 			NCLUSTERS;		// The number of CLUSTERS in the list.
  double 		MAX_CUTOFF ;   		// The maximum interaction cutoff for all clusters.
  vector<string> 	S_MAXIM_INPUT ;		// Input for special smaxim values.
  vector<string> 	S_MINIM_INPUT ; 	// Input for special sminim values.
  vector<CLUSTER> 	VEC;			// A vector containing all clusters ... Need to store cluster pointers so that polymorphism with different cluster types works.
  map<string,int> 	MAP;			// A map from pair types into the VEC index.
  map<int,string> 	MAP_REVERSE;		// A map from the VEC index into the pair types.
  vector<int> 		INT_MAP;
  vector<vector<int>> 	PAIR_INDICES;		// Store indices for ordering of pair properties (Cheby powers, s_minim, etc.)  corresponding to a particular atom ordering.
  vector<int> 		INT_MAP_REVERSE;

  vector<vector<string>> EXCLUDE;

  	////////////////////////
  	//    MEMBER FUNCTIONS
	////////////////////////

  // Build the pair variables for all of the clusters
  void build_pairs(vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP, vector<string> atom_types, vector<int> atom_typeidx);

  // Build all triplets and associated maps for the cluster list.
  int build_all(int cheby_order, vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP,  vector<string> atom_types, vector<int> atomtype_idx);
  
  // Set up the cheby variables (xavg, xdiff, etc)
  void build_cheby_vals(vector<PAIRS> & ATOM_PAIRS);

  // Returns a unique ID number for a cluster of atoms.
  int make_id_int(vector<int>& index);

  void print(bool md_mode);

  // Print minimum distances found for the cluster.
  void print_min_distances();

  // print the force field file header for the cluster list.
  void print_header(ofstream &header, int natoms, int cheby_order);

  // Print the cutoff function parameters controlling the CLUSTER_LIST.
  void print_fcut();

  // Print the cutoff function parameters to the force field header file.
  void print_fcut_header(ostream &header);

  // Allocate the cluster list according to the number of clusters and the number of atoms.  Return a string to search for in the params file that describes the cluster list.
  string allocate(int nclusters, int natoms, const vector<PAIRS> &FF_2BODY);

  // Read the excluded interactions from the input stream.
  void read_exclude(class InputListing & CONTENTS, int lineno);

  static string tuplet_name(int natom, bool plural, bool caps);

  // Build the fast integer maps for the cluster list.  Used by the MD code.
  void build_int_maps(vector<string> ATOMTYPE, vector<int> ATOMTYPE_IDX, vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP);

  // Read the maps from the force field file (MD only).
  void read_maps(ifstream& paramfile, string line);

// Read the force field parameters for a cluster list.
  void read_ff_params(ifstream &PARAMFILE, const vector<string>& TMP_ATOMTYPE);


  // Overloaded for different handling in lsq and MD codes

  // Read smaxim, sminim for the cluster list and store the input.
  void read_cutoff_params(class InputListing & CONTENTS, int lineno, string input_type);

  // Read smaxim, sminim for the cluster list and store the input.
   void read_cutoff_params(istream &input, string LINE, string input_type) ;


  // Process cutoff parameter values from stored input.
  void process_cutoff_params(string input_type,vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP) ;
  
  // Update min/max cutoff values. Required for proper handling of special cutoffs
  void update_minmax_cutoffs(const vector<PAIRS> & FF_2BODY); 

  // Parse the force cutoff parameters for a cluster.
  void parse_fcut(string LINE);

  // Print out special 3 and 4 body force parameters.
  void print_special(ofstream &header);

  // Set defaults for force cutoffs.
  void set_default_cutoffs(const vector<PAIRS>& FF_2BODY);

private:
  void build_pairs_loop(int index, vector<int> atom_index, vector<string> ATOM_CHEMS, vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP, int &count, vector<string> atom_types, vector<int> atom_typeidx);

  // Loop across atom types for each atom in the cluster, and build a corresponding loop.
  void build_int_maps_loop(int index, vector<int> atom_index, vector<string> ATOMTYPE, vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP);

  // Determine whether a particular cluster is excluded.
  bool is_excluded(vector<string> atom_names);

};


// These are no longer inherited from CLUSTER so that CLUSTER_LIST works correctly.
typedef QUADRUPLETS QUAD_FF;
typedef TRIPLETS TRIP_FF;



inline double CLUSTER::get_smaxim(string TYPE)	
// Decides whether outer cutoff should be set by 2-body value or cluster value. Returns the cutoff value.
{	
	double VAL = 0.0;
	
	for (int i=0; i< NPAIRS; i++)
	{
	  if(TYPE == ATOM_PAIRS[i]) 
	  {
		 VAL =  S_MAXIM[i];
		 break;
	  }
	}

	return VAL;	
}


inline double CLUSTER::get_sminim(string TYPE) 
// Decides whether outer cutoff should be set by 2-body value or 4-body value. Returns the cutoff value.
{
	double VAL = 0.0;
	
	
	for (int i=0; i<NPAIRS; i++) 
	{
	  if(TYPE == ATOM_PAIRS[i])
	  {
		 VAL =  S_MINIM[i];
		 break;
	  }
	}

	return VAL;	
}


#define _Cluster_h
#endif // ifndef _Cluster_h

