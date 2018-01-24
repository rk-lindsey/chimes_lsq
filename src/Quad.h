
#ifndef _Quad_h  // Protects against double-inclusion.

struct QUADRUPLETS
{
  int    QUADINDX;
  vector<string> ATOM_PAIRS; 	// Contains the six different atom pair string names
  vector<string> ATOM_NAMES ; // Names of each atom in the quad

  FCUT FORCE_CUTOFF;		// "CUBIC" "COSINE" or "SIGMOID" currently supported
	
  int    N_CFG_CONTRIB;		// How many configurations actually contribute to fitting this triplet??
	
  vector<double> MIN_FOUND;	// Testing two cases: 1. 
	
  vector<double> S_MAXIM_4B;	// A unique outer cutoff for 4B interactions... by default, is set to S_MAXIM
  vector<double> S_MINIM_4B;	// Similar for inner cutoff. This is useful when the 2-body 
  // is refit/extrapolated, thus has a s_min lower than the original fitted value
  // Values need to be specified for each contributing pair
  // [0] -> ij, [1] -> ik, [2] -> il, [3] -> jk, [4] -> jl, [5] -> kl 

	
  int N_TRUE_ALLOWED_POWERS;	// How many UNIQUE sets of powers do we have?
  int N_ALLOWED_POWERS;		// How many total sets of powers do we have?
	
  vector<vector<int> > ALLOWED_POWERS;	// This will keep a list of the allowed polynomial powers for each coefficient
  vector<vector<int> > UNIQUE_POWERS ;  // This is a list of unique polynomial powers for each coefficient
  vector<int> EQUIV_INDICIES;	// For each set of allowed powers, what is the index of the first equivalent set? For example, for the set (OO, OH, OH), (1,0,1) and (1,1,0) are is equivalent
  vector<int>	PARAM_INDICIES;	// For each of the set of allowed powers, what would be the index in the FF? for example, for a set of EQUIV_INDICIES {0,0,2,3}, PARAM_INDICIES would be {0, 0, 1, 2}
	
	
QUADRUPLETS(): ATOM_NAMES(4), ATOM_PAIRS(6), S_MAXIM_4B(6), S_MINIM_4B(6), MIN_FOUND(6) // Default constructor
	 {
		N_CFG_CONTRIB =  0;
		MIN_FOUND[0]   = -1;
		MIN_FOUND[1]   = -1;
		MIN_FOUND[2]   = -1;
		MIN_FOUND[3]   = -1;
		MIN_FOUND[4]   = -1;
		MIN_FOUND[5]   = -1;
		N_TRUE_ALLOWED_POWERS = 0;
	 }
		
  void init() ;  // Initialize values to defaults.
  void build(int cheby_4b_order) ; // The the ALLOWED_POWERS, etc. for an interaction.
  void store_permutations(vector<int> &unsorted_powers) ; // Store all the permutations.

private:
  // Recursively permute atom indices for each element type.
  void permute_atom_indices(int idx, vector<string> names, vector<int> &unsorted_powers, 
									 vector<int> perm, int unique_index, int equiv_index) ;
};


struct QUAD_FF : public QUADRUPLETS
{
	vector<double> 	PARAMS;
};

#define _Quad_h
#endif // ifndef _Quad_h

