
#ifndef _Quad_h  // Protects against double-inclusion.


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

	
  int N_TRUE_ALLOWED_POWERS;	// How many UNIQUE sets of powers do we have?
  int N_ALLOWED_POWERS;		// How many total sets of powers do we have?
	
  vector<vector<int> > ALLOWED_POWERS;	// This will keep a list of the allowed polynomial powers for each coefficient
  map<vector<int>,int> ALLOWED_POWERS_MAP ;  // Use a map for searching efficiency of allowed polynomial powers.
  vector<vector<int> > UNIQUE_POWERS ;  // This is a list of unique polynomial powers for each coefficient
  vector<int> EQUIV_INDICES;	// For each set of allowed powers, what is the index of the first equivalent set? For example, for the set (OO, OH, OH), (1,0,1) and (1,1,0) are is equivalent
  vector<int>	PARAM_INDICES;	// For each of the set of allowed powers, what would be the index in the FF? for example, for a set of EQUIV_INDICES {0,0,2,3}, PARAM_INDICES would be {0, 0, 1, 2}
		
  // Pure virtual (overridable) functions.
  void build(int cheby_4b_order) ; // The the ALLOWED_POWERS, etc. for an interaction.
  void store_permutations(vector<int> &unsorted_powers) ; // Store all the permutations.
  void print()  ;  // Print the quad powers and element types.

  // Print special parameters to the header file.
  void print_special(ofstream &header, string QUAD_MAP_REVERSE, string output_mode) ;

  // Print the params file header for a cluster
  void print_header(ofstream &header) ;

  CLUSTER() {} 
  virtual ~CLUSTER() {} 

CLUSTER(int natom, int npair): ATOM_NAMES(natom), ATOM_PAIRS(npair), S_MAXIM(npair), S_MINIM(npair), MIN_FOUND(npair) // Default constructor
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

class TRIPLETS : public CLUSTER
{
public:
	
  vector<vector<vector< int > > > POP_HIST; // Population histogram that s used to set 3B behavior in unsampled regions
	 
  vector<double> NBINS ;				// Number of bins to use for ij, ik, and jk distances when building the population histograms 
  vector<double> BINWS; 				// Binwidths to use for ij, ik, and jk distances when building the population histograms 
	
   TRIPLETS(): CLUSTER(3,3), NBINS(3), BINWS(3)
	 {
	 }

  // Virtual (overridable) functions.
};


class QUADRUPLETS : public CLUSTER
{
	
public:	
QUADRUPLETS():CLUSTER(4,6) { }


~QUADRUPLETS() { } 
		
};


class CLUSTER_LIST
{
public:
  // The number of CLUSTERS in the list.
  int NCLUSTERS ;        

  // A vector containing all clusters.
  // Need to store cluster pointers so that polymorphism with different cluster types works.
  vector<CLUSTER*> VEC ;

  // A map from pair types into the VEC index.
  map<string,int> MAP ;

  // A map from the VEC index into the pair types.
  map<int,string> MAP_REVERSE ;

  map<int,int> INT_MAP ;

  map<int,int> INT_MAP_REVERSE ;

  void build_maps(vector<struct PAIRS> &atom_pairs) ;

  void build_fast_maps(vector<string>& ATOM_CHEMS) ;

  void build_pairs(vector<string> ATOM_CHEMS, vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP) ;

  void link(vector<CLUSTER> &cluster) ;

  void link(vector<TRIPLETS> &cluster) ;

  void link(vector<QUADRUPLETS> &cluster) ;

private:
  void build_maps_loop(int index, vector<int> pair_index, vector<struct PAIRS> &atom_pairs) ;
  void build_fast_maps_loop(int index, vector<int> atom_index, vector<string>& ATOM_CHEMS) ;
  int make_id_int(vector<int>& index) ;
  void build_pairs_loop(int index, vector<int> atom_index, 
								vector<string> ATOM_CHEMS, vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP, int &count) ;
} ;


struct QUAD_FF : public QUADRUPLETS
{
	vector<double> 	PARAMS;
};

#define _Quad_h
#endif // ifndef _Quad_h

