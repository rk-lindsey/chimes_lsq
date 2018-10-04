#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes
#include<string>
#include<limits>	// Help with handling of over/underflow
#include <algorithm> // For specials vector-related functions (i.e. permute)


#ifdef USE_MPI
#include <mpi.h>
#endif 

using namespace std;

#include "functions.h"
#include "util.h"
#include "Cheby.h"

// define for extra output
//#define DEBUG_CLUSTER

void CLUSTER::build(int cheby_order)
// Build a set of interactions for a cluster.
// Figure out the allowed pair powers. Here are some rules and considerations:
//
// 1. Powers start from zero. So, if order is specified to be 2, polynomial powers range
//    from 0 to n-1, NOT 1 to n!
//
// 2. All four atoms must contribute to the interaction.  For instance, a polynomial with only
//    one non-zero polynomial power would not contribute.
//
// 3. Non-uniqueness upon atom permutation must be taken into consideration. 
//    Pairs permutations occur as a consequence of atom permutations.
//    We explicitly build symmetry-related interactions by permuting atom indices.
{



  vector<int> powers(NPAIRS);
					
  // Loops start at zero for a trick to speed up check for powers > 0
  
#ifdef DEBUG_CLUSTER
  if ( RANK == 0 ) 
  {
	 cout << "Building interactions for atoms: ";
	 for ( int j = 0; j < ATOM_NAMES.size(); j++ )
		cout << ATOM_NAMES[j] << " ";
	 cout << endl;
  }
#endif

  // For the present interaction type, generate a list of all possible allowed power sets

  if ( ! EXCLUDED ) 
	 build_loop(0, cheby_order, powers);
				
  N_TRUE_ALLOWED_POWERS = UNIQUE_POWERS.size();
  N_ALLOWED_POWERS      = PARAM_INDICES.size();

#ifdef DEBUG_CLUSTER
  if ( RANK == 0 ) 
  {
	 cout << "Number of unique powers = " << N_TRUE_ALLOWED_POWERS << endl;
	 cout << "Number of allowed powers = " << N_ALLOWED_POWERS << endl;

	 cout << "UNIQUE_POWERS " << endl;
	 for ( int j = 0; j < UNIQUE_POWERS.size(); j++ ) 
	 {
		for ( int k = 0; k < NPAIRS; k++ ) 
		  	cout << UNIQUE_POWERS[j][k] << " "; 
		cout << endl;
	 }

	 cout << "PARAM_INDICES " << endl;
	 for ( int j = 0; j < PARAM_INDICES.size(); j++ ) 
		cout << PARAM_INDICES[j] << endl;
  }
#endif

  if ( PARAM_INDICES.size() != ALLOWED_POWERS.size() ) 
  {
	 EXIT_MSG("Inconsistent size for PARAM_INDICES and ALLOWED_POWERS");
  }

  // Check that the power count is the same for interactions with the same PARAM_INDICES value.
  // If this is not true, the code will need to be generalized to support permutation sets
  // with differing POWER_COUNT values.
  //
  for ( int j = 0, k = 0; j < ALLOWED_POWERS.size(); ) 
  {
	 for ( ; k < ALLOWED_POWERS.size() && PARAM_INDICES[j] == PARAM_INDICES[k]; k++ ) {
		if ( POWER_COUNT[k] != POWER_COUNT[j] ) 
		{
		  cout << "j = " << j << " k = " << k << endl;
		  EXIT_MSG("Assumption of equal power counts was incorrect");
		}
		else 
		{
		 ;
#ifdef DEBUG_CLUSTER		  
		  cout << "Power count of " << j << " and " << k << "are equal to " << POWER_COUNT[j] << endl;
#endif
		}
	 }
	 j = k;
  }
}

void CLUSTER::build_loop(int indx, int cheby_order, vector<int> powers)
// Recursive function implements nested loops over pair powers.
{

  if ( indx < NPAIRS ) 
  {
	 for(int pair_pow=0; pair_pow < cheby_order; pair_pow++)
	 {
			powers[indx] = pair_pow;
			build_loop(indx+1, cheby_order, powers);
	 }
  } 
  else 
  {
	 vector<bool> represented(NATOMS,false);
	 int count = 0;
		
	 for ( int i = 0; i < NATOMS; i++ ) 
	 {
		for ( int j = i + 1; j < NATOMS; j++ ) 
		{
	 		if ( powers[count++] > 0 ) 
			{
		 		represented[i] = true;
		 		represented[j] = true;
			}
		}
	 }

	 
	 // All atoms must be represented with non-zero powers.
	 // Note:  the original 4-body code also had the requirement of 3 non-zero powers.
	 // An interaction of xij * xkl has 4 atoms represented, but only 2 non-zero powers.
	 // I think this should be included in the 4 body interaction.  Note that this term would
	 // not be included in a 3 body interaction. (LEF)
	 for ( int i = 0; i < NATOMS; i++ ) 
	 {
		if ( ! represented[i] ) 
		  return;
	 }
	 // Store all valid atom index permutations 
	 store_permutations(powers);
  }
}

void CLUSTER::set_atom_indices(vector<string> atomtype, vector<int> atomtype_idx)
// Set the ATOM_INDICES based on the given atomtype names and atomtype_idx indices.
{
	for ( int i = 0; i < NATOMS; i++ ) 
  	{
		int j;
		for ( j = 0; j < atomtype.size(); j++ )
			if ( atomtype[j] == ATOM_NAMES[i] ) 
				break;
	 	if ( j == atomtype.size() ) 
			EXIT_MSG("Could not set atom indices for " + ATOM_NAMES[i] );
		 else
			ATOM_INDICES[i] = atomtype_idx[j];
  	}
}





  
void CLUSTER::store_permutations(vector<int> &unsorted_powers) 
{
  	map<vector<int>,int>::iterator it;

  	// Avoid searching the ALLOWED_POWERS vector directly for efficiency. Search the ALLOWED_POWERS_MAP instead.
	
  	it = ALLOWED_POWERS_MAP.find(unsorted_powers);

  if ( it != ALLOWED_POWERS_MAP.end() ) 
  {
	 // This interaction is already handled.
#ifdef DEBUG_CLUSTER
	 if ( RANK == 0 ) 
	 {
		cout << "Found the following powers already stored:\n";
		for ( int i = 0; i < unsorted_powers.size(); i++ ) 
			cout << unsorted_powers[i] << " ";
		cout << endl;
	 }
#endif

	return;
  } 
  else 
  {
#ifdef DEBUG_CLUSTER
	 if ( RANK == 0 ) 
	 {
		cout << "Adding new unique set of powers" << endl;
		for ( int i = 0; i < unsorted_powers.size(); i++ ) 
		  cout << unsorted_powers[i] << " ";
		cout << endl;
	 } 
	 else 
	 {
		EXIT_MSG("Could not find the unsorted powers in the ALLOWED_POWERS");
	 }
#endif
	 // This is the first time that a combination of powers has been found.
	 UNIQUE_POWERS.push_back(unsorted_powers);	 
  }


  // Recursively permute each element and added to the ALLOWED_POWERS.
#ifdef DEBUG_CLUSTER
  if ( RANK == 0 ) 
  {
	 cout << "Permuting atom indices for ";
		for ( int i = 0; i < NATOMS; i++ ) 
	 		cout << ATOM_NAMES[i] << " ";
	 cout << endl;
  }
#endif

  	vector<int> perm(NATOMS);
  	for ( int i = 0; i < NATOMS; i++ ) 
		perm[i] = i;

  int equiv_index = ALLOWED_POWERS.size();
  permute_atom_indices(0, ATOM_NAMES, unsorted_powers, perm, UNIQUE_POWERS.size()-1, equiv_index );

#ifdef DEBUG_CLUSTER  
  cout << "Done permuting atom indices\n\n";
#endif

}








void CLUSTER::permute_atom_indices(int idx, vector<string> names, vector<int> &unsorted_powers, vector<int> perm, int unique_index, int equiv_index) 
{
  	if ( idx < names.size() ) 
	{
	 // Count up the number of identical atoms.
	 int count = 1;
	 for ( int i = idx; i < names.size() - 1; i++ )
	 {
			if ( names[i+1] == names[i] ) 
		 		++count;
			else
		  		break;
	 }
	 do	// Recursive call to generate permutations for identical atoms.
	{

		// DEBUG
#ifdef DEBUG_CLUSTER
		if ( RANK == 0 ) 
		  cout << "Number of identical atoms of element " << names[idx] << " = " << count << endl;
#endif
		permute_atom_indices(idx + count, names, unsorted_powers, perm, unique_index, equiv_index);
	 	} 
		while ( std::next_permutation(perm.data()+idx, perm.data()+idx+count) );
  	} 
	else 
	{
	 // Store the generated permutation.  Just love c++ syntax.... don't you ?
	 
	 vector<vector<int>> pow_mat(NATOMS, vector<int>(NATOMS));
	 vector<int>         perm_powers(NPAIRS);

	 // Generate the permutation transformation matrix.
	 int count = 0;
	 for ( int i = 0; i < NATOMS; i++ ) 
	 {
		pow_mat[i][i] = 0;
		for ( int j = i + 1; j < NATOMS; j++ ) 
		  	pow_mat[i][j] = unsorted_powers[count++];
	 }
	 /// Fill in other side of the diagonal
	 for ( int k = 0; k < NATOMS; k++ ) 
	 {
		for ( int l = 0; l < k; l++ ) 
		  	pow_mat[k][l] = pow_mat[l][k];
	 }

	 // Permute the powers.
	 count = 0;
	 for ( int i = 0; i < NATOMS; i++ ) 
	 {
		pow_mat[i][i] = 0;
		for ( int j = i + 1; j < NATOMS; j++ ) 
		  	perm_powers[count++] = pow_mat[ perm[i] ] [ perm[j] ];  
	 }

#ifdef DEBUG_CLUSTER
	 if ( RANK == 0 ) 
	 {
		cout << "Atom names: ";      for ( int j = 0; j < NATOMS; j++ ) cout << names[ perm[j] ] << " "; cout << endl;
		cout << "Permutation: ";     for ( int j = 0; j < NATOMS; j++ ) cout << perm[j]          << " "; cout << endl;
		cout << "Permuted powers: "; for ( int j = 0; j < NPAIRS; j++ ) cout << perm_powers[j]   << " "; cout << endl;
	 }
#endif
	 map<vector<int>,int>::iterator it;

	 it = ALLOWED_POWERS_MAP.find(perm_powers);

	 if ( it == ALLOWED_POWERS_MAP.end() ) 
	 {

#ifdef DEBUG_CLUSTER
		if ( RANK == 0 ) 
		{
		  cout << "Adding new allowed powers" << endl;
		  for ( int i = 0; i < perm_powers.size(); i++ ) 
			 cout << perm_powers[i] << " ";
		  cout << endl;
		}
#endif

		ALLOWED_POWERS.push_back(perm_powers);
		POWER_COUNT.push_back(1);
		if ( ALLOWED_POWERS.size() != POWER_COUNT.size() )
		 	EXIT_MSG("Inconsistent size for ALLOWED_POWERS and POWER_COUNT\n");

		// Store the index of perm_powers in the ALLOWED_POWERS vector in the ALLOWED_POWERS_MAP.
		ALLOWED_POWERS_MAP.insert( std::pair<vector<int>,int>(perm_powers, 
																				ALLOWED_POWERS.size()-1) );
		PARAM_INDICES.push_back(unique_index);
		EQUIV_INDICES.push_back(equiv_index);
#ifdef DEBUG_CLUSTER
		cout << "Power count = " << POWER_COUNT[POWER_COUNT.size() - 1] << endl;
#endif
	 }
	 else 
	 {
		int idx_itt = it->second;

		++POWER_COUNT[idx_itt];
#ifdef DEBUG_CLUSTER
		cout << "Found existing powers\n";
		for ( int l = 0; l < ALLOWED_POWERS[idx_itt].size(); l++ )
		{
		  cout << ALLOWED_POWERS[idx_itt][l] << " ";
		}
		cout << endl;
		cout << "Power count = " << POWER_COUNT[idx_itt] << endl << endl;
#endif
	 }
  }
}

void CLUSTER::print(bool md_mode) const
// Print parameters for the cluster.  If md_mode is true, print the potential parameters.
{
  if ( RANK == 0 ) 
  {
	 cout << "    Cluster index       : " <<  INDX << endl;
	 cout << "    Atoms in the cluster: ";
	 for ( int j = 0; j < NATOMS; j++ ) 
		cout << ATOM_NAMES[j] << " ";
	 cout << endl;
	 cout << "    Pairs in the cluster: ";
	 for(int j=0; j < NPAIRS; j++)
	 {
		cout  << ATOM_PAIRS[j] << " ";					
	 }
	 cout << endl;

	 if ( EXCLUDED ) 
	 {
		cout << "   Interaction is excluded by the user" << endl;
	 }
	 else 
	 {
		cout<< "    Number of unique sets of powers: " << N_TRUE_ALLOWED_POWERS << " (" << N_ALLOWED_POWERS << " total)..." << endl; 

		cout << "		     index  |  powers  |  equiv index  |  param index ";
		if ( md_mode ) 
		  cout << " | parameter ";
		else
		  cout << " | count ";
		cout << endl;

		cout << "		   ----------------------------------------------------" << endl;					
				
		for(int j=0; j<ALLOWED_POWERS.size(); j++)
		{

		  cout << "		      " << setw(6) << fixed << left << j << " ";
		
		  for(int k=0; k<NPAIRS; k++)
			 cout << " " << setw(2) << fixed << left << ALLOWED_POWERS[j][k] << " ";
						
		  cout << "       " << setw(8) << EQUIV_INDICES[j] << " ";
		  cout << "       " << setw(8) << PARAM_INDICES[j] << " "; 
		  
		  
		  if ( md_mode )
		  	cout << PARAMS[j];
		  else
		  	cout << "       " << setw(8) << POWER_COUNT[j] << " ";

		  cout << endl;
		
		}
	 }
	 cout << endl;
  }
}

void CLUSTER::print_special(ofstream &header, string MAP_REVERSE, string output_mode)
{
  if ( EXCLUDED ) 
	 return;

  if ( output_mode == "S_MINIM" ) 
  {
	 if( SPECIAL_S_MINIM ) 
		header << " " << MAP_REVERSE << " ";
	 for ( int i = 0; i < NPAIRS; i++ ) 
		header << ATOM_PAIRS[i] << " ";

	 header << fixed << setprecision(5);

	 for ( int i = 0; i < NPAIRS; i++ ) 
		header << S_MINIM[i] << " ";

	 header << endl;						
  } 
  else if ( output_mode == "S_MAXIM" ) 
  {
	 if ( SPECIAL_S_MAXIM )
		header << MAP_REVERSE << " ";
	 for ( int i = 0; i < NPAIRS; i++ ) 
		header << ATOM_PAIRS[i] << " ";

	 header << fixed << setprecision(5);

	 for ( int i = 0; i < NPAIRS; i++ ) 
		header << S_MAXIM[i] << " ";

	 header << endl;						
  }						
  else 
  {
	 cout << "Bad special parameter mode found\n";
	 exit(1);
  }
	 
}

void CLUSTER::print_header(ofstream &header)
// Print the header file for the force field definition
{
  header << " INDEX: " << INDX << " ATOMS: ";
  for ( int m = 0; m < NATOMS; m++ ) 
  {
	 header << ATOM_NAMES[m];
	 if ( m < NPAIRS - 1 ) 
		header << " ";
  }
  header << endl;
  header << " PAIRS: ";
  for(int m=0; m<NPAIRS; m++)
  {
	 header << ATOM_PAIRS[m];
	 if(m<NPAIRS-1)
		header << " ";
  }
  if ( EXCLUDED )
  {
	 header << " EXCLUDED:\n";
  } 
  else 
  {
	 header << "  " << N_TRUE_ALLOWED_POWERS << " parameters, " << N_ALLOWED_POWERS << " total parameters "<< endl;	
	 header << "     index  |  powers  |  equiv index  |  param index  " << endl;
	 header << "   ----------------------------------------------------" << endl;	

	 for(int j=0; j<ALLOWED_POWERS.size(); j++)
	 {
		header << "      " << setw(6) << fixed << left << j << " ";
		header << " ";
		for(int m=0; m<NPAIRS; m++)
		  header << setw(2) << fixed << left << ALLOWED_POWERS[j][m] << " ";
		header << "       " << setw(8) << EQUIV_INDICES[j] << " ";
		header << "       " << setw(8) << PARAM_INDICES[j] << endl; 
		
	 }
  }
  header << endl;
}


void CLUSTER::read_ff_params(ifstream &paramfile, const vector<string> &atomtype)
// Read the force field parameters for a cluster.
{
  string LINE;
  stringstream	STREAM_PARSER;
  string TEMP_STR;

  LINE = get_next_line(paramfile); // Blank line
  LINE = get_next_line(paramfile); // "TRIPLETYP PARAMS:"
  LINE = get_next_line(paramfile); // "INDEX: <indx> ATOMS: <a1> <a2> <a3>"
				 
  vector<string> tokens;
  int ntokens = parse_space(LINE, tokens);

  if ( ntokens != NATOMS + 3 ) 
	 EXIT_MSG("Error reading ff: wrong number of arguments: " + LINE);

  if ( tokens[0] != "INDEX:" ) 
	 EXIT_MSG("Error reading ff: INDEX: not found: " + LINE);

  INDX = stoi(tokens[1]);

  for ( int j = 0; j < NATOMS; j++ ) 
  {
	 ATOM_NAMES[j] = tokens[j+3];
  }

  LINE = get_next_line(paramfile);  // PAIRS: <p1> <p2> <p3> UNIQUE: <uniq> TOTAL: <TOTAL>
  
  ntokens = parse_space(LINE, tokens);

  if ( ntokens > NPAIRS ) {
	 if ( tokens[0] != "PAIRS:" ) 
		EXIT_MSG("Error reading ff: PAIRS: not found: " + LINE);

	 for ( int j = 0; j < NPAIRS; j++ )
		ATOM_PAIRS[j] = tokens[j+1];

	 if ( tokens.size() == NPAIRS + 2 && tokens[NPAIRS+1] == "EXCLUDED:" )
	 {
		EXCLUDED = true;
		N_TRUE_ALLOWED_POWERS = 0;
		N_ALLOWED_POWERS = 0;
		return;
	 }
  }
  else if ( ntokens != NPAIRS + 5 )
	 EXIT_MSG("Error reading ff: wrong number of arguments: " + LINE);

  if ( tokens[1+NPAIRS] != "UNIQUE:" )
	 EXIT_MSG("Error reading ff: UNIQUE: not found: " + LINE);

  N_TRUE_ALLOWED_POWERS = stoi(tokens[2+NPAIRS]);

  if ( tokens[3+NPAIRS] != "TOTAL:" ) 
	 EXIT_MSG("Error reading ff: TOTAL: not found: " + LINE);
  
  N_ALLOWED_POWERS = stoi(tokens[4+NPAIRS]);
				
  ALLOWED_POWERS.resize(N_ALLOWED_POWERS);
  for ( int j = 0; j  < N_ALLOWED_POWERS; j++ ) 
	 ALLOWED_POWERS[j].resize(NPAIRS);
				  
  EQUIV_INDICES.resize(N_ALLOWED_POWERS);
  PARAM_INDICES.resize(N_ALLOWED_POWERS);
  PARAMS        .resize(N_ALLOWED_POWERS);
				
  STREAM_PARSER.str("");
  STREAM_PARSER.clear();	

  LINE = get_next_line(paramfile); // header line
  LINE = get_next_line(paramfile);	// dashes line

  for(int j=0; j<N_ALLOWED_POWERS; j++)
  {
	 paramfile >> TEMP_STR;
	 for ( int k = 0; k < NPAIRS; k++ ) 
		paramfile >> ALLOWED_POWERS[j][k];

	 paramfile >> EQUIV_INDICES[j];
	 paramfile >> PARAM_INDICES[j];
	 paramfile >> PARAMS[j];

	 // if ( RANK == 0 ) 
	 // {
	 // 	cout << EQUIV_INDICES[j] << " " << PARAM_INDICES[j] << " " << PARAMS[j] << endl;
	 // }
	 paramfile.ignore();

  } 
  if (RANK==0)
	 cout << "	...Read " << NATOMS << "-body FF params..." << endl;;				
}


void CLUSTER::set_default_smaxim(const vector<PAIRS> & FF_2BODY)
// Decides whether outer cutoff should be set by 2-body value or cluster value. Returns the cutoff value.
{	
	for (int i=0; i< NPAIRS; i++)
	{

	  if ( ! SPECIAL_S_MAXIM )
	  {
		 int j;
		 for ( j = 0; j < FF_2BODY.size(); j++ ) 
		 {
			if(  FF_2BODY[j].PRPR_NM == ATOM_PAIRS[i] )
			{
			  S_MAXIM[i] = FF_2BODY[j].S_MAXIM;
			  break;
			}
		 }
		 if ( j == FF_2BODY.size() ) 
		 {
			EXIT_MSG("A default smaxim could not be set for " + ATOM_PAIRS[i]);
		 }
	  }
	}
}


void CLUSTER::set_default_sminim(const vector<PAIRS> & FF_2BODY)
// Decides whether outer cutoff should be set by 2-body value or cluster value. Returns the cutoff value.
{	
	for (int i=0; i< NPAIRS; i++)
	{
	  if ( ! SPECIAL_S_MINIM )
	  {
		 int j;
		 for ( j = 0; j < FF_2BODY.size(); j++ ) 
		 {
			if(  FF_2BODY[j].PRPR_NM == ATOM_PAIRS[i] )
			{
			  S_MINIM[i] = FF_2BODY[j].S_MINIM;
			  break;
			}
		 }
		 if ( j == FF_2BODY.size() ) 
		 {
			EXIT_MSG("A default sminim could not be set for " + ATOM_PAIRS[i]);
		 }
	  }
	}
}


int CLUSTER::match_atom_type_idx(string atm_typ)
// Return the atom type index of an atom in the cluster with the given name.
{
  int type = -1;

  for ( int j = 0; j < NATOMS; j++ ) 
  {
	 if ( ATOM_NAMES[j] == atm_typ ) 
	 {
		type = ATOM_INDICES[j];
		break;
	 }
  }
  return type;
}

string CLUSTER_LIST::allocate(int nclusters, int natoms, const vector<PAIRS> &FF_2BODY)
// Allocate the cluster list according to the number of clusters and the number
// of atoms.  Return a string to search for in the params file that describes
// the cluster list.
{
  NCLUSTERS = nclusters;
  int npairs = natoms * ( natoms - 1 ) / 2;
  VEC.resize(nclusters, CLUSTER(natoms,npairs) );	
			
  string pair_type = FF_2BODY[0].PAIRTYP;

  string tuplet = tuplet_name(natoms, false, false);
  for ( auto & c: tuplet ) c = toupper(c);

  string search = tuplet + " ";
  search.append(pair_type); // Syntax ok b/c all pairs have same FF type, and 2b and 3b are same pair type
  search.append(" PARAMS");	

  return search;
}

void CLUSTER_LIST::set_default_cutoffs(const vector<PAIRS>& FF_2BODY)
// Set default values for cutoffs.
{

  for (int i=0; i < NCLUSTERS; i++)
  {	
	 VEC[i].set_default_sminim(FF_2BODY);
	 VEC[i].set_default_smaxim(FF_2BODY);
  }	
}

void CLUSTER_LIST::read_maps(ifstream& paramfile, string line)
// Read the maps from the parameter file.
{
  vector<string> tokens;
  int nmaps = 0;

  if ( parse_space(line, tokens) < 2 ) 
	 EXIT_MSG("Not enough tokens while reading maps: " + line);
  else
	 nmaps = stoi(tokens[1]);

  string tuplets = tuplet_name(VEC[0].NATOMS, false, false);

  if (RANK==0)
	 cout << "	Reading  " << nmaps << " " << tuplets << " for mapping" << endl;
			
  for(int i=0; i< nmaps; i++)
  {
	 line = get_next_line(paramfile);
	 int index;
	 string cluster_pairs;
	 if ( parse_space(line,tokens) >= 2 ) 
	 {
		index = stoi(tokens[0]);
		cluster_pairs = tokens[1];

		if (( cluster_pairs.length() != 2 * VEC[0].NPAIRS)  && (cluster_pairs.length() != 4 * VEC[0].NPAIRS) ) // allow for 2-letter element names
		  EXIT_MSG("Incorrect cluster length: " + line);

		if (RANK==0)
		  cout << "	........Reading " << tuplets << ": " << cluster_pairs << " with mapped index: " << index << endl; 
				
		MAP.insert( make_pair(cluster_pairs, index) );
		MAP_REVERSE.insert( make_pair(index, cluster_pairs) );				
		
	 }
	 else
		EXIT_MSG("Not enough tokens while parsing: " + line);
  }

  if (RANK==0)
	 cout << "	...Read FF " << tuplets << " maps..." << endl;					
}


void CLUSTER_LIST::parse_fcut(string LINE)
// Parse the force cutoff parameters for a cluster.
{
  VEC[0].FORCE_CUTOFF.parse_input(LINE);
  VEC[0].FORCE_CUTOFF.BODIEDNESS = VEC[0].NATOMS;

  // Copy all class members.
  for(int i=1; i<VEC.size(); i++)
	 VEC[i].FORCE_CUTOFF = VEC[0].FORCE_CUTOFF;

}

void CLUSTER_LIST::read_ff_params(ifstream &PARAMFILE, const vector<string>& TMP_ATOMTYPE)
// Read the force field parameters for a cluster list.
{
  for(int i=0; i< NCLUSTERS; i++)
  {
	 VEC[i].read_ff_params(PARAMFILE, TMP_ATOMTYPE);
  }
}


void CLUSTER_LIST::process_cutoff_params(string input_type,vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP)
// Process the stored special lsq cutoff parameters (S_MAXIM, S_MINIM).
// Returns the largest outer cutoff or the smallest inner cutoff found.
{

  // Use the data_vec reference so the same code can change minimum and maximum values.
  using doublevec = vector<double>;

  vector<doublevec*> data_vec(NCLUSTERS);
  vector<bool *> special_flag(NCLUSTERS);

	vector<string> *input_ptr ;
	string LINE ;

  if ( VEC.size()<1)
  {
	 cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "ERROR: Special inner cutoffs specified for many-body chebyshev interactions,"  << COUT_STYLE.ENDSTYLE << endl;
	 cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "ERROR: expected number of pair triplets is zero. Did you forget to set "  << COUT_STYLE.ENDSTYLE << endl;
	 cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "ERROR: the \'ATOM PAIR ****PLETS:\' line in the parameter file?"  << COUT_STYLE.ENDSTYLE << endl;
	 exit_run(0);
  }

  int npairs = VEC[0].NPAIRS;
  int natoms = VEC[0].NATOMS;

  if ( input_type == "S_MAXIM" )
  {
		input_ptr = &S_MAXIM_INPUT ;
	 for ( int i = 0; i < NCLUSTERS; i++ ) 
	 {
		data_vec[i] = &(VEC[i].S_MAXIM);
		special_flag[i] = &(VEC[i].SPECIAL_S_MAXIM);
		
		
	 }
  }
  else if ( input_type == "S_MINIM" ) 
  {
		input_ptr = &S_MINIM_INPUT ;
	 for ( int i = 0; i < NCLUSTERS; i++ )
	 {
		data_vec[i] = &(VEC[i].S_MINIM);
		special_flag[i] = &(VEC[i].SPECIAL_S_MINIM);
	 }
  }
  else
  {
	 EXIT_MSG("UNKNOWN INPUT type: " + input_type);
  }

	if ( input_ptr->size() > 0 ) 
	{
		vector<string> tokens;
		LINE = (*input_ptr)[0] ;
		int count = 0 ;

		if ( parse_space(LINE, tokens) < 5 ) 
			EXIT_MSG("Not enough tokens: " + LINE);

		string TEMP_STR = tokens[3];
			
		if(TEMP_STR == "ALL" )
		{				
			for(int i=0; i < NCLUSTERS; i++)
			{
				for ( int j = 0; j < npairs; j++ )
				{
					(*data_vec[i])[j] = stod(tokens[4]);
					*special_flag[i] = true;
				}
			}
			#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	Note: Setting all " << natoms << "-body " + input_type + " values to " <<  (*data_vec[0])[0] << endl;
			#endif				
		}
		else if(TEMP_STR == "SPECIFIC" )
		{
			int n_spec = stoi(tokens[4]);
				
			#if VERBOSITY == 1
			if ( RANK == 0 ) 
				cout << "	Note: Setting specific " << natoms << "-body " + input_type + " values: " << endl;
			#endif	

			string cluster_name;
			vector<string> pair_names(npairs);
			vector<double> cutoff_vals(npairs);

			for(int i=0; i< n_spec; i++)
			{
				++count ;
				LINE = (*input_ptr)[count] ;

				if ( parse_space(LINE, tokens) != 2 * npairs + 1 )
					EXIT_MSG("Wrong number of arguments: " + LINE);

				cluster_name = tokens[0];
				int cluster_index = 0;
				if ( MAP.find(cluster_name) != MAP.end() ) 
					cluster_index = MAP[cluster_name];
				else
				{
					cout << "Did not find a map for cluster named: " << cluster_name << endl;
					cout << "Possible choices are: " << endl;
					for(int m=0; i<MAP_REVERSE.size(); m++)
						cout << MAP_REVERSE[m] << endl;
		  
		  
					EXIT_MSG("Did not find a map for " + cluster_name);
				}

				for ( int i = 0; i < npairs; i++ ) 
				{
					pair_names[i] = tokens[i+1];
					cutoff_vals[i] = stod(tokens[i + npairs + 1]);

					try
					{
						pair_names[i] = ATOM_PAIRS[ PAIR_MAP[ pair_names[i] ] ].PRPR_NM;
					}
					catch(...)
					{
						cout << "ERROR: Unknown pair: " + LINE;
					}
					for ( int j = 0; j < npairs; j++ ) 
					{
						// Set the cutoff to be the same for all matching pairs.
						// This enforces consistency between pair cutoffs with the same name.
						if ( pair_names[i] == VEC[cluster_index].ATOM_PAIRS[j] )
							(*data_vec[cluster_index])[j] = cutoff_vals[i];
						*special_flag[cluster_index] = true;
					}
				}
			}
			if ( RANK == 0 ) 
			{
				for ( int i = 0; i < NCLUSTERS; i++ ) 
				{
					cout << "	";
					for ( int j = 0; j < natoms; j++ ) 
					{
						cout << VEC[i].ATOM_NAMES[j] << " ";
					}
					for ( int j = 0; j < npairs; j++ ) 
					{
						cout << VEC[i].ATOM_PAIRS[j] << " " << (*data_vec[i])[j] << " ";
					}
					cout << endl;
				}
				cout << endl;
			}
		}
	}
  double val = 0;
  if ( input_type == "S_MAXIM" ) 
  {
	 val = -1.0;
	 for ( int i = 0; i < NCLUSTERS; i++ ) 
	 {
		for ( int j = 0; j < npairs; j++ ) 
		{
		  if ( VEC[i].S_MAXIM[j] > val ) 
			 val = VEC[i].S_MAXIM[j];
		}
	 }
	 if ( RANK == 0 ) 
		cout << "     Max " << natoms << "-body cutoff = " << val << endl << endl;
	 MAX_CUTOFF = val ;
  }
  else if ( input_type == "S_MINIM" ) 
  {
  	 val = 1.0e10;
	 for ( int i = 0; i < NCLUSTERS; i++ ) 
	 {
		for ( int j = 0; j < npairs; j++ ) 
		{
		  if ( VEC[i].S_MINIM[j] < val ) 
			 val = VEC[i].S_MINIM[j];
		}
	 }
	 if ( RANK == 0 ) 
		cout << "     Min " << natoms << "-body cutoff = " << val << endl;
  }
}


void CLUSTER_LIST::read_cutoff_params(istream &input, string LINE, string input_type)
// Read the lsq cutoff parameters (S_MAXIM, S_MINIM) and store for later processing by process_cutoff_params.
{
	vector<string> *line_ptr ;

  if ( input_type == "S_MAXIM" )
  {
		line_ptr = &S_MAXIM_INPUT ;
  }
  else if ( input_type == "S_MINIM" ) 
  {
		line_ptr = &S_MINIM_INPUT ;
  }
  else
  {
	 EXIT_MSG("UNKNOWN INPUT type: " + input_type + " " + LINE);
  }
	line_ptr->push_back(LINE) ;

  vector<string> tokens;

  if ( parse_space(LINE, tokens) < 5 ) 
	 EXIT_MSG("Not enough tokens: " + LINE);

  string TEMP_STR = tokens[3];
			
  if ( TEMP_STR == "SPECIFIC" )
  {
	 int n_spec = stoi(tokens[4]);
				
	 for(int i=0; i< n_spec; i++)
	 {
		LINE = get_next_line(input);
		line_ptr->push_back(LINE) ;
	 }
	}
}

void CLUSTER_LIST::build_int_maps(vector<string> ATOMTYPE, vector<int> ATOMTYPE_IDX, vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP)
// Build the fast integer maps for the cluster list.  Used by the MD code.
{
  // Find the dimension of the map vector.
  int NATOMS = VEC[0].NATOMS;

  int dim = 1;
  for ( int i = 0; i < NATOMS; i++ ) 
	 dim *= MAX_ATOM_TYPES;
  
  INT_MAP.resize(dim, -1);
  PAIR_INDICES.resize(dim);

  INT_MAP_REVERSE.resize(NCLUSTERS, -1);

  string tuplet = tuplet_name(NATOMS, false, true);
  if(RANK==0)
	 cout << endl << "	" << tuplet << " maps:" << endl;

  vector<int> atom_index(NATOMS);

  for ( int i = 0; i < VEC.size(); i++ ) 
	 VEC[i].set_atom_indices(ATOMTYPE, ATOMTYPE_IDX);
	 
  build_int_maps_loop(0, atom_index, ATOMTYPE, ATOM_PAIRS, PAIR_MAP);

  set_default_cutoffs(ATOM_PAIRS);
  for ( int i = 0; i < VEC.size(); i++ ) 
  {
	 VEC[i].set_cheby_vals(ATOM_PAIRS);
  }
  
}

void CLUSTER_LIST::build_int_maps_loop(int index, vector<int> atom_index, vector<string> ATOMTYPE, vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP)
// Loop across atom types for each atom in the cluster, and build a corresponding loop.
{

  int NATMTYP = ATOMTYPE.size();
  int NATOMS = VEC[0].NATOMS;
  int NPAIRS = VEC[0].NPAIRS;

  if ( index < NATOMS ) 
  {
	 for ( int j = 0; j < NATMTYP; j++ ) 
	 {
		atom_index[index] = j;
		build_int_maps_loop(index+1, atom_index, ATOMTYPE, ATOM_PAIRS, PAIR_MAP);
	 }
  } 
  else 
  {
	 vector<string> pairs(NPAIRS);

	 int count = 0;
	 // Build up permuted pair names.
	 for ( int j = 0; j < NATOMS; j++ ) 
	 {
		for ( int k = j + 1; k < NATOMS; k++ ) 
		{
		  pairs[count] = ATOMTYPE[atom_index[j]] + ATOMTYPE[atom_index[k]];
		  pairs[count] = ATOM_PAIRS[ PAIR_MAP[ pairs[count] ] ].PRPR_NM;
		  count++;
		}
	 }

	 // Sort each entry in descending order.
	 //sort(pairs.begin(), pairs.end());
	 //reverse(pairs.begin(), pairs.end() );

	 // Get the integer ID for this triplet.
	 int idx1 = make_id_int(atom_index);

	 // Get the string for the triplet.
	 string int_map_str = "";
	 for ( int j = 0; j < NPAIRS; j++ ) 
	 {
		int_map_str += pairs[j];
	 }

#ifdef DEBUG_CLUSTER
	 cout << "MAP[ " << int_map_str << " = " << MAP[int_map_str] << endl;
#endif
	 
	 if ( MAP.find(int_map_str) != MAP.end() )
	 {
		INT_MAP[idx1] = MAP[int_map_str];

		// Store the mapping of pair indices for later use in force evaluation.
		PAIR_INDICES[idx1].resize(NPAIRS);


		if(MAP[int_map_str]>=0)
		VEC[INT_MAP[idx1]].map_indices_int(atom_index, PAIR_INDICES[idx1]);
	 }
	 else
		EXIT_MSG("Map for " + int_map_str + " was not found\n");

	 // MAP[int_map_str] is -1 for excluded interaction.
	 if ( INT_MAP[idx1] >= 0 ) 
		INT_MAP_REVERSE[MAP[int_map_str]] = idx1;

	 if(RANK == 0)
	 {
		string tuplet = tuplet_name(NATOMS, false, true);
		cout << "		";
		cout<< "Atom type idxs: ";
		for ( int j = 0; j < NATOMS; j++ ) 
		{
		  cout<< fixed << setw(2) << right << atom_index[j];
		}
		cout<< " " << tuplet << " name: "           << setw(12) << right << int_map_str;
		if(MAP[int_map_str]>=0)
			cout<< " Explicit index: " << setw(4) << right << idx1;
		else
			cout<< " Explicit index: " << "EXCLUDED";
		if(MAP[int_map_str]>=0)
			cout<< " Unique index: "   << setw(4) << right << INT_MAP[idx1] << endl;
		else
			cout<< " Unique index: "   << "EXCLUDED" << endl;
	 }
  }
}

int CLUSTER_LIST::build_all(int cheby_order, vector<PAIRS> & ATOM_PAIRS, map<string,int> &PAIR_MAP, vector<string> atom_types, vector<int> atomtype_idx)
// Build all triplets and associated maps for the cluster list.
{
  build_pairs(ATOM_PAIRS, PAIR_MAP, atom_types, atomtype_idx);
			  
  for ( int i = 0; i < VEC.size(); i++ ) 
	 VEC[i].build(cheby_order);

  set_default_cutoffs(ATOM_PAIRS);
  process_cutoff_params("S_MAXIM", ATOM_PAIRS, PAIR_MAP) ;
  process_cutoff_params("S_MINIM", ATOM_PAIRS, PAIR_MAP) ;

  return( VEC.size() );
}


int CLUSTER_LIST::make_id_int(vector<int>& index)
// Returns a unique ID number for a cluster of atoms.
{
  int sum = 0;
  int weight = 1;
  for ( int i = 0; i < index.size(); i++ ) {
	 sum += weight * (index[i]+1);
	 weight *= MAX_ATOM_TYPES;
  }
  return(sum);
}

void CLUSTER_LIST::build_cheby_vals(vector<PAIRS> & ATOM_PAIRS)
{
  for ( int i = 0; i < VEC.size(); i++ ) 
	 VEC[i].set_cheby_vals(ATOM_PAIRS);
}

void CLUSTER_LIST::build_pairs(vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP, vector<string> atom_types,
										 vector<int> atom_typeidx)
// Build the pair variables for all of the clusters
{
			
  // First, extract the atom types
  vector<string>ATOM_CHEMS;

  NCLUSTERS = VEC.size();
  for(int p=0; p<ATOM_PAIRS.size(); p++)
  {
	 string TMP_CHEM = ATOM_PAIRS[p].ATM1TYP;
	 bool   IN_LIST  = false;
	 
	 for(int a=0; a<ATOM_CHEMS.size(); a++)
	 {
		if(ATOM_CHEMS[a] == TMP_CHEM)
		  IN_LIST = true;
	 }
	 
	 if(!IN_LIST)
		ATOM_CHEMS.push_back(TMP_CHEM);
  }
				
  int natoms = VEC[0].NATOMS;
  int npairs = VEC[0].NPAIRS;

  //sort( ATOM_CHEMS.begin(), ATOM_CHEMS.end() );

  vector<int> atom_index(natoms);
  int count = 0;
  
  
  int dim = 1;
  for ( int i = 0; i < natoms; i++ ) 
	 dim *= MAX_ATOM_TYPES;
  
  INT_MAP.resize(dim, -1);
  INT_MAP_REVERSE.resize(NCLUSTERS, -1);
  PAIR_INDICES.resize(dim, vector<int>(npairs,-1));

  build_pairs_loop(0, atom_index, ATOM_CHEMS, ATOM_PAIRS, PAIR_MAP, count,
						 atom_types, atom_typeidx);

}



void CLUSTER_LIST::build_pairs_loop(int index, vector<int> atom_index, 
												vector<string> ATOM_CHEMS, vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP, int &count,
												vector<string> atom_types, vector<int> atom_typeidx)
// Loop over atom types in building the CLUSTER ATOM_PAIRS list.
{
  int NATMTYP = ATOM_CHEMS.size();
  int natoms = VEC[0].NATOMS;
  int npairs = VEC[0].NPAIRS;

  if ( index < atom_index.size() ) {
	 // Keep on looping.
	 for ( int i = 0; i < NATMTYP; i++ ) 
	 {
		atom_index[index] = i;
		build_pairs_loop(index+1, atom_index, ATOM_CHEMS, ATOM_PAIRS, PAIR_MAP, count,
							  atom_types, atom_typeidx);
	 }
  }
  else 
  {
	 // Perform a calculation based on the index values.
	 VEC[count].INDX = count;	// Index for current triplet type
								
	 // Save the names of each atom type in the pair.
	 int count2 = 0;
	 vector<string> pair_names(npairs);
	 vector<string> pair_names_unsorted(npairs);
	 vector<string> atom_names(natoms);
	 //vector<int> atom_index_rev(natoms);

	 // Sort the atom names, not the pair names.  All interactions
	 // must be based on atom names.
	 for ( int i = 0; i < natoms; i++ ) {
		atom_names[i] = ATOM_CHEMS[ atom_index[i] ];
	 }
	 sort(atom_names.begin(), atom_names.end() );

	 // Now get the reverse index:  the ith sorted atom name is the jth ATOM type.
	 // for ( int i = 0; i < natoms; i++ ) {
	 // 	atom_index_rev[i] = -1; 
	 // 	for ( int j = 0; j < NATMTYP; j++ ) {
	 // 	  if ( ATOM_CHEMS[j] == atom_names[i] ) {
	 // 		 atom_index_rev[i] = j;
	 // 		 break;
	 // 	  }
	 // 	}
	 // }
	 for ( int i = 0; i < natoms; i++ ) {
		for ( int j = i + 1; j < natoms; j++ ) {
		  pair_names[count2] = atom_names[i] + atom_names[j];
		  pair_names_unsorted[count2] = ATOM_CHEMS[ atom_index[i] ] + ATOM_CHEMS[ atom_index[j] ];
		  // Get the proper (ordered) name for each pair type.
		  pair_names[count2] = ATOM_PAIRS[ PAIR_MAP[ pair_names[count2] ] ].PRPR_NM;
		  pair_names_unsorted[count2] = ATOM_PAIRS[ PAIR_MAP[ pair_names_unsorted[count2] ] ].PRPR_NM;
		  count2++;
		}
	 }

	 // See if this is a new pair interaction.
	 bool new_pair = true;
	 int map_index = -1;
	 for ( int i = 0; i < count; i++ ) {
		int k;
		for ( k = 0; k < npairs; k++ ) {
		  if ( VEC[i].ATOM_PAIRS[k] != pair_names[k] ) {
			 break;
		  }
		}
		if ( k == npairs ) {
		  // Found a match.
		  new_pair = false;
		  map_index = i;
		  break;
		}
	 }

	 bool excluded = is_excluded(atom_names);
		
	 if ( new_pair ) 
	 {
		// This is a new set of pairs
		if(RANK==0)
		{							
		  cout << endl << "	Made the following " << tuplet_name(natoms, true, false) << ": " << count << " ";
		  int npairs = VEC[0].NPAIRS;
		  for(int m=0; m< npairs; m++) 
			 cout << VEC[count].ATOM_PAIRS[m] << " ";
		  cout << endl << endl;
		}

		if ( count >= VEC.size() ) {
		  VEC.resize(VEC.size()+1);
		  if ( RANK == 0 ) {
			 cout << "Resizing CLUSTER_LIST VEC to " << VEC.size() << endl;
		  }
		}
		count2 = 0;

#ifdef DEBUG_CLUSTER		
		if (RANK == 0 )
		  cout << "ADDING ATOM_PAIRS = ";
#endif
		  
		for ( int i = 0; i < natoms; i++ ) {
		  VEC[count].ATOM_NAMES[i] = atom_names[i];
		  int j;
		  for ( j = 0; j < atom_types.size(); j++ ) 
		  {
			 if ( VEC[count].ATOM_NAMES[i] == atom_types[j] )
				break;
		  }
		  if ( j == atom_types.size() )
			 EXIT_MSG("Could not find an element type for " + VEC[count].ATOM_NAMES[i]);
		  else
		  {
			 VEC[count].ATOM_INDICES[i] = atom_typeidx[j];
#ifdef DEBUG_CLUSTER
			 cout << "The index of cluster atom " << i << "is" << VEC[count].ATOM_INDICES[i] << endl;
#endif
		  }
			 
		  for ( int j = i + 1; j < natoms; j++ ) {
			 VEC[count].ATOM_PAIRS[count2] = pair_names[count2];
#ifdef DEBUG_CLUSTER			 
			 if ( RANK == 0 ) cout << VEC[count].ATOM_PAIRS[count2] << " ";
#endif
			 count2++;
		  }
		}
#ifdef DEBUG_CLUSTER
		if ( RANK == 0 ) cout << endl;
#endif
		
		if ( ! excluded )
		  map_index = count;
		else
		  VEC[count].EXCLUDED = true;

		count++;	
	 }

	 // Set up a mapping between pair name and cluster index, and vice versa.
	 string all_pairs = "";
	 for ( int i = 0; i < pair_names_unsorted.size(); i++ )
		all_pairs += pair_names_unsorted[i];

	 if ( excluded )
		map_index = -1;
	 
	 MAP.insert(make_pair(all_pairs,map_index));

#ifdef DEBUG_CLUSTER
	 if ( RANK == 0 ) cout << "MAP[ " << all_pairs << "] = " << map_index << endl;
#endif
	 
	 MAP_REVERSE.insert(make_pair(map_index,all_pairs));

	 // Construct a unique integer ID for the atom set and store it in the INT_MAP.
	 //vector<int> atom_index_2 = atom_index;
	 //sort( atom_index_2.begin(), atom_index_2.end() );
	 //reverse( atom_index_2.begin(), atom_index_2.end() );

	 // Make a unique ID for this ordering of atoms.
	 int atom_id_int = make_id_int(atom_index);
	 //int atom_id_int = make_id_int(atom_index_2);
	 INT_MAP[atom_id_int] = map_index;

	 // Create a mapping for pair properties like Cheby powers.
	 PAIR_INDICES[atom_id_int].resize(npairs);

	 // Negative map index for excluded interaction.
	 if ( map_index >= 0 ) {
		VEC[map_index].map_indices_int(atom_index, PAIR_INDICES[atom_id_int]);
		INT_MAP_REVERSE[map_index] = atom_id_int;
	 }

	 if(RANK == 0)
	 {
		cout << "		";
		cout<< "Atom type idxs: ";

		for ( int i = 0; i < natoms; i++ ) {
		  cout<< fixed << setw(2) << right << atom_index[i];
		}
		
		string tuplet = tuplet_name(natoms, false, true);
		cout<< "  " << tuplet << " name: "           << setw(12) << right << all_pairs;
		tuplet = tuplet_name(natoms, false, false);
		cout<< " Unique " << tuplet <<  " ID " << atom_id_int << " index: "   << setw(4)  << right << INT_MAP[atom_id_int] << endl;
	 }
  } // End of calculation.
}

bool CLUSTER_LIST::is_excluded(vector<string> atom_names)
// Returns true if the cluster list is excluded from the interaction.
{
  // Pass by value of atom_names makes the sorting safe.
  sort( atom_names.begin(), atom_names.end() );

  for ( int j = 0; j < EXCLUDE.size(); j++ ) 
  {
	 int k;
	 for ( k = 0; k < VEC[0].NATOMS; k++ )
	 {
		if ( EXCLUDE[j][k] != atom_names[k] ) 
		  break;
	 }
	 // See if we found a match from the exclude list.
	 if ( k == VEC[0].NATOMS ) 
		return(true);
  }

  return(false);
}


void CLUSTER_LIST::read_exclude(istream &input, string line)
// Read the excluded interactions from the input stream.
{
  int nexclude;

  vector<string> tokens;

  parse_space(line, tokens);
  if ( tokens.size() < 4 ) 
	 EXIT_MSG("Wrong number of parameters in EXCLUDE command: " + line);

  nexclude = stoi(tokens[3]);

  for(int i=0; i< nexclude; i++)
  {
	 line = get_next_line(input);
	 vector<string> elements;
	 
	 parse_space(line, elements);
	 if ( elements.size() != VEC[0].NATOMS )
		EXIT_MSG("Wrong number of atoms in exclude command: " + line);

	 sort(elements.begin(), elements.end() );
	 EXCLUDE.push_back(elements);
  }
			
}

void CLUSTER_LIST::print(bool md_mode)
{
  if ( RANK == 0 ) 
  {
	 string tuplet = tuplet_name(VEC[0].NATOMS,true,false);

	 cout << endl << "	The following unique " << tuplet << " of atoms have been identified:" << endl;

	 if ( EXCLUDE.size() > 0 )
		cout << "	Note: The following types have been removed, if present: " << endl;
					
	 for(int i=0;i<EXCLUDE.size(); i++)
	 {
		for ( int j = 0; j < EXCLUDE[i].size(); j++ )
		{
		  cout << "		" << EXCLUDE[i][j];
		}
		cout << endl;
	 }
	 cout << endl;

	 for(int i=0;i < VEC.size(); i++)
	 {
		VEC[i].print(md_mode);
	 }
  }
}

void CLUSTER_LIST::print_fcut()
// Print the cutoff function parameters controlling the CLUSTER_LIST.
{

  if ( NCLUSTERS <= 0 ) 
	 return;

  cout << "	Using the following fcut style for " << VEC[0].NATOMS << "B Chebyshev interactions: ";

  VEC[0].FORCE_CUTOFF.print_params();
}

void CLUSTER_LIST::print_min_distances()
{
  string tuplet = tuplet_name(VEC[0].NATOMS, false, false);

  if ( RANK == 0 ) cout << "	Minimum distances between atoms " << tuplet << " pairs: (Angstr.)" << endl;
			
  for (int k=0; k< VEC.size(); k++) 
  {
	 if ( VEC[k].EXCLUDED ) 
		continue;

	 vector<double> min(VEC[k].NPAIRS,2.0e10);
			
#ifdef USE_MPI
			
	 for (int m=0; m< VEC[k].NPAIRS; m++)
	 {
		MPI_Reduce(&(VEC[k].MIN_FOUND[m]), &(min[m]), 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	 }
#else	
	 for(int m=0; m< VEC[k].NPAIRS; m++)
		min[m] = VEC[k].MIN_FOUND[m];
#endif
					
	 if ( RANK == 0 )
	 {
		cout << "		" << k << "	";	
		  for(int m=0; m < VEC[k].NPAIRS; m++)
			 cout << VEC[k].ATOM_PAIRS[m] << " ";
		  for(int m=0; m < VEC[k].NPAIRS; m++)
		  {
			 if ( min[m] < 1.0e9 ) 
				cout << min[m] << " ";
			 else
				cout << "NOT FOUND  ";
		  }
		  cout << endl;
		}
  }
			
  if ( RANK == 0 ) cout << "	Total number of configurations contributing to each "
								<< tuplet << " type:" << endl;
		
  for (int k=0; k<VEC.size(); k++)
  {
	 int sum = 0;
				
#ifdef USE_MPI
	 MPI_Reduce(&(VEC[k].N_CFG_CONTRIB), &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	 sum = VEC[k].N_CFG_CONTRIB;
#endif	
								
	 if ( RANK == 0 ) 
	 {
		cout << "		" << k << "	";	
		for(int m=0; m< VEC[0].NPAIRS; m++)
		  cout << VEC[k].ATOM_PAIRS[m] << " ";
		if ( ! VEC[k].EXCLUDED )
		  cout << sum << endl;
		else
		  cout << "excluded\n";
	 }	
  }
}


void CLUSTER_LIST::print_header(ofstream &header, int natoms, int cheby_order)
// print the force field file header for the cluster list.
// The number of atoms in the cluster is specified in case the
// cluster is not initialized (used).
{
  string tuplet = tuplet_name(natoms, true, false);
  for ( auto & c: tuplet ) c = toupper(c);

  header << endl << "ATOM PAIR " << tuplet << ": ";
  if( cheby_order == 0 || NCLUSTERS == 0 )
	 header << 0 << endl;
  else
  {
	 header << VEC.size() << endl << endl;
	 for(int i=0;i < VEC.size(); i++)
	 {
		VEC[i].print_header(header);
	 }	 
  }
}

void CLUSTER_LIST::print_special(ofstream &header)
// Print out special 3 and 4 body force parameters.
{

  if ( NCLUSTERS == 0 ) 
	 return;
	 

  int FOUND_SPECIAL = 0;
  int NATOMS = VEC[0].NATOMS;

  for(int i=0; i<VEC.size(); i++)
  {
	 if( VEC[i].SPECIAL_S_MINIM && !VEC[i].EXCLUDED ) 
		FOUND_SPECIAL++;
  }
	
  if(FOUND_SPECIAL>0)
  {
	 header << endl << "SPECIAL " << NATOMS << "B S_MINIM: SPECIFIC " << FOUND_SPECIAL << endl;
		
	 for(int i=0; i<VEC.size(); i++)
		if(VEC[i].SPECIAL_S_MINIM)
		  VEC[i].print_special(header, MAP_REVERSE[i], "S_MINIM");
  }
		
  FOUND_SPECIAL = 0;
	
  for(int i=0; i<VEC.size(); i++)
	 if(VEC[i].SPECIAL_S_MAXIM  && !VEC[i].EXCLUDED) 
		FOUND_SPECIAL++;
	
  if(FOUND_SPECIAL>0)
  {
	 header << endl << "SPECIAL " << NATOMS << "B S_MAXIM: SPECIFIC " << FOUND_SPECIAL << endl;
		
	 for(int i=0; i<VEC.size(); i++)
		if(VEC[i].SPECIAL_S_MAXIM) 
		  VEC[i].print_special(header, MAP_REVERSE[i], "S_MAXIM");
  }	
}

string CLUSTER_LIST::tuplet_name(int natom, bool plural, bool caps)
{
  string tuplet;

  switch(natom) 
  {
  case 1:
	 tuplet = "singlet";
	 break;
  case 2:
	 tuplet = "pair";
	 break;
  case 3:
	 tuplet = "triplet";
	 break;
  case 4:
	 tuplet = "quadruplet";
	 break;
  case 5:
	 tuplet = "quintuplet";
	 break;
  case 6:
	 tuplet = "sextuplet";
	 break;
  case 7:
	 tuplet = "septuplet";
	 break;
  case 8:
	 tuplet = "octuplet";
	 break;
  case 9:
	 tuplet = "nontuplet";
	 break;
  default:
	 tuplet = to_string(natom) + "-tuplet";
	 break;
  };

  if ( plural ) 
	 tuplet += "s";

  if ( caps )
	 tuplet[0] = toupper(tuplet[0]);

  return tuplet;
}

bool CLUSTER::init_histogram(vector<PAIRS> &pairs, map<string,int>& pair_map)
// Sets up the histogram for CLUSTER.  Returns true on success,
// false otherwise.
{
		
  double tmp_max, tmp_min;
		
  for ( int j = 0; j < NPAIRS; j++ ) 
  {
	 NBINS[j] = pairs[ pair_map[ ATOM_PAIRS[0]] ].NBINS[j];
	 if(NBINS[j] == 0 ) 
	 {
#if VERBOSITY == 1
		if ( RANK == 0 ) cout << "Found at least one 3b-population histogram nbins = 0. Will not do ANY histogramming. " << endl << endl << endl;
#endif
		
		return false;
	 }
	 tmp_max = S_MAXIM[j];
	 tmp_min = S_MINIM[j];
	 if(tmp_min == -1)
		tmp_min = pairs[ pair_map[ ATOM_PAIRS[j]] ].S_MINIM;
	 if(tmp_max == -1)
		tmp_max = pairs[ pair_map[ ATOM_PAIRS[j]] ].S_MAXIM;
	 BINWS[0] = (tmp_max - tmp_min)/NBINS[0];
  }
		
		
  // Note:  pop_hist is now implemented as a map to allow arbitrary dimensional histograms.
  // No allocation of pop_hist is necessary.  

  if ( RANK == 0 ) 
  {
	 cout << "	" << INDX << " " << fixed << setprecision(1);
	 for ( int j = 0; j < NPAIRS; j++ ) 
	 {
		cout << NBINS[j] << " ";
	 }
  }

#if VERBOSITY == 1
  if ( RANK == 0 ) 
  {
	 cout << "...Initial histogram setup complete!" << endl << endl;
}
#endif
  return(true);
}

void CLUSTER::increment_histogram(vector<int>& index)
// Increment the histogram with the given index vector.
{

  if ( POP_HIST.find(index) != POP_HIST.end() ) 
	 POP_HIST[index]++;
  else
	 POP_HIST.insert(make_pair(index, 1));
}

int CLUSTER::get_histogram(vector<int>& index)
// Get the value of the histogram with the given index vector.  Return 0 if no entry is found.
{

  if ( POP_HIST.find(index) != POP_HIST.end() ) 
	 return(POP_HIST[index]);
  else
	 return(0);
}

void CLUSTER::set_cheby_vals(vector<PAIRS> &FF_2BODY)
// Calculate Chebyshev xmin, xmax, xavg.
{
  for ( int i = 0; i < NPAIRS; i++ ) 
  {
	 int j;
	 for ( j = 0; j < FF_2BODY.size(); j++ ) 
	 {
		if ( ATOM_PAIRS[i] == FF_2BODY[j].PRPR_NM ) 
		  break;
	 }
	 double lambda = 0.0;
	 if ( j < FF_2BODY.size() )
		lambda = FF_2BODY[j].LAMBDA;
	 else
		EXIT_MSG("In function set_cheby_vals, could not find a match for atom pair " + ATOM_PAIRS[i]);

	 Cheby::set_cheby_params(S_MINIM[i], S_MAXIM[i], lambda, FF_2BODY[0].CHEBY_TYPE, X_MINIM[i], X_MAXIM[i], X_DIFF[i], X_AVG[i]);
  }
}


void CLUSTER::map_indices_int(vector<int> & atom_type_idx, vector<int> & pair_map) 
// Matches the allowed powers to the ij, ik, il... type pairs formed from the atoms  ai, aj, ak, al
// This is based on reordering atoms, not pairs.  Pairs can't be interchanged independently, but atoms can.
// "pair_index" serves as the map between the powers in the cluster and the powers between the atoms.
{
  using Pair_int_int = pair<int,int>;

  static vector <Pair_int_int> 	ATOM_TYPE_AND_INDEX;	// [a1/a2/a3/a4 atom pair chemistry][index of pair (from 1-6)]
  //vector <Pair_int_int> 	sorted_index(natoms);	// Sorted to match atom_type_idx.
  static vector<bool> used;
  static vector<int> rev_index;
  static vector<vector<int>> pair_index;
  //static bool called_before = false;

  int npairs = NPAIRS;
  int natoms = NATOMS;

  //if ( ! called_before )
  {
	 // called_before = true;
	 
	 pair_index.resize(natoms);
	 used      .resize(natoms);
	 rev_index .resize(natoms);
	 
	 ATOM_TYPE_AND_INDEX.resize(natoms);

	 for ( int j = 0; j < natoms; j++ ) 
		pair_index[j].resize(natoms);
  }
	
  if ( pair_map.size() != npairs )
  {
  	cout << pair_map.size() << " " << npairs << endl;
	 EXIT_MSG("Wrong dimension for pair map");
  }

  if ( atom_type_idx.size() != natoms ) 
	 EXIT_MSG("Wrong dimension for natoms");

  vector<int> &indices = ATOM_INDICES;
	  
  int count = 0;
  for ( int j = 0; j < natoms; j++ ) 
  {
	 pair_index[j][j] = -1;
	 
	 for ( int k = j + 1; k < natoms; k++ ) 
	 {
		pair_index[j][k] = count;
		pair_index[k][j] = count;
		count++;
	 }
  }

  for(int m = 0; m< natoms; m++)
  { 				
	 ATOM_TYPE_AND_INDEX[m].first  = atom_type_idx[m];					// 
	 ATOM_TYPE_AND_INDEX[m].second = m;	
  }

  // Match the permutation of the atom indices in atom_type_idx to
  // those in the cluster indices[]
  fill(used.begin(), used.end(), false);
  fill(rev_index.begin(), rev_index.end(), -1);

  for ( int m = 0; m < natoms; m++ )
	 for ( int n = 0; n < natoms; n++ ) 
	 {
		if ( ATOM_TYPE_AND_INDEX[m].first == indices[n] && ! used[n] ) 
		{
		  //sorted_index[n] = ATOM_TYPE_AND_INDEX[m];
		  rev_index[ATOM_TYPE_AND_INDEX[m].second] = n;
		  used[n] = true;
		  break;
		}
	 }

  
  for ( int j = 0; j < natoms; j++ ) 
	 if ( rev_index[j] < 0 ) 
		EXIT_MSG("Bad reverse index");
		 
  int m = 0;
  for(int i = 0; i < natoms; i++)
  {
	 int ii = rev_index[i];
	 for ( int j = i + 1; j < natoms; j++ ) 
	 {

		int jj = rev_index[j];
		pair_map[m] = pair_index[ii][jj];

		if ( pair_map[m] < 0 || pair_map[m] > npairs - 1 ) 
		  EXIT_MSG("Permutation power map has a bad value");
		++m;
	 }
  }

  // Double check uniqueness of each value in pair_map.
  for ( int m = 0; m < npairs; m++ ) 
  {
	 for ( int m1 = 0; m1 < npairs; m1++ ) 
	 {
		if ( m1 == m ) 
			continue;
		if ( pair_map[m1] == pair_map[m] )
		  EXIT_MSG("Permutation power map had a repeated value" );
	 }
  }
}

