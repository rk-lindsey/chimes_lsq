#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes
#include<string>
#include<limits>	// Help with handling of over/underflow
#include <algorithm> // For specials vector-related functions (i.e. permute)

using namespace std;

#include "functions.h"

typedef pair<string,int> test;
static bool check_pairs(test a, test b)
{
	return (a.first == b.first);
}


void CLUSTER::build(int cheby_order)
// Build a set of interactions for a quad.
// Figure out the allowed pair quadruplet powers. Here are some rules and considerations:
//
// 1. Powers start from zero. So, if order is specified to be 2, polynomial powers range
//    from 0 to n-1, NOT 1 to n!
//
// 2. At least three pairs must have non-zero powers for the interaction to truly correspond
//    to 4-body interactions
//
// 3. Non-uniqueness upon atom permutation must be taken into consideration. 
//    Pairs permutations occur as a consequence of atom permutations.
{
  vector<int> powers(NPAIRS);
					
  // Loops start at zero for a trick to speed up check for powers > 0


  build_loop(0, cheby_order, powers) ;
				
  N_TRUE_ALLOWED_POWERS = UNIQUE_POWERS.size() ;
  N_ALLOWED_POWERS = PARAM_INDICES.size();
}

void CLUSTER::build_loop(int indx, int cheby_order, vector<int> powers)
// Recursive function implements nested loops over pair powers.
{

  if ( indx < NPAIRS ) 
  {
	 for(int pair_pow=0; pair_pow < cheby_order; pair_pow++)
	 {
		powers[indx] = pair_pow ;
		build_loop(indx+1, cheby_order, powers) ;
	 }

  } 
  else 
  {
	 int THRESHOLD = 0;
	 vector<bool> represented(NATOMS,false) ;
	 int count = 0 ;
	 for ( int i = 0 ; i < NATOMS ; i++ ) 
	 {
		for ( int j = i + 1 ; j < NATOMS ; j++ ) 
		{
		  if ( powers[count++] > 0 ) {
			 represented[i] = true ;
			 represented[j] = true ;
		  }
		}
	 }

	 
	 // All atoms must be represented with non-zero powers.
	 // Note:  the original 4-body code also had the requirement of 3 non-zero powers.
	 // An interaction of xij * xkl has 4 atoms represented, but only 2 non-zero powers.
	 // I think this should be included in the 4 body interaction.  Note that this term would
	 // not be included in a 3 body interaction. (LEF)
	 for ( int i = 0 ; i < NATOMS ; i++ ) 
	 {
		if ( ! represented[i] ) 
		  return ;
	 }
	 // Store all valid atom index permutations 
	 store_permutations(powers) ;
  }
}


void CLUSTER::store_permutations(vector<int> &unsorted_powers) 
{
  map<vector<int>,int>::iterator it ;

  // Avoid searching the ALLOWED_POWERS vector directly for efficiency.
  // Search the ALLOWED_POWERS_MAP instead.
  it = ALLOWED_POWERS_MAP.find(unsorted_powers) ;

  if ( it != ALLOWED_POWERS_MAP.end() ) 
  {
	 // This interaction is already handled.
	 return ;
  } 
  else 
  {
	 UNIQUE_POWERS.push_back(unsorted_powers) ;	 
  }

  // Recursively permute each element and added to the ALLOWED_POWERS.
#if( VERBOSITY > 1)
  if ( RANK == 0 ) 
  {
	 cout << "Permuting atom indices for " ;
	 for ( int i = 0 ; i < NATOMS ; i++ ) cout << ATOM_NAMES[i] << " " ;
	 cout << endl ;
  }
#endif

  vector<int> perm(NATOMS) ;
  for ( int i = 0 ; i < NATOMS ; i++ ) perm[i] = i ;

  int equiv_index = ALLOWED_POWERS.size() ;
  permute_atom_indices(0, ATOM_NAMES, unsorted_powers, perm, UNIQUE_POWERS.size()-1, equiv_index ) ;

}

void CLUSTER::permute_atom_indices(int idx, vector<string> names, vector<int> &unsorted_powers, 
											  vector<int> perm, int unique_index, int equiv_index) 
{

  if ( idx < names.size() ) {
	 // Count up the number of identical atoms.
	 int count = 1 ;
	 for ( int i = idx ; i < names.size() - 1 ; i++ )
	 {
		if ( names[i+1] == names[i] ) {
		  ++count ;
		}
		else {
		  break ;
		}
	 }
	 do {
		// Recursive call to generate permutations for identical atoms.
		permute_atom_indices(idx + count, names, unsorted_powers, perm,
									unique_index, equiv_index) ;
	 } while ( std::next_permutation(perm.data() +idx, perm.data() + idx+count) ) ;
  } else {
	 // Store the generated permutation.
	 // Just love c++ syntax.... don't you ?
	 vector<vector<int>> pow_mat(NATOMS, vector<int>(NATOMS)) ;
	 vector<int> perm_powers(NPAIRS) ;

	 // Generate the permutation transformation matrix.
	 int count = 0 ;
	 for ( int i = 0 ; i < NATOMS ; i++ ) {
		pow_mat[i][i] = 0 ;
		for ( int j = i + 1 ; j < NATOMS ; j++ ) {
		  pow_mat[i][j] = unsorted_powers[count++] ;
		}
	 }
	 /// Fill in other side of the diagonal
	 for ( int k = 0 ; k < NATOMS ; k++ ) {
		for ( int l = 0 ; l < k ; l++ ) {
		  pow_mat[k][l] = pow_mat[l][k] ;
		}
	 }

	 // Permute the powers.
	 count = 0 ;
	 for ( int i = 0 ; i < NATOMS ; i++ ) {
		pow_mat[i][i] = 0 ;
		for ( int j = i + 1 ; j < NATOMS ; j++ ) {
		  perm_powers[count++] = pow_mat[ perm[i] ] [ perm[j] ] ;
		}
	 }

#if(VERBOSITY > 1)
	 if ( RANK == 0 ) 
	 {
		cout << "Permutation: " 
		for ( int j = 0 ; j < NATOMS ; j++ ) {
		  cout << perm[j] << " " ;
		}
		cout << endl ;

		cout << "Permuted powers: " ;
		for ( int j = 0 ; j < NPAIRS ; j++ ) {
		  cout << perm_powers[j] << " " ;
		}
		cout << endl ;
	 }
#endif
	 map<vector<int>,int>::iterator it ;

	 it = ALLOWED_POWERS_MAP.find(perm_powers) ;

	 if ( it == ALLOWED_POWERS_MAP.end() ) {
		ALLOWED_POWERS.push_back(perm_powers) ;
		ALLOWED_POWERS_MAP.insert( std::pair<vector<int>,int>(perm_powers, 1) ) ;
		PARAM_INDICES.push_back(unique_index) ;
		EQUIV_INDICES.push_back(equiv_index) ;
	 } 
  }
}

void CLUSTER::print()
{
  if ( RANK == 0 ) 
  {
	 cout << "		" <<  INDX;
	 for(int j=0; j < NPAIRS ; j++)
		cout  << "  " << ATOM_PAIRS[j];					
	 cout<< ": Number of unique sets of powers: " << N_TRUE_ALLOWED_POWERS << " (" << N_ALLOWED_POWERS << " total)..." << endl; 

	 cout << "		     index  |  powers  |  equiv index  |  param index  " << endl;
	 cout << "		   ----------------------------------------------------" << endl;					
				
	 for(int j=0; j<ALLOWED_POWERS.size(); j++)
	 {

		cout << "		      " << setw(6) << fixed << left << j << " ";
		
		for(int k=0; k<NPAIRS; k++)
		  cout << " " << setw(2) << fixed << left << ALLOWED_POWERS[j][k] << " ";
						
		cout << "       " << setw(8) << EQUIV_INDICES[j] << " ";
		cout << "       " << setw(8) << PARAM_INDICES[j] << endl; 
		
	 }
  }
}

void CLUSTER::print_special(ofstream &header, string QUAD_MAP_REVERSE, string output_mode)
{
  if ( output_mode == "S_MINIM" ) 
  {
	 if( S_MINIM[0] >= 0) 
		header << INDX << " " << QUAD_MAP_REVERSE << " " ;
	 for ( int i = 0 ; i < NPAIRS ; i++ ) 
		header << ATOM_PAIRS[i] << " " ;

	 header << fixed << setprecision(5) ;

	 for ( int i = 0 ; i < NPAIRS ; i++ ) 
		header << S_MINIM[i] << " " ;

	 header << endl;						
  } 
  else if ( output_mode == "S_MAXIM" ) 
  {
	 if ( S_MAXIM[0] >= 0)
		header << INDX << " " << QUAD_MAP_REVERSE << " " ;
	 for ( int i = 0 ; i < NPAIRS ; i++ ) 
		header << ATOM_PAIRS[i] << " " ;

	 header << fixed << setprecision(5) ;

	 for ( int i = 0 ; i < NPAIRS ; i++ ) 
		header << S_MAXIM[i] << " " ;

	 header << endl;						
  }						
  else 
  {
	 cout << "Bad special parameter mode found\n" ;
	 exit(1) ;
  }
	 
}

void CLUSTER::print_header(ofstream &header)
{
  header << INDX << "  ";
  for(int m=0; m<NPAIRS; m++)
  {
	 header << ATOM_PAIRS[m];
	 if(m<NPAIRS-1)
		header << " ";
  }
  header << ": " << N_TRUE_ALLOWED_POWERS << " parameters, " << N_ALLOWED_POWERS << " total parameters "<< endl;	
	
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

  header << endl;
}

void CLUSTER_LIST::link(vector<CLUSTER>& cluster_vec) 
// Link a CLUSTER_LIST to an allocated vector of clusters.
{
  VEC.resize(cluster_vec.size() ) ;
  
  NCLUSTERS = VEC.size() ;

  for ( int j = 0 ; j < VEC.size() ; j++ ) {
	 VEC[j] = &cluster_vec[j] ;
  }
}

void CLUSTER_LIST::link(vector<QUADRUPLETS>& quad_vec) 
// Link a CLUSTER_LIST to an allocated vector of quadruplets
{
  VEC.resize(quad_vec.size() ) ;

  NCLUSTERS = VEC.size() ;

  for ( int j = 0 ; j < VEC.size() ; j++ ) {
	 VEC[j] = &quad_vec[j] ;
  }
}


void CLUSTER_LIST::link(vector<TRIPLETS>& trip_vec) 
// Link a CLUSTER_LIST to an allocated vector of triplets
{
  VEC.resize(trip_vec.size() ) ;

  NCLUSTERS = VEC.size() ;

  for ( int j = 0 ; j < VEC.size() ; j++ ) {
	 VEC[j] = &trip_vec[j] ;
  }
}

  
void CLUSTER_LIST::build_maps(vector<struct PAIRS> &atom_pairs)
{
  int npair = atom_pairs.size() ;

  int cluster_npair = VEC[0]->NPAIRS ;
  vector<int> pair_index(cluster_npair) ;

  build_maps_loop(0, pair_index, atom_pairs) ;
				
  MAP_REVERSE.clear();

  map<string, int>::iterator it;

  for(it = MAP.begin(); it != MAP.end(); it++)
	 MAP_REVERSE.insert(make_pair(it->second,it->first));
}

void CLUSTER_LIST::build_maps_loop(int index, vector<int> pair_index, vector<struct PAIRS> &atom_pairs)
{
  int cluster_npair = VEC[0]->NPAIRS ;

  if ( index < cluster_npair ) 
  {
	 for(int i=0; i< atom_pairs.size() ; i++) 
	 {
		pair_index[index] = i ;
		build_maps_loop(index+1, pair_index, atom_pairs) ;
	 }
  } 
  else
  {
	 // Done looping.  Actually build the loop
	 vector<string> temp_str(cluster_npair);	
	 bool real_cluster ;
	 string full_tmp_str;	

	 for ( int j = 0 ; j < cluster_npair ; j++ ) 
	 {
		temp_str[j] = atom_pairs[ pair_index[j] ].ATM1TYP; 
		temp_str[j].append(atom_pairs[ pair_index[j] ].ATM2TYP) ;
	 }										
										
	 for(int p=0; p < NCLUSTERS ; p++)
	 {
		sort(temp_str.begin(),temp_str.end());
		do
		{
		  real_cluster = true;
												
		  full_tmp_str = "";

		  for (int s = 0 ; s< cluster_npair ; s++)
		  {
			 if(temp_str[s] != VEC[p]->ATOM_PAIRS[s])
				real_cluster = false;
			 
			 full_tmp_str += temp_str[s];
		  }
												
		  if( real_cluster )
		  {													
			 MAP.insert(make_pair(full_tmp_str,p));	
			 break;
		  }
		} while (next_permutation(temp_str.begin(),temp_str.end()));
	 }
  }
}
				
void CLUSTER_LIST::build_fast_maps(vector<string>& ATOM_CHEMS)
{
  int natoms = VEC[0]->NATOMS ;

  if(RANK == 0) 
  {
	 if ( natoms == 4 ) 
		cout << endl << "Generating fast maps for atom quadruplets: " << endl;
	 else if ( natoms == 3 ) 
		cout << endl << "Generating fast maps for atom triplets: " << endl;
	 else if ( natoms == 5 ) 
		cout << endl << "Generating fast maps for atom quintuplets: " << endl;
				
  string 	TMP_QUAD_SIXLET = "";
  string 	TMP_QUAD_PAIR   = "";
				
  int 	ATOM_QUAD_ID_INT;
  int		CURR_QUAD_IDX = 0;

  vector<int> TMP_QUAD_SET(4);

  vector<int> atom_index(natoms) ;

  build_fast_maps_loop(0, atom_index, ATOM_CHEMS) ;
  }
}

void CLUSTER_LIST::build_fast_maps_loop(int index, vector<int> atom_index, vector<string>& ATOM_CHEMS)
{
  // Create a list of each posible combination of atom quadruplets, i through
  // l are given in ascending order
  int natoms = VEC[0]->NATOMS ;
  
  if ( index < natoms ) 
  {
	 for ( int i = 0 ; i < ATOM_CHEMS.size() ; i++ ) {
		atom_index[index] = i ;
		build_fast_maps_loop(index+1, atom_index, ATOM_CHEMS) ;
	 }
  }
  else
  {
	 string tmp_pair = "" ;
	 vector<int> tmp_set(natoms) ;

	 for ( int i = 0 ; i < natoms ; i++ ) {
		for ( int j = i + 1 ; j < natoms ; j++ ) {
			 // Determine the corresponding quadruplet force field type

		  tmp_pair += ATOM_CHEMS[ atom_index[i] ] + ATOM_CHEMS[ atom_index[j] ] ;
		}
		tmp_set[i] = atom_index[i] ;
	 }
	
	 // Always construct ATOM_QUAD_ID_INT lookup key based on atom types in decending order
	 sort(tmp_set.begin(), tmp_set.end());
	 reverse(tmp_set.begin(), tmp_set.end());

	 int atom_id_int = make_id_int(tmp_set) ;
	 
	 INT_MAP.insert(make_pair(atom_id_int,MAP[tmp_pair]));	
	 INT_MAP_REVERSE.insert(make_pair(MAP[tmp_pair], atom_id_int));	
								
	 if(RANK == 0)
	 {
		cout << "		";
		cout<< "Atom type idxs: ";

		for ( int i = 0 ; i < natoms ; i++ ) {
		  cout<< fixed << setw(2) << right << atom_index[i] ;
		}
		
		if ( natoms == 4 ) 
		{
		  cout<< " Quadruplet name: "           << setw(12) << right << tmp_pair ;
		  cout<< " Unique quadruplet index: "   << setw(4)  << right << INT_MAP[atom_id_int] << endl;
		} 
		else if ( natoms == 3 ) 
		{
		  cout<< " Triplet name: "           << setw(12) << right << tmp_pair ;
		  cout<< " Unique triplet index: "   << setw(4)  << right << INT_MAP[atom_id_int] << endl;
		}
		else if ( natoms == 5 ) 
		{
		  cout<< " Quintuplet name: "           << setw(12) << right << tmp_pair ;
		  cout<< " Unique quintuplet index: "   << setw(4)  << right << INT_MAP[atom_id_int] << endl;
		}
	 }
  }
}


int CLUSTER_LIST::make_id_int(vector<int>& index)
// Returns a unique ID number for a cluster of atoms.
{
  int sum = 0 ;
  int weight = 1 ;
  for ( int i = 0 ; i < index.size() ; i++ ) {
	 sum += weight * (index[i]+1) ;
	 weight *= MAX_ATOM_TYPES ;
  }
  return(sum) ;
}

void CLUSTER_LIST::build_pairs(vector<string> ATOM_CHEMS, vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP)
// Build the pair variables for all of the clusters
{
			
  int TEMP_INT = 0;				// Will hold pair quadruplet index
  
  int natoms = VEC[0]->NATOMS ;
  vector<int> atom_index(natoms) ;
  int count = 0 ;

  build_pairs_loop(0, atom_index, ATOM_CHEMS, ATOM_PAIRS, PAIR_MAP, count) ;
}

void CLUSTER_LIST::build_pairs_loop(int index, vector<int> atom_index, 
												vector<string> ATOM_CHEMS, vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP, int &count)
// Loop over atom types in building the CLUSTER ATOM_PAIRS list.
{
  int NATMTYP = ATOM_CHEMS.size() ;

  if ( index < NATMTYP ) {
	 for ( int i = 0 ; i < NATMTYP ; i++ ) 
	 {
		atom_index[index] = i ;
		build_pairs_loop(index+1, atom_index, ATOM_CHEMS, ATOM_PAIRS, PAIR_MAP, count) ;
	 }
  }
  else 
  {
	 VEC[count]->INDX = count ;	// Index for current triplet type
								
	 // Save the names of each atom type in the pair.
	 int natoms = VEC[0]->NATOMS ;
	 int count2 = 0 ;
	 for ( int i = 0 ; i < natoms ; i++ ) {
		VEC[count]->ATOM_NAMES[i] = ATOM_CHEMS[ atom_index[i] ] ;
		for ( int j = i + 1 ; j < natoms ; j++ ) {
		  VEC[count]->ATOM_PAIRS[count2] = ATOM_CHEMS[ atom_index[i] ] 
			 + ATOM_CHEMS[ atom_index[j] ] ;
		  VEC[count]->ATOM_PAIRS[count2] = ATOM_PAIRS[ PAIR_MAP[ VEC[count]->ATOM_PAIRS[count2] ] ].PRPR_NM ;
		  count2++ ;
		}
	 }

	 if(RANK==0)
	 {							
		if ( natoms == 4 ) 
		{
		  cout << "Made the following quadruplets: " << count << " ";
		} else if ( natoms == 3 ) 
		{
		  cout << "Made the following triplets: " << count << " ";
		} else if ( natoms == 5 )
		{
		  cout << "Made the following quintuplets: " << count << " ";
		}
		int npairs = VEC[0]->NPAIRS ;
		for(int m=0; m< npairs; m++) 
		  cout << VEC[count]->ATOM_PAIRS[m] << " ";
		cout << endl;		
	 }
	 count++;	
  }			
}
