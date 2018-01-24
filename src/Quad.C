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


void QUADRUPLETS::init()
// Set initial values for the QUADRUPLETS struct.
{
  for(int j=0; j<6; j++)
  {
	 S_MINIM_4B[j] = -1;
	 S_MAXIM_4B[j] = -1;
  }
				
  FORCE_CUTOFF.TYPE = FCUT_TYPE::CUBIC;	
}	

void QUADRUPLETS::build(int cheby_4b_order)
// Build a set of interactions for a quad.
{
  vector<int> STORED_SORTED_POWERS_EQVS;
  vector<int> UNSORTED_POWERS(6);
					
  // The point of this loop is to determine true unique power sets. Whether or not a given set of 6 powers is 
  // unique depends on the number of unique pairs within the corresponding set of pair types. We will determine
  // this first.
					
  // The data types below will be used to determine uniqueness of power sets
					
  typedef pair<string,int>			PAIR_TYPE_AND_INDEX;
  vector <PAIR_TYPE_AND_INDEX>	PAIR_TYPE_AND_INDEX_VEC(6);	// Used to determine uniqueness of power sets
  vector <PAIR_TYPE_AND_INDEX>	PAIR_TYPE_AND_INDEX_TMP(6);	// A temporary data structure
					
  for(int m=0; m<6; m++)
  {					
	 PAIR_TYPE_AND_INDEX_VEC[m].first  = ATOM_PAIRS[m];
	 PAIR_TYPE_AND_INDEX_VEC[m].second = m;	
  }

  sort (PAIR_TYPE_AND_INDEX_VEC.begin(), PAIR_TYPE_AND_INDEX_VEC.end());	// Sort the vector contents... automatically does on the basis of the .first element, preserving "link" between .first and .second
					
  // Since we don't actually want to erase duplicates from the PAIR_TYPE_AND_INDEX_VEC vector, make a temporary copy to use for unique counting
					
  copy(PAIR_TYPE_AND_INDEX_VEC.begin(), PAIR_TYPE_AND_INDEX_VEC.end(), PAIR_TYPE_AND_INDEX_TMP.begin());
				
  PAIR_TYPE_AND_INDEX_TMP.erase( unique( PAIR_TYPE_AND_INDEX_TMP.begin(), PAIR_TYPE_AND_INDEX_TMP.end(),check_pairs ), PAIR_TYPE_AND_INDEX_TMP.end() );	// Removes duplicated from SORTED vector, based only on the .first element (i.e. uses our check_pairs function)

  //int N_UNIQUE_PAIRS	= PAIR_TYPE_AND_INDEX_VEC.size();
  int N_UNIQUE_PAIRS	 = PAIR_TYPE_AND_INDEX_TMP.size();
					
  // Loops start at zero for a trick to speed up check for powers > 0
					
  int THRESHOLD = 0;
  bool FOUND_I, FOUND_J, FOUND_K, FOUND_L;
					
  for(int pair1_pow=0; pair1_pow<cheby_4b_order; pair1_pow++)
  {
	 for(int pair2_pow=0; pair2_pow<cheby_4b_order; pair2_pow++)
	 {
		for(int pair3_pow=0; pair3_pow<cheby_4b_order; pair3_pow++)
		{
		  for(int pair4_pow=0; pair4_pow<cheby_4b_order; pair4_pow++)
		  {
			 for(int pair5_pow=0; pair5_pow<cheby_4b_order; pair5_pow++)
			 {
				for(int pair6_pow=0; pair6_pow<cheby_4b_order; pair6_pow++)
				{
				  // Before we go any further, make sure that:
				  // 1. at least 3 powers are greater than zero
				  // 2. all 4 atoms (i.e. i, j, k, and l) are represented in the non-zero powers
											
				  THRESHOLD = 0;
											
				  FOUND_I = FOUND_J = FOUND_K = FOUND_L = false;										
											
				  if (pair1_pow == 0)
					 THRESHOLD++;
				  else
					 FOUND_I = FOUND_J = true;
				  if (pair2_pow == 0)
					 THRESHOLD++;
				  else
					 FOUND_I = FOUND_K = true;
				  if (pair3_pow == 0)
					 THRESHOLD++;
				  else
					 FOUND_I = FOUND_L = true;
				  if (pair4_pow == 0)
					 THRESHOLD++;
				  else
					 FOUND_J = FOUND_K = true;
				  if (pair5_pow == 0)
					 THRESHOLD++;
				  else
					 FOUND_J = FOUND_L = true;
				  if (pair6_pow == 0)
					 THRESHOLD++;
				  else
					 FOUND_K = FOUND_L = true;
											
				  if(THRESHOLD > 3) // Then this is not a valid set of powers. Move on to the next set
					 continue;
											
				  if(!(FOUND_I && FOUND_J && FOUND_K && FOUND_L))	// Then this is not a valid set of powers. Move on to the next set
					 continue;
											
				  UNSORTED_POWERS[0] =	pair1_pow;
				  UNSORTED_POWERS[1] =	pair2_pow;
				  UNSORTED_POWERS[2] =	pair3_pow;
				  UNSORTED_POWERS[3] =	pair4_pow;
				  UNSORTED_POWERS[4] =	pair5_pow;
				  UNSORTED_POWERS[5] =	pair6_pow;
											
				  // Store all valid sixlets and permutations 
				  store_permutations(UNSORTED_POWERS) ;
				}
			 }
		  }
		}
	 }
  }
				
  N_TRUE_ALLOWED_POWERS = UNIQUE_POWERS.size() ;
  N_ALLOWED_POWERS = PARAM_INDICIES.size();
}


void QUADRUPLETS::store_permutations(vector<int> &unsorted_powers) 
{
  vector<vector<int>>::iterator it ;

  it = find( ALLOWED_POWERS.begin(), ALLOWED_POWERS.end(), unsorted_powers) ;

  if ( it != ALLOWED_POWERS.end() ) 
  {
	 // This interaction is already handled.
	 return ;
  } 
  else 
  {
	 UNIQUE_POWERS.push_back(unsorted_powers) ;	 
  }

  // Recursively permute each element and added to the ALLOWED_POWERS.
  cout << "Permuting atom indices for " ;
  for ( int i = 0 ; i < 4 ; i++ ) cout << ATOM_NAMES[i] << " " ;
  cout << endl ;

  vector<int> perm(4) ;
  for ( int i = 0 ; i < 4 ; i++ ) perm[i] = i ;

  int equiv_index = it - ALLOWED_POWERS.begin() ;
  permute_atom_indices(0, ATOM_NAMES, unsorted_powers, perm, UNIQUE_POWERS.size()-1, equiv_index ) ;

}

void QUADRUPLETS::permute_atom_indices(int idx, vector<string> names, vector<int> &unsorted_powers, 
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
		int pow_mat[4][4] ;
		vector<int> perm_powers(6) ;

		// Generate the permutation transformation matrix.
		pow_mat[0][0] = pow_mat[1][1] = pow_mat[2][2] = pow_mat[3][3] = 0 ;
		pow_mat[0][1] = unsorted_powers[0] ;
		pow_mat[0][2] = unsorted_powers[1] ;
		pow_mat[0][3] = unsorted_powers[2] ;
		pow_mat[1][2] = unsorted_powers[3] ;
		pow_mat[1][3] = unsorted_powers[4] ;
		pow_mat[2][3] = unsorted_powers[5] ;

		/// Fill in other side of the diagonal
		for ( int k = 0 ; k < 4 ; k++ ) {
		  for ( int l = 0 ; l < k ; l++ ) {
			 pow_mat[k][l] = pow_mat[l][k] ;
		  }
		}

		
		// Permute the powers.
		perm_powers[0] = pow_mat[ perm[0] ][ perm[1] ] ;
		perm_powers[1] = pow_mat[ perm[0] ][ perm[2] ] ;
		perm_powers[2] = pow_mat[ perm[0] ][ perm[3] ] ;
		perm_powers[3] = pow_mat[ perm[1] ][ perm[2] ] ;
		perm_powers[4] = pow_mat[ perm[1] ][ perm[3] ] ;
		perm_powers[5] = pow_mat[ perm[2] ][ perm[3] ] ;

		cout << "Permutation: " << perm[0] << " " << perm[1] << " "
			  << perm[2] << " " << perm[3] << endl ;

		cout << "Permuted powers: " ;
		for ( int j = 0 ; j < 6 ; j++ ) {
		  cout << perm_powers[j] << " " ;
		}
		cout << endl ;

		vector<vector<int>>::iterator it ;

		it = find( ALLOWED_POWERS.begin(), ALLOWED_POWERS.end(), perm_powers) ;

		if ( it == ALLOWED_POWERS.end() ) {
		  ALLOWED_POWERS.push_back(perm_powers) ;
		  PARAM_INDICIES.push_back(unique_index) ;
		  EQUIV_INDICIES.push_back(equiv_index) ;
		} 
  }
}

