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

void QUADRUPLETS::print()
{
  cout << "		" <<  QUADINDX;
  for(int j=0; j<6; j++)
	 cout  << "  " << ATOM_PAIRS[j];					
  cout<< ": Number of unique sets of powers: " << N_TRUE_ALLOWED_POWERS << " (" << N_ALLOWED_POWERS << " total)..." << endl; 

  cout << "		     index  |  powers  |  equiv index  |  param index  " << endl;
  cout << "		   ----------------------------------------------------" << endl;					
				
  for(int j=0; j<ALLOWED_POWERS.size(); j++)
  {

	 cout << "		      " << setw(6) << fixed << left << j << " ";
						
	 for(int k=0; k<6; k++)
		cout << " " << setw(2) << fixed << left << ALLOWED_POWERS[j][k] << " ";
						
	 cout << "       " << setw(8) << EQUIV_INDICIES[j] << " ";
	 cout << "       " << setw(8) << PARAM_INDICIES[j] << endl; 
					
  }
}

void build_quad_maps(map<string,int> &QUAD_MAP, map<int,string> &QUAD_MAP_REVERSE, vector<struct PAIRS> &ATOM_PAIRS, int NPAIR,
                     vector<QUADRUPLETS> PAIR_QUADRUPLETS, int NQUAD)
{
  bool REAL_QUAD;
			
  vector<string> TEMP_STR(6);	
  string FULL_TEMP_STR;	
				
  for(int i=0; i<NPAIR; i++)
  {
	 for(int j=0; j<NPAIR; j++)
	 {
		for(int k=0; k<NPAIR; k++)
		{
		  for(int l=0; l<NPAIR; l++)
		  {
			 for(int m=0; m<NPAIR; m++)
			 {
				for(int n=0; n<NPAIR; n++)
				{
				  TEMP_STR[0] = ATOM_PAIRS[i].ATM1TYP; TEMP_STR[0].append(ATOM_PAIRS[i].ATM2TYP);
				  TEMP_STR[1] = ATOM_PAIRS[j].ATM1TYP; TEMP_STR[1].append(ATOM_PAIRS[j].ATM2TYP);
				  TEMP_STR[2] = ATOM_PAIRS[k].ATM1TYP; TEMP_STR[2].append(ATOM_PAIRS[k].ATM2TYP);
				  TEMP_STR[3] = ATOM_PAIRS[l].ATM1TYP; TEMP_STR[3].append(ATOM_PAIRS[l].ATM2TYP);
				  TEMP_STR[4] = ATOM_PAIRS[m].ATM1TYP; TEMP_STR[4].append(ATOM_PAIRS[m].ATM2TYP);
				  TEMP_STR[5] = ATOM_PAIRS[n].ATM1TYP; TEMP_STR[5].append(ATOM_PAIRS[n].ATM2TYP);
										
										

				  // Compare every possible permutation of these pairs with those for the quadruplet ff types
						
						
				  for(int p=0; p<NQUAD; p++)
				  {
					 sort(TEMP_STR.begin(),TEMP_STR.end());
					 do
					 {
						REAL_QUAD = true;
												
						FULL_TEMP_STR = "";

						for (int s=0; s<6; s++)
						{
						  if(TEMP_STR[s] != PAIR_QUADRUPLETS[p].ATOM_PAIRS[s])
							 REAL_QUAD = false;
													
						  FULL_TEMP_STR += TEMP_STR[s];
						}
												
												
												
						if(REAL_QUAD)
						{													
						  QUAD_MAP        .insert(make_pair(FULL_TEMP_STR,p));	
						  QUAD_MAP_REVERSE.insert(make_pair(p,FULL_TEMP_STR));														
						  break;
						}

												
					 } while (next_permutation(TEMP_STR.begin(),TEMP_STR.end()));
/* // THIS IS THE CULPRIT											
	if(!REAL_QUAD)
	{
	QUAD_MAP        .insert(make_pair(FULL_TEMP_STR,-1*p-1));
	}
*/											
											
				  }
				}
			 }
		  }
		}
	 }
  }
				
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Remove excluded types (NOT CURRENTLY SUPPORTED)
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////					

				
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set up quadruplet REVERSE maps... 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////	
				
  QUAD_MAP_REVERSE.clear();

  map<string, int>::iterator it;

  for(it = QUAD_MAP.begin(); it != QUAD_MAP.end(); it++)
	 QUAD_MAP_REVERSE.insert(make_pair(it->second,it->first));
}


void build_fast_quad_maps(map<string,int> QUAD_MAP, map<int,int> INT_QUAD_MAP, map<int,int> INT_QUAD_MAP_REVERSE, vector<string> ATOM_CHEMS)
{
  if(RANK == 0)
	 cout << endl << "Generating fast maps for atom quadruplets: " << endl;
				
  string 	TMP_QUAD_SIXLET = "";
  string 	TMP_QUAD_PAIR   = "";
				
  int 	ATOM_QUAD_ID_INT;
  int		CURR_QUAD_IDX = 0;

  vector<int> TMP_QUAD_SET(4);

  // Create a list of each posible combination of atom quadruplets, i through
  // l are given in ascending order
				
  for(int i=0; i<ATOM_CHEMS.size(); i++)
  {
	 for(int j=i; j<ATOM_CHEMS.size(); j++)
	 {
		for(int k=j; k<ATOM_CHEMS.size(); k++)
		{
		  for(int l=k; l<ATOM_CHEMS.size(); l++)
		  {
			 // A given set of i-through-l defines a unique set of 4 atoms.
			 // Determine the corresponding quadruplet force field type
								
			 TMP_QUAD_SIXLET = "";
			 TMP_QUAD_PAIR   = "";
								
			 TMP_QUAD_PAIR    = ATOM_CHEMS[i] + ATOM_CHEMS[j]; 	// Define an ij pair
			 TMP_QUAD_SIXLET += TMP_QUAD_PAIR;					// Append pair to string
								
			 TMP_QUAD_PAIR    = ATOM_CHEMS[i] + ATOM_CHEMS[k]; 	// Define an ij pair
			 TMP_QUAD_SIXLET += TMP_QUAD_PAIR;					// Append pair to string
								
			 TMP_QUAD_PAIR    = ATOM_CHEMS[i] + ATOM_CHEMS[l]; 	// Define an ij pair
			 TMP_QUAD_SIXLET += TMP_QUAD_PAIR;					// Append pair to string
								
			 TMP_QUAD_PAIR    = ATOM_CHEMS[j] + ATOM_CHEMS[k]; 	// Define an ij pair
			 TMP_QUAD_SIXLET += TMP_QUAD_PAIR;					// Append pair to string
								
			 TMP_QUAD_PAIR    = ATOM_CHEMS[j] + ATOM_CHEMS[l]; 	// Define an ij pair
			 TMP_QUAD_SIXLET += TMP_QUAD_PAIR;					// Append pair to string
								
			 TMP_QUAD_PAIR    = ATOM_CHEMS[k] + ATOM_CHEMS[l]; 	// Define an ij pair
			 TMP_QUAD_SIXLET += TMP_QUAD_PAIR;					// Append pair to string
								
			 TMP_QUAD_SET[0] = i;
			 TMP_QUAD_SET[1] = j;
			 TMP_QUAD_SET[2] = k;
			 TMP_QUAD_SET[3] = l;
								
			 // Always construct ATOM_QUAD_ID_INT lookup key based on atom types in decending order
								
			 sort   (TMP_QUAD_SET.begin(), TMP_QUAD_SET.end());
			 reverse(TMP_QUAD_SET.begin(), TMP_QUAD_SET.end());

			 ATOM_QUAD_ID_INT = make_quad_id_int(TMP_QUAD_SET[0], TMP_QUAD_SET[1], TMP_QUAD_SET[2], TMP_QUAD_SET[3]) ;

			 INT_QUAD_MAP        .insert(make_pair(ATOM_QUAD_ID_INT,QUAD_MAP[TMP_QUAD_SIXLET]));	
			 INT_QUAD_MAP_REVERSE.insert(make_pair(QUAD_MAP[TMP_QUAD_SIXLET], ATOM_QUAD_ID_INT));	
								
			 if(RANK == 0)
			 {
				cout << "		";
				cout<< "Atom type idxs: ";
				cout<< fixed << setw(2) << right << i;
				cout<< fixed << setw(2) << right << j;
				cout<< fixed << setw(2) << right << k;
				cout<< fixed << setw(2) << right << l;
				cout<< " Quadruplet name: "           << setw(12) << right << TMP_QUAD_SIXLET;
				cout<< " Unique quadruplet index: "   << setw(4)  << right << INT_QUAD_MAP[ATOM_QUAD_ID_INT] << endl;
			 }
								
			 CURR_QUAD_IDX++;

		  }
		}
	 }
  }

}

void build_quad_pairs(vector<QUADRUPLETS>& PAIR_QUADRUPLETS, int NATMTYP, vector<string> ATOM_CHEMS, vector<PAIRS> ATOM_PAIRS, map<string,int> PAIR_MAP)
// Build the pair variables for all of the quads.
{
			
  int TEMP_INT = 0;				// Will hold pair quadruplet index
			
  for(int i=0; i<NATMTYP; i++)
  {
	 for(int j=i; j<NATMTYP; j++)
	 {
		for(int k=j; k<NATMTYP; k++)
		{
		  for(int l=k; l<NATMTYP; l++)
		  {
			 PAIR_QUADRUPLETS[TEMP_INT].QUADINDX = TEMP_INT;	// Index for current triplet type
								
			 // Save the names of each atom type in the pair.

			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_NAMES[0] = ATOM_CHEMS[i] ;
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_NAMES[1] = ATOM_CHEMS[j] ;
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_NAMES[2] = ATOM_CHEMS[k] ;
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_NAMES[3] = ATOM_CHEMS[l] ;

			 // Save the first atom type in each pair
								
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[0] = ATOM_CHEMS[i]; // ij
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[1] = ATOM_CHEMS[i]; // ik
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[2] = ATOM_CHEMS[i]; // il
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[3] = ATOM_CHEMS[j]; // jk
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[4] = ATOM_CHEMS[j]; // jl
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[5] = ATOM_CHEMS[k]; // kl
								
			 // Save the second atom type in each pair
								
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[0] += ATOM_CHEMS[j]; // ij
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[1] += ATOM_CHEMS[k]; // ik
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[2] += ATOM_CHEMS[l]; // il
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[3] += ATOM_CHEMS[k]; // jk
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[4] += ATOM_CHEMS[l]; // jl
			 PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[5] += ATOM_CHEMS[l]; // kl								
								
			 // Now save the "proper" (ordered) name of the pair	

								 
			 for(int m=0; m<6; m++) 
				PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[m] = ATOM_PAIRS[ PAIR_MAP[ PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[m]] ].PRPR_NM;	


			 if(RANK==0)
			 {							
				cout << "Made the following quadruplets: " << TEMP_INT << " ";
				for(int m=0; m<6; m++) 
				  cout << PAIR_QUADRUPLETS[TEMP_INT].ATOM_PAIRS[m] << " ";
				cout << endl;		
			 }
			 TEMP_INT++;	
		  }					
		}
	 }
  }			
}
