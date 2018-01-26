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
						  QUAD_MAP.insert(make_pair(FULL_TEMP_STR,p));	
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
			 PAIR_QUADRUPLETS[TEMP_INT].INDX = TEMP_INT;	// Index for current triplet type
								
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
