// File containing functions for Chebyshev interaction model.

#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes
#include<string>
#include<limits>	// Help with handling of over/underflow
#include<algorithm> // Used for sorting, etc.
#include "functions.h"
#include "util.h"
#include "Cheby.h"

#ifdef USE_MPI
	#include <mpi.h>
#endif

using namespace std;

#define DEBUG_CHEBY

//////////////////////////////////////////
// Cheby transformation functions
//////////////////////////////////////////

// Note: "Inline" tells the compiler to replace function calls with the actual contents of the function
//       to enhance efficiency.. "undoes" modularity at runtime to speed things up.


inline double Cheby::fix_val(double x)
//Takes care of cheby xformed dist behavior outside of allowed range
{
	// If the chebyshev x value is out of range, set it to a limiting value (-1 or 1)
	// This is done purely to maintain numerical stability.  
	// For 2-body interactions there is an additional repulsion for r < rmin, and
	//   a cutoff r > rmax.
	// For 3-body interactions there should be an inner cutoff function for r < rmin.
	//   and an outer cutoff function for r > rmax.
	
	if ( x < -1.0)
		return -1.0;
	else if ( x > 1.0 )
		return 1.0;								
	else
	  return x;  
}

inline void Cheby::transform(double rlen, double x_diff, double x_avg, double lambda, Cheby_trans cheby_type, double & x, double &exprlen)
// Does the actual cheby distance transformation with precalculated limits.										
// Calculates the transformed variable (x) and the exponential term (exprlen) for Morse interactions.
{
	// Given the atomic distance rlen, the fitting minimum and maximum, the Morse lambda variable,
	//	and the type of Chebyshev approximant, calculate:
	// the chebyshev variable x, 
	switch ( cheby_type ) 
	{
		case Cheby_trans::MORSE:
	  		exprlen = exp(-rlen/lambda) ;
	  		x = (exprlen-x_avg)/x_diff;					// pair distances in morse space, normalized to fit over [-1,1]
	  		break ;
		case Cheby_trans::INVRSE_R:
	  		exprlen = 0.0 ;
	  		x     = (1.0/rlen-x_avg) / x_diff;			// pair distances in r^-1 space, normalized to fit over [-1,1]
	  		break ;
		case Cheby_trans::NONE:
	  		exprlen = 0.0 ;
	  		x = (rlen-x_avg) / x_diff;	
	  		break ;
		default:
	  		cout << "ERROR: Undefined CHBTYPE: " << endl ;
	  		cout << "       Excepted values are \"DEFAULT\", \"INVRSE_R\", or \"MORSE\". " << endl;
	  		exit_run(1);
	}
}

inline int Cheby::get_pair_index(int a1, int a2, const vector<int> &atomtype_idx, int natmtyp, 
											const vector<int> &parent)
// Return the pair type index corresponding to atom a1 and atom a2
{
  int curr_pair_type_idx ;

#ifndef LINK_LAMMPS
  int idx = atomtype_idx[parent[a1]]*natmtyp + atomtype_idx[parent[a2]] ;
  if ( idx < 0 || idx >= INT_PAIR_MAP.size() ) 
  {
	 cout << "Index out of range" << endl ;
  }
  curr_pair_type_idx =  INT_PAIR_MAP[idx] ;
#else
  curr_pair_type_idx =  INT_PAIR_MAP[(atomtype_idx[a1]-1)*natmtyp + (atomtype_idx[SYSTEM.PARENT[a2]]-1)];
#endif
  
  return curr_pair_type_idx ;
}

void Cheby::set_cheby_params(double s_minim, double s_maxim, double lambda, Cheby_trans cheby_type, double &xmin, double &xmax, double &xdiff, double &xavg)
// Calculate Cheby-parameters that do not depend on the interatomic distance.
{
	// Chebyshev polynomials are only defined on the range [-1,1], so transorm the pair distance
	// in a way that allows it to fall along that range. Options are:
	//
	// x = 1/pair_dist				// Inverse r, 
	// x = exp(pair_dist/lambda)	// Morse-type
	// x = pair_dist				// default type
	// 
	// All types are normalized by s_min to s_max range to fall along [-1,1]	
	
	switch ( cheby_type ) 
	{
		case Cheby_trans::MORSE:
			xmin  = exp(-s_maxim/lambda); 
			xmax  = exp(-s_minim/lambda); 
			break ;
			
		case Cheby_trans::INVRSE_R:
	  		xmin = 1.0 / s_maxim ;
	  		xmax = 1.0 / s_minim ;
	  		break ;
			
		case Cheby_trans::NONE:
	  		xmin = s_minim ;
	  		xmax = s_maxim ;
	  		break ;
			
		default:
	  		cout << "ERROR: Undefined CHBTYPE: " << endl ;
	  		cout << "       Excepted values are \"DEFAULT\", \"INVRSE_R\", or \"MORSE\". " << endl;
	  		exit_run(1);
	}
	
	xavg  = 0.5 * (xmin + xmax);				// midpoint of possible pair distances in morse space
	xdiff = 0.5 * (xmax - xmin);				// width of possible pair distances in morse space
}


inline double cheby_var_deriv(double xdiff, double rlen, double lambda, Cheby_trans cheby_type, double exprlen)
// Calculate the derivative of the cheby variable x with respect to rlen.
// exprlen in the exponential in the morse term.
{
  double dx_dr;
	
  switch ( cheby_type ) 
  {
  case Cheby_trans::MORSE:
	 dx_dr =  (-exprlen/lambda)/xdiff;
	 break ;
  case Cheby_trans::INVRSE_R:
	 dx_dr = -1.0/(rlen * rlen * xdiff);
	 break ;
  case Cheby_trans::NONE:
	 dx_dr = 1.0 / xdiff;
	 break ;
  default:
	 dx_dr = 0 ;
	 cout << "Error: bad cheby_type: " << endl;
	 exit_run(1);
  }
  return dx_dr;
}

Cheby_trans Cheby::get_trans_type(string cheby_type) 
// Return the transformation type given an input string.
{
  if ( cheby_type == "MORSE" ) 
	 return Cheby_trans::MORSE ;
  else if ( cheby_type == "INVRSE_R" ) 
	 return Cheby_trans::INVRSE_R ;
  else if ( cheby_type == "DEFAULT" )
	 return Cheby_trans::NONE ;
  else
	 EXIT_MSG("Bad Cheby transformation type: " + cheby_type) ;

  // Not reached.
  return Cheby_trans::NONE ;
}

string Cheby::get_trans_string(Cheby_trans trans) 
// Return the string corresponding to a transformation type.
{
  switch ( trans ) 
  {
  case Cheby_trans::MORSE:
	 return "MORSE" ;
	 break ;
  case Cheby_trans::INVRSE_R:
	 return "INVRSE_R" ;
	 break ;
  case Cheby_trans::NONE:
	 return "NONE" ;
	 break ;
  default:
	 EXIT_MSG("Bad Cheby transformation variable") ;
  }
  // Not reached.
  return("Unknown") ;
}


void Cheby::set_polys(int index, double *Tn, double *Tnd, const double rlen, double x_diff, double x_avg, double SNUM)
// Sets the value of the Chebyshev polynomials (Tn) and their derivatives (Tnd).  Tnd is the derivative
// with respect to the interatomic distance, not the transformed distance (x).
{
  double x = 0 ;
  double exprlen = 0 ;

  // Do the Cheby distance transformation

  PAIRS & ff_2body = FF_2BODY[index] ;

  transform(
	 rlen, 
	 x_diff, 
	 x_avg, 
	 ff_2body.LAMBDA, 
	 ff_2body.CHEBY_TYPE, 
	 x,
	 exprlen) ;
		
#if CHECK_CHEBY_RANGE == 1	
  x = fix_val(x) ;

  // Now change the range, if the user requested
	
  x = x*DERIV_CONST + ff_2body.CHEBY_RANGE_LOW - -1.0*DERIV_CONST;
		
  // Sanity check
	
  if ( x < ff_2body.CHEBY_RANGE_LOW || x > ff_2body.CHEBY_RANGE_HIGH )
  {
	 cout << "ERROR: transformed (3B) x falls outside user-defined range." << endl;
	 cout << "x: " << x << endl;
	 cout << "high/low: " << ff_2body.CHEBY_RANGE_HIGH << " " << ff_2body.CHEBY_RANGE_LOW  << endl;
	 exit_run(0);
  }	
		
#endif
		
  // Generate Chebyshev polynomials by recursion. 
  // 
  // What we're doing here. Want to fit using Cheby polynomials of the 1st kinD[i]. "T_n(x)."
  // We need to calculate the derivative of these polynomials.
  // Derivatives are defined through use of Cheby polynomials of the 2nd kind "U_n(x)", as:
  //
  // d/dx[ T_n(x) = n * U_n-1(x)] 
  // 
  // So we need to first set up the 1st-kind polynomials ("Tn[]")
  // Then, to compute the derivatives ("Tnd[]"), first set equal to the 2nd-kind, then multiply by n to get the der's
	 
  // First two 1st-kind Chebys:
		
  Tn[0] = 1.0;
  Tn[1] = x;
	
  // Start the derivative setup. Set the first two 1st-kind Cheby's equal to the first two of the 2nd-kind

  Tnd[0] = 1.0;
  Tnd[1] = 2.0 * x;
	
  // Use recursion to set up the higher n-value Tn and Tnd's

  for ( int i = 2; i <= SNUM; i++ ) 
  {
	 Tn[i]  = 2.0 * x *  Tn[i-1] -  Tn[i-2];
	 Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2];
  }
	
  // Now multiply by n to convert Tnd's to actual derivatives of Tn

  double dx_dr = DERIV_CONST*cheby_var_deriv(x_diff, rlen, ff_2body.LAMBDA, ff_2body.CHEBY_TYPE, exprlen);

  for ( int i = SNUM; i >= 1; i-- ) 
	 Tnd[i] = i * dx_dr * Tnd[i-1];

  Tnd[0] = 0.0;

}

inline void Cheby::set_3b_powers(const TRIPLETS & FF_3BODY, const vector<int> &pair_index, int POWER_SET,
											int & pow_ij, int & pow_ik, int & pow_jk ) 
// Matches the allowed powers to the ij. ik, jk type pairs formed from the atom triplet ai, aj, ak 
// given the index transformation given in pair_index.
{
  pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET][pair_index[0]];
  pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET][pair_index[1]];
  pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET][pair_index[2]];
}

//////////////////////////////////////////
// 4B Cheby functions
//////////////////////////////////////////


void Cheby::map_indices(CLUSTER & cluster, vector<string> & atom_type, vector<int> & pair_map) 
// Matches the allowed powers to the ij, ik, il... type pairs formed from the atom sixlet ai, aj, ak, al
// This is based on reordering atoms, not pairs.  Pairs can't be interchanged independently, but atoms can.
// "pair_index" serves as the map between the powers in the cluster and the powers between the atoms.
//
// Not currently used.  map_indices_int is used instead.
{
  using PAIR_STR_INT = pair<string,int> ;
	
  int npairs = cluster.NPAIRS ;
  int natoms = cluster.NATOMS ;

#ifdef DEBUG_CHEBY
  if ( pair_map.size() != npairs )
	 EXIT_MSG("Wrong dimension for pair map") ;

  if ( atom_type.size() != natoms ) 
	 EXIT_MSG("Wrong dimension for natoms") ;
#endif

  bool called_before = false ;
  static vector <PAIR_STR_INT> 	ATOM_TYPE_AND_INDEX;
  static vector<vector<int>> pair_index ;
  static vector<int> rev_index ;

  if ( ! called_before ) 
  {
	 called_before = true ;
	 ATOM_TYPE_AND_INDEX.resize(natoms) ;
	 pair_index.resize(natoms) ;
	 rev_index.resize(natoms) ;
	 for ( int j = 0 ; j < natoms ; j++ ) 
	 {
		pair_index[j].resize(natoms) ;
	 }
	  
	 int count = 0 ;
	 for ( int j = 0 ; j < natoms ; j++ ) 
	 {
		pair_index[j][j] = -1 ;
		for ( int k = j + 1 ; k < natoms ; k++ ) 
		{
		  pair_index[j][k] = count ;
		  pair_index[k][j] = count ;
		  count++ ;
		}
	 }
  }

  for(int m = 0 ; m< natoms ; m++)
  { 				
	 ATOM_TYPE_AND_INDEX[m].first  = atom_type[m];					// 
	 ATOM_TYPE_AND_INDEX[m].second = m;	
  }

  // Code currently assumes ordering of atom names.
  for(int m = 0 ; m < natoms - 1 ; m++) {
	 if ( cluster.ATOM_NAMES[m] > cluster.ATOM_NAMES[m+1] )
		EXIT_MSG("ATOM NAMES are out of order") ;
  }

  sort (ATOM_TYPE_AND_INDEX.begin(), ATOM_TYPE_AND_INDEX.end());	// Sort the vector contents... automatically does on the basis of the .first element, preserving "link" between .first and .second

  for ( int j = 0 ; j < natoms ; j++ ) 
  {
	 for ( int i = 0 ; i < natoms ; i++ ) 
	 {
		if ( ATOM_TYPE_AND_INDEX[i].second == j ) 
		{
		  rev_index[j] = i ;
		  break ;
		}
	 }
  }

#ifdef DEBUG_CHEBY
  for ( int j = 0 ; j < natoms ; j++ ) 
  {
	 if ( rev_index[j] < 0 ) 
		EXIT_MSG("Bad reverse index") ;
  }
#endif
		 
  int m = 0 ;
  for(int i = 0 ; i < natoms ; i++)
  {
	 int ii = rev_index[i] ;
	 for ( int j = i + 1 ; j < natoms ; j++ ) 
	 {

		int jj = rev_index[j] ;
		pair_map[m] = pair_index[ii][jj] ;

#ifdef DEBUG_CHEBY
		if ( pair_map[m] < 0 || pair_map[m] > npairs - 1 ) 
		  EXIT_MSG("Permutation power map has a bad value") ;
#endif
		++m ;
	 }
  }

#ifdef DEBUG_CHEBY
  // Double check uniqueness of each value in pair_map.
  for ( int m = 0 ; m < npairs ; m++ ) 
  {
	 for ( int m1 = 0 ; m1 < npairs ; m1++ ) 
	 {
		if ( m1 == m ) continue ;
		if ( pair_map[m1] == pair_map[m] )
		  EXIT_MSG("Permutation power map had a repeated value" ) ;
	 }
  }
#endif
}


void Cheby::map_indices_int(CLUSTER & cluster, vector<int> & atom_type_idx, vector<int> & pair_map) 
// Matches the allowed powers to the ij, ik, il... type pairs formed from the atoms  ai, aj, ak, al
// This is based on reordering atoms, not pairs.  Pairs can't be interchanged independently, but atoms can.
// "pair_index" serves as the map between the powers in the cluster and the powers between the atoms.
{
  using Pair_int_int = pair<int,int> ;

  static vector <Pair_int_int> 	ATOM_TYPE_AND_INDEX;	// [a1/a2/a3/a4 atom pair chemistry][index of pair (from 1-6)]
  //vector <Pair_int_int> 	sorted_index(natoms);	// Sorted to match atom_type_idx.
  static vector<bool> used ;
  static vector<int> rev_index ;
  static vector<vector<int>> pair_index ;
  static bool called_before = false ;

  int npairs = cluster.NPAIRS ;
  int natoms = cluster.NATOMS ;

  if ( ! called_before )
  {
	 called_before = true ;
	 pair_index.resize(natoms) ;
	 used.resize(natoms) ;
	 rev_index.resize(natoms) ;
	 ATOM_TYPE_AND_INDEX.resize(natoms) ;

	 for ( int j = 0 ; j < natoms ; j++ ) 
	 {
		pair_index[j].resize(natoms) ;
	 }
  }
	
#ifdef DEBUG_CHEBY
  if ( pair_map.size() != npairs )
	 EXIT_MSG("Wrong dimension for pair map") ;

  if ( atom_type_idx.size() != natoms ) 
	 EXIT_MSG("Wrong dimension for natoms") ;
#endif
  vector<int> &indices = cluster.ATOM_INDICES ;
	  
  int count = 0 ;
  for ( int j = 0 ; j < natoms ; j++ ) 
  {
	 pair_index[j][j] = -1 ;
	 for ( int k = j + 1 ; k < natoms ; k++ ) 
	 {
		pair_index[j][k] = count ;
		pair_index[k][j] = count ;
		count++ ;
	 }
  }

  for(int m = 0 ; m< natoms ; m++)
  { 				
	 ATOM_TYPE_AND_INDEX[m].first  = atom_type_idx[m];					// 
	 ATOM_TYPE_AND_INDEX[m].second = m;	
  }

  // Match the permutation of the atom indices in atom_type_idx to
  // those in the cluster indices[]
  fill(used.begin(), used.end(), false) ;
#ifdef DEBUG_CHEBY
  fill(rev_index.begin(), rev_index.end(), -1) ;
#endif

  for ( int m = 0 ; m < natoms ; m++ )
	 for ( int n = 0 ; n < natoms ; n++ ) 
	 {
		if ( ATOM_TYPE_AND_INDEX[m].first == indices[n] && ! used[n] ) 
		{
		  //sorted_index[n] = ATOM_TYPE_AND_INDEX[m] ;
		  rev_index[ATOM_TYPE_AND_INDEX[m].second] = n ;
		  used[n] = true ;
		  break ;
		}
	 }

  // Find the reverse of the sorted index.
  // vector<int> rev_index(natoms, -1) ;
  // for ( int j = 0 ; j < natoms ; j++ ) 
  // {
  // 	 for ( int i = 0 ; i < natoms ; i++ ) 
  // 	 {
  // 		if ( sorted_index[i].second == j ) 
  // 		{
  // 		  rev_index[j] = i ;
  // 		  break ;
  // 		}
  // 	 }
  // }

#ifdef DEBUG_CHEBY
  for ( int j = 0 ; j < natoms ; j++ ) 
  {
	 if ( rev_index[j] < 0 ) 
		EXIT_MSG("Bad reverse index") ;
  }
#endif
		 
  int m = 0 ;
  for(int i = 0 ; i < natoms ; i++)
  {
	 int ii = rev_index[i] ;
	 for ( int j = i + 1 ; j < natoms ; j++ ) 
	 {

		int jj = rev_index[j] ;
		pair_map[m] = pair_index[ii][jj] ;

#ifdef DEBUG_CHEBY
		if ( pair_map[m] < 0 || pair_map[m] > npairs - 1 ) 
		  EXIT_MSG("Permutation power map has a bad value") ;
#endif
		++m ;
	 }
  }

  // Double check uniqueness of each value in pair_map.
#ifdef DEBUG_CHEBY
  for ( int m = 0 ; m < npairs ; m++ ) 
  {
	 for ( int m1 = 0 ; m1 < npairs ; m1++ ) 
	 {
		if ( m1 == m ) continue ;
		if ( pair_map[m1] == pair_map[m] )
		  EXIT_MSG("Permutation power map had a repeated value" ) ;
	 }
  }
#endif
}

void Cheby::Deriv_2B(vector<vector <XYZ > > & FRAME_A_MATRIX)
 // Calculate derivatives of the forces wrt the Chebyshev parameters. Stores minimum distance between a pair of atoms in minD[i].
 {
	 XYZ RAB; 		// Replaces  Rab[3];
	 double rlen;
	 int vstart;
	 static double *Tn, *Tnd;
	 static bool called_before = false;

	 double fcut; 
	 double fcutderiv; 				
	 double deriv;
	 double tmp_doub; 	

	 double inv_vol = 1.0 / (SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z) ;

	 if ( ! called_before ) 
	 {
		 called_before = true;
		 int dim = 0;

		 for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			 if (FF_2BODY[i].SNUM > dim ) 
				 dim = FF_2BODY[i].SNUM;	 

		 dim++;
		 Tn   = new double [dim];
		 Tnd  = new double [dim];

	 }

	 // Main loop for Chebyshev terms:

	 string TEMP_STR;
	 int curr_pair_type_idx;

	 // Set up for layering

	 int fidx_a2;
	 int a2start, a2end, a2;

	 int MATR_SIZE = FRAME_A_MATRIX.size();

	 for(int a1=0;a1<SYSTEM.ATOMS;a1++)		// Double sum over atom pairs
	 {
		 a2start = 0;
		 a2end   = NEIGHBOR_LIST.LIST[a1].size();

		 for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		 {			
			 a2 = NEIGHBOR_LIST.LIST[a1][a2idx];		

			 curr_pair_type_idx = get_pair_index(a1, a2, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
																	  SYSTEM.PARENT) ;

			 //calculate vstart: (index for populating OO, OH, or HH column block of A).

			 vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			 // Get pair distance

			 rlen = get_dist(SYSTEM, RAB, a1, a2);	// Updates RAB!

			 if ( (rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST))	
				 FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

			 if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM and rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
			 {
				 FF_2BODY[curr_pair_type_idx].N_CFG_CONTRIB++;

				 // Do the distance transformation
				 double x_diff = FF_2BODY[curr_pair_type_idx].X_DIFF ;
				 double x_avg  = FF_2BODY[curr_pair_type_idx].X_AVG ;
				 set_polys(curr_pair_type_idx, Tn, Tnd, rlen, x_diff, x_avg, FF_2BODY[curr_pair_type_idx].SNUM ) ;

				 // fcut and fcutderv are the cutoff functions (1-r/rcut)**3 and its
				 // derivative -3 (1-r/rcut)**2/rcut.  This ensures that
				 // the force goes to 0 as r goes to rcut.
				 // This is not a penalty function in the usual sense.  
				 // I don't see any reason to have a scaling on the cutoff.

				 // That will simply multiply all the forces by a constant,
				 // which will be canceled out during the force matching process.
				 // (LEF)

				 // fcut and fcutderv are the form that the penalty func and its derivative for the morse-type pair distance transformation

				 FF_2BODY[curr_pair_type_idx].FORCE_CUTOFF.get_fcut(fcut, fcutderiv, rlen, 0, 
																					 FF_2BODY[curr_pair_type_idx].S_MAXIM);

				 // Compute part of the derivative
				 // NOTE: All these extra terms are coming from:
				 //
				 // 1. Chain rule to account for transformation from morse-type pair distance to x
				 // 2. Product rule coming from pair distance dependence of fcut, the penalty function

				 fidx_a2 = SYSTEM.PARENT[a2];

				 for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
				 {
					 tmp_doub = (fcut * Tnd[i+1] + fcutderiv * Tn[i+1] );

					 // Finally, account for the x, y, and z unit vectors

					 deriv = tmp_doub * RAB.X / rlen;
					 FRAME_A_MATRIX[a1     ][vstart+i].X += deriv;
					 FRAME_A_MATRIX[fidx_a2][vstart+i].X -= deriv;

					 deriv = tmp_doub * RAB.Y / rlen; 
					 FRAME_A_MATRIX[a1     ][vstart+i].Y += deriv;
					 FRAME_A_MATRIX[fidx_a2][vstart+i].Y -= deriv;

					 deriv = tmp_doub * RAB.Z / rlen;
					 FRAME_A_MATRIX[a1     ][vstart+i].Z += deriv;
					 FRAME_A_MATRIX[fidx_a2][vstart+i].Z -= deriv;

					 if (CONTROLS.FIT_STRESS)
					 {
						 FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].X -= tmp_doub * RAB.X * RAB.X / rlen;
						 FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Y -= tmp_doub * RAB.Y * RAB.Y / rlen;
						 FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Z -= tmp_doub * RAB.Z * RAB.Z / rlen;	
					 }

					 else if (CONTROLS.FIT_STRESS_ALL)
					 {
						 FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].X -= tmp_doub * RAB.X * RAB.X / rlen;	// xx
						 FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Y -= tmp_doub * RAB.X * RAB.Y / rlen;	// xy
						 FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Z -= tmp_doub * RAB.X * RAB.Z / rlen;	// xz

						 FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].X -= tmp_doub * RAB.Y * RAB.X / rlen;	// yx
						 FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Y -= tmp_doub * RAB.Y * RAB.Y / rlen;	// yy
						 FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Z -= tmp_doub * RAB.Y * RAB.Z / rlen;	// yz

						 FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].X -= tmp_doub * RAB.Z * RAB.X / rlen;	// zx
						 FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Y -= tmp_doub * RAB.Z * RAB.Y / rlen;	// zy
						 FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Z -= tmp_doub * RAB.Z * RAB.Z / rlen;	// zz
					 }

					 if(CONTROLS.FIT_ENER) 
					 {
						 FRAME_A_MATRIX[MATR_SIZE-1][vstart+i].X += fcut * Tn[i+1];
						 FRAME_A_MATRIX[MATR_SIZE-1][vstart+i].Y += fcut * Tn[i+1];
						 FRAME_A_MATRIX[MATR_SIZE-1][vstart+i].Z += fcut * Tn[i+1];
					 }
				 }
			 }

		 }
	 }

	 if (CONTROLS.FIT_STRESS)
	 {


		 for ( int i = 0; i < CONTROLS.TOT_SNUM; i++ ) 
		 {
			 FRAME_A_MATRIX[SYSTEM.ATOMS][i].X *= inv_vol;
			 FRAME_A_MATRIX[SYSTEM.ATOMS][i].Y *= inv_vol;
			 FRAME_A_MATRIX[SYSTEM.ATOMS][i].Z *= inv_vol;	 
		 }
	 }
	 else if (CONTROLS.FIT_STRESS_ALL)
	 {		
		 for ( int i = 0; i < CONTROLS.TOT_SNUM; i++ ) 
				{
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][i].Z *= inv_vol;	 
			
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][i].Z *= inv_vol;
			
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][i].Z *= inv_vol;
		}
		
	}

  return;
}

void Cheby::Deriv_3B(vector<vector <XYZ > > & FRAME_A_MATRIX, CLUSTER_LIST &TRIPS) 
							
// Calculate derivatives of the forces wrt the 3-body Chebyshev parameters. 
{
	// This three body interaction stems from: C_n^ij *  C_n^ik * C_n^jk * T_n(x_ij) * T_n(x_ik) * T_n(x_jk)
	//	The logic:
	//	+ Run a triple loop over all atoms in the system.
	//	+ Compute C_ij, C_ik, and C_jk coeffiecients independently as you would do for a normal 2 body 

	XYZ RAB_IJ;
	XYZ RAB_IK;
	XYZ RAB_JK; 		

	int ij_bin;
	int ik_bin;
	int jk_bin;
	
	double rlen_ij,  rlen_ik,  rlen_jk;
	double rlen_ij_dummy, rlen_ik_dummy, rlen_jk_dummy;
	int vstart;
    static int n_2b_cheby_terms, n_3b_cheby_terms;
	static double *Tn_ij,  *Tn_ik,  *Tn_jk;
	static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
	static bool called_before = false;
	
	static int pow_ij, pow_ik, pow_jk;

	double fcut_ij,  fcut_ik,  fcut_jk; 			
	double deriv_ij, deriv_ik, deriv_jk;
	double fcutderiv_ij, fcutderiv_ik, fcutderiv_jk; 	
	double force_wo_coeff_ij, force_wo_coeff_ik, force_wo_coeff_jk;
	
	string TEMP_STR;
	int curr_triple_type_index;
	int curr_pair_type_idx_ij;
	int curr_pair_type_idx_ik;
	int curr_pair_type_idx_jk;
	int row_offset;	
	
	double S_MAXIM_IJ, S_MAXIM_IK, S_MAXIM_JK;
	double S_MINIM_IJ, S_MINIM_IK, S_MINIM_JK;
	
	double inv_vol = 1.0 / (SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z) ;

	vector<CLUSTER> &PAIR_TRIPLETS = TRIPS.VEC ;

	bool FORCE_IS_ZERO_IJ, FORCE_IS_ZERO_IK, FORCE_IS_ZERO_JK;

	vector<int> pair_index(3) ;
	vector<double> x_diff(3), x_avg(3) ;
	vector<int> atom_type_index(3) ;

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;

		
		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
		{
			if (FF_2BODY[i].SNUM_3B_CHEBY > dim ) 
				dim = FF_2BODY[i].SNUM_3B_CHEBY;	
			
			n_2b_cheby_terms += FF_2BODY[i].SNUM;
		}
		for ( int i = 0; i < PAIR_TRIPLETS.size(); i++ ) 
			n_3b_cheby_terms += PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS;;
		
		dim++;
		
		Tn_ij   = new double [dim];
		Tn_ik   = new double [dim];
		Tn_jk   = new double [dim];

		Tnd_ij  = new double [dim];
		Tnd_ik  = new double [dim];
		Tnd_jk  = new double [dim];
		
	}

	// Set up for layering

	int fidx_a2, fidx_a3;
	
	// Set up for MPI
	
	int a1start, a1end;	

	//divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
	a1start = 0;
	a1end = SYSTEM.ATOMS-1;

	// Set up for neighbor lists
	
	int a2start, a2end, a2;
	int a3start, a3end, a3;
	
	int MATR_SIZE = FRAME_A_MATRIX.size();
	
	for(int a1=a1start; a1<=a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev -1)
	{
		a2start = 0;
	
		// Use a special neighbor list for 3 body interations.
		a2end   = NEIGHBOR_LIST.LIST_3B[a1].size();

		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = NEIGHBOR_LIST.LIST_3B[a1][a2idx];

			// Get a3 as a neighbor of a1 to avoid
			// creating neighbor lists for ghost atoms.
			
			a3start = 0;
			a3end   = NEIGHBOR_LIST.LIST_3B[a1].size();
			

			for(int a3idx=a3start; a3idx<a3end; a3idx++)	
			{			
				a3 = NEIGHBOR_LIST.LIST_3B[a1][a3idx];
				
				if ( a3 == a2 || SYSTEM.PARENT[a2] > SYSTEM.PARENT[a3] ) 
					continue;

				curr_pair_type_idx_ij =  get_pair_index(a1, a2, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,SYSTEM.PARENT) ;
				curr_pair_type_idx_ik =  get_pair_index(a1, a3, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,SYSTEM.PARENT) ;
				curr_pair_type_idx_jk =  get_pair_index(a2, a3, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,SYSTEM.PARENT) ;

				atom_type_index[0] = SYSTEM.get_atomtype_idx(a1) ;
				atom_type_index[1] = SYSTEM.get_atomtype_idx(a2) ;
				atom_type_index[2] = SYSTEM.get_atomtype_idx(a3) ;

				int tidx = TRIPS.make_id_int(atom_type_index) ;
				curr_triple_type_index = TRIPS.INT_MAP[tidx];
				
				// If this type has been excluded, then skip to the next iteration of the loop

				if(curr_triple_type_index<0) 
				{
				  //cout << "Interaction " << TEMP_STR << " is excluded\n" ;
				  continue;
				}

				for ( int j = 0 ; j < 3 ; j++ ) 
				  pair_index[j] = TRIPS.PAIR_INDICES[tidx][j] ;

				rlen_ij = get_dist(SYSTEM, RAB_IJ, a1, a2);	// Updates RAB!
				rlen_ik = get_dist(SYSTEM, RAB_IK, a1, a3);	// Updates RAB!
				rlen_jk = get_dist(SYSTEM, RAB_JK, a2, a3);	// Updates RAB!

				S_MAXIM_IJ = PAIR_TRIPLETS[curr_triple_type_index].S_MAXIM[pair_index[0]] ;
				S_MAXIM_IK = PAIR_TRIPLETS[curr_triple_type_index].S_MAXIM[pair_index[1]] ;
				S_MAXIM_JK = PAIR_TRIPLETS[curr_triple_type_index].S_MAXIM[pair_index[2]] ;
				
				S_MINIM_IJ = PAIR_TRIPLETS[curr_triple_type_index].S_MINIM[pair_index[0]] ;
				S_MINIM_IK = PAIR_TRIPLETS[curr_triple_type_index].S_MINIM[pair_index[1]] ;
				S_MINIM_JK = PAIR_TRIPLETS[curr_triple_type_index].S_MINIM[pair_index[2]] ;
				

				
				// Before doing any polynomial/coeff set up, make sure that all ij, ik, and jk distances are 
				// within the allowed range.
				// Unlike the 2-body Cheby, extrapolation/refitting to handle behavior outside of fitting regime is not straightforward.
				
				FORCE_IS_ZERO_IJ = FORCE_IS_ZERO_IK = FORCE_IS_ZERO_JK = false;		
/*				
				cout << "CURR_DISTANCES: " << rlen_ij << " " << rlen_ik << " " << rlen_jk << "		" 
				     << PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_ij, S_MINIM_IJ, S_MAXIM_IJ) << " " 
				     << PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_ik, S_MINIM_IK, S_MAXIM_IK) << " "
				     << PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_jk, S_MINIM_JK, S_MAXIM_JK) << endl;
*/				
				if( PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_ij, S_MINIM_IJ, S_MAXIM_IJ))
				{
					if( PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_ik, S_MINIM_IK, S_MAXIM_IK))
					{
						if( PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_jk, S_MINIM_JK, S_MAXIM_JK))
						{		

							// Populate the appropriate histogram
							
							if(PAIR_TRIPLETS[curr_triple_type_index].NBINS[0]>0)
							{
								
								ij_bin = -10;
								ik_bin = -10;
								jk_bin = -10;

								// Sync up ij ik jk with ATOMPAIR1 ATOMPAIR2  ATOMPAIR3 and NBINS[0] NBINS[1] NBINS[2]
							
								if(FF_2BODY[curr_pair_type_idx_ij].PRPR_NM == PAIR_TRIPLETS[curr_triple_type_index].ATOM_PAIRS[0])
								{
									ij_bin = int(ceil((rlen_ij-S_MINIM_IJ)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[0])-1.0);	
									
									if(FF_2BODY[curr_pair_type_idx_ik].PRPR_NM == PAIR_TRIPLETS[curr_triple_type_index].ATOM_PAIRS[1])
									{
										ik_bin = int(ceil((rlen_ik-S_MINIM_IK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[1])-1.0);
										jk_bin = int(ceil((rlen_jk-S_MINIM_JK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[2])-1.0);
									}
									else
									{
										ik_bin = int(ceil((rlen_ik-S_MINIM_IK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[2])-1.0);
										jk_bin = int(ceil((rlen_jk-S_MINIM_JK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[1])-1.0);
									}
								}
								
								if(FF_2BODY[curr_pair_type_idx_ij].PRPR_NM == PAIR_TRIPLETS[curr_triple_type_index].ATOM_PAIRS[1])
								{
									ij_bin = int(ceil((rlen_ij-S_MINIM_IJ)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[1])-1.0);	
									
									if(FF_2BODY[curr_pair_type_idx_ik].PRPR_NM == PAIR_TRIPLETS[curr_triple_type_index].ATOM_PAIRS[0])
									{
										ik_bin = int(ceil((rlen_ik-S_MINIM_IK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[0])-1.0);
										jk_bin = int(ceil((rlen_jk-S_MINIM_JK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[2])-1.0);
									}
									else
									{
										ik_bin = int(ceil((rlen_ik-S_MINIM_IK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[2])-1.0);
										jk_bin = int(ceil((rlen_jk-S_MINIM_JK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[0])-1.0);
									}
								}
								
								if(FF_2BODY[curr_pair_type_idx_ij].PRPR_NM == PAIR_TRIPLETS[curr_triple_type_index].ATOM_PAIRS[2])
								{
									ij_bin = int(ceil((rlen_ij-S_MINIM_IJ)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[2])-1.0);	
									
									if(FF_2BODY[curr_pair_type_idx_ik].PRPR_NM == PAIR_TRIPLETS[curr_triple_type_index].ATOM_PAIRS[0])
									{
										ik_bin = int(ceil((rlen_ik-S_MINIM_IK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[0])-1.0);
										jk_bin = int(ceil((rlen_jk-S_MINIM_JK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[1])-1.0);
									}
									else
									{
										ik_bin = int(ceil((rlen_ik-S_MINIM_IK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[1])-1.0);
										jk_bin = int(ceil((rlen_jk-S_MINIM_JK)/PAIR_TRIPLETS[curr_triple_type_index].BINWS[0])-1.0);
									}
								}								

								if(ij_bin == -1)
									ij_bin++;
								if(ik_bin == -1)
									ik_bin++;
								if(jk_bin == -1)
									jk_bin++;
								
								if(ij_bin<0 || ik_bin<0 || jk_bin<0)
								{
									cout << "ERROR: bad bin (1)" << endl;
									exit_run(0);
								}
								
								if(ij_bin>=PAIR_TRIPLETS[curr_triple_type_index].NBINS[0] ||
								   ik_bin>=PAIR_TRIPLETS[curr_triple_type_index].NBINS[1] ||
								   jk_bin>=PAIR_TRIPLETS[curr_triple_type_index].NBINS[2]   )
								{
									cout << "ERROR: bad bin (2)" << endl;
									exit_run(0);
								}
							
								vector<int> bin_vecs(3) ;
								bin_vecs[0] = ij_bin ;
								bin_vecs[1] = ik_bin ;
								bin_vecs[2] = jk_bin ;
								PAIR_TRIPLETS[curr_triple_type_index].increment_histogram(bin_vecs) ;
							}
							
						
							// For all types, if r < rcut, then the potential is constant, thus the force  must be zero.
							// Additional, the potential is then taken to be the potential at r_cut.
							
							rlen_ij_dummy = rlen_ij;
							rlen_ik_dummy = rlen_ik;
							rlen_jk_dummy = rlen_jk;
							
							if(rlen_ij < S_MINIM_IJ) 
							{
								rlen_ij_dummy = S_MINIM_IJ;
								FORCE_IS_ZERO_IJ = true;
							}
							
							if(rlen_ik < S_MINIM_IK)
							{
								rlen_ik_dummy = S_MINIM_IK;
								FORCE_IS_ZERO_IK = true;
							}
							
							if(rlen_jk < S_MINIM_JK)
							{
								rlen_jk_dummy = S_MINIM_JK;
								FORCE_IS_ZERO_JK = true;
							}
		
							// Everything is within allowed ranges.
							
							// Track the minimum triplet distances for each given pair
							
							if (PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[0] == -1) 	// Then this is our first check. Just set all equal to current distances
							{
								PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[0]] = rlen_ij;
								PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[1]] = rlen_ik;
								PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[2]] = rlen_jk;
							}
							
							// Case 2: If any distance is smaller than a previous distance
							
							else 
							{
								if (rlen_ij<PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[0]])
									PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[0]] = rlen_ij;
								
								if (rlen_ik<PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[1]])
									PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[1]] = rlen_ik;
								
								if (rlen_jk<PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[2]])
									PAIR_TRIPLETS[curr_triple_type_index].MIN_FOUND[pair_index[2]] = rlen_jk;
									
							}
				
							// Add this to the number of configs contributing to a fit for this triplet type
							
							PAIR_TRIPLETS[curr_triple_type_index].N_CFG_CONTRIB++;

							// Begin setting up the derivative calculation

							// Set up the polynomials
			
							for ( int jj = 0 ; jj < 3 ; jj++ ) 
							{
							  x_avg [jj] = PAIR_TRIPLETS[curr_triple_type_index].X_AVG [pair_index[jj]] ;
							  x_diff[jj] = PAIR_TRIPLETS[curr_triple_type_index].X_DIFF[pair_index[jj]] ;
							}							
							
							set_polys(curr_pair_type_idx_ij, Tn_ij, Tnd_ij, rlen_ij_dummy, x_diff[0], x_avg[0], FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY);			
							set_polys(curr_pair_type_idx_ik, Tn_ik, Tnd_ik, rlen_ik_dummy, x_diff[1], x_avg[1], FF_2BODY[curr_pair_type_idx_ik].SNUM_3B_CHEBY);			
							set_polys(curr_pair_type_idx_jk, Tn_jk, Tnd_jk, rlen_jk_dummy, x_diff[2], x_avg[2], FF_2BODY[curr_pair_type_idx_jk].SNUM_3B_CHEBY);			

							// At this point we've completed all pre-calculations needed to populate the A matrix. Now we need to figure out 
							// where within the matrix to put the data, and to do so. 

							// Note: This syntax is safe since there is only one possible SNUM_3B_CHEBY value for all interactions

							vstart = n_2b_cheby_terms;
			
							for (int i=0; i<curr_triple_type_index; i++)
								vstart += PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS;						
							
							PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_ij, fcutderiv_ij, rlen_ij, S_MINIM_IJ, S_MAXIM_IJ);
							PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_ik, fcutderiv_ik, rlen_ik, S_MINIM_IK, S_MAXIM_IK);
							PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_jk, fcutderiv_jk, rlen_jk, S_MINIM_JK, S_MAXIM_JK);	
/*							
							cout << PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.to_string() << " " << fcut_ij << " " << fcutderiv_ij << " " << rlen_ij << " " << S_MAXIM_IJ << endl;
							cout << PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.to_string() << " " << fcut_ik << " " << fcutderiv_ik << " " << rlen_ik << " " << S_MAXIM_IK << endl;
							cout << PAIR_TRIPLETS[curr_triple_type_index].FORCE_CUTOFF.to_string() << " " << fcut_jk << " " << fcutderiv_jk << " " << rlen_jk << " " << S_MAXIM_JK << endl;
							cout << endl;
*/							
		
							/////////////////////////////////////////////////////////////////////
							/////////////////////////////////////////////////////////////////////
							// Consider special restrictions on allowed triplet types and powers
							/////////////////////////////////////////////////////////////////////
							/////////////////////////////////////////////////////////////////////
			
							row_offset = 0;
			
							// --- THE KEY HERE IS TO UNDERSTAND THAT THE IJ, IK, AND JK HERE IS BASED ON ATOM PAIRS, AND DOESN'T NECESSARILY MATCH THE TRIPLET'S EXPECTED ORDER!
			
							
							fidx_a2 = SYSTEM.PARENT[a2];
							fidx_a3 = SYSTEM.PARENT[a3];

							vector<int> pair_idx(3) ;

							for(int i=0; i<PAIR_TRIPLETS[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
							{
							    row_offset = PAIR_TRIPLETS[curr_triple_type_index].PARAM_INDICES[i];
								
								 set_3b_powers(PAIR_TRIPLETS[curr_triple_type_index], pair_index, i,
													pow_ij, pow_ik, pow_jk) ;

							    deriv_ij =  fcut_ij * Tnd_ij[pow_ij] + fcutderiv_ij * Tn_ij[pow_ij];
							    deriv_ik =  fcut_ik * Tnd_ik[pow_ik] + fcutderiv_ik * Tn_ik[pow_ik];
							    deriv_jk =  fcut_jk * Tnd_jk[pow_jk] + fcutderiv_jk * Tn_jk[pow_jk];	
								
								if(FORCE_IS_ZERO_IJ)
									force_wo_coeff_ij = 0;
								else
									force_wo_coeff_ij = deriv_ij * fcut_ik * fcut_jk * Tn_ik[pow_ik] * Tn_jk[pow_jk];
								
								if(FORCE_IS_ZERO_IK)
									force_wo_coeff_ik = 0;
								else
									force_wo_coeff_ik = deriv_ik * fcut_ij * fcut_jk * Tn_ij[pow_ij] * Tn_jk[pow_jk];
								
								if(FORCE_IS_ZERO_JK)
									force_wo_coeff_jk = 0;
								else
									force_wo_coeff_jk = deriv_jk * fcut_ij * fcut_ik * Tn_ij[pow_ij] * Tn_ik[pow_ik];
						
							    // ij pairs

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].X += force_wo_coeff_ij * RAB_IJ.X / rlen_ij;
							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].X -= force_wo_coeff_ij * RAB_IJ.X / rlen_ij;

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].Y += force_wo_coeff_ij * RAB_IJ.Y / rlen_ij;
							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Y -= force_wo_coeff_ij * RAB_IJ.Y / rlen_ij;

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].Z += force_wo_coeff_ij * RAB_IJ.Z / rlen_ij;
							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Z -= force_wo_coeff_ij * RAB_IJ.Z / rlen_ij;	


							    // ik pairs

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].X += force_wo_coeff_ik * RAB_IK.X / rlen_ik;
							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].X -= force_wo_coeff_ik * RAB_IK.X / rlen_ik;

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].Y += force_wo_coeff_ik * RAB_IK.Y / rlen_ik;
							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Y -= force_wo_coeff_ik * RAB_IK.Y / rlen_ik;

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].Z += force_wo_coeff_ik * RAB_IK.Z / rlen_ik;
							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Z -= force_wo_coeff_ik * RAB_IK.Z / rlen_ik;

							    // jk pairs

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].X += force_wo_coeff_jk * RAB_JK.X / rlen_jk;
							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].X -= force_wo_coeff_jk * RAB_JK.X / rlen_jk;

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Y += force_wo_coeff_jk * RAB_JK.Y / rlen_jk;
							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Y -= force_wo_coeff_jk * RAB_JK.Y / rlen_jk;

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Z += force_wo_coeff_jk * RAB_JK.Z / rlen_jk;
							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Z -= force_wo_coeff_jk * RAB_JK.Z / rlen_jk;

								if (CONTROLS.FIT_STRESS)
								{
								    // ij pairs

								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].X -= force_wo_coeff_ij * RAB_IJ.X * RAB_IJ.X / rlen_ij;
								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].Y -= force_wo_coeff_ij * RAB_IJ.Y * RAB_IJ.Y / rlen_ij;
								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].Z -= force_wo_coeff_ij * RAB_IJ.Z * RAB_IJ.Z / rlen_ij;	


								    // ik pairs

								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].X -= force_wo_coeff_ik * RAB_IK.X * RAB_IK.X / rlen_ik;
								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].Y -= force_wo_coeff_ik * RAB_IK.Y * RAB_IK.Y / rlen_ik;
								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].Z -= force_wo_coeff_ik * RAB_IK.Z * RAB_IK.Z / rlen_ik;

								    // jk pairs

								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].X -= force_wo_coeff_jk * RAB_JK.X * RAB_JK.X / rlen_jk;
								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].Y -= force_wo_coeff_jk * RAB_JK.Y * RAB_JK.Y / rlen_jk;
								    FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].Z -= force_wo_coeff_jk * RAB_JK.Z * RAB_JK.Z / rlen_jk;									
								}
								
								else if (CONTROLS.FIT_STRESS_ALL)
								{
								    // ij pairs: 

								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].X -= force_wo_coeff_ij * RAB_IJ.X * RAB_IJ.X / rlen_ij;
								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].Y -= force_wo_coeff_ij * RAB_IJ.X * RAB_IJ.Y / rlen_ij;
								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].Z -= force_wo_coeff_ij * RAB_IJ.X * RAB_IJ.Z / rlen_ij;	
									
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].X -= force_wo_coeff_ij * RAB_IJ.Y * RAB_IJ.X / rlen_ij;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].Y -= force_wo_coeff_ij * RAB_IJ.Y * RAB_IJ.Y / rlen_ij;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].Z -= force_wo_coeff_ij * RAB_IJ.Y * RAB_IJ.Z / rlen_ij;	
									
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].X -= force_wo_coeff_ij * RAB_IJ.Z * RAB_IJ.X / rlen_ij;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].Y -= force_wo_coeff_ij * RAB_IJ.Z * RAB_IJ.Y / rlen_ij;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].Z -= force_wo_coeff_ij * RAB_IJ.Z * RAB_IJ.Z / rlen_ij;
									
								    // ik pairs

								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].X -= force_wo_coeff_ik * RAB_IK.X * RAB_IK.X / rlen_ik;
								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].Y -= force_wo_coeff_ik * RAB_IK.X * RAB_IK.Y / rlen_ik;
								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].Z -= force_wo_coeff_ik * RAB_IK.X * RAB_IK.Z / rlen_ik;
									
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].X -= force_wo_coeff_ik * RAB_IK.Y * RAB_IK.X / rlen_ik;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].Y -= force_wo_coeff_ik * RAB_IK.Y * RAB_IK.Y / rlen_ik;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].Z -= force_wo_coeff_ik * RAB_IK.Y * RAB_IK.Z / rlen_ik;
									
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].X -= force_wo_coeff_ik * RAB_IK.Z * RAB_IK.X / rlen_ik;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].Y -= force_wo_coeff_ik * RAB_IK.Z * RAB_IK.Y / rlen_ik;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].Z -= force_wo_coeff_ik * RAB_IK.Z * RAB_IK.Z / rlen_ik;	
									
								    // jk pairs

								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].X -= force_wo_coeff_jk * RAB_JK.X * RAB_JK.X / rlen_jk;
								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].Y -= force_wo_coeff_jk * RAB_JK.X * RAB_JK.Y / rlen_jk;
								    FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].Z -= force_wo_coeff_jk * RAB_JK.X * RAB_JK.Z / rlen_jk;			
									
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].X -= force_wo_coeff_jk * RAB_JK.Y * RAB_JK.X / rlen_jk;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].Y -= force_wo_coeff_jk * RAB_JK.Y * RAB_JK.Y / rlen_jk;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].Z -= force_wo_coeff_jk * RAB_JK.Y * RAB_JK.Z / rlen_jk;		
									
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].X -= force_wo_coeff_jk * RAB_JK.Z * RAB_JK.X / rlen_jk;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].Y -= force_wo_coeff_jk * RAB_JK.Z * RAB_JK.Y / rlen_jk;
								    FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].Z -= force_wo_coeff_jk * RAB_JK.Z * RAB_JK.Z / rlen_jk;																																				
									
										
								}
								
									

								if(CONTROLS.FIT_ENER) 
								{
									FRAME_A_MATRIX[MATR_SIZE-1][vstart+row_offset].X += fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk];
									FRAME_A_MATRIX[MATR_SIZE-1][vstart+row_offset].Y += fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk];
									FRAME_A_MATRIX[MATR_SIZE-1][vstart+row_offset].Z += fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk];
								}
								
								
							}
						} // end if rlen_jk within cutoffs...
					} // end if rlen_ik within cutoffs...	
				} // end third loop over atoms							
			}
		}	
	}

	
	if (CONTROLS.FIT_STRESS)
	{
		for(int i=0; i<CONTROLS.NUM_3B_CHEBY; i++) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS][n_2b_cheby_terms+i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS][n_2b_cheby_terms+i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS][n_2b_cheby_terms+i].Z *= inv_vol;	
		}
	}
	
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		for(int i=0; i<CONTROLS.NUM_3B_CHEBY; i++) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][n_2b_cheby_terms+i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][n_2b_cheby_terms+i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][n_2b_cheby_terms+i].Z *= inv_vol;	
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][n_2b_cheby_terms+i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][n_2b_cheby_terms+i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][n_2b_cheby_terms+i].Z *= inv_vol;
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][n_2b_cheby_terms+i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][n_2b_cheby_terms+i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][n_2b_cheby_terms+i].Z *= inv_vol;
		}
	}	
}

void Cheby::Deriv_4B(vector<vector <XYZ > > & FRAME_A_MATRIX, int n_3b_cheby_terms, CLUSTER_LIST& QUADS)		
// Calculate derivatives of the forces wrt the 3-body Chebyshev parameters. 
{
	// BECKY: DON'T FORGET TO UPDATE THIS DESCRIPTION FOR 4-BODY INTERACTIONS
	// Note: 6-element vectors contain data in the form: ij, ik, il, jk, jl, kl
	
	// This three body interaction stems from: C_n^ij *  C_n^ik * C_n^jk * T_n(x_ij) * T_n(x_ik) * T_n(x_jk)
	//	The logic:
	//	+ Run a triple loop over all atoms in the system.
	//	+ Compute C_ij, C_ik, and C_jk coeffiecients independently as you would do for a normal 2 body 

	vector<XYZ> RVEC(6);	// Replaces RVEC_IJ, RVEC_IK...
	vector<XYZ> RAB (6);	// Replaces RAB_IJ, RAB_IK...
	
	vector<double> rlen(6);			// Replaces rlen_ij, rlen_ik...
	vector<double> rlen_dummy(6);	// Replaces rlen_ij_dummy, rlen_ik_dummy

	int vstart;
	static int n_2b_cheby_terms, n_4b_cheby_terms;
	static double *Tn_ij,  *Tn_ik,  *Tn_il,  *Tn_jk,  *Tn_jl,  *Tn_kl;
	static double *Tnd_ij, *Tnd_ik, *Tnd_il, *Tnd_jk, *Tnd_jl, *Tnd_kl;
	static bool called_before = false;
	
	static vector<int> powers(6);	// replaces pow_ij, pow_ik, pow_jk;
	vector<double> fcut0(6);		// replaces fcut0_ij, fcut0_ik, fcut0_jk; 
	vector<double> fcut(6);			// replaces cut_ij,  fcut_ik,  fcut_jk;
	vector<double> fcut_deriv(6);	// replaces fcutderiv_ij, fcutderiv_ik, fcutderiv_jk; 
	vector<double> deriv(6);		// replaces deriv_ij, deriv_ik, deriv_jk;
	vector<double> force_wo_coeff(6);	// replaces force_wo_coeff_ij, force_wo_coeff_ik, force_wo_coeff_jk;

	
	static string TEMP_STR;
	vector<int> atom_type_index(4);		// Index of type of atoms in the quad cluster.
	static int  curr_quad_type_index;
	vector<int> curr_pair_type_idx(6);	// replaces curr_pair_type_idx_ij, etc
	static int row_offset;	
	
	vector<double> S_MAXIM(6);	// replaces S_MAXIM_IJ, S_MAXIM_IK, S_MAXIM_JK;
	vector<double> S_MINIM(6);	// replaces S_MINIM_IJ, S_MINIM_IK, S_MINIM_JK;
	vector<double> x_avg(6) ;
	vector<double> x_diff(6) ;

	double inv_vol = 1.0 / ( SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z ) ;
	
	vector<bool> FORCE_IS_ZERO(6);	// replaces FORCE_IS_ZERO_IJ, FORCE_IS_ZERO_IK, FORCE_IS_ZERO_JK;
	
	double TMP_ENER;
	
	int ATOM_QUAD_ID_INT;
	vector<int>TMP_QUAD_SET(4);
	vector<int> pow_map(6);

	vector<QUADRUPLETS>& PAIR_QUADRUPLETS = QUADS.VEC ;

	if (!called_before) 
	{
		called_before = true;
		int dim = 0;
		n_2b_cheby_terms = 0;
		n_4b_cheby_terms = 0;
		

		
		for (int i=0; i<FF_2BODY.size(); i++) 
		{
			if (FF_2BODY[i].SNUM_4B_CHEBY > dim ) 
				dim = FF_2BODY[i].SNUM_4B_CHEBY;	
			
			n_2b_cheby_terms += FF_2BODY[i].SNUM;
		}
		
		for (int i=0; i<PAIR_QUADRUPLETS.size(); i++) 
			n_4b_cheby_terms += PAIR_QUADRUPLETS[i].N_TRUE_ALLOWED_POWERS;
		
		dim++;
		
		Tn_ij   = new double [dim];
		Tn_ik   = new double [dim];
		Tn_il   = new double [dim];
		Tn_jk   = new double [dim];
		Tn_jl   = new double [dim];
		Tn_kl   = new double [dim];

		Tnd_ij  = new double [dim];
		Tnd_ik  = new double [dim];
		Tnd_il  = new double [dim];
		Tnd_jk  = new double [dim];
		Tnd_jl  = new double [dim];
		Tnd_kl  = new double [dim];
		
	}

	// Set up for layering

	int fidx_a2, fidx_a3, fidx_a4;
	
	// Set up for MPI
	
	int a1start, a1end;	

	//divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
	a1start = 0;
	a1end = SYSTEM.ATOMS-1;

	// Set up for neighbor lists
	
	int a2start, a2end, a2;
	int a3start, a3end, a3;
	int a4start, a4end, a4;	
	
	int MATR_SIZE = FRAME_A_MATRIX.size();
	
	for(int a1=a1start; a1<=a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev -1)
	{
		a2start = 0;
		a2end   = NEIGHBOR_LIST.LIST_4B[a1].size();	// Borrow the special neighbor list for 3 body interations.

		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = NEIGHBOR_LIST.LIST_4B[a1][a2idx];

			// Get a3 as a neighbor of a1 to avoid creating neighbor lists for ghost atoms, do the same for a4
			
			a3start = 0;
			a3end   = NEIGHBOR_LIST.LIST_4B[a1].size();

			for(int a3idx=a3start; a3idx<a3end; a3idx++)	
			{			
				a3 = NEIGHBOR_LIST.LIST_4B[a1][a3idx];
				
				if ( a3 == a2 || SYSTEM.PARENT[a2] > SYSTEM.PARENT[a3] ) 
					continue;
				
				a4start = 0;
				a4end   = NEIGHBOR_LIST.LIST_4B[a1].size();
				
				for(int a4idx=a4start; a4idx<a4end; a4idx++)	
				{			
					a4 = NEIGHBOR_LIST.LIST_4B[a1][a4idx];
				
					// Already checked a2 == a3 and SYSTEM.PARENT[a2] > SYSTEM.PARENT[a3] in a3idx loop.
					if ( a2 == a4  || a3 == a4 || SYSTEM.PARENT[a3] > SYSTEM.PARENT[a4] ) 
						continue;
				
					// Determine the pair types and the triplet type
	
					curr_pair_type_idx[0] = get_pair_index(a1, a2, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP, SYSTEM.PARENT) ;			
					curr_pair_type_idx[1] = get_pair_index(a1, a3, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP, SYSTEM.PARENT) ;
					curr_pair_type_idx[2] = get_pair_index(a1, a4, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP, SYSTEM.PARENT) ;
					curr_pair_type_idx[3] = get_pair_index(a2, a3, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP, SYSTEM.PARENT) ;
					curr_pair_type_idx[4] = get_pair_index(a2, a4, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP, SYSTEM.PARENT) ;
					curr_pair_type_idx[5] = get_pair_index(a3, a4, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP, SYSTEM.PARENT) ;

					fidx_a2 = SYSTEM.PARENT[a2];
					fidx_a3 = SYSTEM.PARENT[a3];
					fidx_a4 = SYSTEM.PARENT[a4];
				
					atom_type_index[0] = SYSTEM.get_atomtype_idx(a1) ;
					atom_type_index[1] = SYSTEM.get_atomtype_idx(a2) ;
					atom_type_index[2] = SYSTEM.get_atomtype_idx(a3) ;
					atom_type_index[3] = SYSTEM.get_atomtype_idx(a4) ;
			
					// Always construct ATOM_QUAD_ID_INT lookup key based on atom types in decending order
			
					//sort   (atom_type_index.begin(), atom_type_index.end());
					//reverse(atom_type_index.begin(), atom_type_index.end());
			
					ATOM_QUAD_ID_INT       = QUADS.make_id_int(atom_type_index) ;

					curr_quad_type_index = QUADS.INT_MAP[ATOM_QUAD_ID_INT];
					
					// If this type has been excluded, then skip to the next iteration of the loop
					if(curr_quad_type_index<0)
						continue;
					
					// Get the atom distances

					rlen[0] = get_dist(SYSTEM, RAB[0], a1, a2);	// Updates RAB!
					rlen[1] = get_dist(SYSTEM, RAB[1], a1, a3);	// Updates RAB!
					rlen[2] = get_dist(SYSTEM, RAB[2], a1, a4);	// Updates RAB!
					rlen[3] = get_dist(SYSTEM, RAB[3], a2, a3);	// Updates RAB!
					rlen[4] = get_dist(SYSTEM, RAB[4], a2, a4);	// Updates RAB!
					rlen[5] = get_dist(SYSTEM, RAB[5], a3, a4);	// Updates RAB!
					
					// Determine the inner and outer cutoffs for each pair type in the quadruplet

					//SET_4B_CHEBY_POWERS(PAIR_QUADRUPLETS[curr_quad_type_index],ATOM_TYPE, pow_map);					
					// map_indices(PAIR_QUADRUPLETS[curr_quad_type_index],ATOM_TYPE, pow_map);					

					// map_indices_int(PAIR_QUADRUPLETS[curr_quad_type_index],atom_type_idx, pow_map);					
					for (int f=0; f<6; f++)
					{
					  pow_map[f] = QUADS.PAIR_INDICES[ATOM_QUAD_ID_INT][f] ;
					  S_MAXIM[f] = PAIR_QUADRUPLETS[curr_quad_type_index].S_MAXIM[pow_map[f]] ;
					  S_MINIM[f] = PAIR_QUADRUPLETS[curr_quad_type_index].S_MINIM[pow_map[f]] ;
					  x_diff [f] = PAIR_QUADRUPLETS[curr_quad_type_index].X_DIFF [pow_map[f]] ;
					  x_avg  [f] = PAIR_QUADRUPLETS[curr_quad_type_index].X_AVG  [pow_map[f]] ;
					}
						
					// Before doing any polynomial/coeff set up, make sure that all ij, ik, and jk distances are within the allowed range.
					// Unlike the 2-body Cheby, extrapolation/refitting to handle behavior outside of fitting regime is not straightforward.
					
					if( !PAIR_QUADRUPLETS[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[0], S_MINIM[0], S_MAXIM[0]))
						continue;
					if( !PAIR_QUADRUPLETS[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[1], S_MINIM[1], S_MAXIM[1]))
						continue;
					if( !PAIR_QUADRUPLETS[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[2], S_MINIM[2], S_MAXIM[2]))
						continue;
					if( !PAIR_QUADRUPLETS[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[3], S_MINIM[3], S_MAXIM[3]))
						continue;
					if( !PAIR_QUADRUPLETS[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[4], S_MINIM[4], S_MAXIM[4]))
						continue;
					if( !PAIR_QUADRUPLETS[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[5], S_MINIM[5], S_MAXIM[5]))
						continue;			
					
					// At this point, all distances are within allowed ranges. We can now proceed to the force derivative calculation
					
					// For all types, if r < rcut, then the potential is constant, thus the force  must be zero.
					// Additionally, the potential is then taken to be the potential at r_cut.
					
					for (int f=0; f<6; f++)
					{
						if(rlen[f] < S_MINIM[f])
						{
							rlen_dummy[f] = S_MINIM[f];
							FORCE_IS_ZERO[f] = true;
						} else {
						  rlen_dummy[f] = rlen[f];
						  FORCE_IS_ZERO[f] = false;
						}
					}

					// Track the minimum quadruplet distances for each given pair
					
					if (PAIR_QUADRUPLETS[curr_quad_type_index].MIN_FOUND[0] == -1) 	// Then this is our first check. Just set all equal to current distances
					{
						for (int f=0; f<6; f++)
							PAIR_QUADRUPLETS[curr_quad_type_index].MIN_FOUND[pow_map[f]] = rlen[f];
					}

					else // Case 2: If any distance is smaller than a previous distance
					{
						for (int f=0; f<6; f++)
						{
							if (rlen[f]<PAIR_QUADRUPLETS[curr_quad_type_index].MIN_FOUND[pow_map[f]])
								PAIR_QUADRUPLETS[curr_quad_type_index].MIN_FOUND[pow_map[f]] = rlen[f];
						}
					}
		
					// Add this to the number of configs contributing to a fit for this triplet type
					
					PAIR_QUADRUPLETS[curr_quad_type_index].N_CFG_CONTRIB++;

					// Begin setting up the derivative calculation

					// Set up the polynomials
	
					set_polys(curr_pair_type_idx[0], Tn_ij, Tnd_ij, rlen_dummy[0], x_diff[0], x_avg[0], FF_2BODY[curr_pair_type_idx[0]].SNUM_4B_CHEBY);
					set_polys(curr_pair_type_idx[1], Tn_ik, Tnd_ik, rlen_dummy[1], x_diff[1], x_avg[1], FF_2BODY[curr_pair_type_idx[1]].SNUM_4B_CHEBY);
					set_polys(curr_pair_type_idx[2], Tn_il, Tnd_il, rlen_dummy[2], x_diff[2], x_avg[2], FF_2BODY[curr_pair_type_idx[2]].SNUM_4B_CHEBY);
					set_polys(curr_pair_type_idx[3], Tn_jk, Tnd_jk, rlen_dummy[3], x_diff[3], x_avg[3], FF_2BODY[curr_pair_type_idx[3]].SNUM_4B_CHEBY);
					set_polys(curr_pair_type_idx[4], Tn_jl, Tnd_jl, rlen_dummy[4], x_diff[4], x_avg[4], FF_2BODY[curr_pair_type_idx[4]].SNUM_4B_CHEBY);
					set_polys(curr_pair_type_idx[5], Tn_kl, Tnd_kl, rlen_dummy[5], x_diff[5], x_avg[5], FF_2BODY[curr_pair_type_idx[5]].SNUM_4B_CHEBY);

					// At this point we've completed all pre-calculations needed to populate the A matrix. Now we need to figure out 
					// where within the matrix to put the data, and to do so. 

					// Note: This syntax is safe since there is only one possible SNUM_3B_CHEBY value for all interactions

					vstart = n_2b_cheby_terms + n_3b_cheby_terms;
	
					for (int i=0; i<curr_quad_type_index; i++)
						vstart += PAIR_QUADRUPLETS[i].N_TRUE_ALLOWED_POWERS;	

					for (int f=0; f<6; f++)
						PAIR_QUADRUPLETS[curr_quad_type_index].FORCE_CUTOFF.get_fcut(fcut[f], fcut_deriv[f], rlen[f], S_MINIM[f], S_MAXIM[f]);
	
					/////////////////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////////////////////
					// Consider special restrictions on allowed quadruplet types and powers
					/////////////////////////////////////////////////////////////////////
					/////////////////////////////////////////////////////////////////////

					row_offset = 0;
	
					// --- THE KEY HERE IS TO UNDERSTAND THAT THE IJ, IK, AND JK HERE IS BASED ON ATOM PAIRS, AND DOESN'T NECESSARILY MATCH THE QUAD'S EXPECTED ORDER!
	
					for(int i=0; i<PAIR_QUADRUPLETS[curr_quad_type_index].N_ALLOWED_POWERS; i++) 
					{
					    row_offset = PAIR_QUADRUPLETS[curr_quad_type_index].PARAM_INDICES[i];
						
						for (int f=0; f<6; f++)	
							powers[f] = PAIR_QUADRUPLETS[curr_quad_type_index].ALLOWED_POWERS[i][pow_map[f]];

					    	 deriv[0] =  fcut[0] * Tnd_ij[powers[0]] + fcut_deriv[0] * Tn_ij[powers[0]];
						 deriv[1] =  fcut[1] * Tnd_ik[powers[1]] + fcut_deriv[1] * Tn_ik[powers[1]];
						 deriv[2] =  fcut[2] * Tnd_il[powers[2]] + fcut_deriv[2] * Tn_il[powers[2]];
						 deriv[3] =  fcut[3] * Tnd_jk[powers[3]] + fcut_deriv[3] * Tn_jk[powers[3]];
						 deriv[4] =  fcut[4] * Tnd_jl[powers[4]] + fcut_deriv[4] * Tn_jl[powers[4]];
						 deriv[5] =  fcut[5] * Tnd_kl[powers[5]] + fcut_deriv[5] * Tn_kl[powers[5]];

						force_wo_coeff[0] = deriv[0] * fcut[1] * fcut[2] * fcut[3] * fcut[4] * fcut[5]  * Tn_ik[powers[1]]  * Tn_il[powers[2]]  * Tn_jk[powers[3]]  * Tn_jl[powers[4]]  * Tn_kl[powers[5]];
						force_wo_coeff[1] = deriv[1] * fcut[0] * fcut[2] * fcut[3] * fcut[4] * fcut[5]  * Tn_ij[powers[0]]  * Tn_il[powers[2]]  * Tn_jk[powers[3]]  * Tn_jl[powers[4]]  * Tn_kl[powers[5]];
						force_wo_coeff[2] = deriv[2] * fcut[0] * fcut[1] * fcut[3] * fcut[4] * fcut[5]  * Tn_ij[powers[0]]  * Tn_ik[powers[1]]  * Tn_jk[powers[3]]  * Tn_jl[powers[4]]  * Tn_kl[powers[5]];
						force_wo_coeff[3] = deriv[3] * fcut[0] * fcut[1] * fcut[2] * fcut[4] * fcut[5]  * Tn_ij[powers[0]]  * Tn_ik[powers[1]]  * Tn_il[powers[2]]  * Tn_jl[powers[4]]  * Tn_kl[powers[5]];
						force_wo_coeff[4] = deriv[4] * fcut[0] * fcut[1] * fcut[2] * fcut[3] * fcut[5]  * Tn_ij[powers[0]]  * Tn_ik[powers[1]]  * Tn_il[powers[2]]  * Tn_jk[powers[3]]  * Tn_kl[powers[5]];
						force_wo_coeff[5] = deriv[5] * fcut[0] * fcut[1] * fcut[2] * fcut[3] * fcut[4]  * Tn_ij[powers[0]]  * Tn_ik[powers[1]]  * Tn_il[powers[2]]  * Tn_jk[powers[3]]  * Tn_jl[powers[4]];
#if(0)
#endif
					    // ij pairs

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].X += force_wo_coeff[0] * RAB[0].X / rlen[0];
					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].X -= force_wo_coeff[0] * RAB[0].X / rlen[0];

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].Y += force_wo_coeff[0] * RAB[0].Y / rlen[0];
					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Y -= force_wo_coeff[0] * RAB[0].Y / rlen[0];

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].Z += force_wo_coeff[0] * RAB[0].Z / rlen[0];
					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Z -= force_wo_coeff[0] * RAB[0].Z / rlen[0];	


					    // ik pairs

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].X += force_wo_coeff[1] * RAB[1].X / rlen[1];
					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].X -= force_wo_coeff[1] * RAB[1].X / rlen[1];

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].Y += force_wo_coeff[1] * RAB[1].Y / rlen[1];
					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Y -= force_wo_coeff[1] * RAB[1].Y / rlen[1];

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].Z += force_wo_coeff[1] * RAB[1].Z / rlen[1];
					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Z -= force_wo_coeff[1] * RAB[1].Z / rlen[1];
						
					    // il pairs

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].X += force_wo_coeff[2] * RAB[2].X / rlen[2];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].X -= force_wo_coeff[2] * RAB[2].X / rlen[2];

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].Y += force_wo_coeff[2] * RAB[2].Y / rlen[2];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].Y -= force_wo_coeff[2] * RAB[2].Y / rlen[2];

					    FRAME_A_MATRIX[a1     ][vstart+row_offset].Z += force_wo_coeff[2] * RAB[2].Z / rlen[2];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].Z -= force_wo_coeff[2] * RAB[2].Z / rlen[2];

					    // jk pairs

					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].X += force_wo_coeff[3] * RAB[3].X / rlen[3];
					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].X -= force_wo_coeff[3] * RAB[3].X / rlen[3];

					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Y += force_wo_coeff[3] * RAB[3].Y / rlen[3];
					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Y -= force_wo_coeff[3] * RAB[3].Y / rlen[3];

					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Z += force_wo_coeff[3] * RAB[3].Z / rlen[3];
					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Z -= force_wo_coeff[3] * RAB[3].Z / rlen[3];
						
					    // jl pairs

					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].X += force_wo_coeff[4] * RAB[4].X / rlen[4];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].X -= force_wo_coeff[4] * RAB[4].X / rlen[4];

					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Y += force_wo_coeff[4] * RAB[4].Y / rlen[4];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].Y -= force_wo_coeff[4] * RAB[4].Y / rlen[4];

					    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Z += force_wo_coeff[4] * RAB[4].Z / rlen[4];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].Z -= force_wo_coeff[4] * RAB[4].Z / rlen[4];
						
					    // kl pairs

					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].X += force_wo_coeff[5] * RAB[5].X / rlen[5];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].X -= force_wo_coeff[5] * RAB[5].X / rlen[5];

					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Y += force_wo_coeff[5] * RAB[5].Y / rlen[5];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].Y -= force_wo_coeff[5] * RAB[5].Y / rlen[5];

					    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Z += force_wo_coeff[5] * RAB[5].Z / rlen[5];
					    FRAME_A_MATRIX[fidx_a4][vstart+row_offset].Z -= force_wo_coeff[5] * RAB[5].Z / rlen[5];

						if (CONTROLS.FIT_STRESS)
						{
							for (int f=0; f<6; f++)
							{
								FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].X -= force_wo_coeff[f] * RAB[f].X * RAB[f].X / rlen[f];
							    	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].Y -= force_wo_coeff[f] * RAB[f].Y * RAB[f].Y / rlen[f];
								FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+row_offset].Z -= force_wo_coeff[f] * RAB[f].Z * RAB[f].Z / rlen[f];	
							}								
						}
						
						else if (CONTROLS.FIT_STRESS_ALL)
						{
							for (int f=0; f<6; f++)
							{
								FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].X -= force_wo_coeff[f] * RAB[f].X * RAB[f].X / rlen[f];
							    	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].Y -= force_wo_coeff[f] * RAB[f].X * RAB[f].Y / rlen[f];
								FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+row_offset].Z -= force_wo_coeff[f] * RAB[f].X * RAB[f].Z / rlen[f];	
							
								FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].X -= force_wo_coeff[f] * RAB[f].Y * RAB[f].X / rlen[f];
							    	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].Y -= force_wo_coeff[f] * RAB[f].Y * RAB[f].Y / rlen[f];
								FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+row_offset].Z -= force_wo_coeff[f] * RAB[f].Y * RAB[f].Z / rlen[f];	
							
								FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].X -= force_wo_coeff[f] * RAB[f].Z * RAB[f].X / rlen[f];
							    	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].Y -= force_wo_coeff[f] * RAB[f].Z * RAB[f].Y / rlen[f];
								FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+row_offset].Z -= force_wo_coeff[f] * RAB[f].Z * RAB[f].Z / rlen[f];
							}	
						}
						
						if(CONTROLS.FIT_ENER) 
						{
							TMP_ENER  = fcut[0] 
							          * fcut[1] 
								  * fcut[2] 
								  * fcut[3] 
								  * fcut[4] 
								  * fcut[5];
								  
							TMP_ENER *=  Tn_ij[powers[0]] 
							           * Tn_ik[powers[1]] 
								   * Tn_il[powers[2]] 
								   * Tn_jk[powers[3]] 
								   * Tn_jl[powers[4]] 
								   * Tn_kl[powers[5]];
							
							FRAME_A_MATRIX[MATR_SIZE-1][vstart+row_offset].X += TMP_ENER;
							FRAME_A_MATRIX[MATR_SIZE-1][vstart+row_offset].Y += TMP_ENER;
							FRAME_A_MATRIX[MATR_SIZE-1][vstart+row_offset].Z += TMP_ENER;
						}
					}
					
				}	// End loop over 4th atom							
			}	// End loop over 3rd atom
		}	// End loop over 2nd atom
	}	// End loop over 1st atom

	
	if (CONTROLS.FIT_STRESS)
	{
		for(int i=0; i<CONTROLS.NUM_4B_CHEBY; i++) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS][n_2b_cheby_terms + n_3b_cheby_terms+i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS][n_2b_cheby_terms + n_3b_cheby_terms+i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS][n_2b_cheby_terms + n_3b_cheby_terms+i].Z *= inv_vol;	
		}
	}
	
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		for(int i=0; i<CONTROLS.NUM_4B_CHEBY; i++) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][n_2b_cheby_terms + n_3b_cheby_terms+i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][n_2b_cheby_terms + n_3b_cheby_terms+i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][n_2b_cheby_terms + n_3b_cheby_terms+i].Z *= inv_vol;	
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][n_2b_cheby_terms + n_3b_cheby_terms+i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][n_2b_cheby_terms + n_3b_cheby_terms+i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][n_2b_cheby_terms + n_3b_cheby_terms+i].Z *= inv_vol;
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][n_2b_cheby_terms + n_3b_cheby_terms+i].X *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][n_2b_cheby_terms + n_3b_cheby_terms+i].Y *= inv_vol;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][n_2b_cheby_terms + n_3b_cheby_terms+i].Z *= inv_vol;
		}
	}	
}

void Cheby::Force_all(CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS) 
// Calculate short-range forces using a Chebyshev polynomial expansion. Can use morse variables similar to the work of Bowman.
{
  // Variables exclusive to 2-body

  double xdiff_2b, xavg_2b;
  static double *Tn, *Tnd;
  static double fcut_2b, fcutderiv_2b, deriv;
  static double rpenalty, Vpenalty;
  static bool called_before = false ;

  // 3-body variables that may be shared with 2-body 
	 		
  XYZ RAB_IJ;
	
#if FORCECHECK == 1	
  static vector<XYZ> FORCE_3B;	// Equivalent of f3b 	
  static ofstream FILE_FORCE_3B;
#endif
	
  string TEMP_STR;
  int fidx_a2 ;
	
  ////////////////////////////////////////////////////////////////////////////////////////
	

  // A penalty function is added to the potential for r + penalty_dist < smin[ipair]
  // All pairs have the same penalty scale and distance
	 
  const double penalty_scale  = FF_2BODY[0].PENALTY_SCALE;	// 1.0e8;
  const double penalty_dist   = FF_2BODY[0].PENALTY_DIST;  	// 0.01;

  if ( ! called_before ) 
  {
	 called_before = true;
		
	 // Set up 2-body polynomials
		
	 int dim = 0;
		
	 for (int i=0; i<FF_2BODY.size(); i++)  
		if (FF_2BODY[i].SNUM > dim ) 
		  dim = FF_2BODY[i].SNUM;	 
		
	 dim++;
	 Tn   = new double [dim];
	 Tnd  = new double [dim];
		
#if FORCECHECK

	 FORCE_3B.resize(SYSTEM.ATOMS);

	 for( int i=0; i<SYSTEM.ATOMS; i++)
		FORCE_3B[i].X = FORCE_3B[i].Y = FORCE_3B[i].Z = 0;

	 FILE_FORCE_3B.open("3b_results.dat");

#endif 
	 // i.e the width of the default cheby range
  }

  // Main loop for Chebyshev terms:
	
  // Set up for MPI
	
  int a1start, a1end;	

#ifndef LINK_LAMMPS
  divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
#else
  a1start = SYSTEM.MY_ATOMS_START;
  a1end   = SYSTEM.MY_ATOMS_START+SYSTEM.MY_ATOMS-1;
#endif	

  // Set up for neighbor lists
	
  int a2start, a2end, a2;
	
  int BAD_CONFIG_1_FOUND = 0; // 0 == false, 1+ == true
  int BAD_CONFIG_2_FOUND = 0; // 0 == false, 1+ == true
	
  /////////////////////////////////////////////
  // EVALUATE THE 2-BODY INTERACTIONS
  /////////////////////////////////////////////
	
  if(FF_2BODY[0].SNUM>0)
  {
	 for(int a1=a1start; a1<=a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev -1)
	 {	
		a2start = 0;
		a2end   = NEIGHBOR_LIST.LIST[a1].size();
			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{
		  a2 = NEIGHBOR_LIST.LIST[a1][a2idx];			
			
		  int curr_pair_type_idx_ij =  get_pair_index(a1, a2, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP, SYSTEM.PARENT) ;
		
		  double rlen_ij = get_dist(SYSTEM, RAB_IJ, a1, a2);	// Updates RAB!
								
		  if(rlen_ij < FF_2BODY[curr_pair_type_idx_ij].S_MAXIM)	// We want to evaluate the penalty function when r < rmin (LEF) .. Assumes 3b inner cutoff is never shorter than 2b's
		  {	
		
			 /////////////////////////////////////////////
			 // EVALUATE THE 2-BODY INTERACTIONS
			 /////////////////////////////////////////////
			
			 // Make sure our newly transformed distance falls in defined range for Cheby polynomials and change the range, if the user requested
			 // Generate Chebyshev polynomials. 
							 
			 xdiff_2b = FF_2BODY[curr_pair_type_idx_ij].X_DIFF ;
			 xavg_2b  = FF_2BODY[curr_pair_type_idx_ij].X_AVG ;
			 set_polys(curr_pair_type_idx_ij, Tn, Tnd, rlen_ij, xdiff_2b, xavg_2b, FF_2BODY[curr_pair_type_idx_ij].SNUM);
				
			 FF_2BODY[curr_pair_type_idx_ij].FORCE_CUTOFF.get_fcut(fcut_2b, fcutderiv_2b, rlen_ij, FF_2BODY[curr_pair_type_idx_ij].S_MINIM,
																					 FF_2BODY[curr_pair_type_idx_ij].S_MAXIM);

			 fidx_a2 = SYSTEM.PARENT[a2];
				
			 for ( int i = 0; i < FF_2BODY[curr_pair_type_idx_ij].SNUM; i++ ) 
			 {
				double coeff                = FF_2BODY[curr_pair_type_idx_ij].PARAMS[i]; // This is the Cheby FF param for the given power
				SYSTEM.TOT_POT_ENER += coeff * fcut_2b * Tn[i+1];
				deriv                = (fcut_2b * Tnd[i+1] + fcutderiv_2b * Tn[i+1]);
				SYSTEM.PRESSURE_XYZ -= coeff * deriv * rlen_ij;				
					
				SYSTEM.PRESSURE_TENSORS_XYZ.X -= coeff * deriv * RAB_IJ.X * RAB_IJ.X / rlen_ij;
				SYSTEM.PRESSURE_TENSORS_XYZ.Y -= coeff * deriv * RAB_IJ.Y * RAB_IJ.Y / rlen_ij;
				SYSTEM.PRESSURE_TENSORS_XYZ.Z -= coeff * deriv * RAB_IJ.Z * RAB_IJ.Z / rlen_ij;
					
				SYSTEM.ACCEL[a1].X += coeff * deriv * RAB_IJ.X / rlen_ij;
				SYSTEM.ACCEL[a1].Y += coeff * deriv * RAB_IJ.Y / rlen_ij;
				SYSTEM.ACCEL[a1].Z += coeff * deriv * RAB_IJ.Z / rlen_ij;

				SYSTEM.ACCEL[fidx_a2].X -= coeff * deriv * RAB_IJ.X / rlen_ij;
				SYSTEM.ACCEL[fidx_a2].Y -= coeff * deriv * RAB_IJ.Y / rlen_ij;
				SYSTEM.ACCEL[fidx_a2].Z -= coeff * deriv * RAB_IJ.Z / rlen_ij;

			 }
								
			 // Add penalty for very short distances(less than smin + penalty_dist), where the fit FF may be unphysical (preserve conservation of E).

			 if ( rlen_ij - penalty_dist < FF_2BODY[curr_pair_type_idx_ij].S_MINIM ) 
				rpenalty = FF_2BODY[curr_pair_type_idx_ij].S_MINIM + penalty_dist - rlen_ij;
			 else 
				rpenalty = 0.0;		

			 if ( rpenalty > 0.0 ) 
			 {
			 	if(rlen_ij < (FF_2BODY[curr_pair_type_idx_ij].S_MINIM+penalty_dist)) // Then we've found a config that should be useful for self-consistent fitting
				{
					BAD_CONFIG_2_FOUND++;
					
					if(rlen_ij < FF_2BODY[curr_pair_type_idx_ij].S_MINIM) // Then we've found a config that should be useful for self-consistent fitting
				  		BAD_CONFIG_1_FOUND++;
				}
				
				Vpenalty = 0.0;
					
				if (isatty(fileno(stdout)))
				  cout << COUT_STYLE.BOLD << COUT_STYLE.MAGENTA << "Warning: (Step " << CONTROLS.STEP << ")Adding penalty in 2B Cheby calc, r < rmin+penalty_dist " << fixed << rlen_ij << " " << FF_2BODY[curr_pair_type_idx_ij].S_MINIM+penalty_dist << " " << TEMP_STR << " " << a1 << " " << a2 << COUT_STYLE.ENDSTYLE << endl;
				else
				  cout << "Warning: (Step " << CONTROLS.STEP << ") Adding penalty in 2B Cheby calc, r < rmin+penalty_dist " << fixed << rlen_ij << " " << FF_2BODY[curr_pair_type_idx_ij].S_MINIM+penalty_dist << " " << TEMP_STR << " " << a1 << " " << a2 << endl;
					
				SYSTEM.ACCEL[fidx_a2].X += 3.0 * rpenalty * rpenalty * penalty_scale * RAB_IJ.X / rlen_ij;
				SYSTEM.ACCEL[fidx_a2].Y += 3.0 * rpenalty * rpenalty * penalty_scale * RAB_IJ.Y / rlen_ij;
				SYSTEM.ACCEL[fidx_a2].Z += 3.0 * rpenalty * rpenalty * penalty_scale * RAB_IJ.Z / rlen_ij;
					
				SYSTEM.ACCEL[a1].X -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB_IJ.X / rlen_ij;
				SYSTEM.ACCEL[a1].Y -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB_IJ.Y / rlen_ij;
				SYSTEM.ACCEL[a1].Z -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB_IJ.Z / rlen_ij;							
					
				Vpenalty = rpenalty * rpenalty * rpenalty * penalty_scale;
				SYSTEM.TOT_POT_ENER += Vpenalty;
				cout << "	...Penalty potential = "<< Vpenalty << endl;
			 }					
		  }			
		}
		
	 }	// If 2-body interaction.
  }	// If 2-body interaction.

  /////////////////////////////////////////////
  // EVALUATE THE 3-BODY INTERACTIONS
  /////////////////////////////////////////////

  if(FF_2BODY[0].SNUM_3B_CHEBY>0)
  { 
	 Force_3B(TRIPS) ;
  }  // If 3-body interaction.
	
	
  /////////////////////////////////////////////
  // EVALUATE THE 4-BODY INTERACTIONS
  /////////////////////////////////////////////

  if(FF_2BODY[0].SNUM_4B_CHEBY>0)
  {
	 Force_4B(QUADS) ;
  }
	
  // Check if a truly bad configuration was found, and if so, print it out
	
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &BAD_CONFIG_1_FOUND,1,MPI_INT, MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &BAD_CONFIG_2_FOUND,1,MPI_INT, MPI_SUM,MPI_COMM_WORLD);
#endif

  // traj file for BAD_CONFIG_1 should only contain configs where rij<rcutin
  // and BAD_CONFIG_2 should only contain configs where rcutin < rij < rcutin+dp.
  
  if  ((BAD_CONFIG_1_FOUND>0) && CONTROLS.PRINT_BAD_CFGS)
  	PRINT_CONFIG(SYSTEM, CONTROLS,1);
   else if ((BAD_CONFIG_2_FOUND>0) && CONTROLS.PRINT_BAD_CFGS)
	PRINT_CONFIG(SYSTEM, CONTROLS,2);

//cout << "COUNTED INTERACTIONS: " << COUNTED_INTERACTIONS << endl;	
} 


void Cheby::Force_3B(CLUSTER_LIST &TRIPS)
// Evaluate the 3-Body Force.
{

  int i_start, i_end;
  vector<TRIP_FF> & FF_3BODY = TRIPS.VEC ;

		
  double rlen_ij,  rlen_ik,  rlen_jk;
  double rlen_ij_dummy, rlen_ik_dummy, rlen_jk_dummy;
	
  static double *Tn_ij,  *Tn_ik,  *Tn_jk;
  static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
  static bool    called_before = false;
	
  int pow_ij, pow_ik, pow_jk;
			  
  double fcut_ij,  fcut_ik,  fcut_jk; 			
  double deriv_ij, deriv_ik, deriv_jk;
  double force_ij, force_ik, force_jk;
  double fcutderiv_ij, fcutderiv_ik, fcutderiv_jk; 		
  int curr_triple_type_index;
  int curr_pair_type_idx_ij;
  int curr_pair_type_idx_ik;
  int curr_pair_type_idx_jk;
  double coeff;
  vector<int> atom_type_idx(3) ;	
  vector<double> s_maxim(3), s_minim(3), x_avg(3), x_diff(3) ;
  vector<int> pair_index(3) ;

  int fidx_a2, fidx_a3 ;

  if ( ! called_before ) 
  {
		
	 // Set up 3-body polynomials
		 
	 int dim = 0;

	 for (int i=0; i<FF_2BODY.size(); i++) 
		if (FF_2BODY[i].SNUM_3B_CHEBY > dim ) 
		  dim = FF_2BODY[i].SNUM_3B_CHEBY;	

	 dim++;
	 Tn_ij   = new double [dim];
	 Tn_ik   = new double [dim];
	 Tn_jk   = new double [dim];

	 Tnd_ij  = new double [dim];
	 Tnd_ik  = new double [dim];
	 Tnd_jk  = new double [dim]; 

	 called_before = true ;
  }

#ifndef LINK_LAMMPS
  divide_atoms(i_start, i_end, NEIGHBOR_LIST.LIST_3B_INT.size());	
#else
  // Not sure if this is correct for LAMMPS.  Please check ! (Larry) -- it looks correct (RKL)
  i_start = 0;
  i_end = NEIGHBOR_LIST.LIST_3B_INT.size() - 1;
  a1start = SYSTEM.MY_ATOMS_START;
  a1end   = SYSTEM.MY_ATOMS_START+SYSTEM.MY_ATOMS-1;
#endif	
		
  int a1;
  int INTERACTIONS = 0;

  // Loop over a1, a2, a3 interaction triples, not atoms
  for ( int ii = i_start; ii <= i_end; ii++ ) 
  {

	 a1 = NEIGHBOR_LIST.LIST_3B_INT[ii].a1;

#ifdef LINK_LAMMPS
	 if ( a1 < a1start || a1 > a1end ) continue;
#endif
			
	 int a2 = NEIGHBOR_LIST.LIST_3B_INT[ii].a2;
	 int a3 = NEIGHBOR_LIST.LIST_3B_INT[ii].a3;

	 if ( a3 == a2 || SYSTEM.PARENT[a2] > SYSTEM.PARENT[a3] ) {
		cout << "Bad pair found " << a3 << a2 << endl;
		cout << "Parents " << SYSTEM.PARENT[a2] << SYSTEM.PARENT[a3] << endl;
	 }
			
	 vector<int> atom_type_index(3) ;

	 atom_type_index[0] = SYSTEM.get_atomtype_idx(a1) ;
	 atom_type_index[1] = SYSTEM.get_atomtype_idx(a2) ;
	 atom_type_index[2] = SYSTEM.get_atomtype_idx(a3) ;

	 curr_pair_type_idx_ij =  get_pair_index(a1, a2, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
														  SYSTEM.PARENT) ;
	 curr_pair_type_idx_ik =  get_pair_index(a1, a3, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
														  SYSTEM.PARENT) ;
	 curr_pair_type_idx_jk =  get_pair_index(a2, a3, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
														  SYSTEM.PARENT) ;

	 int tidx = TRIPS.make_id_int(atom_type_index) ;
	 curr_triple_type_index = TRIPS.INT_MAP[tidx];
					
	 if(curr_triple_type_index<0) 
	 {
		//cout << "Excluded interaction: " << tidx << endl ;
		continue;
	 }

	 XYZ RAB_IJ, RAB_IK, RAB_JK ;
	 rlen_ij = get_dist(SYSTEM, RAB_IJ, a1, a2);	// Updates RAB!
	 rlen_ik = get_dist(SYSTEM, RAB_IK, a1, a3);	// Updates RAB!
	 rlen_jk = get_dist(SYSTEM, RAB_JK, a2, a3);	// Updates RAB!

	 for ( int j = 0 ; j < 3 ; j++ ) 
	 {
		pair_index[j] = TRIPS.PAIR_INDICES[tidx][j] ;
		s_maxim[j] = FF_3BODY[curr_triple_type_index].S_MAXIM[pair_index[j]] ;
		s_minim[j] = FF_3BODY[curr_triple_type_index].S_MINIM[pair_index[j]] ;
		x_diff[j]  = FF_3BODY[curr_triple_type_index].X_DIFF[pair_index[j]] ;
		x_avg[j]	  = FF_3BODY[curr_triple_type_index].X_AVG[pair_index[j]] ;
	 }
					
	 // Before doing any polynomial/coeff set up, make sure that all ij, ik, and jk distances are 
	 // within the allowed range.

	 if ( FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_ij, s_minim[0], s_maxim[0]) )
	 {
		if ( FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_ik, s_minim[1], s_maxim[1]))
		{
		  if ( FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.PROCEED(rlen_jk, s_minim[2], s_maxim[2]) )
		  {											
			 // Everything is within allowed ranges. Begin setting up the force calculation
			  
			 INTERACTIONS++;
			 //cout << a1 << " " << a2 << " " << a3 << " " << endl;

			 rlen_ij_dummy = rlen_ij;
			 rlen_ik_dummy = rlen_ik;
			 rlen_jk_dummy = rlen_jk;
			  
			 if(rlen_ij < s_minim[0])
				rlen_ij_dummy = s_minim[0] ;
			  
			 if(rlen_ik < s_minim[1])
				rlen_ik_dummy = s_minim[1] ;
							
			 if(rlen_jk < s_minim[2])
				rlen_jk_dummy = s_minim[2] ;

			 // Set up the polynomials
			  
			 set_polys(curr_pair_type_idx_ij, Tn_ij, Tnd_ij, rlen_ij_dummy, x_diff[0], x_avg[0], FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY);
			 set_polys(curr_pair_type_idx_ik, Tn_ik, Tnd_ik, rlen_ik_dummy, x_diff[1], x_avg[1], FF_2BODY[curr_pair_type_idx_ik].SNUM_3B_CHEBY);
			 set_polys(curr_pair_type_idx_jk, Tn_jk, Tnd_jk, rlen_jk_dummy, x_diff[2], x_avg[2], FF_2BODY[curr_pair_type_idx_jk].SNUM_3B_CHEBY);
												
			 // Apply the FF

			 // Set up the smoothing functions

			 FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_ij, fcutderiv_ij, rlen_ij, s_minim[0], s_maxim[0]);
			 FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_ik, fcutderiv_ik, rlen_ik, s_minim[1], s_maxim[1]);
			 FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_jk, fcutderiv_jk, rlen_jk, s_minim[2], s_maxim[2]);

			 // Now compute the forces for each set of allowed powers for pairs ij, ik, and jk		
			 // Keep in mind that the order in which allowed powers are stored may not match the
			 // ordering of pairs resulting from the present atom triplet. Thus, we need to order
			 // the stored powers properly before applying the FF.
			  
			 // When doing a compare force calculation, make sure that forces on 
			 // replicate atoms are attributed to the parent atoms
	
			 //fidx_a1 = SYSTEM.PARENT[a1];
			 fidx_a2 = SYSTEM.PARENT[a2];
			 fidx_a3 = SYSTEM.PARENT[a3];
							
			 for(int i=0; i<FF_3BODY[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
			 {
				set_3b_powers(FF_3BODY[curr_triple_type_index], pair_index, i,
								  pow_ij, pow_ik, pow_jk) ;
			      
				coeff = FF_3BODY[curr_triple_type_index].PARAMS[i];
			      
				SYSTEM.TOT_POT_ENER += coeff * fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk]; 

				deriv_ij  = fcut_ij * Tnd_ij[pow_ij] + fcutderiv_ij * Tn_ij [pow_ij];
				deriv_ik  = fcut_ik * Tnd_ik[pow_ik] + fcutderiv_ik * Tn_ik [pow_ik];
				deriv_jk  = fcut_jk * Tnd_jk[pow_jk] + fcutderiv_jk * Tn_jk [pow_jk];
									
				force_ij  = coeff * deriv_ij * fcut_ik * fcut_jk * Tn_ik [pow_ik] * Tn_jk [pow_jk];
				force_ik  = coeff * deriv_ik * fcut_ij * fcut_jk * Tn_ij [pow_ij] * Tn_jk [pow_jk];
				force_jk  = coeff * deriv_jk * fcut_ij * fcut_ik * Tn_ij [pow_ij] * Tn_ik [pow_ik];
							
				SYSTEM.PRESSURE_XYZ    -= force_ij * rlen_ij;
				SYSTEM.PRESSURE_XYZ    -= force_ik * rlen_ik;
				SYSTEM.PRESSURE_XYZ    -= force_jk * rlen_jk;
							
				force_ij /= rlen_ij;
				force_ik /= rlen_ik;
				force_jk /= rlen_jk;

				SYSTEM.PRESSURE_TENSORS_XYZ.X -= force_ij * RAB_IJ.X * RAB_IJ.X ;
				SYSTEM.PRESSURE_TENSORS_XYZ.Y -= force_ij * RAB_IJ.Y * RAB_IJ.Y ;
				SYSTEM.PRESSURE_TENSORS_XYZ.Z -= force_ij * RAB_IJ.Z * RAB_IJ.Z ;
							
				SYSTEM.PRESSURE_TENSORS_XYZ.X -= force_ik * RAB_IK.X * RAB_IK.X ;
				SYSTEM.PRESSURE_TENSORS_XYZ.Y -= force_ik * RAB_IK.Y * RAB_IK.Y ;
				SYSTEM.PRESSURE_TENSORS_XYZ.Z -= force_ik * RAB_IK.Z * RAB_IK.Z ;
							
				SYSTEM.PRESSURE_TENSORS_XYZ.X -= force_jk * RAB_JK.X * RAB_JK.X ;
				SYSTEM.PRESSURE_TENSORS_XYZ.Y -= force_jk * RAB_JK.Y * RAB_JK.Y ;
				SYSTEM.PRESSURE_TENSORS_XYZ.Z -= force_jk * RAB_JK.Z * RAB_JK.Z ;
							

				// Apply forces to ij pair
			      
				SYSTEM.ACCEL[a1]     .X += force_ij * RAB_IJ.X;
				SYSTEM.ACCEL[a1]     .Y += force_ij * RAB_IJ.Y;
				SYSTEM.ACCEL[a1]     .Z += force_ij * RAB_IJ.Z;
			
				SYSTEM.ACCEL[fidx_a2].X -= force_ij * RAB_IJ.X;
				SYSTEM.ACCEL[fidx_a2].Y -= force_ij * RAB_IJ.Y;
				SYSTEM.ACCEL[fidx_a2].Z -= force_ij * RAB_IJ.Z;
			
				// Apply forces to ik pair
			      
				SYSTEM.ACCEL[a1]     .X += force_ik * RAB_IK.X;
				SYSTEM.ACCEL[a1]     .Y += force_ik * RAB_IK.Y;
				SYSTEM.ACCEL[a1]     .Z += force_ik * RAB_IK.Z;	
			
				SYSTEM.ACCEL[fidx_a3].X -= force_ik * RAB_IK.X;
				SYSTEM.ACCEL[fidx_a3].Y -= force_ik * RAB_IK.Y;
				SYSTEM.ACCEL[fidx_a3].Z -= force_ik * RAB_IK.Z;	
			
				// Apply forces to jk pair
			      
				SYSTEM.ACCEL[fidx_a2].X += force_jk * RAB_JK.X;
				SYSTEM.ACCEL[fidx_a2].Y += force_jk * RAB_JK.Y;
				SYSTEM.ACCEL[fidx_a2].Z += force_jk * RAB_JK.Z;
			
				SYSTEM.ACCEL[fidx_a3].X -= force_jk * RAB_JK.X;
				SYSTEM.ACCEL[fidx_a3].Y -= force_jk * RAB_JK.Y;
				SYSTEM.ACCEL[fidx_a3].Z -= force_jk * RAB_JK.Z;	

			      
#if FORCECHECK
			      
				// Apply forces to ij pair
			      
				FORCE_3B[a1]     .X += force_ij * RAB_IJ.X;
				FORCE_3B[a1]     .Y += force_ij * RAB_IJ.Y;
				FORCE_3B[a1]     .Z += force_ij * RAB_IJ.Z;
			      
				FORCE_3B[fidx_a2].X -= force_ij * RAB_IJ.X;
				FORCE_3B[fidx_a2].Y -= force_ij * RAB_IJ.Y;
				FORCE_3B[fidx_a2].Z -= force_ij * RAB_IJ.Z;
			
				// Apply forces to ik pair
			      
				FORCE_3B[a1]     .X += force_ik * RAB_IK.X;
				FORCE_3B[a1]     .Y += force_ik * RAB_IK.Y;
				FORCE_3B[a1]     .Z += force_ik * RAB_IK.Z;	
			
				FORCE_3B[fidx_a3].X -= force_ik * RAB_IK.X;
				FORCE_3B[fidx_a3].Y -= force_ik * RAB_IK.Y;
				FORCE_3B[fidx_a3].Z -= force_ik * RAB_IK.Z;	
			
				// Apply forces to jk pair
			      
				FORCE_3B[fidx_a2].X += force_jk * RAB_JK.X;
				FORCE_3B[fidx_a2].Y += force_jk * RAB_JK.Y;
				FORCE_3B[fidx_a2].Z += force_jk * RAB_JK.Z;
			      
				FORCE_3B[fidx_a3].X -= force_jk * RAB_JK.X;
				FORCE_3B[fidx_a3].Y -= force_jk * RAB_JK.Y;
				FORCE_3B[fidx_a3].Z -= force_jk * RAB_JK.Z;										
#endif								
			 }	
		  }	
		}				
	 }		
  } 	// Loop over interactions.
}


void Cheby::Force_4B(CLUSTER_LIST &QUADS)
{
  // Prepare iterators for outermost loop
		
  int i_start, i_end;
		
  static bool called_before = false ;

  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
	
  // 4-BODY VARIABLES
	
  ////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////
	
  vector<XYZ> RVEC(6);	// Replaces RVEC_IJ, RVEC_IK...
  vector<XYZ> RAB (6);	// Replaces RAB_IJ, RAB_IK...
	
  vector<double> rlen(6);			// Replaces rlen_ij, rlen_ik...
  vector<double> rlen_dummy(6);	// Replaces rlen_ij_dummy, rlen_ik_dummy

  static double *Tn_4b_ij,  *Tn_4b_ik,  *Tn_4b_il,  *Tn_4b_jk,  *Tn_4b_jl,  *Tn_4b_kl;
  static double *Tnd_4b_ij, *Tnd_4b_ik, *Tnd_4b_il, *Tnd_4b_jk, *Tnd_4b_jl, *Tnd_4b_kl;
	
  vector<int> powers(6);	// replaces pow_ij, pow_ik, pow_jk;
  vector<double> x_diff(6), x_avg(6);		// replaces xdiff_ij, xdiff_ik, xdiff_jk; 
  vector<double> fcut_4b(6);		// replaces cut_ij,  fcut_ik,  fcut_jk;
  vector<double> fcut_deriv_4b(6);// replaces fcutderiv_ij, fcutderiv_ik, fcutderiv_jk; 
  vector<double> deriv_4b(6);		// replaces deriv_ij, deriv_ik, deriv_jk;
  vector<double> force_4b(6);		// replaces force

	
//	static string TEMP_STR;
  int curr_quad_type_index;
  vector<int> curr_pair_type_idx(6);	// replaces curr_pair_type_idx_ij, etc
  vector<double> S_MAXIM(6);	// replaces S_MAXIM_IJ, S_MAXIM_IK, S_MAXIM_JK;
  vector<double> S_MINIM(6);	// replaces S_MINIM_IJ, S_MINIM_IK, S_MINIM_JK;
	 
  int quad_id_int;
  vector<int> pow_map(6);
  vector<int> atom_type_index(4);		// Index of type of atoms in the quad cluster.

  vector<CLUSTER>& FF_4BODY = QUADS.VEC ;

#ifndef LINK_LAMMPS
  divide_atoms(i_start, i_end, NEIGHBOR_LIST.LIST_4B_INT.size());	
#else
  i_start = 0;
  i_end = NEIGHBOR_LIST.LIST_4B_INT.size() - 1;
  a1start = SYSTEM.MY_ATOMS_START;
  a1end   = SYSTEM.MY_ATOMS_START+SYSTEM.MY_ATOMS-1;
#endif	

  if ( ! called_before ) 
  {


	 // Set up 4-body polynomials
		
	 int dim = 0 ;

	 for (int i=0; i<FF_2BODY.size(); i++) 
		if (FF_2BODY[i].SNUM_4B_CHEBY > dim ) 
		  dim = FF_2BODY[i].SNUM_4B_CHEBY;	

	 dim++;		
	 Tn_4b_ij   = new double [dim];
	 Tn_4b_ik   = new double [dim];
	 Tn_4b_il   = new double [dim];
	 Tn_4b_jk   = new double [dim];
	 Tn_4b_jl   = new double [dim];
	 Tn_4b_kl   = new double [dim];

	 Tnd_4b_ij  = new double [dim];
	 Tnd_4b_ik  = new double [dim];
	 Tnd_4b_il  = new double [dim];
	 Tnd_4b_jk  = new double [dim];
	 Tnd_4b_jl  = new double [dim];
	 Tnd_4b_kl  = new double [dim];
	 
	 called_before = true ;
  }
		
  ////////////////////////////////////////////////////////////////////////////////////////

  int a1;
		
  // Loop over a1, a2, a3, a4 interaction quadruplets, not atoms

  for ( int ii = i_start; ii <= i_end; ii++ ) 
  {
  
	 a1 = NEIGHBOR_LIST.LIST_4B_INT[ii].a1;

#ifdef LINK_LAMMPS
	 if ( a1 < a1start || a1 > a1end ) continue;
#endif
			
	 int a2 = NEIGHBOR_LIST.LIST_4B_INT[ii].a2;
	 int a3 = NEIGHBOR_LIST.LIST_4B_INT[ii].a3;
	 int a4 = NEIGHBOR_LIST.LIST_4B_INT[ii].a4;

	 if ( a3 == a2 || a4 == a2 || a3 == a4 ) 
	 {
		cout << "Bad pair found " << a1 << ", " << a2 << ", " << a3 << ", " << a4 << endl;
		cout << "Parents " << SYSTEM.PARENT[a1] << ", " << SYSTEM.PARENT[a2] << ", " <<  SYSTEM.PARENT[a3] << ", " <<  SYSTEM.PARENT[a4] << endl;
	 }
	 if (SYSTEM.PARENT[a2] > SYSTEM.PARENT[a3] || SYSTEM.PARENT[a2] > SYSTEM.PARENT[a4] || SYSTEM.PARENT[a3] > SYSTEM.PARENT[a4]) 
	 {
		cout << "Bad pair found " << a1 << ", " << a2 << ", " << a3 << ", " << a4 << endl;
		cout << "Parents " << SYSTEM.PARENT[a1] << ", " << SYSTEM.PARENT[a2] << ", " <<  SYSTEM.PARENT[a3] << ", " <<  SYSTEM.PARENT[a4] << endl;
	 }			
			
	 int fidx_a2 = SYSTEM.PARENT[a2];
	 int fidx_a3 = SYSTEM.PARENT[a3];
	 int fidx_a4 = SYSTEM.PARENT[a4];
			
	 curr_pair_type_idx[0] =  get_pair_index(a1, a2, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
																		SYSTEM.PARENT) ;
	 curr_pair_type_idx[1] =  get_pair_index(a1, a3, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
																		SYSTEM.PARENT) ;
	 curr_pair_type_idx[2] =  get_pair_index(a1, a4, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
																		SYSTEM.PARENT) ;
	 curr_pair_type_idx[3] =  get_pair_index(a2, a3, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
																		SYSTEM.PARENT) ;
	 curr_pair_type_idx[4] =  get_pair_index(a2, a4, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
																		SYSTEM.PARENT) ;
	 curr_pair_type_idx[5] =  get_pair_index(a3, a4, SYSTEM.ATOMTYPE_IDX, CONTROLS.NATMTYP,
																		SYSTEM.PARENT) ;
	 atom_type_index[0] = SYSTEM.get_atomtype_idx(a1) ;
	 atom_type_index[1] = SYSTEM.get_atomtype_idx(fidx_a2) ;
	 atom_type_index[2] = SYSTEM.get_atomtype_idx(fidx_a3) ;
	 atom_type_index[3] = SYSTEM.get_atomtype_idx(fidx_a4) ;
			
	 // Always construct quad_id_int lookup key based on atom types in decending order

	 //sort   (atom_type_index.begin(), atom_type_index.end());
	 //reverse(atom_type_index.begin(), atom_type_index.end());
			
	 quad_id_int         = QUADS.make_id_int(atom_type_index) ;
	 curr_quad_type_index = QUADS.INT_MAP[quad_id_int];

	 if(curr_quad_type_index<0)
	 {
		// cout << "Skipping QUAD " << quad_id_int << endl ;
		continue;
	 }

	 // Get the atom distances

	 rlen[0] = get_dist(SYSTEM, RAB[0], a1, a2);	// Updates RAB!
	 rlen[1] = get_dist(SYSTEM, RAB[1], a1, a3);	// Updates RAB!
	 rlen[2] = get_dist(SYSTEM, RAB[2], a1, a4);	// Updates RAB!
	 rlen[3] = get_dist(SYSTEM, RAB[3], a2, a3);	// Updates RAB!
	 rlen[4] = get_dist(SYSTEM, RAB[4], a2, a4);	// Updates RAB!
	 rlen[5] = get_dist(SYSTEM, RAB[5], a3, a4);	// Updates RAB!
			
	 // Determine the inner and outer cutoffs for each pair type in the quadruplet
		
	 for (int f=0; f<6; f++)
	 {
		int j ;
		pow_map[f] = QUADS.PAIR_INDICES[quad_id_int][f] ;
		j = pow_map[f] ;
		S_MAXIM[f] = FF_4BODY[curr_quad_type_index].S_MAXIM[j] ;
		S_MINIM[f] = FF_4BODY[curr_quad_type_index].S_MINIM[j] ;
		x_diff[f] = FF_4BODY[curr_quad_type_index].X_DIFF[j] ;
		x_avg[f]  = FF_4BODY[curr_quad_type_index].X_AVG[j] ;
	 }

	 // Before doing any polynomial/coeff set up, make sure that all ij, ik, and jk distances are within the allowed range.
	 // Unlike the 2-body Cheby, extrapolation/refitting to handle behavior outside of fitting regime is not straightforward.
			
	 if( !FF_4BODY[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[0], S_MINIM[0], S_MAXIM[0]))
		continue;
	 if( !FF_4BODY[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[1], S_MINIM[1], S_MAXIM[1]))
		continue;
	 if( !FF_4BODY[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[2], S_MINIM[2], S_MAXIM[2]))
		continue;
	 if( !FF_4BODY[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[3], S_MINIM[3], S_MAXIM[3]))
		continue;
	 if( !FF_4BODY[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[4], S_MINIM[4], S_MAXIM[4]))
		continue;
	 if( !FF_4BODY[curr_quad_type_index].FORCE_CUTOFF.PROCEED(rlen[5], S_MINIM[5], S_MAXIM[5]))
		continue;	

// At this point, all distances are within allowed ranges. We can now proceed to the force derivative calculation
			
	 // For all types, if r < rcut, then the potential is constant, thus the force  must be zero.
	 // Additionally, the potential is then taken to be the potential at r_cut.

	 for (int f=0; f<6; f++)
	 {
		rlen_dummy[f] = rlen[f];
				
		if(rlen[f] < S_MINIM[f])
		  rlen_dummy[f] = S_MINIM[f];
	 }

	 
	 // Set up the polynomials

	 set_polys(curr_pair_type_idx[0], Tn_4b_ij, Tnd_4b_ij, rlen_dummy[0], x_diff[0], x_avg[0], FF_2BODY[curr_pair_type_idx[0]].SNUM_4B_CHEBY);
	 set_polys(curr_pair_type_idx[1], Tn_4b_ik, Tnd_4b_ik, rlen_dummy[1], x_diff[1], x_avg[1], FF_2BODY[curr_pair_type_idx[1]].SNUM_4B_CHEBY);
	 set_polys(curr_pair_type_idx[2], Tn_4b_il, Tnd_4b_il, rlen_dummy[2], x_diff[2], x_avg[2], FF_2BODY[curr_pair_type_idx[2]].SNUM_4B_CHEBY);
	 set_polys(curr_pair_type_idx[3], Tn_4b_jk, Tnd_4b_jk, rlen_dummy[3], x_diff[3], x_avg[3], FF_2BODY[curr_pair_type_idx[3]].SNUM_4B_CHEBY);
	 set_polys(curr_pair_type_idx[4], Tn_4b_jl, Tnd_4b_jl, rlen_dummy[4], x_diff[4], x_avg[4], FF_2BODY[curr_pair_type_idx[4]].SNUM_4B_CHEBY);
	 set_polys(curr_pair_type_idx[5], Tn_4b_kl, Tnd_4b_kl, rlen_dummy[5], x_diff[5], x_avg[5], FF_2BODY[curr_pair_type_idx[5]].SNUM_4B_CHEBY);

	 // Set up the smoothing functions
			
	 for (int f=0; f<6; f++)
		FF_4BODY[curr_quad_type_index].FORCE_CUTOFF.get_fcut(fcut_4b[f], fcut_deriv_4b[f], rlen[f], S_MINIM[f], S_MAXIM[f]);
			
	 // Set up terms for derivatives
			
	 for(int i=0; i<FF_4BODY[curr_quad_type_index].N_ALLOWED_POWERS; i++) 
	 {
		double coeff = FF_4BODY[curr_quad_type_index].PARAMS[i];

		for ( int f = 0 ; f < 6 ; f++ ) 
		  powers[f] = FF_4BODY[curr_quad_type_index].ALLOWED_POWERS[i][pow_map[f]] ;


		SYSTEM.TOT_POT_ENER += coeff 
		                     * fcut_4b[0] 
				     * fcut_4b[1] 
				     * fcut_4b[2] 
				     * fcut_4b[3] 
				     * fcut_4b[4] 
				     * fcut_4b[5] 
		  * Tn_4b_ij[powers[0]] 
		  * Tn_4b_ik[powers[1]] 
		  * Tn_4b_il[powers[2]] 
		  * Tn_4b_jk[powers[3]] 
		  * Tn_4b_jl[powers[4]] 
		  * Tn_4b_kl[powers[5]]; 			

		deriv_4b[0] = fcut_4b[0] * Tnd_4b_ij[powers[0]] + fcut_deriv_4b[0] * Tn_4b_ij[powers[0]];
		deriv_4b[1] = fcut_4b[1] * Tnd_4b_ik[powers[1]] + fcut_deriv_4b[1] * Tn_4b_ik[powers[1]];
		deriv_4b[2] = fcut_4b[2] * Tnd_4b_il[powers[2]] + fcut_deriv_4b[2] * Tn_4b_il[powers[2]];
		deriv_4b[3] = fcut_4b[3] * Tnd_4b_jk[powers[3]] + fcut_deriv_4b[3] * Tn_4b_jk[powers[3]];
		deriv_4b[4] = fcut_4b[4] * Tnd_4b_jl[powers[4]] + fcut_deriv_4b[4] * Tn_4b_jl[powers[4]];
		deriv_4b[5] = fcut_4b[5] * Tnd_4b_kl[powers[5]] + fcut_deriv_4b[5] * Tn_4b_kl[powers[5]];
				
		force_4b[0]  = coeff * deriv_4b[0] * fcut_4b[1] * fcut_4b[2] * fcut_4b[3] * fcut_4b[4] * fcut_4b[5] * Tn_4b_ik[powers[1]]  * Tn_4b_il[powers[2]]  * Tn_4b_jk[powers[3]]  * Tn_4b_jl[powers[4]]  * Tn_4b_kl[powers[5]];
		force_4b[1]  = coeff * deriv_4b[1] * fcut_4b[0] * fcut_4b[2] * fcut_4b[3] * fcut_4b[4] * fcut_4b[5] * Tn_4b_ij[powers[0]]  * Tn_4b_il[powers[2]]  * Tn_4b_jk[powers[3]]  * Tn_4b_jl[powers[4]]  * Tn_4b_kl[powers[5]];
		force_4b[2]  = coeff * deriv_4b[2] * fcut_4b[0] * fcut_4b[1] * fcut_4b[3] * fcut_4b[4] * fcut_4b[5] * Tn_4b_ij[powers[0]]  * Tn_4b_ik[powers[1]]  * Tn_4b_jk[powers[3]]  * Tn_4b_jl[powers[4]]  * Tn_4b_kl[powers[5]];
		force_4b[3]  = coeff * deriv_4b[3] * fcut_4b[0] * fcut_4b[1] * fcut_4b[2] * fcut_4b[4] * fcut_4b[5] * Tn_4b_ij[powers[0]]  * Tn_4b_ik[powers[1]]  * Tn_4b_il[powers[2]]  * Tn_4b_jl[powers[4]]  * Tn_4b_kl[powers[5]];
		force_4b[4]  = coeff * deriv_4b[4] * fcut_4b[0] * fcut_4b[1] * fcut_4b[2] * fcut_4b[3] * fcut_4b[5] * Tn_4b_ij[powers[0]]  * Tn_4b_ik[powers[1]]  * Tn_4b_il[powers[2]]  * Tn_4b_jk[powers[3]]  * Tn_4b_kl[powers[5]];
		force_4b[5]  = coeff * deriv_4b[5] * fcut_4b[0] * fcut_4b[1] * fcut_4b[2] * fcut_4b[3] * fcut_4b[4] * Tn_4b_ij[powers[0]]  * Tn_4b_ik[powers[1]]  * Tn_4b_il[powers[2]]  * Tn_4b_jk[powers[3]]  * Tn_4b_jl[powers[4]];

#if(0)
		if ( a1 == 1 || fidx_a2 == 1 || fidx_a3 == 1 || fidx_a4 == 1 )
		{
		  cout << "Target atom found\n" ;

		  cout << std::scientific ;
		  cout.precision(10) ;
		  cout << "A1 = " << a1 << " A2 = " << fidx_a2 << " A3 = " << fidx_a3 << " A4 = " << fidx_a4 << endl ;
		  cout << "rlen_dummy: " ;
		  for ( int ifl = 0 ; ifl < 6 ; ifl++ ) 
		  {
			 cout << " " << rlen_dummy[ifl] ;
		  }
		  cout << endl ;
		  cout << "Force_wo_coeff: " ;
		  for ( int ifl = 0 ; ifl < 6 ; ifl++ ) 
		  {
			 if ( fabs(coeff) > 0.0 ) 
			 {
				cout << " " << force_4b[ifl] / coeff ;
			 } 
			 else
			 {
				cout << " " << 0.0  ;
			 }
		  }
		  cout << endl ;		  

		  cout << "Deriv_4b: " ;
		  for ( int ifl = 0 ; ifl < 6 ; ifl++ ) 
		  {
			 cout << " " << deriv_4b[ifl] ;
		  }
		  cout << endl ;		  

		  cout << "Tn_il " << Tn_4b_il[1] << endl ;
		  cout << "Tnd_il " << Tnd_4b_il[1] << endl ;

		  for ( int ifl = 0 ; ifl < 6 ; ifl++ ) 
		  {
			 cout << "Index = " << ifl << endl ;
			 cout << "Fcut " << fcut_4b[ifl] << endl ;
			 cout << "Fcut_deriv " << fcut_deriv_4b[ifl] << endl ;
			 cout << "Powers " << powers[ifl]] << endl ;
			 cout << "Pow_map " << pow_map[ifl] << endl ;
		  }
		  cout << "Allowed_powers" ;
		  for ( int ifl = 0 ; ifl < 6 ; ifl++ ) 
			 cout << FF_4BODY[curr_quad_type_index].ALLOWED_POWERS[i][ifl] << " " ;
		  cout << endl ;
		  cout.unsetf(ios_base::scientific) ;
		}
#endif

		for(int j=0; j<6; j++)
		{
		  SYSTEM.PRESSURE_XYZ -= force_4b[j] * rlen[j];

		  force_4b[j] /= rlen[j];
					
		  SYSTEM.PRESSURE_TENSORS_XYZ.X -= force_4b[j] * RAB[j].X * RAB[j].X ;
		  SYSTEM.PRESSURE_TENSORS_XYZ.Y -= force_4b[j] * RAB[j].Y * RAB[j].Y ;
		  SYSTEM.PRESSURE_TENSORS_XYZ.Z -= force_4b[j] * RAB[j].Z * RAB[j].Z ;
		}

		// Apply forces to ij pair
				
		SYSTEM.ACCEL[a1]     .X += force_4b[0] * RAB[0].X;
		SYSTEM.ACCEL[a1]     .Y += force_4b[0] * RAB[0].Y;
		SYSTEM.ACCEL[a1]     .Z += force_4b[0] * RAB[0].Z;

		SYSTEM.ACCEL[fidx_a2].X -= force_4b[0] * RAB[0].X;
		SYSTEM.ACCEL[fidx_a2].Y -= force_4b[0] * RAB[0].Y;
		SYSTEM.ACCEL[fidx_a2].Z -= force_4b[0] * RAB[0].Z;				
				
		// Apply forces to ik pair
				
		SYSTEM.ACCEL[a1]     .X += force_4b[1] * RAB[1].X;
		SYSTEM.ACCEL[a1]     .Y += force_4b[1] * RAB[1].Y;
		SYSTEM.ACCEL[a1]     .Z += force_4b[1] * RAB[1].Z;

		SYSTEM.ACCEL[fidx_a3].X -= force_4b[1] * RAB[1].X;
		SYSTEM.ACCEL[fidx_a3].Y -= force_4b[1] * RAB[1].Y;
		SYSTEM.ACCEL[fidx_a3].Z -= force_4b[1] * RAB[1].Z;						
				
		// Apply forces to il pair
				
		SYSTEM.ACCEL[a1]     .X += force_4b[2] * RAB[2].X;
		SYSTEM.ACCEL[a1]     .Y += force_4b[2] * RAB[2].Y;
		SYSTEM.ACCEL[a1]     .Z += force_4b[2] * RAB[2].Z;

		SYSTEM.ACCEL[fidx_a4].X -= force_4b[2] * RAB[2].X;
		SYSTEM.ACCEL[fidx_a4].Y -= force_4b[2] * RAB[2].Y;
		SYSTEM.ACCEL[fidx_a4].Z -= force_4b[2] * RAB[2].Z;					
				
		// Apply forces to jk pair
				
		SYSTEM.ACCEL[fidx_a2].X += force_4b[3] * RAB[3].X;
		SYSTEM.ACCEL[fidx_a2].Y += force_4b[3] * RAB[3].Y;
		SYSTEM.ACCEL[fidx_a2].Z += force_4b[3] * RAB[3].Z;

		SYSTEM.ACCEL[fidx_a3].X -= force_4b[3] * RAB[3].X;
		SYSTEM.ACCEL[fidx_a3].Y -= force_4b[3] * RAB[3].Y;
		SYSTEM.ACCEL[fidx_a3].Z -= force_4b[3] * RAB[3].Z;					
				
		// Apply forces to jl pair

		SYSTEM.ACCEL[fidx_a2].X += force_4b[4] * RAB[4].X;
		SYSTEM.ACCEL[fidx_a2].Y += force_4b[4] * RAB[4].Y;
		SYSTEM.ACCEL[fidx_a2].Z += force_4b[4] * RAB[4].Z;

		SYSTEM.ACCEL[fidx_a4].X -= force_4b[4] * RAB[4].X;
		SYSTEM.ACCEL[fidx_a4].Y -= force_4b[4] * RAB[4].Y;
		SYSTEM.ACCEL[fidx_a4].Z -= force_4b[4] * RAB[4].Z;		
											
		// Apply forces to kl pair
				
		SYSTEM.ACCEL[fidx_a3].X += force_4b[5] * RAB[5].X;
		SYSTEM.ACCEL[fidx_a3].Y += force_4b[5] * RAB[5].Y;
		SYSTEM.ACCEL[fidx_a3].Z += force_4b[5] * RAB[5].Z;

		SYSTEM.ACCEL[fidx_a4].X -= force_4b[5] * RAB[5].X;
		SYSTEM.ACCEL[fidx_a4].Y -= force_4b[5] * RAB[5].Y;
		SYSTEM.ACCEL[fidx_a4].Z -= force_4b[5] * RAB[5].Z;					

	 }	
  }		
}	// If 4-body intereaction.

////////////////////////////////////////////////////////////
//
// FUNCTIONS -- PRINTING OF POTENTIAL ENERGY SURFACE
//
////////////////////////////////////////////////////////////

void Cheby::Print_2B(int ij, string PAIR_NAME, bool INCLUDE_FCUT, bool INCLUDE_CHARGES, bool INCLUDE_PENALTY, string FILE_TAG)
// Generating pair distance scans for the 2-b potential.
// pair distances will range from smin to smax, incremented by sdelta
{

	double rlen;
	static double *Tn, *Tnd;
	static bool called_before = false;
	
	double rpenalty;
	double coeff;			  
	double fcut0; 
	double fcut; 				
	double Vpenalty;

	string OUTFILE = "2b_Cheby_Pot-";
	OUTFILE.append(PAIR_NAME);
	
	if(FILE_TAG != "")
		OUTFILE.append("_for_3B");
	
	OUTFILE.append(".dat");
	ofstream OUTFILE_2B_POT;
	OUTFILE_2B_POT.open(OUTFILE.data());	
	
	SCAN_FILE_2B = OUTFILE;
	
	// A penalty function is added to the potential for r + penalty_dist < smin[ipair]
	// All pairs have the same penalty scale and distance
	const double penalty_scale  = FF_2BODY[0].PENALTY_SCALE;	// 1.0e8;
	const double penalty_dist   = FF_2BODY[0].PENALTY_DIST;  	// 0.01;

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;
		
		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			if (FF_2BODY[i].SNUM > dim ) 
				dim = FF_2BODY[i].SNUM;	 
		
		dim++;
		Tn   = new double [dim];
		Tnd   = new double [dim];
		
	}
  
	// Main loop for Chebyshev terms:
	
	int n_ij = (FF_2BODY[ij].S_MAXIM - FF_2BODY[ij].S_MINIM)/FF_2BODY[ij].S_DELTA;
	double tempx;

	for (int a=1; a<n_ij; a++)
	{
		tempx = 0;
		
		rlen = FF_2BODY[ij].S_MINIM + a * FF_2BODY[ij].S_DELTA;

		if(rlen > FF_2BODY[ij].S_MINIM and rlen < FF_2BODY[ij].S_MAXIM)	
		{
			// Apply a penalty for distances less than smin - penalty_dist.
			
			rpenalty = 0.0;
			
			if ( rlen - penalty_dist < FF_2BODY[ij].S_MINIM ) 
				rpenalty = FF_2BODY[ij].S_MINIM + penalty_dist - rlen;

			// Calculate Cheby polys.

			set_polys(ij, Tn, Tnd, rlen, FF_2BODY[ij].X_DIFF,
						 FF_2BODY[ij].X_AVG, FF_2BODY[ij].SNUM) ; 
	
			// Now compute the force/potential... Coulomb 
			// Ewald sum not required for a scan-type ("cluster") calculation
			// ...If charges are zero, nothing will be added to tempx
			
			if(INCLUDE_CHARGES)
				tempx += FF_2BODY[ij].ATM1CHG * FF_2BODY[ij].ATM2CHG / rlen;

			// Now compute the force/potential... Cheby

			fcut0 = (1.0 - rlen/FF_2BODY[ij].S_MAXIM);
			fcut      = pow(fcut0, FF_2BODY[ij].FORCE_CUTOFF.POWER );
									 
			for ( int i = 0; i < FF_2BODY[ij].SNUM; i++ ) 
			{
				coeff = FF_2BODY[ij].PARAMS[i]; // This is the Cheby FF paramfor the given power

				if(INCLUDE_FCUT)
					tempx += coeff * fcut * Tn[i+1];
				else
					tempx += coeff * Tn[i+1];
			}
			// Add penalty for very short distances, where the fit FF may be unphysical (preserve conservation of E).

			if(INCLUDE_PENALTY)
			{
				if ( rpenalty > 0.0 ) 
				{
					Vpenalty = 0.0;
					Vpenalty = rpenalty * rpenalty * rpenalty * penalty_scale;
					tempx += Vpenalty;
				}
			}



			OUTFILE_2B_POT << rlen << " " << tempx << endl;						
		} 

	}
	
	OUTFILE_2B_POT.close();
	
	return;
}  

// NEW
void Cheby::Print_3B(CLUSTER_LIST &TRIPS, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, 
							int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan)	
	// Print heat map slices for 2+3-body cheby potentials
{
	
  // Variables exclusive to 2-body

  double x_diff2b, x_avg2b ;
  static double *Tn, *Tnd;
  static double fcut_2b;
	
  // 3-body variables that may be shared with 2-body 
	 		
  double rlen_ij,  rlen_ik,  rlen_jk;
	
  static double *Tn_ij,  *Tn_ik,  *Tn_jk;
  static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
  static bool    called_before = false;
	
  static int pow_ij, pow_ik, pow_jk;
			  
  double fcut_ij,  fcut_ik,  fcut_jk; 			
	
  static string TEMP_STR;
  static int curr_triple_type_index;
  int curr_pair_type_idx_ij;
  int curr_pair_type_idx_ik;
  int curr_pair_type_idx_jk;
  vector<int> atom_type_idx(3) ;

	
  double S_MAXIM_IJ, S_MAXIM_IK, S_MAXIM_JK;
  double S_MINIM_IJ, S_MINIM_IK, S_MINIM_JK;
	
  vector<TRIP_FF> & FF_3BODY = TRIPS.VEC ;
  map<string,int> & TRIAD_MAP = TRIPS.MAP ;
  vector<int> pair_index(3) ;
  vector<double> x_diff(3), x_avg(3) ;

  ////////////////////////////////////////////////////////////////////////////////////////
	

  // A penalty function is added to the potential for r + penalty_dist < smin[ipair]
  // All pairs have the same penalty scale and distance
	 
  if ( ! called_before ) 
  {
	 called_before = true;
		
	 // Set up 2-body polynomials
		
	 int dim = 0;
		
	 for ( int i = 0; i < FF_2BODY.size(); i++ ) 
		if (FF_2BODY[i].SNUM > dim ) 
		  dim = FF_2BODY[i].SNUM;	 
		
	 dim++;
	 Tn   = new double [dim];
	 Tnd  = new double [dim];
		
	 // Set up 3-body polynomials
		 
	 dim = 0;

	 for ( int i = 0; i < FF_2BODY.size(); i++ ) 
		if (FF_2BODY[i].SNUM_3B_CHEBY > dim ) 
		  dim = FF_2BODY[i].SNUM_3B_CHEBY;	
			
	 dim++;

	 Tn_ij   = new double [dim];
	 Tn_ik   = new double [dim];
	 Tn_jk   = new double [dim];

	 Tnd_ij  = new double [dim];
	 Tnd_ik  = new double [dim];
	 Tnd_jk  = new double [dim]; 
	
  }
	
  double TMP_DOUB;

  // Determine the FF type for the given triplet

  TEMP_STR =      FF_2BODY[ij].PRPR_NM;
  TEMP_STR.append(FF_2BODY[ik].PRPR_NM);	
  TEMP_STR.append(FF_2BODY[jk].PRPR_NM);	

  curr_pair_type_idx_ij = ij;
  curr_pair_type_idx_ik = ik;
  curr_pair_type_idx_jk = jk;

  // Set the 3-body type
	
  curr_triple_type_index =  TRIAD_MAP[TEMP_STR];

  atom_type_idx[0] = FF_3BODY[curr_triple_type_index].match_atom_type_idx(ATM_TYP_1) ;
  atom_type_idx[1] = FF_3BODY[curr_triple_type_index].match_atom_type_idx(ATM_TYP_2) ;
  atom_type_idx[2] = FF_3BODY[curr_triple_type_index].match_atom_type_idx(ATM_TYP_3) ;

  map_indices_int(FF_3BODY[curr_triple_type_index], atom_type_idx, pair_index) ;

  S_MAXIM_IJ = FF_3BODY[curr_triple_type_index].S_MAXIM[pair_index[0]] ;
  S_MAXIM_IK = FF_3BODY[curr_triple_type_index].S_MAXIM[pair_index[1]] ;
  S_MAXIM_JK = FF_3BODY[curr_triple_type_index].S_MAXIM[pair_index[2]] ;
				
  S_MINIM_IJ = FF_3BODY[curr_triple_type_index].S_MINIM[pair_index[0]] ;
  S_MINIM_IK = FF_3BODY[curr_triple_type_index].S_MINIM[pair_index[1]] ;
  S_MINIM_JK = FF_3BODY[curr_triple_type_index].S_MINIM[pair_index[2]] ;

  // Set the number of slices and steps
	
  int SLICES				= 10;
  int STEPS				= 1000;
	
  double POT_ENER;

	
  cout << "NOTE!!! AUTOMATICALLY INCLUDING 2B PENALTY FUNCTIONS IF 3B TYPE IS CUBIC!!!" << endl;
		
  // ASSUMPTIONS: REQUESTED RMIN RMAX VALUES ALL WITHIN **** 3-BODY'S **** FF RANGES!
	
  for(int i=0; i<=SLICES; i++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev -1)
  {	
	 // Build the output file
		
	 string OUTFILE = "Ternary-3b_Cheby_Pot-";
	 OUTFILE.append(TEMP_STR);
	 OUTFILE.append("_scan-");
	 OUTFILE.append(to_string(static_cast<long long>(i+1)));
	 OUTFILE.append(".dat");
		
	 ofstream OUTFILE_3B_POT;
	 OUTFILE_3B_POT.open(OUTFILE.data());
				
	 // Set the ij distance 
		
	 rlen_ij = S_MINIM_IJ + i * 1.0; //FF_2BODY[ij].S_DELTA*50;
		
	 OUTFILE_3B_POT << "# ij dist = " << fixed << setprecision(2) << rlen_ij << endl;
	
	 if (rlen_ij > FF_2BODY[ij].S_MAXIM)
		break;
		
		
	 for(int j=0; j<STEPS; j++)	
	 {			
		rlen_ik = S_MINIM_IK + j * FF_2BODY[ik].S_DELTA;
			
		if (rlen_ik > FF_2BODY[ik].S_MAXIM)
		  break;
			
		for(int k=0; k<STEPS; k++)
		{			
				
		  rlen_jk = S_MINIM_JK + k * FF_2BODY[jk].S_DELTA;	
				
		  if (rlen_jk > FF_2BODY[jk].S_MAXIM)
			 break;
				
		  // Initialize the potential energy
		
		  POT_ENER = 0;
				
		  /////////////////////////////////////////////
		  // EVALUATE THE 2-BODY INTERACTIONS
		  /////////////////////////////////////////////
						 
		  // ij		 
						 
		  x_diff2b = FF_2BODY[curr_pair_type_idx_ij].X_DIFF ;
		  x_avg2b = FF_2BODY[curr_pair_type_idx_ij].X_AVG ;
		  set_polys(curr_pair_type_idx_ij, Tn, Tnd, rlen_ij, x_diff2b, x_avg2b, FF_2BODY[curr_pair_type_idx_ij].SNUM);
			
		  FF_2BODY[curr_pair_type_idx_ij].FORCE_CUTOFF.get_fcut(fcut_2b, TMP_DOUB, rlen_ij, 0.0, 
																				  FF_2BODY[curr_pair_type_idx_ij].S_MAXIM);

		  for ( int m = 0; m < FF_2BODY[curr_pair_type_idx_ij].SNUM; m++ ) 
			 POT_ENER += FF_2BODY[curr_pair_type_idx_ij].PARAMS[m] * fcut_2b * Tn[m+1];
				
		  if(rlen_ij - FF_2BODY[0].PENALTY_DIST < FF_2BODY[curr_pair_type_idx_ij].S_MINIM )
			 POT_ENER += pow((FF_2BODY[curr_pair_type_idx_ij].S_MINIM + FF_2BODY[0].PENALTY_DIST - rlen_ij),3.0)*FF_2BODY[0].PENALTY_SCALE;
				
		  // ik

		  x_diff2b = FF_2BODY[curr_pair_type_idx_ik].X_DIFF ;
		  x_avg2b = FF_2BODY[curr_pair_type_idx_ik].X_AVG ;
		  set_polys(curr_pair_type_idx_ik, Tn, Tnd, rlen_ik, x_diff2b, x_avg2b, FF_2BODY[curr_pair_type_idx_ik].SNUM);
			
		  FF_2BODY[curr_pair_type_idx_ik].FORCE_CUTOFF.get_fcut(fcut_2b, TMP_DOUB, rlen_ik, 0.0, FF_2BODY[curr_pair_type_idx_ik].S_MAXIM);

		  for ( int m = 0; m < FF_2BODY[curr_pair_type_idx_ik].SNUM; m++ ) 
			 POT_ENER += FF_2BODY[curr_pair_type_idx_ik].PARAMS[m] * fcut_2b * Tn[m+1];
				
		  if(rlen_ik - FF_2BODY[0].PENALTY_DIST < FF_2BODY[curr_pair_type_idx_ik].S_MINIM )
			 POT_ENER += pow((FF_2BODY[curr_pair_type_idx_ik].S_MINIM + FF_2BODY[0].PENALTY_DIST - rlen_ik),3.0)*FF_2BODY[0].PENALTY_SCALE;
				
		  // jk

		  x_diff2b = FF_2BODY[curr_pair_type_idx_jk].X_DIFF ;
		  x_avg2b = FF_2BODY[curr_pair_type_idx_jk].X_AVG ;				
		  set_polys(curr_pair_type_idx_jk, Tn, Tnd, rlen_jk, x_diff2b, x_avg2b, FF_2BODY[curr_pair_type_idx_jk].SNUM);
			
		  FF_2BODY[curr_pair_type_idx_jk].FORCE_CUTOFF.get_fcut(fcut_2b, TMP_DOUB, rlen_jk, 0.0, FF_2BODY[curr_pair_type_idx_jk].S_MAXIM);

		  for ( int m = 0; m < FF_2BODY[curr_pair_type_idx_jk].SNUM; m++ ) 
			 POT_ENER += FF_2BODY[curr_pair_type_idx_jk].PARAMS[m] * fcut_2b * Tn[m+1];
				
		  if(rlen_jk - FF_2BODY[0].PENALTY_DIST < FF_2BODY[curr_pair_type_idx_jk].S_MINIM )
			 POT_ENER += pow((FF_2BODY[curr_pair_type_idx_jk].S_MINIM + FF_2BODY[0].PENALTY_DIST - rlen_jk),3.0)*FF_2BODY[0].PENALTY_SCALE;
				
		  /////////////////////////////////////////////
		  // EVALUATE THE 3-BODY INTERACTIONS
		  /////////////////////////////////////////////

		  // Set up the polynomials

		  for ( int m = 0 ; m < 3 ; m++ ) 
		  {
			 x_avg[m] = FF_3BODY[curr_triple_type_index].X_AVG[pair_index[m]] ;
			 x_diff[m] = FF_3BODY[curr_triple_type_index].X_DIFF[pair_index[m]] ;
		  }

		  set_polys(curr_pair_type_idx_ij, Tn_ij, Tnd_ij, rlen_ij, x_diff[0], x_avg[0], FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY);
		  set_polys(curr_pair_type_idx_ik, Tn_ik, Tnd_ik, rlen_ik, x_diff[1], x_avg[1], FF_2BODY[curr_pair_type_idx_ik].SNUM_3B_CHEBY);
		  set_polys(curr_pair_type_idx_jk, Tn_jk, Tnd_jk, rlen_jk, x_diff[2], x_avg[2], FF_2BODY[curr_pair_type_idx_jk].SNUM_3B_CHEBY);
								
		  // Set up the penalty functions

		  FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_ij, TMP_DOUB, rlen_ij, S_MINIM_IJ, S_MAXIM_IJ);
		  FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_ik, TMP_DOUB, rlen_ik, S_MINIM_IK, S_MAXIM_IK);
		  FF_3BODY[curr_triple_type_index].FORCE_CUTOFF.get_fcut(fcut_jk, TMP_DOUB, rlen_jk, S_MINIM_JK, S_MAXIM_JK);
			
		  for(int i=0; i<FF_3BODY[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
		  {
			 set_3b_powers(FF_3BODY[curr_triple_type_index], pair_index, i,
								pow_ij, pow_ik, pow_jk) ;
			 if( rlen_ij<S_MAXIM_IJ && rlen_ik<S_MAXIM_IK && rlen_jk<S_MAXIM_JK )
				POT_ENER += FF_3BODY[curr_triple_type_index].PARAMS[i] * fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk]; 	
		  }
			
		  OUTFILE_3B_POT <<  rlen_ik << "     " <<  rlen_jk << "     " <<  POT_ENER << endl;	
				
		} 
			
		OUTFILE_3B_POT << endl;
	
	 }
		
	 OUTFILE_3B_POT.close();
  }
	
  return;
} 
