#ifndef _HELPERS_
#define _HELPERS_

#include<stdio.h>
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<vector>
#include<string.h>
#include <map>

#include "functions.h"

using namespace std;

static const double ke      = 332.0637157615209;//converter between electron units and Stillinger units for Charge*Charge.
static const double Hartree = 627.50961; // 1 Hartree in kcal/mol.
static const double Kb      = 0.001987; // Boltzmann constant in kcal/mol-K.
static const double Tfs     = 48.888;   // Internal time unit in fs.
static const double GPa     = 6.9479;



/*
enum Sr_pair_t 
{
  CHEBYSHEV, DFTBPOLY, SPLINE, INVERSE_R, LJ, and STILLINGER
};
*/

struct PAIRS
{
	
	int    PAIRIDX;			// Index for the specific pair type (triplet type)
	string PRPR_NM;			// The order-specific atom pair. For example, OH, never HO, or something similiar "proper name"
	string PAIRTYP;			// Allowed values are CHEBYSHEV, DFTBPOLY, SPLINE, INVERSE_R, LJ, and STILLINGER
	string ATM1TYP;			// Atom chemistry (i.e. C, H, N, O...)
	string ATM2TYP;
	double ATM1CHG;			// Atom partial charge... used when charges are fixed
	double ATM2CHG;
	double ATM1MAS;			// Atomic mass (i.e. 12.011 for C)
	double ATM2MAS;	
	double S_MINIM;			// Minimum allowed pair distance for fitting
	double S_MAXIM;			// Maximum allowed pair distance for fitting
	double S_DELTA;			// Fitting "grid" spacing (width)

	int    SNUM;			// Number of fitting parameters for pair ... WHY WOULD THIS BE DIFFERENT FOR DIFFERENT PAIR TYPES?***
	int    SNUM_3B_CHEBY;	// Number of fitting parameters for pair ... WHY WOULD YOU NEED BOTH SNUM AND THIS SPECIAL CHEBY ONE?***
	
	string CHEBY_TYPE;		// Are distances transformed into inverse-r type or morse-type distances?... or not at all? (default)
	
	double LAMBDA;			
	double MIN_FOUND_DIST;	//  Minimum distance between pairs
	
	bool   USE_OVRPRMS;		// Should overbonding even be computed pair type
	string OVER_TO_ATM;		// Which atom is overbonding *to* being defined for... for example, overbonding to oxygen
	vector<double> OVRPRMS;	// [0] = P_OVERB; [1] = R_0_VAL; [2] = P_1_VAL; [3] = P_2_VAL; [4] = LAMBDA6
	
	PAIRS():OVRPRMS(5){}	// Just a constructor to allow the size of the OVRPRMS vector to be pre-specified
	
	// SPECIAL CASE: When 3B cheby potential is being used.
	//					Use as example atom types O and H
	//					Would have pairs { [O,O], [O,H], [H,H] } 
	// 						.. ATM1TYP for pair 2 would be O, ATM2TYP for pair 2 would be H
	//					Would have triplets: { [OO,OO,OO], [OO,OO,HH], ...	} 
	// 						.. ATM1TYP for pair 2 would be OO, ATM2TYP for pair 2 would be OO, ATM3TYP would be HH
	//					Don't care about charges or mass for triplet's, or several other variables.
	
	// *** IF THESE ARE ALWAYS THE SAME FOR EACH PAIR TYPE, WE SHOULD CREATE A DATA STRUCTURE CALLED PAIR TYPE WHICH WOULD STORE INFO LIKE
	//     THE TYPE (I.E. SPLINE, CHEBY...), AND INFO ON THE NUMBER OF PARAMETERS
};

struct XYZ
{
	double X;
	double Y;
	double Z;
};

struct XYZ_INT
{
	int X;
	int Y;
	int Z;
};

struct TRIPLETS
{
	int    TRIPINDX;
	string ATMPAIR1;
	string ATMPAIR2;
	string ATMPAIR3;
	
	int N_TRUE_ALLOWED_POWERS;	// How many UNIQUE sets of powers do we have?
	int N_ALLOWED_POWERS;	// How many total sets of powers do we have?
	
	vector<XYZ_INT> ALLOWED_POWERS;	// This will keep a list of the allowed polynomial powers for each coefficient
	vector<int>		EQUIV_INDICIES;	// For each set of allowed powers, what is the index of the first equivalent set? For example, for the set (OO, OH, OH), (1,0,1) and (1,1,0) are is equivalent
	vector<int>		PARAM_INDICIES;	// For each of the set of allowed powers, what would be the index in the FF? for example, for a set of EQUIV_INDICIES {0,0,2,3}, PARAM_INDICIES would be {0, 0, 1, 2}
//	vector<int>		POWER_DUPLICATES;	// This will tell the multiplicity of each set of powers. For example, for the set (OO, OH, OH), (1,0,1) has a multiplicity of 2, since (1,1,0 is equivalent) -- No longer used
};


struct FRAME
{
	int ATOMS;
	XYZ BOXDIM;
	
	vector<string> ATOMTYPE;
	vector<XYZ>    COORDS;
	vector<double> CHARGES;
	vector<XYZ>    FORCES;
};


// FUNCTION UPDATED
void ZCalc_Deriv (vector<PAIRS> & ATOM_PAIRS,  vector<TRIPLETS> PAIR_TRIPLETS, FRAME & FRAME_TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<vector <XYZ > > & FRAME_A_MATRIX, vector<vector <XYZ > > & FRAME_COULOMB_FORCES, const int nlayers, bool if_3b_cheby, map<string,int> & PAIR_MAP,  map<string,int> & TRIAD_MAP);

// FUNCTION UPDATED
void SubtractCoordForces (FRAME & TRAJECTORY, bool calc_deriv, vector<XYZ> & P_OVER_FORCES,  vector<PAIRS> & ATOM_PAIRS, vector<int> & ATOM_PAIR_TYPES_ALL, map<string,int> & PAIR_MAP);

// FUNCTION UPDATED
void ZCalc_Ewald (FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP);

// FUNCTION UPDATED
void optimal_ewald_params(double accuracy, double V, int nat, double &alpha, double & rc, int & kc, double & r_acc, double & k_acc);

// FUNCTION UPDATED
void ZCalc_Ewald_Deriv (FRAME & FRAME_TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector <vector <XYZ > > & FRAME_COULOMB_FORCES, map<string,int> & PAIR_MAP);

// FUNCTION UPDATED
void ZCalc_3B_Cheby_Deriv(FRAME & TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers);


#endif

