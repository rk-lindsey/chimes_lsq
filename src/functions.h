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
	
	int    PAIRIDX;			// Index for the specific pair type
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
	
	double LAMBDA;			// FUNCTIONALITY HASN'T BEEN CODED YET, BUT NEEDS TO BE
	double MIN_FOUND_DIST;	// FUNCTIONALITY HASN'T BEEN CODED YET, BUT NEEDS TO BE -- Minimum distance between pairs
	
	bool   USE_OVRPRMS;		// Should overbonding even be computed pair type
	string OVER_TO_ATM;		// Which atom is overbonding *to* being defined for... for example, overbonding to oxygen
	vector<double> OVRPRMS;	// [0] = P_OVERB; [1] = R_0_VAL; [2] = P_1_VAL; [3] = P_2_VAL; [4] = LAMBDA6
	
	PAIRS():OVRPRMS(5){}	// Just a constructor to allow the size of the OVRPRMS vector to be pre-specified
	
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
void ZCalc_Deriv (vector<PAIRS> & ATOM_PAIRS, FRAME & FRAME_TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<vector <XYZ > > & FRAME_A_MATRIX, vector<vector <XYZ > > & FRAME_COULOMB_FORCES, const int nlayers, bool if_3b_cheby, map<string,int> & PAIR_MAP);

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

