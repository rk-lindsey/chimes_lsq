#ifndef VERBOSITY
	#define VERBOSITY 1
#endif

#ifndef FORCECHECK
	#define FORCECHECK 0
#endif

#define FULL	// distance transformation such that pair dist falls on range -1 to 1 
#define LEFT	// distance transformation such that pair dist falls on range -1 to 0
#define RIGHT	// distance transformation such that pair dist falls on range  0 to 1 

#ifndef XFORMTYPE
	#define FULL
#endif 

#define GNUPLOT 1
#define MATLAB  2
#define PYTHON  3

#ifndef PESFORMAT
	#define PESFORMAT GNUPLOT
#endif

#ifndef _HELPERS_
#define _HELPERS_

#include<stdio.h>
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<assert.h>
#include<map>


using namespace std;

static const double ke      = 332.0637157615209;	// Converter between electron units and Stillinger units for Charge*Charge.
static const double Hartree = 627.50961; 			// 1 Hartree in kcal/mol.
static const double Kb      = 0.001987; 			// Boltzmann constant in kcal/mol-K.
static const double Tfs     = 48.888;   			// Internal time unit in fs.
static const double GPa     = 6.9479;				// Unit conversion factor to GPa.

extern string FULL_FILE_3B;	// The 4D PES for 3B FF
extern string SCAN_FILE_3B;	// The 2D PES scans for 3B
extern string SCAN_FILE_2B;	// The 2D PES scans for 2B

struct MD_JOB_CONTROL
{
	int STEP;					// Tracks the current md step
	
	///////////////////////////////////////////////
	// Variables read in from the parameter file:
	///////////////////////////////////////////////
	
	bool   FIT_COUL;			// Replaces fit_coul... If true, take charges from spline parameters.
	bool   USE_COULOMB;	 		// Replaces if_coulomb... If true, calculate Coulomb forces.
    bool   USE_OVERCOORD;		// Replaces if_overcoord... If true, calculate ReaxFF-like overcoordination term.
    bool   FIT_POVER;			// Replaces fit_pover... If true, find linear overcoordination parameter from least-squares fitting. -- this needs to be updated for new handling
    bool   USE_3B_CHEBY;		// Replaces if_3b_cheby... If true, calculate 3-Body Chebyshev interaction.
    
	
	///////////////////////////////////////////////
	// Variables read in from main md code input file
	///////////////////////////////////////////////
	
	// "General control variables"

	int    SEED;				// Replaces rand_seed...
	double TEMPERATURE;			// Replaces TempMD
	bool   PLOT_PES;			// Plot out the potential energy surfaces? IF so, nothing else will be done.
	bool   COMPARE_FORCE;		// Replaces if_read_force... If TRUE, read in read in a set of forces from a file for comparison with those computed by this code
	string COMPARE_FILE;		// Name of the file that holds the forces
	double DELTA_T;				// Relaces deltat
	double DELTA_T_FS;			// Replaces deltat_fs... in femtoseconds
	int	   N_MD_STEPS;			// Replaces nsteps...
	int    N_LAYERS;			// Replaces nlayers... Number of periodic images to use. Usually, no need to go past 1 ...(depends only on spline interaction cutoff).
    string PARAM_FILE;			// Replaces params_file...
    string COORD_FILE;			// Replaces xyz_file...
	bool   SELF_CONSIST;		// Is this part of a self-consistent DFT MD --> FIT --> MM MD --> CYCLE type calculation?
	double NVT_CONV_CUT;		// What is the cutoff for "conservation"... Will kill the program if the current temperature and set temperature differ by more than this fraction
								// ie if | T_set - T_curr | >  NVT_CONV_CUT * T_set, end with error message.. defaults to 0.1 (10%)

	// "Simulation options"
	
	bool   INIT_VEL;				// Replaces if_init_vel... If true, initialize velocities.
	bool   USE_HOOVER_THRMOSTAT;	// Replaces if_hoover... Use a Nose-Hoover thermostat?	
	double FREQ_UPDATE_THERMOSTAT;	// Replaces scale_freq and thoover_fs... it's usage depends on whether USE_HOOVER_THERMOSTAT is true or false.. will be cast as int where required
	bool   USE_NUMERICAL_PRESS;	// Replaces num_pressure... Whether to calculate pressures by finite difference.
	 
	// "Output control" 

	int    FREQ_DFTB_GEN;		// Replaces gen_freq... How often to write the gen file.
	bool   PRINT_VELOC;			// If true, write out the velocities 
	int    FREQ_VELOC;
	int    FREQ_ENER;			// Replaces energy_freq... How often to output energy
	bool   PRINT_FORCE;			// Replaces if_output_force... If TRUE, write out calculated forces.	
	int    FREQ_FORCE;			// How often to print the forces	
	int    SELF_CONSIST_FREQ;	// How frequently to print POSCAR file
	bool   WRAP_COORDS;			// Should coordinates be wrapped?
};

struct NOSE_HOOVER
{
	double VISCO;		// Replaces zeta...	Hoover viscosity
	double VISCO1;		// Replaces zeta1... not sure what this is/does
	double VISCO_DOT0;	// Replaces zeta_dot0.. not sure what this is/does
	double VISCO_DOT1;	// Replaces zeta_dot1.. not sure what this is/does
	double COORD;		// Replaces s_hoover...	Hoover coordinate
	double TIME;		// Replaces thoover and thoover_fs...	In fs
	double N_DOF;		// Replaces ndof...	# degrees of freedom
	double CHRG;		// Replaces QHoover... not sure what this is/does
	double VSCALEH;		// Replaces vscaleh..  not sure what this is/does


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
	
	double 	TEMPERATURE;		// This is the RUNNING temperature, not the set temperature!
	double 	PRESSURE;			// This is the RUNNING pressure, not the set pressure!
	double 	PRESSURE_XYZ;		// This is pxyz from the code.. not sure what it specifically corresponds to...
	XYZ		PRESSURE_TENSORS;	// These are the RUNNING pressure tensors, no the set pressure tensors!
	double	TOT_POT_ENER;		// Replaces VTOT
	
	vector<string> 	ATOMTYPE;
	vector<int> 	ATOMTYPE_IDX;	// Only used for dftbgen file printing
	vector<XYZ>		COORDS;
	vector<double> 	CHARGES;
	vector<double> 	MASS;
	vector<XYZ>		FORCES;
	vector<XYZ>		ACCEL;
	vector<XYZ>		TMP_EWALD;	// Holds temporary ewald accels/forces
	vector<XYZ>		VELOCITY;
	vector<XYZ>		VELOCITY_NEW;
	
};

struct PAIRS	// NEEDS UPDATING
{
	
	int    PAIRIDX;			// Index for the specific pair type (triplet type)
	string PRPR_NM;			// The order-specific atom pair. For example, OH, never HO, or something similiar "proper name"
	string PAIRTYP;			// Allowed values are CHEBYSHEV, DFTBPOLY, SPLINE, INVERSE_R, LJ, and STILLINGER
	string ATM1TYP;			// Atom chemistry (i.e. C, H, N, O...)
	string ATM2TYP;
	double ATM1CHG;			// Atom partial charge... used when charges are fixed
	double ATM2CHG;
	string CHRGSGN;			// Should the fitted charge on a given atom be negative or positive?
	double ATM1MAS;			// Atomic mass (i.e. 12.011 for C)
	double ATM2MAS;	
	double S_MINIM;			// Minimum allowed pair distance for fitting
	double S_MAXIM;			// Maximum allowed pair distance for fitting
	double S_DELTA;			// Fitting "grid" spacing (width)

	int    SNUM;			// Number of fitting parameters for pair ... WHY WOULD THIS BE DIFFERENT FOR DIFFERENT PAIR TYPES?***
	int    SNUM_3B_CHEBY;	// Number of fitting parameters for pair ... WHY WOULD YOU NEED BOTH SNUM AND THIS SPECIAL CHEBY ONE?***
	
	string CHEBY_TYPE;		// Are distances transformed into inverse-r type or morse-type distances?... or not at all? (default)
	double PENALTY_SCALE;	// For 2B Cheby potentials... "a" in vpenalty = a*(smin-penalty_dist-rlen)^3 ... default value is 1.0e8
	double PENALTY_DIST;	// For 2B Cheby potentials... "penalty_dist" in vpenalty = a*(smin-penalty_dist-rlen)^3 ... default value is 0.01
	double CUBIC_SCALE;		// Factor to multiply to the cubic penalty function, (1-rlen/smax)^3... default value is 1
	
	double LAMBDA;			// Morse lambda for CHEBYSHEV type pairs
	double MIN_FOUND_DIST;	// Minimum distance between pairs
	
	bool   USE_OVRPRMS;		// Should overbonding even be computed pair type
	string OVER_TO_ATM;		// Which atom is overbonding *to* being defined for... for example, overbonding to oxygen
	vector<double> OVRPRMS;	// [0] = P_OVERB; [1] = R_0_VAL; [2] = P_1_VAL; [3] = P_2_VAL; [4] = LAMBDA6
	
	PAIRS():OVRPRMS(5){}	// Just a constructor to allow the size of the OVRPRMS vector to be pre-specified
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

struct PAIR_FF : public PAIRS
{
	vector<double> 	PARAMS;
	vector<double>	POT_PARAMS;		// Used by splines to compute pressure by integrating spline eq's
	double			PAIR_CHRG;
};

struct TRIP_FF : public TRIPLETS
{
	vector<double> 	PARAMS;
};


struct PES_PLOTS
{
	int N_PLOTS;
	vector<string> PES_TYPES;
	vector<int> NBODY;
	vector<int> TYPE_INDEX;
	string TYPE_1;
	string TYPE_2;
	string TYPE_3;
	bool INCLUDE_2B;
	
	// Variables used for scans of 3b potential
	
	int N_SCAN;
	
	vector<int>PARENT_TYPE;
	
	vector<int> FIX_PAIR_1;
	vector<int> FIX_PAIR_2;
	
	vector<double> FIX_VAL_1;
	vector<double> FIX_VAL_2;
	
	vector<int> SCAN_PAIR;
	
	vector<string> SCAN_TYPE;
	
	// Variables used to add 2b to 3b
	
	vector<string> 	SEARCH_STRING_2B;
	vector<XYZ_INT>	IJ_IK_JK_TYPE;
		
	
	PES_PLOTS():N_PLOTS(0), N_SCAN(0), INCLUDE_2B(0){}
		
};

struct CHARGE_CONSTRAINT
{
	vector<string> PAIRTYPE;		// This will be read in from the input file
	vector<int>	   PAIRTYPE_IDX;	// This will be set using the map
	vector<double> CONSTRAINTS;
	double		   FORCE;
};

// For MPI calculations.  Number of processors and index of current processor.
extern int NPROCS ;
extern int RANK ;

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- DERIVATIVE CALCULATION
//
//////////////////////////////////////////



// FUNCTION UPDATED
void ZCalc_Deriv (vector<PAIRS> & ATOM_PAIRS,  vector<TRIPLETS> PAIR_TRIPLETS, FRAME & FRAME_TRAJECTORY, vector<vector <XYZ > > & FRAME_A_MATRIX, vector<vector <XYZ > > & FRAME_COULOMB_FORCES, const int nlayers, bool if_3b_cheby, map<string,int> & PAIR_MAP,  map<string,int> & TRIAD_MAP);

// FUNCTION UPDATED
void SubtractCoordForces (FRAME & TRAJECTORY, bool calc_deriv, vector<XYZ> & P_OVER_FORCES,  vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP);

// FUNCTION UPDATED -- OVERLOADING THE FUNCTION.. DIFFERENT INPUT REQUIRED DEPENDING ON WHETHER FUNCTION IS CALLED FROM MD
//                     PROGRAM OR LSQ FITTING PROGRAM...
void ZCalc_Ewald (FRAME & TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP);	// LSQ
void ZCalc_Ewald(FRAME & TRAJECTORY);	// MD

// FUNCTION UPDATED
void optimal_ewald_params(double accuracy, double V, int nat, double &alpha, double & rc, int & kc, double & r_acc, double & k_acc);

// FUNCTION UPDATED
void ZCalc_Ewald_Deriv (FRAME & FRAME_TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, vector <vector <XYZ > > & FRAME_COULOMB_FORCES, map<string,int> & PAIR_MAP);

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- FORCE CALCULATION
//
//////////////////////////////////////////

void ZCalc(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP);

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- PRINTING OF POTENTIAL ENERGY SURFACE
//
//////////////////////////////////////////

void Print_Cheby(vector<PAIR_FF> & FF_2BODY, int ij, string PAIR_NAME, string FILE_TAG);
void Print_3B_Cheby(MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk);
void Print_3B_Cheby_Scan(MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan);

void divide_atoms(int &a1start, int &a1end, int atoms) ;

#endif

