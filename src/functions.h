#ifndef VERBOSITY
	#define VERBOSITY 1
#endif

#ifndef FORCECHECK
	#define FORCECHECK 0
#endif

#ifndef FPENALTY_POWER
	#define FPENALTY_POWER 3.0
#endif

#ifndef WARN
	#define WARN TRUE
#endif 

#define GNUPLOT 1
#define MATLAB  2
#define PYTHON  3

#ifndef PESFORMAT
	#define PESFORMAT GNUPLOT
#endif

#ifndef CHECK_CHEBY_RANGE 
	#define CHECK_CHEBY_RANGE 1	// (true)
#endif

#ifndef _HELPERS_
#define _HELPERS_

// For "original" MD code

#include<stdio.h>
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<assert.h>
#include<map>

// For the LAMMPS version of the MD Code

#if defined(USE_MPI) && defined(LINK_LAMMPS)

	#define  QUOTE_(x) #x
	#define  QUOTE(x) QUOTE_(x)

	#include "lmppath.h"

	#include QUOTE(LMPPATH/src/library.h)
	#include QUOTE(LMPPATH/src/fix.h)
	#include QUOTE(LMPPATH/src/fix_external.h)  

	#include QUOTE(LMPPATH/src/lammps.h)
	#include QUOTE(LMPPATH/src/input.h)
	#include QUOTE(LMPPATH/src/modify.h)
#endif


using namespace std;

// Unit converters 

static const double ke      = 332.0637157615209;	// Converter between electron units and Stillinger units for Charge*Charge.
static const double Hartree = 627.50961; 			// 1 Hartree in kcal/mol.
static const double Kb      = 0.001987; 			// Boltzmann constant in kcal/mol-K.
static const double Tfs     = 48.888;   			// Internal time unit in fs.
static const double GPa     = 6.9479;				// Unit conversion factor... 1 kcal/mol/A^3 = (this constant) GPa
static const double atm     = GPa*9869.23266716;    // stress conversion to atm (for LAMMPS).
static const double GPa2atm = 9869.23266716;		// x_GPa * GPa2atm = x_atm
static const double pi		= 3.14159265359;		

// Global variables declared as externs in functions.h, and declared in functions.C -- general

extern string FULL_FILE_3B;	// The 4D PES for 3B FF
extern string SCAN_FILE_3B;	// The 2D PES scans for 3B
extern string SCAN_FILE_2B;	// The 2D PES scans for 2B

// Global variables declared as externs in functions.h, and declared in functions.C -- MPI calculations.   
 
extern int NPROCS;		// Number of processors
extern int RANK;		// Index of current processor

struct JOB_CONTROL
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
	double PRESSURE;			// Introduced for NPT
	string ENSEMBLE;			// NVE, NPT, NVT
	bool   PLOT_PES;			// Plot out the potential energy surfaces? IF so, nothing else will be done.
	bool   COMPARE_FORCE;		// Replaces if_read_force... If TRUE, read in read in a set of forces from a file for comparison with those computed by this code
	bool   SUBTRACT_FORCE;		// Read frame, compute forces based on parameter file, print out frame where FF forces have been subtracted from input forces
	string COMPARE_FILE;		// Name of the file that holds the forces
	double DELTA_T;				// Relaces deltat
	double DELTA_T_FS;			// Replaces deltat_fs... in femtoseconds
	int	   N_MD_STEPS;			// Replaces nsteps...
	int    N_LAYERS;			// Replaces nlayers... Number of periodic images to use. Usually, no need to go past 1 ...(depends only on spline interaction cutoff).
	double SCALE_SYSTEM_BY;		// Amount to scale boxlengths/coordinates by
    string PARAM_FILE;			// Replaces params_file...
    vector<string> COORD_FILE;	// Replaces xyz_file... Can be a list of files, which will be assembled together along the z-axis to form a single cell
	bool   SELF_CONSIST;		// Is this part of a self-consistent DFT MD --> FIT --> MM MD --> CYCLE type calculation?
	double NVT_CONV_CUT;		// What is the cutoff for "conservation"... Will kill the program if the current temperature and set temperature differ by more than this fraction
								// ie if | T_set - T_curr | >  NVT_CONV_CUT * T_set, end with error message.. defaults to 0.1 (10%)

	// "Simulation options"

	bool   INIT_VEL;				// Replaces if_init_vel... If true, initialize velocities.
	bool   USE_HOOVER_THRMOSTAT;	// Replaces if_hoover... Use a Nose-Hoover thermostat?	
	int    FREEZE_IDX_START;		// First atom in continuous range to freeze.. .counting starts from 0... value of -1 indicates no atoms are to be frozen
	int    FREEZE_IDX_STOP;			// Last atom in continuous range to freeze 
	double FREQ_UPDATE_THERMOSTAT;	// Replaces scale_freq and thoover_fs... it's usage depends on whether USE_HOOVER_THERMOSTAT is true or false.. will be cast as int where required
	double FREQ_UPDATE_BAROSTAT;	// Barostat time constant... defaults to 1000
	bool   USE_NUMERICAL_PRESS;		// Replaces num_pressure... Whether to calculate pressures by finite difference.
	 
	// "Output control" 

	int    FREQ_DFTB_GEN;		// Replaces gen_freq... How often to write the gen file.
	int    FREQ_BACKUP;       // How often to write backup files for restart.
	bool   PRINT_VELOC;			// If true, write out the velocities 
	int    FREQ_VELOC;
	int    FREQ_ENER;			// Replaces energy_freq... How often to output energy
	bool   PRINT_FORCE;			// Replaces if_output_force... If TRUE, write out calculated forces.	
	int    FREQ_FORCE;			// How often to print the forces	
	int    SELF_CONSIST_FREQ;	// How frequently to print POSCAR file
	bool   WRAP_COORDS;			// Should coordinates be wrapped?
	
	// Controls for how to construct the initial system, if desired
	
	bool   BUILD;
	string BUILD_TYPE;
	string BUILD_FILE;
	double BUILD_BOXL;
	int    BUILD_NMOLEC;
	string BUILD_ATOM;
		
	///////////////////////////////////////////////
	// Variables exclusive to the LSQ code
	///////////////////////////////////////////////
	
	bool IS_LSQ;				// Is this for an lsq run or actual md?
	bool FIT_STRESS;			// Should stress tensors be included in the fit?
	bool FIT_ENER;				// Should the total frame energy be included in the fit?
	bool CALL_EWALD;			// Should ewald subroutines be called?
	bool USE_POVER;			// Should overbonding information be printed to the header file?
	
	int		NFRAMES;			// Number of frames in the movie file
	int		NATMTYP;			// How many atom types are in the trajectory?
	int		CHEBY_ORDER;		// Order of Chebyshev polynomial if used... set to 8 for DFTB Erep polynomial
	int		CHEBY_3B_ORDER;		// how many polynomials for 3b cheby?
	int		NUM_3B_CHEBY;		// How many parameters are associated with cheby order CHEBY_3B_ORDER?
	int		INVR_PARAMS;		// currently uses 19 parameters per pair type
	int 	TOT_SNUM;			// total number of force field parameters
	int 	TOT_SHORT_RANGE;	// Number of short tranged FF params... i.e. not Ewald
	
	bool	COUL_CONSV;			// If true, constraints will be applied to charge fitting to try to maintain consistency
	bool	IF_SUBTRACT_COORD;	// If true, subtract overcoordination forces.
	bool	IF_SUBTRACT_COUL;	// If true, subtract Coulombic forces (for use with fixed charges).
	bool	USE_PARTIAL_CHARGES;// Will there be any charges in the system?

	string	CHEBY_TYPE;			// How will distance be transformed?
	string	INFILE;				// Input trajectory file
	
	JOB_CONTROL():N_LAYERS(0), WRAP_COORDS(false),IF_SUBTRACT_COORD(false),IF_SUBTRACT_COUL(false),FIT_COUL(false),USE_PARTIAL_CHARGES(false),COUL_CONSV(false),FIT_POVER(false),USE_3B_CHEBY(false),TOT_SNUM(0){}
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
    int ATOMS;                 		// Just the parent atoms.
    int ALL_ATOMS;          	   	// All atoms, including ghosts. 
	
	int MY_ATOMS;					// Used for lammps linking. Specify how many atoms in SYS the process owns
	int MY_ATOMS_START;				// Used for lammps linking. Specify what index along SYS starts the process' atoms
	XYZ BOXDIM;						// Dimenions of the primitive box.
	XYZ WRAPDIM;              		// Dimensions to use in wrapping coordinates for forces.  
	XYZ STRESS_TENSORS;
	
	double QM_POT_ENER;				// This is the potential energy of the QM calculation!
	
	double 	TEMPERATURE;			// This is the RUNNING temperature, not the set temperature!
	double 	PRESSURE;				// This is the RUNNING pressure, not the set pressure!
	double 	AVG_TEMPERATURE;		// Only used for velocity scaling-type thermostats
	double 	PRESSURE_XYZ;			// This is the running pressure sans the ideal gas term
	XYZ		PRESSURE_TENSORS_XYZ;	// These are the RUNNING pressure tensors, no the set pressure tensors! ...sans the ideal gas term
	XYZ		PRESSURE_TENSORS;		// Adds in the ideal gas term
	double	TOT_POT_ENER;			// Replaces VTOT
	
	vector<int>		PARENT;
	vector<XYZ_INT> LAYER_IDX;
	vector<string> 	ATOMTYPE;
	vector<int> 	ATOMTYPE_IDX;	// Only used for dftbgen and LAMMPS file printing
	vector<XYZ>		COORDS;
	vector<XYZ>    ALL_COORDS;  	// Coordinates of atoms + ghosts used for force evaluation.
	vector<double> 	CHARGES;
	vector<double> 	MASS;
	vector<XYZ>		FORCES;
	vector<XYZ>		ACCEL;
	vector<XYZ>		TMP_EWALD;		// Holds temporary ewald accels/forces
	vector<XYZ>		VELOCITY;
	vector<XYZ>		VELOCITY_NEW;
	vector<XYZ>		VELOCITY_ITER;
	
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
	double OLD_S_MINIM;		// Used for smoothing functions when 3B and 2B inner cutoffs differ. This should be the inner cutoffs used for the actual SVD
	double S_MINIM;			// Minimum allowed pair distance for fitting
	double S_MAXIM;			// Maximum allowed pair distance for fitting
	double S_DELTA;			// Fitting "grid" spacing (width)

	int N_CFG_CONTRIB;		// How many configurations actually contribute to fitting this pair??

	double CHEBY_RANGE_LOW;	//	When fitting to Cheby polynomials, pair distances are typically transformed to exist defined on range -1 to 1 (where Cheby poly's are defined), 
	double CHEBY_RANGE_HIGH;//  but under certain circumstances, it can be advantagous to only fit some sub range (i.e. -1 to 0). These variables define the transformation range
							//  and arr, by default, set to -1 and 1 for low and high, respectively.

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
	
	XYZ NBINS;				// Number of bins to use for ij, ik, and jk distances when building the 3B population histograms 
	
	PAIRS():OVRPRMS(5),N_CFG_CONTRIB(0){}	// Just a constructor to allow the size of the OVRPRMS vector to be pre-specified
};

struct TRIPLETS
{
	int    TRIPINDX;
	string ATMPAIR1;
	string ATMPAIR2;
	string ATMPAIR3;
	
	string FCUT_TYPE;		// "CUBIC" "COSINE" or "SIGMOID" currently supported
	double FCUT_STEEPNESS;
	double FCUT_OFFSET; 
	double FCUT_HEIGHT;
	
	int    N_CFG_CONTRIB;	// How many configurations actually contribute to fitting this triplet??
	
	XYZ MIN_FOUND;			// Testing two cases: 1. 
	
	XYZ S_MAXIM_3B;			// A unique outer cutoff for 3B interactions... by default, is set to S_MAXIM
	XYZ S_MINIM_3B;			// Similar for inner cutoff. This is useful when the 2-body 
							// is refit/extrapolated, thus has a s_min lower than the original fitted value
							// Values need to be specified for each contributing pair
							// X -> IJ, Y -> IK,  -> JK 
	
	vector<vector<vector< int > > > POP_HIST; // Population histogram that s used to set 3B behavior in unsampled regions
	 
	XYZ NBINS;				// Number of bins to use for ij, ik, and jk distances when building the population histograms 
	XYZ BINWS; 				// Binwidths to use for ij, ik, and jk distances when building the population histograms 
	
	int N_TRUE_ALLOWED_POWERS;	// How many UNIQUE sets of powers do we have?
	int N_ALLOWED_POWERS;		// How many total sets of powers do we have?
	
	vector<XYZ_INT> ALLOWED_POWERS;	// This will keep a list of the allowed polynomial powers for each coefficient
	vector<int>		EQUIV_INDICIES;	// For each set of allowed powers, what is the index of the first equivalent set? For example, for the set (OO, OH, OH), (1,0,1) and (1,1,0) are is equivalent
	vector<int>		PARAM_INDICIES;	// For each of the set of allowed powers, what would be the index in the FF? for example, for a set of EQUIV_INDICIES {0,0,2,3}, PARAM_INDICIES would be {0, 0, 1, 2}
	
	TRIPLETS()	// Default constructor
	{
		N_CFG_CONTRIB =  0;
		MIN_FOUND.X   = -1;
		MIN_FOUND.X   = -1;
		MIN_FOUND.X   = -1;
		
	}
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
	bool INCLUDE_FCUT;
	bool INCLUDE_CHARGES;
	bool INCLUDE_PENALTY;
	bool DO_4D;	// Do we want to include the 4D 3-body scan data??
	
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
		
	
	PES_PLOTS():N_PLOTS(0), N_SCAN(0), INCLUDE_2B(0), DO_4D(1){}
		
};

struct CHARGE_CONSTRAINT
{
	vector<string> PAIRTYPE;		// This will be read in from the input file
	vector<int>	   PAIRTYPE_IDX;	// This will be set using the map
	vector<double> CONSTRAINTS;
	double		   FORCE;
};

struct ANSI_COLORS
{
	string MAGENTA;  
	string BLUE;	
	string GREEN;	
	string PINK;	
	string RED;	
	string BOLD;	
	string UNDERLINE;
	string ENDSTYLE; 
};

static const ANSI_COLORS COUT_STYLE  = 
{
	"\033[95m",    // MAGENTA  
	"\033[94m",    // BLUE     
	"\033[92m",    // GREEN    
	"\033[93m",    // PINK     
	"\033[91m",    // RED	   
	"\033[1m ",    // BOLD     
	"\033[4m ",    // UNDERLINE
	"\033[0m ",    // ENDSTYLE 

};

/*

class CELLS
{
	public:
		
		int MAX_CUTOFF;		// What is the maximum force field cutoff?
		int TOTAL_CELLS;	// What is the total number of cells we have?
		bool USE_CELLS;		// Don't use cells if our box is very small
		
		void INITIALIZE(FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
		
		CELLS();
		~CELLS();
	
	private:
		
		XYZ WIDTHS;			// What are the x, y, and z widths of each cell?
		XYZ_INT N_CELLS;	// What is the number of cells in the x, y, and z direction?
		
		vector<XYZ>			LOWER_BOUNDS;	// What is the lower-bound allowed distance for each cell?
		vector<XYZ>			UPPER_BOUNDS;
		vector<int>			MY_CELL;		// What cell does each atom belong to?
		vector<vector<int>> NEARBY;			// What are the immediate neighbors of each cell?	
};

void CELLS::INITIALIZE(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, int MAX_FF_CUTOFF)
{
	// Set the max cutoff
	
	MAX_CUTOFF = MAX_FF_CUTOFF;
	
	// Set the possible number of cells in x, y, and z, and the total resulting
	// number of cells
	
	N_CELLS.X = floor(SYSTEM.BOXDIM.X/MAX_CUTOFF);
	N_CELLS.Y = floor(SYSTEM.BOXDIM.Y/MAX_CUTOFF);
	N_CELLS.Z = floor(SYSTEM.BOXDIM.Z/MAX_CUTOFF);
	
	TOTAL_CELLS = N_CELLS.X * N_CELLS.Y * N_CELLS.Z;
	
	// Set the width of each box in x, y, and z

	WIDTHS.X = SYSTEM.BOXDIM.X/N_CELLS.X;
	WIDTHS.Y = SYSTEM.BOXDIM.Y/N_CELLS.Y;
	WIDTHS.Z = SYSTEM.BOXDIM.Z/N_CELLS.Z;
	
	// Set the upper and lower bounds for each box. We'll
	// use >= lower, < upper
	
	XYZ TMP_VALS;
	
	TMP_VALS.X = 0;
	TMP_VALS.Y = 0;
	TMP_VALS.Z = 0;
	
	LOWER_BOUNDS.resize(TOTAL_CELLS);
	UPPER_BOUNDS.resize(TOTAL_CELLS);
	
	// Temporary variables to help set nearest neighbor cells
	
	vector<XYZ> MIDPOINT(TOTAL_CELLS);	
	NEARBY.resize(TOTAL_CELLS);
	
	int CURRENT_CELL = 0;
	
	for(int i=0; i<N_CELLS.X; i++)
	{
		for(int j=0; j<N_CELLS.Y; j++)
		{
			for(int k=0; k<N_CELLS.Z; k++)
			{
				LOWER_BOUNDS[CURRENT_CELL].X = TMP_VALS.X;
				LOWER_BOUNDS[CURRENT_CELL].Y = TMP_VALS.Y;
				LOWER_BOUNDS[CURRENT_CELL].Z = TMP_VALS.Z;
				
				TMP_VALS.X += (i+1)*WIDTHS.X;
				TMP_VALS.Y += (j+1)*WIDTHS.Y;
				TMP_VALS.Z += (k+1)*WIDTHS.Z;
				
				UPPER_BOUNDS[CURRENT_CELL].X = TMP_VALS.X;
				UPPER_BOUNDS[CURRENT_CELL].Y = TMP_VALS.Y;
				UPPER_BOUNDS[CURRENT_CELL].Z = TMP_VALS.Z;
				
				MIDPOINT[CURRENT_CELL].X = 0.5*(UPPER_BOUNDS.X + LOWER_BOUNDS.X);
				MIDPOINT[CURRENT_CELL].Y = 0.5*(UPPER_BOUNDS.Y + LOWER_BOUNDS.Y);
				MIDPOINT[CURRENT_CELL].Z = 0.5*(UPPER_BOUNDS.Z + LOWER_BOUNDS.Z);

				CURRENT_CELL++;			
			}
		}
	}
	
	// Determine the neighboring cells for each cell
	// Remember, PBC applies here! ... Also, include
	// "self" as a neighbor, because we want to compute interactions 
	// for atoms within the cell
	
	XYZ TMP_D;
	double len_diag = sqrt(WIDTHS.X*WIDTHS.X + WIDTHS.Y*WIDTHS.Y + WIDTHS.Z*WIDTHS.Z);
	
	for(int i=0; i<TOTAL_CELLS; i++)
	{
		for (int j=i; j<TOTAL_CELLS; j++)
		{
			TMP_D.X = MIDPOINT[j].X - MIDPOINT[i].X;
			TMP_D.Y = MIDPOINT[j].Y - MIDPOINT[i].Y;
			TMP_D.Z = MIDPOINT[j].Z - MIDPOINT[i].Z;
			
			TMP_D.X -= floor( 0.5 + TMP_D.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
			TMP_D.Y -= floor( 0.5 + TMP_D.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
			TMP_D.Z -= floor( 0.5 + TMP_D.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;
			
			if( (sqrt(TMP_D.X*TMP_D.X + TMP_D.Y*TMP_D.Y + TMP_D.Z*TMP_D.Z)) < len_diag ) // Then its a neighbor.. add both lists
			{
				NEARBY[i].push_back(j);
				NEARBY[j].push_back(i);
			}
		}
	}
	
	
	
	
	
}

*/


class NEIGHBORS
{
	private:

		bool   FIRST_CALL;					// Is this the first call? if so, need to build initial list
		bool   SECOND_CALL;					// Is this the second call? If so, pick the padding distance.
		double RCUT_PADDING;				// Neighborlist cutoff is r_max + rcut_padding
		double DISPLACEMENT;
		double SAFETY;                 		// Safety factor in calculating neighbors.
		
	public:

		bool   USE;							// Do we even want to use a neighbor list?
		double CURR_VEL;
		double MAX_VEL;
		double MAX_CUTOFF;					// The maximum of all force field outer cutoffs (r_max and s_max)
		double MAX_CUTOFF_3B;
		double EWALD_CUTOFF;           		// The cutoff for Ewald interactions.
		double UPDATE_FREQ;            		// Target update frequency.
		
		vector<vector<int> > LIST;			// The actual (2B) neighbor list. Of size [atoms][neighbors]
		vector<vector<int> > LIST_EWALD;	// The Ewald neighbor list. Of size [atoms][neighbors]
		vector<vector<int> > LIST_UNORDERED;// All neighbors of particle i with i not equal to j.
		vector<vector<int> > LIST_3B;		// The 3B neighbor list (3B interactions likely have a shorter cutoff)

		void UPDATE_LIST(FRAME & SYSTEM, JOB_CONTROL & CONTROLS);	// Will check if lists need updating, and will call DO_UPDATE do so if need be
		void UPDATE_LIST(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, bool FORCE);
		void INITIALIZE(FRAME & SYSTEM);
		void INITIALIZE(FRAME & SYSTEM, double & PAD);
		void DO_UPDATE (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);	// Builds and/or updates neighbor list
		NEIGHBORS();
		~NEIGHBORS();
};

class CONSTRAINT
{
	// NOTE: If velocity scaling is requested, it is taKIN_ENERn care of OUTSIDE of this currently
	//... Eventually all temperature and pressure trackng should be done in here
	//... possible also energetic tracking
	
	private:
		
		double THERM_POSIT_T;		// thermostat position (chi)
		double THERM_VELOC_T;		// Thermostat velocity (chi-dot)... was "visco"
		double THERM_INERT_T;		// Thermostat inertia term ("G")
		double THERM_INERT_Q;		// Thermostat inertia ("Q")
		double THERM_POSIT_0;
		double THERM_INERT_0;
		double THERM_VELOC_0;
			
		double BAROS_POSIT_T;		// ("epsilon")
		double BAROS_VELOC_T;		// ("epsilon-dot")
		double BAROS_FORCE_T;		// ("F-epsilon")
		double BAROS_INERT_W;
		double BAROS_POSIT_0;		// ("epsilon")
		double BAROS_FORCE_0;
		double BAROS_VELOC_0;	
		

		//Barostat variables

		double BAROS_SCALE;
		double VOLUME_0;
		double VOLUME_T;
		
		// Berendsen barostat variables
		
		double BEREND_MU;		// Berendsen scaling factor
		XYZ    BEREND_ANI_MU;	// Berendsen anisotropic scaling factor
		double BEREND_KP;		// Berendsen compressibility divided by time constant... 
		
		// Berendsen thermostat variables ... borrows berend kp
		
		double BEREND_ETA;	// scaling factor
		double BEREND_TAU;	// time constant
		
		// Thermostat variables

		double TIME;		// Hoover time in fs
		double TIME_BARO;	// Barostat time in fs
		double N_DOF;		// # degrees of freedom
		double VSCALEH;		// Replaces 
		double KIN_ENER;	// Kinetic energy
			
	public:
		
		string STYLE;		// NPT, NVT-MTK, NVT-SCALE, NVE

		void INITIALIZE          (string IN_STYLE, JOB_CONTROL & CONTROLS, int ATOMS); 
		void UPDATE_COORDS       (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
		void UPDATE_VELOCS_HALF_1(FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
		void UPDATE_VELOCS_HALF_2(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);
		void SCALE_VELOCITIES    (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
		void UPDATE_TEMPERATURE  (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
		double CONSERVED_QUANT   (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
		
		CONSTRAINT();
		~CONSTRAINT();
		
};








//////////////////////////////////////////
//
//	FUNCTION HEADERS -- MPI
//
//////////////////////////////////////////

void divide_atoms(int &a1start, int &a1end, int atoms);


//////////////////////////////////////////
//
//	FUNCTION HEADERS -- DERIVATIVE CALCULATION
//
//////////////////////////////////////////

void ZCalc_Deriv (JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS,  vector<TRIPLETS> & PAIR_TRIPLETS, FRAME & FRAME_TRAJECTORY, vector<vector <XYZ > > & FRAME_A_MATRIX, vector<vector <XYZ > > & FRAME_COULOMB_FORCES, const int nlayers, bool if_3b_cheby, map<string,int> & PAIR_MAP,  map<string,int> & TRIAD_MAP, NEIGHBORS & NEIGHBOR_LIST);

void SubtractCoordForces (FRAME & TRAJECTORY, bool calc_deriv, vector<XYZ> & P_OVER_FORCES,  vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

void SubtractEwaldForces(FRAME &SYSTEM, NEIGHBORS &NEIGHBOR_LIST, JOB_CONTROL &CONTROLS);

void ZCalc_Ewald(FRAME & TRAJECTORY, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);

void optimal_ewald_params(double accuracy, int nat, double &alpha, double & rc, int & kc, double & r_acc, double & k_acc, XYZ boxdim);

void ZCalc_Ewald_Deriv(FRAME & FRAME_TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, vector <vector <XYZ > > & FRAME_COULOMB_FORCES, map<string,int> & PAIR_MAP,NEIGHBORS & NEIGHBOR_LIST, JOB_CONTROL & CONTROLS);

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- FORCE CALCULATION
//
//////////////////////////////////////////

void   ZCalc                    (FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, NEIGHBORS & NEIGHBOR_LIST);
void   ZCalc_3B_Cheby_Deriv_HIST(JOB_CONTROL & CONTROLS, vector<PAIRS> & FF_2BODY, vector<TRIPLETS> & PAIR_TRIPLETS, vector<vector <vector< XYZ > > > & A_MATRIX, map<string,int> PAIR_MAP, map<string,int> TRIAD_MAP);		
double get_dist                 (const FRAME & SYSTEM, XYZ & RAB, int a1, int a2);

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- PRINTING OF POTENTIAL ENERGY SURFACE
//
//////////////////////////////////////////

void Print_Cheby(vector<PAIR_FF> & FF_2BODY, int ij, string PAIR_NAME, bool INCLUDE_FCUT, bool INCLUDE_CHARGES, bool INCLUDE_PENALTY, string FILE_TAG);
void Print_3B_Cheby(JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk);
void Print_3B_Cheby_Scan(JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan);
void Print_Ternary_Cheby_Scan(JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan);


//////////////////////////////////////////
//
//	FUNCTION HEADERS -- ASSORTED
//
//////////////////////////////////////////

double kinetic_energy(FRAME & SYSTEM, JOB_CONTROL & CONTROLS);					// Overloaded.. compute differentely if for main or new velocities
double kinetic_energy(FRAME & SYSTEM, string TYPE, JOB_CONTROL & CONTROLS);		// Overloaded.. compute differentely if for main or new velocities

void sync_layers(FRAME &SYSTEM, JOB_CONTROL &CONTROLS) ;
void build_layers(FRAME &SYSTEM, JOB_CONTROL &CONTROLS) ;

void enable_fp_exceptions();
void exit_run(int val);



#endif

