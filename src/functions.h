#ifndef VERBOSITY
	#define VERBOSITY 1
#endif

#ifndef FORCECHECK
	#define FORCECHECK 0
#endif

#ifndef WARN
	#define WARN TRUE
#endif 

#ifndef _HELPERS_
#define _HELPERS_

// For "original" MD code

#include<stdio.h>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<algorithm>
#include<assert.h>
#include<map>

using namespace std;

// Maximum number of atom types and powers.

#define MAX_ATOM_TYPES   10 
#define MAX_ATOM_TYPES2 (MAX_ATOM_TYPES  * MAX_ATOM_TYPES)
#define MAX_ATOM_TYPES3 (MAX_ATOM_TYPES2 * MAX_ATOM_TYPES) 
#define MAX_ATOM_TYPES4 (MAX_ATOM_TYPES3 * MAX_ATOM_TYPES) 

// Unit converters 
static const double ke      = 332.0637157615209;// Converter between electron units and Stillinger units for Charge*Charge.
static const double Hartree = 627.50961;	// 1 Hartree in kcal/mol.
static const double Bohr    = 1.889725989;     // 1 Angstrom in Bohr.
static const double Kb      = 0.0019872067 ;		// Boltzmann constant in kcal/mol-K. (Copied from LAMMPS).
static const double Tfs     = 48.888;		// Internal time unit in fs.
static const double GPa     = 6.9476955;		// Unit conversion factor... kcal/mol/A^3 * (this constant) ==> GPa
static const double atm     = GPa*9869.23266716;// stress conversion to atm (for LAMMPS).
static const double GPa2atm = 9869.23266716;	// x_GPa * GPa2atm = x_atm
static const double pi      = 3.14159265359;
static const double mpers   = 0.02045482828; 	// conversion from internal velocity units (kcal^1/2/g^1/2) to m/s... v_internal_units * this const --> m/s ... note: 1 kcal/g = 4184000 m^2/s^2		

// Global variables declared as externs in functions.h, and declared in functions.C -- general

extern string FULL_FILE_3B;	// The 4D PES for 3B FF
extern string SCAN_FILE_3B;	// The 2D PES scans for 3B
extern string SCAN_FILE_2B;	// The 2D PES scans for 2B

// Global variables declared as externs in functions.h, and declared in functions.C -- MPI calculations.   
 
extern int NPROCS;		// Number of processors
extern int RANK;		// Index of current processor

// Enumerated classes.

enum class Cheby_trans	// Supported variable transformations.
{
	MORSE, 
	INVRSE_R, 
	NONE
};

// Algorithms to fix chebyshev polynomials when r < s_minim.
enum class Cheby_fix {
	ZERO_DERIV,
	CONSTANT_DERIV,
	SMOOTH
} ;

// Include Chimes files here.

#include "Fcut.h"
#include "Cluster.h"
#include "A_Matrix.h"

class JOB_CONTROL
{
public:
	int STEP;	// Tracks the current md step
	int NATMTYP;	// How many atom types are in the trajectory?
	
	///////////////////////////////////////////////
	// Variables read in from the parameter file:
	///////////////////////////////////////////////
	
	bool   FIT_COUL;	      // Replaces fit_coul... If true, take charges from spline parameters.
	bool   USE_COULOMB;	      // Replaces if_coulomb... If true, calculate Coulomb forces.
	bool   USE_3B_CHEBY;	      // Replaces if_3b_cheby... If true, calculate 3-Body Chebyshev interaction.
	bool   USE_4B_CHEBY;	      // If true, calculate 4-Body Chebyshev interaction.
	vector<string> ATOMTYPES;     // A list of atom types	      
        
	
	///////////////////////////////////////////////
	// Variables read in from main md code input file
	///////////////////////////////////////////////

	// "General control variables"

	int    SEED;		      // Replaces rand_seed...
	double TEMPERATURE;	      // Replaces TempMD
	double PRESSURE;	      // Introduced for NPT
	string ENSEMBLE;	      // NVE, NPT, NVT
	bool   CHECK_FORCE;	      // If true, numerically check forces from derivatives of energy.
	bool   COMPARE_FORCE;	      // Replaces if_read_force... If TRUE, read in read in a set of forces from a file for comparison with those computed by this code

	bool   SUBTRACT_FORCE;        // Read frame, compute forces based on parameter file, print out frame where FF forces have been subtracted from input forces
	string COMPARE_FILE;	      // Name of the file that holds the forces
	double DELTA_T; 	      // Relaces deltat
	double DELTA_T_FS;	      // Replaces deltat_fs... in femtoseconds
	int    N_MD_STEPS;	      // Replaces nsteps...
	int    N_LAYERS;	      // Replaces nlayers... Number of periodic images to use. Usually, no need to go past 1 ...(depends only on spline interaction cutoff).
	int    REAL_REPLICATES;       // Number of real replicate layers to create. Only goes in the positive direction, so 1 replicate corresponds to 8*attoms
	double SCALE_SYSTEM_BY;       // Amount to scale boxlengths/coordinates by
	string PARAM_FILE;	      // Replaces params_file...
	bool   SERIAL_CHIMES ;    // Use serial chimes calculator ?
	vector<string> COORD_FILE;    // Replaces xyz_file... Can be a list of files, which will be assembled together along the z-axis to form a single cell
	bool   SELF_CONSIST;	      // Is this part of a self-consistent DFT MD --> FIT --> MM MD --> CYCLE type calculation?
	double NVT_CONV_CUT;	      // What is the cutoff for "conservation"... Will kill the program if the current temperature and set temperature differ by more than this fraction
	// ie if | T_set - T_curr | >  NVT_CONV_CUT * T_set, end with error message.. defaults to 0.1 (10%)

	// "Simulation options"

	bool   INIT_VEL;	      // Replaces if_init_vel... If true, initialize velocities.
	bool   USE_HOOVER_THRMOSTAT;  // Replaces if_hoover... Use a Nose-Hoover thermostat?  
	int    FREEZE_IDX_START;      // First atom in continuous range to freeze.. .counting starts from 0... value of -1 indicates no atoms are to be frozen
	int    FREEZE_IDX_STOP;       // Last atom in continuous range to freeze 
	double FREQ_UPDATE_THERMOSTAT;// Replaces scale_freq and thoover_fs... it's usage depends on whether USE_HOOVER_THERMOSTAT is true or false.. will be cast as int where required
	double FREQ_UPDATE_BAROSTAT;  // Barostat time constant... defaults to 1000
	bool   USE_NUMERICAL_PRESS;   // Replaces num_pressure... Whether to calculate pressures by finite difference.
	bool   USE_NUMERICAL_STRESS;   // Whether to calculate the stress tensor by finite difference.	

	// For penalty-function related exit

	double PENALTY_THRESH;        // [0,1]; If a penalty potential more than PENALTY_THRESH*<current per-atom e_cons value> is produced, kill the simulation... it's a runaway
	double IO_ECONS_VAL;	      // The conserved quantity output to screen
	bool USE_KILL_LEN;	      // If true, and if a cheby type interaction, will kill the simulation if pair distances below FF_2BODY[i].KILLLEN are found

	// "Output control" 

	bool   INCLUDE_ATOM_OFFSETS;  // If true, single-atom contirubtions are added to all output energies. Default is false.
	int    FREQ_DFTB_GEN;	      // Replaces gen_freq... How often to write the gen file.
	string TRAJ_FORMAT;	      // .gen, .xyzf, or .lammps (currently)
	bool   SPLIT_FILES ;	      // If TRUE, do not concatenate A matrix files for LSQ.
	bool   HIERARCHICAL_FIT;      // If true, allows 2-body and 1-body interctions to be excluded from fitting
	int    FREQ_BACKUP;	      // How often to write backup files for restart.
	bool   PRINT_VELOC;	      // If true, write out the velocities 
	bool   RESTART; 	      // If true, read a restart file.
	int    FREQ_VELOC;
	int    FREQ_ENER;	      // Replaces energy_freq... How often to output energy
	bool   PRINT_FORCE;	      // Replaces if_output_force... If TRUE, write out calculated forces.
	bool   PRINT_ENERGY_STRESS ;  // If TRUE, add a header with potential energy and configurational stress
	                              // to force output file.
	bool   PRINT_BAD_CFGS;        // Print any config where r < r_cut
	int    FREQ_FORCE;	      // How often to print the forces        
	int    SELF_CONSIST_FREQ;     // How frequently to print POSCAR file
	bool   WRAP_COORDS;	      // Should coordinates be wrapped?
	bool   FORDFTB;	              // Write a special output file for DFTB+ to read in? (default = false)
	int    SKIP_FRAMES ;       // Should parallel processing of LSQ frames skip (>= 1) or be contiguous (0) ?
	
	// Controls for how to construct the initial system, if desired

	bool   BUILD;
	string BUILD_TYPE;
	string BUILD_FILE;
	double BUILD_BOXL;
	int    BUILD_NMOLEC;
	string BUILD_ATOM;

	Cheby_fix cheby_fix_type ;
	double cheby_smooth_distance ;
	
	///////////////////////////////////////////////
	// Variables exclusive to the LSQ code
	///////////////////////////////////////////////

	bool IS_LSQ;		      // Is this for an lsq run or actual md?
	bool FIT_STRESS;	      // Should stress tensors be included in the fit? --> This is ONLY for the diagonal components, xx, yy, zz
	bool FIT_STRESS_ALL;	      // Should stress tensors be included in the fit? --> This is ONLY for ALL components, xx, xy, xz ... zz 
	int  NSTRESS;		      // Only fit stresses for first NSTRESS frames of trajectory
	bool FIT_ENER;  	      // Should the total frame energy be included in the fit?
	bool FIT_ENER_EVER ;	      // Is energy ever included in the fit ?
	int  NENER;
	bool CALL_EWALD;	      // Should ewald subroutines be called?

	int   NFRAMES;  	      // Number of frames in the movie file
	int   CHEBY_ORDER;	      // Order of Chebyshev polynomial if used... set to 8 for DFTB Erep polynomial
	int   CHEBY_3B_ORDER;	      // how many polynomials for 3b cheby?
	int   CHEBY_4B_ORDER;	      // how many polynomials for 4b cheby?
	int   NUM_3B_CHEBY;	      // How many parameters are associated with cheby order CHEBY_3B_ORDER?
	int   NUM_4B_CHEBY;	      // How many parameters are associated with cheby order CHEBY_4B_ORDER?
	int   INVR_PARAMS;	      // currently uses 19 parameters per pair type
	int   TOT_SNUM; 		      // total number of 2-body force field parameters
	int   TOT_SHORT_RANGE;        // Number of short ranged FF params... i.e. not Ewald
	int   TOT_ALL_PARAMS ;        // Total number of LSQ fitting parameters.

	bool  COUL_CONSV;	      // If true, constraints will be applied to charge fitting to try to maintain consistency
	bool  IF_SUBTRACT_COORD;      // If true, subtract overcoordination forces.
	bool  IF_SUBTRACT_COUL;       // If true, subtract Coulombic forces (for use with fixed charges).
	bool  USE_PARTIAL_CHARGES;    // Will there be any charges in the system?

	Cheby_trans CHEBY_TYPE;       // How will distance be transformed?
	vector<string> INFILE;        // Input trajectory file
	vector<int> INFILE_FRAMES;    // How many frames should we read from each file?
	vector<string> INFILE_FORCE_FLAGS;	// Specify a string to prepend to force  labels in b-labels.txt file... default is an empty string
	vector<string> INFILE_STRESS_FLAGS;	// Specify a string to prepend to stress labels in b-labels.txt file... default is an empty string
	vector<string> INFILE_ENERGY_FLAGS;	// Specify a string to prepend to energy labels in b-labels.txt file... default is an empty string

	// These variables are temporary fixes - soon we will have a class to do all input reading

	string FCUT_LINE;
  

	// Constructor... MD values are set in the read_input function in chimes_md.C
	
	JOB_CONTROL(): FIT_COUL(false), 
				   USE_3B_CHEBY(false), 
				   USE_4B_CHEBY(false), 
				   N_LAYERS(0), 
				   WRAP_COORDS(false),
				   TOT_SNUM(0), 
				   COUL_CONSV(false), 
				   IF_SUBTRACT_COORD(false),
				   IF_SUBTRACT_COUL(false),
				   USE_PARTIAL_CHARGES(false),
				   PRINT_ENERGY_STRESS(false)
	{
		NFRAMES         = 0;	// Number of frames in the movie file
		CHEBY_ORDER     = 0;	// Order of Chebyshev polynomial if used... set to 8 for DFTB Erep polynomial
		CHEBY_TYPE      = Cheby_trans::NONE ;
		CHEBY_3B_ORDER  = 0;   
		CHEBY_4B_ORDER  = 0;	// how many polynomials for 4b cheby?
		NUM_3B_CHEBY    = 0;	// How many parameters are associated with cheby order CHEBY_3B_ORDER?
		NUM_4B_CHEBY    = 0;	// How many parameters are associated with cheby order CHEBY_4B_ORDER?
		TOT_SNUM        = 0;	// total number of force field parameters
		TOT_SHORT_RANGE = 0;	// Number of short tranged FF params... i.e. not Ewald

		CHECK_FORCE  = false;
		USE_3B_CHEBY = false;	// Replaces if_3b_cheby... If true, calculate 3-Body Chebyshev interaction.
		USE_4B_CHEBY = false;	//If true, calculate 4-Body Chebyshev interaction.
		SPLIT_FILES  = false ;
		TOT_ALL_PARAMS = 0 ;
		SERIAL_CHIMES = false ;
		USE_KILL_LEN = false;
		//IO_ECONS_VAL = 0.0;
		SKIP_FRAMES = 0 ;
		
		FCUT_LINE = "CUBIC";
		FIT_ENER_EVER = false ;
			

		// Default is the same as the ChIMES calculator.
		cheby_fix_type = Cheby_fix::SMOOTH ;
		cheby_smooth_distance = 0.01 ;
		
		COMPARE_FORCE     = false;	// is this variable really necessary for LSQ?
		CALL_EWALD        = false;
		FIT_ENER          = false;
		FIT_STRESS        = false;
		FIT_STRESS_ALL    = false;
		NSTRESS           = -1;
		NENER             = -1;
		FORDFTB           = false;
		
	}
	void LSQ_SETUP(int npairs, int no_atom_types) ; // Set up JOB_CONTROL for LSQ calculation.
};

class BOX
{
	 /*

		 Notes for usage with LAMMPS linking: 	Box origin must be at 0,0,0, and all cell vector components must be positive 

	 */
	
public:
		
	 // General cell properties
		
	 
	 double VOL;				// Cell volume
		
	 // Cell geometery... if IS_ORTHO is true, ONLY CELL_AX, CELL_BY, and CELL_CZ are ever modified

	 double CELL_AX, CELL_AY, CELL_AZ;	// Cell vectors (transpose gives h-mat)
	 double CELL_BX, CELL_BY, CELL_BZ;
	 double CELL_CX, CELL_CY, CELL_CZ;

	 vector<double> HMAT;
	 vector<double> INVR_HMAT;
		
	 double LATCON_A,  LATCON_B, LATCON_C;	// Cell lattice constants
	 double LAT_ALPHA, LAT_BETA, LAT_GAMMA;	// Cell lattice angles
		
	 double CELL_LX, CELL_LY, CELL_LZ;	// Cell vector lengths (for LAMMPS)
		
	 double XY, XZ, YZ;			// Cell tilt angles (for LAMMPS)

	double xlo, ylo, zlo ;      // LAMMPS-style cell extents.
	double xhi, yhi, zhi ;
		
	 double EXTENT_X, EXTENT_Y, EXTENT_Z;	// Cell extent in the x, y, and z directions

	 bool   IS_ORTHO;			// Is this an orthorhombic cell? If so, we'll use the more computationally efficient operations

	 bool IS_VARIABLE ;     // Does this box vary in time ?

	 BOX();
	 BOX(const BOX & COPY_FROM);
	 ~BOX();
		
	 void READ_BOX();					// Reads box geometery (orth or non-orth) from .xyz(f) file
	 void WRITE_BOX(int LAYERS);
				
	 void   UPDATE_INVER_CELL();				// Computes the inverse h-mat from the h-mat
	 void   UPDATE_LAT_VALUES(); 
	 void   UPDATE_EXTENT();
	 double UPDATE_VOLUME();					// Computes cell volume and saves to VOL
	 void   UPDATE_CELL();					// Updates all cell geometry variables, assuming the h-matrix has been externally modified
		
	 //void PREPARE_COORDS();					// Properly rotates/aligns/shifts coordinate origin to 0,0,0

   // Converts LAMMPS  (lx,ly,lz) = (xhi-xlo,yhi-ylo,zhi-zlo) and tilt factors (xy,xz,yz) to lattice constants and angles alpha, beta, and gamma	 
	 void UNLAMMPSIFY(double lx, double ly, double lz, double xy, double xz, double yz); 

	 void LAMMPSIFY(); // ??		

   // Wraps an atom in box ... Saves wrapped coordinates to UN_WRAPPED_ATOM	 
	 void WRAP_ATOM(XYZ & UNWRAPPED_ATOM, XYZ_INT & WRAP_IDX, bool UPDATE_WRAPDIM);												

   // Wraps an atom in box ... Saves wrapped coordinates to WRAPPED_ATOM or UN_WRAPPED_ATOM
	 void WRAP_ATOM(XYZ & UNWRAPPED_ATOM, XYZ & WRAPPED_ATOM, XYZ_INT & WRAP_IDX, bool UPDATE_WRAPDIM);			

   // Creates replicate atoms (ghost or real) ... Saves replicate coordinates in REFERENCE_ATOM	 
	 void LAYER_ATOM(XYZ & REFERENCE_ATOM, XYZ_INT & LAYER_INDEX);				

	// Creates replicate atoms (ghost or real) ... Saves replicate coordinates in REFERENCE_ATOM or LAYERED_ATOM	 
	 void LAYER_ATOM(XYZ & REFERENCE_ATOM, XYZ_INT & LAYER_INDEX, XYZ & LAYERED_ATOM);

   // Returns "true" if rcut is less than half the shortest cell vector (general for triclinic systems)	 
	 bool IS_RCUT_SAFE(double CUTOFF, int LAYERS);					
		
	 void SCALE_BY_FACTOR(double FACTOR);

	void SCALE_BY_MATRIX(const vector<vector<double>> &FACTOR, bool SCALE_ATOMS, XYZ & ATOM)	;

   // Multiplies cell vectors times a scalar and shifts atoms accordingly, if requested
	 void SCALE_BY_FACTOR(double FACTOR, bool SCALE_ATOMS, XYZ & ATOM);		

   // Updates RAB to contain distance vectors	 
	 void GET_DISTANCE(const XYZ & ATOM1, const XYZ & ATOM2, XYZ & RAB, bool USE_MIC);

};

class FRAME
{
public:
	 int ATOMS;             		// Just the parent atoms.
	 int ALL_ATOMS;         	   	// All atoms, including ghosts. 

		int MY_ATOMS;			// Used for lammps linking. Specify how many atoms in SYS the process owns
		int MY_ATOMS_START;		// Used for lammps linking. Specify what index along SYS starts the process' atoms
	 BOX BOXDIM;			// Dimenions of the primitive box.
	 XYZ STRESS_TENSORS;		// Only used for the diagonal components, xx, yy, zz
	 XYZ STRESS_TENSORS_X;		// Used when all tensor components are requested ... used primarily for the lsq code. When used
	 XYZ STRESS_TENSORS_Y;		// by the MD code, stores stress tensors for comparative purposes (e.g. from DFT)
	 XYZ STRESS_TENSORS_Z;

	 double         QM_POT_ENER;		// This is the potential energy of the QM calculation!
	 vector<double> QM_ENERGY_OFFSET;    	// This is the energy offset between MD and QM energy, as determined by lsq[2].py.
	 vector<double> QM_POT_ENER_PER_ATOM;	// And this is for each atom in the frame, from QM
	 vector<int>    NATOMS_OF_TYPE;	    	// How many atoms of each type there are

	 double 	TEMPERATURE;				// This is the RUNNING temperature, not the set temperature!
	 double 	PRESSURE;					// This is the RUNNING pressure, not the set pressure!
	 double 	AVG_TEMPERATURE;			// Only used for velocity scaling-type thermostats
	 double 	PRESSURE_XYZ;				// This is the running pressure sans the ideal gas term
	 //XYZ	PRESSURE_TENSORS_XYZ;			// These are the RUNNING pressure tensors, no the set pressure tensors! ...sans the ideal gas term
	 vector<XYZ> PRESSURE_TENSORS_XYZ_ALL; 		// These are the RUNNING pressure tensors, no the set pressure tensors! ...sans the ideal gas term ... includes off-diagonals
	 //XYZ	PRESSURE_TENSORS;			// Adds in the ideal gas term
	 vector<XYZ> PRESSURE_TENSORS_ALL;		// Adds in the ideal gas term ... includes off-diagonals
	 double	TOT_POT_ENER;				// Replaces VTOT

	 vector<int> PARENT;
	 vector<XYZ_INT> LAYER_IDX;
	 vector<XYZ_INT> WRAP_IDX;  	// Index of box wrapping for this atom.
	 vector<string> ATOMTYPE;
	 vector<int> 	ATOMTYPE_IDX;	// Only used for dftbgen and LAMMPStrj file printing
	 vector<XYZ> COORDS;         // Current step coordinates.
	 vector<XYZ> COORDS0;        // Prior step coordinates.	 
	 vector<XYZ>     ALL_COORDS;  	// Coordinates of atoms + ghosts used for force evaluation.
	 vector<double> 	CHARGES;
	 vector<double> 	MASS;
	 vector<XYZ>	FORCES;
	 vector<XYZ>     ACCEL;
	 vector<XYZ>     TMP_EWALD;	// Holds temporary ewald accels/forces
	 vector<XYZ>     VELOCITY;
	 vector<XYZ>	VELOCITY_NEW;
	 vector<XYZ>     VELOCITY_ITER;

	 // Update ghost atom positions.

	 void 		update_ghost(int n_layers, bool UPDATE_WRAPDIM);
	 inline int 	get_atomtype_idx(int atom);

	 void SET_NATOMS_OF_TYPE();
	 void READ_XYZF(ifstream &TRAJ_INPUT, const JOB_CONTROL &CONTROLS, const vector<PAIRS> &ATOM_PAIRS, const vector<string> &TMP_ATOMTYPE, int i);
	 void build_layers(int N_LAYERS) ;
};

struct CHARGE_CONSTRAINT
{
	vector<string> PAIRTYPE;		// This will be read in from the input file
	vector<int>	   PAIRTYPE_IDX;	// This will be set using the map
	vector<double> CONSTRAINTS;
	double		   FORCE;
};
	
class INTERACTION_3B
// A 3-body interaction.
{
	// Indices ordered such that i < parent[j] < parent[k].
	// a1 is a non-ghost atom.
	
	public:
		int a1;  // Atom 1.
		int a2;  // Atom 2.
		int a3;  // Atom 3.
};

class INTERACTION_4B
// A 4-body interaction.
{
	// Indices ordered such that i < parent[j] < parent[k] < parent[l].
	// a1 is a non-ghost atom.
	
	public:
		int a1;  // Atom 1.
		int a2;  // Atom 2.
		int a3;  // Atom 3.
		int a4;  // Atom 4.
};
  
class NEIGHBORS
{
private:

	 bool   FIRST_CALL;						// Is this the first call? if so, need to build initial list
	 bool   SECOND_CALL;						// Is this the second call? If so, pick the padding distance.
	 double DISPLACEMENT;
	 double SAFETY;                 					// Safety factor in calculating neighbors.
		
	 void FIX_LAYERS(FRAME & SYSTEM, JOB_CONTROL & CONTROLS);	// Updates ghost atoms based on pbc-wrapped real atoms
	 void DO_UPDATE_SMALL (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);	// Builds and/or updates neighbor list
	 void DO_UPDATE_BIG   (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);	// Builds and/or updates neighbor list
	 void UPDATE_3B_INTERACTION(FRAME & SYSTEM, JOB_CONTROL &CONTROLS);  // Update 3-Body interaction list.
	 void UPDATE_4B_INTERACTION(FRAME & SYSTEM, JOB_CONTROL &CONTROLS);  // Update 4-Body interaction list.

public:

	 bool   UPDATE_WITH_BIG;			// Should we update our neighbor list with DO_UPDATE_BIG? If false, uses DO_UPDATE_SMALL
	 double RCUT_PADDING;			// Neighborlist cutoff is r_max + rcut_padding
	 bool   USE;				// Do we even want to use a neighbor list?
	 double CURR_VEL;
	 double MAX_VEL;
	 double MAX_COORD_STEP;	 // The maximum step of any coordinate.
	 double MAX_CUTOFF;			// The maximum of all force field outer cutoffs (r_max and s_max)
	 double MAX_CUTOFF_3B;
	 double MAX_CUTOFF_4B;
	 double EWALD_CUTOFF;           		// The cutoff for Ewald interactions.
	 double UPDATE_FREQ;            		// Target update frequency.

	 vector<double> PERM_SCALE ;       // Scaling factor for self-interactions.
	 
	 vector<vector<int> > LIST;		// The actual (2B) neighbor list. Of size [atoms][neighbors]
	 vector<vector<int> > LIST_EWALD;	// The Ewald neighbor list. Of size [atoms][neighbors]
	 vector<vector<int> > LIST_UNORDERED;	// All neighbors of particle i with i not equal to j.
	 vector<vector<int> > LIST_3B;		// The 3B neighbor list (3B interactions likely have a shorter cutoff)
	 vector<vector<int> > LIST_4B;		// The 3B neighbor list (3B interactions likely have a shorter cutoff)

	 vector<INTERACTION_3B> LIST_3B_INT;    // A flat list of all 3-body interactions.
	 vector<INTERACTION_4B> LIST_4B_INT;    // A flat list of all 3-body interactions.

	 void UPDATE_LIST(FRAME & SYSTEM, JOB_CONTROL & CONTROLS);	// Will check if lists need updating, and will call DO_UPDATE do so if need be
	 //void UPDATE_LIST(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, bool FORCE);
	 void INITIALIZE(FRAME & SYSTEM);
	 void INITIALIZE_MD(FRAME & SYSTEM, JOB_CONTROL &CONTROLS) ;
	 void INITIALIZE(FRAME & SYSTEM, double & PAD);
	 void DO_UPDATE (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);	// Builds and/or updates neighbor list

	 double MAX_ALL_CUTOFFS() ;
	 
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
	 double THERM_FORCE_T;		// Thermostat force ("G")
	 double THERM_INERT_Q;		// Thermostat inertia ("Q")
	 double THERM_POSIT_0;
	 double THERM_FORCE_0;
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

	 void WRITE(ofstream &STREAM);  // Write variables to file.
	 void READ(ifstream &STREAM);   // Read variables from file.

	 void INITIALIZE          (string IN_STYLE, JOB_CONTROL & CONTROLS, int ATOMS); 
		
	 void INIT_VEL (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);	// Use box Muller to initialize velocities
	 void CHECK_VEL(FRAME & SYSTEM, JOB_CONTROL & CONTROLS); //Test/fix velocity center of mass
		
	 XYZ CENTER_OF_MASS(const FRAME &SYSTEM) ;
 	 void UPDATE_COORDS       (FRAME & SYSTEM, JOB_CONTROL & CONTROLS, NEIGHBORS &NEIGHBORS);
	 void UPDATE_VELOCS_HALF_1(FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
	 void UPDATE_VELOCS_HALF_2(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);
	 void SCALE_VELOCITIES    (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
	 void UPDATE_TEMPERATURE  (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
	 double CONSERVED_QUANT   (FRAME & SYSTEM, JOB_CONTROL & CONTROLS);
		
	 CONSTRAINT();
	 ~CONSTRAINT();
		
};


class THERMO_AVG 
{
	 // Time average of thermodynamic properties.
	
public:
		
	 double	TEMP_SUM;
	 double	PRESS_SUM;
	 double VOLUME_SUM ;
	 double PV_SUM ;
	 XYZ 	STRESS_TENSOR_SUM;
	 vector<XYZ> STRESS_TENSOR_SUM_ALL;
		
	 void WRITE(ofstream &fout);
	 void READ(ifstream &fin);
		
	 THERMO_AVG() 
			{
				 // Zero the initial value of all sums.
		
				 TEMP_SUM  = 0.0;
				 PRESS_SUM = 0.0;
				 VOLUME_SUM = 0.0 ;
				 PV_SUM = 0.0 ;
				 
				 STRESS_TENSOR_SUM.X = 0.0;
				 STRESS_TENSOR_SUM.Y = 0.0;
				 STRESS_TENSOR_SUM.Z = 0.0;
			}
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


void ZCalc_Deriv (JOB_CONTROL & CONTROLS, vector<PAIRS> & FF_2BODY,  CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, FRAME & FRAME_SYSTEM, A_MAT & A_MATRIX,
									map<string,int> & PAIR_MAP,  vector<int> &INT_PAIR_MAP, NEIGHBORS &NEIGHBOR_LIST) ;

void SubtractEwaldForces (FRAME &SYSTEM, NEIGHBORS &NEIGHBOR_LIST, JOB_CONTROL &CONTROLS);

void ZCalc_Ewald         (FRAME & TRAJECTORY, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);

void optimal_ewald_params(double accuracy, int nat, double &alpha, double & rc, int & kc, double & r_acc, double & k_acc, BOX boxdim);

void ZCalc_Ewald_Deriv(FRAME & FRAME_TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, A_MAT & A_MATRIX, map<string,int> & PAIR_MAP,NEIGHBORS & NEIGHBOR_LIST, JOB_CONTROL & CONTROLS) ;

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- FORCE CALCULATION
//
//////////////////////////////////////////

void   ZCalc(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> &PAIR_MAP, vector<int> &INT_PAIR_MAP, CLUSTER_LIST& TRIPS, CLUSTER_LIST &QUADS,  NEIGHBORS & NEIGHBOR_LIST);

//////////////////////////////////////////
// Distance calculation and smoothing functions
//////////////////////////////////////////

inline double get_dist(FRAME & SYSTEM, XYZ & RAB, int a1, int a2)
// Calculates distance as a2 - a1... This function modifies RAB!
{
	if(SYSTEM.ATOMS == SYSTEM.ALL_ATOMS) // Then we're not using ghost atoms. Need to use MIC
		SYSTEM.BOXDIM.GET_DISTANCE(SYSTEM.COORDS[a1], SYSTEM.COORDS[a2], RAB, true);
	else
		SYSTEM.BOXDIM.GET_DISTANCE(SYSTEM.ALL_COORDS[a1], SYSTEM.ALL_COORDS[a2], RAB, false);

	
	//cout << "	" << sqrt(RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z) << endl;
	
	return sqrt(RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z);
}


inline int FRAME::get_atomtype_idx(int atom)
{
  		return(ATOMTYPE_IDX[atom]);
}

//////////////////////////////////////////
// Trig functions for box operations
//////////////////////////////////////////

double VECTOR_MAGNITUDE(vector<double> & vec);
double VECTOR_ANGLE(vector<double> & v1, vector<double> & v2);



//////////////////////////////////////////
//
//	FUNCTION HEADERS -- ASSORTED
//
//////////////////////////////////////////

void OPEN_TRAJFILE(ifstream & TRAJ_INPUT, vector<string> & INFILE, int FILE_IDX);

double kinetic_energy(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<XYZ> &Ktensor);			// Overloaded.. compute differentely if for main or new velocities
double kinetic_energy(FRAME & SYSTEM, string TYPE, JOB_CONTROL & CONTROLS);	// Overloaded.. compute differentely if for main or new velocities

void build_layers      (FRAME &SYSTEM, JOB_CONTROL &CONTROLS);
void build_real_replicates(FRAME &SYSTEM, const JOB_CONTROL &CONTROLS);

void numerical_pressure(const FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, CLUSTER_LIST & TRIPS,  CLUSTER_LIST &QUADS, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST,double & PE_1, double & PE_2, double & dV);
void numerical_stress(const FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY,
											CLUSTER_LIST & TRIPS,  CLUSTER_LIST &QUADS, map<string,int> & PAIR_MAP,
											vector<int> &INT_PAIR_MAP,
											NEIGHBORS & NEIGHBOR_LIST,double & PE_1, double &PE_2, 
											int it1, int it2) ;
void numerical_stress_all(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY,
													CLUSTER_LIST & TRIPS,  CLUSTER_LIST &QUADS, map<string,int> & PAIR_MAP,
													vector<int> &INT_PAIR_MAP,
													NEIGHBORS & NEIGHBOR_LIST) ;
void check_forces(FRAME& SYSTEM, JOB_CONTROL &CONTROLS, vector<PAIR_FF> &FF_2BODY, map<string,int>& PAIR_MAP, vector<int> &INT_PAIR_MAP,  CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, NEIGHBORS &NEIGHBOR_LIST);
void build_int_pair_map(int natmtyp, const vector<string> &atomtype, const vector<int> &atomtype_idx, map<string,int> &pair_map, vector<int> &int_pair_map);
void PRINT_CONFIG(FRAME &SYSTEM, JOB_CONTROL & CONTROLS, int type);
void check_charges(FRAME &SYSTEM, vector<double>& TMP_CHARGES, const vector<string>& TMP_ATOMTYPE, vector<PAIR_FF> &FF_2BODY, int NATMTYP);
void parse_fcut_input(string line, vector<PAIR_FF>& FF_2BODY, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS) ;

#ifdef USE_MPI
void sync_position(vector<XYZ>& coord_vec, NEIGHBORS & neigh_list, vector<XYZ>& velocity_vec,
									 int atoms, bool sync_vel, BOX &BOXDIM) ;

void sum_forces(vector<XYZ>& accel_vec, int atoms, double &pot_energy, double &pressure,
								double &tens_xx,
								double &tens_xy,
								double &tens_xz,
								double &tens_yx,
								double &tens_yy,
								double &tens_yz,
								double &tens_zx,
								double &tens_zy,
								double &tens_zz);
#endif

#endif
