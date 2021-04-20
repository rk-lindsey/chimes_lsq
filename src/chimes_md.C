// USAGE IS: ./a.out < inputfile.in


// Setup for use with MPI

#ifdef USE_MPI
	#include <mpi.h>
#endif

#undef LOG_FORCE
#undef LOG_POS

// Setup with typical headers

#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<cstring>
#include<string>
#include<sstream>
#include<map>
#include<algorithm> // For specials vector-related functions (i.e. permute)

// Used to detect whether i/o is going to terminal or is piped... 
// will help us decide whether to use ANSI color codes

#include<unistd.h>	

// Include our own custom header

#include "functions.h"
#include "Cheby.h"
#include "util.h"
#include "io_styles.h"
#include "input.h"

// For LAMMPS linking

#if defined(USE_MPI) && defined(LINK_LAMMPS)
	#include "lmppath.h"
	using namespace LAMMPS_NS;
#endif
	
using namespace std;	
	
// Define function headers -- general

static void read_input        (string & INFILE, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);
static void read_atom_types(ifstream &PARAMFILE, JOB_CONTROL &CONTROLS, int &NATMTYP, vector<string>& TMP_ATOMTYPE, vector<int>& TMP_NATOMTYPE, vector<int>& TMP_ATOMTYPEIDX, vector<double>& TMP_CHARGES,  vector<double>& TMP_MASS, vector<int> &TMP_SIGN) ;
static void read_ff_params(ifstream &PARAMFILE, JOB_CONTROL &CONTROLS, vector<PAIR_FF>& FF_2BODY, CLUSTER_LIST& TRIPS, CLUSTER_LIST &QUADS, map<string,int> &PAIR_MAP, NEIGHBORS &NEIGHBOR_LIST, FRAME& SYSTEM, int NATMTYP, const vector<string>& TMP_ATOMTYPE, const vector<int>& TMP_ATOMTYPEIDX, vector<double> &TMP_CHARGES, vector<double> &TMP_MASS, const vector<int>& TMP_SIGN, map<int,string>& PAIR_MAP_REVERSE) ;
static void print_ff_summary(const vector<PAIR_FF> &FF_2BODY, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, const JOB_CONTROL &CONTROLS) ;

// Define function headers -- MPI

static void write_xyzv         (FRAME &SYSTEM, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL, THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBOR_LIST, string filename, bool restart);
static void read_restart_params(ifstream &COORDFILE, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL, THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBORS, FRAME &SYSTEM) ;
static void parse_ff_controls  (string &LINE, ifstream &PARAMFILE, JOB_CONTROL &CONTROLS ) ;
static void read_coord_file(int index, JOB_CONTROL &CONTROLS, FRAME &SYSTEM, ifstream &CMPR_FORCEFILE) ;
static void subtract_force(FRAME &SYSTEM, JOB_CONTROL &CONTROLS) ;
static void print_for_dftbplus(FRAME &SYSTEM, JOB_CONTROL &CONTROLS);

// Define function headers -- LAMMPS linking. Note both house_md and lammps need to be compiled for mpi use

#if defined(USE_MPI) && defined(LINK_LAMMPS)

	// Helper functions

	void Write_Lammps_datafile (const FRAME & SYSTEM, const int & NATMTYP, const vector<double> & TMP_MASS);
	void Write_Lammps_inputfile(const JOB_CONTROL & CONTROLS); 

	// callback function for LAMMPS simulation

	void md_callback(void *, bigint, int, int *, double **, double **);

#endif 


// Global variables declared as externs in functions.h, and declared in functions.C -- general

string FULL_FILE_3B;		
string SCAN_FILE_3B;
string SCAN_FILE_2B;

// Global variables declared as externs in functions.h, and declared in functions.C -- MPI calculations.   
 
int NPROCS;		// Number of processors
int RANK;		// Index of current processor

WRITE_TRAJ BAD_CONFIGS_1("XYZ","BAD_1"); // Configs where r_ij < r_cut,in 
WRITE_TRAJ BAD_CONFIGS_2("XYZ","BAD_2"); // Configs where r_ij < r_cut,in +d_penalty
WRITE_TRAJ BAD_CONFIGS_3("XYZ","BAD_3"); // All other configs, but only printed when (CONTROLS.FREQ_DFTB_GEN>0) && ((CONTROLS.STEP+1) % CONTROLS.FREQ_DFTB_GEN == 0)


// Variables that are defined locally for house_md which need to be global for LAMMPS linking

CLUSTER_LIST TRIPS;                   // Holds all 3-body parameters
CLUSTER_LIST QUADS;                   // Holds all 4-body parameters
 
#if defined(USE_MPI) && defined(LINK_LAMMPS)

	// Gloabl variables needed by callback function for LAMMPS  
	
	int TOTAL_ATOMS;
	vector<double>	LMP_CHARGE;		// Global data object to hold charges for the LAMMPS input files     
	vector<string>  TMP_ATOMTYPE;		// Will be used by lammps to map lammps atom types (int's) back to chemical symbols
	JOB_CONTROL	CONTROLS;		// Declare the data object that will hold the main simulation control variables
	map<string,int> PAIR_MAP;		// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
  	vector<int> INT_PAIR_MAP ;		
	map<int,string> PAIR_MAP_REVERSE; 	// Input/output the resverse of PAIR_MAP
	
	vector<PAIR_FF> FF_2BODY;		// Holds all 2-body parameters

	NEIGHBORS NEIGHBOR_LIST;		// Declare the class that will handle the neighbor list

#endif
	
	
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//
//						BEGIN MAIN
//
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////



int main(int argc, char* argv[])
{
  // Set up MPI if requested, otherwise, run in serial
	
#ifdef USE_MPI
  MPI_Init     (&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
  MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
		
  if (RANK==0)
	 cout << "Code compiled in MPI mode."; 
#else
  NPROCS = 1;
  RANK   = 0;
		
  if (RANK==0)
	 cout << "Code compiled in serial mode."; 
#endif

  if (RANK==0)
	 cout <<" Will run on " << NPROCS << " processor(s)." << endl;


#ifdef ENABLE_FP_EXCEPT
	enable_fp_exceptions() ;
#endif

  ////////////////////////////////////////////////////////////
  // Define/initialize important variables
  ////////////////////////////////////////////////////////////
  
  string INFILE = argv[1];
	
  FRAME      SYSTEM;			// Declare the data object that will hold the system coordinates, velocities, accelrations, etc.
  CONSTRAINT ENSEMBLE_CONTROL;		// Declare the class that will handle integration/constraints
  THERMO_AVG AVG_DATA ;

  AVG_DATA.STRESS_TENSOR_SUM_ALL.resize(3);
  
  for(int i=0;i<3; i++)
  {
  	  AVG_DATA.STRESS_TENSOR_SUM_ALL[i].X = 0.0;
  	  AVG_DATA.STRESS_TENSOR_SUM_ALL[i].Y = 0.0;
  	  AVG_DATA.STRESS_TENSOR_SUM_ALL[i].Z = 0.0;
  }

	
  // These are data objects that are normally defined locally, but need to be global for LAMMPS linking
	
#if not defined(LINK_LAMMPS)
	
  JOB_CONTROL CONTROLS;			// Declare the data object that will hold the main simulation control variables
  CONTROLS.IS_LSQ            = false ; 
  NEIGHBORS   NEIGHBOR_LIST;		// Declare the class that will handle the neighbor list
	
  // Data objects to hold coefficients for different force field types, and for FF printing (if requested)

  vector<PAIR_FF> FF_2BODY;		// Holds all 2-body parameters
		
  // Define the mapping variables that let us figure out which FF params to use for a given pair/triplet of pairs
		
  map<string,int> PAIR_MAP;		// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
  map<int,string> PAIR_MAP_REVERSE; 	// Input/output the resverse of PAIR_MAP
  vector<int> INT_PAIR_MAP ;		

  map<string,int> TRIAD_MAP	    = TRIPS.MAP ;
  map<int,string> TRIAD_MAP_REVERSE = TRIPS.MAP_REVERSE ;
  map<string,int> QUAD_MAP	    = QUADS.MAP ;
  map<int,string> QUAD_MAP_REVERSE  = QUADS.MAP_REVERSE ;


// For parameter file parsing
		
  vector<string> TMP_ATOMTYPE;

#endif 

  double Ktot;
  double Vol;
  SYSTEM.AVG_TEMPERATURE = 0.0;
		
  double dens_mol;
  double dens_mass;
  double TEMP_MASS = 0.0;
	
  string  LINE;
  string  TEMP_STR;
  int     NATMTYP = 0;
	
	
  vector<int>	TMP_SIGN;
  vector<int>	TMP_NATOMTYPE;
  vector<int> 	TMP_ATOMTYPEIDX;
  vector<double>TMP_CHARGES;
  vector<double>TMP_MASS;
  stringstream	STREAM_PARSER;
	
   string FIRST_EXT;

  ifstream CMPR_FORCEFILE;	// Holds the forces that were read in for comparison purposes
  ofstream OUT_VELOCFILE;	// Holds the velocities that are computed and are to be printed out
  ofstream STATISTICS;		// Holds a copy of the "Step      Time          Ktot/N          Vtot/N ..." part of the output file.. uninterrupted by warnings
	
  string EXTENSION;		// Holds the input xyz/xyzf file extension, so we know how many fields to read on each atom line
	
  ifstream PARAMFILE;
  ifstream COORDFILE;

  read_input(INFILE, CONTROLS, NEIGHBOR_LIST);		// Populate object with user defined values
	
  cout.precision(15);		// Set output precision
	

	
  ////////////////////////////////////////////////////////////
  // Setup optional output files 
  ////////////////////////////////////////////////////////////

  // trajectory

  WRITE_TRAJ GENFILE(CONTROLS.TRAJ_FORMAT,"STANDARD");	// This is our usual output file

  // force 
  
  WRITE_TRAJ FORCEFILE;

  if ( CONTROLS.PRINT_FORCE ) 
		FORCEFILE.INIT("XYZF_FORCE","FORCE");	// This is the forcout* files

  // velocity

  if ( CONTROLS.PRINT_VELOC ) 
	 OUT_VELOCFILE.open("velocout.txt");	
	 

  ////////////////////////////////////////////////////////////
  // Open all input files and count the total number of expected atoms
  ////////////////////////////////////////////////////////////
	
  if(!CONTROLS.BUILD)
  {
	 SYSTEM.ATOMS = 0;
	
	 for (int i=0; i<CONTROLS.COORD_FILE.size(); i++)
	 {
		COORDFILE.open(CONTROLS.COORD_FILE[i].data());
		if(!COORDFILE.is_open())
		{
		  cout << "ERROR: Cannot open coordinate file: " << CONTROLS.COORD_FILE[i] << endl;
		  exit_run(0);
		}
		else
		{
		  COORDFILE >> LINE;
		  SYSTEM.ATOMS += int(atof(LINE.data()));
		  COORDFILE.close();
		}
	 }
	

	 ////////////////////////////////////////////////////////////
	 // Set up the system data object... Do so by assuming a 
	 // constant number of atoms in the system (almost always
	 // true for md)...
	 ////////////////////////////////////////////////////////////
	
	 if (RANK==0)
		cout << "Setting up the data objects..." << endl;
		
	 SYSTEM.ATOMTYPE     .resize(SYSTEM.ATOMS);
	 SYSTEM.COORDS       .resize(SYSTEM.ATOMS);
	 SYSTEM.CHARGES      .resize(SYSTEM.ATOMS);
	 SYSTEM.MASS         .resize(SYSTEM.ATOMS);
	 SYSTEM.FORCES       .resize(SYSTEM.ATOMS);
	 SYSTEM.ACCEL        .resize(SYSTEM.ATOMS);
	 SYSTEM.VELOCITY     .resize(SYSTEM.ATOMS);
	 SYSTEM.VELOCITY_NEW .resize(SYSTEM.ATOMS);
	 SYSTEM.VELOCITY_ITER.resize(SYSTEM.ATOMS);
	 SYSTEM.ATOMTYPE_IDX .resize(SYSTEM.ATOMS);

	 if (RANK==0)
		cout << "   ...setup complete" << endl << endl;

  }

  ////////////////////////////////////////////////////////////
  // Make sure options are compatible
  ////////////////////////////////////////////////////////////

  if(CONTROLS.COORD_FILE.size() > 1)
  {
 	 if (CONTROLS.RESTART)
	 {
		cout << "ERROR: # CRDFILE # \"CAT\" option cannot be used when # VELINIT # is RESTART." << endl;
		exit_run(0);
	 } 
	 if (CONTROLS.COMPARE_FORCE)
	 {
		cout << "ERROR: # CRDFILE # \"CAT\" option cannot be used when # CMPRFRC # is true." << endl;
		exit_run(0);
	 }
	 if (CONTROLS.SUBTRACT_FORCE)
	 {
		cout << "ERROR: # CRDFILE # \"CAT\" option cannot be used when # SUBTFRC # is true." << endl;
		exit_run(0);
	 }
	 if (CONTROLS.N_LAYERS>0)
	 {
		cout << "ERROR: # CRDFILE # \"CAT\" option cannot be used when # NLAYERS # is greater than zero." << endl;
		exit_run(0);
	 }
  }

  ////////////////////////////////////////////////////////////
  // Read input file input.xyz, where box dims are on the info line:
  ////////////////////////////////////////////////////////////
	
  if(CONTROLS.BUILD)
  {
	 SYSTEM.BOXDIM.CELL_AX = SYSTEM.BOXDIM.CELL_BY = SYSTEM.BOXDIM.CELL_CZ = CONTROLS.BUILD_BOXL;
	 
	 SYSTEM.BOXDIM.UPDATE_CELL();
		
	 // Get the coordinates for the molecule to be inserted
		
	 vector<XYZ>	BASE;
	 XYZ		TMP;
	 int		TMP_NAT;
	 string		TMP_STR;
	 vector<string>	TMP_ATP;

	 if (RANK==0)
		cout << "Building initial coordinates and forces..." << endl;
			
		
	 if(CONTROLS.BUILD_TYPE == "ATOMIC")
	 {
		if (RANK==0)
		  cout << "	...Building for an atomic system." << endl;
			
		TMP.X = TMP.Y = TMP.Z = 0.0;
		BASE.push_back(TMP);
		TMP_ATP.push_back(CONTROLS.BUILD_ATOM);
	 }
	 else
	 {
		if (RANK==0)
		  cout << "	...Building for a molecular system." << endl;
			
		COORDFILE.open(CONTROLS.BUILD_FILE.data());
			
		COORDFILE >> TMP_NAT;					// Get the number of atoms
		getline(COORDFILE,LINE);				// Ignore the comment line
			
		for (int i=0; i<TMP_NAT; i++)
		{
		  COORDFILE >> TMP_STR;
		  COORDFILE >> TMP.X;
		  COORDFILE >> TMP.Y;
		  COORDFILE >> TMP.Z;
		  BASE.push_back(TMP);
		  TMP_ATP.push_back(TMP_STR);
		}
	 }
		
		
	 SYSTEM.ATOMS = CONTROLS.BUILD_NMOLEC*BASE.size();
		
	 cout << "Setting up atoms: " << SYSTEM.ATOMS << " " << CONTROLS.BUILD_NMOLEC << " " << BASE.size() << endl;
		
	 SYSTEM.ATOMTYPE     .resize(SYSTEM.ATOMS);
	 SYSTEM.COORDS       .resize(SYSTEM.ATOMS);
	 SYSTEM.CHARGES      .resize(SYSTEM.ATOMS);
	 SYSTEM.MASS         .resize(SYSTEM.ATOMS);
	 SYSTEM.FORCES       .resize(SYSTEM.ATOMS);
	 SYSTEM.ACCEL        .resize(SYSTEM.ATOMS);
	 SYSTEM.VELOCITY     .resize(SYSTEM.ATOMS);
	 SYSTEM.VELOCITY_NEW .resize(SYSTEM.ATOMS);
	 SYSTEM.VELOCITY_ITER.resize(SYSTEM.ATOMS);
	 SYSTEM.ATOMTYPE_IDX .resize(SYSTEM.ATOMS);
		
	 // Place the molecules on a staggered cubic lattice
		
	 int GRIDPOINTS = ceil(pow(CONTROLS.BUILD_NMOLEC,1.0/3.0)); 
	 int GRIDSPACE  = ceil(SYSTEM.BOXDIM.CELL_AX/GRIDPOINTS);
		
	 if (RANK==0)
		cout << "	...Using " << GRIDPOINTS << " staggered gridpoints with spacing " << GRIDSPACE << "Angstr." << endl;
			
	 int CURR_ATOM = 0;
		
	 for (int i=0; i<GRIDPOINTS; i++)
	 {
		for (int j=0; j<GRIDPOINTS; j++)
		{
		  for (int k=0; k<GRIDPOINTS; k++)
		  {
			 for (int c=0; c<BASE.size(); c++)
			 {
				TMP.X = BASE[c].X + i*GRIDSPACE;
				TMP.Y = BASE[c].Y + j*GRIDSPACE;
						
				if(j%2 == 0)
				  TMP.Z = BASE[c].Z + k*GRIDSPACE;
				else
				  TMP.Z = BASE[c].Z + k*GRIDSPACE + 0.5*GRIDSPACE;
						
				SYSTEM.ATOMTYPE[CURR_ATOM] = TMP_ATP[c];
						
				SYSTEM.COORDS[CURR_ATOM].X = TMP.X;
				SYSTEM.COORDS[CURR_ATOM].Y = TMP.Y;
				SYSTEM.COORDS[CURR_ATOM].Z = TMP.Z;
						
				// Prepare velocities
				SYSTEM.VELOCITY[CURR_ATOM].X = 0;
				SYSTEM.VELOCITY[CURR_ATOM].Y = 0;
				SYSTEM.VELOCITY[CURR_ATOM].Z = 0;
		
				// Prepare forces
				SYSTEM.FORCES[CURR_ATOM].X = 0;
				SYSTEM.FORCES[CURR_ATOM].Y = 0;
				SYSTEM.FORCES[CURR_ATOM].Z = 0;
		
				// Prepare accelerations
				SYSTEM.ACCEL[CURR_ATOM].X = 0;
				SYSTEM.ACCEL[CURR_ATOM].Y = 0;
				SYSTEM.ACCEL[CURR_ATOM].Z = 0;	
						
				CURR_ATOM++;
						
				if (CURR_ATOM == CONTROLS.BUILD_NMOLEC)
				  break;
			 }
					
			 if (CURR_ATOM == CONTROLS.BUILD_NMOLEC)
				break;
		  }
				
		  if (CURR_ATOM == CONTROLS.BUILD_NMOLEC)
			 break;
		}
			
		if (CURR_ATOM == CONTROLS.BUILD_NMOLEC)
		  break;
	 }
		
	 if (RANK==0)
		cout << "Build complete.." << endl << endl;;
			
  }	
  else
  {
	 for (int i=0; i<CONTROLS.COORD_FILE.size(); i++)
	 {
		 if ((i>0)&& (!SYSTEM.BOXDIM.IS_ORTHO))
			 EXIT_MSG("ERROR: Multiple files can only be concatenated when orthorhombic cells are used.");
		 
		 read_coord_file(i, CONTROLS, SYSTEM, CMPR_FORCEFILE) ;
	 }
	 
  }
  ////////////////////////////////////////////////////////////
  // Figure out atom charges and masses, based on parameter 
  // file
  ////////////////////////////////////////////////////////////
	
  if(RANK==0) 
	 cout << "Reading atom info from parameter file..." << endl; 
	
  // Read in the possible atom types and their features	

  read_atom_types(PARAMFILE, CONTROLS, NATMTYP, TMP_ATOMTYPE, TMP_NATOMTYPE, TMP_ATOMTYPEIDX, TMP_CHARGES, TMP_MASS, TMP_SIGN) ;

  // Assign atom features to atoms in SYSTEM data object, and the PAIR_FF object
  // Figure out how many atoms of each type we have... this can be useful for
  // printing VASP POSCAR files

  
  for(int a=0; a<SYSTEM.ATOMS;a++)
  {
  
  	for(int j=0; j<TMP_ATOMTYPE.size(); j++)
		if (SYSTEM.ATOMTYPE[a] == TMP_ATOMTYPE[j])
			TMP_NATOMTYPE[j]++;
	int i ;
	for( i=0; i<NATMTYP; i++)
	{
		if(SYSTEM.ATOMTYPE[a] == TMP_ATOMTYPE[i])
		{
				
		  SYSTEM.CHARGES[a]      = TMP_CHARGES[i];
		  SYSTEM.MASS[a]         = TMP_MASS[i];
		  SYSTEM.ATOMTYPE_IDX[a] = TMP_ATOMTYPEIDX[i];
		  break;
		}			
	}
	if ( i == NATMTYP ) 
	{
		cout << "A definition was not found for atom " << a << endl ;
		exit_run(1) ;
	}
  }

  if(RANK==0)
	 cout << "   ...read complete" << endl << endl;
	
  ////////////////////////////////////////////////////////////
  // Explicitly account for replicate images.. 
  //
  // In contrast to the old approach,  only add images in (+) direction
  // So 0 layers returns the original cell (i.e. NATOMS atoms)
  //    1 layer  returns a 8*NATOMS atoms
  //    2 layers returns  27*NATOMS atoms, and so on.
  //
  // These atoms are explicitly added to the system, so its generally
  // a good idea to use neighbor lists to cut down on cost
  ////////////////////////////////////////////////////////////	
	
  if(CONTROLS.REAL_REPLICATES>0 )
  {
  	if (!SYSTEM.BOXDIM.IS_ORTHO)
		EXIT_MSG("ERROR: Real replicates construction only implemented for orthorhombic systems.");
  
	 build_real_replicates(SYSTEM, CONTROLS) ;
  }
	
  // New layer handling: ghost atoms.. but this should be called whether or not ghost atoms are requested
  // (Required for neighbor lists, which are now ALWAYS used)

  SYSTEM.build_layers(CONTROLS.N_LAYERS);
  
  if ( (CONTROLS.N_LAYERS > 0) && (RANK == 0) )	// Then ghost atoms are used 
  {
	cout << "	Real atoms:                 " << SYSTEM.ATOMS     << endl;
	cout << "	Total atoms (ghost):        " << SYSTEM.ALL_ATOMS << endl;
	cout <<  "	cell vectors (a)            " << SYSTEM.BOXDIM.CELL_AX << " " << SYSTEM.BOXDIM.CELL_AY << " " << SYSTEM.BOXDIM.CELL_AZ << endl;
	cout <<  "	cell vectors (b)            " << SYSTEM.BOXDIM.CELL_BX << " " << SYSTEM.BOXDIM.CELL_BY << " " << SYSTEM.BOXDIM.CELL_BZ << endl;
	cout <<  "	cell vectors (c)            " << SYSTEM.BOXDIM.CELL_CX << " " << SYSTEM.BOXDIM.CELL_CY << " " << SYSTEM.BOXDIM.CELL_CZ << endl;
	cout <<  "	Layers:                     " << CONTROLS.N_LAYERS << endl;
	cout <<  "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_AX * (2*CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_AY * (2*CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_AZ * (2*CONTROLS.N_LAYERS +1) << endl;
	cout <<  "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_BX * (2*CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_BY * (2*CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_BZ * (2*CONTROLS.N_LAYERS +1) << endl;
	cout <<  "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_CX * (2*CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_CY * (2*CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_CZ * (2*CONTROLS.N_LAYERS +1) << endl;   
			
	 if(CONTROLS.WRAP_COORDS)
	 {
		cout << "WARNING: Coordinate wrapping not supported for ghost atom use. Turning option off" << endl;
		CONTROLS.WRAP_COORDS = false;
	 }
  }
  else if(RANK==0) // No ghost atoms.
	 cout << "WARNING: Ghost atoms/implicit layers are NOT being used." << endl;
	
  ////////////////////////////////////////////////////////////
  // Initialize velocities, if requested. Use the box Muller
  // method, setup up thermostat if requested, verify velocity
  // (momentum) COM is ~ zero
  ////////////////////////////////////////////////////////////

  if(RANK==0)
	 cout << "Setting up thermostats/barostats/velocities..." << endl; 

  // Set up the Nose-Hoover thermostat
	
  ENSEMBLE_CONTROL.INITIALIZE(CONTROLS.ENSEMBLE, CONTROLS, SYSTEM.ATOMS);
  
  if((ENSEMBLE_CONTROL.STYLE== "NPT-BEREND")||(ENSEMBLE_CONTROL.STYLE== "NPT-BEREND-ANISO")||(ENSEMBLE_CONTROL.STYLE== "NPT-MTK"))
  {
  	if (!SYSTEM.BOXDIM.IS_ORTHO)
		EXIT_MSG("ERROR: NPT-BEREND, NPT-BEREND-ANISO, and NPT-MTK only compatible with orthorhombic cells. Try LMP-NPT.");
}

  if ( CONTROLS.RESTART ) 
	 read_restart_params(COORDFILE, CONTROLS, ENSEMBLE_CONTROL, AVG_DATA, NEIGHBOR_LIST, SYSTEM) ;

  Vol = SYSTEM.BOXDIM.VOL;

  if ( CONTROLS.INIT_VEL ) // Use box Muller to initialize velocities
  {
	if(RANK == 0)
		cout << "Initializing velocities via box Muller..." << endl; 
	
	 ENSEMBLE_CONTROL.INIT_VEL(SYSTEM,CONTROLS);
  }
	
  if(RANK==0)
	 cout << "   ...setup complete" << endl << endl;
	
  if(RANK==0)
	 cout << "Running velocity sanity checks..." << endl; 
  
  // Test/fix velocity center of mass here

  ENSEMBLE_CONTROL.CHECK_VEL(SYSTEM,CONTROLS);
  
  ////////////////////////////////////////////////////////////
  // Print some info on density, etc. 
  ////////////////////////////////////////////////////////////  

  for(int a=0; a<SYSTEM.ATOMS; a++)
	 TEMP_MASS  += SYSTEM.MASS[a];

  dens_mol  = (SYSTEM.ATOMS * 1.0e24) / (6.0221e23 * Vol);
  dens_mass = (TEMP_MASS/Vol)*(1e24/6.0221e23);

  if(RANK==0)
  {
	 cout << "	Total Mass     = " << fixed << setprecision(4) << right << TEMP_MASS << " au"         << endl;
	 cout << "	Volume         = " << fixed << setprecision(4) << right << Vol       << " Ang^3"     << endl;
	 cout << "	Number density = " << fixed << setprecision(4) << right << dens_mol  << " mol atm/cm^3" << endl;
	 cout << "	Mass density   = " << fixed << setprecision(4) << right << dens_mass << " g/cm^3"       << endl;
	
	 cout << "   ...checks complete" << endl << endl;
  }

  ////////////////////////////////////////////////////////////
  // Begin setup of force field data objects...
  ////////////////////////////////////////////////////////////  

  read_ff_params(PARAMFILE, CONTROLS, FF_2BODY, TRIPS, QUADS, PAIR_MAP, NEIGHBOR_LIST, SYSTEM, NATMTYP, TMP_ATOMTYPE, TMP_ATOMTYPEIDX,TMP_CHARGES, TMP_MASS, TMP_SIGN, PAIR_MAP_REVERSE) ;
	
  // Set the atom charges
	
  for(int a=0; a<FF_2BODY.size(); a++)
  {
	 for(int i=0; i<NATMTYP; i++)
	 {	
		if(FF_2BODY[a].ATM1TYP == TMP_ATOMTYPE[i])
			FF_2BODY[a].ATM1CHG = TMP_CHARGES[i];
		
		if(FF_2BODY[a].ATM2TYP == TMP_ATOMTYPE[i])
			FF_2BODY[a].ATM2CHG = TMP_CHARGES[i];	
	 }	
  }

  ////////////////////////////////////////////////////////////
  // Print a summary of the force field
  ////////////////////////////////////////////////////////////  	
	
  if(RANK==0)
	 print_ff_summary(FF_2BODY, TRIPS, QUADS, CONTROLS) ;

	
  ////////////////////////////////////////////////////////////
  // Set up the neighbor list
  ////////////////////////////////////////////////////////////

  if(RANK == 0)
  	cout << "Initializing the neighbor list..." << endl;

  if (CONTROLS.USE_3B_CHEBY)
  {
   	 TRIPS.update_minmax_cutoffs(FF_2BODY);
  	 NEIGHBOR_LIST.MAX_CUTOFF_3B = TRIPS.MAX_CUTOFF;
  }
  if (CONTROLS.USE_4B_CHEBY)	       
  {
  	QUADS.update_minmax_cutoffs(FF_2BODY);
  	NEIGHBOR_LIST.MAX_CUTOFF_4B = QUADS.MAX_CUTOFF;
  }
  	
  NEIGHBOR_LIST.INITIALIZE_MD(SYSTEM);
  NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);
  
	
	
  ////////////////////////////////////////////////////////////
  // Rebuild maps into a faster data structure
  ////////////////////////////////////////////////////////////
	
	
  // Build the 2-body  and 3-body int-atomtype-ff-maps

  string int_map_str, int_map_3b_str;

  // Sanity checks on fast maps: 
	
  if(RANK==0)
  {
	 cout << endl << "Generating fast maps... " << endl;
	 cout << "	Pair maps:" << endl;
  }
	
  // Build 2-body fast maps

  build_int_pair_map(NATMTYP, TMP_ATOMTYPE, TMP_ATOMTYPEIDX, PAIR_MAP, INT_PAIR_MAP) ;
	
  // Build 3-body fast maps
  if ( FF_2BODY[0].SNUM_3B_CHEBY > 0 ) 
  {
	 TRIPS.build_int_maps(TMP_ATOMTYPE, TMP_ATOMTYPEIDX, FF_2BODY, PAIR_MAP) ;
	 //TRIPS.print(true) ;
  }
		
  // Build 4-body fast maps
  if(FF_2BODY[0].SNUM_4B_CHEBY > 0)
  {
	 QUADS.build_int_maps(TMP_ATOMTYPE, TMP_ATOMTYPEIDX, FF_2BODY, PAIR_MAP) ;
	 //QUADS.print(true) ;
  }


  Cheby cheby{CONTROLS, SYSTEM, NEIGHBOR_LIST, FF_2BODY, INT_PAIR_MAP} ;


  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  // 			LINK TO LAMMPS
  //
  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////

#if defined(USE_MPI) && defined(LINK_LAMMPS)

  if(RANK==0)
	 cout <<	"Running LAMMPS linked code! " << endl;
	
  // Note: All lammps stuff is handled by the first proc, regardless of whether serial or not, however
  //       force calculations (local) are run w/ nmpi procs
		
  // Have LAMMPS handle Ewald sums; turn off here. 
  CONTROLS.USE_COULOMB = false;
		
  // Set the total number of atoms in the system 
		
  TOTAL_ATOMS = SYSTEM.ATOMS;
		
  if(RANK==0)	// Only one rank does file writing
  {
	 // Generate the LAMMPS input file based ENTIRELY on the house_md input. 
	 // The user does not need to provide anything other than the usual chimes_md input!!!!

	 Write_Lammps_datafile (SYSTEM, NATMTYP, TMP_MASS);
	 Write_Lammps_inputfile(CONTROLS);
  }

  // Create an instance of lammps

  LAMMPS *lmp = new LAMMPS(0,NULL,MPI_COMM_WORLD);
		
  // Feed lammps the input file
	
  lmp->input->file("in.lammps");

  // set up the callback function in lammps
	 	
  int ifix = lmp->modify->find_fix_by_style("external");
  FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
  fix->set_callback(md_callback, lmp);

  // Actually run lammps... this is where the time step and simulation length are specified
  // Note:  LAMMPS output is automatically saved to log.lammps
	
  stringstream ss1,ss2;
  string tmpstr;
  const char *command;
  ss1 << "timestep " << CONTROLS.DELTA_T_FS;
  tmpstr = ss1.str();     
  command = tmpstr.c_str();
  if(RANK==0)
	 cout << command << endl;
  lmp->input->one(command);

  if      ((CONTROLS.ENSEMBLE != "LMP-MIN-BOX-ISO") && (CONTROLS.ENSEMBLE != "LMP-MIN-BOX-ANISO") && (CONTROLS.ENSEMBLE != "LMP-MIN-BOX-TRI") && (CONTROLS.ENSEMBLE != "LMP-MIN"))
	 ss2 << "run " << CONTROLS.N_MD_STEPS;
  else if ((CONTROLS.ENSEMBLE == "LMP-MIN-BOX-ISO") || (CONTROLS.ENSEMBLE == "LMP-MIN-BOX-ANISO") || (CONTROLS.ENSEMBLE == "LMP-MIN-BOX-TRI") || (CONTROLS.ENSEMBLE == "LMP-MIN"))
	 ss2 << "minimize " << CONTROLS.MIN_E_CONVG_CRIT << " " << CONTROLS.MIN_F_CONVG_CRIT << " " << CONTROLS.MIN_MAX_ITER << " " << CONTROLS.MIN_MAX_EVAL;
			
  tmpstr = ss2.str();    
		
  command = tmpstr.c_str();
  if(RANK==0)
	 cout << command << endl;
  lmp->input->one(command);
	
  if(RANK==0)
	 cout << "LAMMPS run finished." << endl;

  // exit code here; no need to run simulation code. But first, wait for all the other processes
  //(i.e. those running lammps) to get here.
		
  MPI_Barrier(MPI_COMM_WORLD);
		
  // Delete the lammps pointer
		
  delete lmp;
		
  return 0; 

#endif


  ////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////
  //
  //			START THE SIMULATION
  //
  ////////////////////////////////////////////////////////////  
  ////////////////////////////////////////////////////////////

  STATISTICS.open("md_statistics.out");
  
  if (RANK==0)
	 cout << "BEGIN SIMULATION:" << endl;
	
  int FIRST_STEP = 0;
	
  if(CONTROLS.RESTART)
	 FIRST_STEP = CONTROLS.STEP;
	
  for(CONTROLS.STEP=FIRST_STEP; CONTROLS.STEP<CONTROLS.N_MD_STEPS; CONTROLS.STEP++)	//start Big Loop here.
  {
	 ////////////////////////////////////////////////////////////
	 // Do first half of coordinate/velocity updating
	 ////////////////////////////////////////////////////////////		

	 if(CONTROLS.STEP>FIRST_STEP && RANK==0)	
	 {
		ENSEMBLE_CONTROL.UPDATE_COORDS(SYSTEM, CONTROLS);	// Update coordinates and ghost atoms
			
		if(CONTROLS.WRAP_COORDS)				// Wrap the coordinates:
		{
		 	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
				SYSTEM.BOXDIM.WRAP_ATOM(SYSTEM.COORDS[a1], SYSTEM.WRAP_IDX[a1], false);
		} 

		ENSEMBLE_CONTROL.UPDATE_VELOCS_HALF_1(SYSTEM, CONTROLS);// Update first half of velocity and max velocity for neighbor lists:		
	 }
		
#ifdef USE_MPI
	 sync_position(SYSTEM.COORDS    , NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ATOMS    , true);	// Sync the main coords. Don't need to sync veloc since only proc 1 handles constraints
	 sync_position(SYSTEM.ALL_COORDS, NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ALL_ATOMS, false);	// Sync the main coords. Don't need to sync veloc since only proc 1 handles constraints
	 MPI_Bcast(&NEIGHBOR_LIST.MAX_VEL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);				// Sync the maximum velocites so all procs use the right padding to generate their neigh lists
#endif

	 NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);

	 ////////////////////////////////////////////////////////////
	 // Calculate acceleration
	 ////////////////////////////////////////////////////////////

	 // Do the actual force calculation

	 ZCalc(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);

#ifdef USE_MPI
	// FOR MPI: Synchronize forces, energy, and pressure.
	 sum_forces(SYSTEM.ACCEL, SYSTEM.ATOMS, SYSTEM.TOT_POT_ENER, SYSTEM.PRESSURE_XYZ,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Y,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Z,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Z,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].X,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y,
	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z);
#endif

	if ( (RANK==0)&&(CONTROLS.FORDFTB ) )
		print_for_dftbplus(SYSTEM, CONTROLS);

		
	 if ( CONTROLS.CHECK_FORCE ) 
		check_forces(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST) ;
		
	// Set the econs value
	
	 ////////////////////////////////////////////////////////////
	 // Print some info on new forces, compare to input forces if 
	 // requested
	 ////////////////////////////////////////////////////////////


	 if ( CONTROLS.PRINT_FORCE && (CONTROLS.STEP+1)%CONTROLS.FREQ_FORCE == 0 && RANK == 0 ) 
		FORCEFILE.PRINT_FRAME(CONTROLS,SYSTEM);

	 if ( (CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE) ) 
	 {
		 subtract_force(SYSTEM, CONTROLS) ;
		 normal_exit() ;
	 }
	 else if (RANK == 0)	// Print main simulation header 
	 {
		////////////////////////////////////////////////////////////
		// Print the header for mains simulation output
		////////////////////////////////////////////////////////////		

		if ( CONTROLS.STEP == FIRST_STEP ) 
		{
		  printf("%8s %9s %15s %15s %15s %15s %15s", "Step", "Time", "Ktot/N", "Vtot/N", "Etot/N", "T", "P");
				
		  STATISTICS << "# Step	Time	Ktot/N	Vtot/N	Etot/N	T	P";
	  
		  if ( ENSEMBLE_CONTROL.STYLE == "NVT-MTK" || ENSEMBLE_CONTROL.STYLE == "NPT") 
		  {
			 printf(" %15s\n", "Econs/N");
			 STATISTICS << "	Econs/N" << endl;
		  }
		  else 
		  {
			 printf("\n");
			 STATISTICS << endl;
		  }
			
		  printf("%8s %9s %15s %15s %15s %15s %15s", " ", "(fs)", "(kcal/mol)", "(kcal/mol)", "(kcal/mol)", "(K)", "(GPa)");
		  STATISTICS << "#	(fs)	(kcal/mol)	(kcal/mol)	(kcal/mol)	(K)	(GPa)";
		  if ( ENSEMBLE_CONTROL.STYLE == "NVT-MTK" || ENSEMBLE_CONTROL.STYLE == "NPT") 
		  {
			 printf(" %15s\n", "(kcal/mol)");
			 STATISTICS << "	(kcal/mol)" << endl;
		  }
		  else 
		  {
			 printf("\n");
			 STATISTICS << endl;
		  }

		  std::cout.flush();
		}			
	 }
		
	 ////////////////////////////////////////////////////////////
	 // Do some thermostatting and statistics updating/output (2nd 1/2 v updates)
	 ////////////////////////////////////////////////////////////
		
	 if (RANK == 0)
	 {
		////////////////////////////////////////////////////////////
		//Convert forces to acceleration:
		////////////////////////////////////////////////////////////
		
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
		  SYSTEM.ACCEL[a1].X /= SYSTEM.MASS[a1];
		  SYSTEM.ACCEL[a1].Y /= SYSTEM.MASS[a1];
		  SYSTEM.ACCEL[a1].Z /= SYSTEM.MASS[a1];
		}
	
		////////////////////////////////////////////////////////////
		// Do second half of coordinate/velocity updating
		////////////////////////////////////////////////////////////

		if(CONTROLS.STEP>FIRST_STEP  && RANK==0)	
		  ENSEMBLE_CONTROL.UPDATE_VELOCS_HALF_2(SYSTEM, CONTROLS, NEIGHBOR_LIST);	//update second half of velocity
	 }	
		
	 ////////////////////////////////////////////////////////////
	 // Update temperature and pressure
	 ////////////////////////////////////////////////////////////

	 if (RANK == 0)
	 {
		////////////////////////////////////////////////////////////
		// Store statistics on the average simulation temperature
		// to allow for thermostat verification
		////////////////////////////////////////////////////////////

		Ktot = kinetic_energy(SYSTEM, CONTROLS);

		ENSEMBLE_CONTROL.UPDATE_TEMPERATURE(SYSTEM, CONTROLS);

		// Exit with an error if the set and block temperatures differ by more than CONTROLS.NVT_CONV_CUT
		
		AVG_DATA.TEMP_SUM += SYSTEM.TEMPERATURE;
		
		if (SYSTEM.TEMPERATURE/CONTROLS.TEMPERATURE > (1+CONTROLS.NVT_CONV_CUT) * CONTROLS.TEMPERATURE)
		{
		  cout << "ERROR: System too hot!" << endl;
		  cout << "T Set: " << CONTROLS.TEMPERATURE << endl;
		  cout << "T sys: " << SYSTEM  .TEMPERATURE << endl;
		  cout << "Step : " << CONTROLS.STEP << endl;
		  exit_run(10);
		}
		else if (CONTROLS.TEMPERATURE/SYSTEM.TEMPERATURE > (1+CONTROLS.NVT_CONV_CUT) * CONTROLS.TEMPERATURE)
		{
		  cout << "ERROR: System too cold!" << endl;
		  cout << "T Set: " << CONTROLS.TEMPERATURE << endl;
		  cout << "T sys: " << SYSTEM  .TEMPERATURE << endl;
		  cout << "Step : " << CONTROLS.STEP << endl;
		  exit_run(-10);
		}			
	 }	
	 ////////////////////////////////////////////////////////////
	 // If requested, compute pressure numerically, and accumulate
	 // statistics
	 ////////////////////////////////////////////////////////////

	 if ( CONTROLS.USE_NUMERICAL_PRESS ) 
	 {
		double PE_1, PE_2, dV; 
			
		numerical_pressure(SYSTEM, CONTROLS, FF_2BODY, TRIPS, QUADS, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST, PE_1, PE_2, dV);
			
#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &PE_1,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &PE_2,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
#endif
			
		SYSTEM.PRESSURE_XYZ = -1.0*(PE_2 - PE_1)/dV;
	 }
	 
	 
	// Set the io econs value, sync across procs
	if(RANK == 0)
	{
	 	if ( ENSEMBLE_CONTROL.STYLE == "NVT-MTK" || ENSEMBLE_CONTROL.STYLE == "NPT-MTK") 
			CONTROLS.IO_ECONS_VAL = (Ktot + SYSTEM.TOT_POT_ENER + ENSEMBLE_CONTROL.CONSERVED_QUANT(SYSTEM, CONTROLS)) / SYSTEM.ATOMS;
		else
			CONTROLS.IO_ECONS_VAL = (Ktot + SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS;
	}
	#ifdef USE_MPI
		MPI_Bcast(&CONTROLS.IO_ECONS_VAL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);					// Sync the maximum velocites so all procs use the right padding to generate their neigh lists
	#endif
	 
	 
	 	 
	 if(RANK==0)
	 {
		SYSTEM.PRESSURE = (SYSTEM.PRESSURE_XYZ + 2.0 * Ktot / (3.0 * Vol)) * GPa;	// GPa = Unit conversion factor to GPa.
			
		AVG_DATA.PRESS_SUM += SYSTEM.PRESSURE;

		SYSTEM.PRESSURE_TENSORS_ALL[0].X = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X + 2.0 / 3.0 * Ktot / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[0].Y = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Y + 2.0 / 3.0 * Ktot / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[0].Z = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Z + 2.0 / 3.0 * Ktot / Vol) * GPa;
		
		SYSTEM.PRESSURE_TENSORS_ALL[1].X = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X + 2.0 / 3.0 * Ktot / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[1].Y = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y + 2.0 / 3.0 * Ktot / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[1].Z = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Z + 2.0 / 3.0 * Ktot / Vol) * GPa;
		
		SYSTEM.PRESSURE_TENSORS_ALL[2].X = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].X + 2.0 / 3.0 * Ktot / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[2].Y = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y + 2.0 / 3.0 * Ktot / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[2].Z = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z + 2.0 / 3.0 * Ktot / Vol) * GPa;

		AVG_DATA.STRESS_TENSOR_SUM_ALL[0].X += SYSTEM.PRESSURE_TENSORS_ALL[0].X;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[0].Y += SYSTEM.PRESSURE_TENSORS_ALL[0].Y;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[0].Z += SYSTEM.PRESSURE_TENSORS_ALL[0].Z;
		
		AVG_DATA.STRESS_TENSOR_SUM_ALL[1].X += SYSTEM.PRESSURE_TENSORS_ALL[1].X;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Y += SYSTEM.PRESSURE_TENSORS_ALL[1].Y;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Z += SYSTEM.PRESSURE_TENSORS_ALL[1].Z;
		
		AVG_DATA.STRESS_TENSOR_SUM_ALL[2].X += SYSTEM.PRESSURE_TENSORS_ALL[2].X;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[2].Y += SYSTEM.PRESSURE_TENSORS_ALL[2].Y;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[2].Z += SYSTEM.PRESSURE_TENSORS_ALL[2].Z;
				
		////////////////////////////////////////////////////////////
		// Periodically print simulation output
		////////////////////////////////////////////////////////////
		
		if ( (CONTROLS.STEP+1) % CONTROLS.FREQ_ENER == 0 && RANK == 0 ) 
		{
		 	printf("%8d %9.2f %15.7f %15.7f %15.7f %15.1f %15.8f", 
			CONTROLS.STEP+1, 
			(CONTROLS.STEP+1)*CONTROLS.DELTA_T_FS, 
			Ktot/SYSTEM.ATOMS,
			SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS,
			(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS,
			SYSTEM.TEMPERATURE, 
			SYSTEM.PRESSURE);
			
			STATISTICS 
			<< CONTROLS.STEP+1
			<< "     " 
			<< (CONTROLS.STEP+1)*CONTROLS.DELTA_T_FS
			<< "  " << Ktot/SYSTEM.ATOMS
			<< "      " <<SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS
			<< "        "<<(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS
			<< " " <<SYSTEM.TEMPERATURE<< "      " 
			<< SYSTEM.PRESSURE;
			
			// Print the econs value

		  	printf("%15.7f\n", CONTROLS.IO_ECONS_VAL);
		  	STATISTICS << "        " << CONTROLS.IO_ECONS_VAL << endl;
			
		  cout.flush();
		}	
		
		////////////////////////////////////////////////////////////
		// If requested, scale the velocities
		////////////////////////////////////////////////////////////	
			
		if((ENSEMBLE_CONTROL.STYLE=="NVT-SCALE") && ((CONTROLS.STEP+1) % int(CONTROLS.FREQ_UPDATE_THERMOSTAT) == 0) )
			ENSEMBLE_CONTROL.SCALE_VELOCITIES(SYSTEM, CONTROLS);
		else
			SYSTEM.AVG_TEMPERATURE += SYSTEM.TEMPERATURE;
	 }
		
	 if((ENSEMBLE_CONTROL.STYLE=="NVT-SCALE") && ((CONTROLS.STEP+1) % int(CONTROLS.FREQ_UPDATE_THERMOSTAT) == 0) )
	 {
#ifdef USE_MPI
		sync_position(SYSTEM.COORDS    , NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ATOMS    , true);	// Sync the main coords. Don't need to sync veloc since only proc 1 handles constraints
		sync_position(SYSTEM.ALL_COORDS, NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ALL_ATOMS, false);	// Sync the main coords. Don't need to sync veloc since only proc 1 handles constraints
		MPI_Bcast(&NEIGHBOR_LIST.MAX_VEL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);					// Sync the maximum velocites so all procs use the right padding to generate their neigh lists
#endif
	 }		
	
	 ////////////////////////////////////////////////////////////
	 // If requested, write the dftbgen output file
	 ////////////////////////////////////////////////////////////
		
	 if ( (CONTROLS.FREQ_BACKUP > 0 ) && (CONTROLS.STEP+1) % CONTROLS.FREQ_BACKUP == 0 && RANK == 0) 
	 {
		rename("restart.xyzv", "restart.bak") ;
		write_xyzv(SYSTEM, CONTROLS, ENSEMBLE_CONTROL, AVG_DATA, NEIGHBOR_LIST, "restart.xyzv", true) ;
	 }

	 if ( (CONTROLS.FREQ_DFTB_GEN>0) && ((CONTROLS.STEP+1) % CONTROLS.FREQ_DFTB_GEN == 0) && RANK == 0) 
	 	GENFILE.PRINT_FRAME(CONTROLS,SYSTEM);

		
	 ////////////////////////////////////////////////////////////
	 // If requested, write the velocities
	 ////////////////////////////////////////////////////////////
		
	 if (CONTROLS.PRINT_VELOC && ((CONTROLS.STEP+1) % CONTROLS.FREQ_VELOC == 0) && RANK == 0) 
	 {
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
		  OUT_VELOCFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.VELOCITY[a1].X << endl;
		  OUT_VELOCFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.VELOCITY[a1].Y << endl;
		  OUT_VELOCFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.VELOCITY[a1].Z << endl;
		}	
	 }		
		
	 ////////////////////////////////////////////////////////////
	 ////////////////////////////////////////////////////////////
	 // 
	 // 				END LOOP OVER MD STEPS
	 // 
	 ////////////////////////////////////////////////////////////
	 ////////////////////////////////////////////////////////////
	   
  }//End big loop here.
  
  if (RANK==0)	
  {
	 cout << "END SIMULATION" << endl;

	 TEMP_MASS = 0.0;
	
	 for ( int a = 0; a < SYSTEM.ATOMS; a++ ) 
		TEMP_MASS  += SYSTEM.MASS[a];


	 cout << "	Average temperature over run = " << fixed << setprecision(4) << right << AVG_DATA.TEMP_SUM  / CONTROLS.N_MD_STEPS << " K"   << endl;
	 cout << "	Average pressure    over run = " << fixed << setprecision(4) << right << AVG_DATA.PRESS_SUM / CONTROLS.N_MD_STEPS << " GPa" << endl;
		
	// Why is this only for chebyshev?
	 if( FF_2BODY[0].PAIRTYP == "CHEBYSHEV")
	 {	 
		double Pavg = (AVG_DATA.STRESS_TENSOR_SUM_ALL[0].X + AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Y + AVG_DATA.STRESS_TENSOR_SUM_ALL[2].Z)/3.0/ CONTROLS.N_MD_STEPS ;
		cout << "	Pressures from diagonal stress tensors over run: " << Pavg << endl;
		cout << "	Average stress tensors over run: " << endl;
		cout << "		sigma_xx: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[0].X/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_yy: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Y/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_zz: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[2].Z/CONTROLS.N_MD_STEPS << " GPa" << endl; 
		cout << "		sigma_xy: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[0].Y/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_xz: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[0].Z/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_yz: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Z/CONTROLS.N_MD_STEPS << " GPa" << endl;	   
	 }

	 // Write the final configuration to file.
	 write_xyzv(SYSTEM, CONTROLS, ENSEMBLE_CONTROL, AVG_DATA, NEIGHBOR_LIST, "output.xyz", false);
			
	 STATISTICS.close();
  }

  // MPI -- End our setup
	
	normal_exit() ;
}       


static void read_input(	string & INFILE,
			JOB_CONTROL & CONTROLS,
			NEIGHBORS & NEIGHBOR_LIST)// UPDATED
{
	if (RANK==0)
		cout << endl << "Reading the simulation control input file..." << endl;
		
	// Instantiate the input class

	INPUT MD_INPUT(INFILE);
	
	if(RANK==0)
		cout << "Will read from file: " << INFILE << endl;
	
	// Read the LSQ input file
	
	MD_INPUT.PARSE_INFILE_MD(CONTROLS, NEIGHBOR_LIST);
}

static void print_for_dftbplus(FRAME &SYSTEM, JOB_CONTROL &CONTROLS)
{
	int PRINT_WIDTH     = 21; // Use 21 for testing
	int PRINT_PRECISION = 14; // Use 14 for testing
	
	ofstream dftbplus_file;
	dftbplus_file.open("dftbplus_data.dat");
	
	double offset = 0;
	
	SYSTEM.SET_NATOMS_OF_TYPE();
	
	for (int i=0; i< SYSTEM.QM_ENERGY_OFFSET.size(); i++)
		offset +=  SYSTEM.QM_ENERGY_OFFSET[i]* SYSTEM.NATOMS_OF_TYPE[i];
		
	dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.TOT_POT_ENER+offset << endl;
	
	dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X  << endl; // xx
	dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y  << endl; // yy
	dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z  << endl; // zz
	
dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X  << endl; // xy
	dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y  << endl; // xz
	dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[4].Z  << endl; // yz	
	
	for ( int ia = 0; ia < SYSTEM.ATOMS; ia++ ) 
	{
		dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.ACCEL[ia].X << endl;
		dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.ACCEL[ia].Y << endl;
		dftbplus_file << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.ACCEL[ia].Z << endl;
	}

	dftbplus_file.close();
	
	
}

static void write_xyzv(FRAME &SYSTEM, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL,
							  THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBOR_LIST, string filename, bool restart)
// Output final xyz position in the same format as input.xyz for restarting.
// Includes velocities in addition to positions.
// If restart is true, output information for job restart, including forces.
// Output full precision to aid in testing and restart.
{
	int PRINT_WIDTH     = 21; // Use 21 for testing
	int PRINT_PRECISION = 14; // Use 14 for testing
	
	// Check if the file already exists
	
	string BACKUP = filename;
	
	ifstream FILE_EXISTS(filename);
	if(FILE_EXISTS)
	{
		BACKUP.append(".bak");
		rename(filename.c_str(),BACKUP.c_str());		
	}
	FILE_EXISTS.close();
	
	ofstream fxyz;
	fxyz.open(filename);

	// Print out the file
	
	fxyz << fixed << setw(5) << setprecision(0) << SYSTEM.ATOMS << endl;
	
	if (SYSTEM.BOXDIM.IS_ORTHO)
	{
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_AX << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_BY << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_CZ;		
	}
	else
	{
		fxyz << "NON_ORTHO ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_AX << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_AY << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_AZ << " ";	
		
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_BX << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_BY << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_BZ << " ";
		
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_CX << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_CY << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.CELL_CZ;				
	}

	
	fxyz << endl;

	
	for ( int ia = 0; ia < SYSTEM.ATOMS; ia++ ) 
	{
		XYZ tmp = SYSTEM.COORDS[ia] ;
		
		if ( CONTROLS.WRAP_COORDS ) 	// Wrap into the primitive cell
			SYSTEM.BOXDIM.WRAP_ATOM(tmp, SYSTEM.WRAP_IDX[ia], false);
		
		fxyz << setw(2) << SYSTEM.ATOMTYPE[ia] << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << tmp.X << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << tmp.Y << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << tmp.Z << "    ";

		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.VELOCITY[ia].X << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.VELOCITY[ia].Y << " ";
		fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.VELOCITY[ia].Z << " ";

		if ( restart ) {
			fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.ACCEL[ia].X << " ";
			fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.ACCEL[ia].Y << " ";
			fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.ACCEL[ia].Z << endl;
		} else {
			fxyz << endl ;
		}
	}

	if ( restart ) 
	{
		fxyz << CONTROLS.STEP + 1 << endl ;
		fxyz << NEIGHBOR_LIST.RCUT_PADDING << endl ;
		ENSEMBLE_CONTROL.WRITE(fxyz) ;
		AVG_DATA.WRITE(fxyz) ;
	}		
	fxyz.close();	

}

static void read_restart_params(ifstream &COORDFILE, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL, THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBOR_LIST, FRAME &SYSTEM)
{
	vector<string> tokens;
	int ntokens;
	string LINE;

	if(COORDFILE.is_open())
		COORDFILE.close();
		
	COORDFILE.open(CONTROLS.COORD_FILE[0].data());
		
	bool FIRST = false;
	
	while (true)
	{
		getline(COORDFILE,LINE);
		ntokens = parse_space(LINE,tokens);

		if (ntokens == 1)
			if (!FIRST)
				FIRST = true;
			else
				break;
	}
	
	CONTROLS.STEP = stoi(tokens[0]);
	
	getline(COORDFILE,LINE);
	parse_space(LINE,tokens);
	
	NEIGHBOR_LIST.RCUT_PADDING = stod(tokens[0]);

	ENSEMBLE_CONTROL.READ(COORDFILE) ;
	AVG_DATA.READ(COORDFILE) ;
}


	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	//
	//			LAMMPS CALLBACK FUNCTION AND UTILS
	//
	////////////////////////////////////////////////////////////  
	////////////////////////////////////////////////////////////

#if defined(USE_MPI) && defined(LINK_LAMMPS)

	void Write_Lammps_datafile(const FRAME & SYSTEM, const int & NATMTYP, const vector<double> & TMP_MASS)
	
	// Generate an input file to sent to lammps
	// note: ensemble, timestep, and some other info are hardcoded
	// note: Velocities could be input here, after positions.
	// Format:
	// Velocities
	//
	// atom_index vx vy vz [units: Angstroms/fs, for lammps 'units real']
	// etc. 	
	{		
		ofstream LMPFILE;
		LMPFILE.open("data.lammps");
	
		LMPFILE << "# Position data file" << endl << endl;
	
		LMPFILE << SYSTEM.ATOMS << " atoms" << endl;
		LMPFILE << NATMTYP << " atom types" << endl << endl;
	
		LMPFILE << "0 " << SYSTEM.BOXDIM.CELL_LX << " xlo xhi" << endl;
		LMPFILE << "0 " << SYSTEM.BOXDIM.CELL_LY << " ylo yhi" << endl;
		LMPFILE << "0 " << SYSTEM.BOXDIM.CELL_LZ << " zlo zhi" << endl << endl;  
	
		// LAMMPS allow for a triclinic simulation cell; skew parameters included here.
		LMPFILE << SYSTEM.BOXDIM.XY << " " << SYSTEM.BOXDIM.XZ << " " << SYSTEM.BOXDIM.YZ << " xy xz yz" << endl;  
		LMPFILE << endl;  
	
		LMPFILE << "Masses" << " " << endl << endl;
		for (int i = 1; i <= NATMTYP; i++)
			LMPFILE << i << " " << TMP_MASS[i-1] << endl;
	
		LMPFILE << endl;  
	
		LMPFILE << "Atoms" << endl << endl;  
		for (int i = 0; i < SYSTEM.ATOMS; i++)
			LMPFILE << i+1 << " " << SYSTEM.ATOMTYPE_IDX[i]+1 << " " << SYSTEM.CHARGES[i] << " " << SYSTEM.COORDS[i].X << " " << SYSTEM.COORDS[i].Y << " " << SYSTEM.COORDS[i].Z << endl;  
		
		LMPFILE.close(); 

	}

	void Write_Lammps_inputfile(const JOB_CONTROL & CONTROLS)
	// Eventually we should link the pair_style cutoff distance with whatever the max forcefield cutoff is
	{
		ofstream LMPINFILE;
		LMPINFILE.open("in.lammps");
	
		// These inputs shouldn't *need* to change, so we'll hard-code them

		LMPINFILE << "units            real" << endl;
		LMPINFILE << "kspace_style     ewald 1e-4" << endl;
		LMPINFILE << "atom_style       charge" << endl;
		LMPINFILE << "atom_modify      map array" << endl;
		LMPINFILE << "atom_modify      sort 0 0.0" << endl;
		LMPINFILE << "read_data        data.lammps" << endl;
		LMPINFILE << "neighbor         1.0 bin" << endl;
		LMPINFILE << "neigh_modify     delay 0 every 1 check no page 100000 one 10000" << endl << endl;
		
		// Set the pairstyle and coefficients
	
		LMPINFILE << "pair_style       coul/long 20.0" << endl;	
		LMPINFILE << "pair_coeff       * *"  << endl << endl;
	
		// Initialize the velocity via the temperature (only used if this ISNT a minimization job)
		
		if ((CONTROLS.ENSEMBLE != "LMP-MIN-BOX-ISO") && (CONTROLS.ENSEMBLE != "LMP-MIN-BOX-ANISO") && (CONTROLS.ENSEMBLE != "LMP-MIN-BOX-TRI") && (CONTROLS.ENSEMBLE != "LMP-MIN"))
		{
			LMPINFILE << "velocity         all create " << CONTROLS.TEMPERATURE << " " << CONTROLS.SEED << endl << endl;
	
			if      (CONTROLS.ENSEMBLE == "LMP-NVE")
				LMPINFILE << "fix       1 all nve " << endl;
			else if (CONTROLS.ENSEMBLE == "LMP-NVT")
				LMPINFILE << "fix       1 all nvt temp " << CONTROLS.TEMPERATURE << " " << CONTROLS.TEMPERATURE << " " << CONTROLS.FREQ_UPDATE_THERMOSTAT << endl; 
			else if (CONTROLS.ENSEMBLE == "LMP-NPT")
				LMPINFILE << "fix       1 all npt temp " << CONTROLS.TEMPERATURE << " " << CONTROLS.TEMPERATURE << " " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " iso " << CONTROLS.PRESSURE*GPa2atm << " " << CONTROLS.PRESSURE*GPa2atm <<  " " <<  CONTROLS.FREQ_UPDATE_BAROSTAT << endl;
			else
			{
				cout << "ERROR: Unrecognized ensemble for LAMMPS linking." << endl;
				cout << endl;
				cout << "...Currently supported options are \"LMP-NVE,\" \"LMP-NVT,\" \"LMP-NPT\" \"LMP-MIN-BOX-ISO,\" \"LMP-MIN-BOX-ANISO,\" \"LMP-MIN-BOX-TRI,\" and \"LMP-MIN,\"" << endl;
				exit_run(0);
			}
		}
		else
		{
			if      (CONTROLS.ENSEMBLE == "LMP-MIN-BOX-ISO")
				LMPINFILE << "fix       1 all box/relax iso   "  << CONTROLS.PRESSURE*GPa2atm << endl;
			else if (CONTROLS.ENSEMBLE == "LMP-MIN-BOX-ANISO")
				LMPINFILE << "fix       1 all box/relax aniso "  << CONTROLS.PRESSURE*GPa2atm << endl;
			else if (CONTROLS.ENSEMBLE == "LMP-MIN-BOX-TRI")
				LMPINFILE << "fix       1 all box/relax tri   "  << CONTROLS.PRESSURE*GPa2atm << endl;
			//else if (CONTROLS.ENSEMBLE != "LMP-MIN")
			//	...No fix needed for a plain minimaztion
		}		LMPINFILE << endl;
	
		// I *think* this line converts the potential energy in the thermodynamic output section to the potential energy + thermostat/barostat contributions?
	
		LMPINFILE << "fix_modify       1 energy yes" << endl;
	
		// Not sure of the exact meaning of this syntax, but I think this is saying use the callback function to compute potential energies and forces?
	
		LMPINFILE << "fix              ext all external pf/callback 1 1" << endl;
	
		// Again, adding contributions from the callback function to the PE?
	
		LMPINFILE << "fix_modify       ext energy yes" << endl << endl;
	
		// Custom thermo output ... Make this an optional read-in string from the house_md input file eventually
	
		LMPINFILE << "thermo_style     custom step temp etotal ke pe lx ly lz pxx pyy pzz press vol density" << endl;
		LMPINFILE << "thermo_modify    format float %20.15g flush yes " << endl << endl;
	
		// Print frequencies/formats 
	
		LMPINFILE << "thermo           " << CONTROLS.FREQ_ENER << endl << endl;;
		
		// Dump a plain xyz format and in a lammps format that has boxlengths
							
		LMPINFILE << "dump             dump_1 all xyz    " << CONTROLS.FREQ_DFTB_GEN << " traj.xyz " << endl; 	// name of dump, atom group, type of dump, dump frequency, name of dumped file
		LMPINFILE << "dump             dump_2 all custom " << CONTROLS.FREQ_DFTB_GEN << " traj.lammpstrj element xu yu zu " << endl; 	// name of dump, atom group, type of dump, dump frequency, name of dumped file
	
		// Don't do any spacial ordering of molecules
	
		LMPINFILE << "atom_modify      sort 0 0.0 " << endl;
	
		// Write a restart file every 1000 steps ...Printing alternates between .a and .b
	
		LMPINFILE << "restart          100 restart.a restart.b" << endl;
		
		// Other dump options that I'm ignoring for now
		//
		// dump            1 all custom 5 pos.xyz id type q x y z
		// dump            2 all custom 5 vel.xyz id type q vx vy vz
	
	}



	void md_callback(void *ptr, bigint ntimestep, int nlocal, int *id, double **x, double **f) 
	// Callback to house_md with atom IDs and coords from each proc.
	// Invoke house_md to compute forces, load them into f for LAMMPS to use.
	// "f" can be NULL if proc owns no atoms
	{
	
		class LAMMPS *lmp = (class LAMMPS *) ptr;

		// Get information about the system that LAMMPS has built based on the lammps input file?
	    
		double boxxlo = *((double *) lammps_extract_global(lmp,"boxxlo"));
		double boxylo = *((double *) lammps_extract_global(lmp,"boxylo"));
		double boxzlo = *((double *) lammps_extract_global(lmp,"boxzlo"));
	    
		double boxxhi = *((double *) lammps_extract_global(lmp,"boxxhi"));
		double boxyhi = *((double *) lammps_extract_global(lmp,"boxyhi"));
		double boxzhi = *((double *) lammps_extract_global(lmp,"boxzhi"));
		
		double boxxy  = *((double *) lammps_extract_global(lmp,"xy"));
		double boxxz  = *((double *) lammps_extract_global(lmp,"xz"));
		double boxyz  = *((double *) lammps_extract_global(lmp,"yz"));

		double xprd = (boxxhi-boxxlo);
		double yprd = (boxyhi-boxylo);
		double zprd = (boxzhi-boxzlo);
  

		double *mass   = (double *) lammps_extract_atom(lmp,"mass"  );
		double *virial = (double *) lammps_extract_atom(lmp,"virial");
	    
		int    *type   = (int *)    lammps_extract_atom(lmp,"type"  );
		
		// Now we need to deal with the following:
		// Lammps only passes a subset of atoms to each house_md process.
		// If we push those atom subsets into ZCalc, we will get nonsense
		// because the double loop will only go over that subset of atoms
		//
		// Instead, what we do here is to construct the entire system into
		// SYS, broadcast it to all processors, and then run as usual.

		
		// Let all processors know how many atoms each processor "has"
		
		vector<int> NATOMS_PER_PROC(NPROCS);	// Figure out how many atoms each processor has
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Gather(&nlocal, 1, MPI_INT, &NATOMS_PER_PROC.front(), 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast (&NATOMS_PER_PROC.front(),NPROCS,MPI_INT,0,MPI_COMM_WORLD);

		// Create SYS with size of total number of atoms in the system
		
		int START = 0;

		if(RANK>0)
			for(int i=0; i<RANK; i++)
				START += NATOMS_PER_PROC[i]; 

		// Declare and set up the data object that will hold the system coordinates, velocities, accelrations, etc
		// For the callback force calculation
		
		
		FRAME SYS; 

		SYS.ATOMS      = TOTAL_ATOMS;//nlocal;
		
		SYS.BOXDIM.UNLAMMPSIFY(xprd, yprd, zprd, boxxy, boxxz, boxyz);
		
		SYS.ATOMTYPE    .resize(TOTAL_ATOMS);
		SYS.COORDS      .resize(TOTAL_ATOMS);
		SYS.CHARGES     .resize(TOTAL_ATOMS);
		SYS.MASS        .resize(TOTAL_ATOMS);
		SYS.FORCES      .resize(TOTAL_ATOMS);
		SYS.ACCEL       .resize(TOTAL_ATOMS);
		SYS.VELOCITY    .resize(TOTAL_ATOMS);
		SYS.VELOCITY_NEW.resize(TOTAL_ATOMS);
		SYS.ATOMTYPE_IDX.resize(TOTAL_ATOMS);
	
		for (int iter = 0; iter < nlocal; iter++) 
		{

			int i=START+iter;

  	      		SYS.CHARGES[i] = LMP_CHARGE[type[iter]-1];
		
			SYS.COORDS[i].X = x[iter][0];
			SYS.COORDS[i].Y = x[iter][1];
			SYS.COORDS[i].Z = x[iter][2];
		  
  	      		SYS.MASS[i] = mass[type[iter]];
	  
			f[iter][0] = 0;
			f[iter][1] = 0;
			f[iter][2] = 0;
	  
			SYS.ATOMTYPE_IDX[i] = type[iter];
			SYS.ATOMTYPE    [i] = TMP_ATOMTYPE[type[iter]-1];
		  
			// Wrap the coordinates
			
  	   	   	SYS.BOXDIM.WRAP_ATOM(SYS.COORDS[i], SYS.WRAP_IDX[i], true);
			
			// Prepare forces
			
			SYS.FORCES[i].X = 0;
			SYS.FORCES[i].Y = 0;
			SYS.FORCES[i].Z = 0;
		  
			// Prepare accelerations
		  
			SYS.ACCEL[i].X = 0;
			SYS.ACCEL[i].Y = 0;
			SYS.ACCEL[i].Z = 0;	
		}
	    
		// Now we need to sync SYS.COORDS across all processors.. 
		// Note, masses and velocities have not been updated because they 
		// do not factor into anything in house_md but kinetic energy calculations, which 
		// are not used in conjunction with LAMMPS

	    
		vector<int> STARTS          (NPROCS);
		vector<int> NCOORDS_PER_PROC(NPROCS);
		vector<int> STARTS_COORDS   (NPROCS);
	    
		SYS.MY_ATOMS       = nlocal;
		SYS.MY_ATOMS_START = START; 

		for(int i=0; i<NPROCS; i++)
		{
			STARTS[i]           = 0;
			STARTS_COORDS[i]    = 0;
			NCOORDS_PER_PROC[i] = 3*NATOMS_PER_PROC[i];	
			
			for(int j=0; j<i; j++)
			{
				STARTS[i]        += NATOMS_PER_PROC[j];
				STARTS_COORDS[i] += 3*NATOMS_PER_PROC[j];
			}
		}
	    
	    
		double *coord = (double *) SYS.COORDS.data();
		double *force = (double *) SYS.FORCES.data();
		double *accel = (double *) SYS.ACCEL .data();
		
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &SYS.CHARGES.front(),      &NATOMS_PER_PROC.front(),       &STARTS.front(),        MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &SYS.ATOMTYPE_IDX.front(), &NATOMS_PER_PROC.front(),       &STARTS.front(),        MPI_INT,    MPI_COMM_WORLD);
		MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, coord,                     &NCOORDS_PER_PROC.front(),      &STARTS_COORDS.front(), MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, force,                     &NCOORDS_PER_PROC.front(),      &STARTS_COORDS.front(), MPI_DOUBLE, MPI_COMM_WORLD);
		MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, accel,                     &NCOORDS_PER_PROC.front(),      &STARTS_COORDS.front(), MPI_DOUBLE, MPI_COMM_WORLD);
		
		for(int i=0; i<TOTAL_ATOMS; i++)
			SYS.ATOMTYPE[i] = TMP_ATOMTYPE[SYS.ATOMTYPE_IDX[i]-1];
		
		// Now we need to build the ghost atoms/neighbor lists
		
		SYS.build_layers(CONTROLS.N_LAYERS) ;
		
		// Use very little padding because we will update neighbor list for every frame.
		
		double PADDING = 0.01;
		
		NEIGHBOR_LIST.INITIALIZE(SYS, PADDING);	// Doesn't care about velocity
		NEIGHBOR_LIST.DO_UPDATE(SYS, CONTROLS);	
		
		// Do the actual force calculation using our MD code

		ZCalc(SYS, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);

		// Now we need to sync up all the forces, potential energy, and tensors with LAMMPS
		//... but first, we need to sum the potential energy computed by each process, since 
		// we're saving it directly to LAMMPS' PE which corresponds to the sum over all processes.
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Allreduce(MPI_IN_PLACE, &SYS.TOT_POT_ENER,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		
		// Also, need to add up the forces (accel) to account for f[i]+=val, f[j]-=val, for when 
		// j is on a different proc than i
		
		MPI_Allreduce(MPI_IN_PLACE, accel, 3*SYS.ATOMS, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
		// Give LAMMPS the computed forces
		
		for (int iter = 0; iter < nlocal; iter++)
		{
			int i=START+iter;
		
			f[iter][0] = SYS.ACCEL[i].X;
			f[iter][1] = SYS.ACCEL[i].Y;
			f[iter][2] = SYS.ACCEL[i].Z;
		}
		
		// Integrate?
		
		int ifix = lmp->modify->find_fix_by_style("external");
		FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
		fix->set_energy(SYS.TOT_POT_ENER);
		
		// stress tensor is multiplied by (atm/vol) in LAMMPS; multiply out volume here
		// virial points to LAMMPS array for stress tensor
		// XY, YZ, and XZ components need to be included at a later date.
	    
		// Set the stresses based on our MD code's calculation
	    
		double vol = SYS.BOXDIM.VOL;

		virial[0] =  SYS.PRESSURE_TENSORS_XYZ.X*vol;
		virial[1] =  SYS.PRESSURE_TENSORS_XYZ.Y*vol;
		virial[2] =  SYS.PRESSURE_TENSORS_XYZ.Z*vol;
	    
		virial[3] = SYS.PRESSURE_TENSORS_XYZ[0].Y*vol; // XY
		virial[4] = SYS.PRESSURE_TENSORS_XYZ[0].Z*vol; // XZ
		virial[5] = SYS.PRESSURE_TENSORS_XYZ[1].Z*vol; // YZ
    
	}
	
#endif

static void read_atom_types(ifstream &PARAMFILE, JOB_CONTROL &CONTROLS, int &NATMTYP, vector<string>& TMP_ATOMTYPE,
									 vector<int>& TMP_NATOMTYPE, vector<int>& TMP_ATOMTYPEIDX, vector<double>& TMP_CHARGES, 
									 vector<double>& TMP_MASS, vector<int> &TMP_SIGN)
// Read in the atom types.  Set properties for each atom type.
{

  bool FOUND_END = false ;
  string 	LINE;
  stringstream	STREAM_PARSER;
  string TEMP_STR ;

  PARAMFILE.open(CONTROLS.PARAM_FILE.data());
  if(!PARAMFILE.is_open())
  {
	 cout << "ERROR: Cannot open paramter file: " << CONTROLS.PARAM_FILE << endl;
	 exit_run(0);
  }
	
  FOUND_END = false;
  NATMTYP = 0;
	
  while (FOUND_END == false)
  {
	 getline(PARAMFILE,LINE);
	
	 if(LINE.find("# PAIRIDX #") != string::npos)
	 {
		FOUND_END = true;
			
		if(NATMTYP == 0)
		{
		  cout << "ERROR: ATOM TYPES section not found in parameter file: " << CONTROLS.PARAM_FILE << endl;
		  exit_run(0);
		}
			
		PARAMFILE.close();
		break;
	 }
		
	 else if(LINE.find("ATOM TYPES:") != string::npos)
	 {
		STREAM_PARSER.str(LINE);
		STREAM_PARSER >> TEMP_STR >> TEMP_STR >> NATMTYP;
		CONTROLS.NATMTYP = NATMTYP;
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();
		
		CONTROLS.ATOMTYPES.resize(NATMTYP);
		
		if ( NATMTYP > MAX_ATOM_TYPES ) 
		{
		  cout << "ERROR: TOO MANY ATOM TYPES DEFINED\n" ;
		  exit_run(0) ;
		}
  
		// We need a globally-defined data object to hold charges for passing into LAMMPS
#if defined(USE_MPI) && defined(LINK_LAMMPS)
		LMP_CHARGE         .resize(NATMTYP);
#endif
			
		TMP_ATOMTYPE   .resize(NATMTYP);
		TMP_NATOMTYPE  .resize(NATMTYP);
		TMP_ATOMTYPEIDX.resize(NATMTYP);
		TMP_CHARGES    .resize(NATMTYP);
		TMP_MASS       .resize(NATMTYP);	
		TMP_SIGN	   .resize(NATMTYP);
			
		for(int i=0; i<TMP_NATOMTYPE.size(); i++)
		  TMP_NATOMTYPE[i] = 0;
			
		if(RANK==0)
		  cout << "	Read " << NATMTYP << " atom types:" << endl;			
	 }	
		
	 // Quickly figure out if we're fitting charges... note this process gets repeated more 
	 // comprehensively below.. Here it's just needed to figure out charge signs
		
	 else if(LINE.find("USECOUL: ") != string::npos) 
	 {
		STREAM_PARSER.str(LINE);
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER.str(""); 
		STREAM_PARSER.clear();

		PARAMFILE >> TEMP_STR >> TEMP_STR;
		if (TEMP_STR=="true"  || TEMP_STR=="True"  || TEMP_STR=="TRUE"  || TEMP_STR == "T" || TEMP_STR == "t")
		  CONTROLS.FIT_COUL = true;
		else
		  CONTROLS.FIT_COUL = false;

		PARAMFILE.ignore();
	 }				
		
	 else if(LINE.find("# TYPEIDX #") != string::npos)
	 {
		for (int i=0; i<NATMTYP; i++)
		{
		  getline(PARAMFILE,LINE);
		  STREAM_PARSER.str(LINE);
				
		  STREAM_PARSER >> TEMP_STR;
		  TMP_ATOMTYPEIDX[i] = stoi(TEMP_STR) ;
		  STREAM_PARSER >> TMP_ATOMTYPE[i];
		  CONTROLS.ATOMTYPES[i] = TMP_ATOMTYPE[i];
		  if(CONTROLS.FIT_COUL)
			 STREAM_PARSER >> TEMP_STR;
		  else
			 STREAM_PARSER >> TMP_CHARGES[i];
		  STREAM_PARSER >> TMP_MASS[i];
				
		  STREAM_PARSER.str("");
		  STREAM_PARSER.clear();
				
		  if(RANK==0)
			 cout << "		" 	<< setw(5) << left << i << " " << setw(2) << left << TMP_ATOMTYPE[i] << ", q (e): ";
				
		  if(CONTROLS.FIT_COUL)
		  {
			 if(TEMP_STR == "+")
			 {
				if(RANK==0)
				  cout << "POSITIVE";
				TMP_SIGN[i] = 1;
			 }
			 else
			 {
				if(RANK==0)
				  cout << "NEGATIVE";
				TMP_SIGN[i] = -1;
			 }
		  }
		  else
			 if(RANK==0)
				cout << setw(6) << fixed << setprecision(3) << right << TMP_CHARGES[i];
		  if(RANK==0)
			 cout  << ", mass (amu): " << setw(8) << fixed << setprecision(4) << right << TMP_MASS[i] << endl;
		}
	 }	
  }
}

static void read_ff_params(	ifstream & PARAMFILE, 
			   	JOB_CONTROL &CONTROLS, 
			   	vector<PAIR_FF> & FF_2BODY,
				CLUSTER_LIST& TRIPS, 
				CLUSTER_LIST &QUADS, 
				map<string,int> &PAIR_MAP, 
				NEIGHBORS &NEIGHBOR_LIST, 
				FRAME& SYSTEM, 
				int NATMTYP,
                           	const vector<string>& TMP_ATOMTYPE, 
				const vector<int>& TMP_ATOMTYPEIDX,
				vector<double> &TMP_CHARGES, 
				vector<double> &TMP_MASS,
				const vector<int>& TMP_SIGN, 
				map<int,string>& PAIR_MAP_REVERSE)
				
// Read in the force field parameters.
{
  bool   	FOUND_END = false;
  string 	LINE;
  stringstream	STREAM_PARSER;
  string TEMP_STR ;

  PARAMFILE.open(CONTROLS.PARAM_FILE.data());
  if(!PARAMFILE.is_open())
  {
	 cout << "ERROR: Cannot open parameter file: " << CONTROLS.PARAM_FILE << endl;
	 exit_run(0);
  }
	
  if(RANK==0)
	 cout << endl << "Reading force field parameters..." << endl;
	
  FOUND_END = false;

  string  TEMP_SEARCH_2B = "some default text";
  string  TEMP_SEARCH_3B = "some default text";
  string  TEMP_SEARCH_4B = "some default text";
  string	TEMP_TYPE;
  int     NO_PAIRS, NO_TRIPS, NO_QUADS;
  int		TMP_TERMS1{0}, TMP_TERMS2{0}, TMP_TERMS3{0};
  double	TMP_LOW  = -1;
  double 	TMP_HIGH =  1;
  const double eps_charge = 1.0e-04 ;

  while (FOUND_END == false)
  {
	 LINE = get_next_line(PARAMFILE) ;

	 // Break out of loop

	 if(LINE.find("ENDFILE") != string::npos)
	 {			
		// Rewind so we can set the special 3-body cutoffs
			
		PARAMFILE.seekg(0);
			
		while (FOUND_END == false)		
		{
		  LINE = get_next_line(PARAMFILE) ;

		  if(LINE.find("ENDFILE") != string::npos)
			 break;	
				
		  if(LINE.find("FCUT TYPE:") != string::npos )
		  {
				// Unified parsing of fcut for MD and LSQ code.
			
				LINE = LINE.substr(11) ; // Chop off "FCUT TYPE: ".
				parse_fcut_input(LINE,FF_2BODY,TRIPS,QUADS) ;
		  }
				
		  else if(LINE.find("SPECIAL 3B S_MINIM:") != string::npos)
		  {
				TRIPS.read_cutoff_params(PARAMFILE, LINE, "S_MINIM") ;
				TRIPS.process_cutoff_params("S_MINIM", FF_2BODY, PAIR_MAP) ;
		  }

		  else if(LINE.find("SPECIAL 3B S_MAXIM:") != string::npos)
		  {
				TRIPS.read_cutoff_params(PARAMFILE, LINE, "S_MAXIM") ;
				TRIPS.process_cutoff_params("S_MAXIM", FF_2BODY, PAIR_MAP) ;
//NEIGHBOR_LIST.MAX_CUTOFF_3B = TRIPS.MAX_CUTOFF ;

		  }
				
		  else if(LINE.find("SPECIAL 4B S_MINIM:") != string::npos)
		  {
				QUADS.read_cutoff_params(PARAMFILE, LINE, "S_MINIM") ;
				QUADS.process_cutoff_params("S_MINIM",FF_2BODY, PAIR_MAP) ;
		  }
				
		  else if(LINE.find("SPECIAL 4B S_MAXIM:") != string::npos)
		  {
				QUADS.read_cutoff_params(PARAMFILE, LINE, "S_MAXIM") ;
				QUADS.process_cutoff_params("S_MAXIM",FF_2BODY, PAIR_MAP) ;
				//NEIGHBOR_LIST.MAX_CUTOFF_4B = QUADS.MAX_CUTOFF ;
		  }
		  else if ( LINE.find("NO ENERGY OFFSETS:") != string::npos)
		  {
				vector<string> tokens;
				int ntokens ;
				ntokens = parse_space(LINE,tokens) ;
				
				// Determine the number of offsets (atom types)
				SYSTEM.QM_ENERGY_OFFSET.resize(stoi(tokens[3])); 
				
				for (int i=0; i<SYSTEM.QM_ENERGY_OFFSET.size(); i++)
				{
					LINE = get_next_line(PARAMFILE);
					ntokens = parse_space(LINE,tokens);
				
					if ( ntokens == 4 )
					{
						SYSTEM.QM_ENERGY_OFFSET[i] = stod(tokens[3]) ;
						if ( RANK == 0 ) cout << "QM Energy Offset " << i << " = " << SYSTEM.QM_ENERGY_OFFSET[i] << endl ;
					}
					else
						EXIT_MSG("The ENERGY OFFSET: command in the params file could not be read for offset no", i) ;				
				}
			}
		}

		if(RANK==0)
		{
		  cout << "   ...read complete." << endl << endl;
		  cout << "Notes on simulation: " << endl;
		  
		  cout << "	Using the following fcut style for 2B Chebyshev interactions: ";
		  FF_2BODY[0].FORCE_CUTOFF.print_params();	
				
		  if( CONTROLS.USE_3B_CHEBY )
			 TRIPS.print_fcut() ;

		  if ( CONTROLS.USE_4B_CHEBY )
			 QUADS.print_fcut() ;

		}

		if(CONTROLS.USE_COULOMB)
		{
		  if(RANK==0)
		  {
			 cout << "	Ewald summations will be used ";
			 if(CONTROLS.FIT_COUL)
				cout << "and charges will be taken from fit values" << endl;
			 else
				cout << "and charges will be taken from user-specified values" << endl;
		  }

		}
		else
		{
		  if(RANK==0)
			 cout << "	Electrostatics will not be computed." << endl;
		}
			
		if(RANK==0)
		  cout << endl;  
		
		if(SYSTEM.QM_ENERGY_OFFSET.size() == 0)
		{
			SYSTEM.QM_ENERGY_OFFSET.resize(NO_PAIRS);
			
			for(int i=0; i<NO_PAIRS; i++)
				SYSTEM.QM_ENERGY_OFFSET[i] = 0.0;
		}
			
		PARAMFILE.close();
		break;
	 }
		
	 // Determine what parameters we're actually reading
		
	 else if(LINE.find("USECOUL: ") != string::npos)
	 {
		parse_ff_controls(LINE, PARAMFILE, CONTROLS) ;
	 }
		
	 // Determine the pair type and corresponding orders, etc
		
	 else if(LINE.find("PAIRTYP: ") != string::npos)
	 {

		if(RANK==0)
		  cout << "Attempting to read pairtype..." << endl;
			
		STREAM_PARSER.str(LINE);
			
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> TEMP_TYPE;
			
		if(TEMP_TYPE =="CHEBYSHEV")
		{
			  
		  vector<string> tokens ;
		  string buf ;

		  // Tokenize the input.
		  while ( STREAM_PARSER >> buf ) 
			 tokens.push_back(buf) ;

		  TMP_TERMS1 = 0 ;
		  TMP_TERMS2 = 0 ;
		  TMP_TERMS3 = 0 ;
		  if ( tokens.size() > 0 )
			 TMP_TERMS1 = stoi(tokens[0]) ;
		  if ( tokens.size() > 1 )
			 TMP_TERMS2 = stoi(tokens[1]) ;
		  if ( tokens.size() > 2 )
			 TMP_TERMS3 = stoi(tokens[2]) ;

		  if ( tokens.size() > 3 )
		  {
			 TMP_LOW = stoi(tokens[3]) ;
			 if( TMP_LOW < -1.0 ||  TMP_LOW > +1.0 )
			 {
				cout << "ERROR: CHEBY_RANGE_LOW must be betwee -1 and 1" << endl;
				exit_run(0);
			 }
		  }
				
				
		  if ( tokens.size() > 4 )
		  {
			 TMP_HIGH = stoi(tokens[4]) ;
			 if( TMP_HIGH < -1.0 ||  TMP_HIGH > +1.0 )
			 {
				cout << "ERROR: CHEBY_RANGE_HIGH must be betwee -1 and 1" << endl;
				exit_run(0);
			 }
		  }

		  if(TMP_LOW > TMP_HIGH || TMP_LOW == TMP_HIGH)
		  {
			 cout << "ERROR: CHEBY_RANGE_LOW must be smaller than CHEBY_RANGE_HIGH" << endl;
			 exit_run(0);
		  }

		}			
			
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();	
			
		if(RANK==0)
		  cout << "	...Read FF interaction type..." << endl;	
	 }

	 // Read in pair potential info

	 else if(LINE.find("ATOM PAIRS: ") != string::npos)
	 {	
		// Determine the number of pairs

		STREAM_PARSER.str(LINE);
			
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> NO_PAIRS;
			
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();	
			
		// Resize pair parameter object
			
		FF_2BODY.resize(NO_PAIRS);		
			
		// Set any defaults for that pair style
			
		if(TEMP_TYPE == "CHEBYSHEV")
		{
		  for(int i=0; i<NO_PAIRS; i++)
		  {
			 FF_2BODY[i].PENALTY_DIST     = 0.01;
			 FF_2BODY[i].PENALTY_SCALE    = 1.0e4;
			 FF_2BODY[i].CHEBY_RANGE_HIGH = TMP_HIGH;
			 FF_2BODY[i].CHEBY_RANGE_LOW  = TMP_LOW;
		  }

		}
			
		// Read in general pair parameters
			
		LINE = get_next_line(PARAMFILE) ;
		LINE = get_next_line(PARAMFILE) ;

		for(int i=0; i<NO_PAIRS; i++)
		{
		  FF_2BODY[i].PAIRTYP = TEMP_TYPE;
	
		  PARAMFILE >> FF_2BODY[i].PAIRIDX;
		  PARAMFILE >> FF_2BODY[i].ATM1TYP;
		  PARAMFILE >> FF_2BODY[i].ATM2TYP;	
		  PARAMFILE >> FF_2BODY[i].S_MINIM;	
		  PARAMFILE >> FF_2BODY[i].S_MAXIM;

		  // Set the atom type index.
		  FF_2BODY[i].ATM1TYPE_IDX = -1 ;
		  FF_2BODY[i].ATM2TYPE_IDX = -1 ;
		  for ( int j = 0 ; j < TMP_ATOMTYPE.size() ; j++ ) 
		  {
			 if ( FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] ) 
				FF_2BODY[i].ATM1TYPE_IDX = TMP_ATOMTYPEIDX[j] ;
			 if ( FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j] ) 
				FF_2BODY[i].ATM2TYPE_IDX = TMP_ATOMTYPEIDX[j] ;
		  }
		  if ( FF_2BODY[i].ATM1TYPE_IDX < 0 )
			 EXIT_MSG("While reading the FF file, could not find a match for atom 1 in atom pair " + FF_2BODY[i].ATM1TYP ) ;

		  if ( FF_2BODY[i].ATM2TYPE_IDX < 0 )
			 EXIT_MSG("While reading the FF file, could not find a match for atom 2 in atom pair " + FF_2BODY[i].ATM2TYP ) ;
		  
		  if(FF_2BODY[i].S_MAXIM > NEIGHBOR_LIST.MAX_CUTOFF)
		  {
			 NEIGHBOR_LIST.MAX_CUTOFF    = FF_2BODY[i].S_MAXIM;
		  }

		if ((! SYSTEM.BOXDIM.IS_RCUT_SAFE(FF_2BODY[i].S_MAXIM, CONTROLS.N_LAYERS)))  
		{
			if (isatty(fileno(stdout)) && RANK == 0)
			{
				#if WARN == TRUE
					cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "WARNING: ";
				#else
					cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "Error: ";
				#endif
				
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "Outer cutoff greater than half of at least one layered cell vector at least one box length: "  << FF_2BODY[i].S_MAXIM <<COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Frame:                      " << i << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Pair type:                  " << FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	cell vectors (a)            " << SYSTEM.BOXDIM.CELL_AX << " " << SYSTEM.BOXDIM.CELL_AY << " " << SYSTEM.BOXDIM.CELL_AZ << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	cell vectors (b)            " << SYSTEM.BOXDIM.CELL_BX << " " << SYSTEM.BOXDIM.CELL_BY << " " << SYSTEM.BOXDIM.CELL_BZ << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	cell vectors (c)            " << SYSTEM.BOXDIM.CELL_CX << " " << SYSTEM.BOXDIM.CELL_CY << " " << SYSTEM.BOXDIM.CELL_CZ << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Layers:                     " << CONTROLS.N_LAYERS << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_AX * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_AY * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_AZ * (CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_BX * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_BY * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_BZ * (CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_CX * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_CY * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_CZ * (CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
			
				#if ! WARN
					exit_run(0);
				#endif
			}
			else if ( RANK == 0 ) 
			{
				#if WARN == TRUE
					cout << "WARNING: ";
				#else
					cout << "Error: ";
				#endif
			
				cout <<  "Outer cutoff greater than half of at least one layered cell vector at least one box length: "  << FF_2BODY[i].S_MAXIM <<COUT_STYLE.ENDSTYLE << endl;
				cout <<  "	Frame:                      " << i << COUT_STYLE.ENDSTYLE << endl;
				cout <<  "	Pair type:                  " << FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
				cout <<  "	cell vectors (a)            " << SYSTEM.BOXDIM.CELL_AX << " " << SYSTEM.BOXDIM.CELL_AY << " " << SYSTEM.BOXDIM.CELL_AZ << endl;
				cout <<  "	cell vectors (b)            " << SYSTEM.BOXDIM.CELL_BX << " " << SYSTEM.BOXDIM.CELL_BY << " " << SYSTEM.BOXDIM.CELL_BZ << endl;
				cout <<  "	cell vectors (c)            " << SYSTEM.BOXDIM.CELL_CX << " " << SYSTEM.BOXDIM.CELL_CY << " " << SYSTEM.BOXDIM.CELL_CZ << endl;
				cout <<  "	Layers:                     " << CONTROLS.N_LAYERS << endl;
				cout <<  "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_AX * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_AY * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_AZ * (CONTROLS.N_LAYERS +1) << endl;
				cout <<  "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_BX * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_BY * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_BZ * (CONTROLS.N_LAYERS +1) << endl;
				cout <<  "	Effective cell vectors (a): " << SYSTEM.BOXDIM.CELL_CX * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_CY * (CONTROLS.N_LAYERS +1) << " " << SYSTEM.BOXDIM.CELL_CZ * (CONTROLS.N_LAYERS +1) << endl;	
			
				#if ! WARN
					exit_run(0);
				#endif										
			}
		}
			 				
		  FF_2BODY[i].PRPR_NM = FF_2BODY[i].ATM1TYP;
		  FF_2BODY[i].PRPR_NM.append(FF_2BODY[i].ATM2TYP);	

		  if(TEMP_TYPE =="CHEBYSHEV")
		  {
			 string cheby_type ;
			 PARAMFILE >> cheby_type ;

			 FF_2BODY[i].CHEBY_TYPE = Cheby::get_trans_type(cheby_type) ;

			 if(FF_2BODY[i].CHEBY_TYPE == Cheby_trans::MORSE)
				PARAMFILE >> FF_2BODY[i].LAMBDA;	
					
			 FF_2BODY[i].SNUM          = TMP_TERMS1;
			 FF_2BODY[i].SNUM_3B_CHEBY = TMP_TERMS2;
			 FF_2BODY[i].SNUM_4B_CHEBY = TMP_TERMS3;
					
			 // Set Cheby parameters based on S_MINIM and S_MAXIM.
			 FF_2BODY[i].set_cheby_vals() ;

			 // Setup force cutoff type for 2-body
					
			 FF_2BODY[0].FORCE_CUTOFF.set_type("CUBIC");
					
			 // Copy all class members.
			 for(int i=1; i<FF_2BODY.size(); i++)
				FF_2BODY[i].FORCE_CUTOFF = FF_2BODY[0].FORCE_CUTOFF;										
		  }
		  else // Unknown type
		  {
			 cout << "ERROR: Unknown type: " << TEMP_TYPE << endl; 
		  }				
		}

		TEMP_SEARCH_2B = "PAIR ";
		TEMP_SEARCH_2B.append(FF_2BODY[0].PAIRTYP); // Syntax ok b/c all pairs have same FF type
		TEMP_SEARCH_2B.append(" PARAMS");	
			
		if(RANK==0)
		  cout << "	...Read general FF params..." << endl;
	 }
		
	 // Read any special controls for the pair potential

	 else if(LINE.find("PAIR CHEBYSHEV PENALTY DIST: ") != string::npos)
	 {
		STREAM_PARSER.str(LINE);
		STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR;
		for(int i=0; i<NO_PAIRS; i++) 
		{
		  FF_2BODY[i].PENALTY_DIST = stod(TEMP_STR) ;
		}
		if ( RANK == 0 ) 
		  cout << "PAIR CHEBYSHEV PENALTY DIST: " << FF_2BODY[0].PENALTY_DIST << endl ;

		STREAM_PARSER.str("");
		STREAM_PARSER.clear();	
	 }
		
	 else if(LINE.find("PAIR CHEBYSHEV PENALTY SCALING: ") != string::npos)
	 {
		STREAM_PARSER.str(LINE);
		STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR;
		for(int i=0; i<NO_PAIRS; i++)
		{
		  FF_2BODY[i].PENALTY_SCALE = stod(TEMP_STR) ;
		}
		if ( RANK == 0 ) 
		  cout << "PAIR CHEBYSHEV PENALTY SCALING: " << FF_2BODY[0].PENALTY_SCALE << endl ;
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();	
	 }
				
	 // Setup triplets
		
	 else if(LINE.find("ATOM PAIR TRIPLETS: ") != string::npos)
	 {	
		// Determine the number of pairs

		STREAM_PARSER.str(LINE);
			
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> NO_TRIPS;
			
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();	
			
		// Resize pair parameter object
			
		TEMP_SEARCH_3B = TRIPS.allocate(NO_TRIPS,3,FF_2BODY) ;
			
		if(RANK==0)
		  cout << "	...Read FF triplet specifications..." << endl;
			
	 }
		
	 // Setup quadruplets
		
	 else if(LINE.find("ATOM PAIR QUADRUPLETS: ") != string::npos)
	 {	
		// Determine the number of pairs

		STREAM_PARSER.str(LINE);
			
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> NO_QUADS;
			
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();	
			
		// Resize pair parameter object

		TEMP_SEARCH_4B = QUADS.allocate(NO_QUADS, 4, FF_2BODY) ;
			
		if(RANK==0)
		  cout << "	...Read FF quadruplet specifications..." << endl;
			
	 }
		
	 // Read in pair parameters
		
	 else if(LINE.find(TEMP_SEARCH_2B) != string::npos) // "PAIR <PAIRTYPE> PARAMS"
	 {	

		// Read in the specific atom pair parameters
			
		if(RANK==0)
		  cout << "	...Reading all remaining force field parameters..." << endl;

		for(int i=0; i<NO_PAIRS; i++)
		{
		  LINE = get_next_line(PARAMFILE) ; // Blank line
		  LINE = get_next_line(PARAMFILE) ; // "PAIRTYPE PARAMS: <index> <atom1> <atom2>"	
		  LINE = get_next_line(PARAMFILE) ; // Blank line
				
		  FF_2BODY[i].PARAMS.resize(FF_2BODY[i].SNUM);
				
		  for(int j=0; j<FF_2BODY[i].SNUM; j++)
		  {
			 // Read the short range potential parameters
					
			 PARAMFILE >> TEMP_STR;
			 PARAMFILE >> FF_2BODY[i].PARAMS[j];
		  }
								
		  // If applicable, read the charge parameter
				
		  if(CONTROLS.FIT_COUL)
		  {
			 PARAMFILE >> TEMP_STR >> TEMP_STR >> TEMP_STR >> FF_2BODY[i].PAIR_CHRG;
			 // NOTE: Fit charges are in stillinger units.. convert to e:
			 FF_2BODY[i].PAIR_CHRG /= ke;
		  }
		}	
			
		// At this point we have all the info needed to set the charges...
		//
		// Case 1: FITCOUL is false. Take charges from the individual charges specified by the user
		//        ("ATOM TYPES/TYPEIDX" section of the parameter file)	
		//
		// Case 2: FITCOUL is true. Determine individual charges from charges of type Qxx

			
		if(!CONTROLS.FIT_COUL)
		{
		  for(int i=0; i<NO_PAIRS; i++)
		  {
			 for(int j=0; j<NATMTYP; j++)
			 {
				if( FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] )
				  FF_2BODY[i].ATM1CHG =  TMP_CHARGES [j];
					
				if( FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j] )
				  FF_2BODY[i].ATM2CHG =  TMP_CHARGES [j];
			 }	
				
			 FF_2BODY[i].PAIR_CHRG = FF_2BODY[i].ATM1CHG*FF_2BODY[i].ATM2CHG;
				
		  }
		}
		else
		{				
		  if ( RANK == 0 )
			 cout << "		Re-setting individual atom charges based on pair charges :" << endl;
				
		  for(int j=0; j<NATMTYP; j++)
		  {
			 for(int i=0; i<NO_PAIRS; i++)
			 {
				if( (FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] && FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j]) 
					 || (FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j] && FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j]) )
				{
				  TMP_CHARGES[j] = sqrt(fabs(FF_2BODY[i].PAIR_CHRG))*TMP_SIGN[j];
							
#if defined(USE_MPI) && defined(LINK_LAMMPS)
				  LMP_CHARGE[j] = TMP_CHARGES[j]; // save charges to global variable for LAMMPS
#endif
					
				  break;
				}
			 }
		  }
		  // Consistency check for charges.
		  for(int i=0; i<NO_PAIRS; i++)
		  {
			 for ( int j = 0 ; j < NATMTYP ; j++ ) 
			 {
				if ( FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] ) 
				{
				  for ( int k = 0 ; k < NATMTYP ; k++ ) 
				  {
					 if ( FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[k] ) 
					 {
						if ( fabs(TMP_CHARGES[j]*TMP_CHARGES[k] - FF_2BODY[i].PAIR_CHRG)  > eps_charge )
						{
						  EXIT_MSG("Error: Inconsistent pair charges for pair: " + FF_2BODY[i].ATM1TYP + FF_2BODY[i].ATM2TYP) ;
						}
					 }
				  }
				}
			 }
		  }
		  for(int a=0; a<SYSTEM.ATOMS;a++)
		  {
		  	for(int i=0; i<NATMTYP; i++)
		  	{
		  	  if(SYSTEM.ATOMTYPE[a] == TMP_ATOMTYPE[i])
		  	  {
		  		 SYSTEM.CHARGES[a] = TMP_CHARGES[i];
		  		 break;
		  	  }			
		  	}
		  }
		  check_charges(SYSTEM, TMP_CHARGES, TMP_ATOMTYPE, FF_2BODY, NATMTYP) ;
			
		  if(RANK==0)
		  {
			 cout << "Atom properties:" << endl ;
			 for(int j=0; j<NATMTYP; j++)
			 {
				cout << "		"<<	j << "     "
					  << setw(2) << left << TMP_ATOMTYPE[j] << ", q (e): " 
					  << setw(6) << fixed << setprecision(3) << right << TMP_CHARGES[j] << ", mass (amu): " 
					  << setw(8) << fixed << setprecision(4) << right << TMP_MASS[j] << endl;
			 }
		  }
		}
			
		if(RANK==0)			
		  cout << "	...Read 2-body FF params..." << endl;
	 }		
		
	 // Read 3B parameters
		
	 else if( (LINE.find(TEMP_SEARCH_3B) != string::npos) && CONTROLS.USE_3B_CHEBY) // "PAIR <PAIRTYPE> PARAMS"
	 {	
		TRIPS.read_ff_params(PARAMFILE, TMP_ATOMTYPE) ;
	 }
		
	 // Read 4B parameters
		
	 else if( (LINE.find(TEMP_SEARCH_4B) != string::npos) && CONTROLS.USE_4B_CHEBY) // "PAIR <PAIRTYPE> PARAMS"
	 {	
		QUADS.read_ff_params(PARAMFILE, TMP_ATOMTYPE) ;
	 }		
				
	 // Read the pair maps
		
	 else if(LINE.find("PAIRMAPS: ") != string::npos)
	 {
		STREAM_PARSER.str(LINE);			
		STREAM_PARSER >> TEMP_STR;
		STREAM_PARSER >> TMP_TERMS1;	
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();	
			
		if (RANK==0)
		  cout << "	Reading  " << TMP_TERMS1 << " pairs for mapping" << endl;
			
		for(int i=0; i<TMP_TERMS1; i++)
		{
		  PARAMFILE >> TMP_TERMS2;
		  PARAMFILE >> TEMP_TYPE;
				
		  if (RANK==0)
			 cout << "	........Reading pair: " << TEMP_TYPE << " with mapped index: " << TMP_TERMS2 << endl; 
				
		  PAIR_MAP.insert(make_pair(TEMP_TYPE,TMP_TERMS2));
		  PAIR_MAP_REVERSE.insert(make_pair(TMP_TERMS2,TEMP_TYPE));				
					
		}					
						
		if (RANK==0)
		  cout << "	...Read FF pairmaps..." << endl << endl;					
	 }	
		
	 // Read the triplet maps
		
	 else if(LINE.find("TRIPMAPS: ") != string::npos)
	 {
		TRIPS.read_maps(PARAMFILE, LINE) ;
	 }	
		
	 // Read the quadruplet maps

	 else if(LINE.find("QUADMAPS: ") != string::npos)
	 {
		QUADS.read_maps(PARAMFILE, LINE) ;
	 }
  }	
}

static void parse_ff_controls(string &LINE, ifstream &PARAMFILE, JOB_CONTROL &CONTROLS )
// Parse force field CONTROL flags.
{
  vector<string> tokens ;

  while ( parse_space(LINE, tokens) > 1 ) {
	 if ( tokens[0] == "USECOUL:" ) 
		CONTROLS.USE_COULOMB = is_true(tokens[1]) ;
	 else if ( tokens[0] == "FITCOUL:" ) 
		CONTROLS.FIT_COUL = is_true(tokens[1]) ;
	 else if ( tokens[0] == "USE3BCH:" )
		CONTROLS.USE_3B_CHEBY = is_true(tokens[1]) ;
	 else if ( tokens[0] == "USE4BCH:" )
		CONTROLS.USE_4B_CHEBY = is_true(tokens[1]) ;

	 getline(PARAMFILE, LINE) ;
	 if ( ! PARAMFILE.good() )
	 {
		cout << "Error reading job controls in paramfile\n" ;
		exit_run(0) ;
	 }
  }

  if(RANK==0)
  {
	 cout << "		...Compute electrostatics?      " << boolalpha << CONTROLS.USE_COULOMB << endl;
	 cout << "		...Use fit charges?             " << boolalpha << CONTROLS.FIT_COUL << endl;
	 cout << "		...Use 3-body Cheby params?     " << boolalpha << CONTROLS.USE_3B_CHEBY << endl;
	 cout << "		...Use 4-body Cheby params?     " << boolalpha << CONTROLS.USE_4B_CHEBY << endl;
			
	 cout << "	...Read FF controls..." << endl;	
  }

}

static void print_ff_summary(const vector<PAIR_FF> &FF_2BODY, CLUSTER_LIST& TRIPS,
									  CLUSTER_LIST& QUADS, const JOB_CONTROL &CONTROLS)
// Print out a summary of the force field.
{

  cout << "Force field summary: " << endl;
	
  cout << "	Pair type and number of parameters per pair: " << endl;
  cout << "		" << FF_2BODY[0].PAIRTYP << " " << FF_2BODY[0].SNUM << endl;
	
  if(FF_2BODY[0].PAIRTYP == "CHEBYSHEV")
  {
	 cout << "	Performing pair distance transformations over range: " << endl;
	 cout << "		" << FF_2BODY[0].CHEBY_RANGE_LOW << " to " << FF_2BODY[0].CHEBY_RANGE_HIGH << endl; 
  }
	
  cout  << endl;
	
  cout << "	Interaction parameters for each pair type: " << endl;

  for(int i=0; i<FF_2BODY.size(); i++)
  {
		
	 cout << "		pair type...smin...smax...sdelta...";
	 if(FF_2BODY[0].PAIRTYP == "CHEBYSHEV")
		cout << "cheby type";
	 if(FF_2BODY[i].CHEBY_TYPE == Cheby_trans::MORSE)
		cout << "...cheby lambda";
	 if(FF_2BODY[0].PAIRTYP == "CHEBYSHEV")
		cout << "...penalty dist...penalty scaling...cubic scaling...killlen";
	 cout << endl;
				
	 cout << "		" << FF_2BODY[i].PRPR_NM << " ";
		
	 cout << FF_2BODY[i].S_MINIM << " ";
	 cout << FF_2BODY[i].S_MAXIM << " ";
		
	 if(FF_2BODY[i].PAIRTYP == "CHEBYSHEV")
	 {
		string cheby_type = Cheby::get_trans_string(FF_2BODY[i].CHEBY_TYPE) ;
		cout << cheby_type << " ";
			
		if(FF_2BODY[i].CHEBY_TYPE == Cheby_trans::MORSE )
		  cout << FF_2BODY[i].LAMBDA << " ";
		cout << FF_2BODY[i].PENALTY_DIST << " " << scientific << FF_2BODY[i].PENALTY_SCALE<< " ";
		cout << FF_2BODY[i].KILLLEN;
	 }
		
	 cout << endl;
		
	 cout << "		Parameters: " << endl;
		
	 for(int j=0; j<FF_2BODY[i].SNUM; j++)
		cout << "		" << j << " " << FF_2BODY[i].PARAMS[j] << endl;
		
	 cout << endl; 
  }
	
  if(FF_2BODY[0].SNUM_3B_CHEBY > 0)
  {
	 TRIPS.print(true) ;
  }
		
  if(FF_2BODY[0].SNUM_4B_CHEBY > 0)
  {
	 QUADS.print(true) ;
  }

  if(CONTROLS.FIT_COUL)
  {
	 cout << "	Fitted charges read from parameter file:" << endl;
		
	 for(int i=0; i<FF_2BODY.size(); i++)
		cout << "		" << FF_2BODY[i].PRPR_NM << " " << FF_2BODY[i].PAIR_CHRG << " (" << FF_2BODY[i].PAIR_CHRG*ke << ")" << endl;
	 cout << endl;
  }	
  

  cout << endl;
}


static void read_coord_file(int index, JOB_CONTROL &CONTROLS, FRAME &SYSTEM, ifstream &CMPR_FORCEFILE)
// Read coordinates,  and optionally velocities and forces.  There is support for more than one coordinate
// input file, given by the index.
{
	ifstream 	COORDFILE;
	int 		TEMP_INT;
	BOX 		TMP_BOX;
	stringstream	STREAM_PARSER;
	string 		FIRST_EXT;
	string  	LINE;
	vector<string>	tokens;

	COORDFILE.open(CONTROLS.COORD_FILE[index].data());
	
	COORDFILE >> TEMP_INT;	// Store number of atoms in file in a temp var
	
	COORDFILE >> LINE;
	
	if(LINE == "NON_ORTHO")
	{
		SYSTEM.BOXDIM.IS_ORTHO = false;
		TMP_BOX.IS_ORTHO = false;
		COORDFILE >> TMP_BOX.CELL_AX;
	}
		
	else
		TMP_BOX.CELL_AX = stof(LINE);
	
	
	if (SYSTEM.BOXDIM.IS_ORTHO)
		COORDFILE >> TMP_BOX.CELL_BY >> TMP_BOX.CELL_CZ; // TMP_BOX.CELL_AX >> 
	else
	{
		COORDFILE >>                    TMP_BOX.CELL_AY >> TMP_BOX.CELL_AZ; // TMP_BOX.CELL_AX >> 
		COORDFILE >> TMP_BOX.CELL_BX >> TMP_BOX.CELL_BY >> TMP_BOX.CELL_BZ;
		COORDFILE >> TMP_BOX.CELL_CX >> TMP_BOX.CELL_CY >> TMP_BOX.CELL_CZ;
	}
	
	TMP_BOX.UPDATE_CELL();
	
	if(CONTROLS.FIT_STRESS || CONTROLS.FIT_STRESS_ALL)                                                                                           
		COORDFILE >>  SYSTEM.PRESSURE_TENSORS_ALL[0].X >>  SYSTEM.PRESSURE_TENSORS_ALL[1].Y >>  SYSTEM.PRESSURE_TENSORS_ALL[2].Z;    
	if (CONTROLS.FIT_STRESS_ALL) 
		COORDFILE >>  SYSTEM.PRESSURE_TENSORS_ALL[0].Y >>  SYSTEM.PRESSURE_TENSORS_ALL[0].Z >>  SYSTEM.PRESSURE_TENSORS_ALL[1].Z;

	if(CONTROLS.FIT_ENER)
		COORDFILE >> SYSTEM.QM_POT_ENER;   

	////////////////////////////////////////////////////////////
	// Read in the initial system coordinates, and if requested,
	// initial forces from separate file (i.e. not from .xyzf)
	// ... We also need to figure out the file extension so we
	// know how many fields to read on each atom lne
	////////////////////////////////////////////////////////////

	if (RANK==0)
		cout << "Reading initial coordinates and forces..." << endl;

	if ( CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE) 
	{
		if (RANK==0)
			cout << "      Opening " << CONTROLS.COMPARE_FILE.data() << " to read forces for comparison\n";
		CMPR_FORCEFILE.open(CONTROLS.COMPARE_FILE.data());

		if(!CMPR_FORCEFILE.is_open())
		{
			cout << "ERROR: Cannot open force input file: " << CONTROLS.COMPARE_FILE << endl;
			exit_run(0);
		}
	}
	
	string EXTENSION = CONTROLS.COORD_FILE[index].substr(CONTROLS.COORD_FILE[index].find_last_of(".")+1);

	if(index==0)
		FIRST_EXT = EXTENSION;
		
	if(EXTENSION != FIRST_EXT)
	{
		cout << "ERROR: Extensions for all input coordinate files must match. " << endl;
		cout << "	Found extension "<< FIRST_EXT << " for coordinate file 0 " << endl;
		cout << "	and " << EXTENSION << " for coordinate file " << index << endl;
		exit_run(0);
	}
	
	if (RANK==0)
	{
		cout << "     ...Read the following coordinate file extension: " << EXTENSION << endl;
		cout << "     ...Read the following number of atoms: " << TEMP_INT << endl;
		cout << "     ...Read box dimensions: " << endl;
		TMP_BOX.WRITE_BOX(CONTROLS.N_LAYERS);
		
		if(CONTROLS.FIT_STRESS || CONTROLS.FIT_STRESS_ALL)
			cout << "	...Read xx, yy, & zz stress tensors: " << SYSTEM.PRESSURE_TENSORS_ALL[0].X << " " << SYSTEM.PRESSURE_TENSORS_ALL[1].Y << " " << SYSTEM.PRESSURE_TENSORS_ALL[2].Z << endl;
		if(CONTROLS.FIT_STRESS_ALL)
			cout << "	...Read xy, xz, & yz stress tensors: " << SYSTEM.PRESSURE_TENSORS_ALL[0].Y << " " << SYSTEM.PRESSURE_TENSORS_ALL[0].Z << " " << SYSTEM.PRESSURE_TENSORS_ALL[1].Z << endl;
				
		if(CONTROLS.FIT_ENER)  													
			cout << "	...Read potential energy: " << SYSTEM.QM_POT_ENER << endl;					       
	}
	
	getline(COORDFILE,LINE);

	int     CURR_ATOM    = -1 ;
	XYZ_INT TMP_WRAP_IDX = {0,0,0};

	for(int a=0; a<TEMP_INT;a++)
	{		
		CURR_ATOM++;

		getline(COORDFILE,LINE);

		STREAM_PARSER.str(LINE);
		
		STREAM_PARSER >> SYSTEM.ATOMTYPE[CURR_ATOM];

		// Read the coordinates

		STREAM_PARSER >> SYSTEM.COORDS[CURR_ATOM].X >> SYSTEM.COORDS[CURR_ATOM].Y >> SYSTEM.COORDS[CURR_ATOM].Z;

		// Wrap the coordinates, shift along Z

		TMP_BOX.WRAP_ATOM(SYSTEM.COORDS[CURR_ATOM], TMP_WRAP_IDX, true);
		
		SYSTEM.COORDS[CURR_ATOM].Z += SYSTEM.BOXDIM.CELL_CZ;

		// Prepare velocities
		SYSTEM.VELOCITY[CURR_ATOM].X = 0;
		SYSTEM.VELOCITY[CURR_ATOM].Y = 0;
		SYSTEM.VELOCITY[CURR_ATOM].Z = 0;
		
		// Prepare forces
		SYSTEM.FORCES[CURR_ATOM].X = 0;
		SYSTEM.FORCES[CURR_ATOM].Y = 0;
		SYSTEM.FORCES[CURR_ATOM].Z = 0;
		
		// Prepare accelerations
		SYSTEM.ACCEL[CURR_ATOM].X = 0;
		SYSTEM.ACCEL[CURR_ATOM].Y = 0;
		SYSTEM.ACCEL[CURR_ATOM].Z = 0;		
		

		if ( CONTROLS.RESTART ) 
		{
			if(a==0 && RANK==0)
				cout << "	...Reading positions, velocities, and forces from a restart file " << endl ;

			STREAM_PARSER >> SYSTEM.VELOCITY[CURR_ATOM].X >> SYSTEM.VELOCITY[CURR_ATOM].Y  >> SYSTEM.VELOCITY[CURR_ATOM].Z;	
			STREAM_PARSER >> SYSTEM.ACCEL   [CURR_ATOM].X >> SYSTEM.ACCEL   [CURR_ATOM].Y  >> SYSTEM.ACCEL   [CURR_ATOM].Z;	
		} 
		else if ( CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE ) // Reading positions from *.xyzf for force testing
		{	
			if(EXTENSION == "xyzf")	// Ignore forces in .xyzf file
			{
				if(a==0 && RANK==0)
				{
				  cout << "	...Ignoring last three fields of atom info lines in xyzf file... " << endl;
				  cout << "	...will read from specified force file instead: " << CONTROLS.COMPARE_FILE << endl;
				}

				//COORDFILE >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;			
			}
					
			// Read forces from separate force file
			CMPR_FORCEFILE >> SYSTEM.FORCES[CURR_ATOM].X >> SYSTEM.FORCES[CURR_ATOM].Y >> SYSTEM.FORCES[CURR_ATOM].Z;
			
		}
		else if (!CONTROLS.INIT_VEL) // Reading positions from *.xyz
		{			
			if(EXTENSION == ".xyz")
			{
				cout << "ERROR: Input file requests velocities to be read in. " << endl;
				cout << "Expected .xyzf file, was given .xyzf file." << endl;
				exit_run(0);
			}
			
			// Read in velocities instead of forces...
			// Velocities must be stored so that the code can be restarted 
			// Maybe we should use a different file extension when velocities are stored. (LEF)		
			
			if(a==0 && RANK==0)
				cout << "	...Reading velocities from last three fields of atom info lines in xyzf file... " << endl;	

			STREAM_PARSER >> SYSTEM.VELOCITY[CURR_ATOM].X >> SYSTEM.VELOCITY[CURR_ATOM].Y >> SYSTEM.VELOCITY[CURR_ATOM].Z;
		}
		else if(EXTENSION == "xyzf")
		{
			if(a==0 && RANK==0)
				cout << "	...Ignoring last three fields of atom info lines in xyzf file... " << endl;
			//COORDFILE >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;
		}
		
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();
	}

	// Read in stress components for comparison to current values.
	if ( CONTROLS.FIT_STRESS ) 
	{
		// Units read in are kcal/mol-Ang.^3.
		CMPR_FORCEFILE >> SYSTEM.STRESS_TENSORS.X >> SYSTEM.STRESS_TENSORS.Y >> SYSTEM.STRESS_TENSORS.Z ;

		if ( ! CMPR_FORCEFILE.good() ) 
			EXIT_MSG("Reading force comparson file failed") ;
	} else if ( CONTROLS.FIT_STRESS_ALL ) {
				// Units read in are kcal/mol-Ang.^3.
		CMPR_FORCEFILE >> SYSTEM.STRESS_TENSORS_X.X >> SYSTEM.STRESS_TENSORS_X.Y >> SYSTEM.STRESS_TENSORS_X.Z ;
		CMPR_FORCEFILE >> SYSTEM.STRESS_TENSORS_Y.X >> SYSTEM.STRESS_TENSORS_Y.Y >> SYSTEM.STRESS_TENSORS_Y.Z ;
		CMPR_FORCEFILE >> SYSTEM.STRESS_TENSORS_Z.X >> SYSTEM.STRESS_TENSORS_Z.Y >> SYSTEM.STRESS_TENSORS_Z.Z ;

		if ( ! CMPR_FORCEFILE.good() ) 
			EXIT_MSG("Reading force comparson file failed") ;
	}
	if ( CONTROLS.FIT_ENER )
	{
		vector<double> energy(3) ;
		CMPR_FORCEFILE >> energy[0] ;
		CMPR_FORCEFILE >> energy[1] ;
		CMPR_FORCEFILE >> energy[2] ;

		if ( ! CMPR_FORCEFILE.good() ) 
			EXIT_MSG("Reading force comparson file failed") ;
		else if ( fabs(energy[0]-energy[1]) > 1.0e-08 || fabs(energy[0]-energy[2]) > 1.0e-03 )
			EXIT_MSG("Inconsistent energy values in force comparison file") ;

		SYSTEM.QM_POT_ENER = energy[0] ;
	}
	// Input configurations are added sequentially along z.
		
	SYSTEM.BOXDIM.CELL_AX  = TMP_BOX.CELL_AX;
	SYSTEM.BOXDIM.CELL_AY  = TMP_BOX.CELL_AY;
	SYSTEM.BOXDIM.CELL_AZ  = TMP_BOX.CELL_AZ;
	
	SYSTEM.BOXDIM.CELL_BX  = TMP_BOX.CELL_BX;
	SYSTEM.BOXDIM.CELL_BY  = TMP_BOX.CELL_BY;
	SYSTEM.BOXDIM.CELL_BZ  = TMP_BOX.CELL_BZ;

	SYSTEM.BOXDIM.CELL_CX  = TMP_BOX.CELL_CX;
	SYSTEM.BOXDIM.CELL_CY  = TMP_BOX.CELL_CY;				
	SYSTEM.BOXDIM.CELL_CZ += TMP_BOX.CELL_CZ;
		
	if(CONTROLS.SCALE_SYSTEM_BY != 1.0)
	{
		if (! SYSTEM.BOXDIM.IS_ORTHO)
			EXIT_MSG("ERROR: Box scaling for non-orthorhombic cells has not yet been implemented.");
		
		for(int a1=0; a1<SYSTEM.ATOMS; a1++)
		{
			SYSTEM.COORDS[a1].X *= CONTROLS.SCALE_SYSTEM_BY;
			SYSTEM.COORDS[a1].Y *= CONTROLS.SCALE_SYSTEM_BY;
			SYSTEM.COORDS[a1].Z *= CONTROLS.SCALE_SYSTEM_BY;
		}
			
		SYSTEM.BOXDIM.CELL_AX *= CONTROLS.SCALE_SYSTEM_BY;
		SYSTEM.BOXDIM.CELL_BY *= CONTROLS.SCALE_SYSTEM_BY;
		SYSTEM.BOXDIM.CELL_CZ *= CONTROLS.SCALE_SYSTEM_BY;
	}
	
	SYSTEM.BOXDIM.UPDATE_CELL();	
		
	if (RANK==0)
	{
		cout <<  "     ...Updated simulation box dimensions: " << endl;
		SYSTEM.BOXDIM.WRITE_BOX(CONTROLS.N_LAYERS);
	}

	if ( ! CONTROLS.RESTART ) 
	{
		COORDFILE.close();
		COORDFILE.clear();
	}
	
	if ( CONTROLS.COMPARE_FORCE ) 
		CMPR_FORCEFILE.close();
	
	if(RANK==0)
		cout << "   ...read complete for file " << CONTROLS.COORD_FILE[index] << endl << endl;	
}

static void subtract_force(FRAME &SYSTEM, JOB_CONTROL &CONTROLS)
{
	// To print the modified frame, used if CONTROLS.SUBTRACT_FORCE true
			
	ofstream FORCE_SUBTRACTED_OUTPUT;
	string   FORCE_SUBTRACTED_FILE   = CONTROLS.COORD_FILE[0];
			
	int END = SYSTEM.ATOMS;

	if ( RANK != 0 ) 
		return ;

	if(CONTROLS.SUBTRACT_FORCE)
	{
		FORCE_SUBTRACTED_FILE  .append("_forces_subtracted.xyz");
		FORCE_SUBTRACTED_OUTPUT.open(FORCE_SUBTRACTED_FILE.data());
				
		FORCE_SUBTRACTED_OUTPUT << END << endl;
		
		if (SYSTEM.BOXDIM.IS_ORTHO)
		{
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_AX << " " << SYSTEM.BOXDIM.CELL_BY << " " << SYSTEM.BOXDIM.CELL_CZ << " ";
		}
		else
		{
			FORCE_SUBTRACTED_OUTPUT << "NON_ORTHO ";
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_AX << " ";
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_AY << " ";
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_AZ << " ";	    
		
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_BX << " ";
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_BY << " ";
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_BZ << " ";
		
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_CX << " ";
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_CY << " ";
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.CELL_CZ << " ";	    
		}
				
		// NOTE:  PRESSURE_TENSORS_XYZ hold the potential energy contribution to the current stress tensor.
		// PRESSURE_TENSOR holds potential + kinetic energy to the current stress tensor in GPa units.
		// PRESSURE_TENSORS_XYZ should be used here.
		
		if(CONTROLS.FIT_STRESS)
			FORCE_SUBTRACTED_OUTPUT << 	  SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X - SYSTEM.STRESS_TENSORS.X 
						<< " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y - SYSTEM.STRESS_TENSORS.Y 
						<< " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z - SYSTEM.STRESS_TENSORS.Z;
		else if ( CONTROLS.FIT_STRESS_ALL )						
			FORCE_SUBTRACTED_OUTPUT  << 	   SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X - SYSTEM.STRESS_TENSORS_X.X 
						 << " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Y - SYSTEM.STRESS_TENSORS_X.Y 
						 << " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Z - SYSTEM.STRESS_TENSORS_X.Z
						 << " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X - SYSTEM.STRESS_TENSORS_Y.X 
						 << " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y - SYSTEM.STRESS_TENSORS_Y.Y 
						 << " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Z - SYSTEM.STRESS_TENSORS_Y.Z
						 << " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].X - SYSTEM.STRESS_TENSORS_Z.X 
						 << " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y - SYSTEM.STRESS_TENSORS_Z.Y 
						 << " " << SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z - SYSTEM.STRESS_TENSORS_Z.Z;
						
						
		if(CONTROLS.FIT_ENER)
		{
			//FORCE_SUBTRACTED_OUTPUT << SYSTEM.TOT_POT_ENER - SYSTEM.QM_POT_ENER + SYSTEM.QM_ENERGY_OFFSET ;
		
			// MAJOR ASSUMPTION: atom type orders are identical to the LSQ code
		
			double EDIFF = SYSTEM.TOT_POT_ENER - SYSTEM.QM_POT_ENER;
			
			SYSTEM.SET_NATOMS_OF_TYPE();
			
			for(int i=0; i<SYSTEM.QM_ENERGY_OFFSET.size(); i++)
				EDIFF += SYSTEM.QM_ENERGY_OFFSET[i]*SYSTEM.NATOMS_OF_TYPE[i];
		
			FORCE_SUBTRACTED_OUTPUT << EDIFF;
		}
					
		FORCE_SUBTRACTED_OUTPUT << endl;
	}

	// Check against read-in forces for code verification... Note, ferr is initialized to zero.
	double ferr = 0.0;
		
	for(int a1=0;a1<END;a1++)
	{
		ferr += (SYSTEM.ACCEL[a1].X - SYSTEM.FORCES[a1].X) * (SYSTEM.ACCEL[a1].X - SYSTEM.FORCES[a1].X);
		ferr += (SYSTEM.ACCEL[a1].Y - SYSTEM.FORCES[a1].Y) * (SYSTEM.ACCEL[a1].Y - SYSTEM.FORCES[a1].Y);
		ferr += (SYSTEM.ACCEL[a1].Z - SYSTEM.FORCES[a1].Z) * (SYSTEM.ACCEL[a1].Z - SYSTEM.FORCES[a1].Z);
				
		// Before printing force file with current ff forces subtracted, convert from simulation units (kca/mol/Ang)
		// to Hartree/bohr

		if ( CONTROLS.SUBTRACT_FORCE ) 
		{
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.ATOMTYPE[a1] << "	"
						<< SYSTEM.COORDS[a1].X << "   " << SYSTEM.COORDS[a1].Y << "   " << SYSTEM.COORDS[a1].Z << "   "						
						<< (SYSTEM.FORCES[a1].X - SYSTEM.ACCEL[a1].X)/(Hartree*Bohr) << " "
						<< (SYSTEM.FORCES[a1].Y - SYSTEM.ACCEL[a1].Y)/(Hartree*Bohr) << " "
						<< (SYSTEM.FORCES[a1].Z - SYSTEM.ACCEL[a1].Z)/(Hartree*Bohr) << endl; 			     
		}
	}
			
	if(CONTROLS.SUBTRACT_FORCE)
		FORCE_SUBTRACTED_OUTPUT.close();
			
	ferr = sqrt(ferr/3.0/END);
	cout << "RMS force error = " << fixed << setprecision(6) << ferr << endl;

	ferr = 0.0 ;
	if ( CONTROLS.FIT_STRESS )
	{
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X - SYSTEM.STRESS_TENSORS.X) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X - SYSTEM.STRESS_TENSORS.X) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y - SYSTEM.STRESS_TENSORS.Y) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y - SYSTEM.STRESS_TENSORS.Y) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z - SYSTEM.STRESS_TENSORS.Z) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z - SYSTEM.STRESS_TENSORS.Z) ;

		ferr = sqrt(ferr/3.0) ;
		ferr *= GPa ;
		cout << "RMS stress error = " << ferr << " GPa " << endl ;
	}
	else if ( CONTROLS.FIT_STRESS_ALL )
	{
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X - SYSTEM.STRESS_TENSORS_X.X) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X - SYSTEM.STRESS_TENSORS_X.X) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Y - SYSTEM.STRESS_TENSORS_X.Y) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Y - SYSTEM.STRESS_TENSORS_X.Y) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Z - SYSTEM.STRESS_TENSORS_X.Z) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Z - SYSTEM.STRESS_TENSORS_X.Z) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X - SYSTEM.STRESS_TENSORS_Y.X) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X - SYSTEM.STRESS_TENSORS_Y.X) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y - SYSTEM.STRESS_TENSORS_Y.Y) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y - SYSTEM.STRESS_TENSORS_Y.Y) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Z - SYSTEM.STRESS_TENSORS_Y.Z) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Z - SYSTEM.STRESS_TENSORS_Y.Z) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].X - SYSTEM.STRESS_TENSORS_Z.X) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].X - SYSTEM.STRESS_TENSORS_Z.X) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y - SYSTEM.STRESS_TENSORS_Z.Y) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y - SYSTEM.STRESS_TENSORS_Z.Y) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z - SYSTEM.STRESS_TENSORS_Z.Z) * (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z - SYSTEM.STRESS_TENSORS_Z.Z) ;

		ferr = sqrt(ferr/3.0) ;
		ferr *= GPa ;
		cout << "RMS stress error = " << ferr << " GPa " << endl ;
	}

	if ( CONTROLS.FIT_ENER )
	{
		//ferr = fabs(SYSTEM.TOT_POT_ENER - SYSTEM.QM_POT_ENER + SYSTEM.QM_ENERGY_OFFSET)/SYSTEM.ATOMS ;

		// MAJOR ASSUMPTION: atom type orders are identical to the LSQ code
		
		ferr = SYSTEM.TOT_POT_ENER - SYSTEM.QM_POT_ENER;
		
		SYSTEM.SET_NATOMS_OF_TYPE();
		
		for(int i=0; i<SYSTEM.QM_ENERGY_OFFSET.size(); i++)
			ferr += SYSTEM.QM_ENERGY_OFFSET[i]*SYSTEM.NATOMS_OF_TYPE[i];
		
		ferr = fabs(ferr/SYSTEM.ATOMS);
		
		cout << "Absolute energy error = " << ferr << " kcal/mol/atom \n" ;
	}
}
