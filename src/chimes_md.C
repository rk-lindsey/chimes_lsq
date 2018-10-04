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

// For LAMMPS linking

#if defined(USE_MPI) && defined(LINK_LAMMPS)
	#include "lmppath.h"
	using namespace LAMMPS_NS;
#endif
	
using namespace std;	
	
// Define function headers -- general

static void read_input        (JOB_CONTROL & CONTROLS, PES_PLOTS & FF_PLOTS, NEIGHBORS & NEIGHBOR_LIST);
static void read_atom_types(ifstream &PARAMFILE, JOB_CONTROL &CONTROLS, int &NATMTYP, vector<string>& TMP_ATOMTYPE, vector<int>& TMP_NATOMTYPE, vector<int>& TMP_ATOMTYPEIDX, vector<double>& TMP_CHARGES,  vector<double>& TMP_MASS, vector<int> &TMP_SIGN) ;
static void read_ff_params(ifstream &PARAMFILE, JOB_CONTROL &CONTROLS, vector<PAIR_FF>& FF_2BODY, CLUSTER_LIST& TRIPS, CLUSTER_LIST &QUADS, map<string,int> &PAIR_MAP, NEIGHBORS &NEIGHBOR_LIST, FRAME& SYSTEM, int N_PLOTS, int NATMTYP, const vector<string>& TMP_ATOMTYPE, const vector<int>& TMP_ATOMTYPEIDX, vector<double> &TMP_CHARGES, vector<double> &TMP_MASS, const vector<int>& TMP_SIGN, map<int,string>& PAIR_MAP_REVERSE) ;
static void print_ff_summary(const vector<PAIR_FF> &FF_2BODY, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, const JOB_CONTROL &CONTROLS) ;
static void print_pes(PES_PLOTS &FF_PLOTS, vector<PAIRS> &FF_2BODY, map<string,int> &PAIR_MAP, map<int,string> &PAIR_MAP_REVERSE, CLUSTER_LIST& TRIPS, CLUSTER_LIST& QUADS, JOB_CONTROL &CONTROLS, Cheby &cheby) ;

// Define function headers -- MPI

static void write_xyzv         (FRAME &SYSTEM, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL, THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBOR_LIST, string filename, bool restart);
static void read_restart_params(ifstream &COORDFILE, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL, THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBORS, FRAME &SYSTEM) ;
static void parse_ff_controls  (string &LINE, ifstream &PARAMFILE, JOB_CONTROL &CONTROLS ) ;
static void read_coord_file(int index, JOB_CONTROL &CONTROLS, FRAME &SYSTEM, ifstream &CMPR_FORCEFILE) ;
static void subtract_force(FRAME &SYSTEM, JOB_CONTROL &CONTROLS) ;

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

// Variables that are defined locally for house_md which need to be global for LAMMPS linking
 
#if defined(USE_MPI) && defined(LINK_LAMMPS)

	// Gloabl variables needed by callback function for LAMMPS  
	
	int TOTAL_ATOMS;
	vector<double>	LMP_CHARGE;		// Global data object to hold charges for the LAMMPS input files     
	vector<string>  TMP_ATOMTYPE;		// Will be used by lammps to map lammps atom types (int's) back to chemical symbols
	JOB_CONTROL	CONTROLS;		// Declare the data object that will hold the main simulation control variables
	map<string,int> PAIR_MAP;		// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
  	vector<int> INT_PAIR_MAP ;		
	map<int,string> PAIR_MAP_REVERSE; 	// Input/output the resverse of PAIR_MAP
	map<string,int>& TRIAD_MAP;		// Input is three of any ATOM_PAIRS.PRPR_NM pairs, output is a triplet type index
	map<int,string>& TRIAD_MAP_REVERSE;	// Input/output is the reverse of TRIAD_MAP
	map<string,int>& QUAD_MAP;		// Input is 6 of any ATOM_PAIRS.PRPR_NM pairs, output is a quadruplet type index
	map<int,string>& QUAD_MAP_REVERSE;	// Input/output is the reverse of QUAD_MAP
	
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
	
  PES_PLOTS FF_PLOTS;			// Declare the data object that will help set up PES plots
  FRAME      SYSTEM;			// Declare the data object that will hold the system coordinates, velocities, accelrations, etc.
  CONSTRAINT ENSEMBLE_CONTROL;		// Declare the class that will handle integration/constraints
  THERMO_AVG AVG_DATA ;
	
  // These are data objects that are normally defined locally, but need to be global for LAMMPS linking
	
#if not defined(LINK_LAMMPS)
	
  JOB_CONTROL CONTROLS;			// Declare the data object that will hold the main simulation control variables
  NEIGHBORS   NEIGHBOR_LIST;		// Declare the class that will handle the neighbor list
	
  // Data objects to hold coefficients for different force field types, and for FF printing (if requested)

  vector<PAIR_FF> FF_2BODY;		// Holds all 2-body parameters

  CLUSTER_LIST TRIPS;			// Holds all 3-body parameters
  CLUSTER_LIST QUADS;      		// Holds all 4-body parameters
		
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

  read_input(CONTROLS, FF_PLOTS,NEIGHBOR_LIST);		// Populate object with user defined values
	
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
  // Hop to PES printing, if requested 
  ////////////////////////////////////////////////////////////
	
  if (FF_PLOTS.N_PLOTS > 0)
	 goto FF_SETUP_1; 

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
	
  SYSTEM.BOXDIM.X = SYSTEM.BOXDIM.Y = SYSTEM.BOXDIM.Z = 0;
	
  if(CONTROLS.BUILD)
  {
	 SYSTEM.BOXDIM.X = SYSTEM.BOXDIM.Y = SYSTEM.BOXDIM.Z = CONTROLS.BUILD_BOXL;
		
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
	 int GRIDSPACE  = ceil(SYSTEM.BOXDIM.X/GRIDPOINTS);
		
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
		 read_coord_file(i, CONTROLS, SYSTEM, CMPR_FORCEFILE) ;
	 }
	}
  ////////////////////////////////////////////////////////////
  // Figure out atom charges and masses, based on parameter 
  // file
  ////////////////////////////////////////////////////////////
	 
FF_SETUP_1: 
	
  if(RANK==0) 
	 cout << "Reading atom info from parameter file..." << endl; 
	
  // Read in the possible atom types and their features	

  read_atom_types(PARAMFILE, CONTROLS, NATMTYP, TMP_ATOMTYPE, TMP_NATOMTYPE, TMP_ATOMTYPEIDX, TMP_CHARGES, TMP_MASS, TMP_SIGN) ;
	
  if (FF_PLOTS.N_PLOTS > 0)
	 goto FF_SETUP_2; 
	

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
	 build_real_replicates(SYSTEM, CONTROLS) ;
	
  // New layer handling: ghost atoms.. but this should be called whether or not ghost atoms are requested
  // (Required for neighbor lists, which are now ALWAYS used)


  SYSTEM.WRAPDIM.X = SYSTEM.BOXDIM.X * (2*CONTROLS.N_LAYERS + 1);
  SYSTEM.WRAPDIM.Y = SYSTEM.BOXDIM.Y * (2*CONTROLS.N_LAYERS + 1);
  SYSTEM.WRAPDIM.Z = SYSTEM.BOXDIM.Z * (2*CONTROLS.N_LAYERS + 1);

  build_layers(SYSTEM, CONTROLS);
	
  if ( (CONTROLS.N_LAYERS > 0) && (RANK == 0) )	// Then ghost atoms are used 
  {
	 cout << "	Real atoms:                   " << SYSTEM.ATOMS     << endl;
	 cout << "	Total atoms (ghost):          " << SYSTEM.ALL_ATOMS << endl;
	 cout << "	Real box dimesntions:         " << SYSTEM.BOXDIM.X  << " " << SYSTEM.BOXDIM.Y  << " " << SYSTEM.BOXDIM.Z  << endl;
	 cout << "	Total box dimensions (ghost): " << SYSTEM.WRAPDIM.X << " " << SYSTEM.WRAPDIM.Y << " " << SYSTEM.WRAPDIM.Z << endl << endl;
			
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

  if ( CONTROLS.RESTART ) 
	 read_restart_params(COORDFILE, CONTROLS, ENSEMBLE_CONTROL, AVG_DATA, NEIGHBOR_LIST, SYSTEM) ;

  Vol = SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;

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

FF_SETUP_2:

  read_ff_params(PARAMFILE, CONTROLS, FF_2BODY, TRIPS, QUADS, PAIR_MAP, NEIGHBOR_LIST, SYSTEM, FF_PLOTS.N_PLOTS, NATMTYP, TMP_ATOMTYPE, TMP_ATOMTYPEIDX,TMP_CHARGES, TMP_MASS, TMP_SIGN, PAIR_MAP_REVERSE) ;
	
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

  if(NEIGHBOR_LIST.USE && !CONTROLS.PLOT_PES)
  {
	 if(RANK == 0)
		cout << "Initializing the neighbor list..." << endl;
		
	 NEIGHBOR_LIST.INITIALIZE_MD(SYSTEM);
	 NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);
  }
	
	
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


  ////////////////////////////////////////////////////////////
  // If PES printing is requested, do so here and then 
  // exit the program
  ////////////////////////////////////////////////////////////  	
	
  Cheby cheby{CONTROLS, SYSTEM, NEIGHBOR_LIST, FF_2BODY, INT_PAIR_MAP} ;

  if(RANK==0)
  {
	 if (FF_PLOTS.N_PLOTS > 0)
	 {
		print_pes(FF_PLOTS, FF_2BODY, PAIR_MAP, PAIR_MAP_REVERSE, TRIPS, QUADS, CONTROLS, cheby) ;
		exit_run(0);
	 }
  }	

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
	
  for(CONTROLS.STEP=FIRST_STEP;CONTROLS.STEP<CONTROLS.N_MD_STEPS;CONTROLS.STEP++)	//start Big Loop here.
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
		 	{
			 	SYSTEM.COORDS[a1].X -= floor(SYSTEM.COORDS[a1].X / SYSTEM.BOXDIM.X) * SYSTEM.BOXDIM.X;
			 	SYSTEM.COORDS[a1].Y -= floor(SYSTEM.COORDS[a1].Y / SYSTEM.BOXDIM.Y) * SYSTEM.BOXDIM.Y;
			 	SYSTEM.COORDS[a1].Z -= floor(SYSTEM.COORDS[a1].Z / SYSTEM.BOXDIM.Z) * SYSTEM.BOXDIM.Z;
		 	}	
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
	 sum_forces(SYSTEM.ACCEL, SYSTEM.ATOMS, SYSTEM.TOT_POT_ENER, SYSTEM.PRESSURE_XYZ, SYSTEM.PRESSURE_TENSORS_XYZ.X, SYSTEM.PRESSURE_TENSORS_XYZ.Y, SYSTEM.PRESSURE_TENSORS_XYZ.Z);
#endif
		
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
			
		SYSTEM.PRESSURE_TENSORS.X = (SYSTEM.PRESSURE_TENSORS_XYZ.X + 2.0 / 3.0 * Ktot / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS.Y = (SYSTEM.PRESSURE_TENSORS_XYZ.Y + 2.0 / 3.0 * Ktot / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS.Z = (SYSTEM.PRESSURE_TENSORS_XYZ.Z + 2.0 / 3.0 * Ktot / Vol) * GPa;
				
		AVG_DATA.STRESS_TENSOR_SUM.X += SYSTEM.PRESSURE_TENSORS.X;
		AVG_DATA.STRESS_TENSOR_SUM.Y += SYSTEM.PRESSURE_TENSORS.Y;
		AVG_DATA.STRESS_TENSOR_SUM.Z += SYSTEM.PRESSURE_TENSORS.Z;
		
		////////////////////////////////////////////////////////////
		// Periodically print simulation output
		////////////////////////////////////////////////////////////
		
		if ( (CONTROLS.STEP+1) % CONTROLS.FREQ_ENER == 0 && RANK == 0 ) 
		{
		 	printf("%8d %9.2f %15.7f %15.7f %15.7f %15.1f %15.8f", CONTROLS.STEP+1, (CONTROLS.STEP+1)*CONTROLS.DELTA_T_FS, Ktot/SYSTEM.ATOMS,SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS,(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS,SYSTEM.TEMPERATURE, SYSTEM.PRESSURE);
			STATISTICS << CONTROLS.STEP+1<< "     " << (CONTROLS.STEP+1)*CONTROLS.DELTA_T_FS<< "  " << Ktot/SYSTEM.ATOMS<< "      " <<SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS<< "        " <<(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS<< " " <<SYSTEM.TEMPERATURE<< "      " << SYSTEM.PRESSURE;
			
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
		double Pavg = (AVG_DATA.STRESS_TENSOR_SUM.X + AVG_DATA.STRESS_TENSOR_SUM.Y + AVG_DATA.STRESS_TENSOR_SUM.Z)/3.0/ CONTROLS.N_MD_STEPS ;
		cout << "	Pressures from diagonal stress tensors over run: " << Pavg << endl;
		cout << "	Average diagonal stress tensors over run: " << endl;
		cout << "		sigma_xx: " << AVG_DATA.STRESS_TENSOR_SUM.X/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_yy: " << AVG_DATA.STRESS_TENSOR_SUM.Y/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_zz: " << AVG_DATA.STRESS_TENSOR_SUM.Z/CONTROLS.N_MD_STEPS << " GPa" << endl;     
	 }

	 // Write the final configuration to file.
	 write_xyzv(SYSTEM, CONTROLS, ENSEMBLE_CONTROL, AVG_DATA, NEIGHBOR_LIST, "output.xyz", false);
			
	 STATISTICS.close();
  }

  // MPI -- End our setup
	
	normal_exit() ;
}       


static void read_input(JOB_CONTROL & CONTROLS, PES_PLOTS & FF_PLOTS, NEIGHBORS & NEIGHBOR_LIST) 				// UPDATED
{
	if (RANK==0)
		cout << endl << "Reading the simulation control input file..." << endl;
	
	// Define some variable to help with reading
	
	bool   		FOUND_END = false;
	string 		LINE;
	string		TEMP_STR,TEMP_STR_2;
	stringstream	LINE_PARSER;
	int		ADD_TO_NPLOTS = 0;
	vector<string> tokens ;
	int ntokens ;

	// Set some defaults
	
	CONTROLS.IS_LSQ           = false;
	CONTROLS.SELF_CONSIST     = false;
	CONTROLS.SUBTRACT_FORCE   = false;
	CONTROLS.PLOT_PES         = false;
	CONTROLS.WRAP_COORDS      = true;
	CONTROLS.PRINT_VELOC      = false;
	CONTROLS.NVT_CONV_CUT     = 0.10;
	CONTROLS.FREEZE_IDX_START = -1; 
	CONTROLS.FREEZE_IDX_STOP  = -1; 
	CONTROLS.SCALE_SYSTEM_BY  = 1.0;
	CONTROLS.BUILD            = false;
	CONTROLS.FIT_STRESS       = false;
	CONTROLS.FIT_ENER         = false;
	
	CONTROLS.RESTART          = false;
	
	FF_PLOTS.INCLUDE_FCUT     = true;
	FF_PLOTS.INCLUDE_CHARGES  = true;
	FF_PLOTS.INCLUDE_PENALTY  = true;
	FF_PLOTS.DO_4D            = false;
	
	CONTROLS.FREQ_BACKUP     	= 100;
	CONTROLS.FREQ_UPDATE_THERMOSTAT = -1.0;
	CONTROLS.USE_HOOVER_THRMOSTAT	= false;
	CONTROLS.REAL_REPLICATES        = 0;
	NEIGHBOR_LIST.USE               = true;
	
	CONTROLS.PRINT_BAD_CFGS         = false;
	CONTROLS.TRAJ_FORMAT		= "GEN";
	
	CONTROLS.PENALTY_THRESH = -1.0;	// Default is to not enforce a max-allowed penalty
	CONTROLS.IO_ECONS_VAL   =  0.0;
	
	
	
	
	// Begin reading
	
	while (FOUND_END == false)
	{
		getline(cin,LINE);
		if ( cin.good() ) 
		{
			// Remove comments.
			strip_comments(LINE) ;

			ntokens = parse_space(LINE, tokens) ;
			if ( ntokens == 0 ) continue ;

		}
		else
			EXIT_MSG("Input file terminated without an # ENDFILE # command.") ;

		// Break out of the loop
		if     (LINE.find("# ENDFILE #") != string::npos)
		{
			FOUND_END = true;
			
			if( (CONTROLS.N_LAYERS>0) && (!NEIGHBOR_LIST.USE) )
			{
				if (RANK == 0)
				{
					if (isatty(fileno(stdout)))
						cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "WARNING: Use of neighbor lists HIGHLY reccommended when NLAYERS > 0!" << COUT_STYLE.ENDSTYLE << endl;

					else
						cout << "WARNING: Use of neighbor lists HIGHLY reccommended when NLAYERS > 0!" << endl;
					
					//#if WARN == TRUE
					//	exit_run(0);
					//#endif
				}
			}
			
			if (RANK==0)
				cout << "   ...read complete." << endl << endl;
			
			break;
		}

		// Variables for printing out PES for a given parameter file
		
		else if(LINE.find("# PLOTPES #") != string::npos)
		{
			if (RANK==0)
				cout << "	# PLOTPES #: true... will only plot PES's" << endl;	
			
			getline(cin,LINE);
			LINE_PARSER.str(LINE);
			LINE_PARSER >> LINE;
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				CONTROLS.PLOT_PES = true;
				
				// Read in the number of PES plots to produce

				LINE_PARSER >> FF_PLOTS.N_PLOTS;
				
				// Read in the parameter file
				
				LINE_PARSER >> CONTROLS.PARAM_FILE;		
				
				while(LINE_PARSER >> TEMP_STR)
				{
					if(TEMP_STR=="Exclude"  || TEMP_STR=="exclude"  || TEMP_STR=="EXCLUDE")
					{
						if(LINE_PARSER >> TEMP_STR)
						{
							if(TEMP_STR=="Fcut"  || TEMP_STR=="fcut"  || TEMP_STR=="FCUT")
							{
								FF_PLOTS.INCLUDE_FCUT = false;

								if (RANK==0)
								{
									if (isatty(fileno(stdout)))
										cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "		Excluding cubic scaling (fcut)" << COUT_STYLE.ENDSTYLE << endl;	
									else
										cout << "		Excluding cubic scaling (fcut)" << endl;
								}
							}
							else if (TEMP_STR=="Charges"  || TEMP_STR=="charges"  || TEMP_STR=="CHARGES")
							{
								FF_PLOTS.INCLUDE_CHARGES = false;
								if (RANK==0)
								{
									if (isatty(fileno(stdout)))
										cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "		Excluding charges" << COUT_STYLE.ENDSTYLE << endl;	
									else
										cout << "		Excluding charges" << endl;
								}
							}
							else if (TEMP_STR=="Penalty"  || TEMP_STR=="penalty"  || TEMP_STR=="PENALTY")
							{
								FF_PLOTS.INCLUDE_PENALTY = false;
								if (RANK==0)
								{
									if (isatty(fileno(stdout)))
										cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "		Excluding penalty for 2b (only) scan" << COUT_STYLE.ENDSTYLE << endl;	
									else
										cout << "		Excluding penalty for 2b (only) scan" << endl;
								}
							}
						}
					}
					
					if(TEMP_STR=="No"  || TEMP_STR=="no"  || TEMP_STR=="NO")
					{
						if(LINE_PARSER >> TEMP_STR)
						{
							if(TEMP_STR=="4D"  || TEMP_STR=="4d")
							{
								FF_PLOTS.DO_4D = false;
								if (RANK==0)
								{
									if (isatty(fileno(stdout)))
										cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "		Will not print the 4D 3-body data" << COUT_STYLE.ENDSTYLE << endl;	
									else
										cout << "		Will not print the 4D 3-body data" << endl;
								}
							}
						}
					}
				}
				if(FF_PLOTS.INCLUDE_FCUT)
				{
					if (isatty(fileno(stdout)))
						cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "		Including cubic scaling (fcut)" << COUT_STYLE.ENDSTYLE << endl;	
					else
						cout << "		Including cubic scaling (fcut)" << endl;						
				}
				if(FF_PLOTS.INCLUDE_CHARGES)
				{
					if (isatty(fileno(stdout)))
						cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "		Including charges" << COUT_STYLE.ENDSTYLE << endl;	
					else
						cout << "		Including charges" << endl;						
				}
				if(FF_PLOTS.INCLUDE_PENALTY)
				{
					if (isatty(fileno(stdout)))
						cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "		Including penalty for 2b (only) scan" << COUT_STYLE.ENDSTYLE << endl;	
					else
						cout << "		Including penalty for 2b (only) scan" << endl;						
				}
				LINE_PARSER.str("");
				LINE_PARSER.clear();

				// Read in the search string in a way that makes spacing not matter

				for(int i=0; i<FF_PLOTS.N_PLOTS; i++)
				{
					TEMP_STR = "";
					
					getline(cin,LINE);
					
					LINE_PARSER.str("");
					LINE_PARSER.clear();					
					LINE_PARSER.str(LINE);	
					
					if (RANK==0)
						cout << endl << "	Processing line: " << LINE << endl;				
					
					for(int j=0; j<3; j++)
					{
						LINE_PARSER >> LINE;

						if(j==0)
						{
							if(LINE=="PAIRTYPE")
								FF_PLOTS.NBODY.push_back(2);
							else if(LINE=="TRIPLETTYPE")
								FF_PLOTS.NBODY.push_back(3);
							else
							{
								cout << "ERROR: Unrecognized interaction type: " << LINE << endl;
								cout << "       Allowed values are PAIRTYPE and TRIPLETTYPE" << endl;
								exit_run(0);
							}
						}
						 
						TEMP_STR.append(LINE);
						TEMP_STR.append(" ");							

					}

					FF_PLOTS.TYPE_INDEX.push_back(atoi(LINE.c_str()));
					
					FF_PLOTS.PES_TYPES.push_back(TEMP_STR);

				}

				// Exit.. don't need anything else.
				
				FOUND_END = true;
				
				if (RANK==0)
					cout << "   ...read complete." << endl << endl;
				
				FF_PLOTS.N_PLOTS += ADD_TO_NPLOTS;
				break;				
			}
			else
				if (RANK==0)
					cout << "	# PLOTPES #: false" << endl;	
		}
		
		// For self-consistent fitting type runs
		
		else if(LINE.find("# SLFCNST #") != string::npos)
		{
			cout << "NOTE: Feature \"SLFCNST\" is no longer supported. Ignoring." << endl;
			cin >> LINE; 
			cin.ignore();
		}

		// "General control variables"
		
		else if(LINE.find("# RNDSEED #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.SEED = int(atof(LINE.data()));
			if (RANK==0)
				cout << "	# RNDSEED #: " << CONTROLS.SEED << endl;	
		}	
		
		else if(LINE.find("# TEMPERA #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.TEMPERATURE = double(atof(LINE.data()));
			if (RANK==0)
				cout << "	# TEMPERA #: " << CONTROLS.TEMPERATURE << " K" << endl;	
		}	
		
		else if(LINE.find("# PRESSUR #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.PRESSURE = double(atof(LINE.data()));
			if (RANK==0)
				cout << "	# PRESSUR #: " << CONTROLS.PRESSURE << " GPa" << endl;	
		}
		
		else if(LINE.find("# CONVCUT #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.NVT_CONV_CUT = double(atof(LINE.data()));
			if (RANK==0)
				cout << "	# CONVCUT #: " << CONTROLS.NVT_CONV_CUT*100 << " % of set T" << endl;	
		}	

		else if(LINE.find("# CMPRFRC #") != string::npos)
		{
			getline(cin,LINE);
			LINE_PARSER.str(LINE);
			LINE_PARSER >> LINE;
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				CONTROLS.COMPARE_FORCE = true;
				
				// Read in the name of the file

				LINE_PARSER >> TEMP_STR;
				LINE_PARSER.str("");
				LINE_PARSER.clear();
				
				// Make sure user actually entered a name
				
				if(LINE == TEMP_STR)
				{
					cout << "ERROR: If # CMPRFRC # is true, a force file must be specified on the same line." << endl; 
					cout << "       Example: true comparison_forces.txt" << endl; 
					exit_run(0);
				}
				else
					CONTROLS.COMPARE_FILE = TEMP_STR;
			}
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				CONTROLS.COMPARE_FORCE = false;
			else
			{
				cout << "ERROR: # CMPRFRC # must be specified as true or false." << endl;
				exit_run(1);	
			}					
				
			if(CONTROLS.COMPARE_FORCE)
			{
				if (RANK==0)
					cout << "	# CMPRFRC #: true... will only do 1 md step." << endl;	
				
				// If CONTROLS.COMPARE_FORCE is true, then only do one step,
				// because all we're trying to do is to see if md computed
				// forces match the expected forces

				CONTROLS.N_MD_STEPS = 1;
			}
			else
				if (RANK==0)
					cout << "	# CMPRFRC #: false" << endl;	
		}

		else if(LINE.find("# CHECKFRC #") != string::npos)
		{
			getline(cin,LINE);
			vector<string> tokens ;

			parse_space(LINE, tokens) ;

			if ( is_true(tokens[0]) ) 
			{
			  CONTROLS.CHECK_FORCE = true;
			  if ( RANK == 0 )
				 cout << "\t# CHECKFRC #: true" << endl ;
			}
			else 
			{
			  CONTROLS.CHECK_FORCE = false ;
			  if ( RANK == 0 )
				 cout << "\t# CHECKFRC #: false" << endl ;
			}

		}

		else if(LINE.find("# SUBTFRC #") != string::npos)
		{
			getline(cin,LINE);
			
			LINE_PARSER.str("");
			LINE_PARSER.clear();
			
			LINE_PARSER.str(LINE);
			LINE_PARSER >> LINE;
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				CONTROLS.SUBTRACT_FORCE = true;
				
				// Read in the name of the file

				LINE_PARSER >> TEMP_STR;
				
				// Should we also be subtracting off stresses and/or energy?
				
				CONTROLS.FIT_ENER = CONTROLS.FIT_STRESS = false;
				
				if( bool(LINE_PARSER >> LINE))
				{
					if(LINE == "SUBTR_ENERGY")
						CONTROLS.FIT_ENER = true;
					else if(LINE == "SUBTR_STRESS")
						CONTROLS.FIT_STRESS = true;
					else
					{
						cout << "ERROR: Unrecognized SUBTFRC option " << LINE << ". Allowed values are SUBTR_ENERGY and SUBTR_STRESS.";
						exit(0);
					}
				}
				
				if( bool(LINE_PARSER >> LINE))
				{
					if(LINE == "SUBTR_ENERGY")
						CONTROLS.FIT_ENER = true;
					else if(LINE == "SUBTR_STRESS")
						CONTROLS.FIT_STRESS = true;
					else
					{
						cout << "ERROR: Unrecognized SUBTFRC option " << LINE << ". Allowed values are SUBTR_ENERGY and SUBTR_STRESS.";
						exit(0);
					}
				}
				
				LINE_PARSER.str("");
				LINE_PARSER.clear();
				
				// Make sure user actually entered a name
				
				if(LINE == TEMP_STR)
				{
					cout << "ERROR: If # SUBTFRC # is true, a force file must be specified on the same line." << endl; 
					cout << "       Example: true dft_forces.txt" << endl; 
					exit_run(0);
				}
				else
					CONTROLS.COMPARE_FILE = TEMP_STR;
			}
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				CONTROLS.SUBTRACT_FORCE = false;
			else
			{
				cout << "ERROR: # SUBTFRC # must be specified as true or false." << endl;
				exit_run(1);	
			}					
				
			if(CONTROLS.SUBTRACT_FORCE)
			{
				if (RANK==0)
				{
					cout << "	# SUBTFRC #: true... will only do 1 md step." << endl;	
					
					if(CONTROLS.FIT_ENER)
						cout << "		...Will subtract energies." << endl;	
					if(CONTROLS.FIT_STRESS)
						cout << "		...Will subtract stress tensors." << endl;	
				}
				// If CONTROLS.COMPARE_FORCE is true, then only do one step,
				// because all we're trying to do is to see if md computed
				// forces match the expected forces

				CONTROLS.N_MD_STEPS = 1;
			}
			else
				if (RANK==0)
					cout << "	# SUBTFRC #: false" << endl;	
		}

		else if(LINE.find("# TIMESTP #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.DELTA_T_FS = double(atof(LINE.data()));
			if (RANK==0)
				cout << "	# TIMESTP #: " << CONTROLS.DELTA_T_FS << " fs" << endl;
			
			// Now convert timestep to simulation units
			CONTROLS.DELTA_T = CONTROLS.DELTA_T_FS/Tfs;
		}	
		
		else if(LINE.find("# N_MDSTP #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.N_MD_STEPS = int(atof(LINE.data()));
			if (RANK==0)
				cout << "	# N_MDSTP #: " << CONTROLS.N_MD_STEPS << endl;	
		}
		
		else if(LINE.find("# PENTHRS #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.PENALTY_THRESH = double(atof(LINE.data()));
			if (RANK==0)
				cout << "	# PENTHRS #: " << CONTROLS.PENALTY_THRESH << endl;	
		}		

		else if(LINE.find("# NLAYERS #") != string::npos)
		{
			LINE_PARSER.str("");
			LINE_PARSER.clear();
			
			getline(cin,LINE);
			
			LINE_PARSER.str(LINE);
			LINE_PARSER >> LINE;
				
			CONTROLS.N_LAYERS = int(atof(LINE.data()));
			if (RANK==0)
				cout << "	# NLAYERS #: " << CONTROLS.N_LAYERS << endl;	
			
			if(LINE_PARSER >> LINE)
			{
				if(LINE == "REPLICATE")
				{
					LINE_PARSER >> LINE;
					CONTROLS.REAL_REPLICATES = int(atof(LINE.data()));
					if (RANK==0)
						cout << "	             ... Creating " << CONTROLS.REAL_REPLICATES << " real replicates before ghost atom layering." << endl;	
				}
			}
		}	
		
		else if(LINE.find("# USENEIG #") != string::npos)
		{
			cin >> LINE; 
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				NEIGHBOR_LIST.USE = true;

				if (RANK==0 && NEIGHBOR_LIST.USE )
					cout << "	# USENEIG #: true "<< endl;	
			}
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f") 
			{
				// Everybody is considered a neighbor.
				NEIGHBOR_LIST.USE = false ;
				
				if (RANK==0 && NEIGHBOR_LIST.USE )
					cout << "	# USENEIG #: false "<< endl;	
			}
			cin.ignore();
		}	
		
		else if(LINE.find("# PRMFILE #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.PARAM_FILE = LINE;
			if (RANK==0)
				cout << "	# PRMFILE #: " << CONTROLS.PARAM_FILE << endl;	
		}	
		
		else if(LINE.find("# CRDFILE #") != string::npos)
		{
			LINE_PARSER.str("");
			LINE_PARSER.clear();
			
			getline(cin,LINE);
			
			LINE_PARSER.str(LINE);
			LINE_PARSER >> LINE;
			
			if(LINE == "CAT") // Then we'll be catenating several files together
			{
				cout << "	# CRDFILE #: Creating a cell from multiple input files: ";
				int NO_FILES;
				LINE_PARSER >> NO_FILES;
				cout << NO_FILES  << endl;
				
				for(int i=0; i<NO_FILES; i++)
				{
					LINE_PARSER >> LINE;
					CONTROLS.COORD_FILE.push_back(LINE);
					cout << "		" << CONTROLS.COORD_FILE[i] << endl;
				}	
				
				cout << "	Note: Assumes files have the same x and y box dimensions!" << endl;			
			}
			if(LINE == "SCALE") // Then coordinates/boxlengths will be scaled
			{
				LINE_PARSER >> CONTROLS.SCALE_SYSTEM_BY;
				
				LINE_PARSER >> LINE;
				CONTROLS.COORD_FILE.push_back(LINE);
				
				if (RANK==0)
					cout << "	# CRDFILE #: " << CONTROLS.COORD_FILE[0] << endl;
			}
			if(LINE == "INITIALIZE")
			{
				CONTROLS.BUILD = true;
				
				LINE_PARSER >> CONTROLS.BUILD_TYPE;
				
				if (CONTROLS.BUILD_TYPE == "MOLECULAR")
					LINE_PARSER >>  CONTROLS.BUILD_FILE;
					
				else if (CONTROLS.BUILD_TYPE == "ATOMIC")
					LINE_PARSER >>  CONTROLS.BUILD_ATOM;
				
				LINE_PARSER >> LINE;
				
				if(LINE != "BOXL")
					cout << "ERROR: Expected \'BOXL <value> NMOLEC <value>." << endl;
					
				LINE_PARSER >> CONTROLS.BUILD_BOXL;
				
				LINE_PARSER >> LINE;
				
				if(LINE != "NMOLEC")
					cout << "ERROR: Expected \'BOXL <value> NMOLEC <value>." << endl;

				LINE_PARSER >> CONTROLS.BUILD_NMOLEC;	
				
			}
			else
			{
				CONTROLS.COORD_FILE.push_back(LINE);
				
				if (RANK==0)
					cout << "	# CRDFILE #: " << CONTROLS.COORD_FILE[0] << endl;	
			}
			
			LINE_PARSER.str("");
			LINE_PARSER.clear();
		}			
		
		// " Simulation options"
		
		else if(LINE.find("# VELINIT #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			CONTROLS.INIT_VEL = false;
			
			if (LINE=="READ")
				CONTROLS.INIT_VEL = false;
			else if (LINE=="GEN")
				CONTROLS.INIT_VEL = true;
			else if ( LINE == "RESTART" ) 
				CONTROLS.RESTART = true ;
			else
			{
				cout << "ERROR: # VELINIT # must be specified as READ or GEN." << endl;
				exit_run(1);	
			}			
			
			if (RANK==0)
			{
				if(!CONTROLS.INIT_VEL && CONTROLS.BUILD)
				{
					cout << "ERROR: # VELINIT # must be \'GEN\' when # CRDFILE # \'INITIALIZE\' option is used." << endl;
					exit_run(0);
				}
			}
			

			if (RANK==0)
			{
				if(CONTROLS.INIT_VEL)
					cout << "	# VELINIT #: GEN ... generating velocites via box Muller" << endl;	
				else if(CONTROLS.RESTART)
					cout << "	# VELINIT #: RESTART ... reading velocities from coordinate (restart) file" << endl;	
				else if(!CONTROLS.INIT_VEL && !CONTROLS.RESTART)	
					cout << "	# VELINIT #: READ ... reading velocities from coordinate file" << endl;	
			}
		}	
				
		else if(LINE.find("# CONSRNT #") != string::npos)
		{
			
			LINE_PARSER.str("");
			LINE_PARSER.clear();
			
			getline(cin,LINE);
			
			LINE_PARSER.str(LINE);
			LINE_PARSER >> CONTROLS.ENSEMBLE;	// Which ensemble do we have?
				
			
			if( CONTROLS.ENSEMBLE != "NVT-SCALE" && CONTROLS.ENSEMBLE != "NVT-MTK"    && CONTROLS.ENSEMBLE != "NVE" && 
				CONTROLS.ENSEMBLE != "NPT-MTK"   && CONTROLS.ENSEMBLE != "NPT-BEREND" && CONTROLS.ENSEMBLE != "NVT-BEREND" &&
				CONTROLS.ENSEMBLE != "NPT-BEREND-ANISO" && 
				CONTROLS.ENSEMBLE != "LMP-NVE"	 && CONTROLS.ENSEMBLE != "LMP-NVT"	  && CONTROLS.ENSEMBLE != "LMP-NPT" &&
				CONTROLS.ENSEMBLE != "LMP-MIN-BOX-ISO" &&CONTROLS.ENSEMBLE != "LMP-MIN-BOX-ANISO" && CONTROLS.ENSEMBLE != "LMP-MIN-BOX-TRI" && CONTROLS.ENSEMBLE != "LMP-MIN")
			{
				cout << "ERROR: Unrecognized ensemble: " << CONTROLS.ENSEMBLE << endl;
				cout << "       Options are: NVE, NVT-SCALE, NVT-MTK, NPT-BEREND, or NPT-MTK for plain house_md, " << endl;
				cout << "       and LMP-NVE, LMP-NVT, LMP-NPT, LMP-MIN-BOX-ISO,  LMP-MIN-BOX-ANISO, LMP-MIN-BOX-TRI, and MP-MIN for md linked to LAMMPS" << endl;
				exit_run(0);
			}
			
			
			
			
			
			if(CONTROLS.ENSEMBLE == "LMP-NVE" || CONTROLS.ENSEMBLE == "LMP-NVT" || CONTROLS.ENSEMBLE == "LMP-NPT")
			{
				if(CONTROLS.ENSEMBLE == "LMP-NVT" || CONTROLS.ENSEMBLE == "LMP-NPT")
					LINE_PARSER >> CONTROLS.FREQ_UPDATE_THERMOSTAT;
	
				if(CONTROLS.ENSEMBLE == "LMP-NPT")
					LINE_PARSER >> CONTROLS.FREQ_UPDATE_BAROSTAT;

			}
			else if(CONTROLS.ENSEMBLE == "LMP-MIN-BOX-ISO" || CONTROLS.ENSEMBLE == "LMP-MIN-BOX-ANISO" || CONTROLS.ENSEMBLE == "LMP-MIN-BOX-TRI" || CONTROLS.ENSEMBLE == "LMP-MIN")
			{
				LINE_PARSER >> CONTROLS.MIN_E_CONVG_CRIT;
				LINE_PARSER >> CONTROLS.MIN_F_CONVG_CRIT;
				LINE_PARSER >> CONTROLS.MIN_MAX_ITER;
				LINE_PARSER >> CONTROLS.MIN_MAX_EVAL;
			}
			else
			{
				if(CONTROLS.ENSEMBLE != "NVE")
				{
					// Keep reading..
				
					LINE_PARSER >> LINE;
				
					if (CONTROLS.ENSEMBLE == "NVT-MTK" || CONTROLS.ENSEMBLE == "NPT-MTK")
					{
						CONTROLS.USE_HOOVER_THRMOSTAT = true;
						LINE_PARSER >> LINE;
						CONTROLS.FREQ_UPDATE_THERMOSTAT = double(atof(LINE.data()));
						if (RANK==0)
							cout << "	# CONSRNT #: HOOVER... a Hoover time of " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " will be used." << endl;

						if(CONTROLS.ENSEMBLE == "NPT-MTK")
						{
							if( bool(LINE_PARSER >> LINE))
							{
								CONTROLS.FREQ_UPDATE_BAROSTAT = double(atof(LINE.data()));
							}
							else
							{
								CONTROLS.FREQ_UPDATE_BAROSTAT = 1000;		
							}
							if (RANK==0)
								cout << "	                   ... and barostat will use  " << CONTROLS.FREQ_UPDATE_BAROSTAT << "." << endl;	
						}
	
					}
					else if (CONTROLS.ENSEMBLE == "NVT-SCALE" || CONTROLS.ENSEMBLE == "NPT-BEREND" || CONTROLS.ENSEMBLE == "NVT-BEREND" || CONTROLS.ENSEMBLE == "NPT-BEREND-ANISO")
					{
						CONTROLS.USE_HOOVER_THRMOSTAT = false;
						LINE_PARSER >> LINE;
				
						CONTROLS.FREQ_UPDATE_THERMOSTAT = double(atof(LINE.data()));

						if (CONTROLS.ENSEMBLE == "NPT-BEREND" || CONTROLS.ENSEMBLE == "NPT-BEREND-ANISO")
						{
							if( bool(LINE_PARSER >> LINE))
							{
								LINE_PARSER >> LINE; // ignores the "COORDSCALE" command.
								CONTROLS.FREQ_UPDATE_BAROSTAT = double(atof(LINE.data()));
							}
							else
							{
								CONTROLS.FREQ_UPDATE_BAROSTAT = 1000;		
							}
						}
				
						if (RANK==0)
							cout << "	# CONSRNT #: " << CONTROLS.ENSEMBLE << "... Velocities will be scaled every " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " MD steps." << endl;	
					}
					else
					{
						cout << "ERROR: # CONSRNT # must be specified as HOOVER or VELSCALE, and a Hoover time or velocity " << endl;
						cout << "       scaling frequency must be specified in line. " << endl;
						cout << "       Example: HOOVER 10 "<< endl; 
						exit_run(1);	
					}	
				}
				

			
				if(LINE_PARSER >> LINE)
				{
					if(LINE == "FREEZE")
					{
						cout << "		...Atoms over the following range will be frozen will be frozen (indexed from 0): "; 
					
						LINE_PARSER >> LINE;
						CONTROLS.FREEZE_IDX_START = int(atof(LINE.data()));
						cout << CONTROLS.FREEZE_IDX_START << " - ";
					
						LINE_PARSER >> LINE;
						CONTROLS.FREEZE_IDX_STOP = int(atof(LINE.data()));
						cout << CONTROLS.FREEZE_IDX_STOP << endl;
					
						LINE_PARSER.str("");
						LINE_PARSER.clear();
					}
				}
			}
		}			
		
		else if(LINE.find("# PRSCALC #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if (LINE=="ANALYTICAL")
			{
				CONTROLS.USE_NUMERICAL_PRESS = false;
				if (RANK==0)
					cout << "	# PRSCALC #: ANALYTICAL" << endl;	
			}
			else if (LINE=="NUMERICAL")
			{
				CONTROLS.USE_NUMERICAL_PRESS = true;
				if (RANK==0)
					cout << "	# PRSCALC #: NUMERICAL" << endl;	
			}
			else
			{
				cout << "ERROR: # PRSCALC # must be specified as ANALYTICAL or NUMERICAL." << endl;
				exit_run(1);	
			}				
		}	
				
		// "Output control"
		
		else if(LINE.find("# WRPCRDS #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				CONTROLS.WRAP_COORDS = true;
				if (RANK==0)
					cout << "	# WRPCRDS #: true" << endl;	
			}
			else
			{
				CONTROLS.WRAP_COORDS = false;
				if (RANK==0)
					cout << "	# WRPCRDS #: false" << endl;							
			}
		}			
		
		else if( (LINE.find("# FRQDFTB #") != string::npos) || (LINE.find("# FRQTRAJ #") != string::npos) )
		{
			cin >> LINE; cin.ignore();
			CONTROLS.FREQ_DFTB_GEN = int(atof(LINE.data()));
			if (RANK==0)
				cout << "	# FRQTRAJ #: " << CONTROLS.FREQ_DFTB_GEN << endl;	
			
			if (CONTROLS.DELTA_T_FS > 0 && CONTROLS.N_MD_STEPS > 0 && RANK==0)
			{
				cout << "		... printing every " << CONTROLS.FREQ_DFTB_GEN*CONTROLS.DELTA_T_FS << " fs, " << endl;
				cout << "		... printing       " << CONTROLS.N_MD_STEPS/CONTROLS.FREQ_DFTB_GEN << " frames. " << endl; 
			}
		}
		
		
		else if( LINE.find("# TRAJEXT #") != string::npos)
		{
			cin >> CONTROLS.TRAJ_FORMAT; cin.ignore();

			if (RANK==0)
				cout << "	# TRAJEXT #: " << CONTROLS.TRAJ_FORMAT << endl;	
		}		

		else if(LINE.find("# FRQENER #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.FREQ_ENER = int(atof(LINE.data()));
			if (RANK==0)
				cout << "	# FRQENER #: " << CONTROLS.FREQ_ENER << endl;	
		}		

		else if(LINE.find("# PRNTFRC #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				CONTROLS.PRINT_FORCE = true;
				cin >> LINE; cin.ignore();
				
				if(LINE == "FRQDFTB")
					CONTROLS.FREQ_FORCE = CONTROLS.FREQ_DFTB_GEN;
				else
					CONTROLS.FREQ_FORCE = int(atof(LINE.data()));
				
				if (RANK==0)
					cout << "	# PRNTFRC #: true ...and will be printed every " << CONTROLS.FREQ_FORCE << " frames."<< endl;	
				
				if (CONTROLS.DELTA_T_FS > 0 && CONTROLS.N_MD_STEPS > 0 && RANK == 0)
				{
					cout << "		... printing every " << CONTROLS.FREQ_FORCE*CONTROLS.DELTA_T_FS << " fs, " << endl;
					cout << "		... printing " << CONTROLS.N_MD_STEPS/CONTROLS.FREQ_FORCE << " frames. " << endl; 
				}
			}
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
			{
				CONTROLS.PRINT_FORCE = false;
				if (RANK==0)
					cout << "	# PRNTFRC #: false" << endl;	
			}
			else
			{
				cout << "ERROR: # PRNTFRC # must be specified as true or false." << endl;
				exit_run(1);	
			}								
		}	
		
		else if(LINE.find("# PRNTBAD #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				CONTROLS.PRINT_BAD_CFGS = true;
				//cin >> LINE; cin.ignore();
				
				if (RANK==0)
					cout << "	# PRNTBAD #: true ...will print out all bad configs (r < r_cut and r_cut < r < [r_cut + d_penalty])" << endl;	

			}
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
			{
				CONTROLS.PRINT_BAD_CFGS = false;
				if (RANK==0)
					cout << "	# PRNTBAD #: false" << endl;	
			}
			cin.ignore();							
		}		
		
		else if(LINE.find("# PRNTVEL #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				CONTROLS.PRINT_VELOC = true;
				cin >> LINE; cin.ignore();
				
				if(LINE == "FRQDFTB")
					CONTROLS.FREQ_VELOC = CONTROLS.FREQ_DFTB_GEN;
				else
					CONTROLS.FREQ_VELOC = int(atof(LINE.data()));
				
				if (RANK==0)
					cout << "	# PRNTVEL #: true ...and will be printed every " << CONTROLS.FREQ_VELOC << " frames."<< endl;	
				
				if (CONTROLS.DELTA_T_FS > 0 && CONTROLS.N_MD_STEPS > 0 && RANK == 0)
				{
					cout << "		... printing every " << CONTROLS.FREQ_VELOC*CONTROLS.DELTA_T_FS << " fs, " << endl;
					cout << "		... printing " << CONTROLS.N_MD_STEPS/CONTROLS.FREQ_VELOC << " frames. " << endl; 
				}
				
			}
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
			{
				CONTROLS.FREQ_VELOC = false;
				if (RANK==0)
					cout << "	# PRNTVEL #: false" << endl;	
			}
			else
			{
				cout << "ERROR: # PRNTVEL # must be specified as true or false." << endl;
				exit_run(1);	
			}								
		}			
		else if(LINE.find("# FITSTRS #") != string::npos)
		{
			cin >> LINE;
			
			if (LINE=="first"  || LINE=="First"  || LINE=="FIRST")
			{
				CONTROLS.FIT_STRESS = true;
				cin >> CONTROLS.NSTRESS;
				cin.ignore();
			}			
			else
			{
				if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
					CONTROLS.FIT_STRESS = true;	
				else if (LINE=="all"  || LINE=="All"  || LINE=="ALL"  || LINE == "A" || LINE == "a")
					CONTROLS.FIT_STRESS_ALL = true;
				else if(LINE=="firstall"  || LINE=="FirstAll"  || LINE=="FIRSTALL")
				{
					CONTROLS.FIT_STRESS_ALL = true;
					cin >> CONTROLS.NSTRESS;
					cin.ignore();
				}	
				else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
					CONTROLS.FIT_STRESS = false;
				else
				{
					cout << endl << "ERROR: # FITSTRS # must be specified as true or false." << endl;
					exit(1);	
				}
				
				cin.ignore();
			}
#if VERBOSITY == 1
			
			if ( RANK == 0 ) 
			{
				cout << "	# FITSTRS #: ";
			
				if (CONTROLS.FIT_STRESS)
					cout << "true" << endl;							
				else if (CONTROLS.FIT_STRESS_ALL)
					cout << "true ...will fit to all tensor components" << endl;	
				else
					cout << "false" << endl;
					
				if(CONTROLS.NSTRESS>0)
					cout << "    			 ...will only fit tensors for first " << CONTROLS.NSTRESS << " frames." << endl;
			}
				
#endif
		}
		
		else if(LINE.find("# FITENER #") != string::npos)
		{
			cin >> LINE; 
			
			if (LINE=="first"  || LINE=="First"  || LINE=="FIRST")
			{
				CONTROLS.FIT_ENER = true;
				cin >> CONTROLS.NENER;
				cin.ignore();
			}
			else
			{

				if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
					CONTROLS.FIT_ENER = true;
				else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
					CONTROLS.FIT_ENER = false;
				else
				{
					cout << endl << "ERROR: # FITENER # must be specified as true or false." << endl;
					exit(1);	
				}
				
				cin.ignore();	
			}
			
#if VERBOSITY == 1
			if ( RANK == 0 )
			{
				cout << "	# FITENER #: ";
			
				if (CONTROLS.FIT_ENER)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
					
				if(CONTROLS.NENER>0)
					cout << "    			 ...will only fit energies for first " << CONTROLS.NENER << " frames." << endl;
								
				
			}
#endif
		}
		else if ( RANK == 0 ) 
		{
			cout << "WARNING:  The following input line was ignored:\n" ;
			cout << " " << LINE << endl ;
		}
	}
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
	fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.X << " ";
	fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.Y << " ";
	fxyz << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.Z << endl;
	
	for ( int ia = 0; ia < SYSTEM.ATOMS; ia++ ) 
	{
		XYZ tmp = SYSTEM.COORDS[ia] ;
		
		if ( CONTROLS.WRAP_COORDS ) 	// Wrap into the primitive cell
		{
			
			tmp.X -= floor(tmp.X/SYSTEM.BOXDIM.X)*SYSTEM.BOXDIM.X;
			tmp.Y -= floor(tmp.Y/SYSTEM.BOXDIM.Y)*SYSTEM.BOXDIM.Y;
			tmp.Z -= floor(tmp.Z/SYSTEM.BOXDIM.Z)*SYSTEM.BOXDIM.Z;
		}

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

static void read_restart_params(ifstream &COORDFILE, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL,
										  THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBOR_LIST, FRAME &SYSTEM)
{
	COORDFILE >> CONTROLS.STEP ;
	COORDFILE >> NEIGHBOR_LIST.RCUT_PADDING ;
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
	
		LMPFILE << "0 " << SYSTEM.BOXDIM.X << " xlo xhi" << endl;
		LMPFILE << "0 " << SYSTEM.BOXDIM.Y << " ylo yhi" << endl;
		LMPFILE << "0 " << SYSTEM.BOXDIM.Z << " zlo zhi" << endl << endl;  
	
		// LAMMPS allow for a triclinic simulation cell; skew parameters included here.
		LMPFILE << "0.00000000 0.00000000 0.00000000 xy xz yz" << endl;  
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
		
		SYS.BOXDIM.X   = xprd;
		SYS.BOXDIM.Y   = yprd;
		SYS.BOXDIM.Z   = zprd;
		
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
			
  	   	   	SYS.COORDS[i].X -= floor(SYS.COORDS[i].X/SYS.BOXDIM.X)*SYS.BOXDIM.X;
			SYS.COORDS[i].Y -= floor(SYS.COORDS[i].Y/SYS.BOXDIM.Y)*SYS.BOXDIM.Y;
			SYS.COORDS[i].Z -= floor(SYS.COORDS[i].Z/SYS.BOXDIM.Z)*SYS.BOXDIM.Z;
			
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
		
		build_layers(SYS, CONTROLS) ;

		SYS.WRAPDIM.X = SYS.BOXDIM.X * (2*CONTROLS.N_LAYERS + 1);
		SYS.WRAPDIM.Y = SYS.BOXDIM.Y * (2*CONTROLS.N_LAYERS + 1);
		SYS.WRAPDIM.Z = SYS.BOXDIM.Z * (2*CONTROLS.N_LAYERS + 1);
		
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
	    
		double vol = SYS.BOXDIM.X * SYS.BOXDIM.Y * SYS.BOXDIM.Z;

		virial[0] =  SYS.PRESSURE_TENSORS_XYZ.X*vol;
		virial[1] =  SYS.PRESSURE_TENSORS_XYZ.Y*vol;
		virial[2] =  SYS.PRESSURE_TENSORS_XYZ.Z*vol;
	    
		virial[3] = 0.0; // XY
		virial[4] = 0.0; // XZ
		virial[5] = 0.0; // YZ
    
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

static void read_ff_params(ifstream &PARAMFILE, JOB_CONTROL &CONTROLS, vector<PAIR_FF>& FF_2BODY,
									CLUSTER_LIST& TRIPS, CLUSTER_LIST &QUADS, map<string,int> &PAIR_MAP, 
									NEIGHBORS &NEIGHBOR_LIST, FRAME& SYSTEM, int N_PLOTS, int NATMTYP,
                           const vector<string>& TMP_ATOMTYPE, const vector<int>& TMP_ATOMTYPEIDX,
									vector<double> &TMP_CHARGES, vector<double> &TMP_MASS,
									const vector<int>& TMP_SIGN, map<int,string>& PAIR_MAP_REVERSE)
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
				NEIGHBOR_LIST.MAX_CUTOFF_3B = TRIPS.MAX_CUTOFF ;

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
				NEIGHBOR_LIST.MAX_CUTOFF_4B = QUADS.MAX_CUTOFF ;
		  }
		}

		if(RANK==0)
		{
		  cout << "   ...read complete." << endl << endl;
		  cout << "Notes on simulation: " << endl;
		  cout << "	Using fpenalty power " << FPENALTY_POWER << endl;
		  
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
			
		if(CONTROLS.USE_OVERCOORD)
		{
		  if(RANK==0)
		  {
			 cout << "	Overbonding contributions will be considered." << endl;
				
			 if(CONTROLS.FIT_POVER)
				cout << "		p-over will be read from fit parameters." << endl;
			 else
				cout << "		p-over will be read from specified parameters." << endl;
		  }
		}
			
		if(RANK==0)
		  cout << endl;
			
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
			
		if(TEMP_TYPE == "DFTBPOLY" || TEMP_TYPE == "INVRSE_R")
		{
		  STREAM_PARSER >> TMP_TERMS1;
		}
		else if(TEMP_TYPE =="CHEBYSHEV")
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
			 FF_2BODY[i].CUBIC_SCALE      = 1.0;
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
		  for ( int j = 0 ; j < TMP_ATOMTYPE.size() ; j++ ) {
			 if ( FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] ) 
				FF_2BODY[i].ATM1TYPE_IDX = TMP_ATOMTYPEIDX[j] ;
			 if ( FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j] ) 
				FF_2BODY[i].ATM2TYPE_IDX = TMP_ATOMTYPEIDX[j] ;
		  }
		  if ( FF_2BODY[i].ATM1TYPE_IDX < 0 )
			 EXIT_MSG("Could not find a match for " + FF_2BODY[i].ATM1TYP ) ;

		  if ( FF_2BODY[i].ATM2TYPE_IDX < 0 )
			 EXIT_MSG("Could not find a match for " + FF_2BODY[i].ATM2TYP ) ;
		  
		  if(FF_2BODY[i].S_MAXIM > NEIGHBOR_LIST.MAX_CUTOFF)
		  {
			 NEIGHBOR_LIST.MAX_CUTOFF    = FF_2BODY[i].S_MAXIM;
			 NEIGHBOR_LIST.MAX_CUTOFF_3B = FF_2BODY[i].S_MAXIM;
			 NEIGHBOR_LIST.MAX_CUTOFF_4B = FF_2BODY[i].S_MAXIM;
		  }
				 	
		  PARAMFILE >> FF_2BODY[i].S_DELTA;
				
		  if((N_PLOTS == 0) &&
			  (  FF_2BODY[i].S_MAXIM > 0.5* SYSTEM.WRAPDIM.X
				  || FF_2BODY[i].S_MAXIM > 0.5* SYSTEM.WRAPDIM.Y
				  || FF_2BODY[i].S_MAXIM > 0.5* SYSTEM.WRAPDIM.Z ) )
		  {
#if WARN == TRUE
			 if(RANK==0)
			 {
				if (isatty(fileno(stdout)))
				{
				  cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "WARNING: Outer cutoff greater than half at least one box length" << COUT_STYLE.ENDSTYLE << endl;
				  cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Pair type " <<FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;
				  cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD <<  SYSTEM.BOXDIM.X << " " << SYSTEM.BOXDIM.Y << " " << SYSTEM.BOXDIM.Z << COUT_STYLE.ENDSTYLE << endl;
				}
				else
				{
				  cout << "WARNING: Outer cutoff greater than half at least one box length" << endl;
				  cout << "	Pair type " <<FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM2TYP  << endl;
				  cout << " " <<  SYSTEM.BOXDIM.X << " " << SYSTEM.BOXDIM.Y << " " << SYSTEM.BOXDIM.Z << endl;
				}								
			 }

#else
			 if (isatty(fileno(stdout)))
			 {
				cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "ERROR: Outer cutoff greater than half at least one box length" << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Pair type " <<FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;
				exit_run(0);
			 }
			 else
			 {
				cout << "ERROR: Outer cutoff greater than half at least one box length" << endl;
				cout << "	Pair type " <<FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM2TYP << endl;
				exit_run(0);							
			 }
#endif
		  }

				
		  FF_2BODY[i].PRPR_NM = FF_2BODY[i].ATM1TYP;
		  FF_2BODY[i].PRPR_NM.append(FF_2BODY[i].ATM2TYP);	

		  if(TEMP_TYPE == "DFTBPOLY" || TEMP_TYPE == "INVRSE_R")
		  {
			 FF_2BODY[i].SNUM = TMP_TERMS1;
		  }
		  else if(TEMP_TYPE =="CHEBYSHEV")
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
			 FF_2BODY[0].FORCE_CUTOFF.BODIEDNESS = 2 ;
					
			 // Copy all class members.
			 for(int i=1; i<FF_2BODY.size(); i++)
				FF_2BODY[i].FORCE_CUTOFF = FF_2BODY[0].FORCE_CUTOFF;										
		  }
		  else if(TEMP_TYPE =="LJ")
		  {
			 FF_2BODY[i].SNUM = 2 ;
		  }
		  else if ( TEMP_TYPE == "CLUSTER" )
		  {
			 FF_2BODY[i].SNUM = 0 ;
		  }
		  else // Splines type
		  {
			 FF_2BODY[i].SNUM = (2+floor((FF_2BODY[i].S_MAXIM - FF_2BODY[i].S_MINIM)/FF_2BODY[i].S_DELTA))*2;
		  }				
		}

		if(CONTROLS.USE_OVERCOORD)
		{		
		  LINE = get_next_line(PARAMFILE) ;
		  LINE = get_next_line(PARAMFILE) ;
		  LINE = get_next_line(PARAMFILE) ;

		  for(int i=0; i<NO_PAIRS; i++)
		  {
			 PARAMFILE >> TEMP_STR >> TEMP_STR >> TEMP_STR;	
			 PARAMFILE >> FF_2BODY[i].USE_OVRPRMS;	
			 PARAMFILE >> FF_2BODY[i].OVER_TO_ATM;				
			 PARAMFILE >> FF_2BODY[i].OVRPRMS[0];	
			 PARAMFILE >> FF_2BODY[i].OVRPRMS[1];
			 PARAMFILE >> FF_2BODY[i].OVRPRMS[2];
			 PARAMFILE >> FF_2BODY[i].OVRPRMS[3];
			 PARAMFILE >> FF_2BODY[i].OVRPRMS[4];	
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
		
	 else if(LINE.find("PAIR CHEBYSHEV CUBIC SCALING: ") != string::npos)
	 {
		STREAM_PARSER.str(LINE);
		STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR;
		for(int i=0; i<NO_PAIRS; i++)
		  FF_2BODY[i].CUBIC_SCALE = double(atof(TEMP_STR.data()));
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
				
		  // If this is a spline type, we will over-write parameter coefficients with 
		  // a value less than 1.

		  int 	STOP_FILL_IDX = -1;
		  double 	slope          = -10.0;
				
		  if(FF_2BODY[i].PAIRTYP == "SPLINE")
		  {
			 if(i==0 && RANK ==0)
			 {
				cout << "		Will use a simple linear force model for poorly sampled positions" << endl;
				cout << "			i.e. for |spline coefficients| < 1.0 at close separation distance..." << endl;					
			 }
					
			 STOP_FILL_IDX = -1;
										
			 for(int j=0; j<FF_2BODY[i].SNUM; j++)
			 {
				if(fabs(FF_2BODY[i].PARAMS[j])>1.0)
				{
				  STOP_FILL_IDX = j;
				  if(RANK==0)
					 cout << "			...Generating linear params for pair idx " << i << " for params 0 through " << STOP_FILL_IDX-1 << endl;
				  break;
				}
			 }
								
			 if(STOP_FILL_IDX != -1)
			 {
				for(int j=0; j<STOP_FILL_IDX; j+=2)
				{
				  FF_2BODY[i].PARAMS[j  ] = slope * (STOP_FILL_IDX-j) + FF_2BODY[i].PARAMS[STOP_FILL_IDX];
				  FF_2BODY[i].PARAMS[j+1] = slope/FF_2BODY[i].S_DELTA;	
				  if(RANK==0)
					 cout << "			   " << j << " " << FF_2BODY[i].PARAMS[j] << endl << "			   " << FF_2BODY[i].PARAMS[j+1] << endl;
				}
			 }
					
			 // Compute the integral of the spline equation for use in analytical pressure
			 // calculations...
			 //  We are computing an integral, so we will only have half as many points
					
					
			 FF_2BODY[i].POT_PARAMS.resize(FF_2BODY[i].SNUM/2);
					
			 for(int j=0; j<FF_2BODY[i].POT_PARAMS.size(); j++) 
				FF_2BODY[i].POT_PARAMS[j] = 0;
					
			 for(int j=FF_2BODY[i].SNUM/2-2; j>=0; j--) 
			 {
				FF_2BODY[i].POT_PARAMS[j] = FF_2BODY[i].POT_PARAMS[j+1] - FF_2BODY[i].S_DELTA *
				  (FF_2BODY[i].PARAMS[j*2]/2 + FF_2BODY[i].S_DELTA * FF_2BODY[i].PARAMS[j*2+1]/12 + FF_2BODY[i].PARAMS[j*2+2]/2 - FF_2BODY[i].S_DELTA * FF_2BODY[i].PARAMS[j*2+3]/12);
			 }
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

				  //FF_2BODY[i].ATM1CHG = TMP_CHARGES[j];
				  //FF_2BODY[i].ATM2CHG = TMP_CHARGES[j];
					
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
		  if(!CONTROLS.PLOT_PES)
		  {
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
		
	 // Read in the fit overbonding parameter
		
	 else if(LINE.find("P OVER: ") != string::npos)
	 {
		if(CONTROLS.USE_OVERCOORD && CONTROLS.FIT_POVER)
		{
		  STREAM_PARSER.str(LINE);
		  STREAM_PARSER >> TEMP_STR >> TEMP_STR;
		  STREAM_PARSER >> FF_2BODY[0].OVRPRMS[0];
		  STREAM_PARSER.str("");
		  STREAM_PARSER.clear();	
			
		  for(int i=1; i<NO_PAIRS; i++)
			 FF_2BODY[i].OVRPRMS[0] = FF_2BODY[0].OVRPRMS[0];				
		}
			
		if (RANK==0)
		  cout << "	...Read ReaxFF overbonding FF params..." << endl;
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
	 else if ( tokens[0] == "USEPOVR:" )
		CONTROLS.USE_OVERCOORD = is_true(tokens[1]) ;
	 else if ( tokens[0] == "FITPOVR:" )
		CONTROLS.FIT_POVER = is_true(tokens[1]) ;
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
	 cout << "		...Compute ReaxFF overbonding?  " << boolalpha << CONTROLS.USE_OVERCOORD << endl;
	 cout << "		...Use fit overbonding param?   " << boolalpha << CONTROLS.FIT_POVER << endl;
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
		cout << "...penalty dist...penalty scaling...cubic scaling";
	 cout << endl;
				
	 cout << "		" << FF_2BODY[i].PRPR_NM << " ";
		
	 cout << FF_2BODY[i].S_MINIM << " ";
	 cout << FF_2BODY[i].S_MAXIM << " ";
	 cout << FF_2BODY[i].S_DELTA << " ";
		
	 if(FF_2BODY[i].PAIRTYP == "CHEBYSHEV")
	 {
		string cheby_type = Cheby::get_trans_string(FF_2BODY[i].CHEBY_TYPE) ;
		cout << cheby_type << " ";
			
		if(FF_2BODY[i].CHEBY_TYPE == Cheby_trans::MORSE )
		  cout << FF_2BODY[i].LAMBDA << " ";
		cout << FF_2BODY[i].PENALTY_DIST << " " << scientific << FF_2BODY[i].PENALTY_SCALE<< " ";
		cout << FF_2BODY[i].CUBIC_SCALE;
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
	
  if(CONTROLS.USE_OVERCOORD)
  {
	 cout << "	ReaxFF over-coordination parameters: " << endl;
		
	 cout << "		pair...pover...r0...p1...p2...lambda6" << endl;
		
	 for(int i=0; i<FF_2BODY.size(); i++)
	 {
		cout << "		" << FF_2BODY[i].PRPR_NM << " ";
		cout << FF_2BODY[i].OVRPRMS[0] << " ";
		cout << FF_2BODY[i].OVRPRMS[1] << " ";
		cout << FF_2BODY[i].OVRPRMS[2] << " ";
		cout << FF_2BODY[i].OVRPRMS[3] << " ";
		cout << FF_2BODY[i].OVRPRMS[4] << endl;
	 }
  }	

  cout << endl;
}

static void print_pes(PES_PLOTS &FF_PLOTS, vector<PAIRS> &FF_2BODY, map<string,int> &PAIR_MAP, map<int,string> &PAIR_MAP_REVERSE, CLUSTER_LIST& TRIPS, CLUSTER_LIST& QUADS, JOB_CONTROL &CONTROLS, Cheby &cheby)
// Print the potential energy surface. The quads argument is unused for now.
{

  vector<CLUSTER>& FF_3BODY = TRIPS.VEC ;
  
  int ij, ik, jk;
  string ATM_TYP_1, ATM_TYP_2, ATM_TYP_3;


  ifstream FULL_INFILE_3B;
  ifstream SCAN_INFILE_3B;
  ifstream SCAN_INFILE_2B;
		
  ofstream FULL_OUTFILE_3B;
  ofstream SCAN_OUTFILE_3B;
  ofstream SCAN_OUTFILE_2B;			

  vector<double> IJ_DIST_3B;
  vector<double> IK_DIST_3B;
  vector<double> JK_DIST_3B;
  vector<double> PES_VAL_3B;
		
  vector<double> IJ_DIST_2B;
  vector<double> IK_DIST_2B;
  vector<double> JK_DIST_2B;
  vector<double> PES_VAL_2B_IJ;
  vector<double> PES_VAL_2B_IK;
  vector<double> PES_VAL_2B_JK;
		
  int idx_3b = 0;
  int scan_2b_idx = 0;
		
  for (int i=0; i<FF_PLOTS.N_PLOTS; i++)
  {
	 cout << "	Printing force field PES for: " << FF_PLOTS.PES_TYPES[i] << endl;
			
	 // Figure out whether this is a 2- or 3-body interaction
	 // Figure out the index of the interaction type
	 // If it is 3-body, figure out:
	 // 1. the individual atom types
	 // 2. the pair type indicies
			
		
	 IJ_DIST_3B.clear();
	 IK_DIST_3B.clear();
	 JK_DIST_3B.clear();
	 PES_VAL_3B.clear();		
			
	 IJ_DIST_2B.clear();
	 IK_DIST_2B.clear();
	 JK_DIST_2B.clear();
			
	 PES_VAL_2B_IJ.clear();
	 PES_VAL_2B_IK.clear();
	 PES_VAL_2B_JK.clear();
								
			
	 if(FF_PLOTS.NBODY[i] == 3)	// Eventually also add a check that tye 3B type is Chebyshev...
	 {
		ij =  PAIR_MAP[FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATOM_PAIRS[0]];
		ik =  PAIR_MAP[FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATOM_PAIRS[1]];
		jk =  PAIR_MAP[FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATOM_PAIRS[2]];
			
		ATM_TYP_1 = FF_2BODY[ij].ATM1TYP;
		ATM_TYP_2 = FF_2BODY[jk].ATM1TYP;
		ATM_TYP_3 = FF_2BODY[jk].ATM2TYP;
				
		cout << "	Will work with pair types: " << FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATOM_PAIRS[0] << " " << FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATOM_PAIRS[1] << " " << FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATOM_PAIRS[2] << endl;
		cout << "	and atom types:            " << ATM_TYP_1 << " " << ATM_TYP_2 << " " << ATM_TYP_3 << endl;
				
				
		if(FF_PLOTS.DO_4D)
		{
		  cout << "ERROR: Functionality depreciated. Code needs updating " << endl; 
		  exit_run(0);
		}

		cheby.Print_3B(TRIPS, ATM_TYP_1, ATM_TYP_2, ATM_TYP_3, ij, ik, jk, FF_PLOTS, i);	
				
		scan_2b_idx++;				
	 }
	 else if(FF_PLOTS.NBODY[i] == 2)
	 {

		cout << "	Will work with pair types: " << PAIR_MAP_REVERSE[FF_PLOTS.TYPE_INDEX[i]] << endl;
		cout << "	and atom types:            " << FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM2TYP << endl;
					
		cheby.Print_2B(FF_PLOTS.TYPE_INDEX[i], PAIR_MAP_REVERSE[FF_PLOTS.TYPE_INDEX[i]], FF_PLOTS.INCLUDE_FCUT, FF_PLOTS.INCLUDE_CHARGES, FF_PLOTS.INCLUDE_PENALTY, "");
	 }
	 else
	 {
		cout << "Functionality not programmed yet." << endl;
	 }

	 cout << "	...Force field PES printing complete." << endl;
		
	 idx_3b++;
			
  }
}


static void read_coord_file(int index, JOB_CONTROL &CONTROLS, FRAME &SYSTEM, ifstream &CMPR_FORCEFILE)
// Read coordinates,  and optionally velocities and forces.  There is support for more than one coordinate
// input file, given by the index.
{
  ifstream coordfile;
	int natoms ;
	XYZ tmp_box ;
  stringstream	stream_parser;
	string first_ext;
  string  line;

	coordfile.open(CONTROLS.COORD_FILE[index].data());
	
	coordfile >> natoms;	// Store number of atoms in file in a temp var
	coordfile >> tmp_box.X >> tmp_box.Y >> tmp_box.Z;

	if(CONTROLS.FIT_STRESS)                                                                                           
		coordfile >>  SYSTEM.PRESSURE_TENSORS.X >>  SYSTEM.PRESSURE_TENSORS.Y >>  SYSTEM.PRESSURE_TENSORS.Z;      

	if(CONTROLS.FIT_ENER)
		coordfile >> SYSTEM.QM_POT_ENER;   

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
		first_ext = EXTENSION;
		
	if(EXTENSION != first_ext)
	{
		cout << "ERROR: Extensions for all input coordinate files must match. " << endl;
		cout << "	Found extension "<< first_ext << " for coordinate file 0 " << endl;
		cout << "	and " << EXTENSION << " for coordinate file " << index << endl;
		exit_run(0);
	}
	
	if (RANK==0)
	{
		cout << "     ...Read the following coordinate file extension: " << EXTENSION << endl;
		cout << "   ...Read the following number of atoms: " << natoms << endl;
		cout << "     ...Read box dimensions: " << tmp_box.X << " " << tmp_box.Y << " " << tmp_box.Z << endl;
		
		if(CONTROLS.FIT_STRESS)
			cout << "	...Read stress tensors: " << SYSTEM.PRESSURE_TENSORS.X << " " << SYSTEM.PRESSURE_TENSORS.Y << " " << SYSTEM.PRESSURE_TENSORS.Z << endl;
		if(CONTROLS.FIT_ENER)  													
			cout << "	...Read potential energy: " << SYSTEM.QM_POT_ENER << endl;					       
		
	}
	
	getline(coordfile,line);

	int curr_atom = -1 ;

	for(int a=0; a<natoms;a++)
	{
			
		curr_atom++;

		getline(coordfile,line);

		stream_parser.str(line);
		
		stream_parser >> SYSTEM.ATOMTYPE[curr_atom];

		// Read the coordinates

		stream_parser >> SYSTEM.COORDS[curr_atom].X >> SYSTEM.COORDS[curr_atom].Y >> SYSTEM.COORDS[curr_atom].Z;

		// Wrap the coordinates, shift along Z

		SYSTEM.COORDS[curr_atom].X -= floor(SYSTEM.COORDS[curr_atom].X/tmp_box.X)*tmp_box.X;
		SYSTEM.COORDS[curr_atom].Y -= floor(SYSTEM.COORDS[curr_atom].Y/tmp_box.Y)*tmp_box.Y;
		SYSTEM.COORDS[curr_atom].Z -= floor(SYSTEM.COORDS[curr_atom].Z/tmp_box.Z)*tmp_box.Z;
		SYSTEM.COORDS[curr_atom].Z += SYSTEM.BOXDIM.Z;

		// Prepare velocities
		SYSTEM.VELOCITY[curr_atom].X = 0;
		SYSTEM.VELOCITY[curr_atom].Y = 0;
		SYSTEM.VELOCITY[curr_atom].Z = 0;
		
		// Prepare forces
		SYSTEM.FORCES[curr_atom].X = 0;
		SYSTEM.FORCES[curr_atom].Y = 0;
		SYSTEM.FORCES[curr_atom].Z = 0;
		
		// Prepare accelerations
		SYSTEM.ACCEL[curr_atom].X = 0;
		SYSTEM.ACCEL[curr_atom].Y = 0;
		SYSTEM.ACCEL[curr_atom].Z = 0;		
		

		if ( CONTROLS.RESTART ) 
		{
			if(a==0 && RANK==0)
				cout << "	...Reading positions, velocities, and forces from a restart file " << endl ;

			stream_parser >> SYSTEM.VELOCITY[curr_atom].X >> SYSTEM.VELOCITY[curr_atom].Y  >> SYSTEM.VELOCITY[curr_atom].Z;	
			stream_parser >> SYSTEM.ACCEL   [curr_atom].X >> SYSTEM.ACCEL   [curr_atom].Y  >> SYSTEM.ACCEL   [curr_atom].Z;	
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

				//coordfile >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;			
			}
					
			// Read forces from separate force file
			CMPR_FORCEFILE >> SYSTEM.FORCES[curr_atom].X >> SYSTEM.FORCES[curr_atom].Y >> SYSTEM.FORCES[curr_atom].Z;
			
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

			stream_parser >> SYSTEM.VELOCITY[curr_atom].X >> SYSTEM.VELOCITY[curr_atom].Y >> SYSTEM.VELOCITY[curr_atom].Z;
		}
		else if(EXTENSION == "xyzf")
		{
			if(a==0 && RANK==0)
				cout << "	...Ignoring last three fields of atom info lines in xyzf file... " << endl;
			//coordfile >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;
		}
		
		stream_parser.str("");
		stream_parser.clear();	
	}

	// Read in stress components for comparison to current values.
	if ( CONTROLS.FIT_STRESS ) 
	{
		// Units read in are kcal/mol-Ang.^3.
		CMPR_FORCEFILE >> SYSTEM.STRESS_TENSORS.X >> SYSTEM.STRESS_TENSORS.Y >> SYSTEM.STRESS_TENSORS.Z ;

		if ( ! CMPR_FORCEFILE.good() ) 
			EXIT_MSG("Reading force comparson file failed") ;
	}

	// Input configurations are added sequentially along z.
		
	SYSTEM.BOXDIM.X  = tmp_box.X;
	SYSTEM.BOXDIM.Y  = tmp_box.Y;
	SYSTEM.BOXDIM.Z += tmp_box.Z;
		
	if(CONTROLS.SCALE_SYSTEM_BY != 1.0)
	{
		for(int a1=0; a1<SYSTEM.ATOMS; a1++)
		{
			SYSTEM.COORDS[a1].X *= CONTROLS.SCALE_SYSTEM_BY;
			SYSTEM.COORDS[a1].Y *= CONTROLS.SCALE_SYSTEM_BY;
			SYSTEM.COORDS[a1].Z *= CONTROLS.SCALE_SYSTEM_BY;
		}
			
		SYSTEM.BOXDIM.X *= CONTROLS.SCALE_SYSTEM_BY;
		SYSTEM.BOXDIM.Y *= CONTROLS.SCALE_SYSTEM_BY;
		SYSTEM.BOXDIM.Z *= CONTROLS.SCALE_SYSTEM_BY;
	}
		
	if (RANK==0)
		cout << "	...Updated simulation box dimensions: " << SYSTEM.BOXDIM.X << " " << SYSTEM.BOXDIM.Y << " " << SYSTEM.BOXDIM.Z << endl;

	if ( ! CONTROLS.RESTART ) 
	{
		coordfile.close();
		coordfile.clear();
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
		FORCE_SUBTRACTED_OUTPUT << SYSTEM.BOXDIM.X << " " << SYSTEM.BOXDIM.Y << " " << SYSTEM.BOXDIM.Z << " ";
				
		// NOTE:  PRESSURE_TENSORS_XYZ hold the potential energy contribution to the current stress tensor.
		// PRESSURE_TENSOR holds potential + kinetic energy to the current stress tensor in GPa units.
		// PRESSURE_TENSORS_XYZ should be used here.
		if(CONTROLS.FIT_STRESS)
			FORCE_SUBTRACTED_OUTPUT << 	SYSTEM.PRESSURE_TENSORS_XYZ.X - SYSTEM.STRESS_TENSORS.X 
															<< " " << SYSTEM.PRESSURE_TENSORS_XYZ.Y - SYSTEM.STRESS_TENSORS.Y 
															<< " " << SYSTEM.PRESSURE_TENSORS_XYZ.Z - SYSTEM.STRESS_TENSORS.Z;
		if(CONTROLS.FIT_ENER)
			FORCE_SUBTRACTED_OUTPUT << SYSTEM.TOT_POT_ENER - SYSTEM.QM_POT_ENER;
					
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
	if ( CONTROLS.FIT_STRESS ) {
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ.X - SYSTEM.STRESS_TENSORS.X) * (SYSTEM.PRESSURE_TENSORS_XYZ.X - SYSTEM.STRESS_TENSORS.X) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ.Y - SYSTEM.STRESS_TENSORS.Y) * (SYSTEM.PRESSURE_TENSORS_XYZ.Y - SYSTEM.STRESS_TENSORS.Y) ;
		ferr += (SYSTEM.PRESSURE_TENSORS_XYZ.Z - SYSTEM.STRESS_TENSORS.Z) * (SYSTEM.PRESSURE_TENSORS_XYZ.Z - SYSTEM.STRESS_TENSORS.Z) ;

		ferr = sqrt(ferr/3.0) ;
		ferr *= GPa ;
		cout << "RMS stress error = " << ferr << " GPa " << endl ;
	}
}
