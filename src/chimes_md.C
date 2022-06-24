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
	
using namespace std;	
	
// Define function headers -- general
static void read_input        (string & INFILE, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST);
static void read_atom_types(ifstream &PARAMFILE, JOB_CONTROL &CONTROLS, int &NATMTYP, vector<string>& TMP_ATOMTYPE, vector<int>& TMP_NATOMTYPE, vector<int>& TMP_ATOMTYPEIDX, vector<double>& TMP_CHARGES,  vector<double>& TMP_MASS, vector<int> &TMP_SIGN) ;
static void read_ff_params(ifstream &PARAMFILE, JOB_CONTROL &CONTROLS, vector<PAIR_FF>& FF_2BODY, CLUSTER_LIST& TRIPS, CLUSTER_LIST &QUADS, map<string,int> &PAIR_MAP, NEIGHBORS &NEIGHBOR_LIST, FRAME& SYSTEM, int NATMTYP, const vector<string>& TMP_ATOMTYPE, const vector<int>& TMP_ATOMTYPEIDX, vector<double> &TMP_CHARGES, vector<double> &TMP_MASS, const vector<int>& TMP_SIGN, map<int,string>& PAIR_MAP_REVERSE) ;
static void print_ff_summary(const vector<PAIR_FF> &FF_2BODY, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, const JOB_CONTROL &CONTROLS) ;
static void final_output(FRAME &SYSTEM, THERMO_AVG &AVG_DATA, JOB_CONTROL &CONTROLS, NEIGHBORS &NEIGHBOR_LIST, ofstream &STATISTICS,
												 CONSTRAINT &ENSEMBLE_CONTROL) ;

// Define function headers -- MPI

static void write_xyzv(FRAME &SYSTEM, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL,
											 THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBOR_LIST, string filename, bool restart);
static void read_restart_params(ifstream &COORDFILE, JOB_CONTROL &CONTROLS, CONSTRAINT &ENSEMBLE_CONTROL, THERMO_AVG &AVG_DATA, NEIGHBORS &NEIGHBORS, FRAME &SYSTEM) ;
static void parse_ff_controls  (string &LINE, ifstream &PARAMFILE, JOB_CONTROL &CONTROLS ) ;
static void read_coord_file(int index, JOB_CONTROL &CONTROLS, FRAME &SYSTEM, ifstream &CMPR_FORCEFILE) ;
static void subtract_force(FRAME &SYSTEM, JOB_CONTROL &CONTROLS) ;
static void print_for_dftbplus(FRAME &SYSTEM, JOB_CONTROL &CONTROLS);

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


// Could be defined locally now that LAMMPS contributions have been removed

CLUSTER_LIST TRIPS;                   // Holds all 3-body parameters
CLUSTER_LIST QUADS;                   // Holds all 4-body parameters
 
	
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

  double Ktot ;
	vector<XYZ> Ktensor(3,{0.0,0.0,0.0}) ;
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

  const int STAT_BUF_SZ = 256 ; // Max length of a statistics output line.
	
  cout.precision(15);		// Set output precision
	

	
  ////////////////////////////////////////////////////////////
  // Setup optional output files 
  ////////////////////////////////////////////////////////////

  // trajectory

  WRITE_TRAJ GENFILE(CONTROLS.TRAJ_FORMAT,"STANDARD");	// This is our usual output file

  // force 
  
  WRITE_TRAJ FORCEFILE;

  if ( CONTROLS.PRINT_FORCE ) 
	  FORCEFILE.INIT("XYZF_FORCE","FORCE",CONTROLS.PRINT_ENERGY_STRESS);	// This is the forcout* files

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
	 SYSTEM.COORDS0      .resize(SYSTEM.ATOMS);	 
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
  
  if(  (ENSEMBLE_CONTROL.STYLE == "NPT-BEREND")
		 ||(ENSEMBLE_CONTROL.STYLE == "NPT-BEREND-ANISO")
		 ||(ENSEMBLE_CONTROL.STYLE == "NPT-MTK") )
  {
		 SYSTEM.BOXDIM.IS_VARIABLE = true ;
		 if (!SYSTEM.BOXDIM.IS_ORTHO)
				EXIT_MSG("ERROR: NPT-BEREND, NPT-BEREND-ANISO, and NPT-MTK only compatible with orthorhombic cells.");
	}

  if ( CONTROLS.RESTART ) 
	 read_restart_params(COORDFILE, CONTROLS, ENSEMBLE_CONTROL, AVG_DATA, NEIGHBOR_LIST, SYSTEM) ;

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

  Vol = SYSTEM.BOXDIM.VOL;
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
  	
  NEIGHBOR_LIST.INITIALIZE_MD(SYSTEM,CONTROLS);
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
  //			START THE SIMULATION
  //
  ////////////////////////////////////////////////////////////  
  ////////////////////////////////////////////////////////////

  STATISTICS.open("md_statistics.out");
	STATISTICS.precision(6) ;
  
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
			ENSEMBLE_CONTROL.UPDATE_COORDS(SYSTEM, CONTROLS, NEIGHBOR_LIST);	// Update coordinates and ghost atoms
			
		if(CONTROLS.WRAP_COORDS)				// Wrap the coordinates:
		{
		 	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
				SYSTEM.BOXDIM.WRAP_ATOM(SYSTEM.COORDS[a1], SYSTEM.WRAP_IDX[a1], false);
		} 

		ENSEMBLE_CONTROL.UPDATE_VELOCS_HALF_1(SYSTEM, CONTROLS);// Update first half of velocity and max velocity for neighbor lists:		
	 }
		
#ifdef USE_MPI
	 sync_position(SYSTEM.COORDS    , NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ATOMS    , true,
								 SYSTEM.BOXDIM);

	 // Its faster to recalculate the ghost atoms than to communicate them with MPI.
	 if ( RANK != 0 ) SYSTEM.update_ghost(CONTROLS.N_LAYERS, false) ;

	 //sync_position(SYSTEM.ALL_COORDS, NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ALL_ATOMS, false,
	 //SYSTEM.BOXDIM);	
	 //MPI_Bcast(&NEIGHBOR_LIST.MAX_VEL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);				// Sync the maximum velocites so all procs use the right padding to generate their neigh lists
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

			char stat_buf[STAT_BUF_SZ] ;
			snprintf(stat_buf, STAT_BUF_SZ, "%8s %14s %14s %14s %14s %14s %14s ",
							 "# Step", "Time",  "Ktot/N",  "Vtot/N", 	"Etot/N", "T", "P");

			STATISTICS << stat_buf ;
	  
		  if ( ENSEMBLE_CONTROL.STYLE == "NVT-MTK" ) 
		  {
			 printf(" %15s\n", "Econs/N");

			 snprintf(stat_buf, STAT_BUF_SZ, "%14s %14s\n", "Econs/N", "P_conf") ;
			 STATISTICS << stat_buf ;
		  }
			else if ( ENSEMBLE_CONTROL.STYLE == "NPT-MTK" )
			{
				printf("%15s %15s\n", "Econs/N", "Volume");
				 
				snprintf(stat_buf, STAT_BUF_SZ, "%14s %14s %14s\n", "Econs/N", "P_conf", "Volume") ;
				 STATISTICS << stat_buf ;
			}
			else if ( ENSEMBLE_CONTROL.STYLE == "NPT-BEREND" || ENSEMBLE_CONTROL.STYLE == "NPT-BEREND-ANISO" )
			{
				printf("%15s\n", "Volume");
				 
				snprintf(stat_buf, STAT_BUF_SZ, "%14s %14s\n", "P_conf", "Volume") ;
				 STATISTICS << stat_buf ;
			}
		  else 
		  {
			 printf("\n");
			 snprintf(stat_buf, STAT_BUF_SZ, "%14s\n", "P_conf") ;			 
			 STATISTICS << stat_buf ;
		  }
			
		  printf("%8s %9s %15s %15s %15s %15s %15s", " ", "(fs)", "(kcal/mol)", "(kcal/mol)",
						 "(kcal/mol)", "(K)", "(GPa)");

			snprintf(stat_buf, STAT_BUF_SZ, "%8s %14s %14s %14s %14s %14s %14s ",
							 "# ", "(fs)", "(kcal/mol)", "(kcal/mol)", "(kcal/mol)", "(K)", "(GPa)");

			STATISTICS << stat_buf ;
			
		  if ( ENSEMBLE_CONTROL.STYLE == "NVT-MTK" ) 
				{
					printf(" %15s\n", "(kcal/mol)");

					snprintf(stat_buf, STAT_BUF_SZ, "%14s %14s\n", "(kcal/mol)", "(GPa)") ;			 
					STATISTICS << stat_buf ;
				}
			else if ( ENSEMBLE_CONTROL.STYLE == "NPT-MTK" )
				{
					printf(" %15s %15s\n", "(kcal/mol)", "(Ang.^3)");
					
					snprintf(stat_buf, STAT_BUF_SZ, "%14s %14s %14s\n", "(kcal/mol)", "(GPa)", "(Ang.^3)") ;			 
					STATISTICS << stat_buf ;
				}
			else if ( ENSEMBLE_CONTROL.STYLE == "NPT-BEREND" || ENSEMBLE_CONTROL.STYLE == "NPT-BEREND-ANISO" )
				{
					printf(" %15s\n", "(Ang.^3)");
					
					snprintf(stat_buf, STAT_BUF_SZ, "%14s %14s\n", "(GPa)", "(Ang.^3)") ;			 
					STATISTICS << stat_buf ;
				}			
		  else 
				{
					printf("\n");
					snprintf(stat_buf, STAT_BUF_SZ, "%14s\n", "(GPa)") ;			 			 
					STATISTICS << stat_buf ;
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

		 Ktot = kinetic_energy(SYSTEM, CONTROLS, Ktensor);

		 ENSEMBLE_CONTROL.UPDATE_TEMPERATURE(SYSTEM, CONTROLS);

		 // Exit with an error if the set and block temperatures differ by more than CONTROLS.NVT_CONV_CUT
		
		 AVG_DATA.TEMP_SUM += SYSTEM.TEMPERATURE;

		 if ( CONTROLS.NVT_CONV_CUT > 0.0 )
		 {
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
		double pe1_sum = 0.0, pe2_sum = 0.0 ;
		MPI_Allreduce(&PE_1,&pe1_sum,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		PE_1 = pe1_sum ;
		MPI_Allreduce(&PE_2,&pe2_sum,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		PE_2 = pe2_sum ;
#endif
			
		SYSTEM.PRESSURE_XYZ = -1.0*(PE_2 - PE_1)/dV;
	 }

	 if ( CONTROLS.USE_NUMERICAL_STRESS )
	 {
		 numerical_stress_all(SYSTEM, CONTROLS, FF_2BODY, TRIPS, QUADS, PAIR_MAP, INT_PAIR_MAP,
													NEIGHBOR_LIST) ;
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
	 
	 
	Vol = SYSTEM.BOXDIM.VOL;
	 	 
	 if(RANK==0)
	 {
	
		SYSTEM.PRESSURE = (SYSTEM.PRESSURE_XYZ + 2.0 * Ktot / (3.0 * Vol)) * GPa;	// GPa = Unit conversion factor to GPa.

		// Use the kinetic energy tensor for the instantaneous pressure tensor.
		SYSTEM.PRESSURE_TENSORS_ALL[0].X = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X + 2.0 * Ktensor[0].X / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[0].Y = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Y + 2.0 * Ktensor[0].Y / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[0].Z = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Z + 2.0 * Ktensor[0].Z / Vol) * GPa;
		
		SYSTEM.PRESSURE_TENSORS_ALL[1].X = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X + 2.0 * Ktensor[1].X / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[1].Y = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y + 2.0 * Ktensor[1].Y / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[1].Z = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Z + 2.0 * Ktensor[1].Z / Vol) * GPa;
		
		SYSTEM.PRESSURE_TENSORS_ALL[2].X = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].X + 2.0 * Ktensor[2].X / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[2].Y = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y + 2.0 * Ktensor[2].Y / Vol) * GPa;
		SYSTEM.PRESSURE_TENSORS_ALL[2].Z = (SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z + 2.0 * Ktensor[2].Z / Vol) * GPa;

		AVG_DATA.PRESS_SUM += SYSTEM.PRESSURE;
		AVG_DATA.PV_SUM += Vol * SYSTEM.PRESSURE/GPa ;

		AVG_DATA.STRESS_TENSOR_SUM_ALL[0].X += SYSTEM.PRESSURE_TENSORS_ALL[0].X;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[0].Y += SYSTEM.PRESSURE_TENSORS_ALL[0].Y;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[0].Z += SYSTEM.PRESSURE_TENSORS_ALL[0].Z;
		
		AVG_DATA.STRESS_TENSOR_SUM_ALL[1].X += SYSTEM.PRESSURE_TENSORS_ALL[1].X;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Y += SYSTEM.PRESSURE_TENSORS_ALL[1].Y;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Z += SYSTEM.PRESSURE_TENSORS_ALL[1].Z;
		
		AVG_DATA.STRESS_TENSOR_SUM_ALL[2].X += SYSTEM.PRESSURE_TENSORS_ALL[2].X;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[2].Y += SYSTEM.PRESSURE_TENSORS_ALL[2].Y;
		AVG_DATA.STRESS_TENSOR_SUM_ALL[2].Z += SYSTEM.PRESSURE_TENSORS_ALL[2].Z;
				
		AVG_DATA.VOLUME_SUM += Vol ;

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

			char stat_buf[STAT_BUF_SZ] ;

			snprintf(stat_buf, STAT_BUF_SZ, "%8d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e ",
					 CONTROLS.STEP+1,
					 (CONTROLS.STEP+1)*CONTROLS.DELTA_T_FS,
					 Ktot/SYSTEM.ATOMS,
					 SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS,
					 (Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS,
					 SYSTEM.TEMPERATURE,
					 SYSTEM.PRESSURE) ;
		
			STATISTICS << stat_buf ;


			// Print the econs value
			if ( ENSEMBLE_CONTROL.STYLE == "NVT-MTK" ) 
			{
				snprintf(stat_buf, STAT_BUF_SZ, "%14.7e %14.7e ", CONTROLS.IO_ECONS_VAL, SYSTEM.PRESSURE_XYZ * GPa) ;
				STATISTICS << stat_buf ;
			}
			else if ( ENSEMBLE_CONTROL.STYLE == "NPT-MTK" )
			{
				snprintf(stat_buf, STAT_BUF_SZ, "%14.7e %14.7e %14.7e ", CONTROLS.IO_ECONS_VAL, SYSTEM.PRESSURE_XYZ * GPa, SYSTEM.BOXDIM.VOL) ;
				STATISTICS << stat_buf ;
			}
			else if ( ENSEMBLE_CONTROL.STYLE == "NPT-BEREND" || ENSEMBLE_CONTROL.STYLE == "NPT-BEREND-ANISO" )
			{
				snprintf(stat_buf, STAT_BUF_SZ, "%14.7e %14.7e ", SYSTEM.PRESSURE_XYZ * GPa, SYSTEM.BOXDIM.VOL) ;
				STATISTICS << stat_buf ;
			}			
			else 
			{
				snprintf(stat_buf, STAT_BUF_SZ, "%14.7e ", SYSTEM.PRESSURE_XYZ * GPa) ;				
				STATISTICS << stat_buf ;
			}

			if ( ENSEMBLE_CONTROL.STYLE == "NVT-MTK" || ENSEMBLE_CONTROL.STYLE == "NPT-MKT" )
			{
				printf("%15.7f", CONTROLS.IO_ECONS_VAL);
			} 
			else if ( ENSEMBLE_CONTROL.STYLE == "NPT-MTK"
					  || ENSEMBLE_CONTROL.STYLE == "NPT-BEREND"
					  || ENSEMBLE_CONTROL.STYLE == "NPT-BEREND-ANISO" )
			{
				printf(" %15.7f", SYSTEM.BOXDIM.VOL) ;
			}

			printf("\n") ;
			STATISTICS << endl ;
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
			sync_position(SYSTEM.COORDS    , NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ATOMS, true, SYSTEM.BOXDIM);	// Sync the main coords. Don't need to sync veloc since only proc 1 handles constraints
			sync_position(SYSTEM.ALL_COORDS, NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ALL_ATOMS, false, SYSTEM.BOXDIM);	// Sync the main coords. Don't need to sync veloc since only proc 1 handles constraints

			//MPI_Bcast(&NEIGHBOR_LIST.MAX_VEL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);					// Sync the maximum velocites so all procs use the right padding to generate their neigh lists
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

	final_output(SYSTEM, AVG_DATA, CONTROLS, NEIGHBOR_LIST, STATISTICS, ENSEMBLE_CONTROL) ;
	
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
				if ( ntokens < 4 )
					EXIT_MSG("MISSING NUMBER OF ENERGY OFFSETS\n") ;

				int num_offsets = stoi(tokens[3]) ;

				if ( num_offsets != CONTROLS.ATOMTYPES.size() )
					EXIT_MSG("AN ENERGY OFFSET IS REQUIRED FOR EVERY ATOM TYPE") ;
					
				SYSTEM.QM_ENERGY_OFFSET.resize(num_offsets) ;
				
				for (int i=0; i< num_offsets ; i++)
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
			SYSTEM.QM_ENERGY_OFFSET.resize(CONTROLS.ATOMTYPES.size(),0.0) ;
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

			 // Test for older style input with S_DELTA.
			 stringstream s_test(cheby_type) ;
			 double d_test ;

			 s_test >> d_test ;
			 if ( ! s_test.fail() )
			 {
				 // Numeric input found where cheby type should be.
				 if ( RANK == 0 ) cout << "Detected S_DELTA specification in Cheby pair parameters (not used)\n" ;
				 PARAMFILE >> cheby_type ;
                 						 
			 }
						 

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
		  else if(TEMP_TYPE =="LJ")
		  {
			  FF_2BODY[i].SNUM = 2 ;
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

	if ( CONTROLS.SERIAL_CHIMES )
	{
		if ( RANK == 0 )
			cout << "Will read serial chimes calculator file from " << CONTROLS.PARAM_FILE ;
		 
		FF_2BODY[0].chimes.init_chimesFF(CONTROLS.PARAM_FILE,RANK) ;
		FF_2BODY[0].PAIRTYP = "SERIAL_CHIMES" ;
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
		TMP_BOX.CELL_AX = stod(LINE);
	
	
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
	{
			SYSTEM.PRESSURE_TENSORS_ALL.resize(3) ;
			COORDFILE >>  SYSTEM.PRESSURE_TENSORS_ALL[0].X >>  SYSTEM.PRESSURE_TENSORS_ALL[1].Y >>  SYSTEM.PRESSURE_TENSORS_ALL[2].Z;
	}
	
	if (CONTROLS.FIT_STRESS_ALL)
  {
			SYSTEM.PRESSURE_TENSORS_ALL.resize(3) ;		
			COORDFILE >>  SYSTEM.PRESSURE_TENSORS_ALL[0].Y >>  SYSTEM.PRESSURE_TENSORS_ALL[0].Z >>  SYSTEM.PRESSURE_TENSORS_ALL[1].Z;
	}

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

		SYSTEM.PRESSURE_TENSORS_ALL.resize(3);		
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

		if ( ! CONTROLS.INCLUDE_ATOM_OFFSETS )
		{
			SYSTEM.SET_NATOMS_OF_TYPE();
		
			for(int i=0; i<SYSTEM.QM_ENERGY_OFFSET.size(); i++)
				ferr += SYSTEM.QM_ENERGY_OFFSET[i]*SYSTEM.NATOMS_OF_TYPE[i];
		
		}
		ferr = fabs(ferr/SYSTEM.ATOMS);
		cout << "Absolute energy error = " << ferr << " kcal/mol/atom \n" ;
	}
}

static void final_output(FRAME &SYSTEM, THERMO_AVG &AVG_DATA, JOB_CONTROL &CONTROLS, NEIGHBORS &NEIGHBOR_LIST, ofstream &STATISTICS,
												 CONSTRAINT &ENSEMBLE_CONTROL)
// Final output at the end of an MD simulation.
{
  if (RANK==0)	
  {
		cout << "END SIMULATION" << endl;

		double TEMP_MASS = 0.0;
	
		for ( int a = 0; a < SYSTEM.ATOMS; a++ ) 
			TEMP_MASS  += SYSTEM.MASS[a];


		cout << "	Average temperature over run = " << fixed << setprecision(4) << right << AVG_DATA.TEMP_SUM  / CONTROLS.N_MD_STEPS << " K"   << endl;
		cout << "	Average pressure    over run = " << fixed << setprecision(4) << right << AVG_DATA.PRESS_SUM / CONTROLS.N_MD_STEPS << " GPa" << endl;
		
		// Stress tensors are not calculated for all potentials.
		double Pavg = (AVG_DATA.STRESS_TENSOR_SUM_ALL[0].X + AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Y + AVG_DATA.STRESS_TENSOR_SUM_ALL[2].Z)/3.0/ CONTROLS.N_MD_STEPS ;
		cout << "	Pressures from diagonal stress tensors over run: " << Pavg << endl;
		cout << "	Average stress tensors over run: " << endl;
		cout << "		sigma_xx: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[0].X/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_yy: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Y/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_zz: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[2].Z/CONTROLS.N_MD_STEPS << " GPa" << endl; 
		cout << "		sigma_xy: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[0].Y/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_xz: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[0].Z/CONTROLS.N_MD_STEPS << " GPa" << endl;
		cout << "		sigma_yz: " << AVG_DATA.STRESS_TENSOR_SUM_ALL[1].Z/CONTROLS.N_MD_STEPS << " GPa" << endl;	   

		if ( SYSTEM.BOXDIM.IS_VARIABLE )
			{
				cout << "        Average volume over run = " << fixed << setprecision(4) << right << AVG_DATA.VOLUME_SUM / CONTROLS.N_MD_STEPS << " Ang.^3" << endl ;
				cout << "        Average PV over run     = " << fixed << setprecision(4) << right << AVG_DATA.PV_SUM / CONTROLS.N_MD_STEPS  << " kcal/mol "  
						 << endl ;
				if ( ENSEMBLE_CONTROL.STYLE == "NPT-MTK" )
					{
						// Allows a check on pressure and volume fluctuation magnitude, which should be correctly reproduced
						// by NPT-NTK algorithm.
						// See Martyna, Tobias, Klein JCP 101, 4177(1994) Appendix A.					 
						cout << "        Average of PV predicted by virial theorem = " << fixed << setprecision(4) << right <<
							(CONTROLS.PRESSURE / GPa) * AVG_DATA.VOLUME_SUM / CONTROLS.N_MD_STEPS - Kb * CONTROLS.TEMPERATURE
								 << " kcal/mol" << endl ;
					}
			}
	 
		// Write the final configuration to file.
		write_xyzv(SYSTEM, CONTROLS, ENSEMBLE_CONTROL, AVG_DATA, NEIGHBOR_LIST, "output.xyz", false);
			
		STATISTICS.close();
	} // if ( RANK == 0 )
}
