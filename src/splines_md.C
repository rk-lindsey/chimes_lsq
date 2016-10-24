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
#include<string.h>
#include<sstream>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes

// Include our own custom header

#include "functions.h"

// No fancy namespaces used here

using namespace std;

// Define function headers -- general

static void   read_input        (MD_JOB_CONTROL & CONTROLS, PES_PLOTS & FF_PLOTS, NEIGHBORS & NEIGHBOR_LIST);	// UPDATED
static double kinetic_energy    (FRAME & SYSTEM);				// UPDATED
static double kinetic_energy    (FRAME & SYSTEM, string TYPE);	// UPDATED
double        numerical_pressure(const FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, NEIGHBORS & NEIGHBOR_LIST);

// Define function headers -- MPI

static void sum_forces(vector<XYZ>& accel_vec, int atoms, double &pot_energy, double &pressure);
static void sync_position(vector<XYZ>& coord_vec, NEIGHBORS & neigh_list, vector<XYZ>& velocity_vec, int atoms, bool sync_vel);

// Global variables declared as externs in functions.h, and declared in functions.C -- general

string FULL_FILE_3B;		
string SCAN_FILE_3B;
string SCAN_FILE_2B;

// Global variables declared as externs in functions.h, and declared in functions.C -- MPI calculations.   
 
int NPROCS;		// Number of processors
int RANK;		// Index of current processor


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


	////////////////////////////////////////////////////////////
    // Define/initialize important variables
	////////////////////////////////////////////////////////////

	MD_JOB_CONTROL CONTROLS;			// Declare the data object that will hold the main simulation control variables
	PES_PLOTS FF_PLOTS;					// Declare the data object that will help set up PES plots
	
	// Data objects to hold coefficients for different force field types, and for FF printing (if requested)
	 
	vector<PAIR_FF> FF_2BODY;	// Holds all 2-body parameters
	vector<TRIP_FF> FF_3BODY; 	// Holds all 3-body parameters

	// Define the mapping variables that let us figure out which FF params to use for a given pair/triplet of pairs
	
	map<string,int> PAIR_MAP;			// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
	map<int,string> PAIR_MAP_REVERSE; 	// Input/output the resverse of PAIR_MAP
	map<string,int> TRIAD_MAP;			// Input is three of any ATOM_PAIRS.PRPR_NM pairs, output is a triplet type index
	map<int,string> TRIAD_MAP_REVERSE;	// Input/output is the reverse of TRIAD_MAP

	FRAME     SYSTEM;					// Declare the data object that will hold the system coordinates, velocities, accelrations, etc.
	NEIGHBORS NEIGHBOR_LIST;			// Declare the class that will handle the neighbor list
	
	double Ktot;
	double Vol;
    double avg_temp  = 0.0;
    double temp_sum  = 0.0;
    double press_sum = 0.0;
	
	double ferr = 0.0;
		
    double dens_mol;
    double dens_mass;
	
	double vscale;

	bool   	FOUND_END = false;
	string 	LINE;
	string  TEMP_STR;
	int		NATMTYP = 0;
	
	vector<string> 	TMP_ATOMTYPE;
	vector<int>		TMP_SIGN;
	vector<int>		TMP_NATOMTYPE;
	vector<int> 	TMP_ATOMTYPEIDX;
	vector<double> 	TMP_CHARGES;
	vector<double> 	TMP_MASS;
	stringstream	STREAM_PARSER;
	
	XYZ TEMP_XYZ;
	int TEMP_IDX;
	XYZ_INT TEMP_LAYER;
	
	ofstream GENFILE;			// Holds dftbgen info output.. whatever that is
	ifstream CMPR_FORCEFILE;	// Holds the forces that were read in for comparison purposes
	ofstream OUT_FORCEFILE;		// Holds the forces that are computed and are to be printed out
	ofstream OUT_FORCELABL;		// Holds the forces that are computed and are to be printed out.. has atom labels
	ofstream OUT_VELOCFILE;		// Holds the velocities that are computed and are to be printed out
	ofstream STATISTICS;		// Holds a copy of the "Step      Time          Ktot/N          Vtot/N ..." part of the output file.. uninterrupted by warnings
	
	string EXTENSION;			// Holds the input xyz/xyzf file extension, so we know how many fields to read on each atom line
	
	ifstream PARAMFILE;
	ifstream COORDFILE;


	
	read_input(CONTROLS, FF_PLOTS,NEIGHBOR_LIST);		// Populate object with user defined values
	
	cout.precision(15);			// Set output precision
	
	////////////////////////////////////////////////////////////
    // Hop to PES printing, if requested 
	////////////////////////////////////////////////////////////
	
	if (FF_PLOTS.N_PLOTS > 0)
		goto FF_SETUP_1; 
	
	////////////////////////////////////////////////////////////
	// Setup optional output files 
	////////////////////////////////////////////////////////////

	// dftbgen
	
    if ( CONTROLS.FREQ_DFTB_GEN > 0 ) 
		GENFILE.open("traj.gen");
	
	// force 

    if ( CONTROLS.PRINT_FORCE ) 
	{
		OUT_FORCEFILE.open("forceout.txt");	
		OUT_FORCELABL.open("forceout-labeled.txt");		
	}

	// velocity

    if ( CONTROLS.PRINT_VELOC ) 
		OUT_VELOCFILE.open("velocout.txt");	


	////////////////////////////////////////////////////////////
    // Read input file input.xyz, where box dims are on the info line:
	////////////////////////////////////////////////////////////

	COORDFILE.open(CONTROLS.COORD_FILE.data());
	if(!COORDFILE.is_open())
	{
		cout << "ERROR: Cannot open coordinate file: " << CONTROLS.COORD_FILE << endl;
		exit_run(0);
	}
	
	////////////////////////////////////////////////////////////
	// Set up the system data object... Do so by assuming a 
	// constant number of atoms in the system (almost always
	// true for md)...
	////////////////////////////////////////////////////////////
	
	if (RANK==0)
		cout << "Setting up the data objects..." << endl;

	COORDFILE >> SYSTEM.ATOMS;	
	COORDFILE >> SYSTEM.BOXDIM.X >> SYSTEM.BOXDIM.Y >> SYSTEM.BOXDIM.Z;
	
	SYSTEM.ATOMTYPE    .resize(SYSTEM.ATOMS);
	SYSTEM.COORDS      .resize(SYSTEM.ATOMS);
	SYSTEM.CHARGES     .resize(SYSTEM.ATOMS);
	SYSTEM.MASS        .resize(SYSTEM.ATOMS);
	SYSTEM.FORCES      .resize(SYSTEM.ATOMS);
	SYSTEM.ACCEL       .resize(SYSTEM.ATOMS);
	SYSTEM.VELOCITY    .resize(SYSTEM.ATOMS);
	SYSTEM.VELOCITY_NEW.resize(SYSTEM.ATOMS);
	SYSTEM.ATOMTYPE_IDX.resize(SYSTEM.ATOMS);

	if (RANK==0)
		cout << "   ...setup complete" << endl << endl;

	////////////////////////////////////////////////////////////
	// Read in the initial system coordinates, and if requested,
	// initial forces from separate file (i.e. not from .xyzf)
	// ... We also need to figure out the file extension so we
	// know how many fields to read on each atom lne
	////////////////////////////////////////////////////////////

	if (RANK==0)
		cout << "Reading initial coordinates and forces..." << endl;

    if ( CONTROLS.COMPARE_FORCE ) 
    {
		if (RANK==0)
			cout << "	Opening " << CONTROLS.COMPARE_FILE.data() << " to read forces for comparison\n";
		CMPR_FORCEFILE.open(CONTROLS.COMPARE_FILE.data());

		if(!CMPR_FORCEFILE.is_open())
		{
			cout << "ERROR: Cannot open force input file: " << CONTROLS.COMPARE_FILE << endl;
			exit_run(0);
		}

    }
	
	EXTENSION = CONTROLS.COORD_FILE.substr(CONTROLS.COORD_FILE.find_last_of(".")+1);
	
	if (RANK==0)
	{
		cout << "	...Read the following coordinate file extension: " << EXTENSION << endl;
    	cout << "	...Read box dimensions: " << SYSTEM.BOXDIM.X << " " << SYSTEM.BOXDIM.Y << " " << SYSTEM.BOXDIM.Z << endl;
	}
	
	getline(COORDFILE,LINE);
	
    for(int a=0; a<SYSTEM.ATOMS;a++)
	{
		getline(COORDFILE,LINE);

		STREAM_PARSER.str(LINE);
		
        STREAM_PARSER >> SYSTEM.ATOMTYPE[a];

		// Read the coordinates

		STREAM_PARSER >> SYSTEM.COORDS[a].X >> SYSTEM.COORDS[a].Y >> SYSTEM.COORDS[a].Z;
		
		// Wrap the coordinates
		SYSTEM.COORDS[a].X -= floor(SYSTEM.COORDS[a].X/SYSTEM.BOXDIM.X)*SYSTEM.BOXDIM.X;
		SYSTEM.COORDS[a].Y -= floor(SYSTEM.COORDS[a].Y/SYSTEM.BOXDIM.Y)*SYSTEM.BOXDIM.Y;
		SYSTEM.COORDS[a].Z -= floor(SYSTEM.COORDS[a].Z/SYSTEM.BOXDIM.Z)*SYSTEM.BOXDIM.Z;

		// Prepare velocities
		SYSTEM.VELOCITY[a].X = 0;
		SYSTEM.VELOCITY[a].Y = 0;
		SYSTEM.VELOCITY[a].Z = 0;
		
		// Prepare forces
		SYSTEM.FORCES[a].X = 0;
		SYSTEM.FORCES[a].Y = 0;
		SYSTEM.FORCES[a].Z = 0;
		
		// Prepare accelerations
		SYSTEM.ACCEL[a].X = 0;
		SYSTEM.ACCEL[a].Y = 0;
		SYSTEM.ACCEL[a].Z = 0;		
		
        if ( CONTROLS.COMPARE_FORCE ) // Reading positions from *.xyzf for force testing
		{	
			if(EXTENSION == "xyzf")	// Ignore forces in .xyzf file
			{
				if(a==0 && RANK==0)
				{
					cout << "	...Ignoring last three fields of atom info lines in xyzf file... " << endl;
					cout << "	...will read from specified force file instead." << endl;
				}

				//COORDFILE >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;			
			}
					
			// Read forces from separate force file
			CMPR_FORCEFILE >> SYSTEM.FORCES[a].X >> SYSTEM.FORCES[a].Y >> SYSTEM.FORCES[a].Z;
		}
		
        if(!CONTROLS.INIT_VEL) // Reading positions from *.xyz
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

			STREAM_PARSER >> SYSTEM.VELOCITY[a].X >> SYSTEM.VELOCITY[a].Y >> SYSTEM.VELOCITY[a].Z;
		}
		else
		{
			if(EXTENSION == "xyzf")
			{
				if(a==0 && RANK==0)
					cout << "	...Ignoring last three fields of atom info lines in xyzf file... " << endl;
				//COORDFILE >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;
			}
		}
		
		
		STREAM_PARSER.str("");
		STREAM_PARSER.clear();	
	}
	
    COORDFILE.close();
    COORDFILE.clear();
	
	if ( CONTROLS.COMPARE_FORCE ) 
		CMPR_FORCEFILE.close();
	
	if(RANK==0)
		cout << "   ...read complete" << endl << endl;	
	
	////////////////////////////////////////////////////////////
	// Figure out atom charges and masses, based on parameter 
	// file
	////////////////////////////////////////////////////////////
	 
	FF_SETUP_1: 
	
	if(RANK==0) 
		cout << "Reading atom info from parameter file..." << endl; 
	
	// Read in the possible atom types and thier features	
	
	PARAMFILE.open(CONTROLS.PARAM_FILE.data());
	if(!PARAMFILE.is_open())
	{
		cout << "ERROR: Cannot open coordinate file: " << CONTROLS.PARAM_FILE << endl;
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
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();
			
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
				TMP_ATOMTYPEIDX[i] = int(atof(TEMP_STR.data()));
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
	
	if (FF_PLOTS.N_PLOTS > 0)

		goto FF_SETUP_2; 
	

	// Assign atom features to atoms in SYSTEM data object, and the PAIR_FF object
	
	for(int a=0; a<SYSTEM.ATOMS;a++)
	{
		for(int i=0; i<NATMTYP; i++)
		{
			if(SYSTEM.ATOMTYPE[a] == TMP_ATOMTYPE[i])
			{
				
				SYSTEM.CHARGES[a]      = TMP_CHARGES[i];
				SYSTEM.MASS[a]         = TMP_MASS[i];
				SYSTEM.ATOMTYPE_IDX[a] = TMP_ATOMTYPEIDX[i];
				break;
			}			
		}
	}

	if(RANK==0)
		cout << "   ...read complete" << endl << endl;

	////////////////////////////////////////////////////////////
	// Explicitly account for layers.. 
	//
	// In contrast to the old approach, we only add images in (+) direction
	// So 0 layers returns the original cell (i.e. NATOMS atoms)
	//    1 layer  returns a 8*NATOMS atoms
	//    2 layers returns  27*NATOMS atoms, and so on.
	//
	// These atoms are explicitly added to the system, so its generally
	// a good idea to use neighbor lists to cut down on cost
	////////////////////////////////////////////////////////////	

	if(CONTROLS.N_LAYERS>0 )
	{	
		if(RANK == 0)
			cout << "Building " << CONTROLS.N_LAYERS << " layer(s)..." << endl;
	
		TEMP_IDX = SYSTEM.ATOMS;
		
		SYSTEM.PARENT   .resize(SYSTEM.ATOMS);
		SYSTEM.LAYER_IDX.resize(SYSTEM.ATOMS);
			
		// Create coordinates for the layer atoms. layer elements do not include 0, 0, 0, which is the main cell

		for(int a1=0; a1<SYSTEM.ATOMS; a1++)
		{			
			for(int n1=0; n1<=CONTROLS.N_LAYERS; n1++)
			{
				for(int n2=0; n2<=CONTROLS.N_LAYERS; n2++)
				{
					for(int n3=0; n3<=CONTROLS.N_LAYERS; n3++)
					{	
						if ((n1 == 0) && (n2 == 0) && (n3 == 0) )
						{
							SYSTEM.PARENT[a1] = -1;
							SYSTEM.LAYER_IDX[a1].X = n1;
							SYSTEM.LAYER_IDX[a1].Y = n2;
							SYSTEM.LAYER_IDX[a1].Z = n3;
						}
						else
						{											
							TEMP_XYZ.X = SYSTEM.COORDS.at(a1).X + n1 * SYSTEM.BOXDIM.X;
							TEMP_XYZ.Y = SYSTEM.COORDS.at(a1).Y + n2 * SYSTEM.BOXDIM.Y;
							TEMP_XYZ.Z = SYSTEM.COORDS.at(a1).Z + n3 * SYSTEM.BOXDIM.Z;
							
							TEMP_LAYER.X = n1;
							TEMP_LAYER.Y = n2;
							TEMP_LAYER.Z = n3;

							SYSTEM.COORDS       .push_back(TEMP_XYZ);
							SYSTEM.LAYER_IDX    .push_back(TEMP_LAYER);
							SYSTEM.ATOMTYPE     .push_back(SYSTEM.ATOMTYPE    .at(a1));
							SYSTEM.ATOMTYPE_IDX .push_back(SYSTEM.ATOMTYPE_IDX.at(a1));
							SYSTEM.CHARGES      .push_back(SYSTEM.CHARGES     .at(a1));
							SYSTEM.MASS         .push_back(SYSTEM.MASS        .at(a1));	
							SYSTEM.VELOCITY     .push_back(SYSTEM.VELOCITY    .at(a1));
							SYSTEM.PARENT       .push_back(a1);
							
							TEMP_IDX++;
						}
					}
				}
			}
		}
		
		SYSTEM.ATOMS = TEMP_IDX;

		SYSTEM.BOXDIM.X *= (CONTROLS.N_LAYERS + 1);
		SYSTEM.BOXDIM.Y *= (CONTROLS.N_LAYERS + 1);
		SYSTEM.BOXDIM.Z *= (CONTROLS.N_LAYERS + 1);

		SYSTEM.FORCES      .resize(SYSTEM.ATOMS);
		SYSTEM.ACCEL       .resize(SYSTEM.ATOMS);
		SYSTEM.VELOCITY_NEW.resize(SYSTEM.ATOMS);
		
		for(int a1=0; a1<SYSTEM.ATOMS; a1++)
		{
			SYSTEM.ACCEL[a1].X = 0;
			SYSTEM.ACCEL[a1].Y = 0;
			SYSTEM.ACCEL[a1].Z = 0;
		}
		
		if(!CONTROLS.COMPARE_FORCE)
		{
			for(int a1=0; a1<SYSTEM.ATOMS; a1++)
			{
			
				SYSTEM.FORCES[a1].X = 0;
				SYSTEM.FORCES[a1].Y = 0;
				SYSTEM.FORCES[a1].Z = 0;
			}
		}

		
		if(RANK == 0)
		{
			cout << "	New total atoms:    " << SYSTEM.ATOMS << endl;
			cout << "	New box dimensions: " << SYSTEM.BOXDIM.X << " " << SYSTEM.BOXDIM.Y << " " << SYSTEM.BOXDIM.Z << endl << endl;
		}	
	}
	
	////////////////////////////////////////////////////////////
	// Initialize velocities, if requested. Use the box Muller
	// method, setup up thermostat if requested
	////////////////////////////////////////////////////////////

	if(RANK==0)
		cout << "Setting up thermostats/velocities..." << endl; 

	// Set up the Nose-Hoover thermostat
	
	NOSE_HOOVER THERMOSTAT;
	
	if(CONTROLS.USE_HOOVER_THRMOSTAT)
	{
		THERMOSTAT.TIME  = CONTROLS.FREQ_UPDATE_THERMOSTAT;
		THERMOSTAT.VISCO = 0.0;
		THERMOSTAT.COORD = 0.0;
		THERMOSTAT.N_DOF = 3*SYSTEM.ATOMS - 3;
		THERMOSTAT.CHRG  = THERMOSTAT.N_DOF * Kb * CONTROLS.TEMPERATURE * THERMOSTAT.TIME * THERMOSTAT.TIME / Tfs / Tfs; // Need to convert from fs to sim units
		
	}
	else
	{
		THERMOSTAT.N_DOF = 3*SYSTEM.ATOMS - 3;		
	}
	
	
	// Use box Muller to initialize velocities

	Vol = SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;

	if ( CONTROLS.INIT_VEL ) 
	{
		// Declare some task-specific variables
		
		double x1, x2 , y1 ,y2;
		double sigma;
		int    counter = 0;
		
		srand(CONTROLS.SEED);
		
		for(int a=0; a<SYSTEM.ATOMS; a++)
		{
			
			sigma = sqrt(CONTROLS.TEMPERATURE * Kb / SYSTEM.MASS[a] );

			// Do for x...
			
			x1    = double(rand())/double(RAND_MAX);
			x2    = double(rand())/double(RAND_MAX);
			
			if ( (x1 < 0.0 || x1 > 1.0) || ( x2 < 0.0 || x2 > 1.0 ) )
			{
				cout << "Bad random variable" << endl;
				exit_run(1);
			}
			
			y1 = sqrt(-2.0 * log(x1)) * cos(2.0 * M_PI * x2);
			y2 = sqrt(-2.0 * log(x1)) * sin(2.0 * M_PI * x2);
			y1 *= sigma;
			y2 *= sigma;			

			if ( counter % 2 == 0 ) // Use either y2 or y1 here for maximum fun.
				SYSTEM.VELOCITY[a].X = y2; 
			else 
				SYSTEM.VELOCITY[a].X = y1;
			 
			counter++;
			 
			 
			// Do for y...
			
			x1    = double(rand())/double(RAND_MAX);
			x2    = double(rand())/double(RAND_MAX);
			
			if ( (x1 < 0.0 || x1 > 1.0) || ( x2 < 0.0 || x2 > 1.0 ) )
			{
				cout << "Bad random variable" << endl;
				exit_run(1);
			}
			
			y1 = sqrt(-2.0 * log(x1)) * cos(2.0 * M_PI * x2);
			y2 = sqrt(-2.0 * log(x1)) * sin(2.0 * M_PI * x2);
			y1 *= sigma;
			y2 *= sigma;			

			if ( counter % 2 == 0 ) // Use either y2 or y1 here for maximum fun.
				SYSTEM.VELOCITY[a].Y = y2; 
			else 
				SYSTEM.VELOCITY[a].Y = y1;
			 
			counter++;		
				
			 
			// Do for z...
			
			x1    = double(rand())/double(RAND_MAX);
			x2    = double(rand())/double(RAND_MAX);
			
			if ( (x1 < 0.0 || x1 > 1.0) || ( x2 < 0.0 || x2 > 1.0 ) )
			{
				cout << "Bad random variable" << endl;
				exit_run(1);
			}

			y1 = sqrt(-2.0 * log(x1)) * cos(2.0 * M_PI * x2);
			y2 = sqrt(-2.0 * log(x1)) * sin(2.0 * M_PI * x2);
			y1 *= sigma;
			y2 *= sigma;			

			if ( counter % 2 == 0 ) // Use either y2 or y1 here for maximum fun.
				SYSTEM.VELOCITY[a].Z = y2; 
			else 
				SYSTEM.VELOCITY[a].Z = y1;

			counter++;				 
		}		
	}
	
	if(RANK==0)
		cout << "   ...setup complete" << endl << endl;
	
	
	////////////////////////////////////////////////////////////
	// Setup and test system velocity center of mass... For
	// consevation of momentum, this number should zero (or
	// VERY close)
	////////////////////////////////////////////////////////////  
	// ... Also, figure out how many atoms of each type we have
	// ... this becomes useful for printing VASP POSCAR files
	// //////////////////////////////////////////////////////////// 

	if(RANK==0)
		cout << "Running velocity sanity checks..." << endl; 
  
	// Test/fix velocity center of mass here
	
	XYZ    TEMP_VEL;
	double TEMP_MASS;
	
	TEMP_VEL.X = 0;
	TEMP_VEL.Y = 0;
	TEMP_VEL.Z = 0;
	TEMP_MASS  = 0;

	for(int a=0; a<SYSTEM.ATOMS; a++)
    {
		TEMP_VEL.X += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].X;
		TEMP_VEL.Y += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].Y;
		TEMP_VEL.Z += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].Z;
		TEMP_MASS  += SYSTEM.MASS[a];
		
		for(int j=0; j<TMP_ATOMTYPE.size(); j++)
			if (SYSTEM.ATOMTYPE[a] == TMP_ATOMTYPE[j])
					TMP_NATOMTYPE[j]++;
    }
	
	// Check our velocity center of mass.. hopefully this is (or is very close to) zero
	
	if(RANK==0)  
		cout	<< "	Initial velocity center of mass: (" 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.X/TEMP_MASS << ", " 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Y/TEMP_MASS << ", " 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Z/TEMP_MASS << ") "
				<< endl;
  
	// In case it isn't, correct velocities to make it so: 
	 
	for(int a=0; a<SYSTEM.ATOMS; a++)
    {
		SYSTEM.VELOCITY[a].X -= TEMP_VEL.X/TEMP_MASS;
		SYSTEM.VELOCITY[a].Y -= TEMP_VEL.Y/TEMP_MASS;
		SYSTEM.VELOCITY[a].Z -= TEMP_VEL.Z/TEMP_MASS;
    }

	TEMP_VEL.X = 0;
	TEMP_VEL.Y = 0;
	TEMP_VEL.Z = 0;
	TEMP_MASS  = 0;

	// Now run a sanity check to make sure the new velocity center of mass is ~0
	
	for(int a=0; a<SYSTEM.ATOMS; a++)
    {
		TEMP_VEL.X += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].X;
		TEMP_VEL.Y += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].Y;
		TEMP_VEL.Z += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].Z;
		TEMP_MASS  += SYSTEM.MASS[a];
    }
	
	if(RANK==0)
		cout	<< "	Final   velocity center of mass: (" 
			<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.X/TEMP_MASS << ", " 
			<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Y/TEMP_MASS << ", " 
			<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Z/TEMP_MASS << ") "
			<< endl;
  
	
  	// Print some info on density, etc.

    dens_mol = (SYSTEM.ATOMS * 1.0e24) / (6.0221e23 * Vol);
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

	PARAMFILE.open(CONTROLS.PARAM_FILE.data());
	if(!PARAMFILE.is_open())
	{
		cout << "ERROR: Cannot open coordinate file: " << CONTROLS.PARAM_FILE << endl;
		exit_run(0);
	}
	
	if(RANK==0)
		cout << endl << "Reading force field parameters..." << endl;
	
	FOUND_END = false;

	string  TEMP_SEARCH_2B = "some default text";
	string  TEMP_SEARCH_3B = "some default text";
	string	TEMP_TYPE;
	int     NO_PAIRS, NO_TRIPS;
	int		TMP_TERMS1, TMP_TERMS2;
	double	TMP_LOW  = -1;
	double 	TMP_HIGH =  1;
	
	while (FOUND_END == false)
	{
		getline(PARAMFILE,LINE);

		// Break out of loop

		if(LINE.find("ENDFILE") != string::npos)
		{
			
			// Rewind so we can set the special 3-body cutoffs
			
			PARAMFILE.seekg(0);
			
			while (FOUND_END == false)		
			{
				getline(PARAMFILE,LINE);
				
				if(LINE.find("ENDFILE") != string::npos)
					break;		
				
				else if(LINE.find("SPECIAL 3B S_MINIM:") != string::npos)
				{
					STREAM_PARSER.str(LINE);
	
					STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TMP_TERMS1;
	
					#if VERBOSITY == 1
						cout << "	Note: Setting specific 3-body r_min values: " << endl;
					#endif	

					double TMP_VAL;
					
					string  TMP_IJ,  TMP_IK,  TMP_JK;	
					string TARG_IJ, TARG_IK, TARG_JK;
		
					for(int i=0; i<TMP_TERMS1; i++)
					{
						getline(PARAMFILE,LINE);
						
						STREAM_PARSER.str("");
						STREAM_PARSER.clear();
						
						STREAM_PARSER.str(LINE);
						STREAM_PARSER >> TEMP_STR;	// Just and index for the pair type
						STREAM_PARSER >> TEMP_STR;	// Which 3-body type is it?

						STREAM_PARSER >> TMP_IJ;	// What is the IJ?
						STREAM_PARSER >> TMP_IK;	// What is the IK?
						STREAM_PARSER >> TMP_JK;	// What is the JK?

						// Check that triplet pair types are correct
		
						try
						{
							TMP_IJ = FF_2BODY[ PAIR_MAP[ TMP_IJ ] ].PRPR_NM;
						}
						catch(...)
						{
							cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
							cout << "		Triplet type:              " << TEMP_STR << endl;
							cout << "		First distance, pair type: " << TMP_IJ << endl;
						}
						try
						{
							TMP_IK = FF_2BODY[ PAIR_MAP[ TMP_IK ] ].PRPR_NM;
						}
						catch(...)
						{
							cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
							cout << "		Triplet type:              " << TEMP_STR << endl;
							cout << "		First distance, pair type: " << TMP_IK << endl;
						}
						try
						{
							TMP_JK = FF_2BODY[ PAIR_MAP[ TMP_JK ] ].PRPR_NM;
						}
						catch(...)
						{
							cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
							cout << "		Triplet type:              " << TEMP_STR << endl;
							cout << "		First distance, pair type: " << TMP_JK << endl;
						}
		
						TARG_IJ = FF_2BODY[ PAIR_MAP[ FF_3BODY[TRIAD_MAP[TEMP_STR]].ATMPAIR1 ] ].PRPR_NM;
						TARG_IK = FF_2BODY[ PAIR_MAP[ FF_3BODY[TRIAD_MAP[TEMP_STR]].ATMPAIR2 ] ].PRPR_NM;
						TARG_JK = FF_2BODY[ PAIR_MAP[ FF_3BODY[TRIAD_MAP[TEMP_STR]].ATMPAIR3 ] ].PRPR_NM;
	
						// Read the first inner cutoff

						STREAM_PARSER >> TMP_VAL;
		
						if      ( (TMP_IJ == TARG_IJ) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X = TMP_VAL;
		
						else if ( (TMP_IJ == TARG_IK) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y = TMP_VAL;
		
						else if ( (TMP_IJ == TARG_JK) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z = TMP_VAL;


						// Read the second inner cutoff

						STREAM_PARSER >> TMP_VAL;
		
						if      ( (TMP_IK == TARG_IJ) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X = TMP_VAL;
		
						else if ( (TMP_IK == TARG_IK) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y = TMP_VAL;
		
						else if ( (TMP_IK == TARG_JK) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z = TMP_VAL;

		
						// Read the third inner cutoff

						STREAM_PARSER >> TMP_VAL;
		
						if      ( (TMP_JK == TARG_IJ) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X = TMP_VAL;
		
						else if ( (TMP_JK == TARG_IK) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y = TMP_VAL;
		
						else if ( (TMP_JK == TARG_JK) && (FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z == -1) )
							FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z = TMP_VAL;

						#if VERBOSITY == 1
							cout << "		" << TEMP_STR << "( " <<  TARG_IJ << ", " << TARG_IK << ", " << TARG_JK << "): " 
								              << FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X << ", "
										 	  << FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y << ", "
											  << FF_3BODY[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z << endl;
						#endif	
					}
					
					STREAM_PARSER.str("");
					STREAM_PARSER.clear();
	
				}
				
				else if(LINE.find("SPECIAL 3B S_MAXIM: ") != string::npos)
				{
					STREAM_PARSER.str(LINE);			
					STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TMP_TERMS1;
					STREAM_PARSER.str("");
					STREAM_PARSER.clear();	
			
					for(int i=0;i<TMP_TERMS1; i++)
					{
						getline(PARAMFILE,LINE);
	
						STREAM_PARSER.str(LINE);
						STREAM_PARSER >> TMP_TERMS2 >> TEMP_STR;	// We could use the integer, but this way is probably safer
						STREAM_PARSER >> FF_3BODY[TMP_TERMS2].S_MAXIM_3B;
						
						if(FF_3BODY[TMP_TERMS2].S_MAXIM_3B > NEIGHBOR_LIST.MAX_CUTOFF_3B)
							 NEIGHBOR_LIST.MAX_CUTOFF_3B = FF_3BODY[TMP_TERMS2].S_MAXIM_3B;
						
						STREAM_PARSER.str("");
						STREAM_PARSER.clear();	
					}

				}
				
				
			}
			
			FOUND_END = true;
			
			if(RANK==0)
			{
				cout << "   ...read complete." << endl << endl;
				cout << "Notes on simulation: " << endl;
				cout << "	Using fpenalty power " << FPENALTY_POWER << endl;
				
				
				
				int FOUND_SPECIAL = 0;
	
				for(int i=0; i<FF_3BODY.size(); i++)
					if(FF_3BODY[i].S_MAXIM_3B != -1)
						FOUND_SPECIAL++;
	
				if(FOUND_SPECIAL<FF_3BODY.size())
				{
					cout << "	Special 3-body cutoffs defined: " << endl;
		
					for(int i=0; i<FF_3BODY.size(); i++)
						if(FF_3BODY[i].S_MAXIM_3B != -1)
							cout << "		" <<  i << " " << TRIAD_MAP_REVERSE[i] << " " << FF_3BODY[i].S_MAXIM_3B << endl;
				}
	
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
			STREAM_PARSER.str(LINE);
			STREAM_PARSER >> TEMP_STR;
			STREAM_PARSER >> TEMP_STR;
			STREAM_PARSER.str(""); 
			STREAM_PARSER.clear();
				
			if (TEMP_STR=="true"  || TEMP_STR=="True"  || TEMP_STR=="TRUE"  || TEMP_STR == "T" || TEMP_STR == "t")
				CONTROLS.USE_COULOMB = true;
			else
				CONTROLS.USE_COULOMB = false;

			PARAMFILE >> TEMP_STR >> TEMP_STR;
			if (TEMP_STR=="true"  || TEMP_STR=="True"  || TEMP_STR=="TRUE"  || TEMP_STR == "T" || TEMP_STR == "t")
				CONTROLS.FIT_COUL = true;
			else
				CONTROLS.FIT_COUL = false;
			
			PARAMFILE >> TEMP_STR >> TEMP_STR;
			if (TEMP_STR=="true"  || TEMP_STR=="True"  || TEMP_STR=="TRUE"  || TEMP_STR == "T" || TEMP_STR == "t")
				CONTROLS.USE_OVERCOORD = true;
			else
				CONTROLS.USE_OVERCOORD = false;
			
			PARAMFILE >> TEMP_STR >> TEMP_STR;
			if (TEMP_STR=="true"  || TEMP_STR=="True"  || TEMP_STR=="TRUE"  || TEMP_STR == "T" || TEMP_STR == "t")
				CONTROLS.FIT_POVER = true;
			else
				CONTROLS.FIT_POVER = false;
			
			PARAMFILE >> TEMP_STR >> TEMP_STR;
			if (TEMP_STR=="true"  || TEMP_STR=="True"  || TEMP_STR=="TRUE"  || TEMP_STR == "T" || TEMP_STR == "t")
				CONTROLS.USE_3B_CHEBY = true;
			else
				CONTROLS.USE_3B_CHEBY = false;
			if(RANK==0)
			{
				cout << "		...Compute electrostatics?      " << boolalpha << CONTROLS.USE_COULOMB << endl;
				cout << "		...Use fit charges?             " << boolalpha << CONTROLS.FIT_COUL << endl;
				cout << "		...Compute ReaxFF overbonding?  " << boolalpha << CONTROLS.USE_OVERCOORD << endl;
				cout << "		...Use fit overbonding param?   " << boolalpha << CONTROLS.FIT_POVER << endl;
				cout << "		...Use 3-body Cheby params?     " << boolalpha << CONTROLS.USE_3B_CHEBY << endl;
			
				cout << "	...Read FF controls..." << endl;	
			}

			
			PARAMFILE.ignore();
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
				STREAM_PARSER >> TMP_TERMS1;
				STREAM_PARSER >> TMP_TERMS2;

				if (STREAM_PARSER >>  TMP_LOW)
				{
					if( TMP_LOW < -1.0 ||  TMP_LOW > +1.0 )
					{
						cout << "ERROR: CHEBY_RANGE_LOW must be betwee -1 and 1" << endl;
						exit_run(0);
					}
				}
				
				
				if (STREAM_PARSER >>  TMP_HIGH)
				{
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
					//FF_2BODY[i].PENALTY_SCALE    = 1.0e4;
					FF_2BODY[i].CUBIC_SCALE      = 1.0;
					FF_2BODY[i].CHEBY_RANGE_HIGH = TMP_HIGH;
					FF_2BODY[i].CHEBY_RANGE_LOW  = TMP_LOW;
				}

			}
			
			// Read in general pair parameters
			
			getline(PARAMFILE,LINE);
			getline(PARAMFILE,LINE);
			
			for(int i=0; i<NO_PAIRS; i++)
			{
				FF_2BODY[i].PAIRTYP = TEMP_TYPE;

				PARAMFILE >> FF_2BODY[i].PAIRIDX;
				PARAMFILE >> FF_2BODY[i].ATM1TYP;
				PARAMFILE >> FF_2BODY[i].ATM2TYP;	
				PARAMFILE >> FF_2BODY[i].S_MINIM;	
				PARAMFILE >> FF_2BODY[i].S_MAXIM;
				
				if(FF_2BODY[i].S_MAXIM > NEIGHBOR_LIST.MAX_CUTOFF)
				{
					 NEIGHBOR_LIST.MAX_CUTOFF    = FF_2BODY[i].S_MAXIM;
					 NEIGHBOR_LIST.MAX_CUTOFF_3B = FF_2BODY[i].S_MAXIM;
				}
				 	
				PARAMFILE >> FF_2BODY[i].S_DELTA;
				
				if((FF_PLOTS.N_PLOTS == 0) &&
				   (  FF_2BODY[i].S_MAXIM > 0.5* SYSTEM.BOXDIM.X
				   || FF_2BODY[i].S_MAXIM > 0.5* SYSTEM.BOXDIM.Y
				   || FF_2BODY[i].S_MAXIM > 0.5* SYSTEM.BOXDIM.Z ) )
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
					PARAMFILE >> FF_2BODY[i].CHEBY_TYPE;	// How does the user want distance transformed?

					if(FF_2BODY[i].CHEBY_TYPE == "MORSE")
						PARAMFILE >> FF_2BODY[i].LAMBDA;	
					
					FF_2BODY[i].SNUM          = TMP_TERMS1;
					FF_2BODY[i].SNUM_3B_CHEBY = TMP_TERMS2;
				}
				else if(TEMP_TYPE =="LJ")
				{
					FF_2BODY[i].SNUM = 2;
				}
				else // Splines type
				{
					FF_2BODY[i].SNUM = (2+floor((FF_2BODY[i].S_MAXIM - FF_2BODY[i].S_MINIM)/FF_2BODY[i].S_DELTA))*2;
				}				
			}
			
			if(CONTROLS.USE_OVERCOORD)
			{		
				getline(PARAMFILE,LINE);
				getline(PARAMFILE,LINE);
				getline(PARAMFILE,LINE);
								
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
				FF_2BODY[i].PENALTY_DIST = double(atof(TEMP_STR.data()));
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();	
		}
		
		else if(LINE.find("PAIR CHEBYSHEV PENALTY SCALING: ") != string::npos)
		{
			STREAM_PARSER.str(LINE);
			STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR;
			for(int i=0; i<NO_PAIRS; i++)
				FF_2BODY[i].PENALTY_SCALE = double(atof(TEMP_STR.data()));
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
			
			FF_3BODY.resize(NO_TRIPS);	
			
			for (int i=0; i<NO_TRIPS; i++)
			{	
				FF_3BODY[i].S_MINIM_3B.X = -1;
				FF_3BODY[i].S_MINIM_3B.Y = -1;
				FF_3BODY[i].S_MINIM_3B.Z = -1;
			}	

			TEMP_SEARCH_3B = "TRIPLET ";
			TEMP_SEARCH_3B.append(FF_2BODY[0].PAIRTYP); // Syntax ok b/c all pairs have same FF type, and 2b and 3b are same pair type
			TEMP_SEARCH_3B.append(" PARAMS");	
			
			if(RANK==0)
				cout << "	...Read FF triplet specifications..." << endl;
			
		}
		
		// Read in pair parameters
		
		else if(LINE.find(TEMP_SEARCH_2B) != string::npos) // "PAIR <PAIRTYPE> PARAMS"
		{	

			// Read in the specific atom pair parameters
			
			if(RANK==0)
				cout << "	...Reading all remaining force field parameters..." << endl;

			for(int i=0; i<NO_PAIRS; i++)
			{
				getline(PARAMFILE,LINE);	// Blank line
				getline(PARAMFILE,LINE);	// "PAIRTYPE PARAMS: <index> <atom1> <atom2>"	
				getline(PARAMFILE,LINE);	// Blank line
				
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
				for(int j=0; j<NATMTYP; j++)
				{
					for(int i=0; i<NO_PAIRS; i++)
					{
						if( (FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] && FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j]) || (FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j] && FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j]) )
						{
							TMP_CHARGES[j] = sqrt(fabs(FF_2BODY[i].PAIR_CHRG))*TMP_SIGN[j];
							FF_2BODY[i].ATM1CHG = TMP_CHARGES[j];
							FF_2BODY[i].ATM2CHG = TMP_CHARGES[j];
					
							break;
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

			
				if(RANK==0)
				{
					cout << "		Re-setting individual atom charges based on pair charges :" << endl;
				
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
			for(int i=0; i<NO_TRIPS; i++)
			{
				getline(PARAMFILE,LINE);	// Blank line
				getline(PARAMFILE,LINE);	// "TRIPLETYP PARAMS: <index> <pair1> <pair2> <pair3>"
				 
				STREAM_PARSER.str(LINE);
				
				STREAM_PARSER >> TEMP_STR  >> TEMP_STR;				
				STREAM_PARSER >> FF_3BODY[i].TRIPINDX;				
				STREAM_PARSER >> FF_3BODY[i].ATMPAIR1;
				STREAM_PARSER >> FF_3BODY[i].ATMPAIR2;
				STREAM_PARSER >> FF_3BODY[i].ATMPAIR3;	
				
				// Get rid of the colon at the end of the atom pair:
				
				FF_3BODY[i].ATMPAIR3 = FF_3BODY[i].ATMPAIR3.substr(0,FF_3BODY[i].ATMPAIR3.length()-1);
							
				STREAM_PARSER >> TEMP_STR;
				STREAM_PARSER >> FF_3BODY[i].N_TRUE_ALLOWED_POWERS;
				STREAM_PARSER >> TEMP_STR;
				STREAM_PARSER >> FF_3BODY[i].N_ALLOWED_POWERS;
				
				FF_3BODY[i].ALLOWED_POWERS.resize(FF_3BODY[i].N_ALLOWED_POWERS);
				FF_3BODY[i].EQUIV_INDICIES.resize(FF_3BODY[i].N_ALLOWED_POWERS);
				FF_3BODY[i].PARAM_INDICIES.resize(FF_3BODY[i].N_ALLOWED_POWERS);
				FF_3BODY[i].PARAMS        .resize(FF_3BODY[i].N_ALLOWED_POWERS);
				
				STREAM_PARSER.str("");
				STREAM_PARSER.clear();	

				getline(PARAMFILE,LINE);	// Blank line
				getline(PARAMFILE,LINE);	// header line

				
				getline(PARAMFILE,LINE);	// dashes line
				 
				for(int j=0; j<FF_3BODY[i].N_ALLOWED_POWERS; j++)
				{
					PARAMFILE >> TEMP_STR;
					PARAMFILE >> FF_3BODY[i].ALLOWED_POWERS[j].X;
					PARAMFILE >> FF_3BODY[i].ALLOWED_POWERS[j].Y;
					PARAMFILE >> FF_3BODY[i].ALLOWED_POWERS[j].Z;
					PARAMFILE >> FF_3BODY[i].EQUIV_INDICIES[j];
					PARAMFILE >> FF_3BODY[i].PARAM_INDICIES[j];
					PARAMFILE >> FF_3BODY[i].PARAMS        [j];
					PARAMFILE.ignore();

				} 
			}
			
			if (RANK==0)
				cout << "	...Read 3-body FF params..." << endl << endl;;				
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
		
		// Read the pair maps
		
		else if(LINE.find("TRIPMAPS: ") != string::npos)
		{
			STREAM_PARSER.str(LINE);			
			STREAM_PARSER >> TEMP_STR;
			STREAM_PARSER >> TMP_TERMS1;	
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();	
			
			if (RANK==0)
				cout << "	Reading  " << TMP_TERMS1 << " triplets for mapping" << endl;
			
			for(int i=0; i<TMP_TERMS1; i++)
			{
				PARAMFILE >> TMP_TERMS2;
				PARAMFILE >> TEMP_TYPE;
				
				if (RANK==0)
					cout << "	........Reading triplet: " << TEMP_TYPE << " with mapped index: " << TMP_TERMS2 << endl; 
				
				TRIAD_MAP.insert(make_pair(TEMP_TYPE,TMP_TERMS2));
				TRIAD_MAP_REVERSE.insert(make_pair(TMP_TERMS2,TEMP_TYPE));				
					
			}
			if (RANK==0)
				cout << "	...Read FF triadmaps..." << endl;					
		}	

	}	
	
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
	
	// Free up some memory.. swap the contents of currend vectors with a vector w/ no assigned mem

	vector<double>().swap(TMP_MASS);	
	
	////////////////////////////////////////////////////////////
	// Print a summary of the force field
	////////////////////////////////////////////////////////////  	
	
	if(RANK==0)
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
			if(FF_2BODY[i].CHEBY_TYPE == "MORSE")
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
				cout << FF_2BODY[i].CHEBY_TYPE << " ";
			
				if(FF_2BODY[i].CHEBY_TYPE == "MORSE")
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
			cout << "	Triplet type and polynomial order per triplet: " << endl;
			cout << "		" << FF_2BODY[0].PAIRTYP << " " << FF_2BODY[0].SNUM_3B_CHEBY << endl << endl;	
			
			cout << "	Interaction parameters for each triplet type: " << endl;

			for(int i=0;i<FF_3BODY.size(); i++)
			{
				cout << "	" << FF_3BODY[i].TRIPINDX << "  " << FF_3BODY[i].ATMPAIR1 << " " << FF_3BODY[i].ATMPAIR2 << " " << FF_3BODY[i].ATMPAIR3 << ": ";
				cout << FF_3BODY[i].N_TRUE_ALLOWED_POWERS << " parameters, " << FF_3BODY[i].N_ALLOWED_POWERS << " total parameters "<< endl;	
				cout << "	     index  |  powers  |  equiv index  |  param index  | parameter " << endl;
				cout << "	   --------------------------------------------------------------------" << endl;	

				for(int j=0; j<FF_3BODY[i].ALLOWED_POWERS.size(); j++)
				{
					cout << "	      " << setw(6) << fixed << left << j << " ";
					cout << " " << setw(2) << fixed << left << FF_3BODY[i].ALLOWED_POWERS[j].X  << " ";
					cout << " " << setw(2) << fixed << left << FF_3BODY[i].ALLOWED_POWERS[j].Y  << " ";
					cout << " " << setw(2) << fixed << left << FF_3BODY[i].ALLOWED_POWERS[j].Z  << " ";
					cout << "       " << setw(8) << FF_3BODY[i].EQUIV_INDICIES[j] << " ";
					cout << "       " << setw(8) << FF_3BODY[i].PARAM_INDICIES[j] << " "; 
					cout << "       " << setw(8) << FF_3BODY[i].PARAMS[j] << endl; 
	
				}

				cout << endl;
			}	 	
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
	
	////////////////////////////////////////////////////////////
	// Set up the neighbor list
	////////////////////////////////////////////////////////////

	if(NEIGHBOR_LIST.USE)
	{
		if(RANK == 0)
			cout << "Initializing the neighbor list..." << endl;
		
		NEIGHBOR_LIST.INITIALIZE(SYSTEM);
		NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);
	}
	

	////////////////////////////////////////////////////////////
	// If PES printing is requested, do so here and then 
	// exit the program
	////////////////////////////////////////////////////////////  	
	
	// For now, just do this in serial
	
	if(RANK==0)
	{
		if (FF_PLOTS.N_PLOTS > 0)
		{
			if(FF_PLOTS.INCLUDE_2B)
			{
				if (isatty(fileno(stdout)))
					cout << COUT_STYLE.BOLD << COUT_STYLE.MAGENTA << "WARNING: 2B+3B only used for 2D 3B scans, not 4D 3B scans." << COUT_STYLE.ENDSTYLE << endl;
				else
					cout << "WARNING: 2B+3B only used for 2D 3B scans, not 4D 3B scans." << endl;
			}
			
			int ij, ik, jk;
			string ATM_TYP_1, ATM_TYP_2, ATM_TYP_3;
//			int ADD_2B_TO_3B_IDX = 0;
		
//			int IJ_STEP = 0;
//			int IK_STEP = 0;
//			int JK_STEP = 0;
		
			ifstream FULL_INFILE_3B;
			ifstream SCAN_INFILE_3B;
			ifstream SCAN_INFILE_2B;
		
			ofstream FULL_OUTFILE_3B;
			ofstream SCAN_OUTFILE_3B;
			ofstream SCAN_OUTFILE_2B;			

			XYZ 	TMP_DISTS;
			double 	TMP_PES;
		
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
					ij =  PAIR_MAP[FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATMPAIR1];
					ik =  PAIR_MAP[FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATMPAIR2];
					jk =  PAIR_MAP[FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATMPAIR3];
			
					ATM_TYP_1 = FF_2BODY[ij].ATM1TYP;
					ATM_TYP_2 = FF_2BODY[jk].ATM1TYP;
					ATM_TYP_3 = FF_2BODY[jk].ATM2TYP;
				
					cout << "	Will work with pair types: " << FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATMPAIR1 << " " << FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATMPAIR2 << " " << FF_3BODY[FF_PLOTS.TYPE_INDEX[i]].ATMPAIR3 << endl;
					cout << "	and atom types:            " << ATM_TYP_1 << " " << ATM_TYP_2 << " " << ATM_TYP_3 << endl;
				
				
					if(FF_PLOTS.DO_4D)
						Print_3B_Cheby(CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP, ATM_TYP_1, ATM_TYP_2, ATM_TYP_3, ij, ik, jk);

					cout << "	Now printing for n scans: " << FF_PLOTS.N_SCAN << endl;
				
					for(int j=0; j<FF_PLOTS.N_SCAN; j++)
					{
						IJ_DIST_3B.clear();
						IK_DIST_3B.clear();
						JK_DIST_3B.clear();
						PES_VAL_3B.clear();	
					
						if(FF_PLOTS.PARENT_TYPE[j] == FF_PLOTS.TYPE_INDEX[i] )
						{						
							Print_3B_Cheby_Scan(CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP, ATM_TYP_1, ATM_TYP_2, ATM_TYP_3, ij, ik, jk, FF_PLOTS, j);	
						
							// If 2+3b interactions requested, include them
							if(FF_PLOTS.INCLUDE_2B)	
							{
								// We don't want to double-count ewald interactions, so turn off charges for this part (they've already been calcuated in the 3B scan)
							
								for(int a=0; a<FF_2BODY.size(); a++)
								{
										FF_2BODY[a].ATM1CHG = 0;
										FF_2BODY[a].ATM2CHG = 0;	
								}
							
								// Read in the printed 3b interactions
							
								SCAN_INFILE_3B.open(SCAN_FILE_3B.data());
					
								if(!SCAN_INFILE_3B.is_open())
								{
									cout << "ERROR-1: Cannot open file " << SCAN_FILE_3B << " for 2B addition." << endl;
									exit_run(0);
								}
					
								cout << "	Reading in 3B scans from file: "	<< SCAN_FILE_3B << endl;						
					
								while(!SCAN_INFILE_3B.eof())
								{
									SCAN_INFILE_3B >> TMP_DISTS.X;
									SCAN_INFILE_3B >> TMP_DISTS.Y;
									SCAN_INFILE_3B >> TMP_DISTS.Z;
									SCAN_INFILE_3B >> TMP_PES;
						
									IJ_DIST_3B.push_back(TMP_DISTS.X);
									IK_DIST_3B.push_back(TMP_DISTS.Y);
									JK_DIST_3B.push_back(TMP_DISTS.Z);
						
									PES_VAL_3B.push_back(TMP_PES);						
								}
					
								SCAN_INFILE_3B.close();

								// Compute/save/read the corresponding 2b interactions -- IJ
						
								Print_Cheby(FF_2BODY, FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].X, PAIR_MAP_REVERSE[FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].X], FF_PLOTS.INCLUDE_FCUT, FF_PLOTS.INCLUDE_CHARGES, "for_3b");						
								cout << "		Reading in 2B scans from file (1): "	<< SCAN_FILE_2B << endl;					
						
								SCAN_INFILE_2B.open(SCAN_FILE_2B.data());

								if(!SCAN_INFILE_2B.is_open())
								{
									cout << "ERROR-2: Cannot open file " << SCAN_FILE_2B << " for 2B addition." << endl;
									exit_run(0);
								}								
							
								while(!SCAN_INFILE_2B.eof())
								{
									SCAN_INFILE_2B >> TMP_DISTS.X;
									SCAN_INFILE_2B >> TMP_PES;
								
									IJ_DIST_2B.push_back(TMP_DISTS.X);
									PES_VAL_2B_IJ.push_back(TMP_PES);							
								}
							
								SCAN_INFILE_2B.close();								
							
								// Compute/save/read the corresponding 2b interactions -- IK
							
								Print_Cheby(FF_2BODY,  FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].Y, PAIR_MAP_REVERSE[FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].Y], FF_PLOTS.INCLUDE_FCUT, FF_PLOTS.INCLUDE_CHARGES, "for_3b");
							
								cout << "		Reading in 2B scans from file (2): "	<< SCAN_FILE_2B << endl;					
						
								SCAN_INFILE_2B.open(SCAN_FILE_2B.data());

								if(!SCAN_INFILE_2B.is_open())
								{
									cout << "ERROR-2: Cannot open file " << SCAN_FILE_2B << " for 2B addition." << endl;
									exit_run(0);
								}								
							
								while(!SCAN_INFILE_2B.eof())
								{
									SCAN_INFILE_2B >> TMP_DISTS.X;
									SCAN_INFILE_2B >> TMP_PES;
								
									IK_DIST_2B.push_back(TMP_DISTS.X);
									PES_VAL_2B_IK.push_back(TMP_PES);							
								}
							
								SCAN_INFILE_2B.close();								
							
								// Compute/save/read the corresponding 2b interactions -- JK

								Print_Cheby(FF_2BODY,  FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].Z, PAIR_MAP_REVERSE[FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].Z], FF_PLOTS.INCLUDE_FCUT, FF_PLOTS.INCLUDE_CHARGES, "for_3b");
							
								cout << "		Reading in 2B scans from file (3): "	<< SCAN_FILE_2B << endl;					
						
								SCAN_INFILE_2B.open(SCAN_FILE_2B.data());

								if(!SCAN_INFILE_2B.is_open())
								{
									cout << "ERROR-2: Cannot open file " << SCAN_FILE_2B << " for 2B addition." << endl;
									exit_run(0);
								}								
								
								while(!SCAN_INFILE_2B.eof())
								{
									SCAN_INFILE_2B >> TMP_DISTS.X;
									SCAN_INFILE_2B >> TMP_PES;
								
									JK_DIST_2B.push_back(TMP_DISTS.X);
									PES_VAL_2B_JK.push_back(TMP_PES);						
								}
							
								SCAN_INFILE_2B.close();		

								// Now add the 2B values to the 3B PES
					
					
								cout << PES_VAL_2B_IJ.size() << endl;
								cout << PES_VAL_2B_IK.size() << endl;
								cout << PES_VAL_2B_JK.size() << endl;
					
								for(int a=0; a<PES_VAL_3B.size(); a++)
								{
									// Check if ij distances the same between the 2b ij type and the 3b scan
									
									for(int b=0; b<PES_VAL_2B_IJ.size(); b++)
									{
										if(IJ_DIST_2B[b] == IJ_DIST_3B[a])
										{					
											PES_VAL_3B[a] += PES_VAL_2B_IJ[b];
											break;
										}
									}
									
									// Check if ik distances the same between the 2b ik type and the 3b scan
									
									for(int c=0; c<PES_VAL_2B_IK.size(); c++)
									{
										if(IK_DIST_2B[c] == IK_DIST_3B[a])
										{							
											PES_VAL_3B[a] += PES_VAL_2B_IK[c];
											break;
										}
									}	
									
									// Check if jk distances the same between the 2b jk type and the 3b scan
									
									for(int d=0; d<PES_VAL_2B_JK.size(); d++)
									{
										if(JK_DIST_2B[d] == JK_DIST_3B[a])
										{							
											PES_VAL_3B[a] += PES_VAL_2B_JK[d];
											break;
										}							
									}												
								}	

								// Finally, re-print the 3B PES file
					
								SCAN_FILE_3B.append("+2b.dat");
					
								SCAN_OUTFILE_3B.open(SCAN_FILE_3B.data());
					
								if(!SCAN_OUTFILE_3B.is_open())
								{
									cout << "ERROR-3: Cannot open file " << SCAN_FILE_3B << " for 2B addition." << endl;
									exit_run(0);
								}					
					
								cout << "		...Printing 2B+3B scan to file: " << SCAN_FILE_3B << endl;
								
								for(int a=0; a<PES_VAL_3B.size()-1; a++)
									SCAN_OUTFILE_3B << IJ_DIST_3B[a] << " " << IK_DIST_3B[a] << " " << JK_DIST_3B[a] << " " << PES_VAL_3B[a] << endl;		
					
								SCAN_OUTFILE_3B.close();							
							
								// Add the charges back in for the next round of 3B calculations
							
								for(int a=0; a<FF_2BODY.size(); a++)
								{
									for(int x=0; x<NATMTYP; x++)
									{
										if(FF_2BODY[a].ATM1TYP == TMP_ATOMTYPE[x])
											FF_2BODY[a].ATM1CHG = TMP_CHARGES[x];
		
										if(FF_2BODY[a].ATM2TYP == TMP_ATOMTYPE[x])
											FF_2BODY[a].ATM2CHG = TMP_CHARGES[x];
									}		
								}
							
							}
						}	
					}	
				
					scan_2b_idx++;				
				}
				else if(FF_PLOTS.NBODY[i] == 2)
				{

					cout << "	Will work with pair types: " << PAIR_MAP_REVERSE[FF_PLOTS.TYPE_INDEX[i]] << endl;
					cout << "	and atom types:            " << FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM2TYP << endl;
					
					Print_Cheby(FF_2BODY, FF_PLOTS.TYPE_INDEX[i], PAIR_MAP_REVERSE[FF_PLOTS.TYPE_INDEX[i]], FF_PLOTS.INCLUDE_FCUT, FF_PLOTS.INCLUDE_CHARGES, "");
				}
				else
				{
					cout << "Functionality not programmed yet." << endl;
				}

				cout << "	...Force field PES printing complete." << endl;
		
			idx_3b++;
			
			}
		
			exit_run(0);
		}
	}

  	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	//
	//					START THE SIMULATION
	//
	////////////////////////////////////////////////////////////  
	////////////////////////////////////////////////////////////

	STATISTICS.open("md_statistics.out");

	if (RANK==0)
		cout << "BEGIN SIMULATION:" << endl;
	
	double ke; // A temporary variable used when updating thermostatting 
	
	for(CONTROLS.STEP=0;CONTROLS.STEP<CONTROLS.N_MD_STEPS;CONTROLS.STEP++)	//start Big Loop here.
    {

	  	////////////////////////////////////////////////////////////
		// Do first half of coordinate/velocity updating
		////////////////////////////////////////////////////////////
		
		if(CONTROLS.STEP>0)	
		{
			if (RANK==0)
			{
				// Update coords

				for(int a1=0;a1<SYSTEM.ATOMS;a1++)
				{			
					SYSTEM.COORDS[a1].X += SYSTEM.VELOCITY[a1].X * CONTROLS.DELTA_T + 0.5*SYSTEM.ACCEL[a1].X * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
					SYSTEM.COORDS[a1].Y += SYSTEM.VELOCITY[a1].Y * CONTROLS.DELTA_T + 0.5*SYSTEM.ACCEL[a1].Y * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
					SYSTEM.COORDS[a1].Z += SYSTEM.VELOCITY[a1].Z * CONTROLS.DELTA_T + 0.5*SYSTEM.ACCEL[a1].Z * CONTROLS.DELTA_T * CONTROLS.DELTA_T;

				}

				// Apply thermostatting to coords
			
				if ( CONTROLS.USE_HOOVER_THRMOSTAT ) 
				{
					for(int a1=0;a1<SYSTEM.ATOMS;a1++)
					{					
						SYSTEM.COORDS[a1].X -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].X * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
						SYSTEM.COORDS[a1].Y -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].Y * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
						SYSTEM.COORDS[a1].Z -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].Z * CONTROLS.DELTA_T * CONTROLS.DELTA_T;					
					}				

					ke = kinetic_energy(SYSTEM);
				
					THERMOSTAT.VISCO_DOT0 = ( 2.0 * ke -  THERMOSTAT.N_DOF * Kb * CONTROLS.TEMPERATURE ) / ( THERMOSTAT.CHRG);
					THERMOSTAT.COORD     += THERMOSTAT.VISCO * CONTROLS.DELTA_T + 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * THERMOSTAT.VISCO_DOT0;
				}

				// Wrap the coordinates:
			
				if(CONTROLS.WRAP_COORDS)
				{
					for(int a1=0;a1<SYSTEM.ATOMS;a1++)
					{
						SYSTEM.COORDS[a1].X -= floor(SYSTEM.COORDS[a1].X / SYSTEM.BOXDIM.X) * SYSTEM.BOXDIM.X;
						SYSTEM.COORDS[a1].Y -= floor(SYSTEM.COORDS[a1].Y / SYSTEM.BOXDIM.Y) * SYSTEM.BOXDIM.Y;
						SYSTEM.COORDS[a1].Z -= floor(SYSTEM.COORDS[a1].Z / SYSTEM.BOXDIM.Z) * SYSTEM.BOXDIM.Z;
					}	
				}
			
				// Update first half of velocity (i.e. using a(t))... Also update maximum
				// velocity for neighbor lists:
			
				for(int a1=0;a1<SYSTEM.ATOMS;a1++)
				{
					SYSTEM.VELOCITY[a1].X += 0.5 * SYSTEM.ACCEL[a1].X * CONTROLS.DELTA_T;
					SYSTEM.VELOCITY[a1].Y += 0.5 * SYSTEM.ACCEL[a1].Y * CONTROLS.DELTA_T;
					SYSTEM.VELOCITY[a1].Z += 0.5 * SYSTEM.ACCEL[a1].Z * CONTROLS.DELTA_T;
					
					if(NEIGHBOR_LIST.USE)
					{
					
						NEIGHBOR_LIST.CURR_VEL 
						      = sqrt(SYSTEM.VELOCITY[a1].X*SYSTEM.VELOCITY[a1].X 
						           + SYSTEM.VELOCITY[a1].Y*SYSTEM.VELOCITY[a1].Y 
				      	           + SYSTEM.VELOCITY[a1].Z*SYSTEM.VELOCITY[a1].Z);
					
						if(NEIGHBOR_LIST.CURR_VEL > NEIGHBOR_LIST.MAX_VEL)
							NEIGHBOR_LIST.MAX_VEL = NEIGHBOR_LIST.CURR_VEL;
					}
					
				
				}
			
				// Apply thermostatting to velocities... Also update maximum
				// velocity for neighbor lists:

				if ( CONTROLS.USE_HOOVER_THRMOSTAT ) 
				{
					for(int a1=0;a1<SYSTEM.ATOMS;a1++)
					{
						SYSTEM.VELOCITY[a1].X -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].X * CONTROLS.DELTA_T;
						SYSTEM.VELOCITY[a1].Y -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].Y * CONTROLS.DELTA_T;
						SYSTEM.VELOCITY[a1].Z -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].Z * CONTROLS.DELTA_T;
						
						if(NEIGHBOR_LIST.USE)
						{
						
							NEIGHBOR_LIST.CURR_VEL 
							      = sqrt(SYSTEM.VELOCITY[a1].X*SYSTEM.VELOCITY[a1].X 
							           + SYSTEM.VELOCITY[a1].Y*SYSTEM.VELOCITY[a1].Y 
					      	           + SYSTEM.VELOCITY[a1].Z*SYSTEM.VELOCITY[a1].Z);
					
							if(NEIGHBOR_LIST.CURR_VEL > NEIGHBOR_LIST.MAX_VEL)
								NEIGHBOR_LIST.MAX_VEL = NEIGHBOR_LIST.CURR_VEL;
						}
					}
				}			
			}
		}
		
		cout.precision(14);			// Set output precision
    
      	////////////////////////////////////////////////////////////
		// Calculate acceleration
		////////////////////////////////////////////////////////////


		// FOR MPI:		Broadcast the position and optionally velocity from the root to all other processes. -- WHY ISN'T THIS IS A PRECOMPILER DIRECTIVE?
		sync_position(SYSTEM.COORDS, NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ATOMS, true);
		
		
		////////////////////////////////////////////////////////////
		// Check if neighbor list update needed/do updating
		////////////////////////////////////////////////////////////
		 
		if(CONTROLS.STEP>0)	
		{
			if(NEIGHBOR_LIST.USE)// && (RANK==0))
			{
				#ifdef USE_MPI
					MPI_Bcast(&NEIGHBOR_LIST.MAX_VEL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
				#endif
			
				NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);
			}
		
		}

		// this function calculates the spline and electrostatic forces. --- this should probably be removed
		ZCalc(SYSTEM, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP, NEIGHBOR_LIST);

		// FOR MPI:		Synchronize forces, energy, and pressure.
		
		#ifdef USE_MPI
			sum_forces(SYSTEM.ACCEL, SYSTEM.ATOMS, SYSTEM.TOT_POT_ENER, SYSTEM.PRESSURE_XYZ);
		#endif

		////////////////////////////////////////////////////////////
		// Print some info on new forces, compare to input forces if
		// requested
		////////////////////////////////////////////////////////////


	  	if ( CONTROLS.PRINT_FORCE && (CONTROLS.STEP+1)%CONTROLS.FREQ_FORCE == 0 && RANK == 0 ) 
	  	{
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
			{
				OUT_FORCEFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].X << endl;
				OUT_FORCEFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].Y << endl;
				OUT_FORCEFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].Z << endl;
				
				OUT_FORCELABL <<  SYSTEM.ATOMTYPE[a1] << " " << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].X << endl;
				OUT_FORCELABL <<  SYSTEM.ATOMTYPE[a1] << " " << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].Y << endl;
				OUT_FORCELABL <<  SYSTEM.ATOMTYPE[a1] << " " << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].Z << endl;
			}
		}

		if ( CONTROLS.COMPARE_FORCE && RANK == 0 ) 
		{
			
			int END = SYSTEM.ATOMS;
			
			if(CONTROLS.N_LAYERS>0)
				END = SYSTEM.ATOMS/pow(CONTROLS.N_LAYERS+1,3.0);
		

			// Check against read-in forces for code verification... Note, ferr is initialized to zero.
		
			for(int a1=0;a1<END;a1++)
			{
				ferr += (SYSTEM.ACCEL[a1].X - SYSTEM.FORCES[a1].X) * (SYSTEM.ACCEL[a1].X - SYSTEM.FORCES[a1].X);
				ferr += (SYSTEM.ACCEL[a1].Y - SYSTEM.FORCES[a1].Y) * (SYSTEM.ACCEL[a1].Y - SYSTEM.FORCES[a1].Y);
				ferr += (SYSTEM.ACCEL[a1].Z - SYSTEM.FORCES[a1].Z) * (SYSTEM.ACCEL[a1].Z - SYSTEM.FORCES[a1].Z);
			}
			
			ferr = sqrt(ferr/3.0/END);
			cout << "RMS force error = " << fixed << setprecision(6) << ferr << endl;
			
			return 0;
		}
		else if (RANK == 0)	// Print main simulation header 
		{
		  	////////////////////////////////////////////////////////////
			// Print the header for mains simulation output
			////////////////////////////////////////////////////////////		

			if ( CONTROLS.STEP == 0 ) 
			{
				printf("%8s %9s %15s %15s %15s %15s %15s", "Step", "Time", "Ktot/N", "Vtot/N", "Etot/N", "T", "P");
				
				STATISTICS << "# Step	Time	Ktot/N	Vtot/N	Etot/N	T	P";
	  
				if ( CONTROLS.USE_HOOVER_THRMOSTAT) 
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
				if ( CONTROLS.USE_HOOVER_THRMOSTAT ) 
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
		// Do some thermostatting and statistics updating/output
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

			if(CONTROLS.STEP>0)
			{
				//update second half of velocity
			
				if ( !CONTROLS.USE_HOOVER_THRMOSTAT ) 
				{
					// Update second half of velocity (i.e. using a(t+dt)):
				 
					for(int a1=0;a1<SYSTEM.ATOMS;a1++)
					{
						SYSTEM.VELOCITY[a1].X += 0.5*SYSTEM.ACCEL[a1].X * CONTROLS.DELTA_T;
						SYSTEM.VELOCITY[a1].Y += 0.5*SYSTEM.ACCEL[a1].Y * CONTROLS.DELTA_T;
						SYSTEM.VELOCITY[a1].Z += 0.5*SYSTEM.ACCEL[a1].Z * CONTROLS.DELTA_T;
					}
				} 
				else 
				{
					ke = kinetic_energy(SYSTEM);
				
					THERMOSTAT.VISCO_DOT1 = ( 2.0 * ke - THERMOSTAT.N_DOF * Kb * CONTROLS.TEMPERATURE ) / ( THERMOSTAT.CHRG);
					THERMOSTAT.VISCO1     = 0.0;

					// Iterative determination of velocities
					// See Martyna, Tobias, Klein JCP 101, 4177(1994) Appendix D.
				
					for ( int itr = 0; itr < 10; itr++ ) 
					{
						THERMOSTAT.VISCO1 = THERMOSTAT.VISCO + (THERMOSTAT.VISCO_DOT0 + THERMOSTAT.VISCO_DOT1) * 0.5 * CONTROLS.DELTA_T;
						THERMOSTAT.VSCALEH = 1.0 + 0.5 * CONTROLS.DELTA_T * THERMOSTAT.VISCO1;
					
						for(int a1=0;a1<SYSTEM.ATOMS;a1++)
						{
							SYSTEM.VELOCITY_NEW[a1].X = (SYSTEM.VELOCITY[a1].X + 0.5*SYSTEM.ACCEL[a1].X*CONTROLS.DELTA_T) / THERMOSTAT.VSCALEH;
							SYSTEM.VELOCITY_NEW[a1].Y = (SYSTEM.VELOCITY[a1].Y + 0.5*SYSTEM.ACCEL[a1].Y*CONTROLS.DELTA_T) / THERMOSTAT.VSCALEH;
							SYSTEM.VELOCITY_NEW[a1].Z = (SYSTEM.VELOCITY[a1].Z + 0.5*SYSTEM.ACCEL[a1].Z*CONTROLS.DELTA_T) / THERMOSTAT.VSCALEH;
						}

						ke = kinetic_energy(SYSTEM,"NEW");
						THERMOSTAT.VISCO_DOT1 = ( 2.0 * ke - THERMOSTAT.N_DOF * Kb * CONTROLS.TEMPERATURE ) / ( THERMOSTAT.CHRG);

					}
				
					for(int a1=0;a1<SYSTEM.ATOMS;a1++)
					{
						SYSTEM.VELOCITY[a1].X = SYSTEM.VELOCITY_NEW[a1].X;
						SYSTEM.VELOCITY[a1].Y = SYSTEM.VELOCITY_NEW[a1].Y;
						SYSTEM.VELOCITY[a1].Z = SYSTEM.VELOCITY_NEW[a1].Z;
					}
		  
					THERMOSTAT.VISCO = THERMOSTAT.VISCO1;
				}
			}

		  	////////////////////////////////////////////////////////////
			// Store statistics on the average simulation temperature
			// to allow for thermostat verification
			////////////////////////////////////////////////////////////

			Ktot = kinetic_energy(SYSTEM);	//calculate kinetic energy for scaling:
		
			SYSTEM.TEMPERATURE = 2.0 * Ktot / (THERMOSTAT.N_DOF * Kb);

			temp_sum += SYSTEM.TEMPERATURE;
		
			// Exit with an error if the set and block temperatures differ by more than CONTROLS.NVT_CONV_CUT
		
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
		
		  	////////////////////////////////////////////////////////////
			// If requested, compute pressure numerically, and accumulate
			// statistics
			////////////////////////////////////////////////////////////

			if ( CONTROLS.USE_NUMERICAL_PRESS ) 
			{
				if((CONTROLS.STEP == 0) && (NPROCS>1))
				{
					cout << "WARNING: Numerical pressure computed using only one processor...";
					cout << "This has potential to be a bottleneck!" << endl;
				}
				
				SYSTEM.PRESSURE_XYZ = numerical_pressure(SYSTEM, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP, NEIGHBOR_LIST);
			}

			SYSTEM.PRESSURE = (SYSTEM.PRESSURE_XYZ + 2.0 * Ktot / (3.0 * Vol)) * GPa;	// GPa = Unit conversion factor to GPa.
			press_sum += SYSTEM.PRESSURE;
		
		  	////////////////////////////////////////////////////////////
			// Periodically print simulation output
			////////////////////////////////////////////////////////////
		
			if ( (CONTROLS.STEP+1) % CONTROLS.FREQ_ENER == 0 && RANK == 0 ) 
			{
			
				printf("%8d %9.2f %15.7f %15.7f %15.7f %15.1f %15.3f", CONTROLS.STEP+1, (CONTROLS.STEP+1)*CONTROLS.DELTA_T_FS, Ktot/SYSTEM.ATOMS,SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS,(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS,SYSTEM.TEMPERATURE, SYSTEM.PRESSURE);
				STATISTICS << CONTROLS.STEP+1<< "	" << (CONTROLS.STEP+1)*CONTROLS.DELTA_T_FS<< "	" << Ktot/SYSTEM.ATOMS<< "	" <<SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS<< "	" <<(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS<< "	" <<SYSTEM.TEMPERATURE<< "	" << SYSTEM.PRESSURE;

				if ( CONTROLS.USE_HOOVER_THRMOSTAT ) 
				{
					printf("%15.7f\n", (Ktot + SYSTEM.TOT_POT_ENER + 0.5 * THERMOSTAT.VISCO * THERMOSTAT.VISCO * THERMOSTAT.CHRG + THERMOSTAT.N_DOF * Kb * CONTROLS.TEMPERATURE * THERMOSTAT.COORD) / SYSTEM.ATOMS);
					STATISTICS << "	" << (Ktot + SYSTEM.TOT_POT_ENER + 0.5 * THERMOSTAT.VISCO * THERMOSTAT.VISCO * THERMOSTAT.CHRG + THERMOSTAT.N_DOF * Kb * CONTROLS.TEMPERATURE * THERMOSTAT.COORD) / SYSTEM.ATOMS << endl;
				}
				else 
				{
					printf("%15.7f\n\n", (Ktot + SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS);
					STATISTICS << "	" << (Ktot + SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS << endl;
				}

				cout.flush();
			}	
		
		  	////////////////////////////////////////////////////////////
			// If requested, scale the velocities
			////////////////////////////////////////////////////////////	

			if ( (!CONTROLS.USE_HOOVER_THRMOSTAT) && ((CONTROLS.STEP+1) % int(CONTROLS.FREQ_UPDATE_THERMOSTAT) == 0) )
			{
				avg_temp +=  SYSTEM.TEMPERATURE;
				avg_temp /= int(CONTROLS.FREQ_UPDATE_THERMOSTAT);
				vscale    = sqrt(CONTROLS.TEMPERATURE/avg_temp);

				if (RANK==0)
				{
					cout << "Average temperature     = " << avg_temp << endl;
					cout << "Velocity scaling factor = " << vscale << endl;	
				}
			
				for(int a1=0;a1<SYSTEM.ATOMS;a1++)
				{
					SYSTEM.VELOCITY[a1].X *= vscale;
					SYSTEM.VELOCITY[a1].Y *= vscale;
					SYSTEM.VELOCITY[a1].Z *= vscale;
				}
			
				avg_temp = 0.0;
			}
			else
			  avg_temp += SYSTEM.TEMPERATURE;
		}
		
		// MPI -- Broadcast the position and optionally velocity from the root to all other processes.
		sync_position(SYSTEM.COORDS, NEIGHBOR_LIST, SYSTEM.VELOCITY, SYSTEM.ATOMS, true);


		
	  	////////////////////////////////////////////////////////////
		// If requested, write the dftbgen output file
		////////////////////////////////////////////////////////////
		
		if ( (CONTROLS.FREQ_DFTB_GEN>0) && ((CONTROLS.STEP+1) % CONTROLS.FREQ_DFTB_GEN == 0) && RANK == 0) 
		{
			GENFILE << setw(5) << right << SYSTEM.ATOMS << " S #Step " << CONTROLS.STEP+1 << " Time " << (CONTROLS.STEP+1) * CONTROLS.DELTA_T_FS << " (fs) Temp " << SYSTEM.TEMPERATURE << " (k)" << endl;

			for(int i=0; i<NATMTYP; i++)				// Replaces GENFILE << "O H" << endl;
				GENFILE << TMP_ATOMTYPE[i] << " ";
			GENFILE << endl;
			
			for (int a1=0; a1<SYSTEM.ATOMS; a1++) 
			{
				
				GENFILE << right << setw(4) << a1+1 << " " << setw(2) << SYSTEM.ATOMTYPE_IDX[a1]+1 << " " 
				       << fixed << setprecision(5) << setw(8) << SYSTEM.COORDS[a1].X << " "
					   << fixed << setprecision(5) << setw(8) << SYSTEM.COORDS[a1].Y << " " 	
					   << fixed << setprecision(5) << setw(8) << SYSTEM.COORDS[a1].Z << endl;    	
			}
		   
			GENFILE << fixed << setprecision(5) << setw(8) << 0.0 << " "
					<< fixed << setprecision(5) << setw(8) << 0.0 << " "
					<< fixed << setprecision(5) << setw(8) << 0.0 << endl;
			
			GENFILE << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.X << " "
					<< fixed << setprecision(5) << setw(8) << 0.0 << " "
					<< fixed << setprecision(5) << setw(8) << 0.0 << endl;
			
			GENFILE << fixed << setprecision(5) << setw(8) << 0.0 << " "
					<< fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.Y << " "
					<< fixed << setprecision(5) << setw(8) << 0.0 << endl;
			
			GENFILE << fixed << setprecision(5) << setw(8) << 0.0 << " "
					<< fixed << setprecision(5) << setw(8) << 0.0 << " "
					<< fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.Z << endl;
		}
		
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
		// If requested, print out the VASP POSCAR file for self-consistent fitting
		////////////////////////////////////////////////////////////
	
		if(CONTROLS.SELF_CONSIST && ((CONTROLS.STEP+1) % CONTROLS.SELF_CONSIST_FREQ == 0) && RANK == 0)
		{
			// Format is:
			// Title
			// Scaling factor (i.e. coords = boxlength * scaling_factor * scaled_cooords)
			// Cell vectors (here we assume orthorhombic)
			// Number of atoms of each type
			// "Direct"
			// Scaled coords, in the order of atoms listed above
			// 
			// In order to eventually set atoms back to original order, store a map
			// ... map format is: <MD idx> <VASP idx> < atom type>
	
			int INSTANCE = int((CONTROLS.STEP+1) / CONTROLS.SELF_CONSIST_FREQ);
			stringstream val_to_str;
			
			val_to_str << INSTANCE;
			val_to_str >> TEMP_STR;
			val_to_str.str("");
			val_to_str.clear();
			
			string POSNAME = "POSCAR_";
			POSNAME.append(TEMP_STR);
			POSNAME.append(".mm_md");
			
			ofstream POSCAR_MAPPER;
			static int VASP_IDX;
			
			if(INSTANCE == 1)	// Then this is the first time we're printing.. Need to write the map file
			{
				VASP_IDX = 0;
				POSCAR_MAPPER.open("POSCAR_MAPPER.mm_md");
			}	

			ofstream POSCAR;
			POSCAR.open(POSNAME.data());
	
			POSCAR << "MM MD SNAPSHOT - NATOMS: " << SYSTEM.ATOMS << endl;
	
			POSCAR << "1.0" << endl;
	
			POSCAR << SYSTEM.BOXDIM.X << " 0.0 0.0" << endl;
			POSCAR << "0.0 " << SYSTEM.BOXDIM.Y << " 0.0" <<  endl;
			POSCAR << "0.0 0.0 " << SYSTEM.BOXDIM.Z << endl;
	
			for(int i=0; i<TMP_NATOMTYPE.size(); i++)
				POSCAR << TMP_NATOMTYPE[i] << " ";
			POSCAR << endl;
	
			POSCAR << "Direct" << endl;
	
			for(int i=0; i<TMP_NATOMTYPE.size(); i++)
			{
				for(int j=0; j<SYSTEM.ATOMS; j++)
				{
					if(TMP_ATOMTYPE[i] == SYSTEM.ATOMTYPE[j])
					{
						POSCAR << SYSTEM.COORDS[j].X/SYSTEM.BOXDIM.X << " " 
							   << SYSTEM.COORDS[j].Y/SYSTEM.BOXDIM.Y << " " 
							   << SYSTEM.COORDS[j].Z/SYSTEM.BOXDIM.Z << endl;
					
						if(INSTANCE == 1)
							POSCAR_MAPPER << j << " " << VASP_IDX++ << " " << SYSTEM.ATOMTYPE[j] << endl;
					}
				}
			}

			POSCAR.close();
			if(INSTANCE == 1)
				POSCAR_MAPPER.close();
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
	
		if(OUT_FORCEFILE.is_open())	// Close the force files 
		{
			OUT_FORCEFILE.close();
			OUT_FORCELABL.close();
		}

		TEMP_MASS = 0.0;
	
		for ( int a = 0; a < SYSTEM.ATOMS; a++ ) 
			TEMP_MASS  += SYSTEM.MASS[a];


			cout << "	Average temperature over run = " << fixed << setprecision(4) << right << temp_sum  / CONTROLS.N_MD_STEPS << " K"   << endl;
			cout << "	Average pressure    over run = " << fixed << setprecision(4) << right << press_sum / CONTROLS.N_MD_STEPS << " GPa" << endl;
	


		// Output final xyz position in the same format as input.xyz for restarting.
		ofstream fxyz;
		fxyz.open("output.xyz");

		fxyz << fixed << setw(5) << setprecision(0) << SYSTEM.ATOMS << endl;
		fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.BOXDIM.X << " ";
		fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.BOXDIM.Y << " ";
		fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.BOXDIM.Z << endl;;
	
		// Note the conversion of our string to a c-style string... printf can't handle c++ strings,
		// since it is a C method.
		// 
	 	// I think there is a way to get at the underlying C char* in a C++ String if you like
	 	// printf better (LEF).	

		for ( int ia = 0; ia < SYSTEM.ATOMS; ia++ ) 
		{
			fxyz << setw(2) << SYSTEM.ATOMTYPE[ia] << " ";
			fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.COORDS[ia].X << " ";
			fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.COORDS[ia].Y << " ";
			fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.COORDS[ia].Z << "    ";
			fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.VELOCITY[ia].X << " ";
			fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.VELOCITY[ia].Y << " ";
			fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.VELOCITY[ia].Z << endl;		
		}
	
		fxyz.close();
		STATISTICS.close();

	}
	
	
	// MPI -- End our setup
	
	#ifdef USE_MPI
		MPI_Finalize();
	#endif
  
}       


static void read_input(MD_JOB_CONTROL & CONTROLS, PES_PLOTS & FF_PLOTS, NEIGHBORS & NEIGHBOR_LIST) 				// UPDATED
{
	if (RANK==0)
		cout << endl << "Reading the simulation control input file..." << endl;
	
	// Define some variable to help with reading
	
	bool   			FOUND_END = false;
	string 			LINE;
	string			TEMP_STR,TEMP_STR_2;
	double			TEMP_DOUB;
	int				TEMP_INT;
	int				TMP_NSCAN;
	stringstream	LINE_PARSER;
	int				ITEM_NO;
	int				ADD_TO_NPLOTS = 0;
	
	// Set some defaults
	
	CONTROLS.IS_LSQ       = false;
	CONTROLS.SELF_CONSIST = false;
	CONTROLS.PLOT_PES     = false;
	CONTROLS.WRAP_COORDS  = true;
	CONTROLS.PRINT_VELOC  = false;
	CONTROLS.NVT_CONV_CUT = 0.10;
	
	FF_PLOTS.INCLUDE_FCUT    = true;
	FF_PLOTS.INCLUDE_CHARGES = true;
	
	// Begin reading
	
	while (FOUND_END == false)
	{
		getline(cin,LINE);
		
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
					
						FF_PLOTS.INCLUDE_FCUT = true;
				}
				if(FF_PLOTS.INCLUDE_CHARGES)
				{
					if (isatty(fileno(stdout)))
						cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "		Including charges" << COUT_STYLE.ENDSTYLE << endl;	
					else
						cout << "		Including charges" << endl;						
					
						FF_PLOTS.INCLUDE_FCUT = true;
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
					
					// Are we doing any 3b scans?
					
					LINE_PARSER >> LINE;

					TMP_NSCAN = 0;

					if (LINE=="scan"  || LINE=="SCAN"  || LINE=="Scan")
					{
						LINE_PARSER >> LINE;
 						
						FF_PLOTS.N_SCAN += int(atof(LINE.c_str()));
						TMP_NSCAN = int(atof(LINE.c_str()));
						if (RANK==0)
							cout << "	Found " << TMP_NSCAN << " scans to run for FF type index " << i << endl;
						
					}
					
					LINE_PARSER >> TEMP_STR_2;

					if(TEMP_STR_2 == "INCLUDE")
					{
						LINE_PARSER >> TEMP_STR_2;
						
						if(TEMP_STR_2 == "2B")
						{
							FF_PLOTS.INCLUDE_2B = true;

							LINE_PARSER >> TEMP_STR_2;
							TEMP_STR_2.append(" ");
							LINE_PARSER >> LINE;
							TEMP_STR_2.append(LINE);	// Should be "PAIRTYPE_PARAMS"
							
							if(TEMP_STR_2 != "PAIRTYPE PARAMS")
							{
								cout << "ERROR: 2B Addition specifications incorrect. Format is:" << endl;
								cout << "PAIRTYPE PARAMS IJ <2b_type_index> IK <2b_type_index> JK <2b_type_index>" << endl;
								exit_run(0);
							}
							
							FF_PLOTS.SEARCH_STRING_2B.push_back(TEMP_STR_2);

							XYZ_INT TMP_FF_TYPES;
							
							LINE_PARSER >> TEMP_STR_2  >> TEMP_STR_2;
							TMP_FF_TYPES.X = int(atof(TEMP_STR_2.data()));
							LINE_PARSER >> TEMP_STR_2  >> TEMP_STR_2;
							TMP_FF_TYPES.Y = int(atof(TEMP_STR_2.data()));
							LINE_PARSER >> TEMP_STR_2  >> TEMP_STR_2;
							TMP_FF_TYPES.Z = int(atof(TEMP_STR_2.data()));
							
							FF_PLOTS.IJ_IK_JK_TYPE.push_back(TMP_FF_TYPES);	
							
							if (RANK==0)
								cout << "		...Adding contributions from ij, ik, and jk pair type indicies: " << TMP_FF_TYPES.X << " " << TMP_FF_TYPES.Y << " " << TMP_FF_TYPES.Z << endl;
						}
						else
						{
							cout << "ERROR: Expected \"2B\" " << endl;
							exit_run(0);
						}
						
						LINE_PARSER.str("");
						LINE_PARSER.clear();
					}	

					LINE_PARSER.str("");
					LINE_PARSER.clear();
					
					for(int j=0; j<TMP_NSCAN; j++)
					{
						cin >> LINE >> LINE; 

						if(LINE == "ij" || LINE == "IJ")
								TEMP_INT = 1;
						else if(LINE == "ik" || LINE == "IK")
							TEMP_INT = 2;
						else if(LINE == "jk" || LINE == "JK")
							TEMP_INT = 3;	
						else
						{
							cout << "ERROR: Unrecognized pair type for fixing: " << LINE << endl;
							cout << "       Allowed types are IJ, IK, an JK" << endl;
							cout << "       (PES FF type index " << i << ", scan index " << j << endl;
							exit_run(0);
						}						
						
						cin >> TEMP_DOUB;
						
						FF_PLOTS.FIX_PAIR_1.push_back(TEMP_INT);
						FF_PLOTS.FIX_VAL_1 .push_back(TEMP_DOUB);
						
						cin >> LINE; 

						if(LINE == "ij" || LINE == "IJ")
								TEMP_INT = 1;
						else if(LINE == "ik" || LINE == "IK")
							TEMP_INT = 2;
						else if(LINE == "jk" || LINE == "JK")
							TEMP_INT = 3;						
						else
						{
							cout << "ERROR: Unrecognized pair type for fixing: " << LINE << endl;
							cout << "       Allowed types are IJ, IK, an JK" << endl;
							cout << "       (PES FF type index " << i << ", scan index " << j << endl;
							exit_run(0);
						}							
						cin >> TEMP_DOUB;
						
						FF_PLOTS.FIX_PAIR_2.push_back(TEMP_INT);
						FF_PLOTS.FIX_VAL_2 .push_back(TEMP_DOUB);						
						
						cin >> LINE >> LINE;
						
						if(LINE == "ij" || LINE == "IJ")
								TEMP_INT = 1;
						else if(LINE == "ik" || LINE == "IK")
							TEMP_INT = 2;
						else if(LINE == "jk" || LINE == "JK")
							TEMP_INT = 3;	
						else
						{
							cout << "ERROR: Unrecognized pair type for scanning: " << LINE << endl;
							cout << "       Allowed types are IJ, IK, an JK" << endl;
							cout << "       (PES FF type index " << i << ", scan index " << j << endl;
							exit_run(0);
						}	
						
						ITEM_NO = FF_PLOTS.FIX_PAIR_1.size() - 1;
											
						if(FF_PLOTS.FIX_PAIR_1[ITEM_NO] == FF_PLOTS.FIX_PAIR_2[ITEM_NO])	
						{
								cout << "ERROR: Fixed pair types must be unique." << endl;
								cout << "       (PES FF type index " << i << ", scan index " << ITEM_NO << endl;
								cout << "       (PES FF type index " << i << ", scan index " << ITEM_NO << endl;
								cout << "      " << FF_PLOTS.FIX_PAIR_1[ITEM_NO] << endl;
								cout << "      " << FF_PLOTS.FIX_PAIR_2[ITEM_NO] << endl;
								cout << "      " << TEMP_INT << endl;								
								exit_run(0);
						}					
						if((FF_PLOTS.FIX_PAIR_1[ITEM_NO] == TEMP_INT)	|| (FF_PLOTS.FIX_PAIR_2[ITEM_NO] == TEMP_INT))
						{
							cout << "ERROR: Fixed and scanned pair types cannot be the same." << endl;
							cout << "       (PES FF type index " << i << ", scan index " << ITEM_NO << endl;
							cout << "      " << FF_PLOTS.FIX_PAIR_1[ITEM_NO] << endl;
							cout << "      " << FF_PLOTS.FIX_PAIR_2[ITEM_NO] << endl;
							cout << "      " << TEMP_INT << endl;
							exit_run(0);
						}
						
						FF_PLOTS.SCAN_PAIR.push_back(TEMP_INT);
						FF_PLOTS.PARENT_TYPE.push_back(FF_PLOTS.TYPE_INDEX[i]);
						
						cin.ignore();
						
						if (RANK==0)
							cout << "	Scanning triplet pair index " << TEMP_INT 
							 << " and holding "                <<  FF_PLOTS.FIX_PAIR_1[j] 
						     << " and "                        <<  FF_PLOTS.FIX_PAIR_2[j] 
							 << " fixed at "                   << FF_PLOTS.FIX_VAL_1[j] 
							 << " and "                        << FF_PLOTS.FIX_VAL_2[j] 
							 << " respectively."               << endl;
					}
					


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
			cin >> LINE; 
			
			if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
			{
				CONTROLS.SELF_CONSIST = true;
				cin >> CONTROLS.SELF_CONSIST_FREQ;
				if (RANK==0)
					cout << "	# SLFCNST #: true... will print POSCAR file every " << CONTROLS.SELF_CONSIST_FREQ << " md steps "<< endl;	
			}
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
		
		else if(LINE.find("# NLAYERS #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.N_LAYERS = int(atof(LINE.data()));
			if (RANK==0)
				cout << "	# NLAYERS #: " << CONTROLS.N_LAYERS << endl;	
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
			cin >> LINE; cin.ignore();
			CONTROLS.COORD_FILE = LINE;
			if (RANK==0)
				cout << "	# CRDFILE #: " << CONTROLS.COORD_FILE << endl;	
		}			
		
		// " Simulation options"
		
		else if(LINE.find("# VELINIT #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if (LINE=="READ")
				CONTROLS.INIT_VEL = false;
			else if (LINE=="GEN")
				CONTROLS.INIT_VEL = true;
			else
			{
				cout << "ERROR: # VELINIT # must be specified as READ or GEN." << endl;
				exit_run(1);	
			}			

			if (RANK==0)
			{
				if(CONTROLS.INIT_VEL)
					cout << "	# VELINIT #: GEN ... generating velocites via box Muller" << endl;	
				else	
					cout << "	# VELINIT #: READ ... reading velocities from coordinate file" << endl;	
			}
		}	
				
		else if(LINE.find("# THRMOST #") != string::npos)
		{
			cin >> LINE; 
			
			if (LINE=="HOOVER")
			{
				CONTROLS.USE_HOOVER_THRMOSTAT = true;
				cin >> LINE; cin.ignore();
				CONTROLS.FREQ_UPDATE_THERMOSTAT = double(atof(LINE.data()));
				if (RANK==0)
					cout << "	# THRMOST #: HOOVER... a Hoover time of " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " will be used." << endl;	
			}
			else if (LINE=="VELSCALE")
			{
				CONTROLS.USE_HOOVER_THRMOSTAT = false;
				cin >> LINE; cin.ignore();	
				CONTROLS.FREQ_UPDATE_THERMOSTAT = double(atof(LINE.data()));
				if (RANK==0)
					cout << "	# THRMOST #: VELSCALE... Velocities will be scaled every " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " MD steps." << endl;	
			}
			else
			{
				cout << "ERROR: # THRMOST # must be specified as HOOVER or VELSCALE, and a Hoover time or velocity " << endl;
				cout << "       scaling frequency must be specified in line. " << endl;
				cout << "       Example: HOOVER 10 "<< endl; 
				exit_run(1);	
			}				

		}			
		
		else if(LINE.find("# PRESSUR #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if (LINE=="ANALYTICAL")
			{
				CONTROLS.USE_NUMERICAL_PRESS = false;
				if (RANK==0)
					cout << "	# PRESSUR #: ANALYTICAL" << endl;	
			}
			else if (LINE=="NUMERICAL")
			{
				CONTROLS.USE_NUMERICAL_PRESS = true;
				if (RANK==0)
					cout << "	# PRESSUR #: NUMERICAL" << endl;	
			}
			else
			{
				cout << "ERROR: # NUMPRES # must be specified as ANALYTICAL or NUMERICAL." << endl;
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
		
		else if(LINE.find("# FRQDFTB #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.FREQ_DFTB_GEN = int(atof(LINE.data()));
			if (RANK==0)
				cout << "	# FRQDFTB #: " << CONTROLS.FREQ_DFTB_GEN << endl;	
			
			if (CONTROLS.DELTA_T_FS > 0 && CONTROLS.N_MD_STEPS > 0 && RANK==0)
			{
				cout << "		... printing every " << CONTROLS.FREQ_DFTB_GEN*CONTROLS.DELTA_T_FS << " fs, " << endl;
				cout << "		... printing " << CONTROLS.N_MD_STEPS/CONTROLS.FREQ_DFTB_GEN << " frames. " << endl; 
			}
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
	}
}
static double kinetic_energy(FRAME & SYSTEM)					// UPDATED -- Overloaded.. compute differentely if for main or new velocities
{
	double Ktot = 0.0;
	
	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
			Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY[a1].X * SYSTEM.VELOCITY[a1].X;
			Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY[a1].Y * SYSTEM.VELOCITY[a1].Y;
			Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY[a1].Z * SYSTEM.VELOCITY[a1].Z;
	}

  return(Ktot);
}
static double kinetic_energy(FRAME & SYSTEM, string TYPE)		// UPDATED -- Overloaded.. compute differentely if for main or new velocities
{
	double Ktot = 0.0;
	
	if(TYPE == "NEW")
	{
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
				Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY_NEW[a1].X * SYSTEM.VELOCITY_NEW[a1].X;
				Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY_NEW[a1].Y * SYSTEM.VELOCITY_NEW[a1].Y;
				Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY_NEW[a1].Z * SYSTEM.VELOCITY_NEW[a1].Z;
		}		
	}
	else
	{
		cout << "ERROR: Requested ke type not understood. Check code." << endl;
		exit_run(0);
	}


  return(Ktot);
}

// MPI -- Related functions -- See headers at top of file

static void sum_forces(vector<XYZ>& accel_vec, int atoms, double &pot_energy, double &pressure)
// Add up forces, potential energy, and pressure from all processes.  
{
	#ifdef USE_MPI

		double *accel_sum = new double [3 * atoms];			// Makes a 1d array of 3*atoms doubles
		
		// Makes 1 1d array of 3(X,Y,Z)*atoms doubles... the "(double*)" recasts the XYZ as a double array, 
		// which flattens it out so [0] = vec[0].X, [1] = vec[0].Y, ... [5] = vec[1].Z and so on.
		
		double *accel    = (double*) accel_vec.data();

		if (sizeof(XYZ) != 3*sizeof(double)) 
		{
			// But then what should be done??
			
			printf("Error: compiler padding in XYZ structure detected\n");
			exit_run(1);
		}

		if ( RANK == 0 ) 
		{
			for (int i=0; i<3*atoms; i++) 
				accel_sum[i] = 0.0;
		}

		#ifdef LOG_FORCE

			char buf[20];
			sprintf(buf, "%d.%d", RANK,NPROCS);
			string force_out = string("force.") + string(buf) + string(".out");
			ofstream outf;
			outf.open(force_out.c_str());
			outf.precision(15);

			for (int i=0; i<3*atoms; i++)
				outf << i << " " << accel[i] << " " << endl;

			outf.close();
			
		#endif

		MPI_Reduce(accel, accel_sum, 3*atoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		#ifdef LOG_FORCE

			ofstream outall;
			char buf2[20];
			sprintf(buf2, "force.all.%d.out", NPROCS);
			printf("OPENING: %s\n", buf2);

			outall.open(buf2);
			outall.precision(15);
			for(int i=0; i<3*atoms; i++)
				outall << i << " " << accel_sum[i] << " " << endl;

			outall.close();
			
		#endif
			
		double petot = 0.0;
		double pe    = pot_energy;

		MPI_Reduce(&pot_energy, &petot, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

		pot_energy = petot;

		double pressure_tot = 0.0;
		double press        = pressure;

		MPI_Reduce(&press, &pressure_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		pressure = pressure_tot;
  
		if ( RANK == 0 ) 
			for ( int i = 0; i < 3 * atoms; i++ ) 
				accel[i] = accel_sum[i];
		
		delete [] accel_sum;

	#endif

}


static void sync_position(vector<XYZ>& coord_vec, NEIGHBORS & neigh_list, vector<XYZ>& velocity_vec, int atoms, bool sync_vel) 
// Broadcast the position, neighborlists,  and optionally the velocity to all nodes.
// Velocity-broadcast is only necessary for velocity-dependent forces (not currently implemented).
{
	#ifdef USE_MPI
	
		// Convert out vectors to arrays so they are MPI friendly
	
		double *coord = (double *) coord_vec.data();

		if ( sizeof(XYZ) != 3*sizeof(double) ) 
		{
			printf("Error: compiler padding in XYZ structure detected\n");
			exit_run(1);
		}

		MPI_Bcast(coord, 3*atoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		#ifdef LOG_POS
		
			char buf[20];
			sprintf(buf, "%d.%d", RANK,NPROCS);
			string pos_out = string("pos.") + string(buf) + string(".out");
			ofstream outx;
			outx.open(pos_out.c_str());
			outx.precision(15);

			for (int i=0; i<3*atoms; i++)
				outx << i << " " << coord[i] << " " << endl;

			outx.close();
		
		#endif

			if ( sync_vel )
			{
				double *velocity = (double *) velocity_vec.data();
				MPI_Bcast(velocity, 3 * atoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			}
			
	#endif
}

void exit_run(int value)
// Call this instead of exit(1) to properly terminate all MPI processes.
{
	#ifdef USE_MPI
		MPI_Abort(MPI_COMM_WORLD,value);
	#else
		exit(value);
	#endif

}
