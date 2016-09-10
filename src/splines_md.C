// USAGE IS: ./a.out < inputfile.in


#include<iostream>
#include<iomanip>
#include<fstream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<sstream>

#include "functions.h"


using namespace std;



static void   read_input        (MD_JOB_CONTROL & CONTROLS, PES_PLOTS & FF_PLOTS);	// UPDATED
static double kinetic_energy    (FRAME & SYSTEM);				// UPDATED
static double kinetic_energy    (FRAME & SYSTEM, string TYPE);	// UPDATED
double        numerical_pressure(const FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP) ;


string FULL_FILE_3B;		// Global variables declared as externs in functions.h, and declared in functions.C
string SCAN_FILE_3B;
string SCAN_FILE_2B;


int main(int argc, char* argv[])
{
	
	////////////////////////////////////////////////////////////
    // Define/initialize important variables
	////////////////////////////////////////////////////////////

	MD_JOB_CONTROL CONTROLS;			// Declare the data object that will hold the main simulation control variables
	PES_PLOTS FF_PLOTS;					// Declare the data object that will help set up PES plots
	read_input(CONTROLS, FF_PLOTS);		// Populate object with user defined values
	
	// Data objects to hold coefficients for different force field types, and for FF printing (if requested)
	 
	vector<PAIR_FF> FF_2BODY;	// Holds all 2-body parameters
	vector<TRIP_FF> FF_3BODY; 	// Holds all 3-body parameters

	// Define the mapping variables that let us figure out which FF params to use for a given pair/triplet of pairs
	
	map<string,int> PAIR_MAP;			// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
	map<int,string> PAIR_MAP_REVERSE; 	// Input/output the resverse of PAIR_MAP
	map<string,int> TRIAD_MAP;			// Input is three of any ATOM_PAIRS.PRPR_NM pairs, output is a triplet type index
	map<int,string> TRIAD_MAP_REVERSE;	// Input/output is the reverse of TRIAD_MAP

	FRAME SYSTEM;						// Declare the data object that will hold the system coordinates, velocities, accelrations, etc.
	
	double Ktot;
	double Vol;
    double avg_temp  = 0.0;
    double zeta_dot0 = 0.0;
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
	vector<int> 	TMP_ATOMTYPEIDX;
	vector<double> 	TMP_CHARGES;
	vector<double> 	TMP_MASS;
	stringstream	STREAM_PARSER;
	
	ofstream GENFILE;			// Holds dftbgen info output.. whatever that is
	ifstream CMPR_FORCEFILE;	// Holds the forces that were read in for comparison purposes
	ofstream OUT_FORCEFILE;		// Holds the forces that are computed and are to be printed out
	ofstream OUT_FORCELABL;		// Holds the forces that are computed and are to be printed out.. has atom labels
	ofstream STATISTICS;		// Holds a copy of the "Step      Time          Ktot/N          Vtot/N ..." part of the output file.. uninterrupted by warnings
	
	string EXTENSION;			// Holds the input xyz/xyzf file extension, so we know how many fields to read on each atom line
	
	ifstream PARAMFILE;
	ifstream COORDFILE;
	
	cout.precision(15);			// Set output precision
	
	////////////////////////////////////////////////////////////
    // Hop to PES printing, if requested 
	////////////////////////////////////////////////////////////
	
	if (FF_PLOTS.N_PLOTS > 0)
		goto FF_SETUP_1; 
	
	////////////////////////////////////////////////////////////
	// Setup an a dftb gen output file, if user has requested it
	////////////////////////////////////////////////////////////

    if ( CONTROLS.FREQ_DFTB_GEN > 0 ) 
		GENFILE.open("traj.gen");
	
	////////////////////////////////////////////////////////////
	// Setup an a force output file, if user has requested it	
	////////////////////////////////////////////////////////////

    if ( CONTROLS.PRINT_FORCE ) 
	{
		OUT_FORCEFILE.open("forceout.txt");	
		OUT_FORCELABL.open("forceout-labeled.txt");		
	}

	////////////////////////////////////////////////////////////
    // Read input file input.xyz, where box dims are on the info line:
	////////////////////////////////////////////////////////////

	COORDFILE.open(CONTROLS.COORD_FILE.data());
	if(!COORDFILE.is_open())
	{
		cout << "ERROR: Cannot open coordinate file: " << CONTROLS.COORD_FILE << endl;
		exit(0);
	}
	
	////////////////////////////////////////////////////////////
	// Set up the system data object... Do so by assuming a 
	// constant number of atoms in the system (almost always
	// true for md)...
	////////////////////////////////////////////////////////////
	
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

	cout << "   ...setup complete" << endl << endl;

	////////////////////////////////////////////////////////////
	// Read in the initial system coordinates, and if requested,
	// initial forces from separate file (i.e. not from .xyzf)
	// ... We also need to figure out the file extension so we
	// know how many fields to read on each atom lne
	////////////////////////////////////////////////////////////


	cout << "Reading initial coordinates and forces..." << endl;

    if ( CONTROLS.COMPARE_FORCE ) 
    {
      cout << "Opening " << CONTROLS.COMPARE_FILE.data() << " to read forces for comparison\n";
      CMPR_FORCEFILE.open(CONTROLS.COMPARE_FILE.data());

	  if(!CMPR_FORCEFILE.is_open())
	  {
		  cout << "ERROR: Cannot open force input file: " << CONTROLS.COMPARE_FILE << endl;
		  exit(0);
	  }

    }
	
	EXTENSION = CONTROLS.COORD_FILE.substr(CONTROLS.COORD_FILE.find_last_of(".")+1);
	cout << "	...Read the following coordinate file extension: " << EXTENSION << endl;
    cout << "	...Read box dimensions: " << SYSTEM.BOXDIM.X << " " << SYSTEM.BOXDIM.Y << " " << SYSTEM.BOXDIM.Z << endl;
	
    for(int a=0; a<SYSTEM.ATOMS ;a++)
	{
        COORDFILE >> SYSTEM.ATOMTYPE[a];
	
		// Read the coordinates
		COORDFILE >> SYSTEM.COORDS[a].X >> SYSTEM.COORDS[a].Y >> SYSTEM.COORDS[a].Z;
		
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
				if(a==0)
				{
					cout << "	...Ignoring last three fields of atom info lines in xyzf file... " << endl;
					cout << "	...will read from specified force file instead." << endl;
				}

				COORDFILE >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;			
			}
					
			// Read forces from separate force file
			CMPR_FORCEFILE >> SYSTEM.FORCES[a].X >> SYSTEM.FORCES[a].Y >> SYSTEM.FORCES[a].Z;
		}
        else if(!CONTROLS.INIT_VEL) // Reading positions from *.xyz
		{
			if(EXTENSION == ".xyz")
			{
				cout << "ERROR: Input file requests velocities to be read in. " << endl;
				cout << "Expected .xyzf file, was given .xyzf file." << endl;
				exit(0);
			}
			
			// Read in velocities instead of forces... I guess this is the format of .xyz files for this code? 
			COORDFILE >> SYSTEM.VELOCITY[a].X >> SYSTEM.VELOCITY[a].Y >> SYSTEM.VELOCITY[a].Z;
		}
		else
		{
			if(EXTENSION == "xyzf")
			{
				if(a==0)
					cout << "	...Ignoring last three fields of atom info lines in xyzf file... " << endl;
				COORDFILE >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;
			}
		}
	}
	
    COORDFILE.close();
    COORDFILE.clear();
	
	if ( CONTROLS.COMPARE_FORCE ) 
		CMPR_FORCEFILE.close();
	
	cout << "   ...read complete" << endl << endl;	
	
	////////////////////////////////////////////////////////////
	// Figure out atom charges and masses, based on parameter 
	// file
	////////////////////////////////////////////////////////////
	 
	FF_SETUP_1: 
	 
	cout << "Reading atom info from parameter file..." << endl; 
	
	// Read in the possible atom types and thier features	
	
	PARAMFILE.open(CONTROLS.PARAM_FILE.data());
	if(!PARAMFILE.is_open())
	{
		cout << "ERROR: Cannot open coordinate file: " << CONTROLS.PARAM_FILE << endl;
		exit(0);
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
				exit(0);
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
			TMP_ATOMTYPEIDX.resize(NATMTYP);
			TMP_CHARGES    .resize(NATMTYP);
			TMP_MASS       .resize(NATMTYP);	
			TMP_SIGN	   .resize(NATMTYP);
			
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
				
				cout << "		" 	<< setw(5) << left << i << " " 
									<< setw(2) << left << TMP_ATOMTYPE[i] << ", q (e): ";
				
				if(CONTROLS.FIT_COUL)
				{
					if(TEMP_STR == "+")
					{
						cout << "POSITIVE";
						TMP_SIGN[i] = 1;
					}
					else
					{
						cout << "NEGATIVE";
						TMP_SIGN[i] = -1;
					}
				}
				else
					cout << setw(6) << fixed << setprecision(3) << right << TMP_CHARGES[i];
				cout  << ", mass (amu): " << setw(8) << fixed << setprecision(4) << right << TMP_MASS[i] << endl;
			}
		}	
	}
	
	// Assign atom features to atoms in SYSTEM data object, and the PAIR_FF object
	
    for(int a=0; a<SYSTEM.ATOMS ;a++)
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

	if (FF_PLOTS.N_PLOTS > 0)

		goto FF_SETUP_2; 


	cout << "   ...read complete" << endl << endl;

	////////////////////////////////////////////////////////////
	// Initialize velocities, if requested. Use the box Muller
	// method, setup up thermostat if requested
	////////////////////////////////////////////////////////////

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
				exit(1);
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
				exit(1);
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
				exit(1);
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
	
	cout << "   ...setup complete" << endl << endl;

	////////////////////////////////////////////////////////////
	// Setup and test system velocity center of mass... For
	// consevation of momentum, this number should zero (or
	// VERY close)
	////////////////////////////////////////////////////////////  

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
    }
	
	// Check our velocity center of mass.. hopefully this is (or is very close to) zero
	  
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
	
	cout	<< "	Final   velocity center of mass: (" 
			<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.X/TEMP_MASS << ", " 
			<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Y/TEMP_MASS << ", " 
			<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Z/TEMP_MASS << ") "
			<< endl;
  
	
  	// Print some info on density, etc.

    dens_mol = (SYSTEM.ATOMS * 1.0e24) / (6.0221e23 * Vol);
    dens_mass = (TEMP_MASS/Vol)*(1e24/6.0221e23);

    cout << "	Total Mass     = " << fixed << setprecision(4) << right << TEMP_MASS << " au"         << endl;
    cout << "	Volume         = " << fixed << setprecision(4) << right << Vol       << " Ang^3"     << endl;
    cout << "	Number density = " << fixed << setprecision(4) << right << dens_mol  << " mol atm/cm^3" << endl;
    cout << "	Mass density   = " << fixed << setprecision(4) << right << dens_mass << " g/cm^3"       << endl;
	
	cout << "   ...checks complete" << endl << endl;

	////////////////////////////////////////////////////////////
	// Begin setup of force field data objects...
	////////////////////////////////////////////////////////////  

	FF_SETUP_2:

	PARAMFILE.open(CONTROLS.PARAM_FILE.data());
	if(!PARAMFILE.is_open())
	{
		cout << "ERROR: Cannot open coordinate file: " << CONTROLS.PARAM_FILE << endl;
		exit(0);
	}
	
	cout << endl << "Reading force field parameters..." << endl;
	
	FOUND_END = false;

	string  TEMP_SEARCH_2B = "some default text";
	string  TEMP_SEARCH_3B = "some default text";
	string	TEMP_TYPE;
	int     NO_PAIRS, NO_TRIPS;
	int		TMP_TERMS1, TMP_TERMS2;
	
	while (FOUND_END == false)
	{
		getline(PARAMFILE,LINE);

		// Break out of loop

		if(LINE.find("ENDFILE") != string::npos)
		{
			FOUND_END = true;
			
			cout << "   ...read complete." << endl << endl;
			
			cout << "Notes on simulation: " << endl;
			
			if(CONTROLS.USE_COULOMB)
			{
				cout << "	Ewald summations will be used ";
				if(CONTROLS.FIT_COUL)
					cout << "and charges will be taken from fit values" << endl;
				else
					cout << "and charges will be taken from user-specified values" << endl;
			}
			else
			{
				cout << "	Electrostatics will not be computed." << endl;
			}
			
			if(CONTROLS.USE_OVERCOORD)
			{
				cout << "	Overbonding contributions will be considered." << endl;
				
				if(CONTROLS.FIT_POVER)
					cout << "		p-over will be read from fit parameters." << endl;
				else
					cout << "		p-over will be read from specified parameters." << endl;
			}
			
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
			
			cout << "		...Compute electrostatics?      " << boolalpha << CONTROLS.USE_COULOMB << endl;
			cout << "		...Use fit charges?             " << boolalpha << CONTROLS.FIT_COUL << endl;
			cout << "		...Compute ReaxFF overbonding?  " << boolalpha << CONTROLS.USE_OVERCOORD << endl;
			cout << "		...Use fit overbonding param?   " << boolalpha << CONTROLS.FIT_POVER << endl;
			cout << "		...Use 3-body Cheby params?     " << boolalpha << CONTROLS.USE_3B_CHEBY << endl;
			
			cout << "	...Read FF controls..." << endl;
			
			PARAMFILE.ignore();
/*			
			if(CONTROLS.FIT_COUL)
			{
				cout << endl;
				cout << "	********************* WARNING ********************* " << endl;
				cout << "	Use FIT_COUL option only supported for H2O systems!" << endl;
				cout << "	********************* WARNING ********************* " << endl;
				cout << endl;
			}
*/			

		}
		
		// Determine the pair type and corresponding orders, etc
		
		else if(LINE.find("PAIRTYP: ") != string::npos)
		{
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
			}			
			
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();	
			
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
					FF_2BODY[i].PENALTY_DIST  = 0.01;
			
				for(int i=0; i<NO_PAIRS; i++)
					FF_2BODY[i].PENALTY_SCALE = 1.0e8;
				
				for(int i=0; i<NO_PAIRS; i++)
					FF_2BODY[i].CUBIC_SCALE   = 1.0;
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
				PARAMFILE >> FF_2BODY[i].S_DELTA;
				
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

			TEMP_SEARCH_3B = "TRIPLET ";
			TEMP_SEARCH_3B.append(FF_2BODY[0].PAIRTYP); // Syntax ok b/c all pairs have same FF type, and 2b and 3b are same pair type
			TEMP_SEARCH_3B.append(" PARAMS");	
			
			cout << "	...Read FF triplet specifications..." << endl;
			
		}
		
		// Read in pair parameters
		
		else if(LINE.find(TEMP_SEARCH_2B) != string::npos) // "PAIR <PAIRTYPE> PARAMS"
		{	

			// Read in the specific atom pair parameters
			
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
					if(i==0)
					{
						cout << "		Will use a simple linear force model for poorly sampled positions" << endl;
						cout << "			i.e. for |spline coefficients| < 1.0 at close separation distance..." << endl;					
					}
					
					STOP_FILL_IDX = -1;
										
					for(int j=0; j<FF_2BODY[i].SNUM; j++)
					{
						if(abs(FF_2BODY[i].PARAMS[j])>1.0)
						{
							STOP_FILL_IDX = j;
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
							
							cout << "			   " << j << " " << FF_2BODY[i].PARAMS[j] << endl << "			   " << FF_2BODY[i].PARAMS[j+1] << endl;
						}
					}
					
					// Compute the integral of the spline equation for use in analytical pressure
					// calculations...
					//  We are computing an integral, so we will only have half as many points
					
					
					FF_2BODY[i].POT_PARAMS.resize(FF_2BODY[i].SNUM/2);
					
					for(int j=0; j<FF_2BODY[i].POT_PARAMS.size(); j++) 
						FF_2BODY[i].POT_PARAMS[j] = 0 ;
					
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
/*					 
					// Now actually set the charges for all the atoms
					
					double TMP_qHH;
					double TMP_qOO;
					 
					for(int j=0; j<FF_2BODY.size(); j++)
						if(FF_2BODY[j].PRPR_NM == "HH")
							TMP_qHH = FF_2BODY[j].PAIR_CHRG;
					 
					TMP_qHH = sqrt(TMP_qHH);
					TMP_qOO = -2.0*TMP_qHH;		
				
					/*
					cout << "Will use the following atom charges: " << endl;
					cout << "O: " << TMP_qOO << endl;
					cout << "H: " << TMP_qHH << endl;	
					*/
/*					
					if(!CONTROLS.PLOT_PES)		
					{
						for(int i=0; i<SYSTEM.ATOMS; i++)
						{
							if(SYSTEM.ATOMTYPE[i] == "O")
								SYSTEM.CHARGES[i] = TMP_qOO;
							else if(SYSTEM.ATOMTYPE[i] == "H")
								SYSTEM.CHARGES[i] = TMP_qHH;
							else
							{
								cout << "ERROR: FIT_COUL functionality only valid for H2O-type systems." << endl;
								exit(0);
							}
						}
					}
*/					
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
						if( FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] && FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j] || FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[j] && FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] )
						{
							TMP_CHARGES[j] = sqrt(FF_2BODY[i].PAIR_CHRG)*TMP_SIGN[j];
							FF_2BODY[i].ATM1CHG = TMP_CHARGES[j];
							FF_2BODY[i].ATM2CHG = TMP_CHARGES[j];
						
							break;
						}
					}
				}
				
			    for(int a=0; a<SYSTEM.ATOMS ;a++)
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
				
				cout << "		Re-setting individual atom charges based on pair charges :" << endl;
				for(int j=0; j<NATMTYP; j++)
				{
					cout << "		"<<	j << "     "
									<< setw(2) << left << TMP_ATOMTYPE[j] << ", q (e): " 
									<< setw(6) << fixed << setprecision(3) << right << TMP_CHARGES[j] << ", mass (amu): " 
									<< setw(8) << fixed << setprecision(4) << right << TMP_MASS[j] << endl;
				}
				 
			}
			
			
			
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
			
			cout << "	Reading  " << TMP_TERMS1 << " pairs for mapping" << endl;
			
			for(int i=0; i<TMP_TERMS1; i++)
			{
				PARAMFILE >> TMP_TERMS2;
				PARAMFILE >> TEMP_TYPE;
				
				cout << "	........Reading pair: " << TEMP_TYPE << " with mapped index: " << TMP_TERMS2 << endl; 
				
				PAIR_MAP.insert(make_pair(TEMP_TYPE,TMP_TERMS2));
				PAIR_MAP_REVERSE.insert(make_pair(TMP_TERMS2,TEMP_TYPE));				
					
			}			
			
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
			
			cout << "	Reading  " << TMP_TERMS1 << " triplets for mapping" << endl;
			
			for(int i=0; i<TMP_TERMS1; i++)
			{
				PARAMFILE >> TMP_TERMS2;
				PARAMFILE >> TEMP_TYPE;
				
				cout << "	........Reading triplet: " << TEMP_TYPE << " with mapped index: " << TMP_TERMS2 << endl; 
				
				TRIAD_MAP.insert(make_pair(TEMP_TYPE,TMP_TERMS2));
				TRIAD_MAP_REVERSE.insert(make_pair(TMP_TERMS2,TEMP_TYPE));				
					
			}
			
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

//	vector<string>().swap(TMP_ATOMTYPE); 
//	vector<double>().swap(TMP_CHARGES); 
	vector<double>().swap(TMP_MASS);	
	
	////////////////////////////////////////////////////////////
	// Print a summary of the force field
	////////////////////////////////////////////////////////////  	
	
	cout << "Force field summary: " << endl;
	
	cout << "	Pair type and number of parameters per pair: " << endl;
	cout << "		" << FF_2BODY[0].PAIRTYP << " " << FF_2BODY[0].SNUM << endl << endl;
	
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
			cout << "		" << FF_2BODY[i].PRPR_NM << " " << FF_2BODY[i].PAIR_CHRG << endl;;
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
	
	////////////////////////////////////////////////////////////
	// If PES printing is requested, do so here and then 
	// exit the program
	////////////////////////////////////////////////////////////  	
	
	if (FF_PLOTS.N_PLOTS > 0)
	{
		if(FF_PLOTS.INCLUDE_2B)
			cout << "WARNING: 2B+3B only used for 2D 3B scans, not 4D 3B scans." << endl;
		
		
		int ij, ik, jk;
		string ATM_TYP_1, ATM_TYP_2, ATM_TYP_3;
		int ADD_2B_TO_3B_IDX = 0;
		
		int IJ_STEP = 0;
		int IK_STEP = 0;
		int JK_STEP = 0;
		
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
								exit(0);
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
						
							Print_Cheby(FF_2BODY, FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].X, PAIR_MAP_REVERSE[FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].X],"for_3b");						
							cout << "		Reading in 2B scans from file (1): "	<< SCAN_FILE_2B << endl;					
						
							SCAN_INFILE_2B.open(SCAN_FILE_2B.data());

							if(!SCAN_INFILE_2B.is_open())
							{
								cout << "ERROR-2: Cannot open file " << SCAN_FILE_2B << " for 2B addition." << endl;
								exit(0);
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
							
							Print_Cheby(FF_2BODY,  FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].Y, PAIR_MAP_REVERSE[FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].Y],"for_3b");
							
							cout << "		Reading in 2B scans from file (2): "	<< SCAN_FILE_2B << endl;					
						
							SCAN_INFILE_2B.open(SCAN_FILE_2B.data());

							if(!SCAN_INFILE_2B.is_open())
							{
								cout << "ERROR-2: Cannot open file " << SCAN_FILE_2B << " for 2B addition." << endl;
								exit(0);
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

							Print_Cheby(FF_2BODY,  FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].Z, PAIR_MAP_REVERSE[FF_PLOTS.IJ_IK_JK_TYPE[scan_2b_idx].Z],"for_3b");
							
							cout << "		Reading in 2B scans from file (3): "	<< SCAN_FILE_2B << endl;					
						
							SCAN_INFILE_2B.open(SCAN_FILE_2B.data());

							if(!SCAN_INFILE_2B.is_open())
							{
								cout << "ERROR-2: Cannot open file " << SCAN_FILE_2B << " for 2B addition." << endl;
								exit(0);
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
					
							for(int a=0; a<PES_VAL_3B.size(); a++)
							{
								for(int b=0; b<PES_VAL_2B_IJ.size(); b++)
								{
									if(IJ_DIST_2B[b] == IJ_DIST_3B[a])
									{					
										PES_VAL_3B[a] += PES_VAL_2B_IJ[b];
										break;
									}
								}
								for(int c=0; c<PES_VAL_2B_IK.size(); c++)
								{
									if(IK_DIST_2B[c] == IK_DIST_3B[a])
									{							
										PES_VAL_3B[a] += PES_VAL_2B_IJ[c];
										break;
									}
								}	
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
								exit(0);
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
				cout << "	and atom types:            " << FF_2BODY[i].ATM1TYP << " " << FF_2BODY[i].ATM1TYP << endl;
					
				Print_Cheby(FF_2BODY, FF_PLOTS.TYPE_INDEX[i], PAIR_MAP_REVERSE[FF_PLOTS.TYPE_INDEX[i]],"");
			}
			else
			{
				cout << "Functionality not programmed yet." << endl;
			}

			cout << "	...Force field PES printing complete." << endl;
		
		idx_3b++;
			
		}
		
		exit(0);
	}
	

  	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	//
	//					START THE SIMULATION
	//
	////////////////////////////////////////////////////////////  
	////////////////////////////////////////////////////////////

	STATISTICS.open("md_statistics.out");

	cout << "BEGIN SIMULATION:" << endl;
	
	double ke; // A temporary variable used when updating thermostatting 
	
	for(int A=0;A<CONTROLS.N_MD_STEPS;A++)	//start Big Loop here.
    {
      
	  	////////////////////////////////////////////////////////////
		// Do first half of coordinate/velocity updating
		////////////////////////////////////////////////////////////
		
		if(A>0)	
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
			
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
			{
				SYSTEM.COORDS[a1].X -= floor(SYSTEM.COORDS[a1].X / SYSTEM.BOXDIM.X) * SYSTEM.BOXDIM.X;
				SYSTEM.COORDS[a1].Y -= floor(SYSTEM.COORDS[a1].Y / SYSTEM.BOXDIM.Y) * SYSTEM.BOXDIM.Y;
				SYSTEM.COORDS[a1].Z -= floor(SYSTEM.COORDS[a1].Z / SYSTEM.BOXDIM.Z) * SYSTEM.BOXDIM.Z;
			}	

			// Update first half of velocity (i.e. using a(t)):
			
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
			{
				SYSTEM.VELOCITY[a1].X += 0.5 * SYSTEM.ACCEL[a1].X * CONTROLS.DELTA_T;
				SYSTEM.VELOCITY[a1].Y += 0.5 * SYSTEM.ACCEL[a1].Y * CONTROLS.DELTA_T;
				SYSTEM.VELOCITY[a1].Z += 0.5 * SYSTEM.ACCEL[a1].Z * CONTROLS.DELTA_T;
				
			}
			
			// Apply thermostatting to velocities

			if ( CONTROLS.USE_HOOVER_THRMOSTAT ) 
			{
				for(int a1=0;a1<SYSTEM.ATOMS;a1++)
				{
					SYSTEM.VELOCITY[a1].X -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].X * CONTROLS.DELTA_T;
					SYSTEM.VELOCITY[a1].Y -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].Y * CONTROLS.DELTA_T;
					SYSTEM.VELOCITY[a1].Z -= 0.5*THERMOSTAT.VISCO * SYSTEM.VELOCITY[a1].Z * CONTROLS.DELTA_T;
				}
			}
		}      
		cout.precision(14);			// Set output precision
    
      	////////////////////////////////////////////////////////////
		// Calculate acceleration
		////////////////////////////////////////////////////////////

		// this function calculates the spline and electrostatic forces. --- this should probably be removed
		ZCalc(SYSTEM, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP);

		////////////////////////////////////////////////////////////
		// Print some info on new forces, compare to input forces if
		// requested
		////////////////////////////////////////////////////////////

	  	if ( CONTROLS.PRINT_FORCE && (A+1)%CONTROLS.FREQ_FORCE == 0 ) 
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

		if ( CONTROLS.COMPARE_FORCE ) 
		{
			// Check against read-in forces for code verification... Note, ferr is initialized to zero.
			
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
			{
					ferr += (SYSTEM.ACCEL[a1].X - SYSTEM.FORCES[a1].X) * (SYSTEM.ACCEL[a1].X - SYSTEM.FORCES[a1].X);
					ferr += (SYSTEM.ACCEL[a1].Y - SYSTEM.FORCES[a1].Y) * (SYSTEM.ACCEL[a1].Y - SYSTEM.FORCES[a1].Y);
					ferr += (SYSTEM.ACCEL[a1].Z - SYSTEM.FORCES[a1].Z) * (SYSTEM.ACCEL[a1].Z - SYSTEM.FORCES[a1].Z);
			}
			
			ferr /= 3.0 * SYSTEM.ATOMS;
			ferr = sqrt(ferr);
			
			printf("RMS force error = %13.6e\n", ferr);
			exit(0);
		}
		else	// Print main simulation header 
		{
		  	////////////////////////////////////////////////////////////
			// Print the header for mains simulation output
			////////////////////////////////////////////////////////////		

			if ( A == 0 ) 
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

		if(A>0)
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
		
	  	////////////////////////////////////////////////////////////
		// If requested, compute pressure numerically, and accumulate
		// statistics
		////////////////////////////////////////////////////////////
	  
		if ( CONTROLS.USE_NUMERICAL_PRESS ) 
			SYSTEM.PRESSURE_XYZ = numerical_pressure(SYSTEM, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP);

		SYSTEM.PRESSURE = (SYSTEM.PRESSURE_XYZ + 2.0 * Ktot / (3.0 * Vol)) * GPa;	// GPa = Unit conversion factor to GPa.
		press_sum += SYSTEM.PRESSURE;
		
		
	  	////////////////////////////////////////////////////////////
		// Periodically print simulation output
		////////////////////////////////////////////////////////////
		
		if ( (A+1) % CONTROLS.FREQ_ENER == 0 ) 
		{
			
			printf("%8d %9.2f %15.7f %15.7f %15.7f %15.1f %15.3f", A+1, (A+1)*CONTROLS.DELTA_T_FS, Ktot/SYSTEM.ATOMS,SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS,(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS,SYSTEM.TEMPERATURE, SYSTEM.PRESSURE);
			STATISTICS << A+1<< "	" << (A+1)*CONTROLS.DELTA_T_FS<< "	" << Ktot/SYSTEM.ATOMS<< "	" <<SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS<< "	" <<(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS<< "	" <<SYSTEM.TEMPERATURE<< "	" << SYSTEM.PRESSURE;

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

		if ( (!CONTROLS.USE_HOOVER_THRMOSTAT) && ((A+1) % int(CONTROLS.FREQ_UPDATE_THERMOSTAT) == 0) )
		{
			avg_temp +=  SYSTEM.TEMPERATURE;
			avg_temp /= int(CONTROLS.FREQ_UPDATE_THERMOSTAT);
			vscale    = sqrt(CONTROLS.TEMPERATURE/avg_temp);

			cout << "Average temperature     = " << avg_temp << endl;
			cout << "Velocity scaling factor = " << vscale << endl;
			
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


	  	////////////////////////////////////////////////////////////
		// If requested, write the dftbgen output file
		////////////////////////////////////////////////////////////
		

		if ( (CONTROLS.FREQ_DFTB_GEN>0) && ((A+1) % CONTROLS.FREQ_DFTB_GEN == 0) ) 
		{
			GENFILE << setw(5) << right << SYSTEM.ATOMS << " S #Step " << A+1 << " Time " << (A+1) * CONTROLS.DELTA_T_FS << " (fs) Temp " << SYSTEM.TEMPERATURE << " (k)" << endl;

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
	  	////////////////////////////////////////////////////////////
	  	// 
		// 				END LOOP OVER MD STEPS
		// 
		////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////
   
	   
    }//End big loop here.
    
	cout << "END SIMULATION" << endl;

	if(OUT_FORCEFILE.is_open())
	{
		OUT_FORCEFILE.close();
		OUT_FORCELABL.close();
	}
	
	// WHY ON EARTH WOULD WE NEED TO RE-COMPUTE MASS??? WHAT COULD POSSIBLY
	// HAPPEN TO CAUSE THE MASS TO CHANGE DURING A CLOSED-SYSTEM MD SIMULATION?!?!?!?!?!?!??!?!?!?!?!??!?!?!
	
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


static void read_input(MD_JOB_CONTROL & CONTROLS, PES_PLOTS & FF_PLOTS) 				// UPDATED
{
	cout << endl << "Reading the simulation control input file..." << endl;
	
	bool   			FOUND_END = false;
	string 			LINE;
	string			TEMP_STR,TEMP_STR_2;
	double			TEMP_DOUB;
	int				TEMP_INT;
	int				TMP_NSCAN;
	stringstream	LINE_PARSER;
	int				ITEM_NO;
	int				ADD_TO_NPLOTS = 0;
	
	
	while (FOUND_END == false)
	{
		getline(cin,LINE);
		
		// Break out of the loop
		
		if     (LINE.find("# ENDFILE #") != string::npos)
		{
			FOUND_END = true;
			cout << "   ...read complete." << endl << endl;
			break;
		}
		
		// Variables for printing out PES for a given parameter file
		
		else if(LINE.find("# PLOTPES #") != string::npos)
		{
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
								exit(0);
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
								exit(0);
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
							
							cout << "		...Adding contributions from ij, ik, and jk pair type indicies: " << TMP_FF_TYPES.X << " " << TMP_FF_TYPES.Y << " " << TMP_FF_TYPES.Z << endl;
						}
						else
						{
							cout << "ERROR: Expected \"2B\" " << endl;
							exit(0);
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
							exit(0);
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
							exit(0);
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
							exit(0);
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
								exit(0);
						}					
						if((FF_PLOTS.FIX_PAIR_1[ITEM_NO] == TEMP_INT)	|| (FF_PLOTS.FIX_PAIR_2[ITEM_NO] == TEMP_INT))
						{
							cout << "ERROR: Fixed and scanned pair types cannot be the same." << endl;
							cout << "       (PES FF type index " << i << ", scan index " << ITEM_NO << endl;
							cout << "      " << FF_PLOTS.FIX_PAIR_1[ITEM_NO] << endl;
							cout << "      " << FF_PLOTS.FIX_PAIR_2[ITEM_NO] << endl;
							cout << "      " << TEMP_INT << endl;
							exit(0);
						}
						
						FF_PLOTS.SCAN_PAIR.push_back(TEMP_INT);
						FF_PLOTS.PARENT_TYPE.push_back(FF_PLOTS.TYPE_INDEX[i]);
						
						cin.ignore();
						
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
				cout << "   ...read complete." << endl << endl;
				FF_PLOTS.N_PLOTS += ADD_TO_NPLOTS;
				break;				
			}
			else
				cout << "	# PLOTPES #: false" << endl;	
		}	
		
		
		// "General control variables"
		
		else if(LINE.find("# RNDSEED #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.SEED = int(atof(LINE.data()));
			cout << "	# RNDSEED #: " << CONTROLS.SEED << endl;	
		}	
		
		else if(LINE.find("# TEMPERA #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.TEMPERATURE = double(atof(LINE.data()));
			cout << "	# TEMPERA #: " << CONTROLS.TEMPERATURE << " K" << endl;	
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
				}
				else
					CONTROLS.COMPARE_FILE = TEMP_STR;
			}
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				CONTROLS.COMPARE_FORCE = false;
			else
			{
				cout << "ERROR: # CMPRFRC # must be specified as true or false." << endl;
				exit(1);	
			}					
				
			if(CONTROLS.COMPARE_FORCE)
			{
				cout << "	# CMPRFRC #: true... will only do 1 md step." << endl;	
				
				// If CONTROLS.COMPARE_FORCE is true, then only do one step,
				// because all we're trying to do is to see if md computed
				// forces match the expected forces

				CONTROLS.N_MD_STEPS = 1;
			}
			else
				cout << "	# CMPRFRC #: false" << endl;	
		}		
			
		else if(LINE.find("# TIMESTP #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.DELTA_T_FS = double(atof(LINE.data()));
			cout << "	# TIMESTP #: " << CONTROLS.DELTA_T_FS << " fs" << endl;
			
			// Now convert timestep to simulation units
			CONTROLS.DELTA_T = CONTROLS.DELTA_T_FS/Tfs;
		}	
		
		else if(LINE.find("# N_MDSTP #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.N_MD_STEPS = int(atof(LINE.data()));
			cout << "	# N_MDSTP #: " << CONTROLS.N_MD_STEPS << endl;	
		}	
		
		else if(LINE.find("# NLAYERS #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.N_LAYERS = int(atof(LINE.data()));
			cout << "	# NLAYERS #: " << CONTROLS.N_LAYERS << endl;	
		}	
		
		else if(LINE.find("# PRMFILE #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.PARAM_FILE = LINE;
			cout << "	# PRMFILE #: " << CONTROLS.PARAM_FILE << endl;	
		}	
		
		else if(LINE.find("# CRDFILE #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.COORD_FILE = LINE;
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
				exit(1);	
			}			

			if(CONTROLS.INIT_VEL)
				cout << "	# VELINIT #: GEN ... generating velocites via box Muller" << endl;	
			else
				cout << "	# VELINIT #: READ ... reading velocities from coordinate file" << endl;	
		}	
				
		else if(LINE.find("# THRMOST #") != string::npos)
		{
			cin >> LINE; 
			
			if (LINE=="HOOVER")
			{
				CONTROLS.USE_HOOVER_THRMOSTAT = true;
				cin >> LINE; cin.ignore();
				CONTROLS.FREQ_UPDATE_THERMOSTAT = double(atof(LINE.data()));
				cout << "	# THRMOST #: HOOVER... a Hoover time of " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " will be used." << endl;	
			}
			else if (LINE=="VELSCALE")
			{
				CONTROLS.USE_HOOVER_THRMOSTAT = false;
				cin >> LINE; cin.ignore();	
				CONTROLS.FREQ_UPDATE_THERMOSTAT = double(atof(LINE.data()));
				cout << "	# THRMOST #: VELSCALE... Velocities will be scaled every " << CONTROLS.FREQ_UPDATE_THERMOSTAT << " MD steps." << endl;	
			}
			else
			{
				cout << "ERROR: # THRMOST # must be specified as HOOVER or VELSCALE, and a Hoover time or velocity " << endl;
				cout << "       scaling frequency must be specified in line. " << endl;
				cout << "       Example: HOOVER 10 "<< endl; 
				exit(1);	
			}				

		}			
		
		else if(LINE.find("# PRESSUR #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if (LINE=="ANALYTICAL")
			{
				CONTROLS.USE_NUMERICAL_PRESS = false;
				cout << "	# PRESSUR #: ANALYTICAL" << endl;	
			}
			else if (LINE=="NUMERICAL")
			{
				CONTROLS.USE_NUMERICAL_PRESS = true;
				cout << "	# PRESSUR #: NUMERICAL" << endl;	
			}
			else
			{
				cout << "ERROR: # NUMPRES # must be specified as ANALYTICAL or NUMERICAL." << endl;
				exit(1);	
			}				
		}	
				
		// "Output control"
		
		else if(LINE.find("# FRQDFTB #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.FREQ_DFTB_GEN = int(atof(LINE.data()));
			cout << "	# FRQDFTB #: " << CONTROLS.FREQ_DFTB_GEN << endl;	
		}

		else if(LINE.find("# FRQENER #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.FREQ_ENER = int(atof(LINE.data()));
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
				
				cout << "	# PRNTFRC #: true ...and will be printed every " << CONTROLS.FREQ_FORCE << " frames."<< endl;	
//				cout << "	# PRNTFRC #: true" << endl;	
			}
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
			{
				CONTROLS.PRINT_FORCE = false;
				cout << "	# PRNTFRC #: false" << endl;	
			}
			else
			{
				cout << "ERROR: # PRNTFRC # must be specified as true or false." << endl;
				exit(1);	
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
	}


  return(Ktot);
}




