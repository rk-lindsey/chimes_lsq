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

#ifndef VERBOSITY
	#define VERBOSITY 1
#endif

static void   read_input        (MD_JOB_CONTROL & CONTROLS);	// UPDATED
static double kinetic_energy    (FRAME & SYSTEM);				// UPDATED
static double kinetic_energy    (FRAME & SYSTEM, string TYPE);	// UPDATED
double        numerical_pressure(const FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP) ;

int main(int argc, char* argv[])
{
	
	////////////////////////////////////////////////////////////
    // Define/initialize important variables
	////////////////////////////////////////////////////////////

	MD_JOB_CONTROL CONTROLS;	// Declare the data object that will hold the main simulation control variables
	read_input(CONTROLS);		// Populate object with user defined values
	 
	// Data objects to hold coefficients for different force field types.
	 
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
	
	ofstream GENFILE;			// Holds dftbgen info output.. whatever that is
	ifstream CMPR_FORCEFILE;	// Holds the forces that were read in for comparison purposes
	ofstream OUT_FORCEFILE;		// Holds the forces that are computed and are to be printed out
	
	cout.precision(15);			// Set output precision
	
	
	
	////////////////////////////////////////////////////////////
	// Setup an a dftb gen output file, if user has requested it
	////////////////////////////////////////////////////////////

    if ( CONTROLS.FREQ_DFTB_GEN > 0 ) 
		GENFILE.open("traj.gen");
	
	////////////////////////////////////////////////////////////
	// Setup an a force output file, if user has requested it	
	////////////////////////////////////////////////////////////

    if ( CONTROLS.PRINT_FORCE ) 
		OUT_FORCEFILE.open("forceout.txt");	

	////////////////////////////////////////////////////////////
    // Read input file input.xyz, where box dims are on the info line:
	////////////////////////////////////////////////////////////

    ifstream COORDFILE;
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

	cout << "   ...setup complete" << endl << endl;

	////////////////////////////////////////////////////////////
	// Read in the initial system coordinates, and if requested,
	// initial forces from separate file (i.e. not from .xyzf)
	////////////////////////////////////////////////////////////


	cout << "Reading initial coordinates and forces..." << endl;

	string TEMP_STR;
	
    if ( CONTROLS.COMPARE_FORCE ) 
    {
      cout << "Opening force.txt to read forces for comparison\n";
      CMPR_FORCEFILE.open(CONTROLS.COMPARE_FILE.data());

	  if(!CMPR_FORCEFILE.is_open())
	  {
		  cout << "ERROR: Cannot open force input file: " << CONTROLS.COMPARE_FILE << endl;
		  exit(0);
	  }

    }

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
			// Ignore forces in .xyzf file
			COORDFILE >> TEMP_STR           >> TEMP_STR           >> TEMP_STR;
			
			// Read forces from separate force file
			CMPR_FORCEFILE >> SYSTEM.FORCES[a].X >> SYSTEM.FORCES[a].Y >> SYSTEM.FORCES[a].Z;
		}
        else // Reading positions from *.xyz
		{
			// Read in velocities instead of forces... I guess this is the format of .xyz files for this code? 
			COORDFILE >> SYSTEM.VELOCITY[a].X >> SYSTEM.VELOCITY[a].Y >> SYSTEM.VELOCITY[a].Z;
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
	 
	cout << "Reading atom info from parameter file..." << endl; 
	
	// Read in the possible atom types and thier features	
	
    ifstream PARAMFILE;
	PARAMFILE.open(CONTROLS.PARAM_FILE.data());
	if(!PARAMFILE.is_open())
	{
		cout << "ERROR: Cannot open coordinate file: " << CONTROLS.PARAM_FILE << endl;
		exit(0);
	}
	
	bool   	FOUND_END = false;
	string 	LINE;
	
	int				NATMTYP = 0;
	vector<string> 	TMP_ATOMTYPE;
	vector<double> 	TMP_CHARGES;
	vector<double> 	TMP_MASS;
	stringstream	STREAM_PARSER;
	
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
			
			TMP_ATOMTYPE.resize(NATMTYP);
			TMP_CHARGES .resize(NATMTYP);
			TMP_MASS    .resize(NATMTYP);	
			
			cout << "	Read " << NATMTYP << " atom types:" << endl;			
		}			
		
		else if(LINE.find("# TYPEIDX #") != string::npos)
		{
			for (int i=0; i<NATMTYP; i++)
			{
				getline(PARAMFILE,LINE);
				STREAM_PARSER.str(LINE);
				
				STREAM_PARSER >> TEMP_STR;
				STREAM_PARSER >> TMP_ATOMTYPE[i];
				STREAM_PARSER >> TMP_CHARGES[i];
				STREAM_PARSER >> TMP_MASS[i];
				
				STREAM_PARSER.str("");
				STREAM_PARSER.clear();
				
				cout << "		" 	<< setw(5) << left << TEMP_STR << " " 
									<< setw(2) << left << TMP_ATOMTYPE[i] << ", q (e): " 
									<< setw(6) << fixed << setprecision(3) << right << TMP_CHARGES[i] << ", mass (amu): " 
									<< setw(8) << fixed << setprecision(4) << right << TMP_MASS[i] << endl;
			}
		}	
	}
	
	// Assign atom features to atoms in SYSTEM data object
	
    for(int a=0; a<SYSTEM.ATOMS ;a++)
	{
		for(int i=0; i<NATMTYP; i++)
		{
			if(SYSTEM.ATOMTYPE[a] == TMP_ATOMTYPE[i])
			{
				SYSTEM.CHARGES[a] = TMP_CHARGES[i];
				SYSTEM.MASS[a]    = TMP_MASS[i];
				break;
			}			
		}		
	}
	
	// Free up some memory.. swap the contents of currend vectors with a vector w/ no assigned mem
	
	vector<string>().swap(TMP_ATOMTYPE);
	vector<double>().swap(TMP_CHARGES);
	vector<double>().swap(TMP_MASS);	
	
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
			
			cout << "		...Compute electrostatics?      " << boolalpha << CONTROLS.FIT_COUL << endl;
			cout << "		...Use fit charges?             " << boolalpha << CONTROLS.USE_COULOMB << endl;
			cout << "		...Compute ReaxFF overbonding?  " << boolalpha << CONTROLS.USE_OVERCOORD << endl;
			cout << "		...Use fit overbonding param?   " << boolalpha << CONTROLS.FIT_POVER << endl;
			cout << "		...Use 3-body Cheby params?     " << boolalpha << CONTROLS.USE_3B_CHEBY << endl;
			
			cout << "	...Read FF controls..." << endl;
			
			PARAMFILE.ignore();
			
			if(CONTROLS.FIT_COUL)
			{
				cout << endl;
				cout << "	********************* WARNING ********************* " << endl;
				cout << "	Use FIT_COUL option only supported for H2O systems!" << endl;
				cout << "	********************* WARNING ********************* " << endl;
				cout << endl;
			}

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
					 
					// Now actually set the charges for all the atoms
					
					double TMP_qHH;
					double TMP_qOO;
					 
					for(int i=0; i<FF_2BODY.size(); i++)
						if(FF_2BODY[i].PRPR_NM == "HH")
							TMP_qHH = FF_2BODY[i].PAIR_CHRG;
					 
					TMP_qHH = sqrt(TMP_qHH);
					TMP_qOO = -2.0*TMP_qHH;		
					
					/*
					cout << "Will use the following atom charges: " << endl;
					cout << "O: " << TMP_qOO << endl;
					cout << "H: " << TMP_qHH << endl;	
					*/		
					 
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
			cout << "cheby type...cheby lambda";
		cout << endl;
				
		cout << "		" << FF_2BODY[i].PRPR_NM << " ";
		
		cout << FF_2BODY[i].S_MINIM << " ";
		cout << FF_2BODY[i].S_MAXIM << " ";
		cout << FF_2BODY[i].S_DELTA << " ";
		
		if(FF_2BODY[i].PAIRTYP == "CHEBYSHEV")
		{
			cout << FF_2BODY[i].CHEBY_TYPE << " ";
			cout << FF_2BODY[i].LAMBDA << " ";
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
	////////////////////////////////////////////////////////////
	//
	//					START THE SIMULATION
	//
	////////////////////////////////////////////////////////////  
	////////////////////////////////////////////////////////////

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

	  	if ( CONTROLS.PRINT_FORCE ) 
	  	{
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
			{
				OUT_FORCEFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].X << endl;
				OUT_FORCEFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].Y << endl;
				OUT_FORCEFILE << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[a1].Z << endl;
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
	  
				if ( CONTROLS.USE_HOOVER_THRMOSTAT) 
					printf(" %15s\n", "Econs/N");
				else 
					printf("\n");
			
				printf("%8s %9s %15s %15s %15s %15s %15s", " ", "(fs)", "(kcal/mol)", "(kcal/mol)", "(kcal/mol)", "(K)", "(GPa)");
			
				if ( CONTROLS.USE_HOOVER_THRMOSTAT ) 
					printf(" %15s\n", "(kcal/mol)");
				else 
					printf("\n");

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
			// Commented out b/c i need to define Ptot...
			printf("%8d %9.2f %15.7f %15.7f %15.7f %15.1f %15.3f", A+1, (A+1)*CONTROLS.DELTA_T_FS, Ktot/SYSTEM.ATOMS,SYSTEM.TOT_POT_ENER/SYSTEM.ATOMS,(Ktot+SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS,SYSTEM.TEMPERATURE, SYSTEM.PRESSURE);

			if ( CONTROLS.USE_HOOVER_THRMOSTAT ) 
				printf("%15.7f\n", (Ktot + SYSTEM.TOT_POT_ENER + 0.5 * THERMOSTAT.VISCO * THERMOSTAT.VISCO * THERMOSTAT.CHRG + THERMOSTAT.N_DOF * Kb * CONTROLS.TEMPERATURE * THERMOSTAT.COORD) / SYSTEM.ATOMS);
			else 
				printf("%15.7f\n\n", (Ktot + SYSTEM.TOT_POT_ENER)/SYSTEM.ATOMS);

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

			// This is the part that is particular to H2O... 
			
			GENFILE << "O H" << endl;
			
			int iele;

			for (int a1=0; a1<SYSTEM.ATOMS; a1++) 
			{
				if      ( SYSTEM.ATOMTYPE[a1] == "O" ) 
					iele = 1;
				else if ( SYSTEM.ATOMTYPE[a1] == "H" ) 
					iele = 2;
				else    
				{
					cout << "Error: Option # FMTDFTB # only valid for H2O systems!" << endl; 
					exit(1);
				}
				
				GENFILE << right << setw(4) << a1+1 << " " << setw(2) << iele << " " 
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
		OUT_FORCEFILE.close();
	
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
	fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.BOXDIM.X << " ";
	fxyz << fixed << setw(8) << setprecision(5) << SYSTEM.BOXDIM.X << endl;;
	
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
  
}       


static void read_input(MD_JOB_CONTROL & CONTROLS) 				// UPDATED
{
	cout << endl << "Reading the simulation control input file..." << endl;
	
	bool   			FOUND_END = false;
	string 			LINE;
	string			TEMP_STR;
	stringstream	LINE_PARSER;
	
	
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
				cout << "	# PRNTFRC #: true" << endl;	
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




