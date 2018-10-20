// System headers

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iomanip>
#include <map>
#include <sstream>
#include <unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes
#include <algorithm> // For specials vector-related functions (i.e. permute)

// User-defined headers

#include "functions.h"
#include "util.h"
#include "Cheby.h"
#include "io_styles.h"
#include "A_Matrix.h"

using namespace std;

#ifndef VERBOSITY 
	#define VERBOSITY 1 
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif 

static void read_lsq_input(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST &TRIPS, CLUSTER_LIST & QUADS, map<string,int> & PAIR_MAP,
													 vector<int> &INT_PAIR_MAP, vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, NEIGHBORS & NEIGHBOR_LIST,
													 vector<int>& TMP_ATOMTYPEIDX, vector<string> &TMP_ATOMTYPE);

static void print_bond_stats(vector<PAIRS> &ATOM_PAIRS, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, bool use_3b_cheby, bool use_4b_cheby);

static void build_clusters(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST &TRIPS, CLUSTER_LIST & QUADS, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST, vector<int>& TMP_ATOMTYPEIDX, vector<string>& TMP_ATOMTYPE);
static void print_param_header(JOB_CONTROL &CONTROLS, vector<PAIRS> &ATOM_PAIRS,
															  CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS) ;
static void print_map_file(map<string,int> PAIR_MAP, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS) ;
static int process_frame(A_MAT &A_MATRIX, JOB_CONTROL &CONTROLS, FRAME &SYSTEM, vector<PAIRS> &ATOM_PAIRS,
												 map<string,int> &PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS &NEIGHBOR_LIST,
												 CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, vector<CHARGE_CONSTRAINT> &CHARGE_CONSTRAINTS,
												 int i, int istart)  ;

// Print the params.header file.
// Global variables declared as externs in functions.h, and declared in functions.C

string FULL_FILE_3B;		
string SCAN_FILE_3B;
string SCAN_FILE_2B;
ofstream BAD_CONFIGS_1;	// Configs where r_ij < r_cut,in 
ofstream BAD_CONFIGS_2;	// Configs where r_ij < r_cut,in +d_penalty

// Global variables declared as externs in functions.h, and declared in functions.C -- MPI calculations.   
 
int NPROCS;		// Number of processors
int RANK;		// Index of current processor

// For 4-body interactions, these are used for both the lsq and md parts:

ofstream BAD_CONFIGS;

int main(int argc, char* argv[])
{
	
	// Set up MPI if requested, otherwise, run in serial 

	#ifdef USE_MPI
		MPI_Init     (&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &NPROCS);
		MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
	#else
		NPROCS = 1;
		RANK   = 0;
	#endif
	
	#ifdef ENABLE_FP_EXCEPT
		enable_fp_exceptions();
	#endif
	
	//////////////////////////////////////////////////
	//
	// Job parameters: Declare and set defaults
	//
	//////////////////////////////////////////////////
	
	vector<PAIRS> 		ATOM_PAIRS;		// Will store relevant info regarding atom interaction pair types.. 
	CLUSTER_LIST 		TRIPS;			// Triplet interaction: generic cluster interface.
	CLUSTER_LIST 		QUADS;			// Quadruplet interaction: generic cluster interface.
	
	FRAME		SYSTEM ;			// Stores the trajectory information... box dimensions, atom types, coordinates, and forces
	NEIGHBORS        	NEIGHBOR_LIST;		// Declare the class that will handle the neighbor list
	map<string,int> PAIR_MAP;			// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
	
	vector<CHARGE_CONSTRAINT> CHARGE_CONSTRAINTS;	// Specifies how we constrain charge fitting

	JOB_CONTROL 	CONTROLS;			// Will hold job controls shared by both lsq and md
	CONTROLS.IS_LSQ            = true; 

	
	vector<int>	INT_PAIR_MAP;

	vector<int>	ATOM_TYPE_IDX;			// Used to construct the quadruplet type index
	vector<string>	ATOM_TYPE;			// Used to construct the quadruplet type index
	
	
	//////////////////////////////////////////////////
	//
	// Read and print input to screen
	//
	//////////////////////////////////////////////////
	
	if ( RANK == 0 ) 
		cout << endl << "Reading input file..." << endl;
		
	if ( RANK == 0 ) 
	{
		// Delete temporary files if they exist.
		system("rm -f A.[0-9]*.txt");
		system("rm -f b.[0-9]*.txt");
	}
	
	#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD) ;	// Make sure the files have been removed before they are re-created by other processes.
	#endif

	read_lsq_input(CONTROLS, ATOM_PAIRS, TRIPS, QUADS, PAIR_MAP, INT_PAIR_MAP,  CHARGE_CONSTRAINTS, NEIGHBOR_LIST, ATOM_TYPE_IDX, ATOM_TYPE);

	// Build many-body interaction clusters if necessary.
	build_clusters(CONTROLS, ATOM_PAIRS, TRIPS, QUADS, PAIR_MAP, NEIGHBOR_LIST, ATOM_TYPE_IDX, ATOM_TYPE);

	if ( RANK == 0 ) 
		cout << "...input file read successful: " << endl << endl;
	
	if ( NPROCS > CONTROLS.NFRAMES ) NPROCS = CONTROLS.NFRAMES ;	// Don't use unnecessary processors.


	//////////////////////////////////////////////////
	//
	// Set up force and force derivative vectors
	//
	//////////////////////////////////////////////////	
		 		

	A_MAT A_MATRIX ; // Declare and initialize A-matrix object
	

	// Figure out necessary dimensions for the force/force derivative vectors
	
	CONTROLS.TOT_SNUM     = 0; 
	CONTROLS.NUM_3B_CHEBY = 0;
	CONTROLS.NUM_4B_CHEBY = 0;

	for (int i=0; i<ATOM_PAIRS.size(); i++)
	{
		if ( (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) || (ATOM_PAIRS[i].PAIRTYP == "DFTBPOLY") )	
		{
			ATOM_PAIRS[i].SNUM = CONTROLS.CHEBY_ORDER;
			
			if (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" )
			{
				ATOM_PAIRS[i].SNUM_3B_CHEBY = CONTROLS.CHEBY_3B_ORDER;
				ATOM_PAIRS[i].SNUM_4B_CHEBY = CONTROLS.CHEBY_4B_ORDER;
				ATOM_PAIRS[i].CHEBY_TYPE    = CONTROLS.CHEBY_TYPE;
			}
				
		}
		
		else if (ATOM_PAIRS[i].PAIRTYP == "INVRSE_R") 
			ATOM_PAIRS[i].SNUM = CONTROLS.INVR_PARAMS;

		else // Spline
			ATOM_PAIRS[i].SNUM = (2+floor((ATOM_PAIRS[i].S_MAXIM - ATOM_PAIRS[i].S_MINIM)/ATOM_PAIRS[i].S_DELTA))*2; //2 is for p0/m0/p1/m1.. 
	
		CONTROLS.TOT_SNUM += ATOM_PAIRS[i].SNUM;
	}

	
	if ( RANK == 0 ) 
		cout << "The number of two-body non-coulomb parameters is: " << CONTROLS.TOT_SNUM <<  endl;

	
	if ( ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV" && CONTROLS.CHEBY_3B_ORDER  > 0 ) // All atoms types must use same potential, so check if 3b order is greater than zero 
	{
		CONTROLS.NUM_3B_CHEBY = 0;
		
		for(int i=0; i< TRIPS.VEC.size(); i++)
			CONTROLS.NUM_3B_CHEBY += TRIPS.VEC[i].N_TRUE_ALLOWED_POWERS;

		if ( RANK == 0 ) 
			cout << "The number of three-body Chebyshev parameters is: " << CONTROLS.NUM_3B_CHEBY << endl;

	}
	
	if ( ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV" && CONTROLS.CHEBY_4B_ORDER  > 0 ) // All atoms types must use same potential, so check if 3b order is greater than zero 
	{
		CONTROLS.NUM_4B_CHEBY = 0;
		
		for(int i=0; i< QUADS.VEC.size(); i++)
			CONTROLS.NUM_4B_CHEBY += QUADS.VEC[i].N_TRUE_ALLOWED_POWERS;

		if ( RANK == 0 ) 
			cout << "The number of four-body  Chebyshev parameters is: " << CONTROLS.NUM_4B_CHEBY << endl;
	}	

	CONTROLS.TOT_SHORT_RANGE = CONTROLS.TOT_SNUM + CONTROLS.NUM_3B_CHEBY + CONTROLS.NUM_4B_CHEBY;	

	//////////////////////////////////////////////////
	//
	// Setup variables for reading the trajectory file(s)
	//
	//////////////////////////////////////////////////
	
	cout.precision(16);				// WE SHOULD AT LEAST MOVE THIS TO SOMEWHERE MORE REASONABLE.. LIKE THE SECTION WHERE THE OUTPUT IS ACTUALLY PRINTED
	
	int istart, iend;
	
	divide_atoms(istart, iend, CONTROLS.NFRAMES);	// Each processor only calculates certain frames.  Calculate which frames to use.
	
	if((CONTROLS.FIT_STRESS  || CONTROLS.FIT_STRESS_ALL) && CONTROLS.NSTRESS == -1)
		CONTROLS.NSTRESS = CONTROLS.NFRAMES;
		
	if((CONTROLS.FIT_ENER || CONTROLS.FIT_ENER_PER_ATOM) && CONTROLS.NENER == -1)
		CONTROLS.NENER = CONTROLS.NFRAMES;		
	
	CONTROLS.FIT_ENER_EVER = CONTROLS.FIT_ENER || CONTROLS.FIT_ENER_PER_ATOM ;	// Is energy ever fit ?

	A_MATRIX.OPEN_FILES() ;

	int total_forces = 0 ;
	
	//////////////////////////////////////////////////
	//	
	// Begin processing the trajecory
	//
	//////////////////////////////////////////////////
	 
	int FILE_IDX = 0;		// Index of traj file in CONTROLS.INFILE vector
	ifstream TRAJ_INPUT;

	if (CONTROLS.INFILE.size() == 1)
		CONTROLS.INFILE_FRAMES.push_back(CONTROLS.NFRAMES);

	OPEN_TRAJFILE(TRAJ_INPUT, CONTROLS.INFILE, FILE_IDX);
	
	if ( RANK == 0 ) 
		cout << endl << "Succefully opened the trajectory file..." << endl;

	if ( RANK == 0 ) 
		cout << "Setting up the matrices for A, Coulomb forces, and overbonding..." << endl;
		
	int OFFSET = CONTROLS.INFILE_FRAMES[FILE_IDX];

	for (int i=0; i<CONTROLS.NFRAMES; i++)
	{
		if( (i+1) > OFFSET) // We've reached the end of file. Move on to the next one (if there is a next one)
		{
			//if ( RANK == 0 ) // Debug statement
			//	cout << "Closing file " << CONTROLS.INFILE[FILE_IDX] << " and opening file " << CONTROLS.INFILE[FILE_IDX+1]  << " to read frame " << i+1 << endl;

			FILE_IDX++;
			OPEN_TRAJFILE(TRAJ_INPUT, CONTROLS.INFILE, FILE_IDX); // Closes current file (if open), opens file FILE_IDX
			
			OFFSET += CONTROLS.INFILE_FRAMES[FILE_IDX];			
		}
		
		SYSTEM.READ_XYZF(TRAJ_INPUT, CONTROLS, ATOM_PAIRS, ATOM_TYPE, i) ;

		if ( i >= istart && i <= iend )
			total_forces += process_frame(A_MATRIX, CONTROLS, SYSTEM, ATOM_PAIRS, PAIR_MAP,INT_PAIR_MAP, NEIGHBOR_LIST, TRIPS, QUADS, CHARGE_CONSTRAINTS, i, istart) ;
	}

	// Add the number of force entries in the A matrix across all processes.
	
	int total_forces_all = 0 ;
	
	#ifdef USE_MPI
		MPI_Reduce(&total_forces, &total_forces_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) ;
	#else
		total_forces_all = total_forces ;
	#endif
	
	// What's below is only used to determine whether to add a col of ones to 
	// the A-matrix (i.e. "are ANY energies included in the fit?")

	if ( RANK == 0 ) 
	{
		cout << "...matrix population complete: "  << endl << endl;
		cout << "Printing matrices..." << endl;
		cout << "	...A matrix length (forces): " << total_forces_all << endl << endl ;

		if (CONTROLS.FIT_ENER_EVER)
		{
			cout << "Per-atom type energies will be provided at the end of the parameter file." << endl;
			cout << "Values will correspond to the following atom types, in order: " << endl;
			
			for(int a=0; a<A_MATRIX.NO_ATOM_TYPES; a++)
				cout << "	" << A_MATRIX.ATOM_TYPES[a] << " " << A_MATRIX.NO_ATOMS_OF_TYPE[a] << endl;		
		}
	}

	A_MATRIX.PRINT_CONSTRAINTS(CONTROLS, CHARGE_CONSTRAINTS, ATOM_PAIRS.size() ) ;
	A_MATRIX.CLEANUP_FILES() ;
		
	
 	 //////////////////////////////////////////////////
	//
	// Print out the params file header
	//
	//////////////////////////////////////////////////	   

	print_param_header(CONTROLS, ATOM_PAIRS, TRIPS, QUADS) ;
	
	//////////////////////////////////////////////////
	//
	// Print out the pair/triplet type  map file
	//
	//////////////////////////////////////////////////	  	

	print_map_file(PAIR_MAP, TRIPS, QUADS) ;

	//////////////////////////////////////////////////
	//
	// Print out bonding statistics.
	//
	//////////////////////////////////////////////////	  

	print_bond_stats(ATOM_PAIRS, TRIPS, QUADS, CONTROLS.USE_3B_CHEBY, CONTROLS.USE_4B_CHEBY);

	  
	return 0;		  
}




	//////////////////////////////////////////////////
	//
    // Function definitions
	//
	//////////////////////////////////////////////////



// Read program input from the file "splines_ls.in".
static void read_lsq_input(JOB_CONTROL & CONTROLS, 
                           vector<PAIRS> & ATOM_PAIRS, 
                           CLUSTER_LIST &TRIPS, CLUSTER_LIST & QUADS, map<string,int> & PAIR_MAP, 
                           vector<int> &INT_PAIR_MAP,
                           vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, 
                           NEIGHBORS & NEIGHBOR_LIST,
                           vector<int>& ATOM_TYPE_IDX,
                           vector<string>& ATOM_TYPE)
{
	bool   FOUND_END = false;
	string LINE;
	string TEMP_STR;
	PAIRS  TEMP_PAIR;
	int    TEMP_INT;
	string TEMP_TYPE;
	int    NPAIR;
	int    NTRIP;
	int    NQUAD;
	double SUM_OF_CHARGES = 0;
	stringstream	STREAM_PARSER;
	int ntokens ;
	vector<string> tokens ;
	
	// Set some defaults
	
	double TMP_CHEBY_RANGE_LOW  = -1;
	double TMP_CHEBY_RANGE_HIGH =  1;
	
	NEIGHBOR_LIST.USE	   = true;
	
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

		//if ( RANK == 0 ) cout << "CURRENT LINE IS: " << LINE << endl ;

		if(LINE.find("# ENDFILE #") != string::npos)
		{
			// Run a few checks to make sure logic is correct
		 
			if(CONTROLS.IF_SUBTRACT_COORD && CONTROLS.FIT_POVER)
			{
				cout << "LOGIC ERROR: Problem with code logic. Both fit_pover and ifsubtract_coord cannot be true." << endl;
				cout << "             if_subtract_coord should only be true if overbonding parameters have been " << endl;
				cout << "             specified, and FITPOVR set false." << endl;
				exit(0);
			}
		
			if(CONTROLS.IF_SUBTRACT_COUL && CONTROLS.FIT_COUL)
			{
				cout << "LOGIC ERROR: Problem with code logic. Both fit_coul and ifsubtract_coul cannot be true." << endl;
				cout << "             ifsubtract_coul should only be true if non-zero charges have been specified " << endl;
				cout << "             and FITCOUL set false." << endl;
				exit(0);
			}

			if (CONTROLS.IF_SUBTRACT_COORD && RANK == 0)
			{
				cout << "Special feature: " << endl;
				cout << " Will subtract contributions from user-specified overbonding " << endl;
				cout << " parameters before generating A-matrix." << endl << endl;	
			}					

			if(CONTROLS.USE_PARTIAL_CHARGES && RANK == 0)
			{
				cout << "Special feature: " << endl;
				cout << " Will subtract contributions stemming from user-specified " << endl;
				cout << " charges before generating A-matrix" << endl << endl;	
			}	

			if(CONTROLS.WRAP_COORDS && CONTROLS.N_LAYERS > 0 )
			{
				if ( RANK == 0 ) 
					cout << "WARNING: Coordinate wrapping not supported for ghost atom use. Turning option off" << endl;
				CONTROLS.WRAP_COORDS = false;
			}

			if (CONTROLS.USE_PARTIAL_CHARGES || CONTROLS.FIT_COUL)
				CONTROLS.CALL_EWALD = true;	
				
			if((CONTROLS.FIT_STRESS || CONTROLS.FIT_STRESS_ALL) && CONTROLS.CALL_EWALD)
			{
				EXIT_MSG("ERROR: Inclusion of stress tensors currently not compatible with use of ZCalc_Ewald_Deriv.") ;
			}

			FOUND_END = true;
			break;
				
		}
		
		else if(LINE.find("# WRAPTRJ #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.WRAP_COORDS = bool(LINE.data());
			
			if ( RANK == 0 ) 
				cout << "	# WRAPTRJ #: " << CONTROLS.WRAP_COORDS << endl;	
		}
		
		else if(LINE.find("# TRJFILE #") != string::npos)
		{
			getline(cin, LINE); 
			
			vector<string> tokens;
			
			parse_space(LINE,tokens);
			
			if(tokens.size() == 1) // Then we're reading in a single trajectory file
			{
				CONTROLS.INFILE.push_back(tokens[0]);
				
				if ( RANK == 0 ) 
					cout << "	# TRJFILE #: " << CONTROLS.INFILE[0] << endl;	
			}	
			else if (tokens.size() == 2)
			{
				if (tokens[0] != "MULTI")
				{
					cout << "ERROR: Unknown syntax encoutered reading \"# TRJFILE #\". Exiting" << endl;
					exit_run(0);
				}
				
				else
				{
					if ( RANK == 0 ) 
						cout << "	# TRJFILE #: MULTI " << tokens[1] << endl;
						
					// Read the files (later, will add option to skip first N frames, and to process every M frames)
					
					ifstream MULTI;
					MULTI.open(tokens[1]);
					
					if (!MULTI.is_open())
					{
						cout << "ERROR: Cannot open MULTI file: " << tokens[1] << ". Exiting." << endl;
						exit_run(0);
					}
					
					if ( RANK == 0 ) 
						cout << "		...Will read from the following files: " << endl;
					
					getline(MULTI, LINE);
					
					parse_space(LINE,tokens);
					
					if(tokens.size() != 1)
					{
						cout << "ERROR: Expected to read <nfiles>, got: " << LINE << endl;
						exit_run(0);
					}
					
					int NFILES = stoi(tokens[0]);
					
					for (int i=0; i<NFILES; i++)
					{
						getline(MULTI, LINE);

						parse_space(LINE,tokens);
						
						if(tokens.size() != 2)
						{
							cout << "ERROR: Expected to read <nframes> <filename>, got: " << LINE << endl;
							exit_run(0);
						}
											
						
						CONTROLS.INFILE_FRAMES.push_back(stoi(tokens[0]));
						CONTROLS.INFILE       .push_back(      tokens[1]);
						
						if ( RANK == 0 ) 
							cout << "		-   " 
							     << CONTROLS.INFILE[i]        << " with " 
							     << CONTROLS.INFILE_FRAMES[i] << " frames." << endl;
					}
					
					MULTI.close();
				}
			}
		}			
		
		else if(LINE.find("# NFRAMES #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.NFRAMES = int(atof(LINE.data()));

			if ( RANK == 0 ) 
				cout << "	# NFRAMES #: " << CONTROLS.NFRAMES << endl;			
		}
		
		else if(LINE.find("# NLAYERS #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.N_LAYERS = int(atof(LINE.data()));
			
			if ( RANK == 0 ) 	
				cout << "	# NLAYERS #: " << CONTROLS.N_LAYERS << endl;
		}
		
		else if(LINE.find("# FITCOUL #") != string::npos)
		{
			cin >> LINE; cin.ignore();

			if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
				CONTROLS.FIT_COUL = true;
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				CONTROLS.FIT_COUL = false;
			else
			{
				cout << endl << "ERROR: # FITCOUL # must be specified as true or false." << endl;
				exit(1);	
			}	
			
			if ( RANK == 0 ) 
			{
				cout << "	# FITCOUL #: ";
			
				if (CONTROLS.FIT_COUL)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
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
		}
		
		else if(LINE.find("# FITEATM #") != string::npos)
		{
			cin >> LINE; 
			
			if (LINE=="first"  || LINE=="First"  || LINE=="FIRST")
			{
				CONTROLS.FIT_ENER_PER_ATOM = true;
				cin >> CONTROLS.NENER;
				cin.ignore();
			}
			else
			{

				if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
					CONTROLS.FIT_ENER_PER_ATOM = true;
				else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
					CONTROLS.FIT_ENER_PER_ATOM = false;
				else
				{
					cout << endl << "ERROR: # FITEATM # must be specified as true or false." << endl;
					exit(1);	
				}
				
				cin.ignore();	
			}
			
			if ( RANK == 0 )
			{
				cout << "	# FITEATM #: ";
			
				if (CONTROLS.FIT_ENER_PER_ATOM)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
					
				if(CONTROLS.NENER>0)
					cout << "    			 ...will only fit energies for first " << CONTROLS.NENER << " frames." << endl;
			}
		}		

		else if(LINE.find("# FITPOVR #") != string::npos)
		{
			cin >> LINE; cin.ignore();

			if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
				CONTROLS.FIT_POVER = true;
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				CONTROLS.FIT_POVER = false;
			else
			{
				cout << endl << "ERROR: # FITPOVR # must be specified as true or false." << endl;
				exit(1);	
			}	
			
			if ( RANK == 0 ) 
			{
				cout << "	# FITPOVR #: ";
			
				if (CONTROLS.FIT_POVER)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
			}
		}
		
		else if(LINE.find("# PAIRTYP #") != string::npos)
		{
			cin >> LINE; // cin.ignore();
			
			TEMP_TYPE = LINE;

			if      (LINE != "SPLINE" && LINE != "CHEBYSHEV" && LINE != "DFTBPOLY" && LINE != "INVRSE_R") // These are not supported:  && LINE != "LJ" && LINE != "STILLIN")
			{
				cout << endl;
				cout << "ERROR: Unrecognized pair type. Acceptable options are:" << endl;
				cout << "SPLINE"   << endl;
				cout << "CHEBYSHEV"    << endl;
				cout << "DFTBPOLY" << endl;
				cout << "INVRSE_R" << endl;
//				cout << "LJ"       << endl;
//				cout << "STILLIN"  << endl;
				exit(1);
				
			}
			
			if ( RANK == 0 )
			{
				cout << "	# PAIRTYP #: " << LINE;
				
				if(TEMP_TYPE != "DFTBPOLY")
					cout << " ....NOTE: Forces reported in units of kcal/(mol.A), potential energy in kcal/mol." << endl;
				else
					cout << " ....NOTE: Forces reported in atomic units." << endl;
			}
				
			if(TEMP_TYPE == "INVRSE_R")
			{
				cin >> LINE;
				//cin.ignore();
				CONTROLS.INVR_PARAMS = stoi(LINE);
					
				if ( RANK == 0 ) 
					cout << "	             " << "Will use the following number of inverse-r params: " << CONTROLS.INVR_PARAMS << endl;					
			}

			if(TEMP_TYPE == "DFTBPOLY" || TEMP_TYPE == "CHEBYSHEV")
			{
				getline(cin,LINE);
				
				STREAM_PARSER.str(LINE);
				
				STREAM_PARSER >> CONTROLS.CHEBY_ORDER;				
				
				if ( RANK == 0 ) 
					cout << "	             " << "Will use 2-body order: " << CONTROLS.CHEBY_ORDER << endl;
				
				if(TEMP_TYPE == "CHEBYSHEV")
				{
					STREAM_PARSER >> CONTROLS.CHEBY_3B_ORDER;

					if ( RANK == 0 ) 
						cout << "	             " << "Will use 3-body order: " << CONTROLS.CHEBY_3B_ORDER << endl;
						
					if(CONTROLS.CHEBY_3B_ORDER>0)
						CONTROLS.USE_3B_CHEBY = true;
					
					STREAM_PARSER >> CONTROLS.CHEBY_4B_ORDER;

					if ( RANK == 0 ) 
						cout << "	             " << "Will use 4-body order: " << CONTROLS.CHEBY_4B_ORDER << endl;
						
					if(CONTROLS.CHEBY_4B_ORDER>0)
						CONTROLS.USE_4B_CHEBY = true;
					
					if (STREAM_PARSER >>  TMP_CHEBY_RANGE_LOW)
					{
						if(TMP_CHEBY_RANGE_LOW < -1.0 || TMP_CHEBY_RANGE_LOW > +1.0 )
						{
							cout << "ERROR: TMP_CHEBY_RANGE_LOW must be betwee -1 and 1" << endl;
							exit(0);
						}
					}
					if (STREAM_PARSER >>  TMP_CHEBY_RANGE_HIGH)
					{
						if(TMP_CHEBY_RANGE_HIGH < -1.0 || TMP_CHEBY_RANGE_HIGH > +1.0 )
						{
							cout << "ERROR: TMP_CHEBY_RANGE_LOW must be betwee -1 and 1" << endl;
							exit(0);
						}
					}
					if(TMP_CHEBY_RANGE_HIGH < TMP_CHEBY_RANGE_LOW)
					{
						cout << "ERROR: TMP_CHEBY_RANGE_HIGH must be greater than TMP_CHEBY_RANGE_LOW" << endl;
						exit(0);
					}					
					
					if ( RANK == 0 ) 
						cout << "	             " << "Will transform Chebyshev pair distances to range " << TMP_CHEBY_RANGE_LOW << " to " << TMP_CHEBY_RANGE_HIGH << endl;
				}
				else
					
				cin.ignore();
				
				STREAM_PARSER.str("");
				STREAM_PARSER.clear();		
				
			}
			else
				 cin.ignore();
			
			// Do some error checks
	
			if(CONTROLS.USE_3B_CHEBY && CONTROLS.FIT_POVER)
			{
				cout << "ERROR: Overbonding is not compatible with 3-body Chebyshev potentials." << endl;
				cout << "       Set # FITPOVR # false." << endl;
				exit(0);				
			}

			
		}

		else if( LINE.find("# CHBTYPE #") != string::npos )
		{

		  	cin >> LINE ;
			
			if ( TEMP_TYPE == "CHEBYSHEV" ) 
			{
				vector<string> tokens ;
				
				if ( parse_space(LINE,tokens) >= 1 ) 
				{
					CONTROLS.CHEBY_TYPE = Cheby::get_trans_type(tokens[0]);
			
					if ( RANK == 0 ) 
						cout << "	# CHBTYPE #: " << tokens[0] << endl;	
				}		
				else 
					EXIT_MSG("BAD CHBTYPE" + LINE) ;
			} 
				
			else if ( RANK == 0 ) 
			{
				cout << "Warning: CHBTYPE given for a non-Chebyshev pair type (ignored)" ;
			}

		}		
		
		/////////////////////////////////////////////////////////////////////
		// Read the topology part. For now, ignoring index and atom types.. 
		// Assuming given as OO, OH, HH, as code expects... 
		// will need to be fixed later.
		/////////////////////////////////////////////////////////////////////
		
		else if(LINE.find("EXCLUDE 3B INTERACTION:")!= string::npos)
		{
		  TRIPS.read_exclude(cin, LINE);
		}	

		else if(LINE.find("EXCLUDE 4B INTERACTION:")!= string::npos)
		{
		  QUADS.read_exclude(cin, LINE);
		}	
		
		else if(LINE.find("# NATMTYP #")!= string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.NATMTYP = int(atof(LINE.data()));
			
			if ( RANK == 0 ) 
				cout << "	# NATMTYP #: " << CONTROLS.NATMTYP << endl;	
			
			// Set up pairs
			
			NPAIR = CONTROLS.NATMTYP*(CONTROLS.NATMTYP+1)/2;
			ATOM_PAIRS.resize(NPAIR);
			
			// Set the default cheby range and the fcut type
			
			if (TEMP_TYPE == "CHEBYSHEV")
			{
				for(int i=0; i<NPAIR; i++)
				{
					ATOM_PAIRS[i].CHEBY_RANGE_LOW   = TMP_CHEBY_RANGE_LOW;
					ATOM_PAIRS[i].CHEBY_RANGE_HIGH  = TMP_CHEBY_RANGE_HIGH;
					ATOM_PAIRS[i].FORCE_CUTOFF.TYPE = FCUT_TYPE::CUBIC;
					ATOM_PAIRS[i].PAIRTYP           = TEMP_TYPE;
				}
			}	
			
			// Set up triplets
			
			NTRIP = factorial(CONTROLS.NATMTYP+3-1)/factorial(3)/factorial(CONTROLS.NATMTYP-1);
			TRIPS.allocate(NTRIP, 3, ATOM_PAIRS);
			
			// Set up quadruplets
			
			NQUAD = factorial(CONTROLS.NATMTYP+4-1)/factorial(4)/factorial(CONTROLS.NATMTYP-1);
			QUADS.allocate(NQUAD, 4, ATOM_PAIRS);
		}

		else if(LINE.find("# TYPEIDX #")!= string::npos)
		{
			
			if ( RANK == 0 ) 
				cout << "	# TYPEIDX #    # ATM_TYP #    # ATMCHRG #    # ATMMASS #" << endl;
			
			// Figure out the number of non-unique pairs

			TEMP_INT = CONTROLS.NATMTYP*CONTROLS.NATMTYP;
			
			SUM_OF_CHARGES = 0;
			
			ATOM_TYPE_IDX.resize(CONTROLS.NATMTYP);
			ATOM_TYPE   .resize(CONTROLS.NATMTYP);
			
			for(int i=0; i<CONTROLS.NATMTYP; i++)
			{
				
				// Set the first atom pair types to be of type OO, HH, CC, etc...
				
				ATOM_PAIRS[i].PAIRTYP    = TEMP_TYPE;
				ATOM_PAIRS[i].PAIRIDX	 = i; 
				ATOM_PAIRS[i].CHEBY_TYPE = CONTROLS.CHEBY_TYPE;
				
				cin >> LINE;
				ATOM_TYPE_IDX[i] = stoi(LINE) - 1;
				ATOM_PAIRS[i].ATM1TYPE_IDX = ATOM_TYPE_IDX[i];

				if ( ATOM_TYPE_IDX[i] < 0 || ATOM_TYPE_IDX[i] >= CONTROLS.NATMTYP )
				  EXIT_MSG("Bad atom index: " + LINE );

				cin >> ATOM_PAIRS[i].ATM1TYP >> LINE;
				
				ATOM_TYPE[i]    = ATOM_PAIRS[i].ATM1TYP;
				
				if(!CONTROLS.FIT_COUL)
				{
					ATOM_PAIRS[i].ATM1CHG = double(atof(LINE.data()));
				}
				else
				{
					ATOM_PAIRS[i].CHRGSGN = LINE;
					ATOM_PAIRS[i].ATM1CHG = 0.0;
				}

				SUM_OF_CHARGES += abs(ATOM_PAIRS[i].ATM1CHG);

				cin >> LINE; cin.ignore();
				ATOM_PAIRS[i].ATM1MAS = double(atof(LINE.data()));
				
				ATOM_PAIRS[i].ATM2TYP = ATOM_PAIRS[i].ATM1TYP;
				ATOM_PAIRS[i].ATM2TYPE_IDX = ATOM_PAIRS[i].ATM1TYPE_IDX;

				ATOM_PAIRS[i].ATM2CHG = ATOM_PAIRS[i].ATM1CHG;
				ATOM_PAIRS[i].ATM2MAS = ATOM_PAIRS[i].ATM1MAS;
				
				ATOM_PAIRS[i].PRPR_NM =      ATOM_PAIRS[i].ATM1TYP;
				ATOM_PAIRS[i].PRPR_NM.append(ATOM_PAIRS[i].ATM2TYP);

				if ( RANK == 0 ) 
				{
					cout << " 	" << setw(15) << left << i+1 
						 << setw(15) << left << ATOM_PAIRS[i].ATM1TYP;
					if(!CONTROLS.FIT_COUL) 
						cout << setw(15) << left << ATOM_PAIRS[i].ATM1CHG;
					else
						cout << ATOM_PAIRS[i].CHRGSGN << "		";
					
					cout << setw(15) << left << ATOM_PAIRS[i].ATM1MAS << endl;
				}
			}
			
			if(!CONTROLS.FIT_COUL)
			{
				if (SUM_OF_CHARGES>0)
				{
					CONTROLS.USE_PARTIAL_CHARGES = true;
					CONTROLS.IF_SUBTRACT_COUL    = true;
					
				}
				else
					CONTROLS.USE_PARTIAL_CHARGES = false;				
			}
			
			// Set up all possible unique pair types
			
			TEMP_INT = CONTROLS.NATMTYP;
			
			for(int i=0; i<CONTROLS.NATMTYP; i++)
			{
				for(int j=i+1; j<CONTROLS.NATMTYP; j++)
				{	
					ATOM_PAIRS[TEMP_INT].PAIRTYP    = ATOM_PAIRS[i].PAIRTYP;
					ATOM_PAIRS[TEMP_INT].PAIRIDX    = TEMP_INT; 
					ATOM_PAIRS[TEMP_INT].CHEBY_TYPE = ATOM_PAIRS[i].CHEBY_TYPE;
														
					ATOM_PAIRS[TEMP_INT].ATM1TYP = ATOM_PAIRS[i].ATM1TYP;
					ATOM_PAIRS[TEMP_INT].ATM1TYPE_IDX = ATOM_PAIRS[i].ATM1TYPE_IDX;

					ATOM_PAIRS[TEMP_INT].ATM2TYP = ATOM_PAIRS[j].ATM1TYP;
					ATOM_PAIRS[TEMP_INT].ATM2TYPE_IDX = ATOM_PAIRS[j].ATM1TYPE_IDX;

					ATOM_PAIRS[TEMP_INT].ATM1CHG = ATOM_PAIRS[i].ATM1CHG;
					ATOM_PAIRS[TEMP_INT].ATM2CHG = ATOM_PAIRS[j].ATM1CHG;
					
					ATOM_PAIRS[TEMP_INT].ATM1MAS = ATOM_PAIRS[i].ATM1MAS;
					ATOM_PAIRS[TEMP_INT].ATM2MAS = ATOM_PAIRS[j].ATM1MAS;	
					
					ATOM_PAIRS[TEMP_INT].PRPR_NM =      ATOM_PAIRS[TEMP_INT].ATM1TYP;
					ATOM_PAIRS[TEMP_INT].PRPR_NM.append(ATOM_PAIRS[TEMP_INT].ATM2TYP);											
					
					TEMP_INT++;
				}	
			}

			// Set up the maps to account for non-unique pairs
			
			for(int i=0; i<CONTROLS.NATMTYP; i++)
			{
				for(int j=0; j<CONTROLS.NATMTYP; j++)
				{
					TEMP_STR = ATOM_PAIRS[i].ATM1TYP;
					TEMP_STR.append(ATOM_PAIRS[j].ATM1TYP);
					
					for(int k=0;k<ATOM_PAIRS.size(); k++)
					{
						if((ATOM_PAIRS[i].ATM1TYP == ATOM_PAIRS[k].ATM1TYP && ATOM_PAIRS[j].ATM1TYP == ATOM_PAIRS[k].ATM2TYP)
						 ||(ATOM_PAIRS[j].ATM1TYP == ATOM_PAIRS[k].ATM1TYP && ATOM_PAIRS[i].ATM1TYP == ATOM_PAIRS[k].ATM2TYP))
						{
							TEMP_STR = ATOM_PAIRS[i].ATM1TYP;
							TEMP_STR.append(ATOM_PAIRS[j].ATM2TYP);
							PAIR_MAP.insert(make_pair(TEMP_STR,k));	// Maps the true pair index, k, to the string formed by joining the two atom types.
						}
					}
				}
			}
			
			//cout << "Made the following pairs: " << endl;
			//for(map<string,int>::iterator i=PAIR_MAP.begin(); i!=PAIR_MAP.end(); i++)
			//cout << i->second << " " << i->first << endl;			
					
			if ( RANK == 0 ) 
			{
				cout << endl;
				cout << "	The following unique pair types have been identified:" << endl;
				for(int i=0;i<NPAIR; i++)
					cout << "		" << ATOM_PAIRS[i].PAIRIDX << "  " << ATOM_PAIRS[i].ATM1TYP << " " << ATOM_PAIRS[i].ATM2TYP << endl;
			}
					

		}

		else if(LINE.find("# PAIRIDX #")!= string::npos) // Read the topology part. For now, ignoring index and atom types..
		{
			for(int i=0; i<NPAIR; i++)
			{
				cin >> LINE >> TEMP_PAIR.ATM1TYP >> TEMP_PAIR.ATM2TYP;
				
				TEMP_STR = TEMP_PAIR.ATM1TYP;
				TEMP_STR.append(TEMP_PAIR.ATM2TYP);
				TEMP_INT = PAIR_MAP[TEMP_STR];
				
				ATOM_PAIRS[TEMP_INT].PAIRIDX = TEMP_INT;
				
				ATOM_PAIRS[TEMP_INT].CUBIC_SCALE = 1.0;	// Set the default value

				cin >> LINE;
				ATOM_PAIRS[TEMP_INT].S_MINIM = double(atof(LINE.data()));
					
				cin >> LINE;
				ATOM_PAIRS[TEMP_INT].S_MAXIM = double(atof(LINE.data()));
					
				cin >> LINE; cin.ignore();
				ATOM_PAIRS[TEMP_INT].S_DELTA = double(atof(LINE.data()));
				
				cin >> LINE; cin.ignore();
				ATOM_PAIRS[TEMP_INT].LAMBDA = double(atof(LINE.data()));	
				
				ATOM_PAIRS[TEMP_INT].MIN_FOUND_DIST = 	1.0e10;	// Set an initial minimum distance	
				
				if( ATOM_PAIRS[TEMP_INT].S_MAXIM > NEIGHBOR_LIST.MAX_CUTOFF)
				{
					NEIGHBOR_LIST.MAX_CUTOFF    = ATOM_PAIRS[TEMP_INT].S_MAXIM;
					NEIGHBOR_LIST.MAX_CUTOFF_3B = ATOM_PAIRS[TEMP_INT].S_MAXIM;
					NEIGHBOR_LIST.MAX_CUTOFF_4B = ATOM_PAIRS[TEMP_INT].S_MAXIM;
				}	
				
				cin >> LINE; cin.ignore();
				
				if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
				{
					cin >> ATOM_PAIRS[TEMP_INT].OVER_TO_ATM;
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[0];
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[1];
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[2];
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[3];
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[4];
					
					ATOM_PAIRS[TEMP_INT].USE_OVRPRMS = true;
				
				}
				else
				{
					ATOM_PAIRS[TEMP_INT].USE_OVRPRMS = false;
				}		
				
				// Check if 3B bin widths are specified
				
				if(!cin.eof())	
				{
					cin >> ATOM_PAIRS[TEMP_INT].NBINS[0];
					cin >> ATOM_PAIRS[TEMP_INT].NBINS[1];
					cin >> ATOM_PAIRS[TEMP_INT].NBINS[2];
				}
				else
				{
					ATOM_PAIRS[TEMP_INT].NBINS[0] = 0;
					ATOM_PAIRS[TEMP_INT].NBINS[1] = 0;
					ATOM_PAIRS[TEMP_INT].NBINS[2] = 0;
				}	
			}
		
			if ( RANK == 0 ) 
			{
				cout << "	# PAIRIDX #     ";
				cout << "# ATM_TY1 #     ";
				cout << "# ATM_TY1 #     ";
				cout << "# S_MINIM #     ";
				cout << "# S_MAXIM #     ";
				cout << "# S_DELTA #     ";
				cout << "# MORSE_LAMBDA #";	
				cout << " # USEOVRP #     " << endl;
			}
				
			bool PRINT_OVR = false;
			
			for(int i=0; i<NPAIR; i++)
			{
				if ( RANK == 0 ) 
				{
					cout << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
						 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
						 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
						 << setw(16) << left << ATOM_PAIRS[i].S_MINIM
						 << setw(16) << left << ATOM_PAIRS[i].S_MAXIM							 
						 << setw(16) << left << ATOM_PAIRS[i].S_DELTA
						 << setw(16) << left << ATOM_PAIRS[i].LAMBDA << " " 
						 << setw(16) << left << ATOM_PAIRS[i].USE_OVRPRMS << endl;
				}
					 
				if(ATOM_PAIRS[i].USE_OVRPRMS)
					PRINT_OVR = true;
			}
		
			// If overbonding parameters are provided, and they are not requested to be fit, 
			// subtract thier contribution before generating A matrix
			
			if(PRINT_OVR && !CONTROLS.FIT_POVER)
				CONTROLS.IF_SUBTRACT_COORD = true;

			if(PRINT_OVR)
			{
				if ( RANK == 0 ) 
				{
					cout << "	# PAIRIDX #     ";
					cout << "# ATM_TY1 #     ";
					cout << "# ATM_TY1 #     ";				
					cout << "# P_OVERB #     ";
					cout << "# R_0_VAL #     ";
					cout << "# P_1_VAL #     ";
					cout << "# P_2_VAL #     ";
					cout << "# LAMBDA6 #" << endl;						
				}

				for(int i=0; i<NPAIR; i++)
				{					
					if ( RANK == 0 ) 
					{
						cout << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
							 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
							 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[0] 
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[1]
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[2] 
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[3]
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[4] << endl;	
					}												
				}											
			}
			
			if ( RANK == 0 )
			{
				cout << endl;
				cout << "	Read the following number of ij ik jk bins for pairs: " << endl;
			
				for(int i=0; i<NPAIR; i++)
					cout << "		" << ATOM_PAIRS[i].PRPR_NM  << ": " << ATOM_PAIRS[i].NBINS[0] << " " << ATOM_PAIRS[i].NBINS[1] << " " << ATOM_PAIRS[i].NBINS[2] << endl;					
				cout << endl;
			}

			build_int_pair_map(CONTROLS.NATMTYP, ATOM_TYPE, ATOM_TYPE_IDX, PAIR_MAP, INT_PAIR_MAP);

			for ( int i = 0; i < ATOM_PAIRS.size(); i++ ) 
				if ( ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) 
					ATOM_PAIRS[i].set_cheby_vals();

		}
		
		else if(LINE.find("PAIR CHEBYSHEV CUBIC SCALING")!= string::npos)
		{
			STREAM_PARSER.str(LINE);
			STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR;
			for(int i=0; i<NPAIR; i++)
				ATOM_PAIRS[i].CUBIC_SCALE = double(atof(TEMP_STR.data()));
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();	
		}	
		
		
		else if(LINE.find("CHARGE CONSTRAINTS") != string::npos)
		{		
			if ( RANK == 0 ) 
				cout << endl << "	Attempting to read " << NPAIR-1 << " charge constraints...:" << endl; 
			
			CHARGE_CONSTRAINTS.resize(NPAIR-1);
			for(int i=0; i<NPAIR-1; i++)
			{
				// Read the atom pair types
				
				for(int j=0; j<NPAIR; j++)
				{
					cin >> LINE; 
					CHARGE_CONSTRAINTS[i].PAIRTYPE.push_back(LINE);
					
					CHARGE_CONSTRAINTS[i].PAIRTYPE_IDX.push_back(PAIR_MAP[LINE]);
				}
				
				// Read the constraints 
				
				for(int j=0; j<NPAIR; j++)
				{
					cin >> LINE; 
					CHARGE_CONSTRAINTS[i].CONSTRAINTS.push_back(double(atof(LINE.data())));
				}
				
				// Read the associated force
				cin >> CHARGE_CONSTRAINTS[i].FORCE;				
				cin.ignore();
		
				if ( RANK == 0 ) 
				{
					cout << "		" << i+1 << "	 ";
					for(int j=0; j<NPAIR; j++)
						cout << CHARGE_CONSTRAINTS[i].PAIRTYPE[j] << " (" << CHARGE_CONSTRAINTS[i].PAIRTYPE_IDX[j] << ") ";
					for(int j=0; j<NPAIR; j++)
						cout << CHARGE_CONSTRAINTS[i].CONSTRAINTS[j] << " ";
					cout << CHARGE_CONSTRAINTS[i].FORCE << endl;	
				}	
					
			}
			
			if ( RANK == 0 ) 
				cout << endl;
		}
		
		else if(LINE.find("SPECIAL 3B S_MAXIM:") != string::npos)
		{
		  TRIPS.read_cutoff_params(cin, LINE, "S_MAXIM");
		}
		
		else if(LINE.find("SPECIAL 3B S_MINIM:") != string::npos)
		{
		  TRIPS.read_cutoff_params(cin, LINE, "S_MINIM");
		}
		
		else if(LINE.find("SPECIAL 4B S_MAXIM:") != string::npos)
		{
		  QUADS.read_cutoff_params(cin, LINE, "S_MAXIM");
		}

		else if(LINE.find("SPECIAL 4B S_MINIM:") != string::npos)
		{
		  QUADS.read_cutoff_params(cin, LINE, "S_MINIM");
		}
		
		else if ( (LINE.find("# FCUTTYP #") != string::npos) )
		{
			getline(cin,CONTROLS.FCUT_LINE) ;
			if ( RANK == 0 ) 
			{
				cout << "# FCUTTYP #: " << CONTROLS.FCUT_LINE << "      ... for all Chebyshev interactions" << endl ;
			}		
			/*
			// Unified MD/LSQ parsing of force cutoff.
			getline(cin,TEMP_TYPE) ;
			if ( RANK == 0 ) 
			{
				parse_space(TEMP_TYPE, tokens) ;
				cout << "# FCUTTYP #: " << tokens[0] << "      ... for all Chebyshev interactions" << endl ;
			}
			parse_fcut_input(TEMP_TYPE, ATOM_PAIRS, TRIPS, QUADS) ;
			//cin.ignore();
			*/	
		}
		else if ( RANK == 0 ) 
		{
			cout << "WARNING:  The following input line was ignored:\n" ;
			cout << " " << LINE << endl ;
		}
		
	}	
			
	if ( RANK == 0 ) 
		cout << endl << "Note: Will use cubic scaling of: " << ATOM_PAIRS[0].CUBIC_SCALE << endl << endl; // All types use same scaling
}

static void print_bond_stats(vector<PAIRS> &ATOM_PAIRS, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, bool use_3b_cheby, bool use_4b_cheby)
// Print statistics on bonding.
{

	if ( RANK == 0 ) 
		cout << "	Minimum distances between atoms: (Angstr.)" << endl;

	int NPAIR = ATOM_PAIRS.size();

		for (int k=0; k<NPAIR; k++)
		{
			double sum = 0.0;
#ifdef USE_MPI
			MPI_Reduce(&(ATOM_PAIRS[k].MIN_FOUND_DIST), &sum, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
#else
			sum = ATOM_PAIRS[k].MIN_FOUND_DIST;
#endif
			if ( RANK == 0 ) 
			{
				cout << "		" << k << "	" << ATOM_PAIRS[k].ATM1TYP << 
					" " << ATOM_PAIRS[k].ATM2TYP << "	" << fixed << setprecision(3) << sum << endl;

			}

		}

		if ( RANK == 0 ) 
			cout << "	Total number of configurations contributing to each pair type:" << endl;
		
		for (int k = 0; k < NPAIR; k++) 
		{
			int sum = 0.0;
#ifdef USE_MPI
			MPI_Reduce(&(ATOM_PAIRS[k].N_CFG_CONTRIB), &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
			sum = ATOM_PAIRS[k].N_CFG_CONTRIB;
#endif
			if ( RANK == 0 ) 
				cout << "		" << k << "	" 					
				     << ATOM_PAIRS[k].ATM1TYP << " " 
				     << ATOM_PAIRS[k].ATM2TYP << " " 
				     << sum << endl;
		}

		if( use_3b_cheby ) 
		{
		  TRIPS.print_min_distances();
		}

		if( use_4b_cheby ) 
		{
		  QUADS.print_min_distances();
		}
			
		if ( RANK == 0 ) 
			cout << "...matrix printing complete: " << endl << endl;
}

static void build_clusters(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST &TRIPS, CLUSTER_LIST & QUADS, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST, vector<int>& ATOM_TYPE_IDX, vector<string>& ATOM_TYPE)
// Build the many-body interaction clusters after the input has been read.
{
	if(CONTROLS.USE_3B_CHEBY)
	{
		// Generate unique triplets
		TRIPS.build_all(CONTROLS.CHEBY_3B_ORDER, ATOM_PAIRS, PAIR_MAP,ATOM_TYPE,ATOM_TYPE_IDX); 
		TRIPS.print(false);
		NEIGHBOR_LIST.MAX_CUTOFF_3B = TRIPS.MAX_CUTOFF;
	}	
			
	if(CONTROLS.USE_4B_CHEBY)
	{
		// Generate unique quadruplets and thier corresponding sets of powers
		QUADS.build_all(CONTROLS.CHEBY_4B_ORDER, ATOM_PAIRS, PAIR_MAP, ATOM_TYPE, ATOM_TYPE_IDX);
		QUADS.print(false);
		NEIGHBOR_LIST.MAX_CUTOFF_4B = QUADS.MAX_CUTOFF;
	}

	// Set up the Cheby variables
	if(CONTROLS.USE_3B_CHEBY)
		TRIPS.build_cheby_vals(ATOM_PAIRS);

	if(CONTROLS.USE_4B_CHEBY)
		QUADS.build_cheby_vals(ATOM_PAIRS);
	
	parse_fcut_input(CONTROLS.FCUT_LINE, ATOM_PAIRS, TRIPS, QUADS) ;					
}

static void print_param_header(JOB_CONTROL &CONTROLS, vector<PAIRS> &ATOM_PAIRS,
															 CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS)
// Print the params.header file.
{
	ofstream header;
	header.open("params.header");

	bool USE_POVER ;
	
	//////////////////////////////////////////////////	   
	// THE NEW WAY
	//////////////////////////////////////////////////	

	if(!CONTROLS.FIT_COUL && !CONTROLS.USE_PARTIAL_CHARGES)
		header << "USECOUL: false" << endl;
	else
		header << "USECOUL: true" << endl;

	if(CONTROLS.FIT_COUL)	
		header << "FITCOUL: true" << endl;
	else
		header << "FITCOUL: false" << endl;
	
	USE_POVER = false;
	for(int i=0; i<ATOM_PAIRS.size(); i++)
		if(ATOM_PAIRS[i].USE_OVRPRMS)
			USE_POVER = true;
	
	if(USE_POVER)
		header << "USEPOVR: true" << endl;
	else
		header << "USEPOVR: false" << endl;
	
	if(CONTROLS.FIT_POVER)
		header << "FITPOVR: true" << endl;
	else
		header << "FITPOVR: false" << endl;
	
	if(CONTROLS.USE_3B_CHEBY)
		header << "USE3BCH: true" << endl;
	else
		header << "USE3BCH: false" << endl;
	
	if(CONTROLS.USE_4B_CHEBY)
		header << "USE4BCH: true" << endl;
	else
		header << "USE4BCH: false" << endl;	
	
	header << endl << "PAIRTYP: " << ATOM_PAIRS[0].PAIRTYP << " ";
	
	if     (ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV")
		header << " " << ATOM_PAIRS[0].SNUM << " " << ATOM_PAIRS[0].SNUM_3B_CHEBY << " " << ATOM_PAIRS[0].SNUM_4B_CHEBY << " " << ATOM_PAIRS[0].CHEBY_RANGE_LOW << " " << ATOM_PAIRS[0].CHEBY_RANGE_HIGH << endl;
	else if(ATOM_PAIRS[0].PAIRTYP == "DFTBPOLY")
		header << " " << ATOM_PAIRS[0].SNUM << endl;	
	else if (ATOM_PAIRS[0].PAIRTYP == "INVRSE_R")
		header << " " << CONTROLS.INVR_PARAMS << endl;
	else
		header << endl;
	
	header << endl << "ATOM TYPES: " << CONTROLS.NATMTYP << endl << endl;
	header << "# TYPEIDX #	# ATM_TYP #	# ATMCHRG #	# ATMMASS #" << endl;
	
	
	if(CONTROLS.FIT_COUL)
	{
		for(int i=0; i<CONTROLS.NATMTYP; i++)
			header << i << "		" << ATOM_PAIRS[i].ATM1TYP << "		" << ATOM_PAIRS[i].CHRGSGN << "		" << ATOM_PAIRS[i].ATM1MAS << endl;
	}
	else
	{
		for(int i=0; i<CONTROLS.NATMTYP; i++)
			header << i << "		" << ATOM_PAIRS[i].ATM1TYP << "		" << ATOM_PAIRS[i].ATM1CHG << "		" << ATOM_PAIRS[i].ATM1MAS << endl;		
	}
	
	int NPAIR =  ATOM_PAIRS.size();	
	
	header << endl << "ATOM PAIRS: " << NPAIR << endl << endl;
	
	header << "# PAIRIDX #	";
	header << "# ATM_TY1 #	";
	header << "# ATM_TY1 #	";
	header << "# S_MINIM #	";
	header << "# S_MAXIM #	";
	header << "# S_DELTA #	";
	header << "# CHBDIST #	";	// how pair distance is transformed in cheby calc
	header << "# MORSE_LAMBDA #" << endl;

	for(int i=0; i<NPAIR; i++)
	{

	  string chtype = Cheby::get_trans_string(ATOM_PAIRS[i].CHEBY_TYPE);
		header << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
					 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
					 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
					 << setw(16) << left << ATOM_PAIRS[i].S_MINIM
					 << setw(16) << left << ATOM_PAIRS[i].S_MAXIM							 
					 << setw(16) << left << ATOM_PAIRS[i].S_DELTA
					 << setw(16) << left << chtype;
		if(ATOM_PAIRS[i].CHEBY_TYPE == Cheby_trans::MORSE )
			header << setw(16) << left << ATOM_PAIRS[i].LAMBDA << endl; 
		else
			header << endl;
	}
	

	if(USE_POVER)
	{
		header << endl;
 		header << "# PAIRIDX #	";
 		header << "# ATM_TY1 #	";
 		header << "# ATM_TY1 #	";				
		header << "# USEOVRP #	";
		header << "# TO_ATOM #	";
		header << "# P_OVERB #	";
 		header << "# R_0_VAL #	";
 		header << "# P_1_VAL #	";
 		header << "# P_2_VAL #	";
 		header << "# LAMBDA6 #" << endl;	
			
		for(int i=0; i<NPAIR; i++)
		{
			header << "	"
						 << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
						 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
						 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
						 << setw(16) << left << ATOM_PAIRS[i].USE_OVRPRMS;
			if(ATOM_PAIRS[i].USE_OVRPRMS)
				header	<< setw(16) << left << ATOM_PAIRS[i].OVER_TO_ATM;
			else
				header	<< setw(16) << left << "NONE";
			
			header
			  << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[0] 
			  << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[1]
			  << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[2] 
			  << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[3]
			  << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[4] << endl;	
		}		
	}

	// Quads and triplets both have the same cutoff function parameters.
	// Print out only once.
	if ( ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV" ) 
		ATOM_PAIRS[0].FORCE_CUTOFF.print_header(header) ;
		
	if(ATOM_PAIRS[0].CUBIC_SCALE != 1.0)
		header << endl << "PAIR CHEBYSHEV CUBIC SCALING: " << ATOM_PAIRS[0].CUBIC_SCALE << endl;
	
	// Print out special cutoffs 
	TRIPS.print_special(header);
	QUADS.print_special(header);
	 

	// Print out cluster parameters into the header.
	TRIPS.print_header(header,3,CONTROLS.CHEBY_3B_ORDER);
	QUADS.print_header(header,4,CONTROLS.CHEBY_4B_ORDER);

	header << endl;
	header.close();
}

static void print_map_file(map<string,int> PAIR_MAP, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS)
{
			
	ofstream MAPFILE;
	MAPFILE.open("ff_groups.map");
	
	MAPFILE << "PAIRMAPS: " << PAIR_MAP.size() << endl;
	
	for(map<string,int>::iterator i=PAIR_MAP.begin(); i!=PAIR_MAP.end(); i++)
	{
		MAPFILE << i->second << " " << i->first << endl;
	}

	if(TRIPS.MAP.size() > 0)
	{
		MAPFILE << endl;
		
		MAPFILE << "TRIPMAPS: " << TRIPS.MAP.size() << endl;
		
		for(map<string,int>::iterator i=TRIPS.MAP.begin(); i!=TRIPS.MAP.end(); i++)
			MAPFILE << i->second << " " << i->first << endl;
	}
	
	if(QUADS.MAP.size() > 0)
	{
		MAPFILE << endl;
		
		MAPFILE << "QUADMAPS: " << QUADS.MAP.size() << endl;
		
		for(map<string,int>::iterator i=QUADS.MAP.begin(); i!=QUADS.MAP.end(); i++)
			MAPFILE << i->second << " " << i->first << endl;
	}
	
	MAPFILE.close();

}

static int process_frame(A_MAT &A_MATRIX, JOB_CONTROL &CONTROLS, FRAME &SYSTEM, vector<PAIRS> &ATOM_PAIRS,
												 map<string,int> &PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS &NEIGHBOR_LIST,
												 CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, vector<CHARGE_CONSTRAINT> &CHARGE_CONSTRAINTS,
												 int i, int istart) 
// Process the ith frame of the A_MATRIX.
{

	double 	NEIGHBOR_PADDING = 0.3;

	// NOTE: WE CONTINUALLY RE-USE THE 0th entry of A_MATRIX TO SAVE MEMORY.
	A_MATRIX.INITIALIZE(CONTROLS, SYSTEM, ATOM_PAIRS.size() ) ;

	if (RANK == 0 && i == istart )
	{
		cout << "...matrix setup complete: " << endl << endl;
		cout << "...Populating the matrices for A, Coulomb forces, and overbonding..." << endl << endl;
	}

	// Only include stress tensor data for first NSTRESS frames..
	if(i >= CONTROLS.NSTRESS)
	{
		CONTROLS.FIT_STRESS     = false;	
		CONTROLS.FIT_STRESS_ALL = false;
	}
		
	if(i >= CONTROLS.NENER)
	{
		CONTROLS.FIT_ENER          = false;	
		//CONTROLS.FIT_ENER_PER_ATOM = false;
	}		
	
	// This output is specific to the number of processors.
		
	if(NPROCS==1)
		cout << "	Processing frame: " << setw(5) << i+1 << " of: " << CONTROLS.NFRAMES << endl;

	// Use very little padding because we will update neighbor list for every frame.
		
	NEIGHBOR_LIST.INITIALIZE(SYSTEM, NEIGHBOR_PADDING);
	NEIGHBOR_LIST.DO_UPDATE (SYSTEM, CONTROLS);		

	ZCalc_Deriv(CONTROLS, ATOM_PAIRS, TRIPS, QUADS, SYSTEM, A_MATRIX, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST);
		
	if ( CONTROLS.IF_SUBTRACT_COORD ) // Subtract over-coordination forces from force to be output.
	{
		//SubtractCoordForces(SYSTEM, false, i, A_MATRIX,  ATOM_PAIRS, PAIR_MAP, NEIGHBOR_LIST, true);	
		cout << "Feature deprecated - exiting." << endl;
		exit_run(0);
	}
		
	if (CONTROLS.IF_SUBTRACT_COUL) 
		SubtractEwaldForces(SYSTEM, NEIGHBOR_LIST, CONTROLS);

	if ( CONTROLS.FIT_POVER )	// Fit the overcoordination parameter.
	{
		SubtractCoordForces(SYSTEM, true, A_MATRIX, ATOM_PAIRS, PAIR_MAP, NEIGHBOR_LIST, true);			
		//cout << "Feature deprecated - exiting." << endl;
		//exit_run(0);
	}

	A_MATRIX.PRINT_FRAME(CONTROLS, SYSTEM, ATOM_PAIRS, CHARGE_CONSTRAINTS, i) ;
	int total_forces = 3 * A_MATRIX.FORCES.size() ;
	return total_forces ;
}
