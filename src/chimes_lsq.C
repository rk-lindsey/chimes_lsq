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
#include "input.h"

using namespace std;

#ifndef VERBOSITY 
	#define VERBOSITY 1 
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif 

static void read_lsq_input(	string & INFILE,
				JOB_CONTROL & CONTROLS, 
				vector<PAIRS> & ATOM_PAIRS, 
				CLUSTER_LIST &TRIPS, 
				CLUSTER_LIST & QUADS, 
				map<string,int> & PAIR_MAP,
				vector<int> &INT_PAIR_MAP,
				vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, 
				NEIGHBORS & NEIGHBOR_LIST,
				vector<int>& TMP_ATOMTYPEIDX, 
				vector<string> &TMP_ATOMTYPE);

static void print_bond_stats(	vector<PAIRS> &ATOM_PAIRS, 
				CLUSTER_LIST &TRIPS, 
				CLUSTER_LIST &QUADS, 
				bool use_3b_cheby, 
				bool use_4b_cheby);

static void build_clusters(	JOB_CONTROL & CONTROLS, 
				vector<PAIRS> & ATOM_PAIRS, 
				CLUSTER_LIST &TRIPS, 
				CLUSTER_LIST & QUADS, 
				map<string,int> & PAIR_MAP, 
				NEIGHBORS & NEIGHBOR_LIST, 
				vector<int>& TMP_ATOMTYPEIDX, 
				vector<string>& TMP_ATOMTYPE);
				
static void print_param_header(	JOB_CONTROL &CONTROLS, 
				vector<PAIRS> &ATOM_PAIRS,
				CLUSTER_LIST &TRIPS, 
				CLUSTER_LIST &QUADS) ;
				
static void print_map_file(	map<string,int> PAIR_MAP, 
				CLUSTER_LIST &TRIPS, 
				CLUSTER_LIST &QUADS) ;
				
static int process_frame(	A_MAT &A_MATRIX, 
				JOB_CONTROL &CONTROLS, 
				FRAME &SYSTEM, 
				vector<PAIRS> &ATOM_PAIRS,
				map<string,int> &PAIR_MAP, 
				vector<int> &INT_PAIR_MAP, 
				NEIGHBORS &NEIGHBOR_LIST,
				CLUSTER_LIST &TRIPS, 
				CLUSTER_LIST &QUADS,
				 vector<CHARGE_CONSTRAINT> &CHARGE_CONSTRAINTS,
				int i, 
				int istart)  ;

// Print the params.header file.
// Global variables declared as externs in functions.h, and declared in functions.C

string FULL_FILE_3B;		
string SCAN_FILE_3B;
string SCAN_FILE_2B;
ofstream BAD_CONFIGS_1;	// Configs where r_ij < r_cut,in 
ofstream BAD_CONFIGS_2;	// Configs where r_ij < r_cut,in +d_penalty
ofstream BAD_CONFIGS_3; // All other configs, but only printed when (CONTROLS.FREQ_DFTB_GEN>0) && ((CONTROLS.STEP+1) % CONTROLS.FREQ_DFTB_GEN == 0)

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
	
	if ( argc != 2 )
		EXIT_MSG("Usage: chimes_lsq input_file > <output_file") ;
			
	string INFILE   = argv[1];
	
	vector<PAIRS> 	ATOM_PAIRS;			// Will store relevant info regarding atom interaction pair types.. 
	CLUSTER_LIST 	TRIPS;				// Triplet interaction: generic cluster interface.
	CLUSTER_LIST 	QUADS;				// Quadruplet interaction: generic cluster interface.
	
	FRAME		SYSTEM ;			// Stores the trajectory information... box dimensions, atom types, coordinates, and forces
	NEIGHBORS       NEIGHBOR_LIST;			// Declare the class that will handle the neighbor list
	map<string,int> PAIR_MAP;			// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
	
	vector<CHARGE_CONSTRAINT> CHARGE_CONSTRAINTS;	// Specifies how we constrain charge fitting

	JOB_CONTROL 	CONTROLS;			// Will hold job controls shared by both lsq and md
	CONTROLS.IS_LSQ  = true; 

	
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

	read_lsq_input(INFILE, CONTROLS, ATOM_PAIRS, TRIPS, QUADS, PAIR_MAP, INT_PAIR_MAP, CHARGE_CONSTRAINTS, NEIGHBOR_LIST, ATOM_TYPE_IDX, ATOM_TYPE);

	// Build many-body interaction clusters if necessary.
	build_clusters(CONTROLS, ATOM_PAIRS, TRIPS, QUADS, PAIR_MAP, NEIGHBOR_LIST, ATOM_TYPE_IDX, ATOM_TYPE);

	if ( RANK == 0 ) 
		cout << "...input file read successful: " << endl << endl;
	
	// if ( NPROCS > CONTROLS.NFRAMES ) NPROCS = CONTROLS.NFRAMES ;	// Don't use unnecessary processors.

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

	CONTROLS.LSQ_SETUP( ATOM_PAIRS.size(), ATOM_TYPE.size() ) ;

	//////////////////////////////////////////////////
	//
	// Setup variables for reading the trajectory file(s)
	//
	//////////////////////////////////////////////////
	
	cout.precision(16);				// WE SHOULD AT LEAST MOVE THIS TO SOMEWHERE MORE REASONABLE.. LIKE THE SECTION WHERE THE OUTPUT IS ACTUALLY PRINTED
	
	int istart, iend;
	
	divide_atoms(istart, iend, CONTROLS.NFRAMES);	// Each processor only calculates certain frames.  Calculate which frames to use.

	//////////////////////////////////////////////////
	//	
	// Begin processing the trajecory
	//
	//////////////////////////////////////////////////
	 
	int FILE_IDX = 0;		// Index of traj file in CONTROLS.INFILE vector
	ifstream TRAJ_INPUT;

	OPEN_TRAJFILE(TRAJ_INPUT, CONTROLS.INFILE, FILE_IDX);
	
	if ( RANK == 0 ) 
		cout << endl << "Succefully opened the trajectory file..." << endl;

	if ( RANK == 0 ) 
		cout << "Setting up the matrices for A, Coulomb forces, and overbonding..." << endl;
		
	int OFFSET = CONTROLS.INFILE_FRAMES[FILE_IDX];

	A_MATRIX.OPEN_FILES(CONTROLS) ;
	int total_forces = 0 ;
	
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
	A_MATRIX.CLEANUP_FILES(CONTROLS.SPLIT_FILES) ;
	
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
static void read_lsq_input(string & INFILE,
			   JOB_CONTROL & CONTROLS, 
                           vector<PAIRS> & ATOM_PAIRS, 
                           CLUSTER_LIST &TRIPS, 
			   CLUSTER_LIST & QUADS, 
			   map<string,int> & PAIR_MAP,
			   vector<int> &INT_PAIR_MAP, 
                           vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, 
                           NEIGHBORS & NEIGHBOR_LIST,
                           vector<int>& ATOM_TYPE_IDX,
                           vector<string>& ATOM_TYPE)
{

	// Instantiate the input class

	INPUT LSQ_INPUT(INFILE);
	
	if (RANK==0)
		cout << "Will read from file: " << INFILE << endl;
	
	// Read the LSQ input file
	
	LSQ_INPUT.PARSE_INFILE_LSQ(	CONTROLS, 
					ATOM_PAIRS, 
					TRIPS, 
					QUADS, 
					PAIR_MAP, 
					INT_PAIR_MAP,
					CHARGE_CONSTRAINTS, 
					NEIGHBOR_LIST,
					ATOM_TYPE_IDX,
					ATOM_TYPE);

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

static int process_frame(	A_MAT &A_MATRIX, 
				JOB_CONTROL &CONTROLS, 
				FRAME &SYSTEM, 
				vector<PAIRS> &ATOM_PAIRS,
				map<string,int> &PAIR_MAP, 
				vector<int> &INT_PAIR_MAP, 
				NEIGHBORS &NEIGHBOR_LIST,
				CLUSTER_LIST &TRIPS, 
				CLUSTER_LIST &QUADS, 
				vector<CHARGE_CONSTRAINT> &CHARGE_CONSTRAINTS,
				int i, 
				int istart) 
// Process the ith frame of the A_MATRIX.
{

	double 	NEIGHBOR_PADDING = 0.3;

	// NOTE: WE CONTINUALLY RE-USE THE 0th entry of A_MATRIX TO SAVE MEMORY.
	A_MATRIX.INITIALIZE(CONTROLS, SYSTEM, ATOM_PAIRS.size(),ATOM_PAIRS) ;

	if (RANK == 0 && i == istart )
	{
		cout << "...matrix setup complete: " << endl << endl;
		cout << "...Populating the matrices for A, Coulomb forces, and overbonding..." << endl << endl;
	}
	
		
	bool DUMMY_FIT_STRESS	     = CONTROLS.FIT_STRESS;    
	bool DUMMY_FIT_STRESS_ALL    = CONTROLS.FIT_STRESS_ALL;

	bool DUMMY_FIT_ENER	     = CONTROLS.FIT_ENER;	 
	//bool DUMMY_FIT_ENER_PER_ATOM = CONTROLS.FIT_ENER_PER_ATOM;

	// Only include stress tensor data for first NSTRESS frames..
	
	if((CONTROLS.NSTRESS != -1) && (i >= CONTROLS.NSTRESS))
	{
		CONTROLS.FIT_STRESS     = false;	
		CONTROLS.FIT_STRESS_ALL = false;
	}
		
	if((CONTROLS.NENER != -1) && (i >= CONTROLS.NENER))
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
		SubtractCoordForces(SYSTEM, true, A_MATRIX, ATOM_PAIRS, PAIR_MAP, NEIGHBOR_LIST, true);			

	
	CONTROLS.FIT_STRESS        = DUMMY_FIT_STRESS; 
	CONTROLS.FIT_STRESS_ALL    = DUMMY_FIT_STRESS_ALL;

	CONTROLS.FIT_ENER          = DUMMY_FIT_ENER;	
	//CONTROLS.FIT_ENER_PER_ATOM = DUMMY_FIT_ENER_PER_ATOM;	


	CONTROLS.FIT_ENER_EVER = (CONTROLS.FIT_ENER || CONTROLS.FIT_ENER_PER_ATOM);

	A_MATRIX.PRINT_FRAME(CONTROLS, SYSTEM, ATOM_PAIRS, CHARGE_CONSTRAINTS, i) ;
	int total_forces = 3 * A_MATRIX.FORCES.size() ;
	return total_forces ;
}
