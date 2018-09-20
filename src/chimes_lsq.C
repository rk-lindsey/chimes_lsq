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

using namespace std;

#ifndef VERBOSITY 
	#define VERBOSITY 1 
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif 

static void read_lsq_input(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, CLUSTER_LIST &TRIPS, CLUSTER_LIST & QUADS, map<string,int> & PAIR_MAP, map<int,string> & PAIR_MAP_REVERSE, vector<int> &INT_PAIR_MAP, vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, NEIGHBORS & NEIGHBOR_LIST, vector<int>& TMP_ATOMTYPEIDX, vector<string> &TMP_ATOMTYPE);

static void print_bond_stats(vector<PAIRS> &ATOM_PAIRS, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, bool use_3b_cheby, bool use_4b_cheby) ;

static void add_col_of_ones(int item, bool DO_ENER, bool DO_STRESS, bool DO_STRESS_ALL, int NATOMS, ofstream & OUTFILE);

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

ofstream BAD_CONFIGS ;

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
	
	
	//////////////////////////////////////////////////
	//
	// Job parameters: Declare and set defaults
	//
	//////////////////////////////////////////////////
	
	// MAJOR ASSUMPTIONS: 
	// 
	// 1. Number of atoms does not change over frames ( should always be true for MD)
	// 2. Atom ordering in trajectory file does not change over frames ( should also always be true)
	
	vector<PAIRS> 		ATOM_PAIRS;			// Will store relevant info regarding atom interaction pair types.. 
	CLUSTER_LIST TRIPS ;                // Triplet interaction: generic cluster interface.
	CLUSTER_LIST QUADS ;                // Quadruplet interaction: generic cluster interface.
	
	vector<FRAME> 		TRAJECTORY;			// Stores the trajectory information... box dimensions, atom types, coordinates, and forces
	NEIGHBORS        	 NEIGHBOR_LIST;		// Declare the class that will handle the neighbor list
	vector<int>   		ATOM_PAIR_TYPES;	// Fore use in double loop over atom pairs. Index corresponds to the overall double loop count. 
											// Provides an index that rells you the atom pair's ATOM_PAIRS's type index.. THIS IS FOR LOOPS OF TYPE
											// 	for(int a1=0;a1<nat-1;a1++)	for(int a2=a1+1;a2<nat;a2++)
	vector<int>   ATOM_PAIR_TYPES_ALL;		// THIS IS FOR LOOPS OF TYPE for(int ai=0; ai<nat; ai++), for(int ak=0; ak<nat; ak++)
	
	map<string,int> PAIR_MAP;			// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
	map<int,string> PAIR_MAP_REVERSE; 	// Input/output the resverse of PAIR_MAP
	
	vector<CHARGE_CONSTRAINT> CHARGE_CONSTRAINTS;	// Specifies how we constrain charge fitting

	JOB_CONTROL 	CONTROLS;			// Will hold job controls shared by both lsq and md
	
	vector<int>	INT_PAIR_MAP;

	vector<int>   TMP_ATOMTYPEIDX;			// Used to construct the quadruplet type index
	vector<string> TMP_ATOMTYPE;			// Used to construct the quadruplet type index

	//////////////////////////////////////////////////
	//
	// Read and print input to screen
	//
	//////////////////////////////////////////////////
	
	#if VERBOSITY == 1
		if ( RANK == 0 ) 
		{
			cout << endl << "Reading input file..." << endl;
	
			// Delete temporary files if they exist.
			system("rm -f A.[0-9]*.txt") ;
			system("rm -f b.[0-9]*.txt") ;
		}
	#endif

	read_lsq_input(CONTROLS, ATOM_PAIRS, TRIPS, QUADS, PAIR_MAP, PAIR_MAP_REVERSE, INT_PAIR_MAP,  CHARGE_CONSTRAINTS, NEIGHBOR_LIST, TMP_ATOMTYPEIDX, TMP_ATOMTYPE);
 

	if((CONTROLS.FIT_STRESS || CONTROLS.FIT_STRESS_ALL) && CONTROLS.CALL_EWALD)
	{
		cout << "ERROR: Inclusion of stress tensors currently not compatible with use of ZCalc_Ewald_Deriv." << endl;
		exit(0);
	}

	#if VERBOSITY == 1
		if ( RANK == 0 ) cout << "...input file read successful: " << endl << endl;
	#endif


	//////////////////////////////////////////////////
	//
	// Set up force and force derivative vectors
	//
	//////////////////////////////////////////////////	
		 		
		
	// Quick note on the structure of A_MATRIX's last column:
	// Suppose we have 3 pair types, O--O, O--H, and H--H
	// Suppose also that the pair type potentials have 4 parameters each
	// When laid out as a vector, any given final column in A_MATRIX will be arranged as:
	// 
	// { (Param-1_O--O), (Param-2_O--O), (Param-3_O--O), (Param-4_O--O), 
	//   (Param-1_O--H), (Param-2_O--H), (Param-3_O--H), (Param-4_O--H), 
	//   (Param-1_H--H), (Param-2_H--H), (Param-3_H--H), (Param-4_H--H)	}	
		
	vector<vector <vector <XYZ > > > A_MATRIX;	    	// [ # frames ] [ # atoms ]      [ #fitting parameters ]  .. replaces Als
	vector<vector <vector <XYZ > > > COULOMB_FORCES;	// [ # frames ] [ # pair types ] [ # atoms ]              .. replaces Coul_xx
	vector<vector <XYZ > >           P_OVER_FORCES;  	// [ # frames ] [ # atoms ]                               .. replaces Pover
	
	A_MATRIX      .resize(CONTROLS.NFRAMES);
	COULOMB_FORCES.resize(CONTROLS.NFRAMES);
	P_OVER_FORCES .resize(CONTROLS.NFRAMES);

	// Figure out necessary dimensions for the force/force derivative vectors
	
	CONTROLS.TOT_SNUM = 0; 
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
			}
				
        	}
		
		else if (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" )	// Set the distance transformation type
			ATOM_PAIRS[i].CHEBY_TYPE          = CONTROLS.CHEBY_TYPE;
			
		else if (ATOM_PAIRS[i].PAIRTYP == "INVRSE_R") 
			ATOM_PAIRS[i].SNUM = CONTROLS.INVR_PARAMS;

		else // Spline
			ATOM_PAIRS[i].SNUM = (2+floor((ATOM_PAIRS[i].S_MAXIM - ATOM_PAIRS[i].S_MINIM)/ATOM_PAIRS[i].S_DELTA))*2; //2 is for p0/m0/p1/m1.. 
	
		CONTROLS.TOT_SNUM += ATOM_PAIRS[i].SNUM;
	}
	
	#if VERBOSITY == 1		
		if ( RANK == 0 ) cout << "The number of two-body non-coulomb parameters is: " << CONTROLS.TOT_SNUM <<  endl;
	#endif

	
	if ( ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV" && CONTROLS.CHEBY_3B_ORDER  > 0 ) // All atoms types must use same potential, so check if 3b order is greater than zero 
	{
		CONTROLS.NUM_3B_CHEBY = 0;
		
		for(int i=0; i< TRIPS.VEC.size(); i++)
			CONTROLS.NUM_3B_CHEBY += TRIPS.VEC[i].N_TRUE_ALLOWED_POWERS;

	#if VERBOSITY == 1
		if ( RANK == 0 ) cout << "The number of three-body Chebyshev parameters is: " << CONTROLS.NUM_3B_CHEBY << endl;
	#endif
	}
	
	if ( ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV" && CONTROLS.CHEBY_4B_ORDER  > 0 ) // All atoms types must use same potential, so check if 3b order is greater than zero 
	{
		CONTROLS.NUM_4B_CHEBY = 0;
		
		for(int i=0; i< QUADS.VEC.size(); i++)
			CONTROLS.NUM_4B_CHEBY += QUADS.VEC[i].N_TRUE_ALLOWED_POWERS;

	#if VERBOSITY == 1
		if ( RANK == 0 ) cout << "The number of four-body  Chebyshev parameters is: " << CONTROLS.NUM_4B_CHEBY << endl;
	#endif
	}	

	CONTROLS.TOT_SHORT_RANGE = CONTROLS.TOT_SNUM + CONTROLS.NUM_3B_CHEBY + CONTROLS.NUM_4B_CHEBY;	

	//////////////////////////////////////////////////
	//
	// Read and store MD trajectory (coords and forces)
	//
	//////////////////////////////////////////////////
	
	// Start off with a warning for the MD code... This 
	// won't be necessary once Larry's ghost atom fix 
	// is implemented in the MD code
	
	// Read in the trajectory
					 
	ifstream TRAJ_INPUT;
	TRAJ_INPUT.open(CONTROLS.INFILE.data());
	
	if(!TRAJ_INPUT.is_open())
	{
		cout << "ERROR: Cannot open trajectory file: " << CONTROLS.INFILE << endl;
		exit(1);
	}
	
	TRAJECTORY.resize(CONTROLS.NFRAMES);
	
	#if VERBOSITY == 1
		if ( RANK == 0 ) cout << endl << "Reading in the trajectory file..." << endl;
	#endif
		
	// Prepare for layering
		
	vector<XYZ_INT> LAYER_MATRX;

	for (int i=0; i<CONTROLS.NFRAMES; i++)
	{
		// Read in line with the number of atoms
		
		TRAJ_INPUT >> TRAJECTORY[i].ATOMS;
		
		// Read in line with box dimenstions
		
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.X;
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.Y;
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.Z;

		if((CONTROLS.NSTRESS < 0) || (i<CONTROLS.NSTRESS))
		{
			if(CONTROLS.FIT_STRESS)
			{
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS.X;
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS.Y;
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS.Z;
			}
			else if(CONTROLS.FIT_STRESS_ALL)	// Expects as:xx yy zz xy xz yz // old: xx, xy, xz, yy, yx, yz, zx, zy, zz
			{
				// Read only the "upper" deviatoric components
		
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_X.X;
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_Y.Y;
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_Z.Z;
					
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_X.Y;
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_X.Z;
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_Y.Z;
					
			}
		}
		if((CONTROLS.NENER < 0) || (i<CONTROLS.NENER))
			if(CONTROLS.FIT_ENER) // We're actually fitting to relative *differences* in energy
				TRAJ_INPUT >> TRAJECTORY[i].QM_POT_ENER;

	
		// Check that outer cutoffs do not exceed half of the boxlength
		// with consideration of layering
		
		for(int j=0; j<ATOM_PAIRS.size(); j++)
		{
                        if( (  ATOM_PAIRS[j].S_MAXIM > 0.5* TRAJECTORY[i].BOXDIM.X * (2*CONTROLS.N_LAYERS +1) || ATOM_PAIRS[j].S_MAXIM > 0.5* TRAJECTORY[i].BOXDIM.Y * (2*CONTROLS.N_LAYERS +1) || ATOM_PAIRS[j].S_MAXIM > 0.5* TRAJECTORY[i].BOXDIM.Z * (2*CONTROLS.N_LAYERS +1) ))
			{
					#if WARN == TRUE
						if (isatty(fileno(stdout)) && RANK == 0)
						{
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "WARNING: Outer cutoff greater than half at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Frame:                " << i << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Pair type:            " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Boxlengths:           " << TRAJECTORY[i].BOXDIM.X << " " << TRAJECTORY[i].BOXDIM.Y << " " << TRAJECTORY[i].BOXDIM.Z << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Layers:               " << CONTROLS.N_LAYERS << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Effective boxlengths: " << TRAJECTORY[i].BOXDIM.X * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Y * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Z * (CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
						}
						else if ( RANK == 0 ) 
						{
							cout << "WARNING: Outer cutoff greater than half at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<endl;
							cout << "	Frame:                " << i << COUT_STYLE.ENDSTYLE << endl;
							cout << "	Pair type:            " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
							cout << "	Boxlengths:           " << TRAJECTORY[i].BOXDIM.X << " " << TRAJECTORY[i].BOXDIM.Y << " " << TRAJECTORY[i].BOXDIM.Z << endl;
							cout << "	Layers:               " << CONTROLS.N_LAYERS << endl;
							cout << "	Effective boxlengths: " << TRAJECTORY[i].BOXDIM.X * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Y * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Z * (CONTROLS.N_LAYERS +1) << endl;							
						}

					#else
						if (isatty(fileno(stdout)) && RANK == 0)
						{
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "ERROR: Outer cutoff greater than half at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Frame:                " << i << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Pair type:            " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Boxlengths:           " << TRAJECTORY[i].BOXDIM.X << " " << TRAJECTORY[i].BOXDIM.Y << " " << TRAJECTORY[i].BOXDIM.Z << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Layers:               " << CONTROLS.N_LAYERS << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Effective boxlengths: " << TRAJECTORY[i].BOXDIM.X * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Y * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Z * (CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
							exit(0);
						}
						else if ( RANK == 0 ) 
						{
							cout << "ERROR: Outer cutoff greater than half at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<endl;
							cout << "	Frame:                " << i << COUT_STYLE.ENDSTYLE << endl;
							cout << "	Pair type:            " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << endl;		
							cout << "	Boxlengths:           " << TRAJECTORY[i].BOXDIM.X << " " << TRAJECTORY[i].BOXDIM.Y << " " << TRAJECTORY[i].BOXDIM.Z << endl;
							cout << "	Layers:               " << CONTROLS.N_LAYERS << endl;
							cout << "	Effective boxlengths: " << TRAJECTORY[i].BOXDIM.X * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Y * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Z * (CONTROLS.N_LAYERS +1) << endl;
							exit(0);
						}

					#endif
			}
		}
		
		// Setup the trajectory-holding data object
		
		TRAJECTORY[i].ATOMTYPE    .resize(TRAJECTORY[i].ATOMS);		
		TRAJECTORY[i].COORDS      .resize(TRAJECTORY[i].ATOMS);
		TRAJECTORY[i].FORCES      .resize(TRAJECTORY[i].ATOMS);	// Use for read-in forces.
		TRAJECTORY[i].ACCEL       .resize(TRAJECTORY[i].ATOMS); // Use for calculated forces in ZCalc_Ewald.
		TRAJECTORY[i].CHARGES     .resize(TRAJECTORY[i].ATOMS);
		TRAJECTORY[i].ATOMTYPE_IDX.resize(TRAJECTORY[i].ATOMS);
			
		// Read trajectory, convert to proper units, and apply PBC
			
		for (int j=0; j<TRAJECTORY[i].ATOMS; j++)
		{
			TRAJ_INPUT >> TRAJECTORY[i].ATOMTYPE[j];
			
			for(int k=0; k<TMP_ATOMTYPE.size(); k++)
				if(TRAJECTORY[i].ATOMTYPE[j] == TMP_ATOMTYPE[k])
					TRAJECTORY[i].ATOMTYPE_IDX[j] = k;
			
			TRAJ_INPUT >> TRAJECTORY[i].COORDS[j].X;
			TRAJ_INPUT >> TRAJECTORY[i].COORDS[j].Y;
			TRAJ_INPUT >> TRAJECTORY[i].COORDS[j].Z;
			
			TRAJ_INPUT >> TRAJECTORY[i].FORCES[j].X;
			TRAJ_INPUT >> TRAJECTORY[i].FORCES[j].Y;
			TRAJ_INPUT >> TRAJECTORY[i].FORCES[j].Z;

			if (ATOM_PAIRS[0].PAIRTYP != "DFTBPOLY") // Convert forces to kcal/mol/Angs (Stillinger's units) ... Note, all atom pairs must be of the same type, so using 0 index is ok.
			{
				// Assume units are in Hartree/bohr
				
				TRAJECTORY[i].FORCES[j].X *= 627.50960803*1.889725989;
				TRAJECTORY[i].FORCES[j].Y *= 627.50960803*1.889725989;
				TRAJECTORY[i].FORCES[j].Z *= 627.50960803*1.889725989;
				
				
				// Assume units are in eV/A
				
				/*
				TRAJECTORY[i].FORCES[j].X *= 23.0609;
				TRAJECTORY[i].FORCES[j].Y *= 23.0609;
				TRAJECTORY[i].FORCES[j].Z *= 23.0609;
				*/
			}
						
			if(CONTROLS.WRAP_COORDS)	// Apply PBC (for cases of unwrapped coordinates)
			{
				TRAJECTORY[i].COORDS[j].X -= floor(TRAJECTORY[i].COORDS[j].X/TRAJECTORY[i].BOXDIM.X)*TRAJECTORY[i].BOXDIM.X;
				TRAJECTORY[i].COORDS[j].Y -= floor(TRAJECTORY[i].COORDS[j].Y/TRAJECTORY[i].BOXDIM.Y)*TRAJECTORY[i].BOXDIM.Y;
				TRAJECTORY[i].COORDS[j].Z -= floor(TRAJECTORY[i].COORDS[j].Z/TRAJECTORY[i].BOXDIM.Z)*TRAJECTORY[i].BOXDIM.Z;
			}			
			
			// Assign atom charges.
			if ( CONTROLS.IF_SUBTRACT_COUL ) 
				for(int ii=0; ii< CONTROLS.NATMTYP; ii++)
					if( TRAJECTORY[i].ATOMTYPE[j] == ATOM_PAIRS[ii].ATM1TYP )
					{
						TRAJECTORY[i].CHARGES[j] = ATOM_PAIRS[ii].ATM1CHG;
						break;								
					}
		}
		
		// If layering requested, replicate the system

		build_layers(TRAJECTORY[i], CONTROLS) ;
		
		if(i==0)
		{
			if ( (CONTROLS.N_LAYERS > 0) )	// Then ghost atoms are used 
			{
				if ( RANK == 0 ) 
				{
					cout << "	Reporting outcome of layering for first frame ONLY: " << endl;
					cout << "	Real atoms:                   " << TRAJECTORY[i].ATOMS << endl;
					cout << "	Total atoms (ghost):          " << TRAJECTORY[i].ALL_ATOMS << endl;
					cout << "	Real box dimesntions:         " << TRAJECTORY[i].BOXDIM.X << " " << TRAJECTORY[i].BOXDIM.Y << " " << TRAJECTORY[i].BOXDIM.Z << endl;
					cout << "	Total box dimensions (ghost): " << TRAJECTORY[i].BOXDIM.X * (2*CONTROLS.N_LAYERS + 1) << " " << TRAJECTORY[i].BOXDIM.Y * (2*CONTROLS.N_LAYERS + 1) << " " << TRAJECTORY[i].BOXDIM.Z * (2*CONTROLS.N_LAYERS + 1) << endl << endl;
				}
			
				if(CONTROLS.WRAP_COORDS)
				{
					if ( RANK == 0 ) cout << "WARNING: Coordinate wrapping not supported for ghost atom use. Turning option off" << endl;
					CONTROLS.WRAP_COORDS = false;
				}
			}
			else if ( RANK==0 ) // No ghost atoms.
				cout << "WARNING: Ghost atoms/implicit layers are NOT being used." << endl;
		}
	}
	
	#if VERBOSITY == 1
		if ( RANK == 0 ) cout << "...trajectory file read successful: " << endl << endl;
	#endif

	//////////////////////////////////////////////////
	//
	// Setup A matrix, b vector, for least-squares:
	//
	//////////////////////////////////////////////////

	// Setup the A and Coulomb "matricies"
	
	#if VERBOSITY == 1
	if ( RANK == 0 ) cout << "Setting up the matricies for A, Coulomb forces, and overbonding..." << endl;
	#endif

	int A_SIZE;	// Do we iterate over atoms, or over atoms+1, etc  (where the latter allows us to account for inclusion of stress tensors and/or frame energy)

	for (int f=0; f<CONTROLS.NFRAMES; f++)
	{	
		// Set up the A "matrix"
		
		A_SIZE = TRAJECTORY[f].ATOMS;

		if (CONTROLS.FIT_STRESS)
			A_SIZE++;
		else if (CONTROLS.FIT_STRESS_ALL)
			A_SIZE +=3;

		if (CONTROLS.FIT_ENER)
			A_SIZE++; 

		A_MATRIX[f].resize(A_SIZE);

		for (int i=0; i<A_SIZE; i++)
		{
			A_MATRIX[f][i].resize(CONTROLS.TOT_SHORT_RANGE); 

			for (int j=0; j<CONTROLS.TOT_SHORT_RANGE; j++)
			{		
				A_MATRIX[f][i][j].X = 0.0;
				A_MATRIX[f][i][j].Y = 0.0;
				A_MATRIX[f][i][j].Z = 0.0;
			}
		}

		// Setup the Coulomb force "matrix"
	 
		COULOMB_FORCES[f].resize(ATOM_PAIRS.size());
		P_OVER_FORCES[f] .resize(A_SIZE);
		
		for (int i=0; i<ATOM_PAIRS.size(); i++)
		{
			COULOMB_FORCES[f][i].resize(A_SIZE);		
		
			for (int j=0; j<A_SIZE; j++)
			{			
				COULOMB_FORCES[f][i][j].X = 0;
				COULOMB_FORCES[f][i][j].Y = 0;
				COULOMB_FORCES[f][i][j].Z = 0;	

				if (i==0)	
				{
					P_OVER_FORCES[f][j].X = 0;
					P_OVER_FORCES[f][j].Y = 0;
					P_OVER_FORCES[f][j].Z = 0;						
				}
			}
		}			
	}

	#if VERBOSITY == 1
		if (RANK == 0 ) cout << "...matrix setup complete: " << endl << endl;
	#endif
		
	//////////////////////////////////////////////////
	//
	// Set up the 3B histograms
	//
	//////////////////////////////////////////////////	
	
	if(CONTROLS.USE_3B_CHEBY)
	{
		#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "Setting up the histograms to handle sparse 3b-distance populations..." << endl;
		#endif
	
		for (int i=0; i<TRIPS.VEC.size(); i++)
	     		if ( ! TRIPS.VEC[i].init_histogram(ATOM_PAIRS, PAIR_MAP) ) 
			 break ;
	}
		

	//////////////////////////////////////////////////
	//
	// Generate A matrix, b vector, for least-squares:
	//
	//////////////////////////////////////////////////

	cout.precision(16);	// WE SHOULD AT LEAST MOVE THIS TO SOMEWHERE MORE REASONABLE.. LIKE THE SECTION WHERE THE OUTPUT IS ACTUALLY PRINTED
	
	#if VERBOSITY == 1
		if ( RANK == 0 ) cout << "...Populating the matricies for A, Coulomb forces, and overbonding..." << endl << endl;
	#endif
		

	double 	NEIGHBOR_PADDING = 0.3;

	int istart, iend ;

	// Each processor only calculates certain frames.
	divide_atoms(istart, iend, CONTROLS.NFRAMES) ;
	
	if((CONTROLS.FIT_STRESS  || CONTROLS.FIT_STRESS_ALL) && CONTROLS.NSTRESS == -1)
		CONTROLS.NSTRESS = CONTROLS.NFRAMES;
		
	if(CONTROLS.FIT_ENER && CONTROLS.NENER == -1)
		CONTROLS.NENER = CONTROLS.NFRAMES;		
	
	bool DO_STRESS     = CONTROLS.FIT_STRESS;    
	bool DO_STRESS_ALL = CONTROLS.FIT_STRESS_ALL;
	bool DO_ENER       = CONTROLS.FIT_ENER;

	for(int i= istart ; i <= iend ; i++)
	{
		// Only include stress tensor data for first NSTRESS frames..
		
		if(i >= CONTROLS.NSTRESS)
		{
			CONTROLS.FIT_STRESS     = false;	
			CONTROLS.FIT_STRESS_ALL = false;
		}
		
		if(i >= CONTROLS.NENER)
		{
			CONTROLS.FIT_ENER     = false;	
		}		
	
		// This output is specific to the number of processors.
		
		if(NPROCS==1)
			cout << "	Processing frame: " << setw(5) << i+1 << " of: " << CONTROLS.NFRAMES << endl;

		// Use very little padding because we will update neighbor list for every frame.
		
		NEIGHBOR_LIST.INITIALIZE(TRAJECTORY[i], NEIGHBOR_PADDING);
		NEIGHBOR_LIST.DO_UPDATE(TRAJECTORY[i], CONTROLS);		

		ZCalc_Deriv(CONTROLS, ATOM_PAIRS, TRIPS, QUADS, TRAJECTORY[i], A_MATRIX[i], COULOMB_FORCES[i], CONTROLS.N_LAYERS, CONTROLS.USE_3B_CHEBY, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST);

		if ( CONTROLS.IF_SUBTRACT_COORD ) // Subtract over-coordination forces from force to be output.
			SubtractCoordForces(TRAJECTORY[i], false, P_OVER_FORCES[i],  ATOM_PAIRS, PAIR_MAP, NEIGHBOR_LIST, true);	
		
		if (CONTROLS.IF_SUBTRACT_COUL) 
			SubtractEwaldForces(TRAJECTORY[i], NEIGHBOR_LIST, CONTROLS) ;

		if ( CONTROLS.FIT_POVER )	// Fit the overcoordination parameter.
			SubtractCoordForces(TRAJECTORY[i], true, P_OVER_FORCES[i],  ATOM_PAIRS, PAIR_MAP, NEIGHBOR_LIST, true);	
			
	}

	// Because we need to know later whether stress/energy data were included:
	
	CONTROLS.FIT_STRESS     = DO_STRESS;    
	CONTROLS.FIT_STRESS_ALL = DO_STRESS_ALL;
	CONTROLS.FIT_ENER       = DO_ENER;   
	
	
	#if VERBOSITY == 1
		if ( RANK == 0 ) 
		{
			cout << "...matrix population complete: "  << endl << endl;
			cout << "Printing matricies..." << endl;
		}
	#endif
	
	//////////////////////////////////////////////////
	//
	// Print out least squares matrix:
	// Ax=b.
	// A printed to "A.txt", b printed to "b.txt"
	// use these files for python SVD routine.
	//
	//////////////////////////////////////////////////
	
	// Lets recap. What does any given row of the A file look like???
	// Rows of A have the following order.
	// A's first dimension is the number of frames.  2nd dimension is the number of atoms.
	// 3rd dimension is the number of parameters.
	//
	// parameters in A are ordered as follows:
	// For x forces:
	//	short-range 2-body and 3-body interaction
	//      charge pair parameters (if any).
	//      linear over-coordination parameter (if used)
	// 
	// Same pattern repeated for y and z forces.
	// charge constraints.
	// (LEF)
		
	char nameA[20] ;
	char nameB[20] ;
	char nameBlab[20] ;

	// Label output files by the processor rank
	sprintf(nameA, "A.%04d.txt", RANK) ;
	sprintf(nameB, "b.%04d.txt", RANK) ;
	sprintf(nameBlab, "b-labeled.%04d.txt", RANK) ;

	ofstream fileA(nameA);
	ofstream fileb(nameB);
	ofstream fileb_labeled(nameBlab);

	fileA.precision(16);	//  Reduced precision to 6 for code testing.
	fileA << std::scientific ;

	fileb.precision(16);	//  Usual precision set to 16.
	fileb << std::scientific ;

	if ( RANK == 0 ) cout << "	...A matrix length: " << A_MATRIX.size() << endl << endl;
	  
	// Print only the assigned frames.
	for(int N= istart ; N <= iend ;N++)
	{
		for(int a=0;a<A_MATRIX[N].size();a++)
		{	
			// Check if we need to exclude some tensor data from the A and b text files.
			if((CONTROLS.FIT_STRESS  || CONTROLS.FIT_STRESS_ALL) && N >= CONTROLS.NSTRESS)
			{
				if((a==TRAJECTORY[N].ATOMS) || (a==TRAJECTORY[N].ATOMS+1) ||  (a==TRAJECTORY[N].ATOMS+2))
					continue;
			}
		
			// Check if we need to exclude some energy data from the A and b text files.
			if(CONTROLS.FIT_ENER && N >= CONTROLS.NENER)
			{
				if(a==A_MATRIX[N].size()-1)
					continue;
			}		
		
			
			// Print Afile: .../////////////// -- For X
		  
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
				fileA << A_MATRIX[N][a][n].X  << "   ";
			  
			if ( CONTROLS.FIT_COUL ) 
				for(int i=0; i<COULOMB_FORCES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << COULOMB_FORCES[N][i][a].X << "   ";

			if ( CONTROLS.FIT_POVER ) 
				fileA << " " << P_OVER_FORCES[N][a].X;
			
			add_col_of_ones(a,DO_ENER, DO_STRESS, DO_STRESS_ALL, TRAJECTORY[N].ATOMS, fileA);
			  
			fileA << endl;	
		  
			// Print Afile: .../////////////// -- For Y
		  
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
				fileA << A_MATRIX[N][a][n].Y  << "   ";
			  
			if ( CONTROLS.FIT_COUL ) 
				for(int i=0; i<COULOMB_FORCES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << COULOMB_FORCES[N][i][a].Y << "   ";
			  
			if ( CONTROLS.FIT_POVER ) 
				fileA << " " << P_OVER_FORCES[N][a].Y;
			  
			add_col_of_ones(a,DO_ENER, DO_STRESS, DO_STRESS_ALL, TRAJECTORY[N].ATOMS, fileA);
			  
			fileA << endl;	


			// Print Afile: .../////////////// -- For Z
		  
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
				fileA << A_MATRIX[N][a][n].Z  << "   ";
			  
			if ( CONTROLS.FIT_COUL ) 
				for(int i=0; i<COULOMB_FORCES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << COULOMB_FORCES[N][i][a].Z << "   ";
			  
			if ( CONTROLS.FIT_POVER ) 
				fileA << " " << P_OVER_FORCES[N][a].Z;
			  
			add_col_of_ones(a,DO_ENER, DO_STRESS, DO_STRESS_ALL, TRAJECTORY[N].ATOMS, fileA);
			  
			fileA << endl;		
			
			// Print Bfile: ...
			
			if(N>CONTROLS.NFRAMES-1) // then this is the 3b histogram stuff. All b values (energies) should be zero. -- assumes we're using 3b!!
			{
				if(TRIPS.VEC[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::CUBIC)
				{	
					fileb << 1500.0 << endl;
					fileb << 1500.0 << endl;
					fileb << 1500.0 << endl;
				}
				else		
				{
					fileb << 0.0 << endl;
					fileb << 0.0 << endl;
					fileb << 0.0 << endl;
				}
			}
			else
			{
			
				if(a<TRAJECTORY[N].ATOMS)
				{
					fileb << TRAJECTORY[N].FORCES[a].X << endl;
					fileb << TRAJECTORY[N].FORCES[a].Y << endl;
					fileb << TRAJECTORY[N].FORCES[a].Z << endl;
			
					fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].X << endl;
					fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].Y << endl;
					fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].Z << endl;
				}
				else if ((a==TRAJECTORY[N].ATOMS)&& (CONTROLS.FIT_STRESS))
				{
					// Convert from GPa to internal units to match A-matrix elements

					fileb << TRAJECTORY[N].STRESS_TENSORS.X/GPa << endl;
					fileb << TRAJECTORY[N].STRESS_TENSORS.Y/GPa << endl;
					fileb << TRAJECTORY[N].STRESS_TENSORS.Z/GPa << endl;
			
					fileb_labeled << "-1 " <<  TRAJECTORY[N].STRESS_TENSORS.X/GPa << endl;
					fileb_labeled << "-1 " <<  TRAJECTORY[N].STRESS_TENSORS.Y/GPa << endl;
					fileb_labeled << "-1 " <<  TRAJECTORY[N].STRESS_TENSORS.Z/GPa << endl;
				}
				else if ((a>=TRAJECTORY[N].ATOMS)&& (a<(TRAJECTORY[N].ATOMS+3)) && CONTROLS.FIT_STRESS_ALL)
				{
					if(a==TRAJECTORY[N].ATOMS)
					{
						// Account for the symmetry of the off-diagonal (deviatoric) components
				
						fileb << TRAJECTORY[N].STRESS_TENSORS_X.X/GPa << endl;
						fileb << TRAJECTORY[N].STRESS_TENSORS_X.Y/GPa << endl;
						fileb << TRAJECTORY[N].STRESS_TENSORS_X.Z/GPa << endl;
					
						fileb << TRAJECTORY[N].STRESS_TENSORS_X.Y/GPa << endl; // Symmetry - this is just Y.X
						fileb << TRAJECTORY[N].STRESS_TENSORS_Y.Y/GPa << endl;
						fileb << TRAJECTORY[N].STRESS_TENSORS_Y.Z/GPa << endl;
					
						fileb << TRAJECTORY[N].STRESS_TENSORS_X.Z/GPa << endl; // Symmetry - this is just Z.X
						fileb << TRAJECTORY[N].STRESS_TENSORS_Y.Z/GPa << endl; // Symmetry - this is just Z.Y
						fileb << TRAJECTORY[N].STRESS_TENSORS_Z.Z/GPa << endl;
					
						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_X.X/GPa << endl;
						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_X.Y/GPa << endl;
						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_X.Z/GPa << endl;

						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_X.Y/GPa << endl; // Symmetry - this is just Y.X
						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_Y.Y/GPa << endl;
						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_Y.Z/GPa << endl;

						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_X.Z/GPa << endl; // Symmetry - this is just Z.X
						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_Y.Z/GPa << endl; // Symmetry - this is just Z.Y
						fileb_labeled << "-1 " << TRAJECTORY[N].STRESS_TENSORS_Z.Z/GPa << endl;
					}

				}
				else
				{			
					fileb << TRAJECTORY[N].QM_POT_ENER << endl;
					fileb << TRAJECTORY[N].QM_POT_ENER << endl;
					fileb << TRAJECTORY[N].QM_POT_ENER << endl;
			
					fileb_labeled << "+1 " <<  TRAJECTORY[N].QM_POT_ENER << endl;
					fileb_labeled << "+1 " <<  TRAJECTORY[N].QM_POT_ENER << endl;
					fileb_labeled << "+1 " <<  TRAJECTORY[N].QM_POT_ENER << endl;	
				}	
			}
		}  	
    }
	
	// A and B file charge constraints: generalized
	// 
	// The order that charge constraints are printed needs to match
	// the order atom pairs are expected...

	if ( CONTROLS.FIT_COUL 
		  && CHARGE_CONSTRAINTS.size()>0 
		  && RANK == NPROCS - 1) 
	{
		for(int i=0; i<CHARGE_CONSTRAINTS.size(); i++)
		{
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << "0.0 ";
			
			for(int j=0; j<ATOM_PAIRS.size(); j++)
				for(int k=0; k<CHARGE_CONSTRAINTS.size()+1; k++) // +1 because we n_constr = npairs-1
					if(CHARGE_CONSTRAINTS[i].PAIRTYPE_IDX[k] == j)
						fileA << CHARGE_CONSTRAINTS[i].CONSTRAINTS[k] << " ";
			
			if ( CONTROLS.FIT_POVER ) 
				fileA  << " 0.0 ";
			
			fileA << endl;	
			
			fileb << CHARGE_CONSTRAINTS[i].FORCE << endl;	
		}		
	}
	
	fileA.close();
	fileb.close();
	fileb_labeled.close();

	#ifdef USE_MPI
		// Make sure that every process has closed its files.
		MPI_Barrier(MPI_COMM_WORLD) ;
	#endif

	if ( RANK == 0 ) {
		// Serialize into a single A and b file for now.
		// Could make the SVD program read multiple files.
		system("cat A.[0-9]*.txt > A.txt") ;
		system("rm A.[0-9]*.txt") ;
		system("cat b.[0-9]*.txt > b.txt") ;
		system("rm b.[0-9]*.txt") ;
		system("cat b-labeled.[0-9]*.txt > b-labeled.txt") ;
		system("rm b-labeled.[0-9]*.txt") ;
	}

	  	  
	//////////////////////////////////////////////////
	//
	// Print out the params file header
	//
	//////////////////////////////////////////////////	   
	
	ofstream header;
	header.open("params.header");
	
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
	
	CONTROLS.USE_POVER = false;
	for(int i=0; i<ATOM_PAIRS.size(); i++)
		if(ATOM_PAIRS[i].USE_OVRPRMS)
			CONTROLS.USE_POVER = true;
	
	if(CONTROLS.USE_POVER)
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

	  string chtype = Cheby::get_trans_string(ATOM_PAIRS[i].CHEBY_TYPE) ;
		header << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
			 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
			 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
			 << setw(16) << left << ATOM_PAIRS[i].S_MINIM
			 << setw(16) << left << ATOM_PAIRS[i].S_MAXIM							 
			 << setw(16) << left << ATOM_PAIRS[i].S_DELTA
				 << setw(16) << left << chtype ;
		if(ATOM_PAIRS[i].CHEBY_TYPE == Cheby_trans::MORSE )
			header << setw(16) << left << ATOM_PAIRS[i].LAMBDA << endl; 
		else
			header << endl;
	}
	

	if(CONTROLS.USE_POVER)
	{
		header << endl;
 		header << "!# PAIRIDX #	";
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
	if(ATOM_PAIRS[0].SNUM_3B_CHEBY> 0 || ATOM_PAIRS[0].SNUM_4B_CHEBY> 0)
	  TRIPS.print_fcut_header(header) ;
		
	if(ATOM_PAIRS[0].CUBIC_SCALE != 1.0)
		header << endl << "PAIR CHEBYSHEV CUBIC SCALING: " << ATOM_PAIRS[0].CUBIC_SCALE << endl;
	
	// Print out special cutoffs 
	TRIPS.print_special(header) ;
	QUADS.print_special(header) ;
	 

	// Print out cluster parameters into the header.
	TRIPS.print_header(header,3,CONTROLS.CHEBY_3B_ORDER) ;
	QUADS.print_header(header,4,CONTROLS.CHEBY_4B_ORDER) ;

	header << endl;
	header.close();
	
	//////////////////////////////////////////////////
	//
	// Print out the pair/triplet type  map file
	//
	//////////////////////////////////////////////////	  	

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

	//////////////////////////////////////////////////
	//
	// Print out bonding statistics.
	//
	//////////////////////////////////////////////////	  

	#if VERBOSITY == 1
	print_bond_stats(ATOM_PAIRS, TRIPS, QUADS, CONTROLS.USE_3B_CHEBY, CONTROLS.USE_4B_CHEBY) ;
	#endif
	  
return 0;		  
}




	//////////////////////////////////////////////////
	//
    // Function definitions
	//
	//////////////////////////////////////////////////



// Read program input from the file "splines_ls.in".
static void read_lsq_input(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, 
									CLUSTER_LIST &TRIPS, CLUSTER_LIST & QUADS, map<string,int> & PAIR_MAP, 
									map<int,string> & PAIR_MAP_REVERSE, vector<int> &INT_PAIR_MAP,
									vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, 
									NEIGHBORS & NEIGHBOR_LIST,
									vector<int>& TMP_ATOMTYPEIDX,
									vector<string>& TMP_ATOMTYPE)
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
	
	// Set some defaults
	
	double TMP_CHEBY_RANGE_LOW  = -1;
	double TMP_CHEBY_RANGE_HIGH =  1;
	
	CONTROLS.COMPARE_FORCE = true;	// is this variable really necessary for LSQ?
	CONTROLS.IS_LSQ        = true; 
	CONTROLS.CALL_EWALD    = false;
	CONTROLS.FIT_ENER      = false;
	CONTROLS.FIT_STRESS    = false;
	CONTROLS.FIT_STRESS_ALL= false;
	CONTROLS.NSTRESS       = -1;
	CONTROLS.NENER         = -1;
	
	NEIGHBOR_LIST.USE	   = true;
	
	while (FOUND_END == false)
	{
		getline(cin,LINE);


		if(LINE.find("# ENDFILE #") != string::npos)
		{
		
			// Set up the Cheby variables
			
			if(CONTROLS.USE_3B_CHEBY)
				TRIPS.build_cheby_vals(ATOM_PAIRS);

			if(CONTROLS.USE_4B_CHEBY)
				QUADS.build_cheby_vals(ATOM_PAIRS);
			
		
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
			

#if VERBOSITY == 1

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
			
#endif
				
				
			if (CONTROLS.USE_PARTIAL_CHARGES || CONTROLS.FIT_COUL)
				CONTROLS.CALL_EWALD = true;	
				
			FOUND_END = true;
			break;
				
		}
		
		else if(LINE.find("# WRAPTRJ #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.WRAP_COORDS = bool(LINE.data());
#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	# WRAPTRJ #: " << CONTROLS.WRAP_COORDS << endl;	
#endif
		}
		
		else if(LINE.find("# TRJFILE #") != string::npos)
		{
			cin >> CONTROLS.INFILE; cin.ignore();

#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	# TRJFILE #: " << CONTROLS.INFILE << endl;	
#endif
		}			
		
		else if(LINE.find("# NFRAMES #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.NFRAMES = int(atof(LINE.data()));
			
#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	# NFRAMES #: " << CONTROLS.NFRAMES << endl;			
#endif
		}
		
		else if(LINE.find("# NLAYERS #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.N_LAYERS = int(atof(LINE.data()));
			
#if VERBOSITY == 1
			if ( RANK == 0 ) 	cout << "	# NLAYERS #: " << CONTROLS.N_LAYERS << endl;
#endif
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
			
#if VERBOSITY == 1
			
			if ( RANK == 0 ) 
			{
				cout << "	# FITCOUL #: ";
			
				if (CONTROLS.FIT_COUL)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
			}
				
#endif
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
			
#if VERBOSITY == 1
			if ( RANK == 0 ) 
			{
				cout << "	# FITPOVR #: ";
			
				if (CONTROLS.FIT_POVER)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
			}
#endif
		}
		
		else if(LINE.find("# PAIRTYP #") != string::npos)
		{
			cin >> LINE; //cin.ignore();
			
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
			
#if VERBOSITY == 1
			if ( RANK == 0 )
			{
				cout << "	# PAIRTYP #: " << LINE;
				
				if(TEMP_TYPE != "DFTBPOLY")
					cout << " ....NOTE: Forces reported in units of kcal/(mol.A), potential energy in kcal/mol." << endl;
				else
					cout << " ....NOTE: Forces reported in atomic units." << endl;
			}
#endif
				
			if(TEMP_TYPE == "INVRSE_R")
			{
				cin >> LINE;
				cin.ignore();
				CONTROLS.INVR_PARAMS = int(atof(LINE.data()));
					
#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "	             " << "Will use the following number of inverse-r params: " << CONTROLS.INVR_PARAMS << endl;
#endif
					
			}

			if(TEMP_TYPE == "DFTBPOLY" || TEMP_TYPE == "CHEBYSHEV")
			{
				getline(cin,LINE);
				
				STREAM_PARSER.str(LINE);
				
				STREAM_PARSER >> CONTROLS.CHEBY_ORDER;				
				
#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "	             " << "Will use 2-body order: " << CONTROLS.CHEBY_ORDER << endl;
#endif
				
				if(TEMP_TYPE == "CHEBYSHEV")
				{
					STREAM_PARSER >> CONTROLS.CHEBY_3B_ORDER;

#if VERBOSITY == 1
					if ( RANK == 0 ) cout << "	             " << "Will use 3-body order: " << CONTROLS.CHEBY_3B_ORDER << endl;
#endif
						
					if(CONTROLS.CHEBY_3B_ORDER>0)
						CONTROLS.USE_3B_CHEBY = true;
					
					STREAM_PARSER >> CONTROLS.CHEBY_4B_ORDER;

#if VERBOSITY == 1
					if ( RANK == 0 ) cout << "	             " << "Will use 4-body order: " << CONTROLS.CHEBY_4B_ORDER << endl;
#endif
						
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
					
#if VERBOSITY == 1
					if ( RANK == 0 ) cout << "	             " << "Will transform Chebyshev pair distances to range " << TMP_CHEBY_RANGE_LOW << " to " << TMP_CHEBY_RANGE_HIGH << endl;
#endif
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

		else if( (TEMP_TYPE == "CHEBYSHEV") && (LINE.find("# CHBTYPE #") != string::npos))
		{

		  cin >> LINE ;
		  vector<string> tokens ;
		  if ( parse_space(LINE,tokens) >= 1 ) 
				CONTROLS.CHEBY_TYPE = Cheby::get_trans_type(tokens[0]);
		  else 
				EXIT_MSG("BAD CHBTYPE" + LINE) ;

#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	# CHBTYPE #: " << tokens[0] << endl;	
#endif
		}		
		
		/////////////////////////////////////////////////////////////////////
		// Read the topology part. For now, ignoring index and atom types.. 
		// Assuming given as OO, OH, HH, as code expects... 
		// will need to be fixed later.
		/////////////////////////////////////////////////////////////////////
		
		else if(LINE.find("EXCLUDE 3B INTERACTION:")!= string::npos)
		{
		  TRIPS.read_exclude(cin, LINE) ;
		}	

		else if(LINE.find("EXCLUDE 4B INTERACTION:")!= string::npos)
		{
		  QUADS.read_exclude(cin, LINE) ;
		}	
		
		else if(LINE.find("# NATMTYP #")!= string::npos)
		{
			cin >> LINE; cin.ignore();
			CONTROLS.NATMTYP = int(atof(LINE.data()));
			
#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	# NATMTYP #: " << CONTROLS.NATMTYP << endl;	
#endif
			
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
					ATOM_PAIRS[i].PAIRTYP           = TEMP_TYPE ;
				}
			}	
			
			// Set up triplets
			
			NTRIP = factorial(CONTROLS.NATMTYP+3-1)/factorial(3)/factorial(CONTROLS.NATMTYP-1);
			TRIPS.allocate(NTRIP, 3, ATOM_PAIRS) ;
			
			// Set up quadruplets
			
			NQUAD = factorial(CONTROLS.NATMTYP+4-1)/factorial(4)/factorial(CONTROLS.NATMTYP-1);
			QUADS.allocate(NQUAD, 4, ATOM_PAIRS) ;
		}

		else if(LINE.find("# TYPEIDX #")!= string::npos)
		{
			
#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	# TYPEIDX #    # ATM_TYP #    # ATMCHRG #    # ATMMASS #" << endl;
#endif
			
			// Figure out the number of non-unique pairs

			TEMP_INT = CONTROLS.NATMTYP*CONTROLS.NATMTYP;
			
			SUM_OF_CHARGES = 0;
			
			TMP_ATOMTYPEIDX.resize(CONTROLS.NATMTYP);
			TMP_ATOMTYPE   .resize(CONTROLS.NATMTYP);
			
			for(int i=0; i<CONTROLS.NATMTYP; i++)
			{
				
				// Set the first atom pair types to be of type OO, HH, CC, etc...
				
				ATOM_PAIRS[i].PAIRTYP    = TEMP_TYPE;
				ATOM_PAIRS[i].PAIRIDX	 = i; 
				ATOM_PAIRS[i].CHEBY_TYPE = CONTROLS.CHEBY_TYPE;
				
				cin >> LINE ;
				TMP_ATOMTYPEIDX[i] = stoi(LINE) - 1 ;
				ATOM_PAIRS[i].ATM1TYPE_IDX = TMP_ATOMTYPEIDX[i] ;

				if ( TMP_ATOMTYPEIDX[i] < 0 || TMP_ATOMTYPEIDX[i] >= CONTROLS.NATMTYP )
				  EXIT_MSG("Bad atom index: " + LINE ) ;

				cin >> ATOM_PAIRS[i].ATM1TYP >> LINE;
				
				TMP_ATOMTYPE[i]    = ATOM_PAIRS[i].ATM1TYP;
				
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
				ATOM_PAIRS[i].ATM2TYPE_IDX = ATOM_PAIRS[i].ATM1TYPE_IDX ;

				ATOM_PAIRS[i].ATM2CHG = ATOM_PAIRS[i].ATM1CHG;
				ATOM_PAIRS[i].ATM2MAS = ATOM_PAIRS[i].ATM1MAS;
				
				ATOM_PAIRS[i].PRPR_NM =      ATOM_PAIRS[i].ATM1TYP;
				ATOM_PAIRS[i].PRPR_NM.append(ATOM_PAIRS[i].ATM2TYP);
				
#if VERBOSITY == 1
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
#endif
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
					ATOM_PAIRS[TEMP_INT].ATM1TYPE_IDX = ATOM_PAIRS[i].ATM1TYPE_IDX ;

					ATOM_PAIRS[TEMP_INT].ATM2TYP = ATOM_PAIRS[j].ATM1TYP;
					ATOM_PAIRS[TEMP_INT].ATM2TYPE_IDX = ATOM_PAIRS[j].ATM1TYPE_IDX ;

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
							PAIR_MAP_REVERSE.insert(make_pair(k,TEMP_STR));
						}
					}
				}
			}
			
			//cout << "Made the following pairs: " << endl;
			//for(map<string,int>::iterator i=PAIR_MAP.begin(); i!=PAIR_MAP.end(); i++)
			//cout << i->second << " " << i->first << endl;			

#if VERBOSITY == 1						
			if ( RANK == 0 ) 
			{
				cout << endl;
				cout << "	The following unique pair types have been identified:" << endl;
				for(int i=0;i<NPAIR; i++)
					cout << "		" << ATOM_PAIRS[i].PAIRIDX << "  " << ATOM_PAIRS[i].ATM1TYP << " " << ATOM_PAIRS[i].ATM2TYP << endl;
			}
#endif
					

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

#if VERBOSITY == 1			
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
#endif
				
			bool PRINT_OVR = false;
			
			for(int i=0; i<NPAIR; i++)
			{
#if VERBOSITY == 1	
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
#endif
					 
				if(ATOM_PAIRS[i].USE_OVRPRMS)
					PRINT_OVR = true;
			}
		
			// If overbonding parameters are provided, and they are not requested to be fit, 
			// subtract thier contribution before generating A matrix
			
			if(PRINT_OVR && !CONTROLS.FIT_POVER)
				CONTROLS.IF_SUBTRACT_COORD = true;

			if(PRINT_OVR)
			{
#if VERBOSITY == 1
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
#endif

				for(int i=0; i<NPAIR; i++)
				{
#if VERBOSITY == 1					
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
#endif													
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

			build_int_pair_map(CONTROLS.NATMTYP, TMP_ATOMTYPE, TMP_ATOMTYPEIDX, PAIR_MAP, INT_PAIR_MAP) ;

			for ( int i = 0 ; i < ATOM_PAIRS.size() ; i++ ) {
				if ( ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) 
					ATOM_PAIRS[i].set_cheby_vals() ;
			}
			if(CONTROLS.USE_3B_CHEBY)
			{
				// Generate unique triplets
				TRIPS.build_all(CONTROLS.CHEBY_3B_ORDER, ATOM_PAIRS, PAIR_MAP,TMP_ATOMTYPE,TMP_ATOMTYPEIDX) ; 
				TRIPS.print(false) ;
			}	
			
			if(CONTROLS.USE_4B_CHEBY)	// WORKS
			{
				// Generate unique quadruplets and thier corresponding sets of powers
				QUADS.build_all(CONTROLS.CHEBY_4B_ORDER, ATOM_PAIRS, PAIR_MAP, TMP_ATOMTYPE, TMP_ATOMTYPEIDX) ;
				QUADS.print(false) ;
			}			

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
#if VERBOSITY == 1			
			if ( RANK == 0 ) cout << endl << "	Attempting to read " << NPAIR-1 << " charge constraints...:" << endl; 
#endif	
			
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

#if VERBOSITY == 1			
				if ( RANK == 0 ) 
				{
					cout << "		" << i+1 << "	 ";
					for(int j=0; j<NPAIR; j++)
						cout << CHARGE_CONSTRAINTS[i].PAIRTYPE[j] << " (" << CHARGE_CONSTRAINTS[i].PAIRTYPE_IDX[j] << ") ";
					for(int j=0; j<NPAIR; j++)
						cout << CHARGE_CONSTRAINTS[i].CONSTRAINTS[j] << " ";
					cout << CHARGE_CONSTRAINTS[i].FORCE << endl;	
				}
#endif	
					
			}
			
			if ( RANK == 0 ) cout << endl;
		}
		
		else if(LINE.find("SPECIAL 3B S_MAXIM:") != string::npos)
		{
		  NEIGHBOR_LIST.MAX_CUTOFF_3B = TRIPS.read_cutoff_params(cin, LINE, "S_MAXIM", ATOM_PAIRS, PAIR_MAP) ;
		}
		
		else if(LINE.find("SPECIAL 3B S_MINIM:") != string::npos)
		{
		  TRIPS.read_cutoff_params(cin, LINE, "S_MINIM", ATOM_PAIRS, PAIR_MAP) ;
		}
		
		else if(LINE.find("SPECIAL 4B S_MAXIM:") != string::npos)
		{
		  NEIGHBOR_LIST.MAX_CUTOFF_4B = QUADS.read_cutoff_params(cin, LINE, "S_MAXIM", ATOM_PAIRS, PAIR_MAP) ;
		}

		else if(LINE.find("SPECIAL 4B S_MINIM:") != string::npos)
		{
		  QUADS.read_cutoff_params(cin, LINE, "S_MINIM", ATOM_PAIRS, PAIR_MAP) ;
		}
		
		else if ( (LINE.find("# FCUTTYP #") != string::npos) && ((CONTROLS.CHEBY_3B_ORDER >0)||(CONTROLS.CHEBY_4B_ORDER>0)))
		{

			cin >> TEMP_TYPE;
			
			for(int i=0; i<ATOM_PAIRS.size(); i++)
				ATOM_PAIRS[i].FORCE_CUTOFF.set_type(TEMP_TYPE);
			
			if (CONTROLS.USE_3B_CHEBY)
				for(int i=0; i<TRIPS.VEC.size(); i++)
					TRIPS.VEC[i].FORCE_CUTOFF.set_type(TEMP_TYPE);
					
			if (CONTROLS.USE_4B_CHEBY)
				for(int i=0; i<QUADS.VEC.size(); i++)
					QUADS.VEC[i].FORCE_CUTOFF.set_type(TEMP_TYPE);

#if VERBOSITY == 1
			if ( RANK == 0 ) 
				cout << "	# FCUTTYP #: " << TEMP_TYPE << "	... for all Chebyshev interactions" << endl;	
#endif
			
			double TMP_STEEPNESS = 0.0;
			double TMP_OFFSET    = 0.0;
			double TMP_HEIGHT    = 0.0;
			
			if(TEMP_TYPE=="TERSOFF")
			{
				cin >> TMP_OFFSET;
				
#if VERBOSITY == 1
				if ( RANK == 0 ) 
					cout << "			With an offset of: " << fixed << setprecision(3) << TMP_OFFSET*100.0 << "% of the pair outer cutoffs" << endl;      
#endif			
			}
			
			if(TEMP_TYPE=="SIGMOID" || TEMP_TYPE=="CUBSIG" || TEMP_TYPE=="CUBESTRETCH" || TEMP_TYPE == "SIGFLT")
			{
				cin >> TMP_STEEPNESS;
				cin >> TMP_OFFSET;
				
				if(TEMP_TYPE == "SIGFLT")
					cin >> TMP_HEIGHT;
					
#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "			With steepness, offset, and height of: "
															<< fixed << setprecision(3) 
															<< TMP_STEEPNESS << ",   " 
															<< TMP_OFFSET    << ", and  " 
															<< TMP_HEIGHT    << endl;   
#endif					
			}
			cin.ignore();
			
			
			if (ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV")
			{
				for(int i=0; i<ATOM_PAIRS.size(); i++)
				{
					ATOM_PAIRS[i].FORCE_CUTOFF.BODIEDNESS = 2;
					ATOM_PAIRS[i].FORCE_CUTOFF.STEEPNESS  = TMP_STEEPNESS; 
					ATOM_PAIRS[i].FORCE_CUTOFF.OFFSET     = TMP_OFFSET;
					ATOM_PAIRS[i].FORCE_CUTOFF.HEIGHT     = TMP_HEIGHT;
				}
			}
			if (CONTROLS.USE_3B_CHEBY)
			{
				for(int i=0; i<TRIPS.VEC.size(); i++)
				{
					TRIPS.VEC[i].FORCE_CUTOFF.BODIEDNESS = 3;
					TRIPS.VEC[i].FORCE_CUTOFF.STEEPNESS  = TMP_STEEPNESS; 
					TRIPS.VEC[i].FORCE_CUTOFF.OFFSET     = TMP_OFFSET;
					TRIPS.VEC[i].FORCE_CUTOFF.HEIGHT     = TMP_HEIGHT;
				}			
			}
			if (CONTROLS.USE_4B_CHEBY)
			{
				for(int i=0; i<QUADS.VEC.size(); i++)
				{
					QUADS.VEC[i].FORCE_CUTOFF.BODIEDNESS = 4;
					QUADS.VEC[i].FORCE_CUTOFF.STEEPNESS  = TMP_STEEPNESS; 
					QUADS.VEC[i].FORCE_CUTOFF.OFFSET     = TMP_OFFSET;
					QUADS.VEC[i].FORCE_CUTOFF.HEIGHT     = TMP_HEIGHT;
				}			
			}

		}	
		
		
	}	
	
#if VERBOSITY == 1			
	if ( RANK == 0 ) cout << endl << "Note: Will use cubic scaling of: " << ATOM_PAIRS[0].CUBIC_SCALE << endl << endl;; // All types use same scaling
#endif	
}

static void print_bond_stats(vector<PAIRS> &ATOM_PAIRS, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, bool use_3b_cheby, bool use_4b_cheby)
// Print statistics on bonding.
{

	if ( RANK == 0 ) cout << "	Minimum distances between atoms: (Angstr.)" << endl;

	int NPAIR = ATOM_PAIRS.size() ;

		for (int k=0; k<NPAIR; k++)
		{
			double sum = 0.0 ;
#ifdef USE_MPI
			MPI_Reduce(&(ATOM_PAIRS[k].MIN_FOUND_DIST), &sum, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
#else
			sum = ATOM_PAIRS[k].MIN_FOUND_DIST ;
#endif
			if ( RANK == 0 ) 
			{
				cout << "		" << k << "	" << ATOM_PAIRS[k].ATM1TYP << 
					" " << ATOM_PAIRS[k].ATM2TYP << "	" << fixed << setprecision(3) << sum << endl;

			}

		}

		if ( RANK == 0 ) cout << "	Total number of configurations contributing to each pair type:" << endl;
		
		for (int k = 0 ; k < NPAIR ; k++) 
		{
			int sum = 0.0 ;
#ifdef USE_MPI
			MPI_Reduce(&(ATOM_PAIRS[k].N_CFG_CONTRIB), &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
			sum = ATOM_PAIRS[k].N_CFG_CONTRIB ;
#endif
			if ( RANK == 0 ) cout << "		" << k << "	" 					
										 << ATOM_PAIRS[k].ATM1TYP << " " 
										 << ATOM_PAIRS[k].ATM2TYP << " " 
										 << sum << endl;
		}

		if( use_3b_cheby ) 
		{
		  TRIPS.print_min_distances() ;
		}

		if( use_4b_cheby ) 
		{
		  QUADS.print_min_distances() ;
		}
			
		if ( RANK == 0 ) cout << "...matrix printing complete: " << endl << endl;
}


static void add_col_of_ones(int item, bool DO_ENER, bool DO_STRESS, bool DO_STRESS_ALL, int NATOMS, ofstream & OUTFILE)
{
	// Determine if: Energies are being included in the fit
	// If so, is this A-matrix row is for a force or stress (then print an additional a 0.0)
	// or if it is for an energy (then print an additional 1.0)?
	
	if (DO_ENER)
	{
		if(DO_STRESS)
		{
			if(item<NATOMS+1)
				OUTFILE << " 0.0";
			else
				OUTFILE << " 1.0";
		}
		else if(DO_STRESS_ALL)
		{
			if(item<NATOMS+3)
				OUTFILE << " 0.0";
			else
				OUTFILE << " 1.0";					
		}
		else
		{
			if(item<NATOMS)
				OUTFILE << " 0.0";
			else
				OUTFILE << " 1.0";
		}
			
	}
}
