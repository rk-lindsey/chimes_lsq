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


using namespace std;

#ifndef VERBOSITY 
	#define VERBOSITY 1 
#endif

// For now, keep the lsq code in serial.

#ifdef USE_MPI
#include <mpi.h>
#endif 

static void read_lsq_input(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, vector<TRIPLETS> & PAIR_TRIPLETS, vector<QUADRUPLETS> & PAIR_QUADRUPLETS, map<string,int> & PAIR_MAP, map<int,string> & PAIR_MAP_REVERSE, map<string,int> & TRIAD_MAP, map<int,string> & TRIAD_MAP_REVERSE, map<string,int> & QUAD_MAP, map<int,string> & QUAD_MAP_REVERSE, vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, NEIGHBORS & NEIGHBOR_LIST);

static void print_bond_stats(vector<PAIRS> &ATOM_PAIRS, vector<TRIPLETS> &PAIR_TRIPLETS, vector<QUADRUPLETS> &PAIR_QUADRUPLETS, bool use_3b_cheby, bool use_4b_cheby);

void store_quad_permutations(QUADRUPLETS &quad, vector<int> &unsorted_powers)  ;
void permute_atom_indices(int idx, vector<string> names, QUADRUPLETS &quad, vector<int> &unsorted_powers, vector<int> perm, int unique_index,
								  int equiv_index) ;
bool operator==(const vector<int>& lhs, const vector<int>& rhs)  ;

int factorial(int input)
{
	int result = 1;
	for(int i=input; i>0; i--)
		result *= i;
	
	return result;
}

string FULL_FILE_3B;		// Global variables declared as externs in functions.h, and declared in functions.C
string SCAN_FILE_3B;
string SCAN_FILE_2B;

// Global variables declared as externs in functions.h, and declared in functions.C -- MPI calculations.   
 
int NPROCS;		// Number of processors
int RANK;		// Index of current processor

// Define my new integer maps.. Only needed here to prevent compiler errors related 
// to functions.C... only used in MD part right now

vector<int>	INT_PAIR_MAP;
vector<int>	INT_TRIAD_MAP;

// For 4-body interactions, these are used for both the lsq and md parts:

map<int,int>   INT_QUAD_MAP;			// maps for collections of 4 atoms
map<int,int>   INT_QUAD_MAP_REVERSE;	// corresponding reverse-maps
vector<int>    TMP_ATOMTYPEIDX;			// Used to construct the quadruplet type index
vector<string> TMP_ATOMTYPE;			// Used to construct the quadruplet type index

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
	vector<TRIPLETS> 	PAIR_TRIPLETS;		// Will store relevant info regarding atom interaction triplet types.. i.e. { [OO,OO,OO], [OO,HH,HH], ... }
	vector<QUADRUPLETS> PAIR_QUADRUPLETS;	// Will store relevant info regarding atom interaction quadruplet types.. note that a set of 4 atoms are described by 6 pairs
	vector<FRAME> 		TRAJECTORY;			// Stores the trajectory information... box dimensions, atom types, coordinates, and forces
	NEIGHBORS        	 NEIGHBOR_LIST;		// Declare the class that will handle the neighbor list
	vector<int>   		ATOM_PAIR_TYPES;	// Fore use in double loop over atom pairs. Index corresponds to the overall double loop count. 
											// Provides an index that rells you the atom pair's ATOM_PAIRS's type index.. THIS IS FOR LOOPS OF TYPE
											// 	for(int a1=0;a1<nat-1;a1++)	for(int a2=a1+1;a2<nat;a2++)
	vector<int>   ATOM_PAIR_TYPES_ALL;		// THIS IS FOR LOOPS OF TYPE for(int ai=0; ai<nat; ai++), for(int ak=0; ak<nat; ak++)
	
	map<string,int> PAIR_MAP;			// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
	map<int,string> PAIR_MAP_REVERSE; 	// Input/output the resverse of PAIR_MAP
	
	map<string,int> TRIAD_MAP;			// Input is three of any ATOM_PAIRS.PRPR_NM pairs, output is a triplet type index
	map<int,string> TRIAD_MAP_REVERSE;	// Input/output is the reverse of TRIAD_MAP
	
	map<string,int> QUAD_MAP;			// maps for collections of 4 atoms
	map<int,string> QUAD_MAP_REVERSE;	// corresponding reverse-maps
	
	vector<CHARGE_CONSTRAINT> CHARGE_CONSTRAINTS;	// Specifies how we constrain charge fitting

	JOB_CONTROL 	CONTROLS;			// Will hold job controls shared by both lsq and md

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

	read_lsq_input(CONTROLS, ATOM_PAIRS, PAIR_TRIPLETS, PAIR_QUADRUPLETS, PAIR_MAP, PAIR_MAP_REVERSE, TRIAD_MAP, TRIAD_MAP_REVERSE, QUAD_MAP, QUAD_MAP_REVERSE, CHARGE_CONSTRAINTS, NEIGHBOR_LIST);

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
			ATOM_PAIRS[i].SNUM          = CONTROLS.CHEBY_ORDER;
			
			if (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" )
			{
				ATOM_PAIRS[i].SNUM_3B_CHEBY = CONTROLS.CHEBY_3B_ORDER;
				ATOM_PAIRS[i].SNUM_4B_CHEBY = CONTROLS.CHEBY_4B_ORDER;
			}
				
        }
		
		else if (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) 	// Set the distance transformation type
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
		
		for(int i=0; i<PAIR_TRIPLETS.size(); i++)
			CONTROLS.NUM_3B_CHEBY += PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS;

#if VERBOSITY == 1
		if ( RANK == 0 ) cout << "The number of three-body Chebyshev parameters is: " << CONTROLS.NUM_3B_CHEBY << endl;
#endif
	}
	
	if ( ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV" && CONTROLS.CHEBY_4B_ORDER  > 0 ) // All atoms types must use same potential, so check if 3b order is greater than zero 
	{
		CONTROLS.NUM_4B_CHEBY = 0;
		
		for(int i=0; i<PAIR_QUADRUPLETS.size(); i++)
			CONTROLS.NUM_4B_CHEBY += PAIR_QUADRUPLETS[i].N_TRUE_ALLOWED_POWERS;

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
	
	XYZ MAX_RMIN;
	MAX_RMIN.X = MAX_RMIN.Y = MAX_RMIN.Z = 0;
	
	for(int i=0; i<PAIR_TRIPLETS.size(); i++)
	{
		if(PAIR_TRIPLETS[i].S_MINIM_3B.X > MAX_RMIN.X)
			MAX_RMIN.X = PAIR_TRIPLETS[i].S_MINIM_3B.X;
		if(PAIR_TRIPLETS[i].S_MINIM_3B.Y > MAX_RMIN.Y)
			MAX_RMIN.Y = PAIR_TRIPLETS[i].S_MINIM_3B.Y;
		if(PAIR_TRIPLETS[i].S_MINIM_3B.Z > MAX_RMIN.Z) 
			MAX_RMIN.Z = PAIR_TRIPLETS[i].S_MINIM_3B.Z;
	}
			
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
	int             N_LAYER_ELEM = 0;		

	for (int i=0; i<CONTROLS.NFRAMES; i++)
	{
		// Read in line with the number of atoms
		
		TRAJ_INPUT >> TRAJECTORY[i].ATOMS;
		
		// Read in line with box dimenstions
		
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.X;
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.Y;
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.Z;
		
		bool IFANY = false;
	
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

				/*
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_Y.X;
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_Z.X;
				TRAJ_INPUT >> TRAJECTORY[i].STRESS_TENSORS_Z.Y;
				*/
				
				
		}
		if(CONTROLS.FIT_ENER) // We're actually fitting to relative *differences* in energy
		{
			TRAJ_INPUT >> TRAJECTORY[i].QM_POT_ENER;
			
			if(i>0)
				TRAJECTORY[i].QM_POT_ENER -= TRAJECTORY[0].QM_POT_ENER;
		}
			
		// Check that outer cutoffs do not exceed half of the boxlength
		// with consideration of layering
		
		for(int j=0; j<ATOM_PAIRS.size(); j++)
		{
			if( (  ATOM_PAIRS[j].S_MAXIM > 0.5* TRAJECTORY[i].BOXDIM.X * (2*CONTROLS.N_LAYERS +1)
				|| ATOM_PAIRS[j].S_MAXIM > 0.5* TRAJECTORY[i].BOXDIM.Y * (2*CONTROLS.N_LAYERS +1)
				|| ATOM_PAIRS[j].S_MAXIM > 0.5* TRAJECTORY[i].BOXDIM.Z * (2*CONTROLS.N_LAYERS +1) ))
			{
					#if WARN == TRUE
						if (isatty(fileno(stdout)) && RANK == 0)
						{
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "WARNING: Outer cutoff greater than half at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Pair type:            " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Boxlengths:           " << TRAJECTORY[i].BOXDIM.X << " " << TRAJECTORY[i].BOXDIM.Y << " " << TRAJECTORY[i].BOXDIM.Z << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Layers:               " << CONTROLS.N_LAYERS << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Effective boxlengths: " << TRAJECTORY[i].BOXDIM.X * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Y * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Z * (CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
						}
						else if ( RANK == 0 ) 
						{
							cout << "WARNING: Outer cutoff greater than half at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<endl;
							cout << "	Pair type:            " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
							cout << "	Boxlengths:           " << TRAJECTORY[i].BOXDIM.X << " " << TRAJECTORY[i].BOXDIM.Y << " " << TRAJECTORY[i].BOXDIM.Z << endl;
							cout << "	Layers:               " << CONTROLS.N_LAYERS << endl;
							cout << "	Effective boxlengths: " << TRAJECTORY[i].BOXDIM.X * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Y * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Z * (CONTROLS.N_LAYERS +1) << endl;
														
						}

					#else
						if (isatty(fileno(stdout)) && RANK == 0)
						{
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "ERROR: Outer cutoff greater than half at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Pair type:            " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Boxlengths:           " << TRAJECTORY[i].BOXDIM.X << " " << TRAJECTORY[i].BOXDIM.Y << " " << TRAJECTORY[i].BOXDIM.Z << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Layers:               " << CONTROLS.N_LAYERS << COUT_STYLE.ENDSTYLE << endl;
							cout << COUT_STYLE.RED << COUT_STYLE.BOLD << "	Effective boxlengths: " << TRAJECTORY[i].BOXDIM.X * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Y * (CONTROLS.N_LAYERS +1) << " " << TRAJECTORY[i].BOXDIM.Z * (CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
							
							exit(0);
						}
						else if ( RANK == 0 ) 
						{
							cout << "ERROR: Outer cutoff greater than half at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<endl;
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
	
	TRAJECTORY[0].QM_POT_ENER = 0;
	
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
		
		double tmp_max, tmp_min;
	
		bool ANY_ZERO = false;
		
		for (int i=0; i<PAIR_TRIPLETS.size(); i++)
		{

			PAIR_TRIPLETS[i].NBINS.X = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR1] ].NBINS.X;
			PAIR_TRIPLETS[i].NBINS.Y = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR1] ].NBINS.Y;
			PAIR_TRIPLETS[i].NBINS.Z = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR1] ].NBINS.Z;
		
			if(PAIR_TRIPLETS[i].NBINS.X == 0 || PAIR_TRIPLETS[i].NBINS.Y == 0 || PAIR_TRIPLETS[i].NBINS.X == 0)
			{
				ANY_ZERO = true;
			
				#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "Found at least one 3b-population histogram nbins = 0. Will not do ANY histogramming. " << endl;
				#endif
				
				break;
			}
		
			tmp_max = PAIR_TRIPLETS[i].S_MAXIM_3B.X;
			tmp_min = PAIR_TRIPLETS[i].S_MINIM_3B.X;
			if(tmp_min == -1)
				tmp_min = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR1] ].S_MINIM;
			if(tmp_max == -1)
				tmp_max = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR1] ].S_MAXIM;
			PAIR_TRIPLETS[i].BINWS.X = (tmp_max - tmp_min)/PAIR_TRIPLETS[i].NBINS.X;
		
			tmp_max = PAIR_TRIPLETS[i].S_MAXIM_3B.Y;
			tmp_min = PAIR_TRIPLETS[i].S_MINIM_3B.Y;
			if(tmp_min == -1)
				tmp_min = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR2] ].S_MINIM;
			if(tmp_max == -1)
				tmp_max = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR2] ].S_MAXIM;
			PAIR_TRIPLETS[i].BINWS.Y = (tmp_max - tmp_min)/PAIR_TRIPLETS[i].NBINS.Y;
		
			tmp_max = PAIR_TRIPLETS[i].S_MAXIM_3B.Z;
			tmp_min = PAIR_TRIPLETS[i].S_MINIM_3B.Z;
			if(tmp_min == -1)
				tmp_min = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR3] ].S_MINIM;
			if(tmp_max == -1)
				tmp_max = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[i].ATMPAIR3] ].S_MAXIM;
			PAIR_TRIPLETS[i].BINWS.Z = (tmp_max - tmp_min)/PAIR_TRIPLETS[i].NBINS.Z;
		
			PAIR_TRIPLETS[i].POP_HIST.resize(PAIR_TRIPLETS[i].NBINS.X);
		
			for(int x=0; x<PAIR_TRIPLETS[i].NBINS.X; x++)
			{
				PAIR_TRIPLETS[i].POP_HIST[x].resize(PAIR_TRIPLETS[i].NBINS.Y);
			
				for(int y=0; y<PAIR_TRIPLETS[i].NBINS.Y; y++)
				{
					PAIR_TRIPLETS[i].POP_HIST[x][y].resize(PAIR_TRIPLETS[i].NBINS.Z);

					for(int z=0; z<PAIR_TRIPLETS[i].NBINS.Z; z++)
					{
						PAIR_TRIPLETS[i].POP_HIST[x][y][z] = 0;
					}
				
				}
			}
			
			if ( RANK == 0 ) cout << "	" << i << " " << fixed << setprecision(1) << PAIR_TRIPLETS[i].NBINS.X << " " << PAIR_TRIPLETS[i].NBINS.Y << " " << PAIR_TRIPLETS[i].NBINS.Z << endl;
		}	
		
		#if VERBOSITY == 1
		if ( RANK == 0 ) 
		{
			if(!ANY_ZERO)
				cout << "...Initial histogram setup complete!" << endl << endl;
			else
				cout << endl << endl;
		}
		#endif
	}
		

	//////////////////////////////////////////////////
	//
	// Generate A matrix, b vector, for least-squares:
	//
	//////////////////////////////////////////////////

	// cout.precision(10);		// WE SHOULD AT LEAST MOVE THIS TO SOMEWHERE MORE REASONABLE.. LIKE THE SECTION WHERE THE OUTPUT IS ACTUALLY PRINTED
	cout.precision(16);
	
	#if VERBOSITY == 1
	if ( RANK == 0 ) cout << "...Populating the matricies for A, Coulomb forces, and overbonding..." << endl << endl;
	#endif
		

	int 	FRAME_MATRX_SIZE; 
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
		
		if(i >= CONTROLS.FIT_ENER)
		{
			CONTROLS.FIT_STRESS     = false;	
		}		
	
		// This output is specific to the number of processors.
		
		if(NPROCS==1)
			cout << "	Processing frame: " << setw(5) << i+1 << " of: " << CONTROLS.NFRAMES << endl;

		// Use very little padding because we will update neighbor list for every frame.
		
		NEIGHBOR_LIST.INITIALIZE(TRAJECTORY[i], NEIGHBOR_PADDING);
		NEIGHBOR_LIST.DO_UPDATE(TRAJECTORY[i], CONTROLS);		

		ZCalc_Deriv(CONTROLS, ATOM_PAIRS, PAIR_TRIPLETS, PAIR_QUADRUPLETS, TRAJECTORY[i], A_MATRIX[i], COULOMB_FORCES[i], CONTROLS.N_LAYERS, CONTROLS.USE_3B_CHEBY, PAIR_MAP, TRIAD_MAP, INT_QUAD_MAP, NEIGHBOR_LIST);
		
		if ( CONTROLS.IF_SUBTRACT_COORD ) // Subtract over-coordination forces from force to be output.
			SubtractCoordForces(TRAJECTORY[i], false, P_OVER_FORCES[i],  ATOM_PAIRS, PAIR_MAP, NEIGHBOR_LIST, true);	
		
		if (CONTROLS.IF_SUBTRACT_COUL) 
			SubtractEwaldForces(TRAJECTORY[i], NEIGHBOR_LIST, CONTROLS) ;

		if ( CONTROLS.FIT_POVER )	// Fit the overcoordination parameter.
			SubtractCoordForces(TRAJECTORY[i], true, P_OVER_FORCES[i],  ATOM_PAIRS, PAIR_MAP, NEIGHBOR_LIST, true);	
		
		// If we're including the energy in the fit, we need to subtract off the A matrix entries from the first frame from all the rest
		if((CONTROLS.FIT_ENER) &&(i>0))
		{
			FRAME_MATRX_SIZE = A_MATRIX[i].size()-1;
				
			for (int j=0; j<A_MATRIX[i][FRAME_MATRX_SIZE].size(); j++)
			{
				A_MATRIX[i][FRAME_MATRX_SIZE][j].X -= A_MATRIX[0][FRAME_MATRX_SIZE][j].X;
				A_MATRIX[i][FRAME_MATRX_SIZE][j].Y -= A_MATRIX[0][FRAME_MATRX_SIZE][j].Y;
				A_MATRIX[i][FRAME_MATRX_SIZE][j].Z -= A_MATRIX[0][FRAME_MATRX_SIZE][j].Z;
			}
		}	
	}

	// Because we need to know later whether stress/energy data were included:
	
	CONTROLS.FIT_STRESS     = DO_STRESS;    
	CONTROLS.FIT_STRESS_ALL = DO_STRESS_ALL;
	
	CONTROLS.FIT_ENER       = DO_ENER;   
	
	// If we're including the energy in the fit, we need to subtract off the A matrix entries from the first frame from itself too!
	if(CONTROLS.FIT_ENER)
	{
		FRAME_MATRX_SIZE = A_MATRIX[0].size()-1;
				
		for (int j=0; j<A_MATRIX[0][FRAME_MATRX_SIZE].size(); j++)
		{
			A_MATRIX[0][FRAME_MATRX_SIZE][j].X -= A_MATRIX[0][FRAME_MATRX_SIZE][j].X;
			A_MATRIX[0][FRAME_MATRX_SIZE][j].Y -= A_MATRIX[0][FRAME_MATRX_SIZE][j].Y;
			A_MATRIX[0][FRAME_MATRX_SIZE][j].Z -= A_MATRIX[0][FRAME_MATRX_SIZE][j].Z;
		}
	} 
	
	if(CONTROLS.USE_3B_CHEBY && PAIR_TRIPLETS[0].NBINS.X>0) // Set the 3b-population-histogram based constraints 
	{
		if ( RANK == 0 ) cout << "Setting constraints based on 3b-population histogram constraints " << endl << endl;
		ZCalc_3B_Cheby_Deriv_HIST(CONTROLS, ATOM_PAIRS, PAIR_TRIPLETS, A_MATRIX, PAIR_MAP, TRIAD_MAP);	
	}
	
	
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
	fileb.precision(16);	//  Usual precision set to 16.
	
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
			  
			fileA << endl;	
		  
			// Print Afile: .../////////////// -- For Y
		  
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
				fileA << A_MATRIX[N][a][n].Y  << "   ";
			  
			if ( CONTROLS.FIT_COUL ) 
				for(int i=0; i<COULOMB_FORCES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << COULOMB_FORCES[N][i][a].Y << "   ";
			  
			if ( CONTROLS.FIT_POVER ) 
				fileA << " " << P_OVER_FORCES[N][a].Y;
			  
			fileA << endl;	


			// Print Afile: .../////////////// -- For Z
		  
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
				fileA << A_MATRIX[N][a][n].Z  << "   ";
			  
			if ( CONTROLS.FIT_COUL ) 
				for(int i=0; i<COULOMB_FORCES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << COULOMB_FORCES[N][i][a].Z << "   ";
			  
			if ( CONTROLS.FIT_POVER ) 
				fileA << " " << P_OVER_FORCES[N][a].Z;
			  
			fileA << endl;		
			
			// Print Bfile: ...
			
			if(N>CONTROLS.NFRAMES-1) // then this is the 3b histogram stuff. All b values (energies) should be zero. -- assumes we're using 3b!!
			{
				if(PAIR_TRIPLETS[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::CUBIC)
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
				else if (CONTROLS.FIT_STRESS_ALL)
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
					}
				
				
					/*
				
					fileb << TRAJECTORY[N].STRESS_TENSORS_X.X/GPa << endl;
					fileb << TRAJECTORY[N].STRESS_TENSORS_X.Y/GPa << endl;
					fileb << TRAJECTORY[N].STRESS_TENSORS_X.Z/GPa << endl;
					
					fileb << TRAJECTORY[N].STRESS_TENSORS_Y.X/GPa << endl;
					fileb << TRAJECTORY[N].STRESS_TENSORS_Y.Y/GPa << endl;
					fileb << TRAJECTORY[N].STRESS_TENSORS_Y.Z/GPa << endl;
					
					fileb << TRAJECTORY[N].STRESS_TENSORS_Z.X/GPa << endl;
					fileb << TRAJECTORY[N].STRESS_TENSORS_Z.Y/GPa << endl;
					fileb << TRAJECTORY[N].STRESS_TENSORS_Z.Z/GPa << endl;
					
					*/
				}
				else
				{			
					fileb << TRAJECTORY[N].QM_POT_ENER << endl;
					fileb << TRAJECTORY[N].QM_POT_ENER << endl;
					fileb << TRAJECTORY[N].QM_POT_ENER << endl;
			
					fileb_labeled << "-1 " <<  TRAJECTORY[N].QM_POT_ENER << endl;
					fileb_labeled << "-1 " <<  TRAJECTORY[N].QM_POT_ENER << endl;
					fileb_labeled << "-1 " <<  TRAJECTORY[N].QM_POT_ENER << endl;	
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
	
	bool PRINT_OVR = false;
	int NPAIR =  ATOM_PAIRS.size();	
	
	header << endl << "ATOM PAIRS: " << NPAIR << endl << endl;
	
	header << "!# PAIRIDX #	";
	header << "# ATM_TY1 #	";
	header << "# ATM_TY1 #	";
	header << "# S_MINIM #	";
	header << "# S_MAXIM #	";
	header << "# S_DELTA #	";
	header << "# CHBDIST #	";	// how pair distance is transformed in cheby calc
	header << "# MORSE_LAMBDA #" << endl;

	for(int i=0; i<NPAIR; i++)
	{

		header << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
			 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
			 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
			 << setw(16) << left << ATOM_PAIRS[i].S_MINIM
			 << setw(16) << left << ATOM_PAIRS[i].S_MAXIM							 
			 << setw(16) << left << ATOM_PAIRS[i].S_DELTA
			 << setw(16) << left << ATOM_PAIRS[i].CHEBY_TYPE;
		if(ATOM_PAIRS[i].CHEBY_TYPE == "MORSE")
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
	
	if(ATOM_PAIRS[0].SNUM_3B_CHEBY> 0 || ATOM_PAIRS[0].SNUM_4B_CHEBY> 0)
	{
		header << endl << "FCUT TYPE: " << PAIR_TRIPLETS[0].FORCE_CUTOFF.to_string();
		
		if (PAIR_TRIPLETS[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::SIGMOID || PAIR_TRIPLETS[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::CUBSIG || PAIR_TRIPLETS[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::CUBESTRETCH || PAIR_TRIPLETS[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::SIGFLT)
			header << " " << PAIR_TRIPLETS[0].FORCE_CUTOFF.STEEPNESS << " " << PAIR_TRIPLETS[0].FORCE_CUTOFF.OFFSET;
		if(PAIR_TRIPLETS[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::SIGFLT)
			header << " " << PAIR_TRIPLETS[0].FORCE_CUTOFF.HEIGHT;

		 header << endl;
	 }
		
	if(ATOM_PAIRS[0].CUBIC_SCALE != 1.0)
		header << endl << "PAIR CHEBYSHEV CUBIC SCALING: " << ATOM_PAIRS[0].CUBIC_SCALE << endl;
	
	
	// Print out special cutoffs for 3-body interaction
	
	int FOUND_SPECIAL = 0;

	for(int i=0; i<PAIR_TRIPLETS.size(); i++)
	{
		if(PAIR_TRIPLETS[i].S_MINIM_3B.X >= 0)
			FOUND_SPECIAL++;
	}
	
	if(FOUND_SPECIAL>0)
	{
		header << endl << "SPECIAL 3B S_MINIM: SPECIFIC " << FOUND_SPECIAL << endl;
		
		for(int i=0; i<PAIR_TRIPLETS.size(); i++)
			if(PAIR_TRIPLETS[i].S_MINIM_3B.X >= 0)
				header << i << " " << TRIAD_MAP_REVERSE[i] << " " 
					<< PAIR_TRIPLETS[i].ATMPAIR1 << " " 
					<< PAIR_TRIPLETS[i].ATMPAIR2 << " " 
					<< PAIR_TRIPLETS[i].ATMPAIR3 << " " 
					<< fixed << setprecision(5) 
		            << PAIR_TRIPLETS[i].S_MINIM_3B.X << " "
				 	<< PAIR_TRIPLETS[i].S_MINIM_3B.Y << " "
					<< PAIR_TRIPLETS[i].S_MINIM_3B.Z << endl;						
	}
		
	FOUND_SPECIAL = 0;
	
	for(int i=0; i<PAIR_TRIPLETS.size(); i++)
		if(PAIR_TRIPLETS[i].S_MAXIM_3B.X >= 0)
			FOUND_SPECIAL++;
	
	if(FOUND_SPECIAL>0)
	{
		header << endl << "SPECIAL 3B S_MAXIM: SPECIFIC " << FOUND_SPECIAL << endl;
		
		for(int i=0; i<PAIR_TRIPLETS.size(); i++)
			if(PAIR_TRIPLETS[i].S_MAXIM_3B.X >= 0)
				header << i << " " << TRIAD_MAP_REVERSE[i] << " " 
					<< PAIR_TRIPLETS[i].ATMPAIR1 << " " 
					<< PAIR_TRIPLETS[i].ATMPAIR2 << " " 
					<< PAIR_TRIPLETS[i].ATMPAIR3 << " " 
					<< fixed << setprecision(5) 
		            << PAIR_TRIPLETS[i].S_MAXIM_3B.X << " "
				 	<< PAIR_TRIPLETS[i].S_MAXIM_3B.Y << " "
					<< PAIR_TRIPLETS[i].S_MAXIM_3B.Z << endl;						
	}	

	// Print out special cutoffs for 4-body interactions
	
	FOUND_SPECIAL = 0;

	for(int i=0; i<PAIR_QUADRUPLETS.size(); i++)
	{
		if(PAIR_QUADRUPLETS[i].S_MINIM[0] >= 0)
			FOUND_SPECIAL++;
	}
	
	if(FOUND_SPECIAL>0)
	{
		header << endl << "SPECIAL 4B S_MINIM: SPECIFIC " << FOUND_SPECIAL << endl;
		
		for(int i=0; i<PAIR_QUADRUPLETS.size(); i++)
		  PAIR_QUADRUPLETS[i].print_special(header, QUAD_MAP_REVERSE[i], "S_MINIM") ;
	}
		
	FOUND_SPECIAL = 0;
	
	for(int i=0; i<PAIR_QUADRUPLETS.size(); i++)
		if(PAIR_QUADRUPLETS[i].S_MAXIM[0] >= 0)
			FOUND_SPECIAL++;
	
	if(FOUND_SPECIAL>0)
	{
		header << endl << "SPECIAL 4B S_MAXIM: SPECIFIC " << FOUND_SPECIAL << endl;
		
		for(int i=0; i<PAIR_QUADRUPLETS.size(); i++)
		  PAIR_QUADRUPLETS[i].print_special(header, QUAD_MAP_REVERSE[i], "S_MAXIM" ) ;

	}	
	
	 
	if(!CONTROLS.USE_3B_CHEBY)
		header << endl << "ATOM PAIR TRIPLETS: " << 0 << endl;
	else
	{
		header << endl << "ATOM PAIR TRIPLETS: " << PAIR_TRIPLETS.size() << endl << endl;	

		for(int i=0;i<PAIR_TRIPLETS.size(); i++)
		{
			header << "" << PAIR_TRIPLETS[i].TRIPINDX << "  " << PAIR_TRIPLETS[i].ATMPAIR1 << " " << PAIR_TRIPLETS[i].ATMPAIR2 << " " << PAIR_TRIPLETS[i].ATMPAIR3 << ": ";
			header << PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS << " parameters, " << PAIR_TRIPLETS[i].N_ALLOWED_POWERS << " total parameters "<< endl;	
			header << "     index  |  powers  |  equiv index  |  param index  " << endl;
			header << "   ----------------------------------------------------" << endl;	

			for(int j=0; j<PAIR_TRIPLETS[i].ALLOWED_POWERS.size(); j++)
			{
				header << "      " << setw(6) << fixed << left << j << " ";
				header << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].X  << " ";
				header << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Y  << " ";
				header << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Z  << " ";
				header << "       " << setw(8) << PAIR_TRIPLETS[i].EQUIV_INDICES[j] << " ";
				header << "       " << setw(8) << PAIR_TRIPLETS[i].PARAM_INDICES[j] << endl; 
	
			}

			header << endl;
		}	 
		
	}
	
	if(!CONTROLS.USE_4B_CHEBY)
		header << "ATOM PAIR QUADRUPLETS: " << 0 << endl << endl;
	else
	{
		header << "ATOM PAIR QUADRUPLETS: " << PAIR_QUADRUPLETS.size() << endl << endl;	
		
		
		for(int i=0;i<PAIR_QUADRUPLETS.size(); i++)
		{
		  PAIR_QUADRUPLETS[i].print_header(header) ;
		}

	}

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

	if(TRIAD_MAP.size() > 0)
	{
		MAPFILE << endl;
		
		MAPFILE << "TRIPMAPS: " << TRIAD_MAP.size() << endl;
		
		for(map<string,int>::iterator i=TRIAD_MAP.begin(); i!=TRIAD_MAP.end(); i++)
			MAPFILE << i->second << " " << i->first << endl;
	}
	
	if(QUAD_MAP.size() > 0)
	{
		MAPFILE << endl;
		
		MAPFILE << "QUADMAPS: " << QUAD_MAP.size() << endl;
		
		for(map<string,int>::iterator i=QUAD_MAP.begin(); i!=QUAD_MAP.end(); i++)
			MAPFILE << i->second << " " << i->first << endl;
	}
	
	MAPFILE.close();

	//////////////////////////////////////////////////
	//
	// Print out bonding statistics.
	//
	//////////////////////////////////////////////////	  

	#if VERBOSITY == 1
	print_bond_stats(ATOM_PAIRS, PAIR_TRIPLETS, PAIR_QUADRUPLETS, CONTROLS.USE_3B_CHEBY, CONTROLS.USE_4B_CHEBY) ;
	#endif
	  
return 0;		  
}




	//////////////////////////////////////////////////
	//
    // Function definitions
	//
	//////////////////////////////////////////////////



// Read program input from the file "splines_ls.in".
static void read_lsq_input(JOB_CONTROL & CONTROLS, vector<PAIRS> & ATOM_PAIRS, vector<TRIPLETS> & PAIR_TRIPLETS, vector<QUADRUPLETS> & PAIR_QUADRUPLETS, map<string,int> & PAIR_MAP, map<int,string> & PAIR_MAP_REVERSE, map<string,int> & TRIAD_MAP, map<int,string> & TRIAD_MAP_REVERSE, map<string,int> & QUAD_MAP, map<int,string> & QUAD_MAP_REVERSE, vector<CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, NEIGHBORS & NEIGHBOR_LIST)
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
	
	vector<string> EXCLUDE_3B;
	
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
			if ( RANK == 0 ) 	cout << "	# N_LAYERS #: " << CONTROLS.N_LAYERS << endl;
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
				cin >> LINE;
				cin.ignore();
			}
			else
			{
				if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
					CONTROLS.FIT_STRESS = true;	
				else if (LINE=="all"  || LINE=="All"  || LINE=="ALL"  || LINE == "A" || LINE == "a")
					CONTROLS.FIT_STRESS_ALL = true;	
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
				cin >> LINE;
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
/*			
			if(CONTROLS.USE_3B_CHEBY && CONTROLS.N_LAYERS > 1)
			{
				cout << "ERROR: Use of layers is not supported with 3-body Chebyshev potentials." << endl;
				cout << "       Set # N_LAYERS # to 1." << endl;
				exit(0); 
			}
*/			
			if(CONTROLS.USE_3B_CHEBY && CONTROLS.FIT_POVER)
			{
				cout << "ERROR: Overbonding is not compatible with 3-body Chebyshev potentials." << endl;
				cout << "       Set # FITPOVR # false." << endl;
				exit(0);				
			}

			
		}

		else if( (TEMP_TYPE == "CHEBYSHEV") && (LINE.find("# CHBTYPE #") != string::npos))
		{
			cin >> LINE; cin.ignore();
			CONTROLS.CHEBY_TYPE = LINE;

			#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	# CHBTYPE #: " << CONTROLS.CHEBY_TYPE << endl;	
			#endif
		}		
		
		/////////////////////////////////////////////////////////////////////
		// Read the topology part. For now, ignoring index and atom types.. 
		// Assuming given as OO, OH, HH, as code expects... 
		// will need to be fixed later.
		/////////////////////////////////////////////////////////////////////
		
		else if(LINE.find("EXCLUDE 3B INTERACTION:")!= string::npos)
		{
			STREAM_PARSER.str(LINE);
			STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_INT;

			for(int i=0; i<TEMP_INT; i++)
			{
				cin >> LINE;
				EXCLUDE_3B.push_back(LINE);
			}
			
			// Sort the vector into ascending order
			
			sort (EXCLUDE_3B.begin(), EXCLUDE_3B.end());
			
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();	
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
					ATOM_PAIRS[i].CHEBY_RANGE_LOW  = TMP_CHEBY_RANGE_LOW;
					ATOM_PAIRS[i].CHEBY_RANGE_HIGH = TMP_CHEBY_RANGE_HIGH;
					ATOM_PAIRS[i].FORCE_CUTOFF.TYPE = FCUT_TYPE::CUBIC;
				}
			}	
			
			// Set up triplets
			
			NTRIP = factorial(CONTROLS.NATMTYP+3-1)/factorial(3)/factorial(CONTROLS.NATMTYP-1);
	
			// Account for excluded types:
			
			PAIR_TRIPLETS.resize(NTRIP);
			
			for (int i=0; i<NTRIP; i++)
			{	
				PAIR_TRIPLETS[i].S_MINIM_3B.X = -1;
				PAIR_TRIPLETS[i].S_MINIM_3B.Y = -1;
				PAIR_TRIPLETS[i].S_MINIM_3B.Z = -1;
				
				PAIR_TRIPLETS[i].S_MAXIM_3B.X = -1;
				PAIR_TRIPLETS[i].S_MAXIM_3B.Y = -1;
				PAIR_TRIPLETS[i].S_MAXIM_3B.Z = -1;
				
				PAIR_TRIPLETS[i].FORCE_CUTOFF.TYPE = FCUT_TYPE::CUBIC;
			}	

			// Set up quadruplets
			
			NQUAD = factorial(CONTROLS.NATMTYP+4-1)/factorial(4)/factorial(CONTROLS.NATMTYP-1);
			
			PAIR_QUADRUPLETS.resize(NQUAD);		
			
			for (int i=0; i<NQUAD; i++)
			{	
			  PAIR_QUADRUPLETS[i].init() ;
			}
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
				
				cin >> LINE >> ATOM_PAIRS[i].ATM1TYP >> LINE;
				
				//TMP_ATOMTYPEIDX = i;
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
					ATOM_PAIRS[TEMP_INT].ATM2TYP = ATOM_PAIRS[j].ATM1TYP;
					
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
			
			cout << "Made the following pairs: " << endl;
			for(map<string,int>::iterator i=PAIR_MAP.begin(); i!=PAIR_MAP.end(); i++)
				cout << i->second << " " << i->first << endl;			
			
				
			
			if(CONTROLS.USE_3B_CHEBY)
			{
				// Generate unique triplets
			
				TRIPLETS TRIP_ATOMS;		// Each item is an atom type
				TEMP_INT = 0;				// Will hold pair triplet index
			
				for(int i=0; i<CONTROLS.NATMTYP; i++)
				{
					for(int j=i; j<CONTROLS.NATMTYP; j++)
					{
						for(int k=j; k<CONTROLS.NATMTYP; k++)
						{
							// Get the three atoms that will define the 3-body interaction
						
							TRIP_ATOMS.ATMPAIR1 = ATOM_PAIRS[i].ATM1TYP;	// The first N_ATMTYP pairs are of type AA, BB, CC ... N_ATMTYP
							TRIP_ATOMS.ATMPAIR2 = ATOM_PAIRS[j].ATM1TYP;
							TRIP_ATOMS.ATMPAIR3 = ATOM_PAIRS[k].ATM1TYP;

							// Construct the triplet atom pairs from those atoms
						
							PAIR_TRIPLETS[TEMP_INT].TRIPINDX = TEMP_INT;
						
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR1 = TRIP_ATOMS.ATMPAIR1;	// ij
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR2 = TRIP_ATOMS.ATMPAIR1;	// ik
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR3 = TRIP_ATOMS.ATMPAIR2;	// jk

							PAIR_TRIPLETS[TEMP_INT].ATMPAIR1.append(TRIP_ATOMS.ATMPAIR2);	// ij
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR2.append(TRIP_ATOMS.ATMPAIR3);	// ik
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR3.append(TRIP_ATOMS.ATMPAIR3);	// jk		

							// Now save the "proper" (ordered) name of the pair			

							PAIR_TRIPLETS[TEMP_INT].ATMPAIR1 = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TEMP_INT].ATMPAIR1] ].PRPR_NM;
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR2 = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TEMP_INT].ATMPAIR2] ].PRPR_NM;
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR3 = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TEMP_INT].ATMPAIR3] ].PRPR_NM;
							
//cout << "ADDED TRIPLET: " << TEMP_INT << " " << PAIR_TRIPLETS[TEMP_INT].ATMPAIR1 << " " << PAIR_TRIPLETS[TEMP_INT].ATMPAIR2 << " " << PAIR_TRIPLETS[TEMP_INT].ATMPAIR3 << endl;	
if(RANK==0)
{							
	cout << "Made the following triplets: ";
	cout <<TEMP_INT << " " << PAIR_TRIPLETS[TEMP_INT].ATMPAIR1 << " " << PAIR_TRIPLETS[TEMP_INT].ATMPAIR2 << " " << PAIR_TRIPLETS[TEMP_INT].ATMPAIR3 << endl;		
}
	

							TEMP_INT++;						
						}
					}
				}
				
							
			
				// Figure out the allowed pair triplet powers. Here are some rules and considerations:
				//
				// 1. Powers start from zero. So, if order is specified to be 2, polynomial powers range
				//    from 0 to n-1, NOT 1 to n!
				//
				// 2. At least two pairs must have non-zero powers for the interaction to truly correspond
				//    to 3-body interactions
				//
				// 3. Non-uniqueness upon power sorting must be taken into consideration. For example, for
				//    the pair (OO, OH, OH), powers of (1,1,0) are identical to (1,0,1)
				// 
				// NOTE: We need to also take into consideration the corresponding parameter multiplicities.
			 
				XYZ_INT UNSORTED_POWERS;
				XYZ_INT SORTED_POWERS;
			
				vector<XYZ_INT> STORED_SORTED_POWERS;//(CONTROLS.CHEBY_3B_ORDER*CONTROLS.CHEBY_3B_ORDER*CONTROLS.CHEBY_3B_ORDER);		// Make this the max possible size... it will be destroyed later anyway.
			
				int TOP, BOT;
			
				bool STORED = false;
				int  STORED_IDX;
				int  RUNNING_IDX = 0;
			
				for(int i=0; i<NTRIP; i++)
				{
					vector<int> STORED_SORTED_POWERS_EQVS;
					
					for(int pair1_pow=0; pair1_pow<CONTROLS.CHEBY_3B_ORDER; pair1_pow++)
					{
						for(int pair2_pow=0; pair2_pow<CONTROLS.CHEBY_3B_ORDER; pair2_pow++)
						{
							for(int pair3_pow=0; pair3_pow<CONTROLS.CHEBY_3B_ORDER; pair3_pow++)
							{
								// Check number 1: Are at least two powers greater than 0?
								if( (pair1_pow>0 && pair2_pow>0) || (pair1_pow>0 && pair3_pow>0) || (pair2_pow>0 && pair3_pow>0) )
								{								
									UNSORTED_POWERS.X = pair1_pow;
									UNSORTED_POWERS.Y = pair2_pow;
									UNSORTED_POWERS.Z = pair3_pow;		
								
									// Store all triplet powers that meet the above criteria.
								
									PAIR_TRIPLETS[i].ALLOWED_POWERS.push_back(UNSORTED_POWERS);
								
									// Now, we need to figure out which of these sets of powers is truly unique. This will depend on a few things.
	
									// Case 1: each pair in the current pair triplet is unique ... multiplicities will be one for each item
								
									if( (PAIR_TRIPLETS[i].ATMPAIR1 != PAIR_TRIPLETS[i].ATMPAIR2) 
									 && (PAIR_TRIPLETS[i].ATMPAIR1 != PAIR_TRIPLETS[i].ATMPAIR3) 
									 && (PAIR_TRIPLETS[i].ATMPAIR2 != PAIR_TRIPLETS[i].ATMPAIR3) )
										PAIR_TRIPLETS[i].EQUIV_INDICES.push_back(PAIR_TRIPLETS[i].ALLOWED_POWERS.size()-1);

								
									// Case 2: each pair in the current pair triplet is identical... multiplicity should be 3 when all three powers are not the same
								 
									else if( (PAIR_TRIPLETS[i].ATMPAIR1 == PAIR_TRIPLETS[i].ATMPAIR2) && (PAIR_TRIPLETS[i].ATMPAIR1 == PAIR_TRIPLETS[i].ATMPAIR3))
									{
										// Sort the powers in ascending order

										if (pair1_pow>pair2_pow)
										{
											TOP = pair1_pow;
											BOT = pair2_pow;
										}
										else
										{
											TOP = pair2_pow;
											BOT = pair1_pow;
										}
										if(pair3_pow>TOP)
										{
											SORTED_POWERS.X = pair3_pow;
											SORTED_POWERS.Y = TOP;
											SORTED_POWERS.Z = BOT;										
										}
										else if(pair3_pow>BOT)
										{
											SORTED_POWERS.X = TOP;
											SORTED_POWERS.Y = pair3_pow;
											SORTED_POWERS.Z = BOT;											
										}
										else
										{
											SORTED_POWERS.X = TOP;
											SORTED_POWERS.Y = BOT;
											SORTED_POWERS.Z = pair3_pow;											
										}

										// Check if sorted powers have already been saved
									
										STORED = false;
									
										for(int j=0; j<STORED_SORTED_POWERS.size(); j++)
										{
											if( (STORED_SORTED_POWERS[j].X == SORTED_POWERS.X) &&  (STORED_SORTED_POWERS[j].Y == SORTED_POWERS.Y) &&  (STORED_SORTED_POWERS[j].Z == SORTED_POWERS.Z))
											{
												STORED = true;
												STORED_IDX = j;
												break;
											}										
										}
									
										// Save them, if they have not been already
									
										if(!STORED)
										{
											STORED_SORTED_POWERS           .push_back(SORTED_POWERS);
											STORED_SORTED_POWERS_EQVS      .push_back(PAIR_TRIPLETS[i].ALLOWED_POWERS.size()-1);
											PAIR_TRIPLETS[i].EQUIV_INDICES.push_back(PAIR_TRIPLETS[i].ALLOWED_POWERS.size()-1);	// The current power	
																			
										}
										else
										{
											PAIR_TRIPLETS[i].EQUIV_INDICES.push_back(STORED_SORTED_POWERS_EQVS[STORED_IDX]);											
										}
										STORED = false;

									} 
								
									// Case 3: two pairs in the current pair triplet are identical... multiplicities will be 2 whenever power differs on the two identical pairs
								
									else
									{
										// Sort the powers as if the pair types were arranged A != B == C

										// Case 3a.: A == B != C
									
										if(PAIR_TRIPLETS[i].ATMPAIR1 == PAIR_TRIPLETS[i].ATMPAIR2)
										{
											SORTED_POWERS.X = pair3_pow;
											SORTED_POWERS.Y = pair1_pow;
											SORTED_POWERS.Z = pair2_pow;										
										}
									
										// Case 3b.: A != B == C
									
										if(PAIR_TRIPLETS[i].ATMPAIR2 == PAIR_TRIPLETS[i].ATMPAIR3)
										{
											SORTED_POWERS.X = pair1_pow;
											SORTED_POWERS.Y = pair2_pow;
											SORTED_POWERS.Z = pair3_pow;										
										}			
									
										// Case 3c.: A == C != B
									
										if(PAIR_TRIPLETS[i].ATMPAIR1 == PAIR_TRIPLETS[i].ATMPAIR3)
										{
											SORTED_POWERS.X = pair2_pow;
											SORTED_POWERS.Y = pair1_pow;
											SORTED_POWERS.Z = pair3_pow;										
										}		
									
										// Now arrange the last two powers in ascending order Z>Y
									
										if(SORTED_POWERS.Y > SORTED_POWERS.Z)
										{
											BOT = SORTED_POWERS.Z;
											SORTED_POWERS.Z = SORTED_POWERS.Y;
											SORTED_POWERS.Y = BOT;
										}
									
										// Finally, check whether these sorted powers have already been saved
									
										STORED = false;
									
										for(int j=0; j<STORED_SORTED_POWERS.size(); j++)
										{
											if( (STORED_SORTED_POWERS[j].X == SORTED_POWERS.X) &&  (STORED_SORTED_POWERS[j].Y == SORTED_POWERS.Y) &&  (STORED_SORTED_POWERS[j].Z == SORTED_POWERS.Z))
											{
												STORED = true;
												STORED_IDX = j;
												break;
											}
										}
									
										// Save them, if they have not been already
									
										if(!STORED)
										{
											STORED_SORTED_POWERS           .push_back(SORTED_POWERS);
											STORED_SORTED_POWERS_EQVS      .push_back(PAIR_TRIPLETS[i].ALLOWED_POWERS.size()-1);
											PAIR_TRIPLETS[i].EQUIV_INDICES.push_back(PAIR_TRIPLETS[i].ALLOWED_POWERS.size()-1);	
										}
										else
										{										
											PAIR_TRIPLETS[i].EQUIV_INDICES.push_back(STORED_SORTED_POWERS_EQVS[STORED_IDX]);
										}									
										STORED = false;
					
									}
								}
							}
						}
					}
				
					STORED_SORTED_POWERS.clear();
					vector<XYZ_INT>().swap(STORED_SORTED_POWERS);
				
					// Now all that's left to do is set the force field index for each set of powers
				 
					PAIR_TRIPLETS[i].PARAM_INDICES.resize(PAIR_TRIPLETS[i].EQUIV_INDICES.size());
				
					PAIR_TRIPLETS[i].PARAM_INDICES[0] = 0;
					
					bool FOUND_EQV;
					int  USE_SET = 0;
					int  MAX_SET = 0;
					
					for(int set1=1; set1<PAIR_TRIPLETS[i].EQUIV_INDICES.size(); set1++)
					{
						FOUND_EQV = false;
						
						for(int set2=0; set2<set1; set2++)
						{
							if(PAIR_TRIPLETS[i].EQUIV_INDICES[set1] == PAIR_TRIPLETS[i].EQUIV_INDICES[set2])
							{
								FOUND_EQV = true;
								USE_SET   = set2;
								break;
							}
						}
						
						if(FOUND_EQV)
							PAIR_TRIPLETS[i].PARAM_INDICES[set1] = PAIR_TRIPLETS[i].PARAM_INDICES[USE_SET];
						else
						{
							MAX_SET++;
							PAIR_TRIPLETS[i].PARAM_INDICES[set1] = MAX_SET;	
						}					

					}

					PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS = PAIR_TRIPLETS[i].PARAM_INDICES[PAIR_TRIPLETS[i].PARAM_INDICES.size()-1]+1;
					PAIR_TRIPLETS[i].N_ALLOWED_POWERS = PAIR_TRIPLETS[i].PARAM_INDICES.size();
				}

				// Set up triplet maps... Account for cases where triplet type is meaningless by setting mapped index to -1
					
				bool REAL_TRIPLET = false;
			
				string TEMP_STR_A, TEMP_STR_B, TEMP_STR_C;
				
			
/*			
				cout << "SANITY CHECK: These are your triplets:" << endl;
				for(int m=0; m<NTRIP; m++)
				{
					cout << m << " " << PAIR_TRIPLETS[m].ATMPAIR1 << " " << PAIR_TRIPLETS[m].ATMPAIR2 << " " << PAIR_TRIPLETS[m].ATMPAIR3 << endl;
				}
*/			
				
				for(int i=0; i<NPAIR; i++)
				{
					for(int j=0; j<NPAIR; j++)
					{
						for(int k=0; k<NPAIR; k++)
						{	
							TEMP_STR_A = ATOM_PAIRS[i].ATM1TYP;
							TEMP_STR_A.append(ATOM_PAIRS[i].ATM2TYP);
							
							TEMP_STR_B = ATOM_PAIRS[j].ATM1TYP;
							TEMP_STR_B.append(ATOM_PAIRS[j].ATM2TYP);	
						
							TEMP_STR_C = ATOM_PAIRS[k].ATM1TYP;
							TEMP_STR_C.append(ATOM_PAIRS[k].ATM2TYP);	
						
							REAL_TRIPLET = false;	
						
							for(int m=0; m<NTRIP; m++)
							{
								if ((TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR3) )
									
								{
									TEMP_STR = TEMP_STR_A;
									TEMP_STR.append(TEMP_STR_B);	
									TEMP_STR.append(TEMP_STR_C);

									TRIAD_MAP        .insert(make_pair(TEMP_STR,m));	
									TRIAD_MAP_REVERSE.insert(make_pair(m,TEMP_STR));
									REAL_TRIPLET = true;									
								
								}
								if(!REAL_TRIPLET)	// Then this is not an allowed triplet type!
									TRIAD_MAP        .insert(make_pair(TEMP_STR,-1*m-1));	
							}
						}
					}
				}
				
				

				
				
				// Now that we've got the maps and the 3b ff structures created, we can go back an remove triplet types
				// that the user requested to exclude


				TEMP_INT = NTRIP;
				vector<int>EXCL_IDX;
				bool FOUND = false;
				map<string, int>::iterator it2, it2a,itrem;
				map<string, int>::iterator it, ita, itb;
				ita = TRIAD_MAP.begin();
				itb = TRIAD_MAP.end();
				advance(itb,-1);
				int TARGET;
				
/*				
				//	Sanity check	
				cout << "YOUR OLD MAPS: " << endl;	
				for(it = TRIAD_MAP.begin(); it != TRIAD_MAP.end(); it++)
					cout <<"		" << it->first << " : " << it->second << endl;
*/				
				
				// Start by making all removed type indicies negative
				// ALSO, THIS IS WHERE WE POP OFF ELEMENTS OF OUR PAIR TRIPLET VECTOR
				 
				vector<int> POPPED; 

				for (unsigned j = EXCLUDE_3B.size(); j-- > 0; ) // Since we're popping off by index, iterate over vector (ascending sorted) in reverse
					PAIR_TRIPLETS.erase (PAIR_TRIPLETS.begin() + TRIAD_MAP[EXCLUDE_3B[j]]);
				
				for(int j=0; j<EXCLUDE_3B.size(); j++)
				{
					EXCL_IDX.push_back(TRIAD_MAP[EXCLUDE_3B[j]]);

					NTRIP--;
					
					for(it = TRIAD_MAP.begin(); it != TRIAD_MAP.end(); it++)				
					{
						if(it->second == EXCL_IDX[j])	// Then we need to exclude it!
						{		
							POPPED.push_back(it->second);
							it->second = -1*it->second - 1;	
						}
					}				
				}	
				
				// Sort the popoff list in ascending order
				
				sort (POPPED.begin(), POPPED.end());
				
/*				
				//	Sanity check	
				cout << "YOUR ~~ MAPS: " << endl;	
				for(it = TRIAD_MAP.begin(); it != TRIAD_MAP.end(); it++)
					cout <<"		" << it->first << " : " << it->second << endl;
				
*/
								
				for(int i=0; i<POPPED.size(); i++)
				{
					for(it = TRIAD_MAP.begin(); it != TRIAD_MAP.end(); it++)
					{
						if(it->second>POPPED[i])
							it->second -= 1;
					}
				}
/*				
				//	Sanity check	
				cout << "YOUR NEW MAPS: " << endl;	
				for(it = TRIAD_MAP.begin(); it != TRIAD_MAP.end(); it++)
					cout <<"		" << it->first << " : " << it->second << endl;
*/				
				
				// Finally, rebuild the reverse maps
				
				TRIAD_MAP_REVERSE.clear();

				for(it = TRIAD_MAP.begin(); it != TRIAD_MAP.end(); it++)
					TRIAD_MAP_REVERSE.insert(make_pair(it->second,it->first));

/*				
				//	Sanity check	
							
				cout << "YOUR NEW REV MAPS: " << endl;

				map<int, string>::iterator itc;
				
				for(itc = TRIAD_MAP_REVERSE.begin(); itc != TRIAD_MAP_REVERSE.end(); itc++)
					cout <<"		" << itc->first << " : " << itc->second << endl;
*/				
/*				
				// Sanity check
				
				cout << "Triplet types (force field): " << endl;
				for(int i=0;i<NTRIP; i++)
					cout << "		" << PAIR_TRIPLETS[i].TRIPINDX << "  " << PAIR_TRIPLETS[i].ATMPAIR1 << " " << PAIR_TRIPLETS[i].ATMPAIR2 << " " << PAIR_TRIPLETS[i].ATMPAIR3 << endl;
*/
			}
								
			if(CONTROLS.USE_4B_CHEBY)
			{
				// First, extract the atom types

				vector<string>ATOM_CHEMS;

				for(int p=0; p<ATOM_PAIRS.size(); p++)
				{
					string TMP_CHEM = ATOM_PAIRS[p].ATM1TYP;
					bool   IN_LIST  = false;

					for(int a=0; a<ATOM_CHEMS.size(); a++)
					{
						if(ATOM_CHEMS[a] == TMP_CHEM)
							IN_LIST = true;
					}
	
					if(!IN_LIST)
						ATOM_CHEMS.push_back(TMP_CHEM);
				}
				
				
				//////////////////////////////////////////////////////////////////////
				// Generate unique quadruplets and thier corresponding sets of powers
				//////////////////////////////////////////////////////////////////////

				build_quad_pairs(PAIR_QUADRUPLETS, CONTROLS.NATMTYP, ATOM_CHEMS, ATOM_PAIRS, PAIR_MAP) ;
			
				for(int i=0; i<NQUAD; i++)
				{
				  PAIR_QUADRUPLETS[i].build(CONTROLS.CHEBY_4B_ORDER) ;
				}

				///////////////////////////////////////////////////////////////////////////////////////////////////////////
				// Set up quadruplet maps... Account for cases where quadruplet type is meaningless by setting mapped index to -1
				///////////////////////////////////////////////////////////////////////////////////////////////////////////					


				build_quad_maps(QUAD_MAP, QUAD_MAP_REVERSE, ATOM_PAIRS, NPAIR, PAIR_QUADRUPLETS, NQUAD) ;

				//////////////////////////////////////////////////////////////////////
				// Since there are so many pairs in an atom quadruplet, we'll use "fast" maps here.
				// Generate "fast" maps for 4-body interactions. Assumes never more than MAX_ATOM_TYPES
				// atom types in a simulation
				//////////////////////////////////////////////////////////////////////
				build_fast_quad_maps(QUAD_MAP, INT_QUAD_MAP, INT_QUAD_MAP_REVERSE, ATOM_CHEMS) ;
			}  // CONTROLS.USE_4B_CHEBY

#if VERBOSITY == 1						
			if ( RANK == 0 ) 
			{
				cout << endl;
				cout << "	The following unique pair types have been identified:" << endl;
				for(int i=0;i<NPAIR; i++)
						cout << "		" << ATOM_PAIRS[i].PAIRIDX << "  " << ATOM_PAIRS[i].ATM1TYP << " " << ATOM_PAIRS[i].ATM2TYP << endl;
			}
#endif
					
			if(CONTROLS.USE_3B_CHEBY)
			{
				#if VERBOSITY == 1	
				if ( RANK == 0 ) 
				{
					cout << "	The following unique triplets of pair types and thier allowed pair polynomial powers have been identified:" << endl;
					cout << "	Note: The following types have been removed, if present: " << endl;
					
					for(int i=0;i<EXCLUDE_3B.size(); i++)
						cout << "		" << EXCLUDE_3B[i] << endl;
					cout << "	" << endl;
					for(int i=0;i<NTRIP; i++)
					{
						cout << "		" << PAIR_TRIPLETS[i].TRIPINDX << "  " << PAIR_TRIPLETS[i].ATMPAIR1 << " " << PAIR_TRIPLETS[i].ATMPAIR2 << " " << PAIR_TRIPLETS[i].ATMPAIR3 << ":";
						cout << " Number of unique sets of powers: " << PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS << " (" << PAIR_TRIPLETS[i].N_ALLOWED_POWERS << " total)..." << endl;	
						cout << "		     index  |  powers  |  equiv index  |  param index  " << endl;
						cout << "		   ----------------------------------------------------" << endl;					
					
						for(int j=0; j<PAIR_TRIPLETS[i].ALLOWED_POWERS.size(); j++)
						{
							//						cout << "		   " << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].X << " " << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Y << " " << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Z << endl;		
						 
										

							cout << "		      " << setw(6) << fixed << left << j << " ";
							cout << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].X  << " ";
							cout << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Y  << " ";
							cout << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Z  << " ";
							cout << "       " << setw(8) << PAIR_TRIPLETS[i].EQUIV_INDICES[j] << " ";
							cout << "       " << setw(8) << PAIR_TRIPLETS[i].PARAM_INDICES[j] << endl; 
						
						}

					}
				}
				#endif	
			}	

			if(CONTROLS.USE_4B_CHEBY)
			{
				#if VERBOSITY == 1	
				if ( RANK == 0 ) 
				{
					cout << "	The following unique six-lets of pair types and thier allowed pair polynomial powers have been identified:" << endl;

					for(int i=0;i<NQUAD; i++)
					{
					  PAIR_QUADRUPLETS[i].print() ;
					}
				}
				#endif		
				
			} // CONTROLS.USE_4B_CHEBY
			
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
					cin >> ATOM_PAIRS[TEMP_INT].NBINS.X;
					cin >> ATOM_PAIRS[TEMP_INT].NBINS.Y;
					cin >> ATOM_PAIRS[TEMP_INT].NBINS.Z;
				}
				else
				{
					ATOM_PAIRS[TEMP_INT].NBINS.X = 0;
					ATOM_PAIRS[TEMP_INT].NBINS.Y = 0;
					ATOM_PAIRS[TEMP_INT].NBINS.Z = 0;
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
					cout << "		" << ATOM_PAIRS[i].PRPR_NM  << ": " << ATOM_PAIRS[i].NBINS.X << " " << ATOM_PAIRS[i].NBINS.Y << " " << ATOM_PAIRS[i].NBINS.Z << endl;					
				cout << endl;
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
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();	
			
			STREAM_PARSER.str(LINE);
			
			STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR;
			
			if(TEMP_STR == "ALL" && NTRIP >0 )
			{				
				STREAM_PARSER >> PAIR_TRIPLETS[0].S_MAXIM_3B.X;
				
				if(PAIR_TRIPLETS[0].S_MAXIM_3B.X>NEIGHBOR_LIST.MAX_CUTOFF_3B)
					NEIGHBOR_LIST.MAX_CUTOFF_3B = PAIR_TRIPLETS[0].S_MAXIM_3B.X;
				
				for(int i=0; i<NTRIP; i++)
				{
					PAIR_TRIPLETS[i].S_MAXIM_3B.X = PAIR_TRIPLETS[0].S_MAXIM_3B.X;
					PAIR_TRIPLETS[i].S_MAXIM_3B.Y = PAIR_TRIPLETS[0].S_MAXIM_3B.X;
					PAIR_TRIPLETS[i].S_MAXIM_3B.Z = PAIR_TRIPLETS[0].S_MAXIM_3B.X;
				}
				
				#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "	Note: Setting all 3-body r_max values to " <<  PAIR_TRIPLETS[0].S_MAXIM_3B.X << endl;
				#endif				
			}
			else if(TEMP_STR == "SPECIFIC" && NTRIP >0 )
			{
				STREAM_PARSER >> TEMP_INT;
				STREAM_PARSER.str("");
				STREAM_PARSER.clear();

				
				#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "	Note: Setting specific 3-body r_max values: " << endl;
				#endif	

				double TMP_VAL;
				string  TMP_IJ,  TMP_IK,  TMP_JK;	
				string TARG_IJ, TARG_IK, TARG_JK;
					
				
				for(int i=0; i<TEMP_INT; i++)
				{
					getline(cin,LINE);
					STREAM_PARSER.str(LINE);
					STREAM_PARSER >> TEMP_STR;	// Which 3-body type is it?
					
					
					STREAM_PARSER >> TMP_IJ;	// What is the IJ?
					STREAM_PARSER >> TMP_IK;	// What is the IK?
					STREAM_PARSER >> TMP_JK;	// What is the JK?
					
					// Check that triplet pair types are correct
					
					try
					{
						TMP_IJ = ATOM_PAIRS[ PAIR_MAP[ TMP_IJ ] ].PRPR_NM;
					}
					catch(...)
					{
						cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
						cout << "		Triplet type:              " << TEMP_STR << endl;
						cout << "		First distance, pair type: " << TMP_IJ << endl;
					}
					try
					{
						TMP_IK = ATOM_PAIRS[ PAIR_MAP[ TMP_IK ] ].PRPR_NM;
					}
					catch(...)
					{
						cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
						cout << "		Triplet type:              " << TEMP_STR << endl;
						cout << "		First distance, pair type: " << TMP_IK << endl;
					}
					try
					{
						TMP_JK = ATOM_PAIRS[ PAIR_MAP[ TMP_JK ] ].PRPR_NM;
					}
					catch(...)
					{
						cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
						cout << "		Triplet type:              " << TEMP_STR << endl;
						cout << "		First distance, pair type: " << TMP_JK << endl;
					}
					
					TARG_IJ = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].ATMPAIR1 ] ].PRPR_NM;
					TARG_IK = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].ATMPAIR2 ] ].PRPR_NM;
					TARG_JK = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].ATMPAIR3 ] ].PRPR_NM;
					
					// Read the first inner cutoff

					STREAM_PARSER >> TMP_VAL;
					
					if      ( (TMP_IJ == TARG_IJ) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.X == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.X = TMP_VAL;
					
					else if ( (TMP_IJ == TARG_IK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Y == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Y = TMP_VAL;
					
					else if ( (TMP_IJ == TARG_JK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Z == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Z = TMP_VAL;


					// Read the second inner cutoff

					STREAM_PARSER >> TMP_VAL;
					
					if      ( (TMP_IK == TARG_IJ) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.X == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.X = TMP_VAL;
					
					else if ( (TMP_IK == TARG_IK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Y == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Y = TMP_VAL;
					
					else if ( (TMP_IK == TARG_JK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Z == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Z = TMP_VAL;

					
					// Read the third inner cutoff

					STREAM_PARSER >> TMP_VAL;
					
					if      ( (TMP_JK == TARG_IJ) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.X == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.X = TMP_VAL;
					
					else if ( (TMP_JK == TARG_IK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Y == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Y = TMP_VAL;
					
					else if ( (TMP_JK == TARG_JK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Z == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Z = TMP_VAL;
					
					
					if(PAIR_TRIPLETS[0].S_MAXIM_3B.X > NEIGHBOR_LIST.MAX_CUTOFF_3B)
						NEIGHBOR_LIST.MAX_CUTOFF_3B = PAIR_TRIPLETS[0].S_MAXIM_3B.X;
					if(PAIR_TRIPLETS[0].S_MAXIM_3B.Y > NEIGHBOR_LIST.MAX_CUTOFF_3B)
						NEIGHBOR_LIST.MAX_CUTOFF_3B = PAIR_TRIPLETS[0].S_MAXIM_3B.Y;
					if(PAIR_TRIPLETS[0].S_MAXIM_3B.Z > NEIGHBOR_LIST.MAX_CUTOFF_3B)
						NEIGHBOR_LIST.MAX_CUTOFF_3B = PAIR_TRIPLETS[0].S_MAXIM_3B.Z;
					
					
					PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B;
					STREAM_PARSER.str("");
					STREAM_PARSER.clear();
					
					#if VERBOSITY == 1
					if ( RANK == 0 ) 
					{
						cout << "		" << TEMP_STR << "(" <<  TARG_IJ << ", " << TARG_IK << ", " << TARG_JK << "): " 
							              << PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.X << ", "
									 	  << PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Y << ", "
										  << PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MAXIM_3B.Z << endl;
					}
					#endif	
				}
			}
		}
		
		else if(LINE.find("SPECIAL 3B S_MINIM:") != string::npos)
		{
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();	
			
			STREAM_PARSER.str(LINE);
			
			STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR;
			
			if(TEMP_STR == "ALL" && NTRIP >0 )
			{				
				STREAM_PARSER >> PAIR_TRIPLETS[0].S_MINIM_3B.X;
				
				for(int i=0; i<NTRIP; i++)
				{
					PAIR_TRIPLETS[i].S_MINIM_3B.X = PAIR_TRIPLETS[0].S_MINIM_3B.X;
					PAIR_TRIPLETS[i].S_MINIM_3B.Y = PAIR_TRIPLETS[0].S_MINIM_3B.X;
					PAIR_TRIPLETS[i].S_MINIM_3B.Z = PAIR_TRIPLETS[0].S_MINIM_3B.X;
				}
				
				#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "	Note: Setting all 3-body r_min values to " <<  PAIR_TRIPLETS[0].S_MINIM_3B.X << endl;
				#endif				
			}
			else if(TEMP_STR == "SPECIFIC" && NTRIP >0 )
			{
				STREAM_PARSER >> TEMP_INT;
				STREAM_PARSER.str("");
				STREAM_PARSER.clear();

				
				#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "	Note: Setting specific 3-body r_min values: " << endl;
				#endif	

				double TMP_VAL;
				string  TMP_IJ,  TMP_IK,  TMP_JK;	
				string TARG_IJ, TARG_IK, TARG_JK;
					
				
				for(int i=0; i<TEMP_INT; i++)
				{
					getline(cin,LINE);
					STREAM_PARSER.str(LINE);
					STREAM_PARSER >> TEMP_STR;	// Which 3-body type is it?
					
					STREAM_PARSER >> TMP_IJ;	// What is the IJ?
					STREAM_PARSER >> TMP_IK;	// What is the IK?
					STREAM_PARSER >> TMP_JK;	// What is the JK?
					
					// Check that triplet pair types are correct
					
					try
					{
						TMP_IJ = ATOM_PAIRS[ PAIR_MAP[ TMP_IJ ] ].PRPR_NM;
					}
					catch(...)
					{
						cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
						cout << "		Triplet type:              " << TEMP_STR << endl;
						cout << "		First distance, pair type: " << TMP_IJ << endl;
					}
					try
					{
						TMP_IK = ATOM_PAIRS[ PAIR_MAP[ TMP_IK ] ].PRPR_NM;
					}
					catch(...)
					{
						cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
						cout << "		Triplet type:              " << TEMP_STR << endl;
						cout << "		First distance, pair type: " << TMP_IK << endl;
					}
					try
					{
						TMP_JK = ATOM_PAIRS[ PAIR_MAP[ TMP_JK ] ].PRPR_NM;
					}
					catch(...)
					{
						cout << "ERROR: Unknown triplet pair for special inner cutoff." << endl;
						cout << "		Triplet type:              " << TEMP_STR << endl;
						cout << "		First distance, pair type: " << TMP_JK << endl;
					}
					
					TARG_IJ = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].ATMPAIR1 ] ].PRPR_NM;
					TARG_IK = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].ATMPAIR2 ] ].PRPR_NM;
					TARG_JK = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].ATMPAIR3 ] ].PRPR_NM;
					
					// Read the first inner cutoff

					STREAM_PARSER >> TMP_VAL;
					
					if      ( (TMP_IJ == TARG_IJ) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X = TMP_VAL;
					
					else if ( (TMP_IJ == TARG_IK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y = TMP_VAL;
					
					else if ( (TMP_IJ == TARG_JK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z = TMP_VAL;


					// Read the second inner cutoff

					STREAM_PARSER >> TMP_VAL;
					
					if      ( (TMP_IK == TARG_IJ) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X = TMP_VAL;
					
					else if ( (TMP_IK == TARG_IK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y = TMP_VAL;
					
					else if ( (TMP_IK == TARG_JK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z = TMP_VAL;

					
					// Read the third inner cutoff

					STREAM_PARSER >> TMP_VAL;
					
					if      ( (TMP_JK == TARG_IJ) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X = TMP_VAL;
					
					else if ( (TMP_JK == TARG_IK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y = TMP_VAL;
					
					else if ( (TMP_JK == TARG_JK) && (PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z == -1) )
						PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z = TMP_VAL;
					
					
					
					PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B;
					STREAM_PARSER.str("");
					STREAM_PARSER.clear();
					
					#if VERBOSITY == 1
					if ( RANK == 0 ) 
						cout << "		" << TEMP_STR << "(" <<  TARG_IJ << ", " << TARG_IK << ", " << TARG_JK << "): " 
							  << setw(10) << fixed << right << setprecision(4) << PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.X << ", "
							  << setw(10) << fixed << right << setprecision(4) << PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Y << ", "
							  << setw(10) << fixed << right << setprecision(4) << PAIR_TRIPLETS[TRIAD_MAP[TEMP_STR]].S_MINIM_3B.Z << endl;
					#endif	
				}
			}
		}
		
		else if(LINE.find("SPECIAL 4B S_MAXIM:") != string::npos)
		{
			STREAM_PARSER.str("");
			STREAM_PARSER.clear();	
			
			STREAM_PARSER.str(LINE);
			
			STREAM_PARSER >> TEMP_STR >> TEMP_STR >> TEMP_STR >> TEMP_STR;
			
			if(TEMP_STR == "ALL" && NQUAD >0 )
			{				
				STREAM_PARSER >> PAIR_QUADRUPLETS[0].S_MAXIM[0];
				
				if(PAIR_QUADRUPLETS[0].S_MAXIM[0] > NEIGHBOR_LIST.MAX_CUTOFF_4B)
					NEIGHBOR_LIST.MAX_CUTOFF_4B = PAIR_QUADRUPLETS[0].S_MAXIM[0];

				for(int i=0; i<NQUAD; i++)
				{
					for(int j=0; j<6; j++)
						PAIR_QUADRUPLETS[i].S_MAXIM[j] = PAIR_QUADRUPLETS[0].S_MAXIM[0];
				}
				
				#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "	Note: Setting all 4-body r_max values to " <<  PAIR_QUADRUPLETS[0].S_MAXIM[0] << endl;
				#endif				
			}
			else if(TEMP_STR == "SPECIFIC" && NQUAD >0 )
			{
				STREAM_PARSER >> TEMP_INT;
				STREAM_PARSER.str("");
				STREAM_PARSER.clear();
			
				
				#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "	Note: Setting specific 4-body r_max values: " << endl;
				#endif	

				double TMP_VAL;
				vector<string>TMP_PAIRTYPS(6);
				vector<double>TMP_CUTOFFS(6);
				vector<string>TRG_PAIRTYPS(6);
				
				for(int i=0; i<TEMP_INT; i++)
				{
					getline(cin,LINE);

					STREAM_PARSER.str(LINE);
					STREAM_PARSER >> TEMP_STR;	// Which 3-body type is it?

					// Figure out the pair types for the quadruplet being read, and figure out the expected order
					// of those pair types for the quadruplet force field type
					
					for(int j=0; j<6; j++)
					{
						STREAM_PARSER >> TMP_PAIRTYPS[j];	// What are the pair types (strings)
						
						// Check that sixlet pair types are correct
						
						try
						{
							TMP_PAIRTYPS[j] = ATOM_PAIRS[ PAIR_MAP[ TMP_PAIRTYPS[j] ] ].PRPR_NM;
						}
						catch(...)
						{
							cout << "ERROR: Unknown quad pair for special inner cutoff." << endl;
							cout << "		Quadruplet type:           " << TEMP_STR        << endl;
							cout << "		First distance, pair type: " << TMP_PAIRTYPS[j] << endl;
						}
						
						TRG_PAIRTYPS[j] = ATOM_PAIRS[ PAIR_MAP[ PAIR_QUADRUPLETS[QUAD_MAP[TEMP_STR]].ATOM_PAIRS[j] ] ].PRPR_NM;

					}
					
					// Read the cutoffs
					
					for(int j=0; j<6; j++)
						STREAM_PARSER >> TMP_CUTOFFS[j];
					
					// Store cutoffs
					
					for(int j=0; j<6; j++)	// Iterate over TMP pair types
					{
						for(int k=0; k<6; k++)	// Iterate over TRG pair types
							if( (TMP_PAIRTYPS[j] == TRG_PAIRTYPS[k]) && (PAIR_QUADRUPLETS[QUAD_MAP[TEMP_STR]].S_MAXIM[k] == -1))
								PAIR_QUADRUPLETS[QUAD_MAP[TEMP_STR]].S_MAXIM[k] = TMP_CUTOFFS[j];
						
						// Reset the 4b max cutoff value
						
						if(TMP_CUTOFFS[j] > NEIGHBOR_LIST.MAX_CUTOFF_4B)
							NEIGHBOR_LIST.MAX_CUTOFF_4B = TMP_CUTOFFS[j];
					}
					
					// Print new cutoffs
					
					#if VERBOSITY == 1
					if(RANK==0)
					{
						cout << "		" << TEMP_STR << "(";
						
						for (int j=0; j<5; j++)
							cout << TRG_PAIRTYPS[j] << ", ";
						
						cout << TRG_PAIRTYPS[5] << "): ";
							
						for (int j=0; j<5; j++)
							cout << PAIR_QUADRUPLETS[QUAD_MAP[TEMP_STR]].S_MAXIM[j] << ", ";
						
						cout << PAIR_QUADRUPLETS[QUAD_MAP[TEMP_STR]].S_MAXIM[5] << endl;	
					}
					#endif	
					
					STREAM_PARSER.str("");
					STREAM_PARSER.clear();
				}
			}

		}
		else if(LINE.find("SPECIAL 4B S_MINIM:") != string::npos)
		{
			cout << "ERROR: Functionality not programmed for 4-body yet." << endl;
			exit(0);
		}
		
		else if (LINE.find("# FCUTTYP #") != string::npos)
		{
			cin >> TEMP_TYPE;
			
			for(int i=0; i<ATOM_PAIRS.size(); i++)
				ATOM_PAIRS[i].FORCE_CUTOFF.set_type(TEMP_TYPE);
			
			for(int i=0; i<PAIR_TRIPLETS.size(); i++)
				PAIR_TRIPLETS[i].FORCE_CUTOFF.set_type(TEMP_TYPE);
			
			for(int i=0; i<PAIR_QUADRUPLETS.size(); i++)
				PAIR_QUADRUPLETS[i].FORCE_CUTOFF.set_type(TEMP_TYPE);
			
			#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	# FCUTTYP #: " << TEMP_TYPE << "	... for 3-body and 4-body Chebyshev interactions" << endl;	
			#endif
			
			if(TEMP_TYPE=="SIGMOID" || TEMP_TYPE=="CUBSIG" || TEMP_TYPE=="CUBESTRETCH" || TEMP_TYPE == "SIGFLT")
			{
				cin >> PAIR_TRIPLETS[0].FORCE_CUTOFF.STEEPNESS;
				cin >> PAIR_TRIPLETS[0].FORCE_CUTOFF.OFFSET;
				
				if(PAIR_TRIPLETS[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::SIGFLT)
					cin >> PAIR_TRIPLETS[0].FORCE_CUTOFF.HEIGHT;
				
				#if VERBOSITY == 1
				if ( RANK == 0 ) cout << "			With steepness, offset, and height of: " 
											 << fixed << setprecision(3) 
											 << PAIR_TRIPLETS[0].FORCE_CUTOFF.STEEPNESS 
											 << ",	" << PAIR_TRIPLETS[0].FORCE_CUTOFF.OFFSET 
											 << ", and  " << PAIR_TRIPLETS[0].FORCE_CUTOFF.HEIGHT << endl;	
				#endif					
			}
			cin.ignore();
			
			PAIR_TRIPLETS[0].FORCE_CUTOFF.BODIEDNESS = 3 ;

			for(int i=1; i<PAIR_TRIPLETS.size(); i++)	// Copy all class elements
			{
//				PAIR_TRIPLETS[i].FORCE_CUTOFF = PAIR_TRIPLETS[0].FORCE_CUTOFF;
				PAIR_TRIPLETS[i].FORCE_CUTOFF.BODIEDNESS = PAIR_TRIPLETS[0].FORCE_CUTOFF.BODIEDNESS;
			}
			
			if(CONTROLS.USE_4B_CHEBY)
			{
//				PAIR_QUADRUPLETS[0].FORCE_CUTOFF           =  PAIR_TRIPLETS[0].FORCE_CUTOFF;
				PAIR_QUADRUPLETS[0].FORCE_CUTOFF.STEEPNESS =  PAIR_TRIPLETS[0].FORCE_CUTOFF.STEEPNESS;
				PAIR_QUADRUPLETS[0].FORCE_CUTOFF.OFFSET    =  PAIR_TRIPLETS[0].FORCE_CUTOFF.OFFSET;
			
				if(PAIR_QUADRUPLETS[0].FORCE_CUTOFF.TYPE == FCUT_TYPE::SIGFLT)
					PAIR_QUADRUPLETS[0].FORCE_CUTOFF.HEIGHT = PAIR_TRIPLETS[0].FORCE_CUTOFF.HEIGHT;
			
				PAIR_QUADRUPLETS[0].FORCE_CUTOFF.BODIEDNESS = 4 ;

				for(int i=1; i<PAIR_QUADRUPLETS.size(); i++)	// Copy all class elements
				{
//					PAIR_QUADRUPLETS[i].FORCE_CUTOFF = PAIR_QUADRUPLETS[0].FORCE_CUTOFF;
					PAIR_QUADRUPLETS[i].FORCE_CUTOFF.BODIEDNESS = PAIR_QUADRUPLETS[0].FORCE_CUTOFF.BODIEDNESS;
				}
			}

		}	
		
		
	}	
	
	#if VERBOSITY == 1			
	if ( RANK == 0 ) cout << endl << "Note: Will use cubic scaling of: " << ATOM_PAIRS[0].CUBIC_SCALE << endl << endl;; // All types use same scaling
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

static void print_bond_stats(vector<PAIRS> &ATOM_PAIRS, vector<TRIPLETS> &PAIR_TRIPLETS, vector<QUADRUPLETS> &PAIR_QUADRUPLETS, bool use_3b_cheby, bool use_4b_cheby)
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
			if ( RANK == 0 ) cout << "	Minimum distances between atoms triplet pairs: (Angstr.)" << endl;
			
			for (int k=0; k<PAIR_TRIPLETS.size(); k++) 
			{
				XYZ sum = {0.0, 0.0, 0.0} ;
#ifdef USE_MPI
				MPI_Reduce(&(PAIR_TRIPLETS[k].MIN_FOUND.X), &(sum.X), 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
				MPI_Reduce(&(PAIR_TRIPLETS[k].MIN_FOUND.Y), &(sum.Y), 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
				MPI_Reduce(&(PAIR_TRIPLETS[k].MIN_FOUND.Z), &(sum.Z), 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
#else				
				sum = PAIR_TRIPLETS[k].MIN_FOUND ;
#endif
				if ( RANK == 0 ) cout << "		" << k << "	" 					
											 << PAIR_TRIPLETS[k].ATMPAIR1    << " " 
											 << PAIR_TRIPLETS[k].ATMPAIR2    << " " 
											 << PAIR_TRIPLETS[k].ATMPAIR3    << " "
											 /* 
											 << sum.X << " " 
											 << sum.Y << " " 
											 << sum.Z << endl;
											 */
											 << PAIR_TRIPLETS[k].MIN_FOUND.X << " " 
											 << PAIR_TRIPLETS[k].MIN_FOUND.Y << " " 
											 << PAIR_TRIPLETS[k].MIN_FOUND.Z << endl;
			
			}
			if ( RANK == 0 ) cout << "	Total number of configurations contributing to each triplet type:" << endl;
		
			for (int k=0; k<PAIR_TRIPLETS.size(); k++)
			{
				int sum = 0 ;
#ifdef USE_MPI
				MPI_Reduce(&(PAIR_TRIPLETS[k].N_CFG_CONTRIB), &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
				sum = PAIR_TRIPLETS[k].N_CFG_CONTRIB ;
#endif				
				if ( RANK == 0 ) cout << "		" << k << "	" 					
											 << PAIR_TRIPLETS[k].ATMPAIR1 << " " 
											 << PAIR_TRIPLETS[k].ATMPAIR2 << " " 
											 << PAIR_TRIPLETS[k].ATMPAIR3 << " " 
											 << sum << endl;
			}
		}
		
		if( use_4b_cheby ) 
		{
			
			if ( RANK == 0 ) cout << "	Minimum distances between atoms quadruplet pairs: (Angstr.)" << endl;
			
			for (int k=0; k<PAIR_TRIPLETS.size(); k++) 
			{
				vector<double> sum(6,0.0);
				
				#ifdef USE_MPI
			
					for(int m=0; m<6; m++)
						MPI_Reduce(&(PAIR_QUADRUPLETS[k].MIN_FOUND[m]), &(sum[m]), 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
				#else	
					for(int m=0; m<6; m++)
						sum[m] = PAIR_QUADRUPLETS[k].MIN_FOUND[m];
				#endif
					
				if ( RANK == 0 ) 
				{
					cout << "		" << k << "	" ;	
					for(int m=0; m<6; m++)
						cout << PAIR_QUADRUPLETS[k].ATOM_PAIRS[m] << " ";
					for(int m=0; m<6; m++)
						cout << sum[m] << " ";
					cout << endl;
				}
			}
			
			if ( RANK == 0 ) cout << "	Total number of configurations contributing to each quadruplet type:" << endl;
		
			for (int k=0; k<PAIR_QUADRUPLETS.size(); k++)
			{
				int sum = 0 ;
				
				#ifdef USE_MPI
					MPI_Reduce(&(PAIR_QUADRUPLETS[k].N_CFG_CONTRIB), &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
				#else
					sum = PAIR_QUADRUPLETS[k].N_CFG_CONTRIB ;
				#endif	
								
				if ( RANK == 0 ) 
				{
					cout << "		" << k << "	" ;	
					for(int m=0; m<6; m++)
						cout << PAIR_QUADRUPLETS[k].ATOM_PAIRS[m] << " ";
					cout << sum << endl;
				}	
			}
		}

		if ( RANK == 0 ) cout << "...matrix printing complete: " << endl << endl;

}
  
bool operator==(const vector<int>& lhs, const vector<int>& rhs) 
{
  if ( lhs.size() != rhs.size() ) return false ;

  for ( int i = 0 ; i < lhs.size() ; i++ ) 
  {
	 if ( lhs[i] != rhs[i] ) return false ;
  }
  return true ;
}
