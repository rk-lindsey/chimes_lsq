
#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes
#include<string>
#include<limits>	// Help with handling of over/underflow
#include<algorithm> // Used for sorting, etc.
#include "functions.h"
#include "util.h"
#include "Cheby.h"

#ifdef USE_MPI
	#include <mpi.h>
#endif

using namespace std;


extern 	vector<int>	INT_PAIR_MAP;
extern	vector<int>	INT_TRIAD_MAP;	

//////////////////////////////////////////
//
//	SMALL UTILITY FUNCTION
//
//////////////////////////////////////////


//////////////////////////////////////////
// Functions to manage layers
//////////////////////////////////////////

void build_layers(FRAME &SYSTEM, JOB_CONTROL &CONTROLS)
{
	
	
	int     TEMP_IDX;
	XYZ     TEMP_XYZ;
	XYZ_INT TEMP_LAYER;
	
	SYSTEM.PARENT    .resize(SYSTEM.ATOMS);
	SYSTEM.LAYER_IDX .resize(SYSTEM.ATOMS);
	SYSTEM.ALL_COORDS.resize(SYSTEM.ATOMS);
	SYSTEM.WRAP_IDX  .resize(SYSTEM.ATOMS);
	
	// Set the first ATOMS atoms of ALL_COORDS equivalent to the "real" COORDS
	// Then wrap those ALL_COORDS into primitive cell
	
	SYSTEM.ALL_ATOMS = SYSTEM.ATOMS;	// Default setting: for zero layers
	
	for (int a1=0; a1<SYSTEM.ATOMS; a1++) 
	{
		SYSTEM.WRAP_IDX[a1].X = floor(SYSTEM.COORDS[a1].X/SYSTEM.BOXDIM.X);
		SYSTEM.WRAP_IDX[a1].Y = floor(SYSTEM.COORDS[a1].Y/SYSTEM.BOXDIM.Y);
		SYSTEM.WRAP_IDX[a1].Z = floor(SYSTEM.COORDS[a1].Z/SYSTEM.BOXDIM.Z);

		SYSTEM.ALL_COORDS[a1].X = SYSTEM.COORDS[a1].X - SYSTEM.WRAP_IDX[a1].X * SYSTEM.BOXDIM.X;
		SYSTEM.ALL_COORDS[a1].Y = SYSTEM.COORDS[a1].Y - SYSTEM.WRAP_IDX[a1].Y * SYSTEM.BOXDIM.Y;
		SYSTEM.ALL_COORDS[a1].Z = SYSTEM.COORDS[a1].Z - SYSTEM.WRAP_IDX[a1].Z * SYSTEM.BOXDIM.Z;
		
		SYSTEM.PARENT    [a1] = a1;
		SYSTEM.LAYER_IDX [a1].X = SYSTEM.LAYER_IDX [a1].Y = SYSTEM.LAYER_IDX [a1].Z = 0;
	}
	
	if(CONTROLS.N_LAYERS>0 )
	{	
		TEMP_IDX = SYSTEM.ATOMS;	

		for(int n1 = -CONTROLS.N_LAYERS; n1<=CONTROLS.N_LAYERS; n1++)
		{
			for(int n2 = -CONTROLS.N_LAYERS; n2<=CONTROLS.N_LAYERS; n2++)
			{
				for(int n3 = -CONTROLS.N_LAYERS; n3<=CONTROLS.N_LAYERS; n3++)
				{	
					if (n1 == 0 && n2 == 0 && n3 == 0 ) 
						continue;
					else
					{		
						for(int a1=0; a1<SYSTEM.ATOMS; a1++)
						{
							TEMP_XYZ.X = SYSTEM.ALL_COORDS[a1].X + n1 * SYSTEM.BOXDIM.X;
							TEMP_XYZ.Y = SYSTEM.ALL_COORDS[a1].Y + n2 * SYSTEM.BOXDIM.Y;
							TEMP_XYZ.Z = SYSTEM.ALL_COORDS[a1].Z + n3 * SYSTEM.BOXDIM.Z;
				
							TEMP_LAYER.X = n1;
							TEMP_LAYER.Y = n2;
							TEMP_LAYER.Z = n3;

							SYSTEM.ALL_COORDS	.push_back(TEMP_XYZ);
							SYSTEM.LAYER_IDX    .push_back(TEMP_LAYER);
							SYSTEM.ATOMTYPE     .push_back(SYSTEM.ATOMTYPE    [a1]);
							SYSTEM.CHARGES      .push_back(SYSTEM.CHARGES     [a1]);
							SYSTEM.PARENT       .push_back(a1);
				
							TEMP_IDX++;
						}
					}
				}
			}
		}
		SYSTEM.ALL_ATOMS = TEMP_IDX;
	}		
}


//////////////////////////////////////////
// Distance calculation and smoothing functions
//////////////////////////////////////////

double get_dist(const FRAME & SYSTEM, XYZ & RAB, int a1, int a2)
// Calculates distance as a2 - a1... This function modifies RAB!
{
	if(SYSTEM.ATOMS == SYSTEM.ALL_ATOMS) // Then we're not using ghost atoms. Need to use MIC
	{
		RAB.X = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
		RAB.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
		RAB.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
		
		RAB.X -= floor( 0.5 + RAB.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
		RAB.Y -= floor( 0.5 + RAB.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		RAB.Z -= floor( 0.5 + RAB.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;
	}
	else
	{
		RAB.X = SYSTEM.ALL_COORDS[a2].X - SYSTEM.ALL_COORDS[a1].X;
		RAB.Y = SYSTEM.ALL_COORDS[a2].Y - SYSTEM.ALL_COORDS[a1].Y;
		RAB.Z = SYSTEM.ALL_COORDS[a2].Z - SYSTEM.ALL_COORDS[a1].Z;
	}
	return sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );
}

//////////////////////////////////////////
// MPI compatibility functions
//////////////////////////////////////////
 
void divide_atoms(int &a1start, int &a1end, int atoms) 
{
	int procs_used;

	// Deal gracefully with more tasks than processors.
	if ( NPROCS <= atoms ) 
		procs_used = NPROCS;
	else
		procs_used = atoms;

	// Use ceil so the last process always has fewer tasks than the other
	// This improves load balancing.
	a1start = ceil( (double) RANK * atoms / procs_used);

	if ( RANK > atoms ) 
	{
		a1start = atoms + 1;
		a1end = atoms - 1;
	} else if ( RANK == procs_used - 1 ) {
		// End of the list.
		a1end = atoms - 1;
	} else {
		// Next starting value - 1 .
		a1end   = ceil( (double) (RANK+1) * atoms / procs_used ) - 1;
		if ( a1end > atoms - 1 ) 
			a1end = atoms - 1;
	}

	//cout << "DIVIDING ATOMS: RANK : " << RANK << " " << atoms << " " << a1start << ":" << a1end << endl;
}

//////////////////////////////////////////
// Kinetic energy functions
//////////////////////////////////////////

double kinetic_energy(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)					// UPDATED -- Overloaded.. compute differentely if for main or new velocities
{
	double Ktot = 0.0;
	
	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
		// Don't account for frozen atoms
		
		if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))
			continue;
		
		Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY[a1].X * SYSTEM.VELOCITY[a1].X;
		Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY[a1].Y * SYSTEM.VELOCITY[a1].Y;
		Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY[a1].Z * SYSTEM.VELOCITY[a1].Z;
	}

  return(Ktot);
}
double kinetic_energy(FRAME & SYSTEM, string TYPE, JOB_CONTROL & CONTROLS)		// UPDATED -- Overloaded.. compute differentely if for main or new velocities
{
	double Ktot = 0.0;
	
	if(TYPE == "NEW")
	{
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
			// Don't account for frozen atoms
		
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))
				continue;
			
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


//////////////////////////////////////////
// Pressure functions
//////////////////////////////////////////

void REPLICATE_SYSTEM(const FRAME & SYSTEM, FRAME & REPLICATE)
{
	// Set up for MPI
	
	int a1start, a1end, a2start, a2end;
	/*
	divide_atoms(a1start, a1end, SYSTEM.ATOMS);
	divide_atoms(a2start, a2end, SYSTEM.ALL_ATOMS);
	*/
	
	a1start = 0;
	a2start = 0;
	a1end   = REPLICATE.ATOMS;
	a2end   = REPLICATE.ALL_ATOMS;

	REPLICATE.ATOMS		= SYSTEM.ATOMS;
	REPLICATE.ALL_ATOMS	= SYSTEM.ALL_ATOMS;

	REPLICATE.BOXDIM.X	= SYSTEM.BOXDIM.X;
	REPLICATE.BOXDIM.Y	= SYSTEM.BOXDIM.Y;
	REPLICATE.BOXDIM.Z	= SYSTEM.BOXDIM.Z;

	REPLICATE.TOT_POT_ENER 			= SYSTEM.TOT_POT_ENER;
	
	REPLICATE.FORCES      .resize(REPLICATE.ATOMS);
	REPLICATE.ACCEL       .resize(REPLICATE.ATOMS);
	REPLICATE.COORDS      .resize(REPLICATE.ATOMS);

	for(int i=a1start; i<a1end; i++)
	{
		REPLICATE.COORDS[i].X = SYSTEM.COORDS[i].X;
		REPLICATE.COORDS[i].Y = SYSTEM.COORDS[i].Y;
		REPLICATE.COORDS[i].Z = SYSTEM.COORDS[i].Z;		
	}

	REPLICATE.ATOMTYPE    .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.ATOMTYPE_IDX.resize(REPLICATE.ALL_ATOMS);
	REPLICATE.PARENT      .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.ALL_COORDS  .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.CHARGES     .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.MASS        .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.LAYER_IDX   .resize(REPLICATE.ALL_ATOMS);

	for(int i=a2start; i<a2end; i++)
	{
		REPLICATE.ATOMTYPE    [i] = SYSTEM.ATOMTYPE[i];
		REPLICATE.ATOMTYPE_IDX[i] = SYSTEM.ATOMTYPE_IDX[i];
		REPLICATE.PARENT      [i] = SYSTEM.PARENT[i];
		REPLICATE.CHARGES     [i] = SYSTEM.CHARGES[i];
		REPLICATE.MASS        [i] = SYSTEM.MASS[i];

		REPLICATE.LAYER_IDX[i].X = SYSTEM.LAYER_IDX[i].X;
		REPLICATE.LAYER_IDX[i].Y = SYSTEM.LAYER_IDX[i].X;
		REPLICATE.LAYER_IDX[i].Z = SYSTEM.LAYER_IDX[i].X;		

		REPLICATE.ALL_COORDS[i].X = SYSTEM.ALL_COORDS[i].X;
		REPLICATE.ALL_COORDS[i].Y = SYSTEM.ALL_COORDS[i].Y;
		REPLICATE.ALL_COORDS[i].Z = SYSTEM.ALL_COORDS[i].Z;
	}
}

void numerical_pressure(const FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY,
								CLUSTER_LIST & TRIPS,  CLUSTER_LIST &QUADS, map<string,int> & PAIR_MAP,
								NEIGHBORS & NEIGHBOR_LIST,double & PE_1, double & PE_2, double & dV)
// Evaluates the configurational part of the pressure numerically by -dU/dV.
// Essentially, we are taking the system and expanding/contracting it a bit (lscale) to get the change in potential energy 
{
	const double eps = 1.0e-04;
	double lscale;
	double Vtot1, Vtot2;
	double Vol1, Vol2;

	// Make a copy of the system
	
	static FRAME REPLICATE;
		
	REPLICATE_SYSTEM(SYSTEM, REPLICATE);		// Can we make this one of those "if ! called_before sorts of variables?"

	// Expand coords by a bit (lscale)

	lscale = 1.0  + eps;
	
	for (int j=0; j < REPLICATE.ATOMS; j++ ) 
	{
	 	REPLICATE.COORDS[j].X = lscale * SYSTEM.COORDS[j].X;
		REPLICATE.COORDS[j].Y = lscale * SYSTEM.COORDS[j].Y;
		REPLICATE.COORDS[j].Z = lscale * SYSTEM.COORDS[j].Z;
	}

	for (int j=0; j < REPLICATE.ALL_ATOMS; j++ ) 
	{
	 	REPLICATE.ALL_COORDS[j].X = lscale * SYSTEM.ALL_COORDS[j].X;
		REPLICATE.ALL_COORDS[j].Y = lscale * SYSTEM.ALL_COORDS[j].Y;
		REPLICATE.ALL_COORDS[j].Z = lscale * SYSTEM.ALL_COORDS[j].Z;
	}
	
	REPLICATE.BOXDIM.X = lscale * SYSTEM.BOXDIM.X;
	REPLICATE.BOXDIM.Y = lscale * SYSTEM.BOXDIM.Y;
	REPLICATE.BOXDIM.Z = lscale * SYSTEM.BOXDIM.Z;

	// Compute/store new total potential energy and volume
	
	Vol1 = REPLICATE.BOXDIM.X * REPLICATE.BOXDIM.Y * REPLICATE.BOXDIM.Z;
	
	REPLICATE.TOT_POT_ENER = 0;
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, TRIPS, QUADS, PAIR_MAP, NEIGHBOR_LIST);
	
	Vtot1 = REPLICATE.TOT_POT_ENER;
	
	// Contract coords by a bit

	lscale = 1.0  - eps;
	
	for (int j=0; j < REPLICATE.ATOMS; j++ ) 
	{
	 	REPLICATE.COORDS[j].X = lscale * SYSTEM.COORDS[j].X;
		REPLICATE.COORDS[j].Y = lscale * SYSTEM.COORDS[j].Y;
		REPLICATE.COORDS[j].Z = lscale * SYSTEM.COORDS[j].Z;
	}

	for (int j=0; j < REPLICATE.ALL_ATOMS; j++ ) 
	{
	 	REPLICATE.ALL_COORDS[j].X = lscale * SYSTEM.ALL_COORDS[j].X;
		REPLICATE.ALL_COORDS[j].Y = lscale * SYSTEM.ALL_COORDS[j].Y;
		REPLICATE.ALL_COORDS[j].Z = lscale * SYSTEM.ALL_COORDS[j].Z;
	}
	
	REPLICATE.BOXDIM.X = lscale * SYSTEM.BOXDIM.X;
	REPLICATE.BOXDIM.Y = lscale * SYSTEM.BOXDIM.Y;
	REPLICATE.BOXDIM.Z = lscale * SYSTEM.BOXDIM.Z;
	
	// Compute/store new total potential energy and volume
	
	Vol2 = REPLICATE.BOXDIM.X * REPLICATE.BOXDIM.Y * REPLICATE.BOXDIM.Z;
	
	REPLICATE.TOT_POT_ENER = 0;
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, TRIPS, QUADS, PAIR_MAP, NEIGHBOR_LIST);
	
	Vtot2 = REPLICATE.TOT_POT_ENER;

	// compute (return) pressure 
	
	//return -(Vtot2 - Vtot1)/(Vol2 - Vol1);
	
	PE_1 = Vtot1;
	PE_2 = Vtot2;
	dV   = Vol2 - Vol1;
}


//////////////////////////////////////////
// Bad configuration tracking functions
//////////////////////////////////////////


void PRINT_CONFIG(FRAME &SYSTEM, JOB_CONTROL & CONTROLS)
// Output the current configuration
{
	int PRINT_WIDTH     = 21; // Use 21 for testing
	int PRINT_PRECISION = 14; // Use 14 for testing
	
	// Print out the file
	
	BAD_CONFIGS << fixed << setw(5) << setprecision(0) << SYSTEM.ATOMS << endl;
	BAD_CONFIGS << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.X << " ";
	BAD_CONFIGS << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.Y << " ";
	BAD_CONFIGS << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << SYSTEM.BOXDIM.Z << endl;
	
	for ( int ia = 0; ia < SYSTEM.ATOMS; ia++ ) 
	{
		XYZ tmp = SYSTEM.COORDS[ia];
		
		if ( CONTROLS.WRAP_COORDS ) 	// Wrap into the primitive cell
		{
			
			tmp.X -= floor(tmp.X/SYSTEM.BOXDIM.X)*SYSTEM.BOXDIM.X;
			tmp.Y -= floor(tmp.Y/SYSTEM.BOXDIM.Y)*SYSTEM.BOXDIM.Y;
			tmp.Z -= floor(tmp.Z/SYSTEM.BOXDIM.Z)*SYSTEM.BOXDIM.Z;
		}

		BAD_CONFIGS << setw(2) << SYSTEM.ATOMTYPE[ia] << " ";
		BAD_CONFIGS << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << tmp.X << " ";
		BAD_CONFIGS << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << tmp.Y << " ";
		BAD_CONFIGS << scientific << setw(PRINT_WIDTH) << setprecision(PRINT_PRECISION) << tmp.Z << "    ";

		BAD_CONFIGS << endl;
	}
}


//////////////////////////////////////////
//
//	FUNCTION HEADERS -- DERIVATIVE CALCULATION
//
//////////////////////////////////////////

static void ZCalc_Spline_Deriv  (JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

static void ZCalc_Poly_Deriv    (JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);	

static void ZCalc_InvR_Deriv    (JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);


//////////////////////////////////////////
//
//	FUNCTION HEADERS -- FORCE CALCULATION
//
//////////////////////////////////////////

static void ZCalc_Lj(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

static void ZCalc_Spline(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

static void ZCalcSR_Over(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

////////////////////////////////////////////////////////////
//
// FUNCTIONS -- DERIVATIVE CALCULATION
//
////////////////////////////////////////////////////////////

void ZCalc_Deriv (JOB_CONTROL & CONTROLS, vector<PAIRS> & FF_2BODY,  CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, 
						FRAME & FRAME_SYSTEM, vector<vector <XYZ > > & FRAME_A_MATRIX, 
						vector<vector <XYZ > > & FRAME_COULOMB_FORCES, const int nlayers, bool if_3b_cheby, 
						map<string,int> & PAIR_MAP,  NEIGHBORS &NEIGHBOR_LIST)
// Controls which functions are used to calculate derivatives
{
	// Check for control option compatability:
	
  vector<TRIPLETS> & PAIR_TRIPLETS = TRIPS.VEC ;
  map<string,int> &TRIAD_MAP = TRIPS.MAP ;

  map<int,int> &INT_QUAD_MAP = QUADS.INT_MAP ;


	if((CONTROLS.FIT_ENER) && (FF_2BODY[0].PAIRTYP == "SPLINE"))
	{
		cout << "ERROR: Energy fititng has not been implemented for potential type " << FF_2BODY[0].PAIRTYP << endl;
		exit_run(0);
	}
	
	if(CONTROLS.FIT_STRESS || CONTROLS.FIT_STRESS_ALL || CONTROLS.FIT_ENER)
	{
		if(CONTROLS.IF_SUBTRACT_COORD || CONTROLS.IF_SUBTRACT_COUL || CONTROLS.FIT_POVER)
		{
			cout << "ERROR: The following options must be false when including energies or stresses in fit: " << endl;
			cout << "CONTROLS.FIT_STRESS     " << endl;
			cout << "CONTROLS.FIT_STRESS_ALL" << endl;
			cout << "CONTROLS.FIT_ENER      " << endl;
		}
	}
	
	
  	// Will ewald calculations be needed? DFTB doesn't use the "default"
    // Ewald calculation because it has its own special way of dealing with charges
  
    if (FF_2BODY[0].PAIRTYP == "DFTBPOLY" && CONTROLS.FIT_COUL) 
	{
		cout << "ERROR: Turn of FITCOUL and set charges to zero when using pairtype " << FF_2BODY[0].PAIRTYP << endl;
		exit_run(0);
	}	

	if ( CONTROLS.FIT_COUL ) 
		ZCalc_Ewald_Deriv(FRAME_SYSTEM, FF_2BODY, FRAME_COULOMB_FORCES, PAIR_MAP, NEIGHBOR_LIST, CONTROLS);	
	
    if ( FF_2BODY[0].PAIRTYP == "SPLINE" )
		 ZCalc_Spline_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP, NEIGHBOR_LIST);
	
	else if ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" )
	{
		// Only enter if 2B are requested. For example, skip if user wants to fit ONLY 3B cheby
		// i.e. PAIRTYP: CHEBYSHEV  0 6 or similar
		
		if ( FF_2BODY[0].SNUM > 0)
			ZCalc_Cheby_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP, NEIGHBOR_LIST);
	
		if (if_3b_cheby)
			ZCalc_3B_Cheby_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, PAIR_TRIPLETS, FRAME_A_MATRIX,  nlayers, PAIR_MAP, TRIAD_MAP, NEIGHBOR_LIST);	
			
		if (CONTROLS.USE_4B_CHEBY)
		  ZCalc_4B_Cheby_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, TRIPS, QUADS, FRAME_A_MATRIX,  nlayers, PAIR_MAP, INT_QUAD_MAP, NEIGHBOR_LIST);		
	}			

    else if ( FF_2BODY[0].PAIRTYP == "DFTBPOLY" )	
		 ZCalc_Poly_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP, NEIGHBOR_LIST);

    else if ( FF_2BODY[0].PAIRTYP == "INVRSE_R" )	
		 ZCalc_InvR_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP, NEIGHBOR_LIST);

    else 
    {
		cout << "Error: bad pairtype in ZCalc_Deriv\n";
		exit(1);
    }		
}	

// FUNCTION UPDATED
static void ZCalc_Spline_Deriv(JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)		   
{
	// Original comment: Calculate derivatives of the forces wrt the spline parameters. Stores minimum distance between a pair of atoms in minD[i].
	// New Note: This doesn't actually calcuate any derivatives.. it is just populating A with the cubic hermite basis polynomials needed for fitting
	// 3rd note:  The hermite basis polynomials stored in FRAME_A_MATRIX are equal to derivatives of the force with respect to the spline coefficients..
	//            This is what all of the ZCalc_XX_Deriv functions find.
	//            Since the force is a linear function of the spline coefficients, you are left with just the cubic hermite basis polynomials.
	//            (LEF)
	
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double t;
	double h00,h10,h01,h11;
	int    k0;
	double x,x0;
	int    vstart,kstart;
	string TEMP_STR;
	
	double VOL = SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;

	// Main loop for SPLINE terms:
	
	int curr_pair_type_idx;
	
	// Set up for layering

	int fidx_a2;
	int a2start, a2end, a2;

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)		// Double sum over atom pairs
	{
		a2start = 0;
		a2end   = NEIGHBOR_LIST.LIST[a1].size();
		
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = NEIGHBOR_LIST.LIST[a1][a2idx];
			
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
			//calculate vstart:

			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Get pair distance

			rlen = get_dist(SYSTEM, RAB, a1, a2);	// Updates RAB!

			if ( (rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST) && (a1 !=a2) )	
				FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;
		   
		   
			if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM and rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM) //spline term calculated w/cutoff:
			{
				
				// Setting up for Hermite cubic quadrature:
				
			 	// This part: figure out which distance bin ("unit") that rlen falls into. Use the value of that bin "t" (rather than rlen ("x") ) to 
				// populate the matrix table. Use the lower and upper bounds of the bin to deterime the two constraining data points for the interpolation
				// polynomial fitting
				 
				 
				k0 = int( floor( (rlen - FF_2BODY[curr_pair_type_idx].S_MINIM) / FF_2BODY[curr_pair_type_idx].S_DELTA ) ); // binned location of distance along allowed range

			   
				if ( k0 >= FF_2BODY[curr_pair_type_idx].SNUM/2 - 1 ) // Keep from overrunning table.
				{
					cout << "Table overrun: rlen = " << rlen << endl;
					continue;
				}

			   	x = rlen;

			   	x0 = FF_2BODY[curr_pair_type_idx].S_MINIM + FF_2BODY[curr_pair_type_idx].S_DELTA * k0;
			   	t  = (x-x0) / FF_2BODY[curr_pair_type_idx].S_DELTA;

			   	// for classical cubic quadrature, you would only change h00, h10, h01, h11 polynomials!
			   
			  	// Setup the basis functions for cubic hermite polynomial interpolation...
			   	//
			   	// See: https://www3.nD[i].edu/~coast/jjwteach/www/www/30125/pdfnotes/lecture5_9v14.pdf
			   	// Look for "Each basis function is a third degree polynomial" 
			   	//
			   	// The fitting polynomial, g(x) = f(x_0)*alpha_0 + f(x_1)*alpha_1 + f'(x_0)*beta_0 + f'(x_1)*beta_1
			   	// But we're using t instead of x, and t_0 and t_1 are the two termini of the distance bin.
			   	// i.e. if x = 5.54, delta = 0.1, fall into a bin of value 5.55 with endpoints 5.5 and 5.6.
		   
			   	h00  =  2 * t * t * t - 3 * t * t + 1;			// alpha_0.. basis function for first point
			   
			   	h10  =  t * t * t  -2 * t * t + t;				// beta_0..  basis function for the first derivative at first point
			   	h10 *=  FF_2BODY[curr_pair_type_idx].S_DELTA;	// derivative terms have extra factor.
			   
			   	h01  = -2 * t * t * t + 3 * t * t;				// alpha_1.. basis function for second point
			   
			   	h11  =  t * t * t - t * t;						// beta_1..  basis function for first derivative at second point
			   	h11 *=  FF_2BODY[curr_pair_type_idx].S_DELTA;	// derivative terms have extra factor.

			   	kstart = k0 * 2;								//for each distance unit you have 2 entries. ...I'm guessing these are the start and end of the distance unit bin
			    
			    // Populate the A matrix with these values, muliplied times the xyz unit vectors to get directionality
				//
				// Note:
				//
				// vstart: takes care of index w/r/t possible pair types (i.e. OO, OH, HH...)
				// kstart: takes care of index w/r/t location along allowable pair distances
				// {0..3}: takes care of whether value corresponds to h00, h10, etc.
				//
				// So for a given item in a "Block" of A might correspond to, for example, type H-O, distance 5.0 along range 0.6 - 6.0 \AA, and basis function h01.	

				fidx_a2 = SYSTEM.PARENT[a2];
				
			   	FRAME_A_MATRIX[a1][vstart+kstart+0].X += h00 * RAB.X / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+1].X += h10 * RAB.X / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+2].X += h01 * RAB.X / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+3].X += h11 * RAB.X / rlen;
			  
			   	FRAME_A_MATRIX[a1][vstart+kstart+0].Y += h00 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+1].Y += h10 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+2].Y += h01 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+3].Y += h11 * RAB.Y / rlen;
				  
			   	FRAME_A_MATRIX[a1][vstart+kstart+0].Z += h00 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+1].Z += h10 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+2].Z += h01 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[a1][vstart+kstart+3].Z += h11 * RAB.Z / rlen;
				
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+0].X -= h00 * RAB.X / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+1].X -= h10 * RAB.X / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+2].X -= h01 * RAB.X / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+3].X -= h11 * RAB.X / rlen;
			  
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+0].Y -= h00 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+1].Y -= h10 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+2].Y -= h01 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+3].Y -= h11 * RAB.Y / rlen;
				  
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+0].Z -= h00 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+1].Z -= h10 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+2].Z -= h01 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+3].Z -= h11 * RAB.Z / rlen;
				
				if (CONTROLS.FIT_STRESS)
				{
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+0].X -= h00 * RAB.X * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+1].X -= h10 * RAB.X * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+2].X -= h01 * RAB.X * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+3].X -= h11 * RAB.X * RAB.X / rlen;
			  
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+0].Y -= h00 * RAB.Y * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+1].Y -= h10 * RAB.Y * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+2].Y -= h01 * RAB.Y * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+3].Y -= h11 * RAB.Y * RAB.Y / rlen;
				  
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+0].Z -= h00 * RAB.Z * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+1].Z -= h10 * RAB.Z * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+2].Z -= h01 * RAB.Z * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+kstart+3].Z -= h11 * RAB.Z * RAB.Z / rlen;					
				}
				else if (CONTROLS.FIT_STRESS_ALL)
				{
					// xx, xy, and xz
					
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+0].X -= h00 * RAB.X * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+1].X -= h10 * RAB.X * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+2].X -= h01 * RAB.X * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+3].X -= h11 * RAB.X * RAB.X / rlen;
			  
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+0].Y -= h00 * RAB.X * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+1].Y -= h10 * RAB.X * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+2].Y -= h01 * RAB.X * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+3].Y -= h11 * RAB.X * RAB.Y / rlen;
				  
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+0].Z -= h00 * RAB.X * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+1].Z -= h10 * RAB.X * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+2].Z -= h01 * RAB.X * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+kstart+3].Z -= h11 * RAB.X * RAB.Z / rlen;	
					
					
					// yx, yy, and yz
					
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+0].X -= h00 * RAB.Y * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+1].X -= h10 * RAB.Y * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+2].X -= h01 * RAB.Y * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+3].X -= h11 * RAB.Y * RAB.X / rlen;
			  
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+0].Y -= h00 * RAB.Y * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+1].Y -= h10 * RAB.Y * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+2].Y -= h01 * RAB.Y * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+3].Y -= h11 * RAB.Y * RAB.Y / rlen;
				  
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+0].Z -= h00 * RAB.Y * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+1].Z -= h10 * RAB.Y * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+2].Z -= h01 * RAB.Y * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+kstart+3].Z -= h11 * RAB.Y * RAB.Z / rlen;	
					
					
					// yx, yy, and yz
					
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+0].X -= h00 * RAB.Z * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+1].X -= h10 * RAB.Z * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+2].X -= h01 * RAB.Z * RAB.X / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+3].X -= h11 * RAB.Z * RAB.X / rlen;
			  
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+0].Y -= h00 * RAB.Z * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+1].Y -= h10 * RAB.Z * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+2].Y -= h01 * RAB.Z * RAB.Y / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+3].Y -= h11 * RAB.Z * RAB.Y / rlen;
				  
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+0].Z -= h00 * RAB.Z * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+1].Z -= h10 * RAB.Z * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+2].Z -= h01 * RAB.Z * RAB.Z / rlen;
				   	FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+kstart+3].Z -= h11 * RAB.Z * RAB.Z / rlen;	
				}							
				
			}
		}
	}
	
	if (CONTROLS.FIT_STRESS)
	{
		for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Z /= VOL;	
		}
	}
	
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Z /= VOL;	
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Z /= VOL;
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Z /= VOL;
		}
	}
	
	return;
}

// FUNCTION UPDATED
static void ZCalc_InvR_Deriv(JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)		
// Calculate derivatives of the forces wrt to inverse pair distance to various powers. Stores minimum distance between a pair of atoms in minD[i].
{
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	int vstart;
	int curr_pair_type_idx;
	string TEMP_STR;

	double fc;
	double dfc;
	double rfac;
	
	double VOL = SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;
	
	// Set up for layering

	int fidx_a2;
	int a2start, a2end, a2;
	int a1start, a1end;	

	//divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
	a1start = 0;
	a1end = SYSTEM.ATOMS-1;

	for(int a1=a1start; a1 <= a1end;a1++)		// Double sum over atom pairs
	{
		a2start = 0;
		a2end   = NEIGHBOR_LIST.LIST[a1].size();
			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{		
			a2 = NEIGHBOR_LIST.LIST[a1][a2idx];

			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
		
			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Get pair distance

			rlen = get_dist(SYSTEM, RAB, a1, a2);	// Updates RAB!

			if ( (rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST) && (a1 !=a2) )	// spline term calculated w/cutoff:
				FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

			if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM && rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
			{
				// Calculate the penalty function (fc) and its derivative
				
				fc = (rlen - FF_2BODY[curr_pair_type_idx].S_MAXIM);
				dfc = fc;
				
				fc = fc*fc*fc;							
				dfc = 3*dfc*dfc;
				
				fidx_a2 = SYSTEM.PARENT[a2];
			
				for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
				{
					rfac = ( (i+2) / pow(rlen,i+3) )*fc;
					rfac -= (1/pow(rlen,i+2))*dfc;

					FRAME_A_MATRIX[a1][vstart+i].X += rfac * RAB.X/rlen;
					FRAME_A_MATRIX[a1][vstart+i].Y += rfac * RAB.Y/rlen;
					FRAME_A_MATRIX[a1][vstart+i].Z += rfac * RAB.Z/rlen;
						
					FRAME_A_MATRIX[fidx_a2][vstart+i].X -= rfac * RAB.X/rlen;
					FRAME_A_MATRIX[fidx_a2][vstart+i].Y -= rfac * RAB.Y/rlen;
					FRAME_A_MATRIX[fidx_a2][vstart+i].Z -= rfac * RAB.Z/rlen;		
					
					if (CONTROLS.FIT_STRESS)
					{
						FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].X -= rfac * RAB.X * RAB.X / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Y -= rfac * RAB.Y * RAB.Y / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Z -= rfac * RAB.Z * RAB.Z / rlen;	
					}
					
					else if (CONTROLS.FIT_STRESS_ALL)
					{
						FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].X -= rfac * RAB.X * RAB.X / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Y -= rfac * RAB.X * RAB.Y / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Z -= rfac * RAB.X * RAB.Z / rlen;	
						
						FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].X -= rfac * RAB.Y * RAB.X / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Y -= rfac * RAB.Y * RAB.Y / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Z -= rfac * RAB.Y * RAB.Z / rlen;	
						
						FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].X -= rfac * RAB.Z * RAB.X / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Y -= rfac * RAB.Z * RAB.Y / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Z -= rfac * RAB.Z * RAB.Z / rlen;	
					}						
				}
			}
		}
	}
	
	if (CONTROLS.FIT_STRESS)
	{
		for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Z /= VOL;	
		}
	}
	
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Z /= VOL;	
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Z /= VOL;
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Z /= VOL;
		}
	}
	
	return;

}

// FUNCTION UPDATED
static void ZCalc_Poly_Deriv(JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)	
// Calculate derivatives of the forces wrt the DFTB Erep parameters.
{	
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double x;
	int vstart;

	int curr_pair_type_idx;
	string TEMP_STR;

	static bool called_before = false;
	const  double autoang = 0.5291772488820865;	// Conversion factor fro angstroms to au
	static double *rc;
	
	int dim;
	double rfac;
	
	double VOL = SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;

	if ( ! called_before ) 
	{
		called_before = true;
		dim = 0;
		
		rc = new double [FF_2BODY.size()];
		
		for ( int i=0; i<FF_2BODY.size(); i++ ) 
		{
			// Convert cutoff (S_MAXIM) from units of Angstrom to bohr (atomic units), (rc)
			
			rc[i] = FF_2BODY[i].S_MAXIM/autoang;
			if ( RANK == 0 ) cout << "	rc[" << i << "] = " 
										 << fixed << setprecision(3) 
										 << rc[i] << " bohr (" 
										 << fixed << setprecision(3) 
										 << FF_2BODY[i].S_MAXIM << " Angstr.)" 
										 << endl;
			
			if ( FF_2BODY[i].SNUM > dim ) 
				dim = FF_2BODY[i].SNUM;
		}
		dim++;
	
	}

	// main loop for ninth order polynomial terms:
	
	// Set up for layering

	int fidx_a2;
	int a2start, a2end, a2;
	int a1start, a1end;	


	//divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
	a1start = 0;
	a1end = SYSTEM.ATOMS-1;

	
	for(int a1= a1start;a1<= a1end; a1++)		// Double sum over atom pairs
	{
		a2start = 0;

		a2end   = NEIGHBOR_LIST.LIST[a1].size();
			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = NEIGHBOR_LIST.LIST[a1][a2idx];	

			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
		
			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Get pair distance

			rlen = get_dist(SYSTEM, RAB, a1, a2);	// Updates RAB!
						
			if ( (rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST) && (a1 !=a2) )
				FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

			if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM && rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
			{
				// calculate binning, convert all distances to au from angstroms 
				x = rlen/autoang;
				
				fidx_a2 = SYSTEM.PARENT[a2];

				for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
				{
					rfac = -(i+2)*pow((rc[curr_pair_type_idx]-x),i+1);

					FRAME_A_MATRIX[a1][vstart+i].X += rfac * RAB.X/rlen;
					FRAME_A_MATRIX[a1][vstart+i].Y += rfac * RAB.Y/rlen;
					FRAME_A_MATRIX[a1][vstart+i].Z += rfac * RAB.Z/rlen;
						
					FRAME_A_MATRIX[fidx_a2][vstart+i].X -= rfac * RAB.X/rlen;
					FRAME_A_MATRIX[fidx_a2][vstart+i].Y -= rfac * RAB.Y/rlen;
					FRAME_A_MATRIX[fidx_a2][vstart+i].Z -= rfac * RAB.Z/rlen;
					
					if (CONTROLS.FIT_STRESS)
					{
						FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].X -= 29421.02407027691 * rfac * RAB.X * RAB.X / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Y -= 29421.02407027691 * rfac * RAB.Y * RAB.Y / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Z -= 29421.02407027691 * rfac * RAB.Z * RAB.Z / rlen;	
					}
					
					else if (CONTROLS.FIT_STRESS_ALL)
					{
						FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].X -= 29421.02407027691 * rfac * RAB.X * RAB.X / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Y -= 29421.02407027691 * rfac * RAB.X * RAB.Y / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Z -= 29421.02407027691 * rfac * RAB.X * RAB.Z / rlen;	
						
						FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].X -= 29421.02407027691 * rfac * RAB.Y * RAB.X / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Y -= 29421.02407027691 * rfac * RAB.Y * RAB.Y / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Z -= 29421.02407027691 * rfac * RAB.Y * RAB.Z / rlen;	
						
						FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].X -= 29421.02407027691 * rfac * RAB.Z * RAB.X / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Y -= 29421.02407027691 * rfac * RAB.Z * RAB.Y / rlen;
						FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Z -= 29421.02407027691 * rfac * RAB.Z * RAB.Z / rlen;	
					}
				}
			}

		}
	}
	
	if (CONTROLS.FIT_STRESS)
	{
		for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS][vstart+i].Z /= VOL;	
		}
	}
	
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
		{
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS  ][vstart+i].Z /= VOL;	
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+1][vstart+i].Z /= VOL;
		
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].X /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Y /= VOL;
			FRAME_A_MATRIX[SYSTEM.ATOMS+2][vstart+i].Z /= VOL;
		}
	}
	
	
	return;

}

// FUNCTION UPDATED
void SubtractCoordForces(FRAME & SYSTEM, bool calc_deriv, vector<XYZ> & P_OVER_FORCES,  
								 vector<PAIRS> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST,
								 bool lsq_mode)
// lsq_mode is true if this is a least-squares calculation.
{
	// this function subtracts the ReaxFF over-coordination term to re-fit 
	// splines/charges iteratively for self-consistence.
	// If calc_deriv is true, the derivative of the force wrt the magnitude of the
	// 3-body interaction is placed in Fderiv.  Otherwise, the 3-body force is subtracted
	// from the total forces given in force.  
	 

	// Lucas paper draft parameters.	// Chenoweth values for p1, p2, r0, lambda6 p(ovun2)
	//									// pover = 50.0;	= FF_2BODY[j].OVRPRMS[0];
	// p1  	  =	-2.5881042987450e-01;	// -0.0657;			= FF_2BODY[j].OVRPRMS[1];
	// r0 	  =	 9.6000387075695e-01;	//  5.0451;			= FF_2BODY[j].OVRPRMS[2];
	// p2 	  =	 3.8995379237250e+00;	//  1.0165;			= FF_2BODY[j].OVRPRMS[3];
	// lambda6=	-8.9;  					// -3.6141;			= FF_2BODY[j].OVRPRMS[4];

	XYZ RVEC; 		// Replaces Rvec[3];
	double Vtot = 0;
	
	vector<XYZ> SForce(SYSTEM.ATOMS);	// This is the force ( negative of dEover ) (LEF)
	vector<XYZ> dEover(SYSTEM.ATOMS);	// Over bonding energy derivative (force) w/r/t delta ("S")
	vector<XYZ> dFover(SYSTEM.ATOMS);	// Over bonding force derivative w/r/t pover
	
	// THREE-BODY POTENTIAL:
	
	double Eover = 0;			// Overbonding energy
	double rik,tempr,temps;	
	double S[SYSTEM.ATOMS];				// Difference between oxygen bond order and valency (Called delta in the paper)

	int curr_pair_type_idx; 
	string TEMP_STR;
	
	// Set up some variables
	
	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
		dEover[a1].X = 0.0;
		dEover[a1].Y = 0.0;
		dEover[a1].Z = 0.0;
		
		dFover[a1].X = 0.0;
		dFover[a1].Y = 0.0;
		dFover[a1].Z = 0.0;		
		
		SForce[a1].X = 0;
		SForce[a1].Y = 0;
		SForce[a1].Z = 0;
		
		S[a1] = 0.0;				
	}
		
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Calculate the overbonding energy... Derivative needs to be in separate double loop b/c depends
	// on a term calculated with energies ("S", which is called delta in the paper)
	//
	////////////////////////////////////////////////////////////////////////////////////////////////
			
	// Set up for MPI
	
	int aistart, aiend, ak;	

	if (  lsq_mode ) {
		aistart = 0;
		aiend = SYSTEM.ATOMS-1;
	} else {
		divide_atoms(aistart, aiend, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
	}
			
	for(int ai=aistart; ai <= aiend; ai++) 
	{
		temps = 0.0;
	  
		int a2start = 0;
		int a2end = NEIGHBOR_LIST.LIST_UNORDERED[ai].size();

		for(int akidx = a2start; akidx < a2end; akidx++ )
		{			
			ak = NEIGHBOR_LIST.LIST_UNORDERED[ai][akidx];			
			// TWO-BODY PART... ONLY CARES ABOUT LOOPS OVER AI AND AK
		
			TEMP_STR = SYSTEM.ATOMTYPE[ai];
			TEMP_STR.append(SYSTEM.ATOMTYPE[ak]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
			
			// THE ITEMS COMMENTED OUT IN THIS LOOP SHOULD BE PUT BACK ONCE CODE COMPARISON WITH LUCAS' COMPLETE!!!!!
				
			if(FF_2BODY[curr_pair_type_idx].USE_OVRPRMS && (SYSTEM.ATOMTYPE[ai] == FF_2BODY[curr_pair_type_idx].OVER_TO_ATM)) // Then we should have a defined pair type
			{
				// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
	   
				rik = get_dist(SYSTEM, RVEC, ak, ai);
			
			// Calculate the O--H bond order as defined by ReaxFF
			
				temps +=  exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[2]
				*pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3])); 
			}
		}
		
		S[ai]  = temps - 2.0; // This is the reaxff delta (diff between oxygen valency and bond order) 
		Eover += FF_2BODY[curr_pair_type_idx].OVRPRMS[0] * S[ai] * 1.0/(1.0+exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai]));	

	}

	Vtot += Eover;
	
	///////////////////////////////////////////////////////
	//
	// Now calculate the forces (derivative)
	//
	///////////////////////////////////////////////////////

	// Here, the first loop is over all ai, but we don't do anything unless ai is an oxygen type, since overbonding forces in this program
	// are defined to be between oxygen an either oxygen or H 
	// Second loop is over all atoms != ai
	// For force between a pair of atoms, recall that vectors are defined from ai to ak

	for(int ai=0; ai<SYSTEM.ATOMS; ai++)
	{
		int a2start = 0;
		int a2end = NEIGHBOR_LIST.LIST_UNORDERED[ai].size();

		for(int akidx = a2start; akidx < a2end; akidx++ )
		{	
			ak = NEIGHBOR_LIST.LIST_UNORDERED[ai][akidx];
				
			TEMP_STR = SYSTEM.ATOMTYPE[ai];
			TEMP_STR.append(SYSTEM.ATOMTYPE[ak]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
			
			// If it's a pair type with defined over params
			// and the pair isn't comprised of the exact same atom
			// and ai is the "to" atom (i.e. O in an O--H pair)
		
			if((FF_2BODY[curr_pair_type_idx].USE_OVRPRMS && ai != ak && SYSTEM.ATOMTYPE[ai] == FF_2BODY[curr_pair_type_idx].OVER_TO_ATM))
			{

				rik = get_dist(SYSTEM, RVEC, ai, ak);			
				
				// Calculate the derivative of Eover w/r/t delta. Terms emerge from quotient rule combined with chain rule.

				tempr = 1.0/(1+exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai])) - FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai]*exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai])
					        / pow(1.0+exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai]),2);

				tempr *= FF_2BODY[curr_pair_type_idx].OVRPRMS[2]*pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3]) *
					     FF_2BODY[curr_pair_type_idx].OVRPRMS[3]*exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[2] *
						 pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3]))/rik;	


				dEover[ai].X -= tempr*RVEC.X/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				dEover[ai].Y -= tempr*RVEC.Y/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				dEover[ai].Z -= tempr*RVEC.Z/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];						

				int fidx_ak = SYSTEM.PARENT[ak];

				dEover[fidx_ak].X += tempr*RVEC.X/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				dEover[fidx_ak].Y += tempr*RVEC.Y/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				dEover[fidx_ak].Z += tempr*RVEC.Z/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];	

				dFover[ai].X -= tempr*RVEC.X/rik;
				dFover[ai].Y -= tempr*RVEC.Y/rik;
				dFover[ai].Z -= tempr*RVEC.Z/rik;						

				dFover[fidx_ak].X += tempr*RVEC.X/rik;
				dFover[fidx_ak].Y += tempr*RVEC.Y/rik;
				dFover[fidx_ak].Z += tempr*RVEC.Z/rik;					
			}	
		}
	}

    for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
        SForce[a1].X -= dEover[a1].X;
		SForce[a1].Y -= dEover[a1].Y;
		SForce[a1].Z -= dEover[a1].Z;	
	}

    if ( !calc_deriv) // subtract CoorD[i]. force:
	{ 
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
	        SYSTEM.FORCES[a1].X -= SForce[a1].X;
			SYSTEM.FORCES[a1].Y -= SForce[a1].Y;
			SYSTEM.FORCES[a1].Z -= SForce[a1].Z;	
		}
    }
    else // Then we're fitting pover... Calculate derivative of 3-body force wrt pover.
	{
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
	        P_OVER_FORCES[a1].X -= dFover[a1].X;
			P_OVER_FORCES[a1].Y -= dFover[a1].Y;
			P_OVER_FORCES[a1].Z -= dFover[a1].Z;		
		}
	} // VERIFIED FROM HERE UP
}


////////////////////////////////////////////////////////////
//
// FUNCTIONS -- FORCE CALCULATION
//
////////////////////////////////////////////////////////////
 

void ZCalc(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, 
			  map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
{
  
	for(int a=0;a<SYSTEM.ATOMS;a++)
	{
		SYSTEM.ACCEL[a].X = 0;
		SYSTEM.ACCEL[a].Y = 0;
		SYSTEM.ACCEL[a].Z = 0;
	}

	SYSTEM.TOT_POT_ENER = 0;
	SYSTEM.PRESSURE_XYZ = 0;
	
	// XY, YZ, and XZ need to be included at a future date ...
	
	SYSTEM.PRESSURE_TENSORS_XYZ.X = 0;
	SYSTEM.PRESSURE_TENSORS_XYZ.Y = 0;
	SYSTEM.PRESSURE_TENSORS_XYZ.Z = 0;

	if      ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" ) 
		ZCalc_Cheby_ALL(SYSTEM, CONTROLS, FF_2BODY, TRIPS, QUADS, PAIR_MAP, NEIGHBOR_LIST);
	else if ( FF_2BODY[0].PAIRTYP == "LJ" ) 
		ZCalc_Lj(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, NEIGHBOR_LIST);
	
	else if ( FF_2BODY[0].PAIRTYP == "SPLINE" ) 
		ZCalc_Spline(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, NEIGHBOR_LIST);	
	
    else 
    {
		cout << "Error: bad pairtype in ZCalc: " << FF_2BODY[0].PAIRTYP << endl;
		exit_run(1);
    }	
	
	if ( CONTROLS.USE_COULOMB ) 
		ZCalc_Ewald(SYSTEM, CONTROLS, NEIGHBOR_LIST);
		
	if ( CONTROLS.USE_OVERCOORD ) 	
		ZCalcSR_Over(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, NEIGHBOR_LIST);

	// FUNCTIONS THAT NEED UPDATING:

	/* 

	else if ( pair_type == INVERSE_R ) 	
		ZCalc_SR_Analytic(Coord,Lbc, Latcons,nlayers,nat,smin,smax,snum, SForce,Vtot, Pxyz, params);

	else if ( pair_type == STILLINGER )  // SAVE THIS FOR SECOND TO LAST FOR SIMILAR REASONS  
		ZCalc_Stillinger(Coord,Lbc, Latcons,nlayers,nat,smax, SForce,Vtot,Pxyz);

	else 
		EXIT_MSG("Error: Unknown pair type", pair_type)
	*/
	
	SYSTEM.PRESSURE_XYZ           /= 3.0 *  SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;
	SYSTEM.PRESSURE_TENSORS_XYZ.X /=        SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;
	SYSTEM.PRESSURE_TENSORS_XYZ.Y /=        SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;
	SYSTEM.PRESSURE_TENSORS_XYZ.Z /=        SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;

  return;
}

static void ZCalc_Lj(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
// Calculate LJ interaction.. first parameter is epsilon, second parameter is sigma. ...eventually SMAX should be used for the pair distance cutoff value...
{
  XYZ	RVEC ;
	double	rlen_mi;
	int	curr_pair_type_idx;
	double	fac;
	string	TEMP_STR;
	
	// Set up for MPI
	
	int a1start, a1end;
	int a2start, a2end;
	int fidx_a2;

	#ifndef LINK_LAMMPS
		divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
	#else
		a1start = SYSTEM.MY_ATOMS_START;
		a1end   = SYSTEM.MY_ATOMS_START+SYSTEM.MY_ATOMS-1;
	#endif

	for(int a1=a1start;a1 <= a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS 
	{
		a2start = 0;
		a2end   = NEIGHBOR_LIST.LIST[a1].size();
		
		for(int a2idx =a2start; a2idx < a2end;a2idx++)
		{			
			int a2 = NEIGHBOR_LIST.LIST[a1][a2idx];
					
			#ifndef LINK_LAMMPS
				curr_pair_type_idx =  INT_PAIR_MAP[SYSTEM.ATOMTYPE_IDX[a1]*CONTROLS.NATMTYP + SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]];
			#else
				curr_pair_type_idx =  INT_PAIR_MAP[(SYSTEM.ATOMTYPE_IDX[a1]-1)*CONTROLS.NATMTYP + (SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]-1)];
			#endif

			// pair interaction cutoff distance.
			double rcutoff = FF_2BODY[curr_pair_type_idx].S_MAXIM;

			rlen_mi = get_dist(SYSTEM, RVEC, a1, a2);
	
			if ( rlen_mi < FF_2BODY[curr_pair_type_idx].PARAMS[1]/2.2) 
				EXIT_MSG("Error: close approach", rlen_mi);

			else if ( rlen_mi < rcutoff ) 
			{
				SYSTEM.TOT_POT_ENER += 4.0 * FF_2BODY[curr_pair_type_idx].PARAMS[0] * ( pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,12.0) - pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,6.0) );
				fac = 4.0 * FF_2BODY[curr_pair_type_idx].PARAMS[0] * ( -12.0 * pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,14.0) + 6.0 *    pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,8.0) );
				fac *= 1.0 / ( FF_2BODY[curr_pair_type_idx].PARAMS[1] * FF_2BODY[curr_pair_type_idx].PARAMS[1] );		

				SYSTEM.PRESSURE_XYZ -= fac * (rlen_mi*rlen_mi);
				
				SYSTEM.PRESSURE_TENSORS_XYZ.X -= fac * rlen_mi * RVEC.X * RVEC.X / rlen_mi;
				SYSTEM.PRESSURE_TENSORS_XYZ.Y -= fac * rlen_mi * RVEC.Y * RVEC.Y / rlen_mi;
				SYSTEM.PRESSURE_TENSORS_XYZ.Z -= fac * rlen_mi * RVEC.Z * RVEC.Z / rlen_mi;
	
				SYSTEM.ACCEL[a1].X += RVEC.X*fac;
				SYSTEM.ACCEL[a1].Y += RVEC.Y*fac;
				SYSTEM.ACCEL[a1].Z += RVEC.Z*fac;

				fidx_a2 = SYSTEM.PARENT[a2];

				SYSTEM.ACCEL[fidx_a2].X -= RVEC.X*fac;
				SYSTEM.ACCEL[fidx_a2].Y -= RVEC.Y*fac;
				SYSTEM.ACCEL[fidx_a2].Z -= RVEC.Z*fac;
			}			
		}
	}
}

static void ZCalc_Spline(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
// Calculate spline forces.
{
	XYZ RAB;
	double rlen,rlen2;
	double tempx;
	double S_r;
	
	int k0;
	double x, x0;
	double t, t2, t3, t4;;
	double h00,h10,h01,h11;
	double i00,i10,i01,i11;
	int kstart;
	
	string TEMP_STR;
	int curr_pair_type_idx;

	// Set up for MPI
	
	int a1start, a1end, fidx;
	
	#ifndef LINK_LAMMPS
		divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
	#else
		a1start = SYSTEM.MY_ATOMS_START;
		a1end   = SYSTEM.MY_ATOMS_START+SYSTEM.MY_ATOMS-1;
	#endif

	// Set up for neighbor list
	
	int a2start, a2end, a2;

	for(int a1=a1start;a1 <= a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev-1) 
	{
		a2start = 0;
		a2end   = NEIGHBOR_LIST.LIST[a1].size();
		
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{	
			a2 = NEIGHBOR_LIST.LIST[a1][a2idx];

			rlen = get_dist(SYSTEM, RAB, a1, a2);	// Updates RAB!
								
			#ifndef LINK_LAMMPS
				curr_pair_type_idx =  INT_PAIR_MAP[SYSTEM.ATOMTYPE_IDX[a1]*CONTROLS.NATMTYP + SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]];
			#else
				curr_pair_type_idx =  INT_PAIR_MAP[(SYSTEM.ATOMTYPE_IDX[a1]-1)*CONTROLS.NATMTYP + (SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]-1)];
			#endif
				
			if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM && rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
			{			

				tempx = 0;
				rlen2 = rlen;		

			  	k0  = int(floor((rlen2-FF_2BODY[curr_pair_type_idx].S_MINIM)/FF_2BODY[curr_pair_type_idx].S_DELTA));
			  	x   = rlen2;
			  	x0  = FF_2BODY[curr_pair_type_idx].S_MINIM+FF_2BODY[curr_pair_type_idx].S_DELTA*k0;
			  	t   = (x-x0)/FF_2BODY[curr_pair_type_idx].S_DELTA;

			  	if ( t > 1.000001 || t < -0.000001 ) 
				  	EXIT_MSG("Error: Bad t in spline calculation.");

				t2 = t * t;
				t3 = t2 * t;
				t4 = t3 * t;

				h00 =  2*t3 -3*t2 + 1.0;
				h10 =    t3 -2*t2 + t;
				h01 = -2*t3 +3*t2;
				h11 =    t3 -  t2;
				h10 *= FF_2BODY[curr_pair_type_idx].S_DELTA; 		//derivative terms have extra factor.
				h11 *= FF_2BODY[curr_pair_type_idx].S_DELTA;
				
				kstart = k0*2;	
			  	    
				if ( kstart > FF_2BODY[curr_pair_type_idx].SNUM ) 
				{
					cout << "Error: kstart too large " << kstart << " " << rlen << " " << FF_2BODY[curr_pair_type_idx].S_MINIM << " " << FF_2BODY[curr_pair_type_idx].SNUM << endl;
					exit_run(1);
				}
			  
				S_r =
					h00 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+0]+//lookout!!!
					h10 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+1]+
					h01 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+2]+
					h11 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+3];

				//do integrating here to get potential energy:
				i00 = 0.5 - (0.5*t4 - t3 + t);
				i10 = 1.0/12.0 - (0.25*t4- (2.0/3.0)*t3 + 0.5*t2);
				i01 = 0.5 - (-0.5*t4 + t3);
				i11 = -1.0/12.0 - (0.25*t4 - t3/3.0);
				i10 *= FF_2BODY[curr_pair_type_idx].S_DELTA;
				i11 *= FF_2BODY[curr_pair_type_idx].S_DELTA;

				// TEST: Linear interpolation.
				// i00 = (0.5 - (t - t * t / 2.0));
				// i01 = (0.5 - t * t / 2.0);
				// i10 = 0.0;
				// i11 = 0.0;

				tempx =  -FF_2BODY[curr_pair_type_idx].S_DELTA * (
					i00 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+0]+
					i10 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+1]+
					i01 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+2]+
					i11 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+3] );
			
				// Use the integral of the spline for pressure calculation
				
				if ( kstart/2 + 1 < FF_2BODY[curr_pair_type_idx].SNUM / 2 )
					tempx += FF_2BODY[curr_pair_type_idx].POT_PARAMS[kstart/2+1];
				
				fidx = SYSTEM.PARENT[a2];

				SYSTEM.ACCEL[a1].X += S_r * RAB.X/rlen;
				SYSTEM.ACCEL[a1].Y += S_r * RAB.Y/rlen;
				SYSTEM.ACCEL[a1].Z += S_r * RAB.Z/rlen;

				SYSTEM.ACCEL[fidx].X -= S_r * RAB.X/rlen;
				SYSTEM.ACCEL[fidx].Y -= S_r * RAB.Y/rlen;
				SYSTEM.ACCEL[fidx].Z -= S_r * RAB.Z/rlen;
	   
				SYSTEM.PRESSURE_XYZ -= S_r * rlen;
				
				SYSTEM.PRESSURE_TENSORS_XYZ.X -= S_r * RAB.X * RAB.X / rlen;
				SYSTEM.PRESSURE_TENSORS_XYZ.Y -= S_r * RAB.Y * RAB.Y / rlen;
				SYSTEM.PRESSURE_TENSORS_XYZ.Z -= S_r * RAB.Z * RAB.Z / rlen;
				
				SYSTEM.TOT_POT_ENER += tempx;

			}//rlen
		}
	}
  return;
}

// UPDATED AND VERIFIED
static void ZCalcSR_Over(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
// Calculate short-range overcoordination forces... See SubtractCoordForces for more info
{

	// Notes:
	//
  	// these are short-ranged cutoffs for analytical potential.
  	// e.g. for r<oo_cut, force(r)=force(oo_cut). The values are 
  	// at least 0.2 Angs below any sampled by MD; they are there 
  	// in case trajectories reach them.
	//
  	// Many-body overcoordination term:
  	// The rest of this function can be commented out and the MD 
  	// will still run. The ReaxFF terms have been fitted using
  	// Powell minimization (separate routines, but very easy to do).
	
	XYZ RVEC, RAB;

	// MANY-BODY POTENTIAL:
	double Eover = 0;
	double rik,tempr,temps;
	double powrik;
	double Sexp2;
	
	vector<double> S     (SYSTEM.ATOMS);
	vector<double> Sexp  (SYSTEM.ATOMS);
	vector<XYZ>    dEover(SYSTEM.ATOMS);
	
	int curr_pair_type_idx; 
	string TEMP_STR;
	
	bool SAFE = false;

	// Set up for MPI
	
	int aistart, aiend;	
	
	#ifndef LINK_LAMMPS
		divide_atoms(aistart, aiend, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
	#else
		aistart = SYSTEM.MY_ATOMS_START;
		aiend   = SYSTEM.MY_ATOMS_START+SYSTEM.MY_ATOMS-1;
	#endif

	for(int ai=aistart;ai <= aiend; ai++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS
	{
		temps=0.0;
		S[ai] = 0.0;
		Sexp[ai] = 0.0;

		for(int i=0; i<FF_2BODY.size(); i++)
			if(SYSTEM.ATOMTYPE[ai] == FF_2BODY[i].OVER_TO_ATM)
				SAFE = true;
		
		if(SAFE)	
		{
			int a2start = 0;
			int a2end = NEIGHBOR_LIST.LIST_UNORDERED[ai].size();

			for(int akidx = a2start; akidx < a2end; akidx++ )
			{
				int ak = NEIGHBOR_LIST.LIST_UNORDERED[ai][akidx];

				#ifndef LINK_LAMMPS
					curr_pair_type_idx =  INT_PAIR_MAP[SYSTEM.ATOMTYPE_IDX[ai]*CONTROLS.NATMTYP + SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[ak]]];
				#else
					curr_pair_type_idx =  INT_PAIR_MAP[(SYSTEM.ATOMTYPE_IDX[ai]-1)*CONTROLS.NATMTYP + (SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[ak]]-1)];
				#endif
					
				if(FF_2BODY[curr_pair_type_idx].USE_OVRPRMS)  // Then we should have a defined pair type
				{									
					rik = get_dist(SYSTEM, RAB, ai, ak);
					temps += exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[2] * pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3])); 
				}				
			}
			
			S[ai]    = temps-2.0;
			Sexp[ai] = exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai]);
			Eover   += FF_2BODY[curr_pair_type_idx].OVRPRMS[0] * S[ai] * 1.0/(1.0+Sexp[ai]);			
		}
		
		SAFE = false;
	}
	
	SYSTEM.TOT_POT_ENER += Eover;
	
	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		dEover[a1].X = dEover[a1].Y = dEover[a1].Z = 0;

	for(int ai=aistart;ai <= aiend; ai++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS
	{
		int a2start = 0;
		int a2end = NEIGHBOR_LIST.LIST_UNORDERED[ai].size();

		for(int akidx = a2start; akidx < a2end; akidx++ )
		{	
			int ak = NEIGHBOR_LIST.LIST_UNORDERED[ai][akidx];

			curr_pair_type_idx =  INT_PAIR_MAP[SYSTEM.ATOMTYPE_IDX[ai]*CONTROLS.NATMTYP + SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[ak]]];

			if((FF_2BODY[curr_pair_type_idx].USE_OVRPRMS && SYSTEM.ATOMTYPE[ai] == FF_2BODY[curr_pair_type_idx].OVER_TO_ATM))
			{

				rik = get_dist(SYSTEM, RVEC, ai, ak);

				// Calculate the derivative of Eover w/r/t delta. Terms emerge from quotient rule combined with chain rule.

				powrik = pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3]);
				Sexp2  = (1.0 + Sexp[ai]) * (1.0 + Sexp[ai]);

				tempr  = FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				tempr *= 1.0/(1+Sexp[ai]) - FF_2BODY[curr_pair_type_idx].OVRPRMS[4] * S[ai] * Sexp[ai] / Sexp2;
				tempr *= FF_2BODY[curr_pair_type_idx].OVRPRMS[2] * powrik * FF_2BODY[curr_pair_type_idx].OVRPRMS[3] * exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[2] * powrik) / rik;				

				dEover[ai].X -= tempr*RVEC.X/rik;
				dEover[ai].Y -= tempr*RVEC.Y/rik;
				dEover[ai].Z -= tempr*RVEC.Z/rik;						
					

				int akp = SYSTEM.PARENT[ak];

				dEover[akp].X += tempr*RVEC.X/rik;
				dEover[akp].Y += tempr*RVEC.Y/rik;
				dEover[akp].Z += tempr*RVEC.Z/rik;		
			
				SYSTEM.PRESSURE_XYZ -= tempr * rik;		
				SYSTEM.PRESSURE_TENSORS_XYZ.X -= tempr * (RVEC.X*RVEC.X)/rik;
				SYSTEM.PRESSURE_TENSORS_XYZ.Y -= tempr * (RVEC.Y*RVEC.Y)/rik;
				SYSTEM.PRESSURE_TENSORS_XYZ.Z -= tempr * (RVEC.Z*RVEC.Z)/rik;				
			}
		}
	}

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
		SYSTEM.ACCEL[a1].X -= dEover[a1].X;
		SYSTEM.ACCEL[a1].Y -= dEover[a1].Y;
		SYSTEM.ACCEL[a1].Z -= dEover[a1].Z;
	}
		
}


