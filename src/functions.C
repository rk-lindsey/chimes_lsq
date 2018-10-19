
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
#include "io_styles.h"
#include "A_Matrix.h"

#ifdef USE_MPI
	#include <mpi.h>
#endif

using namespace std;


//////////////////////////////////////////
//
//	SMALL UTILITY FUNCTION
//
//////////////////////////////////////////


//////////////////////////////////////////
// Functions to manage layers
//////////////////////////////////////////


void build_real_replicates(FRAME &SYSTEM, const JOB_CONTROL &CONTROLS)
// Create real (non-ghost) replicates of the input atoms.
{
  if(RANK == 0)
	 cout << "Building " << CONTROLS.REAL_REPLICATES << " replicates..." << endl;
	
  int TEMP_IDX = SYSTEM.ATOMS;
  XYZ TEMP_XYZ{0.0, 0.0, 0.0} ;
		
  SYSTEM.PARENT   .resize(SYSTEM.ATOMS);
  SYSTEM.LAYER_IDX.resize(SYSTEM.ATOMS);
			
  // Create coordinates for the layer atoms. layer elements do not include 0, 0, 0, which is the main cell

  for(int a1=0; a1<SYSTEM.ATOMS; a1++)
  {			
	 for(int n1=0; n1<=CONTROLS.REAL_REPLICATES; n1++)
	 {
		for(int n2=0; n2<=CONTROLS.REAL_REPLICATES; n2++)
		{
		  for(int n3=0; n3<=CONTROLS.REAL_REPLICATES; n3++)
		  {	
			 if ((n1 == 0) && (n2 == 0) && (n3 == 0) )
				SYSTEM.PARENT[a1] = a1;
			 else
			 {											
				TEMP_XYZ.X = SYSTEM.COORDS.at(a1).X + n1 * SYSTEM.BOXDIM.X;
				TEMP_XYZ.Y = SYSTEM.COORDS.at(a1).Y + n2 * SYSTEM.BOXDIM.Y;
				TEMP_XYZ.Z = SYSTEM.COORDS.at(a1).Z + n3 * SYSTEM.BOXDIM.Z;

				SYSTEM.COORDS       .push_back(TEMP_XYZ);
				SYSTEM.ATOMTYPE     .push_back(SYSTEM.ATOMTYPE     .at(a1));
				SYSTEM.ATOMTYPE_IDX .push_back(SYSTEM.ATOMTYPE_IDX .at(a1));
				SYSTEM.CHARGES      .push_back(SYSTEM.CHARGES      .at(a1));
				SYSTEM.MASS         .push_back(SYSTEM.MASS         .at(a1));	
				SYSTEM.VELOCITY     .push_back(SYSTEM.VELOCITY     .at(a1));
				SYSTEM.VELOCITY_ITER.push_back(SYSTEM.VELOCITY_ITER.at(a1));

				TEMP_IDX++;
							
				SYSTEM.PARENT.push_back(TEMP_IDX);
			 }
		  }
		}
	 }
  }
		
  SYSTEM.ATOMS = TEMP_IDX;

  SYSTEM.BOXDIM.X *= (CONTROLS.REAL_REPLICATES + 1);
  SYSTEM.BOXDIM.Y *= (CONTROLS.REAL_REPLICATES + 1);
  SYSTEM.BOXDIM.Z *= (CONTROLS.REAL_REPLICATES + 1);

  SYSTEM.FORCES      .resize(SYSTEM.ATOMS);
  SYSTEM.ACCEL       .resize(SYSTEM.ATOMS);
  SYSTEM.VELOCITY_NEW.resize(SYSTEM.ATOMS);
		
  for(int a1=0; a1<SYSTEM.ATOMS; a1++)
  {			
	 SYSTEM.ACCEL[a1].X = 0;
	 SYSTEM.ACCEL[a1].Y = 0;
	 SYSTEM.ACCEL[a1].Z = 0;
  }
		
  if(!CONTROLS.COMPARE_FORCE && !CONTROLS.SUBTRACT_FORCE)
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


//////////////////////////////////////////
// MPI compatibility functions
//////////////////////////////////////////
 
void divide_atoms(int &a1start, int &a1end, int atoms) 
{
	int procs_used;

	// Deal with no tasks to perform.
	if ( atoms <= 0 ) 
	{

	  a1start = 1 ;
	  a1end = 0 ;
	  return ;
	}

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
	
	REPLICATE.ATOMS		= SYSTEM.ATOMS;
	REPLICATE.ALL_ATOMS	= SYSTEM.ALL_ATOMS;

	a1start = 0;
	a2start = 0;
	a1end   = REPLICATE.ATOMS;
	a2end   = REPLICATE.ALL_ATOMS;

	REPLICATE.MY_ATOMS   = SYSTEM.MY_ATOMS ;
	REPLICATE.MY_ATOMS_START = SYSTEM.MY_ATOMS_START ;
	
	REPLICATE.BOXDIM.X	= SYSTEM.BOXDIM.X;
	REPLICATE.BOXDIM.Y	= SYSTEM.BOXDIM.Y;
	REPLICATE.BOXDIM.Z	= SYSTEM.BOXDIM.Z;

	REPLICATE.WRAPDIM.X	= SYSTEM.WRAPDIM.X;
	REPLICATE.WRAPDIM.Y	= SYSTEM.WRAPDIM.Y;
	REPLICATE.WRAPDIM.Z	= SYSTEM.WRAPDIM.Z;

	REPLICATE.TOT_POT_ENER 			= SYSTEM.TOT_POT_ENER;
	
	REPLICATE.FORCES      .resize(REPLICATE.ATOMS);
	REPLICATE.ACCEL       .resize(REPLICATE.ATOMS);
	REPLICATE.COORDS      .resize(REPLICATE.ATOMS);
	REPLICATE.WRAP_IDX    .resize(REPLICATE.ATOMS) ;
	REPLICATE.MASS        .resize(REPLICATE.ATOMS);

	for(int i=a1start; i<a1end; i++)
	{
		REPLICATE.FORCES[i] = SYSTEM.FORCES[i] ;
		REPLICATE.ACCEL[i]  = SYSTEM.ACCEL[i] ;

		REPLICATE.COORDS[i].X = SYSTEM.COORDS[i].X;
		REPLICATE.COORDS[i].Y = SYSTEM.COORDS[i].Y;
		REPLICATE.COORDS[i].Z = SYSTEM.COORDS[i].Z;		

		REPLICATE.WRAP_IDX[i].X = SYSTEM.WRAP_IDX[i].X ;
		REPLICATE.WRAP_IDX[i].Y = SYSTEM.WRAP_IDX[i].Y ;
		REPLICATE.WRAP_IDX[i].Z = SYSTEM.WRAP_IDX[i].Z ;

		REPLICATE.MASS        [i] = SYSTEM.MASS[i];

	}

	REPLICATE.ATOMTYPE    .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.ATOMTYPE_IDX.resize(REPLICATE.ALL_ATOMS);
	REPLICATE.PARENT      .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.ALL_COORDS  .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.CHARGES     .resize(REPLICATE.ALL_ATOMS);
	REPLICATE.LAYER_IDX   .resize(REPLICATE.ALL_ATOMS);

	for(int i=a2start; i<a2end; i++)
	{
		REPLICATE.ATOMTYPE    [i] = SYSTEM.ATOMTYPE[i];
		REPLICATE.ATOMTYPE_IDX[i] = SYSTEM.ATOMTYPE_IDX[i];
		REPLICATE.PARENT      [i] = SYSTEM.PARENT[i];
		REPLICATE.CHARGES     [i] = SYSTEM.CHARGES[i];

		REPLICATE.LAYER_IDX[i].X = SYSTEM.LAYER_IDX[i].X;
		REPLICATE.LAYER_IDX[i].Y = SYSTEM.LAYER_IDX[i].X;
		REPLICATE.LAYER_IDX[i].Z = SYSTEM.LAYER_IDX[i].X;		

		REPLICATE.ALL_COORDS[i].X = SYSTEM.ALL_COORDS[i].X;
		REPLICATE.ALL_COORDS[i].Y = SYSTEM.ALL_COORDS[i].Y;
		REPLICATE.ALL_COORDS[i].Z = SYSTEM.ALL_COORDS[i].Z;
	}
	
	REPLICATE.TOT_POT_ENER = SYSTEM.TOT_POT_ENER ;
	REPLICATE.PRESSURE     = SYSTEM.PRESSURE ;
	REPLICATE.PRESSURE_XYZ = SYSTEM.PRESSURE_XYZ ;
	REPLICATE.PRESSURE_TENSORS_XYZ = SYSTEM.PRESSURE_TENSORS_XYZ ;
	REPLICATE.PRESSURE_TENSORS = SYSTEM.PRESSURE_TENSORS ;

	REPLICATE.TEMPERATURE  = SYSTEM.TEMPERATURE ;
	REPLICATE.AVG_TEMPERATURE = SYSTEM.AVG_TEMPERATURE ;
	REPLICATE.QM_POT_ENER = SYSTEM.QM_POT_ENER ;

	REPLICATE.STRESS_TENSORS = SYSTEM.STRESS_TENSORS ;
	REPLICATE.STRESS_TENSORS_X = SYSTEM.STRESS_TENSORS_X ;
	REPLICATE.STRESS_TENSORS_Y = SYSTEM.STRESS_TENSORS_Y ;
	REPLICATE.STRESS_TENSORS_Z = SYSTEM.STRESS_TENSORS_Z ;
	
}

void numerical_pressure(const FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY,
								CLUSTER_LIST & TRIPS,  CLUSTER_LIST &QUADS, map<string,int> & PAIR_MAP,
								vector<int> &INT_PAIR_MAP,
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
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP,
			TRIPS, QUADS, NEIGHBOR_LIST);
	
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
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP,
			TRIPS, QUADS, NEIGHBOR_LIST);
	
	Vtot2 = REPLICATE.TOT_POT_ENER;

	// compute (return) pressure 
	
	//return -(Vtot2 - Vtot1)/(Vol2 - Vol1);
	
	PE_1 = Vtot1;
	PE_2 = Vtot2;
	dV   = Vol2 - Vol1;
}


void check_forces(FRAME& SYSTEM, JOB_CONTROL &CONTROLS, vector<PAIR_FF> &FF_2BODY, 
						map<string,int>& PAIR_MAP, vector<int> &INT_PAIR_MAP, 
						CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, NEIGHBORS &NEIGHBOR_LIST)
// Check the forces by finite derivative of the energy.  This will be computationally
// expensive, but an important check.
{
  vector<XYZ> coords(SYSTEM.ATOMS) ;
  vector<XYZ> forces(SYSTEM.ATOMS) ;
  const double eps = 1.0e-06 ;
  const double pass = 1.0e-04 ;


  for ( int a1 = 0 ; a1 < SYSTEM.ATOMS ; a1++ ) 
  {
	 coords[a1] = SYSTEM.COORDS[a1] ;
	 forces[a1] = SYSTEM.ACCEL[a1] ;
  }

  for ( int a1 = 0 ; a1 < SYSTEM.ATOMS ; a1++ ) {
	 for ( int j = 0 ; j < 3 ; j++ ) {

		// Change position of atom.
		if ( j == 0 ) 
		  SYSTEM.COORDS[a1].X += eps ;
		else if ( j == 1 ) 
		  SYSTEM.COORDS[a1].Y += eps ;
		else if ( j == 2 )
		  SYSTEM.COORDS[a1].Z += eps ;

		SYSTEM.update_ghost(CONTROLS.N_LAYERS) ;

		// Recalculate the forces
		NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);

		ZCalc(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);

		// FOR MPI:		Synchronize forces, energy, and pressure.
#ifdef USE_MPI
		sum_forces(SYSTEM.ACCEL, SYSTEM.ATOMS, SYSTEM.TOT_POT_ENER, SYSTEM.PRESSURE_XYZ, SYSTEM.PRESSURE_TENSORS_XYZ.X, SYSTEM.PRESSURE_TENSORS_XYZ.Y, SYSTEM.PRESSURE_TENSORS_XYZ.Z);
#endif

		double energy1 = SYSTEM.TOT_POT_ENER ;

		// Move the atom again.
		if ( j == 0 ) 
		  SYSTEM.COORDS[a1].X = coords[a1].X - eps ;
		else if ( j == 1 ) 
		  SYSTEM.COORDS[a1].Y = coords[a1].Y - eps ;
		else if ( j == 2 )
		  SYSTEM.COORDS[a1].Z = coords[a1].Z - eps ;

		SYSTEM.update_ghost(CONTROLS.N_LAYERS) ;

		NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);

		ZCalc(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);

		// FOR MPI:		Synchronize forces, energy, and pressure.
#ifdef USE_MPI
		sum_forces(SYSTEM.ACCEL, SYSTEM.ATOMS, SYSTEM.TOT_POT_ENER, SYSTEM.PRESSURE_XYZ, SYSTEM.PRESSURE_TENSORS_XYZ.X, SYSTEM.PRESSURE_TENSORS_XYZ.Y, SYSTEM.PRESSURE_TENSORS_XYZ.Z);
#endif

 		// Use symmetric difference for higher accuracy in numerical derivative.
		double fcheck = (SYSTEM.TOT_POT_ENER - energy1) / (2.0 * eps) ;

		double diff ;
		if ( j == 0 )
		  diff = fabs( (fcheck - forces[a1].X) / (1.0 + fabs(forces[a1].X) ) ) ;
		else if ( j == 1 ) 
		  diff = fabs( (fcheck - forces[a1].Y) / (1.0 + fabs(forces[a1].Y)) ) ;
		else if ( j == 2 ) 
		  diff = fabs( (fcheck - forces[a1].Z) / (1.0 + fabs(forces[a1].Z)) ) ;
		
		if ( RANK == 0 ) 
		{
      		if ( diff > pass ) 
		      {
		         cout << "Failed force check for atom " << a1 << " coordinate " << j ;
					cout << " Error = " << setprecision(6) << setw(10) << diff << endl ;
	         }
				else 
				{
		         cout << "Passed force check for atom " << a1 << " coordinate " << j << endl ;
	         }
       }

		SYSTEM.COORDS[a1].X = coords[a1].X ;
		SYSTEM.COORDS[a1].Y = coords[a1].Y ;
		SYSTEM.COORDS[a1].Z = coords[a1].Z ;
	 }
  }
  // Reset original values.
  for ( int a1 = 0 ; a1 < SYSTEM.ATOMS ; a1++ ) 
  {
	 SYSTEM.COORDS[a1] = coords[a1] ;
	 SYSTEM.ACCEL[a1] = forces[a1] ;
  }
  SYSTEM.update_ghost(CONTROLS.N_LAYERS) ;

}


//////////////////////////////////////////
// Bad configuration tracking functions
//////////////////////////////////////////
/*

void PRINT_CONFIG(FRAME &SYSTEM, JOB_CONTROL & CONTROLS, int type)
// Output the current configuration
{
	int PRINT_WIDTH     = 21; // Use 21 for testing
	int PRINT_PRECISION = 14; // Use 14 for testing
	
	if(type == 1)
	{
		ofstream &BAD_CONFIGS = BAD_CONFIGS_1;
	
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
	else if (type == 2)
	{
		ofstream &BAD_CONFIGS = BAD_CONFIGS_2;
	
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
	else
	{
		cout << "ERROR: Unknown print type. See PRINT_CONFIG()." << endl;
		exit(1);
	}	
	
}
*/

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- DERIVATIVE CALCULATION
//
//////////////////////////////////////////

static void ZCalc_Spline_Deriv  (JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, A_MAT & A_MATRIX,  map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

static void ZCalc_Poly_Deriv    (JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, A_MAT & A_MATRIX,  map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);	

static void ZCalc_InvR_Deriv    (JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, A_MAT & A_MATRIX,  map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- FORCE CALCULATION
//
//////////////////////////////////////////

static void ZCalc_Lj(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST) ;

static void ZCalc_Spline(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST) ;

static void ZCalcSR_Over(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

static void ZCalc_Cluster(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST) ;


////////////////////////////////////////////////////////////
//
// FUNCTIONS -- DERIVATIVE CALCULATION
//
////////////////////////////////////////////////////////////

void ZCalc_Deriv (JOB_CONTROL & CONTROLS, vector<PAIRS> & FF_2BODY,  CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, FRAME & FRAME_SYSTEM, A_MAT & A_MATRIX,
									map<string,int> & PAIR_MAP,  vector<int> &INT_PAIR_MAP, NEIGHBORS &NEIGHBOR_LIST)
// Controls which functions are used to calculate derivatives
{
	// Check for control option compatability:
	
  	vector<TRIPLETS> & PAIR_TRIPLETS = TRIPS.VEC ;
  
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
		ZCalc_Ewald_Deriv(FRAME_SYSTEM, FF_2BODY, A_MATRIX, PAIR_MAP, NEIGHBOR_LIST, CONTROLS);	
	
    	if ( FF_2BODY[0].PAIRTYP == "SPLINE" )
		 ZCalc_Spline_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, A_MATRIX, PAIR_MAP, NEIGHBOR_LIST);
	
	else if ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" )
	{
		// Only enter if 2B are requested. For example, skip if user wants to fit ONLY 3B cheby
		// i.e. PAIRTYP: CHEBYSHEV  0 6 or similar

	  Cheby cheby{CONTROLS,FRAME_SYSTEM,NEIGHBOR_LIST,FF_2BODY,INT_PAIR_MAP} ;

	  if ( FF_2BODY[0].SNUM > 0)
		 cheby.Deriv_2B(A_MATRIX) ;
	
	  if (CONTROLS.USE_3B_CHEBY)
		 cheby.Deriv_3B(A_MATRIX,TRIPS) ;
			
		if (CONTROLS.USE_4B_CHEBY) 
		{
			int n_3b_cheby_terms = 0 ;
			for (int i=0; i<PAIR_TRIPLETS.size(); i++) 
				n_3b_cheby_terms += PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS;

		 	cheby.Deriv_4B(A_MATRIX, n_3b_cheby_terms, QUADS) ;
		}
	}			

    	else if ( FF_2BODY[0].PAIRTYP == "DFTBPOLY" )	
		 ZCalc_Poly_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, A_MATRIX, PAIR_MAP, NEIGHBOR_LIST);

    	else if ( FF_2BODY[0].PAIRTYP == "INVRSE_R" )	
		 ZCalc_InvR_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, A_MATRIX, PAIR_MAP, NEIGHBOR_LIST);

    	else 
    {
		cout << "Error: bad pairtype in ZCalc_Deriv\n";
		exit(1);
    }		
}	

// FUNCTION UPDATED
static void ZCalc_Spline_Deriv(JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, A_MAT & A_MATRIX,  map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)		   
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
				
			   	A_MATRIX.FORCES[a1][vstart+kstart+0].X += h00 * RAB.X / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+1].X += h10 * RAB.X / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+2].X += h01 * RAB.X / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+3].X += h11 * RAB.X / rlen;
			  
			   	A_MATRIX.FORCES[a1][vstart+kstart+0].Y += h00 * RAB.Y / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+1].Y += h10 * RAB.Y / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+2].Y += h01 * RAB.Y / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+3].Y += h11 * RAB.Y / rlen;
				  
			   	A_MATRIX.FORCES[a1][vstart+kstart+0].Z += h00 * RAB.Z / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+1].Z += h10 * RAB.Z / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+2].Z += h01 * RAB.Z / rlen;
			   	A_MATRIX.FORCES[a1][vstart+kstart+3].Z += h11 * RAB.Z / rlen;
				
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+0].X -= h00 * RAB.X / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+1].X -= h10 * RAB.X / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+2].X -= h01 * RAB.X / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+3].X -= h11 * RAB.X / rlen;
			  
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+0].Y -= h00 * RAB.Y / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+1].Y -= h10 * RAB.Y / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+2].Y -= h01 * RAB.Y / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+3].Y -= h11 * RAB.Y / rlen;
				  
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+0].Z -= h00 * RAB.Z / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+1].Z -= h10 * RAB.Z / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+2].Z -= h01 * RAB.Z / rlen;
			   	A_MATRIX.FORCES[fidx_a2][vstart+kstart+3].Z -= h11 * RAB.Z / rlen;
				
				if (CONTROLS.FIT_STRESS)
				{
				   	A_MATRIX.STRESSES[vstart+kstart+0].XX -= h00 * RAB.X * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].XX -= h10 * RAB.X * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].XX -= h01 * RAB.X * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].XX -= h11 * RAB.X * RAB.X / rlen;
			  
				   	A_MATRIX.STRESSES[vstart+kstart+0].YY -= h00 * RAB.Y * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].YY -= h10 * RAB.Y * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].YY -= h01 * RAB.Y * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].YY -= h11 * RAB.Y * RAB.Y / rlen;
				  
				   	A_MATRIX.STRESSES[vstart+kstart+0].ZZ -= h00 * RAB.Z * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].ZZ -= h10 * RAB.Z * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].ZZ -= h01 * RAB.Z * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].ZZ -= h11 * RAB.Z * RAB.Z / rlen;					  
				}
				else if (CONTROLS.FIT_STRESS_ALL)
				{
					// xx, xy, and xz
					
				   	A_MATRIX.STRESSES[vstart+kstart+0].XX -= h00 * RAB.X * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].XX -= h10 * RAB.X * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].XX -= h01 * RAB.X * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].XX -= h11 * RAB.X * RAB.X / rlen;
			  
				   	A_MATRIX.STRESSES[vstart+kstart+0].XY -= h00 * RAB.X * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].XY -= h10 * RAB.X * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].XY -= h01 * RAB.X * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].XY -= h11 * RAB.X * RAB.Y / rlen;
				  
				   	A_MATRIX.STRESSES[vstart+kstart+0].XZ -= h00 * RAB.X * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].XZ -= h10 * RAB.X * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].XZ -= h01 * RAB.X * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].XZ -= h11 * RAB.X * RAB.Z / rlen;	
					
					
					// yx, yy, and yz
					
				   	A_MATRIX.STRESSES[vstart+kstart+0].XY -= h00 * RAB.Y * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].XY -= h10 * RAB.Y * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].XY -= h01 * RAB.Y * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].XY -= h11 * RAB.Y * RAB.X / rlen;
			  
				   	A_MATRIX.STRESSES[vstart+kstart+0].YY -= h00 * RAB.Y * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].YY -= h10 * RAB.Y * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].YY -= h01 * RAB.Y * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].YY -= h11 * RAB.Y * RAB.Y / rlen;
				  
				   	A_MATRIX.STRESSES[vstart+kstart+0].YZ -= h00 * RAB.Y * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].YZ -= h10 * RAB.Y * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].YZ -= h01 * RAB.Y * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].YZ -= h11 * RAB.Y * RAB.Z / rlen;	
					
					
					// yx, yy, and yz
					
				   	A_MATRIX.STRESSES[vstart+kstart+0].XZ -= h00 * RAB.Z * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].XZ -= h10 * RAB.Z * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].XZ -= h01 * RAB.Z * RAB.X / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].XZ -= h11 * RAB.Z * RAB.X / rlen;
			  
				   	A_MATRIX.STRESSES[vstart+kstart+0].YZ -= h00 * RAB.Z * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].YZ -= h10 * RAB.Z * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].YZ -= h01 * RAB.Z * RAB.Y / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].YZ -= h11 * RAB.Z * RAB.Y / rlen;
				  
				   	A_MATRIX.STRESSES[vstart+kstart+0].ZZ -= h00 * RAB.Z * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+1].ZZ -= h10 * RAB.Z * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+2].ZZ -= h01 * RAB.Z * RAB.Z / rlen;
				   	A_MATRIX.STRESSES[vstart+kstart+3].ZZ -= h11 * RAB.Z * RAB.Z / rlen;	
				}							
				
			}
		}
	}
	
	if (CONTROLS.FIT_STRESS)
	{
		for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
		{
			 A_MATRIX.STRESSES[i].XX /= VOL;
			 A_MATRIX.STRESSES[i].YY /= VOL;
			 A_MATRIX.STRESSES[i].ZZ /= VOL;     
		}
	}
	
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
		{
			A_MATRIX.STRESSES[i].XX /= VOL;
			A_MATRIX.STRESSES[i].XY /= VOL;
			A_MATRIX.STRESSES[i].XZ /= VOL;      
		
			A_MATRIX.STRESSES[i].XY /= VOL;
			A_MATRIX.STRESSES[i].YY /= VOL;
			A_MATRIX.STRESSES[i].YZ /= VOL;
		
			A_MATRIX.STRESSES[i].XZ /= VOL;
			A_MATRIX.STRESSES[i].YZ /= VOL;
			A_MATRIX.STRESSES[i].ZZ /= VOL;
		}
	}
	
	return;
}

// FUNCTION UPDATED
static void ZCalc_InvR_Deriv(JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, A_MAT & A_MATRIX,  map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)		
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

					A_MATRIX.FORCES[a1][vstart+i].X += rfac * RAB.X/rlen;
					A_MATRIX.FORCES[a1][vstart+i].Y += rfac * RAB.Y/rlen;
					A_MATRIX.FORCES[a1][vstart+i].Z += rfac * RAB.Z/rlen;
						
					A_MATRIX.FORCES[fidx_a2][vstart+i].X -= rfac * RAB.X/rlen;
					A_MATRIX.FORCES[fidx_a2][vstart+i].Y -= rfac * RAB.Y/rlen;
					A_MATRIX.FORCES[fidx_a2][vstart+i].Z -= rfac * RAB.Z/rlen;		
					
					if (CONTROLS.FIT_STRESS)
					{
						A_MATRIX.STRESSES[vstart+i].XX -= rfac * RAB.X * RAB.X / rlen;
						A_MATRIX.STRESSES[vstart+i].YY -= rfac * RAB.Y * RAB.Y / rlen;
						A_MATRIX.STRESSES[vstart+i].ZZ -= rfac * RAB.Z * RAB.Z / rlen;	     
					}
					
					else if (CONTROLS.FIT_STRESS_ALL)
					{
						A_MATRIX.STRESSES[vstart+i].XX -= rfac * RAB.X * RAB.X / rlen;
						A_MATRIX.STRESSES[vstart+i].XY -= rfac * RAB.X * RAB.Y / rlen;
						A_MATRIX.STRESSES[vstart+i].XZ -= rfac * RAB.X * RAB.Z / rlen;	   
						
						A_MATRIX.STRESSES[vstart+i].XY -= rfac * RAB.Y * RAB.X / rlen;
						A_MATRIX.STRESSES[vstart+i].YY -= rfac * RAB.Y * RAB.Y / rlen;
						A_MATRIX.STRESSES[vstart+i].YZ -= rfac * RAB.Y * RAB.Z / rlen;	   
						
						A_MATRIX.STRESSES[vstart+i].XZ -= rfac * RAB.Z * RAB.X / rlen;
						A_MATRIX.STRESSES[vstart+i].YZ -= rfac * RAB.Z * RAB.Y / rlen;
						A_MATRIX.STRESSES[vstart+i].ZZ -= rfac * RAB.Z * RAB.Z / rlen;	   
					}						
				}
			}
		}
	}
	
	if (CONTROLS.FIT_STRESS)
	{
		for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
		{
			A_MATRIX.STRESSES[i].XX /= VOL;
			A_MATRIX.STRESSES[i].YY /= VOL;
			A_MATRIX.STRESSES[i].ZZ /= VOL;        
		}
	}
	
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
		{
			A_MATRIX.STRESSES[i].XX /= VOL;
			A_MATRIX.STRESSES[i].XY /= VOL;
			A_MATRIX.STRESSES[i].XZ /= VOL;      
		
			A_MATRIX.STRESSES[i].XY /= VOL;
			A_MATRIX.STRESSES[i].YY /= VOL;
			A_MATRIX.STRESSES[i].YZ /= VOL;
		
			A_MATRIX.STRESSES[i].XZ /= VOL;
			A_MATRIX.STRESSES[i].YZ /= VOL;
			A_MATRIX.STRESSES[i].ZZ /= VOL;
		}
	}
	
	return;

}

// FUNCTION UPDATED
static void ZCalc_Poly_Deriv(JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, A_MAT & A_MATRIX, map<string,int> PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)	
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

					A_MATRIX.FORCES[a1][vstart+i].X += rfac * RAB.X/rlen;
					A_MATRIX.FORCES[a1][vstart+i].Y += rfac * RAB.Y/rlen;
					A_MATRIX.FORCES[a1][vstart+i].Z += rfac * RAB.Z/rlen;
						
					A_MATRIX.FORCES[fidx_a2][vstart+i].X -= rfac * RAB.X/rlen;
					A_MATRIX.FORCES[fidx_a2][vstart+i].Y -= rfac * RAB.Y/rlen;
					A_MATRIX.FORCES[fidx_a2][vstart+i].Z -= rfac * RAB.Z/rlen;
					
					if (CONTROLS.FIT_STRESS)
					{
						A_MATRIX.STRESSES[vstart+i].XX -= 29421.02407027691 * rfac * RAB.X * RAB.X / rlen;
						A_MATRIX.STRESSES[vstart+i].YY -= 29421.02407027691 * rfac * RAB.Y * RAB.Y / rlen;
						A_MATRIX.STRESSES[vstart+i].ZZ -= 29421.02407027691 * rfac * RAB.Z * RAB.Z / rlen;    
					}
					
					else if (CONTROLS.FIT_STRESS_ALL)
					{
						A_MATRIX.STRESSES[vstart+i].XX -= 29421.02407027691 * rfac * RAB.X * RAB.X / rlen;
						A_MATRIX.STRESSES[vstart+i].XY -= 29421.02407027691 * rfac * RAB.X * RAB.Y / rlen;
						A_MATRIX.STRESSES[vstart+i].XZ -= 29421.02407027691 * rfac * RAB.X * RAB.Z / rlen;  
						
						A_MATRIX.STRESSES[vstart+i].XY -= 29421.02407027691 * rfac * RAB.Y * RAB.X / rlen;
						A_MATRIX.STRESSES[vstart+i].YY -= 29421.02407027691 * rfac * RAB.Y * RAB.Y / rlen;
						A_MATRIX.STRESSES[vstart+i].YZ -= 29421.02407027691 * rfac * RAB.Y * RAB.Z / rlen;  
						
						A_MATRIX.STRESSES[vstart+i].XZ -= 29421.02407027691 * rfac * RAB.Z * RAB.X / rlen;
						A_MATRIX.STRESSES[vstart+i].YZ -= 29421.02407027691 * rfac * RAB.Z * RAB.Y / rlen;
						A_MATRIX.STRESSES[vstart+i].ZZ -= 29421.02407027691 * rfac * RAB.Z * RAB.Z / rlen;  
					}
				}
			}

		}
	}
	
	if (CONTROLS.FIT_STRESS)
	{
		for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
		{
			A_MATRIX.STRESSES[i].XX /= VOL;
			A_MATRIX.STRESSES[i].YY /= VOL;
			A_MATRIX.STRESSES[i].ZZ /= VOL;        
		}
	}
	
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
		{
			A_MATRIX.STRESSES[i].XX /= VOL;
			A_MATRIX.STRESSES[i].XY /= VOL;
			A_MATRIX.STRESSES[i].XZ /= VOL;      
		
			A_MATRIX.STRESSES[i].XY /= VOL;
			A_MATRIX.STRESSES[i].YY /= VOL;
			A_MATRIX.STRESSES[i].YZ /= VOL;
		
			A_MATRIX.STRESSES[i].XZ /= VOL;
			A_MATRIX.STRESSES[i].YZ /= VOL;
			A_MATRIX.STRESSES[i].ZZ /= VOL;
		}
	}
	
	
	return;

}

// FUNCTION UPDATED
void SubtractCoordForces(FRAME & SYSTEM, bool calc_deriv, A_MAT & A_MATRIX, vector<PAIRS> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST, bool lsq_mode)
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
	        	A_MATRIX.OVERBONDING[a1].X -= dFover[a1].X;
			A_MATRIX.OVERBONDING[a1].Y -= dFover[a1].Y;
			A_MATRIX.OVERBONDING[a1].Z -= dFover[a1].Z;		
		}
	} // VERIFIED FROM HERE UP
}


////////////////////////////////////////////////////////////
//
// FUNCTIONS -- FORCE CALCULATION
//
////////////////////////////////////////////////////////////
 

void ZCalc(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int>& INT_PAIR_MAP,CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, NEIGHBORS & NEIGHBOR_LIST)
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

	if ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" ) 
	{
	  Cheby cheby{CONTROLS, SYSTEM, NEIGHBOR_LIST, FF_2BODY, INT_PAIR_MAP};
	  cheby.Force_all(TRIPS, QUADS);
	}
	else if ( FF_2BODY[0].PAIRTYP == "LJ" ) 
	  ZCalc_Lj(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST);
	
	else if ( FF_2BODY[0].PAIRTYP == "SPLINE" ) 
	  ZCalc_Spline(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST);	
	
	else if ( FF_2BODY[0].PAIRTYP == "CLUSTER" )
	  	  ZCalc_Cluster(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST);	
	else 
    {
		cout << "Error: bad pairtype in ZCalc: " << FF_2BODY[0].PAIRTYP << endl;
		exit_run(1);
    }	
	
	if ( CONTROLS.USE_COULOMB ) 
		ZCalc_Ewald(SYSTEM, CONTROLS, NEIGHBOR_LIST);
		
	if ( CONTROLS.USE_OVERCOORD ) 	
	  ZCalcSR_Over(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST);

	// FUNCTIONS THAT NEED UPDATING:

	/* 

	else if ( pair_type == INVERSE_R ) 	
		ZCalc_SR_Analytic(Coord,Lbc, Latcons,nat,smin,smax,snum, SForce,Vtot, Pxyz, params);

	else if ( pair_type == STILLINGER )  // SAVE THIS FOR SECOND TO LAST FOR SIMILAR REASONS  
		ZCalc_Stillinger(Coord,Lbc, Latcons,nat,smax, SForce,Vtot,Pxyz);

	else 
		EXIT_MSG("Error: Unknown pair type", pair_type)
	*/
	
	SYSTEM.PRESSURE_XYZ           /= 3.0 *  SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;
	SYSTEM.PRESSURE_TENSORS_XYZ.X /=        SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;
	SYSTEM.PRESSURE_TENSORS_XYZ.Y /=        SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;
	SYSTEM.PRESSURE_TENSORS_XYZ.Z /=        SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;

  return;
}


static void ZCalc_Lj(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
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



static void ZCalc_Cluster(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
// Many-body Cluster interaction model used for testing purposes.
{
  XYZ	RVEC ;
  double	rlen_mi;
  double	fac;
  string	TEMP_STR;
	
  // Set up for MPI
	
  int a1start, a1end;
  int a2start, a2end;
  int fidx_a2;
  int curr_pair_type_idx ;
  double neighbors = 0.0 ;

#ifndef LINK_LAMMPS
  divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.
#else
  a1start = SYSTEM.MY_ATOMS_START;
  a1end   = SYSTEM.MY_ATOMS_START+SYSTEM.MY_ATOMS-1;
#endif

  double rcutoff2 = 2.0 ;   // Cutoff of the many-body interaction.
  double poten_pow = 3.0 ;  // Power used in calculating the many-body interaction.
  double cluster_rmin = 0.5 ;  // Prevent large cluster forces as small distances.
  
  for(int a1= a1start;a1 <= a1end; a1++ ) 
  {
	 double poten = 1.0 ;

	 a2start = 0;
	 a2end   = a1end ;

	 neighbors = 0.0 ;		
	 for(int a2idx =a2start; a2idx <= a2end ; a2idx++)
	 {			
		int a2 = a2idx ;
		
		if ( a2 == a1 )
		  continue ;

#ifndef LINK_LAMMPS
		curr_pair_type_idx =  INT_PAIR_MAP[SYSTEM.ATOMTYPE_IDX[a1]*CONTROLS.NATMTYP + SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]];
#else
		curr_pair_type_idx =  INT_PAIR_MAP[(SYSTEM.ATOMTYPE_IDX[a1]-1)*CONTROLS.NATMTYP + (SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]-1)];
#endif
		
		// pair interaction cutoff distance.
		double rcutoff = FF_2BODY[curr_pair_type_idx].S_MAXIM;
		double eps = 0.05 ;
		
		rlen_mi = get_dist(SYSTEM, RVEC, a1, a2);

		if ( rlen_mi < rcutoff )
		{
		  // Lennard-Jones Interaction with cubic cutoff.
		  double poten2 = 4.0 * eps * ( pow(1.0/rlen_mi,12.0) - pow(1.0/rlen_mi,6.0) )
			 * pow(rcutoff - rlen_mi,3.0) ;
			 //* pow(rcutoff-rlen_mi, 3.0) ;

		  SYSTEM.TOT_POT_ENER += poten2 ;
		  fac = 4.0 * eps * ( -12.0 * pow(1.0/rlen_mi,14.0) + 6.0 * pow(1.0/rlen_mi,8.0) )
			 * pow(rcutoff - rlen_mi,3.0) 
			 - 3.0 * poten2 / ((rcutoff - rlen_mi) * rlen_mi); // * pow(rcutoff-rlen_mi, 2.0) / rlen_mi ;
		  
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

		// Many-body cluster term.  Includes a regularization for small r (cluster_rmin)
		// and a force cutoff at rcutoff2.
		if ( rlen_mi < rcutoff2 )
		{
		  poten *= pow(1.0/(rlen_mi+cluster_rmin), poten_pow) * pow(rcutoff2-rlen_mi, 3.0) ;
		}
	 }

	 poten = -160.0 * poten ;
	 SYSTEM.TOT_POT_ENER += poten ;
	 
	 for(int a2idx =a2start; a2idx <= a2end; a2idx++)
	 {			
		int a2 = a2idx ;

		if ( a2 == a1 )
		  continue ;
		
#ifndef LINK_LAMMPS
		curr_pair_type_idx =  INT_PAIR_MAP[SYSTEM.ATOMTYPE_IDX[a1]*CONTROLS.NATMTYP + SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]];
#else
		curr_pair_type_idx =  INT_PAIR_MAP[(SYSTEM.ATOMTYPE_IDX[a1]-1)*CONTROLS.NATMTYP + (SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]-1)];
#endif


		rlen_mi = get_dist(SYSTEM, RVEC, a1, a2);

		cout << "Distance between atom " << a1 << " and " << a2 << " = " << rlen_mi << endl ;

 		if ( rlen_mi < rcutoff2 )
 		{
 		  double term = -poten_pow / ( (rlen_mi+cluster_rmin) * rlen_mi ) -
			 3.0 / ( (rcutoff2 - rlen_mi) * rlen_mi ) ;
		  
 		  fac = poten * term ;
		  
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

		  //cout << "Found neighbor between atom " << a1 << " and " << a2 << endl ;
		  neighbors += 1.0 ;
		}	
	 }
	 cout << "Coordination for atom " << a1 << " = " << neighbors << endl ;
  }
}

static void ZCalc_Spline(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
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
static void ZCalcSR_Over(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
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

void FRAME::SET_NATOMS_OF_TYPE() // setting NATOMS_OF_TYPE
{
	int NO_ATOM_TYPES = 0;
	vector<string> ATOM_TYPES;

	vector<string>::iterator it;

	for (int i=0; i<ATOMS; i++)
	{
		// Get the location of the current atom type in the "ATOM_TYPES" list
		// ...if it doesn't exist, add it
		
		it = find(ATOM_TYPES.begin(), ATOM_TYPES.end(), ATOMTYPE[i]);
		
		if (it == ATOM_TYPES.end()) // Then the atom type hasn't been added yet
		{
			NO_ATOM_TYPES++;
			ATOM_TYPES.push_back(ATOMTYPE[i]);
			NATOMS_OF_TYPE.push_back(1);
		}
		else
		{
			NATOMS_OF_TYPE[distance(ATOM_TYPES.begin(), it)]++;
		}
	}
}

void check_charges(FRAME &SYSTEM, vector<double>& TMP_CHARGES, const vector<string>& TMP_ATOMTYPE, vector<PAIR_FF> &FF_2BODY, int NATMTYP)
// Check the charges and adjust values to enforce charge neutrality if necessary.
{
  double eps_charge = 1.0e-03 ;

  double total_charge = 0.0 ;
  for(int a=0; a<SYSTEM.ATOMS;a++)
	 total_charge += SYSTEM.CHARGES[a] ;

  if ( total_charge > eps_charge * SYSTEM.ATOMS ) 
	 EXIT_MSG("System total charge is too large: " + to_string(total_charge) ) ;

  if ( RANK == 0 ) 
	 cout << "      Total system charge = " << total_charge << endl ;

  // // Enforce exact charge neutrality.
  vector<int> count(SYSTEM.ATOMS, 0) ;
  for(int a=0; a<SYSTEM.ATOMS;a++)
  {
	 SYSTEM.CHARGES[a]  -= total_charge / SYSTEM.ATOMS ;
	 count[SYSTEM.ATOMTYPE_IDX[a]]++ ;
  }

  if ( RANK == 0 ) 
	 cout << "      Charges modified for exact neutrality (e):\n" ;
  for ( int i = 0 ; i < NATMTYP ; i++ ) 
  {
	 TMP_CHARGES[i] -= (total_charge * count[i]) / SYSTEM.ATOMS ;
	 if ( RANK == 0 ) cout << "       " << SYSTEM.ATOMTYPE[i] << " " << fixed << setprecision(9) << TMP_CHARGES[i] << endl ;
#if defined(USE_MPI) && defined(LINK_LAMMPS)
	 LMP_CHARGE[i] = TMP_CHARGES[i]; // save charges to global variable for LAMMPS
#endif

  }
  if ( RANK == 0 )
	 cout << endl ;
  
  if ( RANK == 0 ) 
  {
	 cout << "     \n" ;
	 cout << "     Pair charges modified for exact neutrality:\n" ;
  }

  int NO_PAIRS = FF_2BODY.size() ;

  for(int i=0; i<NO_PAIRS; i++)
  {
	 for(int j=0; j<NATMTYP; j++)
	 {
		for(int k=0; k<NATMTYP; k++)
		{
		  if( FF_2BODY[i].ATM1TYP == TMP_ATOMTYPE[j] && FF_2BODY[i].ATM2TYP == TMP_ATOMTYPE[k] ) {
			 FF_2BODY[i].PAIR_CHRG = TMP_CHARGES[j] * TMP_CHARGES[k] ;
			 if ( RANK == 0 ) 
				cout << "       " << FF_2BODY[i].PRPR_NM << " " << FF_2BODY[i].PAIR_CHRG * ke << endl ;
			 break ;
		  }
		}
	 }
  }
}


void build_int_pair_map(int natmtyp, const vector<string> &atomtype, 
								const vector<int> &atomtype_idx,
								map<string,int> &pair_map, vector<int> &int_pair_map)
// Build an integer-valued mapping between the int map index and values in the pair map.
// This determines which element in ATOM_PAIRS or FF_2BODY a particular pair of atom
// types belongs to.
{
  int_pair_map.resize(natmtyp*natmtyp, -1);

  for(int i=0; i<natmtyp; i++)
  {
	 for (int j=0; j<natmtyp; j++)
	 {
		string int_map_str =      atomtype[i];
		int_map_str.append(atomtype[j]);
		
		int int_map_idx = atomtype_idx[i]*natmtyp + atomtype_idx[j];
		
		int_pair_map[int_map_idx] = pair_map[int_map_str];
		
		if(RANK == 0)
		{
		  // Save cout format state.
		  ofstream init ;
		  init.copyfmt(std::cout) ;

		  cout << "		";
		  cout<< "Atom type idxs: ";
		  cout<< fixed << setw(2) << right << i;
		  cout<< fixed << setw(2) << right << j;
		  cout<< " Pair name: "           << setw(4) << right << int_map_str;
		  cout<< " Explicit pair index: " << setw(4) << right << int_map_idx;
		  cout<< " Unique pair index: "   << setw(4) << right << int_pair_map[int_map_idx] << endl;		

		  // Restore cout format state.
		  std::cout.copyfmt(init) ;
		}

	 }
  }
}

void parse_fcut_input(string line, vector<PAIR_FF>& FF_2BODY, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS)
// Parse the input for the force cutoff for all interaction types.
{
	// Each 2-body interaction needs to be done separately.
	for(int i=0; i<FF_2BODY.size(); i++) 
	{
		FF_2BODY[i].FORCE_CUTOFF.parse_input(line);
		FF_2BODY[i].FORCE_CUTOFF.BODIEDNESS = 2;		
	}
		
	// Handle the many-body part
	//if ( FF_2BODY[0].SNUM_3B_CHEBY>0 ) 
	if ( TRIPS.VEC.size()>0 ) 
		TRIPS.parse_fcut(line) ;

	//if(FF_2BODY[0].SNUM_4B_CHEBY>0)
	if ( QUADS.VEC.size()>0 ) 
		QUADS.parse_fcut(line) ;
}

// MPI -- Related functions -- See headers at top of file

#ifdef USE_MPI
	void sum_forces(vector<XYZ>& accel_vec, int atoms, double &pot_energy, double &pressure, double & tens_x, double & tens_y, double & tens_z)
	// Add up forces, potential energy, and pressure from all processes.  
	{
		// Sum up the potential energy, pressure, tensors , and forces from all processors
	#ifdef USE_MPI	
		double *accel = (double *) accel_vec .data();

		MPI_Allreduce(MPI_IN_PLACE, &pot_energy,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &pressure,  1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &tens_x,    1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &tens_y,    1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, &tens_z,    1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
		MPI_Allreduce(MPI_IN_PLACE, accel, 3*atoms, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	#endif	

	}

	void sync_position(vector<XYZ>& coord_vec, NEIGHBORS & neigh_list, vector<XYZ>& velocity_vec, int atoms, bool sync_vel) 
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

#endif // USE_MPI



