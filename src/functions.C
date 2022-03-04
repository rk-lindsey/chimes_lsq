
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

#include "../imports/chimes_calculator/serial_interface/src/serial_chimes_interface.h"

#ifdef USE_MPI
	#include <mpi.h>
#endif

using namespace std;

static void ZCalc_Serial_Chimes(FRAME &SYSTEM, PAIR_FF &FF_2BODY) ;

//////////////////////////////////////////
//
//	SMALL UTILITY FUNCTION
//
//////////////////////////////////////////

void OPEN_TRAJFILE(ifstream & TRAJ_INPUT, vector<string> & INFILE, int FILE_IDX)
{
	if (TRAJ_INPUT.is_open())
		TRAJ_INPUT.close();
		
	if ( (FILE_IDX+1) > INFILE.size())
		EXIT_MSG("(FILE_IDX+1) > INFILE.size(): ", FILE_IDX) ;	

	TRAJ_INPUT.open(INFILE[FILE_IDX].data());
	
	if(!TRAJ_INPUT.is_open())
		EXIT_MSG("ERROR: Cannot open trajectory file: ", INFILE[FILE_IDX]) ;
}

double VECTOR_MAGNITUDE(vector<double> & vec)
{
	double MAG = 0;
	
	for (int i=0; i<vec.size(); i++)
		MAG += vec[i]*vec[i];
	
	return sqrt(MAG);
	
}
double VECTOR_ANGLE(vector<double> & v1, vector<double> & v2)
{
	double ANG = 0;
	
	if (v1.size() != v2.size())
		EXIT_MSG("ERROR: Cannot compute angle between vectors of different dimension");
		
	for (int i=0; i<v1.size(); i++)
		ANG += v1[i]*v2[i];
	
	ANG /= VECTOR_MAGNITUDE(v1);
	ANG /= VECTOR_MAGNITUDE(v2);
	
	return acos(ANG);
			
}

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
	XYZ_INT TMP_LAYER;
		
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
						TMP_LAYER.X = n1;
						TMP_LAYER.X = n2;
						TMP_LAYER.X = n3;
					
						SYSTEM.BOXDIM.LAYER_ATOM(SYSTEM.COORDS[a1], TMP_LAYER, TEMP_XYZ);

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

	SYSTEM.BOXDIM.CELL_AX *= (CONTROLS.REAL_REPLICATES + 1);
	SYSTEM.BOXDIM.CELL_AY *= (CONTROLS.REAL_REPLICATES + 1);
	SYSTEM.BOXDIM.CELL_AZ *= (CONTROLS.REAL_REPLICATES + 1);
  
	SYSTEM.BOXDIM.CELL_BX *= (CONTROLS.REAL_REPLICATES + 1);
	SYSTEM.BOXDIM.CELL_BY *= (CONTROLS.REAL_REPLICATES + 1);
	SYSTEM.BOXDIM.CELL_BZ *= (CONTROLS.REAL_REPLICATES + 1);
  
	SYSTEM.BOXDIM.CELL_CX *= (CONTROLS.REAL_REPLICATES + 1);
	SYSTEM.BOXDIM.CELL_CY *= (CONTROLS.REAL_REPLICATES + 1);
	SYSTEM.BOXDIM.CELL_CZ *= (CONTROLS.REAL_REPLICATES + 1);
  
	SYSTEM.BOXDIM.UPDATE_CELL();

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
		cout << "	New box dimensions: " << endl;
		cout <<  "	cell vectors (a)            " << SYSTEM.BOXDIM.CELL_AX << " " << SYSTEM.BOXDIM.CELL_AY << " " << SYSTEM.BOXDIM.CELL_AZ << endl;
		cout <<  "	cell vectors (b)            " << SYSTEM.BOXDIM.CELL_BX << " " << SYSTEM.BOXDIM.CELL_BY << " " << SYSTEM.BOXDIM.CELL_BZ << endl;
		cout <<  "	cell vectors (c)            " << SYSTEM.BOXDIM.CELL_CX << " " << SYSTEM.BOXDIM.CELL_CY << " " << SYSTEM.BOXDIM.CELL_CZ << endl;
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
	} 
	else if ( RANK == procs_used - 1 ) 
	{
		// End of the list.
		a1end = atoms - 1;
	} 
	else 
	{
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
	else if(TYPE == "ITER")
	{
		 for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		 {
				// Don't account for frozen atoms
		
				if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))
					 continue;
			
				Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY_ITER[a1].X * SYSTEM.VELOCITY_ITER[a1].X;
				Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY_ITER[a1].Y * SYSTEM.VELOCITY_ITER[a1].Y;
				Ktot += 0.5 * SYSTEM.MASS[a1] * SYSTEM.VELOCITY_ITER[a1].Z * SYSTEM.VELOCITY_ITER[a1].Z;
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

	REPLICATE.ATOMS		= SYSTEM.ATOMS;
	REPLICATE.ALL_ATOMS	= SYSTEM.ALL_ATOMS;

	a1start = 0;
	a2start = 0;
	a1end   = REPLICATE.ATOMS;
	a2end   = REPLICATE.ALL_ATOMS;
	
	REPLICATE.BOXDIM = SYSTEM.BOXDIM;

	REPLICATE.TOT_POT_ENER = SYSTEM.TOT_POT_ENER;
	
	REPLICATE.FORCES      .resize(REPLICATE.ATOMS);
	REPLICATE.ACCEL       .resize(REPLICATE.ATOMS);
	REPLICATE.COORDS      .resize(REPLICATE.ATOMS);
	REPLICATE.WRAP_IDX    .resize(REPLICATE.ATOMS) ;
	REPLICATE.MASS        .resize(REPLICATE.ATOMS);

	// Vector copy.
	REPLICATE.QM_ENERGY_OFFSET = SYSTEM.QM_ENERGY_OFFSET ;
	
	for(int i=a1start; i<a1end; i++)
	{
		REPLICATE.FORCES[i] = SYSTEM.FORCES[i] ;
		REPLICATE.ACCEL[i]  = SYSTEM.ACCEL[i] ;

		REPLICATE.COORDS[i].X = SYSTEM.COORDS[i].X;
		REPLICATE.COORDS[i].Y = SYSTEM.COORDS[i].Y;
		REPLICATE.COORDS[i].Z = SYSTEM.COORDS[i].Z;		

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

		REPLICATE.ALL_COORDS[i].X = SYSTEM.ALL_COORDS[i].X;
		REPLICATE.ALL_COORDS[i].Y = SYSTEM.ALL_COORDS[i].Y;
		REPLICATE.ALL_COORDS[i].Z = SYSTEM.ALL_COORDS[i].Z;
	}
	
	REPLICATE.TOT_POT_ENER 		= SYSTEM.TOT_POT_ENER;
	REPLICATE.PRESSURE     		= SYSTEM.PRESSURE;
	REPLICATE.PRESSURE_XYZ 		= SYSTEM.PRESSURE_XYZ;
	
	REPLICATE.PRESSURE_TENSORS_XYZ_ALL.resize(3);
	REPLICATE.PRESSURE_TENSORS_ALL.resize(3);

	for (int i=0; i<3; i++)
	{
		REPLICATE.PRESSURE_TENSORS_XYZ_ALL[i].X = SYSTEM.PRESSURE_TENSORS_XYZ_ALL[i].X;
		REPLICATE.PRESSURE_TENSORS_XYZ_ALL[i].Y = SYSTEM.PRESSURE_TENSORS_XYZ_ALL[i].Y;
		REPLICATE.PRESSURE_TENSORS_XYZ_ALL[i].Z = SYSTEM.PRESSURE_TENSORS_XYZ_ALL[i].Z;
		
		REPLICATE.PRESSURE_TENSORS_ALL[i].X     = SYSTEM.PRESSURE_TENSORS_ALL[i].X; 
		REPLICATE.PRESSURE_TENSORS_ALL[i].Y     = SYSTEM.PRESSURE_TENSORS_ALL[i].Y; 
		REPLICATE.PRESSURE_TENSORS_ALL[i].Z     = SYSTEM.PRESSURE_TENSORS_ALL[i].Z; 
	}
	REPLICATE.TEMPERATURE  		= SYSTEM.TEMPERATURE;
	REPLICATE.AVG_TEMPERATURE 	= SYSTEM.AVG_TEMPERATURE;
	REPLICATE.QM_POT_ENER 		= SYSTEM.QM_POT_ENER;

	REPLICATE.STRESS_TENSORS 	= SYSTEM.STRESS_TENSORS;
	REPLICATE.STRESS_TENSORS_X 	= SYSTEM.STRESS_TENSORS_X;
	REPLICATE.STRESS_TENSORS_Y 	= SYSTEM.STRESS_TENSORS_Y;
	REPLICATE.STRESS_TENSORS_Z 	= SYSTEM.STRESS_TENSORS_Z;
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
		REPLICATE.BOXDIM.SCALE_BY_FACTOR(lscale, true, REPLICATE.COORDS[j]);

	for (int j=0; j < REPLICATE.ALL_ATOMS; j++ ) 
		REPLICATE.BOXDIM.SCALE_BY_FACTOR(lscale, true, REPLICATE.ALL_COORDS[j]);
	
	REPLICATE.BOXDIM.SCALE_BY_FACTOR(lscale);
	REPLICATE.BOXDIM.UPDATE_CELL();
	

	// Compute/store new total potential energy and volume
	
	Vol1 = REPLICATE.BOXDIM.VOL;
	
	REPLICATE.TOT_POT_ENER = 0;
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);
	
	Vtot1 = REPLICATE.TOT_POT_ENER;
	
	// Contract coords by a bit

	lscale = 1.0  - eps;
	
	REPLICATE_SYSTEM(SYSTEM, REPLICATE);

	for (int j=0; j < REPLICATE.ATOMS; j++ ) 
		REPLICATE.BOXDIM.SCALE_BY_FACTOR(lscale, true, REPLICATE.COORDS[j]);
	
	for (int j=0; j < REPLICATE.ALL_ATOMS; j++ ) 
		REPLICATE.BOXDIM.SCALE_BY_FACTOR(lscale, true, REPLICATE.ALL_COORDS[j]);

	REPLICATE.BOXDIM.SCALE_BY_FACTOR(lscale);
	REPLICATE.BOXDIM.UPDATE_CELL();
	
	
	// Compute/store new total potential energy and volume
	
	Vol2 = REPLICATE.BOXDIM.VOL;
	
	REPLICATE.TOT_POT_ENER = 0;
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);
	
	Vtot2 = REPLICATE.TOT_POT_ENER;

	// compute (return) pressure 
	
	//return -(Vtot2 - Vtot1)/(Vol2 - Vol1);
	
	PE_1 = Vtot1;
	PE_2 = Vtot2;
	dV   = Vol2 - Vol1;
}



void numerical_stress(const FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY,
											CLUSTER_LIST & TRIPS,  CLUSTER_LIST &QUADS, map<string,int> & PAIR_MAP,
											vector<int> &INT_PAIR_MAP,
											NEIGHBORS & NEIGHBOR_LIST,double & PE_1, double &PE_2, 
											int it1, int it2)
// Evaluates the configurational part of the stress numerically by -1/V dU/d epsilon.
// Essentially, we are taking the system and distorting it a bit (lscale) to get the change in potential energy
// Two partial potential energies, PE_1 and PE_2 are returned.  Values need to be summed across MPI processes
// to obtain the final stress tensor component.
// it1 and it2 are the requested tensor indices of the stress.
{
	const double eps = 1.0e-04;
	double Vtot1, Vtot2;
	double Vol1, Vol2;

	// Change this if a non-Ewald coulomb option is added.
	if ( CONTROLS.USE_COULOMB )
		EXIT_MSG("Error: Numerical stress with Ewald is not supported") ;
	
	// Make a copy of the system

	vector<vector<double>> lscale(3, vector<double>(3)) ;


	double V0 = SYSTEM.BOXDIM.VOL ;
	
	static FRAME REPLICATE;
		
	REPLICATE_SYSTEM(SYSTEM, REPLICATE);		// Can we make this one of those "if ! called_before sorts of variables?"

	// Expand coords by a bit (lscale)

	for ( int i = 0 ; i < 3 ; i++ ) {
		for ( int j = 0 ; j < 3 ; j++ ) {
			if ( i != j ) 
				lscale[i][j] = 0.0 ;
		}
		lscale[i][i] = 1.0 ;
	}

	lscale[it1][it2] += eps ;
	//if ( it1 != it2 ) 
	//lscale[it2][it1] += eps ;

	for (int j=0; j < REPLICATE.ATOMS; j++ )		
		REPLICATE.BOXDIM.SCALE_BY_MATRIX(lscale, true, REPLICATE.COORDS[j]);

	for (int j=0; j < REPLICATE.ALL_ATOMS; j++ ) 
		REPLICATE.BOXDIM.SCALE_BY_MATRIX(lscale, true, REPLICATE.ALL_COORDS[j]);
	
	REPLICATE.BOXDIM.SCALE_BY_MATRIX(lscale, false, REPLICATE.COORDS[0]) ;
	REPLICATE.BOXDIM.UPDATE_CELL();
	

	// Compute/store new total potential energy and volume
	
	REPLICATE.TOT_POT_ENER = 0;
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);
	
	Vtot1 = REPLICATE.TOT_POT_ENER;
	
	// Contract coords by a bit

	lscale[it1][it2] -= 2.0 * eps ;

	//if ( it1 != it2 ) 
	//lscale[it2][it1] -= 2.0 * eps ;
	
	REPLICATE_SYSTEM(SYSTEM, REPLICATE);

	for (int j=0; j < REPLICATE.ATOMS; j++ ) 
		REPLICATE.BOXDIM.SCALE_BY_MATRIX(lscale, true, REPLICATE.COORDS[j]);
	
	for (int j=0; j < REPLICATE.ALL_ATOMS; j++ ) 
		REPLICATE.BOXDIM.SCALE_BY_MATRIX(lscale, true, REPLICATE.ALL_COORDS[j]);

	REPLICATE.BOXDIM.SCALE_BY_MATRIX(lscale, false, REPLICATE.COORDS[0]);
	REPLICATE.BOXDIM.UPDATE_CELL();
	
	
	// Compute/store new total potential energy and volume
	
	REPLICATE.TOT_POT_ENER = 0;
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);
	
	Vtot2 = REPLICATE.TOT_POT_ENER;

	// compute contribution to stress.
	
	PE_1 = Vtot1 / (2.0 * eps * V0) ;
	PE_2 = Vtot2 / (2.0 * eps * V0) ;
}

void numerical_stress_all(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY,
													CLUSTER_LIST & TRIPS,  CLUSTER_LIST &QUADS, map<string,int> & PAIR_MAP,
													vector<int> &INT_PAIR_MAP,
													NEIGHBORS & NEIGHBOR_LIST)
// Evaluates the configurational part of the stress numerically by -1/V dU/d epsilon.
{
	double PE_1, PE_2, dV; 

	// Loop over cartesion components.
	for ( int tidx0 = 0 ; tidx0 < 3 ; tidx0++ )
		{
			for ( int tidx1 = 0 ; tidx1 <= tidx0 ; tidx1++ )
				{
		 
					numerical_stress(SYSTEM, CONTROLS, FF_2BODY, TRIPS, QUADS, PAIR_MAP, INT_PAIR_MAP,
													 NEIGHBOR_LIST, PE_1, PE_2, tidx0, tidx1);
			
#ifdef USE_MPI
					double pe1_sum = 0.0, pe2_sum = 0.0 ;
					MPI_Allreduce(&PE_1, &pe1_sum,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
					PE_1 = pe1_sum ;
					MPI_Allreduce(&PE_2, &pe2_sum,1,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
					PE_2 = pe2_sum ;
#endif

					switch ( tidx1 ) {
					case 0:
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[tidx0].X = (PE_2 - PE_1);
						break ;
					case 1:
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[tidx0].Y = (PE_2 - PE_1);
						break ;
					case 2:
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[tidx0].Z = (PE_2 - PE_1);
						break ;
					default:
						cout << "Error: bad tensor index\n" ;
						exit_run(1) ;
					}
					switch ( tidx0 ) {
					case 0:
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[tidx1].X = (PE_2 - PE_1);
						break ;
					case 1:
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[tidx1].Y = (PE_2 - PE_1);
						break ;
					case 2:
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[tidx1].Z = (PE_2 - PE_1);
						break ;
					default:
						cout << "Error: bad tensor index\n" ;
						exit_run(1) ;
					}
				}
		}
}

void check_forces(FRAME& SYSTEM, JOB_CONTROL &CONTROLS, vector<PAIR_FF> &FF_2BODY, map<string,int>& PAIR_MAP, vector<int> &INT_PAIR_MAP, CLUSTER_LIST &TRIPS, CLUSTER_LIST &QUADS, NEIGHBORS &NEIGHBOR_LIST)
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

		SYSTEM.update_ghost(CONTROLS.N_LAYERS, false) ;

		// Recalculate the forces
		NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);

		ZCalc(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);

		// FOR MPI:		Synchronize forces, energy, and pressure.
#ifdef USE_MPI
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

		double energy1 = SYSTEM.TOT_POT_ENER ;

		// Move the atom again.
		if ( j == 0 ) 
		  SYSTEM.COORDS[a1].X = coords[a1].X - eps ;
		else if ( j == 1 ) 
		  SYSTEM.COORDS[a1].Y = coords[a1].Y - eps ;
		else if ( j == 2 )
		  SYSTEM.COORDS[a1].Z = coords[a1].Z - eps ;

		SYSTEM.update_ghost(CONTROLS.N_LAYERS, false) ;

		NEIGHBOR_LIST.UPDATE_LIST(SYSTEM, CONTROLS);

		ZCalc(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, TRIPS, QUADS, NEIGHBOR_LIST);

		// FOR MPI:		Synchronize forces, energy, and pressure.
#ifdef USE_MPI
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
  SYSTEM.update_ghost(CONTROLS.N_LAYERS, false) ;

}

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
	
	if(CONTROLS.FIT_STRESS || CONTROLS.FIT_STRESS_ALL || CONTROLS.FIT_ENER)
	{
		if(CONTROLS.IF_SUBTRACT_COORD || CONTROLS.IF_SUBTRACT_COUL)
		{
			cout << "ERROR: The following options must be false when including energies or stresses in fit: " << endl;
			cout << "CONTROLS.FIT_STRESS     " << endl;
			cout << "CONTROLS.FIT_STRESS_ALL" << endl;
			cout << "CONTROLS.FIT_ENER      " << endl;
		}
	}
	
	if ( CONTROLS.FIT_COUL ) 
		ZCalc_Ewald_Deriv(FRAME_SYSTEM, FF_2BODY, A_MATRIX, PAIR_MAP, NEIGHBOR_LIST, CONTROLS);	
	
	if ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" )
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
    else 
    {
		cout << "Error: bad pairtype in ZCalc_Deriv\n";
		exit(1);
    }		
}	


////////////////////////////////////////////////////////////
//
// FUNCTIONS -- FORCE CALCULATION
//
////////////////////////////////////////////////////////////
 
static void ZCalc_Lj(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, vector<int> &INT_PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST); 

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
	
	SYSTEM.PRESSURE_TENSORS_XYZ_ALL.resize(3);
	SYSTEM.PRESSURE_TENSORS_ALL    .resize(3);
  
	for(int i=0;i<3; i++)
	{
		SYSTEM.PRESSURE_TENSORS_XYZ_ALL[i].X = 0.0;
		SYSTEM.PRESSURE_TENSORS_XYZ_ALL[i].Y = 0.0;
		SYSTEM.PRESSURE_TENSORS_XYZ_ALL[i].Z = 0.0;
		SYSTEM.PRESSURE_TENSORS_ALL    [i].X = 0.0;
		SYSTEM.PRESSURE_TENSORS_ALL    [i].Y = 0.0;
		SYSTEM.PRESSURE_TENSORS_ALL    [i].Z = 0.0;
		
	}
	if ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" ) 
	{

	  if(CONTROLS.INCLUDE_ATOM_OFFSETS)
	  {
		  // Add the per-atom contributions to energy, if requested
	  
	  	int a1start, a1end;
	  	divide_atoms(a1start, a1end, SYSTEM.ATOMS);
	  
	  	for(int a=a1start;a<=a1end;a++)
	  		SYSTEM.TOT_POT_ENER += SYSTEM.QM_ENERGY_OFFSET[ SYSTEM.ATOMTYPE_IDX[a] ];
	  }

	  Cheby cheby{CONTROLS, SYSTEM, NEIGHBOR_LIST, FF_2BODY, INT_PAIR_MAP};
	  cheby.Force_all(TRIPS, QUADS);
	  
	}
	else if ( FF_2BODY[0].PAIRTYP == "LJ" )
	{
	  ZCalc_Lj(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, INT_PAIR_MAP, NEIGHBOR_LIST);
	}
	else if ( FF_2BODY[0].PAIRTYP == "SERIAL_CHIMES" )
	{
		 ZCalc_Serial_Chimes(SYSTEM, FF_2BODY[0]) ;
	}
	else 
	{
		 cout << "Error: bad pairtype in ZCalc: " << FF_2BODY[0].PAIRTYP << endl;
		 exit_run(1);
	}	
	
	if ( CONTROLS.USE_COULOMB ) 
		ZCalc_Ewald(SYSTEM, CONTROLS, NEIGHBOR_LIST);
		

	SYSTEM.PRESSURE_XYZ           /= 3.0 *  SYSTEM.BOXDIM.VOL;
	
	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X /= SYSTEM.BOXDIM.VOL;
	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Y /= SYSTEM.BOXDIM.VOL;
	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Z /= SYSTEM.BOXDIM.VOL;

	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X /= SYSTEM.BOXDIM.VOL;
	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y /= SYSTEM.BOXDIM.VOL;
	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Z /= SYSTEM.BOXDIM.VOL;

	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].X /= SYSTEM.BOXDIM.VOL;
	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y /= SYSTEM.BOXDIM.VOL;
	SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z /= SYSTEM.BOXDIM.VOL;

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

	 double perm_scale = NEIGHBOR_LIST.PERM_SCALE[2] ;
	
	 divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.

	 for(int a1=a1start;a1 <= a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS 
	 {
			a2start = 0;
			a2end   = NEIGHBOR_LIST.LIST[a1].size();
		
			for(int a2idx =a2start; a2idx < a2end;a2idx++)
			{			
				 int a2 = NEIGHBOR_LIST.LIST[a1][a2idx];
					
				 curr_pair_type_idx =  INT_PAIR_MAP[SYSTEM.ATOMTYPE_IDX[a1]*CONTROLS.NATMTYP + SYSTEM.ATOMTYPE_IDX[SYSTEM.PARENT[a2]]];

				 // pair interaction cutoff distance.
				 double rcutoff = FF_2BODY[curr_pair_type_idx].S_MAXIM;

				 rlen_mi = get_dist(SYSTEM, RVEC, a1, a2);


				 if ( rlen_mi < FF_2BODY[curr_pair_type_idx].PARAMS[1]/2.2) 
						EXIT_MSG("Error: close approach", rlen_mi);

				 else if ( rlen_mi < rcutoff ) 
				 {
						SYSTEM.TOT_POT_ENER += perm_scale * 4.0 * FF_2BODY[curr_pair_type_idx].PARAMS[0] * ( pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,12.0) - pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,6.0) );
						fac = perm_scale * 4.0 * FF_2BODY[curr_pair_type_idx].PARAMS[0] * ( -12.0 * pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,14.0) + 6.0 *    pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,8.0) );
						fac *= 1.0 / ( FF_2BODY[curr_pair_type_idx].PARAMS[1] * FF_2BODY[curr_pair_type_idx].PARAMS[1] );		

						SYSTEM.PRESSURE_XYZ -= fac * (rlen_mi*rlen_mi);
				
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X -= fac * rlen_mi * RVEC.X * RVEC.X / rlen_mi;
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y -= fac * rlen_mi * RVEC.Y * RVEC.Y / rlen_mi;
						SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z -= fac * rlen_mi * RVEC.Z * RVEC.Z / rlen_mi;
	
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

void sum_forces(vector<XYZ>& accel_vec, int atoms, double &pot_energy, double &pressure,
								double &tens_xx, double &tens_xy,	double &tens_xz,
								double &tens_yx, double &tens_yy,	double &tens_yz,
								double &tens_zx, double &tens_zy,	double &tens_zz)
	
	// Add up forces, potential energy, and pressure from all processes.  
	{
		// Sum up the potential energy, pressure, tensors , and forces from all processors
#ifdef USE_MPI	

		int ndata = 11 + 3 * atoms ;
		vector<double> buf(ndata), sum(ndata,0.0) ;
		buf[0] = pot_energy ;
		buf[1] = pressure ;
		buf[2] = tens_xx ;
		buf[3] = tens_xy ;
		buf[4] = tens_xz ;
		buf[5] = tens_yx ;
		buf[6] = tens_yy ;
		buf[7] = tens_yz ;
		buf[8] = tens_zx ;
		buf[9] = tens_zy ;
		buf[10] = tens_zz ;

		for ( int j = 0 ; j < atoms ; j++ )
		{
			buf[11+3*j] = accel_vec[j].X ;
			buf[12+3*j] = accel_vec[j].Y ;
			buf[13+3*j] = accel_vec[j].Z ;						
		}

		// memcheck_all complained of memory errors when MPI_IN_PLACE was used,
		// so I switched to explicit buffer allocation (LEF 3/3/22)
		MPI_Allreduce(buf.data(), sum.data(), ndata, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

		pot_energy = sum[0] ;
		pressure = sum[1] ;
		tens_xx = sum[2] ;
		tens_xy = sum[3] ;
		tens_xz = sum[4] ;
		tens_yx = sum[5] ;
		tens_yy = sum[6] ;
		tens_yz = sum[7] ;
		tens_zx = sum[8] ;
		tens_zy = sum[9] ;
		tens_zz = sum[10] ;

		for ( int j = 0 ; j < atoms ; j++ )
		{
			accel_vec[j].X = sum[11+3*j] ;
			accel_vec[j].Y = sum[12+3*j] ;
			accel_vec[j].Z = sum[13+3*j] ;						
		}		

#endif	

	}

#ifdef USE_MPI
void sync_position(vector<XYZ>& coord_vec, NEIGHBORS & neigh_list, vector<XYZ>& velocity_vec,
									 int atoms, bool sync_vel, BOX &BOXDIM) 
	// Broadcast the position, neighborlists,  and optionally the velocity to all nodes.
	// Velocity-broadcast is only necessary for velocity-dependent forces (not currently implemented).
{
	 // Convert out vectors to arrays so they are MPI friendly

	 int buf_sz = 3 * atoms ;

	 if ( sync_vel ) buf_sz *= 2 ;

	 if ( BOXDIM.IS_VARIABLE ) buf_sz += 3 ;

	 // One extra for MAX_COORD_STEP.
	 buf_sz++ ;
	 
	 vector<double> buffer(buf_sz) ;

	 for ( int i = 0 ; i < atoms ; i++ )
	 {
			buffer[3*i]   = coord_vec[i].X ;
			buffer[3*i+1] = coord_vec[i].Y ;
			buffer[3*i+2] = coord_vec[i].Z ;
	 }
	 int count = 3 * atoms ;
	 
	 if ( sync_vel )
	 {
			count += 3 * atoms ;
			for ( int i = 0 ; i < atoms ; i++ )
			{
				 buffer[3*atoms+3*i]   = velocity_vec[i].X ;
				 buffer[3*atoms+3*i+1] = velocity_vec[i].Y ;
				 buffer[3*atoms+3*i+2] = velocity_vec[i].Z ;				 				 
			}
	 }

	 if ( BOXDIM.IS_VARIABLE )
	 {
			if ( BOXDIM.IS_ORTHO )
			{
				 buffer[count]   = BOXDIM.CELL_AX ;
				 buffer[count+1] = BOXDIM.CELL_BY ;
				 buffer[count+2] = BOXDIM.CELL_CZ ;
			}
			else
			{
				 EXIT_MSG("Variable non-orthorhombic boxes are not supported") ;
			}
	 }

	 buffer[buf_sz-1] = neigh_list.MAX_COORD_STEP ;
	 
	 MPI_Bcast(buffer.data(), buf_sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	 // Unpack the buffer
	 if ( RANK != 0 )
	 {
			for ( int i = 0 ; i < atoms ; i++ )
			{
				 coord_vec[i].X = buffer[3*i] ;
				 coord_vec[i].Y = buffer[3*i+1] ;
				 coord_vec[i].Z = buffer[3*i+2] ;
			}
	 
			if ( sync_vel )
			{
				 for ( int i = 0 ; i < atoms ; i++ )
				 {
						velocity_vec[i].X = buffer[3*atoms+3*i] ;
						velocity_vec[i].Y = buffer[3*atoms+3*i+1] ;
						velocity_vec[i].Z = buffer[3*atoms+3*i+2] ;				 				 
				 }
			}

			if ( BOXDIM.IS_VARIABLE )
			{
				 BOXDIM.CELL_AX = buffer[count] ;
				 BOXDIM.CELL_BY = buffer[count+1] ;
				 BOXDIM.CELL_CZ = buffer[count+2] ;
				 BOXDIM.UPDATE_CELL();			
			}

			neigh_list.MAX_COORD_STEP = buffer[buf_sz-1] ;
	 }
	 
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
			
#endif // LOG_POS
}

#endif // USE_MPI


static void ZCalc_Serial_Chimes(FRAME &SYSTEM, PAIR_FF &FF_2BODY)
// Force evaluation using serial Chimes calculator interface.
{
    // Only rank 0 does calculation for serial chimes.
    if ( RANK == 0 ) {
        int natoms = SYSTEM.ATOMS ;
        vector<double> xcrds(natoms) ;
        vector<double> ycrds(natoms) ;
        vector<double> zcrds(natoms) ;
        vector<vector<double>> force(natoms) ;
        vector<double> cell_a(3);
        vector<double> cell_b(3);
        vector<double> cell_c(3);
        vector<double> stress(9,0.0); 
        double energy = 0.0 ;
		 
        vector<string> atom_types(SYSTEM.ATOMTYPE.size()) ;
		 
        for ( int j = 0 ; j < natoms ; j++ )
        {
            xcrds[j] = SYSTEM.COORDS[j].X ;
            ycrds[j] = SYSTEM.COORDS[j].Y ;
            zcrds[j] = SYSTEM.COORDS[j].Z ;								
            force[j].resize(3,0.0);
        }

        for ( int j = 0 ; j < atom_types.size() ; j++ )
        {
            atom_types[j] = SYSTEM.ATOMTYPE[j] ;
        }
		 
        cell_a[0] = SYSTEM.BOXDIM.CELL_AX ;
        cell_a[1] = SYSTEM.BOXDIM.CELL_AY ;
        cell_a[2] = SYSTEM.BOXDIM.CELL_AZ ;		 		 

        cell_b[0] = SYSTEM.BOXDIM.CELL_BX ;
        cell_b[1] = SYSTEM.BOXDIM.CELL_BY ;
        cell_b[2] = SYSTEM.BOXDIM.CELL_BZ ;		 		 

        cell_c[0] = SYSTEM.BOXDIM.CELL_CX ;
        cell_c[1] = SYSTEM.BOXDIM.CELL_CY ;
        cell_c[2] = SYSTEM.BOXDIM.CELL_CZ ;		 		 

        FF_2BODY.chimes.calculate(xcrds, ycrds, zcrds, cell_a, cell_b, cell_c, atom_types, energy, force, stress);

        for ( int j = 0 ; j < natoms ; j++ )
        {
            SYSTEM.ACCEL[j].X = force[j][0] ;
            SYSTEM.ACCEL[j].Y = force[j][1] ;
            SYSTEM.ACCEL[j].Z = force[j][2] ;
        }
        SYSTEM.TOT_POT_ENER = energy ;

        // Unit conversions to CHIMES_MD internals.
        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X = stress[0] * SYSTEM.BOXDIM.VOL ;
        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Y = stress[1] * SYSTEM.BOXDIM.VOL ;
        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].Z = stress[2] * SYSTEM.BOXDIM.VOL ;

        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].X = stress[3] * SYSTEM.BOXDIM.VOL ;
        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y = stress[4] * SYSTEM.BOXDIM.VOL ;
        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Z = stress[5] * SYSTEM.BOXDIM.VOL ;

        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].X = stress[6] * SYSTEM.BOXDIM.VOL ;
        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Y = stress[7] * SYSTEM.BOXDIM.VOL ;
        SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z = stress[8] * SYSTEM.BOXDIM.VOL ;		 		 		 		 
        SYSTEM.PRESSURE_XYZ = SYSTEM.PRESSURE_TENSORS_XYZ_ALL[0].X
            + SYSTEM.PRESSURE_TENSORS_XYZ_ALL[1].Y
            + SYSTEM.PRESSURE_TENSORS_XYZ_ALL[2].Z ;
    }
}


