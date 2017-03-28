#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes

#include "functions.h"

using namespace std;




NEIGHBORS::NEIGHBORS()		// Constructor 
{
	RCUT_PADDING  =  0.3;
	DISPLACEMENT  =  0.0;
	MAX_CUTOFF    =  0.0;
	MAX_CUTOFF_3B =  0.0;
	FIRST_CALL    = true;
	SECOND_CALL   = true;
	USE           = false;
	MAX_VEL       =  0.0;
	CURR_VEL      =  0.0;
	
	// New from Larry
	EWALD_CUTOFF  =  0.0;
	UPDATE_FREQ   =  30.0;
	SAFETY        =  1.01;
	
}
NEIGHBORS::~NEIGHBORS(){}	// Deconstructor

void NEIGHBORS::INITIALIZE(FRAME & SYSTEM)		// (overloaded) class constructor -- if no padding specified, default to 0.3
{
	if(USE)
	{
		LIST          .resize(SYSTEM.ATOMS);
		LIST_EWALD    .resize(SYSTEM.ATOMS);
		LIST_UNORDERED.resize(SYSTEM.ATOMS);
		LIST_3B       .resize(SYSTEM.ALL_ATOMS);
	}
}

void NEIGHBORS::INITIALIZE(FRAME & SYSTEM, double & PAD)	// (overloaded) class constructor -- if padding specified, set to value
{
	INITIALIZE(SYSTEM);
	
	if (USE)
		RCUT_PADDING = PAD;
	else 
		RCUT_PADDING = 1.0e+10;
}

void NEIGHBORS::DO_UPDATE(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{
	XYZ RAB;
	XYZ TMP_BOX;
	
	double MAX  = 0;
	double rlen = 0;
	
	// ALL_COORDS is used for force evaluation. Wrap atoms in ALL_COORDS into primitive cell prior to constructing replicas.
	for ( int a = 0; a < SYSTEM.ATOMS; a++ ) 
	{
		SYSTEM.ALL_COORDS[a].X = SYSTEM.COORDS[a].X - floor(SYSTEM.COORDS[a].X / SYSTEM.BOXDIM.X) * SYSTEM.BOXDIM.X;
		SYSTEM.ALL_COORDS[a].Y = SYSTEM.COORDS[a].Y - floor(SYSTEM.COORDS[a].Y / SYSTEM.BOXDIM.Y) * SYSTEM.BOXDIM.Y;
		SYSTEM.ALL_COORDS[a].Z = SYSTEM.COORDS[a].Z - floor(SYSTEM.COORDS[a].Z / SYSTEM.BOXDIM.Z) * SYSTEM.BOXDIM.Z; 
	}
	
	sync_layers(SYSTEM,CONTROLS);
	
	if(FIRST_CALL)									// Set up the first dimension of the list 
	{
		LIST          .resize(SYSTEM.ATOMS);	
		LIST_EWALD    .resize(SYSTEM.ATOMS);	
		LIST_UNORDERED.resize(SYSTEM.ATOMS);	
		LIST_3B       .resize(SYSTEM.ALL_ATOMS);	
	}

	for(int a1=0; a1<SYSTEM.ATOMS; a1++)
	{
		if(!FIRST_CALL)								// Clear out the second dimension so we can start over again
		{										
			vector<int>().swap(LIST   [a1]);	
			vector<int>().swap(LIST_UNORDERED[a1]);
			vector<int>().swap(LIST_EWALD[a1]);
		}
		
		// Search across a2 < a1 for unordered list
		for (int a2 = 0 ; a2< a1 ; a2++)
		{			
			// Get pair distance

			rlen = get_dist(SYSTEM, RAB, a1, a2);		// Updates RAB!

			// I get the correct result when the list contains all a2 achievable
			// at this point in the function, regardless of update frequency
			// i.e. "if(true)"

			if (rlen < MAX_CUTOFF + RCUT_PADDING)		// Then add it to the atoms neighbor list (2B)
				LIST_UNORDERED[a1].push_back(a2);

		}
		
		// Search across a2 > a1 for both regular and unordered lists.
		for (int a2=a1+1 ; a2<SYSTEM.ALL_ATOMS; a2++)
		{			
			// Get pair distance

			rlen = get_dist(SYSTEM, RAB, a1, a2);		// Updates RAB!

			// I get the correct result when the list contains all a2 achievable
			// at this point in the function, regardless of update frequency
			// i.e. "if(true)"

			if (rlen < (MAX_CUTOFF + RCUT_PADDING) )	// Then add it to the atoms neighbor list (2B)
			{
				if (a1 <= SYSTEM.PARENT[a2])			// Select atoms in neighbor list according to parents.
					LIST[a1].push_back(a2);		

				LIST_UNORDERED[a1].push_back(a2) ;
			}

			if (rlen < (EWALD_CUTOFF + RCUT_PADDING) )	// Then add it to the atoms neighbor list (2B)
				if ( a1 <= SYSTEM.PARENT[a2] ) 			// Select atoms in neighbor list according to parents.
					LIST_EWALD[a1].push_back(a2);		
		}
	}
	
	// Calculate neighbors of images for the 3 body list only.
	for(int a1= SYSTEM.ATOMS ; a1<SYSTEM.ALL_ATOMS; a1++)
	{
		if(!FIRST_CALL)	
			vector<int>().swap(LIST_3B[a1]); // Clear out the second dimension so we can start over again

		// Search across a2 > a1 for both regular and unordered lists.
		for (int a2 = 0 ; a2<SYSTEM.ALL_ATOMS; a2++)
		{			
			if ( a1 == a2 ) 
				continue ;
			
			rlen = get_dist(SYSTEM, RAB, a1, a2);	// Get pair distance

			if(rlen < MAX_CUTOFF_3B + RCUT_PADDING)	
				if ( SYSTEM.PARENT[a1] <= SYSTEM.PARENT[a2] ) 
					LIST_3B[a1].push_back(a2);	// Then add it to the atoms neighbor list (3B)
		}
	}
	
	FIRST_CALL = false;	
}

void NEIGHBORS::UPDATE_LIST(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{	
	if(FIRST_CALL)	// Then start from scratch by cycling through all atom (and layer atom) pairs
	{
		DO_UPDATE(SYSTEM, CONTROLS);
		FIRST_CALL = false;
		
		if (!USE) 
			RCUT_PADDING = 1.0e10 ;
		
		if(RANK==0)
		{
			cout << "NEIGHBOR LIST MAX CUTOFF IS: " << fixed << setprecision(5) << MAX_CUTOFF << endl;
			cout << "USING PADDING: " << fixed << setprecision(3) << RCUT_PADDING << endl;
		}
	}
	/* For debugging
	else if(!USE)	// If we don't want to use neighbor lists, this is the same as updating it each time.
	{
		DO_UPDATE(SYSTEM, CONTROLS);
	}
	*/
	else
	{
		// Note: We keep track of the max velocity in the part of splines_md.C 
		// where the first half of thermostatting is applied

		if( (SECOND_CALL))	// Update the cutoff padding
		{
			if(USE)
				RCUT_PADDING = MAX_VEL * UPDATE_FREQ * CONTROLS.DELTA_T;	// should give a distance in AA

			SECOND_CALL = false;
			
			DO_UPDATE(SYSTEM, CONTROLS);
			
			if(RANK == 0)
				cout << "RANK: " << RANK << " RESET RCUT_PADDING TO: " << fixed << RCUT_PADDING << endl;
		}
		else
		{
			DISPLACEMENT += MAX_VEL * CONTROLS.DELTA_T;
			
			if(DISPLACEMENT>0.5*RCUT_PADDING)
			{
				if (USE) 
					RCUT_PADDING = MAX_VEL * UPDATE_FREQ * CONTROLS.DELTA_T;	// Update padding in case max_vel changed. (LEF).
				
				DO_UPDATE(SYSTEM, CONTROLS);
				
				DISPLACEMENT = 0;
				MAX_VEL      = 0;

				if(RANK == 0)
					cout << "RANK: " << RANK << " UPDATING ON STEP: " << CONTROLS.STEP << ", WITH PADDING: " << fixed << setprecision(3) << RCUT_PADDING <<  endl;
				
			}
		}
	}
}
void NEIGHBORS::UPDATE_LIST(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, bool FORCE)
{
	if(FIRST_CALL || SECOND_CALL)
		UPDATE_LIST(SYSTEM, CONTROLS);
	else
	{
		DISPLACEMENT += MAX_VEL * CONTROLS.DELTA_T;
		
		DO_UPDATE(SYSTEM, CONTROLS);
	
		if(RANK == 0)
			cout << "RANK: " << RANK << " FORCING UPDATING ON STEP: " << CONTROLS.STEP << ", WITH PADDING: " << fixed << setprecision(3) << RCUT_PADDING <<  endl;
	}
}



CONSTRAINT::CONSTRAINT(){}	// Constructor
CONSTRAINT::~CONSTRAINT(){}	// Deconstructor

void CONSTRAINT::INITIALIZE(string IN_STYLE, JOB_CONTROL & CONTROLS, int ATOMS)
{
	STYLE = IN_STYLE;
	
	if(IN_STYLE=="NPT-MTK")							// Attempts to use MTK Thermostat and barostat
		STYLE = IN_STYLE;
	else if(IN_STYLE=="NPT-BEREND")					// Uses position scaling to barostat
		STYLE = IN_STYLE;
	else if(IN_STYLE=="NPT-BEREND-ANISO")			// Uses position scaling to barostat anisotropically, but expects alpha = beta = gamma = 90 deg
		STYLE = IN_STYLE;	
	else if(IN_STYLE=="NVT-BEREND")					// Uses velocity scaling to thermostat
		STYLE = IN_STYLE;
	else if(CONTROLS.USE_HOOVER_THRMOSTAT)			// Uses MTK Thermostat
		STYLE = "NVT-MTK";
	else if(CONTROLS.FREQ_UPDATE_THERMOSTAT > -1.0)	// Trivial velocity scaling
		STYLE = "NVT-SCALE";
	else
		STYLE="NVE";
	
	TIME  = CONTROLS.FREQ_UPDATE_THERMOSTAT;
	TIME_BARO = CONTROLS.FREQ_UPDATE_BAROSTAT;

	if(RANK == 0)
	{
		cout << "	...Configuring constraints for a " << STYLE << " simulation." << endl;
		
		if(STYLE != "NVE")
			cout << "	...Setting thermostat time constant to: " << TIME << endl;
		if((STYLE== "NPT-BEREND")||(STYLE== "NPT-BEREND-ANISO")||(STYLE== "NPT-MTK"))
			cout << "	...Setting barostat time constant to: " << TIME_BARO << endl;
	}

	N_DOF = 3*ATOMS - 3;
	THERM_INERT_Q  = N_DOF * Kb * CONTROLS.TEMPERATURE * TIME * TIME / Tfs / Tfs; // Need to convert from fs to sim units ... omega (Frequency) is inverse time
	BAROS_INERT_W = (N_DOF + 3.0) * Kb * CONTROLS.TEMPERATURE * (TIME_BARO) * (TIME_BARO)/ Tfs / Tfs;
	
	BEREND_KP = TIME_BARO/ Tfs;	// Barostat damping parameter... 
	BEREND_TAU = TIME / Tfs;	// Thermostat damping parameter
	
	THERM_POSIT_T = 0;
	THERM_VELOC_T = 0;
	THERM_INERT_T = 0;
	
	BAROS_POSIT_T = 0;
	BAROS_VELOC_T = 0;
	BAROS_FORCE_T = 0;
	
	THERM_POSIT_0 = 0;
	THERM_VELOC_0 = 0;
	THERM_INERT_0 = 0;
	
	BAROS_POSIT_0 = 0;
	BAROS_VELOC_0 = 0;
	BAROS_FORCE_0 = 0;
	
	VOLUME_0 = 0;
	VOLUME_T = 0;
	
	BAROS_SCALE = 1;

	
}

void CONSTRAINT::UPDATE_COORDS(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{
	THERM_POSIT_0 = THERM_POSIT_T;
	THERM_VELOC_0 = THERM_VELOC_T;
	THERM_INERT_0 = THERM_INERT_T;
	
	BAROS_POSIT_0 = BAROS_POSIT_T;
	BAROS_VELOC_0 = BAROS_VELOC_T;
	BAROS_FORCE_0 = BAROS_FORCE_T;
	
	VOLUME_0 = VOLUME_T;

	///////////////////////////////////////////////
	// Do the un-constrained update of coords... this applies to all styles.
	///////////////////////////////////////////////

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
		if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
			continue;
		
		SYSTEM.COORDS[a1].X += SYSTEM.VELOCITY[a1].X * CONTROLS.DELTA_T + 0.5*SYSTEM.ACCEL[a1].X * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
		SYSTEM.COORDS[a1].Y += SYSTEM.VELOCITY[a1].Y * CONTROLS.DELTA_T + 0.5*SYSTEM.ACCEL[a1].Y * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
		SYSTEM.COORDS[a1].Z += SYSTEM.VELOCITY[a1].Z * CONTROLS.DELTA_T + 0.5*SYSTEM.ACCEL[a1].Z * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
	}
	
	if(STYLE=="NVE" || STYLE=="NVT-SCALE" || STYLE=="NVT-BEREND")
		return;

	///////////////////////////////////////////////
	// If requested, do Berendsen barostatting, isotropically
	///////////////////////////////////////////////
	
	if(STYLE=="NPT-BEREND")
	{
		VOLUME_0 = SYSTEM.BOXDIM.X*SYSTEM.BOXDIM.Y*SYSTEM.BOXDIM.Z;
		BEREND_MU = pow(1.0 - CONTROLS.DELTA_T/BEREND_KP*(CONTROLS.PRESSURE-SYSTEM.PRESSURE)/GPa,1.0/3.0);

		if(BEREND_MU < 0)
		{
			cout << "ERROR: Negative Berendsen scaling factor computed. Increase BAROSCALE." << endl;
			cout << SYSTEM.PRESSURE << " " << CONTROLS.PRESSURE << " " << BEREND_MU << endl;
			exit_run(0);
		}
	
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)	
		{	
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
			{
				// Atom freezing doesn't work for a barostat!
				cout << "ERROR: Atoms cannot be frozen in an NPT ensemble." << endl;
				exit_run(0);
			}
			
			SYSTEM.COORDS[a1].X *= BEREND_MU;	
			SYSTEM.COORDS[a1].Y *= BEREND_MU;
			SYSTEM.COORDS[a1].Z *= BEREND_MU;	
		}
	
		SYSTEM.BOXDIM.X *= BEREND_MU;
		SYSTEM.BOXDIM.Y *= BEREND_MU;
		SYSTEM.BOXDIM.Z *= BEREND_MU;

		VOLUME_T = SYSTEM.BOXDIM.X*SYSTEM.BOXDIM.Y*SYSTEM.BOXDIM.Z;
		
		return;
	}
	
	///////////////////////////////////////////////
	// If requested, do Berendsen barostatting, anisotropically
	///////////////////////////////////////////////
	
	if(STYLE=="NPT-BEREND-ANISO")
	{
		VOLUME_0 = SYSTEM.BOXDIM.X*SYSTEM.BOXDIM.Y*SYSTEM.BOXDIM.Z;
		
		BEREND_ANI_MU.X = pow(1.0 - CONTROLS.DELTA_T/BEREND_KP*(CONTROLS.PRESSURE-SYSTEM.PRESSURE_TENSORS.X)/GPa,1.0/3.0);
		BEREND_ANI_MU.Y = pow(1.0 - CONTROLS.DELTA_T/BEREND_KP*(CONTROLS.PRESSURE-SYSTEM.PRESSURE_TENSORS.Y)/GPa,1.0/3.0);
		BEREND_ANI_MU.Z = pow(1.0 - CONTROLS.DELTA_T/BEREND_KP*(CONTROLS.PRESSURE-SYSTEM.PRESSURE_TENSORS.Z)/GPa,1.0/3.0);
		
	
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)	
		{	
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
			{
				// Atom freezing doesn't work for a barostat!
				cout << "ERROR: Atoms cannot be frozen in an NPT ensemble." << endl;
				exit_run(0);
			}
			
			SYSTEM.COORDS[a1].X *= BEREND_ANI_MU.X;	
			SYSTEM.COORDS[a1].Y *= BEREND_ANI_MU.Y;
			SYSTEM.COORDS[a1].Z *= BEREND_ANI_MU.Z;	
		}
	
		SYSTEM.BOXDIM.X *= BEREND_ANI_MU.X;
		SYSTEM.BOXDIM.Y *= BEREND_ANI_MU.Y;
		SYSTEM.BOXDIM.Z *= BEREND_ANI_MU.Z;

		VOLUME_T = SYSTEM.BOXDIM.X*SYSTEM.BOXDIM.Y*SYSTEM.BOXDIM.Z;
		
		return;
	}	

	///////////////////////////////////////////////
	// Do MTTK thermostatting - NVT-MTK
	///////////////////////////////////////////////
	
	if(STYLE=="NVT-MTK" || STYLE=="NPT-MTK")
	{
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)	
		{	
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
				continue;
						
			SYSTEM.COORDS[a1].X -= 0.5*THERM_VELOC_0 * SYSTEM.VELOCITY[a1].X * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
			SYSTEM.COORDS[a1].Y -= 0.5*THERM_VELOC_0 * SYSTEM.VELOCITY[a1].Y * CONTROLS.DELTA_T * CONTROLS.DELTA_T;
			SYSTEM.COORDS[a1].Z -= 0.5*THERM_VELOC_0 * SYSTEM.VELOCITY[a1].Z * CONTROLS.DELTA_T * CONTROLS.DELTA_T;					
		}
		
		// Update thermostat position

		THERM_POSIT_T = THERM_POSIT_0 + THERM_VELOC_0 * CONTROLS.DELTA_T + 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * THERM_INERT_0;
	}

	///////////////////////////////////////////////
	// Do MTTK barostatting part - NPT-MTK
	///////////////////////////////////////////////
	
	if(STYLE=="NPT-MTK") //  HOOVER BAROSTAT: Not working properly yet!
	{
		// Update barostat position, compute scaling

		BAROS_POSIT_T = BAROS_POSIT_0 + BAROS_VELOC_0*CONTROLS.DELTA_T + 0.5*CONTROLS.DELTA_T*CONTROLS.DELTA_T*BAROS_FORCE_0/BAROS_INERT_W  - 0.5*CONTROLS.DELTA_T*CONTROLS.DELTA_T*BAROS_VELOC_0*THERM_VELOC_0;
		BAROS_SCALE   = exp(BAROS_POSIT_T - BAROS_POSIT_0);

		for(int a1=0;a1<SYSTEM.ATOMS;a1++)	
		{	
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
			{
				// Atom freezing doesn't work for a barostat!
				cout << "ERROR: Atoms cannot be frozen in an NPT ensemble." << endl;
				exit_run(0);
			}
			
			SYSTEM.COORDS[a1].X -= 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].X;
			SYSTEM.COORDS[a1].Y -= 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].Y;
			SYSTEM.COORDS[a1].Z -= 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].Z;	

			
			SYSTEM.COORDS[a1].X *= BAROS_SCALE;	
			SYSTEM.COORDS[a1].Y *= BAROS_SCALE;
			SYSTEM.COORDS[a1].Z *= BAROS_SCALE;				
		}
	
		SYSTEM.BOXDIM.X *= BAROS_SCALE;
		SYSTEM.BOXDIM.Y *= BAROS_SCALE;
		SYSTEM.BOXDIM.Z *= BAROS_SCALE;
	
		VOLUME_T = SYSTEM.BOXDIM.X*SYSTEM.BOXDIM.Y*SYSTEM.BOXDIM.Z;	
	}	
}

void CONSTRAINT::UPDATE_VELOCS_HALF_1(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{
	///////////////////////////////////////////////
	// Do the un-constrained update of velocities... this applies to all styles.
	///////////////////////////////////////////////
	
	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
		if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
			continue;
		
		SYSTEM.VELOCITY_ITER[a1].X = SYSTEM.VELOCITY[a1].X + 0.5 * SYSTEM.ACCEL[a1].X * CONTROLS.DELTA_T;
		SYSTEM.VELOCITY_ITER[a1].Y = SYSTEM.VELOCITY[a1].Y + 0.5 * SYSTEM.ACCEL[a1].Y * CONTROLS.DELTA_T;
		SYSTEM.VELOCITY_ITER[a1].Z = SYSTEM.VELOCITY[a1].Z + 0.5 * SYSTEM.ACCEL[a1].Z * CONTROLS.DELTA_T;
	}
	
	if(STYLE=="NVE" || STYLE=="NVT-SCALE" || STYLE=="NVT-BEREND" || STYLE=="NPT-BEREND" || STYLE=="NPT-BEREND-ANISO" )
		return;

	///////////////////////////////////////////////
	// Do MTTK thermostatting
	///////////////////////////////////////////////
	
	if(STYLE=="NVT-MTK" || STYLE=="NPT-MTK")
	{
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
				continue;
		
			SYSTEM.VELOCITY_ITER[a1].X -= 0.5*THERM_VELOC_0 * SYSTEM.VELOCITY[a1].X * CONTROLS.DELTA_T;
			SYSTEM.VELOCITY_ITER[a1].Y -= 0.5*THERM_VELOC_0 * SYSTEM.VELOCITY[a1].Y * CONTROLS.DELTA_T;
			SYSTEM.VELOCITY_ITER[a1].Z -= 0.5*THERM_VELOC_0 * SYSTEM.VELOCITY[a1].Z * CONTROLS.DELTA_T;
		}
	}

	///////////////////////////////////////////////
	// Do MTTK barostatting
	///////////////////////////////////////////////	
		
	if(STYLE=="NPT-MTK")
	{	
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{		
			SYSTEM.VELOCITY_ITER[a1].X -= 0.5 * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].X * CONTROLS.DELTA_T;
			SYSTEM.VELOCITY_ITER[a1].Y -= 0.5 * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].Y * CONTROLS.DELTA_T;
			SYSTEM.VELOCITY_ITER[a1].Z -= 0.5 * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].Z * CONTROLS.DELTA_T;

			SYSTEM.VELOCITY_ITER[a1].X *= BAROS_SCALE;
			SYSTEM.VELOCITY_ITER[a1].Y *= BAROS_SCALE;
			SYSTEM.VELOCITY_ITER[a1].Z *= BAROS_SCALE;
		}
	}	
}

void CONSTRAINT::UPDATE_VELOCS_HALF_2(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST)
{

	static double dV = CONTROLS.DELTA_T*VOLUME_T; //(VOLUME_T-VOLUME_0);
	static double dP = (SYSTEM.PRESSURE - CONTROLS.PRESSURE)/GPa;

	if(STYLE=="NVE" || STYLE=="NVT-SCALE" || STYLE == "NPT-BEREND" || STYLE == "NVT-BEREND" || STYLE=="NPT-BEREND-ANISO")
	{
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
				continue;

			SYSTEM.VELOCITY[a1].X = SYSTEM.VELOCITY_ITER[a1].X + 0.5*SYSTEM.ACCEL[a1].X * CONTROLS.DELTA_T;
			SYSTEM.VELOCITY[a1].Y = SYSTEM.VELOCITY_ITER[a1].Y + 0.5*SYSTEM.ACCEL[a1].Y * CONTROLS.DELTA_T;
			SYSTEM.VELOCITY[a1].Z = SYSTEM.VELOCITY_ITER[a1].Z + 0.5*SYSTEM.ACCEL[a1].Z * CONTROLS.DELTA_T;
		}
		
		///////////////////////////////////////////////
		// Do Berendsen thermostatting
		///////////////////////////////////////////////

		if(STYLE=="NVT-BEREND" || STYLE=="NPT-BEREND" || STYLE=="NPT-BEREND-ANISO")
		{
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
			{
				if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
					continue;
	
				BEREND_ETA = pow(1.0 + CONTROLS.DELTA_T/BEREND_TAU*(CONTROLS.TEMPERATURE/SYSTEM.TEMPERATURE-1.0),0.5);
		
				SYSTEM.VELOCITY[a1].X *= BEREND_ETA;
				SYSTEM.VELOCITY[a1].Y *= BEREND_ETA;
				SYSTEM.VELOCITY[a1].Z *= BEREND_ETA;
			}
		}
		
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
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
	else
	{
		if(STYLE=="NVT-MTK")	
		{
			KIN_ENER = kinetic_energy(SYSTEM, CONTROLS);
			THERM_INERT_T = ( 2.0 * KIN_ENER - N_DOF * Kb * CONTROLS.TEMPERATURE ) / (THERM_INERT_Q);
		}
		else // NPT-MTK
		{	
			KIN_ENER = kinetic_energy(SYSTEM, CONTROLS);	
			
			BAROS_FORCE_T = 3.0/(N_DOF) * 2.0 * KIN_ENER + dV*dP;
			THERM_INERT_T = ( 2.0 * KIN_ENER + BAROS_INERT_W * BAROS_VELOC_T * BAROS_VELOC_T - (N_DOF+1)*(Kb * CONTROLS.TEMPERATURE) ) / (THERM_INERT_Q);
		}

		for ( int itr = 0; itr < 10; itr++ ) // Iterative determination of velocities... See Martyna, Tobias, Klein JCP 101, 4177(1994) Appendix D.
		{
			if(STYLE=="NVT-MTK")
			{
				THERM_VELOC_T = THERM_VELOC_0 + (THERM_INERT_0 + THERM_INERT_T) * 0.5 * CONTROLS.DELTA_T;
				VSCALEH = 1.0 + 0.5 * CONTROLS.DELTA_T * THERM_VELOC_T;
			}
			else // NPT-MTK
			{
				// Need to solve for barostat veloc, since defined in terms of itself
				
				THERM_VELOC_T = THERM_VELOC_0 + (THERM_INERT_0 + THERM_INERT_T) * 0.5 * CONTROLS.DELTA_T;
				BAROS_VELOC_T = 0.5*CONTROLS.DELTA_T*(BAROS_VELOC_0 + BAROS_FORCE_0/BAROS_INERT_W - BAROS_VELOC_0*THERM_VELOC_0 + BAROS_FORCE_T/BAROS_INERT_W)/(1.0 + 0.5 * CONTROLS.DELTA_T * THERM_VELOC_T);

				VSCALEH = 1.0 + 0.5 * CONTROLS.DELTA_T * (THERM_VELOC_T + (2.0 + 3.0/(N_DOF))*BAROS_VELOC_T);
			}
			
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
			{
				if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
					continue;
				
				SYSTEM.VELOCITY_NEW[a1].X = (SYSTEM.VELOCITY_ITER[a1].X + 0.5*SYSTEM.ACCEL[a1].X*CONTROLS.DELTA_T) / VSCALEH;
				SYSTEM.VELOCITY_NEW[a1].Y = (SYSTEM.VELOCITY_ITER[a1].Y + 0.5*SYSTEM.ACCEL[a1].Y*CONTROLS.DELTA_T) / VSCALEH;
				SYSTEM.VELOCITY_NEW[a1].Z = (SYSTEM.VELOCITY_ITER[a1].Z + 0.5*SYSTEM.ACCEL[a1].Z*CONTROLS.DELTA_T) / VSCALEH;
			}

			if(STYLE=="NVT-MTK")
			{
				KIN_ENER = kinetic_energy(SYSTEM,"NEW", CONTROLS);
				THERM_INERT_T = ( 2.0 * KIN_ENER - N_DOF * Kb * CONTROLS.TEMPERATURE ) / (THERM_INERT_Q);
			}
			else // NPT-MTK
			{
				KIN_ENER = kinetic_energy(SYSTEM,"NEW", CONTROLS);
				BAROS_FORCE_T = 3.0/(N_DOF) * 2.0 * KIN_ENER + dV*dP;
				THERM_INERT_T = ( 2.0 * KIN_ENER + BAROS_INERT_W * BAROS_VELOC_T * BAROS_VELOC_T - (N_DOF+1)*(Kb * CONTROLS.TEMPERATURE) ) / (THERM_INERT_Q);
			}

		}
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
				continue;
			
			SYSTEM.VELOCITY[a1].X = SYSTEM.VELOCITY_NEW[a1].X;
			SYSTEM.VELOCITY[a1].Y = SYSTEM.VELOCITY_NEW[a1].Y;
			SYSTEM.VELOCITY[a1].Z = SYSTEM.VELOCITY_NEW[a1].Z;
			
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

void CONSTRAINT::SCALE_VELOCITIES(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{
	double vscale;
	
	SYSTEM.AVG_TEMPERATURE +=  SYSTEM.TEMPERATURE;
	SYSTEM.AVG_TEMPERATURE /= int(CONTROLS.FREQ_UPDATE_THERMOSTAT);
	vscale    = sqrt(CONTROLS.TEMPERATURE/SYSTEM.AVG_TEMPERATURE);

	if (RANK==0)
	{
		cout << "Average temperature     = " << SYSTEM.AVG_TEMPERATURE << endl;
		cout << "Velocity scaling factor = " << vscale << endl;	
	}

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
		SYSTEM.VELOCITY[a1].X *= vscale;
		SYSTEM.VELOCITY[a1].Y *= vscale;
		SYSTEM.VELOCITY[a1].Z *= vscale;
	}

	SYSTEM.AVG_TEMPERATURE = 0.0;
}

void CONSTRAINT::UPDATE_TEMPERATURE(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{
	static double Ktot;
	Ktot = kinetic_energy(SYSTEM, CONTROLS);	//calculate kinetic energy for scaling:

	SYSTEM.TEMPERATURE = 2.0 * Ktot / (N_DOF * Kb);

}

double CONSTRAINT::CONSERVED_QUANT (FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{
	static double THERM_KE, THERM_PE, BAROS_KE, BAROS_PE;
	static double TEMP_VOL;
	
	if(STYLE=="NVT-MTK")
	{
		THERM_KE = 0.5 * THERM_VELOC_T * THERM_VELOC_T * THERM_INERT_Q;
		THERM_PE = N_DOF * Kb * CONTROLS.TEMPERATURE * THERM_POSIT_T;
		
		return THERM_KE + THERM_PE;
	}
	else // NPT-MTK
	{
		TEMP_VOL = SYSTEM.BOXDIM.X*SYSTEM.BOXDIM.Y*SYSTEM.BOXDIM.Z;
	
		THERM_KE = 0.5 * THERM_VELOC_T * THERM_VELOC_T * THERM_INERT_Q;
		BAROS_KE = 0.5 * BAROS_VELOC_T * BAROS_VELOC_T * BAROS_INERT_W;
	
		THERM_PE = N_DOF * Kb * CONTROLS.TEMPERATURE * THERM_POSIT_T;
		BAROS_PE = CONTROLS.PRESSURE/GPa*TEMP_VOL;
		
		TEMP_VOL = pow(TEMP_VOL,1.0/3.0);
		
//		return THERM_KE/TEMP_VOL/TEMP_VOL + THERM_PE/TEMP_VOL + BAROS_KE + BAROS_PE;
		return THERM_KE + THERM_PE + BAROS_KE + BAROS_PE;
	}

	
}
