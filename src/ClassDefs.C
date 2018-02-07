#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes

#include "functions.h"
#include "util.h"

using namespace std;

NEIGHBORS::NEIGHBORS()		// Constructor 
{
	RCUT_PADDING  =  0.3;
	DISPLACEMENT  =  0.0;
	MAX_CUTOFF    =  0.0;
	MAX_CUTOFF_3B =  0.0;
	MAX_CUTOFF_4B =  0.0;
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
		LIST_4B       .resize(SYSTEM.ALL_ATOMS);
	}
}

void NEIGHBORS::INITIALIZE_MD(FRAME & SYSTEM)		// (overloaded) class constructor -- if no padding specified, default to 0.3
{
	if(USE)
	{
		INITIALIZE(SYSTEM);
		
		MAX_VEL = -1.0;

		for(int i=0; i<SYSTEM.ATOMS; i++)
		{
			CURR_VEL 
			      = sqrt(SYSTEM.VELOCITY[i].X*SYSTEM.VELOCITY[i].X 
			           + SYSTEM.VELOCITY[i].Y*SYSTEM.VELOCITY[i].Y 
	      	           + SYSTEM.VELOCITY[i].Z*SYSTEM.VELOCITY[i].Z);
			
			if(CURR_VEL > MAX_VEL)
				MAX_VEL = CURR_VEL;
		}
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

void NEIGHBORS::FIX_LAYERS(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
// Wrap atoms in ALL_COORDS into primitive cell prior to constructing replicas.	
{
	// Set the first NATOMS of ghost atoms to have the wrapped coordinates of the "real" coords
		
	for (int a=0; a<SYSTEM.ATOMS; a++) 
	{
		SYSTEM.WRAP_IDX[a].X = floor(SYSTEM.COORDS[a].X/SYSTEM.BOXDIM.X);
		SYSTEM.WRAP_IDX[a].Y = floor(SYSTEM.COORDS[a].Y/SYSTEM.BOXDIM.Y);
		SYSTEM.WRAP_IDX[a].Z = floor(SYSTEM.COORDS[a].Z/SYSTEM.BOXDIM.Z);
		
		SYSTEM.ALL_COORDS[a].X = SYSTEM.COORDS[a].X - SYSTEM.WRAP_IDX[a].X * SYSTEM.BOXDIM.X;
		SYSTEM.ALL_COORDS[a].Y = SYSTEM.COORDS[a].Y - SYSTEM.WRAP_IDX[a].Y * SYSTEM.BOXDIM.Y;
		SYSTEM.ALL_COORDS[a].Z = SYSTEM.COORDS[a].Z - SYSTEM.WRAP_IDX[a].Z * SYSTEM.BOXDIM.Z; 
	}
	
	// Build the surrounding "cell's" ghost atoms based on the first NATOMS ghost atoms
	
	if(CONTROLS.N_LAYERS>0 )
	{	
		int TEMP_IDX = SYSTEM.ATOMS;	

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
							SYSTEM.ALL_COORDS[TEMP_IDX].X = SYSTEM.ALL_COORDS[a1].X + n1 * SYSTEM.BOXDIM.X;
							SYSTEM.ALL_COORDS[TEMP_IDX].Y = SYSTEM.ALL_COORDS[a1].Y + n2 * SYSTEM.BOXDIM.Y;
							SYSTEM.ALL_COORDS[TEMP_IDX].Z = SYSTEM.ALL_COORDS[a1].Z + n3 * SYSTEM.BOXDIM.Z;
				
							if(SYSTEM.PARENT[TEMP_IDX] != a1)
							{
								cout << "ERROR: Wrong parent atom in found while updating layers" << endl;
								exit_run(0);
							}
				
							TEMP_IDX++;
						}
					}
				}
			}
		}
		
		if ( TEMP_IDX != SYSTEM.ALL_ATOMS ) 
		{
			printf("Error updating layers\n");
			exit(1);
		}
	}	
}

void NEIGHBORS::DO_UPDATE(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
// Choose algorithm based on system size including ghost atoms.
{
	FIX_LAYERS(SYSTEM, CONTROLS);

	if ( SYSTEM.ALL_ATOMS < 200 ) 
		DO_UPDATE_SMALL(SYSTEM, CONTROLS);
	else 
		DO_UPDATE_BIG(SYSTEM, CONTROLS);

	if ( CONTROLS.USE_3B_CHEBY ) 
	  UPDATE_3B_INTERACTION(SYSTEM, CONTROLS);
	
	if ( CONTROLS.USE_4B_CHEBY ) 
	  UPDATE_4B_INTERACTION(SYSTEM, CONTROLS);
}	
void NEIGHBORS::DO_UPDATE_SMALL(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)	
{

	XYZ RAB;
	
	if(!FIRST_CALL)	// Clear out the second dimension so we can start over again
	{
		for(int a1=0; a1<SYSTEM.ATOMS; a1++)
		{
			LIST          [a1].clear();
			LIST_UNORDERED[a1].clear();
			LIST_EWALD    [a1].clear();

			LIST_3B[a1].clear();
			LIST_4B[a1].clear();
		}
	}

	for(int a1=0; a1<SYSTEM.ATOMS; a1++)
	{
		for (int a2=0; a2<SYSTEM.ALL_ATOMS; a2++)
		{
			if(a2 == a1)
				continue;
			
			double rlen = get_dist(SYSTEM, RAB, a1, a2);
			
			if (rlen < MAX_CUTOFF + RCUT_PADDING)
				LIST_UNORDERED[a1].push_back(a2);	
			
			if ((SYSTEM.PARENT[a2]>=a1) && (a2 > a1))
			{
				if(rlen < (MAX_CUTOFF + RCUT_PADDING)) 			// Select atoms in neighbor list according to parents.
					LIST[a1].push_back(a2);		
				
				if (rlen < (EWALD_CUTOFF + RCUT_PADDING) )	
					LIST_EWALD[a1].push_back(a2);	
				
				if (rlen < MAX_CUTOFF_3B + RCUT_PADDING)	
					LIST_3B[a1].push_back(a2);	
				
				if (rlen < MAX_CUTOFF_4B + RCUT_PADDING)	
					LIST_4B[a1].push_back(a2);		
			}
		}
	}

	if(FIRST_CALL == false)
		SECOND_CALL = false;

	FIRST_CALL = false;	
}


void NEIGHBORS::DO_UPDATE_BIG(FRAME & SYSTEM, JOB_CONTROL & CONTROLS) 
// Order-N Neighbor list update with binning of particles.
{
	XYZ RAB;
	double rlen = 0;
	
	if(FIRST_CALL) // Set up the first dimension of the list 
	{
		LIST          .resize(SYSTEM.ATOMS);	
		LIST_EWALD    .resize(SYSTEM.ATOMS);	
		LIST_UNORDERED.resize(SYSTEM.ATOMS);	
		LIST_3B       .resize(SYSTEM.ATOMS);
		LIST_4B       .resize(SYSTEM.ATOMS);	
	}
	
	// Find maximum distance to search for neighbors.
	
	double SEARCH_DIST = MAX_CUTOFF;
	
	if ( EWALD_CUTOFF > SEARCH_DIST )
		SEARCH_DIST = EWALD_CUTOFF;
	
	if ( MAX_CUTOFF_3B > SEARCH_DIST ) 
		SEARCH_DIST = MAX_CUTOFF_3B;
	
	SEARCH_DIST += RCUT_PADDING;

	XYZ_INT NBINS;
	
	NBINS.X = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.X /SEARCH_DIST) + 2;
	NBINS.Y = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.Y /SEARCH_DIST) + 2;
	NBINS.Z = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.Z /SEARCH_DIST) + 2;
	
	int TOTAL_BINS = NBINS.X * NBINS.Y * NBINS.Z;

	vector<vector<int> > BIN(TOTAL_BINS);
	
	for (int i=0; i<TOTAL_BINS; i++) 
		vector<int>().swap(BIN[i]);
	
	XYZ_INT BIN_IDX;
	
	int FULLNESS = 0;
	
	for ( int a1 = 0; a1 < SYSTEM.ALL_ATOMS; a1++ ) 
	{				
		BIN_IDX.X = floor( (SYSTEM.ALL_COORDS[a1].X + SYSTEM.BOXDIM.X * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX.Y = floor( (SYSTEM.ALL_COORDS[a1].Y + SYSTEM.BOXDIM.Y * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX.Z = floor( (SYSTEM.ALL_COORDS[a1].Z + SYSTEM.BOXDIM.Z * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;

		if ( BIN_IDX.X < 0 || BIN_IDX.Y < 0 || BIN_IDX.Z < 0 ) 
		{
			cout << "Error: negative binning BIN_IDX for atom a1 = " << a1 << endl;
			cout << "Check box lengths in .xyz* file." << endl;
			exit(1);
		}
		
		// Calculate bin BIN_IDX of the atom.
		int ibin = BIN_IDX.X + BIN_IDX.Y * NBINS.X + BIN_IDX.Z * NBINS.X * NBINS.Y;

		if ( ibin >= TOTAL_BINS ) 
		{
			cout << "Error: binning BIN_IDX out of range\n";
			cout << "BIN_IDX.X = " << BIN_IDX.X << "BIN_IDX.Y = " << BIN_IDX.Y << "BIN_IDX.Z = " << BIN_IDX.Z << endl;
			exit(1);
		}
		
		// Push the atom into the bin 
		BIN[ibin].push_back(a1);
		FULLNESS++;
		
			
	}

	for(int a1=0; a1<SYSTEM.ATOMS; a1++)
	{
		if(!FIRST_CALL)
		{
			LIST          [a1].clear();
			LIST_UNORDERED[a1].clear();
			LIST_EWALD    [a1].clear();
			LIST_3B       [a1].clear();
			LIST_4B       [a1].clear();
		}
		
		XYZ_INT BIN_IDX_a1;

		BIN_IDX_a1.X = floor( (SYSTEM.ALL_COORDS[a1].X + SYSTEM.BOXDIM.X * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX_a1.Y = floor( (SYSTEM.ALL_COORDS[a1].Y + SYSTEM.BOXDIM.Y * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX_a1.Z = floor( (SYSTEM.ALL_COORDS[a1].Z + SYSTEM.BOXDIM.Z * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;

		if ( BIN_IDX_a1.X < 1 || BIN_IDX_a1.Y < 1 || BIN_IDX_a1.Z < 1 ) 
		{
			cout << "Error: bad binning BIN_IDX\n";
			cout << "BIN_IDX.X = " << BIN_IDX_a1.X << "BIN_IDX.Y = " << BIN_IDX_a1.Y << "BIN_IDX.Z = " << BIN_IDX_a1.Z << endl;
			exit(1);
		}

		// Loop over relevant bins only, not all atoms.
		
		int ibin, a2, a2end;
		
		for (int i=BIN_IDX_a1.X-1; i<= BIN_IDX_a1.X+1; i++)	//BIN_IDX_a1.X 
		{
			for (int j=BIN_IDX_a1.Y-1; j<=BIN_IDX_a1.Y+1; j++ ) //BIN_IDX_a1.Y
			{
				for (int k=BIN_IDX_a1.Z-1; k<=BIN_IDX_a1.Z+1; k++ ) // BIN_IDX_a1.Z
				{
					ibin = i + j * NBINS.X + k * NBINS.X * NBINS.Y;


					if (ibin >= TOTAL_BINS) 
					{
						cout << "Error: binning BIN_IDX out of range\n";
						cout << "BIN_IDX.X = " << i << "BIN_IDX.Y = " << j << "BIN_IDX.Z = " << k << endl;
						exit(1);
					}

					a2end = BIN[ibin].size();
					
					// Check all atoms in the bin.
 
					for (int a2idx=0; a2idx<a2end; a2idx++) 
					{
						a2 = BIN[ibin][a2idx];

						if ( a2 == a1 ) 
							continue;

						rlen = get_dist(SYSTEM, RAB, a1, a2);

						if (rlen < MAX_CUTOFF + RCUT_PADDING)		
							LIST_UNORDERED[a1].push_back(a2);	
						
						if ( a1 <= SYSTEM.PARENT[a2] ) 
						{
							if (rlen < (MAX_CUTOFF + RCUT_PADDING) )		
								LIST[a1].push_back(a2);		

							if (rlen < (EWALD_CUTOFF + RCUT_PADDING) )
								LIST_EWALD[a1].push_back(a2);		

							if(rlen < MAX_CUTOFF_3B + RCUT_PADDING)	
								LIST_3B[a1].push_back(a2);
							
							if(rlen < MAX_CUTOFF_4B + RCUT_PADDING)	
								LIST_4B[a1].push_back(a2);	
						}	
					}
				}
			}
		}
		

	}

	if(FIRST_CALL == false)
		SECOND_CALL = false;

	FIRST_CALL = false;	
	
}

void NEIGHBORS::UPDATE_LIST(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{	
	if(FIRST_CALL)	// Then start from scratch by cycling through all atom (and layer atom) pairs
	{		
		if (!USE) 
			RCUT_PADDING = 1.0e10;
		
		if(RANK==0)
		{
			cout << "NEIGHBOR LIST MAX CUTOFF IS: " << fixed << setprecision(5) << MAX_CUTOFF << endl;
			cout << "USING PADDING: " << fixed << setprecision(3) << RCUT_PADDING << endl;
		}
		
		DO_UPDATE(SYSTEM, CONTROLS);
		FIRST_CALL = false;	
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
			
			if(RANK == 0)
				cout << "RANK: " << RANK << " RESET RCUT_PADDING TO: " << fixed << setprecision(10) << RCUT_PADDING << endl;
			
			DO_UPDATE(SYSTEM, CONTROLS);
			SECOND_CALL = false;
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
	else if(IN_STYLE=="NVE")					// Uses velocity scaling to thermostat
		STYLE = IN_STYLE;	
	else if(CONTROLS.FREQ_UPDATE_THERMOSTAT > -1.0)	// Trivial velocity scaling
		STYLE = "NVT-SCALE";
	else if(IN_STYLE=="LMP-NVE" || IN_STYLE=="LMP-NVT" || IN_STYLE=="LMP-NPT" ||
	 (CONTROLS.ENSEMBLE == "LMP-MIN-BOX-ISO" || CONTROLS.ENSEMBLE == "LMP-MIN-BOX-ANISO" || CONTROLS.ENSEMBLE == "LMP-MIN-BOX-TRI" || CONTROLS.ENSEMBLE == "LMP-MIN"))
		cout << "	...Configuring constraints for a " << IN_STYLE << " simulation." << endl;
	else
	{
		cout << "ERROR: UNKNOWN CONSTRAINT STYLE" << endl;
		exit_run(0);
	}
	
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

void CONSTRAINT::WRITE(ofstream &output)
// Write parameters for restart.
{
	output << THERM_POSIT_T << endl;
	output << THERM_VELOC_T << endl;
	output << THERM_INERT_T << endl;
	output << THERM_INERT_Q << endl;
	output << THERM_POSIT_0 << endl;
	output << THERM_INERT_0 << endl;
	output << THERM_VELOC_0 << endl;

	output << BAROS_POSIT_T << endl;
	output << BAROS_VELOC_T << endl;
	output << BAROS_FORCE_T << endl;
	output << BAROS_INERT_W << endl;
	output << BAROS_POSIT_0 << endl;
	output << BAROS_FORCE_0 << endl;
	output << BAROS_VELOC_0 << endl;

	output << BAROS_SCALE << endl;
	output << VOLUME_0 << endl;
	output << VOLUME_T << endl;

	output << BEREND_MU << endl;
	output << BEREND_ANI_MU.X << endl;
	output << BEREND_ANI_MU.Y << endl;
	output << BEREND_ANI_MU.Z << endl;
	output << BEREND_KP << endl;

	output << BEREND_ETA << endl;
	output << BEREND_TAU << endl;

	output << TIME << endl;
	output << TIME_BARO << endl;
	output << N_DOF << endl;
	output << VSCALEH << endl;
	output << KIN_ENER << endl;
}


void CONSTRAINT::READ(ifstream &input)
// Write parameters for restart.
{
	input >> THERM_POSIT_T;
	input >> THERM_VELOC_T;
	input >> THERM_INERT_T;
	input >> THERM_INERT_Q;
	input >> THERM_POSIT_0;
	input >> THERM_INERT_0;
	input >> THERM_VELOC_0;

	input >> BAROS_POSIT_T;
	input >> BAROS_VELOC_T;
	input >> BAROS_FORCE_T;
	input >> BAROS_INERT_W;
	input >> BAROS_POSIT_0;
	input >> BAROS_FORCE_0;
	input >> BAROS_VELOC_0;

	input >> BAROS_SCALE ;
	input >> VOLUME_0 ;
	input >> VOLUME_T ;

	input >> BEREND_MU ;
	input >> BEREND_ANI_MU.X ;
	input >> BEREND_ANI_MU.Y ;
	input >> BEREND_ANI_MU.Z ;
	input >> BEREND_KP ;

	input >> BEREND_ETA ;
	input >> BEREND_TAU ;

	input >> TIME ;
	input >> TIME_BARO ;
	input >> N_DOF ;
	input >> VSCALEH ;
	input >> KIN_ENER;
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
	
	///////////////////////////////////////////////
	// Refresh ghost atom positions
	///////////////////////////////////////////////

	// Set the first NATOMS of ghost atoms to have the coordinates of the "real" coords

	for (int a=0; a<SYSTEM.ATOMS; a++) 
	{
		SYSTEM.ALL_COORDS[a].X = SYSTEM.COORDS[a].X - SYSTEM.WRAP_IDX[a].X * SYSTEM.BOXDIM.X;
		SYSTEM.ALL_COORDS[a].Y = SYSTEM.COORDS[a].Y - SYSTEM.WRAP_IDX[a].Y * SYSTEM.BOXDIM.Y;
		SYSTEM.ALL_COORDS[a].Z = SYSTEM.COORDS[a].Z - SYSTEM.WRAP_IDX[a].Z * SYSTEM.BOXDIM.Z; 
	}
	
	// Build the surrounding "cell's" ghost atoms based on the first NATOMS ghost atoms
	
	if(CONTROLS.N_LAYERS>0 )
	{	
		int TEMP_IDX = SYSTEM.ATOMS;	

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
							SYSTEM.ALL_COORDS[TEMP_IDX].X = SYSTEM.ALL_COORDS[a1].X + n1 * SYSTEM.BOXDIM.X;
							SYSTEM.ALL_COORDS[TEMP_IDX].Y = SYSTEM.ALL_COORDS[a1].Y + n2 * SYSTEM.BOXDIM.Y;
							SYSTEM.ALL_COORDS[TEMP_IDX].Z = SYSTEM.ALL_COORDS[a1].Z + n3 * SYSTEM.BOXDIM.Z;
				
							if(SYSTEM.PARENT[TEMP_IDX] != a1)
							{
								cout << "ERROR: Wrong parent atom in found while updating layers" << endl;
								exit_run(0);
							}
				
							TEMP_IDX++;
						}
					}
				}
			}
		}
		
		if ( TEMP_IDX != SYSTEM.ALL_ATOMS ) 
		{
			printf("Error updating layers\n");
			exit(1);
		}
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

		return  THERM_KE + THERM_PE;
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

void NEIGHBORS::UPDATE_3B_INTERACTION(FRAME & SYSTEM, JOB_CONTROL &CONTROLS) 
// Build a list of all 3-body interactions.  This "flat" list parallelizes much
// more efficiently than a nested neighbor list loop.
{
	XYZ RAB;  
	INTERACTION_3B inter;

	LIST_3B_INT.clear();
	for ( int i = 0; i < SYSTEM.ATOMS; i++ ) {
		int ai = i;
		for ( int j = 0; j < LIST_3B[i].size(); j++ ) {
			int aj = LIST_3B[i][j];
			for ( int k = 0; k < LIST_3B[i].size(); k++ ) {
				int ak = LIST_3B[i][k];

				if ( aj == ak || SYSTEM.PARENT[aj] > SYSTEM.PARENT[ak] ) 
					continue;

				// The j-k list is possibly outside of the cutoff, so test it here.
				double rlen = get_dist(SYSTEM, RAB, aj, ak);
	
				if ( rlen < MAX_CUTOFF_3B + RCUT_PADDING ) {
					inter.a1 = ai;
					inter.a2 = aj;
					inter.a3 = ak;
	  
					LIST_3B_INT.push_back(inter);
				}
			}
		}
	}
}

void NEIGHBORS::UPDATE_4B_INTERACTION(FRAME & SYSTEM, JOB_CONTROL &CONTROLS) 
// Build a list of all 4-body interactions.  This "flat" list parallelizes much
// more efficiently than a nested neighbor list loop.
{
	XYZ RAB;  
	INTERACTION_4B inter;
	LIST_4B_INT.clear();
	int ai, aj, ak, al;
	
	for (int i=0; i<SYSTEM.ATOMS; i++)	// Loop over all real atoms
	{
		ai = i;
		
		for (int j=0; j<LIST_3B[i].size(); j++) // Loop over all neighbors of i to get atom j
		{
			aj = LIST_4B[i][j];
			
			for (int k=0; k<LIST_3B[i].size(); k++) // Loop over all neighbors of i to get atom k
			{
				ak = LIST_4B[i][k];
				
				for (int l=0; l<LIST_3B[i].size(); l++) // Loop over all neighbors of i to get atom l
				{
					al = LIST_4B[i][l];
				
					// Check that this is a valid quadruplet

					if (aj == ak || aj == al || ak == al)
						continue;
					
					if(SYSTEM.PARENT[aj] > SYSTEM.PARENT[ak] || SYSTEM.PARENT[aj] > SYSTEM.PARENT[al] || SYSTEM.PARENT[ak] > SYSTEM.PARENT[al])
						continue;
					
					// We know ij, ik, and il distances are within the allowed cutoffs, but we still need to check jk, jl, and kl

					if(get_dist(SYSTEM, RAB, aj, ak) >  MAX_CUTOFF_4B + RCUT_PADDING)
						continue;
					if(get_dist(SYSTEM, RAB, aj, al) >  MAX_CUTOFF_4B + RCUT_PADDING)
						continue;
					if(get_dist(SYSTEM, RAB, ak, al) >  MAX_CUTOFF_4B + RCUT_PADDING)
						continue;
					
					inter.a1 = ai;
					inter.a2 = aj;
					inter.a3 = ak;
					inter.a4 = al;
					
					LIST_4B_INT.push_back(inter);
				}
			}
		}
	}
}

void THERMO_AVG::WRITE(ofstream &fout)
// Write out thermodynamic average properties.
{
	fout << TEMP_SUM << endl;
	fout << PRESS_SUM << endl;
	fout << STRESS_TENSOR_SUM.X << endl;
	fout << STRESS_TENSOR_SUM.Y << endl;
	fout << STRESS_TENSOR_SUM.Z << endl;
}


void THERMO_AVG::READ(ifstream &fin)
// Read in thermodynamic average properties.
{
	fin >> TEMP_SUM ;
	fin >> PRESS_SUM ;
	fin >> STRESS_TENSOR_SUM.X ;
	fin >> STRESS_TENSOR_SUM.Y ;
	fin >> STRESS_TENSOR_SUM.Z ;
}


