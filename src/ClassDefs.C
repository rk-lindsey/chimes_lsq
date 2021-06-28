#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes

#include "functions.h"
#include "util.h"
#include "Cheby.h"

using namespace std;

NEIGHBORS::NEIGHBORS()
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
	MAX_COORD_STEP = 0.0 ;
	CURR_VEL      =  0.0;
	
	// New from Larry
	EWALD_CUTOFF  =  0.0;
	UPDATE_FREQ   =  30.0;
	SAFETY        =  1.01;
	
	// New for triclinic support
	UPDATE_WITH_BIG = true;

	PERM_SCALE.resize(MAX_BODIEDNESS+1) ;
	for ( int j = 0 ; j < MAX_BODIEDNESS + 1 ; j++ ) {
		 PERM_SCALE[j] = 1.0 ;
	}
	
}
NEIGHBORS::~NEIGHBORS(){}	// Deconstructor

void NEIGHBORS::INITIALIZE(FRAME & SYSTEM)		// (overloaded) class constructor -- if no padding specified, default to 0.3
{
	 LIST          .resize(SYSTEM.ATOMS);
	 LIST_EWALD    .resize(SYSTEM.ATOMS);
	 LIST_UNORDERED.resize(SYSTEM.ATOMS);
	 LIST_3B       .resize(SYSTEM.ALL_ATOMS);
	 LIST_4B       .resize(SYSTEM.ALL_ATOMS);
		
	 // UPDATE_WITH_BIG is true by default
	 // If it is already false at this point, it is because the user requested it
	 // a small update may be requested for cell vectors that don't obey our requirements (all positive)

	 double max_cutoff = MAX_ALL_CUTOFFS() ;
	 if (   SYSTEM.BOXDIM.EXTENT_X <= max_cutoff
					|| SYSTEM.BOXDIM.EXTENT_Y <= max_cutoff
					|| SYSTEM.BOXDIM.EXTENT_Z <= max_cutoff ) {
			if ( RANK == 0 ) {
				 cout << "Warning:  system size <= cutoff.  Using explicit permutation / small system neighbor algorithm.\n" ;
				 cout << "Warning:  many-body interactions will be substantially slower.\n" ;
			}
			UPDATE_WITH_BIG = false;
	 }

}

void NEIGHBORS::INITIALIZE_MD(FRAME & SYSTEM, JOB_CONTROL &CONTROLS)		// (overloaded) class constructor -- if no padding specified, default to 0.3
{
	INITIALIZE(SYSTEM);
		
	MAX_VEL = -1.0;

	if(USE)
	{
		for(int i=0; i<SYSTEM.ATOMS; i++)
		{
			CURR_VEL 
			      = sqrt(SYSTEM.VELOCITY[i].X*SYSTEM.VELOCITY[i].X 
			           + SYSTEM.VELOCITY[i].Y*SYSTEM.VELOCITY[i].Y 
			           + SYSTEM.VELOCITY[i].Z*SYSTEM.VELOCITY[i].Z);
			
			if(CURR_VEL > MAX_VEL)
				MAX_VEL = CURR_VEL;
		}
		MAX_COORD_STEP = MAX_VEL * CONTROLS.DELTA_T ;
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
	// Build the surrounding "cell's" ghost atoms based on the first NATOMS ghost atoms
	SYSTEM.update_ghost(CONTROLS.N_LAYERS, true) ;
}

void NEIGHBORS::DO_UPDATE(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
// Choose algorithm based on system size including ghost atoms.
{
	FIX_LAYERS(SYSTEM, CONTROLS);
		
	if (UPDATE_WITH_BIG && USE)
		DO_UPDATE_BIG(SYSTEM, CONTROLS);
	else
		 DO_UPDATE_SMALL(SYSTEM, CONTROLS);

	if ( CONTROLS.USE_3B_CHEBY ) 
	  UPDATE_3B_INTERACTION(SYSTEM, CONTROLS);
	
	if ( CONTROLS.USE_4B_CHEBY ) 
	  UPDATE_4B_INTERACTION(SYSTEM, CONTROLS);
}

double NEIGHBORS::MAX_ALL_CUTOFFS()
// Returns the maximum of all cutoff values
{
	 double val = -1.0 ;

	 if ( MAX_CUTOFF > val )    val = MAX_CUTOFF ;
	 if ( MAX_CUTOFF_3B > val ) val = MAX_CUTOFF_3B ;
	 if ( MAX_CUTOFF_4B > val ) val = MAX_CUTOFF_4B ;
	 if ( EWALD_CUTOFF > val )  val = EWALD_CUTOFF ;
	 return val ;
}

void NEIGHBORS::DO_UPDATE_SMALL(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
// This implements the "unordered" neighbor list convention where i < j and i > j are
// explicitly included.  This is necessary for correct evalution of interactions between
// atom i and its self-image when the cutoff is greater than or equal to the box size
// (small box size limit).
{

	XYZ RAB;

	// Set scaling factors for permutations of unordered interactions.
	// PERM_SCALE[n] is the scaling factor for the n-body interaction.
	PERM_SCALE[0] = 1.0 ;
	PERM_SCALE[1] = 1.0 ;	
	for ( int j = 2 ; j < PERM_SCALE.size() ; j++ ) {
		PERM_SCALE[j] = PERM_SCALE[j-1] / j ;
	}
	
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

			if(rlen < (MAX_CUTOFF + RCUT_PADDING)) 			// Select atoms in neighbor list according to parents.
				LIST[a1].push_back(a2);	

			if (rlen < MAX_CUTOFF_3B + RCUT_PADDING)	
				LIST_3B[a1].push_back(a2);

			if (rlen < MAX_CUTOFF_4B + RCUT_PADDING)	
				LIST_4B[a1].push_back(a2);
			
			if (rlen < (EWALD_CUTOFF + RCUT_PADDING) )	
				 LIST_EWALD[a1].push_back(a2);	
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

	for ( int j = 0 ; j < PERM_SCALE.size() ; j++ ) {
		 PERM_SCALE[j] = 1.0 ;
	}
	
	// Find maximum distance to search for neighbors.
	
	double SEARCH_DIST = MAX_CUTOFF;
	
	if ( EWALD_CUTOFF > SEARCH_DIST )
		SEARCH_DIST = EWALD_CUTOFF;
	
	if ( MAX_CUTOFF_3B > SEARCH_DIST ) 
		SEARCH_DIST = MAX_CUTOFF_3B;
		
	if ( MAX_CUTOFF_4B > SEARCH_DIST ) 
		SEARCH_DIST = MAX_CUTOFF_4B;		
	
	SEARCH_DIST += RCUT_PADDING;

	XYZ_INT NBINS;
	
	NBINS.X = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.EXTENT_X/SEARCH_DIST) + 2;
	NBINS.Y = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.EXTENT_Y/SEARCH_DIST) + 2;
	NBINS.Z = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.EXTENT_Z/SEARCH_DIST) + 2;	
	
	
	/* RKL - no longer used - 082319
	
	NBINS.X = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.X /SEARCH_DIST) + 2;
	NBINS.Y = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.Y /SEARCH_DIST) + 2;
	NBINS.Z = ceil((2 * CONTROLS.N_LAYERS+1) * SYSTEM.BOXDIM.Z /SEARCH_DIST) + 2;
	*/
	
	int TOTAL_BINS = NBINS.X * NBINS.Y * NBINS.Z;

	vector<vector<int> > BIN(TOTAL_BINS);
	
	for (int i=0; i<TOTAL_BINS; i++) 
		vector<int>().swap(BIN[i]);
	
	XYZ_INT BIN_IDX;
	
	int FULLNESS = 0;
	
	for ( int a1 = 0; a1 < SYSTEM.ALL_ATOMS; a1++ ) 
	{	
		
		BIN_IDX.X = floor( (SYSTEM.ALL_COORDS[a1].X + SYSTEM.BOXDIM.EXTENT_X * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX.Y = floor( (SYSTEM.ALL_COORDS[a1].Y + SYSTEM.BOXDIM.EXTENT_Y * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX.Z = floor( (SYSTEM.ALL_COORDS[a1].Z + SYSTEM.BOXDIM.EXTENT_Z * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;

		/* RKL - no longer used - 082319
					
		BIN_IDX.X = floor( (SYSTEM.ALL_COORDS[a1].X + SYSTEM.BOXDIM.X * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX.Y = floor( (SYSTEM.ALL_COORDS[a1].Y + SYSTEM.BOXDIM.Y * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX.Z = floor( (SYSTEM.ALL_COORDS[a1].Z + SYSTEM.BOXDIM.Z * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		*/

		if ( BIN_IDX.X < 0 || BIN_IDX.Y < 0 || BIN_IDX.Z < 0 ) 
		{
			cout << "Error: negative binning BIN_IDX for atom a1 = " << a1 << endl;
			cout << "Atom has coordinates: " << SYSTEM.ALL_COORDS[a1].X << ", " << SYSTEM.ALL_COORDS[a1].Y << ", " << SYSTEM.ALL_COORDS[a1].Z << endl;
			cout << "Atom has bins:        " << BIN_IDX.X << " " << BIN_IDX.Y << " " << BIN_IDX.Z << endl;
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
		
		BIN_IDX_a1.X = floor( (SYSTEM.ALL_COORDS[a1].X + SYSTEM.BOXDIM.EXTENT_X * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX_a1.Y = floor( (SYSTEM.ALL_COORDS[a1].Y + SYSTEM.BOXDIM.EXTENT_Y * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX_a1.Z = floor( (SYSTEM.ALL_COORDS[a1].Z + SYSTEM.BOXDIM.EXTENT_Z * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		
		/* RKL - no longer used - 082319

		BIN_IDX_a1.X = floor( (SYSTEM.ALL_COORDS[a1].X + SYSTEM.BOXDIM.X * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX_a1.Y = floor( (SYSTEM.ALL_COORDS[a1].Y + SYSTEM.BOXDIM.Y * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		BIN_IDX_a1.Z = floor( (SYSTEM.ALL_COORDS[a1].Z + SYSTEM.BOXDIM.Z * CONTROLS.N_LAYERS) / SEARCH_DIST ) + 1;
		*/
		
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
	else
	{
		// Note: We keep track of the max velocity in the part of splines_md.C 
		// where the first half of thermostatting is applied

		if( (SECOND_CALL))	// Update the cutoff padding
		{
			if(USE) {
				// RCUT_PADDING = MAX_VEL * UPDATE_FREQ * CONTROLS.DELTA_T;	// should give a distance in AA
				RCUT_PADDING = MAX_COORD_STEP * UPDATE_FREQ ;
			}
			
			if(RANK == 0)
				cout << "RANK: " << RANK << " RESET RCUT_PADDING TO: " << fixed << setprecision(10) << RCUT_PADDING << endl;
			
			DO_UPDATE(SYSTEM, CONTROLS);
			SECOND_CALL = false;
		}
		else
		{
			// DISPLACEMENT += MAX_VEL * CONTROLS.DELTA_T;
			// For NPT dynamics, the time derivative of the position is not equal to the velocity,
			// so displacements should be based on direct evaluation of position steps.
			//
			DISPLACEMENT += MAX_COORD_STEP ;

			if(DISPLACEMENT>0.5*RCUT_PADDING)
			{
				if (USE) {
					// RCUT_PADDING = MAX_VEL * UPDATE_FREQ * CONTROLS.DELTA_T;	// Update padding in case max_vel changed. (LEF).
					RCUT_PADDING = MAX_COORD_STEP * UPDATE_FREQ ;
				}
				
				DO_UPDATE(SYSTEM, CONTROLS);
				
				DISPLACEMENT = 0.0 ;
				MAX_VEL      = 0.0 ;
				MAX_COORD_STEP = 0.0 ;

#if VERBOSITY >= 1				
				if(RANK == 0)
					cout << " Updating neighbor list on step: " << CONTROLS.STEP << ", with padding: " << fixed << setprecision(3) << RCUT_PADDING <<  endl;
#endif
			}
		}
	}
}

// void NEIGHBORS::UPDATE_LIST(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, bool FORCE)
// {
// 	if(FIRST_CALL || SECOND_CALL)
// 		UPDATE_LIST(SYSTEM, CONTROLS);
// 	else
// 	{
// 		DISPLACEMENT += MAX_VEL * CONTROLS.DELTA_T;
		
// 		DO_UPDATE(SYSTEM, CONTROLS);
	
// 		if(RANK == 0)
// 			cout << "RANK: " << RANK << " FORCING UPDATING ON STEP: " << CONTROLS.STEP << ", WITH PADDING: " << fixed << setprecision(3) << RCUT_PADDING <<  endl;
// 	}
// }



CONSTRAINT::CONSTRAINT()
{
	BEREND_MU = 1.0 ;

	THERM_POSIT_T = 0.0 ;
	THERM_VELOC_T = 0.0 ;
	THERM_FORCE_T = 0.0 ;
	THERM_INERT_Q = 1.0 ;
	THERM_POSIT_0 = 0.0 ;
	THERM_FORCE_0 = 0.0 ;
	THERM_VELOC_0 = 0.0 ;
			
	BAROS_POSIT_T = 0.0 ;
	BAROS_VELOC_T = 0.0 ;
	BAROS_FORCE_T = 0.0 ;
	BAROS_INERT_W = 0.0 ;
	BAROS_POSIT_0 = 0.0 ;
	BAROS_FORCE_0 = 0.0 ;
	BAROS_VELOC_0 = 0.0 ;	
		
	 //Barostat variables

	BAROS_SCALE = 1.0 ;
	VOLUME_0 = 1.0 ;
	VOLUME_T = 1.0 ;
		
	 // Berendsen barostat variables
		
	BEREND_MU = 1.0 ;		
	BEREND_ANI_MU.X = 0.0 ;
	BEREND_ANI_MU.Y = 0.0 ;
	BEREND_ANI_MU.Z = 0.0 ;			
	BEREND_KP = 1.0 ;	
		
	BEREND_ETA = 1.0 ;	
	BEREND_TAU = 1.0 ;	
		
	 // Thermostat variables

	TIME = 0.0 ;		
	TIME_BARO = 0.0 ;
	N_DOF = 0 ;		
	VSCALEH = 1.0 ;		
	KIN_ENER = 1.0 ;	

}	// Constructor
CONSTRAINT::~CONSTRAINT(){}	// Deconstructor

void CONSTRAINT::INITIALIZE(string IN_STYLE, JOB_CONTROL & CONTROLS, int ATOMS)
{
	STYLE = IN_STYLE;
	
	if(IN_STYLE=="NPT-MTK")				
		STYLE = IN_STYLE;
	else if(IN_STYLE=="NPT-BEREND")			// Uses position scaling to barostat
		STYLE = IN_STYLE;
	else if(IN_STYLE=="NPT-BEREND-ANISO")		// Uses position scaling to barostat anisotropically, but expects alpha = beta = gamma = 90 deg
		STYLE = IN_STYLE;	
	else if(IN_STYLE=="NVT-BEREND")			// Uses velocity scaling to thermostat
		STYLE = IN_STYLE;
	else if(CONTROLS.USE_HOOVER_THRMOSTAT)		// Uses MTK Thermostat
		STYLE = "NVT-MTK";
	else if(IN_STYLE=="NVE")			// Uses velocity scaling to thermostat
		STYLE = IN_STYLE;	
	else if(CONTROLS.FREQ_UPDATE_THERMOSTAT > -1.0)	// Trivial velocity scaling
		STYLE = "NVT-SCALE";
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
	THERM_FORCE_T = 0;
	
	BAROS_POSIT_T = 0;
	BAROS_VELOC_T = 0;
	BAROS_FORCE_T = 0;
	
	THERM_POSIT_0 = 0;
	THERM_VELOC_0 = 0;
	THERM_FORCE_0 = 0;
	
	BAROS_POSIT_0 = 0;
	BAROS_VELOC_0 = 0;
	BAROS_FORCE_0 = 0;
	
	VOLUME_0 = 0;
	VOLUME_T = 0;
	
	BAROS_SCALE = 1;
	
}

XYZ CONSTRAINT::CENTER_OF_MASS(const FRAME &SYSTEM)
// Returns the center of mass of the system.
{
	XYZ com{0.0,0.0,0.0} ;

	double total_mass = 0.0 ;
	for ( int a = 0 ; a < SYSTEM.ATOMS ; a++ ) {
		com.X += SYSTEM.COORDS[a].X * SYSTEM.MASS[a] ;
		com.Y += SYSTEM.COORDS[a].Y * SYSTEM.MASS[a] ;
		com.Z += SYSTEM.COORDS[a].Z * SYSTEM.MASS[a] ;

		total_mass += SYSTEM.MASS[a] ;
	}
	com.X /= total_mass ;
	com.Y /= total_mass ;
	com.Z /= total_mass ;

	return(com) ;
}

		
void CONSTRAINT::INIT_VEL (FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{
	// Use box Muller to initialize velocities

	double x1, x2 , y1 ,y2;
	double sigma;
	int    counter = 0;
	       
	srand(CONTROLS.SEED);
	       
	for(int a=0; a<SYSTEM.ATOMS; a++)
	{
	       // Don't account for frozen atoms
		       
	       if((CONTROLS.FREEZE_IDX_START != -1) && ((a<CONTROLS.FREEZE_IDX_START) || (a>CONTROLS.FREEZE_IDX_STOP)))
	       {
		 SYSTEM.VELOCITY[a].X = 0;
		 SYSTEM.VELOCITY[a].Y = 0;
		 SYSTEM.VELOCITY[a].Z = 0;
			       
		 continue;
	       }
		       
	       sigma = sqrt(CONTROLS.TEMPERATURE * Kb / SYSTEM.MASS[a] );

	       // Do for x...
		       
	       x1    = double(rand())/double(RAND_MAX);
	       x2    = double(rand())/double(RAND_MAX);
		       
	       if ( (x1 < 0.0 || x1 > 1.0) || ( x2 < 0.0 || x2 > 1.0 ) )
	       {
		 cout << "Bad random variable" << endl;
		 exit_run(1);
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
		 exit_run(1);
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
		 exit_run(1);
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

void CONSTRAINT::CHECK_VEL(FRAME & SYSTEM, JOB_CONTROL & CONTROLS)
{
	  XYZ    TEMP_VEL;
	  double TEMP_MASS;
		
	  TEMP_VEL.X = 0;
	  TEMP_VEL.Y = 0;
	  TEMP_VEL.Z = 0;
	  TEMP_MASS  = 0;

	  for(int a=0; a<SYSTEM.ATOMS; a++)
	  {
		 // Don't account for frozen atoms
			
		 if((CONTROLS.FREEZE_IDX_START != -1) && ((a<CONTROLS.FREEZE_IDX_START) || (a>CONTROLS.FREEZE_IDX_STOP)))
			continue;
			
		 TEMP_VEL.X += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].X;
		 TEMP_VEL.Y += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].Y;
		 TEMP_VEL.Z += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].Z;
		 TEMP_MASS  += SYSTEM.MASS[a];
	  }
		
	  // Check our velocity center of mass.. hopefully this is (or is very close to) zero
		
	  if(RANK==0)  
		 cout	<< "	Initial velocity center of mass: (" 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.X/TEMP_MASS << ", " 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Y/TEMP_MASS << ", " 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Z/TEMP_MASS << ") "
				<< endl;
	  
	  // In case it isn't, correct velocities to make it so: 
		 
	  for(int a=0; a<SYSTEM.ATOMS; a++)
	  {
		 // Don't account for frozen atoms
			
		 if((CONTROLS.FREEZE_IDX_START != -1) && ((a<CONTROLS.FREEZE_IDX_START) || (a>CONTROLS.FREEZE_IDX_STOP)))
			continue;
			
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
		 // Don't account for frozen atoms
			
		 if((CONTROLS.FREEZE_IDX_START != -1) && ((a<CONTROLS.FREEZE_IDX_START) || (a>CONTROLS.FREEZE_IDX_STOP)))
			continue;
			
		 TEMP_VEL.X += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].X;
		 TEMP_VEL.Y += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].Y;
		 TEMP_VEL.Z += SYSTEM.MASS[a] * SYSTEM.VELOCITY[a].Z;
		 TEMP_MASS  += SYSTEM.MASS[a];
	  }
		
	  if(RANK==0)
		 cout	<< "	Final   velocity center of mass: (" 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.X/TEMP_MASS << ", " 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Y/TEMP_MASS << ", " 
				<< setw(10) << fixed << setprecision(4) << right << TEMP_VEL.Z/TEMP_MASS << ") "
				<< endl;
}


void CONSTRAINT::WRITE(ofstream &output)
// Write parameters for restart.
{
	output << THERM_POSIT_T << endl;
	output << THERM_VELOC_T << endl;
	output << THERM_FORCE_T << endl;
	output << THERM_INERT_Q << endl;
	output << THERM_POSIT_0 << endl;
	output << THERM_FORCE_0 << endl;
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
	input >> THERM_FORCE_T;
	input >> THERM_INERT_Q;
	input >> THERM_POSIT_0;
	input >> THERM_FORCE_0;
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


void CONSTRAINT::UPDATE_COORDS(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, NEIGHBORS &NEIGHBORS)
{
	THERM_POSIT_0 = THERM_POSIT_T;
	THERM_VELOC_0 = THERM_VELOC_T;
	THERM_FORCE_0 = THERM_FORCE_T;
	
	BAROS_POSIT_0 = BAROS_POSIT_T;
	BAROS_VELOC_0 = BAROS_VELOC_T;
	BAROS_FORCE_0 = BAROS_FORCE_T;
	
	VOLUME_0 = VOLUME_T;

	///////////////////////////////////////////////
	// Do the un-constrained update of coords... this applies to all styles.
	///////////////////////////////////////////////

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
	
		SYSTEM.COORDS0[a1].X = SYSTEM.COORDS[a1].X ;
		SYSTEM.COORDS0[a1].Y = SYSTEM.COORDS[a1].Y ;
		SYSTEM.COORDS0[a1].Z = SYSTEM.COORDS[a1].Z ;		

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
		
		VOLUME_0 = SYSTEM.BOXDIM.UPDATE_VOLUME();
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
	
		SYSTEM.BOXDIM.CELL_AX *= BEREND_MU;
		SYSTEM.BOXDIM.CELL_BY *= BEREND_MU;
		SYSTEM.BOXDIM.CELL_CZ *= BEREND_MU;
		
		SYSTEM.BOXDIM.UPDATE_CELL();

		VOLUME_T = SYSTEM.BOXDIM.UPDATE_VOLUME();
	}
	
	///////////////////////////////////////////////
	// If requested, do Berendsen barostatting, anisotropically
	///////////////////////////////////////////////
	
	if(STYLE=="NPT-BEREND-ANISO")
	{
		VOLUME_0 = SYSTEM.BOXDIM.UPDATE_VOLUME();
		
		BEREND_ANI_MU.X = pow(1.0 - CONTROLS.DELTA_T/BEREND_KP*(CONTROLS.PRESSURE-SYSTEM.PRESSURE_TENSORS_ALL[0].X)/GPa,1.0/3.0);
		BEREND_ANI_MU.Y = pow(1.0 - CONTROLS.DELTA_T/BEREND_KP*(CONTROLS.PRESSURE-SYSTEM.PRESSURE_TENSORS_ALL[1].Y)/GPa,1.0/3.0);
		BEREND_ANI_MU.Z = pow(1.0 - CONTROLS.DELTA_T/BEREND_KP*(CONTROLS.PRESSURE-SYSTEM.PRESSURE_TENSORS_ALL[2].Z)/GPa,1.0/3.0);
		
	
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
	
		SYSTEM.BOXDIM.CELL_AX *= BEREND_ANI_MU.X;
		SYSTEM.BOXDIM.CELL_BY *= BEREND_ANI_MU.Y;
		SYSTEM.BOXDIM.CELL_CZ *= BEREND_ANI_MU.Z;

		SYSTEM.BOXDIM.UPDATE_CELL();

		VOLUME_T = SYSTEM.BOXDIM.UPDATE_VOLUME();

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

		THERM_POSIT_T = THERM_POSIT_0 + THERM_VELOC_0 * CONTROLS.DELTA_T + 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * THERM_FORCE_0;
	}

	///////////////////////////////////////////////
	// Do MTTK barostatting part - NPT-MTK
	///////////////////////////////////////////////
	
	if(STYLE=="NPT-MTK") 
	{
		// Update barostat position, compute scaling

		BAROS_POSIT_T = BAROS_POSIT_0 + BAROS_VELOC_0*CONTROLS.DELTA_T
			+ 0.5*CONTROLS.DELTA_T*CONTROLS.DELTA_T*BAROS_FORCE_0/BAROS_INERT_W
			- 0.5*CONTROLS.DELTA_T*CONTROLS.DELTA_T*BAROS_VELOC_0*THERM_VELOC_0;
		
		BAROS_SCALE   = exp(BAROS_POSIT_T - BAROS_POSIT_0);

		XYZ com = CENTER_OF_MASS(SYSTEM) ;

		for(int a1=0;a1<SYSTEM.ATOMS;a1++)	
		{	
			if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))
      // Don't account for frozen atoms
			{
				// Atom freezing doesn't work for a barostat!
				cout << "ERROR: Atoms cannot be frozen in an NPT ensemble." << endl;
				exit_run(0);
			}
			
			SYSTEM.COORDS[a1].X -= 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].X;
			SYSTEM.COORDS[a1].Y -= 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].Y;
			SYSTEM.COORDS[a1].Z -= 0.5 * CONTROLS.DELTA_T * CONTROLS.DELTA_T * (2.0 + 3.0/(N_DOF))*BAROS_VELOC_0 * SYSTEM.VELOCITY[a1].Z;	

			// AVOID MOTION OF THE CENTER OF MASS DUE TO BAROSTAT SCALING OF BOX SIZE.
			SYSTEM.COORDS[a1].X = BAROS_SCALE*(SYSTEM.COORDS[a1].X - com.X) + com.X ;	
			SYSTEM.COORDS[a1].Y = BAROS_SCALE*(SYSTEM.COORDS[a1].Y - com.Y) + com.Y ;	
			SYSTEM.COORDS[a1].Z = BAROS_SCALE*(SYSTEM.COORDS[a1].Z - com.Z) + com.Z ;	;				
		}
	
		SYSTEM.BOXDIM.CELL_AX *= BAROS_SCALE;
		SYSTEM.BOXDIM.CELL_BY *= BAROS_SCALE;
		SYSTEM.BOXDIM.CELL_CZ *= BAROS_SCALE;
	
		SYSTEM.BOXDIM.UPDATE_CELL();

		VOLUME_T = SYSTEM.BOXDIM.VOL ;
	}
	
	///////////////////////////////////////////////
	// Refresh ghost atom positions
	///////////////////////////////////////////////

	// Set the first NATOMS of ghost atoms to have the coordinates of the "real" coords

	SYSTEM.update_ghost(CONTROLS.N_LAYERS, false) ;


	NEIGHBORS.MAX_COORD_STEP = 0.0 ;
	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
		double step = sqrt( (SYSTEM.COORDS0[a1].X-SYSTEM.COORDS[a1].X) *  (SYSTEM.COORDS0[a1].X-SYSTEM.COORDS[a1].X) +
												(SYSTEM.COORDS0[a1].Y-SYSTEM.COORDS[a1].Y) *  (SYSTEM.COORDS0[a1].Y-SYSTEM.COORDS[a1].Y) +
												(SYSTEM.COORDS0[a1].Z-SYSTEM.COORDS[a1].Z) *  (SYSTEM.COORDS0[a1].Z-SYSTEM.COORDS[a1].Z) ) ;

		if ( step > NEIGHBORS.MAX_COORD_STEP ) NEIGHBORS.MAX_COORD_STEP = step ;

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

			SYSTEM.VELOCITY_ITER[a1].X *= BAROS_SCALE ;
			SYSTEM.VELOCITY_ITER[a1].Y *= BAROS_SCALE ;
			SYSTEM.VELOCITY_ITER[a1].Z *= BAROS_SCALE ;			
		}
	}	
}

void CONSTRAINT::UPDATE_VELOCS_HALF_2(FRAME & SYSTEM, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST)
{
	const double tol_iter = 1.0e-16 ;

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
			return ;
		}
	else if( STYLE=="NVT-MTK" )	
		{
			KIN_ENER = kinetic_energy(SYSTEM, CONTROLS);
			THERM_FORCE_T = ( 2.0 * KIN_ENER - N_DOF * Kb * CONTROLS.TEMPERATURE ) / (THERM_INERT_Q);
			THERM_VELOC_T = THERM_VELOC_0 + (THERM_FORCE_0 + THERM_FORCE_T) * 0.5 * CONTROLS.DELTA_T;
			
			double vscale_last = 0.0 ;
			double therm_veloc_last = THERM_VELOC_T ;
			double therm_force_last = THERM_FORCE_T ;
			int itr ;
			int itr_max = 10 ;
			
			// Iterative determination of velocities... See Martyna, Tobias, Klein JCP 101, 4177(1994) Appendix D.
			// For treatment of center of mass, see Melchiotta, Ciccotti, Holian, JCP 105, 346 (1996) and ref. therein.
			double err ;
			for ( itr = 0; itr < itr_max ; itr++ ) 
				{
					err = -1.0 ;
					THERM_VELOC_T = THERM_VELOC_0 + (THERM_FORCE_0 + THERM_FORCE_T) * 0.5 * CONTROLS.DELTA_T;

#if VERBOSITY >= 2 												
					if ( RANK == 0 )
						printf("THERM_VEL = %21.14e THERM_FORCE = %21.14e\n", THERM_VELOC_T, THERM_FORCE_T) ;
#endif												
					 
					VSCALEH = 1.0 + 0.5 * CONTROLS.DELTA_T * THERM_VELOC_T;

					for(int a1=0;a1<SYSTEM.ATOMS;a1++)
						{
							if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
								continue;
				
							SYSTEM.VELOCITY_NEW[a1].X = (SYSTEM.VELOCITY_ITER[a1].X + 0.5*SYSTEM.ACCEL[a1].X*CONTROLS.DELTA_T) / VSCALEH;
							SYSTEM.VELOCITY_NEW[a1].Y = (SYSTEM.VELOCITY_ITER[a1].Y + 0.5*SYSTEM.ACCEL[a1].Y*CONTROLS.DELTA_T) / VSCALEH;
							SYSTEM.VELOCITY_NEW[a1].Z = (SYSTEM.VELOCITY_ITER[a1].Z + 0.5*SYSTEM.ACCEL[a1].Z*CONTROLS.DELTA_T) / VSCALEH;
						}
					 
					KIN_ENER = kinetic_energy(SYSTEM,"NEW", CONTROLS);
					THERM_FORCE_T = ( 2.0 * KIN_ENER - N_DOF * Kb * CONTROLS.TEMPERATURE ) / (THERM_INERT_Q);

					if (fabs(therm_veloc_last - THERM_VELOC_T) > err ) err = fabs(therm_veloc_last - THERM_VELOC_T) ;
					if (fabs(therm_force_last - THERM_FORCE_T) > err ) err = fabs(therm_force_last - THERM_FORCE_T) ;
					if (fabs(vscale_last - VSCALEH) > err )            err = fabs(vscale_last - VSCALEH) ;
					
					//if ( err < tol_iter )
					//break ;

					vscale_last = VSCALEH ;
					therm_force_last = THERM_FORCE_T ;
					therm_veloc_last = THERM_VELOC_T ;
				}

			if ( err > tol_iter && RANK == 0 ) {
				 cout << "Warning: NVT-MTK velocity iteration did not converge.\n" ;
			}
			// if ( itr == itr_max && RANK == 0 ) {
			//    cout << "Warning: NVT-MTK velocity iteration did not converge.\n" ;
			// } else if ( RANK == 0 )
			// {
			// 	 cout << "NVT-MTK converged in " << itr + 1 << " iterations\n" ;
			// 	 printf("Error = %21.14e \n", fabs(VSCALEH - vscale_last) ) ;
			// }

		}
	else if ( STYLE == "NPT-MTK" ) // NPT-MTK
		{
			// Advance velocity without thermostat/barostat
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
				{
					if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
						continue;
				
					SYSTEM.VELOCITY_ITER[a1].X +=  0.5*SYSTEM.ACCEL[a1].X*CONTROLS.DELTA_T ;
					SYSTEM.VELOCITY_ITER[a1].Y +=  0.5*SYSTEM.ACCEL[a1].Y*CONTROLS.DELTA_T ;
					SYSTEM.VELOCITY_ITER[a1].Z +=  0.5*SYSTEM.ACCEL[a1].Z*CONTROLS.DELTA_T ;
				}

			//
			// Generate initial guess for iterative calculation of velocities.
			//
			
			// Update thermostat/barostat.
			KIN_ENER = kinetic_energy(SYSTEM,CONTROLS);

			// Use PRESSURE_XYZ + 2 KIN_ENER / (3V) so that kinetic energy contribution to pressure is updated.
			double dP = SYSTEM.PRESSURE_XYZ + 2.0 * KIN_ENER / (3.0 * VOLUME_T) - CONTROLS.PRESSURE/GPa;

			BAROS_FORCE_T = 3.0/(N_DOF) * 2.0 * KIN_ENER + 3.0 * VOLUME_T * dP;				
			THERM_FORCE_T = ( 2.0 * KIN_ENER + BAROS_INERT_W * BAROS_VELOC_0 * BAROS_VELOC_0
												- (N_DOF + 1.0) * Kb * CONTROLS.TEMPERATURE ) / (THERM_INERT_Q);
			THERM_VELOC_T = THERM_VELOC_0 + (THERM_FORCE_0 + THERM_FORCE_T) * 0.5 * CONTROLS.DELTA_T;

			// USE BAROS_VELOC_0 on RHS of equation for guess.  BAROS_VELOS_T is not yet evaluated.
			BAROS_VELOC_T = BAROS_VELOC_0 + 0.5*CONTROLS.DELTA_T*(BAROS_FORCE_0/BAROS_INERT_W
																														- BAROS_VELOC_0*THERM_VELOC_0
																														+ BAROS_FORCE_T/BAROS_INERT_W
																														- BAROS_VELOC_0*THERM_VELOC_T
																														) ;
			VSCALEH = 1.0 + 0.5 * CONTROLS.DELTA_T * (THERM_VELOC_T + (2.0 + 3.0/(N_DOF))*BAROS_VELOC_T);
			
			// Update velocity with thermostat/barostat
			for(int a1=0;a1<SYSTEM.ATOMS;a1++)
				{
					if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
						continue;
				
					SYSTEM.VELOCITY_NEW[a1].X = SYSTEM.VELOCITY_ITER[a1].X / VSCALEH ;
					SYSTEM.VELOCITY_NEW[a1].Y = SYSTEM.VELOCITY_ITER[a1].Y / VSCALEH ;
					SYSTEM.VELOCITY_NEW[a1].Z = SYSTEM.VELOCITY_ITER[a1].Z / VSCALEH ;
				}

			double baros_veloc_last = BAROS_VELOC_T ;
			double therm_veloc_last = THERM_VELOC_T ;
			double therm_force_last = THERM_FORCE_T ;
			double baros_force_last = BAROS_FORCE_T ;						

			//  Iterative determination of velocities.
			//  I made a conservative choice to always perform the same number of iterations, but
			//  the iteration loop could be terminated after convergence is achieved. (LEF 6/24/21)
			//
			int itr_max = 10 ;
			int itr ;
			double err = -1.0 ;
			
			for ( itr = 0 ; itr < itr_max ; itr++ )
				{



					// Update atomic velocity using new thermostat/barostat/velocity
					 // 
					err = -1.0 ;
					for(int a1=0;a1<SYSTEM.ATOMS;a1++)
						{
							if((CONTROLS.FREEZE_IDX_START != -1) && ((a1<CONTROLS.FREEZE_IDX_START) || (a1>CONTROLS.FREEZE_IDX_STOP)))	// Don't account for frozen atoms
								continue;

							XYZ V_OLD ;
							V_OLD.X = SYSTEM.VELOCITY_NEW[a1].X ;
							V_OLD.Y = SYSTEM.VELOCITY_NEW[a1].Y ;
							V_OLD.Z = SYSTEM.VELOCITY_NEW[a1].Z ;														

#if(0)
							// COMMENTED OUT !! THIS ITERATION METHOD WAS LESS STABLE THAN SCALING BY VSCALEH.
							// FORMALLY, I THINK THEY'RE THE SAME.  (LEF 6/24/21).
							SYSTEM.VELOCITY_NEW[a1].X = SYSTEM.VELOCITY_ITER[a1].X
								 - 0.5 * SYSTEM.VELOCITY_NEW[a1].X * ( (2.0 + 3.0/(N_DOF))*BAROS_VELOC_T 
																											 + THERM_VELOC_T) * CONTROLS.DELTA_T ;

							SYSTEM.VELOCITY_NEW[a1].Y = SYSTEM.VELOCITY_ITER[a1].Y
								 - 0.5 *  SYSTEM.VELOCITY_NEW[a1].Y * ( (2.0 + 3.0/(N_DOF))*BAROS_VELOC_T
																												+ THERM_VELOC_T) * CONTROLS.DELTA_T ;

							SYSTEM.VELOCITY_NEW[a1].Z = SYSTEM.VELOCITY_ITER[a1].Z
								 - 0.5 *  SYSTEM.VELOCITY_NEW[a1].Z * ( (2.0 + 3.0/(N_DOF))*BAROS_VELOC_T *
																												+ THERM_VELOC_T) * CONTROLS.DELTA_T ;
#else
							// SCALE VELOCITIES BY VSCALEH factor, as in NVT.
							SYSTEM.VELOCITY_NEW[a1].X = SYSTEM.VELOCITY_ITER[a1].X / VSCALEH;
							SYSTEM.VELOCITY_NEW[a1].Y = SYSTEM.VELOCITY_ITER[a1].Y / VSCALEH;
							SYSTEM.VELOCITY_NEW[a1].Z = SYSTEM.VELOCITY_ITER[a1].Z / VSCALEH;
#endif
							// Keep track of self-consistency error.
							if ( fabs(V_OLD.X - SYSTEM.VELOCITY_NEW[a1].X) > err ) err =  fabs(V_OLD.X - SYSTEM.VELOCITY_NEW[a1].X) ;
							if ( fabs(V_OLD.Y - SYSTEM.VELOCITY_NEW[a1].Y) > err ) err =  fabs(V_OLD.Y - SYSTEM.VELOCITY_NEW[a1].Y) ;
							if ( fabs(V_OLD.Z - SYSTEM.VELOCITY_NEW[a1].Z) > err ) err =  fabs(V_OLD.Z - SYSTEM.VELOCITY_NEW[a1].Z) ;							
						}

					KIN_ENER = kinetic_energy(SYSTEM,"NEW", CONTROLS);

					// Use PRESSURE_XYZ + 2 * KIN_ENER / (3 V) so that kinetic energy contribution to pressure is updated.
					dP = SYSTEM.PRESSURE_XYZ + 2.0 * KIN_ENER / (3.0 * VOLUME_T) - CONTROLS.PRESSURE/GPa;

					BAROS_FORCE_T = 3.0/(N_DOF) * 2.0 * KIN_ENER + 3.0 * VOLUME_T * dP;
					THERM_FORCE_T = ( 2.0 * KIN_ENER + BAROS_INERT_W * BAROS_VELOC_T * BAROS_VELOC_T
														- (N_DOF + 1.0) * Kb * CONTROLS.TEMPERATURE ) / (THERM_INERT_Q);
					THERM_VELOC_T = THERM_VELOC_0 + (THERM_FORCE_0 + THERM_FORCE_T) * 0.5 * CONTROLS.DELTA_T;
					BAROS_VELOC_T = BAROS_VELOC_0 + 0.5*CONTROLS.DELTA_T*(BAROS_FORCE_0/BAROS_INERT_W
																																- BAROS_VELOC_0*THERM_VELOC_0
																																+ BAROS_FORCE_T/BAROS_INERT_W
																																- BAROS_VELOC_T*THERM_VELOC_T
																																) ;
					VSCALEH = 1.0 + 0.5 * CONTROLS.DELTA_T * (THERM_VELOC_T + (2.0 + 3.0/(N_DOF))*BAROS_VELOC_T);					

					if( fabs(baros_veloc_last - BAROS_VELOC_T) > err ) err = fabs(baros_veloc_last - BAROS_VELOC_T) ;
					if (fabs(therm_veloc_last - THERM_VELOC_T) > err ) err =  fabs(therm_veloc_last - THERM_VELOC_T) ;
					if (fabs(therm_force_last - THERM_FORCE_T) > err ) err = fabs(therm_force_last - THERM_FORCE_T) ;
					if (fabs(baros_force_last - BAROS_FORCE_T) > err ) err = fabs(baros_force_last - BAROS_FORCE_T) ;
					
#if VERBOSITY >= 2 												
					if ( RANK == 0 )
					{
						 printf("BAROS_VEL   = %21.14e THERM_VEL   = %21.14e\n", BAROS_VELOC_T, THERM_VELOC_T) ;
						 printf("BAROS_FORCE = %21.14e THERM_FORCE = %21.14e\n", BAROS_FORCE_T, THERM_FORCE_T) ;
					}
#endif												
					// if ( err < tol_iter ) 
					// 	{
					// 		//#if VERBOSITY >= 2						
					// 		 if ( RANK == 0 )
					// 		 {
					// 				cout << "NPT-MTK converged in " << itr + 1 << " iterations\n" ;
					// 				printf("Error = %21.14e \n", err) ;
					// 		 }
							
					// 		//#endif
					//
					// 		 // break ;
					// 	}

					baros_veloc_last = BAROS_VELOC_T ;
					therm_veloc_last = THERM_VELOC_T ;
					therm_force_last = THERM_FORCE_T ;
					baros_force_last = BAROS_FORCE_T ;										

				}

			if ( err > tol_iter && RANK == 0 ) {
				 cout << "Warning: NPT-MTK iteration did not converge\n" ;
			}
		}
	else
		{
			EXIT_MSG("Error: an unknown constraint style") ;
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
	 double THERM_KE, THERM_PE, BAROS_KE, BAROS_PE;
	 double TEMP_VOL;
	
	 if(STYLE=="NVT-MTK")
	 {
			THERM_KE = 0.5 * THERM_VELOC_T * THERM_VELOC_T * THERM_INERT_Q;
			THERM_PE = N_DOF * Kb * CONTROLS.TEMPERATURE * THERM_POSIT_T;

			return  THERM_KE + THERM_PE;
	 }
	 else if ( STYLE == "NPT-MTK" ) // NPT-MTK
	 {
			TEMP_VOL = SYSTEM.BOXDIM.VOL;
	
			THERM_KE = 0.5 * THERM_VELOC_T * THERM_VELOC_T * THERM_INERT_Q;
			BAROS_KE = 0.5 * BAROS_VELOC_T * BAROS_VELOC_T * BAROS_INERT_W;

			THERM_PE = (N_DOF + 1) * Kb * CONTROLS.TEMPERATURE * THERM_POSIT_T;
			BAROS_PE = CONTROLS.PRESSURE/GPa*TEMP_VOL;

			return THERM_KE + THERM_PE + BAROS_KE + BAROS_PE;
	 } else {
			EXIT_MSG("ERROR: a conserved quantity was requested for an unrecognized constraint style") ;
	 }

	 // Fix compiler warning.  Code not reached.
	 return 0.0 ;
	
}

void NEIGHBORS::UPDATE_3B_INTERACTION(FRAME & SYSTEM, JOB_CONTROL &CONTROLS) 
// Build a list of all 3-body interactions.  This "flat" list parallelizes much
// more efficiently than a nested neighbor list loop.
{
	XYZ RAB;  
	INTERACTION_3B inter;

	LIST_3B_INT.clear();

	for ( int i = 0; i < SYSTEM.ATOMS; i++ ) 
	{
		int ai = i;
		for ( int j = 0; j < LIST_3B[i].size(); j++ ) 
		{
			int aj = LIST_3B[i][j];
			for ( int k = 0; k < LIST_3B[i].size(); k++ ) 
			{
				int ak = LIST_3B[i][k];

				if ( aj == ak )
				{
					continue;
				}
				else if ( PERM_SCALE[3] == 1.0 && SYSTEM.PARENT[aj] > SYSTEM.PARENT[ak] )
				{
					 continue ;
				}
				
				// The j-k list is possibly outside of the cutoff, so test it here.
				double rlen = get_dist(SYSTEM, RAB, aj, ak);
	
				if ( rlen < MAX_CUTOFF_3B + RCUT_PADDING ) 
				{
					inter.a1 = ai;
					inter.a2 = aj;
					inter.a3 = ak;
	  
					LIST_3B_INT.push_back(inter);
				}
			}
		}
	}
#if VERBOSITY >= 1 
	if ( RANK == 0 ) 
	  cout << "Number of 3-body interactions = " << LIST_3B_INT.size() << endl ;
#endif
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
		
		for (int j=0; j<LIST_4B[i].size(); j++) // Loop over all neighbors of i to get atom j
		{
			aj = LIST_4B[i][j];
			
			for (int k=0; k<LIST_4B[i].size(); k++) // Loop over all neighbors of i to get atom k
			{

			  ak = LIST_4B[i][k];
				
			  if (aj == ak )
				 continue;

			  if( PERM_SCALE[4] == 1.0 && SYSTEM.PARENT[aj] > SYSTEM.PARENT[ak] )
				 continue;

			  if( get_dist(SYSTEM, RAB, aj, ak) >  MAX_CUTOFF_4B + RCUT_PADDING)
				 continue;

				for (int l=0; l<LIST_4B[i].size(); l++) // Loop over all neighbors of i to get atom l
				{
					al = LIST_4B[i][l];
				
					// Check that this is a valid quadruplet

					if (aj == al || ak == al)
						continue;
					
					if( PERM_SCALE[4] == 1.0 && (SYSTEM.PARENT[ak] > SYSTEM.PARENT[al] || SYSTEM.PARENT[aj] > SYSTEM.PARENT[al]) )
						continue;
					
					// We know ij, ik, il, and jk distances are within the allowed cutoffs, but we still need to check jl, and kl

					if( get_dist(SYSTEM, RAB, aj, al) >=  MAX_CUTOFF_4B + RCUT_PADDING)
						continue;
					if( get_dist(SYSTEM, RAB, ak, al) >=  MAX_CUTOFF_4B + RCUT_PADDING)
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
#if VERBOSITY >= 1 
	if ( RANK == 0 ) 
	  cout << "Number of 4-body interactions = " << LIST_4B_INT.size() << endl ;
#endif
}

void THERMO_AVG::WRITE(ofstream &fout)
// Write out thermodynamic average properties.
{
	fout << TEMP_SUM << endl;
	fout << PRESS_SUM << endl;
	fout << STRESS_TENSOR_SUM.X << endl;
	fout << STRESS_TENSOR_SUM.Y << endl;
	fout << STRESS_TENSOR_SUM.Z << endl;
	fout << VOLUME_SUM << endl ;
	fout << PV_SUM << endl ;
}


void THERMO_AVG::READ(ifstream &fin)
// Read in thermodynamic average properties.
{
	fin >> TEMP_SUM ;
	fin >> PRESS_SUM ;
	fin >> STRESS_TENSOR_SUM.X ;
	fin >> STRESS_TENSOR_SUM.Y ;
	fin >> STRESS_TENSOR_SUM.Z ;
	fin >> VOLUME_SUM ;
	fin >> PV_SUM ;
}



BOX::BOX()
{
	IS_ORTHO = true;
	IS_VARIABLE = false ;
	
	CELL_AX = 0.0;
	CELL_AY = 0.0;
	CELL_AZ = 0.0;

	CELL_BX = 0.0;
	CELL_BY = 0.0;
	CELL_BZ = 0.0;

	CELL_CX = 0.0;
	CELL_CY = 0.0;
	CELL_CZ = 0.0;
	
	LAT_ALPHA = pi/2.0;
	LAT_BETA  = pi/2.0; 
	LAT_GAMMA = pi/2.0;
	
	HMAT.resize(9);
	INVR_HMAT.resize(9);   	
	
	for (int i=0; i<9; i++)
	{	
		HMAT[i]      = 0;
		INVR_HMAT[i] = 0;
	}
}



//void BOX::COPY_BOX_TO(BOX & TO_BOX)
BOX::BOX(const BOX & COPY_FROM)
{
	IS_ORTHO	 = COPY_FROM.IS_ORTHO;
	VOL		 = COPY_FROM.VOL;

	CELL_AX 	 = COPY_FROM.CELL_AX;   
	CELL_AY 	 = COPY_FROM.CELL_AY;   
	CELL_AZ 	 = COPY_FROM.CELL_AZ;   
	CELL_BX 	 = COPY_FROM.CELL_BX;   
	CELL_BY 	 = COPY_FROM.CELL_BY;   
	CELL_BZ 	 = COPY_FROM.CELL_BZ;   
	CELL_CX 	 = COPY_FROM.CELL_CX;   
	CELL_CY 	 = COPY_FROM.CELL_CY;   
	CELL_CZ 	 = COPY_FROM.CELL_CZ;
	
	HMAT.resize(9);
	INVR_HMAT.resize(9);   

	for (int i=0; i<9; i++)
	{	
		HMAT[i]      = COPY_FROM.HMAT[i];;
		INVR_HMAT[i] = COPY_FROM.INVR_HMAT[i];
	}

	LATCON_A	 = COPY_FROM.LATCON_A;  
	LATCON_B	 = COPY_FROM.LATCON_B;  
	LATCON_C	 = COPY_FROM.LATCON_C;  
	LAT_ALPHA	 = COPY_FROM.LAT_ALPHA;
	LAT_BETA	 = COPY_FROM.LAT_BETA;  
	LAT_GAMMA	 = COPY_FROM.LAT_GAMMA;

	CELL_LX 	 = COPY_FROM.CELL_LX;   
	CELL_LY 	 = COPY_FROM.CELL_LY;   
	CELL_LZ 	 = COPY_FROM.CELL_LZ;   

	XY		 = COPY_FROM.XY;	     
	XZ		 = COPY_FROM.XZ;	     
	YZ		 = COPY_FROM.YZ;	     

	EXTENT_X	 = COPY_FROM.EXTENT_X;
	EXTENT_Y	 = COPY_FROM.EXTENT_Y; 
	EXTENT_Z	 = COPY_FROM.EXTENT_Z;  
}

BOX::~BOX(){}

void BOX::WRITE_BOX(int LAYERS)
{
	cout << "       Box information (base units: Angstroms and radians): " << endl;
	cout << "       Orthorhombic:                   " << IS_ORTHO << endl;
	cout << "       Cell vectors (a)                " << CELL_AX	 << " " << CELL_AY	<< " " << CELL_AZ      << endl;
	cout << "       Cell vectors (b)                " << CELL_BX	 << " " << CELL_BY	     << " " << CELL_BZ      << endl;
	cout << "       Cell vectors (c)                " << CELL_CX	 << " " << CELL_CY	     << " " << CELL_CZ      << endl;
	cout << "       Layers:                         " << LAYERS << endl;
	cout << "       Ghost cell vectors (a) 	        " << CELL_AX * (2*LAYERS +1) << " " << CELL_AY * (2*LAYERS +1) << " " << CELL_AZ * (2*LAYERS +1) << endl;
	cout << "       Ghost cell vectors (a) 	        " << CELL_BX * (2*LAYERS +1) << " " << CELL_BY * (2*LAYERS +1) << " " << CELL_BZ * (2*LAYERS +1) << endl;
	cout << "       Ghost cell vectors (a) 	        " << CELL_CX * (2*LAYERS +1) << " " << CELL_CY * (2*LAYERS +1) << " " << CELL_CZ * (2*LAYERS +1) << endl;
	cout << "       Invr. cell vectors (a) 	        " << INVR_HMAT[0] << " " << INVR_HMAT[3] << " " << INVR_HMAT[6] << endl;
	cout << "       Invr. cell vectors (b) 	        " << INVR_HMAT[1] << " " << INVR_HMAT[4] << " " << INVR_HMAT[7] << endl;
	cout << "       Invr. cell vectors (c) 	        " << INVR_HMAT[2] << " " << INVR_HMAT[5] << " " << INVR_HMAT[8] << endl;
	cout << "       Lat. constants (a,b,c) 	        " << LATCON_A	 << " " << LATCON_B	<< " " << LATCON_C     << endl; 	   
	cout << "       Lat. angles (alpha,beta,gamma)  " << LAT_ALPHA    << " " << LAT_BETA     << " " << LAT_GAMMA    << endl;
	cout << "       Cell extents (x,y,z)            " << EXTENT_X     << " " << EXTENT_Y     << " " << EXTENT_Z	  << endl;
	cout << "       LMP cell lengths (lx,ly,lz)     " << CELL_LX      << " " << CELL_LY	   << " " << CELL_LZ	  << endl;
	cout << "       LMP cell tilts facs. (xy,xz,yz) " << XY	    << " " << XZ	   << " " << YZ 	  << endl;
	cout << "       Cell volume                     " << VOL << endl;
	cout << endl;
}



void BOX::UPDATE_INVER_CELL()
{
	// Compute determinant of the HMAT... if 0, no inverse exists
	
	
	HMAT[0] = CELL_AX;
	HMAT[1] = CELL_BX;
	HMAT[2] = CELL_CX;
	
	HMAT[3] = CELL_AY;
	HMAT[4] = CELL_BY;
	HMAT[5] = CELL_CY;
	
	HMAT[6] = CELL_AZ;
	HMAT[7] = CELL_BZ;
	HMAT[8] = CELL_CZ;	
	
	double HMAT_det = HMAT[0] * (HMAT[4]*HMAT[8] - HMAT[5]*HMAT[7])
	                - HMAT[1] * (HMAT[3]*HMAT[8] - HMAT[5]*HMAT[6])
			+ HMAT[2] * (HMAT[3]*HMAT[7] - HMAT[4]*HMAT[6]);
	
	if (HMAT_det == 0)
		EXIT_MSG("ERROR: H-matrix determinant of zero computed in INVERT_HMAT");

	// Compute the adjugate matrix and transpose
	
	vector<double> TMP(9);
	
	TMP[0] =      (HMAT[4]*HMAT[8] - HMAT[5]*HMAT[7]);
	TMP[1] = -1 * (HMAT[3]*HMAT[8] - HMAT[5]*HMAT[6]);
	TMP[2] =      (HMAT[3]*HMAT[7] - HMAT[4]*HMAT[6]);
	
	TMP[3] = -1 * (HMAT[1]*HMAT[8] - HMAT[2]*HMAT[7]);
	TMP[4] =      (HMAT[0]*HMAT[8] - HMAT[2]*HMAT[6]);
	TMP[5] = -1 * (HMAT[0]*HMAT[7] - HMAT[1]*HMAT[6]);
	
	TMP[6] =      (HMAT[1]*HMAT[5] - HMAT[2]*HMAT[4]);
	TMP[7] = -1 * (HMAT[0]*HMAT[5] - HMAT[2]*HMAT[3]);
	TMP[8] =      (HMAT[0]*HMAT[4] - HMAT[1]*HMAT[3]);
	
	INVR_HMAT[0] = TMP[0];
	INVR_HMAT[1] = TMP[3];
	INVR_HMAT[2] = TMP[6];

	INVR_HMAT[3] = TMP[1];
	INVR_HMAT[4] = TMP[4];
	INVR_HMAT[5] = TMP[7];

	INVR_HMAT[6] = TMP[2];
	INVR_HMAT[7] = TMP[5];
	INVR_HMAT[8] = TMP[8];
	
	// Save to the inverse H-Matrix
	
	INVR_HMAT[0] /= HMAT_det;
	INVR_HMAT[1] /= HMAT_det;
	INVR_HMAT[2] /= HMAT_det;
	
	INVR_HMAT[3] /= HMAT_det;
	INVR_HMAT[4] /= HMAT_det;
	INVR_HMAT[5] /= HMAT_det;
	
	INVR_HMAT[6] /= HMAT_det;
	INVR_HMAT[7] /= HMAT_det;
	INVR_HMAT[8] /= HMAT_det;	
}

void BOX::UPDATE_LAT_VALUES()
{
	vector<double> LAT_VEC, LAT_VEC_1, LAT_VEC_2;
	
	LAT_VEC = {CELL_AX, CELL_AY, CELL_AZ};
	LATCON_A = VECTOR_MAGNITUDE(LAT_VEC);
	
	LAT_VEC = {CELL_BX, CELL_BY, CELL_BZ};
	LATCON_B = VECTOR_MAGNITUDE(LAT_VEC);
	
	LAT_VEC = {CELL_CX, CELL_CY, CELL_CZ};
	LATCON_C = VECTOR_MAGNITUDE(LAT_VEC);
	
	
	LAT_VEC_1 = {CELL_BX, CELL_BY, CELL_BZ};
	LAT_VEC_2 = {CELL_CX, CELL_CY, CELL_CZ};
	LAT_ALPHA = VECTOR_ANGLE(LAT_VEC_1, LAT_VEC_2);
	
	LAT_VEC_1 = {CELL_CX, CELL_CY, CELL_CZ};
	LAT_VEC_2 = {CELL_AX, CELL_AY, CELL_AZ}; 
	LAT_BETA  = VECTOR_ANGLE(LAT_VEC_1, LAT_VEC_2);
	
	LAT_VEC_1 = {CELL_AX, CELL_AY, CELL_AZ};
	LAT_VEC_2 = {CELL_BX, CELL_BY, CELL_BZ}; 
	LAT_GAMMA = VECTOR_ANGLE(LAT_VEC_1, LAT_VEC_2);
}

void BOX::UPDATE_EXTENT()
{
	// See LAMMPS manual for formula
	
	
	double tmp;
	    
	double xlo = 0.0;
	if (XY < xlo)
		xlo = XY;
	if (XZ < xlo)
		xlo = XZ;
	if(XY+XZ < xlo)
		xlo = XY+XZ;
		
	double ylo = 0.0;
	if (YZ< ylo)
		ylo = YZ;
	double zlo = 0.0;
	
	double xhi = CELL_AX;
	tmp = 0.0;
	if (XY > tmp)
		tmp = XY;
	if (XZ > tmp)
		tmp = XZ;
	if (XY+XZ> tmp)
		tmp = XY+XZ;
	xhi += tmp;
	
	double yhi = CELL_BY;
	if(YZ > 0.0)
		yhi += YZ;
		
	double zhi = CELL_CZ;
	
	
	
	EXTENT_X = xhi - xlo;
	EXTENT_Y = yhi - ylo;
	EXTENT_Z = zhi - zlo;			
}

void BOX::UPDATE_CELL() // Assumes CELL_* values have been set, either orthorhombically or triclinically
{
	// Build the h-mat
	
	HMAT[0] = CELL_AX;
	HMAT[1] = CELL_BX;
	HMAT[2] = CELL_CX;
	
	HMAT[3] = CELL_AY;
	HMAT[4] = CELL_BY;
	HMAT[5] = CELL_CY;
	
	HMAT[6] = CELL_AZ;
	HMAT[7] = CELL_BZ;
	HMAT[8] = CELL_CZ;

	UPDATE_INVER_CELL();	
	UPDATE_LAT_VALUES();
	LAMMPSIFY();
	//UNLAMMPSIFY(CELL_LX, CELL_LY, CELL_LZ, XY, XZ, YZ);
	//UPDATE_INVER_CELL();
	//UPDATE_LAT_VALUES();	
	UPDATE_EXTENT();
	UPDATE_VOLUME();
}

void BOX::LAMMPSIFY()
{
	/* Here's the ~LAMMPsian~ meaning of all these geometric variables 
	
	See: https://lammps.sandia.gov/doc/Howto_triclinic.html

	(y)
	^
	|     EXTENT_X (lammps xlo_bound)
	||---------------|
	|      LATCON_A
	|    ----------->
	|   . . . . . . .        _ _ yhi
	|  . TRICLINIC .          |
	| . * SYSTEM  .           |  CELL_LY (lammps ly)
	|. . . . . . .______(x)  _|_ ylo
	|   *        |
	xlo *        xhi
	*   *
	*---*
	 ^xy
	
	
	XY, XZ, and YZ are tilt factors (distance, not an angle) in each plane and equal zero for an orthorhombic system
	
	*/

	
	CELL_LX = LATCON_A;
	
	XY      = LATCON_B * cos (LAT_GAMMA);
	
	if (abs(XY) < 1E-12)
		XY = 0.0;
	
	XZ      = LATCON_C * cos (LAT_BETA );
	
	if (abs(XZ) < 1E-12)
		XZ = 0.0;
	
	CELL_LY = sqrt(LATCON_B*LATCON_B - XY*XY);
	
	YZ      = (LATCON_B * LATCON_C * cos(LAT_ALPHA) - XY*XZ)/CELL_LY;
	
	if (abs(YZ) < 1E-12)
		YZ = 0.0;	
	
	CELL_LZ = sqrt( LATCON_C*LATCON_C - XZ*XZ -YZ*YZ );
}

void BOX::UNLAMMPSIFY(double lx, double ly, double lz, double xy, double xz, double yz)
{
	CELL_LX = lx;
	CELL_LY = ly;
	CELL_LZ = lz;
	
	LATCON_A = lx;
	LATCON_B = sqrt(ly*ly + xy*xy);
	LATCON_C = sqrt(lz*lz + xz*xz +yz*yz);
	
	LAT_ALPHA = acos( (xy*xz + ly*yz) / (LATCON_B*LATCON_C) );
	LAT_BETA  = acos(xz/LATCON_C);
	LAT_GAMMA = acos(xy/LATCON_B);
	
	CELL_AX = LATCON_A;
	CELL_AY = 0;
	CELL_AZ = 0;
		
	CELL_BX = LATCON_B*cos(LAT_GAMMA);
	CELL_BY = LATCON_B*sin(LAT_GAMMA);
	CELL_BZ = 0;	
	
	CELL_CX = LATCON_C*cos(LAT_BETA);
	CELL_CY = (LATCON_B*LATCON_C*cos(LAT_ALPHA) - CELL_BX*CELL_CX)/CELL_BY;
	CELL_CZ = sqrt(LATCON_C*LATCON_C + CELL_CX*CELL_CX + CELL_CY*CELL_CY);	
	
	UPDATE_EXTENT();	
	
}

void BOX::WRAP_ATOM(XYZ & UNWRAPPED_ATOM, XYZ_INT & WRAP_IDX, bool UPDATE_WRAPDIM) // Saves wrapped coordinates to UN_WRAPPED_ATOM
{
	XYZ TMP_ATOM;
	
	TMP_ATOM.X = UNWRAPPED_ATOM.X;
	TMP_ATOM.Y = UNWRAPPED_ATOM.Y;
	TMP_ATOM.Z = UNWRAPPED_ATOM.Z;

	WRAP_ATOM(TMP_ATOM, UNWRAPPED_ATOM, WRAP_IDX, UPDATE_WRAPDIM);
}

void BOX::WRAP_ATOM(XYZ & UNWRAPPED_ATOM, XYZ & WRAPPED_ATOM,  XYZ_INT & WRAP_IDX, bool UPDATE_WRAPDIM)	// Saves wrapped coordinates to WRAPPED_ATOM
{

	if (IS_ORTHO) // Then we can do this the fast/cheap way
	{
	
		if (UPDATE_WRAPDIM)
		{
			WRAP_IDX.X = floor(UNWRAPPED_ATOM.X/CELL_LX);
			WRAP_IDX.Y = floor(UNWRAPPED_ATOM.Y/CELL_LY);
			WRAP_IDX.Z = floor(UNWRAPPED_ATOM.Z/CELL_LZ);
		}
		
		WRAPPED_ATOM.X = UNWRAPPED_ATOM.X - WRAP_IDX.X * CELL_LX;
		WRAPPED_ATOM.Y = UNWRAPPED_ATOM.Y - WRAP_IDX.Y * CELL_LY;
		WRAPPED_ATOM.Z = UNWRAPPED_ATOM.Z - WRAP_IDX.Z * CELL_LZ;	   
	}
	else	// We do this the hard way: transform the "system" (single atom in the box) to an orthorhombic system, wrap, then undo the transformation
	{
		// Transform the "system" (single atom in the box) to an orthorhombic system
				
		XYZ TMP_ATOM;
		
		TMP_ATOM.X = INVR_HMAT[0]*UNWRAPPED_ATOM.X + INVR_HMAT[1]*UNWRAPPED_ATOM.Y + INVR_HMAT[2]*UNWRAPPED_ATOM.Z;
		TMP_ATOM.Y = INVR_HMAT[3]*UNWRAPPED_ATOM.X + INVR_HMAT[4]*UNWRAPPED_ATOM.Y + INVR_HMAT[5]*UNWRAPPED_ATOM.Z;
		TMP_ATOM.Z = INVR_HMAT[6]*UNWRAPPED_ATOM.X + INVR_HMAT[7]*UNWRAPPED_ATOM.Y + INVR_HMAT[8]*UNWRAPPED_ATOM.Z;
		
		// Wrap the atom in the transformed box
		
		if (UPDATE_WRAPDIM)
		{
			WRAP_IDX.X = floor(TMP_ATOM.X); // floor(UNWRAPPED_ATOM.X/CELL_LX);
			WRAP_IDX.Y = floor(TMP_ATOM.Y); // floor(UNWRAPPED_ATOM.Y/CELL_LY);
			WRAP_IDX.Z = floor(TMP_ATOM.Z); // floor(UNWRAPPED_ATOM.Z/CELL_LZ);
		}		
		
		TMP_ATOM.X -= WRAP_IDX.X; // round(TMP_ATOM.X);
		TMP_ATOM.Y -= WRAP_IDX.Y; // round(TMP_ATOM.Y);
	    	TMP_ATOM.Z -= WRAP_IDX.Z; // round(TMP_ATOM.Z);
	
		// Undo the transformation
		
		WRAPPED_ATOM.X = HMAT[0]*TMP_ATOM.X + HMAT[1]*TMP_ATOM.Y + HMAT[2]*TMP_ATOM.Z;
		WRAPPED_ATOM.Y = HMAT[3]*TMP_ATOM.X + HMAT[4]*TMP_ATOM.Y + HMAT[5]*TMP_ATOM.Z;
		WRAPPED_ATOM.Z = HMAT[6]*TMP_ATOM.X + HMAT[7]*TMP_ATOM.Y + HMAT[8]*TMP_ATOM.Z;
	}
}

void BOX::LAYER_ATOM(XYZ & REFERENCE_ATOM, XYZ_INT & LAYER_INDEX)	// In place (saves to reference atom)
{
	XYZ TMP_ATOM;
	
	TMP_ATOM.X = REFERENCE_ATOM.X;
	TMP_ATOM.Y = REFERENCE_ATOM.Y;
	TMP_ATOM.Z = REFERENCE_ATOM.Z;

	LAYER_ATOM(TMP_ATOM, LAYER_INDEX, REFERENCE_ATOM);
}

void BOX::LAYER_ATOM(XYZ & REFERENCE_ATOM, XYZ_INT & LAYER_INDEX, XYZ & LAYERED_ATOM) // Saves to layered atom
{
	// Apply shift along cell vectors a, b, and c
	
	if (IS_ORTHO) // Then we can do this the fast/cheap way
	{
		LAYERED_ATOM.X = REFERENCE_ATOM.X + LAYER_INDEX.X*CELL_AX;
		LAYERED_ATOM.Y = REFERENCE_ATOM.Y + LAYER_INDEX.Y*CELL_BY;
		LAYERED_ATOM.Z = REFERENCE_ATOM.Z + LAYER_INDEX.Z*CELL_CZ;	    
	}	
	else
	{
		XYZ TMP_ATOM;

		// Transform the "system" (single atom in the box) to an orthorhombic system
		
		TMP_ATOM.X = INVR_HMAT[0]*REFERENCE_ATOM.X + INVR_HMAT[1]*REFERENCE_ATOM.Y + INVR_HMAT[2]*REFERENCE_ATOM.Z;
		TMP_ATOM.Y = INVR_HMAT[3]*REFERENCE_ATOM.X + INVR_HMAT[4]*REFERENCE_ATOM.Y + INVR_HMAT[5]*REFERENCE_ATOM.Z;
		TMP_ATOM.Z = INVR_HMAT[6]*REFERENCE_ATOM.X + INVR_HMAT[7]*REFERENCE_ATOM.Y + INVR_HMAT[8]*REFERENCE_ATOM.Z;	
		
		// Layer in the reduced coordinate system 	
		
		TMP_ATOM.X += LAYER_INDEX.X;
		TMP_ATOM.Y += LAYER_INDEX.Y;
		TMP_ATOM.Z += LAYER_INDEX.Z;
		
		// Undo the transformation
		
		LAYERED_ATOM.X = HMAT[0]*TMP_ATOM.X + HMAT[1]*TMP_ATOM.Y + HMAT[2]*TMP_ATOM.Z;
		LAYERED_ATOM.Y = HMAT[3]*TMP_ATOM.X + HMAT[4]*TMP_ATOM.Y + HMAT[5]*TMP_ATOM.Z;
		LAYERED_ATOM.Z = HMAT[6]*TMP_ATOM.X + HMAT[7]*TMP_ATOM.Y + HMAT[8]*TMP_ATOM.Z;  	
	}
}

bool BOX::IS_RCUT_SAFE(double CUTOFF, int LAYERS)
{
	vector<double> LAT_VEC;
	
	LAT_VEC = {CELL_AX, CELL_AY, CELL_AZ};
	if ( CUTOFF >= 0.5*VECTOR_MAGNITUDE(LAT_VEC)*(2*LAYERS+1))
	{
		cout << "AX: " << VECTOR_MAGNITUDE(LAT_VEC)*(2*LAYERS+1) << " cut: " << CUTOFF << endl;
		return false;
	}
	
	LAT_VEC = {CELL_BX, CELL_BY, CELL_BZ};
	if ( CUTOFF >= 0.5*VECTOR_MAGNITUDE(LAT_VEC)*(2*LAYERS+1))
	{
		cout << "BY: " << VECTOR_MAGNITUDE(LAT_VEC)*(2*LAYERS+1) << " cut: " << CUTOFF << endl;
		return false;
	}
	
	LAT_VEC = {CELL_CX, CELL_CY, CELL_CZ};
	if ( CUTOFF >= 0.5*VECTOR_MAGNITUDE(LAT_VEC)*(2*LAYERS+1))
	{
		cout << "CZ: " << VECTOR_MAGNITUDE(LAT_VEC)*(2*LAYERS+1) << " cut: " << CUTOFF << endl;
		return false;
	}
	
	return true;
}

double BOX::UPDATE_VOLUME()
{
	if (IS_ORTHO)
	{
		VOL = CELL_AX * CELL_BY * CELL_CZ;
	}
	else
	{
		VOL  = 1;
		VOL += 2*cos(LAT_ALPHA)*cos(LAT_BETA)*cos(LAT_GAMMA);
		VOL -= cos(LAT_ALPHA)*cos(LAT_ALPHA);
		VOL -= cos(LAT_BETA)*cos(LAT_BETA);
		VOL -= cos(LAT_GAMMA)*cos(LAT_GAMMA);
		
		VOL = LATCON_A*LATCON_B*LATCON_C * sqrt(VOL);
	}
	return VOL;
}

void BOX::SCALE_BY_FACTOR(double FACTOR)
{
	static XYZ DUMMY = {0.0, 0.0, 0.0};
	
	SCALE_BY_FACTOR(FACTOR, false, DUMMY);
}


void BOX::SCALE_BY_FACTOR(double FACTOR, bool SCALE_ATOMS, XYZ & ATOM)
{
	// if SCALE_ATOMS is true, ONLY updates COORDS and ALL_COORDS 
	
	if (SCALE_ATOMS)
	{
		XYZ TMP_ATOM;
		
		TMP_ATOM.X = ATOM.X;
		TMP_ATOM.Y = ATOM.Y;
		TMP_ATOM.Z = ATOM.Z;
	
		// Transform the "system" (single atom in the box) to an orthorhombic system
				
		ATOM.X = INVR_HMAT[0]*TMP_ATOM.X + INVR_HMAT[1]*TMP_ATOM.Y + INVR_HMAT[2]*TMP_ATOM.Z;
		ATOM.Y = INVR_HMAT[3]*TMP_ATOM.X + INVR_HMAT[4]*TMP_ATOM.Y + INVR_HMAT[5]*TMP_ATOM.Z;
		ATOM.Z = INVR_HMAT[6]*TMP_ATOM.X + INVR_HMAT[7]*TMP_ATOM.Y + INVR_HMAT[8]*TMP_ATOM.Z; 		    
		
		// transform back to the scaled system 
		
		TMP_ATOM.X = FACTOR*HMAT[0]*ATOM.X + FACTOR*HMAT[1]*ATOM.Y + FACTOR*HMAT[2]*ATOM.Z;
    		TMP_ATOM.Y = FACTOR*HMAT[3]*ATOM.X + FACTOR*HMAT[4]*ATOM.Y + FACTOR*HMAT[5]*ATOM.Z;
		TMP_ATOM.Z = FACTOR*HMAT[6]*ATOM.X + FACTOR*HMAT[7]*ATOM.Y + FACTOR*HMAT[8]*ATOM.Z;
		
		ATOM.X = TMP_ATOM.X;
		ATOM.Y = TMP_ATOM.Y;
		ATOM.Z = TMP_ATOM.Z;
	}
	else
	{
		CELL_AX *= FACTOR;
		CELL_AY *= FACTOR;
		CELL_AZ *= FACTOR;
		
		CELL_BX *= FACTOR;
		CELL_BY *= FACTOR;
		CELL_BZ *= FACTOR;
		
		CELL_CX *= FACTOR;
		CELL_CY *= FACTOR;
		CELL_CZ *= FACTOR;
	}
}

void BOX::GET_DISTANCE(const XYZ & ATOM1, const XYZ & ATOM2, XYZ & RAB, bool USE_MIC)
{

	// Calculates distance as a2 - a1... This function modifies RAB!

	// Convert distances to scaled units
	
	XYZ ATOM1_SCALED, ATOM2_SCALED, TMP_RAB;
	
	ATOM1_SCALED.X = INVR_HMAT[0]*ATOM1.X + INVR_HMAT[1]*ATOM1.Y + INVR_HMAT[2]*ATOM1.Z;
	ATOM1_SCALED.Y = INVR_HMAT[3]*ATOM1.X + INVR_HMAT[4]*ATOM1.Y + INVR_HMAT[5]*ATOM1.Z;
	ATOM1_SCALED.Z = INVR_HMAT[6]*ATOM1.X + INVR_HMAT[7]*ATOM1.Y + INVR_HMAT[8]*ATOM1.Z;
	
	ATOM2_SCALED.X = INVR_HMAT[0]*ATOM2.X + INVR_HMAT[1]*ATOM2.Y + INVR_HMAT[2]*ATOM2.Z;
	ATOM2_SCALED.Y = INVR_HMAT[3]*ATOM2.X + INVR_HMAT[4]*ATOM2.Y + INVR_HMAT[5]*ATOM2.Z;
	ATOM2_SCALED.Z = INVR_HMAT[6]*ATOM2.X + INVR_HMAT[7]*ATOM2.Y + INVR_HMAT[8]*ATOM2.Z;
	
	// Get the (scaled) distance and apply MIC if required
	
	TMP_RAB.X = ATOM2_SCALED.X - ATOM1_SCALED.X;
	TMP_RAB.Y = ATOM2_SCALED.Y - ATOM1_SCALED.Y;
	TMP_RAB.Z = ATOM2_SCALED.Z - ATOM1_SCALED.Z;


	if (USE_MIC)
	{
		TMP_RAB.X -= round( TMP_RAB.X );
		TMP_RAB.Y -= round( TMP_RAB.Y );
		TMP_RAB.Z -= round( TMP_RAB.Z );		
	}
	
	// Convert back to standard units 

	
	RAB.X = HMAT[0]*TMP_RAB.X + HMAT[1]*TMP_RAB.Y + HMAT[2]*TMP_RAB.Z;
	RAB.Y = HMAT[3]*TMP_RAB.X + HMAT[4]*TMP_RAB.Y + HMAT[5]*TMP_RAB.Z;
	RAB.Z = HMAT[6]*TMP_RAB.X + HMAT[7]*TMP_RAB.Y + HMAT[8]*TMP_RAB.Z;
}


void FRAME::update_ghost(int n_layers, bool UPDATE_WRAPDIM)
// Update the ghost atoms using the given number of layers.
{

	for (int a=0; a<ATOMS; a++) 
		BOXDIM.WRAP_ATOM(COORDS[a], ALL_COORDS[a],WRAP_IDX[a], UPDATE_WRAPDIM);		
	
	// Build the surrounding "cell's" ghost atoms based on the first NATOMS ghost atoms
	
	XYZ_INT TEMP_LAYER;
	
	if(n_layers>0 )
	{	
		int TEMP_IDX = ATOMS;	

		for(TEMP_LAYER.X = -n_layers; TEMP_LAYER.X<=n_layers; TEMP_LAYER.X++)
		{
			for(TEMP_LAYER.Y = -n_layers; TEMP_LAYER.Y<=n_layers; TEMP_LAYER.Y++)
			{
				for(TEMP_LAYER.Z = -n_layers; TEMP_LAYER.Z<=n_layers; TEMP_LAYER.Z++)
				{	
					if (TEMP_LAYER.X == 0 && TEMP_LAYER.Y == 0 && TEMP_LAYER.Z == 0 ) 
						continue;
					else
					{
						for(int a1=0; a1<ATOMS; a1++)
						{
							BOXDIM.LAYER_ATOM(ALL_COORDS[a1], TEMP_LAYER, ALL_COORDS[TEMP_IDX]);
				
							if(PARENT[TEMP_IDX] != a1)
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
		
		if ( TEMP_IDX != ALL_ATOMS ) 
		{
			printf("Error updating layers\n");
			exit(1);
		}
	}	
}

void FRAME::READ_XYZF(ifstream &TRAJ_INPUT, const JOB_CONTROL &CONTROLS, const vector<PAIRS> &ATOM_PAIRS, const vector<string> &TMP_ATOMTYPE, int i)
// Read values from the xyzf file into the FRAME.
{
	TRAJ_INPUT >> ATOMS;
		
	// Read in line with box dimenstions
	string header;
	vector<string> tokens;

	// Read line twice to get through newline from last input.
	
	std::getline(TRAJ_INPUT, header);
	std::getline(TRAJ_INPUT, header);

	int ntokens = parse_space(header, tokens);

	// Make sure we at least have boxlengths
	
	if (tokens[0] == "NON_ORTHO")
		BOXDIM.IS_ORTHO = false;

	if ( ntokens >= 3 ) 
	{
		if (BOXDIM.IS_ORTHO)
		{
			BOXDIM.CELL_AX = stod(tokens[0]);
			BOXDIM.CELL_BY = stod(tokens[1]);
			BOXDIM.CELL_CZ = stod(tokens[2]);
			
			BOXDIM.UPDATE_CELL();
		}
		else
		{
			BOXDIM.CELL_AX = stod(tokens[1]); BOXDIM.CELL_AY = stod(tokens[2]); BOXDIM.CELL_AZ = stod(tokens[3]); 
			BOXDIM.CELL_BX = stod(tokens[4]); BOXDIM.CELL_BY = stod(tokens[5]); BOXDIM.CELL_BZ = stod(tokens[6]);
			BOXDIM.CELL_CX = stod(tokens[7]); BOXDIM.CELL_CY = stod(tokens[8]); BOXDIM.CELL_CZ = stod(tokens[9]);
			
			BOXDIM.UPDATE_CELL();
		}
	} 
	else 
	{
		cout << "Error:  Reading frame " << i << endl; // << " of file " << CONTROLS.INFILE << endl;
		cout << "        Missing box length components." << endl;
		cout << "        See offending line below." << endl;
		cout << header << endl;
		exit(1);			
	}
		
	// If requested, check for stress tensors... keep in mind that the 
	// number of tensors and the number of frames with tensors can vary

	if((CONTROLS.NSTRESS < 0) || (i<CONTROLS.NSTRESS))
	{
		if(CONTROLS.FIT_STRESS)
		{
			if (BOXDIM.IS_ORTHO && (ntokens >= 6))
			{
				STRESS_TENSORS.X = stod(tokens[3]);
				STRESS_TENSORS.Y = stod(tokens[4]);
				STRESS_TENSORS.Z = stod(tokens[5]);
			} 
			else if ( (!BOXDIM.IS_ORTHO) && (ntokens >= 13))
			{
				STRESS_TENSORS.X = stod(tokens[10 ]);
				STRESS_TENSORS.Y = stod(tokens[11]);
				STRESS_TENSORS.Z = stod(tokens[12]);
			} 			
			else 
			{
				cout << "Error:  Reading frame " << i << endl; // << " of file " << CONTROLS.INFILE << endl;
				cout << "        Missing diagonal stress tensor components." << endl;
				cout << "        See offending line below." << endl;
				cout << header << endl;
				exit(1);
			}
		}
		else if(CONTROLS.FIT_STRESS_ALL)	// Expects as:xx yy zz xy xz yz // old: xx, xy, xz, yy, yx, yz, zx, zy, zz
		{
			// Read only the "upper" deviatoric components
		
			if (BOXDIM.IS_ORTHO && ( ntokens >= 9 ))
			{
				STRESS_TENSORS_X.X = stod(tokens[3]);
				STRESS_TENSORS_Y.Y = stod(tokens[4]);
				STRESS_TENSORS_Z.Z = stod(tokens[5]);
					
				STRESS_TENSORS_X.Y = stod(tokens[6]);
				STRESS_TENSORS_X.Z = stod(tokens[7]);
				STRESS_TENSORS_Y.Z = stod(tokens[8]);
			} 
			else if ( (!BOXDIM.IS_ORTHO) && (ntokens >= 16))
			{
				STRESS_TENSORS_X.X = stod(tokens[10 ]);
				STRESS_TENSORS_Y.Y = stod(tokens[11]);
				STRESS_TENSORS_Z.Z = stod(tokens[12]);
					
				STRESS_TENSORS_X.Y = stod(tokens[13]);
				STRESS_TENSORS_X.Z = stod(tokens[14]);
				STRESS_TENSORS_Y.Z = stod(tokens[15]);			
			}
			else 
			{
				cout << "Error:  Reading frame " << i  << endl; // << " of file " << CONTROLS.INFILE << endl;
				cout << "        Missing full stress tensor components." << endl;
				cout << "        See offending line below." << endl;
				cout << header << endl;
				exit(1);
			}
		}
	}
	if((CONTROLS.NENER < 0) || (i<CONTROLS.NENER))
		if(CONTROLS.FIT_ENER) // We're fitting to the absolute energy, + an offset (column of 1's at end of A-matrix)
			QM_POT_ENER = stod(tokens[tokens.size()-1]);

	
	// Check that outer cutoffs do not exceed half of the boxlength
	// with consideration of layering
		
	for(int j=0; j<ATOM_PAIRS.size(); j++)
	{
		
		if (! BOXDIM.IS_RCUT_SAFE(ATOM_PAIRS[j].S_MAXIM, CONTROLS.N_LAYERS))
		{
			
			if (isatty(fileno(stdout)) && RANK == 0)
			{
				#if WARN == TRUE
					cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "WARNING: ";
				#else
					cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "Error: ";
				#endif
					
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "Outer cutoff greater than half of at least one layered cell vector at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Frame:                      " << i << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Pair type:                  " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	cell vectors (a)            " << BOXDIM.CELL_AX << " " << BOXDIM.CELL_AY << " " << BOXDIM.CELL_AZ << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	cell vectors (b)            " << BOXDIM.CELL_BX << " " << BOXDIM.CELL_BY << " " << BOXDIM.CELL_BZ << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	cell vectors (c)            " << BOXDIM.CELL_CX << " " << BOXDIM.CELL_CY << " " << BOXDIM.CELL_CZ << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Extent (x,y,z)              " << BOXDIM.EXTENT_X << " " << BOXDIM.EXTENT_Y << " " << BOXDIM.EXTENT_Z << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Layers:                     " << CONTROLS.N_LAYERS << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Effective cell vectors (a): " << BOXDIM.CELL_AX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_AY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_AZ * (2*CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Effective cell vectors (a): " << BOXDIM.CELL_BX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_BY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_BZ * (2*CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;
				cout << COUT_STYLE.MAGENTA << COUT_STYLE.BOLD << "	Effective cell vectors (a): " << BOXDIM.CELL_CX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_CY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_CZ * (2*CONTROLS.N_LAYERS +1) << COUT_STYLE.ENDSTYLE << endl;

				
				
				#if WARN == FALSE
					exit_run(0);
				#endif
			}
			else if ( RANK == 0 ) 
			{
				#if WARN == TRUE
					cout << "WARNING: ";
				#else
					cout << "Error: ";
				#endif
				
				cout <<  "Outer cutoff greater than half of at least one layered cell vector at least one box length: "  << ATOM_PAIRS[j].S_MAXIM <<COUT_STYLE.ENDSTYLE << endl;
				cout <<  "	Frame:                      " << i << COUT_STYLE.ENDSTYLE << endl;
				cout <<  "	Pair type:                  " << ATOM_PAIRS[j].ATM1TYP << " " << ATOM_PAIRS[j].ATM2TYP << COUT_STYLE.ENDSTYLE << endl;		
				cout <<  "	cell vectors (a)            " << BOXDIM.CELL_AX << " " << BOXDIM.CELL_AY << " " << BOXDIM.CELL_AZ << endl;
				cout <<  "	cell vectors (b)            " << BOXDIM.CELL_BX << " " << BOXDIM.CELL_BY << " " << BOXDIM.CELL_BZ << endl;
				cout <<  "	cell vectors (c)            " << BOXDIM.CELL_CX << " " << BOXDIM.CELL_CY << " " << BOXDIM.CELL_CZ << endl;
				cout <<  "	Extent (x,y,z)              " << BOXDIM.EXTENT_X << " " << BOXDIM.EXTENT_Y << " " << BOXDIM.EXTENT_Z << endl;				
				cout <<  "	Layers:                     " << CONTROLS.N_LAYERS << endl;
				cout <<  "	Effective cell vectors (a): " << BOXDIM.CELL_AX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_AY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_AZ * (2*CONTROLS.N_LAYERS +1) << endl;
				cout <<  "	Effective cell vectors (a): " << BOXDIM.CELL_BX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_BY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_BZ * (2*CONTROLS.N_LAYERS +1) << endl;
				cout <<  "	Effective cell vectors (a): " << BOXDIM.CELL_CX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_CY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_CZ * (2*CONTROLS.N_LAYERS +1) << endl;	
				
				#if WARN == FALSE
					exit_run(0);
				#endif										
			}
		}
	}
		
	// Setup the trajectory-holding data object
		
	ATOMTYPE    .resize(ATOMS);		
	COORDS      .resize(ATOMS);
	FORCES      .resize(ATOMS);	// Use for read-in forces.
	ACCEL       .resize(ATOMS); // Use for calculated forces in ZCalc_Ewald.
	CHARGES     .resize(ATOMS);
	ATOMTYPE_IDX.resize(ATOMS);
			
	// Read trajectory, convert to proper units, and apply PBC
			
	for (int j=0; j<ATOMS; j++)
	{
		TRAJ_INPUT >> ATOMTYPE[j];

		int k ;
		for( k=0; k<TMP_ATOMTYPE.size(); k++)
		{
			if(ATOMTYPE[j] == TMP_ATOMTYPE[k])
			{
				ATOMTYPE_IDX[j] = k;
				break ;
			}
		}

		if ( k == TMP_ATOMTYPE.size() )
			EXIT_MSG("Did not recognize atom type " + ATOMTYPE[j]) ;
		
		TRAJ_INPUT >> COORDS[j].X;
		TRAJ_INPUT >> COORDS[j].Y;
		TRAJ_INPUT >> COORDS[j].Z;
			
		TRAJ_INPUT >> FORCES[j].X;
		TRAJ_INPUT >> FORCES[j].Y;
		TRAJ_INPUT >> FORCES[j].Z;

		// Convert forces from atomic (H/B) to kcal/mol/Angs (Stillinger's units) ... Note, all atom pairs must be of the same type, so using 0 index is ok.
				
		FORCES[j].X *= 627.50960803*1.889725989;
		FORCES[j].Y *= 627.50960803*1.889725989;
		FORCES[j].Z *= 627.50960803*1.889725989;

						
		if(CONTROLS.WRAP_COORDS)	// Apply PBC (for cases of unwrapped coordinates)
		{
			BOXDIM.WRAP_ATOM(COORDS[j], WRAP_IDX[j], true);
		}			
			
		// Assign atom charges.
		if ( CONTROLS.IF_SUBTRACT_COUL ) 
			for(int ii=0; ii< CONTROLS.NATMTYP; ii++)
				if( ATOMTYPE[j] == ATOM_PAIRS[ii].ATM1TYP )
				{
					CHARGES[j] = ATOM_PAIRS[ii].ATM1CHG;
					break;								
				}
			
	}
		
	// If layering requested, replicate the system

	build_layers(CONTROLS.N_LAYERS);
		
	if(i==0)
	{
		if ( (CONTROLS.N_LAYERS > 0) )	// Then ghost atoms are used 
		{
			if ( RANK == 0 ) 
			{
				cout << "	Reporting outcome of layering for first frame ONLY: " << endl;
				cout << "	Real atoms:                   " << ATOMS << endl;
				cout << "	Total atoms (ghost):          " << ALL_ATOMS << endl;
				std::streamsize p = cout.precision() ;
				cout.precision(6) ;
				cout <<  "	cell vectors (a)            " << BOXDIM.CELL_AX << " " << BOXDIM.CELL_AY << " " << BOXDIM.CELL_AZ << endl;
				cout <<  "	cell vectors (b)            " << BOXDIM.CELL_BX << " " << BOXDIM.CELL_BY << " " << BOXDIM.CELL_BZ << endl;
				cout <<  "	cell vectors (c)            " << BOXDIM.CELL_CX << " " << BOXDIM.CELL_CY << " " << BOXDIM.CELL_CZ << endl;
				cout <<  "	cell volume (A^3)           " << BOXDIM.VOL << endl; 
				cout <<  "	Extent (x,y,z)              " << BOXDIM.EXTENT_X << " " << BOXDIM.EXTENT_Y << " " << BOXDIM.EXTENT_Z << endl;				
				cout <<  "	Layers:                     " << CONTROLS.N_LAYERS << endl;
				cout <<  "	Effective cell vectors (a): " << BOXDIM.CELL_AX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_AY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_AZ * (2*CONTROLS.N_LAYERS +1) << endl;
				cout <<  "	Effective cell vectors (a): " << BOXDIM.CELL_BX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_BY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_BZ * (2*CONTROLS.N_LAYERS +1) << endl;
				cout <<  "	Effective cell vectors (a): " << BOXDIM.CELL_CX * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_CY * (2*CONTROLS.N_LAYERS +1) << " " << BOXDIM.CELL_CZ * (2*CONTROLS.N_LAYERS +1) << endl;	
				cout.precision(p) ;
			}
		}
		else if ( RANK==0 ) // No ghost atoms.
			cout << "WARNING: Ghost atoms/implicit layers are NOT being used." << endl;
	}
}

void FRAME::build_layers(int N_LAYERS)
// Build layers around the primitive unit cell for use in ghost atom calculations.
{
	
	
	int     TEMP_IDX;
	XYZ     TEMP_XYZ;
	XYZ_INT TEMP_LAYER;
	
	ALL_COORDS   .resize(ATOMS);
	LAYER_IDX    .resize(ATOMS);
	WRAP_IDX     .resize(ATOMS);
	ATOMTYPE     .resize(ATOMS) ;
	ATOMTYPE_IDX .resize(ATOMS) ;
	CHARGES      .resize(ATOMS) ;
	PARENT       .resize(ATOMS);

	// Set the first ATOMS atoms of ALL_COORDS equivalent to the "real" COORDS
	// Then wrap those ALL_COORDS into primitive cell
	
	ALL_ATOMS = ATOMS;	// Default setting: for zero layers
		
	for (int a1=0; a1<ATOMS; a1++) 
	{
		BOXDIM.WRAP_ATOM(COORDS[a1], ALL_COORDS[a1], WRAP_IDX[a1], true);
		
		PARENT    [a1] = a1;
		LAYER_IDX [a1].X = LAYER_IDX [a1].Y = LAYER_IDX [a1].Z = 0;
	}

	if(N_LAYERS>0 )
	{	
		TEMP_IDX = ATOMS;

		for(int n1 = -N_LAYERS; n1<=N_LAYERS; n1++)
		{
			for(int n2 = -N_LAYERS; n2<=N_LAYERS; n2++)
			{
				for(int n3 = -N_LAYERS; n3<=N_LAYERS; n3++)
				{	
					if (n1 == 0 && n2 == 0 && n3 == 0 ) 
						continue;
					else
					{		
						for(int a1=0; a1<ATOMS; a1++)
						{
							TEMP_LAYER.X = n1;
							TEMP_LAYER.Y = n2;
							TEMP_LAYER.Z = n3;
							
							BOXDIM.LAYER_ATOM(ALL_COORDS[a1], TEMP_LAYER, TEMP_XYZ);

							ALL_COORDS   .push_back(TEMP_XYZ);
							LAYER_IDX    .push_back(TEMP_LAYER);
							ATOMTYPE     .push_back(ATOMTYPE    [a1]);
							ATOMTYPE_IDX .push_back(ATOMTYPE_IDX[a1]);
							CHARGES      .push_back(CHARGES     [a1]);
							PARENT       .push_back(a1);
				
							TEMP_IDX++;
						}
					}
				}
			}
		}
		ALL_ATOMS = TEMP_IDX;
	}
}


void PAIRS::set_cheby_vals()
// Calculate Chebyshev xmin, xmax, xavg.
{
  Cheby::set_cheby_params(S_MINIM, S_MAXIM, LAMBDA, CHEBY_TYPE, X_MINIM, X_MAXIM, X_DIFF, X_AVG) ;
}

void JOB_CONTROL::LSQ_SETUP(int npairs, int no_atom_types)
// Setup the JOB_CONTROL structure based on inputs parsed for LSQ calculations.
{
		TOT_SHORT_RANGE = TOT_SNUM + NUM_3B_CHEBY + NUM_4B_CHEBY;
			
		// Keep track of the total number of lsq parameters.
		TOT_ALL_PARAMS = TOT_SHORT_RANGE ;

		if ( FIT_COUL )
			TOT_ALL_PARAMS += npairs ;

		if ( FIT_ENER_EVER || FIT_ENER )
			TOT_ALL_PARAMS += no_atom_types ;
	
		if((FIT_STRESS  || FIT_STRESS_ALL) && NSTRESS == -1)
			NSTRESS = NFRAMES;
		
		if(FIT_ENER && NENER == -1)
			NENER = NFRAMES;		
	
		FIT_ENER_EVER = FIT_ENER;	// Is energy ever fit ?

		if (INFILE.size() == 1)
			INFILE_FRAMES.push_back(NFRAMES);
}
