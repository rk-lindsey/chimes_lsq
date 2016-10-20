#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes

#include "functions.h"

using namespace std;




NEIGHBORS::NEIGHBORS()		// (overloaded) class constructor -- if no padding specified, default to 0.3
{
	RCUT_PADDING  = 0.3;
	DISPLACEMENT  = 0.0;
	MAX_CUTOFF    = 0.0;
	MAX_CUTOFF_3B = 0.0;
	FIRST_CALL    = true;
	SECOND_CALL   = true;
	USE           = false;
	
}
NEIGHBORS::~NEIGHBORS(){}	// Deconstructor

void NEIGHBORS::INITIALIZE(FRAME & SYSTEM)		// (overloaded) class constructor -- if no padding specified, default to 0.3
{
	if(USE)
	{
		LIST   .resize(SYSTEM.ATOMS);
		LIST_3B.resize(SYSTEM.ATOMS);
	}
}

void NEIGHBORS::INITIALIZE(FRAME & SYSTEM, double & PAD)	// (overloaded) class constructor -- if padding specified, set to value
{
	INITIALIZE(SYSTEM);
	RCUT_PADDING = PAD;
}

void NEIGHBORS::DO_UPDATE(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS)
{

	XYZ RAB;
	XYZ TMP_BOX;
	
	double MAX = 0;
	double rlen = 0;
	
	if(FIRST_CALL)									// Set up the first dimension of the list 
	{
		LIST   .resize(SYSTEM.ATOMS);	
		LIST_3B.resize(SYSTEM.ATOMS);	
	}

	for(int a1=0; a1<SYSTEM.ATOMS; a1++)
	{
		if(!FIRST_CALL)								// Clear out the second dimension so we can start over again
		{										
			vector<int>().swap(LIST   [a1]);	
			vector<int>().swap(LIST_3B[a1]);
		}
		
		for (int a2=a1+1; a2<SYSTEM.ATOMS; a2++)		// for (int a2=0; a2<SYSTEM.COORDS.size(); a2++)
		{			
			// Get pair distance

			rlen = get_dist(SYSTEM, CONTROLS, RAB, a1, a2);		// Updates RAB!

			// I get the correct result when the list contains all a2 achievable
			// at this point in the function, regardless of update frequency
			// i.e. "if(true)"

			if (rlen < MAX_CUTOFF + RCUT_PADDING)		// Then add it to the atoms neighbor list (2B)
			{
				LIST[a1].push_back(a2);		
				
				if(rlen < MAX_CUTOFF_3B + RCUT_PADDING)	// Then add it to the atoms neighbor list (3B)
					LIST_3B[a1].push_back(a2);		
			}
		}
	}
}

void NEIGHBORS::UPDATE_LIST(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS)
{	
	if(FIRST_CALL)	// Then start from scratch by cycling through all atom (and layer atom) pairs
	{
		DO_UPDATE(SYSTEM, CONTROLS);
		FIRST_CALL = false;	
		
		if(RANK==0)
		{
			cout << "NEIGHBOR LIST MAX CUTOFF IS: " << MAX_CUTOFF << endl;
			cout << "USING PADDING: " << RCUT_PADDING << endl;
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
			RCUT_PADDING = MAX_VEL * 30.0 * CONTROLS.DELTA_T;	// should give a distance in AA

			SECOND_CALL = false;
			
			DO_UPDATE(SYSTEM, CONTROLS);
			
			if(RANK == 0)
				cout << "RANK: " << RANK << " RESET RCUT_PADDING TO: " << fixed << RCUT_PADDING << endl;
		}
		else
		{

			DISPLACEMENT += MAX_VEL * CONTROLS.DELTA_T;
			if(DISPLACEMENT > 0.5*RCUT_PADDING)	// Then update the neighbor list	
			{
				DO_UPDATE(SYSTEM, CONTROLS);
				
				DISPLACEMENT = 0;
				MAX_VEL      = 0;

				if(RANK == 0)
					cout << "RANK: " << RANK << " UPDATING ON STEP: " << CONTROLS.STEP << ", WITH PADDING: " << fixed << setprecision(3) << RCUT_PADDING <<  endl;
			}
		}
	}
}
