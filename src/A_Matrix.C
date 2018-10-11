#include<vector>
#include<algorithm>
#include<iostream>

using namespace std;

#include "util.h"
#include "A_Matrix.h"

A_MAT::A_MAT(int NFRAMES)
{
	// Set up A-matrix
	
	MY_SIZE = NFRAMES;
	
	FORCES        .resize(MY_SIZE); 
	CHARGES       .resize(MY_SIZE);
	OVERBONDING   .resize(MY_SIZE);
	STRESSES      .resize(MY_SIZE);
	FRAME_ENERGIES.resize(MY_SIZE);
	ATOM_ENERGIES .resize(MY_SIZE);
	
}
A_MAT::~A_MAT(){}

void A_MAT::INITIALIZE_NATOMS  (int ATOMS, vector<string> & FRAME_ATOMTYPES)
{
	// Goal: Determine how many atoms of each type are present

	NO_ATOM_TYPES = 0;

	vector<string>::iterator it;

	for (int i=0; i<ATOMS; i++)
	{
		// Get the location of the current atom type in the "ATOM_TYPES" list
		// ...if it doesn't exist, add it
		
		it = find(ATOM_TYPES.begin(), ATOM_TYPES.end(), FRAME_ATOMTYPES[i]);


		if (it == ATOM_TYPES.end()) // Then the atom type hasn't been added yet
		{
			NO_ATOM_TYPES++;
			ATOM_TYPES.push_back(FRAME_ATOMTYPES[i]);
			NO_ATOMS_OF_TYPE.push_back(1);
		}
		else
		{
			NO_ATOMS_OF_TYPE[distance(ATOM_TYPES.begin(), it)]++;
		}
	}
}


void A_MAT::INITIALIZE_FORCES(int FRAME, int ATOMS, int NPARAM)
{
	FORCES[FRAME].resize(ATOMS);
	
	for (int i=0; i<ATOMS; i++)
	{
		FORCES[FRAME][i].resize(NPARAM);
	
		for (int j=0; j<NPARAM; j++)
		{
			FORCES[FRAME][i][j].X = 0.0;
			FORCES[FRAME][i][j].Y = 0.0;
			FORCES[FRAME][i][j].Z = 0.0;
		}
	}
}

void A_MAT::INITIALIZE_ENERGIES(int FRAME, int ATOMS,int PARAMS, bool FRAME_ENER, bool ATOM_ENER)
{
	if (FRAME_ENER)
	{
		FRAME_ENERGIES[FRAME].resize(PARAMS);
		
		for (int j=0; j<PARAMS; j++)
			FRAME_ENERGIES[FRAME][j] = 0.0;
	}
	else if (ATOM_ENER)
	{
		ATOM_ENERGIES[FRAME].resize(ATOMS);

		for (int i=0; i<ATOMS; i++)
		{
			ATOM_ENERGIES[FRAME][i].resize(PARAMS);
		
			for (int j=0; j<PARAMS; j++)
				ATOM_ENERGIES[FRAME][i][j] = 0.0;
		}
	}
}	
	
void A_MAT::INITIALIZE_STRESSES(int FRAME, int PARAMS, bool DIAG_STRESS, bool ALL_STRESS)
{	
	if(DIAG_STRESS || ALL_STRESS)
	{
		STRESSES[FRAME].resize(PARAMS);
		
		for (int j=0; j<PARAMS; j++)
		{
			STRESSES[FRAME][j].XX = 0.0;
			STRESSES[FRAME][j].YY = 0.0;
			STRESSES[FRAME][j].ZZ = 0.0;
			STRESSES[FRAME][j].XY = 0.0;
			STRESSES[FRAME][j].XZ = 0.0;
			STRESSES[FRAME][j].YZ = 0.0;
		}
	}	
}
	
void A_MAT::INITIALIZE_OVERBOND(int FRAME, int ATOMS)
{
	OVERBONDING[FRAME].resize(ATOMS);
	
	for (int i=0; i<ATOMS; i++)
	{
		OVERBONDING[FRAME][i].X = 0;
		OVERBONDING[FRAME][i].Y = 0;
		OVERBONDING[FRAME][i].Z = 0;						
	}			
}

void A_MAT::INITIALIZE_CHARGES (int FRAME, int FF_PAIRS,int ATOMS)
{
	CHARGES[FRAME].resize(FF_PAIRS);

	for (int i=0; i<FF_PAIRS; i++)
	{
		CHARGES[FRAME][i].resize(ATOMS);		
	
		for (int j=0; j<ATOMS; j++)
		{			
			CHARGES[FRAME][i][j].X = 0;
			CHARGES[FRAME][i][j].Y = 0;
			CHARGES[FRAME][i][j].Z = 0;	
		}
	}
}
