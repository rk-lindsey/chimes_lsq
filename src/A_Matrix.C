#include<vector>
#include<algorithm>
#include<iostream>

using namespace std;
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "util.h"
#include "functions.h"
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

void A_MAT::PRINT_ALL(const struct JOB_CONTROL &CONTROLS,
											const vector<class FRAME> &TRAJECTORY,
											const vector<class PAIRS> & ATOM_PAIRS,
											const vector<struct CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS,
											int istart, int iend)
{
		
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
		
	char nameA[20];
	char nameB[20];
	char nameBlab[20];

	// Label output files by the processor rank
	sprintf(nameA, "A.%04d.txt", RANK);
	sprintf(nameB, "b.%04d.txt", RANK);
	sprintf(nameBlab, "b-labeled.%04d.txt", RANK);

	ofstream fileA(nameA);
	ofstream fileb(nameB);
	ofstream fileb_labeled(nameBlab);

	fileA.precision(16);	//  Reduced precision to 6 for code testing.
	fileA << std::scientific;

	fileb.precision(16);	//  Usual precision set to 16.
	fileb << std::scientific;

	if ( RANK == 0 ) cout << "	...A matrix length (forces): " << FORCES.size()*FORCES[0].size()*3 << endl << endl;

	bool DO_ENER       = CONTROLS.FIT_ENER;
	  
	// Print only the assigned frames.
	for(int N= istart; N <= iend;N++) // Loop over frames
	{
		for(int a=0;a<FORCES[N].size();a++) // Loop over atoms
		{	
			// Print Afile: .../////////////// -- For X
		  
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
				fileA << FORCES[N][a][n].X  << "   ";
			if ( CONTROLS.FIT_COUL ) 
				for(int i=0; i<CHARGES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << CHARGES[N][i][a].X << "   ";
			if ( CONTROLS.FIT_POVER ) 
				fileA << " " << OVERBONDING[N][a].X;

			add_col_of_ones("FORCE", DO_ENER, fileA);			  

			fileA << endl;	
		  
			// Print Afile: .../////////////// -- For Y
		  
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
				fileA << FORCES[N][a][n].Y  << "   ";
			if ( CONTROLS.FIT_COUL ) 
				for(int i=0; i<CHARGES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << CHARGES[N][i][a].Y << "   ";
			if ( CONTROLS.FIT_POVER ) 
				fileA << " " << OVERBONDING[N][a].Y;
			add_col_of_ones("FORCE", DO_ENER, fileA);				  
			fileA << endl;	


			// Print Afile: .../////////////// -- For Z
		  
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
				fileA << FORCES[N][a][n].Z  << "   ";
			if ( CONTROLS.FIT_COUL ) 
				for(int i=0; i<CHARGES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << CHARGES[N][i][a].Z << "   ";
			if ( CONTROLS.FIT_POVER ) 
				fileA << " " << OVERBONDING[N][a].Z;
			add_col_of_ones("FORCE", DO_ENER, fileA);				  
			fileA << endl;		
			
			// Print Bfile: ...
			
			{
				fileb << TRAJECTORY[N].FORCES[a].X << endl;
				fileb << TRAJECTORY[N].FORCES[a].Y << endl;
				fileb << TRAJECTORY[N].FORCES[a].Z << endl;
			
				fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].X << endl;
				fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].Y << endl;
				fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].Z << endl;

			}
		} // End loop over atoms
		
		// Output the stresses and energies 
		
		if (CONTROLS.FIT_STRESS)
		{

			// Check if we need to exclude some tensor data from the A and b text files.
			if( N >= CONTROLS.NSTRESS)
				continue;

			// Output A.txt 
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].XX << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);
			fileA << endl;	
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].YY << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].ZZ << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);
			fileA << endl;	
			
			
			// Convert from GPa to internal units to match A-matrix elements

			fileb << TRAJECTORY[N].STRESS_TENSORS.X/GPa << endl;
			fileb << TRAJECTORY[N].STRESS_TENSORS.Y/GPa << endl;
			fileb << TRAJECTORY[N].STRESS_TENSORS.Z/GPa << endl;
			
			fileb_labeled << "s_xx " <<  TRAJECTORY[N].STRESS_TENSORS.X/GPa << endl;
			fileb_labeled << "s_yy " <<  TRAJECTORY[N].STRESS_TENSORS.Y/GPa << endl;
			fileb_labeled << "s_zz " <<  TRAJECTORY[N].STRESS_TENSORS.Z/GPa << endl;
		
		}
		else if (CONTROLS.FIT_STRESS_ALL)
		{
			// Check if we need to exclude some tensor data from the A and b text files.
			if( N >= CONTROLS.NSTRESS)
				continue;
						
			// Output A.txt
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].XX << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].XY << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].XZ << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;	
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].XY << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;	
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].YY << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;	
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].YZ << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;
			
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].XZ << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);
			fileA << endl;
						
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].YZ << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << STRESSES[N][n].ZZ << " ";
			add_col_of_ones("STRESS", DO_ENER, fileA);	
			fileA << endl;		

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
			
			// Convert from GPa to internal units to match A-matrix elements
					
			fileb_labeled << "s_xx " << TRAJECTORY[N].STRESS_TENSORS_X.X/GPa << endl;
			fileb_labeled << "s_xy " << TRAJECTORY[N].STRESS_TENSORS_X.Y/GPa << endl;
			fileb_labeled << "s_xz " << TRAJECTORY[N].STRESS_TENSORS_X.Z/GPa << endl;	
	
			fileb_labeled << "s_yx " << TRAJECTORY[N].STRESS_TENSORS_X.Y/GPa << endl; // Symmetry - this is just Y.X
			fileb_labeled << "s_yy " << TRAJECTORY[N].STRESS_TENSORS_Y.Y/GPa << endl;
			fileb_labeled << "s_yz " << TRAJECTORY[N].STRESS_TENSORS_Y.Z/GPa << endl;

			fileb_labeled << "s_zx " << TRAJECTORY[N].STRESS_TENSORS_X.Z/GPa << endl; // Symmetry - this is just Z.X
			fileb_labeled << "s_zy " << TRAJECTORY[N].STRESS_TENSORS_Y.Z/GPa << endl; // Symmetry - this is just Z.Y
			fileb_labeled << "s_zz " << TRAJECTORY[N].STRESS_TENSORS_Z.Z/GPa << endl;
		
		}
		if(CONTROLS.FIT_ENER)
		{
			// Check if we need to exclude some energy data from the A and b text files.
			if(N >= CONTROLS.NENER)
				continue;
		
			// Output A.txt 
			
			for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << FRAME_ENERGIES[N][n] << " ";
			add_col_of_ones("ENERGY", DO_ENER, fileA);				
			fileA << endl;
			
			for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << FRAME_ENERGIES[N][n] << " ";
			add_col_of_ones("ENERGY", DO_ENER, fileA);				
			fileA << endl;
			
			for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << FRAME_ENERGIES[N][n] << " ";
			add_col_of_ones("ENERGY", DO_ENER, fileA);				
			fileA << endl;						
			
			// Output b.txt stuff
			
			fileb                  << TRAJECTORY[N].QM_POT_ENER << endl;
			fileb_labeled << "+1 " << TRAJECTORY[N].QM_POT_ENER << endl;
			
			fileb                  << TRAJECTORY[N].QM_POT_ENER << endl;
			fileb_labeled << "+1 " << TRAJECTORY[N].QM_POT_ENER << endl;
			
			fileb                  << TRAJECTORY[N].QM_POT_ENER << endl;
			fileb_labeled << "+1 " << TRAJECTORY[N].QM_POT_ENER << endl;						
		}
		else if(CONTROLS.FIT_ENER_PER_ATOM)
		{
			// Check if we need to exclude some energy data from the A and b text files.
			if(N >= CONTROLS.NENER)
				continue;
		
			// Output A.txt 
			
			for(int a=0; a<ATOM_ENERGIES[N].size(); a++)
			{
				for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
					fileA << ATOM_ENERGIES[N][a][n] << " ";
				add_col_of_ones("ENERGY", DO_ENER, fileA);
				fileA << endl;

				// Output b.txt stuff
				
				fileb                  << TRAJECTORY[N].QM_POT_ENER_PER_ATOM[a] << endl;
				fileb_labeled << "+1 " << TRAJECTORY[N].QM_POT_ENER_PER_ATOM[a] << endl;
			}
		}  	
    }
	
	// A and B file charge constraints: generalized
	// 
	// The order that charge constraints are printed needs to match
	// the order atom pairs are expected...


	if ( CONTROLS.FIT_COUL && (CHARGE_CONSTRAINTS.size()>0) && (RANK == NPROCS - 1)) 
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
		MPI_Barrier(MPI_COMM_WORLD);
	#endif

	if ( RANK == 0 ) {
		// Serialize into a single A and b file for now.
		// Could make the SVD program read multiple files.
		system("cat A.[0-9]*.txt > A.txt");
		system("rm A.[0-9]*.txt");
		system("cat b.[0-9]*.txt > b.txt");
		system("rm b.[0-9]*.txt");
		system("cat b-labeled.[0-9]*.txt > b-labeled.txt");
		system("rm b-labeled.[0-9]*.txt");
	}

}


void A_MAT::add_col_of_ones(string item, bool DO_ENER, ofstream & OUTFILE)
{
	// Determine if: Energies are being included in the fit
	// If so, is this A-matrix row is for a force or stress (then print an additional a 0.0)
	// or if it is for an energy (then print an additional 1.0)?
	
	if (DO_ENER)
	{
		if( (item == "FORCE") || (item == "STRESS") )	// OUTFILE << " 0.0";	
			for(int i=0; i<NO_ATOM_TYPES; i++)
				OUTFILE << " " << "0.0";
		else // OUTFILE << " 1.0";
			for(int i=0; i<NO_ATOM_TYPES; i++)
				OUTFILE << " " << NO_ATOMS_OF_TYPE[i];			
	}
}
