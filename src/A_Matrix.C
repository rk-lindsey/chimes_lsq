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


A_MAT::A_MAT(): FORCES(), CHARGES(), OVERBONDING(), STRESSES(), FRAME_ENERGIES(), ATOM_ENERGIES()
{
	// Set up A-matrix
	
//	FORCES        .resize(MY_SIZE); 
//	CHARGES       .resize(MY_SIZE);
//	OVERBONDING   .resize(MY_SIZE);
//	STRESSES      .resize(MY_SIZE);
//	FRAME_ENERGIES.resize(MY_SIZE);
//	ATOM_ENERGIES .resize(MY_SIZE);
	
}
A_MAT::~A_MAT(){}

void A_MAT::INITIALIZE(JOB_CONTROL &CONTROLS, FRAME& SYSTEM, int NPAIRS)
// Set up the A "matrix"
{
	INITIALIZE_NATOMS  (SYSTEM.ATOMS,SYSTEM.ATOMTYPE);
		
	INITIALIZE_FORCES  (SYSTEM.ATOMS,CONTROLS.TOT_SHORT_RANGE);
	INITIALIZE_ENERGIES(SYSTEM.ATOMS,CONTROLS.TOT_SHORT_RANGE, CONTROLS.FIT_ENER, CONTROLS.FIT_ENER_PER_ATOM);
	INITIALIZE_STRESSES(CONTROLS.TOT_SHORT_RANGE, CONTROLS.FIT_STRESS, CONTROLS.FIT_STRESS_ALL);
	INITIALIZE_OVERBOND(SYSTEM.ATOMS);
	INITIALIZE_CHARGES (NPAIRS,SYSTEM.ATOMS);		
}

void A_MAT::INITIALIZE_NATOMS  (int ATOMS, vector<string> & FRAME_ATOMTYPES)
{
	// Goal: Determine how many atoms of each type are present

	NO_ATOM_TYPES = 0;
	NO_ATOMS_OF_TYPE.resize(0) ;

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


void A_MAT::INITIALIZE_FORCES(int ATOMS, int NPARAM)
{
	FORCES.resize(ATOMS);
	
	for (int i=0; i<ATOMS; i++)
	{
		FORCES[i].resize(NPARAM);
	
		for (int j=0; j<NPARAM; j++)
		{
			FORCES[i][j].X = 0.0;
			FORCES[i][j].Y = 0.0;
			FORCES[i][j].Z = 0.0;
		}
	}
}

void A_MAT::INITIALIZE_ENERGIES(int ATOMS,int PARAMS, bool FRAME_ENER, bool ATOM_ENER)
{
	if (FRAME_ENER)
	{
		FRAME_ENERGIES.resize(PARAMS);
		
		for (int j=0; j<PARAMS; j++)
			FRAME_ENERGIES[j] = 0.0;
	}
	else if (ATOM_ENER)
	{
		ATOM_ENERGIES.resize(ATOMS);

		for (int i=0; i<ATOMS; i++)
		{
			ATOM_ENERGIES[i].resize(PARAMS);
		
			for (int j=0; j<PARAMS; j++)
				ATOM_ENERGIES[i][j] = 0.0;
		}
	}
}	
	
void A_MAT::INITIALIZE_STRESSES(int PARAMS, bool DIAG_STRESS, bool ALL_STRESS)
{	
	if(DIAG_STRESS || ALL_STRESS)
	{
		STRESSES.resize(PARAMS);
		
		for (int j=0; j<PARAMS; j++)
		{
			STRESSES[j].XX = 0.0;
			STRESSES[j].YY = 0.0;
			STRESSES[j].ZZ = 0.0;
			STRESSES[j].XY = 0.0;
			STRESSES[j].XZ = 0.0;
			STRESSES[j].YZ = 0.0;
		}
	}	
}
	
void A_MAT::INITIALIZE_OVERBOND(int ATOMS)
{
	OVERBONDING.resize(ATOMS);
	
	for (int i=0; i<ATOMS; i++)
	{
		OVERBONDING[i].X = 0;
		OVERBONDING[i].Y = 0;
		OVERBONDING[i].Z = 0;						
	}			
}

void A_MAT::INITIALIZE_CHARGES (int FF_PAIRS,int ATOMS)
{
	CHARGES.resize(FF_PAIRS);

	for (int i=0; i<FF_PAIRS; i++)
	{
		CHARGES[i].resize(ATOMS);		
	
		for (int j=0; j<ATOMS; j++)
		{			
			CHARGES[i][j].X = 0;
			CHARGES[i][j].Y = 0;
			CHARGES[i][j].Z = 0;	
		}
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


void A_MAT::PRINT_FRAME(const struct JOB_CONTROL &CONTROLS,
												const class FRAME &SYSTEM,
												const vector<class PAIRS> & ATOM_PAIRS,
												const vector<struct CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS,
	                      int N)
// Print one frame of the A matrix.
{
	bool DO_ENER       = CONTROLS.FIT_ENER_EVER ;

	if ( ! fileb.is_open() ) {
		EXIT_MSG("FILEB was not open") ;
	}
	
	for(int a=0;a<FORCES.size();a++) // Loop over atoms
	{	
		// Print Afile: .../////////////// -- For X
		  
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
			fileA << FORCES[a][n].X  << "   ";
		if ( CONTROLS.FIT_COUL ) 
			for(int i=0; i<CHARGES.size(); i++) // Loop over pair types, i.e. OO, OH, HH
				fileA << CHARGES[i][a].X << "   ";
		if ( CONTROLS.FIT_POVER ) 
			fileA << " " << OVERBONDING[a].X;

		add_col_of_ones("FORCE", DO_ENER, fileA);			  

		fileA << endl;	
		  
		// Print Afile: .../////////////// -- For Y
		  
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
			fileA << FORCES[a][n].Y  << "   ";
		if ( CONTROLS.FIT_COUL ) 
			for(int i=0; i<CHARGES.size(); i++) // Loop over pair types, i.e. OO, OH, HH
				fileA << CHARGES[i][a].Y << "   ";
		if ( CONTROLS.FIT_POVER ) 
			fileA << " " << OVERBONDING[a].Y;
		add_col_of_ones("FORCE", DO_ENER, fileA);				  
		fileA << endl;	


		// Print Afile: .../////////////// -- For Z
		  
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
			fileA << FORCES[a][n].Z  << "   ";
		if ( CONTROLS.FIT_COUL ) 
			for(int i=0; i<CHARGES.size(); i++) // Loop over pair types, i.e. OO, OH, HH
				fileA << CHARGES[i][a].Z << "   ";
		if ( CONTROLS.FIT_POVER ) 
			fileA << " " << OVERBONDING[a].Z;
		add_col_of_ones("FORCE", DO_ENER, fileA);				  
		fileA << endl;		
			
		// Print Bfile: ...
			
		{
			fileb << SYSTEM.FORCES[a].X << endl;
			fileb << SYSTEM.FORCES[a].Y << endl;
			fileb << SYSTEM.FORCES[a].Z << endl;
			fileb.flush() ;
			
			fileb_labeled << SYSTEM.ATOMTYPE[a] << " " <<  SYSTEM.FORCES[a].X << endl;
			fileb_labeled << SYSTEM.ATOMTYPE[a] << " " <<  SYSTEM.FORCES[a].Y << endl;
			fileb_labeled << SYSTEM.ATOMTYPE[a] << " " <<  SYSTEM.FORCES[a].Z << endl;

		}
	}
	if (CONTROLS.FIT_STRESS)
	{

		// Check if we need to exclude some tensor data from the A and b text files.
		if( N >= CONTROLS.NSTRESS)
			return ;

		// Output A.txt 
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].XX << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);
		fileA << endl;	
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].YY << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].ZZ << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);
		fileA << endl;	
			
			
		// Convert from GPa to internal units to match A-matrix elements

		fileb << SYSTEM.STRESS_TENSORS.X/GPa << endl;
		fileb << SYSTEM.STRESS_TENSORS.Y/GPa << endl;
		fileb << SYSTEM.STRESS_TENSORS.Z/GPa << endl;
			
		fileb_labeled << "s_xx " <<  SYSTEM.STRESS_TENSORS.X/GPa << endl;
		fileb_labeled << "s_yy " <<  SYSTEM.STRESS_TENSORS.Y/GPa << endl;
		fileb_labeled << "s_zz " <<  SYSTEM.STRESS_TENSORS.Z/GPa << endl;
		
	}
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		// Check if we need to exclude some tensor data from the A and b text files.
		if( N >= CONTROLS.NSTRESS)
			return ;
						
		// Output A.txt
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].XX << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].XY << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].XZ << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;	
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].XY << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;	
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].YY << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;	
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].YZ << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].XZ << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);
		fileA << endl;
						
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].YZ << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << STRESSES[n].ZZ << " ";
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		fileA << endl;		

		// Account for the symmetry of the off-diagonal (deviatoric) components
			
		fileb << SYSTEM.STRESS_TENSORS_X.X/GPa << endl;
		fileb << SYSTEM.STRESS_TENSORS_X.Y/GPa << endl;
		fileb << SYSTEM.STRESS_TENSORS_X.Z/GPa << endl;
			
		fileb << SYSTEM.STRESS_TENSORS_X.Y/GPa << endl; // Symmetry - this is just Y.X
		fileb << SYSTEM.STRESS_TENSORS_Y.Y/GPa << endl;
		fileb << SYSTEM.STRESS_TENSORS_Y.Z/GPa << endl;
			
		fileb << SYSTEM.STRESS_TENSORS_X.Z/GPa << endl; // Symmetry - this is just Z.X
		fileb << SYSTEM.STRESS_TENSORS_Y.Z/GPa << endl; // Symmetry - this is just Z.Y
		fileb << SYSTEM.STRESS_TENSORS_Z.Z/GPa << endl;			
			
		// Convert from GPa to internal units to match A-matrix elements
					
		fileb_labeled << "s_xx " << SYSTEM.STRESS_TENSORS_X.X/GPa << endl;
		fileb_labeled << "s_xy " << SYSTEM.STRESS_TENSORS_X.Y/GPa << endl;
		fileb_labeled << "s_xz " << SYSTEM.STRESS_TENSORS_X.Z/GPa << endl;	
	
		fileb_labeled << "s_yx " << SYSTEM.STRESS_TENSORS_X.Y/GPa << endl; // Symmetry - this is just Y.X
		fileb_labeled << "s_yy " << SYSTEM.STRESS_TENSORS_Y.Y/GPa << endl;
		fileb_labeled << "s_yz " << SYSTEM.STRESS_TENSORS_Y.Z/GPa << endl;

		fileb_labeled << "s_zx " << SYSTEM.STRESS_TENSORS_X.Z/GPa << endl; // Symmetry - this is just Z.X
		fileb_labeled << "s_zy " << SYSTEM.STRESS_TENSORS_Y.Z/GPa << endl; // Symmetry - this is just Z.Y
		fileb_labeled << "s_zz " << SYSTEM.STRESS_TENSORS_Z.Z/GPa << endl;
		
	}
	if(CONTROLS.FIT_ENER)
	{
		// Check if we need to exclude some energy data from the A and b text files.
		if(N >= CONTROLS.NENER)
			return ;
		
		// Output A.txt 
			
		for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << FRAME_ENERGIES[n] << " ";
		add_col_of_ones("ENERGY", DO_ENER, fileA);				
		fileA << endl;
			
		for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << FRAME_ENERGIES[n] << " ";
		add_col_of_ones("ENERGY", DO_ENER, fileA);				
		fileA << endl;
			
		for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
			fileA << FRAME_ENERGIES[n] << " ";
		add_col_of_ones("ENERGY", DO_ENER, fileA);				
		fileA << endl;						
			
		// Output b.txt stuff
			
		fileb                  << SYSTEM.QM_POT_ENER << endl;
		fileb_labeled << "+1 " << SYSTEM.QM_POT_ENER << endl;
			
		fileb                  << SYSTEM.QM_POT_ENER << endl;
		fileb_labeled << "+1 " << SYSTEM.QM_POT_ENER << endl;
			
		fileb                  << SYSTEM.QM_POT_ENER << endl;
		fileb_labeled << "+1 " << SYSTEM.QM_POT_ENER << endl;						
	}
	else if(CONTROLS.FIT_ENER_PER_ATOM)
	{
		// Check if we need to exclude some energy data from the A and b text files.
		if(N >= CONTROLS.NENER)
			return ;
		
		// Output A.txt 
			
		for(int a=0; a<ATOM_ENERGIES.size(); a++)
		{
			for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << ATOM_ENERGIES[a][n] << " ";
			add_col_of_ones("ENERGY", DO_ENER, fileA);
			fileA << endl;

			// Output b.txt stuff
				
			fileb                  << SYSTEM.QM_POT_ENER_PER_ATOM[a] << endl;
			fileb_labeled << "+1 " << SYSTEM.QM_POT_ENER_PER_ATOM[a] << endl;
		}
	}
	fileA.flush() ;
	fileb.flush() ;
	fileb_labeled.flush() ;

	if ( ! fileA.good() )
		EXIT_MSG("Error in A file") ;

	if ( ! fileb.good() )
		EXIT_MSG("Error in b file") ;

	if ( ! fileb.good() )
		EXIT_MSG("Error in b_labeled file") ;
	
}

void A_MAT::PRINT_CONSTRAINTS(const struct JOB_CONTROL &CONTROLS,
															const vector<struct CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS,
	                            int NPAIRS)
// Print out any charge constraints to the A matrix and b files.
{
	
	if ( CONTROLS.FIT_COUL && (CHARGE_CONSTRAINTS.size()>0) && (RANK == NPROCS - 1)) 
	{
		for(int i=0; i<CHARGE_CONSTRAINTS.size(); i++)
		{
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << "0.0 ";
			
			for(int j=0; j< NPAIRS; j++)
				for(int k=0; k<CHARGE_CONSTRAINTS.size()+1; k++) // +1 because we n_constr = npairs-1
					if(CHARGE_CONSTRAINTS[i].PAIRTYPE_IDX[k] == j)
						fileA << CHARGE_CONSTRAINTS[i].CONSTRAINTS[k] << " ";
			
			if ( CONTROLS.FIT_POVER ) 
				fileA  << " 0.0 ";
			
			fileA << endl;	
			
			fileb << CHARGE_CONSTRAINTS[i].FORCE << endl;	
		}		
	}
}

void A_MAT::CLEANUP_FILES()
// Close and clean up the output files.
{
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
	
void A_MAT::OPEN_FILES()
{
		
	char nameA[80];
	char nameB[80];
	char nameBlab[80];

	// Label output files by the processor rank
	sprintf(nameA, "A.%04d.txt", RANK);
	sprintf(nameB, "b.%04d.txt", RANK);
	sprintf(nameBlab, "b-labeled.%04d.txt", RANK);

	fileA.open(nameA);
	fileb.open(nameB);
	fileb_labeled.open(nameBlab);

	if ( ! fileA.good() || ! fileA.is_open() )
		EXIT_MSG(string("Could not open") + nameA) ;

	if ( ! fileb.good() || ! fileb.is_open() )
		EXIT_MSG(string("Could not open") + nameB) ;

	if ( ! fileb_labeled.good() || ! fileb_labeled.is_open() ) 
		EXIT_MSG(string("Could not open") + nameBlab) ;


//	cout << "Opened " << nameB << endl ;
//	cout << "Opened " << nameBlab << endl ;
//	cout << "Opened " << nameA << endl ;

	fileA.precision(16);	//  Reduced precision to 6 for code testing.
	fileA << std::scientific;

	fileb.precision(16);	//  Usual precision set to 16.
	fileb << std::scientific;

}
