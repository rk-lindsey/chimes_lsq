#include<vector>
#include<algorithm>
#include<iostream>
#include<cmath>

using namespace std;
#ifdef USE_MPI
#include <mpi.h>
#endif

#include "util.h"
#include "functions.h"
#include "A_Matrix.h"


A_MAT::A_MAT(): FORCES(), STRESSES(), FRAME_ENERGIES(), ATOM_ENERGIES(), CHARGES()
{
	// Set up A-matrix
	
	data_count  = 0;
	param_count = 0;
}

A_MAT::~A_MAT(){}

void A_MAT::INITIALIZE(JOB_CONTROL &CONTROLS, FRAME& SYSTEM, int NPAIRS, vector<PAIRS> & ATOM_PAIRS)
// Set up the A "matrix"
{
	// Clear out class variables 

	NO_ATOM_TYPES = 0;	 
	ATOM_TYPES	.clear();	 
	NO_ATOMS_OF_TYPE.clear();
	
	FORCES		.clear(); 	 
	STRESSES	.clear();	 
	FRAME_ENERGIES	.clear();	 
	ATOM_ENERGIES	.clear();	 
	CHARGES		.clear();	         
	
	// Initialize everything

	INITIALIZE_NATOMS  (SYSTEM.ATOMS,SYSTEM.ATOMTYPE, ATOM_PAIRS);
		
	INITIALIZE_FORCES  (SYSTEM.ATOMS,CONTROLS.TOT_SHORT_RANGE);
	INITIALIZE_ENERGIES(SYSTEM.ATOMS,CONTROLS.TOT_SHORT_RANGE, CONTROLS.FIT_ENER);
	INITIALIZE_STRESSES(CONTROLS.TOT_SHORT_RANGE, CONTROLS.FIT_STRESS, CONTROLS.FIT_STRESS_ALL);
	INITIALIZE_CHARGES (NPAIRS,SYSTEM.ATOMS);
}

void A_MAT::INITIALIZE_NATOMS  (int ATOMS, vector<string> & FRAME_ATOMTYPES, vector<PAIRS> & ATOM_PAIRS)
{
	
	// Goal: Determine how many atoms of each type are present

	// First task: Identify the exepcted atom types
	
	ATOM_TYPES      .clear();
	NO_ATOMS_OF_TYPE.clear();
	
	for (int j=0; j<ATOM_PAIRS.size(); j++)
	{
		if (ATOM_PAIRS[j].ATM1TYP == ATOM_PAIRS[j].ATM2TYP )
		{
			ATOM_TYPES.push_back(ATOM_PAIRS[j].ATM1TYP);
			NO_ATOMS_OF_TYPE.push_back(0);
		}
	}
	
	NO_ATOM_TYPES = NO_ATOMS_OF_TYPE.size();
	
	// Second task: Count up how many of each type appear in the trajectory frames 

	for (int i=0; i<ATOMS; i++)
	{
		for (int j=0; j<NO_ATOM_TYPES; j++)
		{
			if (ATOM_TYPES[j] == FRAME_ATOMTYPES[i])
				NO_ATOMS_OF_TYPE[j]++;
		}
	}
	
	/*
	//For debugging:
	
	cout << "IINITIALIZATION RESULTS: " << endl;
	
	for (int i=0; i<NO_ATOM_TYPES; i++)
		cout << "	+ " <<  ATOM_TYPES[i] << " " << NO_ATOMS_OF_TYPE[i] << endl;
	*/
	
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

void A_MAT::INITIALIZE_ENERGIES(int ATOMS,int PARAMS, bool FRAME_ENER)
{
	FRAME_ENERGIES.resize(PARAMS);
		
	for (int j=0; j<PARAMS; j++)
		FRAME_ENERGIES[j] = 0.0;
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
		for(int i=0; i<NO_ATOM_TYPES; i++)
		{	
			if (DO_EXCLUDE_1B)
				if(find(EXCLUDE_1B.begin(), EXCLUDE_1B.end(), i) != EXCLUDE_1B.end())
					continue;
		
			if( (item == "FORCE") || (item == "STRESS") )	// OUTFILE << " 0.0";	
				OUTFILE << " " << "0.0";
			else // OUTFILE << " 1.0";
				OUTFILE << " " << NO_ATOMS_OF_TYPE[i];	
		}		
	}
}

void A_MAT::write_natoms(ofstream & OUTFILE)
{
	// For each line that is printed to A.txt/b.txt, write the corresponding
	// number of atoms in that frame. This is used by the active learning driver
	// for weighting
	
	int natoms = 0;
	
	for(int i=0; i<NO_ATOM_TYPES; i++)
		natoms +=  NO_ATOMS_OF_TYPE[i];
				
	OUTFILE << natoms << endl;			
	
}

bool A_MAT::skip_2b(int n)
{
	int tidx,snum;
	bool exclude = false;
	

	for (int m=0; m<EXCLUDE_2B.size(); m++)
	{
		tidx =   EXCLUDE_2B[m];
		snum = N_EXCLUDE_2B[m];

		if( (n>=tidx*snum) && (n<(tidx+1)*snum) )
			exclude = true;

	}
return exclude;	
}

void A_MAT::PRINT_FRAME(	const struct JOB_CONTROL &CONTROLS,
				const class FRAME &SYSTEM,
				const vector<class PAIRS> & ATOM_PAIRS,
				const vector<struct CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS,
				int N)
// Print one frame of the A matrix.
{
	// Determine which trajectory file this frame came from
	
	int my_file     = 0;
	int frame_check = 0;
	
	for(int i=0; i<CONTROLS.INFILE.size(); i++)
	{
		frame_check += CONTROLS.INFILE_FRAMES[i];
	
		if (N < frame_check)
			break;
			
		my_file++;
	}


	bool DO_ENER       = CONTROLS.FIT_ENER_EVER ;

	if ( ! fileb.is_open() )
		EXIT_MSG("FILEB was not open");

	for(int a=0;a<FORCES.size();a++) // Loop over atoms
	{	
		// Print Afile: .../////////////// -- For X
		  
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;
			fileA << FORCES[a][n].X  << "   ";
		}
		if ( CONTROLS.FIT_COUL ) 
			for(int i=0; i<CHARGES.size(); i++) // Loop over pair types, i.e. OO, OH, HH
				fileA << CHARGES[i][a].X << "   ";

		add_col_of_ones("FORCE", DO_ENER, fileA);
		write_natoms(filena);			  

		fileA << endl;	
		  
		// Print Afile: .../////////////// -- For Y
		  
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			
			fileA << FORCES[a][n].Y  << "   ";
		}
		if ( CONTROLS.FIT_COUL ) 
			for(int i=0; i<CHARGES.size(); i++) // Loop over pair types, i.e. OO, OH, HH
				fileA << CHARGES[i][a].Y << "   ";
		add_col_of_ones("FORCE", DO_ENER, fileA);
		write_natoms(filena);				  
		fileA << endl;	


		// Print Afile: .../////////////// -- For Z
		  
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)	// Afile
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;
							
			fileA << FORCES[a][n].Z  << "   ";
		}
		if ( CONTROLS.FIT_COUL ) 
			for(int i=0; i<CHARGES.size(); i++) // Loop over pair types, i.e. OO, OH, HH
				fileA << CHARGES[i][a].Z << "   ";
		add_col_of_ones("FORCE", DO_ENER, fileA);
		write_natoms(filena);				  
		fileA << endl;		
			
		// Print Bfile: ...
			
		{
			fileb << SYSTEM.FORCES[a].X << endl;
			fileb << SYSTEM.FORCES[a].Y << endl;
			fileb << SYSTEM.FORCES[a].Z << endl;
			data_count += 3 ;
			
			fileb.flush() ;
			
			fileb_labeled << CONTROLS.INFILE_FORCE_FLAGS[my_file] << SYSTEM.ATOMTYPE[a] << " " <<  SYSTEM.FORCES[a].X << endl;
			fileb_labeled << CONTROLS.INFILE_FORCE_FLAGS[my_file] << SYSTEM.ATOMTYPE[a] << " " <<  SYSTEM.FORCES[a].Y << endl;
			fileb_labeled << CONTROLS.INFILE_FORCE_FLAGS[my_file] << SYSTEM.ATOMTYPE[a] << " " <<  SYSTEM.FORCES[a].Z << endl;

		}
	}
	if (CONTROLS.FIT_STRESS)
	{

		// Check if we need to exclude some tensor data from the A and b text files.
		if( N < CONTROLS.NSTRESS)
                {

		// Output A.txt 
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;
			fileA << STRESSES[n].XX << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);
		fileA << endl;	
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << STRESSES[n].YY << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);	
		write_natoms(filena);
		fileA << endl;
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << STRESSES[n].ZZ << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);
		fileA << endl;	
			
			
		// Convert from GPa to internal units to match A-matrix elements

		fileb << SYSTEM.STRESS_TENSORS.X/GPa << endl;
		fileb << SYSTEM.STRESS_TENSORS.Y/GPa << endl;
		fileb << SYSTEM.STRESS_TENSORS.Z/GPa << endl;
		data_count += 3 ;
			
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_xx " <<  SYSTEM.STRESS_TENSORS.X/GPa << endl;
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_yy " <<  SYSTEM.STRESS_TENSORS.Y/GPa << endl;
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_zz " <<  SYSTEM.STRESS_TENSORS.Z/GPa << endl;
		}
	}
	else if (CONTROLS.FIT_STRESS_ALL)
	{
		// Check if we need to exclude some tensor data from the A and b text files.
		if( N < CONTROLS.NSTRESS)
                {
						
		// Output A.txt
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << STRESSES[n].XX << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);	
		fileA << endl;
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << STRESSES[n].XY << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);	
		fileA << endl;
		
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << STRESSES[n].XZ << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);	
		fileA << endl;	
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;
			fileA << STRESSES[n].XY << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);	
		fileA << endl;	
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;
			fileA << STRESSES[n].YY << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);	
		fileA << endl;	
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;
			fileA << STRESSES[n].YZ << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);	
		fileA << endl;
			
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;
			fileA << STRESSES[n].XZ << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);
		fileA << endl;
						
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;
			fileA << STRESSES[n].YZ << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);	
		fileA << endl;
		
		for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << STRESSES[n].ZZ << " ";
		}
		add_col_of_ones("STRESS", DO_ENER, fileA);
		write_natoms(filena);	
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
		data_count += 9 ;
			
		// Convert from GPa to internal units to match A-matrix elements
					
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_xx " << SYSTEM.STRESS_TENSORS_X.X/GPa << endl;
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_xy " << SYSTEM.STRESS_TENSORS_X.Y/GPa << endl;
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_xz " << SYSTEM.STRESS_TENSORS_X.Z/GPa << endl;      
	
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_yx " << SYSTEM.STRESS_TENSORS_X.Y/GPa << endl; // Symmetry - this is just Y.X
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_yy " << SYSTEM.STRESS_TENSORS_Y.Y/GPa << endl;
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_yz " << SYSTEM.STRESS_TENSORS_Y.Z/GPa << endl;

		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_zx " << SYSTEM.STRESS_TENSORS_X.Z/GPa << endl; // Symmetry - this is just Z.X
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_zy " << SYSTEM.STRESS_TENSORS_Y.Z/GPa << endl; // Symmetry - this is just Z.Y
		fileb_labeled << CONTROLS.INFILE_STRESS_FLAGS[my_file] << "s_zz " << SYSTEM.STRESS_TENSORS_Z.Z/GPa << endl;
		}
	}
	if(CONTROLS.FIT_ENER)
	{
		// Check if we need to exclude some energy data from the A and b text files.
		if(N < CONTROLS.NENER)
		{
		// Output A.txt 
			
		for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << FRAME_ENERGIES[n] << " ";
		}
		add_col_of_ones("ENERGY", DO_ENER, fileA);	
		write_natoms(filena);			
		fileA << endl;
			
		for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << FRAME_ENERGIES[n] << " ";
		}
		add_col_of_ones("ENERGY", DO_ENER, fileA);				
		write_natoms(filena);
		fileA << endl;
			
		for(int n=0; n<CONTROLS.TOT_SHORT_RANGE; n++)
		{
			if(CONTROLS.HIERARCHICAL_FIT)
				if(skip_2b(n))
					continue;		
			fileA << FRAME_ENERGIES[n] << " ";
		}
		add_col_of_ones("ENERGY", DO_ENER, fileA);				
		write_natoms(filena);
		fileA << endl;						
			
		// Output b.txt stuff
			
		fileb                  << SYSTEM.QM_POT_ENER << endl;
		fileb_labeled << CONTROLS.INFILE_ENERGY_FLAGS[my_file] << "+1 " << SYSTEM.QM_POT_ENER << endl;
			
		fileb                  << SYSTEM.QM_POT_ENER << endl;
		fileb_labeled << CONTROLS.INFILE_ENERGY_FLAGS[my_file] << "+1 " << SYSTEM.QM_POT_ENER << endl;
			
		fileb                  << SYSTEM.QM_POT_ENER << endl;
		fileb_labeled << CONTROLS.INFILE_ENERGY_FLAGS[my_file] << "+1 " << SYSTEM.QM_POT_ENER << endl;
		data_count += 3 ;
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

void A_MAT::PRINT_CONSTRAINTS(	const struct JOB_CONTROL &CONTROLS,
				const vector<struct CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS,
	                        int NPAIRS)
// Print out any charge constraints to the A matrix and b files.
{
	int print_rank = NPROCS - 1 ;
	
	if ( CONTROLS.SPLIT_FILES ) // Keep files for A split for convenient parallel processing.
	{
		// Write out dimensions.
		vector<int> all_data_count(NPROCS) ;
		
		// Get the total number of data entries (all_data_count)
		
		#ifdef USE_MPI
			MPI_Allgather(&data_count, 1, MPI_INT, all_data_count.data(), 1, MPI_INT, MPI_COMM_WORLD) ;
		#else
			all_data_count[0] = data_count ;
		#endif
		
		int j ;
		for (j=0; j<NPROCS; j++)
		{
			if (all_data_count[j] == 0)
			{
				print_rank = j ;
				break ;
			}
		}
		if ( j == NPROCS ) 
			print_rank = NPROCS - 1 ;
	}
	else
	{
		print_rank = NPROCS - 1 ;
	}
		
	if ( CONTROLS.FIT_COUL && RANK == print_rank )
	{
		for(int i=0; i<CHARGE_CONSTRAINTS.size(); i++)
		{
			for(int n=0; n < CONTROLS.TOT_SHORT_RANGE; n++)
				fileA << "0.0 ";
			
			for(int j=0; j< NPAIRS; j++)
				for(int k=0; k<CHARGE_CONSTRAINTS.size()+1; k++) // +1 because we n_constr = npairs-1
					if(CHARGE_CONSTRAINTS[i].PAIRTYPE_IDX[k] == j)
						fileA << CHARGE_CONSTRAINTS[i].CONSTRAINTS[k] << " ";
			
			fileA << endl;	
			
			fileb << CHARGE_CONSTRAINTS[i].FORCE << endl;
			data_count++ ;
		}		
	}
}

void A_MAT::CLEANUP_FILES(bool SPLIT_FILES)
// Close and clean up the output files.
{
	fileA.close();
	fileb.close();
	fileb_labeled.close();
	filena.close();

	// Make sure that every process has closed its files.

#ifdef USE_MPI
		MPI_Barrier(MPI_COMM_WORLD);
#endif

	if ( RANK == 0 ) 
	{
		system("cat b.[0-9]*.txt > b.txt");
		system("rm  b.[0-9]*.txt");
		
		system("cat b-labeled.[0-9]*.txt > b-labeled.txt");
		system("rm  b-labeled.[0-9]*.txt");
		
		system("cat natoms.[0-9]*.txt > natoms.txt");
		system("rm  natoms.[0-9]*.txt");		

		if ( ! SPLIT_FILES ) 
		{
			// Serialize into a single A
			// Could make the SVD program read multiple files.
			system("cat A.[0-9]*.txt > A.txt");
			system("rm A.[0-9]*.txt");
		}
	}

			
	vector<int> all_data_count(NPROCS) ;
		
	// Get the total number of data entries (all_data_count)
		
#ifdef USE_MPI
	MPI_Allgather(&data_count, 1, MPI_INT, all_data_count.data(), 1, MPI_INT, MPI_COMM_WORLD) ;
#else
	all_data_count[0] = data_count ;
#endif
		
	int start=0, end=0, total=0;
		
	for (int i=0; i<RANK; i++) 
		start += all_data_count[i] ;

	end = start + all_data_count[RANK] - 1 ;
		
	for (int i=0; i<NPROCS; i++) 
		total += all_data_count[i] ;

	if ( SPLIT_FILES ) // Keep files for A split for convenient parallel processing.
	{
		// Write out dimensions.


		// Write the number of columns, row to start, end, and total #of rows to the dimension file.
		if ( all_data_count[RANK] > 0 )
		{
			char name[80] ;
			sprintf(name, "dim.%04d.txt", RANK) ;
			ofstream out ;
			out.open(name) ;
			if ( ! out.is_open()  )
				EXIT_MSG("Could not open " + string(name)) ;
			out << param_count-N_EXCL_LT_3B << " " << start << " " << end << " " << total << "\n" ;
			out.close() ;
		} 
		else 
		{
			// If there is no data, the A file was not used.  Delete it.
			char name[80] ;
			sprintf(name, "A.%04d.txt", RANK) ;
			remove(name) ;
		}
	}
	else {
		// Write dimensions of A and b matrix out for fitting programs.
		ofstream out("dim.txt") ;
		if ( ! out.is_open()  ) 
			EXIT_MSG("Could not open dim.txt") ;
		out << param_count-N_EXCL_LT_3B << " " << total << endl ;
		out.close() ;
	}

#ifdef USE_MPI
  	// Add a barrier to prevent possible timeouts.
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	
}
	
void A_MAT::OPEN_FILES(const JOB_CONTROL &CONTROLS)
{
		
	char nameA[80];
	char nameB[80];
	char nameBlab[80];
	char namena[80];

	// Label output files by the processor rank
	sprintf(nameA, "A.%04d.txt", RANK);
	sprintf(nameB, "b.%04d.txt", RANK);
	sprintf(nameBlab, "b-labeled.%04d.txt", RANK);
	sprintf(namena, "natoms.%04d.txt", RANK);

	fileA.open(nameA);
	fileb.open(nameB);
	fileb_labeled.open(nameBlab);
	filena.open(namena);

	if ( ! fileA.good() || ! fileA.is_open() )
		EXIT_MSG(string("Could not open ") + nameA) ;

	if ( ! fileb.good() || ! fileb.is_open() )
		EXIT_MSG(string("Could not open ") + nameB) ;

	if ( ! fileb_labeled.good() || ! fileb_labeled.is_open() ) 
		EXIT_MSG(string("Could not open ") + nameBlab) ;
		
	if ( ! filena.good() || ! filena.is_open() )
		EXIT_MSG(string("Could not open ") + namena) ;		


	fileA.precision(16);	//  Reduced precision to 6 for code testing.
	fileA << std::scientific;

	fileb.precision(16);	//  Usual precision set to 16.
	fileb << std::scientific;

	param_count = CONTROLS.TOT_ALL_PARAMS ;

}
