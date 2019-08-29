/*

Functions for printing trajectories to different file formats.

Current supported types include:
1. .gen
2. .xyz(f)
3. .lammpstrj ()
4. .pdb

... Needs access to:

1. SYSTEM (coords, forces, velocities, atom types, boxdim, etc)
2. CONTROLS (step, frequency, timestep)

*/

#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes
#include<string>

using namespace std;

#include "functions.h"
#include "util.h"
#include "io_styles.h"


// Constructor/deconstructor

WRITE_TRAJ::WRITE_TRAJ()
{
	FIRST_CALL = true;
}

WRITE_TRAJ::WRITE_TRAJ(string CONTENTS_STR) // Default extension is .gen
{
	string EXTENSION_STR = "GEN";
	
	INIT(EXTENSION_STR, CONTENTS_STR);
}

WRITE_TRAJ::WRITE_TRAJ(string EXTENSION_STR, string CONTENTS_STR)
{
	INIT(EXTENSION_STR, CONTENTS_STR);
}

WRITE_TRAJ::~WRITE_TRAJ()
{
	TRAJFILE.close();
	
	if (TRAJFRCF.is_open())
		TRAJFRCF.close();
	
	if (TRAJFRCL.is_open())
		TRAJFRCL.close();
}

// Initializer called by constructor

void WRITE_TRAJ::INIT(string EXTENSION_STR, string CONTENTS_STR)
{
	FIRST_CALL = true;

	SET_EXTENSION(EXTENSION_STR);	
	SET_CONTENTS (CONTENTS_STR);
	SET_FILENAME();
	
	TRAJFILE.open(FILENAME.data());

	if (!TRAJFILE.is_open())
	{
		cout << "Failed to open trajfile " << FILENAME << " for writing!" << endl;	
		exit_run(0);
	}

	if (EXTENSION == TRAJ_EXT::XYZF_FORCE)
	{ 
		TRAJFRCL.open("forceout-labeled.txt");
		
		if (!TRAJFRCL.is_open())
		{
			cout << " Failed to open associated file forceout-labeled.txt for writing!" << endl;
			exit_run(0);
		}
		
		TRAJFRCF.open("forceout.txt");
		
		if (!TRAJFRCL.is_open())
		{
			cout << " Failed to open associated file forceout.txt for writing!" << endl;	
			exit_run(0);
		}	
	}
}

//  Setup the enum types for file extension and file type

void WRITE_TRAJ::SET_EXTENSION(string EXTENSION_STR)
{	
	if(EXTENSION_STR == "GEN")
		EXTENSION = TRAJ_EXT::GEN;
	else if(EXTENSION_STR == "XYZ")
		EXTENSION = TRAJ_EXT::XYZ;
	else if(EXTENSION_STR == "XYZF_FORCE")
		EXTENSION = TRAJ_EXT::XYZF_FORCE;		
	else if(EXTENSION_STR == "LAMMPSTRJ")
		EXTENSION = TRAJ_EXT::LAMMPSTRJ;
	else if(EXTENSION_STR == "PDB")
	{
		EXTENSION = TRAJ_EXT::PDB;	
		cout << "ERROR: File type .pdb has not been implemented yet." << endl;
		exit_run(0);
	}
	else
	{
		cout << "ERROR: Unknown extension type for output trajectory." << endl;
		cout << "Allowed types are GEN, XYZ, XYZF_FORCE, and LAMMPSTRJ." << endl;
		exit_run(0);	
	}
}

void WRITE_TRAJ::SET_CONTENTS(string CONTENTS_STR)
{	
	if(CONTENTS_STR == "STANDARD")
		CONTENTS = TRAJ_TYPE::STANDARD;
	else if(CONTENTS_STR == "BAD_1")
		CONTENTS = TRAJ_TYPE::BAD_1;
	else if(CONTENTS_STR == "BAD_2")
		CONTENTS = TRAJ_TYPE::BAD_2;
	else if(CONTENTS_STR == "BAD_3")
		CONTENTS = TRAJ_TYPE::BAD_3;		
	else if(CONTENTS_STR == "FORCE")
		CONTENTS = TRAJ_TYPE::FORCE;
	else
	{
		cout << "ERROR: Unknown contents type for output trajectory." << endl;
		cout << "Allowed types are STANDARD, BAD_1, BAD_2, BAD_3 and FORCE." << endl;
		exit_run(0);	
	}	
}

void WRITE_TRAJ::SET_FILENAME()
{
	switch (CONTENTS)
	{
		case TRAJ_TYPE::STANDARD:
			FILENAME = "traj";
			break;
		case TRAJ_TYPE::BAD_1:
			FILENAME = "traj_bad_r.lt.rin";
			break;
		case TRAJ_TYPE::BAD_2:
			FILENAME = "traj_bad_r.lt.rin+dp";
			break;
		case TRAJ_TYPE::BAD_3:
			FILENAME = "traj_bad_r.ge.rin+dp_dftbfrq";
			break;
		case TRAJ_TYPE::FORCE:
			FILENAME = "forceout";
			break;			
		default:
			cout << "Error: Unknown traj type " << RETURN_CONTENTS() << " In WRITE_TRAJ::SET_FILENAME()" << endl;
			exit(1);		
	}
	
	switch (EXTENSION)
	{
		case TRAJ_EXT::GEN:
			FILENAME += ".gen";
			break;
	case TRAJ_EXT::XYZ:
		FILENAME += ".xyz";
			break;
		case TRAJ_EXT::XYZF_FORCE:
			FILENAME += ".xyzf";
			break;			
		case TRAJ_EXT::LAMMPSTRJ:
			FILENAME += ".lammpstrj";
			break;
		case TRAJ_EXT::PDB:
			FILENAME += ".pdb";
			break;	
		default:
			cout << "Error: Unknown extension type " << RETURN_EXTENSION() << " In WRITE_TRAJ::SET_FILENAME()" << endl;
			exit(1);		
	}
}

// Set up enum-to-string mapping functions

string WRITE_TRAJ::RETURN_CONTENTS()
{
	string RESULT;
	
	switch (CONTENTS)
	{
		case TRAJ_TYPE::STANDARD:
			RESULT = "STANDARD";
			break;
		case TRAJ_TYPE::BAD_1:
			RESULT = "BAD_1";
			break;
		case TRAJ_TYPE::BAD_2:
			RESULT = "BAD_2";
			break;
		case TRAJ_TYPE::BAD_3:
			RESULT = "BAD_3";
			break;			
		case TRAJ_TYPE::FORCE:
			RESULT = "FORCE";
			break;			
		default:
			RESULT = "Unknown!";		
	}
		
	return RESULT;	
}

string WRITE_TRAJ::RETURN_EXTENSION()
{
	
	string RESULT;
	
	switch (EXTENSION)
	{
		case TRAJ_EXT::GEN:
			RESULT = "GEN";
			break;
		case TRAJ_EXT::XYZ:
			RESULT = "XYZ";
			break;
		case TRAJ_EXT::XYZF_FORCE:
			RESULT = "XYZF_FORCE";
			break;			
		case TRAJ_EXT::LAMMPSTRJ:
			RESULT = "LAMMPSTRJ";
			break;
		case TRAJ_EXT::PDB:
			RESULT = "PDB";
			break;	
		default:
			RESULT = "Unknown!";	
	}
	
	return RESULT;		
}

// Define the public callable print function


void WRITE_TRAJ::PRINT_FRAME(JOB_CONTROL & CONTROLS, FRAME & SYSTEM)
{

	if (FIRST_CALL)
	{
		FIRST_CALL = false;
		
		ATOMTYPS.resize(CONTROLS.NATMTYP);
		
		for (int i=0; i<CONTROLS.NATMTYP; i++)
			ATOMTYPS[i] = CONTROLS.ATOMTYPES[i];

		bool FOUND = false;
		
		for (int i=0; i<SYSTEM.ATOMS; i++) 
		{
			FOUND = false;
		
			for (int j=0; j<ATOMTYPS.size(); j++)
			{
				if (ATOMTYPS[j] == SYSTEM.ATOMTYPE[i])
				{
					FOUND = true;
					break;
				}
			}
				
			if (!FOUND)
			{
				cout << "ERROR: Unknown atom type " << SYSTEM.ATOMTYPE[i] << endl;
				cout << "See WRITE_TRAJ::PRINT_FRAME function." << endl;
			
				exit_run(0);
			}
		}
	}

	
	switch (EXTENSION)
	{
		case TRAJ_EXT::GEN:
			PRINT_FRAME_GEN(CONTROLS,SYSTEM);
			break;
  		case TRAJ_EXT::XYZ:
			PRINT_FRAME_XYZ(CONTROLS,SYSTEM);
			break;
		case TRAJ_EXT::XYZF_FORCE:
			PRINT_FRAME_XYZF_FORCE(CONTROLS,SYSTEM);
			break;			
		case TRAJ_EXT::LAMMPSTRJ:
			PRINT_FRAME_LAMMPSTRJ(CONTROLS,SYSTEM);
			break;
		case TRAJ_EXT::PDB:
			cout << "Error in WRITE_TRAJ::PRINT_FRAME: PDB not implemented yet!" << endl;
			exit_run(0);
			break;	
		default:
			cout << "Error in WRITE_TRAJ::PRINT_FRAME" << endl;
			exit(1);	
	}

}

// Define the private print functions

void WRITE_TRAJ::PRINT_FRAME_GEN(JOB_CONTROL & CONTROLS, FRAME & SYSTEM)
{
	TRAJFILE << setw(5) << right << SYSTEM.ATOMS << " S #Step " << CONTROLS.STEP+1 << " Time " << (CONTROLS.STEP+1) * CONTROLS.DELTA_T_FS << " (fs) Temp " << SYSTEM.TEMPERATURE << " (k)" << endl;

	for(int i=0; i<CONTROLS.NATMTYP; i++)
		TRAJFILE << ATOMTYPS[i] << " ";
	TRAJFILE << endl;
		
	for (int i=0; i<SYSTEM.ATOMS; i++) 
	{
			
		XYZ tmp = SYSTEM.COORDS[i] ;

		if ( CONTROLS.WRAP_COORDS ) 
			SYSTEM.BOXDIM.WRAP_ATOM(tmp, SYSTEM.WRAP_IDX[i], false);	// Wrap into the primitive cell
			
		TRAJFILE << right << setw(4) << i+1 << " " << setw(2) << SYSTEM.ATOMTYPE_IDX[i]+1 << " " 
				<< fixed << setprecision(5) << setw(8) << tmp.X << " "
				<< fixed << setprecision(5) << setw(8) << tmp.Y << " " 	
				<< fixed << setprecision(5) << setw(8) << tmp.Z << endl;    	
	}
	   
	TRAJFILE << fixed << setprecision(5) << setw(8) << 0.0 << " "
			 << fixed << setprecision(5) << setw(8) << 0.0 << " "
			 << fixed << setprecision(5) << setw(8) << 0.0 << endl;
		
	TRAJFILE << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_AX << " "
			 << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_AY << " "
			 << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_AZ << endl;
		
	TRAJFILE << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_BX << " "
			 << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_BY << " "
			 << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_BZ << endl;
		
	TRAJFILE << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_CX << " "
			 << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_CY << " "
			 << fixed << setprecision(5) << setw(8) << SYSTEM.BOXDIM.CELL_CZ << endl;
}

void WRITE_TRAJ::PRINT_FRAME_XYZ(JOB_CONTROL & CONTROLS, FRAME & SYSTEM)
{
	TRAJFILE << SYSTEM.ATOMS << endl;
	
	if (SYSTEM.BOXDIM.IS_ORTHO)
		TRAJFILE << SYSTEM.BOXDIM.CELL_AX << " " << SYSTEM.BOXDIM.CELL_BY << " " << SYSTEM.BOXDIM.CELL_CZ << endl;
	else
		EXIT_MSG("ERROR: PRINT_FRAME_XYZF_FORCE for non-orthorhombic cells not yet implemented");	
		
	for (int a1=0; a1<SYSTEM.ATOMS; a1++) 
	{
		XYZ tmp = SYSTEM.COORDS[a1] ;

		if ( CONTROLS.WRAP_COORDS ) 
			SYSTEM.BOXDIM.WRAP_ATOM(tmp, SYSTEM.WRAP_IDX[a1], false);	// Wrap into the primitive cell
			
		TRAJFILE << right << setw(4) << SYSTEM.ATOMTYPE[a1] << " "  
				<< fixed << setprecision(5) << setw(15) << tmp.X << " "
				<< fixed << setprecision(5) << setw(15) << tmp.Y << " " 	
				<< fixed << setprecision(5) << setw(15) << tmp.Z << endl;    	
	}
}

void WRITE_TRAJ::PRINT_FRAME_LAMMPSTRJ(JOB_CONTROL & CONTROLS, FRAME & SYSTEM)
{
	TRAJFILE << "ITEM: TIMESTEP" << endl;
	TRAJFILE << (CONTROLS.STEP+1) * CONTROLS.DELTA_T_FS << endl;
	TRAJFILE << "ITEM: NUMBER OF ATOMS" << endl;
	TRAJFILE << SYSTEM.ATOMS << endl;
	TRAJFILE << "ITEM: BOX BOUNDS xy xz yz xy xz yz" << endl;
	TRAJFILE << "0.0 " << SYSTEM.BOXDIM.CELL_LX << " " << SYSTEM.BOXDIM.XY << endl;
	TRAJFILE << "0.0 " << SYSTEM.BOXDIM.CELL_LY << " " << SYSTEM.BOXDIM.XZ << endl;
	TRAJFILE << "0.0 " << SYSTEM.BOXDIM.CELL_LZ << " " << SYSTEM.BOXDIM.YZ << endl;
	
	
	// FYI, can get fancy and control what gets printed to the lammps-format files using:
	// CONTROLS.PRINT_VELOC
	// CONTROLS.PRINT_FORCE 
	// ... but for now, just print coords and velocs
	
	TRAJFILE << "ITEM: ATOMS id type element xu yu zu vx vy vz" << endl;
	
	for (int i=0; i<SYSTEM.ATOMS; i++) 
	{
		TRAJFILE 
			<< right << setw(20) << i+1 
			<< right << setw(4 ) << SYSTEM.ATOMTYPE_IDX[i]+1 
			<< right << setw(4 ) << SYSTEM.ATOMTYPE[i]
			<< fixed << setprecision(5) << setw(15) << SYSTEM.COORDS[i].X
			<< fixed << setprecision(5) << setw(15) << SYSTEM.COORDS[i].Y
			<< fixed << setprecision(5) << setw(15) << SYSTEM.COORDS[i].Z
			<< fixed << setprecision(5) << setw(15) << SYSTEM.VELOCITY[i].X
			<< fixed << setprecision(5) << setw(15) << SYSTEM.VELOCITY[i].Y
			<< fixed << setprecision(5) << setw(15) << SYSTEM.VELOCITY[i].Z
			<< endl;
	}	
}

void WRITE_TRAJ::PRINT_FRAME_XYZF_FORCE(JOB_CONTROL & CONTROLS, FRAME & SYSTEM)
{
	// this style output allows to test force-matching against an MD run.
		
	double fconv = Hartree*Bohr;	// convert xyzf forces to Hartree per Bohr

	TRAJFILE.precision(10);	

	TRAJFILE << SYSTEM.ATOMS << endl;
	
	if (SYSTEM.BOXDIM.IS_ORTHO)
		TRAJFILE << SYSTEM.BOXDIM.CELL_AX << " " << SYSTEM.BOXDIM.CELL_BY << " " << SYSTEM.BOXDIM.CELL_CZ << endl;
	else
		EXIT_MSG("ERROR: PRINT_FRAME_XYZF_FORCE for non-orthorhombic cells not yet implemented");

	for(int i=0;i<SYSTEM.ATOMS;i++)
	{
	
	 	TRAJFILE << right << setw(4) << SYSTEM.ATOMTYPE[i] 
	        	 << fixed << setprecision(5) << setw(15) << " " << SYSTEM.COORDS[i].X 
			 << fixed << setprecision(5) << setw(15) << " " << SYSTEM.COORDS[i].Y
			 << fixed << setprecision(5) << setw(15) << " " << SYSTEM.COORDS[i].Z
			 << fixed << setprecision(5) << setw(15) << " " << SYSTEM.ACCEL[i].X / fconv
			 << fixed << setprecision(5) << setw(15) << " " << SYSTEM.ACCEL[i].Y / fconv
			 << fixed << setprecision(5) << setw(15) << " " << SYSTEM.ACCEL[i].Z / fconv << endl ;

		TRAJFRCF << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[i].X << endl;
	  	TRAJFRCF << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[i].Y << endl;
	  	TRAJFRCF << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[i].Z << endl;
			
	  	TRAJFRCL <<  SYSTEM.ATOMTYPE[i] << " " << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[i].X << endl;
	  	TRAJFRCL <<  SYSTEM.ATOMTYPE[i] << " " << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[i].Y << endl;
	  	TRAJFRCL <<  SYSTEM.ATOMTYPE[i] << " " << fixed << setw(13) << setprecision(6) << scientific << SYSTEM.ACCEL[i].Z << endl;

	}
}

