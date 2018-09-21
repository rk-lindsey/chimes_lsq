#ifndef _IO_STYLES_H
#define _IO_STYLES_H


#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes
#include<string>

using namespace std;

#include "functions.h"

enum class TRAJ_EXT 
{
	GEN,
	XYZ,
	XYZF_FORCE,
	LAMMPSTRJ,
	PDB		
};

enum class TRAJ_TYPE 
{
	STANDARD,	// Standard trajectory file
	BAD_1,		// Configs where r_ij < r_cut,in 
	BAD_2,		// Configs where r_ij < r_cut,in +d_penalty	
	FORCE		// Special: for tracking forces
};

class WRITE_TRAJ
{
	public:
		
		// Constructor is overloaded.
		
		WRITE_TRAJ();
		WRITE_TRAJ(string CONTENTS_STR);	
		WRITE_TRAJ(string EXTENSION_STR, string CONTENTS_STR);
		void INIT(string EXTENSION_STR, string CONTENTS_STR);
		~WRITE_TRAJ();
		
		void PRINT_FRAME(JOB_CONTROL & CONTROLS, FRAME & SYSTEM);

	private:
	
		TRAJ_EXT 	EXTENSION;	// File extension
		TRAJ_TYPE 	CONTENTS;	// File type (contents)
		string		FILENAME;	// Name for the trajectory file
		ofstream	TRAJFILE;	// The atual file object
		ofstream	TRAJFRCF;	// A file of forces for the traj - only used for TRAJ_TYPE == "XYZF_FORCE"
		ofstream	TRAJFRCL;	// Same as TRAJFRCF, but forces are labeled by atom
				
		vector<string>	ATOMTYPS;	// A list of atom types. Used for .gen files
		
		void SET_EXTENSION(string EXTENSION_STR);
		void SET_CONTENTS (string CONTENTS_STR);
		void SET_FILENAME();
		
		string RETURN_EXTENSION();
		string RETURN_CONTENTS();
		
		void PRINT_FRAME_GEN		(JOB_CONTROL & CONTROLS, FRAME & SYSTEM);
		void PRINT_FRAME_XYZ		(JOB_CONTROL & CONTROLS, FRAME & SYSTEM);
		void PRINT_FRAME_XYZF_FORCE	(JOB_CONTROL & CONTROLS, FRAME & SYSTEM);
		void PRINT_FRAME_LAMMPSTRJ	(JOB_CONTROL & CONTROLS, FRAME & SYSTEM);
		void PRINT_FRAME_PDB		(JOB_CONTROL & CONTROLS, FRAME & SYSTEM);
};

#endif
