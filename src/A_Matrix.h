#ifndef _A_MATRIX_H
#define _A_MATRIX_H

#include <vector>

using namespace std;

#include "util.h"

struct STENSOR
{
	double XX;
	double YY;
	double ZZ;
	double XY;
	double XZ;
	double YZ;
};

class A_MAT
{

	// Quick note on the structure of A_MATRIX's last column:
	// Suppose we have 3 pair types, O--O, O--H, and H--H
	// Suppose also that the pair type potentials have 4 parameters each
	// When laid out as a vector, any given final column in A_MATRIX will be arranged as:
	// 
	// { (Param-1_O--O), (Param-2_O--O), (Param-3_O--O), (Param-4_O--O), 
	//   (Param-1_O--H), (Param-2_O--H), (Param-3_O--H), (Param-4_O--H), 
	//   (Param-1_H--H), (Param-2_H--H), (Param-3_H--H), (Param-4_H--H)	}
	
	
	public:
	
	int                     NO_ATOM_TYPES;	 // How many atom types are in the frame
	vector<string>          ATOM_TYPES;	 // What are their chemical symbols
	vector<int>             NO_ATOMS_OF_TYPE;// How many atoms of each type are there?
	
	vector<vector<XYZ> >	FORCES; 	// originally "A_MATRIX"       ... [#frames][#atoms][#fittingparameters]
	vector<STENSOR>         STRESSES;	// new; stresses       ... originally tacked on to "A_MATRIX" ... [#frames][#fittingparameters]
	vector<double>          FRAME_ENERGIES;	// new; frame energies ... originally tacked on to "A_MATRIX" ... [#frames][#fittingparameters]
	vector<vector<double> > ATOM_ENERGIES;	// new; per-atom energies      ... [#frames][#atoms] [#fittingparameters]
	
	// Depreciated for now (on the "to-do" list): 
	
	vector<vector<XYZ> >   CHARGES;	        // originally "COULOMB_FORCES" ... [#frames][#pairtypes][#atoms]

	ofstream fileA, fileb, fileb_labeled, filena;
	
	bool			DO_EXCLUDE_1B;	// Are 1-body interactions being excluded?
	bool			DO_EXCLUDE_2B;	// Are 2-body interactions being excluded?
	
	vector<int>		EXCLUDE_1B;	// Indices for atom type to exclude from per-atom energy offset fitting column (Needed for hierarchical fitting support
	vector<int>		EXCLUDE_2B;	// indices for in-memory A-matrix positions of 2-body parameter types to exclude from printing
	vector<int>		N_EXCLUDE_2B;	// Number of parameters for each excluded 2b
	int			N_EXCL_LT_3B;	// Total number of A-matrix columns excluded due to EXCLUDE_1B and EXCLUDE_2B (needed to property print dim file
	

	A_MAT();
	~A_MAT();
	
	void INITIALIZE_NATOMS  (int ATOMS, vector<string> & FRAME_ATOMTYPES, vector<PAIRS> & ATOM_PAIRS);
	void INITIALIZE_FORCES  (int ATOMS, int PARAMS);
	void INITIALIZE_ENERGIES(int ATOMS, int PARAMS, bool FRAME_ENER);
	void INITIALIZE_STRESSES(int PARAMS, bool DIAG_STRESS, bool ALL_STRESS);
	void INITIALIZE_CHARGES (int FF_PAIRS,int ATOMS);
	void PRINT_FRAME(const struct JOB_CONTROL &CONTROLS, const class FRAME &SYSTEM, const vector<class PAIRS> & ATOM_PAIRS, const vector<struct CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, int N);
	void PRINT_CONSTRAINTS(const struct JOB_CONTROL &CONTROLS,const vector<struct CHARGE_CONSTRAINT> & CHARGE_CONSTRAINTS, int NPAIRS);
	void CLEANUP_FILES(bool SPLIT_FILES);
	void OPEN_FILES(const JOB_CONTROL &CONTROLS);
	void INITIALIZE(JOB_CONTROL &CONTROLS, FRAME& SYSTEM, int NPAIRS, vector<PAIRS> & ATOM_PAIRS);
	bool skip_2b(int n);
	
	private:
	
	void add_col_of_ones(string item, bool DO_ENER, ofstream & OUTFILE);
	void write_natoms(ofstream & OUTFILE);
	int data_count;
	int param_count;
	
};


#endif
	
	
	
	
	
	
	
