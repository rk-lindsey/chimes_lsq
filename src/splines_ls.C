// System headers

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iomanip>
#include <map>

// User-defined headers

#include "functions.h"


using namespace std;

#ifndef VERBOSITY 
	#define VERBOSITY 1 
#endif
 

static void read_lsq_input(string & INFILE, int & nframes, int & nlayers, bool & fit_coul, bool & coul_consv, bool & if_subtract_coord, bool & if_subtract_coul, bool & fit_pover, int & cheby_order, string & cheby_type, int & cheby_3b_order, int & invr_parms,
 int &  NATMTYP, bool & if_3b_cheby, vector<PAIRS> & ATOM_PAIRS, vector<TRIPLETS> & PAIR_TRIPLETS, bool & WRAPCOORDS, map<string,int> & PAIR_MAP, map<int,string> & PAIR_MAP_REVERSE, map<string,int> & TRIAD_MAP, map<int,string> & TRIAD_MAP_REVERSE, bool & use_partial_charges );

int factorial(int input)
{
	int result = 1;
	for(int i=input; i>0; i--)
		result *= i;
	
	return result;
}

string FULL_FILE_3B;		// Global variables declared as externs in functions.h, and declared in functions.C
string SCAN_FILE_3B;
string SCAN_FILE_2B;

int main()
{
	//////////////////////////////////////////////////
	//
	// Job parameters: Declare and set defaults
	//
	//////////////////////////////////////////////////
	
	
	// BELOW: VARIABLES INTRODUCED BY BECKY TO ALLOW EXTENSION TO MULTIPLE ATOM TYPES
	
	vector<PAIRS> ATOM_PAIRS;		// Will store relevant info regarding atom interaction pair types.. 
	vector<TRIPLETS> PAIR_TRIPLETS;	// Will store relevant info regarding atom interaction triplet types.. i.e. { [OO,OO,OO], [OO,HH,HH], ... }
	vector<FRAME> TRAJECTORY;		// Stores the trajectory information... box dimensions, atom types, coordinates, and forces
	vector<int>   ATOM_PAIR_TYPES;	// Fore use in double loop over atom pairs. Index corresponds to the overall double loop count. 
									// Provides an index that rells you the atom pair's ATOM_PAIRS's type index.. THIS IS FOR LOOPS OF TYPE
									// 	for(int a1=0;a1<nat-1;a1++)	for(int a2=a1+1;a2<nat;a2++)
	vector<int>   ATOM_PAIR_TYPES_ALL;	// THIS IS FOR LOOPS OF TYPE for(int ai=0; ai<nat; ai++), for(int ak=0; ak<nat; ak++)
	
	map<string,int> PAIR_MAP;			// Input is two of any atom type contained in xyzf file, in any order, output is a pair type index
	map<int,string> PAIR_MAP_REVERSE; 	// Input/output the resverse of PAIR_MAP
	map<string,int> TRIAD_MAP;			// Input is three of any ATOM_PAIRS.PRPR_NM pairs, output is a triplet type index
	map<int,string> TRIAD_MAP_REVERSE;	// Input/output is the reverse of TRIAD_MAP
	
	// MAJOR ASSUMPTIONS: 
	// 1. Number of atoms does not change over frames ( should always be true for MD)
	// 2. Atom ordering in trajectory file does not change over frames ( should also always be true)
	
	
	// BELOW: VARIABLES MORE-OR-LESS ORIGINAL TO THE VERSION I FIRST HAD ACCESS TO 
	
	string INFILE;					// Input trajectory file
	bool WRAPCOORDS = false;
	
	int nframes;					// Number of frames in the movie file	
	int nlayers=1;					// supercells adjacent to central cell. 1 is enough because of 8-Ang. limit to spline potential; 
									// electrostatic energy doesn't depend on this number. Setting to 2 gives identical result.
  
	bool ifsubtract_coord = false;	// If true, subtract overcoordination forces.
	bool ifsubtract_coul = false;	// If true, subtract Coulombic forces (for use with fixed charges). 
	bool fit_coul = false;			// If true, fit coulomb parameters.
	bool use_partial_charges = true; // Will there be any charges in the system?
	bool fit_pover = false;			// If true, fit overcoordination parameters.
	bool coul_consv = false;		// If true, constraints will be applied to charge fitting to try to maintain consistency

	int cheby_order = 8; 			// Order of Chebyshev polynomial if used... set to 8 for DFTB Erep polynomial
	string cheby_type;				// How will distance be transformed?
	int cheby_3b_order = 2;			// how many polynomials for 3b cheby?
	int num_cheby_3b = 0;
	
	int invr_parms = 4;				// currently uses 19 parameters per pair type
 	int NATMTYP = 0;

	bool if_3b_cheby = false; 		// ************* WHY ISN'T THIS ITS OWN PAIR_TYPE????!?!?!?!?!?!?!

	//////////////////////////////////////////////////
	//
	// Read and print input to screen
	//
	//////////////////////////////////////////////////
	
	#if VERBOSITY == 1
		cout << endl << "Reading input file..." << endl;
	#endif

	read_lsq_input(INFILE, nframes, nlayers, fit_coul, coul_consv, ifsubtract_coord, ifsubtract_coul, fit_pover, cheby_order, cheby_type, cheby_3b_order, invr_parms, NATMTYP, if_3b_cheby, ATOM_PAIRS, PAIR_TRIPLETS, WRAPCOORDS, PAIR_MAP, PAIR_MAP_REVERSE, TRIAD_MAP, TRIAD_MAP_REVERSE, use_partial_charges );

	#if VERBOSITY == 1
		cout << "...input file read successful: " << endl << endl;
	#endif


	//////////////////////////////////////////////////
	//
	// Set up force and force derivatie vectors
	//
	//////////////////////////////////////////////////	
		 		
		
	// Quick note on the structure of A_MATRIX's last column:
	// Suppose we have 3 pair types, O--O, O--H, and H--H
	// Suppose also that the pair type potentials have 4 parameters each
	// When laid out as a vector, any given final column in A_MATRIX will be arranged as:
	// 
	// { (Param-1_O--O), (Param-2_O--O), (Param-3_O--O), (Param-4_O--O), 
	//   (Param-1_O--H), (Param-2_O--H), (Param-3_O--H), (Param-4_O--H), 
	//   (Param-1_H--H), (Param-2_H--H), (Param-3_H--H), (Param-4_H--H)	}	
		
	vector<vector <vector <XYZ > > > A_MATRIX;	    	// [ # frames ] [ # atoms ]      [ #fitting parameters ]  .. replaces Als
	vector<vector <vector <XYZ > > > COULOMB_FORCES;	// [ # frames ] [ # pair types ] [ # atoms ]              .. replaces Coul_xx
	vector<vector <XYZ > >           P_OVER_FORCES;  	// [ # frames ] [ # atoms ]                               .. replaces Pover
	
	A_MATRIX      .resize(nframes);
	COULOMB_FORCES.resize(nframes);
	P_OVER_FORCES.resize(nframes);

	// Figure out necessary dimensions for the force/force derivative vectors
	
	int tot_snum = 0; 

	for (int i=0; i<ATOM_PAIRS.size(); i++)
	{
		if ( (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) || (ATOM_PAIRS[i].PAIRTYP == "DFTBPOLY") )	
        {
			ATOM_PAIRS[i].SNUM          = cheby_order;
			
			if (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" )
				ATOM_PAIRS[i].SNUM_3B_CHEBY = cheby_3b_order;
        }
		
		else if (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) 	// Set the distance transformation type
			ATOM_PAIRS[i].CHEBY_TYPE          = cheby_type;
			
		else if (ATOM_PAIRS[i].PAIRTYP == "INVRSE_R") 
			ATOM_PAIRS[i].SNUM = invr_parms;

		else // Spline
			ATOM_PAIRS[i].SNUM = (2+floor((ATOM_PAIRS[i].S_MAXIM - ATOM_PAIRS[i].S_MINIM)/ATOM_PAIRS[i].S_DELTA))*2; //2 is for p0/m0/p1/m1.. 
	
		tot_snum += ATOM_PAIRS[i].SNUM;
	}
	
	#if VERBOSITY == 1		
		cout << "The number of two-body non-coulomb parameters is: " << tot_snum <<  endl;
	#endif



	if ( ATOM_PAIRS[0].SNUM_3B_CHEBY > 0 ) // All atoms types must use same potential, so check if 3b order is greater than zero 
	{
		num_cheby_3b = 0;
		
		for(int i=0; i<PAIR_TRIPLETS.size(); i++)
			num_cheby_3b += PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS;

		#if VERBOSITY == 1
			cout << "The number of three-body Chebyshev parameters is: " << num_cheby_3b << endl;
		#endif
	}

	int tot_short_range = tot_snum + num_cheby_3b;	
	

	//////////////////////////////////////////////////
	//
	// Read and store MD trajectory (coords and forces)
	//
	//////////////////////////////////////////////////
		
	// Read in the trajectory
					 
	ifstream TRAJ_INPUT;
	TRAJ_INPUT.open(INFILE.data());
	
	if(!TRAJ_INPUT.is_open())
	{
		cout << "ERROR: Cannot open trajectory file: " << INFILE << endl;
		exit(1);
	}
	
	TRAJECTORY.resize(nframes);
	
	#if VERBOSITY == 1
		cout << endl << "Reading in the trajectory file..." << endl;
	#endif

	for (int i=0; i<nframes; i++)
	{
		// Read in line with the number of atoms
		
		TRAJ_INPUT >> TRAJECTORY[i].ATOMS;
		
		// Read in line with box dimenstions
		
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.X;
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.Y;
		TRAJ_INPUT >> TRAJECTORY[i].BOXDIM.Z;
		
		// Setup the trajectory-holding data object
		
		TRAJECTORY[i].ATOMTYPE.resize(TRAJECTORY[i].ATOMS);		
		TRAJECTORY[i].COORDS  .resize(TRAJECTORY[i].ATOMS);
		TRAJECTORY[i].FORCES  .resize(TRAJECTORY[i].ATOMS);
		TRAJECTORY[i].CHARGES .resize(TRAJECTORY[i].ATOMS);

		// Read trajectory, convert to proper units, and apply PBC
		
		for (int j=0; j<TRAJECTORY[i].ATOMS; j++)
		{
			TRAJ_INPUT >> TRAJECTORY[i].ATOMTYPE[j];
			
			TRAJ_INPUT >> TRAJECTORY[i].COORDS[j].X;
			TRAJ_INPUT >> TRAJECTORY[i].COORDS[j].Y;
			TRAJ_INPUT >> TRAJECTORY[i].COORDS[j].Z;
			
			TRAJ_INPUT >> TRAJECTORY[i].FORCES[j].X;
			TRAJ_INPUT >> TRAJECTORY[i].FORCES[j].Y;
			TRAJ_INPUT >> TRAJECTORY[i].FORCES[j].Z;

			if (ATOM_PAIRS[0].PAIRTYP != "DFTBPOLY") // Convert forces from Hartree/bohr to kcal/mol/Angs (Stillinger's units) ... Note, all atom pairs must be of the same type, so using 0 index is ok.
			{
				TRAJECTORY[i].FORCES[j].X *= 627.50960803*1.889725989;
				TRAJECTORY[i].FORCES[j].Y *= 627.50960803*1.889725989;
				TRAJECTORY[i].FORCES[j].Z *= 627.50960803*1.889725989;
			}
						
			if(WRAPCOORDS)	// Apply PBC (for cases of unwrapped coordinates)
			{
				TRAJECTORY[i].COORDS[j].X -= floor(TRAJECTORY[i].COORDS[j].X/TRAJECTORY[i].BOXDIM.X)*TRAJECTORY[i].BOXDIM.X;
				TRAJECTORY[i].COORDS[j].Y -= floor(TRAJECTORY[i].COORDS[j].Y/TRAJECTORY[i].BOXDIM.Y)*TRAJECTORY[i].BOXDIM.Y;
				TRAJECTORY[i].COORDS[j].Z -= floor(TRAJECTORY[i].COORDS[j].Z/TRAJECTORY[i].BOXDIM.Z)*TRAJECTORY[i].BOXDIM.Z;
			} 
		}
	
	}
	
	#if VERBOSITY == 1
		cout << "...trajectory file read successful: " << endl << endl;
	#endif

	//////////////////////////////////////////////////
	//
	// Setup A matrix, b vector, for least-squares:
	//
	//////////////////////////////////////////////////

	// Setup the A and Coulomb "matricies"
	
	#if VERBOSITY == 1
		cout << "Setting up the matricies for A, Coulomb forces, and overbonding..." << endl;
	#endif

	for (int f=0; f<nframes; f++)
	{	
		// Set up the A "matrix"
	
		A_MATRIX[f].resize(TRAJECTORY[f].ATOMS);

		for (int i=0; i<TRAJECTORY[f].ATOMS; i++)
		{
			A_MATRIX[f][i].resize(tot_short_range); 

			for (int j=0; j<tot_short_range; j++)
			{		
				A_MATRIX[f][i][j].X = 0;
				A_MATRIX[f][i][j].Y = 0;
				A_MATRIX[f][i][j].Z = 0;
			}
		}
		
		// Setup the Coulomb force "matrix"
	 
		COULOMB_FORCES[f].resize(ATOM_PAIRS.size());
		P_OVER_FORCES[f] .resize(TRAJECTORY[f].ATOMS);
		
		for (int i=0; i<ATOM_PAIRS.size(); i++)
		{
			COULOMB_FORCES[f][i].resize(TRAJECTORY[f].ATOMS);		
		
			for (int j=0; j<TRAJECTORY[f].ATOMS; j++)
			{			
				COULOMB_FORCES[f][i][j].X = 0;
				COULOMB_FORCES[f][i][j].Y = 0;
				COULOMB_FORCES[f][i][j].Z = 0;	
				
				if (i==0)	
				{
					P_OVER_FORCES[f][j].X = 0;
					P_OVER_FORCES[f][j].Y = 0;
					P_OVER_FORCES[f][j].Z = 0;						
				}
			}
		}			
	}

	#if VERBOSITY == 1
		cout << "...matrix setup complete: " << endl << endl;
	#endif

	//////////////////////////////////////////////////
	//
	// Generate A matrix, b vector, for least-squares:
	//
	//////////////////////////////////////////////////

	// cout.precision(10);		// WE SHOULD AT LEAST MOVE THIS TO SOMEWHERE MORE REASONABLE.. LIKE THE SECTION WHERE THE OUTPUT IS ACTUALLY PRINTED
	cout.precision(16);
	
	#if VERBOSITY == 1
		cout << "...Populating the matricies for A, Coulomb forces, and overbonding..." << endl;
	#endif

	for(int i=0; i<nframes; i++)
	{
		// Run a few checks to make sure logic is correct
		 
		if(ifsubtract_coord && fit_pover)
		{
			cout << "LOGIC ERROR: Problem with code logic. Both fit_pover and ifsubtract_coord cannot be true." << endl;
			cout << "             if_subtract_coord should only be true if overbonding parameters have been " << endl;
			cout << "             specified, and FITPOVR set false." << endl;
		}
		
		if(ifsubtract_coul && fit_coul)
		{
			cout << "LOGIC ERROR: Problem with code logic. Both fit_coul and ifsubtract_coul cannot be true." << endl;
			cout << "             ifsubtract_coul should only be true if non-zero charges have been specified " << endl;
			cout << "             and FITCOUL set false." << endl;
		}
		
		ZCalc_Deriv(ATOM_PAIRS, PAIR_TRIPLETS, TRAJECTORY[i], A_MATRIX[i], COULOMB_FORCES[i], nlayers, if_3b_cheby, PAIR_MAP, TRIAD_MAP);
	
		if ( ifsubtract_coord ) // Subtract over-coordination forces from force to be output.
			SubtractCoordForces(TRAJECTORY[i], false, P_OVER_FORCES[i],  ATOM_PAIRS, PAIR_MAP);	
		
		if (ifsubtract_coul) 
			ZCalc_Ewald(TRAJECTORY[i], ATOM_PAIRS, PAIR_MAP);

		if ( fit_pover )	// Fit the overcoordination parameter.
			SubtractCoordForces(TRAJECTORY[i], true, P_OVER_FORCES[i],  ATOM_PAIRS, PAIR_MAP);	
		
			
	}	
	
	
	#if VERBOSITY == 1
		cout << "...matrix population complete: "  << endl << endl;
		cout << "Printing matricies..." << endl;
	#endif
	
	//////////////////////////////////////////////////
	//
	// Print out least squares matrix:
	// Ax=b.
	// A printed to "A.txt", b printed to "b.txt"
	// use these files for python SVD routine.
	//
	//////////////////////////////////////////////////
	
	// Lets recap. What does any given row of the A file look like???

	ofstream fileA("A.txt");
	ofstream fileb("b.txt");
	ofstream fileb_labeled("b-labeled.txt");
  
	fileA.precision(16);	//  Reduced precision to 6 for code testing.
	fileb.precision(16);	//  Usual precision set to 16.
	  
	for(int N=0;N<nframes;N++)
    {
	
		for(int a=0;a<TRAJECTORY[N].ATOMS;a++)
		{		
			// Print Afile: .../////////////// -- For X
		  
			for(int n=0; n < tot_short_range; n++)	// Afile
				fileA << A_MATRIX[N][a][n].X  << "   ";
			  
			if ( fit_coul ) 
				for(int i=0; i<COULOMB_FORCES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << COULOMB_FORCES[N][i][a].X << "   ";
			  
			if ( fit_pover ) 
				fileA << " " << P_OVER_FORCES[N][a].X;
			  
			fileA << endl;	
		  
			// Print Afile: .../////////////// -- For Y
		  
			for(int n=0; n < tot_short_range; n++)	// Afile
				fileA << A_MATRIX[N][a][n].Y  << "   ";
			  
			if ( fit_coul ) 
				for(int i=0; i<COULOMB_FORCES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << COULOMB_FORCES[N][i][a].Y << "   ";
			  
			if ( fit_pover ) 
				fileA << " " << P_OVER_FORCES[N][a].Y;
			  
			fileA << endl;	


			// Print Afile: .../////////////// -- For Z
		  
			for(int n=0; n < tot_short_range; n++)	// Afile
				fileA << A_MATRIX[N][a][n].Z  << "   ";
			  
			if ( fit_coul ) 
				for(int i=0; i<COULOMB_FORCES[N].size(); i++) // Loop over pair types, i.e. OO, OH, HH
					fileA << COULOMB_FORCES[N][i][a].Z << "   ";
			  
			if ( fit_pover ) 
				fileA << " " << P_OVER_FORCES[N][a].Z;
			  
			fileA << endl;		
			
			// Print Bfile: ...

			fileb << TRAJECTORY[N].FORCES[a].X << endl;
			fileb << TRAJECTORY[N].FORCES[a].Y << endl;
			fileb << TRAJECTORY[N].FORCES[a].Z << endl;
			
			fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].X << endl;
			fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].Y << endl;
			fileb_labeled << TRAJECTORY[N].ATOMTYPE[a] << " " <<  TRAJECTORY[N].FORCES[a].Z << endl;	 	
			
		  
		}  	
    }

	// A and B file charge conservation (q00 + 4 qOH + qHH = 0 ) .. This is specific for water systems
	  
	if ( fit_coul && coul_consv) 
	{
		#if VERBOSITY == 1
			cout << endl;
			cout << "*********************************************************************************" << endl;
			cout << "*********************************************************************************" << endl;
			cout << endl;
			cout << "WARNING: CHARGE CONSERVATION HAS NOT BEEN GENERALIZED. " << endl;
			cout << "         RESULTS WILL BE INCORRECT IF: " << endl;
			cout << "              1. Atom type other than O and H are used" << endl;
			cout << "              2. Atoms are given in order O and H in the input file, respectively" << endl;
			cout << endl;
			cout << "*********************************************************************************" << endl;
			cout << "*********************************************************************************" << endl;
			cout << endl;
		#endif
		
		// A file...
	  
		for(int n=0; n < tot_short_range; n++)
			fileA << "0.0 ";

		fileA << "1000.0 4000.0 4000.0"; // OO, HH, OH
	  		  
		if ( fit_pover ) 
			fileA << " 0.0 ";

		fileA << endl;
	  
		for(int n=0; n < tot_short_range; n++)	//Afile charge consistency (qOO - 4 qHH = 0)
			fileA << "0.0 ";
	  
		fileA << "1000.0 -4000.0 0000.0"; // OO, HH, OH		
	  
		if ( fit_pover ) 
			fileA  << " 0.0 ";
  
		fileA << endl;	 
		  
		// B file...
		  
		fileb << "0.0" << endl;	
		fileb << "0.0" << endl;
		
//		fileb_labeled << "0.0" << endl;	
//		fileb_labeled << "0.0" << endl;	
		 
	} 	
	   
	fileA.close();
	fileb.close();
	fileb_labeled.close();
	  	  
	//////////////////////////////////////////////////
	//
	// Print out the params file header
	//
	//////////////////////////////////////////////////	   
	
	ofstream header;
	header.open("params.header");
	
	//////////////////////////////////////////////////	   
	// THE NEW WAY
	//////////////////////////////////////////////////	
	
cout << "checks: " << 	fit_coul <<  " " << use_partial_charges << endl;

	if(!fit_coul && !use_partial_charges)
		header << "USECOUL: false" << endl;
	else
		header << "USECOUL: true" << endl;

	if(fit_coul)	
		header << "FITCOUL: true" << endl;
	else
		header << "FITCOUL: false" << endl;
	
	bool use_pover = false;
	for(int i=0; i<ATOM_PAIRS.size(); i++)
		if(ATOM_PAIRS[i].USE_OVRPRMS)
			use_pover = true;
	
	if(use_pover)
		header << "USEPOVR: true" << endl;
	else
		header << "USEPOVR: false" << endl;
	
	if(fit_pover)
		header << "FITPOVR: true" << endl;
	else
		header << "FITPOVR: false" << endl;
	
	if(if_3b_cheby)
		header << "USE3BCH: true" << endl;
	else
		header << "USE3BCH: false" << endl;
	
	header << endl << "PAIRTYP: " << ATOM_PAIRS[0].PAIRTYP << " ";
	
	if     (ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV")
		header << " " << ATOM_PAIRS[0].SNUM << " " << ATOM_PAIRS[0].SNUM_3B_CHEBY << endl;
	else if(ATOM_PAIRS[0].PAIRTYP == "DFTBPOLY")
		header << " " << ATOM_PAIRS[0].SNUM << endl;	
	else if (ATOM_PAIRS[0].PAIRTYP == "INVRSE_R")
		header << " " << invr_parms << endl;
	else
		header << endl;
	
	header << endl << "ATOM TYPES: " << NATMTYP << endl << endl;
	header << "# TYPEIDX #	# ATM_TYP #	# ATMCHRG #	# ATMMASS #" << endl;
	
	for(int i=0; i<NATMTYP; i++)
		header << i << "		" << ATOM_PAIRS[i].ATM1TYP << "		" << ATOM_PAIRS[i].ATM1CHG << "		" << ATOM_PAIRS[i].ATM1MAS << endl;
	
	bool PRINT_OVR = false;
	int NPAIR =  ATOM_PAIRS.size();	
	
	header << endl << "ATOM PAIRS: " << NPAIR << endl << endl;
	
	header << "!# PAIRIDX #	";
	header << "# ATM_TY1 #	";
	header << "# ATM_TY1 #	";
	header << "# S_MINIM #	";
	header << "# S_MAXIM #	";
	header << "# S_DELTA #	";
	header << "# CHBDIST #	";	// how pair distance is transformed in cheby calc
	header << "# MORSE_LAMBDA #" << endl;

	for(int i=0; i<NPAIR; i++)
	{

		header << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
			 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
			 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
			 << setw(16) << left << ATOM_PAIRS[i].S_MINIM
			 << setw(16) << left << ATOM_PAIRS[i].S_MAXIM							 
			 << setw(16) << left << ATOM_PAIRS[i].S_DELTA
			 << setw(16) << left << ATOM_PAIRS[i].CHEBY_TYPE;
		if(ATOM_PAIRS[i].CHEBY_TYPE == "MORSE")
			header << setw(16) << left << ATOM_PAIRS[i].LAMBDA << endl; 
		else
			header << endl;
	}
	

	if(use_pover)
	{
		header << endl;
 		header << "!# PAIRIDX #	";
 		header << "# ATM_TY1 #	";
 		header << "# ATM_TY1 #	";				
		header << "# USEOVRP #	";
		header << "# TO_ATOM #	";
		header << "# P_OVERB #	";
 		header << "# R_0_VAL #	";
 		header << "# P_1_VAL #	";
 		header << "# P_2_VAL #	";
 		header << "# LAMBDA6 #" << endl;	
			
		for(int i=0; i<NPAIR; i++)
		{
			header << "	"
				 << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
				 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
				 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
				 << setw(16) << left << ATOM_PAIRS[i].USE_OVRPRMS;
			if(ATOM_PAIRS[i].USE_OVRPRMS)
				header	<< setw(16) << left << ATOM_PAIRS[i].OVER_TO_ATM;
			else
				header	<< setw(16) << left << "NONE";
				header
				 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[0] 
				 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[1]
				 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[2] 
				 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[3]
				 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[4] << endl;	
		}		
			
	}
	 
	if(!if_3b_cheby)
	{
		header << endl << "ATOM PAIR TRIPLETS: " << 0 << endl << endl;
	}
	else
	{
		header << endl << "ATOM PAIR TRIPLETS: " << PAIR_TRIPLETS.size() << endl << endl;	

		for(int i=0;i<PAIR_TRIPLETS.size(); i++)
		{
			header << "" << PAIR_TRIPLETS[i].TRIPINDX << "  " << PAIR_TRIPLETS[i].ATMPAIR1 << " " << PAIR_TRIPLETS[i].ATMPAIR2 << " " << PAIR_TRIPLETS[i].ATMPAIR3 << ": ";
			header << PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS << " parameters, " << PAIR_TRIPLETS[i].N_ALLOWED_POWERS << " total parameters "<< endl;	
			header << "     index  |  powers  |  equiv index  |  param index  " << endl;
			header << "   ----------------------------------------------------" << endl;	

			for(int j=0; j<PAIR_TRIPLETS[i].ALLOWED_POWERS.size(); j++)
			{
				header << "      " << setw(6) << fixed << left << j << " ";
				header << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].X  << " ";
				header << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Y  << " ";
				header << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Z  << " ";
				header << "       " << setw(8) << PAIR_TRIPLETS[i].EQUIV_INDICIES[j] << " ";
				header << "       " << setw(8) << PAIR_TRIPLETS[i].PARAM_INDICIES[j] << endl; 
	
			}

			header << endl;
		}	 
		
	}

	header << endl;

	header.close();
	
	//////////////////////////////////////////////////
	//
	// Print out the pair/triplet type  map file
	//
	//////////////////////////////////////////////////	  	

	ofstream MAPFILE;
	MAPFILE.open("ff_groups.map");
	
	MAPFILE << "PAIRMAPS: " << PAIR_MAP.size() << endl;
	
	for(map<string,int>::iterator i=PAIR_MAP.begin(); i!=PAIR_MAP.end(); i++)
		MAPFILE << i->second << " " << i->first << endl;

	if(TRIAD_MAP.size() > 0)
	{
		MAPFILE << endl;
		
		MAPFILE << "TRIPMAPS: " << TRIAD_MAP.size() << endl;
		
		for(map<string,int>::iterator i=TRIAD_MAP.begin(); i!=TRIAD_MAP.end(); i++)
			MAPFILE << i->second << " " << i->first << endl;
	}
	MAPFILE.close();

	//////////////////////////////////////////////////
	//
	// Print out the params file header
	//
	//////////////////////////////////////////////////	  

	#if VERBOSITY == 1
		cout << "	Minimum distances between atoms: (Angstr.)" << endl;
 
		for (int k=0; k<NPAIR; k++)
			cout << "		" << k << "	" << ATOM_PAIRS[k].ATM1TYP << " " << ATOM_PAIRS[k].ATM2TYP << "	" << fixed << setprecision(3) << ATOM_PAIRS[k].MIN_FOUND_DIST << endl;

		cout << "...matrix printing complete: " << endl << endl;
	#endif
	
	  
return 0;		  
}




	//////////////////////////////////////////////////
	//
    // Function definitions
	//
	//////////////////////////////////////////////////



// Read program input from the file "splines_ls.in".
static void read_lsq_input(string & INFILE, int & nframes, int & nlayers, bool & fit_coul, bool & coul_consv, bool & if_subtract_coord, bool & if_subtract_coul, bool & fit_pover, int & cheby_order, string & cheby_type, int & cheby_3b_order, int & invr_parms,
 int &  NATMTYP, bool & if_3b_cheby, vector<PAIRS> & ATOM_PAIRS, vector<TRIPLETS> & PAIR_TRIPLETS, bool & WRAPCOORDS, map<string,int> & PAIR_MAP, map<int,string> & PAIR_MAP_REVERSE, map<string,int> & TRIAD_MAP, map<int,string> & TRIAD_MAP_REVERSE, bool & use_partial_charges )
{
	bool   FOUND_END = false;
	string LINE;
	string TEMP_STR;
	PAIRS  TEMP_PAIR;
	int    TEMP_INT;
	string TEMP_TYPE;
	int    NPAIR;
	int    NTRIP;
	double SUM_OF_CHARGES = 0;
	
	while (FOUND_END == false)
	{
		getline(cin,LINE);


		if(LINE.find("# ENDFILE #") != string::npos)
		{

			#if VERBOSITY == 1
			
				if (if_subtract_coord)
				{
					cout << "Special feature: " << endl;
					cout << " Will subtract contributions from user-specified overbonding " << endl;
					cout << " parameters before generating A-matrix." << endl << endl;	
				}					

				if(use_partial_charges)
				{
					cout << "Special feature: " << endl;
					cout << " Will subtract contributions stemming from user-specified " << endl;
					cout << " charges before generating A-matrix" << endl << endl;	
				}	
			
			#endif
				
			FOUND_END = true;
			break;
				
		}
		
		else if(LINE.find("# WRAPTRJ #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			WRAPCOORDS = bool(LINE.data());
			#if VERBOSITY == 1
				cout << "	# WRAPTRJ #: " << WRAPCOORDS << endl;	
			#endif
		}
		
		else if(LINE.find("# TRJFILE #") != string::npos)
		{
			cin >> INFILE; cin.ignore();

			#if VERBOSITY == 1
				cout << "	# TRJFILE #: " << INFILE << endl;	
			#endif
		}			
		
		else if(LINE.find("# NFRAMES #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			nframes = int(atof(LINE.data()));
			
			#if VERBOSITY == 1
				cout << "	# NFRAMES #: " << nframes << endl;			
			#endif
		}
		
		else if(LINE.find("# NLAYERS #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			nlayers = int(atof(LINE.data()));
			
			#if VERBOSITY == 1
				cout << "	# NLAYERS #: " << nlayers << endl;
			#endif
		}
		
		else if(LINE.find("# FITCOUL #") != string::npos)
		{
			cin >> LINE; cin.ignore();

			if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
				fit_coul = true;
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				fit_coul = false;
			else
			{
				cout << endl << "ERROR: # FITCOUL # must be specified as true or false." << endl;
				exit(1);	
			}	
			
			#if VERBOSITY == 1
			
				cout << "	# FITCOUL #: ";
			
				if (fit_coul)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
				
			#endif
		}
		else if(LINE.find("# CNSCOUL #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
			if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
				coul_consv = true;
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				coul_consv = false;
			else
			{
				cout << endl << "ERROR: # FITCOUL # must be specified as true or false." << endl;
				exit(1);	
			}	
			
			#if VERBOSITY == 1
						
				cout << "	# CNSCOUL #: ";
			
				if (coul_consv)
					cout << "true" << endl;				
				else
					cout << "false" << endl;							
			#endif
		}
		else if(LINE.find("# FITPOVR #") != string::npos)
		{
			cin >> LINE; cin.ignore();

			if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
				fit_pover = true;
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				fit_pover = false;
			else
			{
				cout << endl << "ERROR: # FITPOVR # must be specified as true or false." << endl;
				exit(1);	
			}	
			
			#if VERBOSITY == 1
			
				cout << "	# FITPOVR #: ";
			
				if (fit_pover)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
			#endif
		}
		
		else if(LINE.find("# PAIRTYP #") != string::npos)
		{
			cin >> LINE; //cin.ignore();
			
			TEMP_TYPE = LINE;

			if      (LINE != "SPLINE" && LINE != "CHEBYSHEV" && LINE != "DFTBPOLY" && LINE != "INVRSE_R") // I don't think these are supported still:  && LINE != "LJ" && LINE != "STILLIN")
			{
				cout << endl;
				cout << "ERROR: Unrecognized pair type. Acceptable options are:" << endl;
				cout << "SPLINE"   << endl;
				cout << "CHEBYSHEV"    << endl;
				cout << "DFTBPOLY" << endl;
				cout << "INVRSE_R" << endl;
//				cout << "LJ"       << endl;
//				cout << "STILLIN"  << endl;
				exit(1);
				
			}
			
			#if VERBOSITY == 1
			
				cout << "	# PAIRTYP #: " << LINE;
				
				if(TEMP_TYPE != "DFTBPOLY")
					cout << " ....NOTE: Forces reported in units of kcal/(mol.A), potential energy in kcal/mol." << endl;
				else
					cout << " ....NOTE: Forces reported in atomic units." << endl;
			#endif
				
				if(TEMP_TYPE == "INVRSE_R")
				{
					cin >> LINE;
					cin.ignore();
					invr_parms = int(atof(LINE.data()));
					
					#if VERBOSITY == 1
						cout << "	             " << "Will use the following number of inverse-r params: " << invr_parms << endl;
					#endif
					
				}

			if(TEMP_TYPE == "DFTBPOLY" || TEMP_TYPE == "CHEBYSHEV")
			{
				cin >> LINE;
				cheby_order = int(atof(LINE.data()));
				
				
				#if VERBOSITY == 1
					cout << "	             " << "Will use 2-body order: " << cheby_order << endl;
				#endif
				
				if(TEMP_TYPE == "CHEBYSHEV")
				{
					cin >> LINE;
					cheby_3b_order = int(atof(LINE.data()));
					cin.ignore();

					#if VERBOSITY == 1
						cout << "	             " << "Will use 3-body order: " << cheby_3b_order << endl;
					#endif
						
					if(cheby_3b_order>0)
						if_3b_cheby = true;
				}
				else
					cin.ignore();
				
			}
			else
				 cin.ignore();
			
			// Do some error checks
			
			if(if_3b_cheby && nlayers > 1)
			{
				cout << "ERROR: Use of layers is not supported with 3-body Chebyshev potentials." << endl;
				cout << "       Set # NLAYERS # to 1." << endl;
				exit(0); 
			}
			if(if_3b_cheby && fit_pover)
			{
				cout << "ERROR: Overbonding is not compatible with 3-body Chebyshev potentials." << endl;
				cout << "       Set # FITPOVR # false." << endl;
				exit(0);				
			}

			
		}
/*		
		else if(LINE.find("# SUBCRDS #") != string::npos)
		{
			cin >> LINE; cin.ignore();

			if      (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
				if_subtract_coord = true;
			else if (LINE=="false" || LINE=="False" || LINE=="FALSE" || LINE == "F" || LINE == "f")
				if_subtract_coord = false;
			else
			{
				cout << "ERROR: # SUBCRDS # must be specified as true or false." << endl;
				exit(1);	
			}	
			
			#if VERBOSITY == 1
			
				cout << "	# SUBCRDS #: ";
			
				if (if_subtract_coord)
					cout << "true" << endl;				
				else
					cout << "false" << endl;
			#endif

		}
*/
		else if( (TEMP_TYPE == "CHEBYSHEV") && (LINE.find("# CHBTYPE #") != string::npos))
		{
			cin >> LINE; cin.ignore();
			cheby_type = LINE;

			#if VERBOSITY == 1
				cout << "	# CHBTYPE #: " << cheby_type << endl;	
			#endif
		}		
		
		/////////////////////////////////////////////////////////////////////
		// Read the topology part. For now, ignoring index and atom types.. 
		// Assuming given as OO, OH, HH, as code expects... 
		// will need to be fixed later.
		/////////////////////////////////////////////////////////////////////
		
		else if(LINE.find("# NATMTYP #")!= string::npos)
		{
			cin >> LINE; cin.ignore();
			NATMTYP = int(atof(LINE.data()));
			
			#if VERBOSITY == 1
				cout << "	# NATMTYP #: " << NATMTYP << endl;	
			#endif
			
			NPAIR = NATMTYP*(NATMTYP+1)/2;
			ATOM_PAIRS.resize(NPAIR);
			
			NTRIP = factorial(NATMTYP+3-1)/factorial(3)/factorial(NATMTYP-1);
			PAIR_TRIPLETS.resize(NTRIP);
			
		}

		else if(LINE.find("# TYPEIDX #")!= string::npos)
		{
			
			#if VERBOSITY == 1
				cout << "	# TYPEIDX #    # ATM_TYP #    # ATMCHRG #    # ATMMASS #" << endl;
			#endif
			
			// Figure out the number of non-unique pairs

			TEMP_INT = NATMTYP*NATMTYP;
			
			SUM_OF_CHARGES = 0;
			
			for(int i=0; i<NATMTYP; i++)
			{
				
				// Set the first atom pair types to be of type OO, HH, CC, etc...
				
				ATOM_PAIRS[i].PAIRTYP    = TEMP_TYPE;
				ATOM_PAIRS[i].PAIRIDX	 = i; 
				ATOM_PAIRS[i].CHEBY_TYPE = cheby_type;
				
				cin >> LINE >> ATOM_PAIRS[i].ATM1TYP >> LINE;
				ATOM_PAIRS[i].ATM1CHG = double(atof(LINE.data()));
				
				SUM_OF_CHARGES += ATOM_PAIRS[i].ATM1CHG;

				cin >> LINE; cin.ignore();
				ATOM_PAIRS[i].ATM1MAS = double(atof(LINE.data()));
				
				ATOM_PAIRS[i].ATM2TYP = ATOM_PAIRS[i].ATM1TYP;
				ATOM_PAIRS[i].ATM2CHG = ATOM_PAIRS[i].ATM1CHG;
				ATOM_PAIRS[i].ATM2MAS = ATOM_PAIRS[i].ATM1MAS;
				
				ATOM_PAIRS[i].PRPR_NM =      ATOM_PAIRS[i].ATM1TYP;
				ATOM_PAIRS[i].PRPR_NM.append(ATOM_PAIRS[i].ATM2TYP);
				
				#if VERBOSITY == 1
					cout << " 	" << setw(15) << left << i+1 
						 << setw(15) << left << ATOM_PAIRS[i].ATM1TYP 
						 << setw(15) << left << ATOM_PAIRS[i].ATM1CHG
						 << setw(15) << left << ATOM_PAIRS[i].ATM1MAS << endl;
				#endif
			}
			
			if(!fit_coul)
			{
				if (SUM_OF_CHARGES>0)
				{
					use_partial_charges = true;
					if_subtract_coul = true;
				}
				else
					use_partial_charges = false;				
			}
			
			// Set up all possible unique pair types
			
			TEMP_INT = NATMTYP;
			
			for(int i=0; i<NATMTYP; i++)
			{
				for(int j=i+1; j<NATMTYP; j++)
				{	
					ATOM_PAIRS[TEMP_INT].PAIRTYP = ATOM_PAIRS[i].PAIRTYP;
					ATOM_PAIRS[TEMP_INT].PAIRIDX = TEMP_INT; 
					ATOM_PAIRS[TEMP_INT].CHEBY_TYPE = ATOM_PAIRS[i].CHEBY_TYPE;
														
					ATOM_PAIRS[TEMP_INT].ATM1TYP = ATOM_PAIRS[i].ATM1TYP;
					ATOM_PAIRS[TEMP_INT].ATM2TYP = ATOM_PAIRS[j].ATM1TYP;
					
					ATOM_PAIRS[TEMP_INT].ATM1CHG = ATOM_PAIRS[i].ATM1CHG;
					ATOM_PAIRS[TEMP_INT].ATM2CHG = ATOM_PAIRS[j].ATM1CHG;
					
					ATOM_PAIRS[TEMP_INT].ATM1MAS = ATOM_PAIRS[i].ATM1MAS;
					ATOM_PAIRS[TEMP_INT].ATM2MAS = ATOM_PAIRS[j].ATM1MAS;	
					
					ATOM_PAIRS[TEMP_INT].PRPR_NM =      ATOM_PAIRS[TEMP_INT].ATM1TYP;
					ATOM_PAIRS[TEMP_INT].PRPR_NM.append(ATOM_PAIRS[TEMP_INT].ATM2TYP);											
					
					TEMP_INT++;
				}	
			}

			// Set up the maps to account for non-unique pairs
			
			for(int i=0; i<NATMTYP; i++)
			{
				for(int j=0; j<NATMTYP; j++)
				{
					TEMP_STR = ATOM_PAIRS[i].ATM1TYP;
					TEMP_STR.append(ATOM_PAIRS[j].ATM1TYP);
					
					for(int k=0;k<ATOM_PAIRS.size(); k++)
					{
						if((ATOM_PAIRS[i].ATM1TYP == ATOM_PAIRS[k].ATM1TYP && ATOM_PAIRS[j].ATM1TYP == ATOM_PAIRS[k].ATM2TYP)
						 ||(ATOM_PAIRS[j].ATM1TYP == ATOM_PAIRS[k].ATM1TYP && ATOM_PAIRS[i].ATM1TYP == ATOM_PAIRS[k].ATM2TYP))
						{
							TEMP_STR = ATOM_PAIRS[i].ATM1TYP;
							TEMP_STR.append(ATOM_PAIRS[j].ATM2TYP);
							PAIR_MAP.insert(make_pair(TEMP_STR,k));	// Maps the true pair index, k, to the string formed by joining the two atom types.
							PAIR_MAP_REVERSE.insert(make_pair(k,TEMP_STR));
						}
					}
				}
			}
			
			if(if_3b_cheby)
			{
				// Generate unique triplets
			
				TRIPLETS TRIP_ATOMS;		// Each item is an atom type
				TEMP_INT = 0;				// Will hold pair triplet index
			
				for(int i=0; i<NATMTYP; i++)
				{
					for(int j=i; j<NATMTYP; j++)
					{
						for(int k=j; k<NATMTYP; k++)
						{
							// Get the three atoms that will define the 3-body interaction
						
							TRIP_ATOMS.ATMPAIR1 = ATOM_PAIRS[i].ATM1TYP;	// The first N_ATMTYP pairs are of type AA, BB, CC ... N_ATMTYP
							TRIP_ATOMS.ATMPAIR2 = ATOM_PAIRS[j].ATM1TYP;
							TRIP_ATOMS.ATMPAIR3 = ATOM_PAIRS[k].ATM1TYP;

							// Construct the triplet atom pairs from those atoms
						
							PAIR_TRIPLETS[TEMP_INT].TRIPINDX = TEMP_INT;
						
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR1 = TRIP_ATOMS.ATMPAIR1;	// ij
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR2 = TRIP_ATOMS.ATMPAIR1;	// ik
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR3 = TRIP_ATOMS.ATMPAIR2;	// jk

							PAIR_TRIPLETS[TEMP_INT].ATMPAIR1.append(TRIP_ATOMS.ATMPAIR2);	// ij
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR2.append(TRIP_ATOMS.ATMPAIR3);	// ik
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR3.append(TRIP_ATOMS.ATMPAIR3);	// jk		

							// Now save the "proper" (ordered) name of the pair			

							PAIR_TRIPLETS[TEMP_INT].ATMPAIR1 = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TEMP_INT].ATMPAIR1] ].PRPR_NM;
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR2 = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TEMP_INT].ATMPAIR2] ].PRPR_NM;
							PAIR_TRIPLETS[TEMP_INT].ATMPAIR3 = ATOM_PAIRS[ PAIR_MAP[ PAIR_TRIPLETS[TEMP_INT].ATMPAIR3] ].PRPR_NM;

							TEMP_INT++;						
						}
					}
				}			
			
				// Figure out the allowed pair triplet powers. Here are some rules and considerations:
				//
				// 1. Powers start from zero. So, if order is specified to be 2, polynomial powers range
				//    from 0 to n-1, NOT 1 to n!
				//
				// 2. At least two pairs must have non-zero powers for the interaction to truly correspond
				//    to 3-body interactions
				//
				// 3. Non-uniqueness upon power sorting must be taken into consideration. For example, for
				//    the pair (OO, OH, OH), powers of (1,1,0) are identical to (1,0,1)
				// 
				// NOTE: We need to also take into consideration the corresponding parameter multiplicities.
			 
				XYZ_INT UNSORTED_POWERS;
				XYZ_INT SORTED_POWERS;
			
				vector<XYZ_INT> STORED_SORTED_POWERS;//(cheby_3b_order*cheby_3b_order*cheby_3b_order);		// Make this the max possible size... it will be destroyed later anyway.
			
				int TOP, BOT;
				int ITEMS = -1; 
			
				bool STORED = false;
				int  STORED_IDX;
				int  RUNNING_IDX = 0;
			
				for(int i=0; i<NTRIP; i++)
				{
					ITEMS = 0;
					for(int pair1_pow=0; pair1_pow<cheby_3b_order; pair1_pow++)
					{
						for(int pair2_pow=0; pair2_pow<cheby_3b_order; pair2_pow++)
						{
							for(int pair3_pow=0; pair3_pow<cheby_3b_order; pair3_pow++)
							{
								// Check number 1: Are at least two powers greater than 0?
								if( (pair1_pow>0 && pair2_pow>0) || (pair1_pow>0 && pair3_pow>0) || (pair2_pow>0 && pair3_pow>0) )
								{								
									UNSORTED_POWERS.X = pair1_pow;
									UNSORTED_POWERS.Y = pair2_pow;
									UNSORTED_POWERS.Z = pair3_pow;		
								
									// Store all triplet powers that meet the above criteria.
								
									PAIR_TRIPLETS[i].ALLOWED_POWERS.push_back(UNSORTED_POWERS);
								
									// Now, we need to figure out which of these sets of powers is truly unique. This will depend on a few things.
	
									// Case 1: each pair in the current pair triplet is unique ... multiplicities will be one for each item
								
									if( (PAIR_TRIPLETS[i].ATMPAIR1 != PAIR_TRIPLETS[i].ATMPAIR2) && (PAIR_TRIPLETS[i].ATMPAIR1 != PAIR_TRIPLETS[i].ATMPAIR3) && (PAIR_TRIPLETS[i].ATMPAIR2 != PAIR_TRIPLETS[i].ATMPAIR3) )
										PAIR_TRIPLETS[i].EQUIV_INDICIES.push_back(PAIR_TRIPLETS[i].ALLOWED_POWERS.size()-1);

								
									// Case 2: each pair in the current pair triplet is identical... multiplicity should be 3 when all three powers are not the same
								 
									else if( (PAIR_TRIPLETS[i].ATMPAIR1 == PAIR_TRIPLETS[i].ATMPAIR2) && (PAIR_TRIPLETS[i].ATMPAIR1 == PAIR_TRIPLETS[i].ATMPAIR3))
									{
										// Sort the powers in ascending order

										if (pair1_pow>pair2_pow)
										{
											TOP = pair1_pow;
											BOT = pair2_pow;
										}
										else
										{
											TOP = pair2_pow;
											BOT = pair1_pow;
										}
										if(pair3_pow>TOP)
										{
											SORTED_POWERS.X = pair3_pow;
											SORTED_POWERS.Y = TOP;
											SORTED_POWERS.Z = BOT;										
										}
										else if(pair3_pow>BOT)
										{
											SORTED_POWERS.X = TOP;
											SORTED_POWERS.Y = pair3_pow;
											SORTED_POWERS.Z = BOT;											
										}
										else
										{
											SORTED_POWERS.X = TOP;
											SORTED_POWERS.Y = BOT;
											SORTED_POWERS.Z = pair3_pow;											
										}

										// Check if sorted powers have already been saved
									
										STORED = false;
									
										for(int j=0; j<STORED_SORTED_POWERS.size(); j++)
										{
											if( (STORED_SORTED_POWERS[j].X == SORTED_POWERS.X) &&  (STORED_SORTED_POWERS[j].Y == SORTED_POWERS.Y) &&  (STORED_SORTED_POWERS[j].Z == SORTED_POWERS.Z))
											{
												STORED = true;
												STORED_IDX = j;
												break;
											}										
										}
									
										// Save them, if they have not been already
									
										if(!STORED)
										{
											STORED_SORTED_POWERS           .push_back(SORTED_POWERS);
											PAIR_TRIPLETS[i].EQUIV_INDICIES.push_back(PAIR_TRIPLETS[i].ALLOWED_POWERS.size()-1);									
										}
										else
										{
											PAIR_TRIPLETS[i].EQUIV_INDICIES.push_back(STORED_IDX);
										}
										STORED = false;

									} 
								
									// Case 3: two pairs in the current pair triplet are identical... multiplicities will be 2 whenever power differs on the two identical pairs
								
									else
									{
										// Sort the powers as if the pair types were arranged A != B == C

										// Case 3a.: A == B != C
									
										if(PAIR_TRIPLETS[i].ATMPAIR1 == PAIR_TRIPLETS[i].ATMPAIR2)
										{
											SORTED_POWERS.X = pair3_pow;
											SORTED_POWERS.Y = pair1_pow;
											SORTED_POWERS.Z = pair2_pow;										
										}
									
										// Case 3b.: A != B == C
									
										if(PAIR_TRIPLETS[i].ATMPAIR2 == PAIR_TRIPLETS[i].ATMPAIR3)
										{
											SORTED_POWERS.X = pair1_pow;
											SORTED_POWERS.Y = pair2_pow;
											SORTED_POWERS.Z = pair3_pow;										
										}			
									
										// Case 3c.: A == C != B
									
										if(PAIR_TRIPLETS[i].ATMPAIR1 == PAIR_TRIPLETS[i].ATMPAIR3)
										{
											SORTED_POWERS.X = pair2_pow;
											SORTED_POWERS.Y = pair1_pow;
											SORTED_POWERS.Z = pair3_pow;										
										}		
									
										// Now arrange the last two powers in ascending order Z>Y
									
										if(SORTED_POWERS.Y > SORTED_POWERS.Z)
										{
											BOT = SORTED_POWERS.Z;
											SORTED_POWERS.Z = SORTED_POWERS.Y;
											SORTED_POWERS.Y = BOT;
										}
									
										// Finally, check whether these sorted powers have already been saved
									
										STORED = false;
									
										for(int j=0; j<STORED_SORTED_POWERS.size(); j++)
										{
											if( (STORED_SORTED_POWERS[j].X == SORTED_POWERS.X) &&  (STORED_SORTED_POWERS[j].Y == SORTED_POWERS.Y) &&  (STORED_SORTED_POWERS[j].Z == SORTED_POWERS.Z))
											{
												STORED = true;
												STORED_IDX = j;
												break;
											}
										}
									
										// Save them, if they have not been already
									
										if(!STORED)
										{
											STORED_SORTED_POWERS           .push_back(SORTED_POWERS);
											PAIR_TRIPLETS[i].EQUIV_INDICIES.push_back(PAIR_TRIPLETS[i].ALLOWED_POWERS.size()-1);	
										}
										else
										{
											PAIR_TRIPLETS[i].EQUIV_INDICIES.push_back(STORED_IDX);
										}									
										STORED = false;
					
									}
								}
							}
						}
					}
				
					STORED_SORTED_POWERS.clear();
					vector<XYZ_INT>().swap(STORED_SORTED_POWERS);
				
					// Now all that's left to do is set the force field index for each set of powers
				 
					PAIR_TRIPLETS[i].PARAM_INDICIES.resize(PAIR_TRIPLETS[i].EQUIV_INDICIES.size());
				
					PAIR_TRIPLETS[i].PARAM_INDICIES[0] = 0;
					
					bool FOUND_EQV;
					int  USE_SET = 0;
					int  MAX_SET = 0;
					
					for(int set1=1; set1<PAIR_TRIPLETS[i].EQUIV_INDICIES.size(); set1++)
					{
						FOUND_EQV = false;
						
						for(int set2=0; set2<set1; set2++)
						{
							if(PAIR_TRIPLETS[i].EQUIV_INDICIES[set1] == PAIR_TRIPLETS[i].EQUIV_INDICIES[set2])
							{
								FOUND_EQV = true;
								USE_SET   = set2;
								break;
							}
						}
						
						if(FOUND_EQV)
							PAIR_TRIPLETS[i].PARAM_INDICIES[set1] = PAIR_TRIPLETS[i].PARAM_INDICIES[USE_SET];
						else
						{
							MAX_SET++;
							PAIR_TRIPLETS[i].PARAM_INDICIES[set1] = MAX_SET;	
						}					

					}

/*	OLD WAY				
					for(int set=1; set<PAIR_TRIPLETS[i].EQUIV_INDICIES.size(); set++)
					{
						if(PAIR_TRIPLETS[i].EQUIV_INDICIES[set] != PAIR_TRIPLETS[i].EQUIV_INDICIES[set-1])
							PAIR_TRIPLETS[i].PARAM_INDICIES[set] = PAIR_TRIPLETS[i].PARAM_INDICIES[set-1]+1;
						else
							PAIR_TRIPLETS[i].PARAM_INDICIES[set] = PAIR_TRIPLETS[i].PARAM_INDICIES[set-1];
					}
*/				
					PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS = PAIR_TRIPLETS[i].PARAM_INDICIES[PAIR_TRIPLETS[i].PARAM_INDICIES.size()-1]+1;
					PAIR_TRIPLETS[i].N_ALLOWED_POWERS = PAIR_TRIPLETS[i].PARAM_INDICIES.size();
				}

				// Set up triplet maps... Account for cases where triplet type is meaningless by setting mapped index to -1
					
				bool REAL_TRIPLET = false;
			
				string TEMP_STR_A, TEMP_STR_B, TEMP_STR_C;
			
				for(int i=0; i<NPAIR; i++)
				{
					for(int j=0; j<NPAIR; j++)
					{
						for(int k=0; k<NPAIR; k++)
						{	
							TEMP_STR_A = ATOM_PAIRS[i].ATM1TYP;
							TEMP_STR_A.append(ATOM_PAIRS[i].ATM2TYP);
							
							TEMP_STR_B = ATOM_PAIRS[j].ATM1TYP;
							TEMP_STR_B.append(ATOM_PAIRS[j].ATM2TYP);	
						
							TEMP_STR_C = ATOM_PAIRS[k].ATM1TYP;
							TEMP_STR_C.append(ATOM_PAIRS[k].ATM2TYP);	
						
							REAL_TRIPLET = false;	
						
							for(int m=0; m<NTRIP; m++)
							{
								if ((TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR3) ||
									(TEMP_STR_C == PAIR_TRIPLETS[m].ATMPAIR1) && (TEMP_STR_B == PAIR_TRIPLETS[m].ATMPAIR2) && (TEMP_STR_A == PAIR_TRIPLETS[m].ATMPAIR3) )
									
								{
									TEMP_STR = TEMP_STR_A;
									TEMP_STR.append(TEMP_STR_B);	
									TEMP_STR.append(TEMP_STR_C);

									TRIAD_MAP        .insert(make_pair(TEMP_STR,m));	
									TRIAD_MAP_REVERSE.insert(make_pair(m,TEMP_STR));
									REAL_TRIPLET = true;									
								
								}
								if(!REAL_TRIPLET)	// Then this is not an allowed triplet type!
									TRIAD_MAP        .insert(make_pair(TEMP_STR,-1));	
							}
						}
					}
				}
			}
					
			#if VERBOSITY == 1						
				cout << endl;
				cout << "	The following unique pair types have been identified:" << endl;
				for(int i=0;i<NPAIR; i++)
						cout << "		" << ATOM_PAIRS[i].PAIRIDX << "  " << ATOM_PAIRS[i].ATM1TYP << " " << ATOM_PAIRS[i].ATM2TYP << endl;
				#endif
					
			if(if_3b_cheby)
			{
				#if VERBOSITY == 1	
					cout << "	The following unique triplets of pair types and thier allowed pair polynomial powers have been identified:" << endl;
					for(int i=0;i<NTRIP; i++)
					{
						cout << "		" << PAIR_TRIPLETS[i].TRIPINDX << "  " << PAIR_TRIPLETS[i].ATMPAIR1 << " " << PAIR_TRIPLETS[i].ATMPAIR2 << " " << PAIR_TRIPLETS[i].ATMPAIR3 << ":";
						cout << " Number of unique sets of powers: " << PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS << " (" << PAIR_TRIPLETS[i].N_ALLOWED_POWERS << " total)..." << endl;	
						cout << "		     index  |  powers  |  equiv index  |  param index  " << endl;
						cout << "		   ----------------------------------------------------" << endl;					
					
						for(int j=0; j<PAIR_TRIPLETS[i].ALLOWED_POWERS.size(); j++)
						{
							//						cout << "		   " << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].X << " " << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Y << " " << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Z << endl;		
						 
										

							cout << "		      " << setw(6) << fixed << left << j << " ";
							cout << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].X  << " ";
							cout << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Y  << " ";
							cout << " " << setw(2) << fixed << left << PAIR_TRIPLETS[i].ALLOWED_POWERS[j].Z  << " ";
							cout << "       " << setw(8) << PAIR_TRIPLETS[i].EQUIV_INDICIES[j] << " ";
							cout << "       " << setw(8) << PAIR_TRIPLETS[i].PARAM_INDICIES[j] << endl; 
						
						}

					}
				#endif	
			}			
						
		}
		
		
		else if(LINE.find("# PAIRIDX #")!= string::npos) // Read the topology part. For now, ignoring index and atom types..
		{
			for(int i=0; i<NPAIR; i++)
			{
				cin >> LINE >> TEMP_PAIR.ATM1TYP >> TEMP_PAIR.ATM2TYP;
				
				TEMP_STR = TEMP_PAIR.ATM1TYP;
				TEMP_STR.append(TEMP_PAIR.ATM2TYP);
				TEMP_INT = PAIR_MAP[TEMP_STR];
				
				ATOM_PAIRS[TEMP_INT].PAIRIDX = TEMP_INT;

				cin >> LINE;
				ATOM_PAIRS[TEMP_INT].S_MINIM = double(atof(LINE.data()));
					
				cin >> LINE;
				ATOM_PAIRS[TEMP_INT].S_MAXIM = double(atof(LINE.data()));
					
				cin >> LINE; cin.ignore();
				ATOM_PAIRS[TEMP_INT].S_DELTA = double(atof(LINE.data()));
				
				cin >> LINE; cin.ignore();
				ATOM_PAIRS[TEMP_INT].LAMBDA = double(atof(LINE.data()));	
				
				ATOM_PAIRS[TEMP_INT].MIN_FOUND_DIST = 	1.0e10;	// Set an initial minimum distance		
				
				cin >> LINE; cin.ignore();
				
				if (LINE=="true"  || LINE=="True"  || LINE=="TRUE"  || LINE == "T" || LINE == "t")
				{
					cin >> ATOM_PAIRS[TEMP_INT].OVER_TO_ATM;
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[0];
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[1];
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[2];
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[3];
					cin >> ATOM_PAIRS[TEMP_INT].OVRPRMS[4];
					
					ATOM_PAIRS[TEMP_INT].USE_OVRPRMS = true;
				
				}
				else
				{
						ATOM_PAIRS[TEMP_INT].USE_OVRPRMS = false;
				}									
		
			}
			
			#if VERBOSITY == 1			
				cout << "	# PAIRIDX #     ";
				cout << "# ATM_TY1 #     ";
				cout << "# ATM_TY1 #     ";
				cout << "# S_MINIM #     ";
				cout << "# S_MAXIM #     ";
				cout << "# S_DELTA #     ";
				cout << "# MORSE_LAMBDA #";	
				cout << " # USEOVRP #     " << endl;
			#endif
				
			bool PRINT_OVR = false;
			
			for(int i=0; i<NPAIR; i++)
			{
				#if VERBOSITY == 1	
					cout << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
						 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
						 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
						 << setw(16) << left << ATOM_PAIRS[i].S_MINIM
						 << setw(16) << left << ATOM_PAIRS[i].S_MAXIM							 
						 << setw(16) << left << ATOM_PAIRS[i].S_DELTA
						 << setw(16) << left << ATOM_PAIRS[i].LAMBDA << " " 
						 << setw(16) << left << ATOM_PAIRS[i].USE_OVRPRMS << endl;
				#endif
					 
					 if(ATOM_PAIRS[i].USE_OVRPRMS)
						 PRINT_OVR = true;
			}
			
			if(PRINT_OVR && !fit_pover)
				if_subtract_coord = true;

			if(PRINT_OVR)
			{
				#if VERBOSITY == 1
					cout << "	# PAIRIDX #     ";
					cout << "# ATM_TY1 #     ";
					cout << "# ATM_TY1 #     ";				
					cout << "# P_OVERB #     ";
					cout << "# R_0_VAL #     ";
					cout << "# P_1_VAL #     ";
					cout << "# P_2_VAL #     ";
					cout << "# LAMBDA6 #" << endl;						
				#endif

				for(int i=0; i<NPAIR; i++)
				{
					#if VERBOSITY == 1					
						cout << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
							 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
							 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[0] 
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[1]
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[2] 
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[3]
							 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[4] << endl;	
					#endif													
				}											
			}
		}		
	
	}	
}
