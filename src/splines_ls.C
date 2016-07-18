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

static void read_lsq_input(string & INFILE, int & nframes, int & nlayers, bool & fit_coul, bool & coul_consv, bool & if_subtract_coord, bool & if_subtract_coul, bool & fit_pover, int & cheby_order, string & cheby_type, int & cheby_3b_order, int & invr_parms, int &  NATMTYP, bool & if_3b_cheby, vector<PAIRS> & ATOM_PAIRS, bool & WRAPCOORDS, map<string,int> & PAIR_MAP);


int main()
{
	//////////////////////////////////////////////////
	//
	// Job parameters: Declare and set defaults
	//
	//////////////////////////////////////////////////
	
	
	// BELOW: VARIABLES INTRODUCED BY BECKY TO ALLOW EXTENSION TO MULTIPLE ATOM TYPES
	// 
	vector<PAIRS> ATOM_PAIRS;		// Will store relevant info regarding atom interaction pair types.. THIS IS FOR LOOPS OF TYPE
	vector<FRAME> TRAJECTORY;		// Stores the trajectory information... box dimensions, atom types, coordinates, and forces
	vector<int>   ATOM_PAIR_TYPES;	// Fore use in double loop over atom pairs. Index corresponds to the overall double loop count. 
									// Provides an index that rells you the atom pair's ATOM_PAIRS's type index.. THIS IS FOR LOOPS OF TYPE
									// 	for(int a1=0;a1<nat-1;a1++)	for(int a2=a1+1;a2<nat;a2++)
	vector<int>   ATOM_PAIR_TYPES_ALL;	// THIS IS FOR LOOPS OF TYPE for(int ai=0; ai<nat; ai++), for(int ak=0; ak<nat; ak++)
	
	map<string,int> PAIR_MAP;
	map<string,int> TRIAD_MAP;		// We'll do this once we get to the 3b stuff...
	
	 
	
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
	bool fit_pover = false;			// If true, fit overcoordination parameters.
	bool coul_consv = false;		// If true, constraints will be applied to charge fitting to try to maintain consistency

	int cheby_order = 8; 			// Order of Chebyshev polynomial if used... set to 8 for DFTB Erep polynomial
	string cheby_type;				// How will distance be transformed?
	int cheby_3b_order = 2;
	int invr_parms = 4;				// currently uses 19 parameters per pair type
 	int NATMTYP = 0;

	bool if_3b_cheby = false; // ************* WHY ISN'T THIS ITS OWN PAIR_TYPE????!?!?!?!?!?!?!

	//////////////////////////////////////////////////
	//
	// Read and print input to screen
	//
	//////////////////////////////////////////////////
	
	cout << endl << "Reading input file..." << endl;

	read_lsq_input(INFILE, nframes, nlayers, fit_coul, coul_consv, ifsubtract_coord, ifsubtract_coul, fit_pover, cheby_order, cheby_type, cheby_3b_order, invr_parms, NATMTYP, if_3b_cheby, ATOM_PAIRS, WRAPCOORDS,PAIR_MAP);
		
	cout << "...input file read successful: " << endl << endl;
	
	
		
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
		
	vector<vector <vector <XYZ > > > A_MATRIX;		// [ # frames ] [ # atoms ]      [ #fitting parameters ]  .. replaces Als
	vector<vector <vector <XYZ > > > COULOMB_FORCES;	// [ # frames ] [ # pair types ] [ # atoms ]              .. replaces Coul_xx
	vector<vector <XYZ > >        P_OVER_FORCES;	// [ # frames ] [ # atoms ]                               .. replaces Pover
	
	A_MATRIX      .resize(nframes);
	COULOMB_FORCES.resize(nframes);
	P_OVER_FORCES.resize(nframes);

	// Figure out necessary dimensions for the force/force derivative vectors
	
	int tot_snum = 0; 

	for (int i=0; i<ATOM_PAIRS.size(); i++)
	{
		if ( (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) || (ATOM_PAIRS[i].PAIRTYP == "DFTBPOLY") )	// WHY DOES DFTBPOLY USE CHEBY? WHY DOES IT USE --3B-- CHEBY??
        {
			ATOM_PAIRS[i].SNUM          = cheby_order;
//			ATOM_PAIRS[i].SNUM_3B_CHEBY = cheby_3b_order;
        }
		
		else if (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV" ) 	// Set the distance transformation type
			ATOM_PAIRS[i].CHEBY_TYPE          = cheby_type;
			
		else if (ATOM_PAIRS[i].PAIRTYP == "INVRSE_R") 
			ATOM_PAIRS[i].SNUM = invr_parms;

		else 
			ATOM_PAIRS[i].SNUM = (2+floor((ATOM_PAIRS[i].S_MAXIM - ATOM_PAIRS[i].S_MINIM)/ATOM_PAIRS[i].S_DELTA))*2; //2 is for p0/m0/p1/m1.. WHAT??
	
		tot_snum += ATOM_PAIRS[i].SNUM;
	}
		
	cout << "The number of two-body non-coulomb parameters is " << tot_snum <<  endl;


/*	 I don't know what to do with this yet... 
	if ( if_3b_cheby ) 
	{
		num_cheby_3b = count_cheby_3b_params(snum_3b_cheby);
		cout << "Number of 3-Body Chebyshev parameter =" << num_cheby_3b << endl;
	}

	int tot_short_range = tot_snum + num_cheby_3b;	
*/		
	
	int tot_short_range = tot_snum; //	int tot_short_range = tot_snum + num_cheby_3b;	
	
	//////////////////////////////////////////////////
	//
	// Read and store MD trajectory (coords and forces)
	//
	//////////////////////////////////////////////////
		
	// Read in the trajectory
					 
	ifstream TRAJ_INPUT;
	TRAJ_INPUT.open(INFILE.data());
	
	TRAJECTORY.resize(nframes);
	
	cout << endl << "Reading in the trajectory file..." << endl;
	
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
/* DON'T NEED THIS ANYMORE SINCE MAPS WERE INTRODUCED!!!		
		// Set up the data object that will define the atom 
		// pair type for all atom pairs.. Main input fie HAS 
		// been read in at this point
		//
		// IMPORTANT: LOOPING METHOD MATTERS! EXPECTS:
		//
		// for(int a1=0;a1<nat-1;a1++)		// Double sum over atom pairs
		// 		for(int a2=a1+1;a2<nat;a2++)		
		
		if(i==0)
		{			
			for(int j=0; j<TRAJECTORY[i].ATOMS-1; j++)
			{
				for (int k=j+1; k<TRAJECTORY[i].ATOMS; k++)
				{
					for(int m=0; m<ATOM_PAIRS.size(); m++)
					{						
						if (TRAJECTORY[i].ATOMTYPE[j] == ATOM_PAIRS[m].ATM1TYP && TRAJECTORY[i].ATOMTYPE[k] == ATOM_PAIRS[m].ATM2TYP)
						{
							ATOM_PAIR_TYPES.push_back(m);
							TRAJECTORY[i].CHARGES[i] = ATOM_PAIRS[m].ATM1CHG;
							TRAJECTORY[i].CHARGES[k] = ATOM_PAIRS[m].ATM1CHG;
						}

						else if (TRAJECTORY[i].ATOMTYPE[k] == ATOM_PAIRS[m].ATM1TYP && TRAJECTORY[i].ATOMTYPE[j] == ATOM_PAIRS[m].ATM2TYP)
						{
							ATOM_PAIR_TYPES.push_back(m);
							TRAJECTORY[i].CHARGES[k] = ATOM_PAIRS[m].ATM1CHG;
							TRAJECTORY[i].CHARGES[j] = ATOM_PAIRS[m].ATM1CHG;
						}
					}
				}
			}
			bool MET;
			
			for(int j=0; j<TRAJECTORY[i].ATOMS; j++)
			{
				for(int k=0; k<TRAJECTORY[i].ATOMS; k++)
				{
					MET = false;
					
					for(int m=0; m<ATOM_PAIRS.size(); m++)
					{					
						if ((TRAJECTORY[i].ATOMTYPE[j] == ATOM_PAIRS[m].ATM1TYP && TRAJECTORY[i].ATOMTYPE[k] == ATOM_PAIRS[m].ATM2TYP) ||
							(TRAJECTORY[i].ATOMTYPE[k] == ATOM_PAIRS[m].ATM1TYP && TRAJECTORY[i].ATOMTYPE[j] == ATOM_PAIRS[m].ATM2TYP))
						{
							ATOM_PAIR_TYPES_ALL.push_back(m);
							MET = true;
							break;				
						}
					}
					
					if (!MET)
					{
							cout << "WARNING: No pair type found for atoms indexed " << j << " and " << k << endl;
							cout << "         Atom types " << TRAJECTORY[i].ATOMTYPE[j] << " and " <<  TRAJECTORY[i].ATOMTYPE[k] << endl;
							cout << "         Possible pairs are:" << endl;
							for(int m=0; m<ATOM_PAIRS.size(); m++)
								cout << "         " << ATOM_PAIRS[m].ATM1TYP << " " << ATOM_PAIRS[m].ATM2TYP << endl;
					}			
				}
			}		
		}		
*/		
		
	}
	
	cout << "...trajectory file read successful: " << endl << endl;

	//////////////////////////////////////////////////
	//
	// Setup A matrix, b vector, for least-squares:
	//
	//////////////////////////////////////////////////

	// Setup the A and Coulomb "matricies"
	
	cout << "Setting up the matricies for A, Coulomb forces, and overbonding..." << endl;
	
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

	cout << "...matrix setup complete: " << endl << endl;


	//////////////////////////////////////////////////
	//
	// Generate A matrix, b vector, for least-squares:
	//
	//////////////////////////////////////////////////
	
	
	// cout.precision(10);		// WE SHOULD AT LEAST MOVE THIS TO SOMEWHERE MORE REASONABLE.. LIKE THE SECTION WHERE THE OUTPUT IS ACTUALLY PRINTED
	cout.precision(16);
	
	cout << "...Populating the matricies for A, Coulomb forces, and overbonding..." << endl;
	
	for(int i=0; i<nframes; i++)
	{
		ZCalc_Deriv(ATOM_PAIRS,TRAJECTORY[i], ATOM_PAIR_TYPES, A_MATRIX[i], COULOMB_FORCES[i], nlayers, if_3b_cheby, PAIR_MAP);

		if ( ifsubtract_coord ) // Subtract over-coordination forces from force to be output.
			SubtractCoordForces(TRAJECTORY[i], false, P_OVER_FORCES[i],  ATOM_PAIRS, ATOM_PAIR_TYPES_ALL, PAIR_MAP);	
		
		if (ifsubtract_coul) 
			ZCalc_Ewald(TRAJECTORY[i], ATOM_PAIR_TYPES, ATOM_PAIRS, PAIR_MAP);

		if ( fit_pover )	// Fit the overcoordination parameter.
			SubtractCoordForces(TRAJECTORY[i], true, P_OVER_FORCES[i],  ATOM_PAIRS, ATOM_PAIR_TYPES_ALL, PAIR_MAP);	
		
			
	}	
	
	cout << "...matrix population complete: "  << endl << endl;
	
	cout << "Printing matricies..." << endl;
	
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
		  
		}  	
    }
  
	// I DON'T UNDERSTAND THIS PART... WHAT IS IS DOING? WHY IS IT PRINTED AT THE END OF THE A FILE?
	// WHY IS IT REPEATED INVERSELY FOR A?
	//
	// A and B file charge conservation (q00 + 4 qOH + qHH = 0 )
	  
	if ( fit_coul && coul_consv) 
	{
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
	} 	
	   
	fileA.close();
	fileb.close();
	  	  
	//////////////////////////////////////////////////
	//
	// Print out the params file header
	//
	//////////////////////////////////////////////////	   
	  
	ofstream header;
	header.open("params.header");
  
	int NPAIR =  ATOM_PAIRS.size();
  
	header << "npair " << NPAIR << endl;
	
	for (int i=0; i<NPAIR; i++) 
	{
		header << fixed << setw(8) << setprecision(5) << ATOM_PAIRS[i].S_MINIM << fixed << setw(8) << setprecision(5) << ATOM_PAIRS[i].S_MAXIM;
		  
		if (ATOM_PAIRS[i].PAIRTYP == "INVRSE_R") 
			header << " " << ATOM_PAIRS[i].SNUM << endl;
		else if (( ATOM_PAIRS[i].PAIRTYP != "CHEBYSHEV" ) && ( ATOM_PAIRS[i].PAIRTYP != "INVRSE_R")) 
		{
			header << fixed << setw(8) << setprecision(5) << ATOM_PAIRS[i].S_DELTA;
			header << " " << ATOM_PAIRS[i].SNUM << endl;
		}
		else if (ATOM_PAIRS[i].PAIRTYP == "CHEBYSHEV") 
		{
			header << fixed << setw(8) << setprecision(5) << ATOM_PAIRS[i].LAMBDA;
			 
			if ( if_3b_cheby ) 
				header << " " << ATOM_PAIRS[i].SNUM << " " << ATOM_PAIRS[i].SNUM_3B_CHEBY << endl;
			else 
				header << " " << ATOM_PAIRS[i].SNUM << endl;	 	
		}
	}

	if ( fit_pover || ifsubtract_coord ) 
	{
		header << 5 << endl;
		  
				for (int i=0; i<NPAIR; i++) 
					if(ATOM_PAIRS[i].USE_OVRPRMS)
						for ( int k=0; k<5; k++ )
							header << scientific << setw(21) << setprecision(13) <<  ATOM_PAIRS[i].OVRPRMS[k] << endl;
		

				for ( int k=0; k<5; k++ )
					header << scientific << setw(21) << setprecision(13) <<  ATOM_PAIRS[0].OVRPRMS[k] << endl;
		
	}
	  
	header.close();

	cout << "	Minimum distances between atoms: (Angstr.)" << endl;
	  
	for (int k=0; k<NPAIR; k++)
		cout << "		" << k << "	" << ATOM_PAIRS[k].ATM1TYP << " " << ATOM_PAIRS[k].ATM2TYP << "	" << fixed << setprecision(3) << ATOM_PAIRS[k].MIN_FOUND_DIST << endl;
	  
	  
	cout << "...matrix printing complete: " << endl << endl;
	
	  
return 0;		  
}




	//////////////////////////////////////////////////
	//
    // Function definitions
	//
	//////////////////////////////////////////////////



// Read program input from the file "splines_ls.in".
static void read_lsq_input(string & INFILE, int & nframes, int & nlayers, bool & fit_coul, bool & coul_consv, bool & if_subtract_coord, bool & if_subtract_coul, bool & fit_pover, int & cheby_order, string & cheby_type, int & cheby_3b_order, int & invr_parms, int &  NATMTYP, bool & if_3b_cheby, vector<PAIRS> & ATOM_PAIRS, bool & WRAPCOORDS, map<string,int> & PAIR_MAP)		   
{
	bool   FOUND_END = false;
	string LINE;
	string TEMP_STR;
	PAIRS  TEMP_PAIR;
	int    TEMP_INT;
	string TEMP_TYPE;
	int    NPAIR;
	
	
	while (FOUND_END == false)
	{
		getline(cin,LINE);


		if(LINE.find("# ENDFILE #") != string::npos)
		{
			FOUND_END = true;
			break;
		}
		
		else if(LINE.find("# WRAPTRJ #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			WRAPCOORDS = bool(LINE.data());
			cout << "	# WRAPTRJ #: " << WRAPCOORDS << endl;	
		}
		
		else if(LINE.find("# TRJFILE #") != string::npos)
		{
			cin >> INFILE; cin.ignore();
			cout << "	# TRJFILE #: " << INFILE << endl;	
		}			
		
		else if(LINE.find("# NFRAMES #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			nframes = int(atof(LINE.data()));
			cout << "	# NFRAMES #: " << nframes << endl;			
		}
		
		else if(LINE.find("# NLAYERS #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			nlayers = int(atof(LINE.data()));
			cout << "	# NLAYERS #: " << nlayers << endl;
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
			
			cout << "	# FITCOUL #: ";
			
			if (fit_coul)
				cout << "true" << endl;				
			else
				cout << "false" << endl;
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
			
			cout << "	# CNSCOUL #: ";
			
			if (coul_consv)
				cout << "true" << endl;				
			else
				cout << "false" << endl;							
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
			
			cout << "	# FITPOVR #: ";
			
			if (fit_pover)
				cout << "true" << endl;				
			else
				cout << "false" << endl;
		}
		
		else if(LINE.find("# PAIRTYP #") != string::npos)
		{
			cin >> LINE; cin.ignore();
			
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
			
			cout << "	# PAIRTYP #: " << LINE;
			
			if(LINE != "DFTBPOLY")
				cout << " ....NOTE: Forces reported in units of kcal/mol." << endl;
			else
				cout << " ....NOTE: Forces reported in atomic units." << endl;

			
		}
		
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
			
			cout << "	# SUBCRDS #: ";
			
			if (if_subtract_coord)
				cout << "true" << endl;				
			else
				cout << "false" << endl;

		}
			
		else if( ((TEMP_TYPE == "CHEBYSHEV")||(TEMP_TYPE == "DFTBPOLY")) && (LINE.find("# POLORDR #") != string::npos))
		{
			cin >> LINE; cin.ignore();

			cheby_order = int(atof(LINE.data()));
			cout << "	# POLORDR #: " << cheby_order << endl;	
		}
		
		else if( (TEMP_TYPE == "CHEBYSHEV") && (LINE.find("# CHBTYPE #") != string::npos))
		{
			cin >> LINE; cin.ignore();
			cheby_type = LINE;

			cout << "	# CHBTYPE #: " << cheby_type << endl;	
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
			
			cout << "	# NATMTYP #: " << NATMTYP << endl;	
			NPAIR = NATMTYP*(NATMTYP+1)/2;
			ATOM_PAIRS.resize(NPAIR);
			
		}

		else if(LINE.find("# TYPEIDX #")!= string::npos)
		{
			cout << "	# TYPEIDX #    # ATM_TYP #    # ATMCHRG #    # ATMMASS #" << endl;
			
			// Figure out the number of non-unique pairs

			TEMP_INT = NATMTYP*NATMTYP;
			
			for(int i=0; i<NATMTYP; i++)
			{
				
				// Set the first atom pair types to be of type OO, HH, CC, etc...
				
				ATOM_PAIRS[i].PAIRTYP    = TEMP_TYPE;
				ATOM_PAIRS[i].CHEBY_TYPE = cheby_type;
				
				cin >> LINE >> ATOM_PAIRS[i].ATM1TYP >> LINE;
				ATOM_PAIRS[i].ATM1CHG = double(atof(LINE.data()));

				cin >> LINE; cin.ignore();
				ATOM_PAIRS[i].ATM1MAS = double(atof(LINE.data()));
				
				ATOM_PAIRS[i].ATM2TYP = ATOM_PAIRS[i].ATM1TYP;
				ATOM_PAIRS[i].ATM2CHG = ATOM_PAIRS[i].ATM1CHG;
				ATOM_PAIRS[i].ATM2MAS = ATOM_PAIRS[i].ATM1MAS;
				
				cout << " 	" << setw(15) << left << i+1 
					 << setw(15) << left << ATOM_PAIRS[i].ATM1TYP 
				     << setw(15) << left << ATOM_PAIRS[i].ATM1CHG
				     << setw(15) << left << ATOM_PAIRS[i].ATM1MAS << endl;
			}
			
			// Set up all possible unique pair types
			
			TEMP_INT = NATMTYP;
			
			for(int i=0; i<NATMTYP; i++)
			{
				for(int j=i+1; j<NATMTYP; j++)
				{	
					ATOM_PAIRS[TEMP_INT].PAIRTYP = ATOM_PAIRS[i].PAIRTYP;
					
					ATOM_PAIRS[TEMP_INT].CHEBY_TYPE = ATOM_PAIRS[i].CHEBY_TYPE;
														
					ATOM_PAIRS[TEMP_INT].ATM1TYP = ATOM_PAIRS[i].ATM1TYP;
					ATOM_PAIRS[TEMP_INT].ATM2TYP = ATOM_PAIRS[j].ATM1TYP;
					
					ATOM_PAIRS[TEMP_INT].ATM1CHG = ATOM_PAIRS[i].ATM1CHG;
					ATOM_PAIRS[TEMP_INT].ATM2CHG = ATOM_PAIRS[j].ATM1CHG;
					
					ATOM_PAIRS[TEMP_INT].ATM1MAS = ATOM_PAIRS[i].ATM1MAS;
					ATOM_PAIRS[TEMP_INT].ATM2MAS = ATOM_PAIRS[j].ATM1MAS;												
					
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
						}
							 
						
					}
				}
			}
/*			
			// Sanity check
			
			cout << "SANITY CHECK " << endl;			
			for(int i=0; i<NATMTYP; i++)
			{
				for(int j=0; j<NATMTYP; j++)
				{
					TEMP_STR = ATOM_PAIRS[i].ATM1TYP;
					TEMP_STR.append(ATOM_PAIRS[j].ATM2TYP);										
					cout << ATOM_PAIRS[i].ATM1TYP << " " << ATOM_PAIRS[j].ATM1TYP << " " << TEMP_STR << " " <<  PAIR_MAP[TEMP_STR] << endl;
				}
			}											
*/											
			
						
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
			
			cout << "	# PAIRIDX #     ";
			cout << "# ATM_TY1 #     ";
			cout << "# ATM_TY1 #     ";
			cout << "# S_MINIM #     ";
			cout << "# S_MAXIM #     ";
			cout << "# S_DELTA #     ";
			cout << "# MORSE_LAMBDA #";	
			cout << " # USEOVRP #     " << endl;
			
			bool PRINT_OVR = false;
			
			for(int i=0; i<NPAIR; i++)
			{
				cout << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
					 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
					 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
					 << setw(16) << left << ATOM_PAIRS[i].S_MINIM
					 << setw(16) << left << ATOM_PAIRS[i].S_MAXIM							 
					 << setw(16) << left << ATOM_PAIRS[i].S_DELTA
					 << setw(16) << left << ATOM_PAIRS[i].LAMBDA << " " 
					 << setw(16) << left << ATOM_PAIRS[i].USE_OVRPRMS << endl;
					 
					 if(ATOM_PAIRS[i].USE_OVRPRMS)
						 PRINT_OVR = true;
			}

			if(PRINT_OVR)
			{
				cout << "	# PAIRIDX #     ";
				cout << "# ATM_TY1 #     ";
				cout << "# ATM_TY1 #     ";				
				cout << "# P_OVERB #     ";
				cout << "# R_0_VAL #     ";
				cout << "# P_1_VAL #     ";
				cout << "# P_2_VAL #     ";
				cout << "# LAMBDA6 #" << endl;						

				for(int i=0; i<NPAIR; i++)
				{
					cout << "	" << setw(16) << left << ATOM_PAIRS[i].PAIRIDX 
						 << setw(16) << left << ATOM_PAIRS[i].ATM1TYP
						 << setw(16) << left << ATOM_PAIRS[i].ATM2TYP 
						 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[0] 
						 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[1]
						 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[2] 
						 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[3]
						 << setw(16) << left << ATOM_PAIRS[i].OVRPRMS[4] << endl;														
				}											
			}
		}		
	}	
}
