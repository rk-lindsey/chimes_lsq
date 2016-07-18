#include<iomanip>

#include "functions.h"


using namespace std;

//////////////////////////////////////////
//
//	FUNCTION HEADERS
//
//////////////////////////////////////////

static void ZCalc_Spline_Deriv(FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);

static void ZCalc_Cheby_Deriv (FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);

static void ZCalc_Poly_Deriv  (FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);	

static void ZCalc_InvR_Deriv  (FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);

////////////////////////////////////////////////////////////
//
// FUNCTIONS THAT SET UP A TO HOLD/COMPUTE DERIVATIVES
//
////////////////////////////////////////////////////////////

// FUNCTION UPDATED
void ZCalc_Deriv(vector<PAIRS> & ATOM_PAIRS, FRAME & FRAME_TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<vector <XYZ > > & FRAME_A_MATRIX, vector<vector <XYZ > > & FRAME_COULOMB_FORCES, const int nlayers, 	bool if_3b_cheby, map<string,int> & PAIR_MAP)
		
	// NEW VARIABLES TO ADD (to be passed in): snum_tot 
	// IS P_OVER_FORCES EVEN USED HERE??? .. I don't think so.. removing from header
	// vector<XYZ <XYZ > > & FRAME_P_OVER_FORCES
	// IS SR_PAIR_T EVEN USED HERE? .. I don't think so.. removing from header
	// 	Sr_pair_t pair_type,
{
  	// Will ewald calculations be needed? DFTB doesn't use the "default"
    // Ewald calculation because it has its own special way of dealing with charges
		
	bool if_ewald ;
  
    if (ATOM_PAIRS[0].PAIRTYP != "DFTBPOLY") 
  	  if_ewald = true;
    else
  	  if_ewald = false;

	if ( if_ewald ) 					// FUNCTION UPDATED
		ZCalc_Ewald_Deriv(FRAME_TRAJECTORY, ATOM_PAIR_TYPES, ATOM_PAIRS, FRAME_COULOMB_FORCES, PAIR_MAP);	
	
    if ( ATOM_PAIRS[0].PAIRTYP == "SPLINE" )			// FUNCTION UPDATED
		ZCalc_Spline_Deriv(FRAME_TRAJECTORY, ATOM_PAIR_TYPES, ATOM_PAIRS, FRAME_A_MATRIX, nlayers, PAIR_MAP);
	
	else if ( ATOM_PAIRS[0].PAIRTYP == "CHEBYSHEV" )	// FUNCTION UPDATED
		ZCalc_Cheby_Deriv(FRAME_TRAJECTORY, ATOM_PAIR_TYPES, ATOM_PAIRS, FRAME_A_MATRIX, nlayers, PAIR_MAP);

    else if ( ATOM_PAIRS[0].PAIRTYP == "DFTBPOLY" )	// FUNCTION UPDATED
		ZCalc_Poly_Deriv(FRAME_TRAJECTORY, ATOM_PAIR_TYPES, ATOM_PAIRS, FRAME_A_MATRIX, nlayers, PAIR_MAP);

    else if ( ATOM_PAIRS[0].PAIRTYP == "INVRSE_R" )	// FUNCTION UPDATED
		ZCalc_InvR_Deriv(FRAME_TRAJECTORY, ATOM_PAIR_TYPES, ATOM_PAIRS, FRAME_A_MATRIX, nlayers, PAIR_MAP);
		
	
		
    else 
    {
		cout << "Error: bad pairtype in ZCalc_Deriv\n" ;
		exit(1) ;
    }
  
/* Commenting out for now, because this will all be updated when Larry's edits are encorperated...  
  
    if (if_3b_cheby)
		ZCalc_3B_Cheby_Deriv(FRAME_TRAJECTORY, ATOM_PAIRS, FRAME_A_MATRIX, nlayers);
//		ZCalc_3B_Cheby_Deriv(FRAME & TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, vector<XYZ <XYZ > > & FRAME_A_MATRIX, const int nlayers);
*/	
 	
	
		
}	


// FUNCTION UPDATED
static void ZCalc_Spline_Deriv(FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)		   
// Original comment: Calculate derivatives of the forces wrt the spline parameters. Stores minimum distance between a pair of atoms in mind.
// New Note: This doesn't actually calcuate any derivatives.. it is just populating A with the cubic hermite basis polynomials needed for fitting
{
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double t;
	double h00,h10,h01,h11;
	int    k0;
	double x,x0;
	int    vstart,kstart;
	string TEMP_STR;

	// Main loop for SPLINE terms:
	
	int curr_pair_type_idx;

	for(int a1=0;a1<TRAJECTORY.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<TRAJECTORY.ATOMS;a2++)
		{
			TEMP_STR = TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(TRAJECTORY.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
			//calculate vstart:

			vstart = curr_pair_type_idx * ATOM_PAIRS[curr_pair_type_idx].SNUM;
			 
			// Start with minimum image convention.  Use layers to access larger distances if desired.	
		   
			RVEC.X  = TRAJECTORY.COORDS[a2].X - TRAJECTORY.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/TRAJECTORY.BOXDIM.X )  * TRAJECTORY.BOXDIM.X;
		   
			RVEC.Y = TRAJECTORY.COORDS[a2].Y - TRAJECTORY.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/TRAJECTORY.BOXDIM.Y )  * TRAJECTORY.BOXDIM.Y;
		   
			RVEC.Z = TRAJECTORY.COORDS[a2].Z - TRAJECTORY.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/TRAJECTORY.BOXDIM.Z )  * TRAJECTORY.BOXDIM.Z;

      
			for(int n1=-1*nlayers;n1<nlayers+1;n1++)
			{
				for(int n2=-1*nlayers;n2<nlayers+1;n2++)
				{
					for(int n3=-1*nlayers;n3<nlayers+1;n3++)
					{
						RAB.X = RVEC.X + n1 * TRAJECTORY.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * TRAJECTORY.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * TRAJECTORY.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );

						if ( rlen < ATOM_PAIRS[curr_pair_type_idx].MIN_FOUND_DIST ) 
							ATOM_PAIRS[curr_pair_type_idx].MIN_FOUND_DIST = rlen;
					   
					   
						if(rlen > ATOM_PAIRS[curr_pair_type_idx].S_MINIM and rlen < ATOM_PAIRS[curr_pair_type_idx].S_MAXIM) //spline term calculated w/cutoff:
						{
							
							// Setting up for Hermite cubic quadrature:
							
						 	// This part: figure out which distance bin ("unit") that rlen falls into. Use the value of that bin "t" (rather than rlen ("x") ) to 
							// populate the matrix table. Use the lower and upper bounds of the bin to deterime the two constraining data points for the interpolation
							// polynomial fitting
							 
							 
							k0 = int( floor( (rlen - ATOM_PAIRS[curr_pair_type_idx].S_MINIM) / ATOM_PAIRS[curr_pair_type_idx].S_DELTA ) ); // binned location of distance along allowed range

						   
							if ( k0 >= ATOM_PAIRS[curr_pair_type_idx].SNUM/2 - 1 ) // Keep from overrunning table.
							{
								cout << "Table overrun: rlen = " << rlen << endl ;
								continue ;
							}

 						   	x = rlen;
	    
 						   	x0 = ATOM_PAIRS[curr_pair_type_idx].S_MINIM + ATOM_PAIRS[curr_pair_type_idx].S_DELTA * k0;
 						   	t  = (x-x0) / ATOM_PAIRS[curr_pair_type_idx].S_DELTA;

						   	// for classical cubic quadrature, you would only change h00, h10, h01, h11 polynomials!
						   
						  	// Setup the basis functions for cubic hermite polynomial interpolation...
						   	//
						   	// See: https://www3.nd.edu/~coast/jjwteach/www/www/30125/pdfnotes/lecture5_9v14.pdf
						   	// Look for "Each basis function is a third degree polynomial" 
						   	//
						   	// The fitting polynomial, g(x) = f(x_0)*alpha_0 + f(x_1)*alpha_1 + f'(x_0)*beta_0 + f'(x_1)*beta_1
						   	// But we're using t instead of x, and t_0 and t_1 are the two termini of the distance bin.
						   	// i.e. if x = 5.54, delta = 0.1, fall into a bin of value 5.55 with endpoints 5.5 and 5.6.
					   
						   	h00  =  2 * t * t * t - 3 * t * t + 1;			// alpha_0.. basis function for first point
						   
						   	h10  =  t * t * t  -2 * t * t + t;				// beta_0..  basis function for the first derivative at first point
						   	h10 *=  ATOM_PAIRS[curr_pair_type_idx].S_DELTA;	// derivative terms have extra factor.
						   
						   	h01  = -2 * t * t * t + 3 * t * t;				// alpha_1.. basis function for second point
						   
						   	h11  =  t * t * t - t * t;						// beta_1..  basis function for first derivative at second point
						   	h11 *=  ATOM_PAIRS[curr_pair_type_idx].S_DELTA;	// derivative terms have extra factor.

						   	kstart = k0 * 2;								//for each distance unit you have 2 entries. ...I'm guessing these are the start and end of the distance unit bin
						    
						    // Populate the A matrix with these values, muliplied times the xyz unit vectors to get directionality
							//
							// Note:
							//
							// vstart: takes care of index w/r/t possible pair types (i.e. OO, OH, HH...)
							// kstart: takes care of index w/r/t location along allowable pair distances
							// {0..3}: takes care of whether value corresponds to h00, h10, etc.
							//
							// So for a given item in a "Block" of A might correspond to, for example, type H-O, distance 5.0 along range 0.6 - 6.0 \AA, and basis function h01.
							
						   	FRAME_A_MATRIX[a1][vstart+kstart+0].X += h00 * RAB.X / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+1].X += h10 * RAB.X / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+2].X += h01 * RAB.X / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+3].X += h11 * RAB.X / rlen;
						  
						   	FRAME_A_MATRIX[a1][vstart+kstart+0].Y += h00 * RAB.Y / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+1].Y += h10 * RAB.Y / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+2].Y += h01 * RAB.Y / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+3].Y += h11 * RAB.Y / rlen;
							  
						   	FRAME_A_MATRIX[a1][vstart+kstart+0].Z += h00 * RAB.Z / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+1].Z += h10 * RAB.Z / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+2].Z += h01 * RAB.Z / rlen;
						   	FRAME_A_MATRIX[a1][vstart+kstart+3].Z += h11 * RAB.Z / rlen;
							
						   	FRAME_A_MATRIX[a2][vstart+kstart+0].X -= h00 * RAB.X / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+1].X -= h10 * RAB.X / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+2].X -= h01 * RAB.X / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+3].X -= h11 * RAB.X / rlen;
						  
						   	FRAME_A_MATRIX[a2][vstart+kstart+0].Y -= h00 * RAB.Y / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+1].Y -= h10 * RAB.Y / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+2].Y -= h01 * RAB.Y / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+3].Y -= h11 * RAB.Y / rlen;
							  
						   	FRAME_A_MATRIX[a2][vstart+kstart+0].Z -= h00 * RAB.Z / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+1].Z -= h10 * RAB.Z / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+2].Z -= h01 * RAB.Z / rlen;
						   	FRAME_A_MATRIX[a2][vstart+kstart+3].Z -= h11 * RAB.Z / rlen;							
	    
						   }
					   }
				   }
			   }
		   }
	   }
	return;
}

// FUNCTION UPDATED
static void ZCalc_Cheby_Deriv(FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)	
// Calculate derivatives of the forces wrt the Chebyshev parameters. Stores minimum distance between a pair of atoms in mind.
// For some reason, derivatives are only calculated for morse-type cheby polynomials. Other forms don't actually use derivatives...
{
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double x;
	int vstart;
	static double *Tn, *Tnd;
	static bool called_before = false;
	
	double xavg, xdiff;
	double xmin; 
	double xmax; 
	double rlen3; 				  
	double fcut0; 
	double fcut; 
	double fcutderiv; 				
	double deriv;
	double tmp_doub; 	

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;
		
		for ( int i = 0; i < ATOM_PAIRS.size(); i++ ) 
			if (ATOM_PAIRS[i].SNUM > dim ) 
				dim = ATOM_PAIRS[i].SNUM;	 // Will probably get a complaint b/c not static...
		
		dim++;
		Tn   = new double [dim];
		Tnd  = new double [dim];
	}

	// Main loop for Chebyshev terms:
	
	string TEMP_STR;
	int curr_pair_type_idx;

	for(int a1=0;a1<TRAJECTORY.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<TRAJECTORY.ATOMS;a2++)
		{
			
			
			TEMP_STR = TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(TRAJECTORY.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
			//calculate vstart:
		
			vstart = curr_pair_type_idx * ATOM_PAIRS[curr_pair_type_idx].SNUM;

			// Start with minimum image convention.  Use layers to access larger distances if desired.	
		   
			RVEC.X  = TRAJECTORY.COORDS[a2].X - TRAJECTORY.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/TRAJECTORY.BOXDIM.X )  * TRAJECTORY.BOXDIM.X;
		   
			RVEC.Y = TRAJECTORY.COORDS[a2].Y - TRAJECTORY.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/TRAJECTORY.BOXDIM.Y )  * TRAJECTORY.BOXDIM.Y;
		   
			RVEC.Z = TRAJECTORY.COORDS[a2].Z - TRAJECTORY.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/TRAJECTORY.BOXDIM.Z )  * TRAJECTORY.BOXDIM.Z;
      
			for(int n1=-1*nlayers;n1<nlayers+1;n1++)
			{
				for(int n2=-1*nlayers;n2<nlayers+1;n2++)
				{
					for(int n3=-1*nlayers;n3<nlayers+1;n3++)
					{
						RAB.X = RVEC.X + n1 * TRAJECTORY.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * TRAJECTORY.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * TRAJECTORY.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );

						if ( rlen < ATOM_PAIRS[curr_pair_type_idx].MIN_FOUND_DIST )	// spline term calculated w/cutoff:
							ATOM_PAIRS[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

						if(rlen > ATOM_PAIRS[curr_pair_type_idx].S_MINIM and rlen < ATOM_PAIRS[curr_pair_type_idx].S_MAXIM)
						{
							// Chebyshev polynomials are only defined on the range [-1,1], so transorm the pair distance
							// in a way that allows it to fall along that range. Options are:
							//
							// x = 1/pair_dist				// Inverse r, 
							// x = exp(pair_dist/lambda)	// Morse-type
							// x = pair_dist				// default type
							// 
							// All types are normalized by s_min to s_max range to fall along [-1,1]
														
							
							if ( ATOM_PAIRS[curr_pair_type_idx].CHEBY_TYPE == "INVRSE_R" ) 
							{
								xavg  =  0.5 * (1.0/ATOM_PAIRS[curr_pair_type_idx].S_MINIM + 1.0/ATOM_PAIRS[curr_pair_type_idx].S_MAXIM); // midpoint of possible pair distances in r^-1 space
								xdiff =  0.5 * (1.0/ATOM_PAIRS[curr_pair_type_idx].S_MINIM - 1.0/ATOM_PAIRS[curr_pair_type_idx].S_MAXIM); // width of possible pair distances in r^-1 space
								x     = (1.0/rlen-xavg) / xdiff;																		  // pair distances in r^-1 space, normalized to fit over [-1,1]
							} 
							else if ( ATOM_PAIRS[curr_pair_type_idx].CHEBY_TYPE == "MORSE" ) 
							{
								xmin  = exp(-ATOM_PAIRS[curr_pair_type_idx].S_MAXIM/ATOM_PAIRS[curr_pair_type_idx].LAMBDA); 
								xmax  = exp(-ATOM_PAIRS[curr_pair_type_idx].S_MINIM/ATOM_PAIRS[curr_pair_type_idx].LAMBDA); 
								xavg  = 0.5 * (xmin + xmax);																// midpoint of possible pair distances in morse space
								xdiff = 0.5 * (xmax - xmin);																// width of possible pair distances in morse space
								x = (exp(-rlen/ATOM_PAIRS[curr_pair_type_idx].LAMBDA)-xavg)/xdiff;							// pair distances in morse space, normalized to fit over [-1,1]
							}
							else if (ATOM_PAIRS[curr_pair_type_idx].CHEBY_TYPE == "DEFAULT")
							{
								xavg  = 0.5 * (ATOM_PAIRS[curr_pair_type_idx].S_MINIM + ATOM_PAIRS[curr_pair_type_idx].S_MAXIM); // midpoint of possible pair distances
								xdiff = 0.5 * (ATOM_PAIRS[curr_pair_type_idx].S_MAXIM - ATOM_PAIRS[curr_pair_type_idx].S_MINIM); // width of possible pair distances
								x = (rlen-xavg) / xdiff;																		 // pair distances, normalized to fit over [-1,1]
							}
							else
							{
								cout << "ERROR: Undefined CHBTYPE: " << ATOM_PAIRS[curr_pair_type_idx].CHEBY_TYPE << endl;
								cout << "       Excepted values are \"DEFAULT\", \"INVRSE_R\", or \"MORSE\". " << endl;
								exit(1);
							}

							if ( x < -1.0 ) // Why are we setting it as -1? shouldn't something else be done here instead? ...since we're outside the range of the fit??
							{
								cout << "Warning:  r < rmin" << endl;
								x = -1.0;
							}
							else if ( x > 1.0 ) // Added a warning...Do we need to turn off the interaction, because it's beyond rcut?
							{
								cout << "Warning:  r < rmax" << endl;								
								x = 1.0;								
							}

							// Generate Chebyshev polynomials by recursion. 
							// 
							// What we're doing here. Want to fit using Cheby polynomials of the 1st kind. "T_n(x)."
							// We need to calculate the derivative of these polynomials.
							// Derivatives are defined through use of Cheby polynomials of the 2nd kind "U_n(x)", as:
							//
							// d/dx[ T_n(x) = n * U_n-1(x)] 
							// 
							// So we need to first set up the 1st-kind polynomials ("Tn[]")
							// Then, to compute the derivatives ("Tnd[]"), first set equal to the 2nd-kind, then multiply by n to get the der's
						  
							// First two 1st-kind Chebys:
						  
							Tn[0] = 1.0;
							Tn[1] = x;
							
							// Start the derivative setup. Set the first two 1st-kind Cheby's equal to the first two of the 2nd-kind
						  
							Tnd[0] = 1.0;
							Tnd[1] = 2.0 * x;
							
							// Use recursion to set up the higher n-value Tn and Tnd's
						  
							for ( int i = 2; i <= ATOM_PAIRS[curr_pair_type_idx].SNUM; i++ ) 
							{
								Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2];
								Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2];
							}
							
							// Now multiply by n to convert Tnd's to actual derivatives of Tn

							for ( int i = ATOM_PAIRS[curr_pair_type_idx].SNUM; i >= 1; i-- ) 
								Tnd[i] = i * Tnd[i-1];

							Tnd[0] = 0.0;
							
							// Define/setup a penalty function to deal with the cases where the pair distance is close to rmin.. to discourage
							// the trajectory from heading towards though poorly sampled regions of the PES.

							rlen3 = rlen * rlen * rlen;
							
							if ( ATOM_PAIRS[curr_pair_type_idx].CHEBY_TYPE == "MORSE" ) // IS THE MORSE TYPE THE ONLY ONE THAT USES DERIVATIVES?? WHY??
							{
								// fcut and fcutderv are the form that the penalty func and its derivative for the morse-type pair distance transformation
								
								fcut0     = (1.0 - rlen/ATOM_PAIRS[curr_pair_type_idx].S_MAXIM);
								fcut      = fcut0 * fcut0 * fcut0;
								fcutderiv = -3.0 * fcut0 * fcut0 / ATOM_PAIRS[curr_pair_type_idx].S_MAXIM;
							
								for ( int i = 0; i < ATOM_PAIRS[curr_pair_type_idx].SNUM; i++ ) 
								{
									
									// NOTE: All these extra terms are coming from:
									//
									// 1. Chain rule to account for transformation from morse-type pair distance to x
									// 2. Product rule coming from pair distance dependence of fcut, the penalty function
									
									tmp_doub = (fcut * Tnd[i+1] *(-exp(-rlen/ATOM_PAIRS[curr_pair_type_idx].LAMBDA)/ATOM_PAIRS[curr_pair_type_idx].LAMBDA)/xdiff + fcutderiv * Tn[i+1] );
									
									// Finally, account for the x, y, and z unit vectors
									
									deriv = tmp_doub * RAB.X / rlen;
									FRAME_A_MATRIX[a1][vstart+i].X += deriv;
									FRAME_A_MATRIX[a2][vstart+i].X -= deriv;
		
									deriv = tmp_doub * RAB.Y / rlen; 
									FRAME_A_MATRIX[a1][vstart+i].Y += deriv;
									FRAME_A_MATRIX[a2][vstart+i].Y -= deriv;
																		
									deriv = tmp_doub * RAB.Z / rlen;
									FRAME_A_MATRIX[a1][vstart+i].Z += deriv;
									FRAME_A_MATRIX[a2][vstart+i].Z -= deriv;

								}
							}
							else 
							{
								for ( int i = 0; i < ATOM_PAIRS[curr_pair_type_idx].SNUM; i++ ) 
								{
									FRAME_A_MATRIX[a1][vstart+i].X += (ATOM_PAIRS[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.X/rlen3;
									FRAME_A_MATRIX[a1][vstart+i].Y += (ATOM_PAIRS[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.Y/rlen3;
									FRAME_A_MATRIX[a1][vstart+i].Z += (ATOM_PAIRS[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.Z/rlen3;
									
									FRAME_A_MATRIX[a2][vstart+i].X -= (ATOM_PAIRS[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.X/rlen3;
									FRAME_A_MATRIX[a2][vstart+i].Y -= (ATOM_PAIRS[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.Y/rlen3;
									FRAME_A_MATRIX[a2][vstart+i].Z -= (ATOM_PAIRS[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.Z/rlen3;

								}
							}
						}
					}
				}
		  	}
			
						
		}
	}
  return;
}

// FUNCTION UPDATED
static void ZCalc_InvR_Deriv(FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)		
// Calculate derivatives of the forces wrt to inverse pair distance to various powers. Stores minimum distance between a pair of atoms in mind.
{
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double x;
	int vstart;
	int curr_pair_type_idx;
	string TEMP_STR;

	double fc;
	double dfc;
	double rfac;

	for(int a1=0;a1<TRAJECTORY.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<TRAJECTORY.ATOMS;a2++)
		{
			TEMP_STR = TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(TRAJECTORY.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
			//calculate vstart:
		
			vstart = curr_pair_type_idx * ATOM_PAIRS[curr_pair_type_idx].SNUM;

			// Start with minimum image convention.  Use layers to access larger distances if desired.	
		   
			RVEC.X  = TRAJECTORY.COORDS[a2].X - TRAJECTORY.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/TRAJECTORY.BOXDIM.X )  * TRAJECTORY.BOXDIM.X;
		   
			RVEC.Y = TRAJECTORY.COORDS[a2].Y - TRAJECTORY.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/TRAJECTORY.BOXDIM.Y )  * TRAJECTORY.BOXDIM.Y;
		   
			RVEC.Z = TRAJECTORY.COORDS[a2].Z - TRAJECTORY.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/TRAJECTORY.BOXDIM.Z )  * TRAJECTORY.BOXDIM.Z;
      
			for(int n1=-1*nlayers;n1<nlayers+1;n1++)
			{
				for(int n2=-1*nlayers;n2<nlayers+1;n2++)
				{
					for(int n3=-1*nlayers;n3<nlayers+1;n3++)
					{
						RAB.X = RVEC.X + n1 * TRAJECTORY.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * TRAJECTORY.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * TRAJECTORY.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );

						if ( rlen < ATOM_PAIRS[curr_pair_type_idx].MIN_FOUND_DIST )	// spline term calculated w/cutoff:
							ATOM_PAIRS[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

						if(rlen > ATOM_PAIRS[curr_pair_type_idx].S_MINIM && rlen < ATOM_PAIRS[curr_pair_type_idx].S_MAXIM)
						{
							x = rlen;
							
							// Calculate the penalty function (fc) and its derivative
							
							fc = (rlen - ATOM_PAIRS[curr_pair_type_idx].S_MAXIM);
							dfc = fc;
							
							fc = fc*fc*fc;							
							dfc = 3*dfc*dfc;
						
							for (int i=0; i<ATOM_PAIRS[curr_pair_type_idx].SNUM; i++) 
							{
								rfac = ( (i+2) / pow(rlen,i+3) )*fc;
								rfac -= (1/pow(rlen,i+2))*dfc;

								FRAME_A_MATRIX[a1][vstart+i].X += rfac * RAB.X/rlen;
								FRAME_A_MATRIX[a1][vstart+i].Y += rfac * RAB.Y/rlen;
								FRAME_A_MATRIX[a1][vstart+i].Z += rfac * RAB.Z/rlen;
									
								FRAME_A_MATRIX[a2][vstart+i].X -= rfac * RAB.X/rlen;
								FRAME_A_MATRIX[a2][vstart+i].Y -= rfac * RAB.Y/rlen;
								FRAME_A_MATRIX[a2][vstart+i].Z -= rfac * RAB.Z/rlen;								
      						}
						}
			 	   }
				}
			}
		}
	}
	return;

}

// FUNCTION UPDATED
static void ZCalc_Poly_Deriv(FRAME & TRAJECTORY, vector<int> & ATOM_PAIR_TYPES, vector<PAIRS> & ATOM_PAIRS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)	
// Calculate derivatives of the forces wrt the DFTB Erep parameters.
{
	
	
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double x;
	int vstart;

	int curr_pair_type_idx;
	string TEMP_STR;
	
	
	
	static bool called_before = false;
	const double autoang = 0.5291772488820865;	// Conversion factor fro angstroms to au
	static double *rc;//static_cast<int>(ATOM_PAIRS.size())];
	
	int dim;
	double rfac;

	if ( ! called_before ) 
	{
		called_before = true;
		dim = 0;
		
		rc = new double [ATOM_PAIRS.size()];
		
		for ( int i=0; i<ATOM_PAIRS.size(); i++ ) 
		{
			// Convert cutoff (S_MAXIM) from units of Angstrom to bohr (atomic units), (rc)
			rc[i] = ATOM_PAIRS[i].S_MAXIM/autoang;
//			printf("rc[%d] = %lf bohr (%lf Angstr.)\n",i,rc[i], ATOM_PAIRS[i].S_MAXIM);
			cout << "	rc[" << i << "] = " << fixed << setprecision(3) << rc[i] << " bohr (" << fixed << setprecision(3) << ATOM_PAIRS[i].S_MAXIM << " Angstr.)" << endl;
			
			if ( ATOM_PAIRS[i].SNUM > dim ) 
				dim = ATOM_PAIRS[i].SNUM;
		}
		dim++;
	
	}

	// main loop for ninth order polynomial terms:
	
	for(int a1=0;a1<TRAJECTORY.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<TRAJECTORY.ATOMS;a2++)
		{
			TEMP_STR = TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(TRAJECTORY.ATOMTYPE[a2]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
			
//			cout << "SANITY CHECK--XX -- " << TRAJECTORY.ATOMTYPE[a1] << " " << TRAJECTORY.ATOMTYPE[a2] << " " << curr_pair_type_idx << " " << 
//				    ATOM_PAIRS[curr_pair_type_idx].ATM1TYP << " " << ATOM_PAIRS[curr_pair_type_idx].ATM2TYP << endl;
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
		
			vstart = curr_pair_type_idx * ATOM_PAIRS[curr_pair_type_idx].SNUM;
			
//			cout << "PAIR TYPE IS:  "  << curr_pair_type_idx << " " << TRAJECTORY.ATOMTYPE[a1] <<  " " << TRAJECTORY.ATOMTYPE[a2] << " " << endl;
//			cout << "AND VSTART IS: " << vstart << endl << endl;

			// Start with minimum image convention.  Use layers to access larger distances if desired.	
		   
			RVEC.X  = TRAJECTORY.COORDS[a2].X - TRAJECTORY.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/TRAJECTORY.BOXDIM.X )  * TRAJECTORY.BOXDIM.X;
		   
			RVEC.Y = TRAJECTORY.COORDS[a2].Y - TRAJECTORY.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/TRAJECTORY.BOXDIM.Y )  * TRAJECTORY.BOXDIM.Y;
		   
			RVEC.Z = TRAJECTORY.COORDS[a2].Z - TRAJECTORY.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/TRAJECTORY.BOXDIM.Z )  * TRAJECTORY.BOXDIM.Z;


			for(int n1=-1*nlayers;n1<nlayers+1;n1++)
			{
				for(int n2=-1*nlayers;n2<nlayers+1;n2++)
				{
			    	for(int n3=-1*nlayers;n3<nlayers+1;n3++)
					{
						RAB.X = RVEC.X + n1 * TRAJECTORY.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * TRAJECTORY.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * TRAJECTORY.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );
						
						if ( rlen < ATOM_PAIRS[curr_pair_type_idx].MIN_FOUND_DIST )	// spline term calculated w/cutoff:
							ATOM_PAIRS[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

						if(rlen > ATOM_PAIRS[curr_pair_type_idx].S_MINIM && rlen < ATOM_PAIRS[curr_pair_type_idx].S_MAXIM)
						{
							// calculate binning, convert all distances to au from angstroms 
							x = rlen/autoang;

							for (int i=0; i<ATOM_PAIRS[curr_pair_type_idx].SNUM; i++) 
							{
								
								rfac = -(i+2)*pow((rc[curr_pair_type_idx]-x),i+1);

								FRAME_A_MATRIX[a1][vstart+i].X += rfac * RAB.X/rlen;
								FRAME_A_MATRIX[a1][vstart+i].Y += rfac * RAB.Y/rlen;
								FRAME_A_MATRIX[a1][vstart+i].Z += rfac * RAB.Z/rlen;
									
								FRAME_A_MATRIX[a2][vstart+i].X -= rfac * RAB.X/rlen;
								FRAME_A_MATRIX[a2][vstart+i].Y -= rfac * RAB.Y/rlen;
								FRAME_A_MATRIX[a2][vstart+i].Z -= rfac * RAB.Z/rlen;
							}
						}
					}
				}		
			}
		}
	}
	return;

}


// FUNCTION UPDATED
void SubtractCoordForces(FRAME & TRAJECTORY, bool calc_deriv, vector<XYZ> & P_OVER_FORCES,  vector<PAIRS> & ATOM_PAIRS, vector<int> & ATOM_PAIR_TYPES_ALL, map<string,int> & PAIR_MAP)
{
	// this function subtracts the ReaxFF over-coordination term to re-fit 
	// splines/charges iteratively for self-consistence.
	// If calc_deriv is true, the derivative of the force wrt the magnitude of the
	// 3-body interaction is placed in Fderiv.  Otherwise, the 3-body force is subtracted
	// from the total forces given in force.  
	 

	// Lucas paper draft parameters.	// Chenoweth values for p1, p2, r0, lambda6 p(ovun2)
	//									// pover = 50.0;	= ATOM_PAIRS[j].OVRPRMS[0];
	// p1  	  =	-2.5881042987450e-01;	// -0.0657;			= ATOM_PAIRS[j].OVRPRMS[1];
	// r0 	  =	 9.6000387075695e-01;	//  5.0451;			= ATOM_PAIRS[j].OVRPRMS[2];
	// p2 	  =	 3.8995379237250e+00;	//  1.0165;			= ATOM_PAIRS[j].OVRPRMS[3];
	// lambda6=	-8.9;  					// -3.6141;			= ATOM_PAIRS[j].OVRPRMS[4];

	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double Vtot = 0;
	
	vector<XYZ> SForce(TRAJECTORY.ATOMS);	// WHAT IS THIS?
	vector<XYZ> dEover(TRAJECTORY.ATOMS);	// Over bonding energy derivative (force) w/r/t delta ("S")
	vector<XYZ> dFover(TRAJECTORY.ATOMS);	// Over bonding force derivative w/r/t pover
	
	// THREE-BODY POTENTIAL:
	
	double Eover = 0;			// Overbonding energy
	double p1,p2,r0,pover,lambda6;
	double rik,tempr,temps;	
	double S[TRAJECTORY.ATOMS];				// Difference between oxygen bond order and valency (Called delta in the paper)
	 
	int i_pair; 
	int curr_pair_type_idx; 
	string TEMP_STR;
	
	// Set up some variables
	
	for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)
	{
		dEover[a1].X = 0.0;
		dEover[a1].Y = 0.0;
		dEover[a1].Z = 0.0;
		
		dFover[a1].X = 0.0;
		dFover[a1].Y = 0.0;
		dFover[a1].Z = 0.0;		
		
		SForce[a1].X = 0;
		SForce[a1].Y = 0;
		SForce[a1].Z = 0;
		
		S[a1] = 0.0;				
	}
		
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Calculate the overbonding energy... Derivative needs to be in separate double loop b/c depends
	// on a term calculated with energies ("S", which is called delta in the paper)
	//
	////////////////////////////////////////////////////////////////////////////////////////////////
			
	for(int ai=0;ai<TRAJECTORY.ATOMS;ai++) // Calculates Eover... I'm assuming this is primarily used for MD, which is why I don't have it
	{
		temps = 0.0;
	  
		for(int ak=0;ak<TRAJECTORY.ATOMS;ak++)
		{			
			if (ai != ak)
			{
				// TWO-BODY PART... ONLY CARES ABOUT LOOPS OVER AI AND AK
			
				TEMP_STR = TRAJECTORY.ATOMTYPE[ai];
				TEMP_STR.append(TRAJECTORY.ATOMTYPE[ak]);
				curr_pair_type_idx = PAIR_MAP[TEMP_STR];
				
				// THE ITEMS COMMENTED OUT IN THIS LOOP SHOULD BE PUT BACK ONCE CODE COMPARISON WITH LUCAS' COMPLETE!!!!!
					
				if(ATOM_PAIRS[curr_pair_type_idx].USE_OVRPRMS && (TRAJECTORY.ATOMTYPE[ai] == ATOM_PAIRS[curr_pair_type_idx].OVER_TO_ATM)) // Then we should have a defined pair type
				{
					// Start with minimum image convention.  Use layers to access larger distances if desired.	
		   
					RVEC.X  = TRAJECTORY.COORDS[ak].X - TRAJECTORY.COORDS[ai].X;
					RVEC.X -= floor(0.5 + RVEC.X/TRAJECTORY.BOXDIM.X )  * TRAJECTORY.BOXDIM.X;

					RVEC.Y = TRAJECTORY.COORDS[ak].Y - TRAJECTORY.COORDS[ai].Y;
					RVEC.Y -= floor(0.5 +  RVEC.Y/TRAJECTORY.BOXDIM.Y )  * TRAJECTORY.BOXDIM.Y;
		   
					RVEC.Z = TRAJECTORY.COORDS[ak].Z - TRAJECTORY.COORDS[ai].Z;
					RVEC.Z -= floor(0.5 +  RVEC.Z/TRAJECTORY.BOXDIM.Z )  * TRAJECTORY.BOXDIM.Z;

				
					rik = sqrt( RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z);	
				
				// Calculate the O--H bond order as defined by ReaxFF
				// 
					temps +=  exp(ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[2]
					*pow(rik/ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[1],ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[3])); 
				}
			}
		}
		
		S[ai]  = temps - 2.0; // This is the reaxff delta (diff between oxygen valency and bond order) 
		Eover += ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[0] * S[ai] * 1.0/(1.0+exp(ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[4]*S[ai]));	

	}

	Vtot += Eover;
	
	///////////////////////////////////////////////////////
	//
	// Now calculate the forces (derivative)
	//
	///////////////////////////////////////////////////////

	// Here, the first loop is over all ai, but we don't do anything unless ai is an oxygen type, since overbonding forces in this program
	// are defined to be between oxygen an either oxygen or H 
	// Second loop is over all atoms != ai
	// For force between a pair of atoms, recall that vectors are defined from ai to ak

	for(int ai=0; ai<TRAJECTORY.ATOMS; ai++)
	{

		for(int ak=0; ak<TRAJECTORY.ATOMS; ak++)
		{	
			TEMP_STR = TRAJECTORY.ATOMTYPE[ai];
			TEMP_STR.append(TRAJECTORY.ATOMTYPE[ak]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
			
			// If it's a pair type with defined over params
			// and the pair isn't comprised of the exact same atom
			// and ai is the "to" atom (i.e. O in an O--H pair)
		
			if((ATOM_PAIRS[curr_pair_type_idx].USE_OVRPRMS && ai != ak && TRAJECTORY.ATOMTYPE[ai] == ATOM_PAIRS[curr_pair_type_idx].OVER_TO_ATM))
			{

				RVEC.X  = TRAJECTORY.COORDS[ak].X - TRAJECTORY.COORDS[ai].X;
				RVEC.X -= floor(0.5 + RVEC.X/TRAJECTORY.BOXDIM.X )  * TRAJECTORY.BOXDIM.X;
		   
				RVEC.Y = TRAJECTORY.COORDS[ak].Y - TRAJECTORY.COORDS[ai].Y;
				RVEC.Y -= floor(0.5 + RVEC.Y/TRAJECTORY.BOXDIM.Y )  * TRAJECTORY.BOXDIM.Y;
		   
				RVEC.Z = TRAJECTORY.COORDS[ak].Z - TRAJECTORY.COORDS[ai].Z;
				RVEC.Z -= floor(0.5 + RVEC.Z/TRAJECTORY.BOXDIM.Z )  * TRAJECTORY.BOXDIM.Z;
				
				rik = sqrt( RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z);				
				
				// Calculate the derivative of Eover w/r/t delta. Terms emerge from quotient rule combined with chain rule.

				tempr = 1.0/(1+exp(ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[4]*S[ai])) - ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[4]*S[ai]*exp(ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[4]*S[ai])
					        / pow(1.0+exp(ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[4]*S[ai]),2);

				tempr *= ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[2]*pow(rik/ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[1],ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[3]) *
					     ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[3]*exp(ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[2] *
						 pow(rik/ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[1],ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[3]))/rik;	


				dEover[ai].X -= tempr*RVEC.X/rik*ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[0];
				dEover[ai].Y -= tempr*RVEC.Y/rik*ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[0];
				dEover[ai].Z -= tempr*RVEC.Z/rik*ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[0];						

				dEover[ak].X += tempr*RVEC.X/rik*ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[0];
				dEover[ak].Y += tempr*RVEC.Y/rik*ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[0];
				dEover[ak].Z += tempr*RVEC.Z/rik*ATOM_PAIRS[curr_pair_type_idx].OVRPRMS[0];	
				
				dFover[ai].X -= tempr*RVEC.X/rik;
				dFover[ai].Y -= tempr*RVEC.Y/rik;
				dFover[ai].Z -= tempr*RVEC.Z/rik;						

				dFover[ak].X += tempr*RVEC.X/rik;
				dFover[ak].Y += tempr*RVEC.Y/rik;
				dFover[ak].Z += tempr*RVEC.Z/rik;					


			}	

		}


	}

    for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)
	{
        SForce[a1].X -= dEover[a1].X;
		SForce[a1].Y -= dEover[a1].Y;
		SForce[a1].Z -= dEover[a1].Z;	
	}

    if ( !calc_deriv) // subtract Coord. force:
	{ 
		for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)
		{
	        TRAJECTORY.FORCES[a1].X -= SForce[a1].X;
			TRAJECTORY.FORCES[a1].Y -= SForce[a1].Y;
			TRAJECTORY.FORCES[a1].Z -= SForce[a1].Z;	
		}
    }
    else // Then we're fitting pover... Calculate derivative of 3-body force wrt pover.
	{
		for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)
		{
	        P_OVER_FORCES[a1].X -= dFover[a1].X;
			P_OVER_FORCES[a1].Y -= dFover[a1].Y;
			P_OVER_FORCES[a1].Z -= dFover[a1].Z;		
		}
	} // VERIFIED FROM HERE UP
}





