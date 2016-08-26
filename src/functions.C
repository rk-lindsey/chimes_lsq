#include<iomanip>
#include<iostream>
#include<fstream>
#include "functions.h"

using namespace std;

//////////////////////////////////////////
//
//	SMALL UTILITY FUNCTION
//
//////////////////////////////////////////

// Overloaded error message exit functions

void EXIT_MSG(string EXIT_STRING) // Single error message
{
	cout << EXIT_STRING << endl;
	exit(0);
}

void EXIT_MSG(string EXIT_STRING, string EXIT_VAR) // error message: var
{
	cout << EXIT_STRING << ": " << EXIT_VAR << endl;
	exit(0);
}

void EXIT_MSG(string EXIT_STRING, double EXIT_VAR) // error message: var
{
	cout << EXIT_STRING << ": " << EXIT_VAR << endl;
	exit(0);
}


void SET_3B_CHEBY_POLYS( PAIRS & FF_2BODY, double *Tn, double *Tnd, const double rlen, double & xdiff)
// Sets the value of the Chebyshev polynomials (Tn) and thier derivatives (Tnd)
{
	double xmin, xmax, xavg;
	double x;
		
	
	if ( FF_2BODY.CHEBY_TYPE == "INVRSE_R" ) 
	{
		xavg  =  0.5 * (1.0/FF_2BODY.S_MINIM + 1.0/FF_2BODY.S_MAXIM); 
		xdiff =  0.5 * (1.0/FF_2BODY.S_MINIM - 1.0/FF_2BODY.S_MAXIM); 
		x     = (1.0/rlen-xavg) / xdiff;	
	} 
	else if ( FF_2BODY.CHEBY_TYPE == "MORSE" ) 
	{
		xmin  = exp(-FF_2BODY.S_MAXIM/FF_2BODY.LAMBDA); 
		xmax  = exp(-FF_2BODY.S_MINIM/FF_2BODY.LAMBDA); 
		xavg  = 0.5 * (xmin + xmax);
		xdiff = 0.5 * (xmax - xmin);
		x = (exp(-rlen/FF_2BODY.LAMBDA)-xavg)/xdiff;
	}
	else if (FF_2BODY.CHEBY_TYPE == "DEFAULT")
	{
		xavg  = 0.5 * (FF_2BODY.S_MINIM + FF_2BODY.S_MAXIM);
		xdiff = 0.5 * (FF_2BODY.S_MAXIM - FF_2BODY.S_MINIM);
		x = (rlen-xavg) / xdiff;	
	}
	else
	{
		cout << "ERROR: Undefined CHBTYPE: " << FF_2BODY.PRPR_NM << endl;
		cout << "       Excepted values are \"DEFAULT\", \"INVRSE_R\", or \"MORSE\". " << endl;
		exit(1);
	}

	if ( x < -1.0 ) // Why are we setting it as -1? shouldn't something else be done here instead? ...since we're outside the range of the fit??
	{
		cout << "Warning: In 3B Cheby transformation, r < rmin " << FF_2BODY.CHEBY_TYPE << endl;
		x = -1.0;
	}
	else if ( x > 1.0 ) // Added a warning...Do we need to turn off the interaction, because it's beyond rcut?
	{
		cout << "Warning:  r > rmax" << endl;								
		x = 1.0;								
	}

	Tn[0] = 1.0;
	Tn[1] = x;

	Tnd[0] = 1.0;
	Tnd[1] = 2.0 * x;

	for ( int i = 2; i <= FF_2BODY.SNUM_3B_CHEBY; i++ ) 
	{
		Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2];
		Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2];
	}

	for ( int i = FF_2BODY.SNUM_3B_CHEBY; i >= 1; i-- ) 
		Tnd[i] = i * Tnd[i-1];

	Tnd[0] = 0.0;
}
//void SET_3B_CHEBY_POLYS( vector<PAIR_FF> & FF_2BODY, const int curr_pair_type_idx, double *Tn, double *Tnd, const double rlen, double & xdiff)

// FUNCTION IS OVERLOADED
void SET_3B_CHEBY_POWERS(vector<PAIR_FF> & FF_2BODY, TRIPLETS & FF_3BODY, map<string,int> & PAIR_MAP,  int & pow_ij, int & pow_ik, int & pow_jk, string PAIR_TYPE_IJ, string PAIR_TYPE_IK, string PAIR_TYPE_JK, int POWER_SET) // LSQ version
// Matches the allowed powers to the ij. ik, jk type pairs formed from the atom triplet ai, aj, ak 
{
	if    (FF_2BODY[ PAIR_MAP[PAIR_TYPE_IJ] ].PRPR_NM == FF_3BODY.ATMPAIR1)
	{
		if    (FF_2BODY[ PAIR_MAP[PAIR_TYPE_IK] ].PRPR_NM == FF_3BODY.ATMPAIR2)
		{
			// Then jk = z
			
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
		}
		else
		{
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
		}

	}
	else if(FF_2BODY[ PAIR_MAP[PAIR_TYPE_IJ] ].PRPR_NM == FF_3BODY.ATMPAIR2)
	{
		if    (FF_2BODY[ PAIR_MAP[PAIR_TYPE_IK] ].PRPR_NM == FF_3BODY.ATMPAIR1)
		{
			// Then jk = z
			
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
		}
		else
		{
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
		}										
	}
	else // PAIR_TYPE_IJ matches alloweD[i].z
	{
		if    (FF_2BODY[ PAIR_MAP[PAIR_TYPE_IK] ].PRPR_NM == FF_3BODY.ATMPAIR1)
		{
			// Then jk = z
			
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
		}
		else
		{
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
		}										
	}
	
}

void SET_3B_CHEBY_POWERS(vector<PAIRS> & FF_2BODY, TRIPLETS & FF_3BODY, map<string,int> & PAIR_MAP,  int & pow_ij, int & pow_ik, int & pow_jk, string PAIR_TYPE_IJ, string PAIR_TYPE_IK, string PAIR_TYPE_JK, int POWER_SET) // MD version
// Matches the allowed powers to the ij. ik, jk type pairs formed from the atom triplet ai, aj, ak 
{
	if    (FF_2BODY[ PAIR_MAP[PAIR_TYPE_IJ] ].PRPR_NM == FF_3BODY.ATMPAIR1)
	{
		if    (FF_2BODY[ PAIR_MAP[PAIR_TYPE_IK] ].PRPR_NM == FF_3BODY.ATMPAIR2)
		{
			// Then jk = z
			
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
		}
		else
		{
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
		}

	}
	else if(FF_2BODY[ PAIR_MAP[PAIR_TYPE_IJ] ].PRPR_NM == FF_3BODY.ATMPAIR2)
	{
		if    (FF_2BODY[ PAIR_MAP[PAIR_TYPE_IK] ].PRPR_NM == FF_3BODY.ATMPAIR1)
		{
			// Then jk = z
			
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
		}
		else
		{
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
		}										
	}
	else // PAIR_TYPE_IJ matches alloweD[i].z
	{
		if    (FF_2BODY[ PAIR_MAP[PAIR_TYPE_IK] ].PRPR_NM == FF_3BODY.ATMPAIR1)
		{
			// Then jk = z
			
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
		}
		else
		{
			pow_ij = FF_3BODY.ALLOWED_POWERS[POWER_SET].Z;
			pow_ik = FF_3BODY.ALLOWED_POWERS[POWER_SET].Y;
			pow_jk = FF_3BODY.ALLOWED_POWERS[POWER_SET].X;
		}										
	}
	
}

void REPLICATE_SYSTEM(const FRAME & SYSTEM, FRAME & REPLICATE)
{
	REPLICATE.ATOMS    				= SYSTEM.ATOMS;
	REPLICATE.BOXDIM.X 				= SYSTEM.BOXDIM.X;
	REPLICATE.BOXDIM.Y 				= SYSTEM.BOXDIM.Y;
	REPLICATE.BOXDIM.Z 				= SYSTEM.BOXDIM.Z;
	REPLICATE.TEMPERATURE			= SYSTEM.TEMPERATURE;
	REPLICATE.PRESSURE 				= SYSTEM.PRESSURE;
	REPLICATE.PRESSURE_XYZ 			= SYSTEM.PRESSURE_XYZ;
	REPLICATE.PRESSURE_TENSORS.X	= SYSTEM.PRESSURE_TENSORS.X;
	REPLICATE.PRESSURE_TENSORS.Y 	= SYSTEM.PRESSURE_TENSORS.Y;
	REPLICATE.PRESSURE_TENSORS.Z	= SYSTEM.PRESSURE_TENSORS.Z;
	REPLICATE.TOT_POT_ENER 			= SYSTEM.TOT_POT_ENER;
	
	REPLICATE.COORDS      .resize(REPLICATE.ATOMS);
	REPLICATE.CHARGES     .resize(REPLICATE.ATOMS);
	REPLICATE.MASS        .resize(REPLICATE.ATOMS);
	REPLICATE.FORCES      .resize(REPLICATE.ATOMS);
	REPLICATE.ACCEL       .resize(REPLICATE.ATOMS);
	REPLICATE.TMP_EWALD   .resize(REPLICATE.ATOMS);
	REPLICATE.VELOCITY    .resize(REPLICATE.ATOMS);
	REPLICATE.VELOCITY_NEW.resize(REPLICATE.ATOMS);
	
	for(int i=0; i<REPLICATE.ATOMS; i++)
	{
		REPLICATE.ATOMTYPE			= SYSTEM.ATOMTYPE;
		
		REPLICATE.COORDS[i].X = SYSTEM.COORDS[i].X;
		REPLICATE.COORDS[i].Y = SYSTEM.COORDS[i].Y;
		REPLICATE.COORDS[i].Z = SYSTEM.COORDS[i].Z;
		
		REPLICATE.CHARGES[i] = SYSTEM.CHARGES[i];
		
		REPLICATE.MASS[i] = SYSTEM.MASS[i];
		
		REPLICATE.FORCES[i].X = SYSTEM.FORCES[i].X;
		REPLICATE.FORCES[i].Y = SYSTEM.FORCES[i].Y;
		REPLICATE.FORCES[i].Z = SYSTEM.FORCES[i].Z;
		
		REPLICATE.ACCEL[i].X = SYSTEM.ACCEL[i].X;
		REPLICATE.ACCEL[i].Y = SYSTEM.ACCEL[i].Y;
		REPLICATE.ACCEL[i].Z = SYSTEM.ACCEL[i].Z;
		
		REPLICATE.VELOCITY[i].X = SYSTEM.VELOCITY[i].X;
		REPLICATE.VELOCITY[i].Y = SYSTEM.VELOCITY[i].Y;
		REPLICATE.VELOCITY[i].Z = SYSTEM.VELOCITY[i].Z;
		
		REPLICATE.VELOCITY_NEW[i].X = SYSTEM.VELOCITY_NEW[i].X;
		REPLICATE.VELOCITY_NEW[i].Y = SYSTEM.VELOCITY_NEW[i].Y;
		REPLICATE.VELOCITY_NEW[i].Z = SYSTEM.VELOCITY_NEW[i].Z;
	}
}

double numerical_pressure(const FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP)   
// Evaluates the configurational part of the pressure numerically by -dU/dV.
// Essentially, we are taking the system and expanding/contracting it a bit (lscale) to get the change in potential energy 
{

	const double eps = 1.0e-04;
	double lscale;
	double Vtot1, Vtot2;
	double Vol1, Vol2;

	// Make a copy of the system
	
	static FRAME REPLICATE;
	
	REPLICATE_SYSTEM(SYSTEM, REPLICATE);		// Can we make this one of those "if ! called_before sorts of variables?"

	// Expand coords by a bit (lscale)

	lscale = 1.0  + eps;
	
	for (int j=0; j<REPLICATE.ATOMS; j++ ) 
	{
	 	REPLICATE.COORDS[j].X = lscale * SYSTEM.COORDS[j].X;
		REPLICATE.COORDS[j].Y = lscale * SYSTEM.COORDS[j].Y;
		REPLICATE.COORDS[j].Z = lscale * SYSTEM.COORDS[j].Z;
	}
	
	REPLICATE.BOXDIM.X = lscale * SYSTEM.BOXDIM.X;
	REPLICATE.BOXDIM.Y = lscale * SYSTEM.BOXDIM.Y;
	REPLICATE.BOXDIM.Z = lscale * SYSTEM.BOXDIM.Z;

	// Compute/store new total potential energy and volume
	
	Vol1 = REPLICATE.BOXDIM.X * REPLICATE.BOXDIM.Y * REPLICATE.BOXDIM.Z;
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP);
	
	Vtot1 = REPLICATE.TOT_POT_ENER;
	
	// Contract coords by a bit

	lscale = 1.0  - eps;
	
	for (int j=0; j<REPLICATE.ATOMS; j++ ) 
	{
	 	REPLICATE.COORDS[j].X = lscale * SYSTEM.COORDS[j].X;
		REPLICATE.COORDS[j].Y = lscale * SYSTEM.COORDS[j].Y;
		REPLICATE.COORDS[j].Z = lscale * SYSTEM.COORDS[j].Z;
	}
	
	REPLICATE.BOXDIM.X = lscale * SYSTEM.BOXDIM.X;
	REPLICATE.BOXDIM.Y = lscale * SYSTEM.BOXDIM.Y;
	REPLICATE.BOXDIM.Z = lscale * SYSTEM.BOXDIM.Z;

	// Compute/store new total potential energy and volume
	
	Vol2 = REPLICATE.BOXDIM.X * REPLICATE.BOXDIM.Y * REPLICATE.BOXDIM.Z;
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP);
	
	Vtot2 = REPLICATE.TOT_POT_ENER;

	// compute (return) pressure 
	
	return -(Vtot2 - Vtot1)/(Vol2 - Vol1);

}


//////////////////////////////////////////
//
//	FUNCTION HEADERS -- DERIVATIVE CALCULATION
//
//////////////////////////////////////////

static void ZCalc_Spline_Deriv  (FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);

static void ZCalc_Cheby_Deriv   (FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);

static void ZCalc_3B_Cheby_Deriv(FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<TRIPLETS> & PAIR_TRIPLETS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, map<string,int> TRIAD_MAP);

static void ZCalc_Poly_Deriv    (FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);	

static void ZCalc_InvR_Deriv    (FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);



//////////////////////////////////////////
//
//	FUNCTION HEADERS -- FORCE CALCULATION
//
//////////////////////////////////////////

static void ZCalc_Cheby(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP);

static void ZCalc_3B_Cheby(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP);

static void ZCalc_Lj(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP);

static void ZCalc_Spline(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP);

static void ZCalcSR_Over(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP);

/*

static void ZCalc_SR_Analytic(double **Coord,const char *Lbc, double *Latcons,const int nlayers, const int nat,const double *smin,const double *smax, const int *snum, double **SForce,double& Vtot,double& Pxyz, double *params);
			 
static void ZCalc_Stillinger(double **Coord, const char *Lbc, double *Latcons,const int nlayers, const int nat, const double *smax, double **SForce,double& Vtot,double& Pxyz);

*/


////////////////////////////////////////////////////////////
//
// FUNCTIONS -- DERIVATIVE CALCULATION
//
////////////////////////////////////////////////////////////


// FUNCTION UPDATED
void ZCalc_Deriv (vector<PAIRS> & FF_2BODY,  vector<TRIPLETS> PAIR_TRIPLETS, FRAME & FRAME_SYSTEM, vector<vector <XYZ > > & FRAME_A_MATRIX, vector<vector <XYZ > > & FRAME_COULOMB_FORCES, const int nlayers, bool if_3b_cheby, map<string,int> & PAIR_MAP,  map<string,int> & TRIAD_MAP)
// Controls which functions are used to calculate derivatives
{
  	// Will ewald calculations be needed? DFTB doesn't use the "default"
    // Ewald calculation because it has its own special way of dealing with charges
		
	bool if_ewald;
  
    if (FF_2BODY[0].PAIRTYP != "DFTBPOLY") 
  	  if_ewald = true;
    else
  	  if_ewald = false;

	if ( if_ewald ) 					// FUNCTION UPDATED
		ZCalc_Ewald_Deriv(FRAME_SYSTEM, FF_2BODY, FRAME_COULOMB_FORCES, PAIR_MAP);	
	
    if ( FF_2BODY[0].PAIRTYP == "SPLINE" )			// FUNCTION UPDATED
		ZCalc_Spline_Deriv(FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP);
	
	else if ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" )	// FUNCTION UPDATED
	{
		ZCalc_Cheby_Deriv(FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP);
	
		if (if_3b_cheby)
			ZCalc_3B_Cheby_Deriv(FRAME_SYSTEM, FF_2BODY, PAIR_TRIPLETS, FRAME_A_MATRIX, nlayers, PAIR_MAP, TRIAD_MAP);		
	}			

    else if ( FF_2BODY[0].PAIRTYP == "DFTBPOLY" )	// FUNCTION UPDATED
		ZCalc_Poly_Deriv(FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP);

    else if ( FF_2BODY[0].PAIRTYP == "INVRSE_R" )	// FUNCTION UPDATED
		ZCalc_InvR_Deriv(FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP);
		
	
		
    else 
    {
		cout << "Error: bad pairtype in ZCalc_Deriv\n";
		exit(1);
    }		
}	

// FUNCTION UPDATED
static void ZCalc_Spline_Deriv(FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)		   
// Original comment: Calculate derivatives of the forces wrt the spline parameters. Stores minimum distance between a pair of atoms in minD[i].
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

	for(int a1=0;a1<SYSTEM.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<SYSTEM.ATOMS;a2++)
		{
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
			//calculate vstart:

			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;
			 
			// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
		   
			RVEC.X  = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
		   
			RVEC.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		   
			RVEC.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;

      
			for(int n1=-1*nlayers;n1<nlayers+1;n1++)
			{
				for(int n2=-1*nlayers;n2<nlayers+1;n2++)
				{
					for(int n3=-1*nlayers;n3<nlayers+1;n3++)
					{
						RAB.X = RVEC.X + n1 * SYSTEM.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * SYSTEM.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * SYSTEM.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );

						if ( rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST ) 
							FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;
					   
					   
						if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM and rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM) //spline term calculated w/cutoff:
						{
							
							// Setting up for Hermite cubic quadrature:
							
						 	// This part: figure out which distance bin ("unit") that rlen falls into. Use the value of that bin "t" (rather than rlen ("x") ) to 
							// populate the matrix table. Use the lower and upper bounds of the bin to deterime the two constraining data points for the interpolation
							// polynomial fitting
							 
							 
							k0 = int( floor( (rlen - FF_2BODY[curr_pair_type_idx].S_MINIM) / FF_2BODY[curr_pair_type_idx].S_DELTA ) ); // binned location of distance along allowed range

						   
							if ( k0 >= FF_2BODY[curr_pair_type_idx].SNUM/2 - 1 ) // Keep from overrunning table.
							{
								cout << "Table overrun: rlen = " << rlen << endl;
								continue;
							}

 						   	x = rlen;
	    
 						   	x0 = FF_2BODY[curr_pair_type_idx].S_MINIM + FF_2BODY[curr_pair_type_idx].S_DELTA * k0;
 						   	t  = (x-x0) / FF_2BODY[curr_pair_type_idx].S_DELTA;

						   	// for classical cubic quadrature, you would only change h00, h10, h01, h11 polynomials!
						   
						  	// Setup the basis functions for cubic hermite polynomial interpolation...
						   	//
						   	// See: https://www3.nD[i].edu/~coast/jjwteach/www/www/30125/pdfnotes/lecture5_9v14.pdf
						   	// Look for "Each basis function is a third degree polynomial" 
						   	//
						   	// The fitting polynomial, g(x) = f(x_0)*alpha_0 + f(x_1)*alpha_1 + f'(x_0)*beta_0 + f'(x_1)*beta_1
						   	// But we're using t instead of x, and t_0 and t_1 are the two termini of the distance bin.
						   	// i.e. if x = 5.54, delta = 0.1, fall into a bin of value 5.55 with endpoints 5.5 and 5.6.
					   
						   	h00  =  2 * t * t * t - 3 * t * t + 1;			// alpha_0.. basis function for first point
						   
						   	h10  =  t * t * t  -2 * t * t + t;				// beta_0..  basis function for the first derivative at first point
						   	h10 *=  FF_2BODY[curr_pair_type_idx].S_DELTA;	// derivative terms have extra factor.
						   
						   	h01  = -2 * t * t * t + 3 * t * t;				// alpha_1.. basis function for second point
						   
						   	h11  =  t * t * t - t * t;						// beta_1..  basis function for first derivative at second point
						   	h11 *=  FF_2BODY[curr_pair_type_idx].S_DELTA;	// derivative terms have extra factor.

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
static void ZCalc_Cheby_Deriv(FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)	
// Calculate derivatives of the forces wrt the Chebyshev parameters. Stores minimum distance between a pair of atoms in minD[i].
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
		
		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			if (FF_2BODY[i].SNUM > dim ) 
				dim = FF_2BODY[i].SNUM;	 
		
		dim++;
		Tn   = new double [dim];
		Tnd  = new double [dim];
	}

	// Main loop for Chebyshev terms:
	
	string TEMP_STR;
	int curr_pair_type_idx;

	for(int a1=0;a1<SYSTEM.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<SYSTEM.ATOMS;a2++)
		{
			
			
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
		
			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
		   
			RVEC.X  = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
		   
			RVEC.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		   
			RVEC.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;
      
			for(int n1=-1*nlayers;n1<nlayers+1;n1++)
			{
				for(int n2=-1*nlayers;n2<nlayers+1;n2++)
				{
					for(int n3=-1*nlayers;n3<nlayers+1;n3++)
					{
						RAB.X = RVEC.X + n1 * SYSTEM.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * SYSTEM.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * SYSTEM.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );

						if ( rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST )	
							FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

						if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM and rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
						{
							// Chebyshev polynomials are only defined on the range [-1,1], so transorm the pair distance
							// in a way that allows it to fall along that range. Options are:
							//
							// x = 1/pair_dist				// Inverse r, 
							// x = exp(pair_dist/lambda)	// Morse-type
							// x = pair_dist				// default type
							// 
							// All types are normalized by s_min to s_max range to fall along [-1,1]
														
							
							if ( FF_2BODY[curr_pair_type_idx].CHEBY_TYPE == "INVRSE_R" ) 
							{
								xavg  =  0.5 * (1.0/FF_2BODY[curr_pair_type_idx].S_MINIM + 1.0/FF_2BODY[curr_pair_type_idx].S_MAXIM); // midpoint of possible pair distances in r^-1 space
								xdiff =  0.5 * (1.0/FF_2BODY[curr_pair_type_idx].S_MINIM - 1.0/FF_2BODY[curr_pair_type_idx].S_MAXIM); // width of possible pair distances in r^-1 space
								x     = (1.0/rlen-xavg) / xdiff;																		  // pair distances in r^-1 space, normalized to fit over [-1,1]
							} 
							else if ( FF_2BODY[curr_pair_type_idx].CHEBY_TYPE == "MORSE" ) 
							{
								xmin  = exp(-FF_2BODY[curr_pair_type_idx].S_MAXIM/FF_2BODY[curr_pair_type_idx].LAMBDA); 
								xmax  = exp(-FF_2BODY[curr_pair_type_idx].S_MINIM/FF_2BODY[curr_pair_type_idx].LAMBDA); 
								xavg  = 0.5 * (xmin + xmax);																// midpoint of possible pair distances in morse space
								xdiff = 0.5 * (xmax - xmin);																// width of possible pair distances in morse space
								x = (exp(-rlen/FF_2BODY[curr_pair_type_idx].LAMBDA)-xavg)/xdiff;							// pair distances in morse space, normalized to fit over [-1,1]
							}
							else if (FF_2BODY[curr_pair_type_idx].CHEBY_TYPE == "DEFAULT")
							{
								xavg  = 0.5 * (FF_2BODY[curr_pair_type_idx].S_MINIM + FF_2BODY[curr_pair_type_idx].S_MAXIM); // midpoint of possible pair distances
								xdiff = 0.5 * (FF_2BODY[curr_pair_type_idx].S_MAXIM - FF_2BODY[curr_pair_type_idx].S_MINIM); // width of possible pair distances
								x = (rlen-xavg) / xdiff;																		 // pair distances, normalized to fit over [-1,1]
							}
							else
							{
								cout << "ERROR: Undefined CHBTYPE: " << FF_2BODY[curr_pair_type_idx].CHEBY_TYPE << endl;
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
								cout << "Warning:  r > rmax" << endl;								
								x = 1.0;								
							}

							// Generate Chebyshev polynomials by recursion. 
							// 
							// What we're doing here. Want to fit using Cheby polynomials of the 1st kinD[i]. "T_n(x)."
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
						  
							for ( int i = 2; i <= FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
							{
								Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2];
								Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2];
							}
							
							// Now multiply by n to convert Tnd's to actual derivatives of Tn

							for ( int i = FF_2BODY[curr_pair_type_idx].SNUM; i >= 1; i-- ) 
								Tnd[i] = i * Tnd[i-1];

							Tnd[0] = 0.0;
							
							// Define/setup a penalty function to deal with the cases where the pair distance is close to rmin.. to discourage
							// the SYSTEM from heading towards though poorly sampled regions of the PES.

							rlen3 = rlen * rlen * rlen;
							
							if ( FF_2BODY[curr_pair_type_idx].CHEBY_TYPE == "MORSE" ) // IS THE MORSE TYPE THE ONLY ONE THAT USES DERIVATIVES?? WHY??
							{
								// fcut and fcutderv are the form that the penalty func and its derivative for the morse-type pair distance transformation
								
								fcut0     = (1.0 - rlen/FF_2BODY[curr_pair_type_idx].S_MAXIM);
								fcut      = fcut0 * fcut0 * fcut0;
								fcutderiv = -3.0 * fcut0 * fcut0 / FF_2BODY[curr_pair_type_idx].S_MAXIM;
							
								for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
								{
									
									// NOTE: All these extra terms are coming from:
									//
									// 1. Chain rule to account for transformation from morse-type pair distance to x
									// 2. Product rule coming from pair distance dependence of fcut, the penalty function
									
									tmp_doub = (fcut * Tnd[i+1] *(-exp(-rlen/FF_2BODY[curr_pair_type_idx].LAMBDA)/FF_2BODY[curr_pair_type_idx].LAMBDA)/xdiff + fcutderiv * Tn[i+1] );
									
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
								for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
								{
									FRAME_A_MATRIX[a1][vstart+i].X += (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.X/rlen3;
									FRAME_A_MATRIX[a1][vstart+i].Y += (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.Y/rlen3;
									FRAME_A_MATRIX[a1][vstart+i].Z += (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.Z/rlen3;
									
									FRAME_A_MATRIX[a2][vstart+i].X -= (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.X/rlen3;
									FRAME_A_MATRIX[a2][vstart+i].Y -= (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.Y/rlen3;
									FRAME_A_MATRIX[a2][vstart+i].Z -= (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * RAB.Z/rlen3;

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

static void ZCalc_3B_Cheby_Deriv(FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<TRIPLETS> & PAIR_TRIPLETS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, map<string,int> TRIAD_MAP)		
// Calculate derivatives of the forces wrt the 3-body Chebyshev parameters. 
// This three body interaction stems from: C_n^ij *  C_n^ik * C_n^jk * T_n(x_ij) * T_n(x_ik) * T_n(x_jk)
{
/*
	The logic:
	+ Run a tripe loop over all atoms in the system.
	+ Compute C_ij, C_ik, and C_jk coeffiecients independently as you would do for a normal 2 body 
*/

	XYZ RVEC_IJ;
	XYZ RVEC_IK;
	XYZ RVEC_JK;
	 		
	XYZ RAB_IJ;
	XYZ RAB_IK;
	XYZ RAB_JK; 		

	double rlen_ij,  rlen_ik,  rlen_jk;
	double rlen3_ij, rlen3_ik, rlen3_jk;
	int vstart;
    static int n_2b_cheby_terms, n_3b_cheby_terms;
	static double *Tn_ij,  *Tn_ik,  *Tn_jk;
	static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
	static bool called_before = false;
	
	static int pow_ij, pow_ik, pow_jk;

	
	double x, xmin, xmax, xavg, xdiff;
	double xdiff_ij, xdiff_ik, xdiff_jk; 
	double rlen3;				  
	double fcut0, fcut0_ij, fcut0_ik, fcut0_jk; 
	double fcut,  fcut_ij,  fcut_ik,  fcut_jk; 			
	double deriv, deriv_ij, deriv_ik, deriv_jk;
	double fcutderiv, fcutderiv_ij, fcutderiv_ik, fcutderiv_jk; 	
	double tmp_doub; 	
	
	static string TEMP_STR;
	static string PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK;
	static int curr_triple_type_index;
	static int curr_pair_type_idx_ij;
	static int curr_pair_type_idx_ik;
	static int curr_pair_type_idx_jk;
	static int row_offset;	
	
	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;

		
		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
		{
			if (FF_2BODY[i].SNUM_3B_CHEBY > dim ) 
				dim = FF_2BODY[i].SNUM_3B_CHEBY;	
			
			n_2b_cheby_terms += FF_2BODY[i].SNUM;
		}
		for ( int i = 0; i < PAIR_TRIPLETS.size(); i++ ) 
			n_3b_cheby_terms += PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS;;
		
		dim++;
		
		Tn_ij   = new double [dim];
		Tn_ik   = new double [dim];
		Tn_jk   = new double [dim];

		Tnd_ij  = new double [dim];
		Tnd_ik  = new double [dim];
		Tnd_jk  = new double [dim];
	}

	// Main loop for Chebyshev terms:

	for(int a1=0;a1<SYSTEM.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		
		for(int a2=a1+1;a2<SYSTEM.ATOMS;a2++)
		{	
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);		
			PAIR_TYPE_IJ = TEMP_STR;					
			curr_pair_type_idx_ij = PAIR_MAP[TEMP_STR];
			
			// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
			// NOTE: LAYERS ARE DISABLED FOR 3B Cheby. Loops become unwieldly... 

			RVEC_IJ.X = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC_IJ.X -= floor( 0.5 + RVEC_IJ.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
			
			RVEC_IJ.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC_IJ.Y -= floor( 0.5 + RVEC_IJ.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
						
			RVEC_IJ.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			RVEC_IJ.Z -= floor( 0.5 + RVEC_IJ.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;			
			
			RAB_IJ.X = RVEC_IJ.X;
			RAB_IJ.Y = RVEC_IJ.Y;
			RAB_IJ.Z = RVEC_IJ.Z;

			rlen_ij = sqrt( RAB_IJ.X*RAB_IJ.X + RAB_IJ.Y*RAB_IJ.Y + RAB_IJ.Z*RAB_IJ.Z );			

			if(rlen_ij > FF_2BODY[curr_pair_type_idx_ij].S_MINIM and rlen_ij < FF_2BODY[curr_pair_type_idx_ij].S_MAXIM)
			{		
				for(int a3=a2+1;a3<SYSTEM.ATOMS;a3++)
				{
					TEMP_STR = SYSTEM.ATOMTYPE[a1];
					TEMP_STR.append(SYSTEM.ATOMTYPE[a3]);	
					PAIR_TYPE_IK = TEMP_STR;							
					curr_pair_type_idx_ik = PAIR_MAP[TEMP_STR];	

					RVEC_IK.X = SYSTEM.COORDS[a3].X - SYSTEM.COORDS[a1].X;
					RVEC_IK.X -= floor( 0.5 + RVEC_IK.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;	
					
					RVEC_IK.Y = SYSTEM.COORDS[a3].Y - SYSTEM.COORDS[a1].Y;
					RVEC_IK.Y -= floor( 0.5 + RVEC_IK.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;	
					
					RVEC_IK.Z = SYSTEM.COORDS[a3].Z - SYSTEM.COORDS[a1].Z;
					RVEC_IK.Z -= floor( 0.5 + RVEC_IK.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;		
	
					RAB_IK.X = RVEC_IK.X;
					RAB_IK.Y = RVEC_IK.Y;
					RAB_IK.Z = RVEC_IK.Z;										
			
					rlen_ik = sqrt( RAB_IK.X*RAB_IK.X + RAB_IK.Y*RAB_IK.Y + RAB_IK.Z*RAB_IK.Z );	

					if(rlen_ik > FF_2BODY[curr_pair_type_idx_ik].S_MINIM and rlen_ik < FF_2BODY[curr_pair_type_idx_ik].S_MAXIM)
					{	
						
						TEMP_STR = SYSTEM.ATOMTYPE[a2];
						TEMP_STR.append(SYSTEM.ATOMTYPE[a3]);	
						PAIR_TYPE_JK = TEMP_STR;							
						curr_pair_type_idx_jk = PAIR_MAP[TEMP_STR];					
				
						RVEC_JK.X = SYSTEM.COORDS[a3].X - SYSTEM.COORDS[a2].X;
						RVEC_JK.X -= floor( 0.5 + RVEC_JK.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
						
						RVEC_JK.Y = SYSTEM.COORDS[a3].Y - SYSTEM.COORDS[a2].Y;
						RVEC_JK.Y -= floor( 0.5 + RVEC_JK.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
						
						RVEC_JK.Z = SYSTEM.COORDS[a3].Z - SYSTEM.COORDS[a2].Z;
						RVEC_JK.Z -= floor( 0.5 + RVEC_JK.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;	
		
						RAB_JK.X = RVEC_JK.X;
						RAB_JK.Y = RVEC_JK.Y;
						RAB_JK.Z = RVEC_JK.Z;										
				
						rlen_jk = sqrt( RAB_JK.X*RAB_JK.X + RAB_JK.Y*RAB_JK.Y + RAB_JK.Z*RAB_JK.Z );							
						
						if(rlen_jk > FF_2BODY[curr_pair_type_idx_jk].S_MINIM and rlen_jk < FF_2BODY[curr_pair_type_idx_jk].S_MAXIM)												
						{				
				
							// Everything is within allowed ranges. Begin setting up the derivative calculation
				
							// Set up the polynomials
							
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_ij], Tn_ij, Tnd_ij, rlen_ij, xdiff_ij);
							
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_ik], Tn_ik, Tnd_ik, rlen_ik, xdiff_ik);
							
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_jk], Tn_jk, Tnd_jk, rlen_jk, xdiff_jk);				
				
														
							// At this point we've completed all pre-calculations needed to populate the A matrix. Now we need to figure out 
							// where within the matrix to put the data, and to do so. 
							
							// Determine the FF type for the given triplet
							
							TEMP_STR =      FF_2BODY[curr_pair_type_idx_ij].PRPR_NM;
							TEMP_STR.append(FF_2BODY[curr_pair_type_idx_ik].PRPR_NM);	
							TEMP_STR.append(FF_2BODY[curr_pair_type_idx_jk].PRPR_NM);	
							
							curr_triple_type_index = TRIAD_MAP[TEMP_STR];

							// Note: This syntax is safe since there is only one possible SNUM_3B_CHEBY value for all interactions

							vstart = n_2b_cheby_terms;
							
							for (int i=0; i<curr_triple_type_index; i++)
								vstart += PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS;			


							if ( FF_2BODY[curr_pair_type_idx_ij].CHEBY_TYPE == "MORSE" ) 
							{
								// Set the pentaly function values
								
								fcut0_ij     = (1.0 - rlen_ij/FF_2BODY[curr_pair_type_idx_ij].S_MAXIM);
								fcut_ij      = fcut0_ij * fcut0_ij * fcut0_ij;
								fcutderiv_ij = -3.0 * fcut0_ij * fcut0_ij / FF_2BODY[curr_pair_type_idx_ij].S_MAXIM;		
								
								fcut0_ik     = (1.0 - rlen_ik/FF_2BODY[curr_pair_type_idx_ik].S_MAXIM);
								fcut_ik      = fcut0_ik * fcut0_ik * fcut0_ik;
								fcutderiv_ik = -3.0 * fcut0_ik * fcut0_ik / FF_2BODY[curr_pair_type_idx_ik].S_MAXIM;	
								
								fcut0_jk     = (1.0 - rlen_jk/FF_2BODY[curr_pair_type_idx_jk].S_MAXIM);
								fcut_jk      = fcut0_jk * fcut0_jk * fcut0_jk;
								fcutderiv_jk = -3.0 * fcut0_jk * fcut0_jk / FF_2BODY[curr_pair_type_idx_jk].S_MAXIM;		
								
								
								/////////////////////////////////////////////////////////////////////
								/////////////////////////////////////////////////////////////////////
								// Consider special restrictions on allowed triplet types and powers
								/////////////////////////////////////////////////////////////////////
								/////////////////////////////////////////////////////////////////////
								
								row_offset = 0;
								
								// --- THE KEY HERE IS TO UNDERSTAND THAT THE IJ, IK, AND JK HERE IS BASED ON ATOM PAIRS, AND DOESN'T NECESSARILY MATCH THE TRIPLET'S EXPECTED ORDER!
								
								for(int i=0; i<PAIR_TRIPLETS[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
								{
									row_offset = PAIR_TRIPLETS[curr_triple_type_index].PARAM_INDICIES[i];
									
									SET_3B_CHEBY_POWERS(FF_2BODY, PAIR_TRIPLETS[curr_triple_type_index], PAIR_MAP,  pow_ij, pow_ik, pow_jk, PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK, i);

									deriv_ij =  fcut_ij * Tnd_ij[pow_ij]
												* (-exp(-rlen_ij/FF_2BODY[curr_pair_type_idx_ij].LAMBDA)/FF_2BODY[curr_pair_type_idx_ij].LAMBDA)/xdiff_ij + fcutderiv_ij
												* Tn_ij [pow_ij];
									
									deriv_ik =  fcut_ik * Tnd_ik[pow_ik] 
												* (-exp(-rlen_ik/FF_2BODY[curr_pair_type_idx_ik].LAMBDA)/FF_2BODY[curr_pair_type_idx_ik].LAMBDA)/xdiff_ik + fcutderiv_ik
												* Tn_ik[pow_ik];
									
									deriv_jk =  fcut_jk * Tnd_jk[pow_jk]
												* (-exp(-rlen_jk/FF_2BODY[curr_pair_type_idx_jk].LAMBDA)/FF_2BODY[curr_pair_type_idx_jk].LAMBDA)/xdiff_jk + fcutderiv_jk
												* Tn_jk[pow_jk];	
																		
									// ij pairs
									
									FRAME_A_MATRIX[a1][vstart+row_offset].X += deriv_ij * fcut_ik * fcut_jk * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.X / rlen_ij;
									FRAME_A_MATRIX[a2][vstart+row_offset].X -= deriv_ij * fcut_ik * fcut_jk * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.X / rlen_ij;
							
									FRAME_A_MATRIX[a1][vstart+row_offset].Y += deriv_ij * fcut_ik * fcut_jk * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.Y / rlen_ij;
									FRAME_A_MATRIX[a2][vstart+row_offset].Y -= deriv_ij * fcut_ik * fcut_jk * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.Y / rlen_ij;
							
									FRAME_A_MATRIX[a1][vstart+row_offset].Z += deriv_ij * fcut_ik * fcut_jk * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.Z / rlen_ij;
									FRAME_A_MATRIX[a2][vstart+row_offset].Z -= deriv_ij * fcut_ik * fcut_jk * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.Z / rlen_ij;	
							
							
									// ik pairs
							
									FRAME_A_MATRIX[a1][vstart+row_offset].X += deriv_ik * fcut_ij * fcut_jk * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.X / rlen_ik;
									FRAME_A_MATRIX[a3][vstart+row_offset].X -= deriv_ik * fcut_ij * fcut_jk * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.X / rlen_ik;
							
									FRAME_A_MATRIX[a1][vstart+row_offset].Y += deriv_ik * fcut_ij * fcut_jk * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.Y / rlen_ik;
									FRAME_A_MATRIX[a3][vstart+row_offset].Y -= deriv_ik * fcut_ij * fcut_jk * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.Y / rlen_ik;
							
									FRAME_A_MATRIX[a1][vstart+row_offset].Z += deriv_ik * fcut_ij * fcut_jk * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.Z / rlen_ik;
									FRAME_A_MATRIX[a3][vstart+row_offset].Z -= deriv_ik * fcut_ij * fcut_jk * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.Z / rlen_ik;
									
									// jk pairs
									
									FRAME_A_MATRIX[a2][vstart+row_offset].X += deriv_jk * fcut_ij * fcut_ik * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.X / rlen_jk;
									FRAME_A_MATRIX[a3][vstart+row_offset].X -= deriv_jk * fcut_ij * fcut_ik * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.X / rlen_jk;
							
									FRAME_A_MATRIX[a2][vstart+row_offset].Y += deriv_jk * fcut_ij * fcut_ik * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.Y / rlen_jk;
									FRAME_A_MATRIX[a3][vstart+row_offset].Y -= deriv_jk * fcut_ij * fcut_ik * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.Y / rlen_jk;
							
									FRAME_A_MATRIX[a2][vstart+row_offset].Z += deriv_jk * fcut_ij * fcut_ik * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.Z / rlen_jk;
									FRAME_A_MATRIX[a3][vstart+row_offset].Z -= deriv_jk * fcut_ij * fcut_ik * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.Z / rlen_jk;
									
								}
							}
							else 
							{
								
								for ( int i = 0; i < FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY; i++ )
								{
									for ( int j = 0; j < FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY; j++ )
									{
										for ( int k = 0; k < FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY; k++ )
										{
											
											row_offset  = FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY*FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY*i; 
											row_offset += FF_2BODY[curr_pair_type_idx_ij].SNUM_3B_CHEBY*j;
											row_offset += k;											
											
											FRAME_A_MATRIX[a1][vstart+row_offset].X += (FF_2BODY[curr_pair_type_idx_ij].S_MAXIM - rlen_ij) * Tn_ij[i] * RAB_IJ.X/rlen3_ij;
											FRAME_A_MATRIX[a1][vstart+row_offset].Y += (FF_2BODY[curr_pair_type_idx_ij].S_MAXIM - rlen_ij) * Tn_ij[i] * RAB_IJ.Y/rlen3_ij;
											FRAME_A_MATRIX[a1][vstart+row_offset].Z += (FF_2BODY[curr_pair_type_idx_ij].S_MAXIM - rlen_ij) * Tn_ij[i] * RAB_IJ.Z/rlen3_ij;
									
											FRAME_A_MATRIX[a2][vstart+row_offset].X -= (FF_2BODY[curr_pair_type_idx_ij].S_MAXIM - rlen_ij) * Tn_ij[i] * RAB_IJ.X/rlen3_ij;
											FRAME_A_MATRIX[a2][vstart+row_offset].Y -= (FF_2BODY[curr_pair_type_idx_ij].S_MAXIM - rlen_ij) * Tn_ij[i] * RAB_IJ.Y/rlen3_ij;
											FRAME_A_MATRIX[a2][vstart+row_offset].Z -= (FF_2BODY[curr_pair_type_idx_ij].S_MAXIM - rlen_ij) * Tn_ij[i] * RAB_IJ.Z/rlen3_ij;		
											
											
											FRAME_A_MATRIX[a1][vstart+row_offset].X += (FF_2BODY[curr_pair_type_idx_ik].S_MAXIM - rlen_ik) * Tn_ik[i] * RAB_IK.X/rlen3_ik;
											FRAME_A_MATRIX[a1][vstart+row_offset].Y += (FF_2BODY[curr_pair_type_idx_ik].S_MAXIM - rlen_ik) * Tn_ik[i] * RAB_IK.Y/rlen3_ik;
											FRAME_A_MATRIX[a1][vstart+row_offset].Z += (FF_2BODY[curr_pair_type_idx_ik].S_MAXIM - rlen_ik) * Tn_ik[i] * RAB_IK.Z/rlen3_ik;
									
											FRAME_A_MATRIX[a2][vstart+row_offset].X -= (FF_2BODY[curr_pair_type_idx_ik].S_MAXIM - rlen_ik) * Tn_ik[i] * RAB_IK.X/rlen3_ik;
											FRAME_A_MATRIX[a2][vstart+row_offset].Y -= (FF_2BODY[curr_pair_type_idx_ik].S_MAXIM - rlen_ik) * Tn_ik[i] * RAB_IK.Y/rlen3_ik;
											FRAME_A_MATRIX[a2][vstart+row_offset].Z -= (FF_2BODY[curr_pair_type_idx_ik].S_MAXIM - rlen_ik) * Tn_ik[i] * RAB_IK.Z/rlen3_ik;	
											
											
											FRAME_A_MATRIX[a1][vstart+row_offset].X += (FF_2BODY[curr_pair_type_idx_jk].S_MAXIM - rlen_jk) * Tn_jk[i] * RAB_JK.X/rlen3_jk;
											FRAME_A_MATRIX[a1][vstart+row_offset].Y += (FF_2BODY[curr_pair_type_idx_jk].S_MAXIM - rlen_jk) * Tn_jk[i] * RAB_JK.Y/rlen3_jk;
											FRAME_A_MATRIX[a1][vstart+row_offset].Z += (FF_2BODY[curr_pair_type_idx_jk].S_MAXIM - rlen_jk) * Tn_jk[i] * RAB_JK.Z/rlen3_jk;
									
											FRAME_A_MATRIX[a2][vstart+row_offset].X -= (FF_2BODY[curr_pair_type_idx_jk].S_MAXIM - rlen_jk) * Tn_jk[i] * RAB_JK.X/rlen3_jk;
											FRAME_A_MATRIX[a2][vstart+row_offset].Y -= (FF_2BODY[curr_pair_type_idx_jk].S_MAXIM - rlen_jk) * Tn_jk[i] * RAB_JK.Y/rlen3_jk;
											FRAME_A_MATRIX[a2][vstart+row_offset].Z -= (FF_2BODY[curr_pair_type_idx_jk].S_MAXIM - rlen_jk) * Tn_jk[i] * RAB_JK.Z/rlen3_jk;
											
										}
									}
								}									

							}
														
						} // end if rlen_jk within cutoffs...
					} // end if rlen_ik within cutoffs...	
	
				} // end third loop over atoms			
			} // end if rlen_ik within cutoffs...
		}
	}
}

// FUNCTION UPDATED
static void ZCalc_InvR_Deriv(FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)		
// Calculate derivatives of the forces wrt to inverse pair distance to various powers. Stores minimum distance between a pair of atoms in minD[i].
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

	for(int a1=0;a1<SYSTEM.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<SYSTEM.ATOMS;a2++)
		{
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
			//calculate vstart:
		
			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
		   
			RVEC.X  = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
		   
			RVEC.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		   
			RVEC.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;
      
			for(int n1=-1*nlayers;n1<nlayers+1;n1++)
			{
				for(int n2=-1*nlayers;n2<nlayers+1;n2++)
				{
					for(int n3=-1*nlayers;n3<nlayers+1;n3++)
					{
						RAB.X = RVEC.X + n1 * SYSTEM.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * SYSTEM.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * SYSTEM.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );

						if ( rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST )	// spline term calculated w/cutoff:
							FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

						if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM && rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
						{
							x = rlen;
							
							// Calculate the penalty function (fc) and its derivative
							
							fc = (rlen - FF_2BODY[curr_pair_type_idx].S_MAXIM);
							dfc = fc;
							
							fc = fc*fc*fc;							
							dfc = 3*dfc*dfc;
						
							for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
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
static void ZCalc_Poly_Deriv(FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)	
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
	static double *rc;//static_cast<int>(FF_2BODY.size())];
	
	int dim;
	double rfac;

	if ( ! called_before ) 
	{
		called_before = true;
		dim = 0;
		
		rc = new double [FF_2BODY.size()];
		
		for ( int i=0; i<FF_2BODY.size(); i++ ) 
		{
			// Convert cutoff (S_MAXIM) from units of Angstrom to bohr (atomic units), (rc)
			
			rc[i] = FF_2BODY[i].S_MAXIM/autoang;
			cout << "	rc[" << i << "] = " << fixed << setprecision(3) << rc[i] << " bohr (" << fixed << setprecision(3) << FF_2BODY[i].S_MAXIM << " Angstr.)" << endl;
			
			if ( FF_2BODY[i].SNUM > dim ) 
				dim = FF_2BODY[i].SNUM;
		}
		dim++;
	
	}

	// main loop for ninth order polynomial terms:
	
	for(int a1=0;a1<SYSTEM.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<SYSTEM.ATOMS;a2++)
		{
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
		
			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
		   
			RVEC.X  = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
		   
			RVEC.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		   
			RVEC.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;


			for(int n1=-1*nlayers;n1<nlayers+1;n1++)
			{
				for(int n2=-1*nlayers;n2<nlayers+1;n2++)
				{
			    	for(int n3=-1*nlayers;n3<nlayers+1;n3++)
					{
						RAB.X = RVEC.X + n1 * SYSTEM.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * SYSTEM.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * SYSTEM.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );
						
						if ( rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST )	// spline term calculated w/cutoff:
							FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

						if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM && rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
						{
							// calculate binning, convert all distances to au from angstroms 
							x = rlen/autoang;

							for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
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
void SubtractCoordForces(FRAME & SYSTEM, bool calc_deriv, vector<XYZ> & P_OVER_FORCES,  vector<PAIRS> & FF_2BODY, map<string,int> & PAIR_MAP)
{
	// this function subtracts the ReaxFF over-coordination term to re-fit 
	// splines/charges iteratively for self-consistence.
	// If calc_deriv is true, the derivative of the force wrt the magnitude of the
	// 3-body interaction is placed in Fderiv.  Otherwise, the 3-body force is subtracted
	// from the total forces given in force.  
	 

	// Lucas paper draft parameters.	// Chenoweth values for p1, p2, r0, lambda6 p(ovun2)
	//									// pover = 50.0;	= FF_2BODY[j].OVRPRMS[0];
	// p1  	  =	-2.5881042987450e-01;	// -0.0657;			= FF_2BODY[j].OVRPRMS[1];
	// r0 	  =	 9.6000387075695e-01;	//  5.0451;			= FF_2BODY[j].OVRPRMS[2];
	// p2 	  =	 3.8995379237250e+00;	//  1.0165;			= FF_2BODY[j].OVRPRMS[3];
	// lambda6=	-8.9;  					// -3.6141;			= FF_2BODY[j].OVRPRMS[4];

	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double Vtot = 0;
	
	vector<XYZ> SForce(SYSTEM.ATOMS);	// WHAT IS THIS?
	vector<XYZ> dEover(SYSTEM.ATOMS);	// Over bonding energy derivative (force) w/r/t delta ("S")
	vector<XYZ> dFover(SYSTEM.ATOMS);	// Over bonding force derivative w/r/t pover
	
	// THREE-BODY POTENTIAL:
	
	double Eover = 0;			// Overbonding energy
	double p1,p2,r0,pover,lambda6;
	double rik,tempr,temps;	
	double S[SYSTEM.ATOMS];				// Difference between oxygen bond order and valency (Called delta in the paper)
	 
	int i_pair; 
	int curr_pair_type_idx; 
	string TEMP_STR;
	
	// Set up some variables
	
	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
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
			
	for(int ai=0;ai<SYSTEM.ATOMS;ai++) // Calculates Eover... I'm assuming this is primarily used for MD, which is why I don't have it
	{
		temps = 0.0;
	  
		for(int ak=0;ak<SYSTEM.ATOMS;ak++)
		{			
			if (ai != ak)
			{
				// TWO-BODY PART... ONLY CARES ABOUT LOOPS OVER AI AND AK
			
				TEMP_STR = SYSTEM.ATOMTYPE[ai];
				TEMP_STR.append(SYSTEM.ATOMTYPE[ak]);
				curr_pair_type_idx = PAIR_MAP[TEMP_STR];
				
				// THE ITEMS COMMENTED OUT IN THIS LOOP SHOULD BE PUT BACK ONCE CODE COMPARISON WITH LUCAS' COMPLETE!!!!!
					
				if(FF_2BODY[curr_pair_type_idx].USE_OVRPRMS && (SYSTEM.ATOMTYPE[ai] == FF_2BODY[curr_pair_type_idx].OVER_TO_ATM)) // Then we should have a defined pair type
				{
					// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
		   
					RVEC.X  = SYSTEM.COORDS[ak].X - SYSTEM.COORDS[ai].X;
					RVEC.X -= floor(0.5 + RVEC.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;

					RVEC.Y = SYSTEM.COORDS[ak].Y - SYSTEM.COORDS[ai].Y;
					RVEC.Y -= floor(0.5 +  RVEC.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		   
					RVEC.Z = SYSTEM.COORDS[ak].Z - SYSTEM.COORDS[ai].Z;
					RVEC.Z -= floor(0.5 +  RVEC.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;

				
					rik = sqrt( RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z);	
				
				// Calculate the O--H bond order as defined by ReaxFF
				// 
					temps +=  exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[2]
					*pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3])); 
				}
			}
		}
		
		S[ai]  = temps - 2.0; // This is the reaxff delta (diff between oxygen valency and bond order) 
		Eover += FF_2BODY[curr_pair_type_idx].OVRPRMS[0] * S[ai] * 1.0/(1.0+exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai]));	

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

	for(int ai=0; ai<SYSTEM.ATOMS; ai++)
	{

		for(int ak=0; ak<SYSTEM.ATOMS; ak++)
		{	
			TEMP_STR = SYSTEM.ATOMTYPE[ai];
			TEMP_STR.append(SYSTEM.ATOMTYPE[ak]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
			
			// If it's a pair type with defined over params
			// and the pair isn't comprised of the exact same atom
			// and ai is the "to" atom (i.e. O in an O--H pair)
		
			if((FF_2BODY[curr_pair_type_idx].USE_OVRPRMS && ai != ak && SYSTEM.ATOMTYPE[ai] == FF_2BODY[curr_pair_type_idx].OVER_TO_ATM))
			{

				RVEC.X  = SYSTEM.COORDS[ak].X - SYSTEM.COORDS[ai].X;
				RVEC.X -= floor(0.5 + RVEC.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
		   
				RVEC.Y = SYSTEM.COORDS[ak].Y - SYSTEM.COORDS[ai].Y;
				RVEC.Y -= floor(0.5 + RVEC.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		   
				RVEC.Z = SYSTEM.COORDS[ak].Z - SYSTEM.COORDS[ai].Z;
				RVEC.Z -= floor(0.5 + RVEC.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;
				
				rik = sqrt( RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z);				
				
				// Calculate the derivative of Eover w/r/t delta. Terms emerge from quotient rule combined with chain rule.

				tempr = 1.0/(1+exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai])) - FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai]*exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai])
					        / pow(1.0+exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai]),2);

				tempr *= FF_2BODY[curr_pair_type_idx].OVRPRMS[2]*pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3]) *
					     FF_2BODY[curr_pair_type_idx].OVRPRMS[3]*exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[2] *
						 pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3]))/rik;	


				dEover[ai].X -= tempr*RVEC.X/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				dEover[ai].Y -= tempr*RVEC.Y/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				dEover[ai].Z -= tempr*RVEC.Z/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];						

				dEover[ak].X += tempr*RVEC.X/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				dEover[ak].Y += tempr*RVEC.Y/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				dEover[ak].Z += tempr*RVEC.Z/rik*FF_2BODY[curr_pair_type_idx].OVRPRMS[0];	
				
				dFover[ai].X -= tempr*RVEC.X/rik;
				dFover[ai].Y -= tempr*RVEC.Y/rik;
				dFover[ai].Z -= tempr*RVEC.Z/rik;						

				dFover[ak].X += tempr*RVEC.X/rik;
				dFover[ak].Y += tempr*RVEC.Y/rik;
				dFover[ak].Z += tempr*RVEC.Z/rik;					


			}	

		}


	}

    for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
        SForce[a1].X -= dEover[a1].X;
		SForce[a1].Y -= dEover[a1].Y;
		SForce[a1].Z -= dEover[a1].Z;	
	}

    if ( !calc_deriv) // subtract CoorD[i]. force:
	{ 
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
	        SYSTEM.FORCES[a1].X -= SForce[a1].X;
			SYSTEM.FORCES[a1].Y -= SForce[a1].Y;
			SYSTEM.FORCES[a1].Z -= SForce[a1].Z;	
		}
    }
    else // Then we're fitting pover... Calculate derivative of 3-body force wrt pover.
	{
		for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		{
	        P_OVER_FORCES[a1].X -= dFover[a1].X;
			P_OVER_FORCES[a1].Y -= dFover[a1].Y;
			P_OVER_FORCES[a1].Z -= dFover[a1].Z;		
		}
	} // VERIFIED FROM HERE UP
}


////////////////////////////////////////////////////////////
//
// FUNCTIONS -- FORCE CALCULATION
//
////////////////////////////////////////////////////////////
 

void ZCalc(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP)
{
	{
		static bool called_before = false;
	
		for(int a=0;a<SYSTEM.ATOMS;a++)
		{
			SYSTEM.ACCEL[a].X = 0;
			SYSTEM.ACCEL[a].Y = 0;
			SYSTEM.ACCEL[a].Z = 0;
		}
	
		SYSTEM.TOT_POT_ENER = 0;
		SYSTEM.PRESSURE_XYZ = 0;

		if      ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" ) 
		{
			ZCalc_Cheby(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP);
			  
			if(FF_2BODY[0].SNUM_3B_CHEBY > 0)
				ZCalc_3B_Cheby(SYSTEM, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP);
		}
		
		else if ( FF_2BODY[0].PAIRTYP == "LJ" ) 
			ZCalc_Lj(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP);
		
		else if ( FF_2BODY[0].PAIRTYP == "SPLINE" ) 
			ZCalc_Spline(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP);	
		
		if ( CONTROLS.USE_COULOMB ) 
			ZCalc_Ewald(SYSTEM);
		

		if ( CONTROLS.USE_OVERCOORD ) 
	      ZCalcSR_Over(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP);

		// FUNCTIONS THAT NEED UPDATING:

		/* 

	
		else if ( pair_type == INVERSE_R ) 	
			ZCalc_SR_Analytic(Coord,Lbc, Latcons,nlayers,nat,smin,smax,snum, SForce,Vtot, Pxyz, params);
    
		else if ( pair_type == STILLINGER )  // SAVE THIS FOR SECOND TO LAST FOR SIMILAR REASONS  
			ZCalc_Stillinger(Coord,Lbc, Latcons,nlayers,nat,smax, SForce,Vtot,Pxyz);
     
	
		// WHAT ABOUT DFTBPOLY? ... DOES THAT JUST MEAN WE WERE FITTING STUFF TO USE IN DFTB?

		else 
			EXIT_MSG("Error: Unknown pair type", pair_type)
	*/
	
		SYSTEM.PRESSURE_XYZ /= 3.0 * SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;


	  return;
	} 
}

// UPDATED AND VERIFIED
static void ZCalc_Cheby(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP)
// Calculate short-range forces using a Chebyshev polynomial expansion. Can use morse variables similar to the work of Bowman.
{
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double x;
	double deriv;
	int vstart;
	static double *Tn, *Tnd;
	static bool called_before = false;
	
	double xavg, xdiff;
	double xmin; 
	double xmax; 
	double rlen3; 	
	double coeff;			  
	double fcut0; 
	double fcut; 
	double fcutderiv; 				
	double tmp_doub; 
	double Vpenalty;

	// A penalty function is added to the potential for r + penalty_dist < smin[ipair]
	const double penalty_scale = 1.0e8;
	const double penalty_dist = 0.01;

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;
		
		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			if (FF_2BODY[i].SNUM > dim ) 
				dim = FF_2BODY[i].SNUM;	 
		
		dim++;
		Tn   = new double [dim];
		Tnd  = new double [dim];
	}
  
	// Main loop for Chebyshev terms:
	
	string TEMP_STR;
	int curr_pair_type_idx;

	for(int a1=0;a1<SYSTEM.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		for(int a2=a1+1;a2<SYSTEM.ATOMS;a2++)
		{
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];

			// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
		   
			RVEC.X  = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC.X -= floor( 0.5 + RVEC.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
		   
			RVEC.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC.Y -= floor( 0.5 + RVEC.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		   
			RVEC.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			RVEC.Z -= floor( 0.5 + RVEC.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;
      
			for(int n1=-1*CONTROLS.N_LAYERS;n1<CONTROLS.N_LAYERS+1;n1++)
			{
				for(int n2=-1*CONTROLS.N_LAYERS;n2<CONTROLS.N_LAYERS+1;n2++)
				{
					for(int n3=-1*CONTROLS.N_LAYERS;n3<CONTROLS.N_LAYERS+1;n3++)
					{
						RAB.X = RVEC.X + n1 * SYSTEM.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * SYSTEM.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * SYSTEM.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );
		
						if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM and rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
						{
							double xavg, xdiff, rpenalty;

							// Apply a penalty for distances less than smin - penalty_dist.
							
							if ( rlen - penalty_dist < FF_2BODY[curr_pair_type_idx].S_MINIM ) 
							{
								rpenalty = FF_2BODY[curr_pair_type_idx].S_MINIM + penalty_dist - rlen;
							} 
							else 
							{
								rpenalty = 0.0;
							}
							
							
							// Chebyshev polynomials are only defined on the range [-1,1], so transorm the pair distance
							// in a way that allows it to fall along that range. Options are:
							//
							// x = 1/pair_dist				// Inverse r, 
							// x = exp(pair_dist/lambda)	// Morse-type
							// x = pair_dist				// default type
							// 
							// All types are normalized by s_min to s_max range to fall along [-1,1]
							
							if ( FF_2BODY[curr_pair_type_idx].CHEBY_TYPE == "INVRSE_R" ) 
							{
								xavg  =  0.5 * (1.0/FF_2BODY[curr_pair_type_idx].S_MINIM + 1.0/FF_2BODY[curr_pair_type_idx].S_MAXIM); // midpoint of possible pair distances in r^-1 space
								xdiff =  0.5 * (1.0/FF_2BODY[curr_pair_type_idx].S_MINIM - 1.0/FF_2BODY[curr_pair_type_idx].S_MAXIM); // width of possible pair distances in r^-1 space
								x     = (1.0/rlen-xavg) / xdiff;																	  // pair distances in r^-1 space, normalized to fit over [-1,1]
							} 
							else if ( FF_2BODY[curr_pair_type_idx].CHEBY_TYPE == "MORSE" ) 
							{
								xmin  = exp(-FF_2BODY[curr_pair_type_idx].S_MAXIM/FF_2BODY[curr_pair_type_idx].LAMBDA); 
								xmax  = exp(-FF_2BODY[curr_pair_type_idx].S_MINIM/FF_2BODY[curr_pair_type_idx].LAMBDA); 
								xavg  = 0.5 * (xmin + xmax);																// midpoint of possible pair distances in morse space
								xdiff = 0.5 * (xmax - xmin);																// width of possible pair distances in morse space
								x = (exp(-rlen/FF_2BODY[curr_pair_type_idx].LAMBDA)-xavg)/xdiff;							// pair distances in morse space, normalized to fit over [-1,1]
							}
							else if (FF_2BODY[curr_pair_type_idx].CHEBY_TYPE == "DEFAULT")
							{
								xavg  = 0.5 * (FF_2BODY[curr_pair_type_idx].S_MINIM + FF_2BODY[curr_pair_type_idx].S_MAXIM); // midpoint of possible pair distances
								xdiff = 0.5 * (FF_2BODY[curr_pair_type_idx].S_MAXIM - FF_2BODY[curr_pair_type_idx].S_MINIM); // width of possible pair distances
								x = (rlen-xavg) / xdiff;																	 // pair distances, normalized to fit over [-1,1]
							}
							else
							{
								cout << "ERROR: Undefined CHBTYPE: " << FF_2BODY[curr_pair_type_idx].CHEBY_TYPE << endl;
								cout << "       Excepted values are \"DEFAULT\", \"INVRSE_R\", or \"MORSE\". " << endl;
								cout << "       Check the parameter input file." << endl;
								exit(1);
							}

							// Make sure our newly transformed distance falls between the bound of -1 to 1 allowed to Cheby polynomials

							if ( x < -1.0 ) 
							{
								cout << "Warning: In 2B Cheby transformation, r < rmin " << TEMP_STR << endl;
								x = -1.0;
							}
							if ( x > 1.0 ) 
							{
								x = 1.0;
							}
							
							
							// Generate Chebyshev polynomials by recursion. 
							// 
							// What we're doing here. Want to fit using Cheby polynomials of the 1st kinD[i]. "T_n(x)."
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
						  
							for ( int i = 2; i <= FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
							{
								Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2];
								Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2];
							}
							
							// Now multiply by n to convert Tnd's to actual derivatives of Tn

							for ( int i = FF_2BODY[curr_pair_type_idx].SNUM; i >= 1; i-- ) 
								Tnd[i] = i * Tnd[i-1];

							Tnd[0] = 0.0;
							
							
							// Now compute the force/potential

							if ( FF_2BODY[curr_pair_type_idx].CHEBY_TYPE == "MORSE" ) 
							{
								fcut0 = (1.0 - rlen/FF_2BODY[curr_pair_type_idx].S_MAXIM);
								fcut = fcut0 * fcut0 * fcut0;
								fcutderiv = -3.0 * fcut0 * fcut0 / FF_2BODY[curr_pair_type_idx].S_MAXIM;
								
								
								for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
								{
									coeff                = FF_2BODY[curr_pair_type_idx].PARAMS[i]; // This is the Cheby FF paramfor the given power
									SYSTEM.TOT_POT_ENER += coeff * fcut * Tn[i+1];
									deriv                = (fcut * Tnd[i+1] *(-exp(-rlen/FF_2BODY[curr_pair_type_idx].LAMBDA)/FF_2BODY[curr_pair_type_idx].LAMBDA)/xdiff + fcutderiv * Tn[i+1]);
									SYSTEM.PRESSURE_XYZ -= coeff * deriv * rlen;

									SYSTEM.ACCEL[a1].X += coeff * deriv * RAB.X / rlen;
									SYSTEM.ACCEL[a1].Y += coeff * deriv * RAB.Y / rlen;
									SYSTEM.ACCEL[a1].Z += coeff * deriv * RAB.Z / rlen;
									
									SYSTEM.ACCEL[a2].X -= coeff * deriv * RAB.X / rlen;
									SYSTEM.ACCEL[a2].Y -= coeff * deriv * RAB.Y / rlen;
									SYSTEM.ACCEL[a2].Z -= coeff * deriv * RAB.Z / rlen;
	
								}
								// Add penalty for very short distances, where the fit FF may be unphysical (preserve conservation of E).

								if ( rpenalty > 0.0 ) 
								{
									Vpenalty = 0.0;
									cout << "Warning: Adding penalty in 2B Cheby calc, r < rmin " << rlen << " " << FF_2BODY[curr_pair_type_idx].S_MINIM << " " << TEMP_STR << endl;
									
									SYSTEM.ACCEL[a2].X += 3.0 * rpenalty * rpenalty * penalty_scale * RAB.X / rlen;
									SYSTEM.ACCEL[a2].Y += 3.0 * rpenalty * rpenalty * penalty_scale * RAB.Y / rlen;
									SYSTEM.ACCEL[a2].Z += 3.0 * rpenalty * rpenalty * penalty_scale * RAB.Z / rlen;
									
									SYSTEM.ACCEL[a1].X -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB.X / rlen;
									SYSTEM.ACCEL[a1].Y -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB.Y / rlen;
									SYSTEM.ACCEL[a1].Z -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB.Z / rlen;								
									
									Vpenalty = rpenalty * rpenalty * rpenalty * penalty_scale;
									SYSTEM.TOT_POT_ENER += Vpenalty;
									cout << "Penalty potential = "<< Vpenalty << endl;
								}
							}
							else 
							{
								rlen3 = rlen * rlen * rlen;
								
								for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++) 
								{
									coeff = FF_2BODY[curr_pair_type_idx].PARAMS[i]; // This is the Cheby FF paramfor the given power
									
									// POTENTIAL NOT YET IMPLEMENTED[i].
									SYSTEM.PRESSURE_XYZ += coeff * (FF_2BODY[curr_pair_type_idx].S_MAXIM-rlen) * Tn[i] / rlen;
	
									SYSTEM.ACCEL[a1].X += coeff * (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * SYSTEM.COORDS[a1].X / rlen3;
									SYSTEM.ACCEL[a1].Y += coeff * (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * SYSTEM.COORDS[a1].Y / rlen3;
									SYSTEM.ACCEL[a1].Z += coeff * (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * SYSTEM.COORDS[a1].Z / rlen3;
								
									SYSTEM.ACCEL[a2].X -= coeff * (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * SYSTEM.COORDS[a1].X / rlen3;
									SYSTEM.ACCEL[a2].Y -= coeff * (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * SYSTEM.COORDS[a1].Y / rlen3;
									SYSTEM.ACCEL[a2].Z -= coeff * (FF_2BODY[curr_pair_type_idx].S_MAXIM - rlen) * Tn[i] * SYSTEM.COORDS[a1].Z / rlen3;

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

// UPDATED AND VERIFIED
static void ZCalc_3B_Cheby(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP)
// Calculate 3-body short-range forces using a Chebyshev polynomial expansion.
// Can use morse variables similar to the work of Bowman.
{
	XYZ RVEC_IJ;
	XYZ RVEC_IK;
	XYZ RVEC_JK;
	 		
	XYZ RAB_IJ;
	XYZ RAB_IK;
	XYZ RAB_JK;
	
	#if FORCECHECK == 1	
		static vector<XYZ> FORCE_3B;	// Equivalent of f3b 	
		static ofstream FILE_FORCE_3B;
	#endif
		
	double rlen_ij,  rlen_ik,  rlen_jk;
	double rlen3_ij, rlen3_ik, rlen3_jk;
	
	static double *Tn_ij,  *Tn_ik,  *Tn_jk;
	static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
	static bool called_before = false;
	
	static int pow_ij, pow_ik, pow_jk;

	double rlen3;				  
	double fcut0, fcut0_ij, fcut0_ik, fcut0_jk; 
	double fcut,  fcut_ij,  fcut_ik,  fcut_jk; 			
	double deriv, deriv_ij, deriv_ik, deriv_jk;
	double        force_ij, force_ik, force_jk;
	double        xdiff_ij, xdiff_ik, xdiff_jk;
	double fcutderiv, fcutderiv_ij, fcutderiv_ik, fcutderiv_jk; 	
	double tmp_doub, tempx; 	
	
	static string TEMP_STR;
	static string PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK;
	static int curr_triple_type_index;
	static int curr_pair_type_idx_ij;
	static int curr_pair_type_idx_ik;
	static int curr_pair_type_idx_jk;
	
	double coeff;

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;

		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			if (FF_2BODY[i].SNUM_3B_CHEBY > dim ) 
				dim = FF_2BODY[i].SNUM_3B_CHEBY;	
		dim++;
	
		Tn_ij   = new double [dim];
		Tn_ik   = new double [dim];
		Tn_jk   = new double [dim];

		Tnd_ij  = new double [dim];
		Tnd_ik  = new double [dim];
		Tnd_jk  = new double [dim];
		
		#if FORCECHECK
		
			FORCE_3B.resize(SYSTEM.ATOMS);
		
			for( int i=0; i<SYSTEM.ATOMS; i++)
				FORCE_3B[i].X = FORCE_3B[i].Y = FORCE_3B[i].Z = 0;
		
			FILE_FORCE_3B.open("3b_results.dat");
		
		#endif 
		
	}
	
	tempx = 0;
  
	// Main loop for Chebyshev terms:

	for(int a1=0;a1<SYSTEM.ATOMS-1;a1++)		// Double sum over atom pairs
	{
		
		for(int a2=a1+1;a2<SYSTEM.ATOMS;a2++)
		{	
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);		
			PAIR_TYPE_IJ = TEMP_STR;					
			curr_pair_type_idx_ij = PAIR_MAP[TEMP_STR];
			
			// Start with minimum image convention.  Use layers to access larger distances if desireD[i].	
			// NOTE: LAYERS ARE DISABLED FOR 3B Cheby. Loops become unwieldly... 

			RVEC_IJ.X = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC_IJ.X -= floor( 0.5 + RVEC_IJ.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
			
			RVEC_IJ.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC_IJ.Y -= floor( 0.5 + RVEC_IJ.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
						
			RVEC_IJ.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			RVEC_IJ.Z -= floor( 0.5 + RVEC_IJ.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;			
			
			RAB_IJ.X = RVEC_IJ.X;
			RAB_IJ.Y = RVEC_IJ.Y;
			RAB_IJ.Z = RVEC_IJ.Z;

			rlen_ij = sqrt( RAB_IJ.X*RAB_IJ.X + RAB_IJ.Y*RAB_IJ.Y + RAB_IJ.Z*RAB_IJ.Z );			

			// Before doing any polynomial/coeff set up, make sure that all ij, ik, and jk distances are 
			// within the allowed range.
			
			if(rlen_ij > FF_2BODY[curr_pair_type_idx_ij].S_MINIM and rlen_ij < FF_2BODY[curr_pair_type_idx_ij].S_MAXIM)
			{
				for(int a3=a2+1;a3<SYSTEM.ATOMS;a3++)	// Move on to ik and jk pairs
				{
					TEMP_STR = SYSTEM.ATOMTYPE[a1];
					TEMP_STR.append(SYSTEM.ATOMTYPE[a3]);	
					PAIR_TYPE_IK = TEMP_STR;							
					curr_pair_type_idx_ik = PAIR_MAP[TEMP_STR];	

					RVEC_IK.X = SYSTEM.COORDS[a3].X - SYSTEM.COORDS[a1].X;
					RVEC_IK.X -= floor( 0.5 + RVEC_IK.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;	
					
					RVEC_IK.Y = SYSTEM.COORDS[a3].Y - SYSTEM.COORDS[a1].Y;
					RVEC_IK.Y -= floor( 0.5 + RVEC_IK.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;	
					
					RVEC_IK.Z = SYSTEM.COORDS[a3].Z - SYSTEM.COORDS[a1].Z;
					RVEC_IK.Z -= floor( 0.5 + RVEC_IK.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;		
	
					RAB_IK.X = RVEC_IK.X;
					RAB_IK.Y = RVEC_IK.Y;
					RAB_IK.Z = RVEC_IK.Z;										
			
					rlen_ik = sqrt( RAB_IK.X*RAB_IK.X + RAB_IK.Y*RAB_IK.Y + RAB_IK.Z*RAB_IK.Z );	

					if(rlen_ik > FF_2BODY[curr_pair_type_idx_ik].S_MINIM and rlen_ik < FF_2BODY[curr_pair_type_idx_ik].S_MAXIM)
					{
						TEMP_STR = SYSTEM.ATOMTYPE[a2];
						TEMP_STR.append(SYSTEM.ATOMTYPE[a3]);	
						PAIR_TYPE_JK = TEMP_STR;							
						curr_pair_type_idx_jk = PAIR_MAP[TEMP_STR];					
				
						RVEC_JK.X = SYSTEM.COORDS[a3].X - SYSTEM.COORDS[a2].X;
						RVEC_JK.X -= floor( 0.5 + RVEC_JK.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
						
						RVEC_JK.Y = SYSTEM.COORDS[a3].Y - SYSTEM.COORDS[a2].Y;
						RVEC_JK.Y -= floor( 0.5 + RVEC_JK.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
						
						RVEC_JK.Z = SYSTEM.COORDS[a3].Z - SYSTEM.COORDS[a2].Z;
						RVEC_JK.Z -= floor( 0.5 + RVEC_JK.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;	
		
						RAB_JK.X = RVEC_JK.X;
						RAB_JK.Y = RVEC_JK.Y;
						RAB_JK.Z = RVEC_JK.Z;										
				
						rlen_jk = sqrt( RAB_JK.X*RAB_JK.X + RAB_JK.Y*RAB_JK.Y + RAB_JK.Z*RAB_JK.Z );							
						
						if(rlen_jk > FF_2BODY[curr_pair_type_idx_jk].S_MINIM and rlen_jk < FF_2BODY[curr_pair_type_idx_jk].S_MAXIM)												
						{
							// Everything is within allowed ranges. Begin setting up the force calculation
				
							// Set up the polynomials
							
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_ij], Tn_ij, Tnd_ij, rlen_ij, xdiff_ij);
							
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_ik], Tn_ik, Tnd_ik, rlen_ik, xdiff_ik);
							
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_jk], Tn_jk, Tnd_jk, rlen_jk, xdiff_jk);
							
							// Apply the FF

							if ( FF_2BODY[curr_pair_type_idx_ij].CHEBY_TYPE == "MORSE" )
							{	
								// Set up the penalty functions
								// 
								fcut0_ij     = (1.0 - rlen_ij/FF_2BODY[curr_pair_type_idx_ij].S_MAXIM);
								fcut_ij      = fcut0_ij * fcut0_ij * fcut0_ij;
								fcutderiv_ij = -3.0 * fcut0_ij * fcut0_ij / FF_2BODY[curr_pair_type_idx_ij].S_MAXIM;	
								
								fcut0_ik     = (1.0 - rlen_ik/FF_2BODY[curr_pair_type_idx_ik].S_MAXIM);
								fcut_ik      = fcut0_ik * fcut0_ik * fcut0_ik;
								fcutderiv_ik = -3.0 * fcut0_ik * fcut0_ik / FF_2BODY[curr_pair_type_idx_ik].S_MAXIM;	
								
								fcut0_jk     = (1.0 - rlen_jk/FF_2BODY[curr_pair_type_idx_jk].S_MAXIM);
								fcut_jk      = fcut0_jk * fcut0_jk * fcut0_jk;
								fcutderiv_jk = -3.0 * fcut0_jk * fcut0_jk / FF_2BODY[curr_pair_type_idx_jk].S_MAXIM;	
								
								// Determine the FF type for the given triplet
							
								TEMP_STR =      FF_2BODY[curr_pair_type_idx_ij].PRPR_NM;
								TEMP_STR.append(FF_2BODY[curr_pair_type_idx_ik].PRPR_NM);	
								TEMP_STR.append(FF_2BODY[curr_pair_type_idx_jk].PRPR_NM);	

								curr_triple_type_index = TRIAD_MAP[TEMP_STR];
								
								// Error check: Certain triplets are impossible...
								
								if(curr_triple_type_index < 0)
								{
									cout << "ERROR: Impossible atom triplet found: " << TEMP_STR << endl;
									cout << "       Check the parameter file." << endl;
									exit(0);
								}
								
								// Now compute the forces for each set of allowed powers for pairs ij, ik, and jk		
								// Keep in mind that the order in which allowed powers are stored may not match the
								// ordering of pairs resulting from the present atom triplet. Thus, we need to order
								// the stored powers properly before applying the FF.
								
								for(int i=0; i<FF_3BODY[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
								{
									SET_3B_CHEBY_POWERS(FF_2BODY, FF_3BODY[curr_triple_type_index], PAIR_MAP,  pow_ij, pow_ik, pow_jk, PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK, i);
									
									coeff = FF_3BODY[curr_triple_type_index].PARAMS[i];
									
									tempx += coeff * fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk];	
									
									deriv_ij = fcut_ij * Tnd_ij[pow_ij]
										     * (-exp(-rlen_ij/FF_2BODY[curr_pair_type_idx_ij].LAMBDA)/FF_2BODY[curr_pair_type_idx_ij].LAMBDA)/xdiff_ij + fcutderiv_ij
										     * Tn_ij [pow_ij];
									
									deriv_ik = fcut_ik * Tnd_ik[pow_ik]
										     * (-exp(-rlen_ik/FF_2BODY[curr_pair_type_idx_ik].LAMBDA)/FF_2BODY[curr_pair_type_idx_ik].LAMBDA)/xdiff_ik + fcutderiv_ik
										     * Tn_ik [pow_ik];
									
									deriv_jk = fcut_jk * Tnd_jk[pow_jk]
										     * (-exp(-rlen_jk/FF_2BODY[curr_pair_type_idx_jk].LAMBDA)/FF_2BODY[curr_pair_type_idx_jk].LAMBDA)/xdiff_jk + fcutderiv_jk
										     * Tn_jk [pow_jk];
									
									force_ij  = coeff * deriv_ij * fcut_ik * fcut_jk * Tn_ik [pow_ik] * Tn_jk [pow_jk];
									SYSTEM.PRESSURE_XYZ    -= force_ij * rlen_ij;
									force_ij /= rlen_ij;
									
									force_ik  = coeff * deriv_ik * fcut_ij * fcut_jk * Tn_ij [pow_ij] * Tn_jk [pow_jk];
									SYSTEM.PRESSURE_XYZ    -= force_ik * rlen_ik;
									force_ik /= rlen_ik;
									
									force_jk  = coeff * deriv_jk * fcut_ij * fcut_ik * Tn_ij [pow_ij] * Tn_ik [pow_ik];
									SYSTEM.PRESSURE_XYZ    -= force_jk * rlen_jk;
									force_jk /= rlen_jk;
								
									// Apply forces to ij pair
									
									SYSTEM.ACCEL[a1].X += force_ij * RAB_IJ.X;
									SYSTEM.ACCEL[a1].Y += force_ij * RAB_IJ.Y;
									SYSTEM.ACCEL[a1].Z += force_ij * RAB_IJ.Z;
									
									SYSTEM.ACCEL[a2].X -= force_ij * RAB_IJ.X;
									SYSTEM.ACCEL[a2].Y -= force_ij * RAB_IJ.Y;
									SYSTEM.ACCEL[a2].Z -= force_ij * RAB_IJ.Z;
									
									// Apply forces to ik pair
									
									SYSTEM.ACCEL[a1].X += force_ik * RAB_IK.X;
									SYSTEM.ACCEL[a1].Y += force_ik * RAB_IK.Y;
									SYSTEM.ACCEL[a1].Z += force_ik * RAB_IK.Z;	
									
									SYSTEM.ACCEL[a3].X -= force_ik * RAB_IK.X;
									SYSTEM.ACCEL[a3].Y -= force_ik * RAB_IK.Y;
									SYSTEM.ACCEL[a3].Z -= force_ik * RAB_IK.Z;	
									
									// Apply forces to jk pair
									
									SYSTEM.ACCEL[a2].X += force_jk * RAB_JK.X;
									SYSTEM.ACCEL[a2].Y += force_jk * RAB_JK.Y;
									SYSTEM.ACCEL[a2].Z += force_jk * RAB_JK.Z;
									
									SYSTEM.ACCEL[a3].X -= force_jk * RAB_JK.X;
									SYSTEM.ACCEL[a3].Y -= force_jk * RAB_JK.Y;
									SYSTEM.ACCEL[a3].Z -= force_jk * RAB_JK.Z;										
									
									#if FORCECHECK
									
										// Apply forces to ij pair
									
										FORCE_3B[a1].X += force_ij * RAB_IJ.X;
										FORCE_3B[a1].Y += force_ij * RAB_IJ.Y;
										FORCE_3B[a1].Z += force_ij * RAB_IJ.Z;
									
										FORCE_3B[a2].X -= force_ij * RAB_IJ.X;
										FORCE_3B[a2].Y -= force_ij * RAB_IJ.Y;
										FORCE_3B[a2].Z -= force_ij * RAB_IJ.Z;
									
										// Apply forces to ik pair
									
										FORCE_3B[a1].X += force_ik * RAB_IK.X;
										FORCE_3B[a1].Y += force_ik * RAB_IK.Y;
										FORCE_3B[a1].Z += force_ik * RAB_IK.Z;	
									
										FORCE_3B[a3].X -= force_ik * RAB_IK.X;
										FORCE_3B[a3].Y -= force_ik * RAB_IK.Y;
										FORCE_3B[a3].Z -= force_ik * RAB_IK.Z;	
									
										// Apply forces to jk pair
									
										FORCE_3B[a2].X += force_jk * RAB_JK.X;
										FORCE_3B[a2].Y += force_jk * RAB_JK.Y;
										FORCE_3B[a2].Z += force_jk * RAB_JK.Z;
									
										FORCE_3B[a3].X -= force_jk * RAB_JK.X;
										FORCE_3B[a3].Y -= force_jk * RAB_JK.Y;
										FORCE_3B[a3].Z -= force_jk * RAB_JK.Z;										
									
									#endif
									
																									
								}	
									
							}
							else
							{
								cout << "ERROR: Functionality not complete for non-Morse-type 3B Chebyshev potentials." << endl;
								exit(0);
							}						
						}	
					}				
				}
			}
		}
    }
	
	#if FORCECHECK
	
		FILE_FORCE_3B << "e3b = " << left << fixed << setprecision(16) << setw(16) << tempx << endl;
	
		for(int i=0; i<SYSTEM.ATOMS; i++)
		{
			FILE_FORCE_3B << left << fixed << setprecision(16) << setw(16) << FORCE_3B[i].X << endl;
			FILE_FORCE_3B << left << fixed << setprecision(16) << setw(16) << FORCE_3B[i].Y << endl;
			FILE_FORCE_3B << left << fixed << setprecision(16) << setw(16) << FORCE_3B[i].Z << endl;
		}
	
	#endif 
	
	SYSTEM.TOT_POT_ENER += tempx;

  return;
}

// UPDATED AND VERIFIED
static void ZCalc_Lj(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP)
// Calculate LJ interaction.. first parameter is epsilon, second parameter is sigma.
{
	// For the LJ potential, the first param for a given pair type is epsilon, and the second, sigma
	// ...eventually SMAX should be used for the pair distance cutoff value...
	
	XYZ		RVEC, RAB; 
	double	rlen_mi;
	int		curr_pair_type_idx;
	double	fac;
	string	TEMP_STR;

	for(int a1=0;a1<SYSTEM.ATOMS;a1++) 
	{
		for(int a2=0;a2<a1;a2++)
		{			
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];

			RVEC.X = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			
			RAB.X = RVEC.X - floor(0.5 + RVEC.X/SYSTEM.BOXDIM.X) * SYSTEM.BOXDIM.X;
			RAB.Y = RVEC.Y - floor(0.5 + RVEC.Y/SYSTEM.BOXDIM.Y) * SYSTEM.BOXDIM.Y;
			RAB.Z = RVEC.Z - floor(0.5 + RVEC.Z/SYSTEM.BOXDIM.Z) * SYSTEM.BOXDIM.Z;

			RVEC.X = RAB.X;
			RVEC.Y = RAB.Y;
			RVEC.Z = RAB.Z;
			
	  	  	rlen_mi = sqrt( RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z );
			
			if ( rlen_mi < 0.5 ) // FF_2BODY[curr_pair_type_idx].PARAMS[1]/2.2)
				EXIT_MSG("Error: close approach", rlen_mi);

			SYSTEM.TOT_POT_ENER += 4.0 * FF_2BODY[curr_pair_type_idx].PARAMS[0] * ( pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,12.0) - pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,6.0) );

			fac = 4.0 * FF_2BODY[curr_pair_type_idx].PARAMS[0] * ( 
				    -12.0 * pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi,14.0) 
				    + 6.0 * pow(FF_2BODY[curr_pair_type_idx].PARAMS[1]/rlen_mi, 8.0) );
			
			fac *= 1.0 / ( FF_2BODY[curr_pair_type_idx].PARAMS[1] * FF_2BODY[curr_pair_type_idx].PARAMS[1] );		
	  
			SYSTEM.PRESSURE_XYZ -= fac * (rlen_mi*rlen_mi);
			
			SYSTEM.ACCEL[a1].X += RVEC.X*fac;
			SYSTEM.ACCEL[a1].Y += RVEC.Y*fac;
			SYSTEM.ACCEL[a1].Z += RVEC.Z*fac;
			
			SYSTEM.ACCEL[a2].X -= RVEC.X*fac;
			SYSTEM.ACCEL[a2].Y -= RVEC.Y*fac;
			SYSTEM.ACCEL[a2].Z -= RVEC.Z*fac;
		}
	}
}

// UPDATED AND VERIFIED
static void ZCalc_Spline(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP)
// Calculate spline forces.
{
	XYZ RVEC, RAB;
	double rlen,rlen2;
	double tempx;
	double S_r;
	int vstart;
	
	int k0;
	double x, x0;
	double t, t2, t3, t4;;
	double h00,h10,h01,h11;
	double i00,i10,i01,i11;
	int kstart;
	
	string TEMP_STR;
	int curr_pair_type_idx;

	for(int a1=0;a1<SYSTEM.ATOMS-1;a1++) 
	{
		for(int a2=a1+1;a2<SYSTEM.ATOMS;a2++)
		{			
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
			
			RVEC.X = SYSTEM.COORDS[a2].X - SYSTEM.COORDS[a1].X;
			RVEC.Y = SYSTEM.COORDS[a2].Y - SYSTEM.COORDS[a1].Y;
			RVEC.Z = SYSTEM.COORDS[a2].Z - SYSTEM.COORDS[a1].Z;
			
			RAB.X = RVEC.X - floor(0.5 + RVEC.X/SYSTEM.BOXDIM.X) * SYSTEM.BOXDIM.X;
			RAB.Y = RVEC.Y - floor(0.5 + RVEC.Y/SYSTEM.BOXDIM.Y) * SYSTEM.BOXDIM.Y;
			RAB.Z = RVEC.Z - floor(0.5 + RVEC.Z/SYSTEM.BOXDIM.Z) * SYSTEM.BOXDIM.Z;
			
			for(int n1=-1*CONTROLS.N_LAYERS;n1<CONTROLS.N_LAYERS+1;n1++)
			{
				for(int n2=-1*CONTROLS.N_LAYERS;n2<CONTROLS.N_LAYERS+1;n2++)
				{
					for(int n3=-1*CONTROLS.N_LAYERS;n3<CONTROLS.N_LAYERS+1;n3++)
					{
						RAB.X = RVEC.X + n1 * SYSTEM.BOXDIM.X;
						RAB.Y = RVEC.Y + n2 * SYSTEM.BOXDIM.Y;
						RAB.Z = RVEC.Z + n3 * SYSTEM.BOXDIM.Z;

						rlen = sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );
		
						if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM && rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
						{			

							tempx = 0;
							rlen2 = rlen;		

						  	k0  = int(floor((rlen2-FF_2BODY[curr_pair_type_idx].S_MINIM)/FF_2BODY[curr_pair_type_idx].S_DELTA));
						  	x   = rlen2;
						  	x0  = FF_2BODY[curr_pair_type_idx].S_MINIM+FF_2BODY[curr_pair_type_idx].S_DELTA*k0;
						  	t   = (x-x0)/FF_2BODY[curr_pair_type_idx].S_DELTA;

						  	if ( t > 1.000001 || t < -0.000001 ) 
							  	EXIT_MSG("Error: Bad t in spline calculation.");

							t2 = t * t;
							t3 = t2 * t;
							t4 = t3 * t;

							h00 =  2*t3 -3*t2 + 1.0;
							h10 =    t3 -2*t2 + t;
							h01 = -2*t3 +3*t2;
							h11 =    t3 -  t2;
							h10 *= FF_2BODY[curr_pair_type_idx].S_DELTA; 		//derivative terms have extra factor.
							h11 *= FF_2BODY[curr_pair_type_idx].S_DELTA;
							
							kstart = k0*2;	
						  	    
							if ( kstart > FF_2BODY[curr_pair_type_idx].SNUM ) 
							{
								cout << "Error: kstart too large " << kstart << " " << rlen << " " << FF_2BODY[curr_pair_type_idx].S_MINIM << " " << FF_2BODY[curr_pair_type_idx].SNUM << endl;
								exit(1);
							}
						  
							S_r =
								h00 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+0]+//lookout!!!
								h10 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+1]+
								h01 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+2]+
								h11 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+3];

							//do integrating here to get potential energy:
							i00 = 0.5 - (0.5*t4 - t3 + t);
							i10 = 1.0/12.0 - (0.25*t4- (2.0/3.0)*t3 + 0.5*t2);
							i01 = 0.5 - (-0.5*t4 + t3);
							i11 = -1.0/12.0 - (0.25*t4 - t3/3.0);
							i10 *= FF_2BODY[curr_pair_type_idx].S_DELTA;
							i11 *= FF_2BODY[curr_pair_type_idx].S_DELTA;

							// TEST: Linear interpolation.
							// i00 = (0.5 - (t - t * t / 2.0));
							// i01 = (0.5 - t * t / 2.0);
							// i10 = 0.0;
							// i11 = 0.0;

							tempx =  -FF_2BODY[curr_pair_type_idx].S_DELTA * (
								i00 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+0]+
								i10 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+1]+
								i01 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+2]+
								i11 * FF_2BODY[curr_pair_type_idx].PARAMS[kstart+3] );
						
							// Use the integral of the spline for pressure calculation
							
							if ( kstart/2 + 1 < FF_2BODY[curr_pair_type_idx].SNUM / 2 )
								tempx += FF_2BODY[curr_pair_type_idx].POT_PARAMS[kstart/2+1];

						SYSTEM.ACCEL[a1].X += S_r * RAB.X/rlen;
						SYSTEM.ACCEL[a1].Y += S_r * RAB.Y/rlen;
						SYSTEM.ACCEL[a1].Z += S_r * RAB.Z/rlen;
		
						SYSTEM.ACCEL[a2].X -= S_r * RAB.X/rlen;
						SYSTEM.ACCEL[a2].Y -= S_r * RAB.Y/rlen;
						SYSTEM.ACCEL[a2].Z -= S_r * RAB.Z/rlen;
				   
						SYSTEM.PRESSURE_XYZ -= S_r * rlen;
						SYSTEM.TOT_POT_ENER += tempx;

						}//rlen
					}
				}
			}
		}
	}
  return;
}

// UPDATED AND VERIFIED
static void ZCalcSR_Over(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP)
// Calculate short-range overcoordination forces... See SubtractCoordForces for more info
{
	
	

	
	
	// Notes:
	//
  	// these are short-ranged cutoffs for analytical potential.
  	// e.g. for r<oo_cut, force(r)=force(oo_cut). The values are 
  	// at least 0.2 Angs below any sampled by MD; they are there 
  	// in case trajectories reach them.
	//
  	// Many-body overcoordination term:
  	// The rest of this function can be commented out and the MD 
  	// will still run. The ReaxFF terms have been fitted using
  	// Powell minimization (separate routines, but very easy to do).
	
	XYZ RVEC, RAB;

	// MANY-BODY POTENTIAL:
	double Eover = 0;
	double pover,p1,p2,r0,lambda6;
	double rik,tempr,temps;
	double powrik;
	double Sexp2;
	
	vector<double> S     (SYSTEM.ATOMS);
	vector<double> Sexp  (SYSTEM.ATOMS);
	vector<XYZ>    dEover(SYSTEM.ATOMS);
	
	int curr_pair_type_idx; 
	string TEMP_STR;
	
	bool SAFE = false;

	for(int ai=0;ai<SYSTEM.ATOMS;ai++)
	{
		temps=0.0;
		S[ai] = 0.0;
		Sexp[ai] = 0.0;

		for(int i=0; i<FF_2BODY.size(); i++)
			if(SYSTEM.ATOMTYPE[ai] == FF_2BODY[i].OVER_TO_ATM)
				SAFE = true;
		
		if(SAFE)	
		{
			for(int ak=0;ak<SYSTEM.ATOMS;ak++)
			{
				if(ai!=ak)
				{			
					TEMP_STR = SYSTEM.ATOMTYPE[ai];
					TEMP_STR.append(SYSTEM.ATOMTYPE[ak]);
					curr_pair_type_idx = PAIR_MAP[TEMP_STR];	
		
					if(FF_2BODY[curr_pair_type_idx].USE_OVRPRMS)  // Then we should have a defined pair type
					{									
						RVEC.X = SYSTEM.COORDS[ak].X - SYSTEM.COORDS[ai].X;
						RVEC.Y = SYSTEM.COORDS[ak].Y - SYSTEM.COORDS[ai].Y;
						RVEC.Z = SYSTEM.COORDS[ak].Z - SYSTEM.COORDS[ai].Z;

						// Short-range interaction, so use minimum image convention.

						RAB.X = RVEC.X - floor(0.5+RVEC.X/SYSTEM.BOXDIM.X)*SYSTEM.BOXDIM.X;
						RAB.Y = RVEC.Y - floor(0.5+RVEC.Y/SYSTEM.BOXDIM.Y)*SYSTEM.BOXDIM.Y;
						RAB.Z = RVEC.Z - floor(0.5+RVEC.Z/SYSTEM.BOXDIM.Z)*SYSTEM.BOXDIM.Z;

						rik    = sqrt(RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z);
						temps += exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[2]
						       * pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3])); 
					}				
				}
			}	
			S[ai]    = temps-2.0;
			Sexp[ai] = exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[4]*S[ai]);
			Eover   += FF_2BODY[curr_pair_type_idx].OVRPRMS[0] * S[ai] * 1.0/(1.0+Sexp[ai]);			
		}
		
		SAFE = false;
	}
	
	SYSTEM.TOT_POT_ENER += Eover;
	
	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
		dEover[a1].X = dEover[a1].Y = dEover[a1].Z = 0;


	for(int ai=0; ai<SYSTEM.ATOMS; ai++)
	{
		for(int ak=0; ak<SYSTEM.ATOMS; ak++)
		{	
			TEMP_STR = SYSTEM.ATOMTYPE[ai];
			TEMP_STR.append(SYSTEM.ATOMTYPE[ak]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];

			if((FF_2BODY[curr_pair_type_idx].USE_OVRPRMS && ai != ak && SYSTEM.ATOMTYPE[ai] == FF_2BODY[curr_pair_type_idx].OVER_TO_ATM))
			{

				RVEC.X  = SYSTEM.COORDS[ak].X - SYSTEM.COORDS[ai].X;
				RVEC.X -= floor(0.5 + RVEC.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
	   
				RVEC.Y = SYSTEM.COORDS[ak].Y - SYSTEM.COORDS[ai].Y;
				RVEC.Y -= floor(0.5 + RVEC.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
	   
				RVEC.Z = SYSTEM.COORDS[ak].Z - SYSTEM.COORDS[ai].Z;
				RVEC.Z -= floor(0.5 + RVEC.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;
			
				rik = sqrt( RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z);				
			
				// Calculate the derivative of Eover w/r/t delta. Terms emerge from quotient rule combined with chain rule.

				powrik = pow(rik/FF_2BODY[curr_pair_type_idx].OVRPRMS[1],FF_2BODY[curr_pair_type_idx].OVRPRMS[3]);
				Sexp2  = (1.0 + Sexp[ai]) * (1.0 + Sexp[ai]);

				tempr  = FF_2BODY[curr_pair_type_idx].OVRPRMS[0];
				tempr *= 1.0/(1+Sexp[ai]) - FF_2BODY[curr_pair_type_idx].OVRPRMS[4] * S[ai] * Sexp[ai] / Sexp2;
				tempr *= FF_2BODY[curr_pair_type_idx].OVRPRMS[2] * powrik * FF_2BODY[curr_pair_type_idx].OVRPRMS[3] * exp(FF_2BODY[curr_pair_type_idx].OVRPRMS[2] * powrik) / rik;				

				dEover[ai].X -= tempr*RVEC.X/rik;
				dEover[ai].Y -= tempr*RVEC.Y/rik;
				dEover[ai].Z -= tempr*RVEC.Z/rik;						

				
				dEover[ak].X += tempr*RVEC.X/rik;
				dEover[ak].Y += tempr*RVEC.Y/rik;
				dEover[ak].Z += tempr*RVEC.Z/rik;		
				
				SYSTEM.PRESSURE_XYZ -= tempr * rik;							
			}
		}
	}

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)
	{
		SYSTEM.ACCEL[a1].X -= dEover[a1].X;
		SYSTEM.ACCEL[a1].Y -= dEover[a1].Y;
		SYSTEM.ACCEL[a1].Z -= dEover[a1].Z;
	}
		
}


////////////////////////////////////////////////////////////
//
// FUNCTIONS -- PRINTING OF POTENTIAL ENERGY SURFACE
//
////////////////////////////////////////////////////////////

void Print_Cheby(vector<PAIR_FF> & FF_2BODY, int ij, string PAIR_NAME, string FILE_TAG)
// Generating pair distance scans for the 2-b potential.
// pair distances will range from smin to smax, incremented by sdelta
{
	
//	string SCAN_FILE_2B;	// Declared as a global (external) variable in functions.h
	
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double x;
	double deriv;
	int vstart;
	static double *Tn, *Tnd;
	static bool called_before = false;
	
	double xavg, xdiff,rpenalty;
	double xmin; 
	double xmax; 
	double rlen3; 	
	double coeff;			  
	double fcut0; 
	double fcut; 
	double fcutderiv; 				
	double tmp_doub; 
	double Vpenalty;
	
	
	string OUTFILE = "2b_Cheby_Pot-";
	OUTFILE.append(PAIR_NAME);
	
	if(FILE_TAG != "")
		OUTFILE.append("_for_3B");
	
	OUTFILE.append(".dat");
	ofstream OUTFILE_2B_POT;
	OUTFILE_2B_POT.open(OUTFILE.data());	
	
	SCAN_FILE_2B = OUTFILE;
	

	// A penalty function is added to the potential for r + penalty_dist < smin[ipair]
	const double penalty_scale = 1.0e8;
	const double penalty_dist = 0.01;

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;
		
		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			if (FF_2BODY[i].SNUM > dim ) 
				dim = FF_2BODY[i].SNUM;	 
		
		dim++;
		Tn   = new double [dim];
		Tnd  = new double [dim];
	}
  
	// Main loop for Chebyshev terms:
	
	int n_ij = (FF_2BODY[ij].S_MAXIM - FF_2BODY[ij].S_MINIM)/FF_2BODY[ij].S_DELTA;
	
	double tempx;

	for (int a=1; a<n_ij; a++)
	{
		tempx = 0;

		rlen = FF_2BODY[ij].S_MINIM + a * FF_2BODY[ij].S_DELTA;

		if(rlen > FF_2BODY[ij].S_MINIM and rlen < FF_2BODY[ij].S_MAXIM)
		{
			// Apply a penalty for distances less than smin - penalty_dist.
			
			rpenalty = 0.0;
			
			if ( rlen - penalty_dist < FF_2BODY[ij].S_MINIM ) 
				rpenalty = FF_2BODY[ij].S_MINIM + penalty_dist - rlen;

			// Chebyshev polynomials are only defined on the range [-1,1], so transorm the pair distance
			// in a way that allows it to fall along that range. Options are:
			//
			// x = 1/pair_dist				// Inverse r, 
			// x = exp(pair_dist/lambda)	// Morse-type
			// x = pair_dist				// default type
			// 
			// All types are normalized by s_min to s_max range to fall along [-1,1]
			
			if ( FF_2BODY[ij].CHEBY_TYPE == "INVRSE_R" ) 
			{
				cout << "Error in Print_2B_Cheby: Functionality not yet programmed." << endl;
				exit(0);
												  // pair distances in r^-1 space, normalized to fit over [-1,1]
			} 
			else if ( FF_2BODY[ij].CHEBY_TYPE == "MORSE" ) 
			{
				xmin  = exp(-FF_2BODY[ij].S_MAXIM/FF_2BODY[ij].LAMBDA); 
				xmax  = exp(-FF_2BODY[ij].S_MINIM/FF_2BODY[ij].LAMBDA); 
				xavg  = 0.5 * (xmin + xmax);												// midpoint of possible pair distances in morse space
				xdiff = 0.5 * (xmax - xmin);												// width of possible pair distances in morse space
				x = (exp(-rlen/FF_2BODY[ij].LAMBDA)-xavg)/xdiff;							// pair distances in morse space, normalized to fit over [-1,1]
			}
			else if (FF_2BODY[ij].CHEBY_TYPE == "DEFAULT")
			{
				cout << "Error in Print_2B_Cheby: Functionality not yet programmed." << endl;
				exit(0);
												 // pair distances, normalized to fit over [-1,1]
			}
			else
			{
				cout << "ERROR: Undefined CHBTYPE: " << FF_2BODY[ij].CHEBY_TYPE << endl;
				cout << "       Excepted values are \"DEFAULT\", \"INVRSE_R\", or \"MORSE\". " << endl;
				cout << "       Check the parameter input file." << endl;
				exit(1);
			}

			// Make sure our newly transformed distance falls between the bound of -1 to 1 allowed to Cheby polynomials

			if ( x < -1.0 ) 
			{
				cout << "Warning:  r < rmin\n";
				x = -1.0;
			}
			if ( x > 1.0 ) 
			{
				x = 1.0;
			}
			
			
			// Generate Chebyshev polynomials by recursion. 
			// 
			// What we're doing here. Want to fit using Cheby polynomials of the 1st kinD[i]. "T_n(x)."
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
		  
			for ( int i = 2; i <= FF_2BODY[ij].SNUM; i++ ) 
			{
				Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2];
				Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2];
			}
			
			// Now multiply by n to convert Tnd's to actual derivatives of Tn

			for ( int i = FF_2BODY[ij].SNUM; i >= 1; i-- ) 
				Tnd[i] = i * Tnd[i-1];

			Tnd[0] = 0.0;
			
			
			// Now compute the force/potential

			if ( FF_2BODY[ij].CHEBY_TYPE == "MORSE" ) 
			{
				fcut0 = (1.0 - rlen/FF_2BODY[ij].S_MAXIM);
				fcut = fcut0 * fcut0 * fcut0;
				fcutderiv = -3.0 * fcut0 * fcut0 / FF_2BODY[ij].S_MAXIM;
				
				
				for ( int i = 0; i < FF_2BODY[ij].SNUM; i++ ) 
				{
					coeff                = FF_2BODY[ij].PARAMS[i]; // This is the Cheby FF paramfor the given power
					tempx += coeff * fcut * Tn[i+1];
				}
				// Add penalty for very short distances, where the fit FF may be unphysical (preserve conservation of E).

				if ( rpenalty > 0.0 ) 
				{
					Vpenalty = 0.0;
					Vpenalty = rpenalty * rpenalty * rpenalty * penalty_scale;
					tempx += Vpenalty;
				}
				
				OUTFILE_2B_POT << rlen << " " << tempx << endl;
			}
			else 
			{
				cout << "Error in Print_2B_Cheby: Functionality not yet programmed." << endl;
				exit(0);
				rlen3 = rlen * rlen * rlen;
			}							
		} 

	}
	
	OUTFILE_2B_POT.close();
	
	return;
}  

void Print_3B_Cheby(MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk)
// Generating 3d heatmaps to draw a 3-b potential.
// pair distances will range from smin to smax, incremented by sdelta
{
//	string FULL_FILE_3B;	// Declared as a global (external) variable in functions.h
	
	static double *Tn_ij,  *Tn_ik,  *Tn_jk;
	static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
	static bool   called_before = false;
	double        xdiff_ij, xdiff_ik, xdiff_jk;
	double        tempx;

	double fcut0, fcut0_ij, fcut0_ik, fcut0_jk; 
	double fcut,  fcut_ij,  fcut_ik,  fcut_jk; 	
	static string TEMP_STR;

	static int curr_triple_type_index;
	static int pow_ij, pow_ik, pow_jk;

	static string PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK;

	double RLEN_IJ;
	double RLEN_IK;
	double RLEN_JK;
	
	
	
	PAIR_TYPE_IJ = ATM_TYP_1;
	PAIR_TYPE_IK = ATM_TYP_1;
	PAIR_TYPE_JK = ATM_TYP_2;

	PAIR_TYPE_IJ.append(ATM_TYP_2);
	PAIR_TYPE_IK.append(ATM_TYP_3);
	PAIR_TYPE_JK.append(ATM_TYP_3);
	
	// Determine the FF type for the given triplet

	TEMP_STR =      FF_2BODY[ij].PRPR_NM;
	TEMP_STR.append(FF_2BODY[ik].PRPR_NM);	
	TEMP_STR.append(FF_2BODY[jk].PRPR_NM);	

	curr_triple_type_index = TRIAD_MAP[TEMP_STR];	
	
	string OUTFILE = "3b_Cheby_Pot-";
	OUTFILE.append(TEMP_STR);
	OUTFILE.append(".dat");
	ofstream OUTFILE_3B_POT;
	OUTFILE_3B_POT.open(OUTFILE.data());
	
	FULL_FILE_3B = OUTFILE;
	
	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;

		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			if (FF_2BODY[i].SNUM_3B_CHEBY > dim ) 
				dim = FF_2BODY[i].SNUM_3B_CHEBY;	
		dim++;

		Tn_ij   = new double [dim];
		Tn_ik   = new double [dim];
		Tn_jk   = new double [dim];

		Tnd_ij  = new double [dim];
		Tnd_ik  = new double [dim];
		Tnd_jk  = new double [dim];
	
	}

	tempx = 0;

	int n_ij = (FF_2BODY[ij].S_MAXIM - FF_2BODY[ij].S_MINIM)/FF_2BODY[ij].S_DELTA;
	int n_ik = (FF_2BODY[ik].S_MAXIM - FF_2BODY[ik].S_MINIM)/FF_2BODY[ik].S_DELTA;
	int n_jk = (FF_2BODY[jk].S_MAXIM - FF_2BODY[jk].S_MINIM)/FF_2BODY[jk].S_DELTA;


	int outer, middle, inner;
	
	outer  = n_ij;
	middle = n_ik;
	inner  = n_jk;
	

	for (int a=0; a<outer; a++)
	{	
		for (int b=0; b<middle; b++)
		{
			for (int c=0; c<inner; c++)
			{

				RLEN_IJ = FF_2BODY[ij].S_MINIM + a * FF_2BODY[ij].S_DELTA;
				RLEN_IK = FF_2BODY[ik].S_MINIM + b * FF_2BODY[ik].S_DELTA;
				RLEN_JK = FF_2BODY[jk].S_MINIM + c * FF_2BODY[jk].S_DELTA;

			
				SET_3B_CHEBY_POLYS(FF_2BODY[ij], Tn_ij, Tnd_ij, RLEN_IJ, xdiff_ij);
			
				SET_3B_CHEBY_POLYS(FF_2BODY[ik], Tn_ik, Tnd_ik, RLEN_IK, xdiff_ik);
			
				SET_3B_CHEBY_POLYS(FF_2BODY[jk], Tn_jk, Tnd_jk, RLEN_JK, xdiff_jk);		
			
				if ( FF_2BODY[ij].CHEBY_TYPE == "MORSE" )
				{	

					fcut0_ij     = (1.0 - RLEN_IJ/FF_2BODY[ij].S_MAXIM);
					fcut_ij      = fcut0_ij * fcut0_ij * fcut0_ij;
				
					fcut0_ik     = (1.0 - RLEN_IK/FF_2BODY[ik].S_MAXIM);
					fcut_ik      = fcut0_ik * fcut0_ik * fcut0_ik;	
				
					fcut0_jk     = (1.0 - RLEN_JK/FF_2BODY[ik].S_MAXIM);
					fcut_jk      = fcut0_jk * fcut0_jk * fcut0_jk;	
				

				
					// Error check: Certain triplets are impossible...
				
					if(curr_triple_type_index < 0)
					{
						cout << "ERROR: Impossible atom triplet found: " << TEMP_STR << endl;
						cout << "       Check the parameter file." << endl;
						exit(0);
					}
				
					// Now compute the forces for each set of allowed powers for pairs ij, ik, and jk		
					// Keep in mind that the order in which allowed powers are stored may not match the
					// ordering of pairs resulting from the present atom triplet. Thus, we need to order
					// the stored powers properly before applying the FF.
					 
					tempx = 0;
				
					for(int i=0; i<FF_3BODY[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
					{
						SET_3B_CHEBY_POWERS(FF_2BODY, FF_3BODY[curr_triple_type_index], PAIR_MAP,  pow_ij, pow_ik, pow_jk, PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK, i);
					
						tempx += FF_3BODY[curr_triple_type_index].PARAMS[i] * fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk];	
					}
					
					#if PESFORMAT == MATLAB
						OUTFILE_3B_POT_MATL_X  << 	RLEN_IJ << " ";
						OUTFILE_3B_POT_MATL_Y  << 	RLEN_IK << " ";
						OUTFILE_3B_POT_MATL_Z  << 	RLEN_JK << " ";
						OUTFILE_3B_POT_MATL_VAL<< 	tempx   << " ";	
						
						if(c > inner-1)
						{
							OUTFILE_3B_POT_MATL_X  << 	RLEN_IJ << ", ";
							OUTFILE_3B_POT_MATL_Y  << 	RLEN_IK << ", ";
							OUTFILE_3B_POT_MATL_Z  << 	RLEN_JK << ", ";
							OUTFILE_3B_POT_MATL_VAL<< 	tempx   << ", ";	
						}
					#else
						OUTFILE_3B_POT << 	RLEN_IJ << " " << RLEN_IK << " " << RLEN_JK << " " << tempx << endl;
					#endif
					
					
				}
				else
				{
					cout << "ERROR: Functionality not complete for non-Morse-type 3B Chebyshev potentials." << endl;
					exit(0);
				}			
			}
			
			#if PESFORMAT == GNUPLOT
				OUTFILE_3B_POT << endl;
			#endif				
		}
	}

	OUTFILE_3B_POT.close();
}

void Print_3B_Cheby_Scan(MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, string & ATM_TYP_1, string & ATM_TYP_2, string & ATM_TYP_3, int ij, int ik, int jk, PES_PLOTS & FF_PLOTS, int scan)
// Draws scans of 3-b potential when a single pair distance is perturbed.
// Not perturbed distances are fixed at a user-specified value
// pair distances will range from smin to smax, incremented by sdelta
{
	
	//string SCAN_FILE_3B;	// Declared as a global (external) variable in functions.h
		
	static double *Tn_ij,  *Tn_ik,  *Tn_jk;
	static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
	static bool   called_before = false;
	double        xdiff_ij, xdiff_ik, xdiff_jk;
	double        tempx;

	double fcut0, fcut0_ij, fcut0_ik, fcut0_jk; 
	double fcut,  fcut_ij,  fcut_ik,  fcut_jk; 	
	static string TEMP_STR;

	static int curr_triple_type_index;
	static int pow_ij, pow_ik, pow_jk;

	static string PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK;

	double RLEN_IJ;
	double RLEN_IK;
	double RLEN_JK;

	PAIR_TYPE_IJ = ATM_TYP_1;
	PAIR_TYPE_IK = ATM_TYP_1;
	PAIR_TYPE_JK = ATM_TYP_2;

	PAIR_TYPE_IJ.append(ATM_TYP_2);
	PAIR_TYPE_IK.append(ATM_TYP_3);
	PAIR_TYPE_JK.append(ATM_TYP_3);
	
	// Determine the FF type for the given triplet

	TEMP_STR =      FF_2BODY[ij].PRPR_NM;
	TEMP_STR.append(FF_2BODY[ik].PRPR_NM);	
	TEMP_STR.append(FF_2BODY[jk].PRPR_NM);	

	curr_triple_type_index = TRIAD_MAP[TEMP_STR];	
	
	string OUTFILE = "3b_Cheby_Pot-";
	OUTFILE.append(TEMP_STR);
	OUTFILE.append("_scan-");
		
	if      (FF_PLOTS.SCAN_PAIR[scan] == 1)
		OUTFILE.append(PAIR_TYPE_IJ);
	else if (FF_PLOTS.SCAN_PAIR[scan] == 2)
		OUTFILE.append(PAIR_TYPE_IK);
	else
		OUTFILE.append(PAIR_TYPE_JK);

	OUTFILE.append(".dat");
	ofstream OUTFILE_3B_POT;
	OUTFILE_3B_POT.open(OUTFILE.data());
	
	SCAN_FILE_3B = OUTFILE;
	
//	cout << "******** WRITING TO FILE:  " << SCAN_FILE_3B << endl;

	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;

		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			if (FF_2BODY[i].SNUM_3B_CHEBY > dim ) 
				dim = FF_2BODY[i].SNUM_3B_CHEBY;	
		dim++;

		Tn_ij   = new double [dim];
		Tn_ik   = new double [dim];
		Tn_jk   = new double [dim];

		Tnd_ij  = new double [dim];
		Tnd_ik  = new double [dim];
		Tnd_jk  = new double [dim];
	
	}

	tempx = 0;
	
	// Set the fixed distances
	
	if      (FF_PLOTS.FIX_PAIR_1[scan] == 1)
		RLEN_IJ = FF_PLOTS.FIX_VAL_1[scan];
	
	else if (FF_PLOTS.FIX_PAIR_1[scan] == 2)
		RLEN_IK = FF_PLOTS.FIX_VAL_1[scan];
	
	else
		RLEN_JK = FF_PLOTS.FIX_VAL_1[scan];
	
	
	if      (FF_PLOTS.FIX_PAIR_2[scan] == 1)
		RLEN_IJ = FF_PLOTS.FIX_VAL_2[scan];
	
	else if (FF_PLOTS.FIX_PAIR_2[scan] == 2)
		RLEN_IK = FF_PLOTS.FIX_VAL_2[scan];
	
	else
		RLEN_JK = FF_PLOTS.FIX_VAL_2[scan];	
	
	// Set the scan pair
	
	int n;
	
	if      (FF_PLOTS.SCAN_PAIR[scan] == 1)
		n = (FF_2BODY[ij].S_MAXIM - FF_2BODY[ij].S_MINIM)/FF_2BODY[ij].S_DELTA;
	
	else if (FF_PLOTS.SCAN_PAIR[scan] == 2)
		n = (FF_2BODY[ik].S_MAXIM - FF_2BODY[ik].S_MINIM)/FF_2BODY[ik].S_DELTA;
	
	else
		n = (FF_2BODY[jk].S_MAXIM - FF_2BODY[jk].S_MINIM)/FF_2BODY[jk].S_DELTA;
	
	// Run the scan	

	for (int a=1; a<n; a++)
	{

		if      (FF_PLOTS.SCAN_PAIR[scan] == 1)
			RLEN_IJ = FF_2BODY[ij].S_MINIM + a * FF_2BODY[ij].S_DELTA;
		else if (FF_PLOTS.SCAN_PAIR[scan] == 2)
			RLEN_IK = FF_2BODY[ik].S_MINIM + a * FF_2BODY[ik].S_DELTA;
		else
			RLEN_JK = FF_2BODY[jk].S_MINIM + a * FF_2BODY[jk].S_DELTA;		

			
		SET_3B_CHEBY_POLYS(FF_2BODY[ij], Tn_ij, Tnd_ij, RLEN_IJ, xdiff_ij);
	
		SET_3B_CHEBY_POLYS(FF_2BODY[ik], Tn_ik, Tnd_ik, RLEN_IK, xdiff_ik);
	
		SET_3B_CHEBY_POLYS(FF_2BODY[jk], Tn_jk, Tnd_jk, RLEN_JK, xdiff_jk);		
	
		if ( FF_2BODY[ij].CHEBY_TYPE == "MORSE" )
		{	

			fcut0_ij     = (1.0 - RLEN_IJ/FF_2BODY[ij].S_MAXIM);
			fcut_ij      = fcut0_ij * fcut0_ij * fcut0_ij;
		
			fcut0_ik     = (1.0 - RLEN_IK/FF_2BODY[ik].S_MAXIM);
			fcut_ik      = fcut0_ik * fcut0_ik * fcut0_ik;	
		
			fcut0_jk     = (1.0 - RLEN_JK/FF_2BODY[ik].S_MAXIM);
			fcut_jk      = fcut0_jk * fcut0_jk * fcut0_jk;	
		

		
			// Error check: Certain triplets are impossible...
		
			if(curr_triple_type_index < 0)
			{
				cout << "ERROR: Impossible atom triplet found: " << TEMP_STR << endl;
				cout << "       Check the parameter file." << endl;
				exit(0);
			}
		
			// Now compute the forces for each set of allowed powers for pairs ij, ik, and jk		
			// Keep in mind that the order in which allowed powers are stored may not match the
			// ordering of pairs resulting from the present atom triplet. Thus, we need to order
			// the stored powers properly before applying the FF.
			 
			tempx = 0;
			
			// If requested, include the 2-body energy contributions
		
			for(int i=0; i<FF_3BODY[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
			{
				SET_3B_CHEBY_POWERS(FF_2BODY, FF_3BODY[curr_triple_type_index], PAIR_MAP,  pow_ij, pow_ik, pow_jk, PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK, i);
			
				tempx += FF_3BODY[curr_triple_type_index].PARAMS[i] * fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk];	
			}
			
			OUTFILE_3B_POT << 	RLEN_IJ << " " << RLEN_IK << " " << RLEN_JK << " " << tempx << endl;
		}
		else
		{
			cout << "ERROR: Functionality not complete for non-Morse-type 3B Chebyshev potentials." << endl;
			exit(0);
		}			
			
	}
	
	OUTFILE_3B_POT.close();
}
 