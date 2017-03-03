#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes

#include "functions.h"

using namespace std;

//////////////////////////////////////////
//
//	SMALL UTILITY FUNCTION
//
//////////////////////////////////////////


//////////////////////////////////////////
// Overloaded error message exit functions
//////////////////////////////////////////
 
void EXIT_MSG(string EXIT_STRING) // Single error message
{
	cout << EXIT_STRING << endl;
	exit_run(0);
}

void EXIT_MSG(string EXIT_STRING, string EXIT_VAR) // error message: var
{
	cout << EXIT_STRING << ": " << EXIT_VAR << endl;
	exit_run(0);
}

void EXIT_MSG(string EXIT_STRING, double EXIT_VAR) // error message: var
{
	cout << EXIT_STRING << ": " << EXIT_VAR << endl;
	exit_run(0);
}

//////////////////////////////////////////
// Distance calculation
//////////////////////////////////////////

double get_dist(const FRAME & SYSTEM, const MD_JOB_CONTROL & CONTROLS, XYZ & RAB, int a1, int a2)
// Calculates distance as a2 - a1... This fuction modifies RAB!
{
	// Get pair distance

	RAB.X = SYSTEM.COORDS.at(a2).X - SYSTEM.COORDS.at(a1).X;
	RAB.Y = SYSTEM.COORDS.at(a2).Y - SYSTEM.COORDS.at(a1).Y;
	RAB.Z = SYSTEM.COORDS.at(a2).Z - SYSTEM.COORDS.at(a1).Z;
	
	// Apply MIC ... this method is also compatbile with coords
	// that aren't wrapped
	
	if(((!CONTROLS.COMPARE_FORCE)&&(!CONTROLS.SUBTRACT_FORCE))||(CONTROLS.N_LAYERS==0))
	{	
		// If this is a normal MD run or an LSQ run with zero layers...
		
		RAB.X -= floor( 0.5 + RAB.X/SYSTEM.BOXDIM.X )  * SYSTEM.BOXDIM.X;
		RAB.Y -= floor( 0.5 + RAB.Y/SYSTEM.BOXDIM.Y )  * SYSTEM.BOXDIM.Y;
		RAB.Z -= floor( 0.5 + RAB.Z/SYSTEM.BOXDIM.Z )  * SYSTEM.BOXDIM.Z;
	}
	else 
	{
		// When doing a compare-force calculation, apply pbc in terms 
		// of parent atoms and the primitive cell
		
		XYZ PRIM_BOX;
		XYZ PAR_DIST;
		
		if(SYSTEM.PARENT[a2]==-1)
		{
			PAR_DIST.X = RAB.X;
			PAR_DIST.Y = RAB.Y;
			PAR_DIST.Z = RAB.Z;
		}
		else
		{
			PAR_DIST.X = SYSTEM.COORDS[SYSTEM.PARENT[a2]].X - SYSTEM.COORDS.at(a1).X;
			PAR_DIST.Y = SYSTEM.COORDS[SYSTEM.PARENT[a2]].Y - SYSTEM.COORDS.at(a1).Y;
			PAR_DIST.Z = SYSTEM.COORDS[SYSTEM.PARENT[a2]].Z - SYSTEM.COORDS.at(a1).Z;	
		}

		if(!CONTROLS.IS_LSQ)
		{
			PRIM_BOX.X = SYSTEM.BOXDIM.X/(CONTROLS.N_LAYERS+1);
			PRIM_BOX.Y = SYSTEM.BOXDIM.Y/(CONTROLS.N_LAYERS+1);
			PRIM_BOX.Z = SYSTEM.BOXDIM.Z/(CONTROLS.N_LAYERS+1);
		}
		else
		{
			// If this is for LSQ, the box dimensions have NOT been
			// replicated
			
			PRIM_BOX.X = SYSTEM.BOXDIM.X;
			PRIM_BOX.Y = SYSTEM.BOXDIM.Y;
			PRIM_BOX.Z = SYSTEM.BOXDIM.Z;
		}



		RAB.X -= floor( 0.5 + PAR_DIST.X/PRIM_BOX.X )  * PRIM_BOX.X;
		RAB.Y -= floor( 0.5 + PAR_DIST.Y/PRIM_BOX.Y )  * PRIM_BOX.Y;
		RAB.Z -= floor( 0.5 + PAR_DIST.Z/PRIM_BOX.Z )  * PRIM_BOX.Z;
	}

	return sqrt( RAB.X*RAB.X + RAB.Y*RAB.Y + RAB.Z*RAB.Z );
}

//////////////////////////////////////////
// Cheby transformation functions
//////////////////////////////////////////

// Note: "Inline" tells the compiler to replace function calls with the actual contents of the function
//       to enhance efficiency.. "undoes" modularity at runtime to speed things up.

static inline double fix_cheby_val(double x, bool inverse_order)
//Takes care of cheby xformed dist behavior outside of allowed range
{
	// If the chebyshev x value is out of range, set it to a limiting value (-1 or 1)
	// This is done purely to maintain numerical stability.  
	// For 2-body interactions there is an additional repulsion for r < rmin, and
	//   a cutoff r > rmax.
	// For 3-body interactions there should be an inner cutoff function for r < rmin.
	//   and an outer cutoff function for r > rmax.
	// If inverse_order is true, x is a decreasing function of r.
	
	
	if ( x < -1.0)
	{

		// Only use these for debugging...
/*		
		if ( inverse_order )
		{
			if (isatty(fileno(stdout)))
				cout << COUT_STYLE.BOLD << COUT_STYLE.MAGENTA << "Warning: r > rmax  " << COUT_STYLE.ENDSTYLE << endl;
			else
				cout << "Warning: r > rmax  " << endl;
		}
		else
		{
			if (isatty(fileno(stdout)))
				cout << COUT_STYLE.BOLD << COUT_STYLE.MAGENTA << "Warning: r < rmin " << COUT_STYLE.ENDSTYLE << endl;
			else
				cout << "Warning: r < rmin " << endl;
		}
*/	
		x = -1.0;
	}
	else if ( x > 1.0 )
    {
		// Added a warning...Do we need to turn off the interaction, because it's beyond rcut?
   	 	// The cutoff function (handled elsewhere) we do that automatially (LEF)

		// Only use these for debugging... t
/*			
		if ( inverse_order )
		{
			if (isatty(fileno(stdout)))
				cout << COUT_STYLE.BOLD << COUT_STYLE.MAGENTA << "Warning: r < rmin  " << COUT_STYLE.ENDSTYLE << endl;
			else
				cout << "Warning: r < rmin  " << endl;
		}
		else
		{
			if (isatty(fileno(stdout)))
				cout << COUT_STYLE.BOLD << COUT_STYLE.MAGENTA << "Warning: r > rmax " << COUT_STYLE.ENDSTYLE << endl;
			else
				cout << "Warning: r > rmax " << endl;
		}
*/		
		x = 1.0;								
	}
	return x;
}

static inline void   cheby_var(double rlen, double s_minim, double s_maxim, double lambda, const string & cheby_type, double & x, double & xdiff, bool & inverse_order, double & exprlen)
// Does the actual cheby distance transformation										
{
	// Given the atomic distance rlen, the fitting minimum and maximum, the Morse lambda variable,
	//	and the type of Chebyshev approximant, calculate:
	// the chebyshev variable x, 
	// the difference between maximum and minimum variable randes (xdiff)
	// a flag indicating the order of max/min variables (whether the chebyshev variable
	// is increasing or decreasing with rlen) (inverse_order)
	// the exponential of the len variable if appropriate.
	//
	// Chebyshev polynomials are only defined on the range [-1,1], so transorm the pair distance
	// in a way that allows it to fall along that range. Options are:
	//
	// x = 1/pair_dist				// Inverse r, 
	// x = exp(pair_dist/lambda)	// Morse-type
	// x = pair_dist				// default type
	// 
	// All types are normalized by s_min to s_max range to fall along [-1,1]	
	
	
	double xavg, xmin, xmax ;

	if ( cheby_type == "INVRSE_R" ) 
	{
		xavg  =  0.5 * (1.0/s_minim + 1.0/s_maxim); // midpoint of possible pair distances in r^-1 space
		xdiff =  0.5 * (1.0/s_minim - 1.0/s_maxim); // width of possible pair distances in r^-1 space
		x     = (1.0/rlen-xavg) / xdiff;			// pair distances in r^-1 space, normalized to fit over [-1,1]
		inverse_order = true;
		exprlen = 0.0 ;
	} 
	else if ( cheby_type == "MORSE" ) 
	{
		xmin  = exp(-s_maxim/lambda); 
		xmax  = exp(-s_minim/lambda); 
		xavg  = 0.5 * (xmin + xmax);				// midpoint of possible pair distances in morse space
		xdiff = 0.5 * (xmax - xmin);				// width of possible pair distances in morse space
		exprlen = exp(-rlen/lambda);
		x = (exprlen-xavg)/xdiff;					// pair distances in morse space, normalized to fit over [-1,1]
		inverse_order = true;
	}
	else if (cheby_type == "DEFAULT")
	{
		xavg  = 0.5 * (s_minim + s_maxim); 			// midpoint of possible pair distances
		xdiff = 0.5 * (s_maxim - s_minim); 			// width of possible pair distances
		x = (rlen-xavg) / xdiff;	
		exprlen = 0.0 ;
		inverse_order = false;
	}
	else
	{
		cout << "ERROR: Undefined CHBTYPE: " << cheby_type << endl;
		cout << "       Excepted values are \"DEFAULT\", \"INVRSE_R\", or \"MORSE\". " << endl;
		exit_run(1);
	}
}

static inline double cheby_var_deriv(double xdiff, double rlen, double lambda, string& cheby_type)
// Calculate the derivative of the cheby variable x with respect to rlen.
{
	double dx_dr;
	
	if ( cheby_type == "MORSE" )
			//dx_dr =  -exp(-rlen/lambda) / (lambda * xdiff);
	        dx_dr =  (-exp(-rlen/lambda)/lambda)/xdiff;

	else if ( cheby_type == "INVRSE_R" ) 
		dx_dr = -1.0/(rlen * rlen * xdiff);
    
	else if ( cheby_type  == "DEFAULT" )
		dx_dr = 1.0 / xdiff;
    
	else
    {
		cout << "Error: bad cheby_type: " << cheby_type << endl;
		exit_run(1);
    }
	
	return dx_dr;
}


//////////////////////////////////////////
// MPI compatibility functions
//////////////////////////////////////////
 
void divide_atoms(int &a1start, int &a1end, int atoms) 
{
	a1start = RANK * atoms / NPROCS;

	if ( NPROCS > atoms ) 
	{
		cout << "Error: number of processors > number of atoms: " << endl;
		cout << NPROCS << " " << atoms << endl;
		exit_run(1);
	}
	  
	if ( RANK == NPROCS - 1 ) 
		a1end = atoms - 1;

	else 
		a1end   = ((RANK+1) * atoms / NPROCS) - 1;

	//  cout << "DIVIDING ATOMS: RANK : " << RANK << " " << a1start << ":" << a1end << endl;
}


//////////////////////////////////////////
// 3B Cheby functions
//////////////////////////////////////////
 
void SET_3B_CHEBY_POLYS( PAIRS & FF_2BODY, double *Tn, double *Tnd, const double rlen, double & xdiff, double SMAX, double SMIN)
// Sets the value of the Chebyshev polynomials (Tn) and thier derivatives (Tnd)
{
	double x, exprlen;
	bool   inverse_order = true;  // Inverse order determines whether x is an increasing (false) or decreasing (true) function of r (LEF).
	double CHEBY_DERIV_CONST;	// Accounts for when cheby range is changed from -1:1 to x:y

	CHEBY_DERIV_CONST = FF_2BODY.CHEBY_RANGE_HIGH - FF_2BODY.CHEBY_RANGE_LOW;	// Ranges should be the same for all types
	CHEBY_DERIV_CONST /= 2.0; // i.e the width of the default cheby range 
	
	// Do the Cheby distance transformation
	
	cheby_var(
		rlen, 
		SMIN, 
		SMAX, 
		FF_2BODY.LAMBDA, 
		FF_2BODY.CHEBY_TYPE, 
		x, 
		xdiff, 
		inverse_order, 
		exprlen);
		
	#if CHECK_CHEBY_RANGE == 1	
	
		if ( x < -1.0 || x > 1.0 ) 
			x = fix_cheby_val(x, inverse_order);

		// Now change the range, if the user requested
	
		x = x*CHEBY_DERIV_CONST + FF_2BODY.CHEBY_RANGE_LOW - -1.0*CHEBY_DERIV_CONST;
	
		// Sanity check
	
		if ( x < FF_2BODY.CHEBY_RANGE_LOW || x > FF_2BODY.CHEBY_RANGE_HIGH )
		{
			cout << "ERROR: transformed (3B) x falls outside user-defined range." << endl;
			cout << "x: " << x << endl;
			cout << "high/low: " << FF_2BODY.CHEBY_RANGE_HIGH << " " << FF_2BODY.CHEBY_RANGE_LOW  << endl;
			exit_run(0);
		}	
		
	#endif
	
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

double SET_SMAXIM(PAIRS & FF_2BODY, TRIPLETS & PAIR_TRIPLETS, string TYPE)
// Used for 3-body interactions. Decides whether outer cutoff should be set by 2-body value or 3-body value. Returns the cutoff value.
{
	
	/*
	if(PAIR_TRIPLETS.S_MAXIM_3B == -1)
		return FF_2BODY.S_MAXIM;
	else
		return PAIR_TRIPLETS.S_MAXIM_3B;
	*/
	
	
	double VAL;
	
	if(PAIR_TRIPLETS.S_MAXIM_3B.X == -1)
		VAL =  FF_2BODY.S_MAXIM;
	else
	{
		if      (TYPE == PAIR_TRIPLETS.ATMPAIR1)
			VAL =  PAIR_TRIPLETS.S_MAXIM_3B.X;
		else if (TYPE == PAIR_TRIPLETS.ATMPAIR2)
			VAL =  PAIR_TRIPLETS.S_MAXIM_3B.Y;
		else if (TYPE == PAIR_TRIPLETS.ATMPAIR3)
			VAL =  PAIR_TRIPLETS.S_MAXIM_3B.Z;
	}
	return VAL;	
}

double SET_SMINIM(PAIRS & FF_2BODY, TRIPLETS & PAIR_TRIPLETS, string TYPE)
// Used for 3-body interactions. Decides whether outer cutoff should be set by 2-body value or 3-body value. Returns the cutoff value.
{
	double VAL;
	
	if(PAIR_TRIPLETS.S_MINIM_3B.X == -1)
		VAL =  FF_2BODY.S_MINIM;
	else
	{
		if      (TYPE == PAIR_TRIPLETS.ATMPAIR1)
			VAL =  PAIR_TRIPLETS.S_MINIM_3B.X;
		else if (TYPE == PAIR_TRIPLETS.ATMPAIR2)
			VAL =  PAIR_TRIPLETS.S_MINIM_3B.Y;
		else if (TYPE == PAIR_TRIPLETS.ATMPAIR3)
			VAL =  PAIR_TRIPLETS.S_MINIM_3B.Z;
	}
	return VAL;
}

//////////////////////////////////////////
// Pressure functions
//////////////////////////////////////////

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

double numerical_pressure(const FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, NEIGHBORS & NEIGHBOR_LIST)   
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
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP, NEIGHBOR_LIST);
	
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
	
	ZCalc(REPLICATE, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP, NEIGHBOR_LIST);
	
	Vtot2 = REPLICATE.TOT_POT_ENER;

	// compute (return) pressure 
	
	return -(Vtot2 - Vtot1)/(Vol2 - Vol1);

}

//////////////////////////////////////////
//
//	FUNCTION HEADERS -- DERIVATIVE CALCULATION
//
//////////////////////////////////////////

static void ZCalc_Spline_Deriv  (MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);

static void ZCalc_Cheby_Deriv   (MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);

static void ZCalc_3B_Cheby_Deriv(MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<TRIPLETS> & PAIR_TRIPLETS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, map<string,int> TRIAD_MAP);

static void ZCalc_Poly_Deriv    (MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);	

static void ZCalc_InvR_Deriv    (MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP);


//////////////////////////////////////////
//
//	FUNCTION HEADERS -- FORCE CALCULATION
//
//////////////////////////////////////////

static void ZCalc_Cheby(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

static void ZCalc_3B_Cheby(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, NEIGHBORS & NEIGHBOR_LIST);	

static void ZCalc_Lj(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP);

static void ZCalc_Spline(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST);

static void ZCalcSR_Over(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP);

////////////////////////////////////////////////////////////
//
// FUNCTIONS -- DERIVATIVE CALCULATION
//
////////////////////////////////////////////////////////////


// FUNCTION UPDATED
void ZCalc_Deriv (MD_JOB_CONTROL & CONTROLS, vector<PAIRS> & FF_2BODY,  vector<TRIPLETS> & PAIR_TRIPLETS, FRAME & FRAME_SYSTEM, vector<vector <XYZ > > & FRAME_A_MATRIX, vector<vector <XYZ > > & FRAME_COULOMB_FORCES, const int nlayers, bool if_3b_cheby, map<string,int> & PAIR_MAP,  map<string,int> & TRIAD_MAP)
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
		ZCalc_Spline_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP);
	
	else if ( FF_2BODY[0].PAIRTYP == "CHEBYSHEV" )	// FUNCTION UPDATED
	{
		// Only enter if 2B are requested. For example, skip if user wants to fit ONLY 3B cheby
		// i.e. PAIRTYP: CHEBYSHEV  0 6 or similar
		//
		
		if ( FF_2BODY[0].SNUM > 0)
			ZCalc_Cheby_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP);
	
		if (if_3b_cheby)
			ZCalc_3B_Cheby_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, PAIR_TRIPLETS, FRAME_A_MATRIX, nlayers, PAIR_MAP, TRIAD_MAP);		
	}			

    else if ( FF_2BODY[0].PAIRTYP == "DFTBPOLY" )	// FUNCTION UPDATED
		ZCalc_Poly_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP);

    else if ( FF_2BODY[0].PAIRTYP == "INVRSE_R" )	// FUNCTION UPDATED
		ZCalc_InvR_Deriv(CONTROLS, FRAME_SYSTEM, FF_2BODY, FRAME_A_MATRIX, nlayers, PAIR_MAP);

    else 
    {
		cout << "Error: bad pairtype in ZCalc_Deriv\n";
		exit(1);
    }		
}	

// FUNCTION UPDATED
static void ZCalc_Spline_Deriv(MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)		   
{
	// Original comment: Calculate derivatives of the forces wrt the spline parameters. Stores minimum distance between a pair of atoms in minD[i].
	// New Note: This doesn't actually calcuate any derivatives.. it is just populating A with the cubic hermite basis polynomials needed for fitting
	// 3rd note:  The hermite basis polynomials stored in FRAME_A_MATRIX are equal to derivatives of the force with respect to the spline coefficients..
	//            This is what all of the ZCalc_XX_Deriv functions find.
	//            Since the force is a linear function of the spline coefficients, you are left with just the cubic hermite basis polynomials.
	//            (LEF)
	
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
	
	// Set up for layering

	int fidx_a2;
	int a2start, a2end, a2;
	int TOTAL_REP_ATOMS = SYSTEM.COORDS.size();

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)		// Double sum over atom pairs
	{
		a2start = a1+1;
		a2end   = TOTAL_REP_ATOMS;

		
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = a2idx;
		
			if(CONTROLS.N_LAYERS>0)
				if(SYSTEM.PARENT[a2] != -1)
					if(SYSTEM.PARENT[a2]<a1)
						continue;
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
			//calculate vstart:

			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Get pair distance

			rlen = get_dist(SYSTEM, CONTROLS, RAB, a1, a2);	// Updates RAB!

			if ( (rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST) && (a1 !=a2) )	
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
				
				
				fidx_a2 = a2;
	
				if(CONTROLS.N_LAYERS>0)
					if(SYSTEM.PARENT[a2] != -1)
							fidx_a2 = SYSTEM.PARENT[a2];
				
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
				
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+0].X -= h00 * RAB.X / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+1].X -= h10 * RAB.X / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+2].X -= h01 * RAB.X / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+3].X -= h11 * RAB.X / rlen;
			  
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+0].Y -= h00 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+1].Y -= h10 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+2].Y -= h01 * RAB.Y / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+3].Y -= h11 * RAB.Y / rlen;
				  
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+0].Z -= h00 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+1].Z -= h10 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+2].Z -= h01 * RAB.Z / rlen;
			   	FRAME_A_MATRIX[fidx_a2][vstart+kstart+3].Z -= h11 * RAB.Z / rlen;							
				
			}
		}
	}
	return;
}

// FUNCTION UPDATED
static void ZCalc_Cheby_Deriv(MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)	
// Calculate derivatives of the forces wrt the Chebyshev parameters. Stores minimum distance between a pair of atoms in minD[i].
{
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double x;
	int vstart;
	static double *Tn, *Tnd;
	static bool called_before = false;

	double xdiff, exprlen;
	double fcut0; 
	double fcut; 
	double fcutderiv; 				
	double deriv;
	double tmp_doub; 	
	double dx_dr;
	
	static double CHEBY_DERIV_CONST;	// Accounts for when cheby range is changed from -1:1 to x:y

	bool inverse_order;
	
	const double fcut_power = 
    #ifndef FPENALTY_POWER
		3.0;
	#else
		FPENALTY_POWER;
    #endif

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
		
		CHEBY_DERIV_CONST = FF_2BODY[0].CHEBY_RANGE_HIGH - FF_2BODY[0].CHEBY_RANGE_LOW;	// Ranges should be the same for all types
		CHEBY_DERIV_CONST /= 2.0; // i.e the width of the default cheby range
	}

	// Main loop for Chebyshev terms:
	
	string TEMP_STR;
	int curr_pair_type_idx;
	
	// Set up for layering

	int fidx_a2;
	int a2start, a2end, a2;
	int TOTAL_REP_ATOMS = SYSTEM.COORDS.size();


	for(int a1=0;a1<SYSTEM.ATOMS;a1++)		// Double sum over atom pairs
	{
		a2start = a1+1;
		a2end   = TOTAL_REP_ATOMS;

			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = a2idx;
			
			if(CONTROLS.N_LAYERS>0)
				if(SYSTEM.PARENT[a2] != -1)
					if(SYSTEM.PARENT[a2]<a1)
						continue;

			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
		
			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Get pair distance

			rlen = get_dist(SYSTEM, CONTROLS, RAB, a1, a2);	// Updates RAB!

			if ( (rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST))	
				FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;
			
			

			if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM and rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
			{
				FF_2BODY[curr_pair_type_idx].N_CFG_CONTRIB++;
			
				// Do the distance transformation
				
				cheby_var(rlen, 
							FF_2BODY[curr_pair_type_idx].S_MINIM,
							FF_2BODY[curr_pair_type_idx].S_MAXIM,
							FF_2BODY[curr_pair_type_idx].LAMBDA,
							FF_2BODY[curr_pair_type_idx].CHEBY_TYPE,
							x, xdiff, inverse_order, exprlen);
							
				#if CHECK_CHEBY_RANGE == 1				
				
			   		// Why are we setting it as -1? shouldn't something else be done here instead? ...since we're outside the range of the fit??
			    	// This is mostly defensive programming against bad behavior of the Cheby polynomial
			    	// outside of the fitting range.  If r < rmin, the repulsive potential should take over.
			    	// If r > rmax, the cutoff function will set the total force to 0. (LEF)							

			  		if ( x < -1.0 || x > 1.0 )
						x = fix_cheby_val(x, inverse_order);
				
					// Now change the range, if the user requested

					x = x*CHEBY_DERIV_CONST + FF_2BODY[0].CHEBY_RANGE_LOW - -1.0*CHEBY_DERIV_CONST;							
				
					// Sanity check
				
					if ( x < FF_2BODY[0].CHEBY_RANGE_LOW || x > FF_2BODY[0].CHEBY_RANGE_HIGH )
					{
						cout << "ERROR: transformed x falls outside user-defined range." << endl;
						cout << "x: " << x << endl;
						cout << "high/low: " << FF_2BODY[0].CHEBY_RANGE_HIGH << " " << FF_2BODY[0].CHEBY_RANGE_LOW  << endl;
						exit(0);
					}
				#endif

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
				// the SYSTEM from heading towards though poorly sampled regions of the PES, but also to ensure that the PES goes to
				// zero at rcut
				
			  	// fcut and fcutderv are the cutoff functions (1-r/rcut)**3 and its
			  	// derivative -3 (1-r/rcut)**2/rcut.  This ensures that
			 	// the force goes to 0 as r goes to rcut.
			  	// This is not a penalty function in the usual sense.  
			  	// I don't see any reason to have a scaling on the cutoff.
			  
			  	// That will simply multiply all the forces by a constant,
			  	// which will be canceled out during the force matching process.
			 	// (LEF)
				
				// fcut and fcutderv are the form that the penalty func and its derivative for the morse-type pair distance transformation
				
				fcut0     = (1.0 - rlen/FF_2BODY[curr_pair_type_idx].S_MAXIM);
				
				fcut      = pow(fcut0, fcut_power);
				fcutderiv = pow(fcut0, fcut_power-1);
				
				fcutderiv *= -1.0 * fcut_power / FF_2BODY[curr_pair_type_idx].S_MAXIM;
				
				// Compute part of the derivative
				
			 	// Added missing code for other cheby types (LEF).
				// NOTE: All these extra terms are coming from:
				//
				// 1. Chain rule to account for transformation from morse-type pair distance to x
				// 2. Product rule coming from pair distance dependence of fcut, the penalty function
				
			 	dx_dr = CHEBY_DERIV_CONST*cheby_var_deriv(xdiff, rlen, FF_2BODY[curr_pair_type_idx].LAMBDA, FF_2BODY[curr_pair_type_idx].CHEBY_TYPE);
				
				fidx_a2 = a2;
				
				if(CONTROLS.N_LAYERS>0)
					if(SYSTEM.PARENT[a2] != -1)
							fidx_a2 = SYSTEM.PARENT[a2];



				for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
				{
					tmp_doub = (fcut * Tnd[i+1] * dx_dr + fcutderiv * Tn[i+1] );
					//tmp_doub = (fcut * Tnd[i+1] *(-exp(-rlen/FF_2BODY[curr_pair_type_idx].LAMBDA)/FF_2BODY[curr_pair_type_idx].LAMBDA)/xdiff + fcutderiv * Tn[i+1] );
					
					// Finally, account for the x, y, and z unit vectors
					
					deriv = tmp_doub * RAB.X / rlen;
					FRAME_A_MATRIX[a1     ][vstart+i].X += deriv;
					FRAME_A_MATRIX[fidx_a2][vstart+i].X -= deriv;

					deriv = tmp_doub * RAB.Y / rlen; 
					FRAME_A_MATRIX[a1     ][vstart+i].Y += deriv;
					FRAME_A_MATRIX[fidx_a2][vstart+i].Y -= deriv;
														
					deriv = tmp_doub * RAB.Z / rlen;
					FRAME_A_MATRIX[a1     ][vstart+i].Z += deriv;
					FRAME_A_MATRIX[fidx_a2][vstart+i].Z -= deriv;

				}
			}

		}
	}

  return;
}

static void ZCalc_3B_Cheby_Deriv(MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<TRIPLETS> & PAIR_TRIPLETS, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP, map<string,int> TRIAD_MAP)		
// Calculate derivatives of the forces wrt the 3-body Chebyshev parameters. 
{
	// This three body interaction stems from: C_n^ij *  C_n^ik * C_n^jk * T_n(x_ij) * T_n(x_ik) * T_n(x_jk)
	//	The logic:
	//	+ Run a triple loop over all atoms in the system.
	//	+ Compute C_ij, C_ik, and C_jk coeffiecients independently as you would do for a normal 2 body 

	XYZ RVEC_IJ;
	XYZ RVEC_IK;
	XYZ RVEC_JK;
	 		
	XYZ RAB_IJ;
	XYZ RAB_IK;
	XYZ RAB_JK; 		

	double rlen_ij,  rlen_ik,  rlen_jk;
	int vstart;
    static int n_2b_cheby_terms, n_3b_cheby_terms;
	static double *Tn_ij,  *Tn_ik,  *Tn_jk;
	static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
	static bool called_before = false;
	
	static int pow_ij, pow_ik, pow_jk;
	double dx_dr_ij, dx_dr_jk, dx_dr_ik;
	double xdiff_ij, xdiff_ik, xdiff_jk; 
	double fcut0_ij, fcut0_ik, fcut0_jk; 
	double fcut_ij,  fcut_ik,  fcut_jk; 			
	double deriv_ij, deriv_ik, deriv_jk;
	double fcutderiv_ij, fcutderiv_ik, fcutderiv_jk; 	
	
	static string TEMP_STR;
	static string PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK;
	static int curr_triple_type_index;
	static int curr_pair_type_idx_ij;
	static int curr_pair_type_idx_ik;
	static int curr_pair_type_idx_jk;
	static int row_offset;	
	
	static int  N_LAYER_IJK;
	static vector<XYZ_INT> LAYER_IJ;
	static vector<XYZ_INT> LAYER_IK;
	static vector<XYZ_INT> LAYER_JK;
	
	static double CHEBY_DERIV_CONST;	// Accounts for when cheby range is changed from -1:1 to x:y

	double S_MAXIM_IJ, S_MAXIM_IK, S_MAXIM_JK;
	double S_MINIM_IJ, S_MINIM_IK, S_MINIM_JK;

	const double fcut_power = 
    #ifndef FPENALTY_POWER
		3.0;
	#else
		FPENALTY_POWER;
    #endif
	
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
		
		CHEBY_DERIV_CONST = FF_2BODY[0].CHEBY_RANGE_HIGH - FF_2BODY[0].CHEBY_RANGE_LOW;	// Ranges should be the same for all types
		CHEBY_DERIV_CONST /= 2.0; // i.e the width of the default cheby range		

	}

	// Set up for layering

	int fidx_a2, fidx_a3;
	
	// Set up for MPI
	
	int a1start, a1end;	

	divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.

	// Set up for neighbor lists
	
	int a2start, a2end, a2;
	int a3start, a3end, a3;
	
	int TOTAL_REP_ATOMS = SYSTEM.COORDS.size();

	for(int a1=a1start; a1<=a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev -1)
	{
		a2start = a1+1;
		a2end   = TOTAL_REP_ATOMS;

			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = a2idx;

			a3start = a2+1;
			a3end   = TOTAL_REP_ATOMS;
			

			if(CONTROLS.N_LAYERS>0)
				if(SYSTEM.PARENT[a2] != -1)
					if(SYSTEM.PARENT[a2]<a1)
						continue;

			for(int a3idx=a3start; a3idx<a3end; a3idx++)	
			{			
				a3 = a3idx;
				

				if(CONTROLS.N_LAYERS>0)
					if(SYSTEM.PARENT[a3] != -1)
						if(SYSTEM.PARENT[a3]<SYSTEM.PARENT[a2])
							continue;

				
				TEMP_STR = SYSTEM.ATOMTYPE[a1];
				TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);		
				PAIR_TYPE_IJ = TEMP_STR;					
				curr_pair_type_idx_ij = PAIR_MAP[TEMP_STR];
	
				TEMP_STR = SYSTEM.ATOMTYPE[a1];
				TEMP_STR.append(SYSTEM.ATOMTYPE[a3]);	
				PAIR_TYPE_IK = TEMP_STR;							
				curr_pair_type_idx_ik = PAIR_MAP[TEMP_STR];	
			
				TEMP_STR = SYSTEM.ATOMTYPE[a2];
				TEMP_STR.append(SYSTEM.ATOMTYPE[a3]);	
				PAIR_TYPE_JK = TEMP_STR;							
				curr_pair_type_idx_jk = PAIR_MAP[TEMP_STR];		
		
				// Determine the FF type for the given triplet

				TEMP_STR =      FF_2BODY[curr_pair_type_idx_ij].PRPR_NM;
				TEMP_STR.append(FF_2BODY[curr_pair_type_idx_ik].PRPR_NM);	
				TEMP_STR.append(FF_2BODY[curr_pair_type_idx_jk].PRPR_NM);	

				curr_triple_type_index = TRIAD_MAP[TEMP_STR];	
				
				// If this type has been excluded, then skip to the next iteration of the loop
				
				if(curr_triple_type_index<0)
					continue;
				
				rlen_ij = get_dist(SYSTEM, CONTROLS, RAB_IJ, a1, a2);	// Updates RAB!
				rlen_ik = get_dist(SYSTEM, CONTROLS, RAB_IK, a1, a3);	// Updates RAB!
				rlen_jk = get_dist(SYSTEM, CONTROLS, RAB_JK, a2, a3);	// Updates RAB!
	
				S_MAXIM_IJ = SET_SMAXIM(FF_2BODY[curr_pair_type_idx_ij], PAIR_TRIPLETS[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ij].PRPR_NM);
				S_MAXIM_IK = SET_SMAXIM(FF_2BODY[curr_pair_type_idx_ik], PAIR_TRIPLETS[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ik].PRPR_NM);
				S_MAXIM_JK = SET_SMAXIM(FF_2BODY[curr_pair_type_idx_jk], PAIR_TRIPLETS[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ij].PRPR_NM);
				
				S_MINIM_IJ = SET_SMINIM(FF_2BODY[curr_pair_type_idx_ij], PAIR_TRIPLETS[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ij].PRPR_NM);
				S_MINIM_IK = SET_SMINIM(FF_2BODY[curr_pair_type_idx_ik], PAIR_TRIPLETS[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ik].PRPR_NM);
				S_MINIM_JK = SET_SMINIM(FF_2BODY[curr_pair_type_idx_jk], PAIR_TRIPLETS[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_jk].PRPR_NM);
				
				// Before doing any polynomial/coeff set up, make sure that all ij, ik, and jk distances are 
				// within the allowed range.
				// Unlike the 2-body Cheby, we only evaluate this for r > s_minim.  We
				// might want an inner cutoff function to ensure energy conservation (LEF)
				// To clarify, this refers to the evaluation in the MD functions

				if(rlen_ij > S_MINIM_IJ && rlen_ij < S_MAXIM_IJ)
				{		
					if(rlen_ik > S_MINIM_IK && rlen_ik < S_MAXIM_IK)
					{	
						if(rlen_jk > S_MINIM_JK && rlen_jk < S_MAXIM_JK)											
						{		
							
							
							// Everything is within allowed ranges.
							
							// Add this to the number of configs contributing to a fit for this triplet type
							
							PAIR_TRIPLETS[curr_triple_type_index].N_CFG_CONTRIB++;

							// Begin setting up the derivative calculation

							// Set up the polynomials
			
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_ij], Tn_ij, Tnd_ij, rlen_ij, xdiff_ij, S_MAXIM_IJ, S_MINIM_IJ);
			
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_ik], Tn_ik, Tnd_ik, rlen_ik, xdiff_ik, S_MAXIM_IK, S_MINIM_IK);
			
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_jk], Tn_jk, Tnd_jk, rlen_jk, xdiff_jk, S_MAXIM_JK, S_MINIM_JK);				

							// At this point we've completed all pre-calculations needed to populate the A matrix. Now we need to figure out 
							// where within the matrix to put the data, and to do so. 

							// Note: This syntax is safe since there is only one possible SNUM_3B_CHEBY value for all interactions

							vstart = n_2b_cheby_terms;
			
							for (int i=0; i<curr_triple_type_index; i++)
								vstart += PAIR_TRIPLETS[i].N_TRUE_ALLOWED_POWERS;			

							// Set the pentaly function values
			
							fcut0_ij     = (1.0 - rlen_ij/S_MAXIM_IJ);
							fcut_ij      = pow(fcut0_ij, fcut_power);
							fcutderiv_ij = pow(fcut0_ij, fcut_power-1);
							fcutderiv_ij *= -1.0 * fcut_power / S_MAXIM_IJ;		
			
							fcut0_ik     = (1.0 - rlen_ik/S_MAXIM_IK);
							fcut_ik      = pow(fcut0_ik, fcut_power);
							fcutderiv_ik = pow(fcut0_ik, fcut_power-1);
							fcutderiv_ik *= -1.0 * fcut_power /S_MAXIM_IK;	
			
							fcut0_jk     = (1.0 - rlen_jk/S_MAXIM_JK);
							fcut_jk      = pow(fcut0_jk, fcut_power);
							fcutderiv_jk = pow(fcut0_jk, fcut_power-1);
							fcutderiv_jk *= -1.0 * fcut_power / S_MAXIM_JK;		
			
							/////////////////////////////////////////////////////////////////////
							/////////////////////////////////////////////////////////////////////
							// Consider special restrictions on allowed triplet types and powers
							/////////////////////////////////////////////////////////////////////
							/////////////////////////////////////////////////////////////////////
			
							row_offset = 0;
			
							// --- THE KEY HERE IS TO UNDERSTAND THAT THE IJ, IK, AND JK HERE IS BASED ON ATOM PAIRS, AND DOESN'T NECESSARILY MATCH THE TRIPLET'S EXPECTED ORDER!
			
							dx_dr_ij = CHEBY_DERIV_CONST*cheby_var_deriv(xdiff_ij, rlen_ij, FF_2BODY[curr_pair_type_idx_ij].LAMBDA, FF_2BODY[curr_pair_type_idx_ij].CHEBY_TYPE);

							dx_dr_jk = CHEBY_DERIV_CONST*cheby_var_deriv(xdiff_jk, rlen_jk, FF_2BODY[curr_pair_type_idx_jk].LAMBDA, FF_2BODY[curr_pair_type_idx_jk].CHEBY_TYPE);

							dx_dr_ik = CHEBY_DERIV_CONST*cheby_var_deriv(xdiff_ik, rlen_ik, FF_2BODY[curr_pair_type_idx_ik].LAMBDA, FF_2BODY[curr_pair_type_idx_ik].CHEBY_TYPE);
							
							fidx_a2 = a2;
							fidx_a3 = a3;
							
							if(CONTROLS.N_LAYERS>0)
							{
								if(SYSTEM.PARENT[a2] != -1)
										fidx_a2 = SYSTEM.PARENT[a2];
								
								if(SYSTEM.PARENT[a3] != -1)
										fidx_a3 = SYSTEM.PARENT[a3];
							}

							for(int i=0; i<PAIR_TRIPLETS[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
							{
							    row_offset = PAIR_TRIPLETS[curr_triple_type_index].PARAM_INDICIES[i];

							    SET_3B_CHEBY_POWERS(FF_2BODY, PAIR_TRIPLETS[curr_triple_type_index], PAIR_MAP, pow_ij, pow_ik, pow_jk, PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK, i);

							    deriv_ij =  fcut_ij * Tnd_ij[pow_ij] * dx_dr_ij + fcutderiv_ij * Tn_ij[pow_ij];

							    deriv_ik =  fcut_ik * Tnd_ik[pow_ik] * dx_dr_ik + fcutderiv_ik * Tn_ik[pow_ik];

							    deriv_jk =  fcut_jk * Tnd_jk[pow_jk] * dx_dr_jk + fcutderiv_jk * Tn_jk[pow_jk];	
						
							    // ij pairs

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].X += deriv_ij * fcut_ik * fcut_jk
							      * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.X / rlen_ij;

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].X -= deriv_ij * fcut_ik * fcut_jk
							      * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.X / rlen_ij;

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].Y += deriv_ij * fcut_ik * fcut_jk
							      * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.Y / rlen_ij;

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Y -= deriv_ij * fcut_ik * fcut_jk
							      * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.Y / rlen_ij;

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].Z += deriv_ij * fcut_ik * fcut_jk
							      * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.Z / rlen_ij;

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Z -= deriv_ij * fcut_ik * fcut_jk
							      * Tn_ik[pow_ik] * Tn_jk[pow_jk] * RAB_IJ.Z / rlen_ij;	


							    // ik pairs

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].X += deriv_ik * fcut_ij * fcut_jk
							      * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.X / rlen_ik;

							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].X -= deriv_ik * fcut_ij * fcut_jk
							      * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.X / rlen_ik;

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].Y += deriv_ik * fcut_ij * fcut_jk
							      * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.Y / rlen_ik;

							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Y -= deriv_ik * fcut_ij * fcut_jk
							      * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.Y / rlen_ik;

							    FRAME_A_MATRIX[a1     ][vstart+row_offset].Z += deriv_ik * fcut_ij * fcut_jk
							      * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.Z / rlen_ik;

							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Z -= deriv_ik * fcut_ij * fcut_jk
							      * Tn_ij[pow_ij] * Tn_jk[pow_jk] * RAB_IK.Z / rlen_ik;

							    // jk pairs

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].X += deriv_jk * fcut_ij * fcut_ik
							      * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.X / rlen_jk;

							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].X -= deriv_jk * fcut_ij * fcut_ik
							      * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.X / rlen_jk;

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Y += deriv_jk * fcut_ij * fcut_ik
							      * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.Y / rlen_jk;

							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Y -= deriv_jk * fcut_ij * fcut_ik
							      * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.Y / rlen_jk;

							    FRAME_A_MATRIX[fidx_a2][vstart+row_offset].Z += deriv_jk * fcut_ij
							      * fcut_ik * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.Z / rlen_jk;

							    FRAME_A_MATRIX[fidx_a3][vstart+row_offset].Z -= deriv_jk * fcut_ij * fcut_ik
							      * Tn_ij[pow_ij] * Tn_ik[pow_ik] * RAB_JK.Z / rlen_jk;

							}
		
						} // end if rlen_jk within cutoffs...
					} // end if rlen_ik within cutoffs...	

				} // end third loop over atoms					
			}
		}
	}
}

// FUNCTION UPDATED
static void ZCalc_InvR_Deriv(MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)		
// Calculate derivatives of the forces wrt to inverse pair distance to various powers. Stores minimum distance between a pair of atoms in minD[i].
{
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	int vstart;
	int curr_pair_type_idx;
	string TEMP_STR;

	double fc;
	double dfc;
	double rfac;

	// Set up for layering

	int fidx_a2;
	int a2start, a2end, a2;
	int TOTAL_REP_ATOMS = SYSTEM.COORDS.size();

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)		// Double sum over atom pairs
	{
		a2start = a1+1;
		a2end   = TOTAL_REP_ATOMS;

			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = a2idx;
			
			if(CONTROLS.N_LAYERS>0)
				if(SYSTEM.PARENT[a2] != -1)
					if(SYSTEM.PARENT[a2]<a1)
						continue;

			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
		
			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Get pair distance

			rlen = get_dist(SYSTEM, CONTROLS, RAB, a1, a2);	// Updates RAB!

			if ( (rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST) && (a1 !=a2) )	// spline term calculated w/cutoff:
				FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

			if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM && rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
			{
				// Calculate the penalty function (fc) and its derivative
				
				fc = (rlen - FF_2BODY[curr_pair_type_idx].S_MAXIM);
				dfc = fc;
				
				fc = fc*fc*fc;							
				dfc = 3*dfc*dfc;
				
				fidx_a2 = a2;
	
				if(CONTROLS.N_LAYERS>0)
					if(SYSTEM.PARENT[a2] != -1)
							fidx_a2 = SYSTEM.PARENT[a2];
			
				for (int i=0; i<FF_2BODY[curr_pair_type_idx].SNUM; i++) 
				{
					rfac = ( (i+2) / pow(rlen,i+3) )*fc;
					rfac -= (1/pow(rlen,i+2))*dfc;

					FRAME_A_MATRIX[a1][vstart+i].X += rfac * RAB.X/rlen;
					FRAME_A_MATRIX[a1][vstart+i].Y += rfac * RAB.Y/rlen;
					FRAME_A_MATRIX[a1][vstart+i].Z += rfac * RAB.Z/rlen;
						
					FRAME_A_MATRIX[fidx_a2][vstart+i].X -= rfac * RAB.X/rlen;
					FRAME_A_MATRIX[fidx_a2][vstart+i].Y -= rfac * RAB.Y/rlen;
					FRAME_A_MATRIX[fidx_a2][vstart+i].Z -= rfac * RAB.Z/rlen;								
				}
			}
		}
	}
	return;

}

// FUNCTION UPDATED
static void ZCalc_Poly_Deriv(MD_JOB_CONTROL & CONTROLS, FRAME & SYSTEM, vector<PAIRS> & FF_2BODY, vector<vector <XYZ > > & FRAME_A_MATRIX, const int nlayers, map<string,int> PAIR_MAP)	
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
	static double *rc;
	
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
	
	// Set up for layering

	int fidx_a2;
	int a2start, a2end, a2;
	int TOTAL_REP_ATOMS = SYSTEM.COORDS.size();

	for(int a1=0;a1<SYSTEM.ATOMS;a1++)		// Double sum over atom pairs
	{
		a2start = a1+1;
		a2end   = TOTAL_REP_ATOMS;

			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			a2 = a2idx;
			
			if(CONTROLS.N_LAYERS>0)
				if(SYSTEM.PARENT[a2] != -1)
					if(SYSTEM.PARENT[a2]<a1)
						continue;

			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];
		  
			//calculate vstart: (index for populating OO, OH, or HH column block of A).
		
			vstart = curr_pair_type_idx * FF_2BODY[curr_pair_type_idx].SNUM;

			// Get pair distance

			rlen = get_dist(SYSTEM, CONTROLS, RAB, a1, a2);	// Updates RAB!
						
			if ( (rlen < FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST) && (a1 !=a2) )
				FF_2BODY[curr_pair_type_idx].MIN_FOUND_DIST = rlen;

			if(rlen > FF_2BODY[curr_pair_type_idx].S_MINIM && rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)
			{
				// calculate binning, convert all distances to au from angstroms 
				x = rlen/autoang;
				
				fidx_a2 = a2;
				
				if(CONTROLS.N_LAYERS>0)
					if(SYSTEM.PARENT[a2] != -1)
							fidx_a2 = SYSTEM.PARENT[a2];

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
	double Vtot = 0;
	
	vector<XYZ> SForce(SYSTEM.ATOMS);	// This is the force ( negative of dEover ) (LEF)
	vector<XYZ> dEover(SYSTEM.ATOMS);	// Over bonding energy derivative (force) w/r/t delta ("S")
	vector<XYZ> dFover(SYSTEM.ATOMS);	// Over bonding force derivative w/r/t pover
	
	// THREE-BODY POTENTIAL:
	
	double Eover = 0;			// Overbonding energy
	double rik,tempr,temps;	
	double S[SYSTEM.ATOMS];				// Difference between oxygen bond order and valency (Called delta in the paper)

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
 

void ZCalc(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, NEIGHBORS & NEIGHBOR_LIST)
{
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
		ZCalc_Cheby(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, NEIGHBOR_LIST);
		  
		if(FF_2BODY[0].SNUM_3B_CHEBY > 0)
			ZCalc_3B_Cheby(SYSTEM, CONTROLS, FF_2BODY, FF_3BODY, PAIR_MAP, TRIAD_MAP, NEIGHBOR_LIST);	;
	}
	
	else if ( FF_2BODY[0].PAIRTYP == "LJ" ) 
		ZCalc_Lj(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP);
	
	else if ( FF_2BODY[0].PAIRTYP == "SPLINE" ) 
		ZCalc_Spline(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP, NEIGHBOR_LIST);	
	
    else 
    {
		cout << "Error: bad pairtype in ZCalc: " << FF_2BODY[0].PAIRTYP << endl;
		exit(1);
    }	
	
	if ( CONTROLS.USE_COULOMB ) 
		ZCalc_Ewald(SYSTEM, CONTROLS, NEIGHBOR_LIST);

	if ( CONTROLS.USE_OVERCOORD ) 
      ZCalcSR_Over(SYSTEM, CONTROLS, FF_2BODY, PAIR_MAP);
	


	// FUNCTIONS THAT NEED UPDATING:

	/* 

	else if ( pair_type == INVERSE_R ) 	
		ZCalc_SR_Analytic(Coord,Lbc, Latcons,nlayers,nat,smin,smax,snum, SForce,Vtot, Pxyz, params);

	else if ( pair_type == STILLINGER )  // SAVE THIS FOR SECOND TO LAST FOR SIMILAR REASONS  
		ZCalc_Stillinger(Coord,Lbc, Latcons,nlayers,nat,smax, SForce,Vtot,Pxyz);

	else 
		EXIT_MSG("Error: Unknown pair type", pair_type)
	*/

	SYSTEM.PRESSURE_XYZ /= 3.0 * SYSTEM.BOXDIM.X * SYSTEM.BOXDIM.Y * SYSTEM.BOXDIM.Z;

  return;
}

// UPDATED AND VERIFIED
static void ZCalc_Cheby(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
// Calculate short-range forces using a Chebyshev polynomial expansion. Can use morse variables similar to the work of Bowman.
{
	XYZ RVEC; 		// Replaces Rvec[3];
	XYZ RAB; 		// Replaces  Rab[3];
	double rlen;
	double x;
	double deriv;
	static double *Tn, *Tnd;
	static bool called_before = false;

	double exprlen;	
	double coeff;			  
	double fcut0; 
	double fcut; 
	double fcutderiv; 				
	double Vpenalty;

	// A penalty function is added to the potential for r + penalty_dist < smin[ipair]
	// All pairs have the same penalty scale and distance
	const double penalty_scale  = FF_2BODY[0].PENALTY_SCALE;	// 1.0e8;
	const double penalty_dist   = FF_2BODY[0].PENALTY_DIST;  	// 0.01;
	
	bool inverse_order;
	double dx_dr;
	
	static double CHEBY_DERIV_CONST;	// Accounts for when cheby range is changed from -1:1 to x:y
		
	const double fcut_power = 
    #ifndef FPENALTY_POWER
		3.0;
	#else
		FPENALTY_POWER;
    #endif

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
		
		CHEBY_DERIV_CONST = FF_2BODY[0].CHEBY_RANGE_HIGH - FF_2BODY[0].CHEBY_RANGE_LOW;	// Ranges should be the same for all types
		CHEBY_DERIV_CONST /= 2.0; 														// i.e the width of the default cheby range
	}
  
	// Main loop for Chebyshev terms:
	
	string 	TEMP_STR;
	int 	curr_pair_type_idx;
	XYZ 	TMP_BOX;
	
	// Set up for MPI
	
	int a1start, a1end, a2;	
	int a2start, a2end;
	int fidx;

	// Divide atoms on a per-processor basis.

	if ((CONTROLS.COMPARE_FORCE && (CONTROLS.N_LAYERS>0)) || (CONTROLS.SUBTRACT_FORCE && (CONTROLS.N_LAYERS>0)))
		divide_atoms(a1start, a1end, SYSTEM.ATOMS/pow(CONTROLS.N_LAYERS+1,3.0));
	else 
		divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.

	for(int a1=a1start; a1<=a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev -1)
	{	
		if(NEIGHBOR_LIST.USE)
		{
			a2start = 0;
			a2end   = NEIGHBOR_LIST.LIST[a1].size();
		}
		else
		{
			a2start = a1+1;
			a2end   = SYSTEM.ATOMS;
		}
			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			if(NEIGHBOR_LIST.USE)
				a2 = NEIGHBOR_LIST.LIST[a1][a2idx];
			else
				a2 = a2idx;
			
			// If a2 is NOT a primitive atom, exclude a1/a2 interactions
			// when a2's parent is larger than a1 (which must be primitive)
			
			if (CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE)
				if(CONTROLS.N_LAYERS>0)
					if(SYSTEM.PARENT[a2] != -1)
						if(SYSTEM.PARENT[a2]<a1)
							continue;
			
			rlen = get_dist(SYSTEM, CONTROLS, RAB, a1, a2);	// Updates RAB!

			TEMP_STR =      SYSTEM.ATOMTYPE.at(a1);
			TEMP_STR.append(SYSTEM.ATOMTYPE.at(a2));
						
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];			
			
			if(rlen < FF_2BODY[curr_pair_type_idx].S_MAXIM)	// We want to evaluate the penalty function when r < rmin (LEF)
			{		
				double xdiff, rpenalty;

				// Apply a penalty for distances less than smin + penalty_dist.
				
				if ( rlen - penalty_dist < FF_2BODY[curr_pair_type_idx].S_MINIM ) 
					rpenalty = FF_2BODY[curr_pair_type_idx].S_MINIM + penalty_dist - rlen;
				else 
					rpenalty = 0.0;							
				
				// Do the cheby distance transfomration
				
				cheby_var(rlen, 
							FF_2BODY[curr_pair_type_idx].S_MINIM,
							FF_2BODY[curr_pair_type_idx].S_MAXIM,
							FF_2BODY[curr_pair_type_idx].LAMBDA,
							FF_2BODY[curr_pair_type_idx].CHEBY_TYPE,
							x, xdiff, inverse_order, exprlen);
							 
				#if CHECK_CHEBY_RANGE == 1	

					// Make sure our newly transformed distance falls in defined range for Cheby polynomials
					if ( x < -1.0 || x > 1.0 )
					{
						x = fix_cheby_val(x, inverse_order );
					
						if (isatty(fileno(stdout)))
							cout << COUT_STYLE.BOLD << COUT_STYLE.MAGENTA << "Warning: (Step " << CONTROLS.STEP << ") In 2B Cheby transformation, r outside of allowed range. " << TEMP_STR << COUT_STYLE.ENDSTYLE << endl;
						else
							cout << "Warning: (Step " << CONTROLS.STEP << ") In 2B Cheby transformation, r outside of allowed range. " << TEMP_STR << endl;
					}	
				
					// Now change the range, if the user requested

					x = x*CHEBY_DERIV_CONST + FF_2BODY[0].CHEBY_RANGE_LOW - -1.0*CHEBY_DERIV_CONST;

					// Sanity check

					if ( x < FF_2BODY[0].CHEBY_RANGE_LOW || x > FF_2BODY[0].CHEBY_RANGE_HIGH )
					{
						cout << "ERROR: transformed x falls outside user-defined range." << endl;
						cout << "x: " << x << endl;
						cout << "high/low: " << FF_2BODY[0].CHEBY_RANGE_HIGH << " " << FF_2BODY[0].CHEBY_RANGE_LOW  << endl;
						exit_run(0);
					}
				#endif						
				
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

				// Smoothing function to force PES to zero at outer cutoff

				fcut0 = (1.0 - rlen/FF_2BODY[curr_pair_type_idx].S_MAXIM);
				fcut      = pow(fcut0, fcut_power);
				fcutderiv = pow(fcut0,fcut_power-1);
				
				fcutderiv *= -1.0 * fcut_power / FF_2BODY[curr_pair_type_idx].S_MAXIM;

				dx_dr = CHEBY_DERIV_CONST*cheby_var_deriv(xdiff, rlen, FF_2BODY[curr_pair_type_idx].LAMBDA, FF_2BODY[curr_pair_type_idx].CHEBY_TYPE);

				fidx = a2;
				
				// When doing a compare force calculation, make sure that forces on 
				// replicate atoms are attributed to the parent atoms
				
				if (CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE)
					if(CONTROLS.N_LAYERS>0)
						if(SYSTEM.PARENT[a2] != -1)
								fidx = SYSTEM.PARENT[a2];

				for ( int i = 0; i < FF_2BODY[curr_pair_type_idx].SNUM; i++ ) 
				{
					coeff                = FF_2BODY[curr_pair_type_idx].PARAMS[i]; // This is the Cheby FF param for the given power
					SYSTEM.TOT_POT_ENER += coeff * fcut * Tn[i+1];

					deriv                = (fcut * Tnd[i+1] * dx_dr + fcutderiv * Tn[i+1]);
					SYSTEM.PRESSURE_XYZ -= coeff * deriv * rlen;

					SYSTEM.ACCEL[a1].X += coeff * deriv * RAB.X / rlen;
					SYSTEM.ACCEL[a1].Y += coeff * deriv * RAB.Y / rlen;
					SYSTEM.ACCEL[a1].Z += coeff * deriv * RAB.Z / rlen;
					
					SYSTEM.ACCEL[fidx].X -= coeff * deriv * RAB.X / rlen;
					SYSTEM.ACCEL[fidx].Y -= coeff * deriv * RAB.Y / rlen;
					SYSTEM.ACCEL[fidx].Z -= coeff * deriv * RAB.Z / rlen;

				}
				// Add penalty for very short distances, where the fit FF may be unphysical (preserve conservation of E).

				if ( rpenalty > 0.0 ) 
				{
					Vpenalty = 0.0;
					
					if (isatty(fileno(stdout)))
						cout << COUT_STYLE.BOLD << COUT_STYLE.MAGENTA << "Warning: (Step " << CONTROLS.STEP << ")Adding penalty in 2B Cheby calc, r < rmin+penalty_dist " << rlen << " " << FF_2BODY[curr_pair_type_idx].S_MINIM+penalty_dist << " " << TEMP_STR << COUT_STYLE.ENDSTYLE << endl;
					else
						cout << "Warning: (Step " << CONTROLS.STEP << ") Adding penalty in 2B Cheby calc, r < rmin+penalty_dist " << rlen << " " << FF_2BODY[curr_pair_type_idx].S_MINIM+penalty_dist << " " << TEMP_STR << endl;
					

							
					SYSTEM.ACCEL[fidx].X += 3.0 * rpenalty * rpenalty * penalty_scale * RAB.X / rlen;
					SYSTEM.ACCEL[fidx].Y += 3.0 * rpenalty * rpenalty * penalty_scale * RAB.Y / rlen;
					SYSTEM.ACCEL[fidx].Z += 3.0 * rpenalty * rpenalty * penalty_scale * RAB.Z / rlen;
					
					SYSTEM.ACCEL[a1].X -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB.X / rlen;
					SYSTEM.ACCEL[a1].Y -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB.Y / rlen;
					SYSTEM.ACCEL[a1].Z -= 3.0 * rpenalty * rpenalty * penalty_scale * RAB.Z / rlen;								
					
					Vpenalty = rpenalty * rpenalty * rpenalty * penalty_scale;
					SYSTEM.TOT_POT_ENER += Vpenalty;
					cout << "	...Penalty potential = "<< Vpenalty << endl;
				}						
			} 
		}
	}
	return;
} 

// UPDATED AND VERIFIED
static void ZCalc_3B_Cheby(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, vector<TRIP_FF> & FF_3BODY, map<string,int> & PAIR_MAP, map<string,int> & TRIAD_MAP, NEIGHBORS & NEIGHBOR_LIST)
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
	
	static double *Tn_ij,  *Tn_ik,  *Tn_jk;
	static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
	static bool called_before = false;
	
	static int pow_ij, pow_ik, pow_jk;
			  
	double fcut0_ij, fcut0_ik, fcut0_jk; 
	double fcut_ij,  fcut_ik,  fcut_jk; 			
	double deriv_ij, deriv_ik, deriv_jk;
	double force_ij, force_ik, force_jk;
	double xdiff_ij, xdiff_ik, xdiff_jk;
	double fcutderiv_ij, fcutderiv_ik, fcutderiv_jk; 	
	double tempx; 	
	
	static string TEMP_STR;
	static string PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK;
	static int curr_triple_type_index;
	static int curr_pair_type_idx_ij;
	static int curr_pair_type_idx_ik;
	static int curr_pair_type_idx_jk;
	
	int fidx_a2, fidx_a3;
	
	double dx_dr_ij, dx_dr_ik, dx_dr_jk;	
	
	static double CHEBY_DERIV_CONST;	// Accounts for when cheby range is changed from -1:1 to x:y
		
	double coeff;
	
	double S_MAXIM_IJ, S_MAXIM_IK, S_MAXIM_JK;
	double S_MINIM_IJ, S_MINIM_IK, S_MINIM_JK;
	
	double DUMMY_RIJ;
	double DUMMY_RIK;
	double DUMMY_RJK;

	double SWITCH_IJ;
	double SWITCH_IK;
	double SWITCH_JK;
	
	double DSWITCH_IJ;
	double DSWITCH_IK;
	double DSWITCH_JK;

	double SWITCH_SMOOTHNESS = 2.0;
	
	bool PENALTY_KICKIN = false;

	// A penalty function is added to the potential for r + penalty_dist < smin[ipair]
	// All pairs have the same penalty scale and distance
	const double penalty_scale  = FF_2BODY[0].PENALTY_SCALE;	// 1.0e8;
	const double penalty_dist   = FF_2BODY[0].PENALTY_DIST;  	// 0.01;

	const double fcut_power = 
    #ifndef FPENALTY_POWER
		3.0;
	#else
		FPENALTY_POWER;
    #endif

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
			
			CHEBY_DERIV_CONST = FF_2BODY[0].CHEBY_RANGE_HIGH - FF_2BODY[0].CHEBY_RANGE_LOW;	// Ranges should be the same for all types
			CHEBY_DERIV_CONST /= 2.0; // i.e the width of the default cheby range
		
	}
	
	tempx = 0;
  
	// Main loop for Chebyshev terms:


	// Set up for MPI
	
	int a1start, a1end;	

	if ((CONTROLS.COMPARE_FORCE && (CONTROLS.N_LAYERS>0)) || (CONTROLS.SUBTRACT_FORCE && (CONTROLS.N_LAYERS>0)))
		divide_atoms(a1start, a1end, SYSTEM.ATOMS/pow(CONTROLS.N_LAYERS+1,3.0));
	else 
		divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.

	// Set up for neighbor lists
	
	int a2start, a2end, a2;
	int a3start, a3end, a3;

	for(int a1=a1start; a1<=a1end; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev -1)
	{		
		if(NEIGHBOR_LIST.USE)
		{
			a2start = 0;
			a2end   = NEIGHBOR_LIST.LIST[a1].size();
		}
		else
		{
			a2start = a1+1;
			a2end   = SYSTEM.ATOMS;
		}
			
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{			
			if(NEIGHBOR_LIST.USE)
				a2 = NEIGHBOR_LIST.LIST[a1][a2idx];
			else
				a2 = a2idx;

			
			if(NEIGHBOR_LIST.USE)
			{
				a3start = 0;
				a3end   = NEIGHBOR_LIST.LIST[a2].size();
			}
			else
			{
				a3start = a2+1;
				a3end   = SYSTEM.ATOMS;
			}
			
			if (CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE)
				if(CONTROLS.N_LAYERS>0)
					if(SYSTEM.PARENT[a2] != -1)
						if(SYSTEM.PARENT[a2]<a1)
							continue;
			
			for(int a3idx=a3start; a3idx<a3end; a3idx++)	
			{			
				if(NEIGHBOR_LIST.USE)
					a3 = NEIGHBOR_LIST.LIST[a2][a3idx];
				else
					a3 = a3idx;
				
				if (CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE)
					if(CONTROLS.N_LAYERS>0)
						if(SYSTEM.PARENT[a3] != -1)
							if(SYSTEM.PARENT[a3]<SYSTEM.PARENT[a2])
								continue;

				
				TEMP_STR = SYSTEM.ATOMTYPE[a1];
				TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);		
				PAIR_TYPE_IJ = TEMP_STR;					
				curr_pair_type_idx_ij = PAIR_MAP[TEMP_STR];
	
				TEMP_STR = SYSTEM.ATOMTYPE[a1];
				TEMP_STR.append(SYSTEM.ATOMTYPE[a3]);	
				PAIR_TYPE_IK = TEMP_STR;							
				curr_pair_type_idx_ik = PAIR_MAP[TEMP_STR];	
			
				TEMP_STR = SYSTEM.ATOMTYPE[a2];
				TEMP_STR.append(SYSTEM.ATOMTYPE[a3]);	
				PAIR_TYPE_JK = TEMP_STR;							
				curr_pair_type_idx_jk = PAIR_MAP[TEMP_STR];		
		
				// Determine the FF type for the given triplet

				TEMP_STR =      FF_2BODY[curr_pair_type_idx_ij].PRPR_NM;
				TEMP_STR.append(FF_2BODY[curr_pair_type_idx_ik].PRPR_NM);	
				TEMP_STR.append(FF_2BODY[curr_pair_type_idx_jk].PRPR_NM);	

				curr_triple_type_index = TRIAD_MAP[TEMP_STR];	

				if(curr_triple_type_index<0)
					continue;
				
									
				rlen_ij = get_dist(SYSTEM, CONTROLS, RAB_IJ, a1, a2);	// Updates RAB!
				rlen_ik = get_dist(SYSTEM, CONTROLS, RAB_IK, a1, a3);	// Updates RAB!
				rlen_jk = get_dist(SYSTEM, CONTROLS, RAB_JK, a2, a3);	// Updates RAB!

				S_MAXIM_IJ = SET_SMAXIM(FF_2BODY[curr_pair_type_idx_ij], FF_3BODY[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ij].PRPR_NM);
				S_MAXIM_IK = SET_SMAXIM(FF_2BODY[curr_pair_type_idx_ik], FF_3BODY[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ik].PRPR_NM);
				S_MAXIM_JK = SET_SMAXIM(FF_2BODY[curr_pair_type_idx_jk], FF_3BODY[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_jk].PRPR_NM);
				
				S_MINIM_IJ = SET_SMINIM(FF_2BODY[curr_pair_type_idx_ij], FF_3BODY[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ij].PRPR_NM);
				S_MINIM_IK = SET_SMINIM(FF_2BODY[curr_pair_type_idx_ik], FF_3BODY[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_ik].PRPR_NM);
				S_MINIM_JK = SET_SMINIM(FF_2BODY[curr_pair_type_idx_jk], FF_3BODY[curr_triple_type_index],FF_2BODY[curr_pair_type_idx_jk].PRPR_NM);


				// Before doing any polynomial/coeff set up, make sure that all ij, ik, and jk distances are 
				// within the allowed range.

				if((rlen_ij < S_MAXIM_IJ) && (rlen_ij > FF_2BODY[curr_pair_type_idx_ij].S_MINIM))
				{			
					if((rlen_ik < S_MAXIM_IK) && (rlen_ik > FF_2BODY[curr_pair_type_idx_ik].S_MINIM))	
					{
						if((rlen_jk < S_MAXIM_JK) && (rlen_jk > FF_2BODY[curr_pair_type_idx_jk].S_MINIM))
						{
							
							PENALTY_KICKIN = false;
							
							// Set up variables for switching
							
							DUMMY_RIJ = rlen_ij;
							DUMMY_RIK = rlen_ik;
							DUMMY_RJK = rlen_jk;
							
							SWITCH_IJ = 0.0;
							SWITCH_IK = 0.0;
							SWITCH_JK = 0.0;
		
							DSWITCH_IJ = 0.0;
							DSWITCH_IK = 0.0;
							DSWITCH_JK = 0.0;
		
							
							if(rlen_ij<S_MINIM_IJ)
								DUMMY_RIJ = S_MINIM_IJ;
							
							if(rlen_ik<S_MINIM_IK)
								DUMMY_RIK = S_MINIM_IK;
							
							if(rlen_jk<S_MINIM_JK)
								DUMMY_RJK = S_MINIM_JK;

							if (rlen_ij < S_MINIM_IJ )
							{
								PENALTY_KICKIN = true;
								SWITCH_IJ  = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_IJ - rlen_ij));
								
								DSWITCH_IJ = -1.0 * penalty_scale * SWITCH_SMOOTHNESS * pow(cosh(SWITCH_SMOOTHNESS*(S_MINIM_IJ - rlen_ij)),-2.0);
								DSWITCH_IJ /= rlen_ij;

								if(rlen_ij < S_MINIM_IJ)
									DUMMY_RIJ =  S_MINIM_IJ;
							}
							if(rlen_ik < S_MINIM_IK )
							{
								PENALTY_KICKIN = true;
								
								SWITCH_IK  = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_IK - rlen_ik));
								
								DSWITCH_IK = -1.0 * penalty_scale * SWITCH_SMOOTHNESS * pow(cosh(SWITCH_SMOOTHNESS*(S_MINIM_IK - rlen_ik)),-2.0); 
								DSWITCH_IK /= rlen_ik;

								if(rlen_ik < S_MINIM_IK)
									DUMMY_RIK =  S_MINIM_IK;
							}
							if (rlen_jk < S_MINIM_JK )
							{
								PENALTY_KICKIN = true;
								
								SWITCH_JK  = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_JK - rlen_jk));
								
								DSWITCH_JK = -1.0 * penalty_scale * SWITCH_SMOOTHNESS * pow(cosh(SWITCH_SMOOTHNESS*(S_MINIM_JK - rlen_jk)),-2.0);
								DSWITCH_JK /= rlen_jk;

								if(rlen_jk < S_MINIM_JK)
									DUMMY_RJK =  S_MINIM_JK;
							}
			
							
							// Everything is within allowed ranges. Begin setting up the force calculation

							// Set up the polynomials
			
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_ij], Tn_ij, Tnd_ij, DUMMY_RIJ, xdiff_ij, S_MAXIM_IJ, S_MINIM_IJ);
			
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_ik], Tn_ik, Tnd_ik, DUMMY_RIK, xdiff_ik, S_MAXIM_IK, S_MINIM_IK);
			
							SET_3B_CHEBY_POLYS(FF_2BODY[curr_pair_type_idx_jk], Tn_jk, Tnd_jk, DUMMY_RJK, xdiff_jk, S_MAXIM_JK, S_MINIM_JK);
			
							// Apply the FF

							// Set up the penalty functions
			
							fcut0_ij     = (1.0 - rlen_ij/S_MAXIM_IJ);
							fcut_ij      = pow(fcut0_ij, fcut_power);
							fcutderiv_ij = pow(fcut0_ij, fcut_power-1);
							fcutderiv_ij *= -1.0 * fcut_power / S_MAXIM_IJ;		
			
							fcut0_ik     = (1.0 - rlen_ik/S_MAXIM_IK);
							fcut_ik      = pow(fcut0_ik, fcut_power);
							fcutderiv_ik = pow(fcut0_ik, fcut_power-1);
							fcutderiv_ik *= -1.0 * fcut_power / S_MAXIM_IK;	
			
							fcut0_jk     = (1.0 - rlen_jk/S_MAXIM_JK);
							fcut_jk      = pow(fcut0_jk, fcut_power);
							fcutderiv_jk = pow(fcut0_jk, fcut_power-1);
							fcutderiv_jk *= -1.0 * fcut_power / S_MAXIM_JK;	
			
							// Set up terms for derivatives
			
							dx_dr_ij = CHEBY_DERIV_CONST*cheby_var_deriv(xdiff_ij, DUMMY_RIJ, FF_2BODY[curr_pair_type_idx_ij].LAMBDA, FF_2BODY[curr_pair_type_idx_ij].CHEBY_TYPE);
							
							dx_dr_ik = CHEBY_DERIV_CONST*cheby_var_deriv(xdiff_ik, DUMMY_RIK, FF_2BODY[curr_pair_type_idx_ik].LAMBDA, FF_2BODY[curr_pair_type_idx_ik].CHEBY_TYPE);

							dx_dr_jk = CHEBY_DERIV_CONST*cheby_var_deriv(xdiff_jk, DUMMY_RJK, FF_2BODY[curr_pair_type_idx_jk].LAMBDA, FF_2BODY[curr_pair_type_idx_jk].CHEBY_TYPE);

							// Error check: Certain triplets are impossible...
			
							if(curr_triple_type_index < 0)
							{
								cout << "ERROR: Impossible atom triplet found: " << TEMP_STR << endl;
								cout << "       Check the parameter file." << endl;
								exit_run(0);
							}
			
							// Now compute the forces for each set of allowed powers for pairs ij, ik, and jk		
							// Keep in mind that the order in which allowed powers are stored may not match the
							// ordering of pairs resulting from the present atom triplet. Thus, we need to order
							// the stored powers properly before applying the FF.
							
							// When doing a compare force calculation, make sure that forces on 
							// replicate atoms are attributed to the parent atoms
			
							fidx_a2 = a2;
							fidx_a3 = a3;
							
							if (CONTROLS.COMPARE_FORCE)
							{
								if(CONTROLS.N_LAYERS>0)
								{
									if(SYSTEM.PARENT[a2] != -1)
											fidx_a2 = SYSTEM.PARENT[a2];
									
									if(SYSTEM.PARENT[a3] != -1)
											fidx_a3 = SYSTEM.PARENT[a3];
								}
							}
			
							for(int i=0; i<FF_3BODY[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
							{
								SET_3B_CHEBY_POWERS(FF_2BODY, FF_3BODY[curr_triple_type_index], PAIR_MAP,  pow_ij, pow_ik, pow_jk, PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK, i);
				
								coeff = FF_3BODY[curr_triple_type_index].PARAMS[i];
				
								tempx += coeff * fcut_ij * fcut_ik * fcut_jk * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk]; 

								deriv_ij  = fcut_ij * Tnd_ij[pow_ij] * dx_dr_ij + fcutderiv_ij * Tn_ij [pow_ij];
								deriv_ik  = fcut_ik * Tnd_ik[pow_ik] * dx_dr_ik + fcutderiv_ik * Tn_ik [pow_ik];
								deriv_jk  = fcut_jk * Tnd_jk[pow_jk] * dx_dr_jk + fcutderiv_jk * Tn_jk [pow_jk];
							
								force_ij  = coeff * deriv_ij * fcut_ik * fcut_jk * Tn_ik [pow_ik] * Tn_jk [pow_jk];
								force_ik  = coeff * deriv_ik * fcut_ij * fcut_jk * Tn_ij [pow_ij] * Tn_jk [pow_jk];
								force_jk  = coeff * deriv_jk * fcut_ij * fcut_ik * Tn_ij [pow_ij] * Tn_ik [pow_ik];
								
								SYSTEM.PRESSURE_XYZ    -= force_ij * rlen_ij;
								SYSTEM.PRESSURE_XYZ    -= force_ik * rlen_ik;
								SYSTEM.PRESSURE_XYZ    -= force_jk * rlen_jk;
								
								force_ij /= rlen_ij;
								force_ik /= rlen_ik;
								force_jk /= rlen_jk;

								// Apply forces to ij pair
				
								SYSTEM.ACCEL[a1]     .X += force_ij * RAB_IJ.X;
								SYSTEM.ACCEL[a1]     .Y += force_ij * RAB_IJ.Y;
								SYSTEM.ACCEL[a1]     .Z += force_ij * RAB_IJ.Z;
				
								SYSTEM.ACCEL[fidx_a2].X -= force_ij * RAB_IJ.X;
								SYSTEM.ACCEL[fidx_a2].Y -= force_ij * RAB_IJ.Y;
								SYSTEM.ACCEL[fidx_a2].Z -= force_ij * RAB_IJ.Z;
				
								// Apply forces to ik pair
				
								SYSTEM.ACCEL[a1]     .X += force_ik * RAB_IK.X;
								SYSTEM.ACCEL[a1]     .Y += force_ik * RAB_IK.Y;
								SYSTEM.ACCEL[a1]     .Z += force_ik * RAB_IK.Z;	
				
								SYSTEM.ACCEL[fidx_a3].X -= force_ik * RAB_IK.X;
								SYSTEM.ACCEL[fidx_a3].Y -= force_ik * RAB_IK.Y;
								SYSTEM.ACCEL[fidx_a3].Z -= force_ik  * RAB_IK.Z;	
				
								// Apply forces to jk pair
				
								SYSTEM.ACCEL[fidx_a2].X += force_jk * RAB_JK.X;
								SYSTEM.ACCEL[fidx_a2].Y += force_jk * RAB_JK.Y;
								SYSTEM.ACCEL[fidx_a2].Z += force_jk * RAB_JK.Z;
				
								SYSTEM.ACCEL[a3]     .X -= force_jk * RAB_JK.X;
								SYSTEM.ACCEL[a3]     .Y -= force_jk * RAB_JK.Y;
								SYSTEM.ACCEL[a3]     .Z -= force_jk * RAB_JK.Z;			
								
								
								
/*								
								///////////
								PE_NORMAL += tempx;

								F_NORMAL.X = force_jk * RAB_JK.X;
								F_NORMAL.Y = force_jk * RAB_JK.Y;
								F_NORMAL.Z = force_jk * RAB_JK.Z;
*/
															
				
								#if FORCECHECK
				
									// Apply forces to ij pair
				
									FORCE_3B[a1]     .X += force_ij * RAB_IJ.X;
									FORCE_3B[a1]     .Y += force_ij * RAB_IJ.Y;
									FORCE_3B[a1]     .Z += force_ij * RAB_IJ.Z;
				
									FORCE_3B[fidx_a2].X -= force_ij * RAB_IJ.X;
									FORCE_3B[fidx_a2].Y -= force_ij * RAB_IJ.Y;
									FORCE_3B[fidx_a2].Z -= force_ij * RAB_IJ.Z;
				
									// Apply forces to ik pair
				
									FORCE_3B[a1]     .X += force_ik * RAB_IK.X;
									FORCE_3B[a1]     .Y += force_ik * RAB_IK.Y;
									FORCE_3B[a1]     .Z += force_ik  * RAB_IK.Z;	
				
									FORCE_3B[fidx_a3].X -= force_ik * RAB_IK.X;
									FORCE_3B[fidx_a3].Y -= force_ik * RAB_IK.Y;
									FORCE_3B[fidx_a3].Z -= force_ik * RAB_IK.Z;	
				
									// Apply forces to jk pair
				
									FORCE_3B[fidx_a2].X += force_jk * RAB_JK.X;
									FORCE_3B[fidx_a2].Y += force_jk * RAB_JK.Y;
									FORCE_3B[fidx_a2].Z += force_jk * RAB_JK.Z;
				
									FORCE_3B[fidx_a3].X -= force_jk * RAB_JK.X;
									FORCE_3B[fidx_a3].Y -= force_jk * RAB_JK.Y;
									FORCE_3B[fidx_a3].Z -= force_jk * RAB_JK.Z;										
				
								#endif
											
							}	
							
							if(PENALTY_KICKIN)
							{
								tempx += SWITCH_IJ + SWITCH_IK + SWITCH_JK;
								
								// Apply forces to ij pair
				
								SYSTEM.ACCEL[a1]     .X += DSWITCH_IJ * RAB_IJ.X;
								SYSTEM.ACCEL[a1]     .Y += DSWITCH_IJ * RAB_IJ.Y;
								SYSTEM.ACCEL[a1]     .Z += DSWITCH_IJ * RAB_IJ.Z;
				
								SYSTEM.ACCEL[fidx_a2].X -= DSWITCH_IJ * RAB_IJ.X;
								SYSTEM.ACCEL[fidx_a2].Y -= DSWITCH_IJ * RAB_IJ.Y;
								SYSTEM.ACCEL[fidx_a2].Z -= DSWITCH_IJ * RAB_IJ.Z;
				
								// Apply forces to ik pair
				
								SYSTEM.ACCEL[a1]     .X += DSWITCH_IK * RAB_IK.X;
								SYSTEM.ACCEL[a1]     .Y += DSWITCH_IK * RAB_IK.Y;
								SYSTEM.ACCEL[a1]     .Z += DSWITCH_IK * RAB_IK.Z;	
				
								SYSTEM.ACCEL[fidx_a3].X -= DSWITCH_IK * RAB_IK.X;
								SYSTEM.ACCEL[fidx_a3].Y -= DSWITCH_IK * RAB_IK.Y;
								SYSTEM.ACCEL[fidx_a3].Z -= DSWITCH_IK * RAB_IK.Z;	
				
								// Apply forces to jk pair
				
								SYSTEM.ACCEL[fidx_a2].X += DSWITCH_JK * RAB_JK.X;
								SYSTEM.ACCEL[fidx_a2].Y += DSWITCH_JK * RAB_JK.Y;
								SYSTEM.ACCEL[fidx_a2].Z += DSWITCH_JK * RAB_JK.Z;
				
								SYSTEM.ACCEL[a3]     .X -= DSWITCH_JK * RAB_JK.X;
								SYSTEM.ACCEL[a3]     .Y -= DSWITCH_JK * RAB_JK.Y;
								SYSTEM.ACCEL[a3]     .Z -= DSWITCH_JK * RAB_JK.Z;	
								
								SYSTEM.TOT_POT_ENER += SWITCH_IJ + SWITCH_IK + SWITCH_JK;	

								#if FORCECHECK
				
									// Apply forces to ij pair
				
									FORCE_3B[a1]     .X += DSWITCH_IJ * RAB_IJ.X;
									FORCE_3B[a1]     .Y += DSWITCH_IJ * RAB_IJ.Y;
									FORCE_3B[a1]     .Z += DSWITCH_IJ * RAB_IJ.Z;
				
									FORCE_3B[fidx_a2].X -= DSWITCH_IJ * RAB_IJ.X;
									FORCE_3B[fidx_a2].Y -= DSWITCH_IJ * RAB_IJ.Y;
									FORCE_3B[fidx_a2].Z -= DSWITCH_IJ * RAB_IJ.Z;
				
									// Apply forces to ik pair
				
									FORCE_3B[a1]     .X += DSWITCH_IK * RAB_IK.X;
									FORCE_3B[a1]     .Y += DSWITCH_IK * RAB_IK.Y;
									FORCE_3B[a1]     .Z += DSWITCH_IK * RAB_IK.Z;	
				
									FORCE_3B[fidx_a3].X -= DSWITCH_IK * RAB_IK.X;
									FORCE_3B[fidx_a3].Y -= DSWITCH_IK * RAB_IK.Y;
									FORCE_3B[fidx_a3].Z -= DSWITCH_IK * RAB_IK.Z;	
				
									// Apply forces to jk pair
				
									FORCE_3B[fidx_a2].X += DSWITCH_JK * RAB_JK.X;
									FORCE_3B[fidx_a2].Y += DSWITCH_JK * RAB_JK.Y;
									FORCE_3B[fidx_a2].Z += DSWITCH_JK * RAB_JK.Z;
				
									FORCE_3B[fidx_a3].X -= DSWITCH_JK * RAB_JK.X;
									FORCE_3B[fidx_a3].Y -= DSWITCH_JK * RAB_JK.Y;
									FORCE_3B[fidx_a3].Z -= DSWITCH_JK * RAB_JK.Z;										
				
								#endif
								
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
	
	// Set up for MPI
	
	int a1start, a1end ;
	divide_atoms(a1start, a1end, SYSTEM.ATOMS) ;	// Divide atoms on a per-processor basis.

	for(int a1=a1start ;a1 <= a1end ; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS 
	{
		for(int a2=0;a2<a1;a2++)
		{			
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];

			// pair interaction cutoff distance.
			double rcutoff = FF_2BODY[curr_pair_type_idx].S_MAXIM ;

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
			
			if ( rlen_mi < FF_2BODY[curr_pair_type_idx].PARAMS[1]/2.2) 
				{				
					EXIT_MSG("Error: close approach", rlen_mi);
				}
			else if ( rlen_mi < rcutoff ) {
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
}

// UPDATED AND VERIFIED
static void ZCalc_Spline(FRAME & SYSTEM, MD_JOB_CONTROL & CONTROLS, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, NEIGHBORS & NEIGHBOR_LIST)
// Calculate spline forces.
{
	XYZ RVEC, RAB;
	double rlen,rlen2;
	double tempx;
	double S_r;
	
	int k0;
	double x, x0;
	double t, t2, t3, t4;;
	double h00,h10,h01,h11;
	double i00,i10,i01,i11;
	int kstart;
	
	string TEMP_STR;
	int curr_pair_type_idx;

	// Set up for MPI
	
	int a1start, a1end, fidx;
	
	if ((CONTROLS.COMPARE_FORCE && (CONTROLS.N_LAYERS>0)) || (CONTROLS.SUBTRACT_FORCE && (CONTROLS.N_LAYERS>0)))
		divide_atoms(a1start, a1end, SYSTEM.ATOMS/pow(CONTROLS.N_LAYERS+1,3.0));
	else 
		divide_atoms(a1start, a1end, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.

	// Set up for neighbor list
	
	int a2start, a2end, a2;

	for(int a1=a1start ;a1 <= a1end ; a1++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS (prev-1) 
	{
		if(NEIGHBOR_LIST.USE)
		{
			a2start = 0;
			a2end   = NEIGHBOR_LIST.LIST[a1].size();
		}
		else
		{
			a2start = a1+1;
			a2end   = SYSTEM.ATOMS;
		}
		
		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{	
			if(NEIGHBOR_LIST.USE)
				a2 = NEIGHBOR_LIST.LIST[a1][a2idx];
			else
				a2 = a2idx;
			
			// If a2 is NOT a primitive atom, exclude a1/a2 interactions
			// when a2's parent is larger than a1 (which must be primitive)
			
			if (CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE)
				if(CONTROLS.N_LAYERS>0)
					if(SYSTEM.PARENT[a2] != -1)
						if(SYSTEM.PARENT[a2]<a1)
							continue;

			rlen = get_dist(SYSTEM, CONTROLS, RAB, a1, a2);	// Updates RAB!
			
			TEMP_STR = SYSTEM.ATOMTYPE[a1];
			TEMP_STR.append(SYSTEM.ATOMTYPE[a2]);
							
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];

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
					exit_run(1);
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
				
				fidx = a2;
				
				// When doing a compare force calculation, make sure that forces on 
				// replicate atoms are attributed to the parent atoms
				
				if (CONTROLS.COMPARE_FORCE || CONTROLS.SUBTRACT_FORCE)
					if(CONTROLS.N_LAYERS>0)
						if(SYSTEM.PARENT[a2] != -1)
								fidx = SYSTEM.PARENT[a2];

				SYSTEM.ACCEL[a1].X += S_r * RAB.X/rlen;
				SYSTEM.ACCEL[a1].Y += S_r * RAB.Y/rlen;
				SYSTEM.ACCEL[a1].Z += S_r * RAB.Z/rlen;

				SYSTEM.ACCEL[fidx].X -= S_r * RAB.X/rlen;
				SYSTEM.ACCEL[fidx].Y -= S_r * RAB.Y/rlen;
				SYSTEM.ACCEL[fidx].Z -= S_r * RAB.Z/rlen;
	   
				SYSTEM.PRESSURE_XYZ -= S_r * rlen;
				SYSTEM.TOT_POT_ENER += tempx;

			}//rlen
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
	double rik,tempr,temps;
	double powrik;
	double Sexp2;
	
	vector<double> S     (SYSTEM.ATOMS);
	vector<double> Sexp  (SYSTEM.ATOMS);
	vector<XYZ>    dEover(SYSTEM.ATOMS);
	
	int curr_pair_type_idx; 
	string TEMP_STR;
	
	bool SAFE = false;

	// Set up for MPI
	
	int aistart, aiend;	
	divide_atoms(aistart, aiend, SYSTEM.ATOMS);	// Divide atoms on a per-processor basis.

	for(int ai=aistart;ai <= aiend; ai++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS
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

	for(int ai=aistart;ai <= aiend; ai++)		// Double sum over atom pairs -- MPI'd over SYSTEM.ATOMS
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

void Print_Cheby(vector<PAIR_FF> & FF_2BODY, int ij, string PAIR_NAME, bool INCLUDE_FCUT, bool INCLUDE_CHARGES, string FILE_TAG)
// Generating pair distance scans for the 2-b potential.
// pair distances will range from smin to smax, incremented by sdelta
{
	
//	string SCAN_FILE_2B;	// Declared as a global (external) variable in functions.h
	
	double rlen;
	double x;
	static double *Tn;
	static bool called_before = false;
	
	double xdiff,rpenalty;
	double exprlen;	
	double coeff;			  
	double fcut0; 
	double fcut; 				
	double Vpenalty;
	
	bool inverse_order;

	static double CHEBY_DERIV_CONST;	// Accounts for when cheby range is changed from -1:1 to x:y
	
	string OUTFILE = "2b_Cheby_Pot-";
	OUTFILE.append(PAIR_NAME);
	
	if(FILE_TAG != "")
		OUTFILE.append("_for_3B");
	
	OUTFILE.append(".dat");
	ofstream OUTFILE_2B_POT;
	OUTFILE_2B_POT.open(OUTFILE.data());	
	
	SCAN_FILE_2B = OUTFILE;
	
	// A penalty function is added to the potential for r + penalty_dist < smin[ipair]
	// All pairs have the same penalty scale and distance
	const double penalty_scale  = FF_2BODY[0].PENALTY_SCALE;	// 1.0e8;
	const double penalty_dist   = FF_2BODY[0].PENALTY_DIST;  	// 0.01;

	const double fcut_power = 
    #ifndef FPENALTY_POWER
		3.0;
	#else
		FPENALTY_POWER;
    #endif


	if ( ! called_before ) 
	{
		called_before = true;
		int dim = 0;
		
		for ( int i = 0; i < FF_2BODY.size(); i++ ) 
			if (FF_2BODY[i].SNUM > dim ) 
				dim = FF_2BODY[i].SNUM;	 
		
		dim++;
		Tn   = new double [dim];
		
		CHEBY_DERIV_CONST = FF_2BODY[0].CHEBY_RANGE_HIGH - FF_2BODY[0].CHEBY_RANGE_LOW;	// Ranges should be the same for all types
		CHEBY_DERIV_CONST /= 2.0; // i.e the width of the default cheby range
	}
  
	// Main loop for Chebyshev terms:
	
	int n_ij = (FF_2BODY[ij].S_MAXIM - FF_2BODY[ij].S_MINIM)/FF_2BODY[ij].S_DELTA;
	
	double tempx;

	for (int a=1; a<n_ij; a++)
	{
		tempx = 0;
		
		double SWITCH, DSWITCH;
		double SWITCH_SMOOTHNESS = 0.1;
	
		SWITCH  = 1.0;
		DSWITCH = 1.0;
		
		rlen = FF_2BODY[ij].S_MINIM + a * FF_2BODY[ij].S_DELTA;
		
		// Check if this is for a 2+3B scan, and whether 2B has been extrapolated from old values
		if (FF_2BODY[ij].S_MINIM != FF_2BODY[ij].OLD_S_MINIM && (FILE_TAG != "") )
		{
			//SWITCH_IJ  = 1.0 - 0.5 + 0.5 * tanh( (RLEN_IJ-S_MINIM_IJ_3B)/SWITCH_SMOOTHNESS );
			SWITCH  =  0.5 + 0.5 * tanh( (rlen-FF_2BODY[ij].OLD_S_MINIM)/SWITCH_SMOOTHNESS );
		}

		if(rlen > FF_2BODY[ij].S_MINIM and rlen < FF_2BODY[ij].S_MAXIM)
		{
			// Apply a penalty for distances less than smin - penalty_dist.
			
			rpenalty = 0.0;
			
			if ( rlen - penalty_dist < FF_2BODY[ij].S_MINIM ) 
				rpenalty = FF_2BODY[ij].S_MINIM + penalty_dist - rlen;

			// Do the distance transfomration

			cheby_var(rlen, 
						FF_2BODY[ij].S_MINIM,
						FF_2BODY[ij].S_MAXIM,
						FF_2BODY[ij].LAMBDA,
						FF_2BODY[ij].CHEBY_TYPE,
						x, xdiff, inverse_order, exprlen);
						 
			#if CHECK_CHEBY_RANGE == 1			 

				// Make sure our newly transformed distance falls between the bound of -1 to 1 allowed to Cheby polynomials
						 
				if ( x < -1.0 ||  x > 1.0)
					x = fix_cheby_val(x, inverse_order );
			
				// Now change the range, if the user requested
			
				x = x*CHEBY_DERIV_CONST + FF_2BODY[0].CHEBY_RANGE_LOW - -1.0*CHEBY_DERIV_CONST;

				// Sanity check
	
				if ( x < FF_2BODY[0].CHEBY_RANGE_LOW || x > FF_2BODY[0].CHEBY_RANGE_HIGH )
				{	
					cout << "ERROR: transformed x falls outside user-defined range." << endl;
					cout << "x: " << x << endl;
					cout << "high/low: " << FF_2BODY[0].CHEBY_RANGE_HIGH << " " << FF_2BODY[0].CHEBY_RANGE_LOW  << endl;
					exit_run(0);
				}
			
			#endif

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
			
			// Use recursion to set up the higher n-value Tn and Tnd's
		  
			for ( int i = 2; i <= FF_2BODY[ij].SNUM; i++ ) 
				Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2];
	
			// Now compute the force/potential... Coulomb 
			// Ewald sum not required for a scan-type ("cluster") calculation
			// ...If charges are zero, nothing will be added to tempx
			
			if(INCLUDE_CHARGES)
				tempx += FF_2BODY[ij].ATM1CHG * FF_2BODY[ij].ATM2CHG / rlen;

			// Now compute the force/potential... Cheby

			fcut0 = (1.0 - rlen/FF_2BODY[ij].S_MAXIM);
			fcut      = pow(fcut0, fcut_power);
									 
			for ( int i = 0; i < FF_2BODY[ij].SNUM; i++ ) 
			{
				coeff = FF_2BODY[ij].PARAMS[i]; // This is the Cheby FF paramfor the given power

				if(INCLUDE_FCUT)
					tempx += coeff * fcut * Tn[i+1];
				else
					tempx += coeff * Tn[i+1];
			}
			// Add penalty for very short distances, where the fit FF may be unphysical (preserve conservation of E).

			if ( rpenalty > 0.0 ) 
			{
				Vpenalty = 0.0;
				Vpenalty = rpenalty * rpenalty * rpenalty * penalty_scale;
				tempx += Vpenalty;
			}
			
//			tempx *= SWITCH;
			
			OUTFILE_2B_POT << rlen << " " << tempx << endl;						
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

	double fcut0_ij, fcut0_ik, fcut0_jk; 
	double fcut_ij,  fcut_ik,  fcut_jk; 	
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
	
	//const double fpenalty_scale = FF_2BODY[0].CUBIC_SCALE;		// 1.0;
	

	
	static double CHEBY_DERIV_CONST;	// Accounts for when cheby range is changed from -1:1 to x:y 
	
	
	
	
	const double fcut_power = 
    #ifndef FPENALTY_POWER
		3.0;
	#else
		FPENALTY_POWER;
    #endif
	
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
		
		CHEBY_DERIV_CONST = FF_2BODY[0].CHEBY_RANGE_HIGH - FF_2BODY[0].CHEBY_RANGE_LOW;	// Ranges should be the same for all types
		CHEBY_DERIV_CONST /= 2.0; // i.e the width of the default cheby range
	
	}


	double S_MAXIM_IJ = SET_SMAXIM(FF_2BODY[ij], FF_3BODY[curr_triple_type_index],FF_2BODY[ij].PRPR_NM);
	double S_MAXIM_IK = SET_SMAXIM(FF_2BODY[ik], FF_3BODY[curr_triple_type_index],FF_2BODY[ik].PRPR_NM);
	double S_MAXIM_JK = SET_SMAXIM(FF_2BODY[jk], FF_3BODY[curr_triple_type_index],FF_2BODY[jk].PRPR_NM);
	
	double S_MINIM_IJ = SET_SMINIM(FF_2BODY[ij], FF_3BODY[curr_triple_type_index],FF_2BODY[ij].PRPR_NM);
	double S_MINIM_IK = SET_SMINIM(FF_2BODY[ik], FF_3BODY[curr_triple_type_index],FF_2BODY[ik].PRPR_NM);
	double S_MINIM_JK = SET_SMINIM(FF_2BODY[jk], FF_3BODY[curr_triple_type_index],FF_2BODY[jk].PRPR_NM);

	tempx = 0;

	int n_ij = (S_MAXIM_IJ - S_MINIM_IJ)/FF_2BODY[ij].S_DELTA;
	int n_ik = (S_MAXIM_IK - S_MINIM_IK)/FF_2BODY[ik].S_DELTA;
	int n_jk = (S_MAXIM_JK - S_MINIM_JK)/FF_2BODY[jk].S_DELTA;


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

				RLEN_IJ = S_MINIM_IJ + a * FF_2BODY[ij].S_DELTA;
				RLEN_IK = S_MINIM_IK + b * FF_2BODY[ik].S_DELTA;
				RLEN_JK = S_MINIM_JK + c * FF_2BODY[jk].S_DELTA;
			
				SET_3B_CHEBY_POLYS(FF_2BODY[ij], Tn_ij, Tnd_ij, RLEN_IJ, xdiff_ij, S_MAXIM_IJ, S_MINIM_IJ);
			
				SET_3B_CHEBY_POLYS(FF_2BODY[ik], Tn_ik, Tnd_ik, RLEN_IK, xdiff_ik, S_MAXIM_IK, S_MINIM_IK);
			
				SET_3B_CHEBY_POLYS(FF_2BODY[jk], Tn_jk, Tnd_jk, RLEN_JK, xdiff_jk, S_MAXIM_JK, S_MINIM_JK);		
			
				fcut0_ij     = (1.0 - RLEN_IJ/S_MAXIM_IJ);
				fcut_ij      =  pow(fcut0_ij, fcut_power);
			
				fcut0_ik     = (1.0 - RLEN_IK/S_MAXIM_IK);
				fcut_ik      = pow(fcut0_ik, fcut_power);
			
				fcut0_jk     = (1.0 - RLEN_JK/S_MAXIM_JK);
				fcut_jk      = pow(fcut0_jk, fcut_power);

				// Error check: Certain triplets are impossible...
			
				if(curr_triple_type_index < 0)
				{
					cout << "ERROR: Impossible atom triplet requested: " << TEMP_STR << endl;
					cout << "       Check the parameter file." << endl;
					exit_run(0);
				}
			
				// Now compute the forces for each set of allowed powers for pairs ij, ik, and jk		
				// Keep in mind that the order in which allowed powers are stored may not match the
				// ordering of pairs resulting from the present atom triplet. Thus, we need to order
				// the stored powers properly before applying the FF.
				 
				tempx = 0;
				
				// Coulombic contributions... IF THEY ARE COMPUTED HERE, THEY SHOULDN'T BE COMPUTED IN THE 2-BODY CALCULATIONS!
				// ...Scan logic in the splines_md.C main function control this.
				
				tempx += FF_2BODY[ij].ATM1CHG * FF_2BODY[ij].ATM2CHG / RLEN_IJ;
				tempx += FF_2BODY[ik].ATM1CHG * FF_2BODY[ik].ATM2CHG / RLEN_IK;
				tempx += FF_2BODY[jk].ATM1CHG * FF_2BODY[jk].ATM2CHG / RLEN_JK;
			
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
			
	static double *Tn_ij,  *Tn_ik,  *Tn_jk;
	static double *Tnd_ij, *Tnd_ik, *Tnd_jk;
	static bool   called_before = false;
	double        xdiff_ij, xdiff_ik, xdiff_jk;
	double        tempx;

	double fcut0_ij, fcut0_ik, fcut0_jk; 
	double fcut_ij,  fcut_ik,  fcut_jk; 	
	static string TEMP_STR;
	
	double inner_fcut0_ij, inner_fcut0_ik, inner_fcut0_jk;
	double inner_fcut_ij,  inner_fcut_ik,  inner_fcut_jk;

	static int curr_triple_type_index;
	static int pow_ij, pow_ik, pow_jk;

	static string PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK;

	double RLEN_IJ;
	double RLEN_IK;
	double RLEN_JK;
	
	double RPEN_IJ;
	double RPEN_IK;
	double RPEN_JK;

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
	
	const double fpenalty_scale = FF_2BODY[0].CUBIC_SCALE;		// 1.0;
	// A penalty function is added to the potential for r + penalty_dist < smin[ipair]
	// All pairs have the same penalty scale and distance
	const double penalty_scale  = FF_2BODY[0].PENALTY_SCALE;	// 1.0e8;
	const double penalty_dist   = FF_2BODY[0].PENALTY_DIST;  	// 0.01;
	
	
	static double CHEBY_DERIV_CONST;	// Accounts for when cheby range is changed from -1:1 to x:y
	
	const double fcut_power = 
    #ifndef FPENALTY_POWER
		3.0;
	#else
		FPENALTY_POWER;
    #endif
	
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
		
		CHEBY_DERIV_CONST = FF_2BODY[0].CHEBY_RANGE_HIGH - FF_2BODY[0].CHEBY_RANGE_LOW;	// Ranges should be the same for all types
		CHEBY_DERIV_CONST /= 2.0; // i.e the width of the default cheby range
	
	}

	double S_MAXIM_IJ = SET_SMAXIM(FF_2BODY[ij], FF_3BODY[curr_triple_type_index],FF_2BODY[ij].PRPR_NM);
	double S_MAXIM_IK = SET_SMAXIM(FF_2BODY[ik], FF_3BODY[curr_triple_type_index],FF_2BODY[ik].PRPR_NM);
	double S_MAXIM_JK = SET_SMAXIM(FF_2BODY[jk], FF_3BODY[curr_triple_type_index],FF_2BODY[jk].PRPR_NM);
	
	double S_MINIM_IJ_3B = SET_SMINIM(FF_2BODY[ij], FF_3BODY[curr_triple_type_index],FF_2BODY[ij].PRPR_NM);
	double S_MINIM_IK_3B = SET_SMINIM(FF_2BODY[ik], FF_3BODY[curr_triple_type_index],FF_2BODY[ik].PRPR_NM);
	double S_MINIM_JK_3B = SET_SMINIM(FF_2BODY[jk], FF_3BODY[curr_triple_type_index],FF_2BODY[jk].PRPR_NM);
	
	double S_MINIM_IJ_2B = FF_2BODY[ij].S_MINIM;
	double S_MINIM_IK_2B = FF_2BODY[ik].S_MINIM;
	double S_MINIM_JK_2B = FF_2BODY[jk].S_MINIM;
	
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
		n = (S_MAXIM_IJ - S_MINIM_IJ_2B)/FF_2BODY[ij].S_DELTA;
	
	else if (FF_PLOTS.SCAN_PAIR[scan] == 2)
		n = (S_MAXIM_IK - S_MINIM_IK_2B)/FF_2BODY[ik].S_DELTA;
	
	else
		n = (S_MAXIM_JK - S_MINIM_JK_2B)/FF_2BODY[jk].S_DELTA;
	
	// Run the scan	

	for (int a=1; a<n; a++)
	{

		if      (FF_PLOTS.SCAN_PAIR[scan] == 1)
			RLEN_IJ = S_MINIM_IJ_2B + a * FF_2BODY[ij].S_DELTA;
		else if (FF_PLOTS.SCAN_PAIR[scan] == 2)
			RLEN_IK = S_MINIM_IK_2B + a * FF_2BODY[ik].S_DELTA;
		else
			RLEN_JK = S_MINIM_JK_2B + a * FF_2BODY[jk].S_DELTA;		
		
		// Set up for switching function... 
		
		tempx = 0;
		
		double SWITCH_IJ,  SWITCH_IK,  SWITCH_JK;
		double DSWITCH_IJ, DSWITCH_IK, DSWITCH_JK;
		double SWITCH_ENDPOINT;
		double SWITCH_SMOOTHNESS = 2.0; // = 0.1;
		
		SWITCH_IJ = 0.0;
		SWITCH_IK = 0.0;
		SWITCH_JK = 0.0;
		
		DSWITCH_IJ = 0.0;
		DSWITCH_IK = 0.0;
		DSWITCH_JK = 0.0;
		
		
		RPEN_IJ = 0.0;
		RPEN_IK = 0.0;
		RPEN_JK = 0.0;
		
		if (  (FF_PLOTS.SCAN_PAIR[scan] == 1) && (RLEN_IJ < S_MINIM_IJ_3B )) 
		{
			RPEN_IJ   = S_MINIM_IJ_3B + penalty_dist - RLEN_IJ;
			//SWITCH_IJ = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_IJ_3B - RLEN_IJ)/S_MINIM_IJ_3B);
			SWITCH_IJ = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_IJ_3B - RLEN_IJ));
			
			if(RLEN_IJ < S_MINIM_IJ_3B)
				RLEN_IJ =  S_MINIM_IJ_3B;
		}
		if(RLEN_IJ < FF_2BODY[ij].S_MINIM )
		{
			cout << "ERROR: IJ pair is fixed below 2b smin value!" << endl;
			cout << "rlen, 2b_smin: " << RLEN_IJ  << " " << FF_2BODY[ij].S_MINIM << endl;
			exit(0);
		}
		
		if (  (FF_PLOTS.SCAN_PAIR[scan] == 2) && (RLEN_IK < S_MINIM_IK_3B )) 
		{
			RPEN_IK   = S_MINIM_IK_3B + penalty_dist - RLEN_IK;
			//SWITCH_IK = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_IK_3B - RLEN_IK)/S_MINIM_IK_3B);
			SWITCH_IK = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_IK_3B - RLEN_IK));
			
			if(RLEN_IK < S_MINIM_IK_3B)
				RLEN_IK =  S_MINIM_IK_3B;
		}
		if(RLEN_IK < FF_2BODY[ik].S_MINIM )
		{
			cout << "ERROR: IK pair is fixed below 2b smin value!" << endl;
			cout << "rlen, 2b_smin: " << RLEN_IK  << " " << FF_2BODY[ik].S_MINIM << endl;
			exit(0);
		}
		
		if (  (FF_PLOTS.SCAN_PAIR[scan] == 3) && (RLEN_JK < S_MINIM_JK_3B )) 
		{
			RPEN_JK   = S_MINIM_JK_3B + penalty_dist - RLEN_JK;
			//SWITCH_JK = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_JK_3B - RLEN_JK)/S_MINIM_JK_3B);
			SWITCH_JK = penalty_scale * tanh( SWITCH_SMOOTHNESS*(S_MINIM_JK_3B - RLEN_JK));
			
			if(RLEN_JK < S_MINIM_JK_3B)
				RLEN_JK =  S_MINIM_JK_3B;
		}
		if(RLEN_JK < FF_2BODY[jk].S_MINIM )
		{
			cout << "ERROR: JK pair is fixed below 2b smin value!" << endl;
			cout << "rlen, 2b_smin: " << RLEN_JK  << " " << FF_2BODY[jk].S_MINIM << endl;
			exit(0);
		}		

		SET_3B_CHEBY_POLYS(FF_2BODY[ij], Tn_ij, Tnd_ij, RLEN_IJ, xdiff_ij, S_MAXIM_IJ, S_MINIM_IJ_3B);
		SET_3B_CHEBY_POLYS(FF_2BODY[ik], Tn_ik, Tnd_ik, RLEN_IK, xdiff_ik, S_MAXIM_IK, S_MINIM_IK_3B);
		SET_3B_CHEBY_POLYS(FF_2BODY[jk], Tn_jk, Tnd_jk, RLEN_JK, xdiff_jk, S_MAXIM_JK, S_MINIM_JK_3B);		
	
		fcut0_ij     = (1.0 - RLEN_IJ/S_MAXIM_IJ);
		fcut_ij      = pow(fcut0_ij, fcut_power);
	
		fcut0_ik     = (1.0 - RLEN_IK/S_MAXIM_IK);
		fcut_ik      = pow(fcut0_ik, fcut_power);
	
		fcut0_jk     = (1.0 - RLEN_JK/S_MAXIM_JK);
		fcut_jk      = pow(fcut0_jk, fcut_power);
	
		fcut_ij *= fpenalty_scale;
		fcut_ik *= fpenalty_scale;
		fcut_jk *= fpenalty_scale;

		// Error check: Certain triplets are impossible...
	
		if(curr_triple_type_index < 0)
		{
			cout << "ERROR: Impossible atom triplet found: " << TEMP_STR << endl;
			cout << "       Check the parameter file." << endl;
			exit_run(0);
		}
	
		// Now compute the forces for each set of allowed powers for pairs ij, ik, and jk		
		// Keep in mind that the order in which allowed powers are stored may not match the
		// ordering of pairs resulting from the present atom triplet. Thus, we need to order
		// the stored powers properly before applying the FF.

		// Coulombic contributions...IF THEY ARE COMPUTED HERE, THEY SHOULDN'T BE COMPUTED IN THE 2-BODY CALCULATIONS!
		// ...Scan logic in the splines_md.C main function control this.
		
		tempx += FF_2BODY[ij].ATM1CHG * FF_2BODY[ij].ATM2CHG / RLEN_IJ;
		tempx += FF_2BODY[ik].ATM1CHG * FF_2BODY[ik].ATM2CHG / RLEN_IK;
		tempx += FF_2BODY[jk].ATM1CHG * FF_2BODY[jk].ATM2CHG / RLEN_JK;			

		for(int i=0; i<FF_3BODY[curr_triple_type_index].N_ALLOWED_POWERS; i++) 
		{
			SET_3B_CHEBY_POWERS(FF_2BODY, FF_3BODY[curr_triple_type_index], PAIR_MAP,  pow_ij, pow_ik, pow_jk, PAIR_TYPE_IJ, PAIR_TYPE_IK, PAIR_TYPE_JK, i);
		
			tempx += FF_3BODY[curr_triple_type_index].PARAMS[i] 
				   * fcut_ij   * fcut_ik   * fcut_jk 
				   * Tn_ij[pow_ij] * Tn_ik[pow_ik] * Tn_jk[pow_jk];		
		}

		tempx += SWITCH_IJ + SWITCH_IK + SWITCH_JK;
		
		// Make sure we're using the right rlen for printing
		
		if      (FF_PLOTS.SCAN_PAIR[scan] == 1)
			RLEN_IJ = S_MINIM_IJ_2B + a * FF_2BODY[ij].S_DELTA;
		else if (FF_PLOTS.SCAN_PAIR[scan] == 2)
			RLEN_IK = S_MINIM_IK_2B + a * FF_2BODY[ik].S_DELTA;
		else
			RLEN_JK = S_MINIM_JK_2B + a * FF_2BODY[jk].S_DELTA;	
				
		OUTFILE_3B_POT << 	RLEN_IJ << " " << RLEN_IK << " " << RLEN_JK << " " << tempx <<endl;						
			
	}
	
	OUTFILE_3B_POT.close();
}
