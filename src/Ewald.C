#include<iomanip>
#include "functions.h"

#define EWALD_ACCURACY 1.0e-07

using namespace std;

//////////////////////////////////////////
//
//	FUNCTION HEADERS
//
//////////////////////////////////////////



static double add_sines    (int kx, int ky, int kz, vector<XYZ> & SIN_XYZ, vector<XYZ> & COS_XYZ);
static void   generate_trig(vector<XYZ> & SIN_XYZ, vector<XYZ> & COS_XYZ, XYZ & RVEC, XYZ & BOXDIM, int kmax);
static void   Ewald_K_Space_New(double alphasq, int k_cut, FRAME & TRAJECTORY, double & UCoul);	
// Function Overloaded .. Different input required depending on whether used by splines_ls or splines_md -- THIS PROBABLY DOESN'T NEED TO BE OVERLOADED
static void   Ewald_K_Space_Orig(double alphasq, FRAME & TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP, double & Vtot);	// LSQ
static void   Ewald_K_Space_Orig(double alphasq, FRAME & TRAJECTORY, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, double & Vtot);	// MD

//////////////////////////////////////////
//
//	FUNCTION DEFINITIONS
//
//////////////////////////////////////////

/* ... I DONT THINK THIS NEEDS OVERLOADING... TEST ON LSQ CODE WITH MD FUNCTION!

static void Ewald_K_Space_New(double alphasq, int k_cut, FRAME & TRAJECTORY)				// LSQ version
// Calculate Ewald K-space components.  Use a rearrangement of the usual Ewald
// expression to generate an order-N evaluation.   See A. Y. Toukmaji et. al,
// Comp. Phys. Comm. 95, 73-92 (1996).
{
	double Volume   = TRAJECTORY.BOXDIM.X * TRAJECTORY.BOXDIM.Y * TRAJECTORY.BOXDIM.Z;
	const double PI = M_PI;
  
	//set up Ewald Coulomb parameters:
  
	double        tempk,tempy;
	static double alpha = 0.0;
	const int     maxk =10000;
	int           kmax =10;
	int           ksq;
	double        Kfac;
	XYZ			  R_K;//        rkx,rky,rkz;
	double        rksq;
	static XYZ    LAST_BOXDIMS;

	//set up Kfac storage vec:
 
	static bool     called_before = false;
	static int      totk;
	static vector <double> KFAC_V(maxk);
	static vector <XYZ>    K_V(maxk);
	vector <double> sin_array(TRAJECTORY.ATOMS);
	vector <double> cos_array(TRAJECTORY.ATOMS);


	if (!called_before) 
	{
		called_before = true;
		LAST_BOXDIMS.X = LAST_BOXDIMS.Y = LAST_BOXDIMS.Z = 0.0;

	}

	if ( TRAJECTORY.BOXDIM.X != LAST_BOXDIMS.X || TRAJECTORY.BOXDIM.Y != LAST_BOXDIMS.Y || TRAJECTORY.BOXDIM.Z != LAST_BOXDIMS.Z ) 
	{
		// Update K factors when box dimensions change.
	
		double min_vfac = 1.0e8;
		totk = 0;
		alpha = sqrt(alphasq);
		int ksqmax = k_cut * k_cut;
		kmax = k_cut;

		for(int kx=-kmax; kx<=kmax; kx++) 
		{
			for(int ky=-kmax; ky<=kmax; ky++)
			{
				for(int kz=-kmax; kz<=kmax; kz++)
				{
					ksq=kx*kx + ky*ky + kz*kz;
					
					if(ksq!=0 && ksq<ksqmax)
					{						
						K_V[totk].X = (2.0*PI/TRAJECTORY.BOXDIM.X)*kx;
						K_V[totk].Y = (2.0*PI/TRAJECTORY.BOXDIM.Y)*ky;
						K_V[totk].Z = (2.0*PI/TRAJECTORY.BOXDIM.Z)*kz;
						
						rksq = K_V[totk].X*K_V[totk].X + K_V[totk].Y*K_V[totk].Y + K_V[totk].Z*K_V[totk].Z;
						
						// Note:  The original version of the Ewald evaluator did not work 
						// when the lattice constants were different.
						
						KFAC_V[totk] = exp(-rksq/(4.0*alphasq))/rksq;
						
						if ( KFAC_V[totk] < min_vfac ) 
							min_vfac = KFAC_V[totk];
						
						totk++;
						
						if(totk>maxk)
						{
							cout << "	totk = " << totk << " greater than maxk = " << maxk << endl;
							exit(1);
						}
					}
				}
			}
		}
		
		if ( LAST_BOXDIMS.X == 0.0 ) 
		{
			#if VERBOSITY == 1
				cout << "	Cutting off K-vectors with magnitude less than "  << k_cut << endl;
				cout << "	Number of Ewald K-vectors = " << totk << endl;
				cout << "	K-space accuracy estimate = " << min_vfac << endl;
			#endif
				
				cout << "	Cutting off K-vectors with magnitude less than "  << k_cut << endl;
				cout << "	Number of Ewald K-vectors = " << totk << endl;
				cout << "	K-space accuracy estimate = " << min_vfac << endl;
		}
		
		LAST_BOXDIMS.X = TRAJECTORY.BOXDIM.X;
		LAST_BOXDIMS.Y = TRAJECTORY.BOXDIM.Y;
		LAST_BOXDIMS.Z = TRAJECTORY.BOXDIM.Z;
	}

	double UCoul = 0;

	// Main loop Ewald Coulomb energy/forces: K-Space loop.
	
	tempk=0;
	
	double tempsin = 0.0;
	double tempcos = 0.0;
	double kdotr;
	double tempd;

	for(int ik=-0; ik<totk; ik++)
	{
		tempsin = 0.0;
		tempcos = 0.0;

		for(int a1=0; a1<TRAJECTORY.ATOMS; a1++) 
		{
			kdotr = TRAJECTORY.COORDS[a1].X * K_V[ik].X + TRAJECTORY.COORDS[a1].Y * K_V[ik].Y + TRAJECTORY.COORDS[a1].Z * K_V[ik].Z;	
			sin_array[a1] = sin(kdotr);																									
			cos_array[a1] = cos(kdotr);						
			tempcos += TRAJECTORY.CHARGES[a1] * cos_array[a1];																			
			tempsin += TRAJECTORY.CHARGES[a1] * sin_array[a1];																			

		}
		
		tempk += KFAC_V[ik]* (tempcos*tempcos + tempsin*tempsin);
		
		for(int a1=0; a1<TRAJECTORY.ATOMS; a1++) 
		{
			tempd = -4.0 * PI * KFAC_V[ik] * ke * TRAJECTORY.CHARGES[a1] * ( tempcos * sin_array[a1] - tempsin * cos_array[a1] ) / Volume;
		
			TRAJECTORY.TMP_EWALD[a1].X += K_V[ik].X * tempd;
			TRAJECTORY.TMP_EWALD[a1].Y += K_V[ik].Y * tempd;
			TRAJECTORY.TMP_EWALD[a1].Z += K_V[ik].Z * tempd;
		}
	}
	
	tempy=(2*PI/Volume)*tempk;
	UCoul+=ke*tempy;

	// Constant energy term.
	for ( int a1=0; a1<TRAJECTORY.ATOMS; a1++ ) 
		UCoul -= ke * alpha /sqrt(PI) * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[a1];


	// End K-Space loop.
	
	// printf("K-space part of Ewald Energy = %13.6e\n", UCoul);

	TRAJECTORY.TOT_POT_ENER += UCoul;
	
	vector<double>().swap(sin_array);
	vector<double>().swap(cos_array);
	
}
*/

static void Ewald_K_Space_New(double alphasq, int k_cut, FRAME & TRAJECTORY, double & UCoul)				// MD version
// Calculate Ewald K-space components.  Use a rearrangement of the usual Ewald
// expression to generate an order-N evaluation.   See A. Y. Toukmaji et. al,
// Comp. Phys. Comm. 95, 73-92 (1996).
{
	double Volume   = TRAJECTORY.BOXDIM.X * TRAJECTORY.BOXDIM.Y * TRAJECTORY.BOXDIM.Z;
	const double PI = M_PI;
  
	//set up Ewald Coulomb parameters:
  
	double        tempk,tempy;
	static double alpha = 0.0;
	const int     maxk =10000;
	int           kmax =10;
	int           ksq;
	double        Kfac;
	XYZ			  R_K;//        rkx,rky,rkz;
	double        rksq;
	static XYZ    LAST_BOXDIMS;

	//set up Kfac storage vec:
 
	static bool     called_before = false;
	static int      totk;
	static vector <double> KFAC_V(maxk);
	static vector <XYZ>    K_V(maxk);
	vector <double> sin_array(TRAJECTORY.ATOMS);
	vector <double> cos_array(TRAJECTORY.ATOMS);


	if (!called_before) 
	{
		called_before = true;
		LAST_BOXDIMS.X = LAST_BOXDIMS.Y = LAST_BOXDIMS.Z = 0.0;
		
		string TEMP_STR;
		
	}

	if ( TRAJECTORY.BOXDIM.X != LAST_BOXDIMS.X || TRAJECTORY.BOXDIM.Y != LAST_BOXDIMS.Y || TRAJECTORY.BOXDIM.Z != LAST_BOXDIMS.Z ) 
	{
		// Update K factors when box dimensions change.
	
		double min_vfac = 1.0e8;
		totk = 0;
		alpha = sqrt(alphasq);
		int ksqmax = k_cut * k_cut;
		kmax = k_cut;

		for(int kx=-kmax; kx<=kmax; kx++) 
		{
			for(int ky=-kmax; ky<=kmax; ky++)
			{
				for(int kz=-kmax; kz<=kmax; kz++)
				{
					ksq=kx*kx + ky*ky + kz*kz;
					
					if(ksq!=0 && ksq<ksqmax)
					{						
						K_V[totk].X = (2.0*PI/TRAJECTORY.BOXDIM.X)*kx;
						K_V[totk].Y = (2.0*PI/TRAJECTORY.BOXDIM.Y)*ky;
						K_V[totk].Z = (2.0*PI/TRAJECTORY.BOXDIM.Z)*kz;
						
						rksq = K_V[totk].X*K_V[totk].X + K_V[totk].Y*K_V[totk].Y + K_V[totk].Z*K_V[totk].Z;
						
						// Note:  The original version of the Ewald evaluator did not work 
						// when the lattice constants were different.
						
						KFAC_V[totk] = exp(-rksq/(4.0*alphasq))/rksq;
						
						if ( KFAC_V[totk] < min_vfac ) 
							min_vfac = KFAC_V[totk];
						
						totk++;
						
						if(totk>maxk)
						{
							cout << "	totk = " << totk << " greater than maxk = " << maxk << endl;
							exit(1);
						}
					}
				}
			}
		}
		
		if ( LAST_BOXDIMS.X == 0.0 ) 
		{
			#if VERBOSITY == 1
				cout << "	Cutting off K-vectors with magnitude less than "  << k_cut << endl;
				cout << "	Number of Ewald K-vectors = " << totk << endl;
				cout << "	K-space accuracy estimate = " << min_vfac << endl << endl;;
			#endif

		}
		
		LAST_BOXDIMS.X = TRAJECTORY.BOXDIM.X;
		LAST_BOXDIMS.Y = TRAJECTORY.BOXDIM.Y;
		LAST_BOXDIMS.Z = TRAJECTORY.BOXDIM.Z;
	}

	UCoul = 0;

	// Main loop Ewald Coulomb energy/forces: K-Space loop.
	
	tempk=0;
	
	double tempsin = 0.0;
	double tempcos = 0.0;
	double kdotr;
	double tempd;

	for(int ik=-0; ik<totk; ik++)
	{
		tempsin = 0.0;
		tempcos = 0.0;

		for(int a1=0; a1<TRAJECTORY.ATOMS; a1++) 
		{
			kdotr = TRAJECTORY.COORDS[a1].X * K_V[ik].X + TRAJECTORY.COORDS[a1].Y * K_V[ik].Y + TRAJECTORY.COORDS[a1].Z * K_V[ik].Z;	
			sin_array[a1] = sin(kdotr);																									
			cos_array[a1] = cos(kdotr);		
			
			tempcos += TRAJECTORY.CHARGES[a1] * cos_array[a1];																			
			tempsin += TRAJECTORY.CHARGES[a1] * sin_array[a1];		
		}
		
		tempk += KFAC_V[ik]* (tempcos*tempcos + tempsin*tempsin);
		
		for(int a1=0; a1<TRAJECTORY.ATOMS; a1++) 
		{
			tempd = -4.0 * PI * KFAC_V[ik] * ke * TRAJECTORY.CHARGES[a1] * ( tempcos * sin_array[a1] - tempsin * cos_array[a1] ) / Volume;
		
			TRAJECTORY.TMP_EWALD[a1].X += K_V[ik].X * tempd;
			TRAJECTORY.TMP_EWALD[a1].Y += K_V[ik].Y * tempd;
			TRAJECTORY.TMP_EWALD[a1].Z += K_V[ik].Z * tempd;
		}
	}
	
	tempy=(2*PI/Volume)*tempk;
	UCoul+=ke*tempy;

	// Constant energy term.
	for ( int a1=0; a1<TRAJECTORY.ATOMS; a1++ ) 
		UCoul -= ke * alpha /sqrt(PI) * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[a1];

	// End K-Space loop.
	
	// printf("K-space part of Ewald Energy = %13.6e\n", UCoul);

	TRAJECTORY.TOT_POT_ENER += UCoul;
	
	vector<double>().swap(sin_array);
	vector<double>().swap(cos_array);
	
}

// FUNCTION UPDATED -- OVERLOADING THE FUNCTION.. DIFFERENT INPUT REQUIRED DEPENDING ON WHETHER FUNCTION IS CALLED FROM MD
// ... I don't think Ewald_K_Space_Orig actually ever get used... these have NOT been updated, and likely have logic errors.
static void Ewald_K_Space_Orig(double alphasq, FRAME & TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, map<string,int> & PAIR_MAP, double & Vtot)	// LSQ version
{
	XYZ RVEC;
	double rlen;
	double tempx;
	double Volume = TRAJECTORY.BOXDIM.X *  TRAJECTORY.BOXDIM.Y *  TRAJECTORY.BOXDIM.Z;
	const double PI=3.14159265359;

	cout.precision(15);

	double			tempk,tempd,tempy;
	const double	alpha  = sqrt(alphasq);
	const int		maxk   = 1000000;
	const int 		kmax   = 10;
	const int 		ksqmax = 50;
	int 			totk;
	int 			ksq;
	double 			Kfac;
	XYZ				RK;
	double 			rksq;
   

	double UCoul = 0;
	totk=0;

	string TEMP_STR;
	int curr_pair_type_idx; 	

	// Main loop Ewald Coulomb energy/forces:
	
	for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)
	{
		for(int a2=0;a2<TRAJECTORY.ATOMS;a2++)
		{
			tempy=0.0;
			
			TEMP_STR = TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(TRAJECTORY.ATOMTYPE[a2]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];		
			
			RVEC.X = TRAJECTORY.COORDS[a2].X - TRAJECTORY.COORDS[a1].X; 
			RVEC.Y = TRAJECTORY.COORDS[a2].Y - TRAJECTORY.COORDS[a1].Y; 
			RVEC.Z = TRAJECTORY.COORDS[a2].Z - TRAJECTORY.COORDS[a1].Z; 

	
			rlen = sqrt(RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z);



			if(a1==a2)
				tempy-=(2*alpha)/sqrt(PI);

			tempk=0;
			totk=0;

			for(int kx=0; kx<=kmax; kx++)
			{
				RK.X = (2.0*PI/TRAJECTORY.BOXDIM.X)*kx;
				
				for(int ky=-kmax; ky<=kmax; ky++)
				{
					RK.Y = (2.0*PI/TRAJECTORY.BOXDIM.Y)*ky;
					
					for(int kz=-kmax; kz<=kmax; kz++)
					{
						RK.Z = (2.0*PI/TRAJECTORY.BOXDIM.Z)*kz;
						
						ksq = kx*kx + ky*ky + kz*kz;
		    
						if(ksq!=0 && ksq<ksqmax)
						{
							totk++;
							
							if(totk > maxk)
							{
								cout << "totk = " << totk << " greater than maxk = " << maxk << endl;
								exit(1);
							}
							rksq = RK.X*RK.X + RK.Y*RK.Y + RK.Z*RK.Z;
							Kfac = exp(-1.0*rksq/(4.0*alpha*alpha))/rksq;
		
							//this is because we want to get a 1/2 speedup by using only 0<kmax for x-coordinate.
							if(kx == 0)
								tempk +=     Kfac * cos( RK.X*RVEC.X + RK.Y*RVEC.Y + RK.Z*RVEC.Z );
							else
								tempk += 2 * Kfac * cos( RK.X*RVEC.X + RK.Y*RVEC.Y + RK.Z*RVEC.Z );
					 
							tempd = sin( RK.X*RVEC.X + RK.Y*RVEC.Y + RK.Z*RVEC.Z );
								
							
							tempx  = Kfac * tempd * ATOM_PAIRS[curr_pair_type_idx].ATM1CHG * ATOM_PAIRS[curr_pair_type_idx].ATM2CHG * 0.5 * (4*PI/Volume);
							tempx *= ke;
							
							// Signs are now flipped, because in the previous version, the TRAJECTORY lines below would be like
							//
							// FCoul[a1] += ...
							// FCoul[a2] -= ...
							//
							// And then at the end of the loop
							//
							// SForce[a1] -= Fcoul
							//
							// Where SForce is then 

							TRAJECTORY.FORCES[a1].X += RK.X*tempx;
							TRAJECTORY.FORCES[a1].Y += RK.Y*tempx;
							TRAJECTORY.FORCES[a1].Z += RK.Z*tempx;
		
							TRAJECTORY.FORCES[a2].X -= RK.X*tempx;
							TRAJECTORY.FORCES[a2].Y -= RK.Y*tempx;
							TRAJECTORY.FORCES[a2].Z -= RK.Z*tempx;

							if(kx != 0)
							{
								TRAJECTORY.FORCES[a1].X += RK.X*tempx;
								TRAJECTORY.FORCES[a1].Y += RK.Y*tempx;
								TRAJECTORY.FORCES[a1].Z += RK.Z*tempx;
		
								TRAJECTORY.FORCES[a2].X -= RK.X*tempx;
								TRAJECTORY.FORCES[a2].Y -= RK.Y*tempx;
								TRAJECTORY.FORCES[a2].Z -= RK.Z*tempx;
							}
		    			}
		 			}
	     		}
			}

			tempy += (4*PI/Volume)*tempk;
			UCoul += ke * 0.5 * ATOM_PAIRS[curr_pair_type_idx].ATM1CHG * ATOM_PAIRS[curr_pair_type_idx].ATM2CHG * tempy;
		} //end main loop Ewald Coulomb forces:
	}
  Vtot+=UCoul;

  return;
}
static void Ewald_K_Space_Orig(double alphasq, FRAME & TRAJECTORY, vector<PAIR_FF> & FF_2BODY, map<string,int> & PAIR_MAP, double & Vtot)	// MD version
{
	XYZ RVEC;
	double rlen;
	double tempx;
	double Volume = TRAJECTORY.BOXDIM.X *  TRAJECTORY.BOXDIM.Y *  TRAJECTORY.BOXDIM.Z;
	const double PI=3.14159265359;

	cout.precision(15);

	double			tempk,tempd,tempy;
	const double	alpha  = sqrt(alphasq);
	const int		maxk   = 1000000;
	const int 		kmax   = 10;
	const int 		ksqmax = 50;
	int 			totk;
	int 			ksq;
	double 			Kfac;
	XYZ				RK;
	double 			rksq;
   

	double UCoul = 0;
	totk=0;

	/*
	string TEMP_STR;
	int curr_pair_type_idx; 
	*/
	

	// Main loop Ewald Coulomb energy/forces:
	
	for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)
	{
		for(int a2=0;a2<TRAJECTORY.ATOMS;a2++)
		{
			tempy=0.0;
/*			
			TEMP_STR = TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(TRAJECTORY.ATOMTYPE[a2]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];			
*/			
			RVEC.X = TRAJECTORY.COORDS[a2].X - TRAJECTORY.COORDS[a1].X; 
			RVEC.Y = TRAJECTORY.COORDS[a2].Y - TRAJECTORY.COORDS[a1].Y; 
			RVEC.Z = TRAJECTORY.COORDS[a2].Z - TRAJECTORY.COORDS[a1].Z; 

	
			rlen = sqrt(RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z);



			if(a1==a2)
				tempy-=(2*alpha)/sqrt(PI);

			tempk=0;
			totk=0;

			for(int kx=0; kx<=kmax; kx++)
			{
				RK.X = (2.0*PI/TRAJECTORY.BOXDIM.X)*kx;
				
				for(int ky=-kmax; ky<=kmax; ky++)
				{
					RK.Y = (2.0*PI/TRAJECTORY.BOXDIM.Y)*ky;
					
					for(int kz=-kmax; kz<=kmax; kz++)
					{
						RK.Z = (2.0*PI/TRAJECTORY.BOXDIM.Z)*kz;
						
						ksq = kx*kx + ky*ky + kz*kz;
		    
						if(ksq!=0 && ksq<ksqmax)
						{
							totk++;
							
							if(totk > maxk)
							{
								cout << "totk = " << totk << " greater than maxk = " << maxk << endl;
								exit(1);
							}
							rksq = RK.X*RK.X + RK.Y*RK.Y + RK.Z*RK.Z;
							Kfac = exp(-1.0*rksq/(4.0*alpha*alpha))/rksq;
							
							//this is because we want to get a 1/2 speedup by using only 0<kmax for x-coordinate.
							if(kx == 0)
								tempk +=     Kfac * cos( RK.X*RVEC.X + RK.Y*RVEC.Y + RK.Z*RVEC.Z );
							else
								tempk += 2 * Kfac * cos( RK.X*RVEC.X + RK.Y*RVEC.Y + RK.Z*RVEC.Z );
					 
							tempd = sin( RK.X*RVEC.X + RK.Y*RVEC.Y + RK.Z*RVEC.Z );
								
							
							tempx  = Kfac * tempd * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[a2] * 0.5 * (4*PI/Volume);
							tempx *= ke;

							TRAJECTORY.FORCES[a1].X += RK.X*tempx;
							TRAJECTORY.FORCES[a1].Y += RK.Y*tempx;
							TRAJECTORY.FORCES[a1].Z += RK.Z*tempx;
		
							TRAJECTORY.FORCES[a2].X -= RK.X*tempx;
							TRAJECTORY.FORCES[a2].Y -= RK.Y*tempx;
							TRAJECTORY.FORCES[a2].Z -= RK.Z*tempx;

							if(kx != 0)
							{
								TRAJECTORY.FORCES[a1].X += RK.X*tempx;
								TRAJECTORY.FORCES[a1].Y += RK.Y*tempx;
								TRAJECTORY.FORCES[a1].Z += RK.Z*tempx;
		
								TRAJECTORY.FORCES[a2].X -= RK.X*tempx;
								TRAJECTORY.FORCES[a2].Y -= RK.Y*tempx;
								TRAJECTORY.FORCES[a2].Z -= RK.Z*tempx;
							}
		    			}
		 			}
	     		}
			}

			tempy += (4*PI/Volume)*tempk;
			UCoul += ke * 0.5 * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[a2] * tempy;
		} //end main loop Ewald Coulomb forces:
	}
  Vtot+=UCoul;

  return;
}


// FUNCTION UPDATED -- OVERLOADING THE FUNCTION.. DIFFERENT INPUT REQUIRED DEPENDING ON WHETHER FUNCTION IS CALLED FROM MD
//                     PROGRAM OR LSQ FITTING PROGRAM...
void ZCalc_Ewald(FRAME & TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, map<string, int>  & PAIR_MAP)	// LSQ version
// Calculate Ewald interactions.
{
  XYZ RVEC; 		// Replaces Rvec[3];
  double tempx;

  const double PI = M_PI;
  const double SQRT_PI = sqrt(PI);

  double 		tempd,tempy;
  double 		dx,dy,dz;
  double 		rlen_mi;
  static double alphasq;
  int 			totk = 0;
  static double r_cut;
  static int	k_cut;
  bool 			if_old_k_space = false;
  static bool 	called_before = false;
  static double alpha;
  const double 	accuracy = EWALD_ACCURACY;
  double		UCoul = 0;
  double 		erfc_val;
  
  string TEMP_STR;


	if ( ! called_before ) 
	{
		double vol = TRAJECTORY.BOXDIM.X * TRAJECTORY.BOXDIM.Y * TRAJECTORY.BOXDIM.Z;
		double r_acc, k_acc;
		(accuracy, vol, TRAJECTORY.ATOMS, alpha, r_cut, k_cut,r_acc, k_acc);	// NEEDS UPDATING!!!
		alphasq = alpha * alpha;
      
		#if VERBOSITY == 1
			printf("Requested Ewald accuracy  = %13.6e\n", accuracy);
			printf("Ewald alpha               = %13.6e\n", alpha);
			printf("R-Space Ewald cutoff      = %13.6e\n", r_cut);
			printf("R-Space accuracy estimate = %13.6e\n", r_acc);
			printf("K-Space accuracy estimate = %13.6e\n", k_acc);
		#endif

		called_before = true;
	}

	// Main loop Ewald Coulomb energy/forces:
	
	int curr_pair_type_idx;	

	for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)		// Double sum over atom pairs
	{
		for(int a2=0;a2<a1;a2++)
		{
			tempy=0.0;
			
			TEMP_STR = TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(TRAJECTORY.ATOMTYPE[a2]);
			curr_pair_type_idx = PAIR_MAP[TEMP_STR];		
	
			RVEC.X = TRAJECTORY.COORDS[a2].X - TRAJECTORY.COORDS[a1].X;
			RVEC.Y = TRAJECTORY.COORDS[a2].Y - TRAJECTORY.COORDS[a1].Y;
			RVEC.Z = TRAJECTORY.COORDS[a2].Z - TRAJECTORY.COORDS[a1].Z;
	
			// THIS IS THE REAL-SPACE LOOP IN MIC:
	
			RVEC.X = RVEC.X - floor( 0.5 + RVEC.X / TRAJECTORY.BOXDIM.X ) * TRAJECTORY.BOXDIM.X;
			RVEC.Y = RVEC.Y - floor( 0.5 + RVEC.Y / TRAJECTORY.BOXDIM.Y ) * TRAJECTORY.BOXDIM.Y;
			RVEC.Z = RVEC.Z - floor( 0.5 + RVEC.Z / TRAJECTORY.BOXDIM.Z ) * TRAJECTORY.BOXDIM.Z;

			rlen_mi = RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z;
	  
			if ( rlen_mi < r_cut * r_cut ) 
			{
				rlen_mi = sqrt(rlen_mi);

				erfc_val = erfc(alpha * rlen_mi);
	    
				tempy += erfc_val/rlen_mi;
	    
				tempd = ( (-2.0/SQRT_PI) * exp(-1.0*alphasq*rlen_mi*rlen_mi) * rlen_mi*alpha - erfc_val ) / (rlen_mi*rlen_mi*rlen_mi);

				tempx =- tempd * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[a1];
				tempx *= ke;
				
		        TRAJECTORY.FORCES[a1].X += RVEC.X * tempx;
				TRAJECTORY.FORCES[a1].Y += RVEC.X * tempx;
				TRAJECTORY.FORCES[a1].Z += RVEC.X * tempx;
				
		        TRAJECTORY.FORCES[a2].X -= RVEC.X * tempx;
				TRAJECTORY.FORCES[a2].Y -= RVEC.Y * tempx;
				TRAJECTORY.FORCES[a2].Z -= RVEC.Z * tempx;								

				UCoul += ke * ATOM_PAIRS[curr_pair_type_idx].ATM1CHG * ATOM_PAIRS[curr_pair_type_idx].ATM2CHG * tempy;
			} 
			else 
			{
				// printf("Cutoff exclusion found\n");
			}
		}
	}

	// K-Space loops.
	
	if ( if_old_k_space ) 
		Ewald_K_Space_Orig(alphasq, TRAJECTORY, ATOM_PAIRS, PAIR_MAP, UCoul);
	else 
		Ewald_K_Space_New(alphasq, k_cut, TRAJECTORY, UCoul);

	return;

}

void ZCalc_Ewald(FRAME & TRAJECTORY)	// MD version
// Calculate Ewald interactions. 
{

  XYZ RVEC; 		// Replaces Rvec[3];
  double tempx;

  const double PI = M_PI;
  const double SQRT_PI = sqrt(PI);

  double 		tempd,tempy;
  double 		dx,dy,dz;
  double 		rlen_mi;
  static double alphasq;
  int 			totk = 0;
  static double r_cut;
  static int	k_cut;
  bool 			if_old_k_space = false;
  static bool 	called_before = false;
  static double alpha;
  const double 	accuracy = EWALD_ACCURACY;
  double		UCoul = 0;
  double		TMP_UCoul;
  double 		erfc_val;
  
  string TEMP_STR;

	if ( ! called_before ) 
	{
		double vol = TRAJECTORY.BOXDIM.X * TRAJECTORY.BOXDIM.Y * TRAJECTORY.BOXDIM.Z;
		double r_acc, k_acc;
		optimal_ewald_params(accuracy, vol, TRAJECTORY.ATOMS, alpha, r_cut, k_cut,r_acc, k_acc);	// NEEDS UPDATING!!!
		alphasq = alpha * alpha;

		#if VERBOSITY == 1
		
			cout << endl;
			cout << "Ewald parameters: " << endl;
			printf("	Requested Ewald accuracy  = %13.6e\n", accuracy);
			printf("	Ewald alpha               = %13.6e\n", alpha);
			printf("	R-Space Ewald cutoff      = %13.6e\n", r_cut);
			printf("	R-Space accuracy estimate = %13.6e\n", r_acc);
			printf("	K-Space accuracy estimate = %13.6e\n", k_acc);
			cout << endl;
			
		#endif
		

		called_before = true;
	}

	TRAJECTORY.TMP_EWALD.resize(TRAJECTORY.ATOMS);

	for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)
	{
		TRAJECTORY.TMP_EWALD[a1].X = 0;
		TRAJECTORY.TMP_EWALD[a1].Y = 0;
		TRAJECTORY.TMP_EWALD[a1].Z = 0;
	}

	// Main loop Ewald Coulomb energy/forces:
	
	int    curr_pair_type_idx;
	
	for(int a1=0;a1<TRAJECTORY.ATOMS;a1++)		// Double sum over atom pairs
	{
		for(int a2=0;a2<a1;a2++)
		{
			tempy=0.0;
	
			RVEC.X = TRAJECTORY.COORDS[a2].X - TRAJECTORY.COORDS[a1].X;
			RVEC.Y = TRAJECTORY.COORDS[a2].Y - TRAJECTORY.COORDS[a1].Y;
			RVEC.Z = TRAJECTORY.COORDS[a2].Z - TRAJECTORY.COORDS[a1].Z;
	
			// THIS IS THE REAL-SPACE LOOP IN MIC:
	
			RVEC.X = RVEC.X - floor( 0.5 + RVEC.X / TRAJECTORY.BOXDIM.X ) * TRAJECTORY.BOXDIM.X;
			RVEC.Y = RVEC.Y - floor( 0.5 + RVEC.Y / TRAJECTORY.BOXDIM.Y ) * TRAJECTORY.BOXDIM.Y;
			RVEC.Z = RVEC.Z - floor( 0.5 + RVEC.Z / TRAJECTORY.BOXDIM.Z ) * TRAJECTORY.BOXDIM.Z;

	      
			rlen_mi = RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z;
	  
			if ( rlen_mi < r_cut * r_cut ) 
			{
				rlen_mi = sqrt(rlen_mi);

				erfc_val = erfc(alpha * rlen_mi);
	    
				tempy += erfc_val/rlen_mi;
	    
				tempd = ( (-2.0/SQRT_PI) * exp(-1.0*alphasq*rlen_mi*rlen_mi) * rlen_mi*alpha - erfc_val ) / (rlen_mi*rlen_mi*rlen_mi);
				
				tempx =- tempd * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[a2];	

				tempx *= ke;

		        TRAJECTORY.TMP_EWALD[a1].X += RVEC.X * tempx;
				TRAJECTORY.TMP_EWALD[a1].Y += RVEC.Y * tempx;
				TRAJECTORY.TMP_EWALD[a1].Z += RVEC.Z * tempx;
				
		        TRAJECTORY.TMP_EWALD[a2].X -= RVEC.X * tempx;
				TRAJECTORY.TMP_EWALD[a2].Y -= RVEC.Y * tempx;
				TRAJECTORY.TMP_EWALD[a2].Z -= RVEC.Z * tempx;	
				

				UCoul += ke * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[a2] * tempy;

			} 
			else 
			{
				// printf("Cutoff exclusion found\n");
			}
		}
	}

	// K-Space loops.
	
	Ewald_K_Space_New(alphasq, k_cut, TRAJECTORY, TMP_UCoul);

	// Update potential energy..
	TRAJECTORY.TOT_POT_ENER += UCoul;

	for(int a1=0;a1<TRAJECTORY.ATOMS;a1++) 
	{
		TRAJECTORY.ACCEL[a1].X -= TRAJECTORY.TMP_EWALD[a1].X;
		TRAJECTORY.ACCEL[a1].Y -= TRAJECTORY.TMP_EWALD[a1].Y;
		TRAJECTORY.ACCEL[a1].Z -= TRAJECTORY.TMP_EWALD[a1].Z;
	}
	
	// Use special relation between energy and pressure of a coulomb system.
	// See Hummer et. al, JCP 109, 2791 (1998) eq. 15.
	
	TRAJECTORY.PRESSURE_XYZ += UCoul+TMP_UCoul;

	return;

}

// FUNCTION UPDATED
void optimal_ewald_params(double accuracy, double V, int nat, double &alpha, double & rc, int & kc, double & r_acc, double & k_acc)
// Calculate optimal Ewald parameters as suggested by Fincham, Mol. Sim. 1, 1-9 (1994) and the Moldy code manual.
{
  double p;
  
  double effort_ratio = 5.5;  // Ratio of time for real vs. fourier term.

  p = -log(accuracy);

  alpha = sqrt(M_PI) * pow(effort_ratio * nat/(V*V), 1.0/6.0);

  rc = 0.9 * sqrt(p)/alpha;
  double rkc = 1.7 * alpha * sqrt(p);

  double lavg = pow(V,1.0/3.0);
  kc = ceil(rkc * lavg / (2.0 * M_PI ));

  // From DLPOLY2 manual.
  k_acc = exp(-rkc * rkc / (4.0 * (alpha*alpha)) ) / (rkc * rkc);

  r_acc = erfc(alpha * rc) / rc;

}

// FUNCTION UPDATED
void ZCalc_Ewald_Deriv(FRAME & FRAME_TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, vector <vector <XYZ > > & FRAME_COULOMB_FORCES, map<string,int> & PAIR_MAP)
{
	XYZ RVEC; // Replaces  Rvec[3];
	double rlen;

	// Main loop Ewald Coulomb:
  
    string			TEMP_STR;
	int				i_pair;
	const  int 		maxk = 10000;			// Max number of kspace vectors
	static int 		kmax;					// ?
	const  int 		ksqmax = 50;			// ?
	static double 	alphasq = 0.7;			// ?
	const  double 	PI = 3.14159265359;
	static int    	totk = 0;
	int    			ksq;
	double 			Kfac;
	XYZ    			R_K; 
	double 			rksq;
	double 			alpha = sqrt(alphasq);
	static double 	r_cut;								// DO WE EVER USE STRAIGHT UP R_CUT?? IF NOT, WHY NOT JUST PASS R_CUT_SQUARED FROM THE OPTIMAL PARAMS FUNCTION??
	static double        	r_cut_squared;				// so we don't need to do an r_cut*r_cut calculation inside of a loop
	const  double 	accuracy = EWALD_ACCURACY;
	double 			tempd, tempd2, tempd3, tempd4;
	
	double 		Volume = FRAME_TRAJECTORY.BOXDIM.X * FRAME_TRAJECTORY.BOXDIM.Y * FRAME_TRAJECTORY.BOXDIM.Z;
	static bool called_before = false;
	XYZ 		D_XYZ; 		// replaces dx,dy,dz
	double 		rlen_mi;
	static vector<XYZ> SIN_XYZ; 	// Replaces sinx,y,z
	static vector<XYZ> COS_XYZ; 	// Replaces cosx,y,z
	double 		ke=1.0; //332.0637157615209;//this is the unit conversion
	//to achieve charges in nice electron units. 
	//currently we apply this conversion at MD-level, not here.

	static double *Kfac_v;
	static vector<XYZ_INT> K_V(maxk); // Replaces *kx_v, *ky_v, *kz_v;

	if (! called_before ) // NOTE: This is modifying STATIC variables.. meaning they exist even after the block exits
	{
		double r_acc, k_acc;

		optimal_ewald_params(accuracy, Volume, FRAME_TRAJECTORY.ATOMS, alpha, r_cut, kmax,r_acc, k_acc);	// NEEDS UPDATING!!!
		alphasq = alpha * alpha;
		r_cut_squared = r_cut*r_cut;
		
		#if VERBOSITY == 1
			printf("\tEwald_Deriv:\n");
			printf("\tR-Space Ewald cutoff      = %13.6e\n", r_cut);
	   	 	printf("\tR-Space accuracy estimate = %13.6e\n", r_acc);
	   	 	printf("\tK-space accuracy estimate = %13.6e\n", k_acc);		
		#endif

	    called_before = true;

    
	    Kfac_v = new double [maxk];	// DECLARING STUFF LIKE THIS IN A FUNCTION A MILLION TIMES WILL KILL EFFICIENCY... CAN WE MAKE THESE GLOBAL-ISH?		
		SIN_XYZ.resize(maxk); // Replaces sinx, siny, sinz
		COS_XYZ.resize(maxk);
  
		totk = 0;

	    double min_vfac = 1.0e8;
	
		for(int kx=0;kx<=kmax;kx++) 
		{
			for(int ky=0;ky<=kmax;ky++)
			{
				for(int kz=0;kz<=kmax;kz++)
				{
					ksq = kx*kx + ky*ky + kz*kz;
				  
					if(ksq!=0 and ksq<ksqmax)
					{

						R_K.X = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.X ) * kx;
						K_V[totk].X = kx;
	
						R_K.Y = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.Y ) * ky;
						K_V[totk].Y = ky;

						R_K.Z = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.Z ) * kz;
						K_V[totk].Z = kz;
					
						rksq = R_K.X*R_K.X + R_K.Y*R_K.Y + R_K.Z*R_K.Z;

    					// Note:  The original version of the Ewald evaluator did not work 
    					// when the lattice constants were different.
    					 
    					Kfac_v[totk] = exp(-rksq/(4.0*alphasq))/rksq;
				
						if ( Kfac_v[totk] < min_vfac ) 
							min_vfac = Kfac_v[totk];
						
						totk++;
				
						if(totk>maxk)
						{
							cout << "totk = " << totk << " greater than maxk = " << maxk << endl;
							exit(1);
						}
					}
				}
			}
		}
		
		#if VERBOSITY == 1
			cout << "	Number of Ewald K-vectors = " << totk << endl;
		#endif
	}
	
	

	
    for(int a1=0;a1<FRAME_TRAJECTORY.ATOMS;a1++) //Ewald real-space sum.
    {
		for(int a2=0;a2<a1;a2++)
		{
		  
			TEMP_STR = FRAME_TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(FRAME_TRAJECTORY.ATOMTYPE[a2]);
							
			i_pair = PAIR_MAP[TEMP_STR];		  
		  
			RVEC.X = FRAME_TRAJECTORY.COORDS[a2].X - FRAME_TRAJECTORY.COORDS[a1].X;
			RVEC.Y = FRAME_TRAJECTORY.COORDS[a2].Y - FRAME_TRAJECTORY.COORDS[a1].Y;
			RVEC.Z = FRAME_TRAJECTORY.COORDS[a2].Z - FRAME_TRAJECTORY.COORDS[a1].Z;
	  
			rlen = sqrt(RVEC.X*RVEC.X + RVEC.Y*RVEC.Y + RVEC.Z*RVEC.Z);

			D_XYZ.X = RVEC.X;
			D_XYZ.Y = RVEC.Y;
			D_XYZ.Z = RVEC.Z;
			
			//real-space sum done in M.I.C.

			D_XYZ.X = D_XYZ.X - floor( 0.5 + D_XYZ.X / FRAME_TRAJECTORY.BOXDIM.X ) * FRAME_TRAJECTORY.BOXDIM.X;	
			D_XYZ.Y = D_XYZ.Y - floor( 0.5 + D_XYZ.Y / FRAME_TRAJECTORY.BOXDIM.Y ) * FRAME_TRAJECTORY.BOXDIM.Y;
			D_XYZ.Z = D_XYZ.Z - floor( 0.5 + D_XYZ.Z / FRAME_TRAJECTORY.BOXDIM.Z ) * FRAME_TRAJECTORY.BOXDIM.Z;

			rlen_mi = D_XYZ.X*D_XYZ.X + D_XYZ.Y*D_XYZ.Y + D_XYZ.Z*D_XYZ.Z; 

			if ( rlen_mi < r_cut_squared ) 
			{
				rlen_mi = sqrt(rlen_mi);
			
				//Note that for least-squares problem, the terms are Force/(Q[a1]*Q[a2]) because Q*Q is the linear multiplier.
		    
				tempd = ( (-2.0/sqrt(PI))*exp(-1*alpha*alpha*rlen_mi*rlen_mi)*rlen_mi*alpha - erfc(alpha*rlen_mi))/(rlen_mi*rlen_mi); //tempd is the derivate of tempy--force.
				tempd *= -1.0;
				tempd *= 2.0;  // Sum a1 > a2, not all a2
				tempd *= ke;
				tempd *= 0.5 / rlen_mi;
				
		
				FRAME_COULOMB_FORCES[i_pair][a1].X += tempd * D_XYZ.X;//populate terms according to x,y,z component.
				FRAME_COULOMB_FORCES[i_pair][a1].Y += tempd * D_XYZ.Y;
				FRAME_COULOMB_FORCES[i_pair][a1].Z += tempd * D_XYZ.Z;
				
				FRAME_COULOMB_FORCES[i_pair][a2].X -= tempd * D_XYZ.X;//populate terms according to x,y,z component.
				FRAME_COULOMB_FORCES[i_pair][a2].Y -= tempd * D_XYZ.Y;
				FRAME_COULOMB_FORCES[i_pair][a2].Z -= tempd * D_XYZ.Z;					
			}
		}
	}
	
	int kx, ky, kz;
	
    for(int a1=0;a1<FRAME_TRAJECTORY.ATOMS;a1++) //Ewald K-space sum.	// -- this is where the slow down occurs
    {
		for(int a2=0 ; a2<a1 ;a2++)
		{
			TEMP_STR = FRAME_TRAJECTORY.ATOMTYPE[a1];
			TEMP_STR.append(FRAME_TRAJECTORY.ATOMTYPE[a2]);
							
			i_pair = PAIR_MAP[TEMP_STR];
				
			RVEC.X = FRAME_TRAJECTORY.COORDS[a2].X - FRAME_TRAJECTORY.COORDS[a1].X;
			RVEC.Y = FRAME_TRAJECTORY.COORDS[a2].Y - FRAME_TRAJECTORY.COORDS[a1].Y;
			RVEC.Z = FRAME_TRAJECTORY.COORDS[a2].Z - FRAME_TRAJECTORY.COORDS[a1].Z;

			// Evaluate sin factors for this pair of atoms. -- NEEDS UPDATING TO HANDLE SIN_XYZ, COS_XYZ, RVEC and FRAME_TRAJECTORY.BOXDIM !!!!
			generate_trig(SIN_XYZ, COS_XYZ, RVEC, FRAME_TRAJECTORY.BOXDIM, kmax);
			
			// Sum over all k vectors.
			for ( int ik = 0; ik < totk; ik++ ) 
			{
				R_K.X = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.X ) *  K_V[ik].X;
				R_K.Y = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.Y ) *  K_V[ik].Y;
				R_K.Z = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.Z ) *  K_V[ik].Z;

				Kfac =  Kfac_v[ik];//exp(-1.0*rksq/(4.0*alpha*alpha))/rksq;

				tempd = add_sines(K_V[ik].X, K_V[ik].Y, K_V[ik].Z, SIN_XYZ, COS_XYZ);

				tempd *= ke;
				tempd *= 2.0; // Sum a1 > a2, not all a2.
				tempd *= Kfac * 0.5*(4*PI/Volume);
				tempd3 = 0.0;

				if (  K_V[ik].Y > 0 ) 
				{
					tempd2 = add_sines(K_V[ik].X, -K_V[ik].Y, K_V[ik].Z, SIN_XYZ, COS_XYZ);

					tempd2 *= ke;
					tempd2 *= 2.0; // Sum a1 > a2, not all a2.
					tempd2 *= Kfac * 0.5*(4*PI/Volume);

					if (  K_V[ik].Z > 0 ) 
					{
						tempd3 = add_sines(K_V[ik].X, -K_V[ik].Y, -K_V[ik].Z, SIN_XYZ, COS_XYZ);
						tempd3 *= ke;
						tempd3 *= 2.0; // Sum a1 > a2, not all a2.
						tempd3 *= Kfac * 0.5*(4*PI/Volume);
					}
					else 
						tempd3 = 0.0;
					 
				} 
				else 
					tempd2 = 0.0;


				if (  K_V[ik].Z > 0 ) 
				{
					tempd4 = add_sines(K_V[ik].X, K_V[ik].Y, -K_V[ik].Z, SIN_XYZ, COS_XYZ);
					tempd4 *= ke;
					tempd4 *= 2.0; // Sum a1 > a2, not all a2.
					tempd4 *= Kfac * 0.5*(4*PI/Volume);
				}
				else 
					tempd4 = 0.0;
	  

				if (  K_V[ik].X > 0 ) 
				{
					tempd  *= 2.0; // Sum kx >= 0, not all kx;
					tempd2 *= 2.0; 
					tempd3 *= 2.0; 
					tempd4 *= 2.0; 
				}
				
//cout << tempd << " " << tempd2 << " " << tempd3 << " " << tempd4 << endl;

				FRAME_COULOMB_FORCES[i_pair][a1].X += ( tempd + tempd2 + tempd3 + tempd4 ) * R_K.X;
				FRAME_COULOMB_FORCES[i_pair][a1].Y += ( tempd - tempd2 - tempd3 + tempd4 ) * R_K.Y;
				FRAME_COULOMB_FORCES[i_pair][a1].Z += ( tempd + tempd2 - tempd3 - tempd4 ) * R_K.Z;
		
				FRAME_COULOMB_FORCES[i_pair][a2].X -= ( tempd + tempd2 + tempd3 + tempd4 ) * R_K.X;
				FRAME_COULOMB_FORCES[i_pair][a2].Y -= ( tempd - tempd2 - tempd3 + tempd4 ) * R_K.Y;
				FRAME_COULOMB_FORCES[i_pair][a2].Z -= ( tempd + tempd2 - tempd3 - tempd4 ) * R_K.Z;	

												
		  	}		
		}
	}
	
	
	
	for(int a1=0;a1<FRAME_TRAJECTORY.ATOMS;a1++)
	{
		for(int i_pair=0; i_pair<ATOM_PAIRS.size(); i_pair++)
		{				
			FRAME_COULOMB_FORCES[i_pair][a1].X *= -1.0;
			FRAME_COULOMB_FORCES[i_pair][a1].Y *= -1.0;
			FRAME_COULOMB_FORCES[i_pair][a1].Z *= -1.0;		
		}
	}
	return;
}

// FUNCTION UPDATED
static double add_sines(int kx, int ky, int kz, vector<XYZ> & SIN_XYZ, vector<XYZ> & COS_XYZ) // Add precomputed sines to get overall factor.
{
	
	static double SINYZ = 0;
	static double COSYZ = 0;
	
	XYZ SIN_K, COS_K;

	if ( kx >= 0 ) 
	{
		SIN_K.X = SIN_XYZ[kx].X;
		COS_K.X = COS_XYZ[kx].X;
	}
	else 
	{
		SIN_K.X = -SIN_XYZ[-kx].X;
		COS_K.X =  COS_XYZ[-kx].X;
	}
	if ( ky >= 0 ) 
	{
		SIN_K.Y = SIN_XYZ[ky].Y;
		COS_K.Y = COS_XYZ[ky].Y;
    }
	else 
    {
		SIN_K.Y = -SIN_XYZ[-ky].Y;
		COS_K.Y =  COS_XYZ[-ky].Y;
    }
	if ( kz >= 0 ) 
	{
		SIN_K.Z = SIN_XYZ[kz].Z;
		COS_K.Z = COS_XYZ[kz].Z;
	}
	else 
    {
		SIN_K.Z = -SIN_XYZ[-kz].Z;
		COS_K.Z =  COS_XYZ[-kz].Z;
    }

	SINYZ = SIN_K.Y * COS_K.Z + COS_K.Y * SIN_K.Z;
	COSYZ = COS_K.Y * COS_K.Z - SIN_K.Y * SIN_K.Z;
	
	return(SIN_K.X * COSYZ + COS_K.X * SINYZ);
}

// FUNCTION UPDATED
static void  generate_trig(vector<XYZ> & SIN_XYZ, vector<XYZ> & COS_XYZ, XYZ & RVEC, XYZ & BOXDIM, int kmax)
{
	SIN_XYZ.resize(kmax);
	COS_XYZ.resize(kmax);
	
	SIN_XYZ[0].X = SIN_XYZ[0].Y = SIN_XYZ[0].Z = 0.0;
	COS_XYZ[0].X = COS_XYZ[0].Y = COS_XYZ[0].Z = 1.0;

	SIN_XYZ[1].X = sin( 2.0 * M_PI * RVEC.X / BOXDIM.X );
	SIN_XYZ[1].Y = sin( 2.0 * M_PI * RVEC.Y / BOXDIM.Y );
	SIN_XYZ[1].Z = sin( 2.0 * M_PI * RVEC.Z / BOXDIM.Z );

	COS_XYZ[1].X = cos( 2.0 * M_PI * RVEC.X / BOXDIM.X );
	COS_XYZ[1].Y = cos( 2.0 * M_PI * RVEC.Y / BOXDIM.Y );
	COS_XYZ[1].Z = cos( 2.0 * M_PI * RVEC.Z / BOXDIM.Z );
  
	// Use angle addition formula recursively.
	
	for ( int ik = 2; ik <= kmax; ik++ ) 
	{
		SIN_XYZ[ik].X= SIN_XYZ[ik-1].X* COS_XYZ[1].X+ COS_XYZ[ik-1].X* SIN_XYZ[1].X;
		SIN_XYZ[ik].Y= SIN_XYZ[ik-1].Y* COS_XYZ[1].Y+ COS_XYZ[ik-1].Y* SIN_XYZ[1].Y;
		SIN_XYZ[ik].Z= SIN_XYZ[ik-1].Z* COS_XYZ[1].Z+ COS_XYZ[ik-1].Z* SIN_XYZ[1].Z;

		COS_XYZ[ik].X= COS_XYZ[ik-1].X* COS_XYZ[1].X- SIN_XYZ[ik-1].X* SIN_XYZ[1].X;
		COS_XYZ[ik].Y= COS_XYZ[ik-1].Y* COS_XYZ[1].Y- SIN_XYZ[ik-1].Y* SIN_XYZ[1].Y;
		COS_XYZ[ik].Z= COS_XYZ[ik-1].Z* COS_XYZ[1].Z- SIN_XYZ[ik-1].Z* SIN_XYZ[1].Z;
	}
}
