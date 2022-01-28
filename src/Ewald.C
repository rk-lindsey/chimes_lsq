#include<iomanip>
#include "functions.h"
#include "util.h"
#include "A_Matrix.h"

#define EWALD_ACCURACY 1.0e-06

using namespace std;

//////////////////////////////////////////
//
//	FUNCTION HEADERS
//
//////////////////////////////////////////

static double add_sines    (int kx, int ky, int kz, vector<XYZ> & SIN_XYZ, vector<XYZ> & COS_XYZ);
static void   generate_trig(vector<XYZ> & SIN_XYZ, vector<XYZ> & COS_XYZ, XYZ & RVEC, BOX & BOXDIM, int kmax);

static void Ewald_K_Space_New(double alphasq, int k_cut, FRAME & TRAJECTORY, double & UCoul, int PRIM_ATOMS, BOX & PRIM_BOX, bool lsq_mode); // MD compare force version


//////////////////////////////////////////
//
//	FUNCTION DEFINITIONS
//
//////////////////////////////////////////

static void Ewald_K_Space_New(double alphasq, int k_cut, FRAME & TRAJECTORY, double & UCoul, int PRIM_ATOMS, BOX & PRIM_BOX, bool lsq_mode)
{
// Calculate Ewald K-space components.  Use a rearrangement of the usual Ewald
// expression to generate an order-N evaluation.   See A. Y. Toukmaji et. al,
// Comp. Phys. Comm. 95, 73-92 (1996).
	
	double Volume   = PRIM_BOX.CELL_AX * PRIM_BOX.CELL_BY * PRIM_BOX.CELL_CZ;
	const double PI = M_PI;
  
	//set up Ewald Coulomb parameters:
  
	double        tempk,tempy;
	static double alpha = 0.0;
	const int     maxk =10000;
	int           kmax =10;
	int           ksq;
	double        rksq;
	static XYZ    LAST_BOXDIMS;

	//set up Kfac storage vec:
 
	static bool     called_before = false;
	static int      totk;
	static vector <double> KFAC_V(maxk);
	static vector <XYZ>    K_V(maxk);
	vector <double> sin_array(PRIM_ATOMS);
	vector <double> cos_array(PRIM_ATOMS);


	if ((!called_before) || (lsq_mode && NPROCS>1)) 
	{
		called_before = true;
		LAST_BOXDIMS.X = LAST_BOXDIMS.Y = LAST_BOXDIMS.Z = 0.0;		
	}

	if (PRIM_BOX.CELL_AX != LAST_BOXDIMS.X || PRIM_BOX.CELL_BY != LAST_BOXDIMS.Y || PRIM_BOX.CELL_CZ != LAST_BOXDIMS.Z ) 
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
						K_V[totk].X = (2.0*PI/PRIM_BOX.CELL_AX)*kx;
						K_V[totk].Y = (2.0*PI/PRIM_BOX.CELL_BY)*ky;
						K_V[totk].Z = (2.0*PI/PRIM_BOX.CELL_CZ)*kz;
						
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
							exit_run(1);
						}
					}
				}
			}
		}
		
		if ( LAST_BOXDIMS.X == 0.0  && RANK == 0 ) 
		{
			#if VERBOSITY == 1
				if(RANK==0)
				{
					cout << "	Cutting off K-vectors with magnitude less than "  << k_cut << endl;
					cout << "	Number of Ewald K-vectors = " << totk << endl;
					cout << "	K-space accuracy estimate = " << min_vfac << endl << endl;
				}
			#endif

		}
		
		LAST_BOXDIMS.X = PRIM_BOX.CELL_AX;
		LAST_BOXDIMS.Y = PRIM_BOX.CELL_BY;
		LAST_BOXDIMS.Z = PRIM_BOX.CELL_CZ;
	}

	UCoul = 0;

	// Main loop Ewald Coulomb energy/forces: K-Space loop.
	
	tempk=0;
	
	double tempsin = 0.0;
	double tempcos = 0.0;
	double kdotr;
	double tempd;

	int ikstart, ikend;
	
	// In this case, we divide work over k vectors.  The code is the same as dividing over atoms. -- MPI

	ikstart = RANK;	// 1 k-vector per process -- default for serial case
	ikend   = RANK;


	if ( lsq_mode ) 
	{
		ikstart = 0 ;
		ikend = totk - 1 ;
	}
	else
	{
		if ( NPROCS <= totk ) 
			divide_atoms(ikstart, ikend, totk);
	}


	XYZ B ;
	double kweight ;
	for(int ik= ikstart; ik <= ikend; ik++)	// Loop is MPI'd over totk
	{
		tempsin = 0.0;
		tempcos = 0.0;

		for(int a1=0; a1<PRIM_ATOMS; a1++) 
		{
			kdotr = TRAJECTORY.COORDS[a1].X * K_V[ik].X + TRAJECTORY.COORDS[a1].Y * K_V[ik].Y + TRAJECTORY.COORDS[a1].Z * K_V[ik].Z;	
			sin_array[a1] = sin(kdotr);																									
			cos_array[a1] = cos(kdotr);		
			
			tempcos += TRAJECTORY.CHARGES[a1] * cos_array[a1];																			
			tempsin += TRAJECTORY.CHARGES[a1] * sin_array[a1];		
		}
		
		kweight = KFAC_V[ik]* (tempcos*tempcos + tempsin*tempsin); 
		tempk += kweight ;

		// K space part of the pressure tensor.
		// See Heyes, PRB, 49, 755(1994), eq. 22.

		rksq = K_V[ik].X * K_V[ik].X + K_V[ik].Y * K_V[ik].Y + K_V[ik].Z * K_V[ik].Z ;

		// Diagonal part of B tensor from Heyes.
		B.X = 1.0 - 2.0 * K_V[ik].X * K_V[ik].X / rksq - 0.5 * K_V[ik].X * K_V[ik].X / alphasq ;
		B.Y = 1.0 - 2.0 * K_V[ik].Y * K_V[ik].Y / rksq - 0.5 * K_V[ik].Y * K_V[ik].Y / alphasq ;
		B.Z = 1.0 - 2.0 * K_V[ik].Z * K_V[ik].Z / rksq - 0.5 * K_V[ik].Z * K_V[ik].Z / alphasq ;

		TRAJECTORY.PRESSURE_TENSORS_XYZ_ALL[0].X += (2*PI*ke/Volume)*kweight * B.X ;
		TRAJECTORY.PRESSURE_TENSORS_XYZ_ALL[1].Y += (2*PI*ke/Volume)*kweight * B.Y ;
		TRAJECTORY.PRESSURE_TENSORS_XYZ_ALL[2].Z += (2*PI*ke/Volume)*kweight * B.Z ;
		
		for(int a1=0; a1<PRIM_ATOMS; a1++) 
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
	
	// Set up for MPI

	int a1start, a1end;
	divide_atoms(a1start, a1end, PRIM_ATOMS);

	for ( int a1=a1start; a1 <=  a1end; a1++ ) // Loop is MPI'd over TRAJECTORY.ATOMS
		UCoul -= ke * alpha /sqrt(PI) * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[a1];

	// End K-Space loop.
	
	// printf("K-space part of Ewald Energy = %13.6e\n", UCoul);

	TRAJECTORY.TOT_POT_ENER += UCoul;
	
	vector<double>().swap(sin_array);
	vector<double>().swap(cos_array);
	
}

void ZCalc_Ewald(FRAME & TRAJECTORY, JOB_CONTROL & CONTROLS, NEIGHBORS & NEIGHBOR_LIST)	// MD version
{
// Calculate Ewald interactions... 
	
  XYZ RVEC; 		// Replaces Rvec[3];
  double tempx;

  const double PI = M_PI;
  const double SQRT_PI = sqrt(PI);

  double 		tempd,tempy;
  double 		rlen_mi;
  static double alphasq;
  static double r_cut;
  static int	k_cut;
  static bool 	called_before = false;
  static double alpha;
  const double 	accuracy = EWALD_ACCURACY;
  double		UCoul = 0;
  double		TMP_UCoul;
  double 		erfc_val;
  int a2start, a2end, a2;
  bool lsq_mode ;
  string TEMP_STR;

	BOX PRIM_BOX;
	int PRIM_ATOMS;
	
	// Primitive box is the same as the original box in this implementation of layers.
	PRIM_BOX.CELL_AX = TRAJECTORY.BOXDIM.CELL_AX;
	PRIM_BOX.CELL_BY = TRAJECTORY.BOXDIM.CELL_BY;
	PRIM_BOX.CELL_CZ = TRAJECTORY.BOXDIM.CELL_CZ;
		
	PRIM_ATOMS = TRAJECTORY.ATOMS;

	if ( CONTROLS.NFRAMES > 0 ) 
		lsq_mode = true ;
	else 
		lsq_mode = false ;

	if ( ! called_before ) 
	{
		double r_acc, k_acc;
		optimal_ewald_params(accuracy, PRIM_ATOMS, alpha, r_cut, k_cut,r_acc, k_acc, PRIM_BOX);	

		// Update the neighbor list based on the Ewald cutoff.
		NEIGHBOR_LIST.EWALD_CUTOFF = r_cut;
		NEIGHBOR_LIST.DO_UPDATE(TRAJECTORY, CONTROLS);

		alphasq = alpha * alpha;

		#if VERBOSITY == 1
			
			if ( RANK == 0 ) 
			{
				cout << endl;
				cout << "Ewald parameters: " << endl;
				printf("	Requested Ewald accuracy  = %13.6e\n", accuracy);
				printf("	Ewald alpha               = %13.6e\n", alpha);
				printf("	R-Space Ewald cutoff      = %13.6e\n", r_cut);
				printf("	R-Space accuracy estimate = %13.6e\n", r_acc);
				printf("	K-Space accuracy estimate = %13.6e\n", k_acc);
				cout << endl;
			}

		#endif
			
		//	comment out line below for LAMMPS calls, in case of flexible simulation cells
		called_before = true;

	}

	TRAJECTORY.TMP_EWALD.resize(PRIM_ATOMS);
	
	// Set up for MPI
	
	int a1start, a1end;
	
	// Set up for a least squares calculation... if statement takes care of MD
	
	a1start = 0;
	a1end = PRIM_ATOMS-1;

	if ( ! lsq_mode ) // MD Calculation.
		divide_atoms(a1start, a1end, PRIM_ATOMS);

	

	for(int a1=0;a1<PRIM_ATOMS;a1++)
	{
		TRAJECTORY.TMP_EWALD[a1].X = 0;
		TRAJECTORY.TMP_EWALD[a1].Y = 0;
		TRAJECTORY.TMP_EWALD[a1].Z = 0;
	}

	// Scaling for permutations of atoms when the cutoff is >= box size.
	double perm_scale = NEIGHBOR_LIST.PERM_SCALE[2] ;

	// Main loop Ewald Coulomb energy/forces:
	for(int a1 = a1start; a1 <= a1end; a1++)	// Double sum over atom pairs (outer loop is MPI'd over TRAJECTORY.ATOMS)
	{
		a2start = 0;
		a2end   = NEIGHBOR_LIST.LIST_EWALD[a1].size();

		for(int a2idx=a2start; a2idx<a2end; a2idx++)	
		{
			a2 = NEIGHBOR_LIST.LIST_EWALD[a1][a2idx];
			int fidx_a2 = TRAJECTORY.PARENT[a2];
			
			tempy=0.0;
	
			// THIS IS THE REAL-SPACE LOOP IN MIC:
	
			rlen_mi = get_dist(TRAJECTORY, RVEC, a1, a2);

			if ( rlen_mi < r_cut ) 
			{
				erfc_val = erfc(alpha * rlen_mi);
	    
				tempy += erfc_val/rlen_mi;
	    
				tempd = ( (-2.0/SQRT_PI) * exp(-1.0*alphasq*rlen_mi*rlen_mi) * rlen_mi*alpha - erfc_val ) / (rlen_mi*rlen_mi*rlen_mi);
				
				tempx =- tempd * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[fidx_a2];	

				tempx *= ke;

				tempx *= perm_scale ;
				tempy *= perm_scale ;
				
				TRAJECTORY.TMP_EWALD[a1].X += RVEC.X * tempx;
				TRAJECTORY.TMP_EWALD[a1].Y += RVEC.Y * tempx;
				TRAJECTORY.TMP_EWALD[a1].Z += RVEC.Z * tempx;

				// Real space part of the pressure tensor.
				// See Heyes, PRB, 49, 755(1994), eq. 22.
				TRAJECTORY.PRESSURE_TENSORS_XYZ_ALL[0].X += tempx * RVEC.X * RVEC.X ;
				TRAJECTORY.PRESSURE_TENSORS_XYZ_ALL[1].Y += tempx * RVEC.Y * RVEC.Y ;
				TRAJECTORY.PRESSURE_TENSORS_XYZ_ALL[2].Z += tempx * RVEC.Z * RVEC.Z ;
				
				TRAJECTORY.TMP_EWALD[fidx_a2].X -= RVEC.X * tempx;
				TRAJECTORY.TMP_EWALD[fidx_a2].Y -= RVEC.Y * tempx;
				TRAJECTORY.TMP_EWALD[fidx_a2].Z -= RVEC.Z * tempx;	
				
				UCoul += ke * TRAJECTORY.CHARGES[a1] * TRAJECTORY.CHARGES[fidx_a2] * tempy;
			} 
		}
	}

	// K-Space loops.

	Ewald_K_Space_New(alphasq, k_cut, TRAJECTORY, TMP_UCoul, PRIM_ATOMS, PRIM_BOX, lsq_mode);

	// Update potential energy..
	TRAJECTORY.TOT_POT_ENER += UCoul;

	for(int a1=0;a1<PRIM_ATOMS;a1++) 
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

void optimal_ewald_params(double accuracy, int nat, double &alpha, double & rc, int & kc, double & r_acc, double & k_acc, BOX boxdim)
{
	// Calculate optimal Ewald parameters as suggested by Fincham, Mol. Sim. 1, 1-9 (1994) and the Moldy code manual.
    double p;
  
    double effort_ratio = 5.5;  // Ratio of time for real vs. fourier term.
    double V;
    double min_boxdim = boxdim.CELL_AX;

    // balance_factor and accuracy_factor determined so that K_acc and r_acc are approximately equal
    // to the target accuray.
    double balance_factor = 1.1;
    double accuracy_factor = 0.8;

    if ( boxdim.CELL_BY < min_boxdim ) 
  	  min_boxdim = boxdim.CELL_BY;
    if ( boxdim.CELL_CZ < min_boxdim ) 
  	  min_boxdim = boxdim.CELL_CZ;

    V = boxdim.CELL_AX * boxdim.CELL_BY * boxdim.CELL_CZ;

    p = -log(accuracy) * accuracy_factor;
    alpha = sqrt(M_PI) * pow(effort_ratio * nat/(V*V), 1.0/6.0);

    // Try optimal cutoff estimate.
    rc = balance_factor * sqrt(p)/alpha;

    if ( rc > 0.5 * min_boxdim ) {
  	  // Optimal cutoff is too big.  Reduce cutoff.
  	  rc = 0.5 * min_boxdim;
  	  alpha = sqrt(p) / rc;
    }

    // rkc is the K cutoff with units 1/r.
    double rkc = 2.0 * alpha * sqrt(p) / balance_factor;
    double lavg = pow(V,1.0/3.0);

    // kc is the interger k cutoff.
    kc = ceil(rkc * lavg / (2.0 * M_PI));

    // Re-adjust rkc for error estimates to integer cutoff value.
    rkc = kc * 2.0 * M_PI / lavg;

    // From Fincham, eq. 1-2.
    k_acc = exp(-rkc * rkc / (4.0 * (alpha*alpha)) ) / (rkc * rkc);
  
    r_acc = erfc(alpha * rc) / rc;
}

void ZCalc_Ewald_Deriv(FRAME & FRAME_TRAJECTORY, vector<PAIRS> & ATOM_PAIRS, A_MAT & A_MATRIX, map<string,int> & PAIR_MAP,NEIGHBORS & NEIGHBOR_LIST, JOB_CONTROL & CONTROLS)
{
	 XYZ RVEC; // Replaces  Rvec[3];

	 // Main loop Ewald Coulomb:
  
	 string			TEMP_STR;
	 int			i_pair;
	 const  int 		maxk     = 10000;	// Max number of kspace vectors
	 static int 		kmax;					
	 const  int 		ksqmax   = 50;			
	 static double 		alphasq  = 0.7;			
	 const  double 		PI       = 3.14159265359;
	 static int    		totk     = 0;
	 int    			ksq;
	 double 			Kfac;
	 XYZ    			R_K; 
	 double 			rksq;
	 double 			alpha    = sqrt(alphasq);
	 static double 		r_cut;						
	 const  double 		accuracy = EWALD_ACCURACY;
	 double 			tempd, tempd2, tempd3, tempd4;
	
	 double 			Volume;
	 XYZ 			D_XYZ; 						
	 double 			rlen_mi;
	 static vector<XYZ>	SIN_XYZ; 				
	 static vector<XYZ>	COS_XYZ; 				
	 double 			ke = 1.0;	// this is the unit conversion to achieve charges in nice electron units. currently we apply this conversion at MD-level, not here.
	 static double 		*Kfac_v;
	 static vector<XYZ_INT>	K_V(maxk); 
	 static XYZ    		LAST_BOXDIMS;
	 int			a2start, a2end, a2;
	
	 Volume = FRAME_TRAJECTORY.BOXDIM.VOL;
	
	 bool BOX_CHANGED = false;
	
	 if (LAST_BOXDIMS.X != FRAME_TRAJECTORY.BOXDIM.CELL_AX  || LAST_BOXDIMS.Y != FRAME_TRAJECTORY.BOXDIM.CELL_BY  || LAST_BOXDIMS.Z != FRAME_TRAJECTORY.BOXDIM.CELL_CZ)
			BOX_CHANGED = true;
	 if(NPROCS>1)
			BOX_CHANGED = true;
			
	 if (BOX_CHANGED) // NOTE: This is modifying STATIC variables.. meaning they exist even after the block exits
	 {
			double r_acc, k_acc;

			LAST_BOXDIMS.X = FRAME_TRAJECTORY.BOXDIM.CELL_AX ;
			LAST_BOXDIMS.Y = FRAME_TRAJECTORY.BOXDIM.CELL_BY ;
			LAST_BOXDIMS.Z = FRAME_TRAJECTORY.BOXDIM.CELL_CZ ;

			optimal_ewald_params(accuracy, FRAME_TRAJECTORY.ATOMS, alpha, r_cut, kmax,r_acc, k_acc, FRAME_TRAJECTORY.BOXDIM);	

			// Update the neighbor list based on the Ewald cutoff.
			NEIGHBOR_LIST.EWALD_CUTOFF = r_cut ;
			NEIGHBOR_LIST.DO_UPDATE(FRAME_TRAJECTORY,CONTROLS) ;
		
			alphasq = alpha * alpha;
		
#if VERBOSITY == 1
			if ( RANK == 0 ) 
			{
				 printf("\tEwald_Deriv:\n");
				 printf("\tR-Space Ewald cutoff      = %13.6e\n", r_cut);
				 printf("\tR-Space accuracy estimate = %13.6e\n", r_acc);
				 printf("\tK-space accuracy estimate = %13.6e\n", k_acc);		
			}
#endif

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

									R_K.X = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.CELL_AX ) * kx;
									K_V[totk].X = kx;
	
									R_K.Y = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.CELL_BY ) * ky;
									K_V[totk].Y = ky;

									R_K.Z = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.CELL_CZ ) * kz;
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
										 exit_run(1);
									}
							 }
						}
				 }
			}
		
#if VERBOSITY == 1
			if ( RANK == 0 ) cout << "	Number of Ewald K-vectors = " << totk << endl;
#endif
	 }

	 // Scaling for permutations of atoms when the cutoff is >= box size.
	 double perm_scale = NEIGHBOR_LIST.PERM_SCALE[2] ;
	
	 for(int a1=0;a1<FRAME_TRAJECTORY.ATOMS;a1++) //Ewald real-space sum.
	 {
			a2start = 0;
			a2end   = NEIGHBOR_LIST.LIST_EWALD[a1].size();

			for(int a2idx=a2start; a2idx<a2end; a2idx++)	
			{
				 a2 = NEIGHBOR_LIST.LIST_EWALD[a1][a2idx] ;
				 int fidx_a2 = FRAME_TRAJECTORY.PARENT[a2] ;
		  
				 TEMP_STR = FRAME_TRAJECTORY.ATOMTYPE[a1];
				 TEMP_STR.append(FRAME_TRAJECTORY.ATOMTYPE[fidx_a2]);
							
				 i_pair = PAIR_MAP[TEMP_STR];		  

				 rlen_mi = get_dist(FRAME_TRAJECTORY, RVEC, a1, a2) ;

				 D_XYZ.X = RVEC.X;
				 D_XYZ.Y = RVEC.Y;
				 D_XYZ.Z = RVEC.Z;
			

				 if ( rlen_mi < r_cut ) 
				 {
						//Note that for least-squares problem, the terms are Force/(Q[a1]*Q[a2]) because Q*Q is the linear multiplier.
		    
						tempd = ( (-2.0/sqrt(PI))*exp(-1*alpha*alpha*rlen_mi*rlen_mi)*rlen_mi*alpha - erfc(alpha*rlen_mi))/(rlen_mi*rlen_mi); //tempd is the derivative of tempy--force.
						tempd *= -1.0;
						tempd *= 2.0;  // Sum a1 > a2, not all a2
						tempd *= ke;
						tempd *= 0.5 / rlen_mi;

						tempd *= perm_scale ;
				
						A_MATRIX.CHARGES[i_pair][a1].X += tempd * D_XYZ.X;//populate terms according to x,y,z component.
						A_MATRIX.CHARGES[i_pair][a1].Y += tempd * D_XYZ.Y;
						A_MATRIX.CHARGES[i_pair][a1].Z += tempd * D_XYZ.Z;
				
						A_MATRIX.CHARGES[i_pair][fidx_a2].X -= tempd * D_XYZ.X;//populate terms according to x,y,z component.
						A_MATRIX.CHARGES[i_pair][fidx_a2].Y -= tempd * D_XYZ.Y;
						A_MATRIX.CHARGES[i_pair][fidx_a2].Z -= tempd * D_XYZ.Z;					
				 }
			}
	 }
	
	 for(int a1=0;a1<FRAME_TRAJECTORY.ATOMS;a1++) //Ewald K-space sum.	// -- this is where the slow down occurs
	 {
			for(int a2=0; a2<a1;a2++)
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
						R_K.X = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.CELL_AX ) *  K_V[ik].X;
						R_K.Y = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.CELL_BY ) *  K_V[ik].Y;
						R_K.Z = ( 2.0 * PI / FRAME_TRAJECTORY.BOXDIM.CELL_CZ ) *  K_V[ik].Z;

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

						A_MATRIX.CHARGES[i_pair][a1].X += ( tempd + tempd2 + tempd3 + tempd4 ) * R_K.X;
						A_MATRIX.CHARGES[i_pair][a1].Y += ( tempd - tempd2 - tempd3 + tempd4 ) * R_K.Y;
						A_MATRIX.CHARGES[i_pair][a1].Z += ( tempd + tempd2 - tempd3 - tempd4 ) * R_K.Z;
		
						A_MATRIX.CHARGES[i_pair][a2].X -= ( tempd + tempd2 + tempd3 + tempd4 ) * R_K.X;
						A_MATRIX.CHARGES[i_pair][a2].Y -= ( tempd - tempd2 - tempd3 + tempd4 ) * R_K.Y;
						A_MATRIX.CHARGES[i_pair][a2].Z -= ( tempd + tempd2 - tempd3 - tempd4 ) * R_K.Z;	

												
				 }		
			}
	 }	
	
	 for(int a1=0;a1<FRAME_TRAJECTORY.ATOMS;a1++)
	 {
			for(int i_pair=0; i_pair<ATOM_PAIRS.size(); i_pair++)
			{				
				 A_MATRIX.CHARGES[i_pair][a1].X *= -1.0;
				 A_MATRIX.CHARGES[i_pair][a1].Y *= -1.0;
				 A_MATRIX.CHARGES[i_pair][a1].Z *= -1.0;		
			}
	 }
	 return;
}

static double add_sines(int kx, int ky, int kz, vector<XYZ> & SIN_XYZ, vector<XYZ> & COS_XYZ) 
{
// Add precomputed sines to get overall factor.	
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

static void  generate_trig(vector<XYZ> & SIN_XYZ, vector<XYZ> & COS_XYZ, XYZ & RVEC, BOX & BOXDIM, int kmax)
{
	SIN_XYZ.resize(kmax);
	COS_XYZ.resize(kmax);
	
	SIN_XYZ[0].X = SIN_XYZ[0].Y = SIN_XYZ[0].Z = 0.0;
	COS_XYZ[0].X = COS_XYZ[0].Y = COS_XYZ[0].Z = 1.0;

	SIN_XYZ[1].X = sin( 2.0 * M_PI * RVEC.X / BOXDIM.CELL_AX );
	SIN_XYZ[1].Y = sin( 2.0 * M_PI * RVEC.Y / BOXDIM.CELL_BY );
	SIN_XYZ[1].Z = sin( 2.0 * M_PI * RVEC.Z / BOXDIM.CELL_CZ );

	COS_XYZ[1].X = cos( 2.0 * M_PI * RVEC.X / BOXDIM.CELL_AX );
	COS_XYZ[1].Y = cos( 2.0 * M_PI * RVEC.Y / BOXDIM.CELL_BY );
	COS_XYZ[1].Z = cos( 2.0 * M_PI * RVEC.Z / BOXDIM.CELL_CZ );
  
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

void SubtractEwaldForces(FRAME & SYSTEM, NEIGHBORS & NEIGHBOR_LIST, JOB_CONTROL & CONTROLS)
// Subtract Coulomb Forces from SYSTEM.FORCES.  Used in fitting a potential with pre-determined charges.
{
	for(int a=0;a<SYSTEM.ATOMS;a++)
	{
		SYSTEM.ACCEL[a].X = 0;
		SYSTEM.ACCEL[a].Y = 0;
		SYSTEM.ACCEL[a].Z = 0;
	}

	ZCalc_Ewald(SYSTEM, CONTROLS, NEIGHBOR_LIST);

	for(int a=0;a<SYSTEM.ATOMS;a++)
	{
		SYSTEM.FORCES[a].X -= SYSTEM.ACCEL[a].X;
		SYSTEM.FORCES[a].Y -= SYSTEM.ACCEL[a].Y;
		SYSTEM.FORCES[a].Z -= SYSTEM.ACCEL[a].Z;
	}
}
