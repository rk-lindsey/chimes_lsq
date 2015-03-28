#include "functions.h"

#define EWALD_ACCURACY 1.0e-07

static void Ewald_K_Space_New(double alphasq, int k_cut,
			      double **Coord, double *Q, double *Latcons, 
			      int nat, double **FCoul,double& Vtot) ;

static void Ewald_K_Space_Orig(double alphasq, double **Coord,string *Lb, double *Latcons,
			       const int nat,double **SForce,double& Vtot) ;

void ZCalc_Ewald(double **Coord, string *Lb, double *Q, double *Latcons,const int nlayers,
	   const int nat,const double smin,const double smax,
	   const double sdelta,const int snum, 
	   double *params,double **SForce,double& Vtot,double& Pxyz)
// Calculate Ewald interactions.
{
  double Rvec[3];
  double tempx ;

  const double PI = M_PI ;
  const double SQRT_PI = sqrt(PI) ;
  //set up Ewald Coulomb parameters:
  double **FCoul ;
  
  FCoul = new double* [nat] ;
  for(int a=0;a<nat;a++)
    FCoul[a] = new double[3] ;

  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      FCoul[a1][c]=0.0 ;

  double tempd,tempy;
  double dx,dy,dz;
  double rlen_mi;
  static double alphasq ;
  int totk;
  static double r_cut ;
  static int k_cut ;
  bool if_old_k_space = false ;
  static bool called_before = false ;
  static double alpha ;
  const double accuracy = EWALD_ACCURACY ;

  if ( ! called_before ) {
    double vol = (Latcons[0] * Latcons[1] * Latcons[2]) ;
    double r_acc, k_acc ;

    optimal_ewald_params(accuracy, vol, nat, alpha,
			 r_cut, k_cut, r_acc, k_acc) ;
    alphasq = alpha * alpha ;
      
    printf("Requested Ewald accuracy  = %13.6e\n", accuracy) ;
    printf("Ewald alpha               = %13.6e\n", alpha) ;
    printf("R-Space Ewald cutoff      = %13.6e\n", r_cut) ;
    printf("R-Space accuracy estimate = %13.6e\n", 
	   r_acc) ;
    printf("K-Space accuracy estimate = %13.6e\n", 
	   k_acc) ;

    called_before = true ;
  }

  double UCoul;
  UCoul=0;
  totk=0;
  
  ////main loop Ewald Coulomb energy/forces:

  for(int a1=0;a1<nat;a1++) 
    {
      for(int a2=0; a2< a1;a2++)
	{
	  tempy=0.0;
	
	  for(int c=0;c<3;c++)
	    Rvec[c]=Coord[a2][c]-Coord[a1][c];
	
	  ////THIS IS THE REAL-SPACE LOOP IN MIC:
	
	  dx=Rvec[0];
	  dy=Rvec[1];
	  dz=Rvec[2];
	    
	  dx=dx-floor(0.5+dx/Latcons[0])*Latcons[0];
	  dy=dy-floor(0.5+dy/Latcons[1])*Latcons[1];
	  dz=dz-floor(0.5+dz/Latcons[2])*Latcons[2];
	      
	  rlen_mi=dx*dx+dy*dy+dz*dz;
	  
	  if ( rlen_mi < r_cut * r_cut ) 
	    {

	      rlen_mi = sqrt(rlen_mi) ;

	      double erfc_val = erfc(alpha * rlen_mi) ;
	    
	      tempy+=erfc_val/rlen_mi;
	    
	      tempd=( (-2.0/SQRT_PI)*exp(-1.0*alphasq*rlen_mi*rlen_mi)*rlen_mi*alpha - erfc_val)/(rlen_mi*rlen_mi*rlen_mi);

	      tempx=-tempd*Q[a1]*Q[a2];
	      tempx*=ke;

	      FCoul[a1][0]+=dx*tempx;
	      FCoul[a1][1]+=dy*tempx;
	      FCoul[a1][2]+=dz*tempx;
	    
	      FCoul[a2][0]-=dx*tempx;
	      FCoul[a2][1]-=dy*tempx;
	      FCoul[a2][2]-=dz*tempx;

	      UCoul+=ke*Q[a1]*Q[a2]*tempy;
	    } 
	  else 
	    {
	      // printf("Cutoff exclusion found\n") ;
	    }
	}
    }

  // printf("R-space part of coulomb energy = %13.6e\n", UCoul) ;

  // K-Space loop.
  if ( if_old_k_space ) {
    Ewald_K_Space_Orig(alphasq,Coord,Lb, Latcons,nat,SForce,UCoul) ;
  } else {
    Ewald_K_Space_New(alphasq, k_cut, Coord, Q, Latcons, nat, FCoul, UCoul) ;
  }

  // Update potential energy
  Vtot+=UCoul;

  ////Add the Ewald Coulomb forces:
  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      SForce[a1][c]-=FCoul[a1][c];

  /// Use special relation between energy and pressure of a coulomb system.
  /// See Hummer et. al, JCP 109, 2791 (1998) eq. 15.
  Pxyz += UCoul ;

  return;
}


static void Ewald_K_Space_New(double alphasq, int k_cut,
			      double **Coord, double *Q, double *Latcons, int nat,
			      double **FCoul,double& Vtot)
// Calculate Ewald K-space components.  Use a rearrangement of the usual Ewald
// expression to generate an order-N evaluation.   See A. Y. Toukmaji et. al,
// Comp. Phys. Comm. 95, 73-92 (1996).
{
  double Volume=Latcons[0]*Latcons[1]*Latcons[2];

  const double PI= M_PI ;
  //set up Ewald Coulomb parameters:
  double tempk,tempy;
  static double alpha = 0.0 ;
  const int maxk=10000;
  int kmax=10;
  int ksq;
  double Kfac;
  double rkx,rky,rkz;
  double rksq;
  static double last_latcons[3] = {0.0,0.0,0.0} ;

  // printf("Value of cutoff = %13.6e\n", erfc(alpha * r_cut) ) ;

  //set up Kfac storage vec:
  static bool called_before = false ;
  static double *Kfac_v ;
  static double *kx_v, *ky_v, *kz_v ;
  static int totk ;

  double *sin_array = new double [nat] ;
  double *cos_array = new double [nat] ;

  if ( ! called_before ) 
    {
      called_before = true ;

      Kfac_v = new double [maxk] ;
      kx_v   = new double [maxk] ;
      ky_v   = new double [maxk] ;
      kz_v   = new double [maxk] ;

    }

  if ( Latcons[0] != last_latcons[0] || Latcons[1] != last_latcons[1] ||
       Latcons[2] != last_latcons[2] ) 
    {
    // Update K factors when box dimensions change.
    //
      double min_vfac = 1.0e8 ;
      totk = 0 ;
      alpha = sqrt(alphasq) ;
      int ksqmax = k_cut * k_cut ;
      kmax = k_cut ;

      for(int kx=-kmax;kx<=kmax;kx++) 
	  for(int ky=-kmax;ky<=kmax;ky++)
	    for(int kz=-kmax;kz<=kmax;kz++)
	      {
		ksq=kx*kx+ky*ky+kz*kz;
		if(ksq!=0 and ksq<ksqmax)
		  {
		    rkx=(2.0*PI/Latcons[0])*kx;
		    kx_v[totk] = rkx ;

		    rky=(2.0*PI/Latcons[1])*ky;
		    ky_v[totk] = rky ;

		    rkz=(2.0*PI/Latcons[2])*kz;
		    kz_v[totk] = rkz ;

		    rksq=rkx*rkx+rky*rky+rkz*rkz;
		    // Note:  The original version of the Ewald evaluator did not work 
		    // when the lattice constants were different.
		    Kfac_v[totk]=exp(-rksq/(4.0*alphasq))/rksq;
		    if ( Kfac_v[totk] < min_vfac ) min_vfac = Kfac_v[totk] ;
		    totk++ ;
		    if(totk>maxk)
		      {
			cout<<"totk="<<totk<<" greater than maxk="<<maxk<<endl;
			exit(1);
		      }
		  }
	      }
      if ( last_latcons[0] == 0.0 ) {
	cout << "Cutting off K-vectors with magnitude less than " 
	     << k_cut << endl ;
	cout << "Number of Ewald K-vectors = " << totk << endl ;
	cout << "K-space accuracy estimate = " << min_vfac << endl ;
      }
      //for(int a1=0;a1<nat;a1++) {
      //printf("Q[%d] = %13.7e\n",a1, Q[a1]) ;
      //}
      for ( int i = 0 ; i < 3 ; i++ ) 
	{
	  last_latcons[i] = Latcons[i] ;
	}
    }

  double UCoul;
  UCoul=0;
  
  ////main loop Ewald Coulomb energy/forces:

  // K-Space loop.
  tempk=0;
  for(int ik=-0;ik<totk;ik++)
      {
	rkx = kx_v[ik] ;
	rky = ky_v[ik] ;
	rkz = kz_v[ik] ;


	Kfac=Kfac_v[ik];//exp(-1.0*rksq/(4.0*alpha*alpha))/rksq;

	double tempsin=0.0 ;
	double tempcos=0.0 ;

	for(int a1=0;a1<nat;a1++) 
	  {
	    double kdotr = Coord[a1][0] * rkx + Coord[a1][1] * rky + Coord[a1][2] * rkz ;
	    sin_array[a1] = sin(kdotr) ;
	    cos_array[a1] = cos(kdotr) ;
	    tempcos += Q[a1] * cos_array[a1] ;
	    tempsin += Q[a1] * sin_array[a1] ;
	  }
	tempk += Kfac * (tempcos * tempcos + tempsin * tempsin) ;
	for(int a1=0;a1<nat;a1++) 
	  {
	    double tempd = -4.0 * PI * Kfac * ke * Q[a1] * ( tempcos * sin_array[a1] - tempsin * cos_array[a1] ) / Volume ;
	    FCoul[a1][0] += rkx * tempd ;
	    FCoul[a1][1] += rky * tempd ;
	    FCoul[a1][2] += rkz * tempd ;
	  }
      }
  tempy=(2*PI/Volume)*tempk;
  UCoul+=ke*tempy;

  // Constant energy term.
  for ( int a1 = 0 ; a1 < nat ; a1++ ) 
    {
      UCoul -= ke * alpha /sqrt(PI) * Q[a1] * Q[a1] ;
    }

  ////end K-Space loop.
  ////end main loop Ewald Coulomb forces:
  // printf("K-space part of Ewald Energy = %13.6e\n", UCoul) ;
  Vtot+=UCoul;

  delete sin_array ;
  delete cos_array ;
}




static void Ewald_K_Space_Orig(double alphasq, double **Coord,string *Lb, double *Latcons,
			const int nat,double **SForce,double& Vtot)
{
  double Rvec[3];
  double rlen ;
  double tempx ;
  double Volume=Latcons[0]*Latcons[1]*Latcons[2];
  const double PI=3.14159265359;

  cout.precision(15);
  
  double p[35];
  for(int n=0;n<35;n++)
    p[n]=0.0;

  float pow_rlen2[35];
  float pow_rlen2_smax[35];
  for(int n=0;n<35;n++)
    {
      pow_rlen2[n]=0.0;
      pow_rlen2_smax[n]=0.0;
    }



  /*
  //This defines original Stillinger parameters for running Stillinger-MD.
  int npower;
 int p_oo_len=3;
  double p_oo[p_oo_len];
  int p_oh_len=6;
  double p_oh[p_oh_len];
  int p_hh_len=7;
  double p_hh[p_hh_len];
  
  //linear parameters:
  p_oo[0]=144.538;
  p_oo[1]=23401.9; 
  p_oh[0]=72.269;
  p_oh[1]=2.6677;
  p_oh[3]=6.0;
  p_hh[0]=36.1345;
  p_hh[1]=20.0;  
  p_hh[4]=17.03002;
  //non-linear parameters:
  p_oo[2]=8.3927;
  p_oh[2]=14.97;
  p_oh[4]=5.49305;
  p_oh[5]=2.2;
  p_hh[2]=40.0;
  p_hh[3]=2.0;
  p_hh[5]=7.60626;
  p_hh[6]=1.4525;
 
  //double qoo=144.538;
  //double qoh=72.269;
  //double qhh=36.1345;
  */

  //these are short-ranged cutoffs for analytical potential.
  //e.g. for r<oo_cut, force(r)=force(oo_cut). The values are 
  //at least 0.2 Angs below any sampled by MD; they are there 
  //in case trajectories reach them.

  ////spline charges Stillinger units
  //double qoo=245.794552892;
  //double qoh=110.774212519;
  //double qhh=49.3255746617;

  //spline charges electron units
  // double qoo=0.86035*0.86035;
  // double qoh=0.86035*0.430175;
  // double qhh=0.430175*0.430175;

  double qhh = 4.09174173723 / ke ;
  double qoo = 4.0 * qhh ;
  double qoh = sqrt(qoo * qhh) ; 

  //set up Ewald Coulomb parameters:
  double FCoul[nat][3];
  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      FCoul[a1][c]=0;

  double tempk,tempd,tempy;
  const double alpha=sqrt(alphasq);
  const int maxk=1000000;
  const int kmax=10;
  const int ksqmax=50;
  int totk;
  int ksq;
  double Kfac;
  double rkx,rky,rkz;
  double rksq;
  
  //set up Kfac storage vec:
  double Kfac_v[ksqmax];
  for(int i=0;i<ksqmax;i++)
    Kfac_v[i]=0.0;

  for(int kx=0;kx<=kmax;kx++)
    for(int ky=-kmax;ky<=kmax;ky++)
      for(int kz=-kmax;kz<=kmax;kz++)
	{
	  ksq=kx*kx+ky*ky+kz*kz;
	  if(ksq!=0 and ksq<ksqmax)
	    {
	      rkx=(2.0*PI/Latcons[0])*kx;
	      rky=(2.0*PI/Latcons[1])*ky;
	      rkz=(2.0*PI/Latcons[2])*kz;
	      rksq=rkx*rkx+rky*rky+rkz*rkz;
	      Kfac_v[ksq]=exp(-1.0*rksq/(4.0*alpha*alpha))/rksq;
	    }
	}
  

  double UCoul;
  UCoul=0;
  totk=0;
  
  double Q[nat];

  ////main loop Ewald Coulomb energy/forces:
  for(int a1=0;a1<nat;a1++)
    for(int a2=0;a2<nat;a2++)
      {
	tempy=0.0;
	
	for(int c=0;c<3;c++)
	  Rvec[c]=Coord[a2][c]-Coord[a1][c];
	
	rlen=sqrt(Rvec[0]*Rvec[0]+Rvec[1]*Rvec[1]+Rvec[2]*Rvec[2]);
	
	if(Lb[a1]=="O" and Lb[a2]=="O")
	  {
	    Q[a1]=sqrt(qoo);
	    Q[a2]=sqrt(qoo);
	  }
	else if((Lb[a1]=="O" and Lb[a2]=="H") or (Lb[a1]=="H" and Lb[a2]=="O"))
	  {
	    Q[a1]=sqrt(qoh)*-1.0;
	    Q[a2]=sqrt(qoh);
	  }
	else if(Lb[a1]=="H" and Lb[a2]=="H")
	  {
	    Q[a1]=sqrt(qhh);
	    Q[a2]=sqrt(qhh);
	  }
	printf("Q[%d] = %13.7e Q[%d] = %13.7e\n", a1, Q[a1], a2, Q[a2]) ;


	if(a1==a2)
	  tempy-=(2*alpha)/sqrt(PI);

	tempk=0;
	totk=0;

	for(int kx=0;kx<=kmax;kx++)
	//for(int kx=-kmax;kx<=kmax;kx++)
	{
	    rkx=(2.0*PI/Latcons[0])*kx;
	    for(int ky=-kmax;ky<=kmax;ky++)
	      {
		rky=(2.0*PI/Latcons[1])*ky;
		for(int kz=-kmax;kz<=kmax;kz++)
		  {
		    rkz=(2.0*PI/Latcons[2])*kz;
		    ksq=kx*kx+ky*ky+kz*kz;
		    
		    if(ksq!=0 and ksq<ksqmax)
		      {
			totk++;
			if(totk>maxk)
			  {
			    cout<<"totk="<<totk<<" greater than maxk="<<maxk<<endl;
			    exit(1);
			  }
			rksq=rkx*rkx+rky*rky+rkz*rkz;
			Kfac=Kfac_v[ksq];//exp(-1.0*rksq/(4.0*alpha*alpha))/rksq;
			
			//this is because we want to get a 1/2 speedup by using
			//only 0<kmax for x-coordinate.
			if(kx==0)
			  tempk+=Kfac*cos(rkx*Rvec[0]+rky*Rvec[1]+rkz*Rvec[2]);
			if(kx!=0)
			  tempk+=2*Kfac*cos(rkx*Rvec[0]+rky*Rvec[1]+rkz*Rvec[2]);

			tempd=sin(rkx*Rvec[0]+rky*Rvec[1]+rkz*Rvec[2]);
									
			tempx=Kfac*tempd*Q[a1]*Q[a2]*0.5*(4*PI/Volume);
			tempx*=ke;

			FCoul[a1][0]+=rkx*tempx;
			FCoul[a1][1]+=rky*tempx;
			FCoul[a1][2]+=rkz*tempx;
			
			FCoul[a2][0]-=rkx*tempx;
			FCoul[a2][1]-=rky*tempx;
			FCoul[a2][2]-=rkz*tempx;

			if(kx != 0)
			  {
			    
			    FCoul[a1][0]+=rkx*tempx;
			    FCoul[a1][1]+=rky*tempx;
			    FCoul[a1][2]+=rkz*tempx;
			
			    FCoul[a2][0]-=rkx*tempx;
			    FCoul[a2][1]-=rky*tempx;
			    FCoul[a2][2]-=rkz*tempx;
			    
			  }

		      }
		  }
	      }
	}
	

	tempy+=(4*PI/Volume)*tempk;
	UCoul+=ke*0.5*Q[a1]*Q[a2]*tempy;
      }
  ////end main loop Ewald Coulomb forces:
  // printf("K-space part of Ewald Energy = %13.6e\n", UCoul) ;

  Vtot+=UCoul;

  ////Add the Ewald Coulomb forces:
  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      SForce[a1][c]-=FCoul[a1][c];

  return;
}


void optimal_ewald_params(double accuracy, 
			  double V, int nat, double &alpha, 
			  double &rc, int &kc, double &r_acc,
			  double &k_acc)
// Calculate optimal Ewald parameters as suggested by Fincham,
// Mol. Sim. 1, 1-9 (1994) and the Moldy code manual.
{
  double p ;
  
  double effort_ratio = 5.5 ;  // Ratio of time for real vs. fourier term.

  p = -log(accuracy) ;

  alpha = sqrt(M_PI) * pow(effort_ratio * nat/(V*V), 1.0/6.0) ;

  rc = 0.9 * sqrt(p)/alpha ;
  double rkc = 1.7 * alpha * sqrt(p) ;

  double lavg = pow(V,1.0/3.0) ;
  kc = ceil(rkc * lavg / (2.0 * M_PI )) ;

  // From DLPOLY2 manual.
  k_acc = exp(-rkc * rkc / (4.0 * (alpha*alpha)) ) / (rkc * rkc) ;

  r_acc = erfc(alpha * rc) / rc ;

}


void ZCalc_Ewald_Deriv(double **Coord, string *Lb, 
		       double *Latcons,const int nlayers,
		       const int nat,
		       double **coul_oo,double **coul_oh,double **coul_hh)
// Calculate derivatives of the force wrt the coulomb parameters.
{
  double Rvec[3];
  double rlen;

  ////main loop Ewald Coulomb:
  const int maxk=10000;
  static int kmax ;
  // const int ksqmax=50;
  const int ksqmax=50;
  static double alphasq=0.7;//0.45;
  // const double alphasq=0.45;
  const double PI=3.14159265359;
  static int totk = 0 ;
  int ksq;
  double Kfac;
  double rkx,rky,rkz;
  double rksq;

  double tempd;
  double alpha=sqrt(alphasq);
  static double r_cut ;
  const double accuracy = EWALD_ACCURACY ;

  // const double r_cut = 4.0 / alpha ;

  double Volume=Latcons[0]*Latcons[1]*Latcons[2];
  static bool called_before = false ;
  double dx,dy,dz,rlen_mi;
  static char *Lbc ;


  double ke=1.0;//332.0637157615209;//this is the unit conversion
  //to achieve charges in nice electron units. 
  //currently we apply this conversion at MD-level, not here.

  static double *Kfac_v, *kx_v, *ky_v, *kz_v ;

  if ( ! called_before ) {
    double r_acc, k_acc ;

    optimal_ewald_params(accuracy, Volume, nat, alpha, 
			 r_cut, kmax,r_acc, k_acc) ;
    alphasq = alpha * alpha ;
    printf("R-Space Ewald cutoff      = %13.6e\n", r_cut) ;
    printf("R-Space accuracy estimate = %13.6e\n", r_acc) ;
    printf("K-space accuracy estimate = %13.6e\n", k_acc) ;

    called_before = true ;

    
    Kfac_v = new double [maxk] ;
    kx_v   = new double [maxk] ;
    ky_v   = new double [maxk] ;
    kz_v   = new double [maxk] ;
    totk = 0 ;
    double min_vfac = 1.0e8 ;
    for(int kx=0;kx<=kmax;kx++) 
      for(int ky=-kmax;ky<=kmax;ky++)
	for(int kz=-kmax;kz<=kmax;kz++)
	  {
	    ksq=kx*kx+ky*ky+kz*kz;
	    if(ksq!=0 and ksq<ksqmax)
	      {
		rkx=(2.0*PI/Latcons[0])*kx;
		kx_v[totk] = rkx ;
		
		rky=(2.0*PI/Latcons[1])*ky;
		ky_v[totk] = rky ;

		rkz=(2.0*PI/Latcons[2])*kz;
		kz_v[totk] = rkz ;

		rksq=rkx*rkx+rky*rky+rkz*rkz;
		// Note:  The original version of the Ewald evaluator did not work 
		// when the lattice constants were different.
		Kfac_v[totk]=exp(-rksq/(4.0*alphasq))/rksq;
		if ( Kfac_v[totk] < min_vfac ) min_vfac = Kfac_v[totk] ;
		totk++ ;
		if(totk>maxk)
		  {
		    cout<<"totk="<<totk<<" greater than maxk="<<maxk<<endl;
		    exit(1);
		  }
	      }
	  }
      cout << "Number of Ewald K-vectors = " << totk << endl ;

      Lbc = new char [nat] ;
  }
  
  // Pack element names into single characters for efficiency
  // Recalculate each step in case of atom reordering between frames.
  for(int a=0;a<nat;a++) {
    if ( Lb[a] == "O" ) {
      Lbc[a] = 'O' ;
    } else if ( Lb[a] == "H" ) {
      Lbc[a] = 'H' ;
    } else {
      cout << "Error: unknown element " << Lb[a] << "\n" ;
    }
  }

  for(int a1=0;a1<nat;a1++)
    for(int a2=0;a2<a1;a2++)
      {
	//Ewald real-space sum.

        for(int c=0;c<3;c++)
          Rvec[c]=Coord[a2][c]-Coord[a1][c];
        rlen=sqrt(Rvec[0]*Rvec[0]+Rvec[1]*Rvec[1]+Rvec[2]*Rvec[2]);

	dx=Rvec[0];
	dy=Rvec[1];
	dz=Rvec[2];

	dx=dx-floor(0.5+dx/Latcons[0])*Latcons[0];//real-space sum done in M.I.C.
	dy=dy-floor(0.5+dy/Latcons[1])*Latcons[1];
	dz=dz-floor(0.5+dz/Latcons[2])*Latcons[2];

	rlen_mi=dx*dx+dy*dy+dz*dz;

	if ( rlen_mi < r_cut * r_cut ) 
	  {
	    rlen_mi = sqrt(rlen_mi) ;
	    //tempd is the derivate of tempy--force.
	    tempd=( (-2.0/sqrt(PI))*exp(-1*alpha*alpha*rlen_mi*rlen_mi)*rlen_mi*alpha - erfc(alpha*rlen_mi))/(rlen_mi*rlen_mi);
	    tempd*=-1.0;

	    tempd*=2.0;  // Sum a1 > a2, not all a2
	    tempd*=ke;
	    tempd *= 0.5 / rlen_mi ;
	    //Note that for least-squares problem, the terms are Force/(Q[a1]*Q[a2]) because Q*Q is the linear multiplier.


	    if(Lbc[a1]=='O' and Lbc[a2]=='O')
	      {
		coul_oo[a1][0]+=tempd*dx;//populate terms according to x,y,z component.
		coul_oo[a1][1]+=tempd*dy;
		coul_oo[a1][2]+=tempd*dz;
			                               
		coul_oo[a2][0]-=tempd*dx;
		coul_oo[a2][1]-=tempd*dy;
		coul_oo[a2][2]-=tempd*dz;
	      }
	    else if((Lbc[a1]=='O' and Lbc[a2]=='H') or (Lbc[a1]=='H' and Lbc[a2]=='O'))
	      {
		coul_oh[a1][0]+=tempd*dx;
		coul_oh[a1][1]+=tempd*dy;
		coul_oh[a1][2]+=tempd*dz;

		coul_oh[a2][0]-=tempd*dx;
		coul_oh[a2][1]-=tempd*dy;
		coul_oh[a2][2]-=tempd*dz;

	      }
	    else if(Lbc[a1]=='H' and Lbc[a2]=='H')
	      {
		coul_hh[a1][0]+=tempd*dx;
		coul_hh[a1][1]+=tempd*dy;
		coul_hh[a1][2]+=tempd*dz;

		coul_hh[a2][0]-=tempd*dx;
		coul_hh[a2][1]-=tempd*dy;
		coul_hh[a2][2]-=tempd*dz;

	      }
	  }
      }
  for(int a1=0;a1<nat;a1++)
    for(int a2=0 ; a2< a1 ;a2++)
      {
	//Ewald K-space sum.

        for(int c=0;c<3;c++)
          Rvec[c]=Coord[a2][c]-Coord[a1][c];

	// Sum over all k vectors.
	for ( int ik = 0 ; ik < totk ; ik++ ) {
	  rkx = kx_v[ik] ;
	  rky = ky_v[ik] ;
	  rkz = kz_v[ik] ;
	  Kfac= Kfac_v[ik];//exp(-1.0*rksq/(4.0*alpha*alpha))/rksq;

	  tempd=sin(rkx*Rvec[0]+rky*Rvec[1]+rkz*Rvec[2]);
	  tempd*=ke;
	  tempd*=2.0 ; // Sum a1 > a2, not all a2.
	  tempd *= Kfac * 0.5*(4*PI/Volume);

	  if ( fabs(rkx) > 1.0e-10 ) tempd*=2.0 ; // Sum kx >= 0, not all kx ;

	  if(Lbc[a1]=='O' and Lbc[a2]=='O')
	    {
	      coul_oo[a1][0]+=tempd*rkx;
	      coul_oo[a1][1]+=tempd*rky;
	      coul_oo[a1][2]+=tempd*rkz;
						                                 
	      coul_oo[a2][0]-=tempd*rkx;
	      coul_oo[a2][1]-=tempd*rky;
	      coul_oo[a2][2]-=tempd*rkz;

	    }
	  else if((Lbc[a1]=='O' and Lbc[a2]=='H') or (Lbc[a1]=='H' and Lbc[a2]=='O'))
	    {
	      coul_oh[a1][0]+=tempd*rkx;
	      coul_oh[a1][1]+=tempd*rky;
	      coul_oh[a1][2]+=tempd*rkz;
		    
	      coul_oh[a2][0]-=tempd*rkx;
	      coul_oh[a2][1]-=tempd*rky;
	      coul_oh[a2][2]-=tempd*rkz;

	    }
	  else if(Lbc[a1]=='H' and Lbc[a2]=='H')
	    {
	      coul_hh[a1][0]+=tempd*rkx;
	      coul_hh[a1][1]+=tempd*rky;
	      coul_hh[a1][2]+=tempd*rkz;
		    
	      coul_hh[a2][0]-=tempd*rkx;
	      coul_hh[a2][1]-=tempd*rky;
	      coul_hh[a2][2]-=tempd*rkz;

	    }
	}
      }

  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      {
	coul_oo[a1][c]*=-1.0;
	coul_oh[a1][c]*=-1.0;
	coul_hh[a1][c]*=-1.0;
      }

  return;
}
