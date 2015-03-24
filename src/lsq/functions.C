#include "functions.h"

static void ZCalc_Spline_Deriv(double **Coord,string *Lb, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double smin,
			       const double smax,
			       const double sdelta,const int snum) ;

static void ZCalc_Ewald_Deriv(double **Coord, string *Lb, 
			      double *Latcons,const int nlayers,
			      const int nat,
			      double **coul_oo,double **coul_oh,double **coul_hh) ;
static void ZCalc_Cheby_Deriv(double **Coord,string *Lb, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double smin,
			       const double smax,
			      const double sdelta,const int snum)       ;

void ZCalc_Deriv(double **Coord,string *Lb, 
		 double *Latcons,const int nlayers,
		 const int nat,double ***A,const double smin,const double smax,
		 const double sdelta,const int snum, double **coul_oo,
		 double **coul_oh,double **coul_hh,bool if_cheby)
{
  bool if_ewald = true ;
  bool if_spline = false ;

  if ( if_cheby == false ) if_spline = true ;
  else if_spline = false ;

  for(int a=0;a<nat;a++)
    for(int n=0;n<snum;n++)
      for(int c=0;c<3;c++)
	  A[a][n][c]=0.0;

  for(int a=0;a<nat;a++)
    {
      for(int c=0;c<3;c++) 
	{
	  coul_oo[a][c]=0.0;
	  coul_oh[a][c]=0.0;
	  coul_hh[a][c]=0.0;
	}
    }

  if ( if_ewald ) {
    ZCalc_Ewald_Deriv(Coord, Lb, Latcons, nlayers, nat, 
		      coul_oo, coul_oh, coul_hh) ;
  }

  if ( if_spline ) {
    ZCalc_Spline_Deriv(Coord,Lb,Latcons,nlayers, nat,A,smin,smax,sdelta,snum) ;
  }

  if ( if_cheby ) {    
    ZCalc_Cheby_Deriv(Coord,Lb,Latcons,nlayers, nat,A,smin,smax,sdelta,snum) ;
  }

}

static void ZCalc_Ewald_Deriv(double **Coord, string *Lb, 
			      double *Latcons,const int nlayers,
			      const int nat,
			      double **coul_oo,double **coul_oh,double **coul_hh)
// Calculate derivatives of the force wrt the coulomb parameters.
{
  double Rvec[3];
  double rlen;

  ////main loop Ewald Coulomb:
  const int maxk=10000;
  const int kmax=10;
  // const int ksqmax=50;
  const int ksqmax=50;
  const double alpha=0.7;//0.45;
  // const double alpha=0.45;
  const double PI=3.14159265359;
  static int totk = 0 ;
  int ksq;
  double Kfac;
  double rkx,rky,rkz;
  double rksq;

  double tempd;
  double kappa=sqrt(alpha);
  const double r_cut = 5.0 / kappa ;

  // const double r_cut = 4.0 / kappa ;

  double Volume=Latcons[0]*Latcons[1]*Latcons[2];
  static bool called_before = false ;
  double dx,dy,dz,rlen_mi;
  static char *Lbc ;


  double ke=1.0;//332.0637157615209;//this is the unit conversion
  //to achieve charges in nice electron units. 
  //currently we apply this conversion at MD-level, not here.

  static double *Kfac_v, *kx_v, *ky_v, *kz_v ;

  if ( ! called_before ) {
    printf("R-Space Ewald cutoff      = %13.6e\n", r_cut) ;
    printf("R-Space accuracy estimate = %13.6e\n", erfc(kappa * r_cut) ) ;
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
		Kfac_v[totk]=exp(-rksq/(4.0*alpha))/rksq;
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
      cout << "K-space accuracy estimate = " << min_vfac << endl ;

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
	    tempd=( (-2.0/sqrt(PI))*exp(-1*kappa*kappa*rlen_mi*rlen_mi)*rlen_mi*kappa - erfc(kappa*rlen_mi))/(rlen_mi*rlen_mi);
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
	  Kfac= Kfac_v[ik];//exp(-1.0*rksq/(4.0*kappa*kappa))/rksq;

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


static void ZCalc_Spline_Deriv(double **Coord,string *Lb, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double smin,
			       const double smax,
			       const double sdelta,const int snum)
// Calculate derivatives of the forces wrt the spline parameters.
{
  double Rvec[3];
  double Rab[3];
  double rlen;
  double t;
  double h00,h10,h01,h11;
  int k0;
  double x,x0;
  int vstart,kstart;

  static bool called_before = false ;
  static char *Lbc ;


  if ( ! called_before ) {
    called_before = true ;
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

  ////main loop for SPLINE terms:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	//calculate vstart: (index for populating OO, OH, or HH column block of A).
	if(Lbc[a1]=='O' and Lbc[a2]=='O')
	  vstart=0;
	if((Lbc[a1]=='O' and Lbc[a2]=='H') or (Lbc[a1]=='H' and Lbc[a2]=='O'))
	  vstart=snum/3;
	if(Lbc[a1]=='H' and Lbc[a2]=='H')
	  vstart=2*snum/3;

	// Start with minimum image convention.  Use layers to access larger distances if desired.	
	for(int c=0;c<3;c++) {
	  Rvec[c]=Coord[a2][c]-Coord[a1][c];
	  Rvec[c] -= floor(0.5+Rvec[c]/Latcons[c])*Latcons[c];
	}
      
	for(int n1=-1*nlayers;n1<nlayers+1;n1++)
	  for(int n2=-1*nlayers;n2<nlayers+1;n2++)
	    for(int n3=-1*nlayers;n3<nlayers+1;n3++)
	      {
		Rab[0]=Rvec[0]+n1*Latcons[0];
		Rab[1]=Rvec[1]+n2*Latcons[1];
		Rab[2]=Rvec[2]+n3*Latcons[2];
		
		rlen=sqrt(Rab[0]*Rab[0]+Rab[1]*Rab[1]+Rab[2]*Rab[2]);
		
		//spline term calculated w/cutoff:
		if(rlen>smin and rlen<smax)
		  {

		    //here is the setup for Hermite cubic quadrature:
		    k0=int(floor((rlen-smin)/sdelta));
		    x=rlen;
		    
		    x0=smin+sdelta*k0;
		    t=(x-x0)/sdelta;

		    //for classical cubic quadrature, you would only
		    //change h00, h10, h01, h11 polynomials!
		    h00=2*t*t*t-3*t*t+1;
		    h10=t*t*t-2*t*t+t;
		    h01=-2*t*t*t+3*t*t;
		    h11=t*t*t-t*t;
		    h10 *= sdelta;//derivative terms have extra factor.
		    h11 *= sdelta;

		    kstart=k0*2;//for each distance unit you have 2 entries.
		    
		    for(int c=0;c<3;c++)
		      {
			A[a1][vstart+kstart+0][c] += h00*Rab[c]/rlen;
			A[a1][vstart+kstart+1][c] += h10*Rab[c]/rlen;
			A[a1][vstart+kstart+2][c] += h01*Rab[c]/rlen;
			A[a1][vstart+kstart+3][c] += h11*Rab[c]/rlen;

			A[a2][vstart+kstart+0][c] -= h00*Rab[c]/rlen;
			A[a2][vstart+kstart+1][c] -= h10*Rab[c]/rlen;
			A[a2][vstart+kstart+2][c] -= h01*Rab[c]/rlen;
			A[a2][vstart+kstart+3][c] -= h11*Rab[c]/rlen;


		      }
	    
		  }//rlen

	      }
      }
  return;
}


static void ZCalc_Cheby_Deriv(double **Coord,string *Lb, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double smin,
			       const double smax,
			       const double sdelta,const int snum)
// Calculate derivatives of the forces wrt the Chebyshev parameters.
{
  double Rvec[3];
  double Rab[3];
  double rlen;
  double x;
  int vstart ;
  static double *Tn, *Tnd ;
  static bool called_before = false ;
  static char *Lbc ;
  const double lambda = 1.25 ;

  if ( ! called_before ) {
    called_before = true ;
    Lbc = new char [nat] ;
    Tn   = new double [snum/3+1] ;
    Tnd  = new double [snum/3+1] ;
    Tn[1] = 1 ;
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

  ////main loop for Chebyshev terms:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	//calculate vstart: (index for populating OO, OH, or HH column block of A).
	if(Lbc[a1]=='O' and Lbc[a2]=='O')
	  vstart=0;
	if((Lbc[a1]=='O' and Lbc[a2]=='H') or (Lbc[a1]=='H' and Lbc[a2]=='O'))
	  vstart=snum/3;
	if(Lbc[a1]=='H' and Lbc[a2]=='H')
	  vstart=2*snum/3;

	// Start with minimum image convention.  Use layers to access larger distances if desired.	
	for(int c=0;c<3;c++) {
	  Rvec[c]=Coord[a2][c]-Coord[a1][c];
	  Rvec[c] -= floor(0.5+Rvec[c]/Latcons[c])*Latcons[c];
	}
      
	for(int n1=-1*nlayers;n1<nlayers+1;n1++)
	  for(int n2=-1*nlayers;n2<nlayers+1;n2++)
	    for(int n3=-1*nlayers;n3<nlayers+1;n3++)
	      {
		Rab[0]=Rvec[0]+n1*Latcons[0];
		Rab[1]=Rvec[1]+n2*Latcons[1];
		Rab[2]=Rvec[2]+n3*Latcons[2];
		
		rlen=sqrt(Rab[0]*Rab[0]+Rab[1]*Rab[1]+Rab[2]*Rab[2]);
		
		//spline term calculated w/cutoff:
		bool inverse_r = false ;
		bool morse     = true ;

		if(rlen>smin and rlen<smax)
		  {
		    double xavg, xdiff ;

		    if ( inverse_r ) {
		      xavg = 0.5 * (1.0/smin + 1.0/smax) ;
		      xdiff = 0.5 * (1.0/smin - 1.0/smax) ;
		      x = (1.0/rlen-xavg) / xdiff ;
		    } else if ( morse ) {
		      double xmin = exp(-smax/lambda) ;
		      double xmax = exp(-smin/lambda) ;
		      xavg = 0.5 * (xmin + xmax) ;
		      xdiff = 0.5 * (xmax - xmin) ;
		      x = (exp(-rlen/lambda)-xavg)/xdiff ;
		    }
		    else {
		      xavg = 0.5 * (smin + smax) ;
		      xdiff = 0.5 * (smax - smin) ;
		      x = (rlen-xavg) / xdiff ;
		    }

		    // x = (rlen-smin) / (smax - smin) ;

		    if ( x < -1.0 ) {
		      cout << "Warning:  r < rmin\n" ;
		      x = -1.0 ;
		    }
		    if ( x > 1.0 ) {
		      x = 1.0 ;
		    }
		    // Generate Chebyshev polynomials by recursion.
		    Tn[0] = 1.0 ;
		    Tn[1] = x ;
		    // Tnd is the derivative of the chebyshev polynomial = n * Chebyshev
		    // polynomial of the second type.  First find Cheby of second type (Un)
		    Tnd[0] = 1.0 ;
		    Tnd[1] = 2.0 * x ;
		    for ( int i = 2 ; i <= snum / 3 ; i++ ) {
		      Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2] ;
		      Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2] ;
		    }
		    // Now store dTn/dx = n * U_(n-1)
		    for ( int i = snum/ 3 ; i >= 1 ; i-- ) {
		      Tnd[i] = i * Tnd[i-1] ;
		    }
		    Tnd[0] = 0.0 ;

		    double rlen3 = rlen * rlen * rlen ;
		    if ( morse ) 
		      {
			double fcut0 = (1.0 - rlen/smax) ;
			double fcut = fcut0 * fcut0 * fcut0 ;
			double fcutderiv = -3.0 * fcut0 * fcut0 / smax ;
			for ( int i = 0 ; i < snum / 3 ; i++ ) 
			  {
			    for(int c=0;c<3;c++)
			      {
				double deriv = 
				  (fcut * Tnd[i+1] *(-exp(-rlen/lambda)/lambda)/xdiff + 
				   fcutderiv * Tn[i+1]
				   ) * Rab[c] / rlen ; 
				A[a1][vstart+i][c] += deriv ;
				A[a2][vstart+i][c] -= deriv ;
			      } 
			  }
		      }
		    else 
		      {
			for ( int i = 0 ; i < snum / 3 ; i++ ) 
			  {
			    for(int c=0;c<3;c++) 
			      {
				A[a1][vstart+i][c] += (smax - rlen) * Tn[i]*Rab[c]/rlen3;
				A[a2][vstart+i][c] -= (smax - rlen) * Tn[i]*Rab[c]/rlen3;
			      }
			  }
		      }
		  }//rlen
	      }
      }
  return;
}

//this function subtracts the ReaxFF over-coordination term to re-fit 
//splines/charges iteratively for self-consistence.
// If calc_deriv is true, the derivative of the force wrt the magnitude of the
// 3-body interaction is placed in Fderiv.  Otherwise, the 3-body force is subtracted
// from the total forces given in force.  
void SubtractCoordForces(double **Coord,double **Force,string *Lb, double *Latcons,
			 const int nlayers, const int nat, bool calc_deriv, 
			 double **Fderiv)
{

  double Rvec[3];
  double Rab[3];
  double Vtot;

  double SForce[nat][3];
  for(int a=0;a<nat;a++)
    for(int c=0;c<3;c++)
      SForce[a][c]=0;

  double p[35];
  for(int n=0;n<35;n++)
    p[n]=0.0;

  ////THREE-BODY POTENTIAL:
  double Eover;
  double p1,p2,r0,pover,lambda6;
  double rik,tempr,temps;


  //really good here:
  // p1=  -2.5881042987450e-01   ; 
  // r0=  9.6000387075695e-01    ;
  // p2=  3.8995379237250e+00    ;
  // lambda6=     -8.9                 ;  
 
  // Chenoweth values for p1, p2, r0, lambda6 p(ovun2)
  pover=50.0;
  p1 = -0.0657 ;
  p2 = 5.0451 ;
  r0 = 1.0165 ;
  lambda6 = -3.6141 ;

  Vtot=0.0;
  Eover=0.0;

  double S[nat];
  for(int ai=0;ai<nat;ai++)
    S[ai]=0.0;

  for(int ai=0;ai<nat;ai++)
    if(Lb[ai]=="O")
      {
	temps=0.0;
	for(int ak=0;ak<nat;ak++)
	  {
	    // Minimum image only for overcoordination.
	    for(int c=0;c<3;c++) {
	      Rvec[c]=Coord[ak][c]-Coord[ai][c];
	      Rvec[c] -= floor(0.5+Rvec[c]/Latcons[c])*Latcons[c];
	    }

	    Rab[0]=Rvec[0];
	    Rab[1]=Rvec[1];
	    Rab[2]=Rvec[2];

	    rik=sqrt(Rab[0]*Rab[0]+Rab[1]*Rab[1]+Rab[2]*Rab[2]);
		    		    
	    if(ai==ak)
	      rik=1000.0;
	    
	    temps+=exp(p1*pow(rik/r0,p2));
	  }
	S[ai]=temps-2.0;
	Eover+=pover*S[ai]*1.0/(1.0+exp(lambda6*S[ai]));
      }

  Vtot+=Eover;
  

  double dEover[nat][3];
  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      dEover[a1][c]=0.0;
  
  for(int aj=0;aj<nat;aj++) 
    {
      for(int ai=0;ai<nat;ai++) 
	{
      if(Lb[ai]=="O")
	for(int ak=0;ak<nat;ak++)
	  if(aj==ai or aj==ak)
	    {
	      for(int c=0;c<3;c++) {
		  Rvec[c]=Coord[ak][c]-Coord[ai][c];
		  Rvec[c] -= floor(0.5+Rvec[c]/Latcons[c])*Latcons[c];
	      }		
	      Rab[0]=Rvec[0];
	      Rab[1]=Rvec[1];
	      Rab[2]=Rvec[2];

	      rik=sqrt(Rab[0]*Rab[0]+Rab[1]*Rab[1]+Rab[2]*Rab[2]);


	      tempr=pover;
	      tempr*=1.0/(1+exp(lambda6*S[ai])) - lambda6*S[ai]*exp(lambda6*S[ai])/pow(1.0+exp(lambda6*S[ai]),2);
	      tempr*=p1*pow(rik/r0,p2)*p2*exp(p1*pow(rik/r0,p2))/rik;
	      
	      if(ai==ak) 
		{
		  tempr=0.0;
		  rik=1000.0;
		}

	      for(int c=0;c<3;c++)
		{
		  if(aj==ai)
		    dEover[aj][c]-=tempr*Rab[c]/rik;
		  if(aj==ak)
		    dEover[aj][c]+=tempr*Rab[c]/rik;
		}
	    }
	}
  }

  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      SForce[a1][c] -= dEover[a1][c];


  if ( calc_deriv == false ) {
    //subtract Coord. force:
    for(int a1=0;a1<nat;a1++)
      for(int c=0;c<3;c++)
	Force[a1][c] -= SForce[a1][c];
  }
  else 
    {
      // Calculate derivative of 3-body force wrt pover.
      for(int a1=0;a1<nat;a1++)
	for(int c=0;c<3;c++)
	  Fderiv[a1][c] = -dEover[a1][c]/pover ;
    }
}




