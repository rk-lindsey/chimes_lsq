#include<iostream>
#include<fstream>
#include "functions.h"

static void ZCalc_Lj(double **Coord,double *Latcons,
		     const int nat,double **SForce,double& Vtot,double& Pxyz) ;
static void ZCalc_SR_Analytic(double **Coord,const char *Lbc, double *Latcons,const int nlayers,
			      const int nat,const double *smin,const double *smax, const int *snum,
			      double **SForce,double& Vtot,double& Pxyz, double *params) ;
static void ZCalcSR_Over(double **Coord,const char *Lbc, double *Latcons,
			 const int nat,double **SForce,double& Vtot,double& Pxyz, int n_over,
			 double *over_param) ;
static void ZCalc_Stillinger(double **Coord, const char *Lbc, double *Latcons,const int nlayers,
			     const int nat, const double *smax, 
			     double **SForce,double& Vtot,double& Pxyz) ;
static void ZCalc_Spline(double **Coord,const char *Lbc, double *Latcons,const int nlayers,
			 const int nat,const double *smin,const double *smax,
			 const double *sdelta,const int *snum, 
			 double *params, double *pot_params,
			 double **SForce,double& Vtot,double& Pxyz) ;

static void ZCalc_Cheby(double **Coord,const char *Lbc, double *Latcons,
			const int nlayers,
			const int nat,const double *smin,
			const double *smax,
			const int *snum, 
			double *params, const double *lambda,
			double **SForce, double &Vtot, double &xyz);


static void ZCalc_Spline_Deriv(double **Coord, const char *Lbc, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double *smin,
			       const double *smax,
			       const double *sdelta,const int *snum,double *mind) ;

static void ZCalc_InvR_Deriv(double **Coord,const char *Lbc, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double *smin,
			       const double *smax,
			      const int *snum, double *mind) ;

static void ZCalc_Poly_Deriv(double **Coord,const char *Lbc, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double *smin,
			       const double *smax,
			      const int *snum, double *mind) ;

static void ZCalc_Cheby_Deriv(double **Coord,const char *Lbc, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double *smin,
			       const double *smax,
			      const int *snum, const double *lambda, double *mind) ;


static int pair_table_offset(int ipair, const int *snum) ;

void ZCalc(double **Coord, const char *Lbc, double *Q, double *Latcons,
	   const int nlayers,
	   const int nat,const double *smin,const double *smax,
	   const double *sdelta,const int *snum, 
	   const int *snum_3b_cheby,
	   double *params, double *pot_params, Sr_pair_t pair_type,
	   bool if_coulomb, bool if_overcoord, 
	   bool if_3b_cheby,
	   int n_over,
	   double *over_params, const double *lambda,
	   double **SForce,double& Vtot,double& Pxyz)
// Calculate the force, potential energy, and pressure.
{
  double volume ;

  for(int a=0;a<nat;a++)
    for(int c=0;c<3;c++)
      SForce[a][c]=0;

  Vtot = 0.0;
  Pxyz = 0.0 ;

  if ( pair_type == LJ ) 
    {
      ZCalc_Lj(Coord, Latcons, nat, SForce, Vtot, Pxyz) ;
    } 
  else if ( pair_type == INVERSE_R ) 
    {
      ZCalc_SR_Analytic(Coord,Lbc, Latcons,nlayers,nat,smin,smax,snum, SForce,Vtot,
			Pxyz, params) ;
    } 
  else if ( pair_type == STILLINGER ) 
    {
      ZCalc_Stillinger(Coord,Lbc, Latcons,nlayers,nat,smax, SForce,Vtot,Pxyz) ;
    } 
  else if ( pair_type == CHEBYSHEV ) 
    {
      ZCalc_Cheby(Coord,Lbc,Latcons,nlayers,nat,smin,smax,snum, 
		  params,lambda,SForce,Vtot,Pxyz) ;
    } 
  else if ( pair_type == SPLINE ) 
    {
      ZCalc_Spline(Coord,Lbc,Latcons,nlayers,nat,smin,smax,sdelta,snum, 
		   params,pot_params,SForce,Vtot,Pxyz) ;
    } 
  else 
    {
      cout << "Error: unknown pair type\n" ;
      exit(1) ;
    }

  if ( if_coulomb ) 
    {
      ZCalc_Ewald(Coord, Lbc, Q, Latcons, nat, SForce, Vtot, Pxyz) ;
    }

  if ( if_overcoord ) 
    {
      // Take the overall magnitude of the overcoordination term from the list
      // of linear parameters in params.
      // cout << "POVER = " << pover << endl ;
      ZCalcSR_Over(Coord,Lbc, Latcons,nat,SForce, Vtot, Pxyz,
		   n_over, over_params) ;
    }

  if ( if_3b_cheby ) 
    {
      // 3-body chebyshev polynomial 
      ZCalc_3B_Cheby(Coord, Lbc, Latcons, nat, smin, smax, snum,
		     snum_3b_cheby,
		     params, lambda, SForce, Vtot, Pxyz) ;
    }

  volume = Latcons[0] * Latcons[1] * Latcons[2] ;
  Pxyz /= 3.0 * volume ;

  return;
}


static void ZCalc_Lj(double **Coord, double *Latcons,
		     const int nat,double **SForce,double& Vtot,double& Pxyz)
// Calculate LJ interaction
{
  double Rvec[3];

  const double eps =    1.0 ;
  const double sigma = 1.1 ;

  for(int a1=0;a1<nat;a1++) 
    {
      for(int a2=0;a2<a1;a2++)
	{

	  for(int c=0;c<3;c++) 
	    Rvec[c]=Coord[a2][c]-Coord[a1][c];
	
	  double dx=Rvec[0];
	  double dy=Rvec[1];
	  double dz=Rvec[2];
	    
	  dx=dx-floor(0.5+dx/Latcons[0])*Latcons[0];
	  dy=dy-floor(0.5+dy/Latcons[1])*Latcons[1];
	  dz=dz-floor(0.5+dz/Latcons[2])*Latcons[2];
	  
	  double rlen_mi=sqrt(dx*dx+dy*dy+dz*dz);
	  if ( rlen_mi < 0.5 ) {
	    cout << "Error: close approach\n" ;
	    exit(1) ;
	  }

	  Vtot += 4.0 * eps * ( pow(sigma/rlen_mi,12.0) - pow(sigma/rlen_mi,6.0) ) ;

	  double fac = 4.0 * eps * ( 
				    -12.0 * pow(sigma/rlen_mi,14.0) 
				    + 6.0 * pow(sigma/rlen_mi, 8.0) ) ;
	  fac *= 1.0 / ( sigma * sigma ) ;
	  
	  Pxyz -= fac * (rlen_mi*rlen_mi) ;

	  SForce[a1][0]+=dx*fac ;
	  SForce[a1][1]+=dy*fac ;
	  SForce[a1][2]+=dz*fac ;
			
	  SForce[a2][0]-=dx * fac ;
	  SForce[a2][1]-=dy * fac ;
	  SForce[a2][2]-=dz * fac ;
	}
    }
}


static void ZCalc_SR_Analytic(double **Coord,const char *Lbc, double *Latcons,const int nlayers,
			      const int nat,const double *smin,const double *smax, const int *snum,
			      double **SForce,double& Vtot,double& Pxyz, double *params)
// Short-range analytic forces - Optimized evaluation.
{
  double Rvec[3];
  double Rab[3] ;
  double rlen,rlen2;
  double S_r;
  double tempx ;
  //these are short-ranged cutoffs for analytical potential.
  //e.g. for r<oo_cut, force(r)=force(oo_cut). The values are 
  //at least 0.2 Angs below any sampled by MD; they are there 
  //in case trajectories reach them.

  double oo_cut=1.8;
  double oh_cut=0.6;
  double hh_cut=0.8;

  double pOO[24] ;
  double pOH[24] ;
  double pHH[24] ;
  for (int i = 0; i < NPAIR; i++){
    for (int j = 0; j < snum[i]; j++) { 
      if (i == 0) {
        pOO[j] = -params[j];
      } else if (i == 1) {
        pOH[j] = -params[snum[i-1]+j];
      } else if (i == 2) {
        pHH[j] = -params[snum[i-1] + snum[i-2] +j];
      }
    }
  }

  ////main loop non-Coulomb short-ranged forces:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	int ipair = pair_index(a1, a2, Lbc) ;

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
		
		    rlen=(Rab[0]*Rab[0]+Rab[1]*Rab[1]+Rab[2]*Rab[2]);

		    //spline term calculated w/cutoff:

		    if(rlen > smin[ipair] * smin[ipair] and rlen < smax[ipair] * smax[ipair] )
		      {
			rlen = sqrt(rlen) ;

			S_r=0.0;

			if(Lbc[a1]=='O' and Lbc[a2]=='O')
			  {
			    rlen2=rlen;
			    if(rlen2<oo_cut)
			      {
				rlen2=oo_cut;
			      }
			  }
			else if((Lbc[a1]=='O' and Lbc[a2]=='H') or (Lbc[a1]=='H' and Lbc[a2]=='O'))
			  {
			    rlen2=rlen;
			    if(rlen2<oh_cut)
			      {
				rlen2=oh_cut;
			      }
			  }
			else if(Lbc[a1]=='H' and Lbc[a2]=='H')
			  {
			    rlen2=rlen;
			    if(rlen2<hh_cut)
			      {
				rlen2=hh_cut;
			      }
			  }
			else
			  {
			    cout << "Bad atom labels found" << endl ;
			    exit(1) ;
			  }
		    
			double smaxpow3 = pow(rlen2-smax[ipair],3) ;
			double smaxpow2 = pow(rlen2-smax[ipair],2) ;
			double rleninv  = 1.0 / rlen2 ;

			//ANALYTICAL 
			//this sum over inverse-powers function for the forces/potential energy is currently
			//commented out. It has the form of Koziol and Fried paper. 
			if(Lbc[a1]=='O' and Lbc[a2]=='O')
			  {

			    double rlenpow2 = rleninv * rleninv ;
			    double rlenpow3 = rlenpow2 * rleninv ;
			    tempx=0.0;
			    for(int n=0;n<24;n++)
			      {
				S_r+=pOO[n]*((-1.0*(n+2))*rlenpow3)*smaxpow3 ;
				S_r+=pOO[n]*(3.0*rlenpow2)*smaxpow2 ;
				tempx+=pOO[n]*rlenpow2 ;
				rlenpow2 *= rleninv ;
				rlenpow3 *= rleninv ;
			      }
			    Vtot += tempx * smaxpow3 ;
			  }
		                                      
			if((Lbc[a1]=='O' and Lbc[a2]=='H') or (Lbc[a1]=='H' and Lbc[a2]=='O'))
			  {
			    double rlenpow2 = rleninv * rleninv ;
			    double rlenpow3 = rlenpow2 * rleninv ;
			    tempx=0.0;
			    for(int n=0;n<24;n++)
			      {
				S_r+=pOH[n]*((-1.0*(n+2))*rlenpow3)*smaxpow3 ;
				S_r+=pOH[n]*(3.0*rlenpow2)*smaxpow2 ;
				tempx += pOH[n]*rlenpow2;

				rlenpow2 *= rleninv ;
				rlenpow3 *= rleninv ;
			      }
                              
			    Vtot += tempx * smaxpow3 ;
			  }
                                  
			if(Lbc[a1]=='H' and Lbc[a2]=='H')
			  {
			    double rlenpow2 = rleninv * rleninv ;
			    double rlenpow3 = rlenpow2 * rleninv ;
			    tempx=0.0;
			    for(int n=0;n<24;n++)
			      {
				S_r+=pHH[n]*((-1.0*(n+2))*rlenpow3)*smaxpow3 ;
				S_r+=pHH[n]*(3.0*rlenpow2)*smaxpow2 ;
				tempx+=pHH[n]*rlenpow2;

				rlenpow2 *= rleninv ;
				rlenpow3 *= rleninv ;
			      }
                              
			    Vtot -= tempx * smaxpow3 ;
			  }
		   
			Pxyz -= S_r * rlen ;

			//this last part just assigns the central-forces to the x,y,z components.
			for(int c=0;c<3;c++)
			  {
			    SForce[a1][c] += S_r*Rab[c]/rlen;
			    SForce[a2][c] -= S_r*Rab[c]/rlen;
			  }

		      }//rlen


		  }

      }
  
  return;
}



static void ZCalcSR_Over(double **Coord, const char *Lbc, double *Latcons,
			 const int nat,double **SForce,double& Vtot,double& Pxyz, 
			 int n_over, double *over_param) 
// Calculate short-range overcoordination forces.
{
  double Rvec[3] ;
  
  //these are short-ranged cutoffs for analytical potential.
  //e.g. for r<oo_cut, force(r)=force(oo_cut). The values are 
  //at least 0.2 Angs below any sampled by MD; they are there 
  //in case trajectories reach them.

  ////Many-body overcoordination term:
  ///The rest of this function can be commented out and the MD 
  //will still run. The ReaxFF terms have been fitted using
  //Powell minimization (separate routines, but very easy to do).
  // 

  ////MANY-BODY POTENTIAL:
  double Eover;
  double pover,p1,p2,r0,lambda6;
  double rik,tempr,temps;

  //ReaxFF over-coordination parameters (from non-linear fitting).
  // pover=75.0;
  // p1=  -2.5881042987450e-01   ; 
  // r0=  9.6000387075695e-01    ;
  // p2=  3.8995379237250e+00    ;
  // lambda6=-8.9;

  // Chenoweth values for p1, p2, r0, lambda6 p(ovun2)
  if ( n_over != 5 ) {
    cout << "Error: wrong number of overcoordination parameters\n" ;
    exit(1) ;
  }

  pover   = over_param[0] ;
  r0      = over_param[1] ;
  p1      = over_param[2] ;
  p2      = over_param[3] ;
  lambda6 = over_param[4] ;

  Eover=0.0;

  double S[nat];
  for(int ai=0;ai<nat;ai++)
    S[ai]=0.0;

  double Sexp[nat] ;
  for(int ai=0;ai<nat;ai++)
    Sexp[ai]=0.0;

  for(int ai=0;ai<nat;ai++)
    if(Lbc[ai]=='O')
      {
	temps=0.0;
	for(int ak=0;ak<nat;ak++)
	  {
	    if(ai==ak)
	      continue ;

	    for(int c=0;c<3;c++)
	      Rvec[c]=Coord[ak][c]-Coord[ai][c];

	    // Short-range interaction, so use minimum image convention.
	    double Rab[3] ;

	    Rab[0] = Rvec[0] - floor(0.5+Rvec[0]/Latcons[0])*Latcons[0];
	    Rab[1] = Rvec[1] - floor(0.5+Rvec[1]/Latcons[1])*Latcons[1];
	    Rab[2] = Rvec[2] - floor(0.5+Rvec[2]/Latcons[2])*Latcons[2];

	    rik=sqrt(Rab[0]*Rab[0]+Rab[1]*Rab[1]+Rab[2]*Rab[2]);

	    temps+=exp(p1*pow(rik/r0,p2));
	  }
	S[ai]=temps-2.0;
	// printf("Overcoordination of %d = %13.6e\n", ai, S[ai]) ;
	Sexp[ai] = exp(lambda6*S[ai]) ;
	Eover+=pover*S[ai]*1.0/(1.0+Sexp[ai]);
      }

  Vtot+=Eover;
  

  double dEover[nat][3];
  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      dEover[a1][c]=0.0;
  
  for(int aj=0;aj<nat;aj++)
    for(int ai=0;ai<nat;ai++)
      if(Lbc[ai]=='O')
	for(int ak=0;ak<nat;ak++)
	  if(aj==ai or aj==ak)
	    {
	      if(ai==ak ) 
		{
		  continue ;
		}

	      for(int c=0;c<3;c++)
		Rvec[c]=Coord[ak][c]-Coord[ai][c];

		
	      // Short-range interaction, so use minimum image convention.
	      double Rab[3] ;
	      Rab[0] = Rvec[0] - floor(0.5+Rvec[0]/Latcons[0])*Latcons[0];
	      Rab[1] = Rvec[1] - floor(0.5+Rvec[1]/Latcons[1])*Latcons[1];
	      Rab[2] = Rvec[2] - floor(0.5+Rvec[2]/Latcons[2])*Latcons[2];

	      rik=sqrt(Rab[0]*Rab[0]+Rab[1]*Rab[1]+Rab[2]*Rab[2]);

	      double powrik = pow(rik/r0,p2) ;
	      double Sexp2  = (1.0 + Sexp[ai]) * (1.0 + Sexp[ai]) ;

	      tempr=pover;
	      tempr*=1.0/(1+Sexp[ai]) - lambda6*S[ai]*Sexp[ai]/Sexp2 ;
	      tempr*=p1*powrik*p2*exp(p1*powrik)/rik;
	      
	      if(aj==ai) 
		{
		  Pxyz -= 0.5 * tempr * rik ;
		  for(int c=0;c<3;c++)
		    {
		      dEover[aj][c]-=tempr*Rab[c]/rik;
		    }
		}

	      if(aj==ak) {
		Pxyz -= 0.5 * tempr * rik ;
		for(int c=0;c<3;c++)
		  {
		    dEover[aj][c]+=tempr*Rab[c]/rik;
		  }
	      }
	    }
  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      SForce[a1][c] -= dEover[a1][c];
  
    return;
}

static void ZCalc_Stillinger(double **Coord,const char *Lbc, double *Latcons,const int nlayers,
			     const int nat, const double *smax, 
			     double **SForce,double& Vtot,double& Pxyz)
// Calculate stillinger model forces (no charges).
{
  double Rvec[3], Rab[3] ;
  double rlen ;
  double S_r;
  bool ifcalcsr;

  //This defines original Stillinger parameters for running Stillinger-MD.
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

  ////spline charges Stillinger units
 
  ////main loop non-Coulomb short-ranged forces:
  for(int a1=0;a1<nat-1;a1++) 
    for(int a2=a1+1;a2<nat;a2++)
      {
	int ipair = pair_index(a1,a2,Lbc) ;

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

		// Stillinger potential calculated with cutoff.

		if( rlen < smax[ipair] )
		  {
		    ifcalcsr=true;

		    S_r=0.0;

		    ////STILLINGER
		    
		    if(Lbc[a1]=='O' and Lbc[a2]=='O')
		      {
			S_r += p_oo[1]*p_oo[2]*(-1)/pow(rlen,p_oo[2]+1);
			////Energy:
			Vtot += p_oo[1]/pow(rlen,p_oo[2]);
		      }
		    
		    else if((Lbc[a1]=='O' and Lbc[a2]=='H') or (Lbc[a1]=='H' and Lbc[a2]=='O'))
		      {
			//Second:
			S_r += (-1)*p_oh[1]*p_oh[2]/pow(rlen,p_oh[2]+1);
			//Third:
			S_r += (-1)*p_oh[3]*(-1)*pow(1+exp(p_oh[4]*(rlen-p_oh[5])),-2)*exp(p_oh[4]*(rlen-p_oh[5]))*p_oh[4];
			
			//Energy:
			Vtot +=  p_oh[1]/pow(rlen,p_oh[2]) - p_oh[3]/(1+exp(p_oh[4]*(rlen-p_oh[5])));	
		      }
		    else if(Lbc[a1]=='H' and Lbc[a2]=='H')
		      {
			//Second:
			S_r += p_hh[1]*(-1)*p_hh[2]/(exp(-1*p_hh[2]*(rlen-p_hh[3]))+exp(p_hh[2]*(rlen-p_hh[3]))+2);
			//Third:
			S_r += (-1)*p_hh[4]*exp(-1*p_hh[5]*(rlen-p_hh[6])*(rlen-p_hh[6]))*(-1*p_hh[5])*2*(rlen-p_hh[6]);
			//Energy
			Vtot += p_hh[1]/(1+exp(p_hh[2]*(rlen-p_hh[3]))) - p_hh[4]*exp(-1*p_hh[5]*(rlen-p_hh[6])*(rlen-p_hh[6]));
		      }
		  

		    //this last part just assigns the central-forces to the x,y,z components.
		    Pxyz -= S_r * rlen ;
		    for(int c=0;c<3;c++)
		      {
			SForce[a1][c] += S_r*Rab[c]/rlen;
			SForce[a2][c] -= S_r*Rab[c]/rlen;
		      }
		  }//rlen
	      }
      }
  return;
}


static void ZCalc_Spline(double **Coord,const char *Lbc, double *Latcons,const int nlayers,
			 const int nat,const double *smin,const double *smax,
			 const double *sdelta,const int *snum, 
			 double *params, double *pot_params, 
			 double **SForce,double& Vtot,double& Pxyz)
// Calculate spline forces.
{
  double Rvec[3], Rab[3] ;
  double rlen,rlen2;
  double tempx;
  double S_r;
  int vstart;

  ////main loop non-Coulomb short-ranged forces:
  for(int a1=0;a1<nat-1;a1++) {
    for(int a2=a1+1;a2<nat;a2++)
      {
	//calculate vstart:
	int ipair = pair_index(a1,a2,Lbc) ;
	vstart = pair_table_offset(ipair, snum) ;
	
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

		if(rlen>smin[ipair] and rlen<smax[ipair])
		  {
		    S_r=0.0;

		    rlen2 = rlen ;
		    tempx = spline_pot(smin[ipair], smax[ipair], sdelta[ipair], 
				       rlen2, params, pot_params, snum[ipair], vstart, S_r) ;

		    //this last part just assigns the central-forces to the x,y,z components.
		    for(int c=0;c<3;c++)
		      {
			SForce[a1][c] += S_r*Rab[c]/rlen;
			SForce[a2][c] -= S_r*Rab[c]/rlen;
		      }

		    Pxyz -= S_r * rlen ;
		    Vtot += tempx ;

		  }//rlen
	      }
      }
  }

  return;
}

double spline_pot(double smin, double smax, double sdelta, double rlen2, double *params, double *pot_params, 
		  int snum, int vstart, double &S_r)
// Return the spline potential at the position rlen2.
// S_r is the derivative of the spline potential with respect to position.
{
  int k0 ;
  double x, x0 ;
  double t;
  double h00,h10,h01,h11;
  double i00,i10,i01,i11;
  int kstart ;
  double tempx ;

  if ( rlen2 > smax ) 
    {
      // Cut off potential.
      S_r = 0.0 ;
      return(0.0) ;
    }
  else if ( rlen2 < smin ) 
    {
      rlen2 = smin ;
    }

  k0=int(floor((rlen2-smin)/sdelta));
  x=rlen2;
  x0=smin+sdelta*k0;
  t=(x-x0)/sdelta;
  if ( t > 1.000001 || t < -0.000001 ) {
    cout << "Error: bad t\n" ;
    exit(1) ;
  }

  double t2, t3, t4 ;

  

  t2 = t * t ;
  t3 = t2 * t ;
  t4 = t3 * t ;

  h00 = 2*t3 -3*t2 + 1.0 ;
  h10 =  t3  -2*t2 + t;
  h01= -2*t3 +3*t2 ;
  h11=   t3 - t2 ;
  h10 *= sdelta;//derivative terms have extra factor.
  h11 *= sdelta;
  kstart=k0*2;		    
  if ( kstart > snum ) {
    cout << "Error: kstart too large " << kstart << " " << rlen2 << " " << smin << " " << snum << endl ;
    exit(1) ;
  }
  S_r=
    h00*params[vstart+kstart+0]+//lookout!!!
    h10*params[vstart+kstart+1]+
    h01*params[vstart+kstart+2]+
    h11*params[vstart+kstart+3];
  
  //do integrating here to get potential energy:
  i00=0.5 - (0.5*t4 - t3 + t) ;
  i10=1.0/12.0 - (0.25*t4- (2.0/3.0)*t3 + 0.5*t2) ;
  i01=0.5 - (-0.5*t4 + t3) ;
  i11=-1.0/12.0 - (0.25*t4 - t3/3.0) ;
  i10*=sdelta;
  i11*=sdelta;

  // TEST: Linear interpolation.
  // i00 = (0.5 - (t - t * t / 2.0)) ;
  // i01 = (0.5 - t * t / 2.0) ;
  // i10 = 0.0 ;
  // i11 = 0.0 ;

  tempx=  -sdelta * (
		     i00*params[vstart+kstart+0]+
		     i10*params[vstart+kstart+1]+
		     i01*params[vstart+kstart+2]+
		     i11*params[vstart+kstart+3] ) ;
  if ( kstart/2 + 1 < snum / 2 ) {
    tempx += pot_params[vstart/2+kstart/2+1] ;
  }
  return(tempx) ;
}


static void ZCalc_Cheby(double **Coord,const char *Lbc, double *Latcons,
			const int nlayers,
			const int nat,const double *smin,
			const double *smax,
			const int *snum, 
			double *params, const double *lambda,
			double **SForce, double &Vtot, double &Pxyz)
// Calculate short-range forces using a Chebyshev polynomial expansion.
// Can use morse variables similar to the work of Bowman.
{
  double Rvec[3];
  double Rab[3];
  double rlen;
  double x;
  int vstart ;
  static double *Tn, *Tnd ;
  static bool called_before = false ;
  //spline term calculated w/cutoff:
  const bool inverse_r = false ;
  const bool morse     = true ;
  double tempx;
  //const double lambda = 1.25 ;
  double exprlen ;

  if ( ! called_before ) {
    called_before = true ;
    int dim = 0 ;
    for ( int i = 0 ; i < NPAIR ; i++ ) 
      {
	if ( snum[i] > dim ) 
	  {
	    dim = snum[i] ;
	  }
      }
    dim++ ;
    Tn   = new double [dim] ;
    Tnd  = new double [dim] ;
  }
  
  // Pack element names into single characters for efficiency
  // Recalculate each step in case of atom reordering between frames.

  ////main loop for Chebyshev terms:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	//calculate vstart: (index for populating OO, OH, or HH column block of A).
	int ipair = pair_index(a1,a2,Lbc) ;
	vstart = pair_table_offset(ipair, snum) ;

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
		
		if( rlen < smax[ipair])
		  {
		    double xavg, xdiff ;

		    if ( inverse_r ) {
		      xavg = 0.5 * (1.0/smin[ipair] + 1.0/smax[ipair]) ;
		      xdiff = 0.5 * (1.0/smin[ipair] - 1.0/smax[ipair]) ;
		      x = (1.0/rlen-xavg) / xdiff ;
		    } else if ( morse ) {
		      double xmin = exp(-smax[ipair]/lambda[ipair]) ;
		      double xmax = exp(-smin[ipair]/lambda[ipair]) ;
		      xavg = 0.5 * (xmin + xmax) ;
		      xdiff = 0.5 * (xmax - xmin) ;
		      exprlen = exp(-rlen/lambda[ipair]) ;
		      x = (exprlen-xavg)/xdiff ;
		    }
		    else {
		      xavg = 0.5 * (smin[ipair] + smax[ipair]) ;
		      xdiff = 0.5 * (smax[ipair] - smin[ipair]) ;
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
		    for ( int i = 2 ; i <= snum[ipair] ; i++ ) {
		      Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2] ;
		      Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2] ;
		    }
		    // Now store dTn/dx = n * U_(n-1)
		    for ( int i = snum[ipair] ; i >= 1 ; i-- ) {
		      Tnd[i] = i * Tnd[i-1] ;
		    }
		    Tnd[0] = 0.0 ;

		    if ( morse ) 
		      {
			double fcut0 = (1.0 - rlen/smax[ipair]) ;
			double fcut = fcut0 * fcut0 * fcut0 ;
			double fcutderiv = -3.0 * fcut0 * fcut0 / smax[ipair] ;
			tempx = 0.0 ;
			for ( int i = 0 ; i < snum[ipair] ; i++ ) 
			  {
			    double coeff = params[vstart+i] ;
			    tempx += coeff * fcut * Tn[i+1] ;
			    double deriv = 
			      (fcut * Tnd[i+1] *(-exprlen/lambda[ipair])/xdiff + 
			       fcutderiv * Tn[i+1]
			       ) ;
			    Pxyz -= coeff * deriv * rlen ;
			    for(int c=0;c<3;c++)
			      {
				SForce[a1][c] += coeff * deriv * Rab[c] / rlen ;
				SForce[a2][c] -= coeff * deriv * Rab[c] / rlen ;
			      } 
			  }
		      }
		    else 
		      {
			double rlen3 = rlen * rlen * rlen ;
			for ( int i = 0 ; i < snum[ipair] ; i++ ) 
			  {
			    double coeff = params[vstart+i] ;
			    tempx = 0.0 ;
			    // POTENTIAL NOT YET IMPLEMENTED.
			    Pxyz += coeff * (smax[ipair]-rlen) * Tn[i] / rlen ;
			    {
			      for(int c=0;c<3;c++) 
				{
				  SForce[a1][c] += coeff * (smax[ipair] - rlen) * Tn[i]*Rab[c]/rlen3;
				  SForce[a2][c] -= coeff * (smax[ipair] - rlen) * Tn[i]*Rab[c]/rlen3;
				}
			    }
			  }
		      }
		    Vtot += tempx ;
		  } //rlen
	      }
      }
  return;
}




void ZCalc_Deriv(double **Coord,const char *Lbc,
		 double *Latcons,const int nlayers,
		 const int nat,double ***A,const double *smin,const double *smax,
		 const double *sdelta,
		 const int *snum, 
		 const int *snum_3b_cheby,
		 const double *lambda,
		 double **coul_oo, double **coul_oh,double **coul_hh, Sr_pair_t pair_type,
		 double *mind, bool if_3b_cheby)
{
  bool if_ewald ;
  if (pair_type != DFTBPOLY) {
    if_ewald = true;
  } else {
    if_ewald = false;
  }
  int snum_tot ;
  
  snum_tot = 0 ;
  for (int i = 0 ; i < NPAIR ; i++ ) {
    snum_tot += snum[i] ;
  }
  if ( if_3b_cheby ) {
    snum_tot += count_cheby_3b_params(snum_3b_cheby) ;
  }

  for(int a=0;a<nat;a++)
    for(int n=0;n<snum_tot;n++)
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
    ZCalc_Ewald_Deriv(Coord, Lbc, Latcons, nat, 
		      coul_oo, coul_oh, coul_hh) ;
  }

  if ( pair_type == SPLINE ) {
    ZCalc_Spline_Deriv(Coord,Lbc,Latcons,nlayers, nat,A,smin,smax,sdelta,snum,mind) ;
  } else if ( pair_type == CHEBYSHEV ) {
    ZCalc_Cheby_Deriv(Coord,Lbc,Latcons,nlayers, nat,A,smin,smax,snum,lambda,mind) ;
  } else if ( pair_type == DFTBPOLY ) {
    ZCalc_Poly_Deriv(Coord,Lbc,Latcons,nlayers, nat,A,smin,smax,snum,mind) ;
  } else if ( pair_type == INVERSE_R ) {
    ZCalc_InvR_Deriv(Coord,Lbc,Latcons,nlayers, nat,A,smin,smax,snum,mind) ;
  } else {
    cout << "Error: bad pairtype in ZCalc_Deriv\n" ;
    exit(1) ;
  }
  
  if ( if_3b_cheby ) {
    ZCalc_3B_Cheby_Deriv(Coord, Lbc, Latcons, nat, A, smin, smax, snum, 
			 snum_3b_cheby, lambda) ;
  }

}

static void ZCalc_Spline_Deriv(double **Coord, const char *Lbc, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double *smin,
			       const double *smax,
			       const double *sdelta,const int *snum,double *mind)
// Calculate derivatives of the forces wrt the spline parameters.
// Stores minimum distance between a pair of atoms in mind.
{
  double Rvec[3];
  double Rab[3];
  double rlen;
  double t;
  double h00,h10,h01,h11;
  int k0;
  double x,x0;
  int vstart,kstart;

  ////main loop for SPLINE terms:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	//calculate vstart: (index for populating OO, OH, or HH column block of A).
	//calculate vstart:
	int ipair = pair_index(a1,a2,Lbc) ;
	vstart = pair_table_offset(ipair,snum) ;

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
		
		if ( rlen < mind[ipair] ) mind[ipair] = rlen ;

		//spline term calculated w/cutoff:
		if(rlen > smin[ipair] and rlen < smax[ipair])
		  {

		    //here is the setup for Hermite cubic quadrature:
		    k0=int(floor((rlen-smin[ipair])/sdelta[ipair]));

		    // Keep from overrunning table.
		    if ( k0 >= snum[ipair]/2 - 1 ) {
		      cout << "Table overrun: rlen = " << rlen << endl ;
		      continue ;
		    }

		    x=rlen;
		    
		    x0=smin[ipair]+sdelta[ipair]*k0;
		    t=(x-x0)/sdelta[ipair];

		    //for classical cubic quadrature, you would only
		    //change h00, h10, h01, h11 polynomials!
		    h00=2*t*t*t-3*t*t+1;
		    h10=t*t*t-2*t*t+t;
		    h01=-2*t*t*t+3*t*t;
		    h11=t*t*t-t*t;
		    h10 *= sdelta[ipair];//derivative terms have extra factor.
		    h11 *= sdelta[ipair];

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

static void ZCalc_Poly_Deriv(double **Coord,const char *Lbc, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double *smin,
			       const double *smax,
			      const int *snum, double *mind)
// Calculate derivatives of the forces wrt the DFTB Erep parameters.
{
  double Rvec[3];
  double Rab[3];
  double rlen;
  double x;
  int vstart ;
  static bool called_before = false ;
  const double autoang = 0.5291772488820865;
  static double rc[NPAIR];

  if ( ! called_before ) {
    called_before = true ;
    int dim = 0 ;
    for ( int i = 0 ; i < NPAIR ; i++ ) 
      {
        rc[i] = smax[i]/autoang;
        printf("rc[%d] = %lf\n",i,rc[i]);
	if ( snum[i] > dim ) 
	  {
	    dim = snum[i] ;
	  }
      }
    dim++ ;
  }

  ////main loop for ninth order polynomial terms:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	//calculate vstart: (index for populating OO, OH, or HH column block of A).
	int ipair = pair_index(a1, a2, Lbc) ;
	vstart = pair_table_offset(ipair,snum) ;

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
		
		if ( rlen < mind[ipair] ) mind[ipair] = rlen ;

		if(rlen>smin[ipair] and rlen<smax[ipair])
		  {
		    //calculate binning, convert all distances to au
		    x=rlen/autoang;
	            for ( int i = 0 ; i < snum[ipair]; i++ ) 
		      {
                        double rfac = -(i+2)*pow((rc[ipair]-x),i+1);
		        for(int c=0;c<3;c++)
		          {
			    A[a1][vstart+i][c] += rfac*Rab[c]/rlen;
			    A[a2][vstart+i][c] -= rfac*Rab[c]/rlen;
		          }
                      }
		  }//rlen
	      }
      }
  return;
}

static void ZCalc_InvR_Deriv(double **Coord,const char *Lbc, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double *smin,
			       const double *smax,
			      const int *snum, double *mind)
// Calculate derivatives of the forces wrt the Chebyshev parameters.
// Stores minimum distance between a pair of atoms in mind.
{
  double Rvec[3];
  double Rab[3];
  double rlen;
  double x;
  int vstart ;

  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	//calculate vstart: (index for populating OO, OH, or HH column block of A).
	int ipair = pair_index(a1, a2, Lbc) ;
	vstart = pair_table_offset(ipair,snum) ;

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
		
		if ( rlen < mind[ipair] ) mind[ipair] = rlen ;

		if(rlen>smin[ipair] and rlen<smax[ipair])
		  {
		    x=rlen;
                    double fc = pow((rlen-smax[ipair]),3);
                    double dfc = 3*pow((rlen-smax[ipair]),2);
 	            for ( int i = 0 ; i < snum[ipair]; i++ ) 
		      {
                        double rfac = ((i+2)/pow(rlen,i+3))*fc;
                        rfac -= (1/pow(rlen,i+2))*dfc;
		        for(int c=0;c<3;c++)
		          {
			    A[a1][vstart+i][c] += rfac*Rab[c]/rlen;
			    A[a2][vstart+i][c] -= rfac*Rab[c]/rlen;
		          }
                      }
		  }//rlen
	      }
      }
  return;
}

static void ZCalc_Cheby_Deriv(double **Coord,const char *Lbc, double *Latcons,
			       const int nlayers,
			       const int nat,double ***A,const double *smin,
			       const double *smax,
			      const int *snum, const double *lambda, double *mind)
// Calculate derivatives of the forces wrt the Chebyshev parameters.
// Stores minimum distance between a pair of atoms in mind.
{
  double Rvec[3];
  double Rab[3];
  double rlen;
  double x;
  int vstart ;
  static double *Tn, *Tnd ;
  static bool called_before = false ;

  if ( ! called_before ) {
    called_before = true ;
    int dim = 0 ;
    for ( int i = 0 ; i < NPAIR ; i++ ) 
      {
	if ( snum[i] > dim ) 
	  {
	    dim = snum[i] ;
	  }
      }
    dim++ ;
    Tn   = new double [dim] ;
    Tnd  = new double [dim] ;
  }

  ////main loop for Chebyshev terms:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	//calculate vstart: (index for populating OO, OH, or HH column block of A).
	int ipair = pair_index(a1, a2, Lbc) ;
	vstart = pair_table_offset(ipair,snum) ;

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
		
		if ( rlen < mind[ipair] ) mind[ipair] = rlen ;

		//spline term calculated w/cutoff:
		bool inverse_r = false ;
		bool morse     = true ;

		if(rlen>smin[ipair] and rlen<smax[ipair])
		  {
		    double xavg, xdiff ;

		    if ( inverse_r ) {
		      xavg = 0.5 * (1.0/smin[ipair] + 1.0/smax[ipair]) ;
		      xdiff = 0.5 * (1.0/smin[ipair] - 1.0/smax[ipair]) ;
		      x = (1.0/rlen-xavg) / xdiff ;
		    } else if ( morse ) {
		      double xmin = exp(-smax[ipair]/lambda[ipair]) ;
		      double xmax = exp(-smin[ipair]/lambda[ipair]) ;
		      xavg = 0.5 * (xmin + xmax) ;
		      xdiff = 0.5 * (xmax - xmin) ;
		      x = (exp(-rlen/lambda[ipair])-xavg)/xdiff ;
		    }
		    else {
		      xavg = 0.5 * (smin[ipair] + smax[ipair]) ;
		      xdiff = 0.5 * (smax[ipair] - smin[ipair]) ;
		      x = (rlen-xavg) / xdiff ;
		    }

		    // x = (rlen-smin[ipair]) / (smax[ipair] - smin[ipair]) ;

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
		    for ( int i = 2 ; i <= snum[ipair] ; i++ ) {
		      Tn[i] =  2.0 * x * Tn[i-1] - Tn[i-2] ;
		      Tnd[i] = 2.0 * x * Tnd[i-1] - Tnd[i-2] ;
		    }
		    // Now store dTn/dx = n * U_(n-1)
		    for ( int i = snum[ipair] ; i >= 1 ; i-- ) {
		      Tnd[i] = i * Tnd[i-1] ;
		    }
		    Tnd[0] = 0.0 ;

		    double rlen3 = rlen * rlen * rlen ;
		    if ( morse ) 
		      {
			double fcut0 = (1.0 - rlen/smax[ipair]) ;
			double fcut = fcut0 * fcut0 * fcut0 ;
			double fcutderiv = -3.0 * fcut0 * fcut0 / smax[ipair] ;
			for ( int i = 0 ; i < snum[ipair] ; i++ ) 
			  {
			    for(int c=0;c<3;c++)
			      {
				double deriv = 
				  (fcut * Tnd[i+1] *(-exp(-rlen/lambda[ipair])/lambda[ipair])/xdiff + 
				   fcutderiv * Tn[i+1]
				   ) * Rab[c] / rlen ; 
				A[a1][vstart+i][c] += deriv ;
				A[a2][vstart+i][c] -= deriv ;
			      } 
			  }
		      }
		    else 
		      {
			for ( int i = 0 ; i < snum[ipair] ; i++ ) 
			  {
			    for(int c=0;c<3;c++) 
			      {
				A[a1][vstart+i][c] += (smax[ipair] - rlen) * Tn[i]*Rab[c]/rlen3;
				A[a2][vstart+i][c] -= (smax[ipair] - rlen) * Tn[i]*Rab[c]/rlen3;
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
			 const int nat, bool calc_deriv, 
			 double **Fderiv, int n_over, double *over_param)
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


  // Lucas paper draft parameters.
  // p1=  -2.5881042987450e-01   ; 
  // r0=  9.6000387075695e-01    ;
  // p2=  3.8995379237250e+00    ;
  // lambda6=     -8.9                 ;  
 
  // Chenoweth values for p1, p2, r0, lambda6 p(ovun2)
  // pover=50.0;
  // p1 = -0.0657 ;
  // p2 = 5.0451 ;
  // r0 = 1.0165 ;
  // lambda6 = -3.6141 ;
  
  if ( n_over != 5 ) 
    {
      cout << "Error: wrong number of overcoordination parameters\n" ;
      exit(1) ;
    }
  
  pover = over_param[0] ;
  r0    = over_param[1] ;
  p1    = over_param[2] ;
  p2   =  over_param[3] ;
  lambda6 = over_param[4] ;

  
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
		  if(aj==ai) {
		    dEover[aj][c]-=tempr*Rab[c]/rik;
                  }
		  if(aj==ak) {
		    dEover[aj][c]+=tempr*Rab[c]/rik;
                  }
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



bool parse_tf(char *val, int bufsz, char *line)
// Parse a true or false token argument.
{

  // val[strlen(val)-1] = 0 ;
  //  printf("val = :%s:\n", val) ;
  if ( strncmp(val, "true",bufsz) == 0 ) 
    {
      return(true) ;
    }
  else if ( strncmp(val, "false", bufsz) == 0 ) 
    {
      return(false) ;
    }
  else 
    {
      printf("Error: need a true or false argument: %s\n",
	     line) ;
    }  
  return(false) ;
}


void parse_param_list(double *params, int nparams, const char* name)
// Parse a given number of floating point parameters on a single input line.
{
  char *val ;

  for ( int k = 0 ; k < nparams ; k++ ) 
    {
      val = strtok(NULL, " ") ;
      if ( val == NULL || val[0] == 0 ) 
	{
	  cout << "Error: " << name <<  " parameter " <<
	    k << " not found\n" ;
	  exit(1) ;
	}
      else 
	{
	  sscanf(val, "%lf", &params[k]) ;
	}
    }
}

int pair_index(int a1, int a2, const char *Lbc)
// Returns the index of the pair type corresponding to atoms a1 and a2.
{
  if ( Lbc[a1] == 'O' && Lbc[a2] == 'O' ) 
    {
      return(0) ;
    } 
  else if ( Lbc[a1] == 'O' && Lbc[a2] == 'H' ) 
    {
      return(1) ;
    }
  else if ( Lbc[a1] == 'H' && Lbc[a2] == 'O' ) 
    {
      return(1) ;
    }
  else if ( Lbc[a1] == 'H' && Lbc[a2] == 'H' ) 
    {
      return(2) ;
    }
  else
    {
      cout << "Error: bad atom labels for atoms " << a1 << " and " << a2 << endl ;
      exit(1) ;
    }
}

int atom_index(int a1, const char *Lbc)
// Returns the index of the pair type corresponding to atoms a1 and a2.
{
  if ( Lbc[a1] == 'O' ) 
    {
      return(0) ;
    } 
  else if ( Lbc[a1] == 'H' )
    {
      return(1) ;
    }
  else
    {
      cout << "Error: bad atom label for atom " << a1 << endl ;
      exit(1) ;
    }
}

static int pair_table_offset(int ipair, const int *snum)
// Calculate the index of the parameters for the given interaction pair.
{
  int vstart = 0 ;
  for ( int i = ipair - 1 ; i >= 0 ; i-- ) {
    vstart += snum[i] ;
  }
  return(vstart) ;
}

Sr_pair_t parse_pair_type(const char *val, int bufsz) 
// Parse the pair type from a string argument.
{

  Sr_pair_t pair_type = CHEBYSHEV ;

  if ( strncmp(val, "spline", bufsz) == 0 ) 
    {
      pair_type = SPLINE ;
    } 
  else if ( strncmp(val, "chebyshev", bufsz) == 0 ) 
    {
      pair_type = CHEBYSHEV ;
    }
  else if ( strncmp(val, "inverse_r", bufsz) == 0 ) 
    {
      pair_type = INVERSE_R ;
    }
  else if ( strncmp(val, "lennard_jones", bufsz) == 0 ) 
    {
      pair_type = LJ ;
    }
  else if ( strncmp(val, "stillinger", bufsz) == 0 ) 
    {
      pair_type = STILLINGER ;
    }
  else 
    {
      printf("Error: did not recognize pair_type |%s|\n", val) ;
      exit(1) ;
    }
  return(pair_type) ;
}

bool read_tf_option(ifstream *paramread, const char *option, const char *params_file)
// Read a true/false option from the given file.
{
  const int bufsz = 1024 ;
  char buf[bufsz] ;
  char cmd_arg[bufsz] ;
  bool result = false ;

  int optlen = strlen(option) ;
  paramread->getline(buf, bufsz) ;

  if ( ! strncmp(buf, option, optlen ) ) {
    sscanf(buf+optlen+1, "%s", cmd_arg) ;
    result = parse_tf(cmd_arg, bufsz, buf) ;
    if ( result ) {
      printf("Option %s is set to true\n", option) ;
    } else {
      printf("Option %s is set to false\n", option); 
    }
  } else {
    printf("Error in reading %s: %s not found\n", params_file, option) ;
    printf("Line was: %s\n", buf) ;
    exit(1) ;
  }
  return(result) ;
}
