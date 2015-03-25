#include "functions.h"

static void ZCalc_Lj(double **Coord,string *Lb, double *Latcons,const int nlayers,
	   const int nat,const double smin,const double smax,
	   const double sdelta,const int snum, 
		     double *params,double **SForce,double& Vtot,double& Pxyz) ;
static void ZCalc_SR_Analytic(double **Coord,char *Lbc, double *Latcons,const int nlayers,
			      const int nat,const double smin,const double smax,
			      double **SForce,double& Vtot,double& Pxyz) ;
static void ZCalcSR_Over(double **Coord, char *Lbc, double *Latcons,const int nlayers,
			 const int nat,double **SForce,double& Vtot,double& Pxyz, double pover) ;
static void ZCalc_Stillinger(double **Coord, char *Lbc, double *Latcons,const int nlayers,
			     const int nat, double smax, 
			     double **SForce,double& Vtot,double& Pxyz) ;
static void ZCalc_Spline(double **Coord,string *Lb, double *Latcons,const int nlayers,
			 const int nat,const double smin,const double smax,
			 const double sdelta,const int snum, 
			 double *params, double *pot_params,
			 double **SForce,double& Vtot,double& Pxyz) ;
static void ZCalc_Cheby(double **Coord,string *Lb, double *Latcons,
			const int nlayers,
			const int nat,const double smin,
			const double smax,
			const double sdelta,const int snum, 
			double *params,
			double **SForce, double &Vtot, double &xyz);

void ZCalc(double **Coord, string *Lb, double *Q, double *Latcons,const int nlayers,
	   const int nat,const double smin,const double smax,
	   const double sdelta,const int snum, 
	   double *params, double *pot_params, bool if_cheby,
	   double **SForce,double& Vtot,double& Pxyz)
{
  bool ifcalc_ewald = true ;       // Ewald charge evaluation
  bool ifcalc_over  = true ;      // Short-range overcoordination potential

  // Choose either spline, cheby, or analytic.
  bool ifcalc_spline     = true ;  // Spline pair force model.
  bool ifcalcsr_analytic = false ; // Short-range analytic potential evaluation.

  // These are for testing purposes only.
  bool ifcalc_lj = false ;           // Lennard-Jones potential for testing.
  bool ifcalc_stillinger = false ;   // Stillinger central force model.

  double pover ;
  double volume ;

  char *Lbc ;

#if(0)
  if ( if_cheby ) {
    ifcalc_spline = false ;
  } else {
    ifcalc_spline = true ;
  }
#else
  // DEBUG !
  ifcalc_spline = false ;
  if_cheby  = false ;
  ifcalc_stillinger = true ;
  ifcalc_lj     = false ;
  ifcalc_ewald  = false ;
  ifcalc_over   = false ;
#endif

  // Pack element names into single characters for efficiency
  Lbc = new char [nat] ;
  for(int a=0;a<nat;a++) {
    if ( Lb[a] == "O" ) {
      Lbc[a] = 'O' ;
    } else if ( Lb[a] == "H" ) {
      Lbc[a] = 'H' ;
    } else {
      cout << "Error: unknown element " << Lb[a] << "\n" ;
    }
  }

  for(int a=0;a<nat;a++)
    for(int c=0;c<3;c++)
      SForce[a][c]=0;

  Vtot = 0.0;
  Pxyz = 0.0 ;

  if ( ifcalc_lj ) {
    ZCalc_Lj(Coord, Lb, Latcons,nlayers, nat, smin, smax, sdelta,snum, 
	     params, SForce, Vtot, Pxyz) ;
  }

  if ( ifcalc_ewald ) {
    ZCalc_Ewald(Coord, Lb, Q, Latcons,nlayers, nat, smin, smax, sdelta,snum, 
		params, SForce, Vtot, Pxyz) ;
  }
  if ( ifcalcsr_analytic ) {
    ZCalc_SR_Analytic(Coord,Lbc, Latcons,nlayers,nat,smin,smax,SForce,Vtot,Pxyz) ;
  }
  if ( ifcalc_over ) {
    pover = params[snum+3] ;
    // cout << "POVER = " << pover << endl ;
    ZCalcSR_Over(Coord,Lbc, Latcons,nlayers, nat,SForce, Vtot, Pxyz,pover) ;
  }

  if ( ifcalc_stillinger ) {
    ZCalc_Stillinger(Coord,Lbc, Latcons,nlayers,nat,smax, SForce,Vtot,Pxyz) ;
  }

  if ( ifcalc_spline ) {
    ZCalc_Spline(Coord,Lb,Latcons,nlayers,nat,smin,smax,sdelta,snum, 
		 params,pot_params,SForce,Vtot,Pxyz) ;
  }
  if ( if_cheby ) {
    ZCalc_Cheby(Coord,Lb,Latcons,nlayers,nat,smin,smax,sdelta,snum, 
		 params,SForce,Vtot,Pxyz) ;
  }
  delete Lbc ;

  volume = Latcons[0] * Latcons[1] * Latcons[2] ;

  Pxyz /= 3.0 * volume ;

  return;
}


static void ZCalc_Lj(double **Coord,string *Lb, double *Latcons,const int nlayers,
	   const int nat,const double smin,const double smax,
	   const double sdelta,const int snum, 
	   double *params,double **SForce,double& Vtot,double& Pxyz)
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


static void ZCalc_SR_Analytic(double **Coord,char *Lbc, double *Latcons,const int nlayers,
			      const int nat,const double smin,const double smax,
			      double **SForce,double& Vtot,double& Pxyz)
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

  bool ifnonbond;

  double pOO[24] ;
  
  pOO[0]=  -18617.280354059098    	  ;
  pOO[1]=  692382.49446701177           ;
  pOO[2]=  -11405958.131264212    	  ;
  pOO[3]=  109588187.69652703           ;
  pOO[4]=  -678681472.79666793    	  ;
  pOO[5]=  2813120111.5392013           ;
  pOO[6]=  -7774279189.7940693    	  ;
  pOO[7]=  13490534921.663717           ;
  pOO[8]=  -11633893383.8025      	  ;
  pOO[9]=  -2281518444.80867      	  ;
  pOO[10]= 12083725865.839197           ;
  pOO[11]= 45862665.384337671           ;
  pOO[12]= -11788101437.910231    	  ;
  pOO[13]= -3101264034.3021011    	  ;
  pOO[14]= 10798442968.405247           ;
  pOO[15]= 8739722513.8273964           ;
  pOO[16]= -7388244851.8771982    	  ;
  pOO[17]= -14539511339.557861    	  ;
  pOO[18]= 11471144515.535652           ;
  pOO[19]= 0.0  ;
  pOO[20]= 0.0  ;
  pOO[21]= 0.0  ;
  pOO[22]= 0.0  ;
  pOO[23]= 0.0  ;

  for(int n=0;n<24;n++)
    pOO[n]*=-1.0;

  double pOH[24] ;
			
  pOH[0]=    -723.98843365906146   	 ;
  pOH[1]=    19521.464403077367    	 ;
  pOH[2]=    -239387.19844904722   	 ;
  pOH[3]=    1779971.391880197     	 ;
  pOH[4]=    -9012011.3066581544   	 ;
  pOH[5]=    33045520.790334914    	 ;
  pOH[6]=    -91065426.852874398   	 ;
  pOH[7]=    193063512.18800989    	 ;
  pOH[8]=    -319567836.69655925   	 ;
  pOH[9]=    416549308.66748077    	 ;
  pOH[10]=   -429085671.67467254   	 ;
  pOH[11]=   348884311.45362383    	 ;
  pOH[12]=   -222517723.55597138   	 ;
  pOH[13]=   109947020.4432832     	 ;
  pOH[14]=   -41208686.826725334   	 ;
  pOH[15]=   11317859.397690322    	 ;
  pOH[16]=   -2147264.8830902833   	 ;
  pOH[17]=   251388.67261027446    	 ;
  pOH[18]=   -13681.419073463892     ;
  pOH[19]= 0.0 ;
  pOH[20]= 0.0 ;
  pOH[21]= 0.0 ;
  pOH[22]= 0.0 ;
  pOH[23]= 0.0 ;

  for(int n=0;n<24;n++)
    pOH[n]*=-1;


  double pHH[24] ;


  pHH[0]=  4008.2237405234468     ;	     
  pHH[1]=  -133312.7298043904     ;	     
  pHH[2]=  2029736.1394321355     ;	     
  pHH[3]=  -18810746.964217555    ;	     
  pHH[4]=  119027276.85909875     ;	     
  pHH[5]=  -546712328.86844385    ;	     
  pHH[6]=  1891487808.0979228     ;	     
  pHH[7]=  -5046788818.0488892    ;	     
  pHH[8]=  10541569604.029257     ;	     
  pHH[9]=  -17389390868.432339    ;	     
  pHH[10]= 22737335322.912605     ;	     
  pHH[11]= -23538172961.926105    ;	     
  pHH[12]= 19171310483.05817      ;	     
  pHH[13]= -12131668715.262671    ;	     
  pHH[14]= 5839324812.613658       ;			     
  pHH[15]= -2064811667.1446371     ;			     
  pHH[16]= 505544757.04635936      ;			     
  pHH[17]= -76541863.838855445     ;			     
  pHH[18]= 5397474.2715492537      ;			     
  pHH[19]=  0.0  ;			     
  pHH[20]=  0.0  ;			     
  pHH[21]=  0.0  ;			     
  pHH[22]=  0.0  ;			     
  pHH[23]=  0.0  ;                       

  for(int n=0;n<24;n++)
    pHH[n]*=-1.0;
					     

  ////main loop non-Coulomb short-ranged forces:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++)
      {
	ifnonbond=true;
        if(Lbc[a1]=='O' and Lbc[a2]=='H')
          if(a2<=a1+2)
            ifnonbond=false;
        if(Lbc[a1]=='H' and Lbc[a2]=='H')
          if(a2<=a1+1)
            ifnonbond=false;

	ifnonbond=true;
	if(ifnonbond)
	  {
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

		    if(rlen > smin * smin and rlen < smax * smax)
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
			if((Lbc[a1]=='O' and Lbc[a2]=='H') or (Lbc[a1]=='H' and Lbc[a2]=='O'))
			  {
			    rlen2=rlen;
			    if(rlen2<oh_cut)
			      {
				rlen2=oh_cut;
			      }
			  }
			if(Lbc[a1]=='H' and Lbc[a2]=='H')
			  {
			    rlen2=rlen;
			    if(rlen2<hh_cut)
			      {
				rlen2=hh_cut;
			      }
			  }
		    
			double smaxpow3 = pow(rlen2-smax,3) ;
			double smaxpow2 = pow(rlen2-smax,2) ;
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
      }
  
  return;
}



static void ZCalcSR_Over(double **Coord, char *Lbc, double *Latcons,const int nlayers,
			 const int nat,double **SForce,double& Vtot,double& Pxyz, double pover)
// Calculate short-range overcoordination forces.
{
  double Rvec[3], Rab[3] ;
  
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
  double p1,p2,r0,lambda6;
  double rik,tempr,temps;

  //ReaxFF over-coordination parameters (from non-linear fitting).
  // pover=75.0;
  // p1=  -2.5881042987450e-01   ; 
  // r0=  9.6000387075695e-01    ;
  // p2=  3.8995379237250e+00    ;
  // lambda6=-8.9;

  // Chenoweth values for p1, p2, r0, lambda6 p(ovun2)
  p1 = -0.0657 ;
  p2 = 5.0451 ;
  r0 = 1.0165 ;
  lambda6 = -3.6141 ;

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

static void ZCalc_Stillinger(double **Coord,char *Lbc, double *Latcons,const int nlayers,
			     const int nat, double smax, 
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

		if( rlen < smax )
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


static void ZCalc_Spline(double **Coord,string *Lb, double *Latcons,const int nlayers,
			 const int nat,const double smin,const double smax,
			 const double sdelta,const int snum, 
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
	if(Lb[a1]=="O" and Lb[a2]=="O")
	  vstart=0;
	if((Lb[a1]=="O" and Lb[a2]=="H") or (Lb[a1]=="H" and Lb[a2]=="O"))
	  vstart=snum/3;
	if(Lb[a1]=="H" and Lb[a2]=="H")
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
		    S_r=0.0;

		    rlen2 = rlen ;
		    tempx = spline_pot(smin, smax, sdelta, rlen2, params, pot_params, snum, vstart, S_r) ;

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
  if ( kstart > snum / 3 ) {
    cout << "Error: kstart too large " << kstart << endl ;
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
  if ( kstart/2 + 1 < snum / 6 ) {
    tempx += pot_params[vstart/2+kstart/2+1] ;
  }
    return(tempx) ;
}


static void ZCalc_Cheby(double **Coord,string *Lb, double *Latcons,
			const int nlayers,
			const int nat,const double smin,
			const double smax,
			const double sdelta,const int snum, 
			double *params,
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
  static char *Lbc ;
  //spline term calculated w/cutoff:
  const bool inverse_r = false ;
  const bool morse     = true ;
  double tempx;
  const double lambda = 1.25 ;
  double exprlen ;

  if ( ! called_before ) {
    called_before = true ;
    Lbc = new char [nat] ;
    Tn   = new double [snum/3+1] ;
    Tnd  = new double [snum/3+1] ;
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
		      exprlen = exp(-rlen/lambda) ;
		      x = (exprlen-xavg)/xdiff ;
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

		    if ( morse ) 
		      {
			double fcut0 = (1.0 - rlen/smax) ;
			double fcut = fcut0 * fcut0 * fcut0 ;
			double fcutderiv = -3.0 * fcut0 * fcut0 / smax ;
			tempx = 0.0 ;
			for ( int i = 0 ; i < snum / 3 ; i++ ) 
			  {
			    double coeff = params[vstart+i] ;
			    tempx += coeff * fcut * Tn[i+1] ;
			    double deriv = 
			      (fcut * Tnd[i+1] *(-exprlen/lambda)/xdiff + 
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
			for ( int i = 0 ; i < snum / 3 ; i++ ) 
			  {
			    double coeff = params[vstart+i] ;
			    tempx = 0.0 ;
			    // POTENTIAL NOT YET IMPLEMENTED.
			    Pxyz += coeff * (smax-rlen) * Tn[i] / rlen ;
			    {
			      for(int c=0;c<3;c++) 
				{
				  SForce[a1][c] += coeff * (smax - rlen) * Tn[i]*Rab[c]/rlen3;
				  SForce[a2][c] -= coeff * (smax - rlen) * Tn[i]*Rab[c]/rlen3;
				}
			    }
			  }
		      }
		    Vtot += tempx ;
		  }//rlen
	      }
      }
  return;
}
