#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include "functions.h"
using namespace std;

static void read_input(const char *file_name,
		       double &TempMD, double &deltat_fs, int &nsteps, 
		       int &nlayers, 
		       bool &if_output_force, bool &if_read_force,
		       bool &if_spline_q, bool &if_init_vel, int &rand_seed, 
		       int &gen_freq, int &energy_freq, int &scale_freq,
		       double &thoover_fs, Sr_pair_t &pair_type,
		       bool &if_coulomb,bool &num_pressure, char *params_file,
		       bool &if_overcoord) ;

static void echo_input(double TempMD, double deltat_fs, int nsteps, 
		       int nlayers, 
		       bool if_output_force, bool if_read_force,
		       bool if_spline_q, bool if_init_vel, int rand_seed, 
		       int gen_freq, int energy_freq, int scale_freq,
		       double thoover_fs, Sr_pair_t pair_type,
		       bool if_coulomb, bool num_pressure, const char *params_file,
		       bool if_overcoord) ; 


static double kinetic_energy(double *Mass, double **Vel, int nat) ;

static double 
numerical_pressure(double **Coord, string *Lb, double *Q, double *Latcons,
		   const int nlayers, const int nat,const double smin,
		   const double smax, const double sdelta,const int snum, 
		   double *params, double *pot_params, Sr_pair_t pair_type,
		   bool if_coulomb,bool if_overcoord, int n_over,
		   double *over_param) ;

int main(int argc, char* argv[])
{
  ////Job parameters:
  /// These are defaults.  Values are read from spline_md.in.

  //temperature
  double TempMD=2000.0;

  //MD timestep in fs.
  double deltat_fs=0.125;

  // Number of MD steps
  int nsteps=80000 ;

  // Number of periodic images to use. Usually, no need to go past 1 
  // (depends only on spline interaction cutoff).
  int nlayers=1;

  // If TRUE, write out calculated forces.
  bool if_output_force = false ;

  // If TRUE, read an xyzf file and compare to forces in xyzf file.
  bool if_read_force = false ;

  // If true, take charges from spline parameters.
  bool if_spline_q = true ;

  // If true, initialize velocities.
  bool if_init_vel = false ;

  // Random number generator seed.
  int rand_seed = 123457 ;

  // If true, output coordinates in DFTB gen format for use with molanal.
  bool output_gen = false ;

  // If true, use hoover-thermostat
  bool if_hoover = false ;

  // If true, calculate Coulomb forces.
  bool if_coulomb = true ;

  // If true, calculate ReaxFF-like overcoordination term.
  bool if_overcoord = true ;

  // How often to output energy
  int energy_freq = 10 ;

  // How often to scale velocities (0 if never).
  int scale_freq = 0 ;  

  // How often to write the gen file.
  int gen_freq = 10 ;

  // Whether to calculate pressures by finite difference.
  bool num_pressure = false ;

  char params_file[1024] = "params.txt" ;

  double over_param[MAXOVERP] =   {
    // Lucas's parameters for overcoordination
    // 75.0,                    // pover
    // 9.6000387075695e-01,     // r0
    // -2.5881042987450e-01,    // p1
    // 3.8995379237250e+00,     // p2
    // -8.9                     // lambda6
    //
    // Chenoweth values
    50.0,      // pover
    1.0165,    // r0
    -0.0657,   // p1
    5.0451,    // p2
    -3.6141    // lambda6
  } ;

  int n_over = 5 ;

  double thoover_fs = 0.0 ;

  Sr_pair_t pair_type = CHEBYSHEV ;

  if ( argc != 2 ) 
    {
      cout << "Usage:  splines_md <input_file>\n" ;
      exit(1) ;
    }

  read_input(argv[1],TempMD, deltat_fs, nsteps, nlayers, if_output_force, if_read_force,
	     if_spline_q, if_init_vel, rand_seed, gen_freq, energy_freq, scale_freq,
	     thoover_fs,pair_type,if_coulomb,num_pressure,params_file,
	     if_overcoord) ;

  echo_input(TempMD, deltat_fs, nsteps, nlayers, if_output_force, if_read_force,
	     if_spline_q, if_init_vel, rand_seed, gen_freq, energy_freq, scale_freq,
	     thoover_fs,pair_type,if_coulomb,num_pressure,params_file,
	     if_overcoord) ;

  if ( thoover_fs > 0.0 ) if_hoover = true ;
  if ( gen_freq > 0 ) output_gen = true ;

  // END OF INPUT

  const double deltat=deltat_fs/Tfs ;//time units for E=kcal/mol,mass=amu,distance=Angs
  
  cout.precision(15);

  ////Read input file input.xyz:
  ////xyz file in Angs with box dimensions in infoline.
  ifstream fileread ;
  if ( if_read_force ) 
    {
      fileread.open("input.xyzf") ;
    } 
  else 
    {
      fileread.open("input.xyz");
    }
  int nat;
  double Latcons[3];
  fileread >> nat;

  fileread >>Latcons[0]>>Latcons[1]>>Latcons[2];
  string *Lb ;
  double **Coord;
  double **Accel;
  double **Fread ;
  double **Vel;
  double **Velnew;
  double Q[nat];  // Charge on each atom.

  // Hoover viscosity.
  double zeta = 0.0 ;

  // Hoover coordinate.
  double s_hoover = 0.0 ; 

  Lb = new string[nat] ;
  Coord=new double *[nat];
  Accel=new double *[nat];
  Fread=new double *[nat] ;
  Vel  = new double *[nat] ;
  Velnew  = new double *[nat] ;

  FILE *fgen ;

  if ( output_gen ) 
    {
      fgen = fopen("traj.gen", "w") ;
      if ( fgen == NULL ) 
	{
	  cout << "Error: could not open traj.xyz\n" ;
	  exit(1) ;
	}
    }
  else
    {
      fgen = NULL ;
    }

  if ( if_read_force ) nsteps = 1 ;

  double thoover = thoover_fs / Tfs ;
  int ndof = 3 * nat - 3 ;
  double QHoover = ndof * Kb * TempMD * thoover * thoover ;

  if ( if_hoover ) {
    // Don't scale velocities when hoover thermostat is on.
    scale_freq = 0 ;
  }

  for(int a=0;a<nat;a++)
    {
      Coord[a] = new double[3];
      Accel[a] = new double[3];
      Fread[a] = new double[3] ;
      Vel[a]   = new double[3] ;
      Velnew[a]   = new double[3] ;
      for(int c=0;c<3;c++)
	{
	  Coord[a][c]=0.0;
	  Accel[a][c]=0.0;
	  Fread[a][c] = 0.0 ;
	  Vel[a][c]   = 0.0 ;
	  Velnew[a][c]   = 0.0 ;
	}
    }

  for(int a=0;a<nat;a++)
    {
      fileread>>Lb[a];
      for(int c=0;c<3;c++)
	fileread>>Coord[a][c];

      if ( if_read_force ) 
	{
	  for(int c=0;c<3;c++) 
	    {
	      fileread>>Fread[a][c];
	      Vel[a][c] = 0.0 ;
	    }
	}
      else
	{
	  for(int c=0;c<3;c++) 
	    {
	      fileread>>Vel[a][c];
	      Fread[a][c] = 0.0 ;
	    }
	}
    }

  //convert OpenMX velocities:
  //for(int a=0;a<nat;a++)
  //for(int c=0;c<3;c++)
  //Vel[a][c] *= 0.00048071145295036654*0.545;


  double Mass[nat];

  for(int a=0;a<nat;a++)
    {
      if(Lb[a]=="O") 
	{
	  Mass[a]=15.999;
	}
      else if(Lb[a]=="H") 
	{
	  Mass[a]=1.008;
	}
      else
	{
	  cout<<"Atom type "<<Lb[a]<<" not supported."<<endl;
	  exit(1);
	}
    }


  ////Wrap the coordinates:
  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      {
	Coord[a1][c] -= floor(Coord[a1][c]/Latcons[c])*Latcons[c];
      }


  ////Define some variables.
  double Ktot,Vtot,Pxyz;

  double Vol = Latcons[0] * Latcons[1] * Latcons[2] ;


  ////Generate random initial velocities w/Box Muller:
  //double Vel[nat][3];
  if ( if_init_vel ) {
    double x1,x2,y1,y2;
    double sigma;//=sqrt(TempMD*0.1);
    for(int a=0;a<nat;a++)
      for(int c=0;c<3;c++)
	Vel[a][c]=0.0;
    srand(rand_seed);
    for(int a=0;a<3*nat;a++)
      {
	int iatm = a/3 ;
	int ixyz = a % 3 ;

	sigma=sqrt(TempMD*Kb/Mass[iatm]) ;
	x1=double(rand())/double(RAND_MAX);
	x2=double(rand())/double(RAND_MAX);
	if ( x1 < 0.0 || x1 > 1.0 ) {
	  cout << "Bad random variable\n" ;
	  exit(1) ;
	}
	if ( x2 < 0.0 || x2 > 1.0 ) {
	  cout << "Bad random variable\n" ;
	  exit(1) ;
	}
	y1 = sqrt(-2.0 * log(x1)) * cos(2.0 * M_PI * x2) ;
	y2 = sqrt(-2.0 * log(x1)) * sin(2.0 * M_PI * x2) ;
	y1=y1*sigma ;
	y2=y2*sigma ;

	if ( ixyz < 0 || ixyz > 2 ) {
	  cout << "Bad xyz index\n" ;
	  exit(1) ;
	}
	
	// Use either y2 or y1 here for maximum fun.
	if ( a % 2 == 0 ) 
	  {
	    Vel[iatm][ixyz]=y2 ;
	  } 
	else 
	  {
	    Vel[iatm][ixyz]=y1 ;
	  }
      }
  }
  ////test/fix v.c.o.m. here
  double tempx=0.0;
  double tempy=0.0;
  double tempz=0.0;
  double tempm=0.0 ;

  for(int a=0;a<nat;a++)
    {
      tempx+=Mass[a]*Vel[a][0];
      tempy+=Mass[a]*Vel[a][1];
      tempz+=Mass[a]*Vel[a][2];
      tempm+=Mass[a];
    }
	  
  cout<<"v.c.o.m. = "<<tempx/tempm<<"   "<<tempy/tempm<<"   "<<tempz/tempm<<endl;
  
  //remove:
  for(int a=0;a<nat;a++)
    {
      Vel[a][0] -= tempx/tempm;
      Vel[a][1] -= tempy/tempm;
      Vel[a][2] -= tempz/tempm;
    }

  tempx=0.0;
  tempy=0.0;
  tempz=0.0;
  for(int a=0;a<nat;a++)
    {
      tempx+=Mass[a]*Vel[a][0];
      tempy+=Mass[a]*Vel[a][1];
      tempz+=Mass[a]*Vel[a][2];
      
      tempm+=Mass[a];
    }

  cout<<"v.c.o.m. = "<<tempx/tempm<<"   "<<tempy/tempm<<"   "<<tempz/tempm<<endl;
  
  //exit(1);
  ////end test/fix v.c.o.m. here


 
 ////define Spline parameters:
  fileread.close();
  fileread.clear();
  double smin=0.0;
  double smax=8.0;
  double sdelta=0.10;
  ifstream paramread(params_file);

  if ( paramread.fail() ) 
    {
      cout << "Could not open " << params_file << endl ;
      exit(1) ;
    }
    
  //spline parameters are in params.txt which is outputted from spline fitting and SVD programs.

  int tempint;

  // Read smin, smax, sdelta from params file instead of using built-in values.
  char buf[1024] ;
  do {
    paramread.getline(buf,1024) ;
  } while (buf[0] == '#' ) ;
  sscanf(buf,"%lf",&smin) ;
  
  int snum ;

  if ( pair_type == SPLINE ) 
    {
      paramread >> smax >> sdelta ;
      cout << "Spline minimum = " << smin << endl ;
      cout << "Spline maximum = " << smax << endl ;
      cout << "Spline step    = " << sdelta << endl ;
      snum=(2+floor((smax-smin)/sdelta))*2*3;//2 is for p0/m0/p1/m1, and 3 is for oo/oh/hh.
    }
  else if ( pair_type == CHEBYSHEV ) {
    paramread >> smax >> snum ;
    cout << "Cheby minimum = " << smin << endl ;
    cout << "Cheby maximum = " << smax << endl ;
    cout << "Cheby order    = " << snum << endl ;
    sdelta = 0.0 ;
    snum *= 3 ;
  } else {
    snum = 0 ;
    smax = Latcons[0] / 2.0 ;
    smin = 0.0 ;
  }
  if ( if_overcoord ) 
    {
      // Read overcoordination parameters from the params file.
      int n_over_read ;
      paramread >> n_over_read ;
      if ( n_over_read > 0 ) 
	{
	  printf("Reading overcoordination parameters from %s\n",
		 params_file) ;
	  n_over = n_over_read ;
	  for ( int k = 0 ; k < n_over ; k++ ) 
	    {
	      paramread >> over_param[k] ;
	      printf("Overcoordination parameter %d = %13.6e\n", k,
		     over_param[k]) ;
	    }
	}
    }

  double params[snum+4];
  double pot_params[snum/2+1] ;
  for(int i=0;i<snum+4;i++)
    params[i]=0.0;

  if ( pair_type == CHEBYSHEV or pair_type == SPLINE ) {
    for(int n=0;n<snum+4;n++)
      {
	paramread >> tempint >> params[n];
	if ( paramread.eof() ) 
	  {
	    cout << "Error reading params.txt\n" ;
	    exit(1) ;
	  }
	if ( tempint != n ) 
	  {
	    cout << "Error reading params.txt: index mismatch\n" ;
	    exit(1) ;
	  }
      }
  }

  if ( pair_type == SPLINE ) 
    {
      // Use a simple linear force model for positions not sampled.
      for ( int n = 0 ; n < 3 ; n++ ) {
	int ntable = snum / 3 ;     // Number of entries per atom pair.
	int nstart = n * ntable ;   // Start of current table.
	int nbegin = 0 ;            // Where populated table entries begin.

	for ( int idx = nstart ; idx < nstart + ntable ; idx += 2 ) {
	  if ( fabs(params[idx]) > 1.0 ) {
	    nbegin = idx ;
	    printf("Filling in table %d for values less than %d\n",
		   n, nbegin) ;
	    break ;
	  }
	}
	const double slope = -10.0 ;
	for ( int idx = nstart ; idx < nbegin ; idx += 2 ) {
	  params[idx] = slope * ( nbegin - idx) + params[nbegin] ;
	  params[idx+1] = slope / sdelta ;
	}    
      }
      for ( int j = 0 ; j < snum / 2 ; j++ ) {
	pot_params[j] = -11111 ;
      }

      // Calculate integral of the spline for the potential.
      for ( int j = snum / 6 - 1 ; j < snum / 2 ; j += snum / 6 ) {
	pot_params[j] = 0.0 ;
	for ( int n = j - 1 ; n >= j - snum / 6 + 1 ; n-- ) {
	  pot_params[n] = pot_params[n+1] - sdelta * ( params[2*n]   /2.0 + sdelta * params[2*n+1]/12.0
						       + params[2*n+2] /2.0 - sdelta * params[2*n+3]/12.0) ;
	}
      }
    } // if_cheby == FALSE .
    

  // Charges on H atom for potential models.
  double q_oo=0.0, q_oh=0.0, q_hh=0.0, q_spline = 0.0, q_stillinger=0.0;
  
  if ( if_spline_q ) {
    q_oo = params[snum]  / ke ;
    q_oh = params[snum+1] / ke ;
    q_hh = params[snum+2] / ke ;

    q_spline = sqrt(q_hh) ;
    printf("Q[H] from params file = %12.5f e\n", q_spline) ;
    printf("Charge equation 1 error = %12.5f\n", 
	   q_oo + 4.0 * q_oh + 4.0 * q_hh) ;

  } else {
    q_stillinger = sqrt(36.1345/ke) ;  // stillinger units
    printf("Built-in Q[H] = %12.5f e\n", q_stillinger) ;
  }



  double qsum = 0.0 ;
  for(int a=0;a<nat;a++)
    {
      if(Lb[a]=="O") 
	{
	  if ( if_spline_q ) {
	    Q[a] = -2.0 * q_spline ;
	  } else {
	    Q[a] = -2.0 * q_stillinger ;
	  }
	}
      else if(Lb[a]=="H") 
	{
	  if ( if_spline_q ) {
	    Q[a] = q_spline ;
	  } else {
	    Q[a] = q_stillinger ;
	  }
	}
      else
	{
	  cout<<"Atom type "<<Lb[a]<<" not supported."<<endl;
	  exit(1);
	}
      qsum += Q[a] ;
    }
  printf("Charge sum per atom = %13.5e\n", 
	 qsum / nat) ;


#if(0)
  // Output for testing the potential calculation.
  FILE *potout = fopen("potparams.txt", "w") ;
  for ( int j = 0 ; j < snum / 2 ; j++ ) {
    fprintf(potout, "%13.6e %13.6e\n", j*sdelta, pot_params[j]) ;
  }
  fclose(potout) ;

  FILE *fout = fopen("potvals.txt", "w") ;
  for ( int j = 0 ; j < 12 * snum / 3 ; j++ ) {
    double r = j * sdelta / 12.0  ;
    double S_r = 0.0 ;
    if ( r < smax ) {
      double spot = spline_pot(smin, smax, sdelta, r, params, pot_params, snum, 0, S_r) ;
      fprintf(fout, "%13.6e %13.6e\n", r, spot) ;
    } else {
      break ;
    }
  }
  fclose(fout) ;
#endif

  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      Accel[a1][c]=0;

  double dens_mol = (nat * 1.0e24) / (6.0221e23 * Vol) ;
  double dens_mass = tempm * dens_mol / nat ;

  printf("Volume                       = %8.5f Ang.^3\n", Vol) ;
  printf("Number density               = %8.5f mol atm/cc\n", dens_mol) ;
  printf("Mass density                 = %8.5f g/cc\n", dens_mass) ;


  double avg_temp = 0.0 ;
  double zeta_dot0 ;
  double temp_sum = 0.0 ;
  double press_sum = 0.0 ;

  for(int A=0;A<nsteps;A++)
    {//start Big Loop here.
      
      //update coordinates
      if(A>0)
	{
	  for(int a1=0;a1<nat;a1++)
	    for(int c=0;c<3;c++)
	      Coord[a1][c] += Vel[a1][c]*deltat+0.5*Accel[a1][c]*deltat*deltat;

	  if ( if_hoover ) 
	    {
	      for(int a1=0;a1<nat;a1++)
		for(int c=0;c<3;c++)
		  Coord[a1][c] -= 0.5 * zeta * Vel[a1][c] * deltat * deltat;
	      double ke = kinetic_energy(Mass, Vel, nat) ;
	      zeta_dot0 = ( 2.0 * ke - ndof * Kb * TempMD ) / ( QHoover) ;
	      s_hoover += zeta * deltat + 0.5 * deltat * deltat * zeta_dot0 ;
	    }
	

	  ////Wrap the coordinates:
	  for(int a1=0;a1<nat;a1++)
	    for(int c=0;c<3;c++)
	      {
		Coord[a1][c] -= floor(Coord[a1][c]/Latcons[c])*Latcons[c];
	      }
	  

	  //update first half of velocity:
	  for(int a1=0;a1<nat;a1++)
	    for(int c=0;c<3;c++)
	      Vel[a1][c]+=0.5*Accel[a1][c]*deltat;//updated with a(t)

	  if ( if_hoover ) 
	    {
	      for(int a1=0;a1<nat;a1++)
		for(int c=0;c<3;c++)
		  Vel[a1][c]-=0.5 * zeta * Vel[a1][c] * deltat;
	    }
	}      
      
      
      //Calculate the acceleration vector (in this case forces from FM):
      for(int a1=0;a1<nat;a1++)
	for(int c=0;c<3;c++)
	  Accel[a1][c]=0;

      //this function calculates the spline and electrostatic forces.
      ZCalc(Coord,Lb,Q,Latcons,nlayers,nat,smin,smax,sdelta,snum,params,pot_params,
	    pair_type,if_coulomb,if_overcoord,n_over,over_param,Accel,Vtot,Pxyz);

      if ( if_output_force ) 
	{
	  FILE *forcefile ;

	  forcefile = fopen("forceout.txt", "w") ;
	  for(int a1=0;a1<nat;a1++)
	    for(int c=0;c<3;c++)
	      fprintf(forcefile, "%13.6e\n", Accel[a1][c]) ;

	}

      if ( if_read_force ) 
	{
	  // Check against read-in forces for code verification.
	  double ferr = 0.0 ;
	  for(int a1=0;a1<nat;a1++)
	    for(int c=0;c<3;c++)
	      ferr += (Accel[a1][c]-Fread[a1][c]) *
		(Accel[a1][c] - Fread[a1][c]) ;
	  ferr /= 3.0 * nat ;
	  ferr = sqrt(ferr) ;
	  printf("RMS force error = %13.6e\n", ferr) ;
	}

      //Convert forces to acceleration:
      for(int a1=0;a1<nat;a1++)
	for(int c=0;c<3;c++)
	  Accel[a1][c] /= Mass[a1];

      if(A>0)
	{
	  //update second half of velocity
	  if ( if_hoover == false ) 
	    {
	      for(int a1=0;a1<nat;a1++)
		for(int c=0;c<3;c++)
		  Vel[a1][c]+=0.5*Accel[a1][c]*deltat;//update velocity with a(t+dt)
	    } 
	  else 
	    {
	      double ke = kinetic_energy(Mass, Vel, nat) ;
	      double zeta_dot1 = ( 2.0 * ke - ndof * Kb * TempMD ) / ( QHoover) ;
	      double zeta1 = 0.0 ;
	      // Iterative determination of velocities
	      // See Martyna, Tobias, Klein JCP 101, 4177(1994) Appendix D.
	      for ( int itr = 0 ; itr < 10 ; itr++ ) {
		zeta1 = zeta + (zeta_dot0 + zeta_dot1) * 0.5 * deltat ;
		double vscaleh = 1.0 + 0.5 * deltat * zeta1 ;
		for(int a1=0;a1<nat;a1++)
		  for(int c=0;c<3;c++)
		    Velnew[a1][c] = (Vel[a1][c] + 0.5*Accel[a1][c]*deltat) / vscaleh ;
		ke = kinetic_energy(Mass, Velnew, nat) ;
		zeta_dot1 = ( 2.0 * ke - ndof * Kb * TempMD ) / ( QHoover) ;
	      }
	      for(int a1=0;a1<nat;a1++)
		for(int c=0;c<3;c++)
		  Vel[a1][c] = Velnew[a1][c] ;
	      zeta = zeta1 ;
	    }
      	}


      //calculate kinetic energy for scaling:
      Ktot=kinetic_energy(Mass, Vel, nat) ;

      // Subtract 3.0 to account for zeroing the center of mass velocity.
      double temperature = 2.0 * Ktot / (ndof * Kb) ;
      temp_sum += temperature ;

      if ( num_pressure ) 
	{
	  Pxyz = numerical_pressure(Coord,Lb, Q, Latcons,nlayers, nat,smin,
				    smax, sdelta,snum, params, pot_params, 
				    pair_type,if_coulomb,if_overcoord,
				    n_over, over_param) ;
	}
      double Ptot = Pxyz + 2.0 * Ktot / (3.0 * Vol) ;
      // Unit conversion factor to GPa.
      Ptot *= GPa ;

      press_sum += Ptot ;

      if ( A == 0 ) 
	{
	  printf("%8s %9s %15s %15s %15s %15s %15s",
		 "Step", "Time", "Ktot/N", "Vtot/N", "Etot/N", 
		 "T", "P") ;
	  if ( if_hoover ) {
	    printf(" %15s\n", "Econs/N") ;
	  } else {
	    printf("\n") ;
	  }
	  printf("%8s %9s %15s %15s %15s %15s %15s",
		 " ", "(fs)", "(kcal/mol)", "(kcal/mol)", "(kcal/mol)", 
		 "(K)", "(GPa)") ;
	  if ( if_hoover ) {
	    printf(" %15s\n", "(kcal/mol)") ;
	  } else {
	    printf("\n") ;
	  }
	}
      if ( (A +1) % energy_freq == 0 ) {
	printf("%8d %9.2f %15.7f %15.7f %15.7f %15.1f %15.3f",
	       A+1, (A+1)*deltat_fs, Ktot/nat,Vtot/nat,(Ktot+Vtot)/nat, 
	       temperature, Ptot) ;
	if ( if_hoover ) {
	  double E_hoover = Ktot + Vtot + 0.5 * zeta * zeta * QHoover +
	    ndof * Kb * TempMD * s_hoover ;
	  printf("%15.7f\n", E_hoover/nat) ;
	} else {
	  printf("\n") ;
	}
      } 

      if ( scale_freq > 0 && (A+1) % scale_freq == 0 ) 
	// Scale velocities.
	{
	  avg_temp += temperature ;
	  avg_temp /= scale_freq ;

	  cout << "Average temperature = " << avg_temp << endl ;

	  double vscale = sqrt(TempMD/avg_temp) ;
	  avg_temp = 0.0 ;

	  cout << "Velocity scaling factor = " << vscale << endl ;
	  for(int a1=0;a1<nat;a1++)
	    for(int c=0;c<3;c++)
	      Vel[a1][c] *= vscale ;
	}
      else {
	avg_temp += temperature ;
      }


      if ( output_gen && (A+1) % gen_freq == 0 ) {
	// Write coordinates to gen file.
	fprintf(fgen, "%5d S #Step %d Time %7.4f (fs) Temp %7.4f (K)\n", nat, A+1, (A+1) * deltat_fs,
		temperature) ;
	// Special to H2O:
	fprintf(fgen, "O H\n") ;
	
	
	for ( int a1 = 0 ; a1 < nat ; a1++ ) 
	  {
	    int iele ;
	    if ( Lb[a1] == "O" ) 
	      {
		iele = 1 ;
	      } else if ( Lb[a1] == "H" ) 
	      {
		iele = 2 ;
	      } else 
	      {
		cout << "Error: bad element label\n" ;
		exit(1) ;
	      }

	    fprintf(fgen, "%4d %2d %8.5f %8.5f %8.5f\n", 
		    a1+1, iele, Coord[a1][0], Coord[a1][1], Coord[a1][2]) ;
	  }
	fprintf(fgen, "%8.5f %8.5f %8.5f\n", 0.0, 0.0, 0.0) ;
	fprintf(fgen, "%8.5f %8.5f %8.5f\n", Latcons[0], 0.0, 0.0) ;
	fprintf(fgen, "%8.5f %8.5f %8.5f\n",        0.0, Latcons[1], 0.0) ;
	fprintf(fgen, "%8.5f %8.5f %8.5f\n",        0.0, 0.0, Latcons[2]) ;
      }
      
    }//End big loop here.

  tempm = 0.0 ;
  for ( int ia = 0 ; ia < nat ; ia++ ) 
    {
      tempm += Mass[ia] ;
    }

  printf("Average temperature over run = %8.2f K\n", temp_sum / nsteps) ;
  printf("Average pressure over run    = %8.4f GPa\n", press_sum / nsteps) ;

  // Output final xyz position in the same format as input.xyz for restarting.
  FILE *fxyz ;
  fxyz = fopen("output.xyz", "w") ;
  if ( fxyz == NULL ) 
    {
      cout << "Error: could not open output.xyz\n" ;
      exit(1) ;
    }
  fprintf(fxyz, "%5d\n", nat) ;
  fprintf(fxyz, "%8.5f %8.5f %8.5f\n", Latcons[0], Latcons[1], Latcons[2]) ;
  for ( int ia = 0 ; ia < nat ; ia++ ) {
    fprintf(fxyz, "%2s %8.5f %8.5f %8.5f    %8.5f %8.5f %8.5f\n",
	    Lb[ia].c_str(), Coord[ia][0], Coord[ia][1], Coord[ia][2],
	    Vel[ia][0], Vel[ia][1], Vel[ia][2]) ;
  }

  //cout<<endl<<endl;
  return 0;    
}       


static void read_input(const char *file_name,
		       double &TempMD, double &deltat_fs, int &nsteps, 
		       int &nlayers, 
		       bool &if_output_force, bool &if_read_force,
		       bool &if_spline_q, bool &if_init_vel, int &rand_seed, 
		       int &gen_freq, int &energy_freq, int &scale_freq,
		       double &thoover_fs, Sr_pair_t &pair_type,
		       bool &if_coulomb, bool &num_pressure, char *params_file,
		       bool &if_overcoord) 
// Read program input from the file "splines_md.in".
{
  FILE *fin ;
  const int bufsz = 1024 ;
  char buf[bufsz] ;

  fin = fopen(file_name, "r") ;
  if ( fin == NULL ) 
    {
      cout << "Error: could not open " << file_name << endl ;
      exit(1) ;
    }
  
  while ( char *line = fgets(buf, bufsz, fin ) ) 
    {
      if ( line[0] == '#' ) {
	continue ;
      }
      // Get rid of new line character.
      // int iend = strlen(line) - 1 ;
      // if ( iend > 0 && line[iend] == '\n' ) {
      // line[iend] = 0 ;
      // }

      // printf("Line = %s\n", line) ;
      char *name = strtok(line, " ") ;
      char *val = strtok(NULL, " ") ;

      // Get rid of new line character.
      int iend = strlen(val) - 1 ;
      if ( iend > 0 && val[iend] == '\n' ) {
	val[iend] = 0 ;
      }
      if ( strncmp(name,"temperature",bufsz) == 0 ) 
	{
	  sscanf(val, "%lf", &TempMD) ;
	  //	  printf("Temperature = %7.4f\n", TempMD) ;	
	}
      else if ( strncmp(name,"deltat",bufsz) == 0 ) 
	{
	  sscanf(val, "%lf", &deltat_fs) ;
	  //	  printf("Time step = %7.4f\n", deltat_fs) ;	
	}
      else if ( strncmp(name,"nsteps",bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &nsteps) ;
	  //	  printf("Number of steps = %d\n", nsteps) ;	
	}
      else if ( strncmp(name,"nlayers",bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &nlayers) ;
	  //	  printf("Number of layers = %d\n", nlayers) ;
	}
      else if ( strncmp(name,"output_force",bufsz) == 0 ) 
	{
	  if_output_force = parse_tf(val, bufsz,line) ;
	  //	  printf("Output_force is %s\n", 
	  //		 if_output_force?"true":"false") ;
	}
      else if ( strncmp(name,"coulomb",bufsz) == 0 ) 
	{
	  if_coulomb = parse_tf(val, bufsz,line) ;
	}
      else if ( strncmp(name,"overcoord",bufsz) == 0 ) 
	{
	  if_overcoord = parse_tf(val, bufsz,line) ;
	}
      else if ( strncmp(name,"read_force",bufsz) == 0 ) 
	{
	  if_read_force = parse_tf(val, bufsz,line) ;
	  //	  printf("Read_force is %s\n", 
	  //		 if_read_force?"true":"false") ;
	}
      else if ( strncmp(name,"num_pressure",bufsz) == 0 ) 
	{
	  num_pressure = parse_tf(val, bufsz,line) ;
	}
      else if ( strncmp(name,"fit_coulomb",bufsz) == 0 ) 
	{
	  if_spline_q = parse_tf(val, bufsz,line) ;
	  //	  printf("spline_q is %s\n", 
	  //		 if_spline_q?"true":"false") ;
	}
      else if ( strncmp(name,"init_vel",bufsz) == 0 ) 
	{
	  if_init_vel = parse_tf(val, bufsz,line) ;
	  //	  printf("init_vel is %s\n", 
	  //		 if_init_vel?"true":"false") ;
	}
      else if ( strncmp(name,"rand_seed",bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &rand_seed) ;
	  //	  printf("Random number seed = %d\n", rand_seed) ;	
	}
      else if ( strncmp(name,"gen_freq",bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &gen_freq) ;
	}
      else if ( strncmp(name,"energy_freq",bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &energy_freq) ;
	  //	  printf("Energy print frequency = %d\n", energy_freq) ;	
	}
      else if ( strncmp(name,"scale_freq",bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &scale_freq) ;
	}
      else if ( strncmp(name,"hoover_time",bufsz) == 0 ) 
	{
	  sscanf(val, "%lf", &thoover_fs) ;
	}
      else if ( strncmp(name,"params_file",bufsz) == 0 ) 
	{
	  strcpy(params_file,val) ;
	}
      else if ( strncmp(name, "pair_type", bufsz) == 0 ) 
	{
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
	}
      else 
	{
	  printf("Error: did not recognize option %s\n", name) ;
	  exit(1) ;
	}
    }

  if ( ferror(fin) ) 
    {
      printf("Error while reading spline_md.in\n") ;
      exit(1) ;
    }
}

static void echo_input(double TempMD, double deltat_fs, int nsteps, 
		       int nlayers, 
		       bool if_output_force, bool if_read_force,
		       bool if_spline_q, bool if_init_vel, int rand_seed, 
		       int gen_freq, int energy_freq, int scale_freq,
		       double thoover_fs, Sr_pair_t pair_type,
		       bool if_coulomb,bool num_pressure,
		       const char* params_file, bool if_overcoord) 
{
  printf("JOB PARAMETERS:\n\n") ;

  printf("Temperature is %8.5f\n", TempMD) ;
  printf("Time step is %8.5f fs\n", deltat_fs) ;
  printf("Number of steps is %d\n", nsteps) ;
  printf("Number of layers is %d\n", nlayers) ;

  if ( pair_type == CHEBYSHEV ) {
      printf("Chebyshev short-range pair forces will be used\n") ;
  } 
  else if ( pair_type == SPLINE ) {
    printf("Spline short-range pair forces will be used\n") ;
  }
  else if ( pair_type == INVERSE_R ) {
    printf("Inverse power short-range pair forces will be used\n") ;
  }
  else if ( pair_type == LJ ) {
    printf("Lennard-Jones short-range pair forces will be used\n") ;
  } 
  else if ( pair_type == STILLINGER ) {
    printf("Stillinger short-range pair forces will be used\n") ;
  } else {
    printf("Error: unknown short-range pair interactions\n") ;
    exit(1) ;
  }

  printf("Force parameters will be read from %s\n", params_file) ;

  if ( if_coulomb ) {
    printf("Coulomb forces will be calculated with Ewald sums\n") ;
  }
  else {
    printf("Coulomb forces will not be calculated\n") ;
  }

  if ( if_overcoord ) {
    printf("ReaxFF overcoordination term will be calculated\n") ;
  } else {
    printf("ReaxFF overcoordination term will NOT be calculated\n") ;
  }

  if ( if_output_force ) {
    printf("Forces will be output (testing)\n") ;
  }     

  if ( if_read_force ) {
    printf("Forces will be read for testing\n") ;
  }

  if ( if_spline_q ) 
    {
      printf("Charges will be taken from params.txt\n") ;
    }
  else
    {
      printf("Built-in charges will be used\n") ;
    }

  if ( num_pressure ) 
    {
      printf("Pressures will be calculated by finite difference of the energy\n") ;
    }
  else
    {
      printf("Pressures will be calculated using the virial\n") ;
    }

  if ( if_init_vel ) 
    {
      printf("Velocities will be initialized\n") ;
    }
  else 
    {
      printf("Velocities will be read from input.xyz\n") ;
    }

    
  printf("Random number seed = %d\n", rand_seed) ;

  if ( gen_freq > 0 ) 
    {
      printf("Trajectory will be written to traj.gen every %d steps\n",
	     gen_freq) ;
    }
  
  printf("Energies will be written every %d steps\n", energy_freq) ;

  if ( scale_freq > 0 ) 
    {
      printf("Velocities will be scaled every %d steps\n",
	     scale_freq) ;
    }

  if ( thoover_fs > 0.0 ) 
    {
      printf("Nose-Hoover thermostat time = %8.5f fs\n", thoover_fs) ;
    }
    
  printf("\nEND OF JOB PARAMETERS\n\n") ;

}

  
static double kinetic_energy(double *Mass, double **Vel, int nat) 
{
  double Ktot = 0.0 ;
  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      Ktot+=0.5*Mass[a1]*Vel[a1][c]*Vel[a1][c];

  return(Ktot) ;
}

static double 
numerical_pressure(double **Coord, string *Lb, double *Q, double *Latcons,
		   const int nlayers, const int nat,const double smin,
		   const double smax, const double sdelta,const int snum, 
		   double *params, double *pot_params, Sr_pair_t pair_type,
		   bool if_coulomb,bool if_overcoord, int n_over,
		   double *over_param) 
// Evaluates the configurational part of the pressure numerically by -dU/dV.
{
  double **Coord1 ;
  double **Accel ;
  const double eps = 1.0e-04 ;
  double lscale ;
  double Latcons1[3] ;
  double Vtot1, Vtot2 ;
  double Vol1, Vol2 ;
  double Pxyz ;

  Coord1= new double* [nat];
  Accel = new double* [nat] ;
  for ( int j = 0 ; j < nat ; j++ ) {
    Coord1[j] = new double [3] ;
    Accel[j]  = new double [3] ;
  }

  lscale = 1.0  + eps ;
  for ( int j = 0 ; j < nat ; j++ ) 
    {
      Coord1[j][0] = lscale * Coord[j][0] ;
      Coord1[j][1] = lscale * Coord[j][1] ;
      Coord1[j][2] = lscale * Coord[j][2] ;
    }

  Latcons1[0] = Latcons[0] * lscale ;
  Latcons1[1] = Latcons[1] * lscale ;
  Latcons1[2] = Latcons[2] * lscale ;

  Vol1 = Latcons1[0] * Latcons1[1] * Latcons1[2] ;

  ZCalc(Coord1,Lb,Q,Latcons1,nlayers,nat,smin,smax,sdelta,snum,params,
	pot_params,pair_type,if_coulomb,if_overcoord,n_over,over_param,
	Accel,Vtot1,Pxyz);
  
  lscale = 1.0 - eps ;

  for ( int j = 0 ; j < nat ; j++ ) 
    {
      Coord1[j][0] = lscale * Coord[j][0] ;
      Coord1[j][1] = lscale * Coord[j][1] ;
      Coord1[j][2] = lscale * Coord[j][2] ;
    }

  Latcons1[0] = Latcons[0] * lscale ;
  Latcons1[1] = Latcons[1] * lscale ;
  Latcons1[2] = Latcons[2] * lscale ;

  Vol2 = Latcons1[0] * Latcons1[1] * Latcons1[2] ;

  ZCalc(Coord1,Lb,Q,Latcons1,nlayers,nat,smin,smax,sdelta,snum,params,
	pot_params,pair_type,if_coulomb,if_overcoord,n_over,over_param,
	Accel,Vtot2,Pxyz);

  double result = -(Vtot2 - Vtot1)/(Vol2 - Vol1) ;

  for ( int j = 0 ; j < nat ; j++ ) {
    delete[] Coord1[j] ;
    delete[] Accel[j]  ;
  }
  delete[] Coord1 ;
  delete[] Accel ;

  // cout << "Numerical pressure = " << result << endl ;

  return(result) ;
}
