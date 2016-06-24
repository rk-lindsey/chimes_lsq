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
		       bool &if_init_vel, int &rand_seed, 
		       int &gen_freq, int &energy_freq, int &scale_freq,
		       double &thoover_fs, 
		       bool &num_pressure, char *params_file, char *xyz_file) ;

static void echo_input(double TempMD, double deltat_fs, int nsteps, 
		       int nlayers, 
		       bool if_output_force, bool if_read_force,
		       bool fit_coul, bool if_init_vel, int rand_seed, 
		       int gen_freq, int energy_freq, int scale_freq,
		       double thoover_fs, Sr_pair_t pair_type,
		       bool if_coulomb, bool num_pressure, const char *params_file,
		       bool if_overcoord, bool if_3b_cheby, const char *xyz_file) ; 


static double kinetic_energy(double *Mass, double **Vel, int nat) ;

static double 
numerical_pressure(double **Coord, const char *Lbc, double *Q, double *Latcons,
		   const int nlayers, const int nat,const double *smin,
		   const double *smax, const double *sdelta,const int *snum, 
		   const int *snum_3b_cheby,
		   double *params, double *pot_params, Sr_pair_t pair_type,
		   bool if_coulomb,bool if_overcoord, bool if_3b_cheby,
		   int n_over,
		   double *over_param, const double *lambda)  ;

static double *Parse_Params_File(char *params_file, double *smin, double *smax, double *sdelta,
				 double *lambda, int *snum, int *snum_3b_cheby, int &n_over,
				 double *over_param, const double *Latcons,
				 int &tot_snum, bool &if_3b_cheby, bool &if_coulomb, bool &fit_coul,
				 int &num_cheby_3b, bool &if_overcoord, bool &fit_pover,
				 Sr_pair_t &pair_type) ;

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
  bool fit_coul = true ;

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

  // If true, find linear overcoordination parameter from least-squares fitting.
  bool fit_pover = true ;

  // If true, calculate 3-Body Chebyshev interaction.
  bool if_3b_cheby = false ;

  // How often to output energy
  int energy_freq = 10 ;

  // How often to scale velocities (0 if never).
  int scale_freq = 0 ;  

  // How often to write the gen file.
  int gen_freq = 10 ;

  // Whether to calculate pressures by finite difference.
  bool num_pressure = false ;

  char params_file[1024] = "params.txt" ;

  char xyz_file[1024] = "input.xyz" ;

  double over_param[MAXOVERP] =   {
    // Lucas's parameters for overcoordination
     75.0,                    // pover
     9.6000387075695e-01,     // r0
     -2.5881042987450e-01,    // p1
     3.8995379237250e+00,     // p2
     -8.9                     // lambda6
    //
    // Chenoweth values
    //50.0,      // pover
    //1.0165,    // r0
    //-0.0657,   // p1
    //5.0451,    // p2
    //-3.6141    // lambda6
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
	     if_init_vel, rand_seed, gen_freq, energy_freq, scale_freq,
	     thoover_fs, num_pressure,params_file,xyz_file) ;

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
      cout << "Reading positions from input.xyzf for force testing\n" ;
      fileread.open("input.xyzf")  ;
      //      if ( fileread.error() ) {
      //cout << "Could not open input.xyzf\n" << endl ;
      //exit(1) ;
      // }
    } 
  else 
    {
      cout << "Reading positions from " << xyz_file << endl  ;
      fileread.open(xyz_file)  ;
      //if ( fileread.error() ) {
      //cout << "Could not open input.xyz\n" << endl ;
      //exit(1) ;
      //}
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

  ifstream forceread ;
  if ( if_read_force ) {
    cout << "Opening force.txt to read forces for comparison\n" ;
    forceread.open("force.txt") ;
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
	      double junk ;
	      // Read past forces in the input.xyzf file
	      fileread >> junk ;

	      // Read forces from force.txt.
	      forceread >>Fread[a][c];
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
	  cout<<"Atom type "<<Lb[a]<< a << " not supported."<<endl;
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
  tempm = 0.0;
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
  double smin[NPAIR];
  double smax[NPAIR];
  double sdelta[NPAIR];
  double lambda[NPAIR];
  int snum[NPAIR] ;
  int snum_3b_cheby[NPAIR] ;
  int tot_snum = 0;
  int num_cheby_3b = 0 ;
  double *params ;

  params = Parse_Params_File(params_file, smin, smax, sdelta, lambda, snum, snum_3b_cheby, n_over,
			     over_param, Latcons, tot_snum, if_3b_cheby, if_coulomb, fit_coul,
			     num_cheby_3b, if_overcoord, fit_pover, pair_type) ;

  echo_input(TempMD, deltat_fs, nsteps, nlayers, if_output_force, if_read_force,
	     fit_coul, if_init_vel, rand_seed, gen_freq, energy_freq, scale_freq,
	     thoover_fs,pair_type,if_coulomb,num_pressure,params_file,
	     if_overcoord, if_3b_cheby,xyz_file) ;


  double pot_params[tot_snum/2+1] ;

  if ( pair_type == SPLINE ) 
    {
      // Use a simple linear force model for positions not sampled.
      int nstart = 0;
      for ( int n = 0 ; n < NPAIR ; n++ ) {
	int ntable = snum[n];     // Number of entries per atom pair.
        if (n > 0) {
	  nstart += snum[n-1];   // Start of current table.
        }
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
	  params[idx+1] = slope / sdelta[n] ;
	}    
	for ( int j = 0 ; j < snum[n] / 2 ; j++ ) {
	  pot_params[j+nstart/2] = -11111 ;
	}

	// Calculate integral of the spline for the potential.
	int j = nstart/2+snum[n]/2-1 ;
	pot_params[j] = 0.0 ;
	for ( int m = j - 1 ; m >= nstart/2 ; m-- ) {
	  pot_params[m] = pot_params[m+1] - sdelta[n] * ( params[2*m]   /2.0 + sdelta[n]* params[2*m+1]/12.0
							  + params[2*m+2] /2.0 - sdelta[n] * params[2*m+3]/12.0) ;
	}
      } // n < NPAIR

    } // if_cheby == FALSE .
    


  if ( if_coulomb ) 
    {
      // Charges on H atom for potential models.
      double q_oo=0.0, q_oh=0.0, q_hh=0.0, q_spline = 0.0, q_stillinger=0.0;
  
      if ( fit_coul ) {
	q_oo = params[tot_snum+num_cheby_3b]  / ke ;
	q_oh = params[tot_snum+num_cheby_3b+1] / ke ;
	q_hh = params[tot_snum+num_cheby_3b+2] / ke ;

	q_spline = sqrt(q_hh) ;
	printf("Q[H] from params file = %12.5f e\n", q_spline) ;
	printf("Charge equation 1 error = %12.5f\n", 
	       q_oo + 4.0 * q_oh + 4.0 * q_hh) ;
      } else {
	q_stillinger = sqrt(44.23918853232863/ke) ;  // spce charge in stillinger units
	printf("Built-in Q[H-DFT] = %12.5f e\n", q_stillinger) ;
      }

      double qsum = 0.0 ;
      for(int a=0;a<nat;a++)
	{
	  if(Lb[a]=="O") 
	    {
	      if ( fit_coul ) {
		Q[a] = -2.0 * q_spline ;
	      } else {
		Q[a] = -2.0 * q_stillinger ;
	      }
	    }
	  else if(Lb[a]=="H") 
	    {
	      if ( fit_coul ) {
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
    }

  if ( pair_type == SPLINE ) 
    {
      // Output for testing the spline potential calculation.
      // potparams.txt contains the modified spline potential as calculated from forces.
      FILE *potout = fopen("potparams.txt", "w") ;
      int jstart = 0 ;
      for (int i = 0; i < NPAIR; i++) {
	for ( int j = 0 ; j < snum[i] / 2 ; j++ ) {
	  fprintf(potout, "%13.6e %13.6e\n", smin[i]+j*sdelta[i], pot_params[j+jstart]) ;
	}
	jstart += snum[i] / 2 ;
      }
      fclose(potout) ;

      // potvals.txt contains the potential as evaluated with the spline_pot function.
      FILE *fout = fopen("potvals.txt", "w") ;
      jstart = 0 ;
      for (int i = 0; i < NPAIR; i++) {
	for ( int j = 0 ; j < snum[i] ; j++ ) {
	  // r chosen to be "off" the spline array values by starting at 0.0, and using 1/2 steps.
	  double r = j * sdelta[i] / 2.0 ;
	  double S_r = 0.0 ;
	  double spot = spline_pot(smin[i], smax[i], sdelta[i], r, params, pot_params, snum[i], jstart, S_r) ;
	  fprintf(fout, "%13.6e %13.6e\n", r, spot) ;
	}
	jstart += snum[i] ;
      }
      fclose(fout) ;
    }

  for(int a1=0;a1<nat;a1++)
    for(int c=0;c<3;c++)
      Accel[a1][c]=0;

  double dens_mol = (nat * 1.0e24) / (6.0221e23 * Vol) ;
  double dens_mass = (tempm/Vol)*(1e24/6.0221e23) ;
  printf("Total Mass                   = %8.5f au\n",tempm);
  printf("Volume                       = %8.5f Ang.^3\n", Vol) ;
  printf("Number density               = %8.5f mol atm/cc\n", dens_mol) ;
  printf("Mass density                 = %8.5f g/cc\n", dens_mass) ;


  double avg_temp = 0.0 ;
  double zeta_dot0 = 0.0 ;
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
      char *Lbc ;

      Lbc = new char [nat] ;
      for(int a=0;a<nat;a++) {
	if ( Lb[a] == "O" ) {
	  Lbc[a] = 'O' ;
	} else if ( Lb[a] == "H" ) {
	  Lbc[a] = 'H' ;
	} else {
	  cout << "Error: unknown element " << Lbc[a] << "\n" ;
	}
      }

      ZCalc(Coord,Lbc,Q,Latcons,nlayers,nat,smin,smax,sdelta,snum,
	    snum_3b_cheby, params,pot_params,
	    pair_type,if_coulomb,if_overcoord,if_3b_cheby,
	    n_over,over_param,lambda,Accel,Vtot,Pxyz);

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
	  exit(0) ;
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
	  Pxyz = numerical_pressure(Coord,Lbc, Q, Latcons,nlayers, nat,smin,
				    smax, sdelta,snum, 
				    snum_3b_cheby,
				    params, pot_params, 
				    pair_type,if_coulomb,if_overcoord, if_3b_cheby,
				    n_over, over_param, lambda) ;
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
	  std::cout.flush() ;
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
	  double E_tot = Ktot + Vtot ;
	  printf("%15.7f\n", E_tot/nat) ;
	  printf("\n") ;
	}
	std::cout.flush() ;
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
	fflush(fgen) ;
      }

      delete[] Lbc ;
      
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
		       bool &if_init_vel, int &rand_seed, 
		       int &gen_freq, int &energy_freq, int &scale_freq,
		       double &thoover_fs, 
		       bool &num_pressure, char *params_file, char *xyz_file) 
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
      else if ( strncmp(name,"xyz_file",bufsz) == 0 ) 
	{
	  strcpy(xyz_file,val) ;
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
		       bool fit_coul, bool if_init_vel, int rand_seed, 
		       int gen_freq, int energy_freq, int scale_freq,
		       double thoover_fs, Sr_pair_t pair_type,
		       bool if_coulomb,bool num_pressure,
		       const char* params_file, bool if_overcoord, bool if_3b_cheby,
		       const char *xyz_file) 
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
  printf("XYZ coordinates will be read from %s\n", xyz_file) ;

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

  if ( if_3b_cheby ) {
    printf("3-Body Chebyshev polynomial will be calculated\n") ;
  } else {
    printf("3-Body Chebyshevy polynomial will NOT be calculated\n") ;
  }

  if ( if_output_force ) {
    printf("Forces will be output (testing)\n") ;
  }     

  if ( if_read_force ) {
    printf("Forces will be read for testing\n") ;
  }

  if ( fit_coul ) 
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
numerical_pressure(double **Coord, const char *Lbc, double *Q, double *Latcons,
		   const int nlayers, const int nat,const double *smin,
		   const double *smax, const double *sdelta,const int *snum, 
		   const int *snum_3b_cheby,
		   double *params, double *pot_params, Sr_pair_t pair_type,
		   bool if_coulomb,bool if_overcoord, bool if_3b_cheby,
		   int n_over,
		   double *over_param, const double *lambda) 
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

  ZCalc(Coord1,Lbc,Q,Latcons1,nlayers,nat,smin,smax,sdelta,snum,
	snum_3b_cheby, params,
	pot_params,pair_type,if_coulomb,if_overcoord,if_3b_cheby,
	n_over,over_param,
	lambda,Accel,Vtot1,Pxyz);
  
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

  ZCalc(Coord1,Lbc,Q,Latcons1,nlayers,nat,smin,smax,sdelta,snum,
	snum_3b_cheby, params,
	pot_params,pair_type,if_coulomb,if_overcoord,if_3b_cheby,
	n_over,over_param,
	lambda,Accel,Vtot2,Pxyz);

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

static double *Parse_Params_File(char *params_file, double *smin, double *smax, double *sdelta,
				 double *lambda, int *snum, int *snum_3b_cheby, int &n_over,
				 double *over_param, const double *Latcons,
				 int &tot_snum, bool &if_3b_cheby, bool &if_coulomb, bool &fit_coul,
				 int &num_cheby_3b, bool &if_overcoord, bool &fit_pover,
				 Sr_pair_t &pair_type)
// Parse a file with potential interaction parameters and settings.
{
  int i;

  ifstream paramread(params_file);

  if ( paramread.fail() ) 
    {
      cout << "Could not open " << params_file << endl ;
      exit(1) ;
    }
    
  //spline parameters are in params.txt which is outputted from spline fitting and SVD programs.

  int tempint;

  // Read smin, smax, sdelta from params file instead of using built-in values.
  const int bufsz = 1024 ;
  char buf[bufsz] ;
  char cmd_arg[bufsz] ;

  do {
    paramread.getline(buf,bufsz) ;
  } while (buf[0] == '#' ) ;

  int npair_read = 0 ;
  sscanf(buf,"npair %d",&npair_read) ;

  if ( npair_read != NPAIR ) {
    cout << "Inconsistent npair found:  compiled npair = " << NPAIR << endl ;
  }
  if_3b_cheby = read_tf_option(&paramread, "3b_cheby", params_file) ;

  paramread.getline(buf, bufsz) ;
  if ( ! strncmp(buf, "pair_type", 9) ) {
    sscanf(buf, "pair_type %s", cmd_arg) ;
    pair_type = parse_pair_type(cmd_arg, bufsz) ;
  } else {
    printf("Error in reading %s: pair type not found\n", params_file) ;
  }

  for (i = 0; i < NPAIR; i++) {
    if ( pair_type == SPLINE ) 
      {
        paramread >> smin[i] >> smax[i] >> sdelta[i] >> snum[i] ;
        cout << "Spline minimum = " << smin[i] << endl ;
        cout << "Spline maximum = " << smax[i] << endl ;
        cout << "Spline step    = " << sdelta[i] << endl ;
	cout << "Number of spline params = " << snum[i] << endl ;
        // snum[i]=(2+floor((smax[i]-smin[i])/sdelta[i]))*2;//2 is for p0/m0/p1/m1
      }
    else if ( pair_type == CHEBYSHEV ) 
      {
	cout << "pair type: " << i << endl;
	paramread >> smin[i] >> smax[i] >> lambda[i] >> snum[i];
	cout << "Cheby minimum = " << smin[i] << endl ;
	cout << "Cheby maximum = " << smax[i]  << endl ;
	cout << "Cheby order    = " << snum[i]  << endl ;
	cout << "Morse lambda value= " << lambda[i]  << endl ;
	if ( if_3b_cheby ) {
	  paramread >> snum_3b_cheby[i] ;
	  cout << "3 Body Cheby order = "<< snum_3b_cheby[i] << endl ;
	}
	sdelta[i] = 0.0 ;
      }
    else if (pair_type == INVERSE_R) 
      {
	cout << "pair type: " << i << endl;
	paramread >> smin[i] >> smax[i] >> snum[i];
	cout << "Inverse R minimum = " << smin[i] << endl ;
	cout << "Inverse R maximum = " << smax[i]  << endl ;
	cout << "Inverse R order    = " << snum[i]  << endl ;
      } else 
      {
	snum[i] = 0 ;
	smax[i] = Latcons[0] / 2.0 ;
	smin[i] = 0.0 ;
      }
    tot_snum += snum[i];
  }
  cout << "Total number of 2-Body potential parameters: " << tot_snum << endl;

  paramread.getline(buf, bufsz) ;

  if_coulomb = read_tf_option(&paramread, "coulomb", params_file) ;
  fit_coul   = read_tf_option(&paramread, "fit_coulomb", params_file) ;

  if_overcoord = read_tf_option(&paramread, "overcoord", params_file) ;

  if ( if_overcoord ) 
    {
      // Read overcoordination parameters from the params file.
      int n_over_read = 0 ;

      paramread.getline(buf, bufsz) ;
      if ( ! strncmp(buf, "nover", 5) ) {
	sscanf(buf, "nover %d", &n_over_read) ;
      }

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
	  // Read final newline.
	  paramread.getline(buf, bufsz) ;
	}
    }
  fit_pover = read_tf_option(&paramread, "fit_pover", params_file) ;

  num_cheby_3b = 0 ;
  if ( if_3b_cheby ) {
    num_cheby_3b = count_cheby_3b_params(snum_3b_cheby) ;
    cout << "Number of 3-Body Chebyshev parameter =" << num_cheby_3b << endl ;
  }

  paramread.getline(buf, bufsz) ;
  if ( strncmp(buf, "least squares", 12) ) {
    printf("Error: least squares parameters not found\n") ;
    printf("Current line = |%s|\n", buf) ;
    exit(1) ;
  }
  // Number of coulomb parameters.
  int ncoulomb ;
  
  // Number of linear fit overcoordination parameters.
  int noverlin ;

  if ( fit_coul == true ) {
    ncoulomb = NPAIR ;
  } else {
    ncoulomb = 0 ;
  }
  if ( fit_pover ) {
    noverlin = 1 ;
  } else {
    noverlin = 0 ;
  }
  double *params = new double[tot_snum+num_cheby_3b+ncoulomb+noverlin];
  for(int i=0;i<tot_snum+num_cheby_3b+ncoulomb+noverlin;i++)
    params[i]=0.0;

  if ( pair_type == CHEBYSHEV or pair_type == SPLINE or pair_type == INVERSE_R ) {
    cout << "Potential parameters read in:\n" ;
    cout << "Number of 2 body parameters = " << tot_snum << endl ;
    cout << "Number of chebyshev 3 body parameters = " << num_cheby_3b << endl ;
    cout << "Fit coulomb parameters = " << ncoulomb << endl ;
    cout << "Fit overcoordination parameters = " << noverlin << endl ;

    for(int n=0;n<tot_snum+num_cheby_3b+ncoulomb+noverlin;n++)
      {
	paramread >> tempint >> params[n];

	if ( tempint != n ) {
	  cout << "Error: parameter index mismatch " << tempint <<  " " 
	       << params[n] << " " << n << endl;
	  exit(1) ;
	}

        cout << tempint << " " << params[n] << endl;
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
  if ( fit_pover ) {
    over_param[0] = params[tot_snum+num_cheby_3b+ncoulomb] ;
    cout << "Linear fit overcoordination parameter re-set to " << over_param[0] << endl ;
  }
  return(params) ;
}
