#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include <string.h>
#include "functions.h"
using namespace std;

static void read_lsq_input(const char *input_filename, int &nlayers, 
			   bool &fit_coul, Sr_pair_t &pair_type, bool &if_subtract_coord,
			   bool &fit_pover, int &cheby_order, double *smin,
			   double *smax, double *sdelta, int &n_over, double *over_param) ;

static void echo_lsq_input(int nlayers, 
			   bool fit_coul, Sr_pair_t pair_type, bool if_subtract_coord,
			   bool fit_pover, int cheby_order, const double *smin,
			   const double *smax, const double *sdelta, int n_over, 
			   const double *over_params)  ;

int main(int argc, char* argv[])
{
  if ( argc != 3 ) {
    cout << "Usage: splines_ls <input filename> <number of frames>\n" ;
    exit(1) ;
  }
  const int nframes=atoi(argv[2]);
  if ( nframes <= 0 ) 
    {
      cout << "Error: a positive number of frames must be specified\n" ;
      exit(1) ;
    }

  ////Job parameters:

  int nlayers=1;//supercells adjacent to central cell.
  //1 is enough because of 8-Ang. limit to spline potential; 
  //electrostatic energy doesn't depend on this number.
  //setting to 2 gives identical result.

  // If true, subtract overcoordination forces.
  bool ifsubtract_coord = false ;  

  // If true, fit coulomb parameters.
  bool fit_coul = true ;

  // If true, fit overcoordination parameters.
  bool fit_pover = true ;

  // Short-range pair potential type.
  Sr_pair_t pair_type ;

  // Order of Chebyshev polynomial if used.
  int cheby_order = 12 ;

  // Number of overcoordination parameters.
  int n_over = 5 ;

  // Values of overcoordination parameters

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
  
  //define spline parameters:
  double smin[NPAIR] ;   // Beginning of spline/fit.
  double smax[NPAIR] ;   // End of spline/fit
  // const double sdelta=0.05;//0.025; //grid spacing
  double sdelta[NPAIR];  // Step size of spline.

  for ( int i = 0 ; i < NPAIR ; i++ ) 
    {
      smax[i] =   6.0 ;
      smin[i] =   0.5 ;
      sdelta[i] = 0.1 ;
    }

  read_lsq_input(argv[1],nlayers, fit_coul, pair_type, ifsubtract_coord, fit_pover, cheby_order,
		 smin, smax, sdelta, n_over, over_param) ;


  echo_lsq_input(nlayers, fit_coul, pair_type, ifsubtract_coord, fit_pover, cheby_order,
		 smin, smax, sdelta, n_over, over_param) ;

  ////Read input file input.xyzf:
  ////xyz format in Angs, three box dimensions in the info line.
  ifstream fileread("input.xyzf");

  std::vector<int> Nat;
  Nat.reserve(nframes);
  std::vector<double*> Latcons;
  Latcons.reserve(nframes);

  std::vector<string*> Lb;
  Lb.reserve(nframes);

  std::vector<double**> Coord;
  Coord.reserve(nframes);
  std::vector<double**> Force;
  Force.reserve(nframes);

  std::vector<double***> Als;//least squares A matrix
  Als.reserve(nframes);

  std::vector<double**> Coul_oo;
  Coul_oo.reserve(nframes);
  std::vector<double**> Coul_oh;
  Coul_oh.reserve(nframes);
  std::vector<double**> Coul_hh;
  Coul_hh.reserve(nframes);

  std::vector<double**> Pover ;
  Pover.reserve(nframes);

  int snum[NPAIR] ;
  int tot_snum ;

  tot_snum = 0 ;
  for ( int i = 0 ; i < NPAIR ; i++ ) 
    {
      if ( pair_type != CHEBYSHEV ) 
	{
	  snum[i]=(2+floor((smax[i]-smin[i])/sdelta[i]))*2;//2 is for p0/m0/p1/m1
	}
      else 
	{
	  snum[i] = cheby_order ;
	}
      tot_snum += snum[i] ;
    }
  cout << "The number of short-range parameters is " << tot_snum << endl ;

  std::vector<double*> Mass;
  Mass.reserve(nframes);

  int tempnat ;
  double *tempLatcons;
  string *tempLb;

  double **tempCoord;
  double **tempForce;
  double ***tempA;
 
  double *tempMass;

  double **tempCoul_oo;
  double **tempCoul_oh;
  double **tempCoul_hh;
  double **tempPover ;

  //Each MD snapshot is pushed back in a std::vector.
  //Coordinates, forces etc. held in pointer arrays.
  for(int N=0;N<nframes;N++)
    {
      fileread >> tempnat;
      tempLatcons=new double[3];
      fileread >> tempLatcons[0] >> tempLatcons[1] >> tempLatcons[2];

      tempLb=new string[tempnat];

      tempCoord=new double *[tempnat];
      tempForce=new double *[tempnat];

      for(int a=0;a<tempnat;a++)
	{
	  tempCoord[a]=new double[3];
	  tempForce[a]=new double[3];
	}
      
      for(int a=0;a<tempnat;a++)
	{
	  fileread>>tempLb[a];
	  for(int c=0;c<3;c++)
	    fileread>>tempCoord[a][c];
	  for(int c=0;c<3;c++)
	    fileread>>tempForce[a][c];
	}


      tempA=new double **[tempnat];
      for(int a=0;a<tempnat;a++)
	{
	  tempA[a]=new double *[tot_snum];
	  for(int n=0;n<tot_snum;n++)
	    {
	      tempA[a][n]=new double[3];
	      for(int c=0;c<3;c++)
		tempA[a][n][c]=0;
	    }
	}

      tempMass=new double[tempnat];
      for(int a=0;a<tempnat;a++)
	{
	  if(tempLb[a]=="O")
	    tempMass[a]=15.9994;
	  else if(tempLb[a]=="H")
	    tempMass[a]=1.00794;
	  else
	    {
	      cerr<<"Atom type "<<tempLb[a]<<" not supported."<<endl;
	      exit(1);
	    }
	}

      tempCoul_oo=new double *[tempnat];
      tempCoul_oh=new double *[tempnat];
      tempCoul_hh=new double *[tempnat];
      tempPover = new double* [tempnat] ;
      for(int a=0;a<tempnat;a++)
	{
	  tempCoul_oo[a]=new double[3];
	  tempCoul_oh[a]=new double[3];
	  tempCoul_hh[a]=new double[3];
	  tempPover[a] = new double [3] ;

	  for(int c=0;c<3;c++)
	    {
	      tempCoul_oo[a][c]=0.0;
	      tempCoul_oh[a][c]=0.0;
	      tempCoul_hh[a][c]=0.0;

	      tempPover[a][c] = 0.0 ;
	    }
	}
      
      Nat.push_back(tempnat);
      Latcons.push_back(tempLatcons);
      Lb.push_back(tempLb);
      Coord.push_back(tempCoord);
      Force.push_back(tempForce);
      Als.push_back(tempA);
      Mass.push_back(tempMass);
      Coul_oo.push_back(tempCoul_oo);
      Coul_oh.push_back(tempCoul_oh);
      Coul_hh.push_back(tempCoul_hh);
      Pover.push_back(tempPover) ;
    }


  ////Wrap coordinates:
  for(int N=0;N<nframes;N++)
    {
      for(int a1=0;a1<Nat[N];a1++)
	for(int c=0;c<3;c++)
	  Coord[N][a1][c] -= floor(Coord[N][a1][c]/Latcons[N][c])*Latcons[N][c];
    }
  ////Convert forces from Hartree/bohr to kcal/Angs (Stillinger's  units).
  for(int N=0;N<nframes;N++)
    for(int a1=0;a1<Nat[N];a1++)
      for(int c=0;c<3;c++)
	Force[N][a1][c] *= 627.50960803*1.889725989;

 
  //Generate A matrix, b vector, for least-squares:

  // cout.precision(10);
  cout.precision(10);

  ////ZCalc_Deriv is the function that calculates elements of A.
  for(int N=0;N<nframes;N++) 
    {
      char *Lbc ;
      
      Lbc = new char [Nat[N]] ;
      for(int a=0;a<Nat[N];a++) {
	if ( Lb[N][a] == "O" ) {
	  Lbc[a] = 'O' ;
	} else if ( Lb[N][a] == "H" ) {
	  Lbc[a] = 'H' ;
	} else {
	  cout << "Error: unknown element " << Lb[N][a] << "\n" ;
	}
      }

      ZCalc_Deriv(Coord[N],Lbc,Latcons[N],nlayers,Nat[N],Als[N],
		  smin,smax,sdelta,snum,Coul_oo[N],Coul_oh[N],Coul_hh[N],
		  pair_type);

      // Subtract over-coordination forces from force to be output.
      if ( ifsubtract_coord ) 
	{
	  SubtractCoordForces(Coord[N], Force[N], Lb[N], Latcons[N], 
			      nlayers, Nat[N], false, Pover[N],
			      n_over, over_param) ;
	}
      if ( fit_pover ) 
	// Fit the overcoordination parameter.
	{
	  SubtractCoordForces(Coord[N], Force[N], Lb[N], Latcons[N], 
			      nlayers, Nat[N], true, Pover[N],
			      n_over, over_param) ;
	}

      delete[] Lbc ;
    }


  ////Print out least squares matrix:
  //Ax=b.
  //A printed to "A.txt", b printed to "b.txt"
  //use these files for python SVD routine.
  ofstream fileA("A.txt");
  ofstream fileb("b.txt");

  //  Reduced precision to 6 for code testing.
  fileA.precision(6);
  fileb.precision(6);

  for(int N=0;N<nframes;N++)
    for(int a=0;a<Nat[N];a++)
      {
	//print Afile:
	for(int c=0;c<3;c++)
	  {
	    //Afile
	    for(int n=0;n<tot_snum;n++) 
	      {
		//	    if ( n % 2 == 1 ) 
		//{
		//Als[N][a][n][c] *= 100.0 ;
		//}
		fileA<<Als[N][a][n][c]<<"   ";
	      }
	  
	    if ( fit_coul ) 
	      {
		fileA<<Coul_oo[N][a][c]<<"   "<<Coul_oh[N][a][c]<<"   "<<Coul_hh[N][a][c];
	      }
	    if ( fit_pover ) 
	      {
		fileA << " " << Pover[N][a][c]  ;
	      }
	    fileA<<endl;
	  }


	  ////print bfile:
	  for(int c=0;c<3;c++)
	    {
	      fileb<<Force[N][a][c]<<endl;
	    }
      }

  //Afile charge conservation (q00 + 4 qOH + qHH = 0 )
  if ( fit_coul ) {
    for(int n=0;n<tot_snum;n++)
      fileA<<"0.0 ";

    fileA<<"1000.0 4000.0 4000.0"; // 0.0"<<endl;
    if ( fit_pover ) 
      {
	fileA << " 0.0 " ;
      }
    fileA << endl ;

    ////bfile charge conservation
    fileb<<"0.0"<<endl;


    //Afile charge consistency (qOO - 4 qHH = 0)
    for(int n=0;n<tot_snum;n++)
      fileA<<"0.0 ";
    fileA<<"1000.0 0000.0 -4000.0" ;
    ////bfile charge consistency
    fileb<<"0.0"<<endl;

    if ( fit_pover ) 
      {
	fileA << " 0.0 " ;
      }
    fileA << endl ;


  } 

  //cout<<endl<<endl;
  FILE *header = fopen("params.header","w") ;
  if ( header == NULL ) 
    {
      cout << "Error: could not write to params.header\n" ;
      exit(1) ;
    }

  fprintf(header, "npair %d\n", NPAIR) ;

  for ( int i = 0 ; i < NPAIR ; i++ ) {
    if ( pair_type != CHEBYSHEV ) {
      fprintf(header, "%8.5f %8.5f %8.5f\n", smin[i], smax[i], sdelta[i]) ;
    } else {
      fprintf(header, "%8.5f %8.5f %d\n", smin[i], smax[i], snum[i]) ;
    }
  }
  
  if ( fit_pover || ifsubtract_coord ) 
    {
      fprintf(header,"%d\n",n_over) ;
      for ( int k = 0 ; k < n_over ; k++ ) 
	{
	  fprintf(header,"%21.13e\n",over_param[k]) ;
	}
    }
  else {
    fprintf(header,"0\n") ;
  }
  fclose(header) ;

  return 0;    
}



static void read_lsq_input(const char *input_filename, int &nlayers, 
			   bool &fit_coul, Sr_pair_t &pair_type, bool &if_subtract_coord,
			   bool &fit_pover, int &cheby_order, double *smin,
			   double *smax, double *sdelta, int &n_over, double *over_params) 
// Read program input from the file "splines_ls.in".
{
  FILE *fin ;
  const int bufsz = 1024 ;
  char buf[bufsz] ;

  fin = fopen(input_filename, "r") ;
  if ( fin == NULL ) 
    {
      cout << "Error: could not open splines_ls.in\n" ; 
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
      if ( val != NULL ) {
	int iend = strlen(val) - 1 ;
	if ( iend > 0 && val[iend] == '\n' ) {
	  val[iend] = 0 ;
	}
      }
      if ( strncmp(name,"nlayers",bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &nlayers) ;
	}
      else if ( strncmp(name,"cheby_order",bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &cheby_order) ;
	}
      else if ( strncmp(name,"smin",bufsz) == 0 ) 
	{
	  sscanf(val, "%lf", &smin[0]) ;
	  parse_param_list(&smin[1], NPAIR-1, "smin") ;
	}
      else if ( strncmp(name,"smax",bufsz) == 0 ) 
	{
	  sscanf(val, "%lf", &smax[0]) ;
	  parse_param_list(&smax[1], NPAIR-1, "smax") ;
	}
      else if ( strncmp(name,"sdelta",bufsz) == 0 ) 
	{
	  sscanf(val, "%lf", &sdelta[0]) ;
	  parse_param_list(&sdelta[1], NPAIR-1, "sdelta") ;
	}
      else if ( strncmp(name,"fit_coulomb",bufsz) == 0 ) 
	{
	  fit_coul = parse_tf(val, bufsz,line) ;
	}
      else if ( strncmp(name,"subtract_coord",bufsz) == 0 ) 
	{
	  if_subtract_coord = parse_tf(val, bufsz, line) ;
	}
      else if ( strncmp(name,"fit_pover",bufsz) == 0 ) 
	{
	  fit_pover = parse_tf(val, bufsz, line) ;
	}
      else if ( strncmp(name, "over_params", bufsz) == 0 ) 
	{
	  sscanf(val, "%d", &n_over) ;
	  if ( val <= 0 ) 
	    {
	      cout << "Error: number of overcoordination params <= 0\n" ;
	      exit(1) ;
	    }
	  else if ( n_over > MAXOVERP ) 
	    {
	      cout << "Error: maximum number of overcoordination params exceeded\n" ;
	      exit(1) ;
	    }
	  parse_param_list(over_params, n_over, "overcoordination") ;
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
#if(0)
	  // NOT SUPPORTED FOR LEAST-SQUARES FITTING
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
#endif
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
      printf("Error while reading splines_ls.in\n") ;
      exit(1) ;
    }
}



static void echo_lsq_input(int nlayers, 
			   bool fit_coul, Sr_pair_t pair_type, bool if_subtract_coord,
			   bool fit_pover, int cheby_order, const double *smin,
			   const double *smax, const double *sdelta, int n_over, 
			   const double *over_params) 
// Read program input from the file "splines_ls.in".
{
  printf("Number of layers used in force evaluation = %d\n", nlayers) ;

  if ( pair_type == CHEBYSHEV ) 
    {
      printf("Chebyshev pair interaction used\n") ;
      printf("Number of chebyshev terms = %d\n", cheby_order) ;
    }
  else if ( pair_type == SPLINE ) 
    {
      printf("Spline pair interaction used\n") ;

      for ( int i = 0 ; i < NPAIR ; i++ ) {
	printf("Spline point spacing for pair %d = %13.6e\n", i, sdelta[i]) ;
      }
    }
  else if ( pair_type == INVERSE_R ) 
    {
      printf("Inverse power pair interaction used\n") ;
    }
  else if ( pair_type == LJ ) 
    {
      printf("Lennard-Jones pair interaction used\n") ;
    }
  else if ( pair_type == STILLINGER ) 
    {
      printf("Stillinger pair interaction used\n") ;
    }

  for ( int i = 0 ; i < NPAIR ; i++ ) {
    printf("Pair interaction minimum distance for pair %d = %13.6e\n", i, smin[i]) ;
    printf("Pair interaction maximum distance for pair %d = %13.6e\n", i, smax[i]) ;
  }

  if ( fit_coul ) 
    {
      printf("Coulomb parameters will be determined by fitting\n") ;
    }
  else 
    {
      printf("Coulomb interactions will not be considered\n") ;
    }
  
  if ( if_subtract_coord ) 
    {
      printf("Overcoordination forces will be subtracted before fitting\n") ;
    }
  else if ( fit_pover ) 
    {
      printf("Overcoordination linear parameter will be fit to forces\n") ;
    }

  if ( if_subtract_coord ) 
    {
      printf("Linear overcoordination parameter pover = %13.6e\n", over_params[0] ) ;
    }
  if ( if_subtract_coord || fit_pover ) 
    {
      printf("Overcoordination r0 = %13.6e\n", over_params[1]) ;
      printf("                 p1 = %13.6e\n", over_params[2]) ;
      printf("                 p2 = %13.6e\n", over_params[3]) ;
      printf("             lambda = %13.6e\n", over_params[4]) ;
    }
}
