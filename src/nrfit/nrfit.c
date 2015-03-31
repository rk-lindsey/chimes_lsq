//
// 
// nrfit is a driver for non-linear optimization of a general
// objective function.  Most of the algorithms are derived from 
// Numerical Recipes in C.  The objective function is implemented
// through a system call.  The system call can be a perl script, etc.
// that executes a program and reports a value to be optimized.
//
// Usage:  nrfit -option option_file < fit.in >& fit.out
// 
//
//
// Parameters that do not depend on the optimization algorithm are
// read from standard input (fit.in in the above example).
//
// Line 1 :  The name of the command that evaluates the objective function.
// Line 2 :  An input file which specifies parameters that the program in 
//           line 1 will read.
// Line 3 :  An output file which specifies parameters that the program 
//           in line 1 will write.  The output file should have a single 
//           line, which is the value of the objective function.
// Line  4 : The number of parameters to optimize (N).
// Lines 5 : 4 + N : Initial values of the parameters.
//
// The options control which optimization methods are used.  A
// global optimization method may be followed by a local optimization
// method by specifying two options on the command line.
//
// The options are as follows:
//  -p : Powell's local optimization method.  See Numerical Recipes.
//  -d file : Davidson-Fletcher-Powell variable metric local optimization method.  
//            See Numerical Recipes.
//             File contents:
//             Line 1:   Optimization function value tolerance.
//
//  -b : The Downhill Simplex (amoeba) local optimization method.  See Numerical Recipes.
//  -s file : Simulated annealing global optimization method.
//             File contents:
//             Line 1 :  Number of cycles
//             Line 2 :  Initial temperature. Final temperature.
//                       The initial temperature should be comparable
//                       to the value of the function.  The final temperature
//                       should be 2-3 orders of magnitude less than the
//                       initial temperature.
//             Line 3:3+N-1 :  Maximum MC step (1 line per parameter). 
//
//  -a  file : Amoeba simulated annealing method.  See Numerical Recipes.
//             File contents:  
//             Line 1: Initial temperature.  Set to a value somewhat less
//                     than a typical value of the objective function.
//             Line 2: Number of blocks:  Set to 10 or greater.  The temperature
//                       is reduced by a factor of 2 each block.
//	       Line 3: Number of steps in a block.
//  -r file : Record to record global optimization method.
//            See G. Dueck, "New optimization hueristics: The great
//            deluge algorithm and the record to record travel",
//            J. Comp. Phys., 104, 86-92(1993).
//
//             File contents:
//             Line 1: Number of cycles to perform.  1 cycle is N trials,
//                     where N is the number of variables to optimize.
//             Line 2 : N + 1: Maximum MC step size.  1 line per variable.
//             Line 2+N :  The deviation parameter.  This is an allowed
//                         fractional tolerance from the best value 
//                         obtained so far.
//
//
//  -g file : Great deluge algorithm global optimization method.
//            See G. Dueck, "New optimization hueristics: The great
//            deluge algorithm and the record to record travel",
//            J. Comp. Phys., 104, 86-92(1993).
//            The algorithm as modified in "Comparison of the performance
//            of modern heuristics for combinatorial optimization of
//            real data", Computers Ops. Res., 20, 687-695(1995) is
//            implemented.  The change in the water level is proportional
//            to the decrease in the objective function.
//            File contents:
//             Line 1: Number of cycles to perform.  1 cycle is N trials,
//                     where N is the number of variables to optimize.
//             Line 2: Maximum MC step size.  1 line per variable.
//             Line 2+N :  The water level change parameter.  Larger
//                         values will converge more rapidly to 
//                         a local minimum.  Smaller values will wander more.

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "dnr.h" 
#include "dnrutil.h" 

typedef double (*doublefn)(double*) ;
double objective_fn(double*) ;
// double min(double x, double y) { return ( x < y ? x : y ) ; }

int GDA, SA, POW, RT, DFP, ASA, AMOEBA;    // commandline options
char infilename[128], outfilename[128], command[128] ;
const  int BLOCK = 8;
int ncoeff;

void equilibrate(int ncycles, double* step, double* rate, 
                 double beta, doublefn objective_fn, double* coeff,
                 int count, int ncoeff) ;
double anneal(double* coeff, double* maxstep, int ncycles,double beta0
	      ,double beta,doublefn objective_fn, int ncoeff)  ;
double  deluge (double* coeff,double* maxstep, int gcycles, 
		doublefn objective_fn, double *UPcoeff, int ncoeff) ;
void gequilibrate(int ncycles, double* step, double* rate, 
		  double *WATER, doublefn objective_fn, double* coeff, 
		  double *UPcoeff, int count, int ncoeff) ;
double record(double* coeff,double* maxstep, int gcycles, 
	      doublefn objective_fn, double deviation, 
	      double *bestcoeff, int ncoeff) ;
void requilibrate(int ncycles, double* step, double* rate, 
		  doublefn objective_fn, double* coeff, double deviation,
		  double *minpotential, double *bestcoeff,
		  int count, int ncoeff) ;

#define MAXCOEFF 100

int main (int argc, char **argv) 
{
     
  int option_char;		// commmandline option
  int i ;
  int ncycles, gcycles;
  double beta0, beta, deviation, UP, *UPcoeff, dfptol ;
  char *rfile, *gfile, *sfile, *dfile, *asafile ; 
  double asa_t ;
  int asa_cyc, asa_block ;
  double coeff[MAXCOEFF] ;
  double bestcoeff[MAXCOEFF] ;
  double maxstep[MAXCOEFF], gmaxstep[MAXCOEFF] ;
  FILE *infile ;
  const char *option_string = "pg:s:r:d:a:b" ;
  int junk ;

  UPcoeff = &UP;     

  while ((option_char = getopt(argc, argv, option_string)) != -1 ) {
    switch (option_char) {
    case 'p' : POW = 1; break;
    case 'g' : GDA = 1; gfile = optarg; break;
    case 's' : SA = 1;  sfile = optarg; break;
    case 'a' : ASA = 1 ; asafile = optarg ; break ;
    case 'r' : RT = 1;  rfile = optarg; break; 
    case 'd' : DFP = 1; dfile = optarg; break ;
    case 'b' : AMOEBA = 1; break ;
    default :
      fprintf(stderr,"Bad option_char\n") ;
      exit(1) ;
    } 
  }

  scanf("%s",command) ;
  scanf("%s",infilename) ;
  scanf("%s",outfilename) ;
  printf("Evaluation command = %s\n" ,command) ;
  printf("Evaluation input file = %s\n",infilename) ;
  printf("Evaluation output file = %s\n",outfilename ) ;

  scanf("%d", &ncoeff) ;
     
  printf("Number of parameters to fit = %d\n", ncoeff) ;
  printf("Initial guess for parameters \n" ) ;
  for ( i = 0 ; i < ncoeff ; i ++ ) {
    scanf("%lf",&coeff[i]) ;
    printf("%11.4e\n",coeff[i]) ;
  }
  
  if ( ASA ) {		// Input for Amoeba Simulated Annealing 

    infile = fopen(asafile,"r") ;
    if ( infile == NULL ) {
      printf("Error: could not open %s\n", asafile) ;
      exit(1) ;
    }
    printf("Amoeba Simulated Annealing (ASA) chosen\n") ;

    fscanf(infile,"%lf",&asa_t) ;     // Initial temperature
    printf("ASA Temperature = %11.4e\n", asa_t) ;
    fscanf(infile,"%d",&asa_cyc) ;   // Number of annealing cycles
    printf("Number of ASA blocks = %d\n", asa_cyc) ;

    fscanf(infile,"%d",&asa_block) ; // Number of steps in a block.
    printf("Number of steps in a block = %d\n", asa_block) ;

    if ( feof(infile) ) {
      printf("Error: could not read all parameters from %s\n", asafile) ;
      exit(1) ;
    }
    fclose(infile) ;
  }				// End of Simulated Annealing Input

  if ( SA ) {		// Input for Simulated Annealing 

    printf("Simulated Annealing Chosen\n") ;
    infile = fopen(sfile,"r") ;
    if ( infile == NULL ) {
      printf("Error: could not open %s\n", sfile) ;
      exit(1) ;
    }

    fscanf(infile,"%d",&ncycles) ;
    printf("Number of cycles = %d\n",ncycles) ;

    fscanf(infile,"%lf %lf", &beta0, &beta) ;

    printf("Initial temperature = %f\n", beta0 ) ;
    printf("Final temperature   = %f\n",beta ) ;
           
    /* Transform to inverse temperature */

    if ( beta0 < 1.0e-20 ) {
      beta0 = 1.0e-20 ;
    }
    if ( beta < 1.0e-20 ) {
      beta = 1.0e-20 ;
    }
    
    beta0 = 1.0 / beta0 ;
    beta  = 1.0 / beta ;

    for ( i = 0 ; i< ncoeff ; i++ ) {
      fscanf(infile,"%lf", &maxstep[i]) ;
    }
    printf( "Maximum MC Step size \n" ) ;
    for ( i = 0 ; i< ncoeff ; i++ ) {
      printf("%lf\n", maxstep[i] ) ;
      if ( maxstep[i] < 0 ) {
	printf("Error : Negative MC step size\n") ;
	exit(1) ;
      }
    }
    fclose(infile) ;      
  }				// End of Simulated Annealing Input

     
  if (GDA) {
	  
    infile = fopen(gfile,"r") ;
    if ( infile == NULL ) {
      printf("Error: could not open %s\n", gfile) ;
      exit(1) ;
    }

    fscanf(infile,"%d",&gcycles) ;	//       Great Deluge Algorithm Input 
    printf("Number of cycles = %d\n",gcycles) ;

    for ( i = 0 ; i< ncoeff ; i++ ) {
      fscanf(infile,"%lf", &gmaxstep[i]) ;
    }
    printf("Maximum MC Step size \n") ;
    for (i= 0 ; i < ncoeff; i++) { 
      printf("%f\n", gmaxstep[i]) ;
      if ( gmaxstep[i] < 0) {  
	printf(" Error : Negative MC step size \n") ;
	exit(1);
      }
    }       
	                    
    fscanf(infile,"%lf",UPcoeff) ;
    printf(" UP = %f\n", *UPcoeff) ;

    
    fclose(infile) ;

  } 
    
  if (RT) {			// Record to Record 
    infile = fopen(rfile,"r");
    if ( infile == NULL ) {
      printf("Error: could not open %s\n", rfile) ;
      exit(1) ;
    }

    fscanf(infile,"%d",&gcycles) ;	
    printf("Number of cycles = %d\n", gcycles) ;

    for ( i = 0 ; i< ncoeff ; i++ ) {
      fscanf(infile,"%lf", &gmaxstep[i]) ;
    }
    printf("Maximum MC Step size \n") ;
    for (i= 0 ; i < ncoeff; i++) { 
      printf("%f\n", gmaxstep[i]) ;
      if ( gmaxstep[i] < 0) {  
	printf(" Error : Negative MC step size \n") ;
	exit(1);
      }
    }       
    fscanf(infile,"%lf",&deviation) ;
    printf( "Deviation = %f\n", deviation) ;

    fclose(infile) ;
  }    

  if (DFP ) {			// DFP optimization
    infile = fopen(dfile,"r");

    if ( infile == NULL ) {
      printf("Error: could not open DFP input\n") ;
      exit(1) ;
    }
    fscanf(infile,"%lf",&dfptol) ;           // Tolerance in dfpmin
    printf( " DFP tolerance = %f\n",dfptol) ;
    fclose(infile) ;
  }    

  if ( ferror(stdin) ) {
    printf("End of file reached prematurely \n" ) ;
    exit(1) ;
  }
  scanf("%d", &junk) ;
  if ( ! feof(stdin) ) {
      printf( "Too many lines in input file \n") ;
      exit(1) ;    
  }

  double fret = 1.0 ;    
  /**** First do simulated annealing ****/
  if ( SA) {
    FILE *coeff_file ;

    fret = anneal(coeff,maxstep, ncycles, beta0, beta, &objective_fn, ncoeff) ;

    coeff_file = fopen("scoeff.out","w");
    printf("Final parameter values\n") ;
    for ( i = 0 ; i < ncoeff ; i++ ) {
      printf("%13.6e\n",coeff[i]) ;
      fprintf(coeff_file, "%13.6e SA\n", coeff[i])  ;
    }
    printf("Value after annealing = %f\n",fret) ;
  }


  //   Great Deluge Algorithm here 
  if (GDA) {
    FILE *coeff_file ;

    coeff_file = fopen("gcoeff.out","w");
    
    printf("Performing Great Deluge Algorithm Optimization\n") ;

    fret = deluge (coeff, gmaxstep,  gcycles, objective_fn, 
		   UPcoeff, ncoeff);

    for ( i = 0 ; i < ncoeff ; i++ ) {
      printf("%f\n", coeff[i]) ;
      fprintf(coeff_file, "%f GDA\n", coeff[i]) ;
    }
		 
    printf(" Value after Great Deluge Algorithm = %lf\n",fret) ;
  }   

  // Record to Record Minimization

  if (RT) {
    FILE *coeff_file ;

    coeff_file = fopen("rcoeff.out","w");

    fret = record(coeff, gmaxstep, gcycles, objective_fn,
		  deviation, bestcoeff, ncoeff);

    printf("Final parameter values found:\n") ;
    for ( i = 0 ; i < ncoeff ; i++ ) {
      printf("%13.7e\n", bestcoeff[i]) ;
      fprintf(coeff_file, "%13.7e RT\n", bestcoeff[i]) ;
    }
    printf(" Value after Record to Record Minimization = %f\n",fret ) ;
  }

  if ( DFP ) {
    // DFP variable metric minimization algorithm.
    int iter = 0 ;
    FILE *coeff_file ;

    coeff_file = fopen("dcoeff.out","w") ;
    printf("Performing DFP local optimization\n") ;

    dfpmin(coeff-1, ncoeff, dfptol, &iter, &fret, objective_fn) ;
    
    printf("Value after DFPMIN = %lf\n",fret) ;
    printf("Number of iterations = %d\n",iter) ;

    for ( i = 0 ; i < ncoeff ; i++ ) {
      printf("%le\n", coeff[i]) ;
      fprintf(coeff_file, "%le P\n",coeff[i] ) ;
    }
  }
  if ( ASA )  {
    // Amoeba simulated annealing.
    double yb ;
    int iter = 0 ;
    double **xi, *y, *p ;
    FILE *coeff_file;

    coeff_file = fopen("asacoeff.out","w") ;          
    int j ;
	  
    xi = matrix(1,ncoeff+1,1,ncoeff) ;
    y  = vector(1,ncoeff) ;
    p = vector(1,ncoeff+1) ;

    for ( j = 1 ; j <= ncoeff ; j++  ) {
      xi[ncoeff+1][j] = coeff[j-1] ;
    }
    y[ncoeff+1] = objective_fn(xi[ncoeff+1]) ;          
    for ( i = 1 ; i <= ncoeff ; i++ ) {
      for ( j = 1 ; j <= ncoeff ; j++  ) {
	xi[i][j] = coeff[j-1] ;
      }
      if ( coeff[i-1] != 0.0 ) {
	xi[i][i] += 0.01 * coeff[i-1] ;
      } else {
	xi[i][i] += 1.0e-05 ;
      }
      y[i] = objective_fn(xi[i]) ;
    }
    yb = 1.0e10 ;
    for ( i = 1 ; i <= asa_cyc ; i++ ) {
      yb = 1.0e12 ;
      iter = asa_block ;
      printf("ASA Temperature = %12.4e\n", asa_t) ;

      amebsa(xi, y, ncoeff, p, &yb, 1.0e-05, &objective_fn,
	     &iter, asa_t) ;
      printf("Best ASA error = %f; current error = %f\n", yb, y[1]) ;
      asa_t /= 2.0 ;
      printf("Current temperature = %12.5e\n", asa_t) ;

    }
    j = 1 ;
    for ( i = 2 ; i <= ncoeff + 1 ; i++ ) {
      if ( y[i] < y[j] ) {
	j = i ;
      }
    }
    printf( "Value after amoeba simulated annealing = %f\n", y[j]) ;
    printf("Final parameter values:\n") ;
    for ( i = 1 ; i <= ncoeff ; i++ ) {
      printf("%13.6e\n",xi[j][i]) ;
      fprintf(coeff_file,"%13.6e A\n",xi[j][i]) ;
    }
  }


  /**** Then clean up with powell's method ****/
  /** Old code using Powell's method **/
  fret = 1.0 ;
  if (POW)  {
    double ftol = 1.0e-04 ;
    int iter = 0 ;	  
    double **xi ;
    FILE *coeff_file; 
    int j ;

    printf("Performing Powell local optimization\n") ;

    coeff_file = fopen("pcoeff.out","w") ;

    xi = matrix(1,ncoeff,1,ncoeff) ;
    for ( i = 1 ; i <= ncoeff ; i++ ) {
      for ( j = 1 ; j <= ncoeff ; j++  ) {
	xi[i][j] = 0.0 ;
      }
      if ( coeff[i-1] != 0.0 ) {
	xi[i][i] = 0.01 * coeff[i-1] ;
      } else {
	xi[i][i] = 1.0e-05 ;
      }
    }
    powell(coeff-1,xi,ncoeff,ftol,&iter,&fret,&objective_fn) ;
    //  coeff-1    initial starting point
    //  xi               initial matrix dimensions ncoeff * ncoeff (logical)
    // ftol  failure to decrease more than this amount in 1 iteration means done
    // &iter             number of iterations taken
    // &fret             the returned function value at coeff.vec()-1
    printf( "Value after powell = %lf\n",fret) ;
    printf( "Number of iterations = %d\n",iter ) ;

    printf("Final parameter values:\n") ;
    for ( i = 0 ; i < ncoeff ; i++ ) {
      printf( "%13.6e\n", coeff[i]) ;
      fprintf(coeff_file, "%13.6e P\n", coeff[i]) ;
    }
  }

  if ( AMOEBA )  {
    // Amoeba Downhill simplex optimization
    int iter = 0 ;
    double **xi, *y ;
    FILE *coeff_file;

    printf("Performing amoeba (simplex) local optimization\n") ;

    coeff_file = fopen("bcoeff.out","w") ;          
    int j ;
	  
    xi = matrix(1,ncoeff+1,1,ncoeff) ;
    y  = vector(1,ncoeff+1) ;

    for ( j = 1 ; j <= ncoeff ; j++  ) {
      xi[ncoeff+1][j] = coeff[j-1] ;
    }
    y[ncoeff+1] = objective_fn(xi[ncoeff+1]) ;          
    for ( i = 1 ; i <= ncoeff ; i++ ) {
      for ( j = 1 ; j <= ncoeff ; j++  ) {
	xi[i][j] = coeff[j-1] ;
      }
      if ( coeff[i-1] != 0.0 ) {
	xi[i][i] += 0.01 * coeff[i-1] ;
      } else {
	xi[i][i] += 1.0e-05 ;
      }
      y[i] = objective_fn(xi[i]) ;
    }
    amoeba(xi, y, ncoeff, 1.0e-05, &objective_fn, &iter) ;
    j = 1 ;
    for ( i = 2 ; i <= ncoeff + 1 ; i++ ) {
      if ( y[i] < y[j] ) {
	j = i ;
      }
    }
    printf( "Value after amoeba optimization = %f\n", y[j]) ;
    printf(  "Number of iterations = %d\n",iter) ;

    for ( i = 1 ; i <= ncoeff ; i++ ) {
      printf("%le\n",xi[j][i]) ;
      fprintf(coeff_file,"%le B\n",xi[j][i]) ;
    }
  }

  return(0) ;
}

double anneal(double* coeff, double* maxstep, int ncycles,double beta0
	      ,double beta,doublefn objective_fn, int ncoeff) 
/** This is a driver for simulated annealing optimization. **/
{
    
     double rate[MAXCOEFF] ;
     double step[MAXCOEFF] ;
     int nstage ;
     int i ;
     int j ;
     int nextra ;
     double retval ;

     for ( j = 0 ; j < ncoeff ; j++ ) {
       rate[j] = 0.0 ;
       step[j] = maxstep[j] / 10.0 ;
     }

     const int NBLOCK = BLOCK ;	/* Block size used in adjusting MC step
				   size */

     nstage = ncycles / NBLOCK - 1 ;

     for ( j = 0 ; j < ncycles / NBLOCK ; j++ ) {
       double betac ;
       if ( j > 0 ) {
	 betac = beta0 * pow(beta/beta0,(1.0*j)/(nstage-1)) ;
       } else {
	 betac = beta0 ;
       }
       printf("Current temperature = %12.5e\n", 1.0/betac) ;
	  
       printf("Current MC step size:\n") ;
       for ( i = 0 ; i < ncoeff ; i++ ) {
	 printf("%12.5e\n", step[i]) ;
       }

       equilibrate(NBLOCK, step, rate, betac, objective_fn, coeff, j,ncoeff) ;

       printf("Acceptance rate:\n") ;
       for ( i = 0 ; i < ncoeff ; i++ ) {
	 printf("%12.5e\n", rate[i]) ;
       }

       /*** Adjust monte carlo step size to get a good acceptance rate ***/
       for ( i = 0 ; i < ncoeff ; i++ ) {
	 if ( rate[i] < 0.4 ) {
	   step[i]   /= 1.5 ;
	 } else if ( rate[i] > 0.6 && step[i] < maxstep[i] ) {
	   step[i]   *= 1.5 ;
	 }
       }
     }
     nextra = ncycles - (ncycles / NBLOCK  ) * NBLOCK ;
     equilibrate(nextra,step,rate,beta,objective_fn,coeff, ncycles/NBLOCK,ncoeff) ;
 
     retval = objective_fn(coeff-1) ;
     return(retval) ;
}

void equilibrate(int ncycles, double* step, double* rate, 
                 double beta, doublefn objective_fn, double* coeff,
                 int count, int ncoeff) 
/**
  Performs ncycles Monte Carlo cycles of simulated annealing.
  beta is the inverse temperature.
**/
{
  int i, j ;
  double oldpoten, oldcoeff, poten, dpoten ;

  for ( i = 0 ; i < ncoeff ; i++ ) {
    rate[i] = 0.0 ;
  }

  if ( ncycles == 0 ) {
    return ;
  }
     
  oldpoten = (*objective_fn)(coeff-1) ;

  for ( i = 0 ; i < ncycles ;  i++ ) {
    for ( j = 0 ; j < ncoeff ; j++ ) {
      oldcoeff = coeff[j] ;
      coeff[j] += step[j] * (drand48()  - 0.5) ;
      poten =  (*objective_fn)(coeff-1) ;
      dpoten = poten - oldpoten ;
	       
      if ( dpoten < 0 ) {
	rate[j]++ ;
	oldpoten = poten ;
      }
      else if ( exp(- beta * dpoten) > drand48() ) {
	rate[j]++ ;
	oldpoten = poten ;
      }
      else {
	coeff[j] = oldcoeff ;
      }
        
    }
    printf("Cycle %d Value %f\n", 
	   ncycles * count + i , (objective_fn)(coeff-1)) ;
  }
  for ( j = 0 ; j < ncoeff ; j++ ) {
    rate[j] /= 1.0 * ncycles ;
  }
}

double  deluge (double* coeff,double* maxstep, int gcycles, 
                 doublefn objective_fn, double *UPcoeff, int ncoeff)
/** Implements the great deluge algorithm, which is a form of
    simplified simulated annealing.  A notional water level is
    used to accept moves.  If the objective function is less than
    the water level, the move is accepted.  The water level is slowly
    decreased as lower values of the objective function are found.

    Some of the parameter names may seem odd, because the original
    GDA was developed for maximization.  I have changed the signs
    of the variables, but not their names, for the minimization problem.
    That means the water level goes down with time, not up as in a deluge.
**/
{
   
     double rate[MAXCOEFF] ;
     double step[MAXCOEFF] ;
     double *WATER, WATER_LEVEL;
     int j ;
     const int NBLOCK = BLOCK ;	
     int i ;
     int nextra ;

     WATER_LEVEL = (*objective_fn)(coeff-1);
     WATER = &WATER_LEVEL;

     for ( j = 0 ; j < ncoeff ; j++ ) {
       rate[j] = 0.0 ;
       step[j] = maxstep[j] / 10.0 ;
     }


     /* Block size used in adjusting MC step size */


     for ( j = 0 ; j < gcycles / NBLOCK ; j++ ) {
	  gequilibrate(NBLOCK, step, rate, WATER,objective_fn, coeff, UPcoeff,
		       j, ncoeff);
	  /*** Adjust monte carlo step size to get a good acceptance rate ***/
	  printf("Acceptance rates:\n") ;
	  for ( i = 0 ; i < ncoeff ; i++ ) {
	    printf("%12.4e\n", rate[i]) ;
	    if ( rate[i] < 0.4 ) {
	      step[i]   /= 1.5 ;
	    } else if ( rate[i] > 0.6 && step[i] < maxstep[i] ) {
	      step[i]   *= 1.5 ;
	    }
	  }
     }
     nextra = gcycles - (gcycles / NBLOCK  ) * NBLOCK; 
     gequilibrate(nextra, step, rate, WATER, objective_fn, coeff, UPcoeff,
                  gcycles / NBLOCK, ncoeff ) ;
     return (objective_fn(coeff-1)) ;
}

void gequilibrate(int ncycles, double* step, double* rate, 
		  double *WATER, doublefn objective_fn, double* coeff, 
		  double *UPcoeff, int count, int ncoeff) 
/**
  Performs Great Deluge Algorithm Monte Carlo.
**/
{
  int j;
  int i ;
  double oldcoeff ;
  double poten ;
  double dpoten ;
  double UP ;

  for ( j = 0 ; j < ncoeff ; j++ ) {
    rate[j]= 0.0 ;
  }
  
  if ( ncycles == 0 ) {
    return ;
  }
     
  for ( i = 0 ; i < ncycles ; i++ ) {
    for ( j = 0 ; j < ncoeff ; j++ ) {
      oldcoeff = coeff[j] ;
      coeff[j] += step[j] * (drand48()  - 0.5) ;
      poten =  (*objective_fn)(coeff-1) ;
      dpoten = poten - *WATER ;
      
      if ( dpoten < 0 ) {
	UP =  *UPcoeff * ( *WATER - poten);                
	//		     double UP = *UPcoeff;
	
	rate[j]++ ;
	*WATER -= UP ;
	printf("Water Level = %12.5e\n", *WATER) ;
      }
      else {
	//                    double UP = *UPcoeff * ( oldpoten - *WATER);
	//		      *WATER -= UP ;
	coeff[j] = oldcoeff ;
      }
    }
    printf(" Cycle %d Value %11.5f\n",count * ncycles + i,(*objective_fn)(coeff-1)) ;
  }
  for ( j = 0 ; j < ncoeff ; j++ ) {
    rate[j] /= 1.0 * ncycles ;
  }
}

double record(double* coeff,double* maxstep, int gcycles, 
	      doublefn objective_fn, double deviation, 
	      double *bestcoeff, int ncoeff)
/** Implements global optimization through the record to record
    algorithm.  All steps with objective function less than
    the running record + deviation * record are accepted.
**/
{
   
  int i, j ;

  double rate[MAXCOEFF] ;
  double step[MAXCOEFF] ;
  double minpotential ; 
  double *min;
  const int NBLOCK = BLOCK ;	
  int nextra ;

  minpotential =  (*objective_fn)(coeff-1);

     
  min = &minpotential;
 
  for ( j = 0 ; j < ncoeff ; j++ ) {
    rate[j] = 0.0 ;
    step[j] = maxstep[j] / 10.0 ;
  }


  /* Block size used in adjusting MC step size */

  for ( j = 0 ; j < gcycles / NBLOCK ; j++ ) {
    printf("Step sizes for variables:\n") ;
    for ( i = 0 ; i < ncoeff ; i++ ) {
      printf("%7.4f\n", step[i]) ;
    }
    requilibrate(NBLOCK, step,rate,objective_fn,coeff,deviation, min, 
		 bestcoeff, j, ncoeff);
    /*** Adjust monte carlo step size to get a good acceptance rate ***/
    printf("Acceptance rates for variables:\n") ;
    for ( i = 0 ; i < ncoeff ; i++ ) {
      printf("%11.5f\n", rate[i]) ;

      if ( rate[i] < 0.4 ) {
	step[i]   /= 1.5 ;
      } else if ( rate[i] > 0.6 && step[i] < maxstep[i] ) {
	step[i]   *= 1.5 ;
      }
    }
  }
  nextra = gcycles - (gcycles / NBLOCK  ) * NBLOCK; 
  requilibrate(nextra, step, rate, objective_fn, coeff, deviation, min,
	       bestcoeff, gcycles / NBLOCK, ncoeff);
  return(minpotential) ;
}

void requilibrate(int ncycles, double* step, double* rate, 
		  doublefn objective_fn, double* coeff, double deviation,
		  double *minpotential, double *bestcoeff,
		  int count, int ncoeff) 
/**
  Performs ncycles Monte Carlo cycles.
  beta is the inverse temperature.
**/
{
  int j ;
  int i ;
  double oldcoeff ;
  double poten ;
  double devi ;
  double dpoten ;

  for ( j = 0 ; j < ncoeff ; j++ ) {
    rate[j] = 0.0 ;
  }

  if ( ncycles == 0 ) {
    return ;
  }
     
  for ( j = 0 ; j < ncoeff ; j++ ) {
    bestcoeff[j] = coeff[j] ;
  }
  for ( i = 0 ; i < ncycles ; i++ ) {
    for ( j = 0 ; j < ncoeff ; j++ ) {
      oldcoeff = coeff[j] ;
      coeff[j] += step[j] * (drand48()  - 0.5) ;
      poten =  (*objective_fn)(coeff-1) ;
      devi = *minpotential * deviation;
      dpoten = poten - ( *minpotential + devi);
	       	      
      if ( dpoten < 0 ) {
	rate[j]++ ;
	if (poten < *minpotential) {
	  *minpotential = poten;
	  bestcoeff[j] = coeff[j] ;
	}
      }
      else 
	coeff[j] = oldcoeff ;
                      
        
    }
    printf("Cycle %d Value %12.5e\n", ncycles * count + i,
	   (*objective_fn)(coeff-1) ) ;
  }
  for ( j = 0 ; j < ncoeff ; j++ ) {
    rate[j] /= 1.0 * ncycles ;
  }
  printf("Best objective function value = %12.5e\n", *minpotential) ;
  printf("Best parameter values:\n") ;
  for ( j = 0 ; j < ncoeff ; j++ ) {
    printf("%12.5e\n", bestcoeff[j]) ;
  }
}

double objective_fn(double * coeff ) 
{
  double result = 0 ;
  FILE *infile, *outfile ;
  int i ;

  infile = fopen(infilename, "w") ;
  if ( infile == NULL ) {
    printf( "Error: could not open file %s\n", infilename) ;
    exit(1) ;
  }

  for ( i = 1 ; i <= ncoeff ; i++ ) {
    printf("Parameter %d = %21.13e\n",i, coeff[i]) ;
    fprintf(infile,"%21.13e\n",coeff[i]) ;
  }
  if ( ferror(infile) ) {
    printf("Error writing new input file\n" ) ;
    exit(1) ;
  }
  fclose(infile) ;

  system(command) ;

  outfile = fopen(outfilename, "r") ;
  // Read results of running the program.
  fscanf(outfile,"%lf",&result) ;

  if ( ferror(outfile) ) {
    printf( "Objective_Fn: Error reading file\n" ) ;
    exit(1) ;
  }
  
  fscanf(outfile,"%d",&i) ;
  if ( ! feof(outfile) ) {
    printf("Objective_Fn: too many numbers in file\n") ;
    exit(1) ;
  }
  fclose(outfile) ;
  
  printf("   Objective function = %21.13e\n", result) ;
  fflush(stdout) ;

  return(result) ;

}
