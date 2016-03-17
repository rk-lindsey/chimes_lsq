#include "functions.h"
//
// 3-Body Chebyshev morse interaction.
//

static void find_pair_cheby(double Rab[3], double *Tn, double *Tnd, double &rlen, double &exprlen,
			    double &xdiff,
			    int a1, int a2, const char *Lbc, double **Coord, const double *Latcons, 
			    const double *smin, const double *smax, 
			    const int *snum, const double *lambda) ;
static double Cheby_3B_Coeff(int a1, int a2, int a3, int n12, int n23, int n13,
			     const char *Lbc, const double *params,
			     const int *snum)  ;
void sort_triple(int &ipair12, int &ipair23, int &ipair13, 
		 int &n12, int &n23, int &n13) ;
static void cubic_cutoff(double &fcut, double &dfcut, double rlen, double smax) ;

void ZCalc_3B_Cheby(double **Coord,const char *Lbc, double *Latcons,
			   const int nat,const double *smin,
			   const double *smax,
			   const int *snum, 
			   double *params, const double *lambda,
			   double **SForce, double &Vtot, double &Pxyz)
// Calculate 3-body short-range forces using a Chebyshev polynomial expansion.
// Can use morse variables similar to the work of Bowman.
{
  double R12[3], R23[3], R13[3];
  double rlen12, rlen23, rlen13 ;
  static double *Tn12, *Tnd12 ;
  static double *Tn23, *Tnd23 ;
  static double *Tn13, *Tnd13 ;
  static bool called_before = false ;
  //spline term calculated w/cutoff:
  double tempx = 0.0;
  //const double lambda = 1.25 ;
  double exprlen12, exprlen23, exprlen13 ;
  double xdiff12, xdiff23, xdiff13 ;


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
    Tn12   = new double [dim] ;
    Tnd12  = new double [dim] ;
    Tn23   = new double [dim] ;
    Tnd23  = new double [dim] ;
    Tn13   = new double [dim] ;
    Tnd13  = new double [dim] ;
  }
  
  ////main loop for 3-body Chebyshev terms:
  for(int a1=0;a1<nat-1;a1++)
    for(int a2=a1+1;a2<nat;a2++) 
      {
	double fcut12, dfcut12 ;
	int ipair12 = pair_index(a1,a2,Lbc) ;
	find_pair_cheby(R12,Tn12, Tnd12, rlen12, exprlen12, xdiff12,
			a1, a2, Lbc, Coord, Latcons, smin, smax, snum, lambda) ;
	
	cubic_cutoff(fcut12, dfcut12, rlen12, smax[ipair12]) ;
	if ( fcut12 == 0.0 ) 
	  continue ;

	for(int a3=a2+1;a3<nat;a3++)
	  {
	    int ipair23 = pair_index(a2,a3,Lbc) ;
	    int ipair13 = pair_index(a1,a3,Lbc) ;

	    find_pair_cheby(R23, Tn23, Tnd23, rlen23, exprlen23, xdiff23,
			    a2, a3, Lbc, Coord, Latcons, smin, smax, snum, lambda) ;
	    find_pair_cheby(R13, Tn13, Tnd13, rlen13, exprlen13, xdiff13,
			    a1, a3, Lbc, Coord, Latcons, smin, smax, snum, lambda) ;

	    for ( int i = 0 ; i < snum[ipair12] ; i++ ) 
	      for ( int j = 0 ; j < snum[ipair23] ; j++ ) 
		for ( int k = 0 ; k < snum[ipair13] ; k++ ) 
		  {
		    double coeff = Cheby_3B_Coeff(a1, a2, a3, i, j, k, Lbc, params, snum) ;

		    double fcut23, dfcut23 ;
		    double fcut13, dfcut13 ;

		    
		    cubic_cutoff(fcut23, dfcut23, rlen23, smax[ipair23]) ;
		    if ( fcut23 == 0.0 ) 
		      continue ;

		    cubic_cutoff(fcut13, dfcut13, rlen13, smax[ipair13]) ;
		    if ( fcut13 == 0.0 ) 
		      continue ;

		    tempx += coeff * fcut12 * fcut23 * fcut13 * Tn12[i] * Tn23[j] * Tn13[k] ;
		    
		    double deriv12 = 
		      fcut12 * Tnd12[i] *(-exprlen12/lambda[ipair12])/xdiff12 +
		      dfcut12 * Tn12[i] ;

		    double deriv23 =
		      fcut23 * Tnd23[j] *(-exprlen23/lambda[ipair23])/xdiff23 +
		      dfcut23 * Tn23[j] ;

		    double deriv13 =
		      fcut13 * Tnd13[k] *(-exprlen13/lambda[ipair13])/xdiff13 +
		      dfcut13 * Tn13[k] ;

		    for(int c=0;c<3;c++)
		      {
			SForce[a1][c] += coeff * deriv12 * fcut13 * fcut23 * Tn23[j] * Tn13[k] * R12[c] / rlen12 ;
			SForce[a2][c] -= coeff * deriv12 * fcut13 * fcut23 * Tn23[j] * Tn13[k] * R12[c] / rlen12 ;

			SForce[a2][c] += coeff * deriv23 * fcut12 * fcut13 * Tn12[i] * Tn13[k] * R23[c] / rlen23 ;
			SForce[a3][c] -= coeff * deriv23 * fcut12 * fcut13 * Tn12[i] * Tn13[k] * R23[c] / rlen23 ;

			SForce[a1][c] += coeff * deriv13 * fcut12 * fcut23 * Tn12[i] * Tn23[j] * R13[c] / rlen13 ;
			SForce[a3][c] -= coeff * deriv13 * fcut12 * fcut23 * Tn12[i] * Tn23[j] * R13[c] / rlen13 ;
		      } 
		    // TO DO:  Update pressure.
		  }
	  }
      }
  Vtot += tempx ;

  return;
}

static void find_pair_cheby(double Rab[3], double *Tn, double *Tnd, double &rlen, double &exprlen,
			    double &xdiff,
			    int a1, int a2, const char *Lbc, double **Coord, const double *Latcons, 
			    const double *smin, const double *smax, 
			    const int *snum, const double *lambda) 
// Calculate Chebyshev pair interaction parameters between two atoms.
// Rab: The displacement vector between the atoms
// Tn:  The Chebyshev polynomials
// Tnd: Derivatives of the Chebyshev polynomials
// rlen: Distance between the atos
// exprlen: Morse exponential function
// xdiff: Chebyshev limit variable.
// a1: Atom index 1
// a2: Atom index 2
// Lbc: Atom type (character).
// Coord: Coordinates of the atoms.
// Latcons: Orthorhombic box lattice constants
// smin: Minimum distance for Chebyshev
// smax: Maximum distance for Chebyshev
// snum: Chebyshev order
// lambda: Morse distance scale
{
  double Rvec[3];
  double x ;
  double xavg ;

  int ipair = pair_index(a1, a2, Lbc) ;

  // Start with minimum image convention.  Use layers to access larger distances if desired.	
  for(int c=0;c<3;c++) {
    Rvec[c]=Coord[a2][c]-Coord[a1][c];
    Rvec[c] -= floor(0.5+Rvec[c]/Latcons[c])*Latcons[c];
  }

  Rab[0]=Rvec[0] ;
  Rab[1]=Rvec[1] ;
  Rab[2]=Rvec[2] ;
		
  rlen=sqrt(Rab[0]*Rab[0]+Rab[1]*Rab[1]+Rab[2]*Rab[2]);
		
  double xmin = exp(-smax[ipair]/lambda[ipair]) ;
  double xmax = exp(-smin[ipair]/lambda[ipair]) ;
  xavg = 0.5 * (xmin + xmax) ;
  xdiff = 0.5 * (xmax - xmin) ;
  exprlen = exp(-rlen/lambda[ipair]) ;

  if( rlen < smax[ipair]) 
    {
      // Use morse variables.
      x = (exprlen-xavg)/xdiff ;
	      
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
	      
    } 
  else 
    {
      // Cut off the interaction.
      exprlen = exp(-smax[ipair]/lambda[ipair]) ;
      xdiff = 0.5 * (xmax - xmin) ;
      for ( int i = 0 ; i <= snum[ipair] ; i++ ) {
	Tn[i] = 0.0 ;
	Tnd[i] = 0.0 ;
      }
    }
}


  
static double Cheby_3B_Coeff(int a1, int a2, int a3, int n12, int n23, int n13,
			     const char *Lbc, const double *params,
			     const int *snum) 
// Extract the desired 3-body chebyshev coefficient from the params array.
{
  int ipair12, ipair23, ipair13 ;
  //int ele1, ele2, ele3 ;
  int index ;

  if ( a1 == a2 || a2 == a3 || a1 == a3 ) {
    cout << "Error in Cheby_3B_Coeff: atom indices are identical" << endl ;
    exit(1) ;
  }
    
  ipair12 = pair_index(a1, a2, Lbc) ;
  ipair23 = pair_index(a2, a3, Lbc) ;
  ipair13 = pair_index(a1, a3, Lbc) ;

  //ele1    = atom_index(a1, Lbc) ;
  //ele2    = atom_index(a2, Lbc) ;
  //ele3    = atom_index(a3, Lbc) ;


  // Move past pair interactions.
  index = 0 ;
  for ( int j = 0 ; j < NPAIR ; j++ ) {
    index += snum[j] ;
  }
  
  // Sort the pairs in ascending order.
  sort_triple(ipair12, ipair23, ipair13, n12, n23, n13) ;

  // Move through all of the pairs sets not used.
  for ( int i12 = 0 ; i12 < ipair12 ; i12++ ) {
    for ( int i23 = i12 ; i23 < ipair23 ; i23++ ) {
      for ( int i13 = i23 ; i13 < ipair13 ; i13++ ) {
	index += snum[i12] * snum[i23] * snum[i13] ;
      }
    }
  }

  // Increment the index within used pair set.
  index += n12 * snum[ipair23] * snum[ipair13] +
    n23 * snum[ipair13] + n13 ;

  return(params[index]) ;
}


void sort_triple(int &ipair12, int &ipair23, int &ipair13, 
		 int &n12, int &n23, int &n13)
// Sort three keys - ipair12, ipair23, ipair13, and corresponding values
// n12, n23, n13, in ascending order.
{
  int t1, t2 ;

  int key[3] = {ipair12, ipair23, ipair13} ;
  int val[3] = {n12, n23, n13} ;

  // Simple bubble sort
  while ( 1 ) {
    int j ;
    for ( j = 0 ; j < 2 ; j++ ) {
      if ( key[j] > key[j+1] ) {
	t1 = key[j] ;
	t2 = val[j] ;

	key[j] = key[j+1] ;
	val[j] = val[j+1] ;

	key[j+1] = t1 ;
	val[j+1] = t2 ;
	break ;
      }
    }
    if ( j == 2 ) {
      break ;
    }
  }

  ipair12 = key[0] ;
  ipair23 = key[1] ;
  ipair13 = key[2] ;

  n12 = val[0] ;
  n23 = val[1] ;
  n13 = val[2] ;

  if ( ipair12 > ipair23 || ipair23 > ipair13 ) {
    printf("Error: sort_triple failed\n") ;
    exit(1) ;
  }
}

static void cubic_cutoff(double &fcut, double &dfcut, double rlen, double smax) 
// Calculate a cubic force cutoff function.
{
  
  if ( rlen < smax ) {
    double fcut0 = (1.0 - rlen/smax) ;
    fcut = fcut0 * fcut0 * fcut0 ;
    dfcut = -3.0 * fcut0 * fcut0 / smax ;
  } else {
    fcut = 0.0 ;
    dfcut = 0.0 ;
  }
}

