#include "functions.h"
//
// 3-Body Chebyshev morse interaction.
//

static void find_pair_cheby(double Rab[3], double *Tn, double *Tnd, double &rlen, double &exprlen,
			    double &xdiff,
			    int a1, int a2, const char *Lbc, double **Coord, const double *Latcons, 
			    const double *smin, const double *smax, 
			    const int *snum, const double *lambda)  ;
static double Cheby_3B_Coeff(int a1, int a2, int a3, int n12, int n23, int n13,
			     const char *Lbc, const double *params,
			     const int *snum, 
			     const int *snum_3b_cheby,
			     int &index, bool return_params) ;
void sort_triple(int &ipair12, int &ipair23, int &ipair13, 
		 int &n12, int &n23, int &n13) ;
static void cubic_cutoff(double &fcut, double &dfcut, double rlen, double smax) ;

void ZCalc_3B_Cheby(double **Coord,const char *Lbc, double *Latcons,
		    const int nat,const double *smin,
		    const double *smax,
  		    const int *snum_3b_cheby,
		    double ******idx_params, const double *lambda,
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
	if ( snum_3b_cheby[i] > dim ) 
	  {
	    dim = snum_3b_cheby[i] ;
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
    {
      int ele1 = atom_index(a1,Lbc) ;
      for(int a2=a1+1;a2<nat;a2++) 
	{
	  double fcut12, dfcut12 ;
	  int ipair12 = pair_index(a1,a2,Lbc) ;
	  int ele2 = atom_index(a2,Lbc) ;
	  find_pair_cheby(R12,Tn12, Tnd12, rlen12, exprlen12, xdiff12,
			  a1, a2, Lbc, Coord, Latcons, smin, smax, snum_3b_cheby, 
			  lambda) ;
	
	  cubic_cutoff(fcut12, dfcut12, rlen12, smax[ipair12]) ;
	  if ( fcut12 == 0.0 ) 
	    continue ;

	  for(int a3=a2+1;a3<nat;a3++)
	    {
	      int ipair23 = pair_index(a2,a3,Lbc) ;
	      int ipair13 = pair_index(a1,a3,Lbc) ;
	      double fcut23, dfcut23 ;
	      double fcut13, dfcut13 ;
	      int ele3 = atom_index(a3, Lbc) ;

	      find_pair_cheby(R23, Tn23, Tnd23, rlen23, exprlen23, xdiff23,
			      a2, a3, Lbc, Coord, Latcons, smin, smax, snum_3b_cheby, 
			      lambda) ;
	      find_pair_cheby(R13, Tn13, Tnd13, rlen13, exprlen13, xdiff13,
			      a1, a3, Lbc, Coord, Latcons, smin, smax, snum_3b_cheby, 
			      lambda) ;

		    
	      cubic_cutoff(fcut23, dfcut23, rlen23, smax[ipair23]) ;
	      if ( fcut23 == 0.0 ) 
		continue ;

	      cubic_cutoff(fcut13, dfcut13, rlen13, smax[ipair13]) ;
	      if ( fcut13 == 0.0 ) 
		continue ;

	      for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
		for ( int j = 0 ; j < snum_3b_cheby[ipair23] ; j++ ) 
		  for ( int k = 0 ; k < snum_3b_cheby[ipair13] ; k++ ) 
		    {
		      double coeff = idx_params[ele1][ele2][ele3][i][j][k] ;

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
    }
  Vtot += tempx ;

  return;
}


double ******Indexed_3B_Cheby_Coeffs(const char *Lbc, 
				  const int nat,
				  const int *snum, 
				  const int *snum_3b_cheby,
				  double *params)

// Generate a "flattened" list of 3-Body Chebyshev Polynomial coefficients to speed up
// indexing in generating 3 body forces.
{
  //spline term calculated w/cutoff:
  int index = 0 ;


  int dim = 0 ;
    for ( int i = 0 ; i < NPAIR ; i++ ) 
      {
	if ( snum_3b_cheby[i] > dim ) 
	  {
	    dim = snum_3b_cheby[i] ;
	  }
      }
    dim++ ;

  double ******flat_params ;
  flat_params = new double*****[NELE] ;
  for ( int i1 = 0 ; i1 < NELE ; i1++ ) {
    flat_params[i1] = new double****[NELE] ;
    for ( int i2 = 0 ; i2 < NELE ; i2++ ) {
      flat_params[i1][i2] = new double***[NELE] ;
      for ( int i3 = 0 ; i3 < NELE ; i3++ ) {
	flat_params[i1][i2][i3] = new double **[dim] ;
	for ( int i4 = 0 ; i4 < dim ; i4++ ) {
	  flat_params[i1][i2][i3][i4] = new double *[dim] ;
	  for ( int i5 = 0 ; i5 < dim ; i5++ ) {
	    flat_params[i1][i2][i3][i4][i5] = new double [dim] ;
	  }
	}
      }
    }
  }

  ////main loop for 3-body Chebyshev terms:
  for(int a1=0;a1<nat-1;a1++) 
    {
      int ele1 = atom_index(a1,Lbc) ;
      for(int a2=a1+1;a2<nat;a2++) 
	{
	  int ele2 = atom_index(a2,Lbc) ;
	  int ipair12 = pair_index(a1,a2,Lbc) ;
	  for(int a3=a2+1;a3<nat;a3++)
	    {
	      int ipair23 = pair_index(a2,a3,Lbc) ;
	      int ipair13 = pair_index(a1,a3,Lbc) ;
	      int ele3 = atom_index(a3,Lbc) ;

	      for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
		for ( int j = 0 ; j < snum_3b_cheby[ipair23] ; j++ ) 
		  for ( int k = 0 ; k < snum_3b_cheby[ipair13] ; k++ ) 
		    {
		      (void) Cheby_3B_Coeff(a1, a2, a3, i, j, k, Lbc, params, 
					    snum, snum_3b_cheby, index, false) ;
		      flat_params[ele1][ele2][ele3][i][j][k] = params[index] ;
		    }
	    }
	}
    }
  return flat_params ;
}



void ZCalc_3B_Cheby_Deriv(double **Coord,const char *Lbc, double *Latcons,
			  const int nat, double ***A,
			  const double *smin,
			  const double *smax,
			  const int *snum, 
			  const int *snum_3b_cheby,
			  const double *lambda)
// Calculate derivatives of the 3-body short-range forces using a Chebyshev polynomial expansion,
//  with respect to multiplicative parameters.
{
  double R12[3], R23[3], R13[3];
  double rlen12, rlen23, rlen13 ;
  static double *Tn12, *Tnd12 ;
  static double *Tn23, *Tnd23 ;
  static double *Tn13, *Tnd13 ;

  static bool called_before = false ;

  double exprlen12, exprlen23, exprlen13 ;
  double xdiff12, xdiff23, xdiff13 ;
  int index = 0 ;

  if ( ! called_before ) {
    called_before = true ;
    int dim = 0 ;
    for ( int i = 0 ; i < NPAIR ; i++ ) 
      {
	if ( snum_3b_cheby[i] > dim ) 
	  {
	    dim = snum_3b_cheby[i] ;
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
			a1, a2, Lbc, Coord, Latcons, smin, smax, 
			snum_3b_cheby,
			lambda) ;
	
	cubic_cutoff(fcut12, dfcut12, rlen12, smax[ipair12]) ;
	if ( fcut12 == 0.0 ) 
	  continue ;

	for(int a3=a2+1;a3<nat;a3++)
	  {
	    int ipair23 = pair_index(a2,a3,Lbc) ;
	    int ipair13 = pair_index(a1,a3,Lbc) ;
	    double fcut23, dfcut23 ;
	    double fcut13, dfcut13 ;

	    find_pair_cheby(R23, Tn23, Tnd23, rlen23, exprlen23, xdiff23,
			    a2, a3, Lbc, Coord, Latcons, smin, smax, 
			    snum_3b_cheby,
			    lambda) ;
	    find_pair_cheby(R13, Tn13, Tnd13, rlen13, exprlen13, xdiff13,
			    a1, a3, Lbc, Coord, Latcons, smin, smax, 
			    snum_3b_cheby,
			    lambda) ;

		    
	    cubic_cutoff(fcut23, dfcut23, rlen23, smax[ipair23]) ;
	    if ( fcut23 == 0.0 ) 
	      continue ;

	    cubic_cutoff(fcut13, dfcut13, rlen13, smax[ipair13]) ;
	    if ( fcut13 == 0.0 ) 
	      continue ;

	    for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
	      for ( int j = 0 ; j < snum_3b_cheby[ipair23] ; j++ ) 
		for ( int k = 0 ; k < snum_3b_cheby[ipair13] ; k++ ) 
		  {
		    (void) Cheby_3B_Coeff(a1, a2, a3, i, j, k, Lbc, NULL, snum, 
					  snum_3b_cheby,
					  index, false) ;


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
			A[a1][index][c] += deriv12 * fcut13 * fcut23 * Tn23[j] * Tn13[k] * R12[c] / rlen12 ;
			A[a2][index][c] -= deriv12 * fcut13 * fcut23 * Tn23[j] * Tn13[k] * R12[c] / rlen12 ;

			A[a2][index][c] += deriv23 * fcut12 * fcut13 * Tn12[i] * Tn13[k] * R23[c] / rlen23 ;
			A[a3][index][c] -= deriv23 * fcut12 * fcut13 * Tn12[i] * Tn13[k] * R23[c] / rlen23 ;

			A[a1][index][c] += deriv13 * fcut12 * fcut23 * Tn12[i] * Tn23[j] * R13[c] / rlen13 ;
			A[a3][index][c] -= deriv13 * fcut12 * fcut23 * Tn12[i] * Tn23[j] * R13[c] / rlen13 ;
		      } 
		  }
	  }
      }

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
// rlen: Distance between the atoms
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
  xdiff = 0.5 * (xmax - xmin) ;

  if( rlen < smax[ipair] && rlen > smin[ipair] ) 
    {
      // Use morse variables.
      xavg = 0.5 * (xmin + xmax) ;
      exprlen = exp(-rlen/lambda[ipair]) ;
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
      exprlen = 1.0 ;
      xdiff = 0.5 * (xmax - xmin) ;
      for ( int i = 0 ; i <= snum[ipair] ; i++ ) {
	Tn[i] = 0.0 ;
	Tnd[i] = 0.0 ;
      }
    }
}


  
static double Cheby_3B_Coeff(int a1, int a2, int a3, int n12, int n23, int n13,
			     const char *Lbc, const double *params,
			     const int *snum, 
			     const int *snum_3b_cheby,
			     int &index, bool return_params) 
// Extract the desired 3-body chebyshev coefficient from the params array.
{
  int ipair12, ipair23, ipair13 ;
  //int ele1, ele2, ele3 ;

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
	index += snum_3b_cheby[i12] * snum_3b_cheby[i23] * snum_3b_cheby[i13] ;
      }
    }
  }

  // Increment the index within used pair set.
  index += n12 * snum_3b_cheby[ipair23] * snum_3b_cheby[ipair13] +
    n23 * snum_3b_cheby[ipair13] + n13 ;

  if ( return_params ) {
    return(params[index]) ;
  } else {
    return(0.0) ;
  }
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
// fcut is the cutoff function.
// dfcut is the first derivative of the cutoff function.
// rlen is the interparticle distance
// smax is the cutoff length.
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



int count_cheby_3b_params(const int *snum)
{
  int num_cheby_3b = 0 ;
  for ( int i12 = 0 ; i12 < NPAIR ; i12++ ) {
    for ( int i23 = i12 ; i23 < NPAIR ; i23++ ) {
      for ( int i13 = i23 ; i13 < NPAIR ; i13++ ) {
	int val = snum[i12] * snum[i23] * snum[i13] ;
	//	printf("3B Params for %d %d %d = %d\n", i12, i23, i13, val) ;
	num_cheby_3b += val ;
      }
    }
  }
  return(num_cheby_3b) ;
}
