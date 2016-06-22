#include "functions.h"
//
// 3-Body Chebyshev morse interaction.
//

/* #define TESTING */

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
static int match_pair_order(int a1, int a2, int a_orig[3], int n12, int n13, int n23) ;
static void sort_pair(int &n, int &m)  ;
static void sort_poly_order(int ipair12, int ipair13, int ipair23, int &n12, int &n13, int &n23) ;

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
		for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
		  for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 

		    {
		      // Require 2 non-zero Cheby orders to create a 3-body interaction.
		      if ( i + j + k < 2 ) 
			continue ;

		      double coeff = idx_params[ele1][ele2][ele3][i][j][k] ;

		      if ( coeff == EMPTY ) {
			printf("Error: empty 3-body parameter found\n") ;
			exit(1) ;
		      }

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

  double ******store_3b_params ;
  store_3b_params = new double*****[NELE] ;
  for ( int i1 = 0 ; i1 < NELE ; i1++ ) {
    store_3b_params[i1] = new double****[NELE] ;
    for ( int i2 = 0 ; i2 < NELE ; i2++ ) {
      store_3b_params[i1][i2] = new double***[NELE] ;
      for ( int i3 = 0 ; i3 < NELE ; i3++ ) {
	store_3b_params[i1][i2][i3] = new double **[dim] ;
	for ( int i4 = 0 ; i4 < dim ; i4++ ) {
	  store_3b_params[i1][i2][i3][i4] = new double *[dim] ;
	  for ( int i5 = 0 ; i5 < dim ; i5++ ) {
	    store_3b_params[i1][i2][i3][i4][i5] = new double [dim] ;
	    for ( int i6 = 0 ; i6 < dim ; i6++ ) {
	      store_3b_params[i1][i2][i3][i4][i5][i6] = EMPTY ;
	    }
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
		for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
		  for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 
		    {
		      if ( i + j + k < 2 ) {
			store_3b_params[ele1][ele2][ele3][i][j][k] = 0.0 ;
		      } else {
			(void) Cheby_3B_Coeff(a1, a2, a3, i, j, k, Lbc, params, 
					      snum, snum_3b_cheby, index, false) ;
			store_3b_params[ele1][ele2][ele3][i][j][k] = params[index] ;
		      }
		    }
	    }
	}
    }

  // Print out parameters.
  printf("3-Body Cheby Parameters\n") ;
  for(int ele1=0; ele1 < NELE ; ele1++) 
    {
      for(int ele2=ele1; ele2 < NELE ; ele2++) 
	{
	  int ipair12 = pair_index_ele(ele1, ele2) ;
	  for(int ele3 = ele2 ; ele3 < NELE ; ele3++)
	    {
	      int ipair23 = pair_index_ele(ele2,ele3) ;
	      int ipair13 = pair_index_ele(ele1,ele3) ;

	      for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
		for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
		  for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 
		    {
		      if ( i + j + k >= 2 ) {
			printf("ele: %d %d %d pow: %d %d %d param: %11.4e\n",
			       ele1, ele2, ele3, i, j, k, 
			       store_3b_params[ele1][ele2][ele3][i][j][k]) ;
		      }
		    }
	    }
	}
    }
  return store_3b_params ;
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
	      for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
		for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 
		  {
		    // Require 2 non-zero Cheby orders to create a 3-body interaction.
		    if ( i + j + k < 2 ) 
		      continue ;

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

  if( rlen < smax[ipair] ) 
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


  
static double Cheby_3B_Coeff(int a1, int a2, int a3, int n12, int n13, int n23,
			     const char *Lbc, const double *params,
			     const int *snum, 
			     const int *snum_3b_cheby, 
			     int &index, bool return_params) 
// Extract the desired 3-body chebyshev coefficient from the params array.
{
  int ipair12, ipair23, ipair13 ;
  int ipair12_orig, ipair23_orig, ipair13_orig ;
  int ele1, ele2, ele3 ;
  int tot_short_range ;
  int tot_2b ;

  if ( a1 == a2 || a2 == a3 || a1 == a3 ) {
    cout << "Error in Cheby_3B_Coeff: atom indices are identical" << endl ;
    exit(1) ;
  }
  if ( n12 + n13 + n23 < 2 ) {
    // Not 3-body.
    cout << "Error in Cheby_3B_Coeff: sum of polynomial powers are too small.\n" ;
  }
    
  ele1    = atom_index(a1, Lbc) ;
  ele2    = atom_index(a2, Lbc) ;
  ele3    = atom_index(a3, Lbc) ;

  ipair12_orig = pair_index_ele(ele1, ele2) ;
  ipair23_orig = pair_index_ele(ele2, ele3) ;
  ipair13_orig = pair_index_ele(ele1, ele3) ;

#ifdef TESTING  
  printf("Input:\n  ele1 = %d ele2 = %d ele3 = %d\n  n12 = %d n13 = %d n23 = %d\n",
	 ele1, ele2, ele3, n12, n13, n23) ;
  printf("  a1 = %d a2 = %d a3 = %d\n", a1, a2, a3) ;
  printf("  ipair12 = %d ipair13 = %d ipair23 = %d\n", ipair12_orig, ipair13_orig, ipair23_orig) ;
#endif

  int ele1_orig, ele2_orig, ele3_orig ;
  int j1, j2, j3 ;

  j1 = j2 = j3 = 0 ;
  ele1_orig = ele1 ;
  ele2_orig = ele2 ;
  ele3_orig = ele3 ;
  
  // Sort element indices in ascending order.
  int a_orig[3] = {a1, a2, a3} ;
  sort_triple(ele1, ele2, ele3, a1, a2, a3) ;

  // Move past pair interactions.
  index = 0 ;
  for ( int j = 0 ; j < NPAIR ; j++ ) {
    index += snum[j] ;
  }
  tot_2b = index ;
  tot_short_range = index ;
  tot_short_range += count_cheby_3b_params(snum_3b_cheby) ;

  // Sort the pairs in ascending order.


  // Move through all of the pairs sets not used.
  for ( int i1 = 0 ; i1 <= ele1 ; i1++ ) {
    for ( int i2 = i1 ; i2 <= ele2 ; i2++ ) {
      for ( int i3 = i2 ; i3 <= ele3 ; i3++ ) {
	int i12 = pair_index_ele(i1, i2) ;
    	int i23 = pair_index_ele(i2, i3) ;
	int i13 = pair_index_ele(i1, i3) ;
	if ( i1 == ele1 && i2 == ele2 && i3 == ele3 ) break ;
#ifdef TESTING
	printf("Skipping over element triple %d %d %d\n", i1, i2, i3) ;
	printf(" Pair index : %d %d %d\n", i12, i13, i23) ;
#endif
	index += snum_3b_cheby[i12] * snum_3b_cheby[i23] * snum_3b_cheby[i13] ;

	// Decrement for cases when sum of Cheby orders < 2 : not 3 body.
	index -= 4 ;
      }
    }
  }

  ipair12 = pair_index_ele(ele1, ele2) ;
  ipair23 = pair_index_ele(ele2, ele3) ;
  ipair13 = pair_index_ele(ele1, ele3) ;

  int n12_new = match_pair_order(a1, a2, a_orig, n12, n13, n23) ;
  int n13_new = match_pair_order(a1, a3, a_orig, n12, n13, n23) ;
  int n23_new = match_pair_order(a2, a3, a_orig, n12, n13, n23) ;

#ifdef TESTING
  printf("Sorted:\n  ele1 = %d ele2 = %d ele3 = %d\n   n12 = %d n13 = %d n23 = %d\n",
	 ele1, ele2, ele3, n12_new, n13_new, n23_new) ;
  printf("  a1 = %d a2 = %d a3 = %d\n", a1, a2, a3) ;
  printf("  ipair12 = %d ipair13 = %d ipair23 = %d\n\n", ipair12, ipair13, ipair23) ;
#endif

  // If the pair index is the same, the polynomial order must be sorted.  This keeps
  // the polynomial invariant to ordering of atoms.
  sort_poly_order(ipair12, ipair13, ipair23, n12_new, n13_new, n23_new) ;

  // Test for special cases.
  if ( n12_new == 0 && n13_new == 1 && n23_new == 1 ) {
    // No additional increment to index needed.
    ;
  } else {
    // Increment the index within used pair set.
    // Decrement by 4 to account for unused indices:
    // These indices correspond to constant or 2-body interactions.
    // 0 0 0 (NOT USED)
    // 0 0 1 (NOT USED)
    // 0 1 1 (USED)
    // 0 1 0 (NOT USED)
    // 1 0 0 (NOT USED)
    // All other parameters are used.
    
    index += n12_new * snum_3b_cheby[ipair13] * snum_3b_cheby[ipair23] +
      n13_new * snum_3b_cheby[ipair23] + n23_new - 4 ;
  }

  assert(index >= tot_2b && index < tot_short_range ) ;

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
// Input: 
//   rlen is the interparticle distance
//   smax is the cutoff length.
// Output:
//   fcut is the cutoff function.
//   dfcut is the first derivative of the cutoff function.
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
// Returns the total number of 3 body Chebyshev parameters.
{
  int num_cheby_3b = 0 ;

  for ( int i1 = 0 ; i1 < NELE ; i1++ ) {
    for ( int i2 = 0 ; i2 <= i1 ; i2++ ) {
      int i12 = pair_index_ele(i1, i2) ;
      for ( int i3 = 0 ; i3 <= i2 ; i3++ ) {
	int i23 = pair_index_ele(i2, i3) ;
	int i13 = pair_index_ele(i1, i3) ;

	num_cheby_3b += snum[i12] * snum[i23] * snum[i13] ;

	// Decrement for cases when sum of Cheby orders < 2 : not 3 body.
	num_cheby_3b -= 4 ;
      }
    }
  }

  return(num_cheby_3b) ;
}

int pair_index_ele(int ele1, int ele2)
// Returns the index of a particular element pair. Element pairs are listed
// in a notional ordered list:  For 3 elements, this would be 
// 00, 01, 02, 11, 12, 22.  This routine picks out the index of a particular pair
// in the list.  The element pairs do not depend on order of the input elements.
{
  if ( ele1 > ele2 ) {
    int tmp = ele1 ;
    ele1 = ele2 ;
    ele2 = tmp ;
  }
  int idx = 0 ;
  int i ;
  for ( i = 0 ; i < ele1 ; i++ ) {
    idx += (NELE-i) ;
  }
  idx += ele2 - ele1 ;
  assert(idx >= 0 && idx < NELE * (NELE+1)/2) ;
  return(idx) ;
}

static int match_pair_order(int a1, int a2, int a_orig[3], int n12, int n13, int n23)
// Return the Chebyshev order (n in T[n]) of the polynomial corresponding to permuted atom
// indices a1 and a2 given original polynomial pair orders n12, n13, and n23 and original
// atom indices a_orig
{
  int i1, i2 ;

  if ( a1 == a_orig[0] ) {
    i1 = 1 ;
  } else if ( a1 == a_orig[1] ) {
    i1 = 2 ;
  } else if ( a1 == a_orig[2] ) {
    i1 = 3 ;
  } else {
    printf("Error in match_pair_order: no match found\n") ;
    exit(1) ;
  }

  if ( a2 == a_orig[0] ) {
    i2 = 1 ;
  } else if ( a2 == a_orig[1] ) {
    i2 = 2 ;
  } else if ( a2 == a_orig[2] ) {
    i2 = 3 ;
  } else {
    printf("Error in match_pair_order: no match found\n") ;
    exit(1) ;
  }

  if ( i1 == 1 ) {
    if ( i2 == 2 ) {
      return(n12) ;
    } else if ( i2 == 3 ) {
      return(n13) ;
    } else {
      printf("Error in match_pair_order\n") ;
      exit(1) ;
    }
  }
  else if ( i1 == 2 ) {
    if ( i2 == 1 ) {
      return(n12) ;
    } else if ( i2 == 3 ) {
      return(n23) ;
    } else {
      printf("Error in match_pair_order\n") ;
      exit(1) ;
    }
  } else if ( i1 == 3 ) {
    if ( i2 == 1 ) {
      return(n13) ;
    } else if ( i2 == 2 ) {
      return(n23) ;
    } else {
      printf("Error in match_pair_order\n") ;
      exit(1) ;
    }
  } else {
    printf("Error in match_pair_order\n") ;
    exit(1) ;
  }

  printf("Error in match_pair_order\n") ;
  exit(1) ;
}

static void sort_poly_order(int ipair12, int ipair13, int ipair23, int &n12, int &n13, int &n23)
// Sort the polynomial in increasing order if the corresponding pair indices are the same.
// This maintains invariance of the polynomial with respect to ordering of coefficients.
{
  if ( ipair12 == ipair13 && ipair12 == ipair23 ) {
    int i, j, k ;

    // Dummy arg to sort_triple.
    i = j = k = 0 ;
    sort_triple(n12, n13, n23, i, j, k) ;
  } else if ( ipair12 == ipair13 ) {
    sort_pair(n12, n13) ;
  } else if ( ipair12 == ipair23 ) {
    sort_pair(n12, n23) ;
  } else if ( ipair13 == ipair23 ) {
    sort_pair(n13, n23) ;
  }
}

static void sort_pair(int &n, int &m) 
// Sort a pair of numbers into increasing order.
{
  if ( n > m ) {
    int tmp = n ;
    n = m ;
    m = tmp ;
  }
}

