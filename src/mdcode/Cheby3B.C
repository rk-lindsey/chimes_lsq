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
static int check_interaction_ordering(int i, int j, int k, int i12, int i13, int i23) ;
static bool is_three_body(int i, int j, int k) ;

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
  double f3b[nat][3];
  bool fcheck = true;

  if (fcheck) {
    memset(f3b, 0, sizeof(f3b));
  }
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
  for(int a1=0; a1<nat; a1++) 
    {
      int ele1 = atom_index(a1,Lbc) ;
      for(int a2=a1+1; a2<nat ; a2++) 
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

	  for(int a3=a2+1; a3<nat ; a3++)
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

              double r12_mag, r23_mag, r13_mag;
              r12_mag = sqrt(R12[0]*R12[0] + R12[1]*R12[1] + R12[2]*R12[2]);
              r23_mag = sqrt(R23[0]*R23[0] + R23[1]*R23[1] + R23[2]*R23[2]);
              r13_mag = sqrt(R13[0]*R13[0] + R13[1]*R13[1] + R13[2]*R13[2]);

	      for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
		for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
		  for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 

		    {
		      // Require 2 non-zero Cheby orders to create a 3-body interaction.
		      if ( ! is_three_body(i,j,k) ) 
			continue ;

		      double coeff = idx_params[ele1][ele2][ele3][i][j][k] ;

		      if ( coeff == EMPTY ) {
			printf("Error: empty 3-body parameter found\n") ;
			exit(1) ;
		      }

		      tempx += coeff * fcut12 * fcut23 * fcut13 * Tn12[i] * Tn23[k] * Tn13[j] ;
		    
		      double deriv12 = 
			fcut12 * Tnd12[i] *(-exprlen12/lambda[ipair12])/xdiff12 +
			dfcut12 * Tn12[i] ;

		      double deriv13 =
			fcut13 * Tnd13[j] *(-exprlen13/lambda[ipair13])/xdiff13 +
			dfcut13 * Tn13[j] ;

		      double deriv23 =
			fcut23 * Tnd23[k] *(-exprlen23/lambda[ipair23])/xdiff23 +
			dfcut23 * Tn23[k] ;
		      
		      double f12 = coeff * deriv12 * fcut13 * fcut23 * Tn23[k] * Tn13[j] ;
		      Pxyz -= f12 * rlen12 ;
		      f12 /= rlen12 ;

		      double f23 = coeff * deriv23 * fcut12 * fcut13 * Tn12[i] * Tn13[j] ;
		      Pxyz -= f23 * rlen23 ;
		      f23 /= rlen23 ;

		      double f13 = coeff * deriv13 * fcut12 * fcut23 * Tn12[i] * Tn23[k] ;
		      Pxyz -= f13 * rlen13 ;
		      f13 /= rlen13 ;

                      // zero-out force if rij < rmin
                      if (r12_mag < smin[ipair12]) {
                        f12 = 0;
                        //printf("Zeroing out f12: r12 = %lf, smin = %lf\n",r12_mag,smin[ipair12]);
                      }
                      if (r23_mag < smin[ipair23]) {
                        f23 = 0;
                        //printf("Zeroing out f23: r23 = %lf, smin = %lf\n",r23_mag,smin[ipair23]);
                      }
                      if (r13_mag < smin[ipair13]) {
                        f13 = 0;
                        //printf("Zeroing out f13: r13 = %lf, smin = %lf\n",r13_mag,smin[ipair13]);
                      }
		      for(int c=0;c<3;c++)
			{
			  SForce[a1][c] += f12 * R12[c] ;
			  SForce[a2][c] -= f12 * R12[c] ; 

			  SForce[a2][c] += f23 * R23[c] ;
			  SForce[a3][c] -= f23 * R23[c] ;

			  SForce[a1][c] += f13 * R13[c] ;
			  SForce[a3][c] -= f13 * R13[c] ;
                          if (fcheck) {
			    f3b[a1][c] += f12 * R12[c] ;
			    f3b[a2][c] -= f12 * R12[c] ; 
			    f3b[a2][c] += f23 * R23[c] ;
			    f3b[a3][c] -= f23 * R23[c] ;
			    f3b[a1][c] += f13 * R13[c] ;
			    f3b[a3][c] -= f13 * R13[c] ;
                          }

			} 
		    }
	    }
	}
    }
  if (fcheck) {
    double e3b = tempx;
    FILE *frs;
    if (called_before) {
      frs = fopen("3b_results.dat","w");
    } else { 
      frs = fopen("3b_results.dat","a");
    }
    fprintf(frs, "e3b = %16.16lf\n",e3b);
    for (int i = 0; i < nat; i++) {
      for (int j = 0; j < 3; j++) {
        fprintf(frs, "%16.16lf\n",f3b[i][j]);
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
				     double *params
				     ) 
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

  // Find an atom for each element type.
  // Need at least 3 examples.
  int *example_atm = new int [ 3 * NELE ] ;
  for ( int j = 0 ; j < NELE ; j++ ) {
    int a1 ;
    int k = 0 ;
    for( a1=0 ; a1 < nat ; a1++) {
      if ( atom_index(a1,Lbc) == j ) {
	example_atm[3*j+k] = a1 ;
	k++ ;
	if ( k == 3 ) break ;
      }
    }
    if ( a1 == nat ) {
      printf("Error: 3 atoms with element index %d could not be found\n", j) ;
    }
  }

  for ( int ele1 = 0 ; ele1 < NELE ; ele1++ ) {
    int a1 = example_atm[3*ele1] ;
    for ( int ele2 = 0 ; ele2 < NELE ; ele2++ ) {
      int a2 = example_atm[3*ele2+1] ;
      int ipair12 = pair_index(a1,a2,Lbc) ;
      for ( int ele3 = 0 ; ele3 < NELE ; ele3++ ) {
	int a3 = example_atm[3*ele3+2] ;
	int ipair23 = pair_index(a2,a3,Lbc) ;
	int ipair13 = pair_index(a1,a3,Lbc) ;
	for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
	  for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
	    for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 
	      {
		if ( ! is_three_body(i,j,k) ) {
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
  printf("3-Body Chebyshev Polynomial Coefficients:\n") ;
  for(int ele1=0; ele1 < NELE ; ele1++) 
    {
      int a1 = example_atm[3*ele1] ;
      for(int ele2=0 ; ele2 < NELE ; ele2++) 
	{
	  int ipair12 = pair_index_ele(ele1, ele2) ;
	  int a2 = example_atm[3*ele2] ;
	  for(int ele3 = 0 ; ele3 < NELE ; ele3++)
	    {
	      int ipair23 = pair_index_ele(ele2,ele3) ;
	      int ipair13 = pair_index_ele(ele1,ele3) ;
	      int a3 = example_atm[3*ele3] ;

	      for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
		for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
		  for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 
		    {
		      if ( is_three_body(i,j,k) ) {
			printf("ele: %c %c %c pair index: %d %d %d: pow: %d %d %d param: %11.4e\n",
			       Lbc[a1], Lbc[a2], Lbc[a3], ipair12, ipair13, ipair23, i, j, k,
			       store_3b_params[ele1][ele2][ele3][i][j][k]) ;
		      }
		    }
	    }
	}
    }
  return store_3b_params ;
  }


int ******Index_3B_Cheby(const char *Lbc, 
			 const int nat,
			 const int *snum, 
			 const int *snum_3b_cheby)
			 
// Generate a precalculated list of 3-Body Chebyshev Polynomial indices to speed up
// indexing in generating 3 body forces force derivatives.
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

  int ******params_index ;
  params_index = new int*****[NELE] ;
  for ( int i1 = 0 ; i1 < NELE ; i1++ ) {
    params_index[i1] = new int****[NELE] ;
    for ( int i2 = 0 ; i2 < NELE ; i2++ ) {
      params_index[i1][i2] = new int***[NELE] ;
      for ( int i3 = 0 ; i3 < NELE ; i3++ ) {
	params_index[i1][i2][i3] = new int **[dim] ;
	for ( int i4 = 0 ; i4 < dim ; i4++ ) {
	  params_index[i1][i2][i3][i4] = new int *[dim] ;
	  for ( int i5 = 0 ; i5 < dim ; i5++ ) {
	    params_index[i1][i2][i3][i4][i5] = new int [dim] ;
	    for ( int i6 = 0 ; i6 < dim ; i6++ ) {
	      params_index[i1][i2][i3][i4][i5][i6] = -1 ;
	    }
	  }
	}
      }
    }
  }

  ////main loop for 3-body Chebyshev terms:

  // Find an atom for each element type.
  // Need at least 3 examples.
  int *example_atm = new int [ 3 * NELE ] ;
  for ( int j = 0 ; j < NELE ; j++ ) {
    int a1 ;
    int k = 0 ;
    for( a1=0 ; a1 < nat ; a1++) {
      if ( atom_index(a1,Lbc) == j ) {
	example_atm[3*j+k] = a1 ;
	k++ ;
	if ( k == 3 ) break ;
      }
    }
    if ( a1 == nat ) {
      printf("Error: 3 atoms with element index %d could not be found\n", j) ;
    }
  }

  for ( int ele1 = 0 ; ele1 < NELE ; ele1++ ) {
    int a1 = example_atm[3*ele1] ;
    for ( int ele2 = 0 ; ele2 < NELE ; ele2++ ) {
      int a2 = example_atm[3*ele2+1] ;
      int ipair12 = pair_index(a1,a2,Lbc) ;
      for ( int ele3 = 0 ; ele3 < NELE ; ele3++ ) {
	int a3 = example_atm[3*ele3+2] ;
	int ipair23 = pair_index(a2,a3,Lbc) ;
	int ipair13 = pair_index(a1,a3,Lbc) ;
	for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
	  for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
	    for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 
	      {
		if ( ! is_three_body(i,j,k) ) {
		  params_index[ele1][ele2][ele3][i][j][k] = -1 ;
		} else {
		  (void) Cheby_3B_Coeff(a1, a2, a3, i, j, k, Lbc, NULL, 
					snum, snum_3b_cheby, index, false) ;
		  params_index[ele1][ele2][ele3][i][j][k] = index ;
		}
	      }
      }
    }
  }
  // Print out parameters.
  printf("3-Body Cheby Parameter Indices\n") ;
  for(int ele1=0; ele1 < NELE ; ele1++) 
    {
      int a1 = example_atm[3*ele1] ;
      for(int ele2=0 ; ele2 < NELE ; ele2++) 
	{
	  int ipair12 = pair_index_ele(ele1, ele2) ;
	  int a2 = example_atm[3*ele2] ;
	  for(int ele3 = 0 ; ele3 < NELE ; ele3++)
	    {
	      int ipair23 = pair_index_ele(ele2,ele3) ;
	      int ipair13 = pair_index_ele(ele1,ele3) ;
	      int a3 = example_atm[3*ele3] ;

	      for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) 
		for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) 
		  for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) 
		    {
		      if ( is_three_body(i,j,k) ) {
			printf("ele: %c %c %c pair index: %d %d %d pow: %d %d %d parameter index: %d\n",
			       Lbc[a1], Lbc[a2], Lbc[a3], ipair12, ipair13, ipair23, i, j, k,
			       params_index[ele1][ele2][ele3][i][j][k]) ;
		      }
		    }
	    }
	}
    }
  return params_index ;
}




void ZCalc_3B_Cheby_Deriv(double **Coord,const char *Lbc, double *Latcons,
			  const int nat, double ***A,
			  const double *smin,
			  const double *smax,
			  const int *snum_3b_cheby,
			  const double *lambda, int ******param_index)
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
  for(int a1=0;a1<nat;a1++) {
    int ele1    = atom_index(a1, Lbc) ;
    for(int a2=a1+1;a2<nat;a2++) 
      {
	double fcut12, dfcut12 ;
	int ipair12 = pair_index(a1,a2,Lbc) ;
	int ele2 = atom_index(a2, Lbc) ;
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
	    int ele3 = atom_index(a3, Lbc) ;
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
		    if ( ! is_three_body(i,j,k) ) 
		      continue ;

		    index = param_index[ele1][ele2][ele3][i][j][k] ;
		    
		    //(void) Cheby_3B_Coeff(a1, a2, a3, i, j, k, Lbc, NULL, snum, 
		    // snum_3b_cheby, index, false) ;

		    double deriv12 = 
		      fcut12 * Tnd12[i] *(-exprlen12/lambda[ipair12])/xdiff12 +
		      dfcut12 * Tn12[i] ;

		    double deriv13 =
		      fcut13 * Tnd13[j] * (-exprlen13/lambda[ipair13])/xdiff13 +
		      dfcut13 * Tn13[j] ;

		    double deriv23 =
		      fcut23 * Tnd23[k] *(-exprlen23/lambda[ipair23])/xdiff23 +
		      dfcut23 * Tn23[k] ;

		    for(int c=0;c<3;c++)
		      {
			A[a1][index][c] += deriv12 * fcut13 * fcut23 * Tn23[k] * Tn13[j] * R12[c] / rlen12 ;
			A[a2][index][c] -= deriv12 * fcut13 * fcut23 * Tn23[k] * Tn13[j] * R12[c] / rlen12 ;

			A[a2][index][c] += deriv23 * fcut12 * fcut13 * Tn12[i] * Tn13[j] * R23[c] / rlen23 ;
			A[a3][index][c] -= deriv23 * fcut12 * fcut13 * Tn12[i] * Tn13[j] * R23[c] / rlen23 ;

			A[a1][index][c] += deriv13 * fcut12 * fcut23 * Tn12[i] * Tn23[k] * R13[c] / rlen13 ;
			A[a3][index][c] -= deriv13 * fcut12 * fcut23 * Tn12[i] * Tn23[k] * R13[c] / rlen13 ;
		      } 
#if(0)
		    cout << a1 << " " << a2 << " " << ipair12 
			 << " " << index << " " << deriv12 << " " << R12[0] / rlen12 << endl ;
		    cout << a1 << " " << a3 << " " << ipair13 
			 << " " << index << " " << deriv13 << " " << R13[0] / rlen13 << endl ;
		    cout << a2 << " " << a3 << " " << ipair23 
			 << " " << index << " " << deriv23 << " " << R23[0] / rlen23 << endl ;
		    cout << endl ;
#endif
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
  if ( ! is_three_body(n12,n13,n23) ) {
    // Not 3-body.
    cout << "Error in Cheby_3B_Coeff: polynomial powers do not correspond to a 3-body interaction\n" ;
    exit(1) ;
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

	if ( i1 == ele1 && i2 == ele2 && i3 == ele3 ) break ;

	int i12 = pair_index_ele(i1, i2) ;
    	int i23 = pair_index_ele(i2, i3) ;
	int i13 = pair_index_ele(i1, i3) ;
#ifdef TESTING
	printf("Skipping over element triple %d %d %d\n", i1, i2, i3) ;
	printf(" Pair index : %d %d %d\n", i12, i13, i23) ;
#endif
	for ( int i = 0 ; i < snum_3b_cheby[i12] ; i++ ) {
	  for ( int j = 0 ; j < snum_3b_cheby[i13] ; j++ ) {
	    for ( int k = 0 ; k < snum_3b_cheby[i23] ; k++ ) {
	      index += check_interaction_ordering(i, j, k, i12, i13, i23) ;
#ifdef TESTING
	      printf("Skip offset = %d %d %d %d\n", i, j, k, index-tot_2b) ;
#endif
	    }
	  }
	}
      }
    }
  }

  ipair12 = pair_index_ele(ele1, ele2) ;
  ipair23 = pair_index_ele(ele2, ele3) ;
  ipair13 = pair_index_ele(ele1, ele3) ;

  int n12_new = match_pair_order(a1, a2, a_orig, n12, n13, n23) ;
  int n13_new = match_pair_order(a1, a3, a_orig, n12, n13, n23) ;
  int n23_new = match_pair_order(a2, a3, a_orig, n12, n13, n23) ;

  // If the pair index is the same, the polynomial order must be sorted.  This keeps
  // the polynomial invariant to ordering of atoms.
  sort_poly_order(ipair12, ipair13, ipair23, n12_new, n13_new, n23_new) ;

  for ( int i = 0 ; i < snum_3b_cheby[ipair12] ; i++ ) {
    for ( int j = 0 ; j < snum_3b_cheby[ipair13] ; j++ ) {
      for ( int k = 0 ; k < snum_3b_cheby[ipair23] ; k++ ) {
	if ( i == n12_new && j == n13_new && k == n23_new ) goto DONE ;
	index += check_interaction_ordering(i, j, k, ipair12, ipair13, ipair23) ;
#ifdef TESTING
	printf("Pow offset = %d %d %d %d\n", i, j, k, index-tot_2b) ;
#endif
      }
    }
  }
 DONE: ;



  assert(index >= tot_2b && index < tot_short_range ) ;
#ifdef TESTING
  printf("Sorted:\n  ele1 = %d ele2 = %d ele3 = %d\n   n12 = %d n13 = %d n23 = %d\n",
	 ele1, ele2, ele3, n12_new, n13_new, n23_new) ;
  printf("  a1 = %d a2 = %d a3 = %d\n", a1, a2, a3) ;
  printf("  ipair12 = %d ipair13 = %d ipair23 = %d\n\n", ipair12, ipair13, ipair23) ;
  printf("  Index = %d 3-B Index = %d\n", index, index - tot_2b) ;
#endif



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
    for ( int i2 = i1 ; i2 < NELE ; i2++ ) {
      int i12 = pair_index_ele(i1, i2) ;
      for ( int i3 = i2 ; i3 < NELE ; i3++ ) {
	int i23 = pair_index_ele(i2, i3) ;
	int i13 = pair_index_ele(i1, i3) ;
	for ( int i = 0 ; i < snum[i12] ; i++ ) {
	  for ( int j = 0 ; j < snum[i13] ; j++ ) {
	    for ( int k = 0 ; k < snum[i23] ; k++ ) {
	      num_cheby_3b += check_interaction_ordering(i, j, k, i12, i13, i23) ;
	    }
	  }
	}
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

static int check_interaction_ordering(int i, int j, int k, int i12, int i13, int i23)
// Check to see if the polynomial powers correspond to standard order.  Return 1 if they
// do, 0 otherwise. i, j, k, are the polynomial powers corresponding to the 12, 13, and
// 23 interactions, respectively.  i12, i13, and i23 are the pair indices for the 
// interactions.
{
  if ( ! is_three_body(i,j,k) ) {
    return(0) ;
  }
  if ( i12 == i13 && i12 == i23 ) {
    // Order all coefficients in ascending order if all elements are the same.
    if ( i <= j && j <= k ) {
      return(1) ;
    } else {
      return(0) ;
    }
  } else if ( i12 == i13 ) {
    // When 2 pair coefficients are identical, the corresponding polynomial 
    // powers must be in ascending order.
    if ( i <= j ) {
      return(1) ;
    } else {
      return(0) ;
    }
  } else if ( i12 == i23 ) {
    if ( i <= k ) {
      return(1) ; ;
    } else {
      return(0) ;
    }
  } else if ( i13 == i23 ) {
    if ( j <= k ) {
      return(1) ; ;
    } else {
      return(0) ;
    }
  } else {
    // All elements are different.  No preferred order of coefficients.
    return(1) ; ;
  }
  // Not reached.
  return(0) ;
}

static bool is_three_body(int i, int j, int k)
// Return TRUE if the polynomial powers i,j, and k correspond to a 3-body interaction.
{
  int count = 0 ;

  // Count up how many polynomials of order 0 there are.
  if ( i == 0 ) {
    count++ ;
  }
  if ( j == 0 ) {
    count++ ;
  }
  if ( k == 0 ) {
    count++ ;
  }
  if ( count > 1 ) {
    return false ;
  }
  return true ;
}
