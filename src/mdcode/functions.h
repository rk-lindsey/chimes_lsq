#ifndef _HELPERS_
#define _HELPERS_


#include<stdio.h>
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <assert.h>
#include "functions.h"
using namespace std;

#define MAXOVERP 5  // Maximum number of overcoordination parameters allowed.
#define EMPTY 1.234e+56  // Value to signify an empty array entry.


const int NELE = 2 ;   // Number of elements.
const int NPAIR = (NELE+1) * NELE / 2 ;  // Number of pair interactions.

static const double ke=332.0637157615209;//converter between electron units and Stillinger units for Charge*Charge.

static const double Hartree = 627.50961 ; // 1 Hartree in kcal/mol.
static const double Kb  = 0.001987 ; // Boltzmann constant in kcal/mol-K.
static const double Tfs = 48.888 ;   // Internal time unit in fs.
static const double GPa = 6.9479 ;

enum Sr_pair_t {
  CHEBYSHEV,
  DFTBPOLY,
  SPLINE,
  INVERSE_R,
  LJ,
  STILLINGER
} ;


void ZCalc(double **Coord, const char *Lbc, double *Q, double *Latcons,
	   const int nlayers,
	   const int nat,const double *smin,const double *smax,
	   const double *sdelta,const int *snum, 
	   const int *snum_3b_cheby,
	   double *params, double *pot_params, Sr_pair_t pair_type,
	   bool if_coulomb, bool if_overcoord, bool if_3b_cheby,
	   int n_over,
	   double *over_params, const double *lambda,
	   double **SForce,double& Vtot,double& Pxyz) ;

void ZCalc_Deriv(double **Coord,const char *Lbc,
		 double *Latcons,const int nlayers,
		 const int nat,double ***A,const double *smin,const double *smax,
		 const double *sdelta, const int *snum, const int *snum_3b_cheby,
		 const double *lambda,
		 double **coul_oo, double **coul_oh,double **coul_hh, Sr_pair_t pair_type,
		 double *mind, bool if_3b_cheby) ;

void SubtractCoordForces(double **Coord,double **Force,string *Lb, double *Latcons,
			 const int nat, bool calc_deriv, 
			 double **Fderiv, int n_over, double *over_param) ;

void ZCalc_Ewald(double **Coord, const char *Lbc, const double *Q, const double *Latcons,
		 const int nat, double **SForce,double& Vtot,double& Pxyz) ;

void ZCalc_Ewald_Orig(double **Coord,string *Lb, const double *Latcons,
		      const int nat,double **SForce,double& Vtot,double& Pxyz) ;

double bondedpot(double **Coord_bonded,double ***I_bonded);
double spline_pot(double smin, double smax, double sdelta, double rlen2, double *params, double *pot_params, int snum, int vstart, double &S_r) ;

bool parse_tf(char *val, int bufsz, char *line) ;

bool read_tf_option(ifstream *paramread, const char *option, const char *params_file) ;

Sr_pair_t parse_pair_type(const char *val, int bufsz)  ;

void optimal_ewald_params(double accuracy, 
			  double V, int nat, double &alpha, 
			  double &rc, int &kc, double &r_acc,
			  double &k_acc) ;


void ZCalc_Ewald_Deriv(double **Coord, const char *Lbc, 
		       double *Latcons,
		       const int nat,
		       double **coul_oo,double **coul_oh,double **coul_hh) ;

void parse_param_list(double *params, int nparams, const char* name) ;

int pair_index(int a1, int a2, const char *Lbc) ; 
int atom_index(int a1, const char *Lbc) ;

void ZCalc_3B_Cheby(double **Coord,const char *Lbc, double *Latcons,
		    const int nat,const double *smin,
		    const double *smax,
  		    const int *snum_3b_cheby,
		    double ******idx_params, const double *lambda,
		    double **SForce, double &Vtot, double &Pxyz) ;

void ZCalc_3B_Cheby_Deriv(double **Coord,const char *Lbc, double *Latcons,
			  const int nat, double ***A,
			  const double *smin,
			  const double *smax,
			  const int *snum, 
			  const int *snum_3b_cheby,
			  const double *lambda, int ******index_params) ;
double ******Indexed_3B_Cheby_Coeffs(const char *Lbc, 
				  const int nat,
				  const int *snum, 
				  const int *snum_3b_cheby,
				     double *params) ;
int ******Index_3B_Cheby(const char *Lbc, 
			const int nat,
			const int *snum, 
			const int *snum_3b_cheby
			) ;

int count_cheby_3b_params(const int *snum) ;

int pair_index_ele(int ele1, int ele2) ;

#endif

