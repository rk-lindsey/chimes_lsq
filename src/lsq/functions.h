#ifndef _HELPERS_
#define _HELPERS_

#include<stdio.h>
#include<iostream>
#include<math.h>
#include<stdlib.h>
#include "functions.h"
using namespace std;


void ZCalc_Deriv(double **Coord,string *Lb, 
		 double *Latcons,const int nlayers,
		 const int nat,double ***A,const double smin,const double smax,
		 const double sdelta,const int snum, double **coul_oo,
		 double **coul_oh,double **coul_hh, bool if_cheby) ;

double bondedpot(double **Coord_bonded,double ***I_bonded);

void SubtractBondedForces(double **Coord,double **Force,string *Lb, double *Latcons,const int nlayers,
                          const int nat,double ***I_bonded);

void SubtractCoulombForces_standard(double **Coord,double **Force,string *Lb, double *Latcons,const int nlayers,
			   const int nat,double ***I_bonded);

void SubtractCoulombForces_Ewald(double **Coord,double **Force,string *Lb, double *Latcons,const int nlayers,
			   const int nat,double ***I_bonded);

void SubtractCoordForces(double **Coord,double **Force,string *Lb, 
			 double *Latcons,const int nlayers, const int nat, 
			 bool calc_deriv, double **Fderiv) ;

#endif

