// Program to generate a spline table for a Lennard-Jones interaction.
#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>
#include<stdlib.h>
#include<time.h>

int main(int argc, char **argv) 
{
  const double smin=0.0;
  const double smax=8.0;
  const double sdelta=0.05;
  const int snum=(1+int((smax-smin)/sdelta));
  double eps = 1.0 ;
  const double sigma = 0.75 ;
  FILE *params ;
  double r ;

  params = fopen("potparams.save", "w") ;

  // 3 copies for OO, OH, HH interaction.
  int l = 0 ;
  for ( int k = 0 ; k < 3 ; k++ ) {
    for ( int j = 0 ; j < snum  ; j++ ) {
      if ( j == 0 ) {
	// Set the force to a constant.
	r = sdelta ;
      } else {
	r = sdelta * j ;
      }
      
      // Potential
      double Vr = 4.0 * eps * ( pow(sigma/r,12.0) - pow(sigma/r,6.0) ) ;
      if ( k == 0 ) {
	fprintf(params,"%13.6e %13.6e\n", r, Vr) ;
      }


      // Force
      double dVdr = 4.0 * eps * ( 
				 -12.0 * pow(sigma/r,13.0) 
				 + 6.0 * pow(sigma/r, 7.0) ) ;
      dVdr /= sigma ;

      
      printf("%3d %21.13e\n", l, dVdr) ;
      l++ ;

      // Derivative of force.
      double d2Vdr2 = 4.0 * eps * ( 
				  12.0 * 13.0 * pow(sigma/r,14.0) 
				 - 6.0 *  7.0 * pow(sigma/r, 8.0) ) ;
      d2Vdr2 /= sigma * sigma ;

      printf("%3d %21.13e\n", l, d2Vdr2) ;
      l++ ;

    }
    // eps *= 2.0 ;
  }
}

    
