/******************************************************************************/
/*

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 20080

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.
*/	 
		
#ifndef ZIGGURAT_H
#define ZIGGURAT_H

float r4_exp ( unsigned long int *jsr, int ke[256], float fe[256], 
  float we[256] );
void r4_exp_setup ( int ke[256], float fe[256], float we[256] );
float r4_nor ( unsigned long int *jsr, int kn[128], float fn[128], 
  float wn[128] );
void r4_nor_setup ( int kn[128], float fn[128], float wn[128] );
float r4_uni ( unsigned long int *jsr );
unsigned long int shr3 ( unsigned long int *jsr );
void timestamp ( void );

#endif
