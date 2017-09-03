// Definition of the FCUT class, controlling force cut-offs.
#include<iomanip>
#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#include<unistd.h>	// Used to detect whether i/o is going to terminal or is piped... will help us decide whether to use ANSI color codes
#include<string>
#include<limits>	// Help with handling of over/underflow

using namespace std;

#include "Fcut.h"
#include "functions.h"

void FCUT::get_fcut(double & fcut, double & fcut_deriv, 
						  const double rlen, const double rmin, const double rmax) 
// Calculates the fcut function and its derivative. Used to force the potential to zero at rcut max, 
// and for 3-body interactions, to zero at rcut min
{	
	static double fcut0;
	static double A, a, a_prime;
	static double B, b, b_prime;
	
	static bool FIRST_PASS = true;

	// Original cubic style smoothing... does not constrian value at rmin
	
	if(TYPE == FCUT_TYPE::CUBESTRETCH )
	{
		cout << "ERROR: TYPE CUBESTRETCH NO LONGER SUPPORTED...Please clean up code!" << endl;
		exit_run(0);
	}
	
	if(BODIEDNESS==2 || (BODIEDNESS==3 && TYPE == FCUT_TYPE::CUBIC))
	{		
		fcut0 = (1.0 - rlen/rmax);
		fcut        = pow(fcut0, POWER);
		fcut_deriv  = pow(fcut0,POWER-1);
		fcut_deriv *= -1.0 * POWER /rmax;

		return;
	}
	else if(BODIEDNESS==3)
	{
		if(TYPE == FCUT_TYPE::COSINE)
		{
			fcut0     = 2.0*pi*( (rlen - rmin) /(rmax - rmin) );
			fcut      = -0.5*cos(fcut0) + 0.5;
			fcut_deriv =  0.5*sin(fcut0) * (2.0*pi/(rmax - rmin));

			return;
		}
		else if(TYPE == FCUT_TYPE::CUBSIG)	// Outer cutoff handled by normal fcut function, inner cutoff is handled by a sigmoid, who is ~1 until r LESS THAN the inner cutoff. does NOT go to zero within rmin<->rmax
		{
			if(FIRST_PASS && RANK==0)
			{
				cout << "		NOTE: Using cubic-sigmoid steepness and offset of " << fixed << setprecision(3) << STEEPNESS << " and " << fixed << setprecision(3) << OFFSET << endl;
				FIRST_PASS = false;
			}
			
			A = exp( (rlen - rmin - OFFSET ) * STEEPNESS);	// A_prime is just STEEPNESS*A
			
			if (A>10E20) // If A is huge, then a = 1 and a' = 0
			{
				a       = 1;
				a_prime = 0;
			}
			else
			{
				a       = -1.0/(1.0+A) + 1;
				a_prime = STEEPNESS*A/(1.0+A)/(1.0+A);
			}
			
			fcut0 = (1.0 - rlen/rmax);
			b        = pow(fcut0, POWER);
			b_prime  = pow(fcut0,POWER-1);
			b_prime *= -1.0 * POWER /rmax;
			
			fcut = a*b;
			
			fcut_deriv = a_prime*b + b_prime*a;

			return;
		}
		
		else if( TYPE == FCUT_TYPE::SIGMOID) // Both outer and inner cutoffs handed by sigmoid functions
		{
			
			/* FOR DEBUGGING 
			
			if(FIRST_PASS && RANK==0)
			{
				cout << "		NOTE: Using sigmoid steepness and offset of " << fixed << setprecision(3) << STEEPNESS << " and " << fixed << setprecision(3) << OFFSET << endl;
				FIRST_PASS = false;
			}
			
			*/
			
			A = exp( (rlen - rmin - OFFSET ) * STEEPNESS);	// A_prime is just STEEPNESS*A
			B = exp( (rlen - rmax + OFFSET ) * STEEPNESS);	// Same for B_prime
			
			if (A>10E20) // If A is huge, then a = 1 and a' = 0
			{
				a       = 1;
				a_prime = 0;
			}
			else
			{
				a       = -1.0/(1.0+A) + 1;
				a_prime = STEEPNESS*A/(1.0+A)/(1.0+A);
			}
				
			if (B>10E20) // If B is huge, then b = 1 and b' = 0
			{
				b       = 1;
				b_prime = 0;
			}
			else
			{
				b       =  1.0/(1.0+B);
				b_prime = -1.0*STEEPNESS*B/(1.0+B)/(1.0+B);
			}

			fcut       = a*b;
			fcut_deriv = a_prime*b + b_prime*a;
			
			if(std::isinf(fcut)||std::isinf(fcut_deriv))
				cout << a << " " << b << " " << a_prime << " " << b_prime << " ++ " << A << " " << B << " -- " << fcut << " " << fcut_deriv << endl;

			return;
		}
		else if( TYPE== FCUT_TYPE::SIGFLT)
		{
			A = exp( (rlen - rmin - OFFSET ) * STEEPNESS);	// A_prime is just STEEPNESS*A
			B = exp( (rlen - rmax + OFFSET ) * STEEPNESS);	// Same for B_prime
			
			if (A>10E20) // If A is huge, then a = 1 and a' = 0
			{
				a       = 1;
				a_prime = 0;
			}
			else
			{
				a       =  HEIGHT/(1.0+A)+1;
				a_prime = -1.0*HEIGHT*STEEPNESS*A/(1.0+A)/(1.0+A);
			}
				
			if (B>10E20) // If B is huge, then b = 1 and b' = 0
			{
				b       = 1;
				b_prime = 0;
			}
			else
			{
				b       =  1.0/(1.0+B);
				b_prime = -1.0*STEEPNESS*B/(1.0+B)/(1.0+B);
			}

			fcut       = a*b;
			fcut_deriv = a_prime*b + b_prime*a;
			
			if(std::isinf(fcut)||std::isinf(fcut_deriv))
				cout << a << " " << b << " " << a_prime << " " << b_prime << " ++ " << A << " " << B << " -- " << fcut << " " << fcut_deriv << endl;

			return;
		}
		else
		{
			cout << "ERROR: UNDEFINED N-BODY PARAMETER!" << endl;
			cout << "       Did you forget to set # FCUTTYP #?" << endl;
			cout << "       Check calls to FCUT::get_fcut." << endl;
			exit(0);
		}
		
	}
	else
	{
		cout << "ERROR: UNDEFINED N-BODY PARAMETER!" << endl;
		cout << "       Check calls to FCUT::get_fcut." << endl;
		exit(0);
	}
}

FCUT::FCUT() 
{
#ifndef FPENALTY_POWER
	POWER = 3.0;
#else
	POWER = FPENALTY_POWER;
#endif

	BODIEDNESS = 2 ;
	
	TYPE = FCUT_TYPE::CUBIC ;

}



void FCUT::set_type(string s) 
// Set the type of the cutoff function.
{
	if ( s == "CUBESTRETCH" ) {
		TYPE = FCUT_TYPE::CUBESTRETCH ;
	} else if ( s == "CUBIC" ) {
		TYPE = FCUT_TYPE::CUBIC ;
	} else if ( s == "COSINE" ) {
		TYPE = FCUT_TYPE::COSINE ;
	} else if ( s == "CUBSIG" ) {
		TYPE = FCUT_TYPE::CUBSIG ;
	} else if ( s == "SIGMOID" ) {
		TYPE = FCUT_TYPE::SIGMOID ;
	} else if ( s == "SIGFLT" ) {
		TYPE = FCUT_TYPE::SIGFLT ;
	} else {
		cout << "Error: unknown cutoff type " << s << endl ;
		exit(1) ;
	}
}

string FCUT::to_string()
// Convert an fcut_type to a string.
{
	string str ;

	switch ( TYPE ) {
	case FCUT_TYPE::CUBESTRETCH:
		str = "CUBESTRETCH" ;
		break ;
	case FCUT_TYPE::CUBIC:
		str = "CUBIC" ;
		break ;
	case FCUT_TYPE::COSINE:
	   str = "COSINE" ;
		break ;
	case FCUT_TYPE::CUBSIG:
		str = "CUBSIG" ;
		break ;
	case FCUT_TYPE::SIGMOID:
		str = "SIGMOID" ;
		break ;
	case FCUT_TYPE::SIGFLT:
		str = "SIGFLT" ;
		break ;
	default:
		cout << "Error: unknown cutoff value " << (int) TYPE << endl ;
		exit(1) ;
	}
	return(str) ;
}


bool FCUT::PROCEED(const double & rlen, const double & rmin, const double & rmax)
// Determines whether rlen is within permissible ranges for the given fcut type... i.e. cubic doesn't care about inner cutoff, while other types do.
{

	if(TYPE == FCUT_TYPE::CUBIC 
		|| TYPE == FCUT_TYPE::CUBSIG
		|| TYPE == FCUT_TYPE::CUBESTRETCH 
		|| TYPE == FCUT_TYPE::SIGFLT)
		if(rlen < rmax)
			return true;
		else
			return false;
	else
		if((rlen < rmax) && (rlen > rmin))
			return true;
		else
			return false;
}
