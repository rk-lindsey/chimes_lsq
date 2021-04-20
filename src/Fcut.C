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
#include "util.h"

void FCUT::get_fcut(double & fcut, double & fcut_deriv, const double rlen, const double rmin, const double rmax) 
// Calculates the fcut function and its derivative. Used to force the potential to zero at rcut max, 
// and for 3-body interactions, to zero at rcut min
{	
	static double fcut0;
	static double fcut0_deriv;
	static double THRESH;
	static double A, a, a_prime;
	static double B, b, b_prime;
	
	static bool FIRST_PASS = true;

	// Original cubic style smoothing... does not constrian value at rmin
	// Cubic cutoff
	if(TYPE == FCUT_TYPE::CUBIC)
	{		
		fcut0 = (1.0 - rlen/rmax);
		fcut        = pow(fcut0, POWER);
		fcut_deriv  = pow(fcut0,POWER-1);
		fcut_deriv *= -1.0 * POWER /rmax;

		return;
	}
	else if(TYPE == FCUT_TYPE::TERSOFF)
	{
		// FYI: For this cutoff type, "OFFSET" is actually a fraction of the outer cutoff
		
		THRESH = rmax-OFFSET*rmax;		
		
		if      (rlen < THRESH)		// Case 1: Our pair distance is less than the fcut kick-in distance
		{
			fcut       = 1.0;
			fcut_deriv = 0.0;
		}
		else if (rlen > rmax)		// Case 2: Our pair distance is greater than the cutoff
		{
			fcut       = 0.0;
			fcut_deriv = 0.0;
		}					
		else				// Case 3: We'll use our modified sin function
		{
			fcut0       = (rlen-THRESH) / (rmax-THRESH) * pi + pi/2.0;
			fcut0_deriv = pi / (rmax - THRESH);
			
			fcut        = 0.5 + 0.5 * sin( fcut0 );
			fcut_deriv  = 0.5 * cos( fcut0 ) * fcut0_deriv; 
		}
		return;	
	}
	else
	{
		cout << "ERROR: UNDEFINED CUTOFF TYPE!"  << endl;
		cout << "       Check calls to FCUT::get_fcut." << endl;
		exit(0);
	}
}

FCUT::FCUT() 
{
	TYPE = FCUT_TYPE::CUBIC;

	// More or less random defaults.

	OFFSET = 0.0;

}



void FCUT::set_type(string s) 
// Set the type of the cutoff function.
{
	else if ( s == "CUBIC" ) 
	{
		TYPE = FCUT_TYPE::CUBIC;
	} 
	else if ( s == "TERSOFF" ) 
	{
	
		TYPE = FCUT_TYPE::TERSOFF;
	} 
	else 
	{
		cout << "Error: unknown cutoff type " << s << endl;
		exit(1);
	}
}

string FCUT::to_string()
// Convert an fcut_type to a string.
{
	string str;

	switch ( TYPE ) {
	case FCUT_TYPE::CUBIC:
		str = "CUBIC";
		break;
	case FCUT_TYPE::TERSOFF:
		str = "TERSOFF";
		break;		
	default:
		cout << "Error: unknown cutoff value " << (int) TYPE << endl;
		exit(1);
	}
	return(str);
}


bool FCUT::PROCEED(const double & rlen, const double & rmin, const double & rmax)
// Determines whether rlen is within permissible ranges for the given fcut type... i.e. cubic doesn't care about inner cutoff, while other types do.
{

	if(TYPE == FCUT_TYPE::CUBIC 
		|| TYPE == FCUT_TYPE::TERSOFF)
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


void FCUT::parse_input(string line)
// Parse the input string and set parameters for the force cutoff
{
	vector<string> tokens;

	int nargs = parse_space(line, tokens);

	validate_num_args(nargs, 1, line);

	set_type(tokens[0]);
	
	if (TYPE == FCUT_TYPE::TERSOFF)
	{
		validate_num_args(nargs, 2, line);
		OFFSET = stod(tokens[1]);
	}
}

void FCUT::print_params()
// Print the force cutoff parameters
{
  cout << to_string() << endl;
  cout << "		...with offset of: " << fixed << setprecision(4) << OFFSET << endl;
}


void FCUT::print_header(ostream &header)
// Print force cutoff function parameters to the header file.
{
  header << endl << "FCUT TYPE: " << to_string();
		
  if(TYPE == FCUT_TYPE::TERSOFF)
	 header << " " << OFFSET;	 

  header << endl;
}
