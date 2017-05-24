// File describing force cut-off functions.
//
#ifndef _FCUT_H
#define _FCUT_H

enum class FCUT_TYPE {
	// Contains allowed force cut-off types.
	CUBESTRETCH,
	CUBIC,
	COSINE,
	CUBSIG,
	SIGMOID,
	SIGFLT
} ;

class FCUT {
public:
	// Convert a string to an FCUT_TYPE.
	static FCUT_TYPE to_val(string s) ;

	// Convert a cut-off function to a string.
	static string to_string(FCUT_TYPE val) ;

	// Evaluate the cut-off function.
	static void get_fcut(int BODIEDNESS, FCUT_TYPE TYPE,
								double & fcut, double & fcut_deriv, const double rlen, 
								const double rmin, const double rmax, const int fcut_power, 
								double OFFSET, double STEEPNESS, double HEIGHT) ;
} ;

#endif
