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

	double STEEPNESS;
	double OFFSET; 
	double HEIGHT;

	// power used in the cutoff function.
	int  POWER ;

	// 2, 3, 4, body interaction.
	int BODIEDNESS ;

	// Type of cutoff-function employed.
	FCUT_TYPE TYPE ;

	// set the type of the cutoff function.
	void set_type(string s) ;

	// Default constructor.
	FCUT() ;

	// Convert a cut-off function to a string.
	string to_string() ;

	// Evaluate the cut-off function.
	void get_fcut(double & fcut, double & fcut_deriv, const double rlen, 
								const double rmin, const double rmax) ;

	// Decide whether to proceed with a pair interaction.
	bool PROCEED(const double & rlen, const double & rmin, const double & rmax) ;

} ;

#endif
