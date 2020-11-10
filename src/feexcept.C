// Floating point exception handling.  This is in it's own file to isolate the _GNU_SOURCE.

#ifdef ENABLE_FP_EXCEPT

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <fenv.h>

void enable_fp_exceptions() 
// Trap floating point exceptions.
{
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW /*| FE_UNDERFLOW*/ ) ;
}

#endif // defined(ENABLE_FP_EXCEPT)

