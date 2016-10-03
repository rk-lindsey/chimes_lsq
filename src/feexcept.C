// Floating point exception handling.  This is in it's own file to isolate the _GNU_SOURCE.
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <fenv.h>

void enable_fp_exceptions() 
{
	feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW /*| FE_UNDERFLOW*/ ) ;
}

