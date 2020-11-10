/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmmin.h
 *
 * Contents: Public interface to the Levenberg-Marquardt core implementation.
 *
 * Author:   Joachim Wuttke 2004-2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 */
 
#ifndef LMMIN_H
#define LMMIN_H

//#ifdef __cplusplus
//extern "C" {
//#endif

#include <float.h>

/** Compact high-level interface. **/

/* Collection of control (input) parameters. */
typedef struct {
    double ftol;      /* relative error desired in the sum of squares. */
    double xtol;      /* relative error between last two approximations. */
    double gtol;      /* orthogonality desired between fvec and its derivs. */
    double epsilon;   /* step used to calculate the jacobian. */
    double stepbound; /* initial bound to steps in the outer loop. */
    int maxcall;      /* maximum number of iterations. */
    int scale_diag;   /* UNDOCUMENTED, TESTWISE automatical diag rescaling? */
    int printflags;   /* OR'ed to produce more noise */
} lm_control_struct;

/* Collection of status (output) parameters. */
typedef struct {
    double fnorm;     /* norm of the residue vector fvec. */
    int nfev;	      /* actual number of iterations. */
    int info;	      /* status of minimization. */
} lm_status_struct;

/* Recommended control parameter settings. */
extern const lm_control_struct lm_control_double;
extern const lm_control_struct lm_control_float;

/* Standard monitoring routine. */
void lm_printout_std( int n_par, const double *par, int m_dat,
                      const void *data, const double *fvec,
                      int printflags, int iflag, int iter, int nfev);

/* Refined calculation of Eucledian norm, typically used in printout routine. */
double lm_enorm( int, const double * );

/* The actual minimization. */
void lmmin( int n_par, double *par, int m_dat, const void *data, 
            void (*evaluate) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info),
            const lm_control_struct *control, lm_status_struct *status,
            void (*printout) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              int printflags, int iflag, int iter, int nfev) );


/** Legacy low-level interface. **/

/* Alternative to lm_minimize, allowing full control, and read-out
   of auxiliary arrays. For usage, see implementation of lm_minimize. */
void lm_lmdif( int m, int n, double *x, double *fvec, double ftol,
	       double xtol, double gtol, int maxfev, double epsfcn,
	       double *diag, int mode, double factor, int *info, int *nfev,
	       double *fjac, int *ipvt, double *qtf, double *wa1,
	       double *wa2, double *wa3, double *wa4,
               void (*evaluate) (const double *par, int m_dat, const void *data,
                                 double *fvec, int *info),
               void (*printout) (int n_par, const double *par, int m_dat,
                                 const void *data, const double *fvec,
                                 int printflags, int iflag, int iter, int nfev),
               int printflags, const void *data );

extern const char *lm_infmsg[];
extern const char *lm_shortmsg[];


/*****************************************************************************/
/*  set numeric constants                                                    */
/*****************************************************************************/

/* machine-dependent constants from float.h */
#define LM_MACHEP     DBL_EPSILON   /* resolution of arithmetic */
#define LM_DWARF      DBL_MIN       /* smallest nonzero number */
#define LM_SQRT_DWARF sqrt(DBL_MIN) /* square should not underflow */
#define LM_SQRT_GIANT sqrt(DBL_MAX) /* square should not overflow */
#define LM_USERTOL    30*LM_MACHEP  /* users are recommended to require this */

/* If the above values do not work, the following seem good for an x86:
 LM_MACHEP     .555e-16
 LM_DWARF      9.9e-324	
 LM_SQRT_DWARF 1.e-160   
 LM_SQRT_GIANT 1.e150 
 LM_USERTOL   1.e-14*/

 /*  The following values should work on any machine:
#define LM_MACHEP     1.2e-16
#define LM_DWARF      1.0e-38
#define LM_SQRT_DWARF 3.834e-20
#define LM_SQRT_GIANT 1.304e19
#define LM_USERTOL    1.e-14*/


/*const lm_control_struct lm_control_double = {
    LM_USERTOL, LM_USERTOL, LM_USERTOL, LM_USERTOL, 1.0, 2000, 1, 0 };

const lm_control_struct lm_control_float = {
    1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.0, 100, 0, 0 };*/



//#ifdef __cplusplus
//}
//#endif

#endif /* LMMIN_H */
