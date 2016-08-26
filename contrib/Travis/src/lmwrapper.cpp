/*****************************************************************************
    TRAVIS - Trajectory Analyzer and Visualizer
    http://www.travis-analyzer.de/

    Copyright (c) 2009-2016 Martin Brehm
                  2012-2016 Martin Thomas

    This file written by Martin Brehm.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#include "lmwrapper.h"
#include <math.h>
#include "tools.h"
#include "globalvar.h"


double g_fZeroWeight;


CLMWrapper::CLMWrapper()
{
}


CLMWrapper::~CLMWrapper()
{
}


void lmcurve_evaluate(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	int z, z2;
	double t;
//	const double eps = 0.001;

	(void)info;  // Suppress "unused parameter" warning

	for (z2=0;z2<m_dat;z2++)
	{
		t = 0;
		for (z=0;z<g_iFitDegree;z++)
		{
			t += par[z*2] * exp(par[z*2+1]*((lmcurve_data_struct*)data)->x[z2]);
/*			if (par[z*2+1] > -eps)
				fvec[z2] += (par[z*2+1]+eps) * 1000.0;
			if (par[z*2] < 0)
				fvec[z2] += fabs(par[z*2]) * 1000.0;*/
		}
		fvec[z2] = ((lmcurve_data_struct*)data)->y[z2] - t;

		if (z2 == 0)
			fvec[z2] *= 100.0;

		for (z=0;z<g_iFitDegree;z++)
			if (par[z*2] < 0)
				fvec[z2] *= 1000.0;
	}
}


void lmcurve_evaluate_expspec(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	int z, z2;
	double t;

	(void)info;  // Suppress "unused parameter" warning

	memcpy(fvec,((lmcurve_data_struct*)data)->y,m_dat*sizeof(double));

//	mprintf("\n#### EVAL");
	for (z=0;z<g_iFitDegree;z++)
	{
//		if (z<10)
//			mprintf("\n    Par%3d: %G.",z,par[z]);
		t = exp(par[z]);
	//	t = par[z];

		for (z2=0;z2<m_dat;z2++)
			fvec[z2] -= t * exp(-g_pExpSpecExpo[z]*((lmcurve_data_struct*)data)->x[z2]);
	}

	fvec[0] *= g_fZeroWeight;

//	mprintf("\nEval: %G;  %G;  %G; ...",fvec[0],fvec[1],fvec[2]);

/*	for (z2=0;z2<m_dat;z2++)
	{
		t = 0;
		for (z=0;z<g_iFitDegree;z++)
			t += exp(par[z]) * exp(-g_pExpSpecExpo[z]*((lmcurve_data_struct*)data)->x[z2]);

		fvec[z2] = ((lmcurve_data_struct*)data)->y[z2] - t;

		if (z2 == 0)
			fvec[z2] *= g_fZeroWeight;
	}*/
}


void lmcurve_evaluate_sd(const double *par, int m_dat, const void *data, double *fvec, int *info)
{
	int z;
	(void)info;  // Suppress "unused parameter" warning

	for (z=0;z<m_dat;z++)
		fvec[z] = ((lmcurve_data_struct*)data)->y[z] - par[0] * (1.0 - exp(par[1] * pow( ((lmcurve_data_struct*)data)->x[z],par[2]) ));
}


double f_PolyExp(double t, const double *p)
{
	double r;
	int z;

	r = 0;
	for (z=0;z<g_iFitDegree;z++)
		r += p[z*2] * exp(p[z*2+1]*t);
	return r;
}


double f_ExpSpec(double t, const double *p)
{
	double r;
	int z;

	r = 0;
	for (z=0;z<g_iFitDegree;z++)
		r += exp(p[z]) * exp(-g_pExpSpecExpo[z]*t);
//		r += p[z] * exp(-g_pExpSpecExpo[z]*t);
	return r;
}


double f_SD(double t, const double *p)
{
	return p[0] * (1.0 - exp(p[1] * pow(t,p[2]) ));
}


void CLMWrapper::Fit_PolyExp(int degree, int n, double *x, double *y, double *par, double &r, double &integral, double *fitcurve, int maxcall)
{
	lm_status_struct status;
	lm_control_struct control;
	int z, z2, i;
	double ya, sse, sst, t;

	g_iFitDegree = degree;

	control.ftol = LM_USERTOL;
	control.xtol = LM_USERTOL;
	control.gtol = LM_USERTOL;
	control.epsilon = LM_USERTOL;
	control.stepbound = 1.0;
	control.maxcall = maxcall;
	control.scale_diag = 0;
	control.printflags = 0; // monitor status (+1) and parameters (+2)

	g_iLMMaxIter = maxcall;

	mprintf("    Performing Levenberg-Marquardt fit with %d-exponential curve...\n",g_iFitDegree);
	mprintf("      f(x) = c1 * exp(e1*x) + c2 * exp(e2*x) + ...\n");
	mprintf("      Starting values:\n");
 	for (z=0;z<g_iFitDegree;z++)
		mprintf("        c%-2d  = % 11.6f      e%-2d = %12.6f  (%12.5f ps)\n",z+1,par[z*2],z+1,par[z*2+1],-1.0/par[z*2+1]);

//	lmcurve_fit(g_iFitDegree*2, par, n, x, y, f_PolyExp, &control, &status );
    lmcurve_data_struct data = { x, y };

	g_bLMFitSilent = false;
	lmmin( g_iFitDegree*2, par, n, (const void*) &data,
           lmcurve_evaluate, &control, &status, lm_printout_std );

	for (z=0;z<g_iFitDegree;z++)
		if (par[z*2] < 0)
			par[z*2] = 0;

	if (fitcurve != NULL)
	{
		for (z=0;z<n;z++)
			fitcurve[z] = f_PolyExp(x[z],par);
	}

	ya = 0;
	for (z=0;z<n;z++)
		ya += y[z];
	ya /= n;

	sse = 0;
	sst = 0;

	for (z=0;z<n;z++)
	{
		sse += (y[z]-f_PolyExp(x[z],par)) * (y[z]-f_PolyExp(x[z],par));
		sst += (y[z]-ya) * (y[z]-ya);
	}

	mprintf("    sse = %.10G,  sst = %.10G,  sse/sst = %.10G\n",sse,sst,sse/sst);

	r = sqrt(1.0 - sse/sst);

//		par[z*2+1] = -(par[z*2+1] * par[z*2+1]);

	// Stack-Sort nach dem groessten Tau
	for (z=0;z<degree-1;z++)
	{
		t = -9e20;
		i = -1;
		for (z2=z;z2<degree;z2++)
		{
			if (-par[z2*2]/par[z2*2+1] > t)
			{
				t = -par[z2*2]/par[z2*2+1];
				i = z2;
			}
		}
		if ((i != -1) && (i != z))
		{
			t = par[z*2];
			par[z*2] = par[i*2];
			par[i*2] = t;
			t = par[z*2+1];
			par[z*2+1] = par[i*2+1];
			par[i*2+1] = t;
		}
	}

	mprintf("      Status: %s (%d cycles).\n",lm_infmsg[status.info],status.nfev);
	for (z=0;z<g_iFitDegree;z++)
	{
		mprintf("        c%-2d  = % 11.6f      e%-2d = % 12.6f",z+1,par[z*2],z+1,par[z*2+1]);
		if ((par[z*2+1] <= 0) && (par[z*2] >= 0))
			mprintf("    tau(%d) = %.8G ps",z+1,-par[z*2]/par[z*2+1]);
		mprintf("\n");
	}
	mprintf("        R    = % 13.9f",r);
	mprintf("  Chi^2 = % 12.8f\n",status.fnorm);
	integral = 0;
	for (z=0;z<g_iFitDegree;z++)
	{
		if (par[z*2+1] > 0)
		{
			mprintf("      Curve integral is infinite (not all exponents <= 0)!\n");
			integral = -1;
			goto _noint;
		}
		integral -= par[z*2]/par[z*2+1];
	}
	mprintf("      Curve integral: tau(total) = %.8G ps.\n",integral);
_noint:
	mprintf("    Fitting done.\n\n");
}


void CLMWrapper::Fit_ExpSpectrum(int degree, double mi, double ma, int n, double *x, double *y, double *par, int fitcurvepoints, double *fitcurvex, double *fitcurve, const char *name, int maxcall, double zeroweight, bool evolve)
{
	lm_status_struct status;
	lm_control_struct control;
	int z, z2, i;
	double ya, sse, sst, r;
//	char buf[256];
	CxString buf;
	FILE *a;

	mprintf("    Performing Lifetime Spectrum fit (max. %d iterations)...\n",maxcall);

	if (evolve)
	{
		mprintf("\n    Allocating %s of memory for evolution buffer...\n",FormatBytes((double)degree*maxcall*2*sizeof(double)));
		try { g_fLSpecEvolveBuf = new double[degree*maxcall*2]; } catch(...) { g_fLSpecEvolveBuf = NULL; }
		if (g_fLSpecEvolveBuf == NULL) NewException((double)degree*maxcall*2*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		for (z=0;z<maxcall*2;z++)
			g_fLSpecEvolveBuf[z*degree] = 0;
	}

	control.ftol = LM_USERTOL;
	control.xtol = LM_USERTOL;
	control.gtol = LM_USERTOL;
	control.epsilon = LM_USERTOL;
	control.stepbound = 1.0;
	control.maxcall = maxcall;
	control.scale_diag = 0;
	control.printflags = 0; // monitor status (+1) and parameters (+2)

//	mprintf("    Epsilon: %G\n",LM_USERTOL);

	g_iLMMaxIter = maxcall;

	g_iFitDegree = degree;

	lmcurve_data_struct data = { x, y };

/*	a = fopen("tmpfit.csv","wt");
	for (z=0;z<n;z++)
		fprintf(a,"%G;  %G\n",x[z],y[z]);
	fclose(a);*/

	g_pExpSpecExpo = new double[degree];
	g_fZeroWeight = zeroweight;

	mprintf("\n    Creating %d exponentials as basis for fitting. List of exponents (in ps^-1):",degree);

	for (z=0;z<degree;z++)
	{
		if ((z%6) == 0)
			mprintf("\n      ");
		g_pExpSpecExpo[z] = 1.0/exp(((double)z/(degree-1))*(log(ma)-log(mi))+log(mi));
		mprintf("%10.5f",g_pExpSpecExpo[z]);
		if (z < degree-1)
			mprintf(", ");
				else mprintf(".");
	}
	mprintf("\n\n");

	g_bLMFitSilent = false;
	g_bLMFitSmooth = false;

	lmmin( degree, par, n, (const void*) &data,
			 lmcurve_evaluate_expspec, &control, &status, lm_printout_std );

	mprintf("    Status: %s (%d cycles).\n",lm_infmsg[status.info],status.nfev);

	ya = 0;
	for (z=0;z<n;z++)
		ya += y[z];
	ya /= n;
	sse = 0;
	sst = 0;
	for (z=0;z<n;z++)
	{
		sse += (y[z]-f_ExpSpec(x[z],par)) * (y[z]-f_ExpSpec(x[z],par));
		sst += (y[z]-ya) * (y[z]-ya);
	}
	r = sqrt(1.0 - sse/sst);
	mprintf("      sse = %.10G,  sst = %.10G,  sse/sst = %.10G\n",sse,sst,sse/sst);
	mprintf("    R value of fit: %.10f\n\n",r);

//	sprintf(buf,"%s.csv",name);
	buf.sprintf("%s.csv",name);
	mprintf("    Saving lifetime spectrum to %s ...\n",(const char*)buf);
	a = OpenFileWrite(buf,true);

	if (evolve)
	{
		fprintf(a,"# Lifetime [ps];  Coefficient;  ln(Coefficient)");
		for (z=0;z<2*maxcall;z++)
		{
			if (g_fLSpecEvolveBuf[z*degree] == 0)
				break;
			fprintf(a,";  ln(C%d)",z+1);
		}
		i = z;
		fprintf(a,"\n");
		for (z=0;z<degree;z++)
		{
			fprintf(a,"%G;  %G;  %G",1.0/g_pExpSpecExpo[z],exp(par[z]),par[z]);
			for (z2=0;z2<i;z2++)
				fprintf(a,";  %G",g_fLSpecEvolveBuf[z2*degree+z]);
			fprintf(a,"\n");
		}
	} else
	{
		fprintf(a,"# Lifetime [ps];  Coefficient;  ln(Coefficient)\n");
		for (z=0;z<degree;z++)
			fprintf(a,"%G;  %G;  %G\n",1.0/g_pExpSpecExpo[z],exp(par[z]),par[z]);
	}

	fclose(a);

/*	if (smooth)
	{
		mprintf("\n    Refining fit with smoothness constraint...\n");

		g_bLMFitSmooth = true;

		lmmin( degree, par, n, (const void*) &data,
				 lmcurve_evaluate_expspec, &control, &status, lm_printout_std );

		ya = 0;
		for (z=0;z<n;z++)
			ya += y[z];
		ya /= n;
		sse = 0;
		sst = 0;
		for (z=0;z<n;z++)
		{
			sse += (y[z]-f_ExpSpec(x[z],par)) * (y[z]-f_ExpSpec(x[z],par));
			sst += (y[z]-ya) * (y[z]-ya);
		}
		r = sqrt(1.0 - sse/sst);
		mprintf("    R value of fit: %.10f\n\n",r);

		sprintf(buf,"%s_smooth.csv",name);
		mprintf("    Saving lifetime spectrum to %s ...\n",buf);
		a = OpenFileWrite(buf,true);
		fprintf(a,"# Exponent [ps];  Coefficient;  ln(Coefficient)\n");
		for (z=0;z<degree;z++)
			fprintf(a,"%G;  %g;  %g\n",g_pExpSpecExpo[z],exp(par[z]),par[z]);
		fclose(a);
	}*/

	g_bLMFitSmooth = false;

	if (fitcurve != NULL)
	{
		for (z=0;z<fitcurvepoints;z++)
			fitcurve[z] = f_ExpSpec(fitcurvex[z],par);
	}

	if (evolve)
		delete[] g_fLSpecEvolveBuf;

	mprintf("\n");
}


void CLMWrapper::Fit_SD(int n, double *x, double *y, double *par, double *fitcurve, double *r, int maxcall)
{
	lm_status_struct status;
	lm_control_struct control;
	int z;
	double ya, sse, sst;

	control.ftol = LM_USERTOL;
	control.xtol = LM_USERTOL;
	control.gtol = LM_USERTOL;
	control.epsilon = LM_USERTOL;
	control.stepbound = 1.0;
	control.maxcall = maxcall;
	control.scale_diag = 0;
	control.printflags = 0; // monitor status (+1) and parameters (+2)

	g_iLMMaxIter = maxcall;

	lmcurve_data_struct data = { x, y };

	g_bLMFitSilent = true;
	lmmin( 3, par, n, (const void*) &data,
			 lmcurve_evaluate_sd, &control, &status, lm_printout_std );
	g_bLMFitSilent = false;

	if (fitcurve != NULL)
	{
		for (z=0;z<n;z++)
			fitcurve[z] = f_SD(x[z],par);
	}

	ya = 0;
	for (z=0;z<n;z++)
		ya += y[z];
	ya /= n;

	sse = 0;
	sst = 0;

	for (z=0;z<n;z++)
	{
		sse += (y[z]-f_SD(x[z],par)) * (y[z]-f_SD(x[z],par));
		sst += (y[z]-ya) * (y[z]-ya);
	}


	*r = sqrt(1.0 - sse/sst);
}


