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


#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include "qr_fact.h"
//#include "defs_and_types.h"


//#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
//#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SIGNF(a, b) ((b) >= 0.0 ? fabsf(a) : -fabsf(a))
//#define SQR(x) ((x)*(x))


/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

 
static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


int ComputeSVD(float **a, int m, int n, float *w, float **v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        eprintf("ComputeSVD(): #rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = (float)((double)a[k][i]/scale);
                    s += ((double)a[k][i] * (double)a[k][i]);
                }
                f = (double)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (float)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k][i] * (double)a[k][j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += (float)(f * (double)a[k][i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = (float)((double)a[k][i]*scale);
            }
        }
        w[i] = (float)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = (float)((double)a[i][k]/scale);
                    s += ((double)a[i][k] * (double)a[i][k]);
                }
                f = (double)a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = (float)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j][k] * (double)a[i][k]);
                        for (k = l; k < n; k++) 
                            a[j][k] += (float)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = (float)((double)a[i][k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j][i] = (float)(((double)a[i][j] / (double)a[i][l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i][k] * (double)v[k][j]);
                    for (k = l; k < n; k++) 
                        v[k][j] += (float)(s * (double)v[k][i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k][i] * (double)a[k][j]);
                    f = (s / (double)a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += (float)(f * (double)a[k][i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = (float)((double)a[j][i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j][nm];
                            z = (double)a[j][i];
                            a[j][nm] = (float)(y * c + z * s);
                            a[j][i] = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (float)(-z);
                    for (j = 0; j < n; j++) 
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                eprintf("ComputeSVD(): No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj][j];
                    z = (double)v[jj][i];
                    v[jj][j] = (float)(x * c + z * s);
                    v[jj][i] = (float)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj][j];
                    z = (double)a[jj][i];
                    a[jj][j] = (float)(y * c + z * s);
                    a[jj][i] = (float)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (float)x;
        }
    }
    free((void*) rv1);
    return(1);
}


/* m - rows, n - cols
   a[zm][zn] = a[zm*n + zn]; */
int ComputeSVD_Flat(float *a, int m, int n, float *w, float *v)
{
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        eprintf("ComputeSVD(): #rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs((double)a[k*n+i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k*n+i] = (float)((double)a[k*n+i]/scale);
                    s += ((double)a[k*n+i] * (double)a[k*n+i]);
                }
                f = (double)a[i*n+i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i*n+i] = (float)(f - g);
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += ((double)a[k*n+i] * (double)a[k*n+j]);
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k*n+j] += (float)(f * (double)a[k*n+i]);
                    }
                }
                for (k = i; k < m; k++) 
                    a[k*n+i] = (float)((double)a[k*n+i]*scale);
            }
        }
        w[i] = (float)(scale * g);
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs((double)a[i*n+k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i*n+k] = (float)((double)a[i*n+k]/scale);
                    s += ((double)a[i*n+k] * (double)a[i*n+k]);
                }
                f = (double)a[i*n+l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i*n+l] = (float)(f - g);
                for (k = l; k < n; k++) 
                    rv1[k] = (double)a[i*n+k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += ((double)a[j*n+k] * (double)a[i*n+k]);
                        for (k = l; k < n; k++) 
                            a[j*n+k] += (float)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++) 
                    a[i*n+k] = (float)((double)a[i*n+k]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((double)w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j*n+i] = (float)(((double)a[i*n+j] / (double)a[i*n+l]) / g);
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += ((double)a[i*n+k] * (double)v[k*n+j]);
                    for (k = l; k < n; k++) 
                        v[k*n+j] += (float)(s * (double)v[k*n+i]);
                }
            }
            for (j = l; j < n; j++) 
                v[i*n+j] = v[j*n+i] = 0.0;
        }
        v[i*n+i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = (double)w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i*n+j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += ((double)a[k*n+i] * (double)a[k*n+j]);
                    f = (s / (double)a[i*n+i]) * g;
                    for (k = i; k < m; k++) 
                        a[k*n+j] += (float)(f * (double)a[k*n+i]);
                }
            }
            for (j = i; j < m; j++) 
                a[j*n+i] = (float)((double)a[j*n+i]*g);
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j*n+i] = 0.0;
        }
        ++a[i*n+i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs((double)w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = (double)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (float)h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = (double)a[j*n+nm];
                            z = (double)a[j*n+i];
                            a[j*n+nm] = (float)(y * c + z * s);
                            a[j*n+i] = (float)(z * c - y * s);
                        }
                    }
                }
            }
            z = (double)w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = (float)(-z);
                    for (j = 0; j < n; j++) 
                        v[j*n+k] = (-v[j*n+k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                eprintf("ComputeSVD(): No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = (double)w[l];
            nm = k - 1;
            y = (double)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = (double)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = (double)v[jj*n+j];
                    z = (double)v[jj*n+i];
                    v[jj*n+j] = (float)(x * c + z * s);
                    v[jj*n+i] = (float)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (float)z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = (double)a[jj*n+j];
                    z = (double)a[jj*n+i];
                    a[jj*n+j] = (float)(y * c + z * s);
                    a[jj*n+i] = (float)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (float)x;
        }
    }
    free((void*) rv1);
    return(1);
}


void ComputePseudoInverse(int m, float *m_in, float *m_out)
{
	float *mat_a, *mat_b, *mat_c, tf;
	int z1, z2, z3;
//	int z;

	mat_a = new float[m*m];
	mat_b = new float[m];
	mat_c = new float[m*m];

	memcpy(mat_a,m_in,m*m*sizeof(float));

	ComputeSVD_Flat(mat_a,m,m,mat_b,mat_c);

/*	mprintf("*** U Matrix:\n");
	mprintf("{ ");
	for (z=0;z<m;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<m;z2++)
		{
			mprintf("%10G",mat_a[z*m+z2]);
			if (z2 < m-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < m-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	mprintf("*** D Matrix:\n");
	mprintf("{ ");
	for (z=0;z<m;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<m;z2++)
		{
			if (z == z2)
				mprintf("%10G",mat_b[z]);
					else mprintf("%10G",0.0);
			if (z2 < m-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < m-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	mprintf("*** V Matrix:\n");
	mprintf("{ ");
	for (z=0;z<m;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<m;z2++)
		{
			mprintf("%10G",mat_c[z*m+z2]);
			if (z2 < m-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < m-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");*/

	for (z1=0;z1<m;z1++)
		if (mat_b[z1] != 0)
			mat_b[z1] = 1.0f / mat_b[z1];

	for (z1=0;z1<m;z1++)
	{
		for (z2=0;z2<m;z2++)
		{
			tf = 0;

			for (z3=0;z3<m;z3++)
				tf += mat_c[z1*m+z3] * mat_b[z3] * mat_a[z2*m+z3];

			m_out[z1*m+z2] = tf;
		}
	}

/*	mprintf("*** Pseudoinverse:\n");
	mprintf("{ ");
	for (z=0;z<m;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<m;z2++)
		{
			mprintf("%10G",m_out[z*m+z2]);
			if (z2 < m-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < m-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");*/

	delete[] mat_a;
	delete[] mat_b;
	delete[] mat_c;
}


void TestPseudoInverse()
{
	float *mat_in, *mat_out;
	int mat_M, z, z2;

	mat_M = 5;

	mat_in = new float[mat_M*mat_M];
	mat_out = new float[mat_M*mat_M];

	for (z=0;z<mat_M;z++)
	{
//		mat_a[z] = new float[mat_M];
//		mat_c[z] = new float[mat_M];
		for (z2=0;z2<mat_M;z2++)
			mat_in[z*mat_M+z2] = ((rand()%20001)-10000)/1000.0f;
	}

	mprintf("*** Input Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_in[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	ComputePseudoInverse(mat_M,mat_in,mat_out);

	mprintf("*** Moore-Penrose Pseudoinverse:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_out[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	delete[] mat_in;
	delete[] mat_out;
}


void TestSVD()
{
	// SVD-Test

	float *mat_a, *mat_b, *mat_c;
	int mat_M, z, z2;

	mat_M = 5;

	mat_a = new float[mat_M*mat_M];
	mat_b = new float[mat_M];
	mat_c = new float[mat_M*mat_M];

	for (z=0;z<mat_M;z++)
	{
//		mat_a[z] = new float[mat_M];
//		mat_c[z] = new float[mat_M];
		for (z2=0;z2<mat_M;z2++)
			mat_a[z*mat_M+z2] = ((rand()%20001)-10000)/1000.0f;
	}

	mprintf("*** Input Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_a[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	ComputeSVD_Flat(mat_a,mat_M,mat_M,mat_b,mat_c);

	mprintf("*** U Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_a[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	mprintf("*** D Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			if (z == z2)
				mprintf("%10G",mat_b[z]);
					else mprintf("%10G",0.0);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	mprintf("*** V Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_c[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");
}


void Solve_LeastSquares_QR(double *a, double *b, double *x, int m, int n)
{
	double **ta;
	int z;

	ta = new double*[m];
	for (z=0;z<m;z++)
	{
		ta[z] = new double[n];
		memcpy(ta[z],&a[z*n],n*sizeof(double));
	}

	QR_least_squares(ta, b, x, m, n);

	for (z=0;z<m;z++)
		delete[] ta[z];
	delete[] ta;
}


void TestLeastSquaresQR()
{
#define M 6
#define N 3

	double A[M*N], b[M], x[N];
	int z, z2;


	A[0*N+0] = 1;   A[0*N+1] = 2;   A[0*N+2] = 3;    b[0] = 3; 
	A[1*N+0] = 4;   A[1*N+1] = 5;   A[1*N+2] = 6;    b[1] = 9; 
	A[2*N+0] = 7;   A[2*N+1] = 8;   A[2*N+2] = 9;    b[2] = 15; 
	A[3*N+0] = 10;  A[3*N+1] = 11;  A[3*N+2] = 12;   b[3] = 22; 
	A[4*N+0] = 13;  A[4*N+1] = 14;  A[4*N+2] = 15;   b[4] = 27; 
	A[5*N+0] = 16;  A[5*N+1] = 17;  A[5*N+2] = -5;   b[5] = 33;

	mprintf("Testing Least-Squares Solver for system:\n\n");

	for (z=0;z<M;z++)
	{
		mprintf(" |");
		for (z2=0;z2<N;z2++)
		{
			mprintf(" %7.3f ",A[z2+z*N]);
		}
		mprintf("|");
		if (z == M/2)
			mprintf(" * ");
				else mprintf("   ");
		if ((z > (M-N)/2) && (z-(M-N)/2 <= N))
			mprintf("| x%d |",(z-(M-N)/2));
				else mprintf("      ");
		if (z == M/2)
			mprintf(" = ");
				else mprintf("   ");
		mprintf("| %7.3f |\n",b[z]);
	}
	mprintf("\n");

	mprintf("Solving...\n");
	Solve_LeastSquares_QR(A,b,x,M,N);

	mprintf("\nResult:\n\n");
	for (z=0;z<N;z++)
		mprintf("  x%d = %G\n",z+1,x[z]);
	mprintf("\n");
}
