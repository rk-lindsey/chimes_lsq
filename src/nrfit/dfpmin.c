#include <stdio.h>
#include <math.h>

#include "dnrutil.h"
#include "dnr.h"

static void   (*dfunc)(double*,double*) ;
static double (*ffunc)(double*) ;
static void dfunc2(double*,double*) ;
static int Np ;

#define ITMAX 200
#define EPS 1.0e-08

void dfpmin(p,n,ftol,iter,fret,func)
double p[],ftol,*fret,(*func)(double*);
int n,*iter;
{
	int j,i,its;
	double fp,fae,fad,fac;
	double *xi,*g,*dg,*hdg,*vector();
	double **hessin,**matrix();
	void linmin(),nrerror(),free_matrix(),free_vector();

	ffunc = func ;
	dfunc = &dfunc2 ;
	Np = n ;

	hessin=matrix(1,n,1,n);
	xi=vector(1,n);
	g=vector(1,n);
	dg=vector(1,n);
	hdg=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,g);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin(p,xi,n,fret,func);

		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			free_vector(hdg,1,n);
			free_vector(dg,1,n);
			free_vector(g,1,n);
			free_vector(xi,1,n);
			free_matrix(hessin,1,n,1,n);
			return;
		}
		fp=(*fret);
		for (i=1;i<=n;i++) dg[i]=g[i];
		*fret=(*func)(p);
		(*dfunc)(p,g);
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
		}
		fac=1.0/fac;
		fad=1.0/fae;
		for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
		for (i=1;i<=n;i++)
			for (j=1;j<=n;j++)
				hessin[i][j] += fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
	nrerror("Too many iterations in DFPMIN");
}

#undef ITMAX
#undef EPS


#define EPS 1.0e-4

static void dfunc2(double *p, double*g )
{
    double z0, z1, temp, h ;
    int i ;

    z0 = (*ffunc)(p) ;
    
    for ( i = 1 ; i <= Np ; i++ ) {
	temp=p[i];
	h=EPS*fabs(temp);
	if (h == 0.0) h=EPS;
	p[i]=temp+h;
	h=p[i]-temp;
	z1 = (*ffunc)(p);
	g[i] = (z1-z0) / h ;
	p[i] = temp ;
    }
}

#undef EPS
#undef NRANSI

/* (C) Copr. 1986-92 Numerical Recipes Software Y229&kH. */

 
 
