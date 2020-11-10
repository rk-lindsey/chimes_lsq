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

#include "xdmatrix3.h"
#include "tools.h"


void CxDMatrix3::Dump() const
{
	BXIN;
	int z/*, z2*/;
	for (z=0;z<3;z++)
	{
/*		for (z2=0;z2<3;z2++)
			if ((GetAt(z2,z) < 0) && (GetAt(z2,z) > -0.000000001))
        			GetAt(z2,z) = 0;*/
		mprintf("( %6.3f %6.3f %6.3f )\n",GetAt(0,z),GetAt(1,z),GetAt(2,z));
	}
	BXOUT;
}


void CxDMatrix3::MatUltra(const CxDVector3 &vec1, const CxDVector3 &vec2)
{
	BXIN;
	CxDMatrix3 mat, mat2;
	CxDVector3 vectemp, vecx, vecy, vec2n, vectemp2;
/*	double mat[3][3], mat2[3][3];
	double vectemp[3], vecx[3], vecy[3], vec2n[3], vectemp2[3];*/
	double a;

#ifdef DEBUG_MATULTRA
	mprintf("\nLinAlg Checker\n\n");
	mprintf("Vektor 1: "); vec1.Dump();
	mprintf("\nVektor 2: "); vec2.Dump();
#endif

	vecx[0] = 1.0f;
	vecx[1] = 0.0f;
	vecx[2] = 0.0f;
	vecy[0] = 0.0f;
	vecy[1] = 1.0f;
	vecy[2] = 0.0f;

#ifdef DEBUG_MATULTRA
	mprintf("\nVektor X: "); vecx.Dump();
	mprintf("\nVektor Y: "); vecy.Dump();
#endif

	// Zu erst vec1 auf X drehen
	vectemp = CrossP(vec1,vecx);
	if (vectemp.GetLength() != 0)
	{
#ifdef DEBUG_MATULTRA
		mprintf("\nDrehachse fuer Vektor1->X: ");
		vectemp.Dump();
#endif
		vectemp.Normalize();
#ifdef DEBUG_MATULTRA
		mprintf("\nNormalisiert: ");
		vectemp.Dump();
#endif
		a = Angle(vec1,vecx);
#ifdef DEBUG_MATULTRA
		mprintf("\nDrehwinkel: %.1f Grad\n",a*180.0/Pi);
#endif
		mat.RotMat(vectemp,a);
	} else
	{
#ifdef DEBUG_MATULTRA
		mprintf("\nDer liegt schon auf der X-Achse!");
#endif
		mat.Unity();
	}
	vec2n = mat*vec2;

#ifdef DEBUG_MATULTRA
	mprintf("Drehmatrix dafuer:\n"); mat.Dump();
	vectemp = mat*vec1;
	mprintf("\nDamit wird Vektor1 zu "); vectemp.Dump();
	mprintf("\nUnd Vektor2 zu "); vec2n.Dump();
#endif

	if ((vec2-vec1).GetLength() < 0.0001)
	{
		*this = mat;
		BXOUT;
		return;
	}

	vectemp2[0] = 0;
	vectemp2[1] = vec2n[1];
	vectemp2[2] = vec2n[2];

#ifdef DEBUG_MATULTRA
	mprintf("\nProjektion von Vektor2n auf YZ-Ebene: ");
	vectemp2.Dump();
#endif

	a = Angle(vectemp2,vecy);
	
#ifdef DEBUG_MATULTRA
	mprintf("\n");
#endif

	if (vectemp2.GetLength()==0)
	{
		eprintf("\nMatUltra: Error caught.\n");
		a = (double)Pi/2.0f;
	}

	if (vec2n[2] > 0)
		a = -a; 
	
#ifdef DEBUG_MATULTRA
	mprintf("\nDrehachse fuer Vektor2->XY: "); 
	vecx.Dump();
	mprintf("\nDrehwinkel: %.1f Grad\n",a*180.0/Pi);
#endif

	mat2.RotMat(vecx,a);

#ifdef DEBUG_MATULTRA
	mprintf("Drehmatrix dafuer:\n"); mat2.Dump();
	vectemp = mat*vec1;
	vectemp2 = mat2*vectemp;
	mprintf("\nDamit wird Vektor1 zu "); vectemp2.Dump();
	vectemp = mat2*vec2n;
	mprintf("\nUnd Vektor2n zu "); vectemp.Dump();
#endif

	*this = mat2*mat;

#ifdef DEBUG_MATULTRA
	mprintf("\nMacht diese Gesamtmatrix:\n"); Dump();
	vectemp = *this*vec1;
	mprintf("\nDamit wird Vektor 1 von "); vec1.Dump(); mprintf(" zu "); vectemp.Dump();
	vectemp = *this*vec2;
	mprintf("\nUnd Vektor 2 wird von "); vec2.Dump(); mprintf(" zu "); vectemp.Dump();
	mprintf("\n");
#endif

	BXOUT;
}


void CxDMatrix3::Invert()
{
	#define a11 GetAt(0,0)
	#define a12 GetAt(0,1)
	#define a13 GetAt(0,2)
	#define a21 GetAt(1,0)
	#define a22 GetAt(1,1)
	#define a23 GetAt(1,2)
	#define a31 GetAt(2,0)
	#define a32 GetAt(2,1)
	#define a33 GetAt(2,2)

	CxDMatrix3 t;
	double det;

	det = a11 * (a33*a22 - a32*a23)
		 - a21 * (a33*a12 - a32*a13)
		 + a31 * (a23*a12 - a22*a13);

	if (det == 0)
	{
		eprintf("CxDMatrix3::Invert(): Failed - matrix is singular.\n");
		return;
	}

	t.GetAt(0,0) =  (a33*a22 - a32*a23);
	t.GetAt(0,1) = -(a33*a12 - a32*a13);
	t.GetAt(0,2) =  (a23*a12 - a22*a13);
	t.GetAt(1,0) = -(a33*a21 - a31*a23);
	t.GetAt(1,1) =  (a33*a11 - a31*a13);
	t.GetAt(1,2) = -(a23*a11 - a21*a13);
	t.GetAt(2,0) =  (a32*a21 - a31*a22);
	t.GetAt(2,1) = -(a32*a11 - a31*a12);
	t.GetAt(2,2) =  (a22*a11 - a21*a12);

	t *= 1.0/det;

	*this = t;
}

