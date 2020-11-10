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

#include "xvector3.h"
#include "tools.h"



void CxVector3::Dump() const
{
	BXIN;
/*	int z;
	for (z=0;z<3;z++)
		if ((m_pData[z] < 0) && (m_pData[z] > -0.000000001))
			m_pData[z] = 0;*/
	mprintf("( %6.3f | %6.3f | %6.3f ) {%.2f}",m_pData[0],m_pData[1],m_pData[2],GetLength());
	BXOUT;
}


void CxVector3::PointRoot(const CxVector3 &vec1, const CxVector3 &vec2, const CxVector3 &point)
{
	BXIN;
	float a, b;
	CxVector3 vn;

	vn = CrossP(vec1,vec2);

	#define ax vec1[0]
	#define ay vec1[1]
	#define az vec1[2]
	#define bx vec2[0]
	#define by vec2[1]
	#define bz vec2[2]
	#define nx vn[0]
	#define ny vn[1]
	#define nz vn[2]
	#define px point[0]
	#define py point[1]
	#define pz point[2]

	a = -((bz*ny*px - by*nz*px - bz*nx*py + bx*nz*py + by*nx*pz - bx*ny*pz)/(-az*by*nx + ay*bz*nx + az*bx*ny - ax*bz*ny - ay*bx*nz + ax*by*nz));
	b = -((az*ny*px - ay*nz*px - az*nx*py + ax*nz*py + ay*nx*pz - ax*ny*pz)/( az*by*nx - ay*bz*nx - az*bx*ny + ax*bz*ny + ay*bx*nz - ax*by*nz));

/*	vec1 *= a;
	vec2 *= b;
	vecadd(vec1,vec2,baseout);*/
	*this = vec1*a + vec2*b;
	BXOUT;
}


float Dihedral(const CxVector3 &vec1, const CxVector3 &vec2, const CxVector3 &norm, bool absolute)
{
	BXIN;
	CxVector3 p1, p2;
	CxVector3 t1, t2;
	float f;

	p1 = CrossP(norm,vec1);
	p2 = CrossP(norm,p1);
	if ((p1.GetLength() == 0) || (p2.GetLength() == 0))
	{
		eprintf("Error in Dihedral.\n");
		return -1.0f;
	}
	t1.PointRoot(p1,p2,vec1);
	t2.PointRoot(p1,p2,vec2);
	f = Angle_Deg(t1,t2);
	if (!absolute)
	{
		if (fabs(Angle_Deg(p1,t2)) > 90.0)
			f = -f;
	}
	BXOUT;
	return f;
}


float Dihedral2(const CxVector3 &vec1, const CxVector3 &vec2, const CxVector3 &norm, bool absolute)
{
	BXIN;
	CxVector3 p1, p2, n;
//	CxVector3 t1, t2;
	float f;

	(void)absolute;

	n = norm;
	n.Normalize();

	p1 = vec1 - n*DotP(vec1,n);
	p2 = vec2 - n*DotP(vec2,n);
	if ((p1.GetLength() == 0) || (p2.GetLength() == 0))
	{
		eprintf("Error in Dihedral.\n");
		return -1.0f;
	}
	f = Angle_Deg(p1,p2);
/*	if (!absolute)
	{
		if (fabs(Angle_Deg(p1,t2)) > 90.0)
			f = -f;
	}*/
	BXOUT;
	return f;
}


CxVector3 PointFromRAD(CxVector3 r1, CxVector3 r2, CxVector3 r3, float r, float a, float d)
{
	CxVector3 res, d1, d2, d3;

	d1 = r2 - r1;
//	mprintf("d1 = "); d1.Dump(); mprintf("\n");

	d2 = r3 - r2;
//	mprintf("d2 = "); d2.Dump(); mprintf("\n");

	d3 = CrossP(d1,d2);
//	mprintf("d3 = "); d3.Dump(); mprintf("\n");

	d1.Normalize();
//	mprintf("d1n = "); d1.Dump(); mprintf("\n");

	d2 = d2 - DotP(d2,d1)*d1;
	d2.Normalize();
//	mprintf("d2n = "); d2.Dump(); mprintf("\n");

	d3.Normalize();
//	mprintf("d3n = "); d3.Dump(); mprintf("\n");

	res = r1 + r * ( (float)cos(a) * d1 + (float)(sin(a) * cos(d)) * d2 + (float)(sin(a) * sin(d)) * d3);
//	mprintf("res = "); res.Dump(); mprintf("\n");

	return res;
}

