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


#include "xdvector3.h"
#include "xvector3.h"
#include "tools.h"


CxDVector3::CxDVector3(const CxVector3 &v)
{
	BXIN;
	m_pData[0] = v.GetAt(0);
	m_pData[1] = v.GetAt(1);
	m_pData[2] = v.GetAt(2);
	BXOUT;
}


void CxDVector3::Dump() const
{
	BXIN;
/*	int z;
	for (z=0;z<3;z++)
		if ((m_pData[z] < 0) && (m_pData[z] > -0.000000001))
			m_pData[z] = 0;*/
	mprintf("( %6.3f | %6.3f | %6.3f ) {%.2f}",m_pData[0],m_pData[1],m_pData[2],GetLength());
	BXOUT;
}


void CxDVector3::PointRoot(const CxDVector3 &vec1, const CxDVector3 &vec2, const CxDVector3 &point)
{
	BXIN;
	double a, b;
	CxDVector3 vn;

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


double Dihedral(const CxDVector3 &vec1, const CxDVector3 &vec2, const CxDVector3 &norm, bool absolute)
{
	BXIN;
	CxDVector3 p1, p2;
	CxDVector3 t1, t2;
	double f;

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


CxDVector3 PointFromRAD(CxDVector3 r1, CxDVector3 r2, CxDVector3 r3, double r, double a, double d)
{
	CxDVector3 res, d1, d2, d3;

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

	res = r1 + r * ( (double)cos(a) * d1 + (double)(sin(a) * cos(d)) * d2 + (double)(sin(a) * sin(d)) * d3);
//	mprintf("res = "); res.Dump(); mprintf("\n");

	return res;
}

