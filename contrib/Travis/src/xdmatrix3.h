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


#ifndef CXDMATRIX3_DEFINED
#define CXDMATRIX3_DEFINED

#include "tools.h"
#include "xobject.h"
#include "xdvector3.h"
#include "backtrace.h"


class CxDMatrix3 : public CxObject
{
public:
	CxDMatrix3() { }

	~CxDMatrix3() { }

/*	CxDMatrix3(CxDMatrix3 m)
	{
		m_pData[0] = m.m_pData[0];
		m_pData[1] = m.m_pData[1];
		m_pData[2] = m.m_pData[2];
		m_pData[3] = m.m_pData[3];
		m_pData[4] = m.m_pData[4];
		m_pData[5] = m.m_pData[5];
		m_pData[6] = m.m_pData[6];
		m_pData[7] = m.m_pData[7];
		m_pData[8] = m.m_pData[8];
	}*/
	
	CxDMatrix3(const CxDMatrix3 &m) : CxObject()
	{
		BXIN;
		m_pData[0] = m.m_pData[0];
		m_pData[1] = m.m_pData[1];
		m_pData[2] = m.m_pData[2];
		m_pData[3] = m.m_pData[3];
		m_pData[4] = m.m_pData[4];
		m_pData[5] = m.m_pData[5];
		m_pData[6] = m.m_pData[6];
		m_pData[7] = m.m_pData[7];
		m_pData[8] = m.m_pData[8];
		BXOUT;
	}

	CxDMatrix3(double f)
	{
		BXIN;
		m_pData[0] = f;
		m_pData[1] = f;
		m_pData[2] = f;
		m_pData[3] = f;
		m_pData[4] = f;
		m_pData[5] = f;
		m_pData[6] = f;
		m_pData[7] = f;
		m_pData[8] = f;
		BXOUT;
	}

	CxDMatrix3(double a, double b, double c, double d, double e, double f, double g, double h, double i)
	{
		BXIN;
		m_pData[0] = a;
		m_pData[1] = b;
		m_pData[2] = c;
		m_pData[3] = d;
		m_pData[4] = e;
		m_pData[5] = f;
		m_pData[6] = g;
		m_pData[7] = h;
		m_pData[8] = i;
		BXOUT;
	}

	double &GetAt(int i)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}

	double &GetAt(int i, int j)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::GetAt(int, int): %d, %d...",i,j);
		#endif
		return m_pData[i*3+j];
	}

	double &operator [] (int i)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}

	double GetAt(int i) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}

	double GetAt(int i, int j) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::GetAt(int, int): %d, %d...",i,j);
		#endif
		return m_pData[i*3+j];
	}

	double operator [] (int i) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}

	void operator *= (const double &f)
	{
		BXIN;
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator *= (double)\n");
		#endif
		m_pData[0] *= f;
		m_pData[1] *= f;
		m_pData[2] *= f;
		m_pData[3] *= f;
		m_pData[4] *= f;
		m_pData[5] *= f;
		m_pData[6] *= f;
		m_pData[7] *= f;
		m_pData[8] *= f;
		BXOUT;
	}

	CxDMatrix3 operator * (const CxDMatrix3 &m) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator * (CxDMatrix3)\n");
		#endif
		return CxDMatrix3(
			GetAt(0,0)*m.GetAt(0,0) + GetAt(1,0)*m.GetAt(0,1) + GetAt(2,0)*m.GetAt(0,2),
			GetAt(0,1)*m.GetAt(0,0) + GetAt(1,1)*m.GetAt(0,1) + GetAt(2,1)*m.GetAt(0,2),
			GetAt(0,2)*m.GetAt(0,0) + GetAt(1,2)*m.GetAt(0,1) + GetAt(2,2)*m.GetAt(0,2),
			GetAt(0,0)*m.GetAt(1,0) + GetAt(1,0)*m.GetAt(1,1) + GetAt(2,0)*m.GetAt(1,2),
			GetAt(0,1)*m.GetAt(1,0) + GetAt(1,1)*m.GetAt(1,1) + GetAt(2,1)*m.GetAt(1,2),
			GetAt(0,2)*m.GetAt(1,0) + GetAt(1,2)*m.GetAt(1,1) + GetAt(2,2)*m.GetAt(1,2),
			GetAt(0,0)*m.GetAt(2,0) + GetAt(1,0)*m.GetAt(2,1) + GetAt(2,0)*m.GetAt(2,2),
			GetAt(0,1)*m.GetAt(2,0) + GetAt(1,1)*m.GetAt(2,1) + GetAt(2,1)*m.GetAt(2,2),
			GetAt(0,2)*m.GetAt(2,0) + GetAt(1,2)*m.GetAt(2,1) + GetAt(2,2)*m.GetAt(2,2)
			);
	}

	CxDVector3 operator * (const CxDVector3 &v) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator * (CxDVector3)\n");
		#endif
		return CxDVector3(
			GetAt(0,0)*v[0] + GetAt(1,0)*v[1] + GetAt(2,0)*v[2],
			GetAt(0,1)*v[0] + GetAt(1,1)*v[1] + GetAt(2,1)*v[2],
			GetAt(0,2)*v[0] + GetAt(1,2)*v[1] + GetAt(2,2)*v[2]);
	}

	void RotMat(const CxDVector3 &vec, const double &a)
	{
		BXIN;
		double ca, sa;
		ca = (double)cos(a);
		sa = (double)sin(a);
		GetAt(0,0) = (double)(vec[0]*vec[0]*(1.0-ca) + ca);
		GetAt(1,0) = (double)(vec[0]*vec[1]*(1.0-ca) - vec[2]*sa);
		GetAt(2,0) = (double)(vec[0]*vec[2]*(1.0-ca) + vec[1]*sa);
		GetAt(0,1) = (double)(vec[1]*vec[0]*(1.0-ca) + vec[2]*sa);
		GetAt(1,1) = (double)(vec[1]*vec[1]*(1.0-ca) + ca);
		GetAt(2,1) = (double)(vec[1]*vec[2]*(1.0-ca) - vec[0]*sa);
		GetAt(0,2) = (double)(vec[2]*vec[0]*(1.0-ca) - vec[1]*sa);
		GetAt(1,2) = (double)(vec[2]*vec[1]*(1.0-ca) + vec[0]*sa);
		GetAt(2,2) = (double)(vec[2]*vec[2]*(1.0-ca) + ca);
		BXOUT;
	}

	inline void Unity()
	{ 
		BXIN;
		GetAt(0,0) = 1;
		GetAt(1,0) = 0;
		GetAt(2,0) = 0;
		GetAt(0,1) = 0;
		GetAt(1,1) = 1;
		GetAt(2,1) = 0;
		GetAt(0,2) = 0;
		GetAt(1,2) = 0;
		GetAt(2,2) = 1;
		BXOUT;
	}

	CxDMatrix3 Transpose()
	{
		return CxDMatrix3(
			GetAt(0,0), GetAt(1,0), GetAt(2,0),
			GetAt(0,1), GetAt(1,1), GetAt(2,1),
			GetAt(0,2), GetAt(1,2), GetAt(2,2));
	}

	void Dump() const;

	void Invert();

	void MatUltra(const CxDVector3 &vec1, const CxDVector3 &vec2);

private:
	double m_pData[9];
};


#endif
