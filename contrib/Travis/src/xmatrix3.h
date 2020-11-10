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


#ifndef CXMATRIX3_DEFINED
#define CXMATRIX3_DEFINED

#include "tools.h"
#include "xobject.h"
#include "xvector3.h"
#include "backtrace.h"


class CxMatrix3 : public CxObject
{
public:
	void Invert();
	float Det();
	CxMatrix3() { }

	~CxMatrix3() { }

/*	CxMatrix3(CxMatrix3 m)
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
	
	CxMatrix3(const CxMatrix3 &m) : CxObject()
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

	CxMatrix3(float f)
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

	CxMatrix3(float a, float b, float c, float d, float e, float f, float g, float h, float i)
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

	float &GetAt(int i)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}

	float &GetAt(int i, int j)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::GetAt(int, int): %d, %d...",i,j);
		#endif
		return m_pData[i*3+j];
	}

	float &operator [] (int i)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}

	float &operator () (int i, int j)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::operator () (int, int): %d, %d\n",i,j);
		#endif
		return GetAt(i,j);
	}

	float GetAt(int i) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}

	float GetAt(int i, int j) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::GetAt(int, int): %d, %d...",i,j);
		#endif
		return m_pData[i*3+j];
	}

	float operator [] (int i) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}

	void operator *= (const float &f)
	{
		BXIN;
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::operator *= (float)\n");
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

	CxMatrix3 operator * (const CxMatrix3 &m) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::operator * (CxMatrix3)\n");
		#endif
		return CxMatrix3(
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

	CxVector3 operator * (const CxVector3 &v) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxMatrix3::operator * (CxVector3)\n");
		#endif
		return CxVector3(
			GetAt(0,0)*v[0] + GetAt(1,0)*v[1] + GetAt(2,0)*v[2],
			GetAt(0,1)*v[0] + GetAt(1,1)*v[1] + GetAt(2,1)*v[2],
			GetAt(0,2)*v[0] + GetAt(1,2)*v[1] + GetAt(2,2)*v[2]);
	}

	void RotMat(const CxVector3 &vec, const float &a)
	{
		BXIN;
		float ca, sa;
		ca = (float)cos(a);
		sa = (float)sin(a);
		GetAt(0,0) = (float)(vec[0]*vec[0]*(1.0f-ca) + ca);
		GetAt(1,0) = (float)(vec[0]*vec[1]*(1.0f-ca) - vec[2]*sa);
		GetAt(2,0) = (float)(vec[0]*vec[2]*(1.0f-ca) + vec[1]*sa);
		GetAt(0,1) = (float)(vec[1]*vec[0]*(1.0f-ca) + vec[2]*sa);
		GetAt(1,1) = (float)(vec[1]*vec[1]*(1.0f-ca) + ca);
		GetAt(2,1) = (float)(vec[1]*vec[2]*(1.0f-ca) - vec[0]*sa);
		GetAt(0,2) = (float)(vec[2]*vec[0]*(1.0f-ca) - vec[1]*sa);
		GetAt(1,2) = (float)(vec[2]*vec[1]*(1.0f-ca) + vec[0]*sa);
		GetAt(2,2) = (float)(vec[2]*vec[2]*(1.0f-ca) + ca);
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

	CxMatrix3 Transpose()
	{
		return CxMatrix3(
			GetAt(0,0), GetAt(1,0), GetAt(2,0),
			GetAt(0,1), GetAt(1,1), GetAt(2,1),
			GetAt(0,2), GetAt(1,2), GetAt(2,2));
	}

	void Dump() const;

	void MatUltra(const CxVector3 &vec1, const CxVector3 &vec2);

private:
	float m_pData[9];
};


#endif
