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


#ifndef CXQUATERNION_DEFINED
#define CXQUATERNION_DEFINED

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


class CxQuaternion : public CxObject
{
public:
	inline CxQuaternion() { }

	inline ~CxQuaternion() { }
	
	inline CxQuaternion(const CxQuaternion &v) : CxObject()
	{
		BXIN;
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
		m_pData[3] = v.m_pData[3];
		BXOUT;
	}

	inline CxQuaternion(double f)
	{
		BXIN;
		m_pData[0] = f;
		m_pData[1] = f;
		m_pData[2] = f;
		m_pData[3] = f;
		BXOUT;
	}

	inline CxQuaternion(double a, double b, double c, double d)
	{
		BXIN;
		m_pData[0] = a;
		m_pData[1] = b;
		m_pData[2] = c;
		m_pData[3] = d;
		BXOUT;
	}

	inline double &GetAt(int i)
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}

	inline double &operator [] (int i)
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}

	inline double GetAt(int i) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
		if (i > 3)
		{
			mprintf("CxQuaternion::GetAt(int): Boundary Error (%d/3)...",i);
			abort();
		}
		#endif
		return m_pData[i];
	}

	inline double operator [] (int i) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}

	inline CxQuaternion operator + (const CxQuaternion &v) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator + (CxQuaternion&)\n");
		#endif
		return CxQuaternion(m_pData[0]+v.m_pData[0],m_pData[1]+v.m_pData[1],m_pData[2]+v.m_pData[2],m_pData[3]+v.m_pData[3]);
	}

	inline CxQuaternion operator - (const CxQuaternion &v) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator - (CxQuaternion&)\n");
		#endif
		return CxQuaternion(m_pData[0]-v.m_pData[0],m_pData[1]-v.m_pData[1],m_pData[2]-v.m_pData[2],m_pData[3]-v.m_pData[3]);
	}

	inline CxQuaternion operator * (const CxQuaternion &v) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator + (CxQuaternion&)\n");
		#endif
		return CxQuaternion(
			m_pData[0]*v.m_pData[0] - m_pData[1]*v.m_pData[1] - m_pData[2]*v.m_pData[2] - m_pData[3]*v.m_pData[3],
			m_pData[0]*v.m_pData[1] + m_pData[1]*v.m_pData[0] + m_pData[2]*v.m_pData[3] - m_pData[3]*v.m_pData[2],
			m_pData[0]*v.m_pData[2] - m_pData[1]*v.m_pData[3] + m_pData[2]*v.m_pData[0] + m_pData[3]*v.m_pData[1],
			m_pData[0]*v.m_pData[3] + m_pData[1]*v.m_pData[2] - m_pData[2]*v.m_pData[1] + m_pData[3]*v.m_pData[0] );
	}

	inline CxQuaternion operator * (const double &f) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator * (double)\n");
		#endif
		return CxQuaternion(m_pData[0]*f,m_pData[1]*f,m_pData[2]*f,m_pData[3]*f);
	}

	inline CxQuaternion operator / (const double &f) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator / (double)\n");
		#endif
		return CxQuaternion(m_pData[0]/f,m_pData[1]/f,m_pData[2]/f,m_pData[3]/f);
	}

	inline void operator += (const CxQuaternion &v)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator += (CxQuaternion&)\n");
		#endif
		m_pData[0] += v.m_pData[0];
		m_pData[1] += v.m_pData[1];
		m_pData[2] += v.m_pData[2];
		m_pData[3] += v.m_pData[3];
		BXOUT;
	}

	inline void operator -= (const CxQuaternion &v)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator -= (CxQuaternion&)\n");
		#endif
		m_pData[0] -= v.m_pData[0];
		m_pData[1] -= v.m_pData[1];
		m_pData[2] -= v.m_pData[2];
		m_pData[3] -= v.m_pData[3];
		BXOUT;
	}

	inline void operator *= (const double &f)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator *= (double)\n");
		#endif
		m_pData[0] *= f;
		m_pData[1] *= f;
		m_pData[2] *= f;
		m_pData[3] *= f;
		BXOUT;
	}

	inline void operator /= (const double &f)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator /= (double)\n");
		#endif
		m_pData[0] /= f;
		m_pData[1] /= f;
		m_pData[2] /= f;
		m_pData[3] /= f;
		BXOUT;
	}

	inline double GetLength() const
	{
		return (double)sqrt(m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2]+m_pData[3]*m_pData[3]);
	}

	inline double GetLengthSqr() const
	{
		return m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2]+m_pData[3]*m_pData[3];
	}

	inline CxQuaternion Conjugate() const
	{
		return CxQuaternion(m_pData[0],-m_pData[1],-m_pData[2],-m_pData[3]);
	}

	inline void BuildRotation(CxDVector3 axis, double angle)
	{
		CxDVector3 tmp;
		tmp = axis;
		tmp.Normalize();
		m_pData[0] = cos(angle/2.0);
		m_pData[1] = sin(angle/2.0) * axis[0];
		m_pData[2] = sin(angle/2.0) * axis[1];
		m_pData[3] = sin(angle/2.0) * axis[2];
	}

	inline CxVector3 Transform(const CxVector3 &v) const
	{
		CxQuaternion out;
		out = this->Conjugate() * CxQuaternion(0,v[0],v[1],v[2]) * (*this);
		return CxVector3((float)out[1],(float)out[2],(float)out[3]);
	}

	inline CxDVector3 Transform(const CxDVector3 &v) const
	{
		CxQuaternion out;
		out = this->Conjugate() * CxQuaternion(0,v[0],v[1],v[2]) * (*this);
		return CxDVector3(out[1],out[2],out[3]);
	}

	inline void Normalize()
	{
		BXIN;
		double l;
		l = GetLength();
		m_pData[0] /= l;
		m_pData[1] /= l;
		m_pData[2] /= l;
		m_pData[3] /= l;
		BXOUT;
	}

	inline void Unity()
	{
		BXIN;
		m_pData[0] = 1.0;
		m_pData[1] = 0;
		m_pData[2] = 0;
		m_pData[3] = 0;
		BXOUT;
	}


private:
	double m_pData[4];
};


#endif
