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


#ifndef CXVECTOR3_DEFINED
#define CXVECTOR3_DEFINED


#include "tools.h"
#include "xobject.h"
#include "backtrace.h"
#include "xdvector3.h"


class CxVector3 : public CxObject
{
public:
	inline CxVector3()
	{
	}


	inline ~CxVector3()
	{
	}


/*	CxVector3(CxVector3 v)
	{
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
	}*/
	

	inline CxVector3(const CxVector3 &v) : CxObject()
	{
		BXIN;
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
		BXOUT;
	}


	inline CxVector3(const CxDVector3 &v)
	{
		BXIN;
		m_pData[0] = (float)v.GetAt(0);
		m_pData[1] = (float)v.GetAt(1);
		m_pData[2] = (float)v.GetAt(2);
		BXOUT;
	}


	inline CxVector3(float f)
	{
		BXIN;
		m_pData[0] = f;
		m_pData[1] = f;
		m_pData[2] = f;
		BXOUT;
	}


	inline CxVector3(float x, float y, float z)
	{
		BXIN;
		m_pData[0] = x;
		m_pData[1] = y;
		m_pData[2] = z;
		BXOUT;
	}


	inline float &GetAt(int i)
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}


	inline float &operator [] (int i)
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	inline float GetAt(int i) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
		if (i > 2)
		{
			mprintf("CxVector3::GetAt(int): Boundary Error (%d/3)...",i);
			abort();
		}
		#endif
		return m_pData[i];
	}


	inline float operator [] (int i) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	inline CxVector3 operator + (const CxVector3 &v) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator + (CxVector3&)\n");
		#endif
		return CxVector3(m_pData[0]+v.m_pData[0],m_pData[1]+v.m_pData[1],m_pData[2]+v.m_pData[2]);
	}


	inline CxVector3 operator - (const CxVector3 &v) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator - (CxVector3&)\n");
		#endif
		return CxVector3(m_pData[0]-v.m_pData[0],m_pData[1]-v.m_pData[1],m_pData[2]-v.m_pData[2]);
	}


	inline CxVector3 operator * (const float &f) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator * (float)\n");
		#endif
		return CxVector3(m_pData[0]*f,m_pData[1]*f,m_pData[2]*f);
	}


	inline CxVector3 operator / (const float &f) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator / (float)\n");
		#endif
		return CxVector3(m_pData[0]/f,m_pData[1]/f,m_pData[2]/f);
	}


	inline void operator += (const CxVector3 &v)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator += (CxVector3&)\n");
		#endif
		m_pData[0] += v.m_pData[0];
		m_pData[1] += v.m_pData[1];
		m_pData[2] += v.m_pData[2];
		BXOUT;
	}


	inline void operator -= (const CxVector3 &v)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator -= (CxVector3&)\n");
		#endif
		m_pData[0] -= v.m_pData[0];
		m_pData[1] -= v.m_pData[1];
		m_pData[2] -= v.m_pData[2];
		BXOUT;
	}


	inline void operator *= (const float &f)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator *= (float)\n");
		#endif
		m_pData[0] *= f;
		m_pData[1] *= f;
		m_pData[2] *= f;
		BXOUT;
	}


	inline void operator /= (const float &f)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxVector3::operator /= (float)\n");
		#endif
		m_pData[0] /= f;
		m_pData[1] /= f;
		m_pData[2] /= f;
		BXOUT;
	}


	inline float GetLength() const
	{
		return (float)sqrt(m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2]);
	}


	inline float GetLengthSqr() const
	{
		return (float)(m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2]);
	}


	inline void Normalize()
	{
		BXIN;
		float l;
		l = GetLength();
		m_pData[0] /= l;
		m_pData[1] /= l;
		m_pData[2] /= l;
		BXOUT;
	}


	inline void Chop(float d)
	{
		if ((m_pData[0] != 0) && (fabsf(m_pData[0]) < d))
			m_pData[0] = 0;
		if ((m_pData[1] != 0) && (fabsf(m_pData[1]) < d))
			m_pData[1] = 0;
		if ((m_pData[2] != 0) && (fabsf(m_pData[2]) < d))
			m_pData[2] = 0;
	}


	void PointRoot(const CxVector3 &vec1, const CxVector3 &vec2, const CxVector3 &point);


	void Dump() const;


private:
	float m_pData[3];
};


float Dihedral(const CxVector3 &vec1, const CxVector3 &vec2, const CxVector3 &norm, bool absolute);


float Dihedral2(const CxVector3 &vec1, const CxVector3 &vec2, const CxVector3 &norm, bool absolute);


inline void Swap(CxVector3 &vec1, CxVector3 &vec2)
{
	BXIN;
	CxVector3 t;
	t = vec1;
	vec1 = vec2;
	vec2 = t;
	BXOUT;
}


inline CxVector3 operator * (const float &f, const CxVector3 &v)
{
	return v*f;
}


inline CxVector3 CrossP(const CxVector3 &vec1, const CxVector3 &vec2)
{
	return CxVector3(vec1[1]*vec2[2] - vec1[2]*vec2[1],
		vec1[2]*vec2[0] - vec1[0]*vec2[2],
		vec1[0]*vec2[1] - vec1[1]*vec2[0]);
}


inline float DotP(const CxVector3 &vec1, const CxVector3 &vec2)
{
	return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
}


inline float Angle(const CxVector3 &vec1, const CxVector3 &vec2)
{
	BXIN;
	float t;

	if ((vec1.GetLength() == 0) || (vec2.GetLength() == 0))
	{
		mprintf("\nAngle(): Warning: Indeterminate angle between ( %f | %f | %f ) and ( %f | %f | %f ).  ",vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2]);
		return 0;
	}

	t = DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0f)
		t = 1.0f;
	if (t < -1.0f)
		t = -1.0f;
	return (float)acos(t);
	BXOUT;
}


inline float Angle_Deg(const CxVector3 &vec1, const CxVector3 &vec2)
{
	BXIN;
	float t;

	if ((vec1.GetLength() == 0) || (vec2.GetLength() == 0))
	{
		mprintf("\nAngle_Deg(): Warning: Indeterminate angle between ( %f | %f | %f ) and ( %f | %f | %f ).  ",vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2]);
		return 0;
	}

	t = DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0f)
		t = 1.0f;
	if (t < -1.0f)
		t = -1.0f;
	t = (float)acos(t);
	BXOUT;
	return (float)fabs(t*180.0f / Pi);
}


inline float VecDist(const CxVector3 &vec1, const CxVector3 &vec2)
{
	return (float)sqrt((vec1[0]-vec2[0])*(vec1[0]-vec2[0]) + (vec1[1]-vec2[1])*(vec1[1]-vec2[1]) + (vec1[2]-vec2[2])*(vec1[2]-vec2[2]));
}


inline CxVector3 Normalize(const CxVector3 &vec)
{
	float tf;

	tf = vec.GetLength();

	return CxVector3(vec[0]/tf,vec[1]/tf,vec[2]/tf);
}


CxVector3 PointFromRAD(CxVector3 r1, CxVector3 r2, CxVector3 r3, float r, float a, float d);


#endif
