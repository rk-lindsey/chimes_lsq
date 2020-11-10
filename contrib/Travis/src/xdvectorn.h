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


#ifndef CXDVECTORN_DEFINED
#define CXDVECTORN_DEFINED

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


class CxDVectorN : public CxObject
{
public:

	inline CxDVectorN()
	{
		m_iDim = 0;
		m_pData = NULL;
	}


	inline ~CxDVectorN()
	{
		if (m_pData != NULL)
		{
			delete[] m_pData;
			m_pData = NULL;
		}
	}


	inline CxDVectorN(const CxDVectorN &v) : CxObject()
	{
		m_iDim = v.m_iDim;
		m_pData = new double[m_iDim];
		memcpy(m_pData,v.m_pData,sizeof(double)*m_iDim);
	}


	inline CxDVectorN & operator = (const CxDVectorN &v)
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator = (const CxDVectorN &)\n");
		#endif
		if (m_pData != NULL)
			delete[] m_pData;
		m_iDim = v.m_iDim;
		m_pData = new double[m_iDim];
		memcpy(m_pData,v.m_pData,sizeof(double)*m_iDim);
		return *this;
	}


	inline CxDVectorN(int i, double f)
	{
		int z;
		m_iDim = i;
		m_pData = new double[m_iDim];
		for (z=0;z<m_iDim;z++)
			m_pData[z] = f;
	}


	inline CxDVectorN(int i)
	{
		m_iDim = i;
		m_pData = new double[m_iDim];
	}


	inline void ZeroVector()
	{
		memset(m_pData,0,sizeof(double)*m_iDim);
	}


	inline int GetDim() const
	{
		return m_iDim;
	}


	inline double &GetAt(int i)
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
		if (i >= m_iDim)
		{
			eprintf("& CxDVectorN::GetAt(int): Boundary Error (%d/%d).\n",i,m_iDim);
			abort();
		}
		#endif
		return m_pData[i];
	}


	inline double &operator [] (int i)
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator [] (int): %d\n",i);
		#endif
		#ifdef DEBUG_ARRAYS
		if ((i < 0) || (i >= m_iDim))
		{
			eprintf("& CxDVectorN::operator [] (int): Boundary Error (%d/%d).\n",i,m_iDim);
			abort();
		}
		#endif
		return m_pData[i];
	}


	inline double GetAt(int i) const
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
		if (i >= m_iDim)
		{
			eprintf("CxDVectorN::GetAt(int): Boundary Error (%d/%d).\n",i,m_iDim);
			abort();
		}
		#endif
		return m_pData[i];
	}


	inline double operator [] (int i) const
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator [] (int): %d\n",i);
		#endif
		#ifdef DEBUG_ARRAYS
		if ((i < 0) || (i >= m_iDim))
		{
			eprintf("CxDVectorN::operator [] (int): Boundary Error (%d/%d).\n",i,m_iDim);
			abort();
		}
		#endif
		return m_pData[i];
	}


	inline CxDVectorN operator + (const CxDVectorN &v) const
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator + (CxDVectorN&)\n");
		#endif
		if (m_iDim != v.m_iDim)
		{
			eprintf("CxDVectorN::operator + (CxDVectorN&): Dimension mismatch (%d vs %d).\n",m_iDim,v.m_iDim);
			abort();
		}
		CxDVectorN r(m_iDim);
		int z;
		for (z=0;z<m_iDim;z++)
			r.m_pData[z] = m_pData[z] + v.m_pData[z];
		return r;
	}


	inline CxDVectorN operator - (const CxDVectorN &v) const
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator - (CxDVectorN&)\n");
		#endif
		if (m_iDim != v.m_iDim)
		{
			eprintf("CxDVectorN::operator - (CxDVectorN&): Dimension mismatch (%d vs %d).\n",m_iDim,v.m_iDim);
			abort();
		}
		CxDVectorN r(m_iDim);
		int z;
		for (z=0;z<m_iDim;z++)
			r.m_pData[z] = m_pData[z] - v.m_pData[z];
		return r;
	}


	inline CxDVectorN operator * (const double &f) const
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator * (double)\n");
		#endif
		CxDVectorN r(m_iDim);
		int z;
		for (z=0;z<m_iDim;z++)
			r.m_pData[z] = m_pData[z] * f;
		return r;
	}


	inline CxDVectorN operator / (const double &f) const
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator / (double)\n");
		#endif
		CxDVectorN r(m_iDim);
		int z;
		for (z=0;z<m_iDim;z++)
			r.m_pData[z] = m_pData[z] / f;
		return r;
	}


	inline void operator += (const CxDVectorN &v)
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator += (CxDVectorN&)\n");
		#endif
		if (m_iDim != v.m_iDim)
		{
			eprintf("CxDVectorN::operator += (CxDVectorN&): Dimension mismatch (%d vs %d).\n",m_iDim,v.m_iDim);
			abort();
		}
		int z;
		for (z=0;z<m_iDim;z++)
			m_pData[z] += v.m_pData[z];
	}


	inline void operator -= (const CxDVectorN &v)
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator -= (CxDVectorN&)\n");
		#endif
		if (m_iDim != v.m_iDim)
		{
			eprintf("CxDVectorN::operator -= (CxDVectorN&): Dimension mismatch (%d vs %d).\n",m_iDim,v.m_iDim);
			abort();
		}
		int z;
		for (z=0;z<m_iDim;z++)
			m_pData[z] -= v.m_pData[z];
	}


	inline void operator *= (const double &f)
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator *= (double)\n");
		#endif
		int z;
		for (z=0;z<m_iDim;z++)
			m_pData[z] *= f;
	}


	inline void operator /= (const double &f)
	{
		#ifdef DEBUG_CDVECTORN
		mprintf("@ CxDVectorN::operator /= (double)\n");
		#endif
		int z;
		for (z=0;z<m_iDim;z++)
			m_pData[z] /= f;
	}


	inline double GetLength() const
	{
		int z;
		double r=0;
		for (z=0;z<m_iDim;z++)
			r += m_pData[z]*m_pData[z];
		return sqrt(r);
	}


	inline double GetLengthSqr() const
	{
		int z;
		double r=0;
		for (z=0;z<m_iDim;z++)
			r += m_pData[z]*m_pData[z];
		return r;
	}


	inline void Normalize()
	{
		BXIN;
		int z;
		double l;
		l = GetLength();
		for (z=0;z<m_iDim;z++)
			m_pData[z] /= l;
		BXOUT;
	}


	void Dump() const;

private:
	double *m_pData;
	int m_iDim;
};


inline CxDVectorN operator * (const double &f, const CxDVectorN &v)
{
	return v*f;
}


inline double DotP(const CxDVectorN &vec1, const CxDVectorN &vec2)
{
	if (vec1.GetDim() != vec2.GetDim())
	{
		eprintf("double DotP(const CxDVectorN &, const CxDVectorN &): Dimension mismatch (%d vs %d).\n",vec1.GetDim(),vec2.GetDim());
		abort();
	}

	int z;
	double f=0;
	for (z=0;z<vec1.GetDim();z++)
		f += vec1[z]*vec2[z];
	return f;
}


inline double Angle(const CxDVectorN &vec1, const CxDVectorN &vec2)
{
	BXIN;
	double t;
	t = DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0f)
		t = 1.0f;
	if (t < -1.0f)
		t = -1.0f;
	return (double)acos(t);
	BXOUT;
}

inline double Angle_Deg(const CxDVectorN &vec1, const CxDVectorN &vec2)
{
	BXIN;
	double t;
	t = DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0f)
		t = 1.0f;
	if (t < -1.0f)
		t = -1.0f;
	t = (double)acos(t);
	BXOUT;
	return (double)fabs(t*180.0f / Pi);
}

inline double VecDist(const CxDVectorN &vec1, const CxDVectorN &vec2)
{
	return (vec1-vec2).GetLength();
}

#endif
