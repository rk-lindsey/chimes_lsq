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

#ifndef _3DF_H
#define _3DF_H

#include "xobject.h"
#include "df.h"
#include "2df.h"
#include "xmemfile.h"
#include "xstring.h"


extern int g_iCubeXStride;
extern int g_iCubeYStride;
extern int g_iCubeZStride;
extern int g_iCubeXMismatch;
extern int g_iCubeYMismatch;
extern int g_iCubeZMismatch;
extern float g_fBoxX;
extern float g_fBoxY;
extern float g_fBoxZ;
extern float g_fCubeXStep;
extern float g_fCubeYStep;
extern float g_fCubeZStep;


template<typename T>
class C3DF : public CxObject
{
public:

	inline void AddToBin(const CxVector3 &vec)
	{
		AddToBin(vec[0],vec[1],vec[2]);
	}


	inline void AddToBin_Sphere(const CxVector3 &vec, double r)
	{
		AddToBin_Sphere(vec[0],vec[1],vec[2],r);
	}


	inline void AddToBin_SphereWrap(const CxVector3 &vec, double r)
	{
		AddToBin_SphereWrap(vec[0],vec[1],vec[2],r);
	}


	inline void AddToBin_Single(double x, double y, double z, T val)
	{
		int ix, iy, iz;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]) || (z < m_fMinVal[2]) || (z > m_fMaxVal[2]))
			return;

		ix = (int)floor((x-m_fMinVal[0])*m_fFac[0]);
		if (ix < 0)
			ix = 0;
		if (ix >= m_iRes[0])
			ix = m_iRes[0] - 1;
		iy = (int)floor((y-m_fMinVal[1])*m_fFac[1]);
		if (iy < 0)
			iy = 0;
		if (iy >= m_iRes[1])
			iy = m_iRes[1] - 1;
		iz = (int)floor((z-m_fMinVal[2])*m_fFac[2]);
		if (iz < 0)
			iz = 0;
		if (iz >= m_iRes[2])
			iz = m_iRes[2] - 1;

		m_pBin[ iz*m_iResXY + iy*m_iRes[0] + ix ] += val;
	}


	inline T GetBinValue_Single(double x, double y, double z) const
	{
		int ix, iy, iz;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]) || (z < m_fMinVal[2]) || (z > m_fMaxVal[2]))
			return -1;

		ix = (int)floor((x-m_fMinVal[0])*m_fFac[0]);
		if (ix < 0)
			ix = 0;
		if (ix >= m_iRes[0])
			ix = m_iRes[0] - 1;
		iy = (int)floor((y-m_fMinVal[1])*m_fFac[1]);
		if (iy < 0)
			iy = 0;
		if (iy >= m_iRes[1])
			iy = m_iRes[1] - 1;
		iz = (int)floor((z-m_fMinVal[2])*m_fFac[2]);
		if (iz < 0)
			iz = 0;
		if (iz >= m_iRes[2])
			iz = m_iRes[2] - 1;

		return m_pBin[ iz*m_iResXY + iy*m_iRes[0] + ix ];
	}


	inline void SetBinValue_Single(double x, double y, double z, T val)
	{
		int ix, iy, iz;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]) || (z < m_fMinVal[2]) || (z > m_fMaxVal[2]))
			return;

		ix = (int)floor((x-m_fMinVal[0])*m_fFac[0]);
		if (ix < 0)
			ix = 0;
		if (ix >= m_iRes[0])
			ix = m_iRes[0] - 1;
		iy = (int)floor((y-m_fMinVal[1])*m_fFac[1]);
		if (iy < 0)
			iy = 0;
		if (iy >= m_iRes[1])
			iy = m_iRes[1] - 1;
		iz = (int)floor((z-m_fMinVal[2])*m_fFac[2]);
		if (iz < 0)
			iz = 0;
		if (iz >= m_iRes[2])
			iz = m_iRes[2] - 1;

		m_pBin[ iz*m_iResXY + iy*m_iRes[0] + ix ] = val;
	}


	double m_fMinVal[3];
	double m_fMaxVal[3];

	CDF *m_pChannels[3];
	C2DF *m_p2DF[3];

	T m_fMinEntry;
	T m_fMaxEntry;
	char *m_sLabelX;
	char *m_sLabelY;
	char *m_sLabelZ;

	double *m_fCountX;
	double *m_fCountY;
	double *m_fCountZ;

	int m_iRes[3];
	int m_iHistogramRes;
	long m_iResXY;
	double m_fBinEntries;
	double m_fSkipEntries;

	T *m_pBin;
	double *m_pHistogram;
	double m_fFac[3];



/******************************************************************************************/

	C3DF()
	{
		m_sLabelX = NULL;
		m_sLabelY = NULL;
		m_sLabelZ = NULL;
		m_pBin = NULL;
		m_fCountX = NULL;
		m_fCountY = NULL;
		m_fCountZ = NULL;
		m_iHistogramRes = 0;
		m_pHistogram = NULL;
	}


	~C3DF()
	{
		if (m_sLabelX != NULL)
		{
			delete[] m_sLabelX;
			m_sLabelX = NULL;
		}
		if (m_sLabelY != NULL)
		{
			delete[] m_sLabelY;
			m_sLabelY = NULL;
		}
		if (m_sLabelZ != NULL)
		{
			delete[] m_sLabelZ;
			m_sLabelZ = NULL;
		}
		if (m_pBin != NULL)
		{
			delete[] m_pBin;
			m_pBin = NULL;
		}
		if (m_fCountX != NULL)
		{
			delete[] m_fCountX;
			m_fCountX = NULL;
		}
		if (m_fCountY != NULL)
		{
			delete[] m_fCountY;
			m_fCountY = NULL;
		}
		if (m_fCountZ != NULL)
		{
			delete[] m_fCountZ;
			m_fCountZ = NULL;
		}
	}


	void Create()
	{
		BTIN;
		int z;
		m_fBinEntries = 0;
		m_fSkipEntries = 0;

		if (m_fCountX != NULL) delete[] m_fCountX;
		try { m_fCountX = new double[m_iRes[0]]; } catch(...) { m_fCountX = NULL; }
		if (m_fCountX == NULL) NewException((double)m_iRes[0]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		if (m_fCountY != NULL) delete[] m_fCountY;
		try { m_fCountY = new double[m_iRes[1]]; } catch(...) { m_fCountY = NULL; }
		if (m_fCountY == NULL) NewException((double)m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		if (m_fCountZ != NULL) delete[] m_fCountZ;
		try { m_fCountZ = new double[m_iRes[2]]; } catch(...) { m_fCountZ = NULL; }
		if (m_fCountZ == NULL) NewException((double)m_iRes[2]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		for (z=0;z<m_iRes[0];z++)
			m_fCountX[z] = 0;
		for (z=0;z<m_iRes[1];z++)
			m_fCountY[z] = 0;
		for (z=0;z<m_iRes[2];z++)
			m_fCountZ[z] = 0;

	//	mprintf("    C3DF: Trying to reserve %s of memory...\n",FormatBytes((double)m_iRes[0]*m_iRes[1]*m_iRes[2]*sizeof(double)));

		if (m_pBin != NULL) delete[] m_pBin;
		try { m_pBin = new T[m_iRes[0]*m_iRes[1]*m_iRes[2]]; } catch(...) { m_pBin = NULL; }
		if (m_pBin == NULL) NewException((double)m_iRes[0]*m_iRes[1]*m_iRes[2]*sizeof(T),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if (m_iHistogramRes != 0)
		{
			if (m_pHistogram != NULL) delete[] m_pHistogram;
			try { m_pHistogram = new double[m_iHistogramRes]; } catch(...) { m_pHistogram = NULL; }
			if (m_pHistogram == NULL) NewException((double)m_iHistogramRes*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}

		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			m_pBin[z] = 0;

		m_fFac[0] = m_iRes[0] / (m_fMaxVal[0]-m_fMinVal[0]);
		m_fFac[1] = m_iRes[1] / (m_fMaxVal[1]-m_fMinVal[1]);
		m_fFac[2] = m_iRes[2] / (m_fMaxVal[2]-m_fMinVal[2]);
		m_iResXY = m_iRes[0] * m_iRes[1];

		BTOUT;
	}


	void AddToBin(double x, double y, double z)
	{
		BXIN;
		T rx, ry, rz;
		int ix, iy, iz;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]) || (z < m_fMinVal[2]) || (z > m_fMaxVal[2]))
		{
			m_fSkipEntries++;
			BXOUT;
			return;
		}
		m_fBinEntries++;

		rx = (T)((x-m_fMinVal[0])*m_fFac[0] - 0.5);
		ix = (int)floor(rx);
		if (ix < 0)
		{
			ix = 0;
			rx = 0;
		} else if (ix > m_iRes[0]-2)
		{
			ix = m_iRes[0]-2;
			rx = 1.0;
		} else
			rx -= ix;

		ry = (T)((y-m_fMinVal[1])*m_fFac[1] - 0.5);
		iy = (int)floor(ry);
		if (iy < 0)
		{
			iy = 0;
			ry = 0;
		} else if (iy > m_iRes[1]-2)
		{
			iy = m_iRes[1]-2;
			ry = 1.0;
		} else
			ry -= iy;

		rz = (T)((z-m_fMinVal[2])*m_fFac[2] - 0.5);
		iz = (int)floor(rz);
		if (iz < 0)
		{
			iz = 0;
			rz = 0;
		} else if (iz > m_iRes[2]-2)
		{
			iz = m_iRes[2]-2;
			rz = 1.0;
		} else
			rz -= iz;

		m_pBin[ iz    * m_iResXY +  iy    * m_iRes[0] + ix    ] += ((T)(1.0)-rx) * ((T)(1.0)-ry) * ((T)(1.0)-rz);
		m_pBin[ iz    * m_iResXY +  iy    * m_iRes[0] + ix + 1] +=           rx  * ((T)(1.0)-ry) * ((T)(1.0)-rz);
		m_pBin[ iz    * m_iResXY + (iy+1) * m_iRes[0] + ix    ] += ((T)(1.0)-rx) *           ry  * ((T)(1.0)-rz);
		m_pBin[ iz    * m_iResXY + (iy+1) * m_iRes[0] + ix + 1] +=           rx  *           ry  * ((T)(1.0)-rz);
		m_pBin[(iz+1) * m_iResXY +  iy    * m_iRes[0] + ix    ] += ((T)(1.0)-rx) * ((T)(1.0)-ry) *           rz;
		m_pBin[(iz+1) * m_iResXY +  iy    * m_iRes[0] + ix + 1] +=           rx  * ((T)(1.0)-ry) *           rz;
		m_pBin[(iz+1) * m_iResXY + (iy+1) * m_iRes[0] + ix    ] += ((T)(1.0)-rx) *           ry  *           rz;
		m_pBin[(iz+1) * m_iResXY + (iy+1) * m_iRes[0] + ix + 1] +=           rx  *           ry  *           rz;

		BXOUT;
	}


	void AddToBin_Sphere(double x, double y, double z, double r)
	{
		BXIN;
		int ix, iy, iz, k, kq, zx, zy, zz, zzq, zyq, tz, ty;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]) || (z < m_fMinVal[2]) || (z > m_fMaxVal[2]))
		{
			m_fSkipEntries++;
			BXOUT;
			return;
		}
		m_fBinEntries++;

		ix = (int)floor((x-m_fMinVal[0])*m_fFac[0]);
		iy = (int)floor((y-m_fMinVal[1])*m_fFac[1]);
		iz = (int)floor((z-m_fMinVal[2])*m_fFac[2]);

		k = (int)ceil(r * m_fFac[0] - 1.0);
		kq = k * k;

		for (zz=-k;zz<=k;zz++)
		{
			if (zz+iz < 0)
				continue;
			if (zz+iz >= m_iRes[2])
				continue;
			tz = (zz+iz) * m_iResXY;
			zzq = zz * zz;
			for (zy=-k;zy<=k;zy++)
			{
				if (zy+iy < 0)
					continue;
				if (zy+iy >= m_iRes[1])
					continue;
				ty = ix + tz + (zy+iy) * m_iRes[0];
				zyq = zzq + zy*zy;
				for (zx=-k;zx<=k;zx++)
				{
					if (zx+ix < 0)
						continue;
					if (zx+ix >= m_iRes[0])
						continue;
					if (zyq + zx*zx <= kq)
						m_pBin[ty+zx]++;
				}
			}
		}

		BXOUT;
	}


	void AddToBin_SphereWrap(double x, double y, double z, double r)
	{
		BXIN;
		int ix, iy, iz, k, kq, zx, zy, zz, zzq, zyq, tz, ty;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]) || (z < m_fMinVal[2]) || (z > m_fMaxVal[2]))
		{
			m_fSkipEntries++;
			BXOUT;
			return;
		}
		m_fBinEntries++;

		ix = (int)floor((x-m_fMinVal[0])*m_fFac[0]);
		iy = (int)floor((y-m_fMinVal[1])*m_fFac[1]);
		iz = (int)floor((z-m_fMinVal[2])*m_fFac[2]);

		k = (int)ceil(r * m_fFac[0] - 1.0);
		kq = k * k;

		for (zz=-k;zz<=k;zz++)
		{
			if (zz+iz < 0)
				tz = (zz+iz+m_iRes[2]) * m_iResXY;
			else if (zz+iz >= m_iRes[2])
				tz = (zz+iz-m_iRes[2]) * m_iResXY;
			else tz = (zz+iz) * m_iResXY;

			zzq = zz * zz;

			for (zy=-k;zy<=k;zy++)
			{
				if (zy+iy < 0)
					ty = tz + (zy+iy+m_iRes[1]) * m_iRes[0];
				else if (zy+iy >= m_iRes[1])
					ty = tz + (zy+iy-m_iRes[1]) * m_iRes[0];
				else ty = tz + (zy+iy) * m_iRes[0];

				zyq = zzq + zy*zy;

				for (zx=-k;zx<=k;zx++)
				{
					if (zyq + zx*zx > kq)
						continue;

					if (zx+ix < 0)
						m_pBin[ty+ix+zx+m_iRes[0]]++;
					else if (zx+ix >= m_iRes[0])
						m_pBin[ty+ix+zx-m_iRes[0]]++;
					else m_pBin[ty+ix+zx]++;
				}
			}
		}

		BXOUT;
	}


	void AddToBin_IntZ(double x, double y, int z)
	{
		BXIN;
		double rx, ry;
		int ix, iy;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]))
		{
			m_fSkipEntries++;
			BXOUT;
			return;
		}
		m_fBinEntries++;

		rx = (x-m_fMinVal[0])*m_fFac[0] - 0.5;
		ix = (int)floor(rx);
		if (ix < 0)
		{
			ix = 0;
			rx = 0;
		} else if (ix > m_iRes[0]-2)
		{
			ix = m_iRes[0]-2;
			rx = 1.0;
		} else
			rx -= ix;

		ry = (y-m_fMinVal[1])*m_fFac[1] - 0.5;
		iy = (int)floor(ry);
		if (iy < 0)
		{
			iy = 0;
			ry = 0;
		} else if (iy > m_iRes[1]-2)
		{
			iy = m_iRes[1]-2;
			ry = 1.0;
		} else
			ry -= iy;

		m_pBin[z * m_iResXY +  iy    * m_iRes[0] + ix    ] += (1.0-rx) * (1.0-ry);
		m_pBin[z * m_iResXY +  iy    * m_iRes[0] + ix + 1] +=      rx  * (1.0-ry);
		m_pBin[z * m_iResXY + (iy+1) * m_iRes[0] + ix    ] += (1.0-rx) *      ry ;
		m_pBin[z * m_iResXY + (iy+1) * m_iRes[0] + ix + 1] +=      rx  *      ry ;

		BXOUT;
	}


	T GetValue(double x, double y, double z)
	{
		BXIN;
		T rx, ry, rz, r;
		int ix, iy, iz;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]) || (z < m_fMinVal[2]) || (z > m_fMaxVal[2]))
		{
			BXOUT;
			return 0;
		}

		rx = (T)((x-m_fMinVal[0])*m_fFac[0] - 0.5);
		ix = (int)floor(rx);
		if (ix < 0)
		{
			ix = 0;
			rx = 0;
		} else if (ix > m_iRes[0]-2)
		{
			ix = m_iRes[0]-2;
			rx = 1.0;
		} else
			rx -= ix;

		ry = (T)((y-m_fMinVal[1])*m_fFac[1] - 0.5);
		iy = (int)floor(ry);
		if (iy < 0)
		{
			iy = 0;
			ry = 0;
		} else if (iy > m_iRes[1]-2)
		{
			iy = m_iRes[1]-2;
			ry = 1.0;
		} else
			ry -= iy;

		rz = (T)((z-m_fMinVal[2])*m_fFac[2] - 0.5);
		iz = (int)floor(rz);
		if (iz < 0)
		{
			iz = 0;
			rz = 0;
		} else if (iz > m_iRes[2]-2)
		{
			iz = m_iRes[2]-2;
			rz = 1.0;
		} else
			rz -= iz;

		r  = m_pBin[ iz    * m_iResXY +  iy    * m_iRes[0] + ix    ] * ((T)(1.0)-rx) * ((T)(1.0)-ry) * ((T)(1.0)-rz);
		r += m_pBin[ iz    * m_iResXY +  iy    * m_iRes[0] + ix + 1] *           rx  * ((T)(1.0)-ry) * ((T)(1.0)-rz);
		r += m_pBin[ iz    * m_iResXY + (iy+1) * m_iRes[0] + ix    ] * ((T)(1.0)-rx) *           ry  * ((T)(1.0)-rz);
		r += m_pBin[ iz    * m_iResXY + (iy+1) * m_iRes[0] + ix + 1] *           rx  *           ry  * ((T)(1.0)-rz);
		r += m_pBin[(iz+1) * m_iResXY +  iy    * m_iRes[0] + ix    ] * ((T)(1.0)-rx) * ((T)(1.0)-ry) *           rz;
		r += m_pBin[(iz+1) * m_iResXY +  iy    * m_iRes[0] + ix + 1] *           rx  * ((T)(1.0)-ry) *           rz;
		r += m_pBin[(iz+1) * m_iResXY + (iy+1) * m_iRes[0] + ix    ] * ((T)(1.0)-rx) *           ry  *           rz;
		r += m_pBin[(iz+1) * m_iResXY + (iy+1) * m_iRes[0] + ix + 1] *           rx  *           ry  *           rz;

		BXOUT;
		return r;
	}


	bool IsZero(double x, double y, double z)
	{
		double rx, ry, rz;
		int ix, iy, iz;

		if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]) || (y < m_fMinVal[1]) || (y > m_fMaxVal[1]) || (z < m_fMinVal[2]) || (z > m_fMaxVal[2]))
			return false;

		rx = (x-m_fMinVal[0])*m_fFac[0] - 0.5;
		ix = (int)floor(rx);
		if (ix < 0)
			ix = 0;
		else if (ix > m_iRes[0]-2)
			ix = m_iRes[0]-2;

		ry = (y-m_fMinVal[1])*m_fFac[1] - 0.5;
		iy = (int)floor(ry);
		if (iy < 0)
			iy = 0;
		else if (iy > m_iRes[1]-2)
			iy = m_iRes[1]-2;

		rz = (z-m_fMinVal[2])*m_fFac[2] - 0.5;
		iz = (int)floor(rz);
		if (iz < 0)
			iz = 0;
		else if (iz > m_iRes[2]-2)
			iz = m_iRes[2]-2;

		if (m_pBin[ iz    * m_iResXY +  iy    * m_iRes[0] + ix    ] != 0) return false;
		if (m_pBin[ iz    * m_iResXY +  iy    * m_iRes[0] + ix + 1] != 0) return false;
		if (m_pBin[ iz    * m_iResXY + (iy+1) * m_iRes[0] + ix    ] != 0) return false;
		if (m_pBin[ iz    * m_iResXY + (iy+1) * m_iRes[0] + ix + 1] != 0) return false;
		if (m_pBin[(iz+1) * m_iResXY +  iy    * m_iRes[0] + ix    ] != 0) return false;
		if (m_pBin[(iz+1) * m_iResXY +  iy    * m_iRes[0] + ix + 1] != 0) return false;
		if (m_pBin[(iz+1) * m_iResXY + (iy+1) * m_iRes[0] + ix    ] != 0) return false;
		if (m_pBin[(iz+1) * m_iResXY + (iy+1) * m_iRes[0] + ix + 1] != 0) return false;

		return true;
	}


	T GetValue(const CxVector3 &vec)
	{
		return GetValue(vec[0],vec[1],vec[2]);
	}


	void SetLabelX(const char *s)
	{
		if (m_sLabelX != NULL)
			delete[] m_sLabelX;

		try { m_sLabelX = new char[strlen(s)+1]; } catch(...) { m_sLabelX = NULL; }
		if (m_sLabelX == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		strcpy(m_sLabelX,s);
	}


	void SetLabelY(const char *s)
	{
		if (m_sLabelY != NULL)
			delete[] m_sLabelY;

		try { m_sLabelY = new char[strlen(s)+1]; } catch(...) { m_sLabelY = NULL; }
		if (m_sLabelY == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		strcpy(m_sLabelY,s);
	}


	void SetLabelZ(const char *s)
	{
		if (m_sLabelZ != NULL)
			delete[] m_sLabelZ;

		try { m_sLabelZ = new char[strlen(s)+1]; } catch(...) { m_sLabelZ = NULL; }
		if (m_sLabelZ == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		strcpy(m_sLabelZ,s);
	}


	void WritePLT(const char *s1, const char *s2, const char *s3, bool sdf)
	{
		BTIN;
		FILE *a;
		unsigned int i;
		int x, y, z;
		float ff;
//		char buf[256];
		CxString buf;
		
//		strcpy(buf,s1);
//		strcat(buf,s2);
//		strcat(buf,s3);
		buf.strcpy(s1);
		buf.strcat(s2);
		buf.strcat(s3);

		a = OpenFileWrite(buf,false);

		i = 3;
		fwrite(&i,4,1,a);
		i = 50;
		fwrite(&i,4,1,a);

		if (sdf)
		{
			i = m_iRes[0];
			fwrite(&i,4,1,a);
			i = m_iRes[1];
			fwrite(&i,4,1,a);
			i = m_iRes[2];
			fwrite(&i,4,1,a);

			ff = (float)m_fMinVal[0]/100.0f;
			fwrite(&ff,4,1,a);
			ff = (float)m_fMaxVal[0]/100.0f;
			fwrite(&ff,4,1,a);

			ff = (float)m_fMinVal[1]/100.0f;
			fwrite(&ff,4,1,a);
			ff = (float)m_fMaxVal[1]/100.0f;
			fwrite(&ff,4,1,a);

			ff = (float)m_fMinVal[2]/100.0f;
			fwrite(&ff,4,1,a);
			ff = (float)m_fMaxVal[2]/100.0f;
			fwrite(&ff,4,1,a);
		} else
		{
			i = m_iRes[0]+2;
			fwrite(&i,4,1,a);
			i = m_iRes[1]+2;
			fwrite(&i,4,1,a);
			i = m_iRes[2]+2;
			fwrite(&i,4,1,a);

			ff = (float)0;
			fwrite(&ff,4,1,a);
			ff = (float)100;
			fwrite(&ff,4,1,a);

			ff = (float)0;
			fwrite(&ff,4,1,a);
			ff = (float)100;
			fwrite(&ff,4,1,a);

			ff = (float)0;
			fwrite(&ff,4,1,a);
			ff = (float)100;
			fwrite(&ff,4,1,a);
		}

		if (sdf)
		{
			for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			{
				ff = (float)m_pBin[z];
				fwrite(&ff,4,1,a);
			}
		} else
		{
			for (z=0;z<m_iRes[2]+2;z++)
			{
				for (y=0;y<m_iRes[1]+2;y++)
				{
					for (x=0;x<m_iRes[0]+2;x++)
					{
						if ((z == 0) || (z == m_iRes[2]+1) || (y == 0) || (y == m_iRes[1]+1) || (x == 0) || (x == m_iRes[0]+1))
							ff = 0;
								else ff = (float)m_pBin[(z-1)*m_iRes[0]*m_iRes[1] + (y-1)*m_iRes[0] + x-1];
						fwrite(&ff,4,1,a);
					}
				}
			}
		}
		fclose(a);
		BTOUT;
	}


	void CorrectRadialDist(int i)
	{
		BTIN;
		int z, z2, z3;
		T d;

		if (i == 0)
		{
			for (z=0;z<m_iRes[0];z++)
			{
				d = (T)(pow(m_fMinVal[0]+(z+1.0)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),3) - pow(m_fMinVal[0]+((double)z)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),3));
	//			mprintf("  %d: %10G  %10G\n",z,d,1.0/d);
				for (z2=0;z2<m_iRes[1];z2++)
				{
					for (z3=0;z3<m_iRes[2];z3++)
					{
						m_pBin[z3*m_iResXY+z2*m_iRes[0]+z] /= d;
					}
				}
			}
		} else if (i == 1)
		{
			for (z2=0;z2<m_iRes[1];z2++)
			{
				d = (T)(pow(m_fMinVal[1]+(z2+1.0)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),3) - pow(m_fMinVal[1]+((double)z2)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),3));
	//			mprintf("  %d: %10G  %10G\n",z2,d,1.0/d);
				for (z=0;z<m_iRes[0];z++)
				{
					for (z3=0;z3<m_iRes[2];z3++)
					{
						m_pBin[z3*m_iResXY+z2*m_iRes[0]+z] /= d;
					}
				}
			}
		} else if (i == 2)
		{
			for (z3=0;z3<m_iRes[2];z3++)
			{
				d = (T)(pow(m_fMinVal[2]+(z3+1.0)/m_iRes[2]*(m_fMaxVal[2]-m_fMinVal[2]),3) - pow(m_fMinVal[2]+((double)z3)/m_iRes[2]*(m_fMaxVal[2]-m_fMinVal[2]),3));
	//			mprintf("  %d: %10G  %10G\n",z3,d,1.0/d);
				for (z=0;z<m_iRes[0];z++)
				{
					for (z2=0;z2<m_iRes[1];z2++)
					{
						m_pBin[z3*m_iResXY+z2*m_iRes[0]+z] /= d;
					}
				}
			}
		}
		BTOUT;
	}


	void CorrectAngle(int i)
	{
		BTIN;
		int z, z2, z3;
		T d;

		if (i == 0)
		{
			for (z=0;z<m_iRes[0];z++)
			{
				d = (T)(cos((m_fMinVal[0]+(double)z*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0])*Pi/180.0) - cos((m_fMinVal[0]+(double)(z+1)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0])*Pi/180.0));
				for (z2=0;z2<m_iRes[1];z2++)
				{
					for (z3=0;z3<m_iRes[2];z3++)
					{
						m_pBin[z3*m_iResXY+z2*m_iRes[0]+z] /= d;
					}
				}
			}
		} else if (i == 1)
		{
			for (z2=0;z2<m_iRes[1];z2++)
			{
				d = (T)(cos((m_fMinVal[1]+(double)z2*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1])*Pi/180.0) - cos((m_fMinVal[1]+(double)(z2+1)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1])*Pi/180.0));
				for (z=0;z<m_iRes[0];z++)
				{
					for (z3=0;z3<m_iRes[2];z3++)
					{
						m_pBin[z3*m_iResXY+z2*m_iRes[0]+z] /= d;
					}
				}
			}
		} else if (i == 2)
		{
			for (z3=0;z3<m_iRes[2];z3++)
			{
				d = (T)(cos((m_fMinVal[2]+(double)z3*(m_fMaxVal[2]-m_fMinVal[2])/m_iRes[2])*Pi/180.0) - cos((m_fMinVal[2]+(double)(z3+1)*(m_fMaxVal[2]-m_fMinVal[2])/m_iRes[2])*Pi/180.0));
				for (z=0;z<m_iRes[0];z++)
				{
					for (z2=0;z2<m_iRes[1];z2++)
					{
						m_pBin[z3*m_iResXY+z2*m_iRes[0]+z] /= d;
					}
				}
			}
		}
		BTOUT;
	}


	double NormalizeBinIntegral(T val)
	{
		BTIN;
		int z;
		T d;

		d = 0;
		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			d += m_pBin[z];
		d = val/d;
		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			m_pBin[z]*=d;
		return d;
		BTOUT;
	}


	void CalcMaxEntry()
	{
		BTIN;
		int z, z2, z3;
		
		m_fMaxEntry = (T)-99999999.0;
		m_fMinEntry = (T)99999999.0;
		for (z=0;z<m_iRes[2];z++)
		{
			for (z2=0;z2<m_iRes[1];z2++)
			{
				for (z3=0;z3<m_iRes[0];z3++)
				{
					if (m_pBin[z*m_iResXY+z2*m_iRes[0]+z3] > m_fMaxEntry)
						m_fMaxEntry = m_pBin[z*m_iResXY+z2*m_iRes[0]+z3];
					if (m_pBin[z*m_iResXY+z2*m_iRes[0]+z3] < m_fMinEntry)
						m_fMinEntry = m_pBin[z*m_iResXY+z2*m_iRes[0]+z3];
				}
			}
		}
		BTOUT;
	}


	void NormalizeBin(T mi, T ma)
	{
		BTIN;
		int z;
		T tmi, tma, d, td;

		tmi = (T)99999999.0f;
		tma = (T)0.0f;
		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
		{
			if (m_pBin[z] < tmi)
				tmi = m_pBin[z];
			if (m_pBin[z] > tma)
				tma = m_pBin[z];
		}
		if (tma-tmi < (T)1E-20f)
			tma += (T)0.00001f;
		d = ma - mi;
		td = tma - tmi;
		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			m_pBin[z] = ((m_pBin[z]-tmi)/td*d)+mi;
		BTOUT;
	}


	void Smooth(int n)
	{
		BTIN;
		int x, y, z;
		int px, py, pz, iy, iz;
		T a, b, r, ty, tz;
		T *tbuf;

		mprintf("    Smoothing 3D distribution (grade %d).",n);

		try { tbuf = new T[m_iRes[0]*m_iRes[1]*m_iRes[2]]; } catch(...) { tbuf = NULL; }
		if (tbuf == NULL) NewException((double)m_iRes[0]*m_iRes[1]*m_iRes[2]*sizeof(T),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		for (z=0;z<m_iRes[2];z++)
		{
			if ((z%(m_iRes[2]/30))==0)
				mprintf(".");
			for (y=0;y<m_iRes[1];y++)
			{
				for (x=0;x<m_iRes[0];x++)
				{
					a = 0;
					b = 0;
					for (pz=z-n;pz<=z+n;pz++)
					{
						if ((pz < 0) || (pz >= m_iRes[2]))
							continue;
						tz = (z-pz)*(z-pz);
						iz = pz*m_iResXY;
						for (py=y-n;py<=y+n;py++)
						{
							if ((py < 0) || (py >= m_iRes[1]))
								continue;
							ty = (y-py)*(y-py);
							iy = py*m_iRes[0];
							for (px=x-n;px<=x+n;px++)
							{
								if ((px < 0) || (px >= m_iRes[0]))
									continue;
								r = 1.0/((x-px)*(x-px)+ty+tz+1.0);
								a += m_pBin[iz+iy+px]*r;
								b += r;
							}
						}
					}
					tbuf[z*m_iResXY+y*m_iRes[0]+x] = a / b;
				}
			}
		}
		delete[] m_pBin;
		m_pBin = tbuf;
		mprintf("Done.\n");
		BTOUT;
	}


	void CopyFrom(C3DF *p)
	{
		BTIN;
		int z;
		m_fBinEntries = p->m_fBinEntries;
		m_fSkipEntries = p->m_fSkipEntries;
		m_iRes[0] = p->m_iRes[0];
		m_iRes[1] = p->m_iRes[1];
		m_iRes[2] = p->m_iRes[2];
		m_iResXY = p->m_iResXY;
		m_iHistogramRes = p->m_iHistogramRes;
		m_fMinVal[0] = p->m_fMinVal[0];
		m_fMaxVal[0] = p->m_fMaxVal[0];
		m_fMinVal[1] = p->m_fMinVal[1];
		m_fMaxVal[1] = p->m_fMaxVal[1];
		m_fMinVal[2] = p->m_fMinVal[2];
		m_fMaxVal[2] = p->m_fMaxVal[2];
		m_fFac[0] = p->m_fFac[0];
		m_fFac[1] = p->m_fFac[1];
		m_fFac[2] = p->m_fFac[2];
		m_fMinEntry = p->m_fMinEntry;
		m_fMaxEntry = p->m_fMaxEntry;

		try { m_pBin = new T[m_iRes[0]*m_iRes[1]*m_iRes[2]]; } catch(...) { m_pBin = NULL; }
		if (m_pBin == NULL) NewException((double)m_iRes[0]*m_iRes[1]*m_iRes[2]*sizeof(T),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	//	for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
	//		m_pBin[z] = p->m_pBin[z];
		memcpy(m_pBin,p->m_pBin,sizeof(T)*m_iRes[0]*m_iRes[1]*m_iRes[2]);

		if (m_iHistogramRes != 0)
		{
			try { m_pHistogram = new double[m_iHistogramRes]; } catch(...) { m_pHistogram = NULL; }
			if (m_pHistogram == NULL) NewException((double)m_iHistogramRes*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z=0;z<m_iHistogramRes;z++)
				m_pHistogram[z] = p->m_pHistogram[z];
		}

		if (p->m_sLabelX != NULL)
		{
			try { m_sLabelX = new char[strlen(p->m_sLabelX)+1]; } catch(...) { m_sLabelX = NULL; }
			if (m_sLabelX == NULL) NewException((double)(strlen(p->m_sLabelX)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(m_sLabelX,p->m_sLabelX);
		}
		if (p->m_sLabelY != NULL)
		{
			try { m_sLabelY = new char[strlen(p->m_sLabelY)+1]; } catch(...) { m_sLabelY = NULL; }
			if (m_sLabelY == NULL) NewException((double)(strlen(p->m_sLabelY)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			strcpy(m_sLabelY,p->m_sLabelY);
		}
		if (p->m_sLabelZ != NULL)
		{
			try { m_sLabelZ = new char[strlen(p->m_sLabelZ)+1]; } catch(...) { m_sLabelZ = NULL; }
			if (m_sLabelZ == NULL) NewException((double)(strlen(p->m_sLabelZ)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			strcpy(m_sLabelZ,p->m_sLabelZ);
		}
		BTOUT;
	}


	void WriteCube(const char *prefix, const char *name, const char *suffix, bool sdf)
	{
		BTIN;
		FILE *a;
		int x, y, z;
//		char buf[256];
		CxString buf;
		
//		strcpy(buf,prefix);
//		strcat(buf,name);
//		strcat(buf,suffix);
		buf.strcpy(prefix);
		buf.strcat(name);
		buf.strcat(suffix);

		a = OpenFileWrite(buf,true);

		mfprintf(a,"GAUSSIAN CUBE FILE; written by TRAVIS\n");
		mfprintf(a,"OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
		if (sdf)
		{
			mfprintf(a,"1 %f %f %f\n",(m_fMinVal[0]+0.5f*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0])/100.0/0.529177249,(m_fMinVal[1]+0.5f*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1])/100.0/0.529177,(m_fMinVal[2]+0.5f*(m_fMaxVal[2]-m_fMinVal[2])/m_iRes[2])/100.0/0.529177);
			mfprintf(a,"%d %f 0.000 0.000\n",m_iRes[0],(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0]/100.0/0.529177249);
			mfprintf(a,"%d 0.000 %f 0.000\n",m_iRes[1],(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1]/100.0/0.529177249);
			mfprintf(a,"%d 0.000 0.000 %f\n",m_iRes[2],(m_fMaxVal[2]-m_fMinVal[2])/m_iRes[2]/100.0/0.529177249);
		} else
		{
			mfprintf(a,"1 0.0 0.0 0.0\n");
			mfprintf(a,"%d %f 0.000 0.000\n",m_iRes[0]+2,100.0/m_iRes[0]/0.529177249);
			mfprintf(a,"%d 0.000 %f 0.000\n",m_iRes[1]+2,100.0/m_iRes[1]/0.529177249);
			mfprintf(a,"%d 0.000 0.000 %f\n",m_iRes[2]+2,100.0/m_iRes[2]/0.529177249);
		}
		mfprintf(a,"0 0 0 0 0\n");

		if (sdf)
		{
			for (x=0;x<m_iRes[0];x++)
				for (y=0;y<m_iRes[1];y++)
					for (z=0;z<m_iRes[2];z++)
					{
						mfprintf(a,"%#.8G ",m_pBin[z*m_iResXY + y*m_iRes[0] + x]);
						if (((x*m_iRes[1]*m_iRes[2] + y*m_iRes[2] + z + 1) % 6) == 0)
							mfprintf(a,"\n");
					}
		} else
		{
			for (x=0;x<m_iRes[0]+2;x++)
				for (y=0;y<m_iRes[1]+2;y++)
					for (z=0;z<m_iRes[2]+2;z++)
					{
						if ((z == 0) || (z == m_iRes[2]+1) || (y == 0) || (y == m_iRes[1]+1) || (x == 0) || (x == m_iRes[0]+1))
							mfprintf(a,"%#.8G ",0.0);
								else mfprintf(a,"%#.8G ",m_pBin[(z-1)*m_iResXY + (y-1)*m_iRes[0] + x-1]);
						if (((x*(m_iRes[1]+2)*(m_iRes[2]+2) + y*(m_iRes[2]+2) + z + 1) % 6) == 0)
							mfprintf(a,"\n");
					}
		}

		fclose(a);
		BTOUT;
	}


	void ReadCube(FILE *a, bool verbose, bool create)
	{
		BTIN;
		int ac, rx, ry, rz, z, ix, iy, iz;
		char buf[256], *p, *q;
		T su, tf;

		fgets(buf,256,a);
		fgets(buf,256,a);

		fgets(buf,256,a);
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;

		ac = atoi(buf);

		fgets(buf,256,a);
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;
		rx = atoi(buf);

		fgets(buf,256,a);
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;
		ry = atoi(buf);

		fgets(buf,256,a);
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;
		rz = atoi(buf);

		for (z=0;z<ac;z++)
			fgets(buf,256,a);

		if (feof(a))
		{
			eprintf("Unexpected end of cube file stream.\n");
			abort();
		}

		if (create)
		{
			m_iRes[0] = rx;
			m_iRes[1] = ry;
			m_iRes[2] = rz;

			m_fMinVal[0] = 0;
			m_fMaxVal[0] = g_fBoxX;
			m_fMinVal[1] = 0;
			m_fMaxVal[1] = g_fBoxY;
			m_fMinVal[2] = 0;
			m_fMaxVal[2] = g_fBoxZ;

			Create();
			mprintf("\n");
		} else
		{
			if ((rx != m_iRes[0]) || (ry != m_iRes[1]) || (rz != m_iRes[2]))
			{
				eprintf("Cube file dimension mismatch (%d-%d, %d-%d, %d-%d).\n",rx,m_iRes[0],ry,m_iRes[1],rz,m_iRes[2]);
				abort();
			}
			for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
				m_pBin[z] = 0;
		}

		if (verbose)
		{
			mprintf("    Reading cube file (resolution %d x %d x %d, %d atoms)...\n",rx,ry,rz,ac);
			mprintf(WHITE,"      [");
		}

		ix = 0;
		iy = 0;
		iz = 0;

		su = 0;
		while (!feof(a))
		{
	_read:
			fgets(buf,256,a);
			if (feof(a))
				break;
			p = buf;

	_next:
			while (*p == ' ')
				p++;
			q = p;
			while ((*p != ' ') && (*p != '\n') && (*p != 0))
				p++;
			if ((p-q) < 8)
				goto _read;

			*p = 0;

			tf = (T)atof(q);
			m_pBin[iz*m_iResXY+iy*m_iRes[0]+ix] = tf;

			if (verbose)
				su += tf;

			iz++;
			if (iz >= rz)
			{
				iz = 0;
				iy++;
				if (iy >= ry)
				{
					iy = 0;
					ix++;
					if (verbose)
						if (fmod(ix,rx/60.0) < 1.0)
							mprintf(WHITE,"#");
				}
			}

			if (ix == rx)
				break;

			p++;
			goto _next;
		}

		if (feof(a))
		{
			eprintf("Unexpected end of cube file stream.\n");
			abort();
		}

		if (verbose)
		{
			mprintf(WHITE,"]");
			mprintf(" done.\n");
			mprintf("    Sum is %f, equals %f electrons.\n",su,su*(g_fBoxX/100.0/0.529177249/rx)*(g_fBoxY/100.0/0.529177249/ry)*(g_fBoxZ/100.0/0.529177249/rz));
		}

		BTOUT;
	}


	void ReadCubeData(FILE *a, bool verbose)
	{
		BTIN;
		int ix, iy, iz;
		char buf[256], *p, *q;
		T su, tf;

		if (verbose)
		{
			mprintf("    Reading cube file (resolution %d x %d x %d)...\n",m_iRes[0],m_iRes[1],m_iRes[2]);
			mprintf(WHITE,"      [");
		}

		ix = 0;
		iy = 0;
		iz = 0;

		su = 0;
		while (!feof(a))
		{
	_read:
			fgets(buf,256,a);
			if (feof(a))
				break;
			p = buf;

	_next:
			while (*p == ' ')
				p++;
			q = p;
			while ((*p != ' ') && (*p != '\n') && (*p != 0))
				p++;
			if ((p-q) < 8)
				goto _read;

			*p = 0;

			tf = (T)atof(q);
			
			m_pBin[iz*m_iResXY+iy*m_iRes[0]+ix] = tf;

			if (verbose)
				su += tf;

			iz += g_iCubeZStride;
			if (iz >= m_iRes[2] + g_iCubeZMismatch) {
				iz = 0;
				iy += g_iCubeYStride;
				if (iy >= m_iRes[1] + g_iCubeYMismatch) {
					iy = 0;
					ix += g_iCubeXStride;
					if (verbose)
						if (fmod(ix,m_iRes[0]/60.0) < 1.0)
							mprintf(WHITE,"#");
				}
			}
			
			if (ix == m_iRes[0] + g_iCubeXMismatch)
				break;
			
	// 		iz++;
	// 		if (iz >= m_iRes[2])
	// 		{
	// 			iz = 0;
	// 			iy++;
	// 			if (iy >= m_iRes[1])
	// 			{
	// 				iy = 0;
	// 				ix++;
	// 				if (verbose)
	// 					if (fmod(ix,m_iRes[0]/60.0) < 1.0)
	// 						mprintf(WHITE,"#");
	// 			}
	// 		}
	// 
	// 		if (ix == m_iRes[0])
	// 			break;

			p++;
			goto _next;
		}
		
		if (g_iCubeXStride > 1) {
			int i;
			for (i = 0; i < m_iRes[2]; i += g_iCubeZStride) {
				int j;
				for (j = 0; j < m_iRes[1]; j += g_iCubeYStride) {
					int k;
					for (k = 0; k < m_iRes[0] - g_iCubeXStride; k += g_iCubeXStride) {
						int l;
						for (l = 1; l < g_iCubeXStride; l++) {
							m_pBin[i*m_iResXY + j*m_iRes[0] + k + l] = m_pBin[i*m_iResXY + j*m_iRes[0] + k] * (g_iCubeXStride - l) / g_iCubeXStride + m_pBin[i*m_iResXY + j*m_iRes[0] + k + g_iCubeXStride] * l / g_iCubeXStride;
						}
					}
					for (k = 1; k < g_iCubeXStride - g_iCubeXMismatch; k++) {
						m_pBin[i*m_iResXY + j*m_iRes[0] + m_iRes[0] - g_iCubeXStride + g_iCubeXMismatch + k] = m_pBin[i*m_iResXY + j*m_iRes[0] + m_iRes[0] - g_iCubeXStride + g_iCubeXMismatch] * (g_iCubeXStride - g_iCubeXMismatch - k) / (g_iCubeXStride - g_iCubeXMismatch) + m_pBin[i*m_iResXY + j*m_iRes[0]] * k / (g_iCubeXStride - g_iCubeXMismatch);
					}
				}
			}
		}
		
		if (g_iCubeYStride > 1) {
			int i;
			for (i = 0; i < m_iRes[2]; i += g_iCubeZStride) {
				int j;
				for (j = 0; j < m_iRes[0]; j++) {
					int k;
					for (k = 0; k < m_iRes[1] - g_iCubeYStride; k += g_iCubeYStride) {
						int l;
						for (l = 1; l < g_iCubeYStride; l++) {
							m_pBin[i*m_iResXY + (k + l)*m_iRes[0] + j] = m_pBin[i*m_iResXY + k*m_iRes[0] + j] * (g_iCubeYStride - l) / g_iCubeYStride + m_pBin[i*m_iResXY + (k + g_iCubeYStride)*m_iRes[0] + j] * l / g_iCubeYStride;
						}
					}
					for (k = 1; k < g_iCubeYStride - g_iCubeYMismatch; k++) {
						m_pBin[i*m_iResXY + (m_iRes[1] - g_iCubeYStride + g_iCubeYMismatch + k)*m_iRes[0] + j] = m_pBin[i*m_iResXY + (m_iRes[1] - g_iCubeYStride + g_iCubeYMismatch)*m_iRes[0] + j] * (g_iCubeYStride - g_iCubeYMismatch - k) / (g_iCubeYStride - g_iCubeYMismatch) + m_pBin[i*m_iResXY + j] * k / (g_iCubeYStride - g_iCubeYMismatch);
					}
				}
			}
		}
		
		if (g_iCubeZStride > 1) {
			int i;
			for (i = 0; i < m_iRes[1]; i++) {
				int j;
				for (j = 0; j < m_iRes[0]; j++) {
					int k;
					for (k = 0; k < m_iRes[2] - g_iCubeZStride; k += g_iCubeZStride) {
						int l;
						for (l = 1; l < g_iCubeZStride; l++) {
							m_pBin[(k + l)*m_iResXY + i*m_iRes[0] + j] = m_pBin[k*m_iResXY + i*m_iRes[0] + j] * (g_iCubeZStride - l) / g_iCubeZStride + m_pBin[(k + g_iCubeZStride)*m_iResXY + i*m_iRes[0] + j] * l / g_iCubeZStride;
						}
					}
					for (k = 1; k < g_iCubeZStride - g_iCubeZMismatch; k++) {
						m_pBin[(m_iRes[2] - g_iCubeZStride + g_iCubeZMismatch + k)*m_iResXY + i*m_iRes[0] + j] = m_pBin[(m_iRes[2] - g_iCubeZStride + g_iCubeZMismatch)*m_iResXY + i*m_iRes[0] + j] * (g_iCubeZStride - g_iCubeZMismatch - k) / (g_iCubeZStride - g_iCubeZMismatch) + m_pBin[i*m_iRes[0] + j] * k / (g_iCubeZStride - g_iCubeZMismatch);
					}
				}
			}
		}

		if (feof(a))
		{
			eprintf("Unexpected end of cube file stream.\n");
			abort();
		}

		if (verbose)
		{
			mprintf(WHITE,"]");
			mprintf(" done.\n");
	// 		mprintf("    Sum is %f, equals %f electrons.\n",su,su*(g_fBoxX/100.0/0.529177249/m_iRes[0])*(g_fBoxY/100.0/0.529177249/m_iRes[1])*(g_fBoxZ/100.0/0.529177249/m_iRes[2]));
			mprintf("    Sum is %f, equals %f electrons.\n",su,su*g_fCubeXStep*g_fCubeYStep*g_fCubeZStep);
		}

		BTOUT;
	}


	void ReadCube(CxMemFile *mf, bool verbose, bool create)
	{
		BTIN;
		int ac, rx, ry, rz, z, ix, iy, iz;
		char buf[256], *p, *q;
		T su, tf;

		mf->Seek(0);
		mf->fgets(buf,256);
		mf->fgets(buf,256);

		mf->fgets(buf,256);
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;

		ac = atoi(buf);

		mf->fgets(buf,256);
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;
		rx = atoi(buf);

		mf->fgets(buf,256);
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;
		ry = atoi(buf);

		mf->fgets(buf,256);
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;
		rz = atoi(buf);

		for (z=0;z<ac;z++)
			mf->fgets(buf,256);

		if (mf->Eof())
		{
			eprintf("Unexpected end of cube file stream.\n");
			abort();
		}

		if (create)
		{
			m_iRes[0] = rx;
			m_iRes[1] = ry;
			m_iRes[2] = rz;

			m_fMinVal[0] = 0;
			m_fMaxVal[0] = g_fBoxX;
			m_fMinVal[1] = 0;
			m_fMaxVal[1] = g_fBoxY;
			m_fMinVal[2] = 0;
			m_fMaxVal[2] = g_fBoxZ;

			Create();
			mprintf("\n");
		} else
		{
			if ((rx != m_iRes[0]) || (ry != m_iRes[1]) || (rz != m_iRes[2]))
			{
				eprintf("Cube file dimension mismatch (%d-%d, %d-%d, %d-%d).\n",rx,m_iRes[0],ry,m_iRes[1],rz,m_iRes[2]);
				abort();
			}
			for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
				m_pBin[z] = 0;
		}

		if (verbose)
		{
			mprintf("    Reading cube file (resolution %d x %d x %d, %d atoms)...\n",rx,ry,rz,ac);
			mprintf(WHITE,"      [");
		}

		ix = 0;
		iy = 0;
		iz = 0;

		su = 0;
		while (!mf->Eof())
		{
	_read:
			mf->fgets(buf,256);
			if (mf->Eof())
				break;
			p = buf;

	_next:
			while (*p == ' ')
				p++;
			q = p;
			while ((*p != ' ') && (*p != '\n') && (*p != 0))
				p++;
			if ((p-q) < 8)
				goto _read;

			*p = 0;

			tf = (T)atof(q);
			m_pBin[iz*m_iResXY+iy*m_iRes[0]+ix] = tf;

			if (verbose)
				su += tf;

			iz++;
			if (iz >= rz)
			{
				iz = 0;
				iy++;
				if (iy >= ry)
				{
					iy = 0;
					ix++;
					if (verbose)
						if (fmod(ix,rx/60.0) < 1.0)
							mprintf(WHITE,"#");
				}
			}

			if (ix == rx)
				break;

			p++;
			goto _next;
		}

		if (mf->Eof())
		{
			eprintf("Unexpected end of cube file stream.\n");
			abort();
		}

		if (verbose)
		{
			mprintf(WHITE,"]");
			mprintf(" done.\n");
			mprintf("    Sum is %f, equals %f electrons.\n",su,su*(g_fBoxX/100.0/0.529177249/rx)*(g_fBoxY/100.0/0.529177249/ry)*(g_fBoxZ/100.0/0.529177249/rz));
		}

		BTOUT;
	}


	void ReadCubeData(CxMemFile *mf, bool verbose)
	{
		BTIN;
		int ix, iy, iz;
		char buf[256], *p, *q;
		T su, tf;

		if (verbose)
		{
			mprintf("    Reading cube file (resolution %d x %d x %d)...\n",m_iRes[0],m_iRes[1],m_iRes[2]);
			mprintf(WHITE,"      [");
		}

		ix = 0;
		iy = 0;
		iz = 0;

		su = 0;
		while (!mf->Eof())
		{
	_read:
			mf->fgets(buf,256);
			if (mf->Eof())
				break;
			p = buf;

	_next:
			while (*p == ' ')
				p++;
			q = p;
			while ((*p != ' ') && (*p != '\n') && (*p != 0))
				p++;
			if ((p-q) < 8)
				goto _read;

			*p = 0;

			tf = (T)atof(q);
			
			m_pBin[iz*m_iResXY+iy*m_iRes[0]+ix] = tf;

			if (verbose)
				su += tf;

			iz += g_iCubeZStride;
			if (iz >= m_iRes[2] + g_iCubeZMismatch) {
				iz = 0;
				iy += g_iCubeYStride;
				if (iy >= m_iRes[1] + g_iCubeYMismatch) {
					iy = 0;
					ix += g_iCubeXStride;
					if (verbose)
						if (fmod(ix,m_iRes[0]/60.0) < 1.0)
							mprintf(WHITE,"#");
				}
			}
			
			if (ix == m_iRes[0] + g_iCubeXMismatch)
				break;
			
	// 		iz++;
	// 		if (iz >= m_iRes[2])
	// 		{
	// 			iz = 0;
	// 			iy++;
	// 			if (iy >= m_iRes[1])
	// 			{
	// 				iy = 0;
	// 				ix++;
	// 				if (verbose)
	// 					if (fmod(ix,m_iRes[0]/60.0) < 1.0)
	// 						mprintf(WHITE,"#");
	// 			}
	// 		}
	// 
	// 		if (ix == m_iRes[0])
	// 			break;

			p++;
			goto _next;
		}
		
		if (g_iCubeXStride > 1) {
			int i;
			for (i = 0; i < m_iRes[2]; i += g_iCubeZStride) {
				int j;
				for (j = 0; j < m_iRes[1]; j += g_iCubeYStride) {
					int k;
					for (k = 0; k < m_iRes[0] - g_iCubeXStride; k += g_iCubeXStride) {
						int l;
						for (l = 1; l < g_iCubeXStride; l++) {
							m_pBin[i*m_iResXY + j*m_iRes[0] + k + l] = m_pBin[i*m_iResXY + j*m_iRes[0] + k] * (g_iCubeXStride - l) / g_iCubeXStride + m_pBin[i*m_iResXY + j*m_iRes[0] + k + g_iCubeXStride] * l / g_iCubeXStride;
						}
					}
					for (k = 1; k < g_iCubeXStride - g_iCubeXMismatch; k++) {
						m_pBin[i*m_iResXY + j*m_iRes[0] + m_iRes[0] - g_iCubeXStride + g_iCubeXMismatch + k] = m_pBin[i*m_iResXY + j*m_iRes[0] + m_iRes[0] - g_iCubeXStride + g_iCubeXMismatch] * (g_iCubeXStride - g_iCubeXMismatch - k) / (g_iCubeXStride - g_iCubeXMismatch) + m_pBin[i*m_iResXY + j*m_iRes[0]] * k / (g_iCubeXStride - g_iCubeXMismatch);
					}
				}
			}
		}
		
		if (g_iCubeYStride > 1) {
			int i;
			for (i = 0; i < m_iRes[2]; i += g_iCubeZStride) {
				int j;
				for (j = 0; j < m_iRes[0]; j++) {
					int k;
					for (k = 0; k < m_iRes[1] - g_iCubeYStride; k += g_iCubeYStride) {
						int l;
						for (l = 1; l < g_iCubeYStride; l++) {
							m_pBin[i*m_iResXY + (k + l)*m_iRes[0] + j] = m_pBin[i*m_iResXY + k*m_iRes[0] + j] * (g_iCubeYStride - l) / g_iCubeYStride + m_pBin[i*m_iResXY + (k + g_iCubeYStride)*m_iRes[0] + j] * l / g_iCubeYStride;
						}
					}
					for (k = 1; k < g_iCubeYStride - g_iCubeYMismatch; k++) {
						m_pBin[i*m_iResXY + (m_iRes[1] - g_iCubeYStride + g_iCubeYMismatch + k)*m_iRes[0] + j] = m_pBin[i*m_iResXY + (m_iRes[1] - g_iCubeYStride + g_iCubeYMismatch)*m_iRes[0] + j] * (g_iCubeYStride - g_iCubeYMismatch - k) / (g_iCubeYStride - g_iCubeYMismatch) + m_pBin[i*m_iResXY + j] * k / (g_iCubeYStride - g_iCubeYMismatch);
					}
				}
			}
		}
		
		if (g_iCubeZStride > 1) {
			int i;
			for (i = 0; i < m_iRes[1]; i++) {
				int j;
				for (j = 0; j < m_iRes[0]; j++) {
					int k;
					for (k = 0; k < m_iRes[2] - g_iCubeZStride; k += g_iCubeZStride) {
						int l;
						for (l = 1; l < g_iCubeZStride; l++) {
							m_pBin[(k + l)*m_iResXY + i*m_iRes[0] + j] = m_pBin[k*m_iResXY + i*m_iRes[0] + j] * (g_iCubeZStride - l) / g_iCubeZStride + m_pBin[(k + g_iCubeZStride)*m_iResXY + i*m_iRes[0] + j] * l / g_iCubeZStride;
						}
					}
					for (k = 1; k < g_iCubeZStride - g_iCubeZMismatch; k++) {
						m_pBin[(m_iRes[2] - g_iCubeZStride + g_iCubeZMismatch + k)*m_iResXY + i*m_iRes[0] + j] = m_pBin[(m_iRes[2] - g_iCubeZStride + g_iCubeZMismatch)*m_iResXY + i*m_iRes[0] + j] * (g_iCubeZStride - g_iCubeZMismatch - k) / (g_iCubeZStride - g_iCubeZMismatch) + m_pBin[i*m_iRes[0] + j] * k / (g_iCubeZStride - g_iCubeZMismatch);
					}
				}
			}
		}

		if (mf->Eof())
		{
			eprintf("Unexpected end of cube file stream.\n");
			abort();
		}

		if (verbose)
		{
			mprintf(WHITE,"]");
			mprintf(" done.\n");
	// 		mprintf("    Sum is %f, equals %f electrons.\n",su,su*(g_fBoxX/100.0/0.529177249/m_iRes[0])*(g_fBoxY/100.0/0.529177249/m_iRes[1])*(g_fBoxZ/100.0/0.529177249/m_iRes[2]));
			mprintf("    Sum is %f, equals %f electrons.\n",su,su*g_fCubeXStep*g_fCubeYStep*g_fCubeZStep);
		}

		BTOUT;
	}


	void MultiplyBin(T d)
	{
		int z;

		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			m_pBin[z] *= d;
	}


	void Invert()
	{
		T f;
		int x, y, z, i;

		f = 0;
		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			if (f < m_pBin[z])
				f = m_pBin[z];
		for (z=0;z<m_iRes[2];z++)
			for (y=0;y<m_iRes[1];y++)
				for (x=0;x<m_iRes[0];x++)
				{
					i = z*m_iResXY+y*m_iRes[0]+x;
					m_pBin[i] = f - m_pBin[i];
				}
	}


	void ClipPlane(int dir, double val)
	{
		int x, y, z, t;

		switch(dir)
		{
			case 0: //x
				t = (int)floor(((val-m_fMinVal[0])/(m_fMaxVal[0]-m_fMinVal[0]))*m_iRes[0]+0.5);
				for (x=t;x<m_iRes[0];x++)
					for (y=0;y<m_iRes[1];y++)
						for (z=0;z<m_iRes[2];z++)
							m_pBin[z*m_iResXY + y*m_iRes[0] + x] = 0;
				break;

			case 1: //y
				t = (int)floor(((val-m_fMinVal[1])/(m_fMaxVal[1]-m_fMinVal[1]))*m_iRes[1]+0.5);
				for (x=0;x<m_iRes[0];x++)
					for (y=t;y<m_iRes[1];y++)
						for (z=0;z<m_iRes[2];z++)
							m_pBin[z*m_iResXY + y*m_iRes[0] + x] = 0;
				break;

			case 2: //z
				t = (int)floor(((val-m_fMinVal[1])/(m_fMaxVal[2]-m_fMinVal[2]))*m_iRes[2]+0.5);
				for (x=0;x<m_iRes[0];x++)
					for (y=0;y<m_iRes[1];y++)
						for (z=t;z<m_iRes[2];z++)
							m_pBin[z*m_iResXY + y*m_iRes[0] + x] = 0;
				break;
		}
	}


	void CalcHistogram()
	{
		int z, z2, t;

		CalcMaxEntry();

		try { m_pHistogram = new double[m_iHistogramRes]; } catch(...) { m_pHistogram = NULL; }
		if (m_pHistogram == NULL) NewException((double)m_iHistogramRes*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<m_iHistogramRes;z++)
			m_pHistogram[z] = 0;
		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
		{
			t = (int)floor((m_pBin[z]-m_fMinEntry)/(m_fMaxEntry-m_fMinEntry)*m_iHistogramRes);
			for (z2=0;z2<t;z2++)
				m_pHistogram[z2]++;
		}
		for (z=0;z<m_iHistogramRes;z++)
			m_pHistogram[z] /= ((double)m_iRes[0])*m_iRes[1]*m_iRes[2];
	}


	void WriteHistogram(const char *prefix, const char *name, const char *suffix)
	{
		BTIN;
		FILE *a;
		int z;
//		char buf[256];
		CxString buf;
		
//		strcpy(buf,prefix);
//		strcat(buf,name);
//		strcat(buf,suffix);
		buf.strcpy(prefix);
		buf.strcat(name);
		buf.strcat(suffix);
		
		a = OpenFileWrite(buf,true);

		for (z=0;z<m_iHistogramRes;z++)
			mfprintf(a," %#.10G;  %#.10G\n",(z+0.5)*(m_fMaxEntry-m_fMinEntry)/m_iHistogramRes,m_pHistogram[z]);
		
		fclose(a);
		BTOUT;
	}


	void CreateSlice(int axis, int v1, int v2, C2DF *cdf)
	{
		int i1, i2, i3;

		if (axis == 0)
		{
			cdf->m_iRes[0] = m_iRes[1];
			cdf->m_fMinVal[0] = m_fMinVal[1];
			cdf->m_fMaxVal[0] = m_fMaxVal[1];
			cdf->SetLabelX(m_sLabelY);

			cdf->m_iRes[1] = m_iRes[2];
			cdf->m_fMinVal[1] = m_fMinVal[2];
			cdf->m_fMaxVal[1] = m_fMaxVal[2];
			cdf->SetLabelY(m_sLabelZ);

			cdf->SetLabelZ("Occurence");

			cdf->Create();
			cdf->m_iSmoothGrade = 0;

			for (i3=0;i3<m_iRes[2];i3++)
			{
				for (i2=0;i2<m_iRes[1];i2++)
				{
					for (i1=v1;i1<=v2;i1++)
					{
						cdf->m_pBin[i3*m_iRes[1]+i2] += m_pBin[i3*m_iResXY+i2*m_iRes[0]+i1];
					}
				}
			}
			cdf->CalcMaxEntry();
		} else if (axis == 1)
		{
			cdf->m_iRes[0] = m_iRes[0];
			cdf->m_fMinVal[0] = m_fMinVal[0];
			cdf->m_fMaxVal[0] = m_fMaxVal[0];
			cdf->SetLabelX(m_sLabelX);

			cdf->m_iRes[1] = m_iRes[2];
			cdf->m_fMinVal[1] = m_fMinVal[2];
			cdf->m_fMaxVal[1] = m_fMaxVal[2];
			cdf->SetLabelY(m_sLabelZ);

			cdf->SetLabelZ("Occurence");

			cdf->Create();
			cdf->m_iSmoothGrade = 0;

			for (i3=0;i3<m_iRes[2];i3++)
			{
				for (i1=0;i1<m_iRes[0];i1++)
				{
					for (i2=v1;i2<=v2;i2++)
					{
						cdf->m_pBin[i3*m_iRes[0]+i1] += m_pBin[i3*m_iResXY+i2*m_iRes[0]+i1];
					}
				}
			}
			cdf->CalcMaxEntry();
		} else if (axis == 2)
		{
			cdf->m_iRes[0] = m_iRes[0];
			cdf->m_fMinVal[0] = m_fMinVal[0];
			cdf->m_fMaxVal[0] = m_fMaxVal[0];
			cdf->SetLabelX(m_sLabelX);

			cdf->m_iRes[1] = m_iRes[1];
			cdf->m_fMinVal[1] = m_fMinVal[1];
			cdf->m_fMaxVal[1] = m_fMaxVal[1];
			cdf->SetLabelY(m_sLabelY);

			cdf->SetLabelZ("Occurence");

			cdf->Create();
			cdf->m_iSmoothGrade = 0;

			for (i1=0;i1<m_iRes[0];i1++)
			{
				for (i2=0;i2<m_iRes[1];i2++)
				{
					for (i3=v1;i3<=v2;i3++)
					{
						cdf->m_pBin[i2*m_iRes[0]+i1] += m_pBin[i3*m_iResXY+i2*m_iRes[0]+i1];
					}
				}
			}
			cdf->CalcMaxEntry();
		}
	}


	void Clear()
	{
		int z;

		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			m_pBin[z] = 0;
	}

};

#endif
