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


#include "statistics.h"
#include "xstring.h"


CStatistics::CStatistics()
{
	m_iSizeX = 0;
	m_iSizeY = 0;
	m_pMin = NULL;
	m_pMax = NULL;
	m_pAvg = NULL;
	m_pCount = NULL;
	m_waXValues.SetName("CStatistics::m_waXValues");
	m_waYValues.SetName("CStatistics::m_waYValues");
}


CStatistics::~CStatistics()
{
	if (m_pMin != NULL)
	{
		delete[] m_pMin;
		m_pMin = NULL;
	}
	if (m_pMax != NULL)
	{
		delete[] m_pMax;
		m_pMax = NULL;
	}
	if (m_pAvg != NULL)
	{
		delete[] m_pAvg;
		m_pAvg = NULL;
	}
	if (m_pCount != NULL)
	{
		delete[] m_pCount;
		m_pCount = NULL;
	}
}


void CStatistics::Init(int x, int y)
{
	int z;

	m_iSizeX = x;
	m_iSizeY = y;

	try { m_pMin = new double[x*y]; } catch(...) { m_pMin = NULL; }
	if (m_pMin == NULL) NewException((double)x*y*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_pMax = new double[x*y]; } catch(...) { m_pMax = NULL; }
	if (m_pMax == NULL) NewException((double)x*y*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_pAvg = new double[x*y]; } catch(...) { m_pAvg = NULL; }
	if (m_pAvg == NULL) NewException((double)x*y*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_pCount = new double[x*y]; } catch(...) { m_pCount = NULL; }
	if (m_pCount == NULL) NewException((double)x*y*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<x*y;z++)
	{
		m_pMin[z] = 9999999999.0;
		m_pMax[z] = 0.0;
		m_pAvg[z] = 0.0;
		m_pCount[z] = 0.0;
	}
}


void CStatistics::Evaluate()
{
	int z;

	for (z=0;z<m_iSizeX*m_iSizeY;z++)
		if (m_pCount[z] != 0)
			m_pAvg[z] /= m_pCount[z];
}


void CStatistics::Write(const char *pre, const char *s, const char *post)
{
	BTIN;
	FILE *a;
	int x, y;
//	char buf[256];
	CxString buf;
	
//	buf[0] = 0;
//	strcpy(buf,pre);
//	strcat(buf,s);
//	strcat(buf,post);

	buf.strcpy(pre);
	buf.strcat(s);
	buf.strcat(post);

	a = OpenFileWrite(buf,true);
	for (x=0;x<m_iSizeX;x++)
	{
		for (y=0;y<m_iSizeY;y++)
		{
			fprintf(a,"%8d; %8d; ",m_waXValues[x]+1,m_waYValues[y]+1);
			if (m_pCount[y*m_iSizeX+x] == 0)
				fprintf(a,"     -    ;      -    ;      -    \n");
			else
				fprintf(a,"%10.3f; %10.3f; %10.3f\n",m_pMin[y*m_iSizeX+x],m_pMax[y*m_iSizeX+x],m_pAvg[y*m_iSizeX+x]);
		}
	}
	fclose(a);
	BTOUT;
}

