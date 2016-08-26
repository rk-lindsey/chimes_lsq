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

#ifndef STATISTICS_H
#define STATISTICS_H

#include "xobject.h"
#include "tools.h"
#include "backtrace.h"
#include "xwordarray.h"


class CStatistics : public CxObject  
{
public:
	CStatistics();
	virtual ~CStatistics();

	void Init(int x, int y);
	void Evaluate();
	void Write(const char *pre, const char *s, const char *post);

	void AddValue(int x, int y, double v)
	{
		int z;
		z = y*m_iSizeX + x;
		m_pCount[z]++;
		m_pAvg[z] += v;
		if (v < m_pMin[z])
			m_pMin[z] = v;
		if (v > m_pMax[z])
			m_pMax[z] = v;
	}

	double *m_pMin;
	double *m_pMax;
	double *m_pAvg;
	double *m_pCount;
	CxWordArray m_waXValues;
	CxWordArray m_waYValues;

	int m_iSizeX, m_iSizeY;

};

#endif

