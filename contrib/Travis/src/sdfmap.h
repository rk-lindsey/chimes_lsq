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

#ifndef SDFMAP_H
#define SDFMAP_H


#include "xobject.h"
#include "xobarray.h"
#include "3df.h"
#include "timestep.h"


class CSDFMap : public CxObject
{
public:
	CSDFMap();
	~CSDFMap();

	void Parse();
	void Create();
	void Process(int rm, CTimeStep *ts);
	void BuildName();
	void Finish();

	int m_iQuantity;
	bool m_bModeMin;
	bool m_bModeMax;
	bool m_bModeAvg;

	CxObArray m_oaAtomGroups;

	int m_iResolution;
	float m_fRadius;
	bool m_bSphere;
	int m_iSphereRadius;

	C3DF<double> *m_pValueMin;
	C3DF<double> *m_pValueMax;
	C3DF<double> *m_pValueAvg;
	C3DF<double> *m_pCount;
	char *m_sName;
	char *m_sSubName;


	inline void AddValue(double x, double y, double z, double v)
	{
		if (m_bModeAvg)
			m_pValueAvg->AddToBin_Single(x,y,z,v);

		if (m_bModeMax)
			if (v > m_pValueMax->GetBinValue_Single(x,y,z))
				m_pValueMax->SetBinValue_Single(x,y,z,v);

		if (m_bModeMin)
		{
			if ((m_pCount->GetBinValue_Single(x,y,z) == 0) || (v < m_pValueMin->GetBinValue_Single(x,y,z)))
				m_pValueMin->SetBinValue_Single(x,y,z,v);
		}

		m_pCount->AddToBin_Single(x,y,z,1.0);
	}

};


#endif


