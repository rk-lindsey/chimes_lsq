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


#ifndef DOMAIN_H
#define DOMAIN_H


#include "xobject.h"
#include "xobarray.h"
#include "xintarray.h"
#include "timestep.h"


class CDomain : public CxObject
{
public:
	CDomain();
	~CDomain();
	void Reset();
	void Assimilate(CDomain *dom);

	bool m_bActive;
	CxIntArray m_iaCells;
	CxIntArray m_iaNeighbors;

	double m_fSurfaceArea;
	double m_fVolume;
	double m_fAVRatio;
	int m_iFaces;
};


class CDomainAnalysis : public CxObject
{
public:
	CDomainAnalysis();
	~CDomainAnalysis();
	void Parse(int i);
	void ProcessStep(CTimeStep *ts);
	void Finish(int i);
	void REC_FuseDomain(int basedom, int dom, CxIntArray *stack);

	bool m_bWriteHistograms;
	CDF *m_pHistoSurface;
	CDF *m_pHistoVolume;
	CDF *m_pHistoAVRatio;
	CDF *m_pHistoFaceCount;
	CDF *m_pHistoCellCount;
	CDF *m_pHistoDomainCount;

	CxObArray m_oaBasePopulation; // Both contain one CAtomGroup object per CMolecule
	CxObArray m_oaDomainSet;

	CxIntArray m_iaBasePopulation;
	CxIntArray m_iaDomainSet;

	CxIntArray m_iaBaseInDomain;

	CxObArray m_oaDomains;

	CxDoubleArray m_faSurfaceAv;
	CxDoubleArray m_faSurfaceMin;
	CxDoubleArray m_faSurfaceMax;
	CxDoubleArray m_faSurfaceSD;

	CxDoubleArray m_faVolumeAv;
	CxDoubleArray m_faVolumeMin;
	CxDoubleArray m_faVolumeMax;
	CxDoubleArray m_faVolumeSD;

	CxDoubleArray m_faAVRatioAv;
	CxDoubleArray m_faAVRatioMin;
	CxDoubleArray m_faAVRatioMax;
	CxDoubleArray m_faAVRatioSD;

	CxDoubleArray m_faFaceCountAv;
	CxDoubleArray m_faFaceCountMin;
	CxDoubleArray m_faFaceCountMax;
	CxDoubleArray m_faFaceCountSD;

	CxDoubleArray m_faCellCountAv;
	CxDoubleArray m_faCellCountMin;
	CxDoubleArray m_faCellCountMax;
	CxDoubleArray m_faCellCountSD;

	CxIntArray m_iaDomainCount;
	CxIntArray m_iaStep;
};


class CDomainEngine : public CxObject
{
public:
	CDomainEngine();
	~CDomainEngine();
	void Parse();
	void ProcessStep(CTimeStep *ts);
	void Finish();

	CxObArray m_oaAnalyses;
};


#endif

