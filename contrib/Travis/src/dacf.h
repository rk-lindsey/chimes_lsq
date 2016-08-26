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

#ifndef DACF_H
#define DACF_H

#include "xobject.h"
#include "xobarray.h"
#include "xfloatarray.h"
#include "nbsearch.h"
#include "moltools.h"
#include "nbexchange.h"


class CDACFSub : public CxObject
{
public:
	void Create(CConditionGroup *c);
	void Parse();
	CDACFSub();
	~CDACFSub();

	bool m_bNewMode;
	bool m_bBorderMode;
	bool m_bDistTrace;
	bool m_bIntermittend;
	float m_fIntGap;
	bool m_bIntTravisStyle;
	bool m_bCorrectEq;

	double m_fEqCounter;

	int m_iRefMol;
	int m_iShowMol;

	char *m_sName;
	void BuildName(const char *n);

	CConditionGroup *m_pCondition;
	CxObArray *m_oaAggregates;
	CxObArray *m_oaNbExPairs;

	C2DF *m_pDLDisp;
	CDF *m_pDACF;
	CDF *m_pDLDF;
	CDF *m_pNDF;
	CDF *m_pDDisp;
	CAF *m_pPairMSD;

	CDF **m_pNbExDF;
	CConditionGroup **m_pNbExConditions;

	CxIntArray *m_piaIntervals;
};


class CDACF : public CxObject
{
public:
	void CreateSubDACFStack(char *s);
	void CreateGridFit2DF(C2DF *df, int degree, bool intermittend);
//	void CreateGridFitDF(CDF *df, int degree, bool intermittend);
	void CalcGridFitParms();
	void UpdateNbEx(int rm, CDACFSub *dacfsub);
	void UpdateDACFSub(int rm, CTimeStep *t, CDACFSub *dacfsub);
	void FinishDACFSub(CTimeStep *t, CDACFSub *dacfsub);
	bool m_bRemoveMaxVel;
	float m_fMaxVel;
	CDACF();
	~CDACF();
	void Parse();
	CxVector3 CalcCenter(CxVec3Array *v, int i1, int i2);

	CxObArray m_oaSubDACFs;
	CConditionGroup *m_pCondition;
	int m_iFirstMol;
	int m_iSecondMol;
	char *m_sName;

	CAtomGroup *m_pCenterAtoms1;
	CAtomGroup *m_pCenterAtoms2;
	CxFloatArray m_faWeight1;
	CxFloatArray m_faWeight2;
	float m_fWeightMol1;

	int m_iLifetimeRes;
	int m_iDACFRes;
	int m_iDisplacementRes;
	float m_fLargestDisplacement;
	float m_fLargestLifetime;
	
	bool m_bDACFGrid;
	int m_iGridMode;
	bool m_bGridCon;
	bool m_bGridInt;
	bool m_bGridIntTravisStyle;
	float m_fGridIntGap;

	CNbExchange *m_pNbExchange;
	bool m_bFitDACF;
	int m_iFitDegreeMin;
	int m_iFitDegreeMax;

	double *m_pFitRMin;
	double *m_pFitRAvg;
	double *m_pFitRMax;

	int m_iGridX;
	int m_iGridY;

	float m_fGridXMin;
	float m_fGridXMax;
	float m_fGridYMin;
	float m_fGridYMax;

	bool m_bLifetimeSpectrum;
	CxObArray m_oaLTSpectra;
};

#endif
