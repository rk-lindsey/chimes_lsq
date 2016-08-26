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

#ifndef NBSEARCH_H
#define NBSEARCH_H

#include "xobject.h"
#include "xobarray.h"
#include "xwordarray.h"
#include "xintarray.h"
#include "xfloatarray.h"
#include "backtrace.h"
#include "df.h"

class CNbSearch;
class CRDF;
class CADF;
class CSingleMolecule;
class CTimeStep;

//#include "moltools.h"


class CNbPair : public CxObject
{
public:
	void ReScan(CNbSearch *parent);
	void CopyFrom(CNbPair *p, CNbSearch *parent);
	void Reset();
	void Create(CNbSearch *parent);
	CNbPair();
	~CNbPair();
	void Check(CNbSearch *parent, CTimeStep *ts, CSingleMolecule *rm, CSingleMolecule *sm);
	void PreCheck(CNbSearch *parent, CTimeStep *ts, CSingleMolecule *rm, CSingleMolecule *sm);

	int m_iNbPosition;
	bool m_bAnyDistPassed;
	bool m_bAnyAnglePassed;
	bool *m_bDistPassed;
	bool *m_bAnglePassed;
	float *m_fDistances;
	float *m_fAngles;
	float m_fMinDist;
	int m_iMinDistIndex;

	CxIntArray m_laDistAtomList;
	CxIntArray m_laAngleAtomList;
};


class CNbSort : public CxObject
{
public:
	CNbSort() { }
	~CNbSort() { }
	float m_fMinDist;
	bool m_bAnyAnglePassed;
	int m_iOM;
};


class CExtendedCondition : public CxObject
{
public:
	void CopyFrom(CExtendedCondition *ec);
	bool Evaluate();
	double m_fX;
	double m_fY;
	double m_fZ;
	bool m_bLarger;
	double m_fA[3];
	double m_fD[3];
	CExtendedCondition();
	~CExtendedCondition();
};


class CNbSearch : public CxObject  
{
public:
	bool CheckExtended(double dist, double angle);
	void ReScan(CSingleMolecule *rm);
	void CopyResults(CNbSearch *p);
	void CopyFrom(CNbSearch *nb);
	void SortNeighbors();
	void Reset();
	void PrintSingle(int om);
	void PrintTable();
	void PrintTable(FILE *a);
	void Parse_OnlyValues();
	void ParseGrid(int rm, int om, int gridmode);
	void Parse(int rm, int om, bool nbana);
	void Create(int obs);
	CNbSearch();
	~CNbSearch();
//	void ScanSingleOM(CSingleMolecule *rm, int om, CTimeStep *t, bool markpassedatoms, bool fold);
	void MarkPassedAtoms(int om);
	void ScanAllOM(CSingleMolecule *rm, CTimeStep *t);
	void PreScanAllOM(CSingleMolecule *rm, CTimeStep *t);

	int m_iRefMol;
	int m_iObsMol;
	bool m_bInactive;
	double m_fCombinationsTotal;
	double m_fCombinationsPassed;
	double m_fMoleculesTotal;
	double m_fMoleculesPassed;
	int m_iAngles;
	int m_iDistances;
	int m_iCombinationsEnabled;
	int m_iNbCountMin;
	int m_iNbCountMax;
	int m_iNumber;

	bool *m_bPassed;
	bool *m_bCombinationMatrix;
	double *m_fCombinationsFound;
	int *m_iCombPassCount;
	int *m_iPassCounter;
//	bool *m_iPassed;
	CRDF *m_pRDF;
	CADF *m_pADF;
	CNbSort *m_pNbSort;

	float *m_pDistances;
	float *m_pAngles;

	CxObArray m_oaNbPairs;

	bool m_bExtendedMode;
	CxObArray m_oaExtendedConditions;
};


class CNbAnalysis : public CxObject  
{
public:
	float m_fBinEntries;
	void AnalyzeStep();
	CNbAnalysis();
	~CNbAnalysis();
	void Parse();
	void BuildName();

	double *m_pDist;
	double *m_pDistMin;
	double *m_pDistMax;
	double *m_pDistAvg;
	double *m_pDistCount;
	char *m_sName;
	CNbSearch *m_pNbSearch;
	CxObArray m_oaDF;
	CxObArray m_oaNPF;
	CDF *m_pNPFCount;
	int m_iResolution;
	float m_fMinDist;
	float m_fMaxDist;
	int m_iNbCount;
//	int m_iMinNbCount;
//	int m_iMaxNbCount;
	int m_iShowMol;
};


class CNbSet : public CxObject
{
public:
	void ResetAlwaysTrue();
	void AddMolecule(int moltype, int mol);
	void Dump();
	void Reset();
	void Scan(CSingleMolecule *rm, CTimeStep *t);
//	CxObArray m_oaMolecules;
	CxObArray m_oaConditionGroups;
	void Parse(int rm);
	CNbSet();
	~CNbSet();
};

#endif
