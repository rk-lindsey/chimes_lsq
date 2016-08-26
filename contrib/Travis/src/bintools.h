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

#ifndef BINTOOLS_H
#define BINTOOLS_H

#include <math.h>
#include <stdio.h>
#include "xvector3.h"
#include "xdoublearray.h"
#include "backtrace.h"

/*class CS6DF : public CxObject
{
public:
	CS6DF();
	~CS6DF();
	void AddToBin(const CxVector3 &vec, const CxVector3 &val);
	double NormalizeBin(double mi, double ma);
	double PPMBin();
	void Write(const char *prefix, const char *name, const char *suffix);
	void Create();
	void WriteHistogram(const char *prefix, const char *name, const char *suffix);
	void CalcHistogram();

//	double m_fRadius;
//	int m_iResolution, m_iResSqr, m_iResTri;

	int	m_iLevelCount;
	double m_fLevelThres[32];

	double *m_pHistogram[32];
//	double *m_pDHistogram[32];
	double *m_pBin[32];
//	double *m_pDBin[3][32];
	double m_fMaxP[32];
//	double m_fMaxDP[32];

	double m_fBinEntries[32];
	double m_fGesBinEntries;
};

class CMSDF : public CxObject
{
public:
	CMSDF() { };
	~CMSDF() { };
	void AddToBinMean(const CxVector3 &vec, const CxVector3 &val);
	void AddToBinMax(const CxVector3 &vec, const CxVector3 &val);
	double NormalizeBin(double mi, double ma);
	void Write(const char *prefix, const char *name, const char *suffix);
	void Create();
	void CalcAvg();
	void WriteHistogram(const char *prefix, const char *name, const char *suffix);
	void CalcHistogram();
	void NormHistogramIntegral();

//	double *m_pDBin[3];
	double m_fMaxAvg;
	double *m_pHistogram;
	
//	double m_fRadius;
//	int m_iResolution, m_iResSqr, m_iResTri;
	double *m_pBin;
	double *m_pRefBin;
	double m_fBinEntries;
};

class CMRDF
{
public:
	CMRDF() { };
	~CMRDF() { };
	void AddToBin(double d, double vec[3]);
	double NormalizeBin(double mi, double ma);
	void Write(char *s);
	void Create();

	int	m_iLevelCount;
	double m_fLevelMin, m_fLevelMax;
	double *m_pDBin[3];

//	double m_fRadius;
//	int m_iResolution;
	double *m_pBin;
	double *m_pRefBin;
	double m_fBinEntries;
};
*/


class CAF : public CxObject
{
public:
	void LinReg(int i1, int i2, double *a0, double *a1, double *r);
	void CalcDeriv(double f);
	CAF();
	~CAF();
	void AddToBin(double x, double y);
	void AddToBin_Index(int i, double y);
	void Write(const char *prefix, const char *name, const char *suffix);
	void Create();
	void BuildAverage();
	
	double m_fMinVal, m_fMaxVal;
	int m_iResolution;
	double *m_pBin;
	double *m_pDBin;
	double m_fBinEntries;
	CxDoubleArray m_faEntries;
};


class CNDF : public CxObject
{
public:
	CNDF() { };
	~CNDF() { };

//	void MultiplyBin(double m);
//	void AddToBin(double x, const CxVector3 &vec);
/*	void AddToBin(double x, int y);
	void AddToBin(double x, double y, double val);
	double NormalizeBin(double mi, double ma);
	double PPMBin();
	void Write(const char *prefix, const char *name, const char *suffix);
	void WriteCSV(const char *prefix, const char *name, const char *suffix);
	void WriteXProjection(const char *prefix, const char *name, const char *suffix);
	void WriteYProjection(const char *prefix, const char *name, const char *suffix);
	void WriteMathematica(const char *prefix, const char *name, const char *suffix);
	void WriteLogMathematica(const char *prefix, const char *name, const char *suffix);
	void WriteVHDFLogMathematica(const char *prefix, const char *name, const char *suffix);
	void Create();
	void AngleCorrect();
	void AngleCorrectX();
	void CorrectRadialDistX();*/

	void AddToBin(double *d);

	int m_iDimensions;
	double *m_fMinVal;
	double *m_fMaxVal;

	double *m_fFac;

	int *m_iRes;

	double *m_pBin;
	double m_fBinEntries;
};


#endif
