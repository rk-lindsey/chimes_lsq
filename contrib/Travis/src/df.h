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

#ifndef DF_H
#define DF_H

#include <math.h>
#include <stdio.h>
#include "xvector3.h"
#include "xdoublearray.h"
#include "backtrace.h"
#include "bintree.h"
#include "grace.h"
#include "lmwrapper.h"

class C2DF;
template<typename T> class C3DF;
class CObservation;


class CDF : public CxObject
{
public:
	void AddFrom(CDF *t);
	int m_iMultiCount;
	void CalcMeanSD();
	void SetAdditionalDatasetLabel(int z, const char *s);
	void SetLabelX(const char *s);
	void SetLabelY(const char *s);
	bool m_bLeft;
	void Fit_PolyExp(int degree, int maxcall);
	void Fit_ExpSpectrum(int res, double mi, double ma, const char *name, int dpoints, int maxcall, bool lindata, double zeroweight, bool evolve);
	void CopyFrom(CDF *p);
	void WriteHistogram(const char *prefix, const char *name, const char *suffix);
	void CalcHistogram();
//	bool m_bRDF;
	void Mirror(float plane);
	void CalcMinMax();
	void ScaleXRange(double fac);
	void CreateCombinedPlot(bool rdf);
	void REC_BinTreeMultiplyBin(double xmin, double xmax, int depth, CBinTree *p, double fac);
	void BinTree_MultiplyBin(double f);
	void REC_BinTreeRadialDist(double xmin, double xmax, int depth, CBinTree *p);
	void BinTree_RadialDist();
	void WriteAdapted(const char *prefix, const char *name, const char *suffix, int mindepth, int maxdepth, double thres, bool rdf);
	void REC_SaveTree(FILE *a, double xmin, double xmax, int depth, int mindepth, int maxdepth, double thres, CBinTree *p, bool rdf);
	void REC_FuseTree(CBinTree *p);
	unsigned long REC_FillBinTree(int pos, int depth, CBinTree *p);
	void PrepareAdapt();
	CDF();
	~CDF();
//	void AddToBin(const CxVector3 &vec);
	void AddToBin(double d);
	void AddToBin(double d, double v);
	void AddToBin_Multi(int i, double d);
	void AddToBin_Multi_Int(int i, int n, double f);
	void AddToBin_Int(int i);
	void AddToBin_Int(int i, double j);
	void AddToBin_Count(int i, int count);
	double NormalizeBin(double mi, double ma);
	void Write(const char *prefix, const char *name, const char *suffix, bool integral);
	void WriteMultiAgr(const char *prefix, const char *name, const char *suffix, const char *title, bool rdf);
	void WriteMultiAgr_Cumulative(const char *prefix, const char *name, const char *suffix, const char *title, bool rdf);
	void WriteMulti(const char *prefix, const char *name, const char *suffix);
	void WriteMulti_Cumulative(const char *prefix, const char *name, const char *suffix);
	void Write_Int(const char *prefix, const char *name, const char *suffix);
	void WriteHenry(const char *prefix, const char *name, const char *suffix);
//	void WriteLog(const char *prefix, const char *name, const char *suffix);
	void Create();
	void CreateMulti(int n);
	void SetLabelMulti(int n, const char *s);
//	void CalcHistogram();
	double NormBinIntegral();
	void NormBinIntegral(double val);
	void NormBinSum(double val);
	void MultiplyBin(double f);
	void SubtractBin(double f);
	void MultiplyIntegral(double f);
//	void WriteHistogram(const char *prefix, const char *name, const char *suffix);
	void AngleCorrect();
	void ZeroBin();

	void Integrate(bool correctradial, double fac);
	void CorrectRadialDist();
	void CorrectRadialDistLong();
	void CorrectLiRadialDist();

	void WriteAgr(const char *prefix, const char *name, const char *suffix, const char *title, bool rdf);

	double GetPercentageRange(double perc);
	
	double m_fMaxP;
	double m_fFac;
//	double m_fMaxDP;
//	double *m_pHistogram;
//	double *m_pDHistogram;
	
//	double *m_pDBin[3];

	double m_fSum;
	double m_fSqSum;
	double m_fMean;
	double m_fSD;
	double m_fMinInput;
	double m_fMaxInput;

	double m_fMinVal, m_fMaxVal;
	double m_fMinEntry, m_fMaxEntry;
	int m_iResolution;
	double *m_pBin;
	double **m_pMultiBin;
	double *m_pIntegral;
	double m_fBinEntries;
	double m_fSkipEntries;
	CBinTree *m_pBinTree;


	bool m_bSaveDist;
	bool m_bCombinedPlot;

	CGrace *m_pCombinedPlot;

	int m_iHistogramRes;
	double *m_pHistogram;

	CLMWrapper *m_pLMWrapper;
	double *m_pCorrCoeff;
	double *m_pFitIntegral;
	double **m_pParameters;

	double **m_pAdditionalSets;
	int m_iAdditionalSets;
	char **m_pAdditionalSetLabels;

	CxObArray m_oaTimeDiffBuf;
	CDF *m_pTimeDiff;
	CDF *m_pTimeDiffAbs;
	CDF *m_pTimeDiffSqr;
	C2DF *m_p3DTimeDiff;
	C2DF *m_p3DTimeDiffAbs;
	C2DF *m_p3DTimeDiffSqr;
	C2DF *m_p3DTimeDiffT;
	C2DF **m_pTimeDiffDistPairs;
	C3DF<double> *m_pTimeDiffDist3DF;

	char *m_sLabelX;
	char *m_sLabelY;
	char **m_sLabelMulti;

	inline void AddToBin_Fast(float d)
	{
		double p;
		int ip;

		if (d > m_fMaxVal)
			return;

		p = d*m_fFac - 0.5;
		ip = (int)floor(p);
		if (ip < 0)
		{
			ip = 0;
			p = 0;
		} else if (ip > m_iResolution-2)
		{
			ip = m_iResolution-2;
			p = 1.0;
		} else
			p -= ip;

		m_pBin[ip    ] += (1.0-p);
		m_pBin[ip + 1] +=      p ;
	}

	inline void AddToBinInt_Fast(int i)
	{
//		if (i >= m_iResolution)
//			abort();

		m_pBin[i] += 1.0;
	}

	inline void AddToBinInt_Fast(int i, double d)
	{
//		if (i >= m_iResolution)
//			abort();

		m_pBin[i] += d;
	}

	inline void AddToBin_Multi_Int_Fast(int i, int n, double f)
	{
		m_pBin[n] += f;
		m_pMultiBin[i][n] += f;
	}

//	double m_fIntegralNorm;
};

#endif

