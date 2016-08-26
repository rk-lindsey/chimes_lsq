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

#ifndef _2DF_H
#define _2DF_H

#include <math.h>
#include <stdio.h>
#include "xvector3.h"
#include "xdoublearray.h"
#include "backtrace.h"
#include "df.h"


class CMathematicaCircle : public CxObject
{
public:
	CMathematicaCircle() { }
	~CMathematicaCircle() { }
	double m_fPosX, m_fPosY, m_fRadius;
	double m_fColorR, m_fColorG, m_fColorB;
};


class C2DF : public CxObject
{
public:
	double CalcCorrelationFactor();
	void NormalizeUniform(double fac);
	void Log();
	void WriteHistogram(const char *prefix, const char *name, const char *suffix);
	void CalcHistogram();
	double GetValue(double x, double y);
	void AddCircle(double x, double y, double r, double cr, double cg, double cb);
	void NormalizeYCount();
	void NormalizeXCount();
	void SwapAxes();
	void Mirror(float plane, int channel);
	void WriteCombinedPlot(const char *prefix, const char *name, const char *suffix);
//	double m_fMathematicaColorOffset;
//	double m_fMathematicaColorScale;
	void Subtract(C2DF *df);
	void CopyFrom(C2DF *df);
	void MakeTensorProduct(C2DF *inp);
	void SetLabelX(const char *s);
	void SetLabelY(const char *s);
	void SetLabelZ(const char *s);
	void CalcMaxEntry();
	void StepY(int y);
	void NormRDF(double n);
	C2DF();
	~C2DF();
	void MultiplyBin(double m);
//	void AddToBin(double x, const CxVector3 &vec);
	void AddToBin(double x, double y);
//	void AddToBin(double x, int y);
	void AddToBin(int x, int y, double val);
	void AddToBin(double x, double y, double val);
	void AddToBin_IntX(int x, double y, double val);
	void AddToBin_IntY(double x, int y, double val);
//	void AddToBin_IntX_fast(int x, double y);
//	void AddToSingleBin(double x, double y, double val);
	void NormalizeBin(double mi, double ma);
	double NormalizeBinIntegral(double val);
	double PPMBin();
	void Write(const char *prefix, const char *name, const char *suffix);
	void WriteGraceBunch(int channel, int graphs, float fac, const char *prefix, const char *name, const char *suffix);
	void WriteCSV(const char *prefix, const char *name, const char *suffix);
	void WriteXProjection(const char *prefix, const char *name, const char *suffix);
	void WriteYProjection(const char *prefix, const char *name, const char *suffix);
	void WriteMathematica(const char *prefix, const char *name, const char *suffix);
	void WriteGnuplotInput(const char *prefix, const char *name, const char *suffix, bool manrange);
	void WriteMathematicaNb(const char *prefix, const char *name, const char *suffix, bool manrange);
	void Create();
	void CorrectRadialDist(int channel);
	void CorrectLiRadialDist(int channel);
	void UnCorrectRadialDist(int channel);
	void CorrectAngle(int channel);
	void UnCorrectAngle(int channel);
//	void AngleCorrectX();

	double m_fMinVal[2];
	double m_fMaxVal[2];

	double m_fFac[2];

	CDF *m_pChannels[2];

	double m_fAspectRatio;
	double m_fMinEntry;
	double m_fMaxEntry;
	char *m_sLabelX;
	char *m_sLabelY;
	char *m_sLabelZ;

	char* m_sName;
	char* m_sShortName;

	double *m_fCountX;
	double *m_fCountY;

	int m_iRes[2];
	int m_iHistogramRes;
	double *m_pBin;
	double *m_pHistogram;
	double m_fBinEntries;
	double m_fSkipEntries;
	unsigned long *m_pStepsY;

	int m_iPlotType;
	int m_iSmoothGrade;
	int m_iInterpolationOrder;
	double m_fPlotExp;
	int m_iExpLegend;
	int m_iColorScale;
	int m_iGPInterpolation;
	bool m_bContourLines;
	int m_iPlotPixel;

	CxObArray m_oaCircles;

	void AddToBin_IntX_fast(int x, double y)
	{ 
		double ry;
		int iy;

		if (y > m_fMaxVal[1])
			return;

		m_fBinEntries++;

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

		m_pBin[ iy    * m_iRes[0] + x] += (1.0-ry);
		m_pBin[(iy+1) * m_iRes[0] + x] +=      ry ;
	}

};


/*inline void C2DF::AddToBin(double x, int y)
{
	BXIN;
	double rx;
	int ix;

	if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;
	rx = ((x-m_fMinVal[0])/(m_fMaxVal[0]-m_fMinVal[0]))*((double)m_iRes[0]-1);
	ix = (int)floor(rx);
	rx -= ix;
	m_pBin[y*m_iRes[0] + ix]            += 1-rx;
	m_pBin[y*m_iRes[0] + ix +1]         += rx;
	m_fCountX[ix] += 1.0-rx;
	m_fCountX[ix+1] += rx;
	m_fCountY[y]++;
	BXOUT;
}*/



#endif
