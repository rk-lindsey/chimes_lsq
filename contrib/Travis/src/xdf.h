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


#ifndef XDF_H
#define XDF_H


#include "xobject.h"
#include "xbytearray.h"
#include "xdoublearray.h"
#include "xfloatarray.h"


class CXDF : public CxObject
{
public:
	void Autocorrelate();
	void InitDeriv();
	CXDF();
	~CXDF();

	void ParseDeriv();

	bool m_bSelf;
	int m_iShowMol;
	int m_iCombinations;
	int m_iDeriv;
	bool m_bDerivAbs;
	bool m_bACF;
	int m_iACFDepth;
	bool m_bACFFFT;
	bool m_bWindowFunction;
	float m_fSpecWaveNumber;
	int m_iMirror;
	int m_iZeroPadding;
	bool m_bACF_DB;
	CxFloatArray **m_pfaACFBuffer;
	CxDoubleArray m_faACF;

	int m_iHistogramRes;
	int m_iResolution;
	CxDoubleArray *m_faData; 
	CxByteArray *m_baDataEnabled; 
	CxDoubleArray **m_pfaDerivBuffer;
	char *m_sName;
	char *m_sShortName;
	char *m_sLabelName;
};

#endif

