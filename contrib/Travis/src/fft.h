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

#ifndef FFT_H
#define FFT_H

#include "xobject.h"
#include "backtrace.h"
#include "tools.h"
#include "kiss_fft.h"

#ifdef USE_FFTW
#include "fftw3.h"
#endif


class CFFT : public CxObject  
{
public:
	CFFT();
	~CFFT();
	int NextFastSize(int i);
//	void PrepareFFT_R2C(int n);
	void PrepareFFT_C2C(int n);
	void PrepareInverseFFT_C2C(int n);
	void DoFFT();

	int m_iSize;
	float *m_pInput;
	float *m_pOutput;

#ifdef USE_FFTW
	fftwf_plan m_pPlan;
#else
	kiss_fft_cfg m_pKISSCfg;
	void KISS_FFT();
#endif
};


#endif
