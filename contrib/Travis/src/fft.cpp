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

#include "fft.h"
#include "tools.h"

#ifdef USE_FFTW

/*************************************************************************
************ FFTW ********************************************************
*************************************************************************/

CFFT::CFFT()
{
	BTIN;
	m_pPlan = NULL;
	BTOUT;
}

CFFT::~CFFT()
{
	BTIN;
	if (m_pPlan != NULL)
	{
		fftwf_destroy_plan(m_pPlan);
		fftwf_free(m_pInput);
		fftwf_free(m_pOutput);
	}
	BTOUT;
}

/*void CFFT::PrepareFFT_R2C(int n)
{
	BTIN;
	int z;
//	m_bInputComplex = false;
//	m_bOutputComplex = true;
	m_iSize = n;
	m_pInput = (float*)fftwf_malloc(sizeof(float) * n);
	for (z=0;z<n;z++)
		m_pInput[z] = 0.0f;
	m_pOutput = (float*)fftwf_malloc(sizeof(fftwf_complex) * n);
	m_pPlan = fftwf_plan_dft_r2c_1d(n, m_pInput, (fftwf_complex*)m_pOutput, FFTW_ESTIMATE);
	BTOUT;
}*/

void CFFT::PrepareFFT_C2C(int n)
{
	BTIN;
	int z;
//	m_bInputComplex = true;
//	m_bOutputComplex = true;
	m_iSize = n;
	m_pInput = (float*)fftwf_malloc(sizeof(fftwf_complex) * n);
	for (z=0;z<n;z++)
		m_pInput[z] = 0.0f;
	m_pOutput = (float*)fftwf_malloc(sizeof(fftwf_complex) * n);
	m_pPlan = fftwf_plan_dft_1d(n, (fftwf_complex*)m_pInput, (fftwf_complex*)m_pOutput, FFTW_FORWARD, FFTW_ESTIMATE);
	BTOUT;
}

void CFFT::PrepareInverseFFT_C2C(int n)
{
	BTIN;
	int z;
//	m_bInputComplex = true;
//	m_bOutputComplex = true;
	m_iSize = n;
	m_pInput = (float*)fftwf_malloc(sizeof(fftwf_complex) * n);
	for (z=0;z<n;z++)
		m_pInput[z] = 0.0f;
	m_pOutput = (float*)fftwf_malloc(sizeof(fftw_complex) * n);
	m_pPlan = fftwf_plan_dft_1d(n, (fftwf_complex*)m_pInput, (fftwf_complex*)m_pOutput, FFTW_BACKWARD, FFTW_ESTIMATE);
	BTOUT;
}

void CFFT::DoFFT()
{
	BTIN;
	fftwf_execute(m_pPlan);
	BTOUT;
}

int CFFT::NextFastSize(int i)
{
	return kiss_fft_next_fast_size(i);
}

#else

/*************************************************************************
************ KISS FFT ****************************************************
*************************************************************************/

CFFT::CFFT()
{
	m_iSize = 0;
	m_pInput = NULL;
	m_pOutput = NULL;
	m_pKISSCfg = NULL;
}

CFFT::~CFFT()
{
	if (m_pKISSCfg != NULL)
		free(m_pKISSCfg);
}

/*void CFFT::PrepareFFT_R2C(int n)
{
	BTIN;
	int z;
	m_iSize = n;
	m_pKISSCfg = kiss_fft_alloc(m_iSize, 0 ,0,0 );
	m_pInput = new float[sizeof(float) * m_iSize * 2];
	for (z=0;z<m_iSize*2;z++)
		m_pInput[z] = 0.0f;
	m_pOutput = new float[sizeof(float) * m_iSize * 2];
//	m_iSign = -1;
	BTOUT;
}*/

void CFFT::PrepareFFT_C2C(int n)
{
	BTIN;
	int z;
	m_iSize = n;
	m_pKISSCfg = kiss_fft_alloc(m_iSize, 0 ,0,0 );

	if (m_pInput != NULL)
		delete[] m_pInput;
	
	try { m_pInput = new float[m_iSize*2]; } catch(...) { m_pInput = NULL; }
	if (m_pInput == NULL) NewException((double)m_iSize*2*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iSize*2;z++)
		m_pInput[z] = 0.0f;

	if (m_pOutput != NULL)
		delete[] m_pOutput;
	
	try { m_pOutput = new float[m_iSize*2]; } catch(...) { m_pOutput = NULL; }
	if (m_pOutput == NULL) NewException((double)m_iSize*2*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
//	m_iSign = 1;
	BTOUT;
}

void CFFT::PrepareInverseFFT_C2C(int n)
{
	BTIN;
	int z;
	m_iSize = n;
	m_pKISSCfg = kiss_fft_alloc(m_iSize, 1 ,0,0 );

	try { m_pInput = new float[m_iSize*2]; } catch(...) { m_pInput = NULL; }
	if (m_pInput == NULL) NewException((double)m_iSize*2*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iSize*2;z++)
		m_pInput[z] = 0.0f;

	try { m_pOutput = new float[m_iSize*2]; } catch(...) { m_pOutput = NULL; }
	if (m_pOutput == NULL) NewException((double)m_iSize*2*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
//	m_iSign = 1;
	BTOUT;
}

void CFFT::DoFFT()
{
	BXIN;
	KISS_FFT();
	BXOUT;
}

void CFFT::KISS_FFT()
{
	kiss_fft(m_pKISSCfg, (kiss_fft_cpx*)m_pInput, (kiss_fft_cpx*)m_pOutput);
}

int CFFT::NextFastSize(int i)
{
	return kiss_fft_next_fast_size(i);
}

#endif


