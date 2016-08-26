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

#include "acf.h"
#include "travis.h"
#include "maintools.h"


CACF::CACF()
{
	BTIN;
	m_pData = NULL;
	m_pCounter = NULL;
	m_iShowMol = -1;
//	m_b2ndDerivative = false;
	m_oaCache.SetName("CACF::m_oaCache");
	m_oaCCRMatrix.SetName("CACF::m_oaCCRMatrix");
	BTOUT;
}


CACF::~CACF()
{
	BTIN;
	if (m_pData != NULL)
		delete[] m_pData;
	if (m_pCounter != NULL)
		delete[] m_pCounter;
	BTOUT;
}


void CACF::Create()
{
	BTIN;
	int z;

	if (m_pData != NULL)
	{
		delete[] m_pData;
		delete[] m_pCounter;
	}

	mprintf("    ACF: Trying to reserve %s of memory...\n",FormatBytes((double)m_iSize*sizeof(double)+(double)m_iSize*sizeof(long)));

	try { m_pData = new double[m_iSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iSize*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_pCounter = new unsigned long[m_iSize]; } catch(...) { m_pCounter = NULL; }
	if (m_pCounter == NULL) NewException((double)m_iSize*sizeof(unsigned long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iSize;z++)
	{
		m_pData[z] = 0.0f;
		m_pCounter[z] = 0;
	}

	if (m_bSpectrum)
	{
		try { m_pSpectrum = new CSpectrum(); } catch(...) { m_pSpectrum = NULL; }
		if (m_pSpectrum == NULL) NewException((double)sizeof(CSpectrum),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pSpectrum->m_iSize = m_iSize;
		m_pSpectrum->m_fWaveNumber = m_fSpecWaveNumber;
	}
	BTOUT;
}


void CACF::Normalize()
{
	BTIN;
	double d;
	int z;

	if (m_pData[0] < 1E-20f)
		return;

	d = m_pCounter[0] / m_pData[0];

	for (z=0;z<m_iSize;z++)
	{
		if (m_pCounter[z] != 0)
			m_pData[z] *= d / m_pCounter[z];
				else m_pData[z] = 0;
	}
	BTOUT;
}


void CACF::Multiply(double f)
{
	BTIN;
//	double d;
	int z;

	for (z=0;z<m_iSize;z++)
		m_pData[z] *= f;
	BTOUT;
}


void CACF::WriteACF(const char *pre, const char *s, const char *post)
{
	BTIN;
	FILE *a;
	int z;
//	char buf[256];
	CxString buf;

//	sprintf(buf,"%s%s%s",pre,s,post);
	buf.sprintf("%s%s%s",pre,s,post);

	a = OpenFileWrite(buf,true);

	for (z=0;z<m_iSize;z++)
		mfprintf(a,"%.1f;  %g\n",z*g_fTimestepLength,m_pData[z]);

	fclose(a);
	BTOUT;
}


void CACF::Window()
{
/*	BTIN;
	int z;
	for (z=0;z<m_iSize;z++)
		m_pData[z] *= (float)pow(sin((float)z/m_iSize*Pi),2);
	BTOUT;*/

	BTIN;
	int z;
	for (z=0;z<m_iSize;z++)
		m_pData[z] *= (float)pow(cos((float)z/m_iSize*Pi),2);
	BTOUT;
}


void CACF::Transform(CFFT *fft)
{
	BTIN;
	int z;
//	memcpy(fft->m_pInput,m_pData,sizeof(float)*(m_iSize/m_iStride));
/*	if (m_bMirror)
	{
		for (z=0;z<(m_iSize+m_iZeroPadding)*2;z++)
		{
			fft->m_pInput[z*2] = 0;
			fft->m_pInput[z*2+1] = 0;
		}
		for (z=0;z<m_iSize/m_iStride;z++)
		{
			fft->m_pInput[(m_iSize+m_iZeroPadding+z)*2] = (float)m_pData[z];
			fft->m_pInput[(m_iSize+m_iZeroPadding+z)*2+1] = 0;
			fft->m_pInput[(m_iSize+m_iZeroPadding-z)*2] = (float)m_pData[z];
			fft->m_pInput[(m_iSize+m_iZeroPadding-z)*2+1] = 0;
		}
	} else
	{
		for (z=0;z<m_iSize+m_iZeroPadding;z++)
		{
			fft->m_pInput[z*2] = 0;
			fft->m_pInput[z*2+1] = 0;
		}*/

/******** MT Neu ************/
/*		for(z = 0; z < m_iZeroPadding / 2; z++)
		{
			fft->m_pInput[z * 2] = 0.0f;
			fft->m_pInput[z * 2 + 1] = 0.0f;
		}

		for(z = 0; z < m_iSize; z++)
		{
			fft->m_pInput[(m_iZeroPadding / 2 + z) * 2] = (float)m_pData[z];
			fft->m_pInput[(m_iZeroPadding / 2 + z) * 2 + 1] = 0.0f;
		}

		for(z = 0; z < m_iZeroPadding / 2; z++)
		{
			fft->m_pInput[(m_iZeroPadding / 2 + m_iSize + z) * 2] = 0.0f;
			fft->m_pInput[(m_iZeroPadding / 2 + m_iSize + z) * 2 + 1] = 0.0f;
		}*/

/******** MT Neu ************/

		for(z = 0; z < m_iSize/2; z++)
		{
			fft->m_pInput[z * 2] = (float)m_pData[z];
			fft->m_pInput[z * 2 + 1] = 0.0f;
		}

		for(z = 0; z < m_iZeroPadding; z++)
		{
			fft->m_pInput[(m_iSize/2 + z) * 2] = 0.0f;
			fft->m_pInput[(m_iSize/2 + z) * 2 + 1] = 0.0f;
		}

		for(z = 0; z < m_iSize/2; z++)
		{
			fft->m_pInput[(m_iSize/2 + m_iZeroPadding + z) * 2] = (float)m_pData[m_iSize/2 + z];
			fft->m_pInput[(m_iSize/2 + m_iZeroPadding + z) * 2 + 1] = 0.0f;
		}

/******** MB Alt ************/
/*		for (z=0;z<m_iSize;z++)
		{
			fft->m_pInput[z*2] = (float)m_pData[z];
			fft->m_pInput[z*2+1] = 0;
		}
		for (z=m_iSize;z<m_iSize+m_iZeroPadding;z++)
		{
			fft->m_pInput[z*2] = 0;
			fft->m_pInput[z*2+1] = 0;
		}*/
//	}
	fft->DoFFT();
//	m_pSpectrum = new CSpectrum();
//	if (m_bMirror)
//		m_pSpectrum->Create((m_iSize+m_iZeroPadding)*2);
		/*	else */ m_pSpectrum->Create(m_iSize+m_iZeroPadding);
	m_pSpectrum->FromComplex((float*)fft->m_pOutput);
	BTOUT;
}


void CACF::Parse()
{
	BTIN;
//	char buf[256];
	CxString buf;
	int ti, z;
	double tf;
	CFFT tfft;

	if (m_iShowMol != -1)
	{
_atoms:
		mprintf("    Which atoms to observe in %s (e.g. \"C1,C3-5,H\")? [all] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		inpprintf("! Which atoms to observe (e.g. \"C1,C3-5,H\")? [all]\n");
		myget(&buf);
		if (buf.GetLength() == 0)
		{
			m_oAtoms.AddAllAtoms((CMolecule*)g_oaMolecules[m_iShowMol],false);
		} else
		{
			if (!m_oAtoms.ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
			{
				eprintf("Wrong input.\n");
				inpprintf("! Wrong input.\n");
				goto _atoms;
			}
		}
		m_iShowAtomGes = m_oAtoms.m_iAtomGes;
	} else
	{
		m_iShowAtomGes = 0;
		for (z=0;z<g_oaAtoms.GetSize();z++)
		{
			if (z == g_iVirtAtomType)
				continue;
			if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius == 0)
			{
				m_bExcludeR0 = AskYesNo("    Neglect atoms excluded from system (e.g. wannier centers) from VACF (y/n)? [yes] ",true);
				mprintf("\n");
				goto _r0done;
			}
		}
		m_bExcludeR0 = false;
		mprintf("    Observing all atoms in the system.\n\n");
_r0done:;
	}
//_acfdepth:

	mprintf(WHITE,"    Hint: ");
	mprintf("The resolution of the ACF may never be higher than the number of processed steps.\n");
	mprintf("          Suggested is 75 percent of the processed steps, but not more than approx. 16384.\n\n");

	if (g_iTrajSteps != -1)
		ti = (int(g_iTrajSteps*0.75)<5120)?int(g_iTrajSteps*0.75):4096;
			else ti = 4096;

	m_iSize = AskUnsignedInteger("    Enter the resolution (=depth) of the velocity ACF (in time steps): [%d] ",ti,ti);

	if (g_bPowerSpec)
		mprintf("\n    This corresponds to a spectral resolution of %.4f cm^-1.\n",33356.41/g_fTimestepLength/2.0/m_iSize);
	
//	if (m_iSize*g_iGesVirtAtomCount*12.0f*sizeof(double)/1024.0f/1024.0f >= 10.0f)
//		if (!AskYesNo("\nThis needs about %.0f MB Memory. Is this OK (y/n)? [yes] ",true,m_iSize*g_iGesVirtAtomCount*4.0f*sizeof(double)/1024.0f/1024.0f))
//			goto _acfdepth;
	if (g_bACFFFT)
	{
		ti = CalcFFTSize(m_iSize,false);
		if (m_iSize != ti)
		{
			mprintf(WHITE,"    The next \"fast\" size for FFT is %d. Using this instead of %d as size.\n\n",ti,m_iSize);
			m_iSize = ti;
		}
	}
//	m_iStride = AskUnsignedInteger("    Take every n-th time step for the Velocity ACF? [1] ",1);
	m_bDerivative = AskYesNo("    Derive the velocity before autocorrelating (y/n)? [no] ",false);
	m_bMassWeight = AskYesNo("    Weight the autocorrelation functions by atomic mass (y/n)? [yes] ",true);
	m_bWindowFunction = AskYesNo("    Apply window function (Cos^2) to Autocorrelation function (y/n)? [yes] ",true);

//	m_bSpectrum = AskYesNo("    Also calculate power spectrum (= Fourier Transform of velocity ACF) (y/n)? [yes] ",true);

	m_bSpectrum = g_bPowerSpec;

	if (m_bSpectrum)
	{
		tf = 33356.41 / g_fTimestepLength / 2.0;
		mprintf("\n    A time step length of %.1f fs allows a spectral range up to %.1f cm^-1.\n\n",g_fTimestepLength,tf);
		m_fSpecWaveNumber = AskRangeFloat("    Calculate spectrum up to which wave number (cm^-1)? [%.1f cm^-1] ",0,tf,(tf<5000.0)?tf:5000.0,(tf<5000.0)?tf:5000.0);
		m_iMirror = AskUnsignedInteger("    No mirroring (0), short-time enhancing (1) or long-time enhancing (2)? [1] ",1);
		m_iZeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ",m_iSize*3,m_iSize*3);

		ti = CalcFFTSize(m_iSize+m_iZeroPadding,false);
		if (m_iSize+m_iZeroPadding != ti)
		{
			mprintf(WHITE,"    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n\n",ti,ti-m_iSize);
			m_iZeroPadding = ti-m_iSize;
		}

		m_iZeroPadding0 = m_iZeroPadding;

		mprintf("    Zero padding increases the spectral resolution to %.4f cm^-1.\n\n",33356.41/g_fTimestepLength/2.0/(m_iSize+m_iZeroPadding));

		m_bACF_DB = AskYesNo("    Convert intensity axis of spectrum to decibel (y/n)? [no] ",false);
                                                                                                                          
		if (g_bAdvanced2)                                                                                          
			m_bDecomposeModes = AskYesNo("    Decompose power spectrum into normal modes (y/n)? [no] ",false); 
				else m_bDecomposeModes = false;                                                            
	}
	BuildName();
	BTOUT;
}


void CACF::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;

	if (m_iShowMol == -1)
		tmp.sprintf("global");
			else tmp.sprintf("%s_%s",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_oAtoms.m_sName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CACF::BuildAtomList(CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z1t, z1a;
	CxIntArray *a1;

	vec->RemoveAll_KeepSize();
	for (z1t=0;z1t<m_oAtoms.m_baAtomType.GetSize();z1t++)
	{
		a1 = (CxIntArray*)m_oAtoms.m_oaAtoms[z1t];
		for (z1a=0;z1a<a1->GetSize();z1a++)
			vec->Add(((CxIntArray*)obs->m_oaAtomOffset[m_oAtoms.m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
	}
	BXOUT;
}


void CACF::NormalizeCached()
{
	BTIN;
	int z;

	for (z=1;z<m_iSize;z++)
		m_pData[z] /= m_pData[0];
	m_pData[0] = 1.0;
	BTOUT;
}


void CACF::MultiplyCached(double f)
{
	BTIN;
	int z;

	for (z=0;z<m_iSize;z++)
		m_pData[z] *= f;
	BTOUT;
}


void CACF::Mirror(int i)
{
/*	int z2;
	double *tda;

	try { tda = new double[m_iSize*2]; } catch(...) { tda = NULL; }
	if (tda == NULL) NewException((double)m_iSize*2*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (i == 1) // Naehenbetont
	{
		for (z2=0;z2<m_iSize;z2++)
		{
			tda[z2] = m_pData[m_iSize-z2-1];
			tda[m_iSize*2-z2-2] = m_pData[m_iSize-z2-1];
		}
		tda[m_iSize*2-1] = 0;
	} else if (i == 2) // Fernbetont
	{
		for (z2=0;z2<m_iSize;z2++)
		{
			tda[z2] = m_pData[z2];
			tda[m_iSize*2-z2-2] = m_pData[z2];
		}
		tda[m_iSize*2-1] = 0;
	}
	delete[] m_pData;
	m_pData = tda;
	m_iSize *= 2;
	m_iZeroPadding *= 2;*/

/******* MT Neu *********/
/*	int z2;
	double *tda;

	try { tda = new double[m_iSize*2-2]; } catch(...) { tda = NULL; }
	if (tda == NULL) NewException((double)m_iSize*2*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (i == 1) // Naehenbetont
	{
		for(z2 = 0; z2 < m_iSize; z2++)
			tda[z2] = m_pData[m_iSize-1-z2];
		for(z2 = 0; z2 < m_iSize - 2; z2++)
			tda[m_iSize+z2] = m_pData[z2+1];
	}
	delete[] m_pData;
	m_pData = tda;
	m_iSize *= 2;
	m_iSize -= 2;
	m_iZeroPadding *= 2;
	m_iZeroPadding += 2;*/

/******* MT Neu *********/
	int z2;
	double *tda;

	try { tda = new double[m_iSize*2-2]; } catch(...) { tda = NULL; }
	if (tda == NULL) NewException((double)m_iSize*2*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (i == 1) // Naehenbetont
	{
		for(z2 = 0; z2 < m_iSize; z2++)
			tda[z2] = m_pData[z2];
		for(z2 = 0; z2 < m_iSize - 2; z2++)
			tda[m_iSize+z2] = m_pData[m_iSize-2-z2];
	}
	delete[] m_pData;
	m_pData = tda;
	m_iSize *= 2;
	m_iSize -= 2;
	m_iZeroPadding *= 2;
	m_iZeroPadding += 2;
}


