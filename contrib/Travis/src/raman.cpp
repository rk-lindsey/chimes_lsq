/*****************************************************************************
    TRAVIS - Trajectory Analyzer and Visualizer
    http://www.travis-analyzer.de/

    Copyright (c) 2009-2016 Martin Brehm
                  2012-2016 Martin Thomas

    This file written by Martin Thomas.

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

#include "raman.h"

#include "acf.h"
#include "df.h"
#include "globalvar.h"
#include "maintools.h"
#include "moltools.h"
#include "timestep.h"
#include "tools.h"
#include "xfloatarray.h"
#include "xobarray.h"

#include <errno.h>
#include <stdio.h>
#include <string.h>

#ifdef TARGET_LINUX
#include <sys/stat.h>
#include <sys/types.h>
#endif

#define BUF_SIZE 4096

static bool g_newRaman;
static char *g_ramanDir;
static bool g_orientAvg;
static float g_fieldStrength;
static int g_step = 0;
static int g_stride;
static char *g_inputTemplate;
static char *g_templateFieldPos;
static char *g_templatePolPos;
static char *g_templateCoordPos;
static char *g_templateStepsPos;

static CxObArray g_ramObserv;

static FILE *g_polFile[3];
static CTimeStep *g_timestep[3];

static FILE *g_reftrajFile;
static int g_steps = 0;

static bool g_ramanCompat = false;


// CRamanDyn::CRamanDyn(int showMol, bool global) {
// 	_global = global;
// 	if(_global) {
// 		m_iShowMol = -1;
// 		m_iMolecules = g_oaSingleMolecules.GetSize();
// 	} else {
// 		m_iShowMol = showMol;
// 		m_iMolecules = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
// 	}
// 
// 	/*** Set Names of dynamic arrays. By M. Brehm */
// 	_dipole0.SetName("CRamanDyn::_dipole0");
// 	char tbuf[256];
// 	for (int tz1=0;tz1<3;tz1++)
// 	{
// 		for (int tz2=0;tz2<3;tz2++)
// 		{
// 			sprintf(tbuf,"CRamanDyn::_polarizability[%d][%d]",tz1,tz2);
// 			_polarizability[tz1][tz2].SetName(tbuf);
// 		}
// 	}
// 	/* End of set names */
// 	
// 	mprintf(YELLOW, ">>> Raman Spectrum >>>\n\n");
// 	
// 	if(_global)
// 		mprintf("    All atoms will be taken.\n\n");
// 	else
// 		mprintf("    All atoms will be taken from the OM %s.\n\n", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
// 	
// 	m_iVecType = 1;
// 	m_iCombinations = 1;
// 	g_bDipole = true;
// 	ParseDipole();
// 	
// 	if(g_iTrajSteps != -1) {
// 		int depth = (int)(g_iTrajSteps / g_stride * 0.75);
// 		if(depth > 4096)
// 			depth = 4096;
// 		m_iDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the ACF (in time steps): [%d] ", depth, depth);
// 	} else {
// 		m_iDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the ACF (in time steps): [2000] ", 2000);
// 	}
// 	
// 	int size = CalcFFTSize(m_iDepth, false);
// 	if(m_iDepth != size) {
// 		mprintf(WHITE,"    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n", size, m_iDepth);
// 		m_iDepth = size;
// 	}
// 	
// 	m_iStride = 1;
// 	m_bSpectrum = true;
// 	
// 	_derivativeOrder = 1;
// 	bool window = true;
// 	bool derive = true;
// 	if(g_bAdvanced2) {
// 		derive = AskYesNo("    Derive the vectors before autocorrelation (y/n)? [yes] ", true);
// 		if(derive)
// 			_derivativeOrder = AskRangeInteger("    Please enter degree of vector derivative (1-2): [1] ", 1, 2, 1);
// 		window = AskYesNo("    Apply window function (Cos^2) to autocorrelation function (y/n)? [yes] ", true);
// 	}
// 	float possibleRange = 33356.41f / g_fTimestepLength / 2.0f / g_stride;
// 	mprintf("    A time step length of %.2f fs with a stride of %d allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, g_stride, possibleRange);
// 	float specWaveNumber = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ", 0, possibleRange, (possibleRange < 5000.0f) ? possibleRange : 5000.0f, (possibleRange < 5000.0f) ? possibleRange : 5000.0f);
// 	int mirror = 1;
// 	int zeroPadding = m_iDepth * 3;
// 	if(g_bAdvanced2) {
// 		mirror = AskRangeInteger("    No mirroring (0) or short-time enhancing (1)? [1] ", 0, 1, 1);
// 		zeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ", m_iDepth * 3, m_iDepth * 3);
// 	}
// 	
// 	size = CalcFFTSize(m_iDepth + zeroPadding, false);
// 	if(m_iDepth + zeroPadding != size) {
// 		mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-m_iDepth);
// 		zeroPadding = size-m_iDepth;
// 	}
// 	
// 	mprintf("    This results in a spectral resolution to %.2f cm^-1.\n\n", 33356.41 / g_fTimestepLength / 2.0f / size);
// 	
// 	try { _isoACF = new CACF(); } catch(...) { _isoACF = NULL; }
// 	if(_isoACF == NULL) NewException((double)sizeof(CACF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	_isoACF->m_iSize = m_iDepth;
// 	_isoACF->m_bSpectrum = true;
// 	_isoACF->m_bDerivative = derive;
// 	_isoACF->m_iDerivative = _derivativeOrder;
// 	_isoACF->m_bWindowFunction = window;
// 	_isoACF->m_fSpecWaveNumber = specWaveNumber;
// 	_isoACF->m_iMirror = mirror;
// 	_isoACF->m_iZeroPadding = zeroPadding;
// 	_isoACF->m_bACF_DB = false;
// 	
// 	try { _anisoACF = new CACF(); } catch(...) { _anisoACF = NULL; }
// 	if(_anisoACF == NULL) NewException((double)sizeof(CACF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	_anisoACF->m_iSize = m_iDepth;
// 	_anisoACF->m_bSpectrum = true;
// 	_anisoACF->m_bDerivative = derive;
// 	_anisoACF->m_iDerivative = _derivativeOrder;
// 	_anisoACF->m_bWindowFunction = window;
// 	_anisoACF->m_fSpecWaveNumber = specWaveNumber;
// 	_anisoACF->m_iMirror = mirror;
// 	_anisoACF->m_iZeroPadding = zeroPadding;
// 	_anisoACF->m_bACF_DB = false;
// 	
// 	BuildName();
// 	
// 	mprintf(YELLOW, "<<< End of Raman Spectrum <<<\n\n");
// }
// 
// CRamanDyn::~CRamanDyn() {
// 	int i, j, k;
// 	
// 	for(i = 0; i < m_iMolecules; i++)
// 		delete (CxFloatArray *)_dipole0[i];
// 	for(i = 0; i < m_iMolecules; i++) {
// 		for(j = 0; j < 3; j++) {
// 			for(k = 0; k < (g_orientAvg ? 3 : 1); k++) {
// 				delete (CxFloatArray *)_polarizability[j][k][i];
// 			}
// 		}
// 	}
// 	delete _isoACF;
// 	delete _anisoACF;
// }
// 
// void CRamanDyn::initialize() {
// 	int i, j, k;
// 	
// 	_isoACF->Create();
// 	_anisoACF->Create();
// 	
// 	if (g_iTrajSteps != -1)
// 		mprintf("    Raman Cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iMolecules * g_iTrajSteps / g_iStride / g_stride * (g_orientAvg ? 9.9 : 3.3) * sizeof(float)));
// 	else
// 		mprintf("    Raman Cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iMolecules * 2000 / g_iStride /g_stride * (g_orientAvg ? 9.9 : 3.3) * sizeof(float)));
// 	for(i = 0; i < m_iMolecules; i++) {
// 		CxVector3 *dipoleVector;
// 		try { dipoleVector = new CxVector3(); } catch(...) { dipoleVector = NULL; }
// 		if(dipoleVector == NULL) NewException((double)sizeof(CxVector3), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		_dipole0.Add(dipoleVector);
// 	}
// 	for(i = 0; i < m_iMolecules; i++) {
// 		for(j = 0; j < 3; j++) {
// 			for(k = 0; k < (g_orientAvg ? 3 : 1); k++) {
// 				CxFloatArray *floatArray;
// 				try { floatArray = new CxFloatArray("CRamanDyn::initialize():floatArray"); } catch(...) { floatArray = NULL; }
// 				if(floatArray == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 				if(g_iTrajSteps != -1) {
// 					floatArray->SetMaxSize((long)(g_iTrajSteps / g_iStride / g_stride * 1.1));
// 					floatArray->SetGrow((long)(g_iTrajSteps / g_iStride / g_stride * 0.1));
// 				} else {
// 					floatArray->SetGrow(10000);
// 				}
// 				_polarizability[j][k].Add(floatArray);
// 			}
// 		}
// 	}
// }
// 
// void CRamanDyn::getDipole0() {
// 	int i;
// 	
// 	for(i = 0; i < m_iMolecules; i++) {
// 		CSingleMolecule *sm;
// 		if(_global)
// 			sm = (CSingleMolecule *)g_oaSingleMolecules[i];
// 		else
// 			sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
// 		*((CxVector3 *)_dipole0[i]) = sm->m_vDipole;
// 	}
// }
// 
// void CRamanDyn::calcPolarizability(int fieldDirection) {
// 	int i;
// 	
// 	for(i = 0; i < m_iMolecules; i++) {
// 		CSingleMolecule *sm;
// 		if(_global)
// 			sm = (CSingleMolecule *)g_oaSingleMolecules[i];
// 		else
// 			sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
// 		for(int j = 0; j < 3; j++) {
// 			((CxFloatArray *)_polarizability[j][fieldDirection][i])->Add(((*((CxVector3 *)_dipole0[i]))[j] - sm->m_vDipole[j]) / g_fieldStrength * 0.393430f); // 0.393430 - Umrechnung von Debye in a.u.
// 		}
// 	}
// }
// 
// void CRamanDyn::finalize() {
// 	int i, j, k, l;
// 	
// 	char filename[BUF_SIZE];
// 	#ifdef TARGET_LINUX
// 	snprintf(filename, BUF_SIZE, "%s/polarizability_%s.dat", g_ramanDir, m_sName);
// 	#else
// 	sprintf(filename, "%s/polarizability_%s.dat", g_ramanDir, m_sName);
// 	#endif
// 	mprintf("    Writing polarizabilities for first molecule to \"%s\"...\n", filename);
// 	FILE *pol_file;
// 	pol_file = OpenFileWrite(filename, false);
// 	for(i = 0; i < ((CxFloatArray *)_polarizability[0][0][0])->GetSize(); i++) {
// 		fprintf(pol_file, "%10.2f", i * g_fTimestepLength * g_iStride * g_stride);
// 		for(j = 0; j < (g_orientAvg ? 3 : 1); j++) {
// 			for(k = 0; k < 3; k++) {
// 				fprintf(pol_file, " %14.8f", (*((CxFloatArray *)_polarizability[k][j][0]))[i]);
// 			}
// 		}
// 		fprintf(pol_file, "\n");
// 	}
// 	fclose(pol_file);
// 	
// 	float step;
// 	switch(_derivativeOrder) {
// 		case 0:
// 			mprintf("    Not deriving polarizabilities.\n");
// 			break;
// 		case 1:
// 			mprintf("    Deriving polarizabilities (1st derivative)...\n");
// 			mprintf(WHITE, "      [");
// 			step = m_iMolecules / 20.0f;
// 			for(i = 0; i < m_iMolecules; i++) {
// 				if(fmodf((float)i, step) < 1.0f)
// 					mprintf(WHITE, "#");
// 				for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - 1; j++) {
// 					for(k = 0; k < 3; k++) {
// 						for(l = 0; l < (g_orientAvg ? 3 : 1); l++) {
// 							(*((CxFloatArray *)_polarizability[k][l][i]))[j] = 0.5f * ((*((CxFloatArray *)_polarizability[k][l][i]))[j+1] - (*((CxFloatArray *)_polarizability[k][l][i]))[j]);
// 						}
// 					}
// 				}
// 			}
// 			mprintf(WHITE, "]\n");
// 			break;
// 		case 2:
// 			mprintf("    Deriving polarizabilities (2nd derivative)...\n");
// 			mprintf(WHITE, "      [");
// 			step = m_iMolecules / 20.0f;
// 			for(i = 0; i < m_iMolecules; i++) {
// 				if(fmodf((float)i, step) < 1.0f)
// 					mprintf(WHITE, "#");
// 				for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - 2; j++) {
// 					for(k = 0; k < 3; k++) {
// 						for(l = 0; l < (g_orientAvg ? 3 : 1); l++) {
// 							(*((CxFloatArray *)_polarizability[k][l][i]))[j] = (*((CxFloatArray *)_polarizability[k][l][i]))[j+2] - 2.0f * (*((CxFloatArray *)_polarizability[k][l][i]))[j+1] + (*((CxFloatArray *)_polarizability[k][l][i]))[j];
// 						}
// 					}
// 				}
// 			}
// 			mprintf(WHITE, "]\n");
// 			break;
// 		default:
// 			mprintf(RED, "Higher derivatives not implemented.\n");
// 			abort();
// 			break;
// 	}
// 	
// 	mprintf("    Processing polarizability tensor components...\n");
// 	step = m_iMolecules / 20.0f;
// 	
// 	CAutoCorrelation *autoCorrelation;
// 	try { autoCorrelation = new CAutoCorrelation(); } catch(...) { autoCorrelation = NULL; }
// 	if(autoCorrelation == NULL) NewException((double)sizeof(CAutoCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	autoCorrelation->Init(((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder, m_iDepth, g_bACFFFT);
// 	
// 	if(g_orientAvg) {
// 		CxFloatArray *isotropyPol, *isotropyACF;
// 		try { isotropyACF = new CxFloatArray("CRamanDyn::finalize():isotropyACF"); } catch(...) { isotropyACF = NULL; }
// 		if(isotropyACF == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		isotropyACF->SetSize(m_iDepth);
// 		try { isotropyPol = new CxFloatArray("CRamanDyn::finalize():isotropyPol"); } catch(...) { isotropyPol = NULL; }
// 		if(isotropyPol == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		isotropyPol->SetSize(((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder);
// 		
// 		mprintf("    Isotropic part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_isoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder; j++) {
// 				(*isotropyPol)[j] = 0.0f;
// 				for(k = 0; k < 3; k++) {
// 					(*isotropyPol)[j] += (*((CxFloatArray *)_polarizability[k][k][i]))[j];
// 				}
// 				(*isotropyPol)[j] /= 3.0f;
// 			}
// 			autoCorrelation->AutoCorrelate(isotropyPol, isotropyACF);
// 			for(j = 0; j < m_iDepth; j++) {
// 				_isoACF->m_pData[j] += (*isotropyACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			if(_global)
// 				_isoACF->m_pData[i] /= g_iSteps / g_stride;
// 			else
// 				_isoACF->m_pData[i] /= g_iSteps / g_stride * m_iMolecules;
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/acf_iso_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/acf_iso_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving ACF as %s...\n", filename);
// 		// 		FILE *acf_file = OpenFileWrite(filename, false);
// 		// 		for(i = 0; i < m_iDepth; i++)
// 		// 			fprintf(acf_file, "%10.2f %12.8f\n", i * g_fTimestepLength * g_iStride * g_stride, _isoACF->m_pData[i]);
// 		// 		fclose(acf_file);
// 		
// 		if(_isoACF->m_iMirror != 0) {
// 			mprintf("    Mirroring ACF...\n");
// 			_isoACF->Mirror(_isoACF->m_iMirror);
// 		}
// 		if(_isoACF->m_bWindowFunction != 0) {
// 			mprintf("    Applying window function to ACF...\n");
// 			_isoACF->Window();
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/acf_iso_%s.mw.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/acf_iso_%s.mw.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving ACF as %s...\n", filename);
// 		// 		acf_file = OpenFileWrite(filename, false);
// 		// 		for(i = 0; i < _isoACF->m_iSize; i++)
// 		// 			fprintf(acf_file, "%10.2f %12.8f\n", i * g_fTimestepLength * g_iStride * g_stride, _isoACF->m_pData[i]);
// 		// 		fclose(acf_file);
// 		
// 		mprintf("    Performing Fourier transformation...\n");
// 		CFFT *fft;
// 		try { fft = new CFFT(); } catch(...) { fft = NULL; }
// 		if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		fft->PrepareFFT_C2C(_isoACF->m_iSize + _isoACF->m_iZeroPadding);
// 		_isoACF->Transform(fft);
// 		delete fft;
// 		_isoACF->m_pSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/spectrum_iso_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/spectrum_iso_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving spectrum as %s...\n", filename);
// 		// 		_isoACF->m_pSpectrum->Write("", filename, "");
// 		
// 		delete isotropyACF;
// 		delete isotropyPol;
// 		
// 		CxFloatArray *anisotropyPol, *anisotropyACF, *anisotropyACFSum;
// 		try { anisotropyPol = new CxFloatArray("CRamanDyn::finalize():anisotropyPol"); } catch(...) { anisotropyPol = NULL; }
// 		if(anisotropyPol == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		anisotropyPol->SetSize(((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder);
// 		try { anisotropyACF = new CxFloatArray("CRamanDyn::finalize():anisotropyACF"); } catch(...) { anisotropyACF = NULL; }
// 		if(anisotropyACF == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		anisotropyACF->SetSize(m_iDepth);
// 		try { anisotropyACFSum = new CxFloatArray("CRamanDyn::finalize():anisotropyACFSum"); } catch(...) { anisotropyACFSum = NULL; }
// 		if(anisotropyACFSum == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		anisotropyACFSum->SetSize(m_iDepth);
// 		
// 		mprintf("    First anisotropic part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_anisoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder; j++) {
// 				(*anisotropyPol)[j] = (*((CxFloatArray *)_polarizability[0][0][i]))[j];
// 				(*anisotropyPol)[j] -= (*((CxFloatArray *)_polarizability[1][1][i]))[j];
// 			}
// 			autoCorrelation->AutoCorrelate(anisotropyPol, anisotropyACF);
// 			for(int j = 0; j < m_iDepth; j++) {
// 				_anisoACF->m_pData[j] += (*anisotropyACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			(*anisotropyACFSum)[i] = 0.5f * _anisoACF->m_pData[i];
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		mprintf("    Second anisotropic part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_anisoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder; j++) {
// 				(*anisotropyPol)[j] = (*((CxFloatArray *)_polarizability[1][1][i]))[j];
// 				(*anisotropyPol)[j] -= (*((CxFloatArray *)_polarizability[2][2][i]))[j];
// 			}
// 			autoCorrelation->AutoCorrelate(anisotropyPol, anisotropyACF);
// 			for(j = 0; j < m_iDepth; j++) {
// 				_anisoACF->m_pData[j] += (*anisotropyACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			(*anisotropyACFSum)[i] += 0.5f * _anisoACF->m_pData[i];
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		mprintf("    Third anisotropic part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_anisoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder; j++) {
// 				(*anisotropyPol)[j] = (*((CxFloatArray *)_polarizability[2][2][i]))[j];
// 				(*anisotropyPol)[j] -= (*((CxFloatArray *)_polarizability[0][0][i]))[j];
// 			}
// 			autoCorrelation->AutoCorrelate(anisotropyPol, anisotropyACF);
// 			for(j = 0; j < m_iDepth; j++) {
// 				_anisoACF->m_pData[j] += (*anisotropyACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			(*anisotropyACFSum)[i] += 0.5f * _anisoACF->m_pData[i];
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		mprintf("    Fourth anisotropic part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_anisoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder; j++) {
// 				(*anisotropyPol)[j] = (*((CxFloatArray *)_polarizability[0][1][i]))[j];
// 				(*anisotropyPol)[j] += (*((CxFloatArray *)_polarizability[1][0][i]))[j];
// 				(*anisotropyPol)[j] *= 0.5f;
// 			}
// 			autoCorrelation->AutoCorrelate(anisotropyPol, anisotropyACF);
// 			for(j = 0; j < m_iDepth; j++) {
// 				_anisoACF->m_pData[j] += (*anisotropyACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			(*anisotropyACFSum)[i] += 3.0f * _anisoACF->m_pData[i];
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		mprintf("    Fifth anisotropic part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_anisoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder; j++) {
// 				(*anisotropyPol)[j] = (*((CxFloatArray *)_polarizability[1][2][i]))[j];
// 				(*anisotropyPol)[j] += (*((CxFloatArray *)_polarizability[2][1][i]))[j];
// 				(*anisotropyPol)[j] *= 0.5f;
// 			}
// 			autoCorrelation->AutoCorrelate(anisotropyPol, anisotropyACF);
// 			for(j = 0; j < m_iDepth; j++) {
// 				_anisoACF->m_pData[j] += (*anisotropyACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			(*anisotropyACFSum)[i] += 3.0f * _anisoACF->m_pData[i];
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		mprintf("    Sixth anisotropic part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_anisoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			for(j = 0; j < ((CxFloatArray *)_polarizability[0][0][0])->GetSize() - _derivativeOrder; j++) {
// 				(*anisotropyPol)[j] = (*((CxFloatArray *)_polarizability[2][0][i]))[j];
// 				(*anisotropyPol)[j] += (*((CxFloatArray *)_polarizability[0][2][i]))[j];
// 				(*anisotropyPol)[j] *= 0.5f;
// 			}
// 			autoCorrelation->AutoCorrelate(anisotropyPol, anisotropyACF);
// 			for(j = 0; j < m_iDepth; j++) {
// 				_anisoACF->m_pData[j] += (*anisotropyACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			(*anisotropyACFSum)[i] += 3.0f * _anisoACF->m_pData[i];
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		for(i = 0; i < m_iDepth; i++) {
// 			if(_global)
// 				_anisoACF->m_pData[i] = (*anisotropyACFSum)[i] / g_iSteps * g_stride;
// 			else
// 				_anisoACF->m_pData[i] = (*anisotropyACFSum)[i] / g_iSteps * g_stride / m_iMolecules;
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/acf_aniso_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/acf_aniso_%s.dat", g_ramanDirm_sName);
// 		// #endif
// 		// 		mprintf("    Saving ACF as %s...\n", filename);
// 		// 		acf_file = OpenFileWrite(filename, false);
// 		// 		for(i = 0; i < m_iDepth; i++)
// 		// 			fprintf(acf_file, "%10.2f %12.8f\n", i * g_fTimestepLength * g_iStride * g_stride, _anisoACF->m_pData[i]);
// 		// 		fclose(acf_file);
// 		
// 		if(_anisoACF->m_iMirror != 0) {
// 			mprintf("    Mirroring ACF...\n");
// 			_anisoACF->Mirror(_anisoACF->m_iMirror);
// 		}
// 		if(_anisoACF->m_bWindowFunction != 0) {
// 			mprintf("    Applying window function to ACF...\n");
// 			_anisoACF->Window();
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/acf_aniso_%s.mw.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/acf_aniso_%s.mw.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving ACF as %s...\n", filename);
// 		// 		acf_file = OpenFileWrite(filename, false);
// 		// 		for(i = 0; i < _anisoACF->m_iSize; i++)
// 		// 			fprintf(acf_file, "%10.2f %12.8f\n", i * g_fTimestepLength * g_iStride * g_stride, _anisoACF->m_pData[i]);
// 		// 		fclose(acf_file);
// 		
// 		mprintf("    Performing Fourier transformation...\n");
// 		try { fft = new CFFT(); } catch(...) { fft = NULL; }
// 		if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		fft->PrepareFFT_C2C(_anisoACF->m_iSize + _anisoACF->m_iZeroPadding);
// 		_anisoACF->Transform(fft);
// 		delete fft;
// 		_anisoACF->m_pSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/spectrum_aniso_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/spectrum_aniso_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving spectrum as %s...\n", filename);
// 		// 		_anisoACF->m_pSpectrum->Write("", filename, "");
// 		
// 		int specSize = _isoACF->m_pSpectrum->m_iSize;
// 		CSpectrum *paraSpectrum, *orthoSpectrum, *sumSpectrum, *depolSpectrum;
// 		try { paraSpectrum = new CSpectrum(); } catch(...) { paraSpectrum = NULL; }
// 		if(paraSpectrum == NULL) NewException((double)sizeof(CSpectrum), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		paraSpectrum->m_fWaveNumber = _isoACF->m_fSpecWaveNumber;
// 		paraSpectrum->Create(specSize);
// 		paraSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		try { orthoSpectrum = new CSpectrum(); } catch(...) { orthoSpectrum = NULL; }
// 		if(orthoSpectrum == NULL) NewException((double)sizeof(CSpectrum), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		orthoSpectrum->m_fWaveNumber = _isoACF->m_fSpecWaveNumber;
// 		orthoSpectrum->Create(specSize);
// 		orthoSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		try { sumSpectrum = new CSpectrum(); } catch(...) { sumSpectrum = NULL; }
// 		if(sumSpectrum == NULL) NewException((double)sizeof(CSpectrum), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		sumSpectrum->m_fWaveNumber = _isoACF->m_fSpecWaveNumber;
// 		sumSpectrum->Create(specSize);
// 		sumSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		try { depolSpectrum = new CSpectrum(); } catch(...) { depolSpectrum = NULL; }
// 		if(depolSpectrum == NULL) NewException((double)sizeof(CSpectrum), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		depolSpectrum->m_fWaveNumber = _isoACF->m_fSpecWaveNumber;
// 		depolSpectrum->Create(specSize);
// 		depolSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		
// 		for(i = 0; i < specSize; i++) {
// 			paraSpectrum->m_pData[i] = _isoACF->m_pSpectrum->m_pData[i] + 4.0f / 45.0f * _anisoACF->m_pSpectrum->m_pData[i];
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/spectrum_para_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/spectrum_para_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving para spectrum as %s...\n", filename);
// 		// 		paraSpectrum->Write("", filename, "");
// 		
// 		for(i = 0; i < specSize; i++) {
// 			orthoSpectrum->m_pData[i] = _anisoACF->m_pSpectrum->m_pData[i] / 15.0f;
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/spectrum_ortho_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/spectrum_ortho_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving ortho spectrum as %s...\n", filename);
// 		// 		orthoSpectrum->Write("", filename, "");
// 		
// 		for(i = 0; i < specSize; i++) {
// 			sumSpectrum->m_pData[i] = paraSpectrum->m_pData[i] + orthoSpectrum->m_pData[i];
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/spectrum_unpol_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/spectrum_unpol_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving unpolarized spectrum as %s...\n", filename);
// 		// 		sumSpectrum->Write("", filename, "");
// 		
// 		for(i = 0; i < specSize; i++) {
// 			float freq = i * paraSpectrum->m_fMaxRWL / paraSpectrum->m_iSize;
// 			float crossSecFactor = powf(g_laser - freq, 4) / freq / (1 - expf(-1.438777f * freq / g_temp));
// 			paraSpectrum->m_pData[i] *= crossSecFactor;
// 			orthoSpectrum->m_pData[i] *= crossSecFactor;
// 			sumSpectrum->m_pData[i] *= crossSecFactor;
// 		}
// 		
// 		for(i = 0; i < specSize; i++) {
// 			depolSpectrum->m_pData[i] = orthoSpectrum->m_pData[i] / paraSpectrum->m_pData[i];
// 		}
// 		
// 		#ifdef TARGET_LINUX
// 		snprintf(filename, BUF_SIZE, "%s/spectrum_para_%s.csv", g_ramanDir, m_sName);
// 		#else
// 		sprintf(filename, "%s/spectrum_para_%s.csv", g_ramanDir, m_sName);
// 		#endif
// 		mprintf("    Saving parallel spectrum as %s...\n", filename);
// 		paraSpectrum->Write("", filename, "");
// 		
// 		#ifdef TARGET_LINUX
// 		snprintf(filename, BUF_SIZE, "%s/spectrum_ortho_%s.csv", g_ramanDir, m_sName);
// 		#else
// 		sprintf(filename, "%s/spectrum_ortho_%s.csv", g_ramanDir, m_sName);
// 		#endif
// 		mprintf("    Saving orthogonal spectrum as %s...\n", filename);
// 		orthoSpectrum->Write("", filename, "");
// 		
// 		#ifdef TARGET_LINUX
// 		snprintf(filename, BUF_SIZE, "%s/spectrum_unpol_%s.csv", g_ramanDir, m_sName);
// 		#else
// 		sprintf(filename, "%s/spectrum_unpol_%s.csv", g_ramanDir, m_sName);
// 		#endif
// 		mprintf("    Saving unpolarized spectrum as %s...\n", filename);
// 		sumSpectrum->Write("", filename, "");
// 		
// 		#ifdef TARGET_LINUX
// 		snprintf(filename, BUF_SIZE, "%s/depol_ratio_%s.csv", g_ramanDir, m_sName);
// 		#else
// 		sprintf(filename, "%s/depol_ratio_%s.csv", g_ramanDir, m_sName);
// 		#endif
// 		mprintf("    Saving depolarization ratio as %s...\n", filename);
// 		depolSpectrum->Write("", filename, "");
// 		
// 		delete anisotropyACFSum;
// 		delete anisotropyACF;
// 		delete anisotropyPol;
// 		
// 		delete paraSpectrum;
// 		delete orthoSpectrum;
// 		delete sumSpectrum;
// 		delete depolSpectrum;
// 	} else {
// 		CxFloatArray *ACF;
// 		try { ACF = new CxFloatArray("CRamanDyn::finalize():ACF"); } catch(...) { ACF = NULL; }
// 		if(ACF == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		ACF->SetSize(m_iDepth);
// 		
// 		mprintf("    Parallel part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_isoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			autoCorrelation->AutoCorrelate((CxFloatArray *)_polarizability[0][0][i], ACF);
// 			for(int j = 0; j < m_iDepth; j++) {
// 				_isoACF->m_pData[j] += (*ACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			if(_global)
// 				_isoACF->m_pData[i] /= g_iSteps / g_stride;
// 			else
// 				_isoACF->m_pData[i] /= g_iSteps / g_stride * m_iMolecules;
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/acf_para_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/acf_para_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving ACF as %s...\n", filename);
// 		// 		FILE *acf_file = OpenFileWrite(filename, false);
// 		// 		for(i = 0; i < m_iDepth; i++)
// 		// 			fprintf(acf_file, "%10.2f %12.8f\n", i * g_fTimestepLength * g_iStride * g_stride, _isoACF->m_pData[i]);
// 		// 		fclose(acf_file);
// 		
// 		if(_isoACF->m_iMirror != 0) {
// 			mprintf("    Mirroring ACF...\n");
// 			_isoACF->Mirror(_isoACF->m_iMirror);
// 		}
// 		if(_isoACF->m_bWindowFunction != 0) {
// 			mprintf("    Applying window function to ACF...\n");
// 			_isoACF->Window();
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/acf_para_%s.mw.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/acf_para_%s.mw.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving ACF as %s...\n", filename);
// 		// 		acf_file = OpenFileWrite(filename, false);
// 		// 		for(i = 0; i < _isoACF->m_iSize; i++)
// 		// 			fprintf(acf_file, "%10.2f %12.8f\n", i * g_fTimestepLength * g_iStride * g_stride, _isoACF->m_pData[i]);
// 		// 		fclose(acf_file);
// 		
// 		mprintf("    Performing Fourier transformation...\n");
// 		CFFT *fft;
// 		try { fft = new CFFT(); } catch(...) { fft = NULL; }
// 		if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		fft->PrepareFFT_C2C(_isoACF->m_iSize + _isoACF->m_iZeroPadding);
// 		_isoACF->Transform(fft);
// 		delete fft;
// 		_isoACF->m_pSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/spectrum_para_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/spectrum_para_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving spectrum as %s...\n", filename);
// 		// 		_isoACF->m_pSpectrum->Write("", filename, "");
// 		
// 		mprintf("    Orthogonal part...\n");
// 		mprintf(WHITE, "      [");
// 		for(i = 0; i < m_iDepth; i++) {
// 			_anisoACF->m_pData[i] = 0.0f;
// 		}
// 		for(i = 0; i < m_iMolecules; i++) {
// 			if(fmodf((float)i, step) < 1.0f)
// 				mprintf(WHITE, "#");
// 			autoCorrelation->AutoCorrelate((CxFloatArray *)_polarizability[1][0][i], ACF);
// 			for(int j = 0; j < m_iDepth; j++) {
// 				_anisoACF->m_pData[j] += (*ACF)[j];
// 			}
// 		}
// 		for(i = 0; i < m_iDepth; i++) {
// 			if(_global)
// 				_anisoACF->m_pData[i] /= g_iSteps / g_stride;
// 			else
// 				_anisoACF->m_pData[i] /= g_iSteps / g_stride * m_iMolecules;
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/acf_ortho_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/acf_ortho_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving ACF as %s...\n", filename);
// 		// 		acf_file = OpenFileWrite(filename, false);
// 		// 		for(i = 0; i < m_iDepth; i++)
// 		// 			fprintf(acf_file, "%10.2f %12.8f\n", i * g_fTimestepLength * g_iStride * g_stride, _anisoACF->m_pData[i]);
// 		// 		fclose(acf_file);
// 		
// 		if(_anisoACF->m_iMirror != 0) {
// 			mprintf("    Mirroring ACF...\n");
// 			_anisoACF->Mirror(_isoACF->m_iMirror);
// 		}
// 		if(_anisoACF->m_bWindowFunction != 0) {
// 			mprintf("    Applying window function to ACF...\n");
// 			_anisoACF->Window();
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/acf_ortho_%s.mw.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/acf_ortho_%s.mw.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving ACF as %s...\n", filename);
// 		// 		acf_file = OpenFileWrite(filename, false);
// 		// 		for(i = 0; i < _isoACF->m_iSize; i++)
// 		// 			fprintf(acf_file, "%10.2f %12.8f\n", i * g_fTimestepLength * g_iStride * g_stride, _anisoACF->m_pData[i]);
// 		// 		fclose(acf_file);
// 		
// 		mprintf("    Performing Fourier transformation...\n");
// 		try { fft = new CFFT(); } catch(...) { fft = NULL; }
// 		if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		fft->PrepareFFT_C2C(_anisoACF->m_iSize + _anisoACF->m_iZeroPadding);
// 		_anisoACF->Transform(fft);
// 		delete fft;
// 		_anisoACF->m_pSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/spectrum_ortho_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/spectrum_ortho_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving spectrum as %s...\n", filename);
// 		// 		_anisoACF->m_pSpectrum->Write("", filename, "");
// 		
// 		delete ACF;
// 		
// 		int specSize = _isoACF->m_pSpectrum->m_iSize;
// 		CSpectrum *sumSpectrum, *depolSpectrum;
// 		try { sumSpectrum = new CSpectrum(); } catch(...) { sumSpectrum = NULL; }
// 		if(sumSpectrum == NULL) NewException((double)sizeof(CSpectrum), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		sumSpectrum->m_fWaveNumber = _isoACF->m_fSpecWaveNumber;
// 		sumSpectrum->Create(specSize);
// 		sumSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		try { depolSpectrum = new CSpectrum(); } catch(...) { depolSpectrum = NULL; }
// 		if(depolSpectrum == NULL) NewException((double)sizeof(CSpectrum), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		depolSpectrum->m_fWaveNumber = _isoACF->m_fSpecWaveNumber;
// 		depolSpectrum->Create(specSize);
// 		depolSpectrum->SetMaxRWL(1e15f/299792458.0f/100.0f/g_fTimestepLength/g_iStride/g_stride);
// 		
// 		for(i = 0; i < specSize; i++) {
// 			sumSpectrum->m_pData[i] = _isoACF->m_pSpectrum->m_pData[i] + _anisoACF->m_pSpectrum->m_pData[i];
// 		}
// 		
// 		// #ifdef TARGET_LINUX
// 		// 		snprintf(filename, BUF_SIZE, "%s/spectrum_unpol_%s.dat", g_ramanDir, m_sName);
// 		// #else
// 		// 		sprintf(filename, "%s/spectrum_unpol_%s.dat", g_ramanDir, m_sName);
// 		// #endif
// 		// 		mprintf("    Saving unpolarized spectrum as %s...\n", filename);
// 		// 		sumSpectrum->Write("", filename, "");
// 		
// 		for(i = 0; i < specSize; i++) {
// 			float freq = i * _isoACF->m_pSpectrum->m_fMaxRWL / _isoACF->m_pSpectrum->m_iSize;
// 			float crossSecFactor = powf(g_laser - freq, 4) / freq / (1 - expf(-1.438777f * freq / g_temp));
// 			_isoACF->m_pSpectrum->m_pData[i] *= crossSecFactor;
// 			_anisoACF->m_pSpectrum->m_pData[i] *= crossSecFactor;
// 			sumSpectrum->m_pData[i] *= crossSecFactor;
// 		}
// 		
// 		for(i = 0; i < specSize; i++) {
// 			depolSpectrum->m_pData[i] = _anisoACF->m_pSpectrum->m_pData[i] / _isoACF->m_pSpectrum->m_pData[i];
// 		}
// 		
// 		#ifdef TARGET_LINUX
// 		snprintf(filename, BUF_SIZE, "%s/spectrum_para_%s.csv", g_ramanDir, m_sName);
// 		#else
// 		sprintf(filename, "%s/spectrum_para_%s.csv", g_ramanDir, m_sName);
// 		#endif
// 		mprintf("    Saving parallel spectrum as %s...\n", filename);
// 		_isoACF->m_pSpectrum->Write("", filename, "");
// 		
// 		#ifdef TARGET_LINUX
// 		snprintf(filename, BUF_SIZE, "%s/spectrum_ortho_%s.csv", g_ramanDir, m_sName);
// 		#else
// 		sprintf(filename, "%s/spectrum_ortho_%s.csv", g_ramanDir, m_sName);
// 		#endif
// 		mprintf("    Saving orthogonal spectrum as %s...\n", filename);
// 		_anisoACF->m_pSpectrum->Write("", filename, "");
// 		
// 		#ifdef TARGET_LINUX
// 		snprintf(filename, BUF_SIZE, "%s/spectrum_unpol_%s.csv", g_ramanDir, m_sName);
// 		#else
// 		sprintf(filename, "%s/spectrum_unpol_%s.csv", g_ramanDir, m_sName);
// 		#endif
// 		mprintf("    Saving unpolarized spectrum as %s...\n", filename);
// 		sumSpectrum->Write("", filename, "");
// 		
// 		#ifdef TARGET_LINUX
// 		snprintf(filename, BUF_SIZE, "%s/depol_ratio_%s.csv", g_ramanDir, m_sName);
// 		#else
// 		sprintf(filename, "%s/depol_ratio_%s.csv", g_ramanDir, m_sName);
// 		#endif
// 		mprintf("    Saving depolarization ratio as %s...\n", filename);
// 		depolSpectrum->Write("", filename, "");
// 		
// 		delete sumSpectrum;
// 		delete depolSpectrum;
// 	}
// 	
// 	delete autoCorrelation;
// }

CRamanObservation::CRamanObservation(bool global) {
	int i;
	if(global) {
		m_iShowMol = -1;
		m_iShowMolCount = g_oaSingleMolecules.GetSize();
		_name = new char[7];
		sprintf(_name, "global");
	} else {
		char buf[BUF_SIZE];
		char buf2[BUF_SIZE];
		size_t remaining = BUF_SIZE;
		if(g_oaMolecules.GetSize() > 1) {
#ifdef TARGET_LINUX
			remaining -= snprintf(buf, remaining, "    Which molecule should be observed (");
#else
			remaining -= sprintf(buf, "    Which molecule should be observed (");
#endif
			for(i = 0; i < g_oaMolecules.GetSize(); i++) {
				if(remaining < 1)
					break;
#ifdef TARGET_LINUX
				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#else
				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;
				if(i < g_oaMolecules.GetSize() - 1) {
#ifdef TARGET_LINUX
					length = snprintf(buf2, remaining, ", ");
#else
					length = sprintf(buf2, ", ");
#endif
					strncat(buf, buf2, remaining - 1);
					remaining -= length;
				}
			}
			strncat(buf, ")? ", remaining - 1);
			m_iShowMol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
		} else {
			m_iShowMol = 0;
			mprintf("    Observing molecule %s.\n", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		}
		m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
		_name = new char[strlen(((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName) + 1];
		strcpy(_name, ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		mprintf("\n");
	}
	
	if (!g_ramanCompat) {
		mprintf("    If the full polarizability tensor is available, you can do orientational averaging.\n    Otherwise, please specify the field polarization used to calculate the polarizability.\n\n");
		int def = 0;
		if (g_iPolarizabilityConf[0] > 0 && g_iPolarizabilityConf[1] > 0 && g_iPolarizabilityConf[2] > 0) {
			def = 1;
		} else {
			if (g_iPolarizabilityConf[0] > 0)
				def = 2;
			else if (g_iPolarizabilityConf[1] > 0)
				def = 3;
			else if (g_iPolarizabilityConf[2] > 0)
				def = 4;
		}
		if (def > 0) {
			m_fieldMode = AskRangeInteger("    Use orientational averaging (1), field along x axis (2), field along y axis (3), or field along z axis (4)? [%d] ", 1, 4, def, def);
		} else {
			m_fieldMode = AskRangeInteger_ND("    Use orientational averaging (1), field along x axis (2), field along y axis (3), or field along z axis (4)? ", 1, 4);
		}
		mprintf("\n    The scattering cross sections are calculated for:\n");
		if (m_fieldMode == 1 || m_fieldMode == 2)
			mprintf("    x polarized incident laser beam propagating along the y axis with detection in z direction.\n");
		else if (m_fieldMode == 3)
			mprintf("    y polarized incident laser beam propagating along the z axis with detection in x direction.\n");
		else if (m_fieldMode == 4)
			mprintf("    z polarized incident laser beam propagating along the x axis with detection in y direction.\n");
		mprintf("\n");
	} else {
		if (g_orientAvg)
			m_fieldMode = 1;
		else
			m_fieldMode = 2;
	}
	
	if(g_iTrajSteps != -1) {
		_correlationDepth = (int)(0.75 * g_iTrajSteps);
		if(_correlationDepth > 4096)
			_correlationDepth = 4096;
		if(g_fTimestepLength * g_stride > 1.0f)
			_correlationDepth = 2048;
		if(g_fTimestepLength * g_stride > 2.0f)
			_correlationDepth = 1024;
		_correlationDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the ACF (in time steps): [%d] ", _correlationDepth, _correlationDepth);
	} else {
		_correlationDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the ACF (in time steps): [256] ", 256);
	}
	int size = CalcFFTSize(_correlationDepth, false);
	if(_correlationDepth != size) {
		mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n", size, _correlationDepth);
		_correlationDepth = size;
	}
	
	if(g_bAdvanced2) {
		_windowFunction = AskRangeInteger("    Window function: cos^2(a*t) (1), exp(-t/a) (2), exp(-(t/a)^2) (3) [1] ", 1, 3, 1);
		if(_windowFunction == 1) {
			mprintf("    The parameter \"a\" is chosen according to the correlation depth.\n");
			_windowFunctionParameter = 0;
		} else if(_windowFunction == 2) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in time steps): [%d] ", _correlationDepth / 4, _correlationDepth / 4);
		} else if(_windowFunction == 3) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in time steps): [%d] ", _correlationDepth / 2, _correlationDepth / 2);
		} else {
			eprintf("This is impossible.\n");
			abort();
		}
	} else {
		_windowFunction = 1;
		_windowFunctionParameter = 0;
	}
	
	if(g_bAdvanced2) {
		_zeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ", _correlationDepth * 3, _correlationDepth * 3);
		size = CalcFFTSize(_correlationDepth + _zeroPadding, false);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
	} else {
		_zeroPadding = _correlationDepth * 3;
		size = CalcFFTSize(_correlationDepth + _zeroPadding, false);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
	}
	
	float possibleRange = 33356.41f / g_fTimestepLength / g_stride / 2.0f;
	_specResolution = possibleRange / (_correlationDepth + _zeroPadding);
	mprintf("    This results in a spectral resolution of %.2f cm^-1.\n", _specResolution);
	mprintf("\n    A time step length of %.2f fs with a stride of %d allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, g_stride, possibleRange);
	float specLimit = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ", 0, possibleRange, (possibleRange < 5000.0f) ? possibleRange : 5000.0f, (possibleRange < 5000.0f) ? possibleRange : 5000.0f);
	_specSize = (int)(specLimit / _specResolution);
	mprintf("\n");
	
	if(g_bAdvanced2) {
		_finiteDifferenceCorrection = AskYesNo("    Apply finite difference correction (y/n)? [yes] ", true);
		mprintf("\n");
	} else {
		_finiteDifferenceCorrection = true;
	}
	
	if(g_bAdvanced2) {
		_saveACF = AskYesNo("    Save autocorrelation functions (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_saveACF = false;
	}
	
	if(g_bAdvanced2 && m_iShowMol == -1) {
		_includeCross = AskYesNo("    Include also cross-correlations (y/n)? [no] ", false);
		if(_includeCross) {
			mprintf(RED, "    This is not implemented.\n");
		}
		mprintf("\n");
	} else {
		_includeCross = false;
	}
	
	{
		_quantumCorrection = 1;
	}
	
	_laserFreq = AskFloat("    Calculate scattering cross section for which laser frequency (cm^-1)? [20000.0] ", 20000.0f);
	_temperature = AskFloat("    Calculate scattering cross section for which temperature (K)? [300.0] ", 300.0f);
	mprintf("\n");
	
	// 	try { _ramanDyn = new CRamanDyn(m_iShowMol, global); } catch(...) { _ramanDyn = NULL; }
// 	if(_ramanDyn == NULL) NewException((double)sizeof(CRamanDyn), __FILE__, __LINE__, __PRETTY_FUNCTION__);

}

CRamanObservation::~CRamanObservation() {
	delete[] _name;

// 	delete _ramanDyn;
}

void CRamanObservation::initialize() {
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	if (g_ramanCompat) {
		_dipoleZero.SetSize(m_iShowMolCount);
		if(g_orientAvg) {
			mprintf("    Polarizability cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * 3.0 * n * sizeof(float)));
		} else {
			mprintf("    Polarizability cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * 9.0 * n * sizeof(float)));
		}
		int i, j, k;
		for(i = 0; i < m_iShowMolCount; i++) {
			for(j = 0; j < 3; j++) {
				for(k = 0; k < (g_orientAvg ? 3 : 1); k++) {
					CxFloatArray *a;
					try { a = new CxFloatArray(); } catch(...) { a = NULL; }
					if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
					a->SetMaxSize(n);
					a->SetGrow(n / 10);
					_polarizabilityCache.Add(a);
				}
			}
		}
	} else {
		mprintf("    Polarizability cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * 9.0 * n * sizeof(float)));
		int i;
		for (i = 0; i < m_iShowMolCount; i++) {
			int j;
			for (j = 0; j < 9; j++) {
				CxFloatArray *a;
				try { a = new CxFloatArray; } catch (...) { a = NULL; }
				if (a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				a->SetMaxSize(n);
				a->SetGrow(n / 10);
				_polarizabilityCache.Add(a);
			}
		}
	}
	
// 	_ramanDyn->initialize();
}

void CRamanObservation::getDipoleZero() {
	int i;
	if(m_iShowMol == -1) {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
			_dipoleZero[i] = sm->m_vDipole;
		}
	} else {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
			_dipoleZero[i] = sm->m_vDipole;
		}
	}
	
// 	_ramanDyn->getDipole0();
}

void CRamanObservation::calcPolarizability(int fieldDirection) {
	int i, j;
	if(m_iShowMol == -1) {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
			for(j = 0; j < 3; j++) {
				((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * fieldDirection + m_iShowMolCount * j + i])->Add((_dipoleZero[i][j] - sm->m_vDipole[j]) / g_fieldStrength * 0.393430f); // Conversion from Debye to a.u.
			}
		}
	} else {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
			for(j = 0; j < 3; j++) {
				((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * fieldDirection + m_iShowMolCount * j + i])->Add((_dipoleZero[i][j] - sm->m_vDipole[j]) / g_fieldStrength * 0.393430f); // Conversion from Debye to a.u.
			}
		}
	}
	
// 	_ramanDyn->calcPolarizability(fieldDirection);
}

void CRamanObservation::process() {
	if (m_iShowMol == -1) {
		int i;
		for (i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
			int j;
			for (j = 0; j < 9; j++) {
				((CxFloatArray *)_polarizabilityCache[m_iShowMolCount * j + i])->Add(sm->m_polarizability[j]);
			}
		}
	} else {
		int i;
		for (i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
			int j;
			for (j = 0; j < 9; j++) {
				((CxFloatArray *)_polarizabilityCache[m_iShowMolCount * j + i])->Add(sm->m_polarizability[j]);
			}
		}
	}
}

void CRamanObservation::finalize() {
	int i, j, k, l;
	CxString filename;
	
	if(_saveACF) {
		if (g_ramanCompat)
			filename.sprintf("%s/polarizability_%s.dat", g_ramanDir, _name);
		else
			filename.sprintf("polarizability_%s.dat", _name);
		mprintf("    Saving polarizabilities for first molecule as %s...\n", (const char *)filename);
		FILE *pol_file;
		pol_file = OpenFileWrite(filename, false);
		for(i = 0; i < ((CxFloatArray *)_polarizabilityCache[0])->GetSize(); i++) {
			fprintf(pol_file, "%.2f;", i * g_fTimestepLength * g_iStride * g_stride);
			if (g_ramanCompat) {
				for(j = 0; j < 3; j++) {
					for(k = 0; k < (g_orientAvg ? 3 : 1); k++) {
						fprintf(pol_file, " %.8G;", ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * k + m_iShowMolCount * j])->GetAt(i));
					}
				}
			} else {
				for (j = 0; j < 9; j++) {
					fprintf(pol_file, " %.8G;", ((CxFloatArray *)_polarizabilityCache[m_iShowMolCount * j])->GetAt(i));
				}
			}
			fprintf(pol_file, "\n");
		}
		fclose(pol_file);
	}
	
	int n = ((CxFloatArray *)_polarizabilityCache[0])->GetSize() - 2;
	float step = (float)m_iShowMolCount / 20.0f;

	if (n < _correlationDepth)
	{
		eprintf("\nError: Autocorrelation depth is %d, but only %d timesteps evaluated.\n",_correlationDepth,n);
		eprintf("       Reduce depth or increase trajectory length.\n\n");
		abort();
	}

	mprintf("    Deriving polarizabilities...\n");
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmodf(i, step) < 1.0f)
			mprintf(WHITE, "#");
		if (g_ramanCompat) {
			for(j = 0; j < 3; j++) {
				for(k = 0; k < (g_orientAvg ? 3 : 1); k++) {
					for(l = 0; l < n; l++) {
						((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * k + m_iShowMolCount * j + i])->GetAt(l) = 0.5f * (((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * k + m_iShowMolCount * j + i])->GetAt(l+2) - ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * k + m_iShowMolCount * j + i])->GetAt(l)) / g_fTimestepLength / g_stride;
					}
				}
			}
		} else {
			for (j = 0; j < 9; j++) {
				for (k = 0; k < n; k++) {
					((CxFloatArray *)_polarizabilityCache[m_iShowMolCount * j + i])->GetAt(k) = 0.5f * (((CxFloatArray *)_polarizabilityCache[m_iShowMolCount * j + i])->GetAt(k + 2) - ((CxFloatArray *)_polarizabilityCache[m_iShowMolCount * j + i])->GetAt(k)) / g_fTimestepLength;
				}
			}
		}
	}
	mprintf(WHITE, "]\n");
	
	mprintf("    Processing polarizability tensor components...\n");
	
	CxFloatArray *acf;
	try { acf = new CxFloatArray(); } catch(...) { acf = NULL; }
	if(acf == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	acf->SetSize(_correlationDepth);
	for(i = 0; i < _correlationDepth; i++) {
		acf->GetAt(i) = 0.0f;
	}
	CAutoCorrelation *ac;
	try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
	if(ac == NULL) NewException((double)sizeof(CAutoCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	ac->Init(n, _correlationDepth, g_bACFFFT);
	CxFloatArray *temp;
	try { temp = new CxFloatArray(); } catch(...) { temp = NULL; }
	if(temp == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp->SetSize(_correlationDepth);
	CxFloatArray *temp2;
	try { temp2 = new CxFloatArray(); } catch(...) { temp2 = NULL; }
	if(temp2 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp2->SetSize(n);
	CxFloatArray *temp3;
	try { temp3 = new CxFloatArray(); } catch(...) { temp3 = NULL; }
	if(temp3 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp3->SetSize(_correlationDepth);
	CFFT *fft;
	try { fft = new CFFT(); } catch(...) { fft = NULL; }
	if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	CxFloatArray *spectrum1;
	try { spectrum1 = new CxFloatArray(); } catch(...) { spectrum1 = NULL; }
	if(spectrum1 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	spectrum1->SetSize(_specSize);
	CxFloatArray *spectrum2;
	try { spectrum2 = new CxFloatArray(); } catch(...) { spectrum2 = NULL; }
	if(spectrum2 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	spectrum2->SetSize(_specSize);
	
// 	if(g_orientAvg) {
	if (m_fieldMode == 1) {
		mprintf("    Isotropic part...\n");
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			for(j = 0; j < n; j++) {
				temp2->GetAt(j) = 0.0f;
				for(k = 0; k < 3; k++) {
					temp2->GetAt(j) += ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * k + m_iShowMolCount * k + i])->GetAt(j);
				}
				temp2->GetAt(j) /= 3.0f;
			}
			ac->AutoCorrelate(temp2, temp);
			for(j = 0; j < _correlationDepth; j++) {
				acf->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				acf->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		
		if(_saveACF) {
			if (g_ramanCompat)
				filename.sprintf("%s/acf_iso_%s.csv", g_ramanDir, _name);
			else
				filename.sprintf("raman_acf_iso_%s.csv", _name);
			mprintf("    Saving autocorrelation function as %s...\n", (const char *)filename);
			FILE *acf_file = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(acf_file, "%.2f; %.10G\n", i * g_fTimestepLength * g_stride, acf->GetAt(i));
			}
			fclose(acf_file);
		}
		
		temp->CopyFrom(acf);
		
		if(_windowFunction == 1) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= powf(cosf((float)i / (temp->GetSize() - 1) / 2.0f * Pi), 2.0f);
			}
		} else if(_windowFunction == 2) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= expf(-(float)i / _windowFunctionParameter);
			}
		} else if(_windowFunction == 3) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= expf(-(float)i * i / _windowFunctionParameter / _windowFunctionParameter);
			}
		} else if(_windowFunction != 0) {
			eprintf("Unknown window function.\n");
			abort();
		}
		
		if(_saveACF) {
			if (g_ramanCompat)
				filename.sprintf("%s/acf_windowed_iso_%s.csv", g_ramanDir, _name);
			else
				filename.sprintf("raman_acf_windowed_iso_%s.csv", _name);
			mprintf("    Saving windowed autocorrelation function as %s...\n", (const char *)filename);
			FILE *acf_file = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(acf_file, "%.2f; %.10G\n", i * g_fTimestepLength * g_stride, acf->GetAt(i));
			}
			fclose(acf_file);
		}
		
		if(_zeroPadding > 0) {
			for(i = 0; i < _zeroPadding; i++) {
				temp->Add(0.0f);
			}
		}
		
		int oldSize = temp->GetSize();
		temp->SetSize(2 * oldSize);
		for(i = 1; i < oldSize; i++) {
			temp->GetAt(oldSize + i) = temp->GetAt(oldSize - i);
		}
		temp->GetAt(oldSize) = 0.0f;
		
		mprintf("    Performing Fourier transformation...\n");
		fft->PrepareFFT_C2C(temp->GetSize());
		for(i = 0; i < temp->GetSize(); i++) {
			fft->m_pInput[2*i] = temp->GetAt(i);
			fft->m_pInput[2*i+1] = 0.0f;
		}
		fft->DoFFT();
		
		for(i = 0; i < _specSize; i++) {
			float freq = i * _specResolution;
			spectrum1->GetAt(i) = 4.160440e-18f * powf(_laserFreq - freq, 4.0f) / freq / (1 - expf(-1.438777f * freq / _temperature)) * fft->m_pOutput[2*i] * g_fTimestepLength * g_stride; // Output in 1e-30*m^2*cm
		}
		spectrum1->GetAt(0) = 0.0f;
		
		if(_finiteDifferenceCorrection) {
			float f = _specResolution * g_fTimestepLength * g_stride * 1.883652e-4f;
			for(i = 1; i < _specSize; i++) {
				spectrum1->GetAt(i) *= powf(f * i / sinf(f * i), 2.0f); // Divide by sinc function to correct finite difference derivation
			}
		}
		
		if(_quantumCorrection != 1) {
			for(i = 0; i < _specSize; i++) {
				float factor = 1.0f;
				if(_quantumCorrection == 2) {
					factor = 2.0f / (1.0f + expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature)) / (1.438777f * _specResolution * i / _quantumCorrectionTemperature) * (1.0f - expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature));
				} else if(_quantumCorrection == 3) {
					factor = expf(0.719388f * _specResolution * i / _quantumCorrectionTemperature) / (1.438777f * _specResolution * i / _quantumCorrectionTemperature) * (1.0f - expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature));
				}
				spectrum1->GetAt(i) *= factor;
			}
		}
		
		mprintf("    First anisotropic part...\n");
		for(i = 0; i < _correlationDepth; i++) {
			temp3->GetAt(i) = 0.0f;
		}
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			for(j = 0; j < n; j++) {
				temp2->GetAt(j) = ((CxFloatArray *)_polarizabilityCache[i])->GetAt(j) - ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount + m_iShowMolCount + i])->GetAt(j);
			}
			ac->AutoCorrelate(temp2, temp);
			for(j = 0; j < _correlationDepth; j++) {
				temp3->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				temp3->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) = 0.5f * temp3->GetAt(i);
		}
		
		mprintf("    Second anisotropic part...\n");
		for(i = 0; i < _correlationDepth; i++) {
			temp3->GetAt(i) = 0.0f;
		}
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			for(j = 0; j < n; j++) {
				temp2->GetAt(j) = ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount + m_iShowMolCount + i])->GetAt(j) - ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * 2 + m_iShowMolCount * 2 + i])->GetAt(j);
			}
			ac->AutoCorrelate(temp2, temp);
			for(j = 0; j < _correlationDepth; j++) {
				temp3->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				temp3->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) += 0.5f * temp3->GetAt(i);
		}
		
		mprintf("    Third anisotropic part...\n");
		for(i = 0; i < _correlationDepth; i++) {
			temp3->GetAt(i) = 0.0f;
		}
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			for(j = 0; j < n; j++) {
				temp2->GetAt(j) = ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * 2 + m_iShowMolCount * 2 + i])->GetAt(j) - ((CxFloatArray *)_polarizabilityCache[i])->GetAt(j);
			}
			ac->AutoCorrelate(temp2, temp);
			for(j = 0; j < _correlationDepth; j++) {
				temp3->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				temp3->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) += 0.5f * temp3->GetAt(i);
		}
		
		mprintf("    Fourth anisotropic part...\n");
		for(i = 0; i < _correlationDepth; i++) {
			temp3->GetAt(i) = 0.0f;
		}
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			for(j = 0; j < n; j++) {
				temp2->GetAt(j) = 0.5f * (((CxFloatArray *)_polarizabilityCache[m_iShowMolCount + i])->GetAt(j) + ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount + i])->GetAt(j));
			}
			ac->AutoCorrelate(temp2, temp);
			for(j = 0; j < _correlationDepth; j++) {
				temp3->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				temp3->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) += 3.0f * temp3->GetAt(i);
		}
		
		mprintf("    Fifth anisotropic part...\n");
		for(i = 0; i < _correlationDepth; i++) {
			temp3->GetAt(i) = 0.0f;
		}
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			for(j = 0; j < n; j++) {
				temp2->GetAt(j) = 0.5f * (((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount + m_iShowMolCount * 2 + i])->GetAt(j) + ((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * 2 + m_iShowMolCount + i])->GetAt(j));
			}
			ac->AutoCorrelate(temp2, temp);
			for(j = 0; j < _correlationDepth; j++) {
				temp3->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				temp3->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) += 3.0f * temp3->GetAt(i);
		}
		
		mprintf("    Sixth anisotropic part...\n");
		for(i = 0; i < _correlationDepth; i++) {
			temp3->GetAt(i) = 0.0f;
		}
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			for(j = 0; j < n; j++) {
				temp2->GetAt(j) = 0.5f * (((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * 2 + i])->GetAt(j) + ((CxFloatArray *)_polarizabilityCache[m_iShowMolCount * 2 + i])->GetAt(j));
			}
			ac->AutoCorrelate(temp2, temp);
			for(j = 0; j < _correlationDepth; j++) {
				temp3->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				temp3->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) += 3.0f * temp3->GetAt(i);
		}
		
		if(_saveACF) {
			if (g_ramanCompat)
				filename.sprintf("%s/acf_aniso_%s.csv", g_ramanDir, _name);
			else
				filename.sprintf("raman_acf_aniso_%s.csv", _name);
			mprintf("    Saving autocorrelation function as %s...\n", (const char *)filename);
			FILE *acf_file = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(acf_file, "%.2f; %.10G\n", i * g_fTimestepLength * g_stride, acf->GetAt(i));
			}
			fclose(acf_file);
		}
		
		temp->CopyFrom(acf);
		
		if(_windowFunction == 1) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= powf(cosf((float)i / (temp->GetSize() - 1) / 2.0f * Pi), 2.0f);
			}
		} else if(_windowFunction == 2) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= expf(-(float)i / _windowFunctionParameter);
			}
		} else if(_windowFunction == 3) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= expf(-(float)i * i / _windowFunctionParameter / _windowFunctionParameter);
			}
		} else if(_windowFunction != 0) {
			eprintf("Unknown window function.\n");
			abort();
		}
		
		if(_saveACF) {
			if (g_ramanCompat)
				filename.sprintf("%s/acf_windowed_aniso_%s.csv", g_ramanDir, _name);
			else
				filename.sprintf("raman_acf_windowed_aniso_%s.csv", _name);
			mprintf("    Saving windowed autocorrelation function as %s...\n", (const char *)filename);
			FILE *acf_file = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(acf_file, "%.2f; %.10G\n", i * g_fTimestepLength * g_stride, acf->GetAt(i));
			}
			fclose(acf_file);
		}
		
		if(_zeroPadding > 0) {
			for(i = 0; i < _zeroPadding; i++) {
				temp->Add(0.0f);
			}
		}
		
		oldSize = temp->GetSize();
		temp->SetSize(2 * oldSize);
		for(i = 1; i < oldSize; i++) {
			temp->GetAt(oldSize + i) = temp->GetAt(oldSize - i);
		}
		temp->GetAt(oldSize) = 0.0f;
		
		mprintf("    Performing Fourier transformation...\n");
		fft->PrepareFFT_C2C(temp->GetSize());
		for(i = 0; i < temp->GetSize(); i++) {
			fft->m_pInput[2*i] = temp->GetAt(i);
			fft->m_pInput[2*i+1] = 0.0f;
		}
		fft->DoFFT();
		
		for(i = 0; i < _specSize; i++) {
			float freq = i * _specResolution;
			spectrum2->GetAt(i) = 4.160440e-18f * powf(_laserFreq - freq, 4.0f) / freq / (1 - expf(-1.438777f * freq / _temperature)) * fft->m_pOutput[2*i] * g_fTimestepLength * g_stride; // Output in 1e-30*m^2*cm
		}
		spectrum2->GetAt(0) = 0.0f;
		
		if(_finiteDifferenceCorrection) {
			float f = _specResolution * g_fTimestepLength * g_stride * 1.883652e-4f;
			for(i = 1; i < _specSize; i++) {
				spectrum2->GetAt(i) *= powf(f * i / sinf(f * i), 2.0f); // Divide by sinc function to correct finite difference derivation
			}
		}
		
		if(_quantumCorrection != 1) {
			for(i = 0; i < _specSize; i++) {
				float factor = 1.0f;
				if(_quantumCorrection == 2) {
					factor = 2.0f / (1.0f + expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature)) / (1.438777f * _specResolution * i / _quantumCorrectionTemperature) * (1.0f - expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature));
				} else if(_quantumCorrection == 3) {
					factor = expf(0.719388f * _specResolution * i / _quantumCorrectionTemperature) / (1.438777f * _specResolution * i / _quantumCorrectionTemperature) * (1.0f - expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature));
				}
				spectrum2->GetAt(i) *= factor;
			}
		}
		
		if (g_ramanCompat)
			filename.sprintf("%s/spectrum_para_%s.csv", g_ramanDir, _name);
		else
			filename.sprintf("raman_spectrum_para_%s.csv", _name);
		mprintf("    Saving parallel spectrum as %s...\n", (const char *)filename);
		FILE *specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (10^-30*K*m^2*cm); Integral (10^-30*K*m^2)\n");

		double integral = 0.0;
		for(i = 0; i < _specSize; i++) {
			integral += ((double)spectrum1->GetAt(i) + 4.0 / 45.0 * (double)spectrum2->GetAt(i)) * _specResolution;
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum1->GetAt(i) + 4.0f / 45.0f * spectrum2->GetAt(i), integral);
		}
		fclose(specFile);
		
		if (g_ramanCompat)
			filename.sprintf("%s/spectrum_ortho_%s.csv", g_ramanDir, _name);
		else
			filename.sprintf("raman_spectrum_ortho_%s.csv", _name);
		mprintf("    Saving orthogonal spectrum as %s...\n", (const char *)filename);
		specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (10^-30*K*m^2*cm); Integral (10^-30*K*m^2)\n");

		integral = 0.0;
		for(i = 0; i < _specSize; i++) {
			integral += (double)spectrum2->GetAt(i) / 15.0 * _specResolution;
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum2->GetAt(i) / 15.0f, integral);
		}
		fclose(specFile);
		
		if (g_ramanCompat)
			filename.sprintf("%s/spectrum_unpol_%s.csv", g_ramanDir, _name);
		else
			filename.sprintf("raman_spectrum_unpol_%s.csv", _name);
		mprintf("    Saving unpolarized spectrum as %s...\n", (const char *)filename);
		specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (10^-30*K*m^2*cm); Integral (10^-30*K*m^2)\n");

		integral = 0.0;
		for(i = 0; i < _specSize; i++) {
			integral += ((double)spectrum1->GetAt(i) + 7.0 / 45.0 * (double)spectrum2->GetAt(i)) * _specResolution;
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum1->GetAt(i) + 7.0f / 45.0f * spectrum2->GetAt(i), integral);
		}
		fclose(specFile);
		
		if (g_ramanCompat)
			filename.sprintf("%s/depol_ratio_%s.csv", g_ramanDir, _name);
		else
			filename.sprintf("raman_depol_ratio_%s.csv", _name);
		mprintf("    Saving depolarization ratio as %s...\n", (const char *)filename);
		specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Depolarization ratio\n");

		for(i = 0; i < _specSize; i++) {
				fprintf(specFile, "%.2f; %.8G\n", _specResolution * i, (spectrum2->GetAt(i) / 15.0f) / (spectrum1->GetAt(i) + 4.0f / 45.0f * spectrum2->GetAt(i)));
		}
		fclose(specFile);
	} else {
		mprintf("    Parallel part...\n");
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			if (g_ramanCompat) {
				ac->AutoCorrelate((CxFloatArray *)_polarizabilityCache[i], temp);
			} else {
				ac->AutoCorrelate((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * (m_fieldMode - 2) + m_iShowMolCount * (m_fieldMode - 2) + i], temp);
			}
			for(j = 0; j < _correlationDepth; j++) {
				acf->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				acf->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		
		if(_saveACF) {
			if (g_ramanCompat)
				filename.sprintf("%s/acf_para_%s.csv", g_ramanDir, _name);
			else
				filename.sprintf("raman_acf_para_%s.csv", _name);
			mprintf("    Saving autocorrelation function as %s...\n", (const char *)filename);
			FILE *acf_file = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(acf_file, "%.2f; %.10G\n", i * g_fTimestepLength * g_stride, acf->GetAt(i));
			}
			fclose(acf_file);
		}
		
		temp->CopyFrom(acf);
		
		if(_windowFunction == 1) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= powf(cosf((float)i / (temp->GetSize() - 1) / 2.0f * Pi), 2.0f);
			}
		} else if(_windowFunction == 2) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= expf(-(float)i / _windowFunctionParameter);
			}
		} else if(_windowFunction == 3) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= expf(-(float)i * i / _windowFunctionParameter / _windowFunctionParameter);
			}
		} else if(_windowFunction != 0) {
			eprintf("Unknown window function.\n");
			abort();
		}
		
		if(_saveACF) {
			if (g_ramanCompat)
				filename.sprintf("%s/acf_windowed_para_%s.csv", g_ramanDir, _name);
			else
				filename.sprintf("raman_acf_windowed_para_%s.csv", _name);
			mprintf("    Saving windowed autocorrelation function as %s...\n", (const char *)filename);
			FILE *acf_file = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(acf_file, "%.2f; %.10G\n", i * g_fTimestepLength * g_stride, acf->GetAt(i));
			}
			fclose(acf_file);
		}
		
		if(_zeroPadding > 0) {
			for(i = 0; i < _zeroPadding; i++) {
				temp->Add(0.0f);
			}
		}
		
		int oldSize = temp->GetSize();
		temp->SetSize(2 * oldSize);
		for(i = 1; i < oldSize; i++) {
			temp->GetAt(oldSize + i) = temp->GetAt(oldSize - i);
		}
		temp->GetAt(oldSize) = 0.0f;
		
		mprintf("    Performing Fourier transformation...\n");
		fft->PrepareFFT_C2C(temp->GetSize());
		for(i = 0; i < temp->GetSize(); i++) {
			fft->m_pInput[2*i] = temp->GetAt(i);
			fft->m_pInput[2*i+1] = 0.0f;
		}
		fft->DoFFT();
		
		for(i = 0; i < _specSize; i++) {
			float freq = i * _specResolution;
			spectrum1->GetAt(i) = 4.160440e-18f * powf(_laserFreq - freq, 4.0f) / freq / (1 - expf(-1.438777f * freq / _temperature)) * fft->m_pOutput[2*i] * g_fTimestepLength * g_stride; // Output in 1e-30*m^2*cm
		}
		spectrum1->GetAt(0) = 0.0f;
		
		if(_finiteDifferenceCorrection) {
			float f = _specResolution * g_fTimestepLength * g_stride * 1.883652e-4f;
			for(i = 1; i < _specSize; i++) {
				spectrum1->GetAt(i) *= powf(f * i / sinf(f * i), 2.0f); // Divide by sinc function to correct finite difference derivation
			}
		}
		
		if(_quantumCorrection != 1) {
			for(i = 0; i < _specSize; i++) {
				float factor = 1.0f;
				if(_quantumCorrection == 2) {
					factor = 2.0f / (1.0f + expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature)) / (1.438777f * _specResolution * i / _quantumCorrectionTemperature) * (1.0f - expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature));
				} else if(_quantumCorrection == 3) {
					factor = expf(0.719388f * _specResolution * i / _quantumCorrectionTemperature) / (1.438777f * _specResolution * i / _quantumCorrectionTemperature) * (1.0f - expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature));
				}
				spectrum1->GetAt(i) *= factor;
			}
		}
		
		if (g_ramanCompat)
			filename.sprintf("%s/spectrum_para_%s.csv", g_ramanDir, _name);
		else
			filename.sprintf("raman_spectrum_para_%s.csv", _name);
		mprintf("    Saving parallel spectrum as %s...\n", (const char *)filename);
		FILE *specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (10^-30*K*m^2*cm); Integral (10^-30*K*m^2)\n");

		double integral = 0.0;
		for(i = 0; i < _specSize; i++) {
			integral += (double)spectrum1->GetAt(i) * _specResolution;
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum1->GetAt(i), integral);
		}
		fclose(specFile);
		
		mprintf("    Orthogonal part...\n");
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) = 0.0f;
		}
		mprintf(WHITE, "     [");
		for(i = 0; i < m_iShowMolCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			if (g_ramanCompat) {
				ac->AutoCorrelate((CxFloatArray *)_polarizabilityCache[m_iShowMolCount + i], temp);
			} else {
				ac->AutoCorrelate((CxFloatArray *)_polarizabilityCache[3*m_iShowMolCount * ((m_fieldMode - 1) % 3) + m_iShowMolCount * (m_fieldMode - 2) + i], temp);
			}
			for(j = 0; j < _correlationDepth; j++) {
				acf->GetAt(j) += temp->GetAt(j);
			}
		}
		mprintf(WHITE, "]\n");
		if(m_iShowMol != -1) {
			for(i = 0; i < _correlationDepth; i++) {
				acf->GetAt(i) /= (float)m_iShowMolCount;
			}
		}
		
		if(_saveACF) {
			if (g_ramanCompat)
				filename.sprintf("%s/acf_ortho_%s.csv", g_ramanDir, _name);
			else
				filename.sprintf("raman_acf_ortho_%s.csv", _name);
			mprintf("    Saving autocorrelation function as %s...\n", (const char *)filename);
			FILE *acf_file = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(acf_file, "%.2f; %.10G\n", i * g_fTimestepLength * g_stride, acf->GetAt(i));
			}
			fclose(acf_file);
		}
		
		temp->CopyFrom(acf);
		
		if(_windowFunction == 1) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= powf(cosf((float)i / (temp->GetSize() - 1) / 2.0f * Pi), 2.0f);
			}
		} else if(_windowFunction == 2) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= expf(-(float)i / _windowFunctionParameter);
			}
		} else if(_windowFunction == 3) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= expf(-(float)i * i / _windowFunctionParameter / _windowFunctionParameter);
			}
		} else if(_windowFunction != 0) {
			eprintf("Unknown window function.\n");
			abort();
		}
		
		if(_saveACF) {
			if (g_ramanCompat)
				filename.sprintf("%s/acf_windowed_ortho_%s.csv", g_ramanDir, _name);
			else
				filename.sprintf("raman_acf_windowed_ortho_%s.csv", _name);
			mprintf("    Saving windowed autocorrelation function as %s...\n", (const char *)filename);
			FILE *acf_file = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(acf_file, "%.2f; %.10G\n", i * g_fTimestepLength * g_stride, acf->GetAt(i));
			}
			fclose(acf_file);
		}
		
		if(_zeroPadding > 0) {
			for(i = 0; i < _zeroPadding; i++) {
				temp->Add(0.0f);
			}
		}
		
		oldSize = temp->GetSize();
		temp->SetSize(2 * oldSize);
		for(i = 1; i < oldSize; i++) {
			temp->GetAt(oldSize + i) = temp->GetAt(oldSize - i);
		}
		temp->GetAt(oldSize) = 0.0f;
		
		mprintf("    Performing Fourier transformation...\n");
		fft->PrepareFFT_C2C(temp->GetSize());
		for(i = 0; i < temp->GetSize(); i++) {
			fft->m_pInput[2*i] = temp->GetAt(i);
			fft->m_pInput[2*i+1] = 0.0f;
		}
		fft->DoFFT();
		
		for(i = 0; i < _specSize; i++) {
			float freq = i * _specResolution;
			spectrum2->GetAt(i) = 4.160440e-18f * powf(_laserFreq - freq, 4.0f) / freq / (1 - expf(-1.438777f * freq / _temperature)) * fft->m_pOutput[2*i] * g_fTimestepLength * g_stride; // Output in 1e-30*m^2*cm
		}
		spectrum2->GetAt(0) = 0.0f;
		
		if(_finiteDifferenceCorrection) {
			float f = _specResolution * g_fTimestepLength * g_stride * 1.883652e-4f;
			for(i = 1; i < _specSize; i++) {
				spectrum2->GetAt(i) *= powf(f * i / sinf(f * i), 2.0f); // Divide by sinc function to correct finite difference derivation
			}
		}
		
		if(_quantumCorrection != 1) {
			for(i = 0; i < _specSize; i++) {
				float factor = 1.0f;
				if(_quantumCorrection == 2) {
					factor = 2.0f / (1.0f + expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature)) / (1.438777f * _specResolution * i / _quantumCorrectionTemperature) * (1.0f - expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature));
				} else if(_quantumCorrection == 3) {
					factor = expf(0.719388f * _specResolution * i / _quantumCorrectionTemperature) / (1.438777f * _specResolution * i / _quantumCorrectionTemperature) * (1.0f - expf(-1.438777f * _specResolution * i / _quantumCorrectionTemperature));
				}
				spectrum2->GetAt(i) *= factor;
			}
		}
		
		if (g_ramanCompat)
			filename.sprintf("%s/spectrum_ortho_%s.csv", g_ramanDir, _name);
		else
			filename.sprintf("raman_spectrum_ortho_%s.csv", _name);
		mprintf("    Saving orthogonal spectrum as %s...\n", (const char *)filename);
		specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (10^-30*K*m^2*cm); Integral (10^-30*K*m^2)\n");

		integral = 0.0;
		for(i = 0; i < _specSize; i++) {
			integral += (double)spectrum2->GetAt(i) * _specResolution;
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum2->GetAt(i), integral);
		}
		fclose(specFile);
		
		if (g_ramanCompat)
			filename.sprintf("%s/spectrum_unpol_%s.csv", g_ramanDir, _name);
		else
			filename.sprintf("raman_spectrum_unpol_%s.csv", _name);
		mprintf("    Saving unpolarized spectrum as %s...\n", (const char *)filename);
		specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (10^-30*K*m^2*cm); Integral (10^-30*K*m^2)\n");

		integral = 0.0;
		for(i = 0; i < _specSize; i++) {
			integral += ((double)spectrum1->GetAt(i) + (double)spectrum2->GetAt(i)) * _specResolution;
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum1->GetAt(i) + spectrum2->GetAt(i), integral);
		}
		fclose(specFile);
		
		if (g_ramanCompat)
			filename.sprintf("%s/depol_ratio_%s.csv", g_ramanDir, _name);
		else
			filename.sprintf("raman_depol_ratio_%s.csv", _name);
		mprintf("    Saving depolarization ratio as %s...\n", (const char *)filename);
		specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Depolarization ratio\n");

		for(i = 0; i < _specSize; i++) {
				fprintf(specFile, "%.2f; %.8G\n", _specResolution * i, spectrum2->GetAt(i) / spectrum1->GetAt(i));
		}
		fclose(specFile);
	}
	
	delete ac;
	delete temp;
	delete temp2;
	delete temp3;
	delete fft;
	delete spectrum1;
	delete spectrum2;
	
// 	_ramanDyn->finalize();
}


#ifdef TARGET_LINUX
static bool parseSettings(FILE *settingsFile) {
	char buf[BUF_SIZE];
	
	if(fgets(buf, BUF_SIZE, settingsFile) == NULL)
		return false;
	int temp = -1;
	if(sscanf(buf, "%d", &temp) < 1)
		return false;
	if(temp == 0)
		g_orientAvg = false;
	else if(temp == 1)
		g_orientAvg = true;
	else
		return false;
	
	if(fgets(buf, BUF_SIZE, settingsFile) == NULL)
		return false;
	if(sscanf(buf, "%f", &g_fieldStrength) < 1)
		return false;
	
	if(fgets(buf, BUF_SIZE, settingsFile) == NULL)
		return false;
	if(sscanf(buf, "%d", &g_stride) < 1)
		return false;
	
	return true;
}
#else
static bool parseSettings(FILE *settingsFile) {
	eprintf("parseSettings(): No Windows implementation available.\n");
	return false;
}
#endif


bool gatherRaman() {
	g_ramanCompat = AskYesNo("    Use Raman compatibility mode (y/n)? [no] ", false);
	
	if (g_ramanCompat) {
		CxString buf;
#ifndef TARGET_LINUX
		mprintf(RED, "Raman calculations are currently possible only under Linux.\n");
		return false;
#else
		mprintf("    To calculate Raman spectra, polarizabilities have to be determined.\n");
		mprintf("    TRAVIS creates CP2K input files for numerical polarizabilities\n    and will process the resulting data in a second run.\n\n");
		
		g_newRaman = AskYesNo("    Do you wish to create new CP2K input files (y) or process existing results (n)? [y] ", true);
	//	char buf[BUF_SIZE];
		if(g_newRaman) {
			AskString("    Please enter a name for the directory to collect the data: [raman] ", &buf, "raman");
		} else {
			AskString("    Please enter the name of the directory to take the data from: [raman] ", &buf, "raman");
		}
		try { g_ramanDir = new char[strlen(buf)+1]; } catch(...) { g_ramanDir = NULL; }
		if(g_ramanDir == NULL) NewException((double)(strlen(buf)+1)*sizeof(char), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		strcpy(g_ramanDir, buf);
		
		if(g_newRaman) {
			g_bKeepOriginalCoords = true;
			if(FileExist(g_ramanDir)) {
				mprintf(RED, "A file or a directory \"%s\" already exists. Please remove it first.\n", g_ramanDir);
				return false;
			}
			if(mkdir(g_ramanDir, S_IRWXU) != 0) {
				mprintf(RED, "Directory \"%s\" could not be created: %s\n", g_ramanDir, strerror(errno));
				return false;
			}
			
			char filename[BUF_SIZE];
			snprintf(filename, BUF_SIZE, "%s/settings.dat", g_ramanDir);
			FILE *settingsFile = fopen(filename, "w");
			if(settingsFile == NULL) {
				mprintf(RED, "Could not open settings file \"%s\": %s\n", filename, strerror(errno));
				return false;
			}
			
			g_orientAvg = AskYesNo("\n    Use orientational averaging? [no] ", false);
			if(g_orientAvg)
				fprintf(settingsFile, "%d\n", 1);
			else
				fprintf(settingsFile, "%d\n", 0);
			
			g_fieldStrength = AskFloat("    Field strength in atomic units [5.0e-4] ", 5.0e-4);
			fprintf(settingsFile, "%.6E\n", g_fieldStrength);
			
			g_stride = AskInteger("    Calculate polarizability for every n-th timestep [1] ", 1);
			fprintf(settingsFile, "%d\n", g_stride);
			
			fclose(settingsFile);
			
			snprintf(filename, BUF_SIZE, "%s/1", g_ramanDir);
			if(mkdir(filename, S_IRWXU) != 0) {
				mprintf(RED, "Directory \"%s\" could not be created: %s\n", filename, strerror(errno));
				return false;
			}
			if(g_orientAvg) {
				snprintf(filename, BUF_SIZE, "%s/2", g_ramanDir);
				if(mkdir(filename, S_IRWXU) != 0) {
					mprintf(RED, "Directory \"%s\" could not be created: %s\n", filename, strerror(errno));
					return false;
				}
				snprintf(filename, BUF_SIZE, "%s/3", g_ramanDir);
				if(mkdir(filename, S_IRWXU) != 0) {
					mprintf(RED, "Directory \"%s\" could not be created: %s\n", filename, strerror(errno));
					return false;
				}
			}
			
			FILE* templateFile;
			snprintf(filename, BUF_SIZE, "%s/template.inp", g_ramanDir);
			templateFile = fopen(filename, "w");
			if(templateFile == NULL) {
				mprintf(RED, "Could not open template file \"%s\": %s\n", filename, strerror(errno));
				return false;
			}
			fprintf(templateFile, "&GLOBAL\n");
			fprintf(templateFile, " PROJECT_NAME polarizability\n");
			fprintf(templateFile, " RUN_TYPE MD\n");
			fprintf(templateFile, " PRINT_LEVEL LOW\n");
			fprintf(templateFile, "&END\n");
			fprintf(templateFile, "&FORCE_EVAL\n");
			fprintf(templateFile, " &DFT\n");
			fprintf(templateFile, "  BASIS_SET_FILE_NAME BASIS_MOLOPT\n");
			fprintf(templateFile, "  POTENTIAL_FILE_NAME POTENTIAL\n");
			fprintf(templateFile, "  &MGRID\n");
			fprintf(templateFile, "   NGRIDS 5\n");
			fprintf(templateFile, "   CUTOFF 280\n");
			fprintf(templateFile, "   REL_CUTOFF 40\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &SCF\n");
			fprintf(templateFile, "   SCF_GUESS ATOMIC\n");
			fprintf(templateFile, "   MAX_SCF 200\n");
			fprintf(templateFile, "   EPS_SCF 1.0E-5\n");
			fprintf(templateFile, "   &OT\n");
			fprintf(templateFile, "    MINIMIZER DIIS\n");
			fprintf(templateFile, "    PRECONDITIONER FULL_SINGLE_INVERSE\n");
			fprintf(templateFile, "   &END\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &XC\n");
			fprintf(templateFile, "   &XC_FUNCTIONAL BLYP\n");
			fprintf(templateFile, "   &END\n");
			fprintf(templateFile, "   &XC_GRID\n");
			fprintf(templateFile, "    XC_SMOOTH_RHO NN10\n");
			fprintf(templateFile, "    XC_DERIV NN10_SMOOTH\n");
			fprintf(templateFile, "   &END\n");
			fprintf(templateFile, "   &VDW_POTENTIAL\n");
			fprintf(templateFile, "    POTENTIAL_TYPE PAIR_POTENTIAL\n");
			fprintf(templateFile, "    &PAIR_POTENTIAL\n");
			fprintf(templateFile, "     TYPE DFTD3\n");
			fprintf(templateFile, "     REFERENCE_FUNCTIONAL BLYP\n");
			fprintf(templateFile, "     PARAMETER_FILE_NAME dftd3.dat\n");
			fprintf(templateFile, "    &END\n");
			fprintf(templateFile, "   &END\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &LOCALIZE\n");
			fprintf(templateFile, "   METHOD CRAZY\n");
			fprintf(templateFile, "   MAX_ITER 2000\n");
			fprintf(templateFile, "   &PRINT\n");
			fprintf(templateFile, "    &WANNIER_CENTERS\n");
			fprintf(templateFile, "     IONS+CENTERS\n");
			fprintf(templateFile, "     FILENAME =polarizability_wannier.xyz\n");
			fprintf(templateFile, "    &END\n");
			fprintf(templateFile, "   &END\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &PERIODIC_EFIELD\n");
			fprintf(templateFile, "   INTENSITY ###!field strength will be put here###\n");
			fprintf(templateFile, "   POLARISATION ###!polarisation vector will be put here###\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, " &END\n");
			fprintf(templateFile, " &SUBSYS\n");
			fprintf(templateFile, "  &CELL\n");
			fprintf(templateFile, "   ABC %.6f %.6f %.6f\n", g_fBoxX / 100.0f, g_fBoxY / 100.0f, g_fBoxZ / 100.0f);
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &COORD\n");
			fprintf(templateFile, "###!coordinates will be put here###\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &KIND H\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q1\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &KIND C\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q4\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &KIND N\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q5\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, "  &KIND O\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q6\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, " &END\n");
			fprintf(templateFile, "&END\n");
			fprintf(templateFile, "&MOTION\n");
			fprintf(templateFile, " &MD\n");
			fprintf(templateFile, "  ENSEMBLE REFTRAJ\n");
			fprintf(templateFile, "  STEPS ###!number of steps will be put here###\n");
			fprintf(templateFile, "  TIMESTEP %f\n", g_fTimestepLength * g_stride);
			fprintf(templateFile, "  &REFTRAJ\n");
			fprintf(templateFile, "   EVAL_ENERGY_FORCES\n");
			fprintf(templateFile, "   TRAJ_FILE_NAME ../polarizability_reftraj.xyz\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, " &END\n");
			fprintf(templateFile, " &PRINT\n");
			fprintf(templateFile, "  &RESTART\n");
			fprintf(templateFile, "   &EACH\n");
			fprintf(templateFile, "    MD 1\n");
			fprintf(templateFile, "   &END\n");
			fprintf(templateFile, "  &END\n");
			fprintf(templateFile, " &END\n");
			fprintf(templateFile, "&END\n");
			fclose(templateFile);
			
			mprintf("\n    An input template for the polarizability calculations has been written to \"%s/template.inp\"\n    Please modify it according to your needs.\n    Press any key when you are finished.\n", g_ramanDir);
			getchar();
			
			if(g_orientAvg) {
				mprintf("    CP2K input files will be created in \"%s/1\", \"%s/2\", and \"%s/3\".\n", g_ramanDir, g_ramanDir, g_ramanDir);
				mprintf("    After TRAVIS has finished, please run them and make sure that the resulting trajectories \"polarizability_wannier.xyz\"\n    are placed in the same directories before you execute TRAVIS for evaluation.\n\n");
			} else {
				mprintf("    A CP2K input file will be created in \"%s/1\".\n", g_ramanDir);
				mprintf("    After TRAVIS has finished, please run it and make sure that the resulting trajectory \"polarizability_wannier.xyz\"\n    is placed in the same directory before you run TRAVIS for evaluation.\n\n");
			}
		} else {
			if(!FileExist(g_ramanDir)) {
				mprintf(RED, "The directory \"%s\" was not found.\n", g_ramanDir);
				return false;
			}
			
			char filename[BUF_SIZE];
			snprintf(filename, BUF_SIZE, "%s/settings.dat", g_ramanDir);
			FILE *settingsFile;
			settingsFile = fopen(filename, "r");
			if(settingsFile == NULL) {
				mprintf(RED, "Could not open settings file \"%s\": %s\n", filename, strerror(errno));
				return false;
			}
			if(!parseSettings(settingsFile)) {
				mprintf(RED, "Could not parse settings file.\n");
				return false;
			}
			
			mprintf("\n");
			if(g_orientAvg)
				mprintf("    Using orientational averaging\n");
			else
				mprintf("    Not using orientational averaging\n");
			mprintf("    Field strength: %.6e a. u.\n", g_fieldStrength);
			mprintf("    Using every %d%s timestep\n\n", g_stride, (g_stride == 1) ? "st" : ((g_stride == 2) ? "nd" : ((g_stride == 3) ? "rd" : "th")));
			
			mprintf("    The scattering cross sections are calculated for an\n    x polarized incident laser beam along the y axis with detection in z direction.\n");
			
			g_bDipole = true;
			ParseDipole();
			
			while(true) {
				mprintf(YELLOW, "\n>>> Raman Observation %d >>>\n\n", g_ramObserv.GetSize() + 1);
				
				CRamanObservation *obs;
				try { obs = new CRamanObservation(); } catch(...) { obs = NULL; }
				if(obs == NULL) NewException((double)sizeof(CRamanObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				g_ramObserv.Add(obs);
				
				mprintf(YELLOW, "<<< End of Raman Observation %d <<<\n\n", g_ramObserv.GetSize());
				
				if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
					break;
				mprintf("\n");
			}
			
			if(AskYesNo("    Compute Raman spectrum of whole system (y/n)? [no] ", false)) {
				mprintf(YELLOW, "\n>>> Global Raman Observation >>>\n\n");
				
				CRamanObservation *obs;
				try { obs = new CRamanObservation(true); } catch(...) { obs = NULL; }
				if(obs == NULL) NewException((double)sizeof(CRamanObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				g_ramObserv.Add(obs);
				
				mprintf(YELLOW, "<<< End of Global Raman Observation <<<\n\n");
			}
		}
#endif
	} else {
		g_stride = 1;
		
		parsePolarizability();
		
		while (true) {
			mprintf(YELLOW, "\n>>> Raman Observation %d >>>\n\n", g_ramObserv.GetSize() + 1);
			
			CRamanObservation *obs;
			try { obs = new CRamanObservation(); } catch(...) { obs = NULL; }
			if (obs == NULL) NewException((double)sizeof(CRamanObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			g_ramObserv.Add(obs);
			
			mprintf(YELLOW, "<<< End of Raman Observation %d <<<\n\n", g_ramObserv.GetSize());
			
			if (!AskYesNo("    Add another observation (y/n)? [no] ", false))
				break;
			mprintf("\n");
		}
		
		if (AskYesNo("    Compute Raman spectrum of whole system (y/n)? [no] ", false)) {
			mprintf(YELLOW, "\n>>> Global Raman Observation >>>\n\n");
			
			CRamanObservation *obs;
			try { obs = new CRamanObservation(true); } catch(...) { obs = NULL; }
			if (obs == NULL) NewException((double)sizeof(CRamanObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			g_ramObserv.Add(obs);
			
			mprintf(YELLOW, "<<< End of Global Raman Observation <<<\n\n");
		}
	}
	
	return true;
}

bool initializeRaman() {
	int i;
	
	if (g_ramanCompat) {
		if(g_newRaman) {
			char filename[BUF_SIZE];
#ifdef TARGET_LINUX
			snprintf(filename, BUF_SIZE, "%s/template.inp", g_ramanDir);
#else
			sprintf(filename, "%s/template.inp", g_ramanDir);
#endif
			FILE *templateFile;
			templateFile = fopen(filename, "r");
			if(templateFile == NULL) {
				mprintf(RED, "Could not open template file \"%s\": %s\n", filename, strerror(errno));
				return false;
			}
			fseek(templateFile, 0L, SEEK_END);
			long length = ftell(templateFile);
			if(length < 0) {
				mprintf(RED, "Could not determine size of template file: %s\n", strerror(errno));
				fclose(templateFile);
				return false;
			}
			
			try { g_inputTemplate = new char[length + 1]; } catch(...) { g_inputTemplate = NULL; }
			if(g_inputTemplate == NULL) NewException((double)(length + 1)*sizeof(char), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			
			rewind(templateFile);
			if((long)fread(g_inputTemplate, sizeof(char), length, templateFile) < length) {
				mprintf(RED, "Could not read template file: %s\n", strerror(errno));
				fclose(templateFile);
				return false;
			}
			g_inputTemplate[length] = '\0';
			
			fclose(templateFile);
			
			g_templateFieldPos = strstr(g_inputTemplate, "###!field strength will be put here###");
			if(g_templateFieldPos == NULL) {
				mprintf(RED, "Position mark for field strength missing in template.\n");
				fclose(templateFile);
				return false;
			}
			g_templatePolPos = strstr(g_inputTemplate, "###!polarisation vector will be put here###");
			if(g_templatePolPos == NULL) {
				mprintf(RED, "Position mark for polarisation vector missing in template.\n");
				fclose(templateFile);
				return false;
			}
			g_templateCoordPos = strstr(g_inputTemplate, "###!coordinates will be put here###");
			if(g_templateCoordPos == NULL) {
				mprintf(RED, "Position mark for coordinates missing in template.\n");
				fclose(templateFile);
				return false;
			}
			g_templateStepsPos = strstr(g_inputTemplate, "###!number of steps will be put here###");
			if(g_templateStepsPos == NULL) {
				mprintf(RED, "Position mark for number of steps missing in template.\n");
				fclose(templateFile);
				return false;
			}
			g_templateFieldPos[0] = '\0';
			g_templatePolPos[0] = '\0';
			g_templateCoordPos[0] = '\0';
			g_templateStepsPos[0] = '\0';
			
#ifdef TARGET_LINUX
			snprintf(filename, BUF_SIZE, "%s/polarizability_reftraj.xyz", g_ramanDir);
#else
			sprintf(filename, "%s/polarizability_reftraj.xyz", g_ramanDir);
#endif
			g_reftrajFile = fopen(filename, "w");
			if(g_reftrajFile == NULL) {
				mprintf(RED, "Could not open reference trajectory \"%s\": %s\n", filename, strerror(errno));
				fclose(templateFile);
				return false;
			}
		} else {
			for(i = 0; i < g_ramObserv.GetSize(); i++) {
				mprintf("Initializing Raman Observation %d...\n", i+1);
				CRamanObservation *obs = (CRamanObservation *)g_ramObserv[i];
				obs->initialize();
			}
			
			char filename[BUF_SIZE];
			for(i = 0; i < (g_orientAvg ? 3 : 1); i++) {
				if(g_iTrajFormat == 5) {
#ifdef TARGET_LINUX
					snprintf(filename, BUF_SIZE, "%s/%d/polarizability_density.cube", g_ramanDir, i + 1);
#else
					sprintf(filename, "%s/%d/polarizability_density.cube", g_ramanDir, i + 1);
#endif
				} else {
#ifdef TARGET_LINUX
					snprintf(filename, BUF_SIZE, "%s/%d/polarizability_wannier.xyz", g_ramanDir, i + 1);
#else
					sprintf(filename, "%s/%d/polarizability_wannier.xyz", g_ramanDir, i + 1);
#endif
				}
				g_polFile[i] = fopen(filename, "r");
				if(g_polFile[i] == NULL) {
					mprintf(RED, "Could not open trajectory \"%s\": %s\n", filename, strerror(errno));
					return false;
				}
			}
			
			for(i = 0; i < (g_orientAvg ? 3 : 1); i++) {
				try { g_timestep[i] = new CTimeStep(); } catch(...) { g_timestep[i] = NULL; }
				if(g_timestep[i] == NULL) NewException((double)sizeof(CTimeStep), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
		}
	} else {
		for (i = 0; i < g_ramObserv.GetSize(); i++) {
			mprintf("Initializing Raman Observation %d...\n", i + 1);
			CRamanObservation *obs = (CRamanObservation *)g_ramObserv[i];
			obs->initialize();
		}
	}
	
	return true;
}

void processRaman(CTimeStep *ts) {
	int i, j;
	if (g_ramanCompat) {
		g_step++;
		if(g_newRaman) {
			if(g_step % g_stride == 0) {
				g_steps++;
				int numAtoms = 0;
				for(i = 0; i < g_iGesAtomCount; i++)
					if(g_waAtomMolIndex[i] != 60000)
						numAtoms++;
				fprintf(g_reftrajFile, "%d\nStep %d\n", numAtoms, g_step);
				for(i = 0; i < g_iGesAtomCount; i++) {
					if(g_waAtomMolIndex[i] == 60000)
						continue;
					fprintf(g_reftrajFile, "%4s %14.10f %14.10f %14.10f\n", ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName, ts->m_vaCoords_Original[i][0] / 100.0f, ts->m_vaCoords_Original[i][1] / 100.0f, ts->m_vaCoords_Original[i][2] / 100.0f);
				}
			}
		} else {
			if(g_step % g_stride == 0) {
				for(i = 0; i < g_ramObserv.GetSize(); i++)
					((CRamanObservation *)g_ramObserv[i])->getDipoleZero();
				for(i = 0; i < (g_orientAvg ? 3 : 1); i++) {
					if(!g_timestep[i]->ReadTimestep(g_polFile[i], false)) {
						mprintf(RED, "Error while reading trajectory for polarizabilities.\n");
						abort();
					}
					if(!g_bSaveCoordsUnchanged) {
						g_timestep[i]->UniteMolecules(false);
						if(g_bRemoveCOM)
							g_timestep[i]->CenterCOM();
					}
					g_timestep[i]->CalcCenters();
					if(g_bWannier)
						g_timestep[i]->ScanWannier(false);
					g_timestep[i]->CalcDipoles(false);
					for(j = 0; j < g_ramObserv.GetSize(); j++) {
						((CRamanObservation *)g_ramObserv[j])->calcPolarizability(i);
					}
				}
			}
		}
	} else {
		ts->CalcPolarizabilities();
		
		for (i = 0; i < g_ramObserv.GetSize(); i++) {
			((CRamanObservation *)g_ramObserv[i])->process();
		}
	}
}

void finalizeRaman() {
	if (g_ramanCompat) {
		int i, j, k, l;
		if(g_newRaman) {
			fclose(g_reftrajFile);
			
			char filename[BUF_SIZE];
			FILE *inputFile;
			for(i = 1; i < (g_orientAvg ? 4 : 2); i++) {
#ifdef TARGET_LINUX
				snprintf(filename, BUF_SIZE, "%s/%d/polarizability.inp", g_ramanDir, i);
#else
				sprintf(filename, "%s/%d/polarizability.inp", g_ramanDir, i);
#endif
				inputFile = fopen(filename, "w");
				if(inputFile == NULL) {
					mprintf(RED, "Could not open input file \"%s\": %s\n", filename, strerror(errno));
					abort();
				}
				
				fprintf(inputFile, "%s", g_inputTemplate);
				for(j = 0; j < 4; j++) {
					char *tempPos = g_inputTemplate;
					for(k = 0; k <= j; k++)
						tempPos = strchr(&tempPos[1], '\0');
					if(tempPos == g_templateFieldPos) {
						fprintf(inputFile, "%.6E", g_fieldStrength);
						fprintf(inputFile, "%s", &g_templateFieldPos[38]);
					} else if(tempPos == g_templatePolPos) {
						if(i == 1)
							fprintf(inputFile, "1.0 0.0 0.0");
						if(i == 2)
							fprintf(inputFile, "0.0 1.0 0.0");
						if(i == 3)
							fprintf(inputFile, "0.0 0.0 1.0");
						fprintf(inputFile, "%s", &g_templatePolPos[43]);
					} else if(tempPos == g_templateCoordPos) {
						CTimeStep *ts = GetTimeStep(0);
						for(l = 0; l < g_iGesAtomCount; l++) {
							if(g_waAtomMolIndex[l] == 60000)
								continue;
							fprintf(inputFile, "   %s %14.10f %14.10f %14.10f\n", ((CAtom *)g_oaAtoms[g_waAtomRealElement[l]])->m_sName, ts->m_vaCoords_Original[l][0] / 100.0f, ts->m_vaCoords_Original[l][1] / 100.0f, ts->m_vaCoords_Original[l][2] / 100.0f);
						}
						fprintf(inputFile, "%s", &g_templateCoordPos[36]);
					} else if(tempPos == g_templateStepsPos) {
						fprintf(inputFile, "%d", g_steps);
						fprintf(inputFile, "%s", &g_templateStepsPos[39]);
					} else {
						mprintf(RED, "Unexpected error while processing the input template.\n");
					}
				}
				fclose(inputFile);
			}
			delete[] g_inputTemplate;
		} else {
			for(i = 0; i < g_ramObserv.GetSize(); i++) {
				mprintf(YELLOW, ">>> Raman Observation %d >>>\n\n", i + 1);
				((CRamanObservation *)g_ramObserv[i])->finalize();
				mprintf(YELLOW, "\n<<< End of Raman Observation %d <<<\n\n", i + 1);
			}
			
			for(i = 0; i < (g_orientAvg ? 3 : 1); i++)
				fclose(g_polFile[i]);
			for(i = 0; i < (g_orientAvg ? 3 : 1); i++)
				delete g_timestep[i];
		}
		delete[] g_ramanDir;
	} else {
		int i;
		for (i = 0; i < g_ramObserv.GetSize(); i++) {
			mprintf(YELLOW, ">>> Raman Observation %d >>>\n\n", i + 1);
			((CRamanObservation *)g_ramObserv[i])->finalize();
			mprintf(YELLOW, "\n<<< End of Raman Observation %d <<<\n\n", i + 1);
		}
	}
}

void parsePolarizability() {
	if (g_bPolarizabilityDefined)
		return;
	
	g_bDipole = true;
	ParseDipole();
	
	mprintf(WHITE, "\n>>> Polarizability Definition >>>\n\n");
	mprintf("    There are the following possibilities to provide polarizabilities:\n");
	mprintf("    (1) Read polarizabilities from an external file\n");
	mprintf("    Calculate polarizabilities by finite differences of dipole moments and\n");
	mprintf("    (2) get dipole moments from Wannier centers\n");
	mprintf("    (3) use Voronoi dipole moments\n");
	mprintf("    (4) load dipole restart files\n");
	mprintf("\n");
	
	while (true) {
		g_iPolarizabilityMode = AskRangeInteger("    Which polarizability mode to set? [2] ", 1, 4, 2);
		
		if (g_iPolarizabilityMode == 1) {
			eprintf("This is not implemented yet.\n");
			continue;
		} else if (g_iPolarizabilityMode == 2) {
			setupWannier();
			break;
		} else if (g_iPolarizabilityMode == 3) {
			if (g_bTegri) {
				if (g_pTetraPak == NULL) {
					try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
					if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					try { g_pTetraPak = new CTetraPak(); } catch(...) { g_pTetraPak = NULL; }
					if (g_pTetraPak == NULL) NewException((double)sizeof(CTetraPak),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					g_pTetraPak->Parse();
				}
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
				continue;
			}
			break;
		} else if (g_iPolarizabilityMode == 4) {
			break;
		} else {
			eprintf("This is impossible.\n");
			continue;
		}
	}
	mprintf("\n");
	
	if (g_iPolarizabilityMode == 1) {
		eprintf("This is not implemented yet.\n");
		abort();
	} else {
		CxString dirname;
		while (true) {
			AskString("    Enter the name of the polarizability directory: [polarizability] ", &dirname, "polarizability");
			if (!FileExist((const char *)dirname)) {
				eprintf("Could not find \"%s\".\n", (const char *)dirname);
				continue;
			}
			break;
		}
		mprintf("\n");
		
		CxString nameString;
		if (g_iPolarizabilityMode == 2)
			nameString = CxString("polarizability_%c%c-wannier.xyz");
		else if (g_iPolarizabilityMode == 3)
			nameString = CxString("polarizability_%c%c-density.cube");
		else if (g_iPolarizabilityMode == 4)
			nameString = CxString("polarizability_%c%c-dipole.restart");
		CxString filename, path;
		filename.sprintf((const char *)nameString, 'x', 'p');
		path.sprintf("%s/%s", (const char *)dirname, (const char *)filename);
		g_iPolarizabilityConf[0] = 0;
		if (FileExist((const char *)path)) {
			g_fPolarizabilityFile[0] = fopen((const char *)path, "r");
			if (g_fPolarizabilityFile[0] != NULL) {
				g_iPolarizabilityConf[0] = 1;
				filename.sprintf((const char *)nameString, 'x', 'n');
				path.sprintf("%s/%s", (const char *)dirname, (const char *)filename);
				if (FileExist((const char *)path)) {
					g_fPolarizabilityFile[1] = fopen((const char *)path, "r");
					if (g_fPolarizabilityFile[1] != NULL) {
						g_iPolarizabilityConf[0] = 2;
					}
				}
			}
		}
		filename.sprintf((const char *)nameString, 'y', 'p');
		path.sprintf("%s/%s", (const char *)dirname, (const char *)filename);
		g_iPolarizabilityConf[1] = 0;
		if (FileExist((const char *)path)) {
			g_fPolarizabilityFile[2] = fopen((const char *)path, "r");
			if (g_fPolarizabilityFile[2] != NULL) {
				g_iPolarizabilityConf[1] = 1;
				filename.sprintf((const char *)nameString, 'y', 'n');
				path.sprintf("%s/%s", (const char *)dirname, (const char *)filename);
				if (FileExist((const char *)path)) {
					g_fPolarizabilityFile[3] = fopen((const char *)path, "r");
					if (g_fPolarizabilityFile[3] != NULL) {
						g_iPolarizabilityConf[1] = 2;
					}
				}
			}
		}
		filename.sprintf((const char *)nameString, 'z', 'p');
		path.sprintf("%s/%s", (const char *)dirname, (const char *)filename);
		g_iPolarizabilityConf[2] = 0;
		if (FileExist((const char *)path)) {
			g_fPolarizabilityFile[4] = fopen((const char *)path, "r");
			if (g_fPolarizabilityFile[4] != NULL) {
				g_iPolarizabilityConf[2] = 1;
				filename.sprintf((const char *)nameString, 'z', 'n');
				path.sprintf("%s/%s", (const char *)dirname, (const char *)filename);
				if (FileExist((const char *)path)) {
					g_fPolarizabilityFile[5] = fopen((const char *)path, "r");
					if (g_fPolarizabilityFile[5] != NULL) {
						g_iPolarizabilityConf[2] = 2;
					}
				}
			}
		}
		
		while (true) {
			mprintf("    The following configuration is set up:\n\n");
			char fieldChar[4] = "xyz";
			int i;
			for (i = 0; i < 3; i++) {
				mprintf("      Electric field along %c axis: ", fieldChar[i]);
				if (g_iPolarizabilityConf[i] == 0)
					mprintf("disabled\n");
				else if (g_iPolarizabilityConf[i] == 1)
					mprintf("forward difference\n");
				else if (g_iPolarizabilityConf[i] == 2)
					mprintf("central difference\n");
				else
					mprintf("UNKNOWN\n");
			}
			mprintf("\n");
			
			if (!AskYesNo("    Change this configuration (y/n)? [no] ", false))
				break;
			
			for (i = 0; i < 3; i++) {
				g_iPolarizabilityConf[i] = AskRangeInteger("    Electric field along %c axis: (1) disable, (2) use forward difference, (3) use central difference? [1] ", 1, 3, 1, fieldChar[i]) - 1;
				if (g_iPolarizabilityConf[i] > 0) {
					while (true) {
						CxString filename;
						FILE *file;
						if (g_iPolarizabilityMode == 2) {
							AskString_ND("    Enter file of Wannier centers with field along positive %c axis: ", &filename, fieldChar[i]);
						} else if (g_iPolarizabilityMode == 3) {
							AskString_ND("    Enter file of electron density with field along positive %c axis: ", &filename, fieldChar[i]);
						} else if (g_iPolarizabilityMode == 4) {
							AskString_ND("    Enter dipole restart file with field along positive %c axis: ", &filename, fieldChar[i]);
						}
						file = fopen((const char *)filename, "r");
						if (file != NULL) {
							if (g_fPolarizabilityFile[2 * i] != NULL)
								fclose(g_fPolarizabilityFile[2 * i]);
							g_fPolarizabilityFile[2 * i] = file;
							break;
						}
						eprintf("Could not open \"%s\": %s.\n", (const char *)filename, strerror(errno));
					}
					if (g_iPolarizabilityConf[i] > 1) {
						while (true) {
							CxString filename;
							FILE *file;
							if (g_iPolarizabilityMode == 2) {
								AskString_ND("    Enter file of Wannier centers with field along negative %c axis: ", &filename, fieldChar[i]);
							} else if (g_iPolarizabilityMode == 3) {
								AskString_ND("    Enter file of electron density with field along negative %c axis: ", &filename, fieldChar[i]);
							} else if (g_iPolarizabilityMode == 4) {
								AskString_ND("    Enter dipole restart file with field along negative %c axis: ", &filename, fieldChar[i]);
							}
							file = fopen((const char *)filename, "r");
							if (file != NULL) {
								if (g_fPolarizabilityFile[2 * i + 1] != NULL)
									fclose(g_fPolarizabilityFile[2 * i + 1]);
								g_fPolarizabilityFile[2 * i + 1] = file;
								break;
							}
							eprintf("Could not open \"%s\": %s.\n", (const char *)filename, strerror(errno));
						}
					}
				}
			}
		}
		mprintf("\n");
		
		if (g_iPolarizabilityMode == 4) {
			int i;
			for (i = 0; i < 3; i++) {
				if (g_iPolarizabilityConf[i] > 0) {
					int numAtoms;
					fread(&numAtoms, sizeof(int), 1, g_fPolarizabilityFile[2 * i]);
					char fieldChar[4] = "xyz";
					if (numAtoms != g_oaSingleMolecules.GetSize()) {
						eprintf("The dipole restart file with field along positive %c axis was written for a different number of molecules.\n", fieldChar[i]);
						abort();
					}
					if (g_iPolarizabilityConf[i] > 1) {
						int numAtoms;
						fread(&numAtoms, sizeof(int), 1, g_fPolarizabilityFile[2 * i + 1]);
						char fieldChar[4] = "xyz";
						if (numAtoms != g_oaSingleMolecules.GetSize()) {
							eprintf("The dipole restart file with field along negative %c axis was written for a different number of molecules.\n", fieldChar[i]);
							abort();
						}
					}
				}
			}
		}
		
		filename.sprintf("%s/settings.dat", (const char *)dirname);
		FILE *settingsFile = fopen((const char *)filename, "r");
		if (settingsFile != NULL) {
			CxString line;
			line.fgets(1024, settingsFile);
			if (sscanf((const char *)line, "%f", &g_fPolarizabilityFieldStrength) < 1)
				g_fPolarizabilityFieldStrength = 1.0f;
			else
				mprintf("    Found electric field strentgh in \"%s/settings.dat\": %g a. u.\n", (const char *)dirname, g_fPolarizabilityFieldStrength);
			fclose(settingsFile);
		}
		
		g_fPolarizabilityFieldStrength = AskFloat("    Electric field strength in atomic units: [%g] ", g_fPolarizabilityFieldStrength, g_fPolarizabilityFieldStrength);
	}
	
	g_bPolarizabilityDefined = true;
	mprintf(WHITE, "\n<<< End of Polarizability Definition <<<\n\n");
}

static CxString g_polarizabilityDir;
static float g_polFieldStrength;
static char *g_polInputTemplate = NULL;
static const char *g_polInputTemplateEnd = NULL;
static char *g_polTemplatePosIdent[2] = { NULL, NULL };
static char *g_polTemplatePosField = NULL;
static char *g_polTemplatePosPol = NULL;
static char *g_polTemplatePosCoord = NULL;
static char *g_polTemplatePosSteps = NULL;
static FILE *g_polReftrajFile = NULL;
static int g_polSteps = 0;

bool gatherPolarizabilityCalc() {

#ifdef TARGET_LINUX
	g_bKeepOriginalCoords = true;
	
	AskString("    Enter a name for the directory to collect the data: [polarizability] ", &g_polarizabilityDir, "polarizability");
	if (FileExist((const char *)g_polarizabilityDir)) {
		eprintf("A file or a directory \"%s\" already exists. Please remove it first.\n", (const char *)g_polarizabilityDir);
		return false;
	}
	if (mkdir((const char *)g_polarizabilityDir, S_IRWXU) != 0) {
		eprintf("Directory \"%s\" could not be created: %s\n", (const char *)g_polarizabilityDir, strerror(errno));
		return false;
	}
	
	mprintf("\n");
	mprintf("    A set of six CP2K input files will be created in the directory \"%s\".\n", (const char *)g_polarizabilityDir);
	mprintf("    These contain the three independent field directions (x, y, z) with positive and negative sign (p, n).\n");
	mprintf("    To get the full polarizability tensor with forward differences, you need to calculate xp, yp, and zp.\n");
	mprintf("    To use (more accurate) central differences, you also need to calculate xn, yn, and zn.\n");
	mprintf("    For isotropic systems, it is usually sufficient to consider only one field direction.\n");
	mprintf("    In this case, you only need xp for forward differences, and also xn for central differences.\n");
	mprintf("\n");
	
	CxString filename;
	filename.sprintf("%s/settings.dat", (const char *)g_polarizabilityDir);
	FILE *settingsFile = fopen((const char *)filename, "w");
	if (settingsFile == NULL) {
		eprintf("Could not open settings file \"%s\": %s\n", (const char *)filename, strerror(errno));
		return false;
	}
	g_polFieldStrength = AskFloat("    Electric field strength in atomic units: [5e-4] ", 5.0e-4f);
	fprintf(settingsFile, "%.6E\n", g_polFieldStrength);
	fclose(settingsFile);
	
	int mode = AskRangeInteger_ND("    Prepare input for Wannier centers (1), cube trajectory (2), or cube streaming (3)? ", 1, 3);
	
	filename.sprintf("%s/template.inp", (const char *)g_polarizabilityDir);
	FILE *templateFile = fopen((const char *)filename, "w");
	if(templateFile == NULL) {
		mprintf(RED, "Could not open template file \"%s\": %s\n", (const char *)filename, strerror(errno));
		return false;
	}
	fprintf(templateFile, "&GLOBAL\n");
	fprintf(templateFile, " PROJECT_NAME polarizability_###!identifier will be put here###\n");
	fprintf(templateFile, " RUN_TYPE MD\n");
	fprintf(templateFile, " PRINT_LEVEL LOW\n");
	fprintf(templateFile, " FFTW_PLAN_TYPE PATIENT\n");
	fprintf(templateFile, "&END\n");
	fprintf(templateFile, "&FORCE_EVAL\n");
	fprintf(templateFile, " &DFT\n");
	fprintf(templateFile, "  BASIS_SET_FILE_NAME BASIS_MOLOPT\n");
	fprintf(templateFile, "  POTENTIAL_FILE_NAME POTENTIAL\n");
	fprintf(templateFile, "  &MGRID\n");
	fprintf(templateFile, "   NGRIDS 5\n");
	fprintf(templateFile, "   CUTOFF 280\n");
	fprintf(templateFile, "   REL_CUTOFF 40\n");
	fprintf(templateFile, "  &END\n");
	fprintf(templateFile, "  &SCF\n");
	fprintf(templateFile, "   SCF_GUESS ATOMIC\n");
	fprintf(templateFile, "   MAX_SCF 200\n");
	fprintf(templateFile, "   EPS_SCF 1.0E-5\n");
	fprintf(templateFile, "   &OT\n");
	fprintf(templateFile, "    MINIMIZER DIIS\n");
	fprintf(templateFile, "    PRECONDITIONER FULL_SINGLE_INVERSE\n");
	fprintf(templateFile, "   &END\n");
	fprintf(templateFile, "  &END\n");
	fprintf(templateFile, "  &XC\n");
	fprintf(templateFile, "   &XC_FUNCTIONAL BLYP\n");
	fprintf(templateFile, "   &END\n");
	fprintf(templateFile, "   &XC_GRID\n");
	fprintf(templateFile, "    XC_SMOOTH_RHO NN10\n");
	fprintf(templateFile, "    XC_DERIV NN10_SMOOTH\n");
	fprintf(templateFile, "   &END\n");
	fprintf(templateFile, "   &VDW_POTENTIAL\n");
	fprintf(templateFile, "    POTENTIAL_TYPE PAIR_POTENTIAL\n");
	fprintf(templateFile, "    &PAIR_POTENTIAL\n");
	fprintf(templateFile, "     TYPE DFTD3\n");
	fprintf(templateFile, "     REFERENCE_FUNCTIONAL BLYP\n");
	fprintf(templateFile, "     PARAMETER_FILE_NAME dftd3.dat\n");
	fprintf(templateFile, "    &END\n");
	fprintf(templateFile, "   &END\n");
	fprintf(templateFile, "  &END\n");
	fprintf(templateFile, "  &PERIODIC_EFIELD\n");
	fprintf(templateFile, "   INTENSITY ###!field strength will be put here###\n");
	fprintf(templateFile, "   POLARISATION ###!polarisation vector will be put here###\n");
	fprintf(templateFile, "  &END\n");
	if (mode == 1) {
		fprintf(templateFile, "  &LOCALIZE\n");
		fprintf(templateFile, "   METHOD CRAZY\n");
		fprintf(templateFile, "   MAX_ITER 2000\n");
		fprintf(templateFile, "   &PRINT\n");
		fprintf(templateFile, "    &WANNIER_CENTERS\n");
		fprintf(templateFile, "     IONS+CENTERS\n");
		fprintf(templateFile, "     FILENAME =polarizability_###!identifier will be put here###-wannier.xyz\n");
		fprintf(templateFile, "    &END\n");
		fprintf(templateFile, "   &END\n");
		fprintf(templateFile, "  &END\n");
	} else if (mode == 2) {
		fprintf(templateFile, "  &PRINT\n");
		fprintf(templateFile, "   &E_DENSITY_CUBE\n");
		fprintf(templateFile, "    FILENAME =polarizability_###!identifier will be put here###-density.cube\n");
		fprintf(templateFile, "    APPEND\n");
		fprintf(templateFile, "    STRIDE 2 2 2\n");
		fprintf(templateFile, "   &END\n");
		fprintf(templateFile, "  &END\n");
	} else if (mode == 3) {
		fprintf(templateFile, "  &PRINT\n");
		fprintf(templateFile, "   &E_DENSITY_CUBE\n");
		fprintf(templateFile, "    STRIDE 2 2 2\n");
		fprintf(templateFile, "   &END\n");
		fprintf(templateFile, "  &END\n");
	}
	fprintf(templateFile, " &END\n");
	fprintf(templateFile, " &SUBSYS\n");
	fprintf(templateFile, "  &CELL\n");
	fprintf(templateFile, "   ABC %.6f %.6f %.6f\n", g_fBoxX / 100.0f, g_fBoxY / 100.0f, g_fBoxZ / 100.0f);
	fprintf(templateFile, "  &END\n");
	fprintf(templateFile, "  &COORD\n");
	fprintf(templateFile, "###!coordinates will be put here###\n");
	fprintf(templateFile, "  &END\n");
	int i;
	for (i = 0; i < g_oaAtoms.GetSize(); i++) {
		if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "H") == 0) {
			fprintf(templateFile, "  &KIND H\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q1\n");
			fprintf(templateFile, "  &END\n");
		} else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "C") == 0) {
			fprintf(templateFile, "  &KIND C\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q4\n");
			fprintf(templateFile, "  &END\n");
		} else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "N") == 0) {
			fprintf(templateFile, "  &KIND N\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q5\n");
			fprintf(templateFile, "  &END\n");
		} else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "O") == 0) {
			fprintf(templateFile, "  &KIND O\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q6\n");
			fprintf(templateFile, "  &END\n");
		} else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "F") == 0) {
			fprintf(templateFile, "  &KIND F\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q7\n");
			fprintf(templateFile, "  &END\n");
		} else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "S") == 0) {
			fprintf(templateFile, "  &KIND S\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q6\n");
			fprintf(templateFile, "  &END\n");
		} else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Cl") == 0) {
			fprintf(templateFile, "  &KIND Cl\n");
			fprintf(templateFile, "   BASIS_SET DZVP-MOLOPT-SR-GTH\n");
			fprintf(templateFile, "   POTENTIAL GTH-BLYP-q7\n");
			fprintf(templateFile, "  &END\n");
		}
	}
	fprintf(templateFile, " &END\n");
	fprintf(templateFile, "&END\n");
	fprintf(templateFile, "&MOTION\n");
	fprintf(templateFile, " &MD\n");
	fprintf(templateFile, "  ENSEMBLE REFTRAJ\n");
	fprintf(templateFile, "  STEPS ###!number of steps will be put here###\n");
	fprintf(templateFile, "  TIMESTEP %f\n", g_fTimestepLength);
	fprintf(templateFile, "  &REFTRAJ\n");
	fprintf(templateFile, "   EVAL_ENERGY_FORCES\n");
	fprintf(templateFile, "   TRAJ_FILE_NAME polarizability_reftraj.xyz\n");
	fprintf(templateFile, "  &END\n");
	fprintf(templateFile, " &END\n");
	fprintf(templateFile, "&END\n");
	fclose(templateFile);
	
	mprintf("\n    An input template for the polarizability calculations has been written to \"%s/template.inp\"\n    Modify it according to your needs and press any key when you are finished.\n", (const char *)g_polarizabilityDir);
	getchar();
	
#else
	eprintf("\nRaman calculations are currently only supported with TARGET_LINUX.\n");
	abort();
#endif

	return true;
}

bool initializePolarizabilityCalc() {
	CxString filename;
	filename.sprintf("%s/template.inp", (const char *)g_polarizabilityDir);
	FILE *templateFile = fopen(filename, "r");
	if (templateFile == NULL) {
		eprintf("Could not open template file \"%s\": %s\n", (const char *)filename, strerror(errno));
		return false;
	}
	
	fseek(templateFile, 0L, SEEK_END);
	long length = ftell(templateFile);
	if (length < 0) {
		eprintf("Could not determine size of template file: %s\n", strerror(errno));
		fclose(templateFile);
		return false;
	}
	try { g_polInputTemplate = new char[length + 1]; } catch(...) { g_polInputTemplate = NULL; }
	if (g_polInputTemplate == NULL) NewException((double)(length + 1) * sizeof(char), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	rewind(templateFile);
	if ((long)fread(g_polInputTemplate, sizeof(char), length, templateFile) < length) {
		eprintf("Could not read template file: %s\n", strerror(errno));
		fclose(templateFile);
		return false;
	}
	g_polInputTemplate[length] = 0;
	g_polInputTemplateEnd = &g_polInputTemplate[length];
	fclose(templateFile);
	
	g_polTemplatePosIdent[0] = strstr(g_polInputTemplate, "###!identifier will be put here###");
	g_polTemplatePosIdent[1] = strstr(g_polTemplatePosIdent[0] + 1, "###!identifier will be put here###");
	g_polTemplatePosField = strstr(g_polInputTemplate, "###!field strength will be put here###");
	if(g_polTemplatePosField == NULL) {
		eprintf("Position mark for field strength missing in template.\n");
		fclose(templateFile);
		return false;
	}
	g_polTemplatePosPol = strstr(g_polInputTemplate, "###!polarisation vector will be put here###");
	if(g_polTemplatePosPol == NULL) {
		eprintf("Position mark for polarisation vector missing in template.\n");
		fclose(templateFile);
		return false;
	}
	g_polTemplatePosCoord = strstr(g_polInputTemplate, "###!coordinates will be put here###");
	if(g_polTemplatePosCoord == NULL) {
		eprintf("Position mark for coordinates missing in template.\n");
		fclose(templateFile);
		return false;
	}
	g_polTemplatePosSteps = strstr(g_polInputTemplate, "###!number of steps will be put here###");
	if(g_polTemplatePosSteps == NULL) {
		eprintf("Position mark for number of steps missing in template.\n");
		fclose(templateFile);
		return false;
	}
	if (g_polTemplatePosIdent[0] != NULL)
		g_polTemplatePosIdent[0][0] = 0;
	if (g_polTemplatePosIdent[1] != NULL)
		g_polTemplatePosIdent[1][0] = 0;
	g_polTemplatePosField[0] = 0;
	g_polTemplatePosPol[0] = 0;
	g_polTemplatePosCoord[0] = 0;
	g_polTemplatePosSteps[0] = 0;
	
	filename.sprintf("%s/polarizability_reftraj.xyz", (const char *)g_polarizabilityDir);
	g_polReftrajFile = fopen((const char *)filename, "w");
	if(g_polReftrajFile == NULL) {
		eprintf("Could not open reference trajectory \"%s\": %s\n", (const char *)filename, strerror(errno));
		fclose(templateFile);
		return false;
	}
		
	return true;
}

void processPolarizabilityCalc(CTimeStep *ts) {
	int numAtoms = 0;
	int i;
	for (i = 0; i < g_iGesAtomCount; i++)
		if (g_waAtomMolIndex[i] != 60000)
			numAtoms++;
	fprintf(g_polReftrajFile, "%d\n", numAtoms);
	fprintf(g_polReftrajFile, "Step %lu\n", g_iSteps);
	for (i = 0; i < g_iGesAtomCount; i++) {
		if (g_waAtomMolIndex[i] == 60000)
			continue;
		fprintf(g_polReftrajFile, "%3s %16.10f %16.10f %16.10f\n", ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName, ts->m_vaCoords_Original[i][0] / 100.0, ts->m_vaCoords_Original[i][1] / 100.0, ts->m_vaCoords_Original[i][2] / 100.0);
	}
	g_polSteps++;
}

void finalizePolarizabilityCalc() {
	fclose(g_polReftrajFile);
	
	char fieldChar[4] = "xyz";
	char dirChar[3] = "pn";
	CxString filename;
	FILE *inputFile;
	int i;
	for (i = 0; i < 3; i++) {
		int j;
		for (j = 0; j < 2; j++) {
			filename.sprintf("%s/polarizability_%c%c.inp", (const char *)g_polarizabilityDir, fieldChar[i], dirChar[j]);
			inputFile = fopen((const char *)filename, "w");
			if (inputFile == NULL) {
				eprintf("Could not open input file \"%s\": %s\n", (const char *)filename, strerror(errno));
				abort();
			}
			
			fprintf(inputFile, "%s", g_polInputTemplate);
			const char *tempPos = strchr(g_polInputTemplate, 0);
			while (tempPos != g_polInputTemplateEnd) {
				if (tempPos == g_polTemplatePosIdent[0] || tempPos == g_polTemplatePosIdent[1]) {
					fprintf(inputFile, "%c%c", fieldChar[i], dirChar[j]);
					tempPos = &tempPos[34];
				} else if (tempPos == g_polTemplatePosField) {
					fprintf(inputFile, "%.6E", g_polFieldStrength);
					tempPos = &tempPos[38];
				} else if (tempPos == g_polTemplatePosPol) {
					if (i == 0)
						fprintf(inputFile, "%c1.0 0.0 0.0", j == 0 ? '+' : '-');
					if (i == 1)
						fprintf(inputFile, "0.0 %c1.0 0.0", j == 0 ? '+' : '-');
					if (i == 2)
						fprintf(inputFile, "0.0 0.0 %c1.0", j == 0 ? '+' : '-');
					tempPos = &tempPos[43];
				} else if (tempPos == g_polTemplatePosCoord) {
					CTimeStep *ts = GetTimeStep(0);
					int k;
					for (k = 0; k < g_iGesAtomCount; k++) {
						if(g_waAtomMolIndex[k] == 60000)
							continue;
						fprintf(inputFile, "   %s %.10f %.10f %.10f\n", ((CAtom *)g_oaAtoms[g_waAtomRealElement[k]])->m_sName, ts->m_vaCoords_Original[k][0] / 100.0, ts->m_vaCoords_Original[k][1] / 100.0, ts->m_vaCoords_Original[k][2] / 100.0);
					}
					tempPos = &tempPos[36];
				} else if (tempPos == g_polTemplatePosSteps) {
					fprintf(inputFile, "%d", g_polSteps);
					tempPos = &tempPos[39];
				}
				fprintf(inputFile, "%s", tempPos);
				tempPos = strchr(tempPos, 0);
			}
			fclose(inputFile);
		}
	}
	
	delete[] g_polInputTemplate;
}

