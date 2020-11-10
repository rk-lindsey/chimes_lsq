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

#include "ir.h"

#include "globalvar.h"
#include "maintools.h"
#include "timestep.h"
#include "xstring.h"

#include <math.h>

#define BUF_SIZE 4096

static CxObArray g_PowerObserv;
// static CxFloatArray g_TemperatureCache;

CPowerObservation::CPowerObservation(bool global) {
	int i;
	CxString buf, buf2;
	
	try { _atoms = new CAtomGroup(); } catch(...) { _atoms = NULL; }
	if(_atoms == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	if(global) {
		m_iShowMol = -1;
		m_iShowMolCount = g_oaSingleMolecules.GetSize();
		_name = new char[7];
		sprintf(_name, "global");
	} else {
//		char buf[BUF_SIZE];
//		char buf2[BUF_SIZE];
//		size_t remaining = BUF_SIZE;
		if(g_oaMolecules.GetSize() > 1) {
/*#ifdef TARGET_LINUX
			remaining -= snprintf(buf, remaining, "    Which molecule should be observed (");
#else
			remaining -= sprintf(buf, "    Which molecule should be observed (");
#endif*/
			buf.sprintf( "    Which molecule should be observed (");
			for(i = 0; i < g_oaMolecules.GetSize(); i++) {

/*				if(remaining < 1)
					break;
#ifdef TARGET_LINUX
				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#else
				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;*/

				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
				buf.strcat(buf2);

				if(i < g_oaMolecules.GetSize() - 1) {
/*#ifdef TARGET_LINUX
					length = snprintf(buf2, remaining, ", ");
#else
					length = sprintf(buf2, ", ");
#endif
					strncat(buf, buf2, remaining - 1);
					remaining -= length;*/

					buf2.sprintf(", ");
					buf.strcat(buf2);
				}
			}
//			strncat(buf, ")? ", remaining - 1);
			buf.strcat(")? ");
			m_iShowMol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
		} else {
			m_iShowMol = 0;
			mprintf("    Observing molecule %s.\n", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		}
		
		while(true) {
			mprintf("    Which atom(s) to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [*] ");
			inpprintf("! Which atom(s) to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [*]\n");
//			char buf[BUF_SIZE];
			myget(&buf);
			if(strlen(buf) == 0) {
				if(!_atoms->ParseAtoms((CMolecule *)g_oaMolecules[m_iShowMol], "*")) {
					eprintf("Weird error.\n");
					continue;
				}
			} else if(!_atoms->ParseAtoms((CMolecule *)g_oaMolecules[m_iShowMol], buf)) {
				continue;
			}
			break;
		}
		mprintf("\n    Observing %d atoms of molecule %s.\n", _atoms->m_iAtomGes, ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
		_name = new char[strlen(((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName) + strlen(_atoms->m_sName) + 4];
		sprintf(_name, "[%s_%s]", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, _atoms->m_sName);
		mprintf("\n");
	}
	
	if(g_iTrajSteps != -1) {
		_correlationDepth = (int)(0.75 * g_iTrajSteps);
		if(_correlationDepth > 4096)
			_correlationDepth = 4096;
		if(g_fTimestepLength > 1.0f)
			_correlationDepth = 2048;
		if(g_fTimestepLength > 2.0f)
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
		mprintf("    Using cos^2 window function; inserting %d zeros for zero padding.\n",_zeroPadding);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
	}
	
	float possibleRange = 33356.41f / g_fTimestepLength / 2.0f;
	_specResolution = possibleRange / (_correlationDepth + _zeroPadding);
	mprintf("    This results in a spectral resolution of %.2f cm^-1.\n", _specResolution);
	mprintf("\n    A time step length of %.2f fs allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, possibleRange);
	float specLimit = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ", 0, possibleRange, (possibleRange < 5000.0f) ? possibleRange : 5000.0f, (possibleRange < 5000.0f) ? possibleRange : 5000.0f);
	_specSize = (int)(specLimit / _specResolution);
	mprintf("\n");
	
	_massWeighting = AskYesNo("    Weight autocorrelation functions by atomic mass (y/n)? [yes] ", true);
	mprintf("\n");
	
	if(g_bAdvanced2) {
		_saveACF = AskYesNo("    Save autocorrelation function (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_saveACF = false;
	}

}


CPowerObservation::~CPowerObservation() {
	delete[] _name;
	delete _atoms;
}


void CPowerObservation::initialize() {
	int i, j, k;
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	if(m_iShowMol == -1) {
		mprintf("    Velocity cache: Trying to allocate %s of memory...\n", FormatBytes((double)g_iGesAtomCount * n * sizeof(CxVector3)));
		for(i = 0; i < g_iGesAtomCount; i++) {
			CxVec3Array *a;
			try { a = new CxVec3Array(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			a->SetMaxSize(n);
			a->SetGrow(n / 10);
			_velocityCache.Add(a);
		}
	} else {
		mprintf("    Velocity cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * _atoms->m_iAtomGes * n * sizeof(CxVector3)));
		for(i = 0; i < m_iShowMolCount * _atoms->m_iAtomGes; i++) {
			CxVec3Array *a;
			try { a = new CxVec3Array(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			a->SetMaxSize(n);
			a->SetGrow(n / 10);
			_velocityCache.Add(a);
		}
	}
	
	if(m_iShowMol == -1) {
		_masses.SetSize(g_iGesAtomCount);
		if(_massWeighting) {
			for(i = 0; i < g_iGesAtomCount; i++) {
				_masses[i] = ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_pElement->m_fMass;
			}
		} else {
			for(i = 0; i < g_iGesAtomCount; i++) {
				_masses[i] = 1.0f;
			}
		}
	} else {
		_masses.SetSize(m_iShowMolCount * _atoms->m_iAtomGes);
		if(_massWeighting) {
			for(i = 0; i < m_iShowMolCount; i++) {
				int n = 0;
				for(j = 0; j < _atoms->m_baRealAtomType.GetSize(); j++) {
					CxIntArray *a = (CxIntArray *)_atoms->m_oaAtoms[j];
					for(k = 0; k < a->GetSize(); k++) {
						_masses[i * _atoms->m_iAtomGes + n] = ((CAtom *)g_oaAtoms[_atoms->m_baRealAtomType[j]])->m_pElement->m_fMass;
						n++;
					}
				}
			}
		} else {
			for(i = 0; i < m_iShowMolCount * _atoms->m_iAtomGes; i++) {
				_masses[i] = 1.0f;
			}
		}
	}
}


void CPowerObservation::process(CTimeStep *ts) {
	int i, j, k;
	if(m_iShowMol == -1) {
		for(i = 0; i < g_iGesAtomCount; i++) {
			((CxVec3Array *)_velocityCache[i])->Add(ts->m_vaVelocities[i]);
		}
	} else {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
			int n = 0;
			for(j = 0; j < _atoms->m_baAtomType.GetSize(); j++) {
				CxIntArray *a = (CxIntArray *)_atoms->m_oaAtoms[j];
				for(k = 0; k < a->GetSize(); k++) {
					((CxVec3Array *)_velocityCache[i * _atoms->m_iAtomGes + n])->Add(ts->m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[_atoms->m_baAtomType[j]])->GetAt(a->GetAt(k))]);
					n++;
				}
			}
		}
	}
}


void CPowerObservation::finalize() {
	int i, j, k, l, n;
	float step;

	n = ((CxVec3Array *)_velocityCache[0])->GetSize();
	if (n < _correlationDepth)
	{
		eprintf("\nError: Autocorrelation depth is %d, but only %d timesteps evaluated.\n",_correlationDepth,n);
		eprintf("       Reduce depth or increase trajectory length.\n\n");
		abort();
	}

	if(m_iShowMol == -1) {
		step = (float)g_iGesAtomCount / 20.0f;
	} else {
		step = (float)m_iShowMolCount * _atoms->m_iAtomGes / 20.0f;
	}
	
	mprintf("    Calculating autocorrelation...\n");
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
	temp->SetSize(n);
	CxFloatArray *temp2;
	try { temp2 = new CxFloatArray(); } catch(...) { temp2 = NULL; }
	if(temp2 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp2->SetSize(_correlationDepth);
	
	mprintf(WHITE, "     [");
	if(m_iShowMol == -1) {
		for(i = 0; i < g_iGesAtomCount; i++) {
			if(fmodf(i, step) < 1.0f)
				mprintf(WHITE, "#");
			for(j = 0; j < 3; j++) {
				for(k = 0; k < n; k++) {
					temp->GetAt(k) = ((CxVec3Array *)_velocityCache[i])->GetAt(k)[j];
				}
				ac->AutoCorrelate(temp, temp2);
				for(k = 0; k < _correlationDepth; k++) {
					acf->GetAt(k) += temp2->GetAt(k) * _masses[i];
				}
			}
		}
	} else {
		for(i = 0; i < m_iShowMolCount; i++) {
			for(j = 0; j < _atoms->m_iAtomGes; j++) {
				if(fmodf(i * _atoms->m_iAtomGes + j, step) < 1.0f)
					mprintf(WHITE, "#");
				for(k = 0; k < 3; k++) {
					for(l = 0; l < n; l++) {
						temp->GetAt(l) = ((CxVec3Array *)_velocityCache[i * _atoms->m_iAtomGes + j])->GetAt(l)[k];
					}
					ac->AutoCorrelate(temp, temp2);
					for(l = 0; l < _correlationDepth; l++) {
						acf->GetAt(l) += temp2->GetAt(l) * _masses[i * _atoms->m_iAtomGes + j];
					}
				}
			}
		}
	}
	mprintf(WHITE, "]\n");
	
	if(m_iShowMol != -1) {
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) /= (float)m_iShowMolCount;
		}
	}
	
	delete ac;
	delete temp2;
	
	if(_saveACF) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "power_acf_%s.csv", _name);
#else
		sprintf(filename, "power_acf_%s.csv", _name);
#endif
		mprintf("    Saving autocorrelation function as %s...\n", filename);
		FILE *acfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(acfFile, "%.2f; %.8G\n", i * g_fTimestepLength, acf->GetAt(i));
		}
		fclose(acfFile);
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
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "power_acf_windowed_%s.csv", _name);
#else
		sprintf(filename, "power_acf_windowed_%s.csv", _name);
#endif
		mprintf("    Saving windowed autocorrelation function as %s...\n", filename);
		FILE *acfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(acfFile, "%.2f; %.8G\n", i * g_fTimestepLength, temp->GetAt(i));
		}
		fclose(acfFile);
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
	
	mprintf("    Fourier transforming autocorrelation function...\n");
	CFFT *fft;
	try { fft = new CFFT(); } catch(...) { fft = NULL; }
	if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	fft->PrepareFFT_C2C(temp->GetSize());
	for(i = 0; i < temp->GetSize(); i++) {
		fft->m_pInput[2*i] = temp->GetAt(i);
		fft->m_pInput[2*i+1] = 0.0f;
	}
	fft->DoFFT();
	
	CxFloatArray *spectrum;
	try { spectrum = new CxFloatArray(); } catch(...) { spectrum = NULL; }
	if(spectrum == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	spectrum->SetSize(_specSize);
	for(i = 0; i < _specSize; i++) {
		spectrum->GetAt(i) = 7.211349e-9f * fft->m_pOutput[2*i] * g_fTimestepLength; // Output in K*cm
	}
	
	char filename[BUF_SIZE];
#ifdef TARGET_LINUX
	snprintf(filename, BUF_SIZE, "power_spectrum_%s.csv", _name);
#else
	sprintf(filename, "power_spectrum_%s.csv", _name);
#endif
	mprintf("    Saving spectrum as %s...\n", filename);
	FILE *specFile = OpenFileWrite(filename, false);

		fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm); Integral (K)\n");

	double integral = 0.0;
	for(i = 0; i < _specSize; i++) {
		integral += (double)spectrum->GetAt(i) * _specResolution;

			fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum->GetAt(i), integral);
	}
	fclose(specFile);
	
	if(m_iShowMol == -1)
		mprintf("\n    Assuming %d degrees of freedom, the average temperature is %.2f K\n", 3 * g_iGesAtomCount, integral / (3.0 * g_iGesAtomCount));
	else
		mprintf("\n    Assuming %d degrees of freedom, the average temperature is %.2f K\n", 3 * _atoms->m_iAtomGes, integral / (3.0 * _atoms->m_iAtomGes));
	
	delete acf;
	delete temp;
	delete fft;
	delete spectrum;
	
	for(i = 0; i < _velocityCache.GetSize(); i++) {
		delete (CxVec3Array *)_velocityCache[i];
	}
}
	

bool gatherPowerSpectrum() {
	g_bUseVelocities = true;
	
	while(true) {
		mprintf(YELLOW, "\n>>> Power Spectrum %d >>>\n\n", g_PowerObserv.GetSize() + 1);
		
		CPowerObservation *obs;
		try { obs = new CPowerObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CPowerObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_PowerObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of Power Spectrum %d <<<\n\n", g_PowerObserv.GetSize());
		
		if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	
	if(AskYesNo("\n    Compute power spectrum of whole system (y/n) [no] ", false)) {
		mprintf(YELLOW, "\n>>> Global Power Spectrum >>>\n\n");
		
		CPowerObservation *obs;
		try { obs = new CPowerObservation(true); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CPowerObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_PowerObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of Global Power Spectrum <<<\n\n");
	}
	
	return true;
}


bool initializePowerSpectrum() {
	int i;
	for(i = 0; i < g_PowerObserv.GetSize(); i++) {
		((CPowerObservation *)g_PowerObserv[i])->initialize();
	}
	
// 	int n;
// 	if(g_iTrajSteps != -1)
// 		n = (int)(1.1 * g_iTrajSteps / g_iStride);
// 	else
// 		n = 10000;
// 	
// 	mprintf("    Temperature cache: Trying to allocate %s of memory...\n", FormatBytes((double)g_iGesAtomCount * n * sizeof(float)));
// 	g_TemperatureCache.SetMaxSize(n);
// 	g_TemperatureCache.SetGrow(n / 10);
	
	return true;
}


void processPowerSpectrum(CTimeStep* ts) {
	int i;
	for(i = 0; i < g_PowerObserv.GetSize(); i++) {
		((CPowerObservation *)g_PowerObserv[i])->process(ts);
	}
	
// 	float temperature = 0.0f;
// 	for (i = 0; i < g_iGesAtomCount; i++) {
// 		temperature += ts->m_vaVelocities[i].GetLengthSqr() * ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_pElement->m_fMass;
// 	}
// 	temperature /= (3.0f * g_iGesAtomCount * 8314.4621f);
// 	g_TemperatureCache.Add(temperature);
}


void finalizePowerSpectrum() {
	int i;
	for(i = 0; i < g_PowerObserv.GetSize(); i++) {
		mprintf(YELLOW, "\n>>> Power Spectrum %d >>>\n\n", i+1);
		((CPowerObservation *)g_PowerObserv[i])->finalize();
		delete (CPowerObservation *)g_PowerObserv[i];
		mprintf(YELLOW, "\n<<< End of Power Spectrum %d <<<\n\n", i+1);
	}
	
// 	CxString filename;
// 	filename.sprintf("temperature.csv");
// 	mprintf("    Saving temperature as %s...\n", (const char *)filename);
// 	FILE *tempFile = OpenFileWrite(filename, false);
// 	fprintf(tempFile, "#Time (fs); Temperature (K)\n");
// 	double average = 0.0;
// 	for(i = 0; i < g_TemperatureCache.GetSize(); i++) {
// 		average += (double)g_TemperatureCache[i];
// 		fprintf(tempFile, "%.2f; %.8G\n", g_fTimestepLength * i, g_TemperatureCache[i]);
// 	}
// 	fclose(tempFile);
// 	average /= g_TemperatureCache.GetSize();
// 	
// 	mprintf("\n    The average temperature is %.2f K\n", average);
}


static CxObArray g_IRObserv;


CIRObservation::CIRObservation(bool global) {
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
	
	if(g_iTrajSteps != -1) {
		_correlationDepth = (int)(0.75 * g_iTrajSteps);
		if(_correlationDepth > 4096)
			_correlationDepth = 4096;
		if(g_fTimestepLength > 1.0f)
			_correlationDepth = 2048;
		if(g_fTimestepLength > 2.0f)
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
		mprintf("    Using cos^2 window function.\n");
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
		mprintf("    Inserting %d zeros for zero padding.\n",_zeroPadding);
	}
	
	float possibleRange = 33356.41f / g_fTimestepLength / 2.0f;
	_specResolution = possibleRange / (_correlationDepth + _zeroPadding);
	mprintf("    This results in a spectral resolution of %.2f cm^-1.\n", _specResolution);
	mprintf("\n    A time step length of %.2f fs allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, possibleRange);
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
		_saveACF = AskYesNo("    Save autocorrelation function (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_saveACF = false;
	}
	
	if(g_bAdvanced2 && m_iShowMol == -1) {
		_includeCross = AskYesNo("    Include also cross-correlations (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_includeCross = false;
	}
	
	
}

CIRObservation::~CIRObservation() {
	delete[] _name;
}

void CIRObservation::initialize() {
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	mprintf("    Moment cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * n * sizeof(CxVector3)));
	int i;
	for(i = 0; i < m_iShowMolCount; i++) {
		CxVec3Array *a;
		try { a = new CxVec3Array(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetMaxSize(n);
		a->SetGrow(n / 10);
		_dipoleCache.Add(a);
	}
}

void CIRObservation::process() {
	int i;
	if(m_iShowMol == -1) {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
			((CxVec3Array *)_dipoleCache[i])->Add(sm->m_vDipole);
		}
	} else {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
			((CxVec3Array *)_dipoleCache[i])->Add(sm->m_vDipole);
		}
	}
}

void CIRObservation::finalize() {
	int n = ((CxVec3Array *)_dipoleCache[0])->GetSize() - 2;
	float step = (float)m_iShowMolCount / 20.0f;

	if (n < _correlationDepth)
	{
		eprintf("\nError: Autocorrelation depth is %d, but only %d timesteps evaluated.\n",_correlationDepth,n);
		eprintf("       Reduce depth or increase trajectory length.\n\n");
		abort();
	}

	mprintf("    Deriving dipole moments...\n");
	mprintf(WHITE, "     [");
	int i, j, k, l;
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmodf(i, step) < 1.0f)
			mprintf(WHITE, "#");
		for(j = 0; j < n; j++) {
			((CxVec3Array *)_dipoleCache[i])->GetAt(j) = 0.5f * (((CxVec3Array *)_dipoleCache[i])->GetAt(j+2) - ((CxVec3Array *)_dipoleCache[i])->GetAt(j)) / g_fTimestepLength;
		}
	}
	mprintf(WHITE, "]\n");
	
	mprintf("    Calculating autocorrelation...\n");
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
	temp->SetSize(n);
	CxFloatArray *temp2;
	try { temp2 = new CxFloatArray(); } catch(...) { temp2 = NULL; }
	if(temp2 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp2->SetSize(n);
	
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmodf(i, step) < 1.0f)
			mprintf(WHITE, "#");
		for(j = 0; j < 3; j++) {
			for(k = 0; k < n; k++) {
				temp->GetAt(k) = ((CxVec3Array *)_dipoleCache[i])->GetAt(k)[j];
			}
			ac->AutoCorrelate(temp, temp2);
			for(k = 0; k < _correlationDepth; k++) {
				acf->GetAt(k) += temp2->GetAt(k);
			}
		}
	}
	mprintf(WHITE, "]\n");
// 	if(m_iShowMol == -1) {
// 		for(i = 0; i < _correlationDepth; i++) {
// 			acf->GetAt(i) /= 3.0f;
// 		}
// 	} else {
// 		for(i = 0; i < _correlationDepth; i++) {
// 			acf->GetAt(i) /= 3.0f * m_iShowMolCount;
// 		}
// 	}
// The 3.0f is included in the final normalization
	if(m_iShowMol != -1) {
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) /= (float)m_iShowMolCount;
		}
	}
	
	delete ac;
	
	CxFloatArray *ccf = NULL;
	if(_includeCross) {
		mprintf("    Calculating cross-correlation...\n");
		try { ccf = new CxFloatArray(); } catch(...) { ccf = NULL; }
		if(ccf == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		ccf->SetSize(_correlationDepth);
		for(i = 0; i < _correlationDepth; i++) {
			ccf->GetAt(i) = 0.0f;
		}
		CCrossCorrelation *cc;
		try { cc = new CCrossCorrelation(); } catch(...) { cc = NULL; }
		if(cc == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		cc->Init(n, _correlationDepth, g_bACFFFT);
		CxFloatArray *temp3;
		try { temp3 = new CxFloatArray(); } catch(...) { temp3 = NULL; }
		if(temp3 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		temp3->SetSize(n);
		temp->SetSize(n);
		temp2->SetSize(n);
		
		mprintf(WHITE, "     [");
		step = (float)m_iShowMolCount * (m_iShowMolCount - 1) / 20.0f;
		int c = 0;
		for(i = 0; i < m_iShowMolCount; i++) {
			for(j = 0; j < m_iShowMolCount; j++) {
				if (i == j)
					continue;
				if(fmodf((float)c++, step) < 1.0f)
					mprintf(WHITE, "#");
				for(k = 0; k < 3; k++) {
					for(l = 0; l < n; l++) {
						temp->GetAt(l) = ((CxVec3Array *)_dipoleCache[i])->GetAt(l)[k];
						temp2->GetAt(l) = ((CxVec3Array *)_dipoleCache[j])->GetAt(l)[k];
					}
					cc->CrossCorrelate(temp, temp2, temp3);
					for(l = 0; l < _correlationDepth; l++) {
						ccf->GetAt(l) += temp3->GetAt(l);
					}
				}
			}
		}
		mprintf(WHITE, "]\n");
// 		for(i = 0; i < _correlationDepth; i++) {
// 			ccf->GetAt(i) /= 3.0f;
// 		}
		
		delete cc;
		delete temp3;
	}
	
	delete temp2;
	
	if(_saveACF) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "ir_acf_%s.csv", _name);
#else
		sprintf(filename, "ir_acf_%s.csv", _name);
#endif
		mprintf("    Saving autocorrelation function as %s...\n", filename);
		FILE *acfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(acfFile, "%.2f; %.10G\n", i * g_fTimestepLength, acf->GetAt(i));
		}
		fclose(acfFile);
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
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "ir_acf_windowed_%s.csv", _name);
#else
		sprintf(filename, "ir_acf_windowed_%s.csv", _name);
#endif
		mprintf("    Saving windowed autocorrelation function as %s...\n", filename);
		FILE *acfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(acfFile, "%.2f; %.10G\n", i * g_fTimestepLength, temp->GetAt(i));
		}
		fclose(acfFile);
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
	
	mprintf("    Fourier transforming autocorrelation function...\n");
	CFFT *fft;
	try { fft = new CFFT(); } catch(...) { fft = NULL; }
	if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	fft->PrepareFFT_C2C(temp->GetSize());
	for(i = 0; i < temp->GetSize(); i++) {
		fft->m_pInput[2*i] = temp->GetAt(i);
		fft->m_pInput[2*i+1] = 0.0f;
	}
	fft->DoFFT();
	
	CxFloatArray *spectrum;
	try { spectrum = new CxFloatArray(); } catch(...) { spectrum = NULL; }
	if(spectrum == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	spectrum->SetSize(_specSize);
	for(i = 0; i < _specSize; i++) {
		spectrum->GetAt(i) = 3047.2310f * fft->m_pOutput[2*i] * g_fTimestepLength; // Output in K*cm*km/mol
	}
	
	if(_finiteDifferenceCorrection) {
		float f = _specResolution * g_fTimestepLength * 1.883652e-4f;
		for(i = 1; i < _specSize; i++) {
			spectrum->GetAt(i) *= powf(f * i / sinf(f * i), 2.0f); // Divide by sinc function to correct finite difference derivation
		}
	}
	
	
	if(_includeCross) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "ir_spectrum_auto_%s.csv", _name);
#else
		sprintf(filename, "ir_spectrum_auto_%s.csv", _name);
#endif
		mprintf("    Saving autocorrelation spectrum as %s...\n", filename);
		FILE *specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");

		double integral = 0.0;
		for(i = 0; i < _specSize; i++) {
			integral += (double)spectrum->GetAt(i) * _specResolution;
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum->GetAt(i), integral);
		}
		fclose(specFile);
		
		if(_saveACF) {
			char filename[BUF_SIZE];
#ifdef TARGET_LINUX
			snprintf(filename, BUF_SIZE, "ir_ccf_%s.csv", _name);
#else
			sprintf(filename, "ir_ccf_%s.csv", _name);
#endif
			mprintf("    Saving cross-correlation function as %s...\n", filename);
			FILE *ccfFile = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(ccfFile, "%.2f; %.10G\n", i * g_fTimestepLength, ccf->GetAt(i));
			}
			fclose(ccfFile);
		}
		
		temp->CopyFrom(ccf);
		
		if(_windowFunction == 1) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= powf(cosf((float)i / temp->GetSize() / 2.0f * Pi), 2.0f);
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
			char filename[BUF_SIZE];
#ifdef TARGET_LINUX
			snprintf(filename, BUF_SIZE, "ir_ccf_windowed_%s.csv", _name);
#else
			sprintf(filename, "ir_ccf_windowed_%s.csv", _name);
#endif
			mprintf("    Saving windowed cross-correlation function as %s...\n", filename);
			FILE *ccfFile = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(ccfFile, "%.2f; %.10G\n", i * g_fTimestepLength, ccf->GetAt(i));
			}
			fclose(ccfFile);
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
		
		mprintf("    Fourier transforming cross-correlation function...\n");
		fft->PrepareFFT_C2C(temp->GetSize());
		for(i = 0; i < temp->GetSize(); i++) {
			fft->m_pInput[2*i] = temp->GetAt(i);
			fft->m_pInput[2*i+1] = 0.0f;
		}
		fft->DoFFT();
		
		CxFloatArray *ccSpectrum;
		try { ccSpectrum = new CxFloatArray(); } catch(...) { ccSpectrum = NULL; }
		if(ccSpectrum == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		ccSpectrum->SetSize(_specSize);
		for(i = 0; i < _specSize; i++) {
			ccSpectrum->GetAt(i) = 3047.2310f * fft->m_pOutput[2*i] * g_fTimestepLength;
		}
		
		
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "ir_spectrum_cross_%s.csv", _name);
#else
		sprintf(filename, "ir_spectrum_cross_%s.csv", _name);
#endif
		mprintf("    Saving cross-correlation spectrum as %s...\n", filename);
		specFile = OpenFileWrite(filename, false);

			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");

		integral = 0.0;
		for(i = 0; i < _specSize; i++) {
			integral += ccSpectrum->GetAt(i) * _specResolution;
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, ccSpectrum->GetAt(i), integral);
		}
		fclose(specFile);
		
		for(i = 0; i < _specSize; i++) {
			spectrum->GetAt(i) += ccSpectrum->GetAt(i);
		}
		
		delete ccSpectrum;
	}
	
	char filename[BUF_SIZE];
#ifdef TARGET_LINUX
	snprintf(filename, BUF_SIZE, "ir_spectrum_%s.csv", _name);
#else
	sprintf(filename, "ir_spectrum_%s.csv", _name);
#endif
	mprintf("    Saving spectrum as %s...\n", filename);
	FILE *specFile = OpenFileWrite(filename, false);

		fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");

	double integral = 0.0;
	for(i = 0; i < _specSize; i++) {
		integral += spectrum->GetAt(i) * _specResolution;
			fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum->GetAt(i), integral);
	}
	fclose(specFile);
	
	delete acf;
	if(_includeCross) {
		delete ccf;
	}
	delete temp;
	delete fft;
	delete spectrum;
	
	for(i = 0; i < _dipoleCache.GetSize(); i++)
		delete (CxVec3Array *)_dipoleCache[i];
}


bool gatherIR() {
	g_bDipole = true;
	ParseDipole();
	
	while(true) {
		mprintf(YELLOW, "\n>>> IR Observation %d >>>\n\n", g_IRObserv.GetSize() + 1);
		
		CIRObservation *obs;
		try { obs = new CIRObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CIRObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_IRObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of IR Observation %d <<<\n\n", g_IRObserv.GetSize());
		
		if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	
	if(AskYesNo("\n    Compute IR spectrum of whole system (y/n) [no] ", false)) {
		mprintf(YELLOW, "\n>>> Global IR Observation >>>\n\n");
		
		CIRObservation *obs;
		try { obs = new CIRObservation(true); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CIRObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_IRObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of Global IR Observation <<<\n\n");
	}
	
	return true;
}

bool initializeIR() {
	int i;
	for(i = 0; i < g_IRObserv.GetSize(); i++) {
		((CIRObservation *)g_IRObserv[i])->initialize();
	}
	return true;
}

void processIR(CTimeStep* ts) {
	(void)ts;
	int i;
	for(i = 0; i < g_IRObserv.GetSize(); i++) {
		((CIRObservation *)g_IRObserv[i])->process();
	}
}

void finalizeIR() {
	int i;
	for(i = 0; i < g_IRObserv.GetSize(); i++) {
		mprintf(YELLOW, "\n>>> IR Observation %d >>>\n\n", i+1);
		((CIRObservation *)g_IRObserv[i])->finalize();
		delete (CIRObservation *)g_IRObserv[i];
		mprintf(YELLOW, "\n<<< End of IR Observation %d <<<\n\n", i+1);
	}
}

static FILE *g_dipoleRestartFile = NULL;

bool gatherDipoleRestart() {
	g_bDipole = true;
	ParseDipole();
	
	mprintf(YELLOW, "\n>>> Dipole Restart File >>>\n\n");
	
	CxString filename;
	AskString("    Enter name of dipole restart file to write: [dipole.restart] ", &filename, "dipole.restart");
	
	g_dipoleRestartFile = OpenFileWrite((const char *)filename, false);
	int numMolecules = g_oaSingleMolecules.GetSize();
	fwrite(&numMolecules, sizeof(int), 1, g_dipoleRestartFile);
	
	mprintf(YELLOW, "\n<<< End of Dipole Restart File <<<\n\n");
	return true;
}

bool initializeDipoleRestart() {
	return true;
}

void processDipoleRestart(CTimeStep *ts) {
	(void)ts;
	int i;
	for (i = 0; i < g_oaSingleMolecules.GetSize(); i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
		int j;
		for (j = 0; j < 3; j++) {
			float val = sm->m_vDipole[j];
			fwrite(&val, sizeof(float), 1, g_dipoleRestartFile);
		}
	}
	fflush(g_dipoleRestartFile);
}

void finalizeDipoleRestart() {
	if (g_dipoleRestartFile != NULL)
		fclose(g_dipoleRestartFile);
}

static CxObArray g_VCDObserv;

CVCDObservation::CVCDObservation(bool global) {
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
		}
		m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
		_name = new char[strlen(((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName) + 1];
		strcpy(_name, ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		mprintf("\n");
	}
	
	if(g_iTrajSteps != -1) {
		_correlationDepth = 0.75 * g_iTrajSteps;
		if(_correlationDepth > 4096)
			_correlationDepth = 4096;
		if(g_fTimestepLength > 1.0f)
			_correlationDepth = 2048;
		if(g_fTimestepLength > 2.0f)
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
		mprintf("    Using cos^2 window function.\n");
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
		mprintf("    Inserting %d zeros for zero padding.\n",_zeroPadding);
	}
	
	float possibleRange = 33356.41f / g_fTimestepLength / 2.0f;
	_specResolution = possibleRange / (_correlationDepth + _zeroPadding);
	mprintf("    This results in a spectral resolution of %.2f cm^-1.\n", _specResolution);
	mprintf("\n    A time step length of %.2f fs allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, possibleRange);
	float specLimit = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ", 0, possibleRange, (possibleRange < 5000.0f) ? possibleRange : 5000.0f, (possibleRange < 5000.0f) ? possibleRange : 5000.0f);
	_specSize = specLimit / _specResolution;
	mprintf("\n");
	
	if(g_bAdvanced2) {
		_smoothWidth = AskUnsignedInteger("    Window width of median filter to smooth magnetic dipole moments: [0] ", 0);
		mprintf("\n");
	} else {
		_smoothWidth = 0;
	}
	
	if(g_bAdvanced2) {
		_saveMoments = AskYesNo("    Save electric and magnetic moments (y/n)? [no] ", false);
	} else {
		_saveMoments = false;
	}
	
	if(g_bAdvanced2) {
		_saveACF = AskYesNo("    Save cross-correlation function (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_saveACF = false;
	}
}

CVCDObservation::~CVCDObservation() {
	delete[] _name;
}

void CVCDObservation::initialize() {
	int n;
	if(g_iTrajSteps != -1)
		n = 1.1 * g_iTrajSteps / g_iStride;
	else
		n = 10000;
	
	mprintf("    Moment cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * 2.0 * n * sizeof(CxVector3)));
	int i;
	for(i = 0; i < m_iShowMolCount; i++) {
		CxVec3Array *a, *b;
		try { a = new CxVec3Array(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		try { b = new CxVec3Array(); } catch(...) { b = NULL; }
		if(b == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetMaxSize(n);
		a->SetGrow(n / 10);
		b->SetMaxSize(n);
		b->SetGrow(n / 10);
		_electricDipoleCache.Add(a);
		_magneticDipoleCache.Add(b);
	}
}

void CVCDObservation::process() {
	int i;
	for(i = 0; i < m_iShowMolCount; i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
		((CxVec3Array *)_electricDipoleCache[i])->Add(sm->m_vDipole);
		((CxVec3Array *)_magneticDipoleCache[i])->Add(sm->m_vMagneticDipole);
	}
}

static int compare_float(const void *f1, const void *f2) {
	if(*(float *)f1 < *(float *)f2)
		return -1;
	if(*(float *)f1 > *(float *)f2)
		return 1;
	return 0;
}

void CVCDObservation::finalize() {
	int i, j, k, l;
	if(_smoothWidth > 0) {
		mprintf("    Smoothing magnetic dipole moments...\n");
		CxVec3Array *tempMoments;
		try { tempMoments = new CxVec3Array(); } catch(...) { tempMoments = NULL; }
		if(tempMoments == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		float *sortArray;
		try { sortArray = new float[2*_smoothWidth+1]; } catch(...) { sortArray = NULL; }
		if(sortArray == NULL) NewException((double)_smoothWidth * sizeof(float), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		for(i = 0; i < m_iShowMolCount; i++) {
			tempMoments->SetSize(((CxVec3Array *)_magneticDipoleCache[i])->GetSize());
			for(j = 0; j < tempMoments->GetSize(); j++) {
				for(k = 0; k < 3; k++) {
					int index = 0;
					for(l = j - _smoothWidth; l <= j + _smoothWidth; l++) {
						if(l < 0)
							continue;
						if(l >= tempMoments->GetSize())
							break;
						sortArray[index] = ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(l)[k];
						index++;
					}
					qsort(sortArray, index, sizeof(float), compare_float);
					if(index % 2 == 0) {
						tempMoments->GetAt(j)[k] = (sortArray[index / 2] + sortArray[index / 2 - 1]) / 2.0f;
					} else {
						tempMoments->GetAt(j)[k] = sortArray[index / 2];
					}
				}
			}
			
// 			for(j = 0; j < 9; j++) {
// 				tempMoments->GetAt(j) = ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j);
// 			}
// 			for(j = 9; j < tempMoments->GetSize() - 9; j++) {
// 				int k;
// 				for(k = 0; k < 3; k++) {
// 					float array[19];
// 					int l;
// 					for(l = j - 9; l <= j + 9; l++) {
// 						array[l - j + 9] = ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(l)[k];
// 					}
// 					qsort(array, 19, sizeof(float), compare_float);
// 					tempMoments->GetAt(j)[k] = array[9];
// 				}
// 			}
// 			for(j = tempMoments->GetSize() - 9; j < tempMoments->GetSize(); j++) {
// 				tempMoments->GetAt(j) = ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j);
// 			}
// 		for(j = 0; j < tempMoments->GetSize(); j++) {
// 			tempMoments->GetAt(j) = CxVector3(0.0f, 0.0f, 0.0f);
// 			int k;
// 			int count = 0;
// 			for(k = j - 10; k <= j + 10; k++) {
// 				if(k < 0)
// 					continue;
// 				if(k >= tempMoments->GetSize())
// 					break;
// 				tempMoments->GetAt(j) += ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(k);
// 				count++;
// 			}
// 			tempMoments->GetAt(j) /= (float)count;
// 		}
// 		for(j = 0; j < 12; j++) {
// 			tempMoments->GetAt(j) = ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j);
// 		}
// 		for(j = 12; j < tempMoments->GetSize() - 12; j++) {
// // 			tempMoments->GetAt(j) = (-3.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-2) + 12.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-1) + 17.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j) + 12.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+1) - 3.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+2)) / 35.0f;
// // 			tempMoments->GetAt(j) = (-21.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-4) + 14.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-3) + 39.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-2) + 54.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-1) + 59.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j) + 54.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+1) + 39.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+2) + 14.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+3) - 21.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+4)) / 231.0f;
// tempMoments->GetAt(j) = (-253.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-12) - 138.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-11) - 33.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-10) + 62.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-9) + 147.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-8) + 222.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-7) + 287.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-6) + 322.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-5) + 387.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-4) + 422.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-3) + 447.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-2) + 462.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j-1) + 467.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j) + 462.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+1) + 447.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+2) +
// 422.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+3) + 387.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+4) + 322.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+5) + 287.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+6) + 222.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+7) + 147.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+8) + 62.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+9) - 33.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+10) - 138.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+11) - 253.0f * ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+12)) / 5175.0f;
// 		}
// 		for(j = tempMoments->GetSize() - 12; j < tempMoments->GetSize(); j++) {
// 			tempMoments->GetAt(j) = ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j);
// 		}
			((CxVec3Array *)_magneticDipoleCache[i])->CopyFrom(tempMoments);
		}
		delete tempMoments;
		delete[] sortArray;
	}
	
	if(_saveMoments) {
		mprintf("    Saving electric and magnetic dipole moments...\n");
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "vcd_moments_%s.csv", _name);
#else
		sprintf(filename, "vcd_moments_%s.csv", _name);
#endif
		FILE *momentsFile = OpenFileWrite(filename, false);
		for(i = 0; i < ((CxVec3Array *)_magneticDipoleCache[0])->GetSize(); i++) {
			fprintf(momentsFile, "%.2f;", g_fTimestepLength * (i + 1));
			for(j = 0; j < m_iShowMolCount; j++) {
				fprintf(momentsFile, " %.10G; %.10G; %.10G; %.10G; %.10G; %.10G;", ((CxVec3Array *)_electricDipoleCache[j])->GetAt(i)[0], ((CxVec3Array *)_electricDipoleCache[j])->GetAt(i)[1], ((CxVec3Array *)_electricDipoleCache[j])->GetAt(i)[2], ((CxVec3Array *)_magneticDipoleCache[j])->GetAt(i)[0], ((CxVec3Array *)_magneticDipoleCache[j])->GetAt(i)[1], ((CxVec3Array *)_magneticDipoleCache[j])->GetAt(i)[2]);
			}
			fprintf(momentsFile, "\n");
		}
		fclose(momentsFile);
	}
	
	mprintf("    Deriving electric dipole moments...\n");
	int n = ((CxVec3Array *)_electricDipoleCache[0])->GetSize() - 2;
	float step = (float)m_iShowMolCount / 20.0f;
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmodf(i, step) < 1.0f)
			mprintf(WHITE, "#");
		for(j = 0; j < n; j++) {
			((CxVec3Array *)_electricDipoleCache[i])->GetAt(j) = 0.5f * (((CxVec3Array *)_electricDipoleCache[i])->GetAt(j+2) - ((CxVec3Array *)_electricDipoleCache[i])->GetAt(j)) / g_fTimestepLength;
		}
	}
	mprintf(WHITE, "]\n");
	mprintf("    Deriving magnetic dipole moments...\n");
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmodf(i, step) < 1.0f)
			mprintf(WHITE, "#");
		for(j = 0; j < n; j++) {
			((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j) = 0.5f * (((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j+2) - ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(j)) / g_fTimestepLength;
		}
	}
	mprintf(WHITE, "]\n");
	
	if(_saveMoments) {
		mprintf("    Saving electric and magnetic dipole moment derivatives...\n");
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "vcd_moments_deriv_%s.csv", _name);
#else
		sprintf(filename, "vcd_moments_deriv_%s.csv", _name);
#endif
		FILE *momentsFile = OpenFileWrite(filename, false);
		for(i = 0; i < n; i++) {
			fprintf(momentsFile, "%.2f;", g_fTimestepLength * (i + 2));
			for(j = 0; j < m_iShowMolCount; j++) {
				fprintf(momentsFile, " %.10G; %.10G; %.10G; %.10G; %.10G; %.10G;", ((CxVec3Array *)_electricDipoleCache[j])->GetAt(i)[0], ((CxVec3Array *)_electricDipoleCache[j])->GetAt(i)[1], ((CxVec3Array *)_electricDipoleCache[j])->GetAt(i)[2], ((CxVec3Array *)_magneticDipoleCache[j])->GetAt(i)[0], ((CxVec3Array *)_magneticDipoleCache[j])->GetAt(i)[1], ((CxVec3Array *)_magneticDipoleCache[j])->GetAt(i)[2]);
			}
			fprintf(momentsFile, "\n");
		}
		fclose(momentsFile);
	}
	
	mprintf("    Calculating cross-correlation...\n");
	CxFloatArray *ccf;
	try { ccf = new CxFloatArray(); } catch(...) { ccf = NULL; }
	if(ccf == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	ccf->SetSize(_correlationDepth);
	for(i = 0; i < _correlationDepth; i++) {
		ccf->GetAt(i) = 0.0f;
	}
	CCrossCorrelation *cc;
	try { cc = new CCrossCorrelation(); } catch(...) { cc = NULL; }
	if(cc == NULL) NewException((double)sizeof(CCrossCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	cc->Init(n, _correlationDepth, g_bACFFFT);
	
	CxFloatArray *temp;
	try { temp = new CxFloatArray(); } catch(...) { temp = NULL; }
	if(temp == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp->SetSize(n);
	CxFloatArray *temp2;
	try { temp2 = new CxFloatArray(); } catch(...) { temp2 = NULL; }
	if(temp2 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp2->SetSize(n);
	CxFloatArray *temp3;
	try { temp3 = new CxFloatArray(); } catch(...) { temp3 = NULL; }
	if(temp3 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp3->SetSize(n);
	
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmodf(i, step) < 1.0f)
			mprintf(WHITE, "#");
		for(j = 0; j < 3; j++) {
			for(k = 0; k < n; k++) {
				temp->GetAt(k) = ((CxVec3Array *)_electricDipoleCache[i])->GetAt(k)[j];
				temp2->GetAt(k) = ((CxVec3Array *)_magneticDipoleCache[i])->GetAt(k)[j];
			}
			cc->CrossCorrelate(temp, temp2, temp3);
			for(k = 0; k < _correlationDepth; k++) {
				ccf->GetAt(k) += temp3->GetAt(k);
			}
			cc->CrossCorrelate(temp2, temp, temp3);
			for(k = 0; k < _correlationDepth; k++) {
				ccf->GetAt(k) -= temp3->GetAt(k);
			}
		}
	}
	mprintf(WHITE, "]\n");
// The 3.0f is included in the final normalization
// 	for(i = 0; i < _correlationDepth; i++) {
// 		ccf->GetAt(i) /= 3.0f * g_iSteps * m_iShowMolCount;
// 	}
	if(m_iShowMol != -1) {
		for(i = 0; i < _correlationDepth; i++) {
			ccf->GetAt(i) /= (float)m_iShowMolCount;
		}
	}
	
	delete cc;
	delete temp2;
	delete temp3;
	
	if(_saveACF) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "vcd_ccf_%s.csv", _name);
#else
		sprintf(filename, "vcd_ccf_%s.csv", _name);
#endif
		mprintf("    Saving cross-correlation function as %s...\n", filename);
		FILE *ccfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(ccfFile, "%.2f; %.10G\n", i * g_fTimestepLength, ccf->GetAt(i));
		}
		fclose(ccfFile);
	}
	
	temp->CopyFrom(ccf);
	
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
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "vcd_ccf_windowed_%s.csv", _name);
#else
		sprintf(filename, "vcd_ccf_windowed_%s.csv", _name);
#endif
		mprintf("    Saving windowed cross-correlation function as %s...\n", filename);
		FILE *ccfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(ccfFile, "%.2f; %.10G\n", i * g_fTimestepLength, temp->GetAt(i));
		}
		fclose(ccfFile);
	}
	
	if(_zeroPadding > 0) {
		for(i = 0; i < _zeroPadding; i++) {
			temp->Add(0.0f);
		}
	}
	
	int oldSize = temp->GetSize();
	temp->SetSize(2 * oldSize);
	for(i = 1; i < oldSize; i++) {
		temp->GetAt(oldSize + i) = -temp->GetAt(oldSize - i);
	}
	temp->GetAt(oldSize) = 0.0f;
	
	mprintf("    Fourier transforming cross-correlation function...\n");
	CFFT *fft;
	try { fft = new CFFT(); } catch(...) { fft = NULL; }
	if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	fft->PrepareFFT_C2C(temp->GetSize());
	for(i = 0; i < temp->GetSize(); i++) {
		fft->m_pInput[2*i] = temp->GetAt(i);
		fft->m_pInput[2*i+1] = 0.0f;
	}
	fft->DoFFT();

	CxFloatArray *spectrum;
	try { spectrum = new CxFloatArray(); } catch(...) { spectrum = NULL; }
	if(spectrum == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	spectrum->SetSize(_specSize);
	for(i = 0; i < _specSize; i++) {
		spectrum->GetAt(i) = 28.260058f * fft->m_pOutput[2*i+1] * g_fTimestepLength; // Output in K*cm*km/mol
	}
	
	if(_finiteDifferenceCorrection) {
		float f = _specResolution * g_fTimestepLength * 1.883652e-4f;
		for(i = 1; i < _specSize; i++) {
			spectrum->GetAt(i) *= powf(f * i / sinf(f * i), 2.0f); // Divide by sinc function to correct finite difference derivation
		}
	}
	
	char filename[BUF_SIZE];
#ifdef TARGET_LINUX
	snprintf(filename, BUF_SIZE, "vcd_spectrum_%s.csv", _name);
#else
	sprintf(filename, "vcd_spectrum_%s.csv", _name);
#endif
	mprintf("    Saving spectrum as %s...\n", filename);
	FILE *specFile = OpenFileWrite(filename, false);
	fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");
	double integral = 0.0;
	for(i = 0; i < _specSize; i++) {
		integral += spectrum->GetAt(i) * _specResolution;
		fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum->GetAt(i), integral);
	}
	fclose(specFile);
	
	delete ccf;
	delete temp;
	delete fft;
	delete spectrum;
	
	for(i = 0; i < _electricDipoleCache.GetSize(); i++) {
		delete (CxVec3Array *)_electricDipoleCache[i];
		delete (CxVec3Array *)_magneticDipoleCache[i];
	}
}

bool gatherVCD() {
	parseMagneticDipole();
	ParseDipole();
	
	while(true) {
		mprintf(YELLOW, "\n>>> VCD Observation %d >>>\n\n", g_VCDObserv.GetSize() + 1);
		
		CVCDObservation *obs;
		try { obs = new CVCDObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CVCDObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_VCDObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of VCD Observation %d <<<\n\n", g_VCDObserv.GetSize());
		
		if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	
	if(AskYesNo("    Compute VCD spectrum of whole system (y/n) [no] ", false)) {
		mprintf(YELLOW, "\n>>> Global VCD Observation >>>\n\n");
		
		CVCDObservation *obs;
		try { obs = new CVCDObservation(true); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CVCDObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_VCDObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of Global VCD Observation <<<\n\n");
	}
	
	return true;
}

bool initializeVCD() {
	int i;
	for(i = 0; i < g_VCDObserv.GetSize(); i++) {
		((CVCDObservation *)g_VCDObserv[i])->initialize();
	}
	return true;
}

void processVCD(CTimeStep *ts) {
	(void)ts;
// 	ts->CalcMagneticDipoles();
	int i;
	for(i = 0; i < g_VCDObserv.GetSize(); i++) {
		((CVCDObservation *)g_VCDObserv[i])->process();
	}
}

void finalizeVCD() {
	int i;
	for(i = 0; i < g_VCDObserv.GetSize(); i++) {
		((CVCDObservation *)g_VCDObserv[i])->finalize();
	}
}

static FILE *g_magneticDipoleRestartFile = NULL;

bool gatherMagneticDipoleRestart() {
	parseMagneticDipole();
	ParseDipole();
	
	mprintf(YELLOW, "\n>>> Magnetic Moment Restart File >>>\n\n");
	
	CxString filename;
	AskString("    Enter name of magnetic moment restart file to write: [magnetic.restart] ", &filename, "magnetic.restart");
	
	g_magneticDipoleRestartFile = OpenFileWrite((const char *)filename, false);
	int numMolecules = g_oaSingleMolecules.GetSize();
	fwrite(&numMolecules, sizeof(int), 1, g_magneticDipoleRestartFile);
	
	mprintf(YELLOW, "\n<<< End of Magnetic Moment Restart File <<<\n\n");
	return true;
}

bool initializeMagneticDipoleRestart() {
	return true;
}

void processMagneticDipoleRestart(CTimeStep *ts) {
	(void)ts;
	int i;
	for (i = 0; i < g_oaSingleMolecules.GetSize(); i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
		int j;
		for (j = 0; j < 3; j++) {
			float val = sm->m_vMagneticDipole[j];
			fwrite(&val, sizeof(float), 1, g_magneticDipoleRestartFile);
		}
	}
	fflush(g_magneticDipoleRestartFile);
}

void finalizeMagneticDipoleRestart() {
	if (g_magneticDipoleRestartFile != NULL)
		fclose(g_magneticDipoleRestartFile);
}

static FILE *g_sortWannierTrajFile;
static int g_sortWannierCount;
static CxVec3Array g_sortWannierHistory;
static CxIntArray g_sortWannierHistoryIndex;
static CxFloatArray g_sortWannierDistanceMatrix;
static CxIntArray g_sortWannierKMStarCol;
static CxIntArray g_sortWannierKMPrimeRow;
static CxFloatArray g_sortWannierKMMinCol;
static CxFloatArray g_sortWannierKMMinRow;

bool gatherSortWannier() {
	g_bKeepOriginalCoords = true;
	g_bWannier = true;
	int watom = -1;
	int i;
	for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
		if(((CAtom *)g_oaAtoms[i])->m_bExclude) {
			watom = i;
			break;
		}
	}
	if(watom == -1)
		eprintf("    Atom label of Wannier centers not found.\n\n");
	else
		mprintf("    Atom type %s is excluded from the system, probably these are the Wannier centers.\n\n", ((CAtom *)g_oaAtoms[watom])->m_sName);
	
	bool ok = false;
	char buf[BUF_SIZE];
	char buf2[BUF_SIZE];
	size_t remaining = BUF_SIZE;
	while(!ok) {
#ifdef TARGET_LINUX
		remaining -= snprintf(buf, remaining, "    Which atom label do the wannier centers have (");
#else
		remaining -= sprintf(buf, "    Which atom label do the wannier centers have (");
#endif
		for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
			if(remaining < 1)
				break;
#ifdef TARGET_LINUX
			size_t length = snprintf(buf2, remaining, "%s", ((CAtom *)g_oaAtoms[i])->m_sName);
#else
			size_t length = sprintf(buf2, "%s", ((CAtom *)g_oaAtoms[i])->m_sName);
#endif
			strncat(buf, buf2, remaining - 1);
			remaining -= length;
			if(i < g_oaAtoms.GetSize() - 2) {
#ifdef TARGET_LINUX
				length = snprintf(buf2, remaining, ", ");
#else
				length = sprintf(buf2, ", ");
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;
			}
		}
#ifdef TARGET_LINUX
		size_t length = snprintf(buf2, remaining, ")? ");
#else
		size_t length = sprintf(buf2, ")? ");
#endif
		strncat(buf, buf2, remaining - 1);
		remaining -= length;
		if(watom != -1) {
#ifdef TARGET_LINUX
			snprintf(buf2, remaining, "[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
#else
			sprintf(buf2, "[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
#endif
			strncat(buf, buf2, remaining - 1);
		}
		
		CxString buf3;
		if(watom == -1)
			AskString_ND(buf, &buf3);
		else
			AskString(buf, &buf3, ((CAtom *)g_oaAtoms[watom])->m_sName);
		
		for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
			if(mystricmp((const char *)buf3, ((CAtom *)g_oaAtoms[i])->m_sName) == 0) {
				g_iWannierAtomType = i;
				ok = true;
				break;
			}
		}
		if(!ok) {
			eprintf("    Wrong input.\n");
			inpprintf("! Wrong input.\n");
		}
	}
	
	g_sortWannierCount = 0;
	for(i = 0; i < g_iGesAtomCount; i++) {
		if(g_baAtomIndex[i] == g_iWannierAtomType)
			g_sortWannierCount++;
	}
	mprintf("\n    There are %d Wannier centers in the first step.\n", g_sortWannierCount);
	
	return true;
}

bool initializeSortWannier() {
	char filename[BUF_SIZE];
#ifdef TARGET_LINUX
	char buf[BUF_SIZE];
	strncpy(buf, g_sInputTraj, BUF_SIZE);
	buf[BUF_SIZE-1] = 0;
	char *p = strrchr(buf, '/');
	if(p == NULL)
		p = buf;
	else
		p++;
	size_t s = strcspn(p, ".");
	if(s > BUF_SIZE - 8)
		s = BUF_SIZE - 8;
	strncpy(filename, p, s);
	filename[s] = 0;
	strcat(filename, "_out.xyz");
#else
	sprintf(filename, "out.xyz");
#endif
	mprintf("    Saving processed trajectory as %s\n\n", filename);
	g_sortWannierTrajFile = OpenFileWrite(filename, false);
	
	mprintf("    Coordinate history: Trying to allocate %s of memory...\n", FormatBytes((double)(g_sortWannierCount * sizeof(int) + g_sortWannierCount * sizeof(CxVector3))));
	g_sortWannierHistory.SetSize(g_sortWannierCount);
	g_sortWannierHistoryIndex.SetSize(g_sortWannierCount);
	
	mprintf("    Hungarian algorithm: Trying to allocate %s of memory...\n", FormatBytes((double)(g_sortWannierCount * g_sortWannierCount * sizeof(float) + 2 * g_sortWannierCount * sizeof(float) + 2 * g_sortWannierCount * sizeof(int))));
	g_sortWannierDistanceMatrix.SetSize(g_sortWannierCount * g_sortWannierCount);
	g_sortWannierKMStarCol.SetSize(g_sortWannierCount);
	g_sortWannierKMPrimeRow.SetSize(g_sortWannierCount);
	g_sortWannierKMMinCol.SetSize(g_sortWannierCount);
	g_sortWannierKMMinRow.SetSize(g_sortWannierCount);
	
	return true;
}

void processSortWannier(CTimeStep *ts) {
	static bool first = true;
	if(first) {
		first = false;
		fprintf(g_sortWannierTrajFile, "%d\n", g_iGesAtomCount);
		if(ts->m_pComment != NULL)
			fprintf(g_sortWannierTrajFile, "%s\n", ts->m_pComment);
		else
			fprintf(g_sortWannierTrajFile, "\n");
		
		int i;
		int pos = 0;
		for(i = 0; i < g_iGesAtomCount; i++) {
			if(g_baAtomIndex[i] == g_iWannierAtomType) {
				g_sortWannierHistory[pos] = ts->m_vaCoords[i];
				g_sortWannierHistoryIndex[pos] = i;
				pos++;
				fprintf(g_sortWannierTrajFile, "%4s %14.10f %14.10f %14.10f\n", ((CAtom *)g_oaAtoms[g_iWannierAtomType])->m_sName, ts->m_vaCoords_Original[i][0] / 100.0f, ts->m_vaCoords_Original[i][1] / 100.0f, ts->m_vaCoords_Original[i][2] / 100.0f);
			} else {
				fprintf(g_sortWannierTrajFile, "%4s %14.10f %14.10f %14.10f\n", ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName, ts->m_vaCoords_Original[i][0] / 100.0f, ts->m_vaCoords_Original[i][1] / 100.0f, ts->m_vaCoords_Original[i][2] / 100.0f);
			}
		}
	} else {
		CxIntArray sortIndex;
		sortIndex.SetSize(g_sortWannierHistory.GetSize());
		
		// Hungarian algorithm as described on http://de.wikipedia.org/wiki/Ungarische_Methode, 07.02.14
		const float thresh = 1.0e-5f;
		int i, j;
		for(i = 0; i < g_sortWannierCount; i++) {
			for(j = 0; j < g_sortWannierCount; j++) {
				CxVector3 dist = g_sortWannierHistory[i] - ts->m_vaCoords[g_sortWannierHistoryIndex[j]];
				if(g_bPeriodicX) {
					while(dist[0] > g_fBoxX / 2.0f) dist[0] -= g_fBoxX;
					while(dist[0] < -g_fBoxX / 2.0f) dist[0] += g_fBoxX;
				}
				if(g_bPeriodicY) {
					while(dist[1] > g_fBoxY / 2.0f) dist[1] -= g_fBoxY;
					while(dist[1] < -g_fBoxY / 2.0f) dist[1] += g_fBoxY;
				}
				if(g_bPeriodicZ) {
					while(dist[2] > g_fBoxZ / 2.0f) dist[2] -= g_fBoxZ;
					while(dist[2] < -g_fBoxZ / 2.0f) dist[2] += g_fBoxZ;
				}
				g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] = dist.GetLengthSqr();
			}
		}
		for(i = 0; i < g_sortWannierCount; i++) {
			g_sortWannierKMStarCol[i] = -1;
		}
		
		for(j = 0; j < g_sortWannierCount; j++) {
			g_sortWannierKMMinCol[j] = g_sortWannierDistanceMatrix[j];
			for(i = 0; i < g_sortWannierCount; i++) {
				if(g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] < g_sortWannierKMMinCol[j]) {
					g_sortWannierKMMinCol[j] = g_sortWannierDistanceMatrix[i*g_sortWannierCount+j];
				}
			}
		}
		for(i = 0; i < g_sortWannierCount; i++) {
			g_sortWannierKMMinRow[i] = g_sortWannierDistanceMatrix[i*g_sortWannierCount] - g_sortWannierKMMinCol[0];
			for(j = 0; j < g_sortWannierCount; j++) {
				if(g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinCol[j] < g_sortWannierKMMinRow[i]) {
					g_sortWannierKMMinRow[i] = g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinCol[j];
				}
			}
		}
		
		while(true) {
			for(i = 0; i < g_sortWannierCount; i++) {
				g_sortWannierKMPrimeRow[i] = -1;
			}
			
			int count = 0;
			for(i = 0; i < g_sortWannierCount; i++) {
				bool star = false;
				for(j = 0; j < g_sortWannierCount; j++) {
					if(g_sortWannierKMStarCol[j] == i) {
						star = true;
						count++;
						break;
					}
				}
				if(star) {
					continue;
				}
				for(j = 0; j < g_sortWannierCount; j++) {
					if(g_sortWannierKMStarCol[j] == -1 && g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinRow[i] - g_sortWannierKMMinCol[j] < thresh) {
						g_sortWannierKMStarCol[j] = i;
						count++;
						break;
					}
				}
			}
			if(count == g_sortWannierCount) {
				break;
			}
			
			int row;
			while(true) {
				float min = 1.0e5f;
				for(i = 0; i < g_sortWannierCount; i++) {
					if(g_sortWannierKMPrimeRow[i] != -1) {
						continue;
					}
					for(j = 0; j < g_sortWannierCount; j++) {
						if((g_sortWannierKMStarCol[j] == -1 || g_sortWannierKMPrimeRow[g_sortWannierKMStarCol[j]] != -1) && g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinRow[i] - g_sortWannierKMMinCol[j] < min) {
							min = g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinRow[i] - g_sortWannierKMMinCol[j];
						}
					}
				}
				
				for(i = 0; i < g_sortWannierCount; i++) {
					if(g_sortWannierKMStarCol[i] == -1 || g_sortWannierKMPrimeRow[g_sortWannierKMStarCol[i]] != -1) {
						g_sortWannierKMMinCol[i] += min;
					}
					if(g_sortWannierKMPrimeRow[i] != -1) {
						g_sortWannierKMMinRow[i] -= min;
					}
				}
				
				row = -1;
				for(i = 0; i < g_sortWannierCount; i++) {
					if(g_sortWannierKMPrimeRow[i] != -1) {
						continue;
					}
					bool found = false;
					for(j = 0; j < g_sortWannierCount; j++) {
						if((g_sortWannierKMStarCol[j] == -1 || g_sortWannierKMPrimeRow[g_sortWannierKMStarCol[j]] != -1) && g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinRow[i] - g_sortWannierKMMinCol[j] < thresh) {
							g_sortWannierKMPrimeRow[i] = j;
							row = i;
							found = true;
							break;
						}
					}
					if(found) {
						break;
					}
				}
				
				bool star = false;
				for(j = 0; j < g_sortWannierCount; j++) {
					if(g_sortWannierKMStarCol[j] == row) {
						star = true;
						break;
					}
				}
				if(!star) {
					break;
				}
			}
			
			while(true) {
				if(g_sortWannierKMStarCol[g_sortWannierKMPrimeRow[row]] == -1) {
					break;
				}
				int tmp = g_sortWannierKMStarCol[g_sortWannierKMPrimeRow[row]];
				g_sortWannierKMStarCol[g_sortWannierKMPrimeRow[row]] = row;
				row = tmp;
			}
			g_sortWannierKMStarCol[g_sortWannierKMPrimeRow[row]] = row;
		}
		
		for(i = 0; i < g_sortWannierHistoryIndex.GetSize(); i++)
			sortIndex[g_sortWannierKMStarCol[i]] = g_sortWannierHistoryIndex[i];
		
		fprintf(g_sortWannierTrajFile, "%d\n", g_iGesAtomCount);
		if(ts->m_pComment != NULL)
			fprintf(g_sortWannierTrajFile, "%s\n", ts->m_pComment);
		else
			fprintf(g_sortWannierTrajFile, "\n");
		
		int pos = 0;
		for(i = 0; i < g_iGesAtomCount; i++) {
			if(g_baAtomIndex[i] == g_iWannierAtomType) {
				g_sortWannierHistory[pos] = ts->m_vaCoords[sortIndex[pos]];
				fprintf(g_sortWannierTrajFile, "%4s %14.10f %14.10f %14.10f\n", ((CAtom *)g_oaAtoms[g_iWannierAtomType])->m_sName, ts->m_vaCoords_Original[sortIndex[pos]][0] / 100.0f, ts->m_vaCoords_Original[sortIndex[pos]][1] / 100.0f, ts->m_vaCoords_Original[sortIndex[pos]][2] / 100.0f);
				pos++;
			} else {
				fprintf(g_sortWannierTrajFile, "%4s %14.10f %14.10f %14.10f\n", ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName, ts->m_vaCoords_Original[i][0] / 100.0f, ts->m_vaCoords_Original[i][1] / 100.0f, ts->m_vaCoords_Original[i][2] / 100.0f);
			}
		}
	}
}

void finalizeSortWannier() {
	fclose(g_sortWannierTrajFile);
}
