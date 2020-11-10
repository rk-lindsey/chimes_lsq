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

#include "pdf.h"

#include "atomgroup.h"
#include "df.h"
#include "globalvar.h"
#include "maintools.h"
#include "xobarray.h"
#include "xvector3.h"
#include "xstring.h"

//#define BUF_SIZE 1024

static CxVector3 g_normalVector;
static CxVector3 g_fixPoint;

static CxObArray g_pdfObserv;

CPDF::CPDF(int showMol) {
//	char buf[256];
	CxString buf;

	m_iShowMol = showMol;
	_showAtomGes = 0;
	
	try { _df = new CDF(); } catch(...) { _df = NULL; }
	if(_df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	mprintf(WHITE, "\n>>> Fixed Plane Density Profile >>>\n\n");
	
	try { _ag = new CAtomGroup(); } catch(...) { _ag = NULL; }
	if(_ag == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	while(true) {
		if(((CMolecule *)g_oaMolecules[m_iShowMol])->m_iAtomGes == 3) {
			mprintf("    %s is only one atom, there is no choice.\n", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
			_ag->Reset();
			_ag->m_pMolecule = (CMolecule *)g_oaMolecules[m_iShowMol];
			_ag->AddAtom(0, 0, false);
			_ag->SortAtoms();
			_ag->BuildName();
		} else {
			mprintf("    Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			inpprintf("! Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			
			myget(&buf);
			if(strlen(buf) == 0) {
				if(!_ag->ParseAtoms((CMolecule *)g_oaMolecules[m_iShowMol], "#2")) {
					eprintf("Weird error.\n");
					continue;
				}
			} else if(!_ag->ParseAtoms((CMolecule *)g_oaMolecules[m_iShowMol], buf)) {
				continue;
			}
		}
		break;
	}
	_showAtomGes += _ag->m_iAtomGes;
	
	ParseDeriv();
	switch(m_iDeriv) {
		case 0:
			_minDist = AskFloat("    Enter the minimal distance of this Density Profile in pm: [%d.0] ", -(float)HalfBox(), -HalfBox());
			_maxDist = AskFloat("    Enter the maximal distance of this Density Profile in pm: [%d.0] ", (float)HalfBox(), HalfBox());
			break;
/*		case 1:
			if(m_bDerivAbs)
				_minDist = AskFloat("    Enter the minimal value of this d1-PDF in pm/ps: [0] ", 0.0f);
			else
				_minDist = AskFloat("    Enter the minimal value of this d1-PDF in pm/ps: [-10.0] ", -10.0f);
			_maxDist = AskFloat("    Enter the maximal value of this d1-RDF in pm/ps: [10.0] ", 10.0f);
			break;
		case 2:
			if(m_bDerivAbs)
				_minDist = AskFloat("    Enter the minimal value of this d2-RDF in pm/ps^2: [0] ", 0.0f);
			else
				_minDist = AskFloat("    Enter the minimal value of this d2-RDF in pm/ps^2: [-10.0] ", -10.0f);
			_maxDist = AskFloat("    Enter the maximal value of this d2-RDF in pm/ps^2: [10.0] ",10.0f);
			break;*/
		default:
			eprintf("Higher derivatives are not implemented.\n");
			abort();
			break;
	}
	m_iResolution = AskUnsignedInteger("    Enter the resolution (bin count) for this Density Profile: [300] ", 300);
	m_iHistogramRes = 0;
	
	if (!g_bBoxNonOrtho) {
		_scaleUniform = AskYesNo("    Scale to uniform density (y) or to nm^(-1) (n)? [no] ", false);
	} else {
		_scaleUniform = false;
	}
	
	buildName();
	m_bSelf = false;
	
	mprintf(WHITE, "\n<<< Enf of Fixed Plane Density Profile <<<\n\n");
}

CPDF::~CPDF() {
	delete _ag;
	delete _df;
	delete[] m_faData;
	delete[] m_sName;
	delete[] m_sShortName;
}

void CPDF::initialize(int showMolCount) {
	mprintf("    Creating Density Profile...\n");
	_df->m_fMinVal = _minDist;
	_df->m_fMaxVal = _maxDist;
	_df->m_iResolution = m_iResolution;
	_df->m_iHistogramRes = m_iHistogramRes;
	_df->SetLabelX("Distance / pm");
	if (_scaleUniform)
		_df->SetLabelY("g(d)");
	else
		_df->SetLabelY("Density / nm^(-1)");
	_df->Create();
	
	try { m_faData = new CxDoubleArray[showMolCount]; } catch(...) { m_faData = NULL; }
	if(m_faData == NULL) NewException((double)showMolCount*sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	// Add timedev here
}

void CPDF::finalize() {
	CxString filename;

	mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n", _df->m_fBinEntries, _df->m_fSkipEntries, ZeroDivide(_df->m_fSkipEntries, _df->m_fBinEntries + _df->m_fSkipEntries) * 100.0f);
	_df->CalcMeanSD();
	mprintf("    Mean value: %10G pm    Standard deviation: %10G pm\n", _df->m_fMean, _df->m_fSD);
	mprintf("    Min. value: %10G pm    Max.value:          %10G pm\n", _df->m_fMinInput, _df->m_fMaxInput);
	
	_df->MultiplyBin(1.0 / g_iSteps);
	_df->Integrate(false, 1.0);
	if (_scaleUniform && !g_bBoxNonOrtho) {
		mprintf("    Scaling Density Profile to uniform density...\n");
		_df->MultiplyBin((g_normalVector[0] * g_normalVector[0] * g_fBoxX + g_normalVector[1] * g_normalVector[1] * g_fBoxY + g_normalVector[2] * g_normalVector[2] * g_fBoxZ) * m_iResolution / (_maxDist-_minDist) / ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize() / _showAtomGes);
	} else {
		mprintf("    Scaling Density Profile to nm^(-1)...\n");
		_df->MultiplyBin(1.0e3 / (_maxDist - _minDist) * m_iResolution);
	}
	
/*	char filename[BUF_SIZE];

#ifdef TARGET_WINDOWS
	_snprintf(filename, BUF_SIZE, "dprof_%s.csv", m_sName);
#elif defined TARGET_LINUX
	snprintf(filename, BUF_SIZE, "dprof_%s.csv", m_sName);
#else
	sprintf(filename, "dprof_%s.csv", m_sName);
#endif*/

	filename.sprintf("dprof_%s.csv", m_sName);

	mprintf("    Saving density profile as \"%s\"...\n", (const char*)filename);
	_df->Write("", filename, "", true);

/*#ifdef TARGET_WINDOWS
	_snprintf(filename, BUF_SIZE, "dprof_%s.agr", m_sName);
#elif defined TARGET_LINUX
	snprintf(filename, BUF_SIZE, "dprof_%s.agr", m_sName);
#else
	sprintf(filename, "dprof_%s.agr", m_sName);
#endif*/

	filename.sprintf("dprof_%s.agr", m_sName);

	mprintf("    Saving density profile AGR as \"%s\"...\n", (const char*)filename);
	_df->WriteAgr("", filename, "", m_sName, true);
}


void CPDF::buildAtomList(CSingleMolecule *sm, CxIntArray *array) {
	int i, j;
	
	array->RemoveAll_KeepSize();
	
	for(i = 0; i < _ag->m_baAtomType.GetSize(); i++) {
		CxIntArray *a = (CxIntArray *)_ag->m_oaAtoms[i];
		for(j = 0; j < a->GetSize(); j++) {
			array->Add(((CxIntArray *)sm->m_oaAtomOffset[_ag->m_baAtomType[i]])->GetAt(a->GetAt(j)));
		}
	}
}


void CPDF::addToDF(float value) {
	_df->AddToBin(value);
}


void CPDF::buildName() {
//	char tmp[BUF_SIZE];
	CxString tmp;
	
	// Check for overflow is missing!!!
/*	tmp[0] = 0;
	if(m_iDeriv != 0) {
		sprintf(tmp, "deriv%d_", m_iDeriv);
	}
	strcat(tmp, "[");
	strcat(tmp, ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
	strcat(tmp, "_");
	strcat(tmp, _ag->m_sName);
	strcat(tmp, "]");*/

	tmp.sprintf("");
	if(m_iDeriv != 0) {
		tmp.sprintf("deriv%d_", m_iDeriv);
	}
	tmp.strcat("[");
	tmp.strcat(((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
	tmp.strcat("_");
	tmp.strcat(_ag->m_sName);
	tmp.strcat("]");
	
	try { m_sShortName = new char[1]; } catch(...) { m_sShortName = NULL; }
	if(m_sShortName == NULL) NewException((double)sizeof(char), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	m_sShortName[0] = 0;
	
	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if(m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	strcpy(m_sName, tmp);
}

CPDFObservation::CPDFObservation() {
	int i;
	CxString buf, buf2;
	
	m_pConditions = NULL;
	m_bTimeDev = false;
	m_bSelf = false;
	m_bOthers = true;
	
	if(g_oaMolecules.GetSize() > 1) {
/*		char buf[BUF_SIZE], buf2[BUF_SIZE];
		size_t remaining = BUF_SIZE;
		
#ifdef TARGET_WINDOWS
		remaining -= _snprintf(buf, remaining, "    Which molecule should be observed (");
#elif defined TARGET_LINUX
		remaining -= snprintf(buf, remaining, "    Which molecule should be observed (");
#else
		remaining -= sprintf(buf, "    Which molecule should be observed (");
#endif*/

		buf.sprintf("    Which molecule should be observed (");
		
		for(i = 0; i < g_oaMolecules.GetSize(); i++) {

/*			if(remaining < 1)
				break;

#ifdef TARGET_WINDOWS
			size_t length = _snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#elif defined TARGET_LINUX
			size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#else
			size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#endif

			strncat(buf, buf2, remaining - 1);
			remaining -= length;*/

			buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
			buf.strcat(buf2);

			if(i < g_oaMolecules.GetSize() - 1) {

/*#ifdef TARGET_WINDOWS
				length = _snprintf(buf2, remaining, ", ");
#elif defined TARGET_LINUX
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
//		strncat(buf, ")? ", remaining - 1);
		buf.strcat(")? ");

		m_iShowMol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
	} else {
		m_iShowMol = 0;
	}
	m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
	m_bObsCertain = false;
	m_bDecompDist = false;
	m_bDecompType = false;
	
	try { _pdf = new CPDF(m_iShowMol); } catch(...) { _pdf = NULL; }
	if(_pdf == NULL) NewException((double)sizeof(CPDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	// Ask for conditions here
	
	// Add observation list here
}

CPDFObservation::~CPDFObservation() {
	delete _pdf;
}

void CPDFObservation::initialize() {
	_pdf->initialize(m_iShowMolCount);
}

void CPDFObservation::process(CTimeStep *ts) {
	int i, j;
	// Process conditions here
	
	for(i = 0; i < m_iShowMolCount; i++) {
		_pdf->m_faData[i].RemoveAll_KeepSize();
	}
	
	for(i = 0; i < m_iShowMolCount; i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
		CxIntArray temp;
		_pdf->buildAtomList(sm, &temp);
		for(j = 0; j < temp.GetSize(); j++) {
			CxVector3 vec = FoldVector(ts->m_vaCoords[temp[j]] - g_fixPoint);
			float val = DotP(vec, g_normalVector);
			_pdf->m_faData[i].Add(val);
		}
	}
	
	for(i = 0; i < m_iShowMolCount; i++) {
		for(j = 0; j < _pdf->m_faData[i].GetSize(); j++) {
			_pdf->addToDF(_pdf->m_faData[i][j]);
		}
	}
}

void CPDFObservation::finalize() {
	mprintf(WHITE, "* Planar Distribution Function\n");
	_pdf->finalize();
}

bool gatherPDF() {
	mprintf(YELLOW, "\n>>> Reference Plane >>>\n\n");
	
	mprintf("    Enter normal vector of reference plane:\n");
	float x, y, z;
	x = AskFloat("    x component [0.0] ", 0.0f);
	y = AskFloat("    y component [0.0] ", 0.0f);
	z = AskFloat("    z component [1.0] ", 1.0f);
	g_normalVector = CxVector3(x, y, z);
	g_normalVector.Normalize();
	mprintf("\n    Normalized normal vector is\n");
	mprintf("        ");
	g_normalVector.Dump();
	mprintf("\n");
	
	mprintf("\n    Enter coordinates of fix point in pm:\n");
	x = AskFloat("    x component [0.0] ", 0.0f);
	y = AskFloat("    y component [0.0] ", 0.0f);
	z = AskFloat("    z component [0.0] ", 0.0f);
	g_fixPoint = CxVector3(x, y, z);
	mprintf("\n    Fix point vector is\n");
	mprintf("        ");
	g_fixPoint.Dump();
	mprintf("\n");
	
	mprintf(YELLOW, "\n<<< End of Reference Plane <<<\n\n");
	
	while(true) {
		mprintf(YELLOW, "\n>>> Density Profile Observation %d >>>\n\n", g_pdfObserv.GetSize() + 1);
		
		CPDFObservation *obs;
		try { obs = new CPDFObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CPDFObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_pdfObserv.Add(obs);
		
		mprintf(YELLOW, "\n<<< End of Density Profile Observation %d <<<\n\n", g_pdfObserv.GetSize());
			
		if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	
	return true;
}

bool initializePDF() {
	int i;
	for(i = 0; i < g_pdfObserv.GetSize(); i++) {
		((CPDFObservation *)g_pdfObserv[i])->initialize();
	}
	return true;
}

void processPDF(CTimeStep *ts) {
	int i;
	for(i = 0; i < g_pdfObserv.GetSize(); i++) {
		((CPDFObservation *)g_pdfObserv[i])->process(ts);
	}
}

void finalizePDF() {
	int i;
	for(i = 0; i < g_pdfObserv.GetSize(); i++) {
		mprintf(YELLOW, "\n>>> Density Profile Observation %d >>>\n\n", i+1);
		((CPDFObservation *)g_pdfObserv[i])->finalize();
		mprintf(YELLOW, "\n<<< End of Density Profile Observation %d <<<\n\n", i+1);
	}
	
	for(i = 0; i < g_pdfObserv.GetSize(); i++) {
		delete (CPDFObservation *)g_pdfObserv[i];
	}
}
