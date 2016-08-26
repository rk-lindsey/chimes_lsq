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

#include "fdf.h"

#include "conversion.h"
#include "globalvar.h"

static CxObArray g_fdfObserv;

CFDFObservation::CFDFObservation() : CObservation() {
	m_pConditions = NULL;
	m_bTimeDev = false;
	m_bSelf = true;
	m_bOthers = false;
	m_iShowMol2 = -1;
	m_bSecondShowMol = false;
	m_bObsCertain = false;
	m_bDecompDist = false;
	m_bDecompType = false;
	
	CxString buf, buf2;
	
	if (g_oaMolecules.GetSize() > 1) {
		buf.sprintf("    Which molecule should be observed (");
		int i;
		for(i = 0; i < g_oaMolecules.GetSize(); i++) {
			buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
			buf.strcat(buf2);
			if(i < g_oaMolecules.GetSize() - 1)
				buf.strcat(", ");
		}
		buf.strcat(")? ");
		m_iShowMol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
	} else {
		m_iShowMol = 0;
		mprintf("    Observing molecule %s.\n\n", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
	}
	m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
	
	while(true) {
		AskString("    Which atom(s) to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ", &buf, "#2");
		if (m_atoms.ParseAtoms((CMolecule *)g_oaMolecules[m_iShowMol], buf))
			break;
	}
	mprintf("\n    Observing %d atoms of molecule %s.\n\n", m_atoms.m_iAtomGes, ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
	
	m_name.sprintf("[%s_%s]", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, m_atoms.m_sName);
	
	m_massWeighting = AskYesNo("    Use forces (y) or omit mass weighting and use accelerations (n)? [yes] ", true);
	mprintf("\n");
	
	if (m_massWeighting) {
		m_fdf.m_fMinVal = AskFloat("    Please enter the minimal force in pN: [0.0] ", 0.0f);
		m_fdf.m_fMaxVal = AskFloat("    Please enter the maximal force in pN: [10000.0] ", 10000.0f);
		m_fdf.m_iResolution = AskInteger("    Please enter the binning resolution: [100] ", 100);
	} else {
		m_fdf.m_fMinVal = AskFloat("    Please enter the minimal acceleration in pm/ps^2: [0.0] ", 0.0f);
		m_fdf.m_fMaxVal = AskFloat("    Please enter the maximal acceleration in pm/ps^2: [100000.0] ", 100000.0f);
		m_fdf.m_iResolution = AskInteger("    Please enter the binning resolution: [100] ", 100);
	}
}

void CFDFObservation::initialize() {
	if (m_massWeighting)
		m_fdf.SetLabelX("Force / pN");
	else
		m_fdf.SetLabelX("Acceleration / pm*ps^-2");
	m_fdf.SetLabelY("Occurrence");
	m_fdf.Create();
	
	m_masses.SetSize(m_atoms.m_iAtomGes);
	if (m_massWeighting) {
		int n = 0;
		int i;
		for(i = 0; i < m_atoms.m_baRealAtomType.GetSize(); i++) {
			CxIntArray *a = (CxIntArray *)m_atoms.m_oaAtoms[i];
			int j;
			for(j = 0; j < a->GetSize(); j++) {
				if (m_atoms.m_baRealAtomType[i] == g_iVirtAtomType)
					m_masses[n] = m_atoms.m_pMolecule->m_fMass;
				else
					m_masses[n] = ((CAtom *)g_oaAtoms[m_atoms.m_baRealAtomType[i]])->m_pElement->m_fMass;
				m_masses[n] *= (float)(1.0E24 * CONST_ATOMIC_MASS_UNIT); // Conversion to obtain final forces in pN
				n++;
			}
		}
	} else {
		int i;
		for (i = 0; i < m_masses.GetSize(); i++) {
			m_masses[i] = 1.0f;
		}
	}
}

void CFDFObservation::process(CTimeStep *ts) {
	int i, j, k;
	for (i = 0; i < m_iShowMolCount; i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
		int n = 0;
		for (j = 0; j < m_atoms.m_baAtomType.GetSize(); j++) {
			CxIntArray *a = (CxIntArray *)m_atoms.m_oaAtoms[j];
			for (k = 0; k < a->GetSize(); k++) {
				m_fdf.AddToBin(ts->m_vaForces[((CxIntArray *)sm->m_oaAtomOffset[m_atoms.m_baAtomType[j]])->GetAt(a->GetAt(k))].GetLength() * m_masses[n]);
				n++;
			}
		}
	}
}

void CFDFObservation::finalize() {
	mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n", m_fdf.m_fBinEntries, m_fdf.m_fSkipEntries, ZeroDivide(m_fdf.m_fSkipEntries, m_fdf.m_fBinEntries + m_fdf.m_fSkipEntries) * 100.0f);
	m_fdf.CalcMeanSD();
	
	CxString unit;
	if (m_massWeighting)
		unit.sprintf("pN");
	else
		unit.sprintf("pm/ps^2");
	mprintf("    Mean value: %10G %s    Standard deviation: %10G %s\n", m_fdf.m_fMean, (const char *)unit, m_fdf.m_fSD, (const char *)unit);
	mprintf("    Min. value: %10G %s    Max.value:          %10G %s\n", m_fdf.m_fMinInput, (const char *)unit, m_fdf.m_fMaxInput, (const char *)unit);
	
	m_fdf.MultiplyBin(1.0 / g_iSteps);
	m_fdf.Integrate(false, 1.0);
	
	CxString filename;
	filename.sprintf("fdf_%s.csv", (const char *)m_name);
	mprintf("    Saving FDF as \"%s\"...\n", (const char *)filename);
	m_fdf.Write("", filename, "", true);
	
	filename.sprintf("fdf_%s.agr", (const char *)m_name);
	mprintf("    Saving FDF AGR as \"%s\"...\n", (const char *)filename);
	m_fdf.WriteAgr("", filename, "", m_name, true);
}

bool gatherFDF() {
	g_bUseForces = true;
	
	while (true) {
		mprintf(YELLOW, "\n>>> FDF Observation %d >>>\n\n", g_fdfObserv.GetSize() + 1);
		
		CFDFObservation *obs;
		try { obs = new CFDFObservation(); } catch(...) { obs = NULL; }
		if (obs == NULL) NewException((double)sizeof(CFDFObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_fdfObserv.Add(obs);
		
		mprintf(YELLOW, "\n<<< End of FDF Observation %d <<<\n\n", g_fdfObserv.GetSize());
		
		if (!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	
	return true;
}

bool initializeFDF() {
	int i;
	for (i = 0; i < g_fdfObserv.GetSize(); i++) {
		((CFDFObservation *)g_fdfObserv[i])->initialize();
	}
	
	return true;
}

void processFDF(CTimeStep *ts) {
	int i;
	for (i = 0; i < g_fdfObserv.GetSize(); i++) {
		((CFDFObservation *)g_fdfObserv[i])->process(ts);
	}
}

void finalizeFDF() {
	int i;
	for (i = 0; i < g_fdfObserv.GetSize(); i++) {
		mprintf(YELLOW, "\n>>> FDF Observation %d >>>\n\n", i + 1);
		((CFDFObservation *)g_fdfObserv[i])->finalize();
		delete (CFDFObservation *)g_fdfObserv[i];
		mprintf(YELLOW, "\n<<< End of FDF Observation %d <<<\n\n", i + 1);
	}
}
