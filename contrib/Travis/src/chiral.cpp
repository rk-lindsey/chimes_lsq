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

#include "chiral.h"

#include "globalvar.h"
#include "maintools.h"
#include "xobarray.h"

//#define BUF_SIZE 4096

static CxObArray g_chiralObserv;
static bool g_fixTraj;
static FILE *g_fixedTrajFile;

CChiralObservation::CChiralObservation() {
	int i, j;
//	char buf[BUF_SIZE], buf2[BUF_SIZE];
//	size_t remaining = BUF_SIZE;
	CxString buf, buf2;
	
	if(g_oaMolecules.GetSize() > 1) {

/*#ifdef TARGET_LINUX
		remaining -= snprintf(buf, remaining, "    Take atoms from which molecule (");
#else
		remaining -= sprintf(buf, "    Take atoms from which molecule (");
#endif*/

		buf.sprintf("    Take atoms from which molecule (");

		for(i = 0; i < g_oaMolecules.GetSize(); i++) {

/*			if(remaining < 1)
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
//		strncat(buf, ")? ", remaining - 1);
		buf.strcat(")? ");

		m_iShowMol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
	} else {
		m_iShowMol = 0;
	}
	m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
	
	while(true) {
//		char buf[BUF_SIZE];
		do {
			mprintf("    Which atom to take from molecule %s? ", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
			inpprintf("! Which atom to take from molecule %s? ", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
			myget(&buf);
			if(strlen(buf) == 0) {
				eprintf("There is no default. Please enter a character string.\n");
				continue;
			}
		} while(!ParseAtom(buf, m_iShowMol, _ty, _rty, _atom));
		
		CMolAtom *ma = NULL;
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]];
		for(i = 0; i < sm->m_oaMolAtoms.GetSize(); i++) {
			ma = (CMolAtom *)sm->m_oaMolAtoms[i];
			if((ma->m_iType == _ty) && (ma->m_iNumber == _atom))
				break;
			ma = NULL;
		}
		if(ma == NULL) {
			eprintf("Weird error.\n");
			abort();
		}
		
		mprintf("    This atom has %d bonds.\n", ma->m_oaBonds.GetSize());
		if(ma->m_oaBonds.GetSize() == 4) {
			mprintf("    The connected atoms are ");
			for(i = 0; i < 4; i++) {
				if(i > 0)
					mprintf(", ");
				mprintf("%s%d", ((CAtom *)g_oaAtoms[sm->m_baAtomIndex[((CMolAtom *)ma->m_oaBonds[i])->m_iType]])->m_sName, ((CMolAtom *)ma->m_oaBonds[i])->m_iNumber+1);
			}
			mprintf(".\n");
			break;
		}
		mprintf(RED, "Please choose an atom with 4 bonds\n");
	}
	
	_molAtoms.SetMaxSize(m_iShowMolCount);
	_orderArrays.SetMaxSize(m_iShowMolCount);
	for(i = 0; i < m_iShowMolCount; i++) {
		CMolAtom *ma = NULL;
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
		for(j = 0; j < sm->m_oaMolAtoms.GetSize(); j++) {
			ma = (CMolAtom *)sm->m_oaMolAtoms[j];
			if((ma->m_iType == _ty) && (ma->m_iNumber == _atom))
				break;
			ma = NULL;
		}
		if(ma == NULL) {
			eprintf("Weird error.\n");
			abort();
		}
		_molAtoms.Add(ma);
		
		CxIntArray *a;
		try { a = new CxIntArray(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetSize(4);
		_orderArrays.Add(a);
		for(j = 0; j < 4; j++)
			a->GetAt(j) = -1;
		double ac;
		int index = 0;
		do {
			ac = 0.0;
			for(j = 0; j < 4; j++)
				if((a->GetAt(j) == -1) && (((CMolAtom *)ma->m_oaBonds[j])->m_fAtomCode > ac))
					ac = ((CMolAtom *)ma->m_oaBonds[j])->m_fAtomCode;
			int minNum = sm->m_oaMolAtoms.GetSize();
			int minIndex;
			while(true) {
				minIndex = 5;
				for(j = 0; j < 4; j++) {
					if((a->GetAt(j) == -1) && (((CMolAtom *)ma->m_oaBonds[j])->m_fAtomCode == ac) && (((CMolAtom *)ma->m_oaBonds[j])->m_iNumber < minNum)) {
						minNum = ((CMolAtom *)ma->m_oaBonds[j])->m_iNumber;
						minIndex = j;
					}
				}
				if(minIndex == 5)
					break;
				a->GetAt(minIndex) = index++;
			}
		} while(ac > 0.0);
	}
	
	mprintf("    The priority of the connected atoms is ");
	for(i = 0; i < 4; i++) {
		int a = ((CxIntArray *)_orderArrays[0])->GetPosition(i);
		if(i > 0)
			mprintf(" > ");
		mprintf("%s%d", ((CAtom *)g_oaAtoms[((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_baAtomIndex[((CMolAtom *)((CMolAtom *)_molAtoms[0])->m_oaBonds[a])->m_iType]])->m_sName, ((CMolAtom *)((CMolAtom *)_molAtoms[0])->m_oaBonds[a])->m_iNumber+1);
	}
	mprintf(".\n");
	
	if(g_fixTraj) {
		_fixTraj = AskYesNo("\n    Adjust chirality for this atom (y/n)? [yes] ", true);
		if(_fixTraj) {
//			char buf[BUF_SIZE], buf2[BUF_SIZE];
//			buf[0] = 0;
//			size_t remaining = BUF_SIZE;
			buf.sprintf("");

			for(i = 0; i < 4; i++) {
				int a = ((CxIntArray *)_orderArrays[0])->GetPosition(i);

/*#ifdef TARGET_LINUX
				size_t length = snprintf(buf2, BUF_SIZE, "%s%d=%d", ((CAtom *)g_oaAtoms[((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_baAtomIndex[((CMolAtom *)((CMolAtom *)_molAtoms[0])->m_oaBonds[a])->m_iType]])->m_sName, ((CMolAtom *)((CMolAtom *)_molAtoms[0])->m_oaBonds[a])->m_iNumber+1, i+1);
#else
				size_t length = sprintf(buf2, "%s%d=%d", ((CAtom *)g_oaAtoms[((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_baAtomIndex[((CMolAtom *)((CMolAtom *)_molAtoms[0])->m_oaBonds[a])->m_iType]])->m_sName, ((CMolAtom *)((CMolAtom *)_molAtoms[0])->m_oaBonds[a])->m_iNumber+1, i+1);
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;*/

				buf2.sprintf("%s%d=%d", ((CAtom *)g_oaAtoms[((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_baAtomIndex[((CMolAtom *)((CMolAtom *)_molAtoms[0])->m_oaBonds[a])->m_iType]])->m_sName, ((CMolAtom *)((CMolAtom *)_molAtoms[0])->m_oaBonds[a])->m_iNumber+1, i+1);
				buf.strcat(buf2);

				if(i < 3) {
/*#ifdef TARGET_LINUX
					length = snprintf(buf2, BUF_SIZE, ", ");
#else
					length = sprintf(buf2, ", ");
#endif
					strncat(buf, buf2, remaining - 1);
					remaining -= length;*/

					buf2.sprintf(", ");
					buf.strcat(buf2);
				}
			}
			mprintf("    Please enter the atoms to swap: (%s)\n", (const char*)buf);
			int a1 = AskRangeInteger_ND("      1st atom: ", 1, 4);
			int a2 = AskRangeInteger_ND("      2nd atom: ", 1, 4);
			
			if(((CMolAtom *)((CMolAtom *)_molAtoms[i])->m_oaBonds[((CxIntArray *)_orderArrays[i])->GetPosition(a1-1)])->m_fAtomCode != ((CMolAtom *)((CMolAtom *)_molAtoms[i])->m_oaBonds[((CxIntArray *)_orderArrays[i])->GetPosition(a2-1)])->m_fAtomCode)
				eprintf("Warning: These atoms are not equivalent!");
			
			_swapOffsets[0].SetMaxSize(m_iShowMolCount);
			_swapOffsets[1].SetMaxSize(m_iShowMolCount);
			for(i = 0; i < m_iShowMolCount; i++) {
				_swapOffsets[0].Add(((CMolAtom *)((CMolAtom *)_molAtoms[i])->m_oaBonds[((CxIntArray *)_orderArrays[i])->GetPosition(a1-1)])->m_iOffset);
				_swapOffsets[1].Add(((CMolAtom *)((CMolAtom *)_molAtoms[i])->m_oaBonds[((CxIntArray *)_orderArrays[i])->GetPosition(a2-1)])->m_iOffset);
			}
			
			_fixToR = AskYesNo("    Adjust chirality to R (y) or S (n)? [yes] ", true);
		}
	} else {
		_fixTraj = false;
	}
}

CChiralObservation::~CChiralObservation() {
	int i;
	for(i = 0; i < m_iShowMolCount; i++) {
		delete (CxIntArray *)_orderArrays[i];
		delete (CxFloatArray *)_timedev[i];
	}
}

void CChiralObservation::initialize() {
	int i;
	
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	mprintf("    Absolute configuration time development: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * n * sizeof(float)));
	for(i = 0; i < m_iShowMolCount; i++) {
		CxFloatArray *a;
		try { a = new CxFloatArray(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetMaxSize(n);
		a->SetGrow((int)(0.1 * n));
		_timedev.Add(a);
	}
}

void CChiralObservation::process(CTimeStep *ts) {
	int i;
	
	for(i = 0; i < m_iShowMolCount; i++) {
		CMolAtom *ma = (CMolAtom *)_molAtoms[i];
		CxIntArray *a = (CxIntArray *)_orderArrays[i];
		CxVector3 vec1, vec2, vec3;
		vec1 = FoldVector(ts->m_vaCoords[((CMolAtom *)ma->m_oaBonds[a->GetPosition(0)])->m_iOffset] - ts->m_vaCoords[((CMolAtom *)ma->m_oaBonds[a->GetPosition(1)])->m_iOffset]);
		vec2 = FoldVector(ts->m_vaCoords[ma->m_iOffset] - ts->m_vaCoords[((CMolAtom *)ma->m_oaBonds[a->GetPosition(2)])->m_iOffset]);
		vec3 = FoldVector(ts->m_vaCoords[((CMolAtom *)ma->m_oaBonds[a->GetPosition(2)])->m_iOffset] - ts->m_vaCoords[((CMolAtom *)ma->m_oaBonds[a->GetPosition(1)])->m_iOffset]);
		float dih = Dihedral(vec1, vec2, vec3, false);
		((CxFloatArray *)_timedev[i])->Add(dih);
	}
}

void CChiralObservation::finalize() {
	int i, j;	
//	char name[BUF_SIZE];
	CxString name;

/*#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "chiral_%s_%s%d.dat", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, ((CAtom *)g_oaAtoms[((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_baAtomIndex[((CMolAtom *)_molAtoms[0])->m_iType]])->m_sName, ((CMolAtom *)_molAtoms[0])->m_iNumber+1);
#else
	sprintf(name, "chiral_%s_%s%d.dat", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, ((CAtom *)g_oaAtoms[((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_baAtomIndex[((CMolAtom *)_molAtoms[0])->m_iType]])->m_sName, ((CMolAtom *)_molAtoms[0])->m_iNumber+1);
#endif*/
	name.sprintf("chiral_%s_%s%d.dat", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, ((CAtom *)g_oaAtoms[((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_baAtomIndex[((CMolAtom *)_molAtoms[0])->m_iType]])->m_sName, ((CMolAtom *)_molAtoms[0])->m_iNumber+1);

	mprintf("    Writing chirality data to %s...\n", (const char*)name);
	FILE *chiralFile = OpenFileWrite(name, false);
	
	for(i = 0; i < ((CxFloatArray *)_timedev[0])->GetSize(); i++) {
		fprintf(chiralFile, "%d", i+1);
		for(j = 0; j < m_iShowMolCount; j++) {
			if(i > 0) {
				if((((CxFloatArray *)_timedev[j])->GetAt(i)) > 0.0f && (((CxFloatArray *)_timedev[j])->GetAt(i-1) <= 0.0f))
					mprintf(RED, "Warning: Molecule %d changed configuration from S to R in step %d.\n", j+1, i+1);
				if((((CxFloatArray *)_timedev[j])->GetAt(i)) <= 0.0f && (((CxFloatArray *)_timedev[j])->GetAt(i-1) > 0.0f))
					mprintf(RED, "Warning: Molecule %d changed configuration from R to S in step %d.\n", j+1, i+1);
			}
			fprintf(chiralFile, " %s", ((CxFloatArray *)_timedev[j])->GetAt(i) > 0.0f ? "R" : "S");
		}
		fprintf(chiralFile, "\n");
	}
	fclose(chiralFile);
}

int CChiralObservation::swapAtom(int offset) {
	if(!_fixTraj)
		return offset;
	int a = _swapOffsets[0].GetPosition(offset);
	if(a != -1) {
		bool lastStepR = ((CxFloatArray *)_timedev[a])->GetAt(((CxFloatArray *)_timedev[a])->GetSize()-1) > 0.0f;
		if(lastStepR != _fixToR)
			return _swapOffsets[1].GetAt(a);
		else
			return offset;
	} else {
		a = _swapOffsets[1].GetPosition(offset);
		if(a != -1) {
			bool lastStepR = ((CxFloatArray *)_timedev[a])->GetAt(((CxFloatArray *)_timedev[a])->GetSize()-1) > 0.0f;
			if(lastStepR != _fixToR)
				return _swapOffsets[0].GetAt(a);
			else
				return offset;
		}
	}
	return offset;
}

bool gatherChiral() {
	g_fixTraj = AskYesNo("    Swap equivalent atoms to have same chirality in all molecules (y/n)? [no] ", false);
	if(g_fixTraj)
		g_bKeepOriginalCoords = true;
	mprintf("\n");
	
	while(true) {
		mprintf(YELLOW, ">>> Chirality Observation %d >>>\n\n", g_chiralObserv.GetSize() + 1);
		
		CChiralObservation *obs;
		try { obs = new CChiralObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CChiralObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_chiralObserv.Add(obs);
		
		mprintf(YELLOW, "\n<<< End of Chirality Observation %d <<<\n\n", g_chiralObserv.GetSize());
		
		if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	mprintf("\n");
	
	return true;
}

bool initializeChiral() {
	int i;
	CxString filename, buf;
	
	if(g_fixTraj) {
//		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
//		char buf[BUF_SIZE];
//		strncpy(buf, g_sInputTraj, BUF_SIZE);
//		buf[BUF_SIZE-1] = 0;
		buf.strcpy(g_sInputTraj);
		char *p = strrchr(buf.GetWritePointer(), '/');
		if (p == NULL)
			p = buf.GetWritePointer();
		else
			p++;
		size_t s = strcspn(p, ".");
		if ((int)s > buf.GetLength() - 8)
			s = (size_t)buf.GetLength() - 8;
//		strncpy(filename, p, s);
		filename.strcpy(p);
		filename(s) = 0;
//		strcat(filename, "_out.xyz");
		filename.strcat("_out.xyz");
#else
//		sprintf(filename, "out.xyz");
		filename.sprintf("out.xyz");
#endif
		mprintf("  Saving processed trajectory as %s\n\n", (const char*)filename);
		g_fixedTrajFile = OpenFileWrite(filename, false);
	}
	
	for(i = 0; i < g_chiralObserv.GetSize(); i++) {
		mprintf("  Initializing Chirality Observation %d...\n", i+1);
		((CChiralObservation *)g_chiralObserv[i])->initialize();
	}
	
	return true;
}

void processChiral(CTimeStep *ts) {
	int i, j;
	for(i = 0; i < g_chiralObserv.GetSize(); i++) {
		((CChiralObservation *)g_chiralObserv[i])->process(ts);
	}
	
	if(g_fixTraj) {
		CxIntArray fixedOrder;
		fixedOrder.SetMaxSize(g_iGesAtomCount);
		for(i = 0; i < g_iGesAtomCount; i++) {
			int si = i;
			for(j = 0; j < g_chiralObserv.GetSize(); j++) {
				si = ((CChiralObservation *)g_chiralObserv[j])->swapAtom(i);
				if(si != i)
					break;
			}
			fixedOrder.Add(si);
		}
		
		fprintf(g_fixedTrajFile, "%d\n", fixedOrder.GetSize());
		if(ts->m_pComment != NULL)
			fprintf(g_fixedTrajFile, "%s\n", ts->m_pComment);
		else
			fprintf(g_fixedTrajFile, "\n");
		for(i = 0; i < fixedOrder.GetSize(); i++) {
			fprintf(g_fixedTrajFile, "%4s %14.10f %14.10f %14.10f\n", ((CAtom *)g_oaAtoms[g_waAtomRealElement[fixedOrder[i]]])->m_sName, ts->m_vaCoords_Original[fixedOrder[i]][0] / 100.0f, ts->m_vaCoords_Original[fixedOrder[i]][1] / 100.0f, ts->m_vaCoords_Original[fixedOrder[i]][2] / 100.0f);
		}
	}
}

void finalizeChiral() {
	int i;
	
	if(g_fixTraj) {
		fclose(g_fixedTrajFile);
	}
	
	for(i = 0; i < g_chiralObserv.GetSize(); i++) {
		mprintf(YELLOW, ">>> Chirality Observation %d >>>\n\n", i+1);
		((CChiralObservation *)g_chiralObserv[i])->finalize();
		mprintf(YELLOW, "\n<<< End of Chirality Observation %d <<<\n\n", i+1);
	}
	for(i = 0; i < g_chiralObserv.GetSize(); i++)
		delete (CChiralObservation *)g_chiralObserv[i];
}
