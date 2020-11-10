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

#include "structurefactor.h"

#include "globalvar.h"
#include "maintools.h"
#include "timestep.h"
#include "xobarray.h"

#define BUF_SIZE 4096

static CxObArray g_isotopes;
static CStructureFactor *g_structureFactor;
static CxObArray g_sFacObserv;

static void createIsotopeList() {
	// Neutron factors from http://www.ncnr.nist.gov/resources/n-lengths/
	// X-Ray factors from http://www.ruppweb.org/new_comp/scattering_factors.htm
	CIsotope *isotope;
	try { isotope = new CIsotope("H", -3.7390f, 0.493f, 0.323f, 0.140f, 0.041f, 10.511f, 26.126f, 3.142f, 57.800f, 0.003f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("1H", -3.7406f, 0.493f, 0.323f, 0.140f, 0.041f, 10.511f, 26.126f, 3.142f, 57.800f, 0.003f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("D", 6.671f, 0.493f, 0.323f, 0.140f, 0.041f, 10.511f, 26.126f, 3.142f, 57.800f, 0.003f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("2H", 6.671f, 0.493f, 0.323f, 0.140f, 0.041f, 10.511f, 26.126f, 3.142f, 57.800f, 0.003f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("3H", 4.792f, 0.493f, 0.323f, 0.140f, 0.041f, 10.511f, 26.126f, 3.142f, 57.800f, 0.003f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("B", 5.30f, 2.055f, 1.333f, 1.098f, 0.707f, 23.219f, 1.021f, 60.350f, 0.140f, -0.193f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("10B", -0.1f, 2.055f, 1.333f, 1.098f, 0.707f, 23.219f, 1.021f, 60.350f, 0.140f, -0.193f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("11B", 6.65f, 2.055f, 1.333f, 1.098f, 0.707f, 23.219f, 1.021f, 60.350f, 0.140f, -0.193f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("C", 6.6460f, 2.310f, 1.020f, 1.589f, 0.865f, 20.844f, 10.208f, 0.569f, 51.651f, 0.216f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("12C", 6.6511f, 2.310f, 1.020f, 1.589f, 0.865f, 20.844f, 10.208f, 0.569f, 51.651f, 0.216f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("13C", 6.19f, 2.310f, 1.020f, 1.589f, 0.865f, 20.844f, 10.208f, 0.569f, 51.651f, 0.216f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("N", 9.36f, 12.213f, 3.132f, 2.013f, 1.166f, 0.006f, 9.893f, 28.997f, 0.583f, -11.529f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("14N", 9.37f, 12.213f, 3.132f, 2.013f, 1.166f, 0.006f, 9.893f, 28.997f, 0.583f, -11.529f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("15N", 6.44f, 12.213f, 3.132f, 2.013f, 1.166f, 0.006f, 9.893f, 28.997f, 0.583f, -11.529f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("O", 5.803f, 3.049f, 2.287f, 1.546f, 0.867f, 13.277f, 5.701f, 0.324f, 32.909f, 0.251f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("16O", 5.803f, 3.049f, 2.287f, 1.546f, 0.867f, 13.277f, 5.701f, 0.324f, 32.909f, 0.251f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("17O", 5.78f, 3.049f, 2.287f, 1.546f, 0.867f, 13.277f, 5.701f, 0.324f, 32.909f, 0.251f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("18O", 5.84f, 3.049f, 2.287f, 1.546f, 0.867f, 13.277f, 5.701f, 0.324f, 32.909f, 0.251f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("F", 5.654f, 3.539f, 2.641f, 1.517f, 1.024f, 10.283f, 4.294f, 0.262f, 26.148f, 0.278f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("S", 2.847f, 6.905f, 5.203f, 1.438f, 1.586f, 1.468f, 22.215f, 0.254f, 56.172f, 0.867f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("32S", 2.804f, 6.905f, 5.203f, 1.438f, 1.586f, 1.468f, 22.215f, 0.254f, 56.172f, 0.867f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("33S", 4.74f, 6.905f, 5.203f, 1.438f, 1.586f, 1.468f, 22.215f, 0.254f, 56.172f, 0.867f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("34S", 3.48f, 6.905f, 5.203f, 1.438f, 1.586f, 1.468f, 22.215f, 0.254f, 56.172f, 0.867f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("36S", 3.1f, 6.905f, 5.203f, 1.438f, 1.586f, 1.468f, 22.215f, 0.254f, 56.172f, 0.867f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("Cl", 9.577f, 11.460f, 7.196f, 6.256f, 1.645f, 0.010f, 1.166f, 18.519f, 47.778f, -9.557f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("35Cl", 11.65f, 11.460f, 7.196f, 6.256f, 1.645f, 0.010f, 1.166f, 18.519f, 47.778f, -9.557f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("37Cl", 3.08f, 11.460f, 7.196f, 6.256f, 1.645f, 0.010f, 1.166f, 18.519f, 47.778f, -9.557f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("Br", 6.795f, 17.179f, 5.236f, 5.638f, 3.985f, 2.172f, 16.580f, 0.261f, 41.433f, -2.956f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("79Br", 6.80f, 17.179f, 5.236f, 5.638f, 3.985f, 2.172f, 16.580f, 0.261f, 41.433f, -2.956f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("81Br", 6.79f, 17.179f, 5.236f, 5.638f, 3.985f, 2.172f, 16.580f, 0.261f, 41.433f, -2.956f); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
}

static void deleteIsotopeList() {
	int i;
	for(i = 0; i < g_isotopes.GetSize(); i++)
		delete (CIsotope *)g_isotopes[i];
}

struct StructureFactorAtomKind: public CxObject
{
	int atomType;
	CxObArray isotopeList;
};

struct StructureFactorMolecule: public CxObject
{
	int moleculeType;
	CxObArray atomKinds;
};

CIsotope::CIsotope(const char *label, float neutronFactor, float cma1, float cma2, float cma3, float cma4, float cmb1, float cmb2, float cmb3, float cmb4, float cmc) {
	try { _label = new char[strlen(label) + 1]; } catch(...) { _label = NULL; }
	if(_label == NULL) NewException((double)sizeof(char) * (strlen(label) + 1), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	strcpy(_label, label);
	_neutronFactor = neutronFactor;
	_cma[0] = cma1;
	_cma[1] = cma2;
	_cma[2] = cma3;
	_cma[3] = cma4;
	_cmb[0] = cmb1;
	_cmb[1] = cmb2;
	_cmb[2] = cmb3;
	_cmb[3] = cmb4;
	_cmc = cmc;
}

CIsotope::~CIsotope() {
	delete[] _label;
}

float CIsotope::xrayFactor(float q) {
	float x = q * 100.0f / 4.0f / Pi;
	float factor = 0.0f;
	int i;
	for(i = 0; i < 4; i++) {
		factor += _cma[i] * expf(-_cmb[i] * x * x);
	}
	factor += _cmc;
	return factor;
}

// CSFacObservation::CSFacObservation(bool global) {
// 	int i, j, k, l;
// //	char buf[BUF_SIZE], buf2[BUF_SIZE];
// 	CxString buf, buf2, temp;
// 
// 	_global = global;
// 	if(_global) {
// 		mprintf("    All atoms will be taken.\n\n");
// 		for(i = 0; i < g_oaMolecules.GetSize(); i++) {
// 			CAtomGroup *ag;
// 			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
// 			if(ag == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			if(!ag->ParseAtoms((CMolecule *)g_oaMolecules[i], "*")) {
// 				eprintf("Weird error.\n");
// 				abort();
// 			}
// 			_agList.Add(ag);
// 		}
// 	} else {
// 		while(true) {
// /*			char buf[BUF_SIZE], buf2[BUF_SIZE];
// 			size_t remaining = BUF_SIZE;
// #ifdef TARGET_LINUX
// 			remaining -= snprintf(buf, remaining, "    Take atoms from which molecule (");
// #else
// 			remaining -= sprintf(buf, "    Take atoms from which molecule (");
// #endif*/
// 			buf.sprintf("    Take atoms from which molecule (");
// 
// 			for(i = 0; i < g_oaMolecules.GetSize(); i++) {
// 
// /*				if(remaining < 1)
// 					break;
// #ifdef TARGET_LINUX
// 				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// #else
// 				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// #endif
// 				strncat(buf, buf2, remaining - 1);
// 				remaining -= length;*/
// 
// 				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// 				buf.strcat(buf2);
// 
// 
// 				if(i < g_oaMolecules.GetSize() - 1) {
// /*#ifdef TARGET_LINUX
// 					length = snprintf(buf2, remaining, ", ");
// #else
// 					length = sprintf(buf2, ", ");
// #endif
// 					strncat(buf, buf2, remaining - 1);
// 					remaining -= length;*/
// 
// 					buf.strcat(", ");
// 				}
// 			}
// //			strncat(buf, ")? ", remaining - 1);
// 			buf.strcat(")? ");
// 			
// 			int mol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
// 			
// 			CAtomGroup *ag;
// 			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
// 			if(ag == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			
// 			while(true) {
// 				if(((CMolecule *)g_oaMolecules[mol])->m_iAtomGes == 3) {
// 					mprintf("    %s is only one atom, there is no choice.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 					ag->Reset();
// 					ag->m_pMolecule = (CMolecule *)g_oaMolecules[mol];
// 					ag->AddAtom(0, 0, false);
// 					ag->SortAtoms();
// 					ag->BuildName();
// 				} else {
// 					mprintf("    Which atom(s) to take from molecule %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ", ((CMolecule*)g_oaMolecules[mol])->m_sName);
// 					inpprintf("! Which atom(s) to take from molecule %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n", ((CMolecule*)g_oaMolecules[mol])->m_sName);
// //					char buf[BUF_SIZE];
// 					myget(&buf);
// 					if(strlen(buf) == 0) {
// 						if(!ag->ParseAtoms((CMolecule *)g_oaMolecules[mol], "#2")) {
// 							eprintf("Weird error.\n");
// 							continue;
// 						}
// 					} else if(!ag->ParseAtoms((CMolecule *)g_oaMolecules[mol], buf)) {
// 						continue;
// 					}
// 				}
// 				break;
// 			}
// 			_agList.Add(ag);
// 			
// 			if(!AskYesNo("    Add atoms from another molecule (y/n)? [no] ", false))
// 				break;
// 		}
// 		mprintf("\n");
// 	}
// 	
// 	if(_global) {
// 		try { _name = new char[strlen("global") + 1]; } catch(...) { _name = NULL; }
// 		if(_name == NULL) NewException((double)sizeof(char) * (strlen("global") + 1), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		sprintf(_name, "global");
// 	} else {
// /*		char temp[BUF_SIZE];
// 		temp[0] = 0;
// 		for(i = 0; i < _agList.GetSize(); i++) {
// 			if(i > 0)
// 				strncat(temp, "_", BUF_SIZE - strlen(temp) - 1);
// 			strncat(temp, ((CAtomGroup *)_agList[i])->m_sName, BUF_SIZE - strlen(temp) - 1);
// 		}*/
// 
// 		temp.sprintf("");
// 		for(i = 0; i < _agList.GetSize(); i++) {
// 			if(i > 0)
// 				temp.strcat("_");
// 			temp.strcat(((CAtomGroup *)_agList[i])->m_sName);
// 		}
// 		
// 		try { _name = new char[strlen(temp) + 3]; } catch(...) { _name = NULL; }
// 		if(_name == NULL) NewException((double)sizeof(char) * (strlen(temp) + 3), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		sprintf(_name, "[%s]", (const char*)temp);
// 	}
// 	
// 	CxObArray isotopeAssignList;
// 	for(i = 0; i < _agList.GetSize(); i++) {
// 		CAtomGroup *ag = (CAtomGroup *)_agList[i];
// 		for(j = 0; j < ag->m_baAtomType.GetSize(); j++) {
// 			CIsotope *isotope = 0;
// 			for(k = 0; k < g_isotopes.GetSize(); k++) {
// 				if(mystricmp(((CAtom *)g_oaAtoms[ag->m_baRealAtomType[j]])->m_sName, ((CIsotope *)g_isotopes[k])->label()) == 0) {
// 					isotope = (CIsotope *)g_isotopes[k];
// 					break;
// 				}
// 			}
// 			if(isotope == 0) {
// 				mprintf(RED, "Isotope data for \"%s\" not found.\n", ((CAtom *)g_oaAtoms[ag->m_baRealAtomType[j]])->m_sName);
// 				abort();
// 			}
// 			for(k = 0; k < ((CxIntArray *)ag->m_oaAtoms[j])->GetSize(); k++) {
// 				isotopeAssignList.Add(isotope);
// 			}
// 		}
// 	}
// 	
// 	if(!AskYesNo("    Use standard atom data (y) or specify isotopes (n)? [yes] ", true)) {
// 		while(true) {
// /*			char buf[BUF_SIZE], buf2[BUF_SIZE];
// 			size_t remaining = BUF_SIZE;
// #ifdef TARGET_LINUX
// 			remaining -= snprintf(buf, remaining, "\n    Change isotopes in which molecule (");
// #else
// 			remaining -= sprintf(buf, "\n    Change isotopes in which molecule (");
// #endif*/
// 			buf.sprintf("\n    Change isotopes in which molecule (");
// 
// 			for(i = 0; i < _agList.GetSize(); i++) {
// 
// /*				if(remaining < 1)
// 					break;
// #ifdef TARGET_LINUX
// 				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)((CAtomGroup *)_agList[i])->m_pMolecule)->m_sName, i+1);
// #else
// 				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)((CAtomGroup *)_agList[i])->m_pMolecule)->m_sName, i+1);
// #endif
// 				strncat(buf, buf2, remaining - 1);
// 				remaining -= length;*/
// 
// 				buf2.sprintf("%s=%d", ((CMolecule *)((CAtomGroup *)_agList[i])->m_pMolecule)->m_sName, i+1);
// 				buf.strcat(buf2);
// 
// 				if(i < _agList.GetSize() - 1) {
// 
// /*#ifdef TARGET_LINUX
// 					length = snprintf(buf2, remaining, ", ");
// #else
// 					length = sprintf(buf2, ", ");
// #endif
// 					strncat(buf, buf2, remaining - 1);
// 					remaining -= length;*/
// 
// 					buf2.sprintf(", ");
// 					buf.strcat(buf2);
// 				}
// 			}
// //			strncat(buf, ")? ", remaining - 1);
// 			buf.strcat(")? ");
// 			
// 			int mol = AskRangeInteger_ND(buf, 1, _agList.GetSize()) - 1;
// 			
// 			int index = 0;
// 			for(i = 0; i < mol; i++) {
// 				CAtomGroup *ag = (CAtomGroup *)_agList[i];
// 				for(j = 0; j < ag->m_baAtomType.GetSize(); j++) {
// 					index += ((CxIntArray *)ag->m_oaAtoms[j])->GetSize();
// 				}
// 			}
// 			
// 			while(true) {
// 				mprintf("\n    The following isotopes are set up:\n");
// 				CAtomGroup *ag = (CAtomGroup *)_agList[mol];
// 				int index1 = index;
// 				for(i = 0; i < ag->m_baAtomType.GetSize(); i++) {
// 					CxIntArray *a = (CxIntArray *)ag->m_oaAtoms[i];
// 					for(j = 0; j < a->GetSize(); j++) {
// 						mprintf("      %s%d: %-4s", ((CAtom *)g_oaAtoms[ag->m_baRealAtomType[i]])->m_sName, a->GetAt(j)+1, ((CIsotope *)isotopeAssignList[index1++])->label());
// 						if((index1 - index) % 5 == 0)
// 							mprintf("\n");
// 					}
// 				}
// 				mprintf("\n\n");
// //				char buf[BUF_SIZE];
// 				unsigned char ca, cb, cc;
// 				int index2;
// 				while(true) {
// 					do {
// 						mprintf("    Change isotope for which atom? [done] ");
// 						inpprintf("! Change isotope for which atom? [done] ");
// 						myget(&buf);
// 						if(strlen(buf) == 0)
// 							break;
// 					} while(!ParseAtom(buf, ((CMolecule *)ag->m_pMolecule)->m_iIndex, ca, cb, cc));
// 					index2 = index;
// 					bool contains = false;
// 					for(i = 0; i < ag->m_baAtomType.GetSize(); i++) {
// 						if(cb == ag->m_baRealAtomType[i]) {
// 							CxIntArray *a = (CxIntArray *)ag->m_oaAtoms[i];
// 							for(j = 0; j < a->GetSize(); j++) {
// 								if(cc == a->GetAt(j)) {
// 									contains = true;
// 									break;
// 								}
// 								index2++;
// 							}
// 							break;
// 						} else {
// 							index2 += ((CxIntArray *)ag->m_oaAtoms[i])->GetSize();
// 						}
// 					}
// 					if(contains)
// 						break;
// 					eprintf("There is no atom \"%s\" in this observation.\n", (const char*)buf);
// 				}
// 				if(strlen(buf) == 0)
// 					break;
// //				char buf2[BUF_SIZE];
// 				while(true) {
// 					mprintf("    Which isotope to set for %s (e.g. \"13C\")? ", (const char*)buf);
// 					inpprintf("! Which isotope to set for %s (e.g. \"13C\")? ", (const char*)buf);
// 					myget(&buf2);
// 					if(strlen(buf2) == 0) {
// 						eprintf("There is no default. Please enter a character string.\n");
// 						continue;
// 					}
// 					int isotopeIndex = -1;
// 					for(i = 0; i < g_isotopes.GetSize(); i++) {
// 						if(mystricmp(((CIsotope *)g_isotopes[i])->label(), buf2) == 0) {
// 							isotopeIndex = i;
// 							break;
// 						}
// 					}
// 					if(isotopeIndex != -1) {
// 						isotopeAssignList[index2] = (CIsotope *)g_isotopes[isotopeIndex];
// 						break;
// 					}
// 					eprintf("Isotope data for \"%s\" not found.\n", (const char*)buf2);
// 				}
// 			}
// 			
// 			if(!AskYesNo("\n    Change isotopes in another molecule (y/n)? [no] ", false))
// 				break;
// 		}
// 		mprintf("\n");
// 	}
// 	
// 	int index = 0;
// 	for(i = 0; i < _agList.GetSize(); i++) {
// 		CAtomGroup *ag = (CAtomGroup *)_agList[i];
// 		CMolecule *mol = ag->m_pMolecule;
// 		for(j = 0; j < ag->m_baAtomType.GetSize(); j++) {
// 			CxIntArray *a = (CxIntArray *)ag->m_oaAtoms[j];
// 			for(k = 0; k < a->GetSize(); k++) {
// 				CIsotope *isotope = (CIsotope *)isotopeAssignList[index++];
// 				int isotopeIndex = -1;
// 				for(l = 0; l < _isotopeTypeList.GetSize(); l++)
// 					if((CIsotope *)_isotopeTypeList[l] == isotope)
// 						isotopeIndex = l;
// 				if(isotopeIndex == -1) {
// 					isotopeIndex = _isotopeTypeList.GetSize();
// 					_isotopeTypeList.Add(isotope);
// 					_isotopeTypeCount.Add(0);
// 				}
// 				for(l = 0; l < mol->m_laSingleMolIndex.GetSize(); l++) {
// 					CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[mol->m_laSingleMolIndex[l]];
// 					_atomIndex.Add(((CxIntArray *)sm->m_oaAtomOffset[ag->m_baAtomType[j]])->GetAt(a->GetAt(k)));
// 					_isotopeList.Add(isotopeIndex);
// 					_isotopeTypeCount[isotopeIndex]++;
// 				}
// 			}
// 		}
// 	}
// 	
// 	_rdfMax = AskFloat("    Enter maximum RDF distance to observe (pm): [%d] ", (float)HalfBox(), HalfBox());
// 	_rdfRes = AskUnsignedInteger("    Enter RDF binning resolution: [%d] ", 2 * HalfBox(), 2 * HalfBox());
// 	_sfacMax = AskFloat("    Enter maximum wave vector modulus (pm^-1): [0.2] ", 0.2f);
// 	_sfacRes = AskUnsignedInteger("    Enter Structure Factor resolution: [1000] ", 1000);
// 	
// 	_normFFac = AskYesNo("\n    Normalize by scattering factors (y/n)? [yes] ", true);
// 	if (_normFFac) {
// 		mprintf("    The following normalization factors are available:\n");
// 		mprintf("    (1) [Sum_i x_i*f_i(q)]^2\n");
// 		mprintf("    (2) Sum_i x_i*[f_i(q)]^2\n");
// 		_normFFacFormula = AskRangeInteger("    Which factor to use? [1] ", 1, 2, 1);
// 	} else {
// 		_normFFacFormula = 0;
// 	}
// 	mprintf("\n");
// 	
// 	if (g_bAdvanced2) {
// 		_sharpening = AskYesNo("    Apply sharpening factor (y/n)? [no] ", false);
// 		if (_sharpening) {
// //			char buf[BUF_SIZE];
// 			while (true) {
// 				mprintf("    Which isotope to use as sharpening atom (e.g. \"N\" or \"14N\")? ", (const char*)buf);
// 				inpprintf("! Which isotope to use as sharpening atom (e.g. \"N\" or \"14N\")? ", (const char*)buf);
// 				myget(&buf);
// 				if(strlen(buf) == 0) {
// 					eprintf("There is no default. Please enter a character string.\n");
// 					continue;
// 				}
// 				int isotopeIndex = -1;
// 				for (i = 0; i < g_isotopes.GetSize(); i++) {
// 					if(mystricmp(((CIsotope *)g_isotopes[i])->label(), buf) == 0) {
// 						isotopeIndex = i;
// 						break;
// 					}
// 				}
// 				if(isotopeIndex != -1) {
// 					_sharpeningIsotope = (CIsotope *)g_isotopes[isotopeIndex];
// 					break;
// 				}
// 				eprintf("Isotope data for \"%s\" not found.\n", (const char*)buf);
// 			}
// 		} else {
// 			_sharpeningIsotope = NULL;
// 		}
// 		mprintf("\n");
// 	} else {
// 		_sharpening = false;
// 		_sharpeningIsotope = NULL;
// 	}
// }
// 
// CSFacObservation::~CSFacObservation() {
// 	int i;
// 	for(i = 0; i < _agList.GetSize(); i++)
// 		delete (CAtomGroup *)_agList[i];
// 	for(i = 0; i < _rdfList.GetSize(); i++)
// 		delete (CDF *)_rdfList[i];
// 	delete[] _name;
// }
// 
// void CSFacObservation::initialize() {
// 	int i;
// 	for(i = 0; i < _isotopeTypeList.GetSize() * (_isotopeTypeList.GetSize() + 1) / 2; i++) {
// 		CDF *df;
// 		try { df = new CDF(); } catch(...) { df = NULL; }
// 		if(df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		df->m_fMinVal = 0.0f;
// 		df->m_fMaxVal = _rdfMax;
// 		df->m_iResolution = _rdfRes;
// 		df->Create();
// 		_rdfList.Add(df);
// 	}
// }
// 
// void CSFacObservation::process(CTimeStep *ts) {
// 	int i, j;
// 	for(i = 0; i < _atomIndex.GetSize(); i++) {
// 		for(j = i+1; j < _atomIndex.GetSize(); j++) {
// 			float dist = FoldedLength(ts->m_vaCoords[_atomIndex[i]] - ts->m_vaCoords[_atomIndex[j]]);
// 			int a = _isotopeList[i];
// 			int b = _isotopeList[j];
// 			if(a < b)
// 				((CDF *)_rdfList[(_isotopeTypeList.GetSize()-1)*a - a*(a-1)/2 + b])->AddToBin(dist);
// 			else
// 				((CDF *)_rdfList[(_isotopeTypeList.GetSize()-1)*b - b*(b-1)/2 + a])->AddToBin(dist);
// 		}
// 	}
// }
// 
// void CSFacObservation::finalize() {
// 	int i, j, k, l;
// 	CDF *rdfTotal;
// 	try { rdfTotal = new CDF(); } catch(...) { rdfTotal = NULL; }
// 	if(rdfTotal == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	rdfTotal->m_fMinVal = 0.0f;
// 	rdfTotal->m_fMaxVal = _rdfMax;
// 	rdfTotal->m_iResolution = _rdfRes;
// 	rdfTotal->Create();
// 	
// 	CDF *sfacTotal;
// 	try { sfacTotal = new CDF(); } catch(...) { sfacTotal = NULL; }
// 	if(sfacTotal == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	sfacTotal->m_fMinVal = 2.0f * Pi / _rdfMax;
// 	sfacTotal->m_fMaxVal = _sfacMax;
// 	sfacTotal->m_iResolution = _sfacRes;
// 	sfacTotal->Create();
// 	
// 	CDF *neutronSFac;
// 	try { neutronSFac = new CDF(); } catch(...) { neutronSFac = NULL; }
// 	if(neutronSFac == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	neutronSFac->m_fMinVal = 2.0f * Pi / _rdfMax;
// 	neutronSFac->m_fMaxVal = _sfacMax;
// 	neutronSFac->m_iResolution = _sfacRes;
// 	neutronSFac->Create();
// 	
// 	CDF *xraySFac;
// 	try { xraySFac = new CDF(); } catch(...) { xraySFac = NULL; }
// 	if(xraySFac == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	xraySFac->m_fMinVal = 2.0f * Pi / _rdfMax;
// 	xraySFac->m_fMaxVal = _sfacMax;
// 	xraySFac->m_iResolution = _sfacRes;
// 	xraySFac->Create();
// 	
// 	for(i = 0; i < _isotopeTypeList.GetSize(); i++) {
// 		for(j = i; j < _isotopeTypeList.GetSize(); j++) {
// 			mprintf("    Processing partial structure factor %s-%s\n", ((CIsotope *)_isotopeTypeList[i])->label(), ((CIsotope *)_isotopeTypeList[j])->label());
// 			CDF *df = (CDF *)_rdfList[(_isotopeTypeList.GetSize()-1)*i - i*(i-1)/2 + j];
// 			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
// 			if(i == j)
// 				df->MultiplyBin(2.0);
// 			df->CorrectRadialDist();
// 			df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0f/3.0f * Pi) / _isotopeTypeCount[i] / _isotopeTypeCount[j]);
// 			if (g_bDoubleBox)
// 				df->MultiplyBin(g_iDoubleBoxFactor);
// 			
// 			char name[BUF_SIZE];
// #ifdef TARGET_LINUX
// 			snprintf(name, BUF_SIZE, "sfac_%s_rdf_%s_%s.csv", _name, ((CIsotope *)_isotopeTypeList[i])->label(), ((CIsotope *)_isotopeTypeList[j])->label());
// #else
// 			sprintf(name, "sfac_%s_rdf_%s_%s.csv", _name, ((CIsotope *)_isotopeTypeList[i])->label(), ((CIsotope *)_isotopeTypeList[j])->label());
// #endif
// 			mprintf("    Writing RDF to %s...\n", name);
// 			df->Write("", name, "", false);
// 			
// 			for(k = 0; k < _rdfRes; k++)
// // 				rdfTotal->m_pBin[k] += df->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * _isotopeTypeCount[i] * _isotopeTypeCount[j] / _atomIndex.GetSize() / _atomIndex.GetSize();
// 				rdfTotal->m_pBin[k] += df->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * _isotopeTypeCount[i] * _isotopeTypeCount[j] / g_iGesAtomCount / g_iGesAtomCount;
// 			
// 			CDF *sfac;
// 			try { sfac = new CDF(); } catch(...) { sfac = NULL; }
// 			if(sfac == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			sfac->m_fMinVal = 2.0f * Pi / _rdfMax;
// 			sfac->m_fMaxVal = _sfacMax;
// 			sfac->m_iResolution = _sfacRes;
// 			sfac->Create();
// 			
// 			for(k = 0; k < _sfacRes; k++) {
// 				float f = 0;
// 				for(l = 0; l < _rdfRes; l++) {
// 					f += ((0.5f + l) / _rdfRes * _rdfMax) * (df->m_pBin[l] - 1.0f) / (2.0f * Pi / _rdfMax + (0.5f + k) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax)) * sinf(((0.5f + l) / _rdfRes * _rdfMax) * (2.0f * Pi / _rdfMax + (0.5f + k) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax)));
// 				}
// 				sfac->m_pBin[k] = f * 4.0f * Pi * g_iGesAtomCount / g_fBoxX / g_fBoxY / g_fBoxZ * _rdfMax / _rdfRes;
// 			}
// 			
// #ifdef TARGET_LINUX
// 			snprintf(name, BUF_SIZE, "sfac_%s_sfac_%s_%s.csv", _name, ((CIsotope *)_isotopeTypeList[i])->label(), ((CIsotope *)_isotopeTypeList[j])->label());
// #else
// 			sprintf(name, "sfac_%s_sfac_%s_%s.csv", _name, ((CIsotope *)_isotopeTypeList[i])->label(), ((CIsotope *)_isotopeTypeList[j])->label());
// #endif
// 			mprintf("    Writing Structure Factor to %s...\n", name);
// 			sfac->Write("", name, "", false);
// 			
// 			for(k = 0; k < _sfacRes; k++) {
// // 				sfacTotal->m_pBin[k] += sfac->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * _isotopeTypeCount[i] * _isotopeTypeCount[j] / _atomIndex.GetSize() / _atomIndex.GetSize();
// 				sfacTotal->m_pBin[k] += sfac->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * _isotopeTypeCount[i] * _isotopeTypeCount[j] / g_iGesAtomCount / g_iGesAtomCount;
// // 				neutronSFac->m_pBin[k] += sfac->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * _isotopeTypeCount[i] * _isotopeTypeCount[j] / _atomIndex.GetSize() / _atomIndex.GetSize() * ((CIsotope *)_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)_isotopeTypeList[j])->neutronFactor();
// 				neutronSFac->m_pBin[k] += sfac->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * _isotopeTypeCount[i] * _isotopeTypeCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)_isotopeTypeList[j])->neutronFactor();
// // 				xraySFac->m_pBin[k] += sfac->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * _isotopeTypeCount[i] * _isotopeTypeCount[j] / _atomIndex.GetSize() / _atomIndex.GetSize() * ((CIsotope *)_isotopeTypeList[i])->xrayFactor(2.0f * Pi / _rdfMax + (0.5f + k) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax)) * ((CIsotope *)_isotopeTypeList[j])->xrayFactor(2.0f * Pi / _rdfMax + (0.5f + k) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax));
// 				xraySFac->m_pBin[k] += sfac->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * _isotopeTypeCount[i] * _isotopeTypeCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)_isotopeTypeList[i])->xrayFactor(2.0f * Pi / _rdfMax + (0.5f + k) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax)) * ((CIsotope *)_isotopeTypeList[j])->xrayFactor(2.0f * Pi / _rdfMax + (0.5f + k) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax));
// 			}
// 			
// 			delete sfac;
// 		}
// 	}
// 	
// 	if(_normFFac) {
// 		float neutronFactor = 0.0f;
// 		for(i = 0; i < _isotopeTypeList.GetSize(); i++) {
// // 			neutronFactor += (float)_isotopeTypeCount[i] / _atomIndex.GetSize() * ((CIsotope *)_isotopeTypeList[i])->neutronFactor();
// 			if (_normFFacFormula == 1) {
// 				neutronFactor += (float)_isotopeTypeCount[i] / g_iGesAtomCount * ((CIsotope *)_isotopeTypeList[i])->neutronFactor();
// 			} else if (_normFFacFormula == 2) {
// 				neutronFactor += (float)_isotopeTypeCount[i] / g_iGesAtomCount * ((CIsotope *)_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)_isotopeTypeList[i])->neutronFactor();
// 			} else {
// 				eprintf("Unknown normalization factor.\n");
// 				abort();
// 			}
// 		}
// 		if (_normFFacFormula == 1) {
// 			neutronSFac->MultiplyBin(1.0f / neutronFactor / neutronFactor);
// 		} else if (_normFFacFormula == 2) {
// 			neutronSFac->MultiplyBin(1.0f / neutronFactor);
// 		} else {
// 			eprintf("Unknown normalization factor.\n");
// 			abort();
// 		}
// 		for(i = 0; i < _sfacRes; i++) {
// 			float xrayFactor = 0.0f;
// 			for(j = 0; j < _isotopeTypeList.GetSize(); j++) {
// // 				xrayFactor += (float)_isotopeTypeCount[j] / _atomIndex.GetSize() * ((CIsotope *)_isotopeTypeList[j])->xrayFactor(2.0f * Pi / _rdfMax + (0.5f + i) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax));
// 				if (_normFFacFormula == 1) {
// 					xrayFactor += (float)_isotopeTypeCount[j] / g_iGesAtomCount * ((CIsotope *)_isotopeTypeList[j])->xrayFactor(2.0f * Pi / _rdfMax + (0.5f + i) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax));
// 				} else if (_normFFacFormula == 2) {
// 					xrayFactor += (float)_isotopeTypeCount[j] / g_iGesAtomCount * ((CIsotope *)_isotopeTypeList[j])->xrayFactor(2.0f * Pi / _rdfMax + (0.5f + i) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax)) * ((CIsotope *)_isotopeTypeList[j])->xrayFactor(2.0f * Pi / _rdfMax + (0.5f + i) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax));
// 				} else {
// 					eprintf("Unknown normalization factor.\n");
// 					abort();
// 				}
// 			}
// 			if (_normFFacFormula == 1) {
// 				xraySFac->m_pBin[i] /= xrayFactor * xrayFactor;
// 			} else if (_normFFacFormula == 2) {
// 				xraySFac->m_pBin[i] /= xrayFactor;
// 			} else {
// 				eprintf("Unknown normalization factor.\n");
// 				abort();
// 			}
// 		}
// 	}
// 	
// 	char name[BUF_SIZE];
// #ifdef TARGET_LINUX
// 	snprintf(name, BUF_SIZE, "sfac_%s_rdf_total.csv", _name);
// #else
// 	sprintf(name, "sfac_%s_rdf_total.csv", _name);
// #endif
// 	mprintf("    Writing total RDF to %s...\n", name);
// 	rdfTotal->Write("", name, "", false);
// 	
// #ifdef TARGET_LINUX
// 	snprintf(name, BUF_SIZE, "sfac_%s_sfac_total.csv", _name);
// #else
// 	sprintf(name, "sfac_%s_sfac_total.csv", _name);
// #endif
// 	mprintf("    Writing total Structure Factor to %s...\n", name);
// 	sfacTotal->Write("", name, "", false);
// 	
// #ifdef TARGET_LINUX
// 	snprintf(name, BUF_SIZE, "sfac_%s_neutron.csv", _name);
// #else
// 	sprintf(name, "sfac_%s_neutron.csv", _name);
// #endif
// 	mprintf("    Writing Neutron Scattering Function to %s...\n", name);
// 	neutronSFac->Write("", name, "", false);
// 	
// #ifdef TARGET_LINUX
// 	snprintf(name, BUF_SIZE, "sfac_%s_xray.csv", _name);
// #else
// 	sprintf(name, "sfac_%s_xray.csv", _name);
// #endif
// 	mprintf("    Writing X-Ray Scattering Function to %s...\n", name);
// 	xraySFac->Write("", name, "", false);
// 	
// 	if (_sharpening) {
// 		for (i = 0; i < _sfacRes; i++) {
// 			float q = 2.0f * Pi / _rdfMax + (0.5f + i) / _sfacRes * (_sfacMax - 2.0f * Pi / _rdfMax);
// 			float neutronFactor = q * expf(-0.01f * q * q);
// 			float xrayFactor = q * expf(-0.01f * q * q) * _sharpeningIsotope->xrayFactor(0.0f) * _sharpeningIsotope->xrayFactor(0.0f) / _sharpeningIsotope->xrayFactor(q) / _sharpeningIsotope->xrayFactor(q);
// 			neutronSFac->m_pBin[i] *= neutronFactor;
// 			xraySFac->m_pBin[i] *= xrayFactor;
// 		}
// 		
// #ifdef TARGET_LINUX
// 		snprintf(name, BUF_SIZE, "sfac_%s_neutron_sharpened.csv", _name);
// #else
// 		sprintf(name, "sfac_%s_neutron.csv", _name);
// #endif
// 		mprintf("    Writing sharpened Neutron Scattering Function to %s...\n", name);
// 		neutronSFac->Write("", name, "", false);
// 		
// #ifdef TARGET_LINUX
// 		snprintf(name, BUF_SIZE, "sfac_%s_xray_sharpened.csv", _name);
// #else
// 		sprintf(name, "sfac_%s_xray.csv", _name);
// #endif
// 		mprintf("    Writing sharpened X-Ray Scattering Function to %s...\n", name);
// 		xraySFac->Write("", name, "", false);
// 	}
// 	
// 	delete rdfTotal;
// 	delete sfacTotal;
// 	delete neutronSFac;
// 	delete xraySFac;
// }

CStructureFactorGroup::CStructureFactorGroup(bool global, CxObArray &isotopeAssignList) {
	int i, j, k, l;
	if (global) {
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			CAtomGroup *ag;
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			if(!ag->ParseAtoms((CMolecule *)g_oaMolecules[i], "*")) {
				eprintf("Weird error.\n");
				abort();
			}
			m_atomGroupList.Add(ag);
		}
		m_sepInterIntra = false;
	} else {
		CxString buf, buf2;
		while (true) {
			buf.sprintf("    Take atoms from which molecule (");
			for (i = 0; i < g_oaMolecules.GetSize(); i++) {
				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
				buf.strcat(buf2);
				if (i < g_oaMolecules.GetSize() - 1) {
					buf.strcat(", ");
				}
			}
			buf.strcat(")? ");
			
			int mol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
			
			CAtomGroup *ag;
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			
			while (true) {
				if (((CMolecule *)g_oaMolecules[mol])->m_iAtomGes == 3) {
					mprintf("    %s is only one atom, there is no choice.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
					ag->Reset();
					ag->m_pMolecule = (CMolecule *)g_oaMolecules[mol];
					ag->AddAtom(0, 0, false);
					ag->SortAtoms();
					ag->BuildName();
				} else {
					AskString("    Which atom(s) to take from molecule %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [*] ", &buf, "*", ((CMolecule*)g_oaMolecules[mol])->m_sName);
					if(!ag->ParseAtoms((CMolecule *)g_oaMolecules[mol], buf)) {
						continue;
					}
					bool containsVirt = false;
					for (i = 0; i < ag->m_baAtomType.GetSize(); i++) {
						if (ag->m_baAtomType[i] == g_iVirtAtomType) {
							containsVirt = true;
							break;
						}
					}
					if (containsVirt) {
						eprintf("The selection must not contain virtual atoms.\n");
						continue;
					}
				}
				break;
			}
			m_atomGroupList.Add(ag);
			
			if (!AskYesNo("\n    Add atoms from another molecule (y/n)? [no] ", false))
				break;
			mprintf("\n");
		}
		mprintf("\n");
		
		m_sepInterIntra = AskYesNo("    Separate intramolecular and intermolecular contributions (y/n)? [no] ", false);
		mprintf("\n");
	}
	
	if (global) {
		m_name.sprintf("total");
	} else {
		CxString temp;
		for (i = 0; i < m_atomGroupList.GetSize(); i++) {
			if (i > 0)
				temp.strcat("__");
			temp.strcat(((CAtomGroup *)m_atomGroupList[i])->m_pMolecule->m_sName);
			temp.strcat("_");
			temp.strcat(((CAtomGroup *)m_atomGroupList[i])->m_sName);
		}
		
		m_name.sprintf("[%s]", (const char *)temp);
	}
	
	for (i = 0; i < m_atomGroupList.GetSize(); i++) {
		CAtomGroup *ag = (CAtomGroup *)m_atomGroupList[i];
		CMolecule *mol = ag->m_pMolecule;
		for (j = 0; j < ag->m_baAtomType.GetSize(); j++) {
			CxIntArray *a = (CxIntArray *)ag->m_oaAtoms[j];
			for (k = 0; k < a->GetSize(); k++) {
				CIsotope *isotope = (CIsotope *)((StructureFactorAtomKind *)((StructureFactorMolecule *)isotopeAssignList.GetAt(mol->m_iIndex))->atomKinds[ag->m_baAtomType[j]])->isotopeList[k];
				int isotopeIndex = -1;
				for (l = 0; l < m_isotopeTypeList.GetSize(); l++) {
					if((CIsotope *)m_isotopeTypeList[l] == isotope) {
						isotopeIndex = l;
						break;
					}
				}
				if (isotopeIndex == -1) {
					isotopeIndex = m_isotopeTypeList.GetSize();
					m_isotopeTypeList.Add(isotope);
					m_isotopeTypeCount.Add(0);
				}
				for (l = 0; l < mol->m_laSingleMolIndex.GetSize(); l++) {
					CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[mol->m_laSingleMolIndex[l]];
					m_atomIndexList.Add(((CxIntArray *)sm->m_oaAtomOffset[ag->m_baAtomType[j]])->GetAt(a->GetAt(k)));
					m_singleMolList.Add(mol->m_laSingleMolIndex[l]);
					m_isotopeList.Add(isotopeIndex);
					m_isotopeTypeCount[isotopeIndex]++;
				}
			}
		}
	}
	
	m_isotopeTypeTotalCount.SetSize(m_isotopeTypeCount.GetSize());
	for (i = 0; i < m_isotopeTypeTotalCount.GetSize(); i++)
		m_isotopeTypeTotalCount[i] = 0;
	
	for (i = 0; i < isotopeAssignList.GetSize(); i++) {
		StructureFactorMolecule *smol = (StructureFactorMolecule *)isotopeAssignList[i];
		for (j = 0; j < smol->atomKinds.GetSize(); j++) {
			StructureFactorAtomKind *satom = (StructureFactorAtomKind *)smol->atomKinds[j];
			for (k = 0; k < satom->isotopeList.GetSize(); k++) {
				CIsotope *isotope = (CIsotope *)satom->isotopeList[k];
				int isotopeIndex = -1;
				for (l = 0; l < m_isotopeTypeList.GetSize(); l++) {
					if ((CIsotope *)m_isotopeTypeList[l] == isotope) {
						isotopeIndex = l;
						break;
					}
				}
				if (isotopeIndex != -1) {
					m_isotopeTypeTotalCount[isotopeIndex] += ((CMolecule *)g_oaMolecules[i])->m_laSingleMolIndex.GetSize();
				}
			}
		}
	}
	
	m_global = global;
}

CStructureFactorGroup::~CStructureFactorGroup() {
	int i;
	for (i = 0; i < m_atomGroupList.GetSize(); i++)
		delete (CAtomGroup *)m_atomGroupList[i];
	for (i = 0; i < m_rdfList.GetSize(); i++)
		delete (CDF *)m_rdfList[i];
	for (i = 0; i < m_rdfIntraList.GetSize(); i++)
		delete (CDF *)m_rdfIntraList[i];
	for (i = 0; i < m_rdfInterList.GetSize(); i++)
		delete (CDF *)m_rdfInterList[i];
}

void CStructureFactorGroup::initialize(float rdfMax, int rdfRes) {
	int i;
	mprintf("    Creating %d RDFs\n", m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
	for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) {
		CDF *df;
		try { df = new CDF(); } catch(...) { df = NULL; }
		if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		df->m_fMinVal = 0.0f;
		df->m_fMaxVal = rdfMax;
		df->m_iResolution = rdfRes;
		df->SetLabelX("Distance / pm");
		df->SetLabelY("g(r)");
		df->Create();
		m_rdfList.Add(df);
	}
	
	if (m_sepInterIntra) {
		mprintf("    Creating %d intramolecular RDFs\n", m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
		for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) {
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0f;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_rdfIntraList.Add(df);
		}
		mprintf("    Creating %d intermolecular RDFs\n", m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
		for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) {
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0f;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_rdfInterList.Add(df);
		}
	}
}

void CStructureFactorGroup::process(CTimeStep *ts) {
	int i, j;
	for (i = 0; i < m_atomIndexList.GetSize(); i++) {
		for (j = i+1; j < m_atomIndexList.GetSize(); j++) {
			float dist = FoldedLength(ts->m_vaCoords[m_atomIndexList[i]] - ts->m_vaCoords[m_atomIndexList[j]]);
			int a = m_isotopeList[i];
			int b = m_isotopeList[j];
			int index;
			if (a < b)
				index = (m_isotopeTypeList.GetSize()-1)*a - a*(a-1)/2 + b;
			else
				index = (m_isotopeTypeList.GetSize()-1)*b - b*(b-1)/2 + a;
			((CDF *)m_rdfList[index])->AddToBin(dist);
			if (m_sepInterIntra) {
				if (m_singleMolList[i] == m_singleMolList[j])
					((CDF *)m_rdfIntraList[index])->AddToBin(dist);
				else
					((CDF *)m_rdfInterList[index])->AddToBin(dist);
			}
		}
	}
}

void CStructureFactorGroup::finalize(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_isotopeTypeList.GetSize(); i++) {
		for (j = i; j < m_isotopeTypeList.GetSize(); j++) {
			mprintf("    Processing partial structure factor %s-%s\n", ((CIsotope *)m_isotopeTypeList[i])->label(), ((CIsotope *)m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfList[(m_isotopeTypeList.GetSize()-1)*i - i*(i-1)/2 + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if(i == j)
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0f / (4.0f/3.0f * Pi) / m_isotopeTypeTotalCount[i] / m_isotopeTypeTotalCount[j]);
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0f/3.0f * Pi) / m_isotopeTypeTotalCount[i] / m_isotopeTypeTotalCount[j]);
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
// 				df->Integrate(true, (4.0f/3.0f * Pi) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / (g_fBoxX * g_fBoxY * g_fBoxZ) / m_isotopeTypeCount[i]);
				filename.sprintf("sfac_%s%s_rdf_%s_%s.csv", m_global ? "" : "self_", (const char *)m_name, (const char *)(((CIsotope *)m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_isotopeTypeList[j])->label()));
				mprintf("    Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			float shift = (float)m_isotopeTypeCount[i] * m_isotopeTypeCount[j] / (m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j]);
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					float r = (0.5f + l) / df->m_fFac;
					f += r * (df->m_pBin[l] - shift) / q * sinf(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount;
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0f);
				else
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_%s%s_sfac_%s_%s.csv", m_global ? "" : "self_", (const char *)m_name, (const char *)((CIsotope *)m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_isotopeTypeList[j])->label());
				mprintf("    Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}

void CStructureFactorGroup::finalizeIntra(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_isotopeTypeList.GetSize(); i++) {
		for (j = i; j < m_isotopeTypeList.GetSize(); j++) {
			mprintf("    Processing intramolecular partial structure factor %s-%s\n", ((CIsotope *)m_isotopeTypeList[i])->label(), ((CIsotope *)m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfIntraList[(m_isotopeTypeList.GetSize()-1)*i - i*(i-1)/2 + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if(i == j)
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0f / (4.0f/3.0f * Pi) / m_isotopeTypeTotalCount[i] / m_isotopeTypeTotalCount[j]);
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0f/3.0f * Pi) / m_isotopeTypeTotalCount[i] / m_isotopeTypeTotalCount[j]);
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_self_%s_rdf_%s_%s_intra.csv", (const char *)m_name, (const char *)(((CIsotope *)m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_isotopeTypeList[j])->label()));
				mprintf("    Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			float shift = 0.0f;
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					float r = (0.5f + l) / df->m_fFac;
					f += r * (df->m_pBin[l] - shift) / q * sinf(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount;
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0f);
				else
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_self_%s_sfac_%s_%s_intra.csv", (const char *)m_name, (const char *)((CIsotope *)m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_isotopeTypeList[j])->label());
				mprintf("    Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}

void CStructureFactorGroup::finalizeInter(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_isotopeTypeList.GetSize(); i++) {
		for (j = i; j < m_isotopeTypeList.GetSize(); j++) {
			mprintf("    Processing intermolecular partial structure factor %s-%s\n", ((CIsotope *)m_isotopeTypeList[i])->label(), ((CIsotope *)m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfInterList[(m_isotopeTypeList.GetSize()-1)*i - i*(i-1)/2 + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if(i == j)
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0f / (4.0f/3.0f * Pi) / m_isotopeTypeTotalCount[i] / m_isotopeTypeTotalCount[j]);
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0f/3.0f * Pi) / m_isotopeTypeTotalCount[i] / m_isotopeTypeTotalCount[j]);
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_self_%s_rdf_%s_%s_inter.csv", (const char *)m_name, (const char *)(((CIsotope *)m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_isotopeTypeList[j])->label()));
				mprintf("    Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			float shift = (float)m_isotopeTypeCount[i] * m_isotopeTypeCount[j] / (m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j]);
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					float r = (0.5f + l) / df->m_fFac;
					f += r * (df->m_pBin[l] - shift) / q * sinf(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount;
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0f : 2.0f) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0f);
				else
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_self_%s_sfac_%s_%s_inter.csv", (const char *)m_name, (const char *)((CIsotope *)m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_isotopeTypeList[j])->label());
				mprintf("    Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}

CStructureFactorCross::CStructureFactorCross(CStructureFactorGroup *sfacGroup1, CStructureFactorGroup *sfacGroup2) {
	m_sfacGroup1 = sfacGroup1;
	m_sfacGroup2 = sfacGroup2;
	m_name.sprintf("%s_%s", (const char *)sfacGroup1->getName(), (const char *)sfacGroup2->getName());
	m_sepInterIntra = false;
}

CStructureFactorCross::~CStructureFactorCross() {
	int i;
	for (i = 0; i < m_rdfList.GetSize(); i++)
		delete (CDF *)m_rdfList[i];
	for (i = 0; i < m_rdfIntraList.GetSize(); i++)
		delete (CDF *)m_rdfIntraList[i];
	for (i = 0; i < m_rdfInterList.GetSize(); i++)
		delete (CDF *)m_rdfInterList[i];
}

void CStructureFactorCross::initialize(float rdfMax, int rdfRes) {
	int i;
	mprintf("    Creating %d RDFs\n", m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize());
	for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize(); i++) {
		CDF *df;
		try { df = new CDF(); } catch(...) { df = NULL; }
		if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		df->m_fMinVal = 0.0f;
		df->m_fMaxVal = rdfMax;
		df->m_iResolution = rdfRes;
		df->SetLabelX("Distance / pm");
		df->SetLabelY("g(r)");
		df->Create();
		m_rdfList.Add(df);
	}
	
	if (m_sepInterIntra) {
		mprintf("    Creating %d intramolecular RDFs\n", m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize());
		for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize(); i++) {
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0f;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_rdfIntraList.Add(df);
		}
		mprintf("    Creating %d intermolecular RDFs\n", m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize());
		for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize(); i++) {
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0f;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_rdfInterList.Add(df);
		}
	}
}

void CStructureFactorCross::process(CTimeStep* ts) {
	int i, j;
	for(i = 0; i < m_sfacGroup1->m_atomIndexList.GetSize(); i++) {
		for(j = 0; j < m_sfacGroup2->m_atomIndexList.GetSize(); j++) {
			float dist = FoldedLength(ts->m_vaCoords[m_sfacGroup1->m_atomIndexList[i]] - ts->m_vaCoords[m_sfacGroup2->m_atomIndexList[j]]);
			int a = m_sfacGroup1->m_isotopeList[i];
			int b = m_sfacGroup2->m_isotopeList[j];
			((CDF *)m_rdfList[a * m_sfacGroup2->m_isotopeTypeList.GetSize() + b])->AddToBin(dist);
			if (m_sepInterIntra) {
				if (m_sfacGroup1->m_singleMolList[i] == m_sfacGroup2->m_singleMolList[j])
					((CDF *)m_rdfIntraList[a * m_sfacGroup2->m_isotopeTypeList.GetSize() + b])->AddToBin(dist);
				else
					((CDF *)m_rdfInterList[a * m_sfacGroup2->m_isotopeTypeList.GetSize() + b])->AddToBin(dist);
			}
		}
	}
}

void CStructureFactorCross::finalize(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize(); i++) {
		for (j = 0; j < m_sfacGroup2->m_isotopeTypeList.GetSize(); j++) {
			mprintf("    Processing partial structure factor %s-%s\n", ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfList[i * m_sfacGroup2->m_isotopeTypeList.GetSize() + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				df->MultiplyBin(2.0f);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0f / (4.0f/3.0f * Pi) / m_sfacGroup1->m_isotopeTypeTotalCount[i] / m_sfacGroup2->m_isotopeTypeTotalCount[j]);
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0f/3.0f * Pi) / m_sfacGroup1->m_isotopeTypeTotalCount[i] / m_sfacGroup2->m_isotopeTypeTotalCount[j]);
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_cross_%s_rdf_%s_%s.csv", (const char *)m_name, (const char *)(((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label()));
				mprintf("    Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			float shift = (float)m_sfacGroup1->m_isotopeTypeCount[i] * m_sfacGroup2->m_isotopeTypeCount[j] / (m_sfacGroup1->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j]);
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				shift *= 2.0f;
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					float r = (0.5f + l) / df->m_fFac;
					f += r * (df->m_pBin[l] - shift) / q * sinf(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = (CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j] ? 1.0f : 2.0f;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup1->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount;
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup1->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup2->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0f);
				else
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_cross_%s_sfac_%s_%s.csv", (const char *)m_name, (const char *)((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
				mprintf("    Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}

void CStructureFactorCross::finalizeIntra(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize(); i++) {
		for (j = 0; j < m_sfacGroup2->m_isotopeTypeList.GetSize(); j++) {
			mprintf("    Processing intramolecular partial structure factor %s-%s\n", ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfIntraList[i * m_sfacGroup2->m_isotopeTypeList.GetSize() + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				df->MultiplyBin(2.0f);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0f / (4.0f/3.0f * Pi) / m_sfacGroup1->m_isotopeTypeTotalCount[i] / m_sfacGroup2->m_isotopeTypeTotalCount[j]);
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0f/3.0f * Pi) / m_sfacGroup1->m_isotopeTypeTotalCount[i] / m_sfacGroup2->m_isotopeTypeTotalCount[j]);
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_cross_%s_rdf_%s_%s_intra.csv", (const char *)m_name, (const char *)(((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label()));
				mprintf("    Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			float shift = 0.0f;
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				shift *= 2.0f;
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					float r = (0.5f + l) / df->m_fFac;
					f += r * (df->m_pBin[l] - shift) / q * sinf(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = (CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j] ? 1.0f : 2.0f;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup1->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount;
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup1->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup2->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0f);
				else
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_cross_%s_sfac_%s_%s_intra.csv", (const char *)m_name, (const char *)((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
				mprintf("    Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}

void CStructureFactorCross::finalizeInter(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize(); i++) {
		for (j = 0; j < m_sfacGroup2->m_isotopeTypeList.GetSize(); j++) {
			mprintf("    Processing intermolecular partial structure factor %s-%s\n", ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfInterList[i * m_sfacGroup2->m_isotopeTypeList.GetSize() + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				df->MultiplyBin(2.0f);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0f / (4.0f/3.0f * Pi) / m_sfacGroup1->m_isotopeTypeTotalCount[i] / m_sfacGroup2->m_isotopeTypeTotalCount[j]);
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0f/3.0f * Pi) / m_sfacGroup1->m_isotopeTypeTotalCount[i] / m_sfacGroup2->m_isotopeTypeTotalCount[j]);
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_cross_%s_rdf_%s_%s_inter.csv", (const char *)m_name, (const char *)(((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label()));
				mprintf("    Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			float shift = (float)m_sfacGroup1->m_isotopeTypeCount[i] * m_sfacGroup2->m_isotopeTypeCount[j] / (m_sfacGroup1->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j]);
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				shift *= 2.0f;
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					float r = (0.5f + l) / df->m_fFac;
					f += r * (df->m_pBin[l] - shift) / q * sinf(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				float q = (sfacTemp->m_fMinVal + (0.5f + k) / sfacTemp->m_fFac) / 1000.0f;
				float f = (CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j] ? 1.0f : 2.0f;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup1->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount;
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup1->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * f * m_sfacGroup2->m_isotopeTypeTotalCount[i] * m_sfacGroup2->m_isotopeTypeTotalCount[j] / g_iGesAtomCount / g_iGesAtomCount * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0f);
				else
					sfacTemp->MultiplyBin(4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_cross_%s_sfac_%s_%s_inter.csv", (const char *)m_name, (const char *)((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
				mprintf("    Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}

CStructureFactor::CStructureFactor() {
	int i, j, k;
	CxString buf, buf2, temp;
	
	CxObArray isotopeAssignList;
	for (i = 0; i < g_oaMolecules.GetSize(); i++) {
		CMolecule *mol = (CMolecule *)g_oaMolecules[i];
		StructureFactorMolecule *smol;
		try { smol = new StructureFactorMolecule; } catch(...) { smol = NULL; }
		if (smol == NULL) NewException((double)sizeof(StructureFactorMolecule), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		smol->moleculeType = i;
		for (j = 0; j < mol->m_baAtomIndex.GetSize(); j++) {
			if (mol->m_baAtomIndex[j] == g_iVirtAtomType)
				continue;
			StructureFactorAtomKind *satom;
			try { satom = new StructureFactorAtomKind; } catch(...) { satom = NULL; }
			if (satom == NULL) NewException((double)sizeof(StructureFactorAtomKind), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			satom->atomType = mol->m_baAtomIndex[j];
			CIsotope *isotope = NULL;
			for (k = 0; k < g_isotopes.GetSize(); k++) {
				if (mystricmp(((CAtom *)g_oaAtoms[mol->m_baAtomIndex[j]])->m_sName, ((CIsotope *)g_isotopes[k])->label()) == 0) {
					isotope = (CIsotope *)g_isotopes[k];
					break;
				}
			}
			if (isotope == NULL) {
				mprintf(RED, "Warning: no isotope data for \"%s\" found. Setting all scattering factors to 0.\n", ((CAtom *)g_oaAtoms[mol->m_baAtomIndex[j]])->m_sName);
				try { isotope = new CIsotope(((CAtom *)g_oaAtoms[mol->m_baAtomIndex[j]])->m_sName, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f); } catch(...) { isotope = NULL; }
				if (isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				g_isotopes.Add(isotope);
			}
			for (k = 0; k < mol->m_waAtomCount[j]; k++) {
				satom->isotopeList.Add(isotope);
			}
			smol->atomKinds.Add(satom);
		}
		isotopeAssignList.Add(smol);
	}
	
	if (g_bAdvanced2) {
		if (!AskYesNo("    Use standard atom data (y) or specify isotopes (n)? [yes] ", true)) {
			while (true) {
				buf.sprintf("\n    Change isotopes in which molecule (");
				for(i = 0; i < g_oaMolecules.GetSize(); i++) {
					buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
					buf.strcat(buf2);
					if(i < g_oaMolecules.GetSize() - 1)
						buf.strcat(", ");
				}
				buf.strcat(")? ");
				int mol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
				CMolecule *m = (CMolecule *)g_oaMolecules[mol];
				
				while (true) {
					mprintf("\n    The following isotopes are set up:\n");
					
					for (i = 0; i < m->m_baAtomIndex.GetSize(); i++) {
						if (m->m_baAtomIndex[i] == g_iVirtAtomType)
							continue;
						for (j = 0; j < m->m_waAtomCount[i]; j++) {
							mprintf("      %s%d: %-4s", ((CAtom *)g_oaAtoms[m->m_baAtomIndex[i]])->m_sName, j+1, ((CIsotope *)((StructureFactorAtomKind *)((StructureFactorMolecule *)isotopeAssignList[mol])->atomKinds[i])->isotopeList[j])->label());
						}
					}
					mprintf("\n\n");
					
					unsigned char ca, cb, cc;
					do {
						AskString("    Change isotope for which atom? [done] ", &buf, "");
						if (buf.GetLength() == 0)
							break;
					} while (!ParseAtom(buf, mol, ca, cb, cc));
					if (buf.GetLength() == 0)
						break;
					
					while (true) {
						AskString_ND("    Which isotope to set for %s (e.g. \"13C\")? ", &buf2, (const char *)buf);
						int isotopeIndex = -1;
						for (i = 0; i < g_isotopes.GetSize(); i++) {
							if(mystricmp(((CIsotope *)g_isotopes[i])->label(), buf2) == 0) {
								isotopeIndex = i;
								break;
							}
						}
						if (isotopeIndex != -1) {
							((StructureFactorAtomKind *)((StructureFactorMolecule *)isotopeAssignList[mol])->atomKinds[ca])->isotopeList[cc] = (CIsotope *)g_isotopes[isotopeIndex];
							break;
						}
						mprintf("\n    No isotope data for \"%s\" found.\n", (const char *)buf2);
						if (AskYesNo("    Enter data for \"%s\" (y/n)? [yes] ", true, (const char *)buf2)) {
							float nf = AskFloat_ND("    Enter neutron scattering factor: ");
							float cma1 = AskFloat_ND("    Enter Cromer-Mann coefficient a1: ");
							float cma2 = AskFloat_ND("    Enter Cromer-Mann coefficient a2: ");
							float cma3 = AskFloat_ND("    Enter Cromer-Mann coefficient a3: ");
							float cma4 = AskFloat_ND("    Enter Cromer-Mann coefficient a4: ");
							float cmb1 = AskFloat_ND("    Enter Cromer-Mann coefficient b1: ");
							float cmb2 = AskFloat_ND("    Enter Cromer-Mann coefficient b2: ");
							float cmb3 = AskFloat_ND("    Enter Cromer-Mann coefficient b3: ");
							float cmb4 = AskFloat_ND("    Enter Cromer-Mann coefficient b4: ");
							float cmc = AskFloat_ND("    Enter Cromer-Mann coefficient c: ");
							CIsotope *isotope;
							try { isotope = new CIsotope(buf2, nf, cma1, cma2, cma3, cma4, cmb1, cmb2, cmb3, cmb4, cmc); } catch(...) { isotope = NULL; }
							if (isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
							g_isotopes.Add(isotope);
							((StructureFactorAtomKind *)((StructureFactorMolecule *)isotopeAssignList[mol])->atomKinds[ca])->isotopeList[cc] = isotope;
							break;
						} else {
							mprintf("\n");
						}
					}
				}
				
				if (!AskYesNo("\n    Change isotopes in another molecule (y/n)? [no] ", false))
					break;
			}
			mprintf("\n");
		} else {
			mprintf("\n");
		}
		
		if (AskYesNo("    Change the predefined scattering factors (y/n)? [no] ", false)) {
			CxObArray isotopeTypeList;
			for (i = 0; i < isotopeAssignList.GetSize(); i++) {
				StructureFactorMolecule *smol = (StructureFactorMolecule *)isotopeAssignList[i];
				for (j = 0; j < smol->atomKinds.GetSize(); j++) {
					StructureFactorAtomKind *satom = (StructureFactorAtomKind *)smol->atomKinds[j];
					for (k = 0; k < satom->isotopeList.GetSize(); k++) {
						CIsotope *isotope = (CIsotope *)satom->isotopeList[k];
						bool found = false;
						int l;
						for (l = 0; l < isotopeTypeList.GetSize(); l++) {
							if ((CIsotope *)isotopeTypeList[l] == isotope) {
								found = true;
								break;
							}
						}
						if (!found) {
							isotopeTypeList.Add(isotope);
						}
					}
				}
			}
			
			while (true) {
				mprintf("\n    The following scattering factors are set up:\n\n");
				mprintf("    Name    neutron       a1       a2       a3       a4       b1       b2       b3       b4        c\n");
				for (i = 0; i < isotopeTypeList.GetSize(); i++) {
					CIsotope *isotope = (CIsotope *)isotopeTypeList[i];
					mprintf("    %-6s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", isotope->_label, isotope->_neutronFactor, isotope->_cma[0], isotope->_cma[1], isotope->_cma[2], isotope->_cma[3], isotope->_cmb[0], isotope->_cmb[1], isotope->_cmb[2], isotope->_cmb[3], isotope->_cmc);
				}
				mprintf("\n");
				CIsotope *isotope = NULL;
				while (true) {
					AskString("    Change scattering factors of which isotope? [done] ", &buf, "");
					if (buf.GetLength() == 0)
						break;
					for (i = 0; i < isotopeTypeList.GetSize(); i++) {
						if (mystricmp(((CIsotope *)isotopeTypeList[i])->label(), buf) == 0) {
							isotope = (CIsotope *)isotopeTypeList[i];
							break;
						}
					}
					if (isotope != NULL)
						break;
					eprintf("The system does not contain isotope \"%s\".\n", (const char *)buf);
				}
				if (buf.GetLength() == 0)
					break;
				isotope->_neutronFactor = AskFloat_ND("    Enter neutron scattering factor: ");
				isotope->_cma[0] = AskFloat_ND("    Enter Cromer-Mann coefficient a1: ");
				isotope->_cma[1] = AskFloat_ND("    Enter Cromer-Mann coefficient a2: ");
				isotope->_cma[2] = AskFloat_ND("    Enter Cromer-Mann coefficient a3: ");
				isotope->_cma[3] = AskFloat_ND("    Enter Cromer-Mann coefficient a4: ");
				isotope->_cmb[0] = AskFloat_ND("    Enter Cromer-Mann coefficient b1: ");
				isotope->_cmb[1] = AskFloat_ND("    Enter Cromer-Mann coefficient b2: ");
				isotope->_cmb[2] = AskFloat_ND("    Enter Cromer-Mann coefficient b3: ");
				isotope->_cmb[3] = AskFloat_ND("    Enter Cromer-Mann coefficient b4: ");
				isotope->_cmc = AskFloat_ND("    Enter Cromer-Mann coefficient c: ");
			}
			mprintf("\n");
		} else {
			mprintf("\n");
		}
	}
	
	try { m_globalSFac = new CStructureFactorGroup(true, isotopeAssignList); } catch(...) { m_globalSFac = NULL; }
	if (m_globalSFac == NULL) NewException((double)sizeof(CStructureFactorGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	if (AskYesNo("    Compute structure factor of specific atom groups (y/n)? [no] ", false)) {
		while (true) {
			mprintf(WHITE, "\n    > Group %d >\n\n", m_sFacGroups.GetSize() + 1);
			CStructureFactorGroup *sfac;
			try { sfac = new CStructureFactorGroup(false, isotopeAssignList); } catch(...) { sfac = NULL; }
			if (sfac == NULL) NewException((double)sizeof(CStructureFactorGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			m_sFacGroups.Add(sfac);
			mprintf(WHITE, "    < End of Group %d <\n\n", m_sFacGroups.GetSize());
			if (!AskYesNo("    Add another group (y/n)? [no] ", false))
				break;
		}
	}
	mprintf("\n");
	
	if (m_sFacGroups.GetSize() > 1) {
		if (AskYesNo("    Compute all cross terms of the defined atom groups (y) or select certain combinations (n)? [yes] ", true)) {
			mprintf("\n");
			int count = 0;
			for (i = 0; i < m_sFacGroups.GetSize(); i++) {
				for (j = i + 1; j < m_sFacGroups.GetSize(); j++) {
					CStructureFactorCross *sfac;
					try { sfac = new CStructureFactorCross((CStructureFactorGroup *)m_sFacGroups[i], (CStructureFactorGroup *)m_sFacGroups[j]); } catch(...) { sfac = NULL; }
					if (sfac == NULL) NewException((double)sizeof(CStructureFactorGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
					m_sFacCrosses.Add(sfac);
					mprintf("    Cross term %d: Group %d - Group %d\n", ++count, i + 1, j + 1);
				}
			}
		} else {
			mprintf("    The following groups are defined:\n\n");
			for (i = 0; i < m_sFacGroups.GetSize(); i++) {
				mprintf("    %3d - %s\n", i + 1, (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			}
			mprintf("\n    Enter each combination as a comma-separated pair of numbers (e.g. 1,1)\n\n");
			CxWordArray wa;
			while (true) {
				wa.RemoveAll_KeepSize();
				AskString("    Enter combination (return=finished): ", &buf, "");
				if (buf.GetLength() == 0)
					break;
				ParseIntList(buf, &wa);
				if (wa.GetSize() != 2) {
					eprintf("Wrong input. Please enter exactly two numbers separated by a comma.\n");
					continue;
				}
				if ((wa[0] < 1) || (wa[1] < 1)) {
					eprintf("Wrong input. Please enter positive numbers.\n");
					continue;
				}
				if ((wa[0] > m_sFacGroups.GetSize()) || (wa[1] > m_sFacGroups.GetSize())) {
					eprintf("Wrong input. There are only %d groups.\n", m_sFacGroups.GetSize());
					continue;
				}
				if (wa[0] == wa[1]) {
					eprintf("Wrong input. Please enter two different numbers.\n");
					continue;
				}
				CStructureFactorCross *sfac;
				try { sfac = new CStructureFactorCross((CStructureFactorGroup *)m_sFacGroups[wa[0] - 1], (CStructureFactorGroup *)m_sFacGroups[wa[1] - 1]); } catch(...) { sfac = NULL; }
				if (sfac == NULL) NewException((double)sizeof(CStructureFactorGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				m_sFacCrosses.Add(sfac);
			}
		}
		mprintf("\n");
	}
	
	if (m_sFacCrosses.GetSize() > 0) {
		if (AskYesNo("    Separate intramolecular and intermolecular contributions in the cross terms (y/n)? [no] ", false)) {
			for (i = 0; i < m_sFacCrosses.GetSize(); i++) {
				((CStructureFactorCross *)m_sFacCrosses[i])->setSepInterIntra(true);
			}
		}
		mprintf("\n");
	}
	
	m_rdfMax = AskFloat("    Enter maximum RDF distance to observe (pm): [%d] ", (float)HalfBox(), HalfBox());
	m_rdfRes = AskUnsignedInteger("    Enter RDF binning resolution: [%d] ", 2 * HalfBox(), 2 * HalfBox());
	m_sfacMax = AskFloat("    Enter maximum wave vector modulus (nm^-1): [100] ", 100.0f);
	m_sfacRes = AskUnsignedInteger("    Enter Structure Factor resolution: [1000] ", 1000);
	
	mprintf("\n    The following normalization factors are available:\n");
	mprintf("    (1) [Sum_i x_i*f_i(q)]^2\n");
	mprintf("    (2) Sum_i x_i*[f_i(q)]^2\n");
	mprintf("    (3) 1 (no normalization)\n\n");
	m_normalization = AskRangeInteger("    Which factor to use? [1] ", 1, 3, 1);
	mprintf("\n");
	
	if (g_bAdvanced2) {
		m_saveIntermediate = AskYesNo("    Save all intermediate data (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		m_saveIntermediate = false;
	}
	
	if (g_bAdvanced2) {
		m_sharpening = AskYesNo("    Apply sharpening factor (y/n)? [no] ", false);
		if (m_sharpening) {
			while (true) {
				AskString_ND("    Which isotope to use as sharpening atom (e.g. \"N\" or \"14N\")? ", &buf);
				int isotopeIndex = -1;
				for (i = 0; i < g_isotopes.GetSize(); i++) {
					if(mystricmp(((CIsotope *)g_isotopes[i])->label(), buf) == 0) {
						isotopeIndex = i;
						break;
					}
				}
				if(isotopeIndex != -1) {
					m_sharpeningIsotope = (CIsotope *)g_isotopes[isotopeIndex];
					break;
				}
				eprintf("Isotope data for \"%s\" not found.\n", (const char *)buf);
			}
		} else {
			m_sharpeningIsotope = NULL;
		}
		mprintf("\n");
	} else {
		m_sharpening = false;
		m_sharpeningIsotope = NULL;
	}
}

CStructureFactor::~CStructureFactor() {
	int i;
	delete m_globalSFac;
	for (i = 0; i < m_sFacGroups.GetSize(); i++)
		delete (CStructureFactorGroup *)m_sFacGroups[i];
	for (i = 0; i < m_sFacCrosses.GetSize(); i++)
		delete (CStructureFactorCross *)m_sFacCrosses[i];
}

void CStructureFactor::initialize() {
	int i;
	mprintf("  Initializing structure factor\n");
	mprintf(WHITE, "\n    > Total >\n\n");
	m_globalSFac->initialize(m_rdfMax, m_rdfRes);
	mprintf(WHITE, "\n    < End of Total <\n");
	for (i = 0; i < m_sFacGroups.GetSize(); i++) {
		mprintf(WHITE, "\n    > Group %d >\n\n", i + 1);
		((CStructureFactorGroup *)m_sFacGroups[i])->initialize(m_rdfMax, m_rdfRes);
		mprintf(WHITE, "\n    < End of Group %d <\n", i + 1);
	}
	for (i = 0; i < m_sFacCrosses.GetSize(); i++) {
		mprintf(WHITE, "\n    > Cross Term %d >\n\n", i + 1);
		((CStructureFactorCross *)m_sFacCrosses[i])->initialize(m_rdfMax, m_rdfRes);
		mprintf(WHITE, "\n    < End of Cross Term %d <\n", i + 1);
	}
}

void CStructureFactor::process(CTimeStep *ts) {
	int i;
	m_globalSFac->process(ts);
	for (i = 0; i < m_sFacGroups.GetSize(); i++) {
		((CStructureFactorGroup *)m_sFacGroups[i])->process(ts);
	}
	for (i = 0; i < m_sFacCrosses.GetSize(); i++) {
		((CStructureFactorCross *)m_sFacCrosses[i])->process(ts);
	}
}

void CStructureFactor::finalize() {
	int i;
	CxString filename;
	CDF *sfac, *xray, *neutron;
	try { sfac = new CDF(); } catch(...) { sfac = NULL; }
	if (sfac == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfac->m_fMinVal = 2.0f * Pi / m_rdfMax * 1000.0f;
	sfac->m_fMaxVal = m_sfacMax;
	sfac->m_iResolution = m_sfacRes;
	sfac->SetLabelX("Wave vector modulus / nm^-1");
	sfac->SetLabelY("Intensity");
	sfac->Create();
	try { xray = new CDF(); } catch(...) { xray = NULL; }
	if (xray == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	xray->m_fMinVal = 2.0f * Pi / m_rdfMax * 1000.0f;
	xray->m_fMaxVal = m_sfacMax;
	xray->m_iResolution = m_sfacRes;
	sfac->SetLabelX("Wave vector modulus / nm^-1");
	sfac->SetLabelY("Intensity");
	xray->Create();
	try { neutron = new CDF(); } catch(...) { neutron = NULL; }
	if (neutron == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	neutron->m_fMinVal = 2.0f * Pi / m_rdfMax * 1000.0f;
	neutron->m_fMaxVal = m_sfacMax;
	neutron->m_iResolution = m_sfacRes;
	sfac->SetLabelX("Wave vector modulus / nm^-1");
	sfac->SetLabelY("Intensity");
	neutron->Create();
	
	mprintf("  Finalizing structure factor\n");
	
	mprintf(WHITE, "\n    > Total >\n\n");
	
	m_globalSFac->finalize(sfac, xray, neutron, m_saveIntermediate);
	if (g_bBoxNonOrtho) {
		float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxVolume * 1000000.0f);
		sfac->MultiplyBin(f);
		xray->MultiplyBin(f);
		neutron->MultiplyBin(f);
	} else {
		float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ);
		sfac->MultiplyBin(f);
		xray->MultiplyBin(f);
		neutron->MultiplyBin(f);
	}
	
	if (m_normalization == 1) {
		float sf = 0.0f;
		for (i = 0; i < m_globalSFac->getIsotopeTypeList().GetSize(); i++) {
			sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(i) / g_iGesAtomCount;
		}
		sfac->MultiplyBin(1.0f / (sf * sf));
		
		for (i = 0; i < xray->m_iResolution; i++) {
			float q = (xray->m_fMinVal + (0.5f + i) / xray->m_fFac) / 1000.0f;
			float xf = 0.0f;
			int j;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->xrayFactor(q);
			}
			xray->m_pBin[i] /= xf * xf;
		}
		
		float nf = 0.0f;
		for (i = 0; i < m_globalSFac->getIsotopeTypeList().GetSize(); i++) {
			nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(i) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(i))->neutronFactor();
		}
		neutron->MultiplyBin(1.0f / (nf * nf));
	} else if (m_normalization == 2) {
		float sf = 0.0f;
		for (i = 0; i < m_globalSFac->getIsotopeTypeList().GetSize(); i++) {
			sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(i) / g_iGesAtomCount;
		}
		sfac->MultiplyBin(1.0f / sf);
		
		for (i = 0; i < xray->m_iResolution; i++) {
			float q = (xray->m_fMinVal + (0.5f + i) / xray->m_fFac) / 1000.0f;
			float xf = 0.0f;
			int j;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->xrayFactor(q) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->xrayFactor(q);
			}
			xray->m_pBin[i] /= xf;
		}
		
		float nf = 0.0f;
		for (i = 0; i < m_globalSFac->getIsotopeTypeList().GetSize(); i++) {
			nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(i) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(i))->neutronFactor() * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(i))->neutronFactor();
		}
		neutron->MultiplyBin(1.0f / nf);
	} else if (m_normalization != 3) {
		eprintf("Weird error. Unknown normalization factor.\n");
		abort();
	}
	
	mprintf("\n");
	filename.sprintf("sfac_total_unweighted.csv");
	mprintf("    Writing unweighted structure factor to %s...\n", (const char *)filename);
	sfac->Write("", filename, "", false);
	filename.sprintf("sfac_total_xray.csv");
	mprintf("    Writing X-ray scattering function to %s...\n", (const char *)filename);
	xray->Write("", filename, "", false);
	filename.sprintf("sfac_total_neutron.csv");
	mprintf("    Writing neutron scattering function to %s...\n", (const char *)filename);
	neutron->Write("", filename, "", false);
	
	if (m_sharpening) {
		for (i = 0; i < xray->m_iResolution; i++) {
			float q = (xray->m_fMinVal + (0.5f + i) / xray->m_fFac) / 1000.0f;
			float nf = q * expf(-0.01f * q * q);
			float xf = q * expf(-0.01f * q * q) * m_sharpeningIsotope->xrayFactor(0.0f) * m_sharpeningIsotope->xrayFactor(0.0f) / m_sharpeningIsotope->xrayFactor(q) / m_sharpeningIsotope->xrayFactor(q);
			neutron->m_pBin[i] *= nf;
			xray->m_pBin[i] *= xf;
		}
		
		filename.sprintf("sfac_total_neutron_sharpened.csv");
		mprintf("    Writing sharpened neutron scattering function to %s...\n", (const char *)filename);
		neutron->Write("", filename, "", false);
		filename.sprintf("sfac_total_xray_sharpened.csv");
		mprintf("    Writing sharpened X-ray scattering function to %s...\n", (const char *)filename);
		xray->Write("", filename, "", false);
	}
	
	mprintf(WHITE, "\n    < End of Total <\n");
	
	for (i = 0; i < m_sFacGroups.GetSize(); i++) {
		mprintf(WHITE, "\n    > Group %d >\n\n", i + 1);
		
		sfac->ZeroBin();
		xray->ZeroBin();
		neutron->ZeroBin();
		
		((CStructureFactorGroup *)m_sFacGroups[i])->finalize(sfac, xray, neutron, m_saveIntermediate);
		if (g_bBoxNonOrtho) {
			float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxVolume * 1000000.0f);
			sfac->MultiplyBin(f);
			xray->MultiplyBin(f);
			neutron->MultiplyBin(f);
		} else {
			float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ);
			sfac->MultiplyBin(f);
			xray->MultiplyBin(f);
			neutron->MultiplyBin(f);
		}
		
		if (m_normalization == 1) {
			float sf = 0.0f;
			int j;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
			}
			sfac->MultiplyBin(1.0f / (sf * sf));
			
			for (j = 0; j < xray->m_iResolution; j++) {
				float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
				float xf = 0.0f;
				int k;
				for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
					xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
				}
				xray->m_pBin[j] /= xf * xf;
			}
			
			float nf = 0.0f;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
			}
			neutron->MultiplyBin(1.0f / (nf * nf));
		} else if (m_normalization == 2) {
			float sf = 0.0f;
			int j;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
			}
			sfac->MultiplyBin(1.0f / sf);
			
			for (j = 0; j < xray->m_iResolution; j++) {
				float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
				float xf = 0.0f;
				int k;
				for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
					xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
				}
				xray->m_pBin[j] /= xf;
			}
			
			float nf = 0.0f;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor() * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
			}
			neutron->MultiplyBin(1.0f / nf);
		} else if (m_normalization != 3) {
			eprintf("Weird error. Unknown normalization factor.\n");
			abort();
		}
		
		mprintf("\n");
		filename.sprintf("sfac_self_%s_unweighted.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
		mprintf("    Writing unweighted structure factor to %s...\n", (const char *)filename);
		sfac->Write("", filename, "", false);
		filename.sprintf("sfac_self_%s_xray.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
		mprintf("    Writing X-ray scattering function to %s...\n", (const char *)filename);
		xray->Write("", filename, "", false);
		filename.sprintf("sfac_self_%s_neutron.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
		mprintf("    Writing neutron scattering function to %s...\n", (const char *)filename);
		neutron->Write("", filename, "", false);
		
		if (m_sharpening) {
			int j;
			for (j = 0; j < xray->m_iResolution; j++) {
				float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
				float nf = q * expf(-0.01f * q * q);
				float xf = q * expf(-0.01f * q * q) * m_sharpeningIsotope->xrayFactor(0.0f) * m_sharpeningIsotope->xrayFactor(0.0f) / m_sharpeningIsotope->xrayFactor(q) / m_sharpeningIsotope->xrayFactor(q);
				neutron->m_pBin[j] *= nf;
				xray->m_pBin[j] *= xf;
			}
			
			filename.sprintf("sfac_self_%s_neutron_sharpened.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			mprintf("    Writing sharpened neutron scattering function to %s...\n", (const char *)filename);
			neutron->Write("", filename, "", false);
			filename.sprintf("sfac_self_%s_xray_sharpened.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			mprintf("    Writing sharpened X-ray scattering function to %s...\n", (const char *)filename);
			xray->Write("", filename, "", false);
		}
		
		if (((CStructureFactorGroup *)m_sFacGroups[i])->sepInterIntra()) {
			sfac->ZeroBin();
			xray->ZeroBin();
			neutron->ZeroBin();
			
			mprintf("\n");
			((CStructureFactorGroup *)m_sFacGroups[i])->finalizeIntra(sfac, xray, neutron, m_saveIntermediate);
			if (g_bBoxNonOrtho) {
				float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxVolume * 1000000.0f);
				sfac->MultiplyBin(f);
				xray->MultiplyBin(f);
				neutron->MultiplyBin(f);
			} else {
				float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ);
				sfac->MultiplyBin(f);
				xray->MultiplyBin(f);
				neutron->MultiplyBin(f);
			}
			
			if (m_normalization == 1) {
				float sf = 0.0f;
				int j;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
				}
				sfac->MultiplyBin(1.0f / (sf * sf));
				
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float xf = 0.0f;
					int k;
					for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
						xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
					}
					xray->m_pBin[j] /= xf * xf;
				}
				
				float nf = 0.0f;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
				}
				neutron->MultiplyBin(1.0f / (nf * nf));
			} else if (m_normalization == 2) {
				float sf = 0.0f;
				int j;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
				}
				sfac->MultiplyBin(1.0f / sf);
				
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float xf = 0.0f;
					int k;
					for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
						xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
					}
					xray->m_pBin[j] /= xf;
				}
				
				float nf = 0.0f;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor() * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
				}
				neutron->MultiplyBin(1.0f / nf);
			} else if (m_normalization != 3) {
				eprintf("Weird error. Unknown normalization factor.\n");
				abort();
			}
			
			mprintf("\n");
			filename.sprintf("sfac_self_%s_unweighted_intra.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			mprintf("    Writing unweighted intramolecular structure factor to %s...\n", (const char *)filename);
			sfac->Write("", filename, "", false);
			filename.sprintf("sfac_self_%s_xray_intra.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			mprintf("    Writing intramolecular X-ray scattering function to %s...\n", (const char *)filename);
			xray->Write("", filename, "", false);
			filename.sprintf("sfac_self_%s_neutron_intra.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			mprintf("    Writing intramolecular neutron scattering function to %s...\n", (const char *)filename);
			neutron->Write("", filename, "", false);
			
			if (m_sharpening) {
				int j;
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float nf = q * expf(-0.01f * q * q);
					float xf = q * expf(-0.01f * q * q) * m_sharpeningIsotope->xrayFactor(0.0f) * m_sharpeningIsotope->xrayFactor(0.0f) / m_sharpeningIsotope->xrayFactor(q) / m_sharpeningIsotope->xrayFactor(q);
					neutron->m_pBin[j] *= nf;
					xray->m_pBin[j] *= xf;
				}
				
				filename.sprintf("sfac_self_%s_neutron_intra_sharpened.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
				mprintf("    Writing intramolecular sharpened neutron scattering function to %s...\n", (const char *)filename);
				neutron->Write("", filename, "", false);
				filename.sprintf("sfac_self_%s_xray_intra_sharpened.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
				mprintf("    Writing intramolecular sharpened X-ray scattering function to %s...\n", (const char *)filename);
				xray->Write("", filename, "", false);
			}
			
			sfac->ZeroBin();
			xray->ZeroBin();
			neutron->ZeroBin();
			
			mprintf("\n");
			((CStructureFactorGroup *)m_sFacGroups[i])->finalizeInter(sfac, xray, neutron, m_saveIntermediate);
			if (g_bBoxNonOrtho) {
				float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxVolume * 1000000.0f);
				sfac->MultiplyBin(f);
				xray->MultiplyBin(f);
				neutron->MultiplyBin(f);
			} else {
				float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ);
				sfac->MultiplyBin(f);
				xray->MultiplyBin(f);
				neutron->MultiplyBin(f);
			}
			
			if (m_normalization == 1) {
				float sf = 0.0f;
				int j;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
				}
				sfac->MultiplyBin(1.0f / (sf * sf));
				
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float xf = 0.0f;
					int k;
					for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
						xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
					}
					xray->m_pBin[j] /= xf * xf;
				}
				
				float nf = 0.0f;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
				}
				neutron->MultiplyBin(1.0f / (nf * nf));
			} else if (m_normalization == 2) {
				float sf = 0.0f;
				int j;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
				}
				sfac->MultiplyBin(1.0f / sf);
				
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float xf = 0.0f;
					int k;
					for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
						xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
					}
					xray->m_pBin[j] /= xf;
				}
				
				float nf = 0.0f;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor() * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
				}
				neutron->MultiplyBin(1.0f / nf);
			} else if (m_normalization != 3) {
				eprintf("Weird error. Unknown normalization factor.\n");
				abort();
			}
			
			mprintf("\n");
			filename.sprintf("sfac_self_%s_unweighted_inter.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			mprintf("    Writing unweighted intermolecular structure factor to %s...\n", (const char *)filename);
			sfac->Write("", filename, "", false);
			filename.sprintf("sfac_self_%s_xray_inter.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			mprintf("    Writing intermolecular X-ray scattering function to %s...\n", (const char *)filename);
			xray->Write("", filename, "", false);
			filename.sprintf("sfac_self_%s_neutron_inter.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			mprintf("    Writing intermolecular neutron scattering function to %s...\n", (const char *)filename);
			neutron->Write("", filename, "", false);
			
			if (m_sharpening) {
				int j;
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float nf = q * expf(-0.01f * q * q);
					float xf = q * expf(-0.01f * q * q) * m_sharpeningIsotope->xrayFactor(0.0f) * m_sharpeningIsotope->xrayFactor(0.0f) / m_sharpeningIsotope->xrayFactor(q) / m_sharpeningIsotope->xrayFactor(q);
					neutron->m_pBin[j] *= nf;
					xray->m_pBin[j] *= xf;
				}
				
				filename.sprintf("sfac_self_%s_neutron_inter_sharpened.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
				mprintf("    Writing intermolecular sharpened neutron scattering function to %s...\n", (const char *)filename);
				neutron->Write("", filename, "", false);
				filename.sprintf("sfac_self_%s_xray_inter_sharpened.csv", (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
				mprintf("    Writing intermolecular sharpened X-ray scattering function to %s...\n", (const char *)filename);
				xray->Write("", filename, "", false);
			}
		}
		
		mprintf(WHITE, "\n    < End of Group %d <\n", i + 1);
	}
	
	for (i = 0; i < m_sFacCrosses.GetSize(); i++) {
		mprintf(WHITE, "\n    > Cross Term %d >\n\n", i + 1);
		
		sfac->ZeroBin();
		xray->ZeroBin();
		neutron->ZeroBin();
		
		((CStructureFactorCross *)m_sFacCrosses[i])->finalize(sfac, xray, neutron, m_saveIntermediate);
		if (g_bBoxNonOrtho) {
			float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxVolume * 1000000.0f);
			sfac->MultiplyBin(f);
			xray->MultiplyBin(f);
			neutron->MultiplyBin(f);
		} else {
			float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ);
			sfac->MultiplyBin(f);
			xray->MultiplyBin(f);
			neutron->MultiplyBin(f);
		}
		
		if (m_normalization == 1) {
			float sf = 0.0f;
			int j;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
			}
			sfac->MultiplyBin(1.0f / (sf * sf));
			
			for (j = 0; j < xray->m_iResolution; j++) {
				float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
				float xf = 0.0f;
				int k;
				for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
					xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
				}
				xray->m_pBin[j] /= xf * xf;
			}
			
			float nf = 0.0f;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
			}
			neutron->MultiplyBin(1.0f / (nf * nf));
		} else if (m_normalization == 2) {
			float sf = 0.0f;
			int j;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
			}
			sfac->MultiplyBin(1.0f / sf);
			
			for (j = 0; j < xray->m_iResolution; j++) {
				float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
				float xf = 0.0f;
				int k;
				for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
					xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
				}
				xray->m_pBin[j] /= xf;
			}
			
			float nf = 0.0f;
			for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
				nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor() * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
			}
			neutron->MultiplyBin(1.0f / nf);
		} else if (m_normalization != 3) {
			eprintf("Weird error. Unknown normalization factor.\n");
			abort();
		}
		
		mprintf("\n");
		filename.sprintf("sfac_cross_%s_unweighted.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
		mprintf("    Writing unweighted structure factor to %s...\n", (const char *)filename);
		sfac->Write("", filename, "", false);
		filename.sprintf("sfac_cross_%s_xray.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
		mprintf("    Writing X-ray scattering function to %s...\n", (const char *)filename);
		xray->Write("", filename, "", false);
		filename.sprintf("sfac_cross_%s_neutron.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
		mprintf("    Writing neutron scattering function to %s...\n", (const char *)filename);
		neutron->Write("", filename, "", false);
		
		if (m_sharpening) {
			int j;
			for (j = 0; j < xray->m_iResolution; j++) {
				float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
				float nf = q * expf(-0.01f * q * q);
				float xf = q * expf(-0.01f * q * q) * m_sharpeningIsotope->xrayFactor(0.0f) * m_sharpeningIsotope->xrayFactor(0.0f) / m_sharpeningIsotope->xrayFactor(q) / m_sharpeningIsotope->xrayFactor(q);
				neutron->m_pBin[j] *= nf;
				xray->m_pBin[j] *= xf;
			}
			
			filename.sprintf("sfac_cross_%s_neutron_sharpened.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
			mprintf("    Writing sharpened neutron scattering function to %s...\n", (const char *)filename);
			neutron->Write("", filename, "", false);
			filename.sprintf("sfac_cross_%s_xray_sharpened.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
			mprintf("    Writing sharpened X-ray scattering function to %s...\n", (const char *)filename);
			xray->Write("", filename, "", false);
		}
		
		if (((CStructureFactorCross *)m_sFacCrosses[i])->sepInterIntra()) {
			sfac->ZeroBin();
			xray->ZeroBin();
			neutron->ZeroBin();
			
			mprintf("\n");
			((CStructureFactorCross *)m_sFacCrosses[i])->finalizeIntra(sfac, xray, neutron, m_saveIntermediate);
			if (g_bBoxNonOrtho) {
				float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxVolume * 1000000.0f);
				sfac->MultiplyBin(f);
				xray->MultiplyBin(f);
				neutron->MultiplyBin(f);
			} else {
				float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ);
				sfac->MultiplyBin(f);
				xray->MultiplyBin(f);
				neutron->MultiplyBin(f);
			}
			
			if (m_normalization == 1) {
				float sf = 0.0f;
				int j;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
				}
				sfac->MultiplyBin(1.0f / (sf * sf));
				
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float xf = 0.0f;
					int k;
					for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
						xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
					}
					xray->m_pBin[j] /= xf * xf;
				}
				
				float nf = 0.0f;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
				}
				neutron->MultiplyBin(1.0f / (nf * nf));
			} else if (m_normalization == 2) {
				float sf = 0.0f;
				int j;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
				}
				sfac->MultiplyBin(1.0f / sf);
				
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float xf = 0.0f;
					int k;
					for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
						xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
					}
					xray->m_pBin[j] /= xf;
				}
				
				float nf = 0.0f;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor() * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
				}
				neutron->MultiplyBin(1.0f / nf);
			} else if (m_normalization != 3) {
				eprintf("Weird error. Unknown normalization factor.\n");
				abort();
			}
			
			mprintf("\n");
			filename.sprintf("sfac_cross_%s_unweighted_intra.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
			mprintf("    Writing unweighted intramolecular structure factor to %s...\n", (const char *)filename);
			sfac->Write("", filename, "", false);
			filename.sprintf("sfac_cross_%s_xray_intra.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
			mprintf("    Writing intramolecular X-ray scattering function to %s...\n", (const char *)filename);
			xray->Write("", filename, "", false);
			filename.sprintf("sfac_cross_%s_neutron_intra.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
			mprintf("    Writing intramolecular neutron scattering function to %s...\n", (const char *)filename);
			neutron->Write("", filename, "", false);
			
			if (m_sharpening) {
				int j;
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float nf = q * expf(-0.01f * q * q);
					float xf = q * expf(-0.01f * q * q) * m_sharpeningIsotope->xrayFactor(0.0f) * m_sharpeningIsotope->xrayFactor(0.0f) / m_sharpeningIsotope->xrayFactor(q) / m_sharpeningIsotope->xrayFactor(q);
					neutron->m_pBin[j] *= nf;
					xray->m_pBin[j] *= xf;
				}
				
				filename.sprintf("sfac_cross_%s_neutron_intra_sharpened.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
				mprintf("    Writing sharpened intramolecular neutron scattering function to %s...\n", (const char *)filename);
				neutron->Write("", filename, "", false);
				filename.sprintf("sfac_cross_%s_xray_intra_sharpened.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
				mprintf("    Writing sharpened intramolcular X-ray scattering function to %s...\n", (const char *)filename);
				xray->Write("", filename, "", false);
			}
			
			sfac->ZeroBin();
			xray->ZeroBin();
			neutron->ZeroBin();
			
			mprintf("\n");
			((CStructureFactorCross *)m_sFacCrosses[i])->finalizeInter(sfac, xray, neutron, m_saveIntermediate);
			if (g_bBoxNonOrtho) {
				float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxVolume * 1000000.0f);
				sfac->MultiplyBin(f);
				xray->MultiplyBin(f);
				neutron->MultiplyBin(f);
			} else {
				float f = 4.0f * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ);
				sfac->MultiplyBin(f);
				xray->MultiplyBin(f);
				neutron->MultiplyBin(f);
			}
			
			if (m_normalization == 1) {
				float sf = 0.0f;
				int j;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
				}
				sfac->MultiplyBin(1.0f / (sf * sf));
				
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float xf = 0.0f;
					int k;
					for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
						xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
					}
					xray->m_pBin[j] /= xf * xf;
				}
				
				float nf = 0.0f;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
				}
				neutron->MultiplyBin(1.0f / (nf * nf));
			} else if (m_normalization == 2) {
				float sf = 0.0f;
				int j;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					sf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount;
				}
				sfac->MultiplyBin(1.0f / sf);
				
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float xf = 0.0f;
					int k;
					for (k = 0; k < m_globalSFac->getIsotopeTypeList().GetSize(); k++) {
						xf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(k) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(k))->xrayFactor(q);
					}
					xray->m_pBin[j] /= xf;
				}
				
				float nf = 0.0f;
				for (j = 0; j < m_globalSFac->getIsotopeTypeList().GetSize(); j++) {
					nf += (float)m_globalSFac->getIsotopeTypeTotalCount().GetAt(j) / g_iGesAtomCount * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor() * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(j))->neutronFactor();
				}
				neutron->MultiplyBin(1.0f / nf);
			} else if (m_normalization != 3) {
				eprintf("Weird error. Unknown normalization factor.\n");
				abort();
			}
			
			mprintf("\n");
			filename.sprintf("sfac_cross_%s_unweighted_inter.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
			mprintf("    Writing unweighted intermolecular structure factor to %s...\n", (const char *)filename);
			sfac->Write("", filename, "", false);
			filename.sprintf("sfac_cross_%s_xray_inter.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
			mprintf("    Writing intermolecular X-ray scattering function to %s...\n", (const char *)filename);
			xray->Write("", filename, "", false);
			filename.sprintf("sfac_cross_%s_neutron_inter.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
			mprintf("    Writing intermolecular neutron scattering function to %s...\n", (const char *)filename);
			neutron->Write("", filename, "", false);
			
			if (m_sharpening) {
				int j;
				for (j = 0; j < xray->m_iResolution; j++) {
					float q = (xray->m_fMinVal + (0.5f + j) / xray->m_fFac) / 1000.0f;
					float nf = q * expf(-0.01f * q * q);
					float xf = q * expf(-0.01f * q * q) * m_sharpeningIsotope->xrayFactor(0.0f) * m_sharpeningIsotope->xrayFactor(0.0f) / m_sharpeningIsotope->xrayFactor(q) / m_sharpeningIsotope->xrayFactor(q);
					neutron->m_pBin[j] *= nf;
					xray->m_pBin[j] *= xf;
				}
				
				filename.sprintf("sfac_cross_%s_neutron_inter_sharpened.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
				mprintf("    Writing sharpened intermolecular neutron scattering function to %s...\n", (const char *)filename);
				neutron->Write("", filename, "", false);
				filename.sprintf("sfac_cross_%s_xray_inter_sharpened.csv", (const char *)((CStructureFactorCross *)m_sFacCrosses[i])->getName());
				mprintf("    Writing sharpened intermolcular X-ray scattering function to %s...\n", (const char *)filename);
				xray->Write("", filename, "", false);
			}
		}
		
		mprintf(WHITE, "\n    < End of Cross Term %d <\n", i + 1);
	}
	
	delete sfac;
	delete xray;
	delete neutron;
}

bool gatherStructureFactor() {
	createIsotopeList();
	
	try { g_structureFactor = new CStructureFactor(); } catch(...) { g_structureFactor = NULL; }
	if (g_structureFactor == NULL) NewException((double)sizeof(CStructureFactor), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
// 	if(AskYesNo("    Compute structure factor of whole system (y/n)? [yes] ", true)) {
// 		mprintf(YELLOW, "\n>>> Structure Factor Observation 1 >>>\n\n");
// 		
// 		CSFacObservation *obs;
// 		try { obs = new CSFacObservation(true); } catch(...) { obs = NULL; }
// 		if(obs == NULL) NewException((double)sizeof(CSFacObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		g_sFacObserv.Add(obs);
// 		
// 		mprintf(YELLOW, "<<< End of Structure Factor Observation 1 <<<\n\n");
// 	}
// 	
// 	if(AskYesNo("    Compute structure factors for specific atoms (y/n)? [no] ", false)) {
// 		while(true) {
// 			mprintf(YELLOW, "\n>>> Structure Factor Observation %d >>>\n\n", g_sFacObserv.GetSize() + 1);
// 			
// 			CSFacObservation *obs;
// 			try { obs = new CSFacObservation(); } catch(...) { obs = NULL; }
// 			if(obs == NULL) NewException((double)sizeof(CSFacObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			g_sFacObserv.Add(obs);
// 			
// 			mprintf(YELLOW, "<<< End of Structure Factor Observation %d <<<\n\n", g_sFacObserv.GetSize());
// 			
// 			if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
// 				break;
// 			mprintf("\n");
// 		}
// 	}
// 	mprintf("\n");
	
	return true;
}

bool initializeStructureFactor() {
	g_structureFactor->initialize();
// 	int i;
// 	for(i = 0; i < g_sFacObserv.GetSize(); i++) {
// 		((CSFacObservation *)g_sFacObserv[i])->initialize();
// 	}
	return true;
}

void processStructureFactor(CTimeStep *ts) {
	g_structureFactor->process(ts);
// 	int i;
// 	for(i = 0; i < g_sFacObserv.GetSize(); i++) {
// 		((CSFacObservation *)g_sFacObserv[i])->process(ts);
// 	}
}

void finalizeStructureFactor() {
	g_structureFactor->finalize();
	delete g_structureFactor;
// 	int i;
// 	for(i = 0; i < g_sFacObserv.GetSize(); i++) {
// 		mprintf(YELLOW, ">>> Structure Factor Observation %d >>>\n\n", i + 1);
// 		((CSFacObservation *)g_sFacObserv[i])->finalize();
// 		mprintf(YELLOW, "\n<<< End of Structure Factor Observation %d <<<\n\n", i + 1);
// 	}
	deleteIsotopeList();
}

// CStructureFactor::CStructureFactor()
// {
// }
// 
// 
// CStructureFactor::~CStructureFactor()
// {
// }
// 
// 
// void CStructureFactor::Create()
// {
// 	int z, z2;
// 	CDF *df;
// 
// 	for (z=0;z<g_oaAtoms.GetSize();z++)
// 	{
// 		if (z == g_iVirtAtomType)
// 			continue;
// 
// 		if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fCoherentNCS == 0)
// 			eprintf("\nError: Coherent Neutron Cross Section for %s not defined. Edit elementdata.cpp.\n",((CAtom*)g_oaAtoms[z])->m_sName);
// 
// 		for (z2=0;z2<g_oaAtoms.GetSize();z2++)
// 		{
// 			if (z2 == g_iVirtAtomType)
// 				continue;
// 
// 			if (z2 < z)
// 			{
// 				m_oaRDFs.Add(NULL);
// 				m_oaSFacs.Add(NULL);
// 			} else
// 			{
// 				df = new CDF();
// 				df->m_fMinVal = 0;
// 				df->m_fMaxVal = m_fRDFRange;
// 				df->m_iResolution = m_iRDFRes;
// 				df->SetLabelX("Distance [pm]");
// 				df->SetLabelY("g(r)");
// 				df->Create();
// 				m_oaRDFs.Add(df);
// 
// 				df = new CDF();
// 				df->m_fMinVal = 0;
// 				df->m_fMaxVal = m_fSQRange;
// 				df->m_iResolution = m_iSQRes;
// 				df->SetLabelX("Wave Vector Modulus [Angstrom^-1]");
// 				df->SetLabelY("S(k)");
// 				df->Create();
// 				m_oaSFacs.Add(df);
// 			}
// 		}
// 	}
// }
// 
// 
// void CStructureFactor::Parse()
// {
// 	mprintf(WHITE,">>> Structure Factor Analysis >>>\n\n");
// 
// 	m_fRDFRange = AskFloat("    Enter maximum RDF distance to observe (in pm): [%d] ",(float)HalfBox(),HalfBox());
// 	m_iRDFRes = AskUnsignedInteger("    Enter RDF binning resolution: [300] ",300);
// 	m_fSQRange = AskFloat("    Enter maximum wave vector modulus (in Angstrom^-1): [50] ",50);
// 	m_iSQRes = AskUnsignedInteger("    Enter Structure Factor binning resolution: [2000] ",2000);
// 
// 	mprintf("\n");
// 
// 	m_bDumpElementRDFs = AskYesNo("    Write out all element RDFs (y/n)? [yes] ",true);
// 	m_bDumpTotalRDF = AskYesNo("    Write out overall RDF (y/n)? [yes] ",true);
// 	m_bDumpElementSFac = AskYesNo("    Write out all element Structure Factor contributions (y/n)? [yes] ",true);
// 
// 	mprintf(WHITE,"\n<<< End of Structure Factor Analysis <<<\n\n");
// }
// 
// 
// void CStructureFactor::ProcessStep(CTimeStep *ts)
// {
// 	int z, z2, t1, t2, i;
// 	float v;
// 
// 	i = g_oaAtoms.GetSize()-1;
// 
// 	for (z=0;z<g_iGesAtomCount;z++)
// 	{
// 		t1 = g_waAtomRealElement[z];
// 		for (z2=z+1;z2<g_iGesAtomCount;z2++)
// 		{
// 			t2 = g_waAtomRealElement[z2];
// 
// 			v = FoldedLength(ts->m_vaCoords[z] - ts->m_vaCoords[z2]);
// 
// 			if (t1 >= t2)
// 			{
// 				((CDF*)m_oaRDFs[t2*i+t1])->AddToBin(v);
// 			} else
// 			{
// 				((CDF*)m_oaRDFs[t1*i+t2])->AddToBin(v);
// 			}
// 		}
// 	}
// }
// 
// 
// void CStructureFactor::Finish()
// {
// 	int z, z2, z3;
// 	CDF *df, *dftotal, *dfstotal;
// 	double fac/*, tf*/;
// 	char buf[256];
// 
// 	dftotal = NULL;
// 
// 	mprintf(WHITE,"* Structure Factor\n");
// 
// 	if (m_bDumpTotalRDF)
// 	{
// 		dftotal = new CDF();
// 		dftotal->m_fMinVal = 0;
// 		dftotal->m_fMaxVal = m_fRDFRange;
// 		dftotal->m_iResolution = m_iRDFRes;
// 		dftotal->SetLabelX("Distance [pm]");
// 		dftotal->SetLabelY("g(r)");
// 		dftotal->Create();
// //		tf = 0;
// 	}
// 
// 	dfstotal = new CDF();
// 	dfstotal->m_fMinVal = 0;
// 	dfstotal->m_fMaxVal = m_fSQRange;
// 	dfstotal->m_iResolution = m_iSQRes;
// 	dfstotal->SetLabelX("Wave Vector Modulus [Angstrom^-1]");
// 	dfstotal->SetLabelY("S(k)");
// 	dfstotal->Create();
// 
// 	for (z=0;z<g_oaAtoms.GetSize();z++)
// 	{
// 		if (z == g_iVirtAtomType)
// 			continue;
// 
// 		for (z2=0;z2<g_oaAtoms.GetSize();z2++)
// 		{
// 			if (z2 == g_iVirtAtomType)
// 				continue;
// 
// 			if (z2 < z)
// 				continue;
// 
// 			df = (CDF*)m_oaRDFs[z*(g_oaAtoms.GetSize()-1)+z2];
// 
// 			mprintf(WHITE,"  Processing Contributions of type %s - %s\n",((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z2])->m_sName);
// 
// 			df->MultiplyBin(1.0 / g_iSteps);
// 			if (z == z2)
// 				df->MultiplyBin(2.0);
// 			df->CorrectRadialDist();
// 			df->MultiplyBin(g_fBoxX*g_fBoxY*g_fBoxZ / (4.0/3.0*Pi) / ((CAtom*)g_oaAtoms[z])->m_iCount / ((CAtom*)g_oaAtoms[z2])->m_iCount);
// 
// 			if (g_bDoubleBox)
// 				df->MultiplyBin(g_iDoubleBoxFactor);
// 
// 			if (m_bDumpElementRDFs)
// 			{
// 				sprintf(buf,"sfac_rdf_%s_%s.csv",((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z2])->m_sName);
// 				mprintf("    Writing RDF to %s ...\n",buf);
// 				df->Write("",buf,"",false);
// 			}
// 
// 			if (m_bDumpTotalRDF)
// 			{
// 				if (z == z2)
// 					fac = 1.0;
// 						else fac = 2.0;
// 				for (z3=0;z3<m_iRDFRes;z3++)
// 					dftotal->m_pBin[z3] += df->m_pBin[z3] * fac * (double)((CAtom*)g_oaAtoms[z])->m_iCount * ((CAtom*)g_oaAtoms[z2])->m_iCount / g_iGesAtomCount / g_iGesAtomCount;
// 	//			tf += fac * (double)((CAtom*)g_oaAtoms[z])->m_iCount * ((CAtom*)g_oaAtoms[z2])->m_iCount / g_iGesAtomCount / g_iGesAtomCount;
// 			}
// 
// 			TransformRDF(df,(CDF*)m_oaSFacs[z*(g_oaAtoms.GetSize()-1)+z2]);
// 
// 			if (z == z2)
// 				fac = 1.0;
// 					else fac = 2.0;
// 			for (z3=0;z3<m_iSQRes;z3++)
// 				dfstotal->m_pBin[z3] += ((CDF*)m_oaSFacs[z*(g_oaAtoms.GetSize()-1)+z2])->m_pBin[z3] * fac * (double)((CAtom*)g_oaAtoms[z])->m_pElement->m_fCoherentNCS * ((CAtom*)g_oaAtoms[z2])->m_pElement->m_fCoherentNCS * ((CAtom*)g_oaAtoms[z])->m_iCount * ((CAtom*)g_oaAtoms[z2])->m_iCount / g_iGesAtomCount / g_iGesAtomCount;
// 
// 			if (m_bDumpElementSFac)
// 			{
// 				sprintf(buf,"sfac_%s_%s.csv",((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z2])->m_sName);
// 				mprintf("    Writing Structure Factor Contribution to %s ...\n",buf);
// 			//	mprintf("    Correction Factor: %G\n",((CDF*)m_oaSFacs[z*(g_oaAtoms.GetSize()-1)+z2])->m_pBin[0]);
// 			/*	tf = 0;
// 				for (z3=0;z3<m_iSQRes;z3++)
// 				{
// 					if (((CDF*)m_oaSFacs[z*(g_oaAtoms.GetSize()-1)+z2])->m_pBin[z3] < tf)
// 						tf = ((CDF*)m_oaSFacs[z*(g_oaAtoms.GetSize()-1)+z2])->m_pBin[z3];
// 				}
// 				if (tf < 0)
// 				{
// 					for (z3=0;z3<m_iSQRes;z3++)
// 						((CDF*)m_oaSFacs[z*(g_oaAtoms.GetSize()-1)+z2])->m_pBin[z3] /= fabs(tf);
// 				}
// 				for (z3=0;z3<m_iSQRes;z3++)
// 					((CDF*)m_oaSFacs[z*(g_oaAtoms.GetSize()-1)+z2])->m_pBin[z3] += 1.0;*/
// 				((CDF*)m_oaSFacs[z*(g_oaAtoms.GetSize()-1)+z2])->Write("",buf,"",false);
// 			}
// 		}
// 	}
// 
// 	if (m_bDumpTotalRDF)
// 	{
// 		mprintf(WHITE,"  * Overall RDF\n");
// 	//	mprintf("    tf = %f.\n",tf);
// 		sprintf(buf,"sfac_rdf_total.csv");
// 		mprintf("    Writing RDF to %s ...\n",buf);
// 		dftotal->Write("",buf,"",false);
// 	}
// 
// 	mprintf(WHITE,"  * Overall S(k)\n");
// 	sprintf(buf,"sfac_total.csv");
// 	mprintf("    Writing S(k) to %s ...\n",buf);
// /*	tf = 0;
// 	for (z3=0;z3<m_iSQRes;z3++)
// 	{
// 		if (dfstotal->m_pBin[z3] < tf)
// 			tf = dfstotal->m_pBin[z3];
// 	}
// 	if (tf < 0)
// 	{
// 		for (z3=0;z3<m_iSQRes;z3++)
// 			dfstotal->m_pBin[z3] /= fabs(tf);
// 	}
// 	for (z3=0;z3<m_iSQRes;z3++)
// 		dfstotal->m_pBin[z3] += 1.0;*/
// 	dfstotal->Write("",buf,"",false);
// }
// 
// 
// void CStructureFactor::TransformRDF(CDF *pin, CDF *pout)
// {
// 	int z, z2;
// 	double tf, tf2;
// 
// 	for (z=0;z<m_iSQRes;z++)
// 	{
// 		tf = 0;
// 
// 		for (z2=0;z2<m_iRDFRes;z2++)
// 		{
// 			tf2 = ((z2+0.5)/m_iRDFRes*m_fRDFRange/100.0) * (pin->m_pBin[z2]-1.0) / ((z+0.5)/m_iSQRes*m_fSQRange) * sin(((z2+0.5)/m_iRDFRes*m_fRDFRange/100.0) * ((z+0.5)/m_iSQRes*m_fSQRange));
// 			tf += tf2;
// 	//		mprintf("\n  %.4fA * %.4f / %.4fA^-1 * sin( %.4fA^-1 * %.4fA ) = %.4f",((z2+0.5)/m_iRDFRes*m_fRDFRange/100.0),(pin->m_pBin[z2]-1.0),((z+0.5)/m_iSQRes*m_fRDFRange/100.0),((z2+0.5)/m_iRDFRes*m_fRDFRange/100.0),((z+0.5)/m_iSQRes*m_fRDFRange/100.0),tf2);
// 		}
// 
// 		pout->m_pBin[z] = tf*4.0*Pi/g_fBoxX/g_fBoxY/g_fBoxZ*1000000.0;
// 	//	mprintf("\n\n## 1 + 4*Pi*%.4f*%.4f = %.4f\n",tf,(1.0/m_iRDFRes*m_fRDFRange/10000.0),pout->m_pBin[z]);
// 	}
// }

