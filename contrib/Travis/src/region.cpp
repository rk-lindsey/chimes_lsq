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
 
#include "region.h" 
 
#include "globalvar.h" 
#include "maintools.h" 
#include "tools.h" 
#include "xobarray.h" 
#include "xvector3.h" 
 
static CxObArray g_regions; 
 
CRegion::CRegion() { 
//	char buf[256]; 
	CxString buf;

	mprintf("    You have to define an orthorhombic part of the simulation box\n    by entering minimum and maximum values along the x, y, and z axes.\n\n"); 
	_xmin = AskFloat("    Minimum x value in pm [0.0] ", 0.0f); 
	_xmax = AskFloat("    Maximum x value in pm [%.1f] ", g_fBoxX, g_fBoxX); 
	_ymin = AskFloat("    Minimum y value in pm [0.0] ", 0.0f); 
	_ymax = AskFloat("    Maximum y value in pm [%.1f] ", g_fBoxY, g_fBoxY); 
	_zmin = AskFloat("    Minimum z value in pm [0.0] ", 0.0f); 
	_zmax = AskFloat("    Maximum z value in pm [%.1f] ", g_fBoxZ, g_fBoxZ); 
	 
	try { _centerAtomTypes = new unsigned char[g_oaMolecules.GetSize()]; } catch(...) { _centerAtomTypes = NULL; } 
	if(_centerAtomTypes == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(unsigned char), __FILE__, __LINE__, __PRETTY_FUNCTION__); 
	try { _centerAtomRealTypes = new unsigned char[g_oaMolecules.GetSize()]; } catch(...) { _centerAtomRealTypes = NULL; } 
	if(_centerAtomRealTypes == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(unsigned char), __FILE__, __LINE__, __PRETTY_FUNCTION__); 
	try { _centerAtoms = new unsigned char[g_oaMolecules.GetSize()]; } catch(...) { _centerAtoms = NULL; } 
	if(_centerAtoms == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(unsigned char), __FILE__, __LINE__, __PRETTY_FUNCTION__); 
	 
	int i; 
	mprintf("\n    By default, the center of mass will be used to decide if a molecule is within the region.\n"); 
	if(AskYesNo("    Change this behavior (y/n)? [no] ", false)) { 
		for(i = 0; i < g_oaMolecules.GetSize(); i++) { 
			while(true) { 
				mprintf("    Which atom to use for %s (e.g. C1)? [#2] ", ((CMolecule*)g_oaMolecules[i])->m_sName); 
				inpprintf("! Which atom to use for %s (e.g. C1)? [#2]\n", ((CMolecule*)g_oaMolecules[i])->m_sName); 
				myget(&buf); 
				if(strlen(buf) == 0) { 
					if(!ParseAtom("#2", i, _centerAtomTypes[i], _centerAtomRealTypes[i], _centerAtoms[i])) { 
						eprintf("Weird error.\n"); 
						abort(); 
					} 
				} else { 
					if(ParseAtom(buf, i, _centerAtomTypes[i], _centerAtomRealTypes[i], _centerAtoms[i])) 
						break; 
				} 
			} 
		} 
	} else { 
		for(i = 0; i < g_oaMolecules.GetSize(); i++) { 
			if(!ParseAtom("#2", i, _centerAtomTypes[i], _centerAtomRealTypes[i], _centerAtoms[i])) { 
				eprintf("Weird error.\n"); 
				abort(); 
			} 
		} 
	} 
} 
 
CRegion::~CRegion() { 
	delete[] _centerAtomTypes; 
	delete[] _centerAtomRealTypes; 
	delete[] _centerAtoms; 
} 
 
bool CRegion::isInRegion(const CxVector3 &vec) { 
	if((vec[0] > _xmin) && (vec[0] < _xmax) && (vec[1] > _ymin) && (vec[1] < _ymax) && (vec[2] > _zmin) && (vec[2] < _zmax)) 
		return true; 
	return false; 
} 
 
bool gatherRegionAnalysis() { 
	mprintf("\n"); 
	mprintf(YELLOW, ">>> Region Definition >>>\n\n"); 
	mprintf("    You can specify as many regions as you want.\n    All parts of the simulation box not contained in these regions will be region 1.\n\n"); 
	 
	while(true) { 
		mprintf(YELLOW, ">>> Region %d >>>\n\n", g_regions.GetSize() + 2); 
		 
		CRegion *region; 
		try { region = new CRegion(); } catch(...) { region = NULL; } 
		if(region == NULL) NewException((double)sizeof(CRegion), __FILE__, __LINE__, __PRETTY_FUNCTION__); 
		g_regions.Add(region); 
		 
		mprintf(YELLOW, "\n<<< End Region %d <<<\n\n", g_regions.GetSize() + 1); 
		 
		if(!AskYesNo("    Add another region (y/n)? [no] ", false)) 
			break; 
		mprintf("\n"); 
	} 
	 
	mprintf(YELLOW, "\n<<< End Region Definition <<<\n\n"); 
	return true; 
} 
 
void processRegionAnalysis(CTimeStep *ts) { 
	int i, j; 
	for(i = 0; i < g_iaSMRegion.GetSize(); i++) { 
		g_iaSMRegion[i] = 1; 
	} 
	 
	for(i = 0; i < g_regions.GetSize(); i++) { 
		CRegion *region = (CRegion *)g_regions[i]; 
		for(j = 0; j < g_oaSingleMolecules.GetSize(); j++) { 
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[j]; 
			CxVector3 position = ts->m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[region->centerAtomType(sm->m_iMolType)])->GetAt(region->centerAtom(sm->m_iMolType))]; 
			if(region->isInRegion(position)) 
				g_iaSMRegion[j] = i + 2; 
		} 
	} 
} 
 
void finalizeRegionAnalysis() { 
	int i; 
	for(i = 0; i < g_regions.GetSize(); i++) { 
		delete g_regions[i]; 
	} 
	g_regions.RemoveAll(); 
} 
