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


#ifndef PLPROJ_H
#define PLPROJ_H


#include "tools.h"
#include "xobject.h"
#include "2df.h"
#include "moltools.h"
#include "xvec3array.h"
#include "xdvec3array.h"
#include "atomgroup.h"
#include "xbytearray.h"


class CPlProj : public CxObject
{
public:
	CPlProj();
	~CPlProj();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();

	int m_iHistogramRes;
	int m_iResolution[2];
	bool m_bIntra;
	int m_iShowMol;
	double m_fMinVal[2];
	double m_fMaxVal[2];
	double m_fSliceBorder[2];
	double m_fAspectRatio;
	double m_fParticleDensity;


	char *m_sName;
	CxVec3Array *m_vaData; 
	C2DF *m_p2DF;
	CAtomGroup m_oAtoms;
	CAtomGroup m_oDrawAtoms;
	bool m_bDrawAtoms;
	bool m_bAverageAtomPos;
	CxDVec3Array m_vaAtomPos;
	int m_iAverageCounter;
	CxByteArray *m_baDataEnabled; 
	int m_iShowAtomGes;

	
};


#endif

