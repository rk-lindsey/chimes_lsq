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

#ifndef REVSDF_H
#define REVSDF_H

#include "xobject.h"
#include "2df.h"
#include "xvec3array.h"
#include "atomgroup.h"
#include "3df.h"


class CSingleMolecule;


class CRevSDF : public CxObject
{
public:
	CRevSDF();
	~CRevSDF();

	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();
	void CreateRevSDF();

	bool m_bIntra;
	int m_iRefOrSec;
	char *m_sName;
	float m_fRadius;
	int m_iResolution;
	bool m_bMirrorY;
	bool m_bMirrorBond;
	int m_iShowMol;
	CxVec3Array *m_vaData; 
	C2DF *m_p2DF;
	C3DF<double> *m_pRevSDF;
	CAtomGroup m_oAtoms;
	double m_fParticleDensity;
	int m_iShowAtomGes;

	double m_fSecondAtomPosX;
	double m_fSecondAtomCount;

	bool m_bCorrectAngle;
	bool m_bCorrectRadial;
	bool m_bCreateRevSDF;
	int m_iRevSDFRes;
	bool m_bDrawAtoms;
	int m_iHistogramRes;
	CxByteArray *m_baDataEnabled; 
};


#endif
