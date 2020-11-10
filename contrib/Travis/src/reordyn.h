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

#ifndef REORDYN_H
#define REORDYN_H

#include "xobject.h"
#include "df.h"
#include "atomgroup.h"
#include "acf.h"


class CSingleMolecule;


class CReorDyn : public CxObject
{
public:
	void Finish_Part2(const char *s, const char *multibuf, int nmol);
	void Finish(const char *multibuf);
	void BuildAtomList(CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();
	void ParseSpec();
	CReorDyn();
	~CReorDyn();

	bool m_bOrtho;
	int m_iVecType; // 0 = Ort, 1 = Dipol, 2 = Geschw., 3 = Kraft
	int m_iCombinations;
	int m_iShowMol;
	char *m_sName;
	char *m_sShortName;
	int m_iDepth;
	int m_iStride;
	int m_iMolecules;
	double *m_pCount;
	bool m_bSpectrum;
	CACF *m_pACF;
	bool m_bCrossCor;

	CxObArray m_oaCache;

	CxObArray m_oaVectors; // Enthaelt CAtomGroups, 3 pro Vektor

	CDF *m_pRDyn;

	bool m_bLifetimeSpectrum;
	CxObArray m_oaLTSpectra;

	bool m_bLegendre2;
};


#endif
