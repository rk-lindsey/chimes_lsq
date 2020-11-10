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

#ifndef ATOMGROUP_H
#define ATOMGROUP_H

#include "xobject.h"
#include "xbytearray.h"
#include "tools.h"

class CMolecule;


class CAtomGroup : public CxObject  
{
public:
	void CopyFrom(CAtomGroup *p);
	void Reset();
	void SortAtoms();
	CAtomGroup();
	~CAtomGroup();

	CMolecule* m_pMolecule;
	CxByteArray m_baAtomType;
	CxByteArray m_baRealAtomType;
	CxByteArray m_baAllAtoms;
	CxObArray m_oaAtoms;       // Enthaelt CxWordArrays
	int m_iAtomGes;
	bool m_bAllAtoms;

	bool ParseAtoms(CMolecule *m, const char *s);
	void AddAtom(int atom, int num, bool all);
	void AddAllAtoms(CMolecule *m, bool virt);
	void BuildName();

	char *m_sName;
};

#endif 
