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

#ifndef NBEXCHANGE_H
#define NBEXCHANGE_H

#include "xobject.h"
#include "xobarray.h"
#include "df.h"
#include "nbsearch.h"

class CConditionGroup;

class CNbExchangePair : public CxObject  
{
public:
	CNbExchangePair() { }
	~CNbExchangePair() { }

	unsigned int m_iShowMol;
	unsigned char m_iStillAlive;
	int m_iStart;
	int m_iEnd;

//	float m_fDistance;
//	unsigned short m_iNeighbor;
};


class CNbExchange : public CxObject  
{
public:
	CNbExchange();
	~CNbExchange();

	bool m_bDifferentCrit;
	CNbSearch *m_pNbSearch;
	int m_iTimeDepth;
	CConditionGroup **m_pNbExConditions;
};

#endif

