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

#ifndef CHIRAL_H
#define CHIRAL_H

#include "moltools.h"
#include "xobarray.h"

class CTimeStep;

class CChiralObservation: public CObservation
{
public:
	CChiralObservation();
	~CChiralObservation();
	
	void initialize();
	void process(CTimeStep *ts);
	void finalize();
	
	int swapAtom(int offset);
	
private:
	unsigned char _ty, _rty, _atom;
	CxObArray _molAtoms;
	CxObArray _orderArrays;
	CxObArray _timedev;
	
	bool _fixTraj;
	CxIntArray _swapOffsets[2];
	bool _fixToR;
};

bool gatherChiral();
bool initializeChiral();
void processChiral(CTimeStep *ts);
void finalizeChiral();

#endif
