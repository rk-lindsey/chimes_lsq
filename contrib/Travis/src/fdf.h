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

#ifndef FDF_H
#define FDF_H

#include "moltools.h"

class CTimeStep;

class CFDFObservation: public CObservation {
public:
	CFDFObservation();
	
	void initialize();
	void process(CTimeStep *ts);
	void finalize();
	
private:
	CAtomGroup m_atoms;
	CxString m_name;
	CDF m_fdf;
	bool m_massWeighting;
	CxFloatArray m_masses;
};

bool gatherFDF();
bool initializeFDF();
void processFDF(CTimeStep *ts);
void finalizeFDF();

#endif
