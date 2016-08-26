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

#ifndef PDF_H
#define PDF_H

#include "moltools.h"
#include "xdf.h"

class CDF;
class CSingleMolecule;
class CTimeStep;
class CxIntArray;


class CPDF: public CXDF
{
public:
	CPDF(int showMol);
	~CPDF();
	
	void initialize(int showMolCount);
	void finalize();
	
	void buildAtomList(CSingleMolecule *mol, CxIntArray *array);
	void addToDF(float value);
	
private:
	CDF *_df;
	CAtomGroup *_ag;
	int _showAtomGes;
	float _minDist, _maxDist;
	bool _scaleUniform;
	
	void buildName();
};


class CPDFObservation: public CObservation
{
public:
	CPDFObservation();
	~CPDFObservation();
	
	void initialize();
	void process(CTimeStep *ts);
	void finalize();
	
private:
	CPDF *_pdf;
};


bool gatherPDF();
bool initializePDF();
void processPDF(CTimeStep *ts);
void finalizePDF();

#endif
