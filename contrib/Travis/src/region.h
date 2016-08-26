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
 
#ifndef REGION_H 
#define REGION_H 
 
#include "xobject.h" 
 
class CTimeStep; 
class CxVector3; 
 

class CRegion: public CxObject 
{ 
public: 
	CRegion(); 
	~CRegion(); 
	 
	unsigned char centerAtomType(int index) const { return _centerAtomTypes[index]; } 
	unsigned char centerAtomRealType(int index) const { return _centerAtomRealTypes[index]; } 
	unsigned char centerAtom(int index) const { return _centerAtoms[index]; } 
	 
	bool isInRegion(const CxVector3 &vector); 
	 
private: 
	float _xmin; 
	float _xmax; 
	float _ymin; 
	float _ymax; 
	float _zmin; 
	float _zmax; 
	 
	unsigned char *_centerAtomTypes; 
	unsigned char *_centerAtomRealTypes; 
	unsigned char *_centerAtoms; 
}; 

 
bool gatherRegionAnalysis(); 
void processRegionAnalysis(CTimeStep *ts); 
void finalizeRegionAnalysis(); 
 
#endif 
