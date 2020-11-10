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

#ifndef IR_H
#define IR_H

#include "moltools.h"

class CTimeStep;

class CPowerObservation: public CObservation
{
public:
	CPowerObservation(bool global = false);
	~CPowerObservation();
	
	void initialize();
	void process(CTimeStep *ts);
	void finalize();
	
private:
	char *_name;
	CAtomGroup *_atoms;
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	int _zeroPadding;
	int _specSize;
	float _specResolution;
	bool _saveACF;
	bool _massWeighting;

	
	CxObArray _velocityCache;
	CxFloatArray _masses;
};

bool gatherPowerSpectrum();
bool initializePowerSpectrum();
void processPowerSpectrum(CTimeStep *ts);
void finalizePowerSpectrum();

class CIRObservation: public CObservation
{
public:
	CIRObservation(bool global = false);
	~CIRObservation();
	
	void initialize();
	void process();
	void finalize();
	
private:
	char *_name;
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	int _zeroPadding;
	int _specSize;
	float _specResolution;
	bool _saveACF;
	bool _includeCross;
	bool _finiteDifferenceCorrection;

	
	CxObArray _dipoleCache;
};

bool gatherIR();
bool initializeIR();
void processIR(CTimeStep *ts);
void finalizeIR();

bool gatherDipoleRestart();
bool initializeDipoleRestart();
void processDipoleRestart(CTimeStep *ts);
void finalizeDipoleRestart();

class CVCDObservation: public CObservation
{
public:
	CVCDObservation(bool global = false);
	~CVCDObservation();
	
	void initialize();
	void process();
	void finalize();
	
private:
	char *_name;
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	int _zeroPadding;
	int _specSize;
	float _specResolution;
	bool _saveMoments;
	bool _saveACF;
	bool _finiteDifferenceCorrection;
	
	int _smoothWidth;
	
	CxObArray _electricDipoleCache;
	CxObArray _magneticDipoleCache;
};

bool gatherVCD();
bool initializeVCD();
void processVCD(CTimeStep *ts);
void finalizeVCD();

bool gatherMagneticDipoleRestart();
bool initializeMagneticDipoleRestart();
void processMagneticDipoleRestart(CTimeStep *ts);
void finalizeMagneticDipoleRestart();

bool gatherSortWannier();
bool initializeSortWannier();
void processSortWannier(CTimeStep *ts);
void finalizeSortWannier();

#endif
