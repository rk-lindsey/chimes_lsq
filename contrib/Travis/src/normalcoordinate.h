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

#ifndef NORMALCOORDINATE_H
#define NORMALCOORDINATE_H

#include "moltools.h"
#include "xfloatarray.h"
#include "xintarray.h"
#include "xobject.h"
#include "xobarray.h"
#include "xvector3.h"

class CTimeStep;
class CxByteArray;
class CxVec3Array;

class CNormalCoordinateObservation;

class CReferenceStructure: public CxObject
{
public:
	CReferenceStructure(int showMol, const char *basename, bool calcIR);
	~CReferenceStructure();
	
	void initialize(const char *basename);
	void processCoordinates(CxVec3Array &coord, CxVector3 &dipole, int showMol);
	void calcVelocities(int showMol);
	float calcMinimumDistance(int showMol, bool useInternals, CxObArray &internals);
	void nextStep();
	float getMinimumDistance(int showMol, int step);
	void setGlobalProbability(int showMol, int step, float prob);
	void finalize(const char *basename);
	
private:
	CTimeStep *_refTimestep;
	CSingleMolecule *_singleMol;
	int _showMolCount;
	int _atomCount;
	int _permutationCount;
	int _numSteps;
	
	bool _writeTransformedTrajectories;
	FILE **_transformedTrajectoryFiles;
	
	CxObArray _permutations;
	CxVector3 _centroid;
	CxObArray _centroidCoords;
	CxFloatArray _weights;
	float _weightSum;
	char _filename[128];
	
	CxObArray _distanceTimedev;
	CxObArray _coordHistory;
	int _historyIndex;
	bool _calcVelocity;
	CxObArray _velocityCache;
	CxObArray _globalProb;
	
	bool _useInternals;
	CxObArray _internals;
	float _gaussianWidth;
	
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	float _specResolution;
	int _specSize;
	int _zeroPadding;
	
	float _convergenceThreshold;
	int _maxIterations;
	
	bool _calcIR;
	CxObArray _dipoleHistory;
	CxObArray _dipoleDerivativeCache;
	
	bool recognizeMolecule(CTimeStep *ts, int showMol, const char *basename);
	void recognizeMoleculeRecursion(unsigned int index, bool *used, int depth, unsigned int *stack, CTimeStep *ts, CxByteArray &baAtomIndex);
	bool recognizeMoleculeBondRange(CTimeStep *ts, int i1, int i2, float *distance, CxByteArray &baAtomIndex);
	void askPermutations(CxIntArray *actions);
	void createPermutationsRecursion(int start, CxIntArray *permutation, CxIntArray *permutationActions);
	
	double offDiagonalNorm(CxObArray &matrix);
	double findRotationAngle(CxObArray &matrix, int i, int j);
	void calcIntegrals(CxObArray &matrix, CxFloatArray *integrals, CxFloatArray *centers);
};

class CNormalCoordinateObservation: public CObservation
{
public:
	CNormalCoordinateObservation();
	~CNormalCoordinateObservation();
	
	void initialize();
	void process(CTimeStep *ts);
	void finalize();
	
private:
	CxObArray _refStructures;
	int _refCount;
	int _numSteps;
	float _gaussianWidth;
	int _correlationDepth;
	
	CxObArray _distanceTimedev;
	bool _useInternals;
	CxObArray _internals;
	
	bool _calcIR;
	
// 	bool _writeTransformedTrajectories;
// 	FILE **_transformedTrajectoryFiles;
// 	CxObArray _distanceTimedev;
// 	CxObArray _coordHistory;
// 	int _historyIndex;
// 	bool _calcVelocity;
// 	CxObArray _velocityCache;
};

bool gatherNormalCoordinate();
bool initializeNormalCoordinate();
void processNormalCoordinate(CTimeStep *ts);
void finalizeNormalCoordinate();

class CEckartReferenceStructure: public CxObject
{
public:
	CEckartReferenceStructure();
	~CEckartReferenceStructure();
	
	void process(CTimeStep *ts);
	const float *rotationMatrix() { return _rotationMatrix; }
	const CxVector3 &translationVector() { return _translationVector; }
	
private:
	CTimeStep *_refTimestep;
	CSingleMolecule *_singleMol;
	char _filename[128];
	int _singleMolIndex;
	CAtomGroup *_mapAtoms;
	CxVector3 _refCentroid;
	CxVec3Array _refCentroidCoord;
	CxFloatArray _weights;
	float _weightSum;
	
	float _rotationMatrix[9];
	CxVector3 _translationVector;
	
	bool recognizeMolecule(CTimeStep *ts);
	void recognizeMoleculeRecursion(unsigned int index, bool *used, int depth, unsigned int *stack, CTimeStep *ts, CxByteArray &baAtomIndex);
	bool recognizeMoleculeBondRange(CTimeStep *ts, int i1, int i2, float *distance, CxByteArray &baAtomIndex);
};

bool gatherEckartTransform();
bool initializeEckartTransform();
void processEckartTransform(CTimeStep *ts);
void finalizeEckartTransform();

#endif
