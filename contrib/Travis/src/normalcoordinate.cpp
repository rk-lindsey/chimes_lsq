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


#include "normalcoordinate.h"

#include "globalvar.h"
#include "linalg.h"
#include "maintools.h"
#include "moltools.h"
#include "timestep.h"
#include "tools.h"
#include "xbytearray.h"
#include "xobarray.h"
#include "xstring.h"

#include <stdio.h>


#define BUF_SIZE 4096


static CxObArray g_normalCoordinateObservations;


static void SVD_3x3(float *input_matrix, float *s_values, float *u_matrix, float *v_matrix) {
	float v[9];
	
	ComputeSVD_Flat(input_matrix, 3, 3, s_values, v);
	
	int sortIndex[3];
	int i, j;
	for(i = 0; i < 3; i++) {
		sortIndex[i] = i;
	}
	for(i = 0; i < 2; i++) {
		for(j = i + 1; j < 3; j++) {
			if(s_values[j] > s_values[i]) {
				float tf = s_values[i];
				s_values[i] = s_values[j];
				s_values[j] = tf;
				int ti = sortIndex[i];
				sortIndex[i] = sortIndex[j];
				sortIndex[j] = ti;
			}
		}
	}
	for(i = 0; i < 3; i++) {
		for(j = 0; j < 3; j++) {
			u_matrix[3*i + j] = input_matrix[3*i + sortIndex[j]];
			v_matrix[3*i + j] = v[3*i + sortIndex[j]];
		}
	}
}

CReferenceStructure::CReferenceStructure(int showMol, const char *basename, bool calcIR) {
	int i, j;
	CxString buf;

	_calcIR = calcIR;
	
	try { _singleMol = new CSingleMolecule(); } catch(...) { _singleMol = NULL; }
	if(_singleMol == NULL) NewException((double)sizeof(CSingleMolecule), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	try { _refTimestep = new CTimeStep(); } catch(...) { _refTimestep = NULL; }
	if(_refTimestep == NULL) NewException((double)sizeof(CTimeStep), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	while(true) {
//		char buf[BUF_SIZE];
		AskString_ND("    Filename for reference structure: ", &buf);
		FILE *refFile = fopen(buf, "r");
		if(refFile == NULL) {
			mprintf(RED, "Could not open reference structure file: %s\n", strerror(errno));
			continue;
		}
		if(!_refTimestep->ReadTimestep(refFile, true)) {
			mprintf(RED, "Error reading reference structure file\n");
			continue;
		}
		strncpy(_filename, buf, 128);
		_filename[127] = 0;
		mprintf("\n    Starting molecule recoginition...\n");
		if(!recognizeMolecule(_refTimestep, showMol, basename)) {
			mprintf(RED, "Molecule recoginition failed\n\n");
			continue;
		}
		break;
	}
	
	if(g_bAdvanced2) {
		_writeTransformedTrajectories = AskYesNo("    Save trajectories after transformation to reference frame for first observed molecule (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_writeTransformedTrajectories = false;
	}
	
	_showMolCount = ((CMolecule *)g_oaMolecules[showMol])->m_laSingleMolIndex.GetSize();
	
	if(_permutationCount > 1) {
		_useInternals = AskYesNo("    Use mass-weighted Cartesians for probabilities (n) or define internal coordinates (y)? [no] ", false);
		if(_useInternals) {
			mprintf("\n    At the moment, only \"simple\" dihedral angles can be defined.\n");
			do {
				CxIntArray *internal;
				try { internal = new CxIntArray(); } catch(...) { internal = NULL; }
				if(internal == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				internal->SetSize(4);
				mprintf("\n");
				for(i = 0; i < 4; i++) {
//					char buf[BUF_SIZE];
					unsigned char parse[3];
					do
						AskString_ND("      Enter the %d. atom (e.g. C7): ", &buf, i+1);
					while(!ParseAtom(buf, showMol, parse[0], parse[1], parse[2]));
					internal->GetAt(i) = -1;
					for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
						CMolAtom *ma = (CMolAtom *)_singleMol->m_oaMolAtoms[j];
						if((ma->m_iType == parse[0]) && (ma->m_iNumber == parse[2])) {
							internal->GetAt(i) = j;
							break;
						}
					}
					if(internal->GetAt(i) == -1) {
						mprintf(RED, "Weird error.\n");
						abort();
					}
				}
				_internals.Add(internal);
			} while(AskYesNo("\n    Add another dihedral angle (y/n)? [no] ", false));
			_gaussianWidth = AskFloat("\n    Gaussian width for probability distribution (deg)? [10.0] ", 10.0f);
		} else {
			_gaussianWidth = AskFloat("    Gaussian width for probability distribution (pm*sqrt(amu))? [1.0] ", 1.0f);
		}
		mprintf("\n");
	} else {
		_useInternals = false;
		_gaussianWidth = 1000.0f;
	}
	
	if(g_iTrajSteps != -1) {
		_correlationDepth = (int)(0.75 * g_iTrajSteps);
		if(_correlationDepth > 4096)
			_correlationDepth = 4096;
		if(g_fTimestepLength > 1.0f)
			_correlationDepth = 2048;
		if(g_fTimestepLength > 2.0f)
			_correlationDepth = 1024;
		_correlationDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the ACF (in time steps): [%d] ", _correlationDepth, _correlationDepth);
	} else {
		_correlationDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the ACF (in time steps): [256] ", 256);
	}
	int size = CalcFFTSize(_correlationDepth, false);
	if(_correlationDepth != size) {
		mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n", size, _correlationDepth);
		_correlationDepth = size;
	}
	
	if(g_bAdvanced2) {
		_windowFunction = AskRangeInteger("    Window function: cos^2(a*t) (1), exp(-t/a) (2), exp(-(t/a)^2) (3) [1] ", 1, 3, 1);
		if(_windowFunction == 1) {
			mprintf("    The parameter \"a\" is chosen according to the correlation depth.\n");
			_windowFunctionParameter = 0;
		} else if(_windowFunction == 2) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in time steps): [%d] ", _correlationDepth / 4, _correlationDepth / 4);
		} else if(_windowFunction == 3) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in time steps): [%d] ", _correlationDepth / 2, _correlationDepth / 2);
		} else {
			eprintf("This is impossible.\n");
			abort();
		}
	} else {
		_windowFunction = 1;
		_windowFunctionParameter = 0;
	}
	
	if(g_bAdvanced2) {
		_zeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ", _correlationDepth * 3, _correlationDepth * 3);
		size = CalcFFTSize(_correlationDepth + _zeroPadding, false);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
	} else {
		_zeroPadding = _correlationDepth * 3;
		size = CalcFFTSize(_correlationDepth + _zeroPadding, false);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
	}
	
	float possibleRange = 33356.41f / g_fTimestepLength / 2.0f;
	_specResolution = possibleRange / (_correlationDepth + _zeroPadding);
	mprintf("    This results in a spectral resolution of %.2f cm^-1.\n", _specResolution);
	mprintf("\n    A time step length of %.2f fs allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, possibleRange);
	float specLimit = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ", 0, possibleRange, (possibleRange < 5000.0f) ? possibleRange : 5000.0f, (possibleRange < 5000.0f) ? possibleRange : 5000.0f);
	_specSize = (int)(specLimit / _specResolution);
	mprintf("\n");
	
	_convergenceThreshold = AskFloat("    Relative convergence threshold for diagonalization: [1.0e-6] ", 1.0e-6f);
	_maxIterations = AskUnsignedInteger("    Maximum number of iterations for diagonalization: [50] ", 50);
	mprintf("\n");
}


CReferenceStructure::~CReferenceStructure() {
	delete _singleMol;
	delete _refTimestep;
	int i;
	for(i = 0; i < _centroidCoords.GetSize(); i++)
		delete (CxVec3Array *)_centroidCoords[i];
	
	for(i = 0; i < _distanceTimedev.GetSize(); i++)
		delete (CxFloatArray *)_distanceTimedev[i];
	for(i = 0; i < _coordHistory.GetSize(); i++)
		delete (CxVec3Array *)_coordHistory[i];
	for(i = 0; i < _velocityCache.GetSize(); i++)
		delete (CxVec3Array *)_velocityCache[i];
	
	for(i = 0; i < _dipoleHistory.GetSize(); i++)
		delete (CxVec3Array *)_dipoleHistory[i];
	for(i = 0; i < _dipoleDerivativeCache.GetSize(); i++)
		delete (CxVec3Array *)_dipoleDerivativeCache[i];
}


void CReferenceStructure::initialize(const char *basename) {
	_numSteps = 0;
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	mprintf("    Distance time development: Trying to allocate %s of memory...\n", FormatBytes((double)_showMolCount * n * _permutationCount * sizeof(float)));
	int i;
	for(i = 0; i < _showMolCount * _permutationCount; i++) {
		CxFloatArray *a;
		try { a = new CxFloatArray(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetMaxSize(n);
		a->SetGrow((int)(0.1 * n));
		_distanceTimedev.Add(a);
	}
	
	mprintf("    Coordinate history: Trying to allocate %s of memory...\n", FormatBytes((double)_showMolCount * _permutationCount * _atomCount * 3 * sizeof(CxVector3)));
	for(i = 0; i < _showMolCount * _permutationCount * _atomCount; i++) {
		CxVec3Array *a;
		try { a = new CxVec3Array(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetSize(3);
		_coordHistory.Add(a);
	}
	_historyIndex = 0;
	_calcVelocity = false;
	if(_calcIR) {
		mprintf("    Dipole history: Trying to allocate %s of memory...\n", FormatBytes((double)_showMolCount * _permutationCount * 3 * sizeof(CxVector3)));
		for(i = 0; i < _showMolCount * _permutationCount; i++) {
			CxVec3Array *a;
			try { a = new CxVec3Array(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			a->SetSize(3);
			_dipoleHistory.Add(a);
		}
	}
	
	mprintf("    Velocity cache: Trying to allocate %s of memory...\n", FormatBytes((double)_showMolCount * _permutationCount * _atomCount * n * sizeof(CxVector3)));
	for(i = 0; i < _showMolCount * _permutationCount * _atomCount; i++) {
		CxVec3Array *a;
		try { a = new CxVec3Array(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetMaxSize(n);
		a->SetGrow((int)(0.1 * n));
		_velocityCache.Add(a);
	}
	if(_calcIR) {
		mprintf("    Dipole derivative cache: Trying to allocate %s of memory...\n", FormatBytes((double)_showMolCount * _permutationCount * n * sizeof(CxVector3)));
		for(i = 0; i < _showMolCount * _permutationCount; i++) {
			CxVec3Array *a;
			try { a = new CxVec3Array(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			a->SetMaxSize(n);
			a->SetGrow((int)(0.1 * n));
			_dipoleDerivativeCache.Add(a);
		}
	}
	
	if(_writeTransformedTrajectories) {
		try { _transformedTrajectoryFiles = new FILE *[_permutationCount]; } catch(...) { _transformedTrajectoryFiles = NULL; }
		if(_transformedTrajectoryFiles == NULL) NewException((double)sizeof(FILE *) * _permutationCount, __FILE__, __LINE__, __PRETTY_FUNCTION__);
		for(i = 0; i < _permutationCount; i++) {
			char name[BUF_SIZE];
#ifdef TARGET_LINUX
			snprintf(name, BUF_SIZE, "%s_p%d_transformed.xyz", basename, i+1);
#else
			sprintf(name, "%s_p%d_transformed.xyz", basename, i+1);
#endif
			_transformedTrajectoryFiles[i] = OpenFileWrite(name, false);
		}
	}
}


void CReferenceStructure::processCoordinates(CxVec3Array &coord, CxVector3 &dipole, int showMol) {
	int i, j, k, l;
	for(i = 0; i < _permutationCount; i++) {
		if(_atomCount != coord.GetSize()) {
			mprintf(RED, "Wrong number of coordinates from trajectory.\n");
			abort();
		}
		
		CxVector3 centroid(0.0f, 0.0f, 0.0f);
		for(j = 0; j < coord.GetSize(); j++) {
			centroid += coord[j] * _weights[j];
		}
		centroid /= _weightSum;
		
		float covarianceMatrix[9] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
		for(j = 0; j < _atomCount; j++) {
			CxVector3 v1 = coord[j] - centroid;
			CxVector3 v2 = ((CxVec3Array *)_centroidCoords[i])->GetAt(j);
			for(k = 0; k < 3; k++) {
				for(l = 0; l < 3; l++) {
					covarianceMatrix[k*3+l] += v1[k] * v2[l] * _weights[j];
				}
			}
		}
		
//************** NEW ********************
// 		float sValues[3];
// 		float vMatrix[9];
// 		ComputeSVD_Flat(covarianceMatrix, 3, 3, sValues, vMatrix);
// 		
// 		float rotationMatrix[9];
// 		for(int j = 0; j < 3; j++) {
// 			for(int k = 0; k < 3; k++) {
// 				rotationMatrix[3*j+k] = 0.0f;
// 				for(int l = 0; l < 3; l++) {
// 					rotationMatrix[3*j+k] += vMatrix[3*j+l] * covarianceMatrix[3*k+l];
// 				}
// 			}
// 		}
// 		
// 		float determinant = rotationMatrix[0]*rotationMatrix[4]*rotationMatrix[8] + rotationMatrix[1]*rotationMatrix[5]*rotationMatrix[6] + rotationMatrix[2]*rotationMatrix[3]*rotationMatrix[7] - rotationMatrix[2]*rotationMatrix[4]*rotationMatrix[6] - rotationMatrix[1]*rotationMatrix[3]*rotationMatrix[8] - rotationMatrix[0]*rotationMatrix[5]*rotationMatrix[7];
// 		if(determinant < 0.0f) {
// // 			mprintf(RED, "det < 0");
// 			for(int j = 0; j < 3; j++) {
// 				for(int k = 0; k < 3; k++) {
// 					rotationMatrix[3*j+k] = 0.0f;
// 					for(int l = 0; l < 3; l++) {
// 						rotationMatrix[3*j+k] += vMatrix[3*j+l] * covarianceMatrix[3*k+l] * (l == 0 ? -1.0f : 1.0f);
// 					}
// 				}
// 			}
// 		}
// 		
// 		mprintf(GREEN, "%f", rotationMatrix[0]*rotationMatrix[4]*rotationMatrix[8] + rotationMatrix[1]*rotationMatrix[5]*rotationMatrix[6] + rotationMatrix[2]*rotationMatrix[3]*rotationMatrix[7] - rotationMatrix[2]*rotationMatrix[4]*rotationMatrix[6] - rotationMatrix[1]*rotationMatrix[3]*rotationMatrix[8] - rotationMatrix[0]*rotationMatrix[5]*rotationMatrix[7]);
//***************************************
//************* NEW 2 *******************
		float sValues[3];
		float uMatrix[9];
		float vMatrix[9];
		SVD_3x3(covarianceMatrix, sValues, uMatrix, vMatrix);
		
		float rotationMatrix[9];
		for(j = 0; j < 3; j++) {
			for(k = 0; k < 3; k++) {
				rotationMatrix[3*j+k] = 0.0f;
				for(l = 0; l < 3; l++) {
					rotationMatrix[3*j+k] += vMatrix[3*j+l] * uMatrix[3*k+l];
				}
			}
		}
		
		float determinant = rotationMatrix[0]*rotationMatrix[4]*rotationMatrix[8] + rotationMatrix[1]*rotationMatrix[5]*rotationMatrix[6] + rotationMatrix[2]*rotationMatrix[3]*rotationMatrix[7] - rotationMatrix[2]*rotationMatrix[4]*rotationMatrix[6] - rotationMatrix[1]*rotationMatrix[3]*rotationMatrix[8] - rotationMatrix[0]*rotationMatrix[5]*rotationMatrix[7];
		if(determinant < 0.0f) {
			for(j = 0; j < 3; j++) {
				for(k = 0; k < 3; k++) {
					rotationMatrix[3*j+k] = 0.0f;
					for(l = 0; l < 3; l++) {
						rotationMatrix[3*j+k] += vMatrix[3*j+l] * uMatrix[3*k+l] * (l == 2 ? -1.0f : 1.0f);
					}
				}
			}
		}
//***************************************
		
// 		float sValues[3];
// 		float uMatrix[9];
// 		float vtMatrix[9];
// 		singularValueDecompose_3x3(covarianceMatrix, sValues, uMatrix, vtMatrix);
// 		
// 		if(vtMatrix[0]*vtMatrix[4]*vtMatrix[8] + vtMatrix[1]*vtMatrix[5]*vtMatrix[6] + vtMatrix[2]*vtMatrix[3]*vtMatrix[7] - vtMatrix[2]*vtMatrix[4]*vtMatrix[6] - vtMatrix[1]*vtMatrix[3]*vtMatrix[8] - vtMatrix[0]*vtMatrix[5]*vtMatrix[7] < 0) {
// 			for(int j = 0; j < 3; j++) {
// 				vtMatrix[6+j] *= -1.0f;
// 			}
// 		}
// 		
// 		float rotationMatrix[9];
// 		for(int j = 0; j < 3; j++) {
// 			for(int k = 0; k < 3; k++) {
// 				rotationMatrix[3*k+j] = 0.0f;
// 				for(int l = 0; l < 3; l++) {
// 					rotationMatrix[3*k+j] += uMatrix[3*j+l] * vtMatrix[3*l+k];
// 				}
// 			}
// 		}
// 		float det = rotationMatrix[0]*rotationMatrix[4]*rotationMatrix[8] + rotationMatrix[1]*rotationMatrix[5]*rotationMatrix[6] + rotationMatrix[2]*rotationMatrix[3]*rotationMatrix[7] - rotationMatrix[2]*rotationMatrix[4]*rotationMatrix[6] - rotationMatrix[1]*rotationMatrix[3]*rotationMatrix[8] - rotationMatrix[0]*rotationMatrix[5]*rotationMatrix[7];
// 		if(det < 0) {
// 			mprintf(RED, "det < 0\n");
// 		}
		
		CxVector3 translation;
		for(j = 0; j < 3; j++) {
			translation[j] = _centroid[j];
			for(k = 0; k < 3; k++) {
				translation[j] -= rotationMatrix[3*j+k] * centroid[k];
			}
		}
		
		float dist = 0.0f;
		for(j = 0; j < _atomCount; j++) {
			CxVector3 newCoord = translation;
			for(k = 0; k < 3; k++) {
				for(l = 0; l < 3; l++) {
					newCoord[k] += rotationMatrix[3*k+l] * coord[j][l];
				}
			}
			((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(j)])->GetAt(_historyIndex) = newCoord;
// 			((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetPosition(j)])->GetAt(_historyIndex) = newCoord;
// 			((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + j])->GetAt(_historyIndex) = newCoord;
			CxVector3 distVec = newCoord - ((CxVec3Array *)_centroidCoords[i])->GetAt(j) - _centroid;
			dist += distVec.GetLengthSqr() * _weights[j];
		}
		if(_useInternals) {
			float idist = 0.0f;
			for(j = 0; j < _internals.GetSize(); j++) {
				float d1, d2;
				CxVector3 vec1, vec2, vec3;
				vec1 = FoldVector(((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(0)) - ((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(1)));
				vec2 = FoldVector(((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(3)) - ((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(2)));
				vec3 = FoldVector(((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(2)) - ((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(1)));
				d1 = Dihedral(vec1, vec2, vec3, false);
				vec1 = FoldVector(((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(0))])->GetAt(_historyIndex) - ((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(1))])->GetAt(_historyIndex));
				vec2 = FoldVector(((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(3))])->GetAt(_historyIndex) - ((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(2))])->GetAt(_historyIndex));
				vec3 = FoldVector(((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(2))])->GetAt(_historyIndex) - ((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)_internals[j])->GetAt(1))])->GetAt(_historyIndex));
				d2 = Dihedral(vec1, vec2, vec3, false);
				float diff = fabsf(d1 - d2);
				if(diff > 180.0f)
					diff = 360.0f - diff;
				idist += diff * diff;
			}
			((CxFloatArray *)_distanceTimedev[showMol*_permutationCount + i])->Add(sqrtf(idist));
		} else {
			((CxFloatArray *)_distanceTimedev[showMol*_permutationCount + i])->Add(sqrtf(dist) / sqrtf(_weightSum));
		}
		
		if(_calcIR) {
			CxVector3 newDipole(0.0f, 0.0f, 0.0f);
			for(j = 0; j < 3; j++) {
				for(k = 0; k < 3; k++) {
					newDipole[j] += rotationMatrix[3*j+k] * dipole[k];
				}
			}
			((CxVec3Array *)_dipoleHistory[showMol*_permutationCount + i])->GetAt(_historyIndex) = newDipole;
		}
		
		if(_writeTransformedTrajectories && (showMol == 0)) {
			fprintf(_transformedTrajectoryFiles[i], "%4d\n\n", _atomCount);
			for(j = 0; j < _atomCount; j++) {
				fprintf(_transformedTrajectoryFiles[i], "%4s %12.6f %12.6f %12.6f\n", ((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType]])->m_sName, ((CxVec3Array *)_coordHistory[i*_atomCount + j])->GetAt(_historyIndex)[0] / 100.0f, ((CxVec3Array *)_coordHistory[i*_atomCount + j])->GetAt(_historyIndex)[1] / 100.0f, ((CxVec3Array *)_coordHistory[i*_atomCount + j])->GetAt(_historyIndex)[2] / 100.0f);
// 				fprintf(_transformedTrajectoryFiles[i], "%4s %12.6f %12.6f %12.6f\n", ((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType]])->m_sName, ((CxVec3Array *)_coordHistory[i*_atomCount + ((CxIntArray *)_permutations[i])->GetPosition(j)])->GetAt(_historyIndex)[0] / 100.0f, ((CxVec3Array *)_coordHistory[i*_atomCount + ((CxIntArray *)_permutations[i])->GetPosition(j)])->GetAt(_historyIndex)[1] / 100.0f, ((CxVec3Array *)_coordHistory[i*_atomCount + ((CxIntArray *)_permutations[i])->GetPosition(j)])->GetAt(_historyIndex)[2] / 100.0f);
			}
		}
	}
}


void CReferenceStructure::calcVelocities(int showMol) {
	int n = (_historyIndex + 1) % 3;
	int i, j;
	for(i = 0; i < _permutationCount; i++) {
		for(j = 0; j < _atomCount; j++) {
			CxVector3 velVec;
			velVec = (((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + j])->GetAt(_historyIndex) - ((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + j])->GetAt(n)) / 2.0f / g_fTimestepLength * 1000.0f;
			((CxVec3Array *)_velocityCache[showMol*_permutationCount*_atomCount + i*_atomCount + j])->Add(velVec);
		}
		
		if(_calcIR) {
			CxVector3 derVec;
			derVec = (((CxVec3Array *)_dipoleHistory[showMol*_permutationCount + i])->GetAt(_historyIndex) - ((CxVec3Array *)_dipoleHistory[showMol*_permutationCount + i])->GetAt(n)) / 2.0f / g_fTimestepLength;
			((CxVec3Array *)_dipoleDerivativeCache[showMol*_permutationCount + i])->Add(derVec);
		}
	}
}


float CReferenceStructure::calcMinimumDistance(int showMol, bool useInternals, CxObArray &internals) {
	float minDist = 0.0f;
	int i, j;
	for(i = 0; i < _permutationCount; i++) {
		if(useInternals) {
			float idist = 0.0f;
			for(j = 0; j < internals.GetSize(); j++) {
				float d1, d2;
				CxVector3 vec1, vec2, vec3;
				vec1 = FoldVector(((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)internals[j])->GetAt(0)) - ((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)internals[j])->GetAt(1)));
				vec2 = FoldVector(((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)internals[j])->GetAt(3)) - ((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)internals[j])->GetAt(2)));
				vec3 = FoldVector(((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)internals[j])->GetAt(2)) - ((CxVec3Array *)_centroidCoords[i])->GetAt(((CxIntArray *)internals[j])->GetAt(1)));
				d1 = Dihedral(vec1, vec2, vec3, false);
				vec1 = FoldVector(((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)internals[j])->GetAt(0))])->GetAt(_historyIndex) - ((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)internals[j])->GetAt(1))])->GetAt(_historyIndex));
				vec2 = FoldVector(((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)internals[j])->GetAt(3))])->GetAt(_historyIndex) - ((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)internals[j])->GetAt(2))])->GetAt(_historyIndex));
				vec3 = FoldVector(((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)internals[j])->GetAt(2))])->GetAt(_historyIndex) - ((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(((CxIntArray *)internals[j])->GetAt(1))])->GetAt(_historyIndex));
				d2 = Dihedral(vec1, vec2, vec3, false);
				float diff = fabsf(d1 - d2);
				if(diff > 180.0f)
					diff = 360.0f - diff;
				idist += diff * diff;
			}
			idist = sqrtf(idist);
			if(i == 0) {
				minDist = idist;
			} else {
				if(idist < minDist) {
					minDist = idist;
				}
			}
		} else {
			float dist = 0.0f;
			for(j = 0; j < _atomCount; j++) {
				CxVector3 distVec = ((CxVec3Array *)_centroidCoords[i])->GetAt(j) - ((CxVec3Array *)_coordHistory[showMol*_permutationCount*_atomCount + i*_atomCount + ((CxIntArray *)_permutations[i])->GetAt(j)])->GetAt(_historyIndex);
				dist += distVec.GetLengthSqr() * _weights[j];
			}
			dist = sqrtf(dist) / sqrtf(_weightSum);
			if(i == 0) {
				minDist = dist;
			} else {
				if(dist < minDist) {
					minDist = dist;
				}
			}
		}
	}
	return minDist;
}


void CReferenceStructure::nextStep() {
	_historyIndex = (_historyIndex + 1) % 3;
	_numSteps++;
}


float CReferenceStructure::getMinimumDistance(int showMol, int step) {
	if(step > _numSteps) {
		mprintf(RED, "Illegal step number.\n");
		abort();
	}
	float minDist = ((CxFloatArray *)_distanceTimedev[showMol*_permutationCount])->GetAt(step);
	int i;
	for(i = 1; i < _permutationCount; i++) {
		if(((CxFloatArray *)_distanceTimedev[showMol*_permutationCount + i])->GetAt(step) < minDist)
			minDist = ((CxFloatArray *)_distanceTimedev[showMol*_permutationCount + i])->GetAt(step);
	}
	return minDist;
}


void CReferenceStructure::setGlobalProbability(int showMol, int step, float prob) {
	if(_globalProb.GetSize() < _showMolCount) {
		int i;
		for(i = 0; i < _showMolCount; i++) {
			CxFloatArray *a;
			try { a = new CxFloatArray(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			a->SetSize(_numSteps);
			_globalProb.Add(a);
		}
	}
	((CxFloatArray *)_globalProb[showMol])->GetAt(step) = prob;
}


void CReferenceStructure::finalize(const char *basename) {
	int i, j, k, l, m, n, o;
	if(_writeTransformedTrajectories) {
		for(i = 0; i < _permutationCount; i++) {
			fclose(_transformedTrajectoryFiles[i]);
		}
		delete[] _transformedTrajectoryFiles;
	}
	
	mprintf("    Saving distance time development...\n");
	char name[BUF_SIZE];
#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "%s_dist_timedev.csv", basename);
#else
	sprintf(name, "%s_dist_timedev.csv", basename);
#endif
	FILE *distanceFile = OpenFileWrite(name, false);
	fprintf(distanceFile, "#Time (fs);");
	for(i = 0; i < _showMolCount; i++) {
		for(j = 0; j < _permutationCount; j++) {
			fprintf(distanceFile, " M%d P%d;", i+1, j+1);
		}
	}
	fprintf(distanceFile, "\n");
	for(i = 0; i < _numSteps; i++) {
		fprintf(distanceFile, "%.2f;", i * g_fTimestepLength * g_iStride);
		for(j = 0; j < _showMolCount; j++) {
			for(k = 0; k < _permutationCount; k++) {
				fprintf(distanceFile, " %.8G;", ((CxFloatArray *)_distanceTimedev[j*_permutationCount + k])->GetAt(i));
			}
		}
		fprintf(distanceFile, "\n");
	}
	fclose(distanceFile);
	
	mprintf("    Saving probability time development...\n");
#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "%s_prob_timedev.csv", basename);
#else
	sprintf(name, "%s_prob_timedev.csv", basename);
#endif
	FILE *probFile = OpenFileWrite(name, false);
	fprintf(probFile, "#Time (fs);");
	for(i = 0; i < _showMolCount; i++) {
		for(j = 0; j < _permutationCount; j++) {
			fprintf(probFile, " M%d P%d;", i+1, j+1);
		}
	}
	fprintf(probFile, "\n");
	CxDoubleArray probSum;
	probSum.SetSize(_showMolCount*_numSteps);
	for(i = 0; i < _numSteps; i++) {
		fprintf(probFile, "%.2f;", i * g_fTimestepLength * g_iStride);
		for(j = 0; j < _showMolCount; j++) {
			CxDoubleArray probs;
			probs.SetSize(_permutationCount);
			probSum[j*_numSteps + i] = 0.0;
			if(_permutationCount == 1) {
				probs[0] = 1.0;
				probSum[j*_numSteps + i] = 1.0;
				fprintf(probFile, " %.8G;", probs[0] / probSum[j*_numSteps + i]);
			} else {
				for(k = 0; k < _permutationCount; k++) {
					probs[k] = exp(-(double)((CxFloatArray *)_distanceTimedev[j*_permutationCount + k])->GetAt(i) * (double)((CxFloatArray *)_distanceTimedev[j*_permutationCount + k])->GetAt(i) / 2.0 / (double)_gaussianWidth / (double)_gaussianWidth);
					probSum[j*_numSteps + i] += probs[k];
				}
				for(k = 0; k < _permutationCount; k++) {
					if(isnan(probs[k] / probSum[j*_numSteps + i])) {
						mprintf(RED, "Invalid probability for permutation %d of molecule %d in step %d.\nMaybe it helps to increase the distribution width.\n", k+1, j+1, i+1);
						abort();
					}
					fprintf(probFile, " %.8G;", probs[k] / probSum[j*_numSteps + i]);
				}
			}
		}
		fprintf(probFile, "\n");
	}
	fclose(probFile);
	
	mprintf("    Calculating cross-correlation matrix...\n");
	CxObArray ccMatrix;
	for(i = 0; i < 3 * _atomCount * (3 * _atomCount + 1) / 2; i++) {
		CxFloatArray *a;
		try { a = new CxFloatArray(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetSize(_correlationDepth);
		for(j = 0; j < _correlationDepth; j++) {
			a->GetAt(j) = 0.0f;
		}
		ccMatrix.Add(a);
	}
	CCrossCorrelation *cc;
	try { cc = new CCrossCorrelation(); } catch(...) { cc = NULL; }
	if(cc == NULL) NewException((double)sizeof(CCrossCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	cc->Init(_numSteps-2, _correlationDepth, g_bACFFFT);
	
	CxFloatArray *temp;
	try { temp = new CxFloatArray(); } catch(...) { temp = NULL; }
	if(temp == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp->SetSize(_numSteps-2);
	
	CxFloatArray *temp2;
	try { temp2 = new CxFloatArray(); } catch(...) { temp2 = NULL; }
	if(temp2 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp2->SetSize(_numSteps-2);
	
	CxFloatArray *temp3;
	try { temp3 = new CxFloatArray(); } catch(...) { temp3 = NULL; }
	if(temp3 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp3->SetSize(_correlationDepth);
	
	mprintf(WHITE, "     [");
	float step = (float)(_showMolCount * _permutationCount * 3*_atomCount * (3*_atomCount+1)/2) / 20.0f;
	int c = 0;
	
	for(i = 0; i < _showMolCount; i++) {
// 	for(int i = 1; i < 2; i++) {
		for(j = 0; j < _permutationCount; j++) {
// 		for(int j = 0; j < 1; j++) {
			for(k = 0; k < _atomCount; k++) {
				for(l = 0; l < 3; l++) {
					for(m = k; m < _atomCount; m++) {
						for(n = ((m == k) ? l : 0); n < 3; n++) {
							if(fmodf((float)c++, step) < 1.0f)
								mprintf(WHITE, "#");
							for(o = 0; o < _numSteps - 2; o++) {
// 								temp->GetAt(o) = expf(-((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) * ((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) / 2.0f / _gaussianWidth / _gaussianWidth) / probSum[i] * ((CxVec3Array *)_velocityCache[i*_permutationCount*_atomCount + j*_atomCount + k])->GetAt(o)[l] * sqrtf(_weights[k]);
// 								temp->GetAt(o) = ((CxFloatArray *)_globalProb[i])->GetAt(o) * expf(-((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) * ((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) / 2.0f / _gaussianWidth / _gaussianWidth) / probSum[i] * ((CxVec3Array *)_velocityCache[i*_permutationCount*_atomCount + j*_atomCount + k])->GetAt(o)[l] * sqrtf(_weights[k]);
// 								mprintf(RED, "%f %f %f\n", exp(-(double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) * (double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) / 2.0 / (double)_gaussianWidth / (double)_gaussianWidth), probSum[i*_numSteps + o+1], (exp(-(double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) * (double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) / 2.0 / (double)_gaussianWidth / (double)_gaussianWidth) / probSum[i*_numSteps + o+1]));
								if(_permutationCount == 1) {
									temp->GetAt(o) = ((CxFloatArray *)_globalProb[i])->GetAt(o+1) * ((CxVec3Array *)_velocityCache[i*_permutationCount*_atomCount + j*_atomCount + k])->GetAt(o)[l] * sqrtf(_weights[k]);
								} else {
									temp->GetAt(o) = ((CxFloatArray *)_globalProb[i])->GetAt(o+1) * (float)(exp(-(double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) * (double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(o+1) / 2.0 / (double)_gaussianWidth / (double)_gaussianWidth) / probSum[i*_numSteps + o+1]) * ((CxVec3Array *)_velocityCache[i*_permutationCount*_atomCount + j*_atomCount + k])->GetAt(o)[l] * sqrtf(_weights[k]);
								}
								temp2->GetAt(o) = ((CxVec3Array *)_velocityCache[i*_permutationCount*_atomCount + j*_atomCount + m])->GetAt(o)[n] * sqrtf(_weights[m]);
							}
							cc->CrossCorrelate(temp, temp2, temp3);
							for(o = 0; o < _correlationDepth; o++) {
								((CxFloatArray *)ccMatrix[(3*_atomCount-1)*(3*k+l) - (3*k+l)*(3*k+l-1)/2 + (3*m+n)])->GetAt(o) += temp3->GetAt(o);
							}
						}
					}
				}
			}
		}
	}
	mprintf(WHITE, "]\n");
	for(k = 0; k < _atomCount; k++) {
		for(l = 0; l < 3; l++) {
			for(m = k; m < _atomCount; m++) {
				for(n = ((m == k) ? l : 0); n < 3; n++) {
					for(o = 0; o < _correlationDepth; o++) {
						((CxFloatArray *)ccMatrix[(3*_atomCount-1)*(3*k+l) - (3*k+l)*(3*k+l-1)/2 + (3*m+n)])->GetAt(o) /= (float)_showMolCount;
					}
				}
			}
		}
	}
	
// 	for(int i = 0; i < 3 * _atomCount; i++) {
// 		for(int j = i; j < 3 * _atomCount; j++) {
// 			char name[BUF_SIZE];
// 			snprintf(name, BUF_SIZE, "ccf_%d_%d.dat", i, j);
// 			FILE *ccfFile = fopen(name, "w");
// 			for(int k = 0; k < _correlationDepth; k++) {
// 				fprintf(ccfFile, "%g %g %g\n", g_fTimestepLength * k, ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k), ((CxVec3Array *)_velocityCache[i])->GetAt(k)[0]);
// 			}
// 			fclose(ccfFile);
// 		}
// 	}
	
	CFFT *fft;
	try { fft = new CFFT(); } catch(...) { fft = NULL; }
	if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	fft->PrepareFFT_C2C(2 * (((CxFloatArray *)ccMatrix[0])->GetSize() + _zeroPadding));
	
	mprintf("    Fourier transforming cross-correlation matrix...\n");
	mprintf(WHITE, "     [");
	step = (float)(3*_atomCount * (3*_atomCount+1)/2) / 20.0f;
	c = 0;
	for(i = 0; i < 3 * _atomCount; i++) {
		for(j = i; j < 3 * _atomCount; j++) {
			if(fmodf((float)c++, step) < 1.0f)
				mprintf(WHITE, "#");
			
			temp->CopyFrom((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j]);
			
			if(_windowFunction == 1) {
				for(k = 0; k < temp->GetSize(); k++) {
					temp->GetAt(k) *= powf(cosf((float)k / (temp->GetSize() - 1) / 2.0f * Pi), 2.0f);
				}
			} else if(_windowFunction == 2) {
				for(k = 0; k < temp->GetSize(); k++) {
					temp->GetAt(k) *= expf(-(float)k / _windowFunctionParameter);
				}
			} else if(_windowFunction == 3) {
				for(k = 0; k < temp->GetSize(); k++) {
					temp->GetAt(k) *= expf(-(float)k * k / _windowFunctionParameter / _windowFunctionParameter);
				}
			} else if(_windowFunction != 0) {
				eprintf("Unknown window function.\n");
				abort();
			}
			
			if(_zeroPadding > 0) {
				for(k = 0; k < _zeroPadding; k++) {
					temp->Add(0.0f);
				}
			}
			
			int oldSize = temp->GetSize();
			temp->SetSize(2 * oldSize);
			for(k = 1; k < oldSize; k++) {
				temp->GetAt(oldSize + k) = temp->GetAt(oldSize - k);
			}
			temp->GetAt(oldSize) = 0.0f;
			
// 			fft->PrepareFFT_C2C(temp->GetSize());
			for(k = 0; k < temp->GetSize(); k++) {
				fft->m_pInput[2*k] = temp->GetAt(k);
				fft->m_pInput[2*k+1] = 0.0f;
			}
			fft->DoFFT();
			((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->SetSize(_specSize);
			for(k = 0; k < _specSize; k++) {
				((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k) = 7.211349e-9f * fft->m_pOutput[2*k] * g_fTimestepLength; // Spectra in K*cm
// 				((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k) = sqrtf(fft->m_pOutput[2*k]*fft->m_pOutput[2*k]+fft->m_pOutput[2*k+1]*fft->m_pOutput[2*k+1]);
			}
		}
	}
	mprintf(WHITE, "]\n");
	
	delete fft;
	
	mprintf("    Saving power spectrum...\n");
#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "%s_power.csv", basename);
#else
	sprintf(name, "%s_power.csv", basename);
#endif
	FILE *powerFile = OpenFileWrite(name, false);
	fprintf(powerFile, "#Wavenumber (cm^-1); Spectrum (K*cm); Integral (K)\n");
	double integral = 0.0;
	for(i = 0; i < _specSize; i++) {
		float value = 0.0f;
		for(j = 0; j < 3*_atomCount; j++) {
			value += ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + j])->GetAt(i);
		}
		integral += (double)value * _specResolution;
		fprintf(powerFile, "%.2f; %.8G; %.14G\n", _specResolution*i, value, integral);
	}
	fclose(powerFile);
	
// 	mprintf("    Saving inital spectra...\n");
// 	snprintf(name, BUF_SIZE, "%s_speci.dat", basename);
// 	FILE *speciFile = OpenFileWrite(name, false);
// 	for(int i = 0; i < _specSize; i++) {
// 		fprintf(speciFile, "%10.4f", _specResolution * i);
// 		for(int j = 0; j < 3*_atomCount; j++) {
// 			fprintf(speciFile, " %20.8f", ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + j])->GetAt(i));
// 		}
// 		fprintf(speciFile, "\n");
// 	}
// 	fclose(speciFile);
// 
	CxFloatArray *transMatrix;
	try { transMatrix = new CxFloatArray(); } catch(...) { transMatrix = NULL; }
	if(transMatrix == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	transMatrix->SetMaxSize(9 * _atomCount * _atomCount);
	for(i = 0; i < 3 * _atomCount; i++) {
		for(j = 0; j < 3 * _atomCount; j++) {
			transMatrix->Add((i == j) ? 1.0f : 0.0f);
		}
	}
	
// 	for(int i = 0; i < 3 * _atomCount; i++) {
// 		for(int j = i; j < 3 * _atomCount; j++) {
// 			char name[BUF_SIZE];
// 			snprintf(name, BUF_SIZE, "speca_%d_%d.dat", i, j);
// 			FILE *specFile = fopen(name, "w");
// 			for(int k = 0; k < _specSize; k++) {
// 				fprintf(specFile, "%g %g\n", _specResolution * k, ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k));
// 			}
// 			fclose(specFile);
// 		}
// 	}
// 	FILE *testFile = fopen("test.out", "w");
// 	for(int i = 0; i < _specSize; i++) {
// 		fprintf(testFile, "%f", _specResolution * i);
// 		for(int j = 0; j < 3 * _atomCount; j++) {
// 			for(int k = j; k < 3 * _atomCount; k++) {
// 				fprintf(testFile, " %f", ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + k])->GetAt(i));
// 			}
// 		}
// 		fprintf(testFile, "\n");
// 	}
// 	fclose(testFile);
	
	mprintf("    Minimizing cross-correlations...\n\n");
	mprintf("     ---------------------------------------------------\n");
	mprintf("     Iteration    Off-diagonal norm               Change\n");
	mprintf("     ---------------------------------------------------\n");
	double norm = offDiagonalNorm(ccMatrix);
	double change = norm;
	mprintf("     %9d %20g %20g\n", 1, norm, change);
	int count = 0;
	temp2->SetSize(3*_atomCount);
	while((fabs(change) > _convergenceThreshold * norm) && (count < _maxIterations)) {
		count++;
		for(i = 0; i < 3*_atomCount; i++) {
			for(j = i+1; j < 3*_atomCount; j++) {
// 				mprintf(RED, "i: %d, j: %d\n", i, j);
				double t = findRotationAngle(ccMatrix, i, j);
				double c = 1.0 / sqrt(t * t + 1.0);
				double s = t * c;
// 				mprintf(RED, "(%d, %d) t: %f\n", i, j, t);
				for(k = 0; k < _specSize; k++) {
// 					mprintf(RED, "k: %d\n", k);
					float temp[3];
					for(l = 0; l < j; l++) {
						if(l < i) {
							temp[0] = c * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + i])->GetAt(k) - s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + j])->GetAt(k);
							temp[1] = s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + i])->GetAt(k) + c * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + j])->GetAt(k);
							((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + i])->GetAt(k) = temp[0];
							((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + j])->GetAt(k) = temp[1];
						} else if(l > i) {
							temp2->GetAt(l) = ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + j])->GetAt(k);
							temp[1] = s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + l])->GetAt(k) + c * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + j])->GetAt(k);
							((CxFloatArray *)ccMatrix[(3*_atomCount-1)*l - l*(l-1)/2 + j])->GetAt(k) = temp[1];
						}
					}
					for(l = i+1; l < 3*_atomCount; l++) {
						if(l < j) {
							temp[0] = c * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + l])->GetAt(k) - s * temp2->GetAt(l);
							((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + l])->GetAt(k) = temp[0];
						} else if(l > j) {
							temp[0] = c * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + l])->GetAt(k) - s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + l])->GetAt(k);
							temp[1] = s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + l])->GetAt(k) + c * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + l])->GetAt(k);
							((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + l])->GetAt(k) = temp[0];
							((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + l])->GetAt(k) = temp[1];
						}
					}
					temp[0] = c * c * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k) + s * s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + j])->GetAt(k) - 2 * c * s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k);
					temp[1] = (c * c - s * s) * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k) + c * s * (((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k) - ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + j])->GetAt(k));
					temp[2] = s * s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k) + c * c * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + j])->GetAt(k) + 2 * c * s * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k);
					((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k) = temp[0];
					((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k) = temp[1];
					((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + j])->GetAt(k) = temp[2];
				}
				float temp[2];
				for(l = 0; l < 3*_atomCount; l++) {
					temp[0] = c * transMatrix->GetAt(i * 3*_atomCount + l) - s * transMatrix->GetAt(j * 3*_atomCount + l);
					temp[1] = s * transMatrix->GetAt(i * 3*_atomCount + l) + c * transMatrix->GetAt(j * 3*_atomCount + l);
					transMatrix->GetAt(i * 3*_atomCount + l) = temp[0];
					transMatrix->GetAt(j * 3*_atomCount + l) = temp[1];
				}
			}
		}
		double newNorm = offDiagonalNorm(ccMatrix);
		change = newNorm - norm;
		norm = newNorm;
		mprintf("     %9d %20g %20g\n", count+1, norm, change);
	}
	mprintf("     ---------------------------------------------------\n");
	if(count == _maxIterations) {
		eprintf("\nMinimization did not converge.\n\n");
		mprintf("Maybe it helps to increase the number of iterations or the convergence threshold.\n");
		mprintf("Will write unconverged results anyway.\n\n"); 
	} else
		mprintf("     Convergence reached after %d iterations.\n\n", count+1);
	
	delete temp;
	delete temp2;
	delete temp3;
	
// 	for(int i = 0; i < 3 * _atomCount; i++) {
// 		for(int j = i; j < 3 * _atomCount; j++) {
// 			char name[BUF_SIZE];
// 			snprintf(name, BUF_SIZE, "specb_%d_%d.dat", i, j);
// 			FILE *specFile = fopen(name, "w");
// 			for(int k = 0; k < _specSize; k++) {
// 				fprintf(specFile, "%g %g\n", _specResolution * k, ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k));
// 			}
// 			fclose(specFile);
// 		}
// 	}
// 	FILE *test2File = fopen("test2.out", "w");
// 	for(int i = 0; i < _specSize; i++) {
// 		fprintf(test2File, "%f", _specResolution * i);
// 		for(int j = 0; j < 3 * _atomCount; j++) {
// 			for(int k = j; k < 3 * _atomCount; k++) {
// 				fprintf(test2File, " %f", ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + k])->GetAt(i));
// 			}
// 		}
// 		fprintf(test2File, "\n");
// 	}
// 	fclose(test2File);
	
	mprintf("    Calculating integrals...\n");
	CxFloatArray *integrals;
	try { integrals = new CxFloatArray(); } catch(...) { integrals = NULL; }
	if(integrals == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	CxFloatArray *centers;
	try { centers = new CxFloatArray(); } catch(...) { centers = NULL; }
	if(centers == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	calcIntegrals(ccMatrix, integrals, centers);
	
	mprintf("    Sorting normal coordinates...\n\n");
	CxIntArray *sortIndex;
	try { sortIndex = new CxIntArray(); } catch(...) { sortIndex = NULL; }
	if(sortIndex == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	for(i = 0; i < 6; i++) {
		int minIndex = 0;
		while(sortIndex->Contains(minIndex) && (minIndex < 3*_atomCount - 1))
			minIndex++;
		float minIntegral = integrals->GetAt(minIndex);
		for(j = minIndex+1; j < 3*_atomCount; j++) {
			if(sortIndex->Contains(j))
				continue;
			if(integrals->GetAt(j) < minIntegral) {
				minIndex = j;
				minIntegral = integrals->GetAt(minIndex);
			}
		}
		sortIndex->Add(minIndex);
	}
	for(i = 6; i < 3*_atomCount; i++) {
		int minIndex = 0;
		while(sortIndex->Contains(minIndex) && (minIndex < 3*_atomCount - 1))
			minIndex++;
		float minCenter = centers->GetAt(minIndex);
		for(j = minIndex+1; j < 3*_atomCount; j++) {
			if(sortIndex->Contains(j))
				continue;
			if(centers->GetAt(j) < minCenter) {
				minIndex = j;
				minCenter = centers->GetAt(minIndex);
			}
		}
		sortIndex->Add(minIndex);
	}
	mprintf("     ----------------------------------------------\n");
	mprintf("     Mode         Integral (K)       Center (cm^-1)\n");
	mprintf("     ----------------------------------------------\n");
	for(i = 0; i < 3*_atomCount; i++) {
		mprintf("     %4d %20g %20.2f\n", i+1, integrals->GetAt(sortIndex->GetAt(i)), centers->GetAt(sortIndex->GetAt(i)));
		if(i == 5)
			mprintf("     ----------------------------------------------\n");
	}
	mprintf("     ----------------------------------------------\n\n");
	
	mprintf("    Saving normal coordinate spectra...\n");
#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "%s_spectra.csv", basename);
#else
	sprintf(name, "%s_spectra.csv", basename);
#endif
	FILE *specFile = OpenFileWrite(name, false);
	fprintf(specFile, "#Wavenumber (cm^-1);");
	for(i = 0; i < 3*_atomCount; i++) {
		fprintf(specFile, " Mode %d (K*cm);", i+1);
	}
	fprintf(specFile, " Trace (K*cm)\n");
	for(i = 0; i < _specSize; i++) {
		fprintf(specFile, "%.2f;", _specResolution * i);
		float value = 0.0f;
		for(j = 0; j < 3*_atomCount; j++) {
			fprintf(specFile, " %.8G;", ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*sortIndex->GetAt(j) - sortIndex->GetAt(j)*(sortIndex->GetAt(j)-1)/2 + sortIndex->GetAt(j)])->GetAt(i));
			value += ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*sortIndex->GetAt(j) - sortIndex->GetAt(j)*(sortIndex->GetAt(j)-1)/2 + sortIndex->GetAt(j)])->GetAt(i);
		}
		fprintf(specFile, " %.8G\n", value);
	}
	fclose(specFile);
	
	mprintf("    Creating Molden file...\n");
#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "%s.molden", basename);
#else
	sprintf(name, "%s.molden", basename);
#endif
	FILE *moldenFile = OpenFileWrite(name, false);
	fprintf(moldenFile, "[Molden Format]\n");
	fprintf(moldenFile, "[FREQ]\n");
	for(i = 0; i < 3*_atomCount; i++) {
		fprintf(moldenFile, "%f\n", i < 6 ? 0.0f : centers->GetAt(sortIndex->GetAt(i)));
	}
	fprintf(moldenFile, "[FR-COORD]\n");
	for(i = 0; i < _atomCount; i++) {
		fprintf(moldenFile, "%s %f %f %f\n", ((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_iType]])->m_sName, ((CxVec3Array *)_centroidCoords[0])->GetAt(i)[0] / 100.0f * 1.889726f, ((CxVec3Array *)_centroidCoords[0])->GetAt(i)[1] / 100.0f * 1.889726f, ((CxVec3Array *)_centroidCoords[0])->GetAt(i)[2] / 100.0f * 1.889726f); 
	}
	fprintf(moldenFile, "[FR-NORM-COORD]\n");
	for(i = 0; i < 3*_atomCount; i++) {
		fprintf(moldenFile, "vibration %d\n", i+1);
		for(j = 0; j < 3*_atomCount; j+=3) {
			float mass = sqrtf(((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[j / 3])->m_iType]])->m_pElement->m_fMass);
			fprintf(moldenFile, "%f %f %f\n", transMatrix->GetAt(sortIndex->GetAt(i) * 3*_atomCount + j) / mass, transMatrix->GetAt(sortIndex->GetAt(i) * 3*_atomCount + j + 1) / mass, transMatrix->GetAt(sortIndex->GetAt(i) * 3*_atomCount + j + 2) / mass);
		}
	}
// 	fprintf(moldenFile, "[INT]\n");
// 	for(int i = 0; i < 3*_atomCount; i++) {
// 		fprintf(moldenFile, "%f\n", integrals->GetAt(i));
// 	}
	fclose(moldenFile);
	
	mprintf("    Creating NMD file...\n");
#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "%s.nmd", basename);
#else
	sprintf(name, "%s.nmd", basename);
#endif
	FILE *nmdFile = OpenFileWrite(name, false);
	fprintf(nmdFile, "atomnames");
	for(i = 0; i < _atomCount; i++) {
		fprintf(nmdFile, " %s", ((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_iType]])->m_sName);
	}
	fprintf(nmdFile, "\n");
	fprintf(nmdFile, "coordinates");
	for(i = 0; i < _atomCount; i++) {
		for(j = 0; j < 3; j++) {
			fprintf(nmdFile, " %f", ((CxVec3Array *)_centroidCoords[0])->GetAt(i)[j] / 100.0f);
		}
	}
	fprintf(nmdFile, "\n");
	for(i = 0; i < 3*_atomCount; i++) {
		fprintf(nmdFile, "mode %d", i+1);
		for(j = 0; j < 3*_atomCount; j+=3) {
			float mass = sqrtf(((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[j / 3])->m_iType]])->m_pElement->m_fMass);
			for(k = 0; k < 3; k++) {
				fprintf(nmdFile, " %f", transMatrix->GetAt(sortIndex->GetAt(i) * 3*_atomCount + j + k) / mass);
			}
		}
		fprintf(nmdFile, "\n");
	}
	fclose(nmdFile);
	
// 	if(_calcIR) {
// 		mprintf("    Calculating dipole auto correlation...\n");
// 		CxFloatArray *dipoleAC;
// 		try { dipoleAC = new CxFloatArray(); } catch(...) { dipoleAC = NULL; }
// 		if(dipoleAC == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		dipoleAC->SetSize(_correlationDepth);
// 		for(int i = 0; i < _correlationDepth; i++)
// 			dipoleAC->GetAt(i) = 0.0f;
// 		
// 		CxFloatArray *tmp;
// 		try { tmp = new CxFloatArray(); } catch(...) { tmp = NULL; }
// 		if(tmp == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		tmp->SetSize(_numSteps-2);
// 		CxFloatArray *tmp2;
// 		try { tmp2 = new CxFloatArray(); } catch(...) { tmp2 = NULL; }
// 		if(tmp2 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		tmp2->SetSize(_numSteps-2);
// 		CxFloatArray *tmp3;
// 		try { tmp3 = new CxFloatArray(); } catch(...) { tmp3 = NULL; }
// 		if(tmp3 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		tmp3->SetSize(_correlationDepth);
// 		
// // 		FILE *testFile = fopen("test.dat", "w");
// // 		for(int i = 0; i < _numSteps - 2; i++) {
// // 			fprintf(testFile, "%d %f %f %f\n", i, ((CxVec3Array *)_dipoleDerivativeCache[0])->GetAt(i)[0], ((CxVec3Array *)_dipoleDerivativeCache[0])->GetAt(i)[1], ((CxVec3Array *)_dipoleDerivativeCache[0])->GetAt(i)[2]);
// // 		}
// // 		fclose(testFile);
// 		
// 		mprintf(WHITE, "     [");
// 		float step = (float)(_showMolCount * _permutationCount * 3) / 20.0f;
// 		int c = 0;
// 		for(int i = 0; i < _showMolCount; i++) {
// 			for(int j = 0; j < _permutationCount; j++) {
// 				for(int k = 0; k < 3; k++) {
// 					if(fmodf((float)c++, step) < 1.0f)
// 						mprintf(WHITE, "#");
// 					for(int l = 0; l < _numSteps - 2; l++) {
// 						tmp->GetAt(l) = ((CxFloatArray *)_globalProb[i])->GetAt(l+1) * (float)(exp(-(double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(l+1) * (double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(l+1) / 2.0 / (double)_gaussianWidth / (double)_gaussianWidth) / probSum[i*_numSteps + l+1]) * ((CxVec3Array *)_dipoleDerivativeCache[i*_permutationCount + j])->GetAt(l)[k];
// 						tmp2->GetAt(l) = ((CxVec3Array *)_dipoleDerivativeCache[i*_permutationCount + j])->GetAt(l)[k];
// 					}
// 					cc->CrossCorrelate(tmp, tmp2, tmp3);
// 					for(int l = 0; l < _correlationDepth; l++) {
// 						dipoleAC->GetAt(l) += tmp3->GetAt(l);
// 					}
// 				}
// 			}
// 		}
// 		mprintf(WHITE, "]\n");
// 		for(int i = 0; i < _correlationDepth; i++) {
// 			dipoleAC->GetAt(i) /= _showMolCount;
// 		}
// 		
// 		mprintf("    Fourier transforming dipole auto correlation...\n");
// 		CFFT *fft;
// 		try { fft = new CFFT(); } catch(...) { fft = NULL; }
// 		if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		CxFloatArray *tmp4;
// 		try { tmp4 = new CxFloatArray(); } catch(...) { tmp4 = NULL; }
// 		if(tmp4 == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		
// 		mprintf(WHITE, "     [");
// 		tmp4->CopyFrom(dipoleAC);
// 		if(_applyWindowFunction) {
// 			for(int i = 0; i < tmp4->GetSize(); i++) {
// 				tmp4->GetAt(i) *= powf(cosf((float)i / tmp4->GetSize() / 2.0f * Pi), 2.0f);
// 			}
// 		}
// 		if(_zeroPadding > 0) {
// 			for(int i = 0; i < _zeroPadding; i++) {
// 				tmp4->Add(0.0f);
// 			}
// 		}
// 		if(_applyMirroring) {
// 			int oldSize = tmp4->GetSize();
// 			tmp4->SetSize(2 * oldSize);
// 			for(int i = 1; i < oldSize; i++) {
// 				tmp4->GetAt(oldSize + i) = tmp4->GetAt(oldSize - i);
// 			}
// 		}
// 		fft->PrepareFFT_C2C(tmp4->GetSize());
// 		for(int i = 0; i < tmp4->GetSize(); i++) {
// 			fft->m_pInput[2*i] = tmp4->GetAt(i);
// 			fft->m_pInput[2*i+1] = 0.0f;
// 		}
// 		fft->DoFFT();
// 		dipoleAC->SetSize(_specSize);
// 		for(int i = 0; i < _specSize; i++) {
// 			dipoleAC->GetAt(i) = fft->m_pOutput[2*i] * 1523615.0f;
// 		}
// 		mprintf(WHITE, "#]\n");
// 		
// 		mprintf("    Saving infrared spectrum...\n");
// 		snprintf(name, BUF_SIZE, "%s_ir.dat", basename);
// 		FILE *irFile = OpenFileWrite(name, false);
// 		for(int i = 0; i < _specSize; i++) {
// 			fprintf(irFile, "%10.4f %20.8f\n", _specResolution*i, dipoleAC->GetAt(i));
// 		}
// 		fclose(irFile);
// 		
// 		mprintf("    Calculating linear least-squares fit...\n");
// 		CxFloatArray *diagonalCoefficients;
// 		try { diagonalCoefficients = new CxFloatArray(); } catch(...) { diagonalCoefficients = NULL; }
// 		if(diagonalCoefficients == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		diagonalCoefficients->SetSize(3*_atomCount);
// 		CxFloatArray *xMatrix;
// 		try { xMatrix = new CxFloatArray(); } catch(...) { xMatrix = NULL; }
// 		if(xMatrix == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		xMatrix->SetSize(3*_atomCount * _specSize);
// 		CxFloatArray *qtMatrix;
// 		try { qtMatrix = new CxFloatArray(); } catch(...) { qtMatrix = NULL; }
// 		if(qtMatrix == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		qtMatrix->SetSize(3*_atomCount * _specSize);
// 		CxFloatArray *rMatrix;
// 		try { rMatrix = new CxFloatArray(); } catch(...) { rMatrix = NULL; }
// 		if(rMatrix == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		rMatrix->SetSize(3*_atomCount * _specSize);
// 		CxFloatArray *qyVector;
// 		try { qyVector = new CxFloatArray(); } catch(...) { qyVector = NULL; }
// 		if(qyVector == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		qyVector->SetSize(3*_atomCount);
// 		
// 		for(int i = 0; i < 3*_atomCount; i++) {
// 			for(int j = 0; j < _specSize; j++) {
// 				xMatrix->GetAt(j * 3*_atomCount + i) = ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*sortIndex->GetAt(i) - sortIndex->GetAt(i)*(sortIndex->GetAt(i)-1)/2 + sortIndex->GetAt(i)])->GetAt(j);
// 			}
// 		}
// // 		qrDecompose(_specSize, 3*_atomCount, xMatrix, qtMatrix, rMatrix);
// 		
// 		for(int i = 0; i < 3*_atomCount; i++) {
// 			qyVector->GetAt(i) = 0.0f;
// 			for(int j = 0; j < _specSize; j++) {
// 				qyVector->GetAt(i) += qtMatrix->GetAt(i * _specSize + j) * dipoleAC->GetAt(j);
// 			}
// 		}
// 		for(int i = 3*_atomCount - 1; i >= 0; i--) {
// 			diagonalCoefficients->GetAt(i) = qyVector->GetAt(i);
// 			for(int j = i+1; j < 3*_atomCount; j++) {
// 				diagonalCoefficients->GetAt(i) -= diagonalCoefficients->GetAt(j) * rMatrix->GetAt(i * 3*_atomCount + j);
// 			}
// 			diagonalCoefficients->GetAt(i) /= rMatrix->GetAt(i * 3*_atomCount + i);
// 		}
// 		
// 		mprintf("\n     -------------------------\n");
// 		mprintf("     Mode   Intensity (km/mol)\n");
// 		mprintf("     -------------------------\n");
// 		for(int i = 0; i < 6; i++) {
// 			mprintf("     %4d %20.6g\n", i+1, diagonalCoefficients->GetAt(i) / 1000.0f);
// 		}
// 		mprintf("     -------------------------\n");
// 		for(int i = 6; i < 3*_atomCount; i++) {
// 			mprintf("     %4d %20.6g\n", i+1, diagonalCoefficients->GetAt(i) / 1000.0f);
// 		}
// 		mprintf("     -------------------------\n\n");
// 		
// 		mprintf("    Saving reconstructed infrared spectrum...\n");
// 		CxFloatArray *irSpec;
// 		try { irSpec = new CxFloatArray(); } catch(...) { irSpec = NULL; }
// 		if(irSpec == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		irSpec->SetSize(_specSize);
// 		for(int i = 0; i < _specSize; i++)
// 			irSpec->GetAt(i) = 0.0f;
// 		for(int i = 0; i < 3*_atomCount; i++) {
// 			for(int j = 0; j < _specSize; j++) {
// 				irSpec->GetAt(j) += diagonalCoefficients->GetAt(i) * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*sortIndex->GetAt(i) - sortIndex->GetAt(i)*(sortIndex->GetAt(i)-1)/2 + sortIndex->GetAt(i)])->GetAt(j);
// 			}
// 		}
// 		snprintf(name, BUF_SIZE, "%s_ir_fit_diagonal.dat", basename);
// 		irFile = OpenFileWrite(name, false);
// 		for(int i = 0; i < _specSize; i++) {
// 			fprintf(irFile, "%10.4f %20.8f\n", _specResolution*i, irSpec->GetAt(i));
// 		}
// 		fclose(irFile);
// 
// 		delete dipoleAC;
// 		
// 		mprintf("    Calculating dipole-velocity cross correlations...\n");
// 		CxObArray cdMatrix;
// 		for(int i = 0; i < 3 * 3*_atomCount; i++) {
// 			CxFloatArray *a;
// 			try { a = new CxFloatArray(); } catch(...) { a = NULL; }
// 			if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			a->SetSize(_correlationDepth);
// 			for(int j = 0; j < _correlationDepth; j++) {
// 				a->GetAt(j) = 0.0f;
// 			}
// 			cdMatrix.Add(a);
// 		}
// 		
// 		mprintf(WHITE, "     [");
// 		step = (float)(_showMolCount * _permutationCount * 3 * 3*_atomCount) / 20.0f;
// 		c = 0;
// 		for(int i = 0; i < _showMolCount; i++) {
// 			for(int j = 0; j < _permutationCount; j++) {
// 				for(int k = 0; k < _atomCount; k++) {
// 					for(int l = 0; l < 3; l++) {
// 						for(int m = 0; m < 3; m++) {
// 							if(fmodf((float)c++, step) < 1.0f)
// 								mprintf(WHITE, "#");
// 							for(int n = 0; n < _numSteps - 2; n++) {
// 								tmp->GetAt(n) = ((CxFloatArray *)_globalProb[i])->GetAt(n+1) * (float)(exp(-(double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(n+1) * (double)((CxFloatArray *)_distanceTimedev[i*_permutationCount + j])->GetAt(n+1) / 2.0 / (double)_gaussianWidth / (double)_gaussianWidth) / probSum[i*_numSteps + n+1]) * ((CxVec3Array *)_dipoleDerivativeCache[i*_permutationCount + j])->GetAt(n)[m];
// 								tmp2->GetAt(n) = ((CxVec3Array *)_velocityCache[i*_permutationCount*_atomCount + j*_atomCount + k])->GetAt(n)[l] * sqrtf(_weights[k]);
// 							}
// 							cc->CrossCorrelate(tmp, tmp2, tmp3);
// 							for(int n = 0; n < _correlationDepth; n++) {
// 								((CxFloatArray *)cdMatrix[m * 3*_atomCount + 3 * k + l])->GetAt(n) += tmp3->GetAt(n);
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 		mprintf(WHITE, "]\n");
// 		for(int i = 0; i < 3*_atomCount; i++) {
// 			for(int j = 0; j < 3; j++) {
// 				for(int k = 0; k < _correlationDepth; k++) {
// 					((CxFloatArray *)cdMatrix[j * 3*_atomCount  + i])->GetAt(k) /= _showMolCount;
// 				}
// 			}
// 		}
// 		
// 		mprintf("    Fourier transforming dipole-velocity cross correlations...\n");
// 		mprintf(WHITE, "     [");
// 		step = (float)(3 * 3*_atomCount) / 20.0f;
// 		c = 0;
// 		for(int i = 0; i < 3*_atomCount; i++) {
// 			for(int j = 0; j < 3; j++) {
// 				if(fmodf((float)c++, step) < 1.0f)
// 					mprintf(WHITE, "#");
// 				tmp4->CopyFrom((CxFloatArray *)cdMatrix[j * 3*_atomCount + i]);
// 				if(_applyWindowFunction) {
// 					for(int k = 0; k < tmp4->GetSize(); k++) {
// 						tmp4->GetAt(k) *= powf(cosf((float)k / tmp4->GetSize() / 2.0f * Pi), 2.0f);
// 					}
// 				}
// 				if(_zeroPadding > 0) {
// 					for(int k = 0; k < _zeroPadding; k++) {
// 						tmp4->Add(0.0f);
// 					}
// 				}
// 				if(_applyMirroring) {
// 					int oldSize = tmp4->GetSize();
// 					tmp4->SetSize(2 * oldSize);
// 					for(int k = 1; k < oldSize; k++) {
// 						tmp4->GetAt(oldSize + k) = tmp4->GetAt(oldSize - k);
// 					}
// 				}
// 				fft->PrepareFFT_C2C(tmp4->GetSize());
// 				for(int k = 0; k < tmp4->GetSize(); k++) {
// 					fft->m_pInput[2*k] = tmp4->GetAt(k);
// 					fft->m_pInput[2*k+1] = 0.0f;
// 				}
// 				fft->DoFFT();
// 				((CxFloatArray *)cdMatrix[j * 3*_atomCount + i])->SetSize(_specSize);
// 				for(int k = 0; k < _specSize; k++) {
// 					((CxFloatArray *)cdMatrix[j * 3*_atomCount + i])->GetAt(k) = fft->m_pOutput[2*k];
// 				}
// 			}
// 		}
// 		mprintf(WHITE, "]\n");
// 		
// 		FILE *testFile = fopen("test.dat", "w");
// 		for(int i = 0; i < _specSize; i++) {
// 			fprintf(testFile, "%10.4f", _specResolution * i);
// 			for(int j = 0; j < 3*_atomCount; j++) {
// 				fprintf(testFile, " %20.8f", ((CxFloatArray *)cdMatrix[j])->GetAt(i));
// 			}
// 			fprintf(testFile, "\n");
// 		}
// 		fclose(testFile);
// 		
// 		delete fft;
// 		
// 		CxObArray cdtMatrix;
// 		for(int i = 0; i < 3 * 3*_atomCount; i++) {
// 			CxFloatArray *a;
// 			try { a = new CxFloatArray(); } catch(...) { a = NULL; }
// 			if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			a->SetSize(_correlationDepth);
// 			for(int j = 0; j < _correlationDepth; j++) {
// 				a->GetAt(j) = 0.0f;
// 			}
// 			cdtMatrix.Add(a);
// 		}
// 		for(int i = 0; i < 3; i++) {
// 			for(int j = 0; j < 3*_atomCount; j++) {
// 				for(int k = 0; k < 3*_atomCount; k++) {
// 					for(int l = 0; l < _specSize; l++) {
// // 						((CxFloatArray *)cdtMatrix[i * 3*_atomCount + j])->GetAt(l) += ((CxFloatArray *)cdMatrix[i * 3*_atomCount + k])->GetAt(l) * transMatrix->GetAt(sortIndex->GetAt(j) * 3*_atomCount + k);
// 						((CxFloatArray *)cdtMatrix[i * 3*_atomCount + j])->GetAt(l) += ((CxFloatArray *)cdMatrix[i * 3*_atomCount + k])->GetAt(l) * transMatrix->GetAt(j * 3*_atomCount + k);
// 					}
// 				}
// 			}
// 		}
// 		
// 		FILE *test2File = fopen("test2.dat", "w");
// 		for(int i = 0; i < _specSize; i++) {
// 			fprintf(test2File, "%10.4f", _specResolution * i);
// 			for(int j = 0; j < 3*_atomCount; j++) {
// 				fprintf(test2File, " %20.8f", ((CxFloatArray *)cdtMatrix[3*_atomCount + j])->GetAt(i));
// 			}
// 			fprintf(test2File, "\n");
// 		}
// 		fclose(test2File);
// 		
// 		mprintf("    Calculating linear least-squares fit...\n");
// 		CxFloatArray *coefficientMatrix;
// 		try { coefficientMatrix = new CxFloatArray(); } catch(...) { coefficientMatrix = NULL; }
// 		if(coefficientMatrix == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		coefficientMatrix->SetSize(3 * 3*_atomCount);
// 		xMatrix->SetSize(3*_atomCount * 3*_atomCount * _specSize);
// 		qtMatrix->SetSize(3*_atomCount * 3*_atomCount * _specSize);
// 		rMatrix->SetSize(3*_atomCount * 3*_atomCount * _specSize);
// 		qyVector->SetSize(3*_atomCount);
// 		
// 		for(int i = 0; i < 3*_atomCount; i++) {
// 			for(int j = 0; j < _specSize; j++) {
// 				for(int k = 0; k < 3*_atomCount; k++) {
// 					if(i <= k)
// 						xMatrix->GetAt((i*_specSize + j) * 3*_atomCount + k) = ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*i - i*(i-1)/2 + k])->GetAt(j);
// // 						xMatrix->GetAt((i*_specSize + j) * 3*_atomCount + k) = ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*sortIndex->GetAt(i) - sortIndex->GetAt(i)*(sortIndex->GetAt(i)-1)/2 + sortIndex->GetAt(k)])->GetAt(j);
// 					else
// 						xMatrix->GetAt((i*_specSize + j) * 3*_atomCount + k) = ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*k - k*(k-1)/2 + i])->GetAt(j);
// // 						xMatrix->GetAt((i*_specSize + j) * 3*_atomCount + k) = ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*sortIndex->GetAt(k) - sortIndex->GetAt(k)*(sortIndex->GetAt(k)-1)/2 + sortIndex->GetAt(i)])->GetAt(j);
// 				}
// 			}
// 		}
// // 		qrDecompose(3*_atomCount * _specSize, 3*_atomCount, xMatrix, qtMatrix, rMatrix);
// 		
// // 		FILE *testFile = fopen("test.dat", "w");
// // 		for(int i = 0; i < 3*_atomCount; i++) {
// // 			for(int j = 0; j < 3*_atomCount; j++) {
// // 				fprintf(testFile, " %f", xMatrix->GetAt(i * 3*_atomCount * _specSize + j));
// // 			}
// // 			fprintf(testFile, "\n");
// // 		}
// // 		fclose(testFile);
// 		
// 		for(int i = 0; i < 3; i++) {
// 			for(int j = 0; j < 3*_atomCount; j++) {
// 				qyVector->GetAt(j) = 0.0f;
// 				for(int k = 0; k < 3*_atomCount; k++) {
// 					for(int l = 0; l < _specSize; l++) {
// 						qyVector->GetAt(j) += qtMatrix->GetAt(j * 3*_atomCount * _specSize + k * _specSize + l) * ((CxFloatArray *)cdtMatrix[i * 3*_atomCount + k])->GetAt(l);
// 					}
// 				}
// 			}
// 			for(int j = 3*_atomCount - 1; j >= 0; j--) {
// 				coefficientMatrix->GetAt(i * 3*_atomCount + j) = qyVector->GetAt(j);
// 				for(int k = j+1; k < 3*_atomCount; k++) {
// 					coefficientMatrix->GetAt(i * 3*_atomCount + j) -= coefficientMatrix->GetAt(i * 3*_atomCount + k) * rMatrix->GetAt(j * 3*_atomCount + k);
// 				}
// 				coefficientMatrix->GetAt(i * 3*_atomCount + j) /= rMatrix->GetAt(j * 3*_atomCount + j);
// 			}
// 		}
// 		
// 		delete xMatrix;
// 		delete qtMatrix;
// 		delete rMatrix;
// 		delete qyVector;
// 		
// 		mprintf("    Saving reconstructed infrared spectrum...\n");
// 		irSpec->SetSize(_specSize);
// 		for(int i = 0; i < _specSize; i++)
// 			irSpec->GetAt(i) = 0.0f;
// 		for(int i = 0; i < 3; i++) {
// 			for(int j = 0; j < 3*_atomCount; j++) {
// 				for(int k = 0; k < 3*_atomCount; k++) {
// 					for(int l = 0; l < _specSize; l++) {
// 						if(j <= k)
// 							irSpec->GetAt(l) += coefficientMatrix->GetAt(i * 3*_atomCount + j) * coefficientMatrix->GetAt(i * 3*_atomCount + k) * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + k])->GetAt(l);
// // 							irSpec->GetAt(l) += coefficientMatrix->GetAt(i * 3*_atomCount + j) * coefficientMatrix->GetAt(i * 3*_atomCount + k) * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*sortIndex->GetAt(j) - sortIndex->GetAt(j)*(sortIndex->GetAt(j)-1)/2 + sortIndex->GetAt(k)])->GetAt(l);
// 						else
// 							irSpec->GetAt(l) += coefficientMatrix->GetAt(i * 3*_atomCount + j) * coefficientMatrix->GetAt(i * 3*_atomCount + k) * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*j - j*(j-1)/2 + k])->GetAt(l);
// // 							irSpec->GetAt(l) += coefficientMatrix->GetAt(i * 3*_atomCount + j) * coefficientMatrix->GetAt(i * 3*_atomCount + k) * ((CxFloatArray *)ccMatrix[(3*_atomCount-1)*sortIndex->GetAt(j) - sortIndex->GetAt(j)*(sortIndex->GetAt(j)-1)/2 + sortIndex->GetAt(k)])->GetAt(l);
// 					}
// 				}
// 			}
// 		}
// 		snprintf(name, BUF_SIZE, "%s_ir_fit_linear.dat", basename);
// 		irFile = OpenFileWrite(name, false);
// 		for(int i = 0; i < _specSize; i++) {
// 			fprintf(irFile, "%10.4f %20.8f\n", _specResolution*i, irSpec->GetAt(i));
// 		}
// 		fclose(irFile);
// 		
// 		delete irSpec;
// 		
// // 		delete temp;
// // 		delete temp2;
// // 		delete temp3;
// 	}
	
	delete cc;
	
	delete integrals;
	delete centers;
	delete sortIndex;
}

bool CReferenceStructure::recognizeMolecule(CTimeStep *ts, int showMol, const char *basename) {
	int i, j;
	bool ok = true;
	unsigned int *stack;
	try { stack = new unsigned int[ts->m_iGesAtomCount]; } catch(...) { stack = NULL; }
	if(stack == NULL) NewException((double)ts->m_iGesAtomCount*sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	bool *used;
	try { used = new bool[ts->m_iGesAtomCount]; } catch(...) { used = NULL; }
	if(used == NULL) NewException((double)ts->m_iGesAtomCount*sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	for(i = 0; (unsigned int)i < ts->m_iGesAtomCount; i++) {
		used[i] = false;
	}
	
	CxByteArray baAtomIndex;
	baAtomIndex.SetSize(ts->m_iGesAtomCount);
	for(i = 0; (unsigned int)i < ts->m_iGesAtomCount; i++) {
//		char buf[64];
		CxString buf;

//		strncpy(buf, (char *)(ts->m_paLabels[i]), 64);
//		buf[63] = 0;

		buf.strcpy((char *)(ts->m_paLabels[i]));
		ReplaceDigits(&buf);

		baAtomIndex[i] = 255;
		for(j = 0; j < g_oaAtoms.GetSize(); j++) {
			if(mystricmp(buf, ((CAtom *)g_oaAtoms[j])->m_sName) == 0) {
				CAtom *a = (CAtom *)g_oaAtoms[j];
				while(a->m_pMergedTo != NULL)
					a = a->m_pMergedTo;
				baAtomIndex[i] = a->m_iIndex;
				break;
			}
		}
		if(baAtomIndex[i] == 255) {
			mprintf(RED, "Atom type \"%s\" not known\n", (const char*)buf);
			ok = false;
		}
	}
	
	if(ok) {
		recognizeMoleculeRecursion(0, used, 0, stack, ts, baAtomIndex);
		for(i = 0; (unsigned int)i < ts->m_iGesAtomCount; i++) {
			if(!used[i]) {
				mprintf(RED, "Some atoms are not connected to the molecule\n");
				ok = false;
				break;
			}
		}
	}
	
	if(ok) {
		mprintf("    Recognized one molecule\n");
		
		mprintf("    Sorting atom types...\n");
		for(i = 0; i < _singleMol->m_baAtomIndex.GetSize() - 1; i++) {
			int a = -1;
			int b = 999;
			for(j = i; j < _singleMol->m_baAtomIndex.GetSize(); j++) {
				if(_singleMol->m_baAtomIndex[j] < b) {
					b = _singleMol->m_baAtomIndex[j];
					a = j;
				}
			}
			if((a != -1) && (a != i)) {
				b = _singleMol->m_baAtomIndex[a];
				_singleMol->m_baAtomIndex[a] = _singleMol->m_baAtomIndex[i];
				_singleMol->m_baAtomIndex[i] = b;
				for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
					if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType == i)
						((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType = a;
					else if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType == a)
						((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType = i;
				}
			}
		}
		
		mprintf("    Setting up bond list...\n");
		for(i = 0; i < _singleMol->m_laBonds.GetSize() / 2; i++) {
			int a = _singleMol->m_laBonds[2*i];
			int b = _singleMol->m_laBonds[2*i+1];
			int c = -1;
			for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
				if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iOffset == a) {
					c = j;
					break;
				}
			}
			if(c == -1) {
				mprintf(RED, "Weird error.\n");
				abort();
			}
			int d = -1;
			for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
				if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iOffset == b) {
					d = j;
					break;
				}
			}
			if(d == -1) {
				mprintf(RED, "Weird error.\n");
				abort();
			}
			((CMolAtom *)_singleMol->m_oaMolAtoms[c])->m_oaBonds.Add((CMolAtom *)_singleMol->m_oaMolAtoms[d]);
			((CMolAtom *)_singleMol->m_oaMolAtoms[d])->m_oaBonds.Add((CMolAtom *)_singleMol->m_oaMolAtoms[c]);
		}
		
		mprintf("    Building atom codes...\n");
		_singleMol->BuildAtomCodes();
		
		mprintf("    Creating topological atom order...\n");
		for(i = 0; i < _singleMol->m_baAtomIndex.GetSize(); i++) {
			CxIntArray *a;
			try { a = new CxIntArray(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			_singleMol->m_oaAtomOffset.Add(a);
			for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
				if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType == i) {
					((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iNumber = a->GetSize();
					a->Add(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iOffset);
				}
			}
		}
		
// 		mprintf("    Creating bond list...\n");
// 		mprintf("    Creating angle list...\n");
		mprintf("    Comparing to global molecule type...\n");
		bool diff = false;
		CSingleMolecule *compMol = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[showMol])->m_laSingleMolIndex[0]];
		if(_singleMol->m_oaMolAtoms.GetSize() != compMol->m_oaMolAtoms.GetSize()) {
			diff = true;
		} else {
			for(i = 0; i < _singleMol->m_oaMolAtoms.GetSize(); i++) {
// 				mprintf(RED, "%f %f\n", ((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_fAtomCode, ((CMolAtom *)compMol->m_oaMolAtoms[i])->m_fAtomCode);
				if(((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_fAtomCode != ((CMolAtom *)compMol->m_oaMolAtoms[i])->m_fAtomCode) {
					diff = true;
					break;
				}
			}
		}
		if(!diff) {
			_singleMol->m_iMolType = showMol;
			mprintf("    The reference structure belongs to molecule type %d (%s)\n", _singleMol->m_iMolType+1, ((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_sName);
			_atomCount = _singleMol->m_oaMolAtoms.GetSize();
		} else {
			ok = false;
			mprintf(RED, "The reference structure does not belong to molecule type %d (%s)\n", _singleMol->m_iMolType+1, ((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_sName);
		}
// 		mprintf("    Refining ring systems...\n");
	}
	mprintf("\n");
	
	if(ok) {
		mprintf("    Looking for equivalent atoms...\n");
		CxIntArray permutationActions;
		askPermutations(&permutationActions);
		
		mprintf("\n    Creating permutations of equivalent atoms...\n");
		
		CxIntArray *nullPerm;
		try { nullPerm = new CxIntArray(); } catch(...) { nullPerm = NULL; }
		if(nullPerm == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		createPermutationsRecursion(0, nullPerm, &permutationActions);
		
		_permutationCount = _permutations.GetSize();
		for(i = 0; i < _permutationCount; i++) {
			CxVec3Array *centroidCoord;
			try { centroidCoord = new CxVec3Array(); } catch(...) { centroidCoord = NULL; }
			if(centroidCoord == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			for(j = 0; j < _atomCount; j++) {
				centroidCoord->Add(_refTimestep->m_vaCoords[((CMolAtom *)_singleMol->m_oaMolAtoms[((CxIntArray *)_permutations[i])->GetAt(j)])->m_iOffset]);
			}
			_centroidCoords.Add(centroidCoord);
		}
		mprintf("    Created %d permutations\n", _permutationCount);
		for(i = 0; i < _permutationCount; i++) {
			mprintf("      %d: ", i+1);
			for(j = 0; j < _atomCount; j++) {
				mprintf(" %d", ((CxIntArray *)_permutations[i])->GetAt(j)+1);
			}
			mprintf("\n");
		}
		mprintf("\n");
		
		if(g_bAdvanced2) {
			if(AskYesNo("    Save xyz files for all permutations (y/n)? [no] ", false)) {
				for(i = 0; i < _centroidCoords.GetSize(); i++) {
					char name[BUF_SIZE];
#ifdef TARGET_LINUX
					snprintf(name, BUF_SIZE, "%s_p%d_reference.xyz", basename, i+1);
#else
					sprintf(name, "%s_p%d_reference.xyz", basename, i+1);
#endif
					FILE *permutationFile = OpenFileWrite(name, false);
					fprintf(permutationFile, "%d\n\n", _singleMol->m_oaMolAtoms.GetSize());
					for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
						fprintf(permutationFile, "%s %f %f %f\n", ((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType]])->m_sName, ((CxVec3Array *)_centroidCoords[i])->GetAt(j)[0] / 100.0f, ((CxVec3Array *)_centroidCoords[i])->GetAt(j)[1] / 100.0f, ((CxVec3Array *)_centroidCoords[i])->GetAt(j)[2] / 100.0f); 
					}
					fclose(permutationFile);
				}
			}
			mprintf("\n");
		}
		
		_weightSum = 0.0f;
		_weights.SetSize(_singleMol->m_oaMolAtoms.GetSize());
		_centroid = CxVector3(0.0f, 0.0f, 0.0f);
		for(i = 0; i < _singleMol->m_oaMolAtoms.GetSize(); i++) {
			_weights[i] = ((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_iType]])->m_pElement->m_fMass;
			_centroid += ((CxVec3Array *)_centroidCoords[0])->GetAt(i) * _weights[i];
			_weightSum += _weights[i];
		}
		_centroid /= _weightSum;
		for(i = 0; i < _centroidCoords.GetSize(); i++) {
			for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
				((CxVec3Array *)_centroidCoords[i])->GetAt(j) -= _centroid;
			}
		}
	}
	
	delete[] used;
	delete[] stack;
	
	return ok;
}

void CReferenceStructure::recognizeMoleculeRecursion(unsigned int index, bool *used, int depth, unsigned int *stack, CTimeStep *ts, CxByteArray &baAtomIndex) {
	int i, j;
	if(g_bVerbose) {
		mprintf("    ");
		for(i = 1; i < depth; i++)
			mprintf("    ");
		if(depth > 0)
			mprintf(WHITE, "\\--");
		mprintf(CYAN, "%s", ((CAtom *)g_oaAtoms[baAtomIndex[index]])->m_sName);
		mprintf("(%d)", index+1);
	}
	
	stack[depth] = index;
	used[index] = true;
	
	bool found = false;
	for(i = 0; i < _singleMol->m_baAtomIndex.GetSize(); i++) {
		if(_singleMol->m_baAtomIndex[i] == baAtomIndex[index]) {
			CMolAtom *molAtom;
			try { molAtom = new CMolAtom(); } catch(...) { molAtom = NULL; }
			if(molAtom == NULL) NewException((double)sizeof(CMolAtom), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			molAtom->m_iOffset = index;
			molAtom->m_iType = i;
			_singleMol->m_oaMolAtoms.Add(molAtom);
			found = true;
			break;
		}
	}
	if(!found) {
		CMolAtom *molAtom;
		try { molAtom = new CMolAtom(); } catch(...) { molAtom = NULL; }
		if(molAtom == NULL) NewException((double)sizeof(CMolAtom), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		molAtom->m_iOffset = index;
		molAtom->m_iType = _singleMol->m_baAtomIndex.GetSize();
		_singleMol->m_oaMolAtoms.Add(molAtom);
		_singleMol->m_baAtomIndex.Add(baAtomIndex[index]);
	}
	
	int numNeighbors = 0;
	int neighborList[MAX_BONDS];
	for(i = 0; (unsigned int)i < ts->m_iGesAtomCount; i++) {
		if((unsigned int)i == index)
			continue;
		if(((CAtom *)g_oaAtoms[baAtomIndex[index]])->m_bExclude)
			continue;
		float distance;
		if(recognizeMoleculeBondRange(ts, index, i, &distance, baAtomIndex)) {
			if(used[i]) {
				if((depth > 0) && ((unsigned int)i != stack[depth-1])) {
					if(_singleMol->m_oaRings.GetSize() < 100) {
						if(g_bVerbose) {
							mprintf(GREEN, " <-- Ring closure: ");
							mprintf("%s(%d)", ((CAtom *)g_oaAtoms[baAtomIndex[stack[depth]]])->m_sName, stack[depth]+1);
						}
						CxIntArray *ring;
						try { ring = new CxIntArray(); } catch(...) { ring = NULL; }
						if(ring == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
						_singleMol->m_oaRings.Add(ring);
						ring->Add(stack[depth]);
						for(j = depth - 1; (stack[j] != (unsigned int)i) && (j >= 0); j--) {
							ring->Add(stack[j]);
							if(g_bVerbose)
								mprintf(" - %s(%d)", ((CAtom *)g_oaAtoms[baAtomIndex[stack[j]]])->m_sName, stack[j]+1);
						}
						ring->Add(i);
						if(g_bVerbose)
							mprintf(" - %s(%d)", ((CAtom *)g_oaAtoms[baAtomIndex[i]])->m_sName, i+1);
					} else {
						mprintf(RED, "More than 100 rings\n");
					}
				}
				continue;
			}
			if(distance < 50.0f) {
				mprintf(RED, "Atoms %s(%d) and %s(%d) are very close to each other: %.4f pm\n", ((CAtom *)g_oaAtoms[baAtomIndex[index]])->m_sName, index + 1, ((CAtom *)g_oaAtoms[baAtomIndex[i]])->m_sName, i + 1, distance);
			}
			neighborList[numNeighbors++] = i;
			if(numNeighbors >= MAX_BONDS) {
				mprintf(RED, "Atom %s(%d) has more than %d bonds\n", ((CAtom *)g_oaAtoms[baAtomIndex[index]])->m_sName, index + 1, MAX_BONDS);
				break;
			}
		}
	}
	
	if(g_bVerbose)
		mprintf("\n");
	
	for(i = 0; i < numNeighbors; i++) {
		_singleMol->m_laBonds.Add(index);
		_singleMol->m_laBonds.Add(neighborList[i]);
		if(!used[neighborList[i]]) {
			recognizeMoleculeRecursion(neighborList[i], used, depth+1, stack, ts, baAtomIndex);
		}
	}
}

bool CReferenceStructure::recognizeMoleculeBondRange(CTimeStep *ts, int i1, int i2, float *distance, CxByteArray &baAtomIndex) {
	float x = ts->m_vaCoords[i1][0] - ts->m_vaCoords[i2][0];
	float y = ts->m_vaCoords[i1][1] - ts->m_vaCoords[i2][1];
	float z = ts->m_vaCoords[i1][2] - ts->m_vaCoords[i2][2];
	
	if(g_bPeriodic)
	{
		if(g_bPeriodicX)
		{
			while(x < -g_fBoxX/2)
				x += g_fBoxX;
			while(x > g_fBoxX/2)
				x -= g_fBoxX;
		}
		if(g_bPeriodicY)
		{
			while(y < -g_fBoxY/2)
				y += g_fBoxY;
			while(y > g_fBoxY/2)
				y -= g_fBoxY;
		}
		if(g_bPeriodicZ)
		{
			while(z < -g_fBoxZ/2)
				z += g_fBoxZ;
			while(z > g_fBoxZ/2)
				z -= g_fBoxZ;
		}
	}
	
	*distance = sqrtf(x*x + y*y + z*z);
	
	return *distance < (((CAtom *)g_oaAtoms[baAtomIndex[i1]])->m_pElement->m_fRadius + ((CAtom *)g_oaAtoms[baAtomIndex[i2]])->m_pElement->m_fRadius) * g_fBondFactor;
}

void CReferenceStructure::askPermutations(CxIntArray *actions) {
	int i;
	int a1 = 0;
	while(a1 < _singleMol->m_oaMolAtoms.GetSize()) {
		int a2 = a1;
		while((a2 < _singleMol->m_oaMolAtoms.GetSize()) && (((CMolAtom *)_singleMol->m_oaMolAtoms[a2])->m_fAtomCode == ((CMolAtom *)_singleMol->m_oaMolAtoms[a1])->m_fAtomCode)) {
			a2++;
		}
		if((a2 - a1) > 1) {
			char buf[BUF_SIZE];
			char buf2[BUF_SIZE];
			size_t remaining = BUF_SIZE;
#ifdef TARGET_LINUX
			remaining -= snprintf(buf, remaining, "Atoms ");
#else
			remaining -= sprintf(buf, "Atoms ");
#endif
			for(i = a1; i < a2; i++) {
				if(remaining < 1)
					break;
#ifdef TARGET_LINUX
				size_t length = snprintf(buf2, remaining, "%s%d ", ((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_iType]])->m_sName, ((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_iNumber + 1);
#else
				size_t length = sprintf(buf2, "%s%d ", ((CAtom *)g_oaAtoms[_singleMol->m_baAtomIndex[((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_iType]])->m_sName, ((CMolAtom *)_singleMol->m_oaMolAtoms[i])->m_iNumber + 1);
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;
				if(i < a2-1) {
#ifdef TARGET_LINUX
					length = snprintf(buf2, remaining, ", ");
#else
					length = sprintf(buf2, ", ");
#endif
					strncat(buf, buf2, remaining - 1);
					remaining -= length;
				}
			}
			strncat(buf, " are equivalent.", remaining - 1);
			mprintf("\n    %s\n", buf);
			int action = AskRangeInteger("    Create no permutation (1), cyclic permutations (2), all permutations (3), or specify certain permutations (4) ? [1] ", 1, 4, 1);
			actions->Add(action);
		}
		a1 = a2;
	}
}

void CReferenceStructure::createPermutationsRecursion(int start, CxIntArray *permutation, CxIntArray *permutationActions) {
	CxString buf;

// 	mprintf(RED, "createPermutationsRecursion (%d):", start);
// 	for(int i = 0; i < permutation->GetSize(); i++) {
// 		mprintf(RED, " %d", permutation->GetAt(i));
// 	}
// 	mprintf(RED, "\n");
	int i, j;
	int a1 = start;
	while(a1 < _singleMol->m_oaMolAtoms.GetSize()) {
		int a2 = a1;
		while((a2 < _singleMol->m_oaMolAtoms.GetSize()) && (((CMolAtom *)_singleMol->m_oaMolAtoms[a2])->m_fAtomCode == ((CMolAtom *)_singleMol->m_oaMolAtoms[a1])->m_fAtomCode)) {
			a2++;
		}
		if((a2 - a1) > 1) {
			CxIntArray newPermutationActions;
			newPermutationActions.CopyFrom(permutationActions);
			int action = newPermutationActions[0];
			newPermutationActions.RemoveAt(0, 1);
			if(action == 1) {
				CxIntArray *newPermutation;
				try { newPermutation = new CxIntArray(); } catch(...) { newPermutation = NULL; }
				if(newPermutation == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				newPermutation->CopyFrom(permutation);
				for(i = a1; i < a2; i++) {
					newPermutation->Add(i);
				}
				createPermutationsRecursion(a2, newPermutation, &newPermutationActions);
				return;
// 				for(int i = a1; i < a2; i++) {
// 					permutation->Add(i);
// 				}
// 				permutationActions->RemoveAt(0, 1);
// 				a1 = a2;
			} else if(action == 2) {
				CxIntArray *newPermutation;
				for(i = 0; i < (a2 - a1); i++) {
					try { newPermutation = new CxIntArray(); } catch(...) { newPermutation = NULL; }
					if(newPermutation == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
					newPermutation->CopyFrom(permutation);
					for(j = a1+i; j < a2; j++)
						newPermutation->Add(j);
					for(j = a1; j < a1+i; j++)
						newPermutation->Add(j);
					createPermutationsRecursion(a2, newPermutation, &newPermutationActions);
				}
				return;
			} else if(action == 3) {
				mprintf(RED, "Not implemented.\n");
				abort();
			} else if(action == 4) {
				CxIntArray *newPermutation;
				while(true) {
					try { newPermutation = new CxIntArray(); } catch(...) { newPermutation = NULL; }
					if(newPermutation == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
					newPermutation->CopyFrom(permutation);
//					char buf[BUF_SIZE];
					const char delim[] = ", ";
					AskString("    Enter permutation: ", &buf, "");
					if(strlen(buf) == 0)
						break;
					char *tok = strtok(buf.GetWritePointer(), delim);
					while(tok != NULL) {
						int n;
						if(sscanf(tok, "%d", &n) == 1) {
							newPermutation->Add(a1+n);
						}
						tok = strtok(NULL, delim);
					}
					createPermutationsRecursion(a2, newPermutation, &newPermutationActions);
				}
				return;
			} else {
				mprintf(RED, "Weird error.\n");
				abort();
			}
		} else {
			permutation->Add(a1);
			a1 = a2;
		}
	}
	_permutations.Add(permutation);
}

double CReferenceStructure::offDiagonalNorm(CxObArray &matrix) {
	double norm = 0.0;
	int i, j, k;
	for(i = 0; i < 3 * _atomCount; i++) {
		for(j = i+1; j < 3 * _atomCount; j++) {
			double integral = 0.0;
			for(k = 0; k < ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetSize(); k++) {
				integral += ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k) * ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + j])->GetAt(k);
			}
// 			mprintf(RED, "Integral (%d, %d): %f\n", i, j, integral);
			integral *= _specResolution;
			norm += integral;
		}
	}
	return sqrt(norm);
}

double CReferenceStructure::findRotationAngle(CxObArray &matrix, int i, int j) {
// 	mprintf(RED, "findRotationAngle (%d, %d)\n", i, j);
	CxFloatArray *aii = (CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i];
	CxFloatArray *aij = (CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + j];
	CxFloatArray *ajj = (CxFloatArray *)matrix[(3*_atomCount-1)*j - j*(j-1)/2 + j];
	double int1 = 0.0, int2 = 0.0, int3 = 0.0;
	
	int k, l;
	for(k = 0; k < aij->GetSize(); k++) {
		int1 += (double)aij->GetAt(k) * (double)aij->GetAt(k);
		int2 += (double)aij->GetAt(k) * ((double)aii->GetAt(k) - (double)ajj->GetAt(k));
		int3 += ((double)aii->GetAt(k) - (double)ajj->GetAt(k)) * ((double)aii->GetAt(k) - (double)ajj->GetAt(k));
	}
	int1 *= (double)_specResolution;
	int2 *= (double)_specResolution;
	int3 *= (double)_specResolution;
	if(fabs(int2) < 1.0e-6) {
		double f = 4.0 * int1 - int3;
		if((fabs(int1) < 1.0) && (fabs(int2) < 1.0) && (fabs(int3) < 1.0))
			return 0.0;
		if(0.5 * f > 0.0)
			return 1.0;
		if(-2.0 * f > 0.0)
			return 0.0;
	} else {
		double f = (4.0 * int1 - int3) / int2;
		double r[2][2];
		double s[2][2];
		for(k = 0; k < 2; k++) {
			for(l = 0; l < 2; l++) {
				r[k][l] = (-f + (k == 0 ? -1.0 : 1.0) * sqrt(16.0 + f * f) + (l == 0 ? -1.0 : 1.0) * sqrt(2.0) * sqrt(16.0 + f * f - (k == 0 ? -1.0 : 1.0) * f * sqrt(16.0 + f * f))) / 4.0;
				s[k][l] = -2.0 * int2 / pow(r[k][l] * r[k][l] + 1.0, 4.0) * (2.0 * r[k][l] * (pow(r[k][l], 4.0) - 14.0 * r[k][l] * r[k][l] + 9.0) + f * (3.0 * pow(r[k][l], 4.0) - 8.0 * r[k][l] * r[k][l] + 1.0));
			}
		}
		for(k = 0; k < 2; k++) {
			for(l = 0; l < 2; l++) {
				if((fabs(r[k][l]) <= 1.0) && (s[k][l] > 0.0))
					return r[k][l];
			}
		}
	}
	
	mprintf(RED, "No suitable rotation angle found. Aborting.\n");
	abort();
}

void CReferenceStructure::calcIntegrals(CxObArray &matrix, CxFloatArray *integrals, CxFloatArray *centers) {
	integrals->SetSize(3*_atomCount);
	centers->SetSize(3*_atomCount);
	int i, k;
	for(i = 0; i < 3*_atomCount; i++) {
		integrals->GetAt(i) = 0.0f;
		centers->GetAt(i) = 0.0f;
		float integral_sqr = 0.0f;
		for(k = 0; k < ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetSize(); k++) {
// 			integrals->GetAt(i) += ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k) * ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k);
			integrals->GetAt(i) += ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k);
			integral_sqr += ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k) * ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k);
			centers->GetAt(i) += ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k) * ((CxFloatArray *)matrix[(3*_atomCount-1)*i - i*(i-1)/2 + i])->GetAt(k) * (float)k;
		}
		integrals->GetAt(i) *= _specResolution;// * 3.605675e-9;
		integral_sqr *= _specResolution;
		centers->GetAt(i) *= _specResolution * _specResolution / integral_sqr;
	}
}

CNormalCoordinateObservation::CNormalCoordinateObservation() {
//	char buf[BUF_SIZE];
//	char buf2[BUF_SIZE];
//	size_t remaining = BUF_SIZE;
	CxString buf, buf2, name;
	int i, j;

	if(g_oaMolecules.GetSize() > 1) {
/*#ifdef TARGET_LINUX
		remaining -= snprintf(buf, remaining, "    Which molecule should be observed (");
#else
		remaining -= sprintf(buf, "    Which molecule should be observed (");
#endif*/

		buf.sprintf("    Which molecule should be observed (");

		for(i = 0; i < g_oaMolecules.GetSize(); i++) {

/*			if(remaining < 1)
				break;
#ifdef TARGET_LINUX
			size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#else
			size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#endif
			strncat(buf, buf2, remaining - 1);
			remaining -= length;*/

			buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
			buf.strcat(buf2);

			if(i < g_oaMolecules.GetSize() - 1) {
/*#ifdef TARGET_LINUX
				length = snprintf(buf2, remaining, ", ");
#else
				length = sprintf(buf2, ", ");
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;*/

				buf2.sprintf(", ");
				buf.strcat(buf2);
			}
		}
//		strncat(buf, ")? ", remaining - 1);
		buf.strcat(")? ");

		m_iShowMol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
	} else {
		m_iShowMol = 0;
	}
	m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
	
// 	_calcIR = AskYesNo("    Calculate also transformed IR spectrum? [no] ", false);
	_calcIR = false;
	if(_calcIR) {
		g_bDipole = true;
		ParseDipole();
	}
	
	mprintf("\n");
	_refCount = 0;
	while(true) {
		_refCount++;
		mprintf(YELLOW, ">>> Reference Structure %d >>>\n\n", _refStructures.GetSize()+1);
		
/*		char name[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(name, BUF_SIZE, "normalcoordinate_%s_r%d", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, _refCount);
#else
		sprintf(name, "normalcoordinate_%s_r%d", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, _refCount);
#endif*/
		name.sprintf("normalcoordinate_%s_r%d", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, _refCount);

		CReferenceStructure *refStruct;
		try { refStruct = new CReferenceStructure(m_iShowMol, name, _calcIR); } catch(...) { refStruct = NULL; }
		if(refStruct == NULL) NewException((double)sizeof(CReferenceStructure), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		_refStructures.Add(refStruct);
		
		mprintf(YELLOW, "<<< End of Reference Structure %d <<<\n\n", _refStructures.GetSize());
		
		if(!AskYesNo("    Add another reference structure (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
// 	_refCount = AskUnsignedInteger("    How many reference structures do you wish to enter? [1] ", 1);
// 	for(int i = 0; i < _refCount; i++) {
// 		mprintf(YELLOW, "\n>>> Reference Structure %d >>>\n\n", _refStructures.GetSize()+1);
// 		
// 		CReferenceStructure *refStruct;
// 		try { refStruct = new CReferenceStructure(m_iShowMol, _calcIR); } catch(...) { refStruct = NULL; }
// 		if(refStruct == NULL) NewException((double)sizeof(CReferenceStructure), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		_refStructures.Add(refStruct);
// 		
// 		mprintf(YELLOW, "<<< End of Reference Structure %d <<<\n", _refStructures.GetSize());
// 	}
	
// 	if(g_bAdvanced2)
// 		_writeTransformedTrajectories = AskYesNo("\n    Save trajectories after transformation to reference frame for first observed molecule (y/n)? [no] ", false);
// 	else
	if(_refCount > 1) {
		_useInternals = AskYesNo("\n    Use mass-weighted Cartesians for probabilities (n) or define internal coordinates (y)? [no] ", false);
		if(_useInternals) {
			mprintf("\n    At the moment, only \"simple\" dihedral angles can be defined.\n");
			do {
				CxIntArray *internal;
				try { internal = new CxIntArray(); } catch(...) { internal = NULL; }
				if(internal == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				internal->SetSize(4);
				mprintf("\n");
				for(i = 0; i < 4; i++) {
//					char buf[BUF_SIZE];
					unsigned char parse[3];
					do
						AskString_ND("      Enter the %d. atom (e.g. C7): ", &buf, i+1);
					while(!ParseAtom(buf, m_iShowMol, parse[0], parse[1], parse[2]));
					internal->GetAt(i) = -1;
					for(j = 0; j < ((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_oaMolAtoms.GetSize(); j++) {
						CMolAtom *ma = (CMolAtom *)((CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]])->m_oaMolAtoms[j];
						if((ma->m_iType == parse[0]) && (ma->m_iNumber == parse[2])) {
							internal->GetAt(i) = j;
							break;
						}
					}
					if(internal->GetAt(i) == -1) {
						mprintf(RED, "Weird error.\n");
						abort();
					}
				}
				_internals.Add(internal);
			} while(AskYesNo("\n    Add another dihedral angle (y/n)? [no] ", false));
			_gaussianWidth = AskFloat("\n    Gaussian width for probability distribution (deg)? [10.0] ", 10.0f);
		} else {
			_gaussianWidth = AskFloat("    Gaussian width for probability distribution (pm*sqrt(amu))? [1.0] ", 1.0f);
		}
	} else {
		_useInternals = false;
		_gaussianWidth = 1000.0f;
	}
	
// 	_gaussianWidth = AskFloat("\n    Gaussian width for probability distribution (pm*sqrt(amu))? [1.0] ", 1.0f);
	
// 	if(g_iTrajSteps != -1) {
// 		_correlationDepth = 0.75 * g_iTrajSteps;
// 		if(_correlationDepth > 4096)
// 			_correlationDepth = 4096;
// 		_correlationDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the ACF (in time steps): [%d] ", _correlationDepth, _correlationDepth);
// 	} else {
// 		_correlationDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the ACF (in time steps): [256] ", 256);
// 	}
// 	int size = CalcFFTSize(_correlationDepth, false);
// 	if(_correlationDepth != size) {
// 		mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n", size, _correlationDepth);
// 		_correlationDepth = size;
// 	}
}

CNormalCoordinateObservation::~CNormalCoordinateObservation() {
	int i;
	for(i = 0; i < _distanceTimedev.GetSize(); i++)
		delete (CxFloatArray *)_distanceTimedev[i];
	for(i = 0; i < _refStructures.GetSize(); i++)
		delete (CReferenceStructure *)_refStructures[i];
}

void CNormalCoordinateObservation::initialize() {
	_numSteps = 0;
// 	if(_writeTransformedTrajectories) {
// 		int numPerm = 0;
// 		for(int i = 0; i < _refStructures.GetSize(); i++) {
// 			numPerm += ((CReferenceStructure *)_refStructures[i])->numPermutations();
// 		}
// 		try { _transformedTrajectoryFiles = new FILE *[numPerm]; } catch(...) { _transformedTrajectoryFiles = NULL; }
// 		if(_transformedTrajectoryFiles == NULL) NewException((double)numPerm * sizeof(FILE *), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		int index = 0;
// 		for(int i = 0; i < _refStructures.GetSize(); i++) {
// 			for(int j = 0; j < ((CReferenceStructure *)_refStructures[i])->numPermutations(); j++) {
// 				char name[BUF_SIZE];
// 				snprintf(name, BUF_SIZE, "reference_%s_%d_%d_transformed.xyz", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, i+1, j+1);
// 				_transformedTrajectoryFiles[index] = OpenFileWrite(name, false);
// 				if(_transformedTrajectoryFiles[index] == NULL) {
// 					mprintf(RED, "Could not open file %s: %s\n", name, strerror(errno));
// 				}
// 				index++;
// 			}
// 		}
// 	}
// 	for(int i = 0; i < m_iShowMolCount * _refStructures.GetSize(); i++) {
// 		CxObArray *a;
// 		try { a = new CxObArray(); } catch(...) { a = NULL; }
// 		if(a == NULL) NewException((double)sizeof(CxObArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		_distanceTimedev.Add(a);
// 	}
// 	for(int i = 0; i < m_iShowMolCount; i++) {
// 		for(int j = 0; j < _refStructures.GetSize(); j++) {
// 			for(int k = 0; k < 3; k++) {
// 				CxObArray *a;
// 				try { a = new CxObArray(); } catch(...) { a = NULL; }
// 				if(a == NULL) NewException((double)sizeof(CxObArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 				for(int l = 0; l < ((CReferenceStructure *)_refStructures[j])->numPermutations(); l++) {
// 					CxVec3Array *b;
// 					try { b = new CxVec3Array(); } catch(...) { b = NULL; }
// 					if(b == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 					a->Add(b);
// 				}
// 				_coordHistory.Add(a);
// 			}
// 		}
// 	}
// 	_historyIndex = 0;
// 	_calcVelocity = false;
// 	for(int i = 0; i < m_iShowMolCount * _refStructures.GetSize(); i++) {
// 		CxObArray *a;
// 		try { a = new CxObArray(); } catch(...) { a = NULL; }
// 		if(a == NULL) NewException((double)sizeof(CxObArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 		_velocityCache.Add(a);
// 	}
	
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	mprintf("    Global distance time development: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * n * _refCount * sizeof(float)));
	int i;
	for(i = 0; i < m_iShowMolCount * _refCount; i++) {
		CxFloatArray *a;
		try { a = new CxFloatArray(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetMaxSize(n);
		a->SetGrow((int)(0.1 * n));
		_distanceTimedev.Add(a);
	}

	for(i = 0; i < _refStructures.GetSize(); i++) {
		mprintf("  Initializing reference structure %d...\n", i+1);
		char name[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(name, BUF_SIZE, "normalcoordinate_%s_r%d", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, i+1);
#else
		sprintf(name, "normalcoordinate_%s_r%d", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, i+1);
#endif
		((CReferenceStructure *)_refStructures[i])->initialize(name);
	}
}

void CNormalCoordinateObservation::process(CTimeStep *ts) {
	int i, j;
	for(i = 0; i < m_iShowMolCount; i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
		CxVec3Array coord;
		coord.SetSize(sm->m_oaMolAtoms.GetSize());
		for(j = 0; j < sm->m_oaMolAtoms.GetSize(); j++) {
			coord[j] = ts->m_vaCoords[((CMolAtom *)sm->m_oaMolAtoms[j])->m_iOffset];
		}
		CxVector3 dipole(0.0f, 0.0f, 0.0f);
		if(_calcIR) {
			dipole = sm->m_vDipole;
		}
		
		for(j = 0; j < _refStructures.GetSize(); j++) {
			((CReferenceStructure *)_refStructures[j])->processCoordinates(coord, dipole, i);
			((CxFloatArray *)_distanceTimedev[i*_refCount + j])->Add(((CReferenceStructure *)_refStructures[j])->calcMinimumDistance(i, _useInternals, _internals));
			if(_numSteps >= 2)
				((CReferenceStructure *)_refStructures[j])->calcVelocities(i);
		}
	}
	for(i = 0; i < _refStructures.GetSize(); i++) {
		((CReferenceStructure *)_refStructures[i])->nextStep();
	}
	_numSteps++;
	
// 	for(int i = 0; i < m_iShowMolCount; i++) {
// 		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
// 		CxVec3Array coord;
// 		coord.SetSize(sm->m_oaMolAtoms.GetSize());
// 		for(int j = 0; j < sm->m_oaMolAtoms.GetSize(); j++) {
// 			coord[j] = ts->m_vaCoords[((CMolAtom *)sm->m_oaMolAtoms[j])->m_iOffset];
// 		}
// 		int index = 0;
// 		for(int j = 0; j < _refStructures.GetSize(); j++) {
// 			CxFloatArray *distances;
// 			try { distances = new CxFloatArray(); } catch(...) { distances = NULL; }
// 			if(distances == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			distances->SetSize(((CReferenceStructure *)_refStructures[j])->numPermutations());
// // 			CxObArray newCoords;
// // 			for(int k = 0; k < ((CReferenceStructure *)_refStructures[j])->numPermutations(); k++) {
// // 				CxVec3Array *a;
// // 				try { a = new CxVec3Array(); } catch(...) { a = NULL; }
// // 				if(a == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// // 				newCoords.Add(a);
// // 			}
// 			((CReferenceStructure *)_refStructures[j])->transformCoordinates(coord, *distances, *((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex]));
// 			if(_writeTransformedTrajectories && (i == 0)) {
// 				for(int k = 0; k < ((CReferenceStructure *)_refStructures[j])->numPermutations(); k++) {
// 					fprintf(_transformedTrajectoryFiles[index], "%6d\n\n", ((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex])->GetAt(k))->GetSize());
// 					for(int l = 0; l < ((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex])->GetAt(k))->GetSize(); l++) {
// 						fprintf(_transformedTrajectoryFiles[index], "%4s %14.8f %14.8f %14.8f\n", ((CAtom *)g_oaAtoms[sm->m_baAtomIndex[((CMolAtom *)sm->m_oaMolAtoms[l])->m_iType]])->m_sName, ((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex])->GetAt(k))->GetAt(l)[0] / 100.0f, ((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex])->GetAt(k))->GetAt(l)[1] / 100.0f, ((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex])->GetAt(k))->GetAt(l)[2] / 100.0f);
// 					}
// 					index++;
// 				}
// 			}
// 			((CxObArray *)_distanceTimedev[i*_refStructures.GetSize() + j])->Add(distances);
// 			if(_calcVelocity) {
// 				CxObArray *vel;
// 				try { vel = new CxObArray(); } catch(...) { vel = NULL; }
// 				if(vel == NULL) NewException((double)sizeof(CxObArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 				for(int k = 0; k < ((CReferenceStructure *)_refStructures[j])->numPermutations(); k++) {
// 					CxVec3Array *vec;
// 					try { vec = new CxVec3Array(); } catch(...) { vec = NULL; }
// 					if(vec == NULL) NewException((double)sizeof(CxVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 					vec->SetSize(((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex])->GetAt(k))->GetSize());
// 					int n = (_historyIndex + 1) % 3;
// 					for(int l = 0; l < ((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex])->GetAt(k))->GetSize(); l++) {
// 						vec->GetAt(l) = (((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+_historyIndex])->GetAt(k))->GetAt(l) - ((CxVec3Array *)((CxObArray *)_coordHistory[3*_refStructures.GetSize()*i+3*j+n])->GetAt(k))->GetAt(l)) / 2.0f / g_fTimestepLength * 1000.0f;
// 					}
// 					vel->Add(vec);
// 				}
// 				((CxObArray *)_velocityCache[i*_refStructures.GetSize() + j])->Add(vel);
// 			}
// // 			for(int k = 0; k < newCoords.GetSize(); k++)
// // 				delete (CxVec3Array *)newCoords[k];
// // 			newCoords.RemoveAll();
// 		}
// 	}
// 	_historyIndex = (_historyIndex + 1) % 3;
// 	if(_historyIndex == 2)
// 		_calcVelocity = true;
// 	_numSteps++;
}

void CNormalCoordinateObservation::finalize() {
	mprintf("    Saving global distance time development...\n");
	char name[BUF_SIZE];
#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "normalcoordinate_%s_global_dist_timedev.csv", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
#else
	sprintf(name, "normalcoordinate_%s_global_dist_timedev.csv", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
#endif
	FILE *distanceFile = OpenFileWrite(name, false);
	fprintf(distanceFile, "#Time (fs);");
	int i, j, k;
	for(i = 0; i < m_iShowMolCount; i++) {
		for(j = 0; j < _refCount; j++) {
			fprintf(distanceFile, " M%d R%d;", i+1, j+1);
		}
	}
	fprintf(distanceFile, "\n");
	for(i = 0; i < _numSteps; i++) {
		fprintf(distanceFile, "%.2f;", i * g_fTimestepLength * g_iStride);
		for(j = 0; j < m_iShowMolCount; j++) {
			for(k = 0; k < _refCount; k++) {
				fprintf(distanceFile, " %.8G;", ((CxFloatArray *)_distanceTimedev[j*_refCount + k])->GetAt(i));
			}
		}
		fprintf(distanceFile, "\n");
	}
	fclose(distanceFile);
	
	mprintf("    Saving global probability time development...\n");
#ifdef TARGET_LINUX
	snprintf(name, BUF_SIZE, "normalcoordinate_%s_global_prob_timedev.csv", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
#else
	sprintf(name, "normalcoordinate_%s_global_prob_timedev.csv", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
#endif
	FILE *probFile = OpenFileWrite(name, false);
	fprintf(probFile, "#Time (fs);");
	for(i = 0; i < m_iShowMolCount; i++) {
		for(j = 0; j < _refCount; j++) {
			fprintf(probFile, " M%d R%d;", i+1, j+1);
		}
	}
	fprintf(probFile, "\n");
	CxDoubleArray probSum;
	probSum.SetSize(m_iShowMolCount*_numSteps);
	for(i = 0; i < _numSteps; i++) {
		fprintf(probFile, "%.2f;", i * g_fTimestepLength * g_iStride);
		for(j = 0; j < m_iShowMolCount; j++) {
			CxDoubleArray probs;
			probs.SetSize(_refCount);
			probSum[j*_numSteps + i] = 0.0;
			if(_refCount == 1) {
				probs[0] = 1.0;
				probSum[j*_numSteps + i] = 1.0;
				((CReferenceStructure *)_refStructures[0])->setGlobalProbability(j, i, (float)(probs[0] / probSum[j*_numSteps + i]));
				fprintf(probFile, " %.8G;", probs[0] / probSum[j*_numSteps + i]);
			} else {
				for(k = 0; k < _refCount; k++) {
					probs[k] = exp(-(double)((CxFloatArray *)_distanceTimedev[j*_refCount + k])->GetAt(i) * (double)((CxFloatArray *)_distanceTimedev[j*_refCount + k])->GetAt(i) / 2.0 / (double)_gaussianWidth / (double)_gaussianWidth);
					probSum[j*_numSteps + i] += probs[k];
				}
				for(k = 0; k < _refCount; k++) {
					if(isnan(probs[k] / probSum[j*_numSteps + i])) {
						mprintf(RED, "Invalid probability for reference %d of molecule %d in step %d.\nMaybe it helps to increase the distribution width.\n", k+1, j+1, i+1);
						abort();
					}
					((CReferenceStructure *)_refStructures[k])->setGlobalProbability(j, i, (float)(probs[k] / probSum[j*_numSteps + i]));
					fprintf(probFile, " %.8G;", probs[k] / probSum[j*_numSteps + i]);
				}
			}
		}
		fprintf(probFile, "\n");
	}
	fclose(probFile);
	
// 	mprintf("    Calculating reference probabilities\n");
// 	char name[BUF_SIZE];
// 	snprintf(name, BUF_SIZE, "normalcoordinate_%s_global_prob_timedev.dat", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
// 	FILE *probFile = OpenFileWrite(name, false);
// 	fprintf(probFile, "# Time    ");
// 	for(int i = 0; i < m_iShowMolCount; i++) {
// 		for(int j = 0; j < _refStructures.GetSize(); j++) {
// 			fprintf(probFile, "  M. %3d,R. %3d", i+1, j+1);
// 		}
// 	}
// 	fprintf(probFile, "\n");
// 	for(int i = 0; i < _numSteps; i++) {
// 		fprintf(probFile, "%10.2f", i * g_fTimestepLength * g_iStride);
// 		for(int j = 0; j < m_iShowMolCount; j++) {
// 			CxFloatArray dist;
// 			dist.SetSize(_refStructures.GetSize());
// 			for(int k = 0; k < _refStructures.GetSize(); k++) {
// 				dist[k] = ((CReferenceStructure *)_refStructures[k])->getMinimumDistance(j, i);
// 			}
// // 			CxFloatArray probs;
// 			CxDoubleArray probs;
// 			probs.SetSize(_refStructures.GetSize());
// // 			float probSum = 0.0f;
// 			double probSum = 0.0;
// 			for(int k = 0; k < _refStructures.GetSize(); k++) {
// // 				probs[k] = expf(-dist[k] * dist[k] / 2.0f / _gaussianWidth / _gaussianWidth);
// 				probs[k] = exp(-(double)dist[k] * (double)dist[k] / 2.0 / (double)_gaussianWidth / (double)_gaussianWidth);
// 				probSum += probs[k];
// 			}
// 			for(int k = 0; k < _refStructures.GetSize(); k++) {
// 				((CReferenceStructure *)_refStructures[k])->setGlobalProbability(j, i, (float)(probs[k] / probSum));
// 				fprintf(probFile, " %14.8f", probs[k] / probSum);
// 			}
// 		}
// 		fprintf(probFile, "\n");
// 	}
// 	fclose(probFile);
	
	for(i = 0; i < _refStructures.GetSize(); i++) {
		mprintf("  Finalizing reference structure %d...\n", i+1);
		char name[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(name, BUF_SIZE, "normalcoordinate_%s_r%d", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, i+1);
#else
		sprintf(name, "normalcoordinate_%s_r%d", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, i+1);
#endif
		((CReferenceStructure *)_refStructures[i])->finalize(name);
	}
	
// 	mprintf("    Saving distance time development...\n");
// 	for(int i = 0; i < _refStructures.GetSize(); i++) {
// 		char name[BUF_SIZE];
// 		snprintf(name, BUF_SIZE, "distance_timedev_reference_%s_%d.dat", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, i+1);
// 		FILE *distanceFile = OpenFileWrite(name, false);
// 		for(int j = 0; j < _numSteps; j++) {
// 			fprintf(distanceFile, "%10.2f", j * g_fTimestepLength * g_iStride);
// 			for(int k = 0; k < m_iShowMolCount; k++) {
// 				for(int l = 0; l < ((CReferenceStructure *)_refStructures[i])->numPermutations(); l++) {
// 					fprintf(distanceFile, " %14.8f", ((CxFloatArray *)((CxObArray *)_distanceTimedev[k*_refStructures.GetSize() + i])->GetAt(j))->GetAt(l));
// 				}
// 			}
// 			fprintf(distanceFile, "\n");
// 		}
// 		fclose(distanceFile);
// 	}
// 	
// 	mprintf("    Saving probability time development...\n");
// 	for(int i = 0; i < _refStructures.GetSize(); i++) {
// 		char name[BUF_SIZE];
// 		snprintf(name, BUF_SIZE, "probability_timedev_reference_%s_%d.dat", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, i+1);
// 		FILE *probFile = OpenFileWrite(name, false);
// 		for(int j = 0; j < _numSteps; j++) {
// 			fprintf(probFile, "%10.2f", j * g_fTimestepLength * g_iStride);
// 			for(int k = 0; k < m_iShowMolCount; k++) {
// 				CxFloatArray probs;
// 				probs.SetSize(((CReferenceStructure *)_refStructures[i])->numPermutations());
// 				float probSum = 0.0f;
// 				for(int l = 0; l < ((CReferenceStructure *)_refStructures[i])->numPermutations(); l++) {
// 					probs[l] = expf(-((CxFloatArray *)((CxObArray *)_distanceTimedev[k*_refStructures.GetSize() + i])->GetAt(j))->GetAt(l) * ((CxFloatArray *)((CxObArray *)_distanceTimedev[k*_refStructures.GetSize() + i])->GetAt(j))->GetAt(l) / 2.0f / _gaussianWidth);
// 					probSum += probs[l];
// 				}
// 				for(int l = 0; l < ((CReferenceStructure *)_refStructures[i])->numPermutations(); l++) {
// 					fprintf(probFile, " %14.8f", probs[l] / probSum);
// 				}
// 			}
// 			fprintf(probFile, "\n");
// 		}
// 		fclose(probFile);
// 	}
// 	
// 	CCrossCorrelation *cc;
// 	try { cc = new CCrossCorrelation(); } catch(...) { cc = NULL; }
// 	if(cc == NULL) NewException((double)sizeof(CCrossCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	cc->Init(_numSteps - 2, _correlationDepth, true);
// 	
// 	for(int i = 0; i < _refStructures.GetSize(); i++) {
// 		mprintf("    Calculating cross correlation matrix for reference structure %d...\n", i+1);
// 		CxObArray ccMatrix;
// 		for(int j = 0; j < 3 * m_iGesAtomCount; j++) {
// 			for(int k = 0; k <= j; k++) {
// 				CxFloatArray *a;
// 				try { a = new CxFloatArray(); } catch(...) { a = NULL; }
// 				if(a == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 				for(int l = 0; l < m_iShowMolCount; l++) {
// 					for(int m = 0; m < ((CReferenceStructure *)_refStructures[i])->numPermutations(); m++) {
// 						
// 					}
// 				}
// 				ccMatrix.Add(a);
// 			}
// 		}
// 	}
// 	
// 	delete cc;
// 	
// 	if(_writeTransformedTrajectories) {
// 		int index = 0;
// 		for(int i = 0; i < _refStructures.GetSize(); i++) {
// 			for(int j = 0; j < ((CReferenceStructure *)_refStructures[i])->numPermutations(); j++) {
// 				fclose(_transformedTrajectoryFiles[index++]);
// 			}
// 		}
// 		delete[] _transformedTrajectoryFiles;
// 	}
// 	for(int i = 0; i < _distanceTimedev.GetSize(); i++) {
// 		for(int j = 0; j < ((CxObArray *)_distanceTimedev[i])->GetSize(); j++) {
// 			delete (CxFloatArray *)((CxObArray *)_distanceTimedev[i])->GetAt(j);
// 		}
// 		delete (CxObArray *)_distanceTimedev[i];
// 	}
// 	for(int i = 0; i < _coordHistory.GetSize(); i++) {
// 		for(int j = 0; j < ((CxObArray *)_coordHistory[i])->GetSize(); j++) {
// 			delete (CxVec3Array *)((CxObArray *)_coordHistory[i])->GetAt(j);
// 		}
// 		delete (CxObArray *)_coordHistory[i];
// 	}
// 	for(int i = 0; i < _velocityCache.GetSize(); i++) {
// 		for(int j = 0; j < ((CxObArray *)_velocityCache[i])->GetSize(); j++) {
// 			for(int k = 0; k < ((CxObArray *)((CxObArray *)_velocityCache[i])->GetAt(j))->GetSize(); k++) {
// 				delete (CxVec3Array *)((CxObArray *)((CxObArray *)_velocityCache[i])->GetAt(j))->GetAt(k);
// 			}
// 			delete (CxObArray *)((CxObArray *)_velocityCache[i])->GetAt(j);
// 		}
// 		delete (CxObArray *)_velocityCache[i];
// 	}
}

bool gatherNormalCoordinate() {
	while(true) {
		mprintf(YELLOW, ">>> Normal Coordinate Observation %d >>>\n\n", g_normalCoordinateObservations.GetSize()+1);
		
		CNormalCoordinateObservation *obs;
		try { obs = new CNormalCoordinateObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CNormalCoordinateObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_normalCoordinateObservations.Add(obs);
		
		mprintf(YELLOW, "\n<<< End of Normal Coordinate Observation %d <<<\n\n", g_normalCoordinateObservations.GetSize());
		
		if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	mprintf("\n");
	
	return true;
}

bool initializeNormalCoordinate() {
	int i;
	for(i = 0; i < g_normalCoordinateObservations.GetSize(); i++) {
		mprintf("Initializing Normal Coordinate Observation %d...\n", i+1);
		((CNormalCoordinateObservation *)g_normalCoordinateObservations[i])->initialize();
	}
	
	return true;
}

void processNormalCoordinate(CTimeStep *ts) {
	int i;
	for(i = 0; i < g_normalCoordinateObservations.GetSize(); i++)
		((CNormalCoordinateObservation *)g_normalCoordinateObservations[i])->process(ts);
}

void finalizeNormalCoordinate() {
	int i;
	for(i = 0; i < g_normalCoordinateObservations.GetSize(); i++) {
		mprintf(YELLOW, ">>> Normal Coordinate Observation %d >>>\n\n", i+1);
		((CNormalCoordinateObservation *)g_normalCoordinateObservations[i])->finalize();
		mprintf(YELLOW, "\n<<< End of Normal Coordinate Observation %d <<<\n", i+1);
	}
	for(i = 0; i < g_normalCoordinateObservations.GetSize(); i++)
		delete (CNormalCoordinateObservation *)g_normalCoordinateObservations[i];
}

CEckartReferenceStructure::CEckartReferenceStructure() {
	CxString buf;

	try { _singleMol = new CSingleMolecule(); } catch(...) { _singleMol = NULL; }
	if(_singleMol == NULL) NewException((double)sizeof(CSingleMolecule), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	try { _refTimestep = new CTimeStep(); } catch(...) { _refTimestep = NULL; }
	if(_refTimestep == NULL) NewException((double)sizeof(CTimeStep), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	try { _mapAtoms = new CAtomGroup(); } catch(...) { _mapAtoms = NULL; }
	if(_mapAtoms == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	while(true) {
//		char buf[BUF_SIZE];
		AskString_ND("    Filename for reference structure: ", &buf);
		FILE *refFile = fopen(buf, "r");
		if(refFile == NULL) {
			mprintf(RED, "Could not open reference structure file: %s\n", strerror(errno));
			continue;
		}
		if(!_refTimestep->ReadTimestep(refFile, true)) {
			mprintf(RED, "Error reading reference structure file\n");
			continue;
		}
		strncpy(_filename, buf, 128);
		_filename[127] = 0;
		mprintf("\n    Starting molecule recoginition...\n");
		if(!recognizeMolecule(_refTimestep)) {
			mprintf(RED, "Molecule recoginition failed\n\n");
			continue;
		}
		break;
	}
	
	if(((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_laSingleMolIndex.GetSize() > 1) {
		_singleMolIndex = AskRangeInteger("    Which representant of %s in the trajectory should be mapped to the reference structure (1-%d)? ", 1, ((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_laSingleMolIndex.GetSize(), 1, ((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_sName, ((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_laSingleMolIndex.GetSize()) - 1;
	} else {
		_singleMolIndex = 0;
		mprintf("    Mapping molecule %s in trajectory to reference structure.\n", ((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_sName);
	}
	mprintf("\n");
	
	while(true) {
		mprintf("    Which atom(s) to map (e.g. \"C1,C3-5,H\", \"*\"=all)? [*] ");
		inpprintf("! Which atom(s) to map (e.g. \"C1,C3-5,H\", \"*\"=all)? [*]\n");
//		char buf[BUF_SIZE];
		myget(&buf);
		if(strlen(buf) == 0) {
			if(!_mapAtoms->ParseAtoms((CMolecule *)g_oaMolecules[_singleMol->m_iMolType], "*")) {
				eprintf("Weird error.\n");
				continue;
			}
		} else if(!_mapAtoms->ParseAtoms((CMolecule *)g_oaMolecules[_singleMol->m_iMolType], buf)) {
			continue;
		}
		break;
	}
	
	mprintf("\n    Mapping %d atoms.\n", _mapAtoms->m_iAtomGes);
	
	_refCentroidCoord.SetSize(_mapAtoms->m_iAtomGes);
	_refCentroid = CxVector3(0.0f, 0.0f, 0.0f);
	_weightSum = 0.0f;
	_weights.SetSize(_mapAtoms->m_iAtomGes);
	int i;
	int k = 0;
	for(i = 0; i < _mapAtoms->m_baAtomType.GetSize(); i++) {
		CxIntArray *a = (CxIntArray *)_mapAtoms->m_oaAtoms[i];
		int j;
		for(j = 0; j < a->GetSize(); j++) {
			_refCentroidCoord[k] = (_refTimestep->m_vaCoords[((CxIntArray *)_singleMol->m_oaAtomOffset[_mapAtoms->m_baAtomType[i]])->GetAt(a->GetAt(j))]);
			_weights[k] = ((CAtom *)g_oaAtoms[_mapAtoms->m_baRealAtomType[i]])->m_pElement->m_fMass;
			_weightSum += _weights[k];
			_refCentroid += _refCentroidCoord[k] * _weights[k];
			k++;
		}
	}
	_refCentroid /= _weightSum;
	for(i = 0; i < _refCentroidCoord.GetSize(); i++) {
		_refCentroidCoord[i] -= _refCentroid;
	}
}

CEckartReferenceStructure::~CEckartReferenceStructure() {
	delete _singleMol;
	delete _refTimestep;
	delete _mapAtoms;
}

void CEckartReferenceStructure::process(CTimeStep *ts) {
	CxVec3Array coord;
	coord.SetSize(_mapAtoms->m_iAtomGes);
	CxVector3 centroid(0.0f, 0.0f, 0.0f);
	CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_laSingleMolIndex[_singleMolIndex]];
	int i, j;
	int k = 0, l;
	for(i = 0; i < _mapAtoms->m_baAtomType.GetSize(); i++) {
		CxIntArray *a = (CxIntArray *)_mapAtoms->m_oaAtoms[i];
		int j;
		for(j = 0; j < a->GetSize(); j++) {
			coord[k] = (ts->m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[_mapAtoms->m_baAtomType[i]])->GetAt(a->GetAt(j))]);
			centroid += coord[k] * _weights[k];
			k++;
		}
	}
	centroid /= _weightSum;
	
	float covarianceMatrix[9] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
	for(i = 0; i < coord.GetSize(); i++) {
		CxVector3 v1 = coord[i] - centroid;
		CxVector3 v2 = _refCentroidCoord[i];
		int k, l;
		for(k = 0; k < 3; k++) {
			for(l = 0; l < 3; l++) {
				covarianceMatrix[k*3+l] += v1[k] * v2[l] * _weights[i];
			}
		}
	}
	
//************** NEW ********************
// 	float sValues[3];
// 	float vMatrix[9];
// 	ComputeSVD_Flat(covarianceMatrix, 3, 3, sValues, vMatrix);
// 	
// 	for(int j = 0; j < 3; j++) {
// 		for(int k = 0; k < 3; k++) {
// 			_rotationMatrix[3*j+k] = 0.0f;
// 			for(int l = 0; l < 3; l++) {
// 				_rotationMatrix[3*j+k] += vMatrix[3*j+l] * covarianceMatrix[3*k+l];
// 			}
// 		}
// 	}
// 	
// 	float determinant = _rotationMatrix[0]*_rotationMatrix[4]*_rotationMatrix[8] + _rotationMatrix[1]*_rotationMatrix[5]*_rotationMatrix[6] + _rotationMatrix[2]*_rotationMatrix[3]*_rotationMatrix[7] - _rotationMatrix[2]*_rotationMatrix[4]*_rotationMatrix[6] - _rotationMatrix[1]*_rotationMatrix[3]*_rotationMatrix[8] - _rotationMatrix[0]*_rotationMatrix[5]*_rotationMatrix[7];
// 	if(determinant < 0.0f) {
// 		for(int j = 0; j < 3; j++) {
// 			for(int k = 0; k < 3; k++) {
// 				_rotationMatrix[3*j+k] = 0.0f;
// 				for(int l = 0; l < 3; l++) {
// 					_rotationMatrix[3*j+k] += vMatrix[3*j+l] * covarianceMatrix[3*k+l] * (l == 0 ? -1.0f : 1.0f);
// 				}
// 			}
// 		}
// 	}
//***************************************
//************* NEW 2 *******************
	float sValues[3];
	float uMatrix[9];
	float vMatrix[9];
	SVD_3x3(covarianceMatrix, sValues, uMatrix, vMatrix);
	
	for(j = 0; j < 3; j++) {
		for(k = 0; k < 3; k++) {
			_rotationMatrix[3*j+k] = 0.0f;
			for(l = 0; l < 3; l++) {
				_rotationMatrix[3*j+k] += vMatrix[3*j+l] * uMatrix[3*k+l];
			}
		}
	}

	float determinant = _rotationMatrix[0]*_rotationMatrix[4]*_rotationMatrix[8] + _rotationMatrix[1]*_rotationMatrix[5]*_rotationMatrix[6] + _rotationMatrix[2]*_rotationMatrix[3]*_rotationMatrix[7] - _rotationMatrix[2]*_rotationMatrix[4]*_rotationMatrix[6] - _rotationMatrix[1]*_rotationMatrix[3]*_rotationMatrix[8] - _rotationMatrix[0]*_rotationMatrix[5]*_rotationMatrix[7];
	if(determinant < 0.0f) {
		for(j = 0; j < 3; j++) {
			for(k = 0; k < 3; k++) {
				_rotationMatrix[3*j+k] = 0.0f;
				for(l = 0; l < 3; l++) {
					_rotationMatrix[3*j+k] += vMatrix[3*j+l] * uMatrix[3*k+l] * (l == 2 ? -1.0f : 1.0f);
				}
			}
		}
	}
//***************************************
// 	float sValues[3];
// 	float uMatrix[9];
// 	float vtMatrix[9];
// 	singularValueDecompose_3x3(covarianceMatrix, sValues, uMatrix, vtMatrix);
// 	
// 	if(vtMatrix[0]*vtMatrix[4]*vtMatrix[8] + vtMatrix[1]*vtMatrix[5]*vtMatrix[6] + vtMatrix[2]*vtMatrix[3]*vtMatrix[7] - vtMatrix[2]*vtMatrix[4]*vtMatrix[6] - vtMatrix[1]*vtMatrix[3]*vtMatrix[8] - vtMatrix[0]*vtMatrix[5]*vtMatrix[7] < 0) {
// 		int j;
// 		for(j = 0; j < 3; j++) {
// 			vtMatrix[6+j] *= -1.0f;
// 		}
// 	}
// 	
// 	int j;
// 	for(j = 0; j < 3; j++) {
// 		int k;
// 		for(k = 0; k < 3; k++) {
// 			_rotationMatrix[3*k+j] = 0.0f;
// 			int l;
// 			for(l = 0; l < 3; l++) {
// 				_rotationMatrix[3*k+j] += uMatrix[3*j+l] * vtMatrix[3*l+k];
// 			}
// 		}
// 	}
	
	for(j = 0; j < 3; j++) {
		_translationVector[j] = _refCentroid[j];
		int k;
		for(k = 0; k < 3; k++) {
			_translationVector[j] -= _rotationMatrix[3*j+k] * centroid[k];
		}
	}
}

bool CEckartReferenceStructure::recognizeMolecule(CTimeStep *ts) {
	bool ok = true;
	unsigned int *stack;
	int i, j;

	try { stack = new unsigned int[ts->m_iGesAtomCount]; } catch(...) { stack = NULL; }
	if(stack == NULL) NewException((double)ts->m_iGesAtomCount*sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	bool *used;
	try { used = new bool[ts->m_iGesAtomCount]; } catch(...) { used = NULL; }
	if(used == NULL) NewException((double)ts->m_iGesAtomCount*sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	for(i = 0; (unsigned int)i < ts->m_iGesAtomCount; i++) {
		used[i] = false;
	}
	
	CxByteArray baAtomIndex;
	baAtomIndex.SetSize(ts->m_iGesAtomCount);
	for(i = 0; (unsigned int)i < ts->m_iGesAtomCount; i++) {
//		char buf[64];
		CxString buf;
//		strncpy(buf, (char *)(ts->m_paLabels[i]), 64);
//		buf[63] = 0;

		buf.strcpy((char *)(ts->m_paLabels[i]));
		ReplaceDigits(&buf);

		baAtomIndex[i] = 255;
		for(j = 0; j < g_oaAtoms.GetSize(); j++) {
			if(mystricmp(buf, ((CAtom *)g_oaAtoms[j])->m_sName) == 0) {
				CAtom *a = (CAtom *)g_oaAtoms[j];
				while(a->m_pMergedTo != NULL)
					a = a->m_pMergedTo;
				baAtomIndex[i] = a->m_iIndex;
				break;
			}
		}
		if(baAtomIndex[i] == 255) {
			mprintf(RED, "Atom type \"%s\" not known\n", (const char*)buf);
			ok = false;
		}
	}
	
	if(ok) {
		recognizeMoleculeRecursion(0, used, 0, stack, ts, baAtomIndex);
		for(i = 0; (unsigned int)i < ts->m_iGesAtomCount; i++) {
			if(!used[i]) {
				mprintf(RED, "Some atoms are not connected to the molecule\n");
				ok = false;
				break;
			}
		}
	}
	
	if(ok) {
		mprintf("    Recognized one molecule\n");
		
		mprintf("    Sorting atom types...\n");
		for(i = 0; i < _singleMol->m_baAtomIndex.GetSize() - 1; i++) {
			int a = -1;
			int b = 999;
			for(j = i; j < _singleMol->m_baAtomIndex.GetSize(); j++) {
				if(_singleMol->m_baAtomIndex[j] < b) {
					b = _singleMol->m_baAtomIndex[j];
					a = j;
				}
			}
			if((a != -1) && (a != i)) {
				b = _singleMol->m_baAtomIndex[a];
				_singleMol->m_baAtomIndex[a] = _singleMol->m_baAtomIndex[i];
				_singleMol->m_baAtomIndex[i] = b;
				for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
					if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType == i)
						((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType = a;
					else if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType == a)
						((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType = i;
				}
			}
		}
		
		mprintf("    Setting up bond list...\n");
		for(i = 0; i < _singleMol->m_laBonds.GetSize() / 2; i++) {
			int a = _singleMol->m_laBonds[2*i];
			int b = _singleMol->m_laBonds[2*i+1];
			int c = -1;
			for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
				if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iOffset == a) {
					c = j;
					break;
				}
			}
			if(c == -1) {
				mprintf(RED, "Weird error.\n");
				abort();
			}
			int d = -1;
			for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
				if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iOffset == b) {
					d = j;
					break;
				}
			}
			if(d == -1) {
				mprintf(RED, "Weird error.\n");
				abort();
			}
			((CMolAtom *)_singleMol->m_oaMolAtoms[c])->m_oaBonds.Add((CMolAtom *)_singleMol->m_oaMolAtoms[d]);
			((CMolAtom *)_singleMol->m_oaMolAtoms[d])->m_oaBonds.Add((CMolAtom *)_singleMol->m_oaMolAtoms[c]);
		}
		
		mprintf("    Building atom codes...\n");
		_singleMol->BuildAtomCodes();
		
		mprintf("    Creating topological atom order...\n");
		for(i = 0; i < _singleMol->m_baAtomIndex.GetSize(); i++) {
			CxIntArray *a;
			try { a = new CxIntArray(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			_singleMol->m_oaAtomOffset.Add(a);
			for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
				if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iType == i) {
					((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iNumber = a->GetSize();
					a->Add(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_iOffset);
				}
			}
		}
		
		mprintf("    Comparing to global molecule types...\n");
		int molType = -1;
		for(i = 0; i < g_oaMolecules.GetSize(); i++) {
			CSingleMolecule *compMol = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[i])->m_laSingleMolIndex[0]];
			if(_singleMol->m_oaMolAtoms.GetSize() != compMol->m_oaMolAtoms.GetSize())
				continue;
			bool diff = false;
			for(j = 0; j < _singleMol->m_oaMolAtoms.GetSize(); j++) {
				if(((CMolAtom *)_singleMol->m_oaMolAtoms[j])->m_fAtomCode != ((CMolAtom *)compMol->m_oaMolAtoms[j])->m_fAtomCode) {
					diff = true;
					break;
				}
			}
			if(!diff) {
				molType = i;
				break;
			}
		}
		if(molType == -1) {
			ok = false;
			mprintf(RED, "The reference structure does not belong to any molecule type\n");
		} else {
			_singleMol->m_iMolType = molType;
			mprintf("    The reference structure belongs to molecule type %d (%s)\n", _singleMol->m_iMolType+1, ((CMolecule *)g_oaMolecules[_singleMol->m_iMolType])->m_sName);
// 			_atomCount = _singleMol->m_oaMolAtoms.GetSize();
		}
	}
	mprintf("\n");
	
	delete[] used;
	delete[] stack;
	
	return ok;
}

void CEckartReferenceStructure::recognizeMoleculeRecursion(unsigned int index, bool *used, int depth, unsigned int *stack, CTimeStep *ts, CxByteArray &baAtomIndex) {
	int i, j;
	if(g_bVerbose) {
		mprintf("    ");
		for(i = 1; i < depth; i++)
			mprintf("    ");
		if(depth > 0)
			mprintf(WHITE, "\\--");
		mprintf(CYAN, "%s", ((CAtom *)g_oaAtoms[baAtomIndex[index]])->m_sName);
		mprintf("(%d)", index+1);
	}
	
	stack[depth] = index;
	used[index] = true;
	
	bool found = false;
	for(i = 0; i < _singleMol->m_baAtomIndex.GetSize(); i++) {
		if(_singleMol->m_baAtomIndex[i] == baAtomIndex[index]) {
			CMolAtom *molAtom;
			try { molAtom = new CMolAtom(); } catch(...) { molAtom = NULL; }
			if(molAtom == NULL) NewException((double)sizeof(CMolAtom), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			molAtom->m_iOffset = index;
			molAtom->m_iType = i;
			_singleMol->m_oaMolAtoms.Add(molAtom);
			found = true;
			break;
		}
	}
	if(!found) {
		CMolAtom *molAtom;
		try { molAtom = new CMolAtom(); } catch(...) { molAtom = NULL; }
		if(molAtom == NULL) NewException((double)sizeof(CMolAtom), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		molAtom->m_iOffset = index;
		molAtom->m_iType = _singleMol->m_baAtomIndex.GetSize();
		_singleMol->m_oaMolAtoms.Add(molAtom);
		_singleMol->m_baAtomIndex.Add(baAtomIndex[index]);
	}
	
	int numNeighbors = 0;
	int neighborList[MAX_BONDS];
	for(i = 0; (unsigned int)i < ts->m_iGesAtomCount; i++) {
		if((unsigned int)i == index)
			continue;
		if(((CAtom *)g_oaAtoms[baAtomIndex[index]])->m_bExclude)
			continue;
		float distance;
		if(recognizeMoleculeBondRange(ts, index, i, &distance, baAtomIndex)) {
			if(used[i]) {
				if((depth > 0) && ((unsigned int)i != stack[depth-1])) {
					if(_singleMol->m_oaRings.GetSize() < 100) {
						if(g_bVerbose) {
							mprintf(GREEN, " <-- Ring closure: ");
							mprintf("%s(%d)", ((CAtom *)g_oaAtoms[baAtomIndex[stack[depth]]])->m_sName, stack[depth]+1);
						}
						CxIntArray *ring;
						try { ring = new CxIntArray(); } catch(...) { ring = NULL; }
						if(ring == NULL) NewException((double)sizeof(CxIntArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
						_singleMol->m_oaRings.Add(ring);
						ring->Add(stack[depth]);
						for(j = depth - 1; (stack[j] != (unsigned int)i) && (j >= 0); j--) {
							ring->Add(stack[j]);
							if(g_bVerbose)
								mprintf(" - %s(%d)", ((CAtom *)g_oaAtoms[baAtomIndex[stack[j]]])->m_sName, stack[j]+1);
						}
						ring->Add(i);
						if(g_bVerbose)
							mprintf(" - %s(%d)", ((CAtom *)g_oaAtoms[baAtomIndex[i]])->m_sName, i+1);
					} else {
						mprintf(RED, "More than 100 rings\n");
					}
				}
				continue;
			}
			if(distance < 50.0f) {
				mprintf(RED, "Atoms %s(%d) and %s(%d) are very close to each other: %.4f pm\n", ((CAtom *)g_oaAtoms[baAtomIndex[index]])->m_sName, index + 1, ((CAtom *)g_oaAtoms[baAtomIndex[i]])->m_sName, i + 1, distance);
			}
			neighborList[numNeighbors++] = i;
			if(numNeighbors >= MAX_BONDS) {
				mprintf(RED, "Atom %s(%d) has more than %d bonds\n", ((CAtom *)g_oaAtoms[baAtomIndex[index]])->m_sName, index + 1, MAX_BONDS);
				break;
			}
		}
	}
	
	if(g_bVerbose)
		mprintf("\n");
	
	for(i = 0; i < numNeighbors; i++) {
		_singleMol->m_laBonds.Add(index);
		_singleMol->m_laBonds.Add(neighborList[i]);
		if(!used[neighborList[i]]) {
			recognizeMoleculeRecursion(neighborList[i], used, depth+1, stack, ts, baAtomIndex);
		}
	}
}

bool CEckartReferenceStructure::recognizeMoleculeBondRange(CTimeStep *ts, int i1, int i2, float *distance, CxByteArray &baAtomIndex) {
	float x = ts->m_vaCoords[i1][0] - ts->m_vaCoords[i2][0];
	float y = ts->m_vaCoords[i1][1] - ts->m_vaCoords[i2][1];
	float z = ts->m_vaCoords[i1][2] - ts->m_vaCoords[i2][2];
	
	if(g_bPeriodic)
	{
		if(g_bPeriodicX)
		{
			while(x < -g_fBoxX/2)
				x += g_fBoxX;
			while(x > g_fBoxX/2)
				x -= g_fBoxX;
		}
		if(g_bPeriodicY)
		{
			while(y < -g_fBoxY/2)
				y += g_fBoxY;
			while(y > g_fBoxY/2)
				y -= g_fBoxY;
		}
		if(g_bPeriodicZ)
		{
			while(z < -g_fBoxZ/2)
				z += g_fBoxZ;
			while(z > g_fBoxZ/2)
				z -= g_fBoxZ;
		}
	}
	
	*distance = sqrtf(x*x + y*y + z*z);
	
	return *distance < (((CAtom *)g_oaAtoms[baAtomIndex[i1]])->m_pElement->m_fRadius + ((CAtom *)g_oaAtoms[baAtomIndex[i2]])->m_pElement->m_fRadius) * g_fBondFactor;
}

static CEckartReferenceStructure *g_eckartReferenceStructure;
static FILE *g_eckartTransformFile;

bool gatherEckartTransform() {
	try { g_eckartReferenceStructure = new CEckartReferenceStructure(); } catch(...) { g_eckartReferenceStructure = NULL; }
	if(g_eckartReferenceStructure == NULL) NewException((double)sizeof(CEckartReferenceStructure), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	if(g_bXYZ4thCol)
		g_bReadChargesFrom4thXYZ = true;
	
	return true;
}

bool initializeEckartTransform() {
	char filename[BUF_SIZE];
#ifdef TARGET_LINUX
	char buf[BUF_SIZE];
	strncpy(buf, g_sInputTraj, BUF_SIZE);
	buf[BUF_SIZE-1] = 0;
	char *p = strrchr(buf, '/');
	if(p == NULL)
		p = buf;
	else
		p++;
	size_t s = strcspn(p, ".");
	if(s > BUF_SIZE - 8)
		s = BUF_SIZE - 8;
	strncpy(filename, p, s);
	filename[s] = 0;
	strcat(filename, "_out.xyz");
#else
	sprintf(filename, "out.xyz");
#endif
	mprintf("    Saving processed trajectory as %s\n\n", filename);
	g_eckartTransformFile = OpenFileWrite(filename, false);
	
	return true;
}

void processEckartTransform(CTimeStep *ts) {
	g_eckartReferenceStructure->process(ts);
	
	fprintf(g_eckartTransformFile, "%d\n", g_iGesAtomCount);
	if(ts->m_pComment != NULL)
		fprintf(g_eckartTransformFile, "%s\n", ts->m_pComment);
	else
		fprintf(g_eckartTransformFile, "\n");
	
	int i;
	for(i = 0; i < g_iGesAtomCount; i++) {
		CxVector3 newCoord = g_eckartReferenceStructure->translationVector();
		int j;
		for(j = 0; j < 3; j++) {
			int k;
			for(k = 0; k < 3; k++) {
				newCoord[j] += g_eckartReferenceStructure->rotationMatrix()[3*j+k] * ts->m_vaCoords[i][k];
			}
		}
		if(g_bReadChargesFrom4thXYZ) {
			fprintf(g_eckartTransformFile, "%2s %10.5f %10.5f %10.5f %10.5f\n", g_waAtomRealElement[i] == 60000 ? "X" : ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName, newCoord[0] / 100.0f, newCoord[1] / 100.0f, newCoord[2] / 100.0f, ts->m_faCharge[i]);
		} else {
			fprintf(g_eckartTransformFile, "%2s %10.5f %10.5f %10.5f\n", g_waAtomRealElement[i] == 60000 ? "X" : ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName, newCoord[0] / 100.0f, newCoord[1] / 100.0f, newCoord[2] / 100.0f);
		}
	}
}

void finalizeEckartTransform() {
	delete g_eckartReferenceStructure;
	fclose(g_eckartTransformFile);
}
