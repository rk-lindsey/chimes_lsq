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

#ifndef VOROWRAPPER_H
#define VOROWRAPPER_H

#include "tools.h"
#include "xobject.h"
#include "voro++.h"
#include "xobarray.h"
#include "df.h"
#include "xmatrix3.h"
#include "xwordarray.h"
#include "xfloatarray.h"
#include "xintarray.h"
#include "xstring.h"
#include "xlongarray.h"


class CTimeStep;
class CVoroWrapper;


class CVoroAtom : public CxObject
{
public:
	void Dump(const char *s);
	void InitAnalyses();
	CVoroAtom();
	~CVoroAtom();

	CVoroWrapper *m_pParent;

	int m_iMolecule;
	int m_iAtomType;
	int m_iRealAtomType;
	int m_iAtom;

	int *m_pNbhTempMol;

	CDF *m_pVolume;
	CDF *m_pSurfaceArea;
	CDF *m_pExposedSurface;
	CDF *m_pExposedSurfacePerc;
	CDF *m_pAVRatio;
	CDF *m_pFaces;
	CDF *m_pFaceOrders;
	CDF *m_pMaxRadius;
	CDF **m_pNbhAtoms;
	CDF **m_pNbhMolecules;
	CDF **m_pNbhDistAtoms;
	CDF **m_pNbhDistMolecules;
};


class CVoroMolecule : public CxObject
{
public:
	void Add(CVoroMolecule *m);
	void Dump(const char *s);
	void InitAnalyses();
	CVoroMolecule();
	~CVoroMolecule();

	CVoroWrapper *m_pParent;

	int m_iMolecule;
	int m_iSingleMol;
	int m_iCenterOffset;

	int *m_pNbhTempMol;
	int *m_pNbhTempAtoms;

	double tempvol;
	double tempsurf;
	int tempfaces;
	double *m_pNbhTempArea;

	CDF *m_pVolume;
	CDF *m_pSurfaceArea;
	CDF *m_pAVRatio;
	CDF *m_pFaces;
	CDF **m_pNbhAtoms;
	CDF **m_pNbhMolecules;
	CDF **m_pNbhDistAtoms;
	CDF **m_pNbhDistMolecules;
	CDF **m_pNbhAreaMolecules;
};




class CVoroWrapper : public CxObject
{
public:
	bool m_bIncludeWritten;
	bool m_bWritePOVMovie;
	CxVector3 m_vPOVFaceBleach;
	float m_fPOVFaceBleach;
	float m_fPOVNbHighlightFac;
	float m_fPOVRayTrans;
	float m_fPOVNbTrans;
	float m_fPOVFaceTrans;
	float m_fPOVEdgeTrans;
	CxFloatArray m_faVoroMetric;
	bool m_bVoroMetricCDF;
	bool m_bVoroMetric;
	C2DF *m_pVoroMetricCDF2;
	C2DF *m_pVoroMetricCDF;
	int m_iPOVThreads;
	int m_iPOVThreadPos;
	FILE *m_fPOVCT;
	double m_fPOVAtomCenter;
	double m_fPOVScale;
	bool m_bPOVMolecules;
	double m_fPOVExplode;
	double m_fPOVClip;
	CxDVector3 m_vPOVBoxClipLow;
	CxDVector3 m_vPOVBoxClipHigh;
	int m_iPOVBoxTotalSteps;
	int m_iPOVBoxTotalFrames;
	int m_iPOVCurrentChoreographyPart;


	double m_fPOV_FPS;
	CxObArray m_oaBoxChoreography;
	double m_fPOVBoxOuterOpacity;
	double m_fPOVZoom;
	int m_iPOVFrameCounter;
	char m_sPOVText[256];
	double m_fPOVAngle;
	CxLongArray m_iaPOVBoxAtomColors;
	FILE **m_fPOVBoxScript;
	bool m_bPOVBox;
	void FinishSurfCover(const char *s);
	void ProcessSurfCover(CTimeStep *ts);
//	void WriteXYZCell(CTimeStep *ts, int mol, int smol);
//	void WriteXYZCell_Start(CTimeStep *ts, const char *s, int mol, int smol);
	void WriteMoleculeInfo(const char *s);
	void WriteAtomInfo(const char *s);
	void WriteNbAtomDistMatrix(const char *s);
	void WriteNbAtomCountMatrix(const char *s);
	void WriteNbMoleculeDistMatrix(const char *s);
	void WriteNbMoleculeAreaMatrix(const char *s);
	void WriteNbMoleculeCountMatrix(const char *s);
	double m_fBoxDens;
	void Finish();
	void Dump(const char *s, CTimeStep *ts);
	void Parse();
	void Init();
	void Init2();
	void Build(CTimeStep *ts);
	CVoroWrapper();
	~CVoroWrapper();
	int m_iBlocksX;
	int m_iBlocksY;
	int m_iBlocksZ;
	container_periodic_poly *m_pContainer;
	CxObArray m_oaVoroAtoms;
	CxObArray m_oaVoroMolecules;
	long *m_pAssignAtoms;
	long *m_pAssignMolecules;
	long *m_pAssignMoleculeTypes;
	double m_fMaxVol;
	double m_fMaxSurf;
	bool *m_pAtomTouched;
	bool m_bVoroStat;

	CxDoubleArray m_faPOVBoxAtomVisible;
	bool m_bWritePOV;
	int m_iPOVMol;
	int m_iPOVSM;
	bool m_bPOVEdges;
	bool m_bPOVFaces;
//	float m_fPOVFaceOpac;
	bool m_bPOVFaceColor;
	bool m_bPOVAtoms;
	bool m_bPOVRot;
	CxMatrix3 m_mPOVMat;
	bool m_bPOVVertices;
	CxString m_sPOVExe;
	CxFloatArray m_faPOVFaceColor;
	CxWordArray m_waPOVFaceColorMol;
	CxWordArray m_waPOVFaceColorElem;
	CxWordArray m_waPOVFaceColorAtom;
	FILE *m_fPOVScript;
	CxVector3 m_vPOVRotInc;
	CxVector3 m_vPOVRotPos;
	bool m_bPOVAtomGrey;
	float m_fPOVCameraDist;
	int m_iPOVResX;
	int m_iPOVResY;
	bool m_bPOVFaceColorRef;
	bool m_bPOVDrawNeighbors;
	bool m_bPOVNbColorFace;
	bool m_bPOVHighlightNbAtoms;

	bool m_bSurfCover;
	int m_iSurfCoverMol;
//	int m_iSurfCoverSM;
	CxWordArray m_waSurfCoverSM;
	CxObArray m_oaSurfCoverData;
	CxWordArray m_waSurfCoverMol;
	CxWordArray m_waSurfCoverElem;
	CxWordArray m_waSurfCoverAtom;
};


#endif


