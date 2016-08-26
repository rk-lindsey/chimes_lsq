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

#ifndef MOLTOOLS_H
#define MOLTOOLS_H

#include <math.h>
#include <stdio.h>
#include "bintools.h"
#include "xobject.h"
#include "xdf.h"
#include "xobarray.h"
#include "xvec3array.h"
#include "xfloatarray.h"
#include "xwordarray.h"
#include "xvector3.h"
#include "xmatrix3.h"
#include "nbsearch.h"
#include "backtrace.h"
#include "statistics.h"
#include "xbytearray.h"
#include "atomgroup.h"
#include "df.h"
#include "2df.h"
#include "3df.h"
//#include "s3df.h"
#include "acf.h"
#include "aggregate.h"
#include "grace.h"
#include "element.h"
#include "reordyn.h"
#include "revsdf.h"
#include "nbexchange.h"

class CDACF;
class CPlProj;
class CTimeStep;
extern int g_iCDFChannels;
//int FindAtom(char *s);


class CAtom : public CxObject
{
public:
	CAtom();
	~CAtom();

	char m_sName[8]; // Elementsymbol
	int m_iCount;    // Wie oft kommt dieses Element in der Simulation vor?
	bool m_bExclude; // Aus der Molekuelerkennung ausschliessen
	float m_fCharge;
	CElement *m_pElement;
	CAtom *m_pMergedTo;
	int m_iIndex;
};


class CVirtualAtom : public CxObject
{
public:
	CVirtualAtom();
	~CVirtualAtom();

	unsigned char m_iMolVirtAtom;
	char m_sLabel[8];
	CAtomGroup m_oCenterAtoms;
	unsigned long m_iIndex;
	unsigned short m_iMolecule;
	unsigned char m_iMode; // 0 - Schwerpunkt, 1 - Bond/Angle/Dihedral, 2 - Dipolvektor, 3 - Geschw.-Vektor, 4 - Kraftvektor
	CxFloatArray m_faWeight;
	float m_fGesWeight;
	unsigned char m_iAtomType[3];
	unsigned char m_iRealAtomType[3];
	unsigned char m_iAtom[3];
	float m_fValues[3];
};


class CMolecule : public CxObject
{
public:
	CMolecule();
	~CMolecule();
	void BuildName();
	int FindAtomInMol(const char *s);
	void Dump();

	CxByteArray m_baAtomIndex; // Welche Atomsorte? -> Verweis auf g_pAtoms
	CxWordArray m_waAtomCount; // Wie viele von dieser Sorte?
	CxIntArray m_laSingleMolIndex; // Verweise auf die ganzen konkreten Molekuele dieser Sorte
	CxIntArray m_laVirtualAtoms;
	CxObArray m_oaNewNumbers; // Enthaelt CxIntArrays mit den neuen Nummern fuer die Atomsorten
	CxObArray m_oaRingAtomTypes; // Enthaelt CxWordArray fuer jeden Ring
	CxObArray m_oaRingAtoms; // Enthaelt CxWordArray fuer jeden Ring
	int m_iAtomGes;       // Gesamtzahl der Atome (inclusive virtueller Atome!)
	int m_iAtomGesNoVirt;
	char *m_sName;     // Name (String)
	int m_iWannierCount;
	float m_fMass;
	int m_iIndex;
	unsigned char m_iDipoleCenterType;   // Atomsorte (fuer Molekuel, nicht global) des Dipolzentrums
	unsigned char m_iDipoleCenterIndex;  // Atomnummer des Dipolzentrums
	float m_fCharge;
	bool m_bPseudo;
	bool m_bChargesAssigned;
	CxObArray m_oaCharges;
	char **m_sFFLabels;
	
	int m_iDipoleMode;
	FILE *m_pDipoleFile;
	int m_iDipoleCommentIndex[3];
	int m_iMagneticDipoleMode;
	FILE *m_pMagneticDipoleFile;
	int m_iMagneticDipoleCommentIndex[3];
	unsigned char m_iMagneticDipoleCenterType;
	unsigned char m_iMagneticDipoleCenterIndex;
	bool m_bPolymer;
};


class CMolAtom : public CxObject
{
public:
	CMolAtom();
	~CMolAtom();
	int m_iMolAtomNumber;
	int m_iType;
	int m_iNumber;
	int m_iOffset;
	double m_fAtomCode;
	double m_fTempAtomCode;
	CxObArray m_oaBonds; // Pointer auf andere CMolAtoms
};


class CSingleMolecule : public CxObject
{
public:
	int m_iMolSMIndex;
	int CountDifferentAtomCodes();
	void BuildAtomCodes();
	CSingleMolecule();
	~CSingleMolecule();
	void Dump();

	CxObArray m_oaBondGroups;
	CxObArray m_oaBonds;
	CxObArray m_oaAngleGroups;
	CxObArray m_oaAngles;

	CxObArray m_oaRings; // Enthaelt ein CxWordArray fuer jeden Ring

	CxIntArray m_laBonds;
	CxIntArray m_laWannier;
	CxObArray m_oaAtomOffset; // Enthaelt ein CxIntArray fuer jede Atomsorte
	CxObArray m_oaMolAtoms; // Enthaelt ein CMolAtom fuer jedes Atom
	int m_iAtomClasses;
	int m_iIndex;
	CxByteArray m_baAtomIndex; 
	int m_iMolType;
	CxVector3 m_vDipole;
	CxVector3 m_vMagneticDipole;
	float m_polarizability[9];
	bool m_bPseudo;
//	CxObArray m_oaAggregates;
};


class CADF : public CXDF
{
public:
	void CopyFrom(CADF *p);
	void ParseCondition(int rm, bool nocrit);
	void ParseConditionGrid(int rm, int gridmode);
	CADF();
	~CADF();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void Parse();
	void ParseCondition_OnlyValues();

	bool m_bOrtho[2];
	int m_iRefOrSec[2][3]; // 0 = Ref, 1 = Obs, 2 = Obs2
	int m_iVecType[2]; // 0 = Ort, 1 = Dipol, 2 = Geschw., 3 = Kraft
	bool m_bStat;
//	bool m_bSameFoot;
	bool m_bMirror;
//	bool m_bSaveAngle;
	FILE **m_fAngle;
	float m_fMinAngle, m_fMaxAngle;
	CxFloatArray m_faMinMaxAngle;
	bool m_bFoldAngle;
	CDF *m_pADF;
	void BuildName();
//	bool m_bOthers;
	bool m_bCosine;
	CxObArray m_oaVectors; // Enthaelt CAtomGroups, 6 pro Vektor-Satz

};


class CDDF : public CXDF
{
public:
	CDDF();
	~CDDF();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();

	bool m_bOrtho[3];
	int m_iRefOrSec[3][3];
	CxObArray m_oaVectors; // Enthaelt CAtomGroups, 9 pro Vektor-Satz
	bool m_bClassical;
	float m_fMinAngle, m_fMaxAngle;
	CDF *m_pDDF;
//	bool m_bSaveAngle;
	FILE **m_fAngle;
//	bool m_bOthers;
	bool m_bCosine;
	bool m_bAbs;
	bool m_bPositive;
	bool m_bSymm;
	bool m_bRotate;
	CxDoubleArray m_faLastData; 
	CxIntArray m_laRotation;
};


class CMSD : public CxObject
{
public:
	void WriteSplit(const char *s);
	void BuildName();
	void Parse();
	CMSD();
	~CMSD();

	int m_iShowAtoms;
	int m_iShowMol;
	int m_iResolution;
	int m_iStride;
	int m_iStride2;
	CAF *m_pMSD;
	CAF **m_pSplitMSD;
	char *m_sName;
	CAtomGroup *m_pAtomGroup;
	CxObArray m_oaCache;
	bool m_bSplit;
	bool m_bTakeX;
	bool m_bTakeY;
	bool m_bTakeZ;
};


class CRDF : public CXDF
{
public:
	bool m_bLongMode;
	void CopyFrom(CRDF *p);
	CRDF();
	~CRDF();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();
	void ParseCondition(int rm, CNbSearch *n, bool nbana);
	void ParseConditionGrid(int rm, CNbSearch *n, int gridmode);
	void ParseCondition_OnlyValues(CNbSearch *n);

	int m_iRefOrSec[2];
	float m_fMinDist, m_fMaxDist;
	CxFloatArray m_faMinMaxDist;
	FILE **m_fDist;
//	CStatistics *m_pDistStat;
	CDF *m_pRDF;
	CxObArray m_oaVectors; // Enthaelt CAtomGroups, 2 pro Vektor-Satz
	int m_iShowAtomGes;
	int m_iRefAtomGes;
	bool m_bAdaptive;
	bool m_bRadialCorrect;
	bool m_bProbDens;
	bool m_bCalcSD;
	int m_iSDBlocks;
	int m_iSDBlockMin;
	int m_iSDBlockMax;
	bool m_bSDVerbose;
	double m_fSDTimesSigma;
};


class CPlDF : public CXDF
{
public:
//	void CopyFrom(CPlDF *p);
	CPlDF();
	~CPlDF();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();

	int m_iRefOrSec[4];
	bool m_bNormal;
	float m_fMinDist, m_fMaxDist;
	CDF *m_pPlDF;
	CxObArray m_oaVectors; // Enthaelt CAtomGroups, 4 pro Vektor-Satz
	int m_iShowAtomGes;
	int m_iRefAtomGes;
};


class CLiDF : public CXDF
{
public:
//	void CopyFrom(CLiDF *p);
	CLiDF();
	~CLiDF();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();

	int m_iRefOrSec[4];
	float m_fMinDist, m_fMaxDist;
	CDF *m_pLiDF;
	CxObArray m_oaVectors; // Enthaelt CAtomGroups, 4 pro Vektor-Satz
	int m_iShowAtomGes;
	int m_iRefAtomGes;
	bool m_bRadialCorrect;
	bool m_bNormal;
};


class CDensDF : public CXDF
{
public:
	void CopyFrom(CRDF *p);
	CDensDF();
	~CDensDF();
	void BuildName();
	void Parse();
	float m_fMinDist, m_fMaxDist;
	unsigned char m_iCenterAtomType;
	unsigned char m_iCenterAtomRealType;
	unsigned char m_iCenterAtom;
	CDF *m_pDensDF;
	bool m_bDensityMass;
	bool *m_pDensityMolSelect;
	CAtomGroup **m_pDensityMolAG;
};


class CVHDF : public CxObject
{
public:
	bool m_bSelf;
	void CorrectCount();
	CVHDF();
	~CVHDF();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();

	char *m_sName;
	char *m_sShortName;
	int m_iRefOrSec[2];
	float m_fMinDist, m_fMaxDist;
	int m_iResolution;
	int m_iDepth;
	int m_iStride;
	C2DF *m_pVHDF;
	double *m_pCount;
	CxObArray m_oaVectors; // Enthaelt CAtomGroups, 2 pro Vektor-Satz
	int m_iShowAtomGes;
	int m_iRefAtomGes;
	int m_iShowMol;
	int m_iCombinations;
	bool m_bRadialCorrect;
	bool m_bSwapAxes;
	int m_iGraceBunchTime;
	int m_iGraceBunchDist;
};


class CVDF : public CXDF
{
public:
	bool m_bScanRange;
	CVDF();
	~CVDF();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();

	int m_iRefOrSec;
	float m_fMinSpeed, m_fMaxSpeed;
	FILE **m_fSpeed;
//	bool m_bSaveSpeed;
	CDF *m_pVDF;
	CAtomGroup m_oAtoms;
	int m_iShowAtomGes;
//	int m_iRefAtomGes;
};


class CDipDF : public CXDF
{
public:
	CDipDF();
	~CDipDF();
	void BuildName();
	void Parse();

	int m_iRefOrSec;
	float m_fDipoleMin, m_fDipoleMax;
//	bool m_bSaveDipole;
	FILE **m_fDipole;
	CDF *m_pDipoleDF;
};


class CSDF : public CxObject
{
public:
	bool m_bVdWSpheres;
	bool m_bCutPlane;
	void CreateCutPlane();
	bool m_bIntra;
	CSDF();
	~CSDF();
	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse(bool voro);

	int m_iRefOrSec;
	char *m_sName;
	float m_fRadius;
	int m_iResolution;
	int m_iHistogramRes;
	bool m_bSDFMirrorXY;
	bool m_bSDFMirrorBisect;
	int m_iShowMol;
	CxVec3Array *m_vaData; 
	CxFloatArray *m_faRadius; 
	C3DF<double> *m_pSDF;
	CAtomGroup m_oAtoms;
	double m_fParticleDensity;
	int m_iShowAtomGes;

	bool m_bClipPlane;
	int m_iClipDirection;
	double m_fClipValue;
	bool m_bInvert;

	double m_fPosCounter;
	double m_fAtom2PosX;
	double m_fAtom3PosX;
	double m_fAtom3PosY;

	int m_iCutPlaneResolution;
	bool m_bCutPlaneShowAtoms;
	C2DF *m_pCutPlane;
	CxByteArray *m_baDataEnabled; 

};


class CCDF : public CxObject
{
public:
	double m_fFactor;
	CCDF();
	~CCDF();

	inline void AddToBin(double *d)
	{
		int z;
		if (m_bDumpDat)
		{
			for (z=0;z<g_iCDFChannels;z++)
				if (z < g_iCDFChannels-1)
					mfprintf(m_fDump,"%f; ",d[z]);
						else mfprintf(m_fDump,"%f\n",d[z]);
		}
		if (g_iCDFChannels == 2)
			m_p2DF->AddToBin(d[0],d[1]);
	}

//	bool *m_bChannelAll; // Soll dieser Kanal nur das aktuelle OM oder alle liefern?
	char *m_sName;
	char *m_sShortName;
	FILE **m_fTimeDev;
	C2DF *m_p2DF;
	C2DF *m_pTensorProduct;
	C2DF *m_pDiffFunction;
	C3DF<double> *m_p3DF;
	CNDF *m_pNDF;
	bool m_bDumpDat;
//	bool m_bTimeDev;
	FILE *m_fDump;
	char m_iNormalize;
	double m_fNormValue;

	int *m_iCombinations;
	int m_iCombinationProd;
	int m_iCombinationsEnabled;
	char *m_pCombineList;

	int *m_iResolution;

	bool m_bGraceBunch;
	int m_iGraceBunchC1;
	int m_iGraceBunchC2;

	bool m_b3DSlices;
	int m_i3DSliceIntervals[3];

	bool m_bAxisDivide;
	bool m_bAxisDivideAll;

	bool m_bTDAnimation;
	int m_iTDASteps;
	int m_iTDAStride;
	int m_iTDATail;
	int m_iTDAResX;
	int m_iTDAResY;
	bool m_bTDATrace; 
	CGrace *m_pTDAPlot;
	int m_iHistogramRes;
};


class CConditionSubGroup : public CxObject
{
public:
	CNbSearch* AddSingleCondition(int rm, int sm, int gridmode);
	void ReScan(CSingleMolecule *rm);
	void CopyResults(CConditionSubGroup *p);
	void Parse_OnlyValues();
	void CopyFrom(CConditionSubGroup *p);
	void Reset();
	int m_iShowMol;
	int m_iCombinations;
	void PrintSingle(int om);
	int m_iNumber;
	void PrintData();
	void PrintData(FILE *a);
	void Parse(int rm, int sm);
//	void ScanNeighborhood(CTimeStep *t, int rm, int mol, bool markpassedatoms);
	void MarkPassedAtoms(int om);
//	void ScanNeighborhoodSingleOM(CTimeStep *t, CSingleMolecule *rm, int om, bool markpassedatoms, bool fold);
	void ScanNeighborhoodAllOM(CTimeStep *t, CSingleMolecule *rm);
	void PreScanNeighborhoodAllOM(CTimeStep *t, CSingleMolecule *rm);
//	void Clear();
	bool Contains(int mol);
	CConditionSubGroup();
	~CConditionSubGroup();

	CxObArray m_oaConditions;
	double m_fPassed;
	double m_fTotal;
	bool *m_bTempPassed;
};


class CConditionGroup : public CxObject
{
public:
	CNbSearch* AddSingleCondition(int rm, int sm, int gridmode);
	void ReScan(CSingleMolecule *rm);
	void CopyResults(CConditionGroup *p);
	void Parse_OnlyValues();
	void CopyFrom(CConditionGroup *p);
	void Reset();
	void PrintTable();
	void PrintSingle(int om);
	void PrintData();
	void PrintData(const char *s);
	CConditionGroup();
	~CConditionGroup();
	void Parse(int rm, int sm);
//	void ScanNeighborhood(CTimeStep *t, int rm, int mol, bool markpassedatoms);
	void MarkPassedAtoms(int om, bool passed);
//	void ScanNeighborhoodSingleOM(CTimeStep *t, CSingleMolecule *rm, int om, bool markpassedatoms, bool fold);
	void ScanNeighborhoodAllOM(CTimeStep *t, CSingleMolecule *rm);
	void PreScanNeighborhoodAllOM(CTimeStep *t, CSingleMolecule *rm);
//	void Clear();
	bool Contains(int mol);

	bool m_bInactive;
	int m_iShowMol;
	int m_iRefMol;
	unsigned long m_iHistoGes;
	bool m_bAnyPassed;
	double m_fTableGes;
	double m_fPassed;
	double m_fTotal;
	bool m_bInvertCondition;

	CxObArray m_oaConditionSubGroups;

	double *m_pTable;
	unsigned long *m_pHistogram;
	bool *m_bAlwaysTrue;
	long *m_iPassCounter;
	long *m_iRMPassCounter;
	long *m_iOMPassCounter;
//	bool m_bNeedNbCountMode;
};


class CObservation : public CxObject
{
public:	
	void CreateTimeDiff(CDF *df, int comb);
	void WriteTimeDiff(CDF *df, const char *anaup, const char *analow, const char *name, const char *multibuf, bool ddf);
	void BuildTimeDiff(CDF *df, bool ddf);
	void ListCDFObservations(int z);
	CObservation();
	~CObservation();

	int m_iShowMol; // In welchem Molekuel beobachten wir Atome?
	int m_iShowMolCount;
	bool m_bSecondShowMol;
	int m_iShowMol2; 
	int m_iShowMol2Count;
	bool m_bExclude1eq2;
	bool m_bSelf;
	bool m_bOthers;
	bool m_bTimeDev;
	bool m_bSaveSeparateFiles;
	bool m_bVelocityRelToRef;
	CxWordArray m_waSaveRefList;
	CxWordArray m_waSaveShowList;

	CxWordArray m_waObsRefList;
	CxWordArray m_waObsShowList;
	CxWordArray m_waObsShow2List;
	bool m_bObsCertain;
	bool m_bDecompDist;

	bool m_bDecompType;
	CxWordArray m_waDecompTypeRefOffs;
	CxWordArray m_waDecompTypeObsOffs;
	CxWordArray m_waDecompTypeRefList;
	CxWordArray m_waDecompTypeObsList;

	bool m_bTimeDiff;
	bool m_b3DTimeDiff;
	int m_iTimeDiffDepth;
	int m_iTimeDiffRes3D;
	int m_iTimeDiffStride3D;
	double m_fTimeDiffMinVal3D;
	double m_fTimeDiffMaxVal3D;
	int m_iTimeDiffDistSteps;
	int m_iTimeDiffDistResX;
	int m_iTimeDiffDistResY;
	double m_fTimeDiffDistMinValX;
	double m_fTimeDiffDistMaxValX;
	double m_fTimeDiffDistMinValY;
	double m_fTimeDiffDistMaxValY;

	bool m_bCombinedGreyMode;
	int m_iCombinedGreyMin;
	int m_iCombinedGreyMax;
	int m_iCombinedGreyShades;

	bool m_bCombinedPlot;

	bool m_bBinOnlyPassedAtoms;
	bool m_bBinOnlyNotPassedAtoms;

	CxIntArray m_iaRMRegions;
	CxIntArray m_iaOM1Regions;
	CxIntArray m_iaOM2Regions;

	CConditionGroup *m_pConditions;
	CConditionGroup *m_pConditionsOM2;

	CReorDyn *m_pRDyn;
	CReorDyn *m_pIRSpec;

	CPlProj  *m_pPlProj;

	CDensDF  *m_pDensityDF;

	CRevSDF  *m_pRevSDF;

	CVHDF    *m_pVHDF;
	CMSD     *m_pMSD;
	CSDF     *m_pSDF;
	CCDF     *m_pCDF;
	CACF     *m_pVACF;
	CACF     *m_pDipACF;
	CADF    **m_pADF;
	CDDF    **m_pDDF;
	CRDF    **m_pRDF;
	CDipDF  **m_pDipDF;
	CVDF    **m_pVDF;
	CPlDF   **m_pPlDF;
	CLiDF   **m_pLiDF;

	CNbAnalysis *m_pNbAnalysis;
	CNbExchange *m_pNbExchange;
	CDACF *m_pDACF;
};


#endif
