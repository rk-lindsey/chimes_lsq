/*****************************************************************************
    TRAVIS - Trajectory Analyzer and Visualizer
    http://www.travis-analyzer.de/

    Copyright (c) 2009-2016 Martin Brehm
                  2012-2016 Martin Thomas

    This file written by Martin Brehm and Martin Thomas.

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

#ifndef GLOBALVAR_H
#define GLOBALVAR_H

#include "asciiart.h"
#include "database.h"
#include "timestep.h"
#include "domain.h"
#include "vorowrapper.h"

#include "tetrapak.h"


class CxString;

extern CxObArray g_oaAtoms;  // Die Atomsorten der Simulation
extern CxObArray g_oaMolecules; // Die Molekuelsorten der Simulation
extern CxObArray g_oaSingleMolecules;
extern CxObArray g_oaObserv;
extern CxObArray g_oaVirtualAtoms;
extern CTimeStep g_TimeStep;
extern char g_iFixMol; // Welches Molekuel fixieren wir?
extern unsigned char g_iFixAtomType[3]; // Welches Element fixieren wir? (3 mal wegen den 3 Fixpunkten)
extern unsigned char g_iFixRealAtomType[3]; // Wie Eben, nur Index auf g_pAtoms
extern unsigned char g_iFixAtom[3]; // Welche Atome in in dem Molekuel fixieren?
extern int g_iCurrentTimeStep, g_iLastTimeStep, g_iNextTimeStep;
extern unsigned long g_iSteps; // Zaehler fuer die schon verarbeiteten Zeitschritte
extern unsigned char g_iBinning;
extern bool g_bFold;
extern bool g_bWarnUnsteady;
extern float g_fUnsteadyLimit;
extern FILE *g_fSaveJustTraj;
extern unsigned char g_iVirtAtomType;
extern float g_fVelLimit, g_fForceLimit;
extern float g_fMSVDFLevel, g_fMSFDFLevel;
extern bool g_bCalcVel, g_bCalcForces;
extern bool g_bManSVDFLevel, g_bManSFDFLevel;
extern bool g_bManSMeanVDFLevel, g_bManSMeanFDFLevel;
extern bool g_bManSMaxVDFLevel, g_bManSMaxFDFLevel;
extern float g_fSVDFLevel[16];
extern float g_fSMeanVDFLevel[16];
extern float g_fSMaxVDFLevel[16];
extern float g_fSFDFLevel[16];
extern float g_fSMeanFDFLevel[16];
extern float g_fSMaxFDFLevel[16];
//extern char g_sRefEnv[256];
extern CxString g_sRefEnv;
extern bool g_bUseVelocities;
extern bool g_bUseForces;
extern float g_fBoxX, g_fBoxY, g_fBoxZ;
extern CxVec3Array g_pRefMol;  // Die Koordinaten des Referenzmolekueles (zum akkumulieren)
extern bool g_bAvg; // Gemitteltes Referenzmolekuel ausgeben
extern bool g_bMiddleAvg; // Soll das Referenzmolekuel wirklich gemittelt werden (true), oder soll einfach das erstbeste Molekuel genommen werden (false)?
extern int g_iGesAtomCount; // Gesamtzahl der Atome pro Zeitschritt - also in der Simulation
extern int g_iGesVirtAtomCount;
extern int g_iMaxStep; // Bis zu welchem Zeitschritt geht die Analyse? -1 fuer alle
extern unsigned char g_iSwapAtoms; // Sollen Atome, die sich frei drehen oder sich vertauschen, zurueckgetauscht werden?
extern unsigned char g_iNormSDF; // 0 - jede SDF einzeln auf 100, 1 - gar nicht normieren
extern unsigned char g_iSVDFLevel, g_iSFDFLevel;
extern bool g_bMegaMat;
extern bool g_bMatOnlyBind;
extern bool g_bVFHisto;
extern unsigned short g_iHistogramRes;
extern FILE *g_fRefTrajec, *g_fRefEnv, *g_fVRDF[32];
extern CTimeStep *g_pTempTimestep;
extern bool g_bFoldAtomwise;
extern float g_fBondFactor;
extern bool g_bSaveJustTraj;
extern bool g_bSaveVirtAtoms;
extern FILE *g_fVFCorr[64];
extern unsigned short g_iVFCorrCount;
extern bool g_bSaveRefWithEnv;
extern bool g_bRefEnvCenter;
extern bool g_bPeriodicX, g_bPeriodicY, g_bPeriodicZ;
extern bool g_bPeriodic;
extern int g_iStride;
extern CAsciiArt g_oAsciiArt;
extern CTimeStep *g_pT2Timestep;
//extern CNbSearch *g_pNbAll;
extern int g_iClusterPos;
extern CxByteArray g_baAtomIndex;
extern bool g_bScanMolecules;
extern CACF *g_pGlobalVACF;
extern CACF *g_pGlobalDipACF;
extern int g_iRefSystemDim;
//extern bool g_bSaveTrajVirt;
extern FILE *g_fPos, *g_fVel, *g_fForce;
extern float g_fTimestepLength;
extern CFFT *g_pFFT;
extern long g_iHighestStepNumber;
extern int g_iDotCounter;
extern bool g_bStepSkipped;
extern bool g_bAbortAnalysis;
extern CxIntArray g_laBondBlackList;
extern int g_iBondBlackListUsed;
extern int g_iFirstStepSkipped;
extern bool g_bGlobalPsycho;
extern int g_iWannierAtomType;
extern float g_fWannierCharge;

extern bool g_bCombinedDF;
extern bool g_bSDF, g_bRDF; 
extern bool g_bVDF, g_bFDF, g_bADF;
extern bool g_bDipole;

extern CxObArray g_oaAnalysisGroups;
extern CxObArray g_oaAnalyses;
extern int g_iStepHistory;
extern CxObArray g_oaTimeSteps;
extern float g_fMaxVel, g_fMaxForce;
extern float g_fLMaxVel, g_fLMaxForce;
extern float g_fLMidVel, g_fLMidForce;
extern bool g_bAsciiArt;

//extern char g_sInputVel[256];
//extern char g_sInputForce[256];
//extern char g_sInputCtrl[256];

extern CxString g_sInputVel;
extern CxString g_sInputForce;
extern CxString g_sInputCtrl;

extern CDatabase *g_pDatabase;

extern CxIntArray g_laAtomSMIndex; // Enthaelt fuer jedes Atom, in welchem SingleMolecule es sich befindet, oder -1 wenn gar nicht
extern CxIntArray g_laAtomSMLocalIndex;
extern CxWordArray g_waAtomElement;
extern CxWordArray g_waAtomMolNumber;
extern CxWordArray g_waAtomMolIndex;
extern CxWordArray g_waAtomRealElement;
extern CxDoubleArray g_faAtomCode;
extern CxDoubleArray g_faVdWRadius;
extern CxWordArray g_waAtomMolUID;

extern CxByteArray g_baAtomPassedCondition;

extern int g_iCDFChannels;
extern bool g_bWannier;

extern int g_iScanMolStep;
extern char *g_sInputTraj;
extern bool g_bCDF;
extern int *g_iObsChannel; // 1 - RDF, 2 - ADF, 3 - DDF, 4 - DipDF, 5 - VDF, 6 - FDF
extern bool g_bDDF;
extern bool g_bRevSDF;
extern bool g_bMSD;
extern bool g_bCutCluster;
extern bool g_bSaveRefEnv;
extern unsigned char g_iSaveRefMol;
extern bool g_bRefEnvFix;
extern int g_iClusterSteps, g_iClusterCount;
extern CxIntArray g_iaClusterSteps, g_iaClusterMol;
extern bool g_bNbAnalysis;
extern bool g_bVHDF;
//extern CxObArray g_oaNbSearches;
extern CNbSet *g_pNbSet;
extern bool g_bVACF;
extern bool g_bVFDF;
extern bool g_bGlobalVACF;
extern bool g_bGlobalDipACF;
extern int g_iSDFSmoothGrade;
extern bool g_bACFWindowFunction;
extern bool g_bSaveVelForce;
extern float g_fVelPercentage, g_fForcePercentage;
//extern bool g_bRefNoVirt;
extern unsigned long g_iSaveGesAtoms;
extern bool g_bUnwrap;
extern CxObArray g_oaSaveMolecules;
extern CxVec3Array g_vaUnwrapArray;
extern int g_iBeginStep;
extern bool g_bSkipDoubleSteps;

extern int g_iNumberPos; // Die wievielte Nummer in der XYZ-Kommentarzeile ist die Schrittzahl?
extern FILE *g_fInput;

extern bool g_bDipACF;
extern bool g_bACF;
extern bool g_bDipDF;

extern bool g_bInputRedirected;

//extern int g_iSDFScale; // 0 = ppm, 1 = pm^-3, 2 = nm^-3, 3 = Rel. to Uniform Density
extern bool g_bSDFUniform;

extern bool g_bDDisp; 
extern bool g_bDLDF; 
extern bool g_bDLDisp; 


extern bool g_bDACF; 
extern bool g_bAggregation;
//extern CAggregation *g_pAggregation;

extern bool g_bKeepUnfoldedCoords;

extern CxObArray g_oaElements;

extern bool g_bTDO;
extern bool g_bTDOEqui;
extern int g_iTDOCount;
extern int g_iTDOStride;
extern int g_iTDOStart;
extern CxIntArray g_laTDOSteps;
extern float g_fTDOBleaching;


extern bool g_bMultiInterval;
extern int g_iMultiIntervalBegin;
extern int g_iMultiIntervalStride;
extern int g_iMultiIntervalLength;
extern CxIntArray g_laMultiIntervalStart;
extern CxIntArray g_laMultiIntervalEnd;

extern bool g_bDoubleBox;
extern int g_iDoubleBoxX;
extern int g_iDoubleBoxY;
extern int g_iDoubleBoxZ;
extern int g_iDoubleBoxFactor;

extern int g_iTrajSteps;

extern bool g_bRDyn;

extern bool g_bBondACF;
extern int g_iBondACFDepth;
extern bool g_bBondACFDebug;
extern bool g_bBondACFNormalize;
extern bool g_bBondACFSymmetrize;
extern bool g_bBondACFWindow;

extern bool g_bACFFFT;

extern bool g_bMSDCacheMode;
extern bool g_bRDynCacheMode;
extern bool g_bVACFCacheMode;

extern char *g_sInputFile;
extern FILE *g_fInputFile;

extern bool g_bSaveCondSnapshot;
extern bool g_bSaveCondWholeBox;
extern FILE *g_fSaveCondFile;
extern int g_iSaveCondCount;

extern bool g_bCond;
//extern CConditionGroup *g_pCondition;

extern bool g_bCombined;

extern unsigned long g_iFastForwardPos;

extern int g_iNbhMode;

extern bool g_bWriteAtomwise;


extern time_t g_iStartTime, g_iEndTime;

extern int g_iTrajFormat;

extern bool g_bScanVelocities;

extern int g_iScanNbhStart;
extern int g_iScanNbhSteps;
extern int g_iScanNbhStride;

extern int g_iScanVelStart;
extern int g_iScanVelSteps;
extern int g_iScanVelStride;


extern bool g_bSilentProgress;

extern char *g_sHomeDir;
extern char *g_sSettingsFile;
extern char *g_sHostName;
extern char *g_sWorkingDir;

extern bool g_bSMode;

extern bool g_bCreateRevSDF;
extern bool g_bNoColor;

extern int g_iColorIntensity;
extern bool g_bNPT;

//extern char g_sNPTFile[256];

extern CxString g_sNPTFile;

extern FILE *g_fNPTFile;

extern bool g_bNbExchange;

extern bool g_bVerbose;
extern bool *g_pUniteTemp;

extern bool g_bCenterZero;

extern bool g_bSaveJustCenter;
extern int g_iSaveJustMol;
extern int g_iSaveJustSM;
extern unsigned char g_iSaveJustAtomType;
extern unsigned char g_iSaveJustRealAtomType;
extern unsigned char g_iSaveJustAtom;

extern bool g_bTimeDiff;
extern bool g_bDeriv;
extern int g_iDerivLast;
extern int g_iDerivCurr;
extern int g_iDerivNext;

extern int g_iCloseAtomCounter;

extern bool g_bVoro;
extern CVoroWrapper *g_pVoroWrapper;

extern bool g_bRemoveCOM;

extern int g_iVoroMemory;

extern bool g_bSaxonize;

extern bool g_bUnknownElements;

extern bool g_bNeedMoleculeWrap;


extern bool g_bAdvanced1;
extern bool g_bAdvanced2;

extern bool g_bVoid;


extern bool g_bDipoleDefined;

extern bool g_bDipolGrimme;


extern bool g_bSaveCoordsUnchanged;

extern bool g_bDens;

extern bool g_bSaveTrajNoRot;

extern bool g_bCheckWrite;

extern bool g_bRaman;
extern bool g_bIRSpec;
extern bool g_bPowerSpec;

extern bool g_bKeepOriginalCoords;

extern bool g_bShowConf;
extern bool g_bWriteConf;
extern char *g_sConfFile;

extern bool g_bPDF;

extern int g_iStrideDetect;

extern float g_fMinPeriodic;


extern bool g_bRegionAnalysis;
extern CxIntArray g_iaSMRegion;

extern bool g_bDumpDipoleVector;
extern bool g_bDumpDipoleAbs;
extern bool g_bDumpDipoleXYZ;
extern CxObArray g_oaDumpDipoleVector; // Contains CxIntAtrrays
extern FILE *g_fDumpDipole;
extern FILE *g_fDumpDipoleXYZ;
extern int g_iDumpDipoleSMCount;
extern int g_iDumpDipoleXYZAtoms;
extern float g_fDumpDipoleScale;

extern bool g_bPlDF;
extern bool g_bLiDF;

extern bool g_bStreamInput;

extern bool g_bGlobalIR;
extern CReorDyn *g_pGlobalIR;

extern bool g_bLMFitSilent;

extern int g_iFitDegree;
extern double *g_pExpSpecExpo;
extern bool g_bLMFitSmooth;

typedef struct {
    const double *x;
    const double *y;
} lmcurve_data_struct;

extern int g_iLMMaxIter;

extern bool g_bUnwrapWannier;
extern bool g_bDipoleRefFixed;

extern double *g_fLSpecEvolveBuf;

extern bool g_bXYZ4thCol;
extern bool g_bXYZComment6Numbers;
extern bool g_bReadChargesFrom4thXYZ;

extern bool g_bPairMSD;

extern bool g_bSFac;
//extern CStructureFactor *g_pSFac;


extern bool g_bEnvWriteDetailedInfo;
extern bool g_bEnvSortNb;
extern bool g_bEnvDisableSortNb;

extern bool g_bShowCredits;

extern bool g_bWriteInputOrder;

extern bool g_bPlProj;

extern bool g_bNormalCoordinate;
extern bool g_bChiral;
extern bool g_bSortWannier;
extern bool g_bVCD;
extern bool g_bEckartTransform;
extern bool g_bPower;
extern bool g_bIR;

extern CTetraPak *g_pTetraPak;

extern bool g_bTegri;

extern CxFloatArray g_faVoronoiRadii;

extern bool g_bBoxNonOrtho;
extern bool g_bFoundNonOrtho;
extern float g_fBoxAngleA;
extern float g_fBoxAngleB;
extern float g_fBoxAngleC;
extern CxMatrix3 g_mBoxToOrtho;
extern CxMatrix3 g_mBoxFromOrtho;
extern float g_fBoxMinDiamA;
extern float g_fBoxMinDiamB;
extern float g_fBoxMinDiamC;
extern float g_fBoxMinDiam;
extern float g_fBoxVolume;
extern bool g_bWriteOrtho;
extern float g_fWriteOrthoFac;

extern bool g_bSDFVoro;
extern int g_iVoroPrintLevel;

extern bool g_bSDFMap;
extern CxObArray g_oaSDFMaps;
extern int g_iSDFMapSmoothGrade;

extern bool g_bCHDF;

extern float g_fCubeXStep;
extern float g_fCubeYStep;
extern float g_fCubeZStep;
extern float g_fCubeXVector[3];
extern float g_fCubeYVector[3];
extern float g_fCubeZVector[3];
extern CxMatrix3 g_mCubeCell;
extern int g_iCubeXStride;
extern int g_iCubeYStride;
extern int g_iCubeZStride;
extern int g_iCubeXMismatch;
extern int g_iCubeYMismatch;
extern int g_iCubeZMismatch;

extern bool g_bCubeTimeDev;
extern FILE *g_fPDESolverInfoFile;
extern float g_fBackgroundDensity;
extern float g_fPDEConvThresh;
extern int g_iPDEMaxIter;

extern bool g_bVoronoiMoments;
extern bool g_bVoroIntEquitable;
extern bool g_bVoroIntegrateCharge;
extern bool g_bVoroIntegrateDipoleMoment;
extern bool g_bVoroIntegrateTotalCurrent;
extern bool g_bVoroIntegrateMagneticMoment;

extern bool g_bVoroSilent;
extern bool g_bDomA;
extern CDomainEngine *g_pDomainEngine;

extern CxString g_sAmberParmFile;



extern bool g_bCubeStream;
extern CxMemFile *g_fCubeMemFile;
extern int g_iCubeMemFileLines;
extern int g_iCubeMemFileSteps;

extern bool g_bDipoleRestart;
extern bool g_bLoadDipoleRestart;
extern FILE *g_fDipoleRestartFile;
extern bool g_bMagneticDipoleRestart;
extern bool g_bLoadMagneticDipoleRestart;
extern FILE *g_fMagneticDipoleRestartFile;


extern bool g_bSetUpPolarizabilityCalc;

extern bool g_bPolarizabilityDefined;
extern int g_iPolarizabilityMode;
extern int g_iPolarizabilityConf[3];
extern FILE *g_fPolarizabilityFile[6];
extern float g_fPolarizabilityFieldStrength;

extern bool g_bLAMMPSCharge;
extern bool g_bReadLAMMPSCharges;

#endif


