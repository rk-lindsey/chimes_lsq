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

#ifndef MAINTOOLS_H
#define MAINTOOLS_H

#include "travis.h"
#include "tools.h"
#include "database.h"
#include "statistics.h"
#include "element.h"
#include "base64.h"


class CAutoCorrelation : public CxObject
{
public:
	CAutoCorrelation();
	~CAutoCorrelation();

	void Init(int input, int depth, bool fft);
	void AutoCorrelate(CxFloatArray *inp, CxFloatArray *outp);

	int m_iInput;
	int m_iDepth;
	int m_iFFTSize;
	bool m_bFFT;
	CFFT *m_pFFT;
	CFFT *m_pFFT2;
};


class CCrossCorrelation : public CxObject                                              
{                                                                                      
public:                                                                                
	CCrossCorrelation();                                                            
	~CCrossCorrelation();                                                           
                                                                                       
	void Init(int input, int depth, bool fft);                                      
	void CrossCorrelate(CxFloatArray *inp1, CxFloatArray *inp2, CxFloatArray *outp);
                                                                                       
	int m_iInput;                                                                   
	int m_iDepth;                                                                   
	int m_iFFTSize;                                                                 
	bool m_bFFT;                                                                    
	CFFT *m_pFFT;                                                                   
	CFFT *m_pFFT2;                                                                  
	CFFT *m_pFFTback;                                                               
};    
                                                                                 

/*void AutoCorrelate(CxFloatArray *inp, CxFloatArray *outp, int depth);
void AutoCorrelate(CxFloatArray *inp, CxFloatArray *outp, int depth, CFFT *fft, CFFT *fft2);*/
CAnalysisGroup* AddAnalysisGroup(const char *name);
void AddAnalysis(CAnalysisGroup* g, const char *name, const char *abbrev);
void InitAnalyses();
void DumpAnalyses();
//void UniteNb();
bool ParseAtom(const char *s, int refmol, unsigned char &ty, unsigned char &rty, unsigned char &atom);
bool ParseRefSystem(int refmol, const char *s, int points);
CTimeStep* GetTimeStep(int i);
CTimeStep** GetTimeStepAddress(int i);
void CalcVelocities();
void CalcVolumetricDataTimeDev();
void CalcCurrentDensity();
void CalcForces();
/*float AtomMass(char *s);
int AtomOrd(char *s);
float AtomRadius(char *s);*/
CElement* FindElement(const char *s, bool quiet);
float GuessBoxSize();
void strtolower(char *s);
void SortAtoms();
bool SetAnalysis(const char *s);
bool ParseFunctions(const char *s);
bool ParsePeriodic(const char *s);
void WriteHeader();
void CommandLineHelp();
bool ParseArgs(int argc, const char *argv[]);
void ParsePassiveArgs(int argc, const char *argv[]);
//void VariablesToDatabase();
//void DatabaseToVariables();
//void WriteDefaultSettings(const char *s);
void CreateDatabaseDefaults();
void LoadSettings();
void InitDatabase();
void RECURSION_BuildCDF(CObservation *o, int channel, int om, CxDoubleArray **data, double *result);
CVirtualAtom* AddVirtualAtom(int mol);
void RemoveAllElements();
void RemoveAllAtoms();
void RemoveAllAnalyses();
void RemoveAllMolecules();
void RemoveAllObservations();
void GetTravisPath();
void ReorderAtoms(int molecule);
void ReorderLikeInput();
void DoubleBoxHelper(unsigned char tpx, unsigned char tpy, unsigned char tpz);
unsigned long GraceColor(int z, double bleach);
void parseCoreCharges();
bool setupWannier();
void ParseDipole();
void parseMagneticDipole();
void DipolGrimme(const char *s);


inline CxVector3 FoldVector(CxVector3 v)
{
	if (g_bBoxNonOrtho)
	{
		CxVector3 w;
		w = g_mBoxToOrtho * v;
		while (w[0] > 0.5)
			w[0] -= 1.0;
		while (w[0] <= -0.5)
			w[0] += 1.0;
		while (w[1] > 0.5)
			w[1] -= 1.0;
		while (w[1] <= -0.5)
			w[1] += 1.0;
		while (w[2] > 0.5)
			w[2] -= 1.0;
		while (w[2] <= -0.5)
			w[2] += 1.0;
		return g_mBoxFromOrtho * w;
	} else
	{
		if (g_bPeriodicX)
		{
			while (v[0] > g_fBoxX/2)
				v[0] -= g_fBoxX;
			while (v[0] <= -g_fBoxX/2)
				v[0] += g_fBoxX;
		}
		if (g_bPeriodicY)
		{
			while (v[1] > g_fBoxY/2)
				v[1] -= g_fBoxY;
			while (v[1] <= -g_fBoxY/2)
				v[1] += g_fBoxY;
		}
		if (g_bPeriodicZ)
		{
			while (v[2] > g_fBoxZ/2)
				v[2] -= g_fBoxZ;
			while (v[2] <= -g_fBoxZ/2)
				v[2] += g_fBoxZ;
		}
	}
	return v;
}


inline float FoldedLength(CxVector3 v)
{
	if (g_bBoxNonOrtho)
	{
		CxVector3 w;
		w = g_mBoxToOrtho * v;
		while (w[0] > 0.5)
			w[0] -= 1.0;
		while (w[0] <= -0.5)
			w[0] += 1.0;
		while (w[1] > 0.5)
			w[1] -= 1.0;
		while (w[1] <= -0.5)
			w[1] += 1.0;
		while (w[2] > 0.5)
			w[2] -= 1.0;
		while (w[2] <= -0.5)
			w[2] += 1.0;
		return (g_mBoxFromOrtho * w).GetLength();
	} else
	{
		if (g_bPeriodicX)
		{
			while (v[0] > g_fBoxX/2)
				v[0] -= g_fBoxX;
			while (v[0] <= -g_fBoxX/2)
				v[0] += g_fBoxX;
		}
		if (g_bPeriodicY)
		{
			while (v[1] > g_fBoxY/2)
				v[1] -= g_fBoxY;
			while (v[1] <= -g_fBoxY/2)
				v[1] += g_fBoxY;
		}
		if (g_bPeriodicZ)
		{
			while (v[2] > g_fBoxZ/2)
				v[2] -= g_fBoxZ;
			while (v[2] <= -g_fBoxZ/2)
				v[2] += g_fBoxZ;
		}
		return v.GetLength();
	}
}


/*inline CxVector3 FoldVector1(CxVector3 v)
{
	while (v[0] >= 0.5)
		v[0] -= 1.0;
	while (v[0] < -0.5)
		v[0] += 1.0;

	while (v[1] >= 0.5)
		v[1] -= 1.0;
	while (v[1] < -0.5)
		v[1] += 1.0;

	while (v[2] >= 0.5)
		v[2] -= 1.0;
	while (v[2] < -0.5)
		v[2] += 1.0;

	return v;
}*/


void BuildAtomIndices();
bool DetermineTrajFormat();
void PrintSMode();
void WriteCredits();
void WriteCredits_Long();
unsigned long CalcFFTSize(unsigned long i, bool silent);

//void FormatTime(unsigned long eta, char *buf);
void FormatTime(unsigned long eta, CxString *buf);

void RenderStructFormulas(int tries);
void RenderFormula(const char *s, int tries);
void InitGlobalVars();

void ParseVoronoiRadii();

void DumpNonOrthoCellData();
void ExtractXYZCellGeometry(const char *s);



#endif
