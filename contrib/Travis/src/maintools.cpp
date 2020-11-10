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

#include "maintools.h"
#include "conversion.h"

#include "interface.h"
#include "bicgstab.h"

#ifdef TARGET_LINUX
#include <unistd.h>
#endif

CAnalysisGroup* AddAnalysisGroup(const char *name)
{
	BTIN;

	CAnalysisGroup *g;
	try { g = new CAnalysisGroup(); } catch(...) { g = NULL; }
	if (g == NULL) NewException((double)sizeof(CAnalysisGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { g->m_sGroupName = new char[strlen(name)+1]; } catch(...) { g->m_sGroupName = NULL; }
	if (g->m_sGroupName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(g->m_sGroupName,name);
	g_oaAnalysisGroups.Add(g);
	BTOUT;
	return g;
}


void AddAnalysis(CAnalysisGroup* g, const char *name, const char *abbrev)
{
	BTIN;

	CAnalysis *a;
	try { a = new CAnalysis(); } catch(...) { a = NULL; }
	if (a == NULL) NewException((double)sizeof(CAnalysis),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { a->m_sName = new char[strlen(name)+1]; } catch(...) { a->m_sName = NULL; }
	if (a->m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { a->m_sAbbrev = new char[strlen(abbrev)+1]; } catch(...) { a->m_sAbbrev = NULL; }
	if (a->m_sAbbrev == NULL) NewException((double)(strlen(abbrev)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(a->m_sName,name);
	strcpy(a->m_sAbbrev,abbrev);
	g_oaAnalyses.Add(a);
	g->m_oaAnalyses.Add(a);
	BTOUT;
}


void InitAnalyses()
{
	BTIN;
	CAnalysisGroup *g;

	/**********************************************************************/
	g = AddAnalysisGroup("Static (time independent) Functions");

	AddAnalysis(g,"Combined Distribution Function","cdf");
	AddAnalysis(g,"Radial Distribution Function","rdf");
	AddAnalysis(g,"Angular Distribution Function","adf");
	AddAnalysis(g,"Dihedral Distribution Function","ddf");
	AddAnalysis(g,"Point-Plane Distance Distribution","pldf");
	AddAnalysis(g,"Point-Line Distance Distribution","lidf");

	AddAnalysis(g,"Plane Projection Distribution","plproj");
	AddAnalysis(g,"Fixed Plane Density Profile","dprof");
	AddAnalysis(g,"Density Distribution Function","dens");
	AddAnalysis(g,"Spatial Distribution Function","sdf");
	AddAnalysis(g,"Pseudo SDF (only 2 ref. atoms)","psdf");
	AddAnalysis(g,"Dipole Distribution Function","dip");
	AddAnalysis(g,"Evaluate structural condition","cond");
	AddAnalysis(g,"Compute Structure Factor","sfac");



	/**********************************************************************/
	g = AddAnalysisGroup("Dynamic (time dependent) Functions");

	AddAnalysis(g,"Velocity Distribution Function","vdf");
	AddAnalysis(g,"Force Distribution Function","fdf");
	AddAnalysis(g,"Mean Square Displacement / Diffusion Coefficients","msd");
	AddAnalysis(g,"Velocity Autocorrelation Functions","acf");
	AddAnalysis(g,"Vector Reorientation Dynamics","rdyn");
	AddAnalysis(g,"Van Hove Correlation Function","vhcf");
	AddAnalysis(g,"Aggregation Functions (DACF, DLDF, DDisp)","aggr");



	/**********************************************************************/
	g = AddAnalysisGroup("Spectroscopic Functions");

	AddAnalysis(g,"Calculate Power Spectrum","power");
	AddAnalysis(g,"Calculate IR Spectrum","ir");
	AddAnalysis(g,"Calculate Raman Spectrum","raman");
	AddAnalysis(g, "Normal coordinate analysis", "nc");
	AddAnalysis(g, "Calculate VCD Spectrum", "vcd");
	AddAnalysis(g, "Save dipole restart file", "drst");
	AddAnalysis(g, "Save magnetic moment restart file", "mrst");
	AddAnalysis(g, "Set up polarizability calculation", "pol");

	/**********************************************************************/
	g = AddAnalysisGroup("Miscellaneous Functions");

	AddAnalysis(g,"Save trajectory of RM environment / TDO Plot","env");
	AddAnalysis(g,"Save processed Trajectory","proc");
	AddAnalysis(g,"Cut Clusters","cut");
	AddAnalysis(g,"Region-specific Analysis","reg");
	AddAnalysis(g,"Chirality Analysis", "chi");
	AddAnalysis(g,"Transform to Eckart Frame", "eck");
	AddAnalysis(g,"Sort Wannier Centers", "swan");
	AddAnalysis(g,"Basic Voronoi Analysis","voro");
	AddAnalysis(g,"Domain Analysis","doma");
	AddAnalysis(g,"Voronoi Integration Functions","vori");
	
	/**********************************************************************/
	BTOUT;
}


void DumpAnalyses()
{
	BTIN;
	int z, z2;
	CAnalysis *a;
	CAnalysisGroup *g;

	for (z=0;z<g_oaAnalysisGroups.GetSize();z++)
	{
		g = (CAnalysisGroup*)g_oaAnalysisGroups[z];
		mprintf(YELLOW," *** %s\n",g->m_sGroupName);
		for (z2=0;z2<g->m_oaAnalyses.GetSize();z2++)
		{
			a = (CAnalysis*)g->m_oaAnalyses[z2];
			mprintf(WHITE," %-7s",a->m_sAbbrev);
			mprintf("- %s\n",a->m_sName);
		}
		mprintf("\n");
	}
	BTOUT;
}


/*void UniteNb()
{
	BTIN;
	int z, z2, z3, z4;
	CNbSearch *nb;
	CxIntArray *w, *w2;

	if (g_pNbAll != NULL)
		delete g_pNbAll;
	g_pNbAll = new CNbSearch();
	g_pNbAll->Create();

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		w2 = (CxIntArray*)g_pNbAll->m_oaScanNeighbors[z];
		for (z2=0;z2<g_oaNbSearches.GetSize();z2++)
		{
			nb = (CNbSearch*)g_oaNbSearches[z2];
			w = (CxIntArray*)nb->m_oaScanNeighbors[z];
			for (z3=0;z3<nb->m_waScanNeighborCount[z];z3++)
			{
				if (!g_bKeepNbCount)
					for (z4=0;z4<g_pNbAll->m_waScanNeighborCount[z];z4++)
					{
						if (w->GetAt(z3) == w2->GetAt(z4))
							goto _unb_found;
					}
				w2->Add(w->GetAt(z3));
				g_pNbAll->m_waScanNeighborCount[z]++;
_unb_found:;
			}
		}
	}
	BTOUT; 
}*/


bool ParseAtom(const char *s, int refmol, unsigned char &ty, unsigned char &rty, unsigned char &atom)
{
	BTIN;
	CMolecule *mol;
	char buf[64];
	const char *r, *t;
	int a, b;

	mol = (CMolecule*)g_oaMolecules[refmol];
	t = s;
	while ((*t == ' ') && (*t != 0))
		t++;
	r = t;
	while ((!isdigit(*r)) && (*r != 0))
		r++;
	if (r-t == 0)
	{
		eprintf("ParseAtom(): Missing Atom Label for Atom.\n");
		BTOUT;
		return false;
	}
	if (r-t >= 64)
	{
		eprintf("Internal Error in ParseAtom(): A, %d >= 64.\n",r-t);
		return false;
	}
	memcpy(buf,t,r-t);
	buf[r-t] = 0;
//	printf("Label \"%s\"\n",buf);
	a = mol->FindAtomInMol(buf);
	if (a == -1)
	{
		eprintf("ParseAtom(): Atom Type \"%s\" not found in Molecule.\n",buf);
		BTOUT;
		return false;
	}
	b = atoi(r);
	if (b != 0)
		b--;
	if (b >= ((CMolecule*)g_oaMolecules[refmol])->m_waAtomCount[a])
	{
		eprintf("ParseAtom(): The Molecule only has %d %s-Atoms.\n",((CMolecule*)g_oaMolecules[refmol])->m_waAtomCount[a],buf);
		BTOUT;
		return false;
	}
//	printf("ParseRefSystem: Atom 1 ist in Molekuel %d Typ %d Nummer %d.\n",refmol,a,b);
	ty = a;
	rty = ((CMolecule*)g_oaMolecules[refmol])->m_baAtomIndex[a];
	atom = b;
	BTOUT;
	return true;
}


bool ParseRefSystem(int refmol, const char *s, int points)
{
	BTIN;
	char buf[256];
	char *p, *q;

	if (points == 1)
	{
		BTOUT;
		return ParseAtom(s,refmol,g_iFixAtomType[0],g_iFixRealAtomType[0],g_iFixAtom[0]);
	}
	strcpy(buf,s);
	if (points == 2)
	{
		p = strchr(buf,',');
		if (p == NULL)
		{
			eprintf("ParseRefSystem(): No comma found.\n");
			BTOUT;
			return false;
		}
		*p = 0;
		p++;
		if  (!ParseAtom(buf,refmol,g_iFixAtomType[0],g_iFixRealAtomType[0],g_iFixAtom[0]))
		{
			BTOUT;
			return false;
		}
		BTOUT;
		return ParseAtom(p,refmol,g_iFixAtomType[1],g_iFixRealAtomType[1],g_iFixAtom[1]);
	}
	if (points == 3)
	{
		p = strchr(buf,',');
		if (p == NULL)
		{
			eprintf("ParseRefSystem(): No comma found.\n");
			BTOUT;
			return false;
		}
		*p = 0;
		p++;
		q = strchr(p,',');
		if (q == NULL)
		{
			eprintf("ParseRefSystem(): No second comma found.\n");
			BTOUT;
			return false;
		}
		*q = 0;
		q++;
		if  (!ParseAtom(buf,refmol,g_iFixAtomType[0],g_iFixRealAtomType[0],g_iFixAtom[0]))
		{
			BTOUT;
			return false;
		}
		if  (!ParseAtom(p,refmol,g_iFixAtomType[1],g_iFixRealAtomType[1],g_iFixAtom[1]))
		{
			BTOUT;
			return false;
		}
		BTOUT;
		return ParseAtom(q,refmol,g_iFixAtomType[2],g_iFixRealAtomType[2],g_iFixAtom[2]);
	}
	BTOUT;
	return false;
}


CTimeStep* GetTimeStep(int i)
{
	BTIN;
	if (g_iCurrentTimeStep != -1)
	{
		if (g_bUseVelocities || g_bUseForces)
			i++;
		if (g_iCurrentTimeStep-i+g_iStepHistory < 0)
		{
			eprintf("GetTimeStep(): Error! i=%d, g_iCurrentTimeStep=%d, g_iStepHistory=%d.\n",i,g_iCurrentTimeStep,g_iStepHistory);
			BTOUT;
			return NULL;
		}
		if (g_iCurrentTimeStep-i >= 0)
		{
			BTOUT;
			return (CTimeStep*)g_oaTimeSteps[g_iCurrentTimeStep-i];
		} else
		{
			BTOUT;
			return (CTimeStep*)g_oaTimeSteps[g_iCurrentTimeStep-i+g_iStepHistory];
		}
	} else 
	{
		BTOUT;
		return NULL;
	}
}


CTimeStep** GetTimeStepAddress(int i)
{
	BTIN;
	if (g_iCurrentTimeStep != -1)
	{
		if (g_bUseVelocities || g_bUseForces)
			i++;
		if (g_iCurrentTimeStep-i >= 0)
		{
			BTOUT;
			return (CTimeStep**)&g_oaTimeSteps[g_iCurrentTimeStep-i];
		} else
		{
			BTOUT;
			return (CTimeStep**)&g_oaTimeSteps[g_iCurrentTimeStep-i+g_iStepHistory];
		}
	} else
	{
		BTOUT;
		return NULL;
	}
}


void CalcVelocities()
{
	BTIN;
	int z;
	float v, imv, ilmv;

	g_fLMidVel = 0;
	imv = g_fMaxVel;
	ilmv = g_fLMaxVel;
	GetTimeStep(0)->m_vaVelocities.SetSize(g_iGesVirtAtomCount);
	for (z=0;z<g_iGesVirtAtomCount;z++)
	{
//		if (z == 0)
//			mprintf("\n(%f - %f) / (2 * %f) * 100000 = %f.",GetTimeStep(1)->m_vaCoords[z][0],GetTimeStep(-1)->m_vaCoords[z][0],g_fTimestepLength,(GetTimeStep(1)->m_vaCoords[z][0] - GetTimeStep(-1)->m_vaCoords[z][0]) / (2.0f*g_fTimestepLength) * 1000.0f);
		
		CxVector3 dist = GetTimeStep(-1)->m_vaCoords[z] - GetTimeStep(1)->m_vaCoords[z];
// 		mprintf(RED, "d: ( %12G | %12G | %12G )\n", dist[0], dist[1], dist[2]);
		if(g_bPeriodicX) {
			while(dist[0] > g_fBoxX / 2.0f) dist[0] -= g_fBoxX;
			while(dist[0] < -g_fBoxX / 2.0f) dist[0] += g_fBoxX;
		}
		if(g_bPeriodicY) {
			while(dist[1] > g_fBoxY / 2.0f) dist[1] -= g_fBoxY;
			while(dist[1] < -g_fBoxY / 2.0f) dist[1] += g_fBoxY;
		}
		if(g_bPeriodicZ) {
			while(dist[2] > g_fBoxZ / 2.0f) dist[2] -= g_fBoxZ;
			while(dist[2] < -g_fBoxZ / 2.0f) dist[2] += g_fBoxZ;
		}
		GetTimeStep(0)->m_vaVelocities[z] = dist * 1000.0f / 2.0f / g_fTimestepLength;
// 		mprintf(RED, "v: ( %12G | %12G | %12G )\n", GetTimeStep(0)->m_vaVelocities[z][0], GetTimeStep(0)->m_vaVelocities[z][1], GetTimeStep(0)->m_vaVelocities[z][2]);
		
// 		GetTimeStep(0)->m_vaVelocities[z] = (GetTimeStep(1)->m_vaCoords[z] - GetTimeStep(-1)->m_vaCoords[z]) / (2.0f*g_fTimestepLength) * 1000.0f;
		v = GetTimeStep(0)->m_vaVelocities[z].GetLength();
		g_fLMidVel += v;
		if (v > imv)
			imv = v;
		if (v > ilmv)
			ilmv = v;
	}
	g_fLMidVel /= g_iGesVirtAtomCount;
	g_fLMaxVel = ilmv;
	if (imv < g_fUnsteadyLimit)
	{
		g_fMaxVel = imv;
	} else if (!g_bWarnUnsteady)
	{
		g_bWarnUnsteady = true;
		mprintf("Warning: Discontinuity in time step %d (maximum velocity = %f > %f pm/ps).\n",(int)g_iSteps,imv,g_fUnsteadyLimit);
	}
	BTOUT; 
}

void CalcVolumetricDataTimeDev() {
	if (GetTimeStep(0)->m_pVolumetricDataTimeDev == NULL)
		return;
	if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin == NULL) {
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] = GetTimeStep(0)->m_pVolumetricData->m_iRes[0];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] = GetTimeStep(0)->m_pVolumetricData->m_iRes[1];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] = GetTimeStep(0)->m_pVolumetricData->m_iRes[2];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMinVal[0] = GetTimeStep(0)->m_pVolumetricData->m_fMinVal[0];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMaxVal[0] = GetTimeStep(0)->m_pVolumetricData->m_fMaxVal[0];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMinVal[1] = GetTimeStep(0)->m_pVolumetricData->m_fMinVal[1];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMaxVal[1] = GetTimeStep(0)->m_pVolumetricData->m_fMaxVal[1];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMinVal[2] = GetTimeStep(0)->m_pVolumetricData->m_fMinVal[2];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMaxVal[2] = GetTimeStep(0)->m_pVolumetricData->m_fMaxVal[2];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->Create();
	}
	
	int i;
	double integral = 0.0;
	int timedevSize = 0;
	for (i = 0; i < GetTimeStep(0)->m_pVolumetricData->m_iRes[0] * GetTimeStep(0)->m_pVolumetricData->m_iRes[1] * GetTimeStep(0)->m_pVolumetricData->m_iRes[2]; i++) {
		if (fabs(GetTimeStep(0)->m_pVolumetricData->m_pBin[i]) > -1.0e-20) {
			GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i] = (GetTimeStep(-1)->m_pVolumetricData->m_pBin[i] - GetTimeStep(1)->m_pVolumetricData->m_pBin[i]) / (2.0 * g_fTimestepLength);
			integral += GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i];
			timedevSize++;
		} else {
			GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i] = 0.0;
		}
	}
// 	mprintf(GREEN, "%g %g\n", integral, integral / timedevSize);
	for (i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2]; i++) {
		if (fabs(GetTimeStep(0)->m_pVolumetricData->m_pBin[i]) > -1.0e-20) {
			GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i] -= integral / timedevSize;
		}
	}
// 	integral = 0.0;
// 	for (i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2]; i++) {
// 		integral += GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i];
// 	}
// 	mprintf(GREEN, "%g\n", integral);
	
// 	mprintf(GREEN, "%#20.10g %#20.10g %#20.10g %#20.10g\n", GetTimeStep(-1)->m_pVolumetricData->m_pBin[10], GetTimeStep(0)->m_pVolumetricData->m_pBin[10], GetTimeStep(1)->m_pVolumetricData->m_pBin[10], GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[10]);
	
// 		FILE *cubeFile = fopen("test.cube", "w");
// 		if (cubeFile == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		fprintf(cubeFile, "\n\n");
// 		fprintf(cubeFile, "%5lu %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
// 		for (int i = 0; i < (int)GetTimeStep(0)->m_iGesAtomCount; i++) {
// 			fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, GetTimeStep(0)->m_vaCoords[i][0] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][1] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][2] / 0.529177249f / 100.0f);
// 		}
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				for (int k = 0; k < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6; k++) {
// 					for (int l = 0; l < 6; l++) {
// 						fprintf(cubeFile, "%13.5E", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]]);
// 					}
// 					fprintf(cubeFile, "\n");
// 				}
// 				if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6 != 0) {
// 					for (int l = 0; l < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6; l++) {
// 						fprintf(cubeFile, "%13.5E", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]]);
// 					}
// 					fprintf(cubeFile, "\n");
// 				}
// 			}
// 		}
// 		
// 		fclose(cubeFile);
}

void CalcCurrentDensity() {
	int res[3];
	res[0] = GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0];
	res[1] = GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1];
	res[2] = GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2];
	
	C3DF<VORI_FLOAT> density;
	density.CopyFrom(GetTimeStep(0)->m_pVolumetricData);
	
	C3DF<VORI_FLOAT> densityGradient[3];
	int i;
	for (i = 0; i < 3; i++) {
		densityGradient[i].m_iRes[0] = density.m_iRes[0];
		densityGradient[i].m_iRes[1] = density.m_iRes[1];
		densityGradient[i].m_iRes[2] = density.m_iRes[2];
		densityGradient[i].m_fMinVal[0] = density.m_fMinVal[0];
		densityGradient[i].m_fMaxVal[0] = density.m_fMaxVal[0];
		densityGradient[i].m_fMinVal[1] = density.m_fMinVal[1];
		densityGradient[i].m_fMaxVal[1] = density.m_fMaxVal[1];
		densityGradient[i].m_fMinVal[2] = density.m_fMinVal[2];
		densityGradient[i].m_fMaxVal[2] = density.m_fMaxVal[2];
		densityGradient[i].Create();
	}
	
	for (i = 0; i < res[2]; i++) {
		int j;
		for (j = 0; j < res[1]; j++) {
			int k;
			for (k = 1; k < res[0] - 1; k++) {
				densityGradient[0].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (density.m_pBin[i * res[0] * res[1] + j * res[0] + k + 1] - density.m_pBin[i * res[0] * res[1] + j * res[0] + k - 1]) / (2.0 * g_fCubeXStep);
			}
			densityGradient[0].m_pBin[i * res[0] * res[1] + j * res[0]] = (density.m_pBin[i * res[0] * res[1] + j * res[0] + 1] - density.m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 1]) / (2.0 * g_fCubeXStep);
			densityGradient[0].m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 1] = (density.m_pBin[i * res[0] * res[1] + j * res[0]] - density.m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 2]) / (2.0 * g_fCubeXStep);
		}
	}
	
	for (i = 0; i < res[2]; i++) {
		int j;
		for (j = 1; j < res[1] - 1; j++) {
			int k;
			for (k = 0; k < res[0]; k++) {
				densityGradient[1].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (density.m_pBin[i * res[0] * res[1] + (j + 1) * res[0] + k] - density.m_pBin[i * res[0] * res[1] + (j - 1) * res[0] + k]) / (2.0 * g_fCubeYStep);
			}
		}
		int k;
		for (k = 0; k < res[0]; k++) {
			densityGradient[1].m_pBin[i * res[0] * res[1] + k] = (density.m_pBin[i * res[0] * res[1] + res[0] + k] - density.m_pBin[i * res[0] * res[1] + (res[1] - 1) * res[0] + k]) / (2.0 * g_fCubeYStep);
			densityGradient[1].m_pBin[i * res[0] * res[1] + (res[1] - 1) * res[0] + k] = (density.m_pBin[i * res[0] * res[1] + k] - density.m_pBin[i * res[0] * res[1] + (res[1] - 2) * res[0] + k]) / (2.0 * g_fCubeYStep);
		}
	}
	
	for (i = 1; i < res[2] - 1; i++) {
		int j;
		for (j = 0; j < res[1]; j++) {
			int k;
			for (k = 0; k < res[0]; k++) {
				densityGradient[2].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (density.m_pBin[(i + 1) * res[0] * res[1] + j * res[0] + k] - density.m_pBin[(i - 1) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
			}
		}
	}
	int j;
	for (j = 0; j < res[1]; j++) {
		int k;
		for (k = 0; k < res[0]; k++) {
			densityGradient[2].m_pBin[j * res[0] + k] = (density.m_pBin[res[0] * res[1] + j * res[0] + k] - density.m_pBin[(res[2] - 1) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
			densityGradient[2].m_pBin[(res[2] - 1) * res[0] * res[1] + j * res[0] + k] = (density.m_pBin[j * res[0] + k] - density.m_pBin[(res[2] - 2) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
		}
	}
// 	mprintf(GREEN, "%#.10g %#.10g %#.10g %#.10g\n", density.m_pBin[1942857], densityGradient[0].m_pBin[1942857], densityGradient[1].m_pBin[1942857], densityGradient[2].m_pBin[1942857]);
	
// 	FILE *cubeFile2 = fopen("test2.cube", "w");
// 	if (cubeFile2 == NULL) {
// 		printf("Could not open output file!\n");
// 		abort();
// 	}
// 	
// 	fprintf(cubeFile2, "\n\n");
// 	fprintf(cubeFile2, "%5lu %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
// 	fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
// 	fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
// 	fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
// 	for (int i = 0; i < (int)GetTimeStep(0)->m_iGesAtomCount; i++) {
// 		fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, GetTimeStep(0)->m_vaCoords[i][0] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][1] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][2] / 0.529177249f / 100.0f);
// 	}
// 	for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 		for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 			for (int k = 0; k < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6; k++) {
// 				for (int l = 0; l < 6; l++) {
// 					fprintf(cubeFile2, "%13.5E", densityGradient[2].m_pBin[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]]);
// 				}
// 				fprintf(cubeFile2, "\n");
// 			}
// 			if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6 != 0) {
// 				for (int l = 0; l < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6; l++) {
// 					fprintf(cubeFile2, "%13.5E", densityGradient[2].m_pBin[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]]);
// 				}
// 				fprintf(cubeFile2, "\n");
// 			}
// 		}
// 	}
// 	
// 	fclose(cubeFile2);
	
// 	CSparseMatrix testMatrix;
// 	testMatrix.setTest();
// 	double solution[9];
// 	for (int i = 0; i < 9; i++)
// 		solution[i] = 0.0;
// 	double rhs[9];
// 	rhs[0] = 1.0;
// 	rhs[1] = 2.0;
// 	rhs[2] = 1.0;
// 	rhs[3] = 2.0;
// 	rhs[4] = 1.0;
// 	rhs[5] = 2.0;
// 	rhs[6] = 1.0;
// 	rhs[7] = 2.0;
// 	rhs[8] = 1.0;
// 	
// 	CCurrentPDESolver::bicgstabl(4, &testMatrix, solution, rhs, 100, 1e-10, stdout);
// 	
// 	return;
	
	for (i = 0; i < res[0] * res[1] * res[2]; i++) {
		density.m_pBin[i] += g_fBackgroundDensity;
	}
	
	CSparseMatrix pdeMatrix;
	CCurrentPDEDiscretizer::discretize(&pdeMatrix, density, densityGradient[0], densityGradient[1], densityGradient[2]);
	
	static double *pdeSolution = NULL;
	if (pdeSolution == NULL) {
		try { pdeSolution = new double[res[0] * res[1] * res[2]]; } catch(...) { pdeSolution = NULL; }
		if (pdeSolution == NULL) NewException((double)res[0] * res[1] * res[2] * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
	}
	
	static bool first = true;
	static double thresh;
	if (first) {
		memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
		thresh = g_fPDEConvThresh * CCurrentPDESolver::calcResidual(&pdeMatrix, pdeSolution, GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin);
		first = false;
	}
	
	if (g_fPDESolverInfoFile != NULL) {
		fprintf(g_fPDESolverInfoFile, "Step %lu\n", g_iSteps - 1);
	}
	double thresh2 = thresh;
	if (!CCurrentPDESolver::bicgstabl(4, &pdeMatrix, pdeSolution, GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, g_fPDESolverInfoFile)) {
		memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
		if (g_fPDESolverInfoFile != NULL) {
			fprintf(g_fPDESolverInfoFile, "Resetting solution guess\n");
		}
		thresh2 = thresh;
		if (!CCurrentPDESolver::bicgstabl(4, &pdeMatrix, pdeSolution, GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, g_fPDESolverInfoFile)) {
			memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
			if (g_fPDESolverInfoFile != NULL) {
				fprintf(g_fPDESolverInfoFile, "Trying different threshold\n");
			}
			thresh2 *= 1.05;
			if (!CCurrentPDESolver::bicgstabl(4, &pdeMatrix, pdeSolution, GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, g_fPDESolverInfoFile)) {
				memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
				if (g_fPDESolverInfoFile != NULL) {
					fprintf(g_fPDESolverInfoFile, "Severe convergence problem! Resetting solution guess\n");
				}
			}
		}
	}
	if (g_fPDESolverInfoFile != NULL) {
		fprintf(g_fPDESolverInfoFile, "\n");
	}
	
	GetTimeStep(0)->m_pCurrentDensity->SetSize(3 * res[0] * res[1] * res[2]);
	
// 	fftwf_complex *fft_data;
// 	fft_data = fftwf_alloc_complex(res[0] * res[1] * res[2]);
// 	fftwf_plan plan1, plan2;
// 	plan1 = fftwf_plan_dft_3d(res[2], res[1], res[0], fft_data, fft_data, FFTW_FORWARD, FFTW_MEASURE);
// 	plan2 = fftwf_plan_dft_3d(res[2], res[1], res[0], fft_data, fft_data, FFTW_BACKWARD, FFTW_MEASURE);
// 	
// 	int i;
// 	for (i = 0; i < res[0] * res[1] * res[2]; i++) {
// 		fft_data[i][0] = GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i];
// 		fft_data[i][1] = 0.0f;
// 	}
// 	fftwf_execute(plan1);
// 	for (i = 0; i < res[2]; i++) {
// 		int j;
// 		for (j = 0; j < res[1]; j++) {
// 			int k;
// 			for (k = 0; k < res[0]; k++) {
// 				if (i == 0 && j == 0 && k == 0) {
// 					fft_data[i * res[1] * res[0] + j * res[0] + k][0] = 0.0f;
// 					fft_data[i * res[1] * res[0] + j * res[0] + k][1] = 0.0f;
// 					continue;
// 				}
// 				float f = (2.0f * cosf(2.0f * (float)Pi / res[2] * i) - 2.0f) / g_fCubeXStep / g_fCubeXStep + (2.0f * cosf(2.0f * (float)Pi / res[1] * j) - 2.0f) / g_fCubeYStep / g_fCubeYStep + (2.0f * cosf(2.0f * (float)Pi / res[0] * k) - 2.0f) / g_fCubeZStep / g_fCubeZStep;
// 				fft_data[i * res[1] * res[0] + j * res[0] + k][0] /= f;
// 				fft_data[i * res[1] * res[0] + j * res[0] + k][1] /= f;
// 			}
// 		}
// 	}
// 	fftwf_execute(plan2);
// 	for (i = 0; i < res[0] * res[1] * res[2]; i++) {
// 		fft_data[i][0] /= res[0] * res[1] * res[2];
// 		fft_data[i][1] /= res[0] * res[1] * res[2];
// 	}
	
// 		FILE *cubeFile2 = fopen("test2.cube", "w");
// 		if (cubeFile2 == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		fprintf(cubeFile2, "\n\n");
// 		fprintf(cubeFile2, "%5lu %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
// 		fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
// 		fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
// 		fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
// 		for (int i = 0; i < (int)GetTimeStep(0)->m_iGesAtomCount; i++) {
// 			fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, GetTimeStep(0)->m_vaCoords[i][0] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][1] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][2] / 0.529177249f / 100.0f);
// 		}
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				for (int k = 0; k < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6; k++) {
// 					for (int l = 0; l < 6; l++) {
// 						fprintf(cubeFile2, "%13.5E", fft_data[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]][0]);
// 					}
// 					fprintf(cubeFile2, "\n");
// 				}
// 				if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6 != 0) {
// 					for (int l = 0; l < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6; l++) {
// 						fprintf(cubeFile2, "%13.5E", fft_data[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]][0]);
// 					}
// 					fprintf(cubeFile2, "\n");
// 				}
// 			}
// 		}
// 		
// 		fclose(cubeFile2);
// 	
// 		FILE *testFile2 = fopen("test2.dat", "w");
// 		if (testFile2 == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				fprintf(testFile2, "%.6f %.6f %.14f\n", i * g_fCubeXStep, j * g_fCubeYStep, fft_data[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 2 * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]][0]);
// 			}
// 			fprintf(testFile2, "\n");
// 		}
// 		
// 		fclose(testFile2);
		
// 	for (i = 0; i < res[2]; i++) {
// 		int j;
// 		for (j = 0; j < res[1]; j++) {
// 			int k;
// 			for (k = 1; k < res[0] - 1; k++) {
// 				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3) = (fft_data[i * res[1] * res[0] + j * res[0] + k + 1][0] - fft_data[i * res[1] * res[0] + j * res[0] + k - 1][0]) / g_fCubeXStep;
// 			}
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3) = (fft_data[i * res[1] * res[0] + j * res[0] + 1][0] - fft_data[i * res[1] * res[0] + j * res[0] + res[0] - 1][0]) / g_fCubeXStep;
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + (res[0] - 1) * 3) = (fft_data[i * res[1] * res[0] + j * res[0]][0] - fft_data[i * res[1] * res[0] + j * res[0] + res[0] - 2][0]) / g_fCubeXStep;
// 		}
// 	}
// 	for (i = 0; i < res[2]; i++) {
// 		int j;
// 		for (j = 1; j < res[1] - 1; j++) {
// 			int k;
// 			for (k = 0; k < res[0]; k++) {
// 				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 1) = (fft_data[i * res[1] * res[0] + (j + 1) * res[0] + k][0] - fft_data[i * res[1] * res[0] + (j - 1) * res[0] + k][0]) / g_fCubeYStep;
// 			}
// 		}
// 		int k;
// 		for (k = 0; k < res[0]; k++) {
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + k * 3 + 1) = (fft_data[i * res[1] * res[0] + res[0] + k][0] - fft_data[i * res[1] * res[0] + (res[1] - 1) * res[0] + k][0]) / g_fCubeYStep;
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + (res[1] - 1) * res[0] * 3 + k * 3 + 1) = (fft_data[i * res[1] * res[0] + k][0] - fft_data[i * res[1] * res[0] + (res[1] - 2) * res[0] + k][0]) / g_fCubeYStep;
// 		}
// 	}
// 	for (i = 1; i < res[2] - 1; i++) {
// 		int j;
// 		for (j = 0; j < res[1]; j++) {
// 			int k;
// 			for (k = 0; k < res[0]; k++) {
// 				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (fft_data[(i + 1) * res[1] * res[0] + j * res[0] + k][0] - fft_data[(i - 1) * res[1] * res[0] + j * res[0] + k][0]) / g_fCubeZStep;
// 			}
// 		}
// 	}
// 	int j;
// 	for (j = 0; j < res[1]; j++) {
// 		int k;
// 		for (k = 0; k < res[0]; k++) {
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(j * res[0] * 3 + k * 3 + 2) = (fft_data[res[1] * res[0] + j * res[0] + k][0] - fft_data[(res[2] - 1) * res[1] * res[0] + j * res[0] + k][0]) / g_fCubeZStep;
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt((res[2] - 1) * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (fft_data[j * res[0] + k][0] - fft_data[(res[2] - 2) * res[1] * res[0] + j * res[0] + k][0]) / g_fCubeZStep;
// 		}
// 	}
// 	
// 	fftwf_destroy_plan(plan1);
// 	fftwf_destroy_plan(plan2);
// 	fftwf_free(fft_data);
	
	for (i = 0; i < res[2]; i++) {
		int j;
		for (j = 0; j < res[1]; j++) {
			int k;
			for (k = 1; k < res[0] - 1; k++) {
				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3) = (pdeSolution[i * res[1] * res[0] + j * res[0] + k + 1] - pdeSolution[i * res[1] * res[0] + j * res[0] + k - 1]) / (2.0 * g_fCubeXStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
			}
			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3) = (pdeSolution[i * res[1] * res[0] + j * res[0] + 1] - pdeSolution[i * res[1] * res[0] + j * res[0] + res[0] - 1]) / (2.0 * g_fCubeXStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0]];
			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + (res[0] - 1) * 3) = (pdeSolution[i * res[1] * res[0] + j * res[0]] - pdeSolution[i * res[1] * res[0] + j * res[0] + res[0] - 2]) / (2.0 * g_fCubeXStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + res[0] - 1];
		}
	}
	for (i = 0; i < res[2]; i++) {
		int j;
		for (j = 1; j < res[1] - 1; j++) {
			int k;
			for (k = 0; k < res[0]; k++) {
				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 1) = (pdeSolution[i * res[1] * res[0] + (j + 1) * res[0] + k] - pdeSolution[i * res[1] * res[0] + (j - 1) * res[0] + k]) / (2.0 * g_fCubeYStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
			}
		}
		int k;
		for (k = 0; k < res[0]; k++) {
			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + k * 3 + 1) = (pdeSolution[i * res[1] * res[0] + res[0] + k] - pdeSolution[i * res[1] * res[0] + (res[1] - 1) * res[0] + k]) / (2.0 * g_fCubeYStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + k];
			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + (res[1] - 1) * res[0] * 3 + k * 3 + 1) = (pdeSolution[i * res[1] * res[0] + k] - pdeSolution[i * res[1] * res[0] + (res[1] - 2) * res[0] + k]) / (2.0 * g_fCubeYStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + (res[1] - 1) * res[0] + k];
		}
	}
	for (i = 1; i < res[2] - 1; i++) {
		int j;
		for (j = 0; j < res[1]; j++) {
			int k;
			for (k = 0; k < res[0]; k++) {
				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (pdeSolution[(i + 1) * res[1] * res[0] + j * res[0] + k] - pdeSolution[(i - 1) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
			}
		}
	}
	for (j = 0; j < res[1]; j++) {
		int k;
		for (k = 0; k < res[0]; k++) {
			GetTimeStep(0)->m_pCurrentDensity->GetAt(j * res[0] * 3 + k * 3 + 2) = (pdeSolution[res[1] * res[0] + j * res[0] + k] - pdeSolution[(res[2] - 1) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[j * res[0] + k];
			GetTimeStep(0)->m_pCurrentDensity->GetAt((res[2] - 1) * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (pdeSolution[j * res[0] + k] - pdeSolution[(res[2] - 2) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[(res[2] - 1) * res[1] * res[0] + j * res[0] + k];
		}
	}
// 	mprintf(GREEN, "%#.10g\n", GetTimeStep(0)->m_pCurrentDensity->GetAt(0));
	
// 		FILE *cubeFile = fopen("test.cube", "w");
// 		if (cubeFile == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		fprintf(cubeFile, "\n\n");
// 		fprintf(cubeFile, "%5lu %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
// 		for (int i = 0; i < (int)GetTimeStep(0)->m_iGesAtomCount; i++) {
// 			fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, GetTimeStep(0)->m_vaCoords[i][0] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][1] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][2] / 0.529177249f / 100.0f);
// 		}
// 		double sum = 0.0;
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				for (int k = 0; k < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6; k++) {
// 					for (int l = 0; l < 6; l++) {
// 						fprintf(cubeFile, "%13.5E", GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * 3));
// 						sum += GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * 3);
// 					}
// 					fprintf(cubeFile, "\n");
// 				}
// 				if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6 != 0) {
// 					for (int l = 0; l < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6; l++) {
// 						fprintf(cubeFile, "%13.5E", GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * 3));
// 						sum += GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * 3);
// 					}
// 					fprintf(cubeFile, "\n");
// 				}
// 			}
// 		}
// 		mprintf(GREEN, "%g\n", sum);
// 		
// 		fclose(cubeFile);
// 		
// 		FILE *testFile = fopen("test.dat", "w");
// 		if (testFile == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				fprintf(testFile, "%.6f %.6f %.14f %.14f\n", i * g_fCubeXStep, j * g_fCubeYStep, GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 2 * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3), GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 2 * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + 1));
// 			}
// 			fprintf(testFile, "\n");
// 		}
// 		
// 		fclose(testFile);
}

void CalcForces()
{
	BTIN;
	int z;
	float f, imf, ilmf;

	g_fLMidForce = 0;
	imf = g_fMaxForce;
	ilmf = g_fLMaxForce;
	GetTimeStep(0)->m_vaForces.SetSize(g_iGesVirtAtomCount);
	for (z=0;z<g_iGesVirtAtomCount;z++)
	{
		GetTimeStep(0)->m_vaForces[z] = (GetTimeStep(-1)->m_vaCoords[z] - 2*GetTimeStep(0)->m_vaCoords[z] + GetTimeStep(1)->m_vaCoords[z])/(float)pow((g_fTimestepLength/1000.0),2);
		f = GetTimeStep(0)->m_vaForces[z].GetLength();
// 		mprintf(GREEN, "%g\n", f);
		g_fLMidForce += f;
		if (f > imf)
			imf = f;
		if (f > ilmf)
			ilmf = f;
	}
	g_fLMidForce /= g_iGesVirtAtomCount;
	g_fLMaxForce = ilmf;
	if (imf < g_fUnsteadyLimit)
	{
		g_fMaxForce = imf;
	} else if (!g_bWarnUnsteady)
	{
		g_bWarnUnsteady = true;
		mprintf("Warning: Discontinuity in time step %d (maximum acceleration = %f > %f pm/ps^2).\n",(int)g_iSteps,imf,g_fUnsteadyLimit);
	}
	BTOUT;
}


/*float AtomMass(char *s)
{
	BTIN;
	int z;
	for (z=0;z<g_oaElements.GetSize();z++)
		if (mystricmp(s,((CElement*)g_oaElements[z])->m_sLabel)==0)
		{
			BTOUT;
			return ((CElement*)g_oaElements[z])->m_fMass;
		}
	eprintf("Atom \"%s\" not found.\n",s);
	BTOUT;
	return -1.0f;
}

int AtomOrd(char *s)
{
	BTIN;
	int z;
	for (z=0;z<g_oaElements.GetSize();z++)
		if (mystricmp(s,((CElement*)g_oaElements[z])->m_sLabel)==0)
		{
			BTOUT;
			return ((CElement*)g_oaElements[z])->m_iOrd;
		}
	eprintf("Atom \"%s\" not found.\n",s);
	BTOUT;
	return -1;
}

float AtomRadius(char *s)
{
	BTIN;
	int z;
	for (z=0;z<g_oaElements.GetSize();z++)
		if (mystricmp(s,((CElement*)g_oaElements[z])->m_sLabel)==0)
		{
			BTOUT;
			return ((CElement*)g_oaElements[z])->m_fRadius;
		}
	mprintf("No Atom Radius found for Atom \"%s\".\n",s);
	BTOUT;
	return 0;
}*/


CElement* FindElement(const char *s, bool quiet)
{
	BTIN;
	int z;

	for (z=0;z<g_oaElements.GetSize();z++)
	{
		if (mystricmp(s,((CElement*)g_oaElements[z])->m_sLabel)==0)
		{
			BTOUT;
			return (CElement*)g_oaElements[z];
		}
	}
	if (!quiet)
		eprintf("No element data found for atom \"%s\".\n",s);
	g_bUnknownElements = true;
	BTOUT;
	return NULL;
}


float GuessBoxSize()
{
	BTIN;
	int z;
	float m, gm;

	gm = 0;
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		if (z == g_iVirtAtomType)
			continue;
		m = ((CAtom*)g_oaAtoms[z])->m_pElement->m_fMass;
/*		if (m < 0) // Dieses Atom ist nicht in der Liste
		{
			if (z != g_iVirtAtomType)
				mprintf("Atom \"%s\" nicht in der Liste.\n",((CAtom*)g_oaAtoms[z])->m_sName);
			continue;
		}*/
		gm += m * ((CAtom*)g_oaAtoms[z])->m_iCount;
	}
	m = (float)pow(gm/6.022f*10,0.33333333f) * 100.0f;  // In pm
	BTOUT;
	return m;
}


void strtolower(char *s)
{
	BTIN;
	char *p;
	p = s;
	while (*p != 0)
	{
		if ((*p >= 'A') && (*p <= 'Z'))
			*p += 32;
		p++;
	}
	BTOUT;
}


void SortAtoms()
{
	BTIN;
	int z, z2;
	int i;
	CAtom *a;

	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		i = z;
		for (z2=z+1;z2<g_oaAtoms.GetSize();z2++)
			if (strcmp(((CAtom*)g_oaAtoms[i])->m_sName,((CAtom*)g_oaAtoms[z2])->m_sName) > 0)
				i = z2;
		if (i != z)
		{
			a = (CAtom*)g_oaAtoms[z];
			g_oaAtoms[z] = g_oaAtoms[i];
			g_oaAtoms[i] = a;
		}
	}
	for (z=0;z<g_oaAtoms.GetSize();z++)
		((CAtom*)g_oaAtoms[z])->m_iIndex = z;
	BTOUT;
}


bool SetAnalysis(const char *s)
{
	BTIN;

	if (mystricmp(s, "chdf") == 0) {
		g_bCHDF = true;
		return true;
	}
	if(mystricmp(s, "chi") == 0) {
		g_bChiral = true;
		return true;
	}
	if(mystricmp(s, "eck") == 0) {
		g_bEckartTransform = true;
		return true;
	}
	if (mystricmp(s, "pol") == 0) {
		g_bSetUpPolarizabilityCalc = true;
		return true;
	}

	if(mystricmp(s, "swan") == 0) {
		g_bSortWannier = true;
		return true;
	}
	if(mystricmp(s, "vcd") == 0) {
		g_bVCD = true;
		return true;
	}
	if (mystricmp(s, "drst") == 0) {
		g_bDipoleRestart = true;
		return true;
	}
	if (mystricmp(s, "mrst") == 0) {
		g_bMagneticDipoleRestart = true;
		return true;
	}

	if(mystricmp(s, "nc") == 0) {
		g_bNormalCoordinate = true;
		return true;
	}
	if (mystricmp(s,"dprof")==0)
	{
/*		if ((!g_bPeriodicX) || (!g_bPeriodicY) || (!g_bPeriodicZ))
		{
			eprintf("\n    Error: Requires XYZ-periodic box!\n\n");
			BTOUT;
			return false;
		}*/
		g_bPDF = true;
		BTOUT;
		return true;
	}

	if (mystricmp(s,"vori")==0)
	{
		g_bTegri = true;
		BTOUT;
		return true;
	}
	
	if (mystricmp(s,"doma")==0)
	{
		g_bDomA = true;
		BTOUT;
		return true;
	}

	if (mystricmp(s,"plproj")==0)
	{
		g_bPlProj = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"pldf")==0)
	{
		g_bPlDF = true;
		BTOUT;
		return true;
	}


	if (mystricmp(s,"lidf")==0)
	{
		g_bLiDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"power")==0)
	{
// 		g_bACF = true;
// 		g_bPowerSpec = true;
		g_bPower = true;
		BTOUT;
		return true;
	}


	if (mystricmp(s,"ir")==0)
	{
// 		g_bRDyn = true;
// 		g_bIRSpec = true;
		g_bIR = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"reg")==0)
	{
		g_bRegionAnalysis = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"raman")==0)
	{
		g_bRaman = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"sfac")==0)
	{
		g_bSFac = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"dens")==0)
	{
		if ((!g_bPeriodicX) || (!g_bPeriodicY) || (!g_bPeriodicZ))
		{
			eprintf("\n    Error: Requires XYZ-periodic box!\n\n");
			BTOUT;
			return false;
		}
		g_bDens = true;
		BTOUT;
		return true;
	}


	if (mystricmp(s,"cond")==0)
	{
		g_bCond = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"sdf")==0)
	{
		g_bSDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"psdf")==0)
	{
		g_bRevSDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"rdf")==0)
	{
		g_bRDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"vdf")==0)
	{
		g_bVDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"fdf")==0)
	{
		g_bFDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"adf")==0)
	{
		g_bADF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"ddf")==0)
	{
		g_bDDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"cdf")==0)
	{
		g_bCDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"msd")==0)
	{
		g_bMSD = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"aggr")==0)
	{
		g_bAggregation = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"nbex")==0)
	{
		g_bNbExchange = true;
		BTOUT;
		return true;
	}
/*	if (mystricmp(s,"ddisp")==0)
	{
		g_bDDisp = true;
		g_bAggregation = true;
		BTOUT;
		return true;
	}*/
	if (mystricmp(s,"acf")==0)
	{
		g_bACF = true;
		BTOUT;
		return true;
	}


	if (mystricmp(s,"nbx")==0)
	{
		g_bNbExchange = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"rdyn")==0)
	{
		g_bRDyn = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"cut")==0)
	{
		g_bCutCluster = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"avg")==0)
	{
		g_bAvg = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"env")==0)
	{
		g_bSaveRefEnv = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"proc")==0)
	{
		g_bSaveJustTraj = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"vhcf")==0)
	{
		g_bVHDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"nbh")==0)
	{
		g_bNbAnalysis = true;
		BTOUT;
		return true;
	}
/*	if (mystricmp(s,"vfd")==0)
	{
		g_bSaveVelForce = true;
		BTOUT;
		return true;
	}*/
	if (mystricmp(s,"dip")==0)
	{
		g_bDipDF = true;
		g_bDipole = true;
		BTOUT;
		return true;
	}

	if (mystricmp(s,"voro")==0)
	{
		if ((!g_bPeriodicX) || (!g_bPeriodicY) || (!g_bPeriodicZ))
		{
			eprintf("\n    Error: Basic Voronoi Analysis currently requires XYZ-periodic box.\n\n");
			BTOUT;
			return false;
		}
		g_bVoro = true;
		BTOUT;
		return true;
	}


	BTOUT;
	return false;
}


bool ParseFunctions(const char *s)
{
	BTIN;
	const char *cp;
	char *p, *q;
	char buf[256];
	int z;
	CAnalysis *a;

	
	g_bSFac = false;
	g_bVoid = false;
	g_bTegri = false;
	g_bVoro = false;
	g_bDomA = false;
	g_bPlProj = false;
	g_bPlDF = false;
	g_bLiDF = false;
	g_bDens = false;
	g_bCond = false;
	g_bSDF = false;
	g_bRDF = false;
	g_bVDF = false;
	g_bFDF = false;
	g_bADF = false;
	g_bAvg = false;
	g_bSaveVelForce = false;
	g_bSaveRefEnv = false;
	g_bRefEnvCenter = false;
	g_bSaveJustTraj = false;
	g_bVHDF = false;
	g_bCutCluster = false;
	g_bNbAnalysis = false;
	g_bVACF = false;
	g_bDipACF = false;
	g_bDDF = false;
	g_bMSD = false;
	g_bDipole = false;
	g_bDipDF = false;
	g_bCDF = false;
	g_bAggregation = false;
	g_bDDisp = false; 
	g_bDACF = false; 
	g_bDLDF = false; 
	g_bRDyn = false;
	g_bNbExchange = false;
	g_bRaman = false;
	g_bIRSpec = false;
	g_bPowerSpec = false;
	g_bRegionAnalysis = false;

	cp = s;
	q = buf;

	while (*cp != 0)
	{
		if (*cp == ' ')
			cp++;
		else {
			*q = *cp;
			cp++;
			q++;
		}
	}
	*q = 0;

	p = buf;
	
	while (true)
	{
		q = strchr(p,',');
		if (q != NULL)
			*q = 0;

		for (z=0;z<g_oaAnalyses.GetSize();z++)
		{
			a = (CAnalysis*)g_oaAnalyses[z];
			if (mystricmp(p,a->m_sAbbrev)==0)
			{
				mprintf(WHITE," %-7s",a->m_sAbbrev);
				mprintf("- %s\n",a->m_sName);
				if (!SetAnalysis(p))
				{
					eprintf("Function cannot be applied to this system: \"%s\"\n",p);
					return false;
				}
				goto _done;
			}
		}
		eprintf("Unknown Function: \"%s\"\n",p);
		return false;
_done:
		if ((q == NULL) || (*(q+1) == 0))
		{
			BTOUT;
			return true;
		}
		p = q+1;
	}
}


bool ParsePeriodic(const char *s)
{
	const char *p;

	BTIN;
	if (strlen(s) == 0)
	{
		g_bPeriodicX = true;
		g_bPeriodicY = true;
		g_bPeriodicZ = true;
		g_bPeriodic = true;
		BTOUT; 
		return true;
	}
	g_bPeriodicX = false;
	g_bPeriodicY = false;
	g_bPeriodicZ = false;

	p = s;
	while (*p != 0)
	{
		if ((*p == 'x') || (*p == 'X'))
			g_bPeriodicX = true;
		else if ((*p == 'y') || (*p == 'Y'))
			g_bPeriodicY = true;
		else if ((*p == 'z') || (*p == 'Z'))
			g_bPeriodicZ = true;
		else if (*p != '0')
		{
			eprintf("Wrong input.\n");
			BTOUT;
			return false;
		}
		p++;
	}
	g_bPeriodic = g_bPeriodicX || g_bPeriodicY || g_bPeriodicZ;
	BTOUT;
	return true;
}


void WriteHeader()
{
	struct tm *today;
	time_t ltime;
	char buf[64];
	unsigned long l, *lp;
	unsigned char *ca, *cb, *cc, *cd;

	BTIN;
	mprintf("\n");  /* http://patorjk.com/software/taag/ "Big Money-SW" */

/*	mprintf(" ________          ______   __     __  __           \n");
	mprintf("/        |        /      \\ /  |   /  |/  |          \n");
	mprintf("########/______  /######  |## |   ## |##/   _______ \n");
	mprintf("   ## | /      \\ ## |__## |## |   ## |/  | /       |\n");
	mprintf("   ## |/######  |##    ## |##  \\ /##/ ## |/#######/ \n");
	mprintf("   ## |## |  ##/ ######## | ##  /##/  ## |##      \\ \n");
	mprintf("   ## |## |      ## |  ## |  ## ##/   ## | ######  |\n");
	mprintf("   ## |## |      ## |  ## |   ###/    ## |/     ##/ \n");
	mprintf("   ##/ ##/       ##/   ##/     #/     ##/ #######/  \n");*/

/*	mprintf("   ________                             __\n");
	mprintf("  /        |                           /  |\n");
	mprintf("  ########/______   ______   __     __ ##/   _______\n");
	mprintf("     ## | /      \\ /      \\ /  \\   /  |/  | /       |\n");
	mprintf("     ## |/######  |######  |##  \\ /##/ ## |/#######/\n");
	mprintf("     ## |## |  ##/ /    ## | ##  /##/  ## |##      \\\n");
	mprintf("     ## |## |     /####### |  ## ##/   ## | ######  |\n");
	mprintf("     ## |## |     ##    ## |   ###/    ## |/     ##/\n");
	mprintf("     ##/ ##/       #######/     #/     ##/ #######/\n");*/

	mprintf(YELLOW,"   ________                                 __\n");
	mprintf(YELLOW,"  /        |                               /  |\n");
	mprintf(YELLOW,"  ########/ ______    ______    __     __  ##/    _______\n");
	mprintf(YELLOW,"     ## |  /      \\  /      \\  /  \\   /  | /  |  /       |\n");
	mprintf(YELLOW,"     ## | /######  | ######  | ##  \\ /##/  ## | /#######/\n");
	mprintf(YELLOW,"     ## | ## |  ##/  /    ## |  ##  /##/   ## | ##      \\\n");
	mprintf(YELLOW,"     ## | ## |      /####### |   ## ##/    ## |  ######  |\n");
	mprintf(YELLOW,"     ## | ## |      ##    ## |    ###/     ## | /     ##/\n");
	mprintf(YELLOW,"     ##/  ##/        #######/      #/      ##/  #######/\n");

	mprintf(WHITE,"\n     TRajectory Analyzer and VISualizer  -  Open-source freeware under GNU GPL v3\n\n");
	mprintf("");
	mprintf("     Copyright (c) Martin Brehm      (2009-2016)\n");
	mprintf("                   Martin Thomas     (2012-2016)\n");
	mprintf("                   Barbara Kirchner  (2009-2016)\n");
	mprintf("                   University of Leipzig / University of Bonn.\n\n");
	mprintf(YELLOW,"     http://www.travis-analyzer.de\n\n");
//	mprintf("     Open-source freeware; Licensed under the GNU General Public License v3.\n\n");
	mprintf("     Please cite:\n");
	mprintf(WHITE,"     M. Brehm and B. Kirchner, J. Chem. Inf. Model. 2011, 51 (8), pp 2007-2023.\n\n");
	mprintf("     There is absolutely no warranty on any results obtained from TRAVIS.\n\n");


	time(&ltime);
	today = localtime(&ltime);
	strcpy(buf,asctime(today));
	buf[strlen(buf)-1] = 0;
	mprintf(WHITE,"  #  ");
	if (g_sHostName != NULL)
		mprintf("Running on %s at %s",g_sHostName,buf);
			else mprintf("Running at %s",buf);

#ifdef TARGET_LINUX
	mprintf(" (PID %d).\n",getpid());
#else
	mprintf(".\n");
#endif

	if (g_sWorkingDir != NULL)
	{
		mprintf(WHITE,"  #  ");
		mprintf("Running in %s\n",g_sWorkingDir);
	}

	mprintf(WHITE,"  #  ");
	mprintf("Source code version: ");
	mprintf("%s.\n",SOURCE_VERSION);
	mprintf(WHITE,"  #  ");
	mprintf("Compiled at ");
	mprintf("%s %s.\n",__DATE__,__TIME__);
#ifdef __VERSION__
	mprintf(WHITE,"  #  ");
	mprintf("Compiler version: %s\n",__VERSION__);
#endif

#ifdef TARGET_WINDOWS
	mprintf(WHITE,"  #  ");
	mprintf("Target platform: Windows\n");
#elif defined(TARGET_LINUX)
	mprintf(WHITE,"  #  ");
	mprintf("Target platform: Linux\n");
#else
	mprintf(WHITE,"  #  ");
	mprintf("Target platform: Generic\n");
#endif
	mprintf(WHITE,"  #  ");
	mprintf("Compile flags: ");
#ifdef DEBUG_BACKTRACE
	mprintf("DEBUG_BACKTRACE ");
#endif
#ifdef DEBUG_EXTENDED_BACKTRACE
	mprintf("DEBUG_EXTENDED_BACKTRACE ");
#endif
#ifdef USE_FFTW
	mprintf("USE_FFTW ");
#endif
#ifdef DEBUG_ARRAYS
	mprintf("DEBUG_ARRAYS ");
#endif
#ifdef DEBUG_COBARRAY
	mprintf("DEBUG_COBARRAY ");
#endif
#ifdef DEBUG_CSTRING
	mprintf("DEBUG_CSTRING ");
#endif
#ifdef DEBUG_CVEC3ARRAY
	mprintf("DEBUG_CVEC3ARRAY ");
#endif
#ifdef DEBUG_DATABASE
	mprintf("DEBUG_DATABASE ");
#endif
	mprintf("\n");

	l = 0xA0B0C0D0;
	lp = &l;
	ca = (unsigned char*)lp;
	cb = ca+1;
	cc = ca+2;
	cd = ca+3;
	mprintf(WHITE,"  #  ");
	mprintf("Machine: int=%db, long=%db, addr=%db, 0xA0B0C0D0=%02X,%02X,%02X,%02X.\n",sizeof(int),sizeof(long),sizeof(void*),*ca,*cb,*cc,*cd);

	if (g_sHomeDir != NULL)
	{
		mprintf(WHITE,"  #  ");
		mprintf("User home: %s\n",g_sHomeDir);
	}

	mprintf(WHITE,"  #  ");
	mprintf("Exe path: %s\n",g_sExeName);
	mprintf(WHITE,"  #  ");

	if ((!IsTTY(stdin)) || (g_sInputFile != NULL))
	{
		if (g_sInputFile != NULL)
			mprintf("Input from %s, ",g_sInputFile);
				else mprintf("Input is redirected, ");
	} else mprintf("Input from terminal, ");

	if (IsTTY(stdout))
		mprintf("Output to terminal\n");
			else mprintf("Output is redirected\n");

	if (g_bStreamInput)
	{
		mprintf(WHITE,"  #  ");
		mprintf("Input trajectory treated as stream\n");
	}

	mprintf("\n");
	mprintf(" >>> Please use a color scheme with dark background or specify \"-nocolor\"! <<<\n\n");
//	mprintf("\nSome general hints:\n  - If you don't Enter something and just press RETURN,\n    the Value in [Square Brackets] will be used as default.\n");
//	mprintf("  - If your input was wrong and you want to jump back,\n    please enter \"$\".\n\n");
	BTOUT;
}


void CommandLineHelp()
{
	BTIN;
	mprintf(WHITE,"    List of supported command line options:\n\n");
	mprintf("      -p <file>       Loads position data from the specified trajectory file.\n");
	mprintf("                      The file format may be *.xyz, *.pdb, *.lmp (Lammps), HISTORY (DLPOLY), or *.prmtop/*.mdcrd (Amber).\n");
	mprintf("      -i <file>       Reads input from the specified text file.\n\n");
	mprintf("      -config <file>  Load the specified configuration file.\n");
	mprintf("      -stream         Treats input trajectory as a stream (e.g. named pipe): No fseek, etc.\n");
	mprintf("      -showconf       Shows a tree structure of the configuration file.\n");
	mprintf("      -writeconf      Writes the default configuration file, including all defines values.\n\n");
	mprintf("      -verbose        Show detailed information about what's going on.\n");
	mprintf("      -nocolor        Executes TRAVIS in monochrome mode.\n");
	mprintf("      -dimcolor       Uses dim instead of bright colors.\n\n");
	mprintf("      -credits        Display a list of persons who contributed to TRAVIS.\n");
	mprintf("      -help, -?       Shows this help.\n");
	mprintf("\n");
	mprintf("    If only one argument is specified, it is assumed to be the name of a trajectory file.\n");
	mprintf("    If argument is specified at all, TRAVIS asks for the trajectory file to open.\n");
	BTOUT;
}


bool ParseArgs(int argc, const char *argv[])
{
	BTIN;
	const char *p;
	char SModeFlag[64];
	int z;
//	bool namegiven;

	UnBase64((unsigned char *)SModeFlag,(const unsigned char *)"LXNjaGlzcw==",12);

	g_bAsciiArt = false;
	z = 1;

//	g_sInputTraj[0] = 0;

//	g_sInputVel[0] = 0;
	g_sInputVel.sprintf("");

//	g_sInputForce[0] = 0;
	g_sInputForce.sprintf("");

//	g_sInputCtrl[0] = 0;
	g_sInputCtrl.sprintf("");

//	namegiven = false;

	while (z < argc)
	{
		if ((memcmp(argv[z],"-?",2)==0) || (memcmp(argv[z],"--?",3)==0) || (memcmp(argv[z],"-help",5)==0) || (memcmp(argv[z],"--help",6)==0))
		{
			CommandLineHelp();
			BTOUT;
			return false;
		}

		if (mystricmp(argv[z],"-p")==0)
		{
//			mprintf("Dateiname: \"%s\"\n",argv[z+1]);
			z++;

			try { g_sInputTraj = new char[strlen(argv[z])+1]; } catch(...) { g_sInputTraj = NULL; }
			if (g_sInputTraj == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			strcpy(g_sInputTraj,argv[z]);
			goto _argend;
		}

		if (mystricmp(argv[z],"-conf")==0)
		{
//			mprintf("Dateiname: \"%s\"\n",argv[z+1]);
			z++;

			try { g_sConfFile = new char[strlen(argv[z])+1]; } catch(...) { g_sConfFile = NULL; }
			if (g_sConfFile == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			strcpy(g_sConfFile,argv[z]);
			goto _argend;
		}

		if (mystricmp(argv[z],"-i")==0)
		{
//			mprintf("Dateiname: \"%s\"\n",argv[z+1]);
			z++;

/*			try { g_sInputFile = new char[strlen(argv[z])+1]; } catch(...) { g_sInputFile = NULL; }
			if (g_sInputFile == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(g_sInputFile,argv[z]);*/
			goto _argend;
		}

		if (mystricmp(argv[z],"-verbose")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-credits")==0)
		{
			goto _argend;
		}


		if (mystricmp(argv[z],"-showconf")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-stream")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-writeconf")==0)
		{
			goto _argend;
		}

/*		if (mystricmp(argv[z],"-v")==0)
		{
			z++;
			strcpy(g_sInputVel,argv[z]);
			goto _argend;
		}

		if (mystricmp(argv[z],"-f")==0)
		{
			z++;
			strcpy(g_sInputForce,argv[z]);
			goto _argend;
		}*/

		if ((mystricmp(argv[z],"-nocolor")==0) || (mystricmp(argv[z],"--nocolor")==0))
		{
			goto _argend;
		}

		if ((mystricmp(argv[z],"-dimcolor")==0) || (mystricmp(argv[z],"--dimcolor")==0))
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-a")==0)
		{
//			mprintf("AA ist an.\n");
			g_bAsciiArt = true;
			goto _argend;
		}

		if (mystricmp(argv[z],"-lsd")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-readdipole")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-sax")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],SModeFlag)==0)
		{
			g_bSMode = true;
			goto _argend;
		}

		if ((argc > 2) || (argv[z][0] == '-'))
		{
			eprintf("Unknown parameter: \"%s\".\n\n",argv[z]);
			CommandLineHelp();
			BTOUT;
			return false;
		}

		p = strrchr(argv[z],'.');
		if (p != NULL)
		{
			p++;
			if ((mystricmp(p,"xyz")==0) || (mystricmp(p,"pdb")==0) || (mystricmp(p,"mol2")==0) || (mystricmp(p,"lmp")==0) || (mystricmp(p,"HISTORY")==0) || (mystricmp(p,"prmtop")==0) || (mystricmp(p,"mdcrd")==0))
			{
				try { g_sInputTraj = new char[strlen(argv[z])+1]; } catch(...) { g_sInputTraj = NULL; }
				if (g_sInputTraj == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				strcpy(g_sInputTraj,argv[z]);
			} else
			{
				eprintf("Unknown parameter: \"%s\".\n\n",argv[z]);
				CommandLineHelp();
				BTOUT;
				return false;
			}

//			try { g_sInputTraj = new char[strlen(argv[z])+1]; } catch(...) { g_sInputTraj = NULL; }
//			if (g_sInputTraj == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

//			strcpy(g_sInputTraj,argv[z]);
		} else
		{

			if (mystricmp(&argv[z][strlen(argv[z])-7],"HISTORY")==0)
			{
				try { g_sInputTraj = new char[strlen(argv[z])+1]; } catch(...) { g_sInputTraj = NULL; }
				if (g_sInputTraj == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				strcpy(g_sInputTraj,argv[z]);
			} else
			{
				eprintf("Unknown parameter: \"%s\".\n\n",argv[z]);
				CommandLineHelp();
				BTOUT;
				return false;
			}

		}
		BTOUT;
		return true;

_argend:
		z++;
	}
	BTOUT;
	return true;
}


void ParsePassiveArgs(int argc, const char *argv[])
{
	BTIN;
	int z;

	g_bGlobalPsycho = false;

	z = 1;
	while (z < argc)
	{
		if (mystricmp(argv[z],"-credits")==0)
			g_bShowCredits = true;

		if (mystricmp(argv[z],"-lsd")==0)
			g_bGlobalPsycho = true;

		if (mystricmp(argv[z],"-readdipole")==0)
			g_bDipolGrimme = true;

		if (mystricmp(argv[z],"-sax")==0)
			g_bSaxonize = true;

		if (mystricmp(argv[z],"-verbose")==0)
			g_bVerbose = true;


		if (mystricmp(argv[z],"-showconf")==0)
			g_bShowConf = true;

		if (mystricmp(argv[z],"-stream")==0)
			g_bStreamInput = true;

		if (mystricmp(argv[z],"-writeconf")==0)
			g_bWriteConf = true;

		if ((mystricmp(argv[z],"-nocolor")==0) || (mystricmp(argv[z],"--nocolor")==0))
			g_bNoColor = true;

		if ((mystricmp(argv[z],"-dimcolor")==0) || (mystricmp(argv[z],"--dimcolor")==0))
			g_iColorIntensity = 2;

		if (mystricmp(argv[z],"-i")==0)
		{
			z++;

			try { g_sInputFile = new char[strlen(argv[z])+1]; } catch(...) { g_sInputFile = NULL; }
			if (g_sInputFile == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(g_sInputFile,argv[z]);
		}

		z++;
	}
	BTOUT;
}


void CreateDatabaseDefaults()
{
	g_pDatabase->AddString("/GLOBAL/TRAVIS_VERSION",SOURCE_VERSION);
	g_pDatabase->AddBool("/GLOBAL/SHOWCREDITS",false);

	g_pDatabase->AddBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA",true);
	g_pDatabase->AddBool("/PLOT2D/FORMATS/WRITE_GNUPLOT",true);
	g_pDatabase->AddBool("/PLOT2D/FORMATS/WRITE_TRIPLES",true);
	g_pDatabase->AddBool("/PLOT2D/FORMATS/WRITE_MATRIX",true);

	g_pDatabase->AddInt("/PLOT2D/DEFAULTS/BIN_RES",100);
	g_pDatabase->AddInt("/PLOT2D/DEFAULTS/IMAGE_RES",1000);
	g_pDatabase->AddFloat("/PLOT2D/DEFAULTS/PLOT_EXP",0.5);
	g_pDatabase->AddInt("/PLOT2D/DEFAULTS/CONTOUR_LINES",30);

	g_pDatabase->AddBool("/PLOT3D/FORMATS/WRITE_CUBE",true);
	g_pDatabase->AddBool("/PLOT3D/FORMATS/WRITE_PLT",true);
	g_pDatabase->AddInt("/PLOT3D/DEFAULTS/BIN_RES",100);

/*	g_pDatabase->AddString("/BLA/BLUBB/PLOEPP/STRING1","String 1 Content");
	g_pDatabase->AddInt("/BLA/BLUBB/PLOEPP/INT1",123456);
	g_pDatabase->AddFloat("/BLA/BLUBB/PLOEPP/FLOAT1",1.23456);
	g_pDatabase->AddBool("/BLA/BLUBB/PLOEPP/BOOL1",true);*/
	//g_pDatabase->AddInt("/BLA/BLUBB/PLOEPP/STRING1",17);
}


void LoadSettings()
{
//	char buf[256];
	CxString buf;
	char sep;
	FILE *a;

#ifdef TARGET_WINDOWS
	sep = '\\';
#else
	sep = '/';
#endif

	if (g_sConfFile != NULL)
	{
		mprintf(WHITE,"    Custom configuration file specified: %s\n",g_sConfFile);
		if (FileExist(g_sConfFile))
		{
			mprintf("    Loading configuration from %s ...\n",g_sConfFile);
			g_pDatabase->ParseInputFile(g_sConfFile);
			return;
		} else
		{
			eprintf("    File does not exist or cannot be read.\n");
			abort();
		}
	}

//	sprintf(buf,"%s%ctravis.conf",g_sWorkingDir,sep);
	buf.sprintf("%s%ctravis.conf",g_sWorkingDir,sep);
	if (FileExist(buf))
	{
		mprintf("    Loading configuration from %s ...\n",(const char*)buf);
		g_pDatabase->ParseInputFile(buf);
		return;
	}

//	sprintf(buf,"%s%c.travis.conf",g_sWorkingDir,sep);
	buf.sprintf("%s%c.travis.conf",g_sWorkingDir,sep);
	if (FileExist(buf))
	{
		mprintf("    Loading configuration from %s ...\n",(const char*)buf);
		g_pDatabase->ParseInputFile(buf);
		return;
	}

//	sprintf(buf,"%s%ctravis.conf",g_sHomeDir,sep);
	buf.sprintf("%s%ctravis.conf",g_sHomeDir,sep);
	if (FileExist(buf))
	{
		mprintf("    Loading configuration from %s ...\n",(const char*)buf);
		g_pDatabase->ParseInputFile(buf);
		return;
	}

//	sprintf(buf,"%s%c.travis.conf",g_sHomeDir,sep);
	buf.sprintf("%s%c.travis.conf",g_sHomeDir,sep);
	if (FileExist(buf))
	{
		mprintf("    Loading configuration from %s ...\n",(const char*)buf);
		g_pDatabase->ParseInputFile(buf);
		return;
	}

	mprintf(WHITE,"    No configuration file found.\n");
	if (g_sHomeDir == NULL)
	{
		eprintf("\n    Could not detect user home directory, writing config file to current directory.\n\n");
//		sprintf(buf,".travis.conf");
		buf.sprintf(".travis.conf");
	} else
//		sprintf(buf,"%s%c.travis.conf",g_sHomeDir,sep);
		buf.sprintf("%s%c.travis.conf",g_sHomeDir,sep);
	mprintf("    Writing default configuration to %s ...\n",(const char*)buf);
	a = fopen(buf,"wt");
	if (a == NULL)
	{
		eprintf("    Cannot open %s for writing.\n",(const char*)buf);
		return;
	}
	fclose(a);

	g_pDatabase->WriteOutputFile(buf);
}


void InitDatabase()
{
//	char buf[256];
	CxString buf;
	char sep;
	FILE *a;

#ifdef TARGET_WINDOWS
	sep = '\\';
#else
	sep = '/';
#endif

	try { g_pDatabase = new CDatabase(); } catch(...) { g_pDatabase = NULL; }
	if (g_pDatabase == NULL) NewException((double)sizeof(CDatabase),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	CreateDatabaseDefaults();

	/******* Interface *********/
	Interface_DefaultConf();

	LoadSettings();
	
//	mprintf("\"%s\" vs \"%s\"\n",g_pDatabase->GetString("/GLOBAL/TRAVIS_VERSION"),SOURCE_VERSION);

	if (strcmp(g_pDatabase->GetString("/GLOBAL/TRAVIS_VERSION"),SOURCE_VERSION) != 0)
	{
		g_pDatabase->SetString("/GLOBAL/TRAVIS_VERSION",SOURCE_VERSION);
//		sprintf(buf,"%s%c.travis.conf",g_sHomeDir,sep);
		buf.sprintf("%s%c.travis.conf",g_sHomeDir,sep);
		mprintf(WHITE,"\n    TRAVIS version has changed.\n");
		if (!g_bWriteConf)
			mprintf("\n    You should write a new configuration file (command line argument \"-writeconf\").\n");
	}

	if (g_bWriteConf)
	{
		if (g_sHomeDir == NULL)
		{
			eprintf("\n    Could not detect user's home directory, writing config file to current directory.\n\n");
//			sprintf(buf,".travis.conf");
			buf.sprintf(".travis.conf");
		} else
//			sprintf(buf,"%s%c.travis.conf",g_sHomeDir,sep);
			buf.sprintf("%s%c.travis.conf",g_sHomeDir,sep);

		mprintf("    Writing new default configuration to %s ...\n",(const char*)buf);
		a = fopen(buf,"wt");
		if (a == NULL)
		{
			eprintf("    Cannot open %s for writing.\n",(const char*)buf);
			goto _writeend;
		}
		fclose(a);
		g_pDatabase->WriteOutputFile(buf);
	}
_writeend:

	if (g_bShowConf)
	{
		mprintf(WHITE,"\n  Output of Database Tree:\n\n");
		g_pDatabase->DumpTree();
	}

	if (g_pDatabase->GetBool("/GLOBAL/SHOWCREDITS"))
		g_bShowCredits = true;

	mprintf("\n");
}


void RECURSION_BuildCDF(CObservation *o, int channel, int om, CxDoubleArray **data, double *result)
{
	BXIN;
	int z/*, z2*/;

	if (channel == g_iCDFChannels)
	{
		o->m_pCDF->AddToBin(result);
		BXOUT;
		return;
	}
/*	if (o->m_pCDF->m_bChannelAll[channel])
	{
		for (z=0;z<o->m_iShowMolCount;z++)
			for (z2=0;z2<data[channel][z].GetSize();z2++)
			{
				result[channel] = data[channel][z][z2];
				RECURSION_BuildCDF(o,channel+1,om,data,result);
			}
	} else*/
	{
		for (z=0;z<data[channel][om].GetSize();z++)
		{
			result[channel] = data[channel][om][z];
			RECURSION_BuildCDF(o,channel+1,om,data,result);
		}
	}
	BXOUT;
}


CVirtualAtom* AddVirtualAtom(int mol)
{
	BTIN;
	int z2;
	CMolecule *m;
	CVirtualAtom *va;
	CSingleMolecule *sm;
	CxIntArray *la;

	m = (CMolecule*)g_oaMolecules[mol];

	try { va = new CVirtualAtom(); } catch(...) { va = NULL; }
	if (va == NULL) NewException((double)sizeof(CVirtualAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	g_oaVirtualAtoms.Add(va);
	va->m_iIndex = g_iGesVirtAtomCount;
	va->m_iMolecule = mol;
	va->m_iMolVirtAtom = (unsigned char)m->m_laVirtualAtoms.GetSize();
	sprintf(va->m_sLabel,"#%d",va->m_iMolVirtAtom+1);
	if (m->m_laVirtualAtoms.GetSize()==0)
	{
		m->m_baAtomIndex.Add(g_iVirtAtomType);
		m->m_waAtomCount.Add(1);
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			sm->m_baAtomIndex.Add(g_iVirtAtomType);

			try { la = new CxIntArray("AddVirtualAtom():sm->m_oaAtomOffset[]"); } catch(...) { la = NULL; }
			if (la == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			sm->m_oaAtomOffset.Add(la);
		}
	} else
	{
		m->m_waAtomCount[m->m_baAtomIndex.GetSize()-1]++;
//		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
//		{
//			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
//			sm->m_iAtomCount[sm->m_iElements-1]++;
//		}
	}
	for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
		((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->Add(g_iGesVirtAtomCount);
		g_iGesVirtAtomCount++;
	}
	m->m_iAtomGes++;
	m->m_laVirtualAtoms.Add((unsigned short)g_oaVirtualAtoms.GetSize()-1);
	BTOUT;
	return va;
}


void AddElement(const char *s, int ord, float mass, float radius, float vdw)
{
	BTIN;
	CElement *e;

	try { e = new CElement(); } catch(...) { e = NULL; }
	if (e == NULL) NewException((double)sizeof(CElement),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { e->m_sLabel = new char[strlen(s)+1]; } catch(...) { e->m_sLabel = NULL; }
	if (e->m_sLabel == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(e->m_sLabel,s);
	e->m_iOrd = ord;
	e->m_fMass = mass;
	e->m_fRadius = radius;
	e->m_fVdWRadius = vdw;
	e->m_fCoherentNCS = 0;
	g_oaElements.Add(e);
	BTOUT;
}


void AddElement(const char *s, int ord, float mass, float radius, float vdw, float ncs)
{
	BTIN;
	CElement *e;

	try { e = new CElement(); } catch(...) { e = NULL; }
	if (e == NULL) NewException((double)sizeof(CElement),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { e->m_sLabel = new char[strlen(s)+1]; } catch(...) { e->m_sLabel = NULL; }
	if (e->m_sLabel == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(e->m_sLabel,s);
	e->m_iOrd = ord;
	e->m_fMass = mass;
	e->m_fRadius = radius;
	e->m_fVdWRadius = vdw;
	e->m_fCoherentNCS = ncs;
	g_oaElements.Add(e);
	BTOUT;
}


void SetElementColor(const char *s, unsigned char r, unsigned char g, unsigned char b, unsigned char br, unsigned char bb, unsigned char bg)
{
	BTIN;
	CElement *e;
	e = FindElement(s,true);
	if (e == NULL)
	{
		eprintf("SetElementColor(): Element \"%s\" not found.\n");
		return;
	}
	e->m_iColorR = r;
	e->m_iColorG = g;
	e->m_iColorB = b;
	e->m_iColorBleachedR = br;
	e->m_iColorBleachedG = bg;
	e->m_iColorBleachedB = bb;
	BTOUT;
}


void RemoveAllElements()
{
	BTIN;
	int z;
	for (z=0;z<g_oaElements.GetSize();z++)
		delete (CElement*)g_oaElements[z];
	g_oaElements.RemoveAll();
	BTOUT;
}


void RemoveAllAtoms()
{
	BTIN;
	int z;
	for (z=0;z<g_oaAtoms.GetSize();z++)
		delete (CAtom*)g_oaAtoms[z];
	g_oaAtoms.RemoveAll();
	BTOUT;
}


void RemoveAllAnalyses()
{
	BTIN;
	int z;
	for (z=0;z<g_oaAnalyses.GetSize();z++)
		delete (CAnalysis*)g_oaAnalyses[z];
	g_oaAnalyses.RemoveAll();
	for (z=0;z<g_oaAnalysisGroups.GetSize();z++)
		delete (CAnalysisGroup*)g_oaAnalysisGroups[z];
	g_oaAnalysisGroups.RemoveAll();
	BTOUT;
}


void RemoveAllMolecules()
{
	BTIN;
	int z;

	for (z=0;z<g_oaMolecules.GetSize();z++)
		delete (CMolecule*)g_oaMolecules[z];
	g_oaMolecules.RemoveAll();
	for (z=0;z<g_oaSingleMolecules.GetSize();z++)
		delete (CSingleMolecule*)g_oaSingleMolecules[z];
	g_oaSingleMolecules.RemoveAll();
	BTOUT;
}


void RemoveAllObservations()
{
	BTIN;
	int z;

	for (z=0;z<g_oaObserv.GetSize();z++)
		delete (CObservation*)g_oaObserv[z];
	g_oaObserv.RemoveAll();
	BTOUT;
}


void GetTravisPath()
{
	BTIN;
	char *p, *q, *tmp, *env;

	if (FileExist(g_sExeName))
		return;
	if ((strchr(g_sExeName,'/')!=NULL) || (strchr(g_sExeName,'\\')!=NULL))
		return;

	tmp = getenv("PATH");
	if (tmp == NULL)
	{
		BTOUT;
		return;
	}

	try { env = new char[strlen(tmp)+1]; } catch(...) { env = NULL; }
	if (env == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(env,tmp);
	p = env;
	while (true)
	{
#ifdef TARGET_WINDOWS
		q = strchr(p,';');
#else
		q = strchr(p,':');
#endif
		if (q != NULL)
			*q = 0;

		try { tmp = new char[strlen(p)+strlen(g_sExeName)+2]; } catch(...) { tmp = NULL; }
		if (tmp == NULL) NewException((double)(strlen(p)+strlen(g_sExeName)+2)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(tmp,p);
#ifdef TARGET_WINDOWS
		strcat(tmp,"\\");
#else
		strcat(tmp,"/");
#endif
		strcat(tmp,g_sExeName);
		if (FileExist(tmp))
		{
			delete[] g_sExeName;

			try { g_sExeName = new char[strlen(tmp)+1]; } catch(...) { g_sExeName = NULL; }
			if (g_sExeName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(g_sExeName,tmp);
			delete[] env;
			BTOUT;
			return;
		}
		delete[] tmp;
		if (q != NULL)
			p = q+1;
				else break;
	}
	delete[] env;
	BTOUT;
}


void ReorderAtoms(int molecule)
{ // Written @ Ballmer Peak
	int z, z2, z3;
	CMolecule *mol;
	CSingleMolecule *sm;
	CxIntArray *twa;

	mol = (CMolecule*)g_oaMolecules[molecule];

	twa = NULL;

	for (z=0;z<mol->m_laSingleMolIndex.GetSize();z++)
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[mol->m_laSingleMolIndex[z]];
//		sm->m_oaTempAtomOffset.SetSize(mol->m_baAtomIndex.GetSize());

		for (z2=0;z2<mol->m_baAtomIndex.GetSize();z2++)
		{
/*			if (sm->m_oaTempAtomOffset[z2] != NULL)
				delete sm->m_oaTempAtomOffset[z2];
			sm->m_oaTempAtomOffset[z2] = new CxIntArray();*/
			if (twa != NULL)
				delete twa;

			try { twa = new CxIntArray("ReorderAtoms():twa"); } catch(...) { twa = NULL; }
			if (twa == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			twa->SetSize(((CxIntArray*)mol->m_oaNewNumbers[z2])->GetSize());
			for (z3=0;z3<((CxIntArray*)mol->m_oaNewNumbers[z2])->GetSize();z3++)
			{
//				ti = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3);
//				(*((CxIntArray*)sm->m_oaAtomOffset[z2]))[z3] = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z3));
				(*twa)[z3] = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z3));
//				mprintf("  (%s%d wird zu %s%d; %d -> %d)\n",((CAtom*)g_oaAtoms[mol->m_baAtomIndex[z2]])->m_sName,z3+1,((CAtom*)g_oaAtoms[mol->m_baAtomIndex[z2]])->m_sName,((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z3)+1,((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3),((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z3)));
			}
			for (z3=0;z3<((CxIntArray*)mol->m_oaNewNumbers[z2])->GetSize();z3++)
				(*((CxIntArray*)sm->m_oaAtomOffset[z2]))[z3] = (*twa)[z3];
		}
	}
	delete twa;
}


void ReorderLikeInput()
{ // Written @ Ballmer Peak
	int z, z2, z3, z4, z5, ti, tl;
	CMolecule *mol;
	CSingleMolecule *sm;

	mprintf("\n");
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		mol = (CMolecule*)g_oaMolecules[z];
		mprintf("Reordering atoms in %s...",mol->m_sName);
		mol->m_oaNewNumbers.SetSize(mol->m_baAtomIndex.GetSize());
		sm = (CSingleMolecule*)g_oaSingleMolecules[mol->m_laSingleMolIndex[0]];
		for (z2=0;z2<mol->m_baAtomIndex.GetSize();z2++)
		{
			if (mol->m_oaNewNumbers[z2] != NULL)
				delete mol->m_oaNewNumbers[z2];

			try { mol->m_oaNewNumbers[z2] = new CxIntArray("ReorderLikeInput():mol->m_oaNewNumbers[z2]"); } catch(...) { mol->m_oaNewNumbers[z2] = NULL; }
			if (mol->m_oaNewNumbers[z2] == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			// F***ing Stack sort is incompatible with Ballmer Peak ^^
			for (z3=0;z3<mol->m_waAtomCount[z2];z3++)
			{
				// Sortiere nach kleinstem Offset
				tl = 9999999;
				ti = -1;
				for (z4=0;z4<mol->m_waAtomCount[z2];z4++)
				{
					for (z5=0;z5<z3;z5++)
						if (((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z5) == z4)
							goto _nextreord; // Der ist schon einsortiert. Keyswapping hier offensichtlich nicht moeglich
					if (((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z4) < tl)
					{
						tl = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z4);
						ti = z4;
					}
_nextreord:;
				}
				((CxIntArray*)mol->m_oaNewNumbers[z2])->Add(ti);
			}
		}
		ReorderAtoms(z);
		mprintf("Done.\n");
	}
	mprintf("\n");
}


unsigned long GraceColor(int z, double bleach)
{
	int r, g, b;

	z = z % 14;
	switch(z)
	{
		case 0: r = 0; g = 0; b = 0; break;
		case 1: r = 0; g = 0; b = 255; break;
		case 2: r = 255; g = 0; b = 0; break;
		case 3: r = 0; g = 255; b = 0; break;
		case 4: r = 255; g = 0; b = 255; break;
		case 5: r = 0; g = 255; b = 255; break;
		case 6: r = 255; g = 255; b = 0; break;
		case 7: r = 128; g = 128; b = 128; break;
		case 8: r = 0; g = 0; b = 128; break;
		case 9: r = 128; g = 0; b = 0; break;
		case 10: r = 0; g = 128; b = 0; break;
		case 11: r = 128; g = 0; b = 128; break;
		case 12: r = 0; g = 128; b = 128; break;
		case 13: r = 128; g = 255; b = 0; break;
		default: r = 0; g = 0; b = 0;
	}
	r = (int)((1.0-bleach)*r + bleach*255.0);
	g = (int)((1.0-bleach)*g + bleach*255.0);
	b = (int)((1.0-bleach)*b + bleach*255.0);
	if (r > 255)
		r = 255;
	if (g > 255)
		g = 255;
	if (b > 255)
		b = 255;
	return r + g*0x100 + b*0x10000;
}


unsigned long CalcFFTSize(unsigned long i, bool silent)
{
	unsigned long t, r, r2, r3, r5, t2, t3, t5;
	int i2, i3, i5, m2, m3, m5;
	bool b;

//	mprintf("*** CalcFFTSize %d ***\n",i);

	m2 = (int)ceil(log((double)i) / log(2.0));
	m3 = (int)ceil(log((double)i) / log(3.0));
	m5 = (int)ceil(log((double)i) / log(5.0));

	r = 2147483647;
	r2 = m2;
	r3 = m3;
	r5 = m5;
//	mprintf("m2=%d, m3=%d, m5=%d, r=%d.\n",m2,m3,m5,r);
	t2 = 1;
	for (i2=0;i2<=m2;i2++)
	{
		t3 = 1;
		for (i3=0;i3<=m3;i3++)
		{
			t5 = 1;
			for (i5=0;i5<=m5;i5++)
			{
				t = t2 * t3 * t5;
	//			mprintf("%d * %d * %d = %d\n",t2,t3,t5,t);
				if ((t >= i) && (t < r))
				{
					r = t;
					r2 = i2;
					r3 = i3;
					r5 = i5;
//					mprintf("2^%d + 3^%d +5^%d = %d.\n",r2,r3,r5,r);
				}
				t5 *= 5;
			}
			t3 *= 3;
		}
		t2 *= 2;
	}
	if (!silent)
	{
		mprintf("\n    CalcFFTSize(): %d = ",i);
		if (r2 != 0)
		{
			mprintf("2^%d",r2);
			b = true;
		}
		if (r3 != 0)
		{
			if (b)
				mprintf(" *");
			mprintf(" 3^%d",r3);
			b = true;
		}
		if (r5 != 0)
		{
			if (b)
				mprintf(" *");
			mprintf(" 5^%d",r5);
		}
		if ((i - r) == 0)
			mprintf(". All prime factors within 2, 3, 5. Fine.\n\n");
				else mprintf(" - %d. Prime factors other than 2, 3, 5 not allowed.\n",r-i);
	}
//	mprintf("*** CalcFFTSize done ***\n");
	return r;
}


void BuildAtomIndices()
{
	int z, z2, z3, z4, z5, ti, ti2;
	CMolecule *m;
	CSingleMolecule *sm;

	g_waAtomMolIndex.SetSize(g_iGesVirtAtomCount);
	g_waAtomMolUID.SetSize(g_iGesVirtAtomCount);
	g_laAtomSMLocalIndex.SetSize(g_iGesVirtAtomCount);
	g_waAtomMolNumber.SetSize(g_iGesVirtAtomCount);
	g_waAtomElement.SetSize(g_iGesVirtAtomCount);
	g_waAtomRealElement.SetSize(g_iGesVirtAtomCount);
	g_faAtomCode.SetSize(g_iGesAtomCount);
	g_faVdWRadius.SetSize(g_iGesAtomCount);

	for (z=0;z<g_iGesVirtAtomCount;z++)
	{
		g_waAtomMolIndex[z] = 60000;
		g_waAtomRealElement[z] = 60000;
	}

	ti2 = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		if (m->m_bPseudo)
			continue;
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			for (z3=0;z3<sm->m_baAtomIndex.GetSize();z3++)
			{
				for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
				{
					ti = ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4);
					if (g_waAtomMolIndex[ti] != 60000)
						eprintf("BuildAtomIndices(): Atom in more than one molecule! z=%d, z2=%d, z3=%d, z4=%d.\n",z,z2,z3,z4);
					g_waAtomMolIndex[ti] = z;
					g_waAtomMolUID[ti] = ti2;
					g_laAtomSMLocalIndex[ti] = z2;
					g_waAtomElement[ti] = z3;
					g_waAtomMolNumber[ti] = z4;
					g_waAtomRealElement[ti] = sm->m_baAtomIndex[z3];

					if (sm->m_baAtomIndex[z3] == g_iVirtAtomType)
						continue;

					g_faVdWRadius[ti] = ((CAtom*)g_oaAtoms[sm->m_baAtomIndex[z3]])->m_pElement->m_fVdWRadius;

					for (z5=0;z5<sm->m_oaMolAtoms.GetSize();z5++)
					{
						if (((CMolAtom*)sm->m_oaMolAtoms[z5])->m_iOffset == ti)
						{
							g_faAtomCode[ti] = ((CMolAtom*)sm->m_oaMolAtoms[z5])->m_fAtomCode;
							goto _found;
						}
					}
					eprintf("BuildAtomIndices(): Atom Code not found (z=%d, z2=%d, z3=%d, z4=%d).\n",z,z2,z3,z4);
_found:;
				}
			}
			ti2++;
		}
	}
	g_baAtomPassedCondition.SetSize(g_iGesVirtAtomCount);
	for (z=0;z<g_iGesVirtAtomCount;z++)
		g_baAtomPassedCondition[z] = 0;
}


bool DetermineTrajFormat()
{
	char *p;

	p = strrchr(g_sInputTraj,'.');
	if (p == NULL)
	{
		if (strlen(g_sInputTraj) >= 7)
		{
			if (mystricmp(&g_sInputTraj[strlen(g_sInputTraj)-7],"HISTORY") == 0)
			{
				g_iTrajFormat = 4;
				return true;
			}
		}
		goto _unk;
	}
	p++;
	if (mystricmp(p,"xyz")==0)
	{
		g_iTrajFormat = 0;
		return true;
	}
	if (mystricmp(p,"pdb")==0)
	{
		g_iTrajFormat = 1;
		return true;
	}
	if (mystricmp(p,"mol2")==0)
	{
		g_iTrajFormat = 2;
		return true;
	}
	if (mystricmp(p,"lmp")==0)
	{
		g_iTrajFormat = 3;
		return true;
	}
	if (mystricmp(p,"cube")==0)
	{
		g_iTrajFormat = 5;
		return true;
	}
	if ((mystricmp(p,"prmtop")==0) || (mystricmp(p,"mdcrd")==0))
	{
		g_iTrajFormat = 6;
		*p = 0; // Delete file extension
		return true;
	}
_unk:
	eprintf("Could not determine trajectory file format.\n\n",p);
	mprintf(WHITE,"The following formats are supported:\n");
	mprintf("  - XYZ trajectories (extension .xyz)\n");
	mprintf("  - PDB trajectories (extension .pdb)\n");
	mprintf("  - LAMMPS trajectories \"dump custom element xu yu zu\" (extension .lmp)\n");
	mprintf("  - DLPOLY trajectories (file name \"HISTORY\", no extension)\n");
	mprintf("  - Amber trajectories (extension .prmtop / .mdcrd)\n");
	mprintf("  - Cube file trajectories (extension .cube)\n");
	mprintf("\n");
	return false;
}


void PrintSMode()
{
	const char *stext[16];
	char buf[256];
	int z;

	stext[0]  = "ICAgU1NTU1NTU1NTU1NTU1NTICAgICAgICAgICAgICAgICAgaGhoaGhoICAgICAgICAgICAgICAgaWlpICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgISEh";
	stext[1]  = "IFNTOjo6Ojo6Ojo6Ojo6Ojo6UyAgICAgICAgICAgICAgICAgaDo6OjpoICAgICAgICAgICAgICBpOjo6aSAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAhITohIQ==";
	stext[2]  = "Uzo6Ojo6U1NTU1NTOjo6Ojo6UyAgICAgICAgICAgICAgICAgaDo6OjpoICAgICAgICAgICAgICAgaWlpICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAhOjo6IQ==";
	stext[3]  = "Uzo6Ojo6UyAgICAgU1NTU1NTUyAgICAgICAgICAgICAgICAgIGg6OjpoICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAhOjo6IQ==";
	stext[4]  = "Uzo6Ojo6UyAgICAgICAgICAgICAgICAgY2NjY2NjY2NjY2NjIGg6OjpoIGhoaGhoICAgICAgIGlpaWlpaSAgICAgIHNzc3Nzc3NzICAgICAgICBzc3Nzc3NzcyAgICAhOjo6IQ==";
	stext[5]  = "Uzo6Ojo6UyAgICAgICAgICAgICAgIGNjOjo6Ojo6Ojo6OjpjIGg6OjpoaDo6Ojo6aGggICAgIGk6Ojo6aSAgICBzczo6Ojo6Ojo6cyAgICAgc3M6Ojo6Ojo6OnMgICAhOjo6IQ==";
	stext[6]  = "IFM6Ojo6U1NTUyAgICAgICAgICAgYzo6Ojo6Ojo6Ojo6OjpjIGg6Ojo6Ojo6Ojo6OjpoaCAgICBpOjo6aSAgc3M6Ojo6Ojo6Ojo6OnMgIHNzOjo6Ojo6Ojo6OjpzICAhOjo6IQ==";
	stext[7]  = "ICBTUzo6Ojo6OlNTU1NTICAgICBjOjo6OjpjY2NjY2M6OjpjIGg6Ojo6OjpoaGg6Ojo6OmggICBpOjo6aSAgczo6Ojo6c3Nzczo6OjpzIHM6Ojo6OnNzc3M6Ojo6cyAhOjo6IQ==";
	stext[8]  = "ICAgIFNTUzo6Ojo6Ojo6U1MgICBjOjo6OmMgICAgIGNjY2NjIGg6Ojo6OmggICBoOjo6OjpoICBpOjo6aSAgIHM6Ojo6cyAgIHNzc3MgICBzOjo6OnMgICBzc3NzICAhOjo6IQ==";
	stext[9]  = "ICAgICAgIFNTU1NTUzo6OjpTICBjOjo6YyAgICAgICAgICAgIGg6Ojo6aCAgICAgaDo6OjpoICBpOjo6aSAgICAgczo6Ojo6cyAgICAgICAgIHM6Ojo6OnMgICAgICAhOjo6IQ==";
	stext[10] = "ICAgICAgICAgICAgUzo6Ojo6UyBjOjo6YyAgICAgICAgICAgIGg6Ojo6aCAgICAgaDo6OjpoICBpOjo6aSAgICAgICBzOjo6OjpzICAgICAgICAgczo6Ojo6cyAgICAhITohIQ==";
	stext[11] = "ICAgICAgICAgICAgUzo6Ojo6UyBjOjo6OmMgICAgIGNjY2NjIGg6Ojo6aCAgICAgaDo6OjpoICBpOjo6aSAgc3NzcyAgICBzOjo6OnMgIHNzc3MgICAgczo6OjpzICAgISEh";
	stext[12] = "U1NTU1NTUyAgICAgUzo6Ojo6UyBjOjo6OjpjY2NjY2M6OjpjIGg6Ojo6aCAgICAgaDo6OjpoIGk6Ojo6Omkgczo6Ojpzc3NzOjo6OjpzIHM6Ojo6c3Nzczo6Ojo6cw==";
	stext[13] = "Uzo6Ojo6OlNTU1NTUzo6Ojo6UyAgYzo6Ojo6Ojo6Ojo6OjpjIGg6Ojo6aCAgICAgaDo6OjpoIGk6Ojo6Omkgczo6Ojo6Ojo6Ojo6OnMgIHM6Ojo6Ojo6Ojo6OjpzICAgISEh";
	stext[14] = "Uzo6Ojo6Ojo6Ojo6Ojo6OlNTICAgIGNjOjo6Ojo6Ojo6OjpjIGg6Ojo6aCAgICAgaDo6OjpoIGk6Ojo6OmkgIHM6Ojo6Ojo6OjpzcyAgICBzOjo6Ojo6Ojo6c3MgICAhITohIQ==";
	stext[15] = "IFNTU1NTU1NTU1NTU1NTUyAgICAgICAgY2NjY2NjY2NjY2NjIGhoaGhoaCAgICAgaGhoaGhoIGlpaWlpaWkgICBzc3Nzc3Nzc3MgICAgICAgc3Nzc3Nzc3NzICAgICAgISEh";

	for (z=0;z<16;z++)
	{
		UnBase64((unsigned char*)buf,(const unsigned char*)stext[z],strlen(stext[z]));
		mprintf(YELLOW,"%s\n",buf);
	}
	mprintf("\n");
}


void WriteCredits()
{
	mprintf("\n");
	if (g_bShowCredits)
	{
		WriteCredits_Long();
	} else
	{
		mprintf(YELLOW,"    Note:");
		mprintf(" To show a list of all persons who contributed to TRAVIS,\n");
		mprintf("          please add \"");
		mprintf(WHITE,"-credits");
		mprintf("\" to your command line arguments, or set the\n");
		mprintf("          variable \"SHOWCREDITS\" to \"TRUE\" in your travis.conf file.\n");
	}
	mprintf("\n");
	mprintf("    Source code from other projects used in TRAVIS:\n");
	mprintf("      - lmfit     from "); mprintf(WHITE,"Joachim Wuttke\n");
	mprintf("      - kiss_fft  from "); mprintf(WHITE,"Mark Borgerding\n");
	mprintf("      - voro++    from "); mprintf(WHITE,"Chris Rycroft\n");
	mprintf("\n");
	mprintf(WHITE,"    http://www.travis-analyzer.de\n\n");
	mprintf(YELLOW,"    Please cite:\n\n");
	mprintf(RED,"  * ");
	mprintf("\"TRAVIS - A Free Analyzer and Visualizer for Monte Carlo and Molecular Dynamics Trajectories\",\n    M. Brehm, B. Kirchner; J. Chem. Inf. Model. 2011, 51 (8), pp 2007-2023.\n\n");
	if (g_bIRSpec || g_bPowerSpec || g_bRaman || g_bPower || g_bIR)
	{
		mprintf(RED,"  * ");
		mprintf("\"Computing vibrational spectra from ab initio molecular dynamics\",\n    M. Thomas, M. Brehm, R. Fligg, P. Voehringer, B. Kirchner; Phys. Chem. Chem. Phys. 2013, 15, pp 6608-6622.\n\n");
	}
	if (g_bNormalCoordinate || g_bEckartTransform) {
		mprintf(RED,"  * ");
		mprintf("\"Simulating the vibrational spectra of ionic liquid systems:\n    1-Ethyl-3-methylimidazolium acetate and its mixtures\",\n    M. Thomas, M. Brehm, O. Holloczki, Z. Kelemen, L. Nyulaszi, T. Pasinszki, B. Kirchner;\n    J. Chem. Phys. 2014, 141, 024510.\n\n");
	}
	if (g_bVoro || g_bDomA) {
		mprintf(RED,"  * ");
		mprintf("\"Domain Analysis in Nanostructured Liquids: A Post-Molecular Dynamics Study at the Example of Ionic Liquids\",\n    M. Brehm, H. Weber, M. Thomas, O. Holloczki, B. Kirchner; ChemPhysChem 2015, 16, pp 3271-3277.\n\n");
	}
	if (g_bSFac) {
		mprintf(RED,"  * ");
		mprintf("\"Triphilic Ionic-Liquid Mixtures: Fluorinated and Non-fluorinated Aprotic Ionic-Liquid Mixtures\",\n    O. Holloczki, M. Macchiagodena, H. Weber, M. Thomas, M. Brehm, A. Stark, O. Russina, A. Triolo, B. Kirchner; ChemPhysChem 2015, 16, pp 3325-3333.\n\n");
	}
	if (g_bTegri) {
		mprintf(RED,"  * ");
		mprintf("\"Voronoi dipole moments for the simulation of bulk phase vibrational spectra\",\n    M. Thomas, M. Brehm, B. Kirchner; Phys. Chem. Chem. Phys. 2015, 17, pp 3207-3213.\n\n");
	}
	if (g_bVCD || g_bCubeTimeDev) {
		mprintf(RED,"  * ");
		mprintf("\"Classical Magnetic Dipole Moments for the Simulation of Vibrational Circular Dichroism by Ab Initio Molecular Dynamics\",\n    M. Thomas, B. Kirchner; J. Phys. Chem. Lett. 2016, 7, pp 509-513\n\n");
	}
}


void WriteCredits_Long()
{
	mprintf(YELLOW,"  >>> Credits <<<\n");
	mprintf(YELLOW,"\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Martin Brehm          "); mprintf("Main Developer (2009 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Martin Thomas         "); mprintf("Main Developer (2012 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Barbara Kirchner      "); mprintf("Group Leader and idea\n\n");

	mprintf(YELLOW,"    > "); mprintf(WHITE,"Marc Bruessel         "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Philipp di Dio        "); mprintf("Atom Parameters and Math Support\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Michael von Domaros   "); mprintf("New Ideas and Linux Repository updating\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Dzmitry Firaha        "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Dorothea Golze        "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Oldamur Holloczki     "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Daniela Kerle         "); mprintf("Name Finding\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Miriam Kohagen        "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Marina Macchiagodena  "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Fred Malberg          "); mprintf("Leading Bug Finder :-)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Matthias Schoeppke    "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Daniele Tedesco       "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Jens Thar             "); mprintf("Fruitful Discussion and Scientific Input\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Henry Weber           "); mprintf("Testing, Debugging and creative Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Stefan Zahn           "); mprintf("Testing and Debugging\n");

/*	mprintf(YELLOW,"  >>> Credits <<<\n");
	mprintf(YELLOW,"\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Martin Brehm "); mprintf("("); mprintf(RED,"*"); mprintf(")    "); mprintf("Main Developer (2009 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Martin Thomas "); mprintf("("); mprintf(RED,"*"); mprintf(")   "); mprintf("Main Developer (2012 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Barbara Kirchner    "); mprintf("Group Leader\n\n");

	mprintf(YELLOW,"    > "); mprintf(WHITE,"Marc Bruessel       "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Philipp di Dio "); mprintf("("); mprintf(RED,"*"); mprintf(")  "); mprintf("Atom Parameters and Math Support\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Dorothea Golze      "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Daniela Kerle       "); mprintf("Name Finding\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Miriam Kohagen      "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Fred Malberg        "); mprintf("Leading Bug Finder :-)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Matthias Schoeppke  "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Jens Thar "); mprintf("("); mprintf(RED,"*"); mprintf(")       Fruitful Discussion and Scientific Input\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Henry Weber "); mprintf("("); mprintf(RED,"*"); mprintf(")     ");mprintf("Testing, Debugging and creative Ideas\n");
	mprintf("\n");
	mprintf("    ("); mprintf(RED,"*"); mprintf(")"); mprintf(" Office 42 - \"answer to all questions\".\n");*/
}


CAutoCorrelation::CAutoCorrelation()
{
	m_iInput = 0;
	m_iDepth = 0;
	m_pFFT = NULL;
	m_pFFT2 = NULL;
}


CAutoCorrelation::~CAutoCorrelation()
{
}


void CAutoCorrelation::Init(int input, int depth, bool fft)
{
	m_iInput = input;
	m_iDepth = depth;
	if (fft)
	{
		m_bFFT = true;
		m_iFFTSize = CalcFFTSize(input,true);

		try { m_pFFT = new CFFT(); } catch(...) { m_pFFT = NULL; }
		if (m_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pFFT->PrepareFFT_C2C(2*m_iFFTSize);

		try { m_pFFT2 = new CFFT(); } catch(...) { m_pFFT2 = NULL; }
		if (m_pFFT2 == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pFFT2->PrepareInverseFFT_C2C(2*m_iFFTSize);
	} else
	{
		m_bFFT = false;
	}
}


void CAutoCorrelation::AutoCorrelate(CxFloatArray *inp, CxFloatArray *outp)
{
	int z, z2;
	double tf;

	outp->SetSize(m_iDepth);

	if (m_bFFT)
	{
		for (z=0;z<m_iInput;z++)
		{
			m_pFFT->m_pInput[z*2] = (*inp)[z];
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++)
		{
			m_pFFT->m_pInput[z*2] = 0;
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		m_pFFT->DoFFT();

		for (z=0;z<m_iFFTSize*2;z++)
		{
			m_pFFT2->m_pInput[z*2] = (float)((m_pFFT->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2] + m_pFFT->m_pOutput[z*2+1]*m_pFFT->m_pOutput[z*2+1]));
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		m_pFFT2->DoFFT();

		for (z=0;z<m_iDepth;z++)
			(*outp)[z] = (float)((double)m_pFFT2->m_pOutput[2*z] / m_iFFTSize / 2.0f / ((double)m_iInput - z));
	} else
	{
		for (z=0;z<m_iDepth;z++) // Tau
		{
			tf = 0;
			for (z2=0;z2<m_iInput-z;z2++)
				tf += (double)(*inp)[z2] * (*inp)[z2+z];
			(*outp)[z] = (float)(tf / (double)(m_iInput-z));
		}
	}
}


/*
void FormatTime(unsigned long eta, char *buf)
{
	char tbuf[256], tbuf2[256];

	if ((eta/60) > 0)
		sprintf(tbuf,"%02lus",eta%60);
			else sprintf(tbuf,"%2lus",eta%60);
	eta /= 60;
	if (eta > 0)
	{
		strcpy(tbuf2,tbuf);
		if ((eta/60) > 0)
			sprintf(tbuf,"%02lum",eta%60);
				else sprintf(tbuf,"%2lum",eta%60);
		strcat(tbuf,tbuf2);
		eta /= 60;
		if (eta > 0)
		{
			strcpy(tbuf2,tbuf);
			if ((eta/60) > 0)
				sprintf(tbuf,"%02luh",eta%24);
					else sprintf(tbuf,"%2luh",eta%24);
			strcat(tbuf,tbuf2);
			eta /= 24;
		}
		if (eta > 0)
		{
			strcpy(tbuf2,tbuf);
			if ((eta/60) > 0)
				sprintf(tbuf,"%02lud",eta);
					else sprintf(tbuf,"%2lud",eta);
			strcat(tbuf,tbuf2);
		}
	}
	strcpy(buf,tbuf);
}
*/


void FormatTime(unsigned long eta, CxString *buf)
{
//	char tbuf[256], tbuf2[256];
	CxString tbuf, tbuf2;

	if (buf == NULL)
		return;

	if ((eta/60) > 0)
//		sprintf(tbuf,"%02lus",eta%60);
		tbuf.sprintf("%02lus",eta%60);
	else
//		sprintf(tbuf,"%2lus",eta%60);
		tbuf.sprintf("%2lus",eta%60);

	eta /= 60;
	if (eta > 0)
	{
//		strcpy(tbuf2,tbuf);
		tbuf2.strcpy(tbuf);

		if ((eta/60) > 0)
//			sprintf(tbuf,"%02lum",eta%60);
			tbuf.sprintf("%02lum",eta%60);
		else
//			sprintf(tbuf,"%2lum",eta%60);
			tbuf.sprintf("%2lum",eta%60);

//		strcat(tbuf,tbuf2);
		tbuf.strcat(tbuf2);

		eta /= 60;
		if (eta > 0)
		{
//			strcpy(tbuf2,tbuf);
			tbuf2.strcpy(tbuf);

			if ((eta/60) > 0)
//				sprintf(tbuf,"%02luh",eta%24);
				tbuf.sprintf("%02luh",eta%24);
			else
//				sprintf(tbuf,"%2luh",eta%24);
				tbuf.sprintf("%2luh",eta%24);

//			strcat(tbuf,tbuf2);
			tbuf.strcat(tbuf2);

			eta /= 24;
		}
		if (eta > 0)
		{
//			strcpy(tbuf2,tbuf);
			tbuf2.strcpy(tbuf);

			if ((eta/60) > 0)
//				sprintf(tbuf,"%02lud",eta);
				tbuf.sprintf("%02lud",eta);
			else
//				sprintf(tbuf,"%2lud",eta);
				tbuf.sprintf("%2lud",eta);

//			strcat(tbuf,tbuf2);
			tbuf.strcat(tbuf2);
		}
	}

//	strcpy(buf,tbuf);
	buf->strcpy(tbuf);
}


void RenderStructFormulas(int tries)
{
	int z, z2, z3, ti;
	bool nohyd;
	double tf;
	FILE *a;
//	char buf[256];
	CxString buf;
	CMolecule *m;
	CSingleMolecule *sm;
	CMolAtom *ma1, *ma2;

	mprintf("\n");
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
		nohyd = false;

_again:
		if (nohyd)
		{
//			sprintf(buf,"mol%d_%s_no_H.dot",z+1,m->m_sName);
			buf.sprintf("mol%d_%s_no_H.dot",z+1,m->m_sName);
			mprintf(WHITE,"    * Molecule %d: %s (without hydrogen atoms)\n",z+1,m->m_sName);
		} else
		{
//			sprintf(buf,"mol%d_%s.dot",z+1,m->m_sName);
			buf.sprintf("mol%d_%s.dot",z+1,m->m_sName);
			mprintf(WHITE,"    * Molecule %d: %s\n",z+1,m->m_sName);
		}
		mprintf("      Writing dot file %s...\n",(const char*)buf);
		a = OpenFileWrite(buf,true);
		mfprintf(a,"graph molecule {\n");
		mfprintf(a,"  graph [pack=true,splines=true,overlap=false];\n");
		mfprintf(a,"  node [shape=none,fontsize=16,fontname=\"Arial\",margin=0,fixedsize=true,height=0.28];\n");
		mfprintf(a,"  edge [style=bold,len=0.70];\n");

		for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++)
		{
			ma1 = (CMolAtom*)sm->m_oaMolAtoms[z2];
			if (nohyd && (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,"H")==0))
				continue;
			mfprintf(a,"  %s%d [label=\"%s%d\",width=%.2f];\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,ma1->m_iNumber+1,((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,ma1->m_iNumber+1,(strlen(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName)+(((ma1->m_iNumber+1)>9)?2:1))*0.17f);
		}

		for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++)
		{
			ma1 = (CMolAtom*)sm->m_oaMolAtoms[z2];
			if (nohyd && (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,"H")==0))
				continue;
			for (z3=0;z3<ma1->m_oaBonds.GetSize();z3++)
			{
				ma2 = (CMolAtom*)ma1->m_oaBonds[z3];
				if (ma2->m_iMolAtomNumber < ma1->m_iMolAtomNumber)
					continue;
				if (nohyd && (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma2->m_iType]])->m_sName,"H")==0))
					continue;
				// Calculate edge weight: C-C is most important
				if ((mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,"C")==0) && (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma2->m_iType]])->m_sName,"C")==0))
				{
					ti = 10000;
					tf = 2.0;
				} else 
				{
					ti = ((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_pElement->m_iOrd * ((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma2->m_iType]])->m_pElement->m_iOrd;
					tf = 1.0;
				}
				if (nohyd)
					ti = 40;
				tf = 0.5 + (((ti>80)?80:ti)/80.0*2.0);
				mfprintf(a,"  %s%d -- %s%d [weight=%d, penwidth=%f];\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,ma1->m_iNumber+1,((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma2->m_iType]])->m_sName,ma2->m_iNumber+1,ti,tf);
			}
		}
		mfprintf(a,"}\n");
		fclose(a);
//		if (nohyd)
//			sprintf(buf,"dot mol%d_%s_no_H.dot -Tpng -omol%d_%s_no_H.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000",z+1,m->m_sName,z+1,m->m_sName);
//				else sprintf(buf,"dot mol%d_%s.dot -Tpng -omol%d_%s.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000",z+1,m->m_sName,z+1,m->m_sName);
		if (nohyd)
//			sprintf(buf,"mol%d_%s_no_H",z+1,m->m_sName);
			buf.sprintf("mol%d_%s_no_H",z+1,m->m_sName);
		else
//			sprintf(buf,"mol%d_%s",z+1,m->m_sName);
			buf.sprintf("mol%d_%s",z+1,m->m_sName);

		RenderFormula(buf,tries);

		if ((!nohyd) && (m->m_iAtomGes > 30))
		{
			nohyd = true;
			goto _again;
		}
	}
	mprintf("    If the command above worked, you can now view the PNG files that have been created.\n\n");
}


void RenderFormula(const char *s, int tries)
{
	int z, zm;
	char buf[256], buf2[256], *p, *q;
	double tf, mi, ma, av;
	FILE *a;
	
//	mprintf("%d Tries:\n",tries);
	zm = -1;
	av = 0;
	ma = 0;
	mi = 1e30;
	mprintf("      Command: dot %s.dot -Tpng -o%s.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=%d\n",s,s,rand());
	mprintf("      Optimizing (%d tries): ",tries);
	mprintf(WHITE,"[");
	for (z=0;z<tries;z++)
	{
		mprintf(WHITE,"#");
		sprintf(buf,"dot %s.dot -Tpng -o%s.%d.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=%d > dot%d.log 2>&1",s,s,z,rand(),z);

//		mprintf("    (%s)\n",buf);

		system(buf);

		sprintf(buf,"dot%d.log",z);
		a = fopen(buf,"rt");
		if (a == NULL)
		{
			eprintf("\nRenderFormula(): Error opening %s for reading.\n",buf);
			eprintf("It seems that GraphViz is not working (try command \"dot\").\n\n");
			return;
		}
		while (!feof(a))
		{
			fgets(buf,256,a);
			if (strstr(buf,"final") != NULL)
			{
				p = buf;
	//			mprintf("    buf=\"%s\".\n",buf);
				while (((*p < '0') || (*p > '9')) && (*p != 0))
					p++;
	//			mprintf("    p=\"%s\".\n",p);
				if (*p == 0)
					continue;
				q = p;
//				mprintf("    *p = %d.\n",*p);
				while (((*p >= '0') && (*p <= '9')) || (*p == '.'))
				{
	//				mprintf("    *p = %d --> p++;\n",*p);
					p++;
				}
	//			mprintf("    p2=\"%s\".\n",p);
				*p = 0;
		//		mprintf("    q=\"%s\".\n",q);
				tf = atof(q);
//				mprintf("    tf = %g.\n",tf);
				if (tf > ma)
					ma = tf;
				if (tf < mi)
				{
					mi = tf;
					zm = z;
				}
				av += tf;
				goto _done;
			}
		}
		eprintf("\nError with GraphViz output. See dot%d.log.\n",z);
		eprintf("It seems like GraphViz is not working (try command \"dot\").\n\n");
		return;
_done:
		fclose(a);
		sprintf(buf,"dot%d.log",z);
		remove(buf);
	}
	mprintf(WHITE,"]\n");
	av /= tries;
	mprintf("      Quality statistics: Min %g, Max %g, Avg %g.\n",mi,ma,av);
	mprintf("      Using best result to produce output file.\n");
	for (z=0;z<tries;z++)
	{
		if (z == zm)
			continue;
		sprintf(buf,"%s.%d.png",s,z);
		if (remove(buf) != 0)
			eprintf("      Error removing \"%s\".\n",buf);
	}
	sprintf(buf,"%s.%d.png",s,zm);
	sprintf(buf2,"%s.png",s);
	if (rename(buf,buf2) != 0)
		eprintf("      Error renaming \"%s\" to \"%s\".\n",buf,buf2);
	mprintf("\n");
}

void parseCoreCharges() {
	static bool coreChargesDefined = false;
	if (coreChargesDefined)
		return;
	
	mprintf(WHITE, "\n  > Core Charge Definition >\n\n");
		mprintf("    Please mind pseudopotentials while entering core charges.\n\n");
		int i;
		for (i = 0; i < g_oaAtoms.GetSize(); i++) {
			if (i == g_iWannierAtomType)
				continue;
			if (i == g_iVirtAtomType)
				continue;
			float def = 0.0f;
			if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "B") == 0) def = 3.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "C") == 0) def = 4.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "N") == 0) def = 5.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "O") == 0) def = 6.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "F") == 0) def = 7.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Si") == 0) def = 4.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "P") == 0) def = 5.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "S") == 0) def = 6.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Cl") == 0) def = 7.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Br") == 0) def = 7.0f;
			else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "I") == 0) def = 7.0f;
			else def = (float)((CAtom *)g_oaAtoms[i])->m_pElement->m_iOrd;
			((CAtom *)g_oaAtoms[i])->m_fCharge = AskFloat("    Enter core charge for atom type %s: [%.1f] ", def, ((CAtom*)g_oaAtoms[i])->m_sName, def);
		}
	mprintf(WHITE, "\n  < End Core Charge Definition <\n");
	
	coreChargesDefined = true;
}

bool setupWannier() {
	if(!g_bWannier) {
		mprintf(WHITE, "\n  > Wannier Centers >\n\n");
		g_bWannier = true;
		int watom = -1;
		int i;
		for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
			if(((CAtom *)g_oaAtoms[i])->m_bExclude) {
				watom = i;
				break;
			}
		}
		if(watom == -1)
			eprintf("    Atom label of Wannier centers not found.\n\n");
		else
			mprintf("    Atom type %s is excluded from the system, probably these are the Wannier centers.\n\n", ((CAtom *)g_oaAtoms[watom])->m_sName);
		
		bool ok = false;
		while(!ok) {
			CxString buf, buf2;
			
			buf.sprintf("    Which atom label do the wannier centers have (");
			for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
				buf2.sprintf("%s", ((CAtom *)g_oaAtoms[i])->m_sName);
				buf.strcat(buf2);
				if(i < g_oaAtoms.GetSize() - 2) {
					buf2.sprintf(", ");
					buf.strcat(buf2);
				}
			}
			buf2.sprintf(")? ");
			buf.strcat(buf2);
			if(watom != -1) {
				buf2.sprintf("[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
				buf.strcat(buf2);
			}
			
			if(watom == -1)
				AskString_ND(buf, &buf2);
			else
				AskString(buf, &buf2, ((CAtom *)g_oaAtoms[watom])->m_sName);
			
			for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
				if(mystricmp(buf2, ((CAtom *)g_oaAtoms[i])->m_sName) == 0) {
					g_iWannierAtomType = i;
					ok = true;
					break;
				}
			}
			if(!ok) {
				eprintf("    Wrong input.\n");
				inpprintf("! Wrong input.\n");
			}
		}
		
		g_fWannierCharge = fabsf(AskFloat("\n    Enter the negative charge of the Wannier centers (without sign): [2.0] ", 2.0f));
		parseCoreCharges();
		mprintf("\n");
		if(!g_TimeStep.ScanWannier(true)) {
			eprintf("\n    Setup of Wannier centers failed.\n");
			return false;
		}
		mprintf(WHITE, "\n  < End Wannier Centers <\n");
	}
	return true;
}

static void setupFluctuatingCharges()
{
	if ((!g_bReadChargesFrom4thXYZ) && (!g_bReadLAMMPSCharges))
	{
		if (g_bXYZ4thCol)
			g_bReadChargesFrom4thXYZ = true;

		if (g_bLAMMPSCharge)
			g_bReadLAMMPSCharges = true;

		mprintf(WHITE, "\n  > Fluctuating Atomic Partial Charges >\n");
		
		FILE *a;
		int z, z2, z3, z4;
		float tf;
		CMolecule *m;
		CSingleMolecule *sm;
		mprintf("\n    Trying to read charges from first time step...\n");
		a = fopen(g_sInputTraj,"rb");
		if (a == NULL)
		{
			eprintf("    Could not open \"%s\" for reading.\n",g_sInputTraj);
			abort();
		}
		
		if (!g_TimeStep.ReadTimestep(a,true))
		{
			eprintf("    Could not read first time step of \"%s\".\n",g_sInputTraj);
			abort();
		}
		
		fclose(a);
		
		mprintf("\n    Uniting molecules which have been broken by wrapping...\n");
		g_TimeStep.UniteMolecules(true);

		g_TimeStep.CalcCenters();
		
		mprintf("\n    Found the following atomic charges:\n\n");
		
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			m->m_bChargesAssigned = true;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
					continue;
				for (z3=0;z3<m->m_waAtomCount[z2];z3++)
				{
					for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
					{
						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
						mprintf("      %s[%d] %s%d:  %8.6f\n",m->m_sName,z4+1,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)]);
					}
				}
			}
		}
		mprintf("\n    Total charges of the molecules:\n\n");
		
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
				tf = 0;
				for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
				{
					if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
						continue;
					
					for (z3=0;z3<m->m_waAtomCount[z2];z3++)
						tf += g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)];
				}
				mprintf("      %s[%d]:  Total charge is %8.6f.\n",m->m_sName,z4+1,tf);
				if (z4 == 0)
					m->m_fCharge = tf;
			}
		}
		mprintf(WHITE, "\n  < End Fluctuating Atomic Partial Charges <\n");
	}
}


static void setupDipoleRestartFile() {
	if (!g_bLoadDipoleRestart) {
		mprintf(WHITE, "\n  > Dipole Restart File >\n\n");
		g_bLoadDipoleRestart = true;
		while (true) {
			CxString filename;
			AskString("    Enter name of dipole restart file: [dipole.restart] ", &filename, "dipole.restart");
			g_fDipoleRestartFile = fopen((const char *)filename, "r");
			if (g_fDipoleRestartFile == NULL) {
				eprintf("    Could not open \"%s\": %s.\n\n", (const char *)filename, strerror(errno));
				continue;
			}
			int numMolecules;
			fread(&numMolecules, sizeof(int), 1, g_fDipoleRestartFile);
			if (numMolecules != g_oaSingleMolecules.GetSize()) {
				eprintf("    This dipole restart file was written for a different number of molecules.\n\n");
				continue;
			}
			break;
		}
		mprintf(WHITE, "\n  < End Dipole Restart File <\n");
	}
}


static void setupMagneticDipoleRestartFile() {
	if (!g_bLoadMagneticDipoleRestart) {
		mprintf(WHITE, "\n  > Magnetic Moment Restart File >\n\n");
		g_bLoadMagneticDipoleRestart = true;
		while (true) {
			CxString filename;
			AskString("    Enter name of magnetic moment restart file: [magnetic.restart] ", &filename, "magnetic.restart");
			g_fMagneticDipoleRestartFile = fopen((const char *)filename, "r");
			if (g_fMagneticDipoleRestartFile == NULL) {
				eprintf("    Could not open \"%s\": %s.\n\n", (const char *)filename, strerror(errno));
				continue;
			}
			int numMolecules;
			fread(&numMolecules, sizeof(int), 1, g_fMagneticDipoleRestartFile);
			if (numMolecules != g_oaSingleMolecules.GetSize()) {
				eprintf("    This magnetic moment restart file was written for a different number of molecules.\n\n");
				continue;
			}
			break;
		}
		mprintf(WHITE, "\n  < End Magnetic Moment Restart File <\n");
	}
}


void ParseDipole()
{
	if (g_bDipoleDefined)
		return;
	
	g_bDipole = true;
	
	if (g_bTegri && (g_pTetraPak == NULL)) {
		try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
		if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { g_pTetraPak = new CTetraPak(); } catch(...) { g_pTetraPak = NULL; }
		if (g_pTetraPak == NULL) NewException((double)sizeof(CTetraPak),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pTetraPak->Parse();
	}
	
	mprintf(WHITE, "\n>>> Dipole Definition >>>\n\n");
	
	if(g_bAdvanced2) {
		g_bDipoleRefFixed = AskYesNo("    Use a box-fixed reference point for dipole calculation (y/n)? [no] ", false);
		mprintf("\n\n");
	} else {
		g_bDipoleRefFixed = false;
	}
	
	mprintf("    There are the following possibilities to provide dipole moments:\n");
	mprintf("    (1) Calculate dipole moments from Wannier centers\n        (These need to be atoms in the trajectory)\n");
	mprintf("    (2) Read dipole moments from an external file\n        (This file is expected to contain one line per timestep\n         with the three dipole vector components separated by spaces)\n");
	mprintf("    (3) Read dipole moments from the trajectory file\n        (These need to be specified in the comment line)\n");
	mprintf("    (4) Calculate dipole moments from fixed atomic partial charges\n");
	mprintf("    (5) Calculate dipole moments from fluctuating atomic partial charges\n        (These need to be in the fourth column of the trajectory)\n");
	mprintf("    (6) Calculate dipole moments from Voronoi partial charges\n        (This needs electron density data)\n");
	mprintf("    (7) Calculate Voronoi dipole moments\n        (This needs electron density data)\n");
	mprintf("    (8) Load a dipole restart file\n");
	mprintf("\n");
	
	mprintf("    The following dipole modes can be set for all molecules at once: 1, 5, 6, 7, 8\n");
	
	while (true) {
		int globalMode = AskRangeInteger("    Which dipole mode to set for all molecules at once? [none] ", 0, 8, 0);
		
		if (globalMode == 0) {
			int i;
			for (i = 0; i < g_oaMolecules.GetSize(); i++)
				((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 0;
			break;
		} else if (globalMode == 1) {
			if (setupWannier()) {
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 1;
				break;
			}
			continue;
		} else if (globalMode == 2) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 3) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 4) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 5) {
			if (g_bXYZ4thCol || g_bLAMMPSCharge) {
				setupFluctuatingCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 5;
				break;
			}
			eprintf("    The trajectory does not have a fourth column.\n");
			continue;
		} else if (globalMode == 6) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				parseCoreCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++) {
					unsigned char rty;
					ParseAtom("#2", i, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterIndex);
					((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 6;
				}
				break;
			}
			eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			continue;
		} else if (globalMode == 7) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				g_bVoroIntegrateDipoleMoment = true;
				parseCoreCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++) {
					unsigned char rty;
					ParseAtom("#2", i, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterIndex);
					((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 7;
				}
				break;
			}
			eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			continue;
		} else if (globalMode == 8) {
			setupDipoleRestartFile();
			int i;
			for (i = 0; i < g_oaMolecules.GetSize(); i++)
				((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 8;
			break;
		}
		eprintf("    This is impossible.\n");
	}
	
	while (true) {
		mprintf("\n    The following dipole modes are set up:\n\n");
		int longest = 0;
		int i;
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			int length = strlen(((CMolecule *)g_oaMolecules[i])->m_sName);
			if (length > longest)
				longest = length;
		}
		CxString buf;
		buf.sprintf("      %%%ds - %%d", longest);
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			mprintf(buf, ((CMolecule *)g_oaMolecules[i])->m_sName, ((CMolecule *)g_oaMolecules[i])->m_iDipoleMode);
			if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 0)
				mprintf(" (none)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 1)
				mprintf(" (Wannier centers)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 2)
				mprintf(" (external file)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 3)
				mprintf(" (trajectory comment line)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 4)
				mprintf(" (fixed atomic partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 5)
				mprintf(" (fluctuating atomic partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 6)
				mprintf(" (Voronoi partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 7)
				mprintf(" (Voronoi dipole moments)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 8)
				mprintf(" (restart file)\n");
			else
				mprintf("\n");
		}
		mprintf("\n");
		
		if (!AskYesNo("    Change dipole mode for a molecule (y/n)? [no] ", false))
			break;
		
		int mol;
		if (g_oaMolecules.GetSize() > 1) {
			buf.sprintf("    Change dipole mode for which molecule (");
			for (i = 0; i < g_oaMolecules.GetSize(); i++) {
				CxString buf2;
				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i + 1);
				buf.strcat(buf2);
				if(i < g_oaMolecules.GetSize() - 1) {
					buf.strcat(", ");
				}
			}
			buf.strcat(")? ");
			
			mol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
		} else {
			mol = 0;
		}
		
		int mode = AskRangeInteger_ND("    Which dipole mode to set for molecule %s? ", 1, 8, ((CMolecule *)g_oaMolecules[mol])->m_sName);
		
		if (mode == 0) {
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 0;
		} else if (mode == 1) {
			if (setupWannier()) {
				((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 1;
			}
		} else if (mode == 2) {
			FILE *dipoleFile = NULL;
			while (true) {
				if(dipoleFile != NULL) {
					fclose(dipoleFile);
					dipoleFile = NULL;
				}
				AskString_ND("    Name of file for dipole moments: ", &buf);
				mprintf("\n    Trying to read first line...\n\n");
				dipoleFile = fopen(buf, "r");
				if(dipoleFile == NULL) {
					eprintf("    Could not open \"%s\".\n\n", (const char*)buf);
					continue;
				}
				if (buf.fgets(1024, dipoleFile) == NULL) {
					eprintf("    Could not read first line.\n\n");
					continue;
				}
				float dipole[3];
				if(sscanf((const char *)buf, "%f %f %f", &dipole[0], &dipole[1], &dipole[2]) < 3) {
					eprintf("    Could not find three floating point numbers in first line.\n\n");
					continue;
				}
				mprintf("    Found the following dipole moment in the first line:\n");
				mprintf("      ( %8G | %8G | %8G ) Debye\n", dipole[0], dipole[1], dipole[2]);
				rewind(dipoleFile);
				break;
			}
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 2;
			((CMolecule *)g_oaMolecules[mol])->m_pDipoleFile = dipoleFile;
			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
		} else if (mode == 3) {
			CxString comment;
			comment.strcpy(g_TimeStep.m_pComment);
			mprintf("    The comment line is \"%s\".\n", (const char *)comment);
			buf.sprintf("");
			const char delim[] = "\t ";
			char *tok = strtok(comment.GetWritePointer(), delim);
			int num = 0;
			while(tok != NULL) {
				float f;
				if(sscanf(tok, "%f", &f) == 1) {
					CxString buf2;
					num++;
					if(num > 1) {
						buf2.sprintf(", %s (%d)", tok, num);
					} else {
						buf2.sprintf("%s (%d)", tok, num);
					}
					buf.strcat(buf2);
				}
				tok = strtok(NULL, delim);
			}
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[0] = AskRangeInteger("    Choose the X dipole component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[1] = AskRangeInteger("    Choose the Y dipole component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[2] = AskRangeInteger("    Choose the Z dipole component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 3;
			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
		} else if (mode == 4) {
			CMolecule *m = (CMolecule*)g_oaMolecules[mol];
			CSingleMolecule *sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			m->m_bChargesAssigned = true;
			int z2, z3, ti;
			CMolAtom *ma;
			CxFloatArray *tfa;
			bool b;
			float tf;
			double td;
			CxString buf2;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				mprintf("\n");
				ma = NULL;
				td = 1.0e50;
				
				try { tfa = new CxFloatArray("ParseDipole():tfa"); } catch(...) { tfa = NULL; }
				if (tfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				m->m_oaCharges.Add(tfa);
				_cnm:
				buf.sprintf("");
				buf2.sprintf("");
				b = false;
				ti = 0;
				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
				{
					if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
						continue;
					
					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
					
					if (ma->m_fAtomCode > td)
					{
						continue;
					}
					if (b)
					{
						if (ma->m_fAtomCode < td)
						{
							continue;
						}
						ti++;
						buf2.sprintf(", %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						buf.strcat(buf2);
						
					} else
					{
						if (ma->m_fAtomCode < td)
						{
							td = ma->m_fAtomCode;
						}
						
						buf.sprintf("    Enter partial atomic charge on %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						b = true;
						ti++;
					}
				}
				if (b)
				{
					buf.strcat(": [0.0] ");
					tf = AskFloat(buf,0);
					for (z3=0;z3<ti;z3++)
						tfa->Add(tf);
					td -= 1.0;
					goto _cnm;
				}
			}
			
			tf = 0;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				for (z3=0;z3<((CxFloatArray*)m->m_oaCharges[z2])->GetSize();z3++)
				{
					mprintf("    %s%d: %7.4f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3));
					tf += ((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3);
				}
			}
			m->m_fCharge = tf;
			mprintf(WHITE,"\n    Molecule %s: Total charge is %.4f.\n\n",m->m_sName,m->m_fCharge);
			
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 4;
		} else if (mode == 5) {
			if (g_bXYZ4thCol || g_bLAMMPSCharge) {
				setupFluctuatingCharges();
				((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 5;
			}
			eprintf("    The trajectory does not have a fourth column.\n");
		} else if (mode == 6) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				parseCoreCharges();
				unsigned char rty;
				ParseAtom("#2", mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex);
				((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 6;
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			}
		} else if (mode == 7) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				g_bVoroIntegrateDipoleMoment = true;
				parseCoreCharges();
				unsigned char rty;
				ParseAtom("#2", mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex);
				((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 7;
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			}
		} else if (mode == 8) {
			setupDipoleRestartFile();
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 8;
		} else {
			eprintf("    This is impossible.\n");
		}
	}
	
	if (!g_bDipoleRefFixed) {
		int i;
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			if (fabsf(((CMolecule *)g_oaMolecules[i])->m_fCharge) > 0.001f) {
				mprintf("\n");
				mprintf("    Molecule %s is charged (%.4f). Input of a dipole reference point is required.\n",((CMolecule *)g_oaMolecules[i])->m_sName, ((CMolecule *)g_oaMolecules[i])->m_fCharge);
				unsigned char rty;
				CxString buf;
				do {
					AskString("    Please enter dipole reference point: [#2] ", &buf, "#2");
				} while (!ParseAtom(buf, i, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterIndex));
			}
		}
	}
	mprintf("\n");
	
// 	CxString buf, buf2, comment;
// 
// //	const int BUF_SIZE = 1024;
// 	
// 	if (g_bDipoleDefined)
// 		return;
// 
// 	
// 	if(g_bAdvanced2) {
// 		g_bDipoleRefFixed = AskYesNo("    Use a box-fixed reference point for dipole calculation (y/n)? [no] ", false);
// 		mprintf("\n");
// 	} else {
// 		g_bDipoleRefFixed = false;
// 	}
// 	
// 	mprintf(WHITE, "\n>>> Dipole Definition >>>\n\n");
// 	mprintf("    There are different ways to provide dipole moments in TRAVIS:\n");
// 	mprintf("    (1) Calculate dipole moments from Wannier centers\n        (These need to be atoms in the trajectory)\n");
// 	mprintf("    (2) Read dipole moments from an external file\n        (This file is expected to contain one line per timestep\n         with the three dipole vector components separated by spaces)\n");
// 	mprintf("    (3) Read dipole moments from the trajectory file\n        (These need to be specified in the comment line)\n");
// 	mprintf("    (4) Calculate dipole moments from fixed atomic partial charges\n");
// 	mprintf("    (5) Calculate dipole moments from fluctuating atomic partial charges\n        (These need to be in the fourth column of the trajectory)\n");
// 	mprintf("    (8) Load dipole moments from a dipole restart file\n");
// 	mprintf("\n");
// 	
// 	while(true) {
// //		char buf[BUF_SIZE];
// //		char buf2[BUF_SIZE];
// //		size_t remaining = BUF_SIZE;
// 		int mol;
// 		if(g_oaMolecules.GetSize() > 1) {
// /*#ifdef TARGET_LINUX
// 			remaining -= snprintf(buf, remaining, "    Add dipole information for which molecule (");
// #else
// 			remaining -= sprintf(buf, "    Add dipole information for which molecule (");
// #endif*/
// 
// 			buf.sprintf("    Add dipole information for which molecule (");
// 
// 			int i;
// 			for(i = 0; i < g_oaMolecules.GetSize(); i++) {
// 
// /*				if(remaining < 1)
// 					break;
// #ifdef TARGET_LINUX
// 				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// #else
// 				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// #endif
// 				strncat(buf, buf2, remaining - 1);
// 				remaining -= length;*/
// 
// 				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// 				buf.strcat(buf2);
// 
// 				if(i < g_oaMolecules.GetSize() - 1) {
// /*#ifdef TARGET_LINUX
// 					length = snprintf(buf2, remaining, ", ");
// #else
// 					length = sprintf(buf2, ", ");
// #endif
// 					strncat(buf, buf2, remaining - 1);
// 					remaining -= length;*/
// 
// 					buf2.sprintf(", ");
// 					buf.strcat(buf2);
// 				}
// 			}
// //			strncat(buf, ")? ", remaining - 1);
// 			buf.strcat(")? ");
// 
// 			mol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
// 		} else {
// 			mol = 0;
// 		}
// 		
// 		int mode;
// 
// 		{
// 			if(g_bXYZ4thCol)
// 				mode = AskRangeInteger("\n    Use Wannier centers (1), external file (2), comment line (3), fixed charges (4), or fluctuating charges (5) for molecule %s? [1] ", 1, 5, 1, ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			else
// 				mode = AskRangeInteger("\n    Use Wannier centers (1), external file (2), comment line (3), or fixed charges (4) for molecule %s? [1] ", 1, 4, 1, ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		}
// 		mprintf("\n");
// 		
// 		if(mode == 1) {
// 			if(!g_bWannier) {
// 				mprintf(WHITE, "  > Wannier Centers >\n\n");
// 				g_bWannier = true;
// 				int watom = -1;
// 				int i;
// 				for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
// 					if(((CAtom *)g_oaAtoms[i])->m_bExclude) {
// 						watom = i;
// 						break;
// 					}
// 				}
// 				if(watom == -1)
// 					eprintf("    Atom label of Wannier centers not found.\n\n");
// 				else
// 					mprintf("    Atom type %s is excluded from the system, probably these are the Wannier centers.\n\n", ((CAtom *)g_oaAtoms[watom])->m_sName);
// 				
// 				bool ok = false;
// 				while(!ok) {
// /*					remaining = BUF_SIZE;
// #ifdef TARGET_LINUX
// 					remaining -= snprintf(buf, remaining, "    Which atom label do the wannier centers have (");
// #else
// 					remaining -= sprintf(buf, "    Which atom label do the wannier centers have (");
// #endif*/
// 
// 					buf.sprintf("    Which atom label do the wannier centers have (");
// 
// 					for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
// 
// /*						if(remaining < 1)
// 							break;
// #ifdef TARGET_LINUX
// 						size_t length = snprintf(buf2, remaining, "%s", ((CAtom *)g_oaAtoms[i])->m_sName);
// #else
// 						size_t length = sprintf(buf2, "%s", ((CAtom *)g_oaAtoms[i])->m_sName);
// #endif
// 						strncat(buf, buf2, remaining - 1);
// 						remaining -= length;*/
// 
// 						buf2.sprintf("%s", ((CAtom *)g_oaAtoms[i])->m_sName);
// 						buf.strcat(buf2);
// 
// 						if(i < g_oaAtoms.GetSize() - 2) {
// /*#ifdef TARGET_LINUX
// 							length = snprintf(buf2, remaining, ", ");
// #else
// 							length = sprintf(buf2, ", ");
// #endif
// 							strncat(buf, buf2, remaining - 1);
// 							remaining -= length;*/
// 
// 							buf2.sprintf(", ");
// 							buf.strcat(buf2);
// 						}
// 					}
// 
// /*#ifdef TARGET_LINUX
// 					size_t length = snprintf(buf2, remaining, ")? ");
// #else
// 					size_t length = sprintf(buf2, ")? ");
// #endif
// 					strncat(buf, buf2, remaining - 1);
// 					remaining -= length;*/
// 
// 					buf2.sprintf(")? ");
// 					buf.strcat(buf2);
// 
// 					if(watom != -1) {
// /*#ifdef TARGET_LINUX
// 						snprintf(buf2, remaining, "[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
// #else
// 						sprintf(buf2, "[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
// #endif
// 						strncat(buf, buf2, remaining - 1);*/
// 
// 						buf2.sprintf("[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
// 						buf.strcat(buf2);
// 					}
// 					
// 					if(watom == -1)
// 						AskString_ND(buf, &buf2);
// 					else
// 						AskString(buf, &buf2, ((CAtom *)g_oaAtoms[watom])->m_sName);
// 					
// 					for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
// 						if(mystricmp(buf2, ((CAtom *)g_oaAtoms[i])->m_sName) == 0) {
// 							g_iWannierAtomType = i;
// 							ok = true;
// 							break;
// 						}
// 					}
// 					if(!ok) {
// 						eprintf("    Wrong input.\n");
// 						inpprintf("! Wrong input.\n");
// 					}
// 				}
// 				
// 				g_fWannierCharge = fabsf(AskFloat("\n    Enter the negative charge of the Wannier centers (without sign): [2.0] ", 2.0f));
// 				for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
// 					if(i == g_iWannierAtomType)
// 						continue;
// 					if(i == g_iVirtAtomType)
// 						continue;
// 					float def;
// 					if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "B") == 0) def = 3.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "C") == 0) def = 4.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "N") == 0) def = 5.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "O") == 0) def = 6.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "F") == 0) def = 7.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Si") == 0) def = 4.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "P") == 0) def = 5.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "S") == 0) def = 6.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Cl") == 0) def = 7.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Br") == 0) def = 7.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "I") == 0) def = 7.0f;
// 					else def = (float)((CAtom *)g_oaAtoms[i])->m_pElement->m_iOrd;
// 					((CAtom *)g_oaAtoms[i])->m_fCharge = AskFloat("    Enter core charge (mind pseudopotentials!) for atom type %s: [%.1f] ", def, ((CAtom*)g_oaAtoms[i])->m_sName, def);
// 				}
// 				mprintf("\n");
// 				if(!g_TimeStep.ScanWannier(true)) {
// 					eprintf("\n    No dipole information for molecule %s added.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 				}
// 				mprintf(WHITE, "\n  < End Wannier Centers <\n\n");
// 			}
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 1;
// 			mprintf("    Set dipole information for %s to Wannier centers.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 2) {
// 			FILE *dipoleFile = NULL;
// 			while(true) {
// 				if(dipoleFile != NULL) {
// 					fclose(dipoleFile);
// 					dipoleFile = NULL;
// 				}
// 				AskString_ND("    Name of file for dipole moments: ", &buf);
// 				mprintf("\n    Trying to read first line...\n\n");
// 				dipoleFile = fopen(buf, "r");
// 				if(dipoleFile == NULL) {
// 					eprintf("    Could not open \"%s\".\n\n", (const char*)buf);
// 					continue;
// 				}
// //				if(fgets(buf, BUF_SIZE, dipoleFile) == NULL) {
// 				if (buf.fgets(1024, dipoleFile) == NULL) {
// 					eprintf("    Could not read first line.\n\n");
// 					continue;
// 				}
// 				float dipole[3];
// 				if(sscanf(buf, "%f %f %f", &dipole[0], &dipole[1], &dipole[2]) < 3) {
// 					eprintf("    Could not find three floating point numbers in first line.\n\n");
// 					continue;
// 				}
// 				mprintf("    Found the following dipole moment in the first line:\n");
// 				mprintf("      ( %8G | %8G | %8G ) Debye\n", dipole[0], dipole[1], dipole[2]);
// 				rewind(dipoleFile);
// 				break;
// 			}
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 2;
// 			mprintf("\n    Set dipole information for %s to external file.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			((CMolecule *)g_oaMolecules[mol])->m_pDipoleFile = dipoleFile;
// 			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
// 				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 3) {
// //			char comment[BUF_SIZE];
// //			strncpy(comment, g_TimeStep.m_pComment, BUF_SIZE);
// //			comment[BUF_SIZE-1] = '\0';
// 			comment.strcpy(g_TimeStep.m_pComment);
// 			mprintf("    The comment line is \"%s\".\n", (const char*)comment);
// //			char buf[BUF_SIZE];
// //			buf[0] = '\0';
// //			size_t remaining = BUF_SIZE;
// 			buf.sprintf("");
// 			const char delim[] = "\t ";
// 			char *tok = strtok(comment.GetWritePointer(), delim);
// 			int num = 0;
// 			while(tok != NULL) {
// 				float f;
// 				if(sscanf(tok, "%f", &f) == 1) {
// 					num++;
// //					char buf2[BUF_SIZE];
// //					size_t length;
// 					if(num > 1) {
// /*#ifdef TARGET_LINUX
// 						length = snprintf(buf2, BUF_SIZE, ", %s (%d)", tok, num);
// #else
// 						length = sprintf(buf2, ", %s (%d)", tok, num);
// #endif*/
// 						buf2.sprintf(", %s (%d)", tok, num);
// 					} else {
// /*#ifdef TARGET_LINUX
// 						length = snprintf(buf2, BUF_SIZE, "%s (%d)", tok, num);
// #else
// 						length = sprintf(buf2, "%s (%d)", tok, num);
// #endif*/
// 						buf2.sprintf("%s (%d)", tok, num);
// 					}
// //					strncat(buf, buf2, remaining - 1);
// 					buf.strcat(buf2);
// /*					remaining -= length;
// 					if(remaining < 1)
// 						break;*/
// 				}
// 				tok = strtok(NULL, delim);
// 			}
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[0] = AskRangeInteger("    Choose the X dipole component: %s [1] ", 1, num, 1, (const char*)buf) - 1;
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[1] = AskRangeInteger("    Choose the Y dipole component: %s [1] ", 1, num, 1, (const char*)buf) - 1;
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[2] = AskRangeInteger("    Choose the Z dipole component: %s [1] ", 1, num, 1, (const char*)buf) - 1;
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 3;
// 			mprintf("\n    Set dipole information for %s to comment line.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
// 				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 4) {
// 			CMolecule *m = (CMolecule*)g_oaMolecules[mol];
// 			CSingleMolecule *sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
// 			m->m_bChargesAssigned = true;
// 			int z2, z3, ti;
// 			CMolAtom *ma;
// 			CxFloatArray *tfa;
// 			bool b;
// 			float tf;
// 			double td;
// 			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 			{
// 				mprintf("\n");
// //					mprintf("%s Anfang.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 				ma = NULL;
// 				td = 1.0e50;
// 
// 				try { tfa = new CxFloatArray("ParseDipole():tfa"); } catch(...) { tfa = NULL; }
// 				if (tfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
// 							
// 				m->m_oaCharges.Add(tfa);
// _cnm:
// //				buf[0] = 0;
// //				buf2[0] = 0;
// 				buf.sprintf("");
// 				buf2.sprintf("");
// 				b = false;
// 				ti = 0;
// //					mprintf("D.\n");
// 				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
// 				{
// //						mprintf("z3=%d.\n",z3);
// 					if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
// 						continue;
// 
// 					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
// 
// //						mprintf("AtomType %s, z3=%d, td=%f, ac=%f.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3,td,ma->m_fAtomCode);
// 
// 					if (ma->m_fAtomCode > td)
// 					{
// //							mprintf("  AtomCode = %f > %f = td. Skip.\n",ma->m_fAtomCode,td);
// 						continue;
// 					}
// 					if (b)
// 					{
// 						if (ma->m_fAtomCode < td)
// 						{
// //								mprintf("  AtomCode = %f < %f = td. Skip.\n",ma->m_fAtomCode,td);
// 							continue;
// 						}
// 						ti++;
// //							mprintf("  Taking %s%d into account.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 
// //						sprintf(buf2,", %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// //						strcat(buf,buf2);
// 
// 						buf2.sprintf(", %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						buf.strcat(buf2);
// 
// 					} else
// 					{
// 						if (ma->m_fAtomCode < td)
// 						{
// //								mprintf("  Changing td from %f to %f.\n",td,ma->m_fAtomCode);
// //								mprintf("  Taking %s%d into account.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 							td = ma->m_fAtomCode;
// 						}
// 
// //						sprintf(buf,"    Enter partial atomic charge on %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						buf.sprintf("    Enter partial atomic charge on %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						b = true;
// 						ti++;
// 					}
// 				}
// //					mprintf("%s Mitte.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 				if (b)
// 				{
// //					strcat(buf,": [0.0] ");
// 					buf.strcat(": [0.0] ");
// 					tf = AskFloat(buf,0);
// 					for (z3=0;z3<ti;z3++)
// 						tfa->Add(tf);
// 					td -= 1.0;
// 					goto _cnm;
// 				}
// //					mprintf("%s Ende.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 			}
// 
// //			mprintf("\n");
// 			tf = 0;
// 			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 			{
// 				for (z3=0;z3<((CxFloatArray*)m->m_oaCharges[z2])->GetSize();z3++)
// 				{
// 					mprintf("    %s%d: %7.4f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3));
// 					tf += ((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3);
// 				}
// 			}
// 			m->m_fCharge = tf;
// 			mprintf(WHITE,"\n    Molecule %s: Total charge is %.4f.\n\n",m->m_sName,m->m_fCharge);
// 			
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 4;
// 			mprintf("\n    Set dipole information for %s to fixed atomic partial charges.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 5) {
// 			if(!g_bReadChargesFrom4thXYZ) {
// 				mprintf(WHITE, "  > Fluctuating Atomic Partial Charges >\n");
// 				g_bReadChargesFrom4thXYZ = true;
// 				
// 				FILE *a;
// 				int z, z2, z3, z4;
// 				float tf;
// 				CMolecule *m;
// 				CSingleMolecule *sm;
// 				mprintf("\n    Trying to read charges from first time step...\n");
// 				a = fopen(g_sInputTraj,"rb");
// 				if (a == NULL)
// 				{
// 					eprintf("    Could not open \"%s\" for reading.\n",g_sInputTraj);
// 					abort();
// 				}
// 
// 				if (!g_TimeStep.ReadTimestep(a,true))
// 				{
// 					eprintf("    Could not read first time step of \"%s\".\n",g_sInputTraj);
// 					abort();
// 				}
// 
// 				fclose(a);
// 
// 				mprintf("\n    Uniting molecules which have been broken by wrapping...\n");
// 				g_TimeStep.UniteMolecules(true);
// 
// 				mprintf("\n    Found the following atomic charges:\n\n");
// 
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					m = (CMolecule*)g_oaMolecules[z];
// 					m->m_bChargesAssigned = true;
// 					for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 					{
// 						if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
// 							continue;
// 						for (z3=0;z3<m->m_waAtomCount[z2];z3++)
// 						{
// 							for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
// 							{
// 								sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
// 								mprintf("      %s[%d] %s%d:  %8.6f\n",m->m_sName,z4+1,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)]);
// 							}
// 						}
// 					}
// 				}
// 				mprintf("\n    Total charges of the molecules:\n\n");
// 
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					m = (CMolecule*)g_oaMolecules[z];
// 					for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
// 					{
// 						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
// 						tf = 0;
// 						for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 						{
// 							if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
// 								continue;
// 
// 							for (z3=0;z3<m->m_waAtomCount[z2];z3++)
// 								tf += g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)];
// 						}
// 						mprintf("      %s[%d]:  Total charge is %8.6f.\n",m->m_sName,z4+1,tf);
// 						if (z4 == 0)
// 							m->m_fCharge = tf;
// 					}
// 				}
// 				mprintf(WHITE, "\n  < End Fluctuating Atomic Partial Charges <\n\n");
// 			}
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 5;
// 			mprintf("    Set dipole information for %s to fluctuating atomic partial charges.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 6) {
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 6;
// 			mprintf("    Set dipole information for %s to Voronoi partial charges.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			unsigned char rty;
// 			ParseAtom("#2", mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex);
// 		} else if (mode == 7) {
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 7;
// 			mprintf("    Set dipole information for %s to Voronoi dipole moments.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			unsigned char rty;
// 			ParseAtom("#2", mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex);
// 		} else if (mode == 8) {
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 8;
// // 			g_bDipoleRestartFile = fopen("dipole.restart", "r");
// 			mprintf("    Using dipole restart file.\n");
// 		} else {
// 			eprintf("Weird error.\n");
// 			abort();
// 		}
// 		
// 		if(!g_bDipoleRefFixed) {
// 			if(fabsf(((CMolecule *)g_oaMolecules[mol])->m_fCharge) > 0.001f) {
// 				mprintf("\n");
// 				mprintf("    The molecule is charged (%.4f). Input of a dipole reference point is required.\n", ((CMolecule *)g_oaMolecules[mol])->m_fCharge);
// 				unsigned char rty;
// 				do {
// 					AskString("    Please enter dipole reference point: [#2] ", &buf, "#2");
// 				} while(!ParseAtom(buf, mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex));
// 			}
// 		}
// 		
// 		if(!(g_oaMolecules.GetSize() > 1) || !AskYesNo("\n    Add dipole information for another molecule (y/n)? [yes] ", true))
// 			break;
// 		mprintf("\n");
// 	}
// 	mprintf("\n");
	
	size_t nameLength = 0;
	int i;
	for(i = 0; i < g_oaMolecules.GetSize(); i++) {
		CMolecule *m = (CMolecule *)g_oaMolecules[i];
		if(m->m_iDipoleMode == 1 || m->m_iDipoleMode == 4 || m->m_iDipoleMode == 5)
			if(strlen(m->m_sName) > nameLength)
				nameLength = strlen(m->m_sName);
	}
	
	if(!g_bDipoleRefFixed) {
		for(i = 0; i < g_oaMolecules.GetSize(); i++) {
			CMolecule *m = (CMolecule *)g_oaMolecules[i];
			if(m->m_iDipoleMode == 1 || m->m_iDipoleMode == 4 || m->m_iDipoleMode == 5) {
				mprintf(WHITE, "  - %s:", m->m_sName);
				size_t j;
				for(j = strlen(m->m_sName); j <= nameLength; j++)
					mprintf(" ");
				mprintf(" Dipole reference point is %s%d.\n", ((CAtom*)g_oaAtoms[m->m_baAtomIndex[m->m_iDipoleCenterType]])->m_sName, m->m_iDipoleCenterIndex + 1);
			}
		}
	}
	
	mprintf("\n    Calculating dipole moments from first time step...\n\n");
	if (g_bTegri) {
		g_pTetraPak->ProcessStep(&g_TimeStep, true);
		mprintf("\n");
	}
	g_TimeStep.CalcDipoles(true);
	for(i = 0; i < g_oaMolecules.GetSize(); i++)
		if(((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 2)
			rewind(((CMolecule *)g_oaMolecules[i])->m_pDipoleFile);
	if (g_bLoadDipoleRestart) {
		fseek(g_fDipoleRestartFile, sizeof(int), SEEK_SET);
	}
	
	for(i = 0; i < g_oaMolecules.GetSize(); i++) {
		CMolecule *m = (CMolecule *)g_oaMolecules[i];
		if(m->m_iDipoleMode != 0) {
			float min = 1.0e6f, max = 0.0f, ave = 0.0f;
			int j;
			for(j = 0; j < m->m_laSingleMolIndex.GetSize(); j++) {
				float val = ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[j]])->m_vDipole.GetLength();
				if(val > max)
					max = val;
				if(val < min)
					min = val;
				ave += val;
			}
			ave /= (float)m->m_laSingleMolIndex.GetSize();
			mprintf(WHITE, "  - %s:", m->m_sName);
			size_t k;
			for(k = strlen(m->m_sName); k <= nameLength; k++)
				mprintf(" ");
			mprintf(" Min. %12G, Max. %12G, Avg. %12G Debye.\n", min, max, ave);
			mprintf("        Cartesian dipole vector of first molecule is ( %8G | %8G | %8G ) Debye.\n\n", ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[0], ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[1], ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[2]);
		}
	}
	mprintf("\n");
	
	if (g_bAdvanced2)
	{
		int z, z2, ti;
		CxIntArray *tia;
//		char buf[256];
		CxString buf;
		if (AskYesNo("    Write out cartesian dipole vectors of molecules in each step (y/n)? [no] ",false))
		{
			g_bDumpDipoleVector = true;
			mprintf("\n");
			ti = 0;
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				if (AskYesNo("      Consider %s molecules for writing out dipole vectors (y/n)? [yes] ",true,((CMolecule*)g_oaMolecules[z])->m_sName))
				{
					try { tia = new CxIntArray("ParseDipole():tia"); } catch(...) { tia = NULL; }
					if (tia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					g_oaDumpDipoleVector.Add(tia);
					AskString("        Which representants to consider (e.g. 1,3-6,8)? [all] ",&buf,"*");
					if (buf[0] == '*')
					{
						for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
							tia->Add(z2+1);
					} else ParseIntList(buf,tia);
					for (z2=0;z2<tia->GetSize();z2++)
						tia->GetAt(z2)--;
					ti += tia->GetSize();
					mprintf("\n");
				} else g_oaDumpDipoleVector.Add(NULL);
			}
			g_iDumpDipoleSMCount = ti;
			g_bDumpDipoleAbs = AskYesNo("    Write out absolute value as 4th column per molecule (y/n)? [yes] ",true);
			g_bDumpDipoleXYZ = AskYesNo("    Write out XYZ trajectory to visualize dipole vectors (y/n)? [no] ",false);

			if (g_bDumpDipoleXYZ)
			{
				g_fDumpDipoleScale = AskFloat("    Please enter dipole scale in trajectory (unit: Angstrom/Debye): [1.0] ",1.0f);
				g_iDumpDipoleXYZAtoms = 0;
				for (z=0;z<g_oaMolecules.GetSize();z++)
				{
					if (g_oaDumpDipoleVector[z] == NULL)
						continue;
					g_iDumpDipoleXYZAtoms += ((CMolecule*)g_oaMolecules[z])->m_iAtomGesNoVirt * ((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();
				}
				g_iDumpDipoleXYZAtoms += g_iDumpDipoleSMCount*2;
			}
		} else g_bDumpDipoleVector = false;
		mprintf("\n");
	}

	if (g_bTegri && (g_iTrajFormat != 5))
	{
		mprintf("    Rewinding density cube file to beginning...\n\n");
		fseek(g_pTetraPak->m_fCubePipe,0,SEEK_SET);
	}
	
	mprintf(WHITE, "<<< End of Dipole Definition <<<\n\n\n");
	
	g_bDipoleDefined = true;
	
//----------------- OLD -----------------------
// 	
// 	int z, z2, z3, z4, ti;
// 	bool b;
// 	char buf[256], buf2[256];
// 	double td, td2, td3;
// 	float tf;
// 	unsigned char rty;
// 	CMolecule *m;
// 	CSingleMolecule *sm;
// 	CMolAtom *ma;
// 	CxFloatArray *tfa;
// 	CxIntArray *tia;
// 	FILE *a;
// 
// 	if (g_bDipoleDefined)
// 		return;
// 
// _dipbeg:
// 	mprintf(WHITE,"\n>>> Dipole definition >>>\n\n");
// 	mprintf("    TRAVIS can calculate dipole moments either from wannier centers (you need\n");
// 	mprintf("    to have those in the trajectory then) or from fixed atomic partial charges.\n\n");
// 	g_bWannier = AskYesNo("    Obtain dipole vectors from wannier centers (y/n)? [yes] ",true);
// 	mprintf("\n");
// 	if (g_bWannier)
// 	{
// _wantype:
// 		ti = -1;
// 		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
// 		{
// 			if (((CAtom*)g_oaAtoms[z])->m_bExclude)
// 			{
// 				ti = z;
// 				break;
// 			}
// 		}
// 
// 		if (ti == -1)
// 			eprintf("    Atom label of wannier centers not found.\n\n");
// 				else mprintf("    Atom type %s is excluded from system, probably the wannier centers.\n\n",((CAtom*)g_oaAtoms[ti])->m_sName);
// 		
// 		sprintf(buf,"    Which atom label do the wannier centers have (");
// 
// 		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
// 		{
// 			if (z < (int)g_oaAtoms.GetSize()-2)
// 			{
// 				sprintf(buf2,"%s, ",((CAtom*)g_oaAtoms[z])->m_sName);
// 				strcat(buf,buf2);
// 			} else 
// 			{
// 				if (ti == -1)
// 					sprintf(buf2,"%s)? ",((CAtom*)g_oaAtoms[z])->m_sName);
// 						else sprintf(buf2,"%s)? [%s] ",((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[ti])->m_sName);
// 				strcat(buf,buf2);
// 			}
// 		}
// 
// 		if (ti == -1)
// 			AskString_ND(buf,buf2);
// 				else AskString(buf,buf2,((CAtom*)g_oaAtoms[ti])->m_sName);
// 
// 		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
// 		{
// 			if (mystricmp(buf2,((CAtom*)g_oaAtoms[z])->m_sName)==0)
// 			{
// 				g_iWannierAtomType = z;
// 				goto _wandone;
// 			}
// 		}
// 		eprintf("    Wrong input.\n");
// 		inpprintf("! Wrong input.\n");
// 		goto _wantype;
// _wandone:
// 		g_fWannierCharge = (float)fabs(AskFloat("    Enter the negative charge of the wannier centers (without sign): [2.0] ",2.0f));
// 		mprintf("\n");
// 		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
// 		{
// 			if (z == g_iWannierAtomType)
// 				continue;
// 			if (z == g_iVirtAtomType)
// 				continue;
// 
// 			// Pseudopotential core charges for frequently used atoms
// 			if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"C") == 0)
// 				z2 = 4;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"N") == 0)
// 				z2 = 5;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"O") == 0)
// 				z2 = 6;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"F") == 0)
// 				z2 = 7;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"P") == 0)
// 				z2 = 5;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"S") == 0)
// 				z2 = 6;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"Cl") == 0)
// 				z2 = 7;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"B") == 0)
// 				z2 = 3;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"Br") == 0)
// 				z2 = 7;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"I") == 0)
// 				z2 = 7;
// 			else z2 = ((CAtom*)g_oaAtoms[z])->m_pElement->m_iOrd;
// 
// 			((CAtom*)g_oaAtoms[z])->m_fCharge = AskFloat("    Enter core charge (mind pseudopotentials!) for atom type %s: [%d.0] ",(float)z2,((CAtom*)g_oaAtoms[z])->m_sName,z2);
// 		}
// 		mprintf("\n");
// 		if (!g_TimeStep.ScanWannier(true))
// 			goto _dipbeg;
// /*		if (g_bAdvanced2)
// 		{
// 			g_bUnwrapWannier = AskYesNo("    Write out trajectory with unwrapped molecules and wannier centers (y/n)? [no] ",false);
// 		} else g_bUnwrapWannier = false;*/
// 	} else // NOT WANNIER
// 	{
// 		if (g_bBetaFeatures && g_bXYZ4thCol)
// 		{
// 			if (AskYesNo("    Read atomic charges from 4th numeric column of XYZ file (y/n)? [yes] ",true))
// 			{
// 				g_bReadChargesFrom4thXYZ = true;
// 
// 				mprintf("\n    Trying to read charges from first time step...\n");
// 				a = fopen(g_sInputTraj,"rb");
// 				if (a == NULL)
// 				{
// 					eprintf("    Could not open \"%s\" for reading.\n",g_sInputTraj);
// 					abort();
// 				}
// 
// 				if (!g_TimeStep.ReadTimestep(a,true))
// 				{
// 					eprintf("    Could not read first time step of \"%s\".\n",g_sInputTraj);
// 					abort();
// 				}
// 
// 				fclose(a);
// 
// 				mprintf("\n    Uniting molecules which have been broken by wrapping...\n");
// 				g_TimeStep.UniteMolecules(true);
// 
// 				mprintf("\n    Found the following atomic charges:\n\n");
// 
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					m = (CMolecule*)g_oaMolecules[z];
// 					m->m_bChargesAssigned = true;
// 					for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 					{
// 						if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
// 							continue;
// 						for (z3=0;z3<m->m_waAtomCount[z2];z3++)
// 						{
// 							for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
// 							{
// 								sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
// 								mprintf("      %s[%d] %s%d:  %8.6f\n",m->m_sName,z4+1,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)]);
// 							}
// 						}
// 					}
// 				}
// 				mprintf("\n    Total charges of the molecules:\n\n");
// 
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					m = (CMolecule*)g_oaMolecules[z];
// 					for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
// 					{
// 						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
// 						tf = 0;
// 						for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 						{
// 							if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
// 								continue;
// 
// 							for (z3=0;z3<m->m_waAtomCount[z2];z3++)
// 								tf += g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)];
// 						}
// 						mprintf("      %s[%d]:  Total charge is %8.6f.\n",m->m_sName,z4+1,tf);
// 						if (z4 == 0)
// 							m->m_fCharge = tf;
// 					}
// 				}
// 				mprintf("\n");
// 
// 				goto _chargedone;
// 			}
// 		}
// 		for (z=0;z<g_oaMolecules.GetSize();z++)
// 		{
// 			m = (CMolecule*)g_oaMolecules[z];
// 			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
// 			mprintf(YELLOW,"\n  * Molecule %s\n\n",m->m_sName);
// 
// 			m->m_bChargesAssigned = AskYesNo("    Enter partial atomic charges for this molecule (y/n)? [yes] ",true);
// 
// 			if (!m->m_bChargesAssigned)
// 			{
// 				mprintf("\n    Dipole moment will be always zero for this molecule.\n");
// 				goto _nocharge;
// 			}
// 
// 			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 			{
// 				mprintf("\n");
// //					mprintf("%s Anfang.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 				ma = NULL;
// 				td = 1.0e50;
// 
// 				try { tfa = new CxFloatArray("ParseDipole():tfa"); } catch(...) { tfa = NULL; }
// 				if (tfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
// 							
// 				m->m_oaCharges.Add(tfa);
// _cnm:
// 				buf[0] = 0;
// 				buf2[0] = 0;
// 				b = false;
// 				ti = 0;
// //					mprintf("D.\n");
// 				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
// 				{
// //						mprintf("z3=%d.\n",z3);
// 					if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
// 						continue;
// 
// 					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
// 
// //						mprintf("AtomType %s, z3=%d, td=%f, ac=%f.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3,td,ma->m_fAtomCode);
// 
// 					if (ma->m_fAtomCode > td)
// 					{
// //							mprintf("  AtomCode = %f > %f = td. Skip.\n",ma->m_fAtomCode,td);
// 						continue;
// 					}
// 					if (b)
// 					{
// 						if (ma->m_fAtomCode < td)
// 						{
// //								mprintf("  AtomCode = %f < %f = td. Skip.\n",ma->m_fAtomCode,td);
// 							continue;
// 						}
// 						ti++;
// //							mprintf("  Taking %s%d into account.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						sprintf(buf2,", %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						strcat(buf,buf2);
// 					} else
// 					{
// 						if (ma->m_fAtomCode < td)
// 						{
// //								mprintf("  Changing td from %f to %f.\n",td,ma->m_fAtomCode);
// //								mprintf("  Taking %s%d into account.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 							td = ma->m_fAtomCode;
// 						}
// 						sprintf(buf,"    Enter partial atomic charge on %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						b = true;
// 						ti++;
// 					}
// 				}
// //					mprintf("%s Mitte.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 				if (b)
// 				{
// 					strcat(buf,": [0.0] ");
// 					tf = AskFloat(buf,0);
// 					for (z3=0;z3<ti;z3++)
// 						tfa->Add(tf);
// 					td -= 1.0;
// 					goto _cnm;
// 				}
// //					mprintf("%s Ende.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 			}
// 
// //			mprintf("\n");
// 			tf = 0;
// 			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 			{
// 				for (z3=0;z3<((CxFloatArray*)m->m_oaCharges[z2])->GetSize();z3++)
// 				{
// 					mprintf("    %s%d: %7.4f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3));
// 					tf += ((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3);
// 				}
// 			}
// 			m->m_fCharge = tf;
// 			mprintf(WHITE,"\n    Molecule %s: Total charge is %.4f.\n\n",m->m_sName,m->m_fCharge);
// _nocharge:;
// 		}
// 	}
// _chargedone:
// 	if (g_bAdvanced2)
// 	{
// 		mprintf("\n");
// 		g_bDipoleRefFixed = AskYesNo("    Use a box-fixed reference point for dipole calculation (y/n)? [no] ",false);
// 		if (g_bDipoleRefFixed)
// 			mprintf("\n    Warning: Absolute dipole moments of charged molecules will be erroneous.\n    Only the derivatives will be useful.\n\n");
// 	} else g_bDipoleRefFixed = false;
// 
// 	if (!g_bDipoleRefFixed)
// 	{
// 		for (z=0;z<g_oaMolecules.GetSize();z++)
// 		{
// 			m = (CMolecule*)g_oaMolecules[z];
// 			if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 				continue;
// 
// 			if (fabs(m->m_fCharge) > 0.001)
// 			{
// 				mprintf("\n");
// 				mprintf("    Molecule %s is charged (%.4f). Input of a dipole reference point is required.\n",m->m_sName,m->m_fCharge);
// 	_diprefagain:
// 				AskString("    Please enter dipole reference point for %s: [#1] ",buf,"#1",m->m_sName);
// 				if (!ParseAtom(buf,z,m->m_iDipoleCenterType,rty,m->m_iDipoleCenterIndex))
// 					goto _diprefagain;
// 			}
// 		}
// 		mprintf("\n");
// 	}
// 
// 	ti = 0;
// 	for (z=0;z<g_oaMolecules.GetSize();z++)
// 	{
// 		m = (CMolecule*)g_oaMolecules[z];
// 		if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 			continue;
// 		if (((int)strlen(m->m_sName)) > ti)
// 			ti = strlen(m->m_sName);
// 	}
// 
// 	if (!g_bDipoleRefFixed)
// 	{
// 		for (z=0;z<g_oaMolecules.GetSize();z++)
// 		{
// 			m = (CMolecule*)g_oaMolecules[z];
// 			if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 				continue;
// 			mprintf(WHITE,"  - %s:",m->m_sName);
// 			for (z2=strlen(m->m_sName);z2<=ti;z2++)
// 				mprintf(" ");
// 			mprintf(" Dipole reference point is %s%d.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[m->m_iDipoleCenterType]])->m_sName,m->m_iDipoleCenterIndex+1);
// 		}
// 		mprintf("\n");
// 	}
// 
// 	mprintf("    Calculating dipole moments from first time step...\n\n");
// 	g_TimeStep.CalcDipoles();
// 
// 	sm = NULL;
// 	for (z=0;z<g_oaMolecules.GetSize();z++)
// 	{
// 		m = (CMolecule*)g_oaMolecules[z];
// 		if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 			continue;
// 		td = 0;
// 		td2 = 1.0e10;
// 		td3 = 0;
// 		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
// 		{
// 			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
// 			tf = sm->m_vDipole.GetLength();
// 			if (tf > td)
// 				td = tf;
// 			if (tf < td2)
// 				td2 = tf;
// 			td3 += tf;
// 		}
// 		td3 /= (double)m->m_laSingleMolIndex.GetSize();
// 		mprintf(WHITE,"  - %s:",m->m_sName);
// 		for (z2=strlen(m->m_sName);z2<=ti;z2++)
// 			mprintf(" ");
// 		mprintf(" Min. %12G, Max. %12G, Avg. %12G Debye.\n",td2,td,td3);
// 		if (m->m_laSingleMolIndex.GetSize() == 1)
// 		{
// 			for (z2=strlen(m->m_sName);z2<=ti;z2++)
// 				mprintf(" ");
// 			mprintf("    Cartesian dipole vector is ( %8G | %8G | %8G ).\n\n",sm->m_vDipole[0],sm->m_vDipole[1],sm->m_vDipole[2]);
// 		}
// 	}
// 
// 	mprintf("\n");
// 	if (g_bAdvanced2)
// 	{
// 		if (AskYesNo("    Write out cartesian dipole vectors of molecules in each step (y/n)? [no] ",false))
// 		{
// 			g_bDumpDipoleVector = true;
// 			mprintf("\n");
// 			ti = 0;
// 			for (z=0;z<g_oaMolecules.GetSize();z++)
// 			{
// 				if (AskYesNo("      Consider %s molecules for writing out dipole vectors (y/n)? [yes] ",true,((CMolecule*)g_oaMolecules[z])->m_sName))
// 				{
// 					try { tia = new CxIntArray("ParseDipole():tia"); } catch(...) { tia = NULL; }
// 					if (tia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
// 					
// 					g_oaDumpDipoleVector.Add(tia);
// 					AskString("        Which representants to consider (e.g. 1,3-6,8)? [all] ",buf,"*");
// 					if (buf[0] == '*')
// 					{
// 						for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
// 							tia->Add(z2+1);
// 					} else ParseIntList(buf,tia);
// 					for (z2=0;z2<tia->GetSize();z2++)
// 						tia->GetAt(z2)--;
// 					ti += tia->GetSize();
// 					mprintf("\n");
// 				} else g_oaDumpDipoleVector.Add(NULL);
// 			}
// 			g_iDumpDipoleSMCount = ti;
// 			g_bDumpDipoleAbs = AskYesNo("    Write out absolute value as 4th column per molecule (y/n)? [yes] ",true);
// 			g_bDumpDipoleXYZ = AskYesNo("    Write out XYZ trajectory to visualize dipole vectors (y/n)? [no] ",false);
// 
// 			if (g_bDumpDipoleXYZ)
// 			{
// 				g_fDumpDipoleScale = AskFloat("    Please enter dipole scale in trajectory (unit: Angstrom/Debye): [1.0] ",1.0f);
// 				g_iDumpDipoleXYZAtoms = 0;
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					if (g_oaDumpDipoleVector[z] == NULL)
// 						continue;
// 					g_iDumpDipoleXYZAtoms += ((CMolecule*)g_oaMolecules[z])->m_iAtomGesNoVirt * ((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();
// 				}
// 				g_iDumpDipoleXYZAtoms += g_iDumpDipoleSMCount*2;
// 			}
// 		} else g_bDumpDipoleVector = false;
// 	}
// 
// 	mprintf(WHITE,"\n<<< End of Dipole definition <<<\n\n");
// 	
// 	g_bDipoleDefined = true;
}

void parseMagneticDipole() {
	static bool magneticDipoleDefined = false;
	if (magneticDipoleDefined)
		return;
	
	g_bUseVelocities = true;
	
	if (g_bTegri && (g_pTetraPak == NULL)) {
		try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
		if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { g_pTetraPak = new CTetraPak(); } catch(...) { g_pTetraPak = NULL; }
		if (g_pTetraPak == NULL) NewException((double)sizeof(CTetraPak),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pTetraPak->Parse();
	}
		
	mprintf(WHITE, "\n>>> Magnetic Moment Definition >>>\n\n");
	
	mprintf("    There are the following possibilities to provide magnetic moments:\n");
	mprintf("    (1) Calculate magnetic moments from Wannier centers\n        (These need to be atoms in the trajectory and have to be sorted, use \"swan\" in the main menu)\n");
	mprintf("    (2) Read magnetic moments from an external file\n        (This file is expected to contain one line per timestep\n         with the three vector components separated by spaces)\n");
	mprintf("    (3) Read magnetic moments from the trajectory file\n        (These need to be specified in the comment line)\n");
	mprintf("    (4) Calculate magnetic moments from fixed atomic partial charges\n");
	mprintf("    (5) Calculate magnetic moments from fluctuating atomic partial charges\n        (These need to be in the fourth column of the trajectory)\n");
	mprintf("    (6) Calculate magnetic moments from Voronoi partial charges\n        (This needs electron density data)\n");
	mprintf("    (7) Calculate Voronoi magnetic moments\n        (This needs electron density data)\n");
	mprintf("    (8) Load a magnetic moment restart file\n");
	mprintf("\n");
	
	mprintf("    The following magnetic moment modes can be set for all molecules at once: 1, 5, 6, 7, 8\n");
	
	while (true) {
		int globalMode = AskRangeInteger("    Which magnetic moment mode to set for all molecules at once? [none] ", 0, 8, 0);
		
		if (globalMode == 0) {
			int i;
			for (i = 0; i < g_oaMolecules.GetSize(); i++)
				((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 0;
			break;
		} else if (globalMode == 1) {
			if (setupWannier()) {
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 1;
				break;
			}
			continue;
		} else if (globalMode == 2) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 3) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 4) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 5) {
			if (g_bXYZ4thCol || g_bLAMMPSCharge) {
				setupFluctuatingCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 5;
				break;
			}
			eprintf("    The trajectory does not have a fourth column.\n");
			continue;
		} else if (globalMode == 6) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				parseCoreCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 6;
				break;
			}
			eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			continue;
		} else if (globalMode == 7) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				g_bVoroIntegrateDipoleMoment = true;
				g_bVoroIntegrateTotalCurrent = true;
				g_bVoroIntegrateMagneticMoment = true;
				g_bCubeTimeDev = true;
				parseCoreCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 7;
				break;
			}
			eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			continue;
		} else if (globalMode == 8) {
			setupMagneticDipoleRestartFile();
			int i;
			for (i = 0; i < g_oaMolecules.GetSize(); i++)
				((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 8;
			break;
		}
		eprintf("    This is impossible.\n");
	}
	
	int i;
	for (i = 0; i < g_oaMolecules.GetSize(); i++) {
		unsigned char rty;
		ParseAtom("#2", i, ((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleCenterIndex);
	}
		
	while (true) {
		mprintf("\n    The following magnetic moment modes are set up:\n\n");
		int longest = 0;
		int i;
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			int length = strlen(((CMolecule *)g_oaMolecules[i])->m_sName);
			if (length > longest)
				longest = length;
		}
		CxString buf;
		buf.sprintf("      %%%ds - %%d", longest);
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			mprintf(buf, ((CMolecule *)g_oaMolecules[i])->m_sName, ((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode);
			if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 0)
				mprintf(" (none)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 1)
				mprintf(" (Wannier centers)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 2)
				mprintf(" (external file)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 3)
				mprintf(" (trajectory comment line)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 4)
				mprintf(" (fixed atomic partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 5)
				mprintf(" (fluctuating atomic partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 6)
				mprintf(" (Voronoi partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 7)
				mprintf(" (Voronoi magnetic moments)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 8)
				mprintf(" (restart file)\n");
			else
				mprintf("\n");
		}
		mprintf("\n");
		
		if (!AskYesNo("    Change magnetic moment mode for a molecule (y/n)? [no] ", false))
			break;
		
		int mol;
		if (g_oaMolecules.GetSize() > 1) {
			buf.sprintf("    Change dipole mode for which molecule (");
			for (i = 0; i < g_oaMolecules.GetSize(); i++) {
				CxString buf2;
				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i + 1);
				buf.strcat(buf2);
				if(i < g_oaMolecules.GetSize() - 1) {
					buf.strcat(", ");
				}
			}
			buf.strcat(")? ");
			
			mol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
		} else {
			mol = 0;
		}
		
		int mode = AskRangeInteger_ND("    Which magnetic moment mode to set for molecule %s? ", 1, 8, ((CMolecule *)g_oaMolecules[mol])->m_sName);
		
		if (mode == 0) {
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 0;
		} else if (mode == 1) {
			if (setupWannier()) {
				((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 1;
			}
		} else if (mode == 2) {
			FILE *dipoleFile = NULL;
			while (true) {
				if(dipoleFile != NULL) {
					fclose(dipoleFile);
					dipoleFile = NULL;
				}
				AskString_ND("    Name of file for magnetic moments: ", &buf);
				mprintf("\n    Trying to read first line...\n\n");
				dipoleFile = fopen(buf, "r");
				if(dipoleFile == NULL) {
					eprintf("    Could not open \"%s\".\n\n", (const char*)buf);
					continue;
				}
				if (buf.fgets(1024, dipoleFile) == NULL) {
					eprintf("    Could not read first line.\n\n");
					continue;
				}
				float dipole[3];
				if(sscanf((const char *)buf, "%f %f %f", &dipole[0], &dipole[1], &dipole[2]) < 3) {
					eprintf("    Could not find three floating point numbers in first line.\n\n");
					continue;
				}
				mprintf("    Found the following magnetic moment in the first line:\n");
				mprintf("      ( %8G | %8G | %8G ) Bohr magneton.\n", dipole[0], dipole[1], dipole[2]);
				rewind(dipoleFile);
				break;
			}
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 2;
			((CMolecule *)g_oaMolecules[mol])->m_pMagneticDipoleFile = dipoleFile;
			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
				eprintf("\n    Warning: There are several molecules %s. The same magnetic moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
		} else if (mode == 3) {
			CxString comment;
			comment.strcpy(g_TimeStep.m_pComment);
			mprintf("    The comment line is \"%s\".\n", (const char *)comment);
			buf.sprintf("");
			const char delim[] = "\t ";
			char *tok = strtok(comment.GetWritePointer(), delim);
			int num = 0;
			while(tok != NULL) {
				float f;
				if(sscanf(tok, "%f", &f) == 1) {
					CxString buf2;
					num++;
					if(num > 1) {
						buf2.sprintf(", %s (%d)", tok, num);
					} else {
						buf2.sprintf("%s (%d)", tok, num);
					}
					buf.strcat(buf2);
				}
				tok = strtok(NULL, delim);
			}
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleCommentIndex[0] = AskRangeInteger("    Choose the X magnetic moment component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleCommentIndex[1] = AskRangeInteger("    Choose the Y magnetic moment component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleCommentIndex[2] = AskRangeInteger("    Choose the Z magnetic moment component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 3;
			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
		} else if (mode == 4) {
			CMolecule *m = (CMolecule*)g_oaMolecules[mol];
			CSingleMolecule *sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			m->m_bChargesAssigned = true;
			int z2, z3, ti;
			CMolAtom *ma;
			CxFloatArray *tfa;
			bool b;
			float tf;
			double td;
			CxString buf2;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				mprintf("\n");
				ma = NULL;
				td = 1.0e50;
				
				try { tfa = new CxFloatArray("ParseDipole():tfa"); } catch(...) { tfa = NULL; }
				if (tfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				m->m_oaCharges.Add(tfa);
				_cnm:
				buf.sprintf("");
				buf2.sprintf("");
				b = false;
				ti = 0;
				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
				{
					if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
						continue;
					
					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
					
					if (ma->m_fAtomCode > td)
					{
						continue;
					}
					if (b)
					{
						if (ma->m_fAtomCode < td)
						{
							continue;
						}
						ti++;
						buf2.sprintf(", %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						buf.strcat(buf2);
						
					} else
					{
						if (ma->m_fAtomCode < td)
						{
							td = ma->m_fAtomCode;
						}
						
						buf.sprintf("    Enter partial atomic charge on %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						b = true;
						ti++;
					}
				}
				if (b)
				{
					buf.strcat(": [0.0] ");
					tf = AskFloat(buf,0);
					for (z3=0;z3<ti;z3++)
						tfa->Add(tf);
					td -= 1.0;
					goto _cnm;
				}
			}
			
			tf = 0;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				for (z3=0;z3<((CxFloatArray*)m->m_oaCharges[z2])->GetSize();z3++)
				{
					mprintf("    %s%d: %7.4f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3));
					tf += ((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3);
				}
			}
			m->m_fCharge = tf;
			mprintf(WHITE,"\n    Molecule %s: Total charge is %.4f.\n\n",m->m_sName,m->m_fCharge);
			
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 4;
		} else if (mode == 5) {
			if (g_bXYZ4thCol || g_bLAMMPSCharge) {
				setupFluctuatingCharges();
				((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 5;
			}
			eprintf("    The trajectory does not have a fourth column.\n");
		} else if (mode == 6) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				parseCoreCharges();
				((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 6;
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			}
		} else if (mode == 7) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				g_bVoroIntegrateDipoleMoment = true;
				g_bVoroIntegrateTotalCurrent = true;
				g_bVoroIntegrateMagneticMoment = true;
				g_bCubeTimeDev = true;
				parseCoreCharges();
				((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 7;
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			}
		} else if (mode == 8) {
			setupMagneticDipoleRestartFile();
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 8;
		} else {
			eprintf("    This is impossible.\n");
		}
	}
	mprintf("\n");
	
	if (g_bCubeTimeDev) {
		g_fBackgroundDensity = AskFloat("    Background density to improve PDE solver convergence in atomic units [1e-6] ", 1.0e-6f);
		g_fPDEConvThresh = AskFloat("    Relative convergence threshold for PDE solver [1e-2] ", 0.01f);
		g_iPDEMaxIter = AskInteger("    Maximum number of iterations in PDE solver [20] ", 20);
		mprintf("\n");
	}
	
	mprintf(WHITE, "<<< End of Magnetic Moment Definition <<<\n\n");
	
	magneticDipoleDefined = true;
}

void DipolGrimme(const char *s)
{
	FILE *a;
	char buf[256], *p, *q;
	const char *sep = " ,;\t";
	int i, ti, z, z3;
	CxFloatArray *fa, *ptfa2, *ptfa3;
	float fx, fy, fz, tf;
	CFFT tfft;
	CAutoCorrelation *ac;
	CReorDyn *m_pRDyn;

	mprintf("\n    This analysis reads dipole vectors from a text file\n");
	mprintf("    and calculates the dipole vector ACF / IR spectrum.\n\n");
	mprintf("    The text file needs to contain the three cartesian components\n");
	mprintf("    of the dipole vector (three numbers per line); one line equals one time step.\n\n");

	if (s == NULL)
	{
		eprintf("Please specify dipole text file as first command line parameter.\n\n");
		return;
	}

	mprintf("    If you are not sure what to enter, use the default values (just press return).\n\n");

	g_fTimestepLength = AskFloat("    Enter the length of one trajectory time step in fs: [0.5] ",0.5f);
	mprintf("\n");

	try { m_pRDyn = new CReorDyn(); } catch(...) { m_pRDyn = NULL; }
	if (m_pRDyn == NULL) NewException((double)sizeof(CReorDyn),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
//	m_pRDyn->m_sName = "dipole";

	try { m_pRDyn->m_pRDyn = new CDF(); } catch(...) { m_pRDyn->m_pRDyn = NULL; }
	if (m_pRDyn->m_pRDyn == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pRDyn->m_pRDyn->m_bLeft = true;
	mprintf(WHITE,"\n*** Vector Reorientation Dynamics ***\n\n");

	m_pRDyn->m_iDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the vector ACF (in time steps): [4096] ",4096);

	mprintf("\n    This corresponds to a spectral resolution of %.4f cm^-1.\n",33356.41/g_fTimestepLength/2.0/m_pRDyn->m_iDepth);

	ti = CalcFFTSize(m_pRDyn->m_iDepth,false);
	if (m_pRDyn->m_iDepth != ti)
	{
		mprintf(WHITE,"\n    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n",tfft.NextFastSize(m_pRDyn->m_iDepth),m_pRDyn->m_iDepth);
		m_pRDyn->m_iDepth = tfft.NextFastSize(m_pRDyn->m_iDepth);
	}

	try { m_pRDyn->m_pACF = new CACF(); } catch(...) { m_pRDyn->m_pACF = NULL; }
	if (m_pRDyn->m_pACF == NULL) NewException((double)sizeof(CACF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pRDyn->m_pACF->m_iSize = m_pRDyn->m_iDepth;
	m_pRDyn->m_pACF->m_bSpectrum = true;

	m_pRDyn->m_pACF->m_bWindowFunction = AskYesNo("    Apply window function (Cos^2) to Autocorrelation function (y/n)? [yes] ",true);

	tf = 33356.41 / g_fTimestepLength / 2.0;
	mprintf("\n    A time step length of %.1f fs allows a spectral range up to %.1f cm^-1.\n\n",g_fTimestepLength,tf);
	m_pRDyn->m_pACF->m_fSpecWaveNumber = AskRangeFloat("    Calculate spectrum up to which wave number (cm^-1)? [%.1f cm^-1] ",0,tf,(tf<5000.0)?tf:5000.0,(tf<5000.0)?tf:5000.0);
	m_pRDyn->m_pACF->m_iMirror = AskUnsignedInteger("    No mirroring (0), short-time enhancing (1) or long-time enhancing (2)? [1] ",1);
	m_pRDyn->m_pACF->m_iZeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ",m_pRDyn->m_iDepth*3,m_pRDyn->m_iDepth*3);
	m_pRDyn->m_pACF->m_iZeroPadding0 = m_pRDyn->m_pACF->m_iZeroPadding;

	ti = CalcFFTSize(m_pRDyn->m_iDepth+m_pRDyn->m_pACF->m_iZeroPadding,false);
	if (m_pRDyn->m_iDepth+m_pRDyn->m_pACF->m_iZeroPadding != ti)
	{
		mprintf(WHITE,"\n    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n",ti,ti-m_pRDyn->m_iDepth);
		m_pRDyn->m_pACF->m_iZeroPadding = ti-m_pRDyn->m_iDepth;
	}

	mprintf("    Zero padding increases the spectral resolution to %.4f cm^-1.\n\n",33356.41/g_fTimestepLength/2.0/(m_pRDyn->m_iDepth+m_pRDyn->m_pACF->m_iZeroPadding));

	m_pRDyn->m_pACF->m_bACF_DB = AskYesNo("    Convert intensity axis of spectrum to decibel (y/n)? [no] ",false);
	m_pRDyn->m_pACF->Create();

	m_pRDyn->m_pRDyn->m_fMinVal = 0;
	m_pRDyn->m_pRDyn->m_fMaxVal = m_pRDyn->m_iDepth * g_fTimestepLength / 1000.0f;
	m_pRDyn->m_pRDyn->m_iResolution = m_pRDyn->m_iDepth;
	m_pRDyn->m_pRDyn->SetLabelX("Tau [ps]");
	m_pRDyn->m_pRDyn->SetLabelY("Vector autocorrelation");
	m_pRDyn->m_pRDyn->Create();

	mprintf("    Trying to open dipole file \"%s\"...\n",s);

	a = fopen(s,"rt");

	if (a == NULL)
	{
		eprintf("Error: Could not open file for reading.\n\n");
		return;
	}

	mprintf("\nReading dipole data:\n");

	try { fa = new CxFloatArray("DipolGrimme():fa"); } catch(...) { fa = NULL; }
	if (fa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	i = 0;
	while (!feof(a))
	{
		fgets(buf,256,a);
		if (feof(a))
		{
			mprintf("\n\nEnd of dipole file reached.");
			break;
		}
		buf[strlen(buf)-1] = 0;

		p = buf;

		while (strchr(sep,*p) != NULL)
			p++;
		q = p;
		while ((strchr(sep,*p) == NULL) && (*p != 0))
			p++;
		*p = 0;
		p++;
		fx = atof(q);

		while (strchr(sep,*p) != NULL)
			p++;
		q = p;
		while ((strchr(sep,*p) == NULL) && (*p != 0))
			p++;
		*p = 0;
		p++;
		fy = atof(q);

		while (strchr(sep,*p) != NULL)
			p++;
		q = p;
		while ((strchr(sep,*p) == NULL) && (*p != 0))
			p++;
		*p = 0;
		fz = atof(q);

//		mprintf("\n( %f | %f | %f)",fx,fy,fz);

		fa->Add(fx);
		fa->Add(fy);
		fa->Add(fz);
	
		if ((i % 50) == 0)
			mprintf("\n%6d ",i);
		mprintf(".");

		i++;
	}
	mprintf("\n\n%d time steps read.\n\n",i);

	fclose(a);

	mprintf("    Autocorrelating cached vectors via FFT...\n");

	try { ptfa2 = new CxFloatArray("DipolGrimme():ptfa2"); } catch(...) { ptfa2 = NULL; }
	if (ptfa2 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	ptfa2->SetSize(i);

	try { ptfa3 = new CxFloatArray("DipolGrimme():ptfa3"); } catch(...) { ptfa3 = NULL; }
	if (ptfa3 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	ptfa3->SetSize(i);

	try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
	if (ac == NULL) NewException((double)sizeof(CAutoCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	ac->Init(i,m_pRDyn->m_iDepth,true);

	/* X */
	for (z3=0;z3<(int)i;z3++)
		(*ptfa2)[z3] = (*fa)[z3*3];

	ac->AutoCorrelate(ptfa2,ptfa3);

	for (z3=0;z3<(int)m_pRDyn->m_iDepth;z3++)
		m_pRDyn->m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3]);

	/* Y */
	for (z3=0;z3<(int)i;z3++)
		(*ptfa2)[z3] = (*fa)[z3*3+1];

	ac->AutoCorrelate(ptfa2,ptfa3);

	for (z3=0;z3<(int)m_pRDyn->m_iDepth;z3++)
		m_pRDyn->m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3]);

	/* Z */
	for (z3=0;z3<(int)i;z3++)
		(*ptfa2)[z3] = (*fa)[z3*3+2];

	ac->AutoCorrelate(ptfa2,ptfa3);

	for (z3=0;z3<(int)m_pRDyn->m_iDepth;z3++)
	{
		m_pRDyn->m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3]);
		m_pRDyn->m_pRDyn->m_fBinEntries += (double)(i-z3) - 3.0;
	}

	delete ac;
	delete ptfa2;
	delete ptfa3;

	mprintf("    Finished.\n\n");
	mprintf("    %.0f bin entries.\n",m_pRDyn->m_pRDyn->m_fBinEntries);
	mprintf("    Starting value is %f.\n",m_pRDyn->m_pRDyn->m_pBin[0]);
	m_pRDyn->m_pRDyn->MultiplyBin(1.0/m_pRDyn->m_pRDyn->m_pBin[0]);
	sprintf(buf,"rdyn_dipole.csv");
	mprintf("    Saving result as %s ...\n",buf);
	m_pRDyn->m_pRDyn->Write("",buf,"",false);

	mprintf("\n    Creating reorientation spectrum:\n");

	for (z=0;z<m_pRDyn->m_iDepth;z++)
 		m_pRDyn->m_pACF->m_pData[z] = m_pRDyn->m_pRDyn->m_pBin[z];

	if (m_pRDyn->m_pACF->m_iMirror != 0)
	{
		mprintf("    Mirroring ACF...\n");
		m_pRDyn->m_pACF->Mirror(m_pRDyn->m_pACF->m_iMirror);
		sprintf(buf,"rdyn_dipole.mirrored.csv");
		mprintf("    Saving mirrored ACF as %s ...\n",buf);
		m_pRDyn->m_pACF->WriteACF("",buf,"");
	}

	if (m_pRDyn->m_pACF->m_bWindowFunction)
	{
		mprintf("    Applying window function to ACF...\n");
		m_pRDyn->m_pACF->Window();
		sprintf(buf,"rdyn_dipole.windowed.csv");
		mprintf("    Saving windowed ACF as %s ...\n",buf);
		m_pRDyn->m_pACF->WriteACF("",buf,"");
	}

	mprintf("    Performing Fourier transformation...\n");

	try { g_pFFT = new CFFT(); } catch(...) { g_pFFT = NULL; }
	if (g_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	g_pFFT->PrepareFFT_C2C(m_pRDyn->m_pACF->m_iSize+m_pRDyn->m_pACF->m_iZeroPadding);

	m_pRDyn->m_pACF->Transform(g_pFFT);
	delete g_pFFT;
	m_pRDyn->m_pACF->m_pSpectrum->SetMaxRWL(1E7f/299.792f/g_fTimestepLength);
	if (m_pRDyn->m_pACF->m_bACF_DB)
	{
		mprintf("    Normalizing spectrum to decibel...\n");
		m_pRDyn->m_pACF->m_pSpectrum->MakeDB();
	}

	sprintf(buf,"spectrum_dipole.csv");
	mprintf("    Saving spectrum as %s ...\n",buf);
	m_pRDyn->m_pACF->m_pSpectrum->Write("",buf,"");

	mprintf("\n\n*** Finished ***\n\n");

	WriteCredits();
}


void InitGlobalVars()
{
	g_bReadLAMMPSCharges = false;
	g_bLAMMPSCharge = false;
	g_iVoroPrintLevel = 0;
	g_bSDFVoro = false;
	g_bDoubleBox = false;
	g_bWriteOrtho = false;
	g_bBoxNonOrtho = false;
	g_bFoundNonOrtho = false;
	g_bWriteInputOrder = false;
	g_bShowCredits = false;
	g_bPairMSD = false;
	g_bStreamInput = false;
	g_sConfFile = NULL;
	g_bShowConf = false;
	g_bWriteConf = false;
	g_bCheckWrite = true;
	g_bSaxonize = false;
	g_bUnknownElements = false;
	g_iStartTime = 0;
	g_bSMode = false;
	g_bNoColor = false;
	g_iColorIntensity = 1;
	g_bNPT = true; // Nur anfaenglich
//	g_sNPTFile[0] = 0;
	g_sNPTFile.sprintf("");
	g_fNPTFile = NULL;
	g_fBoxX = 0;
	g_fBoxY = 0;
	g_fBoxZ = 0;
	g_bVerbose = false;
	g_bNeedMoleculeWrap = false;
	g_bDipoleDefined = false;
	g_bDipolGrimme = false;
	g_fTimestepLength = 0;
	g_sInputTraj = NULL;
	g_sInputFile = NULL;
	g_bKeepOriginalCoords = false;
	g_bAbortAnalysis = true;
	g_bMiddleAvg = true;
	g_bVFHisto = false;
	g_bUseVelocities = false;
	g_bUseForces = false;
	g_bWriteAtomwise = false;
	g_fBondFactor = 1.15f;
	g_fVelPercentage = 0.95f;
	g_fForcePercentage = 0.95f;
	g_fPos = NULL;
	g_fVel = NULL;
	g_fForce = NULL;
	g_iRefSystemDim = 0;
	g_iFixMol = -1;
	g_fInputFile = NULL;
	g_fInput = NULL;
	g_fSaveCondFile = NULL;
	g_bSaveCondSnapshot = false;
	g_bScanVelocities = false;
	g_bSilentProgress = false;
	g_bCreateRevSDF = false;
	g_bDeriv = false;
	g_iStrideDetect = -1;
	g_bGlobalIR = false;
	g_pGlobalIR = NULL;
	g_bLMFitSilent = false;
	g_iLMMaxIter = 200;
	g_bXYZ4thCol = false;
	g_bXYZComment6Numbers = false;
	g_bReadChargesFrom4thXYZ = false;
	g_bEnvWriteDetailedInfo = false;
	g_bEnvSortNb = false;
	g_bEnvDisableSortNb = false;
	g_bSDFMap = false;
	g_iSDFMapSmoothGrade = 0;
	g_bVoroSilent = false;
	g_bVoronoiMoments = true;
	g_iWannierAtomType = -1;
}


CCrossCorrelation::CCrossCorrelation()                                                                                                                                                                                   
{                                                                                                                                                                                                                        
	m_iInput = 0;                                                                                                                                                                                                     
	m_iDepth = 0;                                                                                                                                                                                                     
	m_pFFT = NULL;                                                                                                                                                                                                    
	m_pFFT2 = NULL;                                                                                                                                                                                                 
}                                                                                                                                                                                                                          
                                                                                                                                                                                                                         
                                                                                                                                                                                                                         
CCrossCorrelation::~CCrossCorrelation()                                                                                                                                                                                  
{                                                                                                                                                                                                                        
}                                                                                                                                                                                                                        
                                                                                                                                                                                                                         
                                                                                                                                                                                                                         
void CCrossCorrelation::Init(int input, int depth, bool fft)                                                                                                                                                             
{                                                                                                                                                                                                                        
	m_iInput = input;                                                                                                                                                                                                 
	m_iDepth = depth;                                                                                                                                                                                                 
	if (fft)                                                                                                                                                                                                          
	{                                                                                                                                                                                                                 
		m_bFFT = true;                                                                                                                                                                                            
		m_iFFTSize = CalcFFTSize(input,true);                                                                                                                                                                     
                                                                                                                                                                                                                         
		try { m_pFFT = new CFFT(); } catch(...) { m_pFFT = NULL; }                                                                                                                                                
		if (m_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);                                                                                                             
		                                                                                                                                                                                                          
		m_pFFT->PrepareFFT_C2C(2*m_iFFTSize);                                                                                                                                                                     
                                                                                                                                                                                                                         
		try { m_pFFT2 = new CFFT(); } catch(...) { m_pFFT2 = NULL; }                                                                                                                                              
		if (m_pFFT2 == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);                                                                                                            
		                                                                                                                                                                                                          
		m_pFFT2->PrepareFFT_C2C(2*m_iFFTSize);                                                                                                                                                                    
                                                                                                                                                                                                                         
		try { m_pFFTback = new CFFT(); } catch(...) { m_pFFTback = NULL; }                                                                                                                                        
		if (m_pFFTback == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);                                                                                                         
		                                                                                                                                                                                                          
		m_pFFTback->PrepareInverseFFT_C2C(2*m_iFFTSize);                                                                                                                                                          
	} else                                                                                                                                                                                                            
	{                                                                                                                                                                                                                 
		m_bFFT = false;                                                                                                                                                                                           
	}                                                                                                                                                                                                                 
}                                                                                                                                                                                                                        
                                                                                                                                                                                                                         
                                                                                                                                                                                                                         
void CCrossCorrelation::CrossCorrelate(CxFloatArray *inp1, CxFloatArray *inp2, CxFloatArray *outp)                                                                                                                       
{                                                                                                                                                                                                                        
	int z, z2;                                                                                                                                                                                                        
	double tf;                                                                                                                                                                                                        
                                                                                                                                                                                                                         
	outp->SetSize(m_iDepth);                                                                                                                                                                                          
                                                                                                                                                                                                                         
	if (m_bFFT)                                                                                                                                                                                                       
	{                                                                                                                                                                                                                 
		for (z=0;z<m_iInput;z++)                                                                                                                                                                                  
		{                                                                                                                                                                                                         
			m_pFFT->m_pInput[z*2] = (*inp1)[z];                                                                                                                                                               
			m_pFFT->m_pInput[z*2+1] = 0;                                                                                                                                                                      
		}                                                                                                                                                                                                         
		for (z=m_iInput;z<2*m_iFFTSize;z++)                                                                                                                                                                       
		{                                                                                                                                                                                                         
			m_pFFT->m_pInput[z*2] = 0;                                                                                                                                                                        
			m_pFFT->m_pInput[z*2+1] = 0;                                                                                                                                                                      
		}                                                                                                                                                                                                         
		m_pFFT->DoFFT();                                                                                                                                                                                          
                                                                                                                                                                                                                         
		for (z=0;z<m_iInput;z++)                                                                                                                                                                                  
		{                                                                                                                                                                                                         
			m_pFFT2->m_pInput[z*2] = (*inp2)[z];                                                                                                                                                              
			m_pFFT2->m_pInput[z*2+1] = 0;                                                                                                                                                                     
		}                                                                                                                                                                                                         
		for (z=m_iInput;z<2*m_iFFTSize;z++)                                                                                                                                                                       
		{                                                                                                                                                                                                         
			m_pFFT2->m_pInput[z*2] = 0;                                                                                                                                                                       
			m_pFFT2->m_pInput[z*2+1] = 0;                                                                                                                                                                     
		}                                                                                                                                                                                                         
		m_pFFT2->DoFFT();                                                                                                                                                                                         
                                                                                                                                                                                                                         
		for (z=0;z<m_iFFTSize*2;z++)                                                                                                                                                                              
		{                                                                                                                                                                                                         
			// a1*a2 + b1*b2                                                                                                                                                                                  
			m_pFFTback->m_pInput[z*2]   = (float)((m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2]   + m_pFFT->m_pOutput[z*2+1]*m_pFFT2->m_pOutput[z*2+1]));                                                   
			// a2*b1 - a1*b2                                                                                                                                                                                  
			m_pFFTback->m_pInput[z*2+1] = (float)((-m_pFFT2->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2+1] + m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2+1]));                                                    
		}                                                                                                                                                                                                         
		m_pFFTback->DoFFT();                                                                                                                                                                                      
                                                                                                                                                                                                                         
		for (z=0;z<m_iDepth;z++)                                                                                                                                                                                  
			(*outp)[z] = (float)((double)m_pFFTback->m_pOutput[2*z] / m_iFFTSize / 2.0f / ((double)m_iInput - z));                                                                                            
	} else                                                                                                                                                                                                            
	{                                                                                                                                                                                                                 
		for (z=0;z<m_iDepth;z++) // Tau                                                                                                                                                                           
		{                                                                                                                                                                                                         
			tf = 0;                                                                                                                                                                                           
			for (z2=0;z2<m_iInput-z;z2++)                                                                                                                                                                     
				tf += (double)(*inp1)[z2] * (*inp2)[z2+z];                                                                                                                                                
			(*outp)[z] = (float)(tf / (double)(m_iInput-z));                                                                                                                                                  
		}                                                                                                                                                                                                         
	}                                                                                                                                                                                                                 
}                                                                                                                                                                                                                        
                                                                                                                                                                                                                        

void ParseVoronoiRadii()
{
	int i, z, z2, z3, z4;
	float *f, tf;
	CMolecule *m;

	if (g_faVoronoiRadii.GetSize() != 0)
		return; // Already parsed

	mprintf("\n    When performing the Voronoi decomposition, radii may be assigned to the atoms (\"radical Voronoi tesselation\").\n\n");

	i = AskRangeInteger("    Do not assign radii (0), use covalent radii (1), use Van-der-Waals radii (2), or specify other radii (3)? [0] ",0,3,0);

	if (i == 0)
	{
		g_faVoronoiRadii.SetSize(g_iGesAtomCount);

		for (z=0;z<g_iGesAtomCount;z++)
			g_faVoronoiRadii[z] = 0.5f;
	}

	if (i == 1)
	{
		mprintf("\n    Using the following atom radii:\n");
		for (z=0;z<g_oaAtoms.GetSize();z++)
		{
			if (z == g_iVirtAtomType)
				continue;

			mprintf("      %2s: %5.1f pm.\n",((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius);
			
			if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius < 1.0f)
				mprintf(RED, "Warning: The radius of %s is very small.\n", ((CAtom*)g_oaAtoms[z])->m_sName);
		}

		g_faVoronoiRadii.SetSize(g_iGesAtomCount);

		for (z=0;z<g_iGesAtomCount;z++)
			g_faVoronoiRadii[z] = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fRadius;
	}

	if (i == 2)
	{
		mprintf("\n    Using the following atom radii:\n");
		for (z=0;z<g_oaAtoms.GetSize();z++)
		{
			if (z == g_iVirtAtomType)
				continue;

			mprintf("      - %2s: %5.1f pm.\n",((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z])->m_pElement->m_fVdWRadius);
			
			if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius < 1.0f)
				mprintf(RED, "Warning: The radius of %s is very small.\n", ((CAtom*)g_oaAtoms[z])->m_sName);
		}

		g_faVoronoiRadii.SetSize(g_iGesAtomCount);

		for (z=0;z<g_iGesAtomCount;z++)
			g_faVoronoiRadii[z] = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fVdWRadius;
	}

	if (i == 3)
	{
		mprintf("    Enter 0 as radius to exclude atoms from the Voronoi decomposition.\n\n");

		g_faVoronoiRadii.SetSize(g_iGesAtomCount);

		if (AskYesNo("    Assign atom radii per element type (y) or per atom (n)? [yes] ",true))
		{
			mprintf("\n");
			f = new float[g_oaAtoms.GetSize()];

			for (z=0;z<g_oaAtoms.GetSize();z++)
			{
				if (z == g_iVirtAtomType)
					continue;
		
				f[z] = AskFloat_ND("      Which radius to use for element type %s (in pm)? ",((CAtom*)g_oaAtoms[z])->m_sName);
			}

			for (z=0;z<g_iGesAtomCount;z++)
				g_faVoronoiRadii[z] = f[g_waAtomRealElement[z]];

			delete[] f;
		} else
		{
			mprintf("\n");
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				m = (CMolecule*)g_oaMolecules[z];
				mprintf("    * Molecule %s\n\n",m->m_sName);
				for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
				{
					if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
						continue;
					for (z3=0;z3<m->m_waAtomCount[z2];z3++)
					{
						tf = AskFloat_ND("      Radius for %2s%-2d (in pm): ",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1);
						for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
							g_faVoronoiRadii[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]])->m_oaAtomOffset[z2])->GetAt(z3)] = tf;
					}
				}
				mprintf("\n");
			}
		}
	}

	if (i == 0)
		mprintf("\n    Not using Voronoi radii.\n\n");
	else mprintf("\n    Voronoi radii defined.\n\n");
}


void DumpNonOrthoCellData()
{
	double tf;
	CxVector3 veca, vecb, vecc;

	veca = CxVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
	vecb = CxVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
	vecc = CxVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

	mprintf(WHITE,"    *** Non-orthorhombic cell geometry ***\n");
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Cell vector A: ( %10.4f | %10.4f | %10.4f ), Length %10.3f pm.\n",g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2),veca.GetLength());
	mprintf(WHITE,"    *"); mprintf("   Cell vector B: ( %10.4f | %10.4f | %10.4f ), Length %10.3f pm.\n",g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2),vecb.GetLength());
	mprintf(WHITE,"    *"); mprintf("   Cell vector C: ( %10.4f | %10.4f | %10.4f ), Length %10.3f pm.\n",g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2),vecc.GetLength());
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Cell angle Alpha (between B and C): %8.3f degree.\n",g_fBoxAngleA);
	mprintf(WHITE,"    *"); mprintf("   Cell angle Beta  (between A and C): %8.3f degree.\n",g_fBoxAngleB);
	mprintf(WHITE,"    *"); mprintf("   Cell angle Gamma (between A and B): %8.3f degree.\n",g_fBoxAngleC);
	mprintf(WHITE,"    *\n");

	tf = DotP(CrossP(veca,vecb),vecc)/1000000.0;
	mprintf(WHITE,"    *"); mprintf("   Cell volume:  %.3f Angstrom^3 = %.6f nm^3.\n",tf,tf/1000.0);
	g_fBoxVolume = tf;
	mprintf(WHITE,"    *"); mprintf("   Cell density: %.6f g/cm^3.\n",pow(GuessBoxSize(),3)/tf/1000000.0);
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Orthogonal bounding box of the cell:\n");
	mprintf(WHITE,"    *"); mprintf("     dX = %10.4f pm,  dY = %10.4f pm,  dZ = %10.4f pm.\n",g_fBoxX,g_fBoxY,g_fBoxZ);
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Minimal periodic diameter (smallest distance between atom and its periodic image):\n");
	mprintf(WHITE,"    *"); mprintf("     Min. Diam. = %10.4f pm.\n",g_fBoxMinDiam);
	mprintf(WHITE,"    *"); mprintf("     ( dA = %10.4f pm,  dB = %10.4f pm,  dC = %10.4f pm )\n",g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Determinant of the forward transformation matrix is  %12G.\n",g_mBoxFromOrtho.Det());
	mprintf(WHITE,"    *"); mprintf("   Determinant of the backward transformation matrix is %12G.\n",g_mBoxToOrtho.Det());
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *** End of non-orthorhombic cell geometry ***\n\n");
}


void ExtractXYZCellGeometry(const char *s)
{
	int z;
	char buf[64];
	const char *p, *q;
	double data[6];
	CxVector3 veca, vecb, vecc;

	z = 0;
	p = s;
	while (*p != 0)
	{
//		mprintf("# (A) z=%d \"%s\".\n",z,p);
		while (strchr("0123456789+-.Ee",*p) == NULL)
			p++;
		q = p;
//		mprintf("# (B) z=%d \"%s\".\n",z,p);
		while ((strchr("0123456789+-.Ee",*q) != NULL) && (*q != 0))
			q++;
		if (*(q-1) == 0)
			q--;
//		mprintf("# (C) z=%d q=\"%s\" p=\"%s\".\n",z,q,p);
		if (q > p)
		{
//			mprintf("# q>p, q-p = %d\n",q-p);
//			mprintf("# p=\"%s\".\n",p);
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
//			mprintf("# buf = \"%s\".\n",buf);
			if (atof(buf) != 0)
			{
				data[z] = atof(buf);
				z++;
			}
		}// else mprintf("# q <= p\n");
//		mprintf("# z=%d.\n",z);
		if (*q == 0)
			break;
//		mprintf("# X.\n");
		if (z == 6)
			break;
//		mprintf("# Y.\n");
		p = q;
	}
	if (z < 6)
	{
		eprintf("\nError: ExtractXYZCellGeometry(): Comment line of XYZ trajectory does not contain 6 numbers.\n");
		eprintf("  (Comment line is \"%s\").   ",s);
		return;
	}

	g_fBoxX = data[0] * 100.0; // Angstrom --> pm
	g_fBoxY = data[1] * 100.0;
	g_fBoxZ = data[2] * 100.0;
	g_fBoxAngleA = data[3];
	g_fBoxAngleB = data[4];
	g_fBoxAngleC = data[5];
//	mprintf("    Computing vectors...\n");

	g_mBoxFromOrtho(0,0) = g_fBoxX;
	g_mBoxFromOrtho(0,1) = 0;
	g_mBoxFromOrtho(0,2) = 0;

	g_mBoxFromOrtho(1,0) = g_fBoxY*cos(g_fBoxAngleC*Pi/180.0);
	g_mBoxFromOrtho(1,1) = g_fBoxY*sin(g_fBoxAngleC*Pi/180.0);
	g_mBoxFromOrtho(1,2) = 0;

	g_mBoxFromOrtho(2,0) = g_fBoxZ*cos(g_fBoxAngleB*Pi/180.0);
	g_mBoxFromOrtho(2,1) = (-(g_mBoxFromOrtho(1,0)*g_fBoxZ*cos(g_fBoxAngleB*Pi/180.0)) + g_fBoxY*g_fBoxZ*cos(g_fBoxAngleA*Pi/180.0))/g_mBoxFromOrtho(1,1);
	g_mBoxFromOrtho(2,2) = sqrtf(-((pow(g_fBoxZ,2)*(pow(g_mBoxFromOrtho(1,1),2)*(-1 + pow(cos(g_fBoxAngleB*Pi/180.0),2)) + pow(g_mBoxFromOrtho(1,0)*cos(g_fBoxAngleB*Pi/180.0) - g_fBoxY*cos(g_fBoxAngleA*Pi/180.0),2)))/pow(g_mBoxFromOrtho(1,1),2)));

	if (g_bDoubleBox)
	{
		g_mBoxFromOrtho(0,0) *= g_iDoubleBoxX;
		g_mBoxFromOrtho(1,0) *= g_iDoubleBoxX;
		g_mBoxFromOrtho(2,0) *= g_iDoubleBoxX;

		g_mBoxFromOrtho(0,1) *= g_iDoubleBoxY;
		g_mBoxFromOrtho(1,1) *= g_iDoubleBoxY;
		g_mBoxFromOrtho(2,1) *= g_iDoubleBoxY;

		g_mBoxFromOrtho(0,2) *= g_iDoubleBoxZ;
		g_mBoxFromOrtho(1,2) *= g_iDoubleBoxZ;
		g_mBoxFromOrtho(2,2) *= g_iDoubleBoxZ;
	}

//	mprintf("    Recalculating angles from vectors...\n");

	veca = CxVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
	vecb = CxVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
	vecc = CxVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

	g_fBoxAngleA = acosf(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
	g_fBoxAngleB = acosf(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
	g_fBoxAngleC = acosf(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;

//	mprintf("    Computing inverse of transformation matrix...\n\n");
	g_mBoxToOrtho = CxMatrix3(g_mBoxFromOrtho);
	g_mBoxToOrtho.Invert();

	// Orthogonal bounding box
	g_fBoxX = fabsf(g_mBoxFromOrtho(0,0)) + fabsf(g_mBoxFromOrtho(1,0)) + fabsf(g_mBoxFromOrtho(2,0));
	g_fBoxY = fabsf(g_mBoxFromOrtho(0,1)) + fabsf(g_mBoxFromOrtho(1,1)) + fabsf(g_mBoxFromOrtho(2,1));
	g_fBoxZ = fabsf(g_mBoxFromOrtho(0,2)) + fabsf(g_mBoxFromOrtho(1,2)) + fabsf(g_mBoxFromOrtho(2,2));

	veca = CxVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
	vecb = CxVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
	vecc = CxVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

	// Minimal diameters
	g_fBoxMinDiamA = fabsf(DotP(veca,Normalize(CrossP(vecb,vecc))));
	g_fBoxMinDiamB = fabsf(DotP(vecb,Normalize(CrossP(veca,vecc))));
	g_fBoxMinDiamC = fabsf(DotP(vecc,Normalize(CrossP(veca,vecb))));
	g_fBoxMinDiam = MIN3(g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);
}




