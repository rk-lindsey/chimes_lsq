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


#include "tetrapak.h"

#include "maintools.h"
#include "tools.h"
#include "globalvar.h"
#include "conversion.h"

static char *g_cubeMemFileNameFixed = NULL;
static int g_cubeMemFileIndex = 1;
static int g_cubeMemFileStride = 1;


CTetraPak::CTetraPak()
{
// 	m_atomChargeFile = NULL;
// 	m_molChargeFile = NULL;
}


CTetraPak::~CTetraPak()
{
// 	if (m_atomChargeFile != NULL)
// 		fclose(m_atomChargeFile);
// 	if (m_molChargeFile != NULL)
// 		fclose(m_molChargeFile);
	if (g_cubeMemFileNameFixed != NULL)
		delete[] g_cubeMemFileNameFixed;
}


void CTetraPak::Parse()
{
	CTimeStep *t;
	int z, i, rx, ry, rz;
	double fs;
	bool sanityall;
//	char buf[256];
	CxString buf;
	FILE *a;

	mprintf(WHITE,"\n>>> Voronoi Integration Functions >>>\n\n");

	mprintf("    Initializing Voronoi tesselation...\n");
	g_pVoroWrapper->Init();
	mprintf("\n");

	mprintf("*** Voro: Box density is %f particles / Angstrom^3.\n",g_pVoroWrapper->m_fBoxDens);
	mprintf("*** Voro: Using %d x %d x %d blocks.\n",g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ);

	try { t = new CTimeStep(); } catch(...) { t = NULL; }
	if (t == NULL) NewException((double)sizeof(CTimeStep),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	t->CopyFrom(&g_TimeStep);

	t->FoldAtomsPositive();
	g_pVoroWrapper->Dump("voro.txt",t);
	mprintf("\n");
	mprintf("    Voro++: Using cell memory for %d particles.\n\n",g_iVoroMemory);

//	Removed: Did not work.
//	g_bVoroIntEquitable = AskYesNo("    Use equitable binning for Voronoi integration (y/n)? [no] ",false);

	if (g_TimeStep.m_pVolumetricData != NULL)
	{
		m_bInterpolation = AskYesNo("    Use on-the-fly interpolation of volumetric data (y/n)? [no] ",false);

		if (m_bInterpolation)
		{
			mprintf("\n    Interpolation currently only works for charges, not for moments.\n\n");
//			m_bInterpolationLinear = AskYesNo("    Use linear interpolation (y) or constant interpolation (n)? [yes] ",true);
			m_iInterpolationFactor = AskUnsignedInteger("    Enter interpolation factor: [2] ",2);

			m_iInterpolationFactorCube = m_iInterpolationFactor * m_iInterpolationFactor * m_iInterpolationFactor;

			rx = g_TimeStep.m_pVolumetricData->m_iRes[0];
			ry = g_TimeStep.m_pVolumetricData->m_iRes[1];
			rz = g_TimeStep.m_pVolumetricData->m_iRes[2];
			mprintf("\n    Original grid of %d x %d x %d will be interpolated to %d x %d x %d.\n",rx,ry,rz,rx*m_iInterpolationFactor,ry*m_iInterpolationFactor,rz*m_iInterpolationFactor);
		} else {
			m_iInterpolationFactor = 1;
			m_iInterpolationFactorCube = 1;
		}
	} else {
		m_bInterpolation = false;
		m_iInterpolationFactor = 1;
		m_iInterpolationFactorCube = 1;
	}

	mprintf("\n    Performing sanity check...\n\n");

	if (g_TimeStep.m_pVolumetricData != NULL)
	{
		rx = g_TimeStep.m_pVolumetricData->m_iRes[0];
		ry = g_TimeStep.m_pVolumetricData->m_iRes[1];
		rz = g_TimeStep.m_pVolumetricData->m_iRes[2];
		mprintf("    Input trajectory contains volumetric data. Using grid of %d x %d x %d.\n",rx,ry,rz);
	} else
	{
		rx = 100;
		ry = 100;
		rz = 100;
		mprintf("    Input trajectory does not contain volumetric data. Using grid of 100 x 100 x 100.\n");
	}

	m_p3DF = new C3DF<VORI_FLOAT>();
	m_p3DF->m_iRes[0] = rx;
	m_p3DF->m_iRes[1] = ry;
	m_p3DF->m_iRes[2] = rz;
	m_p3DF->m_fMinVal[0] = 0;
	m_p3DF->m_fMinVal[1] = 0;
	m_p3DF->m_fMinVal[2] = 0;
// 	m_p3DF->m_fMaxVal[0] = g_fBoxX;
// 	m_p3DF->m_fMaxVal[1] = g_fBoxY;
// 	m_p3DF->m_fMaxVal[2] = g_fBoxZ;
// 	m_p3DF->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
// 	m_p3DF->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
// 	m_p3DF->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
// 	m_p3DF->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * rx;
// 	m_p3DF->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * ry;
// 	m_p3DF->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * rz;
	m_p3DF->m_fMaxVal[0] = g_TimeStep.m_pVolumetricData->m_fMaxVal[0];
	m_p3DF->m_fMaxVal[1] = g_TimeStep.m_pVolumetricData->m_fMaxVal[1];
	m_p3DF->m_fMaxVal[2] = g_TimeStep.m_pVolumetricData->m_fMaxVal[2];
	m_p3DF->Create();
	
	try { m_hitCount = new int[rx * ry * rz * m_iInterpolationFactorCube]; } catch (...) { m_hitCount = NULL; }
	if (m_hitCount == NULL) NewException((double)sizeof(int) * rx * ry * rz * m_iInterpolationFactorCube, __FILE__, __LINE__, __PRETTY_FUNCTION__);
	memset(m_hitCount, 0, rx * ry * rz * m_iInterpolationFactorCube * sizeof(int));

	if (g_bAdvanced2)
	{
		mprintf("\n");
		sanityall = AskYesNo("    Perform sanity check for first time step (n) or for all steps in trajectory (y)? [no] ",false);
	} else
		sanityall = false;

	if (sanityall)
	{
		mprintf("\n");
		g_fPos = fopen(g_sInputTraj,"rt");
		if (g_fPos == NULL)
		{
			eprintf("Error. Could not open \"%s\".\n",g_sInputTraj);
			abort();
		}

		i = 0;
		while (true)
		{
			mprintf(YELLOW,"Step %5d:  ",i+1);
			i++;

			if (!g_TimeStep.ReadTimestep(g_fPos,false))
				break;

			g_TimeStep.CalcCenters();

			t->CopyFrom(&g_TimeStep);

			t->FoldAtomsPositive();

			for (z=0;z<rx*ry*rz;z++)
				m_p3DF->m_pBin[z] = 0;
			
			for (z = 0; z < rx * ry * rz * m_iInterpolationFactorCube; z++)
				m_hitCount[z] = 0;

			if (!BuildVoronoi(t,true,true,true))
			{
				eprintf("\n    Sanity check failed.\n\n");
				if (!AskYesNo("    This should not have happened. Continue anyway (y/n)? [yes] ",true))
					abort();
				mprintf("\n");
			}
		}
		mprintf("\n    All done. Leaving.\n");
		exit(0);
	} else
	{
		if (!BuildVoronoi(t,true,true))
		{
			eprintf("\n    Sanity check failed.\n\n");
			if (!AskYesNo("    This should not have happened. Continue anyway (y/n)? [yes] ",true))
				abort();
			mprintf("\n");
		} else
		{
			mprintf("\n    Sanity check done.\n\n");
		}
	}

//	m_p3DF->WritePLT("test.plt","","",true);

	delete m_p3DF;
	
	delete[] m_hitCount;

	if (g_bAdvanced2) {
		if (AskYesNo("    Compute Voronoi charges from a single electron density cube file (y/n)? [no] ",false))
		{
			mprintf("\n");
_n2:
			AskString_ND("    Enter file name of the cube file: ",&buf);

			a = fopen(buf,"rt");
			if (a == NULL)
			{
				eprintf("    Could not open \"%s\" for reading.\n",(const char*)buf);
				goto _n2;
			}

			mprintf("\n");

			m_p3DF = new C3DF<VORI_FLOAT>();
			m_p3DF->ReadCube(a,true,true);

			fclose(a);

			mprintf("\n    Calculating Voronoi Partial Charges...\n");

			m_faCharge.SetSize(g_iGesAtomCount);
			for (z=0;z<g_iGesAtomCount;z++)
				m_faCharge[z] = 0;

			BuildVoronoi(t,true,false);

			mprintf("    Done.\n\n    Results:\n");

			fs = 0;
			for (z=0;z<g_iGesAtomCount;z++)
			{
// 			m_faCharge[z] *= (g_fBoxX/100.0/0.529177249/m_p3DF->m_iRes[0])*(g_fBoxY/100.0/0.529177249/m_p3DF->m_iRes[1])*(g_fBoxZ/100.0/0.529177249/m_p3DF->m_iRes[2]);
				m_faCharge[z] *= g_fCubeXStep*g_fCubeYStep*g_fCubeZStep;
				
				if (mystricmp(((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,"H") == 0)
					m_faCharge[z] += 1.0;
				else if (mystricmp(((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,"C") == 0)
					m_faCharge[z] += 4.0;
				else if (mystricmp(((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,"N") == 0)
					m_faCharge[z] += 5.0;
				else if (mystricmp(((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,"O") == 0)
					m_faCharge[z] += 6.0;
				else
					eprintf("Unknown atom type.\n");

				mprintf("      Atom %4d (%-2s): Charge %10.6f.\n",z+1,((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,m_faCharge[z]);
				fs += m_faCharge[z];
			}

			mprintf("\n    Sum of charges: %10.6f.\n",fs);
		}
		mprintf("\n");
	}
	delete t;

	if (g_iTrajFormat == 5)
		m_bVoronoiCharges = AskYesNo("    Compute Voronoi charges in each step from electron density in cube trajectory (y/n)? [yes] ",true);
	else
		m_bVoronoiCharges = AskYesNo("    Compute Voronoi charges in each step from a stream of electron density cube files (y/n)? [no] ",false);

	if (m_bVoronoiCharges)
	{
		mprintf("\n");
		if (g_iTrajFormat == 5)
		{
			g_bCubeStream = AskYesNo("    Use streaming mode (y/n)? [no] ", false);
			if (g_bCubeStream) {
				try { g_fCubeMemFile = new CxMemFile(); } catch (...) { g_fCubeMemFile = NULL; }
				if (g_fCubeMemFile == NULL) NewException((double)sizeof(g_fCubeMemFile), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				
				CxString buf;
				AskString_ND("    Enter file base name: ", &buf);
				if (g_cubeMemFileNameFixed != NULL)
					delete[] g_cubeMemFileNameFixed;
				try { g_cubeMemFileNameFixed = new char[strlen(buf) + 1]; } catch (...) { g_cubeMemFileNameFixed = NULL; }
				if (g_cubeMemFileNameFixed == NULL) NewException((double)sizeof(char) * (strlen(buf) + 1), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				strcpy(g_cubeMemFileNameFixed, buf);
				g_cubeMemFileIndex = AskUnsignedInteger("    First index to process: [1] ", 1);
				g_cubeMemFileStride = AskUnsignedInteger("    Stride of cube files: [1] ", 1);
				g_iCubeMemFileSteps = AskUnsignedInteger_ND("    Number of cube files to read: ");
				mprintf("\n    Expecting files \"%s%d.cube\", \"%s%d.cube\", ..., \"%s%d.cube\" to appear\n", g_cubeMemFileNameFixed, g_cubeMemFileIndex, g_cubeMemFileNameFixed, g_cubeMemFileIndex + g_cubeMemFileStride, g_cubeMemFileNameFixed, g_cubeMemFileIndex + (g_iCubeMemFileSteps - 1) * g_cubeMemFileStride);
				g_cubeMemFileIndex -= g_cubeMemFileStride;
				
				g_iCubeMemFileLines = (g_TimeStep.m_pVolumetricData->m_iRes[2] + g_iCubeZMismatch) / g_iCubeZStride / 6;
				if (g_TimeStep.m_pVolumetricData->m_iRes[2] % 6 > 0)
					g_iCubeMemFileLines++;
				g_iCubeMemFileLines *= (g_TimeStep.m_pVolumetricData->m_iRes[1] + g_iCubeYMismatch) / g_iCubeYStride * (g_TimeStep.m_pVolumetricData->m_iRes[0] + g_iCubeZMismatch) / g_iCubeZStride;
				g_iCubeMemFileLines += g_iGesAtomCount + 6;
			} else {
				mprintf("    Taking electron density from input cube trajectory.\n");
			}
		} else
		{
_nameagain:
			AskString_ND("    Enter name of the cube file for the volumetric electron density data: ",&buf);
			m_fCubePipe = fopen(buf,"rt");
			if (m_fCubePipe == NULL)
			{
				eprintf("    Could not open \"%s\" for reading.\n");
				goto _nameagain;
			}
		}
		if (g_bVoroIntegrateCharge)
			m_faCharge.SetSize(g_iGesAtomCount);
		if (g_bVoroIntegrateDipoleMoment)
			m_moments.SetSize(g_iGesAtomCount);
		if (g_bVoroIntegrateTotalCurrent)
			m_totalCurrent.SetSize(g_iGesAtomCount);
		if (g_bVoroIntegrateMagneticMoment)
			m_magneticMoments.SetSize(g_iGesAtomCount);
		m_p3DF = NULL;

		
		if (g_bAdvanced2)
			m_saveTotalIntegrals = AskYesNo("\n    Save total cube integrals in each step (y/n)? [no] ", false);
		else
			m_saveTotalIntegrals = false;
		
		if (m_saveTotalIntegrals) {
			m_totalIntegralFile = OpenFileWrite("vori_integrals.csv", false);
		}
	}

	if (g_bSDF)
	{
		mprintf("\n");
		g_bSDFVoro = AskYesNo("    Compute a SDF of Voronoi cells (y/n)? [yes] ",true);

		if (g_bSDFVoro)
		{
		}
	}
	
// 	m_saveAtomCharges = AskYesNo("\n    Save atomic Voronoi charges (y/n)? [no] ", false);
// 	if (m_saveAtomCharges) {
// 		m_atomChargeFile = OpenFileWrite("vori_charge_atoms.csv", false);
// 		fprintf(m_atomChargeFile, "#Time (fs);");
// 		int i, j, k, l;
// 		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
// 			CMolecule *m = (CMolecule *)g_oaMolecules[i];
// 			for (j = 0; j < m->m_laSingleMolIndex.GetSize(); j++) {
// 				CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[j]];
// 				for (k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
// 					if (m->m_baAtomIndex[k] == g_iVirtAtomType)
// 						continue;
// 					for (l = 0; l < ((CxIntArray *)sm->m_oaAtomOffset[k])->GetSize(); l++) {
// 						fprintf(m_atomChargeFile, " %s[%d]-%s%d;", m->m_sName, j + 1, ((CAtom *)g_oaAtoms[m->m_baAtomIndex[k]])->m_sName, l + 1);
// 					}
// 				}
// 			}
// 		}
// 		fprintf(m_atomChargeFile, "\n");
// 	}
// 	
// 	m_saveMolCharges = AskYesNo("    Save molecular Voronoi charges (y/n)? [no] ", false);
// 	if (m_saveMolCharges) {
// 		m_molChargeFile = OpenFileWrite("vori_charge_molecules.csv", false);
// 		fprintf(m_molChargeFile, "#Time (fs);");
// 		int i, j;
// 		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
// 			CMolecule *m = (CMolecule *)g_oaMolecules[i];
// 			for (j = 0; j < m->m_laSingleMolIndex.GetSize(); j++) {
// 				fprintf(m_molChargeFile, " %s[%d];", m->m_sName, j + 1);
// 			}
// 		}
// 		fprintf(m_molChargeFile, "\n");
// 	}

	mprintf(WHITE,"\n<<< End of Voronoi Integration Functions <<<\n\n");
	
// 	if (g_bVCD || g_bMagneticDipoleRestart) {
// 		g_bDipole = true;
// 		g_bUseVelocities = true;
// 		g_bCubeTimeDev = true;
// // 		g_fTimestepLength = 0.1f;
// 	}
}


void CTetraPak::ProcessStep(CTimeStep *ts, bool verbose)
{
	int z;
	double su;
	CTimeStep temp;
	
	double integralFactor;
	if (g_bBoxNonOrtho) {
		integralFactor = g_fCubeXVector[0] * g_fCubeYVector[1] * g_fCubeZVector[2] + g_fCubeXVector[1] * g_fCubeYVector[2] * g_fCubeZVector[0] + g_fCubeXVector[2] * g_fCubeYVector[0] * g_fCubeZVector[1] - g_fCubeXVector[0] * g_fCubeYVector[2] * g_fCubeZVector[1] - g_fCubeXVector[1] * g_fCubeYVector[0] * g_fCubeZVector[2] - g_fCubeXVector[2] * g_fCubeYVector[1] * g_fCubeZVector[0];
	} else {
		integralFactor = g_fCubeXStep * g_fCubeYStep * g_fCubeZStep;
	}

	if (m_bVoronoiCharges)
	{
		temp.CopyFrom(ts);
		temp.FoldAtomsPositive();
		
		if (g_bVoroIntegrateCharge) {
			m_faCharge.SetSize(g_iGesAtomCount);
			ts->m_faCharge.SetSize(g_iGesAtomCount);
		}
		if (g_bVoroIntegrateDipoleMoment) {
			m_moments.SetSize(g_iGesAtomCount);
			ts->m_dipoleMoments.SetSize(g_iGesAtomCount);
		}
		if (g_bVoroIntegrateTotalCurrent) {
			m_totalCurrent.SetSize(g_iGesAtomCount);
			ts->m_totalCurrents.SetSize(g_iGesAtomCount);
		}
		if (g_bVoroIntegrateMagneticMoment) {
			m_magneticMoments.SetSize(g_iGesAtomCount);
			ts->m_magneticDipoleMoments.SetSize(g_iGesAtomCount);
		}
		//		fseek(m_fCubePipe,0,SEEK_SET);

		if (g_iTrajFormat == 5)
		{
			m_p3DF = ts->m_pVolumetricData;
			
			if (m_saveTotalIntegrals || verbose)
			{
				su = 0;
				for (z=0;z<m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2];z++)
					su += m_p3DF->m_pBin[z];
//				double totCharge = -su * g_fCubeXStep*g_fCubeYStep*g_fCubeZStep;
				double totCharge = -su * integralFactor;
				if (g_bVoroIntegrateCharge) {
					int i;
					for (i = 0; i < g_iGesAtomCount; i++) {
						totCharge += ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge;
					}
				}
				double totDip[3] = { 0.0, 0.0, 0.0 };
				if (g_bVoroIntegrateDipoleMoment) {
					int i, j, k;
					for (i = 0; i < m_p3DF->m_iRes[0]; i++) {
//						double x = (double)i * g_fCubeXStep * LEN_AU2PM;
						for (j = 0; j < m_p3DF->m_iRes[1]; j++) {
//							double y = (double)j * g_fCubeYStep * LEN_AU2PM;
							for (k = 0; k < m_p3DF->m_iRes[2]; k++) {
								double x, y, z;
								if (g_bBoxNonOrtho) {
									x = ((double)i * g_fCubeXVector[0] + (double)j * g_fCubeYVector[0] + (double)k * g_fCubeZVector[0]) * LEN_AU2PM;
									y = ((double)i * g_fCubeXVector[1] + (double)j * g_fCubeYVector[1] + (double)k * g_fCubeZVector[1]) * LEN_AU2PM;
									z = ((double)i * g_fCubeXVector[2] + (double)j * g_fCubeYVector[2] + (double)k * g_fCubeZVector[2]) * LEN_AU2PM;
								} else {
									x = (double)i * g_fCubeXStep * LEN_AU2PM;
									y = (double)j * g_fCubeYStep * LEN_AU2PM;
									z = (double)k * g_fCubeZStep * LEN_AU2PM;
								}
//								double z = (double)k * g_fCubeZStep * LEN_AU2PM;
								totDip[0] -= x * m_p3DF->m_pBin[i + j * m_p3DF->m_iRes[0] + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0]];
								totDip[1] -= y * m_p3DF->m_pBin[i + j * m_p3DF->m_iRes[0] + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0]];
								totDip[2] -= z * m_p3DF->m_pBin[i + j * m_p3DF->m_iRes[0] + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0]];
							}
						}
					}
					for (i = 0; i < 3; i++)
//						totDip[i] *= g_fCubeXStep*g_fCubeYStep*g_fCubeZStep;
						totDip[i] *= integralFactor;
					for (i = 0; i < g_iGesAtomCount; i++) {
						for (j = 0; j < 3; j++)
							totDip[j] += ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge * temp.m_vaCoords[i][j];
					}
					for (i = 0; i < 3; i++)
						totDip[i] *= DIP_EPM2DEBYE;
				}
				double totCurrent[3] = { 0.0, 0.0, 0.0 };
				if (g_bVoroIntegrateTotalCurrent && ts->m_pCurrentDensity != NULL) {
					int i, j, k;
					for (i = 0; i < m_p3DF->m_iRes[0]; i++) {
						for (j = 0; j < m_p3DF->m_iRes[1]; j++) {
							for (k = 0; k < m_p3DF->m_iRes[2]; k++) {
								totCurrent[0] += ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3);
								totCurrent[1] += ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3 + 1);
								totCurrent[2] += ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3 + 2);
							}
						}
					}
					for (i = 0; i < 3; i++)
//						totCurrent[i] *= g_fCubeXStep*g_fCubeYStep*g_fCubeZStep * CURR_AUFS2MBPM;
						totCurrent[i] *= integralFactor * CURR_AUFS2MBPM;
					for (i = 0; i < g_iGesAtomCount; i++) {
						for (j = 0; j < 3; j++)
							totCurrent[j] += ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge * temp.m_vaVelocities[i][j] * CURR_EMS2MBPM;
					}
				}
				double totMagMom[3] = { 0.0, 0.0, 0.0 };
				if (g_bVoroIntegrateMagneticMoment && ts->m_pCurrentDensity != NULL) {
					int i, j, k;
					for (i = 0; i < m_p3DF->m_iRes[0]; i++) {
//						double x = (double)i * g_fCubeXStep * LEN_AU2PM;
						for (j = 0; j < m_p3DF->m_iRes[1]; j++) {
//							double y = (double)j * g_fCubeYStep * LEN_AU2PM;
							for (k = 0; k < m_p3DF->m_iRes[2]; k++) {
								double x, y, z;
								if (g_bBoxNonOrtho) {
										x = ((double)i * g_fCubeXVector[0] + (double)j * g_fCubeYVector[0] + (double)k * g_fCubeZVector[0]) * LEN_AU2PM;
										y = ((double)i * g_fCubeXVector[1] + (double)j * g_fCubeYVector[1] + (double)k * g_fCubeZVector[1]) * LEN_AU2PM;
										z = ((double)i * g_fCubeXVector[2] + (double)j * g_fCubeYVector[2] + (double)k * g_fCubeZVector[2]) * LEN_AU2PM;
								} else {
										x = (double)i * g_fCubeXStep * LEN_AU2PM;
										y = (double)j * g_fCubeYStep * LEN_AU2PM;
										z = (double)k * g_fCubeZStep * LEN_AU2PM;
								}
//								double z = (double)k * g_fCubeZStep * LEN_AU2PM;
								totMagMom[0] += y * ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3 + 2);
								totMagMom[0] -= z * ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3 + 1);
								totMagMom[1] += z * ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3);
								totMagMom[1] -= x * ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3 + 2);
								totMagMom[2] += x * ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3 + 1);
								totMagMom[2] -= y * ts->m_pCurrentDensity->GetAt(i * 3 + j * m_p3DF->m_iRes[0] * 3 + k * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[0] * 3);
							}
						}
					}
					for (i = 0; i < 3; i++)
//						totMagMom[i] *= 0.5 * g_fCubeXStep*g_fCubeYStep*g_fCubeZStep * MAG_AUPMFS2MB;
						totMagMom[i] *= 0.5 * integralFactor * MAG_AUPMFS2MB;
					for (i = 0; i < g_iGesAtomCount; i++) {
						totMagMom[0] += 0.5 * ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge * temp.m_vaCoords[i][1] * temp.m_vaVelocities[i][2] * MAG_EPMMS2MB;
						totMagMom[0] -= 0.5 * ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge * temp.m_vaCoords[i][2] * temp.m_vaVelocities[i][1] * MAG_EPMMS2MB;
						totMagMom[1] += 0.5 * ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge * temp.m_vaCoords[i][2] * temp.m_vaVelocities[i][0] * MAG_EPMMS2MB;
						totMagMom[1] -= 0.5 * ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge * temp.m_vaCoords[i][0] * temp.m_vaVelocities[i][2] * MAG_EPMMS2MB;
						totMagMom[2] += 0.5 * ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge * temp.m_vaCoords[i][0] * temp.m_vaVelocities[i][1] * MAG_EPMMS2MB;
						totMagMom[2] -= 0.5 * ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_fCharge * temp.m_vaCoords[i][1] * temp.m_vaVelocities[i][0] * MAG_EPMMS2MB;
					}
				}
				
				if (verbose) {
					mprintf("    Volumetric electron density data read from cube file:\n");
					mprintf("      Cube file has %d x %d x %d grid points.\n",m_p3DF->m_iRes[0],m_p3DF->m_iRes[1],m_p3DF->m_iRes[2]);
	// 				mprintf("      Sum is %.3f, equals %.6f electrons.\n",su,su*(g_fBoxX/100.0/0.529177249/m_p3DF->m_iRes[0])*(g_fBoxY/100.0/0.529177249/m_p3DF->m_iRes[1])*(g_fBoxZ/100.0/0.529177249/m_p3DF->m_iRes[2]));
//					mprintf("      Sum is %.3f, equals %.6f electrons.\n",su,su*g_fCubeXStep*g_fCubeYStep*g_fCubeZStep);
					mprintf("      Sum is %.3f, equals %.6f electrons.\n",su,su*integralFactor);
					if (g_bVoroIntegrateCharge)
						mprintf("      Total charge is %.6f.\n", totCharge);
					if (g_bVoroIntegrateDipoleMoment)
						mprintf("      Total dipole moment is ( %.6f | %.6f | %.6f ) Debye.\n", totDip[0], totDip[1], totDip[2]);
					if (g_bVoroIntegrateTotalCurrent && ts->m_pCurrentDensity != NULL)
						mprintf("      Total current is ( %.6f | %.6f | %.6f ) Bohr magneton * pm^-1.\n", totCurrent[0], totCurrent[1], totCurrent[2]);
					if (g_bVoroIntegrateMagneticMoment && ts->m_pCurrentDensity != NULL)
						mprintf("      Total magnetic moment is ( %.6f | %.6f | %.6f ) Bohr magneton.\n", totMagMom[0], totMagMom[1], totMagMom[2]);
				}
				if (m_saveTotalIntegrals) {
					fprintf(m_totalIntegralFile, "%lu", g_iSteps);
					if (g_bVoroIntegrateCharge)
						fprintf(m_totalIntegralFile, "; %.8G", totCharge);
					if (g_bVoroIntegrateDipoleMoment)
						fprintf(m_totalIntegralFile, "; %.8G; %.8G; %.8G", totDip[0], totDip[1], totDip[2]);
					if (g_bVoroIntegrateTotalCurrent)
						fprintf(m_totalIntegralFile, "; %.8G; %.8G; %.8G", totCurrent[0], totCurrent[1], totCurrent[2]);
					if (g_bVoroIntegrateMagneticMoment)
						fprintf(m_totalIntegralFile, "; %.8G; %.8G; %.8G", totMagMom[0], totMagMom[1], totMagMom[2]);
					fprintf(m_totalIntegralFile, "\n");
				}
			}
		} else
		{
			if (m_p3DF == NULL)
			{
				m_p3DF = new C3DF<VORI_FLOAT>();
				m_p3DF->ReadCube(m_fCubePipe,verbose,true);
			} else m_p3DF->ReadCube(m_fCubePipe,verbose,false);
		}

		for (z=0;z<g_iGesAtomCount;z++)
			m_faCharge[z] = 0;

		if (verbose)
			mprintf("\n    Integrating over electron density grid...\n");

		BuildVoronoi(&temp,verbose,false);
		
		for (z=0;z<g_iGesAtomCount;z++)
		{
// 			m_faCharge[z] *= (g_fBoxX/100.0/0.529177249/m_p3DF->m_iRes[0])*(g_fBoxY/100.0/0.529177249/m_p3DF->m_iRes[1])*(g_fBoxZ/100.0/0.529177249/m_p3DF->m_iRes[2]);
// 			m_faCharge[z] *= g_fCubeXStep*g_fCubeYStep*g_fCubeZStep;
			m_faCharge[z] *= integralFactor;
			
			m_faCharge[z] += ((CAtom *)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge;

			if (g_bVoroIntegrateCharge)
				ts->m_faCharge[z] = m_faCharge[z];

// 			ts->m_dipoleMoments[z] = m_moments[z] * 1000.0 * (g_fBoxX/100.0/0.529177249/m_p3DF->m_iRes[0])*(g_fBoxY/100.0/0.529177249/m_p3DF->m_iRes[1])*(g_fBoxZ/100.0/0.529177249/m_p3DF->m_iRes[2]);

			if (g_bVoroIntegrateDipoleMoment)
// 				ts->m_dipoleMoments[z] = m_moments[z] * 1000.0 * g_fCubeXStep*g_fCubeYStep*g_fCubeZStep; // Conversion to e*pm
				ts->m_dipoleMoments[z] = m_moments[z] * 1000.0 * integralFactor; // Conversion to e*pm
			if (g_bVoroIntegrateTotalCurrent)
// 				ts->m_totalCurrents[z] = m_totalCurrent[z] * g_fCubeXStep * g_fCubeYStep * g_fCubeZStep * CURR_AUFS2MBPM; // Conversion to Bohr magnetons*pm^-1
				ts->m_totalCurrents[z] = m_totalCurrent[z] * integralFactor * CURR_AUFS2MBPM; // Conversion to Bohr magnetons*pm^-1
			if (g_bVoroIntegrateMagneticMoment)
// 				ts->m_magneticDipoleMoments[z] = m_magneticMoments[z] * 0.5 * g_fCubeXStep * g_fCubeYStep * g_fCubeZStep * 1000.0 * MAG_AUPMFS2MB; // Conversion to Bohr magnetons
				ts->m_magneticDipoleMoments[z] = m_magneticMoments[z] * 0.5 * integralFactor * 1000.0 * MAG_AUPMFS2MB; // Conversion to Bohr magnetons
				
// 						mprintf("%s %f\n", ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName, m_faCharge[z]);

		}
		
// 		if (m_saveAtomCharges || m_saveMolCharges) {
// 			if (m_saveAtomCharges)
// 				fprintf(m_atomChargeFile, "%.2f;", g_iSteps * g_fTimestepLength * g_iStride);
// 			if (m_saveMolCharges)
// 				fprintf(m_molChargeFile, "%.2f;", g_iSteps * g_fTimestepLength * g_iStride);
// 			
// 			int i, j, k, l;
// 			for (i = 0; i < g_oaMolecules.GetSize(); i++) {
// 				CMolecule *m = (CMolecule *)g_oaMolecules[i];
// 				for (j = 0; j < m->m_laSingleMolIndex.GetSize(); j++) {
// 					CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[j]];
// 					float sum = 0.0f;
// 					for (k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
// 						if (m->m_baAtomIndex[k] == g_iVirtAtomType)
// 							continue;
// 						for (l = 0; l < ((CxIntArray *)sm->m_oaAtomOffset[k])->GetSize(); l++) {
// 							if (m_saveAtomCharges) {
// 								fprintf(m_atomChargeFile, " %.6f;", m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]);
// 							}
// 							sum += m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)];
// 						}
// 					}
// 					if (m_saveMolCharges) {
// 						fprintf(m_molChargeFile, " %.6f;", sum);
// 					}
// 				}
// 			}
// 			
// 			if (m_saveAtomCharges)
// 				fprintf(m_atomChargeFile, "\n");
// 			if (m_saveMolCharges)
// 				fprintf(m_molChargeFile, "\n");
// 		}
	}
}

bool CTetraPak::BuildVoronoi(CTimeStep *ts, bool verbose, bool sanity, bool sanityall)
{
	voronoicell_neighbor c;
	container_periodic_poly *con;
	int ijk, q, z, z2, fc, id, ti, cc;
	double *pp, tx, ty, tz;
	CxDVector3 vec1;
	vector<int> nb;
	vector<int> fv;
	vector<int> fo;
	vector<double> fare;
	CTetraFace *fa;
	double tfx, tfy, tfz, lastvert[2];
	double mi[2], ma[2];
	bool b;

	b = true;

//	mprintf("Container Size: %f | %f | %f.\n",g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0);

// 	try { con = new container_periodic_poly(g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0,g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ,g_iVoroMemory); } catch(...) { con = NULL; }
// 	mprintf(GREEN, "%.20g %.20g %.20g\n", m_p3DF->m_fMaxVal[0]/1000.0, m_p3DF->m_fMaxVal[1]/1000.0, m_p3DF->m_fMaxVal[2]/1000.0);
	if (g_bBoxNonOrtho) {
		if (g_mCubeCell(0, 1) != 0 || g_mCubeCell(0, 2) != 0 || g_mCubeCell(1, 2) != 0) {
			eprintf("This cell matrix cannot be handled by the Voro++ library.\nThe first cell vector has to be parallel to the x axis and the second cell vector has to be in the xy plane.\n");
			abort();
		}
		try { con = new container_periodic_poly(g_mCubeCell(0, 0) / 1000.0, g_mCubeCell(1, 0) / 1000.0, g_mCubeCell(1, 1) / 1000.0, g_mCubeCell(2, 0) / 1000.0, g_mCubeCell(2, 1) / 1000.0, g_mCubeCell(2, 2) / 1000.0, g_pVoroWrapper->m_iBlocksX, g_pVoroWrapper->m_iBlocksY, g_pVoroWrapper->m_iBlocksZ, g_iVoroMemory); } catch(...) { con = NULL; }
// 		mprintf(GREEN, "%.20g %.20g %.20g %.20g %.20g %.20g\n", g_mCubeCell(0, 0) / 1000.0, g_mCubeCell(1, 0) / 1000.0, g_mCubeCell(1, 1) / 1000.0, g_mCubeCell(2, 0) / 1000.0, g_mCubeCell(2, 1) / 1000.0, g_mCubeCell(2, 2) / 1000.0);
	} else {
		try { con = new container_periodic_poly(m_p3DF->m_fMaxVal[0]/1000.0,0,m_p3DF->m_fMaxVal[1]/1000.0,0,0,m_p3DF->m_fMaxVal[2]/1000.0,g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ,g_iVoroMemory); } catch(...) { con = NULL; }
// 		mprintf(GREEN, "%.20g %.20g %.20g\n", m_p3DF->m_fMaxVal[0]/1000.0, m_p3DF->m_fMaxVal[1]/1000.0, m_p3DF->m_fMaxVal[2]/1000.0);
		// 		mprintf(GREEN, "%.20g %.20g %.20g\n", g_fCubeXStep, g_fCubeYStep, g_fCubeZStep);
	}
	if (con == NULL) NewException((double)sizeof(container_periodic_poly),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<g_iGesAtomCount;z++) {
		con->put(z,ts->m_vaCoords[z][0]/1000.0,ts->m_vaCoords[z][1]/1000.0,ts->m_vaCoords[z][2]/1000.0,g_faVoronoiRadii[z]/1000.0);
// 		mprintf(GREEN, "%d %.20g %.20g %.20g %.20g\n", z+1, ts->m_vaCoords[z][0]/1000.0, ts->m_vaCoords[z][1]/1000.0, ts->m_vaCoords[z][2]/1000.0, g_faVoronoiRadii[z]/1000.0);
	}
	
	c_loop_all_periodic vl(*con);

	tfx = 0;
	tfy = 0;

	if (sanity)
	{
		if (sanityall)
			mprintf("  Sanity check:");
		else
			mprintf("\n    Voronoi decomposition and integration over grid:\n");
	}

	if (verbose)
		mprintf(WHITE,"      [");
	
	m_fRayCount = 0;
	m_fRayHisto[0] = 0;
	m_fRayHisto[1] = 0;
	m_fRayHisto[2] = 0;
	m_fRayHisto[3] = 0;
	cc = 0;
	if (vl.start()) 
	{
		do 
		{
			if (verbose)
				if (fmod(cc,g_iGesAtomCount/60.0) < 1.0)
					mprintf(WHITE,"#");

			if (con->compute_cell(c,vl))
			{
				ijk = vl.ijk;
				q = vl.q;
				pp = con->p[ijk]+con->ps*q;

				id = con->id[ijk][q];

// 				if (id == 103) {
// 					mprintf("\n****************************************************************************************************\n");
// 					mprintf("***   Cell %3d *************************************************************************************\n",id);
// 					mprintf("****************************************************************************************************\n");
// 					mprintf("Center: ( %.20G | %.20G | %.20G ).\n", pp[0], pp[1], pp[2]);
// 					CxDVector3 v(pp[0], pp[1], pp[2]);
// 					v *= 1000.0;
// 					v = g_mBoxToOrtho * v;
// 					mprintf("Center2: ( %.20G | %.20G | %.20G ).\n", v[0], v[1], v[2]);
// 				}

				m_iFaces = c.number_of_faces();
				c.face_vertices(fv);
				c.neighbors(nb);
				c.face_areas(fare);
				c.face_orders(fo);

// 				mprintf("Found %d faces.\n",m_iFaces);

				for (z=m_oaFaces.GetSize();z<m_iFaces;z++)
				{
					fa = new CTetraFace();
					m_oaFaces.Add(fa);
				}

				mi[0] = 1.0E20f;
				mi[1] = 1.0E20f;
				ma[0] = -1.0E20f;
				ma[1] = -1.0E20f;

				ti = 0;
				for (z=0;z<m_iFaces;z++)
				{
// 					if (id == 103) {
// 						mprintf("  Face %3d: Neighbor is %d, %d edges, area %.10f.\n",z+1,nb[z],fo[z],fare[z]);
// 						if (cc < 2)
// 							mprintf("****************** Face %d ***********************\n",z+1);
// 					}

					fc = fv[ti];

					fa = (CTetraFace*)m_oaFaces[z];

					fa->m_fMin[0] = 1.0E20f;
					fa->m_fMin[1] = 1.0E20f;
					fa->m_fMax[0] = -1.0E20f;
					fa->m_fMax[1] = -1.0E20f;

					tx = 0;
					ty = 0;
					tz = 0;
					for (z2=0;z2<fc;z2++)
					{
						tfy = pp[1]+c.pts[fv[ti+z2+1]*3+1]*0.5;
						tfz = pp[2]+c.pts[fv[ti+z2+1]*3+2]*0.5;

						tx += *pp+c.pts[fv[ti+z2+1]*3]*0.5;
						ty += tfy;
						tz += tfz;

// 						if (id == 103) {
// 							CxDVector3 v(*pp+c.pts[fv[ti+z2+1]*3]*0.5,pp[1]+c.pts[fv[ti+z2+1]*3+1]*0.5,pp[2]+c.pts[fv[ti+z2+1]*3+2]*0.5);
// 							v *= 1000.0;
// 							v = g_mBoxToOrtho * v;
// 							mprintf("  Vertex %d: ( %.20G | %.20G | %.20G ).\n",z2+1,v[0], v[1], v[2]);
// 							mprintf("  Vertex %d: ( %.20G | %.20G | %.20G ).\n",z2+1,*pp+c.pts[fv[ti+z2+1]*3]*0.5,pp[1]+c.pts[fv[ti+z2+1]*3+1]*0.5,pp[2]+c.pts[fv[ti+z2+1]*3+2]*0.5);
// 						}
						
						if (g_bBoxNonOrtho) {
							CxDVector3 coord((*pp+c.pts[fv[ti+z2+1]*3]*0.5) * 1000.0, tfy * 1000.0, tfz * 1000.0);
							coord = g_mBoxToOrtho * coord;
							if (fa->m_fMin[0] > coord[1])
								fa->m_fMin[0] = coord[1];
							if (fa->m_fMax[0] < coord[1])
								fa->m_fMax[0] = coord[1];
							if (fa->m_fMin[1] > coord[2])
								fa->m_fMin[1] = coord[2];
							if (fa->m_fMax[1] < coord[2])
								fa->m_fMax[1] = coord[2];
						} else {
							if (fa->m_fMin[0] > tfy)
								fa->m_fMin[0] = tfy;
							if (fa->m_fMax[0] < tfy)
								fa->m_fMax[0] = tfy;
							if (fa->m_fMin[1] > tfz)
								fa->m_fMin[1] = tfz;
							if (fa->m_fMax[1] < tfz)
								fa->m_fMax[1] = tfz;
						}
					}

// 					mprintf("Sum t: ( %20G | %20G | %20G ).\n",tx,ty,tz);

		/*			if ((tx != 0) && (fabs(tx) < 1.0E-14))
						tx = 0;
					if ((ty != 0) && (fabs(ty) < 1.0E-14))
						ty = 0;
					if ((tz != 0) && (fabs(tz) < 1.0E-14))
						tz = 0;*/

		//			mprintf("Chopped t: ( %20G | %20G | %20G ).\n",tx,ty,tz);

					tx /= fc;
					ty /= fc;
					tz /= fc;

// 					mprintf("Divided t: ( %20G | %20G | %20G ).\n",tx,ty,tz);

					if (fa->m_fMin[0] < mi[0])
						mi[0] = fa->m_fMin[0];
					if (fa->m_fMax[0] > ma[0])
						ma[0] = fa->m_fMax[0];
					if (fa->m_fMin[1] < mi[1])
						mi[1] = fa->m_fMin[1];
					if (fa->m_fMax[1] > ma[1])
						ma[1] = fa->m_fMax[1];

					fa->m_vCenter[0] = tx;
					fa->m_vCenter[1] = ty;
					fa->m_vCenter[2] = tz;

					fa->m_vSpan1[0] = *pp+c.pts[fv[ti+1]*3]*0.5;
					fa->m_vSpan1[1] = pp[1]+c.pts[fv[ti+1]*3+1]*0.5;
					fa->m_vSpan1[2] = pp[2]+c.pts[fv[ti+1]*3+2]*0.5;

					vec1[0] = *pp+c.pts[fv[ti+2]*3]*0.5;
					vec1[1] = pp[1]+c.pts[fv[ti+2]*3+1]*0.5;
					vec1[2] = pp[2]+c.pts[fv[ti+2]*3+2]*0.5;

// 					if (id == 103)
// 						mprintf("m_vSpan1: ( %20G | %20G | %20G ).\n",fa->m_vSpan1[0],fa->m_vSpan1[1],fa->m_vSpan1[2]);
// 					mprintf("vec1: ( %20G | %20G | %20G ).\n",vec1[0],vec1[1],vec1[2]);

					fa->m_vSpan1 -= fa->m_vCenter;
					vec1 -= fa->m_vCenter;

// 					if (id == 103)
// 						mprintf("m_vSpan1: ( %20G | %20G | %20G ).\n",fa->m_vSpan1[0],fa->m_vSpan1[1],fa->m_vSpan1[2]);
					
					fa->m_vSpan1.Normalize();
					fa->m_vSpan1.Chop(1.0E-14);
					
// 					if (id == 103)
// 						mprintf("m_vSpan1: ( %20G | %20G | %20G ).\n",fa->m_vSpan1[0],fa->m_vSpan1[1],fa->m_vSpan1[2]);
					
		/*			if ((fa->m_vSpan1[0] != 0) && (fabs(fa->m_vSpan1[0]) < 1.0E-14))
						fa->m_vSpan1[0] = 0;
					if ((fa->m_vSpan1[1] != 0) && (fabs(fa->m_vSpan1[1]) < 1.0E-14))
						fa->m_vSpan1[1] = 0;
					if ((fa->m_vSpan1[2] != 0) && (fabs(fa->m_vSpan1[2]) < 1.0E-14))
						fa->m_vSpan1[2] = 0;*/

// 					mprintf("m_vCenter: ( %20G | %20G | %20G ).\n",fa->m_vCenter[0],fa->m_vCenter[1],fa->m_vCenter[2]);

// 					mprintf("(*) m_vSpan1: ( %20G | %20G | %20G ).\n",fa->m_vSpan1[0],fa->m_vSpan1[1],fa->m_vSpan1[2]);
// 					mprintf("(*) vec1: ( %20G | %20G | %20G ).\n",vec1[0],vec1[1],vec1[2]);

					fa->m_vNormal = CrossP(fa->m_vSpan1,vec1);

// 					if (id == 103)
// 						mprintf("m_vNormal: ( %20G | %20G | %20G ).\n",fa->m_vNormal[0],fa->m_vNormal[1],fa->m_vNormal[2]);
					
					fa->m_vNormal.Normalize();
					fa->m_vNormal.Chop(1.0E-14);
					
// 					if (id == 103)
// 						mprintf("m_vNormal: ( %20G | %20G | %20G ).\n",fa->m_vNormal[0],fa->m_vNormal[1],fa->m_vNormal[2]);
					
// 					mprintf("Normal is ( %20G | %20G | %20G ).\n",fa->m_vNormal[0],fa->m_vNormal[1],fa->m_vNormal[2]);

		/*			if ((fa->m_vNormal[0] != 0) && (fabs(fa->m_vNormal[0]) < 1.0E-14))
						fa->m_vNormal[0] = 0;
					if ((fa->m_vNormal[1] != 0) && (fabs(fa->m_vNormal[1]) < 1.0E-14))
						fa->m_vNormal[1] = 0;
					if ((fa->m_vNormal[2] != 0) && (fabs(fa->m_vNormal[2]) < 1.0E-14))
						fa->m_vNormal[2] = 0;*/

		//			mprintf("Chopped normal is ( %20G | %20G | %20G ).\n",fa->m_vNormal[0],fa->m_vNormal[1],fa->m_vNormal[2]);

					fa->m_vSpan2 = CrossP(fa->m_vSpan1,fa->m_vNormal);

// 					if (id == 103)
// 						mprintf(GREEN, "span2: %.20g %.20g %.20g\n", fa->m_vSpan2[0], fa->m_vSpan2[1], fa->m_vSpan2[2]);
					
					fa->m_vSpan2.Normalize();
 					fa->m_vSpan2.Chop(1.0E-14);
					
// 					if (id == 103)
// 						mprintf(GREEN, "span2: %.20g %.20g %.20g\n", fa->m_vSpan2[0], fa->m_vSpan2[1], fa->m_vSpan2[2]);

		/*			if ((fa->m_vSpan2[0] != 0) && (fabs(fa->m_vSpan2[0]) < 1.0E-14))
						fa->m_vSpan2[0] = 0;
					if ((fa->m_vSpan2[1] != 0) && (fabs(fa->m_vSpan2[1]) < 1.0E-14))
						fa->m_vSpan2[1] = 0;
					if ((fa->m_vSpan2[2] != 0) && (fabs(fa->m_vSpan2[2]) < 1.0E-14))
						fa->m_vSpan2[2] = 0;*/

/*					if (fabs(DotP(CxDVector3(1.0,0,0),fa->m_vNormal)) < 1.0e-10)
					{
						mprintf("## Error: Normal vector in Y-Z plane: %.20G.\n",DotP(CxDVector3(1.0,0,0),fa->m_vNormal));
						mprintf("   Span1:  %.20G  %.20G  %.20G\n",fa->m_vSpan1[0],fa->m_vSpan1[1],fa->m_vSpan1[2]);
						mprintf("   vec1:   %.20G  %.20G  %.20G\n",vec1[0],vec1[1],vec1[2]);
						mprintf("   Normal: %.20G  %.20G  %.20G\n",fa->m_vNormal[0],fa->m_vNormal[1],fa->m_vNormal[2]);
						mprintf("   Span2:  %.20G  %.20G  %.20G\n",fa->m_vSpan2[0],fa->m_vSpan2[1],fa->m_vSpan2[2]);
					} */

// 					fa->m_vSpan1.Normalize();
// 					fa->m_vSpan2.Normalize();
// 					fa->m_vNormal.Normalize();
					
// 					if (id == 103) {
// 						mprintf(GREEN, "span1: %.20g %.20g %.20g\n", fa->m_vSpan1[0], fa->m_vSpan1[1], fa->m_vSpan1[2]);
// 						mprintf(GREEN, "span2: %.20g %.20g %.20g\n", fa->m_vSpan2[0], fa->m_vSpan2[1], fa->m_vSpan2[2]);
// 						mprintf(GREEN, "normal: %.20g %.20g %.20g\n", fa->m_vNormal[0], fa->m_vNormal[1], fa->m_vNormal[2]);
// 					}
					
					fa->m_vaEdges.RemoveAll_KeepSize();
					fa->m_faEdgesLength.RemoveAll_KeepSize();

					#define vec1x (fa->m_vSpan1[0])
					#define vec1y (fa->m_vSpan1[1])
					#define vec1z (fa->m_vSpan1[2])
					#define vec2x (fa->m_vSpan2[0])
					#define vec2y (fa->m_vSpan2[1])
					#define vec2z (fa->m_vSpan2[2])
					#define vecnx (fa->m_vNormal[0])
					#define vecny (fa->m_vNormal[1])
					#define vecnz (fa->m_vNormal[2])

//					tf = -vec1z*vec2y*vecnx + vec1y*vec2z*vecnx + vec1z*vec2x*vecny - vec1x*vec2z*vecny - vec1y*vec2x*vecnz + vec1x*vec2y*vecnz;

//					if (cc < 2)
//						mprintf("    tf=%10.6f\n",tf);


					for (z2=0;z2<fc+1;z2++)
					{
						tx = (*pp+c.pts[fv[ti+(z2%fc)+1]*3]*0.5)-fa->m_vCenter[0];
						ty = (pp[1]+c.pts[fv[ti+(z2%fc)+1]*3+1]*0.5)-fa->m_vCenter[1];
						tz = (pp[2]+c.pts[fv[ti+(z2%fc)+1]*3+2]*0.5)-fa->m_vCenter[2];

						lastvert[0] = tfx;
						lastvert[1] = tfy;

						tfx = - (-tz*vec2y*vecnx + ty*vec2z*vecnx + tz*vec2x*vecny - tx*vec2z*vecny - ty*vec2x*vecnz + tx*vec2y*vecnz) /*/ tf*/;
						tfy = - ( tz*vec1y*vecnx - ty*vec1z*vecnx - tz*vec1x*vecny + tx*vec1z*vecny + ty*vec1x*vecnz - tx*vec1y*vecnz) /*/ tf*/;
//						tfn = - (-tz*vec1y*vec2x + ty*vec1z*vec2x + tz*vec1x*vec2y - tx*vec1z*vec2y - ty*vec1x*vec2z + tx*vec1y*vec2z) /*/ tf*/;

//						if (cc < 2)
//							mprintf("  Vertex %d: a=%10.6f,  b=%10.6f,  c=%10.6f\n",z2+1,tfx,tfy,tfn);

						if (z2 > 0)
						{
							fa->m_vaEdges.Add(CxDVector3(tfy-lastvert[1],-tfx+lastvert[0],-(tfy-lastvert[1])*tfx + (tfx-lastvert[0])*tfy));
							fa->m_faEdgesLength.Add(sqrt((tfx-lastvert[0])*(tfx-lastvert[0]) + (tfy-lastvert[1])*(tfy-lastvert[1])));
// 							mprintf(GREEN, "%.10g %.10g %.10g %.10g %.10g\n", tfx, tfy, lastvert[0], lastvert[1], sqrt((tfx-lastvert[0])*(tfx-lastvert[0]) + (tfy-lastvert[1])*(tfy-lastvert[1])));
// 							mprintf(GREEN, "%.20g %.20g %.20g\n", tfy-lastvert[1],-tfx+lastvert[0],-(tfy-lastvert[1])*tfx + (tfx-lastvert[0])*tfy);

//							if (cc < 2)
//								mprintf("    da=%10.6f,  db=%10.6f,  cv=%1.06f\n",fa->m_vaEdges[fa->m_vaEdges.GetSize()-1][0],fa->m_vaEdges[fa->m_vaEdges.GetSize()-1][1],fa->m_vaEdges[fa->m_vaEdges.GetSize()-1][2]);
						}
					}


/*					if (cc < 2)
					{
						mprintf("\n");
						for (z2=0;z2<30;z2++)
						{
							tfy = fa->m_fMin[0] + (float)z2*(fa->m_fMax[0]-fa->m_fMin[0])/29.0f;
							for (z3=0;z3<60;z3++)
							{
								tfn = fa->m_fMin[1] + (float)z3*(fa->m_fMax[1]-fa->m_fMin[1])/59.0f;

								for (z4=0;z4<fc;z4++)
								{
									ty = (pp[1]+c.pts[fv[ti+z4+1]*3+1]*0.5-g_fBoxY/2000.0);
									tz = (pp[2]+c.pts[fv[ti+z4+1]*3+2]*0.5-g_fBoxZ/2000.0);

									if (((int)((fa->m_vCenter[1]-fa->m_fMin[0])/(fa->m_fMax[0]-fa->m_fMin[0])*29.0+0.5) == z2) &&
										((int)((fa->m_vCenter[2]-fa->m_fMin[1])/(fa->m_fMax[1]-fa->m_fMin[1])*59.0+0.5) == z3))
									{
										mprintf("*");
										goto _done;
									}

									if (((int)((ty-fa->m_fMin[0])/(fa->m_fMax[0]-fa->m_fMin[0])*29.0+0.5) == z2) &&
										((int)((tz-fa->m_fMin[1])/(fa->m_fMax[1]-fa->m_fMin[1])*59.0+0.5) == z3))
									{
										mprintf("+");
										goto _done;
									}
								}

								if (fa->XRay_Hit(tfy,tfn,tfx))
									mprintf("#");
								else mprintf(".");
_done:;
							}
							mprintf("\n");
						}
						mprintf("\n");

					}*/

					ti += fv[ti]+1;
				}

// 				if (magnetic) {
// 					m_totalCurrent[id] = integrateTotalCurrent(m_p3DF, ts->m_pCurrentDensity, mi, ma);
// 					m_magneticMoments[id] = integrateMagneticMoment(m_p3DF, ts->m_pCurrentDensity, mi, ma, pp);
// 				} else {
// 					if (sanity) {
// 						Integrate_Verbose(m_p3DF,mi,ma);
// 					} else {
//         					if (m_bInterpolation)
//         						m_faCharge[id] -= Integrate_Refine(m_p3DF,mi,ma);
//         					else
//         						m_faCharge[id] -= Integrate(m_p3DF,mi,ma);
// 						if (g_bVoronoiMoments)
// 							m_moments[id] = integrateMoment(m_p3DF, mi, ma, pp);
// 					}
// 				}
// 
				double charge;
				CxDVector3 dipoleMoment, totalCurrent, magneticMoment;
				if (sanity) {
// 					Integrate_Verbose(m_p3DF,mi,ma);
					integrateCell(m_p3DF, ts->m_pCurrentDensity, mi, ma, pp, &charge, &dipoleMoment, &totalCurrent, &magneticMoment, true/*, id*/);
				} else {
// 					if (m_bInterpolation) {
// 						m_faCharge[id] -= Integrate_Refine(m_p3DF,mi,ma);
// 					} else {
					integrateCell(m_p3DF, ts->m_pCurrentDensity, mi, ma, pp, &charge, &dipoleMoment, &totalCurrent, &magneticMoment);
					if (g_bVoroIntegrateCharge)
						m_faCharge[id] -= charge;
					if (g_bVoroIntegrateDipoleMoment)
						m_moments[id] = dipoleMoment;
					if (g_bVoroIntegrateTotalCurrent)
						m_totalCurrent[id] = totalCurrent;
					if (g_bVoroIntegrateMagneticMoment)
						m_magneticMoments[id] = magneticMoment;
// 					}
				}
				
/*
//				mprintf("\n");
				for (z2=0;z2<30;z2++)
				{
					tfy = mi[0] + (float)z2*(ma[0]-mi[0])/29.0f;
					for (z3=0;z3<60;z3++)
					{
						tfn = mi[1] + (float)z3*(ma[1]-mi[1])/59.0f;

						ti = 0;
						for (z4=0;z4<faces;z4++)
						{
							if (((CTetraFace*)m_oaFaces[z4])->XRay_Hit(tfy,tfn,tfx))
								ti++;
						}

						if (ti == 0)
							rhi[0]++;
						else if (ti == 1)
							rhi[1]++;
						else if (ti == 2)
							rhi[2]++;
						else rhi[3]++;

						rc++;

//						if ((ti != 0) && (ti != 2))
//							mprintf("\nError: Voronoi cell %d: %d ",id+1,ti);

//						if (ti == 0)
//							mprintf(".");
//						else mprintf("%d",ti);
					}
//					mprintf("\n");
				}
//				mprintf("\n");*/

			}
			cc++;
		} while (vl.inc());

	}

	if (verbose)
	{
		mprintf(WHITE,"]");
		mprintf(" done.\n");
	}

	if (sanity)
	{
		if (sanityall)
		{
			if ((m_fRayHisto[1] != 0) || (m_fRayHisto[3] != 0))
				b = false;

			for (z2=0;z2<m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube;z2++)
			{
				if (m_hitCount[z2] != 1) {
					b = false;
					break;
				}
// 				if (m_p3DF->m_pBin[z2] == 0)
// 				{
// 					b = false;
// 					break;
// 				}
// 				else if (m_p3DF->m_pBin[z2] == 1)
// 					m_fRayHisto[1]++;
// 				else 
// 				{
// 					b = false;
// 					break;
// 				}
			}

			if (!b)
			{
				mprintf("      %d Voronoi cells processed.\n",cc);
				mprintf("\n    %10d rays cast. Intersection histogram:\n",m_fRayCount);
				mprintf("    %10d rays hit   0 times (%10.6f%c).\n",m_fRayHisto[0],(double)m_fRayHisto[0]/m_fRayCount*100.0,'%');
				mprintf("    %10d rays hit   1 times (%10.6f%c).\n",m_fRayHisto[1],(double)m_fRayHisto[1]/m_fRayCount*100.0,'%');
				mprintf("    %10d rays hit   2 times (%10.6f%c).\n",m_fRayHisto[2],(double)m_fRayHisto[2]/m_fRayCount*100.0,'%');
				mprintf("    %10d rays hit > 2 times (%10.6f%c).\n",m_fRayHisto[3],(double)m_fRayHisto[3]/m_fRayCount*100.0,'%');

				mprintf("\n    Grid analysis:\n");
				m_fRayHisto[0] = 0;
				m_fRayHisto[1] = 0;
				m_fRayHisto[2] = 0;
				
				int z2a, z2b, z2c, z2d, z2e, z2f;
				for (z2a = 0; z2a < m_p3DF->m_iRes[0]; z2a++) {
					for (z2b = 0; z2b < m_p3DF->m_iRes[1]; z2b++) {
						for (z2c = 0; z2c < m_p3DF->m_iRes[2]; z2c++) {
							for (z2d = 0; z2d < m_iInterpolationFactor; z2d++) {
								for (z2e = 0; z2e < m_iInterpolationFactor; z2e++) {
									for (z2f = 0; z2f < m_iInterpolationFactor; z2f++) {
										if (m_hitCount[(z2d + z2e * m_iInterpolationFactor + z2f * m_iInterpolationFactor * m_iInterpolationFactor) * m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2] + z2a + z2b * m_p3DF->m_iRes[0] + z2c * m_p3DF->m_iResXY] == 0) {
											m_fRayHisto[0]++;
											b = false;
											mprintf("  Missed: ( %d.%d | %d.%d | %d.%d ).\n", z2a, z2d, z2b, z2e, z2c, z2f);
										} else if (m_hitCount[(z2d + z2e * m_iInterpolationFactor + z2f * m_iInterpolationFactor * m_iInterpolationFactor) * m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2] + z2a + z2b * m_p3DF->m_iRes[0] + z2c * m_p3DF->m_iResXY] == 1) {
											m_fRayHisto[1]++;
										} else {
											m_fRayHisto[2]++;
											b = false;
											mprintf("  Multiple: ( %d.%d | %d.%d | %d.%d ).\n", z2a, z2d, z2b, z2e, z2c, z2f);
										}
									}
								}
							}
						}
					}
				}
				
// 				for (z2=0;z2<m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube;z2++)
// 				{
// 					if (m_hitCount[z2] == 0) {
// 						m_fRayHisto[0]++;
// // 						mprintf("  Missed: ( %d | %d | %d ).\n", z2%m_p3DF->m_iRes[0], (z2 / m_p3DF->m_iRes[0]) % m_p3DF->m_iRes[1], (z2 / m_p3DF->m_iResXY));
// 					} else if (m_hitCount[z2] == 1) {
// 						m_fRayHisto[1]++;
// 					} else {
// 						m_fRayHisto[2]++;
// // 						mprintf("  Multiple: ( %d | %d | %d ).\n", z2%m_p3DF->m_iRes[0], (z2 / m_p3DF->m_iRes[0]) % m_p3DF->m_iRes[1], (z2 / m_p3DF->m_iResXY));
// 					}
// // 					if (m_p3DF->m_pBin[z2] == 0)
// // 					{
// // 						m_fRayHisto[0]++;
// // 						b = false;
// // 		//				mprintf("  Missed: ( %d | %d | %d ).\n",z2%100,(z2/100)%100,(z2/10000));
// // 					}
// // 					else if (m_p3DF->m_pBin[z2] == 1)
// // 						m_fRayHisto[1]++;
// // 					else 
// // 					{
// // 						m_fRayHisto[2]++;
// // 						b = false;
// // 		//				mprintf("  Multiple: ( %d | %d | %d ).\n",z2%100,(z2/100)%100,(z2/10000));
// // 					}
// 				}

				mprintf("      %10d points missed             (%10.6f%c).\n",m_fRayHisto[0],(double)m_fRayHisto[0]/(m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube)*100.0,'%');
				mprintf("      %10d points hit                (%10.6f%c).\n",m_fRayHisto[1],(double)m_fRayHisto[1]/(m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube)*100.0,'%');
				mprintf("      %10d points hit multiple times (%10.6f%c).\n",m_fRayHisto[2],(double)m_fRayHisto[2]/(m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube)*100.0,'%');
			}
		} else
		{
			mprintf("      %d Voronoi cells processed.\n",cc);
			mprintf("\n    %10d rays cast. Intersection histogram:\n",m_fRayCount);
			mprintf("    %10d rays hit   0 times (%10.6f%c).\n",m_fRayHisto[0],(double)m_fRayHisto[0]/m_fRayCount*100.0,'%');
			mprintf("    %10d rays hit   1 times (%10.6f%c).\n",m_fRayHisto[1],(double)m_fRayHisto[1]/m_fRayCount*100.0,'%');
			mprintf("    %10d rays hit   2 times (%10.6f%c).\n",m_fRayHisto[2],(double)m_fRayHisto[2]/m_fRayCount*100.0,'%');
			mprintf("    %10d rays hit > 2 times (%10.6f%c).\n",m_fRayHisto[3],(double)m_fRayHisto[3]/m_fRayCount*100.0,'%');

			if ((m_fRayHisto[1] != 0) || (m_fRayHisto[3] != 0))
				b = false;

			mprintf("\n    Grid analysis:\n");
			m_fRayHisto[0] = 0;
			m_fRayHisto[1] = 0;
			m_fRayHisto[2] = 0;
			
			int z2a, z2b, z2c, z2d, z2e, z2f;
			for (z2a = 0; z2a < m_p3DF->m_iRes[0]; z2a++) {
				for (z2b = 0; z2b < m_p3DF->m_iRes[1]; z2b++) {
					for (z2c = 0; z2c < m_p3DF->m_iRes[2]; z2c++) {
						for (z2d = 0; z2d < m_iInterpolationFactor; z2d++) {
							for (z2e = 0; z2e < m_iInterpolationFactor; z2e++) {
								for (z2f = 0; z2f < m_iInterpolationFactor; z2f++) {
									if (m_hitCount[(z2d + z2e * m_iInterpolationFactor + z2f * m_iInterpolationFactor * m_iInterpolationFactor) * m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2] + z2a + z2b * m_p3DF->m_iRes[0] + z2c * m_p3DF->m_iResXY] == 0) {
										m_fRayHisto[0]++;
										b = false;
										mprintf("  Missed: ( %d.%d | %d.%d | %d.%d ).\n", z2a, z2d, z2b, z2e, z2c, z2f);
									} else if (m_hitCount[(z2d + z2e * m_iInterpolationFactor + z2f * m_iInterpolationFactor * m_iInterpolationFactor) * m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2] + z2a + z2b * m_p3DF->m_iRes[0] + z2c * m_p3DF->m_iResXY] == 1) {
										m_fRayHisto[1]++;
									} else {
										m_fRayHisto[2]++;
										b = false;
										mprintf("  Multiple: ( %d.%d | %d.%d | %d.%d ).\n", z2a, z2d, z2b, z2e, z2c, z2f);
									}
								}
							}
						}
					}
				}
			}
			
// 			for (z2=0;z2<m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube;z2++)
// 			{
// // 				if (m_p3DF->m_pBin[z2] == 0)
// 				if (m_hitCount[z2] == 0)
// 				{
// 					m_fRayHisto[0]++;
// 					b = false;
// 					mprintf("  Missed: ( %d.%d | %d.%d | %d.%d ).\n", (z2 % m_iInterpolationFactorCube) % m_p3DF->m_iRes[0], (z2 / (m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2])) % m_iInterpolationFactor, ((z2 % m_iInterpolationFactorCube) / m_p3DF->m_iRes[0]) % m_p3DF->m_iRes[1], ((z2 / (m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2])) / m_iInterpolationFactor) % m_iInterpolationFactor, ((z2 % m_iInterpolationFactorCube) / m_p3DF->m_iResXY), ((z2 / (m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2])) / m_iInterpolationFactor * m_iInterpolationFactor));
// 				}
// // 				else if (m_p3DF->m_pBin[z2] == 1)
// 				else if (m_hitCount[z2] == 1)
// 					m_fRayHisto[1]++;
// 				else 
// 				{
// 					m_fRayHisto[2]++;
// 					b = false;
// 					mprintf("  Multiple: ( %d.%d | %d.%d | %d.%d ).\n", (z2 % m_iInterpolationFactorCube) % m_p3DF->m_iRes[0], (z2 / (m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2])) % m_iInterpolationFactor, ((z2 % m_iInterpolationFactorCube) / m_p3DF->m_iRes[0]) % m_p3DF->m_iRes[1], ((z2 / (m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2])) / m_iInterpolationFactor) % m_iInterpolationFactor, ((z2 % m_iInterpolationFactorCube) / m_p3DF->m_iResXY), ((z2 / (m_p3DF->m_iRes[0] * m_p3DF->m_iRes[1] * m_p3DF->m_iRes[2])) / m_iInterpolationFactor * m_iInterpolationFactor));
// 				}
// 			}

			mprintf("      %10d points missed             (%10.6f%c).\n",m_fRayHisto[0],(double)m_fRayHisto[0]/(m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube)*100.0,'%');
			mprintf("      %10d points hit                (%10.6f%c).\n",m_fRayHisto[1],(double)m_fRayHisto[1]/(m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube)*100.0,'%');
			mprintf("      %10d points hit multiple times (%10.6f%c).\n",m_fRayHisto[2],(double)m_fRayHisto[2]/(m_p3DF->m_iRes[0]*m_p3DF->m_iRes[1]*m_p3DF->m_iRes[2]*m_iInterpolationFactorCube)*100.0,'%');
		}

	}

	delete con;

	return b;
}


void CTetraPak::BuildVoronoiBuffer(CTimeStep *ts)
{
	voronoicell_neighbor c;
	container_periodic_poly *con;
	int ijk, q, z, z2, fc, id, ti;
	double *pp;
	vector<int> fv;
	CTetraCellBuffer *cbu;
	CTetraFaceBuffer *fbu;

	try { con = new container_periodic_poly(g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0,g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ,g_iVoroMemory); } catch(...) { con = NULL; }
	if (con == NULL) NewException((double)sizeof(container_periodic_poly),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<g_iGesAtomCount;z++)
 		con->put(z,ts->m_vaCoords[z][0]/1000.0,ts->m_vaCoords[z][1]/1000.0,ts->m_vaCoords[z][2]/1000.0,g_faVoronoiRadii[z]/1000.0);
	
	c_loop_all_periodic vl(*con);

	for (z=m_oaVoroBuffer.GetSize();z<g_iGesAtomCount;z++)
	{
		cbu = new CTetraCellBuffer();
		m_oaVoroBuffer.Add(cbu);
	}

	if (vl.start()) 
	{
		do 
		{
			if (con->compute_cell(c,vl))
			{
				ijk = vl.ijk;
				q = vl.q;
				pp = con->p[ijk]+con->ps*q;
				id = con->id[ijk][q];
				c.face_vertices(fv);

				cbu = (CTetraCellBuffer*)m_oaVoroBuffer[id];

				for (z=cbu->m_oaFaces.GetSize();z<c.number_of_faces();z++)
				{
					fbu = new CTetraFaceBuffer();
					cbu->m_oaFaces.Add(fbu);
				}

				cbu->m_iFaces = c.number_of_faces();
				cbu->m_vCenter = CxDVector3(pp[0]*1000.0,pp[1]*1000.0,pp[2]*1000.0);

				ti = 0;
				for (z=0;z<c.number_of_faces();z++)
				{
					fc = fv[ti];
					fbu = (CTetraFaceBuffer*)cbu->m_oaFaces[z];

					fbu->m_vaVertices.RemoveAll_KeepSize();

					for (z2=0;z2<fc;z2++)
						fbu->m_vaVertices.Add(
							CxDVector3(
								(*pp+c.pts[fv[ti+z2+1]*3]*0.5)*1000.0,
								(pp[1]+c.pts[fv[ti+z2+1]*3+1]*0.5)*1000.0,
								(pp[2]+c.pts[fv[ti+z2+1]*3+2]*0.5)*1000.0  )  );

					ti += fv[ti]+1;
				}
			}
		} while (vl.inc());
	}

	delete con;
}



void CTetraPak::Integrate_Verbose(C3DF<VORI_FLOAT> *df, double mi[2], double ma[2])
{
	int z2, z2b, z3, z3b, z4, z4b, ti, ixa, ixb, iya, iyb, iza, izb;
	double xv[10], tfy, tfz, tf;
	double EPS = 1.0E-20;
// 	int xi[10];

	iya = (int)floor(mi[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]+EPS);
	iyb = (int)ceil(ma[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]-EPS);

//	mprintf("Y: %f becomes %d, %f becomes %d.\n",mi[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]+EPS,iya,ma[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]-EPS,iyb);
//	mprintf(" MaxVal %f\n",df->m_fMaxVal[1]);
	
	iza = (int)floor(mi[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]+EPS);
	izb = (int)ceil(ma[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]-EPS);
	
//	mprintf(BLUE,"Integrating Y %d..%d, Z %d..%d.\n",iya,iyb,iza,izb);
//	mprintf(BLUE,"  (range Y %f..%f, Z %f..%f).\n",mi[0],ma[0],mi[1],ma[1]);

	if (iyb-iya >= df->m_iRes[1])
	{
//		mprintf(YELLOW,"\n    Averted calamity (y = %d..%d).",iya,iyb);
//		iyb = iya + df->m_iRes[1] - 1;
//		mprintf(YELLOW," Now %d .. %d.  ",iya,iyb);
	}

	if (izb-iza >= df->m_iRes[2])
	{
//		mprintf(YELLOW,"\n    Averted calamity (z = %d..%d).",iza,izb);
//		izb = iza + df->m_iRes[2] - 1;
//		mprintf(YELLOW," Now %d .. %d.  ",iza,izb);
	}

/*	mprintf("BoxX=%.20f, BoxY=%.20f, BoxZ=%.20f\n",g_fBoxX,g_fBoxY,g_fBoxZ);
	mprintf("Res0=%d, Res1=%d, Res2=%d.\n",df->m_iRes[0],df->m_iRes[1],df->m_iRes[2]);
	mprintf("mi0=%.20f, ma0=%.20f, mi1=%.20f, ma1=%.20f\n",mi[0],ma[0],mi[1],ma[1]);
	mprintf("iya=%d, iyb=%d, iza=%d, izb=%d\n",iya,iyb,iza,izb);*/

	for (z2=iya;z2<=iyb;z2++)
	{
// 		tfy = (float)(z2+0.0)*(g_fBoxY/1000.0)/df->m_iRes[1];
		tfy = (float)z2*(df->m_fMaxVal[1]/1000.0)/df->m_iRes[1];
		
		z2b = z2%df->m_iRes[1];
		if (z2b < 0)
			z2b += df->m_iRes[1];
		z2b *= df->m_iRes[0];

		for (z3=iza;z3<=izb;z3++)
		{
// 			tfz = (float)(z3+0.0)*(g_fBoxZ/1000.0)/df->m_iRes[2];
			tfz = (float)z3*(df->m_fMaxVal[2]/1000.0)/df->m_iRes[2];
			
			z3b = z3%df->m_iRes[2];
			if (z3b < 0)
				z3b += df->m_iRes[2];
			z3b *= df->m_iResXY;
			z3b += z2b;

			ti = 0;
			for (z4=0;z4<m_iFaces;z4++)
			{
				if (((CTetraFace*)m_oaFaces[z4])->XRay_Hit(tfy,tfz,xv[ti]/*,(z3==50)?true:false*/))
				{
	//				mprintf("  %3d | %3d [%2d]: Hit at (%d) %.6f\n",z2,z3,z4,ti,xv[ti]);
// 					xi[ti] = z4;
					ti++;
				} 
	//			else
	//			{
	//				mprintf("  %3d | %3d | %3d: Missed at (%d)\n",z2,z3,z4,ti);
	//			}
			}

	/*		if (ti != 2)
			{
				mprintf("***  %3d | %3d: Hit %d times (",z2,z3,ti);
				for (z4=0;z4<ti;z4++)
				{
					mprintf("%d@%f",xi[z4],xv[z4]);
					if (z4+1 == ti)
						mprintf(")\n");
					else
						mprintf(",");
				}
			}*/

			if (ti == 0)
				m_fRayHisto[0]++;
			else if (ti == 1)
				m_fRayHisto[1]++;
			else if (ti == 2)
				m_fRayHisto[2]++;
			else m_fRayHisto[3]++;

			m_fRayCount++;

			if (ti == 2)
			{
				if (xv[0] > xv[1])
				{
					tf = xv[0];
					xv[0] = xv[1];
					xv[1] = tf;
				}

				ixa = (int)ceil(xv[0] / (df->m_fMaxVal[1]/1000.0) * df->m_iRes[1]);
				ixb = (int)floor(xv[1] / (df->m_fMaxVal[1]/1000.0) * df->m_iRes[1]);

//				if ((z2==23) && (z3==50))
//					mprintf(YELLOW,"Hit 23|50: Range %d - %d\n",ixa,ixb);
				
//				if (ixa != ((int)ceil(xv[0] / (g_fBoxX / 1000.0) * df->m_iRes[0])))
//					mprintf("\n ixa %d <--> %d. ",ixa,(int)ceil(xv[0] / (g_fBoxX / 1000.0) * df->m_iRes[0]));

//				if ((z2==67) && (z3==51))
//					mprintf("\nixa=%3d (%8f) ixb=%3d (%8f):",ixa,xv[0],ixb,xv[1]);

				for (z4=ixa;z4<=ixb;z4++)
				{
					z4b = z4%df->m_iRes[0];
					if (z4b < 0)
						z4b += df->m_iRes[0];


//					if ((z2==67) && (z3==51))
//						mprintf(" %3d",z4b);

					df->m_pBin[z3b+z4b]++;

//					if ((z2==23) && (z3==50))
//						mprintf("%d... %f (%d)\n",z4b,df->m_pBin[z3b+z4b],z3b+z4b);

//					if (z3b+z4b == 502300)
//						mprintf(GREEN,"%d + %d = %d (z2=%d, z3=%d, z4=%d).\n",z3b,z4b,z3b+z4b,z2,z3,z4);
				}
			}
		}
	}
}


// VORI_FLOAT CTetraPak::Integrate(C3DF<VORI_FLOAT> *df, double mi[2], double ma[2])
// {
// 	int z2, z2b, z3, z3b, z4, z4b, ti, ixa, ixb, iya, iyb, iza, izb;
// 	double xv[4], tfy, tfz, tf, rv;
// 
// 	iya = (int)floor(mi[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
// 	iyb = (int)ceil(ma[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
// 	
// 	iza = (int)floor(mi[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
// 	izb = (int)ceil(ma[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
// 
// //	if (iyb-iya >= df->m_iRes[1])
// //		iyb = iya + df->m_iRes[1] - 1;
// 
// //	if (izb-iza >= df->m_iRes[2])
// //		izb = iza + df->m_iRes[2] - 1;
// 	
// 	rv = 0;
// 	for (z2=iya;z2<=iyb;z2++)
// 	{
// 		tfy = (float)z2*(df->m_fMaxVal[1]/1000.0)/df->m_iRes[1];
// 		
// 		z2b = z2%df->m_iRes[1];
// 		if (z2b < 0)
// 			z2b += df->m_iRes[1];
// 		z2b *= df->m_iRes[0];
// 
// 		for (z3=iza;z3<=izb;z3++)
// 		{
// 			tfz = (float)z3*(df->m_fMaxVal[2]/1000.0)/df->m_iRes[2];
// 			
// 			z3b = z3%df->m_iRes[2];
// 			if (z3b < 0)
// 				z3b += df->m_iRes[2];
// 			z3b *= df->m_iResXY;
// 			z3b += z2b;
// 
// 			ti = 0;
// 			for (z4=0;z4<m_iFaces;z4++)
// 				if (((CTetraFace*)m_oaFaces[z4])->XRay_Hit(tfy,tfz,xv[ti]))
// 					ti++;
// 
// 			if (ti == 2)
// 			{
// 				if (xv[0] > xv[1])
// 				{
// 					tf = xv[0];
// 					xv[0] = xv[1];
// 					xv[1] = tf;
// 				}
// 
// 				ixa = (int)ceil(xv[0] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]);
// 				
// 				ixb = (int)floor(xv[1] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]);
// 
// /*				if (g_bVoroIntEquitable)
// 				{
// 					tf = 0.5 * ((double)ixa - (xv[0] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]));
// 
// 					z4b = (ixa-1)%df->m_iRes[0];
// 					if (z4b < 0)
// 						z4b += df->m_iRes[0];
// 					
// 					rv += tf * df->m_pBin[z3b+z4b];
// 					
// 					z4b = ixa%df->m_iRes[0];
// 					if (z4b < 0)
// 						z4b += df->m_iRes[0];
// 
// 					rv += (0.5 + tf) * df->m_pBin[z3b+z4b];
// 
// 					ixa++;
// 
// 					tf = 0.5 * ((xv[1] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]) - (double)ixb);
// 
// 					z4b = ixb%df->m_iRes[0];
// 					if (z4b < 0)
// 						z4b += df->m_iRes[0];
// 					
// 					rv += (0.5 + tf) * df->m_pBin[z3b+z4b];
// 					
// 					z4b = (ixb+1)%df->m_iRes[0];
// 					if (z4b < 0)
// 						z4b += df->m_iRes[0];
// 
// 					rv += tf * df->m_pBin[z3b+z4b];
// 					
// 					ixb--;
// 				}*/
// 
// 				for (z4=ixa;z4<=ixb;z4++)
// 				{
// 					z4b = z4%df->m_iRes[0];
// 					if (z4b < 0)
// 						z4b += df->m_iRes[0];
// 
// 					rv += df->m_pBin[z3b+z4b];
// 				}
// 			}
// 		}
// 	}
// 
// 	return rv;
// }


VORI_FLOAT CTetraPak::Integrate_Refine(C3DF<VORI_FLOAT> *df, double mi[2], double ma[2])
{
	int z2, z2i, z3, z3i, z4, ti, ixa, ixb, iya, iyb, iza, izb;
	double xv[4], tfy, tfz, tf;
	int iz1, iz2, iy1, iy2, ix1, ix2;
	VORI_FLOAT rv;
	VORI_FLOAT cz1, cz2, cy1, cy2, cx1, cx2;
	VORI_FLOAT c11, c12, c21, c22;


	iya = (int)floor(mi[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
	iyb = (int)ceil(ma[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
	
	iza = (int)floor(mi[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
	izb = (int)ceil(ma[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);

//	if (iyb-iya >= df->m_iRes[1])
//		iyb = iya + df->m_iRes[1] - 1;

//	if (izb-iza >= df->m_iRes[2])
//		izb = iza + df->m_iRes[2] - 1;
	
	rv = 0;
	for (z2=iya;z2<=iyb;z2++)
	{		
		iy1 = z2%df->m_iRes[1];
		if (iy1 < 0)
			iy1 += df->m_iRes[1];
		iy1 *= df->m_iRes[0];

		iy2 = (z2+1)%df->m_iRes[1];
		if (iy2 < 0)
			iy2 += df->m_iRes[1];
		iy2 *= df->m_iRes[0];

		for (z3=iza;z3<=izb;z3++)
		{
			iz1 = z3%df->m_iRes[2];
			if (iz1 < 0)
				iz1 += df->m_iRes[2];
			iz1 *= df->m_iResXY;

			iz2 = (z3+1)%df->m_iRes[2];
			if (iz2 < 0)
				iz2 += df->m_iRes[2];
			iz2 *= df->m_iResXY;

			for (z2i=0;z2i<m_iInterpolationFactor;z2i++)
			{
	//			tfy = (z2-0.5+(z2i+0.5)/(double)m_iInterpolationFactor)/df->m_iRes[1]*df->m_fMaxVal[1]/1000.0;
				tfy = (z2+z2i/(double)m_iInterpolationFactor)/df->m_iRes[1]*df->m_fMaxVal[1]/1000.0;

				cy2 = z2i/(VORI_FLOAT)m_iInterpolationFactor;
				cy1 = (VORI_FLOAT)(1.0 - cy2);

				for (z3i=0;z3i<m_iInterpolationFactor;z3i++)
				{
//					tfz = (z3-0.5+(z3i+0.5)/(double)m_iInterpolationFactor)/df->m_iRes[2]*df->m_fMaxVal[2]/1000.0;
					tfz = (z3+z3i/(double)m_iInterpolationFactor)/df->m_iRes[2]*df->m_fMaxVal[2]/1000.0;

					cz2 = z3i/(VORI_FLOAT)m_iInterpolationFactor;
					cz1 = (VORI_FLOAT)(1.0 - cz2);

					c11 = cy1*cz1;
					c12 = cy1*cz2;
					c21 = cy2*cz1;
					c22 = cy2*cz2;

					ti = 0;
					for (z4=0;z4<m_iFaces;z4++)
						if (((CTetraFace*)m_oaFaces[z4])->XRay_Hit(tfy,tfz,xv[ti]))
							ti++;

					if (ti == 2)
					{
						if (xv[0] > xv[1])
						{
							tf = xv[0];
							xv[0] = xv[1];
							xv[1] = tf;
						}

						ixa = (int)ceil(xv[0] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0] * m_iInterpolationFactor);
						ixb = (int)floor(xv[1] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0] * m_iInterpolationFactor);

						for (z4=ixa;z4<=ixb;z4++)
						{
							ix1 = (z4/m_iInterpolationFactor)%df->m_iRes[0];
							if (ix1 < 0)
								ix1 += df->m_iRes[0];

							ix2 = (ix1+1)%df->m_iRes[0];

							cx2 = (z4%m_iInterpolationFactor)/(VORI_FLOAT)m_iInterpolationFactor;
							cx1 = (VORI_FLOAT)(1.0 - cx2);

							rv += cx1*(c11*df->m_pBin[ix1+iy1+iz1]
								+ c12*df->m_pBin[ix1+iy1+iz2]
								+ c21*df->m_pBin[ix1+iy2+iz1]
								+ c22*df->m_pBin[ix1+iy2+iz2]);

							rv += cx2*(c11*df->m_pBin[ix2+iy1+iz1]
								+ c12*df->m_pBin[ix2+iy1+iz2]
								+ c21*df->m_pBin[ix2+iy2+iz1]
								+ c22*df->m_pBin[ix2+iy2+iz2]);
						}
					}
				}
			}
		}
	}

	return rv / (VORI_FLOAT)(m_iInterpolationFactorCube);
}

void CTetraPak::integrateCell(C3DF<VORI_FLOAT> *electronDensity, CxFloatArray *currentDensity, double mi[2], double ma[2], double refPoint[3], double *charge, CxDVector3 *dipoleMoment, CxDVector3 *totalCurrent, CxDVector3 *magneticMoment, bool sanity/*, int debug*/) {
	int z2, z3, z4, ti, ixa, ixb, iya, iyb, iza, izb;
	int ix1, ix2, iy1, iy2, iz1, iz2, z2i, z3i, z4a;
	double xv[4], tfy, tfz, tf;
	VORI_FLOAT cx1, cx2, cy1, cy2, cz1, cz2, c11, c12, c21, c22;
	
	*charge = 0.0;
	*dipoleMoment = CxDVector3(0.0, 0.0, 0.0);
	*totalCurrent = CxDVector3(0.0, 0.0, 0.0);
	*magneticMoment = CxDVector3(0.0, 0.0, 0.0);
	
	if (g_bBoxNonOrtho) {
		iya = (int)floor(mi[0] * electronDensity->m_iRes[1]);
		iyb = (int)ceil(ma[0] * electronDensity->m_iRes[1]);
		
		iza = (int)floor(mi[1] * electronDensity->m_iRes[2]);
		izb = (int)ceil(ma[1] * electronDensity->m_iRes[2]);
	} else {
		iya = (int)floor(mi[0] / (electronDensity->m_fMaxVal[1] / 1000.0) * electronDensity->m_iRes[1]);
		iyb = (int)ceil(ma[0] / (electronDensity->m_fMaxVal[1] / 1000.0) * electronDensity->m_iRes[1]);
		
		iza = (int)floor(mi[1] / (electronDensity->m_fMaxVal[2] / 1000.0) * electronDensity->m_iRes[2]);
		izb = (int)ceil(ma[1] / (electronDensity->m_fMaxVal[2] / 1000.0) * electronDensity->m_iRes[2]);
	}
	
	for (z2 = iya; z2 <= iyb; z2++) {
// 		tfy = (double)z2 * (electronDensity->m_fMaxVal[1] / 1000.0) / electronDensity->m_iRes[1];
		
// 		z2b = z2 % electronDensity->m_iRes[1];
// 		if (z2b < 0)
// 			z2b += electronDensity->m_iRes[1];
// 		z2b *= electronDensity->m_iRes[0];
		
		iy1 = z2 % electronDensity->m_iRes[1];
		if (iy1 < 0)
			iy1 += electronDensity->m_iRes[1];
		iy1 *= electronDensity->m_iRes[0];
		
		iy2 = (z2 + 1) % electronDensity->m_iRes[1];
		if (iy2 < 0)
			iy2 += electronDensity->m_iRes[1];
		iy2 *= electronDensity->m_iRes[0];
		
		for (z3 = iza; z3 <= izb; z3++) {
// 			tfz = (double)z3 * (electronDensity->m_fMaxVal[2] / 1000.0) / electronDensity->m_iRes[2];
			
// 			z3b = z3 % electronDensity->m_iRes[2];
// 			if (z3b < 0)
// 				z3b += electronDensity->m_iRes[2];
// 			z3b *= electronDensity->m_iResXY;
// 			z3b += z2b;
			
			iz1 = z3 % electronDensity->m_iRes[2];
			if (iz1 < 0)
				iz1 += electronDensity->m_iRes[2];
			iz1 *= electronDensity->m_iResXY;
			
			iz2 = (z3 + 1) % electronDensity->m_iRes[2];
			if (iz2 < 0)
				iz2 += electronDensity->m_iRes[2];
			iz2 *= electronDensity->m_iResXY;
			
// 			mprintf(GREEN, "Hallo %d\n", m_iInterpolationFactor);
			for (z2i = 0; z2i < m_iInterpolationFactor; z2i++) {
// 			for (z2i = 0; z2i < 1; z2i++) {
				if (g_bBoxNonOrtho) {
					tfy = ((double)z2 + (double)z2i / m_iInterpolationFactor) / electronDensity->m_iRes[1];
				} else {
					tfy = ((double)z2 + (double)z2i / m_iInterpolationFactor) * (electronDensity->m_fMaxVal[1] / 1000.0) / electronDensity->m_iRes[1];
				}
				
				cy2 = (VORI_FLOAT)z2i / m_iInterpolationFactor;
				cy1 = 1.0 - cy2;
				
				for (z3i = 0; z3i < m_iInterpolationFactor; z3i++) {
// 				for (z3i = 0; z3i < 1; z3i++) {
					if (g_bBoxNonOrtho) {
						tfz = ((double)z3 + (double)z3i / m_iInterpolationFactor) / electronDensity->m_iRes[2];
					} else {
						tfz = ((double)z3 + (double)z3i / m_iInterpolationFactor) * (electronDensity->m_fMaxVal[2] / 1000.0) / electronDensity->m_iRes[2];
					}
					
					cz2 = (VORI_FLOAT)z3i / m_iInterpolationFactor;
					cz1 = 1.0 - cz2;
					
					c11 = cy1 * cz1;
					c12 = cy1 * cz2;
					c21 = cy2 * cz1;
					c22 = cy2 * cz2;
					
					ti = 0;
					for (z4 = 0; z4 < m_iFaces; z4++) {
						if (((CTetraFace *)m_oaFaces[z4])->XRay_Hit(tfy, tfz, xv[ti]/*, debug == 103 && fabs(tfy - 0.16666666666666665741) < 1.0e-8 && fabs(tfz - 0.46239999999999997771) < 1.0e-8*/)) {
// 							if (debug == 103 && fabs(tfy - 0.16666666666666665741) < 1.0e-8 && fabs(tfz - 0.46239999999999997771) < 1.0e-8)
// 								mprintf(GREEN, "Hit %d: %.20g\n", z4+1, xv[ti]);
							ti++;
						}
					}
						
					if (sanity) {
						if (ti == 0)
							m_fRayHisto[0]++;
						else if (ti == 1)
							m_fRayHisto[1]++;
						else if (ti == 2)
							m_fRayHisto[2]++;
						else m_fRayHisto[3]++;
						
						m_fRayCount++;
					}
					
					if (ti == 2) {
						if (xv[0] > xv[1]) {
							tf = xv[0];
							xv[0] = xv[1];
							xv[1] = tf;
						}
						
// 						while (fabs(xv[0] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] - ceil(xv[0] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0])) < VORI_EPS)
// 							xv[0] += VORI_EPS * (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0];
// 						while (fabs(xv[1] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] - ceil(xv[1] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0])) < VORI_EPS)
// 							xv[1] += VORI_EPS * (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0];
						
						if (g_bBoxNonOrtho) {
							while (fabs(xv[0] * electronDensity->m_iRes[0] * m_iInterpolationFactor - ceil(xv[0] * electronDensity->m_iRes[0] * m_iInterpolationFactor)) < VORI_EPS)
								xv[0] += VORI_EPS/* * electronDensity->m_iRes[0] * m_iInterpolationFactor*/;
							while (fabs(xv[1] * electronDensity->m_iRes[0] * m_iInterpolationFactor - ceil(xv[1] * electronDensity->m_iRes[0] * m_iInterpolationFactor)) < VORI_EPS)
								xv[1] += VORI_EPS/* * electronDensity->m_iRes[0] * m_iInterpolationFactor*/;
							
							ixa = (int)ceil(xv[0] * electronDensity->m_iRes[0] * m_iInterpolationFactor);
							ixb = (int)floor(xv[1] * electronDensity->m_iRes[0] * m_iInterpolationFactor);
						} else {
							while (fabs(xv[0] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] * m_iInterpolationFactor - ceil(xv[0] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] * m_iInterpolationFactor)) < VORI_EPS)
								xv[0] += VORI_EPS * (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] * m_iInterpolationFactor;
							while (fabs(xv[1] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] * m_iInterpolationFactor - ceil(xv[1] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] * m_iInterpolationFactor)) < VORI_EPS)
								xv[1] += VORI_EPS * (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] * m_iInterpolationFactor;
							
							ixa = (int)ceil(xv[0] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] * m_iInterpolationFactor);
							ixb = (int)floor(xv[1] / (electronDensity->m_fMaxVal[0] / 1000.0) * electronDensity->m_iRes[0] * m_iInterpolationFactor);
						}
						
						for (z4 = ixa; z4 <= ixb; z4++) {
// 							z4b = z4 % electronDensity->m_iRes[0];
// 							if (z4b < 0)
// 								z4b += electronDensity->m_iRes[0];
							
							ix1 = floori(z4, m_iInterpolationFactor) % electronDensity->m_iRes[0];
							if (ix1 < 0)
								ix1 += electronDensity->m_iRes[0];
							ix2 = (ix1 + 1) % electronDensity->m_iRes[0];
							
							z4a = z4 % m_iInterpolationFactor;
							if (z4a < 0)
								z4a += m_iInterpolationFactor;
							
							cx2 = (VORI_FLOAT)z4a / m_iInterpolationFactor;
							cx1 = 1.0 - cx2;
							
							double x, y, z;
							if (g_bBoxNonOrtho) {
								x = (((double)z4 / m_iInterpolationFactor) / electronDensity->m_iRes[0] * g_mCubeCell(0,0) + tfy * g_mCubeCell(1,0) + tfz * g_mCubeCell(2,0)) / 1000.0;
								y = (((double)z4 / m_iInterpolationFactor) / electronDensity->m_iRes[0] * g_mCubeCell(0,1) + tfy * g_mCubeCell(1,1) + tfz * g_mCubeCell(2,1)) / 1000.0;
								z = (((double)z4 / m_iInterpolationFactor) / electronDensity->m_iRes[0] * g_mCubeCell(0,2) + tfy * g_mCubeCell(1,2) + tfz * g_mCubeCell(2,2)) / 1000.0;
							} else {
								x = ((double)z4 / m_iInterpolationFactor) * (electronDensity->m_fMaxVal[0] / 1000.0) / electronDensity->m_iRes[0];
								y = tfy;
								z = tfz;
							}
// 							if (ix1 == 0 && iy1 == 0 && iz1 == 0) {
// 								mprintf(GREEN, "%d %d %d\n", z4, z2, z3);
// 								mprintf(GREEN, "%.10g %.10g %.10g\n", x, y, z);
// 								mprintf(GREEN, "%d %d %d\n", z4a, z2i, z3i);
// 								mprintf(GREEN, "%d %d\n", (z4a + z2i * m_iInterpolationFactor + z3i * m_iInterpolationFactor * m_iInterpolationFactor) * electronDensity->m_iRes[0] * electronDensity->m_iRes[1] * electronDensity->m_iRes[2] + ix1 + iy1 + iz1, m_hitCount[(z4a + z2i * m_iInterpolationFactor + z3i * m_iInterpolationFactor * m_iInterpolationFactor) * electronDensity->m_iRes[0] * electronDensity->m_iRes[1] * electronDensity->m_iRes[2] + ix1 + iy1 + iz1]);
// 							}
							
							if (sanity) {
// 								mprintf(GREEN, "%d\n", (z4a + z2i * m_iInterpolationFactor + z3i * m_iInterpolationFactor * m_iInterpolationFactor) * electronDensity->m_iRes[0] * electronDensity->m_iRes[1] * electronDensity->m_iRes[2] + ix1 + iy1 + iz1);
								m_hitCount[(z4a + z2i * m_iInterpolationFactor + z3i * m_iInterpolationFactor * m_iInterpolationFactor) * electronDensity->m_iRes[0] * electronDensity->m_iRes[1] * electronDensity->m_iRes[2] + ix1 + iy1 + iz1] += 1;
							} else {
								double ed = (cx1 * (c11 * electronDensity->m_pBin[ix1+iy1+iz1] + c12 * electronDensity->m_pBin[ix1+iy1+iz2] + c21 * electronDensity->m_pBin[ix1+iy2+iz1] + c22 * electronDensity->m_pBin[ix1+iy2+iz2])
								+ cx2 * (c11 * electronDensity->m_pBin[ix2+iy1+iz1] + c12 * electronDensity->m_pBin[ix2+iy1+iz2] + c21 * electronDensity->m_pBin[ix2+iy2+iz1] + c22 * electronDensity->m_pBin[ix2+iy2+iz2])) / m_iInterpolationFactorCube;
								if (g_bVoroIntegrateCharge) {
// 									*charge += electronDensity->m_pBin[z3b+z4b];
									*charge += ed;
								}
								if (g_bVoroIntegrateDipoleMoment) {
// 									(*dipoleMoment)[0] += -(x - refPoint[0]) * electronDensity->m_pBin[z3b+z4b];
// 									(*dipoleMoment)[1] += -(tfy - refPoint[1]) * electronDensity->m_pBin[z3b+z4b];
// 									(*dipoleMoment)[2] += -(tfz - refPoint[2]) * electronDensity->m_pBin[z3b+z4b];
									(*dipoleMoment)[0] += -(x - refPoint[0]) * ed;
									(*dipoleMoment)[1] += -(y - refPoint[1]) * ed;
									(*dipoleMoment)[2] += -(z - refPoint[2]) * ed;
								}
								if (currentDensity != NULL) {
									double cdx = (cx1 * (c11 * currentDensity->GetAt((ix1+iy1+iz1)*3) + c12 * currentDensity->GetAt((ix1+iy1+iz2)*3) + c21 * currentDensity->GetAt((ix1+iy2+iz1)*3) + c22 * currentDensity->GetAt((ix1+iy2+iz2)*3))
									+ cx2 * (c11 * currentDensity->GetAt((ix2+iy1+iz1)*3) + c12 * currentDensity->GetAt((ix2+iy1+iz2)*3) + c21 * currentDensity->GetAt((ix2+iy2+iz1)*3) + c22 * currentDensity->GetAt((ix2+iy2+iz2)*3))) / m_iInterpolationFactorCube;
									double cdy = (cx1 * (c11 * currentDensity->GetAt((ix1+iy1+iz1)*3 + 1) + c12 * currentDensity->GetAt((ix1+iy1+iz2)*3 + 1) + c21 * currentDensity->GetAt((ix1+iy2+iz1)*3 + 1) + c22 * currentDensity->GetAt((ix1+iy2+iz2)*3 + 1))
									+ cx2 * (c11 * currentDensity->GetAt((ix2+iy1+iz1)*3 + 1) + c12 * currentDensity->GetAt((ix2+iy1+iz2)*3 + 1) + c21 * currentDensity->GetAt((ix2+iy2+iz1)*3 + 1) + c22 * currentDensity->GetAt((ix2+iy2+iz2)*3 + 1))) / m_iInterpolationFactorCube;
									double cdz = (cx1 * (c11 * currentDensity->GetAt((ix1+iy1+iz1)*3 + 2) + c12 * currentDensity->GetAt((ix1+iy1+iz2)*3 + 2) + c21 * currentDensity->GetAt((ix1+iy2+iz1)*3 + 2) + c22 * currentDensity->GetAt((ix1+iy2+iz2)*3 + 2))
									+ cx2 * (c11 * currentDensity->GetAt((ix2+iy1+iz1)*3 + 2) + c12 * currentDensity->GetAt((ix2+iy1+iz2)*3 + 2) + c21 * currentDensity->GetAt((ix2+iy2+iz1)*3 + 2) + c22 * currentDensity->GetAt((ix2+iy2+iz2)*3 + 2))) / m_iInterpolationFactorCube;
									if (g_bVoroIntegrateTotalCurrent) {
// 										(*totalCurrent)[0] += currentDensity->GetAt((z3b+z4b)*3);
// 										(*totalCurrent)[1] += currentDensity->GetAt((z3b+z4b)*3 + 1);
// 										(*totalCurrent)[2] += currentDensity->GetAt((z3b+z4b)*3 + 2);
										(*totalCurrent)[0] += cdx;
										(*totalCurrent)[1] += cdy;
										(*totalCurrent)[2] += cdz;
									}
									if (g_bVoroIntegrateMagneticMoment) {
// 										(*magneticMoment)[0] += (tfy - refPoint[1]) * currentDensity->GetAt((z3b+z4b)*3 + 2);
// 										(*magneticMoment)[0] -= (tfz - refPoint[2]) * currentDensity->GetAt((z3b+z4b)*3 + 1);
// 										(*magneticMoment)[1] += (tfz - refPoint[2]) * currentDensity->GetAt((z3b+z4b)*3);
// 										(*magneticMoment)[1] -= (x - refPoint[0]) * currentDensity->GetAt((z3b+z4b)*3 + 2);
// 										(*magneticMoment)[2] += (x - refPoint[0]) * currentDensity->GetAt((z3b+z4b)*3 + 1);
// 										(*magneticMoment)[2] -= (tfy - refPoint[1]) * currentDensity->GetAt((z3b+z4b)*3);
										(*magneticMoment)[0] += (y - refPoint[1]) * cdz;
										(*magneticMoment)[0] -= (z - refPoint[2]) * cdy;
										(*magneticMoment)[1] += (z - refPoint[2]) * cdx;
										(*magneticMoment)[1] -= (x - refPoint[0]) * cdz;
										(*magneticMoment)[2] += (x - refPoint[0]) * cdy;
										(*magneticMoment)[2] -= (y - refPoint[1]) * cdx;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

// CxDVector3 CTetraPak::integrateMoment(C3DF<VORI_FLOAT>* df, double mi[2], double ma[2], double refPoint[3]) {
// 	int z2, z2b, z3, z3b, z4, z4b, ti, ixa, ixb, iya, iyb, iza, izb;
// 	double xv[4], tfy, tfz, tf;
// 	CxDVector3 rv(0.0, 0.0, 0.0);
// 	
// // 	iya = (int)(mi[0]/(g_fBoxY/1000.0)*df->m_iRes[1] - 0.5);
// // 	iya = (int)floor(mi[0]/(g_fBoxY/1000.0)*df->m_iRes[1]);
// 	iya = (int)floor(mi[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
// // 	iyb = (int)(ma[0]/(g_fBoxY/1000.0)*df->m_iRes[1] + 0.5);
// // 	iyb = (int)ceil(ma[0]/(g_fBoxY/1000.0)*df->m_iRes[1]);
// 	iyb = (int)ceil(ma[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
// 	
// // 	iza = (int)(mi[1]/(g_fBoxZ/1000.0)*df->m_iRes[2] - 0.5);
// // 	iza = (int)floor(mi[1]/(g_fBoxZ/1000.0)*df->m_iRes[2]);
// 	iza = (int)floor(mi[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
// // 	izb = (int)(ma[1]/(g_fBoxZ/1000.0)*df->m_iRes[2] + 0.5);
// // 	izb = (int)ceil(ma[1]/(g_fBoxZ/1000.0)*df->m_iRes[2]);
// 	izb = (int)ceil(ma[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
// 
// //	if (iyb-iya >= df->m_iRes[1])
// //		iyb = iya + df->m_iRes[1] - 1;
// 
// //	if (izb-iza >= df->m_iRes[2])
// //		izb = iza + df->m_iRes[2] - 1;
// 	
// 	for (z2=iya;z2<=iyb;z2++)
// 	{
// // 		tfy = (double)z2 * (g_fBoxY / 1000.0) / df->m_iRes[1];
// 		tfy = (double)z2 * (df->m_fMaxVal[1] / 1000.0) / df->m_iRes[1];
// 
// 		z2b = z2%df->m_iRes[1];
// 		if (z2b < 0)
// 			z2b += df->m_iRes[1];
// 		z2b *= df->m_iRes[0];
// 		
// 		for (z3=iza;z3<=izb;z3++)
// 		{
// // 			tfz = (double)z3 * (g_fBoxZ / 1000.0) / df->m_iRes[2];
// 			tfz = (double)z3 * (df->m_fMaxVal[2] / 1000.0) / df->m_iRes[2];
// 
// 			z3b = z3%df->m_iRes[2];
// 			if (z3b < 0)
// 				z3b += df->m_iRes[2];
// 			z3b *= df->m_iResXY;
// 			z3b += z2b;
// 			
// 			ti = 0;
// 			for (z4=0;z4<m_iFaces;z4++)
// 				if (((CTetraFace*)m_oaFaces[z4])->XRay_Hit(tfy,tfz,xv[ti]))
// 					ti++;
// 				
// 			if (ti == 2)
// 			{
// 				if (xv[0] > xv[1])
// 				{
// 					tf = xv[0];
// 					xv[0] = xv[1];
// 					xv[1] = tf;
// 				}
// 				
// // 				if ((xv[0]/(g_fBoxX/1000.0)*df->m_iRes[0]+0.5) < 0)
// // 					ixa = (int)((xv[0]+g_fBoxX/1000.0)/(g_fBoxX/1000.0)*df->m_iRes[0] + 0.5) - df->m_iRes[0];
// // 				else ixa = (int)(xv[0]/(g_fBoxX/1000.0)*df->m_iRes[0] + 0.5);
// // 				ixa = (int)ceil(xv[0] / (g_fBoxX / 1000.0) * df->m_iRes[0]);
// 				ixa = (int)ceil(xv[0] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]);
// 				
// // 				if ((xv[1]/(g_fBoxX/1000.0)*df->m_iRes[0] - 0.5) < 0)
// // 					ixb = (int)((xv[1]+g_fBoxX/1000.0)/(g_fBoxX/1000.0)*df->m_iRes[0] - 0.5) - df->m_iRes[0];
// // 				else ixb = (int)(xv[1]/(g_fBoxX/1000.0)*df->m_iRes[0] - 0.5);
// // 				ixb = (int)floor(xv[1] / (g_fBoxX / 1000.0) * df->m_iRes[0]);
// 				ixb = (int)floor(xv[1] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]);
// 				
// 				for (z4=ixa;z4<=ixb;z4++)
// 				{
// 					z4b = z4%df->m_iRes[0];
// 					if (z4b < 0)
// 						z4b += df->m_iRes[0];
// 					double x = (double)z4 * (df->m_fMaxVal[0] / 1000.0) / df->m_iRes[0];
// 					
// 					rv[0] += -(x - refPoint[0]) * df->m_pBin[z3b+z4b];
// 					rv[1] += -(tfy - refPoint[1]) * df->m_pBin[z3b+z4b];
// 					rv[2] += -(tfz - refPoint[2]) * df->m_pBin[z3b+z4b];
// 				}
// 			}
// 		}
// 	}
// 	
// 	return rv;
// }

// CxDVector3 CTetraPak::integrateTotalCurrent(C3DF<VORI_FLOAT> *df, CxFloatArray *currentDensity, double mi[2], double ma[2]) {
// 	int z2, z2b, z3, z3b, z4, z4b, ti, ixa, ixb, iya, iyb, iza, izb;
// 	double xv[4], tfy, tfz, tf;
// 	CxDVector3 rv(0.0, 0.0, 0.0);
// 	
// 	iya = (int)floor(mi[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
// 	iyb = (int)ceil(ma[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
// 	
// 	iza = (int)floor(mi[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
// 	izb = (int)ceil(ma[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
// 	
// 	for (z2=iya;z2<=iyb;z2++)
// 	{
// 		tfy = (double)z2 * (df->m_fMaxVal[1] / 1000.0) / df->m_iRes[1];
// 		
// 		z2b = z2%df->m_iRes[1];
// 		if (z2b < 0)
// 			z2b += df->m_iRes[1];
// 		z2b *= df->m_iRes[0];
// 		
// 		for (z3=iza;z3<=izb;z3++)
// 		{
// 			tfz = (double)z3 * (df->m_fMaxVal[2] / 1000.0) / df->m_iRes[2];
// 			
// 			z3b = z3%df->m_iRes[2];
// 			if (z3b < 0)
// 				z3b += df->m_iRes[2];
// 			z3b *= df->m_iResXY;
// 			z3b += z2b;
// 			
// 			ti = 0;
// 			for (z4=0;z4<m_iFaces;z4++)
// 				if (((CTetraFace*)m_oaFaces[z4])->XRay_Hit(tfy,tfz,xv[ti]))
// 					ti++;
// 				
// 			if (ti == 2)
// 			{
// 				if (xv[0] > xv[1])
// 				{
// 					tf = xv[0];
// 					xv[0] = xv[1];
// 					xv[1] = tf;
// 				}
// 				
// 				ixa = (int)ceil(xv[0] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]);
// 				ixb = (int)floor(xv[1] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]);
// 				
// 				for (z4=ixa;z4<=ixb;z4++)
// 				{
// 					z4b = z4%df->m_iRes[0];
// 					if (z4b < 0)
// 						z4b += df->m_iRes[0];
// 					
// 					rv[0] += currentDensity->GetAt((z3b+z4b)*3);
// 					rv[1] += currentDensity->GetAt((z3b+z4b)*3 + 1);
// 					rv[2] += currentDensity->GetAt((z3b+z4b)*3 + 2);
// 				}
// 			}
// 		}
// 	}
// 	
// 	return rv;
// }
// 
// CxDVector3 CTetraPak::integrateMagneticMoment(C3DF<VORI_FLOAT> *df, CxFloatArray *currentDensity, double mi[2], double ma[2], double refPoint[3]) {
// 	int z2, z2b, z3, z3b, z4, z4b, ti, ixa, ixb, iya, iyb, iza, izb;
// 	double xv[4], tfy, tfz, tf;
// 	CxDVector3 rv(0.0, 0.0, 0.0);
// 	
// 	iya = (int)floor(mi[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
// 	iyb = (int)ceil(ma[0]/(df->m_fMaxVal[1]/1000.0)*df->m_iRes[1]);
// 	
// 	iza = (int)floor(mi[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
// 	izb = (int)ceil(ma[1]/(df->m_fMaxVal[2]/1000.0)*df->m_iRes[2]);
// 	
// 	for (z2=iya;z2<=iyb;z2++)
// 	{
// 		tfy = (double)z2 * (df->m_fMaxVal[1] / 1000.0) / df->m_iRes[1];
// 		
// 		z2b = z2%df->m_iRes[1];
// 		if (z2b < 0)
// 			z2b += df->m_iRes[1];
// 		z2b *= df->m_iRes[0];
// 		
// 		for (z3=iza;z3<=izb;z3++)
// 		{
// 			tfz = (double)z3 * (df->m_fMaxVal[2] / 1000.0) / df->m_iRes[2];
// 			
// 			z3b = z3%df->m_iRes[2];
// 			if (z3b < 0)
// 				z3b += df->m_iRes[2];
// 			z3b *= df->m_iResXY;
// 			z3b += z2b;
// 			
// 			ti = 0;
// 			for (z4=0;z4<m_iFaces;z4++)
// 				if (((CTetraFace*)m_oaFaces[z4])->XRay_Hit(tfy,tfz,xv[ti]))
// 					ti++;
// 				
// 			if (ti == 2)
// 			{
// 				if (xv[0] > xv[1])
// 				{
// 					tf = xv[0];
// 					xv[0] = xv[1];
// 					xv[1] = tf;
// 				}
// 				
// 				ixa = (int)ceil(xv[0] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]);
// 				ixb = (int)floor(xv[1] / (df->m_fMaxVal[0] / 1000.0) * df->m_iRes[0]);
// 				
// 				for (z4=ixa;z4<=ixb;z4++)
// 				{
// 					z4b = z4%df->m_iRes[0];
// 					if (z4b < 0)
// 						z4b += df->m_iRes[0];
// 					double x = (double)z4 * (df->m_fMaxVal[0] / 1000.0) / df->m_iRes[0];
// 					
// 					rv[0] += (tfy - refPoint[1]) * currentDensity->GetAt((z3b+z4b)*3 + 2);
// 					rv[0] -= (tfz - refPoint[2]) * currentDensity->GetAt((z3b+z4b)*3 + 1);
// 					rv[1] += (tfz - refPoint[2]) * currentDensity->GetAt((z3b+z4b)*3);
// 					rv[1] -= (x - refPoint[0]) * currentDensity->GetAt((z3b+z4b)*3 + 2);
// 					rv[2] += (x - refPoint[0]) * currentDensity->GetAt((z3b+z4b)*3 + 1);
// 					rv[2] += (tfy - refPoint[1]) * currentDensity->GetAt((z3b+z4b)*3);
// 				}
// 			}
// 		}
// 	}
// 	
// 	return rv;
// }


inline bool CTetraFace::XRay_Hit(double py, double pz, double &px/*, bool verbose*/) {
	double tf, a, b, ty, tz;
	int z;
	
	#define lvec1x (m_vSpan1[0])
	#define lvec1y (m_vSpan1[1])
	#define lvec1z (m_vSpan1[2])
	#define lvec2x (m_vSpan2[0])
	#define lvec2y (m_vSpan2[1])
	#define lvec2z (m_vSpan2[2])
	
	if (g_bBoxNonOrtho) {
		
		CxDVector3 veca = CxDVector3(g_mCubeCell(0,0) / 1000.0, g_mCubeCell(0,1) / 1000.0, g_mCubeCell(0,2) / 1000.0);
		CxDVector3 vecb = CxDVector3(g_mCubeCell(1,0) / 1000.0, g_mCubeCell(1,1) / 1000.0, g_mCubeCell(1,2) / 1000.0);
		CxDVector3 vecc = CxDVector3(g_mCubeCell(2,0) / 1000.0, g_mCubeCell(2,1) / 1000.0, g_mCubeCell(2,2) / 1000.0);
		
		double coorda;
		CxDVector3 intersection;
		bool run;
		do {
			run = false;
		
			coorda = DotP(m_vNormal, m_vCenter - py * vecb - pz * vecc) / DotP(m_vNormal, veca);
			
			intersection = coorda * veca + py * vecb + pz * vecc - m_vCenter;
// 			if (verbose) {
// 				mprintf(GREEN, "%26.20g %26.20g %26.20g %26.20g\n", coorda, intersection[0], intersection[1], intersection[2]);
// 			}
			
			a = DotP(intersection, m_vSpan1);
			b = DotP(intersection, m_vSpan2);
// 			if (verbose) {
// 				CxDVector3 v(m_vSpan1);
// 				v *= 1e-5;
// 				v += m_vCenter;
// 				v *= 1000.0;
// 				v = g_mBoxToOrtho * v;
// 				mprintf(GREEN, "Span1: %26.20g %26.20g %26.20g\n", v[0], v[1], v[2]);
// 				v = m_vSpan2;
// 				v *= 1e-5;
// 				v += m_vCenter;
// 				v *= 1000.0;
// 				v = g_mBoxToOrtho * v;
// 				mprintf(GREEN, "Span2: %26.20g %26.20g %26.20g\n", v[0], v[1], v[2]);
// 				v = m_vCenter;
// 				v *= 1000.0;
// 				v = g_mBoxToOrtho * v;
// 				mprintf(GREEN, "Center: %26.20g %26.20g %26.20g\n", v[0], v[1], v[2]);
// 				mprintf(GREEN, "%26.20g %26.20g\n", a, b);
// 			}
			
#ifdef TARGET_WINDOWS
			if (!_finite(a) || !_finite(b))
				return false;
#else
			if (!isfinite(a) || !isfinite(b))
				return false;
#endif
			
			for (z=0;z<m_vaEdges.GetSize();z++) {
				if (m_faEdgesLength[z] < VORI_EPS)
					return false;
				double area = m_vaEdges[z][0]*a + m_vaEdges[z][1]*b + m_vaEdges[z][2];
				double dist = area / m_faEdgesLength[z];
// 				if (verbose)
// 					mprintf(GREEN, "%26.20g\n", dist);
				if (fabs(dist) < VORI_EPS) {
					py += VORI_EPS;
					pz += VORI_EPS;
					run = true;
// 					if (verbose)
// 						mprintf(YELLOW, "run\n");
					break;
				}
				if (area * m_vaEdges[z][2] < 0)
					return false;
			}
		} while (run);
		
		px = coorda;
		return true;
	} else {
		// Hack from 09.10.2014: Ray cannot hit walls which are almost parallel...
		//		if (fabs(m_vNormal[0]) < 1.0E-14)
		//			return false;
		
		ty = py - m_vCenter[1];
		tz = pz - m_vCenter[2];
		
		tf = -lvec1z*lvec2y + lvec1y*lvec2z;
		
		bool run;
		do {
			run = false;
			
			a = (-tz*lvec2y + ty*lvec2z)/tf;
			b = ( tz*lvec1y - ty*lvec1z)/tf;
			
	// 		if (verbose)
	// 			mprintf(GREEN, "%.20g %.20g\n", a, b);
			
#ifdef TARGET_WINDOWS
			if (!_finite(a) || !_finite(b))
				return false;
#else
			if (!isfinite(a) || !isfinite(b))
				return false;
#endif
					
			for (z=0;z<m_vaEdges.GetSize();z++) {
				double area = m_vaEdges[z][0]*a + m_vaEdges[z][1]*b + m_vaEdges[z][2];
				double dist = 2.0 * area / m_faEdgesLength[z];
	// 			if (verbose)
	// 				mprintf(GREEN, "%d %.20g\n", z, dist);
				if (fabs(dist) < VORI_EPS) {
	// 				if (verbose)
	// 					mprintf(GREEN, "repeat %d\n", z);
					ty += VORI_EPS;
					tz += VORI_EPS;
					run = true;
					break;
				}
				if (area * m_vaEdges[z][2] < 0)
	// 			if ((m_vaEdges[z][0]*a + m_vaEdges[z][1]*b + m_vaEdges[z][2])*m_vaEdges[z][2] < 0)
					return false;
			}
		} while (run);
		
		px = a*lvec1x + b*lvec2x + m_vCenter[0];
	// 	if (verbose)
	// 		mprintf(GREEN, "%.20g", m_vCenter[0]);
		
		/*
		* #ifdef TARGET_WINDOWS
		*		if (!_finite(px))
		* #else
		*		if (!isfinite(px))
		* #endif
		*		{
		*			mprintf("\nXRay_Hit(py=%.10f, pz=%.10f) returned %.10f:\n",py,pz,px);
		*			mprintf("  ty=%.20G  tz=%.20G\n",ty,tz);
		*			mprintf("  m_vSpan1: %f  %f  %f\n",m_vSpan1[0],m_vSpan1[1],m_vSpan1[2]);
		*			mprintf("  m_vSpan2: %f  %f  %f\n",m_vSpan2[0],m_vSpan2[1],m_vSpan2[2]);
		*			mprintf("  m_vNormal: %f  %f  %f\n",m_vNormal[0],m_vNormal[1],m_vNormal[2]);
		*			mprintf("  CrossP: %10G\n",CrossP(m_vSpan1,m_vSpan2).GetLength());
		*			mprintf("  m_vCenter: %f  %f  %f\n",m_vCenter[0],m_vCenter[1],m_vCenter[2]);
		*			mprintf("  tf=%.20G,  a=%f,  b=%f.\n",tf,a,b);
		}*/
		
		/*		if (verbose)
		*		{
		*			mprintf("\nXRay_Hit(py=%.10f, pz=%.10f) returned %.10f:\n",py,pz,px);
		*			mprintf("  ty=%.20G  tz=%.20G\n",ty,tz);
		*			mprintf("  m_vSpan1: %f  %f  %f\n",m_vSpan1[0],m_vSpan1[1],m_vSpan1[2]);
		*			mprintf("  m_vSpan2: %f  %f  %f\n",m_vSpan2[0],m_vSpan2[1],m_vSpan2[2]);
		*			mprintf("  m_vNormal: %g  %g  %g\n",m_vNormal[0],m_vNormal[1],m_vNormal[2]);
		*			mprintf("  CrossP: %10G\n",CrossP(m_vSpan1,m_vSpan2).GetLength());
		*			mprintf("  m_vCenter: %f  %f  %f\n",m_vCenter[0],m_vCenter[1],m_vCenter[2]);
		*			mprintf("  tf=%.20G,  a=%f,  b=%f.\n",tf,a,b);
		}*/
		
		return true;
	}
	
}

void printCubeMemFileName(char *name, int length, bool nextFile) {
	if (nextFile)
		g_cubeMemFileIndex += g_cubeMemFileStride;

#ifdef TARGET_WINDOWS
	_snprintf(name, length, "%s%d.cube", g_cubeMemFileNameFixed, g_cubeMemFileIndex);
#elif defined TARGET_LINUX
	snprintf(name, length, "%s%d.cube", g_cubeMemFileNameFixed, g_cubeMemFileIndex);
#else
	sprintf(name, "%s%d.cube", g_cubeMemFileNameFixed, g_cubeMemFileIndex);
#endif
}
