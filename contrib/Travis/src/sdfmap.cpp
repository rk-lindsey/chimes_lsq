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


#include "sdfmap.h"
#include "globalvar.h"


CSDFMap::CSDFMap()
{
	m_pValueMin = NULL;
	m_pValueMax = NULL;
	m_pValueAvg = NULL;
	m_pCount = NULL;
	m_sName = NULL;
	m_sSubName = NULL;
}


CSDFMap::~CSDFMap()
{
	if (m_pValueMin != NULL)
	{
		delete m_pValueMin;
		m_pValueMin = NULL;
	}

	if (m_pValueMax != NULL)
	{
		delete m_pValueMax;
		m_pValueMax = NULL;
	}

	if (m_pValueAvg != NULL)
	{
		delete m_pValueAvg;
		m_pValueAvg = NULL;
	}

	if (m_pCount != NULL)
	{
		delete m_pCount;
		m_pCount = NULL;
	}

	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}

	if (m_sSubName != NULL)
	{
		delete[] m_sSubName;
		m_sSubName = NULL;
	}
}


void CSDFMap::Parse()
{
	int ti, z;
//	char buf[256], buf2[256];
	CxString buf, buf2;
	CAtomGroup *ag;

	mprintf(WHITE,"\n>>> SDF Surface Map >>>\n\n");

	mprintf("    List of scalar quantities which can be mapped to SDF surface:\n");
	mprintf("      1.) Absolute velocity\n");
	mprintf("      2.) Absolute force\n");

	mprintf("\n");

_quantityagain:
	m_iQuantity = AskRangeInteger_ND("    Select quantity to map onto SDF surface (1-2): ",1,2);

	switch(m_iQuantity)
	{
		/********************************************************************************************************/
		case 1:
			g_bUseVelocities = true;
			if (g_fTimestepLength == 0)
			{
				mprintf("\n");
				g_fTimestepLength = AskFloat("    Enter the length of one trajectory time step in fs: [0.5] ",0.5f);
			}

			_vnextmol:
			mprintf("\n");
//			sprintf(buf,"    Evaluate velocity of atoms from which molecule type (");
			buf.sprintf("    Evaluate velocity of atoms from which molecule type (");
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
//				sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//				strcat(buf,buf2);
//				if (z < g_oaMolecules.GetSize()-1)
//					strcat(buf,", ");
				buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
				buf.strcat(buf2);
				if (z < g_oaMolecules.GetSize()-1)
					buf.strcat(", ");
			}
//			strcat(buf,")? ");
			buf.strcat(")? ");

			ti = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize()) - 1;

			for (z=0;z<m_oaAtomGroups.GetSize();z++)
			{
				if (((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_iIndex == ti)
				{
					eprintf("This molecule type is already in use.\n");
					goto _vnextmol;
				}
			}

			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			_vagain:
			AskString("    Which atoms to take into account from %s (*=all)? [*] ",&buf,"*",((CMolecule*)g_oaMolecules[ti])->m_sName);
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[ti],buf))
				goto _vagain;

			m_oaAtomGroups.Add(ag);

			mprintf("\n");
			if (AskYesNo("    Evaluate velocity of atoms from another molecule type (y/n)? [no] ",false))
				goto _vnextmol;

//			sprintf(buf,"velocity");
			buf.sprintf("velocity");
			for (z=0;z<m_oaAtomGroups.GetSize();z++)
			{
//				sprintf(buf2,"_%s_%s",((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_sName,((CAtomGroup*)m_oaAtomGroups[z])->m_sName);
//				strcat(buf,buf2);
				buf2.sprintf("_%s_%s",((CAtomGroup*)m_oaAtomGroups[z])->m_pMolecule->m_sName,((CAtomGroup*)m_oaAtomGroups[z])->m_sName);
				buf.strcat(buf2);
			}

			try { m_sSubName = new char[strlen(buf)+1]; } catch(...) { m_sSubName = NULL; }
			if (m_sSubName == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(m_sSubName,buf);
			break;

		/********************************************************************************************************/
		default:
			eprintf("Not yet implemented.\n");
			goto _quantityagain;
	}

	mprintf("\n");

	m_bModeAvg = AskYesNo("    Create surface map of average values (y/n)? [yes] ",true);
	m_bModeMin = AskYesNo("    Create surface map of minimal values (y/n)? [no] ",false);
	m_bModeMax = AskYesNo("    Create surface map of maximal values (y/n)? [no] ",false);

	mprintf("\n");

	m_bSphere = AskYesNo("    Add each value to a sphere of more than one grid points (y/n)? [no] ",false);

	if (m_bSphere)
		m_iSphereRadius = AskUnsignedInteger("    Enter the radius of this sphere in grid points: [3] ",3);

	mprintf("\n");

	m_fRadius = AskFloat("    Please enter radius of this SDF surface map in pm: [%d.0] ",(float)HalfBox(),HalfBox());

_sdfresagain:
	ti = g_pDatabase->GetInt("/PLOT3D/DEFAULTS/BIN_RES");
	m_iResolution = AskUnsignedInteger("    Please enter binning resolution of this SDF surface map per dimension: [%d] ",ti,ti);

	if (m_iResolution > 300)
	{
		eprintf("\nWarning: ");
		mprintf("This large resolution will consume much RAM (%s)\n",FormatBytes(pow((double)m_iResolution,3)*sizeof(double)));
		mprintf("         and result in very large SDF map output files.\n\n");
		if (!AskYesNo("    Use this resolution (y/n)? [no] ",false))
			goto _sdfresagain;
	}

	BuildName();

	mprintf(WHITE,"\n<<< End of SDF Surface Map <<<\n\n");
}


void CSDFMap::Create()
{
	m_pCount = new C3DF<double>();
	m_pCount->m_fMinVal[0] = -m_fRadius;
	m_pCount->m_fMaxVal[0] = m_fRadius;
	m_pCount->m_fMinVal[1] = -m_fRadius;
	m_pCount->m_fMaxVal[1] = m_fRadius;
	m_pCount->m_fMinVal[2] = -m_fRadius;
	m_pCount->m_fMaxVal[2] = m_fRadius;
	m_pCount->m_iRes[0] = m_iResolution;
	m_pCount->m_iRes[1] = m_iResolution;
	m_pCount->m_iRes[2] = m_iResolution;
	m_pCount->m_iHistogramRes = 0;
	m_pCount->Create();

	if (m_bModeAvg)
	{
		m_pValueAvg = new C3DF<double>();
		m_pValueAvg->m_fMinVal[0] = -m_fRadius;
		m_pValueAvg->m_fMaxVal[0] = m_fRadius;
		m_pValueAvg->m_fMinVal[1] = -m_fRadius;
		m_pValueAvg->m_fMaxVal[1] = m_fRadius;
		m_pValueAvg->m_fMinVal[2] = -m_fRadius;
		m_pValueAvg->m_fMaxVal[2] = m_fRadius;
		m_pValueAvg->m_iRes[0] = m_iResolution;
		m_pValueAvg->m_iRes[1] = m_iResolution;
		m_pValueAvg->m_iRes[2] = m_iResolution;
		m_pValueAvg->m_iHistogramRes = 0;
		m_pValueAvg->Create();
	}

	if (m_bModeMin)
	{
		m_pValueMin = new C3DF<double>();
		m_pValueMin->m_fMinVal[0] = -m_fRadius;
		m_pValueMin->m_fMaxVal[0] = m_fRadius;
		m_pValueMin->m_fMinVal[1] = -m_fRadius;
		m_pValueMin->m_fMaxVal[1] = m_fRadius;
		m_pValueMin->m_fMinVal[2] = -m_fRadius;
		m_pValueMin->m_fMaxVal[2] = m_fRadius;
		m_pValueMin->m_iRes[0] = m_iResolution;
		m_pValueMin->m_iRes[1] = m_iResolution;
		m_pValueMin->m_iRes[2] = m_iResolution;
		m_pValueMin->m_iHistogramRes = 0;
		m_pValueMin->Create();
	}

	if (m_bModeMax)
	{
		m_pValueMax = new C3DF<double>();
		m_pValueMax->m_fMinVal[0] = -m_fRadius;
		m_pValueMax->m_fMaxVal[0] = m_fRadius;
		m_pValueMax->m_fMinVal[1] = -m_fRadius;
		m_pValueMax->m_fMaxVal[1] = m_fRadius;
		m_pValueMax->m_fMinVal[2] = -m_fRadius;
		m_pValueMax->m_fMaxVal[2] = m_fRadius;
		m_pValueMax->m_iRes[0] = m_iResolution;
		m_pValueMax->m_iRes[1] = m_iResolution;
		m_pValueMax->m_iRes[2] = m_iResolution;
		m_pValueMax->m_iHistogramRes = 0;
		m_pValueMax->Create();
	}
}


void CSDFMap::Process(int rm, CTimeStep *ts)
{
	int z, z2, z3, z4, ti;
	double px, py, pz, v;
	CAtomGroup *ag;
	CMolecule *m;
	CSingleMolecule *sm;

	(void)rm; // Keep this: probably other quantities will require knowledge of current reference molecule

	switch(m_iQuantity)
	{
		case 1: // Absolute Velocity
			for (z=0;z<m_oaAtomGroups.GetSize();z++)
			{
				ag = (CAtomGroup*)m_oaAtomGroups[z];
				m = (CMolecule*)g_oaMolecules[ag->m_pMolecule->m_iIndex];
				for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
				{
					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
					for (z3=0;z3<ag->m_baAtomType.GetSize();z3++)
					{
						for (z4=0;z4<((CxIntArray*)ag->m_oaAtoms[z3])->GetSize();z4++)
						{
							ti = ((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z3]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z3])->GetAt(z4));

							px = ts->m_vaCoords[ti][0];
							py = ts->m_vaCoords[ti][1];
							pz = ts->m_vaCoords[ti][2];

							v = ts->m_vaVelocities[ti].GetLength();

							AddValue(px,py,pz,v);
						}
					}
				}
			}
			break;

		default:
			return;
	}
}


void CSDFMap::BuildName()
{
//	char tmp[32768];
	CxString tmp;

	if (m_bSphere)
//		sprintf(tmp,"sdfmap_%s_%s%d%s%d%s%d_%s_sph%d",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,m_sSubName,m_iSphereRadius);
		tmp.sprintf("sdfmap_%s_%s%d%s%d%s%d_%s_sph%d",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,m_sSubName,m_iSphereRadius);
	else
//		sprintf(tmp,"sdfmap_%s_%s%d%s%d%s%d_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,m_sSubName);
		tmp.sprintf("sdfmap_%s_%s%d%s%d%s%d_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,m_sSubName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
}


void CSDFMap::Finish()
{
	int z;
	double mimi, mima, miavg, mico, mami, mama, maavg, maco, avmi, avma, avavg, avco, toco, ze;
	C3DF<double> *tempSDF;
//	char buf[256];
	CxString buf;

	mprintf(YELLOW,"    * Finishing %s ...\n",m_sName);

	mprintf("        Calculating statistics...\n");

	toco = 0;
	ze = 0;

	if (m_bModeMin)
	{
		mimi = 1.0E30;
		mima = -1.0E30;
		miavg = 0;
		mico = 0;
	}

	if (m_bModeMax)
	{
		mami = 1.0E30;
		mama = -1.0E30;
		maavg = 0;
		maco = 0;
	}

	if (m_bModeAvg)
	{
		avmi = 1.0E30;
		avma = -1.0E30;
		avavg = 0;
		avco = 0;
	}

	for (z=0;z<m_iResolution*m_iResolution*m_iResolution;z++)
	{
		if (m_pCount->m_pBin[z] != 0)
		{
			toco += m_pCount->m_pBin[z];

			if (m_bModeAvg)
			{
				m_pValueAvg->m_pBin[z] /= m_pCount->m_pBin[z];
				avco++;
				avavg += m_pValueAvg->m_pBin[z];
				if (m_pValueAvg->m_pBin[z] > avma)
					avma = m_pValueAvg->m_pBin[z];
				if (m_pValueAvg->m_pBin[z] < avmi)
					avmi = m_pValueAvg->m_pBin[z];
			}

			if (m_bModeMin)
			{
				mico++;
				miavg += m_pValueMin->m_pBin[z];
				if (m_pValueMin->m_pBin[z] > mima)
					mima = m_pValueMin->m_pBin[z];
				if (m_pValueMin->m_pBin[z] < mimi)
					mimi = m_pValueMin->m_pBin[z];
			}

			if (m_bModeMax)
			{
				maco++;
				maavg += m_pValueMax->m_pBin[z];
				if (m_pValueMax->m_pBin[z] > mama)
					mama = m_pValueMax->m_pBin[z];
				if (m_pValueMax->m_pBin[z] < mami)
					mami = m_pValueMax->m_pBin[z];
			}

		} else ze++;
	}

	mprintf("        %.0f bin entries in total.\n",toco);
	mprintf("        %.0f bin values are zero (%.3f%c).\n",ze,ze/m_iResolution/m_iResolution/m_iResolution*100.0,'%');

	if (m_bModeAvg)
	{
		mprintf(WHITE,"\n        Average values:\n");
		mprintf("          Min. value %10.5f\n",avmi);
		mprintf("          Max. value %10.5f\n",avma);
		mprintf("          Avg. value %10.5f  (without all zeros %10.5f)\n",avavg/m_iResolution/m_iResolution/m_iResolution,avavg/avco);
	}

	if (m_bModeMin)
	{
		mprintf(WHITE,"\n        Minimal values:\n");
		mprintf("          Min. value %10.5f\n",mimi);
		mprintf("          Max. value %10.5f\n",mima);
		mprintf("          Avg. value %10.5f  (without all zeros %10.5f)\n",miavg/m_iResolution/m_iResolution/m_iResolution,miavg/mico);
	}

	if (m_bModeMax)
	{
		mprintf(WHITE,"\n        Maximal values:\n");
		mprintf("          Min. value %10.5f\n",mami);
		mprintf("          Max. value %10.5f\n",mama);
		mprintf("          Avg. value %10.5f  (without all zeros %10.5f)\n",maavg/m_iResolution/m_iResolution/m_iResolution,maavg/maco);
	}

	if (m_bModeAvg)
	{
		mprintf(WHITE,"\n        Writing average SDF map:\n");
		for (z=0;z<=g_iSDFMapSmoothGrade;z++)
		{
			try { tempSDF = new C3DF<double>(); } catch(...) { tempSDF = NULL; }
			if (tempSDF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			tempSDF->CopyFrom(m_pValueAvg);

			if (z != 0)
				tempSDF->Smooth(z);

			tempSDF->CalcMaxEntry();
			mprintf("        Data range from %.6f to %.6f.\n",tempSDF->m_fMinEntry,tempSDF->m_fMaxEntry);

			if (g_pDatabase->GetBool("/PLOT3D/FORMATS/WRITE_PLT"))
			{
				if (z != 0)
//					sprintf(buf,"_avg.s%d.plt",z);
					buf.sprintf("_avg.s%d.plt",z);
				else
//					sprintf(buf,"_avg.plt");
					buf.sprintf("_avg.plt");

				mprintf("        Saving SDF map as \"%s%s\"...\n",m_sName,(const char*)buf);
				tempSDF->WritePLT("",m_sName,buf,true);
			}

			if (g_pDatabase->GetBool("/PLOT3D/FORMATS/WRITE_CUBE"))
			{
				if (z != 0)
//					sprintf(buf,"_avg.s%d.cube",z);
					buf.sprintf("_avg.s%d.cube",z);
				else
//					sprintf(buf,"_avg.cube");
					buf.sprintf("_avg.cube");

				mprintf("        Saving SDF map as \"%s%s\"...\n",m_sName,(const char*)buf);
				tempSDF->WriteCube("",m_sName,buf,true);
			}

			delete tempSDF;
		}
	}

	if (m_bModeMin)
	{
		mprintf(WHITE,"\n        Writing minimal SDF map:\n");
		for (z=0;z<=g_iSDFMapSmoothGrade;z++)
		{
			try { tempSDF = new C3DF<double>(); } catch(...) { tempSDF = NULL; }
			if (tempSDF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			tempSDF->CopyFrom(m_pValueMin);

			if (z != 0)
				tempSDF->Smooth(z);

			tempSDF->CalcMaxEntry();
			mprintf("        Data range from %.6f to %.6f.\n",tempSDF->m_fMinEntry,tempSDF->m_fMaxEntry);

			if (g_pDatabase->GetBool("/PLOT3D/FORMATS/WRITE_PLT"))
			{
				if (z != 0)
//					sprintf(buf,"_min.s%d.plt",z);
					buf.sprintf("_min.s%d.plt",z);
				else
//					sprintf(buf,"_min.plt");
					buf.sprintf("_min.plt");

				mprintf("        Saving SDF map as \"%s%s\"...\n",m_sName,(const char*)buf);
				tempSDF->WritePLT("",m_sName,buf,true);
			}

			if (g_pDatabase->GetBool("/PLOT3D/FORMATS/WRITE_CUBE"))
			{
				if (z != 0)
//					sprintf(buf,"_min.s%d.cube",z);
					buf.sprintf("_min.s%d.cube",z);
				else
//					sprintf(buf,"_min.cube");
					buf.sprintf("_min.cube");

				mprintf("        Saving SDF map as \"%s%s\"...\n",m_sName,(const char*)buf);
				tempSDF->WriteCube("",m_sName,buf,true);
			}

			delete tempSDF;
		}
	}

	if (m_bModeMax)
	{
		mprintf(WHITE,"\n        Writing maximal SDF map:\n");
		for (z=0;z<=g_iSDFMapSmoothGrade;z++)
		{
			try { tempSDF = new C3DF<double>(); } catch(...) { tempSDF = NULL; }
			if (tempSDF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			tempSDF->CopyFrom(m_pValueMax);

			if (z != 0)
				tempSDF->Smooth(z);

			tempSDF->CalcMaxEntry();
			mprintf("        Data range from %.6f to %.6f.\n",tempSDF->m_fMinEntry,tempSDF->m_fMaxEntry);

			if (g_pDatabase->GetBool("/PLOT3D/FORMATS/WRITE_PLT"))
			{
				if (z != 0)
//					sprintf(buf,"_max.s%d.plt",z);
					buf.sprintf("_max.s%d.plt",z);
				else
//					sprintf(buf,"_max.plt");
					buf.sprintf("_max.plt");

				mprintf("        Saving SDF map as \"%s%s\"...\n",m_sName,(const char*)buf);
				tempSDF->WritePLT("",m_sName,buf,true);
			}

			if (g_pDatabase->GetBool("/PLOT3D/FORMATS/WRITE_CUBE"))
			{
				if (z != 0)
//					sprintf(buf,"_max.s%d.cube",z);
					buf.sprintf("_max.s%d.cube",z);
				else
//					sprintf(buf,"_max.cube");
					buf.sprintf("_max.cube");

				mprintf("        Saving SDF map as \"%s%s\"...\n",m_sName,(const char*)buf);
				tempSDF->WriteCube("",m_sName,buf,true);
			}

			delete tempSDF;
		}
	}
}

