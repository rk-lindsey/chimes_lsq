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

#include "revsdf.h"
#include "xwordarray.h"
#include "moltools.h"
#include "globalvar.h"


CRevSDF::CRevSDF()
{
	m_fSecondAtomPosX = 0;
	m_fSecondAtomCount = 0;
}


CRevSDF::~CRevSDF()
{
}


void CRevSDF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z1t, z1a;
	CxIntArray *a1;

	vec->RemoveAll_KeepSize();
	for (z1t=0;z1t<m_oAtoms.m_baAtomType.GetSize();z1t++)
	{
//		mprintf("(a) BuildAtomList z1t=%d, WA=%X\n",z1t,m_oAtoms.m_oaAtoms[z1t]);
		a1 = (CxIntArray*)m_oAtoms.m_oaAtoms[z1t];
//		mprintf("(b) BuildAtomList wa.GetSize()=%d\n",a1->GetSize());
		for (z1a=0;z1a<a1->GetSize();z1a++)
			if (m_iRefOrSec)
			{
//				mprintf("(c) BuildAtomList z1t=%d, z1a=%d, Type=%d, WA=%X\n",z1t,z1a,m_oAtoms.m_baAtomType[z1t],obs->m_oaAtomOffset[m_oAtoms.m_baAtomType[z1t]]);
				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[m_oAtoms.m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
			} else vec->Add(((CxIntArray*)ref->m_oaAtomOffset[m_oAtoms.m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
	}
	BXOUT;
}


void CRevSDF::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;

	if (m_iRefOrSec)
//		sprintf(tmp,"%s_%s%d%s%d_%s_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_oAtoms.m_sName);
		tmp.sprintf("%s_%s%d%s%d_%s_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_oAtoms.m_sName);
	else
//		sprintf(tmp,"%s_%s%d%s%d_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,m_oAtoms.m_sName);
		tmp.sprintf("%s_%s%d%s%d_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,m_oAtoms.m_sName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CRevSDF::Parse()
{
	BTIN;
//	char buf[256];
	CxString buf;

	m_iShowAtomGes = 0;

	try { m_p2DF = new C2DF(); } catch(...) { m_p2DF = NULL; }
	if (m_p2DF == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Pseudo SDF >>>\n\n");
_sdfatoms:
	if (m_bIntra)
	{
		mprintf("    Observing Atoms in reference molecule %s.\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		m_iRefOrSec = 0;
	} else
	{
		mprintf("    Observing Atoms in observed molecule %s.\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		m_iRefOrSec = 1;
	}
/*	if (m_iShowMol != -1)
		m_iRefOrSec = AskRangeInteger("    Observe atoms in RM %s (0) or in OM %s (1)? [1] ",0,1,1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			else m_iRefOrSec = 0;*/
	mprintf("    Which atoms to observe (e.g. \"C1,C3-5,H\")? [all] ");
	inpprintf("! Which atoms to observe (e.g. \"C1,C3-5,H\")? [all]\n");
	myget(&buf);

	if (strlen(buf) == 0)
	{
		m_oAtoms.AddAllAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec==0)?g_iFixMol:m_iShowMol],false);
	} else if (!m_oAtoms.ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec==0)?g_iFixMol:m_iShowMol],buf))
	{
		eprintf("Wrong input.\n");
		inpprintf("! Wrong input.\n");
		goto _sdfatoms;
	}

	m_iShowAtomGes += m_oAtoms.m_iAtomGes;
	m_fParticleDensity = m_oAtoms.m_iAtomGes * ((CMolecule*)g_oaMolecules[(m_iRefOrSec==0)?g_iFixMol:m_iShowMol])->m_laSingleMolIndex.GetSize() / g_fBoxX / g_fBoxY / g_fBoxZ * 1E9f;

	m_fRadius = AskFloat("    Please enter radius of this Pseudo SDF in pm: [%.0f.0] ",(float)HalfBox(),(float)HalfBox());

	if (g_bPeriodic && (m_fRadius > HalfBox_Exact()+1.0f))
	{
		eprintf("\nWarning: ");
		mprintf("The specified max. radius is larger than half of the smallest periodic cell vector.\n");
		mprintf("         TRAVIS counts every atom only once (central periodic image).\n");
		mprintf("         Expect the analysis to decay to zero for large radii.\n\n");
		AskYesNo("         Acknowledged [yes] ",true);
		mprintf("\n");
	}

	m_iResolution = AskUnsignedInteger("    Please enter binning resolution of this Pseudo SDF per dimension: [100] ",100)+1;
	m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
	m_bMirrorY = AskYesNo("    Force this Pseudo SDF to be mirror-symmetrical to the X axis (y/n)? [no] ",false);
	if (m_bMirrorY)
		m_bMirrorBond = AskYesNo("    Mirror this Pseudo SDF at middle of bond (y) or ref. atom (n)? [yes] ",true);
	m_bCorrectRadial = AskYesNo("    Correct radial distribition (y/n)? [yes] ",true);
	m_bCorrectAngle = AskYesNo("    Apply cone correction (y/n)? [yes] ",true);
	m_bDrawAtoms = AskYesNo("    Show reference atoms in Pseudo SDF plot (y/n)? [yes] ",true);
	m_bCreateRevSDF = AskYesNo("    Create a volumetric Revolution SDF (y/n)? [no] ",false);
	if (m_bCreateRevSDF)
	{
		g_bCreateRevSDF = true;
		m_iRevSDFRes = AskUnsignedInteger("    Enter binning resolution of Revolution SDF per dimension: [100] ",100);
	}
	BuildName();
	mprintf(WHITE,"\n<<< End of Pseudo SDF <<<\n\n");
	BTOUT;
}

void CRevSDF::CreateRevSDF()
{
	int x, y, z;
	CxVector3 v;

	try { m_pRevSDF = new C3DF<double>(); } catch(...) { m_pRevSDF = NULL; }
	if (m_pRevSDF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pRevSDF->m_fMinVal[0] = -m_fRadius;
	m_pRevSDF->m_fMaxVal[0] = m_fRadius;
	m_pRevSDF->m_fMinVal[1] = -m_fRadius;
	m_pRevSDF->m_fMaxVal[1] = m_fRadius;
	m_pRevSDF->m_fMinVal[2] = -m_fRadius;
	m_pRevSDF->m_fMaxVal[2] = m_fRadius;
	m_pRevSDF->m_iRes[0] = m_iRevSDFRes;
	m_pRevSDF->m_iRes[1] = m_iRevSDFRes;
	m_pRevSDF->m_iRes[2] = m_iRevSDFRes;
	m_pRevSDF->Create();

	for (z=0;z<m_pRevSDF->m_iRes[2];z++)
	{
		v[2] = (float)(((double)(z+0.5)/m_pRevSDF->m_iRes[2]*(m_pRevSDF->m_fMaxVal[2]-m_pRevSDF->m_fMinVal[2])) + m_pRevSDF->m_fMinVal[2]);
		for (y=0;y<m_pRevSDF->m_iRes[1];y++)
		{
			v[1] = (float)(((double)(y+0.5)/m_pRevSDF->m_iRes[1]*(m_pRevSDF->m_fMaxVal[1]-m_pRevSDF->m_fMinVal[1])) + m_pRevSDF->m_fMinVal[1]);
			for (x=0;x<m_pRevSDF->m_iRes[0];x++)
			{
				v[0] = (float)(((double)(x+0.5)/m_pRevSDF->m_iRes[0]*(m_pRevSDF->m_fMaxVal[0]-m_pRevSDF->m_fMinVal[0])) + m_pRevSDF->m_fMinVal[0]);
				m_pRevSDF->m_pBin[z*m_pRevSDF->m_iResXY + y*m_pRevSDF->m_iRes[0] + x] = m_p2DF->GetValue(sqrt(v[0]*v[0]+v[2]*v[2]),v[1]);
			}
		}
	}
}

