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


#include "plproj.h"
#include "globalvar.h"


CPlProj::CPlProj()
{
	m_p2DF = NULL;
	m_iAverageCounter = 0;
}


CPlProj::~CPlProj()
{
}


void CPlProj::Parse()
{
	int ti, z;
//	char buf[256], buf2[256];
	CxString buf, buf2;
	double tf;

	try { m_p2DF = new C2DF(); } catch(...) { m_p2DF = NULL; }
	if (m_p2DF == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Plane Projection Distribution Function >>>\n\n");

_plprojatoms:
	if (m_bIntra)
	{
		mprintf("    Observing atoms in molecules of kind %s.\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
	} else
	{
		mprintf("    Observing atoms in reference molecule (%s).\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	}
	mprintf("    Which atoms to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ");
	inpprintf("! Which atoms to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n");
	myget(&buf);
	if (strlen(buf) == 0)
	{
		if (!m_oAtoms.ParseAtoms((CMolecule*)g_oaMolecules[m_bIntra?g_iFixMol:m_iShowMol],"#2"))
		{
			eprintf("Strange error ^^\n");
			inpprintf("! Strange error ^\n");
			goto _plprojatoms;
		}
	} else 
	{
		if (!m_oAtoms.ParseAtoms((CMolecule*)g_oaMolecules[m_bIntra?g_iFixMol:m_iShowMol],buf))
		{
			eprintf("Wrong input.\n");
			inpprintf("! Wrong input.\n");
			goto _plprojatoms;
		}
	}

	m_iShowAtomGes = m_oAtoms.m_iAtomGes;

	if (g_bPeriodicX && g_bPeriodicY && g_bPeriodicZ)
		m_fParticleDensity = m_oAtoms.m_iAtomGes * ((CMolecule*)g_oaMolecules[m_bIntra?g_iFixMol:m_iShowMol])->m_laSingleMolIndex.GetSize() / g_fBoxX / g_fBoxY / g_fBoxZ * 1E9f;

	mprintf("\n");
	m_bDrawAtoms = AskYesNo("    Draw atoms of reference molecule in plot (y/n)? [yes] ",true);
	if (m_bDrawAtoms)
	{
_drawagain:
//		sprintf(buf,"    Draw reference atoms (%s%d, %s%d, %s%d) in plot (y/n)? [yes] ",((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1);
		buf.sprintf("    Draw reference atoms (%s%d, %s%d, %s%d) in plot (y/n)? [yes] ",((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1);
		if (AskYesNo(buf,true))
		{
//			sprintf(buf,"%s%d,%s%d,%s%d",((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1);
			buf.sprintf("%s%d,%s%d,%s%d",((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1);
			AskString("    Which additional atoms to draw (e.g. \"C1,C3-5,H\")? [none] ",&buf2,"");
			if (strlen(buf2) != 0)
			{
//				strcat(buf,",");
//				strcat(buf,buf2);
				buf.strcat(",");
				buf.strcat(buf2);
			}
		} else
		{
			AskString("    Which atoms to draw (e.g. \"C1,C3-5,H\", \"*\"=all)? [none] ",&buf,"");
		}
//		mprintf("---> \"%s\".\n",buf);
		if (strlen(buf) != 0)
		{
			if (!m_oDrawAtoms.ParseAtoms((CMolecule*)g_oaMolecules[g_iFixMol],buf))
			{
				eprintf("Wrong input (\"%s\").\n",(const char*)buf);
				goto _drawagain;
			}
		} else
		{
			m_bDrawAtoms = false;
			mprintf("    Not drawing any atoms from reference molecule.\n\n");
		}
	}

	if (m_bDrawAtoms)
	{
		mprintf("    Drawing %d atoms from reference molecule.\n",m_oDrawAtoms.m_iAtomGes);
		m_vaAtomPos.SetSize(m_oDrawAtoms.m_iAtomGes);
		for (z=0;z<m_oDrawAtoms.m_iAtomGes;z++)
			m_vaAtomPos[z] = CxDVector3(0,0,0);
		m_bAverageAtomPos = AskYesNo("    Determine atom draw positions as average (y) or from one snapshot (n)? [yes] ",true);
		mprintf("\n");
	}


_againx:
	m_fMinVal[0] = AskFloat("    Enter lower border in X direction (in pm): [%d.0] ",(float)-HalfBox(),-HalfBox());
	m_fMaxVal[0] = AskFloat("    Enter upper border in X direction (in pm): [%.1f] ",fabs(m_fMinVal[0]),fabs(m_fMinVal[0]));
	if (m_fMaxVal[0] - m_fMinVal[0] <= 0)
	{
		eprintf("    Error: Upper border needs to be larger than lower border.\n");
		inpprintf("! Wrong input.\n");
		goto _againx;
	}

_againy:
	m_fMinVal[1] = AskFloat("    Enter lower border in Y direction (in pm): [%d.0] ",(float)-HalfBox(),-HalfBox());
	m_fMaxVal[1] = AskFloat("    Enter upper border in Y direction (in pm): [%.1f] ",fabs(m_fMinVal[1]),fabs(m_fMinVal[1]));
	if (m_fMaxVal[1] - m_fMinVal[1] <= 0)
	{
		eprintf("    Error: Upper border needs to be larger than lower border.\n");
		inpprintf("! Wrong input.\n");
		goto _againy;
	}

_againslice:
	m_fSliceBorder[0] = AskFloat("    Enter lower border of slice (parallel shift of reference plane, in pm): [-100.0] ",-100.0,-100.0);
	m_fSliceBorder[1] = AskFloat("    Enter upper border of slice (parallel shift of reference plane, in pm): [%.1f] ",fabs(m_fSliceBorder[0]),fabs(m_fSliceBorder[0]));
	if (m_fSliceBorder[1] - m_fSliceBorder[0] <= 0)
	{
		eprintf("    Error: Upper border needs to be larger than lower border.\n");
		inpprintf("! Wrong input.\n");
		goto _againslice;
	}

	tf = (m_fMaxVal[1]-m_fMinVal[1]) / (m_fMaxVal[0]-m_fMinVal[0]);
	mprintf("\n");
	mprintf("    The \"natural\" aspect ratio of this plot would be %.0fpm / %.0fpm = %.4f.\n\n",m_fMaxVal[1]-m_fMinVal[1],m_fMaxVal[0]-m_fMinVal[0],tf);
	m_fAspectRatio = AskFloat("    Please enter aspect ratio (=height/width) of the plot: [%.4f] ",tf,tf);
	m_p2DF->m_fAspectRatio = m_fAspectRatio;
	mprintf("\n");

	ti = g_pDatabase->GetInt("/PLOT2D/DEFAULTS/BIN_RES");
	m_iResolution[0] = AskUnsignedInteger("    Please enter binning resolution in X direction: [%d] ",ti,ti);
	m_iResolution[1] = AskUnsignedInteger("    Please enter binning resolution in Y direction: [%d] ",m_iResolution[0],m_iResolution[0]);

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
	else m_iHistogramRes = 0;

	BuildName();

	mprintf(WHITE,"\n<<< End of Plane Projection Distribution Function <<<\n\n");
}


void CPlProj::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;

//	tmp[0] = 0;
	tmp.sprintf("");

	if (m_bIntra)
//		sprintf(tmp,"plproj_%s_%s%d%s%d%s%d_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,m_oAtoms.m_sName);
		tmp.sprintf("plproj_%s_%s%d%s%d%s%d_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,m_oAtoms.m_sName);
	else
//		sprintf(tmp,"plproj_%s_%s%d%s%d%s%d_%s_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_oAtoms.m_sName);
		tmp.sprintf("plproj_%s_%s%d%s%d%s%d_%s_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_oAtoms.m_sName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CPlProj::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z1t, z1a;
	CxIntArray *a1;

	vec->RemoveAll_KeepSize();
	for (z1t=0;z1t<m_oAtoms.m_baAtomType.GetSize();z1t++)
	{
		a1 = (CxIntArray*)m_oAtoms.m_oaAtoms[z1t];
		for (z1a=0;z1a<a1->GetSize();z1a++)
		{
			if (m_bIntra)
				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[m_oAtoms.m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
			else 
				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[m_oAtoms.m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
		}
	}
	BXOUT;
}




