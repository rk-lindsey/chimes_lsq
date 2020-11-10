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

#include "dacf.h"
#include "globalvar.h"


CDACFSub::CDACFSub()
{
	m_bDistTrace = false;
	m_bNewMode = false;
	m_fEqCounter = 0;
}


CDACFSub::~CDACFSub()
{
}


CDACF::CDACF()
{
	m_bLifetimeSpectrum = false;
	m_oaSubDACFs.SetName("CDACF::m_oaSubDACFs");
	m_faWeight1.SetName("CDACF::m_faWeight1");
	m_faWeight2.SetName("CDACF::m_faWeight2");
	m_oaLTSpectra.SetName("CDACF::m_oaLTSpectra");
}


CDACF::~CDACF()
{
}


void CDACF::Parse()
{
	int z, z2, ti;
	float tf, tf2;
	CDACFSub *ds;
	CxFloatArray fa;
	CxIntArray wa;
//	char buf[256], buf2[256];
	CxString buf, buf2;

	if (g_oaMolecules.GetSize() > 1)
	{
//		sprintf(buf,"    Enter first molecule for aggregation (");
		buf.sprintf("    Enter first molecule for aggregation (");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
//			sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//			strcat(buf,buf2);
//			if (z < g_oaMolecules.GetSize()-1)
//				strcat(buf,", ");
			buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
			buf.strcat(buf2);
			if (z < g_oaMolecules.GetSize()-1)
				buf.strcat(", ");
		}
//		strcat(buf,"): ");
		buf.strcat("): ");
		m_iFirstMol = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize()) - 1;

//		sprintf(buf,"    Enter second molecule for aggregation (");
		buf.sprintf("    Enter second molecule for aggregation (");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
//			sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//			strcat(buf,buf2);
//			if (z < g_oaMolecules.GetSize()-1)
//				strcat(buf,", ");
			buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
			buf.strcat(buf2);
			if (z < g_oaMolecules.GetSize()-1)
				buf.strcat(", ");
		}
//		strcat(buf,"): ");
		buf.strcat("): ");
		m_iSecondMol = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize()) - 1;
	} else
	{
		m_iFirstMol = 0;
		m_iSecondMol = 0;
	}
	mprintf("    Observing aggregates between %s and %s.\n\n",((CMolecule*)g_oaMolecules[m_iFirstMol])->m_sName,((CMolecule*)g_oaMolecules[m_iSecondMol])->m_sName);

//	sprintf(buf,"%s_%s",((CMolecule*)g_oaMolecules[m_iFirstMol])->m_sName,((CMolecule*)g_oaMolecules[m_iSecondMol])->m_sName);
	buf.sprintf("%s_%s",((CMolecule*)g_oaMolecules[m_iFirstMol])->m_sName,((CMolecule*)g_oaMolecules[m_iSecondMol])->m_sName);

	try { m_sName = new char[strlen(buf)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,buf);

	if (g_bDACF)
	{
		if (g_iTrajSteps != -1)
			m_iDACFRes = AskUnsignedInteger("    Enter the depth (resolution) of the DACF in time steps: [%d] ",int(g_iTrajSteps*0.75),int(g_iTrajSteps*0.75));
				else m_iDACFRes = AskUnsignedInteger("    Enter the depth (resolution) of the DACF in time steps: [5000] ",5000);
		m_bFitDACF = AskYesNo("    Fit poly-exponential functions to DACF and integrate over [0,infinity) (y/n)? [yes] ",true);
		if (m_bFitDACF)
		{
			mprintf("\n    TRAVIS will fit functions with different count of exponential terms to your data.\n");
			mprintf("    This enables you to see how many different processes contribute to your decay.\n\n");
			m_iFitDegreeMin = AskUnsignedInteger("    How many exponential terms to use at least: [1] ",1);
			m_iFitDegreeMax = AskUnsignedInteger("    How many exponential terms to use at most: [4] ",4);
			mprintf("\n");
		}
			m_bLifetimeSpectrum = false;
	}

	if (g_bNbExchange)
	{
		if (g_iTrajSteps != -1)
			m_pNbExchange->m_iTimeDepth = AskUnsignedInteger("    Enter the depth (resolution) of the neighborhood exchange functions in time steps: [%d] ",int(g_iTrajSteps*0.75),int(g_iTrajSteps*0.75));
				else m_pNbExchange->m_iTimeDepth = AskUnsignedInteger("    Enter the depth (resolution) of the neighborhood exchange functions in time steps: [5000] ",5000);
	}

	if (g_bDDisp)
	{
		g_bKeepUnfoldedCoords = true;

		try { m_pCenterAtoms1 = new CAtomGroup(); } catch(...) { m_pCenterAtoms1 = NULL; }
		if (m_pCenterAtoms1 == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		try { m_pCenterAtoms2 = new CAtomGroup(); } catch(...) { m_pCenterAtoms2 = NULL; }
		if (m_pCenterAtoms2 == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		mprintf(YELLOW,"\n*** Dimer Displacement\n\n");
		mprintf("    The \"dimer center\" is the weightened average of 2 atoms from the monomers.\n\n");
_ca1:	mprintf("    Enter atoms from 1st molecule [center of mass]: ");
		inpprintf("! Enter atoms from 1st molecule [center of mass]:\n");
		myget(&buf);
		if (strlen(buf)!=0)
		{
			if (!m_pCenterAtoms1->ParseAtoms((CMolecule*)g_oaMolecules[m_iFirstMol],buf))
			{
				eprintf("Wrong input.\n");
				inpprintf("! Wrong input.\n");
				goto _ca1;
			}
		} else
		{
			if (!m_pCenterAtoms1->ParseAtoms((CMolecule*)g_oaMolecules[m_iFirstMol],"#2"))
			{
				eprintf("Weird error.\n");
				inpprintf("! Weird error.\n");
				abort();
			}
		}
		if (m_pCenterAtoms1->m_iAtomGes > 1)
		{
_wm1again:	ti = AskRangeInteger("    Weight Atoms in 1st mol. uniform (1), by mass (2) or manually (3)? [2]",1,3,2);
			switch(ti)
			{
				case 1:
					for (z=0;z<m_pCenterAtoms1->m_iAtomGes;z++)
						m_faWeight1.Add(1.0f);
					break;
				case 2:
					for (z=0;z<m_pCenterAtoms1->m_baAtomType.GetSize();z++)
						for (z2=0;z2<((CxIntArray*)m_pCenterAtoms1->m_oaAtoms[z])->GetSize();z2++)
							m_faWeight1.Add(((CAtom*)g_oaAtoms[m_pCenterAtoms1->m_baRealAtomType[z]])->m_pElement->m_fMass);
					break;
				case 3:
					for (z=0;z<m_pCenterAtoms1->m_baAtomType.GetSize();z++)
						for (z2=0;z2<((CxIntArray*)m_pCenterAtoms1->m_oaAtoms[z])->GetSize();z2++)
						{
							tf = AskFloat("    Enter weight for %s%d: [1.0] ",1.0f,((CAtom*)g_oaAtoms[m_pCenterAtoms1->m_baRealAtomType[z]])->m_sName,((CxIntArray*)m_pCenterAtoms1->m_oaAtoms[z])->GetAt(z2)+1);
							m_faWeight1.Add(tf);
						}
					break;
				default:
					goto _wm1again;
			}
			tf = 0;
			for (z=0;z<m_faWeight1.GetSize();z++)
				tf += m_faWeight1[z];
			for (z=0;z<m_faWeight1.GetSize();z++)
				m_faWeight1[z] /= tf;
		} else m_faWeight1.Add(1.0f);
_ca2:	mprintf("    Enter atoms from 2nd molecule [center of mass]: ");
		inpprintf("! Enter atoms from 2nd molecule [center of mass]:\n");
		myget(&buf);
		if (strlen(buf)!=0)
		{
			if (!m_pCenterAtoms2->ParseAtoms((CMolecule*)g_oaMolecules[m_iSecondMol],buf))
			{
				eprintf("Wrong input.\n");
				inpprintf("! Wrong input.\n");
				goto _ca2;
			}
		} else
		{
			if (!m_pCenterAtoms2->ParseAtoms((CMolecule*)g_oaMolecules[m_iSecondMol],"#2"))
			{
				eprintf("Weird error.\n");
				inpprintf("! Weird error.\n");
				abort();
			}
		}
		if (m_pCenterAtoms2->m_iAtomGes > 1)
		{
_wm2again:	ti = AskRangeInteger("    Weight atoms in 2nd mol. uniform (1), by mass (2) or manually (3)? [2]",1,3,2);
			switch(ti)
			{
				case 1:
					for (z=0;z<m_pCenterAtoms2->m_iAtomGes;z++)
						m_faWeight2.Add(1.0f);
					break;
				case 2:
					for (z=0;z<m_pCenterAtoms2->m_baAtomType.GetSize();z++)
						for (z2=0;z2<((CxIntArray*)m_pCenterAtoms2->m_oaAtoms[z])->GetSize();z2++)
							m_faWeight2.Add(((CAtom*)g_oaAtoms[m_pCenterAtoms2->m_baRealAtomType[z]])->m_pElement->m_fMass);
					break;
				case 3:
					for (z=0;z<m_pCenterAtoms2->m_baAtomType.GetSize();z++)
						for (z2=0;z2<((CxIntArray*)m_pCenterAtoms2->m_oaAtoms[z])->GetSize();z2++)
						{
							tf = AskFloat("    Enter weight for %s%d: [1.0] ",1.0f,((CAtom*)g_oaAtoms[m_pCenterAtoms2->m_baRealAtomType[z]])->m_sName,((CxIntArray*)m_pCenterAtoms2->m_oaAtoms[z])->GetAt(z2)+1);
							m_faWeight2.Add(tf);
						}
					break;
				default:
					goto _wm2again;
			}
			tf = 0;
			for (z=0;z<m_faWeight2.GetSize();z++)
				tf += m_faWeight2[z];
			for (z=0;z<m_faWeight2.GetSize();z++)
				m_faWeight2[z] /= tf;
		} else m_faWeight2.Add(1.0f);
_mwagain:	
		ti = AskRangeInteger("    Weight 1st to 2nd mol. uniform (1), by mass (2) or manually (3)? [2] ",1,3,2);
		switch(ti)
		{
			case 1:
				m_fWeightMol1 = 0.5f;
				break;
			case 2:
				tf = ((CMolecule*)g_oaMolecules[m_iFirstMol])->m_fMass;
				tf2 = ((CMolecule*)g_oaMolecules[m_iSecondMol])->m_fMass;
				m_fWeightMol1 = tf / (tf+tf2);
				break;
			case 3:
				tf = AskFloat("    Enter weight for 1st molecule: [1.0] ",1.0f);
				tf2 = AskFloat("    Enter weight for 2nd molecule: [1.0] ",1.0f);
				m_fWeightMol1 = tf / (tf+tf2);
				break;
			default:
				goto _mwagain;
		}
		mprintf("    Weight Mol.1 : Mol.2 is %.3f : %.3f.\n\n",m_fWeightMol1,1.0f-m_fWeightMol1);
		m_fLargestDisplacement = AskFloat("    Enter largest displacement in pm [1000.0]: ",1000.0f);
		m_iDisplacementRes = AskUnsignedInteger("    Enter displacement resolution [100]: ",100);
		m_fMaxVel = AskFloat("    Enter maximum dimer center velocity in pm/ps (0 means no max. vel.): [10000.0] ",10000.0f);
		if (m_fMaxVel != 0)
			m_bRemoveMaxVel = AskYesNo("    If Velocity is above max. value: Remove value (y) or only warn (n)? [no] ",false);
	} else m_fMaxVel = 0.0f;

	if (g_bDLDF)
	{
		if (g_bDACF)
		{
			m_fLargestLifetime = AskFloat("    Enter largest lifetime in ps [%.3f]: ",m_iDACFRes*g_fTimestepLength/1000.0,m_iDACFRes*g_fTimestepLength/1000.0);
			m_iLifetimeRes = AskUnsignedInteger("    Enter lifetime binning resolution [100]: ",100);
		} else
		{
			if (g_iTrajSteps != -1)
				m_fLargestLifetime = AskFloat("    Enter largest lifetime in ps [%.3f]: ",g_iTrajSteps*0.75*g_fTimestepLength/1000.0,g_iTrajSteps*0.75*g_fTimestepLength/1000.0);
					else m_fLargestLifetime = AskFloat("    Enter largest lifetime in ps [5.0]: ",5.0f);
			m_iLifetimeRes = AskUnsignedInteger("    Enter lifetime binning resolution [100]: ",int(m_fLargestLifetime*1000.0/g_fTimestepLength),int(m_fLargestLifetime*1000.0/g_fTimestepLength));
		}

		if (m_fLargestLifetime * 1000.0 / g_fTimestepLength < m_iLifetimeRes)
		{
			mprintf(WHITE,"\nWARNING: Largest lifetime of %.3f ps is %d timesteps - less than lifetime resolution %d!\n",m_fLargestLifetime,(int)(m_fLargestLifetime * 1000.0 / g_fTimestepLength),m_iLifetimeRes);
			mprintf(WHITE,"         Changing lifetime resolution down to %d.\n",(int)(m_fLargestLifetime * 1000.0 / g_fTimestepLength));
			m_iLifetimeRes = (int)(m_fLargestLifetime * 1000.0 / g_fTimestepLength);
		}
	}

/*	if (g_bDDisp && g_bDLDF)
		g_bDLDisp = AskYesNo("\n    Create an combined Dimer Lifetime/Displacement-Function (y/n)? [yes] ",true);
			else g_bDLDisp = false;*/

	mprintf("\n");

		m_bDACFGrid = false;

	{
		mprintf("\nNow you have to enter the condition(s) for aggregation.\n\n");

		try { m_pCondition = new CConditionGroup(); } catch(...) { m_pCondition = NULL; }
		if (m_pCondition == NULL) NewException((double)sizeof(CConditionGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pCondition->Parse(m_iFirstMol,m_iSecondMol);

		try { ds = new CDACFSub(); } catch(...) { ds = NULL; }
		if (ds == NULL) NewException((double)sizeof(CDACFSub),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		ds->m_iRefMol = m_iFirstMol;
		ds->m_iShowMol = m_iSecondMol;
		ds->Create(m_pCondition);
		ds->Parse();
		ds->BuildName(m_sName);
		m_oaSubDACFs.Add(ds);
		mprintf("\n");

		while (AskYesNo("    Add another aggregation function with different criterion values (y/n)? [no] ",false))
		{
			try { ds = new CDACFSub(); } catch(...) { ds = NULL; }
			if (ds == NULL) NewException((double)sizeof(CDACFSub),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			ds->m_iRefMol = m_iFirstMol;
			ds->m_iShowMol = m_iSecondMol;
			ds->Create(m_pCondition);
			mprintf("\n");
			ds->m_pCondition->Parse_OnlyValues();
			ds->Parse();
			ds->BuildName(m_sName);
			m_oaSubDACFs.Add(ds);
			mprintf("\n");
		}
	} // END IF NOT GRID
	mprintf(WHITE,"\n<<< End of Aggregation Functions <<<\n");
}


void CDACFSub::BuildName(const char *n)
{
	int z, z2, z3;
//	char buf[256], buf2[256];
	CxString buf, buf2;
	CConditionSubGroup *sg;
	CNbSearch *nb;

//	strcpy(buf,n);
	buf.strcpy(n);
	for (z=0;z<m_pCondition->m_oaConditionSubGroups.GetSize();z++)
	{
		sg = (CConditionSubGroup*)m_pCondition->m_oaConditionSubGroups[z];
		for (z2=0;z2<sg->m_oaConditions.GetSize();z2++)
		{
			nb = (CNbSearch*)sg->m_oaConditions[z2];
			if (nb->m_pRDF != NULL)
			{
				if (nb->m_iNbCountMin == -1)
				{
//					strcat(buf,"_r");
					buf.strcat("_r");
					for (z3=0;z3<nb->m_pRDF->m_faMinMaxDist.GetSize()/2;z3++)
					{
						if (z3 < (nb->m_pRDF->m_faMinMaxDist.GetSize()/2)-1)
//							sprintf(buf2,"%.0f-%.0f,",nb->m_pRDF->m_faMinMaxDist[z3*2],nb->m_pRDF->m_faMinMaxDist[z3*2+1]);
							buf2.sprintf("%.0f-%.0f,",nb->m_pRDF->m_faMinMaxDist[z3*2],nb->m_pRDF->m_faMinMaxDist[z3*2+1]);
								else
//									sprintf(buf2,"%.0f-%.0f",nb->m_pRDF->m_faMinMaxDist[z3*2],nb->m_pRDF->m_faMinMaxDist[z3*2+1]);
									buf2.sprintf("%.0f-%.0f",nb->m_pRDF->m_faMinMaxDist[z3*2],nb->m_pRDF->m_faMinMaxDist[z3*2+1]);
					}
//					strcat(buf,buf2);
					buf.strcat(buf2);
				} else
				{
//					sprintf(buf2,"_n%d-%d",nb->m_iNbCountMin+1,nb->m_iNbCountMax+1);
//					strcat(buf,buf2);
					buf2.sprintf("_n%d-%d",nb->m_iNbCountMin+1,nb->m_iNbCountMax+1);
					buf.strcat(buf2);
				}
			}
			if (nb->m_pADF != NULL)
			{
//				strcat(buf,"_a");
				buf.strcat("_a");
				for (z3=0;z3<nb->m_pADF->m_faMinMaxAngle.GetSize()/2;z3++)
				{
					if (z3 < (nb->m_pADF->m_faMinMaxAngle.GetSize()/2)-1)
//						sprintf(buf2,"%.0f-%.0f,",nb->m_pADF->m_faMinMaxAngle[z3*2],nb->m_pADF->m_faMinMaxAngle[z3*2+1]);
						buf2.sprintf("%.0f-%.0f,",nb->m_pADF->m_faMinMaxAngle[z3*2],nb->m_pADF->m_faMinMaxAngle[z3*2+1]);
							else 
//								sprintf(buf2,"%.0f-%.0f",nb->m_pADF->m_faMinMaxAngle[z3*2],nb->m_pADF->m_faMinMaxAngle[z3*2+1]);
								buf2.sprintf("%.0f-%.0f",nb->m_pADF->m_faMinMaxAngle[z3*2],nb->m_pADF->m_faMinMaxAngle[z3*2+1]);
				}
//				strcat(buf,buf2);
				buf.strcat(buf2);
			}
		}
	}

	if (m_bIntermittend)
	{
		if (m_fIntGap != 0)
//			sprintf(buf2,"_int%.3f",m_fIntGap);
			buf2.sprintf("_int%.3f",m_fIntGap);
				else
//					sprintf(buf2,"_int");
					buf2.sprintf("_int");
//		strcat(buf,buf2);
		buf.strcat(buf2);
		if (!m_bIntTravisStyle)
//			strcat(buf,"_mpstyle");
			buf.strcat("_mpstyle");
	}

	if (m_bDistTrace)
	{
//		strcat(buf,"_trace");
		buf.strcat("_trace");
	}

	if (m_bCorrectEq)
	{
//		strcat(buf,"_correq");
		buf.strcat("_correq");
	}

	if (m_bNewMode)
	{
//		strcat(buf,"_new");
		buf.strcat("_new");
	}

	if (m_bBorderMode)
	{
//		strcat(buf,"_border");
		buf.strcat("_border");
	}

	try { m_sName = new char[strlen(buf)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,buf);
}


void CDACFSub::Parse()
{
	mprintf("\n");

	if (g_bDLDisp)
		m_bDistTrace = AskYesNo("    Add values to Lifetime/Displacement plot at cleavage (n) or full trace (y)? [yes] ",true);

	m_bIntermittend = (AskRangeInteger("    Calculate this function as continuous (0) or intermittent (1) type? [0] ",0,1,0) == 1);

	{
		m_bNewMode = false;
		m_bBorderMode = false;
		m_bCorrectEq = false;
	}

	if (m_bIntermittend)
	{
		if (g_bAdvanced2 && (!m_bNewMode))
			m_bIntTravisStyle = AskYesNo("    Use \"bridged style\" (y) or \"true ACF style\" (n) for intermittent function? [yes] ",true);
				else m_bIntTravisStyle = true;

		if (!m_bNewMode)
			m_fIntGap = AskFloat("    Enter maximum allowed gap for intermittent function in ps (0=infinite)? [0] ",0);
				else m_fIntGap = 0;
	}
}


CxVector3 CDACF::CalcCenter(CxVec3Array *v, int i1, int i2)
{
	CxVector3 v1, v2;
	int z, z2, c;
	CSingleMolecule *s1, *s2;

	s1 = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iFirstMol])->m_laSingleMolIndex[i1]];
	s2 = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex[i2]];

	c = 0;
	v1 = CxVector3(0,0,0);
	for (z=0;z<m_pCenterAtoms1->m_baAtomType.GetSize();z++)
	{
		for (z2=0;z2<((CxIntArray*)m_pCenterAtoms1->m_oaAtoms[z])->GetSize();z2++)
		{
			v1 += v->GetAt(((CxIntArray*)s1->m_oaAtomOffset[z])->GetAt(z2)) * m_faWeight1[c];
			c++;
		}
	}
	c = 0;
	v2 = CxVector3(0,0,0);
	for (z=0;z<m_pCenterAtoms2->m_baAtomType.GetSize();z++)
	{
		for (z2=0;z2<((CxIntArray*)m_pCenterAtoms2->m_oaAtoms[z])->GetSize();z2++)
		{
			v2 += v->GetAt(((CxIntArray*)s2->m_oaAtomOffset[z])->GetAt(z2)) * m_faWeight2[c];
			c++;
		}
	}

	v1 = v1 * m_fWeightMol1 + v2 * (1.0f - m_fWeightMol1);

	return v1;
}


void CDACF::FinishDACFSub(CTimeStep *t, CDACFSub *dacfsub)
{
	int z2, z4, z5, z6, z7, z8, ti, j;
//	CSingleMolecule *smfix;
	double tfs, tf, tf2;
	CAggregate *ag;

	ti = 0;
	for (z2=0;z2<((CMolecule*)g_oaMolecules[m_iFirstMol])->m_laSingleMolIndex.GetSize();z2++)
		ti += dacfsub->m_oaAggregates[z2].GetSize();
	mprintf("    Finishing Functions... (%d pairs pending)\n",ti);
	mprintf(WHITE,"    [");
	tfs = ((CMolecule*)g_oaMolecules[m_iFirstMol])->m_laSingleMolIndex.GetSize() / 60.0;
	for (z2=0;z2<((CMolecule*)g_oaMolecules[m_iFirstMol])->m_laSingleMolIndex.GetSize();z2++)
	{
		if (fmod(z2,tfs) < 1.0)
			mprintf(WHITE,"#");

//		smfix = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iFirstMol])->m_laSingleMolIndex[z2]];

		if (dacfsub->m_bIntermittend)
		{
			j = 0;
			for (z4=0;z4<dacfsub->m_oaAggregates[z2].GetSize();z4++)
			{
				ag = (CAggregate*)dacfsub->m_oaAggregates[z2][z4];
				if (ag->m_iEnd == -1)
					ag->m_iEnd = g_iSteps;

				if ((!dacfsub->m_bIntTravisStyle) && (ag->m_iEnd == (long)g_iSteps))
				{
					ag->m_iEnd = g_iSteps;
					ag->m_laLifeIntervals.Add(ag->m_iStart);
					ag->m_laLifeIntervals.Add(ag->m_iEnd);
				}
				if (g_bDDisp || g_bPairMSD)
				{
					if (ag->m_iEnd == (long)g_iSteps)
						ag->m_vEnd = CalcCenter(&t->m_vaCoords_Unfolded,z2,ag->m_iSingleMol2);

					tf = (ag->m_vEnd-ag->m_vFirstStart).GetLength();

					if (g_bDDisp)
						dacfsub->m_pDDisp->AddToBin(tf);

					if (g_bDLDisp)
						dacfsub->m_pDLDisp->AddToBin((ag->m_iEnd - ag->m_iFirstStart+1)*g_fTimestepLength/1000.0,tf);

					if (g_bPairMSD)
						dacfsub->m_pPairMSD->AddToBin_Index(ag->m_iEnd - ag->m_iFirstStart + 1,tf*tf);

					if (g_bVerbose)
						mprintf("\n## %d: (Int) Aggregate %d.%d lived from %d to %d, end at ( %.2f | %.2f | %.2f ), traveled %.2f pm.",g_iSteps,z2,j,ag->m_iFirstStart,ag->m_iEnd,ag->m_vEnd[0],ag->m_vEnd[1],ag->m_vEnd[2],(ag->m_vEnd-ag->m_vFirstStart).GetLength());
				}

				if (g_bDLDF)
					dacfsub->m_pDLDF->AddToBin((ag->m_iEnd - ag->m_iFirstStart + 1)*g_fTimestepLength/1000.0);

				if (g_bDACF)
				{
					if (!dacfsub->m_bIntTravisStyle)
					{
						for (z5=0;z5<ag->m_laLifeIntervals.GetSize()/2;z5++)
						{
							for (z6=ag->m_laLifeIntervals[z5*2];z6<=ag->m_laLifeIntervals[z5*2+1];z6++)
							{
								for (z7=z5;z7<ag->m_laLifeIntervals.GetSize()/2;z7++)
								{
									dacfsub->m_fEqCounter += ag->m_laLifeIntervals[z7*2+1] - ag->m_laLifeIntervals[z7*2] + 1;

									for (z8=ag->m_laLifeIntervals[z7*2];z8<=ag->m_laLifeIntervals[z7*2+1];z8++)
									{
										if (z8 < z6)
											continue;
										dacfsub->m_pDACF->AddToBin_Int(z8-z6);
									}
								}
							}
						}
					} else
					{
						dacfsub->m_fEqCounter += ag->m_iEnd - ag->m_iFirstStart + 1;

						if (dacfsub->m_bNewMode)
						{
							dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2].Add(ag->m_iFirstStart);
							dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2].Add(ag->m_iEnd);
						} else
						{
							if (m_iDACFRes < (int)(ag->m_iEnd-ag->m_iFirstStart+1))
							{
								ti = m_iDACFRes;
								for (z5=0;z5<ti;z5++)
									dacfsub->m_pDACF->AddToBin_Count(z5,(int)(ag->m_iEnd-ag->m_iFirstStart+1)-m_iDACFRes);
							} else ti = (int)(ag->m_iEnd-ag->m_iFirstStart+1);

							for (z5=0;z5<ti;z5++)
								dacfsub->m_pDACF->AddToBin_Count(z5,ti-z5);
						}
					}
				}

				delete ag;
				dacfsub->m_oaAggregates[z2].RemoveAt_NoShrink(z4,1);
				z4--;
				j++;
			} // END FOR Z4
		} else // IF NOT INTERMITTENT
		{
			for (z4=0;z4<dacfsub->m_oaAggregates[z2].GetSize();z4++)
			{
				ag = (CAggregate*)dacfsub->m_oaAggregates[z2][z4];
				if (ag->m_iEnd == -1)
				{
					ag->m_iEnd = g_iSteps;
					if (g_bDDisp || g_bPairMSD)
					{
						ag->m_vEnd = CalcCenter(&t->m_vaCoords_Unfolded,z2,ag->m_iSingleMol2);
						
						tf = (ag->m_vEnd-ag->m_vStart).GetLength();
						tf2 = (ag->m_iEnd-ag->m_iStart+1)*g_fTimestepLength/1000.0;

						if (m_fMaxVel != 0)
						{
							if (tf/tf2 > m_fMaxVel)
							{
								mprintf("\n    Dimer %d - %d lived from Step %d to %d and moved %.2f pm (%.1f pm/ps) - too fast.",z2+1,ag->m_iSingleMol2+1,ag->m_iStart,g_iSteps,tf,tf/tf2*1000.0f);
								if (m_bRemoveMaxVel)
									goto _maxvel2;
							}
						}

						if (g_bDDisp)
							dacfsub->m_pDDisp->AddToBin(tf);

						if (g_bDLDisp)
							dacfsub->m_pDLDisp->AddToBin(tf2,tf);

						if (g_bPairMSD)
							dacfsub->m_pPairMSD->AddToBin_Index(ag->m_iEnd-ag->m_iStart+1,tf*tf);
					}
	_maxvel2:
					if (g_bDLDF)
						dacfsub->m_pDLDF->AddToBin((ag->m_iEnd - ag->m_iStart+1)*g_fTimestepLength/1000.0);

					if (g_bDACF)
					{
						dacfsub->m_fEqCounter += g_iSteps - ag->m_iStart + 1;

						if (dacfsub->m_bNewMode)
						{
							dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2].Add(ag->m_iStart);
							dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2].Add(ag->m_iEnd);
						} else
						{
							if (m_iDACFRes < (int)(g_iSteps-ag->m_iStart+1))
							{
								ti = m_iDACFRes;
								for (z5=0;z5<ti;z5++)
									dacfsub->m_pDACF->AddToBin_Count(z5,(int)(g_iSteps-ag->m_iStart+1)-m_iDACFRes);
							} else ti = (int)(g_iSteps-ag->m_iStart+1);

							for (z5=0;z5<ti;z5++)
								dacfsub->m_pDACF->AddToBin_Count(z5,ti-z5);
						}
					}

					delete ag;
					dacfsub->m_oaAggregates[z2].RemoveAt_NoShrink(z4,1);
					z4--;
					continue;
				}
			} // END FOR Z4
		} // END IF NOT INTERMITTEND
	} // END FOR Z2
	mprintf(WHITE,"]\n\n");
}


void CDACF::UpdateDACFSub(int rm, CTimeStep *t, CDACFSub *dacfsub)
{
	int z3, z4, z5, z6, z7, z8, ti;
	double tf, tf2;
	CAggregate *ag;
	CxVector3 tv;

	for (z4=0;z4<dacfsub->m_oaAggregates[rm].GetSize();z4++)
		((CAggregate*)dacfsub->m_oaAggregates[rm][z4])->m_bStillAlive = false;

	ti = 0;
	for (z3=0;z3<((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize();z3++) // Alle anderen Molekuele durchgehen
	{
		if ((m_iFirstMol == m_iSecondMol) && (z3 == rm)) // Wir wollen nicht das Referenzmolekuel mitzaehlen
			continue;
		if (dacfsub->m_pCondition->Contains(z3))
		{
			ti++;
			for (z4=0;z4<dacfsub->m_oaAggregates[rm].GetSize();z4++)
			{
				ag = (CAggregate*)dacfsub->m_oaAggregates[rm][z4];
				if (ag->m_iSingleMol2 == z3)
				{
					if (dacfsub->m_bIntermittend && (ag->m_iEnd != -1))
					{
						ag->m_iEnd = -1;
						ag->m_iStart = g_iSteps;
						if (g_bDDisp)
							ag->m_vStart = CalcCenter(&t->m_vaCoords_Unfolded,rm,z3);
					}
					ag->m_bStillAlive = true;
					goto _found;
				}
			}

			try { ag = new CAggregate(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAggregate),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			dacfsub->m_oaAggregates[rm].Add(ag);
			ag->m_iSingleMol2 = z3;
			ag->m_bStillAlive = true;
			ag->m_iStart = g_iSteps;
			ag->m_iFirstStart = g_iSteps;
			ag->m_iEnd = -1;
			if (g_bDDisp || g_bPairMSD)
			{
				ag->m_vFirstStart = CalcCenter(&t->m_vaCoords_Unfolded,rm,z3);
				ag->m_vStart = ag->m_vFirstStart;
				if (g_bVerbose)
					mprintf("\n## %d: Aggregate %d.%d spawned at ( %.2f | %.2f | %.2f ).",g_iSteps,rm,dacfsub->m_oaAggregates[rm].GetSize()-1,ag->m_vFirstStart[0],ag->m_vFirstStart[1],ag->m_vFirstStart[2]);
			}
			
_found:;
		}
	} // END FOR Z3

	dacfsub->m_pNDF->AddToBin_Int(ti);

	if (dacfsub->m_bIntermittend)
	{
		for (z4=0;z4<dacfsub->m_oaAggregates[rm].GetSize();z4++)
		{
			ag = (CAggregate*)dacfsub->m_oaAggregates[rm][z4];
			if ((!ag->m_bStillAlive) && (ag->m_iEnd == -1))
			{
				ag->m_iEnd = g_iSteps;

				// Neu 27.08.
				if (g_bDDisp || g_bPairMSD)
					ag->m_vEnd = CalcCenter(&t->m_vaCoords_Unfolded,rm,ag->m_iSingleMol2);

				if (g_bDACF && (!dacfsub->m_bIntTravisStyle))
				{
					ag->m_laLifeIntervals.Add(ag->m_iStart);
					ag->m_laLifeIntervals.Add(g_iSteps);
				}
			}

			if ((dacfsub->m_fIntGap != 0) && (ag->m_iEnd != -1) && (ag->m_iEnd <= g_iSteps-(dacfsub->m_fIntGap/(g_fTimestepLength/1000.0))))
			{
				if (g_bDDisp || g_bPairMSD)
				{
					tf = (ag->m_vEnd-ag->m_vFirstStart).GetLength();

					if (g_bDDisp)
					{
						// Weg 27.08.
						// ag->m_vEnd = CalcCenter(&t->m_vaCoords_Unfolded,rm,ag->m_iSingleMol2);

						dacfsub->m_pDDisp->AddToBin(tf);
						if (g_bDLDisp)
	//						dacfsub->m_pDLDisp->AddToBin((g_iSteps - ag->m_iFirstStart)*g_fTimestepLength/1000.0,(ag->m_vEnd-ag->m_vFirstStart).GetLength());
							dacfsub->m_pDLDisp->AddToBin((ag->m_iEnd - ag->m_iFirstStart + 1)*g_fTimestepLength/1000.0,tf);
					}

					if (g_bPairMSD)
						dacfsub->m_pPairMSD->AddToBin_Index(ag->m_iEnd - ag->m_iFirstStart + 1,tf*tf);
				}

				if (g_bDLDF)
				{
//					dacfsub->m_pDLDF->AddToBin((g_iSteps - ag->m_iFirstStart)*g_fTimestepLength/1000.0);
					dacfsub->m_pDLDF->AddToBin((ag->m_iEnd - ag->m_iFirstStart + 1)*g_fTimestepLength/1000.0);
				}

				if (g_bDACF)
				{
					if (!dacfsub->m_bIntTravisStyle)
					{
						for (z5=0;z5<ag->m_laLifeIntervals.GetSize()/2;z5++)
						{
							for (z6=ag->m_laLifeIntervals[z5*2];z6<=ag->m_laLifeIntervals[z5*2+1];z6++)
							{
								for (z7=z5;z7<ag->m_laLifeIntervals.GetSize()/2;z7++)
								{
									for (z8=ag->m_laLifeIntervals[z7*2];z8<=ag->m_laLifeIntervals[z7*2+1];z8++)
									{
										if (z8 < z6)
											continue;
										dacfsub->m_pDACF->AddToBin(z8-z6);
									}
									dacfsub->m_fEqCounter += ag->m_laLifeIntervals[z7*2+1] - ag->m_laLifeIntervals[z7*2] + 1;
								}
							}
						}
					} else
					{
						dacfsub->m_fEqCounter += ag->m_iEnd-ag->m_iFirstStart+1;

						if (dacfsub->m_bNewMode)
						{
							dacfsub->m_piaIntervals[rm*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2].Add(ag->m_iFirstStart);
							dacfsub->m_piaIntervals[rm*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2].Add(ag->m_iEnd);
						} else
						{
							if (m_iDACFRes < (int)(ag->m_iEnd-ag->m_iFirstStart+1))
							{
								ti = m_iDACFRes;
								for (z5=0;z5<ti;z5++)
									dacfsub->m_pDACF->AddToBin_Count(z5,(int)(ag->m_iEnd-ag->m_iFirstStart+1)-m_iDACFRes);

							} else ti = (int)(ag->m_iEnd-ag->m_iFirstStart+1);

							for (z5=0;z5<ti;z5++)
								dacfsub->m_pDACF->AddToBin_Count(z5,ti-z5);
						}

/*						if (m_iDACFRes < (int)(g_iSteps-ag->m_iFirstStart))
						{
							ti = m_iDACFRes;
							for (z5=0;z5<ti;z5++)
								dacfsub->m_pDACF->AddToBin_Count(z5,(int)(g_iSteps-ag->m_iFirstStart)-m_iDACFRes);
						} else ti = (int)(g_iSteps-ag->m_iFirstStart);
						for (z5=0;z5<ti;z5++)
							dacfsub->m_pDACF->AddToBin_Count(z5,ti-z5);*/
					}
				}

	//			mprintf("Int: RM=%d, OM=%d, t=%d, FirstStart=%d, End=%d, MaxGap ueberschritten.\n",rm+1,ag->m_iSingleMol2+1,g_iSteps,ag->m_iFirstStart,ag->m_iEnd);
				delete ag;
				dacfsub->m_oaAggregates[rm].RemoveAt_NoShrink(z4,1);
				z4--;
				continue;
			} else if (ag->m_iEnd != -1) // END IF finally dead
			{ // broken, but still within gap
				if ((g_bDDisp && dacfsub->m_bDistTrace) || g_bPairMSD)
				{
					tv = CalcCenter(&t->m_vaCoords_Unfolded,rm,ag->m_iSingleMol2);
					ag->m_faTempTrace.Add(g_iSteps - ag->m_iFirstStart + 1);
					ag->m_faTempTrace.Add((tv-ag->m_vFirstStart).GetLength());
				}
			} else // END IF broken, but still within gap
			{ // Currently together
				if ((g_bDDisp && dacfsub->m_bDistTrace) || g_bPairMSD)
				{
					tv = CalcCenter(&t->m_vaCoords_Unfolded,rm,ag->m_iSingleMol2);
					tf = (tv-ag->m_vFirstStart).GetLength();
					if (g_bDDisp && dacfsub->m_bDistTrace)
					{
						dacfsub->m_pDDisp->AddToBin(tf);
						if (g_bDLDisp)
						{
							dacfsub->m_pDLDisp->AddToBin((g_iSteps - ag->m_iFirstStart + 1)*g_fTimestepLength/1000.0,tf);
							for (z5=0;z5<ag->m_faTempTrace.GetSize();z5+=2)
								dacfsub->m_pDLDisp->AddToBin(ag->m_faTempTrace[z5]*g_fTimestepLength/1000.0,ag->m_faTempTrace[z5+1]);
						}
					}
					if (g_bPairMSD)
					{
						dacfsub->m_pPairMSD->AddToBin_Index(g_iSteps - ag->m_iFirstStart + 1,tf*tf);
						for (z5=0;z5<ag->m_faTempTrace.GetSize();z5+=2)
							dacfsub->m_pPairMSD->AddToBin_Index((int)ag->m_faTempTrace[z5],ag->m_faTempTrace[z5+1]*ag->m_faTempTrace[z5+1]);
					}
					ag->m_faTempTrace.RemoveAll_KeepSize();
				}
				if (g_bVerbose)
				{
					tv = CalcCenter(&t->m_vaCoords_Unfolded,rm,ag->m_iSingleMol2);
					tf = (tv-ag->m_vFirstStart).GetLength();
					mprintf("\n## %d: (Int) Aggregate %d.%d lives since %d at ( %.2f | %.2f | %.2f ), traveled %.2f pm.",g_iSteps,rm,z4,ag->m_iFirstStart,tv[0],tv[1],tv[2],tf);
				}
			}
		} // END FOR Z4
	} else // IF NOT INTERMITTENT
	{
		for (z4=0;z4<dacfsub->m_oaAggregates[rm].GetSize();z4++)
		{
			ag = (CAggregate*)dacfsub->m_oaAggregates[rm][z4];
			if ((!ag->m_bStillAlive) && (ag->m_iEnd == -1))
			{
				ag->m_iEnd = g_iSteps;
				if (g_bDDisp || g_bPairMSD)
				{
					ag->m_vEnd = CalcCenter(&t->m_vaCoords_Unfolded,rm,ag->m_iSingleMol2);

					tf = (ag->m_vEnd-ag->m_vStart).GetLength();
					tf2 = (ag->m_iEnd-ag->m_iStart+1)*g_fTimestepLength/1000.0;

					if (g_bVerbose)
						mprintf("\n## %d: Aggregate %d.%d died at ( %.2f | %.2f | %.2f ), traveled %.2f pm.",g_iSteps,rm,z4,ag->m_vEnd[0],ag->m_vEnd[1],ag->m_vEnd[2],tf);
					
					if (m_fMaxVel != 0)
					{
						if (tf/tf2 > m_fMaxVel)
						{
							mprintf("\nDimer %d - %d lived from Step %d to %d and moved %.2f pm (%.1f pm/ps) - too fast.",rm+1,ag->m_iSingleMol2+1,ag->m_iStart,g_iSteps,tf,tf/tf2*1000.0f);
							if (m_bRemoveMaxVel)
								goto _maxvel;
						}
					}
					if (g_bDDisp)
						dacfsub->m_pDDisp->AddToBin(tf);
					if (g_bDLDisp)
						dacfsub->m_pDLDisp->AddToBin(tf2,tf);
					if (g_bPairMSD)
						dacfsub->m_pPairMSD->AddToBin_Index(ag->m_iEnd-ag->m_iStart+1,tf*tf);
				}
		_maxvel:
				if (g_bDLDF)
					dacfsub->m_pDLDF->AddToBin((g_iSteps - ag->m_iStart + 1)*g_fTimestepLength/1000.0);

				if (g_bDACF)
				{
					dacfsub->m_fEqCounter += ag->m_iEnd-ag->m_iStart+1;

					if (dacfsub->m_bNewMode)
					{
				//		mprintf("%d*%d+%d=%d\n",rm,((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize(),ag->m_iSingleMol2,rm*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2);
						dacfsub->m_piaIntervals[rm*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2].Add(ag->m_iStart);
						dacfsub->m_piaIntervals[rm*((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize()+ag->m_iSingleMol2].Add(ag->m_iEnd);
					} else
					{
						if (m_iDACFRes < (int)(ag->m_iEnd-ag->m_iStart+1))
						{
							ti = m_iDACFRes;
							for (z5=0;z5<ti;z5++)
								dacfsub->m_pDACF->AddToBin_Count(z5,(int)(ag->m_iEnd-ag->m_iStart+1)-m_iDACFRes);
						} else ti = (int)(ag->m_iEnd-ag->m_iStart+1);

						for (z5=0;z5<ti;z5++)
							dacfsub->m_pDACF->AddToBin_Count(z5,ti-z5);
					}
				}

	//			mprintf("Con: RM=%d, OM=%d, t=%d, FirstStart=%d, End=%d.\n",rm+1,ag->m_iSingleMol2+1,g_iSteps,ag->m_iFirstStart,ag->m_iEnd);
				delete ag;
				dacfsub->m_oaAggregates[rm].RemoveAt_NoShrink(z4,1);
				z4--;
				continue;
			} else
			{
				if ((g_bDDisp && dacfsub->m_bDistTrace) || g_bPairMSD)
				{
					tv = CalcCenter(&t->m_vaCoords_Unfolded,rm,ag->m_iSingleMol2);
					tf = (tv-ag->m_vStart).GetLength();
					if (g_bDDisp && dacfsub->m_bDistTrace)
						dacfsub->m_pDDisp->AddToBin(tf);
					if (g_bDLDisp && dacfsub->m_bDistTrace)
						dacfsub->m_pDLDisp->AddToBin((g_iSteps - ag->m_iStart+1)*g_fTimestepLength/1000.0,tf);
					if (g_bPairMSD)
						dacfsub->m_pPairMSD->AddToBin_Index(g_iSteps - ag->m_iStart + 1,tf*tf);
				}
				if (g_bVerbose)
				{
					tv = CalcCenter(&t->m_vaCoords_Unfolded,rm,ag->m_iSingleMol2);
					tf = (tv-ag->m_vStart).GetLength();
					mprintf("\n## %d: Aggregate %d.%d lives at ( %.2f | %.2f | %.2f ), traveled %.2f pm.",g_iSteps,rm,z4,ag->m_vEnd[0],ag->m_vEnd[1],ag->m_vEnd[2],tf);
				}
			}
		} // END FOR Z4
	} // END IF NOT INTERMITTEND
}


void CDACF::UpdateNbEx(int rm, CDACFSub *dacfsub)
{
	int z, z2, ti, ti2;
	CNbExchangePair *nx;
//	CNbPair *np;
//	CNbSearch *ns;

//	ns = (CNbSearch*)((CConditionSubGroup*)dacfsub->m_pCondition->m_oaConditionSubGroups[0])->m_oaConditions[0];
	for (z2=0;z2<dacfsub->m_oaNbExPairs[rm].GetSize();z2++)
		((CNbExchangePair*)dacfsub->m_oaNbExPairs[rm][z2])->m_iStillAlive = 0;
	ti = 0;
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iSecondMol])->m_laSingleMolIndex.GetSize();z++) // Alle anderen Molekuele durchgehen
	{
		if ((m_iFirstMol == m_iSecondMol) && (z == rm)) // Wir wollen nicht das Referenzmolekuel mitzaehlen
			continue;
		if (dacfsub->m_pCondition->Contains(z))
		{
			ti++;
			for (z2=0;z2<dacfsub->m_oaNbExPairs[rm].GetSize();z2++)
			{
				nx = (CNbExchangePair*)dacfsub->m_oaNbExPairs[rm][z2];
				if ((int)(nx->m_iShowMol) == z)
				{
					if (nx->m_iEnd == -1) // Zuletzt noch am leben
					{
						nx->m_iStillAlive = 1;
						goto _found;
					}
				}
			}

			try { nx = new CNbExchangePair(); } catch(...) { nx = NULL; }
			if (nx == NULL) NewException((double)sizeof(CNbExchangePair),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			dacfsub->m_oaNbExPairs[rm].Add(nx);
			nx->m_iShowMol = z;
			nx->m_iStart = g_iSteps;
			nx->m_iEnd = -1;
			nx->m_iStillAlive = 1;
_found:;
		}
	} // END FOR Z

	if (!g_bAggregation)
		dacfsub->m_pNDF->AddToBin_Int(ti);

	for (z=0;z<dacfsub->m_oaAggregates[rm].GetSize();z++)
	{
		nx = (CNbExchangePair*)dacfsub->m_oaNbExPairs[rm][z];

		if ((nx->m_iEnd == -1) && (nx->m_iStillAlive == 0))
			nx->m_iEnd = g_iSteps; // Noch am leben? Nicht mehr am leben.

		if (nx->m_iEnd != -1)
		{
			if (nx->m_iEnd < (long)g_iSteps-m_pNbExchange->m_iTimeDepth)
			{
				dacfsub->m_oaAggregates[rm].RemoveAt(z,1);
				z--;
				continue;
			}
		}

//		np = (CNbPair*)ns->m_oaNbPairs[z];

		if (nx->m_iStart < (long)g_iSteps-m_pNbExchange->m_iTimeDepth)
			ti = g_iSteps-m_pNbExchange->m_iTimeDepth;
				else ti = nx->m_iStart;
		if (nx->m_iEnd == -1)
			ti2 = g_iSteps;
				else ti2 = nx->m_iEnd-1; // Weil es bei m_iEnd bereits erstmalig tot ist.
		for (z2=ti;z2<=ti2;z2++)
		{
			// AddToBin Tau=g_iSteps-z2; Pos=np->m_fMinDist; Nb=np->m_iNbPosition
		}
	} // END FOR Z

}


void CDACFSub::Create(CConditionGroup *c)
{
	try { m_pCondition = new CConditionGroup(); } catch(...) { m_pCondition = NULL; }
	if (m_pCondition == NULL) NewException((double)sizeof(CConditionGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pCondition->CopyFrom(c);
	if (g_bAggregation)
	{
		try { m_oaAggregates = new CxObArray[((CMolecule*)g_oaMolecules[m_iRefMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_oaAggregates = NULL; }
		if (m_oaAggregates == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iRefMol])->m_laSingleMolIndex.GetSize()*sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}
}


void CDACF::CalcGridFitParms()
{
	int z, z2;
	CDACFSub *dacfsub;

	try { m_pFitRMin = new double[m_iFitDegreeMax+1]; } catch(...) { m_pFitRMin = NULL; }
	if (m_pFitRMin == NULL) NewException((double)(m_iFitDegreeMax+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_pFitRAvg = new double[m_iFitDegreeMax+1]; } catch(...) { m_pFitRAvg = NULL; }
	if (m_pFitRAvg == NULL) NewException((double)(m_iFitDegreeMax+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_pFitRMax = new double[m_iFitDegreeMax+1]; } catch(...) { m_pFitRMax = NULL; }
	if (m_pFitRMax == NULL) NewException((double)(m_iFitDegreeMax+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<=m_iFitDegreeMax;z++)
	{
		m_pFitRMin[z] = 1;
		m_pFitRAvg[z] = 0;
		m_pFitRMax[z] = 0;
	}

	for (z=0;z<m_oaSubDACFs.GetSize();z++)
	{
		dacfsub = (CDACFSub*)m_oaSubDACFs[z];
		for (z2=m_iFitDegreeMin;z2<=m_iFitDegreeMax;z2++)
		{
			if (m_pFitRMin[z2] > dacfsub->m_pDACF->m_pCorrCoeff[z2])
				m_pFitRMin[z2] = dacfsub->m_pDACF->m_pCorrCoeff[z2];
			if (m_pFitRMax[z2] < dacfsub->m_pDACF->m_pCorrCoeff[z2])
				m_pFitRMax[z2] = dacfsub->m_pDACF->m_pCorrCoeff[z2];
			m_pFitRAvg[z2] += dacfsub->m_pDACF->m_pCorrCoeff[z2];
		}
	}

	for (z=m_iFitDegreeMin;z<=m_iFitDegreeMax;z++)
		m_pFitRAvg[z] /= m_oaSubDACFs.GetSize();
}


void CDACF::CreateGridFit2DF(C2DF *df, int degree, bool intermittend)
{
	int x, y, z;

	switch(m_iGridMode)
	{
		case 3:
			df->m_iRes[0] = m_iGridX;
			df->m_iRes[1] = m_iGridY;
			df->m_fMinVal[0] = m_fGridXMin;
			df->m_fMaxVal[0] = m_fGridXMax;
			df->m_fMinVal[1] = m_fGridYMin;
			df->m_fMaxVal[1] = m_fGridYMax;
			df->Create();
			df->SetLabelX("Distance [pm]");
			df->SetLabelY("Angle [degree]");

			for (x=0;x<m_iGridX;x++)
			{
				for (y=0;y<m_iGridY;y++)
				{
					if (m_bGridCon && m_bGridInt)
					{
						if (intermittend)
							df->m_pBin[y*m_iGridX+x] = ((CDACFSub*)m_oaSubDACFs[(x*m_iGridY+y)*2+1])->m_pDACF->m_pFitIntegral[degree];
								else df->m_pBin[y*m_iGridX+x] = ((CDACFSub*)m_oaSubDACFs[(x*m_iGridY+y)*2])->m_pDACF->m_pFitIntegral[degree];
					} else df->m_pBin[y*m_iGridX+x] = ((CDACFSub*)m_oaSubDACFs[x*m_iGridY+y])->m_pDACF->m_pFitIntegral[degree];
				}
			}

			df->CalcMaxEntry();
			break;

		case 5:
			df->m_iRes[0] = m_iGridY-m_iGridX+1;
			df->m_iRes[1] = m_iGridY-m_iGridX+1;
			df->m_fMinVal[0] = m_iGridX;
			df->m_fMaxVal[0] = m_iGridY;
			df->m_fMinVal[1] = m_iGridX;
			df->m_fMaxVal[1] = m_iGridY;
			df->Create();
			df->SetLabelX("From n-th neighbor");
			df->SetLabelY("To n-th neighbor");

			z = 0;
			for (x=0;x<m_iGridY-m_iGridX+1;x++)
			{
				for (y=x;y<m_iGridY-m_iGridX+1;y++)
				{
					if (m_bGridCon && m_bGridInt)
					{
						if (intermittend)
							df->m_pBin[y*(m_iGridY-m_iGridX+1)+x] = ((CDACFSub*)m_oaSubDACFs[z*2+1])->m_pDACF->m_pFitIntegral[degree];
								else df->m_pBin[y*(m_iGridY-m_iGridX+1)+x] = ((CDACFSub*)m_oaSubDACFs[z*2])->m_pDACF->m_pFitIntegral[degree];
					} else df->m_pBin[y*(m_iGridY-m_iGridX+1)+x] = ((CDACFSub*)m_oaSubDACFs[z])->m_pDACF->m_pFitIntegral[degree];
					z++;
				}
			}

			df->CalcMaxEntry();
			break;
	}
}

/*void CDACF::CreateGridFitDF(CDF *df, int degree, bool intermittend)
{
	int x, y, z;

	switch(m_iGridMode)
	{
	}
}*/


void CDACF::CreateSubDACFStack(char *s)
{
//	CGrace *g;
	int z, z2;
	FILE *a;
//	char buf[256];
	CxString buf;

//	sprintf(buf,"%s.csv",s);
	buf.sprintf("%s.csv",s);
	mprintf("      Writing DACF grid CSV file \"%s\"...\n",(const char*)buf);
	a = OpenFileWrite(buf,true);
	mfprintf(a,"# Tau [ps]");
	for (z=0;z<m_oaSubDACFs.GetSize();z++)
		mfprintf(a,";  %s",((CDACFSub*)m_oaSubDACFs[z])->m_sName);
	mfprintf(a,"\n");
	for (z=0;z<m_iDACFRes;z++)
	{
		mfprintf(a,"%G",z * g_fTimestepLength / 1000.0);
		for (z2=0;z2<m_oaSubDACFs.GetSize();z2++)
			mfprintf(a,";  %G",((CDACFSub*)m_oaSubDACFs[z2])->m_pDACF->m_pBin[z]);
	}
	fclose(a);

/*	sprintf(buf,"%s.agr",s);
	mprintf("      Writing DACF grid AGR file \"%s\"...\n",buf);
	g = new CGrace();
	g->AddDataset();

	g->SetLabelX("Tau [ps]");
	g->SetLabelY("DACF");
	g->SetTitle("DACF Grid");
	g->SetRangeX(0,m_iDACFRes*g_fTimestepLength/1000.0);
	g->SetRangeY(0,1.0);
	g->MakeTicks();
	g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iDACFRes);
	mprintf("      Preparing set  1: [");
	tfs = m_iResolution / 50.0;
	for (z=0;z<m_iResolution;z++)
	{
		if (m_bLeft)
			x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;

		g->AddXYTupel(0,x,m_pBin[z]);

		if (fmod(z,tfs) < 1)
			mprintf(WHITE,"#");
	}
	mprintf("]\n");
	zi = 0;
	for (z0=0;z0<m_iAdditionalSets;z0++)
	{
		if (m_pAdditionalSets[z0] == NULL)
			continue;
		mprintf("      Preparing set %2d: [",zi+2);
		g->AddDataset();
		g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iResolution*2);
		if (m_pAdditionalSetLabels[z0] != NULL)
			g->SetDatasetName(m_pAdditionalSetLabels[z0]);
		for (z=0;z<m_iResolution;z++)
		{
			if (m_bLeft)
				x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
					else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;
			g->AddXYTupel(zi+1,x,m_pAdditionalSets[z0][z]);

			if (fmod(z,tfs) < 1)
				mprintf(WHITE,"#");
		}
		mprintf("]\n");
		zi++;
	}
	g->WriteAgr(buf);*/
}
