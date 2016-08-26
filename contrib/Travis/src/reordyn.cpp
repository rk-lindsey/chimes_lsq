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


#include "reordyn.h"
#include "travis.h"
#include "maintools.h"


CReorDyn::CReorDyn()
{
	m_iShowMol = -1;
	m_bLifetimeSpectrum = false;
	m_bLegendre2 = false;
	m_oaCache.SetName("CReorDyn::m_oaCache");
	m_oaVectors.SetName("CReorDyn::m_oaVectors"); 
	m_oaLTSpectra.SetName("CReorDyn::m_oaLTSpectra");
}


CReorDyn::~CReorDyn()
{
}


void CReorDyn::Parse()
{
	BTIN;
//	char buf[256];
	CxString buf;
	int z2, ti;
	CAtomGroup *ag;
	float tf;
	CFFT tfft;

	try { m_pRDyn = new CDF(); } catch(...) { m_pRDyn = NULL; }
	if (m_pRDyn == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pRDyn->m_bLeft = true;
	m_bCrossCor = false;

	mprintf(YELLOW,"\n>>> Vector Reorientation Dynamics >>>\n\n");

	mprintf("    All atoms will be taken from the OM %s.\n\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);

	m_iVecType = AskRangeInteger("    Should the vector depict position (1), dipole (2), velocity (3) or force (4)? [1] ",1,4,1) - 1;
	mprintf("\n");

	z2 = 0;
	m_iCombinations = 0;
	m_iMolecules = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
	ti = 1;

	if (m_iVecType == 0) // Position
	{
		m_bOrtho = (AskRangeInteger("    Should the vector connect 2 points (0) or stand perpendicular to 3 points (1)? [0] ",0,1,0) != 0);

		do {
			ti = 1;
			if (z2 != 0)
				mprintf("\n    %d. vector\n\n",z2+1);

			if (m_bOrtho)
			{
_ax1:			mprintf("      Please enter the atom(s) at the base point (e.g. C7): ");
				inpprintf("! Please enter the atom(s) at the base point (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _ax1;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
_ax2:			mprintf("      Please enter the 2nd atom(s) of the normal plane (e.g. C7): ");
				inpprintf("! Please enter the 2nd atom(s) of the normal plane (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _ax2;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
_ax3:			mprintf("      Please enter the 3rd atom(s) of the normal plane (e.g. C7): ");
				inpprintf("! Please enter the 3rd atom(s) of the normal plane (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _ax3;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
			} else // IF ORTHO
			{
_ax4:			mprintf("      Please enter the atom(s) at the base point (e.g. C7): ");
				inpprintf("! Please enter the atom(s) at the base point (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _ax4;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
_ax5:			mprintf("      Please enter the atom(s) at the tip point (e.g. C7): ");
				inpprintf("! Please enter the atom(s) at the tip point (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _ax5;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
				m_oaVectors.Add(NULL);
			} // END IF NOT ORTHO
			z2++;
			m_iCombinations += ti;
		} while (AskYesNo("\n    Enter another set of vectors (y/n)? [no] ",false));
	} else if (m_iVecType == 1) // Dipol
	{
		mprintf("    Taking dipole vector of OM %s.\n\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		m_iCombinations = 1;
		g_bDipole = true;
		ParseDipole();
	} else if (m_iVecType == 2) // Geschwindigkeit
	{
_ax6:		mprintf("      Velocity vector of which atoms to use (e.g. C7)? [#2] ");
		inpprintf("! Velocity vector of which atoms to use (e.g. C7)? [#2]\n");
		myget(&buf);

		try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
		if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if (strlen(buf)==0)
		{
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],"#2"))
			{
				eprintf("Weird error.\n");
				inpprintf("! Weird error.\n");
				abort();
			}
		} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
		{
			eprintf("Wrong input.\n");
			inpprintf("! Wrong input.\n");
			goto _ax6;
		}
		m_oaVectors.Add(ag);
		ti *= ag->m_iAtomGes;
	} else if (m_iVecType == 3) // Kraft
	{
_ax7:		mprintf("      Force vector of which atoms to use (e.g. C7)? [#2] ");
		inpprintf("! Force vector of which atoms to use (e.g. C7)? [#2]\n");
		myget(&buf);

		try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
		if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if (strlen(buf)==0)
		{
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],"#2"))
			{
				eprintf("Weird Error.\n");
				inpprintf("! Weird Error.\n");
				abort();
			}
		} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
		{
			eprintf("Wrong input.\n");
			inpprintf("! Wrong input.\n");
			goto _ax7;
		}
		m_oaVectors.Add(ag);
		ti *= ag->m_iAtomGes;
	}

	mprintf("\n    Observing %d vectors per OM.\n\n",m_iCombinations);

//_depth:

	mprintf(WHITE,"    Hint: ");
	mprintf("The resolution of the ACF may never be higher than the number of processed steps.\n");
	mprintf("          Suggested is up to 75 percent of the processed steps.\n\n");
	if (g_iTrajSteps != -1)
		m_iDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the vector ACF (in time steps): [%d] ",int(g_iTrajSteps*0.75),int(g_iTrajSteps*0.75));
			else m_iDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the vector ACF (in time steps): [10000] ",10000);

/*	if (g_bRDynCacheMode)
	{
		tf = g_iTrajSteps*m_iCombinations*((CMolecule*)g_oaMolecules[m_iShowMol])->m_waSingleMolIndex.GetSize()*3.0f*sizeof(float)/1024.0f/1024.0f;
//		if (tf >= 10.0f)
			if (!AskYesNo("    This will occupy around %.0f MB RAM (for each R.Dyn.). Continue (y/n)? [yes] ",true,tf))
				goto _depth;
	} else
	{
		tf = m_iDepth*g_iGesVirtAtomCount*3.0f*sizeof(double)/1024.0f/1024.0f;
//		if (tf >= 10.0f)
			if (!AskYesNo("    This will occupy around %.0f MB RAM (once, not for each R.Dyn.). Continue (y/n)? [yes] ",true,tf))
				goto _depth;
	}*/
	if (g_bACFFFT)
	{
		ti = CalcFFTSize(m_iDepth,false);
		if (m_iDepth != ti)
		{
			mprintf(WHITE,"\n    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n",tfft.NextFastSize(m_iDepth),m_iDepth);
			m_iDepth = tfft.NextFastSize(m_iDepth);
		}
	} else mprintf("\n");

	m_iStride = AskUnsignedInteger("    Take each n-th time step for temporal axis: [1] ",1);

		m_bLegendre2 = false;

	mprintf("\n");

	m_bSpectrum = AskYesNo("    Calculate reorientation spectrum (FFT of vector ACF) (y/n)? [no] ",false);

	if (m_bSpectrum)
	{
		try { m_pACF = new CACF(); } catch(...) { m_pACF = NULL; }
		if (m_pACF == NULL) NewException((double)sizeof(CACF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pACF->m_iSize = m_iDepth;
		m_pACF->m_bSpectrum = true;

		m_pACF->m_bDerivative = AskYesNo("    Derive the vectors before autocorrelation (y/n)? [yes] ",true);

		if (m_pACF->m_bDerivative)
			m_pACF->m_iDerivative = AskRangeInteger("    Please enter degree of vector derivative (1-6): [1] ",1,6,1);
				else m_pACF->m_iDerivative = 0;

		m_pACF->m_bWindowFunction = AskYesNo("    Apply window function (Cos^2) to autocorrelation function (y/n)? [yes] ",true);

		tf = 33356.41 / g_fTimestepLength / 2.0;
		mprintf("\n    A time step length of %.1f fs allows a spectral range up to %.1f cm^-1.\n\n",g_fTimestepLength,tf);
		m_pACF->m_fSpecWaveNumber = AskRangeFloat("    Calculate spectrum up to which wave number (cm^-1)? [%.1f cm^-1] ",0,tf,(tf<5000.0)?tf:5000.0,(tf<5000.0)?tf:5000.0);
		m_pACF->m_iMirror = AskUnsignedInteger("    No mirroring (0), short-time enhancing (1) or long-time enhancing (2)? [1] ",1);
		m_pACF->m_iZeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ",m_iDepth*3,m_iDepth*3);
		m_pACF->m_iZeroPadding0 = m_pACF->m_iZeroPadding;

		ti = CalcFFTSize(m_iDepth+m_pACF->m_iZeroPadding,false);
		if (m_iDepth+m_pACF->m_iZeroPadding != ti)
		{
			mprintf(WHITE,"\n    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n",ti,ti-m_iDepth);
			m_pACF->m_iZeroPadding = ti-m_iDepth;
		}

		m_pACF->m_bACF_DB = AskYesNo("    Convert intensity axis of spectrum to decibel (y/n)? [no] ",false);
		m_pACF->Create();
	}

	BuildName();

	mprintf(YELLOW,"\n<<< End of Vector Reorientation Dynamics <<<\n\n");
	BTOUT;
}


void CReorDyn::ParseSpec()
{
	BTIN;
	int ti;
	float tf;
	CFFT tfft;

	try { m_pRDyn = new CDF(); } catch(...) { m_pRDyn = NULL; }
	if (m_pRDyn == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pRDyn->m_bLeft = true;

	mprintf(YELLOW,"\n>>> Infrared Spectrum >>>\n\n");

	m_iVecType = 1;

	m_iCombinations = 0;
	ti = 1;

	if (m_iShowMol != -1)
	{
		mprintf("    Taking dipole vector of OM %s.\n\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		m_iMolecules = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
	} else m_iMolecules = g_oaSingleMolecules.GetSize();

	m_iCombinations = 1;
	
	g_bDipole = true;
	ParseDipole();

	mprintf(WHITE,"    Hint: ");
	mprintf("The resolution of the ACF may never be higher than the number of processed steps.\n");
	mprintf("          Suggested is 75 percent of the processed steps, but not more than approx. 16384.\n\n");
	if (g_iTrajSteps != -1)
		ti = (int(g_iTrajSteps*0.75)<5120)?int(g_iTrajSteps*0.75):4096;
			else ti = 4096;

	m_iDepth = AskUnsignedInteger("    Enter the resolution (=depth) of the dipole ACF (in time steps): [%d] ",ti,ti);

	if (g_bACFFFT)
	{
		ti = CalcFFTSize(m_iDepth,false);
		if (m_iDepth != ti)
		{
			mprintf(WHITE,"\n    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n",tfft.NextFastSize(m_iDepth),m_iDepth);
			m_iDepth = tfft.NextFastSize(m_iDepth);
		}
	} else mprintf("\n");

	mprintf("\n    This corresponds to a spectral resolution of %.4f cm^-1.\n",33356.41/g_fTimestepLength/2.0/m_iDepth);

	if (g_bAdvanced2)
		m_iStride = AskUnsignedInteger("    Take each n-th time step for temporal axis: [1] ",1);
			else m_iStride = 1;
	
	m_bSpectrum = true;		

	try { m_pACF = new CACF(); } catch(...) { m_pACF = NULL; }
	if (m_pACF == NULL) NewException((double)sizeof(CACF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pACF->m_iSize = m_iDepth;
	m_pACF->m_bSpectrum = true;

	if (g_bAdvanced2)
	{
		mprintf("\n");
		m_pACF->m_bDerivative = AskYesNo("    Derive the dipole vectors before autocorrelation (y/n)? [yes] ",true);

		if (m_pACF->m_bDerivative)
			m_pACF->m_iDerivative = AskRangeInteger("    Please enter degree of vector derivative (1-6): [1] ",1,6,1);
				else m_pACF->m_iDerivative = 0;

		m_pACF->m_bWindowFunction = AskYesNo("    Apply window function (Cos^2) to autocorrelation function (y/n)? [yes] ",true);
	} else
	{
		m_pACF->m_bWindowFunction = true;
		m_pACF->m_bDerivative = true;
		m_pACF->m_iDerivative = 1;
	}

	tf = 33356.41 / g_fTimestepLength / 2.0;
	mprintf("\n    A time step length of %.1f fs allows a spectral range up to %.1f cm^-1.\n\n",g_fTimestepLength,tf);
	m_pACF->m_fSpecWaveNumber = AskRangeFloat("    Calculate spectrum up to which wave number (cm^-1)? [%.1f cm^-1] ",0,tf,(tf<5000.0)?tf:5000.0,(tf<5000.0)?tf:5000.0);

	if (g_bAdvanced2)
	{
		m_pACF->m_iMirror = AskUnsignedInteger("    No mirroring (0), short-time enhancing (1) or long-time enhancing (2)? [1] ",1);
		m_pACF->m_iZeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ",m_iDepth*3,m_iDepth*3);
	} else
	{
		m_pACF->m_iMirror = 1;
		m_pACF->m_iZeroPadding = m_iDepth*3;
	}

	ti = CalcFFTSize(m_iDepth+m_pACF->m_iZeroPadding,false);
	if (m_iDepth+m_pACF->m_iZeroPadding != ti)
	{
		mprintf(WHITE,"    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n",ti,ti-m_iDepth);
		m_pACF->m_iZeroPadding = ti-m_iDepth;
	}

	mprintf("    Zero padding increases the spectral resolution to %.4f cm^-1.\n\n",33356.41/g_fTimestepLength/2.0/(m_iDepth+m_pACF->m_iZeroPadding));

	m_pACF->m_iZeroPadding0 = m_pACF->m_iZeroPadding;

	m_pACF->m_bACF_DB = AskYesNo("    Convert intensity axis of spectrum to decibel (y/n)? [no] ",false);
	m_pACF->Create();

		m_bCrossCor = false;

	BuildName();

	mprintf(YELLOW,"\n<<< End of Infrared Spectrum <<<\n\n");
	BTOUT;
}


void CReorDyn::BuildName()
{
	BTIN;
	int z2;
//	char tmp[32768];
	CxString tmp;
	CAtomGroup *ag;

//	tmp[0] = 0;
	tmp.sprintf("");

//	strcat(tmp,"[");
	tmp.strcat("[");
	if (m_iVecType == 0) // Position
	{
		for (z2=0;z2<m_oaVectors.GetSize()/3;z2++)
		{
			if (m_bOrtho)
			{
				ag = (CAtomGroup*)m_oaVectors[z2*3];
//				strcat(tmp,ag->m_sName);
//				strcat(tmp,"_");
//				strcat(tmp,((CAtomGroup*)m_oaVectors[z2*3+1])->m_sName);
//				strcat(tmp,"_");
//				strcat(tmp,((CAtomGroup*)m_oaVectors[z2*3+2])->m_sName);
				tmp.strcat(ag->m_sName);
				tmp.strcat("_");
				tmp.strcat(((CAtomGroup*)m_oaVectors[z2*3+1])->m_sName);
				tmp.strcat("_");
				tmp.strcat(((CAtomGroup*)m_oaVectors[z2*3+2])->m_sName);
			} else
			{
				ag = (CAtomGroup*)m_oaVectors[z2*3];
//				strcat(tmp,ag->m_sName);
//				strcat(tmp,"_");
//				strcat(tmp,((CAtomGroup*)m_oaVectors[z2*3+1])->m_sName);
				tmp.strcat(ag->m_sName);
				tmp.strcat("_");
				tmp.strcat(((CAtomGroup*)m_oaVectors[z2*3+1])->m_sName);
			}
			if (z2<(m_oaVectors.GetSize()/3)-1)
//				strcat(tmp,"]_[");
				tmp.strcat("]_[");
		}
	} else if (m_iVecType == 1) // Dipol
	{
//		strcat(tmp,"dip_");
		tmp.strcat("dip_");
		if (m_iShowMol == -1)
//			strcat(tmp,"global");
			tmp.strcat("global");
		else
//			strcat(tmp,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			tmp.strcat(((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	} else if (m_iVecType == 2) // Geschwindigkeit
	{
//		strcat(tmp,"vel_");
//		strcat(tmp,((CAtomGroup*)m_oaVectors[0])->m_sName);
		tmp.strcat("vel_");
		tmp.strcat(((CAtomGroup*)m_oaVectors[0])->m_sName);
	} else if (m_iVecType == 3) // Kraft
	{
//		strcat(tmp,"frc_");
//		strcat(tmp,((CAtomGroup*)m_oaVectors[0])->m_sName);
		tmp.strcat("frc_");
		tmp.strcat(((CAtomGroup*)m_oaVectors[0])->m_sName);
	}
//	strcat(tmp,"]");
	tmp.strcat("]");

	if (m_bLegendre2)
//		strcat(tmp,"_P2");
		tmp.strcat("_P2");

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	if (m_iVecType != 1)
	{
		strcpy(m_sShortName,tmp);
//		sprintf(tmp,"%s_",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
//		strcat(tmp,m_sShortName);
		tmp.sprintf("%s_",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.strcat(m_sShortName);
	}

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CReorDyn::BuildAtomList(CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z, z1t, z1a, z2t, z2a, z3t, z3a;
	CAtomGroup *g1, *g2, *g3;
	CxIntArray *a1, *a2, *a3;

	vec->RemoveAll_KeepSize();
	for (z=0;z<m_oaVectors.GetSize()/3;z++)
	{
		if (m_bOrtho)
		{
			g1 = (CAtomGroup*)m_oaVectors[z*3];
			g2 = (CAtomGroup*)m_oaVectors[z*3+1];
			g3 = (CAtomGroup*)m_oaVectors[z*3+2];
			for (z1t=0;z1t<g1->m_baAtomType.GetSize();z1t++)
			{
				a1 = (CxIntArray*)g1->m_oaAtoms[z1t];
				for (z1a=0;z1a<a1->GetSize();z1a++)
				{
					for (z2t=0;z2t<g2->m_baAtomType.GetSize();z2t++)
					{
						a2 = (CxIntArray*)g2->m_oaAtoms[z2t];
						for (z2a=0;z2a<a2->GetSize();z2a++)
						{
							for (z3t=0;z3t<g3->m_baAtomType.GetSize();z3t++)
							{
								a3 = (CxIntArray*)g3->m_oaAtoms[z3t];
								for (z3a=0;z3a<a3->GetSize();z3a++)
								{
									vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
									vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
									vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
								}
							}
						}
					}
				}
			}
		} else
		{
			g1 = (CAtomGroup*)m_oaVectors[z*3];
			g2 = (CAtomGroup*)m_oaVectors[z*3+1];
			for (z1t=0;z1t<g1->m_baAtomType.GetSize();z1t++)
			{
				a1 = (CxIntArray*)g1->m_oaAtoms[z1t];
				for (z1a=0;z1a<a1->GetSize();z1a++)
				{
					for (z2t=0;z2t<g2->m_baAtomType.GetSize();z2t++)
					{
						a2 = (CxIntArray*)g2->m_oaAtoms[z2t];
						for (z2a=0;z2a<a2->GetSize();z2a++)
						{
							vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
							vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
							vec->Add(0);
						}
					}
				}
			}
		}
	}
	BXOUT;
}


void CReorDyn::Finish(const char *multibuf)
{
	int z2, z3, z4, z5, z6, ti;
	float tfs;
//	char buf[256];
	CxString buf;
	CxFloatArray *ptfa, *ptfab, *ptfa2, *ptfa2b, *ptfa3;
	CAutoCorrelation *ac;
	CCrossCorrelation *ccr;

	ti = 0;

	if (g_bRDynCacheMode)
	{
		if (m_bSpectrum)
		{
			switch(m_pACF->m_iDerivative)
			{
				case 0:
					mprintf("    Not deriving vectors.\n");
					ti = 0;
					break;

				case 1:
					ti = 2;
					mprintf("    Deriving vectors (1st derivative)...\n");
					mprintf(WHITE,"      [");
					tfs = (m_iMolecules*m_iCombinations)/60.0;
					for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
					{
						if (fmod(z2,tfs) < 1.0)
							mprintf(WHITE,"#");
						ptfa = (CxFloatArray*)m_oaCache[z2];
						for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
						{
							(*ptfa)[z3*3]   = 0.5 * ((*ptfa)[(z3+2)*3]   - (*ptfa)[z3*3]);
							(*ptfa)[z3*3+1] = 0.5 * ((*ptfa)[(z3+2)*3+1] - (*ptfa)[z3*3+1]);
							(*ptfa)[z3*3+2] = 0.5 * ((*ptfa)[(z3+2)*3+2] - (*ptfa)[z3*3+2]);
						}
					}
					mprintf(WHITE,"]\n");
					break;

				case 2:
					ti = 2;
					mprintf("    Deriving vectors (2nd derivative)...\n");
					mprintf(WHITE,"      [");
					tfs = (m_iMolecules*m_iCombinations)/60.0;
					for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
					{
						if (fmod(z2,tfs) < 1.0)
							mprintf(WHITE,"#");
						ptfa = (CxFloatArray*)m_oaCache[z2];
						for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
						{
							(*ptfa)[z3*3]   = (*ptfa)[(z3+2)*3]   + (*ptfa)[z3*3]   - 2*(*ptfa)[(z3+1)*3];
							(*ptfa)[z3*3+1] = (*ptfa)[(z3+2)*3+1] + (*ptfa)[z3*3+1] - 2*(*ptfa)[(z3+1)*3+1];
							(*ptfa)[z3*3+2] = (*ptfa)[(z3+2)*3+2] + (*ptfa)[z3*3+2] - 2*(*ptfa)[(z3+1)*3+2];
						}
					}
					mprintf(WHITE,"]\n");
					break;

				case 3:
					ti = 4;
					mprintf("    Deriving vectors (3rd derivative)...\n");
					mprintf(WHITE,"      [");
					tfs = (m_iMolecules*m_iCombinations)/60.0;
					for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
					{
						if (fmod(z2,tfs) < 1.0)
							mprintf(WHITE,"#");
						ptfa = (CxFloatArray*)m_oaCache[z2];
						for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
						{
							(*ptfa)[z3*3]   = 0.5 * ((*ptfa)[(z3+4)*3]   - 2*(*ptfa)[(z3+3)*3]   + 2*(*ptfa)[(z3+1)*3]   - (*ptfa)[z3*3]);
							(*ptfa)[z3*3+1] = 0.5 * ((*ptfa)[(z3+4)*3+1] - 2*(*ptfa)[(z3+3)*3+1] + 2*(*ptfa)[(z3+1)*3+1] - (*ptfa)[z3*3+1]);
							(*ptfa)[z3*3+2] = 0.5 * ((*ptfa)[(z3+4)*3+2] - 2*(*ptfa)[(z3+3)*3+2] + 2*(*ptfa)[(z3+1)*3+2] - (*ptfa)[z3*3+2]);
						}
					}
					mprintf(WHITE,"]\n");
					break;

				case 4:
					ti = 4;
					mprintf("    Deriving vectors (4th derivative)...\n");
					mprintf(WHITE,"      [");
					tfs = (m_iMolecules*m_iCombinations)/60.0;
					for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
					{
						if (fmod(z2,tfs) < 1.0)
							mprintf(WHITE,"#");
						ptfa = (CxFloatArray*)m_oaCache[z2];
						for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
						{
							(*ptfa)[z3*3]   = (*ptfa)[(z3+4)*3]   - 4*(*ptfa)[(z3+3)*3]   + 6*(*ptfa)[(z3+2)*3]   - 4*(*ptfa)[(z3+1)*3]   + (*ptfa)[z3*3];
							(*ptfa)[z3*3+1] = (*ptfa)[(z3+4)*3+1] - 4*(*ptfa)[(z3+3)*3+1] + 6*(*ptfa)[(z3+2)*3+1] - 4*(*ptfa)[(z3+1)*3+1] + (*ptfa)[z3*3+1];
							(*ptfa)[z3*3+2] = (*ptfa)[(z3+4)*3+2] - 4*(*ptfa)[(z3+3)*3+2] + 6*(*ptfa)[(z3+2)*3+2] - 4*(*ptfa)[(z3+1)*3+2] + (*ptfa)[z3*3+2];
						}
					}
					mprintf(WHITE,"]\n");
					break;

				case 5:
					ti = 6;
					mprintf("    Deriving vectors (5th derivative)...\n");
					mprintf(WHITE,"      [");
					tfs = (m_iMolecules*m_iCombinations)/60.0;
					for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
					{
						if (fmod(z2,tfs) < 1.0)
							mprintf(WHITE,"#");
						ptfa = (CxFloatArray*)m_oaCache[z2];
						for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
						{
							(*ptfa)[z3*3]   = 0.5 * ((*ptfa)[(z3+6)*3]   - 4*(*ptfa)[(z3+5)*3]   + 5*(*ptfa)[(z3+4)*3]   - 5*(*ptfa)[(z3+2)*3]   + 4*(*ptfa)[(z3+1)*3]   - (*ptfa)[z3*3]);
							(*ptfa)[z3*3+1] = 0.5 * ((*ptfa)[(z3+6)*3+1] - 4*(*ptfa)[(z3+5)*3+1] + 5*(*ptfa)[(z3+4)*3+1] - 5*(*ptfa)[(z3+2)*3+1] + 4*(*ptfa)[(z3+1)*3+1] - (*ptfa)[z3*3+1]);
							(*ptfa)[z3*3+2] = 0.5 * ((*ptfa)[(z3+6)*3+2] - 4*(*ptfa)[(z3+5)*3+2] + 5*(*ptfa)[(z3+4)*3+2] - 5*(*ptfa)[(z3+2)*3+2] + 4*(*ptfa)[(z3+1)*3+2] - (*ptfa)[z3*3+2]);
						}
					}
					mprintf(WHITE,"]\n");
					break;

				case 6:
					ti = 6;
					mprintf("    Deriving vectors (6th derivative)...\n");
					mprintf(WHITE,"      [");
					tfs = (m_iMolecules*m_iCombinations)/60.0;
					for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
					{
						if (fmod(z2,tfs) < 1.0)
							mprintf(WHITE,"#");
						ptfa = (CxFloatArray*)m_oaCache[z2];
						for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
						{
							(*ptfa)[z3*3]   = (*ptfa)[(z3+6)*3]   - 6*(*ptfa)[(z3+5)*3]   + 15*(*ptfa)[(z3+4)*3]   - 20*(*ptfa)[(z3+3)*3]   + 15*(*ptfa)[(z3+2)*3]   - 6*(*ptfa)[(z3+1)*3]   + (*ptfa)[z3*3];
							(*ptfa)[z3*3+1] = (*ptfa)[(z3+6)*3+1] - 6*(*ptfa)[(z3+5)*3+1] + 15*(*ptfa)[(z3+4)*3+1] - 20*(*ptfa)[(z3+3)*3+1] + 15*(*ptfa)[(z3+2)*3+1] - 6*(*ptfa)[(z3+1)*3+1] + (*ptfa)[z3*3+1];
							(*ptfa)[z3*3+2] = (*ptfa)[(z3+6)*3+2] - 6*(*ptfa)[(z3+5)*3+2] + 15*(*ptfa)[(z3+4)*3+2] - 20*(*ptfa)[(z3+3)*3+2] + 15*(*ptfa)[(z3+2)*3+2] - 6*(*ptfa)[(z3+1)*3+2] + (*ptfa)[z3*3+2];
						}
					}
					mprintf(WHITE,"]\n");
					break;

				default:
					ti = 0;
					eprintf("    Error in CReorDyn::Finish().\n");
			}
		} else // if not spectrum
		{
			ti = 0;
		}

		if (((long)g_iSteps)/g_iStride-ti <= m_iDepth)
		{
			eprintf("Error: Not enough time steps analyzed (less than requested correlation depth %d).\n",m_iDepth);
			goto _irdone;
		}

		mprintf("    Autocorrelating cached vectors...\n");
		mprintf(WHITE,"      [");
		tfs = (m_iMolecules*m_iCombinations)/60.0;

		try { ptfa2 = new CxFloatArray("CReorDyn::Finish():ptfa2"); } catch(...) { ptfa2 = NULL; }
		if (ptfa2 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ptfa2->SetSize(g_iSteps/g_iStride-ti);

		try { ptfa3 = new CxFloatArray("CReorDyn::Finish():ptfa3"); } catch(...) { ptfa3 = NULL; }
		if (ptfa3 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ptfa3->SetSize(g_iSteps/g_iStride-ti);

		try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
		if (ac == NULL) NewException((double)sizeof(CAutoCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ac->Init(g_iSteps/g_iStride-ti,m_iDepth,g_bACFFFT);
		for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
		{
			if (fmod(z2,tfs) < 1.0)
				mprintf(WHITE,"#");
			ptfa = (CxFloatArray*)m_oaCache[z2];

			if (m_bLegendre2) // 2nd Legendre Polynomial of dot product
			{
				/* 3/2 X^2 */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3] * (*ptfa)[z3*3];
				ac->AutoCorrelate(ptfa2,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,1.5*(*ptfa3)[z3*m_iStride]);

				/* 3/2 Y^2 */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3+1] * (*ptfa)[z3*3+1];
				ac->AutoCorrelate(ptfa2,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,1.5*(*ptfa3)[z3*m_iStride]);

				/* 3/2 Z^2 */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3+2] * (*ptfa)[z3*3+2];
				ac->AutoCorrelate(ptfa2,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,1.5*(*ptfa3)[z3*m_iStride]);

				/* 3 X*Y */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3] * (*ptfa)[z3*3+1];
				ac->AutoCorrelate(ptfa2,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,3*(*ptfa3)[z3*m_iStride]);

				/* 3 X*Z */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3] * (*ptfa)[z3*3+2];
				ac->AutoCorrelate(ptfa2,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,3*(*ptfa3)[z3*m_iStride]);

				/* 3 Y*Z */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3+1] * (*ptfa)[z3*3+2];
				ac->AutoCorrelate(ptfa2,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,3*(*ptfa3)[z3*m_iStride]);

				/* -1/2 */
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,-0.5);

				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->m_fBinEntries += (double)(g_iSteps-z3*m_iStride) - 3.0;
			} else // Classical vector autocorrelation
			{
				/* X */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3];
	//					(*ptfa2)[z3] = sqrt((*ptfa)[z3*3]*(*ptfa)[z3*3] + (*ptfa)[z3*3+1]*(*ptfa)[z3*3+1] + (*ptfa)[z3*3+2]*(*ptfa)[z3*3+2]);

				ac->AutoCorrelate(ptfa2,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);

				/* Y */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3+1];
				ac->AutoCorrelate(ptfa2,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);

				/* Z */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
					(*ptfa2)[z3] = (*ptfa)[z3*3+2];
				ac->AutoCorrelate(ptfa2,ptfa3);

				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
				{
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);
					m_pRDyn->m_fBinEntries += (double)(g_iSteps-z3*m_iStride) - 3.0;
				}
			}
		}
		delete ac;
		delete ptfa2;
		delete ptfa3;

		mprintf(WHITE,"]\n");
	} else
	{
		for (z2=0;z2<m_pRDyn->m_iResolution;z2++)
			m_pRDyn->m_pBin[z2] /= m_pCount[z2];
	}

	if (m_iShowMol == -1) // Keine Normalisierung fuer globales Spektrum
		Finish_Part2("",multibuf,1); 
			else Finish_Part2("",multibuf,m_iMolecules*m_iCombinations);

	if (m_bSpectrum && m_bCrossCor)
	{

/*************************************************************************/

		mprintf(WHITE,"\n  * Computing Cross-Correlation \"cross\"\n");

		m_pACF->m_iSize = m_iDepth;
		m_pACF->m_iZeroPadding = m_pACF->m_iZeroPadding0;
		m_pACF->Create();

		mprintf("    Cross-Correlating cached vectors...\n");
		mprintf(WHITE,"      [");
		tfs = (m_iMolecules*m_iCombinations)*(m_iMolecules*m_iCombinations)/60.0;
		m_pRDyn->ZeroBin();

		try { ptfa2 = new CxFloatArray("CReorDyn::Finish():ptfa2"); } catch(...) { ptfa2 = NULL; }
		if (ptfa2 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ptfa2->SetSize(g_iSteps/g_iStride-ti);

		try { ptfa2b = new CxFloatArray("CReorDyn::Finish():ptfa2b"); } catch(...) { ptfa2b = NULL; }
		if (ptfa2b == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ptfa2b->SetSize(g_iSteps/g_iStride-ti);

		try { ptfa3 = new CxFloatArray("CReorDyn::Finish():ptfa3"); } catch(...) { ptfa3 = NULL; }
		if (ptfa3 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ptfa3->SetSize(g_iSteps/g_iStride-ti);

		try { ccr = new CCrossCorrelation(); } catch(...) { ccr = NULL; }
		if (ccr == NULL) NewException((double)sizeof(CCrossCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ccr->Init(g_iSteps/g_iStride-ti,m_iDepth,g_bACFFFT);
		for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
		{
			ptfa = (CxFloatArray*)m_oaCache[z2];

			for (z4=0;z4<m_iMolecules*m_iCombinations;z4++)
			{
				if (fmod(z4+z2*m_iMolecules*m_iCombinations,tfs) < 1.0)
					mprintf(WHITE,"#");

				if (z2 == z4)
					continue;

				ptfab = (CxFloatArray*)m_oaCache[z4];

				/* X */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
				{
					(*ptfa2)[z3]  = (*ptfa)[z3*3];
					(*ptfa2b)[z3] = (*ptfab)[z3*3];
				}
				ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);

				/* Y */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
				{
					(*ptfa2)[z3]  = (*ptfa)[z3*3+1];
					(*ptfa2b)[z3] = (*ptfab)[z3*3+1];
				}
				ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);

				/* Z */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
				{
					(*ptfa2)[z3]  = (*ptfa)[z3*3+2];
					(*ptfa2b)[z3] = (*ptfab)[z3*3+2];
				}
				ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);

				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
				{
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);
					m_pRDyn->m_fBinEntries += (double)(g_iSteps-z3*m_iStride) - 3.0;
				}
			}
		}
		delete ccr;
		delete ptfa2;
		delete ptfa2b;
		delete ptfa3;

		mprintf(WHITE,"]\n");

		Finish_Part2("cross_",multibuf,1);

/*************************************************************************/

		mprintf(WHITE,"\n  * Computing Cross-Correlation \"both\"\n");

		m_pACF->m_iSize = m_iDepth;
		m_pACF->m_iZeroPadding = m_pACF->m_iZeroPadding0;
		m_pACF->Create();

		mprintf("    Cross-Correlating cached vectors...\n");
		mprintf(WHITE,"      [");
		tfs = (m_iMolecules*m_iCombinations)*(m_iMolecules*m_iCombinations)/60.0;
		m_pRDyn->ZeroBin();

		try { ptfa2 = new CxFloatArray("CReorDyn::Finish():ptfa2"); } catch(...) { ptfa2 = NULL; }
		if (ptfa2 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ptfa2->SetSize(g_iSteps/g_iStride-ti);

		try { ptfa2b = new CxFloatArray("CReorDyn::Finish():ptfa2b"); } catch(...) { ptfa2b = NULL; }
		if (ptfa2b == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ptfa2b->SetSize(g_iSteps/g_iStride-ti);

		try { ptfa3 = new CxFloatArray("CReorDyn::Finish():ptfa3"); } catch(...) { ptfa3 = NULL; }
		if (ptfa3 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ptfa3->SetSize(g_iSteps/g_iStride-ti);

		try { ccr = new CCrossCorrelation(); } catch(...) { ccr = NULL; }
		if (ccr == NULL) NewException((double)sizeof(CCrossCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ccr->Init(g_iSteps/g_iStride-ti,m_iDepth,g_bACFFFT);
		for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
		{
			ptfa = (CxFloatArray*)m_oaCache[z2];

			for (z4=0;z4<m_iMolecules*m_iCombinations;z4++)
			{
				if (fmod(z4+z2*m_iMolecules*m_iCombinations,tfs) < 1.0)
					mprintf(WHITE,"#");

				ptfab = (CxFloatArray*)m_oaCache[z4];

				/* X */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
				{
					(*ptfa2)[z3]  = (*ptfa)[z3*3];
					(*ptfa2b)[z3] = (*ptfab)[z3*3];
				}
				ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);

				/* Y */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
				{
					(*ptfa2)[z3]  = (*ptfa)[z3*3+1];
					(*ptfa2b)[z3] = (*ptfab)[z3*3+1];
				}
				ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);
				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);

				/* Z */
				for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
				{
					(*ptfa2)[z3]  = (*ptfa)[z3*3+2];
					(*ptfa2b)[z3] = (*ptfab)[z3*3+2];
				}
				ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);

				for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
				{
					m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);
					m_pRDyn->m_fBinEntries += (double)(g_iSteps-z3*m_iStride) - 3.0;
				}
			}
		}
		delete ccr;
		delete ptfa2;
		delete ptfa2b;
		delete ptfa3;

		mprintf(WHITE,"]\n");

		Finish_Part2("both_",multibuf,1);

/*************************************************************************/

		if (m_iShowMol == -1)
		{
			for (z5=0;z5<g_oaMolecules.GetSize();z5++)
			{
				for (z6=z5;z6<g_oaMolecules.GetSize();z6++)
				{
					mprintf(WHITE,"\n  * Computing Cross-Correlation \"%s - %s\"\n",((CMolecule*)g_oaMolecules[z5])->m_sName,((CMolecule*)g_oaMolecules[z6])->m_sName);
					m_pACF->m_iSize = m_iDepth;
					m_pACF->m_iZeroPadding = m_pACF->m_iZeroPadding0;
					m_pACF->Create();

					mprintf("    Cross-Correlating cached vectors...\n");
					mprintf(WHITE,"      [");
					tfs = (m_iMolecules*m_iCombinations)*(m_iMolecules*m_iCombinations)/60.0;
					m_pRDyn->ZeroBin();

					try { ptfa2 = new CxFloatArray("CReorDyn::Finish():ptfa2"); } catch(...) { ptfa2 = NULL; }
					if (ptfa2 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					ptfa2->SetSize(g_iSteps/g_iStride-ti);

					try { ptfa2b = new CxFloatArray("CReorDyn::Finish():ptfa2b"); } catch(...) { ptfa2b = NULL; }
					if (ptfa2b == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					ptfa2b->SetSize(g_iSteps/g_iStride-ti);

					try { ptfa3 = new CxFloatArray("CReorDyn::Finish():ptfa3"); } catch(...) { ptfa3 = NULL; }
					if (ptfa3 == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					ptfa3->SetSize(g_iSteps/g_iStride-ti);

					try { ccr = new CCrossCorrelation(); } catch(...) { ccr = NULL; }
					if (ccr == NULL) NewException((double)sizeof(CCrossCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					ccr->Init(g_iSteps/g_iStride-ti,m_iDepth,g_bACFFFT);

					for (z2=0;z2<m_iMolecules*m_iCombinations;z2++)
					{
						ptfa = (CxFloatArray*)m_oaCache[z2];

						for (z4=0;z4<m_iMolecules*m_iCombinations;z4++)
						{
							if (fmod(z4+z2*m_iMolecules*m_iCombinations,tfs) < 1.0)
								mprintf(WHITE,"#");

							if (((((CSingleMolecule*)g_oaSingleMolecules[z2])->m_iMolType != z5) || (((CSingleMolecule*)g_oaSingleMolecules[z4])->m_iMolType != z6)) &&
							    ((((CSingleMolecule*)g_oaSingleMolecules[z2])->m_iMolType != z6) || (((CSingleMolecule*)g_oaSingleMolecules[z4])->m_iMolType != z5)))
								continue;

							ptfab = (CxFloatArray*)m_oaCache[z4];

							/* X */
							for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
							{
								(*ptfa2)[z3]  = (*ptfa)[z3*3];
								(*ptfa2b)[z3] = (*ptfab)[z3*3];
							}
							ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);
							for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
								m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);

							/* Y */
							for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
							{
								(*ptfa2)[z3]  = (*ptfa)[z3*3+1];
								(*ptfa2b)[z3] = (*ptfab)[z3*3+1];
							}
							ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);
							for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
								m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);

							/* Z */
							for (z3=0;z3<(int)g_iSteps/g_iStride-ti;z3++)
							{
								(*ptfa2)[z3]  = (*ptfa)[z3*3+2];
								(*ptfa2b)[z3] = (*ptfab)[z3*3+2];
							}
							ccr->CrossCorrelate(ptfa2,ptfa2b,ptfa3);

							for (z3=0;z3<(int)m_iDepth/m_iStride;z3++)
							{
								m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3*m_iStride]);
								m_pRDyn->m_fBinEntries += (double)(g_iSteps-z3*m_iStride) - 3.0;
							}
						}
					}
					delete ccr;
					delete ptfa2;
					delete ptfa2b;
					delete ptfa3;

					mprintf(WHITE,"]\n");

//					sprintf(buf,"cc_%s_%s_",((CMolecule*)g_oaMolecules[z5])->m_sName,((CMolecule*)g_oaMolecules[z6])->m_sName);
					buf.sprintf("cc_%s_%s_",((CMolecule*)g_oaMolecules[z5])->m_sName,((CMolecule*)g_oaMolecules[z6])->m_sName);
					Finish_Part2(buf,multibuf,1);
				}
			}
		}

/*************************************************************************/

	} // END IF SPECTRUM
_irdone:;
}


void CReorDyn::Finish_Part2(const char *s, const char *multibuf, int nmol)
{
//	char buf[256];
	CxString buf;
	int z2;

	mprintf("    %.0f bin entries.\n",m_pRDyn->m_fBinEntries);

	if (m_bSpectrum)
		m_pRDyn->MultiplyBin(1.0/g_iSteps/nmol);
			else m_pRDyn->MultiplyBin(1.0/m_pRDyn->m_pBin[0]);

//	sprintf(buf,"rdyn_%s%s%s.csv",s,m_sName,multibuf);
	buf.sprintf("rdyn_%s%s%s.csv",s,m_sName,multibuf);
	mprintf("    Saving result as %s ...\n",(const char*)buf);
	m_pRDyn->Write("",buf,"",false);


	if (m_bSpectrum)
	{
		mprintf("    Creating reorientation spectrum:\n");

		for (z2=0;z2<m_iDepth;z2++)
 			m_pACF->m_pData[z2] = m_pRDyn->m_pBin[z2];

		if (m_pACF->m_iDerivative != 0)
//			sprintf(buf,"acf_%s%s%s.d%d.csv",s,m_sName,multibuf,m_pACF->m_iDerivative);
			buf.sprintf("acf_%s%s%s.d%d.csv",s,m_sName,multibuf,m_pACF->m_iDerivative);
		else
//			sprintf(buf,"acf_%s%s%s.csv",s,m_sName,multibuf);
			buf.sprintf("acf_%s%s%s.csv",s,m_sName,multibuf);

		mprintf("      Saving ACF as %s ...\n",(const char*)buf);
		m_pACF->WriteACF("",buf,"");

		if (m_pACF->m_iMirror != 0)
		{
			mprintf("      Mirroring ACF...\n");
			m_pACF->Mirror(m_pACF->m_iMirror);

			if (m_pACF->m_bDerivative)
//				sprintf(buf,"acf_%s%s%s.d%d.m.csv",s,m_sName,multibuf,m_pACF->m_iDerivative);
				buf.sprintf("acf_%s%s%s.d%d.m.csv",s,m_sName,multibuf,m_pACF->m_iDerivative);
			else
//				sprintf(buf,"acf_%s%s%s.m.csv",s,m_sName,multibuf);
				buf.sprintf("acf_%s%s%s.m.csv",s,m_sName,multibuf);

			mprintf("      Saving mirrored ACF as %s ...\n",(const char*)buf);
			m_pACF->WriteACF("",buf,"");
		}

		if (m_pACF->m_bWindowFunction)
		{
			mprintf("      Applying window function to ACF...\n");
			m_pACF->Window();

			if (m_pACF->m_bDerivative)
//				sprintf(buf,"acf_%s%s%s.d%d.w.csv",s,m_sName,multibuf,m_pACF->m_iDerivative);
				buf.sprintf("acf_%s%s%s.d%d.w.csv",s,m_sName,multibuf,m_pACF->m_iDerivative);
			else
//				sprintf(buf,"acf_%s%s%s.w.csv",s,m_sName,multibuf);
				buf.sprintf("acf_%s%s%s.w.csv",s,m_sName,multibuf);

			mprintf("      Saving windowed ACF as %s ...\n",(const char*)buf);
			m_pACF->WriteACF("",buf,"");
		}

		mprintf("      Performing fourier transformation...\n");

		try { g_pFFT = new CFFT(); } catch(...) { g_pFFT = NULL; }
		if (g_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pFFT->PrepareFFT_C2C(m_pACF->m_iSize+m_pACF->m_iZeroPadding);
		m_pACF->Transform(g_pFFT);
		delete g_pFFT;
		m_pACF->m_pSpectrum->SetMaxRWL(1E7f/299.792f/g_fTimestepLength/g_iStride);

		if (m_pACF->m_bACF_DB)
		{
			mprintf("      Normalising spectrum to decibel...\n");
			m_pACF->m_pSpectrum->MakeDB();
		}

		if (m_pACF->m_bDerivative)
//			sprintf(buf,"spectrum_%s%s%s.d%d.csv",s,m_sName,multibuf,m_pACF->m_iDerivative);
			buf.sprintf("spectrum_%s%s%s.d%d.csv",s,m_sName,multibuf,m_pACF->m_iDerivative);
		else
//			sprintf(buf,"spectrum_%s%s%s.csv",s,m_sName,multibuf);
			buf.sprintf("spectrum_%s%s%s.csv",s,m_sName,multibuf);

		mprintf("      Saving spectrum as %s ...\n",(const char*)buf);
		m_pACF->m_pSpectrum->Write("",buf,"");
	}
}
