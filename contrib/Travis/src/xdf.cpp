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

#include "xdf.h"
#include "globalvar.h"
#include "maintools.h"


CXDF::CXDF()
{
	m_iDeriv = 0;
	m_bACF = false;
	m_bACFFFT = false;
	m_sName = NULL;
	m_sShortName = NULL;
	m_sLabelName = NULL;
	m_faData = NULL;
	m_baDataEnabled = NULL;
}


CXDF::~CXDF()
{
}


void CXDF::ParseDeriv()
{
//	int ti;

	if (!g_bAdvanced2)
		m_iDeriv = 0;
			else m_iDeriv = AskUnsignedInteger("    Use values (0), 1st time derivative (1), or 2nd time derivative (2) of the values? [0] ",0);

	if (m_iDeriv != 0)
	{
		g_bDeriv = true;
		m_bDerivAbs = AskYesNo("    Take the absolute values of the derivatives (y/n)? [no] ",false);
		m_bACF = AskYesNo("    Autocorrelate derived values (y/n)? [no] ",false);
		if (m_bACF)
		{
			eprintf("\n    This feature is not working yet.\n\n");
			m_bACF = false;
		}
/*		if (m_bACF)
		{
			m_iACFDepth = AskUnsignedInteger("    Enter depth of the autocorrelation in time steps: [10000] ",10000);
			ti = CalcFFTSize(m_iACFDepth,false);
			if (m_iACFDepth != ti)
			{
				mprintf(WHITE,"\nThe next \"fast\" size for autocorrelation is %d. Using this instead of %d as depth.\n",ti,m_iACFDepth);
				m_iACFDepth = ti;
			}

			m_bACFFFT = AskYesNo("    Calculate the fourier transform of the autocorrelation (y/n)? [yes] ",true);
			if (m_bACFFFT)
			{
				m_bWindowFunction = AskYesNo("    Apply window function (Cos^2) to Autocorrelation function (y/n)? [yes] ",true);
				m_fSpecWaveNumber = AskFloat("    Calculate spectrum up to which wave number (cm^-1) (0=full spectrum)? [4000cm^-1] ",4000.0);
				m_iMirror = AskUnsignedInteger("    No mirroring (0), short-time enhancing (1) or long-time enhancing (2)? [0] ",0);
				m_iZeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ",m_iACFDepth*3,m_iACFDepth*3);

				ti = CalcFFTSize(m_iACFDepth+m_iZeroPadding,false);
				if (m_iACFDepth+m_iZeroPadding != ti)
				{
					mprintf(WHITE,"\nThe next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n",ti,ti-m_iACFDepth);
					m_iZeroPadding = ti-m_iACFDepth;
				}

				m_bACF_DB = AskYesNo("    Convert intensity axis of spectrum to decibel (y/n)? [no] ",false);
			}
		}*/
	}
}


void CXDF::InitDeriv()
{
	int z2;
	CxFloatArray *ptfa;

	try { m_pfaDerivBuffer = new CxDoubleArray*[3]; } catch(...) { m_pfaDerivBuffer = NULL; }
	if (m_pfaDerivBuffer == NULL) NewException((double)3*sizeof(CxDoubleArray*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z2=0;z2<3;z2++)
	{
		try { m_pfaDerivBuffer[z2] = new CxDoubleArray("CXDF::m_pfaDerivBuffer[z2]"); } catch(...) { m_pfaDerivBuffer[z2] = NULL; }
		if (m_pfaDerivBuffer[z2] == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		if (m_bSelf)
			m_pfaDerivBuffer[z2]->SetSize(((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*m_iCombinations);
				else m_pfaDerivBuffer[z2]->SetSize(((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*m_iCombinations);
	}

	if (m_bACF)
	{
		if (m_bSelf)
		{
			try { m_pfaACFBuffer = new CxFloatArray*[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*m_iCombinations]; } catch(...) { m_pfaACFBuffer = NULL; }
			if (m_pfaACFBuffer == NULL) NewException((double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*m_iCombinations*sizeof(CxFloatArray*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*m_iCombinations;z2++)
			{
				try { ptfa = new CxFloatArray("CXDF::InitDeriv():ptfa"); } catch(...) { ptfa = NULL; }
				if (ptfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				if (g_iTrajSteps != -1)
				{
					ptfa->SetMaxSize((long)(g_iTrajSteps*1.1));
					ptfa->SetGrow((long)(g_iTrajSteps*0.1));
				} else ptfa->SetGrow(1000);
				m_pfaACFBuffer[z2] = ptfa;
			}
		} else
		{
			try { m_pfaACFBuffer = new CxFloatArray*[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*m_iCombinations]; } catch(...) { m_pfaACFBuffer = NULL; }
			if (m_pfaACFBuffer == NULL) NewException((double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*m_iCombinations*sizeof(CxFloatArray*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*m_iCombinations;z2++)
			{
				try { ptfa = new CxFloatArray("CXDF::InitDeriv():ptfa"); } catch(...) { ptfa = NULL; }
				if (ptfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				if (g_iTrajSteps != -1)
				{
					ptfa->SetMaxSize((long)(g_iTrajSteps*1.1));
					ptfa->SetGrow((long)(g_iTrajSteps*0.1));
				} else ptfa->SetGrow(1000);
				m_pfaACFBuffer[z2] = ptfa;
			}
		}
	}
}


void CXDF::Autocorrelate()
{
	int n, z, z2;
	double tfs;
	CAutoCorrelation *ac;
	CxFloatArray *ptfa;

	if (m_bSelf)
		n = ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*m_iCombinations;
			else n = ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*m_iCombinations;

	tfs = n/60.0;
	mprintf("    Autocorrelating cached values...\n");
	mprintf(WHITE,"      [");

	try { ptfa = new CxFloatArray("CXDF::Autocorrelate():ptfa"); } catch(...) { ptfa = NULL; }
	if (ptfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	ptfa->SetSize(m_pfaACFBuffer[0]->GetSize());

	try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
	if (ac == NULL) NewException((double)sizeof(CAutoCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	ac->Init(m_pfaACFBuffer[0]->GetSize(),m_iACFDepth,true);
	m_faACF.SetSize(m_iACFDepth);
	for (z=0;z<m_iACFDepth;z++)
		m_faACF[z] = 0;
	for (z=0;z<n;z++)
	{
		if (fmod(z,tfs) < 1.0)
			mprintf(WHITE,"#");

		ac->AutoCorrelate(m_pfaACFBuffer[z],ptfa);
		for (z2=0;z2<m_iACFDepth;z2++)
			m_faACF[z2] += (double)(*ptfa)[z2];
	}
	delete ac;
	delete ptfa;
	for (z=1;z<m_iACFDepth;z++)
		m_faACF[z] /= m_faACF[0];
	m_faACF[0] = 1.0f;
	mprintf(WHITE,"]\n");
}
