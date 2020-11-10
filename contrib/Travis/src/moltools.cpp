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

#include "moltools.h"
#include "travis.h"
#include "maintools.h"
#include "plproj.h"



CAtom::CAtom()
{
	m_pElement = NULL;
	m_bExclude = false;
	m_pMergedTo = NULL;
	m_iIndex = -1;
}


CAtom::~CAtom()
{
}


CVirtualAtom::CVirtualAtom()
{
	m_faWeight.SetName("CVirtualAtom::m_faWeight");
}


CVirtualAtom::~CVirtualAtom()
{
}


CMolecule::CMolecule()
{
	m_iWannierCount = 0;
	m_sName = NULL;
	m_fCharge = 0;
	m_bPseudo = false;
	m_bChargesAssigned = false;
	m_bPolymer = false;

	m_baAtomIndex.SetName("CMolecule::m_baAtomIndex");
	m_waAtomCount.SetName("CMolecule::m_waAtomCount");
	m_laSingleMolIndex.SetName("CMolecule::m_laSingleMolIndex");
	m_laVirtualAtoms.SetName("CMolecule::m_laVirtualAtoms");
	m_oaNewNumbers.SetName("CMolecule::m_oaNewNumbers");
	m_oaRingAtomTypes.SetName("CMolecule::m_oaRingAtomTypes");
	m_oaRingAtoms.SetName("CMolecule::m_oaRingAtoms");
	m_oaCharges.SetName("CMolecule::m_oaCharges");
	
	m_iDipoleMode = 0;
	m_pDipoleFile = NULL;
}


CMolecule::~CMolecule()
{
	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}
	
	if(m_pDipoleFile != NULL)
		fclose(m_pDipoleFile);
}


CMolAtom::CMolAtom()
{
	m_oaBonds.SetName("CMolAtom::m_oaBonds");
}


CMolAtom::~CMolAtom()
{
}


CSingleMolecule::CSingleMolecule()
{
	m_bPseudo = false;

	m_oaBondGroups.SetName("CSingleMolecule::m_oaBondGroups");
	m_oaBonds.SetName("CSingleMolecule::m_oaBonds");
	m_oaAngleGroups.SetName("CSingleMolecule::m_oaAngleGroups");
	m_oaAngles.SetName("CSingleMolecule::m_oaAngles");
	m_oaRings.SetName("CSingleMolecule::m_oaRings");
	m_laBonds.SetName("CSingleMolecule::m_laBonds");
	m_laWannier.SetName("CSingleMolecule::m_laWannier");
	m_oaAtomOffset.SetName("CSingleMolecule::m_oaAtomOffset");
	m_oaMolAtoms.SetName("CSingleMolecule::m_oaMolAtoms");
	m_baAtomIndex.SetName("CSingleMolecule::m_baAtomIndex"); 
}


CSingleMolecule::~CSingleMolecule()
{
}


CADF::CADF()
{
	m_pADF = NULL;
	m_oaVectors.SetName("CADF::m_oaVectors");
	m_faACF.SetName("CADF::m_faACF");
	m_faMinMaxAngle.SetName("CADF::m_faMinMaxAngle");
}


CADF::~CADF()
{
}


CDDF::CDDF()
{
	m_bRotate = false;
	m_oaVectors.SetName("CDDF::m_oaVectors");
	m_faLastData.SetName("CDDF::m_faLastData"); 
	m_laRotation.SetName("CDDF::m_laRotation");
	m_faACF.SetName("CDDF::m_faACF");
}


CDDF::~CDDF()
{
}


CMSD::CMSD()
{
	m_oaCache.SetName("CMSD::m_oaCache");
}


CMSD::~CMSD()
{
}


CVHDF::CVHDF()
{
	m_sName = NULL;
	m_sShortName = NULL;
	m_pVHDF = NULL;
	m_oaVectors.SetName("CVHDF::m_oaVectors");
}


CVHDF::~CVHDF()
{
	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}
	if (m_sShortName != NULL)
	{
		delete[] m_sShortName;
		m_sShortName = NULL;
	}
	if (m_pVHDF != NULL)
	{
		delete m_pVHDF;
		m_pVHDF = NULL;
	}
}


CRDF::CRDF()
{
	m_fDist = NULL;
	m_pRDF = NULL;
	m_bProbDens = false;
	m_faMinMaxDist.SetName("CRDF::m_faMinMaxDist");
	m_oaVectors.SetName("CRDF::m_oaVectors");
	m_faACF.SetName("CRDF::m_faACF");
}


CPlDF::CPlDF()
{
	m_pPlDF = NULL;
	m_oaVectors.SetName("CPlDF::m_oaVectors");
	m_faACF.SetName("CPlDF::m_faACF");
}


CLiDF::CLiDF()
{
	m_pLiDF = NULL;
	m_oaVectors.SetName("CLiDF::m_oaVectors");
	m_faACF.SetName("CLiDF::m_faACF");
}


CRDF::~CRDF()
{
	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}
	if (m_sShortName != NULL)
	{
		delete[] m_sShortName;
		m_sShortName = NULL;
	}
	if (m_fDist != NULL)
	{
		delete[] m_fDist;
		m_fDist = NULL;
	}
	if (m_pRDF != NULL)
	{
		delete m_pRDF;
		m_pRDF = NULL;
	}
	if (m_faData != NULL)
	{
		delete[] m_faData;
		m_faData = NULL;
	}
}


CDensDF::CDensDF()
{
	m_pDensDF = NULL;
}


CDensDF::~CDensDF()
{
	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}
	if (m_sShortName != NULL)
	{
		delete[] m_sShortName;
		m_sShortName = NULL;
	}
	if (m_pDensDF != NULL)
	{
		delete m_pDensDF;
		m_pDensDF = NULL;
	}
}


CVDF::CVDF()
{
	m_faACF.SetName("CVDF::m_faACF");
}


CVDF::~CVDF()
{
}


CDipDF::CDipDF()
{
	m_faACF.SetName("CDipDF::m_faACF");
}


CDipDF::~CDipDF()
{
}


CSDF::CSDF()
{
	m_pCutPlane = NULL;
	m_fPosCounter = 0;
	m_fAtom2PosX = 0;
	m_fAtom3PosX = 0;
	m_fAtom3PosY = 0;
}


CSDF::~CSDF()
{
}


CCDF::CCDF()
{
	m_bAxisDivide = false;
}


CCDF::~CCDF()
{
}


CConditionSubGroup::CConditionSubGroup()
{
	m_fPassed = 0;
	m_fTotal = 0;
	m_oaConditions.SetName("CConditionSubGroup::m_oaConditions");
}


CConditionSubGroup::~CConditionSubGroup()
{
}


CConditionGroup::CConditionGroup()
{
	m_pTable = NULL;
	m_bInactive = false;
	m_fPassed = 0;
	m_fTotal = 0;
	m_iPassCounter = NULL;
	m_bAlwaysTrue = NULL;
	m_oaConditionSubGroups.SetName("CConditionGroup::m_oaConditionSubGroups");
}


CConditionGroup::~CConditionGroup()
{
}


bool ContainsDigit(const char *s)
{
	if (strcspn(s,"0123456789") != strlen(s))
		return true;
	return false;
}


/*
void ReplaceDigits(char *s)
{
	char uf[32];
	char *p, *q;
	bool b;

	p = s;
	q = buf;
	b = false;
	while (*p != 0)
	{
		if (b)
		{
			switch(*p)
			{
				case '1': *q = 'a'; break;
				case '2': *q = 'b'; break;
				case '3': *q = 'c'; break;
				case '4': *q = 'd'; break;
				case '5': *q = 'e'; break;
				case '6': *q = 'f'; break;
				case '7': *q = 'g'; break;
				case '8': *q = 'h'; break;
				case '9': *q = 'i'; break;
				case '0': *q = 'z'; break;
				default: *q = *p;
			}
		} else
		{
			switch(*p)
			{
				case '1': b = true; *q = '_'; q++; *q = 'a'; break;
				case '2': b = true; *q = '_'; q++; *q = 'b'; break;
				case '3': b = true; *q = '_'; q++; *q = 'c'; break;
				case '4': b = true; *q = '_'; q++; *q = 'd'; break;
				case '5': b = true; *q = '_'; q++; *q = 'e'; break;
				case '6': b = true; *q = '_'; q++; *q = 'f'; break;
				case '7': b = true; *q = '_'; q++; *q = 'g'; break;
				case '8': b = true; *q = '_'; q++; *q = 'h'; break;
				case '9': b = true; *q = '_'; q++; *q = 'i'; break;
				case '0': b = true; *q = '_'; q++; *q = 'z'; break;
				default: *q = *p;
			}
		}
		p++;
		q++;
	}
	*q = 0;
	strcpy(s,buf);
}
*/


void ReplaceDigits(CxString *s)
{
	char buf[32];
	char *p, *q;
	bool b;

	p = s->GetWritePointer();
	q = buf;
	b = false;
	while (*p != 0)
	{
		if (b)
		{
			switch(*p)
			{
				case '1': *q = 'a'; break;
				case '2': *q = 'b'; break;
				case '3': *q = 'c'; break;
				case '4': *q = 'd'; break;
				case '5': *q = 'e'; break;
				case '6': *q = 'f'; break;
				case '7': *q = 'g'; break;
				case '8': *q = 'h'; break;
				case '9': *q = 'i'; break;
				case '0': *q = 'z'; break;
				default: *q = *p;
			}
		} else
		{
			switch(*p)
			{
				case '1': b = true; *q = '_'; q++; *q = 'a'; break;
				case '2': b = true; *q = '_'; q++; *q = 'b'; break;
				case '3': b = true; *q = '_'; q++; *q = 'c'; break;
				case '4': b = true; *q = '_'; q++; *q = 'd'; break;
				case '5': b = true; *q = '_'; q++; *q = 'e'; break;
				case '6': b = true; *q = '_'; q++; *q = 'f'; break;
				case '7': b = true; *q = '_'; q++; *q = 'g'; break;
				case '8': b = true; *q = '_'; q++; *q = 'h'; break;
				case '9': b = true; *q = '_'; q++; *q = 'i'; break;
				case '0': b = true; *q = '_'; q++; *q = 'z'; break;
				default: *q = *p;
			}
		}
		p++;
		q++;
	}
	*q = 0;
	s->strcpy(buf);
}


void xAddAtom(const char *s)
{
	BTIN;
	int z;
	CAtom *a;
//	char buf[64];
	CxString buf;

	if (g_oaAtoms.GetSize() >= 254)
	{
		eprintf("More than 254 different atom types not supported.\n");
		return;
	}

//	strcpy(buf,s);
	buf.strcpy(s);

	ReplaceDigits(&buf);

//	printf("AddAtom: \"%s\".\n",s);
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		if (mystricmp(buf,((CAtom*)g_oaAtoms[z])->m_sName)==0)
		{
			((CAtom*)g_oaAtoms[z])->m_iCount++;
			BTOUT; 
			return;
		}
	}
	if (ContainsDigit(s))
		eprintf("Digits in element labels not allowed. Renaming %s to %s.\n",s,(const char*)buf);

	try { a = new CAtom(); } catch(...) { a = NULL; }
	if (a == NULL) NewException((double)sizeof(CAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	a->m_iIndex = g_oaAtoms.GetSize();
	strcpy(a->m_sName,buf);
	if (islower(a->m_sName[0]))
		a->m_sName[0] = toupper(a->m_sName[0]);
	if (strlen(a->m_sName) > 1)
		if (isupper(a->m_sName[1]))
			a->m_sName[1] = tolower(a->m_sName[1]);
	a->m_iCount = 1;
/*	if (s[0] != '#')
	{*/
		a->m_pElement = FindElement(buf,false);
		if (a->m_pElement == NULL)
		{
			try { a->m_pElement = new CElement(); } catch(...) { a->m_pElement = NULL; }
			if (a->m_pElement == NULL) NewException((double)sizeof(CElement),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			try { a->m_pElement->m_sLabel = new char[strlen(a->m_sName)+1]; } catch(...) { a->m_pElement->m_sLabel = NULL; }
			if (a->m_pElement->m_sLabel == NULL) NewException((double)(strlen(a->m_sName)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(a->m_pElement->m_sLabel,a->m_sName);
		}
/*		a->m_fMass = AtomMass(s);
		a->m_fRadius = AtomRadius(s);
		a->m_iOrd = AtomOrd(s);*/
//		a->m_fVDWRadius = AtomVDWRadius(s);
/*	} else
	{
		a->m_fMass = 0.0f;
		a->m_fRadius = 0.0f;
		a->m_iOrd = 0;
//		a->m_fVDWRadius = 0.0f;
	}*/
	g_oaAtoms.Add(a);
//	g_pAtoms[g_iElementCount].Offset = offset;
//	printf("Fuege Atom %s an Stelle %d neu hinzu. Der Offset ist %d.\n",s,g_iAtomCount,offset);
//	g_iElementCount++;
	BTOUT; 
}


void CSingleMolecule::Dump()
{
	BTIN;
	int z, z2;
	mprintf("### Single Molecule Dump ###\n");
	mprintf("%d Elemente.\n",m_baAtomIndex.GetSize());
	for (z=0;z<m_baAtomIndex.GetSize();z++)
	{
		mprintf(" * Element %d: %s. %d Vertreter *\n    Atome ",z+1,((CAtom*)g_oaAtoms[m_baAtomIndex[z]])->m_sName,((CxIntArray*)m_oaAtomOffset[z])->GetSize());
		for (z2=0;z2<((CxIntArray*)m_oaAtomOffset[z])->GetSize();z2++)
		{
			mprintf("%d",((CxIntArray*)m_oaAtomOffset[z])->GetAt(z2));
			if (z2 < ((CxIntArray*)m_oaAtomOffset[z])->GetSize()-1)
				mprintf(", ");
		}
		mprintf("\n");
	}
	BTOUT;
}


void CMolecule::Dump()
{
	BTIN;
	int z;
	mprintf("### Molecule Type Dump ###\n");
	mprintf("%d Elemente.\n",m_baAtomIndex.GetSize());
	for (z=0;z<m_baAtomIndex.GetSize();z++)
		mprintf(" * Element %d: %s. %d Vertreter *\n",z+1,((CAtom*)g_oaAtoms[m_baAtomIndex[z]])->m_sName,m_waAtomCount[z]);
	BTOUT;
}


void CADF::BuildName()
{
	BTIN;
	int z, z2;
//	char tmp[32768];
	CxString tmp;
	CAtomGroup *ag;

//	tmp[0] = 0;
	tmp.sprintf("");

	if (m_iDeriv != 0)
//		sprintf(tmp,"deriv%d_",m_iDeriv);
		tmp.sprintf("deriv%d_",m_iDeriv);
	for (z2=0;z2<m_oaVectors.GetSize()/6;z2++)
	{
		for (z=0;z<2;z++)
		{
//			strcat(tmp,"[");
			tmp.strcat("[");
			if (m_iVecType[z] == 0) // Position
			{
				if (m_bOrtho[z])
				{
/*					if (m_bSameFoot && (z == 1))
						ag = (CAtomGroup*)m_oaVectors[z2*6];
							else */ag = (CAtomGroup*)m_oaVectors[z2*6+z*3];

//					strcat(tmp,ag->m_sName);
					tmp.strcat(ag->m_sName);

					if (m_iRefOrSec[z][0])
//						strcat(tmp,"o_");
						tmp.strcat("o_");
					else
//						strcat(tmp,"r_");
						tmp.strcat("r_");

					tmp.strcat(((CAtomGroup*)m_oaVectors[z2*6+z*3+1])->m_sName);

					if (m_iRefOrSec[z][1])
//						strcat(tmp,"o_");
						tmp.strcat("o_");
					else
//						strcat(tmp,"r_");
						tmp.strcat("r_");

//					strcat(tmp,((CAtomGroup*)m_oaVectors[z2*6+z*3+2])->m_sName);
					tmp.strcat(((CAtomGroup*)m_oaVectors[z2*6+z*3+2])->m_sName);

					if (m_iRefOrSec[z][2])
//						strcat(tmp,"o");
						tmp.strcat("o");
					else
//						strcat(tmp,"r");
						tmp.strcat("r");

				} else
				{
			/*		if (m_bSameFoot && (z == 1))
						ag = (CAtomGroup*)m_oaVectors[z2*6];
							else */ag = (CAtomGroup*)m_oaVectors[z2*6+z*3];

//					strcat(tmp,ag->m_sName);
					tmp.strcat(ag->m_sName);

					if (m_iRefOrSec[z][0])
//						strcat(tmp,"o_");
						tmp.strcat("o_");
					else
//						strcat(tmp,"r_");
						tmp.strcat("r_");

//					strcat(tmp,((CAtomGroup*)m_oaVectors[z2*6+z*3+1])->m_sName);
					tmp.strcat(((CAtomGroup*)m_oaVectors[z2*6+z*3+1])->m_sName);

					if (m_iRefOrSec[z][1])
//						strcat(tmp,"o");
						tmp.strcat("o");
					else
//						strcat(tmp,"r");
						tmp.strcat("r");
				}
			} else if (m_iVecType[z] == 1) // Dipol
			{
//				strcat(tmp,"dip_");
//				strcat(tmp,(m_iRefOrSec[z][0]!=0)?((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName:((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
				tmp.strcat("dip_");
				tmp.strcat((m_iRefOrSec[z][0]!=0)?((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName:((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
			} else if (m_iVecType[z] == 2) // Geschwindigkeit
			{
//				strcat(tmp,"vel_");
//				strcat(tmp,((CAtomGroup*)m_oaVectors[z*3])->m_sName);
				tmp.strcat("vel_");
				tmp.strcat(((CAtomGroup*)m_oaVectors[z*3])->m_sName);
			} else if (m_iVecType[z] == 3) // Kraft
			{
//				strcat(tmp,"frc_");
//				strcat(tmp,((CAtomGroup*)m_oaVectors[z*3])->m_sName);
				tmp.strcat("frc_");
				tmp.strcat(((CAtomGroup*)m_oaVectors[z*3])->m_sName);
			}
			if (z == 0)
//				strcat(tmp,"]-");
				tmp.strcat("]-");
			else
//				strcat(tmp,"]");
				tmp.strcat("]");
		}
		if (z2<(m_oaVectors.GetSize()/6)-1)
//			strcat(tmp,"_");
			tmp.strcat("_");
	}

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sShortName,tmp);

	if (m_iShowMol != -1)
//		sprintf(tmp,"%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.sprintf("%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	else
//		sprintf(tmp,"%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);
		tmp.sprintf("%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);

//	strcat(tmp,m_sShortName);
	tmp.strcat(m_sShortName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);


	tmp.sprintf("");
	if (m_iDeriv != 0)
		tmp.sprintf("deriv%d_",m_iDeriv);
	for (z=0;z<2;z++)
	{
		tmp.strcat("[");
		if (m_iVecType[z] == 0) // Position
		{
			if (m_bOrtho[z])
			{
				ag = (CAtomGroup*)m_oaVectors[z*3];
				tmp.strcat(ag->m_sName);
				if (m_iRefOrSec[z][0])
					tmp.strcat("o_");
				else
					tmp.strcat("r_");
				tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+1])->m_sName);
				if (m_iRefOrSec[z][1])
					tmp.strcat("o_");
				else
					tmp.strcat("r_");
				tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+2])->m_sName);
				if (m_iRefOrSec[z][2])
					tmp.strcat("o");
				else
					tmp.strcat("r");
			} else
			{
				ag = (CAtomGroup*)m_oaVectors[z*3];
				tmp.strcat(ag->m_sName);
				if (m_iRefOrSec[z][0])
					tmp.strcat("o_");
				else
					tmp.strcat("r_");
				tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+1])->m_sName);
				if (m_iRefOrSec[z][1])
					tmp.strcat("o");
				else
					tmp.strcat("r");
			}
		} else if (m_iVecType[z] == 1) // Dipol
		{
			tmp.strcat("dip_");
			tmp.strcat((m_iRefOrSec[z][0]!=0)?((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName:((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		} else if (m_iVecType[z] == 2) // Geschwindigkeit
		{
			tmp.strcat("vel_");
			tmp.strcat(((CAtomGroup*)m_oaVectors[z*3])->m_sName);
		} else if (m_iVecType[z] == 3) // Kraft
		{
			tmp.strcat("frc_");
			tmp.strcat(((CAtomGroup*)m_oaVectors[z*3])->m_sName);
		}
		if (z == 0)
			tmp.strcat("]-");
		else
			tmp.strcat("]");
	}
	if (m_oaVectors.GetSize() > 6)
		tmp.strcat(",(...)");
	try { m_sLabelName = new char[strlen(tmp)+1]; } catch(...) { m_sLabelName = NULL; }
	if (m_sLabelName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sLabelName,tmp);

	BTOUT;
}


void CPlDF::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;
	CAtomGroup *ag;

//	tmp[0] = 0;
	tmp.sprintf("");

//	strcat(tmp,"[");
	tmp.strcat("[");
	if (m_bNormal)
	{
		ag = (CAtomGroup*)m_oaVectors[0];
//		strcat(tmp,ag->m_sName);
		tmp.strcat(ag->m_sName);

		if (m_iRefOrSec[0])
//			strcat(tmp,"o_");
			tmp.strcat("o_");
		else
//			strcat(tmp,"r_");
			tmp.strcat("r_");

//		strcat(tmp,((CAtomGroup*)m_oaVectors[1])->m_sName);
		tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);

		if (m_iRefOrSec[1])
//			strcat(tmp,"o");
			tmp.strcat("o");
		else
//			strcat(tmp,"r");
			tmp.strcat("r");
	} else
	{
		ag = (CAtomGroup*)m_oaVectors[0];
//		strcat(tmp,ag->m_sName);
		tmp.strcat(ag->m_sName);

		if (m_iRefOrSec[0])
//			strcat(tmp,"o_");
			tmp.strcat("o_");
		else
//			strcat(tmp,"r_");
			tmp.strcat("r_");

//		strcat(tmp,((CAtomGroup*)m_oaVectors[1])->m_sName);
		tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);

		if (m_iRefOrSec[1])
//			strcat(tmp,"o_");
			tmp.strcat("o_");
		else
//			strcat(tmp,"r_");
			tmp.strcat("r_");

//		strcat(tmp,((CAtomGroup*)m_oaVectors[2])->m_sName);
		tmp.strcat(((CAtomGroup*)m_oaVectors[2])->m_sName);

		if (m_iRefOrSec[2])
//			strcat(tmp,"o");
			tmp.strcat("o");
		else
//			strcat(tmp,"r");
			tmp.strcat("r");
	}

//	strcat(tmp,"]_");
//	strcat(tmp,((CAtomGroup*)m_oaVectors[3])->m_sName);
	tmp.strcat("]_");
	tmp.strcat(((CAtomGroup*)m_oaVectors[3])->m_sName);

	if (m_iRefOrSec[3])
//		strcat(tmp,"o");
		tmp.strcat("o");
	else
//		strcat(tmp,"r");
		tmp.strcat("r");

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sShortName,tmp);

	if (m_iShowMol != -1)
//		sprintf(tmp,"%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.sprintf("%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	else
//		sprintf(tmp,"%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);
		tmp.sprintf("%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);

//	strcat(tmp,m_sShortName);
	tmp.strcat(m_sShortName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);


	tmp.sprintf("[");
	if (m_bNormal)
	{
		ag = (CAtomGroup*)m_oaVectors[0];
		tmp.strcat(ag->m_sName);
		if (m_iRefOrSec[0])
			tmp.strcat("o_");
		else
			tmp.strcat("r_");
		tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);
		if (m_iRefOrSec[1])
			tmp.strcat("o");
		else
			tmp.strcat("r");
	} else
	{
		ag = (CAtomGroup*)m_oaVectors[0];
		tmp.strcat(ag->m_sName);
		if (m_iRefOrSec[0])
			tmp.strcat("o_");
		else
			tmp.strcat("r_");
		tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);
		if (m_iRefOrSec[1])
			tmp.strcat("o_");
		else
			tmp.strcat("r_");
		tmp.strcat(((CAtomGroup*)m_oaVectors[2])->m_sName);
		if (m_iRefOrSec[2])
			tmp.strcat("o");
		else
			tmp.strcat("r");
	}
	tmp.strcat("]_");
	tmp.strcat(((CAtomGroup*)m_oaVectors[3])->m_sName);
	if (m_iRefOrSec[3])
		tmp.strcat("o");
	else
		tmp.strcat("r");
	if (m_oaVectors.GetSize() > 4)
		tmp.strcat(",(...)");
	try { m_sLabelName = new char[strlen(tmp)+1]; } catch(...) { m_sLabelName = NULL; }
	if (m_sLabelName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sLabelName,tmp);

	BTOUT;
}


void CLiDF::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;
	CAtomGroup *ag;

//	tmp[0] = 0;
	tmp.sprintf("");

//	strcat(tmp,"[");
	tmp.strcat("[");
	if (m_bNormal)
	{
		ag = (CAtomGroup*)m_oaVectors[0];
//		strcat(tmp,ag->m_sName);
		tmp.strcat(ag->m_sName);

		if (m_iRefOrSec[0])
//			strcat(tmp,"o_");
			tmp.strcat("o_");
		else
//			strcat(tmp,"r_");
			tmp.strcat("r_");

//		strcat(tmp,((CAtomGroup*)m_oaVectors[1])->m_sName);
		tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);

		if (m_iRefOrSec[1])
//			strcat(tmp,"o_");
			tmp.strcat("o_");
		else
//			strcat(tmp,"r_");
			tmp.strcat("r_");

//		strcat(tmp,((CAtomGroup*)m_oaVectors[2])->m_sName);
		tmp.strcat(((CAtomGroup*)m_oaVectors[2])->m_sName);

		if (m_iRefOrSec[2])
//			strcat(tmp,"o");
			tmp.strcat("o");
		else
//			strcat(tmp,"r");
			tmp.strcat("r");
	} else
	{
		ag = (CAtomGroup*)m_oaVectors[0];
//		strcat(tmp,ag->m_sName);
		tmp.strcat(ag->m_sName);

		if (m_iRefOrSec[0])
//			strcat(tmp,"o_");
			tmp.strcat("o_");
		else
//			strcat(tmp,"r_");
			tmp.strcat("r_");

//		strcat(tmp,((CAtomGroup*)m_oaVectors[1])->m_sName);
		tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);

		if (m_iRefOrSec[1])
//			strcat(tmp,"o");
			tmp.strcat("o");
		else
//			strcat(tmp,"r");
			tmp.strcat("r");
	}
//	strcat(tmp,"]_");
//	strcat(tmp,((CAtomGroup*)m_oaVectors[3])->m_sName);
	tmp.strcat("]_");
	tmp.strcat(((CAtomGroup*)m_oaVectors[3])->m_sName);

	if (m_iRefOrSec[3])
//		strcat(tmp,"o");
		tmp.strcat("o");
	else
//		strcat(tmp,"r");
		tmp.strcat("r");

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sShortName,tmp);

	if (m_iShowMol != -1)
//		sprintf(tmp,"%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.sprintf("%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	else
//		sprintf(tmp,"%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);
		tmp.sprintf("%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);

//	strcat(tmp,m_sShortName);
	tmp.strcat(m_sShortName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);


	tmp.sprintf("[");
	if (m_bNormal)
	{
		ag = (CAtomGroup*)m_oaVectors[0];
		tmp.strcat(ag->m_sName);
		if (m_iRefOrSec[0])
			tmp.strcat("o_");
		else
			tmp.strcat("r_");
		tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);
		if (m_iRefOrSec[1])
			tmp.strcat("o_");
		else
			tmp.strcat("r_");
		tmp.strcat(((CAtomGroup*)m_oaVectors[2])->m_sName);
		if (m_iRefOrSec[2])
			tmp.strcat("o");
		else
			tmp.strcat("r");
	} else
	{
		ag = (CAtomGroup*)m_oaVectors[0];
		tmp.strcat(ag->m_sName);
		if (m_iRefOrSec[0])
			tmp.strcat("o_");
		else
			tmp.strcat("r_");
		tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);
		if (m_iRefOrSec[1])
			tmp.strcat("o");
		else
			tmp.strcat("r");
	}
	tmp.strcat("]_");
	tmp.strcat(((CAtomGroup*)m_oaVectors[3])->m_sName);
	if (m_iRefOrSec[3])
		tmp.strcat("o");
	else
		tmp.strcat("r");
	if (m_oaVectors.GetSize() > 4)
		tmp.strcat(",(...)");
	try { m_sLabelName = new char[strlen(tmp)+1]; } catch(...) { m_sLabelName = NULL; }
	if (m_sLabelName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);	
	strcpy(m_sLabelName,tmp);


	BTOUT;
}


void CDDF::BuildName()
{
	BTIN;
	int z;
//	char tmp[32768];
	CxString tmp;
	CAtomGroup *ag;

//	tmp[0] = 0;
	tmp.sprintf("");

	if (m_iDeriv != 0)
//		sprintf(tmp,"deriv%d_",m_iDeriv);
		tmp.sprintf("deriv%d_",m_iDeriv);

	for (z=0;z<3;z++)
	{
//		strcat(tmp,"[");
		tmp.strcat("[");
		if (m_bOrtho[z])
		{
			ag = (CAtomGroup*)m_oaVectors[z*3];
//			strcat(tmp,ag->m_sName);
			tmp.strcat(ag->m_sName);

			if (m_iRefOrSec[z][0])
//				strcat(tmp,"o_");
				tmp.strcat("o_");
			else
//				strcat(tmp,"r_");
				tmp.strcat("r_");

//			strcat(tmp,((CAtomGroup*)m_oaVectors[z*3+1])->m_sName);
			tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+1])->m_sName);

			if (m_iRefOrSec[z][1])
//				strcat(tmp,"o_");
				tmp.strcat("o_");
			else
//				strcat(tmp,"r_");
				tmp.strcat("r_");

//			strcat(tmp,((CAtomGroup*)m_oaVectors[z*3+2])->m_sName);
			tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+2])->m_sName);

			if (m_iRefOrSec[z][2])
//				strcat(tmp,"o");
				tmp.strcat("o");
			else
//				strcat(tmp,"r");
				tmp.strcat("r");
		} else
		{
			ag = (CAtomGroup*)m_oaVectors[z*3];
//			strcat(tmp,ag->m_sName);
			tmp.strcat(ag->m_sName);

			if (m_iRefOrSec[z][0])
//				strcat(tmp,"o_");
				tmp.strcat("o_");
			else
//				strcat(tmp,"r_");
				tmp.strcat("r_");

//			strcat(tmp,((CAtomGroup*)m_oaVectors[z*3+1])->m_sName);
			tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+1])->m_sName);

			if (m_iRefOrSec[z][1])
//				strcat(tmp,"o");
				tmp.strcat("o");
			else
//				strcat(tmp,"r");
				tmp.strcat("r");
		}

		if (z < 2)
//			strcat(tmp,"]-");
			tmp.strcat("]-");
		else
//			strcat(tmp,"]");
			tmp.strcat("]");
	}

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sShortName,tmp);

	if (m_iShowMol != -1)
//		sprintf(tmp,"%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.sprintf("%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	else
//		sprintf(tmp,"%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);
		tmp.sprintf("%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);

//	strcat(tmp,m_sShortName);
	tmp.strcat(m_sShortName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);


	tmp.sprintf("");
	if (m_iDeriv != 0)
		tmp.sprintf("deriv%d_",m_iDeriv);
	for (z=0;z<3;z++)
	{
		tmp.strcat("[");
		if (m_bOrtho[z])
		{
			ag = (CAtomGroup*)m_oaVectors[z*3];
			tmp.strcat(ag->m_sName);
			if (m_iRefOrSec[z][0])
				tmp.strcat("o_");
			else
				tmp.strcat("r_");
			tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+1])->m_sName);
			if (m_iRefOrSec[z][1])
				tmp.strcat("o_");
			else
				tmp.strcat("r_");
			tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+2])->m_sName);
			if (m_iRefOrSec[z][2])
				tmp.strcat("o");
			else
				tmp.strcat("r");
		} else
		{
			ag = (CAtomGroup*)m_oaVectors[z*3];
			tmp.strcat(ag->m_sName);
			if (m_iRefOrSec[z][0])
				tmp.strcat("o_");
			else
				tmp.strcat("r_");
			tmp.strcat(((CAtomGroup*)m_oaVectors[z*3+1])->m_sName);
			if (m_iRefOrSec[z][1])
				tmp.strcat("o");
			else
				tmp.strcat("r");
		}
		if (z < 2)
			tmp.strcat("]-");
		else
			tmp.strcat("]");
	}
	if (m_oaVectors.GetSize() > 9)
		tmp.strcat(",(...)");
	try { m_sLabelName = new char[strlen(tmp)+1]; } catch(...) { m_sLabelName = NULL; }
	if (m_sLabelName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sLabelName,tmp);

	BTOUT;
}


void CDipDF::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;

//	tmp[0] = 0;
	tmp.sprintf("");

	if (m_iDeriv != 0)
//		sprintf(tmp,"deriv%d_",m_iDeriv);
		tmp.sprintf("deriv%d_",m_iDeriv);

	if (m_iRefOrSec)
//		sprintf(tmp,"%s",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.sprintf("%s",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	else
//		sprintf(tmp,"%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		tmp.sprintf("%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,tmp);

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sShortName,tmp);

	try { m_sLabelName = new char[strlen(tmp)+1]; } catch(...) { m_sLabelName = NULL; }
	if (m_sLabelName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sLabelName,tmp);

	BTOUT;
}


void CADF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z, z1t, z1a, z2t, z2a, z3t, z3a, z4t, z4a, z5t, z5a, z6t, z6a;
	CAtomGroup *g1, *g2, *g3, *g4, *g5, *g6;
	CxIntArray *a1, *a2, *a3, *a4, *a5, *a6;

//	mprintf("Vec: %X\n",vec);

	vec->RemoveAll_KeepSize();
	for (z=0;z<m_oaVectors.GetSize()/6;z++)
	{
//		mprintf("[%d]Vec: %X\n",z,vec);
		if (m_bOrtho[0])
		{
			g1 = (CAtomGroup*)m_oaVectors[z*6];
			g2 = (CAtomGroup*)m_oaVectors[z*6+1];
			g3 = (CAtomGroup*)m_oaVectors[z*6+2];
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
									if (m_bOrtho[1])
									{
										g4 = (CAtomGroup*)m_oaVectors[z*6+3];
										g5 = (CAtomGroup*)m_oaVectors[z*6+4];
										g6 = (CAtomGroup*)m_oaVectors[z*6+5];
										for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
										{
											a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
											for (z4a=0;z4a<a4->GetSize();z4a++)
											{
												for (z5t=0;z5t<g5->m_baAtomType.GetSize();z5t++)
												{
													a5 = (CxIntArray*)g5->m_oaAtoms[z5t];
													for (z5a=0;z5a<a5->GetSize();z5a++)
													{
														for (z6t=0;z6t<g6->m_baAtomType.GetSize();z6t++)
														{
															a6 = (CxIntArray*)g6->m_oaAtoms[z6t];
															for (z6a=0;z6a<a6->GetSize();z6a++)
															{
																if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																	vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																else
																	vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																	vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																else
																	vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																if ((!m_iRefOrSec[0][2]) || (obs == NULL))
																	vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
																else
																	vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));

																if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																	vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																else
																	vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																	vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																else
																	vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																if ((!m_iRefOrSec[1][2]) || (obs == NULL))
																	vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));
																else
																	vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));
															}
														}
													}
												}
											}
										}
									} else
									{
										g4 = (CAtomGroup*)m_oaVectors[z*6+3];
										g5 = (CAtomGroup*)m_oaVectors[z*6+4];
										for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
										{
											a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
											for (z4a=0;z4a<a4->GetSize();z4a++)
											{
												for (z5t=0;z5t<g5->m_baAtomType.GetSize();z5t++)
												{
													a5 = (CxIntArray*)g5->m_oaAtoms[z5t];
													for (z5a=0;z5a<a5->GetSize();z5a++)
													{
														if ((!m_iRefOrSec[0][0]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

														if ((!m_iRefOrSec[0][1]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

														if ((!m_iRefOrSec[0][2]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));

														if ((!m_iRefOrSec[1][0]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

														if ((!m_iRefOrSec[1][1]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

														vec->Add(0);
													}
												}
											}
										}
									} // END IF NOT ORTHO[1]
								}
							}
						}
					}
				}
			}
		} else
		{
			g1 = (CAtomGroup*)m_oaVectors[z*6];
			g2 = (CAtomGroup*)m_oaVectors[z*6+1];
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
							if (m_bOrtho[1])
							{
								g4 = (CAtomGroup*)m_oaVectors[z*6+3];
								g5 = (CAtomGroup*)m_oaVectors[z*6+4];
								g6 = (CAtomGroup*)m_oaVectors[z*6+5];
								for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
								{
									a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
									for (z4a=0;z4a<a4->GetSize();z4a++)
									{
										for (z5t=0;z5t<g5->m_baAtomType.GetSize();z5t++)
										{
											a5 = (CxIntArray*)g5->m_oaAtoms[z5t];
											for (z5a=0;z5a<a5->GetSize();z5a++)
											{
												for (z6t=0;z6t<g6->m_baAtomType.GetSize();z6t++)
												{
													a6 = (CxIntArray*)g6->m_oaAtoms[z6t];
													for (z6a=0;z6a<a6->GetSize();z6a++)
													{
														if ((!m_iRefOrSec[0][0]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

														if ((!m_iRefOrSec[0][1]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

														vec->Add(0);

														if ((!m_iRefOrSec[1][0]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

														if ((!m_iRefOrSec[1][1]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

														if ((!m_iRefOrSec[1][2]) || (obs == NULL))
															vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));
														else
															vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));
													}
												}
											}
										}
									}
								}
							} else
							{
								g4 = (CAtomGroup*)m_oaVectors[z*6+3];
								g5 = (CAtomGroup*)m_oaVectors[z*6+4];
								for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
								{
									a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
									for (z4a=0;z4a<a4->GetSize();z4a++)
									{
										for (z5t=0;z5t<g5->m_baAtomType.GetSize();z5t++)
										{
											a5 = (CxIntArray*)g5->m_oaAtoms[z5t];
											for (z5a=0;z5a<a5->GetSize();z5a++)
											{
											//	mprintf("### BuildAtomList\n");
											//	mprintf("    z1a=%d, z1t=%d, z2a=%d, z2t=%d, z4a=%d, z4t=%d, z5a=%d, z5t=%d.\n",z1a,z1t,z2a,z2t,z4a,z4t,z5a,z5t);
											//	mprintf("    g1=\"%s\", g2=\"%s\", g4=\"%s\", g5=\"%s\".\n",g1->m_sName,g2->m_sName,g4->m_sName,g5->m_sName);

												if ((!m_iRefOrSec[0][0]) || (obs == NULL))
													vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
												else
													vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

												if ((!m_iRefOrSec[0][1]) || (obs == NULL))
													vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
												else
													vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

												vec->Add(0);

												if ((!m_iRefOrSec[1][0]) || (obs == NULL))
													vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
												else
													vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

												if ((!m_iRefOrSec[1][1]) || (obs == NULL))
													vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
												else
													vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

												vec->Add(0);
											}
										}
									}
								}
							} // END IF NOT ORTHO[1]
						}
					}
				}
			}
		}
	}
	BXOUT;
}


/* Mega abartig !!!! */
void CDDF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z, z1t, z1a, z2t, z2a, z3t, z3a, z4t, z4a, z5t, z5a, z6t, z6a, z7t, z7a, z8t, z8a, z9t, z9a;
	CAtomGroup *g1, *g2, *g3, *g4, *g5, *g6, *g7, *g8, *g9;
	CxIntArray *a1, *a2, *a3, *a4, *a5, *a6, *a7, *a8, *a9;

	vec->RemoveAll_KeepSize();
	for (z=0;z<m_oaVectors.GetSize()/9;z++)
	{
		if (m_bOrtho[0])
		{
			g1 = (CAtomGroup*)m_oaVectors[z*9];
			g2 = (CAtomGroup*)m_oaVectors[z*9+1];
			g3 = (CAtomGroup*)m_oaVectors[z*9+2];
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
									if (m_bOrtho[1])
									{
										g4 = (CAtomGroup*)m_oaVectors[z*9+3];
										g5 = (CAtomGroup*)m_oaVectors[z*9+4];
										g6 = (CAtomGroup*)m_oaVectors[z*9+5];
										for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
										{
											a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
											for (z4a=0;z4a<a4->GetSize();z4a++)
											{
												for (z5t=0;z5t<g5->m_baAtomType.GetSize();z5t++)
												{
													a5 = (CxIntArray*)g5->m_oaAtoms[z5t];
													for (z5a=0;z5a<a5->GetSize();z5a++)
													{
														for (z6t=0;z6t<g6->m_baAtomType.GetSize();z6t++)
														{
															a6 = (CxIntArray*)g6->m_oaAtoms[z6t];
															for (z6a=0;z6a<a6->GetSize();z6a++)
															{
																if (m_bOrtho[2])
																{
																	g7 = (CAtomGroup*)m_oaVectors[z*9+6];
																	g8 = (CAtomGroup*)m_oaVectors[z*9+7];
																	g9 = (CAtomGroup*)m_oaVectors[z*9+8];
																	for (z7t=0;z7t<g7->m_baAtomType.GetSize();z7t++)
																	{
																		a7 = (CxIntArray*)g7->m_oaAtoms[z7t];
																		for (z7a=0;z7a<a7->GetSize();z7a++)
																		{
																			for (z8t=0;z8t<g8->m_baAtomType.GetSize();z8t++)
																			{
																				a8 = (CxIntArray*)g8->m_oaAtoms[z8t];
																				for (z8a=0;z8a<a8->GetSize();z8a++)
																				{
																					for (z9t=0;z9t<g9->m_baAtomType.GetSize();z9t++)
																					{
																						a9 = (CxIntArray*)g9->m_oaAtoms[z9t];
																						for (z9a=0;z9a<a9->GetSize();z9a++)
																						{
																							if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																							if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																							if ((!m_iRefOrSec[0][2]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));

																							if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																							if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																							if ((!m_iRefOrSec[1][2]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));

																							if ((!m_iRefOrSec[2][0]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));

																							if ((!m_iRefOrSec[2][1]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));

																							if ((!m_iRefOrSec[2][2]) || (obs == NULL))
																								vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g9->m_baAtomType[z9t]])->GetAt(a9->GetAt(z9a)));
																							else
																								vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g9->m_baAtomType[z9t]])->GetAt(a9->GetAt(z9a)));
																						}
																					}
																				}
																			}
																		}
																	}
																} else // IF NOT ORTHO[2]
																{
																	g7 = (CAtomGroup*)m_oaVectors[z*9+6];
																	g8 = (CAtomGroup*)m_oaVectors[z*9+7];
																	for (z7t=0;z7t<g7->m_baAtomType.GetSize();z7t++)
																	{
																		a7 = (CxIntArray*)g7->m_oaAtoms[z7t];
																		for (z7a=0;z7a<a7->GetSize();z7a++)
																		{
																			for (z8t=0;z8t<g8->m_baAtomType.GetSize();z8t++)
																			{
																				a8 = (CxIntArray*)g8->m_oaAtoms[z8t];
																				for (z8a=0;z8a<a8->GetSize();z8a++)
																				{
																					if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																					if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																					if ((!m_iRefOrSec[0][2]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));

																					if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																					if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																					if ((!m_iRefOrSec[1][2]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));

																					if ((!m_iRefOrSec[2][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));

																					if ((!m_iRefOrSec[2][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));

																					vec->Add(0);
																				}
																			}
																		}
																	}
																} // END IF NOT ORTHO[2]
															}
														}
													}
												}
											}
										}
									} else // IF NOT ORTHO[1]
									{
										g4 = (CAtomGroup*)m_oaVectors[z*9+3];
										g5 = (CAtomGroup*)m_oaVectors[z*9+4];
										for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
										{
											a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
											for (z4a=0;z4a<a4->GetSize();z4a++)
											{
												for (z5t=0;z5t<g5->m_baAtomType.GetSize();z5t++)
												{
													a5 = (CxIntArray*)g5->m_oaAtoms[z5t];
													for (z5a=0;z5a<a5->GetSize();z5a++)
													{
														if (m_bOrtho[2])
														{
															g7 = (CAtomGroup*)m_oaVectors[z*9+6];
															g8 = (CAtomGroup*)m_oaVectors[z*9+7];
															g9 = (CAtomGroup*)m_oaVectors[z*9+8];
															for (z7t=0;z7t<g7->m_baAtomType.GetSize();z7t++)
															{
																a7 = (CxIntArray*)g7->m_oaAtoms[z7t];
																for (z7a=0;z7a<a7->GetSize();z7a++)
																{
																	for (z8t=0;z8t<g8->m_baAtomType.GetSize();z8t++)
																	{
																		a8 = (CxIntArray*)g8->m_oaAtoms[z8t];
																		for (z8a=0;z8a<a8->GetSize();z8a++)
																		{
																			for (z9t=0;z9t<g9->m_baAtomType.GetSize();z9t++)
																			{
																				a9 = (CxIntArray*)g9->m_oaAtoms[z9t];
																				for (z9a=0;z9a<a9->GetSize();z9a++)
																				{
																					if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																					if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																					if ((!m_iRefOrSec[0][2]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));

																					if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																					if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																					vec->Add(0);

																					if ((!m_iRefOrSec[2][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));

																					if ((!m_iRefOrSec[2][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));

																					if ((!m_iRefOrSec[2][2]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g9->m_baAtomType[z9t]])->GetAt(a9->GetAt(z9a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g9->m_baAtomType[z9t]])->GetAt(a9->GetAt(z9a)));
																				}
																			}
																		}
																	}
																}
															}
														} else // IF NOT ORTHO[2]
														{
															g7 = (CAtomGroup*)m_oaVectors[z*9+6];
															g8 = (CAtomGroup*)m_oaVectors[z*9+7];
															for (z7t=0;z7t<g7->m_baAtomType.GetSize();z7t++)
															{
																a7 = (CxIntArray*)g7->m_oaAtoms[z7t];
																for (z7a=0;z7a<a7->GetSize();z7a++)
																{
																	for (z8t=0;z8t<g8->m_baAtomType.GetSize();z8t++)
																	{
																		a8 = (CxIntArray*)g8->m_oaAtoms[z8t];
																		for (z8a=0;z8a<a8->GetSize();z8a++)
																		{
																			if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																			if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																			if ((!m_iRefOrSec[0][2]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));

																			if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																			if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																			vec->Add(0);

																			if ((!m_iRefOrSec[2][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));

																			if ((!m_iRefOrSec[2][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));

																			vec->Add(0);
																		}
																	}
																}
															}
														} // END IF NOT ORTHO[2]
													}
												}
											}
										}
									} // END IF NOT ORTHO[1]
								}
							}
						}
					}
				}
			}
		} else // IF NOT ORTHO[0]
		{
			g1 = (CAtomGroup*)m_oaVectors[z*9];
			g2 = (CAtomGroup*)m_oaVectors[z*9+1];
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
							if (m_bOrtho[1])
							{
								g4 = (CAtomGroup*)m_oaVectors[z*9+3];
								g5 = (CAtomGroup*)m_oaVectors[z*9+4];
								g6 = (CAtomGroup*)m_oaVectors[z*9+5];
								for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
								{
									a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
									for (z4a=0;z4a<a4->GetSize();z4a++)
									{
										for (z5t=0;z5t<g5->m_baAtomType.GetSize();z5t++)
										{
											a5 = (CxIntArray*)g5->m_oaAtoms[z5t];
											for (z5a=0;z5a<a5->GetSize();z5a++)
											{
												for (z6t=0;z6t<g6->m_baAtomType.GetSize();z6t++)
												{
													a6 = (CxIntArray*)g6->m_oaAtoms[z6t];
													for (z6a=0;z6a<a6->GetSize();z6a++)
													{
														if (m_bOrtho[2])
														{
															g7 = (CAtomGroup*)m_oaVectors[z*9+6];
															g8 = (CAtomGroup*)m_oaVectors[z*9+7];
															g9 = (CAtomGroup*)m_oaVectors[z*9+8];
															for (z7t=0;z7t<g7->m_baAtomType.GetSize();z7t++)
															{
																a7 = (CxIntArray*)g7->m_oaAtoms[z7t];
																for (z7a=0;z7a<a7->GetSize();z7a++)
																{
																	for (z8t=0;z8t<g8->m_baAtomType.GetSize();z8t++)
																	{
																		a8 = (CxIntArray*)g8->m_oaAtoms[z8t];
																		for (z8a=0;z8a<a8->GetSize();z8a++)
																		{
																			for (z9t=0;z9t<g9->m_baAtomType.GetSize();z9t++)
																			{
																				a9 = (CxIntArray*)g9->m_oaAtoms[z9t];
																				for (z9a=0;z9a<a9->GetSize();z9a++)
																				{
																					if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																					if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																					vec->Add(0);

																					if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																					if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																					if ((!m_iRefOrSec[1][2]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));

																					if ((!m_iRefOrSec[2][0]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));

																					if ((!m_iRefOrSec[2][1]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));

																					if ((!m_iRefOrSec[2][2]) || (obs == NULL))
																						vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g9->m_baAtomType[z9t]])->GetAt(a9->GetAt(z9a)));
																					else
																						vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g9->m_baAtomType[z9t]])->GetAt(a9->GetAt(z9a)));
																				}
																			}
																		}
																	}
																}
															}
														} else // IF NOT ORTHO[2]
														{
															g7 = (CAtomGroup*)m_oaVectors[z*9+6];
															g8 = (CAtomGroup*)m_oaVectors[z*9+7];
															for (z7t=0;z7t<g7->m_baAtomType.GetSize();z7t++)
															{
																a7 = (CxIntArray*)g7->m_oaAtoms[z7t];
																for (z7a=0;z7a<a7->GetSize();z7a++)
																{
																	for (z8t=0;z8t<g8->m_baAtomType.GetSize();z8t++)
																	{
																		a8 = (CxIntArray*)g8->m_oaAtoms[z8t];
																		for (z8a=0;z8a<a8->GetSize();z8a++)
																		{
																			if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																			if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																			vec->Add(0);

																			if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																			if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																			if ((!m_iRefOrSec[1][2]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g6->m_baAtomType[z6t]])->GetAt(a6->GetAt(z6a)));

																			if ((!m_iRefOrSec[2][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));

																			if ((!m_iRefOrSec[2][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));

																			vec->Add(0);
																		}
																	}
																}
															}
														} // END IF NOT ORTHO[2]
													}
												}
											}
										}
									}
								}
							} else // IF NOT ORTHO[1]
							{
								g4 = (CAtomGroup*)m_oaVectors[z*9+3];
								g5 = (CAtomGroup*)m_oaVectors[z*9+4];
								for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
								{
									a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
									for (z4a=0;z4a<a4->GetSize();z4a++)
									{
										for (z5t=0;z5t<g5->m_baAtomType.GetSize();z5t++)
										{
											a5 = (CxIntArray*)g5->m_oaAtoms[z5t];
											for (z5a=0;z5a<a5->GetSize();z5a++)
											{
												if (m_bOrtho[2])
												{
													g7 = (CAtomGroup*)m_oaVectors[z*9+6];
													g8 = (CAtomGroup*)m_oaVectors[z*9+7];
													g9 = (CAtomGroup*)m_oaVectors[z*9+8];
													for (z7t=0;z7t<g7->m_baAtomType.GetSize();z7t++)
													{
														a7 = (CxIntArray*)g7->m_oaAtoms[z7t];
														for (z7a=0;z7a<a7->GetSize();z7a++)
														{
															for (z8t=0;z8t<g8->m_baAtomType.GetSize();z8t++)
															{
																a8 = (CxIntArray*)g8->m_oaAtoms[z8t];
																for (z8a=0;z8a<a8->GetSize();z8a++)
																{
																	for (z9t=0;z9t<g9->m_baAtomType.GetSize();z9t++)
																	{
																		a9 = (CxIntArray*)g9->m_oaAtoms[z9t];
																		for (z9a=0;z9a<a9->GetSize();z9a++)
																		{
																			if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																			if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																			vec->Add(0);

																			if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																			if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																			vec->Add(0);

																			if ((!m_iRefOrSec[2][0]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));

																			if ((!m_iRefOrSec[2][1]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));

																			if ((!m_iRefOrSec[2][2]) || (obs == NULL))
																				vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g9->m_baAtomType[z9t]])->GetAt(a9->GetAt(z9a)));
																			else
																				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g9->m_baAtomType[z9t]])->GetAt(a9->GetAt(z9a)));
																		}
																	}
																}
															}
														}
													}
												} else // IF NOT ORTHO[2]
												{
													g7 = (CAtomGroup*)m_oaVectors[z*9+6];
													g8 = (CAtomGroup*)m_oaVectors[z*9+7];
													for (z7t=0;z7t<g7->m_baAtomType.GetSize();z7t++)
													{
														a7 = (CxIntArray*)g7->m_oaAtoms[z7t];
														for (z7a=0;z7a<a7->GetSize();z7a++)
														{
															for (z8t=0;z8t<g8->m_baAtomType.GetSize();z8t++)
															{
																a8 = (CxIntArray*)g8->m_oaAtoms[z8t];
																for (z8a=0;z8a<a8->GetSize();z8a++)
																{
																	if ((!m_iRefOrSec[0][0]) || (obs == NULL))
																		vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
																	else
																		vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));

																	if ((!m_iRefOrSec[0][1]) || (obs == NULL))
																		vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
																	else
																		vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));

																	vec->Add(0);

																	if ((!m_iRefOrSec[1][0]) || (obs == NULL))
																		vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
																	else
																		vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));

																	if ((!m_iRefOrSec[1][1]) || (obs == NULL))
																		vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));
																	else
																		vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g5->m_baAtomType[z5t]])->GetAt(a5->GetAt(z5a)));

																	vec->Add(0);

																	if ((!m_iRefOrSec[2][0]) || (obs == NULL))
																		vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));
																	else
																		vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g7->m_baAtomType[z7t]])->GetAt(a7->GetAt(z7a)));

																	if ((!m_iRefOrSec[2][1]) || (obs == NULL))
																		vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));
																	else
																		vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g8->m_baAtomType[z8t]])->GetAt(a8->GetAt(z8a)));

																	vec->Add(0);
																}
															}
														}
													}
												} // END IF NOT ORTHO[2]
											}
										}
									}
								}
							} // END IF NOT ORTHO[1]
						}
					}
				}
			}
		} // END IF NOT ORTHO[0]
	} // END FOR ALL SETS
	BXOUT;
}


void CPlDF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z, z1t, z1a, z2t, z2a, z3t, z3a, z4t, z4a;
	CAtomGroup *g1, *g2, *g3, *g4;
	CxIntArray *a1, *a2, *a3, *a4;

	vec->RemoveAll_KeepSize();
	for (z=0;z<m_oaVectors.GetSize()/4;z++)
	{
		if (m_bNormal)
		{
			g1 = (CAtomGroup*)m_oaVectors[z*4];
			g2 = (CAtomGroup*)m_oaVectors[z*4+1];
			g4 = (CAtomGroup*)m_oaVectors[z*4+3];

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
							for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
							{
								a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
								for (z4a=0;z4a<a4->GetSize();z4a++)
								{
									if ((!m_iRefOrSec[0]) || (obs == NULL))
										vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
											else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
									if ((!m_iRefOrSec[1]) || (obs == NULL))
										vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
											else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
									vec->Add(0);
									if ((!m_iRefOrSec[3]) || (obs == NULL))
										vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
											else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
								}
							}
						}
					}
				}
			}
		} else // IF NOT NORMAL
		{
			g1 = (CAtomGroup*)m_oaVectors[z*4];
			g2 = (CAtomGroup*)m_oaVectors[z*4+1];
			g3 = (CAtomGroup*)m_oaVectors[z*4+2];
			g4 = (CAtomGroup*)m_oaVectors[z*4+3];

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
									for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
									{
										a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
										for (z4a=0;z4a<a4->GetSize();z4a++)
										{
											if ((!m_iRefOrSec[0]) || (obs == NULL))
												vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
													else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
											if ((!m_iRefOrSec[1]) || (obs == NULL))
												vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
													else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
											if ((!m_iRefOrSec[2]) || (obs == NULL))
												vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
													else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
											if ((!m_iRefOrSec[3]) || (obs == NULL))
												vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
													else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
										}
									}
								}
							}
						}
					}
				}
			}
		}
	} // END FOR ALL SETS
	BXOUT;
}


void CLiDF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z, z1t, z1a, z2t, z2a, z3t, z3a, z4t, z4a;
	CAtomGroup *g1, *g2, *g3, *g4;
	CxIntArray *a1, *a2, *a3, *a4;

	vec->RemoveAll_KeepSize();
	for (z=0;z<m_oaVectors.GetSize()/4;z++)
	{
		if (m_bNormal)
		{
			g1 = (CAtomGroup*)m_oaVectors[z*4];
			g2 = (CAtomGroup*)m_oaVectors[z*4+1];
			g3 = (CAtomGroup*)m_oaVectors[z*4+2];
			g4 = (CAtomGroup*)m_oaVectors[z*4+3];

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
									for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
									{
										a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
										for (z4a=0;z4a<a4->GetSize();z4a++)
										{
											if ((!m_iRefOrSec[0]) || (obs == NULL))
												vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
													else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
											if ((!m_iRefOrSec[1]) || (obs == NULL))
												vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
													else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
											if ((!m_iRefOrSec[2]) || (obs == NULL))
												vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
													else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g3->m_baAtomType[z3t]])->GetAt(a3->GetAt(z3a)));
											if ((!m_iRefOrSec[3]) || (obs == NULL))
												vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
													else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
										}
									}
								}
							}
						}
					}
				}
			}
		} else // If not Normal
		{
			g1 = (CAtomGroup*)m_oaVectors[z*4];
			g2 = (CAtomGroup*)m_oaVectors[z*4+1];
			g4 = (CAtomGroup*)m_oaVectors[z*4+3];

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
							for (z4t=0;z4t<g4->m_baAtomType.GetSize();z4t++)
							{
								a4 = (CxIntArray*)g4->m_oaAtoms[z4t];
								for (z4a=0;z4a<a4->GetSize();z4a++)
								{
									if ((!m_iRefOrSec[0]) || (obs == NULL))
										vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
											else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
									if ((!m_iRefOrSec[1]) || (obs == NULL))
										vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
											else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
									vec->Add(0);
									if ((!m_iRefOrSec[3]) || (obs == NULL))
										vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
											else vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g4->m_baAtomType[z4t]])->GetAt(a4->GetAt(z4a)));
								}
							}
						}
					}
				}
			}
		}
	} // END FOR ALL SETS
	BXOUT;
}


void CRDF::BuildName()
{
	BTIN;
	int z;
//	char tmp[32768];
	CxString tmp;

//	tmp[0] = 0;
	tmp.sprintf("");

	if (m_iDeriv != 0)
//		sprintf(tmp,"deriv%d_",m_iDeriv);
		tmp.sprintf("deriv%d_",m_iDeriv);

	for (z=0;z<m_oaVectors.GetSize()/2;z++)
	{
//		strcat(tmp,"[");
//		strcat(tmp,((CAtomGroup*)m_oaVectors[z*2])->m_sName);
		tmp.strcat("[");
		tmp.strcat(((CAtomGroup*)m_oaVectors[z*2])->m_sName);

		if (m_iRefOrSec[0])
//			strcat(tmp,"o_");
			tmp.strcat("o_");
		else
//			strcat(tmp,"r_");
			tmp.strcat("r_");

//		strcat(tmp,((CAtomGroup*)m_oaVectors[z*2+1])->m_sName);
		tmp.strcat(((CAtomGroup*)m_oaVectors[z*2+1])->m_sName);

		if (m_iRefOrSec[1])
//			strcat(tmp,"o");
			tmp.strcat("o");
		else
//			strcat(tmp,"r");
			tmp.strcat("r");

		if (z < (m_oaVectors.GetSize()/2)-1)
//			strcat(tmp,"],");
			tmp.strcat("],");
		else
//			strcat(tmp,"]");
			tmp.strcat("]");
	}

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sShortName,(const char*)tmp);

	if (m_iShowMol != -1)
//		sprintf(tmp,"%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.sprintf("%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	else
//		sprintf(tmp,"%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);
		tmp.sprintf("%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);

//	strcat(tmp,m_sShortName);
	tmp.strcat(m_sShortName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);

	tmp.sprintf("");
	if (m_iDeriv != 0)
		tmp.sprintf("deriv%d_",m_iDeriv);
	tmp.strcat("[");
	tmp.strcat(((CAtomGroup*)m_oaVectors[0])->m_sName);
	if (m_iRefOrSec[0])
		tmp.strcat("o_");
	else
		tmp.strcat("r_");
	tmp.strcat(((CAtomGroup*)m_oaVectors[1])->m_sName);
	if (m_iRefOrSec[1])
		tmp.strcat("o");
	else
		tmp.strcat("r");
	tmp.strcat("]");
	if (m_oaVectors.GetSize() > 2)
		tmp.strcat(",(...)");
	try { m_sLabelName = new char[strlen(tmp)+1]; } catch(...) { m_sLabelName = NULL; }
	if (m_sLabelName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sLabelName,tmp);

	BTOUT;
}


void CDensDF::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;
	int z;

//	sprintf(tmp,"%s_%s%d_%s",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,((CAtom*)g_oaAtoms[m_iCenterAtomRealType])->m_sName,m_iCenterAtom+1,(m_bDensityMass?"mass":"particle"));
	tmp.sprintf("%s_%s%d_%s",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,((CAtom*)g_oaAtoms[m_iCenterAtomRealType])->m_sName,m_iCenterAtom+1,(m_bDensityMass?"mass":"particle"));

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		if (!m_pDensityMolSelect[z])
			continue;

//		strcat(tmp,"_");
//		strcat(tmp,((CMolecule*)g_oaMolecules[z])->m_sName);
		tmp.strcat("_");
		tmp.strcat(((CMolecule*)g_oaMolecules[z])->m_sName);

		if (!m_pDensityMolAG[z]->m_bAllAtoms)
		{
//			strcat(tmp,"_");
//			strcat(tmp,m_pDensityMolAG[z]->m_sName);
			tmp.strcat("_");
			tmp.strcat(m_pDensityMolAG[z]->m_sName);
		}
	}

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CVHDF::BuildName()
{
	BTIN;
	int z;
//	char tmp[32768];
	CxString tmp;

//	tmp[0] = 0;
	tmp.sprintf("");

	for (z=0;z<m_oaVectors.GetSize()/2;z++)
	{
//		strcat(tmp,"[");
//		strcat(tmp,((CAtomGroup*)m_oaVectors[z*2])->m_sName);
		tmp.strcat("[");
		tmp.strcat(((CAtomGroup*)m_oaVectors[z*2])->m_sName);

		if (m_iRefOrSec[0])
//			strcat(tmp,"o_");
			tmp.strcat("o_");
		else
//			strcat(tmp,"r_");
			tmp.strcat("r_");

//		strcat(tmp,((CAtomGroup*)m_oaVectors[z*2+1])->m_sName);
		tmp.strcat(((CAtomGroup*)m_oaVectors[z*2+1])->m_sName);

		if (m_iRefOrSec[1])
//			strcat(tmp,"o");
			tmp.strcat("o");
		else
//			strcat(tmp,"r");
			tmp.strcat("r");

		if (z < (m_oaVectors.GetSize()/2)-1)
//			strcat(tmp,"],");
			tmp.strcat("],");
		else
//			strcat(tmp,"]");
			tmp.strcat("]");
	}

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sShortName,tmp);

	if (m_iShowMol != -1)
//		sprintf(tmp,"%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.sprintf("%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	else
//		sprintf(tmp,"%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);
		tmp.sprintf("%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);

//	strcat(tmp,m_sShortName);
	tmp.strcat(m_sShortName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CRDF::Parse()
{
	BTIN;
//	char buf[1024];
	CxString buf;
	int ti;
	CAtomGroup *ag;

	try { m_pRDF = new CDF(); } catch(...) { m_pRDF = NULL; }
	if (m_pRDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_iShowAtomGes = 0;
	m_iRefAtomGes = 0;
	m_iCombinations = 0;
	mprintf(WHITE,"\n>>> Radial Distribution Function >>>\n\n");
/*	if (m_bSelf)
	{
		m_iRefOrSec[0] = 0;
		m_iRefOrSec[1] = 0;
	} else
	{
		m_iRefOrSec[0] = 0;
		m_iRefOrSec[1] = 1;
	}*/
	if (m_iShowMol != -1)
		m_iRefOrSec[0] = AskRangeInteger("    Take reference atom(s) from RM %s (0) or from OM %s (1)? [0] ",0,1,0,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			else m_iRefOrSec[0] = 0; // Kein OM: Nimm alles aus RM
	if (m_iShowMol != -1)
		m_iRefOrSec[1] = AskRangeInteger("    Take observed atom(s) from RM %s (0) or from OM %s (1)? [1] ",0,1,1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			else m_iRefOrSec[1] = 0; // Kein OM: Nimm alles aus RM
		
_rdfnewset:
	mprintf("\n");

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
_rdfatom1:
	// 1 reales + 2 virtuelle = 3 gesamt
	if (((CMolecule*)g_oaMolecules[(m_iRefOrSec[0])?m_iShowMol:g_iFixMol])->m_iAtomGes == 3)
	{
		mprintf("    %s is only one atom, there is no choice.\n",((CMolecule*)g_oaMolecules[(m_iRefOrSec[0])?m_iShowMol:g_iFixMol])->m_sName);
		ag->Reset();
		ag->m_pMolecule = (CMolecule*)g_oaMolecules[(m_iRefOrSec[0])?m_iShowMol:g_iFixMol];
		ag->AddAtom(0,0,false);
		ag->SortAtoms();
		ag->BuildName();
	} else
	{
		if (m_iRefOrSec[0])
		{
			mprintf("    Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			inpprintf("! Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		} else
		{
			mprintf("    Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
			inpprintf("! Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		}
		myget(&buf);
		if (strlen(buf) == 0)
		{
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec[0])?m_iShowMol:g_iFixMol],"#2"))
			{
				eprintf("Weird error.\n");
				goto _rdfatom1;
			}
		} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec[0])?m_iShowMol:g_iFixMol],buf))
			goto _rdfatom1;
	}
	m_oaVectors.Add(ag);
	m_iRefAtomGes += ag->m_iAtomGes;
	ti = ag->m_iAtomGes;

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
_rdfatom2:
	if (((CMolecule*)g_oaMolecules[(m_iRefOrSec[1])?m_iShowMol:g_iFixMol])->m_iAtomGes == 3)
	{
		mprintf("    %s is only one atom, there is no choice.\n",((CMolecule*)g_oaMolecules[(m_iRefOrSec[1])?m_iShowMol:g_iFixMol])->m_sName);
		ag->Reset();
		ag->m_pMolecule = (CMolecule*)g_oaMolecules[(m_iRefOrSec[1])?m_iShowMol:g_iFixMol];
		ag->AddAtom(0,0,false);
		ag->SortAtoms();
		ag->BuildName();
	} else
	{
		if (m_iRefOrSec[1])
		{
			mprintf("    Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			inpprintf("! Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		} else
		{
			mprintf("    Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
			inpprintf("! Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		}
		myget(&buf);

		try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
		if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if (strlen(buf) == 0)
		{
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec[1])?m_iShowMol:g_iFixMol],"#2"))
			{
				eprintf("Weird error.\n");
				goto _rdfatom2;
			}
		} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec[1])?m_iShowMol:g_iFixMol],buf))
			goto _rdfatom2;
	}
	m_oaVectors.Add(ag);
	m_iShowAtomGes += ag->m_iAtomGes;
	m_iCombinations += ti * ag->m_iAtomGes;

	if (g_bAdvanced2)
		if (AskYesNo("    Add another set of atoms to this (!) RDF (y/n)? [no] ",false))
			goto _rdfnewset;

	mprintf("    This yields in %d combinations.\n\n",m_iCombinations);

	ParseDeriv();

	switch(m_iDeriv)
	{
		case 0:
			m_fMinDist = AskFloat("    Enter the minimal radius of this RDF in pm: [0] ",0.0f);
			m_fMaxDist = AskFloat("    Enter the maximal radius of this RDF in pm: [%d.0] ",(float)HalfBox(),HalfBox());
			break;
		case 1:
			if (m_bDerivAbs)
				m_fMinDist = AskFloat("    Enter the minimal value of this d1-RDF in pm/ps: [0] ",0.0f);
					else m_fMinDist = AskFloat("    Enter the minimal value of this d1-RDF in pm/ps: [-10.0] ",-10.0f);
			m_fMaxDist = AskFloat("    Enter the maximal value of this d1-RDF in pm/ps: [10.0] ",10.0f);
			break;
		case 2:
			if (m_bDerivAbs)
				m_fMinDist = AskFloat("    Enter the minimal value of this d2-RDF in pm/ps^2: [0] ",0.0f);
					else m_fMinDist = AskFloat("    Enter the minimal value of this d2-RDF in pm/ps^2: [-10.0] ",-10.0f);
			m_fMaxDist = AskFloat("    Enter the maximal value of this d2-RDF in pm/ps^2: [10.0] ",10.0f);
			break;
	}

	if (g_bPeriodic && (m_iDeriv == 0) && (m_fMaxDist > HalfBox_Exact()+1.0f))
	{
		eprintf("\nWarning: ");
		mprintf("The specified max. radius is larger than half of the smallest periodic cell vector.\n");
		mprintf("         TRAVIS counts every atom only once (central periodic image).\n");
		mprintf("         Expect the analysis to decay to zero for large radii.\n\n");
		AskYesNo("         Acknowledged [yes] ",true);
		mprintf("\n");
	}

	m_bAdaptive = false/*AskYesNo("    Enter binning resolution (n) or use adaptive binnig (y)? [no] ",false)*/;
	if (!m_bAdaptive)
		m_iResolution = AskUnsignedInteger("    Enter the resolution (bin count) for this RDF: [300] ",/*(int)((m_fMaxDist-m_fMinDist)/10.0f),(int)((m_fMaxDist-m_fMinDist)/10.0f)*/300);
			else m_iResolution = 65536;

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
			else m_iHistogramRes = 0;

	if (m_iDeriv == 0)
	{
		if (m_iShowMol != -1)
			m_bRadialCorrect = AskYesNo("    Correct radial distribution for this RDF (y/n)? [yes] ",true);
				else m_bRadialCorrect = AskYesNo("    Correct radial distribution for this RDF (y/n)? [no] ",false);
	} else m_bRadialCorrect = false;

	if (g_bAdvanced2 && m_bRadialCorrect)
		m_bProbDens = AskYesNo("    Compute occurence in nm^(-3) (y) or rel. to uniform density (n)? [no] ",false);
			else m_bProbDens = false;

		m_bLongMode = false;

	if (m_bCalcSD)
	{
		m_iSDBlocks = AskUnsignedInteger("    How many different block lenghts to use? [100] ",100);
		m_iSDBlockMin = AskUnsignedInteger("    Enter minimal block length in time steps: [1] ",1);
		if (g_iTrajSteps != -1)
			m_iSDBlockMax = AskUnsignedInteger("    Enter maximal block length in time steps: [%d] ",g_iTrajSteps/10,g_iTrajSteps/10);
				else m_iSDBlockMax = AskUnsignedInteger("    Enter maximal block length in time steps: [1000] ",1000);
		m_bSDVerbose = AskYesNo("    Write out correlation length extrapolation fit data for each bin (y/n)? [no] ",false);
		m_fSDTimesSigma = AskFloat("    Use which factor of sigma for confidence range? [3.0] ",3.0);
		mprintf("\n    %.2f sigma leads to a confidence level of %.5f%c.\n",m_fSDTimesSigma,2.0*NormalDistIntegral(m_fSDTimesSigma)-1.0,'%');
	}

	BuildName();
	mprintf(WHITE,"\n<<< End of Radial Distribution Function <<<\n\n");
	BTOUT;
}


void CDensDF::Parse()
{
	BTIN;
//	char buf[1024];
	CxString buf;
	int z;
	CMolecule *m;

	try { m_pDensDF = new CDF(); } catch(...) { m_pDensDF = NULL; }
	if (m_pDensDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Density Distribution Function >>>\n\n");

	m_bDensityMass = AskYesNo("    Observe mass density (y) or particle density (n)? [yes] ",true);

	mprintf("\n    Choose a reference atom around which the density will be analyzed.\n\n");

	// 1 reales + 2 virtuelle = 3 gesamt
	if (((CMolecule*)g_oaMolecules[m_iShowMol])->m_iAtomGesNoVirt == 1)
	{
		mprintf("    %s is only one atom, there is no choice.\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		if (!ParseAtom("#2",m_iShowMol,m_iCenterAtomType,m_iCenterAtomRealType,m_iCenterAtom))
		{
			eprintf("Weird error.\n");
			abort();
		}
	} else
	{
_densatom1:
		mprintf("    Which atom to take from OM %s (e.g. C1)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		inpprintf("! Which atom to take from OM %s (e.g. C1)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		myget(&buf);
		if (strlen(buf) == 0)
		{
			if (!ParseAtom("#2",m_iShowMol,m_iCenterAtomType,m_iCenterAtomRealType,m_iCenterAtom))
			{
				eprintf("Weird error.\n");
				abort();
			}
		} else if (!ParseAtom(buf,m_iShowMol,m_iCenterAtomType,m_iCenterAtomRealType,m_iCenterAtom))
			goto _densatom1;
	}

	try { m_pDensityMolSelect = new bool[g_oaMolecules.GetSize()]; } catch(...) { m_pDensityMolSelect = NULL; }
	if (m_pDensityMolSelect == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_pDensityMolAG = new CAtomGroup*[g_oaMolecules.GetSize()]; } catch(...) { m_pDensityMolAG = NULL; }
	if (m_pDensityMolAG == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(CAtomGroup*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf("\n    Please choose the atoms to observe:\n\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		m_pDensityMolSelect[z] = AskYesNo("    Consider contributions from %s (y/n)? [%s] ",!m->m_bPseudo,((CMolecule*)g_oaMolecules[z])->m_sName,m->m_bPseudo?"no":"yes");
		if (m_pDensityMolSelect[z])
		{
			try { m_pDensityMolAG[z] = new CAtomGroup(); } catch(...) { m_pDensityMolAG[z] = NULL; }
			if (m_pDensityMolAG[z] == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
_densatom2:
			AskString("      Which atoms from %s to observe (e.g. C1-3,C6,N)? [all] ",&buf,"",m->m_sName);
			if (strlen(buf) == 0)
			{
				m_pDensityMolAG[z]->AddAllAtoms(m,false);
			} else if (!m_pDensityMolAG[z]->ParseAtoms(m,buf))
				goto _densatom2;
		}
	}

	mprintf("\n");

	m_fMinDist = AskFloat("    Enter the minimal radius of this Density DF in pm: [0] ",0.0f);
	m_fMaxDist = AskFloat("    Enter the maximal radius of this Density DF in pm: [%d.0] ",(float)HalfBox(),HalfBox());

	m_iResolution = AskUnsignedInteger("    Enter the resolution (bin count) for this Density DF: [300] ",300);

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
			else m_iHistogramRes = 0;

	BuildName();
	mprintf(WHITE,"\n<<< End of Density Distribution Function <<<\n\n");
	BTOUT;
}


void CRDF::ParseCondition(int rm, CNbSearch *n, bool nbana)
{
	BTIN;
//	char buf[1024];
	CxString buf;
	int ti;
	CAtomGroup *ag;

	try { m_pRDF = new CDF(); } catch(...) { m_pRDF = NULL; }
	if (m_pRDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_iShowAtomGes = 0;
	m_iRefAtomGes = 0;
	m_iCombinations = 0;
	m_iRefOrSec[0] = 0;
	m_iRefOrSec[1] = 1;
	mprintf(WHITE,"\n>>> Distance Condition >>>\n");
	
_rdfnewset:
	mprintf("\n");

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

_rdfatom1:
	if (((CMolecule*)g_oaMolecules[rm])->m_iAtomGes == 3)
	{
		mprintf("    %s is only one atom, there is no choice.\n",((CMolecule*)g_oaMolecules[rm])->m_sName);
		ag->Reset();
		ag->m_pMolecule = (CMolecule*)g_oaMolecules[rm];
		ag->AddAtom(0,0,false);
		ag->SortAtoms();
		ag->BuildName();
	} else
	{
		mprintf("    Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[rm])->m_sName);
		inpprintf("! Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[rm])->m_sName);
		myget(&buf);
		if (strlen(buf) == 0)
		{
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[rm],"#2"))
			{
				eprintf("Weird error.\n");
				goto _rdfatom1;
			}
		} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[rm],buf))
			goto _rdfatom1;
	}
	m_oaVectors.Add(ag);
	m_iRefAtomGes += ag->m_iAtomGes;
	ti = ag->m_iAtomGes;

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

_rdfatom2:
	if (((CMolecule*)g_oaMolecules[m_iShowMol])->m_iAtomGes == 3)
	{
		mprintf("    %s is only one atom, there is no choice.\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		ag->Reset();
		ag->m_pMolecule = (CMolecule*)g_oaMolecules[m_iShowMol];
		ag->AddAtom(0,0,false);
		ag->SortAtoms();
		ag->BuildName();
	} else
	{
		mprintf("    Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		inpprintf("! Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		myget(&buf);

		try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
		if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if (strlen(buf) == 0)
		{
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],"#2"))
			{
				eprintf("Weird error.\n");
				goto _rdfatom2;
			}
		} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
			goto _rdfatom2;
	}
	m_oaVectors.Add(ag);
	m_iShowAtomGes += ag->m_iAtomGes;
	m_iCombinations += ti * ag->m_iAtomGes;

	if (AskYesNo("    Enter another set of atoms for this condition (y/n)? [no] ",false))
		goto _rdfnewset;

	if (!nbana)
	{
		if (AskUnsignedInteger("\n    Enter min./max. distance (0) or min./max. nearest neighbor count (1)? [0] ",0)==0)
		{
			g_bEnvDisableSortNb = true;

			do {
				m_faMinMaxDist.Add(AskFloat("    Enter the minimal distance in pm: [0] ",0.0f));
				m_faMinMaxDist.Add(AskFloat("    Enter the maximal distance in pm: [400] ",400.0f));
			} while (AskYesNo("    Enter another distance interval (y/n)? [no] ",false));
			n->m_iNbCountMin = -1;
			n->m_iNbCountMax = -1;
		} else
		{
			n->m_iNbCountMin = AskRangeInteger("    Use next neighbors from the n-th on? [1] ",0,((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(),1)-1;
			n->m_iNbCountMax = AskRangeInteger("    Use next neighbors up to the n-th? [%d] ",n->m_iNbCountMin+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(),n->m_iNbCountMin+1,n->m_iNbCountMin+1)-1;
//			m_faMinMaxDist.Add(0);
//			m_faMinMaxDist.Add(9E20f);
		}
	} else
	{
		m_faMinMaxDist.Add(0);
		m_faMinMaxDist.Add(1.0e30f);
		n->m_iNbCountMin = -2;
		n->m_iNbCountMax = -2;
	}
	mprintf(WHITE,"\n<<< End of Distance Condition <<<\n\n");
	BTOUT;
}


void CRDF::ParseConditionGrid(int rm, CNbSearch *n, int gridmode)
{
	BTIN;
//	char buf[1024];
	CxString buf;
	int ti;
	CAtomGroup *ag;

	try { m_pRDF = new CDF(); } catch(...) { m_pRDF = NULL; }
	if (m_pRDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_iShowAtomGes = 0;
	m_iRefAtomGes = 0;
	m_iCombinations = 0;
	m_iRefOrSec[0] = 0;
	m_iRefOrSec[1] = 1;
	mprintf(WHITE,"\n>>> Distance Condition >>>\n");

_rdfnewset:
	mprintf("\n");

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

_rdfatom1:
	if (((CMolecule*)g_oaMolecules[rm])->m_iAtomGes == 3)
	{
		mprintf("    %s is only one atom, there is no choice.\n",((CMolecule*)g_oaMolecules[rm])->m_sName);
		ag->Reset();
		ag->m_pMolecule = (CMolecule*)g_oaMolecules[rm];
		ag->AddAtom(0,0,false);
		ag->SortAtoms();
		ag->BuildName();
	} else
	{
		mprintf("    Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[rm])->m_sName);
		inpprintf("! Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[rm])->m_sName);
		myget(&buf);
		if (strlen(buf) == 0)
		{
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[rm],"#2"))
			{
				eprintf("Weird error.\n");
				goto _rdfatom1;
			}
		} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[rm],buf))
			goto _rdfatom1;
	}
	m_oaVectors.Add(ag);
	m_iRefAtomGes += ag->m_iAtomGes;
	ti = ag->m_iAtomGes;

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

_rdfatom2:
	if (((CMolecule*)g_oaMolecules[m_iShowMol])->m_iAtomGes == 3)
	{
		mprintf("    %s is only one atom, there is no choice.\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		ag->Reset();
		ag->m_pMolecule = (CMolecule*)g_oaMolecules[m_iShowMol];
		ag->AddAtom(0,0,false);
		ag->SortAtoms();
		ag->BuildName();
	} else
	{
		mprintf("    Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		inpprintf("! Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		myget(&buf);

		try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
		if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if (strlen(buf) == 0)
		{
			if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],"#2"))
			{
				eprintf("Weird error.\n");
				goto _rdfatom2;
			}
		} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
			goto _rdfatom2;
	}
	m_oaVectors.Add(ag);
	m_iShowAtomGes += ag->m_iAtomGes;
	m_iCombinations += ti * ag->m_iAtomGes;

	if (AskYesNo("    Enter another set of atoms (y/n)? [no] ",false))
		goto _rdfnewset;

	if (gridmode == 6)
	{
		if (AskUnsignedInteger("\n    Enter min./max. distance (0) or min./max. nearest neighbor count (1)? [0] ",0)==0)
		{
			m_faMinMaxDist.Add(AskFloat("    Enter the minimal distance in pm: [0] ",0.0f));
			m_faMinMaxDist.Add(AskFloat("    Enter the maximal distance in pm: [400] ",400.0f));
			n->m_iNbCountMin = -1;
			n->m_iNbCountMax = -1;
		} else
		{
			n->m_iNbCountMin = AskRangeInteger("    Use next neighbors from the n-th on? [1] ",0,((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(),1)-1;
			n->m_iNbCountMax = AskRangeInteger("    Use next neighbors up to the n-th? [%d] ",n->m_iNbCountMin+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(),n->m_iNbCountMin+1,n->m_iNbCountMin+1)-1;
//			m_faMinMaxDist.Add(0);
//			m_faMinMaxDist.Add(9E20f);
		}
	} else if ((gridmode == 4) || (gridmode == 5))
	{
		n->m_iNbCountMin = 0;
		n->m_iNbCountMax = 0;
//		m_faMinMaxDist.Add(0);
//		m_faMinMaxDist.Add(9E20f);
	} else if (gridmode == 2)
	{
		m_faMinMaxDist.Add(AskFloat("    Enter the minimal distance in pm: [0] ",0.0f));
		m_faMinMaxDist.Add(AskFloat("    Enter the maximal distance in pm: [400] ",400.0f));
		n->m_iNbCountMin = -1;
		n->m_iNbCountMax = -1;
	} else 
	{
		m_faMinMaxDist.Add(0.0f);
		m_faMinMaxDist.Add(400.0f);
		n->m_iNbCountMin = -1;
		n->m_iNbCountMax = -1;
	}
	mprintf(WHITE,"\n<<< End of Distance Condition <<<\n\n");
	BTOUT;
}


void CRDF::ParseCondition_OnlyValues(CNbSearch *n)
{
	BTIN;
	int z;

	mprintf("    Distance condition between ");
	for (z=0;z<m_oaVectors.GetSize()/2;z++)
	{
		mprintf("%s and %s",((CAtomGroup*)m_oaVectors[z*2])->m_sName,((CAtomGroup*)m_oaVectors[z*2+1])->m_sName);
		if (z < (m_oaVectors.GetSize()/2)-1)
			mprintf(", ");
	}
	mprintf("\n");
	if (n->m_iNbCountMin == -1)
	{
		for (z=0;z<m_faMinMaxDist.GetSize()/2;z++)
		{
			mprintf("    Interval %d:\n",z+1,m_faMinMaxDist.GetSize()/2);
			m_faMinMaxDist[z*2] = AskFloat("      Enter the minimal distance in pm: [0] ",0.0f);
			m_faMinMaxDist[z*2+1] = AskFloat("      Enter the maximal distance in pm: [400] ",400.0f);
		}
	} else
	{
		n->m_iNbCountMin = AskRangeInteger("    Use next neighbors from the n-th on? [1] ",0,((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(),1)-1;
		n->m_iNbCountMax = AskRangeInteger("    Use next neighbors up to the n-th? [%d] ",n->m_iNbCountMin+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(),n->m_iNbCountMin+1,n->m_iNbCountMin+1)-1;
	}

	BTOUT;
}


void CVHDF::Parse()
{
	BTIN;
//	char buf[1024];
	CxString buf;
	int ti;
	CAtomGroup *ag;

	try { m_pVHDF = new C2DF(); } catch(...) { m_pVHDF = NULL; }
	if (m_pVHDF == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_iShowAtomGes = 0;
	m_iRefAtomGes = 0;
	m_iCombinations = 0;
	mprintf(WHITE,"\n>>> Van Hove Correlation Function >>>\n\n");
	if (m_iShowMol != -1)
	{
		m_iRefOrSec[0] = AskRangeInteger("    Take (fixed) reference atom(s) from RM %s (0) or from OM %s (1)? [0] ",0,1,0,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		m_iRefOrSec[1] = AskRangeInteger("    Take (moving) observed atom(s) from RM %s (0) or from OM %s (1)? [1] ",0,1,1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	} else // Kein OM: Nimm alles aus RM
	{
		m_iRefOrSec[0] = 0;
		m_iRefOrSec[1] = 0;
	}
_rdfnewset:
	mprintf("\n");

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

_rdfatom1:
	if (m_iRefOrSec[0])
	{
		mprintf("    Which (fixed) reference atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		inpprintf("! Which (fixed) reference atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	} else
	{
		mprintf("    Which (fixed) reference atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		inpprintf("! Which (fixed) reference atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
	}
	myget(&buf);
	if (strlen(buf) == 0)
	{
		if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec[0])?m_iShowMol:g_iFixMol],"#2"))
			goto _rdfatom1;
	} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec[0])?m_iShowMol:g_iFixMol],buf))
		goto _rdfatom1;
	m_oaVectors.Add(ag);
	m_iRefAtomGes += ag->m_iAtomGes;
	ti = ag->m_iAtomGes;

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

_rdfatom2:
	if (m_iRefOrSec[1])
	{
		mprintf("    Which (moving) observed atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		inpprintf("! Which (moving) observed atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	} else
	{
		mprintf("    Which (moving) observed atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		inpprintf("! Which (moving) observed atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
	}
	myget(&buf);

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (strlen(buf) == 0)
	{
		if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec[1])?m_iShowMol:g_iFixMol],"#2"))
			goto _rdfatom2;
	} else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec[1])?m_iShowMol:g_iFixMol],buf))
		goto _rdfatom2;
	m_oaVectors.Add(ag);
	m_iShowAtomGes += ag->m_iAtomGes;
	m_iCombinations += ti * ag->m_iAtomGes;

	if (g_bAdvanced2)
		if (AskYesNo("    Enter another set of atoms (y/n)? [no] ",false))
			goto _rdfnewset;

	m_fMinDist = 0; //AskFloat("    Enter the minimal radius of this VHCF in pm: [0] ",0.0f);
	m_fMaxDist = AskFloat("    Enter the radius of this VHCF in pm: [%d.0] ",(float)HalfBox(),HalfBox());
	m_iResolution = AskUnsignedInteger("    Enter the radial resolution of this VHCF: [100] ",100);
	m_bRadialCorrect = AskYesNo("    Correct radial distribution for this VHCF (y/n)? [%s] ",!m_bSelf,m_bSelf?"no":"yes");

_depth:
	if (g_iTrajSteps != -1)
		m_iDepth = AskUnsignedInteger("    Enter the temporal resolution (time depth) of this VHCF (in time steps): [%d] ",int(g_iTrajSteps*0.75),int(g_iTrajSteps*0.75));
			else m_iDepth = AskUnsignedInteger("    Enter the temporal resolution (time depth) of this VHCF (in time steps): [5000] ",5000);

	mprintf("\n    This will occupy %s of RAM.\n",FormatBytes((double)m_iDepth*g_iGesVirtAtomCount*3.0*sizeof(float)));
	if (m_iDepth*g_iGesVirtAtomCount*3.0*sizeof(float)/1024.0f/1024.0f >= 1000.0f)
		if (!AskYesNo("    Make sure that enough free RAM is available. Continue (y/n)? [yes] ",true))
			goto _depth;
	mprintf("\n");

_resagain:
	m_iStride = AskUnsignedInteger("    Take each n-th time step for the temporal axis? [%d] ",max(1,m_iDepth/100),max(1,m_iDepth/100));

	mprintf("\n    This results in a plot resolution of %d on the temporal axis.\n",m_iDepth/m_iStride);

	if (m_iDepth/m_iStride > 200)
	{
		mprintf("\n");
		if (!AskYesNo("    The resolution seems quite high, the plot will take much time to render. Contiune (y/n)? [yes] ",true))
			goto _resagain;
	}

	if (g_bAdvanced2)
	{
		mprintf("\n");
		m_bSwapAxes = AskYesNo("    Put distance on X axis and time on Y axis (y) or swap axes (n)? [yes] ",true);
		m_iGraceBunchTime = AskUnsignedInteger("    How many time graphs do you want do draw in the distance grace stack (0=disable)? [10] ",10);
		m_iGraceBunchDist = AskUnsignedInteger("    How many distance graphs do you want do draw in the time grace stack (0=disable)? [10] ",10);
	} else
	{
		m_bSwapAxes = true;
		m_iGraceBunchTime = 0;
		m_iGraceBunchDist = 0;
	}

	BuildName();
	mprintf(WHITE,"\n<<< End of Van Hove Correlation Function <<<\n\n");
	BTOUT;
}


void CADF::Parse()
{
	BTIN;
//	char buf[1024];
	CxString buf;
	int z, z2, ti;
	CAtomGroup *ag;
	float tf;
	bool dip;

	try { m_pADF = new CDF(); } catch(...) { m_pADF = NULL; }
	if (m_pADF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Angular Distribution Function >>>\n\n");
	dip = false;
	for (z=0;z<2;z++)
	{
		m_iVecType[z] = AskRangeInteger("    Should the %d. vector depict position (1), dipole (2), velocity (3) or force (4)? [1] ",1,4,1,z+1) - 1;
		if (m_iVecType[z] == 1)
		{
			g_bDipole = true;
			ParseDipole();
		}
		if (!dip)
		{
			dip = true;
			mprintf(YELLOW,"    Please note: ");
			mprintf("The last virtual atom (e.g., #3) of each molecule is defined as the tip of the\n");
			mprintf("    molecular dipole vector starting from molecular center of mass (i.e., #2). So you can involve the\n");
			mprintf("    dipole vector in any analysis. 1 Debye corresponds to 100 pm length of the vector.\n\n");
		}
	}
	z2 = 0;
	m_iCombinations = 0;
	do {
		ti = 1;
		if (z2 != 0)
			mprintf("\n    %d. set of vectors\n\n",z2+1);
		for (z=0;z<2;z++)
		{
			if (m_iVecType[z] == 0) // Position
			{
/*				if (z2 == 0)
					if ((z == 1) && (m_iVecType[0] == 0) && (m_iVecType[1] == 0))
						m_bSameFoot = AskYesNo("    Should the base points of the 1st and 2nd vectors always be equal (y/n)? [no] ",false);
							else m_bSameFoot = false;*/

				if (z2 == 0)
				{
					mprintf("\n");
					m_bOrtho[z] = (AskRangeInteger("    Should the %d. vector connect 2 points (0) or stand perpendicular to 3 points (1)? [0] ",0,1,0,z+1) != 0);
					mprintf("\n");
				}

				if (m_bOrtho[z])
				{
_ax1:				if (z2 == 0)
					{
						if (m_iShowMol != -1)
							m_iRefOrSec[z][0] = AskRangeInteger("    Take atom(s) at the base point from RM %s (0) or from OM %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else
							m_iRefOrSec[z][0] = 0;
					}
					mprintf("      Please enter the atom(s) at the base point (e.g. C7): ");
					inpprintf("! Please enter the atom(s) at the base point (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						delete ag;
						goto _ax1;
					}
					m_oaVectors.Add(ag);
					ti *= ag->m_iAtomGes;

_ax2:				if (z2 == 0)
					{
						if (m_iShowMol != -1)
							m_iRefOrSec[z][1] = AskRangeInteger("    Take 2nd atom(s) of normal plane from RM %s (0) or from OM %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else
							m_iRefOrSec[z][1] = 0;
					}
					mprintf("      Please enter the 2nd atom(s) of the normal plane (e.g. C7): ");
					inpprintf("! Please enter the 2nd atom(s) of the normal plane (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					if (!ag->ParseAtoms((!m_iRefOrSec[z][1])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						delete ag;
						goto _ax2;
					}
					m_oaVectors.Add(ag);
					ti *= ag->m_iAtomGes;
_ax3:				if (z2 == 0)
					{
						if (m_iShowMol != -1)
							m_iRefOrSec[z][2] = AskRangeInteger("    Take 3rd atom(s) of normal plane from RM %s (0) or from OM %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
								else m_iRefOrSec[z][2] = 0;
					}
					mprintf("      Please enter the 3rd atom(s) of the normal plane (e.g. C7): ");
					inpprintf("! Please enter the 3rd atom(s) of the normal plane (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					if (!ag->ParseAtoms((!m_iRefOrSec[z][2])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						delete ag;
						goto _ax3;
					}
					m_oaVectors.Add(ag);
					ti *= ag->m_iAtomGes;
				} else // IF ORTHO
				{
_ax4:				if (z2 == 0)
					{
						if (m_iShowMol != -1)
							m_iRefOrSec[z][0] = AskRangeInteger("    Take atom(s) at the base point from RM %s (0) or from OM %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else
							m_iRefOrSec[z][0] = 0;
					}
					mprintf("      Please enter the atom(s) at the base point (e.g. C7): ");
					inpprintf("! Please enter the atom(s) at the base point (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						delete ag;
						goto _ax4;
					}
		//			mprintf("## ParseAtoms: \"%s\" --> \"%s\".\n",(const char*)buf,ag->m_sName);
					m_oaVectors.Add(ag);
					ti *= ag->m_iAtomGes;

_ax5:				if (z2 == 0)
					{
						if (m_iShowMol != -1)
							m_iRefOrSec[z][1] = AskRangeInteger("    Take atom(s) at the tip point from RM %s (0) or from OM %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else
							m_iRefOrSec[z][1] = 0;
					}
					mprintf("      Please enter the atom(s) at the tip point (e.g. C7): ");
					inpprintf("! Please enter the atom(s) at the tip point (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					if (!ag->ParseAtoms((!m_iRefOrSec[z][1])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						delete ag;
						goto _ax5;
					}
		//			mprintf("## ParseAtoms: \"%s\" --> \"%s\".\n",(const char*)buf,ag->m_sName);
					m_oaVectors.Add(ag);
					ti *= ag->m_iAtomGes;
					m_oaVectors.Add(NULL);
				} // END IF NOT ORTHO
			} else if (m_iVecType[z] == 1) // Dipol
			{
				mprintf("\n");
				if (m_iShowMol != -1)
					m_iRefOrSec[z][0] = AskRangeInteger("    Take dipole vector from RM %s (0) or from OM %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
				else
					m_iRefOrSec[z][0] = 0;

				m_bOrtho[z] = false;

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],"#2"))
				{
					eprintf("CADF::Parse(): Weird error. \"#2\" not found.\n");
					abort();
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				buf.Format("#%d",((CMolecule*)((!m_iRefOrSec[z][0])?g_oaMolecules[g_iFixMol]:g_oaMolecules[m_iShowMol]))->m_laVirtualAtoms.GetSize());

				if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("CADF::Parse(): Weird error. \"%s\" not found.\n",(const char*)buf);
					abort();
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
				m_oaVectors.Add(NULL);

			} else if (m_iVecType[z] == 2) // Geschwindigkeit
			{
_ax6:			if (m_iShowMol != -1)
					m_iRefOrSec[z][0] = AskRangeInteger("    Take velocity vector from RM %s (0) or from OM %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else m_iRefOrSec[z][0] = 0;
				mprintf("      Velocity vector of which atoms to use (e.g. C7)? [#2] ");
				inpprintf("! Velocity vector of which atoms to use (e.g. C7)? [#2]\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (strlen(buf)==0)
				{
					if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],"#2"))
					{
						eprintf("Weird error.\n");
						inpprintf("! Weird error.\n");
						abort();
					}
				} else if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					delete ag;
					goto _ax6;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
			} else if (m_iVecType[z] == 3) // Kraft
			{
_ax7:			
				if (m_iShowMol != -1)
					m_iRefOrSec[z][0] = AskRangeInteger("    Take force vector from RM %s (0) or from OM %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else m_iRefOrSec[z][0] = 0;
				mprintf("      Force vector of which atoms to use (e.g. C7)? [#2] ");
				inpprintf("! Force vector of which atoms to use (e.g. C7)? [#2]\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (strlen(buf)==0)
				{
					if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],"#2"))
					{
						eprintf("Weird Error.\n");
						inpprintf("! Weird Error.\n");
						abort();
					}
				} else if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					delete ag;
					goto _ax7;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
			}
		} // END FOR 0..1
		z2++;
		m_iCombinations += ti;
	//	mprintf("RefOrSec: %d %d %d %d.\n",m_iRefOrSec[0][0],m_iRefOrSec[0][1],m_iRefOrSec[1][0],m_iRefOrSec[1][1]);
	} while (g_bAdvanced2?AskYesNo("\n    Add another set of vectors to this (!) ADF (y/n)? [no] ",false):false);

	ParseDeriv();

	switch(m_iDeriv)
	{
		case 0:
			mprintf("\n");
_anglowag:
			m_fMinAngle = AskFloat("    Enter minimal angle between the vectors: [0 deg] ",0.0f);
			if (m_fMinAngle < 0)
			{
				eprintf("    Angles < 0 degree are not defined.\n");
				goto _anglowag;
			}
			m_fMaxAngle = AskFloat("    Enter maximal angle between the vectors: [180 deg] ",180.0f);
			break;
		case 1:
			if (m_bDerivAbs)
				m_fMinAngle = AskFloat("    Enter the minimal value of this d1-ADF in deg/ps: [0] ",0.0f);
					else m_fMinAngle = AskFloat("    Enter the minimal value of this d1-ADF in deg/ps: [-10.0] ",-10.0f);
			m_fMaxAngle = AskFloat("    Enter the maximal value of this d1-ADF in deg/ps: [10.0] ",10.0f);
			break;
		case 2:
			if (m_bDerivAbs)
				m_fMinAngle = AskFloat("    Enter the minimal value of this d2-ADF in deg/ps^2: [0] ",0.0f);
					else m_fMinAngle = AskFloat("    Enter the minimal value of this d2-ADF in deg/ps^2: [-10.0] ",-10.0f);
			m_fMaxAngle = AskFloat("    Enter the maximal value of this d2-ADF in deg/ps^2: [10.0] ",10.0f);
			break;
	}

/*	if (m_fMaxAngle <= 90.0f)
		m_bFoldAngle = AskYesNo("    Should angles > 90 deg be \"mirrored\" (180 deg = 0 deg) (y/n)? [yes] ",true);
			else*/ m_bFoldAngle = false;
	m_bCosine = (AskRangeInteger("    Plot ADF against angle (0) or against cosine (1)? [0] ",0,1,0)!=0);
	if (m_bCosine)
	{
		m_fMinAngle = (float)cos(m_fMinAngle/180.0*Pi);
		m_fMaxAngle = (float)cos(m_fMaxAngle/180.0*Pi);
		if (m_fMinAngle > m_fMaxAngle)
		{
			tf = m_fMinAngle;
			m_fMinAngle = m_fMaxAngle;
			m_fMaxAngle = tf;
		}
		mprintf("    The data range is %.2f to %.2f.\n",m_fMinAngle,m_fMaxAngle);
	}
	if (!m_bCosine)
		m_bMirror = AskYesNo("    Force this ADF to be mirror-symmetric to the 90 deg line (y/n)? [no] ",false);
			else m_bMirror = false;
	m_iResolution = AskUnsignedInteger("    Please enter the resolution (bin count) for this ADF: [100] ",100);

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
			else m_iHistogramRes = 0;

	if (m_iShowMol != -1)
		m_bStat = AskYesNo("    Apply cone correction (y/n)? [%c] ",!m_bCosine,(!m_bCosine)?'y':'n');
			else m_bStat = 0;
/*	mprintf("\n    Save temporal development of this ADF (0=nein, 1=ja)? [0] ");
	myget(buf);
	m_bSaveAngle = (atoi(buf)!=0);*/
	BuildName();
	mprintf(WHITE,"\n<<< End of Angular Distribution Function <<<\n\n");
	BTOUT;
}


void CDDF::Parse()
{
	BTIN;
//	char buf[1024];
	CxString buf;
	CAtomGroup *ag;
	int z0, z, z2, i;
	float tf;

	try { m_pDDF = new CDF(); } catch(...) { m_pDDF = NULL; }
	if (m_pDDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Dihedral Distribution Function >>>\n\n");
	m_bClassical = AskYesNo("    Use \"simple\" (y) (4 atoms) or \"generalized\" (n) (3 vectors) Dihedrals? [yes] ",true);

	m_iCombinations = 1;
	if (m_bClassical)
	{
		m_bOrtho[0] = false;
		m_bOrtho[1] = false;
		m_bOrtho[2] = false;
		z0 = m_oaVectors.GetSize();
		for (z=0;z<3;z++)
		{
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			m_oaVectors.Add(ag);

			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			m_oaVectors.Add(ag);
			m_oaVectors.Add(NULL);
		}
		for (z=0;z<4;z++)
		{
_bx:		if (m_iShowMol != -1)
				i = AskRangeInteger("    Take the %d. atom from RM %s (0) or from OM %s (1)? [0] ",0,1,0,z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
					else i = 0;
			switch(z)
			{
				case 0:
					m_iRefOrSec[0][1] = i;
					break;
				case 1:
					m_iRefOrSec[0][0] = i;
					m_iRefOrSec[2][0] = i;
					break;
				case 2:
					m_iRefOrSec[2][1] = i;
					m_iRefOrSec[1][0] = i;
					break;
				case 3:
					m_iRefOrSec[1][1] = i;
					break;
			}
			mprintf("      Enter the %d. atom(s) (e.g. C7): ",z+1);
			inpprintf("! Enter the %d. atom(s) (e.g. C7):\n",z+1);
			myget(&buf);
			switch(z)
			{
				case 0:
					if (!((CAtomGroup*)m_oaVectors[z0+0*3+1])->ParseAtoms((!m_iRefOrSec[0][1])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _bx;
					}
					m_iCombinations *= ((CAtomGroup*)m_oaVectors[z0+0*3+1])->m_iAtomGes;
					break;
				case 1:
					if (!((CAtomGroup*)m_oaVectors[z0+0*3+0])->ParseAtoms((!m_iRefOrSec[0][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _bx;
					}
					if (!((CAtomGroup*)m_oaVectors[z0+2*3+0])->ParseAtoms((!m_iRefOrSec[2][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _bx;
					}
					m_iCombinations *= ((CAtomGroup*)m_oaVectors[z0+0*3+0])->m_iAtomGes;
					m_iCombinations *= ((CAtomGroup*)m_oaVectors[z0+0*3+0])->m_iAtomGes;
					break;
				case 2:
					if (!((CAtomGroup*)m_oaVectors[z0+2*3+1])->ParseAtoms((!m_iRefOrSec[2][1])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _bx;
					}
					if (!((CAtomGroup*)m_oaVectors[z0+1*3+0])->ParseAtoms((!m_iRefOrSec[1][0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _bx;
					}
					m_iCombinations *= ((CAtomGroup*)m_oaVectors[z0+2*3+1])->m_iAtomGes;
					m_iCombinations *= ((CAtomGroup*)m_oaVectors[z0+2*3+1])->m_iAtomGes;
					break;
				case 3:
					if (!((CAtomGroup*)m_oaVectors[z0+1*3+1])->ParseAtoms((!m_iRefOrSec[1][1])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _bx;
					}
					m_iCombinations *= ((CAtomGroup*)m_oaVectors[z0+1*3+1])->m_iAtomGes;
					break;
			}
		}
	} else // NONCLASSIC
	{
		mprintf("\nYou now have to define 3 vectors:\n");
		mprintf("The 1st and 2nd vector are projected onto the normal plane of the 3rd vector.\n");
		mprintf("The angle between the two projected vectors in the plane is evaluated.\n");
		for (z=0;z<3;z++)
		{
			mprintf(WHITE,"\n  * Vector %d\n",z+1);
			m_bOrtho[z] = (AskRangeInteger("    Shall the %d. vector connect 2 points (0) or be orthogonal to a plane (1)? [0] ",0,1,0,z+1) != 0);
			if (m_bOrtho[z])
			{
				for (z2=0;z2<3;z2++)
				{
					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

_by:				if (m_iShowMol != -1)
						m_iRefOrSec[z][z2] = AskRangeInteger("    Take the %d. atom(s) of the plane from RM %s (0) or from OM %s (1)? [0] ",0,1,0,z2+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
							else m_iRefOrSec[z][z2] = 0;
					mprintf("      Please enter the %d. atom(s) (e.g. C7): ",z2+1);
					inpprintf("! Please enter the %d. atom(s) (e.g. C7):\n",z2+1);
					myget(&buf);
					if (!ag->ParseAtoms((!m_iRefOrSec[z][z2])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _by;
					}
					m_iCombinations *= ag->m_iAtomGes;
					m_oaVectors.Add(ag);
				}
			} else // Ortho
			{
				for (z2=0;z2<2;z2++)
				{
					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

_bz:				if (m_iShowMol != -1)
						m_iRefOrSec[z][z2] = AskRangeInteger("    Take the %d. atom(s) of the vector from RM %s (0) or from OM %s (1)? [0] ",0,1,0,z2+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
							else m_iRefOrSec[z][z2] = 0;
					mprintf("      Please enter the %d. atom(s) (e.g. C7): ",z2+1);
					inpprintf("! Please enter the %d. atom(s) (e.g. C7):\n",z2+1);
					myget(&buf);
					if (!ag->ParseAtoms((!m_iRefOrSec[z][z2])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _bz;
					}
					m_iCombinations *= ag->m_iAtomGes;
					m_oaVectors.Add(ag);
				}
				m_oaVectors.Add(NULL);
			} // ENDIF ORTHO
		} // END FOR
	} // END IF NONCLASSIC

/*	mprintf("\n");
	mprintf("    The \"classical\" dihedral angle is defined for a range of 0 ... 180 deg.\n\n");
	m_bAbs = !AskYesNo("    Do you want to extend the range to -180 ... +180 deg (y/n)? [no] ",false);*/

	mprintf("\n");
	mprintf("    Per default, the dihedral angle is defined in a range of -180 ... 180 deg.\n\n");
	m_bPositive = AskYesNo("    Use range of 0 ... 360 deg instead (e.g. -90 deg becomes 270 deg, ...) (y/n)? [no] ",false);
	m_bAbs = false;
	mprintf("\n");

	ParseDeriv();

	switch(m_iDeriv)
	{
		case 0:
			if (m_bPositive)
			{
				m_fMinAngle = AskFloat("\n    Enter minimal dihedral angle to observe:   [0 deg] ",0);
				m_fMaxAngle = AskFloat("    Enter maximal dihedral angle to observe: [360 deg] ",360.0f);
			} else
			{
				m_fMinAngle = AskFloat("\n    Enter minimal dihedral angle to observe: [-180 deg] ",-180.0f);
				m_fMaxAngle = AskFloat("    Enter maximal dihedral angle to observe:  [180 deg] ",180.0f);
			}
			break;
		case 1:
			if (m_bDerivAbs)
				m_fMinAngle = AskFloat("    Enter the minimal value of this d1-DDF in deg/ps: [0] ",0.0f);
					else m_fMinAngle = AskFloat("    Enter the minimal value of this d1-DDF in deg/ps: [-10.0] ",-10.0f);
			m_fMaxAngle = AskFloat("    Enter the maximal value of this d1-DDF in deg/ps: [10.0] ",10.0f);
			break;
		case 2:
			if (m_bDerivAbs)
				m_fMinAngle = AskFloat("    Enter the minimal value of this d2-DDF in deg/ps^2: [0] ",0.0f);
					else m_fMinAngle = AskFloat("    Enter the minimal value of this d2-DDF in deg/ps^2: [-10.0] ",-10.0f);
			m_fMaxAngle = AskFloat("    Enter the maximal value of this d2-DDF in deg/ps^2: [10.0] ",10.0f);
			break;
	}

	if (m_iDeriv == 0)
	{
		m_bCosine = (AskRangeInteger("    Plot DDF against angle (0) or against cosine (1)? [0] ",0,1,0) != 0);
		if (m_bCosine)
		{
			if ((m_fMinAngle <= 0) && (m_fMaxAngle >= 0))
			{
				if (fabsf(m_fMinAngle) > fabsf(m_fMaxAngle))
					m_fMinAngle = cosf(fabsf(m_fMinAngle));
						else m_fMinAngle = cosf(fabsf(m_fMaxAngle));
				m_fMaxAngle = 1.0f;
			} else
			{
				m_fMinAngle = cosf(m_fMinAngle/180.0f*Pi);
				m_fMaxAngle = cosf(m_fMaxAngle/180.0f*Pi);
			}
			if (m_fMinAngle > m_fMaxAngle)
			{
				tf = m_fMinAngle;
				m_fMinAngle = m_fMaxAngle;
				m_fMaxAngle = tf;
			}
			m_bAbs = false;
			m_bSymm = false;
			mprintf("\n    The data range is %.2f to %.2f.\n\n",m_fMinAngle,m_fMaxAngle);
		} else
		{
/*			m_bAbs = AskYesNo("    Use absolute angle values (y) or signed values (n)? [yes] ",true);
			if (!m_bAbs)
			{
				m_fMinAngle = -m_fMaxAngle;
				mprintf("\n    The data range is %.2f to %.2f.\n\n",m_fMinAngle,m_fMaxAngle);
				m_bSymm = AskYesNo("    Force this DDF to be symmetrical to the 0 degree line (y/n)? [no] ",false);
			} else m_bSymm = false;*/
//			if (!m_bAbs)
//			{
//				m_fMinAngle = -m_fMaxAngle;
				mprintf("\n    The data range is %.2f to %.2f.\n\n",m_fMinAngle,m_fMaxAngle);
				if (!m_bPositive)
					m_bSymm = AskYesNo("    Add also bin entries for inverse angles (force symmetry) (y/n)? [no] ",false);
						else m_bSymm = false;
//			} else m_bSymm = false;
		}
	} else 
	{
		m_bCosine = false;
		m_bAbs = false;
		m_bSymm = false;
	}
	m_iResolution = AskUnsignedInteger("    Please enter binning resolution for this DDF: [100] ",100);

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
			else m_iHistogramRes = 0;

/*	mprintf("\n    Save temporal Development of Dihedrals (0=no, 1=yes)? [0] ");
	myget(buf);
	m_bSaveAngle = atoi(buf)!=0;*/
	BuildName();
	mprintf(WHITE,"\n<<< End of Dihedral Distribution Function <<<\n\n");
	BTOUT;
}


void CPlDF::Parse()
{
	BTIN;
//	char buf[1024];
	CxString buf;
	CAtomGroup *ag;
	int z;

	try { m_pPlDF = new CDF(); } catch(...) { m_pPlDF = NULL; }
	if (m_pPlDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Plane Distance Distribution Function >>>\n\n");

	m_bNormal = !AskYesNo("    Define plane via 3 atoms (y) or as normal plane of a vector (n) (y/n)? [yes] ",true);

	m_iCombinations = 1;
	if (m_bNormal)
	{
		try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
		if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
_b2:	
		if (m_iShowMol != -1)
			m_iRefOrSec[0] = AskRangeInteger("    Take the base atom of the vector/plane from RM %s (0) or from OM %s (1)? [0] ",0,1,0,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
				else m_iRefOrSec[0] = 0;
		mprintf("      Please enter the base atom (e.g. C7): ");
		inpprintf("! Please enter the base atom (e.g. C7):\n");
		myget(&buf);
		if (!ag->ParseAtoms((!m_iRefOrSec[0])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
		{
			eprintf("Wrong input.\n");
			inpprintf("! Wrong input.\n");
			goto _b2;
		}
		m_iCombinations *= ag->m_iAtomGes;
		m_oaVectors.Add(ag);

		try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
		if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
_b3:	
		if (m_iShowMol != -1)
			m_iRefOrSec[1] = AskRangeInteger("    Take the atom at the tip of the vector from RM %s (0) or from OM %s (1)? [0] ",0,1,0,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
				else m_iRefOrSec[1] = 0;
		mprintf("      Please enter the tip atom (e.g. C7): ");
		inpprintf("! Please enter the tip atom (e.g. C7):\n");
		myget(&buf);
		if (!ag->ParseAtoms((!m_iRefOrSec[1])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
		{
			eprintf("Wrong input.\n");
			inpprintf("! Wrong input.\n");
			goto _b3;
		}
		m_iCombinations *= ag->m_iAtomGes;
		m_oaVectors.Add(ag);
		m_oaVectors.Add(NULL); // Dummy placeholder
	} else
	{
		for (z=0;z<3;z++)
		{
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
_b1:	
			if (m_iShowMol != -1)
				m_iRefOrSec[z] = AskRangeInteger("    Take the %d. atom(s) of the plane from RM %s (0) or from OM %s (1)? [0] ",0,1,0,z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
					else m_iRefOrSec[z] = 0;
			mprintf("      Please enter the %d. atom (e.g. C7): ",z+1);
			inpprintf("! Please enter the %d. atom (e.g. C7):\n",z+1);
			myget(&buf);
			if (!ag->ParseAtoms((!m_iRefOrSec[z])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
			{
				eprintf("Wrong input.\n");
				inpprintf("! Wrong input.\n");
				goto _b1;
			}
			m_iCombinations *= ag->m_iAtomGes;
			m_oaVectors.Add(ag);
		}
	}
	mprintf("    Plane defined.\n\n");

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
_b4:	
	if (m_iShowMol != -1)
		m_iRefOrSec[3] = AskRangeInteger("    Take the atom(s) to observe from RM %s (0) or from OM %s (1)? [0] ",0,1,0,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			else m_iRefOrSec[3] = 0;
	mprintf("      Please enter the atom(s) to observe (e.g. C7): ");
	inpprintf("! Please enter the atom(s) to observe (e.g. C7):\n");
	myget(&buf);
	if (!ag->ParseAtoms((!m_iRefOrSec[3])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
	{
		eprintf("Wrong input.\n");
		inpprintf("! Wrong input.\n");
		goto _b4;
	}
	m_iCombinations *= ag->m_iAtomGes;
	m_oaVectors.Add(ag);

	m_fMinDist = AskFloat("    Enter the minimal distance to observe (in pm): [-%d.0] ",(float)-HalfBox(),HalfBox());
	m_fMaxDist = AskFloat("    Enter the maximal distance to observe (in pm): [%d.0] ",(float)HalfBox(),HalfBox());

	m_iResolution = AskUnsignedInteger("    Please enter binning resolution for this PlDF: [100] ",100);

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
			else m_iHistogramRes = 0;

	BuildName();
	mprintf(WHITE,"\n<<< End of Plane Distance Distribution Function <<<\n\n");
	BTOUT;
}


void CLiDF::Parse()
{
	BTIN;
//	char buf[1024];
	CxString buf;
	CAtomGroup *ag;
	int z;

	try { m_pLiDF = new CDF(); } catch(...) { m_pLiDF = NULL; }
	if (m_pLiDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Line Distance Distribution Function >>>\n\n");

	m_bNormal = !AskYesNo("    Define line via 2 atoms (y) or as normal vector of a plane (n) (y/n)? [yes] ",true);

	m_iCombinations = 1;
	if (m_bNormal)
	{
		mprintf("\n    The line will pass through the 1st point of the plane.\n\n");
		for (z=0;z<3;z++)
		{
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	_b3:	
			if (m_iShowMol != -1)
				m_iRefOrSec[z] = AskRangeInteger("    Take the %d. atom(s) of the plane from RM %s (0) or from OM %s (1)? [0] ",0,1,0,z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
					else m_iRefOrSec[z] = 0;
			mprintf("      Please enter the %d. atom (e.g. C7): ",z+1);
			inpprintf("! Please enter the %d. atom (e.g. C7):\n",z+1);
			myget(&buf);
			if (!ag->ParseAtoms((!m_iRefOrSec[z])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
			{
				eprintf("Wrong input.\n");
				inpprintf("! Wrong input.\n");
				goto _b3;
			}
			m_iCombinations *= ag->m_iAtomGes;
			m_oaVectors.Add(ag);
		}
	} else
	{
		for (z=0;z<2;z++)
		{
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	_b1:	
			if (m_iShowMol != -1)
				m_iRefOrSec[z] = AskRangeInteger("    Take the %d. atom(s) of the line from RM %s (0) or from OM %s (1)? [0] ",0,1,0,z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
					else m_iRefOrSec[z] = 0;
			mprintf("      Please enter the %d. atom (e.g. C7): ",z+1);
			inpprintf("! Please enter the %d. atom (e.g. C7):\n",z+1);
			myget(&buf);
			if (!ag->ParseAtoms((!m_iRefOrSec[z])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
			{
				eprintf("Wrong input.\n");
				inpprintf("! Wrong input.\n");
				goto _b1;
			}
			m_iCombinations *= ag->m_iAtomGes;
			m_oaVectors.Add(ag);
		}
		m_oaVectors.Add(NULL);
	}

	mprintf("    Line defined.\n\n");

	try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
	if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
_b2:	
	if (m_iShowMol != -1)
		m_iRefOrSec[3] = AskRangeInteger("    Take the atom(s) to observe from RM %s (0) or from OM %s (1)? [0] ",0,1,0,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			else m_iRefOrSec[3] = 0;
	mprintf("      Please enter the atom(s) to observe (e.g. C7): ");
	inpprintf("! Please enter the atom(s) to observe (e.g. C7):\n");
	myget(&buf);
	if (!ag->ParseAtoms((!m_iRefOrSec[3])?(CMolecule*)g_oaMolecules[g_iFixMol]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
	{
		eprintf("Wrong input.\n");
		inpprintf("! Wrong input.\n");
		goto _b2;
	}
	m_iCombinations *= ag->m_iAtomGes;
	m_oaVectors.Add(ag);

	m_fMinDist = AskFloat("    Enter the minimal distance to observe (in pm): [0] ",0.0f);
	m_fMaxDist = AskFloat("    Enter the maximal distance to observe (in pm): [%d.0] ",(float)HalfBox(),HalfBox());

	m_iResolution = AskUnsignedInteger("    Please enter binning resolution for this LiDF: [100] ",100);

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
			else m_iHistogramRes = 0;

	if (m_iShowMol != -1)
		m_bRadialCorrect = AskYesNo("    Correct radial distribution for this LiDF (y/n)? [yes] ",true);
			else m_bRadialCorrect = AskYesNo("    Correct radial distribution for this LiDF (y/n)? [no] ",false);

	BuildName();
	mprintf(WHITE,"\n<<< End of Line Distance Distribution Function <<<\n\n");
	BTOUT;
}


void CDipDF::Parse()
{
	BTIN;
	float tf;
	double td, td2, td3;
	CMolecule *m;
	CSingleMolecule *sm;
	int z;

	g_bDipole = true;

	try { m_pDipoleDF = new CDF(); } catch(...) { m_pDipoleDF = NULL; }
	if (m_pDipoleDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Dipole Distribution Function >>>\n\n");
	if (m_iShowMol != -1)
		m_iRefOrSec = AskRangeInteger("    Observe dipole moment from RM %s (0) or from OM %s (1)? [0] ",0,1,0,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			else m_iRefOrSec = 0;

	ParseDeriv();

	if (m_iRefOrSec != 0)
		m = (CMolecule*)g_oaMolecules[m_iShowMol];
			else m = (CMolecule*)g_oaMolecules[g_iFixMol];

	td = 0;
	td2 = 1.0e10;
	td3 = 0;
	for (z=0;z<m->m_laSingleMolIndex.GetSize();z++)
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z]];
		tf = sm->m_vDipole.GetLength();
		if (tf > td)
			td = tf;
		if (tf < td2)
			td2 = tf;
		td3 += tf;
	}
	td3 /= (double)m->m_laSingleMolIndex.GetSize();

	mprintf("\n    Dipole moment of %s (1st step): Min. %.3f, Max. %.3f, Avg. %.3f Debye.\n\n",m->m_sName,td2,td,td3);

	switch(m_iDeriv)
	{
		case 0:
			m_fDipoleMin = AskFloat("    Enter the lower bound for dipole values (in Debye): [0] ",0.0f);
			m_fDipoleMax = AskFloat("    Enter the upper bound for dipole values (in Debye): [%d] ",float(int(td*2.0)),int(td*2.0));
			break;
		case 1:
			if (m_bDerivAbs)
				m_fDipoleMin = AskFloat("    Enter the minimal value of this d1-DDF in Debye/ps: [0] ",0.0f);
					else m_fDipoleMin = AskFloat("    Enter the minimal value of this d1-DDF in Debye/ps: [-10.0] ",-10.0f);
			m_fDipoleMax = AskFloat("    Enter the maximal value of this d1-DDF in Debye/ps: [10.0] ",10.0f);
			break;
		case 2:
			if (m_bDerivAbs)
				m_fDipoleMin = AskFloat("    Enter the minimal value of this d2-DDF in Debye/ps^2: [0] ",0.0f);
					else m_fDipoleMin = AskFloat("    Enter the minimal value of this d2-DDF in Debye/ps^2: [-10.0] ",-10.0f);
			m_fDipoleMax = AskFloat("    Enter the maximal value of this d2-DDF in Debye/ps^2: [10.0] ",10.0f);
			break;
	}

	m_iResolution = AskUnsignedInteger("    Enter binning resolution for the Dipole DF: [%d] ",int((m_fDipoleMax-m_fDipoleMin)*10.0),int((m_fDipoleMax-m_fDipoleMin)*10.0));

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
			else m_iHistogramRes = 0;

	BuildName();
	mprintf(WHITE,"\n<<< End of Dipole Distribution Function <<<\n\n");
	BTOUT;
}


void CSDF::Parse(bool voro)
{
	BTIN;
	int ti;
//	char buf[1024];
	CxString buf;

	m_iShowAtomGes = 0;

	try { m_pSDF = new C3DF<double>(); } catch(...) { m_pSDF = NULL; }
	if (m_pSDF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Spatial Distribution Function >>>\n\n");
_sdfatoms:
	if (!voro)
	{
		if (m_bIntra)
		{
			mprintf("    Observing atoms in reference molecule %s.\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
			m_iRefOrSec = 0;
		} else
		{
			mprintf("    Observing atoms in observed molecule %s.\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			m_iRefOrSec = 1;
		}
	/*	if (m_iShowMol != -1)
			m_iRefOrSec = AskRangeInteger("    Observe atoms in RM %s (0) or in OM %s (1)? [1] ",0,1,1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
				else m_iRefOrSec = 0;*/
		mprintf("    Which atoms to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ");
		inpprintf("! Which atoms to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n");
		myget(&buf);
		if (strlen(buf) == 0)
		{
			if (!m_oAtoms.ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec==0)?g_iFixMol:m_iShowMol],"#2"))
			{
				eprintf("Strange error ^^\n");
				inpprintf("! Strange error ^\n");
				goto _sdfatoms;
			}
		} else 
		{
			if (!m_oAtoms.ParseAtoms((CMolecule*)g_oaMolecules[(m_iRefOrSec==0)?g_iFixMol:m_iShowMol],buf))
			{
				eprintf("Wrong input.\n");
				inpprintf("! Wrong input.\n");
				goto _sdfatoms;
			}
		}

		m_iShowAtomGes += m_oAtoms.m_iAtomGes;
		m_fParticleDensity = m_oAtoms.m_iAtomGes * ((CMolecule*)g_oaMolecules[(m_iRefOrSec==0)?g_iFixMol:m_iShowMol])->m_laSingleMolIndex.GetSize() / g_fBoxX / g_fBoxY / g_fBoxZ * 1E9f;
	}

	m_fRadius = AskFloat("    Please enter radius of this SDF in pm: [%d.0] ",(float)HalfBox(),HalfBox());

	if (g_bPeriodic && (m_fRadius > HalfBox_Exact()+1.0f))
	{
		eprintf("\nWarning: ");
		mprintf("The specified max. radius is larger than half of the smallest periodic cell vector.\n");
		mprintf("         TRAVIS counts every atom only once (central periodic image).\n");
		mprintf("         Expect the analysis to decay to zero for large radii.\n\n");
		AskYesNo("         Acknowledged [yes] ",true);
		mprintf("\n");
	}

_sdfresagain:
	ti = g_pDatabase->GetInt("/PLOT3D/DEFAULTS/BIN_RES");
	m_iResolution = AskUnsignedInteger("    Please enter binning resolution of this SDF per dimension: [%d] ",ti,ti);

	if (m_iResolution > 300)
	{
		eprintf("\nWarning: ");
		mprintf("This large resolution will consume much RAM (%s)\n",FormatBytes(pow((double)m_iResolution,3)*sizeof(double)));
		mprintf("         and result in very large SDF output files.\n\n");
		if (!AskYesNo("    Use this resolution (y/n)? [no] ",false))
			goto _sdfresagain;
	}

	if (g_bAdvanced2)
	{
		m_bVdWSpheres = false; //AskYesNo("    Process atoms as VdW spheres rather than points (y/n)? [no] ",false);
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
		m_bClipPlane = AskYesNo("    Use a clipping plane for this SDF (y/n)? [no] ",false);
		if (m_bClipPlane)
		{
	_cut:	mprintf("    Should the clipping plane stand perpendicular to the X, Y or Z axis (x/y/z)? [x] ");
			inpprintf("! Should the clipping plane stand perpendicular to the X, Y or Z axis (x/y/z)? [x]\n");
			myget(&buf);
			switch(tolower(buf[0]))
			{
				case 0:
				case 'x':
					m_iClipDirection = 0;
					break;
				case 'y':
					m_iClipDirection = 1;
					break;
				case 'z':
					m_iClipDirection = 2;
					break;
				default:
					eprintf("    Wrong input! Please enter 'x', 'y' or 'z'.\n");
					inpprintf("! Wrong input! Please enter 'x', 'y' or 'z'.\n");
					goto _cut;	
			}
			m_fClipValue = AskFloat("    Please enter distance of clipping plane from origin (in pm): [0] ",0.0f);
		}
		m_bInvert = AskYesNo("    Should this SDF be inverted (y/n)? [no] ",false);
		m_bSDFMirrorXY = AskYesNo("    Force this SDF to be symmetrical to the XY plane (y/n)? [no] ",false);
		m_bSDFMirrorBisect = AskYesNo("    Force this SDF to be symmetrical to the angle bisector (y/n)? [no] ",false);
		m_bCutPlane = AskYesNo("    Add a particle density cut plane to this SDF (y/n)? [no] ",false);
		if (m_bCutPlane)
		{
			m_iCutPlaneResolution = AskUnsignedInteger("    Enter the binning resolution for the cut plane per dimension: [100] ",100);
			m_bCutPlaneShowAtoms = AskYesNo("    Show reference atoms in Pseudo SDF plot (y/n)? [yes] ",true);
		}
	} else
	{
		m_bVdWSpheres = false;
		m_bClipPlane = false;
		m_bInvert = false;
		m_bSDFMirrorXY = false;
		m_bSDFMirrorBisect = false;
		m_bCutPlane = false;
		m_iHistogramRes = 0;
	}

	if (!voro)
		BuildName();

	mprintf(WHITE,"\n<<< End of Spatial Distribution Function <<<\n\n");
	BTOUT;
}


void CVDF::Parse()
{
	BTIN;
//	char buf[1024];
	CxString buf;

	m_iShowAtomGes = 0;

	try { m_pVDF = new CDF(); } catch(...) { m_pVDF = NULL; }
	if (m_pVDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n>>> Velocity Distribution Function >>>\n\n");
_atoms:
	if (g_iFixMol == -1)
	{
		mprintf("    Observing atoms from OM %s.\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		m_iRefOrSec = 1;
	} else if (m_iShowMol == -1)
	{
		mprintf("    Observing atoms from RM %s.\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		m_iRefOrSec = 0;
	} else m_iRefOrSec = AskRangeInteger("    Observe atoms in RM %s (0) or in OM %s (1)? [1]",0,1,1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);

	mprintf("    Which atoms to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ");
	inpprintf("! Which atoms to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n");
	myget(&buf);
	if (strlen(buf) == 0)
	{
		m_oAtoms.AddAllAtoms((CMolecule*)g_oaMolecules[m_iRefOrSec?m_iShowMol:g_iFixMol],false);
		if (!m_oAtoms.ParseAtoms((CMolecule*)g_oaMolecules[m_iRefOrSec?m_iShowMol:g_iFixMol],"#2"))
		{
			eprintf("Strange error ^^\n");
			inpprintf("! Strange error ^^\n");
			goto _atoms;
		}
	} else 
	{
		if (!m_oAtoms.ParseAtoms((CMolecule*)g_oaMolecules[m_iRefOrSec?m_iShowMol:g_iFixMol],buf))
		{
			eprintf("Wrong input.\n");
			inpprintf("! Wrong input.\n");
			goto _atoms;
		}
	}
	m_iCombinations = m_oAtoms.m_iAtomGes;
	m_iShowAtomGes += m_oAtoms.m_iAtomGes;
//	m_fParticleDensity = m_oAtoms.m_iAtomGes * ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize() / g_fBoxX / g_fBoxY / g_fBoxZ;

	ParseDeriv();

	switch(m_iDeriv)
	{
		case 0:
		/*	if (AskYesNo("    Scan for velocity range (y) or enter range manually (n)? [no] ",false))
			{
				m_bScanRange = true;
				g_bScanVelocities = true;
			} else*/
			{
				m_bScanRange = false;
				m_fMinSpeed = 0; //AskFloat("    Please enter the minimal velocity: [0.0 pm/fs] ",0.0f);
				m_fMaxSpeed = AskFloat("    Please enter the maximal velocity: [10000.0 pm/ps] ",10000.0f);
			}
			break;
		case 1:
			if (m_bDerivAbs)
				m_fMinSpeed = AskFloat("    Enter the minimal value of this d1-VDF in pm/ps^2: [0] ",0.0f);
					else m_fMinSpeed = AskFloat("    Enter the minimal value of this d1-VDF in pm/ps^2: [-10.0] ",-10.0f);
			m_fMaxSpeed = AskFloat("    Enter the maximal value of this d1-VDF in pm/ps^2: [10.0] ",10.0f);
			break;
		case 2:
			if (m_bDerivAbs)
				m_fMinSpeed = AskFloat("    Enter the minimal value of this d2-VDF in pm/ps^3: [0] ",0.0f);
					else m_fMinSpeed = AskFloat("    Enter the minimal value of this d2-VDF in pm/ps^3: [-10.0] ",-10.0f);
			m_fMaxSpeed = AskFloat("    Enter the maximal value of this d2-VDF in pm/ps^3: [10.0] ",10.0f);
			break;
	}

	m_iResolution = AskInteger("    Please enter the binning resolution for this VDF: [300] ",300);

	if (g_bAdvanced2)
		m_iHistogramRes = AskUnsignedInteger("    Please enter histogram resolution (0=no histogram): [5000] ",5000);
			else m_iHistogramRes = 0;

/*	mprintf("    Save temporal Development of Velocity (0=no, 1=yes)? [0] ");
	myget(buf);
	m_bSaveSpeed = atoi(buf)!=0;*/
	BuildName();
	mprintf(WHITE,"\n<<< End of Velocity Distribution Function <<<\n\n");
	BTOUT;
}


void CSDF::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;

//	tmp[0] = 0;
	tmp.sprintf("");

	if (m_iRefOrSec)
//		sprintf(tmp,"%s_%s%d%s%d%s%d_%s_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_oAtoms.m_sName);
		tmp.sprintf("%s_%s%d%s%d%s%d_%s_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_oAtoms.m_sName);
	else
//		sprintf(tmp,"%s_%s%d%s%d%s%d_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,m_oAtoms.m_sName);
		tmp.sprintf("%s_%s%d%s%d%s%d_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,m_oAtoms.m_sName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CVDF::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;

	if (m_iDeriv != 0)
//		sprintf(tmp,"deriv%d_%s",m_iDeriv,m_oAtoms.m_sName);
		tmp.sprintf("deriv%d_%s",m_iDeriv,m_oAtoms.m_sName);
	else
//		sprintf(tmp,"%s",m_oAtoms.m_sName);
		tmp.sprintf("%s",m_oAtoms.m_sName);

	try { m_sShortName = new char[strlen(tmp)+1]; } catch(...) { m_sShortName = NULL; }
	if (m_sShortName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sShortName,m_oAtoms.m_sName);

	if (g_iFixMol == -1)
	{
//		sprintf(tmp,"%s_",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		tmp.sprintf("%s_",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	} else
		if (m_iRefOrSec)
//			sprintf(tmp,"%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			tmp.sprintf("%s_%s%d_%s_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		else
//			sprintf(tmp,"%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);
			tmp.sprintf("%s_%s%d_",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1);

//	strcat(tmp,m_sShortName);
	tmp.strcat(m_sShortName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);


	try { m_sLabelName = new char[strlen(tmp)+1]; } catch(...) { m_sLabelName = NULL; }
	if (m_sLabelName == NULL) NewException((double)(strlen(m_sShortName)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sLabelName,m_sShortName);

	BTOUT;
}


void CSDF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
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


void CVDF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z1t, z1a;
	CxIntArray *a1;

	vec->RemoveAll_KeepSize();
	for (z1t=0;z1t<m_oAtoms.m_baAtomType.GetSize();z1t++)
	{
		a1 = (CxIntArray*)m_oAtoms.m_oaAtoms[z1t];
		for (z1a=0;z1a<a1->GetSize();z1a++)
			if (m_iRefOrSec)
				vec->Add(((CxIntArray*)obs->m_oaAtomOffset[m_oAtoms.m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
					else vec->Add(((CxIntArray*)ref->m_oaAtomOffset[m_oAtoms.m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
	}
	BXOUT;
}


void CRDF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z, z1t, z1a, z2t, z2a;
	CAtomGroup *g1, *g2;
	CxIntArray *a1, *a2;

	vec->RemoveAll_KeepSize();
	for (z=0;z<m_oaVectors.GetSize()/2;z++)
	{
		g1 = (CAtomGroup*)m_oaVectors[z*2];
		for (z1t=0;z1t<g1->m_baAtomType.GetSize();z1t++)
		{
			a1 = (CxIntArray*)g1->m_oaAtoms[z1t];
			for (z1a=0;z1a<a1->GetSize();z1a++)
			{
				g2 = (CAtomGroup*)m_oaVectors[z*2+1];
				for (z2t=0;z2t<g2->m_baAtomType.GetSize();z2t++)
				{
					a2 = (CxIntArray*)g2->m_oaAtoms[z2t];
					for (z2a=0;z2a<a2->GetSize();z2a++)
					{
						if (m_iRefOrSec[0])
							vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
								else vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
						if (m_iRefOrSec[1])
							vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
								else vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
			//			mprintf("Vector z=%d, z1t=%d, z1a=%d, z2t=%d, z2a=%d.\n",z,z1t,z1a,z2t,z2a);
					}
				}
			}
		}
	}
	BXOUT;
}


void CVHDF::BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec)
{
	BXIN;
	int z, z1t, z1a, z2t, z2a;
	CAtomGroup *g1, *g2;
	CxIntArray *a1, *a2;

	vec->RemoveAll_KeepSize();
	for (z=0;z<m_oaVectors.GetSize()/2;z++)
	{
		g1 = (CAtomGroup*)m_oaVectors[z*2];
		for (z1t=0;z1t<g1->m_baAtomType.GetSize();z1t++)
		{
			a1 = (CxIntArray*)g1->m_oaAtoms[z1t];
			for (z1a=0;z1a<a1->GetSize();z1a++)
			{
				g2 = (CAtomGroup*)m_oaVectors[z*2+1];
				for (z2t=0;z2t<g2->m_baAtomType.GetSize();z2t++)
				{
					a2 = (CxIntArray*)g2->m_oaAtoms[z2t];
					for (z2a=0;z2a<a2->GetSize();z2a++)
					{
						if (m_iRefOrSec[0])
							vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
								else vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
						if (m_iRefOrSec[1])
							vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
								else vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
					}
				}
			}
		}
	}
	BXOUT;
}


/*int FindAtom(char *s)
{
	BTIN;
	int z;
	for (z=0;z<g_oaAtoms.GetSize();z++)
		if (mystricmp(s,((CAtom*)g_oaAtoms[z])->m_sName)==0)
		{
			BTOUT; 
			return z;
		}
	BTOUT; 
	return -1;
}*/


int CMolecule::FindAtomInMol(const char *s)
{
	BTIN;
	int z;
	for (z=0;z<m_baAtomIndex.GetSize();z++)
		if (mystricmp(s,((CAtom*)g_oaAtoms[m_baAtomIndex[z]])->m_sName)==0)
		{
			BTOUT; 
			return z;
		}
	BTOUT; 
	return -1;
}


void CMolecule::BuildName()
{
	BTIN;
	int z, z2;
//	char buf[64], buf2[1024];
	CxString buf, buf2;

	if (m_bPseudo)
	{
//		buf2[0] = '$';
//		buf2[1] = 0;
		buf2.sprintf("$");
	} else 
//		buf2[0] = 0;
		buf2.sprintf("");

	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		for (z2=0;z2<m_baAtomIndex.GetSize();z2++)
		{
			if (m_baAtomIndex[z2] != z)
				continue;

			if (m_waAtomCount[z2] > 1)
//				sprintf(buf,"%s%d",((CAtom*)g_oaAtoms[z])->m_sName,m_waAtomCount[z2]);
				buf.sprintf("%s%d",((CAtom*)g_oaAtoms[z])->m_sName,m_waAtomCount[z2]);
			else
//				sprintf(buf,"%s",((CAtom*)g_oaAtoms[z])->m_sName);
				buf.sprintf("%s",((CAtom*)g_oaAtoms[z])->m_sName);

//			strcat(buf2,buf);
			buf2.strcat(buf);
		}
	}

	try { m_sName = new char[strlen(buf2)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(buf2)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,buf2);
	BTOUT; 
}


bool CConditionGroup::Contains(int mol)
{
	BXIN;
	int z;

	if (m_bInactive)
		return m_bAlwaysTrue[mol];

	if (m_bInvertCondition)
	{
		if (m_bAlwaysTrue[mol])
		{
			BXOUT;
			return false;
		}
		for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		{
			if (((CConditionSubGroup*)m_oaConditionSubGroups[z])->Contains(mol))
			{
				BXOUT;
				return false;
			}
		}
		BXOUT;
		return true;
	} else
	{
		if (m_bAlwaysTrue[mol])
		{
			BXOUT;
			return true;
		}
		for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		{
			if (((CConditionSubGroup*)m_oaConditionSubGroups[z])->Contains(mol))
			{
				BXOUT;
				return true;
			}
		}
		BXOUT;
		return false;
	}
}


void CConditionGroup::MarkPassedAtoms(int om, bool passed)
{
	BXIN;
	int z;
//	int ti;

//	mprintf("######## OM %d ##########\n",om+1);

	for (z=0;z<g_iGesVirtAtomCount;z++)
		g_baAtomPassedCondition[z] = 110; // 110 heisst: Atom kommt nicht in einer Condition vor

	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->MarkPassedAtoms(om);

//	for (z=0;z<g_iGesVirtAtomCount;z++)
//		mprintf("\nAtom %d: %d.",z+1,g_baAtomPassedCondition[z]);

//	ti = 0;
	if (m_bInvertCondition)
	{
		for (z=0;z<g_iGesVirtAtomCount;z++)
		{
			if (g_baAtomPassedCondition[z] != 100) // Entweder Condition nicht bestanden oder Atom kam gar nicht darin vor
			{
//				ti++;
				g_baAtomPassedCondition[z] = 1;
//				mprintf("%d: Ok\n",z+1);
//				mprintf("\n %d --> Wird genommen.",z+1);
			} else // Condition bestanden
			{
	//			mprintf("Atom %s%d nicht bestanden (%d).\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,g_waAtomMolNumber[z]+1,g_baAtomPassedCondition[z]);
	//			mprintf("\n%d erfuellt --> durchgefallen",z+1);
				g_baAtomPassedCondition[z] = 0;
			}
		}
	} else
	{
		for (z=0;z<g_iGesVirtAtomCount;z++)
		{
			if (passed)
			{
				if (g_baAtomPassedCondition[z] >= 100) // Entweder Condition bestanden oder Atom kam gar nicht darin vor
				{
//					ti++;
					g_baAtomPassedCondition[z] = 1;
				} else // Condition nicht bestanden
				{
		//			mprintf("Atom %s%d nicht bestanden (%d).\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,g_waAtomMolNumber[z]+1,g_baAtomPassedCondition[z]);
					g_baAtomPassedCondition[z] = 0;
				}
			} else
			{
				if (g_baAtomPassedCondition[z] == 100) // Condition bestanden
				{
//					ti++;
//					mprintf("### Atom %s[%d] %s%d bestanden (%d).\n",((CMolecule*)g_oaMolecules[g_waAtomMolIndex[z]])->m_sName,g_laAtomSMLocalIndex[z]+1,((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,g_waAtomMolNumber[z]+1,g_baAtomPassedCondition[z]);
					g_baAtomPassedCondition[z] = 0;
				} else // Entweder Condition nicht bestanden oder Atom kam gar nicht darin vor
				{
//					mprintf("--- Atom %s[%d] %s%d nicht bestanden (%d).\n",((CMolecule*)g_oaMolecules[g_waAtomMolIndex[z]])->m_sName,g_laAtomSMLocalIndex[z]+1,((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,g_waAtomMolNumber[z]+1,g_baAtomPassedCondition[z]);
					g_baAtomPassedCondition[z] = 1;
				}
			}
		}
	}
//	mprintf("\nMarked %d / %d atoms as enabled.",ti,g_iGesVirtAtomCount);

	BXOUT;
}


void CConditionGroup::ScanNeighborhoodAllOM(CTimeStep *t, CSingleMolecule *rm)
{
	BXIN;
	int z, z2, t1, t2, i;

	if (m_bInactive)
		return;

	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->ScanNeighborhoodAllOM(t,rm);

	if (m_oaConditionSubGroups.GetSize() == 2)
	{
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
		{
			t1 = 0;
			t2 = 0;
			for (z2=0;z2<((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_oaConditions.GetSize();z2++)
				t1 += ((CNbSearch*)((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_oaConditions[z2])->m_iCombPassCount[z];
			for (z2=0;z2<((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_oaConditions.GetSize();z2++)
				t2 += ((CNbSearch*)((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_oaConditions[z2])->m_iCombPassCount[z];
			m_pTable[t1*(((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)+t2]++;
			m_fTableGes++;
		}
	}

	i = 0;
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
	{
		m_fTotal++;
		if (m_bInvertCondition)
		{
			for (z2=0;z2<m_oaConditionSubGroups.GetSize();z2++)
				if (((CConditionSubGroup*)m_oaConditionSubGroups[z2])->Contains(z))
					goto _nopass;
			m_fPassed++;
			m_iPassCounter[z]++;
			m_bAnyPassed = true;
		} else
		{
			for (z2=0;z2<m_oaConditionSubGroups.GetSize();z2++)
			{
				if (((CConditionSubGroup*)m_oaConditionSubGroups[z2])->Contains(z))
				{
//					mprintf("passed[%d]++\n",z);
					m_fPassed++;
					m_iPassCounter[z]++;
					m_bAnyPassed = true;
					i++;
					goto _nopass;
				}
			}
		}
_nopass:;
	}
	m_pHistogram[i]++;
	m_iHistoGes++;

	BXOUT;
}


void CConditionGroup::PreScanNeighborhoodAllOM(CTimeStep *t, CSingleMolecule *rm)
{
	BXIN;
	int z;

	if (m_bInactive)
		return;

	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->PreScanNeighborhoodAllOM(t,rm);

	BXOUT;
}


void CConditionGroup::Parse(int rm, int sm)
{
	BTIN;
	CConditionSubGroup *sg;
	int z;

	m_iShowMol = sm;
	m_iRefMol = rm;

	try { m_bAlwaysTrue = new bool[((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_bAlwaysTrue = NULL; }
	if (m_bAlwaysTrue == NULL) NewException((double)((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_iPassCounter = new long[((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iPassCounter = NULL; }
	if (m_iPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_iOMPassCounter = new long[((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iOMPassCounter = NULL; }
	if (m_iOMPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize();z++)
	{
		m_bAlwaysTrue[z] = false;
		m_iPassCounter[z] = 0;
		m_iOMPassCounter[z] = 0;
	}

	try { m_iRMPassCounter = new long[((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iRMPassCounter = NULL; }
	if (m_iRMPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex.GetSize();z++)
		m_iRMPassCounter[z] = 0;

	mprintf("Different sets of conditions are connected with OR (at least one of them has to be fulfilled).\n");
	mprintf("Different conditions within one set are connected with AND (have to be fulfilled at the same time).\n");
	mprintf("If you create 1 or 2 sets of conditions, a table with subcondition count populations will be printed.\n\n");
_newset:
	mprintf(YELLOW,">>> %d. set of conditions >>>\n",m_oaConditionSubGroups.GetSize()+1);

	try { sg = new CConditionSubGroup(); } catch(...) { sg = NULL; }
	if (sg == NULL) NewException((double)sizeof(CConditionSubGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	sg->m_iShowMol = m_iShowMol;
	sg->m_iNumber = m_oaConditionSubGroups.GetSize();
	m_oaConditionSubGroups.Add(sg);
	sg->Parse(rm,sm);
	mprintf(YELLOW,"\n<<< End of %d. set of conditions <<<\n\n",m_oaConditionSubGroups.GetSize());
	if (AskYesNo("    Enter an additional set of conditions (y/n)? [no] ",false))
		goto _newset;
	
	m_bInvertCondition = AskYesNo("    Invert this condition (only add bin entries if failed) (y/n)? [no] ",false);

	if (m_oaConditionSubGroups.GetSize() == 2)
	{
		m_fTableGes = 0;

		try { m_pTable = new double[(((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_iCombinations+1) * (((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)]; } catch(...) { m_pTable = NULL; }
		if (m_pTable == NULL) NewException((double)(((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_iCombinations+1) * (((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<(((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_iCombinations+1) * (((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1);z++)
			m_pTable[z] = 0;
	}

/*	m_bNeedNbCountMode = false;
	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
	{
		sg = (CConditionSubGroup*)m_oaConditionSubGroups[z];
		for (z2=0;z2<sg->m_oaConditions.GetSize();z2++)
			if (((CNbSearch*)sg->m_oaConditions[z2])->m_iNbCountMin != -1)
				m_bNeedNbCountMode = true;
	}*/
	m_iHistoGes = 0;

	try {m_pHistogram = new unsigned long[((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()+1]; } catch(...) { m_pHistogram = NULL; }
	if (m_pHistogram == NULL) NewException((double)((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()+1*sizeof(unsigned long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()+1;z++)
		m_pHistogram[z] = 0;
	BTOUT;
}


void CConditionGroup::PrintData()
{
	int z;

	mprintf(GREEN,"\n>>> Condition Data >>>\n\n");
	mprintf("    %.0f of %.0f molecules passed the conditions (%.2f percent).\n",m_fPassed,m_fTotal,ZeroDivide(m_fPassed,m_fTotal)*100.0f);
	if (m_fPassed == 0)
	{
		mprintf(YELLOW,"\n    No molecules at all passed the condition.\n");
		mprintf("    This is probably not what you wanted.\n");
	} else
	{
		mprintf(YELLOW,"\nList of reference molecules (%s) that passed the conditions:\n\n",((CMolecule*)g_oaMolecules[m_iRefMol])->m_sName);
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iRefMol])->m_laSingleMolIndex.GetSize();z++)
			if (m_iRMPassCounter[z] > 0)
				mprintf("  - %2d: %10.4f percent of the time (%d hits)\n",z+1,((double)m_iRMPassCounter[z])/g_iSteps*100.0,m_iRMPassCounter[z]);

		mprintf(YELLOW,"\nList of observed molecules (%s) that passed the conditions:\n\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
			if (m_iOMPassCounter[z] > 0)
				mprintf("  - %2d: %10.4f percent of the time (%d hits)\n",z+1,((double)m_iOMPassCounter[z])/g_iSteps*100.0,m_iOMPassCounter[z]);
	}
	mprintf(YELLOW,"\nNeighbor count histogram:\n\n");
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()+1;z++)
		if (m_pHistogram[z] != 0)
			mprintf("  - %d Neighbors: %8.4f percent of the time (%d hits).\n",z,((double)m_pHistogram[z])/m_iHistoGes*100.0,m_pHistogram[z]);
	mprintf(YELLOW,"\nListing for all sets of conditions:\n");
	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->PrintData();
	mprintf(GREEN,"\n<<< End of Condition Data <<<\n\n");
}


void CConditionGroup::PrintData(const char *s)
{
	int z;
	FILE *a;

	mprintf(GREEN,"\n>>> Condition Data >>>\n\n");
	mprintf("    %.0f of %.0f molecules passed the conditions (%.2f percent).\n",m_fPassed,m_fTotal,ZeroDivide(m_fPassed,m_fTotal)*100.0f);
	if (m_fPassed == 0)
	{
		mprintf("\n    No molecules at all passed the condition.\n");
		mprintf("    This is probably not what you wanted.\n");
	}
	mprintf("\nNeighbor count histogram:\n\n");
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()+1;z++)
		if (m_pHistogram[z] != 0)
			mprintf("  - %d Neighbors: %8.4f percent of the time (%d hits).\n",z,((double)m_pHistogram[z])/m_iHistoGes*100.0,m_pHistogram[z]);
	mprintf("\n    Saving detailed condition data as \"%s\"...\n",s);

	a = OpenFileWrite(s,true);

	mfprintf(a,"\n>>> Condition Data >>>\n\n");
	mfprintf(a,"    %.0f of %.0f molecules passed the conditions (%.2f percent).\n",m_fPassed,m_fTotal,ZeroDivide(m_fPassed,m_fTotal)*100.0f);
	if (m_fPassed == 0)
	{
		mfprintf(a,"\n    No molecules at all passed the condition.\n");
		mfprintf(a,"    This is probably not what you wanted.\n");
	} else
	{
		mfprintf(a,"\nList of reference molecules (%s) that passed the conditions:\n\n",((CMolecule*)g_oaMolecules[m_iRefMol])->m_sName);
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iRefMol])->m_laSingleMolIndex.GetSize();z++)
			if (m_iRMPassCounter[z] > 0)
				mfprintf(a,"  - %2d: %10.4f percent of the time (%lu hits)\n",z+1,((double)m_iRMPassCounter[z])/g_iSteps*100.0,m_iRMPassCounter[z]);

		mfprintf(a,"\nList of observed molecules (%s) that passed the conditions:\n\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
			if (m_iOMPassCounter[z] > 0)
				mfprintf(a,"  - %2d: %10.4f percent of the time (%lu hits)\n",z+1,((double)m_iOMPassCounter[z])/g_iSteps*100.0,m_iOMPassCounter[z]);
	}
	mfprintf(a,"\nNeighbor count histogram:\n\n");
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()+1;z++)
		if (m_pHistogram[z] != 0)
			mfprintf(a,"  - %d Neighbors: %8.4f percent of the time (%lu hits).\n",z,((double)m_pHistogram[z])/m_iHistoGes*100.0,m_pHistogram[z]);
	mfprintf(a,"\nListing for all sets of conditions:\n");
	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->PrintData(a);
	mfprintf(a,"\n<<< End of Condition Data <<<\n\n");

	fclose(a);
	mprintf(GREEN,"\n<<< End of Condition Data <<<\n\n");
}


bool CConditionSubGroup::Contains(int mol)
{
	BXIN;
	int z;

	for (z=0;z<m_oaConditions.GetSize();z++)
	{
		if (!((CNbSearch*)m_oaConditions[z])->m_bPassed[mol])
		{
			BXOUT;
			return false;
		}
	}
	BXOUT;
	return true;
}


void CConditionSubGroup::MarkPassedAtoms(int om)
{
	int z, z2;

	for (z=0;z<g_iGesVirtAtomCount;z++)
		m_bTempPassed[z] = false;

	for (z=0;z<m_oaConditions.GetSize();z++)
	{
		for (z2=0;z2<g_iGesVirtAtomCount;z2++)
			if (g_baAtomPassedCondition[z2] < 100)
				g_baAtomPassedCondition[z2] = 0;

		((CNbSearch*)m_oaConditions[z])->MarkPassedAtoms(om);

		for (z2=0;z2<g_iGesVirtAtomCount;z2++)
			if (g_baAtomPassedCondition[z2] == 1)
				m_bTempPassed[z2] = true;
	}

	for (z=0;z<g_iGesVirtAtomCount;z++)
		if (m_bTempPassed[z])
			g_baAtomPassedCondition[z] = 100;
}


void CConditionSubGroup::ScanNeighborhoodAllOM(CTimeStep *t, CSingleMolecule *rm)
{
	int z, z2;

	for (z=0;z<m_oaConditions.GetSize();z++)
		((CNbSearch*)m_oaConditions[z])->ScanAllOM(rm,t);

	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
	{
		m_fTotal++;
		for (z2=0;z2<m_oaConditions.GetSize();z2++)
			if (!((CNbSearch*)m_oaConditions[z2])->m_bPassed[z])
				goto _nopass;
		m_fPassed++;
_nopass:;
	}
}


void CConditionSubGroup::PreScanNeighborhoodAllOM(CTimeStep *t, CSingleMolecule *rm)
{
	int z;

	for (z=0;z<m_oaConditions.GetSize();z++)
		((CNbSearch*)m_oaConditions[z])->PreScanAllOM(rm,t);
}


void CConditionSubGroup::Parse(int rm, int sm)
{
	BTIN;
	CNbSearch *n;

	m_iCombinations = 0;
	do
	{
		try { n = new CNbSearch(); } catch(...) { n = NULL; }
		if (n == NULL) NewException((double)sizeof(CNbSearch),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		n->m_iNumber = m_oaConditions.GetSize();
		m_oaConditions.Add(n);
		n->Parse(rm,sm,false);
		m_iCombinations += n->m_iCombinationsEnabled;
	} while (AskYesNo("    Add another condition to this set of conditions (y/n)? [no] ",false));

	try { m_bTempPassed = new bool[g_iGesVirtAtomCount]; } catch(...) { m_bTempPassed = NULL; }
	if (m_bTempPassed == NULL) NewException((double)g_iGesVirtAtomCount*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	BTOUT;
}


void CADF::ParseCondition(int rm, bool nocrit)
{
	BTIN;
//	char buf[1024];
	CxString buf;
	int z, z2, ti;
	CAtomGroup *ag;

	mprintf(WHITE,"\n>>> Angular Condition >>>\n");
	for (z=0;z<2;z++)
		m_iVecType[z] = 0;
	z2 = 0;
	m_iCombinations = 0;
	do {
		ti = 1;
		if (z2 != 0)
			mprintf("\n    %d. set of vectors\n\n",z2+1);
		for (z=0;z<2;z++)
		{
/*			if (z == 1)
				m_bSameFoot = AskYesNo("    Should the base points of the 1st and 2nd vectors always be equal (y/n)? [no] ",false);
					else m_bSameFoot = false;*/
			m_bOrtho[z] = (AskRangeInteger("\n    Should the %d. vector connect two points (0) or stand perpendicular on 3 points (1)? [0] ",0,1,0,z+1) != 0);
			mprintf("\n");
			if (m_bOrtho[z])
			{
/*				if (m_bSameFoot && (z == 1))
				{
					m_iRefOrSec[1][0] = m_iRefOrSec[0][0];
					m_oaVectors.Add(NULL);
				} else*/
				{
_ax1:				if (m_iShowMol != -1)
						m_iRefOrSec[z][0] = AskRangeInteger("    Take atom(s) in base point from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
							else m_iRefOrSec[z][0] = 0;
					mprintf("      Enter atom(s) in the base point (e.g. C7): ");
					inpprintf("! Enter atom(s) in the base point (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						delete ag;
						goto _ax1;
					}
					m_oaVectors.Add(ag);
					ti *= ag->m_iAtomGes;
				}
_ax2:			if (m_iShowMol != -1)
					m_iRefOrSec[z][1] = AskRangeInteger("   Take 2nd atom(s) of normal plane from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else m_iRefOrSec[z][1] = 0;
				mprintf("      Enter 2nd atom(s) of normal plane (e.g. C7): ");
				inpprintf("! Enter 2nd atom(s) of normal plane (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (!ag->ParseAtoms((!m_iRefOrSec[z][1])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					delete ag;
					goto _ax2;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
_ax3:			if (m_iShowMol != -1)
					m_iRefOrSec[z][2] = AskRangeInteger("    Take 3rd atom(s) of normal plane from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else m_iRefOrSec[z][2] = 0;
				mprintf("      Enter 3rd atom(s) of normal plane (e.g. C7): ");
				inpprintf("! Enter 3rd atom(s) of normal plane (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (!ag->ParseAtoms((!m_iRefOrSec[z][2])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					delete ag;
					goto _ax3;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
			} else // IF NOT ORTHO
			{
/*				if (m_bSameFoot && (z == 1))
				{
					m_iRefOrSec[1][0] = m_iRefOrSec[0][0];
					m_oaVectors.Add(NULL);
				} else*/
				{
_ax4:				if (m_iShowMol != -1)
						m_iRefOrSec[z][0] = AskRangeInteger("    Take atom(s) in base point from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
							else m_iRefOrSec[z][0] = 0;
					mprintf("      Enter atom(s) in the base point (e.g. C7): ");
					inpprintf("! Enter atom(s) in the base point (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						delete ag;
						goto _ax4;
					}
					m_oaVectors.Add(ag);
					ti *= ag->m_iAtomGes;
				} // END IF NOT SAMEFOOT
_ax5:			if (m_iShowMol != -1)
					m_iRefOrSec[z][1] = AskRangeInteger("    Take atom(s) in tip point from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else m_iRefOrSec[z][1] = 0;
				mprintf("      Enter atom(s) in the tip point (e.g. C7): ");
				inpprintf("! Enter atom(s) in the tip point (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if (!ag->ParseAtoms((!m_iRefOrSec[z][1])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					delete ag;
					goto _ax5;
				}
				m_oaVectors.Add(ag);
				ti *= ag->m_iAtomGes;
				m_oaVectors.Add(NULL);
			} // END IF NOT ORTHO
		} // END FOR 0..1
		m_iCombinations += ti;
		z2++;
	} while (AskYesNo("\n    Enter another set of vectors (y/n)? [no] ",false));
	if (!nocrit)
	{
		do {
			m_faMinMaxAngle.Add(AskFloat("    Enter minimal angle between the vectors in degree: [0.0] ",0.0f));
			m_faMinMaxAngle.Add(AskFloat("    Enter maximal angle between the vectors in degree: [180.0] ",180.0f));
		} while (AskYesNo("    Enter another angle interval (y/n)? [no] ",false));
	} else
	{
		m_faMinMaxAngle.Add(0.0f);
		m_faMinMaxAngle.Add(180.0f);
	}
	mprintf(WHITE,"\n<<< End of Angular Condition <<<\n\n");
	BTOUT;
}


void CADF::ParseConditionGrid(int rm, int gridmode)
{
	BTIN;
//	char buf[1024];
	CxString buf;
	int z, z2, ti;
	CAtomGroup *ag;

	mprintf(WHITE,"\n>>> Angular Condition >>>\n");
	for (z=0;z<2;z++)
		m_iVecType[z] = 0;
	z2 = 0;
	m_iCombinations = 0;
	do {
		ti = 1;
		if (z2 != 0)
			mprintf("\n    %d. set of vectors\n\n",z2+1);
		for (z=0;z<2;z++)
		{
/*			if (z == 1)
				m_bSameFoot = AskYesNo("    Should the base points of the 1st and 2nd vectors always be equal (y/n)? [no] ",false);
					else m_bSameFoot = false;*/
			m_bOrtho[z] = (AskRangeInteger("\n    Should the %d. vector connect two points (0) or stand perpendicular on 3 points (1)? [0] ",0,1,0,z+1) != 0);
			mprintf("\n");
			if (m_bOrtho[z])
			{
/*				if (m_bSameFoot && (z == 1))
				{
					m_iRefOrSec[1][0] = m_iRefOrSec[0][0];
					m_oaVectors.Add(NULL);
				} else*/
				{
_ax1:				if (m_iShowMol != -1)
						m_iRefOrSec[z][0] = AskRangeInteger("    Take atom(s) in base point from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
							else m_iRefOrSec[z][0] = 0;
					mprintf("      Enter atom(s) in the base point (e.g. C7): ");
					inpprintf("! Enter atom(s) in the base point (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					m_oaVectors.Add(ag);
					if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _ax1;
					}
					ti *= ag->m_iAtomGes;
				}
_ax2:			if (m_iShowMol != -1)
					m_iRefOrSec[z][1] = AskRangeInteger("   Take 2nd atom(s) of normal plane from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else m_iRefOrSec[z][1] = 0;
				mprintf("      Enter 2nd atom(s) of normal plane (e.g. C7): ");
				inpprintf("! Enter 2nd atom(s) of normal plane (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				m_oaVectors.Add(ag);
				if (!ag->ParseAtoms((!m_iRefOrSec[z][1])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _ax2;
				}
				ti *= ag->m_iAtomGes;
_ax3:			if (m_iShowMol != -1)
					m_iRefOrSec[z][2] = AskRangeInteger("    Take 3rd atom(s) of normal plane from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else m_iRefOrSec[z][2] = 0;
				mprintf("      Enter 3rd atom(s) of normal plane (e.g. C7): ");
				inpprintf("! Enter 3rd atom(s) of normal plane (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				m_oaVectors.Add(ag);
				if (!ag->ParseAtoms((!m_iRefOrSec[z][2])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _ax3;
				}
				ti *= ag->m_iAtomGes;
			} else // IF NOT ORTHO
			{
/*				if (m_bSameFoot && (z == 1))
				{
					m_iRefOrSec[1][0] = m_iRefOrSec[0][0];
					m_oaVectors.Add(NULL);
				} else*/
				{
_ax4:				if (m_iShowMol != -1)
						m_iRefOrSec[z][0] = AskRangeInteger("    Take atom(s) in base point from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
							else m_iRefOrSec[z][0] = 0;
					mprintf("      Enter atom(s) in the base point (e.g. C7): ");
					inpprintf("! Enter atom(s) in the base point (e.g. C7):\n");
					myget(&buf);

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					m_oaVectors.Add(ag);
					if (!ag->ParseAtoms((!m_iRefOrSec[z][0])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
					{
						eprintf("Wrong input.\n");
						inpprintf("! Wrong input.\n");
						goto _ax4;
					}
					ti *= ag->m_iAtomGes;
				} // END IF NOT SAMEFOOT
_ax5:			if (m_iShowMol != -1)
					m_iRefOrSec[z][1] = AskRangeInteger("    Take atom(s) in tip point from 1st mol. %s (0) or from 2nd mol. %s (1)? [%d] ",0,1,z,((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,z);
						else m_iRefOrSec[z][1] = 0;
				mprintf("      Enter atom(s) in the tip point (e.g. C7): ");
				inpprintf("! Enter atom(s) in the tip point (e.g. C7):\n");
				myget(&buf);

				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				m_oaVectors.Add(ag);
				if (!ag->ParseAtoms((!m_iRefOrSec[z][1])?(CMolecule*)g_oaMolecules[rm]:(CMolecule*)g_oaMolecules[m_iShowMol],buf))
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _ax5;
				}
				ti *= ag->m_iAtomGes;
				m_oaVectors.Add(NULL);
			} // END IF NOT ORTHO
		} // END FOR 0..1
		m_iCombinations += ti;
		z2++;
	} while (AskYesNo("\n    Enter another set of vectors (y/n)? [no] ",false));
	if ((gridmode == 1) || (gridmode == 6))
	{
		m_faMinMaxAngle.Add(AskFloat("    Enter minimal angle between the vectors in degree: [0.0] ",0.0f));
		m_faMinMaxAngle.Add(AskFloat("    Enter maximal angle between the vectors in degree: [180.0] ",180.0f));
	} else
	{
		m_faMinMaxAngle.Add(0.0f);
		m_faMinMaxAngle.Add(180.0f);
	}
	mprintf(WHITE,"\n<<< End of Angular Condition <<<\n\n");
	BTOUT;
}


void CADF::ParseCondition_OnlyValues()
{
	BTIN;
	int z;

	mprintf("    Angular condition between vectors:\n");
	for (z=0;z<m_oaVectors.GetSize()/6;z++)
	{
		mprintf("      ");
		if (m_oaVectors[z*6+2] != NULL)
			mprintf("normal(%s, %s, %s)",((CAtomGroup*)m_oaVectors[z*6])->m_sName,((CAtomGroup*)m_oaVectors[z*6+1])->m_sName,((CAtomGroup*)m_oaVectors[z*6+2])->m_sName);
				else mprintf("(%s --> %s)",((CAtomGroup*)m_oaVectors[z*6])->m_sName,((CAtomGroup*)m_oaVectors[z*6+1])->m_sName);
		mprintf(" and ");
		if (m_oaVectors[z*6+5] != NULL)
			mprintf("normal(%s, %s, %s)",((CAtomGroup*)m_oaVectors[z*6+3])->m_sName,((CAtomGroup*)m_oaVectors[z*6+4])->m_sName,((CAtomGroup*)m_oaVectors[z*6+5])->m_sName);
				else mprintf("(%s --> %s)",((CAtomGroup*)m_oaVectors[z*6+3])->m_sName,((CAtomGroup*)m_oaVectors[z*6+4])->m_sName);
		mprintf("\n");
	}
	for (z=0;z<m_faMinMaxAngle.GetSize()/2;z++)
	{
		mprintf("    Interval %d:\n",z+1,m_faMinMaxAngle.GetSize()/2);
		m_faMinMaxAngle[z*2] = AskFloat("      Enter minimal angle between the vectors in degree: [0.0] ",0.0f);
		m_faMinMaxAngle[z*2+1] = AskFloat("      Enter maximal angle between the vectors in degree: [180.0] ",180.0f);
	} 
//	mprintf(WHITE,"\n<<< End of Angular Condition <<<\n\n");
	BTOUT;
}


CObservation::CObservation()
{
	m_pConditions = NULL;
	m_pConditionsOM2 = NULL;
	m_pMSD = NULL;
	m_pSDF = NULL;
	m_pCDF = NULL;
	m_pVACF = NULL;
	m_pDipACF = NULL;
	m_pADF = NULL;
	m_pDDF = NULL;
	m_pRDF = NULL;
	m_pDipDF = NULL;
	m_pVDF = NULL;
	m_pPlProj = NULL;
	m_bBinOnlyPassedAtoms = false;
	m_bBinOnlyNotPassedAtoms = false;

	m_waSaveRefList.SetName("CObservation::m_waSaveRefList");
	m_waSaveShowList.SetName("CObservation::m_waSaveShowList");
	m_waObsRefList.SetName("CObservation::m_waObsRefList");
	m_waObsShowList.SetName("CObservation::m_waObsShowList");
	m_waObsShow2List.SetName("CObservation::m_waObsShow2List");
	m_waDecompTypeRefOffs.SetName("CObservation::m_waDecompTypeRefOffs");
	m_waDecompTypeObsOffs.SetName("CObservation::m_waDecompTypeObsOffs");
	m_waDecompTypeRefList.SetName("CObservation::m_waDecompTypeRefList");
	m_waDecompTypeObsList.SetName("CObservation::m_waDecompTypeObsList");
	m_iaRMRegions.SetName("CObservation::m_iaRMRegions");
	m_iaOM1Regions.SetName("CObservation::m_iaOM1Regions");
	m_iaOM2Regions.SetName("CObservation::m_iaOM2Regions");
}


CObservation::~CObservation()
{
	int z;

	if (m_pConditions != NULL)
	{
		delete m_pConditions;
		m_pConditions = NULL;
	}
	if (m_pMSD != NULL)
	{
		delete m_pMSD;
		m_pMSD = NULL;
	}
	if (m_pSDF != NULL)
	{
		delete m_pSDF;
		m_pSDF = NULL;
	}
	if (m_pCDF != NULL)
	{
		delete m_pCDF;
		m_pCDF = NULL;
	}
	if (m_pVACF != NULL)
	{
		delete m_pVACF;
		m_pVACF = NULL;
	}
	if (m_pDipACF != NULL)
	{
		delete m_pDipACF;
		m_pDipACF = NULL;
	}

	if (m_pRDF != NULL)
	{
		for (z=0;z<g_iCDFChannels;z++)
		{
			if (m_pRDF[z] != NULL)
			{
				delete m_pRDF[z];
				m_pRDF[z] = NULL;
			}
		}
		delete[] m_pRDF;
		m_pRDF = NULL;
	}
	if (m_pADF != NULL)
	{
		for (z=0;z<g_iCDFChannels;z++)
		{
			if (m_pADF[z] != NULL)
			{
				delete m_pADF[z];
				m_pADF[z] = NULL;
			}
		}
		delete[] m_pADF;
		m_pADF = NULL;
	}
	if (m_pDDF != NULL)
	{
		for (z=0;z<g_iCDFChannels;z++)
		{
			if (m_pDDF[z] != NULL)
			{
				delete m_pDDF[z];
				m_pDDF[z] = NULL;
			}
		}
		delete[] m_pDDF;
		m_pDDF = NULL;
	}
	if (m_pDipDF != NULL)
	{
		for (z=0;z<g_iCDFChannels;z++)
		{
			if (m_pDipDF[z] != NULL)
			{
				delete m_pDipDF[z];
				m_pDipDF[z] = NULL;
			}
		}
		delete[] m_pDipDF;
		m_pDipDF = NULL;
	}
	if (m_pVDF != NULL)
	{
		for (z=0;z<g_iCDFChannels;z++)
		{
			if (m_pVDF[z] != NULL)
			{
				delete m_pVDF[z];
				m_pVDF[z] = NULL;
			}
		}
		delete[] m_pVDF;
		m_pVDF = NULL;
	}
}


void CSingleMolecule::BuildAtomCodes()
{
	int z, z2, c1, c2, i;
	double ac;
	CMolAtom *ma;

	// Die Anfangswerte der AtomCodes: [Ordnungszahl] * 10.0 + [Zahl der Nicht-Wasserstoff-Bindungen]
	for (z=0;z<m_oaMolAtoms.GetSize();z++)
	{
		ma = (CMolAtom*)m_oaMolAtoms[z];
		ma->m_fAtomCode = 10.0 * ((CAtom*)g_oaAtoms[m_baAtomIndex[ma->m_iType]])->m_pElement->m_fMass;

		i = 0;
		for (z2=0;z2<ma->m_oaBonds.GetSize();z2++)
		{
			// Alle Wasserstoff-Atome ueberspringen
			if (((CAtom*)g_oaAtoms[m_baAtomIndex[((CMolAtom*)ma->m_oaBonds[z2])->m_iType]])->m_pElement->m_fMass < 4.5)
				continue;
			ma->m_fAtomCode++;
			i++;
		}

		if (g_bVerbose)
			mprintf("    Atom %d (%s%d): 10.0 * %.2f + %d = %.2f\n",z+1,((CAtom*)g_oaAtoms[m_baAtomIndex[ma->m_iType]])->m_sName,ma->m_iOffset+1,((CAtom*)g_oaAtoms[m_baAtomIndex[ma->m_iType]])->m_pElement->m_fMass,i,ma->m_fAtomCode);
	}

	i = 0;
	do {
		c1 = CountDifferentAtomCodes();

		if (g_bVerbose)
			mprintf(WHITE,"\n  Cycle %d: %d different atom codes exist.\n\n",i+1,c1);

		for (z=0;z<m_oaMolAtoms.GetSize();z++)
		{
			ma = (CMolAtom*)m_oaMolAtoms[z];
			ma->m_fTempAtomCode = ma->m_fAtomCode * 5.0;

			for (z2=0;z2<ma->m_oaBonds.GetSize();z2++)
				ma->m_fTempAtomCode += ((CMolAtom*)ma->m_oaBonds[z2])->m_fAtomCode;
		}
		for (z=0;z<m_oaMolAtoms.GetSize();z++)
		{
			ma = (CMolAtom*)m_oaMolAtoms[z];
			ma->m_fAtomCode = ma->m_fTempAtomCode;

			if (g_bVerbose)
				mprintf("    Atom %d (%s%d): %.2f\n",z+1,((CAtom*)g_oaAtoms[m_baAtomIndex[ma->m_iType]])->m_sName,ma->m_iOffset+1,ma->m_fAtomCode);
		}
		c2 = CountDifferentAtomCodes();
		i++;
//		mprintf("Iteration %d: %d classes before, %d classes after.\n",i,c1,c2);
	} while (c1 != c2);

	if (g_bVerbose)
		mprintf(WHITE,"\nSorting...\n");

	m_iAtomClasses = c2;
//	mprintf("Finished.\n");
//	mprintf("%d Iterations, %d atom classes found.\n",i,c2);
//	mprintf("Sorting Atom Codes...");

	// Sortieren mittels StackSort
	for (z=0;z<m_oaMolAtoms.GetSize();z++)
	{
		ac = -1;
		i = -1;
		for (z2=z;z2<m_oaMolAtoms.GetSize();z2++)
		{
			if (((CMolAtom*)m_oaMolAtoms[z2])->m_fAtomCode > ac)
			{
				ac = ((CMolAtom*)m_oaMolAtoms[z2])->m_fAtomCode;
				i = z2;
			}
		}
		if (i != -1)
		{
			ma = (CMolAtom*)m_oaMolAtoms[z];
			m_oaMolAtoms[z] = m_oaMolAtoms[i];
			m_oaMolAtoms[i] = ma;
			((CMolAtom*)m_oaMolAtoms[z])->m_iMolAtomNumber = z;
		} else
		{
			eprintf("CSingleMolecule::BuildAtomCodes(): Weird error.\n");
			return;
		}
	}

	if (g_bVerbose)
	{
		for (z=0;z<m_oaMolAtoms.GetSize();z++)
		{
			ma = (CMolAtom*)m_oaMolAtoms[z];
			mprintf("    Atom %d (%s%d): %.2f\n",z+1,((CAtom*)g_oaAtoms[m_baAtomIndex[ma->m_iType]])->m_sName,ma->m_iOffset+1,ma->m_fAtomCode);
		}
		mprintf(WHITE,"Finished.\n\n");
	}
//	mprintf("Finished.\n");
}


int CSingleMolecule::CountDifferentAtomCodes()
{
	int z, z2, i;
	double *d;

	try { d = new double[m_oaMolAtoms.GetSize()]; } catch(...) { d = NULL; }
	if (d == NULL) NewException((double)m_oaMolAtoms.GetSize()*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	i = 0;
	for (z=0;z<m_oaMolAtoms.GetSize();z++)
	{
		for (z2=0;z2<i;z2++)
			if (d[z2] == ((CMolAtom*)m_oaMolAtoms[z])->m_fAtomCode)
				goto _next;
		d[i] = ((CMolAtom*)m_oaMolAtoms[z])->m_fAtomCode;
		i++;
_next:;
	}
	delete[] d;
	return i;
}


void CMSD::Parse()
{
//	char buf[1024];
	CxString buf;
	int ti;
//	float tf;

	BTIN;

	try { m_pMSD = new CAF(); } catch(...) { m_pMSD = NULL; }
	if (m_pMSD == NULL) NewException((double)sizeof(CAF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf(WHITE,"\n*** Mean Square Displacement\n\n");

	try { m_pAtomGroup = new CAtomGroup(); } catch(...) { m_pAtomGroup = NULL; }
	if (m_pAtomGroup == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
_rdfatom1:
	mprintf("    Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	inpprintf("! Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2]\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	myget(&buf);
	if (strlen(buf) == 0)
	{
		if (!m_pAtomGroup->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],"#2"))
		{
			eprintf("Weird error.\n");
			inpprintf("! Weird error.\n");
			goto _rdfatom1;
		}
	} else if (!m_pAtomGroup->ParseAtoms((CMolecule*)g_oaMolecules[m_iShowMol],buf))
	{
		eprintf("Wrong input.\n");
		inpprintf("! Wrong input.\n");
		goto _rdfatom1;
	}

//_depth:
	if (g_iTrajSteps != -1)
		m_iResolution = AskUnsignedInteger("    Enter the resolution (=depth) for this MSD (in time steps): [%d] ",int(g_iTrajSteps*0.75),int(g_iTrajSteps*0.75));
			else m_iResolution = AskUnsignedInteger("    Enter the resolution (=depth) for this MSD (in time steps): [10000] ",10000);

	m_iShowAtoms = m_pAtomGroup->m_iAtomGes;

/*	if (g_iTrajSteps != -1)
	{
		if (g_bMSDCacheMode)
		{
			tf = g_iTrajSteps*m_pAtomGroup->m_iAtomGes*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*3.0f*sizeof(float)/1024.0f/1024.0f;
			if (tf >= 10.0f)
				if (!AskYesNo("    This will occupy %.0f MB RAM (once, not for every MSD!). Continue (y/n)? [yes] ",true,tf))
					goto _depth;
		} else
		{
			tf = m_iResolution*g_iGesVirtAtomCount*12.0f*sizeof(double)/1024.0f/1024.0f;
			if (tf >= 10.0f)
				if (!AskYesNo("    This will occupy %.0f MB RAM (once, not for every MSD!). Continue (y/n)? [yes] ",true,tf))
					goto _depth;
		}
	}*/
/*	if (m_iResolution > g_iMaxMSDDepth)
		g_iMaxMSDDepth = m_iReolution;*/

	ti = m_iResolution / 1000;
	if (ti < 1)
		ti = 1;
	m_iStride = AskUnsignedInteger("    Take every n-th step for the time axis of the MSD: [%d] ",ti,ti);
	m_iStride2 = AskUnsignedInteger("    Shift correlation window n time steps at once: [%d] ",m_iStride,m_iStride);

	if (g_bAdvanced2)
	{
		m_bSplit = AskYesNo("    Decompose this MSD into contributions from each individual molecule (y/n)? [no] ",false);
		if (AskYesNo("    Take into account only certain dimensions for displacement (y/n)? [no] ",false))
		{
			m_bTakeX = AskYesNo("      Take into account contributions along X axis (y/n)? [yes] ",true);
			m_bTakeY = AskYesNo("      Take into account contributions along Y axis (y/n)? [yes] ",true);
			m_bTakeZ = AskYesNo("      Take into account contributions along Z axis (y/n)? [yes] ",true);
		} else
		{
			m_bTakeX = true;
			m_bTakeY = true;
			m_bTakeZ = true;
		}
	} else
	{
		m_bSplit = false;
		m_bTakeX = true;
		m_bTakeY = true;
		m_bTakeZ = true;
	}

	BuildName();
	BTOUT;
}


void CMSD::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;

//	sprintf(tmp,"%s_%s",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_pAtomGroup->m_sName);
	tmp.sprintf("%s_%s",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,m_pAtomGroup->m_sName);

	if ((!m_bTakeX) || (!m_bTakeY) || (!m_bTakeZ))
	{
//		strcat(tmp,"_");
		tmp.strcat("_");

		if (m_bTakeX)
//			strcat(tmp,"X");
			tmp.strcat("X");

		if (m_bTakeY)
//			strcat(tmp,"Y");
			tmp.strcat("Y");

		if (m_bTakeZ)
//			strcat(tmp,"Z");
			tmp.strcat("Z");
	}

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CVHDF::CorrectCount()
{
	int x, y;

	for (y=0;y<m_iResolution;y++)
		for (x=0;x<m_iDepth/m_iStride;x++)
			if (m_pCount[x] != 0)
				m_pVHDF->m_pBin[y*m_iDepth/m_iStride+x] /= m_pCount[x];
}


void CConditionSubGroup::PrintData()
{
	int z;
	double d;

	mprintf(GREEN,"\n*** Data for %d. set of conditions ***\n\n",m_iNumber+1);
	mprintf("    %.0f of %.0f molecules passed this set of conditions (%.4f percent).\n",m_fPassed,m_fTotal,ZeroDivide(m_fPassed,m_fTotal)*100.0f);
	d = 1.0;
	for (z=0;z<m_oaConditions.GetSize();z++)
	{
		mprintf("      - Condition %d: %.0f of %.0f molecules passed (%.4f percent).\n",z+1,((CNbSearch*)m_oaConditions[z])->m_fMoleculesPassed,((CNbSearch*)m_oaConditions[z])->m_fMoleculesTotal,ZeroDivide(((CNbSearch*)m_oaConditions[z])->m_fMoleculesPassed,((CNbSearch*)m_oaConditions[z])->m_fMoleculesTotal)*100.0);
		d *= ZeroDivide(((CNbSearch*)m_oaConditions[z])->m_fMoleculesPassed,((CNbSearch*)m_oaConditions[z])->m_fMoleculesTotal);
	}
	mprintf("\n");
	if (m_oaConditions.GetSize() > 1)
	{
		mprintf("    The product of the condition probabilities is %.4f percent.\n",d*100.0);
		mprintf("    If the conditions would be uncorrelated, this would be also the probability for the set.\n");
		mprintf("    Therefore, the conditions are %s correlated for %.4f percent.\n\n",(d<ZeroDivide(m_fPassed,m_fTotal))?"positively":"negatively",(d<ZeroDivide(m_fPassed,m_fTotal))?(ZeroDivide(ZeroDivide(m_fPassed,m_fTotal),d)-1.0)*100.0:(ZeroDivide(d,ZeroDivide(m_fPassed,m_fTotal))-1.0)*100.0);
	}

	for (z=0;z<m_oaConditions.GetSize();z++)
		((CNbSearch*)m_oaConditions[z])->PrintTable();
}


void CConditionSubGroup::PrintData(FILE *a)
{
	int z;

	mfprintf(a,"\n*** Data for %d. set of conditions ***\n\n",m_iNumber+1);
	mfprintf(a,"    %.0f of %.0f molecules passed this set of conditions (%.4f percent).\n",m_fPassed,m_fTotal,ZeroDivide(m_fPassed,m_fTotal)*100.0f);
	for (z=0;z<m_oaConditions.GetSize();z++)
		((CNbSearch*)m_oaConditions[z])->PrintTable(a);
}


void CConditionSubGroup::PrintSingle(int om)
{
	int z;

	for (z=0;z<m_oaConditions.GetSize();z++)
		((CNbSearch*)m_oaConditions[z])->PrintSingle(om);
}


void CConditionGroup::PrintSingle(int om)
{
	int z;

	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->PrintSingle(om);
}


void CConditionGroup::PrintTable()
{
	int z, z2;
	double tf1, tf2, tf3;

	mprintf(WHITE,"*** Condition Table ***\n\n");
	mprintf("    The rows indicate how many subconditions of condition 1 are fulfilled.\n");
	mprintf("    The columns indicate how many subconditions of condition 2 are fulfilled.\n\n");

	mprintf("     ");
	for (z=0;z<((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1;z++)
		mprintf("| %4d    ",z);
	mprintf("\n");
	mprintf("-----");
	for (z=0;z<((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1;z++)
		mprintf("|---------",z);
	mprintf("\n");
	tf1 = 0;
	tf2 = 0;
	tf3 = 0;
	for (z=0;z<((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_iCombinations+1;z++)
	{
		mprintf(" %3d ",z);
		for (z2=0;z2<((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1;z2++)
		{
			if ((z > 0) && (z2 == 0))
				tf1 += m_pTable[z*(((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)+z2];
			if ((z == 0) && (z2 > 0))
				tf2 += m_pTable[z*(((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)+z2];
			if ((z > 0) && (z2 > 0))
				tf3 += m_pTable[z*(((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)+z2];
			mprintf("| %7.3f ",m_pTable[z*(((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)+z2]/m_fTableGes*100.0);
		}
		mprintf("\n");
	}
	mprintf("\n");
	mprintf("    A RM/OM pair fulfills no condition:      %7.3f percent of the time.\n",m_pTable[0]/m_fTableGes*100.0);
	mprintf("    A RM/OM pair fulfills only condition 1:  %7.3f percent of the time.\n",tf1/m_fTableGes*100.0);
	mprintf("    A RM/OM pair fulfills only condition 2:  %7.3f percent of the time.\n",tf2/m_fTableGes*100.0);
	mprintf("    A RM/OM pair fulfills both conditions:   %7.3f percent of the time.\n",tf3/m_fTableGes*100.0);
	mprintf("                                    Total:   %7.3f percent of the time.\n",100.0);
	mprintf("\n");
	mprintf(WHITE,"*** Condition Table End ***\n\n");
}


void CObservation::ListCDFObservations(int z)
{
	int z2, ti, ti2, ti3;
	CxIntArray tempwa;

	switch(g_iObsChannel[z])
	{
		case 0: // RDF
			if (m_bOthers)
			{
				if (m_bSecondShowMol && (z == 1))
					m_pRDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol2])->m_laSingleMolIndex[0]],&tempwa);
				else
					m_pRDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]],&tempwa);
			} else
				m_pRDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],&tempwa);
			for (z2=0;z2<tempwa.GetSize()/2;z2++)
			{
				ti = tempwa[z2*2];
				ti2 = tempwa[z2*2+1];
				mprintf("  * %2d.) Distance %s%d (%s) <--> %s%d (%s)\n",z2+1,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
			}
			break;

		case 1: // ADF
			if (m_bOthers)
			{
				if (m_bSecondShowMol && (z == 1))
					m_pADF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol2])->m_laSingleMolIndex[0]],&tempwa);
				else
					m_pADF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]],&tempwa);
			} else
				m_pADF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],&tempwa);
			for (z2=0;z2<tempwa.GetSize()/6;z2++)
			{
				mprintf("  * %2d.) Angle ",z2+1);
				if (m_pADF[z]->m_bOrtho[0])
				{
					ti = tempwa[z2*6];
					ti2 = tempwa[z2*6+1];
					ti3 = tempwa[z2*6+2];
					mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)] to ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
				} else
				{
					ti = tempwa[z2*6];
					ti2 = tempwa[z2*6+1];
					mprintf("[%s%d (%s) --> %s%d (%s)] to ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
				}
				if (m_pADF[z]->m_bOrtho[1])
				{
					ti = tempwa[z2*6+3];
					ti2 = tempwa[z2*6+4];
					ti3 = tempwa[z2*6+5];
					mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)]\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
				} else
				{
					ti = tempwa[z2*6+3];
					ti2 = tempwa[z2*6+4];
					mprintf("[%s%d (%s) --> %s%d (%s)]\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
				}
			}
			break;

		case 2: // DDF
			if (m_bOthers)
			{
				if (m_bSecondShowMol && (z == 1))
					m_pDDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol2])->m_laSingleMolIndex[0]],&tempwa);
				else
					m_pDDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]],&tempwa);
			} else
				m_pDDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],&tempwa);
			for (z2=0;z2<tempwa.GetSize()/9;z2++)
			{
				mprintf("  * %2d.) Dihedral Angle ",z2+1);
				if (m_pDDF[z]->m_bOrtho[0])
				{
					ti = tempwa[z2*9];
					ti2 = tempwa[z2*9+1];
					ti3 = tempwa[z2*9+2];
					mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)] - ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
				} else
				{
					ti = tempwa[z2*9];
					ti2 = tempwa[z2*9+1];
					mprintf("[%s%d (%s) --> %s%d (%s)] - ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
				}
				if (m_pDDF[z]->m_bOrtho[1])
				{
					ti = tempwa[z2*9+3];
					ti2 = tempwa[z2*9+4];
					ti3 = tempwa[z2*9+5];
					mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)] - ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
				} else
				{
					ti = tempwa[z2*9+3];
					ti2 = tempwa[z2*9+4];
					mprintf("[%s%d (%s) --> %s%d (%s)] - ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
				}
				if (m_pDDF[z]->m_bOrtho[2])
				{
					ti = tempwa[z2*9+6];
					ti2 = tempwa[z2*9+7];
					ti3 = tempwa[z2*9+8];
					mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)]\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
				} else
				{
					ti = tempwa[z2*9+6];
					ti2 = tempwa[z2*9+7];
					mprintf("[%s%d (%s) --> %s%d (%s)]\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
				}
			}
			break;

		case 3:
			eprintf("(not implemented)\n");
			break;

		case 4:
			eprintf("(not implemented)\n");
			break;

		case 5: // PlDF
			if (m_bOthers)
			{
				if (m_bSecondShowMol && (z == 1))
					m_pPlDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol2])->m_laSingleMolIndex[0]],&tempwa);
				else
					m_pPlDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]],&tempwa);
			} else
				m_pPlDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],&tempwa);
			for (z2=0;z2<tempwa.GetSize()/4;z2++)
			{
				mprintf("  * %2d.) Distance from Plane ",z2+1);
				if (m_pPlDF[z]->m_bNormal)
				{
					ti = tempwa[z2*4];
					ti2 = tempwa[z2*4+1];
					mprintf("[%s%d (%s) --> %s%d (%s)] - ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
				} else
				{
					ti = tempwa[z2*4];
					ti2 = tempwa[z2*4+1];
					ti3 = tempwa[z2*4+2];
					mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)] - ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
				}
				ti = tempwa[z2*4+3];
				mprintf("%s%d (%s)\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName);
			}
			break;

		case 6: // LiDF
			if (m_bOthers)
			{
				if (m_bSecondShowMol && (z == 1))
					m_pLiDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol2])->m_laSingleMolIndex[0]],&tempwa);
				else
					m_pLiDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[0]],&tempwa);
			} else
				m_pLiDF[z]->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]],&tempwa);
			for (z2=0;z2<tempwa.GetSize()/4;z2++)
			{
				mprintf("  * %2d.) Distance from Plane ",z2+1);
				if (m_pLiDF[z]->m_bNormal)
				{
					ti = tempwa[z2*4];
					ti2 = tempwa[z2*4+1];
					ti3 = tempwa[z2*4+2];
					mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)] - ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
				} else
				{
					ti = tempwa[z2*4];
					ti2 = tempwa[z2*4+1];
					mprintf("[%s%d (%s) --> %s%d (%s)] - ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
				}
				ti = tempwa[z2*4+3];
				mprintf("%s%d (%s)\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName);
			}
			break;
	}
}


void CSDF::CreateCutPlane()
{
	int x, y;
	CxVector3 vec;

	try { m_pCutPlane = new C2DF(); } catch(...) { m_pCutPlane = NULL; }
	if (m_pCutPlane == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pCutPlane->m_fMinVal[0] = -m_fRadius;
	m_pCutPlane->m_fMinVal[1] = -m_fRadius;
	m_pCutPlane->m_fMaxVal[0] = m_fRadius;
	m_pCutPlane->m_fMaxVal[1] = m_fRadius;
	m_pCutPlane->m_iRes[0] = m_iCutPlaneResolution;
	m_pCutPlane->m_iRes[1] = m_iCutPlaneResolution;
	m_pCutPlane->SetLabelX("X [pm]");
	m_pCutPlane->SetLabelY("Y [pm]");
	m_pCutPlane->Create();

	vec[2] = 0;
	for (y=0;y<m_iCutPlaneResolution;y++)
	{
		vec[1] = (float)(((y+0.5)/m_pCutPlane->m_iRes[1])*(m_pCutPlane->m_fMaxVal[1]-m_pCutPlane->m_fMinVal[1])+m_pCutPlane->m_fMinVal[1]);
		for (x=0;x<m_iCutPlaneResolution;x++)
		{
			vec[0] = (float)(((x+0.5)/m_pCutPlane->m_iRes[0])*(m_pCutPlane->m_fMaxVal[0]-m_pCutPlane->m_fMinVal[0])+m_pCutPlane->m_fMinVal[0]);
			m_pCutPlane->AddToBin(x,y,m_pSDF->GetValue(vec));
		}
	}
	if (m_bCutPlaneShowAtoms)
	{
		m_pCutPlane->AddCircle(0,0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_fRadius*0.75,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_iColorR/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_iColorG/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_iColorB/255.0);
		m_pCutPlane->AddCircle(m_fAtom2PosX / m_fPosCounter,0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_fRadius*0.75,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_iColorR/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_iColorG/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_iColorB/255.0);
		m_pCutPlane->AddCircle(m_fAtom3PosX / m_fPosCounter,m_fAtom3PosY / m_fPosCounter,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_pElement->m_fRadius*0.75,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_pElement->m_iColorR/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_pElement->m_iColorG/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_pElement->m_iColorB/255.0);
	}
}


void CConditionGroup::Reset()
{
	int z;

	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->Reset();

	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
		m_iPassCounter[z] = 0;
}


void CConditionSubGroup::Reset()
{
	int z;

	for (z=0;z<m_oaConditions.GetSize();z++)
		((CNbSearch*)m_oaConditions[z])->Reset();
}


void CRDF::CopyFrom(CRDF *p)
{
	int z;
	CAtomGroup *ag;

	m_bAdaptive = p->m_bAdaptive;
	m_bRadialCorrect = p->m_bRadialCorrect;
	m_bSelf = p->m_bSelf;
	m_fMinDist = p->m_fMinDist;
	m_fMaxDist = p->m_fMaxDist;
	m_iCombinations = p->m_iCombinations;
	m_iHistogramRes = p->m_iHistogramRes;
	m_iRefAtomGes = p->m_iRefAtomGes;
	m_iRefOrSec[0] = p->m_iRefOrSec[0];
	m_iRefOrSec[1] = p->m_iRefOrSec[1];
	m_iResolution = p->m_iResolution;
	m_iShowAtomGes = p->m_iShowAtomGes;
	m_iShowMol = p->m_iShowMol;

	if (p->m_sName != NULL)
	{
		try { m_sName = new char[strlen(p->m_sName)+1]; } catch(...) { m_sName = NULL; }
		if (m_sName == NULL) NewException((double)(strlen(p->m_sName)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(m_sName,p->m_sName);
	} else m_sName = NULL;

	if (p->m_sShortName != NULL)
	{
		try { m_sShortName = new char[strlen(p->m_sShortName)+1]; } catch(...) { m_sShortName = NULL; }
		if (m_sShortName == NULL) NewException((double)(strlen(p->m_sShortName)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(m_sShortName,p->m_sShortName);
	} else m_sShortName = NULL;

	if (p->m_baDataEnabled != NULL)
	{
		try { m_baDataEnabled = new CxByteArray[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_baDataEnabled = NULL; }
		if (m_baDataEnabled == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
			m_baDataEnabled[z].CopyFrom(&p->m_baDataEnabled[z]);
	}

	if (p->m_faData != NULL)
	{
		try { m_faData = new CxDoubleArray[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_faData = NULL; }
		if (m_faData == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
			m_faData[z].CopyFrom(&p->m_faData[z]);
	}

	for (z=0;z<p->m_oaVectors.GetSize();z++)
	{
		try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
		if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		ag->CopyFrom((CAtomGroup*)p->m_oaVectors[z]);
		m_oaVectors.Add(ag);
	}

	m_faMinMaxDist.CopyFrom(&p->m_faMinMaxDist);

	if (p->m_pRDF != NULL)
	{
		try { m_pRDF = new CDF(); } catch(...) { m_pRDF = NULL; }
		if (m_pRDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pRDF->CopyFrom(p->m_pRDF);
	}
}


void CADF::CopyFrom(CADF *p)
{
	int z, z2;
	CAtomGroup *ag;

	m_bCosine = p->m_bCosine;
	m_bFoldAngle = p->m_bFoldAngle;
	m_bMirror = p->m_bMirror;
	m_bOrtho[0] = p->m_bOrtho[0];
	m_bOrtho[1] = p->m_bOrtho[1];
	m_iVecType[0] = p->m_iVecType[0];
	m_iVecType[1] = p->m_iVecType[1];
	m_bSelf = p->m_bSelf;
	m_bStat = p->m_bStat;
	m_fMinAngle = p->m_fMinAngle;
	m_fMaxAngle = p->m_fMaxAngle;
	m_iCombinations = p->m_iCombinations;
	m_iHistogramRes = p->m_iHistogramRes;
	for (z=0;z<2;z++)
		for (z2=0;z2<3;z2++)
			m_iRefOrSec[z][z2] = p->m_iRefOrSec[z][z2];
	m_iResolution = p->m_iResolution;
	m_iShowMol = p->m_iShowMol;

	if (p->m_sName != NULL)
	{
		try { m_sName = new char[strlen(p->m_sName)+1]; } catch(...) { m_sName = NULL; }
		if (m_sName == NULL) NewException((double)(strlen(p->m_sName)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(m_sName,p->m_sName);
	} else m_sName = NULL;

	if (p->m_sShortName != NULL)
	{
		try { m_sShortName = new char[strlen(p->m_sShortName)+1]; } catch(...) { m_sShortName = NULL; }
		if (m_sShortName == NULL) NewException((double)(strlen(p->m_sShortName)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(m_sShortName,p->m_sShortName);
	} else m_sShortName = NULL;

	if (p->m_baDataEnabled != NULL)
	{
		try { m_baDataEnabled = new CxByteArray[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_baDataEnabled = NULL; }
		if (m_baDataEnabled == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
			m_baDataEnabled[z].CopyFrom(&p->m_baDataEnabled[z]);
	}

	if (p->m_faData != NULL)
	{
		try { m_faData = new CxDoubleArray[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_faData = NULL; }
		if (m_faData == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
			m_faData[z].CopyFrom(&p->m_faData[z]);
	}

	for (z=0;z<p->m_oaVectors.GetSize();z++)
	{
		if (p->m_oaVectors[z] != NULL)
		{
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			ag->CopyFrom((CAtomGroup*)p->m_oaVectors[z]);
			m_oaVectors.Add(ag);
		} else m_oaVectors.Add(NULL);
	}

	m_faMinMaxAngle.CopyFrom(&p->m_faMinMaxAngle);

	if (p->m_pADF != NULL)
	{
		try { m_pADF = new CDF(); } catch(...) { m_pADF = NULL; }
		if (m_pADF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pADF->CopyFrom(p->m_pADF);
	}
}


void CConditionSubGroup::CopyFrom(CConditionSubGroup *p)
{
	int z;
	CNbSearch *nb;

	m_fPassed = p->m_fPassed;
	m_fTotal = p->m_fTotal;
	m_iCombinations = p->m_iCombinations;
	m_iNumber = p->m_iNumber;
	m_iShowMol = p->m_iShowMol;
	if (p->m_bTempPassed != NULL)
	{
		try { m_bTempPassed = new bool[g_iGesVirtAtomCount]; } catch(...) { m_bTempPassed = NULL; }
		if (m_bTempPassed == NULL) NewException((double)g_iGesVirtAtomCount*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_bTempPassed,p->m_bTempPassed,sizeof(bool)*g_iGesVirtAtomCount);
	}
	for (z=0;z<p->m_oaConditions.GetSize();z++)
	{
		try { nb = new CNbSearch(); } catch(...) { nb = NULL; }
		if (nb == NULL) NewException((double)sizeof(CNbSearch),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		nb->CopyFrom((CNbSearch*)p->m_oaConditions[z]);
		m_oaConditions.Add(nb);
	}
}


void CConditionGroup::CopyFrom(CConditionGroup *p)
{
	int z;
	CConditionSubGroup *sg;

	m_bInactive = p->m_bInactive;
	m_iShowMol = p->m_iShowMol;
	m_iRefMol = p->m_iRefMol;
	m_iHistoGes = p->m_iHistoGes;
	m_bAnyPassed = p->m_bAnyPassed;
	m_fTableGes = p->m_fTableGes;
	m_fPassed = p->m_fPassed;
	m_fTotal = p->m_fTotal;
	m_bInvertCondition = p->m_bInvertCondition;

	if (p->m_bAlwaysTrue != NULL)
	{
		try { m_bAlwaysTrue = new bool[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_bAlwaysTrue = NULL; }
		if (m_bAlwaysTrue == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_bAlwaysTrue,p->m_bAlwaysTrue,sizeof(bool)*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize());
	}
	if (p->m_iPassCounter != NULL)
	{
		try { m_iPassCounter = new long[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iPassCounter = NULL; }
		if (m_iPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_iPassCounter,p->m_iPassCounter,sizeof(long)*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize());
	}
	if (p->m_iOMPassCounter != NULL)
	{
		try { m_iOMPassCounter = new long[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iOMPassCounter = NULL; }
		if (m_iOMPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_iOMPassCounter,p->m_iOMPassCounter,sizeof(long)*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize());
	}
	if (p->m_iRMPassCounter != NULL)
	{
		try { m_iRMPassCounter = new long[((CMolecule*)g_oaMolecules[m_iRefMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iRMPassCounter = NULL; }
		if (m_iRMPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iRefMol])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_iRMPassCounter,p->m_iRMPassCounter,sizeof(long)*((CMolecule*)g_oaMolecules[m_iRefMol])->m_laSingleMolIndex.GetSize());
	}
	for (z=0;z<p->m_oaConditionSubGroups.GetSize();z++)
	{
		try { sg = new CConditionSubGroup(); } catch(...) { sg = NULL; }
		if (sg == NULL) NewException((double)sizeof(CConditionSubGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		sg->CopyFrom((CConditionSubGroup*)p->m_oaConditionSubGroups[z]);
		m_oaConditionSubGroups.Add(sg);
	}

	if (p->m_pTable != NULL)
	{
		try { m_pTable = new double[(((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_iCombinations+1) * (((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)]; } catch(...) { m_pTable = NULL; }
		if (m_pTable == NULL) NewException((double)(((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_iCombinations+1) * (((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_pTable,p->m_pTable,sizeof(double)*(((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_iCombinations+1) * (((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1));
	}

	if (p->m_pHistogram != NULL)
	{
		try { m_pHistogram = new unsigned long[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()+1]; } catch(...) { m_pHistogram = NULL; }
		if (m_pHistogram == NULL) NewException((double)(((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()+1)*sizeof(unsigned long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_pHistogram,p->m_pHistogram,sizeof(unsigned long)*(((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()+1));
	}
}


void CConditionGroup::Parse_OnlyValues()
{
	BTIN;
	int z;

	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
	{
		mprintf(YELLOW,">>> %d. set of conditions >>>\n\n",z+1);
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->Parse_OnlyValues();
		mprintf(YELLOW,"<<< End of %d. set of conditions <<<\n\n",z+1);
	}
}


void CConditionSubGroup::Parse_OnlyValues()
{
	BTIN;
	int z;

	for (z=0;z<m_oaConditions.GetSize();z++)
	{
		mprintf(YELLOW,"    >>> %d. condition within this set >>>\n\n",z+1);
		((CNbSearch*)m_oaConditions[z])->Parse_OnlyValues();
		mprintf(YELLOW,"\n    <<< End of %d. condition <<<\n\n",z+1);
	}
	BTOUT;
}


void CConditionGroup::CopyResults(CConditionGroup *p)
{
	int z;

	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->CopyResults((CConditionSubGroup*)p->m_oaConditionSubGroups[z]);
}


void CConditionSubGroup::CopyResults(CConditionSubGroup *p)
{
	int z;

	for (z=0;z<m_oaConditions.GetSize();z++)
		((CNbSearch*)m_oaConditions[z])->CopyResults((CNbSearch*)p->m_oaConditions[z]);
}


void CConditionGroup::ReScan(CSingleMolecule *rm)
{
	BXIN;
	int z, z2, t1, t2, i;

	if (m_bInactive)
		return;

	for (z=0;z<m_oaConditionSubGroups.GetSize();z++)
		((CConditionSubGroup*)m_oaConditionSubGroups[z])->ReScan(rm);

	if (m_oaConditionSubGroups.GetSize() == 2)
	{
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
		{
			t1 = 0;
			t2 = 0;
			for (z2=0;z2<((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_oaConditions.GetSize();z2++)
				t1 += ((CNbSearch*)((CConditionSubGroup*)m_oaConditionSubGroups[0])->m_oaConditions[z2])->m_iCombPassCount[z];
			for (z2=0;z2<((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_oaConditions.GetSize();z2++)
				t2 += ((CNbSearch*)((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_oaConditions[z2])->m_iCombPassCount[z];
			m_pTable[t1*(((CConditionSubGroup*)m_oaConditionSubGroups[1])->m_iCombinations+1)+t2]++;
			m_fTableGes++;
		}
	}

	i = 0;
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
	{
		m_fTotal++;
		if (m_bInvertCondition)
		{
			for (z2=0;z2<m_oaConditionSubGroups.GetSize();z2++)
				if (((CConditionSubGroup*)m_oaConditionSubGroups[z2])->Contains(z))
					goto _nopass;
			m_fPassed++;
			m_iPassCounter[z]++;
			m_bAnyPassed = true;
		} else
		{
			for (z2=0;z2<m_oaConditionSubGroups.GetSize();z2++)
				if (((CConditionSubGroup*)m_oaConditionSubGroups[z2])->Contains(z))
				{
//					mprintf("passed[%d]++\n",z);
					m_fPassed++;
					m_iPassCounter[z]++;
					m_bAnyPassed = true;
					i++;
					goto _nopass;
				}
		}
_nopass:;
	}
	m_pHistogram[i]++;
	m_iHistoGes++;

	BXOUT;
}


void CConditionSubGroup::ReScan(CSingleMolecule *rm)
{
	int z, z2;

	for (z=0;z<m_oaConditions.GetSize();z++)
		((CNbSearch*)m_oaConditions[z])->ReScan(rm);

	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
	{
		m_fTotal++;
		for (z2=0;z2<m_oaConditions.GetSize();z2++)
			if (!((CNbSearch*)m_oaConditions[z2])->m_bPassed[z])
				goto _nopass;
		m_fPassed++;
_nopass:;
	}
}


CNbSearch* CConditionGroup::AddSingleCondition(int rm, int sm, int gridmode)
{
	BTIN;
	CConditionSubGroup *sg;
	int z;

	m_iShowMol = sm;
	m_iRefMol = rm;

	try { m_bAlwaysTrue = new bool[((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_bAlwaysTrue = NULL; }
	if (m_bAlwaysTrue == NULL) NewException((double)((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_iPassCounter = new long[((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iPassCounter = NULL; }
	if (m_iPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_iOMPassCounter = new long[((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iOMPassCounter = NULL; }
	if (m_iOMPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize();z++)
	{
		m_bAlwaysTrue[z] = false;
		m_iPassCounter[z] = 0;
		m_iOMPassCounter[z] = 0;
	}

	try { m_iRMPassCounter = new long[((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iRMPassCounter = NULL; }
	if (m_iRMPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex.GetSize();z++)
		m_iRMPassCounter[z] = 0;

	try { sg = new CConditionSubGroup(); } catch(...) { sg = NULL; }
	if (sg == NULL) NewException((double)sizeof(CConditionSubGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	sg->m_iShowMol = m_iShowMol;
	sg->m_iNumber = m_oaConditionSubGroups.GetSize();
	m_oaConditionSubGroups.Add(sg);
	
	m_bInvertCondition = false;

	m_iHistoGes = 0;

	try { m_pHistogram = new unsigned long[((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()+1]; } catch(...) { m_pHistogram = NULL; }
	if (m_pHistogram == NULL) NewException((double)(((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()+1)*sizeof(unsigned long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<((CMolecule*)g_oaMolecules[sm])->m_laSingleMolIndex.GetSize()+1;z++)
		m_pHistogram[z] = 0;

	BTOUT;
	return sg->AddSingleCondition(rm,sm,gridmode);
}


CNbSearch* CConditionSubGroup::AddSingleCondition(int rm, int sm, int gridmode)
{
	BTIN;
	CNbSearch *n;

	m_iCombinations = 0;

	try { n = new CNbSearch(); } catch(...) { n = NULL; }
	if (n == NULL) NewException((double)sizeof(CNbSearch),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	n->m_iNumber = m_oaConditions.GetSize();
	m_oaConditions.Add(n);
	n->ParseGrid(rm,sm,gridmode);
	m_iCombinations += n->m_iCombinationsEnabled;

	try { m_bTempPassed = new bool[g_iGesVirtAtomCount]; } catch(...) { m_bTempPassed = NULL; }
	if (m_bTempPassed == NULL) NewException((double)g_iGesVirtAtomCount*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	BTOUT;
	return n;
}


void CObservation::BuildTimeDiff(CDF *df, bool ddf)
{
	int z, z2, z3, z4, i;
	CxFloatArray *ptfa;
	double tf, tfa, tfsq, tfs, tf0;
//	char buf[1024];
	CxString buf;

	try { df->m_pTimeDiff = new CDF(); } catch(...) { df->m_pTimeDiff = NULL; }
	if (df->m_pTimeDiff == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	df->m_pTimeDiff->m_iResolution = m_iTimeDiffDepth;
	df->m_pTimeDiff->m_fMinVal = 0;
	df->m_pTimeDiff->m_fMaxVal = m_iTimeDiffDepth * g_fTimestepLength / 1000.0;
	df->m_pTimeDiff->Create();
	df->m_pTimeDiff->SetLabelX("Tau [ps]");
//	sprintf(buf,"Delta ");
//	strcat(buf,df->m_sLabelX);
	buf.sprintf("Delta ");
	buf.strcat(df->m_sLabelX);
	df->m_pTimeDiff->SetLabelY(buf);

	try { df->m_pTimeDiffAbs = new CDF(); } catch(...) { df->m_pTimeDiffAbs = NULL; }
	if (df->m_pTimeDiffAbs == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	df->m_pTimeDiffAbs->m_iResolution = m_iTimeDiffDepth;
	df->m_pTimeDiffAbs->m_fMinVal = 0;
	df->m_pTimeDiffAbs->m_fMaxVal = m_iTimeDiffDepth * g_fTimestepLength / 1000.0;
	df->m_pTimeDiffAbs->Create();
	df->m_pTimeDiffAbs->SetLabelX("Tau [ps]");
	df->m_pTimeDiffAbs->SetLabelY(buf);

//	sprintf(buf,"Delta Square ");
//	strcat(buf,df->m_sLabelX);
	buf.sprintf("Delta Square ");
	buf.strcat(df->m_sLabelX);

	try { df->m_pTimeDiffSqr = new CDF(); } catch(...) { df->m_pTimeDiffSqr = NULL; }
	if (df->m_pTimeDiffSqr == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	df->m_pTimeDiffSqr->m_iResolution = m_iTimeDiffDepth;
	df->m_pTimeDiffSqr->m_fMinVal = 0;
	df->m_pTimeDiffSqr->m_fMaxVal = m_iTimeDiffDepth * g_fTimestepLength / 1000.0;
	df->m_pTimeDiffSqr->Create();
	df->m_pTimeDiffSqr->SetLabelX("Tau [ps]");
	df->m_pTimeDiffSqr->SetLabelY(buf);

	if (m_b3DTimeDiff)
	{
//		sprintf(buf,"Delta ");
//		strcat(buf,df->m_sLabelX);
		buf.sprintf("Delta ");
		buf.strcat(df->m_sLabelX);

		try { df->m_p3DTimeDiff = new C2DF(); } catch(...) { df->m_p3DTimeDiff = NULL; }
		if (df->m_p3DTimeDiff == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		df->m_p3DTimeDiff->m_iRes[0] = m_iTimeDiffDepth/m_iTimeDiffStride3D;
		df->m_p3DTimeDiff->m_iRes[1] = m_iTimeDiffRes3D;
		df->m_p3DTimeDiff->m_fMinVal[0] = 0;
		df->m_p3DTimeDiff->m_fMaxVal[0] = m_iTimeDiffDepth * g_fTimestepLength / 1000.0;
		df->m_p3DTimeDiff->m_fMinVal[1] = m_fTimeDiffMinVal3D;
		df->m_p3DTimeDiff->m_fMaxVal[1] = m_fTimeDiffMaxVal3D;
		df->m_p3DTimeDiff->Create();
		df->m_p3DTimeDiff->SetLabelX("Tau [ps]");
		df->m_p3DTimeDiff->SetLabelY(df->m_sLabelX);
		df->m_p3DTimeDiff->SetLabelZ(buf);

		try { df->m_p3DTimeDiffAbs = new C2DF(); } catch(...) { df->m_p3DTimeDiffAbs = NULL; }
		if (df->m_p3DTimeDiffAbs == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		df->m_p3DTimeDiffAbs->m_iRes[0] = m_iTimeDiffDepth/m_iTimeDiffStride3D;
		df->m_p3DTimeDiffAbs->m_iRes[1] = m_iTimeDiffRes3D;
		df->m_p3DTimeDiffAbs->m_fMinVal[0] = 0;
		df->m_p3DTimeDiffAbs->m_fMaxVal[0] = m_iTimeDiffDepth * g_fTimestepLength / 1000.0;
		df->m_p3DTimeDiffAbs->m_fMinVal[1] = m_fTimeDiffMinVal3D;
		df->m_p3DTimeDiffAbs->m_fMaxVal[1] = m_fTimeDiffMaxVal3D;
		df->m_p3DTimeDiffAbs->Create();
		df->m_p3DTimeDiffAbs->SetLabelX("Tau [ps]");
		df->m_p3DTimeDiffAbs->SetLabelY(df->m_sLabelX);
		df->m_p3DTimeDiffAbs->SetLabelZ(buf);

		if (m_iTimeDiffDistSteps != 0)
		{
			try { df->m_pTimeDiffDistPairs = new C2DF*[m_iTimeDiffDepth/m_iTimeDiffStride3D/m_iTimeDiffDistSteps]; } catch(...) { df->m_pTimeDiffDistPairs = NULL; }
			if (df->m_pTimeDiffDistPairs == NULL) NewException((double)m_iTimeDiffDepth/m_iTimeDiffStride3D/m_iTimeDiffDistSteps*sizeof(C2DF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z=0;z<m_iTimeDiffDepth/m_iTimeDiffStride3D/m_iTimeDiffDistSteps;z++)
			{
				try { df->m_pTimeDiffDistPairs[z] = new C2DF(); } catch(...) { df->m_pTimeDiffDistPairs[z] = NULL; }
				if (df->m_pTimeDiffDistPairs[z] == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				df->m_pTimeDiffDistPairs[z]->m_iRes[0] = m_iTimeDiffDistResX;
				df->m_pTimeDiffDistPairs[z]->m_iRes[1] = m_iTimeDiffDistResY;
				df->m_pTimeDiffDistPairs[z]->m_fMinVal[0] = m_fTimeDiffDistMinValX;
				df->m_pTimeDiffDistPairs[z]->m_fMaxVal[0] = m_fTimeDiffDistMaxValX;
				df->m_pTimeDiffDistPairs[z]->m_fMinVal[1] = m_fTimeDiffDistMinValY;
				df->m_pTimeDiffDistPairs[z]->m_fMaxVal[1] = m_fTimeDiffDistMaxValY;
				df->m_pTimeDiffDistPairs[z]->Create();
				df->m_pTimeDiffDistPairs[z]->SetLabelX(df->m_sLabelX);
				df->m_pTimeDiffDistPairs[z]->SetLabelY(df->m_sLabelX);
			}
		}

		try { df->m_pTimeDiffDist3DF = new C3DF<double>(); } catch(...) { df->m_pTimeDiffDist3DF = NULL; }
		if (df->m_pTimeDiffDist3DF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		df->m_pTimeDiffDist3DF->m_iRes[0] = m_iTimeDiffDistResX;
		df->m_pTimeDiffDist3DF->m_iRes[1] = m_iTimeDiffDistResY;
		df->m_pTimeDiffDist3DF->m_iRes[2] = m_iTimeDiffDepth/m_iTimeDiffStride3D;
		df->m_pTimeDiffDist3DF->m_fMinVal[0] = m_fTimeDiffDistMinValX;
		df->m_pTimeDiffDist3DF->m_fMaxVal[0] = m_fTimeDiffDistMaxValX;
		df->m_pTimeDiffDist3DF->m_fMinVal[1] = m_fTimeDiffDistMinValY;
		df->m_pTimeDiffDist3DF->m_fMaxVal[1] = m_fTimeDiffDistMaxValY;
		df->m_pTimeDiffDist3DF->m_fMinVal[2] = 0;
		df->m_pTimeDiffDist3DF->m_fMaxVal[2] = m_iTimeDiffDepth * g_fTimestepLength / 1000.0;
		df->m_pTimeDiffDist3DF->Create();
		df->m_pTimeDiffDist3DF->SetLabelX(df->m_sLabelX);
		df->m_pTimeDiffDist3DF->SetLabelY(df->m_sLabelX);

//		sprintf(buf,"Delta Square ");
//		strcat(buf,df->m_sLabelX);
		buf.sprintf("Delta Square ");
		buf.strcat(df->m_sLabelX);

		try { df->m_p3DTimeDiffSqr = new C2DF(); } catch(...) { df->m_p3DTimeDiffSqr = NULL; }
		if (df->m_p3DTimeDiffSqr == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		df->m_p3DTimeDiffSqr->m_iRes[0] = m_iTimeDiffDepth/m_iTimeDiffStride3D;
		df->m_p3DTimeDiffSqr->m_iRes[1] = m_iTimeDiffRes3D;
		df->m_p3DTimeDiffSqr->m_fMinVal[0] = 0;
		df->m_p3DTimeDiffSqr->m_fMaxVal[0] = m_iTimeDiffDepth * g_fTimestepLength / 1000.0;
		df->m_p3DTimeDiffSqr->m_fMinVal[1] = m_fTimeDiffMinVal3D;
		df->m_p3DTimeDiffSqr->m_fMaxVal[1] = m_fTimeDiffMaxVal3D;
		df->m_p3DTimeDiffSqr->Create();
		df->m_p3DTimeDiffSqr->SetLabelX("Tau [ps]");
		df->m_p3DTimeDiffSqr->SetLabelY(df->m_sLabelX);
		df->m_p3DTimeDiffSqr->SetLabelZ(buf);

		try { df->m_p3DTimeDiffT = new C2DF(); } catch(...) { df->m_p3DTimeDiffT = NULL; }
		if (df->m_p3DTimeDiffT == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		df->m_p3DTimeDiffT->m_iRes[0] = m_iTimeDiffDepth/m_iTimeDiffStride3D;
		df->m_p3DTimeDiffT->m_iRes[1] = m_iTimeDiffRes3D;
		df->m_p3DTimeDiffT->m_fMinVal[0] = 0;
		df->m_p3DTimeDiffT->m_fMaxVal[0] = m_iTimeDiffDepth * g_fTimestepLength / 1000.0;
		df->m_p3DTimeDiffT->m_fMinVal[1] = m_fTimeDiffMinVal3D;
		df->m_p3DTimeDiffT->m_fMaxVal[1] = m_fTimeDiffMaxVal3D;
		df->m_p3DTimeDiffT->Create();
	}

	for (z2=0;z2<df->m_oaTimeDiffBuf.GetSize();z2++)
	{
		mprintf("      %4d/%d 2D  [",z2+1,df->m_oaTimeDiffBuf.GetSize());
		ptfa = (CxFloatArray*)df->m_oaTimeDiffBuf[z2];
		tfs = m_iTimeDiffDepth/50.0;
		for (z3=0;z3<m_iTimeDiffDepth;z3++) // Das ist das Tau
		{
			if (fmod(z3,tfs) < 1.0)
				mprintf(WHITE,"#");
			tf = 0;
			tfa = 0;
			tfsq = 0;
			for (z4=0;z4<(int)ptfa->GetSize()-z3-1;z4++) // Das ist der Startpunkt
			{
				tf0 = (*ptfa)[z3+z4]-(*ptfa)[z4];
				if (ddf)
				{
					if (tf0 > 180.0)
						tf0 -= 360.0;
					if (tf0 <= -180.0)
						tf0 += 360.0;
				}
				tfa += fabs(tf0);
				tf += tf0;
				tfsq += tf0 * tf0;
			}
			df->m_pTimeDiff->AddToBin_Int(z3,tf/(ptfa->GetSize()-z3-1));
			df->m_pTimeDiffAbs->AddToBin_Int(z3,tfa/(ptfa->GetSize()-z3-1));
			df->m_pTimeDiffSqr->AddToBin_Int(z3,tfsq/(ptfa->GetSize()-z3-1));

			df->m_pTimeDiff->m_fBinEntries += ptfa->GetSize()-z3-2;
			df->m_pTimeDiffAbs->m_fBinEntries += ptfa->GetSize()-z3-2;
			df->m_pTimeDiffSqr->m_fBinEntries += ptfa->GetSize()-z3-2;
		}
		mprintf("]\n");
		if (m_b3DTimeDiff)
		{
			mprintf("      %4d/%d 3D  [",z2+1,df->m_oaTimeDiffBuf.GetSize());
			tfs = m_iTimeDiffDepth/50.0;
			for (z3=0;z3<m_iTimeDiffDepth;z3+=m_iTimeDiffStride3D) // Das ist das Tau
			{
				if (fmod(z3,tfs) < 1.0)
					mprintf(WHITE,"#");
				i = z3/m_iTimeDiffStride3D;
				for (z4=0;z4<(int)ptfa->GetSize()-z3-1;z4++) // Das ist der Startpunkt
				{
					tf0 = (*ptfa)[z3+z4]-(*ptfa)[z4];
					if (ddf)
					{
						if (tf0 > 180.0)
							tf0 -= 360.0;
						if (tf0 <= -180.0)
							tf0 += 360.0;
					}
					df->m_p3DTimeDiff->AddToBin_IntX(i,(*ptfa)[z4],tf0);
					df->m_p3DTimeDiffAbs->AddToBin_IntX(i,(*ptfa)[z4],fabs(tf0));
					df->m_p3DTimeDiffSqr->AddToBin_IntX(i,(*ptfa)[z4],tf0*tf0);
					df->m_p3DTimeDiffT->AddToBin_IntX(i,(*ptfa)[z4],1.0);
					df->m_pTimeDiffDist3DF->AddToBin_IntZ((*ptfa)[z4],(*ptfa)[z3+z4],i);
					if (m_iTimeDiffDistSteps != 0)
						if ((i % m_iTimeDiffDistSteps) == 0)
							df->m_pTimeDiffDistPairs[i/m_iTimeDiffDistSteps]->AddToBin((*ptfa)[z4],(*ptfa)[z3+z4]);
				}
			}
			mprintf("]\n");
		}
	}
	for (z3=0;z3<m_iTimeDiffDepth;z3++) // Das ist das Tau
	{
		df->m_pTimeDiff->m_pBin[z3] /= df->m_oaTimeDiffBuf.GetSize();
		df->m_pTimeDiffAbs->m_pBin[z3] /= df->m_oaTimeDiffBuf.GetSize();
		df->m_pTimeDiffSqr->m_pBin[z3] /= df->m_oaTimeDiffBuf.GetSize();
	}
	if (m_b3DTimeDiff)
	{
		for (z2=0;z2<df->m_p3DTimeDiff->m_iRes[0]*df->m_p3DTimeDiff->m_iRes[1];z2++)
		{
			if (df->m_p3DTimeDiffT->m_pBin[z2] != 0)
			{
				df->m_p3DTimeDiff->m_pBin[z2] /= df->m_p3DTimeDiffT->m_pBin[z2];
				df->m_p3DTimeDiffAbs->m_pBin[z2] /= df->m_p3DTimeDiffT->m_pBin[z2];
				df->m_p3DTimeDiffSqr->m_pBin[z2] /= df->m_p3DTimeDiffT->m_pBin[z2];
			}
		}
		if (m_iTimeDiffDistSteps != 0)
			for (z=0;z<m_iTimeDiffDepth/m_iTimeDiffStride3D/m_iTimeDiffDistSteps;z++)
				df->m_pTimeDiffDistPairs[z]->NormalizeBinIntegral(1000000.0);
	}
}


void CObservation::WriteTimeDiff(CDF *df, const char *anaup, const char *analow, const char *name, const char *multibuf, bool ddf)
{
//	char buf[1024];
	CxString buf;
	int z, z2;
	C3DF<double> *temp3DF;

	mprintf("    Creating temporal difference plot...\n");
	BuildTimeDiff(df,ddf);
	mprintf("      (%.0f bin entries)\n",df->m_pTimeDiff->m_fBinEntries);
//	sprintf(buf,"%s_timediff_%s%s.csv",analow,name,multibuf);
	buf.sprintf("%s_timediff_%s%s.csv",analow,name,multibuf);
	mprintf("      Saving %s temporal difference plot as \"%s\"...\n",anaup,(const char*)buf);
	df->m_pTimeDiff->Write("",buf,"",true);
//	sprintf(buf,"%s_timediff_%s%s.agr",analow,name,multibuf);
	buf.sprintf("%s_timediff_%s%s.agr",analow,name,multibuf);
	mprintf("      Saving %s temporal difference plot AGR file as \"%s\"...\n",anaup,(const char*)buf);
	df->m_pTimeDiff->WriteAgr("",buf,"",name,false);

//	sprintf(buf,"%s_timediff_%s%s_abs.csv",analow,name,multibuf);
	buf.sprintf("%s_timediff_%s%s_abs.csv",analow,name,multibuf);
	mprintf("      Saving %s absolute temporal difference plot as \"%s\"...\n",anaup,(const char*)buf);
	df->m_pTimeDiffAbs->Write("",buf,"",true);
//	sprintf(buf,"%s_timediff_%s%s_abs.agr",analow,name,multibuf);
	buf.sprintf("%s_timediff_%s%s_abs.agr",analow,name,multibuf);
	mprintf("      Saving %s absolute temporal difference plot AGR file as \"%s\"...\n",anaup,(const char*)buf);
	df->m_pTimeDiffAbs->WriteAgr("",buf,"",name,false);

//	sprintf(buf,"%s_timediff_%s%s_sqr.csv",analow,name,multibuf);
	buf.sprintf("%s_timediff_%s%s_sqr.csv",analow,name,multibuf);
	mprintf("      Saving %s squared temporal difference plot as \"%s\"...\n",anaup,(const char*)buf);
	df->m_pTimeDiffSqr->Write("",buf,"",true);
//	sprintf(buf,"%s_timediff_%s%s_sqr.agr",analow,name,multibuf);
	buf.sprintf("%s_timediff_%s%s_sqr.agr",analow,name,multibuf);
	mprintf("      Saving %s squared temporal difference plot AGR file as \"%s\"...\n",anaup,(const char*)buf);
	df->m_pTimeDiffSqr->WriteAgr("",buf,"",name,false);

	if (m_b3DTimeDiff)
	{
		mprintf("    Writing temporal 3D difference plot...\n");
		mprintf("      (%.0f bin entries)\n",df->m_p3DTimeDiff->m_fBinEntries);
//		sprintf(buf,"%s_timediff_%s%s_triples.csv",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_triples.csv",analow,name,multibuf);
		mprintf("      Saving %s temporal difference plot triples as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiff->Write("",buf,"");
//		sprintf(buf,"%s_timediff_%s%s_matrix.csv",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_matrix.csv",analow,name,multibuf);
		mprintf("      Saving %s temporal difference plot matrix as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiff->WriteCSV("",buf,"");
//		sprintf(buf,"%s_timediff_%s%s.nb",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s.nb",analow,name,multibuf);
		mprintf("      Saving %s temporal difference plot Mathematica Notebook as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiff->WriteMathematicaNb("",buf,"",false);
//		sprintf(buf,"%s_timediff_%s%s",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s",analow,name,multibuf);
		mprintf("      Saving %s temporal difference plot Gnuplot Input as \"%s.gp\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiff->WriteGnuplotInput("",buf,"",false);

//		sprintf(buf,"%s_timediff_%s%s_triples_abs.csv",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_triples_abs.csv",analow,name,multibuf);
		mprintf("      Saving %s absolute temporal difference plot triples as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiffAbs->Write("",buf,"");
//		sprintf(buf,"%s_timediff_%s%s_matrix_abs.csv",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_matrix_abs.csv",analow,name,multibuf);
		mprintf("      Saving %s absolute temporal difference plot matrix as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiffAbs->WriteCSV("",buf,"");
//		sprintf(buf,"%s_timediff_%s%s_abs.nb",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_abs.nb",analow,name,multibuf);
		mprintf("      Saving %s absolute temporal difference plot Mathematica Notebook as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiffAbs->WriteMathematicaNb("",buf,"",false);
//		sprintf(buf,"%s_timediff_%s%s_abs",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_abs",analow,name,multibuf);
		mprintf("      Saving %s absolute temporal difference plot Gnuplot Input as \"%s.gp\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiffAbs->WriteGnuplotInput("",buf,"",false);

//		sprintf(buf,"%s_timediff_%s%s_triples_sqr.csv",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_triples_sqr.csv",analow,name,multibuf);
		mprintf("      Saving %s squared temporal difference plot triples as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiffSqr->Write("",buf,"");
//		sprintf(buf,"%s_timediff_%s%s_matrix_sqr.csv",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_matrix_sqr.csv",analow,name,multibuf);
		mprintf("      Saving %s squared temporal difference plot matrix as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiffSqr->WriteCSV("",buf,"");
//		sprintf(buf,"%s_timediff_%s%s_sqr.nb",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_sqr.nb",analow,name,multibuf);
		mprintf("      Saving %s squared temporal difference plot Mathematica Notebook as \"%s\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiffSqr->WriteMathematicaNb("",buf,"",false);
//		sprintf(buf,"%s_timediff_%s%s_sqr",analow,name,multibuf);
		buf.sprintf("%s_timediff_%s%s_sqr",analow,name,multibuf);
		mprintf("      Saving %s squared temporal difference plot Gnuplot Input as \"%s.gp\"...\n",anaup,(const char*)buf);
		df->m_p3DTimeDiffSqr->WriteGnuplotInput("",buf,"",false);

		if (m_iTimeDiffDistSteps != 0)
		{
			for (z=0;z<m_iTimeDiffDepth/m_iTimeDiffStride3D/m_iTimeDiffDistSteps;z++)
			{
//				sprintf(buf,"%s_timediff_beforeafter_%s%s_tau%.3f_triples.csv",analow,name,multibuf,z*m_iTimeDiffStride3D*m_iTimeDiffDistSteps*g_fTimestepLength/1000.0);
				buf.sprintf("%s_timediff_beforeafter_%s%s_tau%.3f_triples.csv",analow,name,multibuf,z*m_iTimeDiffStride3D*m_iTimeDiffDistSteps*g_fTimestepLength/1000.0);
				mprintf("      Saving %s before/after plot triples as \"%s\"...\n",anaup,(const char*)buf);
				df->m_pTimeDiffDistPairs[z]->Write("",buf,"");
//				sprintf(buf,"%s_timediff_beforeafter_%s%s_tau%.3f_matrix.csv",analow,name,multibuf,z*m_iTimeDiffStride3D*m_iTimeDiffDistSteps*g_fTimestepLength/1000.0);
				buf.sprintf("%s_timediff_beforeafter_%s%s_tau%.3f_matrix.csv",analow,name,multibuf,z*m_iTimeDiffStride3D*m_iTimeDiffDistSteps*g_fTimestepLength/1000.0);
				mprintf("      Saving %s before/after plot matrix as \"%s\"...\n",anaup,(const char*)buf);
				df->m_pTimeDiffDistPairs[z]->WriteCSV("",buf,"");
//				sprintf(buf,"%s_timediff_beforeafter_%s%s_tau%.3f.nb",analow,name,multibuf,z*m_iTimeDiffStride3D*m_iTimeDiffDistSteps*g_fTimestepLength/1000.0);
				buf.sprintf("%s_timediff_beforeafter_%s%s_tau%.3f.nb",analow,name,multibuf,z*m_iTimeDiffStride3D*m_iTimeDiffDistSteps*g_fTimestepLength/1000.0);
				mprintf("      Saving %s before/after plot Mathematica Notebook as \"%s\"...\n",anaup,(const char*)buf);
				df->m_pTimeDiffDistPairs[z]->WriteMathematicaNb("",buf,"",false);
//				sprintf(buf,"%s_timediff_beforeafter_%s%s_tau%.3f",analow,name,multibuf,z*m_iTimeDiffStride3D*m_iTimeDiffDistSteps*g_fTimestepLength/1000.0);
				buf.sprintf("%s_timediff_beforeafter_%s%s_tau%.3f",analow,name,multibuf,z*m_iTimeDiffStride3D*m_iTimeDiffDistSteps*g_fTimestepLength/1000.0);
				mprintf("      Saving %s before/after plot Gnuplot Input as \"%s.gp\"...\n",anaup,(const char*)buf);
				df->m_pTimeDiffDistPairs[z]->WriteGnuplotInput("",buf,"",false);
			}
		}
		for (z2=0;z2<=3;z2++)
		{
			try { temp3DF = new C3DF<double>(); } catch(...) { temp3DF = NULL; }
			if (temp3DF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			temp3DF->CopyFrom(df->m_pTimeDiffDist3DF);
			if (z2 != 0)
			{
				temp3DF->Smooth(z2);
//				sprintf(buf,".s%d%s.plt",z2,multibuf);
				buf.sprintf(".s%d%s.plt",z2,multibuf);
			} else
//				sprintf(buf,"%s.plt",multibuf);
				buf.sprintf("%s.plt",multibuf);

			mprintf("    Saving 3D Plot as \"%s%s\"...\n",name,(const char*)buf);
			temp3DF->WritePLT("",name,buf,false);

			if (z2 != 0)
//				sprintf(buf,".s%d%s.cube",z2,multibuf);
				buf.sprintf(".s%d%s.cube",z2,multibuf);
			else
//				sprintf(buf,"%s.cube",multibuf);
				buf.sprintf("%s.cube",multibuf);

			mprintf("    Saving 3D Plot as \"%s%s\"...\n",name,(const char*)buf);
			temp3DF->WriteCube("",name,buf,false);
		}
	}
}


void CObservation::CreateTimeDiff(CDF *df, int comb)
{
	int z2;
	CxFloatArray *ptfa;

	if (m_bSelf)
	{
		for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*comb;z2++)
		{
			try { ptfa = new CxFloatArray("CObservation::CreateTimeDiff():ptfa"); } catch(...) { ptfa = NULL; }
			if (ptfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			if (g_iTrajSteps != -1)
			{
				ptfa->SetMaxSize((long)(g_iTrajSteps*1.1));
				ptfa->SetGrow((long)(g_iTrajSteps*0.1));
			} else ptfa->SetGrow(1000);
			df->m_oaTimeDiffBuf.Add(ptfa);
		}
	} else
	{
		for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*comb;z2++)
		{
			try { ptfa = new CxFloatArray("CObservation::CreateTimeDiff():ptfa"); } catch(...) { ptfa = NULL; }
			if (ptfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			if (g_iTrajSteps != -1)
			{
				ptfa->SetMaxSize((long)(g_iTrajSteps*1.1));
				ptfa->SetGrow((long)(g_iTrajSteps*0.1));
			} else ptfa->SetGrow(1000);
			df->m_oaTimeDiffBuf.Add(ptfa);
		}
	}
}


void CMSD::WriteSplit(const char *s)
{
	FILE *a;
	int z, z2, z3;
	CAF *af;

	a = OpenFileWrite(s,true);
	mfprintf(a,"# tau [ps];  Total MSD [pm^2]");

	for (z=0;z<m_pAtomGroup->m_oaAtoms.GetSize();z++)
		for (z2=0;z2<((CxIntArray*)m_pAtomGroup->m_oaAtoms[z])->GetSize();z2++)
			for (z3=0;z3<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z3++)
			{
				if (m_pAtomGroup->m_baRealAtomType[z] == g_iVirtAtomType)
					mfprintf(a,";  Mol. %d #%d",z3+1,((CxIntArray*)m_pAtomGroup->m_oaAtoms[z])->GetAt(z2)+1);
						else mfprintf(a,";  Mol. %d %s%d",z3+1,((CAtom*)g_oaAtoms[m_pAtomGroup->m_baRealAtomType[z]])->m_sName,((CxIntArray*)m_pAtomGroup->m_oaAtoms[z])->GetAt(z2)+1);
			}
	mfprintf(a,"\n");

	for (z=0;z<m_pMSD->m_iResolution;z++)
	{
		mfprintf(a,"%f;  %f",m_pMSD->m_fMinVal+z*(m_pMSD->m_fMaxVal-m_pMSD->m_fMinVal)/m_pMSD->m_iResolution,m_pMSD->m_pBin[z]);

		for (z2=0;z2<m_pAtomGroup->m_iAtomGes;z2++)
			for (z3=0;z3<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z3++)
			{
				af = m_pSplitMSD[z2*((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()+z3];
				mfprintf(a,";  %f",af->m_pBin[z]);
			}
		mfprintf(a,"\n");
	}

	fclose(a);
}
