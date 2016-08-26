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

#include "nbsearch.h"
#include "travis.h"
#include "maintools.h"


CNbPair::CNbPair()
{
	m_fDistances = NULL;
	m_fAngles = NULL;
	m_bDistPassed = NULL;
	m_bAnglePassed = NULL;
	m_bAnyAnglePassed = false;
	m_bAnyDistPassed = false;
	m_laDistAtomList.SetName("CNbPair::m_laDistAtomList");
	m_laAngleAtomList.SetName("CNbPair::m_laAngleAtomList");
}


CNbPair::~CNbPair()
{
}


void CNbPair::Check(CNbSearch *parent, CTimeStep *ts, CSingleMolecule *rm, CSingleMolecule *sm)
{
	int z, z2;
	CxVector3 vec1, vec2, veca, vecb, vecc;

	if (parent->m_pRDF != NULL)
	{
		m_bAnyDistPassed = false;
		m_fMinDist = 1E30f;
		m_laDistAtomList.RemoveAll_KeepSize();
		parent->m_pRDF->BuildAtomList(rm,sm,&m_laDistAtomList);
		for (z=0;z<m_laDistAtomList.GetSize()/2;z++)
		{
			m_fDistances[z] = FoldedLength(ts->m_vaCoords[m_laDistAtomList[z*2]]-ts->m_vaCoords[m_laDistAtomList[z*2+1]]);
			if (m_fDistances[z] < m_fMinDist)
			{
				m_fMinDist = m_fDistances[z];
				m_iMinDistIndex = z;
			}

			m_bDistPassed[z] = false;

			for (z2=0;z2<parent->m_pRDF->m_faMinMaxDist.GetSize();z2+=2)
			{
				if ((m_fDistances[z] >= parent->m_pRDF->m_faMinMaxDist[z2]) && (m_fDistances[z] <= parent->m_pRDF->m_faMinMaxDist[z2+1]))
				{
					m_bAnyDistPassed = true;
					m_bDistPassed[z] = true;
					break;
				} 
			}
		}
	} else
	{
		m_bAnyDistPassed = true;
		m_bDistPassed[0] = true;
	}

	if (parent->m_pADF != NULL)
	{
		m_bAnyAnglePassed = false;
		m_laAngleAtomList.RemoveAll_KeepSize();
		parent->m_pADF->BuildAtomList(rm,sm,&m_laAngleAtomList);
		for (z=0;z<m_laAngleAtomList.GetSize()/6;z++)
		{
			if (parent->m_pADF->m_bOrtho[0])
			{
				veca = ts->m_vaCoords[m_laAngleAtomList[z*6]];
				vecb = ts->m_vaCoords[m_laAngleAtomList[z*6+1]];
				vecc = ts->m_vaCoords[m_laAngleAtomList[z*6+2]];
				vec1 = CrossP(FoldVector(vecb-veca),FoldVector(vecc-veca));
			} else
			{
				veca = ts->m_vaCoords[m_laAngleAtomList[z*6]];
				vecb = ts->m_vaCoords[m_laAngleAtomList[z*6+1]];
				vec1 = FoldVector(vecb-veca);
			}
			if (parent->m_pADF->m_bOrtho[1])
			{
				veca = ts->m_vaCoords[m_laAngleAtomList[z*6+3]];
				vecb = ts->m_vaCoords[m_laAngleAtomList[z*6+4]];
				vecc = ts->m_vaCoords[m_laAngleAtomList[z*6+5]];
				vec2 = CrossP(FoldVector(vecb-veca),FoldVector(vecc-veca));
			} else
			{
				veca = ts->m_vaCoords[m_laAngleAtomList[z*6+3]];
				vecb = ts->m_vaCoords[m_laAngleAtomList[z*6+4]];
				vec2 = FoldVector(vecb-veca);
			}
			m_fAngles[z] = Angle_Deg(vec1,vec2);
			m_bAnglePassed[z] = false;
			for (z2=0;z2<parent->m_pADF->m_faMinMaxAngle.GetSize();z2+=2)
			{
				if ((m_fAngles[z] >= parent->m_pADF->m_faMinMaxAngle[z2]) && (m_fAngles[z] <= parent->m_pADF->m_faMinMaxAngle[z2+1]))
				{
					m_bAnyAnglePassed = true;
					m_bAnglePassed[z] = true;
					break;
				}
			}
		}
	} else
	{
		m_bAnyAnglePassed = true;
		m_bAnglePassed[0] = true;
	}
//	mprintf("Pair: Angle=%f --> %s.\n",m_fAngles[0],m_bAnglePassed[0]?"passed.":"failed.");
}


void CNbPair::PreCheck(CNbSearch *parent, CTimeStep *ts, CSingleMolecule *rm, CSingleMolecule *sm)
{
	int z, z2;
	CxVector3 vec1, vec2, veca, vecb, vecc;

	if (parent->m_pRDF != NULL)
	{
		m_fMinDist = 1E30f;
		m_laDistAtomList.RemoveAll_KeepSize();
		parent->m_pRDF->BuildAtomList(rm,sm,&m_laDistAtomList);
		for (z=0;z<m_laDistAtomList.GetSize()/2;z++)
		{
			m_fDistances[z] = FoldedLength(ts->m_vaCoords[m_laDistAtomList[z*2]]-ts->m_vaCoords[m_laDistAtomList[z*2+1]]);
			if (m_fDistances[z] < m_fMinDist)
			{
				m_fMinDist = m_fDistances[z];
				m_iMinDistIndex = z;
			}
		}
	} else
	{
		m_bAnyDistPassed = true;
		m_bDistPassed[0] = true;
	}

	if (parent->m_pADF != NULL)
	{
		m_bAnyAnglePassed = false;
		m_laAngleAtomList.RemoveAll_KeepSize();
		parent->m_pADF->BuildAtomList(rm,sm,&m_laAngleAtomList);
		for (z=0;z<m_laAngleAtomList.GetSize()/6;z++)
		{
			if (parent->m_pADF->m_bOrtho[0])
			{
				veca = ts->m_vaCoords[m_laAngleAtomList[z*6]];
				vecb = ts->m_vaCoords[m_laAngleAtomList[z*6+1]];
				vecc = ts->m_vaCoords[m_laAngleAtomList[z*6+2]];
				vec1 = CrossP(FoldVector(vecb-veca),FoldVector(vecc-veca));
			} else
			{
				veca = ts->m_vaCoords[m_laAngleAtomList[z*6]];
				vecb = ts->m_vaCoords[m_laAngleAtomList[z*6+1]];
				vec1 = FoldVector(vecb-veca);
			}
			if (parent->m_pADF->m_bOrtho[1])
			{
				veca = ts->m_vaCoords[m_laAngleAtomList[z*6+3]];
				vecb = ts->m_vaCoords[m_laAngleAtomList[z*6+4]];
				vecc = ts->m_vaCoords[m_laAngleAtomList[z*6+5]];
				vec2 = CrossP(FoldVector(vecb-veca),FoldVector(vecc-veca));
			} else
			{
				veca = ts->m_vaCoords[m_laAngleAtomList[z*6+3]];
				vecb = ts->m_vaCoords[m_laAngleAtomList[z*6+4]];
				vec2 = FoldVector(vecb-veca);
			}
			m_fAngles[z] = Angle_Deg(vec1,vec2);
			m_bAnglePassed[z] = false;
			for (z2=0;z2<parent->m_pADF->m_faMinMaxAngle.GetSize();z2+=2)
			{
				if ((m_fAngles[z] >= parent->m_pADF->m_faMinMaxAngle[z2]) && (m_fAngles[z] <= parent->m_pADF->m_faMinMaxAngle[z2+1]))
				{
					m_bAnyAnglePassed = true;
					m_bAnglePassed[z] = true;
					break;
				}
			}
		}
	} else
	{
		m_bAnyAnglePassed = true;
		m_bAnglePassed[0] = true;
	}
}


CExtendedCondition::CExtendedCondition()
{
}


CExtendedCondition::~CExtendedCondition()
{
}


CNbSearch::CNbSearch()
{
	m_pRDF = NULL;
	m_pADF = NULL;
	m_iNbCountMin = -1;
	m_bInactive = false;
	m_pNbSort = NULL;
	m_pDistances = NULL;
	m_pAngles = NULL;
	m_oaNbPairs.SetName("CNbSearch::m_oaNbPairs");
	m_oaExtendedConditions.SetName("CNbSearch::m_oaExtendedConditions");
}


CNbSearch::~CNbSearch()
{
}


void CNbSearch::MarkPassedAtoms(int om)
{
	int /*z,*/ z2, z3, z4;
	CNbPair *n;

	n = (CNbPair*)m_oaNbPairs[om];

/*	mprintf("****A****\n");

	for (z=0;z<g_iGesVirtAtomCount;z++)
		mprintf("!    %d: %d\n",z+1,g_baAtomPassedCondition[z]);*/

	if (m_pRDF != NULL)
	{
		for (z2=0;z2<m_iDistances;z2++)
			for (z4=z2*2;z4<z2*2+2;z4++)
				if (g_baAtomPassedCondition[n->m_laDistAtomList[z4]] == 110)
					g_baAtomPassedCondition[n->m_laDistAtomList[z4]] = 0;
	}

	if (m_pADF != NULL)
	{
		for (z3=0;z3<m_iAngles;z3++)
		{
			if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6]] == 110)
				g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6]] = 0;
			if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+1]] == 110)
				g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+1]] = 0;
			if (m_pADF->m_bOrtho[0])
				if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+2]] == 110)
					g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+2]] = 0;
			if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+3]] == 110)
				g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+3]] = 0;
			if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+4]] == 110)
				g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+4]] = 0;
			if (m_pADF->m_bOrtho[1])
				if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+5]] == 110)
					g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+5]] = 0;
		}
	}

/*	mprintf("****B****\n");

	for (z=0;z<g_iGesVirtAtomCount;z++)
		mprintf("!    %d: %d\n",z+1,g_baAtomPassedCondition[z]);*/

	for (z2=0;z2<m_iDistances;z2++)
	{
		for (z3=0;z3<m_iAngles;z3++)
		{
			if (n->m_bDistPassed[z2] && n->m_bAnglePassed[z3])
			{
//				mprintf("@ Blubb 0\n");
				if (m_bCombinationMatrix[z2*m_iAngles+z3])
				{
//					mprintf("@ Blubb\n");
					if (m_pRDF)
					{
						for (z4=z2*2;z4<z2*2+2;z4++)
							if (g_baAtomPassedCondition[n->m_laDistAtomList[z4]] == 0)
								g_baAtomPassedCondition[n->m_laDistAtomList[z4]] = 1;
					}

					if (m_pADF)
					{
						if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6]] == 0)
							g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6]] = 1;

						if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+1]] == 0)
							g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+1]] = 1;

						if (m_pADF->m_bOrtho[0])
							if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+2]] == 0)
								g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+2]] = 1;

						if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+3]] == 0)
							g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+3]] = 1;

						if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+4]] == 0)
							g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+4]] = 1;

						if (m_pADF->m_bOrtho[1])
							if (g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+5]] == 0)
								g_baAtomPassedCondition[n->m_laAngleAtomList[z3*6+5]] = 1;
					}
				}
			}
		}
	}

/*	mprintf("****C****\n");

	for (z=0;z<g_iGesVirtAtomCount;z++)
		mprintf("!    %d: %d\n",z+1,g_baAtomPassedCondition[z]);*/
}


void CNbSearch::ScanAllOM(CSingleMolecule *rm, CTimeStep *t)
{
	int z, z2, z3, ti, ti3;
	CNbPair *n;
	bool anypassed;

	m_fMoleculesTotal += ((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();
	m_fCombinationsTotal += ((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize() * m_iDistances * m_iAngles;

	ti3 = 0;
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
	{
		n = (CNbPair*)m_oaNbPairs[z];
		m_bPassed[z] = false;
		anypassed = false;
		m_iCombPassCount[z] = 0;
		if ((rm->m_iMolSMIndex == z) && (rm->m_iMolType == m_iObsMol))
		{
			n->m_fMinDist = 9e20f;
			continue;
		}

		n->Check(this,t,rm,(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex[z]]);

//		mprintf("A\n");
		if (m_iNbCountMin <= -1) // Abstands-Modus
		{
//			mprintf("B\n");
			for (z2=0;z2<m_iDistances;z2++)
			{
				if (n->m_bDistPassed[z2])
				{
					anypassed = true;
					m_fCombinationsFound[(z2+1)*(m_iAngles+1)] += 1.0; //m_iAngles;
				}
			}

			for (z3=0;z3<m_iAngles;z3++)
			{
				if (n->m_bAnglePassed[z3])
				{
					anypassed = true;
					m_fCombinationsFound[z3+1] += 1.0; //m_iDistances;
				}
			}

/*			mprintf("$$$ davor\n");
			for (z2=0;z2<m_iDistances;z2++)
				mprintf("Distance %d: %d\n",z2+1,n->m_bDistPassed[z2]);

			for (z3=0;z3<m_iAngles;z3++)
				mprintf("Angle %d: %d\n",z3+1,n->m_bAnglePassed[z3]);*/

			for (z2=0;z2<m_iDistances;z2++)
			{
				for (z3=0;z3<m_iAngles;z3++)
				{
					if (n->m_bDistPassed[z2] && n->m_bAnglePassed[z3])
					{
//						mprintf("C\n");
						m_fCombinationsFound[(z2+1)*(m_iAngles+1)+z3+1]++;
						if (m_bCombinationMatrix[z2*m_iAngles+z3])
						{
//							mprintf("D\n");
							if (m_bExtendedMode)
							{
//								mprintf("E\n");
								if (CheckExtended(n->m_fDistances[z2],n->m_fAngles[z3]))
								{
		//							mprintf("%f | %f passed: Dist %d Angle %d\n",n->m_fDistances[z2],n->m_fAngles[z3],n->m_bDistPassed[z2],n->m_bAnglePassed[z3]);
									if (!m_bPassed[z])
									{
										ti3++;
										m_bPassed[z] = true;
										m_iPassCounter[z]++;
									}
									m_iCombPassCount[z]++;
									m_fCombinationsPassed++;
								} else
								{
									n->m_bDistPassed[z2] = false;
									n->m_bAnglePassed[z3] = false;
								}
							} else
							{
								if (!m_bPassed[z])
								{
									ti3++;
									m_bPassed[z] = true;
									m_iPassCounter[z]++;
								}
								m_iCombPassCount[z]++;
								m_fCombinationsPassed++;
							}
						}
					}
				}
			}

/*			mprintf("$$$ danach\n");
			for (z2=0;z2<m_iDistances;z2++)
				mprintf("Distance %d: %d\n",z2+1,n->m_bDistPassed[z2]);

			for (z3=0;z3<m_iAngles;z3++)
				mprintf("Angle %d: %d\n",z3+1,n->m_bAnglePassed[z3]);*/

			if (anypassed)
				m_fCombinationsFound[0] += 1.0; //m_iDistances * m_iAngles;
			if (m_bPassed[z])
				m_fMoleculesPassed++;
		} // Ende Abstandsmodus
	}
	if (m_iNbCountMin > -1) // Nachbar-Anzahl-Modus
	{
		SortNeighbors();

		for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
		{
			ti = m_pNbSort[z].m_iOM;
			if ((rm->m_iMolSMIndex == ti) && (rm->m_iMolType == m_iObsMol))
				continue;
			n = (CNbPair*)m_oaNbPairs[ti];
			if ((z >= m_iNbCountMin) && (z <= m_iNbCountMax))
			{
				n->m_bAnyDistPassed = true;
				n->m_bDistPassed[n->m_iMinDistIndex] = true;
				if (n->m_bAnyAnglePassed)
				{
					m_bPassed[ti] = true;
					m_iPassCounter[ti]++;
					m_fMoleculesPassed++;
					m_fCombinationsFound[0] += m_iDistances * m_iAngles;
				} else m_bPassed[ti] = false;
			} else
			{
				n->m_bAnyDistPassed = false;
				m_bPassed[ti] = false;
			}
		}
	}
}


void CNbSearch::PreScanAllOM(CSingleMolecule *rm, CTimeStep *t)
{
	int z;
	CNbPair *n;

	for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
	{
		n = (CNbPair*)m_oaNbPairs[z];

		if ((rm->m_iMolSMIndex == z) && (rm->m_iMolType == m_iObsMol))
			continue;

		n->PreCheck(this,t,rm,(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex[z]]);
	}

	if (m_iNbCountMin > -1) // Nachbar-Anzahl-Modus
		SortNeighbors();
}


void CNbSearch::Create(int obs)
{
	int z;
	CNbPair *n;

	try { m_bPassed = new bool[((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_bPassed = NULL; }
	if (m_bPassed == NULL) NewException((double)((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_iCombPassCount = new int[((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iCombPassCount = NULL; }
	if (m_iCombPassCount == NULL) NewException((double)((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_iPassCounter = new int[((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iPassCounter = NULL; }
	if (m_iPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	if (m_pADF != NULL)
		m_iAngles = m_pADF->m_iCombinations;
			else m_iAngles = 1;
	if (m_pRDF != NULL)
		m_iDistances = m_pRDF->m_iCombinations;
			else m_iDistances = 1;
	m_oaNbPairs.RemoveAll();

	try { m_pDistances = new float[m_iDistances*((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pDistances = NULL; }
	if (m_pDistances == NULL) NewException((double)m_iDistances*((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_pAngles = new float[m_iAngles*((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pAngles = NULL; }
	if (m_pAngles == NULL) NewException((double)m_iAngles*((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize()*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<((CMolecule*)g_oaMolecules[obs])->m_laSingleMolIndex.GetSize();z++)
	{
		m_iPassCounter[z] = 0;

		try { n = new CNbPair(); } catch(...) { n = NULL; }
		if (n == NULL) NewException((double)sizeof(CNbPair),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		n->m_fDistances = &m_pDistances[z*m_iDistances];
		n->m_fAngles = &m_pAngles[z*m_iAngles];
		n->Create(this);
		m_oaNbPairs.Add(n);
	}

	try { m_fCombinationsFound = new double[(m_iAngles+1)*(m_iDistances+1)]; } catch(...) { m_fCombinationsFound = NULL; }
	if (m_fCombinationsFound == NULL) NewException((double)(m_iAngles+1)*(m_iDistances+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<(m_iAngles+1)*(m_iDistances+1);z++)
		m_fCombinationsFound[z] = 0;

	try { m_bCombinationMatrix = new bool[m_iAngles*m_iDistances]; } catch(...) { m_bCombinationMatrix = NULL; }
	if (m_bCombinationMatrix == NULL) NewException((double)m_iAngles*m_iDistances*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_fCombinationsTotal = 0;
	m_fCombinationsPassed = 0;
	m_fMoleculesTotal = 0;
	m_fMoleculesPassed = 0;
}


void CNbSearch::Parse(int rm, int om, bool nbana)
{
	CxIntArray wa, tempwa;
	int z, z2, z3, ti, ti2, ti3;
	CxObArray *oa;
	CExtendedCondition *ec;
//	char buf[256];
	CxString buf;

	mprintf(YELLOW,"\n>>> Condition between %s and %s >>>\n\n",((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[om])->m_sName);
	m_iRefMol = rm;
	m_iObsMol = om;
	if (nbana)
		goto _rdf;
	m_bExtendedMode = AskYesNo("    Use set of linear functions as condition (\"extended mode\") (y/n)? [no] ",false);
	
	if (m_bExtendedMode)
	{
		g_bEnvDisableSortNb = true;

		mprintf(WHITE,"\n    For the extended mode you have to define a distance and an angle.\n\n");

		try { m_pRDF = new CRDF(); } catch(...) { m_pRDF = NULL; }
		if (m_pRDF == NULL) NewException((double)sizeof(CRDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pRDF->m_iShowMol = om;
		m_pRDF->ParseCondition(rm,this,true);

		try {m_pADF = new CADF(); } catch(...) { m_pADF = NULL; }
		if (m_pADF == NULL) NewException((double)sizeof(CADF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pADF->m_iShowMol = om;
		m_pADF->ParseCondition(rm,true);
		mprintf("\n    You can now define linear functions of distance and angle.\n");
		mprintf("    You can create several blocks of these functions.\n");
		mprintf("    The functions in each block are connected with \"AND\" (have to be fulfilled at the same time).\n");
		mprintf("    The blocks of functions are connected with \"OR\" (at least one of them has to be fulfilled).\n\n");
		mprintf("    For each linear function, you have to enter two distance-angle pairs, which define a line.\n");
		mprintf("    A third distange-angle pair controls which side of the line is to be included.\n");

		z = 0;
		z2 = 0;
		while (true)
		{
			mprintf(YELLOW,"\n*** %d. block of functions ***\n",z+1);

			try {	oa = new CxObArray("CNbSearch::Parse():oa"); } catch(...) { oa = NULL; }
			if (oa == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
			m_oaExtendedConditions.Add(oa);
			while (true)
			{
				try {	ec = new CExtendedCondition(); } catch(...) { ec = NULL; }
				if (ec == NULL) NewException((double)sizeof(CExtendedCondition),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				oa->Add(ec);
				mprintf(WHITE,"\n    * %d. function in %d. block\n\n",z2+1,z+1);
_lineagain:
				ec->m_fD[0] = AskFloat_ND("      Line point 1 - enter distance value in pm: ");
				ec->m_fA[0] = AskFloat_ND("      Line point 1 - enter angle value in degree: ");
				ec->m_fD[1] = AskFloat_ND("      Line point 2 - enter distance value in pm: ");
				ec->m_fA[1] = AskFloat_ND("      Line point 2 - enter angle value in degree: ");
				if ((ec->m_fD[1] == ec->m_fD[0]) && (ec->m_fA[1] == ec->m_fA[0]))
				{
					eprintf("\n      The two points need to be different.\n\n");
					goto _lineagain;
				}
_sampleagain:
				ec->m_fD[2] = AskFloat_ND("      Included sample point - enter distance value in pm: ");
				ec->m_fA[2] = AskFloat_ND("      Included sample point - enter angle value in degree: ");
				
				if (!ec->Evaluate())
				{
					eprintf("\n      Sample point may not be located on the line!\n\n");
					goto _sampleagain;
				}
				mprintf("\n      (resulting equation: %G * dist + %G * angle %c %G).\n",ec->m_fX,ec->m_fY,ec->m_bLarger?'>':'<',ec->m_fZ);

				z2++;
				mprintf("\n");
				if (!AskYesNo("    Add another function to this block (y/n)? [no] ",false))
					break;
			}
			z++;
			mprintf("\n");
			if (!AskYesNo("    Add another block of functions (y/n)? [no] ",false))
				break;
		}
		mprintf("\n");
	} else // END IF EXTENDED MODE
	{
_condagain:
		if (AskYesNo("    Do you want to define a distance condition (y/n)? [yes] ",true))
		{
_rdf:
			try {	m_pRDF = new CRDF(); } catch(...) { m_pRDF = NULL; }
			if (m_pRDF == NULL) NewException((double)sizeof(CRDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			m_pRDF->m_iShowMol = om;
			m_pRDF->ParseCondition(rm,this,nbana);
		}
		if (!nbana)
		{
			if (AskYesNo("    Do you want to define an angular condition (y/n)? [%s] ",(m_pRDF == NULL),(m_pRDF == NULL)?"yes":"no"))
			{
				g_bEnvDisableSortNb = true;

				try {	m_pADF = new CADF(); } catch(...) { m_pADF = NULL; }
				if (m_pADF == NULL) NewException((double)sizeof(CADF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				m_pADF->m_iShowMol = om;
				m_pADF->ParseCondition(rm,false);
			}
		}
		if ((m_pRDF == NULL) && (m_pADF == NULL))
		{
			eprintf("Error: You need to specify at least one condition.\n");
			goto _condagain;
		}
	}
	Create(om);
	if ((m_pADF != NULL) && (m_pRDF != NULL))
	{
		mprintf(WHITE,"\n    There are %d distance subconditions:\n",m_iDistances);
		m_pRDF->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[om])->m_laSingleMolIndex[0]],&wa);
		for (z=0;z<wa.GetSize()/2;z++)
		{
			ti = wa[z*2];
			ti2 = wa[z*2+1];
			mprintf("      * %2d.) Distance %s%d (%s) <--> %s%d (%s)\n",z+1,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
		}
		mprintf("\n");
		mprintf(WHITE,"    There are %d angular subconditions:\n",m_iAngles);
		m_pADF->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[om])->m_laSingleMolIndex[0]],&wa);
		for (z=0;z<wa.GetSize()/6;z++)
		{
			mprintf("      * %2d.) Angle ",z+1);
			if (m_pADF->m_bOrtho[0])
			{
				ti = wa[z*6];
				ti2 = wa[z*6+1];
				ti3 = wa[z*6+2];
				mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)] to ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
			} else
			{
				ti = wa[z*6];
				ti2 = wa[z*6+1];
				mprintf("[%s%d (%s) --> %s%d (%s)] to ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
			}
			if (m_pADF->m_bOrtho[1])
			{
				ti = wa[z*6+3];
				ti2 = wa[z*6+4];
				ti3 = wa[z*6+5];
				mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)]\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
			} else
			{
				ti = wa[z*6+3];
				ti2 = wa[z*6+4];
				mprintf("[%s%d (%s) --> %s%d (%s)]\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
			}
		}
		mprintf("\n");
		if (m_iNbCountMin > -1)
		{
			mprintf("    Min./max. next neighbor count specified, combining each with each subcondition.\n");
			for (z=0;z<m_iDistances;z++)
				for (z2=0;z2<m_iAngles;z2++)
					m_bCombinationMatrix[z*m_iAngles+z2] = true;
		} else
		{
			if ((m_iDistances > 1) && (m_iAngles > 1))
			{
				if (m_iDistances == m_iAngles)
				{
					if (AskYesNo("    Combine n-th distance subcondition with n-th angular subcondition (y/n)? [yes] ",true))
					{
						for (z=0;z<m_iDistances;z++)
							for (z2=0;z2<m_iAngles;z2++)
								if (z == z2)
									m_bCombinationMatrix[z*m_iAngles+z2] = true;
										else m_bCombinationMatrix[z*m_iAngles+z2] = false;
					} else goto _notn;
				} else
				{
_notn:
					if (AskYesNo("    Combine all distance subconditions with all angular subconditions (y/n)? [yes] ",true))
					{
						for (z=0;z<m_iDistances;z++)
							for (z2=0;z2<m_iAngles;z2++)
								m_bCombinationMatrix[z*m_iAngles+z2] = true;
					} else
					{
						for (z=0;z<m_iDistances;z++)
							for (z2=0;z2<m_iAngles;z2++)
									m_bCombinationMatrix[z*m_iAngles+z2] = false;

						mprintf("\n    Please enter all the combinations you want to use.\n");
						mprintf("    Enter each combination as a comma-separated 2-tuple (e.g. 1,2).\n\n");
_nextcomb:
						AskString("    Enter combination (return=finished): ",&buf,"");
						if (strlen(buf)==0)
							goto _combdone;
						tempwa.RemoveAll_KeepSize();
						ParseIntList(buf,&tempwa);
						if (tempwa.GetSize() != 2)
						{
							eprintf("    Wrong input, %d instead of 2 values.\n",tempwa.GetSize());
							goto _nextcomb;
						}
						if ((tempwa[0] < 1) || (tempwa[0] > m_iDistances))
						{
							eprintf("    Wrong input, only %d distance subconditions (%d requested).\n",m_iDistances,tempwa[0]);
							goto _nextcomb;
						}
						if ((tempwa[1] < 1) || (tempwa[1] > m_iAngles))
						{
							eprintf("    Wrong input, only %d angle subconditions (%d requested).\n",m_iAngles,tempwa[1]);
							goto _nextcomb;
						}
						tempwa[0]--;
						tempwa[1]--;
						if (m_bCombinationMatrix[tempwa[0]*m_iAngles+tempwa[1]])
						{
							eprintf("    This combination has already been added.\n");
							goto _nextcomb;
						}
						m_bCombinationMatrix[tempwa[0]*m_iAngles+tempwa[1]] = true;
						goto _nextcomb;
_combdone:;
					}
				}
				mprintf("\n");
			} else
			{
				for (z=0;z<m_iDistances;z++)
					for (z2=0;z2<m_iAngles;z2++)
						m_bCombinationMatrix[z*m_iAngles+z2] = true;
			}
		}
	} else
	{
		for (z=0;z<m_iDistances*m_iAngles;z++)
			m_bCombinationMatrix[z] = true;
	}
	m_iCombinationsEnabled = 0;
	for (z=0;z<m_iDistances;z++)
		for (z2=0;z2<m_iAngles;z2++)
			if (m_bCombinationMatrix[z*m_iAngles+z2])
				m_iCombinationsEnabled++;
	if (!nbana)
	{
		mprintf("\n    Using %d combinations of subconditions for this condition:\n\n",m_iCombinationsEnabled);
		z3 = 0;
		for (z=0;z<m_iDistances;z++)
		{
			for (z2=0;z2<m_iAngles;z2++)
			{
				if (m_bCombinationMatrix[z*m_iAngles+z2])
				{
					mprintf("      %2d.) Distance %2d - Angle %2d.\n",z3+1,z+1,z2+1);
					z3++;
				}
			}
		}
	}

	if (m_bExtendedMode)
	{
		mprintf(WHITE,"\n    You have the following extended condition functions:\n\n");
		for (z=0;z<m_oaExtendedConditions.GetSize();z++)
		{
			oa = (CxObArray*)m_oaExtendedConditions[z];
			for (z2=0;z2<oa->GetSize();z2++)
			{
				ec = (CExtendedCondition*)oa->GetAt(z2);
				mprintf("    %G * dist + %G * angle %c %G\n",ec->m_fX,ec->m_fY,ec->m_bLarger?'>':'<',ec->m_fZ);
				if (z2+1 < oa->GetSize())
					mprintf("      AND\n");
			}
			if (z+1 < m_oaExtendedConditions.GetSize())
				mprintf("\n      OR\n\n");
		}
	}
	mprintf(YELLOW,"\n<<< End of this condition <<<\n\n");
}


void CNbSearch::ParseGrid(int rm, int om, int gridmode)
{
	CxIntArray wa, tempwa;
	CxObArray *oa;
	CExtendedCondition *ec;
	int z, z2, ti, ti2, ti3;
//	char buf[256];
	CxString buf;

	mprintf(YELLOW,"\n>>> Condition between %s and %s >>>\n\n",((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[om])->m_sName);
	m_iRefMol = rm;
	m_iObsMol = om;
	if (gridmode == 6)
	{
		m_bExtendedMode = AskYesNo("    Use set of linear functions as condition (\"extended mode\") (y/n)? [no] ",false);

		if (m_bExtendedMode)
		{
			mprintf(WHITE,"\n    For the extended mode you have to define a distance and an angle.\n\n");

			try {	m_pRDF = new CRDF(); } catch(...) { m_pRDF = NULL; }
			if (m_pRDF == NULL) NewException((double)sizeof(CRDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			m_pRDF->m_iShowMol = om;
			m_pRDF->ParseCondition(rm,this,true);

			try {	m_pADF = new CADF(); } catch(...) { m_pADF = NULL; }
			if (m_pADF == NULL) NewException((double)sizeof(CADF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			m_pADF->m_iShowMol = om;
			m_pADF->ParseCondition(rm,true);
			mprintf("\n    You can now define linear functions of distance and angle.\n");
			mprintf("    You can create several blocks of these functions.\n");
			mprintf("    The functions in each block are connected with \"AND\" (have to be fulfilled at the same time).\n");
			mprintf("    The blocks of functions are connected with \"OR\" (at least one of them has to be fulfilled).\n\n");
			mprintf("    For each linear function, you have to enter two distance-angle pairs, which define a line.\n");
			mprintf("    A third distange-angle pair controls which side of the line is to be included.\n");

			z = 0;
			z2 = 0;
			while (true)
			{
				mprintf(YELLOW,"\n*** %d. block of functions ***\n",z+1);

				try {	oa = new CxObArray("CNbSearch::ParseGrid():oa"); } catch(...) { oa = NULL; }
				if (oa == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				m_oaExtendedConditions.Add(oa);
				while (true)
				{
					try {	ec = new CExtendedCondition(); } catch(...) { ec = NULL; }
					if (ec == NULL) NewException((double)sizeof(CExtendedCondition),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					oa->Add(ec);
					mprintf(WHITE,"\n    * %d. function in %d. block\n\n",z2+1,z+1);
	_lineagain:
					ec->m_fD[0] = AskFloat_ND("      Line point 1 - enter distance value in pm: ");
					ec->m_fA[0] = AskFloat_ND("      Line point 1 - enter angle value in degree: ");
					ec->m_fD[1] = AskFloat_ND("      Line point 2 - enter distance value in pm: ");
					ec->m_fA[1] = AskFloat_ND("      Line point 2 - enter angle value in degree: ");
					if ((ec->m_fD[1] == ec->m_fD[0]) && (ec->m_fA[1] == ec->m_fA[0]))
					{
						eprintf("\n      The two points need to be different.\n\n");
						goto _lineagain;
					}
	_sampleagain:
					ec->m_fD[2] = AskFloat_ND("      Included sample point - enter distance value in pm: ");
					ec->m_fA[2] = AskFloat_ND("      Included sample point - enter angle value in degree: ");
					
					if (!ec->Evaluate())
					{
						eprintf("\n      Sample point may not be located on the line!\n\n");
						goto _sampleagain;
					}
					mprintf("\n      (resulting equation: %G * dist + %G * angle %c %G).\n",ec->m_fX,ec->m_fY,ec->m_bLarger?'>':'<',ec->m_fZ);

					z2++;
					mprintf("\n");
					if (!AskYesNo("    Add another function to this block (y/n)? [no] ",false))
						break;
				}
				z++;
				mprintf("\n");
				if (!AskYesNo("    Add another block of functions (y/n)? [no] ",false))
					break;
			}
			mprintf("\n");
			goto _noadf;
		}
	} else m_bExtendedMode = false;

	if ((gridmode == 1) || (gridmode == 3) || (gridmode == 4) || (gridmode == 5) || (gridmode == 6))
	{
		if (gridmode == 6)
			if (!AskYesNo("    Do you want to define a distance condition (y/n)? [yes] ",true))
				goto _nordf;

		try {	m_pRDF = new CRDF(); } catch(...) { m_pRDF = NULL; }
		if (m_pRDF == NULL) NewException((double)sizeof(CRDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pRDF->m_iShowMol = om;
		m_pRDF->ParseConditionGrid(rm,this,gridmode);
_nordf:;
	}
	if ((gridmode == 2) || (gridmode == 3) || (gridmode == 4) || (gridmode == 5) || (gridmode == 6))
	{
		if ((gridmode == 4) || (gridmode == 5) || (gridmode == 6))
			if (!AskYesNo("    Do you want to define an angular condition (y/n)? [no] ",false))
				goto _noadf;

		try {	m_pADF = new CADF(); } catch(...) { m_pADF = NULL; }
		if (m_pADF == NULL) NewException((double)sizeof(CADF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pADF->m_iShowMol = om;
		m_pADF->ParseConditionGrid(rm,gridmode);
_noadf:;
	}
	Create(om);
	if ((m_pADF != NULL) && (m_pRDF != NULL))
	{
		mprintf(WHITE,"\nThere are %d distance subconditions:\n",m_iDistances);
		m_pRDF->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[om])->m_laSingleMolIndex[0]],&wa);
		for (z=0;z<wa.GetSize()/2;z++)
		{
			ti = wa[z*2];
			ti2 = wa[z*2+1];
			mprintf("  * %2d.) Distance %s%d (%s) <--> %s%d (%s)\n",z+1,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
		}
		mprintf("\n");
		mprintf(WHITE,"There are %d angular subconditions:\n",m_iDistances);
		m_pADF->BuildAtomList((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[rm])->m_laSingleMolIndex[0]],(CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[om])->m_laSingleMolIndex[0]],&wa);
		for (z=0;z<wa.GetSize()/6;z++)
		{
			mprintf("  * %2d.) Angle ",z+1);
			if (m_pADF->m_bOrtho[0])
			{
				ti = wa[z*6];
				ti2 = wa[z*6+1];
				ti3 = wa[z*6+2];
				mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)] to ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
			} else
			{
				ti = wa[z*6];
				ti2 = wa[z*6+1];
				mprintf("[%s%d (%s) --> %s%d (%s)] to ",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
			}
			if (m_pADF->m_bOrtho[1])
			{
				ti = wa[z*6+3];
				ti2 = wa[z*6+4];
				ti3 = wa[z*6+5];
				mprintf("[%s%d (%s), %s%d (%s), %s%d (%s)]\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti3]])->m_sName,g_waAtomMolNumber[ti3]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti3]])->m_sName);
			} else
			{
				ti = wa[z*6+3];
				ti2 = wa[z*6+4];
				mprintf("[%s%d (%s) --> %s%d (%s)]\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[ti2]])->m_sName,g_waAtomMolNumber[ti2]+1,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti2]])->m_sName);
			}
		}
		mprintf("\n");
		if (m_iNbCountMin != -1)
		{
			mprintf("    Min./max. next neighbor count specified, combining each with each subcondition.\n");
			for (z=0;z<m_iDistances;z++)
				for (z2=0;z2<m_iAngles;z2++)
					m_bCombinationMatrix[z*m_iAngles+z2] = true;
		} else
		{
			if ((m_iDistances > 1) && (m_iAngles > 1))
			{
				if (m_iDistances == m_iAngles)
				{
					if (AskYesNo("    Combine n-th distance subcondition with n-th angular subcondition (y/n)? [yes] ",true))
					{
						for (z=0;z<m_iDistances;z++)
							for (z2=0;z2<m_iAngles;z2++)
								if (z == z2)
									m_bCombinationMatrix[z*m_iAngles+z2] = true;
										else m_bCombinationMatrix[z*m_iAngles+z2] = false;
					} else goto _notn;
				} else
				{
_notn:
					if (AskYesNo("    Combine all distance subconditions with all angular subconditions (y/n)? [yes] ",true))
					{
						for (z=0;z<m_iDistances;z++)
							for (z2=0;z2<m_iAngles;z2++)
								m_bCombinationMatrix[z*m_iAngles+z2] = true;
					} else
					{
						for (z=0;z<m_iDistances;z++)
							for (z2=0;z2<m_iAngles;z2++)
									m_bCombinationMatrix[z*m_iAngles+z2] = false;

						mprintf("\n    Please enter all the combinations you want to use.\n");
						mprintf("    Enter each combination as a comma-separated 2-tuple (e.g. 1,2).\n\n");
_nextcomb:
						AskString("    Enter combination (return=finished): ",&buf,"");
						if (strlen(buf)==0)
							goto _combdone;
						tempwa.RemoveAll_KeepSize();
						ParseIntList(buf,&tempwa);
						if (tempwa.GetSize() != 2)
						{
							eprintf("    Wrong input, %d instead of 2 values.\n",tempwa.GetSize());
							goto _nextcomb;
						}
						if ((tempwa[0] < 1) || (tempwa[0] > m_iDistances))
						{
							eprintf("    Wrong input, only %d distance subconditions (%d requested).\n",m_iDistances,tempwa[0]);
							goto _nextcomb;
						}
						if ((tempwa[1] < 1) || (tempwa[1] > m_iAngles))
						{
							eprintf("    Wrong input, only %d angle subconditions (%d requested).\n",m_iAngles,tempwa[1]);
							goto _nextcomb;
						}
						tempwa[0]--;
						tempwa[1]--;
						if (m_bCombinationMatrix[tempwa[0]*m_iAngles+tempwa[1]])
						{
							eprintf("    This combination has already been added.\n");
							goto _nextcomb;
						}
						m_bCombinationMatrix[tempwa[0]*m_iAngles+tempwa[1]] = true;
						goto _nextcomb;
_combdone:;
					}
				}
				mprintf("\n");
			} else
			{
				for (z=0;z<m_iDistances;z++)
					for (z2=0;z2<m_iAngles;z2++)
						m_bCombinationMatrix[z*m_iAngles+z2] = true;
			}
		}
	} else
	{
		for (z=0;z<m_iDistances*m_iAngles;z++)
			m_bCombinationMatrix[z] = true;
	}
	m_iCombinationsEnabled = 0;
	for (z=0;z<m_iDistances;z++)
		for (z2=0;z2<m_iAngles;z2++)
			if (m_bCombinationMatrix[z*m_iAngles+z2])
				m_iCombinationsEnabled++;
	mprintf("\nUsing %d combinations of subconditions for this condition.\n",m_iCombinationsEnabled);

	if (m_bExtendedMode)
	{
		mprintf(WHITE,"\n    You have the following extended condition functions:\n\n");
		for (z=0;z<m_oaExtendedConditions.GetSize();z++)
		{
			oa = (CxObArray*)m_oaExtendedConditions[z];
			for (z2=0;z2<oa->GetSize();z2++)
			{
				ec = (CExtendedCondition*)oa->GetAt(z2);
				mprintf("    %G * dist + %G * angle %c %G\n",ec->m_fX,ec->m_fY,ec->m_bLarger?'>':'<',ec->m_fZ);
				if (z2+1 < oa->GetSize())
					mprintf("      AND\n");
			}
			if (z+1 < m_oaExtendedConditions.GetSize())
				mprintf("\n      OR\n\n");
		}
	}

	mprintf(YELLOW,"\n<<< End of this condition <<<\n\n");
}


void CNbSearch::Parse_OnlyValues()
{
//	mprintf(YELLOW,"\n>>> Condition between %s and %s >>>\n\n",((CMolecule*)g_oaMolecules[rm])->m_sName,((CMolecule*)g_oaMolecules[om])->m_sName);
	if (m_pRDF != NULL)
	{
		m_pRDF->ParseCondition_OnlyValues(this);
		if (m_pADF != NULL)
			mprintf("\n");
	}
	if (m_pADF != NULL)
		m_pADF->ParseCondition_OnlyValues();
}


void CNbSearch::PrintTable()
{
	int z, z2;

	mprintf(WHITE,"\n*** Condition %d ***\n\n",m_iNumber+1);
	mprintf("    %.0f molecules total, %.0f passed (%.4f percent).\n",m_fMoleculesTotal,m_fMoleculesPassed,m_fMoleculesPassed/m_fMoleculesTotal*100.0f);
	if (m_iNbCountMin == -1)
	{
		mprintf("    %.0f subconditions total, %.0f passed (%.4f percent).\n\n",m_fCombinationsTotal,m_fCombinationsPassed,m_fCombinationsPassed/m_fCombinationsTotal*100.0f);
		mprintf("    Frequency table for Neighborhood between %s and %s:\n",((CMolecule*)g_oaMolecules[m_iRefMol])->m_sName,((CMolecule*)g_oaMolecules[m_iObsMol])->m_sName);
		mprintf("      (enabled pairs of conditions are printed in color)\n\n");
		mprintf("          | Dist.  |");
		for (z=0;z<m_iDistances;z++)
			mprintf(" %4d   ",z+1);
		mprintf("\n");
		mprintf("    Angle | (all)  |");
		for (z=0;z<m_iDistances;z++)
			mprintf("%7.3f ",m_fCombinationsFound[(z+1)*(m_iAngles+1)]/m_fCombinationsFound[0]*100.0f);
		mprintf("\n");
		mprintf("    ------|--------|");
		for (z=0;z<m_iDistances;z++)
			mprintf("--------");
		mprintf("\n");
		for (z=0;z<m_iAngles;z++)
		{
			mprintf("      %3d |%7.3f |",z+1,m_fCombinationsFound[z+1]/m_fCombinationsFound[0]*100.0f);
			for (z2=0;z2<m_iDistances;z2++)
			{
				if (m_bCombinationMatrix[z2*m_iAngles+z])
					mprintf(GREEN,"%7.3f ",m_fCombinationsFound[(z2+1)*(m_iAngles+1)+z+1]/m_fCombinationsFound[0]*100.0f);
						else mprintf("%7.3f ",m_fCombinationsFound[(z2+1)*(m_iAngles+1)+z+1]/m_fCombinationsFound[0]*100.0f);
			}
			mprintf("\n");
		}
		mprintf("\n");
	}
}


void CNbSearch::PrintTable(FILE *a)
{
	int z, z2;

	mfprintf(a,"\n*** Condition %d ***\n\n",m_iNumber+1);
	mfprintf(a,"    %.0f molecules total, %.0f passed (%.4f percent).\n",m_fMoleculesTotal,m_fMoleculesPassed,m_fMoleculesPassed/m_fMoleculesTotal*100.0f);
	if (m_iNbCountMin == -1)
	{
		mfprintf(a,"    %.0f subconditions total, %.0f passed (%.4f percent).\n\n",m_fCombinationsTotal,m_fCombinationsPassed,m_fCombinationsPassed/m_fCombinationsTotal*100.0f);
		mfprintf(a,"    Frequency table for Neighborhood between %s and %s   \n\n",((CMolecule*)g_oaMolecules[m_iRefMol])->m_sName,((CMolecule*)g_oaMolecules[m_iObsMol])->m_sName);
		mfprintf(a,"          | Dist.  |");
		for (z=0;z<m_iDistances;z++)
			mfprintf(a," %4d   ",z+1);
		mfprintf(a,"\n");
		mfprintf(a,"    Angle | (all)  |");
		for (z=0;z<m_iDistances;z++)
			mfprintf(a," %7.3f ",m_fCombinationsFound[(z+1)*(m_iAngles+1)]/m_fCombinationsFound[0]*100.0f);
		mfprintf(a,"\n");
		mfprintf(a,"    ------|--------|");
		for (z=0;z<m_iDistances;z++)
			mfprintf(a,"--------");
		mfprintf(a,"\n");
		for (z=0;z<m_iAngles;z++)
		{
			mfprintf(a,"      %3d |%7.3f |",z+1,m_fCombinationsFound[z+1]/m_fCombinationsFound[0]*100.0f);
			for (z2=0;z2<m_iDistances;z2++)
				mfprintf(a,"%7.3f ",m_fCombinationsFound[(z2+1)*(m_iAngles+1)+z+1]/m_fCombinationsFound[0]*100.0f);
			mfprintf(a,"\n");
		}
		mfprintf(a,"\n");
	}
}


void CNbPair::Create(CNbSearch *parent)
{
//	m_fDistances = new float[parent->m_iDistances];
	try {	m_bDistPassed = new bool[parent->m_iDistances]; } catch(...) { m_bDistPassed = NULL; }
	if (m_bDistPassed == NULL) NewException((double)parent->m_iDistances*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
//	m_fAngles = new float[parent->m_iAngles];
	try {	m_bAnglePassed = new bool[parent->m_iAngles]; } catch(...) { m_bAnglePassed = NULL; }
	if (m_bAnglePassed == NULL) NewException((double)parent->m_iAngles*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
}


void CNbSearch::PrintSingle(int om)
{
	int z, z4;
	CNbPair *n;

	mprintf("*** Single Condition Dump ***\n");
	n = (CNbPair*)m_oaNbPairs[om];
	for (z=0;z<m_iDistances;z++)
	{
		mprintf("    Distance %d: %f (%s). Atoms: ",z+1,n->m_fDistances[z],n->m_bDistPassed[z]?"passed":"failed");
		for (z4=z*2;z4<z*2+2;z4++)
			mprintf("%d, ",n->m_laDistAtomList[z4]);
		mprintf("\n");
	}
	for (z=0;z<m_iAngles;z++)
	{
		mprintf("    Angle %d: %f (%s). Atoms: ",z+1,((CNbPair*)m_oaNbPairs[om])->m_fAngles[z],((CNbPair*)m_oaNbPairs[om])->m_bAnglePassed[z]?"passed":"failed");
		mprintf("%d, ",n->m_laAngleAtomList[z*6]);
		mprintf("%d, ",n->m_laAngleAtomList[z*6+1]);
		if (m_pADF->m_bOrtho[0])
			mprintf("%d, ",n->m_laAngleAtomList[z*6+2]);
		mprintf("%d, ",n->m_laAngleAtomList[z*6+3]);
		mprintf("%d, ",n->m_laAngleAtomList[z*6+4]);
		if (m_pADF->m_bOrtho[1])
			mprintf("%d, ",n->m_laAngleAtomList[z*6+5]);
		mprintf("\n");
	}
}


void CNbSet::Parse(int rm)
{
	int z;
//	CNbSearch *n;
	CConditionGroup *g;

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		if ((z == rm) && (((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize() == 1))
		{
			m_oaConditionGroups.Add(NULL);
			continue;
		}
		if (!AskYesNo("    Use neighboring %s molecules (y/n)? [yes] ",true,((CMolecule*)g_oaMolecules[z])->m_sName))
		{
			m_oaConditionGroups.Add(NULL);
			continue;
		}
		mprintf(WHITE,"\n>>> Conditions for neighboring %s molecules >>>\n\n",((CMolecule*)g_oaMolecules[z])->m_sName);

		try {	g = new CConditionGroup(); } catch(...) { g = NULL; }
		if (g == NULL) NewException((double)sizeof(CConditionGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g->Parse(rm,z);
		m_oaConditionGroups.Add(g);
		mprintf(WHITE,"\n<<< End of Conditions for neighboring %s molecules <<<\n\n",((CMolecule*)g_oaMolecules[z])->m_sName);
/*		n = new CNbSearch();
		n->Parse(rm,z);
		m_oaMolecules.Add(n);*/
	}
}


void CNbSet::Scan(CSingleMolecule *rm, CTimeStep *t)
{
	int z;

	for (z=0;z<m_oaConditionGroups.GetSize();z++)
		if (m_oaConditionGroups[z]!= NULL)
			if (!((CConditionGroup*)m_oaConditionGroups[z])->m_bInactive)
				((CConditionGroup*)m_oaConditionGroups[z])->ScanNeighborhoodAllOM(t,rm);
}


void CNbSearch::Reset()
{
	int z;

	for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
	{
		m_bPassed[z] = false;
		m_iPassCounter[z] = 0;
		((CNbPair*)m_oaNbPairs[z])->Reset();
	}
	for (z=0;z<(m_iAngles+1)*(m_iDistances+1);z++)
		m_fCombinationsFound[z] = 0;
	m_fCombinationsTotal = 0;
	m_fCombinationsPassed = 0;
	m_fMoleculesTotal = 0;
	m_fMoleculesPassed = 0;
}


void CNbSet::Reset()
{
	int z;

	for (z=0;z<m_oaConditionGroups.GetSize();z++)
	{
		if (m_oaConditionGroups[z] != NULL)
			((CConditionGroup*)m_oaConditionGroups[z])->Reset();
	}
}


void CNbPair::Reset()
{

}


void CNbSet::Dump()
{
	int z, z2;
	CConditionGroup *g;
	CMolecule *m;

	mprintf(WHITE,"Neighborhood Environment:\n\n");
	for (z=0;z<m_oaConditionGroups.GetSize();z++)
	{
		g = (CConditionGroup*)m_oaConditionGroups[z];
		if (g == NULL)
			continue;
/*		if (g->m_bInactive)
			continue;*/
		m = (CMolecule*)g_oaMolecules[z];
		mprintf("  * Molecule %s: ",m->m_sName);
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
			if (g->Contains(z2))
				mprintf("%d ",z2+1);
		mprintf("\n");
	}
}


void CNbSet::AddMolecule(int moltype, int mol)
{
	CConditionGroup *cg;
	int z;

	if (m_oaConditionGroups[moltype] == NULL)
	{
		try {	cg = new CConditionGroup(); } catch(...) { cg = NULL; }
		if (cg == NULL) NewException((double)sizeof(CConditionGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_oaConditionGroups[moltype] = cg;
		cg->m_bInactive = true;
		cg->m_iShowMol = moltype;

		try {	cg->m_bAlwaysTrue = new bool[((CMolecule*)g_oaMolecules[moltype])->m_laSingleMolIndex.GetSize()]; } catch(...) { cg->m_bAlwaysTrue = NULL; }
		if (cg->m_bAlwaysTrue == NULL) NewException((double)((CMolecule*)g_oaMolecules[moltype])->m_laSingleMolIndex.GetSize()*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try {	cg->m_iPassCounter = new long[((CMolecule*)g_oaMolecules[moltype])->m_laSingleMolIndex.GetSize()]; } catch(...) { cg->m_iPassCounter = NULL; }
		if (cg->m_iPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[moltype])->m_laSingleMolIndex.GetSize()*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<((CMolecule*)g_oaMolecules[moltype])->m_laSingleMolIndex.GetSize();z++)
		{
			cg->m_bAlwaysTrue[z] = false;
			cg->m_iPassCounter[z] = 0;
		}
	}
	((CConditionGroup*)m_oaConditionGroups[moltype])->m_bAlwaysTrue[mol] = true;
}


CNbAnalysis::CNbAnalysis()
{
	m_oaDF.SetName("CNbAnalysis::m_oaDF");
	m_oaNPF.SetName("CNbAnalysis::m_oaNPF");
}


CNbAnalysis::~CNbAnalysis()
{
}


void CNbAnalysis::Parse()
{
	int z;

	mprintf(WHITE,"\n>>> Neighborhood Analysis >>>\n\n");
	mprintf("    You now have to define a neighborhood criterion.\n");
	m_fBinEntries = 0;

	try {	m_pNbSearch = new CNbSearch(); } catch(...) { m_pNbSearch = NULL; }
	if (m_pNbSearch == NULL) NewException((double)sizeof(CNbSearch),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pNbSearch->Parse(g_iFixMol,m_iShowMol,true);
	if (g_iFixMol == m_iShowMol)
		m_iNbCount = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()-1;
			else m_iNbCount = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
/*	m_iMinNbCount = AskUnsignedInteger("    From which neighbor on to analyze? [1] ",1)-1;
	if (m_iShowMol == g_iFixMol)
		ti = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()-1;
			else ti = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
	m_iMaxNbCount = AskRangeInteger("    Up to which neighbor to analyze? [%d] ",m_iMinNbCount+1,ti,ti,ti)-1;
*/	m_fMinDist = AskFloat("    Enter minimal distance for neighborhood analysis (pm): [0] ",0);
	m_fMaxDist = AskFloat("    Enter maximal distance for neighborhood analysis (pm): [%d] ",(float)HalfBox()*2.0f,HalfBox()*2);
	m_iResolution = AskUnsignedInteger("    Enter the binning resolution: [200] ",200);
	BuildName();

	try {	m_pDist = new double[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pDist = NULL; }
	if (m_pDist == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try {	m_pDistMin = new double[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pDistMin = NULL; }
	if (m_pDistMin == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try {	m_pDistMax = new double[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pDistMax = NULL; }
	if (m_pDistMax == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try {	m_pDistAvg = new double[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pDistAvg = NULL; }
	if (m_pDistAvg == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try {	m_pDistCount = new double[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pDistCount = NULL; }
	if (m_pDistCount == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
	{
		m_pDistMin[z] = 1e20;
		m_pDistMax[z] = 0;
		m_pDistAvg[z] = 0;
		m_pDistCount[z] = 0;
	}
	mprintf(WHITE,"\n<<< End of Neighborhood Analysis <<<\n");
}


void CNbAnalysis::BuildName()
{
	BTIN;
//	char tmp[32768];
	CxString tmp;

//	sprintf(tmp,"%s_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	tmp.sprintf("%s_%s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
	BTOUT;
}


void CNbAnalysis::AnalyzeStep()
{
	int z, z2;
	float r;

/*	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
		m_pDist[z] = ((CNbPair*)m_pNbSearch->m_oaNbPairs[z])->m_fMinDist;

	qsort(m_pDist,((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(),sizeof(double),compare_double);
*/
/*	mprintf("Sorted: ");
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
		mprintf("%.2f  ",m_pDist[z]);
	mprintf("\n");*/

	m_pNbSearch->SortNeighbors();

	m_fBinEntries += m_iNbCount;

	for (z=0;z<((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
		m_pDist[z] = m_pNbSearch->m_pNbSort[z].m_fMinDist;

	for (z=0;z<m_iNbCount;z++)
	{
		((CDF*)m_oaDF[z])->AddToBin(m_pDist[z]);
		if (m_pDist[z] < m_pDistMin[z])
			m_pDistMin[z] = m_pDist[z];
		if (m_pDist[z] > m_pDistMax[z])
			m_pDistMax[z] = m_pDist[z];
		m_pDistAvg[z] += m_pDist[z];
		m_pDistCount[z]++;
	}

	z2 = 0;
	for (z=0;z<m_iResolution;z++)
	{
		r = (float)(m_fMinDist+(z+0.5)*(m_fMaxDist-m_fMinDist)/m_iResolution);
		while ((m_pDist[z2] <= r) && (z2 < ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize()))
			z2++;
//		z2--;
//		mprintf("z=%d, r=%f, z2=%d\n",z,r,z2);
		((CDF*)m_oaNPF[z2])->AddToBin_Int(z);
		m_pNPFCount->AddToBin_Int(z);
	}
}


int compare_nbsort(const void *arg1, const void *arg2)
{
	if (((CNbSort*)arg1)->m_bAnyAnglePassed && (!((CNbSort*)arg2)->m_bAnyAnglePassed))
		return -1;
	if ((!((CNbSort*)arg1)->m_bAnyAnglePassed) && ((CNbSort*)arg2)->m_bAnyAnglePassed)
		return 1;

	if (((CNbSort*)arg1)->m_fMinDist > ((CNbSort*)arg2)->m_fMinDist)
		return 1;
	else if (((CNbSort*)arg1)->m_fMinDist < ((CNbSort*)arg2)->m_fMinDist)
		return -1;
	else return 0;
}


void CNbSearch::SortNeighbors()
{
	int z;
	CNbPair *n;

	if (m_pNbSort == NULL)
	{
		try {	m_pNbSort = new CNbSort[((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pNbSort = NULL; }
		if (m_pNbSort == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()*sizeof(CNbSort),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}

	for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
	{
		n = (CNbPair*)m_oaNbPairs[z];
		m_pNbSort[z].m_fMinDist = n->m_fMinDist;
		m_pNbSort[z].m_iOM = z;
		m_pNbSort[z].m_bAnyAnglePassed = n->m_bAnyAnglePassed;
	}
	
	qsort(m_pNbSort,((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize(),sizeof(CNbSort),compare_nbsort);

	for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
		((CNbPair*)m_oaNbPairs[m_pNbSort[z].m_iOM])->m_iNbPosition = z;
}


void CNbSet::ResetAlwaysTrue()
{
	CConditionGroup *cg;
	int z, z2;

	for (z=0;z<m_oaConditionGroups.GetSize();z++)
	{
		if (m_oaConditionGroups[z] == NULL)
			continue;
		cg = (CConditionGroup*)m_oaConditionGroups[z];
		for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
		{
			cg->m_bAlwaysTrue[z2] = false;
			cg->m_iPassCounter[z2] = 0;
		}
	}
}


void CNbSearch::CopyFrom(CNbSearch *nb)
{
	int z;
	CNbPair *np;
	CExtendedCondition *ec;

	m_iRefMol = nb->m_iRefMol;
	m_iObsMol = nb->m_iObsMol;
	m_bInactive = nb->m_bInactive;
	m_fCombinationsTotal = nb->m_fCombinationsTotal;
	m_fCombinationsPassed = nb->m_fCombinationsPassed;
	m_fMoleculesTotal = nb->m_fMoleculesTotal;
	m_fMoleculesPassed = nb->m_fMoleculesPassed;
	m_iAngles = nb->m_iAngles;
	m_iDistances = nb->m_iDistances;
	m_iCombinationsEnabled = nb->m_iCombinationsEnabled;
	m_iNbCountMin = nb->m_iNbCountMin;
	m_iNbCountMax = nb->m_iNbCountMax;
	m_iNumber = nb->m_iNumber;

	m_bExtendedMode = nb->m_bExtendedMode;
	m_oaExtendedConditions.RemoveAll_KeepSize();
	for (z=0;z<nb->m_oaExtendedConditions.GetSize();z++)
	{
		try {	ec = new CExtendedCondition(); } catch(...) { ec = NULL; }
		if (ec == NULL) NewException((double)sizeof(CExtendedCondition),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ec->CopyFrom((CExtendedCondition*)nb->m_oaExtendedConditions[z]);
		m_oaExtendedConditions.Add(ec);
	}

	if (nb->m_pRDF != NULL)
	{
		try {	m_pRDF = new CRDF(); } catch(...) { m_pRDF = NULL; }
		if (m_pRDF == NULL) NewException((double)sizeof(CRDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pRDF->CopyFrom(nb->m_pRDF);
	}
	if (nb->m_pADF != NULL)
	{
		try {	m_pADF = new CADF(); } catch(...) { m_pADF = NULL; }
		if (m_pADF == NULL) NewException((double)sizeof(CADF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pADF->CopyFrom(nb->m_pADF);
	}

	if (nb->m_bPassed != NULL)
	{
		try {	m_bPassed = new bool[((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_bPassed = NULL; }
		if (m_bPassed == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_bPassed,nb->m_bPassed,sizeof(bool)*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize());
	}
	if (nb->m_iCombPassCount != NULL)
	{
		try {	m_iCombPassCount = new int[((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iCombPassCount = NULL; }
		if (m_iCombPassCount == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_iCombPassCount,nb->m_iCombPassCount,sizeof(int)*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize());
	}
	if (nb->m_iPassCounter != NULL)
	{
		try {	m_iPassCounter = new int[((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_iPassCounter = NULL; }
		if (m_iPassCounter == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_iPassCounter,nb->m_iPassCounter,sizeof(int)*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize());
	}
	if (nb->m_fCombinationsFound != NULL)
	{
		try {	m_fCombinationsFound = new double[(m_iAngles+1)*(m_iDistances+1)]; } catch(...) { m_fCombinationsFound = NULL; }
		if (m_fCombinationsFound == NULL) NewException((double)(m_iAngles+1)*(m_iDistances+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_fCombinationsFound,nb->m_fCombinationsFound,sizeof(double)*(m_iAngles+1)*(m_iDistances+1));
	}
	if (nb->m_bCombinationMatrix != NULL)
	{
		try {	m_bCombinationMatrix = new bool[m_iAngles*m_iDistances]; } catch(...) { m_bCombinationMatrix = NULL; }
		if (m_bCombinationMatrix == NULL) NewException((double)m_iAngles*m_iDistances*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_bCombinationMatrix,nb->m_bCombinationMatrix,sizeof(bool)*m_iAngles*m_iDistances);
	}
	if (nb->m_pNbSort != NULL)
	{
		try {	m_pNbSort = new CNbSort[((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pNbSort = NULL; }
		if (m_pNbSort == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()*sizeof(CNbSort),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
			m_pNbSort[z] = nb->m_pNbSort[z];
	}

	if (nb->m_pDistances != NULL)
	{
		try {	m_pDistances = new float[m_iDistances*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pDistances = NULL; }
		if (m_pDistances == NULL) NewException((double)m_iDistances*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_pDistances,nb->m_pDistances,sizeof(float)*m_iDistances*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize());
	}

	if (nb->m_pAngles != NULL)
	{
		try {	m_pAngles = new float[m_iAngles*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pAngles = NULL; }
		if (m_pAngles == NULL) NewException((double)m_iAngles*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_pAngles,nb->m_pAngles,sizeof(float)*m_iAngles*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize());
	}
	
	m_oaNbPairs.RemoveAll_KeepSize();
	for (z=0;z<nb->m_oaNbPairs.GetSize();z++)
	{
		try {	np = new CNbPair(); } catch(...) { np = NULL; }
		if (np == NULL) NewException((double)sizeof(CNbPair),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		np->m_fDistances = &m_pDistances[z*m_iDistances];
		np->m_fAngles = &m_pAngles[z*m_iAngles];
		np->CopyFrom((CNbPair*)nb->m_oaNbPairs[z],nb);
		m_oaNbPairs.Add(np);
	}
}


void CNbPair::CopyFrom(CNbPair *p, CNbSearch *parent)
{
	m_iNbPosition = p->m_iNbPosition;
	m_bAnglePassed = p->m_bAnglePassed;
	m_iNbPosition = p->m_iNbPosition;
	m_bAnyDistPassed = p->m_bAnyDistPassed;
	m_bAnyAnglePassed = p->m_bAnyAnglePassed;
	m_fMinDist = p->m_fMinDist;
	m_iMinDistIndex = p->m_iMinDistIndex;

/*	if (p->m_fDistances != NULL)
	{
		m_fDistances = new float[parent->m_iDistances];
		memcpy(m_fDistances,p->m_fDistances,sizeof(float)*parent->m_iDistances);
	}*/
	if (p->m_bDistPassed != NULL)
	{
		try {	m_bDistPassed = new bool[parent->m_iDistances]; } catch(...) { m_bDistPassed = NULL; }
		if (m_bDistPassed == NULL) NewException((double)parent->m_iDistances*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_bDistPassed,p->m_bDistPassed,sizeof(bool)*parent->m_iDistances);
	}
/*	if (p->m_fAngles != NULL)
	{
		m_fAngles = new float[parent->m_iAngles];
		memcpy(m_fAngles,p->m_fAngles,sizeof(float)*parent->m_iAngles);
	}*/
	if (p->m_bAnglePassed != NULL)
	{
		try {m_bAnglePassed = new bool[parent->m_iAngles]; } catch(...) { m_bAnglePassed = NULL; }
		if (m_bAnglePassed == NULL) NewException((double)parent->m_iAngles*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_bAnglePassed,p->m_bAnglePassed,sizeof(bool)*parent->m_iAngles);
	}

	m_laDistAtomList.CopyFrom(&p->m_laDistAtomList);
	m_laAngleAtomList.CopyFrom(&p->m_laAngleAtomList);
}


void CNbSearch::CopyResults(CNbSearch *p)
{
/*	int z, z2*/;

/*	for (z=0;z<m_oaNbPairs.GetSize();z++)
	{
		memcpy(((CNbPair*)m_oaNbPairs[z])->m_fDistances,((CNbPair*)p->m_oaNbPairs[z])->m_fDistances,sizeof(float)*m_iDistances);
		memcpy(((CNbPair*)m_oaNbPairs[z])->m_fAngles,((CNbPair*)p->m_oaNbPairs[z])->m_fAngles,sizeof(float)*m_iAngles);
	}*/
	if (m_iNbCountMin > -1)
	{
		if (m_pNbSort == NULL)
		{
			try { m_pNbSort = new CNbSort[((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { m_pNbSort = NULL; }
			if (m_pNbSort == NULL) NewException((double)((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize()*sizeof(CNbSort),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}
/*		for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
			m_pNbSort[z] = p->m_pNbSort[z];*/
		memcpy(m_pNbSort,p->m_pNbSort,sizeof(CNbSort)*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize());
	} else
	{
		memcpy(m_pDistances,p->m_pDistances,sizeof(float)*m_iDistances*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize());
		memcpy(m_pAngles,p->m_pAngles,sizeof(float)*m_iAngles*((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize());
	}
}


void CNbSearch::ReScan(CSingleMolecule *rm)
{
	int z, z2, z3, ti, ti3;
	CNbPair *n;
	bool anypassed;

	m_fMoleculesTotal += ((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();
	m_fCombinationsTotal += ((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize() * m_iDistances * m_iAngles;

	ti3 = 0;
	for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
	{
		n = (CNbPair*)m_oaNbPairs[z];
		m_bPassed[z] = false;
		anypassed = false;
		m_iCombPassCount[z] = 0;
		if ((rm->m_iMolSMIndex == z) && (rm->m_iMolType == m_iObsMol))
		{
			n->m_fMinDist = 9e20f;
			continue;
		}

		n->ReScan(this);

		if (m_iNbCountMin <= -1) // Abstands-Modus
		{
			for (z2=0;z2<m_iDistances;z2++)
			{
				if (n->m_bDistPassed[z2])
				{
					anypassed = true;
					m_fCombinationsFound[(z2+1)*(m_iAngles+1)] += m_iAngles;
				}
			}

			for (z3=0;z3<m_iAngles;z3++)
			{
				if (n->m_bAnglePassed[z3])
				{
					anypassed = true;
					m_fCombinationsFound[z3+1] += m_iDistances;
				}
			}

			for (z2=0;z2<m_iDistances;z2++)
			{
				for (z3=0;z3<m_iAngles;z3++)
				{
					if (n->m_bDistPassed[z2] && n->m_bAnglePassed[z3])
					{
						m_fCombinationsFound[(z2+1)*(m_iAngles+1)+z3+1]++;
						if (m_bCombinationMatrix[z2*m_iAngles+z3])
						{
							if (m_bExtendedMode)
							{
								if (CheckExtended(n->m_fDistances[z2],n->m_fAngles[z3]))
								{
									if (!m_bPassed[z])
									{
										ti3++;
										m_bPassed[z] = true;
										m_iPassCounter[z]++;
									}
									m_iCombPassCount[z]++;
									m_fCombinationsPassed++;
								} else
								{
									n->m_bDistPassed[z2] = false;
									n->m_bAnglePassed[z3] = false;
								}
							} else
							{
								if (!m_bPassed[z])
								{
									ti3++;
									m_bPassed[z] = true;
									m_iPassCounter[z]++;
								}
								m_iCombPassCount[z]++;
								m_fCombinationsPassed++;
							}
						}
					}
				}
			}

			if (anypassed)
				m_fCombinationsFound[0] += m_iDistances * m_iAngles;
			if (m_bPassed[z])
				m_fMoleculesPassed++;
		} // Ende Abstandsmodus
	}
	if (m_iNbCountMin > -1) // Nachbar-Anzahl-Modus
	{
		for (z=0;z<((CMolecule*)g_oaMolecules[m_iObsMol])->m_laSingleMolIndex.GetSize();z++)
		{
			ti = m_pNbSort[z].m_iOM;
			if ((rm->m_iMolSMIndex == ti) && (rm->m_iMolType == m_iObsMol))
				continue;
			n = (CNbPair*)m_oaNbPairs[ti];
			if ((z >= m_iNbCountMin) && (z <= m_iNbCountMax))
			{
				n->m_bAnyDistPassed = true;
				n->m_bDistPassed[n->m_iMinDistIndex] = true;
				if (n->m_bAnyAnglePassed)
				{
					m_bPassed[ti] = true;
					m_iPassCounter[ti]++;
					m_fMoleculesPassed++;
					m_fCombinationsFound[0] += m_iDistances * m_iAngles;
				} else m_bPassed[ti] = false;
			} else
			{
				n->m_bAnyDistPassed = false;
				m_bPassed[ti] = false;
			}
		}
	}
}


void CNbPair::ReScan(CNbSearch *parent)
{
	int z, z2;

	if (parent->m_pRDF != NULL)
	{
		m_bAnyDistPassed = false;
		m_fMinDist = 1E30f;
		for (z=0;z<parent->m_iDistances;z++)
		{
			if (m_fDistances[z] < m_fMinDist)
			{
				m_fMinDist = m_fDistances[z];
				m_iMinDistIndex = z;
			}

			m_bDistPassed[z] = false;
			for (z2=0;z2<parent->m_pRDF->m_faMinMaxDist.GetSize();z2+=2)
			{
				if ((m_fDistances[z] >= parent->m_pRDF->m_faMinMaxDist[z2]) && (m_fDistances[z] <= parent->m_pRDF->m_faMinMaxDist[z2+1]))
				{
					m_bAnyDistPassed = true;
					m_bDistPassed[z] = true;
					break;
				} 
			}
		}
	} else
	{
		m_bAnyDistPassed = true;
		m_bDistPassed[0] = true;
	}

	if (parent->m_pADF != NULL)
	{
		m_bAnyAnglePassed = false;
		for (z=0;z<parent->m_iAngles;z++)
		{
			m_bAnglePassed[z] = false;
			for (z2=0;z2<parent->m_pADF->m_faMinMaxAngle.GetSize();z2+=2)
			{
				if ((m_fAngles[z] >= parent->m_pADF->m_faMinMaxAngle[z2]) && (m_fAngles[z] <= parent->m_pADF->m_faMinMaxAngle[z2+1]))
				{
					m_bAnyAnglePassed = true;
					m_bAnglePassed[z] = true;
					break;
				}
			}
		}
	} else
	{
		m_bAnyAnglePassed = true;
		m_bAnglePassed[0] = true;
	}
}


bool CNbSearch::CheckExtended(double dist, double angle)
{
	CxObArray *oa;
	CExtendedCondition *ec;
	int z, z2;

//	mprintf("\nChecking Dist=%G, Angle=%G...\n",dist,angle);
	for (z=0;z<m_oaExtendedConditions.GetSize();z++)
	{
		oa = (CxObArray*)m_oaExtendedConditions[z];
		for (z2=0;z2<oa->GetSize();z2++)
		{
			ec = (CExtendedCondition*)(*oa)[z2];
//			mprintf("   Function %d: %G * %G + %G * %G %c %G?  ",z2+1,ec->m_fX,dist,ec->m_fY,angle,ec->m_bLarger?'>':'<',ec->m_fZ);
			if (ec->m_bLarger)
			{
				if (ec->m_fX*dist + ec->m_fY*angle <= ec->m_fZ)
					goto _fail;
//				mprintf("Yes!\n");
			} else
			{
				if (ec->m_fX*dist + ec->m_fY*angle >= ec->m_fZ)
					goto _fail;
//				mprintf("Yes!\n");
			}
		}
//		mprintf("<return true>\n");
//		mprintf("\nBestanden! dist=%G, angle=%G.\n",dist,angle);
		return true;
_fail:;
//	  mprintf("No!\n");
	}
//	mprintf("<return false>\n");
	return false;
}


bool CExtendedCondition::Evaluate()
{
	double t;

	if ((m_fA[0]*m_fD[1] - m_fA[1]*m_fD[0]) == 0)
	{
//		mprintf("A\n");
		m_fZ = 0;
		if (m_fA[1] == m_fA[0])
		{
//			mprintf("X\n");
			m_fX = 0;
			m_fY = 1.0;
		} else
		{
//			mprintf("Y\n");
			m_fX = 1.0;
			m_fY = - (m_fD[1] - m_fD[0]) / (m_fA[1] - m_fA[0]);
		}
	} else
	{
//		mprintf("B\n");
		m_fX = -(m_fA[1] - m_fA[0]) / (m_fA[0]*m_fD[1] - m_fA[1]*m_fD[0]);
		m_fY = (m_fD[1] - m_fD[0]) / (m_fA[0]*m_fD[1] - m_fA[1]*m_fD[0]);
		m_fZ = 1.0;
	}

	t = m_fX * m_fD[2] + m_fY * m_fA[2];

	if (t == m_fZ)
		return false;

	if (t > m_fZ)
		m_bLarger = true;
			else m_bLarger = false;
	return true;
}


void CExtendedCondition::CopyFrom(CExtendedCondition *ec)
{
	int z;

	m_bLarger = ec->m_bLarger;
	m_fX = ec->m_fX;
	m_fY = ec->m_fY;
	m_fZ = ec->m_fZ;
	for (z=0;z<3;z++)
	{
		m_fA[z] = ec->m_fA[z];
		m_fD[z] = ec->m_fD[z];
	}
}


CNbSet::CNbSet()
{
	m_oaConditionGroups.SetName("CNBSet::m_oaConditionGroups");
}


CNbSet::~CNbSet()
{
}


