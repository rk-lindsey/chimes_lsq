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

#include "vorowrapper.h"
#include "voro++.h"
#include "tools.h"
#include "globalvar.h"
#include "maintools.h"


CVoroWrapper::CVoroWrapper()
{
	m_pContainer = NULL;

	m_oaVoroAtoms.SetName("CVoroWrapper::m_oaVoroAtoms");
	m_oaVoroMolecules.SetName("CVoroWrapper::m_oaVoroMolecules");
	m_faPOVFaceColor.SetName("CVoroWrapper::m_faPOVFaceColor");
	m_waPOVFaceColorMol.SetName("CVoroWrapper::m_waPOVFaceColorMol");
	m_waPOVFaceColorElem.SetName("CVoroWrapper::m_waPOVFaceColorElem");
	m_waPOVFaceColorAtom.SetName("CVoroWrapper::m_waPOVFaceColorAtom");
	m_waSurfCoverSM.SetName("CVoroWrapper::m_waSurfCoverSM");
	m_oaSurfCoverData.SetName("CVoroWrapper::m_oaSurfCoverData");
	m_waSurfCoverMol.SetName("CVoroWrapper::m_waSurfCoverMol");
	m_waSurfCoverElem.SetName("CVoroWrapper::m_waSurfCoverElem");
	m_waSurfCoverAtom.SetName("CVoroWrapper::m_waSurfCoverAtom");

	m_fPOVAngle = 90.0f;
	m_sPOVText[0] = 0;
	m_iPOVFrameCounter = 0;
	m_fPOVZoom = 1.0f;
	m_fPOVBoxOuterOpacity = 1.0f;
	m_vPOVBoxClipLow[0] = -1.0f;
	m_vPOVBoxClipLow[1] = -1.0f;
	m_vPOVBoxClipLow[2] = -1.0f;
	m_vPOVBoxClipHigh[0] = 1.0f;
	m_vPOVBoxClipHigh[1] = 1.0f;
	m_vPOVBoxClipHigh[2] = 1.0f;
	m_fPOVAtomCenter = 0.0f;
	m_fPOVScale = 0.95f;
	m_bPOVMolecules = false;
	m_fPOVExplode = 1.0f;
	m_fPOVClip = 1.0f;

	m_fPOVRayTrans = 0.0f;
	m_fPOVNbTrans = 0.0f;
	m_fPOVFaceTrans = 0.0f;
	m_fPOVEdgeTrans = 0.0f;
	m_fPOVNbHighlightFac = 1.3f;
	m_fPOVFaceBleach = 0.5f;
	m_vPOVFaceBleach[0] = 0.5f;
	m_vPOVFaceBleach[1] = 0.5f;
	m_vPOVFaceBleach[2] = 1.0f;

	m_bIncludeWritten = false;

	m_bVoroStat = false;
	m_bSurfCover = false;
	m_bWritePOV = false;

	strcpy(m_sPOVText,"ABCabc");
}


CVoroWrapper::~CVoroWrapper()
{
}


void CVoroWrapper::Build(CTimeStep *ts)
{
	int ijk, q, z, z2, z3, id, i, faces, co;
	CVoroAtom *va;
	CVoroMolecule *vm;
	double /**pp,*/ ta, surf, volume, vg;
	voronoicell_neighbor c;
	vector<int> nb;
	vector<int> fo;
	vector<double> fa;
	int *nbtemp;

	try { nbtemp = new int[m_oaVoroAtoms.GetSize()]; } catch(...) { nbtemp = NULL; }
	if (nbtemp == NULL) NewException((double)m_oaVoroAtoms.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_pContainer = new container_periodic_poly(g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0,m_iBlocksX,m_iBlocksY,m_iBlocksZ,g_iVoroMemory); } catch(...) { m_pContainer = NULL; }
	if (m_pContainer == NULL) NewException((double)sizeof(container_periodic_poly),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<g_iGesAtomCount;z++)
		m_pContainer->put(z,ts->m_vaCoords[z][0]/1000.0,ts->m_vaCoords[z][1]/1000.0,ts->m_vaCoords[z][2]/1000.0,g_faVoronoiRadii[z]/1000.0);

	c_loop_all_periodic vl(*m_pContainer);

	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		vm->tempfaces = 0;
		vm->tempsurf = 0;
		vm->tempvol = 0;
		for (z2=0;z2<m_oaVoroMolecules.GetSize();z2++)
			vm->m_pNbhTempMol[z2] = 0;
		for (z2=0;z2<m_oaVoroMolecules.GetSize();z2++)
			vm->m_pNbhTempArea[z2] = 0.0;
		for (z2=0;z2<g_iGesAtomCount;z2++)
			vm->m_pNbhTempAtoms[z2] = 0;
	}

	for (z=0;z<g_iGesAtomCount;z++)
		m_pAtomTouched[z] = false;

	co = 0;
	vg = 0;
	if (vl.start()) 
	{
		do 
		{
			if (m_pContainer->compute_cell(c,vl))
			{
				co++;

				ijk=vl.ijk;
				q=vl.q;
	//			pp=m_pContainer->p[ijk]+m_pContainer->ps*q;

				id = m_pContainer->id[ijk][q];
				c.face_areas(fa);
				c.face_orders(fo);
				c.neighbors(nb);

				m_pAtomTouched[id] = true;

				va = (CVoroAtom*)m_oaVoroAtoms[m_pAssignAtoms[id]];
				vm = (CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[id]];

				surf = c.surface_area();
				faces = c.number_of_faces();
				volume = c.volume();

				vg += volume;

				va->m_pVolume->AddToBin(volume*1000.0);
				va->m_pSurfaceArea->AddToBin(surf*100.0);
				va->m_pFaces->AddToBin_Int(faces);
				va->m_pMaxRadius->AddToBin(sqrt(c.max_radius_squared())*1000.0);
				va->m_pAVRatio->AddToBin(10.6347231054330961*volume/pow(surf,1.5));

				vm->tempvol += volume*1000.0;

				for (z2=0;z2<m_oaVoroAtoms.GetSize();z2++)
				{
					nbtemp[z2] = 0;
				}

				for (z2=0;z2<m_oaVoroMolecules.GetSize();z2++)
				{
					va->m_pNbhTempMol[z2] = 0;
				}

				ta = 0;
				for (z2=0;z2<faces;z2++)
				{
					va->m_pFaceOrders->AddToBin_Int(fo[z2]);

					if (m_pAssignMolecules[id] != m_pAssignMolecules[nb[z2]])
					{
						nbtemp[m_pAssignAtoms[nb[z2]]]++;
						va->m_pNbhDistAtoms[m_pAssignAtoms[nb[z2]]]->AddToBin(ts->FoldedDistance(id,nb[z2]));

						((CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[nb[z2]]])->m_pNbhDistAtoms[m_pAssignAtoms[id]]->AddToBin(ts->FoldedDistance(id,((CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[nb[z2]]])->m_iCenterOffset));

						if (vm->m_pNbhTempAtoms[nb[z2]] == 0)
						{
							vm->m_pNbhTempAtoms[nb[z2]]++;
							vm->m_pNbhDistAtoms[m_pAssignAtoms[nb[z2]]]->AddToBin(ts->FoldedDistance(id,((CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[nb[z2]]])->m_iCenterOffset));
						}

						if (vm->m_pNbhTempMol[m_pAssignMolecules[nb[z2]]] == 0)
						{
							vm->m_pNbhTempMol[m_pAssignMolecules[nb[z2]]]++;
							vm->m_pNbhDistMolecules[m_pAssignMoleculeTypes[nb[z2]]]->AddToBin(ts->FoldedDistance(((CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[id]])->m_iCenterOffset,((CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[nb[z2]]])->m_iCenterOffset));
						}
						vm->m_pNbhTempArea[m_pAssignMoleculeTypes[nb[z2]]] += fa[z2] * 100.0;

						if (va->m_pNbhTempMol[m_pAssignMolecules[nb[z2]]] == 0)
							va->m_pNbhDistMolecules[m_pAssignMoleculeTypes[nb[z2]]]->AddToBin(ts->FoldedDistance(id,((CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[nb[z2]]])->m_iCenterOffset));

						va->m_pNbhTempMol[m_pAssignMolecules[nb[z2]]]++;
						vm->tempsurf += fa[z2]*100.0;
						vm->tempfaces++;
						ta += fa[z2]*100.0;
					}
				}

				va->m_pExposedSurface->AddToBin(ta);
				va->m_pExposedSurfacePerc->AddToBin(ta/surf);

				for (z2=0;z2<m_oaVoroAtoms.GetSize();z2++)
				{
					va->m_pNbhAtoms[z2]->AddToBin_Int(nbtemp[z2]);
				}

				for (z2=0;z2<g_oaMolecules.GetSize();z2++)
				{
					i = 0;
					for (z3=0;z3<m_oaVoroMolecules.GetSize();z3++)
					{
						if (((CVoroMolecule*)m_oaVoroMolecules[z3])->m_iMolecule != z2)
							continue;
						if (va->m_pNbhTempMol[z3] != 0)
							i++;
					}
					va->m_pNbhMolecules[z2]->AddToBin_Int(i);
				}
			} else
			{
				mprintf("\nWarning: Compute_Cell failed.");
			}
		} while (vl.inc());
	}

	for (z=0;z<g_iGesAtomCount;z++)
	{
		if (!m_pAtomTouched[z])
			mprintf("\nWarning: Atom %s(%d) %s%d not touched by voronoi!",((CMolecule*)g_oaMolecules[g_waAtomMolIndex[z]])->m_sName,g_laAtomSMIndex[z]+1,((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,g_waAtomMolNumber[z]+1);
	}

	if ((co != g_iGesAtomCount) || (fabs(vg*1000.0-g_fBoxX*g_fBoxY*g_fBoxZ/1000000.0)/(g_fBoxX*g_fBoxY*g_fBoxZ/1000000.0) > 0.001))
	{
		mprintf("\nVoro++: Problem. %d/%d atoms touched. V=%f/%f A^3.",co,g_iGesAtomCount,vg*1000.0,g_fBoxX*g_fBoxY*g_fBoxZ/1000000.0);
	}

	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		vm->m_pVolume->AddToBin(vm->tempvol);
		vm->m_pSurfaceArea->AddToBin(vm->tempsurf);
		vm->m_pFaces->AddToBin_Int(vm->tempfaces);
		vm->m_pAVRatio->AddToBin(10.6347231054330961*vm->tempvol/pow(vm->tempsurf,1.5));

		for (z2=0;z2<m_oaVoroAtoms.GetSize();z2++)
		{
			i = 0;
			for (z3=0;z3<g_iGesAtomCount;z3++)
			{
				if (m_pAssignAtoms[z3] != z2)
					continue;
				if (vm->m_pNbhTempAtoms[z3] != 0)
					i++;
			}
			vm->m_pNbhAtoms[z2]->AddToBin_Int(i);
		}

		for (z2=0;z2<g_oaMolecules.GetSize();z2++)
		{
			i = 0;
			for (z3=0;z3<m_oaVoroMolecules.GetSize();z3++)
			{
				if (((CVoroMolecule*)m_oaVoroMolecules[z3])->m_iMolecule != z2)
					continue;
				if (vm->m_pNbhTempMol[z3] != 0)
					i++;
			}
			vm->m_pNbhMolecules[z2]->AddToBin_Int(i);
		}
		
		for (z2=0;z2<g_oaMolecules.GetSize();z2++)
		{
			vm->m_pNbhAreaMolecules[z2]->AddToBin(vm->m_pNbhTempArea[z2]);
		}
	}

	delete[] nbtemp;
	delete m_pContainer;
}


void CVoroWrapper::Init()
{
	int z, z2, z3, z4, i;
	CVoroAtom *va;
	CVoroMolecule *vm;
	CMolecule *m;
	CSingleMolecule *sm;

	ParseVoronoiRadii();

	m_iBlocksX = (int)(pow(g_iGesAtomCount/optimal_particles/g_fBoxX/g_fBoxY/g_fBoxZ,1/3.0)*g_fBoxX)+1;
	m_iBlocksY = (int)(pow(g_iGesAtomCount/optimal_particles/g_fBoxX/g_fBoxY/g_fBoxZ,1/3.0)*g_fBoxY)+1;
	m_iBlocksZ = (int)(pow(g_iGesAtomCount/optimal_particles/g_fBoxX/g_fBoxY/g_fBoxZ,1/3.0)*g_fBoxZ)+1;

	m_fBoxDens = g_iGesAtomCount / g_fBoxX / g_fBoxY / g_fBoxZ * 1000000.0; // Particles / Angstrom^3

	if (m_bVoroStat || m_bSurfCover || m_bWritePOV || g_bVoid || g_bTegri || g_bDomA)
	{
		try { m_pAssignAtoms = new long[g_iGesAtomCount]; } catch(...) { m_pAssignAtoms = NULL; }
		if (m_pAssignAtoms == NULL) NewException((double)g_iGesAtomCount*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { m_pAssignMolecules = new long[g_iGesAtomCount]; } catch(...) { m_pAssignMolecules = NULL; }
		if (m_pAssignMolecules == NULL) NewException((double)g_iGesAtomCount*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { m_pAssignMoleculeTypes = new long[g_iGesAtomCount]; } catch(...) { m_pAssignMoleculeTypes = NULL; }
		if (m_pAssignMoleculeTypes == NULL) NewException((double)g_iGesAtomCount*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		try { m_pAtomTouched = new bool[g_iGesAtomCount]; } catch(...) { m_pAtomTouched = NULL; }
		if (m_pAtomTouched == NULL) NewException((double)g_iGesAtomCount*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		i = 0;
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			for (z2=0;z2<m->m_waAtomCount.GetSize()-1;z2++)
			{
				for (z3=0;z3<m->m_waAtomCount[z2];z3++)
				{
					try { va = new CVoroAtom(); } catch(...) { va = NULL; }
					if (va == NULL) NewException((double)sizeof(CVoroAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					va->m_pParent = this;
					va->m_iAtomType = z2;
					va->m_iAtom = z3;
					va->m_iRealAtomType = m->m_baAtomIndex[z2];
					va->m_iMolecule = z;

					m_oaVoroAtoms.Add(va);

					for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
					{
						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
	//					mprintf("Mol %d, AT %d, A %d, Sm %d, Offset %f, i=%d.\n",z+1,z2+1,z3+1,z4+1,((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3),i);
						m_pAssignAtoms[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)] = i;
					}

					i++;
				}
			}
		}


		i = 0;
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];

				try { vm = new CVoroMolecule(); } catch(...) { vm = NULL; }
				if (vm == NULL) NewException((double)sizeof(CVoroMolecule),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				vm->m_pParent = this;
				vm->m_iMolecule = z;
				vm->m_iSingleMol = z4;
				vm->m_iCenterOffset = ((CxIntArray*)sm->m_oaAtomOffset[sm->m_oaAtomOffset.GetSize()-1])->GetAt(1); // 1 is #2, which is Center of Mass

				m_oaVoroMolecules.Add(vm);

				for (z2=0;z2<m->m_waAtomCount.GetSize()-1;z2++)
				{
					for (z3=0;z3<m->m_waAtomCount[z2];z3++)
					{
						m_pAssignMolecules[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)] = i;
						m_pAssignMoleculeTypes[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)] = z;
					}
				}

				i++;
			}
		}
	}

	g_iVoroMemory = 16;
}


void CVoroWrapper::Init2()
{
	int z;

	if (m_bVoroStat)
	{
		for (z=0;z<m_oaVoroAtoms.GetSize();z++)
		{
			((CVoroAtom*)m_oaVoroAtoms[z])->InitAnalyses();
			try { ((CVoroAtom*)m_oaVoroAtoms[z])->m_pNbhTempMol = new int[m_oaVoroMolecules.GetSize()]; } catch(...) { ((CVoroAtom*)m_oaVoroAtoms[z])->m_pNbhTempMol = NULL; }
			if (((CVoroAtom*)m_oaVoroAtoms[z])->m_pNbhTempMol == NULL) NewException((double)m_oaVoroMolecules.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}

		for (z=0;z<m_oaVoroMolecules.GetSize();z++)
		{
			((CVoroMolecule*)m_oaVoroMolecules[z])->InitAnalyses();

			try { ((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempMol = new int[m_oaVoroMolecules.GetSize()]; } catch(...) { ((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempMol = NULL; }
			if (((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempMol == NULL) NewException((double)m_oaVoroMolecules.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			try { ((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempArea = new double[m_oaVoroMolecules.GetSize()]; } catch(...) { ((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempArea = NULL; }
			if (((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempArea == NULL) NewException((double)m_oaVoroMolecules.GetSize()*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			try { ((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempAtoms = new int[g_iGesAtomCount]; } catch(...) { ((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempAtoms = NULL; }
			if (((CVoroMolecule*)m_oaVoroMolecules[z])->m_pNbhTempAtoms == NULL) NewException((double)g_iGesAtomCount*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}
	}
}


void CVoroWrapper::Parse()
{
	CTimeStep *t;
	CMolecule *m;
	CSingleMolecule *sm;
//	unsigned long col;
//	char buf[256];
	CxString buf;
	int z, z2, z3;
	CxIntArray tia, tia2;
	CxDoubleArray *fa;


	mprintf(WHITE,">>> Voronoi Analysis >>>\n\n");

	m_bVoroStat = AskYesNo("    Compute Voronoi statistics (y/n)? [yes] ",true);

	if (m_bVoroStat)
	{
		mprintf("\n    Please choose the print detail level of the Voronoi statistics:\n");
		mprintf("      0 - Only write min. / max. / average / std.dev. of basic quantities per atom / molecule type\n");
		mprintf("      1 - Write some additional info (e.g., neighborhood probabilities) per atom / molecule type\n");
		mprintf("      2 - Write the histograms of these quantities to additional files\n");
		mprintf("      3 - Also write neighborhood histograms for all pairs of atom / molecule types\n\n");
		mprintf("      Values >= 1 will create a \"voronoi\" subdirectory to host the detailed result files.\n\n");
		g_iVoroPrintLevel = AskRangeInteger("    Enter print detail level of Voronoi statistics (0-3): [0] ",0,3,0);
		mprintf("\n");
	}

	m_bSurfCover = AskYesNo("    Calculate surface coverage of one molecule (y/n)? [no] ",false);

	if (m_bSurfCover)
	{
		if (g_oaMolecules.GetSize() > 1)
		{
			m_iSurfCoverMol = AskRangeInteger_ND("    Use Voronoi cell of which molecule type (1-%d)? ",1,g_oaMolecules.GetSize(),g_oaMolecules.GetSize())-1;
		} else m_iSurfCoverMol = 0;

		m = (CMolecule*)g_oaMolecules[m_iSurfCoverMol];

		if (m->m_laSingleMolIndex.GetSize() > 1)
/*			m_iSurfCoverSM = AskInteger_ND("    Which representants of %s to use (1-%d)? ",m->m_sName,m->m_laSingleMolIndex.GetSize())-1;
				else m_iSurfCoverSM = 0;*/
		{
_surfsmagain:
			AskString("    Which representants of %s to use (allowed 1-%d, eg. 5-7,9)? [all] ",&buf,"",m->m_sName,m->m_laSingleMolIndex.GetSize());
			if (strlen(buf) == 0)
			{
				for (z=0;z<m->m_laSingleMolIndex.GetSize();z++)
					m_waSurfCoverSM.Add(z+1);
			} else ParseIntList(buf,&m_waSurfCoverSM);
			for (z=0;z<m_waSurfCoverSM.GetSize();z++)
			{
				if ((m_waSurfCoverSM[z] < 1) || (m_waSurfCoverSM[z] > m->m_laSingleMolIndex.GetSize()))
				{
					eprintf("Error: %d is not within the allowed range.\n",m_waSurfCoverSM[z]);
					goto _surfsmagain;
				}
				m_waSurfCoverSM[z]--;
			}
		} else m_waSurfCoverSM.Add(0);

		mprintf("\n    Observing voronoi cells of %d molecules.\n\n",m_waSurfCoverSM.GetSize());

		g_iFixMol = m_iSurfCoverMol;

		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			for (z2=0;z2<m->m_waAtomCount.GetSize()-1;z2++)
			{
//				at = (CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]];
				for (z3=0;z3<m->m_waAtomCount[z2];z3++)
				{
					try { fa = new CxDoubleArray("CVoroWrapper::Parse():fa"); } catch(...) { fa = NULL; }
					if (fa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					fa->SetGrow(1000);
					m_oaSurfCoverData.Add(fa);
					m_waSurfCoverMol.Add(z);
					m_waSurfCoverElem.Add(z2);
					m_waSurfCoverAtom.Add(z3);
				}
			}
		}
	}

/***************************************************************************************************************************/

	m_bWritePOV = false;


/******************************************************************************************************************************/

	m_bPOVBox = false;


/***************************************************************************************************************************/

	m_bVoroMetric = false;


	mprintf("\n");
	mprintf("    Initializing Voronoi data structures...\n");
	g_pVoroWrapper->Init();
	mprintf("\n");

	mprintf("*** Voro: Box density is %f particles / Angstrom^3.\n",m_fBoxDens);
	mprintf("*** Voro: Using %d x %d x %d blocks.\n",m_iBlocksX,m_iBlocksY,m_iBlocksZ);
	mprintf("*** Voro: Using cell memory for %d particles.\n",g_iVoroMemory);

	if (m_bVoroStat || m_bSurfCover || m_bWritePOV)
	{
		mprintf(WHITE,"*** Dumping Voronoi Infos to \"voro.txt\".\n");

		try { t = new CTimeStep(); } catch(...) { t = NULL; }
		if (t == NULL) NewException((double)sizeof(CTimeStep),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		t->CopyFrom(&g_TimeStep);

		m = (CMolecule*)g_oaMolecules[0];
		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];

		t->CenterPos(t->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(1)]);
		t->CenterPos(CxVector3(-g_fBoxX/2.0,-g_fBoxY/2.0,-g_fBoxZ/2.0));
		t->FoldAtomsPositive();
		g_pVoroWrapper->Dump("voro.txt",t);
	//	g_pVoroWrapper->WritePOV(t,"voro.pov",0,0);
		delete t;
		g_pVoroWrapper->Init2();
	}

//	AskYesNo("    This is hard-coded for Ulrikes Bmim-Br ^^ [accept] ",false);

	mprintf(WHITE,"\n<<< End of Voronoi Analysis <<<\n\n");
}


void CVoroWrapper::Dump(const char *s, CTimeStep *ts)
{
	int ijk, q, z, z2, ci, id;
	double *pp, cv;
	FILE *a;
	voronoicell_neighbor c;
	vector<int> nb;
	vector<int> fo;
	vector<double> fa;
	CVoroAtom *va, *va2;
	CVoroMolecule *vm, *vm2;

	a = OpenFileWrite(s,true);

	mfprintf(a,"*** Voro: Box density is %f particles / Angstrom^3.\n",m_fBoxDens);
	mfprintf(a,"*** Voro: Using %d x %d x %d blocks.\n",m_iBlocksX,m_iBlocksY,m_iBlocksZ);

//	d0 = g_iGesAtomCount / g_fBoxX / g_fBoxY / g_fBoxZ * 1000000.0;
//	d = 2.35 / pow(d0,1/3.0);

	mfprintf(a,"*** Voro: Creating Container...\n");
	fflush(a);

//	m_pContainer = new container_periodic_poly(1.0,0,1.0,0,0,1.0,m_iBlocksX,m_iBlocksY,m_iBlocksZ,g_iVoroMemory);

	try { m_pContainer = new container_periodic_poly(g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0,m_iBlocksX,m_iBlocksY,m_iBlocksZ,g_iVoroMemory); } catch(...) { m_pContainer = NULL; }
	if (m_pContainer == NULL) NewException((double)sizeof(container_periodic_poly),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	mfprintf(a,"*** Voro: Adding %d particles...\n",g_iGesAtomCount);
	fflush(a);

	for (z=0;z<g_iGesAtomCount;z++)
		m_pContainer->put(z,ts->m_vaCoords[z][0]/1000.0,ts->m_vaCoords[z][1]/1000.0,ts->m_vaCoords[z][2]/1000.0,g_faVoronoiRadii[z]/1000.0);

	m_fMaxVol = 0;
	m_fMaxSurf = 0;

	c_loop_all_periodic vl(*m_pContainer);

	ci = 0;
	cv = 0;
	if (vl.start()) 
	{
		do 
		{
			if (m_pContainer->compute_cell(c,vl))
			{
				ci++;
				ijk=vl.ijk;
				q=vl.q;
				pp=m_pContainer->p[ijk]+m_pContainer->ps*q;
				id = m_pContainer->id[ijk][q];

				va = (CVoroAtom*)m_oaVoroAtoms[m_pAssignAtoms[id]];
				vm = (CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[id]];

				mfprintf(a,"Cell %5d (%s[%d] %s%d); p=(%.2f | %.2f | %.2f); r=%.2f pm; %d Faces; Vol=%.4f A^3, Surface=%.4f A^2; Max.Radius=%.3f pm\n",id,((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,vm->m_iSingleMol+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType])->m_sName,va->m_iAtom+1,pp[0]*1000.0,pp[1]*1000.0,pp[2]*1000.0,pp[3]*1000.0,c.number_of_faces(),c.volume()*1000.0,c.surface_area()*100.0,sqrt(c.max_radius_squared())*1000.0);

				if (c.volume()*1000.0 > m_fMaxVol)
					m_fMaxVol = c.volume()*1000.0;

				cv += c.volume()*1000.0;

				if (c.surface_area()*100.0 > m_fMaxSurf)
					m_fMaxSurf = c.surface_area()*100.0;

				c.face_areas(fa);
				c.face_orders(fo);
				c.neighbors(nb);

				for (z2=0;z2<c.number_of_faces();z2++)
				{
					va2 = (CVoroAtom*)m_oaVoroAtoms[m_pAssignAtoms[nb[z2]]];
					vm2 = (CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[nb[z2]]];
					mfprintf(a,"  - Face %2d: %2d Vertices, Area=%7.4f A^2, Neighbor is %5d (%s[%d] %s%d)\n",z2+1,fo[z2],fa[z2]*100.0,nb[z2],((CMolecule*)g_oaMolecules[va2->m_iMolecule])->m_sName,vm2->m_iSingleMol+1,((CAtom*)g_oaAtoms[va2->m_iRealAtomType])->m_sName,va2->m_iAtom+1);
				}
			}
		} while (vl.inc());
	}
	mprintf("    %d/%d Particles touched.\n",ci,g_iGesAtomCount);
	mprintf("    Volume %f/%f Angstrom^3.\n",cv,g_fBoxX*g_fBoxY*g_fBoxZ/1000000.0);
	m_fMaxVol *= 3.0;
	m_fMaxSurf *= 3.0;
	delete m_pContainer;
	mfprintf(a,"*** The End ***\n");
	fclose(a);
}


void CVoroWrapper::Finish()
{
	int z, z2;
//	char buf[256];
	CxString buf;
	CVoroAtom *va;
	CVoroMolecule *vm, *vmbase;

	if (m_bPOVBox || m_bWritePOV)
	{
		if (m_bWritePOVMovie)
		{
			for (z=0;z<m_iPOVThreads;z++)
				fclose(m_fPOVBoxScript[z]);
		}
		mprintf(YELLOW,"\n    Note: ");
		mprintf("Please use POV-Ray 3.7 to render the pov files.\n");
	}

#ifdef TARGET_LINUX
	if (m_bWritePOV)
	{
		mprintf("\n");
		if (m_bWritePOVMovie)
		{
			if (m_iPOVThreads == 1)
			{
//				sprintf(buf,"chmod 755 POV_Single/render_single");
				buf.sprintf("chmod 755 POV_Single/render_single");
				mprintf("    Executing \"%s\" ...\n",(const char*)buf);
				system(buf);
			} else
			{
				for (z=0;z<m_iPOVThreads;z++)
				{
//					sprintf(buf,"chmod 755 POV_Single/%d/render_single",z+1);
					buf.sprintf("chmod 755 POV_Single/%d/render_single",z+1);
					mprintf("    Executing \"%s\" ...\n",(const char*)buf);
					system(buf);
				}
			}
//			sprintf(buf,"chmod 755 POV_Single/run_mencoder");
			buf.sprintf("chmod 755 POV_Single/run_mencoder");
			mprintf("    Executing \"%s\" ...\n",(const char*)buf);
			system(buf);
		} else
		{
//			sprintf(buf,"chmod 755 render_single");
			buf.sprintf("chmod 755 render_single");
			mprintf("    Executing \"%s\" ...\n",(const char*)buf);
			system(buf);
		}
	}

	if (m_bPOVBox)
	{
		mprintf("\n");
		if (m_bWritePOVMovie)
		{
			if (m_iPOVThreads == 1)
			{
//				sprintf(buf,"chmod 755 POV_Box/render_box");
				buf.sprintf("chmod 755 POV_Box/render_box");
				mprintf("    Executing \"%s\" ...\n",(const char*)buf);
				system(buf);
			} else
			{
				for (z=0;z<m_iPOVThreads;z++)
				{
//					sprintf(buf,"chmod 755 POV_Box/%d/render_box",z+1);
					buf.sprintf("chmod 755 POV_Box/%d/render_box",z+1);
					mprintf("    Executing \"%s\" ...\n",(const char*)buf);
					system(buf);
				}
			}
//			sprintf(buf,"chmod 755 POV_Box/run_mencoder");
			buf.sprintf("chmod 755 POV_Box/run_mencoder");
			mprintf("    Executing \"%s\" ...\n",(const char*)buf);
			system(buf);
		} else
		{
//			sprintf(buf,"chmod 755 render_box");
			buf.sprintf("chmod 755 render_box");
			mprintf("    Executing \"%s\" ...\n",(const char*)buf);
			system(buf);
		}
	}
#endif

	if (m_bVoroMetric)
	{
		if (m_bVoroMetricCDF)
		{
			m_pVoroMetricCDF->WriteMathematicaNb("","vorometric.nb","",false);
			m_pVoroMetricCDF2->WriteMathematicaNb("","vorometric2.nb","",false);
		}
	}

//	if (m_bWritePOV)
//		fclose(m_fPOVScript);

	if (m_bSurfCover)
	{
		mprintf(WHITE,"*** Voronoi Surface Coverage\n\n");
		mprintf("    Writing Surface data to voro_surfcover.csv\n");
		FinishSurfCover("voro_surfcover.csv");
		mprintf("\n");
	}

	if (m_bVoroStat)
	{
		mprintf(WHITE,"*** Voronoi Functions\n\n");

		if (g_iVoroPrintLevel >= 1)
		{
			mprintf("    Creating the sub-directory \"voronoi\"...\n");
			system("mkdir voronoi");
			mprintf("\n");

			mprintf("    Atom Information:\n");
			for (z=0;z<m_oaVoroAtoms.GetSize();z++)
			{
				va = (CVoroAtom*)m_oaVoroAtoms[z];
	//			sprintf(buf,"voronoi/voro_atom_%s_%s%d",((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[va->m_iRealAtomType])->m_sName,va->m_iAtom+1);
				buf.sprintf("voronoi/voro_atom_%s_%s%d",((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[va->m_iRealAtomType])->m_sName,va->m_iAtom+1);
				mprintf("      - Writing %s ...\n",(const char*)buf);
				va->Dump(buf);
			}

			mprintf("\n    Molecule Information:\n");
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				vmbase = NULL;
				for (z2=0;z2<m_oaVoroMolecules.GetSize();z2++)
				{
					vm = (CVoroMolecule*)m_oaVoroMolecules[z2];
					if (vm == NULL)
						continue;
					if (vm->m_iMolecule != z)
						continue;
					if (vmbase == NULL)
					{
						vmbase = vm;
					} else
					{
						vmbase->Add(vm);
						delete vm;
						m_oaVoroMolecules[z2] = NULL;
					}
				}
	//			sprintf(buf,"voronoi/voro_molecule_%s",((CMolecule*)g_oaMolecules[vmbase->m_iMolecule])->m_sName);
				buf.sprintf("voronoi/voro_molecule_%s",((CMolecule*)g_oaMolecules[vmbase->m_iMolecule])->m_sName);
				mprintf("      - Writing %s ...\n",(const char*)buf);
				vmbase->Dump(buf);
			}
			mprintf("\n");
		}

		mprintf("    Mean Atom Neighbor Count Matrix:\n");
		mprintf("      Writing voro_atom_nbcount.csv ...\n");
		WriteNbAtomCountMatrix("voro_atom_nbcount.csv");

		mprintf("\n    Mean Atom Neighbor Distance Matrix:\n");
		mprintf("      Writing voro_atom_nbdist.csv ...\n");
		WriteNbAtomDistMatrix("voro_atom_nbdist.csv");

		mprintf("\n    Mean Molecule Neighbor Count Matrix:\n");
		mprintf("      Writing voro_molecule_nbcount.csv ...\n");
		WriteNbMoleculeCountMatrix("voro_molecule_nbcount.csv");

		mprintf("\n    Mean Molecule Neighbor Distance Matrix:\n");
		mprintf("      Writing voro_molecule_nbdist.csv ...\n");
		WriteNbMoleculeDistMatrix("voro_molecule_nbdist.csv");

		mprintf("\n    Mean Molecule Neighbor Area Matrix:\n");
		mprintf("      Writing voro_molecule_nbarea.csv ...\n");
		WriteNbMoleculeAreaMatrix("voro_molecule_nbarea.csv");
		
		mprintf("\n    Atom information table:\n");
		mprintf("      Writing voro_atoms.csv ...\n");
		WriteAtomInfo("voro_atoms.csv");

		mprintf("\n    Molecule information table:\n");
		mprintf("      Writing voro_molecules.csv ...\n");
		WriteMoleculeInfo("voro_molecules.csv");
		mprintf("\n");
	}
}


CVoroAtom::CVoroAtom()
{
}


CVoroAtom::~CVoroAtom()
{
}


void CVoroAtom::InitAnalyses()
{
	int z;

	try { m_pSurfaceArea = new CDF(); } catch(...) { m_pSurfaceArea = NULL; }
	if (m_pSurfaceArea == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pSurfaceArea->m_fMinVal = 0;
	m_pSurfaceArea->m_fMaxVal = m_pParent->m_fMaxSurf;
	m_pSurfaceArea->m_iResolution = 10000;
	m_pSurfaceArea->SetLabelX("Surface Area [A^2]");
	m_pSurfaceArea->SetLabelY("Occurrence");
	m_pSurfaceArea->Create();

	try { m_pExposedSurface = new CDF(); } catch(...) { m_pExposedSurface = NULL; }
	if (m_pExposedSurface == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pExposedSurface->m_fMinVal = 0;
	m_pExposedSurface->m_fMaxVal = m_pParent->m_fMaxSurf;
	m_pExposedSurface->m_iResolution = 10000;
	m_pExposedSurface->SetLabelX("Exposed Surface Area [A^2]");
	m_pExposedSurface->SetLabelY("Occurrence");
	m_pExposedSurface->Create();

	try { m_pExposedSurfacePerc = new CDF(); } catch(...) { m_pExposedSurfacePerc = NULL; }
	if (m_pExposedSurfacePerc == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pExposedSurfacePerc->m_fMinVal = 0;
	m_pExposedSurfacePerc->m_fMaxVal = 100.0;
	m_pExposedSurfacePerc->m_iResolution = 10000;
	m_pExposedSurfacePerc->SetLabelX("Exposed Surface Percentage");
	m_pExposedSurfacePerc->SetLabelY("Occurrence");
	m_pExposedSurfacePerc->Create();

	try { m_pVolume = new CDF(); } catch(...) { m_pVolume = NULL; }
	if (m_pVolume == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pVolume->m_fMinVal = 0;
	m_pVolume->m_fMaxVal = m_pParent->m_fMaxVol;
	m_pVolume->m_iResolution = 10000;
	m_pVolume->SetLabelX("Volume [A^3]");
	m_pVolume->SetLabelY("Occurrence");
	m_pVolume->Create();

	try { m_pMaxRadius = new CDF(); } catch(...) { m_pMaxRadius = NULL; }
	if (m_pMaxRadius == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pMaxRadius->m_fMinVal = 0;
	m_pMaxRadius->m_fMaxVal = g_fBoxX/4.0;
	m_pMaxRadius->m_iResolution = int(g_fBoxX*0.2);
	m_pMaxRadius->SetLabelX("Max. Radius [pm]");
	m_pMaxRadius->SetLabelY("Occurrence");
	m_pMaxRadius->Create();

	try { m_pAVRatio = new CDF(); } catch(...) { m_pAVRatio = NULL; }
	if (m_pAVRatio == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pAVRatio->m_fMinVal = 0;
	m_pAVRatio->m_fMaxVal = 1.0;
	m_pAVRatio->m_iResolution = 1000;
	m_pAVRatio->SetLabelX("A/V Ratio");
	m_pAVRatio->SetLabelY("Occurrence");
	m_pAVRatio->Create();

	try { m_pFaces = new CDF(); } catch(...) { m_pFaces = NULL; }
	if (m_pFaces == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pFaces->m_fMinVal = 0;
	m_pFaces->m_fMaxVal = 100;
	m_pFaces->m_iResolution = 101;
	m_pFaces->SetLabelX("Face Count");
	m_pFaces->SetLabelY("Occurrence");
	m_pFaces->Create();

	try { m_pFaceOrders = new CDF(); } catch(...) { m_pFaceOrders = NULL; }
	if (m_pFaceOrders == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pFaceOrders->m_fMinVal = 0;
	m_pFaceOrders->m_fMaxVal = 50;
	m_pFaceOrders->m_iResolution = 51;
	m_pFaceOrders->SetLabelX("Face Order");
	m_pFaceOrders->SetLabelY("Occurrence");
	m_pFaceOrders->Create();

	try { m_pNbhAtoms = new CDF*[m_pParent->m_oaVoroAtoms.GetSize()]; } catch(...) { m_pNbhAtoms = NULL; }
	if (m_pNbhAtoms == NULL) NewException((double)m_pParent->m_oaVoroAtoms.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
		try { m_pNbhAtoms[z] = new CDF(); } catch(...) { m_pNbhAtoms[z] = NULL; }
		if (m_pNbhAtoms[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhAtoms[z]->m_fMinVal = 0;
		m_pNbhAtoms[z]->m_fMaxVal = 50;
		m_pNbhAtoms[z]->m_iResolution = 51;
		m_pNbhAtoms[z]->SetLabelX("Neighbor Count");
		m_pNbhAtoms[z]->SetLabelY("Occurrence");
		m_pNbhAtoms[z]->Create();
	}

	try { m_pNbhDistAtoms = new CDF*[m_pParent->m_oaVoroAtoms.GetSize()]; } catch(...) { m_pNbhDistAtoms = NULL; }
	if (m_pNbhDistAtoms == NULL) NewException((double)m_pParent->m_oaVoroAtoms.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
		try { m_pNbhDistAtoms[z] = new CDF(); } catch(...) { m_pNbhDistAtoms[z] = NULL; }
		if (m_pNbhDistAtoms[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhDistAtoms[z]->m_fMinVal = 0;
		m_pNbhDistAtoms[z]->m_fMaxVal = g_fBoxX*2.0;
		m_pNbhDistAtoms[z]->m_iResolution = int(g_fBoxX*0.5);
		m_pNbhDistAtoms[z]->SetLabelX("Neighbor Distance [pm]");
		m_pNbhDistAtoms[z]->SetLabelY("Occurrence");
		m_pNbhDistAtoms[z]->Create();
	}

	try { m_pNbhDistMolecules = new CDF*[g_oaMolecules.GetSize()]; } catch(...) { m_pNbhDistMolecules = NULL; }
	if (m_pNbhDistMolecules == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		try { m_pNbhDistMolecules[z] = new CDF(); } catch(...) { m_pNbhDistMolecules[z] = NULL; }
		if (m_pNbhDistMolecules[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhDistMolecules[z]->m_fMinVal = 0;
		m_pNbhDistMolecules[z]->m_fMaxVal = g_fBoxX*2.0;
		m_pNbhDistMolecules[z]->m_iResolution = int(g_fBoxX*0.5);
		m_pNbhDistMolecules[z]->SetLabelX("Neighbor Distance [pm]");
		m_pNbhDistMolecules[z]->SetLabelY("Occurrence");
		m_pNbhDistMolecules[z]->Create();
	}

	try { m_pNbhMolecules = new CDF*[g_oaMolecules.GetSize()]; } catch(...) { m_pNbhMolecules = NULL; }
	if (m_pNbhMolecules == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		try { m_pNbhMolecules[z] = new CDF(); } catch(...) { m_pNbhMolecules[z] = NULL; }
		if (m_pNbhMolecules[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhMolecules[z]->m_fMinVal = 0;
		m_pNbhMolecules[z]->m_fMaxVal = 50;
		m_pNbhMolecules[z]->m_iResolution = 51;
		m_pNbhMolecules[z]->SetLabelX("Neighbor Count");
		m_pNbhMolecules[z]->SetLabelY("Occurrence");
		m_pNbhMolecules[z]->Create();
	}
}


CVoroMolecule::CVoroMolecule()
{
}


CVoroMolecule::~CVoroMolecule()
{
}


void CVoroMolecule::InitAnalyses()
{
	int z;

	try { m_pSurfaceArea = new CDF(); } catch(...) { m_pSurfaceArea = NULL; }
	if (m_pSurfaceArea == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pSurfaceArea->m_fMinVal = 0;
	m_pSurfaceArea->m_fMaxVal = m_pParent->m_fMaxSurf*((CMolecule*)g_oaMolecules[m_iMolecule])->m_iAtomGes;
	m_pSurfaceArea->m_iResolution = 1000;
	m_pSurfaceArea->SetLabelX("Surface Area [A^2]");
	m_pSurfaceArea->SetLabelY("Occurrence");
	m_pSurfaceArea->Create();

	try { m_pVolume = new CDF(); } catch(...) { m_pVolume = NULL; }
	if (m_pVolume == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pVolume->m_fMinVal = 0;
	m_pVolume->m_fMaxVal = m_pParent->m_fMaxVol*((CMolecule*)g_oaMolecules[m_iMolecule])->m_iAtomGes;
	m_pVolume->m_iResolution = 1000;
	m_pVolume->SetLabelX("Volume [A^3]");
	m_pVolume->SetLabelY("Occurrence");
	m_pVolume->Create();

	try { m_pAVRatio = new CDF(); } catch(...) { m_pAVRatio = NULL; }
	if (m_pAVRatio == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pAVRatio->m_fMinVal = 0;
	m_pAVRatio->m_fMaxVal = 1.0;
	m_pAVRatio->m_iResolution = 1000;
	m_pAVRatio->SetLabelX("A/V Ratio");
	m_pAVRatio->SetLabelY("Occurrence");
	m_pAVRatio->Create();

	try { m_pFaces = new CDF(); } catch(...) { m_pFaces = NULL; }
	if (m_pFaces == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pFaces->m_fMinVal = 0;
	m_pFaces->m_fMaxVal = 100*((CMolecule*)g_oaMolecules[m_iMolecule])->m_iAtomGes;
	m_pFaces->m_iResolution = int(100*((CMolecule*)g_oaMolecules[m_iMolecule])->m_iAtomGes)+1;
	m_pFaces->SetLabelX("Face count");
	m_pFaces->SetLabelY("Occurrence");
	m_pFaces->Create();

	try { m_pNbhAtoms = new CDF*[m_pParent->m_oaVoroAtoms.GetSize()]; } catch(...) { m_pNbhAtoms = NULL; }
	if (m_pNbhAtoms == NULL) NewException((double)m_pParent->m_oaVoroAtoms.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
		try { m_pNbhAtoms[z] = new CDF(); } catch(...) { m_pNbhAtoms[z] = NULL; }
		if (m_pNbhAtoms[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhAtoms[z]->m_fMinVal = 0;
		m_pNbhAtoms[z]->m_fMaxVal = 50;
		m_pNbhAtoms[z]->m_iResolution = 51;
		m_pNbhAtoms[z]->SetLabelX("Neighbor count");
		m_pNbhAtoms[z]->SetLabelY("Occurrence");
		m_pNbhAtoms[z]->Create();
	}

	try { m_pNbhDistAtoms = new CDF*[m_pParent->m_oaVoroAtoms.GetSize()]; } catch(...) { m_pNbhDistAtoms = NULL; }
	if (m_pNbhDistAtoms == NULL) NewException((double)m_pParent->m_oaVoroAtoms.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
		try { m_pNbhDistAtoms[z] = new CDF(); } catch(...) { m_pNbhDistAtoms[z] = NULL; }
		if (m_pNbhDistAtoms[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhDistAtoms[z]->m_fMinVal = 0;
		m_pNbhDistAtoms[z]->m_fMaxVal = g_fBoxX*2.0;
		m_pNbhDistAtoms[z]->m_iResolution = int(g_fBoxX*0.5);
		m_pNbhDistAtoms[z]->SetLabelX("Neighbor distance [pm]");
		m_pNbhDistAtoms[z]->SetLabelY("Occurrence");
		m_pNbhDistAtoms[z]->Create();
	}

	try { m_pNbhDistMolecules = new CDF*[g_oaMolecules.GetSize()]; } catch(...) { m_pNbhDistMolecules = NULL; }
	if (m_pNbhDistMolecules == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		try { m_pNbhDistMolecules[z] = new CDF(); } catch(...) { m_pNbhDistMolecules[z] = NULL; }
		if (m_pNbhDistMolecules[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhDistMolecules[z]->m_fMinVal = 0;
		m_pNbhDistMolecules[z]->m_fMaxVal = g_fBoxX*2.0;
		m_pNbhDistMolecules[z]->m_iResolution = int(g_fBoxX*0.5);
		m_pNbhDistMolecules[z]->SetLabelX("Neighbor distance [pm]");
		m_pNbhDistMolecules[z]->SetLabelY("Occurrence");
		m_pNbhDistMolecules[z]->Create();
	}

	try { m_pNbhAreaMolecules = new CDF*[g_oaMolecules.GetSize()]; } catch(...) { m_pNbhAreaMolecules = NULL; }
	if (m_pNbhAreaMolecules == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		try { m_pNbhAreaMolecules[z] = new CDF(); } catch(...) { m_pNbhAreaMolecules[z] = NULL; }
		if (m_pNbhAreaMolecules[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhAreaMolecules[z]->m_fMinVal = 0;
		m_pNbhAreaMolecules[z]->m_fMaxVal = g_fBoxX*g_fBoxX / 10000.0f;
		m_pNbhAreaMolecules[z]->m_iResolution = int(g_fBoxX*g_fBoxX*0.25 / 10000.0f);
		m_pNbhAreaMolecules[z]->SetLabelX("Neighbor area [A^2]");
		m_pNbhAreaMolecules[z]->SetLabelY("Occurrence");
		m_pNbhAreaMolecules[z]->Create();
	}
	
	try { m_pNbhMolecules = new CDF*[g_oaMolecules.GetSize()]; } catch(...) { m_pNbhMolecules = NULL; }
	if (m_pNbhMolecules == NULL) NewException((double)g_oaMolecules.GetSize()*sizeof(CDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		try { m_pNbhMolecules[z] = new CDF(); } catch(...) { m_pNbhMolecules[z] = NULL; }
		if (m_pNbhMolecules[z] == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pNbhMolecules[z]->m_fMinVal = 0;
		m_pNbhMolecules[z]->m_fMaxVal = 50;
		m_pNbhMolecules[z]->m_iResolution = 51;
		m_pNbhMolecules[z]->SetLabelX("Neighbor count");
		m_pNbhMolecules[z]->SetLabelY("Occurrence");
		m_pNbhMolecules[z]->Create();
	}
}


void CVoroAtom::Dump(const char *s)
{
	FILE *a;
//	char buf[256];
	CxString buf;
	int z, z2;

//	sprintf(buf,"%s.txt",s);
	buf.sprintf("%s.txt",s);
	a = OpenFileWrite(buf,true);

	mfprintf(a,"*** Voronoi Analysis for %s %s%d ***\n\n",((CMolecule*)g_oaMolecules[m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[m_iRealAtomType])->m_sName,m_iAtom+1);

	m_pVolume->CalcMeanSD();
	mfprintf(a,"  * Volume:\n");
	mfprintf(a,"      Avg. %10G Angstrom^3\n",m_pVolume->m_fMean);
	mfprintf(a,"      Min. %10G Angstrom^3\n",m_pVolume->m_fMinInput);
	mfprintf(a,"      Max. %10G Angstrom^3\n",m_pVolume->m_fMaxInput);
	mfprintf(a,"      SD   %10G Angstrom^3\n",m_pVolume->m_fSD);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_volume.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pVolume->m_fBinEntries,m_pVolume->m_fSkipEntries);
		m_pVolume->NormBinSum(100.0);
		m_pVolume->Write("",s,"_volume.csv",true);
	}

	m_pSurfaceArea->CalcMeanSD();
	mfprintf(a,"  * Surface Area:\n");
	mfprintf(a,"      Avg. %10G Angstrom^2\n",m_pSurfaceArea->m_fMean);
	mfprintf(a,"      Min. %10G Angstrom^2\n",m_pSurfaceArea->m_fMinInput);
	mfprintf(a,"      Max. %10G Angstrom^2\n",m_pSurfaceArea->m_fMaxInput);
	mfprintf(a,"      SD   %10G Angstrom^2\n",m_pSurfaceArea->m_fSD);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_surface.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pSurfaceArea->m_fBinEntries,m_pSurfaceArea->m_fSkipEntries);
		m_pSurfaceArea->NormBinSum(100.0);
		m_pSurfaceArea->Write("",s,"_surface.csv",true);
	}

	m_pExposedSurface->CalcMeanSD();
	mfprintf(a,"  * Exposed Surface Area:\n");
	mfprintf(a,"      Avg. %10G Angstrom^2\n",m_pExposedSurface->m_fMean);
	mfprintf(a,"      Min. %10G Angstrom^2\n",m_pExposedSurface->m_fMinInput);
	mfprintf(a,"      Max. %10G Angstrom^2\n",m_pExposedSurface->m_fMaxInput);
	mfprintf(a,"      SD   %10G Angstrom^2\n",m_pExposedSurface->m_fSD);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_exposed_surface.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pExposedSurface->m_fBinEntries,m_pExposedSurface->m_fSkipEntries);
		m_pExposedSurface->NormBinSum(100.0);
		m_pExposedSurface->Write("",s,"_exposed_surface.csv",true);
	}

	m_pExposedSurfacePerc->CalcMeanSD();
	mfprintf(a,"  * Exposed Surface Area Percentage:\n");
	mfprintf(a,"      Avg. %10G percent\n",m_pExposedSurfacePerc->m_fMean);
	mfprintf(a,"      Min. %10G percent\n",m_pExposedSurfacePerc->m_fMinInput);
	mfprintf(a,"      Max. %10G percent\n",m_pExposedSurfacePerc->m_fMaxInput);
	mfprintf(a,"      SD   %10G percent\n",m_pExposedSurfacePerc->m_fSD);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_exposed_surface_perc.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pExposedSurfacePerc->m_fBinEntries,m_pExposedSurfacePerc->m_fSkipEntries);
		m_pExposedSurfacePerc->NormBinSum(100.0);
		m_pExposedSurfacePerc->Write("",s,"_exposed_surface_perc.csv",true);
	}

	m_pMaxRadius->CalcMeanSD();
	mfprintf(a,"  * Maximum Radius:\n");
	mfprintf(a,"      Avg. %.4f pm\n",m_pMaxRadius->m_fMean);
	mfprintf(a,"      Min. %.4f pm\n",m_pMaxRadius->m_fMinInput);
	mfprintf(a,"      Max. %.4f pm\n",m_pMaxRadius->m_fMaxInput);
	mfprintf(a,"      SD   %.4f pm\n\n",m_pMaxRadius->m_fSD);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_maxradius.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pMaxRadius->m_fBinEntries,m_pMaxRadius->m_fSkipEntries);
		m_pMaxRadius->NormBinSum(100.0);
		m_pMaxRadius->Write("",s,"_maxradius.csv",true);
	}

	m_pFaceOrders->CalcMeanSD();
	mfprintf(a,"  * Face Orders:\n");
	mfprintf(a,"      Avg. %.6f Edges\n",m_pFaceOrders->m_fMean);
	mfprintf(a,"      Min. %.0f Edges\n",m_pFaceOrders->m_fMinInput);
	mfprintf(a,"      Max. %.0f Edges\n",m_pFaceOrders->m_fMaxInput);
	mfprintf(a,"      SD   %.6f Edges\n\n",m_pFaceOrders->m_fSD);
	m_pFaceOrders->NormBinSum(100.0);
	for (z=0;z<m_pFaceOrders->m_iResolution;z++)
		if (m_pFaceOrders->m_pBin[z] != 0)
			mfprintf(a,"      - %2d Edges: %10.6f percent\n",z,m_pFaceOrders->m_pBin[z]);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_faceorders.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pFaceOrders->m_fBinEntries,m_pFaceOrders->m_fSkipEntries);
		m_pFaceOrders->Write_Int("",s,"_faceorders.csv");
	}

	m_pFaces->CalcMeanSD();
	mfprintf(a,"  * Faces:\n");
	mfprintf(a,"      Avg. %.6f\n",m_pFaces->m_fMean);
	mfprintf(a,"      Min. %.0f\n",m_pFaces->m_fMinInput);
	mfprintf(a,"      Max. %.0f\n",m_pFaces->m_fMaxInput);
	mfprintf(a,"      SD   %.6f\n\n",m_pFaces->m_fSD);
	m_pFaces->NormBinSum(100.0);
	for (z=0;z<m_pFaces->m_iResolution;z++)
		if (m_pFaces->m_pBin[z] != 0)
			mfprintf(a,"      - %2d Faces: %10.6f percent\n",z,m_pFaces->m_pBin[z]);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_faces.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pFaces->m_fBinEntries,m_pFaces->m_fSkipEntries);
		m_pFaces->Write_Int("",s,"_faces.csv");
	}

	m_pAVRatio->CalcMeanSD();
	mfprintf(a,"  * A/V Ratio:\n");
	mfprintf(a,"      Avg. %.6f\n",m_pAVRatio->m_fMean);
	mfprintf(a,"      Min. %.6f\n",m_pAVRatio->m_fMinInput);
	mfprintf(a,"      Max. %.6f\n",m_pAVRatio->m_fMaxInput);
	mfprintf(a,"      SD   %.6f\n\n",m_pAVRatio->m_fSD);
	if (g_iVoroPrintLevel >= 2)
	{
		m_pAVRatio->NormBinSum(100.0);
		mfprintf(a,"      Saving Histogram as %s_avratio.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pAVRatio->m_fBinEntries,m_pAVRatio->m_fSkipEntries);
		m_pAVRatio->Write("",s,"_avratio.csv",false);
	}

	mfprintf(a,"\n\n###### Neighboring Molecules ######\n\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m_pNbhMolecules[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s molecules:\n",((CMolecule*)g_oaMolecules[z])->m_sName);
		mfprintf(a,"      Avg. %.6f\n",m_pNbhMolecules[z]->m_fMean);
		mfprintf(a,"      Min. %.0f\n",m_pNbhMolecules[z]->m_fMinInput);
		mfprintf(a,"      Max. %.0f\n",m_pNbhMolecules[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.6f\n",m_pNbhMolecules[z]->m_fSD);
		m_pNbhMolecules[z]->NormBinSum(100.0);
		for (z2=0;z2<m_pNbhMolecules[z]->m_iResolution;z2++)
			if (m_pNbhMolecules[z]->m_pBin[z2] != 0)
				mfprintf(a,"      - %2d Molecules: %10.6f percent\n",z2,m_pNbhMolecules[z]->m_pBin[z2]);
		if (g_iVoroPrintLevel >= 3)
		{
	//		sprintf(buf,"_nbcount_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
			buf.sprintf("_nbcount_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhMolecules[z]->m_fBinEntries,m_pNbhMolecules[z]->m_fSkipEntries);
			m_pNbhMolecules[z]->Write_Int("",s,buf);
		}
	}

	mfprintf(a,"\n\n###### Neighboring Molecule Distances ######\n\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m_pNbhDistMolecules[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s molecule distances:\n",((CMolecule*)g_oaMolecules[z])->m_sName);
		mfprintf(a,"      Avg. %.4f pm\n",m_pNbhDistMolecules[z]->m_fMean);
		mfprintf(a,"      Min. %.4f pm\n",m_pNbhDistMolecules[z]->m_fMinInput);
		mfprintf(a,"      Max. %.4f pm\n",m_pNbhDistMolecules[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.4f pm\n",m_pNbhDistMolecules[z]->m_fSD);
		if (g_iVoroPrintLevel >= 3)
		{
			m_pNbhDistMolecules[z]->NormBinSum(100.0);
	//		sprintf(buf,"_nbdist_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
			buf.sprintf("_nbdist_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhDistMolecules[z]->m_fBinEntries,m_pNbhDistMolecules[z]->m_fSkipEntries);
			m_pNbhDistMolecules[z]->Write("",s,buf,true);
		}
	}

	mfprintf(a,"\n\n###### Neighboring Atoms ######\n\n");

	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
//		sprintf(buf,"%s %s%d",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		buf.sprintf("%s %s%d",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		m_pNbhAtoms[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s atoms:\n",(const char*)buf);
		mfprintf(a,"      Avg. %.6f\n",m_pNbhAtoms[z]->m_fMean);
		mfprintf(a,"      Min. %.0f\n",m_pNbhAtoms[z]->m_fMinInput);
		mfprintf(a,"      Max. %.0f\n",m_pNbhAtoms[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.6f\n",m_pNbhAtoms[z]->m_fSD);
		m_pNbhAtoms[z]->NormBinSum(100.0);
		for (z2=0;z2<m_pNbhAtoms[z]->m_iResolution;z2++)
			if (m_pNbhAtoms[z]->m_pBin[z2] != 0)
				mfprintf(a,"      - %2d Atoms: %10.6f percent\n",z2,m_pNbhAtoms[z]->m_pBin[z2]);
		if (g_iVoroPrintLevel >= 3)
		{
	//		sprintf(buf,"_nbcount_atom_%s_%s%d.csv",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
			buf.sprintf("_nbcount_atom_%s_%s%d.csv",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhAtoms[z]->m_fBinEntries,m_pNbhAtoms[z]->m_fSkipEntries);
			m_pNbhAtoms[z]->Write_Int("",s,buf);
		}
	}

	mfprintf(a,"\n\n###### Neighboring Atom Distances ######\n\n");

	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
//		sprintf(buf,"%s %s%d",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		buf.sprintf("%s %s%d",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		m_pNbhDistAtoms[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s atom distances:\n",(const char*)buf);
		mfprintf(a,"      Avg. %.4f pm\n",m_pNbhDistAtoms[z]->m_fMean);
		mfprintf(a,"      Min. %.4f pm\n",m_pNbhDistAtoms[z]->m_fMinInput);
		mfprintf(a,"      Max. %.4f pm\n",m_pNbhDistAtoms[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.4f pm\n",m_pNbhDistAtoms[z]->m_fSD);
		if (g_iVoroPrintLevel >= 3)
		{
			m_pNbhDistAtoms[z]->NormBinSum(100.0);
	//		sprintf(buf,"_nbdist_atom_%s_%s%d.csv",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
			buf.sprintf("_nbdist_atom_%s_%s%d.csv",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhDistAtoms[z]->m_fBinEntries,m_pNbhDistAtoms[z]->m_fSkipEntries);
			m_pNbhDistAtoms[z]->Write("",s,buf,true);
		}
	}

	fclose(a);
}


void CVoroMolecule::Dump(const char *s)
{
	FILE *a;
//	char buf[256];
	CxString buf;
	int z, z2;

//	sprintf(buf,"%s.txt",s);
	buf.sprintf("%s.txt",s);
	a = OpenFileWrite(buf,true);

	mfprintf(a,"*** Voronoi Analysis for Molecule %s ***\n\n",((CMolecule*)g_oaMolecules[m_iMolecule])->m_sName);

	m_pVolume->CalcMeanSD();
	mfprintf(a,"  * Volume:\n");
	mfprintf(a,"      Avg. %10G Angstrom^3\n",m_pVolume->m_fMean);
	mfprintf(a,"      Min. %10G Angstrom^3\n",m_pVolume->m_fMinInput);
	mfprintf(a,"      Max. %10G Angstrom^3\n",m_pVolume->m_fMaxInput);
	mfprintf(a,"      SD   %10G Angstrom^3\n",m_pVolume->m_fSD);

	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_volume.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pVolume->m_fBinEntries,m_pVolume->m_fSkipEntries);
		m_pVolume->NormBinSum(100.0);
		m_pVolume->Write("",s,"_volume.csv",true);
	}

	m_pSurfaceArea->CalcMeanSD();
	mfprintf(a,"  * Surface Area:\n");
	mfprintf(a,"      Avg. %10G Angstrom^2\n",m_pSurfaceArea->m_fMean);
	mfprintf(a,"      Min. %10G Angstrom^2\n",m_pSurfaceArea->m_fMinInput);
	mfprintf(a,"      Max. %10G Angstrom^2\n",m_pSurfaceArea->m_fMaxInput);
	mfprintf(a,"      SD   %10G Angstrom^2\n",m_pSurfaceArea->m_fSD);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_surface.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pSurfaceArea->m_fBinEntries,m_pSurfaceArea->m_fSkipEntries);
		m_pSurfaceArea->NormBinSum(100.0);
		m_pSurfaceArea->Write("",s,"_surface.csv",true);
	}

	m_pAVRatio->CalcMeanSD();
	mfprintf(a,"  * A/V Ratio:\n");
	mfprintf(a,"      Avg. %10G Angstrom^2\n",m_pAVRatio->m_fMean);
	mfprintf(a,"      Min. %10G Angstrom^2\n",m_pAVRatio->m_fMinInput);
	mfprintf(a,"      Max. %10G Angstrom^2\n",m_pAVRatio->m_fMaxInput);
	mfprintf(a,"      SD   %10G Angstrom^2\n",m_pAVRatio->m_fSD);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_avratio.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pAVRatio->m_fBinEntries,m_pAVRatio->m_fSkipEntries);
		m_pAVRatio->NormBinSum(100.0);
		m_pAVRatio->Write("",s,"_avratio.csv",true);
	}

	m_pFaces->CalcMeanSD();
	mfprintf(a,"  * Faces:\n");
	mfprintf(a,"      Avg. %.6f\n",m_pFaces->m_fMean);
	mfprintf(a,"      Min. %.0f\n",m_pFaces->m_fMinInput);
	mfprintf(a,"      Max. %.0f\n",m_pFaces->m_fMaxInput);
	mfprintf(a,"      SD   %.6f\n\n",m_pFaces->m_fSD);
	m_pFaces->NormBinSum(100.0);
	for (z=0;z<m_pFaces->m_iResolution;z++)
		if (m_pFaces->m_pBin[z] != 0)
			mfprintf(a,"      - %2d Faces: %10.6f percent\n",z,m_pFaces->m_pBin[z]);
	if (g_iVoroPrintLevel >= 2)
	{
		mfprintf(a,"      Saving Histogram as %s_faces.csv.\n",s);
		mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pFaces->m_fBinEntries,m_pFaces->m_fSkipEntries);
		m_pFaces->Write_Int("",s,"_faces.csv");
	}

	mfprintf(a,"\n\n###### Neighboring Molecules ######\n\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m_pNbhMolecules[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s molecules:\n",((CMolecule*)g_oaMolecules[z])->m_sName);
		mfprintf(a,"      Avg. %.6f\n",m_pNbhMolecules[z]->m_fMean);
		mfprintf(a,"      Min. %.0f\n",m_pNbhMolecules[z]->m_fMinInput);
		mfprintf(a,"      Max. %.0f\n",m_pNbhMolecules[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.6f\n",m_pNbhMolecules[z]->m_fSD);
		m_pNbhMolecules[z]->NormBinSum(100.0);
		for (z2=0;z2<m_pNbhMolecules[z]->m_iResolution;z2++)
			if (m_pNbhMolecules[z]->m_pBin[z2] != 0)
				mfprintf(a,"      - %2d Molecules: %10.6f percent\n",z2,m_pNbhMolecules[z]->m_pBin[z2]);
//		sprintf(buf,"_nbcount_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
		if (g_iVoroPrintLevel >= 3)
		{
			buf.sprintf("_nbcount_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhMolecules[z]->m_fBinEntries,m_pNbhMolecules[z]->m_fSkipEntries);
			m_pNbhMolecules[z]->Write_Int("",s,buf);
		}
	}

	mfprintf(a,"\n\n###### Neighboring Molecule Areas ######\n\n");
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m_pNbhAreaMolecules[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s molecule areas:\n",((CMolecule*)g_oaMolecules[z])->m_sName);
		mfprintf(a,"      Avg. %.4f\n",m_pNbhAreaMolecules[z]->m_fMean);
		mfprintf(a,"      Min. %.4f\n",m_pNbhAreaMolecules[z]->m_fMinInput);
		mfprintf(a,"      Max. %.4f\n",m_pNbhAreaMolecules[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.4f\n",m_pNbhAreaMolecules[z]->m_fSD);
		if (g_iVoroPrintLevel >= 3)
		{
			m_pNbhAreaMolecules[z]->NormBinSum(100.0);
			buf.sprintf("_nbarea_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhAreaMolecules[z]->m_fBinEntries,m_pNbhAreaMolecules[z]->m_fSkipEntries);
			m_pNbhAreaMolecules[z]->Write("",s,buf,true);
		}
	}
	
	mfprintf(a,"\n\n###### Neighboring Molecule Distances ######\n\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m_pNbhDistMolecules[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s molecule distances:\n",((CMolecule*)g_oaMolecules[z])->m_sName);
		mfprintf(a,"      Avg. %.4f pm\n",m_pNbhDistMolecules[z]->m_fMean);
		mfprintf(a,"      Min. %.4f pm\n",m_pNbhDistMolecules[z]->m_fMinInput);
		mfprintf(a,"      Max. %.4f pm\n",m_pNbhDistMolecules[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.4f pm\n",m_pNbhDistMolecules[z]->m_fSD);
		if (g_iVoroPrintLevel >= 3)
		{
			m_pNbhDistMolecules[z]->NormBinSum(100.0);
	//		sprintf(buf,"_nbdist_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
			buf.sprintf("_nbdist_mol_%s.csv",((CMolecule*)g_oaMolecules[z])->m_sName);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhDistMolecules[z]->m_fBinEntries,m_pNbhDistMolecules[z]->m_fSkipEntries);
			m_pNbhDistMolecules[z]->Write("",s,buf,true);
		}
	}

	mfprintf(a,"\n\n###### Neighboring Atoms ######\n\n");

	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
//		sprintf(buf,"%s %s%d",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		buf.sprintf("%s %s%d",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		m_pNbhAtoms[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s atoms:\n",(const char *)buf);
		mfprintf(a,"      Avg. %.6f\n",m_pNbhAtoms[z]->m_fMean);
		mfprintf(a,"      Min. %.0f\n",m_pNbhAtoms[z]->m_fMinInput);
		mfprintf(a,"      Max. %.0f\n",m_pNbhAtoms[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.6f\n",m_pNbhAtoms[z]->m_fSD);
		m_pNbhAtoms[z]->NormBinSum(100.0);
		for (z2=0;z2<m_pNbhAtoms[z]->m_iResolution;z2++)
			if (m_pNbhAtoms[z]->m_pBin[z2] != 0)
				mfprintf(a,"      - %2d Atoms: %10.6f percent\n",z2,m_pNbhAtoms[z]->m_pBin[z2]);
//		sprintf(buf,"_nbcount_atom_%s_%s%d.csv",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		if (g_iVoroPrintLevel >= 3)
		{
			buf.sprintf("_nbcount_atom_%s_%s%d.csv",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhAtoms[z]->m_fBinEntries,m_pNbhAtoms[z]->m_fSkipEntries);
			m_pNbhAtoms[z]->Write_Int("",s,buf);
		}
	}

	mfprintf(a,"\n\n###### Neighboring Atom Distances ######\n\n");

	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
//		sprintf(buf,"%s %s%d",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		buf.sprintf("%s %s%d",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
		m_pNbhDistAtoms[z]->CalcMeanSD();
		mfprintf(a,"  * Neighboring %s atom distances:\n",(const char*)buf);
		mfprintf(a,"      Avg. %.4f pm\n",m_pNbhDistAtoms[z]->m_fMean);
		mfprintf(a,"      Min. %.4f pm\n",m_pNbhDistAtoms[z]->m_fMinInput);
		mfprintf(a,"      Max. %.4f pm\n",m_pNbhDistAtoms[z]->m_fMaxInput);
		mfprintf(a,"      SD   %.4f pm\n",m_pNbhDistAtoms[z]->m_fSD);
		if (g_iVoroPrintLevel >= 3)
		{
			m_pNbhDistAtoms[z]->NormBinSum(100.0);
	//		sprintf(buf,"_nbdist_atom_%s_%s%d.csv",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
			buf.sprintf("_nbdist_atom_%s_%s%d.csv",((CMolecule*)g_oaMolecules[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iRealAtomType])->m_sName,((CVoroAtom*)m_pParent->m_oaVoroAtoms[z])->m_iAtom+1);
			mfprintf(a,"      Saving Histogram as %s%s.\n",s,(const char*)buf);
			mfprintf(a,"        (%.0f bin entries, %.0f skipped)\n\n",m_pNbhDistAtoms[z]->m_fBinEntries,m_pNbhDistAtoms[z]->m_fSkipEntries);
			m_pNbhDistAtoms[z]->Write("",s,buf,true);
		}
	}

	fclose(a);
}

void CVoroMolecule::Add(CVoroMolecule *m)
{
	int z;

	for (z=0;z<m_pParent->m_oaVoroAtoms.GetSize();z++)
	{
		m_pNbhAtoms[z]->AddFrom(m->m_pNbhAtoms[z]);
		m_pNbhDistAtoms[z]->AddFrom(m->m_pNbhDistAtoms[z]);
	}
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m_pNbhMolecules[z]->AddFrom(m->m_pNbhMolecules[z]);
		m_pNbhDistMolecules[z]->AddFrom(m->m_pNbhDistMolecules[z]);
		m_pNbhAreaMolecules[z]->AddFrom(m->m_pNbhAreaMolecules[z]);
	}
	m_pVolume->AddFrom(m->m_pVolume);
	m_pSurfaceArea->AddFrom(m->m_pSurfaceArea);
	m_pAVRatio->AddFrom(m->m_pAVRatio);
}

void CVoroWrapper::WriteNbAtomCountMatrix(const char *s)
{
	FILE *a;
	int z, z2;
	CVoroAtom *va;

	a = OpenFileWrite(s,true);

	mfprintf(a,"Row=From, Column=To\n");
	mfprintf(a,"*");
	for (z=0;z<m_oaVoroAtoms.GetSize();z++)
	{
		va = (CVoroAtom*)m_oaVoroAtoms[z];
		mfprintf(a,";  %s %s%d",((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[va->m_iRealAtomType])->m_sName,va->m_iAtom+1);
	}
	mfprintf(a,"\n");
	
	for (z=0;z<m_oaVoroAtoms.GetSize();z++)
	{
		va = (CVoroAtom*)m_oaVoroAtoms[z];
		mfprintf(a,"%s %s%d",((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[va->m_iRealAtomType])->m_sName,va->m_iAtom+1);
		for (z2=0;z2<m_oaVoroAtoms.GetSize();z2++)
			mfprintf(a,";  %f",va->m_pNbhAtoms[z2]->m_fMean);
		mfprintf(a,"\n");
	}
	
	fclose(a);
}

void CVoroWrapper::WriteNbAtomDistMatrix(const char *s)
{
	FILE *a;
	int z, z2;
	CVoroAtom *va;

	a = OpenFileWrite(s,true);

	mfprintf(a,"Row=From, Column=To, all numbers in pm\n");
	mfprintf(a,"*");
	for (z=0;z<m_oaVoroAtoms.GetSize();z++)
	{
		va = (CVoroAtom*)m_oaVoroAtoms[z];
		mfprintf(a,";  %s %s%d",((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[va->m_iRealAtomType])->m_sName,va->m_iAtom+1);
	}
	mfprintf(a,"\n");
	
	for (z=0;z<m_oaVoroAtoms.GetSize();z++)
	{
		va = (CVoroAtom*)m_oaVoroAtoms[z];
		mfprintf(a,"%s %s%d",((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[va->m_iRealAtomType])->m_sName,va->m_iAtom+1);
		for (z2=0;z2<m_oaVoroAtoms.GetSize();z2++)
			mfprintf(a,";  %f",va->m_pNbhDistAtoms[z2]->m_fMean);
		mfprintf(a,"\n");
	}
	
	fclose(a);
}

void CVoroWrapper::WriteNbMoleculeCountMatrix(const char *s)
{
	FILE *a;
	int z, z2;
	CVoroMolecule *vm;

	a = OpenFileWrite(s,true);

	mfprintf(a,"Row=From, Column=To\n");
	mfprintf(a,"*");
	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		if (vm == NULL)
			continue;
		mfprintf(a,";  %s",((CMolecule*)g_oaMolecules[vm->m_iMolecule])->m_sName);
	}
	mfprintf(a,"\n");
	
	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		if (vm == NULL)
			continue;
		mfprintf(a,"%s",((CMolecule*)g_oaMolecules[vm->m_iMolecule])->m_sName);
		for (z2=0;z2<m_oaVoroMolecules.GetSize();z2++)
		{
			if (m_oaVoroMolecules[z2] == NULL)
				continue;
			mfprintf(a,";  %f",vm->m_pNbhMolecules[((CVoroMolecule*)m_oaVoroMolecules[z2])->m_iMolecule]->m_fMean);
		}
		mfprintf(a,"\n");
	}
	
	fclose(a);
}

void CVoroWrapper::WriteNbMoleculeDistMatrix(const char *s)
{
	FILE *a;
	int z, z2;
	CVoroMolecule *vm;

	a = OpenFileWrite(s,true);

	mfprintf(a,"Row=From, Column=To, all numbers in pm\n");
	mfprintf(a,"*");
	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		if (vm == NULL)
			continue;
		mfprintf(a,";  %s",((CMolecule*)g_oaMolecules[vm->m_iMolecule])->m_sName);
	}
	mfprintf(a,"\n");
	
	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		if (vm == NULL)
			continue;
		mfprintf(a,"%s",((CMolecule*)g_oaMolecules[vm->m_iMolecule])->m_sName);
		for (z2=0;z2<m_oaVoroMolecules.GetSize();z2++)
		{
			if (m_oaVoroMolecules[z2] == NULL)
				continue;
			mfprintf(a,";  %f",vm->m_pNbhDistMolecules[((CVoroMolecule*)m_oaVoroMolecules[z2])->m_iMolecule]->m_fMean);
		}
		mfprintf(a,"\n");
	}
	
	fclose(a);
}

void CVoroWrapper::WriteNbMoleculeAreaMatrix(const char *s)
{
	FILE *a;
	int z, z2;
	CVoroMolecule *vm;
	
	a = OpenFileWrite(s,true);
	
	mfprintf(a,"Row=From, Column=To, all numbers in A^2\n");
	mfprintf(a,"*");
	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		if (vm == NULL)
			continue;
		mfprintf(a,";  %s",((CMolecule*)g_oaMolecules[vm->m_iMolecule])->m_sName);
	}
	mfprintf(a,"\n");
	
	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		if (vm == NULL)
			continue;
		mfprintf(a,"%s",((CMolecule*)g_oaMolecules[vm->m_iMolecule])->m_sName);
		for (z2=0;z2<m_oaVoroMolecules.GetSize();z2++)
		{
			if (m_oaVoroMolecules[z2] == NULL)
				continue;
			mfprintf(a,";  %f",vm->m_pNbhAreaMolecules[((CVoroMolecule*)m_oaVoroMolecules[z2])->m_iMolecule]->m_fMean);
		}
		mfprintf(a,"\n");
	}
	
	fclose(a);
}

void CVoroWrapper::WriteAtomInfo(const char *s)
{
	FILE *a;
	int z;
	CVoroAtom *va;

	a = OpenFileWrite(s,true);

	mfprintf(a,"Molecule;  Atom;  Mean Volume [A^3];  Mean Surface [A^2];  Mean Exposed Surface [A^2];  Mean Exposed Surface Percentage;  Mean Faces;  Mean Face Orders;  Mean Max Radius [pm];  Mean A/V Ratio\n");

	for (z=0;z<m_oaVoroAtoms.GetSize();z++)
	{
		va = (CVoroAtom*)m_oaVoroAtoms[z];
		mfprintf(a,"%s;  %s%d;  ",((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,((CAtom*)g_oaAtoms[va->m_iRealAtomType])->m_sName,va->m_iAtom+1);
		mfprintf(a,"%f;  ",va->m_pVolume->m_fMean);
		mfprintf(a,"%f;  ",va->m_pSurfaceArea->m_fMean);
		mfprintf(a,"%f;  ",va->m_pExposedSurface->m_fMean);
		mfprintf(a,"%f;  ",va->m_pExposedSurfacePerc->m_fMean);
		mfprintf(a,"%f;  ",va->m_pFaces->m_fMean);
		mfprintf(a,"%f;  ",va->m_pFaceOrders->m_fMean);
		mfprintf(a,"%f;  ",va->m_pMaxRadius->m_fMean);
		mfprintf(a,"%f\n",va->m_pAVRatio->m_fMean);
	}
	
	fclose(a);
}

void CVoroWrapper::WriteMoleculeInfo(const char *s)
{
	FILE *a;
	int z;
	CVoroMolecule *vm;

	a = OpenFileWrite(s,true);

	mfprintf(a,"Molecule;  Mean Volume [A^3];  Mean Surface [A^2];  Mean Faces;  Mean A/V Ratio\n");

	for (z=0;z<m_oaVoroMolecules.GetSize();z++)
	{
		vm = (CVoroMolecule*)m_oaVoroMolecules[z];
		if (vm == NULL)
			continue;
		mfprintf(a,"%s;  ",((CMolecule*)g_oaMolecules[vm->m_iMolecule])->m_sName);
		mfprintf(a,"%f;  ",vm->m_pVolume->m_fMean);
		mfprintf(a,"%f;  ",vm->m_pSurfaceArea->m_fMean);
		mfprintf(a,"%f;  ",vm->m_pFaces->m_fMean);
		mfprintf(a,"%f\n",vm->m_pAVRatio->m_fMean);
	}
	
	fclose(a);
}



void CVoroWrapper::ProcessSurfCover(CTimeStep *ts)
{
	int ijk, q, z, z3, id, faces;
//	CVoroAtom *va;
	CVoroMolecule *vm;
//	CMolecule *m;
//	CSingleMolecule *sm;
	voronoicell_neighbor c;
	vector<int> nb;
	vector<int> fo;
	vector<double> fa;
	CxVector3 vec;
	double atot;
	bool *tb;
	CxDoubleArray *tfa;

//	m = (CMolecule*)g_oaMolecules[m_iSurfCoverMol];
/*	sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[m_iSurfCoverSM]];

	vec = ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[g_iFixAtomType[0]])->GetAt(g_iFixAtom[0])] - CxVector3(g_fBoxX/2.0,g_fBoxY/2.0,g_fBoxZ/2.0);

	ts->CenterPos(vec);*/

	ts->FoldAtomsPositive();

	try { m_pContainer = new container_periodic_poly(g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0,m_iBlocksX,m_iBlocksY,m_iBlocksZ,g_iVoroMemory); } catch(...) { m_pContainer = NULL; }
	if (m_pContainer == NULL) NewException((double)sizeof(container_periodic_poly),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<g_iGesAtomCount;z++)
		m_pContainer->put(z,ts->m_vaCoords[z][0]/1000.0,ts->m_vaCoords[z][1]/1000.0,ts->m_vaCoords[z][2]/1000.0,g_faVoronoiRadii[z]/1000.0);

	c_loop_all_periodic vl(*m_pContainer);

	atot = 0;

	try { tb = new bool[m_waSurfCoverAtom.GetSize()]; } catch(...) { tb = NULL; }
	if (tb == NULL) NewException((double)m_waSurfCoverAtom.GetSize()*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_waSurfCoverAtom.GetSize();z++)
		tb[z] = false;

	if (vl.start()) 
	{
		do 
		{
			if (m_pContainer->compute_cell(c,vl))
			{
				ijk=vl.ijk;
				q=vl.q;

				id = m_pContainer->id[ijk][q];

//				va = (CVoroAtom*)m_oaVoroAtoms[m_pAssignAtoms[id]];
				vm = (CVoroMolecule*)m_oaVoroMolecules[m_pAssignMolecules[id]];

				if ((vm->m_iMolecule == m_iSurfCoverMol) && (/*vm->m_iSingleMol == m_iSurfCoverSM*/m_waSurfCoverSM.Contains(vm->m_iSingleMol)))
				{
					c.neighbors(nb);
					faces = c.number_of_faces();
					c.face_areas(fa);

					for (z=0;z<faces;z++)
					{
						if (m_pAssignMolecules[id] != m_pAssignMolecules[nb[z]])
						{
							for (z3=0;z3<m_waSurfCoverMol.GetSize();z3++)
							{
								if ((g_waAtomMolIndex[nb[z]] == m_waSurfCoverMol[z3]) &&
									(g_waAtomElement[nb[z]] == m_waSurfCoverElem[z3]) &&
									(g_waAtomMolNumber[nb[z]] == m_waSurfCoverAtom[z3]))
								{
									tfa = (CxDoubleArray*)m_oaSurfCoverData[z3];

									if (!tb[z3])
									{
										tb[z3] = true;
										tfa->Add(fa[z]);
									} else (*tfa)[tfa->GetSize()-1] += fa[z];

									atot += fa[z];

									goto _done;
								}
							}
							eprintf("Error.\n");
_done:;
						}
					}
				}
			}
		} while (vl.inc());
	}

	for (z=0;z<m_waSurfCoverAtom.GetSize();z++)
	{
		if (!tb[z])
			((CxDoubleArray*)m_oaSurfCoverData[z])->Add(0);

		if (atot > 0)
			((CxDoubleArray*)m_oaSurfCoverData[z])->GetAt(((CxDoubleArray*)m_oaSurfCoverData[z])->GetSize()-1) /= atot;
	}

	delete m_pContainer;
	delete[] tb;
}


void CVoroWrapper::FinishSurfCover(const char *s)
{
	FILE *a;
	int z, z2;
	CMolecule *m;
	double t;

	a = OpenFileWrite(s,true);

	mfprintf(a,"# Step;  ");

	for (z=0;z<m_waSurfCoverMol.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[m_waSurfCoverMol[z]];
		mfprintf(a,"%s %s%d",m->m_sName,((CAtom*)g_oaAtoms[m->m_baAtomIndex[m_waSurfCoverElem[z]]])->m_sName,m_waSurfCoverAtom[z]+1);
		if (z < m_waSurfCoverMol.GetSize()-1)
			mfprintf(a,";  ");
	}
	mfprintf(a,"\n");

	mfprintf(a,"# Averages;  ");
	for (z2=0;z2<m_waSurfCoverMol.GetSize();z2++)
	{
		t = 0;
		for (z=0;z<((CxDoubleArray*)m_oaSurfCoverData[z2])->GetSize();z++)
			t += ((CxDoubleArray*)m_oaSurfCoverData[z2])->GetAt(z);

		if (t > 0)
			t /= ((CxDoubleArray*)m_oaSurfCoverData[z2])->GetSize();

		mfprintf(a,"%f",t*100.0f);
		if (z2 < m_waSurfCoverMol.GetSize()-1)
			mfprintf(a,";  ");
	}
	mfprintf(a,"\n");

	for (z=0;z<((CxDoubleArray*)m_oaSurfCoverData[0])->GetSize();z++)
	{
		mfprintf(a,"%d;  ",z);
		for (z2=0;z2<m_waSurfCoverMol.GetSize();z2++)
		{
			mfprintf(a,"%f",((CxDoubleArray*)m_oaSurfCoverData[z2])->GetAt(z)*100.0f);
			if (z2 < m_waSurfCoverMol.GetSize()-1)
				mfprintf(a,";  ");
		}
		mfprintf(a,"\n");
	}

	fclose(a);
}




