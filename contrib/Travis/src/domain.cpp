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


#include "domain.h"
#include "tools.h"
#include "timestep.h"
#include "vorowrapper.h"
#include "globalvar.h"


CDomain::CDomain()
{
	Reset();
}


CDomain::~CDomain()
{
}


void CDomain::Reset()
{
	m_iaCells.RemoveAll_KeepSize();
	m_iaNeighbors.RemoveAll_KeepSize();
	m_bActive = true;

	m_iFaces = 0;
	m_fSurfaceArea = 0;
	m_fVolume = 0;
	m_fAVRatio = 0;
}


void CDomain::Assimilate(CDomain *dom)
{
	int z;

	for (z=0;z<dom->m_iaCells.GetSize();z++)
		m_iaCells.Add(dom->m_iaCells[z]);

	dom->m_bActive = false;

	m_fVolume += dom->m_fVolume;
	m_fSurfaceArea += dom->m_fSurfaceArea;
	m_iFaces += dom->m_iFaces;
	m_fAVRatio = 10.6347231054330961*m_fVolume/pow(m_fSurfaceArea,1.5);
}


CDomainAnalysis::CDomainAnalysis()
{
}


CDomainAnalysis::~CDomainAnalysis()
{
}


void CDomainAnalysis::Parse(int i)
{
	int z, z2, z3, z4;
	CAtomGroup *ag;
	CMolecule *m;
	CSingleMolecule *sm;
//	char buf[256];
	CxString buf;

	mprintf(YELLOW,"    >>> Domain Analysis #%d >>>\n\n",i+1);

	mprintf("    First, you have to specify the base population (i.e., all atoms\n");
	mprintf("    that will have a Voronoi cell in the analysis).\n\n");
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];

_ag:
		if (AskYesNo("      Include some atoms of molecule type %d (%s) in base population (y/n)? [yes] ",true,z+1,m->m_sName))
		{
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			AskString("        Which atoms to include (e.g. C1,C5-7, *=all)? [all] ",&buf,"*");
			if (!ag->ParseAtoms(m,buf))
				goto _ag;
			m_oaBasePopulation.Add(ag);
		} else
			m_oaBasePopulation.Add(NULL);
	}

	mprintf("\n    Processing selection...\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		ag = (CAtomGroup*)m_oaBasePopulation[z];
		if (ag == NULL)
			continue;
		m = ag->m_pMolecule;
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			for (z3=0;z3<ag->m_oaAtoms.GetSize();z3++)
			{
				for (z4=0;z4<((CxIntArray*)ag->m_oaAtoms[z3])->GetSize();z4++)
					m_iaBasePopulation.Add(((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z3]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z3])->GetAt(z4)));
			}
		}
	}

	mprintf("      Added %d atoms to base population.\n\n",m_iaBasePopulation.GetSize());

	mprintf("    Now you have to specify which atoms from the base population\n");
	mprintf("    belong to your domain of interest.\n\n");
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];

_ag2:
		if (AskYesNo("      Include some atoms of molecule type %d (%s) in domain (y/n)? [no] ",false,z+1,m->m_sName))
		{
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			AskString("        Which atoms to include (e.g. C1,C5-7, *=all)? [all] ",&buf,"*");
			if (!ag->ParseAtoms(m,buf))
				goto _ag2;
			m_oaDomainSet.Add(ag);
		} else
			m_oaDomainSet.Add(NULL);
	}

	mprintf("\n    Processing selection...\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		ag = (CAtomGroup*)m_oaDomainSet[z];
		if (ag == NULL)
			continue;
		m = ag->m_pMolecule;
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			for (z3=0;z3<ag->m_oaAtoms.GetSize();z3++)
			{
				for (z4=0;z4<((CxIntArray*)ag->m_oaAtoms[z3])->GetSize();z4++)
				{
					m_iaDomainSet.Add(((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z3]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z3])->GetAt(z4)));
				}
			}
		}
	}

	mprintf("      Added %d atoms to domain.\n\n",m_iaDomainSet.GetSize());

	mprintf("    Checking if domain is a subset of base population...\n");

	for (z=0;z<m_iaDomainSet.GetSize();z++)
	{
		for (z2=0;z2<m_iaBasePopulation.GetSize();z2++)
		{
			if (m_iaDomainSet[z] == m_iaBasePopulation[z2])
				goto _done;
		}
		eprintf("Error: Some atoms of domain are not included in base population.\n");
		abort();
_done:;
	}

	m_iaBaseInDomain.SetSize(m_iaBasePopulation.GetSize());
	for (z=0;z<m_iaBasePopulation.GetSize();z++)
	{
		for (z2=0;z2<m_iaDomainSet.GetSize();z2++)
		{
			if (m_iaBasePopulation[z] == m_iaDomainSet[z2])
			{
				m_iaBaseInDomain[z] = 1;
				goto _done2;
			}
		}
		m_iaBaseInDomain[z] = 0;
_done2:;
	}

	m_oaDomains.SetSize(m_iaBasePopulation.GetSize());

	for (z=0;z<m_iaBasePopulation.GetSize();z++)
	{
		if (m_iaBaseInDomain[z])
			m_oaDomains[z] = new CDomain();
		else
			m_oaDomains[z] = NULL;
	}

	mprintf("\n");
	m_bWriteHistograms = AskYesNo("    Write histograms of domain properties (y) or only min./max./average/std.dev.? [yes] ",true);
	if (m_bWriteHistograms)
	{
		mprintf("\n");
		mprintf(WHITE,"    Domain Surface Area Histogram\n");
		m_pHistoSurface = new CDF();
		m_pHistoSurface->m_fMinVal = 0;
		m_pHistoSurface->m_fMaxVal = AskFloat("      Enter max. domain surface area value [A^2]: [10000.0] ",10000);
		m_pHistoSurface->m_iResolution = AskUnsignedInteger("      Enter histogram resolution: [100] ",100);
		m_pHistoSurface->SetLabelX("Domain Surface Area");
		m_pHistoSurface->SetLabelY("Occurrence");
		m_pHistoSurface->Create();

		mprintf("\n");
		mprintf(WHITE,"    Domain Volume Histogram\n");
		m_pHistoVolume = new CDF();
		m_pHistoVolume->m_fMinVal = 0;
		m_pHistoVolume->m_fMaxVal = AskFloat("      Enter max. domain volume value [A^3]: [10000.0] ",10000);
		m_pHistoVolume->m_iResolution = AskUnsignedInteger("      Enter histogram resolution: [100] ",100);
		m_pHistoVolume->SetLabelX("Domain Volume");
		m_pHistoVolume->SetLabelY("Occurrence");
		m_pHistoVolume->Create();

		mprintf("\n");
		mprintf(WHITE,"    Domain Isoperimetric Quotient Histogram\n");
		m_pHistoAVRatio = new CDF();
		m_pHistoAVRatio->m_fMinVal = 0;
		m_pHistoAVRatio->m_fMaxVal = AskFloat("      Enter max. domain isoperimetric quotient value: [1.0] ",1);
		m_pHistoAVRatio->m_iResolution = AskUnsignedInteger("      Enter histogram resolution: [100] ",100);
		m_pHistoAVRatio->SetLabelX("Domain Isoperimetric Quotient");
		m_pHistoAVRatio->SetLabelY("Occurrence");
		m_pHistoAVRatio->Create();

		mprintf("\n");
		mprintf(WHITE,"    Domain Face Count Histogram\n");
		m_pHistoFaceCount = new CDF();
		m_pHistoFaceCount->m_fMinVal = 0;
		m_pHistoFaceCount->m_fMaxVal = AskUnsignedInteger("      Enter max. domain face count: [10000] ",10000);
		m_pHistoFaceCount->m_iResolution = AskUnsignedInteger("      Enter histogram resolution: [100] ",100);
		m_pHistoFaceCount->SetLabelX("Domain Face Count");
		m_pHistoFaceCount->SetLabelY("Occurrence");
		m_pHistoFaceCount->Create();

		mprintf("\n");
		mprintf(WHITE,"    Domain Cell Count Histogram\n");
		m_pHistoCellCount = new CDF();
		m_pHistoCellCount->m_fMinVal = 0;
		m_pHistoCellCount->m_fMaxVal = AskUnsignedInteger("      Enter max. domain cell count: [%d] ",m_iaDomainSet.GetSize(),m_iaDomainSet.GetSize());
		m_pHistoCellCount->m_iResolution = AskUnsignedInteger("      Enter histogram resolution: [100] ",100);
		m_pHistoCellCount->SetLabelX("Domain Cell Count");
		m_pHistoCellCount->SetLabelY("Occurrence");
		m_pHistoCellCount->Create();

		mprintf("\n");
		mprintf(WHITE,"    Domain Count Histogram\n");
		m_pHistoDomainCount = new CDF();
		m_pHistoDomainCount->m_fMinVal = 0;
		m_pHistoDomainCount->m_fMaxVal = AskUnsignedInteger("      Enter max. domain count: [%d] ",m_iaDomainSet.GetSize(),m_iaDomainSet.GetSize());
		m_pHistoDomainCount->m_iResolution = AskUnsignedInteger("      Enter histogram resolution: [100] ",100);
		m_pHistoDomainCount->SetLabelX("Domain Count");
		m_pHistoDomainCount->SetLabelY("Occurrence");
		m_pHistoDomainCount->Create();
	}

	mprintf(YELLOW,"\n    <<< End of Domain Analysis #%d <<<\n\n",i+1);
}


void CDomainAnalysis::ProcessStep(CTimeStep *ts)
{
	voronoicell_neighbor c;
	container_periodic_poly *con;
	int ijk, q, z, z2, faces, id;
	double tf;
	vector<int> nb/*, fo*/;
	vector<double> fa;
	CDomain *domain;
	CxIntArray stack;
	double sumvol, sqvol, sumsurf, sqsurf, sumav, sqav, sumcell, sqcell, sumface, sqface;
	double minvol, maxvol, minsurf, maxsurf, minav, maxav, mincell, maxcell, minface, maxface;

	try { con = new container_periodic_poly(g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0,g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ,g_iVoroMemory); } catch(...) { con = NULL; }
	if (con == NULL) NewException((double)sizeof(container_periodic_poly),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_iaBasePopulation.GetSize();z++)
		con->put(z,ts->m_vaCoords[m_iaBasePopulation[z]][0]/1000.0,ts->m_vaCoords[m_iaBasePopulation[z]][1]/1000.0,ts->m_vaCoords[m_iaBasePopulation[z]][2]/1000.0,g_faVoronoiRadii[m_iaBasePopulation[z]]/1000.0);

	c_loop_all_periodic vl(*con);

	for (z=0;z<m_oaDomains.GetSize();z++)
		if (m_oaDomains[z] != NULL)
			((CDomain*)m_oaDomains[z])->Reset();

	if (vl.start()) 
	{
		do 
		{
			if (con->compute_cell(c,vl))
			{
				ijk=vl.ijk;
				q=vl.q;

				id = con->id[ijk][q];

				if (!m_iaBaseInDomain[id])
					continue;

				c.neighbors(nb);
			//	c.face_vertices(fv);
				c.face_areas(fa);
			//	c.face_orders(fo);

				faces = c.number_of_faces();

				domain = (CDomain*)m_oaDomains[id];

				if (domain == NULL)
					mprintf("Error.\n");

				domain->m_iaCells.Add(id);

				domain->m_fSurfaceArea = 0;
				domain->m_iFaces = 0;

				for (z=0;z<faces;z++)
				{
					if (m_iaBaseInDomain[nb[z]])
					{ // Face leads to other domain member
						domain->m_iaNeighbors.Add(nb[z]);
					} else
					{ // Face leads out of domain --> outer surface
						domain->m_fSurfaceArea += fa[z]*100.0;
						domain->m_iFaces++;
					}
				}

				domain->m_fVolume = c.volume()*1000.0;
				domain->m_fAVRatio = 10.6347231054330961*domain->m_fVolume/pow(domain->m_fSurfaceArea,1.5);
			}
		} while (vl.inc());
	}

	delete con;

	for (z=0;z<m_oaDomains.GetSize();z++)
	{
		if (m_oaDomains[z] == NULL)
			continue;

		if (!((CDomain*)m_oaDomains[z])->m_bActive)
			continue;

	//	mprintf("Starting Recursion for %d (%d)...\n",z,((CDomain*)m_oaDomains[z])->m_iaNeighbors.GetSize());

		stack.RemoveAll_KeepSize();
		REC_FuseDomain(z,z,&stack);
	}

	z2 = 0;
	sumvol = 0;
	sqvol = 0;
	minvol = 1.0e30;
	maxvol = -1.0e30;
	sumsurf = 0;
	sqsurf = 0;
	minsurf = 1.0e30;
	maxsurf = -1.0e30;
	sumav = 0;
	sqav = 0;
	minav = 1.0e30;
	maxav = -1.0e30;
	sumface = 0;
	sqface = 0;
	minface = 1.0e30;
	maxface = -1.0e30;
	sumcell = 0;
	sqcell = 0;
	mincell = 1.0e30;
	maxcell = -1.0e30;

	for (z=0;z<m_oaDomains.GetSize();z++)
	{
		if (m_oaDomains[z] == NULL)
			continue;

		domain = (CDomain*)m_oaDomains[z];

		if (!domain->m_bActive)
			continue;

		tf = domain->m_fVolume;
		sumvol += tf;
		sqvol += tf*tf;
		if (tf < minvol)
			minvol = tf;
		if (tf > maxvol)
			maxvol = tf;

		tf = domain->m_fSurfaceArea;
		sumsurf += tf;
		sqsurf += tf*tf;
		if (tf < minsurf)
			minsurf = tf;
		if (tf > maxsurf)
			maxsurf = tf;

		tf = domain->m_fAVRatio;
		sumav += tf;
		sqav += tf*tf;
		if (tf < minav)
			minav = tf;
		if (tf > maxav)
			maxav = tf;

		tf = domain->m_iFaces;
		sumface += tf;
		sqface += tf*tf;
		if (tf < minface)
			minface = tf;
		if (tf > maxface)
			maxface = tf;

		tf = domain->m_iaCells.GetSize();
		sumcell += tf;
		sqcell += tf*tf;
		if (tf < mincell)
			mincell = tf;
		if (tf > maxcell)
			maxcell = tf;

		if (m_bWriteHistograms)
		{
			m_pHistoVolume->AddToBin(domain->m_fVolume);
			m_pHistoSurface->AddToBin(domain->m_fSurfaceArea);
			m_pHistoAVRatio->AddToBin(domain->m_fAVRatio);
			m_pHistoFaceCount->AddToBin(domain->m_iFaces);
			m_pHistoCellCount->AddToBin(domain->m_iaCells.GetSize());
		}

		z2++;
	}

	if (m_bWriteHistograms)
		m_pHistoDomainCount->AddToBin(z2);

//	mprintf("%d domains active.\n",z2);

	sumvol /= z2;
	m_faVolumeAv.Add(sumvol);
	m_faVolumeSD.Add(sqrt(sqvol/z2 - sumvol*sumvol));
	m_faVolumeMin.Add(minvol);
	m_faVolumeMax.Add(maxvol);

	sumsurf /= z2;
	m_faSurfaceAv.Add(sumsurf);
	m_faSurfaceSD.Add(sqrt(sqsurf/z2 - sumsurf*sumsurf));
	m_faSurfaceMin.Add(minsurf);
	m_faSurfaceMax.Add(maxsurf);

	sumav /= z2;
	m_faAVRatioAv.Add(sumav);
	m_faAVRatioSD.Add(sqrt(sqav/z2 - sumav*sumav));
	m_faAVRatioMin.Add(minav);
	m_faAVRatioMax.Add(maxav);

	sumface /= z2;
	m_faFaceCountAv.Add(sumface);
	m_faFaceCountSD.Add(sqrt(sqface/z2 - sumface*sumface));
	m_faFaceCountMin.Add(minface);
	m_faFaceCountMax.Add(maxface);

	sumcell /= z2;
	m_faCellCountAv.Add(sumcell);
	m_faCellCountSD.Add(sqrt(sqcell/z2 - sumcell*sumcell));
	m_faCellCountMin.Add(mincell);
	m_faCellCountMax.Add(maxcell);

	m_iaDomainCount.Add(z2);

	m_iaStep.Add(g_iSteps);
}


void CDomainAnalysis::REC_FuseDomain(int basedom, int dom, CxIntArray *stack)
{
	int i, z, z2;
	CDomain *pdom;

	stack->Add(dom);

	pdom = (CDomain*)m_oaDomains[dom];

	if (dom != basedom)
	{
	//	mprintf("    Assimilating %d into %d.\n",dom,basedom);
		((CDomain*)m_oaDomains[basedom])->Assimilate(pdom);
	}

	for (z=0;z<pdom->m_iaNeighbors.GetSize();z++)
	{
		i = pdom->m_iaNeighbors[z];
		if (((CDomain*)m_oaDomains[i])->m_bActive)
		{
			for (z2=0;z2<stack->GetSize();z2++)
				if (stack->GetAt(z2) == i)
					goto _next;
			REC_FuseDomain(basedom,i,stack);
		}
_next:;
	}

	stack->Pop_KeepSize();
}


void CDomainAnalysis::Finish(int i)
{
//	char buf[1024];
	CxString buf, buf2;
	double tf;
	int z;
	FILE *a;

//	sprintf(buf,"domain");
	buf.sprintf("domain");

	for (z=0;z<m_oaDomainSet.GetSize();z++)
	{
		if (m_oaDomainSet[z] == NULL)
			continue;
//		strcat(buf,"_");
//		strcat(buf,((CAtomGroup*)m_oaDomainSet[z])->m_pMolecule->m_sName);
//		strcat(buf,"_");
//		strcat(buf,((CAtomGroup*)m_oaDomainSet[z])->m_sName);
		buf.strcat("_");
		buf.strcat(((CAtomGroup*)m_oaDomainSet[z])->m_pMolecule->m_sName);
		buf.strcat("_");
		buf.strcat(((CAtomGroup*)m_oaDomainSet[z])->m_sName);
	}
//	strcat(buf,".csv");

	buf2 = buf;
	buf2.strcat("_statistics.csv");

	mprintf("      Analysis %d: Writing results to %s ...\n",i+1,(const char*)buf);

	a = OpenFileWrite(buf2,true);

	fprintf(a,"#Step;  Domain_count;  Volume_(A^3)_Avg;  Min;  Max;  SD;  Surface_Area_(A^2)_Avg;  Min;  Max;  SD;  Isoperimetric_Quotient_Avg;  Min;  Max;  SD;  Face_Count_Avg;  Min;  Max;  SD;  Cell_Count_Avg;  Min;  Max;  SD\n");

	fprintf(a,"#Simulation_Average;  ");

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_iaDomainCount[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faVolumeAv[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faVolumeMin[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faVolumeMax[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faVolumeSD[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faSurfaceAv[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faSurfaceMin[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faSurfaceMax[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faSurfaceSD[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faAVRatioAv[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faAVRatioMin[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);
	
	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faAVRatioMax[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faAVRatioSD[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faFaceCountAv[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faFaceCountMin[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faFaceCountMax[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faFaceCountSD[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faCellCountAv[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);

	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faCellCountMin[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);
	
	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faCellCountMax[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f;  ",tf);
	
	tf = 0;
	for (z=0;z<m_iaStep.GetSize();z++)
		tf += m_faCellCountSD[z];
	tf /= m_iaStep.GetSize();
	fprintf(a,"%.10f\n",tf);

	for (z=0;z<m_iaDomainCount.GetSize();z++)
		fprintf(a,"%d;  %d;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.10f;  %.0f;  %.0f;  %.10f;  %.10f;  %.0f;  %.0f;  %.10f\n",m_iaStep[z],m_iaDomainCount[z],m_faVolumeAv[z],m_faVolumeMin[z],m_faVolumeMax[z],m_faVolumeSD[z],m_faSurfaceAv[z],m_faSurfaceMin[z],m_faSurfaceMax[z],m_faSurfaceSD[z],m_faAVRatioAv[z],m_faAVRatioMin[z],m_faAVRatioMax[z],m_faAVRatioSD[z],m_faFaceCountAv[z],m_faFaceCountMin[z],m_faFaceCountMax[z],m_faFaceCountSD[z],m_faCellCountAv[z],m_faCellCountMin[z],m_faCellCountMax[z],m_faCellCountSD[z]);

	fclose(a);

	if (m_bWriteHistograms)
	{
		mprintf("      Saving Domain Surface Area Histogram as %s_histo_surface.csv.\n",(const char*)buf);
		mprintf("        (%.0f bin entries, %.0f skipped)\n\n",m_pHistoSurface->m_fBinEntries,m_pHistoSurface->m_fSkipEntries);
		m_pHistoSurface->NormBinSum(100.0);
		m_pHistoSurface->Write("",buf,"_histo_surface.csv",true);

		mprintf("      Saving Domain Volume Histogram as %s_histo_volume.csv.\n",(const char*)buf);
		mprintf("        (%.0f bin entries, %.0f skipped)\n\n",m_pHistoVolume->m_fBinEntries,m_pHistoVolume->m_fSkipEntries);
		m_pHistoVolume->NormBinSum(100.0);
		m_pHistoVolume->Write("",buf,"_histo_volume.csv",true);

		mprintf("      Saving Domain Isoperimetric Quotient Histogram as %s_histo_avratio.csv.\n",(const char*)buf);
		mprintf("        (%.0f bin entries, %.0f skipped)\n\n",m_pHistoAVRatio->m_fBinEntries,m_pHistoAVRatio->m_fSkipEntries);
		m_pHistoAVRatio->NormBinSum(100.0);
		m_pHistoAVRatio->Write("",buf,"_histo_avratio.csv",true);

		mprintf("      Saving Domain Face Count Histogram as %s_histo_faces.csv.\n",(const char*)buf);
		mprintf("        (%.0f bin entries, %.0f skipped)\n\n",m_pHistoFaceCount->m_fBinEntries,m_pHistoFaceCount->m_fSkipEntries);
		m_pHistoFaceCount->NormBinSum(100.0);
		m_pHistoFaceCount->Write("",buf,"_histo_faces.csv",true);

		mprintf("      Saving Domain Cell Count Histogram as %s_histo_cells.csv.\n",(const char*)buf);
		mprintf("        (%.0f bin entries, %.0f skipped)\n\n",m_pHistoCellCount->m_fBinEntries,m_pHistoCellCount->m_fSkipEntries);
		m_pHistoCellCount->NormBinSum(100.0);
		m_pHistoCellCount->Write("",buf,"_histo_cells.csv",true);

		mprintf("      Saving Domain Count Histogram as %s_histo_domaincount.csv.\n",(const char*)buf);
		mprintf("        (%.0f bin entries, %.0f skipped)\n\n",m_pHistoDomainCount->m_fBinEntries,m_pHistoDomainCount->m_fSkipEntries);
		m_pHistoDomainCount->NormBinSum(100.0);
		m_pHistoDomainCount->Write("",buf,"_histo_domaincount.csv",true);
	}
}


CDomainEngine::CDomainEngine()
{
}


CDomainEngine::~CDomainEngine()
{
}


void CDomainEngine::Parse()
{
	CTimeStep *ts;
	CDomainAnalysis *tempdoma;

	try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
	if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	mprintf("    Initializing Voronoi tesselation...\n");
	g_pVoroWrapper->Init();
	mprintf("\n");
	mprintf("*** Voro: Box density is %f particles / Angstrom^3.\n",g_pVoroWrapper->m_fBoxDens);
	mprintf("*** Voro: Using %d x %d x %d blocks.\n",g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ);
	try { ts = new CTimeStep(); } catch(...) { ts = NULL; }
	if (ts == NULL) NewException((double)sizeof(CTimeStep),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	ts->CopyFrom(&g_TimeStep);
	ts->FoldAtomsPositive();
	g_pVoroWrapper->Dump("voro.txt",ts);
	mprintf("\n");
	mprintf("    Voro++: Using cell memory for %d particles.\n\n",g_iVoroMemory);

	do {
		try { tempdoma = new CDomainAnalysis(); } catch(...) { tempdoma = NULL; }
		if (tempdoma == NULL) NewException(sizeof(CDomainAnalysis),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		tempdoma->Parse(m_oaAnalyses.GetSize());
		m_oaAnalyses.Add(tempdoma);
	} while (AskYesNo("    Add another domain analysis (y/n)? [no] ",false));
}


void CDomainEngine::ProcessStep(CTimeStep *ts)
{
	int z;

	for (z=0;z<m_oaAnalyses.GetSize();z++)
		((CDomainAnalysis*)m_oaAnalyses[z])->ProcessStep(ts);
}


void CDomainEngine::Finish()
{
	int z;

	mprintf(YELLOW,"    >>> Finishing Domain Analysis >>>\n");

	for (z=0;z<m_oaAnalyses.GetSize();z++)
		((CDomainAnalysis*)m_oaAnalyses[z])->Finish(z);

	mprintf(YELLOW,"    <<< Finishing Domain Analysis done <<<\n");
}

