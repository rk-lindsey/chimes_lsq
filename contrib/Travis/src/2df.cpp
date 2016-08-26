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

#include "2df.h"
#include "travis.h"


C2DF::C2DF()
{
	m_sLabelX = NULL;
	m_sLabelY = NULL;
	m_sLabelZ = NULL;
	m_pBin = NULL;
	m_pStepsY = NULL;
	m_fCountX = NULL;
	m_fCountY = NULL;
	m_iHistogramRes = 0;

	m_iPlotType = 1;
	m_iSmoothGrade = 1;
	m_iInterpolationOrder = 2;
	m_fPlotExp = g_pDatabase->GetFloat("/PLOT2D/DEFAULTS/PLOT_EXP");
	m_iExpLegend = 1;
	m_iColorScale = 1;
	m_iGPInterpolation = 5;
	m_bContourLines = true;
	m_fAspectRatio = 1.0;

	m_iPlotPixel = g_pDatabase->GetInt("/PLOT2D/DEFAULTS/IMAGE_RES");

	m_oaCircles.SetName("C2DF::m_oaCircles");
}


C2DF::~C2DF()
{
	if (m_sLabelX != NULL)
	{
		delete[] m_sLabelX;
		m_sLabelX = NULL;
	}
	if (m_sLabelY != NULL)
	{
		delete[] m_sLabelY;
		m_sLabelY = NULL;
	}
	if (m_sLabelZ != NULL)
	{
		delete[] m_sLabelZ;
		m_sLabelZ = NULL;
	}
	if (m_pBin != NULL)
	{
		delete[] m_pBin;
		m_pBin = NULL;
	}
	if (m_pStepsY != NULL)
	{
		delete[] m_pStepsY;
		m_pStepsY = NULL;
	}
	if (m_fCountX != NULL)
	{
		delete[] m_fCountX;
		m_fCountX = NULL;
	}
	if (m_fCountY != NULL)
	{
		delete[] m_fCountY;
		m_fCountY = NULL;
	}
}


void C2DF::CorrectAngle(int channel)
{
	BTIN;
	int z, z2;
	double d;

	if (channel == 0)
	{
		for (z=0;z<m_iRes[0];z++)
		{
			d = cos((m_fMinVal[0]+(double)z*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0])*Pi/180.0) - cos((m_fMinVal[0]+(double)(z+1)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0])*Pi/180.0);
			for (z2=0;z2<m_iRes[1];z2++)
			{
				m_pBin[z2*m_iRes[0]+z] /= d;
			}
		}
	} else
	{
		for (z2=0;z2<m_iRes[1];z2++)
		{
			d = cos((m_fMinVal[1]+(double)z2*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1])*Pi/180.0) - cos((m_fMinVal[1]+(double)(z2+1)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1])*Pi/180.0);
			for (z=0;z<m_iRes[0];z++)
			{
				m_pBin[z2*m_iRes[0]+z] /= d;
			}
		}
	}
	BTOUT;
}


void C2DF::UnCorrectAngle(int channel)
{
	BTIN;
	int z, z2;

	if (channel == 0)
	{
		for (z=0;z<m_iRes[0];z++)
			for (z2=0;z2<m_iRes[1];z2++)
				m_pBin[z2*m_iRes[0]+z] *= cos((m_fMinVal[0]+(double)z*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0])*Pi/180.0) - cos((m_fMinVal[0]+(double)(z+1)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0])*Pi/180.0);
	} else
	{
		for (z=0;z<m_iRes[0];z++)
			for (z2=0;z2<m_iRes[1];z2++)
				m_pBin[z2*m_iRes[0]+z] *= cos((m_fMinVal[1]+(double)z2*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1])*Pi/180.0) - cos((m_fMinVal[1]+(double)(z2+1)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1])*Pi/180.0);
	}
	BTOUT;
}


void C2DF::CorrectRadialDist(int channel)
{
	BTIN;
	int z, z2;
	double d;

	if (channel == 0)
	{
		for (z=0;z<m_iRes[0];z++)
		{
			d = pow(m_fMinVal[0]+(z+1.0)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),3) - pow(m_fMinVal[0]+((double)z)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),3);
			for (z2=0;z2<m_iRes[1];z2++)
			{
				m_pBin[z2*m_iRes[0]+z] /= d;
			}
		}
	} else
	{
		for (z2=0;z2<m_iRes[1];z2++)
		{
			d = pow(m_fMinVal[1]+(z2+1.0)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),3) - pow(m_fMinVal[1]+((double)z2)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),3);
			for (z=0;z<m_iRes[0];z++)
			{
				m_pBin[z2*m_iRes[0]+z] /= d;
			}
		}
	}
	BTOUT;
}


void C2DF::CorrectLiRadialDist(int channel)
{
	BTIN;
	int z, z2;
	double d;

	if (channel == 0)
	{
		for (z=0;z<m_iRes[0];z++)
		{
			d = pow(m_fMinVal[0]+(z+1.0)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),2) - pow(m_fMinVal[0]+((double)z)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),2);
			for (z2=0;z2<m_iRes[1];z2++)
			{
				m_pBin[z2*m_iRes[0]+z] /= d;
			}
		}
	} else
	{
		for (z2=0;z2<m_iRes[1];z2++)
		{
			d = pow(m_fMinVal[1]+(z2+1.0)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),2) - pow(m_fMinVal[1]+((double)z2)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),2);
			for (z=0;z<m_iRes[0];z++)
			{
				m_pBin[z2*m_iRes[0]+z] /= d;
			}
		}
	}
	BTOUT;
}


void C2DF::UnCorrectRadialDist(int channel)
{
	BTIN;
	int z, z2;

	if (channel == 0)
	{
		for (z=0;z<m_iRes[0];z++)
			for (z2=0;z2<m_iRes[1];z2++)
				m_pBin[z2*m_iRes[0]+z] *= (double)(pow(m_fMinVal[0]+(z+1.0)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),3) - pow(m_fMinVal[0]+((double)z)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),3));
	} else
	{
		for (z=0;z<m_iRes[0];z++)
			for (z2=0;z2<m_iRes[1];z2++)
				m_pBin[z2*m_iRes[0]+z] *= (double)(pow(m_fMinVal[1]+(z2+1.0)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),3) - pow(m_fMinVal[1]+((double)z2)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),3));
	}
	BTOUT;
}


void C2DF::AddToBin(double x, double y)
{
	BXIN;
	double rx, ry;
	int ix, iy;

	if ((x < m_fMinVal[0]) || (y < m_fMinVal[1]) || (x > m_fMaxVal[0]) || (y > m_fMaxVal[1]))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;

	rx = (x-m_fMinVal[0])*m_fFac[0] - 0.5;
	ix = (int)floor(rx);
	if (ix < 0)
	{
		ix = 0;
		rx = 0;
	} else if (ix > m_iRes[0]-2)
	{
		ix = m_iRes[0]-2;
		rx = 1.0;
	} else
		rx -= ix;

	ry = (y-m_fMinVal[1])*m_fFac[1] - 0.5;
	iy = (int)floor(ry);
	if (iy < 0)
	{
		iy = 0;
		ry = 0;
	} else if (iy > m_iRes[1]-2)
	{
		iy = m_iRes[1]-2;
		ry = 1.0;
	} else
		ry -= iy;

	m_pBin[ iy    * m_iRes[0] + ix    ] += (1.0-rx) * (1.0-ry);
	m_pBin[ iy    * m_iRes[0] + ix + 1] +=      rx  * (1.0-ry);
	m_pBin[(iy+1) * m_iRes[0] + ix    ] += (1.0-rx) *      ry;
	m_pBin[(iy+1) * m_iRes[0] + ix + 1] +=      rx  *      ry;

	m_fCountX[ix  ] += (1.0-rx);
	m_fCountX[ix+1] +=      rx;
	m_fCountY[iy  ] += (1.0-ry);
	m_fCountY[iy+1] +=      ry;
	BXOUT;
}


void C2DF::AddToBin(double x, double y, double val)
{ 
	BXIN;
	double rx, ry;
	int ix, iy;

	if ((x < m_fMinVal[0]) || (y < m_fMinVal[1]) || (x > m_fMaxVal[0]) || (y > m_fMaxVal[1]))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;

	rx = (x-m_fMinVal[0])*m_fFac[0] - 0.5;
	ix = (int)floor(rx);
	if (ix < 0)
	{
		ix = 0;
		rx = 0;
	} else if (ix > m_iRes[0]-2)
	{
		ix = m_iRes[0]-2;
		rx = 1.0;
	} else
		rx -= ix;

	ry = (y-m_fMinVal[1])*m_fFac[1] - 0.5;
	iy = (int)floor(ry);
	if (iy < 0)
	{
		iy = 0;
		ry = 0;
	} else if (iy > m_iRes[1]-2)
	{
		iy = m_iRes[1]-2;
		ry = 1.0;
	} else
		ry -= iy;

	m_pBin[ iy    * m_iRes[0] + ix    ] += (1.0-rx) * (1.0-ry) * val;
	m_pBin[ iy    * m_iRes[0] + ix + 1] +=      rx  * (1.0-ry) * val;
	m_pBin[(iy+1) * m_iRes[0] + ix    ] += (1.0-rx) *      ry  * val;
	m_pBin[(iy+1) * m_iRes[0] + ix + 1] +=      rx  *      ry  * val;

	m_fCountX[ix  ] += (1.0-rx);
	m_fCountX[ix+1] +=      rx;
	m_fCountY[iy  ] += (1.0-ry);
	m_fCountY[iy+1] +=      ry;
	BXOUT;
}


void C2DF::AddToBin_IntX(int x, double y, double val)
{ 
	BXIN;
	double ry;
	int iy;

	if ((y < m_fMinVal[1]) || (y > m_fMaxVal[1]))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;

	ry = (y-m_fMinVal[1])*m_fFac[1] - 0.5;
	iy = (int)floor(ry);
	if (iy < 0)
	{
		iy = 0;
		ry = 0;
	} else if (iy > m_iRes[1]-2)
	{
		iy = m_iRes[1]-2;
		ry = 1.0;
	} else
		ry -= iy;

	m_pBin[ iy    * m_iRes[0] + x] += (1.0-ry) * val;
	m_pBin[(iy+1) * m_iRes[0] + x] +=      ry  * val;

	m_fCountX[x]++;
	m_fCountY[iy  ] += (1.0-ry);
	m_fCountY[iy+1] +=      ry;
	BXOUT;
}


void C2DF::AddToBin_IntY(double x, int y, double val)
{ 
	BXIN;
	double rx;
	int ix;

	if ((x < m_fMinVal[0]) || (x > m_fMaxVal[0]))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;

	rx = (x-m_fMinVal[0])*m_fFac[0] - 0.5;
	ix = (int)floor(rx);
	if (ix < 0)
	{
		ix = 0;
		rx = 0;
	} else if (ix > m_iRes[0]-2)
	{
		ix = m_iRes[0]-2;
		rx = 1.0;
	} else
		rx -= ix;

	m_pBin[ ix    + y * m_iRes[0]] += (1.0-rx) * val;
	m_pBin[(ix+1) + y * m_iRes[0]] +=      rx  * val;

	m_fCountY[y]++;
	m_fCountX[ix  ] += (1.0-rx);
	m_fCountX[ix+1] +=      rx;
	BXOUT;
}



/*void C2DF::AddToSingleBin(double x, double y, double val)
{ 
	BXIN;
	double rx, ry;
	int ix, iy;

	if ((x < m_fMinVal[0]) || (y < m_fMinVal[1]) || (x > m_fMaxVal[0]) || (y > m_fMaxVal[1]))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;
	rx = ((x-m_fMinVal[0])/(m_fMaxVal[0]-m_fMinVal[0]))*((double)m_iRes[0]-1);
	ry = ((y-m_fMinVal[1])/(m_fMaxVal[1]-m_fMinVal[1]))*((double)m_iRes[1]-1);
	ix = (int)(rx+0.5);
	iy = (int)(ry+0.5);
	m_pBin[iy*(m_iRes[0]+1) + ix] += val;
	BXOUT;
}*/


//DEL void C2DF::CorrectRadialDist(bool xdim)
//DEL {
//DEL 	BTIN;
//DEL 	int x, y;
//DEL 	double f;
//DEL 
//DEL 	if (xdim)
//DEL 	{
//DEL 		for (x=0;x<m_iRes[0];x++)
//DEL 		{
//DEL 			f = (double)(1.0f/(pow(m_fMinVal[0]+(x+1.0f)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),3)-pow(m_fMinVal[0]+((double)x)/m_iRes[0]*(m_fMaxVal[0]-m_fMinVal[0]),3)));
//DEL 			for (y=0;y<m_iRes[1];y++)
//DEL 				m_pBin[x+y*(m_iRes[0]+1)] *= f;
//DEL 		}
//DEL 	} else
//DEL 	{
//DEL 		for (y=0;y<m_iRes[1];y++)
//DEL 		{
//DEL 			f = (double)(1.0f/(pow(m_fMinVal[1]+(y+1.0f)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),3)-pow(m_fMinVal[1]+((double)y)/m_iRes[1]*(m_fMaxVal[1]-m_fMinVal[1]),3)));
//DEL 			for (x=0;x<m_iRes[0];y++)
//DEL 				m_pBin[x+y*(m_iRes[0]+1)] *= f;
//DEL 		}
//DEL 	}
//DEL 	BTOUT;
//DEL }


void C2DF::MultiplyBin(double m)
{
	BTIN;
	int z;

	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
		m_pBin[z] *= m;
	BTOUT;
}


double C2DF::PPMBin()
{
	BTIN;
	double f;
	int z;

	f = 0;
	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
	{
		m_pBin[z] *= 1000000.0f / (double)g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize();
		if (f < m_pBin[z])
			f = m_pBin[z];
	}
	BTOUT;
	return f;
}


void C2DF::NormalizeBin(double mi, double ma)
{
	BTIN;
	int z;
	double tmi, tma, d, td;

	tmi = 99999999.0f;
	tma = 0.0f;
	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
	{
		if (m_pBin[z] < tmi)
			tmi = m_pBin[z];
		if (m_pBin[z] > tma)
			tma = m_pBin[z];
	}
	if (tma-tmi < 1E-20f)
		tma += 0.00001f;
	d = ma - mi;
	td = tma - tmi;
	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
		m_pBin[z] = ((m_pBin[z]-tmi)/td*d)+mi;
	BTOUT;
}


double C2DF::NormalizeBinIntegral(double val)
{
	BTIN;
	int z;
	double d;

	d = 0;
	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
		d += m_pBin[z];
	if (d == 0)
	{
		BTOUT;
		return 0.0;
	}
	d = val/d;
	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
		m_pBin[z]*=d;
	return d;
	BTOUT;
}


void C2DF::Create()
{
	BTIN;
	int z;
	m_fBinEntries = 0;
	m_fSkipEntries = 0;

	try { m_pStepsY = new unsigned long[m_iRes[1]]; } catch(...) { m_pStepsY = NULL; }
	if (m_pStepsY == NULL) NewException((double)m_iRes[1]*sizeof(unsigned long),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_fCountX = new double[m_iRes[0]]; } catch(...) { m_fCountX = NULL; }
	if (m_fCountX == NULL) NewException((double)m_iRes[0]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_fCountY = new double[m_iRes[1]]; } catch(...) { m_fCountY = NULL; }
	if (m_fCountY == NULL) NewException((double)m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_iRes[0];z++)
		m_fCountX[z] = 0;

	for (z=0;z<m_iRes[1];z++)
	{
		m_pStepsY[z] = 0;
		m_fCountY[z] = 0;
	}

	try { m_pBin = new double[m_iRes[0]*m_iRes[1]]; } catch(...) { m_fCountY = NULL; }
	if (m_fCountY == NULL) NewException((double)m_iRes[0]*m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
		m_pBin[z] = 0;
	m_fFac[0] = (double)m_iRes[0] / (m_fMaxVal[0]-m_fMinVal[0]);
	m_fFac[1] = (double)m_iRes[1] / (m_fMaxVal[1]-m_fMinVal[1]);
	BTOUT;
}


void C2DF::Write(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z, z2;
//	char buf[32768];
	CxString buf;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);

	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	if (m_sLabelX == NULL)
		mfprintf(a,"# (no label); ");
			else mfprintf(a,"# %s; ",m_sLabelX);
	if (m_sLabelY == NULL)
		mfprintf(a,"(no label); ");
			else mfprintf(a,"%s; ",m_sLabelY);
	if (m_sLabelZ == NULL)
		mfprintf(a,"Occurrence\n");
			else mfprintf(a,"%s\n",m_sLabelZ);
	for (z=0;z<m_iRes[1];z++)
		for (z2=0;z2<m_iRes[0];z2++)
			mfprintf(a,"%#.10G;  %#.10G;  %#.10G\n",m_fMinVal[0]+(z2+0.5)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0],m_fMinVal[1]+(z+0.5)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1],m_pBin[z*m_iRes[0]+z2]);
	fclose(a);
	
	BTOUT;
}


void C2DF::WriteGraceBunch(int channel, int graphs, float fac, const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	int x, y, t;
//	char buf[32768];
	CxString buf;
	CGrace *g;

	try { g = new CGrace(); } catch(...) { g = NULL; }
	if (g == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	g->CurrentGraph()->m_bLegend = true;

	if (channel == 0)
	{
		g->SetLabelX(m_sLabelY);
		if (m_sLabelZ == NULL)
			g->SetLabelY("Occurrence");
				else g->SetLabelY(m_sLabelZ);
		CalcMaxEntry();
		g->SetRangeX(m_fMinVal[1],m_fMaxVal[1]);
		g->SetRangeY(0,m_fMaxEntry*1.1*fac);
		g->MakeTicks();
//		sprintf(buf,"Legend unit is %s",m_sLabelX);
		buf.sprintf("Legend unit is %s",m_sLabelX);
		g->SetSubTitle(buf);

		for (x=0;x<graphs;x++)
		{
			t = (int)((double)x / (graphs-1) * m_iRes[0]);
			g->AddDataset();
//			sprintf(buf,"%.2f",m_fMinVal[0]+t*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0]);
			buf.sprintf("%.2f",m_fMinVal[0]+t*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0]);
			g->SetDatasetName(buf);
			g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iRes[1]*2);
			for (y=0;y<m_iRes[1]-1;y++)
				g->AddXYTupel(m_fMinVal[1]+(y+0.5)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1],m_pBin[y*m_iRes[0]+t]*fac);
		}
	} else
	{
		g->SetLabelX(m_sLabelX);
		if (m_sLabelZ == NULL)
			g->SetLabelY("Occurrence");
				else g->SetLabelY(m_sLabelZ);
		CalcMaxEntry();
		g->SetRangeX(m_fMinVal[0],m_fMaxVal[0]);
		g->SetRangeY(0,m_fMaxEntry*1.1*fac);
		g->MakeTicks();
//		sprintf(buf,"Legend unit is %s",m_sLabelY);
		buf.sprintf("Legend unit is %s",m_sLabelY);
		g->SetSubTitle(buf);

		for (y=0;y<graphs;y++)
		{
			t = (int)((double)y / (graphs-1) * m_iRes[1]);
			g->AddDataset();
//			sprintf(buf,"%.2f",m_fMinVal[1]+t*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1]);
			buf.sprintf("%.2f",m_fMinVal[1]+t*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1]);
			g->SetDatasetName(buf);
			g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iRes[0]*2);
			for (x=0;x<m_iRes[0]-1;x++)
				g->AddXYTupel(m_fMinVal[0]+(x+0.5)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0],m_pBin[t*m_iRes[0]+x]*fac);
		}
	}
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);
	g->WriteAgr(buf,false);
	*strrchr(buf.GetWritePointer(),'.') = 0;
	buf.strcat(".csv");
	g->WriteCSV(buf);

	BTOUT;
}


void C2DF::WriteCSV(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z, z2;
//	char buf[32768];
	CxString buf;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	if (m_sLabelX == NULL)
		mfprintf(a,"(no label) \\ ");
			else mfprintf(a,"%s \\ ",m_sLabelX);
	if (m_sLabelY == NULL)
		mfprintf(a,"(no label)");
			else mfprintf(a,"%s",m_sLabelY);
	for (z=0;z<m_iRes[1];z++)
			mfprintf(a,"; %#.10G",m_fMinVal[1]+(z+0.5)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1]);
	mfprintf(a,"\n");
	for (z2=0;z2<m_iRes[0];z2++)
	{
		mfprintf(a,"%#.10G",m_fMinVal[0]+(z2+0.5)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0]);
		for (z=0;z<m_iRes[1];z++)
			mfprintf(a,"; %#.10G",m_pBin[z*m_iRes[0]+z2]);
		mfprintf(a,"\n");
	}
	fclose(a);
	BTOUT;
}


void C2DF::WriteXProjection(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z, z2;
//	char buf[32768];
	CxString buf;
	double d;
	double *tb;

	try { tb = new double[m_iRes[0]]; } catch(...) { tb = NULL; }
	if (tb == NULL) NewException((double)m_iRes[0]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_iRes[0];z++)
	{
		d = 0;
		for (z2=0;z2<m_iRes[1];z2++)
			d += m_pBin[z2*m_iRes[0]+z];
		d /= m_iRes[1];
		tb[z] = d;
	}
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	for (z=0;z<m_iRes[0];z++)
		mfprintf(a,"%f;  %f\n",m_fMinVal[0]+(z+0.5)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0],tb[z]);
	fclose(a);
	delete[] tb;
	BTOUT;
}


void C2DF::WriteYProjection(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z, z2;
//	char buf[32768];
	CxString buf;
	double d;
	double *tb;

	try { tb = new double[m_iRes[1]]; } catch(...) { tb = NULL; }
	if (tb == NULL) NewException((double)m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_iRes[1];z++)
	{
		d = 0;
		for (z2=0;z2<m_iRes[0];z2++)
			d += m_pBin[z*m_iRes[0]+z2];
		d /= m_iRes[1];
		tb[z] = d;
	}
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	for (z=0;z<m_iRes[1];z++)
		mfprintf(a,"%f;  %f\n",m_fMinVal[1]+(z+0.5)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1],tb[z]);
	fclose(a);
	delete[] tb;
	BTOUT;
}


void C2DF::WriteMathematica(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z, z2, i;
//	char buf[32768];
	CxString buf;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	mfprintf(a,"{ ");
	i = 0;
	for (z=0;z<m_iRes[1];z++)
	{
		for (z2=0;z2<m_iRes[0];z2++)
		{
			mfprintf(a,"{ %#.10f,  %#.10f,  %#.10f } ",m_fMinVal[0]+(z2+0.5)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0],m_fMinVal[1]+(z+0.5)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1],m_pBin[z*m_iRes[0]+z2]);
			i++;
			if ((z2+1 < m_iRes[0]) || (z+1 < m_iRes[1]))
				mfprintf(a,", ");
			if (i % 5 == 0)
				mfprintf(a,"\n");
		}
	}
	mfprintf(a," }\n");
	fclose(a);
	BTOUT;
}


void C2DF::WriteGnuplotInput(const char *prefix, const char *name, const char *suffix, bool manrange)
{
	FILE *a, *b;
	int z, z2;
//	char buf1[32768], buf2[32768], buf[32768];
	CxString buf1, buf2, buf, out, out2;
	int minorx, majorx, minory, majory;

//	sprintf(buf1,"%s%s%s.gp.csv",prefix,name,suffix);
	buf1.sprintf("%s%s%s.gp.csv",prefix,name,suffix);
	a = OpenFileWrite(buf1,true,&out);

	out2 = out;
	out2(out2.GetLength()-7) = 0;

	for (z2=0;z2<m_iRes[0];z2++)
	{
		for (z=0;z<m_iRes[1];z++)
		{
			mfprintf(a,"%.10G",m_pBin[z*m_iRes[0]+z2]);
			if (z < m_iRes[1]-1)
				mfprintf(a,"; ");
		}
		mfprintf(a,"\n");
	}

	fclose(a);

//	sprintf(buf2,"%s%s%s.gp",prefix,name,suffix);
	buf2.sprintf("%s%s%s.gp",prefix,name,suffix);
	b = OpenFileWrite(buf2,true);

	if (!manrange)
		CalcMaxEntry();

	CreateTicks(m_fMinVal[0],m_fMaxVal[0],majorx,minorx);
	CreateTicks(m_fMinVal[1],m_fMaxVal[1],majory,minory);

	mfprintf(b,"# Total plot width in pixel\n");
	mfprintf(b,"s_widthpixel = %d\n",m_iPlotPixel);
	mfprintf(b,"\n");
	mfprintf(b,"# Border width scaling\n");
	mfprintf(b,"s_borderscale = 1.0\n");
	mfprintf(b,"# Grid width scaling\n");
	mfprintf(b,"s_gridscale = 1.0\n");
	mfprintf(b,"# Contour width scaling\n");
	mfprintf(b,"s_contourscale = 1.0\n");
	mfprintf(b,"\n");
	mfprintf(b,"# Title font scaling\n");
	mfprintf(b,"s_titlefontscale = 1.0\n");
	mfprintf(b,"# Axis title font scaling\n");
	mfprintf(b,"s_labelfontscale = 1.0\n");
	mfprintf(b,"# Tics font scaling\n");
	mfprintf(b,"s_ticsfontscale = 1.0\n");
	mfprintf(b,"\n");
	mfprintf(b,"# Aspect ratio\n");
	mfprintf(b,"s_ratio = %f\n",m_fAspectRatio);
	mfprintf(b,"\n");
	mfprintf(b,"# Global title\n");
	mfprintf(b,"s_title = \"\"\n");

	mfprintf(b,"# X axis title\n");
	ProtectCharacters(&buf,m_sLabelX,"_^","\\\\");
	mfprintf(b,"s_xtitle = \"%s\"\n",(const char*)buf);

	mfprintf(b,"# Y axis title\n");
	ProtectCharacters(&buf,m_sLabelY,"_^","\\\\");
	mfprintf(b,"s_ytitle = \"%s\"\n",(const char*)buf);

	mfprintf(b,"\n");
	mfprintf(b,"# X interpolation steps\n");
	mfprintf(b,"s_xipl = %d\n",m_iGPInterpolation);
	mfprintf(b,"# Y interpolation steps\n");
	mfprintf(b,"s_yipl = %d\n",m_iGPInterpolation);
	mfprintf(b,"\n");
	mfprintf(b,"# Draw contour lines (0 = No, 1 = Yes)\n");
	mfprintf(b,"s_contour = %d\n",m_bContourLines?1:0);
	mfprintf(b,"# Number of contour levels\n");
	mfprintf(b,"s_ncontour = %d\n",g_pDatabase->GetInt("/PLOT2D/DEFAULTS/CONTOUR_LINES"));
	mfprintf(b,"\n");
	mfprintf(b,"# Plotting function\n");
	mfprintf(b,"f_plot(x) = x**%G\n",m_fPlotExp);
	mfprintf(b,"\n");
	mfprintf(b,"# Minimum X\n");
	mfprintf(b,"s_xmin = %f\n",m_fMinVal[0]);
	mfprintf(b,"# Maximum X\n");
	mfprintf(b,"s_xmax = %f\n",m_fMaxVal[0]);
	mfprintf(b,"# Major X tics distance\n");
	mfprintf(b,"s_xtics = %f\n",(m_fMaxVal[0]-m_fMinVal[0])/(majorx-1));
	mfprintf(b,"# Minor X tics number\n");
	mfprintf(b,"s_mxtics = %d\n",minorx+1);
	mfprintf(b,"\n");
	mfprintf(b,"# Minimum Y\n");
	mfprintf(b,"s_ymin = %f\n",m_fMinVal[1]);
	mfprintf(b,"# Maximum Y\n");
	mfprintf(b,"s_ymax = %f\n",m_fMaxVal[1]);
	mfprintf(b,"# Major Y tics distance\n");
	mfprintf(b,"s_ytics = %f\n",(m_fMaxVal[1]-m_fMinVal[1])/(majory-1));
	mfprintf(b,"# Minor Y tics number\n");
	mfprintf(b,"s_mytics = %d\n",minory+1);
	mfprintf(b,"\n");
	mfprintf(b,"# Use automatic Z values (0 = No, 1 = Yes)\n");
	mfprintf(b,"s_zauto = 1\n");
	mfprintf(b,"# Minimum Z\n");
	mfprintf(b,"s_zmin = 0\n");
	mfprintf(b,"# Maximum Z\n");
	mfprintf(b,"s_zmax = 1000\n");
	mfprintf(b,"# Major Z tics distance\n");
	mfprintf(b,"s_ztics = 500\n");
	mfprintf(b,"# Minor Z tics number\n");
	mfprintf(b,"s_mztics = 1\n");
	mfprintf(b,"\n");
	mfprintf(b,"# Legend orientation (0 = Horizontal, 1 = Vertical)\n");
	mfprintf(b,"s_lorient = 0\n");
	mfprintf(b,"\n");
	mfprintf(b,"# The red line\n");
	mfprintf(b,"#set arrow 1 from 50,50 to 200,200 nohead back linewidth 0.7 linecolor rgbcolor \"red\"\n");
	mfprintf(b,"\n");
	mfprintf(b,"# Do not change anything below\n");
	mfprintf(b,"\n");
	mfprintf(b,"s_scale = s_widthpixel / 1024.0\n");
	mfprintf(b,"\n");
	mfprintf(b,"set term unknown\n");
	mfprintf(b,"\n");
	mfprintf(b,"splot \"%s\" matrix using ($2*%G+%G):($1*%G+%G):3\n",(const char*)out,(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0],m_fMinVal[0]+(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0]/2.0,(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1],m_fMinVal[1]+(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1]/2.0);
	mfprintf(b,"\n");
	mfprintf(b,"if (s_zauto) s_zmin = GPVAL_Z_MIN; s_zmax = GPVAL_Z_MAX\n");
	mfprintf(b,"\n");
	mfprintf(b,"font_string = sprintf(\"%cs,%cd\", \"Helvetica\", 24*s_scale)\n",'%','%');
	mfprintf(b,"titlefont_string = sprintf(\"%cs,%cd\", \"Helvetica\", 24*1.2*s_scale*s_titlefontscale)\n",'%','%');
	mfprintf(b,"labelfont_string = sprintf(\"%cs,%cd\", \"Helvetica\", 24*s_scale*s_labelfontscale)\n",'%','%');
	mfprintf(b,"ticsfont_string = sprintf(\"%cs,%cd\", \"Helvetica\", 24*s_scale*s_ticsfontscale)\n",'%','%');
	mfprintf(b,"\n");
	mfprintf(b,"set term pngcairo enhanced font font_string linewidth 4.0*s_scale size 1024*s_scale,1024*s_scale\n");
	mfprintf(b,"set output \"%s.png\"\n",(const char*)out2);
	mfprintf(b,"set size ratio s_ratio\n");
	mfprintf(b,"\n");
	mfprintf(b,"set title s_title font titlefont_string\n");
	mfprintf(b,"set xlabel s_xtitle font labelfont_string\n");
	mfprintf(b,"set ylabel s_ytitle font labelfont_string\n");
	mfprintf(b,"unset key\n");
	mfprintf(b,"unset clabel\n");
	mfprintf(b,"\n");
	mfprintf(b,"set pm3d interpolate s_xipl,s_yipl map\n");
	mfprintf(b,"\n");
	mfprintf(b,"set palette model RGB functions gray<=0.0 ? 1.0 : gray<0.2 ? 1.0-(5.0**1.5*gray**1.5)*0.8 : gray<0.4 ? 0.2-(gray-0.2) : gray<0.6 ? 5.0**0.5*(gray-0.4)**0.5 : gray<0.8 ? 1.0 : gray<1.0 ? 1.0 : 1.0, gray<=0.0 ? 1.0 : gray<0.2 ? 1.0-(5.0**1.5*gray**1.5)*0.8 : gray<0.4 ? 0.2+0.8*(5.0**0.75*(gray-0.2)**0.75) : gray<0.6 ? 1.0 : gray<0.8 ? 1.0-5.0*(gray-0.6) : gray<1.0 ? 0.0 : 0.0, gray<=0.0 ? 1.0 : gray<0.2 ? 1.0 : gray<0.4 ? 1.0-5.0**1.33*(gray-0.2)**1.33 : gray<0.6 ? 0.0 : gray<0.8 ? 0.0 : gray<1.0 ? 5.0*(gray-0.8) : 1.0\n");
	mfprintf(b,"\n");
	mfprintf(b,"if (s_contour) set palette maxcolors s_ncontour; else set palette maxcolors 0\n");
	mfprintf(b,"\n");
	mfprintf(b,"set xrange [s_xmin:s_xmax]\n");
	mfprintf(b,"set yrange [s_ymin:s_ymax]\n");
	mfprintf(b,"set zrange [f_plot(s_zmin):f_plot(s_zmax)]\n");
	mfprintf(b,"set cbrange [f_plot(s_zmin):f_plot(s_zmax)]\n");
	mfprintf(b,"unset colorbox\n");
	mfprintf(b,"\n");
	mfprintf(b,"set tics font ticsfont_string\n");
	mfprintf(b,"set xtics s_xtics out scale 0.5*s_borderscale offset 0,0.5\n");
	mfprintf(b,"set mxtics s_mxtics\n");
	mfprintf(b,"set ytics s_ytics out scale 0.5*s_borderscale offset 0.5,0\n");
	mfprintf(b,"set mytics s_mytics\n");
	mfprintf(b,"\n");
	mfprintf(b,"set border linewidth s_borderscale\n");
	mfprintf(b,"\n");
	mfprintf(b,"set multiplot\n");
	mfprintf(b,"\n");
	mfprintf(b,"set origin -0.05,-0.05\n");
	mfprintf(b,"set size 1.15,1.15\n");
	mfprintf(b,"\n");
	mfprintf(b,"splot \"%s\" matrix using ($2*%G+%G):($1*%G+%G):(f_plot($3))\n",(const char*)out,(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0],m_fMinVal[0]+(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0]/2.0,(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1],m_fMinVal[1]+(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1]/2.0);
//	mfprintf(b,"splot \"test5.dat\" matrix using ($2*1.5+50.75):($1*5.0+2.5):(f_plot($3))\n");
	mfprintf(b,"\n");
	mfprintf(b,"set grid xtics mxtics ytics mytics back linestyle -1 linewidth 0.15*s_gridscale linecolor rgbcolor \"gray40\"\n");
	mfprintf(b,"set title \"\"\n");
	mfprintf(b,"set xlabel \"\"\n");
	mfprintf(b,"set ylabel \"\"\n");
	mfprintf(b,"set xtics in format \"\"\n");
	mfprintf(b,"set ytics in format \"\"\n");
	mfprintf(b,"if (s_contour) set contour base; set cntrparam bspline; set cntrparam levels incremental f_plot(s_zmin),1.0*(f_plot(s_zmax)-f_plot(s_zmin))/s_ncontour,f_plot(s_zmax)\n");
	mfprintf(b,"\n");
	mfprintf(b,"splot \"%s\" matrix using ($2*%G+%G):($1*%G+%G):(s_contour ? f_plot($3) : 1/0) with lines linewidth 0.15*s_contourscale linecolor rgbcolor \"black\" nosurface\n",(const char*)out,(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0],m_fMinVal[0]+(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0]/2.0,(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1],m_fMinVal[1]+(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1]/2.0);
	mfprintf(b,"\n");
	mfprintf(b,"unset multiplot\n");
	mfprintf(b,"\n");
	mfprintf(b,"if (s_lorient) t_x = 196*s_scale; else t_x = 896*s_scale\n");
	mfprintf(b,"if (s_lorient) t_y = 896*s_scale; else t_y = 196*s_scale\n");
	mfprintf(b,"\n");
	mfprintf(b,"set term pngcairo enhanced font font_string linewidth 4.0*s_scale size t_x,t_y\n");
	mfprintf(b,"set output \"%s_box.png\"\n",(const char*)out2);
	mfprintf(b,"\n");
	mfprintf(b,"if (s_lorient) t_ratio = 20; else t_ratio = 0.05\n");
	mfprintf(b,"\n");
	mfprintf(b,"set origin -0.13,-0.13\n");
	mfprintf(b,"set size ratio t_ratio 1.25,1.25\n");
	mfprintf(b,"\n");
	mfprintf(b,"unset title\n");
	mfprintf(b,"unset xlabel\n");
	mfprintf(b,"unset ylabel\n");
	mfprintf(b,"\n");
	mfprintf(b,"set grid front\n");
	mfprintf(b,"unset grid\n");
	mfprintf(b,"\n");
	mfprintf(b,"if (s_lorient) set xrange [0:1]; else set xrange [s_zmin:s_zmax]\n");
	mfprintf(b,"if (s_lorient) set yrange [s_zmin:s_zmax]; else set yrange [0:1]\n");
	mfprintf(b,"if (s_lorient) set format y; else set format x\n");
	mfprintf(b,"if (s_zauto) set xtics scale 0.5*s_borderscale offset s_lorient ? -0.5 : 0,s_lorient ? 0 : 0.5 autofreq; else set xtics scale 0.5*s_borderscale offset 0,0.5 s_ztics\n");
	mfprintf(b,"if (s_zauto) set mxtics default; else set mxtics s_mztics\n");
	mfprintf(b,"if (s_zauto) set ytics scale 0.5*s_borderscale offset s_lorient ? -0.5 : 0,s_lorient ? 0 : 0.5 autofreq; else set ytics scale 0.5*s_borderscale offset 0,0.5 s_ztics\n");
	mfprintf(b,"if (s_zauto) set mytics default; else set mytics s_mztics\n");
	mfprintf(b,"if (s_lorient) unset xtics; else unset ytics\n");
	mfprintf(b,"\n");
	mfprintf(b,"if (s_lorient) set isosamples 2,500\n");
	mfprintf(b,"\n");
	mfprintf(b,"if (s_lorient) set view 0,0; set border 15; unset ztics; set size 0.85,1.5; set origin -0.13,-0.25\n");
	mfprintf(b,"\n");
	mfprintf(b,"splot f_plot(s_lorient ? y : x) nocontour, s_contour ? f_plot(s_lorient ? y : x) : 1/0 with lines linewidth 0.15*s_contourscale linecolor rgbcolor \"black\" nosurface\n");

	fclose(b);
}


void C2DF::WriteMathematicaNb(const char *prefix, const char *name, const char *suffix, bool manrange)
{
	BTIN;
	FILE *a;
	int z, z2;
//	char buf[32768];
	CxString buf;
	int minorx, majorx, minory, majory;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	if (!manrange)
		CalcMaxEntry();

	if ((m_iPlotType == 1) && (!m_bContourLines))
		m_iPlotType = 2;

	CreateTicks(m_fMinVal[0],m_fMaxVal[0],majorx,minorx);
	CreateTicks(m_fMinVal[1],m_fMaxVal[1],majory,minory);

	a = OpenFileWrite(buf,true);

	mfprintf(a,"\n");
	mfprintf(a,"Notebook[{\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," StyleBox[\n");
	mfprintf(a,"  RowBox[{\"TRAVIS\", \" \", \"Analysis\"}], \"Section\",\n");
	mfprintf(a,"  Evaluatable->False]], \"Input\",\n");
	mfprintf(a," Evaluatable->False],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"Cell[\"Input Data\", \"Section\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"dat\", \"=\",\n");
	mfprintf(a,"   RowBox[{\"dat2\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"{\", \n");
	mfprintf(a,"     RowBox[{\n");

	for (z=0;z<m_iRes[1];z++)
	{
		for (z2=0;z2<m_iRes[0];z2++)
		{
			mfprintf(a,"     RowBox[{\"{\", \n");
			mfprintf(a,"      RowBox[{\"%#.10f\", \",\", \"%#.10f\", \",\", \"%#.10f\"}], \"}\"}]",m_fMinVal[0]+(z2+0.5)*(m_fMaxVal[0]-m_fMinVal[0])/m_iRes[0],m_fMinVal[1]+(z+0.5)*(m_fMaxVal[1]-m_fMinVal[1])/m_iRes[1],m_pBin[z*m_iRes[0]+z2]);
			if ((z < m_iRes[1]-1) || (z2 < m_iRes[0]-1))
				mfprintf(a,", \",\", \n");
		}
	}
	
	mfprintf(a,"}], \"}\"}]}]}], \n");
	mfprintf(a,"  \";\"}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	
	mfprintf(a,"Cell[\"Internal Routines\", \"Section\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[{\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"prim\", \"=\", \n");
	mfprintf(a,"   RowBox[{\"{\", \"}\"}]}], \";\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"funcX\", \"=\", \n");
	mfprintf(a,"   RowBox[{\"{\", \"}\"}]}], \";\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"funcY\", \"=\", \n");
	mfprintf(a,"   RowBox[{\"{\", \"}\"}]}], \";\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"x\", \"=.\"}], \";\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"s\", \"=.\"}], \";\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\n");

	mfprintf(a,"   RowBox[{\"maxsig\", \"[\",\n");
	mfprintf(a,"    RowBox[{\"x_\", \",\", \"s_\"}], \"]\"}], \"=\",\n");
	mfprintf(a,"   RowBox[{\"N\", \"[\",\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"      RowBox[{\n");
	mfprintf(a,"       RowBox[{\"x\", \"\\[Equal]\", \"0\"}], \",\", \"0\", \",\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\"Ceiling\", \"[\", \n");
	mfprintf(a,"         RowBox[{\"x\", \"/\", \n");
	mfprintf(a,"          RowBox[{\"10\", \"^\", \n");
	mfprintf(a,"           RowBox[{\"Floor\", \"[\", \n");
	mfprintf(a,"            RowBox[{\"N\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"Log\", \"[\", \n");
	mfprintf(a,"               RowBox[{\"10\", \",\", \n");
	mfprintf(a,"                RowBox[{\"Abs\", \"[\", \"x\", \"]\"}]}], \"]\"}], \"+\", \"1\", \"-\", \"s\"}],\n");
	mfprintf(a,"              \"]\"}], \"]\"}]}]}], \"]\"}], \"*\", \n");
	mfprintf(a,"        RowBox[{\"10\", \"^\", \n");
	mfprintf(a,"         RowBox[{\"Floor\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"N\", \"[\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Log\", \"[\", \n");
	mfprintf(a,"             RowBox[{\"10\", \",\", \n");
	mfprintf(a,"              RowBox[{\"Abs\", \"[\", \"x\", \"]\"}]}], \"]\"}], \"+\", \"1\", \"-\", \"s\"}], \n");
	mfprintf(a,"           \"]\"}], \"]\"}]}]}]}], \"]\"}], \"+\", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"Abs\", \"[\", \"x\", \"]\"}], \"/\", \"1000000\"}]}], \"]\"}]}], \n");
	mfprintf(a,"  \";\"}], \"\\n\", \n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minsig\", \"[\", \n");
	mfprintf(a,"    RowBox[{\"x_\", \",\", \"s_\"}], \"]\"}], \"=\", \n");
	mfprintf(a,"   RowBox[{\"N\", \"[\", \n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"      RowBox[{\n");
	mfprintf(a,"       RowBox[{\"x\", \"\\[Equal]\", \"0\"}], \",\", \"0\", \",\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\"Floor\", \"[\", \n");
	mfprintf(a,"         RowBox[{\"x\", \"/\", \n");
	mfprintf(a,"          RowBox[{\"10\", \"^\", \n");
	mfprintf(a,"           RowBox[{\"Floor\", \"[\", \n");
	mfprintf(a,"            RowBox[{\"N\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"Log\", \"[\", \n");
	mfprintf(a,"               RowBox[{\"10\", \",\", \n");
	mfprintf(a,"                RowBox[{\"Abs\", \"[\", \"x\", \"]\"}]}], \"]\"}], \"+\", \"1\", \"-\", \"s\"}],\n");
	mfprintf(a,"              \"]\"}], \"]\"}]}]}], \"]\"}], \"*\", \n");
	mfprintf(a,"        RowBox[{\"10\", \"^\", \n");
	mfprintf(a,"         RowBox[{\"Floor\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"N\", \"[\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Log\", \"[\", \n");
	mfprintf(a,"             RowBox[{\"10\", \",\", \n");
	mfprintf(a,"              RowBox[{\"Abs\", \"[\", \"x\", \"]\"}]}], \"]\"}], \"+\", \"1\", \"-\", \"s\"}], \n");
	mfprintf(a,"           \"]\"}], \"]\"}]}]}]}], \"]\"}], \"-\", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"Abs\", \"[\", \"x\", \"]\"}], \"/\", \"1000000\"}]}], \"]\"}]}], \n");
	mfprintf(a,"  \";\"}], \"\\n\", \n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"prec\", \"[\", \"x_\", \"]\"}], \"=\", \n");
	mfprintf(a,"   RowBox[{\"Max\", \"[\", \n");
	mfprintf(a,"    RowBox[{\"0\", \",\", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"Floor\", \"[\", \n");
	mfprintf(a,"       RowBox[{\"N\", \"[\", \n");
	mfprintf(a,"        RowBox[{\"Log\", \"[\", \n");
	mfprintf(a,"         RowBox[{\"10\", \",\", \n");
	mfprintf(a,"          RowBox[{\"1\", \"/\", \n");
	mfprintf(a,"           RowBox[{\"Abs\", \"[\", \"x\", \"]\"}]}]}], \"]\"}], \"]\"}], \"]\"}], \"+\", \n");
	mfprintf(a,"      \"2\"}]}], \"]\"}]}], \";\"}], \"\\[IndentingNewLine]\", \n");

	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"plot\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"x\", \"=.\"}], \";\", \n");
	mfprintf(a,"      RowBox[{\"y\", \"=.\"}], \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"(*\", \" \", \n");

	mfprintf(a,"       RowBox[{\"Select\", \" \", \"Color\", \" \", \"Function\"}], \" \", \"*)\"}],");
	mfprintf(a,"      \"\\[IndentingNewLine]\",");
	mfprintf(a,"      RowBox[{\"If\", \"[\", ");
	mfprintf(a,"       RowBox[{");
	mfprintf(a,"        RowBox[{\"coloring\", \"\\[Equal]\", \"1\"}], \",\", ");
	mfprintf(a,"        RowBox[{");
	mfprintf(a,"         RowBox[{\"ColFunc\", \"[\", \"x_\", \"]\"}], \"=\", ");
	mfprintf(a,"         RowBox[{\"If\", \"[\", ");
	mfprintf(a,"          RowBox[{");
	mfprintf(a,"           RowBox[{\"x\", \"\\[LessEqual]\", \"0\"}], \",\", ");
	mfprintf(a,"           RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"            RowBox[{\"1\", \",\", \"1\", \",\", \"1\"}], \"]\"}], \",\", ");
	mfprintf(a,"           RowBox[{\"If\", \"[\", ");
	mfprintf(a,"            RowBox[{");
	mfprintf(a,"             RowBox[{\"x\", \"<\", ");
	mfprintf(a,"              RowBox[{\"1\", \"/\", \"5\"}]}], \",\", ");
	mfprintf(a,"             RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"              RowBox[{");
	mfprintf(a,"               RowBox[{\"1\", \"-\", ");
	mfprintf(a,"                RowBox[{");
	mfprintf(a,"                 RowBox[{\"(\", ");
	mfprintf(a,"                  RowBox[{");
	mfprintf(a,"                   RowBox[{\"5\", \"^\", \"1.5\"}], \"*\", ");
	mfprintf(a,"                   RowBox[{\"x\", \"^\", \"1.5\"}]}], \")\"}], \"*\", \"0.8\"}]}], \",\", ");
	mfprintf(a,"               RowBox[{\"1\", \"-\", ");
	mfprintf(a,"                RowBox[{");
	mfprintf(a,"                 RowBox[{\"(\", ");
	mfprintf(a,"                  RowBox[{");
	mfprintf(a,"                   RowBox[{\"5\", \"^\", \"1.5\"}], \"*\", ");
	mfprintf(a,"                   RowBox[{\"x\", \"^\", \"1.5\"}]}], \")\"}], \"*\", \"0.8\"}]}], \",\", ");
	mfprintf(a,"               \"1\"}], \"]\"}], \",\", ");
	mfprintf(a,"             RowBox[{\"If\", \"[\", ");
	mfprintf(a,"              RowBox[{");
	mfprintf(a,"               RowBox[{\"x\", \"<\", ");
	mfprintf(a,"                RowBox[{\"2\", \"/\", \"5\"}]}], \",\", ");
	mfprintf(a,"               RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                RowBox[{");
	mfprintf(a,"                 RowBox[{\"0.2\", \"-\", ");
	mfprintf(a,"                  RowBox[{\"(\", ");
	mfprintf(a,"                   RowBox[{\"x\", \"-\", ");
	mfprintf(a,"                    RowBox[{\"1\", \"/\", \"5\"}]}], \")\"}]}], \",\", ");
	mfprintf(a,"                 RowBox[{\"0.2\", \"+\", ");
	mfprintf(a,"                  RowBox[{\"0.8\", \"*\", ");
	mfprintf(a,"                   RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"5\", \"^\", \"0.75\"}], \"*\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"x\", \"-\", ");
	mfprintf(a,"                    RowBox[{\"1\", \"/\", \"5\"}]}], \")\"}], \"^\", \"0.75\"}]}], ");
	mfprintf(a,"                    \")\"}]}]}], \",\", ");
	mfprintf(a,"                 RowBox[{\"1\", \"-\", ");
	mfprintf(a,"                  RowBox[{");
	mfprintf(a,"                   RowBox[{\"5\", \"^\", \"1.33\"}], \"*\", ");
	mfprintf(a,"                   RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"x\", \"-\", ");
	mfprintf(a,"                    RowBox[{\"1\", \"/\", \"5\"}]}], \")\"}], \"^\", \"1.33\"}]}]}]}], ");
	mfprintf(a,"                \"]\"}], \",\", ");
	mfprintf(a,"               RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                RowBox[{");
	mfprintf(a,"                 RowBox[{\"x\", \"<\", ");
	mfprintf(a,"                  RowBox[{\"3\", \"/\", \"5\"}]}], \",\", ");
	mfprintf(a,"                 RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                  RowBox[{");
	mfprintf(a,"                   RowBox[{");
	mfprintf(a,"                    RowBox[{\"5\", \"^\", \"0.5\"}], \"*\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"x\", \"-\", ");
	mfprintf(a,"                    RowBox[{\"2\", \"/\", \"5\"}]}], \")\"}], \"^\", \"0.5\"}]}], \",\", ");
	mfprintf(a,"                   \"1\", \",\", \"0\"}], \"]\"}], \",\", ");
	mfprintf(a,"                 RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                  RowBox[{");
	mfprintf(a,"                   RowBox[{\"x\", \"<\", ");
	mfprintf(a,"                    RowBox[{\"4\", \"/\", \"5\"}]}], \",\", ");
	mfprintf(a,"                   RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{\"1\", \",\", ");
	mfprintf(a,"                    RowBox[{\"1\", \"-\", ");
	mfprintf(a,"                    RowBox[{\"5\", \"*\", ");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"x\", \"-\", ");
	mfprintf(a,"                    RowBox[{\"3\", \"/\", \"5\"}]}], \")\"}]}]}], \",\", \"0\"}], \"]\"}], ");
	mfprintf(a,"                   \",\", ");
	mfprintf(a,"                   RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"x\", \"<\", \"1\"}], \",\", ");
	mfprintf(a,"                    RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{\"1\", \",\", \"0\", \",\", ");
	mfprintf(a,"                    RowBox[{\"5\", \"*\", ");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"x\", \"-\", ");
	mfprintf(a,"                    RowBox[{\"4\", \"/\", \"5\"}]}], \")\"}]}]}], \"]\"}], \",\", ");
	mfprintf(a,"                    RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{\"1\", \",\", \"0\", \",\", \"1\"}], \"]\"}]}], \"]\"}]}], ");
	mfprintf(a,"                  \"]\"}]}], \"]\"}]}], \"]\"}]}], \"]\"}]}], \"]\"}]}], \",\", ");
	mfprintf(a,"        RowBox[{\"If\", \"[\", ");
	mfprintf(a,"         RowBox[{");
	mfprintf(a,"          RowBox[{\"coloring\", \"\\[Equal]\", \"2\"}], \",\", ");
	mfprintf(a,"          RowBox[{");
	mfprintf(a,"           RowBox[{\"ColFunc\", \"[\", \"x_\", \"]\"}], \"=\", ");
	mfprintf(a,"           RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"            RowBox[{");
	mfprintf(a,"             RowBox[{\"1\", \"-\", \"x\"}], \",\", ");
	mfprintf(a,"             RowBox[{\"1\", \"-\", \"x\"}], \",\", ");
	mfprintf(a,"             RowBox[{\"1\", \"-\", \"x\"}]}], \"]\"}]}], \",\", ");
	mfprintf(a,"          RowBox[{\"If\", \"[\", ");
	mfprintf(a,"           RowBox[{");
	mfprintf(a,"            RowBox[{\"coloring\", \"\\[Equal]\", \"3\"}], \",\", ");
	mfprintf(a,"            RowBox[{");
	mfprintf(a,"             RowBox[{\"ColFunc\", \"[\", \"x_\", \"]\"}], \"=\", ");
	mfprintf(a,"             RowBox[{");
	mfprintf(a,"              RowBox[{\"ColorData\", \"[\", \"\\\"\\<TemperatureMap\\>\\\"\", \"]\"}], \"[\", ");
	mfprintf(a,"              \"x\", \"]\"}]}], \",\", ");
	mfprintf(a,"            RowBox[{");
	mfprintf(a,"             RowBox[{\"ColFunc\", \"[\", \"x_\", \"]\"}], \"=\", ");
	mfprintf(a,"             RowBox[{\"If\", \"[\", ");
	mfprintf(a,"              RowBox[{");
	mfprintf(a,"               RowBox[{\"x\", \"\\[LessEqual]\", \"0\"}], \",\", ");
	mfprintf(a,"               RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                RowBox[{\"0\", \",\", \"0\", \",\", \"0.5\"}], \"]\"}], \",\", ");
	mfprintf(a,"               RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                RowBox[{");
	mfprintf(a,"                 RowBox[{\"x\", \"<\", \"0.1\"}], \",\", ");
	mfprintf(a,"                 RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                  RowBox[{\"0\", \",\", \"0\", \",\", ");
	mfprintf(a,"                   RowBox[{\"0.5\", \"+\", ");
	mfprintf(a,"                    RowBox[{\"5\", \"*\", \"x\"}]}]}], \"]\"}], \",\", ");
	mfprintf(a,"                 RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                  RowBox[{");
	mfprintf(a,"                   RowBox[{\"x\", \"<\", \"0.4\"}], \",\", ");
	mfprintf(a,"                   RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{\"0\", \",\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"x\", \"-\", \"0.1\"}], \")\"}], \"/\", \"0.3\"}], \",\", ");
	mfprintf(a,"                    \"1\"}], \"]\"}], \",\", ");
	mfprintf(a,"                   RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"x\", \"<\", \"0.5\"}], \",\", ");
	mfprintf(a,"                    RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"x\", \"-\", \"0.4\"}], \")\"}], \"/\", \"0.1\"}], \",\", \"1\", ");
	mfprintf(a,"                    \",\", \"1\"}], \"]\"}], \",\", ");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"x\", \"<\", \"0.6\"}], \",\", ");
	mfprintf(a,"                    RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{\"1\", \",\", \"1\", \",\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"0.6\", \"-\", \"x\"}], \")\"}], \"*\", \"10\"}]}], \"]\"}], ");
	mfprintf(a,"                    \",\", ");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"x\", \"<\", \"0.8\"}], \",\", ");
	mfprintf(a,"                    RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{\"1\", \",\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"0.8\", \"-\", \"x\"}], \")\"}], \"/\", \"0.2\"}], \")\"}], ");
	mfprintf(a,"                    \"^\", \"0.5\"}], \",\", \"0\"}], \"]\"}], \",\", ");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"x\", \"<\", \"1\"}], \",\", ");
	mfprintf(a,"                    RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{\"1\", \",\", \"0\", \",\", ");
	mfprintf(a,"                    RowBox[{");
	mfprintf(a,"                    RowBox[{\"(\", ");
	mfprintf(a,"                    RowBox[{\"x\", \"-\", \"0.9\"}], \")\"}], \"*\", \"10\"}]}], \"]\"}], ");
	mfprintf(a,"                    \",\", ");
	mfprintf(a,"                    RowBox[{\"RGBColor\", \"[\", ");
	mfprintf(a,"                    RowBox[{\"1\", \",\", \"0\", \",\", \"1\"}], \"]\"}]}], \"]\"}]}], ");
	mfprintf(a,"                    \"]\"}]}], \"]\"}]}], \"]\"}]}], \"]\"}]}], \"]\"}]}], \"]\"}]}]}], ");
	mfprintf(a,"           \"]\"}]}], \"]\"}]}], \"]\"}], \";\", \"\\[IndentingNewLine]\", ");
	mfprintf(a,"      \"\\[IndentingNewLine]\", ");


	mfprintf(a,"      RowBox[{\"oplotrangeX1\", \"=\", \"%f\"}], \";\", \"\\[IndentingNewLine]\", \n",m_fMinVal[0]);
	mfprintf(a,"      RowBox[{\"oplotrangeY1\", \"=\", \"%f\"}], \";\", \"\\[IndentingNewLine]\", \n",m_fMinVal[1]);
	mfprintf(a,"      RowBox[{\"oplotrangeX2\", \"=\", \"%f\"}], \";\", \"\\[IndentingNewLine]\", \n",m_fMaxVal[0]);
	mfprintf(a,"      RowBox[{\"oplotrangeY2\", \"=\", \"%f\"}], \";\", \"\\[IndentingNewLine]\", \n",m_fMaxVal[1]);
	mfprintf(a,"      RowBox[{\"oplotresX\", \"=\", \"%d\"}], \";\", \"\\[IndentingNewLine]\", \n",m_iRes[0]);
	mfprintf(a,"      RowBox[{\"oplotresY\", \"=\", \"%d\"}], \";\", \"\\[IndentingNewLine]\", \n",m_iRes[1]);
	
	mfprintf(a,"      \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"plotresX\", \"=\", \n");
	mfprintf(a,"       RowBox[{\"Ceiling\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"plotrangeX2\", \"-\", \"plotrangeX1\"}], \")\"}], \"/\", \n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"oplotrangeX2\", \"-\", \"oplotrangeX1\"}], \")\"}]}], \"*\", \n");
	mfprintf(a,"         \"oplotresX\"}], \"]\"}]}], \";\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"plotresY\", \"=\", \n");
	mfprintf(a,"       RowBox[{\"Ceiling\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"plotrangeY2\", \"-\", \"plotrangeY1\"}], \")\"}], \"/\", \n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"oplotrangeY2\", \"-\", \"oplotrangeY1\"}], \")\"}]}], \"*\", \n");
	mfprintf(a,"         \"oplotresY\"}], \"]\"}]}], \";\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"plotresIndX\", \"=\", \n");
	mfprintf(a,"       RowBox[{\"Ceiling\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"plotrangeX1\", \"-\", \"oplotrangeX1\"}], \")\"}], \"/\", \n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"oplotrangeX2\", \"-\", \"oplotrangeX1\"}], \")\"}]}], \"*\", \n");
	mfprintf(a,"         \"oplotresX\"}], \"]\"}]}], \";\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"plotresIndY\", \"=\", \n");
	mfprintf(a,"       RowBox[{\"Ceiling\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"plotrangeY1\", \"-\", \"oplotrangeY1\"}], \")\"}], \"/\", \n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"oplotrangeY2\", \"-\", \"oplotrangeY1\"}], \")\"}]}], \"*\", \n");
	mfprintf(a,"         \"oplotresY\"}], \"]\"}]}], \";\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"       RowBox[{\"Smoothen\", \" \", \"Data\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"dat\", \"=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"ox\", \"=\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"             RowBox[{\"i\", \",\", \"plotresX\"}], \"]\"}], \"+\", \"plotresIndX\"}]}], \n");
	mfprintf(a,"          \";\", \n");
	mfprintf(a,"          RowBox[{\"oy\", \"=\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Quotient\", \"[\", \n");
	mfprintf(a,"             RowBox[{\"i\", \",\", \"plotresX\"}], \"]\"}], \"+\", \"plotresIndY\"}]}], \n");
	mfprintf(a,"          \";\", \n");
	mfprintf(a,"          RowBox[{\"x\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"            RowBox[{\"i\", \",\", \"plotresX\"}], \"]\"}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"y\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"Quotient\", \"[\", \n");
	mfprintf(a,"            RowBox[{\"i\", \",\", \"plotresX\"}], \"]\"}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"plotrangeX1\", \"+\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\"plotrangeX2\", \"-\", \"plotrangeX1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"              RowBox[{\"x\", \"/\", \n");
	mfprintf(a,"               RowBox[{\"(\", \n");
	mfprintf(a,"                RowBox[{\"plotresX\", \"-\", \"1\"}], \")\"}]}]}]}], \",\", \n");
	mfprintf(a,"            RowBox[{\"plotrangeY1\", \"+\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\"plotrangeY2\", \"-\", \"plotrangeY1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"              RowBox[{\"y\", \"/\", \n");
	mfprintf(a,"               RowBox[{\"(\", \n");
	mfprintf(a,"                RowBox[{\"plotresY\", \"-\", \"1\"}], \")\"}]}]}]}], \",\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"s\", \"=\", \"0\"}], \";\", \n");
	mfprintf(a,"             RowBox[{\"t\", \"=\", \"0\"}], \";\", \n");
	mfprintf(a,"             RowBox[{\"For\", \"[\", \n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"ty\", \"=\", \n");
	mfprintf(a,"                RowBox[{\"-\", \"smoothgrade\"}]}], \",\", \n");
	mfprintf(a,"               RowBox[{\"ty\", \"\\[LessEqual]\", \"smoothgrade\"}], \",\", \n");
	mfprintf(a,"               RowBox[{\"ty\", \"++\"}], \",\", \n");
	mfprintf(a,"               RowBox[{\"For\", \"[\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"tx\", \"=\", \n");
	mfprintf(a,"                  RowBox[{\"-\", \"smoothgrade\"}]}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"tx\", \"\\[LessEqual]\", \"smoothgrade\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"tx\", \"++\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"tx\", \"^\", \"2\"}], \"+\", \n");
	mfprintf(a,"                    RowBox[{\"ty\", \"^\", \"2\"}]}], \"\\[LessEqual]\", \n");
	mfprintf(a,"                    RowBox[{\"smoothgrade\", \"^\", \"2\"}]}], \",\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"r\", \"=\", \n");
	mfprintf(a,"                    RowBox[{\"1\", \"/\", \n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"tx\", \"^\", \"2\"}], \"+\", \n");
	mfprintf(a,"                    RowBox[{\"ty\", \"^\", \"2\"}], \"+\", \"1\"}], \")\"}]}]}], \";\", \n");
	mfprintf(a,"                    RowBox[{\"t\", \"+=\", \"r\"}], \";\", \n");
	mfprintf(a,"                    RowBox[{\"s\", \"+=\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat2\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Max\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Min\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ox\", \"+\", \"tx\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"oplotresX\", \"-\", \"1\"}]}], \"]\"}], \",\", \"0\"}], \n");
	mfprintf(a,"                    \"]\"}], \"+\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Max\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Min\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"oy\", \"+\", \"ty\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"oplotresY\", \"-\", \"1\"}]}], \"]\"}], \",\", \"0\"}], \n");
	mfprintf(a,"                    \"]\"}], \"*\", \"oplotresX\"}], \"+\", \"1\"}], \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"*\", \"r\"}]}]}], \",\"}], \n");
	mfprintf(a,"                  \"]\"}]}], \"]\"}]}], \"]\"}], \";\", \n");
	mfprintf(a,"             RowBox[{\"s\", \"/\", \"t\"}]}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"plotresX\", \"*\", \"plotresY\"}], \"-\", \"1\"}]}], \"}\"}]}], \n");
	mfprintf(a,"        \"]\"}]}], \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"       \"Add\", \" \", \"Functions\", \" \", \"to\", \" \", \"primitive\", \" \", \"list\"}], \n");
	mfprintf(a,"       \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"For\", \"[\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\"i\", \"=\", \"1\"}], \",\", \n");
	mfprintf(a,"        RowBox[{\"i\", \"<=\", \n");
	mfprintf(a,"         RowBox[{\"Length\", \"[\", \"funcX\", \"]\"}]}], \",\", \n");
	mfprintf(a,"        RowBox[{\"i\", \"++\"}], \",\", \n");
	mfprintf(a,"        RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"         RowBox[{\"prim\", \",\", \n");
	mfprintf(a,"          RowBox[{\"ListLinePlot\", \"[\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"x\", \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"funcX\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                 RowBox[{\"[\", \"1\", \"]\"}], \"]\"}]}], \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"x\", \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"funcX\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                  RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                 RowBox[{\"[\", \"1\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"funcX\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                  RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                 RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"funcX\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"-\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"funcX\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"1\", \"]\"}], \"]\"}]}], \")\"}], \"/\", \"1000\"}]}], \n");
	mfprintf(a,"               \"}\"}]}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"PlotStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"             RowBox[{\"Evaluate\", \"[\", \n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"funcX\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"               RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}]}]}], \"]\"}]}], \"]\"}]}], \n");
	mfprintf(a,"       \"]\"}], \";\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"For\", \"[\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\"i\", \"=\", \"1\"}], \",\", \n");
	mfprintf(a,"        RowBox[{\"i\", \"<=\", \n");
	mfprintf(a,"         RowBox[{\"Length\", \"[\", \"funcY\", \"]\"}]}], \",\", \n");
	mfprintf(a,"        RowBox[{\"i\", \"++\"}], \",\", \n");
	mfprintf(a,"        RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"         RowBox[{\"prim\", \",\", \n");
	mfprintf(a,"          RowBox[{\"ListLinePlot\", \"[\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"funcY\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                 RowBox[{\"[\", \"1\", \"]\"}], \"]\"}], \",\", \"y\"}], \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"y\", \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"funcY\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                  RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                 RowBox[{\"[\", \"1\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"funcY\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                  RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                 RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"funcY\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"-\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"funcY\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"1\", \"]\"}], \"]\"}]}], \")\"}], \"/\", \"1000\"}]}], \n");
	mfprintf(a,"               \"}\"}]}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"PlotStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"             RowBox[{\"Evaluate\", \"[\", \n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"funcY\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"               RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}]}]}], \"]\"}]}], \"]\"}]}], \n");
	mfprintf(a,"       \"]\"}], \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"tfx\", \":=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"i2\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"tickminX\", \"+\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"i\", \"/\", \n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"majorticksX\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksX\", \"+\", \"1\"}], \")\"}]}], \")\"}]}], \"*\", \n");
	mfprintf(a,"             RowBox[{\"(\", \n");
	mfprintf(a,"              RowBox[{\"tickmaxX\", \"-\", \"tickminX\"}], \")\"}]}]}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\"i2\", \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksX\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ticklabelprefixX\", \"<>\", \n");
	mfprintf(a,"               RowBox[{\"ToString\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"tickprecX\", \"\\[Equal]\", \"0\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"Round\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\"i2\", \",\", \"1\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"PaddedForm\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\"i2\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"20\", \",\", \"tickprecX\"}], \"}\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"NumberPadding\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"\\\"\\<\\>\\\"\", \",\", \"\\\"\\<0\\>\\\"\"}], \"}\"}]}]}], \n");
	mfprintf(a,"                   \"]\"}]}], \"]\"}], \"]\"}], \"<>\", \"ticklabelsuffixX\"}], \",\", \n");
	mfprintf(a,"              \"\\\"\\<\\>\\\"\"}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksX\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"majortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \"majortickshift\"}]}], \"}\"}], \n");
	mfprintf(a,"              \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"minortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \"minortickshift\"}]}], \"}\"}]}],\n");
	mfprintf(a,"              \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"majorticksX\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"minorticksX\", \"+\", \"1\"}], \")\"}]}]}], \"}\"}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"tfxnolabel\", \":=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"i2\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"tickminX\", \"+\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"i\", \"/\", \n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"majorticksX\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksX\", \"+\", \"1\"}], \")\"}]}], \")\"}]}], \"*\", \n");
	mfprintf(a,"             RowBox[{\"(\", \n");
	mfprintf(a,"              RowBox[{\"tickmaxX\", \"-\", \"tickminX\"}], \")\"}]}]}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\"i2\", \",\", \"\\\"\\<\\>\\\"\", \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksX\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"majortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \"majortickshift\"}]}], \"}\"}], \n");
	mfprintf(a,"              \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"minortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \"minortickshift\"}]}], \"}\"}]}],\n");
	mfprintf(a,"              \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"majorticksX\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"minorticksX\", \"+\", \"1\"}], \")\"}]}]}], \"}\"}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"tfy\", \":=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"i2\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"tickminY\", \"+\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"i\", \"/\", \n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"majorticksY\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksY\", \"+\", \"1\"}], \")\"}]}], \")\"}]}], \"*\", \n");
	mfprintf(a,"             RowBox[{\"(\", \n");
	mfprintf(a,"              RowBox[{\"tickmaxY\", \"-\", \"tickminY\"}], \")\"}]}]}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\"i2\", \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksY\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ticklabelprefixY\", \"<>\", \n");
	mfprintf(a,"               RowBox[{\"ToString\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"tickprecY\", \"\\[Equal]\", \"0\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"Round\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\"i2\", \",\", \"1\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"PaddedForm\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\"i2\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"20\", \",\", \"tickprecY\"}], \"}\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"NumberPadding\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"\\\"\\<\\>\\\"\", \",\", \"\\\"\\<0\\>\\\"\"}], \"}\"}]}]}], \n");
	mfprintf(a,"                   \"]\"}]}], \"]\"}], \"]\"}], \"<>\", \"ticklabelsuffixY\"}], \",\", \n");
	mfprintf(a,"              \"\\\"\\<\\>\\\"\"}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksY\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"majortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \"majortickshift\"}]}], \"}\"}], \n");
	mfprintf(a,"              \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"minortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \"minortickshift\"}]}], \"}\"}]}],\n");
	mfprintf(a,"              \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"majorticksY\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"minorticksY\", \"+\", \"1\"}], \")\"}]}]}], \"}\"}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"tfynolabel\", \":=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"i2\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"tickminY\", \"+\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"i\", \"/\", \n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"majorticksY\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksY\", \"+\", \"1\"}], \")\"}]}], \")\"}]}], \"*\", \n");
	mfprintf(a,"             RowBox[{\"(\", \n");
	mfprintf(a,"              RowBox[{\"tickmaxY\", \"-\", \"tickminY\"}], \")\"}]}]}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\"i2\", \",\", \"\\\"\\<\\>\\\"\", \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksY\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"majortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \"majortickshift\"}]}], \"}\"}], \n");
	mfprintf(a,"              \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"minortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \"minortickshift\"}]}], \"}\"}]}],\n");
	mfprintf(a,"              \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"majorticksY\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"minorticksY\", \"+\", \"1\"}], \")\"}]}]}], \"}\"}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"tfz\", \":=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"i2\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"tickminZ\", \"+\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"i\", \"/\", \n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"majorticksZ\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksZ\", \"+\", \"1\"}], \")\"}]}], \")\"}]}], \"*\", \n");
	mfprintf(a,"             RowBox[{\"(\", \n");
	mfprintf(a,"              RowBox[{\"tickmaxZ\", \"-\", \"tickminZ\"}], \")\"}]}]}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\"i2\", \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksZ\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ticklabelprefixZ\", \"<>\", \n");
	mfprintf(a,"               RowBox[{\"ToString\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"tickprecZ\", \"\\[Equal]\", \"0\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"Round\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"i2\", \"^\", \n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"1\", \"/\", \"exp\"}], \")\"}]}], \",\", \"1\"}], \"]\"}], \n");
	mfprintf(a,"                  \",\", \n");
	mfprintf(a,"                  RowBox[{\"PaddedForm\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"i2\", \"^\", \n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"1\", \"/\", \"exp\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"20\", \",\", \"tickprecZ\"}], \"}\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"NumberPadding\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"\\\"\\<\\>\\\"\", \",\", \"\\\"\\<0\\>\\\"\"}], \"}\"}]}]}], \n");
	mfprintf(a,"                   \"]\"}]}], \"]\"}], \"]\"}], \"<>\", \"ticklabelsuffixZ\"}], \",\", \n");
	mfprintf(a,"              \"\\\"\\<\\>\\\"\"}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksZ\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"majortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \"majortickshift\"}]}], \"}\"}], \n");
	mfprintf(a,"              \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"minortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \"minortickshift\"}]}], \"}\"}]}],\n");
	mfprintf(a,"              \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"majorticksZ\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"minorticksZ\", \"+\", \"1\"}], \")\"}]}]}], \"}\"}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"tfznolabel\", \":=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"i2\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"tickminZ\", \"+\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"i\", \"/\", \n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"majorticksZ\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksZ\", \"+\", \"1\"}], \")\"}]}], \")\"}]}], \"*\", \n");
	mfprintf(a,"             RowBox[{\"(\", \n");
	mfprintf(a,"              RowBox[{\"tickmaxZ\", \"-\", \"tickminZ\"}], \")\"}]}]}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\"i2\", \",\", \"\\\"\\<\\>\\\"\", \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksZ\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"majortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"majorticklength\", \"*\", \"majortickshift\"}]}], \"}\"}], \n");
	mfprintf(a,"              \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"minortickshift\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"minorticklength\", \"*\", \"minortickshift\"}]}], \"}\"}]}],\n");
	mfprintf(a,"              \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"majorticksZ\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"minorticksZ\", \"+\", \"1\"}], \")\"}]}]}], \"}\"}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"tfl\", \":=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"i2\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"tickminL\", \"+\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"i\", \"/\", \n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"majorticksL\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksL\", \"+\", \"1\"}], \")\"}]}], \")\"}]}], \"*\", \n");
	mfprintf(a,"             RowBox[{\"(\", \n");
	mfprintf(a,"              RowBox[{\"tickmaxL\", \"-\", \"tickminL\"}], \")\"}]}]}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\"i2\", \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksL\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ToString\", \"[\", \n");
	mfprintf(a,"               RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"tickprecL\", \"\\[Equal]\", \"0\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"Round\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"i2\", \",\", \"1\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"PaddedForm\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"i2\", \",\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"20\", \",\", \"tickprecL\"}], \"}\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"NumberPadding\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"\\\"\\<\\>\\\"\", \",\", \"\\\"\\<0\\>\\\"\"}], \"}\"}]}]}], \n");
	mfprintf(a,"                  \"]\"}]}], \"]\"}], \"]\"}], \",\", \"\\\"\\<\\>\\\"\"}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksL\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"majorticklengthL\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"majortickshiftL\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"majorticklengthL\", \"*\", \"majortickshiftL\"}]}], \"}\"}],\n");
	mfprintf(a,"               \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"minorticklengthL\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"minortickshiftL\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"minorticklengthL\", \"*\", \"minortickshiftL\"}]}], \n");
	mfprintf(a,"               \"}\"}]}], \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"majorticksL\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"minorticksL\", \"+\", \"1\"}], \")\"}]}]}], \"}\"}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"tflnolabel\", \":=\", \n");
	mfprintf(a,"       RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"i2\", \"=\", \n");
	mfprintf(a,"           RowBox[{\"tickminL\", \"+\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"i\", \"/\", \n");
	mfprintf(a,"              RowBox[{\"(\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"majorticksL\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                RowBox[{\"(\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksL\", \"+\", \"1\"}], \")\"}]}], \")\"}]}], \"*\", \n");
	mfprintf(a,"             RowBox[{\"(\", \n");
	mfprintf(a,"              RowBox[{\"tickmaxL\", \"-\", \"tickminL\"}], \")\"}]}]}]}], \";\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\"i2\", \",\", \"\\\"\\<\\>\\\"\", \",\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\n");
	mfprintf(a,"               RowBox[{\"Mod\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"i\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"minorticksL\", \"+\", \"1\"}]}], \"]\"}], \"\\[Equal]\", \n");
	mfprintf(a,"               \"0\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"majorticklengthL\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"majortickshiftL\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"majorticklengthL\", \"*\", \"majortickshiftL\"}]}], \"}\"}],\n");
	mfprintf(a,"               \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"minorticklengthL\", \"*\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\"1\", \"-\", \"minortickshiftL\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"minorticklengthL\", \"*\", \"minortickshiftL\"}]}], \n");
	mfprintf(a,"               \"}\"}]}], \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"i\", \",\", \"0\", \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"majorticksL\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\"(\", \n");
	mfprintf(a,"             RowBox[{\"minorticksL\", \"+\", \"1\"}], \")\"}]}]}], \"}\"}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"       RowBox[{\"Create\", \" \", \"Contour\", \" \", \"Plot\"}], \" \", \"*)\"}], \n");
	mfprintf(a,"      \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"datmin\", \"=\", \n");
	mfprintf(a,"       RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"ValueQ\", \"[\", \"mincolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"         \"mincolorvalue\", \",\", \n");
	mfprintf(a,"         RowBox[{\"Min\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"              RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"             RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\"i\", \",\", \"1\", \",\", \n");
	mfprintf(a,"              RowBox[{\"plotresX\", \"*\", \"plotresY\"}]}], \"}\"}]}], \"]\"}], \n");
	mfprintf(a,"          \"]\"}]}], \"]\"}]}], \";\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"datminexp\", \"=\", \n");
	mfprintf(a,"       RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"ValueQ\", \"[\", \"mincolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"          RowBox[{\"Abs\", \"[\", \n");
	mfprintf(a,"           RowBox[{\"mincolorvalue\", \"^\", \"exp\"}], \"]\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"Sign\", \"[\", \"datmin\", \"]\"}], \"*\", \n");
	mfprintf(a,"          RowBox[{\"Abs\", \"[\", \n");
	mfprintf(a,"           RowBox[{\"datmin\", \"^\", \"exp\"}], \"]\"}]}]}], \"]\"}]}], \";\", \n");
	mfprintf(a,"      \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"datmax\", \"=\", \n");
	mfprintf(a,"       RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"ValueQ\", \"[\", \"maxcolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"         \"maxcolorvalue\", \",\", \n");
	mfprintf(a,"         RowBox[{\"Max\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"              RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"             RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\"i\", \",\", \"1\", \",\", \n");
	mfprintf(a,"              RowBox[{\"plotresX\", \"*\", \"plotresY\"}]}], \"}\"}]}], \"]\"}], \n");
	mfprintf(a,"          \"]\"}]}], \"]\"}]}], \";\", \"\\[IndentingNewLine]\", \n");


	mfprintf(a,"      RowBox[{\"datmaxexp\", \"=\",\n");
	mfprintf(a,"       RowBox[{\"If\", \"[\",\n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"ValueQ\", \"[\", \"maxcolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"Sign\", \"[\", \"maxcolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"Abs\", \"[\", \"maxcolorvalue\", \"]\"}], \"^\", \"exp\"}]}], \",\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"Sign\", \"[\", \"datmax\", \"]\"}], \"*\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"Abs\", \"[\", \"datmax\", \"]\"}], \"^\", \"exp\"}]}]}], \"]\"}]}], \n");
	mfprintf(a,"      \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"      RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\"coloring\", \"\\[Equal]\", \"4\"}], \",\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Abs\", \"[\", \"datmin\", \"]\"}], \">\", \n");
	mfprintf(a,"            RowBox[{\"Abs\", \"[\", \"datmax\", \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"datmax\", \"=\", \n");
	mfprintf(a,"            RowBox[{\"Abs\", \"[\", \"datmin\", \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"datmin\", \"=\", \n");
	mfprintf(a,"            RowBox[{\"-\", \n");
	mfprintf(a,"             RowBox[{\"Abs\", \"[\", \"datmax\", \"]\"}]}]}]}], \"]\"}], \";\", \n");
	mfprintf(a,"         RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Abs\", \"[\", \"datminexp\", \"]\"}], \">\", \n");
	mfprintf(a,"            RowBox[{\"Abs\", \"[\", \"datmaxexp\", \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"datmaxexp\", \"=\", \n");
	mfprintf(a,"            RowBox[{\"Abs\", \"[\", \"datminexp\", \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"datminexp\", \"=\", \n");
	mfprintf(a,"            RowBox[{\"-\", \n");
	mfprintf(a,"             RowBox[{\"Abs\", \"[\", \"datmaxexp\", \"]\"}]}]}]}], \"]\"}]}], \",\"}], \n");
	mfprintf(a,"       \"]\"}], \";\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\",\n");


	mfprintf(a,"      RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\"plottype\", \"\\[Equal]\", \"1\"}], \",\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"        RowBox[{\"Show\", \"[\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"ListContourPlot\", \"[\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                 RowBox[{\"[\", \"1\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                 RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"Sign\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                   RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}], \"*\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"Abs\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}], \"^\", \"exp\"}]}]}], \n");
	mfprintf(a,"               \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"i\", \",\", \"1\", \",\", \n");
	mfprintf(a,"                RowBox[{\"plotresX\", \"*\", \"plotresY\"}]}], \"}\"}]}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"PlotRange\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"{\", \n");
	mfprintf(a,"                RowBox[{\"datminexp\", \",\", \"datmaxexp\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ColorFunction\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"Function\", \"[\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"{\", \"z\", \"}\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"ColFunc\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"z\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"mincolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"-\", \n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}]}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"mincolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"-\", \"datminexp\"}]}], \"]\"}]}], \")\"}], \"/\", \n");
	mfprintf(a,"                   RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"maxcolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"maxcolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"maxcolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}], \"-\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"mincolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"datmaxexp\", \"-\", \"datminexp\"}]}], \"]\"}], \"/\", \n");
	mfprintf(a,"                    \"colorscale\"}], \")\"}]}], \"]\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ColorFunctionScaling\", \"\\[Rule]\", \"False\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"Contours\", \"\\[Rule]\", \"contours\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"PlotRangePadding\", \"\\[Rule]\", \"framemargins\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"FrameStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                RowBox[{\"Thickness\", \"[\", \"framethickness\", \"]\"}], \"]\"}]}], \n");
	mfprintf(a,"              \",\", \n");
	mfprintf(a,"              RowBox[{\"FrameTicks\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"{\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\"tfy\", \",\", \n");
	mfprintf(a,"                   RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"ticksbothsidesY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"labelsbothsidesY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"tfy\", \",\", \"tfynolabel\"}], \"]\"}], \",\", \"None\"}], \"]\"}]}],\n");
	mfprintf(a,"                   \"}\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\"tfx\", \",\", \n");
	mfprintf(a,"                   RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"ticksbothsidesX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"labelsbothsidesX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"tfx\", \",\", \"tfxnolabel\"}], \"]\"}], \",\", \"None\"}], \"]\"}]}],\n");
	mfprintf(a,"                   \"}\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"              RowBox[{\"FrameLabel\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"{\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"xlabel\", \",\", \n");
	mfprintf(a,"                   RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"ylabel\", \",\", \n");
	mfprintf(a,"                   RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                 \"\\\"\\<\\>\\\"\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ContourLabels\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"ValueQ\", \"[\", \"usecontourlabels\", \"]\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"(\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"Text\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"#3\", \"^\", \n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"1\", \"/\", \"exp\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontFamily\", \"\\[Rule]\", \"\\\"\\<Arial\\>\\\"\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"contourlabelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"#1\", \",\", \"#2\"}], \"}\"}]}], \"]\"}], \"&\"}], \")\"}], \n");
	mfprintf(a,"                 \",\", \"None\"}], \"]\"}]}], \",\", \n");
	mfprintf(a,"              RowBox[{\"InterpolationOrder\", \"\\[Rule]\", \"interpolationorder\"}],\n");
	mfprintf(a,"               \",\", \n");
	mfprintf(a,"              RowBox[{\"ContourStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"Thickness\", \"[\", \"contourthickness\", \"]\"}], \",\", \n");
	mfprintf(a,"                 \"contourcolor\", \",\", \n");
	mfprintf(a,"                 RowBox[{\"Opacity\", \"[\", \"contouropacity\", \"]\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"Dashing\", \"[\", \"contourdash\", \"]\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"              RowBox[{\"LabelStyle\", \"\\[Rule]\", \n"); /*********************************************/
	mfprintf(a,"               RowBox[{\"{\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"FontFamily\", \"\\[Rule]\", \"\\\"\\<Arial\\>\\\"\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                  RowBox[{\"ticklabelsize\", \"*\", \n");
	mfprintf(a,"                   RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"              RowBox[{\"AspectRatio\", \"\\[Rule]\", \"aspectratio\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ImageSize\", \"\\[Rule]\", \"plotpixels\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"drawmesh\", \"\\[Equal]\", \"1\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"{\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"Mesh\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"majorticksX\", \"-\", \"2\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"majorticksX\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                    \"minorticksX\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"majorticksY\", \"-\", \"2\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"majorticksY\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                    \"minorticksY\"}]}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"MeshFunctions\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"#1\", \"&\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"#2\", \"&\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"MeshStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Thickness\", \"[\", \"meshthicknessX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Opacity\", \"[\", \"meshopacityX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"meshcolorX\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Thickness\", \"[\", \"meshthicknessY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Opacity\", \"[\", \"meshopacityY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"meshcolorY\"}], \"]\"}]}], \"}\"}]}]}], \"}\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"Mesh\", \"\\[Rule]\", \"None\"}]}], \"]\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"PerformanceGoal\", \"\\[Rule]\", \"\\\"\\<Quality\\>\\\"\"}]}], \n");
	mfprintf(a,"             \"}\"}]}], \"]\"}], \",\", \"\\[IndentingNewLine]\", \"prim\"}], \"]\"}], \",\",\n");
	mfprintf(a,"         \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"        RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"plottype\", \"\\[Equal]\", \"2\"}], \",\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"          RowBox[{\"Show\", \"[\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"ListDensityPlot\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"{\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                   RowBox[{\"[\", \"1\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                   RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"Sign\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}], \"*\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}], \"^\", \"exp\"}]}]}], \n");
	mfprintf(a,"                 \"}\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"{\", \n");
	mfprintf(a,"                 RowBox[{\"i\", \",\", \"1\", \",\", \n");
	mfprintf(a,"                  RowBox[{\"plotresX\", \"*\", \"plotresY\"}]}], \"}\"}]}], \"]\"}], \n");
	mfprintf(a,"              \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"PlotRange\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\"datminexp\", \",\", \"datmaxexp\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"ColorFunction\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"Function\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"{\", \"z\", \"}\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"ColFunc\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"z\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"mincolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"-\", \n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}]}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"mincolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"-\", \"datminexp\"}]}], \"]\"}]}], \")\"}], \"/\", \n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"maxcolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"maxcolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"maxcolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}], \"-\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"mincolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"datmaxexp\", \"-\", \"datminexp\"}]}], \"]\"}], \"/\", \n");
	mfprintf(a,"                    \"colorscale\"}], \")\"}]}], \"]\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"ColorFunctionScaling\", \"\\[Rule]\", \"False\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"PlotRangePadding\", \"\\[Rule]\", \"framemargins\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"FrameStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"Thickness\", \"[\", \"framethickness\", \"]\"}], \"]\"}]}], \n");
	mfprintf(a,"                \",\", \n");
	mfprintf(a,"                RowBox[{\"FrameTicks\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"tfy\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"ticksbothsidesY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"labelsbothsidesY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"tfy\", \",\", \"tfynolabel\"}], \"]\"}], \",\", \"None\"}], \"]\"}]}],\n");
	mfprintf(a,"                     \"}\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"tfx\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"ticksbothsidesX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"labelsbothsidesX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"tfx\", \",\", \"tfxnolabel\"}], \"]\"}], \",\", \"None\"}], \"]\"}]}],\n");
	mfprintf(a,"                     \"}\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"FrameLabel\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"xlabel\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"ylabel\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                   \"\\\"\\<\\>\\\"\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                \"InterpolationOrder\", \"\\[Rule]\", \"interpolationorder\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"LabelStyle\", \"\\[Rule]\", \n");  /*********************************************/
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"FontFamily\", \"\\[Rule]\", \"\\\"\\<Arial\\>\\\"\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"ticklabelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"AspectRatio\", \"\\[Rule]\", \"aspectratio\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"ImageSize\", \"\\[Rule]\", \"plotpixels\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"drawmesh\", \"\\[Equal]\", \"1\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"{\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Mesh\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"majorticksX\", \"-\", \"2\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"majorticksX\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                    \"minorticksX\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"majorticksY\", \"-\", \"2\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"majorticksY\", \"-\", \"1\"}], \")\"}], \"*\", \n");
	mfprintf(a,"                    \"minorticksY\"}]}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"MeshFunctions\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"#1\", \"&\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"#2\", \"&\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"MeshStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Thickness\", \"[\", \"meshthicknessX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Opacity\", \"[\", \"meshopacityX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"meshcolorX\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Thickness\", \"[\", \"meshthicknessY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Opacity\", \"[\", \"meshopacityY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"meshcolorY\"}], \"]\"}]}], \"}\"}]}]}], \"}\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"Mesh\", \"\\[Rule]\", \"None\"}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"PerformanceGoal\", \"\\[Rule]\", \"\\\"\\<Quality\\>\\\"\"}]}], \n");
	mfprintf(a,"               \"}\"}]}], \"]\"}], \",\", \"\\[IndentingNewLine]\", \"prim\"}], \"]\"}], \n");
	mfprintf(a,"          \",\", \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"          RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"plottype\", \"\\[Equal]\", \"3\"}], \",\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"            RowBox[{\"ReliefPlot\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"Sign\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"i\", \"*\", \"plotresX\"}], \"+\", \"j\"}], \"]\"}], \"]\"}], \n");
	mfprintf(a,"                   \"[\", \n");
	mfprintf(a,"                   RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}], \"*\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"Abs\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"i\", \"*\", \"plotresX\"}], \"+\", \"j\"}], \"]\"}], \"]\"}], \n");
	mfprintf(a,"                    \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}], \"^\", \"exp\"}]}], \n");
	mfprintf(a,"                \",\", \n");
	mfprintf(a,"                RowBox[{\"{\", \n");
	mfprintf(a,"                 RowBox[{\"i\", \",\", \"1\", \",\", \n");
	mfprintf(a,"                  RowBox[{\"plotresY\", \"-\", \"1\"}]}], \"}\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"{\", \n");
	mfprintf(a,"                 RowBox[{\"j\", \",\", \"1\", \",\", \n");
	mfprintf(a,"                  RowBox[{\"plotresX\", \"-\", \"1\"}]}], \"}\"}]}], \"]\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"PlotRange\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\"datminexp\", \",\", \"datmaxexp\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"ColorFunction\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"Function\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"{\", \"z\", \"}\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"ColFunc\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"z\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"mincolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"-\", \n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}]}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"mincolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"-\", \"datminexp\"}]}], \"]\"}]}], \")\"}], \"/\", \n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"maxcolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"maxcolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"maxcolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}], \"-\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"mincolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"datmaxexp\", \"-\", \"datminexp\"}]}], \"]\"}], \"/\", \n");
	mfprintf(a,"                    \"colorscale\"}], \")\"}]}], \"]\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"ColorFunctionScaling\", \"\\[Rule]\", \"False\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"PlotRangePadding\", \"\\[Rule]\", \"framemargins\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"FrameStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                  RowBox[{\"Thickness\", \"[\", \"framethickness\", \"]\"}], \"]\"}]}], \n");
	mfprintf(a,"                \",\", \n");
	mfprintf(a,"                RowBox[{\"FrameTicks\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"tfy\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"ticksbothsidesY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"labelsbothsidesY\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"tfy\", \",\", \"tfynolabel\"}], \"]\"}], \",\", \"None\"}], \"]\"}]}],\n");
	mfprintf(a,"                     \"}\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"tfx\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"ticksbothsidesX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"labelsbothsidesX\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"tfx\", \",\", \"tfxnolabel\"}], \"]\"}], \",\", \"None\"}], \"]\"}]}],\n");
	mfprintf(a,"                     \"}\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"FrameLabel\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"xlabel\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"ylabel\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                   \"\\\"\\<\\>\\\"\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"AspectRatio\", \"\\[Rule]\", \"aspectratio\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"ImageSize\", \"\\[Rule]\", \"plotpixels\"}], \",\", \n");
	mfprintf(a,"                RowBox[{\"DataRange\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"plotrangeX1\", \",\", \"plotrangeX2\"}], \"}\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"plotrangeY1\", \",\", \"plotrangeY2\"}], \"}\"}]}], \n");
	mfprintf(a,"                  \"}\"}]}], \",\", \n");
	mfprintf(a,"                RowBox[{\"LabelStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"{\", \n");
	mfprintf(a,"                  RowBox[{\n");
	mfprintf(a,"                   RowBox[{\"FontFamily\", \"\\[Rule]\", \"\\\"\\<Arial\\>\\\"\"}], \",\", \n");
	mfprintf(a,"                   RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"ticklabelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"}\"}]}]}], \n");
	mfprintf(a,"               \"}\"}]}], \"]\"}], \",\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"            \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"plottype\", \"\\[Equal]\", \"4\"}], \",\", \n");
	mfprintf(a,"              \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"              RowBox[{\"ListPlot3D\", \"[\", \n");
	mfprintf(a,"               RowBox[{\n");
	mfprintf(a,"                RowBox[{\"Table\", \"[\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"{\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"1\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"2\", \"]\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"dat\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"i\", \"]\"}], \"]\"}], \"[\", \n");
	mfprintf(a,"                    RowBox[{\"[\", \"3\", \"]\"}], \"]\"}], \"]\"}], \"^\", \"exp\"}]}]}], \n");
	mfprintf(a,"                   \"}\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"{\", \n");
	mfprintf(a,"                   RowBox[{\"i\", \",\", \"1\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"plotresX\", \"*\", \"plotresY\"}]}], \"}\"}]}], \"]\"}], \n");
	mfprintf(a,"                \",\", \n");
	mfprintf(a,"                RowBox[{\"{\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"PlotRange\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"datminexp\", \",\", \"datmaxexp\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"ColorFunction\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"Function\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"x\", \",\", \"y\", \",\", \"z\"}], \"}\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"ColFunc\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\"z\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"mincolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"-\", \n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}]}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"mincolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"-\", \"datminexp\"}]}], \"]\"}]}], \")\"}], \"/\", \n");
	mfprintf(a,"                    RowBox[{\"(\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"maxcolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"maxcolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"maxcolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}], \"-\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Sign\", \"[\", \"mincolorvalue\", \"]\"}], \"*\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Abs\", \"[\", \"mincolorvalue\", \"]\"}], \"^\", \n");
	mfprintf(a,"                    \"exp\"}]}]}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"datmaxexp\", \"-\", \"datminexp\"}]}], \"]\"}], \"/\", \n");
	mfprintf(a,"                    \"colorscale\"}], \")\"}]}], \"]\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"ColorFunctionScaling\", \"\\[Rule]\", \"False\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"PlotRangePadding\", \"\\[Rule]\", \"framemargins\"}], \n");
	mfprintf(a,"                  \",\", \n");
	mfprintf(a,"                  RowBox[{\"InterpolationOrder\", \"\\[Rule]\", \"2\"}], \",\", \n");

	mfprintf(a,"                  RowBox[{\"AspectRatio\", \"\\[Rule]\", \"aspectratio\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"ImageSize\", \"\\[Rule]\", \"plotpixels\"}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"Mesh\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"gridxcount\", \",\", \"gridycount\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"MeshStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Thickness\", \"[\", \"gridxthickness\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Opacity\", \"[\", \"gridxopacity\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"gridxcolor\"}], \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Thickness\", \"[\", \"gridythickness\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Opacity\", \"[\", \"gridyopacity\", \"]\"}], \",\", \n");
	mfprintf(a,"                    \"gridycolor\"}], \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"Boxed\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"boxed\", \"]\"}], \",\", \"True\", \",\", \n");
	mfprintf(a,"                    \"False\"}], \"]\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"BoxStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"boxcolor\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"Thickness\", \"[\", \"boxthickness\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Opacity\", \"[\", \"boxopacity\", \"]\"}]}], \"]\"}]}], \n");
	mfprintf(a,"                  \",\", \n");
	mfprintf(a,"                  RowBox[{\"Ticks\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\"tfx\", \",\", \"tfy\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"labelsz\", \"]\"}], \",\", \"tfz\", \",\", \n");
	mfprintf(a,"                    \"tfznolabel\"}], \"]\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"TicksStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"Thickness\", \"[\", \"axesthickness\", \"]\"}], \"]\"}]}],\n");
	mfprintf(a,"                   \",\", \n");
	mfprintf(a,"                  RowBox[{\"AxesStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"Thickness\", \"[\", \"axesthickness\", \"]\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"AxesEdge\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    \"xaxisposition\", \",\", \"yaxisposition\", \",\", \n");
	mfprintf(a,"                    \"zaxisposition\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"AxesLabel\", \"\\[Rule]\", \n");
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"xlabel\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"ylabel\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\"zlabel\", \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}]}], \n");
	mfprintf(a,"                    \"}\"}]}], \",\", \n");
	mfprintf(a,"                  RowBox[{\"LabelStyle\", \"\\[Rule]\", \n");   /*********************************************/
	mfprintf(a,"                   RowBox[{\"{\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"FontFamily\", \"\\[Rule]\", \"\\\"\\<Arial\\>\\\"\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                    RowBox[{\"ticklabelsize\", \"*\", \n");
	mfprintf(a,"                    RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"}\"}]}]}], \n");
	mfprintf(a,"                 \"}\"}]}], \"]\"}], \"\\[IndentingNewLine]\", \",\"}], \"]\"}]}], \n");
	mfprintf(a,"           \"]\"}]}], \"]\"}]}], \"]\"}]}], \"\\[IndentingNewLine]\", \")\"}]}], \";\"}], \n");
	mfprintf(a,"  \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"Create\", \" \", \"Legend\"}], \" \", \"*)\"}]}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"plotlegend\", \":=\", \n");
	mfprintf(a,"   RowBox[{\"(\", \n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\"x\", \"=.\"}], \";\", \n");
	mfprintf(a,"     RowBox[{\"y\", \"=.\"}], \";\", \n");
	mfprintf(a,"     RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"      RowBox[{\n");
	mfprintf(a,"       RowBox[{\"plottype\", \">\", \"1\"}], \",\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"       \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"       RowBox[{\"DensityPlot\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"explegend\", \"==\", \"1\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Sign\", \"[\", \"x\", \"]\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"Abs\", \"[\", \"x\", \"]\"}], \"^\", \"exp\"}]}], \",\", \"x\"}], \n");
	mfprintf(a,"          \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"x\", \",\", \"datmin\", \",\", \"datmax\"}], \"}\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"y\", \",\", \"0\", \",\", \"1\"}], \"}\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"PlotRange\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"explegend\", \"==\", \"1\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"datminexp\", \",\", \"datmaxexp\"}], \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"datmin\", \",\", \"datmax\"}], \"}\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"ColorFunction\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"Function\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"{\", \"z\", \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ColFunc\", \"[\", \n");
	mfprintf(a,"               RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"explegend\", \"==\", \"1\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"(\", \n");
	mfprintf(a,"                   RowBox[{\"z\", \" \", \"-\", \"datminexp\"}], \")\"}], \"/\", \n");
	mfprintf(a,"                  RowBox[{\"(\", \n");
	mfprintf(a,"                   RowBox[{\"datmaxexp\", \"-\", \"datminexp\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"(\", \n");
	mfprintf(a,"                   RowBox[{\"z\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"mincolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"-\", \"mincolorvalue\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"-\", \"datmin\"}]}], \"]\"}]}], \")\"}], \"/\", \n");
	mfprintf(a,"                  RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"maxcolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"maxcolorvalue\", \"-\", \"mincolorvalue\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"datmax\", \"-\", \"datmin\"}]}], \"]\"}]}]}], \"]\"}], \n");
	mfprintf(a,"               \"]\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"ColorFunctionScaling\", \"\\[Rule]\", \"False\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"Axes\", \"\\[Rule]\", \"False\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"Frame\", \"\\[Rule]\", \"True\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"FrameStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"             RowBox[{\"Thickness\", \"[\", \"framethickness\", \"]\"}], \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"FrameTicks\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"None\", \",\", \"None\"}], \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"tfl\", \",\", \"None\"}], \"}\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"PlotRangePadding\", \"\\[Rule]\", \"None\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"FrameLabel\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"               RowBox[{\"zlabel\", \",\", \n");
	mfprintf(a,"                RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                  RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"              \"\\\"\\<\\>\\\"\", \",\", \"\\\"\\<\\>\\\"\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"LabelStyle\", \"\\[Rule]\", \n");   /*********************************************/
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"FontFamily\", \"\\[Rule]\", \"\\\"\\<Arial\\>\\\"\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"ticklabelsize\", \"*\", \n");
	mfprintf(a,"                RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"AspectRatio\", \"\\[Rule]\", \"0.04\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"ImageSize\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"plotpixels\", \"*\", \n");
	mfprintf(a,"             RowBox[{\"7\", \"/\", \"8\"}]}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"PlotPoints\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\"2000\", \",\", \"2\"}], \"}\"}]}]}], \"}\"}]}], \"]\"}], \",\", \n");
	mfprintf(a,"       \"\\[IndentingNewLine]\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"       RowBox[{\"ContourPlot\", \"[\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"explegend\", \"==\", \"1\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Sign\", \"[\", \"x\", \"]\"}], \"*\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"Abs\", \"[\", \"x\", \"]\"}], \"^\", \"exp\"}]}], \",\", \"x\"}], \n");
	mfprintf(a,"          \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"x\", \",\", \"datmin\", \",\", \"datmax\"}], \"}\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\"y\", \",\", \"0\", \",\", \"1\"}], \"}\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"{\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"PlotRange\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"explegend\", \"==\", \"1\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"datminexp\", \",\", \"datmaxexp\"}], \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"datmin\", \",\", \"datmax\"}], \"}\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"ColorFunction\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"Function\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"{\", \"z\", \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"ColFunc\", \"[\", \n");
	mfprintf(a,"               RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                RowBox[{\n");
	mfprintf(a,"                 RowBox[{\"explegend\", \"==\", \"1\"}], \",\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"(\", \n");
	mfprintf(a,"                   RowBox[{\"z\", \" \", \"-\", \"datminexp\"}], \")\"}], \"/\", \n");
	mfprintf(a,"                  RowBox[{\"(\", \n");
	mfprintf(a,"                   RowBox[{\"datmaxexp\", \"-\", \"datminexp\"}], \")\"}]}], \",\", \n");
	mfprintf(a,"                 RowBox[{\n");
	mfprintf(a,"                  RowBox[{\"(\", \n");
	mfprintf(a,"                   RowBox[{\"z\", \"+\", \n");
	mfprintf(a,"                    RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                    RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"mincolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"-\", \"mincolorvalue\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"-\", \"datmin\"}]}], \"]\"}]}], \")\"}], \"/\", \n");
	mfprintf(a,"                  RowBox[{\"If\", \"[\", \n");
	mfprintf(a,"                   RowBox[{\n");
	mfprintf(a,"                    RowBox[{\"ValueQ\", \"[\", \"maxcolorvalue\", \"]\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"maxcolorvalue\", \"-\", \"mincolorvalue\"}], \",\", \n");
	mfprintf(a,"                    RowBox[{\"datmax\", \"-\", \"datmin\"}]}], \"]\"}]}]}], \"]\"}], \n");
	mfprintf(a,"               \"]\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"ColorFunctionScaling\", \"\\[Rule]\", \"False\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"Axes\", \"\\[Rule]\", \"False\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"Frame\", \"\\[Rule]\", \"True\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"FrameStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"             RowBox[{\"Thickness\", \"[\", \"framethickness\", \"]\"}], \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"Contours\", \"\\[Rule]\", \"contours\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"FrameTicks\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"None\", \",\", \"None\"}], \"}\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"{\", \n");
	mfprintf(a,"               RowBox[{\"tfl\", \",\", \"None\"}], \"}\"}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"PlotRangePadding\", \"\\[Rule]\", \"None\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"ContourStyle\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"Directive\", \"[\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"Thickness\", \"[\", \"contourthickness\", \"]\"}], \",\", \n");
	mfprintf(a,"              \"contourcolor\", \",\", \n");
	mfprintf(a,"              RowBox[{\"Opacity\", \"[\", \"contouropacity\", \"]\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"Dashing\", \"[\", \"contourdash\", \"]\"}]}], \"]\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"FrameLabel\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"Style\", \"[\", \n");
	mfprintf(a,"               RowBox[{\"zlabel\", \",\", \n");
	mfprintf(a,"                RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"                 RowBox[{\"labelsize\", \"*\", \n");
	mfprintf(a,"                  RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"]\"}], \",\", \n");
	mfprintf(a,"              \"\\\"\\<\\>\\\"\", \",\", \"\\\"\\<\\>\\\"\"}], \"}\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"LabelStyle\", \"\\[Rule]\", \n");   /*********************************************/
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\n");
	mfprintf(a,"              RowBox[{\"FontFamily\", \"\\[Rule]\", \"\\\"\\<Arial\\>\\\"\"}], \",\", \n");
	mfprintf(a,"              RowBox[{\"FontSize\", \"\\[Rule]\", \n");
	mfprintf(a,"               RowBox[{\"ticklabelsize\", \"*\", \n");
	mfprintf(a,"                RowBox[{\"plotpixels\", \"/\", \"800\"}]}]}]}], \"}\"}]}], \",\", \n");
	mfprintf(a,"           RowBox[{\"AspectRatio\", \"\\[Rule]\", \"0.04\"}], \",\", \n");
	mfprintf(a,"           RowBox[{\"ImageSize\", \"\\[Rule]\", \n");
	mfprintf(a,"            RowBox[{\"plotpixels\", \"*\", \n");
	mfprintf(a,"             RowBox[{\"7\", \"/\", \"8\"}]}]}]}], \"}\"}]}], \"]\"}]}], \"]\"}]}], \n");
	mfprintf(a,"    \")\"}]}], \";\"}]}], \"Input\" ]\n");
	mfprintf(a,"}, Closed]],\n");



	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Plotting Parameters\", \"Section\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell],\n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Misc. Properties\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = True}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\"Choose\", \" \", \"the\", \" \", \"plot\", \" \", \n");
	mfprintf(a,"      RowBox[{\"type\", \".\", \" \", \"1\"}]}], \" \", \"=\", \" \", \n");
	mfprintf(a,"     RowBox[{\"Contour\", \" \", \"Plot\"}]}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"2\", \" \", \"=\", \" \", \n");
	mfprintf(a,"     RowBox[{\"Density\", \" \", \"Plot\"}]}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"3\", \" \", \"=\", \" \", \n");
	mfprintf(a,"     RowBox[{\"Relief\", \" \", \"Plot\"}]}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"4\", \" \", \"=\", \" \", \n");
	mfprintf(a,"     RowBox[{\"3\", \"D\", \" \", \"Plot\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"plottype\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",m_iPlotType);
	mfprintf(a,"\n");


	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\"Choose\", \" \", \"the\", \" \", \"color\", \" \", \n");
	mfprintf(a,"      RowBox[{\"scale\", \":\", \"1\"}]}], \"=\", \n");
	mfprintf(a,"     RowBox[{\"Travis\", \" \", \"Scale\"}]}], \",\", \n");
	mfprintf(a,"    RowBox[{\"2\", \"=\", \n");
	mfprintf(a,"     RowBox[{\"black\", \" \", \"and\", \" \", \"white\"}]}], \",\", \n");
	mfprintf(a,"    RowBox[{\"3\", \"=\", \n");
	mfprintf(a,"     RowBox[{\"Temperature\", \" \", \"scale\"}]}], \",\", \n");
	mfprintf(a,"    RowBox[{\"4\", \"=\", \n");
	mfprintf(a,"     RowBox[{\"PlusMinus\", \" \", \"scale\"}]}]}], \" \", \"*)\"}], \n");
	mfprintf(a,"  \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"coloring\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",m_iColorScale);

	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Choose\", \" \", \"the\", \" \", \"exponent\", \" \", \"your\", \" \", \"data\", \" \", \n");
	mfprintf(a,"    \"will\", \" \", \"be\", \" \", \"processed\", \" \", \n");
	mfprintf(a,"    RowBox[{\"with\", \".\", \" \", \"Try\"}], \" \", \"values\", \" \", \"between\", \" \", \n");
	mfprintf(a,"    \"0.2\", \" \", \"and\", \" \", \"1.\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"exp\", \"=\", \"%.1f\"}], \";\"}]}]], \"Input\"],\n",m_fPlotExp);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Also\", \" \", \"scale\", \" \", \"the\", \" \", \"legend\", \" \", \"with\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"exponent\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"explegend\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",m_iExpLegend);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Choose\", \" \", \"the\", \" \", \"smoothing\", \" \", \"grade\", \" \", \"for\", \" \", \n");
	mfprintf(a,"     \"your\", \" \", \n");
	mfprintf(a,"     RowBox[{\"data\", \".\", \" \", \"The\"}], \" \", \"higher\", \" \", \"it\", \" \", \"is\"}],\n");
	mfprintf(a,"     \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"the\", \" \", \"longer\", \" \", \"it\", \" \", \n");
	mfprintf(a,"     RowBox[{\"takes\", \".\", \" \", \"0\"}], \" \", \"means\", \" \", \"no\", \" \", \n");
	mfprintf(a,"     RowBox[{\"smoothing\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"smoothgrade\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",m_iSmoothGrade);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"interpolation\", \" \", \"order\", \" \", \"mathematica\", \" \", \"uses\",\n");
	mfprintf(a,"     \" \", \"between\", \" \", \"the\", \" \", \"data\", \" \", \n");
	mfprintf(a,"    RowBox[{\"points\", \".\", \" \", \"0\"}], \" \", \"shows\", \" \", \"the\", \" \", \"data\", \n");
	mfprintf(a,"    \" \", \"points\", \" \", \"like\", \" \", \"they\", \" \", \n");
	mfprintf(a,"    RowBox[{\"are\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"interpolationorder\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",m_iInterpolationOrder);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"Choose\", \" \", \"the\", \" \", \"color\", \" \", \n");
	mfprintf(a,"    RowBox[{\"scaling\", \".\", \" \", \"1\"}], \" \", \"means\", \" \", \"that\", \" \", \"the\",\n");
	mfprintf(a,"     \" \", \"highest\", \" \", \"peak\", \" \", \"only\", \" \", \"gets\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    \"highest\", \" \", \n");
	mfprintf(a,"    RowBox[{\"color\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"colorscale\", \"=\", \"1.0\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"width\", \" \", \"of\", \" \", \"the\", \" \", \"Plot\", \" \", \"in\", \" \", \n");
	mfprintf(a,"    \"pixels\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"plotpixels\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",m_iPlotPixel);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"The\", \" \", \"aspect\", \" \", \"ratio\", \" \", \"of\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    RowBox[{\"plot\", \".\", \" \", \"1\"}], \" \", \"means\", \" \", \n");
	mfprintf(a,"    RowBox[{\"quadratic\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"aspectratio\", \"=\", \"%f\"}], \";\"}]}]], \"Input\"],\n",m_fAspectRatio);

	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"     \"Controls\", \" \", \"if\", \" \", \"a\", \" \", \"mesh\", \" \", \"is\", \" \", \"drawn\", \n");
	mfprintf(a,"      \" \", \"across\", \" \", \"the\", \" \", \n");
	mfprintf(a,"      RowBox[{\"plot\", \".\", \" \", \"1\"}]}], \"=\", \"enable\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"0\", \"=\", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\n");
	mfprintf(a,"       RowBox[{\"disable\", \".\", \" \", \"For\"}], \" \", \"mesh\", \" \", \"appearance\", \n");
	mfprintf(a,"       \" \", \"see\", \" \", \"\\\"\\<Axis Properties\\>\\\"\"}], \" \", \"\\[Rule]\", \" \", \n");
	mfprintf(a,"      \"\\\"\\<Mesh Properties\\>\\\"\"}]}]}], \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"drawmesh\", \"=\", \"1\"}], \";\"}]}]], \"Input\"]\n");

	mfprintf(a,"}, Open  ]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Contour Properties\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Choose\", \" \", \"the\", \" \", \"number\", \" \", \"of\", \" \", \"contours\", \" \", \n");
	mfprintf(a,"    \"that\", \" \", \"will\", \" \", \"be\", \" \", \"displayed\"}], \" \", \"*)\"}], \n");
	mfprintf(a,"  \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"contours\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",g_pDatabase->GetInt("/PLOT2D/DEFAULTS/CONTOUR_LINES"));
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Print\", \" \", \"the\", \" \", \"height\", \" \", \"value\", \" \", \"onto\", \" \", \n");
	mfprintf(a,"     \"each\", \" \", \"contour\", \" \", \n");
	mfprintf(a,"     RowBox[{\"line\", \"?\", \" \", \"If\"}], \" \", \"desired\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"replace\", \" \", \"the\", \" \", \"dot\", \" \", \"by\", \" \", \n");
	mfprintf(a,"     RowBox[{\"\\\"\\<1\\>\\\"\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"usecontourlabels\", \"=.\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"If\", \" \", \"contour\", \" \", \"labels\", \" \", \"are\", \" \", \n");
	mfprintf(a,"    RowBox[{\"used\", \":\", \" \", \n");
	mfprintf(a,"     RowBox[{\"The\", \" \", \"font\", \" \", \"size\", \" \", \"of\", \" \", \"these\", \" \", \n");
	mfprintf(a,"      RowBox[{\"labels\", \".\"}]}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"contourlabelsize\", \"=\", \"10\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"thickness\", \" \", \"of\", \" \", \"the\", \" \", \"contour\", \" \", \n");
	mfprintf(a,"    \"lines\", \" \", \n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\"(\", \n");
	mfprintf(a,"      RowBox[{\n");
	mfprintf(a,"      \"relative\", \" \", \"to\", \" \", \"the\", \" \", \"overall\", \" \", \"plot\", \" \", \n");
	mfprintf(a,"       \"width\"}], \")\"}], \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"contourthickness\", \"=\", \"0.001\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"The\", \" \", \"color\", \" \", \"of\", \" \", \"the\", \" \", \"contour\", \" \", \n");
	mfprintf(a,"     RowBox[{\"lines\", \".\", \" \", \"Enter\"}], \" \", \"three\", \" \", \"values\", \" \", \n");
	mfprintf(a,"     \"between\", \" \", \"0\", \" \", \"and\", \" \", \"1\", \" \", \"for\", \" \", \"red\"}], \",\",\n");
	mfprintf(a,"     \" \", \n");
	mfprintf(a,"    RowBox[{\"green\", \" \", \"and\", \" \", \"blue\", \" \", \"color\", \" \", \n");
	mfprintf(a,"     RowBox[{\"component\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"contourcolor\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"The\", \" \", \"opacity\", \" \", \"of\", \" \", \"the\", \" \", \"contour\", \" \", \n");
	mfprintf(a,"     RowBox[{\"lines\", \".\", \" \", \"1\"}], \" \", \"means\", \" \", \"fully\", \" \", \n");
	mfprintf(a,"     \"opaque\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"0\", \" \", \"means\", \" \", \"totally\", \" \", \"transparent\", \" \", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"(\", \"invisible\", \")\"}], \".\"}]}]}], \" \", \"*)\"}], \n");
	mfprintf(a,"  \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"contouropacity\", \"=\", \"0.7\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"The\", \" \", \"pattern\", \" \", \"of\", \" \", \"the\", \" \", \"contour\", \" \", \n");
	mfprintf(a,"    RowBox[{\"lines\", \".\", \" \", \"Large\"}], \" \", \"values\", \" \", \"like\", \" \", \n");
	mfprintf(a,"    \"1000\", \" \", \"mean\", \" \", \"no\", \" \", \n");
	mfprintf(a,"    RowBox[{\"dashing\", \".\", \" \", \"e\", \".\", \"g\", \".\", \" \", \"0.01\"}], \" \", \n");
	mfprintf(a,"    \"results\", \" \", \"in\", \" \", \"dashed\", \" \", \n");
	mfprintf(a,"    RowBox[{\"lines\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"contourdash\", \"=\", \"1000\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"3D Plot Properties\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Draw\", \" \", \"Tick\", \" \", \"labels\", \" \", \"on\", \" \", \"the\", \" \", \"Z\", \" \", \n");
	mfprintf(a,"    \"axis\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"labelsz\", \"=\", \".\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Draw\", \" \", \"a\", \" \", \"rectangular\", \" \", \"box\", \" \", \"around\", \" \", \n");
	mfprintf(a,"    \"the\", \" \", \"3\", \"D\", \" \", \"plot\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"boxed\", \"=\", \".\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"line\", \" \", \"width\", \" \", \"of\", \" \", \"this\", \" \", \"box\"}], \n");
	mfprintf(a,"   \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"boxthickness\", \"=\", \"0.002\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"The\", \" \", \"color\", \" \", \"of\", \" \", \"this\", \" \", \"box\"}], \" \", \n");
	mfprintf(a,"   \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"boxcolor\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"The\", \" \", \"opacity\", \" \", \"of\", \" \", \"this\", \" \", \"box\"}], \" \", \n");
	mfprintf(a,"   \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"boxopacity\", \"=\", \"1\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"May\", \" \", \"be\", \" \", \"either\", \" \", \n");
	mfprintf(a,"     RowBox[{\"{\", \n");
	mfprintf(a,"      RowBox[{\"y\", \",\", \"z\"}], \"}\"}], \" \", \"with\", \" \", \"y\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"z\", \" \", \"=\", \" \", \n");
	mfprintf(a,"     RowBox[{\"1\", \" \", \"/\", \" \", \n");
	mfprintf(a,"      RowBox[{\"-\", \"1\"}]}]}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"or\", \" \", \"Automatic\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"or\", \" \", \"None\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"xaxisposition\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"{\", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"-\", \"1\"}], \",\", \n");
	mfprintf(a,"      RowBox[{\"-\", \"1\"}]}], \"}\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"May\", \" \", \"be\", \" \", \"either\", \" \", \n");
	mfprintf(a,"     RowBox[{\"{\", \n");
	mfprintf(a,"      RowBox[{\"x\", \",\", \"z\"}], \"}\"}], \" \", \"with\", \" \", \"x\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"z\", \" \", \"=\", \" \", \n");
	mfprintf(a,"     RowBox[{\"1\", \" \", \"/\", \" \", \n");
	mfprintf(a,"      RowBox[{\"-\", \"1\"}]}]}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"or\", \" \", \"Automatic\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"or\", \" \", \"None\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"yaxisposition\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"{\", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"-\", \"1\"}], \",\", \n");
	mfprintf(a,"      RowBox[{\"-\", \"1\"}]}], \"}\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");

	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"May\", \" \", \"be\", \" \", \"either\", \" \", \n");
	mfprintf(a,"     RowBox[{\"{\", \n");
	mfprintf(a,"      RowBox[{\"x\", \",\", \"y\"}], \"}\"}], \" \", \"with\", \" \", \"x\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"y\", \" \", \"=\", \" \", \n");
	mfprintf(a,"     RowBox[{\"1\", \" \", \"/\", \" \", \n");
	mfprintf(a,"      RowBox[{\"-\", \"1\"}]}]}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"or\", \" \", \"Automatic\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"or\", \" \", \"None\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"zaxisposition\", \"=\", \"None\"}], \";\"}]}]], \"Input\"],\n");

	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"thickness\", \" \", \"of\", \" \", \"the\", \" \", \"axes\", \" \", \"and\", \n");
	mfprintf(a,"    \" \", \"tick\", \" \", \"marks\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"axesthickness\", \"=\", \"0.005\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"label\", \" \", \"to\", \" \", \"plot\", \" \", \"on\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    \"Z\", \" \", \"axis\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");

	if (m_sLabelZ == NULL)
		mfprintf(a,"   RowBox[{\"zlabel\", \"=\", \"\\\"\\<Occurrence\\>\\\"\"}], \";\"}]}]], \"Input\"],\n");
	else
		mfprintf(a,"   RowBox[{\"zlabel\", \"=\", \"\\\"\\<%s\\>\\\"\"}], \";\"}]}]], \"Input\"],\n",m_sLabelZ);

	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"How\", \" \", \"many\", \" \", \"grid\", \" \", \"lines\", \" \", \"to\", \" \", \"draw\", \" \",\n");
	mfprintf(a,"     \"in\", \" \", \"X\", \" \", \"direction\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"gridxcount\", \"=\", \"25\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Color\", \" \", \"for\", \" \", \"the\", \" \", \"grid\", \" \", \"lines\", \" \", \"in\", \" \",\n");
	mfprintf(a,"     \"X\", \" \", \"direction\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"gridxcolor\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Thickness\", \" \", \"of\", \" \", \"the\", \" \", \"grid\", \" \", \"lines\", \" \", \"in\", \n");
	mfprintf(a,"    \" \", \"X\", \" \", \"direction\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"gridxthickness\", \"=\", \"0.001\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Opacity\", \" \", \"of\", \" \", \"the\", \" \", \"grid\", \" \", \"lines\", \" \", \"in\", \n");
	mfprintf(a,"    \" \", \"X\", \" \", \"direction\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"gridxopacity\", \"=\", \"1\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"How\", \" \", \"many\", \" \", \"grid\", \" \", \"lines\", \" \", \"to\", \" \", \"draw\", \" \",\n");
	mfprintf(a,"     \"in\", \" \", \"Y\", \" \", \"direction\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"gridycount\", \"=\", \"25\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Color\", \" \", \"for\", \" \", \"the\", \" \", \"grid\", \" \", \"lines\", \" \", \"in\", \" \",\n");
	mfprintf(a,"     \"Y\", \" \", \"direction\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"gridycolor\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Thickness\", \" \", \"of\", \" \", \"the\", \" \", \"grid\", \" \", \"lines\", \" \", \"in\", \n");
	mfprintf(a,"    \" \", \"Y\", \" \", \"direction\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"gridythickness\", \"=\", \"0.001\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Opacity\", \" \", \"of\", \" \", \"the\", \" \", \"grid\", \" \", \"lines\", \" \", \"in\", \n");
	mfprintf(a,"    \" \", \"Y\", \" \", \"direction\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"gridyopacity\", \"=\", \"1\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Plot Range\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");

	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\"(*\", \" \",\n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Here\", \" \", \"you\", \" \", \"can\", \" \", \"\\\"\\<zoom in\\>\\\"\", \" \", \"into\", \" \", \n");
	mfprintf(a,"     \"a\", \" \", \"specified\", \" \", \"area\", \" \", \"of\", \" \", \"your\", \" \", \n");
	mfprintf(a,"     RowBox[{\"plot\", \".\", \" \", \"The\"}], \" \", \"range\", \" \", \"of\", \" \", \"your\", \n");
	mfprintf(a,"     \" \", \"data\", \" \", \"was\", \" \", \"originally\", \" \", \"X1\"}], \"=\", \"%f\"}], \n",m_fMinVal[0]);
	mfprintf(a,"   \",\", \" \", \n");
	mfprintf(a,"   RowBox[{\"Y1\", \"=\", \"%f\"}], \",\", \" \", \n",m_fMinVal[1]);
	mfprintf(a,"   RowBox[{\"X2\", \"=\", \"%f\"}], \",\", \" \", \n",m_fMaxVal[0]);
	mfprintf(a,"   RowBox[{\"Y2\", \"=\", \"%f.\"}]}], \" \", \"*)\"}]], \"Input\"],\n",m_fMaxVal[1]);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"plotrangeX1\", \"=\", \"%f\"}], \";\"}]], \"Input\"],\n",m_fMinVal[0]);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"plotrangeY1\", \"=\", \"%f\"}], \";\"}]], \"Input\"],\n",m_fMinVal[1]);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"plotrangeX2\", \"=\", \"%f\"}], \";\"}]], \"Input\"],\n",m_fMaxVal[0]);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"plotrangeY2\", \"=\", \"%f\"}], \";\"}]], \"Input\"]\n",m_fMaxVal[1]);
	mfprintf(a,"}, Closed]],\n");

	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Manual Color/Contour Scale\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");

	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"     \"Enter\", \" \", \"a\", \" \", \"value\", \" \", \"here\", \" \", \"ONLY\", \" \", \"if\", \n");
	mfprintf(a,"      \" \", \"you\", \" \", \"want\", \" \", \"to\", \" \", \"create\", \" \", \"different\", \n");
	mfprintf(a,"      \" \", \"plots\", \" \", \"with\", \" \", \"the\", \" \", \"same\", \" \", \"color\", \" \", \n");
	mfprintf(a,"      RowBox[{\"scale\", \".\", \"\\[IndentingNewLine]\", \"Then\"}], \" \", \"no\", \" \", \n");
	mfprintf(a,"      \"rescaling\", \" \", \"for\", \" \", \"smoothened\", \" \", \"data\", \" \", \"will\", \n");
	mfprintf(a,"      \" \", \n");
	mfprintf(a,"      RowBox[{\"occur\", \".\", \" \", \"Data\"}], \" \", \"for\", \" \", \"this\", \" \", \n");
	mfprintf(a,"      RowBox[{\"plot\", \":\", \" \", \"maxcolorvalue\"}]}], \" \", \"=\", \" \", \n");
	mfprintf(a,"     \"%f\"}], \";\", \" \", \n",m_fMaxEntry);
	mfprintf(a,"    RowBox[{\"mincolorvalue\", \" \", \"=\", \" \", \n");
	mfprintf(a,"     RowBox[{\"%c\", \"%f\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n",(m_fMinEntry<0)?'-':'+',m_fMinEntry);
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\n");

	if (manrange)
		mfprintf(a,"    RowBox[{\"mincolorvalue\", \"=\", \"%f\"}], \";\"}], \"\\[IndentingNewLine]\", \n",0.0);
	else
		mfprintf(a,"    RowBox[{\"mincolorvalue\", \"=.\"}], \";\"}], \"\\[IndentingNewLine]\", \n");

	mfprintf(a,"   RowBox[{\n");

	if (manrange)
		mfprintf(a,"    RowBox[{\"maxcolorvalue\", \"=\", \"%f\"}], \";\"}]}]}]], \"Input\"],\n",m_fMaxEntry);
	else
		mfprintf(a,"    RowBox[{\"maxcolorvalue\", \"=.\"}], \";\"}]}]}]], \"Input\"]\n");

	mfprintf(a,"\n");

	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Frame Thickness / Font Sizes\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Choose\", \" \", \"the\", \" \", \"font\", \" \", \"size\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"tick\", \" \", \"labels\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticklabelsize\", \"=\", \"28\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Choose\", \" \", \"the\", \" \", \"font\", \" \", \"size\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"axis\", \" \", \"labels\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"labelsize\", \"=\", \"32\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Choose\", \" \", \"the\", \" \", \"thickness\", \" \", \"of\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    \"frame\", \" \", \"and\", \" \", \"the\", \" \", \"axis\", \" \", \"ticks\"}], \" \", \"*)\"}],\n");
	mfprintf(a,"   \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"framethickness\", \"=\", \"0.005\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Define\", \" \", \"if\", \" \", \"space\", \" \", \"should\", \" \", \"be\", \" \", \"left\", \n");
	mfprintf(a,"    \" \", \"between\", \" \", \"the\", \" \", \"plot\", \" \", \"and\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axes\", \".\", \" \", \"Set\"}], \" \", \"this\", \" \", \"to\", \" \", \n");
	mfprintf(a,"    \"\\\"\\<Automatic\\>\\\"\", \" \", \"or\", \" \", \"to\", \" \", \n");
	mfprintf(a,"    RowBox[{\"\\\"\\<0\\>\\\"\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"framemargins\", \"=\", \"0\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Axis Properties\", \"Section\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"General Properties\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"length\", \" \", \"of\", \" \", \"the\", \" \", \"major\", \" \", \"ticks\", \n");
	mfprintf(a,"    \" \", \"on\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axes\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"majorticklength\", \"=\", \"0.03\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"length\", \" \", \"of\", \" \", \"the\", \" \", \"minor\", \" \", \"ticks\", \n");
	mfprintf(a,"    \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"the\", \" \", \"ones\", \" \", \"without\", \" \", \"tick\", \" \", \"labels\"}],\n");
	mfprintf(a,"      \")\"}], \" \", \"on\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axes\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minorticklength\", \"=\", \"0.015\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Controls\", \" \", \"if\", \" \", \"the\", \" \", \"major\", \" \", \"ticks\", \" \", \n");
	mfprintf(a,"     \"point\", \" \", \"to\", \" \", \"the\", \" \", \"inside\", \" \", \"or\", \" \", \n");
	mfprintf(a,"     RowBox[{\"outside\", \".\", \" \", \"0\"}], \" \", \"is\", \" \", \"inside\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"0.5\", \" \", \"is\", \" \", \"symmetrical\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"1\", \" \", \"is\", \" \", \n");
	mfprintf(a,"     RowBox[{\"outside\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"majortickshift\", \"=\", \"0.5\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Controls\", \" \", \"if\", \" \", \"the\", \" \", \"minor\", \" \", \"ticks\", \" \", \n");
	mfprintf(a,"     \"point\", \" \", \"to\", \" \", \"the\", \" \", \"inside\", \" \", \"or\", \" \", \n");
	mfprintf(a,"     RowBox[{\"outside\", \".\", \" \", \"0\"}], \" \", \"is\", \" \", \"inside\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"0.5\", \" \", \"is\", \" \", \"symmetrical\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"1\", \" \", \"is\", \" \", \n");
	mfprintf(a,"     RowBox[{\"outside\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minortickshift\", \"=\", \"0.5\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");

	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Mesh Properties\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]],\n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"thickness\", \" \", \"of\", \" \", \"the\", \" \", \"mesh\", \" \", \n");
	mfprintf(a,"    \"perpendicular\", \" \", \"to\", \" \", \"the\", \" \", \"X\", \" \", \"axis\", \" \", \n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\"(\", \n");
	mfprintf(a,"      RowBox[{\n");
	mfprintf(a,"      \"relative\", \" \", \"to\", \" \", \"the\", \" \", \"overall\", \" \", \"plot\", \" \", \n");
	mfprintf(a,"       \"width\"}], \")\"}], \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"meshthicknessX\", \"=\", \"0.002\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \",\n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"The\", \" \", \"color\", \" \", \"of\", \" \", \"the\", \" \", \"mesh\", \" \", \n");
	mfprintf(a,"     \"perpendicular\", \" \", \"to\", \" \", \"the\", \" \", \"X\", \" \", \n");
	mfprintf(a,"     RowBox[{\"axis\", \".\", \" \", \"Enter\"}], \" \", \"three\", \" \", \"values\", \" \", \n");
	mfprintf(a,"     \"between\", \" \", \"0\", \" \", \"and\", \" \", \"1\", \" \", \"for\", \" \", \"red\"}], \",\",\n");
	mfprintf(a,"     \" \", \n");
	mfprintf(a,"    RowBox[{\"green\", \" \", \"and\", \" \", \"blue\", \" \", \"color\", \" \", \n");
	mfprintf(a,"     RowBox[{\"component\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"meshcolorX\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"The\", \" \", \"opacity\", \" \", \"of\", \" \", \"the\", \" \", \"mesh\", \" \", \n");
	mfprintf(a,"     \"perpendicular\", \" \", \"to\", \" \", \"the\", \" \", \"X\", \" \", \n");
	mfprintf(a,"     RowBox[{\"axis\", \".\", \" \", \"1\"}], \" \", \"means\", \" \", \"fully\", \" \", \n");
	mfprintf(a,"     \"opaque\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"0\", \" \", \"means\", \" \", \"totally\", \" \", \"transparent\", \" \", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"(\", \"invisible\", \")\"}], \".\"}]}]}], \" \", \"*)\"}], \n");
	mfprintf(a,"  \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"meshopacityX\", \"=\", \"0.5\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"thickness\", \" \", \"of\", \" \", \"the\", \" \", \"mesh\", \" \", \n");
	mfprintf(a,"    \"perpendicular\", \" \", \"to\", \" \", \"the\", \" \", \"Y\", \" \", \"axis\", \" \", \n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"     RowBox[{\"(\", \n");
	mfprintf(a,"      RowBox[{\n");
	mfprintf(a,"      \"relative\", \" \", \"to\", \" \", \"the\", \" \", \"overall\", \" \", \"plot\", \" \", \n");
	mfprintf(a,"       \"width\"}], \")\"}], \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"meshthicknessY\", \"=\", \"0.002\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"The\", \" \", \"color\", \" \", \"of\", \" \", \"the\", \" \", \"mesh\", \" \", \n");
	mfprintf(a,"     \"perpendicular\", \" \", \"to\", \" \", \"the\", \" \", \"Y\", \" \", \n");
	mfprintf(a,"     RowBox[{\"axis\", \".\", \" \", \"Enter\"}], \" \", \"three\", \" \", \"values\", \" \", \n");
	mfprintf(a,"     \"between\", \" \", \"0\", \" \", \"and\", \" \", \"1\", \" \", \"for\", \" \", \"red\"}], \",\",\n");
	mfprintf(a,"     \" \", \n");
	mfprintf(a,"    RowBox[{\"green\", \" \", \"and\", \" \", \"blue\", \" \", \"color\", \" \", \n");
	mfprintf(a,"     RowBox[{\"component\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"meshcolorY\", \"=\", \n");
	mfprintf(a,"    RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"The\", \" \", \"opacity\", \" \", \"of\", \" \", \"the\", \" \", \"mesh\", \" \", \n");
	mfprintf(a,"     \"perpendicular\", \" \", \"to\", \" \", \"the\", \" \", \"Y\", \" \", \n");
	mfprintf(a,"     RowBox[{\"axis\", \".\", \" \", \"1\"}], \" \", \"means\", \" \", \"fully\", \" \", \n");
	mfprintf(a,"     \"opaque\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"0\", \" \", \"means\", \" \", \"totally\", \" \", \"transparent\", \" \", \n");
	mfprintf(a,"     RowBox[{\n");
	mfprintf(a,"      RowBox[{\"(\", \"invisible\", \")\"}], \".\"}]}]}], \" \", \"*)\"}], \n");
	mfprintf(a,"  \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"meshopacityY\", \"=\", \"0.5\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed  ]],\n");

	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"X Axis\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"label\", \" \", \"to\", \" \", \"appear\", \" \", \"on\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    \"x\", \" \", \"axis\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"xlabel\", \"=\", \"\\\"\\<%s\\>\\\"\"}], \n",m_sLabelX);
	mfprintf(a,"   \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"lowest\", \" \", \"tick\", \" \", \"value\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"x\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickminX\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"maxsig\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"plotrangeX1\", \",\", \"2\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"highest\", \" \", \"tick\", \" \", \"value\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"x\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickmaxX\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"minsig\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"plotrangeX2\", \",\", \"2\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Print\", \" \", \"this\", \" \", \"number\", \" \", \"of\", \" \", \"digits\", \" \", \n");
	mfprintf(a,"    \"after\", \" \", \"the\", \" \", \"decimal\", \" \", \"point\", \" \", \"of\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"labels\", \".\", \" \", \"0\"}], \" \", \"means\", \" \", \"no\", \" \", \n");
	mfprintf(a,"    \"decimal\", \" \", \n");
	mfprintf(a,"    RowBox[{\"point\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickprecX\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"prec\", \"[\", \"plotrangeX2\", \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Use\", \" \", \"this\", \" \", \"count\", \" \", \"of\", \" \", \"major\", \" \", \"ticks\", \n");
	mfprintf(a,"    \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"with\", \" \", \"tick\", \" \", \"label\"}], \")\"}], \" \", \"on\", \" \", \n");
	mfprintf(a,"    \"the\", \" \", \"x\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"majorticksX\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",majorx);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"How\", \" \", \"many\", \" \", \"minor\", \" \", \"ticks\", \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"without\", \" \", \"tick\", \" \", \"label\"}], \")\"}], \" \", \"to\", \" \", \n");
	mfprintf(a,"    \"show\", \" \", \"PER\", \" \", \"MAJOR\", \" \", \n");
	mfprintf(a,"    RowBox[{\"TICK\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minorticksX\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",minorx);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"A\", \" \", \"character\", \" \", \"string\", \" \", \"to\", \" \", \"stand\", \" \", \n");
	mfprintf(a,"    \"before\", \" \", \"each\", \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"label\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticklabelprefixX\", \"=\", \"\\\"\\<\\>\\\"\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"A\", \" \", \"character\", \" \", \"string\", \" \", \"to\", \" \", \"stand\", \" \", \n");
	mfprintf(a,"    \"after\", \" \", \"each\", \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"label\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticklabelsuffixX\", \"=\", \"\\\"\\<\\>\\\"\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Show\", \" \", \"tick\", \" \", \"marks\", \" \", \"on\", \" \", \"both\", \" \", \"the\", \" \",\n");
	mfprintf(a,"     \"top\", \" \", \"X\", \" \", \"axis\", \" \", \"and\", \" \", \"the\", \" \", \"bottom\", \" \",\n");
	mfprintf(a,"     \"X\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \"?\", \" \", \"Set\"}], \" \", \"to\", \" \", \"\\\"\\<1\\>\\\"\", \" \", \"if\",\n");
	mfprintf(a,"     \" \", \"desired\", \" \", \"or\", \" \", \"to\", \" \", \"\\\"\\<.\\>\\\"\", \" \", \"if\", \" \", \n");
	mfprintf(a,"    RowBox[{\"not\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticksbothsidesX\", \"=\", \"1\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Show\", \" \", \"tick\", \" \", \"labels\", \" \", \"on\", \" \", \"both\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"top\", \" \", \"X\", \" \", \"axis\", \" \", \"and\", \" \", \"the\", \" \", \"bottom\", \n");
	mfprintf(a,"    \" \", \"X\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \"?\", \" \", \"Set\"}], \" \", \"to\", \" \", \"\\\"\\<1\\>\\\"\", \" \", \"if\",\n");
	mfprintf(a,"     \" \", \"desired\", \" \", \"or\", \" \", \"to\", \" \", \"\\\"\\<.\\>\\\"\", \" \", \"if\", \" \", \n");
	mfprintf(a,"    RowBox[{\"not\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"labelsbothsidesX\", \"=.\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Y Axis\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"label\", \" \", \"to\", \" \", \"appear\", \" \", \"on\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    \"Y\", \" \", \"axis\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ylabel\", \"=\", \"\\\"\\<%s\\>\\\"\"}], \n",m_sLabelY);
	mfprintf(a,"   \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"lowest\", \" \", \"tick\", \" \", \"value\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"Y\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickminY\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"maxsig\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"plotrangeY1\", \",\", \"2\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"highest\", \" \", \"tick\", \" \", \"value\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"Y\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickmaxY\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"minsig\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"plotrangeY2\", \",\", \"2\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Print\", \" \", \"this\", \" \", \"number\", \" \", \"of\", \" \", \"digits\", \" \", \n");
	mfprintf(a,"    \"after\", \" \", \"the\", \" \", \"decimal\", \" \", \"point\", \" \", \"of\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"labels\", \".\", \" \", \"0\"}], \" \", \"means\", \" \", \"no\", \" \", \n");
	mfprintf(a,"    \"decimal\", \" \", \n");
	mfprintf(a,"    RowBox[{\"point\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickprecY\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"prec\", \"[\", \"plotrangeY2\", \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Use\", \" \", \"this\", \" \", \"count\", \" \", \"of\", \" \", \"major\", \" \", \"ticks\", \n");
	mfprintf(a,"    \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"with\", \" \", \"tick\", \" \", \"label\"}], \")\"}], \" \", \"on\", \" \", \n");
	mfprintf(a,"    \"the\", \" \", \"Y\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"majorticksY\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",majory);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"How\", \" \", \"many\", \" \", \"minor\", \" \", \"ticks\", \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"without\", \" \", \"tick\", \" \", \"label\"}], \")\"}], \" \", \"to\", \" \", \n");
	mfprintf(a,"    \"show\", \" \", \"PER\", \" \", \"MAJOR\", \" \", \n");
	mfprintf(a,"    RowBox[{\"TICK\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minorticksY\", \"=\", \"%d\"}], \";\"}]}]], \"Input\"],\n",minory);
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"A\", \" \", \"character\", \" \", \"string\", \" \", \"to\", \" \", \"stand\", \" \", \n");
	mfprintf(a,"    \"before\", \" \", \"each\", \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"label\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticklabelprefixY\", \"=\", \"\\\"\\<\\>\\\"\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"A\", \" \", \"character\", \" \", \"string\", \" \", \"to\", \" \", \"stand\", \" \", \n");
	mfprintf(a,"    \"after\", \" \", \"each\", \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"label\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticklabelsuffixY\", \"=\", \"\\\"\\<\\>\\\"\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Show\", \" \", \"tick\", \" \", \"marks\", \" \", \"on\", \" \", \"both\", \" \", \"the\", \" \",\n");
	mfprintf(a,"     \"left\", \" \", \"Y\", \" \", \"axis\", \" \", \"and\", \" \", \"the\", \" \", \"right\", \" \",\n");
	mfprintf(a,"     \"Y\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \"?\", \" \", \"Set\"}], \" \", \"to\", \" \", \"\\\"\\<1\\>\\\"\", \" \", \"if\",\n");
	mfprintf(a,"     \" \", \"desired\", \" \", \"or\", \" \", \"to\", \" \", \"\\\"\\<.\\>\\\"\", \" \", \"if\", \" \", \n");
	mfprintf(a,"    RowBox[{\"not\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticksbothsidesY\", \"=\", \"1\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Show\", \" \", \"tick\", \" \", \"labels\", \" \", \"on\", \" \", \"both\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"left\", \" \", \"Y\", \" \", \"axis\", \" \", \"and\", \" \", \"the\", \" \", \"right\", \n");
	mfprintf(a,"    \" \", \"Y\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \"?\", \" \", \"Set\"}], \" \", \"to\", \" \", \"\\\"\\<1\\>\\\"\", \" \", \"if\",\n");
	mfprintf(a,"     \" \", \"desired\", \" \", \"or\", \" \", \"to\", \" \", \"\\\"\\<.\\>\\\"\", \" \", \"if\", \" \", \n");
	mfprintf(a,"    RowBox[{\"not\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"labelsbothsidesY\", \"=.\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Z Axis\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False},\n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"label\", \" \", \"to\", \" \", \"appear\", \" \", \"on\", \" \", \"the\", \" \", \n");
	mfprintf(a,"    \"Z\", \" \", \"axis\"}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");

	if (m_sLabelZ == NULL)
		mfprintf(a,"   RowBox[{\"zlabel\", \"=\", \"\\\"\\<Occurrence\\>\\\"\"}], \n");
	else
		mfprintf(a,"   RowBox[{\"zlabel\", \"=\", \"\\\"\\<%s\\>\\\"\"}], \n",m_sLabelZ);

	mfprintf(a,"   \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"lowest\", \" \", \"tick\", \" \", \"value\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"Z\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickminZ\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"maxsig\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"datminexp\", \",\", \"2\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"highest\", \" \", \"tick\", \" \", \"value\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"Z\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickmaxZ\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"minsig\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"datmaxexp\", \",\", \"2\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Print\", \" \", \"this\", \" \", \"number\", \" \", \"of\", \" \", \"digits\", \" \", \n");
	mfprintf(a,"    \"after\", \" \", \"the\", \" \", \"decimal\", \" \", \"point\", \" \", \"of\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"labels\", \".\", \" \", \"0\"}], \" \", \"means\", \" \", \"no\", \" \", \n");
	mfprintf(a,"    \"decimal\", \" \", \n");
	mfprintf(a,"    RowBox[{\"point\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickprecZ\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"prec\", \"[\", \"datmaxexp\", \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Use\", \" \", \"this\", \" \", \"count\", \" \", \"of\", \" \", \"major\", \" \", \"ticks\", \n");
	mfprintf(a,"    \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"with\", \" \", \"tick\", \" \", \"label\"}], \")\"}], \" \", \"on\", \" \", \n");
	mfprintf(a,"    \"the\", \" \", \"Z\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"majorticksZ\", \"=\", \"5\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"How\", \" \", \"many\", \" \", \"minor\", \" \", \"ticks\", \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"without\", \" \", \"tick\", \" \", \"label\"}], \")\"}], \" \", \"to\", \" \", \n");
	mfprintf(a,"    \"show\", \" \", \"PER\", \" \", \"MAJOR\", \" \", \n");
	mfprintf(a,"    RowBox[{\"TICK\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minorticksZ\", \"=\", \"3\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"A\", \" \", \"character\", \" \", \"string\", \" \", \"to\", \" \", \"stand\", \" \", \n");
	mfprintf(a,"    \"before\", \" \", \"each\", \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"label\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticklabelprefixZ\", \"=\", \"\\\"\\<\\>\\\"\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"A\", \" \", \"character\", \" \", \"string\", \" \", \"to\", \" \", \"stand\", \" \", \n");
	mfprintf(a,"    \"after\", \" \", \"each\", \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"label\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"ticklabelsuffixZ\", \"=\", \"\\\"\\<\\>\\\"\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Legend Axis\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"lowest\", \" \", \"tick\", \" \", \"value\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"legend\", \" \", \"X\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\", \" \", \"\\\"\\<datmin\\>\\\"\"}], \" \", \"sets\", \" \", \"this\", \n");
	mfprintf(a,"    \" \", \"automatically\", \" \", \"to\", \" \", \"the\", \" \", \"lowest\", \" \", \"plot\", \n");
	mfprintf(a,"    \" \", \n");
	mfprintf(a,"    RowBox[{\"element\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickminL\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"maxsig\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"datmin\", \",\", \"3\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"highest\", \" \", \"tick\", \" \", \"value\", \" \", \"for\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"legend\", \" \", \"Y\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\", \" \", \"\\\"\\<datmax\\>\\\"\"}], \" \", \"sets\", \" \", \"this\", \n");
	mfprintf(a,"    \" \", \"automatically\", \" \", \"to\", \" \", \"the\", \" \", \"highest\", \" \", \"plot\", \n");
	mfprintf(a,"    \" \", \n");
	mfprintf(a,"    RowBox[{\"element\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickmaxL\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"minsig\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"datmax\", \",\", \"3\"}], \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Print\", \" \", \"this\", \" \", \"number\", \" \", \"of\", \" \", \"digits\", \" \", \n");
	mfprintf(a,"    \"after\", \" \", \"the\", \" \", \"decimal\", \" \", \"point\", \" \", \"of\", \" \", \"the\", \n");
	mfprintf(a,"    \" \", \"tick\", \" \", \n");
	mfprintf(a,"    RowBox[{\"labels\", \".\", \" \", \"0\"}], \" \", \"means\", \" \", \"no\", \" \", \n");
	mfprintf(a,"    \"decimal\", \" \", \n");
	mfprintf(a,"    RowBox[{\"point\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"tickprecL\", \":=\", \n");
	mfprintf(a,"    RowBox[{\"prec\", \"[\", \"datmax\", \"]\"}]}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"Use\", \" \", \"this\", \" \", \"count\", \" \", \"of\", \" \", \"major\", \" \", \"ticks\", \n");
	mfprintf(a,"    \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"with\", \" \", \"tick\", \" \", \"label\"}], \")\"}], \" \", \"on\", \" \", \n");
	mfprintf(a,"    \"the\", \" \", \"legend\", \" \", \"X\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"majorticksL\", \"=\", \"5\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"How\", \" \", \"many\", \" \", \"minor\", \" \", \"ticks\", \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"without\", \" \", \"tick\", \" \", \"label\"}], \")\"}], \" \", \"to\", \" \", \n");
	mfprintf(a,"    \"show\", \" \", \"PER\", \" \", \"MAJOR\", \" \", \n");
	mfprintf(a,"    RowBox[{\"TICK\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minorticksL\", \"=\", \"3\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"length\", \" \", \"of\", \" \", \"the\", \" \", \"major\", \" \", \"ticks\", \n");
	mfprintf(a,"    \" \", \"on\", \" \", \"the\", \" \", \"legend\", \" \", \"X\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"majorticklengthL\", \"=\", \"0.015\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"   \"The\", \" \", \"length\", \" \", \"of\", \" \", \"the\", \" \", \"minor\", \" \", \"ticks\", \n");
	mfprintf(a,"    \" \", \n");
	mfprintf(a,"    RowBox[{\"(\", \n");
	mfprintf(a,"     RowBox[{\"the\", \" \", \"ones\", \" \", \"without\", \" \", \"tick\", \" \", \"labels\"}],\n");
	mfprintf(a,"      \")\"}], \" \", \"on\", \" \", \"the\", \" \", \"legend\", \" \", \"X\", \" \", \n");
	mfprintf(a,"    RowBox[{\"axis\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minorticklengthL\", \"=\", \"0.008\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Controls\", \" \", \"if\", \" \", \"the\", \" \", \"major\", \" \", \"ticks\", \" \", \n");
	mfprintf(a,"     \"point\", \" \", \"to\", \" \", \"the\", \" \", \"inside\", \" \", \"or\", \" \", \n");
	mfprintf(a,"     RowBox[{\"outside\", \".\", \" \", \"0\"}], \" \", \"is\", \" \", \"inside\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"0.5\", \" \", \"is\", \" \", \"symmetrical\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"1\", \" \", \"is\", \" \", \n");
	mfprintf(a,"     RowBox[{\"outside\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"majortickshiftL\", \"=\", \"0\"}], \";\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Controls\", \" \", \"if\", \" \", \"the\", \" \", \"minor\", \" \", \"ticks\", \" \", \n");
	mfprintf(a,"     \"point\", \" \", \"to\", \" \", \"the\", \" \", \"inside\", \" \", \"or\", \" \", \n");
	mfprintf(a,"     RowBox[{\"outside\", \".\", \" \", \"0\"}], \" \", \"is\", \" \", \"inside\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"0.5\", \" \", \"is\", \" \", \"symmetrical\"}], \",\", \" \", \n");
	mfprintf(a,"    RowBox[{\"1\", \" \", \"is\", \" \", \n");
	mfprintf(a,"     RowBox[{\"outside\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\n");
	mfprintf(a,"   RowBox[{\"minortickshiftL\", \"=\", \"0\"}], \";\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Lines / Rectangles / Functions\", \"Section\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Function Drawing\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");


	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Here\", \" \", \"some\", \" \", \"functions\", \" \", \"of\", \" \", \"the\", \" \", \"form\",\n");
	mfprintf(a,"      \" \", \"f\", \n");
	mfprintf(a,"     RowBox[{\"(\", \"x\", \")\"}]}], \" \", \"=\", \" \", \n");
	mfprintf(a,"    RowBox[{\"y\", \" \", \"are\", \" \", \n");
	mfprintf(a,"     RowBox[{\"defined\", \".\", \" \", \"Uncomment\"}], \" \", \"them\", \" \", \"for\", \" \",\n");
	mfprintf(a,"      \"testing\", \" \", \"and\", \" \", \"define\", \" \", \"your\", \" \", \"own\", \" \", \n");
	mfprintf(a,"     RowBox[{\"ones\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\"(*\", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"funcX\", \",\", \n");
	mfprintf(a,"      RowBox[{\"{\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"x\", \"/\", \"%f\"}], \")\"}], \"^\", \"2\"}], \"*\", \n",m_fMaxVal[0]);
	mfprintf(a,"         \"%f\"}], \",\", \n",m_fMaxVal[1]);
	mfprintf(a,"        RowBox[{\"{\", \n");
	mfprintf(a,"         RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}], \",\", \n",m_fMinVal[0],m_fMaxVal[0]);
	mfprintf(a,"        RowBox[{\"{\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"Thickness\", \"[\", \"0.01\", \"]\"}], \",\", \n");
	mfprintf(a,"          RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"           RowBox[{\"1\", \",\", \"1\", \",\", \"0\"}], \"]\"}], \",\", \n");
	mfprintf(a,"          RowBox[{\"Dashing\", \"[\", \"0.02\", \"]\"}]}], \"}\"}]}], \"}\"}]}], \"]\"}], \n");
	mfprintf(a,"    \";\"}], \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\"(*\", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"funcX\", \",\", \n");
	mfprintf(a,"      RowBox[{\"{\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"(\", \n");
	mfprintf(a,"           RowBox[{\"x\", \"/\", \"%f\"}], \")\"}], \"^\", \"3\"}], \"*\", \n",m_fMaxVal[0]);
	mfprintf(a,"         \"%f\"}], \",\", \n",m_fMaxVal[1]);
	mfprintf(a,"        RowBox[{\"{\", \n");
	mfprintf(a,"         RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}], \",\", \n",m_fMinVal[0]+(m_fMaxVal[0]-m_fMinVal[0])*0.2,m_fMaxVal[0]-(m_fMaxVal[0]-m_fMinVal[0])*0.2);
	mfprintf(a,"        RowBox[{\"{\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"Thickness\", \"[\", \"0.02\", \"]\"}], \",\", \n");
	mfprintf(a,"          RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"           RowBox[{\"0\", \",\", \"1\", \",\", \"1\"}], \"]\"}]}], \"}\"}]}], \"}\"}]}], \n");
	mfprintf(a,"     \"]\"}], \";\"}], \"*)\"}]}]], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\n");
	mfprintf(a,"    \"Here\", \" \", \"some\", \" \", \"functions\", \" \", \"of\", \" \", \"the\", \" \", \"form\",\n");
	mfprintf(a,"      \" \", \"f\", \n");
	mfprintf(a,"     RowBox[{\"(\", \"y\", \")\"}]}], \" \", \"=\", \" \", \n");
	mfprintf(a,"    RowBox[{\"x\", \" \", \"are\", \" \", \n");
	mfprintf(a,"     RowBox[{\"defined\", \".\", \" \", \"Uncomment\"}], \" \", \"them\", \" \", \"for\", \" \",\n");
	mfprintf(a,"      \"testing\", \" \", \"and\", \" \", \"define\", \" \", \"your\", \" \", \"own\", \" \", \n");
	mfprintf(a,"     RowBox[{\"ones\", \".\"}]}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\"(*\", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"funcY\", \",\", \n");
	mfprintf(a,"      RowBox[{\"{\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"Sin\", \"[\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"y\", \"/\", \"%f\"}], \"*\", \"10\", \"Pi\"}], \"]\"}], \"*\", \n",m_fMaxVal[1]);
	mfprintf(a,"         \"%f\"}], \",\", \n",m_fMaxVal[0]*0.5);
	mfprintf(a,"        RowBox[{\"{\", \n");
	mfprintf(a,"         RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}], \",\", \n",m_fMinVal[1],m_fMaxVal[1]);
	mfprintf(a,"        RowBox[{\"{\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"Thickness\", \"[\", \"0.01\", \"]\"}], \",\", \n");
	mfprintf(a,"          RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"           RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \"}\"}]}], \"}\"}]}], \n");
	mfprintf(a,"     \"]\"}], \";\"}], \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\"(*\", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"funcY\", \",\", \n");
	mfprintf(a,"      RowBox[{\"{\", \n");
	mfprintf(a,"       RowBox[{\n");
	mfprintf(a,"        RowBox[{\"y\", \"*\", \"%f\"}], \",\", \n",m_fMaxVal[0]/m_fMaxVal[1]);
	mfprintf(a,"        RowBox[{\"{\", \n");
	mfprintf(a,"         RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}], \",\", \n",m_fMinVal[1],m_fMaxVal[1]*0.5);
	mfprintf(a,"        RowBox[{\"{\", \n");
	mfprintf(a,"         RowBox[{\n");
	mfprintf(a,"          RowBox[{\"Thickness\", \"[\", \"0.02\", \"]\"}], \",\", \n");
	mfprintf(a,"          RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"           RowBox[{\"1\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \"}\"}]}], \"}\"}]}], \n");
	mfprintf(a,"     \"]\"}], \";\"}], \"*)\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");

	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"Cell[\"Rectangle and Circle Drawing\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"Draw\", \" \", \"Rectangles\", \" \", \"and\", \" \", \n");
	mfprintf(a,"    RowBox[{\"Circles\", \".\", \" \", \"Uncomment\"}], \" \", \"the\", \" \", \"examples\", \n");
	mfprintf(a,"    \" \", \"below\", \" \", \"and\", \" \", \"try\", \" \", \"it\", \" \", \n");
	mfprintf(a,"    RowBox[{\"out\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\"(*\", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"prim\", \",\", \n");
	mfprintf(a,"      RowBox[{\"Graphics\", \"[\", \n");
	mfprintf(a,"       RowBox[{\"{\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"EdgeForm\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Thickness\", \"[\", \"0.006\", \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"             RowBox[{\"1\", \",\", \"0\", \",\", \"1\"}], \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"Dashing\", \"[\", \"0.02\", \"]\"}]}], \"}\"}], \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"1\", \",\", \"0\", \",\", \"0\"}], \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"Opacity\", \"[\", \"0\", \"]\"}], \",\", \n");


	mfprintf(a,"         RowBox[{\"Rectangle\", \"[\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"{\", \n");
	mfprintf(a,"            RowBox[{\n");
	mfprintf(a,"             RowBox[{\"%f\"}], \",\", \n",m_fMinVal[0]);
	mfprintf(a,"             RowBox[{\"%f\"}]}], \"}\"}], \",\", \n",m_fMinVal[1]);
	mfprintf(a,"           RowBox[{\"{\", \n");
	mfprintf(a,"            RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}]}], \"]\"}]}], \n",m_fMaxVal[0]*0.7,m_fMaxVal[1]*0.3);

	mfprintf(a,"        \"}\"}], \"]\"}]}], \"]\"}], \";\", \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"    RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"prim\", \",\", \n");
	mfprintf(a,"      RowBox[{\"Graphics\", \"[\", \n");
	mfprintf(a,"       RowBox[{\"{\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"EdgeForm\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"Thickness\", \"[\", \"0.004\", \"]\"}], \",\", \n");
	mfprintf(a,"            RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"             RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \"}\"}], \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"0\", \",\", \"1\", \",\", \"0\"}], \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"Opacity\", \"[\", \"1\", \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"Disk\", \"[\", \n");
	mfprintf(a,"          RowBox[{\n");
	mfprintf(a,"           RowBox[{\"{\", \n");
	mfprintf(a,"            RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}]}], \"]\"}]}], \"}\"}], \n",m_fMinVal[0]+(m_fMaxVal[0]-m_fMinVal[0])*0.41666,(m_fMaxVal[0]-m_fMinVal[0])*0.0625,(m_fMaxVal[0]-m_fMinVal[0])*0.0625,(m_fMaxVal[1]-m_fMinVal[1])*0.0625);
	mfprintf(a,"       \"]\"}]}], \"]\"}], \";\"}], \"*)\"}]\n");

	if (m_oaCircles.GetSize() > 0)
	{
		mfprintf(a,",\"\\[IndentingNewLine]\", \n");
		mfprintf(a,"RowBox[{\n");
		for (z=0;z<m_oaCircles.GetSize();z++)
		{
			mfprintf(a,"   RowBox[{\n");
			mfprintf(a,"    RowBox[{\"AppendTo\", \"[\", \n");
			mfprintf(a,"     RowBox[{\"prim\", \",\", \n");
			mfprintf(a,"      RowBox[{\"Graphics\", \"[\", \n");
			mfprintf(a,"       RowBox[{\"{\", \n");
			mfprintf(a,"        RowBox[{\n");
			mfprintf(a,"         RowBox[{\"EdgeForm\", \"[\", \n");
			mfprintf(a,"          RowBox[{\"{\", \n");
			mfprintf(a,"           RowBox[{\n");
			mfprintf(a,"            RowBox[{\"Thickness\", \"[\", \"0.002\", \"]\"}], \",\", \n");
			mfprintf(a,"            RowBox[{\"RGBColor\", \"[\", \n");
			mfprintf(a,"             RowBox[{\"0\", \",\", \"0\", \",\", \"0\"}], \"]\"}]}], \"}\"}], \"]\"}], \",\", \n");
			mfprintf(a,"         RowBox[{\"RGBColor\", \"[\", \n");
			mfprintf(a,"          RowBox[{\"%f\", \",\", \"%f\", \",\", \"%f\"}], \"]\"}], \",\", \n",((CMathematicaCircle*)m_oaCircles[z])->m_fColorR,((CMathematicaCircle*)m_oaCircles[z])->m_fColorG,((CMathematicaCircle*)m_oaCircles[z])->m_fColorB);
			mfprintf(a,"         RowBox[{\"Opacity\", \"[\", \"1\", \"]\"}], \",\", \n");
			mfprintf(a,"         RowBox[{\"Disk\", \"[\", \n");
			mfprintf(a,"          RowBox[{\n");
			mfprintf(a,"           RowBox[{\"{\", \n");
			mfprintf(a,"            RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}], \",\", \"%f\"}], \"]\"}]}], \"}\"}], \n",((CMathematicaCircle*)m_oaCircles[z])->m_fPosX,((CMathematicaCircle*)m_oaCircles[z])->m_fPosY,((CMathematicaCircle*)m_oaCircles[z])->m_fRadius);
			mfprintf(a,"       \"]\"}]}], \"]\"}], \";\"}]");
			if (z < m_oaCircles.GetSize()-1)
				mfprintf(a,", \"\\[IndentingNewLine]\", ");
			mfprintf(a,"\n");
		}
		mfprintf(a,"}]\n");
	}
	mfprintf(a,"}]], \"Input\"]\n");
	mfprintf(a,"}, Closed  ]],\n");

	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"Line Drawing\", \"Subsection\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = False},\n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\n");
	mfprintf(a," RowBox[{\n");
	mfprintf(a,"  RowBox[{\"(*\", \" \", \n");
	mfprintf(a,"   RowBox[{\"Draw\", \" \", \n");
	mfprintf(a,"    RowBox[{\"Lines\", \".\", \" \", \"Uncomment\"}], \" \", \"the\", \" \", \"example\", \" \",\n");
	mfprintf(a,"     \"below\", \" \", \"and\", \" \", \"try\", \" \", \"it\", \" \", \n");
	mfprintf(a,"    RowBox[{\"out\", \".\"}]}], \" \", \"*)\"}], \"\\[IndentingNewLine]\", \n");
	mfprintf(a,"  RowBox[{\"(*\", \n");
	mfprintf(a,"   RowBox[{\n");
	mfprintf(a,"    RowBox[{\"AppendTo\", \"[\", \n");
	mfprintf(a,"     RowBox[{\"prim\", \",\", \n");
	mfprintf(a,"      RowBox[{\"Graphics\", \"[\", \n");
	mfprintf(a,"       RowBox[{\"{\", \n");
	mfprintf(a,"        RowBox[{\n");
	mfprintf(a,"         RowBox[{\"Thickness\", \"[\", \"0.01\", \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"RGBColor\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"0\", \",\", \"0\", \",\", \"1\"}], \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"Dashing\", \"[\", \"1000\", \"]\"}], \",\", \n");
	mfprintf(a,"         RowBox[{\"Line\", \"[\", \n");
	mfprintf(a,"          RowBox[{\"{\", \n");
	mfprintf(a,"           RowBox[{\n");
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}], \",\", \n",m_fMinVal[0],m_fMinVal[1]);
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}], \",\", \n",m_fMaxVal[0]*0.4,m_fMaxVal[1]*0.8);
	mfprintf(a,"            RowBox[{\"{\", \n");
	mfprintf(a,"             RowBox[{\"%f\", \",\", \"%f\"}], \"}\"}]}], \"}\"}], \"]\"}]}],\n",m_fMaxVal[0]*0.4,m_fMinVal[1]);
	mfprintf(a,"         \"}\"}], \"]\"}]}], \"]\"}], \";\"}], \"*)\"}]}]], \"Input\"]\n");
	mfprintf(a,"}, Closed]]\n");
	mfprintf(a,"}, Closed]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[CellGroupData[{\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[\"\\<\\\n");
	mfprintf(a,"Plot output. Click \\\"Evaluation -> Evaluate Notebook\\\" to draw/redraw!\\\n");
	mfprintf(a,"\\>\", \"Section\",\n");
	mfprintf(a," CellDingbat->DynamicModuleBox[{$CellContext`state$$ = True}, \n");
	mfprintf(a,"   OpenerBox[\n");
	mfprintf(a,"    Dynamic[$CellContext`state$$, (FrontEndExecute[{\n");
	mfprintf(a,"        FrontEnd`SelectionMove[\n");
	mfprintf(a,"         FrontEnd`ButtonNotebook[], All, ButtonCell], \n");
	mfprintf(a,"        FrontEndToken[\"OpenCloseGroup\"]}]; $CellContext`state$$ = #)& ]], \n");
	mfprintf(a,"   DynamicModuleValues :> {}]],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\"plot\"], \"Input\"],\n");
	mfprintf(a,"\n");
	mfprintf(a,"Cell[BoxData[\"plotlegend\"], \"Input\"]\n");
	mfprintf(a,"}, Open  ]]\n");
	mfprintf(a,"}\n");
	mfprintf(a,"]\n");

	fclose(a);
	BTOUT;
}


void C2DF::NormRDF(double n)
{
	BTIN;
	int x, y;
	double f;

	for (y=0;y<m_iRes[1];y++)
	{
		f = n/(double)m_pStepsY[y];
		for (x=0;x<m_iRes[0];x++)
			m_pBin[x+y*m_iRes[0]] *= f;
	}
	BTOUT;
}


void C2DF::StepY(int y)
{
	BXIN;
	m_pStepsY[y]++;
	BXOUT;
}


void C2DF::CalcMaxEntry()
{
	BTIN;
	int z, z2;
	
	m_fMaxEntry = -99999999.0;
	m_fMinEntry =  99999999.0;
	for (z=0;z<m_iRes[1];z++)
	{
		for (z2=0;z2<m_iRes[0];z2++)
		{
			if (m_pBin[z*m_iRes[0]+z2] > m_fMaxEntry)
				m_fMaxEntry = m_pBin[z*m_iRes[0]+z2];
			if (m_pBin[z*m_iRes[0]+z2] < m_fMinEntry)
				m_fMinEntry = m_pBin[z*m_iRes[0]+z2];
		}
	}
	BTOUT;
}


void C2DF::SetLabelX(const char *s)
{
	if (m_sLabelX != NULL)
		delete[] m_sLabelX;

	try { m_sLabelX = new char[strlen(s)+1]; } catch(...) { m_sLabelX = NULL; }
	if (m_sLabelX == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	strcpy(m_sLabelX,s);
}


void C2DF::SetLabelY(const char *s)
{
	if (m_sLabelY != NULL)
		delete[] m_sLabelY;

	try { m_sLabelY = new char[strlen(s)+1]; } catch(...) { m_sLabelY = NULL; }
	if (m_sLabelY == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	strcpy(m_sLabelY,s);
}


void C2DF::SetLabelZ(const char *s)
{
	if (m_sLabelZ != NULL)
		delete[] m_sLabelZ;

	try { m_sLabelZ = new char[strlen(s)+1]; } catch(...) { m_sLabelZ = NULL; }
	if (m_sLabelZ == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	strcpy(m_sLabelZ,s);
}


void C2DF::MakeTensorProduct(C2DF *inp)
{
	int x, y;
	double tf, *px, *py;

	try { px = new double[m_iRes[0]]; } catch(...) { px = NULL; }
	if (px == NULL) NewException((double)m_iRes[0]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { py = new double[m_iRes[1]]; } catch(...) { py = NULL; }
	if (py == NULL) NewException((double)m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (y=0;y<m_iRes[1];y++)
	{
		tf = 0;
		for (x=0;x<m_iRes[0];x++)
			tf += inp->m_pBin[y*m_iRes[0]+x];
		tf /= (double)m_iRes[0];
		py[y] = tf;
	}
	for (x=0;x<m_iRes[0];x++)
	{
		tf = 0;
		for (y=0;y<m_iRes[1];y++)
			tf += inp->m_pBin[y*m_iRes[0]+x];
		tf /= (double)m_iRes[1];
		px[x] = tf;
	}
	for (x=0;x<m_iRes[0];x++)
	{
		for (y=0;y<m_iRes[1];y++)
			m_pBin[y*m_iRes[0]+x] = px[x] * py[y];
	}
	delete[] py;
	delete[] px;
}


double C2DF::CalcCorrelationFactor()
{
	int x, y;
	double tf, *px, *py, sum, tsum, cf;

	try { px = new double[m_iRes[0]]; } catch(...) { px = NULL; }
	if (px == NULL) NewException((double)m_iRes[0]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { py = new double[m_iRes[1]]; } catch(...) { py = NULL; }
	if (py == NULL) NewException((double)m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	sum = 0;
	for (x=0;x<m_iRes[0]*m_iRes[1];x++)
		sum += m_pBin[x];

	for (y=0;y<m_iRes[1];y++)
	{
		tf = 0;
		for (x=0;x<m_iRes[0];x++)
			tf += m_pBin[y*m_iRes[0]+x];
		tf /= (double)m_iRes[0];
		py[y] = tf;
	}

	for (x=0;x<m_iRes[0];x++)
	{
		tf = 0;
		for (y=0;y<m_iRes[1];y++)
			tf += m_pBin[y*m_iRes[0]+x];
		tf /= (double)m_iRes[1];
		px[x] = tf;
	}

	tsum = 0;
	for (x=0;x<m_iRes[0];x++)
		for (y=0;y<m_iRes[1];y++)
			tsum += px[x] * py[y];

	cf = 0;
	for (x=0;x<m_iRes[0];x++)
		for (y=0;y<m_iRes[1];y++)
			cf += fabs( m_pBin[y*m_iRes[0]+x]/sum - (px[x]*py[y])/tsum );

	delete[] py;
	delete[] px;

	return cf;
}


void C2DF::CopyFrom(C2DF *df)
{
	long z;

	if (df->m_sLabelX != NULL)
	{
		try { m_sLabelX = new char[strlen(df->m_sLabelX)+1]; } catch(...) { m_sLabelX = NULL; }
		if (m_sLabelX == NULL) NewException((double)(strlen(df->m_sLabelX)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(m_sLabelX,df->m_sLabelX);
	}
	if (df->m_sLabelY != NULL)
	{
		try { m_sLabelY = new char[strlen(df->m_sLabelY)+1]; } catch(...) { m_sLabelY = NULL; }
		if (m_sLabelY == NULL) NewException((double)(strlen(df->m_sLabelY)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		strcpy(m_sLabelY,df->m_sLabelY);
	}
	if (df->m_sLabelZ != NULL)
	{
		try { m_sLabelZ = new char[strlen(df->m_sLabelZ)+1]; } catch(...) { m_sLabelZ = NULL; }
		if (m_sLabelZ == NULL) NewException((double)(strlen(df->m_sLabelZ)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		strcpy(m_sLabelZ,df->m_sLabelZ);
	}
	m_iRes[0] = df->m_iRes[0];
	m_iRes[1] = df->m_iRes[1];
	m_fMinVal[0] = df->m_fMinVal[0];
	m_fMinVal[1] = df->m_fMinVal[1];
	m_fMaxVal[0] = df->m_fMaxVal[0];
	m_fMaxVal[1] = df->m_fMaxVal[1];

	try { m_fCountX = new double[m_iRes[0]]; } catch(...) { m_fCountX = NULL; }
	if (m_fCountX == NULL) NewException((double)m_iRes[0]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_fCountY = new double[m_iRes[1]]; } catch(...) { m_fCountY = NULL; }
	if (m_fCountY == NULL) NewException((double)m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_pStepsY = new unsigned long[m_iRes[1]]; } catch(...) { m_pStepsY = NULL; }
	if (m_pStepsY == NULL) NewException((double)m_iRes[0]*sizeof(unsigned long),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iRes[0];z++)
		m_fCountX[z] = df->m_fCountX[z];

	for (z=0;z<m_iRes[1];z++)
	{
		m_fCountY[z] = df->m_fCountY[z];
		m_pStepsY[z] = df->m_pStepsY[z];
	}

	try { m_pBin = new double[m_iRes[0]*m_iRes[1]]; } catch(...) { m_pBin = NULL; }
	if (m_pBin == NULL) NewException((double)m_iRes[0]*m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
		m_pBin[z] = df->m_pBin[z];
	m_fFac[0] = m_iRes[0] / (m_fMaxVal[0]-m_fMinVal[0]);
	m_fFac[1] = m_iRes[1] / (m_fMaxVal[1]-m_fMinVal[1]);
	m_fBinEntries = df->m_fBinEntries;
	m_fSkipEntries = df->m_fSkipEntries;
//	m_fMathematicaColorScale = df->m_fMathematicaColorScale;
//	m_fMathematicaColorOffset = df->m_fMathematicaColorOffset;
	m_pChannels[0] = df->m_pChannels[0];
	m_pChannels[1] = df->m_pChannels[1];
}


void C2DF::Subtract(C2DF *df)
{
	int z;

	if ((m_iRes[0] != df->m_iRes[0]) || (m_iRes[1] != df->m_iRes[1]))
	{
		eprintf("C2DF::Subtract(C2DF*): Resolution of CDFs must be identical.\n");
		return;
	}
	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
		m_pBin[z] -= df->m_pBin[z];
}


void C2DF::WriteCombinedPlot(const char *prefix, const char *name, const char *suffix)
{
	CGrace *comb;
//	char buf[32768];
	CxString buf;
	int z;

	try { comb = new CGrace(); } catch(...) { comb = NULL; }
	if (comb == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pChannels[0]->CalcMinMax();
	m_pChannels[1]->CalcMinMax();

	comb->SetViewport(0.15f,0.8f,0.8f,0.95f);
	comb->CurrentGraph()->m_bTicks = false;
	comb->CurrentGraph()->m_bTickLabels = false;
	comb->CurrentGraph()->m_fFrameWidth = 0.5;
	comb->SetRangeX(m_pChannels[0]->m_fMinVal,m_pChannels[0]->m_fMaxVal);
	comb->SetRangeY(/*m_pChannels[0]->m_fMinEntry*/0.0,m_pChannels[0]->m_fMaxEntry+(m_pChannels[0]->m_fMaxEntry/*-m_pChannels[0]->m_fMinEntry*/)*0.1);
	comb->MakeTicks();
	if (g_iObsChannel[0] == 0) // RDF
		comb->AddLine(m_pChannels[0]->m_fMinVal,1,m_pChannels[0]->m_fMaxVal,1,1,1);
	comb->AddDataset();
	comb->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_pChannels[0]->m_iResolution*2);
	for (z=0;z<m_pChannels[0]->m_iResolution-1;z++)
		comb->AddXYTupel(0,m_pChannels[0]->m_fMinVal+(z+0.5)*(m_pChannels[0]->m_fMaxVal-m_pChannels[0]->m_fMinVal)/m_pChannels[0]->m_iResolution,m_pChannels[0]->m_pBin[z]);

	comb->AddGraph();
	comb->SetViewport(0.8f, 0.15f, 0.95f, 0.8f);
	comb->CurrentGraph()->m_bTicks = false;
	comb->CurrentGraph()->m_bTickLabels = false;
	comb->CurrentGraph()->m_fFrameWidth = 0.5;
	comb->SetRangeX(/*m_pChannels[1]->m_fMinEntry*/0.0,m_pChannels[1]->m_fMaxEntry+(m_pChannels[1]->m_fMaxEntry/*-m_pChannels[1]->m_fMinEntry*/)*0.1);
	comb->SetRangeY(m_pChannels[1]->m_fMinVal,m_pChannels[1]->m_fMaxVal);
	comb->MakeTicks();
	if (g_iObsChannel[1] == 0) // RDF
		comb->AddLine(1,m_pChannels[1]->m_fMinVal,1,m_pChannels[1]->m_fMaxVal,1,1);
	comb->AddDataset();
	comb->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_pChannels[1]->m_iResolution*2);
	for (z=0;z<m_pChannels[1]->m_iResolution-1;z++)
		comb->AddXYTupel(0,m_pChannels[1]->m_pBin[z],m_pChannels[1]->m_fMinVal+(z+0.5)*(m_pChannels[1]->m_fMaxVal-m_pChannels[1]->m_fMinVal)/m_pChannels[1]->m_iResolution);

	comb->AddGraph();
	comb->SetViewport(0.15f, 0.15f, 0.8f, 0.8f);
	comb->CurrentGraph()->m_bTicksBothSidesX = true;
	comb->CurrentGraph()->m_bTicksBothSidesY = true;
	comb->CurrentGraph()->m_bTickInX = false;
	comb->CurrentGraph()->m_bTickInY = false;
	comb->SetLabelX(m_sLabelX);
	comb->SetLabelY(m_sLabelY);
	comb->SetRangeX(m_pChannels[0]->m_fMinVal,m_pChannels[0]->m_fMaxVal);
	comb->SetRangeY(m_pChannels[1]->m_fMinVal,m_pChannels[1]->m_fMaxVal);
	comb->MakeTicks();
	comb->AddDataset();

/*	strcpy(buf,prefix);
	strcat(buf,name);
	strcat(buf,suffix);*/
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);
	comb->WriteAgr(buf,false);
}


void C2DF::Mirror(float plane, int channel)
{
	int x, y, t;
	double *tbin, p;

	try { tbin = new double[m_iRes[0]*m_iRes[1]]; } catch(...) { tbin = NULL; }
	if (tbin == NULL) NewException((double)m_iRes[0]*m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (channel == 0)
	{
		p = (plane-m_fMinVal[0])*m_fFac[0];
		for (x=0;x<m_iRes[0];x++)
		{
			for (y=0;y<m_iRes[1];y++)
			{
				tbin[y*m_iRes[0]+x] = m_pBin[y*m_iRes[0]+x];
				t = (int)(2.0*p - x - 1);
				if ((t >= 0) && (t < m_iRes[0]))
				{
					tbin[y*m_iRes[0]+x] += m_pBin[y*m_iRes[0]+t];
					tbin[y*m_iRes[0]+x] /= 2.0;
				}
			}
		}
	} else if (channel == 1)
	{
		p = (plane-m_fMinVal[1])*m_fFac[1];
		for (x=0;x<m_iRes[0];x++)
		{
			for (y=0;y<m_iRes[1];y++)
			{
				tbin[y*m_iRes[0]+x] = m_pBin[y*m_iRes[0]+x];
				t = (int)(2*p - y - 1);
				if ((t >= 0) && (t < m_iRes[1]))
				{
					tbin[y*m_iRes[0]+x] += m_pBin[t*m_iRes[0]+x];
					tbin[y*m_iRes[0]+x] /= 2.0;
				}
			}
		}
	}

	delete[] m_pBin;
	m_pBin = tbin;
}


void C2DF::SwapAxes()
{
	char *t;
	int x, y;
	double *pd;
	double d;
	int i;
	CDF *df;

	t = m_sLabelX;
	m_sLabelX = m_sLabelY;
	m_sLabelY = t;

	try { pd = new double[m_iRes[0]*m_iRes[1]]; } catch(...) { pd = NULL; }
	if (pd == NULL) NewException((double)m_iRes[0]*m_iRes[1]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (y=0;y<m_iRes[1];y++)
		for (x=0;x<m_iRes[0];x++)
			pd[x*m_iRes[1]+y] = m_pBin[y*m_iRes[0]+x];
	delete[] m_pBin;
	m_pBin = pd;
	d = m_fMinVal[0];
	m_fMinVal[0] = m_fMinVal[1];
	m_fMinVal[1] = d;
	d = m_fMaxVal[0];
	m_fMaxVal[0] = m_fMaxVal[1];
	m_fMaxVal[1] = d;
	i = m_iRes[0];
	m_iRes[0] = m_iRes[1];
	m_iRes[1] = i;
	df = m_pChannels[0];
	m_pChannels[0] = m_pChannels[1];
	m_pChannels[1] = df;
	d = m_fFac[0];
	m_fFac[0] = m_fFac[1];
	m_fFac[1] = d;
}


void C2DF::NormalizeXCount()
{
	int x, y;

	for (x=0;x<m_iRes[0];x++)
		for (y=0;y<m_iRes[1];y++)
			if (m_fCountX[x] != 0)
				m_pBin[y*m_iRes[0]+x] /= m_fCountX[x];
}


void C2DF::NormalizeYCount()
{
	int x, y;

	for (y=0;y<m_iRes[1];y++)
		for (x=0;x<m_iRes[0];x++)
			if (m_fCountY[y] != 0)
				m_pBin[y*m_iRes[0]+x] /= m_fCountY[y];
}


void C2DF::AddCircle(double x, double y, double r, double cr, double cg, double cb)
{
	CMathematicaCircle *c;

	try { c = new CMathematicaCircle(); } catch(...) { c = NULL; }
	if (c == NULL) NewException((double)sizeof(CMathematicaCircle),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	c->m_fPosX = x;
	c->m_fPosY = y;
	c->m_fRadius = r;
	c->m_fColorR = cr;
	c->m_fColorG = cg;
	c->m_fColorB = cb;
	m_oaCircles.Add(c);
}


double C2DF::GetValue(double x, double y)
{
	BXIN;
	double rx, ry, r;
	int ix, iy;

	if ((x < m_fMinVal[0]) || (y < m_fMinVal[1]) || (x > m_fMaxVal[0]) || (y > m_fMaxVal[1]))
	{
		BXOUT;
		return 0;
	}

	rx = (x-m_fMinVal[0])*m_fFac[0] - 0.5;
	ix = (int)floor(rx);
	if (ix < 0)
	{
		ix = 0;
		rx = 0;
	} else if (ix > m_iRes[0]-2)
	{
		ix = m_iRes[0]-2;
		rx = 1.0;
	} else
		rx -= ix;

	ry = (y-m_fMinVal[1])*m_fFac[1] - 0.5;
	iy = (int)floor(ry);
	if (iy < 0)
	{
		iy = 0;
		ry = 0;
	} else if (iy > m_iRes[1]-2)
	{
		iy = m_iRes[1]-2;
		ry = 1.0;
	} else
		ry -= iy;

	r  = m_pBin[ iy    * m_iRes[0] + ix    ] * (1.0-rx) * (1.0-ry);
	r += m_pBin[ iy    * m_iRes[0] + ix + 1] *      rx  * (1.0-ry);
	r += m_pBin[(iy+1) * m_iRes[0] + ix    ] * (1.0-rx) *      ry;
	r += m_pBin[(iy+1) * m_iRes[0] + ix + 1] *      rx  *      ry;

	BXOUT;
	return r;
}


void C2DF::AddToBin(int x, int y, double val)
{
	BXIN;

	m_fBinEntries++;
	m_pBin[y*m_iRes[0] + x] += val;
	m_fCountX[x]++;
	m_fCountY[y]++;
	BXOUT;
}


void C2DF::CalcHistogram()
{
	int z, z2, t;

	CalcMaxEntry();

	try { m_pHistogram = new double[m_iHistogramRes]; } catch(...) { m_pHistogram = NULL; }
	if (m_pHistogram == NULL) NewException((double)m_iHistogramRes*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iHistogramRes;z++)
		m_pHistogram[z] = 0;
	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
	{
		t = (int)floor((m_pBin[z]-m_fMinEntry)/(m_fMaxEntry-m_fMinEntry)*m_iHistogramRes);
		for (z2=0;z2<t;z2++)
			m_pHistogram[z2]++;
	}
	for (z=0;z<m_iHistogramRes;z++)
		m_pHistogram[z] /= ((double)m_iRes[0])*m_iRes[1];
}


void C2DF::WriteHistogram(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z;
//	char buf[32768];
	CxString buf;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);
	
	a = OpenFileWrite(buf,true);

	for (z=0;z<m_iHistogramRes;z++)
		mfprintf(a," %#.10G;  %#.10G\n",(z+0.5)*(m_fMaxEntry-m_fMinEntry)/m_iHistogramRes,m_pHistogram[z]);
	
	fclose(a);
	BTOUT;
}


void C2DF::Log()
{
	int z;

	for (z=0;z<m_iRes[0]*m_iRes[1];z++)
	{
		if (m_pBin[z] >= 1.0E-10)
			m_pBin[z] = log10(m_pBin[z])+10.0;
				else m_pBin[z] = 0;
	}
}


void C2DF::NormalizeUniform(double fac)
{
	int zx, zy;
	double r1, r2, a1, a2, dca, vol;

	for (zy=0;zy<m_iRes[1];zy++)
	{
		a1 = (m_fMinVal[1] + (double)zy     * (m_fMaxVal[1] - m_fMinVal[1]) / m_iRes[1]) * Pi / 180.0;
		a2 = (m_fMinVal[1] + (double)(zy+1) * (m_fMaxVal[1] - m_fMinVal[1]) / m_iRes[1]) * Pi / 180.0;

		dca = (1.0-cos(a2)) - (1.0-cos(a1));

		for (zx=0;zx<m_iRes[0];zx++)
		{
			r1 = m_fMinVal[0] + ((double)zx) / m_iRes[0] * (m_fMaxVal[0] - m_fMinVal[0]);
			r2 = m_fMinVal[0] + (zx+1.0)     / m_iRes[0] * (m_fMaxVal[0] - m_fMinVal[0]);

			vol = 2.0/3.0*Pi * ( pow(r2,3) * dca - pow(r1,3) * dca );

			m_pBin[zy*m_iRes[0]+zx] *= fac / vol;
		}
	}
}
