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

#include "bintools.h"
#include "travis.h"

/*void CMSDF::CalcHistogram()
{
	BTIN;
	int z, z2, t;

	m_fMaxAvg = 0;
	for (z=0;z<g_iSDFResTri;z++)
		if (m_fMaxAvg < m_pBin[z])
			m_fMaxAvg = m_pBin[z];

	m_pHistogram = new double[g_iHistogramRes];
	for (z=0;z<g_iHistogramRes;z++)
	{
		t = 0;
		for (z2=0;z2<g_iSDFResTri;z2++)
		{
			if ((m_pBin[z2] > z*m_fMaxAvg/(g_iHistogramRes-1)) && (m_pBin[z2] <= (z+1)*m_fMaxAvg/(g_iHistogramRes-1)))
				t++;
		}
		m_pHistogram[z] = (double)t;
	}
	BTOUT;
}


void CMSDF::NormHistogramIntegral()
{
	BTIN;
	int z;
	double f;

	f = 0;
	for (z=0;z<g_iHistogramRes;z++)
		f += m_pHistogram[z];

	for (z=0;z<g_iHistogramRes;z++)
		m_pHistogram[z] /= f;
	BTOUT;
}


void CMSDF::WriteHistogram(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z;
	char buf[256];
	
	strcpy(buf,prefix);
	strcat(buf,name);
	strcat(buf,suffix);
	
	a = OpenFileWrite(buf,true);

	for (z=0;z<g_iHistogramRes;z++)
		mfprintf(a,"%10f;  %10f\n",z*m_fMaxAvg/(g_iHistogramRes-1),m_pHistogram[z]);
	
	fclose(a);
	BTOUT;
}*/

/*void CMSDF::CalcAvg()
{
	BTIN;
	int z;

	for (z=0;z<g_iSDFResTri;z++)
	{
		if (m_pRefBin[z] < 1E-20)
			continue;
		m_pBin[z] /= m_pRefBin[z];
//		m_pDBin[0][z] /= m_pRefBin[z];
//		m_pDBin[1][z] /= m_pRefBin[z];
//		m_pDBin[2][z] /= m_pRefBin[z];
	}
	BTOUT;
}*/


/*void CS6DF::CalcHistogram()
{
	BTIN;
	int z, z2, z0, t;

	for (z0=0;z0<m_iLevelCount;z0++)
	{
		m_fMaxP[z0] = 0;
		for (z=0;z<g_iSDFResTri;z++)
			if (m_fMaxP[z0] < m_pBin[z0][z])
				m_fMaxP[z0] = m_pBin[z0][z];

		m_pHistogram[z0] = new double[g_iHistogramRes];
		for (z=0;z<g_iHistogramRes;z++)
		{
			t = 0;
			for (z2=0;z2<g_iSDFResTri;z2++)
				if ((m_pBin[z0][z2] > z*m_fMaxP[z0]/(g_iHistogramRes-1)) && (m_pBin[z0][z2] <= (z+1)*m_fMaxP[z0]/(g_iHistogramRes-1)))
					t++;
			m_pHistogram[z0][z] = (double)t;
		}
	}
	BTOUT;
}*/

/*void CS6DF::WriteHistogram(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z, z0;
	char buf[256], buf2[32];
	
	for (z0=0;z0<m_iLevelCount;z0++)
	{
		strcpy(buf,prefix);
		strcat(buf,name);
		sprintf(buf2,".%d",z0+1);
		strcat(buf,buf2);
		strcat(buf,suffix);
	
		a = fopen(buf,"wt");

		for (z=0;z<g_iHistogramRes;z++)
			mfprintf(a,"%10f;  %10f\n",z*m_fMaxP[z0]/(g_iHistogramRes-1),m_pHistogram[z0][z]);
	
		fclose(a);
	}
	BTOUT;
}*/


/*void CS6DF::AddToBin(const CxVector3 &vec, const CxVector3 &val)
{
	BXIN;
	double x, y, z, tx, ty, tz, d;
	int ix, iy, iz, j;

	if (vec.GetLength() > m_fRadius)
	{
		BXOUT;
		return;
	}
	if (m_bCutSDF)
	{
		switch(m_iCutPlane)
		{
			case 0:
				if (vec[0] > g_fCutValue)
				{
					BXOUT;
					return;
				}
				break;
			case 1:
				if (vec[1] > g_fCutValue)
				{
					BXOUT;
					return;
				}
				break;
			case 2:
				if (vec[2] > g_fCutValue)
				{
					BXOUT;
					return;
				}
				break;
		}
	}
	d = val.GetLength();
//	printf("AddToBin x=%f y=%f z=%f d=%f Level=%d\n",vec[0],vec[1],vec[2],d,m_iLevelCount);
	m_fGesBinEntries++;
	x = ((vec[0]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	y = ((vec[1]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	z = ((vec[2]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	ix = (int)floor(x);
	iy = (int)floor(y);
	iz = (int)floor(z);

	j = 0;
	while (d > m_fLevelThres[j])
	{
		j++;
		if (j >= m_iLevelCount)
			break;
	}
	if (j >= m_iLevelCount)
		j = m_iLevelCount-1;
//	printf("  Ok. x=%f y=%f z=%f d=%f j=%d Level=%d\n",vec[0],vec[1],vec[2],d,j,g_iLevelCount);
	m_fBinEntries[j]++;
//	switch(g_iBinning)
//	{
//		case 1: // Einfach einsortieren
//			m_pBin[j][iz*g_iSDFResSqr + iy*g_iSDFResolution + ix] += 1.0f;
//			break;
//		case 2: // Auf 8 Bins aufspalten
			x -= ix;
			y -= iy;
			z -= iz;
			tx = 1-x;
			ty = 1-y;
			tz = 1-z;
			m_pBin[j][iz*g_iSDFResSqr + iy*g_iSDFResolution + ix]            += tx*ty*tz;
			m_pBin[j][iz*g_iSDFResSqr + iy*g_iSDFResolution + ix +1]         +=  x*ty*tz;
			m_pBin[j][iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix]        += tx* y*tz;
			m_pBin[j][iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1]     +=  x* y*tz;
			m_pBin[j][(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix]        += tx*ty* z;
			m_pBin[j][(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix +1]     +=  x*ty* z;
			m_pBin[j][(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix]    += tx* y* z;
			m_pBin[j][(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1] +=  x* y* z;
//			break;
//		default:
//			printf("CS6DF::AddToBin(): Error!\n");
//			break;
//	}
//	printf("  Fertig.\n");
	BXOUT;
}*/

/*void CMSDF::AddToBinMean(const CxVector3 &vec, const CxVector3 &val)
{
	BXIN;
	double x, y, z, tx, ty, tz, v;
	int ix, iy, iz;

	if (vec.GetLength() > g_fSDFRadius)
	{
		BXOUT;
		return;
	}
	if (g_bCutSDF)
	{
		switch(g_iCutPlane)
		{
			case 0:
				if (vec[0] > g_fCutValue)
				{
					BXOUT;
					return;
				}
				break;
			case 1:
				if (vec[1] > g_fCutValue)
				{	
					BXOUT;
					return;
				}
				break;
			case 2:
				if (vec[2] > g_fCutValue)
				{
					BXOUT;
					return;
				}
				break;
		}
	}
	v = val.GetLength();
	m_fBinEntries++;
	x = ((vec[0]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	y = ((vec[1]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	z = ((vec[2]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	ix = (int)floor(x);
	iy = (int)floor(y);
	iz = (int)floor(z);
//	printf("AddToBin x=%f y=%f z=%f d=%f\n",vec[0],vec[1],vec[2],v);
//	switch(g_iBinning)
//	{
//		case 1: // Einfach einsortieren
//			m_pBin   [iz*g_iSDFResSqr + iy*g_iSDFResolution + ix] += v;
//			m_pRefBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix] += 1.0;
//			break;
//		case 2: // Auf 8 Bins aufspalten
			x -= ix;
			y -= iy;
			z -= iz;
			tx = 1-x;
			ty = 1-y;
			tz = 1-z;
			m_pBin   [iz*g_iSDFResSqr + iy*g_iSDFResolution + ix]            += tx*ty*tz*v;
			m_pBin   [iz*g_iSDFResSqr + iy*g_iSDFResolution + ix +1]         +=  x*ty*tz*v;
			m_pBin   [iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix]        += tx* y*tz*v;
			m_pBin   [iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1]     +=  x* y*tz*v;
			m_pBin   [(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix]        += tx*ty* z*v;
			m_pBin   [(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix +1]     +=  x*ty* z*v;
			m_pBin   [(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix]    += tx* y* z*v;
			m_pBin   [(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1] +=  x* y* z*v;

			m_pRefBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix]            += tx*ty*tz;
			m_pRefBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix +1]         +=  x*ty*tz;
			m_pRefBin[iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix]        += tx* y*tz;
			m_pRefBin[iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1]     +=  x* y*tz;
			m_pRefBin[(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix]        += tx*ty* z;
			m_pRefBin[(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix +1]     +=  x*ty* z;
			m_pRefBin[(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix]    += tx* y* z;
			m_pRefBin[(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1] +=  x* y* z;
//			break;
//		default:
//			printf("CMSDF::AddToBin(): Error!\n");
//			break;
//	}
	BXOUT;
}*/

/*void CMSDF::AddToBinMax(const CxVector3 &vec, const CxVector3 &val)
{
	BXIN;
	double x, y, z, tx, ty, tz, v;
	int ix, iy, iz;

	if (vec.GetLength() > g_fSDFRadius)
	{
		BXOUT;
		return;
	}
	if (g_bCutSDF)
	{
		switch(g_iCutPlane)
		{
			case 0:
				if (vec[0] > g_fCutValue)
				{
					BXOUT;
					return;
				}
				break;
			case 1:
				if (vec[1] > g_fCutValue)
				{
					BXOUT;
					return;
				}
				break;
			case 2:
				if (vec[2] > g_fCutValue)
				{
					BXOUT;
					return;
				}
				break;
		}
	}
	v = val.GetLength();
	m_fBinEntries++;
	x = ((vec[0]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	y = ((vec[1]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	z = ((vec[2]+g_fSDFRadius)/(2*g_fSDFRadius))*((double)g_iSDFResolution-1);
	ix = (int)floor(x);
	iy = (int)floor(y);
	iz = (int)floor(z);
//	printf("AddToBin x=%f y=%f z=%f d=%f\n",vec[0],vec[1],vec[2],v);
//	switch(g_iBinning)
//	{
//		case 1: // Einfach einsortieren
//			if (m_pBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix] < v)
//				m_pBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix] = v;
//			break;
//		case 2: // Auf 8 Bins aufspalten
			x -= ix;
			y -= iy;
			z -= iz;
			tx = 1-x;
			ty = 1-y;
			tz = 1-z;
			if (m_pBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix] < tx*ty*tz*v)
				m_pBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix] = tx*ty*tz*v;
			if (m_pBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix +1] < x*ty*tz*v)
				m_pBin[iz*g_iSDFResSqr + iy*g_iSDFResolution + ix +1] =  x*ty*tz*v;
			if (m_pBin[iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix] < tx* y*tz*v)
				m_pBin[iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix] = tx* y*tz*v;
			if (m_pBin[iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1] < x* y*tz*v)
				m_pBin[iz*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1] =  x* y*tz*v;
			if (m_pBin[(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix] < tx*ty* z*v)
				m_pBin[(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix] = tx*ty* z*v;
			if (m_pBin[(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix +1] < x*ty* z*v)
				m_pBin[(iz+1)*g_iSDFResSqr + iy*g_iSDFResolution + ix +1] =  x*ty* z*v;
			if (m_pBin[(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix] < tx* y* z*v)
				m_pBin[(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix] = tx* y* z*v;
			if (m_pBin[(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1] < x* y* z*v)
				m_pBin[(iz+1)*g_iSDFResSqr + (iy+1)*g_iSDFResolution + ix +1] =  x* y* z*v;
//			break;
//		default:
//			printf("CMSDF::AddToBin(): Error!\n");
//			break;
//	}
	BXOUT;
}*/

/*
void CMRDF::AddToBin(double d, double val[3])
{
	double p, v;
	int ip;

	if (d > g_fRDFRadius)
		return;
	v = length(val);
	m_fBinEntries++;
	p = (d/g_fRDFRadius)*(g_iRDFResolution-1.0);
	ip = floor(p);
	switch(g_iBinning)
	{
		case 1: // Einfach einsortieren
			m_pBin   [ip] +=   v;	
			m_pRefBin[ip] += 1.0;	
			break;
		case 2: // Auf 2 Bins aufspalten
			p -= ip;
			m_pBin   [ip]   += v*(1-p);
			m_pBin   [ip+1] += v*   p;
			m_pRefBin[ip]   += (1-p);
			m_pRefBin[ip+1] +=    p;
			break;
		default:
			printf("CMRDF::AddToBin(): Error!\n");
			break;
	}
}
*/


CAF::CAF()
{
	m_pBin = NULL;
	m_pDBin = NULL;
	m_faEntries.SetName("CAF::m_faEntries");
}


CAF::~CAF()
{
}


void CAF::AddToBin(double x, double y)
{
	BXIN;
	double p;
	int ip;

	if ((x < m_fMinVal) || (x > m_fMaxVal))
	{
		BXOUT;
		return;
	}
	m_fBinEntries++;
	p = ((x-m_fMinVal)/(m_fMaxVal-m_fMinVal))*(m_iResolution-1.0f);
	ip = (int)floor(p);
//	switch(g_iBinning)
//	{
//		case 1: // Einfach einsortieren
//			m_pBin[ip] += y;
//			m_faEntries[ip]++;
//			break;
//		case 2: // Auf 2 Bins aufspalten
			p -= ip;
			m_pBin[ip]   += (1-p)*y;
			m_pBin[ip+1] +=    p *y;
			m_faEntries[ip]   += (1-p);
			m_faEntries[ip+1] +=    p;
//			break;
//		default:
//			printf("CDF::AddToBin(): Error!\n");
//			break;
//	}
	BXOUT;
}

void CAF::AddToBin_Index(int i, double y)
{
	BXIN;

	if ((i < 0) || (i >= m_iResolution))
	{
		BXOUT;
		return;
	}
	m_fBinEntries++;
//	switch(g_iBinning)
//	{
//		case 1: // Einfach einsortieren
//			m_pBin[ip] += y;
//			m_faEntries[ip]++;
//			break;
//		case 2: // Auf 2 Bins aufspalten
			m_pBin[i] += y;
			m_faEntries[i]++;
//			break;
//		default:
//			printf("CDF::AddToBin(): Error!\n");
//			break;
//	}
	BXOUT;
}

void CAF::BuildAverage()
{
	BTIN;
	int z;

	for (z=0;z<m_iResolution;z++)
	{
		if (m_faEntries[z] != 0)
			m_pBin[z] /= m_faEntries[z];
	}
	BTOUT;
}

/*double CS6DF::PPMBin()
{
	BTIN;
	double f;
	int z, z2;

	f = 0;
	for (z2=0;z2<m_iLevelCount;z2++)
	{
		for (z=0;z<g_iSDFResTri;z++)
		{
			m_pBin[z2][z] *= 1000000.0f / (double)g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_waSingleMolIndex.GetSize();
			if (f < m_pBin[z2][z])
				f = m_pBin[z2][z];
		}
	}
	BTOUT;
	return f;
}*/


/*void CMSDF::Create()
{
	BTIN;
	int z;
	m_fBinEntries = 0;
	m_pBin = new double[g_iSDFResTri];
	m_pRefBin = new double[g_iSDFResTri];
//	m_pDBin[0] = new double[g_iSDFResTri];
//	m_pDBin[1] = new double[g_iSDFResTri];
//	m_pDBin[2] = new double[g_iSDFResTri];

	for (z=0;z<g_iSDFResTri;z++)
	{
		m_pBin[z] = 0;
		m_pRefBin[z] = 0;
//		m_pDBin[0][z] = 0;
//		m_pDBin[1][z] = 0;
//		m_pDBin[2][z] = 0;
	}
	m_pHistogram = NULL;
	BTOUT;
}*/

/*void CS6DF::Create()
{
	BTIN;
	int z, z2;
	m_fGesBinEntries = 0;
	for (z=0;z<m_iLevelCount;z++)
	{
		m_fBinEntries[z] = 0;
		m_pBin[z] = new double[g_iSDFResTri];
//		m_pDBin[0][z] = new double[g_iSDFResTri];
//		m_pDBin[1][z] = new double[g_iSDFResTri];
//		m_pDBin[2][z] = new double[g_iSDFResTri];
		for (z2=0;z2<g_iSDFResTri;z2++)
		{
			m_pBin[z][z2] = 0;
//			m_pDBin[0][z][z2] = 0;
//			m_pDBin[1][z][z2] = 0;
//			m_pDBin[2][z][z2] = 0;
		}
		m_pHistogram[z] = NULL;
//		m_pDHistogram[z] = NULL;
	}
	BTOUT;
}*/

void CAF::Create()
{
	BTIN;
	int z;
	m_fBinEntries = 0.0;

	try { m_pBin = new double[m_iResolution]; } catch(...) { m_pBin = NULL; }
	if (m_pBin == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_faEntries.SetSize(m_iResolution);
	for (z=0;z<m_iResolution;z++)
	{
		m_pBin[z] = 0.0f;
		m_faEntries[z] = 0.0;
	}
	BTOUT;
}



/*void CMSDF::Write(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	unsigned int i;
	int z;
	float f;
	char buf[256];
	
	strcpy(buf,prefix);
	strcat(buf,name);
	strcat(buf,suffix);

	a = fopen(buf,"wb");
	if (a == NULL)
	{
		mprintf("\nFehler! Konnte Ausgabedatei %s nicht zum Schreiben oeffnen.\n",buf);
		BTOUT;
		return;
	}
	i = 3;
	fwrite(&i,4,1,a);
	i = 50;
	fwrite(&i,4,1,a);
	i = g_iSDFResolution;
	fwrite(&i,4,1,a);
	fwrite(&i,4,1,a);
	fwrite(&i,4,1,a);
	f = -g_fSDFRadius;
	fwrite(&f,4,1,a);
	f = g_fSDFRadius;
	fwrite(&f,4,1,a);
	f = -g_fSDFRadius;
	fwrite(&f,4,1,a);
	f = g_fSDFRadius;
	fwrite(&f,4,1,a);
	f = -g_fSDFRadius;
	fwrite(&f,4,1,a);
	f = g_fSDFRadius;
	fwrite(&f,4,1,a);
	for (z=0;z<g_iSDFResTri;z++)
	{
		f = (float)m_pBin[z];
		fwrite(&f,4,1,a);
	}
	fclose(a);
	BTOUT;
}*/

/*void CS6DF::Write(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	unsigned int i;
	int z, z0;
	float f;
	char buf[256], buf2[32];

	for (z0=0;z0<m_iLevelCount;z0++)
	{
		strcpy(buf,prefix);
		strcat(buf,name);
		sprintf(buf2,".%d",z0+1);
		strcat(buf,buf2);
		strcat(buf,suffix);

		a = fopen(buf,"wb");
		if (a == NULL)
		{
			mprintf("\nFehler! Konnte Ausgabedatei %s nicht zum Schreiben Oeffnen.\n",buf);
			BTOUT;
			return;
		}
		i = 3;
		fwrite(&i,4,1,a);
		i = 50;
		fwrite(&i,4,1,a);
		i = g_iSDFResolution;
		fwrite(&i,4,1,a);
		fwrite(&i,4,1,a);
		fwrite(&i,4,1,a);
		f = -g_fSDFRadius;
		fwrite(&f,4,1,a);
		f = g_fSDFRadius;
		fwrite(&f,4,1,a);
		f = -g_fSDFRadius;
		fwrite(&f,4,1,a);
		f = g_fSDFRadius;
		fwrite(&f,4,1,a);
		f = -g_fSDFRadius;
		fwrite(&f,4,1,a);
		f = g_fSDFRadius;
		fwrite(&f,4,1,a);
		for (z=0;z<g_iSDFResTri;z++)
		{
			f = (float)m_pBin[z0][z];
			fwrite(&f,4,1,a);
		}
		fclose(a);
	}
	BTOUT;
}*/


void CAF::Write(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z;
//	char buf[256];
	CxString buf;
	
//	buf[0] = 0;
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	mfprintf(a,"# tau [ps];  MSD [pm^2];  Derivative\n");
	if (m_pDBin != NULL)
	{
		for (z=0;z<m_iResolution;z++)
			mfprintf(a,"%f; %f; %f\n",m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution,m_pBin[z],m_pDBin[z]);
	} else
	{
		for (z=0;z<m_iResolution;z++)
			mfprintf(a,"%f; %f\n",m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution,m_pBin[z]);
	}
	fclose(a);
	BTOUT;
}

/*CS6DF::CS6DF()
{
}

CS6DF::~CS6DF()
{
	BTIN;
	int z;
	for (z=0;z<m_iLevelCount;z++)
	{
		delete m_pBin[z];
//		delete m_pDBin[0][z];
//		delete m_pDBin[1][z];
//		delete m_pDBin[2][z];
		if (m_pHistogram[z] != NULL)
			delete m_pHistogram[z];
//		if (m_pDHistogram[z] != NULL)
//			delete m_pDHistogram[z];
	}
	BTOUT;
}*/

void CAF::CalcDeriv(double f)
{
	BTIN;
	int z;

	try { m_pDBin = new double[m_iResolution]; } catch(...) { m_pDBin = NULL; }
	if (m_pDBin == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_iResolution-1;z++)
		m_pDBin[z] = (m_pBin[z+1]-m_pBin[z])*f/((m_fMaxVal-m_fMinVal)/m_iResolution);
	m_pDBin[m_iResolution-1] = 0.0f;
	BTOUT;
}


void CAF::LinReg(int i1, int i2, double *a0, double *a1, double *r)
{
	double x, mx, my, sxx, sxy, sse, sst;
	int z;

	mx = 0;
	my = 0;
	sxx = 0;
	sxy = 0;
	for (z=i1;z<=i2;z++)
	{
		x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
		mx += x;
		my += m_pBin[z];
		sxx += x*x;
		sxy += x*m_pBin[z];
	}
	mx /= i2-i1+1;
	my /= i2-i1+1;

	*a1 = (sxy - (i2-i1+1)*mx*my) / (sxx - (i2-i1+1)*mx*mx);
	*a0 = my - *a1 * mx;

	sse = 0;
	sst = 0;
	for (z=i1;z<=i2;z++)
	{
		x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
		sse += (m_pBin[z] - *a0 - *a1*x) * (m_pBin[z] - *a0 - *a1*x);
		sst += (m_pBin[z]-my) * (m_pBin[z]-my);
	}

	*r = sqrt(1.0 - sse/sst);
}
