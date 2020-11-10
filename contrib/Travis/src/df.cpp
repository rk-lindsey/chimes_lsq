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


#include "df.h"
#include "travis.h"


CDF::CDF()
{
	m_pHistogram = NULL;
	m_fBinEntries = 0;
	m_fSkipEntries = 0;
	m_fSum = 0;
	m_fSqSum = 0;
	m_fMean = 0;
	m_fSD = 0;
	m_fMinInput = 1e50;
	m_fMaxInput = -1e50;
	m_iMultiCount = 0;
	m_iHistogramRes = 0;
	m_pIntegral = NULL;
	m_pBin = 0;
	m_pIntegral = 0;
	m_iResolution = 0; 
	m_bLeft = false;
	m_pAdditionalSets = NULL;
	m_pAdditionalSetLabels = NULL;
	m_iAdditionalSets = 0;
	m_sLabelX = NULL;
	m_sLabelY = NULL;
	m_sLabelMulti = NULL;
	m_oaTimeDiffBuf.SetName("CDF::m_oaTimeDiffBuf");
}

	
CDF::~CDF()
{
}


void CDF::ZeroBin()
{
	int z;

	for (z=0;z<m_iResolution;z++)
		m_pBin[z] = 0;

	m_fBinEntries = 0;
	m_fSkipEntries = 0;
	m_fSum = 0;
	m_fSqSum = 0;
}


void CDF::AddToBin_Int(int i)
{
	double d;
	BXIN;

	d = ((double)i)/m_iResolution*(m_fMaxVal-m_fMinVal)-m_fMinVal;
	m_fSum += d;
	m_fSqSum += d*d;
	if (d < m_fMinInput)
		m_fMinInput = d;
	if (d > m_fMaxInput)
		m_fMaxInput = d;

	if ((i < 0) || (i >= m_iResolution))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;

	m_pBin[i]++;

	BXOUT;
}


void CDF::AddToBin(double d)
{
	BXIN;
	double p;
	int ip;

	m_fSum += d;
	m_fSqSum += d*d;
	if (d < m_fMinInput)
		m_fMinInput = d;
	if (d > m_fMaxInput)
		m_fMaxInput = d;

	if ((d < m_fMinVal) || (d > m_fMaxVal))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;


	p = (d-m_fMinVal)*m_fFac - 0.5;
	ip = (int)floor(p);
	if (ip < 0)
	{
		ip = 0;
		p = 0;
	} else if (ip > m_iResolution-2)
	{
		ip = m_iResolution-2;
		p = 1.0;
	} else
		p -= ip;

	m_pBin[ip    ] += (1.0-p);
	m_pBin[ip + 1] +=      p ;


/*	if (ip < 0)
	{
		m_pBin[0]++;
	} else if (ip > m_iResolution-2)
	{
		m_pBin[m_iResolution-1]++;
	} else
	{
		p -= ip;
		m_pBin[ip]   += (1-p);
		m_pBin[ip+1] +=    p;
	}*/

/*	ip = floor((d-m_fMinVal)*m_fFac);
	if (ip >= m_iResolution)
		ip = m_iResolution -1;
	m_pBin[ip]++;*/

	BXOUT;
}


void CDF::AddToBin(double d, double v)
{
	BXIN;
	double p;
	int ip;

	m_fSum += d;
	m_fSqSum += d*d;
	if (d < m_fMinInput)
		m_fMinInput = d;
	if (d > m_fMaxInput)
		m_fMaxInput = d;

	if ((d < m_fMinVal) || (d > m_fMaxVal))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;

	p = (d-m_fMinVal)*m_fFac - 0.5;
	ip = (int)floor(p);
	if (ip < 0)
	{
		ip = 0;
		p = 0;
	} else if (ip > m_iResolution-2)
	{
		ip = m_iResolution-2;
		p = 1.0;
	} else
		p -= ip;

	m_pBin[ip    ] += (1.0-p) * v;
	m_pBin[ip + 1] +=      p  * v;

	BXOUT;
}


void CDF::AddToBin_Multi(int i, double d)
{
	BXIN;
	double p;
	int ip;

	if (i >= m_iMultiCount)
	{
		eprintf("CDF::AddToBin_Multi(): %d >= %d.\n",i,m_iMultiCount);
		abort();
	}

	m_fSum += d;
	m_fSqSum += d*d;
	if (d < m_fMinInput)
		m_fMinInput = d;
	if (d > m_fMaxInput)
		m_fMaxInput = d;

	if ((d < m_fMinVal) || (d > m_fMaxVal))
	{
		m_fSkipEntries++;
		BXOUT;
		return;
	}
	m_fBinEntries++;

	p = (d-m_fMinVal)*m_fFac - 0.5;
	ip = (int)floor(p);
	if (ip < 0)
	{
		ip = 0;
		p = 0;
	} else if (ip > m_iResolution-2)
	{
		ip = m_iResolution-2;
		p = 1.0;
	} else
		p -= ip;

	m_pBin[ip    ] += (1.0-p);
	m_pBin[ip + 1] +=      p ;
	m_pMultiBin[i][ip    ] += (1.0-p);
	m_pMultiBin[i][ip + 1] +=      p ;


	BXOUT;
}


void CDF::AddToBin_Multi_Int(int i, int n, double f)
{
	BXIN;
	double d;

	if (i >= m_iMultiCount)
	{
		eprintf("CDF::AddToBin_Multi(): %d >= %d.\n",i,m_iMultiCount);
		abort();
	}

	d = ((double)n)/m_iResolution*(m_fMaxVal-m_fMinVal)-m_fMinVal;
	m_fSum += d;
	m_fSqSum += d*d;
	if (d < m_fMinInput)
		m_fMinInput = d;
	if (d > m_fMaxInput)
		m_fMaxInput = d;

	m_fBinEntries++;

	m_pBin[n] += f;
	m_pMultiBin[i][n] += f;

	BXOUT;
}


void CDF::SetAdditionalDatasetLabel(int z, const char *s)
{
	if (m_pAdditionalSetLabels[z] != NULL)
		delete[] m_pAdditionalSetLabels[z];

	try { m_pAdditionalSetLabels[z] = new char[strlen(s)+1]; } catch(...) { m_pAdditionalSetLabels[z] = NULL; }
	if (m_pAdditionalSetLabels[z] == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_pAdditionalSetLabels[z],s);
}


void CDF::SetLabelX(const char *s)
{
	if (m_sLabelX != NULL)
		delete[] m_sLabelX;

	try { m_sLabelX = new char[strlen(s)+1]; } catch(...) { m_sLabelX = NULL; }
	if (m_sLabelX == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sLabelX,s);
}


void CDF::SetLabelY(const char *s)
{
	if (m_sLabelY != NULL)
		delete[] m_sLabelY;

	try { m_sLabelY = new char[strlen(s)+1]; } catch(...) { m_sLabelY = NULL; }
	if (m_sLabelY == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sLabelY,s);
}


void CDF::AddToBin_Int(int i, double j)
{
	double d;
	BXIN;

	d = ((double)i)/m_iResolution*(m_fMaxVal-m_fMinVal)-m_fMinVal;
	m_fSum += d;
	m_fSqSum += d*d;
	if (d < m_fMinInput)
		m_fMinInput = d;
	if (d > m_fMaxInput)
		m_fMaxInput = d;

	if ((i < 0) || (i >= m_iResolution))
	{
		m_fSkipEntries++;
//		mprintf("CDF: %.3f out of Range (%.3f-%.3f).\n",d,m_fMinVal,m_fMaxVal);
		BXOUT;
		return;
	}
	m_fBinEntries++;
	m_pBin[i]+=j;
	BXOUT;
}


void CDF::AddToBin_Count(int i, int count)
{
	double d;
	BXIN;

	d = ((double)i)/m_iResolution*(m_fMaxVal-m_fMinVal)-m_fMinVal;
	m_fSum += d*count;
	m_fSqSum += d*d*count;
	if (d < m_fMinInput)
		m_fMinInput = d;
	if (d > m_fMaxInput)
		m_fMaxInput = d;

	if ((i < 0) || (i >= m_iResolution))
	{
		m_fSkipEntries+=count;
//		mprintf("CDF: %.3f out of Range (%.3f-%.3f).\n",d,m_fMinVal,m_fMaxVal);
		BXOUT;
		return;
	}
	m_fBinEntries+=count;
	m_pBin[i]+=count;
	BXOUT;
}


void CDF::AngleCorrect()
{
	BTIN;
	int z;

	for (z=0;z<m_iResolution;z++)
		m_pBin[z] /= cos((m_fMinVal+(double)z*(m_fMaxVal-m_fMinVal)/m_iResolution)*Pi/180.0) - cos((m_fMinVal+(double)(z+1)*(m_fMaxVal-m_fMinVal)/m_iResolution)*Pi/180.0);
	BTOUT;
}


void CDF::MultiplyBin(double f)
{ 
	BTIN;
	int z, z2;

	for (z=0;z<m_iResolution;z++)
		m_pBin[z] *= f;

	if (m_iMultiCount != 0)
	{
		for (z2=0;z2<m_iMultiCount;z2++)
			for (z=0;z<m_iResolution;z++)
				m_pMultiBin[z2][z] *= f;
	}

	for (z2=0;z2<((m_iAdditionalSets<4)?m_iAdditionalSets:4);z2++)
		for (z=0;z<m_iResolution;z++)
			m_pAdditionalSets[z2][z] *= f;

	BTOUT;
}


void CDF::SubtractBin(double f)
{ 
	BTIN;
	int z/*, z2*/;

	for (z=0;z<m_iResolution;z++)
		m_pBin[z] -= f;

/*	if (m_iMultiCount != 0)
	{
		for (z2=0;z2<m_iMultiCount;z2++)
			for (z=0;z<m_iResolution;z++)
				m_pMultiBin[z2][z] *= f;
	}

	for (z2=0;z2<((m_iAdditionalSets<4)?m_iAdditionalSets:4);z2++)
		for (z=0;z<m_iResolution;z++)
			m_pAdditionalSets[z2][z] *= f;*/

	BTOUT;
}


void CDF::MultiplyIntegral(double f)
{ 
	BTIN;
	int z;

	if (m_pIntegral == NULL)
		return;

	for (z=0;z<m_iResolution;z++)
		m_pIntegral[z] *= f;

	BTOUT;
}


double CDF::NormBinIntegral()
{
	BTIN;
	int z, z2;
	double f;

	if (m_fBinEntries == 0)
		return 0;

	f = 0;
	for (z=0;z<m_iResolution;z++)
		f += m_pBin[z];
	f /= 100.0*m_iResolution;
	for (z=0;z<m_iResolution;z++)
		m_pBin[z] /= f;

	if (m_iMultiCount != 0)
	{
		f = 0;
		for (z2=0;z2<m_iMultiCount;z2++)
			for (z=0;z<m_iResolution;z++)
				f += m_pMultiBin[z2][z];
		f /= 100.0*m_iResolution;
		for (z2=0;z2<m_iMultiCount;z2++)
			for (z=0;z<m_iResolution;z++)
				m_pMultiBin[z2][z] /= f;
	}
	return f;
	BTOUT;
}


void CDF::NormBinIntegral(double val)
{
	BTIN;
	int z, z2;
	double f;

	if (m_fBinEntries == 0)
		return;
	f = 0;
	for (z=0;z<m_iResolution;z++)
		f += m_pBin[z];
	f *= 1.0 / m_iResolution * (m_fMaxVal-m_fMinVal) / val;
	for (z=0;z<m_iResolution;z++)
		m_pBin[z] /= f;

	if (m_iMultiCount != 0)
	{
		f = 0;
		for (z2=0;z2<m_iMultiCount;z2++)
			for (z=0;z<m_iResolution;z++)
				f += m_pMultiBin[z2][z];
		f *= val / m_iResolution * (m_fMaxVal-m_fMinVal);
		for (z2=0;z2<m_iMultiCount;z2++)
			for (z=0;z<m_iResolution;z++)
				m_pMultiBin[z2][z] /= f;
	}
	BTOUT;
}


void CDF::NormBinSum(double val)
{
	BTIN;
	int z;
	double f;

	if (m_fBinEntries == 0)
		return;
	f = 0;
	for (z=0;z<m_iResolution;z++)
		f += m_pBin[z];
	for (z=0;z<m_iResolution;z++)
		m_pBin[z] *= val / f;

	BTOUT;
}


void CDF::CorrectRadialDist()
{
	BTIN;
	int z, z2;

	for (z=0;z<m_iResolution;z++)
		m_pBin[z] /= pow(m_fMinVal+(z+1.0)/m_iResolution*(m_fMaxVal-m_fMinVal),3) - pow(m_fMinVal+((double)z)/m_iResolution*(m_fMaxVal-m_fMinVal),3);

	if (m_iMultiCount != 0)
	{
		for (z2=0;z2<m_iMultiCount;z2++)
			for (z=0;z<m_iResolution;z++)
				m_pMultiBin[z2][z] /= pow(m_fMinVal+(z+1.0)/m_iResolution*(m_fMaxVal-m_fMinVal),3) - pow(m_fMinVal+((double)z)/m_iResolution*(m_fMaxVal-m_fMinVal),3);
	}

	if (m_iAdditionalSets != 0)
	{
		for (z2=0;z2<((m_iAdditionalSets<4)?m_iAdditionalSets:4);z2++)
			for (z=0;z<m_iResolution;z++)
//			{
//				mprintf("z2=%d, z=%d, v=%G, d=%G\n",z2,z,m_pAdditionalSets[z2][z],pow(m_fMinVal+(z+1.0)/m_iResolution*(m_fMaxVal-m_fMinVal),3) - pow(m_fMinVal+((double)z)/m_iResolution*(m_fMaxVal-m_fMinVal),3));
				m_pAdditionalSets[z2][z] /= pow(m_fMinVal+(z+1.0)/m_iResolution*(m_fMaxVal-m_fMinVal),3) - pow(m_fMinVal+((double)z)/m_iResolution*(m_fMaxVal-m_fMinVal),3);
//			}
	}

	BTOUT;
}


double LongMode_F1(double r)
{
	return atan( sqrt(4.0*r*r - 2.0) );
}


double LongMode_F2(double r)
{
	return 8.0*r*atan( 2.0*r*(4.0*r*r-3.0) / ( sqrt(4.0*r*r-2.0) * (4.0*r*r + 1.0) ) );
}


double LongMode_P(double r)
{
	if (2.0*r <= 1.0)
	{
//		return 4.0*Pi*r*r;
		return r*r;
	} else if (2.0*r <= sqrt(2.0))
	{
//		return 2.0*Pi*r*(3.0-4.0*r);
		return 0.5*r*(3.0-4.0*r);
	} else if (2.0*r <= sqrt(3.0))
	{
//		return 2.0*r * ( 3.0*Pi - 12.0*LongMode_F1(r) + LongMode_F2(r) );
		return 0.5/Pi*r * ( 3.0*Pi - 12.0*LongMode_F1(r) + LongMode_F2(r) );
	} else return 1e100;
}


void CDF::CorrectRadialDistLong()
{
	BTIN;
	int z;

	for (z=0;z<m_iResolution;z++)
	{
//		mprintf("Dist %.6G pm, Val %.6G pm, Fac %.6G pm\n",(m_fMinVal+(z+0.5)/m_iResolution*(m_fMaxVal-m_fMinVal)),(m_fMinVal+(z+0.5)/m_iResolution*(m_fMaxVal-m_fMinVal)) / g_fBoxX,LongMode_P( (m_fMinVal+(z+0.5)/m_iResolution*(m_fMaxVal-m_fMinVal)) / g_fBoxX ));
		m_pBin[z] /= 3*((m_fMaxVal-m_fMinVal)/m_iResolution) * LongMode_P( (m_fMinVal+(z+0.5)/m_iResolution*(m_fMaxVal-m_fMinVal)) / g_fBoxX ) * g_fBoxX * g_fBoxX;
	}

	BTOUT;
}


void CDF::CorrectLiRadialDist()
{
	BTIN;
	int z, z2;

	for (z=0;z<m_iResolution;z++)
		m_pBin[z] /= pow(m_fMinVal+(z+1.0)/m_iResolution*(m_fMaxVal-m_fMinVal),2) - pow(m_fMinVal+((double)z)/m_iResolution*(m_fMaxVal-m_fMinVal),2);

	if (m_iMultiCount != 0)
	{
		for (z2=0;z2<m_iMultiCount;z2++)
			for (z=0;z<m_iResolution;z++)
				m_pMultiBin[z2][z] /= pow(m_fMinVal+(z+1.0)/m_iResolution*(m_fMaxVal-m_fMinVal),2) - pow(m_fMinVal+((double)z)/m_iResolution*(m_fMaxVal-m_fMinVal),2);
	}
	BTOUT;
}


double CDF::NormalizeBin(double mi, double ma)
{
	BTIN;
	int z;
	double tmi, tma, d, td;

	if (m_fBinEntries == 0)
		return 0;
	tmi = 99999999.0f;
	tma = 0.0f;
	for (z=0;z<m_iResolution;z++)
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
	for (z=0;z<m_iResolution;z++)
		m_pBin[z] = ((m_pBin[z]-tmi)/td*d)+mi;
	BTOUT;
	return tma;
}


double CDF::GetPercentageRange(double perc)
{
	BTIN;
	int z;
	double s, s2, l;

	if (perc >= 1.0f)
		return m_fMaxVal;
	s = 0;
	for (z=0;z<m_iResolution;z++)
		s += m_pBin[z];
	l = s*perc;
	s2 = 0;
	for (z=0;z<m_iResolution;z++)
	{
		s2 += m_pBin[z];
		if (s2 >= l)
			return m_fMinVal+z*(m_fMaxVal-m_fMinVal)/(double)m_iResolution;
	}
	BTOUT;
	return m_fMaxVal;
}


void CDF::Create()
{
	BTIN;
	int z;
//	m_bRDF = rdf;

	try { m_pBin = new double[m_iResolution]; } catch(...) { m_pBin = NULL; }
	if (m_pBin == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iResolution;z++)
		m_pBin[z] = 0;
	m_fFac = m_iResolution/(m_fMaxVal-m_fMinVal);
//	m_pDHistogram = NULL;
	BTOUT;
}


void CDF::CreateMulti(int n)
{
	BTIN;
	int z, z2;
	
	m_iMultiCount = n;
	m_fBinEntries = 0;
	m_fSkipEntries = 0;

	try { m_pBin = new double[m_iResolution]; } catch(...) { m_pBin = NULL; }
	if (m_pBin == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<m_iResolution;z++)
		m_pBin[z] = 0;

	try { m_sLabelMulti = new char*[n]; } catch(...) { m_sLabelMulti = NULL; }
	if (m_sLabelMulti == NULL) NewException((double)n*sizeof(char*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<n;z++)
		m_sLabelMulti[z] = NULL;

	try { m_pMultiBin = new double*[m_iMultiCount]; } catch(...) { m_pMultiBin = NULL; }
	if (m_pMultiBin == NULL) NewException((double)m_iMultiCount*sizeof(double*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iMultiCount;z++)
	{
		try { m_pMultiBin[z] = new double[m_iResolution]; } catch(...) { m_pMultiBin[z] = NULL; }
		if (m_pMultiBin[z] == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z2=0;z2<m_iResolution;z2++)
			m_pMultiBin[z][z2] = 0;
	}
	m_pHistogram = NULL;
	m_fFac = m_iResolution/(m_fMaxVal-m_fMinVal);
	m_fSum = 0;
	m_fSqSum = 0;
	m_fMean = 0;
	m_fSD = 0;
	m_fMinInput = 1e50;
	m_fMaxInput = -1e50;
	BTOUT;
}


void CDF::Write(const char *prefix, const char *name, const char *suffix, bool integral)
{
	BTIN;
	FILE *a;
	int z, z2;
//	char buf[32768];
	CxString buf;
	double /*d,*/ x;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);

	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	if (m_sLabelX != NULL)
		mfprintf(a,"# %s;  ",m_sLabelX);
			else mfprintf(a,"# (no label);  ");
	if (m_sLabelY != NULL)
		mfprintf(a,"%s",m_sLabelY);
			else mfprintf(a,"(no label)");
	if (integral && (m_pIntegral != NULL))
		mfprintf(a,";  Integral");
	for (z2=0;z2<m_iAdditionalSets;z2++)
	{
		if (m_pAdditionalSetLabels[z2] != NULL)
			mfprintf(a,";  %s",m_pAdditionalSetLabels[z2]);
	}
	mfprintf(a,"\n");
//	d = 0;
	for (z=0;z<m_iResolution;z++)
	{
		if (m_bLeft)
			x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;
//		d += m_pBin[z];
		if (integral && (m_pIntegral != NULL))
			mfprintf(a,"%#.10G;  %#.10G;  %#.10G",x,m_pBin[z],m_pIntegral[z]);
				else mfprintf(a,"%#.10G;  %#.10G",x,m_pBin[z]);
		for (z2=0;z2<m_iAdditionalSets;z2++)
		{
			if (m_pAdditionalSets[z2] != NULL)
				mfprintf(a,"; %#.10G",m_pAdditionalSets[z2][z]);
		}
		mfprintf(a,"\n");
	}
	fclose(a);
	BTOUT;
}


void CDF::WriteMulti(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z, z2;
//	char buf[32768];
	CxString buf;
	double x;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	if (m_sLabelX != NULL)
		mfprintf(a,"# %s",m_sLabelX);
			else mfprintf(a,"# (no X label)");

	if (m_sLabelY != NULL)
		mfprintf(a,";  %s",m_sLabelY);
			else mfprintf(a,";  (no Y label)");

	for (z=0;z<m_iMultiCount;z++)
	{
		if (m_sLabelMulti != NULL)
		{
			if (m_sLabelMulti[z] != NULL)
				mfprintf(a,";  %s",m_sLabelMulti[z]);
					else mfprintf(a,";  (no label)");
		} else mfprintf(a,";  (no label)");
	}

	mfprintf(a,"\n");

	for (z=0;z<m_iResolution;z++)
	{
		if (m_bLeft)
			x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;
		mfprintf(a,"%#.10G",x);

		mfprintf(a,";  %#.10G",m_pBin[z]);

		for (z2=0;z2<m_iMultiCount;z2++)
			mfprintf(a,";  %#.10G",m_pMultiBin[z2][z]);

		mfprintf(a,"\n");
	}
	fclose(a);
	BTOUT;
}


void CDF::WriteMulti_Cumulative(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z, z2;
//	char buf[32768];
	CxString buf;
	double d, x;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);

	if (m_sLabelX != NULL)
		mfprintf(a,"# %s",m_sLabelX);
			else mfprintf(a,"# (no X label)");

	if (m_sLabelY != NULL)
		mfprintf(a,";  %s",m_sLabelY);
			else mfprintf(a,";  (no Y label)");

	for (z=0;z<m_iMultiCount;z++)
	{
		if (m_sLabelMulti != NULL)
		{
			if (m_sLabelMulti[z] != NULL)
				mfprintf(a,";  %s",m_sLabelMulti[z]);
					else mfprintf(a,";  (no label)");
		} else mfprintf(a,";  (no label)");
	}

	mfprintf(a,"\n");

	for (z=0;z<m_iResolution;z++)
	{
		if (m_bLeft)
			x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;
		mfprintf(a,"%#.10G",x);

		mfprintf(a,";  %#.10G",m_pBin[z]);

		d = 0;
		for (z2=0;z2<m_iMultiCount;z2++)
		{
			d += m_pMultiBin[z2][z];
			mfprintf(a,";  %#.10G",d);
		}

		mfprintf(a,"\n");
	}
	fclose(a);
	BTOUT;
}


void CDF::WriteMultiAgr(const char *prefix, const char *name, const char *suffix, const char *title, bool rdf)
{
	CGrace *g;
	int z0, z;
//	char buf[32768];
	CxString buf;
	double x, tfs;

	tfs = 1.0;

//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	try { g = new CGrace(); } catch(...) { g = NULL; }
	if (g == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	g->CurrentGraph()->m_bLegend = true;
	g->AddDataset();

	CalcMinMax();

	if (m_sLabelX != NULL)
		g->SetLabelX(m_sLabelX);
	if (m_sLabelY != NULL)
		g->SetLabelY(m_sLabelY);
	g->SetTitle(title);
	g->SetRangeX(m_fMinVal,m_fMaxVal);
	g->SetRangeY((m_fMinEntry<0.0)?m_fMinEntry:0.0,m_fMaxEntry+(m_fMaxEntry-((m_fMinEntry<0.0)?m_fMinEntry:0.0))*0.1);
	g->MakeTicks();
	g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iResolution*2);
	g->SetDatasetName("Total");
	if (rdf)
		g->AddLine(m_fMinVal,1.0,m_fMaxVal,1.0,1,1,0,0,0);

	if (m_iResolution > 2000)
	{
		mprintf("      Preparing set  1: [");
		tfs = m_iResolution / 50.0;
	}

	for (z=0;z<m_iResolution;z++)
	{
		if (m_bLeft)
			x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;

		g->AddXYTupel(0,x,m_pBin[z]);

		if (m_iResolution > 2000)
			if (fmod(z,tfs) < 1)
				mprintf(WHITE,"#");
	}

	if (m_iResolution > 2000)
		mprintf("]\n");

	for (z0=0;z0<m_iMultiCount;z0++)
	{
		if (m_iResolution > 2000)
			mprintf("      Preparing set %2d: [",z0+2);

		g->AddDataset();
		g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iResolution*2);
		if (m_sLabelMulti != NULL)
			if (m_sLabelMulti[z0] != NULL)
				g->SetDatasetName(m_sLabelMulti[z0]);

		for (z=0;z<m_iResolution;z++)
		{
			if (m_bLeft)
				x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
					else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;
			g->AddXYTupel(z0+1,x,m_pMultiBin[z0][z]);

			if (m_iResolution > 2000)
				if (fmod(z,tfs) < 1)
					mprintf(WHITE,"#");
		}

		if (m_iResolution > 2000)
			mprintf("]\n");
	}
	g->WriteAgr(buf,false);
	delete g;
}


void CDF::WriteMultiAgr_Cumulative(const char *prefix, const char *name, const char *suffix, const char *title, bool rdf)
{
	CGrace *g;
	int z0, z;
//	char buf[32768];
	CxString buf;
	double x, tfs;
	double *td;

//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	try { g = new CGrace(); } catch(...) { g = NULL; }
	if (g == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	g->CurrentGraph()->m_bLegend = true;
	g->CurrentGraph()->m_bInvert = true;
//	g->AddDataset();

	CalcMinMax();

	if (m_sLabelX != NULL)
		g->SetLabelX(m_sLabelX);
	if (m_sLabelY != NULL)
		g->SetLabelY(m_sLabelY);
	g->SetTitle(title);
	g->SetRangeX(m_fMinVal,m_fMaxVal);
	g->SetRangeY((m_fMinEntry<0.0)?m_fMinEntry:0.0,m_fMaxEntry+(m_fMaxEntry-((m_fMinEntry<0.0)?m_fMinEntry:0.0))*0.1);
	g->MakeTicks();
	if (rdf)
		g->AddLine(m_fMinVal,1.0,m_fMaxVal,1.0,1,1,0,0,0);
/*	g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iResolution*2);
	g->SetDatasetName("Total");
	mprintf("      Preparing set  1: [");
	for (z=0;z<m_iResolution;z++)
	{
		if (m_bLeft)
			x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;

		g->AddXYTupel(0,x,m_pBin[z]);

		if (fmod(z,tfs) < 1)
			mprintf(WHITE,"#");
	}
	mprintf("]\n");*/
	tfs = m_iResolution / 50.0;

	try { td = new double[m_iResolution]; } catch(...) { td = NULL; }
	if (td == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iResolution;z++)
		td[z] = 0;
	for (z0=0;z0<m_iMultiCount;z0++)
	{
		g->AddDataset();
		g->CurrentGraph()->CurrentDataset()->m_bFill = true;
//		g->CurrentGraph()->CurrentDataset()->m_fLineWidth = 0;
		g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iResolution*2);
	}
	for (z0=0;z0<m_iMultiCount;z0++)
	{
		if (m_iResolution > 2000)
			mprintf("      Preparing set %2d: [",z0+1);

		if (m_sLabelMulti != NULL)
			if (m_sLabelMulti[z0] != NULL)
				g->SetDatasetName(m_iMultiCount-z0-1,m_sLabelMulti[z0]);
		for (z=0;z<m_iResolution;z++)
		{
			if (m_bLeft)
				x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
					else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;

			td[z] += m_pMultiBin[z0][z];
			g->AddXYTupel(m_iMultiCount-z0-1,x,td[z]);

			if (m_iResolution > 2000)
				if (fmod(z,tfs) < 1)
					mprintf(WHITE,"#");
		}

		if (m_iResolution > 2000)
			mprintf("]\n");
	}
	g->WriteAgr(buf,false);
	delete g;
	delete[] td;
}


void CDF::Write_Int(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z;
//	char buf[32768];
	CxString buf;
	double d;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	if (m_sLabelX != NULL)
		mfprintf(a,"# %s;  ",m_sLabelX);
			else mfprintf(a,"# (no label);  ");
	if (m_sLabelY != NULL)
		mfprintf(a,"%s",m_sLabelY);
			else mfprintf(a,"(no label)");
	mfprintf(a,";  Cumulative Sum\n");
	d = 0;
	for (z=0;z<m_iResolution;z++)
	{
		d += m_pBin[z];
		mfprintf(a,"%d;  %#.10G;  %#.10G\n",z,m_pBin[z],d);
	}
	fclose(a);
	BTOUT;
}


void CDF::WriteAgr(const char *prefix, const char *name, const char *suffix, const char *title, bool rdf)
{
	CGrace *g;
	int z0, z, zi;
//	char buf[32768];
	CxString buf;
	double x, tfs;

	tfs = 1.0;

//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	try { g = new CGrace(); } catch(...) { g = NULL; }
	if (g == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	g->AddDataset();

	CalcMinMax();

	if (m_sLabelX != NULL)
		g->SetLabelX(m_sLabelX);
	if (m_sLabelY != NULL)
		g->SetLabelY(m_sLabelY);

	g->SetTitle(title);
	g->SetRangeX(m_fMinVal,m_fMaxVal);
	g->SetRangeY((m_fMinEntry<0.0)?m_fMinEntry:0.0,m_fMaxEntry+(m_fMaxEntry-((m_fMinEntry<0.0)?m_fMinEntry:0.0))*0.1);

	g->MakeTicks();
	g->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iResolution*2);
	if (rdf)
		g->AddLine(m_fMinVal,1.0,m_fMaxVal,1.0,1,1,0,0,0);

	if (m_iResolution > 2000)
	{
		mprintf("      Preparing set  1: [");
		tfs = m_iResolution / 50.0;
	}

	for (z=0;z<m_iResolution;z++)
	{
		if (m_bLeft)
			x = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;

		g->AddXYTupel(0,x,m_pBin[z]);

		if (m_iResolution > 2000)
			if (fmod(z,tfs) < 1)
				mprintf(WHITE,"#");
	}

	if (m_iResolution > 2000)
		mprintf("]\n");

	zi = 0;
	for (z0=0;z0<m_iAdditionalSets;z0++)
	{
		if (m_pAdditionalSets[z0] == NULL)
			continue;

		if (m_iResolution > 2000)
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

			if (m_iResolution > 2000)
				if (fmod(z,tfs) < 1)
					mprintf(WHITE,"#");
		}
		if (m_iResolution > 2000)
			mprintf("]\n");
		zi++;
	}
	g->WriteAgr(buf,false);
	delete g;
}


void CDF::WriteHenry(const char *prefix, const char *name, const char *suffix)
{
	BTIN;
	FILE *a;
	int z;
//	char buf[32768];
	CxString buf;
	double d;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	d = 0;
	for (z=0;z<m_iResolution-1;z+=10)
	{
		d += m_pBin[z];
		mfprintf(a,"%9.4f;  %f;  %f\n",m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution,m_pBin[z],d);
	}
	fclose(a);
	BTOUT;
}


void CDF::Integrate(bool correctradial, double fac)
{
	int z;
	double d;

	try { m_pIntegral = new double[m_iResolution]; } catch(...) { m_pIntegral = NULL; }
	if (m_pIntegral == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iResolution;z++)
	{
		if (correctradial)
			d = m_pBin[z] * (pow(m_fMinVal+(z+1.0)/m_iResolution*(m_fMaxVal-m_fMinVal),3) - pow(m_fMinVal+((double)z)/m_iResolution*(m_fMaxVal-m_fMinVal),3));
				else d = m_pBin[z];
		if (z == 0)
			m_pIntegral[z] = d * fac;
				else m_pIntegral[z] = m_pIntegral[z-1] + d * fac;
	}
}


void CDF::PrepareAdapt()
{
	try { m_pBinTree = new CBinTree(); } catch(...) { m_pBinTree = NULL; }
	if (m_pBinTree == NULL) NewException((double)sizeof(CBinTree),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pBinTree->m_pParent = NULL;
	m_pBinTree->m_pValue = NULL;
	REC_FillBinTree(0,0,m_pBinTree);
	REC_FuseTree(m_pBinTree);
}


unsigned long CDF::REC_FillBinTree(int pos, int depth, CBinTree *p)
{
	CBinTree *t;
	double *d;

	if (depth < 16)
	{
		try { t = new CBinTree(); } catch(...) { t = NULL; }
		if (t == NULL) NewException((double)sizeof(CBinTree),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		t->m_pValue = NULL;
		t->m_pParent = p;
		p->m_pChildren[0] = t;
		pos = REC_FillBinTree(pos,depth+1,t);

		try { t = new CBinTree(); } catch(...) { t = NULL; }
		if (t == NULL) NewException((double)sizeof(CBinTree),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		t->m_pValue = NULL;
		t->m_pParent = p;
		p->m_pChildren[1] = t;
		pos = REC_FillBinTree(pos,depth+1,t);
		return pos;
	} else
	{
		try { d = new double; } catch(...) { d = NULL; }
		if (d == NULL) NewException((double)sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		*d = m_pBin[pos];
		p->m_pValue = (void*)d;
//		mprintf("FillTree: Pos=%d\n",pos);
		return pos+1;
	}
}


void CDF::REC_FuseTree(CBinTree *p)
{
	double *d, *d1, *d2;

	if (p->m_pChildren[0]->m_pValue == NULL)
		REC_FuseTree(p->m_pChildren[0]);

	if (p->m_pChildren[1]->m_pValue == NULL)
		REC_FuseTree(p->m_pChildren[1]);

	d1 = (double*)p->m_pChildren[0]->m_pValue;
	d2 = (double*)p->m_pChildren[1]->m_pValue;

	try { d = new double; } catch(...) { d = NULL; }
	if (d == NULL) NewException((double)sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*d = *d1 + *d2;
	p->m_pValue = (void*)d;
}


void CDF::REC_SaveTree(FILE *a, double xmin, double xmax, int depth, int mindepth, int maxdepth, double thres, CBinTree *p, bool rdf)
{
//	int z;
	double /*thr,*/ d1, d2, v;

	if (depth == maxdepth)
	{
		v = *(double*)p->m_pValue;
		mfprintf(a,"%9.4f;  %f\n",(xmax+xmin)/2.0,v);
		return;
	}

	d1 = *((double*)p->m_pChildren[0]->m_pValue);
	d2 = *((double*)p->m_pChildren[1]->m_pValue);

	if ((sqrt(pow(thres,2)-pow(d1-d2,2)) < (xmax-xmin)) && (depth >= mindepth))
	{
		v = *(double*)p->m_pValue;
		mfprintf(a,"%9.4f;  %f\n",(xmax+xmin)/2.0,v);
	} else
	{
		REC_SaveTree(a,xmin,xmin+(xmax-xmin)/2.0,depth+1,mindepth,maxdepth,thres,p->m_pChildren[0],rdf);
		REC_SaveTree(a,xmin+(xmax-xmin)/2.0,xmax,depth+1,mindepth,maxdepth,thres,p->m_pChildren[1],rdf);
	}
}


void CDF::WriteAdapted(const char *prefix, const char *name, const char *suffix, int mindepth, int maxdepth, double thres, bool rdf)
{
	BTIN;
	FILE *a;
//	char buf[32768];
	CxString buf;
	
//	strcpy(buf,prefix);
//	strcat(buf,name);
//	strcat(buf,suffix);
	buf.strcpy(prefix);
	buf.strcat(name);
	buf.strcat(suffix);

	a = OpenFileWrite(buf,true);
	REC_SaveTree(a,m_fMinVal,m_fMaxVal,0,mindepth,maxdepth,thres,m_pBinTree,rdf);
	fclose(a);
	BTOUT;
}


void CDF::BinTree_RadialDist()
{
	REC_BinTreeRadialDist(m_fMinVal,m_fMaxVal,0,m_pBinTree);
}


void CDF::REC_BinTreeRadialDist(double xmin, double xmax, int depth, CBinTree *p)
{
	*(double*)p->m_pValue /= pow(xmax,3)-pow(xmin,3);

	if (depth == 16)
		return;

	REC_BinTreeRadialDist(xmin,xmin+(xmax-xmin)/2.0,depth+1,p->m_pChildren[0]);
	REC_BinTreeRadialDist(xmin+(xmax-xmin)/2.0,xmax,depth+1,p->m_pChildren[1]);
}


void CDF::BinTree_MultiplyBin(double f)
{
	REC_BinTreeMultiplyBin(m_fMinVal,m_fMaxVal,0,m_pBinTree,f);
}


void CDF::REC_BinTreeMultiplyBin(double xmin, double xmax, int depth, CBinTree *p, double fac)
{
	*(double*)p->m_pValue *= fac;

	if (depth == 16)
		return;

	REC_BinTreeMultiplyBin(xmin,xmin+(xmax-xmin)/2.0,depth+1,p->m_pChildren[0],fac);
	REC_BinTreeMultiplyBin(xmin+(xmax-xmin)/2.0,xmax,depth+1,p->m_pChildren[1],fac);
}


void CDF::CreateCombinedPlot(bool rdf)
{
	double x, x2, imax, tmx, px;
	int z/*, tpx*/, mi, ma;
//	char buf[64];
	CxString buf;

	m_pCombinedPlot->FindMinMaxVal();
	x = m_pCombinedPlot->CurrentGraph()->m_fMaxValX;
	x2 = 0.25*x;
	m_pCombinedPlot->SetRangeX(0,x+x2+0.025*x);
	m_pCombinedPlot->SetRangeY(m_fMinVal,m_fMaxVal);

	CreateTicks(m_fMinVal,m_fMaxVal,ma,mi);
	m_pCombinedPlot->CurrentGraph()->m_fTickMajorY = (m_fMaxVal-m_fMinVal)/(ma-1);
	m_pCombinedPlot->CurrentGraph()->m_iTickMinorY = mi;
	m_pCombinedPlot->CurrentGraph()->m_iTickPrecY = (int)max(0,1-int(ceil(log10(m_pCombinedPlot->CurrentGraph()->m_fTickMajorY))));

	m_pCombinedPlot->AddDataset();
	m_pCombinedPlot->CurrentGraph()->CurrentDataset()->m_iLineColorIndex = 1;

//	mprintf("A\n");

	imax = 0;
	for (z=0;z<m_iResolution;z++)
		if (m_pBin[z] > imax)
			imax = m_pBin[z];

	// Die Trennlinie zwischen TD und Histogramm
	m_pCombinedPlot->AddLine(x,m_fMinVal,x,m_fMaxVal,2,1);

	// Die Nulllinie
//	m_pCombinedPlot->AddLine(x+x2+0.025*x,m_fMinVal,x+x2+0.025*x,m_fMaxVal,2,1);

	// Die "Eins"-Linie
	if (rdf)
		m_pCombinedPlot->AddLine(x+x2+0.025*x-(1.0/imax*x2),m_fMinVal,x+x2+0.025*x-(1.0/imax*x2),m_fMaxVal,1,1);

	m_pCombinedPlot->CurrentGraph()->m_bTicksBothSidesY = true;
	m_pCombinedPlot->CurrentGraph()->m_bTickLabelsBothSidesY = true;

	x += 0.025*x;

	m_pCombinedPlot->CurrentGraph()->CurrentDataset()->m_faValues.SetMaxSize(m_iResolution*2);
	for (z=0;z<m_iResolution-1;z++)
		m_pCombinedPlot->AddXYTupel(m_pCombinedPlot->CurrentGraph()->m_oaDatasets.GetSize()-1,x+x2-(m_pBin[z]/imax*x2),m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution);

	tmx = majorticks(0,(float)x);
//	tpx = (int)max(0,1-int(ceil(log10(tmx))));

//	mprintf("B\n");

//	mprintf("x=%f, tmx=%f, tpx=%f.\n",x,tmx,tpx);

	px = 0;
	z = 0;
	do {
		if ((z%2)==0)
		{
//			sprintf(buf,"%.1f ps",px);
			buf.sprintf("%.1f ps",px);
			m_pCombinedPlot->AddCustomLabelX(true,px,buf);
		} else
			m_pCombinedPlot->AddCustomLabelX(false,px,"");
		z++;
		px += tmx/2.0;
	} while (px+tmx/2.0 < x);
//	mprintf("C\n");

	m_pCombinedPlot->AddCustomLabelX(true,x+x2,"0");

	if (rdf)
		m_pCombinedPlot->AddCustomLabelX(true,x+x2-(1.0/imax*x2),"1");
/*	sprintf(buf,"%.1f",imax);
	m_pCombinedPlot->AddCustomLabelX(true,x,buf);*/
	for (z=2;z<imax-1;z+=max(1,(int)(imax/4)))
		m_pCombinedPlot->AddCustomLabelX(false,x+x2-((double)z/imax*x2),"");
//	sprintf(buf,"%.0f",floor(imax));
	buf.sprintf("%.0f",floor(imax));
	m_pCombinedPlot->AddCustomLabelX(true,x+x2-(floor(imax)/imax*x2),buf);
}


void CDF::ScaleXRange(double fac)
{
	m_fMinVal *= fac;
	m_fMaxVal *= fac;
}


void CDF::CalcMinMax()
{
	int z, z2;

	m_fMinEntry = 9E99;
	m_fMaxEntry = -9E99;

	for (z=0;z<m_iResolution;z++)
	{
		if (m_pBin[z] > m_fMaxEntry)
			m_fMaxEntry = m_pBin[z];
		if (m_pBin[z] < m_fMinEntry)
			m_fMinEntry = m_pBin[z];
	}

	if (m_iMultiCount != 0)
	{
		for (z2=0;z2<m_iMultiCount;z2++)
		{
			for (z=0;z<m_iResolution;z++)
			{
				if (m_pMultiBin[z2][z] > m_fMaxEntry)
					m_fMaxEntry = m_pMultiBin[z2][z];
				if (m_pMultiBin[z2][z] < m_fMinEntry)
					m_fMinEntry = m_pMultiBin[z2][z];
			}
		}
	}
}


void CDF::CalcHistogram()
{
	int z, z2, t;

	CalcMinMax();

	try { m_pHistogram = new double[m_iHistogramRes]; } catch(...) { m_pHistogram = NULL; }
	if (m_pHistogram == NULL) NewException((double)m_iHistogramRes*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iHistogramRes;z++)
		m_pHistogram[z] = 0;
	for (z=0;z<m_iResolution;z++)
	{
		t = (int)floor((m_pBin[z]-m_fMinEntry)/(m_fMaxEntry-m_fMinEntry)*m_iHistogramRes);
		for (z2=0;z2<t;z2++)
			m_pHistogram[z2]++;
	}
	for (z=0;z<m_iHistogramRes;z++)
		m_pHistogram[z] /= (double)m_iResolution;
}


void CDF::WriteHistogram(const char *prefix, const char *name, const char *suffix)
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


void CDF::CopyFrom(CDF *p)
{
	m_bCombinedPlot = p->m_bCombinedPlot;
	m_bSaveDist = p->m_bSaveDist;
	m_fBinEntries = p->m_fBinEntries;
	m_fFac = p->m_fFac;
	m_fMaxEntry = p->m_fMaxEntry;
	m_fMaxP = p->m_fMaxP;
	m_fMaxVal = p->m_fMaxVal;
	m_fMinEntry = p->m_fMinEntry;
	m_fMinVal = p->m_fMinVal;
	m_fSkipEntries = p->m_fSkipEntries;
	m_iHistogramRes = p->m_iHistogramRes;
	m_iResolution = p->m_iResolution;
	m_fSum = p->m_fSum;
	m_fSqSum = p->m_fSqSum;
	m_fMean = p->m_fMean;
	m_fSD = p->m_fSD;
	m_fMinInput = p->m_fMinInput;
	m_fMaxInput = p->m_fMaxInput;

	if (p->m_pBin != NULL)
	{
		try { m_pBin = new double[m_iResolution]; } catch(...) { m_pBin = NULL; }
		if (m_pBin == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_pBin,p->m_pBin,sizeof(double)*m_iResolution);
	}
	m_pBinTree = NULL;
	m_pCombinedPlot = NULL;

	if (p->m_pHistogram != NULL)
	{
		try { m_pHistogram = new double[m_iHistogramRes]; } catch(...) { m_pHistogram = NULL; }
		if (m_pHistogram == NULL) NewException((double)m_iHistogramRes*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_pHistogram,p->m_pHistogram,sizeof(double)*m_iHistogramRes);
	}

	if (p->m_pIntegral != NULL)
	{
		try { m_pIntegral = new double[m_iResolution]; } catch(...) { m_pIntegral = NULL; }
		if (m_pIntegral == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		memcpy(m_pIntegral,p->m_pIntegral,sizeof(double)*m_iResolution);
	}
}


void CDF::Fit_PolyExp(int degree, int maxcall)
{
	double *x, *y;
	int z;
//	char buf[256];
	CxString buf;

	try { m_pAdditionalSets[degree] = new double[m_iResolution]; } catch(...) { m_pAdditionalSets[degree] = NULL; }
	if (m_pAdditionalSets[degree] == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

//	sprintf(buf,"%d-exp fit",degree);
	buf.sprintf("%d-exp fit",degree);

	try { m_pAdditionalSetLabels[degree] = new char[strlen(buf)+1]; } catch(...) { m_pAdditionalSetLabels[degree] = NULL; }
	if (m_pAdditionalSetLabels[degree] == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_pAdditionalSetLabels[degree],buf);

	try { m_pLMWrapper = new CLMWrapper(); } catch(...) { m_pLMWrapper = NULL; }
	if (m_pLMWrapper == NULL) NewException((double)sizeof(CLMWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { x = new double[m_iResolution]; } catch(...) { x = NULL; }
	if (x == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { y = new double[m_iResolution]; } catch(...) { y = NULL; }
	if (y == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_pParameters[degree] = new double[degree*2]; } catch(...) { m_pParameters[degree] = NULL; }
	if (m_pParameters[degree] == NULL) NewException((double)degree*2*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<m_iResolution;z++)
	{
		if (m_bLeft)
			x[z] = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x[z] = m_fMinVal+(z+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;
		y[z] = m_pBin[z];
	}
	if (degree == 1)
	{
		m_pParameters[1][0] = 1.0;
		if (g_iTrajSteps != -1)
			m_pParameters[degree][1] = -1.0/((g_fTimestepLength/1000.0)*g_iTrajSteps/10.0);
				else m_pParameters[degree][1] = -1.0/((g_fTimestepLength/1000.0)*1000.0);
	} else
	{
		for (z=0;z<degree;z++)
		{
			// Smallest lifetime: [g_fTimestepLength/1000.0] ps
			// Largest lifetime:  [g_iTrajSteps*g_fTimestepLength/1000.0] ps
			m_pParameters[degree][z*2] = 1.0/degree;
			if (g_iTrajSteps != -1)
				m_pParameters[degree][z*2+1] = -1.0/((g_fTimestepLength/1000.0)*pow(10.0,log10(g_iTrajSteps)*double(z)/(degree-1)));
					else m_pParameters[degree][z*2+1] = -1.0/((g_fTimestepLength/1000.0)*pow(10.0,log10(10000.0)*double(z)/(degree-1)));
		}
	}

	m_pLMWrapper->Fit_PolyExp(degree,m_iResolution,x,y,m_pParameters[degree],m_pCorrCoeff[degree],m_pFitIntegral[degree],m_pAdditionalSets[degree],maxcall);

	delete[] x;
	delete[] y;
	delete m_pLMWrapper;
	m_pLMWrapper = NULL;
}


void CDF::Fit_ExpSpectrum(int res, double mi, double ma, const char *name, int dpoints, int maxcall, bool lindata, double zeroweight, bool evolve)
{
	double *x, *y;
	int z;
//	char buf[256];
	CxString buf;
	int *dpmapping;
	double *val, *valx;

	try { dpmapping = new int[dpoints]; } catch(...) { dpmapping = NULL; }
	if (dpmapping == NULL) NewException((double)dpoints*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { val = new double[m_iResolution]; } catch(...) { val = NULL; }
	if (val == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { valx = new double[m_iResolution]; } catch(...) { valx = NULL; }
	if (valx == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_pAdditionalSets = new double*; } catch(...) { m_pAdditionalSets = NULL; }
	if (m_pAdditionalSets == NULL) NewException((double)sizeof(double*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_pAdditionalSetLabels = new char*; } catch(...) { m_pAdditionalSetLabels = NULL; }
	if (m_pAdditionalSetLabels == NULL) NewException((double)sizeof(char*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_iAdditionalSets = 1;

	try { m_pAdditionalSets[0] = new double[m_iResolution]; } catch(...) { m_pAdditionalSets[0] = NULL; }
	if (m_pAdditionalSets[0] == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

//	sprintf(buf,"Lifetime spectrum fit");
	buf.sprintf("Lifetime spectrum fit");

	try { m_pAdditionalSetLabels[0] = new char[strlen(buf)+1]; } catch(...) { m_pAdditionalSetLabels[0] = NULL; }
	if (m_pAdditionalSetLabels[0] == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_pAdditionalSetLabels[0],buf);

	try { m_pLMWrapper = new CLMWrapper(); } catch(...) { m_pLMWrapper = NULL; }
	if (m_pLMWrapper == NULL) NewException((double)sizeof(CLMWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { x = new double[dpoints]; } catch(...) { x = NULL; }
	if (x == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { y = new double[dpoints]; } catch(...) { y = NULL; }
	if (y == NULL) NewException((double)m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_pParameters = new double*; } catch(...) { m_pParameters = NULL; }
	if (m_pParameters == NULL) NewException((double)sizeof(double*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_pParameters[0] = new double[res]; } catch(...) { m_pParameters[0] = NULL; }
	if (m_pParameters[0] == NULL) NewException((double)res*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	if (lindata)
		mprintf("    Reducing %d data points to %d fit points linearly.\n\n",m_iResolution,dpoints);
	else
		mprintf("    Reducing %d data points to %d fit points logarithmically.\n\n",m_iResolution,dpoints);

	mprintf("    Using the following data points:");

	for (z=0;z<dpoints;z++)
	{
		if ((z%10) == 0)
			mprintf("\n      ");

		if (lindata)
			dpmapping[z] = (int)(z/(dpoints-1.0)*(m_iResolution-1));
				else dpmapping[z] = (int)(exp(z/(dpoints-1.0)*log(m_iResolution))+0.5) - 1;

		if (z != 0)
			if (dpmapping[z] <= dpmapping[z-1])
				dpmapping[z] = dpmapping[z-1]+1;

		mprintf("%6d",dpmapping[z]);
		if (z < dpoints-1)
			mprintf(", ");
				else mprintf(".");
	}
	mprintf("\n\n");

	mprintf("    Using a weight of %.2f for the first data point.\n\n",zeroweight);

	for (z=0;z<dpoints;z++)
	{
		if (m_bLeft)
			x[z] = m_fMinVal+dpmapping[z]*(m_fMaxVal-m_fMinVal)/m_iResolution;
				else x[z] = m_fMinVal+(dpmapping[z]+0.5)*(m_fMaxVal-m_fMinVal)/m_iResolution;
		y[z] = m_pBin[dpmapping[z]];
	}

	for (z=0;z<res;z++)
//		m_pParameters[0][z] = 1.0/res;
		m_pParameters[0][z] = log(1.0/res);//+((rand()%10000)-5000.0)/50000.0;

	for (z=0;z<m_iResolution;z++)
		valx[z] = m_fMinVal+z*(m_fMaxVal-m_fMinVal)/(m_iResolution-1);

	m_pLMWrapper->Fit_ExpSpectrum(res,mi,ma,dpoints,x,y,m_pParameters[0],m_iResolution,valx,val,name,maxcall,zeroweight,evolve);

	for (z=0;z<m_iResolution;z++)
		m_pAdditionalSets[0][z] = val[z];

	delete[] val;
	delete[] valx;
	delete[] dpmapping;
	delete[] x;
	delete[] y;
	delete m_pLMWrapper;
	m_pLMWrapper = NULL;
}


void CDF::CalcMeanSD()
{
	if ((m_fBinEntries+m_fSkipEntries) == 0)
	{
		m_fMinInput = 0;
		m_fMaxInput = 0;
		return;
	}
	m_fMean = m_fSum / (m_fBinEntries+m_fSkipEntries);
	m_fSD = sqrt(m_fSqSum / (m_fBinEntries+m_fSkipEntries) - m_fMean*m_fMean);
}


void CDF::AddFrom(CDF *t)
{
	int z;

	if (m_iResolution != t->m_iResolution)
		abort();
	for (z=0;z<m_iResolution;z++)
		m_pBin[z] += t->m_pBin[z];
	if (t->m_fMinEntry < m_fMinEntry)
		m_fMinEntry = t->m_fMinEntry;
	if (t->m_fMaxEntry > m_fMaxEntry)
		m_fMaxEntry = t->m_fMaxEntry;
	m_fBinEntries += t->m_fBinEntries;
	m_fSkipEntries += t->m_fSkipEntries;
	m_fSum += t->m_fSum;
	m_fSqSum += t->m_fSqSum;
}


void CDF::SetLabelMulti(int n, const char *s)
{
	if (n >= m_iMultiCount)
	{
		eprintf("CDF::SetLabelMulti(): %d >= %d.\n",n,m_iMultiCount);
		abort();
	}

	if (m_sLabelMulti[n] != NULL)
		delete[] m_sLabelMulti[n];

	try { m_sLabelMulti[n] = new char[strlen(s)+1]; } catch(...) { m_sLabelMulti[n] = NULL; }
	if (m_sLabelMulti[n] == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	strcpy(m_sLabelMulti[n],s);
}


void CDF::Mirror(float plane)
{
// 	int z, i, j;
// 	float p, f;

	if ((plane <= m_fMinVal) || (plane >= m_fMaxVal))
		return;
	
	int i, t;
	double *tbin, p;
	try { tbin = new double[m_iResolution]; } catch(...) { tbin = NULL; }
	if (tbin == NULL) NewException((double)m_iResolution * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	p = (plane - m_fMinVal) * m_fFac;
	for (i = 0; i < m_iResolution; i++) {
		tbin[i] = m_pBin[i];
		t = (int)(2.0 * p - i - 1);
		if ((t >= 0) && (t < m_iResolution)) {
			tbin[i] += m_pBin[t];
			tbin[i] /= 2.0;
		}
	}
	
	delete[] m_pBin;
	m_pBin = tbin;
	
// 	// Array index of mirror plane
// 	i = (plane-m_fMinVal)*m_fFac;
// 	mprintf(GREEN, "i: %d %f\n", i, m_fFac);
// 
// 	// Traverse array until mirror plane
// 	for (z=0;z<i;z++)
// 	{
// 		// Current position
// 		p = z/m_fFac + m_fMinVal;
// 
// 		// Index of the mirrored point
// 		j = (-p-m_fMinVal)*m_fFac;
// 		mprintf(GREEN, "%d %d\n", z, j);
// 
// 		if ((j < 0) || (j >= m_iResolution))
// 			continue;
// 
// 		// Average of both points
// 		f = (m_pBin[z]+m_pBin[j]) / 2.0f;
// 
// 		// Write average into both array positions
// 		m_pBin[z] = f;
// 		m_pBin[j] = f;
// 	}
}

