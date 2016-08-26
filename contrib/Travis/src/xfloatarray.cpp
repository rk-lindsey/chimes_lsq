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

#include "xfloatarray.h"
#include "tools.h"


#ifndef DEBUG_ARRAYS
#define m_sName (NULL)
#endif


CxFloatArray::CxFloatArray()
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::CxFloatArray()\n");
#endif
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
	m_iGrow = 16;
#ifdef DEBUG_ARRAYS
	m_sName = NULL;
#endif
	BXOUT;
}


CxFloatArray::CxFloatArray(const char *name)
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::CxFloatArray(const char *)\n");
#endif
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
	m_iGrow = 16;
#ifdef DEBUG_ARRAYS
	m_sName = NULL;
#endif
	SetName(name);
	BXOUT;
}

	
CxFloatArray::~CxFloatArray()
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::~CxFloatArray()\n");
#endif
	RemoveAll();
#ifdef DEBUG_ARRAYS
	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}
#endif
	BXOUT;
}

	
CxFloatArray::CxFloatArray(CxFloatArray &o) : CxObject()
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::CxFloatArray(CxFloatArray &)...");
#endif
	unsigned long z;
	m_iSize = o.m_iSize;
	m_iMaxSize = o.m_iMaxSize;
	m_iGrow = o.m_iGrow;

	try { m_pData = new float[m_iMaxSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	for (z=0;z<m_iSize;z++)
		m_pData[z] = o.m_pData[z];
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxFloatArray::CopyFrom(CxFloatArray *o)
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::CopyFrom(CxFloatArray*)...");
#endif
	if (m_iMaxSize != o->m_iMaxSize)
	{
		m_iMaxSize = o->m_iMaxSize;
		if (m_pData != NULL)
			delete[] m_pData;

		try { m_pData = new float[m_iMaxSize]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	}
	m_iGrow = o->m_iGrow;
	m_iSize = o->m_iSize;
	if (o->m_pData != NULL)
		memcpy(m_pData,o->m_pData,m_iSize*sizeof(float));
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}	

	
void CxFloatArray::SetAt(unsigned long pos, float f)
{
	BXIN;
	unsigned long z;
#ifdef DEBUG_CFLOATARRAY
	bool s = false;
	mprintf("@ CxFloatArray::Add(float)...");
#endif
	if (pos >= m_iMaxSize)
	{
#ifdef DEBUG_CFLOATARRAY
		s = true;
		mprintf("\n");
#endif
		SetMaxSize(pos+1);
	}
	m_pData[pos] = f;
	if (pos >= m_iSize)
	{
		for (z=m_iSize;z<pos;z++)
			m_pData[z] = 0.0f;
		m_iSize = pos+1;
	}
#ifdef DEBUG_CFLOATARRAY
	if (s)
		mprintf("@ done.\n");
	else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxFloatArray::Add(float f)
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	bool s = false;
	mprintf("@ CxFloatArray::Add(float)...");
#endif
	if (m_iSize+1 > m_iMaxSize)
	{
#ifdef DEBUG_CFLOATARRAY
		s = true;
		mprintf("\n");
#endif
//		SetMaxSize(m_iMaxSize + m_iGrow);
		if (m_iMaxSize == 0)
			SetMaxSize(m_iMaxSize + m_iGrow);
				else SetMaxSize(m_iMaxSize*2);
	}
	m_pData[m_iSize] = f;
	m_iSize++;
#ifdef DEBUG_CFLOATARRAY
	if (s)
		mprintf("@ done.\n");
			else mprintf("done.\n");
#endif
	BXOUT;
}
	

void CxFloatArray::SetSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::SetSize(int): %d...",i);
#endif
	if (m_iSize == i)
		return;
	float *temp;

	try { temp = new float[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);

	if (m_pData != NULL) {
		if (m_iSize > i)
			memcpy(temp,m_pData,i*sizeof(float));
		else
			memcpy(temp,m_pData,m_iSize*sizeof(float));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxFloatArray::SetMaxSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::SetMaxSize(int): %d...",i);
#endif
	float *temp;

	try { temp = new float[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,m_iSize*sizeof(float));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxFloatArray::SetGrow(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::SetGrow(int): %d\n",i);
#endif
	m_iGrow = i;
	BXOUT;
}
		
	
void CxFloatArray::RemoveAll()
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::RemoveAll():...");
#endif
	if (m_pData != NULL)
		delete[] m_pData;
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxFloatArray::RemoveAll_KeepSize()
{
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::RemoveAll_KeepSize():...");
#endif
	m_iSize = 0;
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxFloatArray::RemoveAt(unsigned long pos, unsigned long count)
{
	BXIN;
	float *temp;
		
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::RemoveAt(int, int): %d, %d...",pos,count);
#endif

	try { temp = new float[m_iSize-count]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize-count)*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(float));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(float));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxFloatArray::InsertAt(float f, unsigned long pos)
{
	BXIN;
	float *temp;
		
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::InsertAt(float, int): %d...");
#endif

	try { temp = new float[m_iSize+1]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+1)*sizeof(float),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(float));
		memcpy(&temp[pos+1],&m_pData[pos],(m_iSize-pos)*sizeof(float));
		delete[] m_pData;
	}
	temp[pos] = f;
	m_pData = temp;
	m_iSize++;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxFloatArray::SetName(const char *name)
{
	(void)name;
#ifdef DEBUG_ARRAYS
	BXIN;
#ifdef DEBUG_CFLOATARRAY
	mprintf("@ CxFloatArray::SetName(const char *): \"%s\"...",name);
#endif
	if (m_sName != NULL)
		delete[] m_sName;
	try { m_sName = new char[strlen(name)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,name);
#ifdef DEBUG_CFLOATARRAY
	mprintf("done.\n");
#endif
	BXOUT;
#endif
}
