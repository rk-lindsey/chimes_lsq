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

#include "xlongarray.h"


#ifndef DEBUG_ARRAYS
#define m_sName (NULL)
#endif


CxLongArray::CxLongArray()
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::CxLongArray()\n");
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


CxLongArray::CxLongArray(const char *name)
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::CxLongArray(const char *)\n");
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

	
CxLongArray::~CxLongArray()
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::~CxLongArray()\n");
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

	
CxLongArray::CxLongArray(CxLongArray &o) : CxObject()
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::CxLongArray(CxLongArray &)...");
#endif
	unsigned long z;
	m_iSize = o.m_iSize;
	m_iMaxSize = o.m_iMaxSize;
	m_iGrow = o.m_iGrow;

	try { m_pData = new long[m_iMaxSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	for (z=0;z<m_iSize;z++)
		m_pData[z] = o.m_pData[z];
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxLongArray::CopyFrom(CxLongArray *o)
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::CopyFrom(CxLongArray*)...");
#endif
//	unsigned long z;
	m_iSize = o->m_iSize;
	m_iMaxSize = o->m_iMaxSize;
	m_iGrow = o->m_iGrow;
	if (m_pData != NULL)
		delete[] m_pData;

	try { m_pData = new long[m_iMaxSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
//	for (z=0;z<m_iSize;z++)
//		m_pData[z] = o->m_pData[z];
	if (o->m_pData != NULL)
		memcpy(m_pData,o->m_pData,m_iSize*sizeof(long));
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}
	

void CxLongArray::Append(CxLongArray *o)
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::Append(CxLongArray*)...");
#endif
//	unsigned long z;
	if (m_iSize+o->m_iSize > m_iMaxSize)
		SetMaxSize(m_iSize+o->m_iSize);
//	for (z=0;z<o->m_iSize;z++)
//		m_pData[m_iSize+z] = o->m_pData[z];
	if (o->m_pData != NULL)
		memcpy(&m_pData[m_iSize],o->m_pData,o->m_iSize*sizeof(long));
	m_iSize += o->m_iSize;
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}	


void CxLongArray::SetAt(unsigned long pos, long f)
{
	BXIN;
//	unsigned long z;
#ifdef DEBUG_CLONGARRAY
	bool s = false;
	mprintf("@ CxLongArray::Add(long)...");
#endif
	if (pos >= m_iMaxSize)
	{
#ifdef DEBUG_CLONGARRAY
		s = true;
		mprintf("\n");
#endif
		SetMaxSize(pos+1);
	}
	m_pData[pos] = f;
	if (pos >= m_iSize)
//	{
//		for (z=m_iSize;z<pos;z++)
//			m_pData[z] = 0;
		m_iSize = pos+1;
//	}
#ifdef DEBUG_CLONGARRAY
	if (s)
		mprintf("@ done.\n");
	else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxLongArray::Add(long f)
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	bool s = false;
	mprintf("@ CxLongArray::Add(long)...");
#endif
	if (m_iSize+1 > m_iMaxSize)
	{
#ifdef DEBUG_CLONGARRAY
		s = true;
		mprintf("\n");
#endif
		SetMaxSize(m_iMaxSize + m_iGrow);
	}
	m_pData[m_iSize] = f;
	m_iSize++;
#ifdef DEBUG_CLONGARRAY
	if (s)
		mprintf("@ done.\n");
			else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxLongArray::SetSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::SetSize(int): %d...",i);
#endif
	if (m_iSize == i)
		return;
	long *temp;

	try { temp = new long[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		if (m_iSize > i)
			memcpy(temp,m_pData,i*sizeof(long));
		else
			memcpy(temp,m_pData,m_iSize*sizeof(long));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxLongArray::SetMaxSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::SetMaxSize(int): %d...",i);
#endif
	long *temp;

	try { temp = new long[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,m_iSize*sizeof(long));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxLongArray::SetGrow(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::SetGrow(int): %d\n",i);
#endif
	m_iGrow = i;
	BXOUT;
}
		
	
void CxLongArray::RemoveAll()
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::RemoveAll():...");
#endif
	if (m_pData != NULL)
		delete[] m_pData;
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxLongArray::RemoveAll_KeepSize()
{
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::RemoveAll_KeepSize():...");
#endif
	m_iSize = 0;
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxLongArray::RemoveAt(unsigned long pos, unsigned long count)
{
	BXIN;
	long *temp;
		
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::RemoveAt(int, int): %d, %d...",pos,count);
#endif

	try { temp = new long[m_iSize-count]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize-count)*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(long));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(long));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxLongArray::RemoveAt_KeepSize(unsigned long pos, unsigned long count)
{
	BXIN;
	long *temp;
		
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::RemoveAt_KeepSize(int, int): %d, %d...",pos,count);
#endif

	try { temp = new long[m_iMaxSize]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)m_iMaxSize*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(long));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(long));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxLongArray::InsertAt(long f, unsigned long pos)
{
	BXIN;
	long *temp;
		
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::InsertAt(long, int): %d...");
#endif

	try { temp = new long[m_iSize+1]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+1)*sizeof(long),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(long));
		memcpy(&temp[pos+1],&m_pData[pos],(m_iSize-pos)*sizeof(long));
		delete[] m_pData;
	}
	temp[pos] = f;
	m_pData = temp;
	m_iSize++;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxLongArray::SetName(const char *name)
{
	(void)name;
#ifdef DEBUG_ARRAYS
	BXIN;
#ifdef DEBUG_CLONGARRAY
	mprintf("@ CxLongArray::SetName(const char *): \"%s\"...",name);
#endif
	if (m_sName != NULL)
		delete[] m_sName;
	try { m_sName = new char[strlen(name)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,name);
#ifdef DEBUG_CLONGARRAY
	mprintf("done.\n");
#endif
	BXOUT;
#endif
}
