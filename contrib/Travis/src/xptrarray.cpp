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

#include "xptrarray.h"


#ifndef DEBUG_ARRAYS
#define m_sName (NULL)
#endif


CxPtrArray::CxPtrArray()
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::CxPtrArray()\n");
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


CxPtrArray::CxPtrArray(const char *name)
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::CxPtrArray(const char *)\n");
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

	
CxPtrArray::~CxPtrArray()
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::~CxPtrArray()\n");
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

	
CxPtrArray::CxPtrArray(CxPtrArray &o) : CxObject()
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::CxPtrArray(CxPtrArray &)...");
#endif
	unsigned long z;
	m_iSize = o.m_iSize;
	m_iMaxSize = o.m_iMaxSize;
	m_iGrow = o.m_iGrow;

	try { m_pData = new void*[m_iMaxSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(void*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	for (z=0;z<m_iSize;z++)
		m_pData[z] = o.m_pData[z];
#ifdef DEBUG_CPTRARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxPtrArray::SetAt(unsigned long pos, void *o)
{
	BXIN;
	unsigned long z;
#ifdef DEBUG_CPTRARRAY
	bool s = false;
	mprintf("@ CxPtrArray::Add(void *)...");
#endif
	if (pos >= m_iMaxSize)
	{
#ifdef DEBUG_CPTRARRAY
		s = true;
		mprintf("\n");
#endif
		SetMaxSize(pos+1);
	}
	m_pData[pos] = o;
	if (pos >= m_iSize)
	{
		for (z=m_iSize;z<pos;z++)
			m_pData[z] = NULL;
		m_iSize = pos+1;
	}
#ifdef DEBUG_CPTRARRAY
	if (s)
		mprintf("@ done.\n");
	else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxPtrArray::Add(void *o)
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	bool s = false;
	mprintf("@ CxPtrArray::Add(void *)...");
#endif
	if (m_iSize+1 > m_iMaxSize)
	{
#ifdef DEBUG_CPTRARRAY
		s = true;
		mprintf("\n");
#endif
		if (m_iGrow < 1)
			m_iGrow = 1;
//		SetMaxSize(m_iMaxSize + m_iGrow);
		if (m_iMaxSize == 0)
			SetMaxSize(m_iMaxSize + m_iGrow);
				else SetMaxSize(m_iMaxSize*2);
	}
	m_pData[m_iSize] = o;
	m_iSize++;
#ifdef DEBUG_CPTRARRAY
	if (s)
		mprintf("@ done.\n");
			else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxPtrArray::SetSize(unsigned long i)
{
	BXIN;
	int z;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::SetSize(int): %d...",i);
#endif
	if (i == m_iSize)
		return;
	void **temp;

	try { temp = new void*[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(void*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		if (m_iSize > i)
			memcpy(temp,m_pData,i*sizeof(void*));
		else
			memcpy(temp,m_pData,m_iSize*sizeof(void*));
		delete[] m_pData;
	}
// Neu und heikel (Ballmer-Peak ^^)
	for (z=m_iSize;z<(int)i;z++)
		temp[z] = NULL;
// Ende neu
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_CPTRARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxPtrArray::SetMaxSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::SetMaxSize(int): %d...",i);
#endif
	if (i == m_iMaxSize)
	{
		BXOUT;
		return;
	}
	void **temp;

	try { temp = new void*[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(void*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,m_iSize*sizeof(void*));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
#ifdef DEBUG_CPTRARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxPtrArray::SetGrow(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::SetGrow(int): %d\n",i);
#endif
	m_iGrow = i;
	BXOUT;
}
		
	
void CxPtrArray::RemoveAll()
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::RemoveAll():...");
#endif
	if (m_pData != NULL)
		delete[] m_pData;
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
#ifdef DEBUG_CPTRARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxPtrArray::RemoveAll_KeepSize()
{
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::RemoveAll_KeepSize():...");
#endif
	m_iSize = 0;
#ifdef DEBUG_CPTRARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxPtrArray::RemoveAt(unsigned long pos, unsigned long count)
{
	BXIN;
	void **temp;
		
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::RemoveAt(int, int): %d, %d...",pos,count);
#endif

	try { temp = new void*[m_iSize-count]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize-count)*sizeof(void*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
//	for (z=pos;z<pos+count;z++)
//		delete m_pData[z];
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(void*));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(void*));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CPTRARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxPtrArray::InsertAt(void *o, unsigned long pos)
{
	BXIN;
	void **temp;
		
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::InsertAt(void *, int): %d...");
#endif

	try { temp = new void*[m_iSize+1]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+1)*sizeof(void*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(void*));
		memcpy(&temp[pos+1],&m_pData[pos],(m_iSize-pos)*sizeof(void*));
		delete[] m_pData;
	}
	temp[pos] = o;
	m_pData = temp;
	m_iSize++;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CPTRARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxPtrArray::SetName(const char *name)
{
	(void)name;
#ifdef DEBUG_ARRAYS
	BXIN;
#ifdef DEBUG_CPTRARRAY
	mprintf("@ CxPtrArray::SetName(const char *): \"%s\"...",name);
#endif
	if (m_sName != NULL)
		delete[] m_sName;
	try { m_sName = new char[strlen(name)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,name);
#ifdef DEBUG_CPTRARRAY
	mprintf("done.\n");
#endif
	BXOUT;
#endif
}
