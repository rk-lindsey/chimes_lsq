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

#include "xintarray.h"


#ifndef DEBUG_ARRAYS
#define m_sName (NULL)
#endif


CxIntArray::CxIntArray()
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::CxIntArray()\n");
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


CxIntArray::CxIntArray(const char *name)
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::CxIntArray(const char *)\n");
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

	
CxIntArray::~CxIntArray()
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::~CxIntArray()\n");
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

	
CxIntArray::CxIntArray(CxIntArray &o) : CxObject()
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::CxIntArray(CxIntArray &)...");
#endif
	unsigned int z;
	m_iSize = o.m_iSize;
	m_iMaxSize = o.m_iMaxSize;
	m_iGrow = o.m_iGrow;

	try { m_pData = new int[m_iMaxSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	for (z=0;z<m_iSize;z++)
		m_pData[z] = o.m_pData[z];
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxIntArray::CopyFrom(CxIntArray *o)
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::CopyFrom(CxIntArray*)...");
#endif
	m_iSize = o->m_iSize;
	m_iMaxSize = o->m_iMaxSize;
	m_iGrow = o->m_iGrow;
	if (m_pData != NULL)
		delete[] m_pData;

	try { m_pData = new int[m_iMaxSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (o->m_pData != NULL)
		memcpy(m_pData,o->m_pData,m_iSize*sizeof(int));
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}
	

void CxIntArray::Append(CxIntArray *o)
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::Append(CxIntArray*)...");
#endif
//	unsigned int z;
	if (m_iSize+o->m_iSize > m_iMaxSize)
		SetMaxSize(m_iSize+o->m_iSize);
//	for (z=0;z<o->m_iSize;z++)
//		m_pData[m_iSize+z] = o->m_pData[z];
	if (o->m_pData != NULL)
		memcpy(&m_pData[m_iSize],o->m_pData,o->m_iSize*sizeof(int));
	m_iSize += o->m_iSize;
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}	


void CxIntArray::SetAt(unsigned int pos, int f)
{
	BXIN;
//	unsigned int z;
#ifdef DEBUG_CINTARRAY
	bool s = false;
	mprintf("@ CxIntArray::Add(long)...");
#endif
	if (pos >= m_iMaxSize)
	{
#ifdef DEBUG_CINTARRAY
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
#ifdef DEBUG_CINTARRAY
	if (s)
		mprintf("@ done.\n");
	else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxIntArray::Add(int f)
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	bool s = false;
	mprintf("@ CxIntArray::Add(long)...");
#endif
//	mprintf("Add %d %X\n",f,this);
	if (m_iSize+1 > m_iMaxSize)
	{
#ifdef DEBUG_CINTARRAY
		s = true;
		mprintf("\n");
#endif
//		mprintf("Grow A.\n");
//		SetMaxSize(m_iMaxSize + m_iGrow);
		if (m_iMaxSize == 0)
			SetMaxSize(m_iMaxSize + m_iGrow);
				else SetMaxSize(m_iMaxSize*2);
//		mprintf("Grow B.\n");
	}
//	mprintf("B %X data=%X size=%d\n",this,m_pData,m_iSize);
	m_pData[m_iSize] = f;
//	mprintf("C %X\n",this);
	m_iSize++;
#ifdef DEBUG_CINTARRAY
	if (s)
		mprintf("@ done.\n");
			else mprintf("done.\n");
#endif
//	mprintf("Add Done %X.\n",this);
	BXOUT;
}

	
void CxIntArray::SetSize(unsigned int i)
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::SetSize(int): %d...",i);
#endif
	if (m_iSize == i)
		return;
	int *temp;

	try { temp = new int[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		if (m_iSize > i)
			memcpy(temp,m_pData,i*sizeof(int));
		else
			memcpy(temp,m_pData,m_iSize*sizeof(int));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxIntArray::SetMaxSize(unsigned int i)
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::SetMaxSize(int): %d...",i);
#endif
	int *temp;

	try { temp = new int[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,m_iSize*sizeof(int));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxIntArray::SetGrow(unsigned int i)
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::SetGrow(int): %d\n",i);
#endif
	m_iGrow = i;
	BXOUT;
}
		
	
void CxIntArray::RemoveAll()
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::RemoveAll():...");
#endif
	if (m_pData != NULL)
		delete[] m_pData;
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxIntArray::RemoveAll_KeepSize()
{
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::RemoveAll_KeepSize():...");
#endif
	m_iSize = 0;
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxIntArray::RemoveAt(unsigned int pos, unsigned int count)
{
	BXIN;
	int *temp;
		
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::RemoveAt(int, int): %d, %d...",pos,count);
#endif

	try { temp = new int[m_iSize-count]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize-count)*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(int));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(int));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxIntArray::RemoveAt_KeepSize(unsigned int pos, unsigned int count)
{
	BXIN;
	int *temp;
		
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::RemoveAt_KeepSize(int, int): %d, %d...",pos,count);
#endif

	try { temp = new int[m_iMaxSize]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)m_iMaxSize*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(int));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(int));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


int CxIntArray::Pop_KeepSize()
{
	if (m_iSize > 0)
	{
		m_iSize--;
		return m_pData[m_iSize];
	} else
		return 0;
}

	
void CxIntArray::InsertAt(int f, unsigned int pos)
{
	BXIN;
	int *temp;
		
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::InsertAt(long, int): %d...");
#endif

	try { temp = new int[m_iSize+1]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+1)*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(int));
		memcpy(&temp[pos+1],&m_pData[pos],(m_iSize-pos)*sizeof(int));
		delete[] m_pData;
	}
	temp[pos] = f;
	m_pData = temp;
	m_iSize++;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxIntArray::SetName(const char *name)
{
	(void)name;
#ifdef DEBUG_ARRAYS
	BXIN;
#ifdef DEBUG_CINTARRAY
	mprintf("@ CxIntArray::SetName(const char *): \"%s\"...",name);
#endif
	if (m_sName != NULL)
		delete[] m_sName;
	try { m_sName = new char[strlen(name)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,name);
#ifdef DEBUG_CINTARRAY
	mprintf("done.\n");
#endif
	BXOUT;
#endif
}
