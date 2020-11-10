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

#include "xwordarray.h"


#ifndef DEBUG_ARRAYS
#define m_sName (NULL)
#endif


CxWordArray::CxWordArray()
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::CxWordArray()\n");
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


CxWordArray::CxWordArray(const char *name)
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::CxWordArray(const char *)\n");
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

	
CxWordArray::~CxWordArray()
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::~CxWordArray()\n");
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

	
CxWordArray::CxWordArray(CxWordArray &o) : CxObject()
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::CxWordArray(CxWordArray &)...");
#endif
	unsigned long z;
	m_iSize = o.m_iSize;
	m_iMaxSize = o.m_iMaxSize;
	m_iGrow = o.m_iGrow;

	try { m_pData = new unsigned short[m_iMaxSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(unsigned short),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	for (z=0;z<m_iSize;z++)
		m_pData[z] = o.m_pData[z];
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxWordArray::CopyFrom(CxWordArray *o)
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::CopyFrom(CxWordArray*)...");
#endif
	if (m_iMaxSize != o->m_iMaxSize)
	{
		m_iMaxSize = o->m_iMaxSize;
		if (m_pData != NULL)
			delete[] m_pData;

		try { m_pData = new unsigned short[m_iMaxSize]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(unsigned short),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	}
	m_iSize = o->m_iSize;
	m_iGrow = o->m_iGrow;
	if (o->m_pData != NULL)
		memcpy(m_pData,o->m_pData,m_iSize*sizeof(unsigned short));
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}	


void CxWordArray::Append(CxWordArray *o)
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::Append(CxWordArray*)...");
#endif
//	unsigned long z;
	if (m_iSize+o->m_iSize > m_iMaxSize)
		SetMaxSize(m_iSize+o->m_iSize);
//	for (z=0;z<o->m_iSize;z++)
//		m_pData[m_iSize+z] = o->m_pData[z];
	if (o->m_pData != NULL)
		memcpy(&m_pData[m_iSize],o->m_pData,o->m_iSize*sizeof(unsigned short));
	m_iSize += o->m_iSize;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}
	

void CxWordArray::SetAt(unsigned long pos, unsigned short f)
{
	BXIN;
//	unsigned long z;
#ifdef DEBUG_CWORDARRAY
	bool s = false;
	mprintf("@ CxWordArray::Add(unsigned short)...");
#endif
	if (pos >= m_iMaxSize)
	{
#ifdef DEBUG_CWORDARRAY
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
#ifdef DEBUG_CWORDARRAY
	if (s)
		mprintf("@ done.\n");
	else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxWordArray::Add(unsigned short f)
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	bool s = false;
	mprintf("@ CxWordArray::Add(unsigned short)...");
#endif
	if (m_iSize+1 > m_iMaxSize)
	{
#ifdef DEBUG_CWORDARRAY
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
#ifdef DEBUG_CWORDARRAY
	if (s)
		mprintf("@ done.\n");
			else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxWordArray::SetSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::SetSize(int): %d...",i);
#endif
	if (m_iSize == i)
		return;
	unsigned short *temp;

	try { temp = new unsigned short[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(unsigned short),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		if (m_iSize > i)
			memcpy(temp,m_pData,i*sizeof(unsigned short));
		else
			memcpy(temp,m_pData,m_iSize*sizeof(unsigned short));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxWordArray::GrowBy(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::GrowBy(unsigned long): %d...",i);
#endif
	if (i == 0)
		return;
	unsigned short *temp;

	try { temp = new unsigned short[m_iSize+i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+i)*sizeof(unsigned short),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,m_iSize*sizeof(unsigned short));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxWordArray::SetMaxSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::SetMaxSize(int): %d...",i);
#endif
	unsigned short *temp;

	try { temp = new unsigned short[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(unsigned short),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,m_iSize*sizeof(unsigned short));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxWordArray::SetGrow(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::SetGrow(int): %d\n",i);
#endif
	m_iGrow = i;
	BXOUT;
}
		
	
void CxWordArray::RemoveAll()
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::RemoveAll():...");
#endif
	if (m_pData != NULL)
		delete[] m_pData;
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxWordArray::RemoveAll_KeepSize()
{
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::RemoveAll_KeepSize():...");
#endif
	m_iSize = 0;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxWordArray::RemoveAt(unsigned long pos, unsigned long count)
{
	BXIN;
	unsigned short *temp;
		
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::RemoveAt(int, int): %d, %d...",pos,count);
#endif

	try { temp = new unsigned short[m_iSize-count]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize-count)*sizeof(unsigned short),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(unsigned short));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(unsigned short));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxWordArray::RemoveAt_KeepSize(unsigned long pos, unsigned long count)
{
	BXIN;
	unsigned short *temp;
		
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::RemoveAt_KeepSize(int, int): %d, %d...",pos,count);
#endif

	try { temp = new unsigned short[m_iMaxSize]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)m_iMaxSize*sizeof(unsigned short),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(unsigned short));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(unsigned short));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxWordArray::InsertAt(unsigned short f, unsigned long pos)
{
	BXIN;
	unsigned short *temp;
		
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::InsertAt(unsigned short, int): %d...");
#endif

	try { temp = new unsigned short[m_iSize+1]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+1)*sizeof(unsigned short),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(unsigned short));
		memcpy(&temp[pos+1],&m_pData[pos],(m_iSize-pos)*sizeof(unsigned short));
		delete[] m_pData;
	}
	temp[pos] = f;
	m_pData = temp;
	m_iSize++;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxWordArray::SetName(const char *name)
{
	(void)name;
#ifdef DEBUG_ARRAYS
	BXIN;
#ifdef DEBUG_CWORDARRAY
	mprintf("@ CxWordArray::SetName(const char *): \"%s\"...",name);
#endif
	if (m_sName != NULL)
		delete[] m_sName;
	try { m_sName = new char[strlen(name)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,name);
#ifdef DEBUG_CWORDARRAY
	mprintf("done.\n");
#endif
	BXOUT;
#endif
}


