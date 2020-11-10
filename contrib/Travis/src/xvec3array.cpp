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

#include "xvec3array.h"


#ifndef DEBUG_ARRAYS
#define m_sName (NULL)
#endif


CxVec3Array::CxVec3Array()
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::CxVec3Array()\n");
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


CxVec3Array::CxVec3Array(const char *name)
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::CxVec3Array(const char *)\n");
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

	
CxVec3Array::~CxVec3Array()
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::~CxVec3Array()\n");
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
	

CxVec3Array::CxVec3Array(CxVec3Array &o) : CxObject()
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::CxVec3Array(CxVec3Array &)...");
#endif
//	unsigned long z;
	m_iSize = o.m_iSize;
	m_iMaxSize = o.m_iMaxSize;
	m_iGrow = o.m_iGrow;
//	m_pData = new CxVector3*[m_iMaxSize];

	try { m_pData = new CxVector3[m_iMaxSize]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(CxVector3),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
/*	for (z=0;z<m_iSize;z++)
		m_pData[z] = new CxVector3(*o.m_pData[z]);
	for (z=m_iSize;z<(int)m_iMaxSize;z++)
		m_pData[z] = new CxVector3();*/
	if (o.m_pData != NULL)
		memcpy(m_pData,o.m_pData,sizeof(CxVector3)*m_iSize);
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxVec3Array::SetAt(unsigned long pos, CxVector3 o)
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	bool s = false;
	mprintf("@ CxVec3Array::Add(CxVector3)...");
#endif
	if (pos >= m_iMaxSize)
	{
#ifdef DEBUG_CVEC3ARRAY
		s = true;
		mprintf("\n");
#endif
		SetMaxSize(pos+1);
	}
//	*m_pData[pos] = *o;
	m_pData[pos] = o;
	if (pos >= m_iSize)
		m_iSize = pos+1;
#ifdef DEBUG_CVEC3ARRAY
	if (s)
		mprintf("@ done.\n");
	else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxVec3Array::Add(CxVector3 o)
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	bool s = false;
	mprintf("@ CxVec3Array::Add(CxVector3)...");
#endif
	if (m_iSize+1 > m_iMaxSize)
	{
#ifdef DEBUG_CVEC3ARRAY
		s = true;
		mprintf("\n");
#endif
//		SetMaxSize(m_iMaxSize + m_iGrow);
		if (m_iMaxSize == 0)
			SetMaxSize(m_iMaxSize + m_iGrow);
				else SetMaxSize(m_iMaxSize*2);
	}
//	*m_pData[m_iSize] = *o;
	m_pData[m_iSize] = o;
	m_iSize++;
#ifdef DEBUG_CVEC3ARRAY
	if (s)
		mprintf("@ done.\n");
			else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxVec3Array::SetSize(unsigned long i)
{
	BXIN;
//	int z;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::SetSize(int): %d...",i);
#endif
	if (m_iSize == i)
		return;
//	CxVector3 **temp;
	CxVector3 *temp;
	if (i == m_iSize)
	{
		BXOUT;
		return;
	}
//	temp = new CxVector3*[i];

	try { temp = new CxVector3[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(CxVector3),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		if (i > m_iSize)
		{
	/*		for (z=0;z<(int)m_iSize;z++)
			{
				temp[z] = new CxVector3(*m_pData[z]);
				delete m_pData[z];
			}
			for (z=m_iSize;z<(int)i;z++)
				temp[z] = new CxVector3();*/
			memcpy(temp,m_pData,sizeof(CxVector3)*m_iSize);
		} else
		{
	/*		for (z=0;z<(int)i;z++)
			{
				temp[z] = new CxVector3(*m_pData[z]);
				delete m_pData[z];
			}
			for (z=i;z<(int)m_iSize;z++)
				delete m_pData[z];*/
			memcpy(temp,m_pData,sizeof(CxVector3)*i);
		}
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxVec3Array::SetMaxSize(unsigned long i)
{
	BXIN;
//	int z;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::SetMaxSize(int): %d...",i);
#endif
//	CxVector3 **temp;
	CxVector3 *temp;
	if (i <= m_iSize)
	{
		eprintf("CxVec3Array::SetMaxSize(%d): Size %d > %d.\n",i,m_iSize,i);
		return;
	}
//	temp = new CxVector3*[i];

	try { temp = new CxVector3[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(CxVector3),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
/*	for (z=0;z<(int)m_iSize;z++)
	{
		temp[z] = new CxVector3(*m_pData[z]);
		delete m_pData[z];
	}
	for (z=m_iSize;z<(int)i;z++)
		temp[z] = new CxVector3();*/
	if (m_pData != NULL) {
		memcpy(temp,m_pData,sizeof(CxVector3)*m_iSize);
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxVec3Array::SetGrow(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::SetGrow(int): %d\n",i);
#endif
	m_iGrow = i;
	BXOUT;
}
		
	
void CxVec3Array::RemoveAll()
{
	BXIN;
//	int z;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::RemoveAll():...");
#endif
/*	for (z=0;z<(int)m_iMaxSize;z++)
		delete m_pData[z];*/
	if (m_pData != NULL)
		delete[] m_pData;
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxVec3Array::RemoveAll_KeepSize()
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::RemoveAll_KeepSize():...");
#endif
	m_iSize = 0;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxVec3Array::RemoveAt(unsigned long pos, unsigned long count)
{
	BXIN;
//	CxVector3 **temp;
	CxVector3 *temp;
//	int z;
		
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::RemoveAt(int, int): %d, %d...",pos,count);
#endif
//	temp = new CxVector3*[m_iSize-count];

	try { temp = new CxVector3[m_iSize-count]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize-count)*sizeof(CxVector3),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
/*	for (z=0;z<(int)pos;z++)
	{
		temp[z] = new CxVector3(*m_pData[z]);
		delete m_pData[z];
	}
	for (z=pos;z<(int)(m_iSize-count);z++)
	{
		temp[z] = new CxVector3(*m_pData[z+count]);
		delete m_pData[z+count];
	}*/
	if (m_pData != NULL) {
		memcpy(temp,m_pData,sizeof(CxVector3)*pos);
		memcpy(&temp[pos],&m_pData[pos+count],sizeof(CxVector3)*(m_iSize-count-pos));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxVec3Array::InsertAt(CxVector3 o, unsigned long pos)
{
	BXIN;
//	int z;
//	CxVector3 **temp;
	CxVector3 *temp;
		
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::InsertAt(CxVector3, int): %d...");
#endif
//	temp = new CxVector3*[m_iSize+1];

	try { temp = new CxVector3[m_iSize+1]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+1)*sizeof(CxVector3),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
/*	for (z=0;z<(int)pos;z++)
	{
		temp[z] = new CxVector3(*m_pData[z]);
		delete m_pData[z];
	}
	for (z=pos;z<(int)m_iSize;z++)
	{
		temp[z+1] = new CxVector3(*m_pData[z]);
		delete m_pData[z];
	}*/
	if (m_pData != NULL) {
		memcpy(temp,m_pData,sizeof(CxVector3)*pos);
		memcpy(&temp[pos+1],&m_pData[pos],sizeof(CxVector3)*(m_iSize-pos));
	//	*temp[pos] = *o;
		temp[pos] = o;
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize++;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxVec3Array::CopyFrom(CxVec3Array *a)
{
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::CopyFrom(CxVec3Array*)...");
#endif
	unsigned long z;
	if (m_iMaxSize != a->m_iMaxSize)
	{
		if (m_pData != NULL)
		{
/*			for (z=0;z<m_iSize;z++)
				delete m_pData[z];*/
			delete[] m_pData;
		}
		m_iMaxSize = a->m_iMaxSize;
//		m_pData = new CxVector3*[m_iMaxSize];

		try { m_pData = new CxVector3[m_iMaxSize]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(CxVector3),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
		
/*		for (z=0;z<m_iMaxSize;z++)
			m_pData[z] = new CxVector3();*/
	}
	m_iSize = a->m_iSize;
	m_iGrow = a->m_iGrow;
	for (z=0;z<m_iSize;z++)
		m_pData[z] = a->m_pData[z];
//		*m_pData[z] = *a->m_pData[z];
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxVec3Array::SetName(const char *name)
{
	(void)name;
#ifdef DEBUG_ARRAYS
	BXIN;
#ifdef DEBUG_CVEC3ARRAY
	mprintf("@ CxVec3Array::SetName(const char *): \"%s\"...",name);
#endif
	if (m_sName != NULL)
		delete[] m_sName;
	try { m_sName = new char[strlen(name)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,name);
#ifdef DEBUG_CVEC3ARRAY
	mprintf("done.\n");
#endif
	BXOUT;
#endif
}
