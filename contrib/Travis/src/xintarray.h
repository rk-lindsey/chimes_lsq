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

#ifndef XINTARRAY_H
#define XINTARRAY_H

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


class CxIntArray : public CxObject
{
public:
	CxIntArray();
	~CxIntArray();
	CxIntArray(const char *name);
	void SetName(const char *name);
	CxIntArray(CxIntArray &o);
	void CopyFrom(CxIntArray *o);
	void Add(int f);
	void Append(CxIntArray *o);
	void SetAt(unsigned int pos, int f);
	void SetSize(unsigned int i);
	void SetMaxSize(unsigned int i);
	void SetGrow(unsigned int i);
	void RemoveAll();
	void RemoveAll_KeepSize();
	void RemoveAt(unsigned int pos, unsigned int count);
	void RemoveAt_KeepSize(unsigned int pos, unsigned int count);
	void InsertAt(int f, unsigned int pos);
	int Pop_KeepSize();

	inline bool Contains(int i)
	{
		BXIN;
		int z;
		for (z=0;z<(int)m_iSize;z++)
			if (m_pData[z] == i)
				return true;
		return false;
		BXOUT;
	}

	inline int GetPosition(int i)
	{
		BXIN;
		int z;
		for (z=0;z<(int)m_iSize;z++)
			if (m_pData[z] == i)
				return z;
		return -1;
		BXOUT;
	}

	inline int &GetAt(unsigned int i)
	{
		BXIN;
#ifdef DEBUG_CINTARRAY
		mprintf("@ CxIntArray::GetAt(int): %d...",i);
#endif
#ifdef DEBUG_ARRAYS
		if (i >= m_iSize)
		{
			if (m_sName != NULL)
				eprintf("CxIntArray \"%s\" Boundary Error (%d/%d).\n",m_sName,i,m_iSize);
					else eprintf("CxIntArray Boundary Error (%d/%d).\n",i,m_iSize);
			abort();
		}
#endif
#ifdef DEBUG_CINTARRAY
		mprintf("done.\n");
#endif
		BXOUT;
		return m_pData[i];
	}

	inline int &operator [] (unsigned int i)
	{
#ifdef DEBUG_CINTARRAY
		mprintf("@ CxIntArray::operator [] (int): %d\n",i);
#endif
		return GetAt(i);
	}

	inline int GetSize()
	{
	#ifdef DEBUG_CINTARRAY
		mprintf("@ CxIntArray::GetSize(): %d\n",m_iSize);
	#endif
		return m_iSize;
	}	
	
private:	
	int *m_pData;
	unsigned int m_iSize;
	unsigned int m_iMaxSize;
	unsigned int m_iGrow;
#ifdef DEBUG_ARRAYS
	char *m_sName;
#endif
};

#endif
