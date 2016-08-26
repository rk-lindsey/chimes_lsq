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

#ifndef XLONGARRAY_H
#define XLONGARRAY_H

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


class CxLongArray : public CxObject
{
public:
	CxLongArray();
	~CxLongArray();
	CxLongArray(CxLongArray &o);
	CxLongArray(const char *name);
	void SetName(const char *name);
	void CopyFrom(CxLongArray *o);
	void Add(long f);
	void Append(CxLongArray *o);
	void SetAt(unsigned long pos, long f);
	void SetSize(unsigned long i);
	void SetMaxSize(unsigned long i);
	void SetGrow(unsigned long i);
	void RemoveAll();
	void RemoveAll_KeepSize();
	void RemoveAt(unsigned long pos, unsigned long count);
	void RemoveAt_KeepSize(unsigned long pos, unsigned long count);
	void InsertAt(long f, unsigned long pos);

	inline bool Contains(long i)
	{
		BXIN;
		int z;
		for (z=0;z<(int)m_iSize;z++)
			if (m_pData[z] == i)
				return true;
		return false;
		BXOUT;
	}

	inline int GetPosition(long i)
	{
		BXIN;
		int z;
		for (z=0;z<(int)m_iSize;z++)
			if (m_pData[z] == i)
				return z;
		return -1;
		BXOUT;
	}

	inline long &GetAt(unsigned long i)
	{
		BXIN;
#ifdef DEBUG_CLONGARRAY
		mprintf("@ CxLongArray::GetAt(int): %d...",i);
#endif
#ifdef DEBUG_ARRAYS
		if (i >= m_iSize)
		{
			if (m_sName != NULL)
				eprintf("CxLongArray \"%s\" Boundary Error (%d/%d).\n",m_sName,i,m_iSize);
					else eprintf("CxLongArray Boundary Error (%d/%d).\n",i,m_iSize);
			abort();
		}
#endif
#ifdef DEBUG_CLONGARRAY
		mprintf("done.\n");
#endif
		BXOUT;
		return m_pData[i];
	}

	inline long &operator [] (unsigned long i)
	{
#ifdef DEBUG_CLONGARRAY
		mprintf("@ CxLongArray::operator [] (int): %d\n",i);
#endif
		return GetAt(i);
	}

	inline int GetSize()
	{
	#ifdef DEBUG_CLONGARRAY
		mprintf("@ CxLongArray::GetSize(): %d\n",m_iSize);
	#endif
		return m_iSize;
	}	
	
private:	
	long *m_pData;
	unsigned long m_iSize;
	unsigned long m_iMaxSize;
	unsigned long m_iGrow;
#ifdef DEBUG_ARRAYS
	char *m_sName;
#endif
};

#endif
