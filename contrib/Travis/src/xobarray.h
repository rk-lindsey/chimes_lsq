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

#ifndef XOBARRAY_H
#define XOBARRAY_H

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


class CxObArray : public CxObject
{
public:
	
	CxObArray();
	CxObArray(CxObArray &o);
	~CxObArray();
	CxObArray(const char *name);
	void SetName(const char *name);
	void Add(CxObject *o);
	void SetAt(unsigned long pos, CxObject *o);


	inline CxObject* &GetAt(unsigned long i)
	{
#ifdef DEBUG_COBARRAY
		mprintf("@ CxObArray::GetAt(int): %d...",i);
#endif
#ifdef DEBUG_ARRAYS
		if (i >= m_iSize)
		{
			if (m_sName != NULL)
				eprintf("CxObArray \"%s\" Boundary Error (%d/%d).\n",m_sName,i,m_iSize);
					else eprintf("CxObArray Boundary Error (%d/%d).\n",i,m_iSize);
			abort();
		}
#endif
#ifdef DEBUG_COBARRAY
		mprintf("done.\n");
#endif
		return m_pData[i];
	}

	inline CxObject* &operator [] (unsigned long i)
	{
#ifdef DEBUG_COBARRAY
		mprintf("@ CxObArray::operator [] (int): %d\n",i);
#endif
		return GetAt(i);
	}

	inline int GetSize()
	{
	#ifdef DEBUG_COBARRAY
		mprintf("@ CxObArray::GetSize(): %d\n",m_iSize);
	#endif
		return m_iSize;
	}	

	void SetSize(unsigned long i);
	void SetSize_NoShrink(unsigned long i);
	void SetMaxSize(unsigned long i);
	void SetGrow(unsigned long i);
	void RemoveAll();
	void RemoveAll_KeepSize();
	void RemoveAt(unsigned long pos, unsigned long count);
	void RemoveAt_NoShrink(unsigned long pos, unsigned long count);
	void InsertAt(CxObject *o, unsigned long pos);
		
private:	
	CxObject **m_pData;
	unsigned long m_iSize;
	unsigned long m_iMaxSize;
	unsigned long m_iGrow;
#ifdef DEBUG_ARRAYS
	char *m_sName;
#endif
};

#endif
