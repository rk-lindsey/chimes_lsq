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


#ifndef XSTRING_H
#define XSTRING_H


#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


#define XVSPRINTF_LINUX(obj, format, pre)                   \
{                                                           \
	va_list XXXparams;                                      \
	int XXXi;                                               \
                                                            \
	va_start(XXXparams, pre);                               \
	XXXi = obj.Format_Internal(format,0,XXXparams);         \
	va_end(XXXparams);                                      \
                                                            \
	if (XXXi != 0)                                          \
	{                                                       \
		va_start(XXXparams, pre);                           \
		obj.Format_Internal(format,XXXi,XXXparams);         \
		va_end(XXXparams);                                  \
	}                                                       \
}


#define XVSPRINTF_WINDOWS(obj, format, pre)                 \
{                                                           \
	va_list XXXparams;                                      \
	int XXXi;                                               \
                                                            \
	XXXi = 0;                                               \
	do {                                                    \
		va_start(XXXparams, pre);                           \
		XXXi = obj.Format_Internal(format,XXXi,XXXparams);  \
		va_end(XXXparams);                                  \
	} while (XXXi != 0);                                    \
}



class CxString : public CxObject
{
public:
	CxString();
	~CxString();
	CxString(const char *s);
	CxString(const CxString &s);
	CxString(const CxString &s1, const CxString &s2);
	operator const char*() const;
	CxString & operator = (const CxString &s);
	CxString operator + (const CxString &s) const;
	bool operator == (const CxString &s) const;
	bool operator != (const CxString &s) const;
	bool operator > (const CxString &s) const;
	bool operator < (const CxString &s) const;
	bool operator >= (const CxString &s) const;
	bool operator <= (const CxString &s) const;
	void operator += (const CxString &s);
//	char& operator [] (int i);
	char& operator () (int i);
	int GetLength() const;
	int FindFirst(char c) const;
	int FindNext(int i, char c) const;
	int FindLast(char c) const;
	CxString Mid(int pos, int count) const;
	CxString Mid(int pos) const;
	CxString Left(int count) const;
	void Format(const char *s, ...);
	void sprintf(const char *s, ...);
//	void vsprintf(const char *s, va_list params);
	void strcat(const char *s);
	void strcpy(const char *s);
	void memcpy(const char *s, int size);
	void Dump();
	void SetBufSize(int i);
	const char* fgets(int len, FILE *a);
	void ToLowercase();
//	void ToUppercase();
	char* GetWritePointer();

	int Format_Internal(const char *s, int length, va_list params);
	
private:

	char *m_pData;
	int m_iBufLen;
};


CxString operator + (const char *s1, const CxString& s2);


#endif
