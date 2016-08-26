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

#include "xstring.h"
#include <string.h>
#include "maintools.h"


// Don't worry, this is only the first attempt. If this buffer is too small,
// larger buffers will be tried out. Strings of arbitrary length can be handled.

#define TEMP_BUF_SIZE (1024)

static char g_pTempBuf[TEMP_BUF_SIZE];


CxString operator + (const char *s1, const CxString& s2)
{
#ifdef DEBUG_CxString
	mprintf("@ operator + (const char *, const CxString &): \"%s\", \"%s\"...",s1,(const char*)s2);
#endif
	return CxString(s1,s2);
}


CxString::CxString()
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::CxString()\n");
#endif
//	m_pData = NULL;
//	m_iBufLen = 0;

	try { m_pData = new char[1]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)1*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_pData[0] = 0;
	m_iBufLen = 1;

	BXOUT;
}

	
CxString::~CxString()
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::~CxString(): \"%s\"...",m_pData);
#endif
	if (m_pData != NULL)
	{
		delete[] m_pData;
		m_pData = NULL;
	}
	m_iBufLen = 0;
#ifdef DEBUG_CxString
	mprintf("done.\n");
#endif
	BXOUT;
}

	
CxString::CxString(const char *s)
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::CxString(const char *): \"%s\"...",s);
#endif

	int i;

	if (s == NULL)
	{
		eprintf("CxString::CxString(const char *): NULL pointer passed.\n");
		abort();
	}

	i = strlen(s);

//	if (i != 0)
//	{
		try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		::strcpy(m_pData,s);

		m_iBufLen = i+1;
//	} else
//	{
//		m_pData = NULL;
//		m_iBufLen = 0;
//	}

#ifdef DEBUG_CxString
	mprintf("done.\n");
#endif
	BXOUT;
}

	
CxString::CxString(const CxString &s) : CxObject()
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::CxString(const CxString &): \"%s\"...",(const char*)s);
#endif

	int i;

	i = s.GetLength();

//	if (i != 0)
//	{
		try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		::strcpy(m_pData,s.m_pData);

		m_iBufLen = i+1;
//	} else
//	{
//		m_pData = NULL;
//		m_iBufLen = 0;
//	}

#ifdef DEBUG_CxString
	mprintf("done.\n");
#endif
	BXOUT;
}


CxString & CxString::operator = (const CxString &s)
{
	int i;

	i = s.GetLength();

	if (i+1 > m_iBufLen)
	{
		if (m_pData != NULL)
			delete[] m_pData;

		try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		::strcpy(m_pData,s.m_pData);

		m_iBufLen = i+1;
	} else
		::strcpy(m_pData,s.m_pData);

	return *this;
}

	
CxString::CxString(const CxString &s1, const CxString &s2)
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::CxString(const CxString &, const CxString &): \"%s\", \"%s\"...",(const char*)s1,(const char*)s2);
#endif

	int i;

	i = s1.GetLength()+s2.GetLength();

	if (i != 0)
	{
		if (m_iBufLen < i+1)
		{
			if (m_pData != NULL)
				delete[] m_pData;

			try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
			if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			m_iBufLen = i+1;

			::strcpy(m_pData,s1.m_pData);
		}

		::strcat(m_pData,s2.m_pData);
		
	} else
	{
//		m_iBufLen = 0;
//		m_pData = NULL;

		try { m_pData = new char[1]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)1*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_pData[0] = 0;
		m_iBufLen = 1;
	}

#ifdef DEBUG_CxString
	mprintf("done.\n");
#endif
	BXOUT;
}

	
CxString::operator const char*() const
{
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator const char* () const: \"%s\"\n",m_pData);
#endif
	return m_pData;
}

	
CxString CxString::operator + (const CxString &s) const
{
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator + (const CxString &) const: \"%s\", \"%s\"\n",m_pData,(const char*)s);
#endif
	return CxString(*this,s);
}

	
bool CxString::operator == (const CxString &s) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator == (const CxString &) const: \"%s\", \"%s\"\n",m_pData,(const char*)s);
#endif
	if (strcmp(m_pData,s.m_pData) == 0)
	{
		BXOUT;
		return true;
	} else
	{
		BXOUT;
		return false;
	}
}

	
bool CxString::operator != (const CxString &s) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator != (const CxString &) const: \"%s\", \"%s\"\n",m_pData,(const char*)s);
#endif
	if (strcmp(m_pData,s.m_pData) == 0)
	{
		BXOUT;
		return false;
	} else
	{
		BXOUT;
		return true;
	}
}

	
bool CxString::operator > (const CxString &s) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator > (const CxString &) const: \"%s\", \"%s\"\n",m_pData,(const char*)s);
#endif
	if (strcmp(m_pData,s.m_pData) > 0)
	{
		BXOUT;
		return true;
	} else
	{
		BXOUT;
		return false;
	}	
}

	
bool CxString::operator < (const CxString &s) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator < (const CxString &) const: \"%s\", \"%s\"\n",m_pData,(const char*)s);
#endif
	if (strcmp(m_pData,s.m_pData) < 0)
	{
		BXOUT;
		return true;
	} else
	{
		BXOUT;
		return false;
	}
}

	
bool CxString::operator >= (const CxString &s) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator >= (const CxString &) const: \"%s\", \"%s\"\n",m_pData,(const char*)s);
#endif
	if (strcmp(m_pData,s.m_pData) >= 0)
	{
		BXOUT;
		return true;
	} else
	{
		BXOUT;
		return false;
	}
}

	
bool CxString::operator <= (const CxString &s) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator <= (const CxString &) const: \"%s\", \"%s\"\n",m_pData,(const char*)s);
#endif
	if (strcmp(m_pData,s.m_pData) <= 0)
	{
		BXOUT;
		return true;
	} else
	{
		BXOUT;
		return false;
	}
}

	
void CxString::operator += (const CxString &s)
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator += (const CxString &): \"%s\", \"%s\"...",m_pData,(const char*)s);
#endif

	strcat((const char*)s);

#ifdef DEBUG_CxString
	mprintf("done.\n");
#endif
	BXOUT;
}
	

/*char& CxString::operator [] (int i)
{
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator [] (int i): \"%s\", %d\n",m_pData,i);
#endif
	if ((i < 0) || (i >= m_iBufLen))
	{
		eprintf("Error: CxString::operator []: Boundary Error (%d/%d).\n",i,m_iBufLen);
		abort();
	}
	return m_pData[i];
}*/


char& CxString::operator () (int i)
{
#ifdef DEBUG_CxString
	mprintf("@ CxString::operator () (int i): \"%s\", %d\n",m_pData,i);
#endif
	if ((i < 0) || (i >= m_iBufLen))
	{
		eprintf("Error: CxString::operator (): Boundary Error (%d/%d).\n",i,m_iBufLen);
		abort();
	}
	return m_pData[i];
}

  
int CxString::GetLength() const
{
#ifdef DEBUG_CxString
	mprintf("@ CxString::GetLength() const: \"%s\", %d\n",m_pData,strlen(m_pData));
#endif
	if (m_pData != NULL)
		return strlen(m_pData);
	else
		return 0;
}

	
int CxString::FindFirst(char c) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::FindFirst(char) const: \"%s\", '%c'\n",m_pData,c);
#endif
	char *p;
	p = strchr(m_pData,c);
	if (p == NULL)
	{
		BXOUT;
		return -1;
	}
	BXOUT;
	return p-m_pData;
}
	
	
int CxString::FindNext(int i, char c) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::FindNext(int, char) const: \"%s\", %d, '%c'\n",m_pData,i,c);
#endif
	char *p;
	p = strchr(&m_pData[i],c);
	if (p == NULL)
	{
		BXOUT;
		return -1;
	}
	BXOUT;
	return p-m_pData;
}	

	
int CxString::FindLast(char c) const
{
	BXIN;
#ifdef DEBUG_CxString
	mprintf("@ CxString::FindLast(char) const: \"%s\", '%c'\n",m_pData,c);
#endif
	char *p;
	p = strrchr(m_pData,c);
	if (p == NULL)
	{
		BXOUT;
		return -1;
	}
	BXOUT;
	return p-m_pData;
}

	
CxString CxString::Mid(int pos, int count) const
{
	BXIN;
	char *buf;
		
	try { buf = new char[count+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(count+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
#ifdef DEBUG_CxString
	mprintf("@ CxString::Mid(int, int) const: \"%s\", %d, %d...",m_pData,pos,count);
#endif
	::memcpy(buf,&m_pData[pos],count);
	buf[count] = 0;
#ifdef DEBUG_CxString
	mprintf("done: \"%s\"\n",buf);
#endif
	CxString s = CxString(buf);
	delete[] buf;
	BXOUT;
	return s;
}

	
CxString CxString::Mid(int pos) const
{
	BXIN;
	char *buf;
		
	try { buf = new char[strlen(&m_pData[pos])+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(strlen(&m_pData[pos])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
#ifdef DEBUG_CxString
	mprintf("@ CxString::Mid(int) const: \"%s\", %d...",m_pData,pos);
#endif
	::memcpy(buf,&m_pData[pos],strlen(&m_pData[pos])+1);
#ifdef DEBUG_CxString
	mprintf("done: \"%s\"\n",buf);
#endif
	CxString s = CxString(buf);
	delete[] buf;
	BXOUT;
	return s;
}

	
CxString CxString::Left(int count) const
{
	BXIN;
	char *buf;

	try { buf = new char[count+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(count+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
#ifdef DEBUG_CxString
	mprintf("@ CxString::Left(int) const: \"%s\", %d...",m_pData,count);
#endif
	::memcpy(buf,m_pData,count);
	buf[count] = 0;
#ifdef DEBUG_CxString
	mprintf("done: \"%s\"\n",buf);
#endif
	CxString s = CxString(buf);
	delete[] buf;
	BXOUT;
	return s;
}

	
void CxString::Dump()
{
	BXIN;
	mprintf(m_pData);
	BXOUT;
}


int CxString::Format_Internal(const char *s, int length, va_list params)
{
	BXIN;
	int i;
	
#ifdef DEBUG_CxString
	mprintf("@ CxString::Format_Internal(const char *s): \"%s\"...",s);
#endif

	if (s == NULL)
		return 0;

#ifdef TARGET_LINUX

	// First attempt: Try if string fits into TEMP_BUF_SIZE
	if (length == 0)
	{
		// Return value is the number of bytes (not including the terminal NULL character)
		// that should have been written to buffer if it was large enough
		i = vsnprintf(g_pTempBuf,TEMP_BUF_SIZE-1,s,params);

		// It fits
		if (i+1 <= TEMP_BUF_SIZE)
		{
			if (i+1 > m_iBufLen)
			{
				if (m_pData != NULL)
					delete[] m_pData;

				try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
				if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				::strcpy(m_pData,g_pTempBuf);

				m_iBufLen = i+1;
			} else
				::strcpy(m_pData,g_pTempBuf);

			return 0;
		} else // It does not fit - return required buffer size
			return i+1;

	} else // 2nd attempt with known buffer size
	{
		if (m_pData != NULL)
			delete[] m_pData;

		try { m_pData = new char[length+1]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)(length+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		vsnprintf(m_pData,length,s,params);

		m_iBufLen = length+1;

		return 0;
	}

#else
	
	// First attempt: Try if string fits into TEMP_BUF_SIZE
	if (length == 0)
	{
		// Accept two different implementations for maximum compatibility:
		// (1) If the buffer is too small, the return value is negative (absolute value without any meaning)
		// (2) Like above in the TARGET_LINUX case
		#ifdef TARGET_WINDOWS
			i = _vsnprintf(g_pTempBuf,TEMP_BUF_SIZE-1,s,params);
		#else
			i = vsnprintf(g_pTempBuf,TEMP_BUF_SIZE-1,s,params);
		#endif

		// It fits
		if ((i >= 0) && (i+1 <= TEMP_BUF_SIZE))
		{
			if (i+1 > m_iBufLen)
			{
				if (m_pData != NULL)
					delete[] m_pData;

				try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
				if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				::strcpy(m_pData,g_pTempBuf);

				m_iBufLen = i+1;
			} else
				::strcpy(m_pData,g_pTempBuf);

			return 0;
		} else // It does not fit - return larger guess for buffer size
			return 2*TEMP_BUF_SIZE;

	} else // n-th attempt
	{
		if (m_pData != NULL)
			delete[] m_pData;

		try { m_pData = new char[length]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)(length)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_iBufLen = length;

		#ifdef TARGET_WINDOWS
			i = _vsnprintf(m_pData,length-1,s,params);
		#else
			i = vsnprintf(m_pData,length-1,s,params);
		#endif

		// It fits
		if ((i >= 0) && (i+1 <= length))
			return 0;
		else // It does not fit - return larger guess for buffer size
			return 2*length;
	}

#endif

#ifdef DEBUG_CxString
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxString::Format(const char *s, ...)
{
#ifdef TARGET_LINUX
	XVSPRINTF_LINUX((*this),s,s);
#else
	XVSPRINTF_WINDOWS((*this),s,s);
#endif

/*	va_list params;
	int i;

#ifdef TARGET_LINUX

	va_start(params, s); 
	i = Format_Internal(s,0,params);
	va_end(params);

	if (i != 0)
	{
		va_start(params, s); 
		Format_Internal(s,i,params);
		va_end(params);
	}

#else

	i = 0;
	do {
		va_start(params, s); 
		printf("Trying buffer size %d...\n",i);
		i = Format_Internal(s,i,params);
		va_end(params);
	} while (i != 0);

#endif*/
}


void CxString::sprintf(const char *s, ...)
{
#ifdef TARGET_LINUX
	XVSPRINTF_LINUX((*this),s,s);
#else
	XVSPRINTF_WINDOWS((*this),s,s);
#endif

/*	va_list params;
	int i;

#ifdef TARGET_LINUX

	va_start(params, s); 
	i = Format_Internal(s,0,params);
	va_end(params);

	if (i != 0)
	{
		va_start(params, s); 
		Format_Internal(s,i,params);
		va_end(params);
	}

#else

	i = 0;
	do {
		va_start(params, s); 
		printf("Trying buffer size %d...\n",i);
		i = Format_Internal(s,i,params);
		va_end(params);
	} while (i != 0);

#endif*/
}


void CxString::strcat(const char *s)
{
	int i, j;
	char *tmp;

	if (s == NULL)
		return;

	j = strlen(s);

	if (j == 0)
		return;

	if (m_pData == NULL)
	{
		*this = CxString(s);
		return;
	}

	i = GetLength() + j;

	if (m_iBufLen >= i+1)
	{
		::strcat(m_pData,s);
	} else
	{
		try { tmp = new char[i+1]; } catch(...) { tmp = NULL; }
		if (tmp == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		::strcpy(tmp,m_pData);

		::strcat(tmp,s);

		m_iBufLen = i+1;

		delete[] m_pData;

		m_pData = tmp;
	}
}


void CxString::strcpy(const char *s)
{
	*this = CxString(s);
}


void CxString::SetBufSize(int i)
{
	char *tmp;

	if (m_iBufLen >= i)
		return;

	try { tmp = new char[i]; } catch(...) { tmp = NULL; }
	if (tmp == NULL) NewException((double)i*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (m_pData != NULL)
	{
		::memcpy(tmp,m_pData,m_iBufLen);
		delete[] m_pData;
	}

	m_iBufLen = i;

	m_pData = tmp;
}


void CxString::memcpy(const char *s, int size)
{
	if (size > m_iBufLen)
		SetBufSize(size);

	::memcpy(m_pData,s,size);
}


const char* CxString::fgets(int len, FILE *a)
{
	if (m_iBufLen <= len)
		SetBufSize(len+1);

	return ::fgets(m_pData,len,a);
}


void CxString::ToLowercase()
{
	if (m_pData != NULL)
		strtolower(m_pData);
}


/*
void CxString::ToUppercase()
{
	if (m_pData != NULL)
		strtoupper(m_pData);
}
*/


char* CxString::GetWritePointer()
{
	return m_pData;
}
