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


#include "xmemfile.h"
#include "tools.h"
#include "xstring.h"


#ifdef TARGET_LINUX
#include <unistd.h>
#endif


CxMemFile::CxMemFile()
{
	m_pBuffer = NULL;
	m_pPointer = NULL;
	m_iBufSize = 0;
}


CxMemFile::CxMemFile(CxMemFile *p)
{
	if (p->m_pBuffer != NULL)
	{
		m_pBuffer = new unsigned char[p->m_iBufSize];
		m_pPointer = p->m_pPointer;
		m_iBufSize = p->m_iBufSize;
		memcpy(m_pBuffer,p->m_pBuffer,m_iBufSize);
	} else
	{
		m_pBuffer = NULL;
		m_pPointer = NULL;
		m_iBufSize = 0;
	}
}


CxMemFile::~CxMemFile()
{
	ReleaseBuffer();
}


bool CxMemFile::SetSize(int i)
{
	unsigned char *p;

	p = new unsigned char[i];

	if (m_pBuffer != NULL)
	{
		if (i > m_iBufSize)
		{
			memcpy(p,m_pBuffer,m_iBufSize);
		} else
		{
			memcpy(p,m_pBuffer,i);
		}

		if (m_pPointer-m_pBuffer >= i)
			m_pPointer = &p[i-1];
		else
			m_pPointer = &p[m_pPointer - m_pBuffer];

		delete[] m_pBuffer;
	}

	m_iBufSize = i;
	m_pBuffer = p;

	return true;
}


void CxMemFile::ReleaseBuffer()
{
	if (m_pBuffer != NULL)
	{
		delete[] m_pBuffer;
		m_pBuffer = NULL;
	}
	m_pPointer = NULL;
	m_iBufSize = 0;
}


bool CxMemFile::Eof()
{
	if (m_pPointer >= &m_pBuffer[m_iBufSize-1])
		return true;
	else
		return false;
}


bool CxMemFile::Match(const char*s)
{
	SkipEmpty();
	if (memcmp(m_pPointer,s,strlen(s)) == 0)
	{
		m_pPointer += strlen(s);
		return true;
	} else
		return false;
}


bool CxMemFile::ReadFile(const char *s, bool text)
{
	FILE *a;
	int i, z;
	char buf[4096];

	if (text)
		a = fopen(s,"rt");
	else
		a = fopen(s,"rb");

	if (a == NULL)
	{
		mprintf("CMemFile::ReadFile(): Could not read \"%s\".\n",s);
		return false;
	}

	fseek(a,0,SEEK_END);
	i = ftell(a);
	fseek(a,0,SEEK_SET);


	mprintf("Reading %d bytes from %s ...\n",i,s);

	m_iBufSize = i+1;
	m_pBuffer = new unsigned char[i+1];
	m_pPointer = m_pBuffer;

	while (!feof(a))
	{
		z = fread(buf,1,4096,a);
		memcpy(m_pPointer,buf,z);
		m_pPointer += z;
		if (z < 4096)
			break;
	}

	m_pBuffer[m_iBufSize-1] = 0;

	m_pPointer = m_pBuffer;

//	mprintf("Done.\n");

	fclose(a);

	return true;
}


void CxMemFile::SkipEmpty()
{
	while ((m_pPointer < &m_pBuffer[m_iBufSize-1]) && ((*m_pPointer == ' ') || (*m_pPointer == '\r') || (*m_pPointer == '\n') || (*m_pPointer == '\t')))
		m_pPointer++;
}


void CxMemFile::ReverseSkipEmpty()
{
	while ((m_pPointer > m_pBuffer) && ((*m_pPointer == ' ') || (*m_pPointer == '\r') || (*m_pPointer == '\n') || (*m_pPointer == '\t')))
	{
//		mprintf("ReverseSkipEmpty: Skipping %d...\n",*((char*)m_pPointer));
		m_pPointer--;
	}
}


void CxMemFile::SkipEmpty(const char *sep)
{
	while ((strchr(sep,*m_pPointer) != 0) && (m_pPointer < &m_pBuffer[m_iBufSize-1]))
		m_pPointer++;
}


void CxMemFile::SkipWord()
{
	while ((m_pPointer < &m_pBuffer[m_iBufSize-1]) && (*m_pPointer != ' ') && (*m_pPointer != '\r') && (*m_pPointer != '\n') && (*m_pPointer != '\t'))
		m_pPointer++;
}


void CxMemFile::Seek(int pos)
{
	m_pPointer = &m_pBuffer[pos];
}


void CxMemFile::Create(int length)
{
	m_pBuffer = new unsigned char[length];
	m_iBufSize = length;
	Seek(0);
}


void CxMemFile::ReadWord(char *buf, int len)
{
	char *p;

	SkipEmpty();

	p = buf;
	while ((m_pPointer < &m_pBuffer[m_iBufSize-1]) && (*m_pPointer != ' ') && (*m_pPointer != '\r') && (*m_pPointer != '\n') && (*m_pPointer != '\t'))
	{
		*p = *m_pPointer;
		m_pPointer++;
		p++;
		if (p-buf+1 >= len)
		{
			mprintf("CxMemFile::ReadWord(): Warning: Buffer overflow prevented (%d).\n",len);
			*p = 0;
			return;
		}
	}
	*p = 0;
}


void CxMemFile::ReadWord(char *buf, int len, const char *sep)
{
	char *p;

	SkipEmpty(sep);

	p = buf;
	while ((m_pPointer < &m_pBuffer[m_iBufSize-1]) &&  (strchr(sep,*m_pPointer) == 0))
	{
		*p = *m_pPointer;
		m_pPointer++;
		p++;
		if (p-buf+1 >= len)
		{
			mprintf("CxMemFile::ReadWord(): Warning: Buffer overflow prevented (%d).\n",len);
			*p = 0;
			return;
		}
	}
	*p = 0;
}


int CxMemFile::scanf(const char *s, void* data)
{
	int i;
	static char buf[4096];
	unsigned char *p;

	while ((*m_pPointer == ' ') || (*m_pPointer == '\r') || (*m_pPointer == '\n'))
	{
		m_pPointer++;
		if (*m_pPointer == 0)
		{
			mprintf("CMemFile::scanf(): Unexpected end of file (1).\n");
			abort();
		}
	}
	p = m_pPointer;

	while ((*m_pPointer != ' ') && (*m_pPointer != '\r') && (*m_pPointer != '\n'))
	{
		m_pPointer++;
		if (*m_pPointer == 0)
		{
			mprintf("CMemFile::scanf(): Unexpected end of file (2).\n");
			abort();
		}
	}

	memcpy(buf,p,m_pPointer-p);
	buf[m_pPointer-p] = 0;

//	printf("Request: %s, Read \"%s\".\n",s,buf);

	i = sscanf(buf,s,data);

	return i;
}


int CxMemFile::fgets(char *buf, int len)
{
// 	unsigned char *p;
// 
// 	mprintf(GREEN, "%p\n", m_pPointer);
// 	p = m_pPointer;
// 
// 	while ((*m_pPointer != ' ') && (*m_pPointer != '\r') && (*m_pPointer != '\n'))
// 	{
// 		m_pPointer++;
// 		if (*m_pPointer == 0)
// 		{
// 			mprintf("CMemFile::scanf(): Unexpected end of file (2).\n");
// 			abort();
// 		}
// 	}
// 
// 	if (m_pPointer-p < len-1)
// 	{
// 		memcpy(buf,p,m_pPointer-p);
// 		buf[m_pPointer-p] = 0;
// 	} else
// 	{
// 		memcpy(buf,p,len-1);
// 		buf[len-1] = 0;
// 	}
// 
// 	return m_pPointer-p;
	
	int count = 0;
	unsigned char *start = m_pPointer;
	
	while ((*m_pPointer != '\n') && (*m_pPointer != 0) && (count < len - 1)) {
		m_pPointer++;
		count++;
	}
	if ((*m_pPointer == '\n') && (count < len - 1)) {
		m_pPointer++;
		count++;
	}
	memcpy(buf, start, count);
	buf[count] = 0;
	
	return count;
}


void CxMemFile::printf(const char *s, ...)
{
//	static char buf[4096];
	CxString buf;
	int i, j;
//	va_list args;
	unsigned char *p;

//	va_start(args,s);
//	vsprintf(buf,s,args);
//	buf.vsprintf(s,args);
//	va_end(args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(buf,s,s);
#else
	XVSPRINTF_WINDOWS(buf,s,s);
#endif


	i = strlen(buf);

	if (m_pPointer-m_pBuffer+i+1 > m_iBufSize)
	{
		p = new unsigned char[m_iBufSize+1024*64];
		if (m_iBufSize != 0)
		{
			j = m_pPointer-m_pBuffer;
			memcpy(p,m_pBuffer,m_pPointer-m_pBuffer+1);
			delete[] m_pBuffer;
			m_pBuffer = p;
			m_pPointer = m_pBuffer+j;
			m_iBufSize += 1024*64;
		} else
		{
			m_pBuffer = p;
			m_pPointer = m_pBuffer;
			m_iBufSize = 1024*64;
		}
	}

	memcpy(m_pPointer,(const char*)buf,i);
	m_pPointer += i;
	*m_pPointer = 0;
}


void CxMemFile::WriteFile(const char *s, bool text)
{
	FILE *a;
	unsigned char *p;
	int i, k;

	if (text)
		a = fopen(s,"wt");
	else
		a = fopen(s,"wb");

	p = m_pBuffer;

	i = m_iBufSize;

	while (i > 0)
	{
		if (i >= 4096)
		{
			k = fwrite(p,1,4096,a);
			p += k;
			i -= k;
		} else
		{
			k = fwrite(p,1,i,a);
			p += k;
			i -= k;
			break;
		}
	}

	mprintf("WriteFile(): %d Bytes written.\n",m_iBufSize);

	fclose(a);
}


bool CxMemFile::ReadFileSuccessive(const char *s, int lines, bool verbose)
{
	FILE *a;
	int i, k, l, z;
	static char buf[4096];
	unsigned char *tc;
	long t0;

	a = fopen(s,"rb");
	while (a == NULL) {
		if (verbose)
			mprintf("ReadFileSuccessive(): Waiting for file \"%s\"...\n", s);
#ifdef TARGET_LINUX
		usleep(1000000); // Wait 1 sec
#else
		eprintf("I have no usleep(), I am slightly in a hurry ;-)\n");
#endif
		a = fopen(s, "rb");
	}

	if (verbose)
		mprintf("ReadFileSuccessive(): Reading file \"%s\", expecting %d lines...\n",s,lines);

// 	if (a == NULL)
// 	{
// 		eprintf("CxMemFile::ReadFileSuccessive(): Could not open %s for reading.\n",s);
// 		return false;
// 	}

	t0 = time(NULL);

	if (m_pBuffer == NULL)
// 		Create(16384);
		Create(lines * 80);
	
	Seek(0);

	i = 0;
	l = 0; // Total line count read
	while (l < lines) // Until all expected lines have been read
	{
		k = fread(buf,1,4096,a);
		i += k;

		if (m_pPointer+k > m_pBuffer+m_iBufSize)
		{
			if (verbose)
				mprintf("  ReadFileSuccessive(): Increasing buffer size to %d bytes.\n",m_iBufSize*2);

			tc = new unsigned char[m_iBufSize*2];
			memcpy(tc,m_pBuffer,m_iBufSize);
			m_pPointer = tc + (m_pPointer - m_pBuffer);
			delete[] m_pBuffer;
			m_pBuffer = tc;
			m_iBufSize *= 2;
		}

		memcpy(m_pPointer,buf,k);
		m_pPointer += k;

		// Count line breaks in block just read
		for (z=0;z<k;z++)
			if (buf[z] == '\n')
				l++;
		
		if (l >= lines)
			break;

		if (feof(a))
		{
			if (time(NULL) > t0+180)
			{
				eprintf("CxMemFile::ReadFileSuccessive(): Timeout exceeded (180 sec) while waiting for more data in file.\n");
				fclose(a);
				return false;
			}

			if (verbose)
				mprintf("  ReadFileSuccessive(): EOF after reading %d bytes. %d/%d lines in total.\n",i,l,lines);
			
			i = 0;

#ifdef TARGET_LINUX
			usleep(1000000); // Wait 1 sec
#else
			eprintf("I have no usleep(), I am slightly in a hurry ;-)\n");
#endif

			clearerr(a); // Reset EOF flag
		}

	}
	fclose(a);

	m_pPointer++;
	*m_pPointer = 0;
	// This is probably smaller than the allocated block, but don't care - size of valid data matters
	m_iBufSize = m_pPointer - m_pBuffer + 1;

	if (verbose)
		mprintf("ReadFileSuccessive(): Reading complete.\n");

	Seek(0);

	return true;
}


