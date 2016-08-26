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


#include "tools.h"
#include "backtrace.h"
#include "xwordarray.h"
#include "travis.h"
#include "xintarray.h"


#ifdef TARGET_WINDOWS
#include <windows.h>
#include <wincon.h>
#include <io.h>
HANDLE g_hConsole;
#endif

#ifdef TARGET_LINUX
#include <unistd.h>
#endif


FILE *g_pLogFile;
jmp_buf g_JumpBuf;


/********* Beginn Quellcode *****************/

#ifdef TARGET_WINDOWS


void InitColor()
{
	BTIN;
#ifdef USE_COLOR
	if (g_bNoColor)
		return;
	g_hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
#endif
	BTOUT;
}


void TextColor(int c)
{
//	BTIN;
#ifdef USE_COLOR
	int flag;

	if (g_bNoColor)
		return;

	flag = 0;
	if ((c & 1) != 0)
		flag |= FOREGROUND_BLUE;
	if ((c & 2) != 0)
		flag |= FOREGROUND_GREEN;
	if ((c & 4) != 0)
		flag |= FOREGROUND_RED;
	if ((c & 8) != 0)
		flag |= FOREGROUND_INTENSITY;
	SetConsoleTextAttribute(g_hConsole, flag);
#endif
//	BTOUT;
}


void TextNormal()
{
//	BTIN;
#ifdef USE_COLOR
	if (g_bNoColor)
		return;
	SetConsoleTextAttribute(g_hConsole, FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);
#endif
//	BTOUT;
}


#elif defined(TARGET_LINUX)


void InitColor()
{
	BTIN;
#ifdef USE_COLOR
#endif
	BTOUT;
}


void TextColor(int c)
{
//	BTIN;
#ifdef USE_COLOR
	if (g_bNoColor)
		return;
	char command[13];
//	sprintf(command, "%c[%d;%d;40m", 0x1B, g_iColorIntensity, c + 30);
	sprintf(command, "%c[%d;%dm", 0x1B, g_iColorIntensity, c + 30);
	printf("%s", command);
#endif
//	BTOUT;
}


void TextNormal()
{
//	BTIN;
#ifdef USE_COLOR
	if (g_bNoColor)
		return;
	char command[13];
//	sprintf(command, "%c[0;37;40m", 0x1B);
	sprintf(command, "%c[0;0;0m", 0x1B);
	printf("%s", command);
#endif
//	BTOUT;
}


#else


void InitColor()
{
}


void TextColor(int c)
{
	(void)c;
}


void TextNormal()
{
}


#endif


char* fgets_bin(char *buf, int n, FILE *a)
{
	char *p;
	char *r;

	r = fgets(buf,n,a);
	p = &buf[strlen(buf)-1];

	while (((*p == '\r') || (*p == '\n')) && (p > buf))
		p--;

	p += 2;
	*p = 0;
	return r;
}


void LoadPosition()
{
	longjmp(g_JumpBuf,1);
}


bool AskYesNo(const char *s, bool def, ...)
{
	BTIN;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, def);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,def);
#else
	XVSPRINTF_WINDOWS(obuf,s,def);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		BTOUT;
		return def;
	}
	if ((mystricmp(buf,"y")==0) || (mystricmp(buf,"yes")==0))
	{
		BTOUT;
		return true;
	}
	if ((mystricmp(buf,"n")==0) || (mystricmp(buf,"no")==0))
	{
		BTOUT;
		return false;
	}
	eprintf("Wrong input. Enter \"y\" or \"n\".\n");
	inpprintf("! Wrong input. Enter \"y\" or \"n\".\n");
	goto _again;
	return false; // Never happens ^^
}


bool AskYesNo_ND(const char *s, ...)
{
	BTIN;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, s);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,s);
#else
	XVSPRINTF_WINDOWS(obuf,s,s);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		eprintf("There is no default value. Enter \"y\" or \"n\".\n");
		inpprintf("! There is no default value. Enter \"y\" or \"n\".\n");
		goto _again;
	}
	if ((mystricmp(buf,"y")==0) || (mystricmp(buf,"yes")==0))
	{
		BTOUT;
		return true;
	}
	if ((mystricmp(buf,"n")==0) || (mystricmp(buf,"no")==0))
	{
		BTOUT;
		return false;
	}
	eprintf("Wrong input. Enter \"y\" or \"n\".\n");
	inpprintf("! Wrong input. Enter \"y\" or \"n\".\n");
	goto _again;
	return false; // Never happens ^^
}


int AskUnsignedInteger(const char *s, int def, ...)
{
	BTIN;
	int i;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, def);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,def);
#else
	XVSPRINTF_WINDOWS(obuf,s,def);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		BTOUT;
		return def;
	}
	i = atoi(buf);
	if (i>0)
	{
		BTOUT;
		return i;
	}
	if (strcmp(buf,"0")==0)
	{
		BTOUT;
		return 0;
	}
	eprintf("Wrong input. Enter a positive integer number.\n");
	inpprintf("! Wrong input. Enter a positive integer number.\n");
	goto _again;
	return 0; // Never happens ^^
}


int AskInteger(const char *s, int def, ...)
{
	BTIN;
	int i;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, def);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,def);
#else
	XVSPRINTF_WINDOWS(obuf,s,def);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		BTOUT;
		return def;
	}
	i = atoi(buf);
	if (i!=0)
	{
		BTOUT;
		return i;
	}
	if (strcmp(buf,"0")==0)
	{
		BTOUT;
		return 0;
	}
	eprintf("Wrong input. Enter a integer number.\n");
	eprintf("! Wrong input. Enter a integer number.\n");
	goto _again;
	return 0; // Never happens ^^
}


int AskUnsignedInteger_ND(const char *s, ...)
{
	BTIN;
	int i;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, s);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,s);
#else
	XVSPRINTF_WINDOWS(obuf,s,s);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		eprintf("There is no default value. Enter a positive integer number.\n");
		inpprintf("! There is no default value. Enter a positive integer number.\n");
		goto _again;
	}
	i = atoi(buf);
	if (i>0)
	{
		BTOUT;
		return i;
	}
	if (strcmp(buf,"0")==0)
	{
		BTOUT;
		return 0;
	}
	eprintf("Wrong input. Enter a positive integer number.\n");
	inpprintf("! Wrong input. Enter a positive integer number.\n");
	goto _again;
	return 0; // Never happens ^^
}


int AskRangeInteger(const char *s, int mini, int maxi, int def, ...)
{
	BTIN;
	int i;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, def);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,def);
#else
	XVSPRINTF_WINDOWS(obuf,s,def);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		BTOUT;
		return def;
	}
	i = atoi(buf);
	if ((i > 0) && (i >= mini) && (i <= maxi))
	{
		BTOUT;
		return i;
	}
	if ((strcmp(buf,"0")==0) && (mini <= 0) && (maxi >= 0))
	{
		BTOUT;
		return 0;
	}
	eprintf("Wrong input. Enter an integer number between %d and %d.\n",mini,maxi);
	inpprintf("! Wrong input. Enter an integer number between %d and %d.\n",mini,maxi);
	goto _again;
	return 0; // Never happens ^^
}


float AskRangeFloat(const char *s, float mini, float maxi, float def, ...)
{
	BTIN;
	float f;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, def);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,def);
#else
	XVSPRINTF_WINDOWS(obuf,s,def);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		BTOUT;
		return def;
	}
	f = (float)atof(buf);
	if ((f > 0) && (f >= mini) && (f <= maxi))
	{
		BTOUT;
		return f;
	}
	if (((strcmp(buf,"0")==0) || (strcmp(buf,"0.0")==0)) && (mini <= 0.0) && (maxi >= 0.0))
	{
		BTOUT;
		return 0;
	}
	eprintf("Wrong input. Enter an floatin point number between %f and %f.\n",mini,maxi);
	inpprintf("! Wrong input. Enter an floatin point number between %f and %f.\n",mini,maxi);
	goto _again;
	return 0; // Never happens ^^
}


int AskRangeInteger_ND(const char *s, int mini, int maxi, ...)
{
	BTIN;
	int i;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, maxi);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,maxi);
#else
	XVSPRINTF_WINDOWS(obuf,s,maxi);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		eprintf("There is no default value. Enter a positive integer number.\n");
		inpprintf("! There is no default value. Enter a positive integer number.\n");
		goto _again;
	}
	i = atoi(buf);
	if ((i > 0) && (i >= mini) && (i <= maxi))
	{
		BTOUT;
		return i;
	}
	if ((strcmp(buf,"0")==0) && (mini <= 0) && (maxi >= 0))
	{
		BTOUT;
		return 0;
	}
	eprintf("Wrong input. Enter an integer number between %d and %d.\n",mini,maxi);
	inpprintf("! Wrong input. Enter an integer number between %d and %d.\n",mini,maxi);
	goto _again;
	return 0; // Never happens ^^
}


/*void AskString_ND(const char *s, char *buf, ...)
{
	BTIN;
	static char obuf[4096], *p;
	va_list args;
	va_start (args, buf);
	vsprintf (obuf,s, args);
	va_end (args);
_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(buf);
//	if (strcmp(buf,"$")==0)
//	{
//		BTOUT;
//		LoadPosition();
//	}
	if (strlen(buf)==0)
	{
		eprintf("There is no default value. Enter a character string.\n");
		inpprintf("! There is no default value. Enter a character string.\n");
		goto _again;
	}
	BTOUT;
	return;
}


void AskString(const char *s, char *buf, const char *def, ...)
{
	BTIN;
	static char obuf[4096], *p;
	va_list args;
	va_start (args, def);
	vsprintf (obuf,s, args);
	va_end (args);
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(buf);
//	if (strcmp(buf,"$")==0)
//	{
//		BTOUT;
//		LoadPosition();
//	}
	if (strlen(buf)==0)
	{
		strcpy(buf,def);
		BTOUT;
		return;
	}
	BTOUT;
	return;
}*/


void AskString_ND(const char *s, CxString *buf, ...)
{
	BTIN;
//	static char obuf[4096];
	CxString obuf;
	const char *p;
//	va_list args;

//	va_start (args, buf);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,buf);
#else
	XVSPRINTF_WINDOWS(obuf,s,buf);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(*buf) == 0)
	{
		eprintf("There is no default value. Enter a character string.\n");
		inpprintf("! There is no default value. Enter a character string.\n");
		goto _again;
	}
	BTOUT;
	return;
}


void AskString(const char *s, CxString *buf, const char *def, ...)
{
	BTIN;
//	static char obuf[4096];
	CxString obuf;
	const char *p;
//	va_list args;

//	va_start (args, def);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,def);
#else
	XVSPRINTF_WINDOWS(obuf,s,def);
#endif

	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(*buf) == 0)
	{
		buf->strcpy(def);
		BTOUT;
		return;
	}
	BTOUT;
	return;
}


float AskFloat_ND(const char *s, ...)
{
	BTIN;
	float d;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, s);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,s);
#else
	XVSPRINTF_WINDOWS(obuf,s,s);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		eprintf("There is no default value. Enter a floating point number.\n");
		inpprintf("! There is no default value. Enter a floating point number.\n");
		goto _again;
	}
	d = (float)atof(buf);
	if (d!=0)
	{
		BTOUT;
		return d;
	}
	if (strcmp(buf,"0")==0)
	{
		BTOUT;
		return 0;
	}
	eprintf("Wrong input. Enter a floating point number.\n");
	inpprintf("! Wrong input. Enter a floating point number.\n");
	goto _again;
	return 0; // Never happens ^^
}


float AskFloat(const char *s, float def, ...)
{
	BTIN;
	float d;
//	static char buf[256], obuf[4096];
	CxString buf, obuf;
	const char *p;
//	va_list args;

//	va_start (args, def);
//	vsprintf (obuf,s, args);
//	obuf.vsprintf(s,args);
//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(obuf,s,def);
#else
	XVSPRINTF_WINDOWS(obuf,s,def);
#endif

_again:
	mprintf(obuf);
	p = obuf;
	while ((*p == ' ') || (*p == '\n'))
		p++;
	inpprintf("! %s\n",p);
	myget(&buf);
/*	if (strcmp(buf,"$")==0)
	{
		BTOUT;
		LoadPosition();
	}*/
	if (strlen(buf)==0)
	{
		BTOUT;
		return def;
	}
	d = (float)atof(buf);
	if (d!=0)
	{
		BTOUT;
		return d;
	}
	if (strcmp(buf,"0")==0)
	{
		BTOUT;
		return 0;
	}
	eprintf("Wrong input. Enter a floating point number.\n");
	inpprintf("! Wrong input. Enter a floating point number.\n");
	goto _again;
	return 0; // Never happens ^^
}


int AskMolecule(const char *s)
{
	int z;
//	static char buf[4096], buf2[256];
	CxString buf, buf2;

//	sprintf(buf,"%s (",s);
	buf.sprintf("%s (",s);
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
//		sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//		strcat(buf,buf2);
//		if (z < g_oaMolecules.GetSize()-1)
//			strcat(buf,", ");
		buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
		buf.strcat(buf2);
		if (z < g_oaMolecules.GetSize()-1)
			buf.strcat(", ");
	}
//	strcat(buf,")? ");
	buf.strcat(")? ");
	return AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize()) - 1;
}


int HalfBox()
{
	BTIN;
	int i;

	i = 999999;
	if (g_bPeriodic)
	{
		if (g_bBoxNonOrtho)
		{
			i = ((int)(g_fBoxMinDiam/200))*100;
		} else
		{
			if (g_bPeriodicX && (g_fBoxX/2.0f < i))
				i = ((int)(g_fBoxX/200))*100;
			if (g_bPeriodicY && (g_fBoxY/2.0f < i))
				i = ((int)(g_fBoxY/200))*100;
			if (g_bPeriodicZ && (g_fBoxZ/2.0f < i))
				i = ((int)(g_fBoxZ/200))*100;
		}
	} else
	{
		if ((g_TimeStep.m_vMax[0]-g_TimeStep.m_vMin[0])/2.0f < i)
			i = ((int)((g_TimeStep.m_vMax[0]-g_TimeStep.m_vMin[0])/200))*100;
		if ((g_TimeStep.m_vMax[1]-g_TimeStep.m_vMin[1])/2.0f < i)
			i = ((int)((g_TimeStep.m_vMax[1]-g_TimeStep.m_vMin[1])/200))*100;
		if ((g_TimeStep.m_vMax[2]-g_TimeStep.m_vMin[2])/2.0f < i)
			i = ((int)((g_TimeStep.m_vMax[2]-g_TimeStep.m_vMin[2])/200))*100;
	}
	BTOUT; 
	return i;
}


float HalfBox_Exact()
{
	float f;

	f = 1.0e30f;
	if (g_bPeriodic)
	{
		if (g_bBoxNonOrtho)
		{
			f = g_fBoxMinDiam;
		} else
		{
			if (g_bPeriodicX && (g_fBoxX/2.0f < f))
				f = g_fBoxX;
			if (g_bPeriodicY && (g_fBoxY/2.0f < f))
				f = g_fBoxY;
			if (g_bPeriodicZ && (g_fBoxZ/2.0f < f))
				f = g_fBoxZ;
		}
	}
	return f;
}


int HalfBoxSq3()
{
	BTIN;
	int i;

	i = 999999;
	if (g_bPeriodic)
	{
		if (g_bPeriodicX && (g_fBoxX/2.0f*sqrt(3.0) < i))
			i = ((int)(g_fBoxX/200*sqrt(3.0)))*100;
		if (g_bPeriodicY && (g_fBoxY/2.0f*sqrt(3.0) < i))
			i = ((int)(g_fBoxY/200*sqrt(3.0)))*100;
		if (g_bPeriodicZ && (g_fBoxZ/2.0f*sqrt(3.0) < i))
			i = ((int)(g_fBoxZ/200*sqrt(3.0)))*100;
	} else
	{
		if ((g_TimeStep.m_vMax[0]-g_TimeStep.m_vMin[0])/2.0f*sqrt(3.0) < i)
			i = ((int)((g_TimeStep.m_vMax[0]-g_TimeStep.m_vMin[0])/200*sqrt(3.0)))*100;
		if ((g_TimeStep.m_vMax[1]-g_TimeStep.m_vMin[1])/2.0f*sqrt(3.0) < i)
			i = ((int)((g_TimeStep.m_vMax[1]-g_TimeStep.m_vMin[1])/200*sqrt(3.0)))*100;
		if ((g_TimeStep.m_vMax[2]-g_TimeStep.m_vMin[2])/2.0f*sqrt(3.0) < i)
			i = ((int)((g_TimeStep.m_vMax[2]-g_TimeStep.m_vMin[2])/200*sqrt(3.0)))*100;
	}
	BTOUT; 
	return i;
}


void mprintf(const char *s, ...)
{
//	static char buffer[32768];
	CxString buffer;
//	va_list args;
	int z, z2, zold;

//	va_start (args, s);

/*#ifdef TARGET_WINDOWS
	_vsnprintf(buffer,32767,s,args);
#elif defined(TARGET_LINUX)
	vsnprintf(buffer,32767,s,args);
#else
	vsprintf(buffer,s,args);
#endif*/

//	buffer.vsprintf(s,args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(buffer,s,s);
#else
	XVSPRINTF_WINDOWS(buffer,s,s);
#endif

//	va_end (args);

	if (g_bSaxonize)
	{
		buffer.SetBufSize(buffer.GetLength()*2);
		saxonize(buffer.GetWritePointer());
	}

	if (g_bGlobalPsycho)
	{
		zold = -1;
		for (z=0;z<(int)strlen(buffer);z++)
		{
			do {
				z2 = rand()%8+1;
			} while (zold == z2);
			zold = z2;
			switch(z2)
			{
				case 0: TextColor(GREY); break;
				case 1: TextColor(BLUE); break;
				case 2: TextColor(GREEN); break;
				case 3: TextColor(CYAN); break;
				case 4: TextColor(RED); break;
				case 5: TextColor(PINK); break;
				case 6: TextColor(YELLOW); break;
				case 7: TextColor(WHITE); break;
				case 8: TextNormal(); break;
			}
			printf("%c",buffer[z]);
		}
		TextNormal();
	} else
		printf("%s",(const char*)buffer);
	fflush(stdout);
	if (g_pLogFile != NULL)
	{
		mfprintf(g_pLogFile,"%s",(const char*)buffer);
		fflush(g_pLogFile);
	}
}


void mfprintf(FILE *a, const char *s, ...)
{
	va_list args;
	int i;
	
	va_start (args, s);
	i = vfprintf(a,s,args);
	if (g_bCheckWrite && ((i < 0) || (ferror(a) != 0)))
	{
#ifdef TARGET_LINUX
		char buf[512], buf2[512];
		snprintf(buf, sizeof(buf), "/proc/self/fd/%d", fileno(a));
		memset(buf2,0,512);
		if (readlink(buf, buf2, sizeof(buf2)) > 0)
		{
			eprintf("\nError writing data to file \"%s\".\n",buf2);
		} else 
		{
			eprintf("\nError writing data to file.\n");
		}
#else
		eprintf("\nError writing data to file.\n");
#endif
		eprintf("vfprintf returned %d. ferror returned %d. Errno=%d.\n",i,ferror(a),errno);
		eprintf("System error message: \"%s\".\n\n",strerror(errno));
		abort();
	}
	va_end (args);
//	fflush(a);
}


void inpprintf(const char *s, ...)
{
	va_list args;

	if (g_fInput != NULL)
	{
		va_start(args, s);
		vfprintf(g_fInput,s,args);
		va_end(args);
		fflush(g_fInput);
	}
}


void saxonize(char *buf)
{
	char tbuf[32768], *p, *q, *d;
	bool num, beg;
	int length;

	strcpy(tbuf,buf);

	p = tbuf;
	q = tbuf;
	d = buf;

	while (*q != 0)
	{
		while (*p == ' ')
		{
			p++;
			*d = ' ';
			d++;
		}
		q = p;
		num = false;
		while ((*q != ' ') && (*q != 0))
		{
			if (isdigit(*q) || (*q == '[') || (*q == '\\') || (*q == '~'))
				num = true;
			q++;
		}
		if (num)
		{
			while (p != q)
			{
				*d = *p;
				d++;
				p++;
			}
			continue;
		}
		length = q-p;
		beg = true;
		while (p != q)
		{
			switch(*p)
			{
				case 'A':
					if (beg && (tolower(*(p+1)) != 'r'))
						*d = UML_AE;
							 else *d = *p;
					break;

				case 'a':
					if (beg && (tolower(*(p+1)) != 'r'))
						*d = UML_ae;
							 else *d = *p;
					break;

				case 'c':
					if (beg)
					{
						if ((tolower(*(p+1)) == 'e') || (tolower(*(p+1)) == 'i'))
							*d = 's';
						else if (tolower(*(p+1)) == 'h') 
						{
							*(d++) = 's';
							*d = 'c';
						} else *d = 'g';
					} else if (p+1 == q)
						*d = 'g';
					else if ((tolower(*(p+1)) == 'a') || (tolower(*(p+1)) == 'o') || (tolower(*(p+1)) == 'u'))
						*d = 'g';
					else if ((tolower(*(p+1)) == 'e') || (tolower(*(p+1)) == 'i'))
						*d = 's';
					else if (tolower(*(p+1)) != 'h')
						*d = 'g';
					else *d = *p;
					break;

				case 'C':
					if (beg)
					{
						if ((tolower(*(p+1)) == 'e') || (tolower(*(p+1)) == 'i'))
							*d = 'S';
						else if (tolower(*(p+1)) == 'h') 
						{
							*(d++) = 'S';
							*d = 'c';
						} else *d = 'G';
					} else if (p+1 == q)
						*d = 'G';
					else if ((tolower(*(p+1)) == 'a') || (tolower(*(p+1)) == 'o') || (tolower(*(p+1)) == 'u'))
						*d = 'G';
					else if ((tolower(*(p+1)) == 'e') || (tolower(*(p+1)) == 'i'))
						*d = 'S';
					else if (tolower(*(p+1)) != 'h')
						*d = 'G';
					else *d = *p;
					break;

				case 'E':
					if (beg)
						*d = UML_AE;
							else *d = *p;
					break;

				case 'e':
					if (beg)
						*d = UML_ae;
							else *d = *p;
					break;

				case 'i':
					if (tolower(*(p+1)) == 'g')
					{
						*(d++) = 'i';
						*(d++) = 's';
						*(d++) = 'c';
						*d = 'h';
						p++;
					} else if ((tolower(*(p+1)) == 'c') && (tolower(*(p+2)) == 'h'))
					{
						*(d++) = 'i';
						*(d++) = 'd';
						*(d++) = 's';
						*(d++) = 'c';
						*d = 'h';
						p+=2;
					} else *d = *p;
					break;

				case 'I':
					if (tolower(*(p+1)) == 'g')
					{
						*(d++) = 'I';
						*(d++) = 'S';
						*(d++) = 'C';
						*d = 'H';
						p++;
					}  else if ((tolower(*(p+1)) == 'c') && (tolower(*(p+2)) == 'h'))
					{
						*(d++) = 'I';
						*(d++) = 'D';
						*(d++) = 'S';
						*(d++) = 'C';
						*d = 'H';
						p+=2;
					} else *d = *p;
					break;

				case 'k':
					*d = 'g';
					break;

				case 'K':
					*d = 'G';
					break;

				case 'p':
					*d = 'b';
					break;

				case 'P':
					*d = 'B';
					break;

				case 'q':
					if (tolower(*(p+1)) == 'u')
					{
						*(d++) = 'g';
						*d = 'w';
						p++;
					} else *d = *p;
					break;

				case 'Q':
					if (tolower(*(p+1)) == 'u')
					{
						*(d++) = 'G';
						*d = 'w';
						p++;
					} else *d = *p;
					break;

				case 's':
					if ((*(p+1) == 't') && beg)
					{
						*(d++) = 's';
						*(d++) = 'c';
						*(d++) = 'h';
						*d = 'd';
						p++;
					} else if ((*(p+1) == 'p') && beg)
					{
						*(d++) = 's';
						*(d++) = 'c';
						*(d++) = 'h';
						*d = 'b';
						p++;
					} else *d = *p;
					break;

				case 'S':
					if ((*(p+1) == 't') && beg)
					{
						*(d++) = 'S';
						*(d++) = 'c';
						*(d++) = 'h';
						*d = 'd';
						p++;
					} else if ((*(p+1) == 'p') && beg)
					{
						*(d++) = 'S';
						*(d++) = 'c';
						*(d++) = 'h';
						*d = 'b';
						p++;
					} else *d = *p;
					break;

				case 't':
					if ((tolower(*(p+1)) == 'i') && ((tolower(*(p+2)) == 'a') || (tolower(*(p+2)) == 'e') ||(tolower(*(p+2)) == 'i') || (tolower(*(p+2)) == 'o') || (tolower(*(p+2)) == 'u')))
					{
						*(d++) = 's';
						*(d++) = 'c';
						*(d++) = 'h';
						*d = 'i';
						p++;
					} else if (beg && (tolower(*(p+1)) == 'h'))
					{
						*d = 'd';
						p++;
					} else *d = 'd';
					break;

				case 'T':
					if ((tolower(*(p+1)) == 'i') && ((tolower(*(p+2)) == 'a') || (tolower(*(p+2)) == 'e') ||(tolower(*(p+2)) == 'i') || (tolower(*(p+2)) == 'o') || (tolower(*(p+2)) == 'u')))
					{
						*(d++) = 'S';
						*(d++) = 'C';
						*(d++) = 'H';
						*d = 'I';
						p++;
					} else if (beg && (tolower(*(p+1)) == 'h'))
					{
						*d = 'D';
						p++;
					} else *d = 'D';
					break;

				case 'u':
					if ((beg) && ((tolower(*(p+1)) == 's') || ((tolower(*(p+1)) == 'n') && (tolower(*(p+2)) == 'i'))))
					{
						*(d++) = 'j';
						*(d++) = 'u';
						*d = 'h';
					} else *d = *p;
					break;

				case 'U':
					if ((beg) && ((tolower(*(p+1)) == 's') || ((tolower(*(p+1)) == 'n') && (tolower(*(p+2)) == 'i'))))
					{
						*(d++) = 'J';
						*(d++) = 'u';
						*d = 'h';
					} else *d = *p;
					break;

				case 'x':
					if (length > 2)
					{
						*(d++) = 'g';
						*d = 's';
					} else *d = *p;
					break;

				case 'X':
					if (length > 2)
					{
						*(d++) = 'G';
						*d = 'S';
					} else *d = *p;
					break;

				case 'y':
					if ((p+1 == q) && (length > 2))
					{
						*(d++) = 'i';
						*(d++) = 'e';
						*d = 'h';
					} else *d = *p;
					break;

				case 'Y':
					if ((p+1 == q) && (length > 2))
					{
						*(d++) = 'I';
						*(d++) = 'E';
						*d = 'H';
					} else *d = *p;
					break;

				default:
					*d = *p;
			}
			d++;
			p++;
			beg = false;
		}
	}
	*d = 0;
}


void mprintf(int color, const char *s, ...)
{
//	char buffer[32768];
	CxString buffer;
//	va_list args;
	int z, z2, zold;
	
//	va_start (args, s);

/*#ifdef TARGET_WINDOWS
	_vsnprintf(buffer,32768,s,args);
#elif defined(TARGET_LINUX)
	vsnprintf(buffer,32768,s,args);
#else
	vsprintf(buffer,s,args);
#endif*/

//	buffer.vsprintf(s,args);

//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(buffer,s,s);
#else
	XVSPRINTF_WINDOWS(buffer,s,s);
#endif

	if (g_bSaxonize)
	{
		buffer.SetBufSize(buffer.GetLength()*2);
		saxonize(buffer.GetWritePointer());
	}

	if (g_bGlobalPsycho)
	{
		zold = -1;
		for (z=0;z<(int)strlen(buffer);z++)
		{
			do {
				z2 = rand()%8+1;
			} while (zold == z2);
			zold = z2;
			switch(z2)
			{
				case 0: TextColor(GREY); break;
				case 1: TextColor(BLUE); break;
				case 2: TextColor(GREEN); break;
				case 3: TextColor(CYAN); break;
				case 4: TextColor(RED); break;
				case 5: TextColor(PINK); break;
				case 6: TextColor(YELLOW); break;
				case 7: TextColor(WHITE); break;
				case 8: TextNormal(); break;
			}
			printf("%c",buffer[z]);
		}
	} else
	{
		TextColor(color);
		printf("%s",(const char*)buffer);
	}
	TextNormal();
	fflush(stdout);
	if (g_pLogFile != NULL)
	{
		mfprintf(g_pLogFile,"%s",(const char*)buffer);
		fflush(g_pLogFile);
	}
}


void eprintf(const char *s, ...)
{
//	char buffer[32768];
	CxString buffer;
//	va_list args;
//	int z, z2, zold;
	
//	va_start (args, s);

/*#ifdef TARGET_WINDOWS
	_vsnprintf(buffer,32768,s,args);
#elif defined(TARGET_LINUX)
	vsnprintf(buffer,32768,s,args);
#else
	vsprintf(buffer,s,args);
#endif*/

//	buffer.vsprintf(s,args);

//	va_end (args);

#ifdef TARGET_LINUX
	XVSPRINTF_LINUX(buffer,s,s);
#else
	XVSPRINTF_WINDOWS(buffer,s,s);
#endif

/*	if (g_bGlobalPsycho)
	{
		zold = -1;
		for (z=0;z<(int)strlen(buffer);z++)
		{
			do {
				z2 = rand()%8+1;
			} while (zold == z2);
			zold = z2;
			switch(z2)
			{
				case 0: TextColor(GREY); break;
				case 1: TextColor(BLUE); break;
				case 2: TextColor(GREEN); break;
				case 3: TextColor(CYAN); break;
				case 4: TextColor(RED); break;
				case 5: TextColor(PINK); break;
				case 6: TextColor(YELLOW); break;
				case 7: TextColor(WHITE); break;
				case 8: TextNormal(); break;
			}
			printf("%c",buffer[z]);
		}
	} else*/
	{
		TextColor(RED);
//		fprintf(stderr,buffer);
		printf("%s",(const char*)buffer);
	}
	TextNormal();
//	fflush(stderr);
	fflush(stdout);
	if (g_pLogFile != NULL)
	{
		fprintf(g_pLogFile,"%s",(const char*)buffer);
		fflush(g_pLogFile);
	}
}


/*
void myget(char *s)
{
	BTIN;
	FILE *t;
	char buf[256];
_beg:
	if (g_fInputFile != NULL)
	{
		if (fgets(s,200,g_fInputFile) == NULL)
		{
			if (feof(g_fInputFile))
			{
				eprintf("\nEnd of input file reached. Switching to keyboard input.\n\n");
				fclose(g_fInputFile);
				g_fInputFile = NULL;
				g_bInputRedirected = false;
				if (AskYesNo("    Do you want to write an input file for this session (y/n)? [yes] ",true))
				{
					mprintf("    Copying contents of old input file to input.txt ...\n");
					FreeFileName("intmp.txt");
					if (!FileExist(g_sInputFile))
					{
						eprintf("Error: Input file \"%s\" is no longer available.\n\n",g_sInputFile);
						eprintf("*** Please note that the input file must not be named \"intmp.txt\". ***\n\n");
						abort();
					}
					rename(g_sInputFile,"intmp.txt");
					t = fopen("intmp.txt","rt");
					if (t == NULL)
					{
						eprintf("Error: Could not rename input file \"%s\" to \"intmp.txt\".\n\n",g_sInputFile);
						abort();
					}
					g_fInput = OpenFileWrite("input.txt",true);
					while (!feof(t))
					{
						fgets(buf,256,t);
						if (feof(t))
							break;
						buf[strlen(buf)-1] = 0;
						fprintf(g_fInput,"%s\n",buf);
					}
					fclose(t);
				} 
				mprintf("Please continue with answering the last question above.\n\n");
				goto _beg;
			} else
			{
				eprintf("\nError while reading input file.\n\n");
				mprintf("\"End Of File\" has not been signaled. Some other weird error ^^\n\n");
				mprintf("TRAVIS has to stop here, sorry.\n\n");
				exit(0);
			}
		}
	} else
	{
		if (fgets(s,200,stdin) == NULL)
		{
			eprintf("\nError while reading input from stdin (probably end of input file reached).\n\n");
			mprintf("Your redirected input file seems to have not enough lines.\n");
			mprintf("The Operating System does not allow console applications with redirected input\n");
			mprintf("to access the keyboard, so you cannot enter the missing lines per hand.\n\n");
			mprintf("TRAVIS has to stop here. I'm sorry. Use the -i flag instead for reading input files.\n\n");
			exit(0);
		}
	}
	if (strlen(s) > 195)
		eprintf(" (Input line too long, truncated to 200 characters) ");
	if (strlen(s) != 0)
		s[strlen(s)-1] = 0;
	if (g_bInputRedirected)
	{
		if (strlen(s)==0)
		{
			mprintf(WHITE,"<RETURN>\n");
			return;
		}
		if (s[0] == '!')
			goto _beg; // Kommentar im Input File
		mprintf(WHITE,"%s\n",s);
	} else
	{
		if (g_pLogFile != NULL)
		{
			fprintf(g_pLogFile,"%s\n",s);
			fflush(g_pLogFile);
		}
		if (g_fInput != NULL)
		{
			fprintf(g_fInput,"%s\n",s);
			fflush(g_fInput);
		}
	}
	BTOUT; 
}
*/


void myget(CxString *st)
{
	BTIN;
	FILE *t;

	st->SetBufSize(4096);

_beg:
	if (g_fInputFile != NULL)
	{
//		if (fgets(st,4095,g_fInputFile) == NULL)
		if (st->fgets(4095,g_fInputFile) == NULL)
		{
			if (feof(g_fInputFile))
			{
				eprintf("\nEnd of input file reached. Switching to keyboard input.\n\n");
				fclose(g_fInputFile);
				g_fInputFile = NULL;
				g_bInputRedirected = false;
				if (AskYesNo("    Do you want to write an input file for this session (y/n)? [yes] ",true))
				{
					mprintf("    Copying contents of old input file to input.txt ...\n");
					FreeFileName("intmp.txt");
					if (!FileExist(g_sInputFile))
					{
						eprintf("Error: Input file \"%s\" is no longer available.\n\n",g_sInputFile);
						eprintf("*** Please note that the input file must not be named \"intmp.txt\". ***\n\n");
						abort();
					}
					rename(g_sInputFile,"intmp.txt");
					t = fopen("intmp.txt","rt");
					if (t == NULL)
					{
						eprintf("Error: Could not rename input file \"%s\" to \"intmp.txt\".\n\n",g_sInputFile);
						abort();
					}
					g_fInput = OpenFileWrite("input.txt",true);
					char buf[4096];
					while (!feof(t))
					{
						fgets(buf,4095,t);
						if (feof(t))
							break;
						buf[strlen(buf)-1] = 0;
						fprintf(g_fInput,"%s\n",buf);
					}
					fclose(t);
				} 
				mprintf("Please continue with answering the last question above.\n\n");
				goto _beg;
			} else
			{
				eprintf("\nError while reading input file.\n\n");
				mprintf("\"End Of File\" has not been signaled. Some other weird error ^^\n\n");
				mprintf("TRAVIS has to stop here, sorry.\n\n");
				exit(0);
			}
		}
	} else
	{
//		if (fgets(st,4095,stdin) == NULL)
		if (st->fgets(4095,stdin) == NULL)
		{
			eprintf("\nError while reading input from stdin (probably end of input file reached).\n\n");
			mprintf("Your redirected input file seems to have not enough lines.\n");
			mprintf("The Operating System does not allow console applications with redirected input\n");
			mprintf("to access the keyboard, so you cannot enter the missing lines per hand.\n\n");
			mprintf("TRAVIS has to stop here. I'm sorry. Use the -i flag instead for reading input files.\n\n");
			exit(0);
		}
	}
	if (st->GetLength() >= 4090)
		eprintf(" (Input line too long, truncated to 4090 characters) ");

	if (st->GetLength() != 0)
		(*st)(st->GetLength()-1) = 0;

	if (g_bInputRedirected)
	{
		if (st->GetLength() == 0)
		{
			mprintf(WHITE,"<RETURN>\n");
			return;
		}
		if ((*st)[0] == '!')
			goto _beg; // Kommentar im Input File
		mprintf(WHITE,"%s\n",(const char*)(*st));
	} else
	{
		if (g_pLogFile != NULL)
		{
			fprintf(g_pLogFile,"%s\n",(const char*)(*st));
			fflush(g_pLogFile);
		}
		if (g_fInput != NULL)
		{
			fprintf(g_fInput,"%s\n",(const char*)(*st));
			fflush(g_fInput);
		}
	}

	BTOUT; 
}


bool FileExist(const char *s)
{
	FILE *a;
	a = fopen(s,"rb");
	if (a == NULL)
		return false;
	fclose(a);
	return true;
}


void FreeFileName(const char *s)
{
	BTIN;
//	char tbuf[256], buf[256];
	CxString tbuf, buf;
	const char *p, *q;
	int z;

	if (FileExist(s))
	{
		p = strrchr(s,'\\');

		if (p == NULL)
			p = strrchr(s,'/');

		if (p != NULL)
		{
			p++;
//			memcpy(tbuf,s,p-s);
//			tbuf[p-s] = 0;
			tbuf.SetBufSize(p-s+1);
			tbuf.memcpy(s,p-s);
			tbuf(p-s) = 0;
		}

		q = strrchr(s,'.');
		z = 1;

		while (true)
		{
			if (p != NULL)
//				sprintf(buf,"%s#%d#%s",tbuf,z,p);
				buf.sprintf("%s#%d#%s",(const char*)tbuf,z,p);
			else
//				sprintf(buf,"#%d#%s",z,s);
				buf.sprintf("#%d#%s",z,s);

			if (!FileExist(buf))
			{
				if (strlen(buf) > 20)
				{
					if (q != NULL)
						mprintf(GREY,"[Renaming existing File %s to #%d#[...]%s ...",s,z,q);
					else
						mprintf(GREY,"[Renaming existing File %s to #%d#[...] ...",s,z);
				} else
					mprintf(GREY,"[Renaming existing File %s to %s ...",s,(const char*)buf);

				if (rename(s,buf)!=0)
				{
					eprintf("Error.]\n");
					eprintf("Could not rename file \"%s\".\n",s);
					eprintf("Make sure you have writing permission and the file is not opened by any external program.\n\n");
					eprintf("Aborting Program Execution.\n");
					BTOUT;
					abort();
				}

				if (FileExist(s))
				{
					eprintf("Error.]\n");
					eprintf("File \"%s\" still existed after attempt to rename.\n",s);
					eprintf("Make sure you have writing permission and the file is not opened by any external program.\n\n");
					eprintf("Aborting Program Execution.\n");
					BTOUT;
					abort();
				}

				mprintf(GREY,"OK.]\n");
				BTOUT;
				return;
			}
			z++;
			if (z > 10000)
			{
				eprintf("Error.]\n");
				eprintf("File %s has more than 10000 backups.\nHint: rm \\#*  ;-)\n\nAborting Program Execution.\n",s);
				BTOUT;
				abort();
			}
		}
	}
	BTOUT;
}


/*void FreeFileName(char *pre, char *s, char *post)
{
	BTIN;
	char buf[256];
	sprintf(buf,"%s%s%s",pre,s,post);
	FreeFileName(buf);
	BTOUT;
}*/


int mystricmp(const char *s1, const char *s2)
{
	BXIN;
	char *buf, *buf2;
	unsigned long z, i;

	try { buf = new char[strlen(s1)+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(strlen(s1)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { buf2 = new char[strlen(s2)+1]; } catch(...) { buf2 = NULL; }
	if (buf2 == NULL) NewException((double)(strlen(s2)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (z=0;z<strlen(s1);z++)
		buf[z] = toupper(s1[z]);
	buf[strlen(s1)]=0;
	for (z=0;z<strlen(s2);z++)
		buf2[z] = toupper(s2[z]);
	buf2[strlen(s2)]=0;
	i = strcmp(buf,buf2);
	delete[] buf;
	delete[] buf2;
	BXOUT;
	return i;
}


void RemoveDoubleBackslash(char *buf)
{
	BXIN;
	char *p, *q;
	
	p = buf;
	q = buf;

	while (*q != 0)
	{
		if ((*q == '\\') && (*(q+1) == '\\'))
			q++;
		*p = *q;
		p++;
		q++;
	}
	*p = 0;
	BXOUT;
}


const char* TreeElement(const char *s)
{
/*#ifdef TARGET_LINUX
	static char buf[16];
	if (strcmp(s,"|")==0)
	{
		sprintf(buf,"%c%c%c",0xe2,0x94,0x83);
		return buf;
	} else if (strcmp(s,"|---")==0)
	{
		sprintf(buf,"%c%c%c%c%c%c%c%c%c%c%c%c",0xe2,0x94,0xa0,0xe2,0x94,0x80,0xe2,0x94,0x80,0xe2,0x94,0x80);
		return buf;
	} else if (strcmp(s,"\\---")==0)
	{
		sprintf(buf,"%c%c%c%c%c%c%c%c%c%c%c%c",0xe2,0x94,0x96,0xe2,0x94,0x80,0xe2,0x94,0x80,0xe2,0x94,0x80);
		return buf;
	} else return "???";
#else*/
	if (strcmp(s,"|")==0)
		return "|";
	else if (strcmp(s,"|--")==0)
		return "|--";
	else if (strcmp(s,"\\--")==0)
		return "\\--";
	else return "???";
//#endif
}


bool ParseIntList(const char *s, CxIntArray *la)
{
	BTIN;
	const char *p, *q;
	char buf[32];
	int i, i2, z;
	bool m;

	i2 = -1;
	p = s;
	i = -1;
	m = false;
	while (*p != 0)
	{
		while ((*p == ' ') && (*p != 0))
			p++;
		q = p;
		if (isdigit(*q))
		{
			while (isdigit(*q))
				q++;
			if (q-p >= 31)
			{
				eprintf("ParseIntList: Entry too long (%d >= 31).\n",q-p);
				return false;
			}
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			if (m)
				i2 = atoi(buf);
			else
				i = atoi(buf);
		} else if (*q == '-')
		{
			if (i == -1)
			{
				eprintf("Error: \"-\" without preceding number.\n");
				BTOUT; 
				return false;
			}
			m = true;
			q++;
		} else if (*q == ',')
		{
			if (m)
			{
				if (i2 == -1)
					abort();
				for (z=i;z<=i2;z++)
					la->Add(z);
			} else
				la->Add(i);
			m = false;
			i = -1;
			i2 = -1;
			q++;
		}
//		_end:
		p = q;
	}
	
	if (m)
	{
		for (z=i;z<=i2;z++)
			la->Add(z);
	} else if (i != -1)
		la->Add(i);
	
	BTOUT; 
	return true;
}


bool ParseIntList(const char *s, CxWordArray *wa)
{
	BTIN;
	const char *p, *q;
	char buf[32];
	int i, i2, z;
	bool m;

	i2 = -1;
	p = s;
	i = -1;
	m = false;
	while (*p != 0)
	{
		while ((*p == ' ') && (*p != 0))
			p++;
		q = p;
		if (isdigit(*q))
		{
			while (isdigit(*q))
				q++;
			if (q-p >= 31)
			{
				eprintf("ParseIntList: Entry too long (%d >= 31).\n",q-p);
				return false;
			}
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			if (m)
				i2 = atoi(buf);
			else
				i = atoi(buf);
		} else if (*q == '-')
		{
			if (i == -1)
			{
				eprintf("Error: \"-\" without preceding number.\n");
				BTOUT; 
				return false;
			}
			m = true;
			q++;
		} else if (*q == ',')
		{
			if (m)
			{
				if (i2 == -1)
					abort();
				for (z=i;z<=i2;z++)
					wa->Add(z);
			} else wa->Add(i);
			m = false;
			i = -1;
			i2 = -1;
			q++;
		}
//		_end:
		p = q;
	}
	
	if (m)
	{
		for (z=i;z<=i2;z++)
			wa->Add(z);
	} else if (i != -1)
		wa->Add(i);
	
	BTOUT; 
	return true;
}


void SortSingleMolAtomTypes()
{
	BTIN;
	int z, z2, z3, z4, i, j;
	CSingleMolecule *s;

	for (z=0;z<g_oaSingleMolecules.GetSize();z++)
	{
//		mprintf("Molekuel %d\n",z+1);
		s = (CSingleMolecule*)g_oaSingleMolecules[z];
		for (z2=0;z2<s->m_baAtomIndex.GetSize()-1;z2++)
		{
//			printf("  Element %d\n",z2+1);
			i = -1;
			j = 999;
			for (z3=z2;z3<s->m_baAtomIndex.GetSize();z3++)
			{
				if (s->m_baAtomIndex[z3] < j)
				{
					j = s->m_baAtomIndex[z3];
					i = z3;
				}
			}
			if ((i != -1) && (i != z2))
			{
//				mprintf("  Vertausche Elemente %d <-> %d...\n",z2+1,i+1);
				j = s->m_baAtomIndex[i];
				s->m_baAtomIndex[i] = s->m_baAtomIndex[z2];
				s->m_baAtomIndex[z2] = j;
				for (z4=0;z4<s->m_oaMolAtoms.GetSize();z4++)
				{
					if (((CMolAtom*)s->m_oaMolAtoms[z4])->m_iType == z2)
						((CMolAtom*)s->m_oaMolAtoms[z4])->m_iType = i;
					else if (((CMolAtom*)s->m_oaMolAtoms[z4])->m_iType == i)
						((CMolAtom*)s->m_oaMolAtoms[z4])->m_iType = z2;
				}
/*				j = s->m_baAtomCount[i];
				s->m_baAtomCount[i] = s->m_baAtomCount[z2];
				s->m_baAtomCount[z2] = j;*/
//				mprintf("  Vertausche %d Offsets...\n",max(s->m_iAtomCount[i],s->m_iAtomCount[z2]));
/*				wa = (CxIntArray*)s->m_oaAtomOffset[i];
				s->m_oaAtomOffset[i] = (CxIntArray*)s->m_oaAtomOffset[z2];
				s->m_oaAtomOffset[z2] = wa;*/
/*				for (z3=0;z3<max(()->GetSize(),((CxIntArray*)s->m_oaAtomOffset[z2])->GetSize());z3++)
				{
//					mprintf("    %d: Vertausche %d <-> %d...",z3+1,s->m_iAtomOffset[z2][z3]+1,s->m_iAtomOffset[i][z3]+1);
					j = ((CxIntArray*)s->m_oaAtomOffset[z2])->GetAt(z3);
					((CxIntArray*)s->m_oaAtomOffset[z2])->SetAt(z3,((CxIntArray*)s->m_oaAtomOffset[i])->GetAt(z3));
					((CxIntArray*)s->m_oaAtomOffset[i])->SetAt(z3,j);
//					mprintf("Ok: %d <-> %d\n",s->m_iAtomOffset[z2][z3]+1,s->m_iAtomOffset[i][z3]+1);
				}*/
			}
		}
//		mprintf("Molekuel zu Ende.\n");
	}
	BTOUT;
}


FILE *OpenFileWrite(const char *buf, bool text, CxString *out)
{
	BTIN;
	FILE *a;
	char tmpbuf[256];

	if (strlen(buf) > 220)
	{
		memcpy(tmpbuf,buf,100);
		tmpbuf[100] = 0;
		strcat(tmpbuf,"...");
		memcpy(&tmpbuf[103],&buf[strlen(buf)-100],100);
		tmpbuf[203] = 0;
		mprintf("[");
		mprintf(RED,"Warning:");
		mprintf(" File name too long (%d), shortening to \"%s\"]\n",strlen(buf),tmpbuf);
	} else
		strcpy(tmpbuf,buf);

	FreeFileName(tmpbuf);
	if (text)
		a = fopen(tmpbuf,"wt");
			else a = fopen(tmpbuf,"wb");
	if (a == NULL)
	{
		eprintf("Could not open file \"%s\" for writing.\nBe sure to have writing permission for this directory.\n\nAborting execution.\n",tmpbuf);
		abort();
	}

	if (out != NULL)
		out->strcpy(tmpbuf);

	BTOUT;
	return a;
}


bool IsTTY(FILE *f)
{
	(void)f;
#ifdef TARGET_WINDOWS
	return (_isatty(_fileno(f))!=0);
#elif defined(TARGET_LINUX)
	return (isatty(fileno(f))!=0);
#else
	return true;
#endif
}


double ZeroDivide(double a, double b)
{
	if ((a != a) || (b != b))
		return 0;
	if (b == 0)
		return 0;
	return a / b;
}


bool isdigit(char c)
{
	if ((c >= '0') && (c <= '9'))
		return true;
			else return false;
}


bool isnumeric(char c)
{
	if (((c >= '0') && (c <= '9')) || (c == '.') || (c == '+') || (c == '-') || (c == 'e') || (c == 'E'))
		return true;
			else return false;
}


char *FormatBytes(double i)
{
	static char buf[64], buf2[64];

	if (i < 0)
	{
		buf2[0] = '-';
		buf2[1] = 0;
		i = fabs(i);
	} else buf2[0] = 0;
	if (i > 1024.0*1024.0*1024.0*2.0)
		sprintf(buf,"%.2f GB",i/1024.0/1024.0/1024.0);
	else if (i > 1024.0*1024.0*2.0)
		sprintf(buf,"%.2f MB",i/1024.0/1024.0);
	else if (i > 1024.0*2.0)
		sprintf(buf,"%.2f KB",i/1024.0);
	else sprintf(buf,"%.0f Bytes",i);
	strcat(buf2,buf);
	return buf2;
}


const char* RemovePath(const char *s)
{
	const char *p;

	p = strrchr(s,'/');
	if (p == NULL)
		p = strrchr(s,'\\');
	if (p == NULL)
		return s;
	return p+1;
}


float dec(float a, float dig)
{
	if (a == 0)
		return 1;
	return (float)pow(10.0,floor(log10(pow(10.0,dig)/a)));
}


void decomp(float &a, float &b, float &c)
{
	c = 1/dec(a,1);
	b = a/c;
}


float maxbound(float a, float r)
{
	return (float)ceil(a*r)/r;
}


float minbound(float a, float r)
{
	return (float)floor(a*r)/r;
}


float majorticks(float lower, float upper)
{
	float t1, t2, t3, r;
	t1 = (upper-lower)/6;
	decomp(t1,t2,t3);
	if (t2 <= 2)
		r = 2;
	else if (t2 <= 5)
		r = 5;
	else r = 10;
	return r * t3;
}


int GetSignificantDigit(double d, int sig)
{
	int t;

//	printf("GetDigit %G, %d\n",d,sig);
	if (d == 0)
		return 0;
	d = fabs(d);
//	printf("A: %G\n",d);
	d /= pow(10.0,floor(log10(d)));
//	printf("B: %G\n",d);
	d *= pow(10.0,sig);
//	printf("C: %G\n",d);
	t = (int)floor(d+0.01);
//	printf("D: %d\n",t);
	t = t % 10;
//	printf("E: %d\n",t);
	return t;
}

//FILE *g_fTicks = NULL;

void CreateTicks(double mi, double ma, int &major, int &minor)
{
	double t, tm, p, p2, px, l, t2, t3;
	int z, z2, z3, zm, z0, i, im, j;
	bool havezero;

	mprintf("       CreateTicks(): Optimizing axis ticks for range %g - %g...\n",mi,ma);
	/* First the major ticks */
	l = 1e30;
	i = -1;
	im = -1;
	for (z=4;z<9;z++)
	{
		for (zm=1;zm<6;zm++)
		{
			havezero = false;
	//		printf("*** Checke %d / %d = %d:\n",z,zm,z-2+(z-1)*zm);
			t = 1.0;
			for (z2=0;z2<z;z2++)
			{
				p = mi + (ma-mi)*z2/(z-1);
		//		printf(" *  %.10G  ",p);
				t2 = 0;
				for (z3=1;z3<6;z3++)
				{
					j = GetSignificantDigit(p,z3);
		//			printf("  <D%d=%d: ",z3+2,j);
					t3 = 0;
					switch(j)
					{
						case 0:
							break;
						case 5:
							if (GetSignificantDigit(p,z3+1) == 0)
							{
								if (z3 < 3)
									t3 = 0 * pow(10.0,z3-1);
										else t3 = 2 * pow(10.0,z3-1);
							} else t3 = 9 * pow(10.0,z3-1);
							break;
						case 2:
						case 4:
						case 6:
						case 8:
							if (GetSignificantDigit(p,z3+1) == 0)
								t3 = 5 * pow(10.0,z3-1);
									else t3 = 9 * pow(10.0,z3-1);
							break;
						default:
							t3 = 9 * pow(10.0,z3-1);
							break;
					}
		//			printf("%G> ",t3);
					t2 += t3;
				}
	//			printf("    Ergebnis %G\n",t2);
				t += t2;
				if (p == 0)
					havezero = true;
			}
			t /= z;

			tm = 1.0;
			for (z2=0;z2<z-1;z2++)
			{
				p = mi + (ma-mi)*z2/(z-1);
				p2 = mi + (ma-mi)*(z2+1)/(z-1);
				for (z0=0;z0<zm;z0++)
				{
					px = p + (p2-p)*(z0+1.0)/(zm+1.0);
		//			printf("      %.10G  ",px);
					t2 = 0;
					for (z3=1;z3<6;z3++)
					{
						j = GetSignificantDigit(px,z3);
	//					printf("  <D%d=%d: ",z3+2,j);
						t3 = 0;
						switch(j)
						{
							case 0:
								break;
							case 5:
								if (GetSignificantDigit(px,z3+1) == 0)
								{
									if (z3 < 3)
										t3 = 0 * pow(10.0,z3-1);
											else t3 = 2 * pow(10.0,z3-1);
								} else t3 = 9 * pow(10.0,z3-1);
								break;
							case 2:
							case 4:
							case 6:
							case 8:
								if (GetSignificantDigit(px,z3+1) == 0)
									t3 = 5 * pow(10.0,z3-1);
										else t3 = 9 * pow(10.0,z3-1);
								break;
							default:
								t3 = 9 * pow(10.0,z3-1);
								break;
						}
		//				printf("%G> ",t3);
						t2 += t3;
					}
	//				printf("    Ergebnis %G\n",t2);
					tm += t2;
				}
			}
			tm /= zm;

			t = 2*t + tm;

			switch(z)
			{
				case 4: t *= 3.0; break;
				case 5: t *= 2.0; break;
				case 6: t *= 1.0; break;
				case 7: t *= 2.0; break;
				case 8: t *= 6.0; break;
			}

			switch(z-2+(z-1)*zm)
			{
				case 1: t *= 2500.0; break;
				case 2: t *= 2500.0; break;
				case 3: t *= 2500.0; break;
				case 4: t *= 2500.0; break;
				case 5: t *= 750.0; break;
				case 6: t *= 750.0; break;
				case 7: t *= 250.0; break;
				case 8: t *= 250.0; break;
				case 9: t *= 100.0; break;
				case 10: t *= 100.0; break;
				case 11: t *= 10.0; break;
				case 12: t *= 10.0; break;
				case 13: t *= 10.0; break;
				case 14: t *= 8.0; break;
				case 15: t *= 5.0; break;
				case 16: t *= 1.0; break;
				case 17: t *= 1.0; break;
				case 18: t *= 1.0; break;
				case 19: t *= 1.0; break;
				case 20: t *= 1.0; break;
				case 21: t *= 1.0; break;
				case 22: t *= 1.0; break;
				case 23: t *= 2.0; break;
				case 24: t *= 2.0; break;
				case 25: t *= 5.0; break;
				case 26: t *= 5.0; break;
				case 27: t *= 5.0; break;
				case 28: t *= 10.0; break;
				case 29: t *= 10.0; break;
				case 30: t *= 10.0; break;
				case 31: t *= 10.0; break;
				case 32: t *= 50.0; break;
				case 33: t *= 50.0; break;
				case 34: t *= 50.0; break;
			}
			if (havezero)
				t -= 10000.0 / ((float)z-2+(z-1)*zm);
	//		printf("## Summe %G\n",t);
			if (t < l)
			{
				l = t;
				i = z;
				im = zm;
			}
		}
	}
//	printf("--> Gewonnen hat %d / %d = %d Ticks mit %G.\n",i,im,i-2+(i-1)*im,l);
	major = i;
	minor = im;

	mprintf("          Using %d major ticks: ",major);
	t = mi;
	for (z=0;z<major;z++)
	{
		mprintf("%g",t);
		if (z<major-1)
			mprintf(", ");
		t += (ma-mi)/(major-1.0);
	}
	mprintf(".\n");
	mprintf("          Using %d minor ticks per major tick: %g, ",minor,mi);
	t = mi;
	for (z=0;z<minor;z++)
	{
		t += (ma-mi)/(major-1.0)/(minor+1.0);
		mprintf("(%g), ",t);
	}
	mprintf("%g.\n",mi+(ma-mi)/(major-1.0));

/*	if (g_fTicks == NULL)
		g_fTicks = fopen("ticks.csv","wt");
	fprintf(g_fTicks,"%f;  %f;  %d;  %d;  %f;  %f",mi,ma,major,minor,l,mi);
	t = mi;
	for (z=0;z<minor;z++)
	{
		t += (ma-mi)/(major-1.0)/(minor+1.0);
		fprintf(g_fTicks,";  (%f)",t);
	}
	t = mi+(ma-mi)/(major-1.0);
	for (z=0;z<major-1;z++)
	{
		fprintf(g_fTicks,"; %f",t);
		t += (ma-mi)/(major-1.0);
	}
	fprintf(g_fTicks,"\n",t);
	fflush(g_fTicks);*/

/*	for (z=0;z<major;z++)
	{
		p = mi + (ma-mi)*z/(major-1);
		p2 = mi + (ma-mi)*(z+1)/(major-1);
		printf(" * %G\n",p);
		if (z+1 < major)
		{
			for (z2=0;z2<minor;z2++)
			{
				px = p + (p2-p)*(z2+1)/(minor+1);
				printf("   %G\n",px);
			}
		}
	}*/
}


double MaxDiff_DoubleArray(double *p, int n)
{
	double m, t;
	int z, z2;
	
	m = 0;
	for (z=0;z<n;z++)
	{
		for (z2=0;z2<n;z2++)
		{
			if (z == z2)
				continue;

			t = fabs(p[z]-p[z2]);
			if (t > m)
				m = t;
		}
	}
	return m;
}


/*void ProtectCharacters(char *dest, const char *src, const char *rep, const char *prot)
{
	const char *p, *r;
	char *q;

	p = src;
	q = dest;

	while (*p != 0)
	{
		if (strchr(rep,*p) != NULL)
		{
			r = prot;
			while (*r != 0)
			{
				*q = *r;
				q++;
				r++;
			}
		}
		*q = *p;
		p++;
		q++;
	}
	*q = 0;
}*/


void ProtectCharacters(CxString *dest, const char *src, const char *rep, const char *prot)
{
	const char *p, *r;
	char *q;
	int i, j;

	j = strlen(prot);
	i = strlen(src);

	p = src;
	while (*p != 0)
	{
		if (strchr(rep,*p) != NULL)
			i += j;
		p++;
	}

	dest->SetBufSize(i+1);

	p = src;
	q = dest->GetWritePointer();

	while (*p != 0)
	{
		if (strchr(rep,*p) != NULL)
		{
			r = prot;
			while (*r != 0)
			{
				*q = *r;
				q++;
				r++;
			}
		}
		*q = *p;
		p++;
		q++;
	}
	*q = 0;
}


// Cumulative density function of a standard normal random variable
// Taken from http://www.johndcook.com/cpp_phi.html
// Original source:  Handbook of Mathematical Functions by Abramowitz and Stegun, formula 7.1.26
double NormalDistIntegral(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}


/* Returns the value of the Gamma function at argument xx.
   Source: Numerical Recipes in C, 2nd Ed., Page 214 */
double GammaLn(double xx)
{
	double x, y, tmp, ser;
	static double cof[6] = { 76.18009172947146, -86.50532032941677,
	24.01409824083091, -1.231739572450155,
	0.1208650973866179e-2, -0.5395239384953e-5 };
	int j;

	y = xx;
	x = xx;
	tmp = x + 5.5;
	tmp -= (x+0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j=0;j<=5;j++)
		ser += cof[j] / ++y;

	return -tmp + log(2.5066282746310005 * ser / x);
}


/* Returns the natural logarithm of the factorial of n.
   Source: Numerical Recipes in C, 2nd Ed., Page 215 */
double FactorialLn(int n)
{
	if (n < 0)
		return -1.0;

	if (n <= 1)
		return 0.0;

	return GammaLn(n+1.0);
}


/* Returns the binomial coefficient "n over k".
   Source: Numerical Recipes in C, 2nd Ed., Page 215 */
double BinomialCoeff(int n, int k)
{
	return floor(0.5+exp(FactorialLn(n)-FactorialLn(k)-FactorialLn(n-k)));
}


void HSV2RGB(float h, float s, float v, float &r, float &g, float &b)
{
	int hi;
	float f, p, q, t;

	h *= 360.0f;

	hi = (int)floorf(h / 60.0f);
	f = (h / 60.0f) - hi;
	p = v * (1.0f - s);
	q = v * (1.0f - s * f);
	t = v * (1.0f - s * (1.0f - f) );

	switch(hi)
	{
		case 0:
		case 6:
			r = v;
			g = t;
			b = p;
			break;
		case 1:
			r = q;
			g = v;
			b = p;
			break;
		case 2:
			r = p;
			g = v;
			b = t;
			break;
		case 3:
			r = p;
			g = q;
			b = v;
			break;
		case 4:
			r = t;
			g = p;
			b = v;
			break;
		case 5:
			r = v;
			g = p;
			b = q;
			break;
	}

	return;
}


