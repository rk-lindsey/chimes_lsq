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


#include "backtrace.h"
#include "travis.h"
#include "maintools.h"

#ifdef TARGET_LINUX
#include <execinfo.h>
#endif


CxObArray g_oaBackTrace;
char *g_sExeName;


void DumpBacktrace()
{
	int z;
	CxTracePoint *p;
	mprintf(YELLOW,">>> Backtrace\n");
	for (z=0;z<g_oaBackTrace.GetSize();z++)
	{
		p = (CxTracePoint*)g_oaBackTrace[z];
		if (p->m_sMessage[0] == 0)
			mprintf("    %d.) %s in %s after Line %d.\n",z+1,p->m_sFunction,p->m_sFile,p->m_iLine);
				else mprintf("    %d.) %s in %s after Line %d: \"%s\"\n",z+1,p->m_sFunction,p->m_sFile,p->m_iLine,p->m_sMessage);
	}
	mprintf(YELLOW,"<<< Backtrace\n");
}


#ifdef TARGET_LINUX


void DumpLinuxBacktrace()
{
	void *buffer[4096];
	char **strings, buf[256], *p;
	int n, z, i, j;
	bool showall;
	FILE *a;
	
//	mprintf(YELLOW,"*** Linux Backtrace ***\n");
	n = backtrace(buffer,4096);
	if (n <= 0)
	{
		eprintf("Could not retrieve stack trace.\n");
		return;
	}
//	mprintf("Stack contains %d Frames.\n",n);
	strings = backtrace_symbols(buffer,n);
	if (strings == NULL)
		eprintf("Could not retrieve debug symbol strings.\n");
/*	for (z=0;z<n;z++)
	{
		mprintf(GREEN," # %2d",n-z);
		mprintf(" [%016p]",buffer[z]);
		if (strings != NULL)
			mprintf(" \"%s\"",strings[z]);
		mprintf("\n");
	}*/
	mprintf(YELLOW,"Decoding line numbers... ");
	mprintf("(\"addr2line <addr> -f -i -C -e %s\")\n",g_sExeName);
	a = fopen(g_sExeName,"rb");
	if (a == NULL)
	{
		eprintf("Could not locate \"%s\" :-(\n",g_sExeName);
	} else
	{
		fclose(a);
		j = 0;
		for (z=0;z<n;z++)
			if ((buffer[z] < (void*)0xF00000) && (buffer[z] > (void*)0x400000))
				j++;
		if (j < 3)
			showall = true;
				else showall = false;
		for (z=0;z<n;z++)
		{
			if (showall || ((buffer[z] < (void*)0xF00000) && (buffer[z] > (void*)0x400000)))
			{
				mprintf(GREEN," # %2d ",n-z);
				sprintf(buf,"addr2line %p -f -i -C -e %s > tmp.txt",buffer[z],g_sExeName);
				system(buf);
				a = fopen("tmp.txt","rt");
				if (a == NULL)
				{
					eprintf("Error. Make sure you have the addr2line tool installed and writing permission to this directory.\n");
					continue;
				}
				mprintf("- [%010p] ",buffer[z]);
				i = 0;
				while (true)
				{
					fgets(buf,256,a);
					if (feof(a))
						break;
					if (strlen(buf) > 2)
					{
						if (i == 2)
						{
							i = 0;
							mprintf("\n                     ");
						}
						buf[strlen(buf)-1] = 0;
						p = strrchr(buf,'/');
						if (p == NULL)
							p = buf;
								else p++;
						mprintf("- %s ",p);
						i++;
					}
				}
				fclose(a);
				mprintf("\n");
			}
		}
		system("rm tmp.txt");
	}
}


#endif


void Crash()
{
	bool b;

	g_bCheckWrite = false;
	UninstallSignalHandler();
	b = g_bGlobalPsycho;
	g_bGlobalPsycho = false;
	mprintf(WHITE,"TRAVIS apparently has crashed :-( Sorry.\n");
	mprintf(WHITE,"Please send the Log File (travis.log) to the developer,\n");
	mprintf(WHITE,"making him able to analyze this error.\n");
	BTDUMP;
#ifdef TARGET_LINUX
	DumpLinuxBacktrace();
#endif
	if (g_bSMode)
	{
		mprintf("\n");
		g_bGlobalPsycho = b;
		PrintSMode();
		g_bGlobalPsycho = false;
		mprintf("\n");
	}
	eprintf("Delivering Control to Operating System.\n");
/*#ifdef TARGET_WINDOWS
	char buf[2];
	gets(buf);
#endif*/
	abort();
}


void CrashInt()
{
	bool b;

	UninstallSignalHandler();
	b = g_bGlobalPsycho;
	g_bGlobalPsycho = false;
	BTDUMP;
#ifdef TARGET_LINUX
	DumpLinuxBacktrace();
#endif
	if (g_bSMode)
	{
		mprintf("\n");
		g_bGlobalPsycho = b;
		PrintSMode();
		g_bGlobalPsycho = false;
		mprintf("\n");
	}
	eprintf("Delivering Control to Operating System.\n");
/*#ifdef TARGET_WINDOWS
	char buf[2];
	gets(buf);
#endif*/
	abort();
}


void SIGNAL_SEGV(int param)
{
	(void)param;  // Suppress "unused parameter" warning
	eprintf("\n\n*** SIGSEGV caught: Segmentation Fault ***\n");
	Crash();
}


void SIGNAL_FPE(int param)
{
	(void)param;  // Suppress "unused parameter" warning
	eprintf("\n\n*** SIGFPE caught: Floating Point Exception ***\n");
	Crash();
}


void SIGNAL_ILL(int param)
{
	(void)param;  // Suppress "unused parameter" warning
	eprintf("\n\n*** SIGILL caught: Illegal Instruction ***\n");
	Crash();
}


void SIGNAL_INT(int param)
{
	bool b;

	(void)param;  // Suppress "unused parameter" warning
	b = g_bGlobalPsycho;
	if (!g_bAbortAnalysis)
	{
		mprintf("\n\n*** SIGINT: Analysis softly interrupted by user ***\n");
		mprintf("    Press CTRL+C again to immediately stop execution.\n");
		g_bAbortAnalysis = true;
		signal(SIGINT,SIGNAL_INT);
		return;
	} else
	{
		g_bGlobalPsycho = false;
		eprintf("\n\n*** SIGINT caught: Hard interrupt by user ***\n");
		eprintf("Stopping Execution.\n");
		g_bGlobalPsycho = b;
		CrashInt();
		exit(0);
	}
}


void SIGNAL_ABRT(int param)
{
	(void)param;  // Suppress "unused parameter" warning
	eprintf("\n\n*** SIGABRT caught: Abnormal Program Termination ***\n");
	Crash();
}


void SIGNAL_TERM(int param)
{
	(void)param;  // Suppress "unused parameter" warning
	eprintf("\n\n*** SIGTERM caught: Terminating ***\n");
	Crash();
}


#ifdef TARGET_LINUX


void SIGNAL_HANGUP(int param)
{
	(void)param;  // Suppress "unused parameter" warning
	eprintf("\n\n*** SIGHUP caught: Trying to ignore hangup ^^ ***\n");
	signal(SIGHUP,SIGNAL_HANGUP);
}


#endif


void InstallSignalHandler()
{
	BTIN;
	signal(SIGSEGV,SIGNAL_SEGV);
	signal(SIGFPE,SIGNAL_FPE);
	signal(SIGILL,SIGNAL_ILL);
	signal(SIGINT,SIGNAL_INT);
	signal(SIGABRT,SIGNAL_ABRT);
	signal(SIGTERM,SIGNAL_TERM);
#ifdef TARGET_LINUX
	signal(SIGHUP,SIGNAL_HANGUP);
#endif
	BTOUT;
}


void UninstallSignalHandler()
{
	BTIN;
	signal(SIGSEGV,SIG_DFL);
	signal(SIGFPE,SIG_DFL);
	signal(SIGILL,SIG_DFL);
	signal(SIGINT,SIG_DFL);
	signal(SIGABRT,SIG_DFL);
	signal(SIGTERM,SIG_DFL);
#ifdef TARGET_LINUX
	signal(SIGHUP,SIG_DFL);
#endif
	BTOUT;
}


void NewException(double d, const char *filename, int line, const char *function)
{
	UninstallSignalHandler();
	g_bGlobalPsycho = false;
	eprintf("\n*** Out of Memory ***\n\n");
	eprintf("Caught an exception while allocating %s of memory.\n",FormatBytes(d));
	mprintf("The exception was raised in %s:%d",filename,line);
	if (function[0] != 0)
		mprintf(" - %s",function);
	mprintf(".\n");
	mprintf("It is likely that TRAVIS requested more RAM than available on this machine.\n");
	if (d >= 4000000000.0)
		mprintf("Apart from that, more than 4 GB cannot be allocated in one single call.\n");
	mprintf("Adjust your input parameters.\n\n");
	BTDUMP;
#ifdef TARGET_LINUX
	DumpLinuxBacktrace();
#endif
	eprintf("Delivering Control to Operating System.\n");
	abort();
}


void NewException(double d, const char *filename, int line, const char *function, const char *info)
{
	UninstallSignalHandler();
	g_bGlobalPsycho = false;
	eprintf("\n*** Out of Memory ***\n\n");
	eprintf("Caught an exception while allocating %s of memory.\n",FormatBytes(d));
	mprintf("The exception was raised in %s:%d",filename,line);
	if (function[0] != 0)
		mprintf(" - %s",function);
	if (info != NULL)
		mprintf(", concerning %s",info);
	mprintf(".\n");
	mprintf("It is likely that TRAVIS requested more RAM than available on this machine.\n");
	if (d >= 4000000000.0)
		mprintf("Apart from that, more than 4 GB cannot be allocated in one single call.\n");
	mprintf("Adjust your input parameters.\n\n");
	BTDUMP;
#ifdef TARGET_LINUX
	DumpLinuxBacktrace();
#endif
	eprintf("Delivering Control to Operating System.\n");
	abort();
}

