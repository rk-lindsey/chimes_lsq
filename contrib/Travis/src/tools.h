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

#ifndef TOOLS_H
#define TOOLS_H

#include "config.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>
#include <wchar.h>
#include <ctype.h>
#include <setjmp.h>
//#include <assert.h>


//#include <sys/time.h>
//#include <termios.h>

class CxString;

/*************************************************************/

#define Pi (3.1415926535897932385)

#define MIN(a,b) (((a)<=(b)) ? (a) : (b))
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))

#define MIN3(a,b,c) MIN(MIN(a,b),c)
#define MAX3(a,b,c) MAX(MAX(a,b),c)

#define MIN4(a,b,c,d) MIN(MIN(a,b),MIN(c,d))
#define MAX4(a,b,c,d) MAX(MAX(a,b),MAX(c,d))

#define SQR(x)   (x)*(x)


#ifndef TARGET_LINUX
#define __PRETTY_FUNCTION__ ""
#endif


#ifdef TARGET_WINDOWS

#define isnan(X) _isnan(X)

#define GREY 8
#define BLUE 9
#define GREEN 10
#define CYAN 11
#define RED 12
#define PINK 13
#define YELLOW 14
#define WHITE 15

#define UML_AE ((unsigned char)142)
#define UML_ae ((unsigned char)132)
#define UML_OE ((unsigned char)153)
#define UML_oe ((unsigned char)148)
#define UML_UE ((unsigned char)154)
#define UML_ue ((unsigned char)129)
#define UML_ss ((unsigned char)225)

#else // Linux or generic

#define GREY 0
#define BLUE 4
#define GREEN 2
#define CYAN 6
#define RED 1
#define PINK 5
#define YELLOW 3
#define WHITE 7

#define UML_AE ((unsigned char)196)
#define UML_ae ((unsigned char)228)
#define UML_OE ((unsigned char)214)
#define UML_oe ((unsigned char)246)
#define UML_UE ((unsigned char)220)
#define UML_ue ((unsigned char)252)
#define UML_ss ((unsigned char)223)

#endif

#include <errno.h>
#include <sys/types.h>

class CxWordArray;
class CxIntArray;

#define max(a,b)  (((a) > (b)) ? (a) : (b))

#define SAVEPOS	mprintf(GREEN,"You may jump back to here by entering \"$\" <---\n\n"); if (setjmp(g_JumpBuf)!=0)	mprintf(GREEN,"\n<--- Going Back\n\n");


extern FILE *g_pLogFile;
extern jmp_buf g_JumpBuf;

char* fgets_bin(char *buf, int n, FILE *a);

void SavePosition();
void LoadPosition();

//void AskString(const char *s, char *buf, const char *def, ...);
//void AskString_ND(const char *s, char *buf, ...);

void AskString(const char *s, CxString *buf, const char *def, ...);
void AskString_ND(const char *s, CxString *buf, ...);

bool AskYesNo(const char *s, bool def, ...);
int AskUnsignedInteger(const char *s, int def, ...);
float AskFloat(const char *s, float def, ...);
bool AskYesNo_ND(const char *s, ...);
int AskUnsignedInteger_ND(const char *s, ...);
float AskFloat_ND(const char *s, ...);
int AskRangeInteger(const char *s, int mini, int maxi, int def, ...);
int AskRangeInteger_ND(const char *s, int mini, int maxi, ...);
float AskRangeFloat(const char *s, float mini, float maxi, float def, ...);
int AskInteger(const char *s, int def, ...);

int AskMolecule(const char *s);

void InitColor();
void TextColor(int fg);
void TextNormal();

int HalfBox();
float HalfBox_Exact();
int HalfBoxSq3();
bool FileExist(const char *s);
void FreeFileName(const char *s);
//void FreeFileName(char *pre, char *s, char *post);
void saxonize(char *buf);
void mprintf(const char *s, ...);
void mfprintf(FILE *a, const char *s, ...);
void mprintf(int color, const char *s, ...);
void inpprintf(const char *s, ...);
void eprintf(const char *s, ...);
//void myget(char *s);
void myget(CxString *st);
int mystricmp(const char *s1, const char *s2);
void RemoveDoubleBackslash(char *buf);
const char* TreeElement(const char *s);
bool ParseIntList(const char *s, CxIntArray *la);
bool ParseIntList(const char *s, CxWordArray *wa);
void SortSingleMolAtomTypes();
void xAddAtom(const char *s);

bool ContainsDigit(const char *s);
//void ReplaceDigits(char *s);
void ReplaceDigits(CxString *s);
const char* RemovePath(const char *s);

FILE *OpenFileWrite(const char *buf, bool text, CxString *out = NULL);
bool IsTTY(FILE *f);
double ZeroDivide(double a, double b);
bool isdigit(char c);
bool isnumeric(char c);
char *FormatBytes(double i);

float dec(float a, float dig);
void decomp(float &a, float &b, float &c);
float maxbound(float a, float r);
float minbound(float a, float r);
float majorticks(float lower, float upper);

int GetSignificantDigit(double d, int sig);
void CreateTicks(double mi, double ma, int &major, int &minor);

double MaxDiff_DoubleArray(double *p, int n);

//void ProtectCharacters(char *dest, const char *src, const char *rep, const char *prot);
void ProtectCharacters(CxString *dest, const char *src, const char *rep, const char *prot);

double NormalDistIntegral(double x);

double GammaLn(double xx);
double FactorialLn(int n);
double BinomialCoeff(int n, int k);

void HSV2RGB(float h, float s, float v, float &r, float &g, float &b);

inline int floori(int a, int b) {
	if (0 < (a ^ b)) {
		return a / b;
	} else {
		if (a % b != 0)
			return a / b - 1;
		else
			return a / b;
	}
}

#endif


