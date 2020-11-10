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

#ifndef ASCIIART_H
#define ASCIIART_H

#include "xobject.h"
#include "xobarray.h"
#include "backtrace.h"


class CAsciiPiece : public CxObject
{
public:
	char GetAt(int x, int y);
	CAsciiPiece() { }
	~CAsciiPiece() { }
	unsigned short m_iWidth;
	unsigned short m_iHeight;
	char *m_pBuf;
};


class CAsciiArt : public CxObject
{
public:
	void NewLine();
	char GetChar();
	int GetWidth();
	int GetHeight();
	int m_iPosX;
	int m_iPosY;
	void SelectOne();
	CAsciiArt();
	~CAsciiArt();
	void AddPiece(const char *s);
	int m_iSelectedPiece;
	CxObArray m_oaAsciiPieces;
};

#endif
