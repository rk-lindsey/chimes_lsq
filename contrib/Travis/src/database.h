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

#ifndef DATABASE_H
#define DATABASE_H

#include "tools.h"
#include "xobject.h"
#include "xstring.h"
#include "xobarray.h"
#include "backtrace.h"


class CDatabaseTable : public CxObject
{
	public:
	
		CDatabaseTable();
		~CDatabaseTable();
		void DumpTable();
		void DumpTreeTable(int depth, unsigned long bitmask, bool last);
		void DumpOutputFile(FILE *a, int depth);
		void AddInt(int x, int y, signed long val);
		void AddDouble(int x, int y, double val);
		void AddString(int x, int y, CxString val);
		void AddBool(int x, int y, bool val);
		
//private:
		unsigned char m_iWidth;
		unsigned short m_iHeight;
		CxObArray m_pValues;
};


class CDatabaseValue : public CxObject
{
public:

	CDatabaseValue();
	~CDatabaseValue();
	void DumpValue();
	void DumpType();
	void DumpOutputFile(FILE *a);
	
//private:
	CxString m_sName;
	void *m_pValue;
	unsigned char m_iType;  /* 0 - nichts, 1 - int, 2 - double, 3 - string, 4 - bool, 5 - array, 6 - table */
};


class CDatabaseNode : public CxObject
{
public:
	
	CDatabaseNode();
	~CDatabaseNode();
	CDatabaseNode(const char *s);
	CDatabaseNode* REC_FindNode(CxString s);
	CDatabaseNode* REC_FindAddNode(CxString s);
	CDatabaseValue* FindValue(CxString s);
	void REC_DumpTree(int depth, unsigned long bitmask, bool last);
	void REC_WriteOutputFile(FILE *a, int depth);
	
	CxString m_sName;
	
//private:
	CxObArray m_oaChildren;
	CxObArray m_oaValues;
	CDatabaseNode *m_pParent;
};


class CDatabase : public CxObject
{
public:
	CDatabase();
	~CDatabase();
	
	void AddNode(CxString s);
	void DumpTree();
	bool WriteOutputFile(CxString s);
	
	int GetElementType(CxString s);
	unsigned long GetInt(CxString s);
	double GetFloat(CxString s);
	const char* GetString(CxString s);
	bool GetBool(CxString s);
	
/*******************************************/
	
	void AddArray(CxString s);
	void AddArrayInt(CxString s, int i);
	void AddArrayFloat(CxString s, double d);
	void AddArrayString(CxString s, CxString st);
	void AddArrayBool(CxString s, bool b);

	int GetArraySize(CxString s);
	int GetArrayElementType(CxString s, int i);

	int GetArrayInt(CxString s, int i);
	double GetArrayFloat(CxString s, int i);
	const char* GetArrayString(CxString s, int i);
	bool GetArrayBool(CxString s, int i);

	void SetArrayInt(CxString s, int p, int i);
	void SetArrayFloat(CxString s, int p, double d);
	void SetArrayString(CxString s, int p, CxString st);
	void SetArrayBool(CxString s, int p, bool b);

/********* Noch zu implementieren *************/

	bool ExistNode(CxString s);
		
	int GetTableWidth(CxString s);
	int GetTableHeight(CxString s);
	int GetTableElementType(CxString s, int x, int y);
	
	unsigned long GetTableInt(CxString s, int x, int y);
	double GetTableFloat(CxString s, int x, int y);
	const char* GetTableString(CxString s, int x, int y);
	bool GetTableBool(CxString s, int x, int y);
	
	void SetTableInt(CxString s, int x, int y, int i);
	void SetTableFloat(CxString s, int x, int y, double d);
	void SetTableString(CxString s, int x, int y, CxString st);
	void SetTableBool(CxString s, int x, int y, bool b);
	
	void AddTable(CxString s);
	
/************************************************/
	
	void AddInt(CxString s, int i);
	void AddFloat(CxString s, double d);
	void AddString(CxString s, CxString s2);
	void AddBool(CxString s, bool val);

	void SetInt(CxString s, int i);
	void SetFloat(CxString s, double d);
	void SetString(CxString s, CxString s2);
	void SetBool(CxString s, bool val);

	void ParseInputFile(const char *s);
	
private:
	CDatabaseNode *m_pRoot;
};

#endif


