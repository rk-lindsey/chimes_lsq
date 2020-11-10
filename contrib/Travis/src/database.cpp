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

#include "database.h"

int treecolor=RED, namecolor=GREEN, valuecolor=YELLOW, typecolor=BLUE, symbolcolor=PINK, nodecolor=CYAN;


CDatabaseTable::CDatabaseTable()
{
	BTIN;
	m_iWidth = 0;
	m_iHeight = 0;
	m_pValues.SetName("CDatabaseTable::m_pValues");
	BTOUT;
}


CDatabaseTable::~CDatabaseTable()
{
	BTIN;
	int x, y;
	for (y=0;y<m_iHeight;y++)
	{
		for (x=0;x<m_iHeight;x++)
			if (((CxObArray*)m_pValues.GetAt(y))->GetAt(x) != NULL)
				delete (CDatabaseValue*)((CxObArray*)m_pValues.GetAt(y))->GetAt(x);
		delete (CxObArray*)m_pValues.GetAt(y);
	}
	m_pValues.RemoveAll();
	BTOUT;
}


void CDatabaseTable::DumpTable()
{
	BTIN;
	int x, y;
	unsigned char *t, tw=0, ttw=0;
//	char buf[256];
	CxString buf;
	CDatabaseValue *v;
	
	try { t = new unsigned char[m_iWidth]; } catch(...) { t = NULL; }
	if (t == NULL) NewException((double)m_iWidth*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (x=0;x<m_iWidth;x++)
	{
		tw = 0;
		for (y=0;y<m_iHeight;y++)
		{
			v = (CDatabaseValue*)((CxObArray*)m_pValues.GetAt(y))->GetAt(x);
			if (v == NULL)
				continue;
			switch(v->m_iType)
			{
				case 1:
					ttw = ((int)log10((double)((signed long)v->m_pValue)))+1;
					if ((signed long)v->m_pValue < 0)
						ttw++;
					break;
					
				case 2:
					ttw = ((int)log10(*((double*)v->m_pValue)))+1;
					if (*((double*)v->m_pValue) < 0)
						ttw++;
					ttw+=5;
					break;
					
				case 3:
					ttw = strlen((char*)v->m_pValue)+2;
					break;
				
				case 4:
					ttw = 5;
					break;
			}
			if (ttw > tw)
				tw = ttw;
		}
		t[x] = tw;
	}
	
	for (y=0;y<m_iHeight;y++)
	{
		for (x=0;x<m_iWidth;x++)
		{
			v = (CDatabaseValue*)((CxObArray*)m_pValues.GetAt(y))->GetAt(x);
			if (v == NULL)
			{
				ttw = (t[x]-1)/2;
				tw = (t[x]-1)-ttw;
//				sprintf(buf,"  [%c%ds-%c%ds]  ",'%',tw,'%',ttw);
				buf.sprintf("  [%c%ds-%c%ds]  ",'%',tw,'%',ttw);
				mprintf(buf,"","");
				continue;
			}
			switch(v->m_iType)
			{
				case 1:
					tw = ((int)log10((double)((signed long)v->m_pValue)))+1;
					if ((signed long)v->m_pValue < 0)
						tw++;
					ttw = t[x]-tw;
//					sprintf(buf,"  [%c%ds%cd]  ",'%',ttw,'%');
					buf.sprintf("  [%c%ds%cd]  ",'%',ttw,'%');
					mprintf(buf,"",(signed long)v->m_pValue);
					break;
				case 2:
					tw = ((int)log10(*((double*)v->m_pValue)))+1;
					if (*((double*)v->m_pValue) < 0)
						tw++;
					tw+=5;
					ttw = t[x]-tw;
//					sprintf(buf,"  [%c%ds%c.4f]  ",'%',ttw,'%');
					buf.sprintf("  [%c%ds%c.4f]  ",'%',ttw,'%');
					mprintf(buf,"",*((double*)v->m_pValue));
					break;
				case 3:
					tw = strlen((char*)v->m_pValue)+2;
					ttw = t[x]-tw;
//					sprintf(buf,"  [\"%cs\"%c%ds]  ",'%','%',ttw);
					buf.sprintf("  [\"%cs\"%c%ds]  ",'%','%',ttw);
					mprintf(buf,(char*)v->m_pValue,"");
					break;
				case 4:
					ttw = (t[x]-5)/2;
					tw = (t[x]-5)-ttw;
					if (v->m_pValue == 0)
//						sprintf(buf,"  [%c%dsFALSE%c%ds]  ",'%',tw,'%',ttw);
						buf.sprintf("  [%c%dsFALSE%c%ds]  ",'%',tw,'%',ttw);
							else
//								sprintf(buf,"  [%c%dsTRUE %c%ds]  ",'%',tw,'%',ttw);
								buf.sprintf("  [%c%dsTRUE %c%ds]  ",'%',tw,'%',ttw);
					mprintf(buf,"","");
					break;
			}
		}
		mprintf("\n");
	}
	delete[] t;
	BTOUT;
}


void CDatabaseTable::DumpTreeTable(int depth, unsigned long bitmask, bool last)
{
	BTIN;
	int x, y, z2;
	unsigned char *t, tw=0, ttw=0;
//	char buf[256];
	CxString buf;
	CDatabaseValue *v;
	
	try { t = new unsigned char[m_iWidth]; } catch(...) { t = NULL; }
	if (t == NULL) NewException((double)m_iWidth*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (x=0;x<m_iWidth;x++)
	{
		tw = 0;
		for (y=0;y<m_iHeight;y++)
		{
			v = (CDatabaseValue*)((CxObArray*)m_pValues.GetAt(y))->GetAt(x);
			if (v == NULL)
				continue;
			switch(v->m_iType)
			{
				case 1:
					ttw = ((int)log10((double)((signed long)v->m_pValue)))+1;
					if ((signed long)v->m_pValue < 0)
						ttw++;
					break;
					
				case 2:
					ttw = ((int)log10(*((double*)v->m_pValue)))+1;
					if (*((double*)v->m_pValue) < 0)
						ttw++;
					ttw+=5;
					break;
					
				case 3:
					ttw = strlen((char*)v->m_pValue)+2;
					break;
				
				case 4:
					ttw = 5;
					break;
			}
			if (ttw > tw)
				tw = ttw;
		}
		t[x] = tw;
	}
	
	for (y=0;y<m_iHeight;y++)
	{
		TextColor(treecolor);
		mprintf("   ");
		for (z2=0;z2<depth;z2++)
			if ((bitmask & (unsigned long)pow(2.0,z2)) != 0)
				mprintf("%s   ",TreeElement("|"));
					else mprintf("    ");
		if (last)
			mprintf("      ");
				else mprintf("%s   ",TreeElement("|"));
		TextColor(valuecolor);
		
		for (x=0;x<m_iWidth;x++)
		{
			v = (CDatabaseValue*)((CxObArray*)m_pValues.GetAt(y))->GetAt(x);
			if (v == NULL)
			{
				ttw = (t[x]-1)/2;
				tw = (t[x]-1)-ttw;
//				sprintf(buf,"%c%ds-%c%ds",'%',tw,'%',ttw);
				buf.sprintf("%c%ds-%c%ds",'%',tw,'%',ttw);
				TextColor(symbolcolor);
				mprintf(" [");
				TextColor(valuecolor);
				mprintf(buf,"","");
				TextColor(symbolcolor);
				mprintf("] ");
				continue;
			}
			switch(v->m_iType)
			{
				case 1:
					tw = ((int)log10((double)((signed long)v->m_pValue)))+1;
					if ((signed long)v->m_pValue < 0)
						tw++;
					ttw = t[x]-tw;
//					sprintf(buf,"%c%ds%cd",'%',ttw,'%');
					buf.sprintf("%c%ds%cd",'%',ttw,'%');
					TextColor(symbolcolor);
					mprintf(" [");
					TextColor(valuecolor);
					mprintf(buf,"",(signed long)v->m_pValue);
					TextColor(symbolcolor);
					mprintf("] ");
					break;
				case 2:
					tw = ((int)log10(*((double*)v->m_pValue)))+1;
					if (*((double*)v->m_pValue) < 0)
						tw++;
					tw+=5;
					ttw = t[x]-tw;
//					sprintf(buf,"%c%ds%c.4f",'%',ttw,'%');
					buf.sprintf("%c%ds%c.4f",'%',ttw,'%');
					TextColor(symbolcolor);
					mprintf(" [");
					TextColor(valuecolor);
					mprintf(buf,"",*((double*)v->m_pValue));
					TextColor(symbolcolor);
					mprintf("] ");
					break;
				case 3:
					tw = strlen((char*)v->m_pValue)+2;
					ttw = t[x]-tw;
//					sprintf(buf,"\"%cs\"%c%ds",'%','%',ttw);
					buf.sprintf("\"%cs\"%c%ds",'%','%',ttw);
					TextColor(symbolcolor);
					mprintf(" [");
					TextColor(valuecolor);
					mprintf(buf,(char*)v->m_pValue,"");
					TextColor(symbolcolor);
					mprintf("] ");
					break;
				case 4:
					ttw = (t[x]-5)/2;
					tw = (t[x]-5)-ttw;
					if (v->m_pValue == 0)
//						sprintf(buf,"%c%dsFALSE%c%ds",'%',tw,'%',ttw);
						buf.sprintf("%c%dsFALSE%c%ds",'%',tw,'%',ttw);
					else
//						sprintf(buf,"%c%dsTRUE %c%ds",'%',tw,'%',ttw);
						buf.sprintf("%c%dsTRUE %c%ds",'%',tw,'%',ttw);
					TextColor(symbolcolor);
					mprintf(" [");
					TextColor(valuecolor);
					mprintf(buf,"","");
					TextColor(symbolcolor);
					mprintf("] ");
					break;
			}
		}
		mprintf("\n");
	}
	delete[] t;
	BTOUT;
}


void CDatabaseTable::DumpOutputFile(FILE *a, int depth)
{
	BTIN;
	int x, y;
	unsigned char *t, tw=0, ttw=0;
//	char buf[256];
	CxString buf;
	CDatabaseValue *v;
	
	try { t = new unsigned char[m_iWidth]; } catch(...) { t = NULL; }
	if (t == NULL) NewException((double)m_iWidth*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for (x=0;x<m_iWidth;x++)
	{
		tw = 0;
		for (y=0;y<m_iHeight;y++)
		{
			v = (CDatabaseValue*)((CxObArray*)m_pValues.GetAt(y))->GetAt(x);
			if (v == NULL)
				continue;
			switch(v->m_iType)
			{
				case 1:
					ttw = ((int)log10((double)((signed long)v->m_pValue)))+1;
					if ((signed long)v->m_pValue < 0)
						ttw++;
					break;
					
				case 2:
					ttw = ((int)log10(*((double*)v->m_pValue)))+1;
					if (*((double*)v->m_pValue) < 0)
						ttw++;
					ttw+=5;
					break;
					
				case 3:
					ttw = strlen((char*)v->m_pValue)+2;
					break;
				
				case 4:
					ttw = 5;
					break;
			}
			if (ttw > tw)
				tw = ttw;
		}
		t[x] = tw;
	}
	
	for (y=0;y<m_iHeight;y++)
	{
		for (x=0;x<depth;x++)
			mfprintf(a,"  ");
		for (x=0;x<m_iWidth;x++)
		{
			v = (CDatabaseValue*)((CxObArray*)m_pValues.GetAt(y))->GetAt(x);
			if (v == NULL)
			{
				ttw = (t[x]-1)/2;
				tw = (t[x]-1)-ttw;
//				sprintf(buf,"  %c%ds-%c%ds  ",'%',tw,'%',ttw);
				buf.sprintf("  %c%ds-%c%ds  ",'%',tw,'%',ttw);
				mfprintf(a,buf,"","");
				continue;
			}
			switch(v->m_iType)
			{
				case 1:
					tw = ((int)log10((double)((signed long)v->m_pValue)))+1;
					if ((signed long)v->m_pValue < 0)
						tw++;
					ttw = t[x]-tw;
//					sprintf(buf,"  %c%ds%cd  ",'%',ttw,'%');
					buf.sprintf("  %c%ds%cd  ",'%',ttw,'%');
					mfprintf(a,buf,"",(signed long)v->m_pValue);
					break;
				case 2:
					tw = ((int)log10(*((double*)v->m_pValue)))+1;
					if (*((double*)v->m_pValue) < 0)
						tw++;
					tw+=5;
					ttw = t[x]-tw;
//					sprintf(buf,"  %c%ds%c.4f  ",'%',ttw,'%');
					buf.sprintf("  %c%ds%c.4f  ",'%',ttw,'%');
					mfprintf(a,buf,"",*((double*)v->m_pValue));
					break;
				case 3:
					tw = strlen((char*)v->m_pValue)+2;
					ttw = t[x]-tw;
//					sprintf(buf,"  \"%cs\"%c%ds  ",'%','%',ttw);
					buf.sprintf("  \"%cs\"%c%ds  ",'%','%',ttw);
					mfprintf(a,buf,(char*)v->m_pValue,"");
					break;
				case 4:
					ttw = (t[x]-5)/2;
					tw = (t[x]-5)-ttw;
					if (v->m_pValue == 0)
//						sprintf(buf,"  %c%dsFALSE%c%ds  ",'%',tw,'%',ttw);
						buf.sprintf("  %c%dsFALSE%c%ds  ",'%',tw,'%',ttw);
					else
//						sprintf(buf,"  %c%dsTRUE %c%ds  ",'%',tw,'%',ttw);
						buf.sprintf("  %c%dsTRUE %c%ds  ",'%',tw,'%',ttw);
					mfprintf(a,buf,"","");
					break;
			}
		}
		if (y+1 < m_iHeight)
			mfprintf(a,"\n");
	}
	delete[] t;
	BTOUT;
}


void CDatabaseTable::AddInt(int x, int y, signed long val)
{
	BTIN;
	CDatabaseValue *v;
	CxObArray *o;
	int z, z2;
	
	if (y >= m_iHeight)
	{
		for (z=m_iHeight;z<=y;z++)
		{
			try { o = new CxObArray("CDatabaseTable::AddInt():o"); } catch(...) { o = NULL; }
			if (o == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			m_pValues.SetAt(z,(CxObject*)o);
			for (z2=0;z2<m_iWidth;z2++)
				o->SetAt(z2,NULL);
		}
		m_iHeight = y+1;
	}
	if (x >= m_iWidth)
	{
		for (z=m_iWidth;z<=x;z++)
			for (z2=0;z2<m_iHeight;z2++)
				((CxObArray*)m_pValues.GetAt(z2))->SetAt(z,(CxObject*)NULL);
		m_iWidth = x+1;
	}

	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	v->m_iType = 1;
	v->m_pValue = (void*)val;
	((CxObArray*)m_pValues.GetAt(y))->SetAt(x,(CxObject*)v);
	BTOUT;
}


void CDatabaseTable::AddDouble(int x, int y, double val)
{
	BTIN;
	CDatabaseValue *v;
	CxObArray *o;
	int z, z2;
	
	if (y >= m_iHeight)
	{
		for (z=m_iHeight;z<=y;z++)
		{
			try { o = new CxObArray("CDatabaseTable::AddDouble():o"); } catch(...) { o = NULL; }
			if (o == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			m_pValues.SetAt(z,(CxObject*)o);
			for (z2=0;z2<m_iWidth;z2++)
				o->SetAt(z2,NULL);
		}
		m_iHeight = y+1;
	}
	if (x >= m_iWidth)
	{
		for (z=m_iWidth;z<=x;z++)
			for (z2=0;z2<m_iHeight;z2++)
				((CxObArray*)m_pValues.GetAt(z2))->SetAt(z,(CxObject*)NULL);
		m_iWidth = x+1;
	}

	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v->m_iType = 2;

	double *t;
	try { t = new double; } catch(...) { t = NULL; }
	if (t == NULL) NewException((double)sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*t = val;
	v->m_pValue = (void*)t;
	((CxObArray*)m_pValues.GetAt(y))->SetAt(x,(CxObject*)v);
	BTOUT;
}


void CDatabaseTable::AddString(int x, int y, CxString val)
{
	BTIN;
	CDatabaseValue *v;
	CxObArray *o;
	int z, z2;
	
	if (y >= m_iHeight)
	{
		for (z=m_iHeight;z<=y;z++)
		{
			try { o = new CxObArray("CDatabaseTable::AddString():o"); } catch(...) { o = NULL; }
			if (o == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			m_pValues.SetAt(z,(CxObject*)o);
			for (z2=0;z2<m_iWidth;z2++)
				o->SetAt(z2,NULL);
		}
		m_iHeight = y+1;
	}
	if (x >= m_iWidth)
	{
		for (z=m_iWidth;z<=x;z++)
			for (z2=0;z2<m_iHeight;z2++)
				((CxObArray*)m_pValues.GetAt(z2))->SetAt(z,(CxObject*)NULL);
		m_iWidth = x+1;
	}

	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v->m_iType = 3;

	char *buf;
	try { buf = new char[val.GetLength()+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(val.GetLength()+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(buf,(const char*)val);
	v->m_pValue = (void*)buf;
	((CxObArray*)m_pValues.GetAt(y))->SetAt(x,(CxObject*)v);
	BTOUT;
}


void CDatabaseTable::AddBool(int x, int y, bool val)
{
	BTIN;
	CDatabaseValue *v;
	CxObArray *o;
	int z, z2;
	
	if (y >= m_iHeight)
	{
		for (z=m_iHeight;z<=y;z++)
		{
			try { o = new CxObArray("CDatabaseTable::AddBool():o"); } catch(...) { o = NULL; }
			if (o == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			m_pValues.SetAt(z,(CxObject*)o);
			for (z2=0;z2<m_iWidth;z2++)
				o->SetAt(z2,NULL);
		}
		m_iHeight = y+1;
	}
	if (x >= m_iWidth)
	{
		for (z=m_iWidth;z<=x;z++)
			for (z2=0;z2<m_iHeight;z2++)
				((CxObArray*)m_pValues.GetAt(z2))->SetAt(z,(CxObject*)NULL);
		m_iWidth = x+1;
	}

	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v->m_iType = 4;
	v->m_pValue = (void*)((unsigned long)(val?1:0));
	((CxObArray*)m_pValues.GetAt(y))->SetAt(x,(CxObject*)v);
	BTOUT;
}


CDatabaseValue::CDatabaseValue()
{
	BTIN;
	m_pValue = NULL;
	BTOUT;
}


CDatabaseValue::~CDatabaseValue()
{
}


void CDatabaseValue::DumpOutputFile(FILE *a)
{
	BTIN;
	int z;
	switch(m_iType)
	{
		case 0:
			mfprintf(a,"<NO VALUE>");
			break;
		case 1:
			mfprintf(a,"%ld",(long)m_pValue);
			break;
		case 2:
			mfprintf(a,"%g",*(double*)m_pValue);
			break;
		case 3:
			mfprintf(a,"\"%s\"",(char*)m_pValue);
			break;
		case 4:
			if (m_pValue == 0)
				mfprintf(a,"FALSE");
					else mfprintf(a,"TRUE");
			break;
		case 5:
			mfprintf(a,"{ ");
			for (z=0;z<((CxObArray*)m_pValue)->GetSize();z++)
			{
				((CDatabaseValue*)((CxObArray*)m_pValue)->GetAt(z))->DumpOutputFile(a);
				if (z < ((CxObArray*)m_pValue)->GetSize()-1)
					mfprintf(a,",");
				mfprintf(a," ");
			}
			mfprintf(a,"}");
			break;
		case 6:
			break;
	}
	BTOUT;
}


void CDatabaseValue::DumpValue()
{
	BTIN;
	int z;
	switch(m_iType)
	{
		case 0:
			mprintf("<NO VALUE>");
			break;
		case 1:
			TextColor(valuecolor);
			mprintf("%ld",m_pValue);
			TextColor(symbolcolor);
			break;
		case 2:
			TextColor(valuecolor);
			mprintf("%g",*(double*)m_pValue);
			TextColor(symbolcolor);
			break;
		case 3:
			TextColor(symbolcolor);
			mprintf("\"");
			TextColor(valuecolor);
			mprintf("%s",m_pValue);
			TextColor(symbolcolor);
			mprintf("\"");
			break;
		case 4:
			TextColor(valuecolor);
			if (m_pValue == 0)
				mprintf("FALSE");
					else mprintf("TRUE");
			TextColor(symbolcolor);
			break;
		case 5:
			TextColor(symbolcolor);
			mprintf("{ ");
			for (z=0;z<((CxObArray*)m_pValue)->GetSize();z++)
			{
				TextColor(valuecolor);
				((CDatabaseValue*)((CxObArray*)m_pValue)->GetAt(z))->DumpValue();
				TextColor(symbolcolor);
				if (z < ((CxObArray*)m_pValue)->GetSize()-1)
					mprintf(",");
				mprintf(" ");
			}
			mprintf("}");
			break;
		case 6:
			TextColor(valuecolor);
			mprintf("$...$");
			TextColor(symbolcolor);
			break;
	}
	BTOUT;
}


void CDatabaseValue::DumpType()
{
	BTIN;
	switch(m_iType)
	{
		case 0:
			mprintf("<INVALID>");
			break;
		case 1:
			mprintf("INT");
			break;
		case 2:
			mprintf("FLOAT");
			break;
		case 3:
			mprintf("STRING");
			break;
		case 4:
			mprintf("BOOL");
			break;
		case 5:
			mprintf("ARRAY");
			break;
		case 6:
			mprintf("TABLE");
			break;
	}
	BTOUT;
}



CDatabaseNode::CDatabaseNode()
{
	/*mprintf("CDatabaseNode::CDatabaseNode()\n");*/
	m_oaChildren.SetName("CDatabaseNode::m_oaChildren");
	m_oaValues.SetName("CDatabaseNode::m_oaValues");
}


CDatabaseNode::~CDatabaseNode()
{
	/*mprintf("CDatabaseNode::~CDatabaseNode()\n");*/
}


CDatabaseNode::CDatabaseNode(const char *s)
{
	BTIN;
	m_sName.Format(s);
	BTOUT;
}


CDatabaseNode* CDatabaseNode::REC_FindNode(CxString s)
{
	BTIN;
	int i, z;
		
#ifdef DEBUG_DATABASE		
	mprintf("CDatabaseNode::FindNode(CxString): Searching for \"%s\" in \"%s\".\n",(const char*)s,(const char*)m_sName);		
#endif		
		
	i = s.FindFirst('/');
	if (i == -1)
		i = s.GetLength();
		
	CxString t = s.Left(i);
		
	for (z=0;z<m_oaChildren.GetSize();z++)
	{
		if (((CDatabaseNode*)m_oaChildren[z])->m_sName == t)
		{
#ifdef DEBUG_DATABASE		
			mprintf("  Node \"%s\" found.\n",(const char*)t);
#endif		
			if (i != s.GetLength())
			{
				BTOUT;
				return ((CDatabaseNode*)m_oaChildren[z])->REC_FindNode(s.Mid(i+1));
			} else
			{
				BTOUT;
				return (CDatabaseNode*)m_oaChildren[z];
			}
		}
	}
#ifdef DEBUG_DATABASE		
	mprintf("  Node \"%s\" not found in \"%s\".\n",(const char*)t,(const char*)m_sName);
#endif		
	if (i != s.GetLength())
	{
		BTOUT;
		return NULL;
	}
	BTOUT;
	return NULL;
}

	
CDatabaseNode* CDatabaseNode::REC_FindAddNode(CxString s)
{
	BTIN;
	int i, z;
		
#ifdef DEBUG_DATABASE		
	mprintf("CDatabaseNode::FindAddNode(CxString): Searching for \"%s\" in \"%s\".\n",(const char*)s,(const char*)m_sName);		
#endif		
		
	i = s.FindFirst('/');
	if (i == -1)
		i = s.GetLength();
	
	CxString t = s.Left(i);
	
	for (z=0;z<m_oaChildren.GetSize();z++)
		if (((CDatabaseNode*)m_oaChildren[z])->m_sName == t)
	{
#ifdef DEBUG_DATABASE		
		mprintf("  Node \"%s\" found.\n",(const char*)t);
#endif		
		if (i != s.GetLength())
		{
			BTOUT;
			return ((CDatabaseNode*)m_oaChildren[z])->REC_FindAddNode(s.Mid(i+1));
		} else
		{
			BTOUT;
			return (CDatabaseNode*)m_oaChildren[z];
		}
	}
#ifdef DEBUG_DATABASE		
	mprintf("  Adding Node \"%s\" to \"%s\".\n",(const char*)t,(const char*)m_sName);
#endif	
	
	CDatabaseNode *p;
	try { p = new CDatabaseNode(t); } catch(...) { p = NULL; }
	if (p == NULL) NewException((double)sizeof(CDatabaseNode),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	p->m_pParent = this;
	m_oaChildren.Add(p);
	if (i != s.GetLength())
	{
		BTOUT;
		return p->REC_FindAddNode(s.Mid(i+1));
	} else
	{
		BTOUT;
		return p;
	}
}


CDatabaseValue* CDatabaseNode::FindValue(CxString s)
{
	BTIN;
	int z;
	
	for (z=0;z<m_oaValues.GetSize();z++)
		if (s == ((CDatabaseValue*)m_oaValues[z])->m_sName)
		{
			BTOUT;
			return (CDatabaseValue*)m_oaValues[z];
		}
	BTOUT;
	return NULL;
}


void CDatabaseNode::REC_DumpTree(int depth, unsigned long bitmask, bool last)
{
	BTIN;
	int z, z2;
	
	TextColor(treecolor);
	mprintf("  ");
	for (z=0;z<depth-1;z++)
	{
		if ((bitmask & (unsigned long)pow(2.0,z)) != 0)
//			mprintf("|   ");
			mprintf("%s  ",TreeElement("|"));
				else mprintf("   ");
	}
	if (depth != 0)
	{
		if (last)
//			mprintf("\\---");
			mprintf("%s",TreeElement("\\--"));
//				else mprintf("|---");
				else mprintf("%s",TreeElement("|--"));
	}
	TextColor(nodecolor);
	mprintf("[%s]\n",(const char*)m_sName);
	
	if ((m_oaValues.GetSize() != 0) || (m_oaChildren.GetSize() != 0))
	{
		TextColor(treecolor);
		mprintf("  ");
		for (z2=0;z2<depth;z2++)
			if ((bitmask & (unsigned long)pow(2.0,z2)) != 0)
//					mprintf("|   ");
				mprintf("%s  ",TreeElement("|"));
					else mprintf("   ");
//			mprintf("|\n");
		mprintf("%s\n",TreeElement("|"));
	}
	
	for (z=0;z<m_oaValues.GetSize();z++)
	{
		TextColor(treecolor);
		mprintf("  ");
		for (z2=0;z2<depth;z2++)
			if ((bitmask & (unsigned long)pow(2.0,z2)) != 0)
//					mprintf("|   ");
				mprintf("%s  ",TreeElement("|"));
					else mprintf("   ");
		if ((z+1 == m_oaValues.GetSize()) && (m_oaChildren.GetSize() == 0))
//				mprintf("\\---");
			mprintf("%s",TreeElement("\\--"));
//					else mprintf("|---");
				else mprintf("%s",TreeElement("|--"));
		TextColor(namecolor);
		mprintf("%s",(const char*)((CDatabaseValue*)m_oaValues[z])->m_sName);
//			if (((CDatabaseValue*)m_oaValues[z])->m_iType != 4)
		{
			TextColor(symbolcolor);
			mprintf(" = ");
			((CDatabaseValue*)m_oaValues[z])->DumpValue();
		}
		TextColor(symbolcolor);
		mprintf(" (");
		TextColor(typecolor);
		((CDatabaseValue*)m_oaValues[z])->DumpType();
		TextColor(symbolcolor);
		mprintf(")\n");
		if (((CDatabaseValue*)m_oaValues[z])->m_iType == 6)
		{
			((CDatabaseTable*)((CDatabaseValue*)m_oaValues[z])->m_pValue)->DumpTreeTable(depth,bitmask,(z+1 == m_oaValues.GetSize()) && (m_oaChildren.GetSize() == 0));
			
			if (z+1 < m_oaValues.GetSize())
			{
				TextColor(treecolor);
				mprintf("  ");
				for (z2=0;z2<depth;z2++)
					if ((bitmask & (unsigned long)pow(2.0,z2)) != 0)
//					mprintf("|   ");
						mprintf("%s  ",TreeElement("|"));
				else mprintf("   ");
//			mprintf("|\n");
				mprintf("%s\n",TreeElement("|"));
			}
		}
		if (z+1 == m_oaValues.GetSize())
		{
			TextColor(treecolor);
			mprintf("  ");
			for (z2=0;z2<depth;z2++)
				if ((bitmask & (unsigned long)pow(2.0,z2)) != 0)
//					mprintf("|   ");
					mprintf("%s  ",TreeElement("|"));
						else mprintf("   ");
			if (m_oaChildren.GetSize() != 0)
//				mprintf("|\n");
				mprintf("%s\n",TreeElement("|"));
					else mprintf("\n");
		}
	}
	for (z=0;z<m_oaChildren.GetSize();z++)
	{
		if (z+1 == m_oaChildren.GetSize())
			bitmask -= (unsigned long)pow(2.0,depth);
		((CDatabaseNode*)m_oaChildren[z])->REC_DumpTree(depth+1,bitmask,z+1==m_oaChildren.GetSize());
	}
	if ((m_oaChildren.GetSize() == 0) && last)
	{
		TextColor(treecolor);
		mprintf("  ");
		for (z2=0;z2<depth;z2++)
			if ((bitmask & (unsigned long)pow(2.0,z2)) != 0)
//				mprintf("|   ");
				mprintf("%s  ",TreeElement("|"));
					else mprintf("   ");
		mprintf("\n");
	}
	TextNormal();
	BTOUT;
}


void CDatabaseNode::REC_WriteOutputFile(FILE *a, int depth)
{
	BTIN;
	int z, z2;
	
	if (depth != 0)
		mfprintf(a,"\n");

	for (z=0;z<depth-1;z++)
		mfprintf(a,"  ");

	if (depth != 0)
		mfprintf(a,"&%s\n",(const char*)m_sName);
	
	for (z=0;z<m_oaValues.GetSize();z++)
	{
		for (z2=0;z2<depth;z2++)
			mfprintf(a,"  ");
		mfprintf(a,"%s",(const char*)((CDatabaseValue*)m_oaValues[z])->m_sName);
//		if (((CDatabaseValue*)m_oaValues[z])->m_iType != 4)
		{
			mfprintf(a," = ");
			if (((CDatabaseValue*)m_oaValues[z])->m_iType == 6)
			{
				mfprintf(a,"$\n");
				((CDatabaseTable*)((CDatabaseValue*)m_oaValues[z])->m_pValue)->DumpOutputFile(a,depth+1);
				mfprintf(a,"$");
			} else ((CDatabaseValue*)m_oaValues[z])->DumpOutputFile(a);
		}
		mfprintf(a,"\n");
	}
	
	for (z=0;z<m_oaChildren.GetSize();z++)
		((CDatabaseNode*)m_oaChildren[z])->REC_WriteOutputFile(a,depth+1);
	
	if (m_oaChildren.GetSize() != 0)
		mfprintf(a,"\n");

	for (z=0;z<depth-1;z++)
		mfprintf(a,"  ");
	if (depth != 0)
		mfprintf(a,"&/%s\n",(const char*)m_sName);
	BTOUT;
}


CDatabase::CDatabase()
{ 
	BTIN;
	try { m_pRoot = new CDatabaseNode("ROOT"); } catch(...) { m_pRoot = NULL; }
	if (m_pRoot == NULL) NewException((double)sizeof(CDatabaseNode),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	BTOUT;
}


CDatabase::~CDatabase()
{
	BTIN;
	delete m_pRoot;
	BTOUT;
}


void CDatabase::AddNode(CxString s)
{
	BTIN;
//	mprintf("CDatabase::AddNode(CxString): \"%s\"\n",(const char*)s);
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddNode(CxString): Wrong format.\n");
		BTOUT;
		return;
	}
	m_pRoot->REC_FindAddNode(s.Mid(1));
	BTOUT;
}


void CDatabase::DumpTree()
{
	BTIN;
	m_pRoot->REC_DumpTree(0,0xFFFFFFFF,false);
	BTOUT;
}


int CDatabase::GetElementType(CxString s)
{
	BTIN;
	if (s[0] != '/')
	{
		mprintf("CDatabase::GetElementType(CxString): Wrong format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		BTOUT;
		return 0;
	}
	BTOUT;
	return n->FindValue(s.Mid(s.FindLast('/')+1))->m_iType;
}


unsigned long CDatabase::GetInt(CxString s)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetInt(CxString): Wrong format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetInt(CxString): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return 0;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetInt(CxString): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 1)
	{
		eprintf("CDatabase::GetInt(CxString): Value \"%s\" is not an integer value.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	BTOUT;
	return (unsigned long)v->m_pValue;
}	


double CDatabase::GetFloat(CxString s)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetFloat(CxString): Wrong format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetFloat(CxString): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return 0;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetFloat(CxString): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 2)
	{
		eprintf("CDatabase::GetFloat(CxString): Value \"%s\" ist not a float value.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	BTOUT;
	return *(double*)v->m_pValue;
}	


int CDatabase::GetArraySize(CxString s)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetArraySize(CxString): Wrong format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetArraySize(CxString): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return 0;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetArraySize(CxString): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::GetArraySize(CxString): Value \"%s\" ist not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	BTOUT;
	return ((CxObArray*)v->m_pValue)->GetSize();
}	


int CDatabase::GetArrayElementType(CxString s, int i)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetArrayElementType(CxString,int): Wrong Format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetArrayElementType(CxString,int): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return 0;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetArrayElementType(CxString,int): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::GetArrayElementType(CxString,int): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	BTOUT;
	return ((CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i))->m_iType;
}	


int CDatabase::GetArrayInt(CxString s, int i)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetArrayInt(CxString,int): Wrong format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetArrayInt(CxString,int): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return 0;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetArrayInt(CxString,int): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::GetArrayInt(CxString,int): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	v2 = (CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i);
	if (v2->m_iType != 1)
	{
		eprintf("CDatabase::GetArrayInt(CxString,int): Value \"%s\"(%d) is not an integer value.\n",(const char*)s.Mid(s.FindLast('/')+1),i);
		BTOUT;
		return 0;
	}
	BTOUT;
	return (long)v2->m_pValue;
}	


double CDatabase::GetArrayFloat(CxString s, int i)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetArrayFloat(CxString,int): Wrong format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetArrayFloat(CxString,int): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return 0;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetArrayFloat(CxString,int): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::GetArrayFloat(CxString,int): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	v2 = (CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i);
	if (v2->m_iType != 2)
	{
		eprintf("CDatabase::GetArrayFloat(CxString,int): Value \"%s\"(%d) is not a floating point value.\n",(const char*)s.Mid(s.FindLast('/')+1),i);
		BTOUT;
		return 0;
	}
	BTOUT;
	return *(double*)v2->m_pValue;
}

	
bool CDatabase::GetArrayBool(CxString s, int i)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetArrayBool(CxString,int): Wrong format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetArrayBool(CxString,int): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return 0;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetArrayBool(CxString,int): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::GetArrayBool(CxString,int): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	v2 = (CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i);
	if (v2->m_iType != 4)
	{
		eprintf("CDatabase::GetArrayBool(CxString,int): Value \"%s\"(%d) is not a boolean value.\n",(const char*)s.Mid(s.FindLast('/')+1),i);
		BTOUT;
		return 0;
	}
	BTOUT;
	return (v2->m_pValue!=NULL);
}	


const char* CDatabase::GetArrayString(CxString s, int i)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetArrayString(CxString,int): Wrong format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetArrayString(CxString,int): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return 0;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetArrayString(CxString,int): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::GetArrayString(CxString,int): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	v2 = (CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i);
	if (v2->m_iType != 3)
	{
		eprintf("CDatabase::GetArrayString(CxString,int): Value \"%s\"(%d) is not a character string.\n",(const char*)s.Mid(s.FindLast('/')+1),i);
		BTOUT;
		return 0;
	}
	BTOUT;
	return (char*)v2->m_pValue;
}	


void CDatabase::SetArrayInt(CxString s, int i, int i2)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::SetArrayInt(CxString,int,int): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::SetArrayInt(CxString,int,int): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::SetArrayInt(CxString,int,int): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::SetArrayInt(CxString,int,int): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	v2 = (CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i);
	if (v2->m_iType != 1)
	{
		eprintf("CDatabase::SetArrayInt(CxString,int,int): Value \"%s\"(%d) is not an integer value.\n",(const char*)s.Mid(s.FindLast('/')+1),i);
		BTOUT;
		return;
	}
	v2->m_pValue = (void*)((unsigned long)i2);
	BTOUT;
}	


void CDatabase::SetArrayFloat(CxString s, int i, double d)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::SetArrayFloat(CxString,int,double): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::SetArrayFloat(CxString,int,double): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::SetArrayFloat(CxString,int,double): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::SetArrayFloat(CxString,int,double): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	v2 = (CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i);
	if (v2->m_iType != 2)
	{
		eprintf("CDatabase::SetArrayFloat(CxString,int,double): Value \"%s\"(%d) is not a floating point value.\n",(const char*)s.Mid(s.FindLast('/')+1),i);
		BTOUT;
		return;
	}
	*(double*)v2->m_pValue = d;
	BTOUT;
}

	
void CDatabase::SetArrayBool(CxString s, int i, bool b)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetArrayBool(CxString,int,bool): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetArrayBool(CxString,int,bool): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetArrayBool(CxString,int,bool): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::GetArrayBool(CxString,int,bool): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	v2 = (CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i);
	if (v2->m_iType != 4)
	{
		eprintf("CDatabase::GetArrayBool(CxString,int,bool): Value \"%s\"(%d) is not a boolean value.\n",(const char*)s.Mid(s.FindLast('/')+1),i);
		BTOUT;
		return;
	}
	v2->m_pValue = (void*)((unsigned long)(b?1:0));
	BTOUT;
}	


void CDatabase::SetArrayString(CxString s, int i, CxString s2)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::SetArrayString(CxString,int,CxString): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::SetArrayString(CxString,int,CxString): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::SetArrayString(CxString,int,CxString): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::SetArrayString(CxString,int,CxString): Value \"%s\" is not an array.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	v2 = (CDatabaseValue*)((CxObArray*)v->m_pValue)->GetAt(i);
	if (v2->m_iType != 3)
	{
		eprintf("CDatabase::SetArrayString(CxString,int,CxString): Value \"%s\"(%d) is not a character string.\n",(const char*)s.Mid(s.FindLast('/')+1),i);
		BTOUT;
		return;
	}
	delete (char*)v2->m_pValue;

	char *buf;
	try { buf = new char[s2.GetLength()+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(s2.GetLength()+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(buf,(const char*)s2);
	v2->m_pValue = (void*)buf;
	BTOUT;
}	


const char* CDatabase::GetString(CxString s)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetString(CxString): Wrong format.\n");
		BTOUT;
		return NULL;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetString(CxString): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return NULL;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetString(CxString): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return 0;
	}
	if (v->m_iType != 3)
	{
		eprintf("CDatabase::GetString(CxString): Value \"%s\" is not a character string.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return NULL;
	}
	BTOUT;
	return (char*)v->m_pValue;
}	


bool CDatabase::GetBool(CxString s)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::GetBool(CxString): Wrong format.\n");
		BTOUT;
		return false;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::GetBool(CxString): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return false;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::GetBool(CxString): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return false;
	}
	if (v->m_iType != 4)
	{
		eprintf("CDatabase::GetBool(CxString): Value \"%s\" is not a boolean value.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return false;
	}
	BTOUT;
	return (v->m_pValue != NULL);
}	


void CDatabase::AddArray(CxString s)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddArray(CxString): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	
	n = m_pRoot->REC_FindAddNode((s.Mid(1)).Left(s.FindLast('/')-1));

	CDatabaseValue *v;
	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	v->m_sName.Format(s.Mid(s.FindLast('/')+1));
	v->m_iType = 5;

	try { v->m_pValue = (void*)new CxObArray("CDatabase::AddArray():v->m_pValue"); } catch(...) { v->m_pValue = NULL; }
	if (v->m_pValue == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	n->m_oaValues.Add(v);
	BTOUT;
}


void CDatabase::AddArrayBool(CxString s, bool b)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddArrayBool(CxString,bool): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;

	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::AddArrayBool(CxString,bool): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::AddArrayBool(CxString,bool): Value \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::AddArrayBool(CxString,bool): Value \"%s\" is not an array.\n",(const char*)s);
		BTOUT;
		return;
	}

	try { v2 = new CDatabaseValue(); } catch(...) { v2 = NULL; }
	if (v2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v2->m_iType = 4;
	v2->m_pValue = (void*)((unsigned long)(b?1:0));
	v2->m_sName.Format("");
	((CxObArray*)v->m_pValue)->Add(v2);
	BTOUT;
}


void CDatabase::AddArrayFloat(CxString s, double d)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddArrayFloat(CxString,double): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;

	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::AddArrayFloat(CxString,double): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::AddArrayFloat(CxString,double): Value \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::AddArrayFloat(CxString,double): Value \"%s\" is not an array.\n",(const char*)s);
		BTOUT;
		return;
	}

	try { v2 = new CDatabaseValue(); } catch(...) { v2 = NULL; }
	if (v2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v2->m_iType = 2;

	double *t;
	try { t = new double; } catch(...) { t = NULL; }
	if (t == NULL) NewException((double)sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*t = d;
	v2->m_pValue = (void*)t;
	v2->m_sName.Format("");
	((CxObArray*)v->m_pValue)->Add(v2);
	BTOUT;
}


void CDatabase::AddArrayInt(CxString s, int i)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddArrayInt(CxString,int): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;

	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::AddArrayInt(CxString,int): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::AddArrayInt(CxString,int): Value \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::AddArrayInt(CxString,int): Value \"%s\" is not an array.\n",(const char*)s);
		BTOUT;
		return;
	}

	try { v2 = new CDatabaseValue(); } catch(...) { v2 = NULL; }
	if (v2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v2->m_iType = 1;
	v2->m_pValue = (void*)((unsigned long)i);
	v2->m_sName.Format("");
	((CxObArray*)v->m_pValue)->Add(v2);
	BTOUT;
}


void CDatabase::AddArrayString(CxString s, CxString s2)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddArrayString(CxString,CxString): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v, *v2;

	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::AddArrayString(CxString,CxString): Node \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::AddArrayString(CxString,CxString): Value \"%s\" not found.\n",(const char*)s);
		BTOUT;
		return;
	}
	if (v->m_iType != 5)
	{
		eprintf("CDatabase::AddArrayString(CxString,CxString): Value \"%s\" is not an array.\n",(const char*)s);
		BTOUT;
		return;
	}

	try { v2 = new CDatabaseValue(); } catch(...) { v2 = NULL; }
	if (v2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v2->m_iType = 3;

	char *buf;
	try { buf = new char[s2.GetLength()+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(s2.GetLength()+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(buf,(const char*)s2);
	v2->m_pValue = (void*)buf;
	v2->m_sName.Format("");
	((CxObArray*)v->m_pValue)->Add(v2);
	BTOUT;
}


void CDatabase::AddInt(CxString s, int i)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddInt(CxString,int): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	
	n = m_pRoot->REC_FindAddNode((s.Mid(1)).Left(s.FindLast('/')-1));

	CDatabaseValue *v;
	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v->m_sName.Format(s.Mid(s.FindLast('/')+1));
	v->m_iType = 1;
	v->m_pValue = (void*)((unsigned long)i);
	n->m_oaValues.Add(v);
	BTOUT;
}


void CDatabase::AddFloat(CxString s, double d)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddFloat(CxString,double): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	
	n = m_pRoot->REC_FindAddNode((s.Mid(1)).Left(s.FindLast('/')-1));

	CDatabaseValue *v;
	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v->m_sName.Format(s.Mid(s.FindLast('/')+1));
	v->m_iType = 2;

	double *t;
	try { t = new double; } catch(...) { t = NULL; }
	if (t == NULL) NewException((double)sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	*t = d;
	v->m_pValue = (void*)t;
	n->m_oaValues.Add(v);
	BTOUT;
}


void CDatabase::AddString(CxString s, CxString s2)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddString(CxString,CxString): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	
	n = m_pRoot->REC_FindAddNode((s.Mid(1)).Left(s.FindLast('/')-1));

	CDatabaseValue *v;
	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	v->m_sName.Format(s.Mid(s.FindLast('/')+1));
	v->m_iType = 3;

	char *buf;
	try { buf = new char[s2.GetLength()+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(s2.GetLength()+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(buf,(const char*)s2);
	v->m_pValue = (void*)buf;
	n->m_oaValues.Add(v);
	BTOUT;
}


void CDatabase::AddBool(CxString s, bool val)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::AddBool(CxString): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	
	n = m_pRoot->REC_FindAddNode((s.Mid(1)).Left(s.FindLast('/')-1));

	CDatabaseValue *v;
	try { v = new CDatabaseValue(); } catch(...) { v = NULL; }
	if (v == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	v->m_sName.Format(s.Mid(s.FindLast('/')+1));
	v->m_iType = 4;
	v->m_pValue = (void*)((unsigned long)(val?1:0));
	n->m_oaValues.Add(v);
	BTOUT;
}


void CDatabase::SetInt(CxString s, int i)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::SetInt(CxString,int): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::SetInt(CxString,int): Node \"%s\" not found.\n",(const char*)(s.Mid(1)).Left(s.FindLast('/')-1));
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::SetInt(CxString,int): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	if (v->m_iType != 1)
	{
		eprintf("CDatabase::SetInt(CxString,int): Value \"%s\" is not an integer value.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	v->m_pValue = (void*)((unsigned long)i);
	BTOUT;
}


void CDatabase::SetFloat(CxString s, double d)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::SetFloat(CxString,double): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::SetFloat(CxString,double): Node \"%s\" not found.\n",(const char*)(s.Mid(1)).Left(s.FindLast('/')-1));
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::SetFloat(CxString,double): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	if (v->m_iType != 2)
	{
		eprintf("CDatabase::SetFloat(CxString,double): Value \"%s\" is not a floating point value.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	*(double*)v->m_pValue = d;
	BTOUT;
}


void CDatabase::SetString(CxString s, CxString s2)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::SetString(CxString,CxString): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::SetString(CxString,CxString): Node \"%s\" not found.\n",(const char*)(s.Mid(1)).Left(s.FindLast('/')-1));
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::SetString(CxString,CxString): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	if (v->m_iType != 3)
	{
		eprintf("CDatabase::SetString(CxString,CxString): Value \"%s\" is not a character string.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	delete[] (char*)v->m_pValue;

	char *buf;
	try { buf = new char[s2.GetLength()+1]; } catch(...) { buf = NULL; }
	if (buf == NULL) NewException((double)(s2.GetLength()+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(buf,(const char*)s2);
	v->m_pValue = (void*)buf;
	BTOUT;
}


void CDatabase::SetBool(CxString s, bool b)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::SetBool(CxString,bool): Wrong format.\n");
		BTOUT;
		return;
	}
	CDatabaseNode *n;
	CDatabaseValue *v;
	n = m_pRoot->REC_FindNode((s.Mid(1)).Left(s.FindLast('/')-1));
	if (n == NULL)
	{
		eprintf("CDatabase::SetBool(CxString,bool): Node \"%s\" not found.\n",(const char*)(s.Mid(1)).Left(s.FindLast('/')-1));
		BTOUT;
		return;
	}
	v = n->FindValue(s.Mid(s.FindLast('/')+1));
	if (v == NULL)
	{
		eprintf("CDatabase::SetBool(CxString,bool): Value \"%s\" not found.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	if (v->m_iType != 4)
	{
		eprintf("CDatabase::SetBool(CxString,bool): Value \"%s\" is not a boolean value.\n",(const char*)s.Mid(s.FindLast('/')+1));
		BTOUT;
		return;
	}
	v->m_pValue = (void*)((unsigned long)(b?1:0));
	BTOUT;
}


bool CDatabase::WriteOutputFile(CxString s)
{
	BTIN;
	FILE *a;
	//a = OpenFileWrite((const char*)s,true);
	a = fopen((const char*)s,"wt");
	if (a == NULL)
	{
		eprintf("Could not open File \"%s\" for writing.\nBe sure to have writing permission for this directory.\n",(const char*)s);
		abort();
	}
	m_pRoot->REC_WriteOutputFile(a,0);
	mfprintf(a,"\n");
	fclose(a);
	BTOUT;
	return true;
}


void CDatabase::ParseInputFile(const char *s)
{
	BTIN;
	FILE *a;
	char *p, *q=NULL, *tc;
	char buf[1024], tbuf[1024], t2buf[1024], buf2[1024];
//	CxString buf, tbuf, t2buf, buf2;
	CDatabaseNode *cnode, *tnode;
	CDatabaseValue *val=NULL, *val2;
	double *td;
	int tabx=0, taby=0;
	bool tstring, section, endsection, name, value, namegiven, eq, tdouble, tint, valuegiven, array, table, linevalue;
	
	cnode = m_pRoot;
	
	a = fopen(s,"rt");
	
	if (a == NULL)
	{
		eprintf("CDatabase::ParseInputFile(char*): Could not open file \"%s\".\n",s);
		BTOUT;
		return;
	}
	
	table = false;
	while (!feof(a))
	{
		if (!fgets(buf,256,a))
			continue;
		
		buf[strlen(buf)-1] = 0;
#ifdef DEBUG_DATABASE
		mprintf("*** \"%s\" ***\n",buf);
#endif
		
		p = buf;
		
		tstring = false;
		section = false;
		endsection = false;
		name = false;
		value = false;
		namegiven = false;
		eq = false;
		tdouble = false;
		tint = false;
		valuegiven = false;
		array = false;
		linevalue = false;
		
		while (true)
		{
			if (*p == '\\')
			{
				if (p > buf)
				{
					if (*(p-1) == '\\')
						goto _nolincont;
				} else if (*p != 0)
					if (*(p+1) == '\\')
						goto _nolincont;
#ifdef DEBUG_DATABASE
				mprintf("Zeilenumbruch!\np=\"%s\"\nq=\"%s\"\n",p,q);
#endif
				*p = 0;
				fgets(buf2,256,a);
				buf2[strlen(buf2)-1] = 0;
#ifdef DEBUG_DATABASE
				mprintf("buf=\"%s\"\nbuf2=\"%s\"\n",buf,buf2);
#endif
				strcat(buf,buf2);
#ifdef DEBUG_DATABASE
				mprintf("Neuer buf=\"%s\"\n",buf);
#endif
			}
_nolincont:
			if (array)
			{
				if ((!tstring) && (*p == '{'))
				{
#ifdef DEBUG_DATABASE
					mprintf("ARRAY ERROR: Nested Arrays not allowed.\n");
#endif
					array = false;
					goto _linedone;
				}
				if ((!tstring) && ((*p == ',') || (*p == '}') || (*p == ' ') || (*p == '\t')))
				{
					if ((*p == ',') || (*p == ' ') ||(*p == '\t'))
					{
						if (value)
						{
#ifdef DEBUG_DATABASE
							mprintf("ARRAY NEXT ELEMENT\n");
#endif
						} else
						{
							q = p+1;
							goto _chardone;
						}
					} else if (*p == '}')
					{
						array = false;
#ifdef DEBUG_DATABASE
						mprintf("ARRAY DONE\n");
#endif
						if (!value)
							goto _linedone;
					}
					if (valuegiven)
					{
						valuegiven = false;
						q = p+1;
						value = false;
						goto _chardone;
					}
					if (tint)
					{
						memcpy(tbuf,q,p-q);
						tbuf[p-q] = 0;

						try { val2 = new CDatabaseValue(); } catch(...) { val2 = NULL; }
						if (val2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

						val2->m_sName.Format(t2buf);
						val2->m_iType = 1;
						val2->m_pValue = (void*)((unsigned long)atoi(tbuf));
						((CxObArray*)val->m_pValue)->Add(val2);
#ifdef DEBUG_DATABASE
						mprintf("ARRAY INT \"%s\" = %d\n",tbuf,val2->m_pValue);
#endif
					} else if (tdouble)
					{
						memcpy(tbuf,q,p-q);
						tbuf[p-q] = 0;

						try { val2 = new CDatabaseValue(); } catch(...) { val2 = NULL; }
						if (val2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

						val2->m_sName.Format(t2buf);
						val2->m_iType = 2;

						try { td = new double; } catch(...) { td = NULL; }
						if (td == NULL) NewException((double)sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						*td = atof(tbuf);
						val2->m_pValue = (void*)td;
						((CxObArray*)val->m_pValue)->Add(val2);
#ifdef DEBUG_DATABASE
						mprintf("ARRAY DOUBLE  \"%s\" = %f\n",tbuf,*td);
#endif
					} else
					{
						memcpy(tbuf,q,p-q);
						tbuf[p-q] = 0;
						if (mystricmp(tbuf,"true")==0)
						{

							try { val2 = new CDatabaseValue(); } catch(...) { val2 = NULL; }
							if (val2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

							val2->m_sName.Format(t2buf);
							val2->m_iType = 4;
							val2->m_pValue = (void*)1;
							((CxObArray*)val->m_pValue)->Add(val2);
#ifdef DEBUG_DATABASE
							mprintf("ARRAY BOOL \"%s\" = TRUE\n",tbuf);
#endif
						} else if (mystricmp(tbuf,"false")==0)
						{

							try { val2 = new CDatabaseValue(); } catch(...) { val2 = NULL; }
							if (val2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

							val2->m_sName.Format(t2buf);
							val2->m_iType = 4;
							val2->m_pValue = (void*)0;
							((CxObArray*)val->m_pValue)->Add(val2);
#ifdef DEBUG_DATABASE
							mprintf("ARRAY BOOL \"%s\" = FALSE\n",tbuf);
#endif
						} else eprintf("CDatabase::ParseInputFile(char*): Unrecognized array element \"%s\".\n",tbuf);
					}
					tint = false;
					tdouble = false;
					value = false;
					q = p+1;
					if (!array)
						goto _linedone;
							else goto _chardone;
				}
				
				if ((!tstring) && (*p == '"'))
				{
					tstring = true;
					q = p+1;
#ifdef DEBUG_DATABASE
					mprintf("ARRAY STRING ANFANG\n");
#endif
					goto _chardone;
				}
				if (tstring && (*p == '"'))
				{
					tstring = false;
					valuegiven = true;
					value = true;
					memcpy(tbuf,q,p-q);
					tbuf[p-q] = 0;
#ifdef DEBUG_DATABASE
					mprintf("ARRAY STRING ENDE: \"%s\".\n",tbuf);
#endif

					try { val2 = new CDatabaseValue(); } catch(...) { val2 = NULL; }
					if (val2 == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					val2->m_sName.Format(t2buf);
					val2->m_iType = 3;

					try { tc = new char[strlen(tbuf)+1]; } catch(...) { tc = NULL; }
					if (tc == NULL) NewException((double)(strlen(tbuf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					strcpy(tc,tbuf);
					RemoveDoubleBackslash(tc);
					val2->m_pValue = (void*)tc;
					((CxObArray*)val->m_pValue)->Add(val2);
					goto _chardone;
				}
					
				if ((!tstring) && (*p >= '0') && (*p <= '9') && (!tdouble))
				{
#ifdef DEBUG_DATABASE
					mprintf("ARRAY INT ANFANG\n");
#endif
					tint = true;
					value = true;
					goto _chardone;
				}
				
				if ((!tstring) && (*p == '.'))
				{
#ifdef DEBUG_DATABASE
					mprintf("ARRAY DOUBLE ANFANG\n");
#endif
					tint = false;
					tdouble = true;
					value = true;
					goto _chardone;
				}
				if (!value)
				{
#ifdef DEBUG_DATABASE
					mprintf("ARRAY BOOL ANFANG\n");
#endif
					value = true;
				}
				goto _chardone;
			} // ARRAY ENDE 
			
			if (table)
			{
				if ((!tstring) && (*p == '$'))
				{
					if (!value)
					{
#ifdef DEBUG_DATABASE
						mprintf("TABLE ENDE\n");
#endif
						table = false;
#ifdef DEBUG_DATABASE
						mprintf("\nHier die Tabelle:\n");
						((CDatabaseTable*)val->m_pValue)->DumpTable();
						mprintf("**ENDE**\n");
#endif
						goto _linedone;
					}
				}
				if ((*p == '\r') || (*p == '\n') || (*p == 0))
				{
					if (!linevalue)
					{
#ifdef DEBUG_DATABASE
						mprintf("TABLE EMPTY LINE\n");
#endif
						goto _linedone;
					} else
					{
#ifdef DEBUG_DATABASE
						mprintf("TABLE LINE END\n");
#endif
						if (!value)
						{
							tabx = 0;
							taby++;
							goto _linedone;
						}
					}
				}
				if ((!tstring) && ((*p == ',') || (*p == ' ') || (*p == '\t') || (*p == 0) || (*p == '\r') || (*p == '\n') || (*p == '$')))
				{
					if ((*p == ',') || (*p == ' ') || (*p == '\t'))
					{
						if (value)
						{
#ifdef DEBUG_DATABASE
							mprintf("TABLE NEXT ELEMENT\n");
#endif
						} else
						{
							q = p+1;
							goto _chardone;
						}
					} 
					if (valuegiven)
					{
						valuegiven = false;
						q = p+1;
						value = false;
						goto _chardone;
					}
					if (tint)
					{
						memcpy(tbuf,q,p-q);
						tbuf[p-q] = 0;
						((CDatabaseTable*)val->m_pValue)->AddInt(tabx,taby,atoi(tbuf));
#ifdef DEBUG_DATABASE
						mprintf("TABLE INT \"%s\" = %d\n",tbuf,atoi(tbuf));
#endif
						tabx++;
					} else if (tdouble)
					{
						memcpy(tbuf,q,p-q);
						tbuf[p-q] = 0;
						((CDatabaseTable*)val->m_pValue)->AddDouble(tabx,taby,atof(tbuf));
#ifdef DEBUG_DATABASE
						mprintf("TABLE DOUBLE  \"%s\" = %f\n",tbuf,atof(tbuf));
#endif
						tabx++;
					} else
					{
						memcpy(tbuf,q,p-q);
						tbuf[p-q] = 0;
						if (mystricmp(tbuf,"true")==0)
						{
							((CDatabaseTable*)val->m_pValue)->AddBool(tabx,taby,true);
#ifdef DEBUG_DATABASE
							mprintf("TABLE BOOL \"%s\" = TRUE\n",tbuf);
#endif
							tabx++;
						} else if (mystricmp(tbuf,"false")==0)
						{
							((CDatabaseTable*)val->m_pValue)->AddBool(tabx,taby,false);
#ifdef DEBUG_DATABASE
							mprintf("TABLE BOOL \"%s\" = FALSE\n",tbuf);
#endif
							tabx++;
						} else if (strcmp(tbuf,"-")==0)
						{
#ifdef DEBUG_DATABASE
							mprintf("TABLE PLACEHOLDER\n");
#endif
							tabx++;
						} else eprintf("CDatabase::ParseInputFile(char*): Unrecognized table element \"%s\".\n",tbuf);
					}
					if (*p == '$')
					{
#ifdef DEBUG_DATABASE
						mprintf("TABLE END\n");
#endif
						table = false;
					}
					tint = false;
					tdouble = false;
					value = false;
					q = p+1;
					if (!table)
						goto _linedone;
					else goto _chardone;
				}
				if ((!tstring) && (*p == '"'))
				{
					tstring = true;
					linevalue = true;
					q = p+1;
#ifdef DEBUG_DATABASE
					mprintf("TABLE STRING ANFANG\n");
#endif
					goto _chardone;
				}
				if (tstring && (*p == '"'))
				{
					tstring = false;
					valuegiven = true;
					value = true;
					memcpy(tbuf,q,p-q);
					tbuf[p-q] = 0;
#ifdef DEBUG_DATABASE
					mprintf("TABLE STRING ENDE: \"%s\".\n",tbuf);
#endif

					try { tc = new char[strlen(tbuf)+1]; } catch(...) { tc = NULL; }
					if (tc == NULL) NewException((double)(strlen(tbuf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					strcpy(tc,tbuf);
					RemoveDoubleBackslash(tc);
					((CDatabaseTable*)val->m_pValue)->AddString(tabx,taby,tc);
					delete[] tc;
					tabx++;
					goto _chardone;
				}
					
				if ((!tstring) && (*p >= '0') && (*p <= '9') && (!tdouble))
				{
#ifdef DEBUG_DATABASE
					mprintf("TABLE INT ANFANG\n");
#endif
					tint = true;
					value = true;
					linevalue = true;
					goto _chardone;
				}
				
				if ((!tstring) && (*p == '.'))
				{
#ifdef DEBUG_DATABASE
					mprintf("TABLE DOUBLE ANFANG\n");
#endif
					tint = false;
					tdouble = true;
					value = true;
					linevalue = true;
					goto _chardone;
				}
				if (!value)
				{
#ifdef DEBUG_DATABASE
					mprintf("TABLE BOOL ANFANG\n");
#endif
					value = true;
					linevalue = true;
				}
				goto _chardone;
			} // TABLE ENDE
	
			if ((tstring && (*p != '"')) || (*p == ' ') || (*p == '\t') || (*p == '#') || ((p > buf) && (*p == '/') && (*(p-1) == '/')) || (*p == '\r') || (*p == '\n') || (*p == 0))
			{
//				mprintf("CHARDONE.\n",*p);
				if (endsection)
				{
					memcpy(tbuf,q,p-q);
					tbuf[p-q] = 0;
#ifdef DEBUG_DATABASE
					mprintf("ENDSECTION \"%s\"\n",tbuf);
#endif
					cnode = cnode->m_pParent;
				} else if (section)
				{
					memcpy(tbuf,q,p-q);
					tbuf[p-q] = 0;
#ifdef DEBUG_DATABASE
					mprintf("SECTION \"%s\"\n",tbuf);
#endif
					if ((tnode=cnode->REC_FindNode(tbuf)) != NULL)
					{
#ifdef DEBUG_DATABASE
						mprintf("SECTION already exists.\n",tbuf);
#endif
					} else
					{
						try { tnode = new CDatabaseNode(); } catch(...) { tnode = NULL; }
						if (tnode == NULL) NewException((double)sizeof(CDatabaseNode),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
#ifdef DEBUG_DATABASE
						mprintf("SECTION created.\n",tbuf);
#endif

						tnode->m_pParent = cnode;
						tnode->m_sName.Format(tbuf);
						cnode->m_oaChildren.Add(tnode);
					}
					cnode = tnode;
				} else if (name)
				{
					name = false;
					namegiven = true;
					memcpy(tbuf,q,p-q);
					tbuf[p-q] = 0;
					strcpy(t2buf,tbuf);
#ifdef DEBUG_DATABASE
					mprintf("NAME \"%s\".\n",tbuf);
#endif
				} else if (value)
				{
					value = false;
					valuegiven = true;
					memcpy(tbuf,q,p-q);
					tbuf[p-q] = 0;
#ifdef DEBUG_DATABASE
					mprintf("VALUE: ");
#endif
					if (tint)
					{
						if ((val=cnode->FindValue(t2buf))!=NULL)
						{
#ifdef DEBUG_DATABASE
							mprintf("VALUE already exists.\n");
#endif
							if (val->m_iType != 1)
							{
								eprintf("CDatabase::ParseInputFile(char*): Redefinition of %s with different type.\n",t2buf);
								abort();
							}
						} else
						{
							try { val = new CDatabaseValue(); } catch(...) { val = NULL; }
							if (val == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

#ifdef DEBUG_DATABASE
							mprintf("VALUE created.\n");
#endif

							val->m_sName.Format(t2buf);
							val->m_iType = 1;
							cnode->m_oaValues.Add(val);
						}

						val->m_pValue = (void*)((unsigned long)atoi(tbuf));
#ifdef DEBUG_DATABASE
						mprintf("INT\n");
#endif
					} else if (tdouble)
					{
						if ((val=cnode->FindValue(t2buf))!=NULL)
						{
#ifdef DEBUG_DATABASE
							mprintf("VALUE already existent.\n");
#endif
							if (val->m_iType != 2)
							{
								eprintf("CDatabase::ParseInputFile(char*): Redefinition of %s with different type.\n",t2buf);
								abort();
							}
						} else
						{
							try { val = new CDatabaseValue(); } catch(...) { val = NULL; }
							if (val == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

#ifdef DEBUG_DATABASE
							mprintf("VALUE created.\n");
#endif

							val->m_sName.Format(t2buf);
							val->m_iType = 2;
							cnode->m_oaValues.Add(val);
						}

						try { td = new double; } catch(...) { td = NULL; }
						if (td == NULL) NewException((double)sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						*td = atof(tbuf);
						val->m_pValue = (void*)td;
#ifdef DEBUG_DATABASE
						mprintf("DOUBLE\n");
#endif
					} else
					{
						if ((val=cnode->FindValue(t2buf))!=NULL)
						{
#ifdef DEBUG_DATABASE
							mprintf("VALUE already existent.\n");
#endif
							if (val->m_iType != 4)
							{
								eprintf("CDatabase::ParseInputFile(char*): Redefinition of %s with different type.\n",t2buf);
								abort();
							}
						} else
						{
							try { val = new CDatabaseValue(); } catch(...) { val = NULL; }
							if (val == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

#ifdef DEBUG_DATABASE
							mprintf("VALUE created.\n");
#endif

							val->m_sName.Format(t2buf);
							val->m_iType = 4;
							cnode->m_oaValues.Add(val);
						}

						if (mystricmp(tbuf,"true")==0)
							val->m_pValue = (void*)1;
						else if (mystricmp(tbuf,"false")==0)
							val->m_pValue = (void*)0;
						else 
						{
							eprintf("CDatabase::ParseInputFile(char*): Unrecognized bool value \"%s\".\n",buf);
							goto _chardone;
						}
#ifdef DEBUG_DATABASE
						mprintf("BOOL\n");
#endif
					}
#ifdef DEBUG_DATABASE
					mprintf(" \"%s\".\n",tbuf);
#endif
				} 
				goto _chardone;
			}
			
			if (tstring && (*p == '"'))
			{
				tstring = false;
				valuegiven = true;
				memcpy(tbuf,q,p-q);
				tbuf[p-q] = 0;
#ifdef DEBUG_DATABASE
				mprintf("VALUE: STRING \"%s\".\n",tbuf);
#endif
				if ((val=cnode->FindValue(t2buf))!=NULL)
				{
#ifdef DEBUG_DATABASE
					mprintf("VALUE already existent.\n");
#endif
					if (val->m_iType != 3)
					{
						eprintf("CDatabase::ParseInputFile(char*): Redefinition of %s with different type.\n",t2buf);
						abort();
					}
				} else
				{
					try { val = new CDatabaseValue(); } catch(...) { val = NULL; }
					if (val == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

#ifdef DEBUG_DATABASE
					mprintf("VALUE created.\n");
#endif

					val->m_sName.Format(t2buf);
					val->m_iType = 3;
					cnode->m_oaValues.Add(val);
				}

				try { tc = new char[strlen(tbuf)+1]; } catch(...) { tc = NULL; }
				if (tc == NULL) NewException((double)(strlen(tbuf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				strcpy(tc,tbuf);
				RemoveDoubleBackslash(tc);
				val->m_pValue = (void*)tc;
			
				goto _chardone;
			}
			
			if (namegiven && (!tstring) && (*p == '"'))
			{
				tstring = true;
				q = p+1;
//				mprintf("STRING\n");
				goto _chardone;
			}
			
			if (namegiven && (!tstring) && (*p == '{'))
			{
				array = true;
				q = p+1;
#ifdef DEBUG_DATABASE
				mprintf("ARRAY START\n");
#endif
				try { val = new CDatabaseValue(); } catch(...) { val = NULL; }
				if (val == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				val->m_sName.Format(t2buf);
				val->m_iType = 5;

				try { val->m_pValue = new CxObArray("CDatabase::ParseInputFile():val->m_pValue"); } catch(...) { val->m_pValue = NULL; }
				if (val->m_pValue == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				cnode->m_oaValues.Add(val);
				goto _chardone;
			}
			
			if (namegiven && (!tstring) && (*p == '$'))
			{
				table = true;
				tabx = 0;
				taby = 0;
				q = p+1;
#ifdef DEBUG_DATABASE
				mprintf("TABLE START\n");
#endif
				try { val = new CDatabaseValue(); } catch(...) { val = NULL; }
				if (val == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				val->m_sName.Format(t2buf);

				try { val->m_pValue = (void*) new CDatabaseTable(); } catch(...) { val->m_pValue = NULL; }
				if (val->m_pValue == NULL) NewException((double)sizeof(CDatabaseTable),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				val->m_iType = 6;
				cnode->m_oaValues.Add(val);
				goto _chardone;
			}
			
			if (*p == '&')
			{
#ifdef DEBUG_DATABASE
				mprintf("SECTION.\n");
#endif
				section = true;
				q = p+1;
				goto _chardone;
			}
			
			if ((*p == '/') && section)
			{
#ifdef DEBUG_DATABASE
				mprintf("ENDSECTION.\n");
#endif
				endsection = true;
				q++;
				goto _chardone;
			}
			
			if (*p == '=')
			{
				if (name)
				{
					name = false;
					namegiven = true;
					memcpy(tbuf,q,p-q);
					tbuf[p-q] = 0;
					strcpy(t2buf,tbuf);
#ifdef DEBUG_DATABASE
					mprintf("NAME \"%s\".\n",tbuf);
#endif
				}
				if (!namegiven)
				{
					eprintf("CDatabase::ParseInputFile(char*): Error: \"=\" without variable name.\n  (\"%s\")\n",buf);
					goto _linedone;
				}
				eq = true;
				goto _chardone;
			}
			
			if ((isalpha(*p) || (*p == '_')) && (!name) && (!namegiven) && (!section))
			{
//				mprintf("NAME\n");
				name = true;
				q = p;
				goto _chardone;
			}
			
			if (namegiven && eq && (!value))
			{
//				mprintf("VALUE\n");
				value = true;
				q = p;
			}
			
			if (namegiven && (!eq) && (!value) && (!name))
			{
#ifdef DEBUG_DATABASE
				mprintf("(1)STATEMENT \"%s\"\n",t2buf);
#endif

				try { val = new CDatabaseValue(); } catch(...) { val = NULL; }
				if (val == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				val->m_sName.Format(t2buf);
				val->m_pValue = (void*)1;
				val->m_iType = 4;
				cnode->m_oaValues.Add(val);
				q = p;
				name = true;
			}	
					
			if (value && (*p >= '0') && (*p <= '9') && (!tdouble))
			{
//				mprintf("INT. double=%d\n",tdouble);
				tint = true;
				goto _chardone;
			}
			
			if (value && (*p == '.'))
			{
//				mprintf("DOUBLE\n");
				tint = false;
				tdouble = true;
				goto _chardone;
			}
			
_chardone:
			if ((*p == '#') || ((p > buf) && (*p == '/') && (*(p-1) == '/')) || (*p == '\r') || (*p == '\n') || (*p == 0))
			{
//				mprintf("LINEDONE.\n",*p);
				if (namegiven && (!valuegiven))
				{
#ifdef DEBUG_DATABASE
					mprintf("(2)STATEMENT \"%s\"\n",t2buf);
#endif

					try { val = new CDatabaseValue(); } catch(...) { val = NULL; }
					if (val == NULL) NewException((double)sizeof(CDatabaseValue),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					val->m_sName.Format(t2buf);
					val->m_pValue = (void*)1;
					val->m_iType = 4;
					cnode->m_oaValues.Add(val);
				}
				goto _linedone;
			}
			p++;
		}
		
_linedone:;
	}
	
	fclose(a);
	BTOUT;
}


bool CDatabase::ExistNode(CxString s)
{
	BTIN;
	if (s[0] != '/')
	{
		eprintf("CDatabase::ExistNode(CxString): Wrong Format.\n");
		BTOUT;
		return 0;
	}
	CDatabaseNode *n;
	
	n = m_pRoot->REC_FindNode(s.Mid(1));
	if (n == NULL)
	{
		BTOUT;
		return false;
	}
	BTOUT;
	return true;
}	
