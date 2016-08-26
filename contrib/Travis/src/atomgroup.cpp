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

#include "atomgroup.h"
#include "moltools.h"
#include "travis.h"


CAtomGroup::CAtomGroup()
{
	m_iAtomGes = 0;
	m_pMolecule = NULL;
	m_sName = NULL;
	m_bAllAtoms = false;

	m_baAtomType.SetName("CAtomGroup::m_baAtomType");
	m_baRealAtomType.SetName("CAtomGroup::m_baRealAtomType");
	m_baAllAtoms.SetName("CAtomGroup::m_baAllAtoms");
	m_oaAtoms.SetName("CAtomGroup::m_oaAtoms");
}


CAtomGroup::~CAtomGroup()
{
	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}
}


void CAtomGroup::AddAtom(int atom, int num, bool all)
{
	BTIN;
	int z, z2;
	CxIntArray *atoms;

	for (z=0;z<m_baAtomType.GetSize();z++)
	{
		if (m_baAtomType[z] == atom)
		{
			atoms = (CxIntArray*)m_oaAtoms[z];

			for (z2=0;z2<atoms->GetSize();z2++)
			{
				if (atoms->GetAt(z2) == num)
				{
					eprintf("CAtomGroup::AddAtom(): Atom %s%d already included, ignoring this.\n",((CAtom*)g_oaAtoms[atom])->m_sName,num+1);
					BTOUT; 
					return;
				}
			}

			atoms->Add(num);
			m_iAtomGes++;
			if (all && (m_baAllAtoms[z]==0))
				m_baAllAtoms[z] = 1;
			BTOUT; 
			return;
		}
	}

	m_baAtomType.Add(atom);
	m_baRealAtomType.Add(m_pMolecule->m_baAtomIndex[atom]);

	try { atoms = new CxIntArray("CAtomGroup::AddAtom():atoms"); } catch(...) { atoms = NULL; }
	if (atoms == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	atoms->Add(num);
	m_oaAtoms.Add(atoms);
	m_baAllAtoms.Add(all?1:0);
	m_iAtomGes++;
	BTOUT; 
}


void CAtomGroup::AddAllAtoms(CMolecule *m, bool virt)
{
	BTIN;
	int z, z2;
	CxIntArray *atoms;

	m_pMolecule = m;
	for (z=0;z<m->m_baAtomIndex.GetSize()-(virt?0:1);z++) // Virtuellen Atome mitnehmen?
	{
		try { atoms = new CxIntArray("CAtomGroup::AddAllAtoms():atoms"); } catch(...) { atoms = NULL; }
		if (atoms == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		for (z2=0;z2<m->m_waAtomCount[z];z2++)
			atoms->Add(z2);
		m_iAtomGes += m->m_waAtomCount[z];
		m_oaAtoms.Add(atoms);
		m_baAtomType.Add(z);
		m_baRealAtomType.Add(m_pMolecule->m_baAtomIndex[z]);
		m_baAllAtoms.Add(1);
	}
	m_bAllAtoms = true;
	BuildName();
	BTOUT; 
}


bool CAtomGroup::ParseAtoms(CMolecule *mol, const char *s)
{
	BTIN;
	const char *p, *q;
	char buf[32];
//	CxString buf;
	int atom, la, i, i2, z;
	bool m, all;
	const char *allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890,-#_ ";

	if ((s[0] == '*') && (s[1] == 0))
	{
		AddAllAtoms(mol,false);
		return true;
	}

	Reset();
	m_pMolecule = mol;
	p = s;
	atom = -1;
	la = -1;
	i = -1;
	i2 = -1;
	m = false;
	all = false;
	while (*p != 0)
	{
		while (*p == ' ')
			p++;
		if (strchr(allowed,*p) == NULL)
		{
			eprintf("Error: Character \"%c\" not allowed.\n",*p);
			BTOUT; 
			return false;
		}
		q = p;
		if (isalpha(*q) || (*q == '_') || (*q == '#'))
		{
			if (m)
			{
				eprintf("Error: Only digit allowed after \"-\".\n",buf);
				BTOUT; 
				return false;
			}
			while (isalpha(*q) || (*q == '_') || (*q == '#'))
				q++;
			if (q-p >= 32)
			{
				eprintf("ParseAtoms(): Internal Error A (%d >= 32).\n",q-p);
				return false;
			}
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			atom = mol->FindAtomInMol(buf);
			if (atom == -1)
			{
				eprintf("Error: Atom \"%s\" not in the molecule.\n",buf);
				BTOUT; 
				return false;
			}
		} else if (isdigit(*q))
		{
			if (atom == -1)
				atom = la;
			if (atom == -1)
			{
				eprintf("Error: Number in the beginning not possible.\n");
				BTOUT; 
				return false;
			}
			while (isdigit(*q)) q++;
			if ((*q != '-') && (*q != ',') && (*q != ' ') && (*q != 0))
			{
				eprintf("Error: Only \",\" or \"-\" may follow after a number.\n");
				BTOUT; 
				return false;
			}
			if (q-p >= 32)
			{
				eprintf("ParseAtoms(): Internal Error B (%d >= 32).\n",q-p);
				return false;
			}
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			if (atoi(buf)-1 >= mol->m_waAtomCount[atom])
			{
				eprintf("Error: Only %d %s atoms in the molecule (requested: %d)\n",mol->m_waAtomCount[atom],((CAtom*)g_oaAtoms[mol->m_baAtomIndex[atom]])->m_sName,atoi(buf));
				BTOUT; 
				return false;
			}
			if (m)
			{
				i2 = atoi(buf)-1;
				if (i2 < i)
				{
					eprintf("Error: Invalid atom range, %d < %d.\n",i2+1,i+1);
					BTOUT; 
					return false;
				}
			} else i = atoi(buf)-1;
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
			if (atom == -1)
			{
				eprintf("Error: Comma without atom.\n");
				BTOUT; 
				return false;
			}
			if (i == -1)
			{
				i = 0;
				i2 = mol->m_waAtomCount[atom]-1;
				all = true;
			} else all = false;
			if (all)
			{
				for (z=i;z<=i2;z++)
					AddAtom(atom,z,true);
			} else if (m)
			{
				for (z=i;z<=i2;z++)
					AddAtom(atom,z,false);
			} else AddAtom(atom,i,false);
			la = atom;
			m = false;
			i = -1;
			i2 = -1;
			atom = -1;
			q++;
		}
//		_end:
		p = q;
	}
	if (atom != -1)
	{
		if (i == -1)
		{
			i = 0;
			i2 = mol->m_waAtomCount[atom]-1;
			all = true;
		} else all = false;
		if (all)
		{
			for (z=i;z<=i2;z++)
				AddAtom(atom,z,true);
		} else if (m)
		{
			for (z=i;z<=i2;z++)
				AddAtom(atom,z,false);
		} else AddAtom(atom,i,false);
	}
	if (m_iAtomGes == 0)
		return false;
	SortAtoms();
	BuildName();
	BTOUT; 
	return true;
}


void CAtomGroup::BuildName()
{
//	char tmp[4096], buf2[32];
	CxString tmp, buf2;
	CxIntArray *atoms;
	int z, z2, z3;

	tmp.sprintf("");
	for (z=0;z<m_baAtomType.GetSize();z++)
	{
		if (m_baAllAtoms[z])
		{
//			strcat(tmp,((CAtom*)g_oaAtoms[m_baRealAtomType[z]])->m_sName);
			tmp.strcat(((CAtom*)g_oaAtoms[m_baRealAtomType[z]])->m_sName);
		} else
		{
			atoms = (CxIntArray*)m_oaAtoms[z];
			for (z2=0;z2<atoms->GetSize();z2++)
			{
				z3 = z2;
				if (z2 < atoms->GetSize()-1)
				{
					while (atoms->GetAt(z3)+1 == atoms->GetAt(z3+1))
					{
						z3++;
						if (z3 >= atoms->GetSize()-1)
							break;
					}
//					if (z3 > z2)
//						z3--;
				}	
				if (z3 > z2)
				{
//					sprintf(buf2,"%s%d-%d",((CAtom*)g_oaAtoms[m_baRealAtomType[z]])->m_sName,atoms->GetAt(z2)+1,atoms->GetAt(z3)+1);
					buf2.sprintf("%s%d-%d",((CAtom*)g_oaAtoms[m_baRealAtomType[z]])->m_sName,atoms->GetAt(z2)+1,atoms->GetAt(z3)+1);
					z2 = z3;
				} else 
					buf2.sprintf("%s%d",((CAtom*)g_oaAtoms[m_baRealAtomType[z]])->m_sName,atoms->GetAt(z2)+1);
//					sprintf(buf2,"%s%d",((CAtom*)g_oaAtoms[m_baRealAtomType[z]])->m_sName,atoms->GetAt(z2)+1);

//				strcat(tmp,buf2);
				tmp.strcat(buf2);
			}
		}
	}

	try { m_sName = new char[strlen(tmp)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(m_sName,tmp);
}


void CAtomGroup::SortAtoms()
{
	CxIntArray *atoms;
	int z, z2, z3, i, j;

	for (z=0;z<m_baAtomType.GetSize();z++)
	{
		if (!m_baAllAtoms[z])
		{
			atoms = (CxIntArray*)m_oaAtoms[z];
			for (z2=0;z2<atoms->GetSize();z2++)
			{
				j = 999999;
				i = -1;
				for (z3=z2;z3<atoms->GetSize();z3++)
				{
					if (atoms->GetAt(z3) < j)
					{
						j = atoms->GetAt(z3);
						i = z3;
					}
				}
				j = atoms->GetAt(i);
				atoms->SetAt(i,atoms->GetAt(z2));
				atoms->SetAt(z2,j);
			}
		}
	}
}


void CAtomGroup::Reset()
{
	int z;

	m_iAtomGes = 0;
	m_baAllAtoms.RemoveAll();
	m_baAtomType.RemoveAll();
	m_baRealAtomType.RemoveAll();
	for (z=0;z<m_oaAtoms.GetSize();z++)
		delete (CxIntArray*)m_oaAtoms[z];
	m_oaAtoms.RemoveAll();
}


void CAtomGroup::CopyFrom(CAtomGroup *p)
{
	int z;
	CxIntArray *wa;

	m_baAllAtoms.CopyFrom(&p->m_baAllAtoms);
	m_baAtomType.CopyFrom(&p->m_baAtomType);
	m_baRealAtomType.CopyFrom(&p->m_baRealAtomType);
	m_iAtomGes = p->m_iAtomGes;
	m_pMolecule = p->m_pMolecule;
	m_bAllAtoms = p->m_bAllAtoms;
	if (p->m_sName != NULL)
	{
		try { m_sName = new char[strlen(p->m_sName)+1]; } catch(...) { m_sName = NULL; }
		if (m_sName == NULL) NewException((double)(strlen(p->m_sName)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		strcpy(m_sName,p->m_sName);
	}
	for (z=0;z<m_oaAtoms.GetSize();z++)
	{
		try { wa = new CxIntArray("CAtomGroup::CopyFrom():wa"); } catch(...) { wa = NULL; }
		if (wa == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		wa->CopyFrom((CxIntArray*)p->m_oaAtoms[z]);
		m_oaAtoms.Add(wa);
	}
}
