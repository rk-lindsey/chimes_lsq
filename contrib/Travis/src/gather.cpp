/*****************************************************************************
    TRAVIS - Trajectory Analyzer and Visualizer
    http://www.travis-analyzer.de/

    Copyright (c) 2009-2016 Martin Brehm
                  2012-2016 Martin Thomas

    This file written by Martin Brehm and Martin Thomas.

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


#include "travis.h"
#include "tools.h"
#include "database.h"
#include "statistics.h"
#include "maintools.h"
#include "dacf.h"
#include "interface.h"
#include "plproj.h"
#include "sdfmap.h"
#include "domain.h"
#include "conversion.h"

bool GatherInfos()
{
	BTIN;
//	char buf[32768], buf2[32768];
	CxString buf, buf2;
	float tf, tf2;
	FILE *a;
	bool tb=false, tb2;
	int z, z0, z2, z3, z4, z5, z6, z7;
	int ti, ti2, ti3;
	unsigned char ty, rty, atom;
	unsigned char ty2, rty2, atom2;
	CVirtualAtom *va;
	CObservation *o;
	CMolecule *mol;
	CSingleMolecule *sm, *sm2;
	CAtomGroup *ag;
	CxWordArray tempwa;
	CMolecule *m, *m2;
	CAtom *at;
	CElement *el;
	CxIntArray *pia, ia;
	CxVector3 veca(0.0f, 0.0f, 0.0f), vecb(0.0f, 0.0f, 0.0f), vecc(0.0f, 0.0f, 0.0f);

	g_bKeepUnfoldedCoords = false;
	g_iDoubleBoxFactor = 1;
	g_bBondACF = false;
	g_bCombined = false;
	g_bRefEnvFix = false;
	g_bTimeDiff = false;
	g_bRemoveCOM = false;
	g_bSaveCoordsUnchanged = false;
	g_bSDFUniform = false;
	a = NULL;
	ti2 = -1;

	mprintf(YELLOW,"\n*** Interactive query of settings ***\n\n");
//	SAVEPOS;


	ti = 0;
	if (g_bUnknownElements)
	{
		tb2 = false;
_unkstart:
		tb = false;
		mprintf(WHITE,"Unrecognized atom types: ");
		for (z=0;z<g_oaAtoms.GetSize();z++)
		{
			if (z == g_iVirtAtomType)
				continue;

			at = (CAtom*)g_oaAtoms[z];
			while (at->m_pMergedTo != NULL)
			{
//				mprintf("%s merged to %s\n",at->m_sName,at->m_pMergedTo->m_sName);
				at = at->m_pMergedTo;
			}
			if (at->m_pElement->m_fRadius != 0)
				continue;
			if (tb)
				mprintf(", ");
			tb = true;
			mprintf("%s",at->m_sName);
			ti++;
			ti2 = z;
		}
		mprintf("\n\n");

		if ((ti == 1) && (mystricmp("X",((CAtom*)g_oaAtoms[ti2])->m_sName) == 0))
		{  // Only Wannier centers unrecognized.
			mprintf("    Possibly, atom type \"X\" are Wannier centers. If so, do not assign data to them.\n\n");
			tb2 = false; 
		} else tb2 = true;

		tb = false;
		if (AskYesNo("    Do you want to assign atom data to them (y/n)? [%s] ",tb2,tb2?"yes":"no"))
		{
			if (tb2)
			{
				if (AskYesNo("    Automatically rename all atoms to corresponding elements (y/n)? [yes] ",true))
				{
					mprintf("\n");
					tb2 = true;
					tb = false;
					for (z=0;z<g_oaAtoms.GetSize();z++)
					{
						if (z == g_iVirtAtomType)
							continue;

						at = (CAtom*)g_oaAtoms[z];
						if (at->m_pElement->m_fRadius != 0)
							continue;

						buf2.SetBufSize(3);
						buf2(0) = at->m_sName[0];
						buf2(1) = at->m_sName[1];
						buf2(2) = 0;
						el = FindElement(buf2,true);
						if (el != NULL)
						{
							if (el->m_fRadius == 0)
								goto _unkX2;
						} else
						{
_unkX2:
							buf2(1) = 0;
							el = FindElement(buf2,true);
							if (el != NULL)
							{
								if (el->m_fRadius == 0)
									goto _unkX3;
							} else
							{
_unkX3:
//								strcpy(buf2,at->m_sName);
								buf2.strcpy(at->m_sName);
							}
						}
						if (el == NULL)
						{
							mprintf("    - Skipping %s... (have no guess for atom type)\n",at->m_sName);
							tb = true;
							continue;
						}

						xAddAtom(buf2);

						mprintf("    - Renaming %s to %s...\n",at->m_sName,(const char*)buf2);
						for (z2=0;z2<g_oaAtoms.GetSize();z2++)
						{
							if (((CAtom*)g_oaAtoms[z2])->m_pElement == el)
							{
								at = (CAtom*)g_oaAtoms[z2];
								el = at->m_pElement;
								goto _unkX1;
							}
						}
						eprintf("    Strange error ^^\n");
_unkX1:;
						for (z2=0;z2<(int)g_TimeStep.m_iGesAtomCount;z2++)
						{
//							strcpy(buf,(char*)g_TimeStep.m_paLabels[z2]);
							buf.strcpy((char*)g_TimeStep.m_paLabels[z2]);
							ReplaceDigits(&buf);
							if (mystricmp(buf,((CAtom*)g_oaAtoms[z])->m_sName) == 0)
							{
								at->m_iCount++;
								((CAtom*)g_oaAtoms[z])->m_iCount--;
							}
						}
//						mprintf("\nmerging %s to %s.\n",((CAtom*)g_oaAtoms[z])->m_sName,at->m_sName);
						((CAtom*)g_oaAtoms[z])->m_pMergedTo = at;
					}
					if (tb)
					{
						mprintf("\n");
						goto _unkstart;
					}
					tb = false;
					goto _unkdone;
				}
			}
			ti = g_oaAtoms.GetSize();
			for (z=0;z<ti;z++)
			{
				if (z == g_iVirtAtomType)
					continue;

				at = (CAtom*)g_oaAtoms[z];
				
				while (at->m_pMergedTo != NULL)
				{
//					mprintf("%s merged to %s\n",at->m_sName,at->m_pMergedTo->m_sName);
					at = at->m_pMergedTo;
				}

				if (at->m_pElement->m_fRadius != 0)
					continue;
				mprintf(WHITE,"\n  * Element %s\n\n",at->m_sName);

				buf2.SetBufSize(3);
				buf2(0) = at->m_sName[0];
				buf2(1) = at->m_sName[1];
				buf2(2) = 0;
				el = FindElement(buf2,true);
				if (el != NULL)
				{
					if (el->m_fRadius == 0)
						goto _unk2;
				} else
				{
_unk2:
					buf2(1) = 0;
					el = FindElement(buf2,true);
					if (el != NULL)
					{
						if (el->m_fRadius == 0)
							goto _unk3;
					} else
					{
_unk3:
//						strcpy(buf2,at->m_sName);
						buf2.strcpy(at->m_sName);
					}
				}

				if (AskYesNo("    Do you want to rename this element (y/n)? [yes] ",true))
				{
_unk5:
					AskString("    Enter new name for %s: [%s] ",&buf,buf2,at->m_sName,(const char*)buf2);
					if (strlen(buf) > 7)
					{
						eprintf("Atom labels may have maximum length of 7 characters.\n");
						goto _unk5;
					}
					if (ContainsDigit(buf))
					{
						eprintf("Digits in element labels not allowed.\n");
						goto _unk5;
					}
					el = FindElement(buf,true);
					if (el != NULL)
					{
						xAddAtom(buf);
						if (AskYesNo("    Element %s is known. Merge %s atoms into %s (y/n)? [yes] ",true,(const char*)buf,at->m_sName,(const char*)buf))
						{
							for (z2=0;z2<g_oaAtoms.GetSize();z2++)
							{
								if (((CAtom*)g_oaAtoms[z2])->m_pElement == el)
								{
									at = (CAtom*)g_oaAtoms[z2];
									el = at->m_pElement;
									goto _unk1;
								}
							}
							eprintf("    Strange error ^^\n");
_unk1:;
						}
					} else
					{
						mprintf("      Adding element %s...",(const char*)buf);

						try { at = new CAtom(); } catch(...) { at = NULL; }
						if (at == NULL) NewException((double)sizeof(CAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						at->m_iIndex = g_oaAtoms.GetSize();
						g_oaAtoms.Add(at);

						try { at->m_pElement = new CElement(); } catch(...) { at->m_pElement = NULL; }
						if (at->m_pElement == NULL) NewException((double)sizeof(CElement),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						el = at->m_pElement;
						sprintf(at->m_sName,"%s",(const char*)buf);
					}
					for (z2=0;z2<(int)g_TimeStep.m_iGesAtomCount;z2++)
					{
//						strcpy(buf,(char*)g_TimeStep.m_paLabels[z2]);
						buf.strcpy((char*)g_TimeStep.m_paLabels[z2]);
						ReplaceDigits(&buf);
						if (mystricmp(buf,((CAtom*)g_oaAtoms[z])->m_sName) == 0)
						{
							at->m_iCount++;
							((CAtom*)g_oaAtoms[z])->m_iCount--;
						}
					}
					((CAtom*)g_oaAtoms[z])->m_pMergedTo = at;
					if (el->m_fRadius != 0)
						continue;
				} // END IF RENAME

				if (AskYesNo("    Do you want to copy the atom parameters from another element (y/n)? [yes] ",true))
				{
_unk4:
					AskString("    Copy atom parameters from which element? [%s] ",&buf,buf2,(const char*)buf2);
					el = FindElement(buf,true);
					if (el == NULL)
					{
						eprintf("This element is not defined.\n");
						goto _unk4;
					}
					if (el->m_fRadius == 0)
					{
						if (!AskYesNo("    This element has a radius of 0 (will form no bonds). Use it anyway (y/n)? [yes] ",true))
							goto _unk4;
					}

					mprintf(WHITE,"\n    A covalent radius of 0 pm hinders atoms of this kind from forming any bonds.\n\n");

					at->m_pElement->m_fRadius = AskFloat("    Please enter covalent radius in pm: [%.2f] ",el->m_fRadius,el->m_fRadius);
					at->m_pElement->m_fMass = AskFloat("    Please enter atom mass in u: [%.2f] ",el->m_fMass,el->m_fMass);
	
					if (!tb && (fabs(at->m_pElement->m_fMass-el->m_fMass) < 0.01))
					{
						tb = true;
						mprintf(RED,"\n    Warning: ");
						mprintf("Elements %s and %s keep their different labels, but will have\n",at->m_sName,(const char*)buf);
						mprintf("             the same mass and therefore also the same atom codes. They will not\n");
						mprintf("             be distinguished in the molecule recognition!\n");
					}
//					mprintf("    Using mass=%.2fu, radius=%.2fpm, ord. number=%d.\n",at->m_pElement->m_fMass,el->m_fRadius,el->m_iOrd);
//					at->m_pElement->CopyData(el);
				} else
				{
					el = FindElement(buf2,true);
					mprintf(WHITE,"\n    A covalent radius of 0 pm hinders atoms of this kind from forming any bonds.\n\n");
					if (el != NULL)
					{
						mprintf("    (Default values are for %s)\n\n",(const char*)buf2);
						at->m_pElement->m_fRadius = AskFloat("    Please enter covalent radius in pm: [%.2f] ",el->m_fRadius,el->m_fRadius);
						at->m_pElement->m_fMass = AskFloat("    Please enter atom mass in u: [%.2f] ",el->m_fMass,el->m_fMass);
					} else
					{
						at->m_pElement->m_fRadius = AskFloat("    Please enter covalent radius in pm: [100] ",100);
						at->m_pElement->m_fMass = AskFloat_ND("    Please enter atom mass in u: ");
					}
				}
			}
_unkdone:
			mprintf("\n");
			mprintf(WHITE,"%d atoms in the system: ",g_iGesAtomCount);

			for (z=0;z<g_oaAtoms.GetSize();z++)
			{
				if (z == g_iVirtAtomType)
					continue;
				if (((CAtom*)g_oaAtoms[z])->m_pMergedTo != NULL)
					continue;
				mprintf("%dx %s",((CAtom*)g_oaAtoms[z])->m_iCount,((CAtom*)g_oaAtoms[z])->m_sName);
				if (((CAtom*)g_oaAtoms[z])->m_pMergedTo != NULL)
					mprintf(" (merged to %s)",((CAtom*)g_oaAtoms[z])->m_pMergedTo->m_sName);
				if (z < (int)g_oaAtoms.GetSize()-1)
					if (z+1 != g_iVirtAtomType)
						mprintf(", ");
			}
			mprintf("\n");
		} // END IF CHANGE UNKNOWN
		mprintf("\n");
	} // END IF UNKNOWN

	ia.RemoveAll();
	tb2 = false;
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		tb = false;
		for (z3=0;z3<ia.GetSize();z3++)
		{
			if (z == ia[z3])
				goto _next;
		}
		if (((CAtom*)g_oaAtoms[z])->m_pMergedTo != NULL)
			continue;
		if (z == g_iVirtAtomType)
			continue;
		for (z2=z+1;z2<g_oaAtoms.GetSize();z2++)
		{
			if (((CAtom*)g_oaAtoms[z2])->m_pMergedTo != NULL)
				continue;
			if (z2 == g_iVirtAtomType)
				continue;
			if (fabs(((CAtom*)g_oaAtoms[z])->m_pElement->m_fMass - ((CAtom*)g_oaAtoms[z2])->m_pElement->m_fMass) < 0.01)
			{
				if (!tb)
				{
					tb = true;
					tb2 = true;
					mprintf(RED,"\n    Warning: ");
					mprintf("The following atoms have the same mass and will therefore not be\n");
					mprintf("             distinguished in the molecule recognition:\n\n");
					mprintf(WHITE,"      Mass=%.2f: %s, %s",((CAtom*)g_oaAtoms[z])->m_pElement->m_fMass,((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z2])->m_sName);
					ia.Add(z);
					ia.Add(z2);
				} else
				{
					mprintf(WHITE,", %s",((CAtom*)g_oaAtoms[z2])->m_sName);
				}
			}
		}
		if (tb)
			mprintf("\n");
_next:;
	}
	if (tb2)
		mprintf("\n");

	mprintf("    The advanced mode includes many additional options which are quite powerful, yet possibly\n");
    mprintf("    weird or seldomly required. This includes support for non-orthorhombic simulation cells,\n");
	mprintf("    NPT ensemble, non-periodic systems, user-definded virtual atoms, pseudomolecules, ...\n\n");

	g_bAdvanced1 = AskYesNo("    Use the advanced mode until the analysis selection menu (y/n)? [no] ",false);

	mprintf("\n");

	if (g_bXYZComment6Numbers)
	{
		mprintf("    TRAVIS can read cell vectors/angles from the comment line of the XYZ file.\n");
		mprintf("    The format needs to be \"a b c alpha beta gamma\". Lengths in Angstrom / Angles in Degree.\n\n");
		mprintf("    The current comment line is \"%s\".\n\n",g_TimeStep.m_pComment);
		if (AskYesNo("    Are the six numbers in the trajectory's comment line cell geometry data in this format (y/n)? [yes] ",true))
		{
			g_bFoundNonOrtho = true;
			mprintf("\n    Extracting cell geometry...\n\n");
			ExtractXYZCellGeometry(g_TimeStep.m_pComment);
		} else
			g_bXYZComment6Numbers = false;
	}

	if ((g_fBoxX != 0) || g_bFoundNonOrtho)
	{
		mprintf("    Found cell geometry data in trajectory file:\n\n");
		if (g_bFoundNonOrtho)
		{
			DumpNonOrthoCellData();
		} else
		{
			mprintf("      X = %.2f pm,   Y = %.2f pm,   Z = %.2f pm\n\n",g_fBoxX,g_fBoxY,g_fBoxZ);
			
			float tf = GuessBoxSize();
			mprintf("    The overall box density is %.6f g/cm^3.\n\n",tf*tf*tf/g_fBoxX/g_fBoxY/g_fBoxZ);
		}

		mprintf("    Assuming a 3D-periodic cell.\n\n");

		if (AskYesNo("    Use these values (y) or enter different values (n)? [yes] ",true))
		{
			if (g_bFoundNonOrtho)
				g_bBoxNonOrtho = true;
			else
				g_bBoxNonOrtho = false;

			g_bPeriodic = true;
			g_bPeriodicX = true;
			g_bPeriodicY = true;
			g_bPeriodicZ = true;

			if (AskYesNo("    Update cell geometry in every time step (i.e., NPT ensemble) (y) or use fixed cell (n)? [yes] ",true))
				g_bNPT = true;
			else
				g_bNPT = false;

			goto _celldefined;
		} else
			g_bXYZComment6Numbers = false;
	}

	if (g_bAdvanced1)
	{
		do
		{
			mprintf("    Use periodic boundary conditions (0=no, x, y, z, xy, xz, xyz)? [xyz] ");
			inpprintf("! Use periodic boundary conditions (0=no, x, y, z, xy, xz, xyz)? [xyz]\n");
			myget(&buf);
		} while (!ParsePeriodic(buf));

//		if (!g_bPeriodic)
//		{
/*			if (AskYesNo("    Enter Box Size anyways (e.g. fo RDFs) (yes/no)? [yes] ",true))
			{
				goto _askbox;
			} else
			{
				g_fBoxX = GuessBoxSize();
				g_fBoxY = GuessBoxSize();
				g_fBoxZ = GuessBoxSize();
				mprintf("\nGuessed box size (for density 1.0 g/cm^3) is %.2f x %.2f x %.2f pm.\n\n",g_fBoxX,g_fBoxX,g_fBoxX);
			}*/
//			g_fBoxX = (g_TimeStep.m_vMax[0]-g_TimeStep.m_vMin[0])*10.0;
//			g_fBoxY = (g_TimeStep.m_vMax[1]-g_TimeStep.m_vMin[1])*10.0;
//			g_fBoxZ = (g_TimeStep.m_vMax[2]-g_TimeStep.m_vMin[2])*10.0;
//			mprintf("\n    Using cell vector of %.2f x %.2f x %.2f pm to \"fake\" non-periodic box.\n\n",g_fBoxX,g_fBoxY,g_fBoxZ);
//		}
	} else
		ParsePeriodic("xyz");

	if (g_bPeriodic)
	{
		tf = GuessBoxSize();
		mprintf("    A cell vector of %.2f pm would result in a density of 1.0 g/cm^3.\n\n",tf);

		if (g_bAdvanced1)
		{
			g_bNPT = AskYesNo("    Use time-dependent cell vector (NPT ensemble) (y/n)? [no] ",false);
			mprintf("\n");
		} else
			g_bNPT = false;

		if (g_bAdvanced1)
		{
			g_bBoxNonOrtho = AskYesNo("    Is the simulation cell non-orthorhombic (cell angles other than 90 deg) (y/n)? [no] ",false);
			mprintf("\n");
		} else
			g_bBoxNonOrtho = false;

		if (g_bNPT)
		{
			if (g_fBoxX != 0)
			{
				mprintf("    Your PDB file contains cell vector information.\n\n");
			} else
			{
				mprintf("    The trajectory does not contain cell vector information.\n");
				mprintf("    (supported: PDB file with \"CRYST1\" section, like written by GROMACS trjconv)\n");
				mprintf("    You have to supply a plain text file which contains the cell vectors for each step.\n");
				mprintf("    It should contain three space-separated floating point numbers (X, Y, Z vector) per line.\n");
				mprintf("    These numbers need to be IN ANGSTROM! (1 angstrom = 100 pm)\n\n");
_nptfileagain:
				AskString_ND("    Please enter cell vector text file name: ",&g_sNPTFile);
				if (!FileExist(g_sNPTFile))
				{
					eprintf("Could not open file \"%s\" for reading.\n",(const char*)g_sNPTFile);
					goto _nptfileagain;
				}
				g_fNPTFile = fopen(g_sNPTFile,"rt");
				g_TimeStep.ReadCellVector(g_fNPTFile);
				fclose(g_fNPTFile);
				mprintf("\n");
			}
			mprintf("    The initial cell vector is ( %.2f pm | %.2f pm | %.2f pm ).\n",g_fBoxX,g_fBoxY,g_fBoxZ);
			mprintf("    The initial box density is %.6f g/cm^3.\n\n",tf*tf*tf/g_fBoxX/g_fBoxY/g_fBoxZ);
		} else if (!g_bBoxNonOrtho)
		{
			if (g_fBoxX != 0)
			{
				mprintf("    Cell vector found in trajectory file:\n");
				mprintf("      X = %.2f pm,   Y = %.2f pm,   Z = %.2f pm\n\n",g_fBoxX,g_fBoxY,g_fBoxZ);
				if (AskYesNo("    Use these values (y/n)? [yes] ",true))
					goto _celldone;
			}
_askbox:
			if (!AskYesNo("    Are the 3 cell vectors of the same size (yes/no)? [yes] ",true))
			{
				if (g_bPeriodicX)
					g_fBoxX = AskFloat_ND("    Enter length of X cell vector in pm: ");
				if (g_bPeriodicY)
					g_fBoxY = AskFloat_ND("    Enter length of Y cell vector in pm: ");
				if (g_bPeriodicZ)
					g_fBoxZ = AskFloat_ND("    Enter length of Z cell vector in pm: ");
			} else
			{
				g_fBoxX = AskFloat_ND("    Enter length of cell vector in pm: ");
				g_fBoxY = g_fBoxX;
				g_fBoxZ = g_fBoxX;
			}
			if ((g_fBoxX <= 0) || (g_fBoxY <= 0) || (g_fBoxZ <= 0))
			{
				eprintf("Cell vector of length <= 0 is not allowed.\n");
				goto _askbox;
			}
_celldone:
			if (g_bPeriodicX && g_bPeriodicY && g_bPeriodicZ)
			{
				mprintf("\n    The box size is %.2f x %.2f x %.2f pm.\n",g_fBoxX,g_fBoxY,g_fBoxZ);
				mprintf("    The overall box density is %.6f g/cm^3.\n",tf*tf*tf/g_fBoxX/g_fBoxY/g_fBoxZ);
				if (tf*tf*tf/g_fBoxX/g_fBoxY/g_fBoxZ > 15.0f)
				{
					mprintf("\n");
					if (!AskYesNo("    The density of your box seems to be very high. Continue (y) or change input (n)? [n] ",false))
					{
						mprintf("\n");
						goto _askbox;
					}
					mprintf("\n");
				}
			}
		} else // Non-orthorhombic cell definition
		{
			if (AskYesNo("    Enter three cell vectors (y) or three edge lengths and angles (n)? [no] ",false))
			{
				mprintf("\n");

				g_mBoxFromOrtho(0,0) = AskFloat_ND("    Enter X component of cell vector A in pm: ");
				g_mBoxFromOrtho(0,1) = AskFloat_ND("    Enter Y component of cell vector A in pm: ");
				g_mBoxFromOrtho(0,2) = AskFloat_ND("    Enter Z component of cell vector A in pm: ");
				mprintf("\n");
				g_mBoxFromOrtho(1,0) = AskFloat_ND("    Enter X component of cell vector B in pm: ");
				g_mBoxFromOrtho(1,1) = AskFloat_ND("    Enter Y component of cell vector B in pm: ");
				g_mBoxFromOrtho(1,2) = AskFloat_ND("    Enter Z component of cell vector B in pm: ");
				mprintf("\n");
				g_mBoxFromOrtho(2,0) = AskFloat_ND("    Enter X component of cell vector C in pm: ");
				g_mBoxFromOrtho(2,1) = AskFloat_ND("    Enter Y component of cell vector C in pm: ");
				g_mBoxFromOrtho(2,2) = AskFloat_ND("    Enter Z component of cell vector C in pm: ");

				veca = CxVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
				vecb = CxVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
				vecc = CxVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

				g_fBoxAngleA = acosf(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
				g_fBoxAngleB = acosf(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
				g_fBoxAngleC = acosf(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;
			} else
			{
				mprintf("\n    Convention: \"Vector A\" is located on the X axis, \"vector B\" in the XY plane.\n\n");
				g_fBoxX = AskFloat_ND("    Enter length of cell vector A in pm: ");
				g_fBoxY = AskFloat_ND("    Enter length of cell vector B in pm: ");
				g_fBoxZ = AskFloat_ND("    Enter length of cell vector C in pm: ");
				mprintf("\n");
				g_fBoxAngleA = AskFloat("    Enter angle Alpha (between B and C) in deg: [90.0] ",90.0);
				g_fBoxAngleB = AskFloat("    Enter angle Beta  (between A and C) in deg: [90.0] ",90.0);
				g_fBoxAngleC = AskFloat("    Enter angle Gamma (between A and B) in deg: [90.0] ",90.0);
				mprintf("\n");
				mprintf("    Computing vectors...\n");

				g_mBoxFromOrtho(0,0) = g_fBoxX;
				g_mBoxFromOrtho(0,1) = 0;
				g_mBoxFromOrtho(0,2) = 0;

				g_mBoxFromOrtho(1,0) = g_fBoxY*cos(g_fBoxAngleC*Pi/180.0);
				g_mBoxFromOrtho(1,1) = g_fBoxY*sin(g_fBoxAngleC*Pi/180.0);
				g_mBoxFromOrtho(1,2) = 0;

				g_mBoxFromOrtho(2,0) = g_fBoxZ*cos(g_fBoxAngleB*Pi/180.0);
				g_mBoxFromOrtho(2,1) = (-(g_mBoxFromOrtho(1,0)*g_fBoxZ*cos(g_fBoxAngleB*Pi/180.0)) + g_fBoxY*g_fBoxZ*cos(g_fBoxAngleA*Pi/180.0))/g_mBoxFromOrtho(1,1);
				g_mBoxFromOrtho(2,2) = sqrtf(-((pow(g_fBoxZ,2)*(pow(g_mBoxFromOrtho(1,1),2)*(-1 + pow(cos(g_fBoxAngleB*Pi/180.0),2)) + pow(g_mBoxFromOrtho(1,0)*cos(g_fBoxAngleB*Pi/180.0) - g_fBoxY*cos(g_fBoxAngleA*Pi/180.0),2)))/pow(g_mBoxFromOrtho(1,1),2)));

				mprintf("    Recalculating angles from vectors...\n");

				veca = CxVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
				vecb = CxVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
				vecc = CxVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

				g_fBoxAngleA = acosf(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
				g_fBoxAngleB = acosf(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
				g_fBoxAngleC = acosf(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;
			}

			mprintf("    Computing inverse of transformation matrix...\n\n");
			g_mBoxToOrtho = CxMatrix3(g_mBoxFromOrtho);
			g_mBoxToOrtho.Invert();

			// Orthogonal bounding box
			g_fBoxX = fabsf(g_mBoxFromOrtho(0,0)) + fabsf(g_mBoxFromOrtho(1,0)) + fabsf(g_mBoxFromOrtho(2,0));
			g_fBoxY = fabsf(g_mBoxFromOrtho(0,1)) + fabsf(g_mBoxFromOrtho(1,1)) + fabsf(g_mBoxFromOrtho(2,1));
			g_fBoxZ = fabsf(g_mBoxFromOrtho(0,2)) + fabsf(g_mBoxFromOrtho(1,2)) + fabsf(g_mBoxFromOrtho(2,2));

			veca = CxVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
			vecb = CxVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
			vecc = CxVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

			// Minimal diameters
			g_fBoxMinDiamA = fabsf(DotP(veca,Normalize(CrossP(vecb,vecc))));
			g_fBoxMinDiamB = fabsf(DotP(vecb,Normalize(CrossP(veca,vecc))));
			g_fBoxMinDiamC = fabsf(DotP(vecc,Normalize(CrossP(veca,vecb))));
			g_fBoxMinDiam = MIN3(g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);

			DumpNonOrthoCellData();
		}
		
		if (g_iTrajFormat == 5) {
			mprintf("\n    Cube file has a resolution of %d x %d x %d.\n",g_TimeStep.m_pVolumetricData->m_iRes[0],g_TimeStep.m_pVolumetricData->m_iRes[1],g_TimeStep.m_pVolumetricData->m_iRes[2]);
			if (g_bBoxNonOrtho) {
				float ang = Angle_Deg(CxVector3(g_mCubeCell(0, 0), g_mCubeCell(0, 1), g_mCubeCell(0, 2)), veca);
				if (ang > 1.0f) {
					eprintf("The direction of cell vector A does not match the grid data in the cube file (the angle is %.2f°).\nPlease check the cell vectors.\n", ang);
					abort();
				}
				ang = Angle_Deg(CxVector3(g_mCubeCell(1, 0), g_mCubeCell(1, 1), g_mCubeCell(1, 2)), vecb);
				if (ang > 1.0f) {
					eprintf("The direction of cell vector B does not match the grid data in the cube file (the angle is %.2f°).\nPlease check the cell vectors.\n", ang);
					abort();
				}
				ang = Angle_Deg(CxVector3(g_mCubeCell(2, 0), g_mCubeCell(2, 1), g_mCubeCell(2, 2)), vecc);
				if (ang > 1.0f) {
					eprintf("The direction of cell vector C does not match the grid data in the cube file (the angle is %.2f°).\nPlease check the cell vectors.\n", ang);
					abort();
				}
			}
			float diff;
			if (g_bBoxNonOrtho) {
				CxVector3 diffvec = CxVector3(g_fCubeXVector[0], g_fCubeXVector[1], g_fCubeXVector[2]) * (float)LEN_AU2PM * g_TimeStep.m_pVolumetricData->m_iRes[0];
				diff = diffvec.GetLength() - veca.GetLength();
			} else {
				diff = g_fBoxX - g_fCubeXStep * LEN_AU2PM * g_TimeStep.m_pVolumetricData->m_iRes[0];
			}
			if (fabsf(diff) > 0.1f) {
				mprintf("\n    The X cell vector does not match the grid data in the cube file (the difference is %.3f pm).\n", diff);
				mprintf("    This is possibly connected to the stride of the cube file.\n");
				g_iCubeXStride = AskInteger("    Which X stride was used to write the cube file? [2] ", 2);
				int i;
				for (i = 1; i < g_iCubeXStride; i++) {
					if (g_bBoxNonOrtho) {
						CxVector3 diffvec = CxVector3(g_fCubeXVector[0], g_fCubeXVector[1], g_fCubeXVector[2]) * (float)LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - i) / g_iCubeXStride;
						diff = CxVector3(g_fCubeXVector[0], g_fCubeXVector[1], g_fCubeXVector[2]).GetLength() - veca.GetLength();
					} else {
						diff = g_fBoxX - g_fCubeXStep * LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - i) / g_iCubeXStride;
					}
					if (fabsf(diff) < 0.1f) {
						mprintf("    Assuming that the original grid had %d points, the data would match (the difference is %.3f pm).\n", g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - i, diff);
						if (AskYesNo("    Take this value (y/n)? [yes] ", true)) {
							g_iCubeXMismatch = i;
							break;
						}
					}
				}
				if (g_iCubeXMismatch == 0) {
					mprintf("    The original number of grid points could not be determined.\n");
					g_iCubeXMismatch = g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - AskRangeInteger_ND("    Please enter the original number of grid points: ", (g_TimeStep.m_pVolumetricData->m_iRes[0] - 1) * g_iCubeXStride + 1, g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride);
				}
// 				if (g_bBoxNonOrtho) {
					g_fCubeXVector[0] /= g_iCubeXStride;
					g_fCubeXVector[1] /= g_iCubeXStride;
					g_fCubeXVector[2] /= g_iCubeXStride;
// 				} else {
					g_fCubeXStep /= g_iCubeXStride;
// 				}
			} else {
				mprintf("\n    X cell vector matches grid data in cube file (difference is %.3f pm).\n", diff);
			}
			
			if (g_bBoxNonOrtho) {
				CxVector3 diffvec = CxVector3(g_fCubeYVector[0], g_fCubeYVector[1], g_fCubeYVector[2]) * (float)LEN_AU2PM * g_TimeStep.m_pVolumetricData->m_iRes[1];
				diff = diffvec.GetLength() - vecb.GetLength();
			} else {
				diff = g_fBoxY - g_fCubeYStep * LEN_AU2PM * g_TimeStep.m_pVolumetricData->m_iRes[1];
			}
			if (fabsf(diff) > 0.1f) {
				mprintf("\n    The Y cell vector does not match the grid data in the cube file (the difference is %.3f pm).\n", diff);
				mprintf("    This is possibly connected to the stride of the cube file.\n");
				g_iCubeYStride = AskInteger("    Which Y stride was used to write the cube file? [2] ", 2);
				int i;
				for (i = 1; i < g_iCubeYStride; i++) {
					if (g_bBoxNonOrtho) {
						CxVector3 diffvec = CxVector3(g_fCubeYVector[0], g_fCubeYVector[1], g_fCubeYVector[2]) * (float)LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - i) / g_iCubeXStride;
						diff = CxVector3(g_fCubeYVector[0], g_fCubeYVector[1], g_fCubeYVector[2]).GetLength() - vecb.GetLength();
					} else {
						diff = g_fBoxY - g_fCubeYStep * LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - i) / g_iCubeYStride;
					}
					if (fabsf(diff) < 0.1f) {
						mprintf("    Assuming that the original grid had %d points, the data would match (the difference is %.3f pm).\n", g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - i, diff);
						if (AskYesNo("    Take this value (y/n)? [yes] ", true)) {
							g_iCubeYMismatch = i;
							break;
						}
					}
				}
				if (g_iCubeYMismatch == 0) {
					mprintf("    The original number of grid points could not be determined.\n");
					g_iCubeYMismatch = g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - AskRangeInteger_ND("    Please enter the original number of grid points: ", (g_TimeStep.m_pVolumetricData->m_iRes[1] - 1) * g_iCubeYStride + 1, g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride);
				}
// 				if (g_bBoxNonOrtho) {
					g_fCubeYVector[0] /= g_iCubeYStride;
					g_fCubeYVector[1] /= g_iCubeYStride;
					g_fCubeYVector[2] /= g_iCubeYStride;
// 				} else {
					g_fCubeYStep /= g_iCubeYStride;
// 				}
			} else {
				mprintf("\n    Y cell vector matches grid data in cube file (difference is %.3f pm).\n", diff);
			}
			
			if (g_bBoxNonOrtho) {
				CxVector3 diffvec = CxVector3(g_fCubeZVector[0], g_fCubeZVector[1], g_fCubeZVector[2]) * (float)LEN_AU2PM * g_TimeStep.m_pVolumetricData->m_iRes[2];
				diff = diffvec.GetLength() - vecc.GetLength();
			} else {
				diff = g_fBoxZ - g_fCubeZStep * LEN_AU2PM * g_TimeStep.m_pVolumetricData->m_iRes[2];
			}
			if (fabsf(diff) > 0.1f) {
				mprintf("\n    The Z cell vector does not match the grid data in the cube file (the difference is %.3f pm).\n", diff);
				mprintf("    This is possibly connected to the stride of the cube file.\n");
				g_iCubeZStride = AskInteger("    Which Z stride was used to write the cube file? [2] ", 2);
				int i;
				for (i = 1; i < g_iCubeZStride; i++) {
					if (g_bBoxNonOrtho) {
						CxVector3 diffvec = CxVector3(g_fCubeZVector[0], g_fCubeZVector[1], g_fCubeZVector[2]) * (float)LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - i) / g_iCubeZStride;
						diff = CxVector3(g_fCubeZVector[0], g_fCubeZVector[1], g_fCubeZVector[2]).GetLength() - vecc.GetLength();
					} else {
						diff = g_fBoxZ - g_fCubeZStep * LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - i) / g_iCubeZStride;
					}
					if (fabsf(diff) < 0.1f) {
						mprintf("    Assuming that the original grid had %d points, the data would match (the difference is %.3f pm).\n", g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - i, diff);
						if (AskYesNo("    Take this value (y/n)? [yes] ", true)) {
							g_iCubeZMismatch = i;
							break;
						}
					}
				}
				if (g_iCubeZMismatch == 0) {
					mprintf("    The original number of grid points could not be determined.\n");
					g_iCubeZMismatch = g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - AskRangeInteger_ND("    Please enter the original number of grid points: ", (g_TimeStep.m_pVolumetricData->m_iRes[2] - 1) * g_iCubeZStride + 1, g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride);
				}
// 				if (g_bBoxNonOrtho) {
					g_fCubeZVector[0] /= g_iCubeZStride;
					g_fCubeZVector[1] /= g_iCubeZStride;
					g_fCubeZVector[2] /= g_iCubeZStride;
// 				} else {
					g_fCubeZStep /= g_iCubeZStride;
// 				}
			} else {
				mprintf("\n    Z cell vector matches grid data in cube file (difference is %.3f pm).\n", diff);
			}
			
			if (g_bAdvanced1) {
				mprintf("\n    You may enter a cube interpolation factor along each axis (1 means no interpolation)\n");
				
				int factor = 0;
				while (true) {
					factor = AskInteger("    X axis: [1] ", 1);
					if (factor >= 1)
						break;
					mprintf(RED, "The factor has to be positive.\n");
				}
				if (factor > 1) {
					g_iCubeXStride *= factor;
					g_iCubeXMismatch *= factor;
					if (g_bBoxNonOrtho) {
						g_fCubeXVector[0] /= g_iCubeXStride;
						g_fCubeXVector[1] /= g_iCubeXStride;
						g_fCubeXVector[2] /= g_iCubeXStride;
					} else {
						g_fCubeXStep /= factor;
					}
				}
				
				while (true) {
					factor = AskInteger("    Y axis: [1] ", 1);
					if (factor >= 1)
						break;
					mprintf(RED, "The factor has to be positive.\n");
				}
				if (factor > 1) {
					g_iCubeYStride *= factor;
					g_iCubeYMismatch *= factor;
					if (g_bBoxNonOrtho) {
						g_fCubeYVector[0] /= g_iCubeYStride;
						g_fCubeYVector[1] /= g_iCubeYStride;
						g_fCubeYVector[2] /= g_iCubeYStride;
					} else {
						g_fCubeYStep /= factor;
					}
				}
				
				while (true) {
					factor = AskInteger("    Z axis: [1] ", 1);
					if (factor >= 1)
						break;
					mprintf(RED, "The factor has to be positive.\n");
				}
				if (factor > 1) {
					g_iCubeZStride *= factor;
					g_iCubeZMismatch *= factor;
					if (g_bBoxNonOrtho) {
						g_fCubeZVector[0] /= g_iCubeZStride;
						g_fCubeZVector[1] /= g_iCubeZStride;
						g_fCubeZVector[2] /= g_iCubeZStride;
					} else {
						g_fCubeZStep /= factor;
					}
				}
			}
			
			if (g_bBoxNonOrtho) {
				g_fCubeXVector[0] = g_mBoxFromOrtho(0, 0) / (g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM;
				g_fCubeXVector[1] = g_mBoxFromOrtho(0, 1) / (g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM;
				g_fCubeXVector[2] = g_mBoxFromOrtho(0, 2) / (g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM;
				g_fCubeYVector[0] = g_mBoxFromOrtho(1, 0) / (g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM;
				g_fCubeYVector[1] = g_mBoxFromOrtho(1, 1) / (g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM;
				g_fCubeYVector[2] = g_mBoxFromOrtho(1, 2) / (g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM;
				g_fCubeZVector[0] = g_mBoxFromOrtho(2, 0) / (g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM;
				g_fCubeZVector[1] = g_mBoxFromOrtho(2, 1) / (g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM;
				g_fCubeZVector[2] = g_mBoxFromOrtho(2, 2) / (g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM;
				int k, l;
				for (k = 0; k < 3; k++) {
					for (l = 0; l < 3; l++) {
						g_mCubeCell(k, l) = g_mBoxFromOrtho(k, l);
					}
				}
				g_TimeStep.m_pVolumetricData->m_fMaxVal[0] = sqrtf(g_mCubeCell(0, 0) * g_mCubeCell(0, 0) + g_mCubeCell(0, 1) * g_mCubeCell(0, 1) + g_mCubeCell(0, 2) * g_mCubeCell(0, 2));
				g_TimeStep.m_pVolumetricData->m_fMaxVal[1] = sqrtf(g_mCubeCell(1, 0) * g_mCubeCell(1, 0) + g_mCubeCell(1, 1) * g_mCubeCell(1, 1) + g_mCubeCell(1, 2) * g_mCubeCell(1, 2));
				g_TimeStep.m_pVolumetricData->m_fMaxVal[2] = sqrtf(g_mCubeCell(2, 0) * g_mCubeCell(2, 0) + g_mCubeCell(2, 1) * g_mCubeCell(2, 1) + g_mCubeCell(2, 2) * g_mCubeCell(2, 2));
			} else {
				if (fabsf(g_fCubeXStep - g_fBoxX / (g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM) > 0.1f || fabsf(g_fCubeYStep - g_fBoxY / (g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM) > 0.1f || fabsf(g_fCubeZStep - g_fBoxZ / (g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM) > 0.1f) {
					eprintf("Error while adapting cube grid step to cell size.\n");
					return false;
				}
				g_fCubeXStep = g_fBoxX / (g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM;
				g_fCubeYStep = g_fBoxY / (g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM;
				g_fCubeZStep = g_fBoxZ / (g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM;
				g_mCubeCell.Unity();
				g_mCubeCell(0, 0) = g_fBoxX;
				g_mCubeCell(1, 1) = g_fBoxY;
				g_mCubeCell(2, 2) = g_fBoxZ;
				g_TimeStep.m_pVolumetricData->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[0] * g_iCubeXStride - g_iCubeXMismatch);
				g_TimeStep.m_pVolumetricData->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[1] * g_iCubeYStride - g_iCubeYMismatch);
				g_TimeStep.m_pVolumetricData->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (g_TimeStep.m_pVolumetricData->m_iRes[2] * g_iCubeZStride - g_iCubeZMismatch);
			}
			
		}

_celldefined:

		if (g_bAdvanced1 && g_bPeriodic)
		{
			mprintf("\n");
			if (AskYesNo("    Should the periodic box be multiplied (y/n)? [no] ",false))
			{
				mprintf("\n");
				g_bDoubleBox = true;

				if (g_bPeriodicX)
					g_iDoubleBoxX = AskUnsignedInteger("    Replicate the box n times in X direction: [2] ",2);
						else g_iDoubleBoxX = 1;

				if (g_bPeriodicY)
					g_iDoubleBoxY = AskUnsignedInteger("    Replicate the box n times in Y direction: [2] ",2);
						else g_iDoubleBoxY = 1;

				if (g_bPeriodicZ)
					g_iDoubleBoxZ = AskUnsignedInteger("    Replicate the box n times in Z direction: [2] ",2);
						else g_iDoubleBoxZ = 1;

				g_iDoubleBoxFactor = g_iDoubleBoxX * g_iDoubleBoxY * g_iDoubleBoxZ;
				g_iGesAtomCount *= g_iDoubleBoxFactor;
				g_iGesVirtAtomCount *= g_iDoubleBoxFactor;

				if (g_bBoxNonOrtho)
				{
					g_mBoxFromOrtho(0,0) *= g_iDoubleBoxX;
					g_mBoxFromOrtho(1,0) *= g_iDoubleBoxX;
					g_mBoxFromOrtho(2,0) *= g_iDoubleBoxX;

					g_mBoxFromOrtho(0,1) *= g_iDoubleBoxY;
					g_mBoxFromOrtho(1,1) *= g_iDoubleBoxY;
					g_mBoxFromOrtho(2,1) *= g_iDoubleBoxY;

					g_mBoxFromOrtho(0,2) *= g_iDoubleBoxZ;
					g_mBoxFromOrtho(1,2) *= g_iDoubleBoxZ;
					g_mBoxFromOrtho(2,2) *= g_iDoubleBoxZ;

					g_mBoxToOrtho = CxMatrix3(g_mBoxFromOrtho);
					g_mBoxToOrtho.Invert();
				} else
				{
					g_fBoxX *= g_iDoubleBoxX;
					g_fBoxY *= g_iDoubleBoxY;
					g_fBoxZ *= g_iDoubleBoxZ;
				}
			} else
				g_bDoubleBox = false;
		} else
			g_bDoubleBox = false;
	} else
		g_bDoubleBox = false;

	try { g_pUniteTemp = new bool[g_iGesAtomCount]; } catch(...) { g_pUniteTemp = NULL; }
	if (g_pUniteTemp == NULL) NewException((double)g_iGesAtomCount*sizeof(bool),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (!g_bPeriodic)
		g_bNPT = false;

	if (g_bPeriodic)
	{
		g_fMinPeriodic = 1e30f;
		if (g_bPeriodicX)
			if (g_fMinPeriodic > g_fBoxX)
				g_fMinPeriodic = g_fBoxX;
		if (g_bPeriodicY)
			if (g_fMinPeriodic > g_fBoxY)
				g_fMinPeriodic = g_fBoxY;
		if (g_bPeriodicZ)
			if (g_fMinPeriodic > g_fBoxZ)
				g_fMinPeriodic = g_fBoxZ;
	}

//	SAVEPOS;
_molbegin:
	if (!g_bStreamInput)
	{
		if (g_bAdvanced1)
		{
			mprintf("\n");
			g_iScanMolStep = AskInteger("    Execute molecule recognition for which time step (-1 = disable)? [0] ",0);
		} else g_iScanMolStep = 0;

		a = fopen(g_sInputTraj,"rb");
		if ((g_bNPT) && (g_sNPTFile[0] != 0))
			g_fNPTFile = fopen(g_sNPTFile,"rt");
		if (g_iScanMolStep > 0)
		{
			mprintf("    Fast-forwarding to step %d...\n",g_iScanMolStep);
			mprintf(WHITE,"      [");
			for (z=0;z<g_iScanMolStep;z++)
			{
				if ((g_bNPT) && (g_sNPTFile[0] != 0))
//					fgets(buf,256,g_fNPTFile);
					buf.fgets(256,g_fNPTFile);
				if (fmod(z,g_iScanMolStep/60.0) < 1.0)
					mprintf(WHITE,"#");
				if (!g_TimeStep.SkipTimestep(a))
				{
					eprintf("Error.\n",g_iScanMolStep);
					goto _skipfail;
				}
			}
			mprintf(WHITE,"]\n");
			g_iFastForwardPos = ftell(a);
			mprintf("  Step %d begins at offset %lu (%.1f MB).\n\n",g_iScanMolStep+1,g_iFastForwardPos,g_iFastForwardPos/1024.0/1024.0);
		}
	_skipfail:;
		g_TimeStep.ReadTimestep(a,true);
		if (g_bNPT)
		{
			if (g_sNPTFile[0] != 0)
				g_TimeStep.ReadCellVector(g_fNPTFile);

			if (g_bBoxNonOrtho)
			{
				mprintf("\n    NPT: Using cell geometry ( %.2f pm | %.2f pm | %.2f pm | %.2f deg | %.2f deg | %.2f deg ) for molecule recognition.\n",sqrt(g_mBoxFromOrtho(0,0)*g_mBoxFromOrtho(0,0)+g_mBoxFromOrtho(0,1)*g_mBoxFromOrtho(0,1)+g_mBoxFromOrtho(0,2)*g_mBoxFromOrtho(0,2)),
					sqrt(g_mBoxFromOrtho(1,0)*g_mBoxFromOrtho(1,0)+g_mBoxFromOrtho(1,1)*g_mBoxFromOrtho(1,1)+g_mBoxFromOrtho(1,2)*g_mBoxFromOrtho(1,2)),
					sqrt(g_mBoxFromOrtho(2,0)*g_mBoxFromOrtho(2,0)+g_mBoxFromOrtho(2,1)*g_mBoxFromOrtho(2,1)+g_mBoxFromOrtho(2,2)*g_mBoxFromOrtho(2,2)),
					g_fBoxAngleA,g_fBoxAngleB,g_fBoxAngleC);
			} else
			{
				mprintf("\n    NPT: Using a cell vector of ( %.2f pm | %.2f pm | %.2f pm ) for molecule recognition.\n",g_fBoxX,g_fBoxY,g_fBoxZ);
			}
		}
	} else
	{
		g_iScanMolStep = 0;
		mprintf("\n    Stream input: Executing molecule recognition from first time step...\n");
	}

	g_iCloseAtomCounter = 0;

	if (!g_TimeStep.ScanMolecules())
		return false;

	if (!g_bStreamInput)
	{
		fclose(a);
		if (g_sNPTFile[0] != 0)
			fclose(g_fNPTFile);
	}

	mprintf(YELLOW,"\n*** The following %d kind%s of molecules ha%s been recognized:\n",g_oaMolecules.GetSize(),(g_oaMolecules.GetSize()>1)?"s":"",(g_oaMolecules.GetSize()>1)?"ve":"s");
	mprintf("    (ordered by molecular mass)\n\n");
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		mprintf(WHITE,"  - Molecule %d: %s ",z+1,((CMolecule*)g_oaMolecules[z])->m_sName);
		mprintf("(%d piece%s, %.2f g/mol)\n",((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize(),(((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize()>1)?"s":"",((CMolecule*)g_oaMolecules[z])->m_fMass);
		sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex[0]];
		mprintf("      (%d noneq. atom%s, %d noneq. bond%s, %d noneq. angle%s)\n",sm->m_iAtomClasses,(sm->m_iAtomClasses==1)?"":"s",sm->m_oaBondGroups.GetSize(),(sm->m_oaBondGroups.GetSize()==1)?"":"s",sm->m_oaAngleGroups.GetSize(),(sm->m_oaAngleGroups.GetSize()==1)?"":"s");
		if (((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms.GetSize() != 0)
		{
			if (((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms.GetSize() == 1)
				mprintf(WHITE,"      Detected %d ring:\n",((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms.GetSize());
			else
				mprintf(WHITE,"      Detected %d rings:\n",((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms.GetSize());

			for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms.GetSize();z2++)
			{
				mprintf(WHITE,"        %2d.) %2d-ring: ",z2+1,((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetSize());
				if (((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetSize() > 16)
				{
					for (z3=0;z3<5;z3++)
					{
						mprintf("%s%d",((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[z])->m_baAtomIndex[((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtomTypes[z2])->GetAt(z3)]])->m_sName,((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetAt(z3)+1);
						if (z3+1 < ((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetSize())
							mprintf(" - ");
					}
					mprintf(RED,"...");
					mprintf(" - ");
					for (z3=((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetSize()-5;z3<((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetSize();z3++)
					{
						mprintf("%s%d",((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[z])->m_baAtomIndex[((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtomTypes[z2])->GetAt(z3)]])->m_sName,((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetAt(z3)+1);
						if (z3+1 < ((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetSize())
							mprintf(" - ");
					}
				} else
				{
					for (z3=0;z3<((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetSize();z3++)
					{
						mprintf("%s%d",((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[z])->m_baAtomIndex[((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtomTypes[z2])->GetAt(z3)]])->m_sName,((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetAt(z3)+1);
						if (z3+1 < ((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[z2])->GetSize())
							mprintf(" - ");
					}
				}
				mprintf("\n");
				if (z2 >= 9)
				{
					mprintf(RED,"\n    Showing only the first 10 rings.\n");
					break;
				}
			}
		}

		mprintf("\n");
	}

	if (g_oaMolecules.GetSize() == 0)
	{
		eprintf("You don't have any atom/molecule in your simulation. This is probably not what you wanted.\n");
		eprintf("Stopping execution.\n");
		return false;
	}

/*	if (g_bVerbose)
	{
		mprintf(WHITE,">>> Output of the molecule tree >>>\n");
		g_TimeStep.PrintMegaTree();
		mprintf(WHITE,"<<< Output of the molecule tree <<<\n\n");
	}*/

_matrixagain:

	if (g_bAdvanced1)
	{
		g_bMegaMat = !AskYesNo("    Show bond matrices only for first representant of each molecule type (y/n)? [yes] ",true);
		g_bMatOnlyBind = AskYesNo("    Show only bonds in the bond matrices (y/n)? [yes] ",true);
	} else
	{
		g_bMegaMat = false;
		g_bMatOnlyBind = true;
	}
			
	mprintf(WHITE,"\nOutput of bond matrices:\n");
	g_TimeStep.PrintMatrix(!g_bMegaMat,g_bMatOnlyBind);

	mprintf("\n");
	if (!AskYesNo("    Accept these molecules (y) or change something (n)? [yes] ",true))
	{
		tb = false;
		mprintf("\n");
		mprintf("    If you want some atom types to form no bonds at all, assign covalent radius 0 to them.\n");
_modagain:
		mprintf(YELLOW,"\n    *** Modify Molecules ***\n\n");
		mprintf("    Your choices:\n\n");

		mprintf("    1.) Change covalent atom radii used for bond recognition\n");
		mprintf("    2.) Break specific bonds\n");
		mprintf("    3.) Rename elements\n");

		mprintf("\n");
		switch(AskRangeInteger("    Please select: [done] ",1,3,0))
		{
			case 1:
				mprintf("\n    These values have been used (covalent radii from literature multiplied with %.2f):\n\n",g_fBondFactor);
				for (z=0;z<g_oaAtoms.GetSize()-1;z++)
					mprintf("      - %-2s    %5.1f pm\n",((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius*g_fBondFactor);
				mprintf("\n");
_arnew:
				AskString("    Which radius do you want to change (RETURN=done)? ",&buf,"");
				if (strlen(buf) == 0)
					goto _ardone;
				for (z=0;z<g_oaAtoms.GetSize()-1;z++)
				{
					if (mystricmp(buf,((CAtom*)g_oaAtoms[z])->m_sName) == 0)
					{
						((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius = AskFloat("    Please enter new bond radius for %s in pm: [%.1f] ",((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius*g_fBondFactor,((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius*g_fBondFactor) / g_fBondFactor;
						if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius == 0)
							mprintf("\n      Atom %s will form no bonds now.\n\n",((CAtom*)g_oaAtoms[z])->m_sName);
						mprintf("\n");
						tb = true;
						goto _arnew;
					}
				}
				eprintf("    Atom \"%s\" is not in the system.\n\n",(const char*)buf);
				goto _arnew;
_ardone:
				mprintf("\n");
				break;

			case 2:
_breakagain:
				mprintf("\n");
				if (g_oaMolecules.GetSize() > 1)
				{
_breakmol:			ti = AskUnsignedInteger_ND("    Break bonds in which molecule (1-%d)? ",g_oaMolecules.GetSize())-1;
					if ((ti < 0) || (ti >= (int)g_oaMolecules.GetSize()))
					{
						eprintf("Wrong Input.\n\n");
						goto _breakmol;
					}
				} else
					ti = 0;

				mol = (CMolecule*)g_oaMolecules[ti];
				tb = true;
				if (AskYesNo("    Break all bonds within this molecule (y/n)? [no] ",false))
				{
					for (z=0;z<mol->m_laSingleMolIndex.GetSize();z++)
					{
						sm = (CSingleMolecule*)g_oaSingleMolecules[mol->m_laSingleMolIndex[z]];
						g_laBondBlackList.Append(&sm->m_laBonds);
					}
				} else
				{
_breaka1:			AskString_ND("    Enter 1st atom of the bond to break (e.g. O2): ",&buf);
					if (!ParseAtom(buf,ti,ty,rty,atom))
						goto _breaka1;
_breaka2:			AskString_ND("    Enter 2nd atom of the bond to break (e.g. O2): ",&buf);
					if (!ParseAtom(buf,ti,ty2,rty2,atom2))
						goto _breaka2;
					for (z=0;z<mol->m_laSingleMolIndex.GetSize();z++)
					{
						sm = (CSingleMolecule*)g_oaSingleMolecules[mol->m_laSingleMolIndex[z]];
						g_laBondBlackList.Add(((CxIntArray*)sm->m_oaAtomOffset[ty])->GetAt(atom));
						g_laBondBlackList.Add(((CxIntArray*)sm->m_oaAtomOffset[ty2])->GetAt(atom2));
					}
				}
				if (AskYesNo("\n    Break another bond (y/n)? [no] ",false))
					goto _breakagain;
				mprintf("\n");
				break;

			case 3:
/*********************************************************************************************************************/
_rennew:
				mprintf("\n    The system contains the following atoms:\n\n    ");
				tb2 = false;
				for (z=0;z<g_oaAtoms.GetSize();z++)
				{
					if (z == g_iVirtAtomType)
						continue;
					if (((CAtom*)g_oaAtoms[z])->m_pMergedTo != NULL)
						continue;
					if (tb2)
						mprintf(", ");
					tb2 = true;
					mprintf("%dx %s",((CAtom*)g_oaAtoms[z])->m_iCount,((CAtom*)g_oaAtoms[z])->m_sName);
				}
				mprintf("\n\n");

_renerr:
				AskString("    Which element to rename: [done] ",&buf,"");

				if (strlen(buf) == 0)
					goto _rendone;

				for (z=0;z<g_oaAtoms.GetSize();z++)
					if (mystricmp(buf,((CAtom*)g_oaAtoms[z])->m_sName) == 0)
						goto _renfound;
				eprintf("\n    Atom \"%s\" not in the system.\n\n",(const char*)buf);
				goto _renerr;
_renfound:
				tb = true;
				at = (CAtom*)g_oaAtoms[z];
				
				AskString_ND("    Enter new name for %s: ",&buf2,(const char*)buf);
				if (strlen(buf2) > 7)
				{
					eprintf("\n    Atom labels may have maximum length of 7 characters.\n\n");
					goto _renfound;
				}
				if (ContainsDigit(buf2))
				{
					eprintf("\n    Digits in element labels not allowed.\n\n");
					goto _renfound;
				}
				el = FindElement(buf2,true);
				if (el != NULL)
				{
					xAddAtom(buf2);
					if (AskYesNo("    Element %s is known. Merge %s atoms into %s (y/n)? [yes] ",true,(const char*)buf2,(const char*)buf,(const char*)buf2))
					{
						for (z2=0;z2<g_oaAtoms.GetSize();z2++)
						{
							if (((CAtom*)g_oaAtoms[z2])->m_pElement == el)
							{
								at = (CAtom*)g_oaAtoms[z2];
								at->m_iCount--;
								el = at->m_pElement;
								goto _ren1;
							}
						}
						eprintf("\n    Strange error ^^\n");
_ren1:;
					}
				} else
				{
					mprintf("\n      Adding element %s...\n\n",(const char*)buf2);

					try { at = new CAtom(); } catch(...) { at = NULL; }
					if (at == NULL) NewException((double)sizeof(CAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					at->m_iIndex = g_oaAtoms.GetSize();
					g_oaAtoms.Add(at);
					at->m_iCount = 0;

					try { at->m_pElement = new CElement(); } catch(...) { at->m_pElement = NULL; }
					if (at->m_pElement == NULL) NewException((double)sizeof(CElement),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					at->m_pElement->m_fRadius = AskFloat("    Please enter covalent radius in pm: [100.0] ",100.0f);
					at->m_pElement->m_fMass = AskFloat("    Please enter atom mass in u: [0] ",0);
					sprintf(at->m_sName,"%s",(const char*)buf2);
				}
				for (z2=0;z2<(int)g_TimeStep.m_iGesAtomCount;z2++)
				{
//					strcpy(buf,(char*)g_TimeStep.m_paLabels[z2]);
					buf.strcpy((char*)g_TimeStep.m_paLabels[z2]);
					ReplaceDigits(&buf);
					if (mystricmp(buf,((CAtom*)g_oaAtoms[z])->m_sName) == 0)
					{
						at->m_iCount++;
						((CAtom*)g_oaAtoms[z])->m_iCount--;
					}
				}
				((CAtom*)g_oaAtoms[z])->m_pMergedTo = at;
				mprintf("\n");
				goto _rennew;
_rendone:;
				mprintf("\n    The system contains the following atoms:\n\n    ");
				tb2 = false;
				for (z=0;z<g_oaAtoms.GetSize();z++)
				{
					if (z == g_iVirtAtomType)
						continue;
					if (((CAtom*)g_oaAtoms[z])->m_pMergedTo != NULL)
						continue;
					if (tb2)
						mprintf(", ");
					tb2 = true;
					mprintf("%dx %s",((CAtom*)g_oaAtoms[z])->m_iCount,((CAtom*)g_oaAtoms[z])->m_sName);
				}
				mprintf("\n");
/*********************************************************************************************************************/
				break;

			default:
				goto _fin;
		}
		goto _modagain;
_fin:
		if (tb)
		{
			mprintf(WHITE,"\n    Going back to molecule recognition with changed settings.\n\n");
			goto _molbegin;
		} else
			mprintf(WHITE,"\n    Nothing was changed, not repeating molecule recognition.\n\n");
	}

	if (g_bPeriodic)
	{
		tb = false;
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			if (((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms.GetSize() != 0)
			{
				if ((((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms.GetSize() > 5) || (((CxIntArray*)((CMolecule*)g_oaMolecules[z])->m_oaRingAtoms[0])->GetSize() > 12))
				{
					if (!tb)
					{
						tb = true;
						mprintf("\n    Some of your molecules are polycyclic. Probably some of them are periodic across the\n");
						mprintf("    cell boundary, i.e., are polymers of infinite extent (like, e.g., a periodic solid lattice).\n");
						mprintf("    Those need to be handled differently (wrapped atom-wise, etc.).\n\n");
					}
					((CMolecule*)g_oaMolecules[z])->m_bPolymer = AskYesNo("      Molecule %d (%s): Is this an infinite polymer/lattice (y/n)? [yes] ",true,z+1,((CMolecule*)g_oaMolecules[z])->m_sName);
				}
			}
		}

		mprintf("\n    Uniting molecules which have been broken by wrapping...\n");
		g_TimeStep.UniteMolecules(true);
	}
	
	mprintf("\n");
	mprintf("    You can create images of the structural formulas with the atom labels for easier\n");
	mprintf("    identification of the atoms (requires installed GraphViz package - see www.graphviz.org).\n\n");

	if (AskYesNo("    Create images of the structural formulas (y/n)? [no] ",false))
	{
		if (g_bAdvanced1)
			ti = AskUnsignedInteger("    How many iterations to perform for formula optimization? [10] ",10);
				else ti = 10;
		RenderStructFormulas(ti);
	}

	mprintf(WHITE,"\n    The atoms are currently ordered by topological priority.\n");

	if (g_bAdvanced1)
	{
		mprintf("\n");

		// This block written @ Ballmer peak ;-P Try it out yourself.
		if (AskYesNo("    Change the atom ordering in some molecule (y/n)? [no] ",false))
		{
			if (AskYesNo("    Order atoms like in the input file instead (y/n)? [no] ",false))
			{
				ReorderLikeInput(); // Crazy shit!! (Ballmer Peak)
			} else
			{
	_reordernext:
//				sprintf(buf,"    Change atom ordering in which of the molecules (");
				buf.sprintf("    Change atom ordering in which of the molecules (");
				for (z=0;z<g_oaMolecules.GetSize();z++)
				{
//					sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//					strcat(buf,buf2);
//					if (z < g_oaMolecules.GetSize()-1)
//						strcat(buf,", ");
					buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
					buf.strcat(buf2);
					if (z < g_oaMolecules.GetSize()-1)
						buf.strcat(", ");
				}
//				strcat(buf,")? ");
				buf.strcat(")? ");
				ti = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize())-1;
				mol = (CMolecule*)g_oaMolecules[ti];
				mol->m_oaNewNumbers.SetSize(mol->m_baAtomIndex.GetSize());
				for (z=0;z<mol->m_baAtomIndex.GetSize();z++)
				{
					mprintf("* Atom type %s\n",((CAtom*)g_oaAtoms[mol->m_baAtomIndex[z]])->m_sName);
					if (mol->m_oaNewNumbers[z] != NULL)
						delete mol->m_oaNewNumbers[z];

					try { mol->m_oaNewNumbers[z] = new CxIntArray("gather():mol->m_oaNewNumbers[z]"); } catch(...) { mol->m_oaNewNumbers[z] = NULL; }
					if (mol->m_oaNewNumbers[z] == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					for (z2=0;z2<mol->m_waAtomCount[z];z2++)
					{
		_renumberagain:
						ti2 = AskRangeInteger("    Which number should %s%d bear? [%d] ",1,mol->m_waAtomCount[z],z2+1,((CAtom*)g_oaAtoms[mol->m_baAtomIndex[z]])->m_sName,z2+1,z2+1);
						for (z3=0;z3<z2;z3++)
							if (((CxIntArray*)mol->m_oaNewNumbers[z])->GetAt(z3) == ti2-1)
							{
								eprintf("This number already was chosen for %s%d!\n",((CAtom*)g_oaAtoms[mol->m_baAtomIndex[z]])->m_sName,z3+1);
								goto _renumberagain;
							}
						(*((CxIntArray*)mol->m_oaNewNumbers[z])).Add(ti2-1);
					}
				}
				mprintf("\n    Reordering atoms...");
				ReorderAtoms(ti);
				mprintf("Done.\n\n");
				if (AskYesNo("    Change atom ordering of another molecule (y/n)? [no] ",false))
					goto _reordernext;
				goto _matrixagain;
			} // End "this block"

			mprintf(WHITE,"\n    Going back to molecule recognition with changed settings.\n\n");
			goto _molbegin;
		} // END IF REORDER
	} // END IF ADVANCED

	mprintf(WHITE,"\n    Defining virtual atom #1 as molecular Center of Geometry:\n");
	for (z0=0;z0<g_oaMolecules.GetSize();z0++)
	{
		mprintf("      - %s...\n",((CMolecule*)g_oaMolecules[z0])->m_sName);
		va = AddVirtualAtom(z0);
		va->m_iMode = 0;
		va->m_oCenterAtoms.AddAllAtoms((CMolecule*)g_oaMolecules[va->m_iMolecule],false);
		va->m_faWeight.SetSize(va->m_oCenterAtoms.m_iAtomGes);
		for (z=0;z<va->m_oCenterAtoms.m_iAtomGes;z++)
			va->m_faWeight[z] = 1.0f;
		va->m_fGesWeight = (float)va->m_oCenterAtoms.m_iAtomGes;
	}

	mprintf(WHITE,"\n    Defining virtual atom #2 as molecular Center of Mass:\n");
	for (z0=0;z0<g_oaMolecules.GetSize();z0++)
	{
		mprintf("      - %s...\n",((CMolecule*)g_oaMolecules[z0])->m_sName);
		va = AddVirtualAtom(z0);
		va->m_iMode = 0;
		va->m_oCenterAtoms.AddAllAtoms((CMolecule*)g_oaMolecules[va->m_iMolecule],false);
		va->m_faWeight.SetSize(va->m_oCenterAtoms.m_iAtomGes);
		z3 = 0;
		tf2 = 0;
		for (z=0;z<va->m_oCenterAtoms.m_baAtomType.GetSize();z++)
		{
			tf = ((CAtom*)g_oaAtoms[va->m_oCenterAtoms.m_baRealAtomType[z]])->m_pElement->m_fMass;
			for (z2=0;z2<((CxIntArray*)va->m_oCenterAtoms.m_oaAtoms[z])->GetSize();z2++)
			{
				va->m_faWeight[z3] = tf;
				tf2 += tf;
				z3++;
			}
		}
		va->m_fGesWeight = tf2;
		if (va->m_fGesWeight == 0)
		{
			eprintf("        Molecule %s has total mass of zero. Defining #2 as Center of Geometry.\n",((CMolecule*)g_oaMolecules[z0])->m_sName);
			for (z=0;z<va->m_oCenterAtoms.m_iAtomGes;z++)
				va->m_faWeight[z] = 1.0f;
			va->m_fGesWeight = (float)va->m_oCenterAtoms.m_iAtomGes;
		}
	}

	for (z0=0;z0<g_oaMolecules.GetSize();z0++)
	{
		m = (CMolecule*)g_oaMolecules[z0];
		if (m->m_oaRingAtoms.GetSize() != 0)
		{
			if (m->m_bPolymer)
			{
				mprintf(WHITE,"\n    Skipping ring centers of %s: Is an infinite polymer/lattice.\n",m->m_sName);
			} else
			{
				mprintf(WHITE,"\n    Defining ring centers in %s:\n",m->m_sName);
				for (z=0;z<m->m_oaRingAtoms.GetSize();z++)
				{
					mprintf("      - Defining #%d as Center of Ring %d: ",z+3,z+1);
					for (z3=0;z3<((CxIntArray*)m->m_oaRingAtoms[z])->GetSize();z3++)
					{
						mprintf("%s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[((CxIntArray*)m->m_oaRingAtomTypes[z])->GetAt(z3)]])->m_sName,((CxIntArray*)m->m_oaRingAtoms[z])->GetAt(z3)+1);
						if (z3+1 < ((CxIntArray*)m->m_oaRingAtoms[z])->GetSize())
							mprintf(" - ");
					}
					mprintf("\n");
					va = AddVirtualAtom(z0);
					va->m_iMode = 0;
					va->m_oCenterAtoms.m_pMolecule = m;
					for (z3=0;z3<((CxIntArray*)m->m_oaRingAtoms[z])->GetSize();z3++)
						va->m_oCenterAtoms.AddAtom(((CxIntArray*)m->m_oaRingAtomTypes[z])->GetAt(z3),((CxIntArray*)m->m_oaRingAtoms[z])->GetAt(z3),false);
					va->m_faWeight.SetSize(va->m_oCenterAtoms.m_iAtomGes);
					for (z3=0;z3<va->m_oCenterAtoms.m_iAtomGes;z3++)
						va->m_faWeight[z3] = 1.0f;
					va->m_fGesWeight = (float)va->m_oCenterAtoms.m_iAtomGes;
				}
			}
		}
	}

	if (g_bAdvanced1)
	{
		mprintf("\n");
		if (AskYesNo("\n    Define additional virtual atoms (y/n)? [no] ",false))
		{
	_vabeg:
//			sprintf(buf,"    To which molecule shall the virtual atom belong (");
			buf.sprintf("    To which molecule shall the virtual atom belong (");
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
//				sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//				strcat(buf,buf2);
//				if (z < g_oaMolecules.GetSize()-1)
//					strcat(buf,", ");
				buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
				buf.strcat(buf2);
				if (z < g_oaMolecules.GetSize()-1)
					buf.strcat(", ");
			}
//			strcat(buf,")? ");
			buf.strcat(")? ");

			z = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize());
			va = AddVirtualAtom(z-1);
			mprintf(WHITE,"\n*** Defining virtual atom #%d in %s\n",va->m_iMolVirtAtom+1,((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName);
			va->m_iMode = AskRangeInteger("    Define v.a. as center (0) or through distance/angle/dihedral (1)? [center] ",0,1,0);
			if (va->m_iMode == 0)
			{
	_vaatoms:
				mprintf("    Which atoms to use for center (e.g. \"C1,C3-5,H\")? [all] ");
				inpprintf("! Which atoms to use for center (e.g. \"C1,C3-5,H\")? [all]\n");
				myget(&buf);
				if (strlen(buf)==0)
					va->m_oCenterAtoms.AddAllAtoms((CMolecule*)g_oaMolecules[va->m_iMolecule],false);
				else if (!va->m_oCenterAtoms.ParseAtoms((CMolecule*)g_oaMolecules[va->m_iMolecule],buf))
					goto _vaatoms;
				va->m_faWeight.SetSize(va->m_oCenterAtoms.m_iAtomGes);
	_vaweight:
				mprintf("    How shall these atoms be weightened (RETURN=equal, *=mass, #=manually)? [equal] ");
				inpprintf("! How shall these atoms be weightened (RETURN=equal, *=mass, #=manually)? [equal]\n");
				myget(&buf);
				if (strlen(buf)==0)
				{
					for (z=0;z<va->m_oCenterAtoms.m_iAtomGes;z++)
						va->m_faWeight[z] = 1.0f;
					va->m_fGesWeight = (float)va->m_oCenterAtoms.m_iAtomGes;
				} else if (buf[0] == '*')
				{
					z3 = 0;
					tf2 = 0;
					for (z=0;z<va->m_oCenterAtoms.m_baAtomType.GetSize();z++)
					{
						tf = ((CAtom*)g_oaAtoms[va->m_oCenterAtoms.m_baRealAtomType[z]])->m_pElement->m_fMass;
						for (z2=0;z2<((CxIntArray*)va->m_oCenterAtoms.m_oaAtoms[z])->GetSize();z2++)
						{
							va->m_faWeight[z3] = tf;
							tf2 += tf;
							z3++;
						}
					}
					va->m_fGesWeight = tf2;
					if (va->m_fGesWeight == 0)
					{
						eprintf("\n    Molecule has total mass of zero. Weighting atoms equally.\n\n");
						for (z=0;z<va->m_oCenterAtoms.m_iAtomGes;z++)
							va->m_faWeight[z] = 1.0f;
						va->m_fGesWeight = (float)va->m_oCenterAtoms.m_iAtomGes;
					}
				} else if (buf[0] == '#')
				{
_vawagain:
					z3 = 0;
					tf2 = 0;
					mprintf("\n");
					for (z=0;z<va->m_oCenterAtoms.m_baAtomType.GetSize();z++)
					{
						for (z2=0;z2<((CxIntArray*)va->m_oCenterAtoms.m_oaAtoms[z])->GetSize();z2++)
						{
							tf = AskFloat("      Enter the weight for %s%d: [1] ",1.0f,((CAtom*)g_oaAtoms[va->m_oCenterAtoms.m_baRealAtomType[z]])->m_pElement->m_sLabel, ((CxIntArray *)va->m_oCenterAtoms.m_oaAtoms[z])->GetAt(z2)+1);
							va->m_faWeight[z3] = tf;
							tf2 += tf;
							z3++;
						}
					}
					mprintf("    The sum of weights is %.4f.\n\n",tf2);
					va->m_fGesWeight = tf2;
					if (va->m_fGesWeight == 0)
					{
						eprintf("Sum of weights may not be zero. Enter the weights again.\n\n");
						goto _vawagain;
					}
				} else
				{
					eprintf("Wrong input.\n");
					inpprintf("! Wrong input.\n");
					goto _vaweight;
				}
			} else
			{
	_vabond:
				mprintf("    Enter 2nd atom for distance #%d- (e.g. C2): ",va->m_iMolVirtAtom+1);
				inpprintf("! Enter 2nd atom for distance #%d- (e.g. C2):\n",va->m_iMolVirtAtom+1);
				myget(&buf);
				if (!ParseAtom(buf,va->m_iMolecule,va->m_iAtomType[0],va->m_iRealAtomType[0],va->m_iAtom[0]))
					goto _vabond;
				if ((va->m_iRealAtomType[0] == g_iVirtAtomType) && (va->m_iAtom[0] == va->m_iMolVirtAtom))
				{
					eprintf("This atom was already chosen.\n");
					goto _vabond;
				}
	_vaangle:
				mprintf("    Enter 3rd atom for angle #%d-%s%d- (e.g. C2)? ",va->m_iMolVirtAtom+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1);
				inpprintf("! Enter 3rd atom for angle #%d-%s%d- (e.g. C2)?\n",va->m_iMolVirtAtom+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1);
				myget(&buf);
				if (!ParseAtom(buf,va->m_iMolecule,va->m_iAtomType[1],va->m_iRealAtomType[1],va->m_iAtom[1]))
					goto _vaangle;
				if (((va->m_iRealAtomType[1] == g_iVirtAtomType) && (va->m_iAtom[1] == va->m_iMolVirtAtom)) || ((va->m_iRealAtomType[1] == va->m_iRealAtomType[0]) && (va->m_iAtom[1] == va->m_iAtom[0])))
				{
					eprintf("This atom was already chosen.\n");
					goto _vaangle;
				}
	_vadihedral:
				mprintf("    Enter 4th atom for dihedral #%d-%s%d-%s%d- (e.g. C2)? ",va->m_iMolVirtAtom+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[1]])->m_sName,va->m_iAtom[1]+1);
				inpprintf("! Enter 4th atom for dihedral #%d-%s%d-%s%d- (e.g. C2)?\n",va->m_iMolVirtAtom+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[1]])->m_sName,va->m_iAtom[1]+1);
				myget(&buf);
				if (!ParseAtom(buf,va->m_iMolecule,va->m_iAtomType[2],va->m_iRealAtomType[2],va->m_iAtom[2]))
					goto _vadihedral;
				if (((va->m_iRealAtomType[2] == g_iVirtAtomType) && (va->m_iAtom[2] == va->m_iMolVirtAtom)) || ((va->m_iRealAtomType[2] == va->m_iRealAtomType[0]) && (va->m_iAtom[2] == va->m_iAtom[0])) || ((va->m_iRealAtomType[2] == va->m_iRealAtomType[1]) && (va->m_iAtom[2] == va->m_iAtom[1])))
				{
					eprintf("This atom was already chosen.\n");
					goto _vadihedral;
				}
		
				va->m_fValues[0] = AskFloat_ND("    Enter distance #%d-%s%d in pm: ",va->m_iMolVirtAtom+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1);
				va->m_fValues[1] = AskFloat_ND("    Enter angle #%d-%s%d-%s%d in degree: ",va->m_iMolVirtAtom+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[1]])->m_sName,va->m_iAtom[1]+1) * (float)Pi / 180.0f;
				va->m_fValues[2] = AskFloat_ND("    Enter dihedral #%d-%s%d-%s%d-%s%d in degree: ",va->m_iMolVirtAtom+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[1]])->m_sName,va->m_iAtom[1]+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[2]])->m_sName,va->m_iAtom[2]+1) * (float)Pi / 180.0f;
			}

			if (AskYesNo("\n    Define further virtual atoms (y/n)? [no] ",false))
				goto _vabeg;
		}

		mprintf("\n");
		if (AskYesNo("    Do you want to define pseudo-molecules (y/n)? [no] ",false))
		{
			tb2 = false;
			ti = 0;
			mprintf("\n    If you need the center of mass or center of geometry of the whole system, you should\n");
			mprintf("    define a pseudo-molecule over the whole system, and use its atoms #1 and #2.\n\n");
			if (AskYesNo("    Define a pseudo-molecule over the whole system (y/n)? [yes] ",true))
			{
				ti++;
/************************************************************************/
				mprintf("\n");
				try { m2 = new CMolecule(); } catch(...) { m2 = NULL; }
				if (m2 == NULL) NewException((double)sizeof(CMolecule),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				try { sm2 = new CSingleMolecule(); } catch(...) { sm2 = NULL; }
				if (sm2 == NULL) NewException((double)sizeof(CSingleMolecule),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				m2->m_bPseudo = true;
				sm2->m_bPseudo = true;
				m2->m_iIndex = g_oaMolecules.GetSize();
				sm2->m_iMolType = g_oaMolecules.GetSize();
				sm2->m_iMolSMIndex = 0;
				m2->m_laSingleMolIndex.Add(g_oaSingleMolecules.GetSize());
				g_oaMolecules.Add(m2);
				g_oaSingleMolecules.Add(sm2);

				mprintf(WHITE,"    >>> Pseudo-molecule %d (= molecule %d) >>>\n\n",ti,g_oaMolecules.GetSize());
				mprintf("    Adding all atoms from all molecules...\n");

				for (z=0;z<g_oaMolecules.GetSize();z++)
				{
					m = (CMolecule*)g_oaMolecules[z];
					if (m->m_bPseudo)
						continue;

					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					ag->AddAllAtoms(m,false);

					ia.RemoveAll();

					for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
						ia.Add(z2);

					for (z2=0;z2<ag->m_baAtomType.GetSize();z2++)
					{
						m2->m_iAtomGes += ((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize();

						for (z4=0;z4<m2->m_baAtomIndex.GetSize();z4++)
						{
							if (m2->m_baAtomIndex[z4] == ag->m_baRealAtomType[z2])
							{
								m2->m_waAtomCount[z4] += ((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize();
								pia = (CxIntArray*)sm2->m_oaAtomOffset[z4];
								goto _pmafound;
							}
						}

						m2->m_baAtomIndex.Add(ag->m_baRealAtomType[z2]);
						sm2->m_baAtomIndex.Add(ag->m_baRealAtomType[z2]);
						m2->m_waAtomCount.Add(((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize());

						try { pia = new CxIntArray("gather():pia"); } catch(...) { pia = NULL; }
						if (pia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

						sm2->m_oaAtomOffset.Add(pia);
_pmafound:
						for (z4=0;z4<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z4++)
						{
							for (z5=0;z5<ia.GetSize();z5++)
							{
								sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[ia[z5]]];
								for (z6=0;z6<sm2->m_oaAtomOffset.GetSize();z6++)
								{
									for (z7=0;z7<((CxIntArray*)sm2->m_oaAtomOffset[z6])->GetSize();z7++)
									{
										if (((CxIntArray*)sm2->m_oaAtomOffset[z6])->GetAt(z7) == ((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)))
										{
											eprintf("Weird error: Atom %s%d from %s[%d] is already in the pseudo-molecule; skipping.\n",((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z2]])->m_sName,((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)+1,m->m_sName,ia[z5]+1);
											m2->m_waAtomCount[m2->m_waAtomCount.GetSize()-1]--;
											abort();
										}
									}
								}
								pia->Add(((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)));
							}
						}
					}
					delete ag;
				}
				mprintf("\n");
/********************************************************************************/
				goto _nopseudomol;
			}
	_pseudomolnext:
			ti++;

			if (!tb2)
			{
				tb2 = true;
				mprintf("\n    You may define custom element labels. This enables keeping apart\n    atoms of the same type from different molecules/positions.\n\n");
				tb = AskYesNo("    Keep the standard element labels (y) or define custom labels (n)? [yes] ",true);
				mprintf("\n");
			}

			try { m2 = new CMolecule(); } catch(...) { m2 = NULL; }
			if (m2 == NULL) NewException((double)sizeof(CMolecule),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			try { sm2 = new CSingleMolecule(); } catch(...) { sm2 = NULL; }
			if (sm2 == NULL) NewException((double)sizeof(CSingleMolecule),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			m2->m_bPseudo = true;
			sm2->m_bPseudo = true;
			m2->m_iIndex = g_oaMolecules.GetSize();
			sm2->m_iMolType = g_oaMolecules.GetSize();
			sm2->m_iMolSMIndex = 0;
			m2->m_laSingleMolIndex.Add(g_oaSingleMolecules.GetSize());
			g_oaMolecules.Add(m2);
			g_oaSingleMolecules.Add(sm2);

			mprintf(WHITE,"    >>> Pseudo-molecule %d (= molecule %d) >>>\n\n",ti,g_oaMolecules.GetSize());
	/*		if (ti == 1)
				mprintf(WHITE,"      You may define custom element labels. This enables keeping apart\n      atoms of the same type from different molecules/positions.\n\n");
			tb = AskYesNo("      Keep the standard element labels (y) or define custom labels (n)? [yes] ",true);
			mprintf("\n");
*/
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				m = (CMolecule*)g_oaMolecules[z];
				if (m->m_bPseudo)
					continue;
				if (AskYesNo("    Use atoms from molecule %d (%s) (y/n)? [no] ",false,z+1,m->m_sName))
				{
					mprintf("\n");
					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (AskYesNo("      Use all atoms (y) from %s or only certain atoms (n)? [no] ",false,m->m_sName))
					{
						ag->AddAllAtoms(m,false);
					} else
					{
	_pseudomolp1:
						AskString_ND("      Please enter atoms from %s to use (e.g. C1-3,H,O4): ",&buf,m->m_sName);
						if (!ag->ParseAtoms(m,buf))
							goto _pseudomolp1;
					}
					ia.RemoveAll();
					if (m->m_laSingleMolIndex.GetSize() > 1)
					{
						if (!AskYesNo("      Add selected atoms from all %s molecules (y), or only from certain molecules (n)? [yes] ",true,m->m_sName))
						{
	_pseudomolil:
							AskString_ND("      Enter the %s molecules to use (range 1-%d, e.g. 1,3-6,8): ",&buf,m->m_sName,m->m_laSingleMolIndex.GetSize());
							ia.RemoveAll();
							if (!ParseIntList(buf,&ia))
								goto _pseudomolil;
							for (z2=0;z2<ia.GetSize();z2++)
							{
								if ((ia[z2] < 1) || (ia[z2] > m->m_laSingleMolIndex.GetSize()))
								{
									eprintf("Invalid number: %d (should be between 1 and %d).\n",ia[z2],m->m_laSingleMolIndex.GetSize());
									goto _pseudomolil;
								}
								ia[z2]--;
							}
						} else
						{
							for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
								ia.Add(z2);
						}
					} else ia.Add(0);

					if (!tb) // Custom Labels
					{
						for (z2=0;z2<ag->m_baAtomType.GetSize();z2++)
						{
	_pmlagain:
							AskString("      How should %s from %s be labeled? [%s] ",&buf,((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z2]])->m_sName,((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z2]])->m_sName,m->m_sName,((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z2]])->m_sName);
							if (ContainsDigit(buf))
							{
								eprintf("Digits in element labels not allowed.\n");
								goto _pmlagain;
							}
							for (z3=0;z3<g_oaAtoms.GetSize();z3++)
							{
								if (mystricmp(buf,((CAtom*)g_oaAtoms[z3])->m_sName) == 0)
								{
									ti2 = z3;
									goto _pseudomold1;
								}
							}
							if (strlen(buf) > 7)
							{
								mprintf("        Element labels may only have up to 7 characters; truncating.\n");
								buf(7) = 0;
							}
							mprintf("        Adding new element label \"%s\"...\n",(const char*)buf);

							try { at = new CAtom(); } catch(...) { at = NULL; }
							if (at == NULL) NewException((double)sizeof(CAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
							
							at->m_iIndex = g_oaAtoms.GetSize();
							at->m_pElement = ((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z2]])->m_pElement;
							memcpy(at->m_sName,buf,8);
							at->m_iCount = ((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize();
							ti2 = g_oaAtoms.GetSize();
							g_oaAtoms.Add(at);
	_pseudomold1:
							m2->m_iAtomGes += ((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize();
							for (z3=0;z3<m2->m_baAtomIndex.GetSize();z3++)
							{
								if (m2->m_baAtomIndex[z3] == ti2)
								{
									m2->m_waAtomCount[z3] += ((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize();
									for (z4=0;z4<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z4++)
									{
										for (z5=0;z5<ia.GetSize();z5++)
										{
											sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[ia[z5]]];
											for (z6=0;z6<sm2->m_oaAtomOffset.GetSize();z6++)
											{
												for (z7=0;z7<((CxIntArray*)sm2->m_oaAtomOffset[z6])->GetSize();z7++)
												{
													if (((CxIntArray*)sm2->m_oaAtomOffset[z6])->GetAt(z7) == ((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)))
													{
														eprintf("Atom %s%d from %s[%d] is already in the pseudo-molecule; skipping.\n",((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z2]])->m_sName,((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)+1,m->m_sName,ia[z5]+1);
														m2->m_waAtomCount[z3]--;
														goto _pseudomolsk1;
													}
												}
											}
							//				mprintf("Molecule %d, SM %d, Type %d, Atom %d: Line %d.\n",z,ia[z5],((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)));
											((CxIntArray*)sm2->m_oaAtomOffset[z3])->Add(((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)));
	_pseudomolsk1:;
										}
									}
									goto _pseudomold2;
								}
							}
							m2->m_baAtomIndex.Add(ti2);
							sm2->m_baAtomIndex.Add(ti2);
							m2->m_waAtomCount.Add(((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize());

							try { pia = new CxIntArray("gather():pia"); } catch(...) { pia = NULL; }
							if (pia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
							
							sm2->m_oaAtomOffset.Add(pia);
							for (z4=0;z4<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z4++)
							{
								for (z5=0;z5<ia.GetSize();z5++)
								{
									sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[ia[z5]]];
									for (z6=0;z6<sm2->m_oaAtomOffset.GetSize();z6++)
									{
										for (z7=0;z7<((CxIntArray*)sm2->m_oaAtomOffset[z6])->GetSize();z7++)
										{
											if (((CxIntArray*)sm2->m_oaAtomOffset[z6])->GetAt(z7) == ((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)))
											{
												eprintf("Atom %s%d from %s[%d] is already in the pseudo-molecule; skipping.\n",((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z2]])->m_sName,((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)+1,m->m_sName,ia[z5]+1);
												m2->m_waAtomCount[m2->m_waAtomCount.GetSize()-1]--;
												goto _pseudomolsk2;
											}
										}
									}
									pia->Add(((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)));
	_pseudomolsk2:;
								}
							}
	_pseudomold2:;
						}
					} else // If not custom labels
					{
						for (z2=0;z2<ag->m_baAtomType.GetSize();z2++)
						{
							for (z4=0;z4<m2->m_baAtomIndex.GetSize();z4++)
							{
								if (m2->m_baAtomIndex[z4] == ag->m_baRealAtomType[z2])
								{
									m2->m_waAtomCount[z4] += ((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize();
									pia = (CxIntArray*)sm2->m_oaAtomOffset[z4];
									goto _pmafound2;
								}
							}

							m2->m_iAtomGes += ((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize();
							m2->m_baAtomIndex.Add(ag->m_baRealAtomType[z2]);
							sm2->m_baAtomIndex.Add(ag->m_baRealAtomType[z2]);
							m2->m_waAtomCount.Add(((CxIntArray*)ag->m_oaAtoms[z2])->GetSize() * ia.GetSize());

							try { pia = new CxIntArray("gather():pia"); } catch(...) { pia = NULL; }
							if (pia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

							sm2->m_oaAtomOffset.Add(pia);
_pmafound2:
							for (z4=0;z4<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z4++)
							{
								for (z5=0;z5<ia.GetSize();z5++)
								{
									sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[ia[z5]]];
									for (z6=0;z6<sm2->m_oaAtomOffset.GetSize();z6++)
									{
										for (z7=0;z7<((CxIntArray*)sm2->m_oaAtomOffset[z6])->GetSize();z7++)
										{
											if (((CxIntArray*)sm2->m_oaAtomOffset[z6])->GetAt(z7) == ((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)))
											{
												eprintf("Atom %s%d from %s[%d] is already in the pseudo-molecule; skipping.\n",((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z2]])->m_sName,((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)+1,m->m_sName,ia[z5]+1);
												m2->m_waAtomCount[m2->m_waAtomCount.GetSize()-1]--;
												goto _pseudomolsk3;
											}
										}
									}
									pia->Add(((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z4)));
	_pseudomolsk3:;
								}
							}
						}
					}
					delete ag;
					if (AskYesNo("      Use another set of atoms from %s (y/n)? [no] ",false,m->m_sName))
					{
						try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
						if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						goto _pseudomolp1;
					}
				}
				mprintf("\n");
			}
_nopseudomol:
			m2->BuildName();
			mprintf(WHITE,"    <<< Pseudo-molecule %d (= molecule %d) defined as %s <<<\n\n",ti,g_oaMolecules.GetSize(),m2->m_sName);
			if (AskYesNo("    Define another pseudo-molecule (y/n)? [no] ",false))
				goto _pseudomolnext;
			mprintf(WHITE,"\n    The following pseudo-molecules have been defined:\n");
			mprintf("      (note: pseudo-molecules start with a $ sign)\n\n");
			ti = 0;
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				m = (CMolecule*)g_oaMolecules[z];
				if (m->m_bPseudo)
				{
					ti++;
					mprintf("      %d.) Molecule %d - %s\n",ti,z+1,m->m_sName);
				}
			}
			mprintf(WHITE,"\n    Defining virtual atom #1 as molecular Center of Geometry:\n");
			for (z0=0;z0<g_oaMolecules.GetSize();z0++)
			{
				m = (CMolecule*)g_oaMolecules[z0];
				if (!m->m_bPseudo)
					continue;
				mprintf("      - %s...\n",m->m_sName);
				va = AddVirtualAtom(z0);
				va->m_iMode = 0;
				va->m_oCenterAtoms.AddAllAtoms((CMolecule*)g_oaMolecules[va->m_iMolecule],false);
				va->m_faWeight.SetSize(va->m_oCenterAtoms.m_iAtomGes);
				for (z=0;z<va->m_oCenterAtoms.m_iAtomGes;z++)
					va->m_faWeight[z] = 1.0f;
				va->m_fGesWeight = (float)va->m_oCenterAtoms.m_iAtomGes;
			}

			mprintf(WHITE,"\n    Defining virtual atom #2 as molecular Center of Mass:\n");
			for (z0=0;z0<g_oaMolecules.GetSize();z0++)
			{
				m = (CMolecule*)g_oaMolecules[z0];
				if (!m->m_bPseudo)
					continue;
				mprintf("      - %s...\n",m->m_sName);
				va = AddVirtualAtom(z0);
				va->m_iMode = 0;
				va->m_oCenterAtoms.AddAllAtoms((CMolecule*)g_oaMolecules[va->m_iMolecule],false);
				va->m_faWeight.SetSize(va->m_oCenterAtoms.m_iAtomGes);
				z3 = 0;
				tf2 = 0;
				for (z=0;z<va->m_oCenterAtoms.m_baAtomType.GetSize();z++)
				{
					tf = ((CAtom*)g_oaAtoms[va->m_oCenterAtoms.m_baRealAtomType[z]])->m_pElement->m_fMass;
					for (z2=0;z2<((CxIntArray*)va->m_oCenterAtoms.m_oaAtoms[z])->GetSize();z2++)
					{
						va->m_faWeight[z3] = tf;
						tf2 += tf;
						z3++;
					}
				}
				va->m_fGesWeight = tf2;
				if (va->m_fGesWeight == 0)
				{
					eprintf("        Molecule %s has total mass of zero. Defining #2 as Center of Geometry.\n",m->m_sName);
					for (z=0;z<va->m_oCenterAtoms.m_iAtomGes;z++)
						va->m_faWeight[z] = 1.0f;
					va->m_fGesWeight = (float)va->m_oCenterAtoms.m_iAtomGes;
				}
			}
		}
	}

	mprintf(WHITE,"\n    Defining molecular dipole vectors:\n");
	mprintf("      The last virtual atom of each molecule is defined as the tip of the dipole vector\n");
	mprintf("      starting in center of mass (#2). 1 Debye corresponds to 100pm of vector length.\n");
	mprintf("      Only useful if dipole/charge data is provided. Otherwise vector has always length of zero.\n\n");
	for (z0=0;z0<g_oaMolecules.GetSize();z0++)
	{
		m = (CMolecule*)g_oaMolecules[z0];
		mprintf("      - Dipole vector #2 --> #%d in %s ...\n",m->m_laVirtualAtoms.GetSize()+1,m->m_sName);
		va = AddVirtualAtom(z0);
		va->m_iMode = 2;
	}

	mprintf(GREEN,"\n>>> %d virtual atoms have been defined: >>>\n",g_oaVirtualAtoms.GetSize());
	for (z0=0;z0<g_oaMolecules.GetSize();z0++)
	{
		for (z=0;z<((CMolecule*)g_oaMolecules[z0])->m_laVirtualAtoms.GetSize();z++)
		{
			va = (CVirtualAtom*)g_oaVirtualAtoms[((CMolecule*)g_oaMolecules[z0])->m_laVirtualAtoms[z]];
			if (va->m_iMolVirtAtom < 2)
			{
				mprintf(WHITE,"\n    #%d in %s: %s.",va->m_iMolVirtAtom+1,((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,(va->m_iMolVirtAtom==0)?"Center of Geometry":"Center of Mass");
				if (va->m_iMolVirtAtom == 1)
					mprintf("\n");
				continue;
			}
			if (va->m_iMode < 2)
			{
				mprintf(WHITE,"\n    #%d in %s - defined through %s\n",va->m_iMolVirtAtom+1,((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName,(va->m_iMode==0)?"center:":"distance/angle/dihedral:");
				if (va->m_iMode == 0)
				{
					if (!((CMolecule*)g_oaMolecules[z0])->m_bPseudo || (z > 1))
					{
						z4 = 0;
						for (z2=0;z2<va->m_oCenterAtoms.m_baAtomType.GetSize();z2++)
						{
							for (z3=0;z3<((CxIntArray*)va->m_oCenterAtoms.m_oaAtoms[z2])->GetSize();z3++)
							{
								mprintf("      - %2s%-3d  Weight %.2f%c\n",((CAtom*)g_oaAtoms[va->m_oCenterAtoms.m_baRealAtomType[z2]])->m_sName,((CxIntArray*)va->m_oCenterAtoms.m_oaAtoms[z2])->GetAt(z3)+1,va->m_faWeight[z4]/va->m_fGesWeight*100.0f,'%');
								z4++;
							}
						}
					} else mprintf("        (pseudo-molecule, skipping atoms)\n");
				} else
				{
					mprintf("      - Distance #%d-%s%d = %.2f pm\n",z+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1,va->m_fValues[0]);
					mprintf("      - Angle #%d-%s%d-%s%d = %.2f degree\n",z+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[1]])->m_sName,va->m_iAtom[1]+1,va->m_fValues[1]*180.0/Pi);
					mprintf("      - Dihedral #%d-%s%d-%s%d-%s%d = %.2f degree\n",z+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[0]])->m_sName,va->m_iAtom[0]+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[1]])->m_sName,va->m_iAtom[1]+1,((CAtom*)g_oaAtoms[va->m_iRealAtomType[2]])->m_sName,va->m_iAtom[2]+1,va->m_fValues[2]*180.0/Pi);
				}
			} else if (va->m_iMode == 2)
			{
				mprintf(WHITE,"\n    #%d in %s - tip of molecular dipole vector starting in #2.\n",va->m_iMolVirtAtom+1,((CMolecule*)g_oaMolecules[va->m_iMolecule])->m_sName);
			}
		}
	}
	mprintf(GREEN,"\n<<< End of virtual atoms <<<\n\n");

	if (g_bPeriodic)
		g_TimeStep.UniteMolecules(false);
	
	g_TimeStep.CalcCenters();
	g_bFoldAtomwise = false;

	if (g_bPeriodic)
		g_TimeStep.FoldMolecules();
	
	BuildAtomIndices();
	//	g_TimeStep.WritePOV("step.pov");

	mprintf(WHITE,"\n>>> List of functions <<<\n\n");
	DumpAnalyses();
	
	mprintf("   (You may specify multiple analyses at once, but\n");
	mprintf("    the safe way is to perform only one kind of analysis at a time.)\n");
_fncinput:
	inpprintf("! Which functions to compute (comma separated)?\n");
	mprintf("\n    Which functions to compute (comma separated)? ");
	myget(&buf);
	mprintf("\n");
	if (!ParseFunctions(buf))
	{
		eprintf("Wrong input.\n");
		inpprintf("! Wrong input.\n");
		goto _fncinput;
	}
	mprintf("\n");

	if (g_bRegionAnalysis)
	{
		if (!(g_bRDF || g_bSDF || g_bADF || g_bDDF || g_bCDF || g_bRevSDF || g_bDens))
		{
			eprintf("    Region-specific analysis needs to be applied to some static analysis (e.g., RDF).\n\n");
			return false;
		}
		g_iaSMRegion.SetSize(g_oaSingleMolecules.GetSize());
		for (z=0;z<g_oaSingleMolecules.GetSize();z++)
			g_iaSMRegion[z] = (rand()%2)+1;
	}



	if (g_bNPT && (g_bACF || g_bMSD))
	{
		mprintf(RED,"    Warning: ");
		mprintf("When a variable cell vector is used, dynamical analyses (like spectra or MSDs)\n");
		mprintf("             will give erroneous results.\n\n");
	}

	if (g_bSaveJustTraj || g_bSaveRefEnv || g_bCutCluster)
		g_bNeedMoleculeWrap = true;

	if (g_bDomA || g_bVoro)
		g_bNeedMoleculeWrap = true;

	if (g_bBoxNonOrtho)
	{
		if (g_bDomA)
		{
			eprintf("    Domain Analysis currently only works with orthorhombic simulation cells.\n\n");
			return false;
		}
		if (g_bVoro)
		{
			eprintf("    Basic Voronoi analysis currently only works with orthorhombic simulation cells.\n\n");
			return false;
		}
	}

	if (g_bTegri)
		g_bNeedMoleculeWrap = true;

	mprintf(WHITE,"    The advanced mode includes some options which are quite powerful,\n    yet possibly weird or seldomly required.\n\n");

	g_bAdvanced2 = AskYesNo("    Use the advanced mode for the main part (y/n)? [no] ",false);
	mprintf("\n");

	if (g_bVDF || g_bRaman || g_bCombined || g_bBondACF || g_bVHDF || g_bACF || g_bUseVelocities || g_bUseForces || g_bMSD || g_bAggregation || g_bDLDisp || g_bDLDF || g_bDACF || g_bRDyn || g_bIRSpec || g_bTimeDiff || g_bDeriv || g_bNormalCoordinate || g_bPower || g_bIR || g_bSetUpPolarizabilityCalc || g_bFDF || g_bVCD || g_bMagneticDipoleRestart)
	{
		g_fTimestepLength = AskFloat("    Enter the length of one trajectory time step in fs: [0.5] ",0.5f);
		mprintf("\n");
	}

	/*********** Interface ***************/
	if (!Interface_BeforeAnalysis())
		return false;

	if (g_bCDF)
	{
		mprintf(YELLOW,"*** Combined Distribution Function\n\n");
		g_iCDFChannels = AskUnsignedInteger("    How many channels should the Combined Distribution Function have? [2] ",2);

		try { g_iObsChannel = new int[g_iCDFChannels]; } catch(...) { g_iObsChannel = NULL; }
		if (g_iObsChannel == NULL) NewException((double)g_iCDFChannels*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		mprintf("    Choose from these functions: rdf, adf, ddf, dip, vdf, pldf, lidf");
		mprintf("\n");
		for (z=0;z<g_iCDFChannels;z++)
		{
_entertype:
			AskString_ND("    Channel %d: Enter function type (e.g. rdf): ",&buf,z+1);
			if (mystricmp(buf,"rdf")==0)
			{
				g_iObsChannel[z] = 0;
				g_bRDF = true;
			} else if (mystricmp(buf,"adf")==0)
			{
				g_iObsChannel[z] = 1;
				g_bADF = true;
			} else if (mystricmp(buf,"ddf")==0)
			{
				g_iObsChannel[z] = 2;
				g_bDDF = true;
			} else if (mystricmp(buf,"dip")==0)
			{
				g_iObsChannel[z] = 3;
				g_bDipDF = true;
			} else if (mystricmp(buf,"vdf")==0)
			{
				g_iObsChannel[z] = 4;
				g_bVDF = true;
			} else if (mystricmp(buf,"pldf")==0)
			{
				g_iObsChannel[z] = 5;
				g_bPlDF = true;
			} else if (mystricmp(buf,"lidf")==0)
			{
				g_iObsChannel[z] = 6;
				g_bLiDF = true;
			} else
			{
				eprintf("Wrong input.\n");
				inpprintf("! Wrong input.\n");
				goto _entertype;
			}
		}
	} else g_iCDFChannels = 1;

	if (g_bDipDF || (g_bDipole && !g_bCHDF))
		ParseDipole();
	
	if (g_bDomA)
	{
		try { g_pDomainEngine = new CDomainEngine(); } catch(...) { g_pDomainEngine = NULL; }
		if (g_pDomainEngine == NULL) NewException(sizeof(CDomainEngine*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		g_pDomainEngine->Parse();
	}

	if (g_bVoro)
	{
		try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
		if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pVoroWrapper->Parse();
	}

	if (g_bTegri && (g_pTetraPak == NULL))
	{
		try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
		if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		try { g_pTetraPak = new CTetraPak(); } catch(...) { g_pTetraPak = NULL; }
		if (g_pTetraPak == NULL) NewException((double)sizeof(CTetraPak),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pTetraPak->Parse();
	}

/*	if (g_bSFac)
	{
		try { g_pSFac = new CStructureFactor(); } catch(...) { g_pSFac = NULL; }
		if (g_pSFac == NULL) NewException((double)sizeof(CStructureFactor),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pSFac->Parse();
	}*/


	if (g_bAggregation)
	{
		mprintf(WHITE,">>> Selection of Aggregation Functions >>>\n\n");

		g_bDACF = AskYesNo("    Compute Dimer Existence Autocorrelation Functions (DACFs) (y/n)? [yes] ",true);
		g_bDLDF = AskYesNo("    Compute Dimer Lifetime Distribution Functions (DLDFs) (y/n)? [yes] ",true);
		g_bDDisp = AskYesNo("    Compute Dimer Displacement Functions / Pair Diffusion (y/n)? [no] ",false);
		if (g_bDDisp)
		{
			g_bPairMSD = AskYesNo("    Compute Pair Mean Square Displacement (Pair Diffusion) (y/n)? [yes] ",true);
			g_bDLDisp = AskYesNo("    Compute combined Dimer Lifetime/Displacement Function (DLDisp) (y/n)? [yes] ",true);
		} else
		{
			g_bDLDisp = false;
			g_bPairMSD = false;
		}

		mprintf(WHITE,"\n<<< Selection of Aggregation Functions <<<\n\n");
	}

	if (g_bNbExchange)
	{
		mprintf(WHITE,">>> Selection of Neighborhood Exchange Functions >>>\n\n");
		mprintf(WHITE,"\n<<< Selection of Neighborhood Exchange Functions <<<\n\n");
	}
	
	if (g_bAdvanced2)
	{
		if (g_bRDyn || g_bIRSpec)
			g_bRDynCacheMode = AskYesNo("    Use RDyn cached mode (do this unless there are problems) (y/n)? [yes] ",true);

		if (g_bMSD)
			g_bMSDCacheMode = AskYesNo("    Use MSD cached mode (do this unless there are problems) (y/n)? [yes] ",true);
	} else
	{
		g_bRDynCacheMode = true;
		g_bMSDCacheMode = true;
	}


	if (g_bACF)
	{
//		mprintf(YELLOW,">>> Select Autocorrelation Functions >>>\n\n");

		if (!g_bPowerSpec)
			g_bVACF = AskYesNo("    Calculate velocity autocorrelation function / power spectra (y/n)? [yes] ",true);
				else g_bVACF = true;

		if (g_bVACF)
			g_bUseVelocities = true;

/*		if (g_bBetaFeatures)
			g_bBondACF = AskYesNo("    Calculate bond vibration spectra (y/n)? [no] ",false);
				else*/ g_bBondACF = false;

/*		g_bDipACF = AskYesNo("    Calculate dipole autocorrelation function / IR spectra (y/n)? [no] ",false);
		if (g_bDipACF)
			g_bDipole = true;*/

//		mprintf(YELLOW,"\n<<< End of Select Autocorrelation Functions <<<\n\n");
	}
	
	if (g_bBondACF)
	{
		mprintf(WHITE,">>> Bond vibration spectra >>>\n\n");
		if (g_iTrajSteps != -1)
			g_iBondACFDepth = AskUnsignedInteger("    Enter bond ACF depth in time steps: [%d] ",g_iTrajSteps/2,g_iTrajSteps/2);
				else g_iBondACFDepth = AskUnsignedInteger("    Enter bond ACF depth in time steps: [%d] ",4096,4096);
		ti = CalcFFTSize(g_iBondACFDepth,false);
		if (g_iBondACFDepth != ti)
		{
			mprintf(WHITE,"\n    The next \"fast\" size for FFT is %d. Using this instead of %d as size.\n",ti,g_iBondACFDepth);
			g_iBondACFDepth = ti;
		}
		g_bBondACFNormalize = AskYesNo("    Normalize all bond ACFs (removes intensity info for the peaks) (y/n)? [no] ",false);
		g_bBondACFSymmetrize = AskYesNo("    Symmetrize all bond ACFs (y/n)? [no] ",false);
		g_bBondACFWindow = AskYesNo("    Apply window function to bond ACFs (y/n)? [yes] ",true);
		g_bBondACFDebug = AskYesNo("    Write out bond ACF debug data (y/n)? [no] ",false);
		mprintf(WHITE,"\n<<< End of Bond vibration spectra <<<\n\n");
	}

	if (g_bRDyn || g_bIRSpec || g_bBondACF || g_bVACF)
	{
		if (g_bAdvanced2)
			g_bACFFFT = AskYesNo("    Use fourier transform for autocorrelation (much faster) (y/n)? [yes] ",true);
				else g_bACFFFT = true;
	} else g_bACFFFT = true;

	if (g_bCDF || g_bPlProj || g_bDipDF || g_bRDF || g_bVHDF || g_bSDF || g_bPlDF || g_bLiDF || 
      g_bRevSDF || g_bADF || g_bDDF || g_bSaveRefEnv || g_bCutCluster || g_bCond || g_bNbAnalysis)
	{
		if (g_oaMolecules.GetSize() > 1)
		{
//			sprintf(buf,"\n    Which of the molecules should be the reference molecule (");
			buf.sprintf("\n    Which of the molecules should be the reference molecule (");
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
//				sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//				strcat(buf,buf2);
//				if (z < g_oaMolecules.GetSize()-1)
//					strcat(buf,", ");
				buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
				buf.strcat(buf2);
				if (z < g_oaMolecules.GetSize()-1)
					buf.strcat(", ");
			}
//			strcat(buf,")? ");
			buf.strcat(")? ");
			g_iFixMol = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize()) - 1;
		} else g_iFixMol = 0;
		mprintf(WHITE,"\n    %s is the reference molecule.\n\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
	} else g_iFixMol = -1;

	if (g_bSaveRefEnv)
	{
		mprintf(YELLOW,">>> Save environment of the reference molecule >>>\n\n");
		g_bSaveRefWithEnv = AskYesNo("    Save the reference molecule itself (y/n)? [yes] ",true);
		g_bRefEnvCenter = AskYesNo("    Center the reference molecule in the box (y/n)? [yes] ",true);
		if (g_bRefEnvCenter)
		{
			mprintf("\n    The first reference atom will be put into the middle of the box.\n\n");
			g_bRefEnvFix = AskYesNo("    Fix rotational freedom of the reference molecule (y/n)? [no] ",false);
			if (g_bRefEnvFix)
			{
				mprintf("\n    The 2nd reference atom will be put onto the positive X axis,\n");
				mprintf("    and the 3rd reference atom into the X-Y plane with positive Y values.\n\n");
			}
			mprintf(WHITE,"    You will be asked for the reference atom(s) lateron.\n\n");
		}
		if (((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() > 1)
			g_iSaveRefMol = AskRangeInteger("    Which representative of the reference molecules to use (1-%d)? [1] ",1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize(),1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()) - 1;
				else g_iSaveRefMol = 0;
		mprintf(WHITE,"\nPlease choose mode for neighborhood detection:\n\n");
		inpprintf("! Please choose mode for neighborhood detection:\n");
		mprintf("  1.) Search neighbors once in the beginning and always show those.\n");
		mprintf("  2.) Search neighbors in every step and show the current neighbors.\n");
		mprintf("      (Warning: Neighbor count differs from step to step - be aware.\n");
		mprintf("  3.) Find the frequentiest neighbors over all steps and always show those.\n\n");
		mprintf(WHITE,"    Please note: ");
		mprintf("Programs such as VMD cannot handle variable atom count \n");
		mprintf("    along a trajectory and will therefore not work properly with mode 2!\n\n");
		g_iNbhMode = AskRangeInteger("    Choice (1-3): [3] ",1,3,3);
		mprintf(YELLOW,"\n<<< End of Save environment of the reference molecule <<<\n\n");
	}

	if (g_bCutCluster)
	{
		mprintf(WHITE,">>> Cut clusters >>>\n\n");
		g_bSaveRefWithEnv = AskYesNo("    Show reference molecule itself (y/n)? [yes] ",true);
		g_bRefEnvCenter = AskYesNo("    Center reference molecule in the box (y/n)? [yes] ",true);
		if (g_bRefEnvCenter)
			g_bRefEnvFix = AskYesNo("    Fix rotational freedom of the reference molecule (y/n)? [no] ",false);
_clustercount:
		g_iClusterCount = AskUnsignedInteger("    How many clusters to create? [100] ",100);
		g_iClusterSteps = AskUnsignedInteger("    How many time steps to use for cluster creation? [%d] ",g_iTrajSteps,g_iTrajSteps);
		mprintf("\nCreating cluster distribution...");
		for (z=0;z<g_iClusterCount;z++)
		{
			if (fmod(z,g_iClusterCount/25.0) < 1.0)
				mprintf(".");
			z3=0;
_clustagain:
			z3++;
			if (z3 > 500)
			{
				eprintf("Error: Too few molecules / time steps for requested cluster count.\n");
				goto _clustercount;
			}
			ti = ((((unsigned long)rand()%16384)+rand()*16384) % g_iClusterSteps) + 1;
			ti2 = rand()%((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize();
			for (z2=0;z2<g_iaClusterSteps.GetSize();z2++)
				if ((g_iaClusterSteps[z2] == ti) && (g_iaClusterMol[z2] == ti2))
					goto _clustagain;
			g_iaClusterSteps.Add(ti);
			g_iaClusterMol.Add(ti2);
		}
		mprintf("Done.\n");
		mprintf("Sorting cluster distribution...");
		for (z=0;z<g_iClusterCount-1;z++)
		{
			if (fmod(z,(g_iClusterCount-1)/25.0) < 1.0)
				mprintf(".");
			ti = 99999999;
			ti2 = 99999999;
			ti3 = -1;
			for (z2=z;z2<g_iClusterCount;z2++)
			{
				if (g_iaClusterSteps[z2] < ti)
				{
					ti = g_iaClusterSteps[z2];
					ti2 = g_iaClusterMol[z2];
					ti3 = z2;
				} else if (g_iaClusterSteps[z2] == ti)
				{
					if (g_iaClusterMol[z2] < ti2)
					{
						ti2 = g_iaClusterMol[z2];
						ti3 = z2;
					}
				}
			}
			if (ti3 == -1)
				abort();
			ti = g_iaClusterSteps[z];
			g_iaClusterSteps[z] = g_iaClusterSteps[ti3];
			g_iaClusterSteps[ti3] = ti;
			ti = g_iaClusterMol[z];
			g_iaClusterMol[z] = g_iaClusterMol[ti3];
			g_iaClusterMol[ti3] = ti;
		}
		mprintf("Done.\n\n");
		mprintf("Writing cluster selection to \"cluster.csv\"...");
		a = OpenFileWrite("cluster.csv",true);
		mfprintf(a,"# Cluster;  Timestep;  Ref. Mol.\n");

		for (z=0;z<g_iClusterCount;z++)
			mfprintf(a,"%d;  %d;  %d\n",z+1,g_iaClusterSteps[z],g_iaClusterMol[z]+1);

		fclose(a);
		mprintf("Done.\n\n");
		mprintf(WHITE,"<<< End of Cut clusters <<<\n\n");
	}

	if (g_bSDF || g_bPlProj || 
      g_bAvg || ((g_bSaveRefEnv || g_bCutCluster) && g_bRefEnvFix)) // Nur fuer SDFs: Welche zwei weiteren Atome fixieren?
	{
		g_iRefSystemDim = 3;
_ref3again:
		mprintf(WHITE,"    You have to choose three reference atoms from the reference molecule %s:\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		mprintf("    The first reference atom will be put into the center.\n");
		mprintf("    The 2nd reference atom will be put onto the positive X axis.\n");
		mprintf("    The 3rd reference atom will be put into the X-Y plane with positive Y values.\n\n");
		mprintf("    Please enter three comma-separated reference atoms (e.g. \"C1,H2,O1\"): ");
		inpprintf("! Please enter three comma-separated reference atoms (e.g. \"C1,H2,O1\"):\n");
		myget(&buf);
		if (!ParseRefSystem(g_iFixMol,buf,3))
			goto _ref3again;
		for (z=1;z<3;z++)
		{
			for (z2=0;z2<z;z2++)
			{
				if ((g_iFixAtomType[z] == g_iFixAtomType[z2]) && (g_iFixAtom[z] == g_iFixAtom[z2]))
				{
					eprintf("Please enter three pairwise different atoms.\n\n");
					inpprintf("! Please enter three pairwise different atoms.\n");
					goto _ref3again;
				}
			}
		}
 		mprintf("\n");
		mprintf(WHITE,"    Reference plane: Fixing in %s the %d. %s-, the %d. %s- and the %d. %s atom.\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[2]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName);
		mprintf("\n");
	} else if (g_bRevSDF)
	{
		g_iRefSystemDim = 2;
_ref2again:
		mprintf("    Please enter two reference atoms (e.g. C1,H3): ");
		inpprintf("! Please enter two reference atoms (e.g. C1,H3):\n");
		myget(&buf);
		if (!ParseRefSystem(g_iFixMol,buf,2))
			goto _ref2again;
		if ((g_iFixAtomType[0] == g_iFixAtomType[1]) && (g_iFixAtom[0] == g_iFixAtom[1]))
		{
			eprintf("Please enter two different atoms.\n\n");
			inpprintf("! Please enter two different atoms.\n");
			goto _ref2again;
		}
 		mprintf("\n");
		mprintf(WHITE,"Reference axis: Fixing in %s the %d. %s and the %d. %s atom.\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName);
		mprintf("\n");
	} else if (g_bVHDF || g_bRDF || g_bADF || g_bPlDF || g_bLiDF || g_bDDF || g_bSaveRefEnv || g_bCutCluster) // Fuer Verteilungsfunktion: Welches Atom im Ursprung
	{
		g_iRefSystemDim = 1;
		if (g_bAdvanced2)
		{
_ref1again:
			mprintf("    Please enter the atom to put into the box center (e.g. C3): [center of mass] ");
			inpprintf("! Please enter the atom to put into the box center (e.g. C3): [center of mass]\n");
			myget(&buf);
			if (strlen(buf)==0)
			{
				if (!ParseRefSystem(g_iFixMol,"#2",1))
				{
					eprintf("Weird error.\n");
					inpprintf("! Weird error.\n");
					abort();
				}
			} else if (!ParseRefSystem(g_iFixMol,buf,1))
				goto _ref1again;
 			mprintf("\n");
			mprintf(WHITE,"Reference atom: Fixing in %s the %d. %s atom.\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName);
			mprintf("\n");
		} else
		{
			if (!ParseRefSystem(g_iFixMol,"#2",1))
			{
				eprintf("Weird error.\n");
				inpprintf("! Weird error.\n");
				abort();
			}
		}
	} else if (g_iFixMol != -1)
	{
		g_iRefSystemDim = 1;
		if (!ParseRefSystem(g_iFixMol,"#2",1))
		{
			eprintf("Weird error.\n");
			abort();
		}
		mprintf("    Reference atom: Fixing in %s the %d. %s atom.\n\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName);
	}

	if (g_bSaveRefEnv || g_bCutCluster)
	{
		mprintf(WHITE,">>> Neighborhood Definition >>>\n\n");

		try { g_pNbSet = new CNbSet(); } catch(...) { g_pNbSet = NULL; }
		if (g_pNbSet == NULL) NewException((double)sizeof(CNbSet),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pNbSet->Parse(g_iFixMol);
		mprintf(WHITE,"\n<<< End of Neighborhood Definition <<<\n\n");

		if (g_bSaveRefEnv)
		{
			if (AskYesNo("    Create a temporal development overlay (TDO) plot (y/n)? [no] ",false))
			{
				g_bTDO = true;
				g_laTDOSteps.RemoveAll();
				if (AskYesNo("    Use equidistant intervals (y) or specify each point (n) (y/n)? [yes] ",true))
				{
					g_bTDOEqui = true;
					g_iTDOCount = AskUnsignedInteger("    Enter the overlay count: [5] ",5);
					g_iTDOStart = AskUnsignedInteger("    Enter the starting timestep: [0] ",0);
					g_iTDOStride = AskUnsignedInteger("    Enter the overlay distance in steps: [1000] ",1000);
					for (z=0;z<g_iTDOCount;z++)
						g_laTDOSteps.Add(g_iTDOStart+z*g_iTDOStride);
				} else
				{
					g_bTDOEqui = false;
					while (true)
					{
						ti = AskUnsignedInteger("    Enter timestep for %dth overlay: [done] ",999999,g_laTDOSteps.GetSize()+1);
						if (ti == 999999)
							break;
						g_laTDOSteps.Add(ti);
					}
				}
				g_fTDOBleaching = AskRangeFloat("    Enter the TDO bleaching grade (0-1): [0.8] ",0.0f,1.0f,0.8f);
			} else g_bTDO = false;
		}
	}

	if (g_bPlProj || g_bAggregation || g_bNbAnalysis || g_bPlDF || g_bLiDF || g_bRDyn || g_bIRSpec || g_bDipDF || g_bMSD || g_bADF || g_bDDF || g_bVACF || g_bDipACF || g_bRevSDF || g_bSDF || g_bVHDF || g_bRDF || g_bDens || g_bVDF || /*g_bFDF ||*/ g_bCond) // Fuer Verteilungsfunktion: Welches Atom beobachten?
	{
		if (g_bVACF)
		{
			if (g_bAdvanced2)
				g_bVACFCacheMode = AskYesNo("    Use VACF cached mode (do this unless there are problems) (y/n)? [yes] ",true);
					else g_bVACFCacheMode = true;

			if (g_bPowerSpec)
			{
				g_bGlobalVACF = AskYesNo("    Compute power spectrum of whole system (y/n)? [yes] ",true);

				if (g_bGlobalVACF)
					if (!AskYesNo("    Compute also power spectra for certain atoms/molecules (y/n)? [no] ",false))
						goto _endobs;
			} else
			{
				g_bGlobalVACF = AskYesNo("    Compute global velocity ACF of whole system (y/n)? [yes] ",true);

				if (g_bGlobalVACF)
					if (!AskYesNo("    Compute also velocity ACFs for certain atoms/molecules (y/n)? [no] ",false))
						goto _endobs;
			}
		}

		if (g_bIRSpec)
		{
			g_bGlobalIR = AskYesNo("    Compute IR spectrum of whole system (y/n)? [yes] ",true);
			if (g_bGlobalIR)
				if (!AskYesNo("    Compute also IR spectra for certain molecule types (y/n)? [no] ",false))
					goto _endobs;
		}

/*		if (g_bDipACF)
		{
			g_bGlobalDipACF = AskYesNo("    Create global dipole ACF of all atoms (y/n)? [yes] ",true);
		}*/
_nextsdf:
		mprintf(YELLOW,"\n>>> Observation %d >>>\n\n",g_oaObserv.GetSize()+1);

		try { o = new CObservation(); } catch(...) { o = NULL; }
		if (o == NULL) NewException((double)sizeof(CObservation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		o->m_pConditions = NULL;
		g_oaObserv.Add(o);
		o->m_bTimeDev = false;
		if (g_bAggregation || g_bNbExchange)
		{
			o->m_bSelf = false;
			o->m_bOthers = false;
		} else if (g_bNbAnalysis)
		{
			o->m_bSelf = false;
			o->m_bOthers = true;
			mprintf("    Performing this observation intermolecular.\n\n");
		} else if (g_iFixMol != -1)
		{
			if (AskRangeInteger("    Perform this observation intramolecular (within the reference molecule) (0) or intermolecular (1)? [1] ",0,1,1) == 1)
			{
				o->m_bSelf = false;
				o->m_bOthers = true;
			} else
			{
				o->m_bSelf = true;
				o->m_bOthers = false;
			}
			mprintf("\n");
		} else 
		{
			o->m_bSelf = false;
			o->m_bOthers = true;
		}
		if (o->m_bOthers)
		{
			if (g_bCDF)
			{
				mprintf(WHITE,"    Please note: ");
				mprintf("Although you specified this observation to be intermolecular, you may of course\n");
				mprintf("                 choose atoms from RM only (or OM only), which will yield an intramolecular quantity.\n\n");
			}
			if (g_bCDF)
				o->m_bSecondShowMol = AskYesNo("    CDF: Perform a three-body analysis (y) or observe one molecule at a time (n)? [no] ",false);
					else o->m_bSecondShowMol = false;
			if (o->m_bSecondShowMol)
			{
				mprintf(WHITE,"\n    You have to select two observed molecules (OMs):\n");
				mprintf("    The first CDF channel observes the 1st OM, the second channel observes the 2nd OM.\n\n");
			}
			if (g_oaMolecules.GetSize() > 1)
			{
_obsmolagain:
				if (o->m_bSecondShowMol)
//					sprintf(buf,"    Which 1st molecule should be observed (");
					buf.sprintf("    Which 1st molecule should be observed (");
				else
//					sprintf(buf,"    Which molecule should be observed (");
					buf.sprintf("    Which molecule should be observed (");

				for (z=0;z<g_oaMolecules.GetSize();z++)
				{
//					sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//					strcat(buf,buf2);
//					if (z < g_oaMolecules.GetSize()-1)
//						strcat(buf,", ");
					buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
					buf.strcat(buf2);
					if (z < g_oaMolecules.GetSize()-1)
						buf.strcat(", ");
				}
//				strcat(buf,")? ");
				buf.strcat(")? ");

				o->m_iShowMol = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize()) - 1;
				if ((((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize() == 1) && (o->m_iShowMol == g_iFixMol))
				{
					eprintf("Error: There is only 1 molecule of %s.\n--> Intermolecular observation between %s and %s not possible.\n\n",((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName);
					goto _obsmolagain;
				}
			} else 
			{
				mprintf("    Only one molecule type, choosing %s as observed molecule (OM).\n",((CMolecule*)g_oaMolecules[0])->m_sName);
				o->m_iShowMol = 0;
			}
			o->m_iShowMolCount = ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();
			if (o->m_bSecondShowMol)
			{
				if (g_oaMolecules.GetSize() > 1)
				{
_obsmol2again:
//					sprintf(buf,"    Which 2nd molecule should be observed (");
					buf.sprintf("    Which 2nd molecule should be observed (");
					for (z=0;z<g_oaMolecules.GetSize();z++)
					{
//						sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//						strcat(buf,buf2);
//						if (z < g_oaMolecules.GetSize()-1)
//							strcat(buf,", ");
						buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
						buf.strcat(buf2);
						if (z < g_oaMolecules.GetSize()-1)
							buf.strcat(", ");
					}
//					strcat(buf,")? ");
					buf.strcat(")? ");

					o->m_iShowMol2 = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize()) - 1;
					if ((((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_laSingleMolIndex.GetSize() == 1) && (o->m_iShowMol2 == g_iFixMol))
					{
						eprintf("Error: There is only 1 molecule of %s.\n  --> Intermolecular observation between %s and %s not possible.\n\n",((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_sName);
						goto _obsmol2again;
					}
				} else o->m_iShowMol2 = 0;
				o->m_iShowMol2Count = ((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_laSingleMolIndex.GetSize();
				if (o->m_iShowMol == o->m_iShowMol2)
					o->m_bExclude1eq2 = AskYesNo("    Exclude the case \"1st OM = 2nd OM\" (y/n)? [yes] ",true);
						else o->m_bExclude1eq2 = false;
			}
		} else 
		{
			o->m_iShowMol = -1;
			o->m_iShowMol2 = -1;
			o->m_bSecondShowMol = false;
			o->m_iShowMolCount = 1;
		}

//		buf[0] = 0;
		buf.sprintf("");

		if ((o->m_bSecondShowMol) && (o->m_iShowMol != o->m_iShowMol2))
//			sprintf(buf," of %s / %s / %s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_sName);
			buf.sprintf(" of %s / %s / %s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_sName);
		else if ((o->m_iShowMol != -1) && (g_iFixMol != -1))
//			sprintf(buf," of %s / %s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName);
			buf.sprintf(" of %s / %s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName);
		else if (g_iFixMol != -1)
//			sprintf(buf," of %s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
			buf.sprintf(" of %s",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);

		if (g_iFixMol != -1)
		{
			if (g_bAdvanced2)
			{
				if (AskYesNo("    Observe only certain molecules%s (y/n)? [no] ",false,(const char*)buf))
				{
					o->m_bObsCertain = true;
					if (!g_bCDF && (g_bRDF || g_bADF || g_bDDF || g_bDipDF))
						o->m_bDecompDist = AskYesNo("    Decompose distribution functions into contributions from molecules (y/n)? [yes] ",true);
					
					if (((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() > 1)
					{
						mprintf("    Which %s molecules (RM) to take into account (e.g. 1,3-7)? [all] ",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
						inpprintf("! Which %s molecules (RM) to take into account (e.g. 1,3-7)? [all]\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
						myget(&buf);
						if (strlen(buf)==0)
						{
							for (z=0;z<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
								o->m_waObsRefList.Add(z);
						} else
						{
							ParseIntList(buf,&o->m_waObsRefList);
							for (z=0;z<o->m_waObsRefList.GetSize();z++)
								o->m_waObsRefList[z]--;
						}
					} else o->m_waObsRefList.Add(0);

				//	mprintf("RefList[0] = %d\n",o->m_waObsRefList[0]);

					if (o->m_iShowMol != -1)
					{
						if (((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize() > 1)
						{
							mprintf("    Which %s molecules (OM) to take into account (e.g. 1,3-7)? [all] ",((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName);
							inpprintf("! Which %s molecules (OM) to take into account (e.g. 1,3-7)? [all]\n",((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName);
							myget(&buf);
							if (strlen(buf)==0)
							{
								for (z=0;z<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
									o->m_waObsShowList.Add(z);
							} else
							{
								ParseIntList(buf,&o->m_waObsShowList);
								for (z=0;z<o->m_waObsShowList.GetSize();z++)
									o->m_waObsShowList[z]--;
							}
						} else o->m_waObsShowList.Add(0);
					} else o->m_waObsShowList.Add(0); // Dummy

					if (o->m_bSecondShowMol)
					{
						if (o->m_iShowMol2 != -1)
						{
							if (((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_laSingleMolIndex.GetSize() > 1)
							{
								mprintf("    Which %s molecules (2nd OM) to take into account (e.g. 1,3-7)? [all] ",((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_sName);
								inpprintf("! Which %s molecules (2nd OM) to take into account (e.g. 1,3-7)? [all]\n",((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_sName);
								myget(&buf);
								if (strlen(buf)==0)
								{
									for (z=0;z<((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_laSingleMolIndex.GetSize();z++)
										o->m_waObsShow2List.Add(z);
								} else
								{
									ParseIntList(buf,&o->m_waObsShow2List);
									for (z=0;z<o->m_waObsShow2List.GetSize();z++)
										o->m_waObsShow2List[z]--;
								}
							} else o->m_waObsShow2List.Add(0);
						} else o->m_waObsShow2List.Add(0); // Dummy
					}
				} else // if not onlysome
				{
					o->m_bObsCertain = false;
					o->m_bDecompDist = false;
				}
			} else // if not advanced
			{
				o->m_bObsCertain = false;
				o->m_bDecompDist = false;
			}
		} else
		{
			o->m_bObsCertain = false;
			o->m_bDecompDist = false;
		}

		o->m_bDecompType = false;

		if ((!o->m_bObsCertain) && (!g_bCDF) && g_bAdvanced2 && (g_bRDF || g_bADF || g_bDDF || g_bDipDF))
		{
			if (AskYesNo("    Decompose this observation into contributions from different elements (y/n)? [no] ",false))
			{
				o->m_bDecompType = true;
			}
		}

		if (g_bRegionAnalysis)
		{
			AskString("    Take reference molecules from which regions? [all] ",&buf,"0");
			ParseIntList(buf,&o->m_iaRMRegions);
	//		for (z6=0;z6<o->m_iaRMRegions.GetSize();z6++)
	//			if (o->m_iaRMRegions[z6] != 0)
	//				o->m_iaRMRegions[z6]--;

			if (o->m_bOthers)
			{
				if (o->m_bSecondShowMol)
				{
					AskString("    Take 1st observed molecules from which regions? [all] ",&buf,"0");
					ParseIntList(buf,&o->m_iaOM1Regions);
				//	for (z6=0;z6<o->m_iaOM1Regions.GetSize();z6++)
				//		if (o->m_iaOM1Regions[z6] != 0)
				//			o->m_iaOM1Regions[z6]--;
					AskString("    Take 2nd observed molecules from which regions? [all] ",&buf,"0");
					ParseIntList(buf,&o->m_iaOM2Regions);
				//	for (z6=0;z6<o->m_iaOM2Regions.GetSize();z6++)
				//		if (o->m_iaOM2Regions[z6] != 0)
				//			o->m_iaOM2Regions[z6]--;
				} else
				{
					AskString("    Take observed molecules from which regions? [all] ",&buf,"0");
					ParseIntList(buf,&o->m_iaOM1Regions);
				//	for (z6=0;z6<o->m_iaOM1Regions.GetSize();z6++)
				//		if (o->m_iaOM1Regions[z6] != 0)
				//			o->m_iaOM1Regions[z6]--;
				}
			}
		}

		if (g_bAggregation || g_bNbExchange)
		{
			try { o->m_pDACF = new CDACF(); } catch(...) { o->m_pDACF = NULL; }
			if (o->m_pDACF == NULL) NewException((double)sizeof(CDACF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pDACF->Parse();
		}

		if (g_bRDyn)
		{
			try { o->m_pRDyn = new CReorDyn(); } catch(...) { o->m_pRDyn = NULL; }
			if (o->m_pRDyn == NULL) NewException((double)sizeof(CReorDyn),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pRDyn->m_iShowMol = o->m_iShowMol;
			o->m_pRDyn->Parse();
		}

		if (g_bIRSpec)
		{
			try { o->m_pIRSpec = new CReorDyn(); } catch(...) { o->m_pIRSpec = NULL; }
			if (o->m_pIRSpec == NULL) NewException((double)sizeof(CReorDyn),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pIRSpec->m_iShowMol = o->m_iShowMol;
			o->m_pIRSpec->ParseSpec();
		}

		if (g_bDens)
		{
			try { o->m_pDensityDF = new CDensDF(); } catch(...) { o->m_pDensityDF = NULL; }
			if (o->m_pDensityDF == NULL) NewException((double)sizeof(CDensDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pDensityDF->m_iShowMol = o->m_iShowMol;
			o->m_pDensityDF->Parse();
		}
			
		if (g_bSDF)
		{
			try { o->m_pSDF = new CSDF(); } catch(...) { o->m_pSDF = NULL; }
			if (o->m_pSDF == NULL) NewException((double)sizeof(CSDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			 
			o->m_pSDF->m_bIntra = o->m_bSelf;
			o->m_pSDF->m_iShowMol = o->m_iShowMol;
			o->m_pSDF->Parse(false);
		}

		if (g_bPlProj)
		{
			try { o->m_pPlProj = new CPlProj(); } catch(...) { o->m_pPlProj = NULL; }
			if (o->m_pPlProj == NULL) NewException((double)sizeof(CPlProj),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			 
			o->m_pPlProj->m_bIntra = o->m_bSelf;
			o->m_pPlProj->m_iShowMol = o->m_iShowMol;
			o->m_pPlProj->Parse();
		}

		if (g_bRevSDF)
		{
			try { o->m_pRevSDF = new CRevSDF();  } catch(...) { o->m_pRevSDF = NULL; }
			if (o->m_pRevSDF == NULL) NewException((double)sizeof(CRevSDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pRevSDF->m_bIntra = o->m_bSelf;
			o->m_pRevSDF->m_iShowMol = o->m_iShowMol;
			o->m_pRevSDF->Parse();
		}

		if (g_bNbAnalysis)
		{
			try { o->m_pNbAnalysis = new CNbAnalysis();  } catch(...) { o->m_pNbAnalysis = NULL; }
			if (o->m_pNbAnalysis == NULL) NewException((double)sizeof(CNbAnalysis),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pNbAnalysis->m_iShowMol = o->m_iShowMol;
			o->m_pNbAnalysis->Parse();
		}

		try { o->m_pRDF = new CRDF*[g_iCDFChannels];  } catch(...) { o->m_pRDF = NULL; }
		if (o->m_pRDF == NULL) NewException((double)g_iCDFChannels*sizeof(CRDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { o->m_pADF = new CADF*[g_iCDFChannels];  } catch(...) { o->m_pADF = NULL; }
		if (o->m_pADF == NULL) NewException((double)g_iCDFChannels*sizeof(CADF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { o->m_pDDF = new CDDF*[g_iCDFChannels];  } catch(...) { o->m_pDDF = NULL; }
		if (o->m_pDDF == NULL) NewException((double)g_iCDFChannels*sizeof(CDDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { o->m_pDipDF = new CDipDF*[g_iCDFChannels];  } catch(...) { o->m_pDipDF = NULL; }
		if (o->m_pDipDF == NULL) NewException((double)g_iCDFChannels*sizeof(CDipDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { o->m_pVDF = new CVDF*[g_iCDFChannels];  } catch(...) { o->m_pVDF = NULL; }
		if (o->m_pVDF == NULL) NewException((double)g_iCDFChannels*sizeof(CVDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		try { o->m_pPlDF = new CPlDF*[g_iCDFChannels];  } catch(...) { o->m_pPlDF = NULL; }
		if (o->m_pPlDF == NULL) NewException((double)g_iCDFChannels*sizeof(CVDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		try { o->m_pLiDF = new CLiDF*[g_iCDFChannels];  } catch(...) { o->m_pLiDF = NULL; }
		if (o->m_pLiDF == NULL) NewException((double)g_iCDFChannels*sizeof(CVDF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		for (z=0;z<g_iCDFChannels;z++)
		{
			o->m_pRDF[z] = NULL;
			o->m_pADF[z] = NULL;
			o->m_pDDF[z] = NULL;
			o->m_pDipDF[z] = NULL;
			o->m_pVDF[z] = NULL;
			o->m_pPlDF[z] = NULL;
			o->m_pLiDF[z] = NULL;
		}

		if (g_bCDF)
		{
			try { o->m_pCDF = new CCDF();  } catch(...) { o->m_pCDF = NULL; }
			if (o->m_pCDF == NULL) NewException((double)sizeof(CCDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			try { o->m_pCDF->m_iCombinations = new int[g_iCDFChannels];  } catch(...) { o->m_pCDF->m_iCombinations = NULL; }
			if (o->m_pCDF->m_iCombinations == NULL) NewException((double)g_iCDFChannels*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			mprintf(BLUE,"\n>>> Combined Distribution Function >>>\n");
//			o->m_pCDF->m_bChannelAll = new bool[g_iCDFChannels];
			o->m_pCDF->m_iCombinationProd = 1;
			for (z=0;z<g_iCDFChannels;z++)
			{
				mprintf(BLUE,"\n### Channel %d ###\n",z+1);
				switch(g_iObsChannel[z])
				{
					case 0:
						try { o->m_pRDF[z] = new CRDF(); } catch(...) { o->m_pRDF[z] = NULL; }
						if (o->m_pRDF[z] == NULL) NewException((double)sizeof(CRDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						if (o->m_bSecondShowMol && (z == 1))
							o->m_pRDF[z]->m_iShowMol = o->m_iShowMol2;
								else o->m_pRDF[z]->m_iShowMol = o->m_iShowMol;
						o->m_pRDF[z]->Parse();
						o->m_pCDF->m_iCombinations[z] = o->m_pRDF[z]->m_iCombinations;
						o->m_pCDF->m_iCombinationProd *= o->m_pRDF[z]->m_iCombinations;
/*						if (o->m_pRDF[z]->m_bSaveDist)
							o->m_bTimeDev = true;*/
						break;

					case 1:
						try { o->m_pADF[z] = new CADF(); } catch(...) { o->m_pADF[z] = NULL; }
						if (o->m_pADF[z] == NULL) NewException((double)sizeof(CADF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						if (o->m_bSecondShowMol && (z == 1))
							o->m_pADF[z]->m_iShowMol = o->m_iShowMol2;
								else o->m_pADF[z]->m_iShowMol = o->m_iShowMol;
						o->m_pADF[z]->Parse();
						o->m_pCDF->m_iCombinations[z] = o->m_pADF[z]->m_iCombinations;
						o->m_pCDF->m_iCombinationProd *= o->m_pADF[z]->m_iCombinations;
/*						if (o->m_pADF[z]->m_bSaveAngle)
							o->m_bTimeDev = true;*/
						break;

					case 2:
						try { o->m_pDDF[z] = new CDDF(); } catch(...) { o->m_pDDF[z] = NULL; }
						if (o->m_pDDF[z] == NULL) NewException((double)sizeof(CDDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						if (o->m_bSecondShowMol && (z == 1))
							o->m_pDDF[z]->m_iShowMol = o->m_iShowMol2;
								else o->m_pDDF[z]->m_iShowMol = o->m_iShowMol;
						o->m_pDDF[z]->Parse();
						o->m_pCDF->m_iCombinations[z] = o->m_pDDF[z]->m_iCombinations;
						o->m_pCDF->m_iCombinationProd *= o->m_pDDF[z]->m_iCombinations;
/*						if (o->m_pDDF[z]->m_bSaveAngle)
							o->m_bTimeDev = true;*/
						break;

					case 3:
						try { o->m_pDipDF[z] = new CDipDF(); } catch(...) { o->m_pDipDF[z] = NULL; }
						if (o->m_pDipDF[z] == NULL) NewException((double)sizeof(CDipDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						if (o->m_bSecondShowMol && (z == 1))
							o->m_pDipDF[z]->m_iShowMol = o->m_iShowMol2;
								else o->m_pDipDF[z]->m_iShowMol = o->m_iShowMol;
						o->m_pDipDF[z]->m_iCombinations = 1;
						o->m_pDipDF[z]->Parse();
						o->m_pCDF->m_iCombinations[z] = 1;
/*						if (o->m_pDipDF[z]->m_bSaveDipole)
							o->m_bTimeDev = true;*/
						break;

					case 4:
						try { o->m_pVDF[z] = new CVDF(); } catch(...) { o->m_pVDF[z] = NULL; }
						if (o->m_pVDF[z] == NULL) NewException((double)sizeof(CVDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						 
						if (o->m_bSecondShowMol && (z == 1))
							o->m_pVDF[z]->m_iShowMol = o->m_iShowMol2;
								else o->m_pVDF[z]->m_iShowMol = o->m_iShowMol;
						o->m_pVDF[z]->Parse();
						o->m_pCDF->m_iCombinations[z] = o->m_pVDF[z]->m_iCombinations;
						o->m_pCDF->m_iCombinationProd *= o->m_pVDF[z]->m_iCombinations;
						g_bUseVelocities = true;
/*						if (o->m_pVDF[z]->m_bSaveSpeed)
							o->m_bTimeDev = true;*/
						break;

					case 5:
						try { o->m_pPlDF[z] = new CPlDF(); } catch(...) { o->m_pPlDF[z] = NULL; }
						if (o->m_pPlDF[z] == NULL) NewException((double)sizeof(CPlDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						if (o->m_bSecondShowMol && (z == 1))
							o->m_pPlDF[z]->m_iShowMol = o->m_iShowMol2;
								else o->m_pPlDF[z]->m_iShowMol = o->m_iShowMol;
						o->m_pPlDF[z]->Parse();
						o->m_pCDF->m_iCombinations[z] = o->m_pPlDF[z]->m_iCombinations;
						o->m_pCDF->m_iCombinationProd *= o->m_pPlDF[z]->m_iCombinations;
						break;

					case 6:
						try { o->m_pLiDF[z] = new CLiDF(); } catch(...) { o->m_pLiDF[z] = NULL; }
						if (o->m_pLiDF[z] == NULL) NewException((double)sizeof(CLiDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						if (o->m_bSecondShowMol && (z == 1))
							o->m_pLiDF[z]->m_iShowMol = o->m_iShowMol2;
								else o->m_pLiDF[z]->m_iShowMol = o->m_iShowMol;
						o->m_pLiDF[z]->Parse();
						o->m_pCDF->m_iCombinations[z] = o->m_pLiDF[z]->m_iCombinations;
						o->m_pCDF->m_iCombinationProd *= o->m_pLiDF[z]->m_iCombinations;
						break;
				}
/*				if (z != 0)
					o->m_pCDF->m_bChannelAll[z] = AskYesNo("    Should this CDF Channel observe all OMs (y) or only the current one (n)? [no] ",false);
						else o->m_pCDF->m_bChannelAll[z] = false;*/
			}
/*			mprintf("    Write out temporal development for this CDF (1=yes,0=no)? [no] ");
			myget(buf);
			o->m_pCDF->m_bTimeDev = (atoi(buf)!=0);
			if (o->m_pCDF->m_bTimeDev)
				o->m_bTimeDev = true;*/

			mprintf("\n");

			try { o->m_pCDF->m_iResolution = new int[g_iCDFChannels]; } catch(...) { o->m_pCDF->m_iResolution = NULL; }
			if (o->m_pCDF->m_iResolution == NULL) NewException((double)g_iCDFChannels*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z=0;z<g_iCDFChannels;z++)
			{
				switch(g_iObsChannel[z])
				{
//					case 0: sprintf(buf,"RDF"); break;
//					case 1: sprintf(buf,"ADF"); break;
//					case 2: sprintf(buf,"DDF"); break;
//					case 3: sprintf(buf,"DipDF"); break;
//					case 4: sprintf(buf,"VDF"); break;
//					case 5: sprintf(buf,"PlDF"); break;
//					case 6: sprintf(buf,"LiDF"); break;
					case 0: buf.sprintf("RDF"); break;
					case 1: buf.sprintf("ADF"); break;
					case 2: buf.sprintf("DDF"); break;
					case 3: buf.sprintf("DipDF"); break;
					case 4: buf.sprintf("VDF"); break;
					case 5: buf.sprintf("PlDF"); break;
					case 6: buf.sprintf("LiDF"); break;
				}
				o->m_pCDF->m_iResolution[z] = AskUnsignedInteger("    Please enter the resolution (bin count) for CDF channel %d (%s): [100] ",100,z+1,(const char*)buf);
			}

			if (g_bAdvanced2)
				o->m_pCDF->m_iHistogramRes = AskUnsignedInteger("    Please enter CDF histogram resolution (0=no histogram): [5000] ",5000);
					else o->m_pCDF->m_iHistogramRes = 0;

			if (g_iCDFChannels == 2)
			{
				try { o->m_pCDF->m_pCombineList = new char[o->m_pCDF->m_iCombinationProd]; } catch(...) { o->m_pCDF->m_pCombineList = NULL; }
				if (o->m_pCDF->m_pCombineList == NULL) NewException((double)o->m_pCDF->m_iCombinationProd*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				if ((o->m_pCDF->m_iCombinations[0] > 1) && (o->m_pCDF->m_iCombinations[1] > 1))
				{
					mprintf("\n");
					for (z=0;z<g_iCDFChannels;z++)
					{
						mprintf(WHITE,"CDF channel %d has the following %d observations:\n",z+1,o->m_pCDF->m_iCombinations[z]);
						o->ListCDFObservations(z);
						mprintf("\n");
					}
					if (o->m_pCDF->m_iCombinations[0] == o->m_pCDF->m_iCombinations[1])
					{
						if (AskYesNo("    Combine n-th with n-th observation? (y/n) [yes] ",true))
						{
							for (z=0;z<o->m_pCDF->m_iCombinations[0];z++)
								for (z2=0;z2<o->m_pCDF->m_iCombinations[1];z2++)
									o->m_pCDF->m_pCombineList[z*o->m_pCDF->m_iCombinations[1]+z2] = (z==z2)?1:0;
							goto _combdone;
						} goto _askcomball;
					} else
					{
_askcomball:
						if (!AskYesNo("    Combine each with each observation (y), or use only some combinations (n)? [yes] ",true))
						{
							for (z=0;z<o->m_pCDF->m_iCombinations[0];z++)
								for (z2=0;z2<o->m_pCDF->m_iCombinations[1];z2++)
									o->m_pCDF->m_pCombineList[z*o->m_pCDF->m_iCombinations[1]+z2] = 0;

							mprintf("\n    Please enter all the combinations you want to observe.\n");
							mprintf("    Enter each combination as a comma-separated %d-tuple (e.g. ",g_iCDFChannels);
							for (z=0;z<g_iCDFChannels;z++)
							{
								mprintf("%d",z+1);
								if (z+1 < g_iCDFChannels)
									mprintf(",");
							}
							mprintf(").\n\n");
_nextcomb:
							AskString("    Enter combination (return=finished): ",&buf,"");
							if (strlen(buf)==0)
								goto _combdone;
							tempwa.RemoveAll_KeepSize();
							ParseIntList(buf,&tempwa);
							if (tempwa.GetSize() != g_iCDFChannels)
							{
								eprintf("    Wrong input, %d instead of %d values.\n",tempwa.GetSize(),g_iCDFChannels);
								goto _nextcomb;
							}
							for (z=0;z<g_iCDFChannels;z++)
							{
								if ((tempwa[z] < 1) || (tempwa[z] > o->m_pCDF->m_iCombinations[z]))
								{
									eprintf("    Wrong input, channel %d has only %d observations (%d requested).\n",z+1,o->m_pCDF->m_iCombinations[z],tempwa[z]);
									goto _nextcomb;
								}
								tempwa[z]--;
							}
							if (o->m_pCDF->m_pCombineList[tempwa[0]*o->m_pCDF->m_iCombinations[1]+tempwa[1]] == 1)
							{
								eprintf("    This combination has already been added.\n");
								goto _nextcomb;
							}
							o->m_pCDF->m_pCombineList[tempwa[0]*o->m_pCDF->m_iCombinations[1]+tempwa[1]] = 1;
							goto _nextcomb;
						} else goto _combineall;
_combdone:;
					}
				} else
				{
_combineall:
					for (z=0;z<o->m_pCDF->m_iCombinationProd;z++)
						o->m_pCDF->m_pCombineList[z] = 1;
				}

				o->m_pCDF->m_iCombinationsEnabled = 0;
				for (z=0;z<o->m_pCDF->m_iCombinationProd;z++)
					if (o->m_pCDF->m_pCombineList[z] != 0)
						o->m_pCDF->m_iCombinationsEnabled++;

				if (o->m_pCDF->m_iCombinationsEnabled == 1)
				{
					mprintf(WHITE,"\n    Using 1 combination for each RM-OM pair.\n\n");
				} else
				{
					mprintf(WHITE,"\n    Using %d combinations for each RM-OM pair:\n\n",o->m_pCDF->m_iCombinationsEnabled);
					z3 = 0;
					for (z=0;z<o->m_pCDF->m_iCombinations[0];z++)
					{
						for (z2=0;z2<o->m_pCDF->m_iCombinations[1];z2++)
						{
							if (o->m_pCDF->m_pCombineList[z*o->m_pCDF->m_iCombinations[1]+z2] == 1)
							{
								mprintf("      %2d.) %2d - %2d\n",z3+1,z+1,z2+1);
								z3++;
							}
						}
					}
					mprintf("\n");
				}

				if (g_bAdvanced2)
				{
					o->m_pCDF->m_bAxisDivide = AskYesNo("    Write out +/- correlation plot for this CDF (y/n)? [no] ",false);

					if (o->m_pCDF->m_bAxisDivide)
						o->m_pCDF->m_bAxisDivideAll = AskYesNo("    Also write tensor product and axis projection quotients (y/n)? [no] ",false);
				} else
				{
					o->m_pCDF->m_bAxisDivide = false;
					o->m_pCDF->m_bAxisDivideAll = false;
				}

				if (g_bAdvanced2)
				{
					o->m_pCDF->m_bGraceBunch = AskYesNo("    Write out grace stack (multiple 2D plots) for this CDF (y/n)? [no] ",false);
					if (o->m_pCDF->m_bGraceBunch)
					{
						o->m_pCDF->m_iGraceBunchC1 = AskUnsignedInteger("    How many graphs do you want do draw in the channel 1 grace stack (0=disable)? [10] ",10);
						o->m_pCDF->m_iGraceBunchC2 = AskUnsignedInteger("    How many graphs do you want do draw in the channel 2 grace stack (0=disable)? [10] ",10);
					}
				} else o->m_pCDF->m_bGraceBunch = false;
			}

			if (g_iCDFChannels == 3)
			{
				mprintf("\n");

				try { o->m_pCDF->m_pCombineList = new char[o->m_pCDF->m_iCombinationProd]; } catch(...) { o->m_pCDF->m_pCombineList = NULL; }
				if (o->m_pCDF->m_pCombineList == NULL) NewException((double)o->m_pCDF->m_iCombinationProd*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				if ((o->m_pCDF->m_iCombinations[0] > 1) || (o->m_pCDF->m_iCombinations[1] > 1) || (o->m_pCDF->m_iCombinations[2] > 1))
				{
					for (z=0;z<g_iCDFChannels;z++)
					{
						mprintf(WHITE,"CDF channel %d has the following %d observations:\n",z+1,o->m_pCDF->m_iCombinations[z]);
						o->ListCDFObservations(z);
						mprintf("\n");
					}
					if (!AskYesNo("    Combine each with each observation (y), or use only some combinations (n)? [yes] ",true))
					{
						for (z=0;z<o->m_pCDF->m_iCombinations[0];z++)
							for (z2=0;z2<o->m_pCDF->m_iCombinations[1];z2++)
								for (z3=0;z3<o->m_pCDF->m_iCombinations[2];z3++)
									o->m_pCDF->m_pCombineList[z*o->m_pCDF->m_iCombinations[1]*o->m_pCDF->m_iCombinations[2]+z2*o->m_pCDF->m_iCombinations[2]+z3] = 0;

						mprintf("\n    Please enter all the combinations you want to observe.\n");
						mprintf("    Enter each combination as a comma-separated %d-tuple (e.g. ",g_iCDFChannels);
						for (z=0;z<g_iCDFChannels;z++)
						{
							mprintf("%d",z+1);
							if (z+1 < g_iCDFChannels)
								mprintf(",");
						}
						mprintf(").\n\n");
_3nextcomb:
						AskString("    Enter combination (return=finished): ",&buf,"");
						if (strlen(buf)==0)
							goto _3combdone;
						tempwa.RemoveAll_KeepSize();
						ParseIntList(buf,&tempwa);
						if (tempwa.GetSize() != g_iCDFChannels)
						{
							eprintf("    Wrong input, %d instead of %d values.\n",tempwa.GetSize(),g_iCDFChannels);
							goto _3nextcomb;
						}
						for (z=0;z<g_iCDFChannels;z++)
						{
							if ((tempwa[z] < 1) || (tempwa[z] > o->m_pCDF->m_iCombinations[z]))
							{
								eprintf("    Wrong input, channel %d has only %d observations (%d requested).\n",z+1,o->m_pCDF->m_iCombinations[z],tempwa[z]);
								goto _3nextcomb;
							}
							tempwa[z]--;
						}
						if (o->m_pCDF->m_pCombineList[tempwa[0]*o->m_pCDF->m_iCombinations[1]*o->m_pCDF->m_iCombinations[2]+tempwa[1]*o->m_pCDF->m_iCombinations[2]+tempwa[2]] == 1)
						{
							eprintf("    This combination has already been added.\n");
							goto _3nextcomb;
						}
						o->m_pCDF->m_pCombineList[tempwa[0]*o->m_pCDF->m_iCombinations[1]*o->m_pCDF->m_iCombinations[2]+tempwa[1]*o->m_pCDF->m_iCombinations[2]+tempwa[2]] = 1;
						goto _3nextcomb;
					} else goto _3combineall;
				} else
				{
_3combineall:
					for (z=0;z<o->m_pCDF->m_iCombinationProd;z++)
						o->m_pCDF->m_pCombineList[z] = 1;
				}
_3combdone:
				o->m_pCDF->m_iCombinationsEnabled = 0;
				for (z=0;z<o->m_pCDF->m_iCombinationProd;z++)
					if (o->m_pCDF->m_pCombineList[z] != 0)
						o->m_pCDF->m_iCombinationsEnabled++;

				if (o->m_pCDF->m_iCombinationsEnabled == 1)
				{
					mprintf(WHITE,"\n    Using 1 combination for each RM-OM pair.\n\n");
				} else
				{
					mprintf(WHITE,"\n    Using %d combinations for each RM-OM pair:\n\n",o->m_pCDF->m_iCombinationsEnabled);
					z4 = 0;
					for (z=0;z<o->m_pCDF->m_iCombinations[0];z++)
						for (z2=0;z2<o->m_pCDF->m_iCombinations[1];z2++)
							for (z3=0;z3<o->m_pCDF->m_iCombinations[2];z3++)
								if (o->m_pCDF->m_pCombineList[z*o->m_pCDF->m_iCombinations[1]*o->m_pCDF->m_iCombinations[2]+z2*o->m_pCDF->m_iCombinations[2]+z3] == 1)
								{
									mprintf("      %2d.) %2d - %2d - %2d\n",z4+1,z+1,z2+1,z3+1);
									z4++;
								}
					mprintf("\n");
				}

				if (g_bAdvanced2)
				{
					o->m_pCDF->m_b3DSlices = AskYesNo("    Write out 2D slices for this CDF (y/n)? [no] ",false);

					if (o->m_pCDF->m_b3DSlices)
						for (z=0;z<3;z++)
							o->m_pCDF->m_i3DSliceIntervals[z] = AskUnsignedInteger("    How many slice intervals to create along channel %d axis (0=skip)? [0] ",0,z+1);
				} else o->m_pCDF->m_b3DSlices = false;
			} // end if channels == 3


			if (g_bAdvanced2)
			{
				o->m_pCDF->m_bDumpDat = AskYesNo("    Write out input tuples (very large!) for this CDF (y/n)? [no] ",false);
				if ((g_iCDFChannels == 2) && (g_iObsChannel[0] == 0) && (g_iObsChannel[1] == 1))
					o->m_pCDF->m_iNormalize = AskRangeInteger("    Normalize uniformly (3), data range (2), integral (1), or do not normalize (0)? [1] ",0,3,1);
						else o->m_pCDF->m_iNormalize = AskRangeInteger("    Normalize data range (2), integral (1), or do not normalize (0)? [1] ",0,2,1);
				if (o->m_pCDF->m_iNormalize == 1)
					o->m_pCDF->m_fNormValue = AskFloat("    Set CDF integral to which value? [1000000] ",1000000.0f);
				if (o->m_pCDF->m_iNormalize == 2)
					o->m_pCDF->m_fNormValue = AskFloat("    Set maximum entry to which value? [100] ",100.0f);
			} else
			{
				o->m_pCDF->m_bDumpDat = false;
				o->m_pCDF->m_iNormalize = 1;
				o->m_pCDF->m_fNormValue = 1000000.0f;
			}

			mprintf(BLUE,"\n<<< End of Combined Distribution Function <<<\n\n");
		} else
		{
			if (g_bRDF)
			{
				try { o->m_pRDF[0] = new CRDF(); } catch(...) { o->m_pRDF[0] = NULL; }
				if (o->m_pRDF[0] == NULL) NewException((double)sizeof(CRDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pRDF[0]->m_iShowMol = o->m_iShowMol;
				o->m_pRDF[0]->Parse();
				o->m_pRDF[0]->m_bSelf = o->m_bSelf;
			}
			if (g_bPlDF)
			{
				try { o->m_pPlDF[0] = new CPlDF(); } catch(...) { o->m_pPlDF[0] = NULL; }
				if (o->m_pPlDF[0] == NULL) NewException((double)sizeof(CPlDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pPlDF[0]->m_iShowMol = o->m_iShowMol;
				o->m_pPlDF[0]->Parse();
				o->m_pPlDF[0]->m_bSelf = o->m_bSelf;
			}
			if (g_bLiDF)
			{
				try { o->m_pLiDF[0] = new CLiDF(); } catch(...) { o->m_pLiDF[0] = NULL; }
				if (o->m_pLiDF[0] == NULL) NewException((double)sizeof(CLiDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pLiDF[0]->m_iShowMol = o->m_iShowMol;
				o->m_pLiDF[0]->Parse();
				o->m_pLiDF[0]->m_bSelf = o->m_bSelf;
			}
			if (g_bVHDF)
			{
				try { o->m_pVHDF = new CVHDF(); } catch(...) { o->m_pVHDF = NULL; }
				if (o->m_pVHDF == NULL) NewException((double)sizeof(CVHDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pVHDF->m_iShowMol = o->m_iShowMol;
				o->m_pVHDF->m_bSelf = o->m_bSelf;
				o->m_pVHDF->Parse();
			}
			if (g_bADF)
			{
				try { o->m_pADF[0] = new CADF(); } catch(...) { o->m_pADF[0] = NULL; }
				if (o->m_pADF[0] == NULL) NewException((double)sizeof(CADF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pADF[0]->m_iShowMol = o->m_iShowMol;
				o->m_pADF[0]->m_bSelf = o->m_bSelf;
				o->m_pADF[0]->Parse();
			}
			if (g_bDDF)
			{
				try { o->m_pDDF[0] = new CDDF(); } catch(...) { o->m_pDDF[0] = NULL; }
				if (o->m_pDDF[0] == NULL) NewException((double)sizeof(CDDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pDDF[0]->m_iShowMol = o->m_iShowMol;
				o->m_pDDF[0]->m_bSelf = o->m_bSelf;
				o->m_pDDF[0]->Parse();
			}
			if (g_bDipDF)
			{
				try { o->m_pDipDF[0] = new CDipDF(); } catch(...) { o->m_pDipDF[0] = NULL; }
				if (o->m_pDipDF[0] == NULL) NewException((double)sizeof(CDipDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pDipDF[0]->m_iShowMol = o->m_iShowMol;
				o->m_pDipDF[0]->m_bSelf = o->m_bSelf;
				o->m_pDipDF[0]->Parse();
			}
			if (g_bVDF)
			{
				try { o->m_pVDF[0] = new CVDF(); } catch(...) { o->m_pVDF[0] = NULL; }
				if (o->m_pVDF[0] == NULL) NewException((double)sizeof(CVDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pVDF[0]->m_iShowMol = o->m_iShowMol;
				o->m_pVDF[0]->m_bSelf = true;
				o->m_pVDF[0]->Parse();
				g_bUseVelocities = true;
			}
			if (g_bMSD)
			{
				try { o->m_pMSD = new CMSD(); } catch(...) { o->m_pMSD = NULL; }
				if (o->m_pMSD == NULL) NewException((double)sizeof(CMSD),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				 
				o->m_pMSD->m_iShowMol = o->m_iShowMol;
				o->m_pMSD->Parse();
			}
		}

		if ((g_bCDF && (g_iCDFChannels == 2)) || g_bRDF || g_bADF || g_bVDF || g_bDDF || g_bDipDF)
		{
			if (!o->m_bSecondShowMol)
			{
				if (AskYesNo("    Save temporal development for this observation (y/n)? [no] ",false))
				{
					mprintf(WHITE,"\n>>> Save temporal development >>>\n\n");
					o->m_bTimeDev = true;

					if (g_fTimestepLength == 0)
					{
						g_fTimestepLength = AskFloat("    Enter the length of one trajectory time step in fs: [0.5] ",0.5f);
						mprintf("\n");
					}

					if (g_bDDF && (!g_bCDF))
					{
						if ((o->m_pDDF[0]->m_iDeriv == 0) && (!o->m_pDDF[0]->m_bCosine) && (!o->m_pDDF[0]->m_bAbs) && (!o->m_pDDF[0]->m_bPositive))
						{
							mprintf("\n    You can include full rotations into the DDF's temporal development.\n");
							mprintf("    This means, if a dihedral makes a full rotation, the value will be 360 deg insteat of 0 deg.\n");
							mprintf("    You may count the number of full rotations during simulation in this way.\n\n");
							o->m_pDDF[0]->m_bRotate = AskYesNo("    Include full rotations into the DDF's temporal development (y/n)? [no] ",false);
						} else o->m_pDDF[0]->m_bRotate = false;
						if (o->m_pDDF[0]->m_bRotate)
						{
							mprintf(WHITE,"\n   Warning: ");
							mprintf("This only works if your trajectory time step is small enough to ensure\n");
							mprintf("             that the dihedral is never changing more than 180 deg within one step!\n\n");
						}
					}
					if (g_bRDF || g_bADF || g_bVDF || g_bDDF || g_bDipDF)
					{
						if (AskYesNo("    Create combined development/histogram-plots for 2D analyses (y/n)? [yes] ",true))
						{
							g_bCombined = true;
							o->m_bCombinedPlot = true;
							if (AskYesNo("    Use grey tones for combined plot (y) or standard colors (n)? [no] ",false))
							{
								o->m_bCombinedGreyMode = true;
								o->m_iCombinedGreyMin = AskRangeInteger("    Darkest shade of grey to use (0=black, 255=white)? [128] ",0,255,128);
								o->m_iCombinedGreyMax = AskRangeInteger("    Lightest shade of grey to use (0=black, 255=white)? [224] ",o->m_iCombinedGreyMin,255,224);
								o->m_iCombinedGreyShades = AskUnsignedInteger("    How many different shades of grey to use? [4] ",4);
							} else o->m_bCombinedGreyMode = false;
						} else o->m_bCombinedPlot = false;
					} else o->m_bCombinedPlot = false;

					if (!g_bVDF && ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() > 1)
					{
						mprintf("    Save development for which rep. of the reference molecule (e.g. 1,3-7)? [1] ");
						inpprintf("! Save development for which rep. of the reference molecule (e.g. 1,3-7)? [1]\n");
						myget(&buf);
						if (strlen(buf)==0)
							o->m_waSaveRefList.Add(0);
						else {
							ParseIntList(buf,&o->m_waSaveRefList);
							for (z=0;z<o->m_waSaveRefList.GetSize();z++)
								o->m_waSaveRefList[z]--;
						}
					} else o->m_waSaveRefList.Add(0);
					if (o->m_iShowMol != -1)
					{
						if (((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize() > 1)
						{
							mprintf("    Save development for which rep. of the observed molecule (e.g. 1,3-7)? [all] ");
							inpprintf("! Save development for which rep. of the observed molecule (e.g. 1,3-7)? [all]\n");
							myget(&buf);
							if (strlen(buf)==0)
							{
								for (z=0;z<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z++)
									o->m_waSaveShowList.Add(z);
							} else
							{
								ParseIntList(buf,&o->m_waSaveShowList);
								for (z=0;z<o->m_waSaveShowList.GetSize();z++)
									o->m_waSaveShowList[z]--;
							}
						} else o->m_waSaveShowList.Add(0);
					} else o->m_waSaveShowList.Add(0); // Dummy
					if ((o->m_waSaveRefList.GetSize() > 1) && (o->m_bOthers))
						o->m_bSaveSeparateFiles = AskYesNo("    Developments for different ref. molecules to same (n) or different (y) files? [yes] ",true);
							else o->m_bSaveSeparateFiles = false;
					if ((g_bCDF) && (g_iCDFChannels == 2))
					{
						if (AskYesNo("    Generate time development animation (TDO) for CDF (y/n)? [no] ",false))
						{
							o->m_pCDF->m_bTDAnimation = true;
							o->m_pCDF->m_iTDASteps = AskUnsignedInteger("    How many steps should the animation have? [1000] ",1000);
							o->m_pCDF->m_iTDAStride = AskUnsignedInteger("    Enter distance between the frames in timesteps: [100] ",50);
							o->m_pCDF->m_iTDATail = AskUnsignedInteger("    Enter length of the tail in timesteps: [100] ",100);
							o->m_pCDF->m_bTDATrace = AskYesNo("    Show trace (y/n)? [yes] ",true);
							o->m_pCDF->m_iTDAResX = AskUnsignedInteger("    Enter width (in pixel) of the TDA images: [640] ",640);
							o->m_pCDF->m_iTDAResY = AskUnsignedInteger("    Enter height (in pixel) of the TDA images: [480] ",480);
						} else o->m_pCDF->m_bTDAnimation = false;
					}
					mprintf("\n    Saving temporal development for reference molecules ");
					for (z=0;z<o->m_waSaveRefList.GetSize();z++)
					{
						if (z < (int)o->m_waSaveRefList.GetSize()-1)
						{
							mprintf("%d, ",o->m_waSaveRefList[z]+1);
							if (((z+1) % 16) == 0)
								mprintf("\n      ");
						} else mprintf("%d",o->m_waSaveRefList[z]+1);
					}
					if (o->m_bOthers)
					{
						mprintf("\n    and for observed molecules ");
						for (z=0;z<o->m_waSaveShowList.GetSize();z++)
						{
							if (z < (int)o->m_waSaveShowList.GetSize()-1)
							{
								mprintf("%d, ",o->m_waSaveShowList[z]+1);
								if (((z+1) % 16) == 0)
									mprintf("\n      ");
							} else mprintf("%d",o->m_waSaveShowList[z]+1);
						}
					}
					mprintf(".\n");
					mprintf(WHITE,"\n<<< End of Save temporal development <<<\n\n");
				} else o->m_bTimeDev = false;

				if (g_bAdvanced2)
				{
					if (AskYesNo("    Create a temporal difference plot for this observation (y/n)? [no] ",false))
					{
						o->m_bTimeDiff = true;
						g_bTimeDiff = true;
						o->m_iTimeDiffDepth = AskUnsignedInteger("    Enter temporal depth of this plot in time steps: [1000] ",1000);
						o->m_b3DTimeDiff = AskYesNo("    Create also 3D temporal difference plots (y/n)? [yes] ",true);
						if (o->m_b3DTimeDiff)
						{
							tf = 0;
							tf2 = 0;
							if (g_bRDF)
							{
								tf = o->m_pRDF[0]->m_fMinDist;
								tf2 = o->m_pRDF[0]->m_fMaxDist;
							}
							if (g_bADF)
							{
								tf = o->m_pADF[0]->m_fMinAngle;
								tf2 = o->m_pADF[0]->m_fMaxAngle;
							}
							if (g_bDDF)
							{
								tf = o->m_pDDF[0]->m_fMinAngle;
								tf2 = o->m_pDDF[0]->m_fMaxAngle;
							}
							if (g_bVDF)
							{
								tf = o->m_pVDF[0]->m_fMinSpeed;
								tf2 = o->m_pVDF[0]->m_fMaxSpeed;
							}
							if (g_bDipDF)
							{
								tf = o->m_pDipDF[0]->m_fDipoleMin;
								tf2 = o->m_pDipDF[0]->m_fDipoleMax;
							}
							o->m_iTimeDiffStride3D = AskUnsignedInteger("    Take every n-th time step for the tau axis: [%d] ",o->m_iTimeDiffDepth/100,o->m_iTimeDiffDepth/100);
							o->m_fTimeDiffMinVal3D = AskFloat("    Enter min. value for the starting point: [%.2f] ",tf,tf);
							o->m_fTimeDiffMaxVal3D = AskFloat("    Enter max. value for the starting point: [%.2f] ",tf2,tf2);
							o->m_iTimeDiffRes3D = AskUnsignedInteger("    Enter binning resolution for the starting value: [100] ",100);
							mprintf("\n    Temporal before/after 3D plot:\n");
							o->m_iTimeDiffDistSteps= AskUnsignedInteger("    Every how many time steps create a 2D slice (0 to disable): [0] ",0);
							o->m_fTimeDiffDistMinValX = AskFloat("    Enter min. value for the starting point: [%.2f] ",tf,tf);
							o->m_fTimeDiffDistMaxValX = AskFloat("    Enter max. value for the starting point: [%.2f] ",tf2,tf2);
							o->m_iTimeDiffDistResX = AskUnsignedInteger("    Enter binning resolution for the X axis: [100] ",100);
							o->m_fTimeDiffDistMinValY = AskFloat("    Enter min. value for the end point: [%.2f] ",tf,tf);
							o->m_fTimeDiffDistMaxValY = AskFloat("    Enter max. value for the end point: [%.2f] ",tf2,tf2);
							o->m_iTimeDiffDistResY = AskUnsignedInteger("    Enter binning resolution for the Y axis: [100] ",100);
						}
						mprintf("\n");
					} else o->m_bTimeDiff = false;
				} else o->m_bTimeDiff = false;
			} else
			{
				o->m_bTimeDev = false;
				o->m_bTimeDiff = false;
			}
		} else
		{
			o->m_bTimeDev = false;
			o->m_bTimeDiff = false;
		}

		if ((g_iFixMol != -1) && (o->m_iShowMol != -1) && (g_bRDyn || g_bDipDF || g_bMSD || g_bADF || g_bPlDF || g_bLiDF || g_bDDF || g_bVACF || g_bDipACF ||g_bRevSDF || g_bSDF || g_bPlProj || g_bVHDF || g_bRDF || g_bVDF || /*g_bFDF ||*/ g_bCond))
		{
			if (g_bCond)
				goto _askcond;
			if (AskYesNo("    Add a condition to this observation (y/n)? [no] ",false))
			{
_askcond:
				mprintf(GREEN,"\n>>> Condition input >>>\n\n");
//				mprintf("You may enter several sets of conditions. They are connected with \"or\",\nonly one of them needs to apply to take a configuration into account.\n\n");
//				mprintf("In each set of conditions, you may enter several conditions. They are connected\nwith \"and\", all of them have to apply to make this set of conditions apply.\n");

				if (o->m_bSecondShowMol)
				{
					if (!AskYesNo("    Add a condition between RM and 1st OM (y/n)? [yes] ",true))
						goto _no1stcond;
					mprintf(WHITE,"\n    *** Input of Condition between RM and 1st OM ***\n\n");
				}

				try { o->m_pConditions = new CConditionGroup(); } catch(...) { o->m_pConditions = NULL; }
				if (o->m_pConditions == NULL) NewException((double)sizeof(CConditionGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				o->m_pConditions->m_iShowMol = o->m_iShowMol;
				o->m_pConditions->Parse(g_iFixMol,o->m_iShowMol);
_no1stcond:

				if (o->m_bSecondShowMol)
				{
					mprintf("\n");
					if (!AskYesNo("    Add a condition between RM and 2nd OM (y/n)? [yes] ",true))
						goto _no2ndcond;
					mprintf(WHITE,"\n    *** Input of Condition between RM and 2nd OM ***\n\n");

					try { o->m_pConditionsOM2 = new CConditionGroup(); } catch(...) { o->m_pConditionsOM2 = NULL; }
					if (o->m_pConditionsOM2 == NULL) NewException((double)sizeof(CConditionGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					o->m_pConditionsOM2->m_iShowMol = o->m_iShowMol2;
					o->m_pConditionsOM2->Parse(g_iFixMol,o->m_iShowMol2);
_no2ndcond:;
				}

				if (o->m_pConditions != NULL)
				{
					if (o->m_pConditions->m_bInvertCondition)
					{
						o->m_bBinOnlyNotPassedAtoms = false;
						o->m_bBinOnlyPassedAtoms = false;
						goto _nobinonly;
					}
				}

				if (o->m_pConditionsOM2 != NULL)
				{
					if (o->m_pConditionsOM2->m_bInvertCondition)
					{
						o->m_bBinOnlyNotPassedAtoms = false;
						o->m_bBinOnlyPassedAtoms = false;
						goto _nobinonly;
					}
				}

				mprintf("\n    Normally, TRAVIS evaluates the condition for each RM-OM pair. If the condition\n");
				mprintf("    is fulfilled, all atoms from this pair are evaluated. The following question\n");
				mprintf("    enables to take into account only exactly the atoms which fulfilled the condition.\n\n");

				mprintf("    If you have a distance condition with \"nearest neighbor mode\", only the atom\n");
				mprintf("    from the OM that is closest to the RM will be taken into account.\n\n");

				o->m_bBinOnlyPassedAtoms = AskYesNo("    Add only atoms to the bin that passed the condition(s) (y/n)? [no] ",false);
/*				if (!o->m_bBinOnlyPassedAtoms)
					o->m_bBinOnlyNotPassedAtoms = AskYesNo("    Add only atoms to the bin that did NOT pass the condition(s) (y/n)? [no] ",false);
						else*/ o->m_bBinOnlyNotPassedAtoms = false;
_nobinonly:

				g_bSaveCondSnapshot = AskYesNo("    Save a snapshot every time the conditions are fulfilled (y/n)? [no] ",false);
				
				if (g_bSaveCondSnapshot)
					g_bNeedMoleculeWrap = true;

				if (g_bSaveCondSnapshot)
				{
					g_bSaveCondWholeBox = AskYesNo("    Save the whole box (y) or only the RM/OM pair (n)? [no] ",false);
					if (g_bSaveCondWholeBox)
						mprintf("\n    The RM that fulfills the condition will be centered (coordinates 0|0|0).\n");
				}

				mprintf(GREEN,"\n<<< End of Condition input <<<\n\n");
			}
		}

		if (g_bVACF)
		{
			mprintf(WHITE,"\n>>> Velocity autocorrelation function >>>\n\n");

			try { o->m_pVACF = new CACF(); } catch(...) { o->m_pVACF = NULL; }
			if (o->m_pVACF == NULL) NewException((double)sizeof(CACF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pVACF->m_iShowMol = o->m_iShowMol;
			o->m_pVACF->Parse();
			mprintf(WHITE,"\n<< End of Velocity autocorrelation function <<<\n\n");
		}

/*		if (g_bDipACF)
		{
			mprintf(WHITE,"\n>>> Dipole moment autocorrelation function >>>\n\n");

			try { o->m_pDipACF = new CACF(); } catch(...) { o->m_pDipACF = NULL; }
			if (o->m_pDipACF == NULL) NewException((double)sizeof(CACF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pDipACF->Parse();
			mprintf(WHITE,"\n<<< End of Dipole moment autocorrelation function <<<\n\n");
		}*/

		if (g_bUseVelocities)
		{
			if (g_iFixMol != -1)
				o->m_bVelocityRelToRef = (AskRangeInteger("    Use absolute velocities (0) or relative to the reference molecule (1)? [0] ",0,1,0) != 0);
					else o->m_bVelocityRelToRef = false;
		}

		mprintf(YELLOW,"\n<<< End of Observation %d <<<\n",g_oaObserv.GetSize());

		if (AskYesNo("\n    Add another observation (y/n)? [no] ",false))
			goto _nextsdf;
		mprintf("\n");

		mprintf(WHITE,">>> Observation List >>>\n\n");
		for (z2=0;z2<g_oaObserv.GetSize();z2++)
		{
			mprintf(YELLOW,"  *** Observation %d\n",z2+1);
			o = (CObservation*)g_oaObserv[z2];

			if (g_bCond)
				mprintf("    - Condition between %s and %s\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName);

			if (g_bAggregation)
			{
				mprintf("    - Aggregation functions - %d value sets:\n",o->m_pDACF->m_oaSubDACFs.GetSize());
				for (z3=0;z3<o->m_pDACF->m_oaSubDACFs.GetSize();z3++)
					mprintf("      > %s\n",((CDACFSub*)o->m_pDACF->m_oaSubDACFs[z3])->m_sName);
			}

			if (g_bNbExchange)
			{
				mprintf("    - Neighborhood exchange functions:\n");
				for (z3=0;z3<o->m_pDACF->m_oaSubDACFs.GetSize();z3++)
					mprintf("      > %s\n",((CDACFSub*)o->m_pDACF->m_oaSubDACFs[z3])->m_sName);
			}

			if (g_bRDyn)
				mprintf("    - Reorientation Dynamics: %s\n",o->m_pRDyn->m_sName);

			if (g_bIRSpec)
				mprintf("    - IR Spectrum: %s\n",o->m_pIRSpec->m_sName);

			if (g_bDens)
				mprintf("    - Density DF: %s\n",o->m_pDensityDF->m_sName);

			if (g_bSDF)
				mprintf("    - SDF: %s\n",o->m_pSDF->m_sName);

			if (g_bPlProj)
				mprintf("    - Plane Projection DF: %s\n",o->m_pPlProj->m_sName);

			if (g_bRevSDF)
				mprintf("    - Pseudo SDF: %s\n",o->m_pRevSDF->m_sName);

			if (g_bNbAnalysis)
				mprintf("    - Neighborhood Analysis: %s\n",o->m_pNbAnalysis->m_sName);

			if (g_bVACF)
				mprintf("    - VACF: %s\n",o->m_pVACF->m_sName);

			if (g_bMSD)
				mprintf("    - MSD: %s\n",o->m_pMSD->m_sName);

			if (g_bVHDF)
				mprintf("    - VHDF: %s\n",o->m_pVHDF->m_sName);

			if (g_bCDF)
			{
				mprintf("    - Combined Distribution Function:\n");
				for (z3=0;z3<g_iCDFChannels;z3++)
				{
					if (o->m_pRDF[z3] != NULL)
						mprintf("      - Channel %d: RDF: %s\n",z3+1,o->m_pRDF[z3]->m_sName);
					if (o->m_pADF[z3] != NULL)
						mprintf("      - Channel %d: ADF: %s\n",z3+1,o->m_pADF[z3]->m_sName);
					if (o->m_pDDF[z3] != NULL)
						mprintf("      - Channel %d: DDF: %s\n",z3+1,o->m_pDDF[z3]->m_sName);
					if (o->m_pDipDF[z3] != NULL)
						mprintf("      - Channel %d: DipDF: %s\n",z3+1,o->m_pDipDF[z3]->m_sName);
					if (o->m_pVDF[z3] != NULL)
						mprintf("      - Channel %d: VDF: %s\n",z3+1,o->m_pVDF[z3]->m_sName);
					if (o->m_pPlDF[z3] != NULL)
						mprintf("      - Channel %d: PlDF: %s\n",z3+1,o->m_pPlDF[z3]->m_sName);
					if (o->m_pLiDF[z3] != NULL)
						mprintf("      - Channel %d: LiDF: %s\n",z3+1,o->m_pLiDF[z3]->m_sName);
				}
			} else
			{
				if (g_bRDF)
					mprintf("    - RDF: %s\n",o->m_pRDF[0]->m_sName);
				if (g_bDipDF)
					mprintf("    - DipDF: %s\n",o->m_pDipDF[0]->m_sName);
				if (g_bADF)
					mprintf("    - ADF: %s\n",o->m_pADF[0]->m_sName);
				if (g_bDDF)
					mprintf("    - DDF: %s\n",o->m_pDDF[0]->m_sName);
				if (g_bVDF)
					mprintf("    - VDF: %s\n",o->m_pVDF[0]->m_sName);
				if (g_bPlDF)
					mprintf("    - PlDF: %s\n",o->m_pPlDF[0]->m_sName);
				if (g_bLiDF)
					mprintf("    - LiDF: %s\n",o->m_pLiDF[0]->m_sName);
			}
		}
		mprintf(WHITE,"\n<<< End of Observation List <<<\n\n");
	}
_endobs:

	if (!Interface_BeforeAnalysis2())
		return false;
	
	if (g_bVACF && g_bGlobalVACF)
	{
		try { g_pGlobalVACF = new CACF(); } catch(...) { g_pGlobalVACF = NULL; }
		if (g_pGlobalVACF == NULL) NewException((double)sizeof(CACF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		mprintf(WHITE,"\n>>> Global velocity autocorrelation function >>>\n\n");
		g_pGlobalVACF->Parse();
		mprintf(WHITE,"\n<<< End of global velocity autocorrelation function <<<\n\n");
	}

	if (g_bIRSpec && g_bGlobalIR)
	{
		try { g_pGlobalIR = new CReorDyn(); } catch(...) { g_pGlobalIR = NULL; }
		if (g_pGlobalIR == NULL) NewException((double)sizeof(CReorDyn),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		mprintf(WHITE,"\n*** Definition of system-wide IR spectrum ***\n\n");
		g_pGlobalIR->ParseSpec();
	}

	if (g_bSDF)
	{
		g_bSDFMap = AskYesNo("    Map a quantity to the surface of the SDF (y/n)? [no] ",false);

		if (g_bSDFMap)
		{
_sdfmapnext:
			CSDFMap *sma;

			sma = new CSDFMap();
			g_oaSDFMaps.Add(sma);

			sma->Parse();

			if (AskYesNo("    Map another quantity to the surface of the SDF (y/n)? [no] ",false))
				goto _sdfmapnext;

			mprintf("\n");
		}
	}

	if (g_bSDF || 
      g_bCreateRevSDF)
	{
		g_iSDFSmoothGrade = AskInteger("    Up to which degree should the SDFs be smoothened (0=not at all)? [2] ",2);

		if (g_bSDFMap)
			g_iSDFMapSmoothGrade = AskInteger("    Up to which degree should the SDF surface maps be smoothened (0=not at all)? [2] ",2);

		if (g_bSDF)
		{
_sdf_norm_again:
			g_bSDFUniform = AskYesNo("    Calculate SDF values in nm^-3 (n) or relative to uniform density (y)? [no] ",false);
			if (g_bSDFUniform && ((!g_bPeriodicX) || (!g_bPeriodicY) || (!g_bPeriodicZ)))
			{
				eprintf("\n    Error: Uniform particle density is only defined for XYZ-periodic systems.\n\n");
				goto _sdf_norm_again;
			}
		}
//			g_iSDFScale = AskRangeInteger("    SDF values in ppm (1), pm^-3 (2), nm^-3 (3) or rel. to uniform density (4)? [4] ",1,4,4) - 1;
//			g_iSDFScale = 2; // nm^-3
	} else if (g_bCDF && (g_iCDFChannels == 3))
		g_iSDFSmoothGrade = AskInteger("    Up to which degree should the 3D CDFs be smoothened (0=not at all)? [2] ",2);

//	if (g_bVDF || g_bSVDF)
//		g_bSaveVelForce = AskYesNo("    Save maximum/averaged velocity time development (y/n)? [yes] ",true);

//	if (g_bVDF || g_bFDF)
//	{
//		g_fVelPercentage = AskFloat("    Which percentage of all velocities/forces to take into account? [95] ",95.0f);
//		g_fForcePercentage = g_fVelPercentage;
//	}
	g_bSaveVelForce = false;

	if ((g_bSDF && g_bAdvanced2) || g_bSaveRefEnv || g_bCutCluster || g_bSaveJustTraj)
	{
		mprintf(WHITE,"\n>>> Coordinate output options >>>\n\n");

		if ((g_bSDF && g_bAdvanced2) || g_bSaveRefEnv || g_bCutCluster)
		{
			if (g_bAdvanced2)
				g_bWriteAtomwise = (AskRangeInteger("    Sort output coordinates by molecules (0) or by element types (1)? [molecules] ",0,1,0)!=0);
					else g_bWriteAtomwise = false;
		}

		g_bSaveVirtAtoms = AskYesNo("    Save also virtual atoms (y/n)? [no] ",false);

		if (g_bSaveVirtAtoms)
		{
			if (AskYesNo("    Ref.Env.: Use alias names (instead of #) for virtual atoms (y/n)? [no] ",false))
			{
				for (z2=0;z2<g_oaMolecules.GetSize();z2++)
				{
					CMolecule *m = (CMolecule*)g_oaMolecules[z2];
					for (z=0;z<m->m_baAtomIndex.GetSize();z++)
					{
						if (m->m_baAtomIndex[z] != g_iVirtAtomType)
							continue;
						for (z3=0;z3<m->m_waAtomCount[z];z3++)
						{
							mprintf("    How to name atom #%d in %s? [#%d] ",z3+1,m->m_sName,z3+1);
							inpprintf("! How to name atom #%d in %s? [#%d]\n",z3+1,m->m_sName,z3+1);
							myget(&buf);
							if (strlen(buf) != 0)
								strcpy(((CVirtualAtom*)g_oaVirtualAtoms[m->m_laVirtualAtoms[z3]])->m_sLabel,buf);
						}
					}
				}
			}
		}

		if (g_bSaveRefEnv || g_bCutCluster)
			g_bEnvWriteDetailedInfo = AskYesNo("    Write detailed (long) comment lines into output xyz file (y/n)? [yes] ",true);
				else g_bEnvWriteDetailedInfo = false;

		if ((g_bSaveRefEnv || g_bCutCluster) && !g_bEnvDisableSortNb)
			g_bEnvSortNb = AskYesNo("    Sort molecules according to distance from reference molecule (y/n)? [yes] ",true);
				else g_bEnvSortNb = false;

		mprintf(WHITE,"\n<<< End of Coordinate output options <<<\n\n");
	} else
	{
		g_bWriteAtomwise = false;
		g_bSaveVirtAtoms = false;
	}

	if (g_bSaveRefEnv || g_bCutCluster)
		g_bCenterZero = AskYesNo("    Put the center of the system to (0|0|0) (y) or to (x/2|y/2|z/2) (n)? [yes] ",true);

	g_bMiddleAvg = false;
	g_iSwapAtoms = 0;

	if (g_bSaveJustTraj)
	{
		mprintf(WHITE,"\n>>> Process Trajectory >>>\n\n");

		if (g_bAdvanced2)
			g_bSaveTrajNoRot = AskYesNo("    Remove angular momentum from trajectory (y/n)? [no] ",false);
				else g_bSaveTrajNoRot = false;

		if (g_bSaveTrajNoRot)
		{
			mprintf("\n    This centers the mass center of the system. All atoms of the system are saved.\n\n");
			g_bSaveCoordsUnchanged = false;
			g_bSaveJustCenter = false;
			if(AskYesNo("    Put the center of the system to (0|0|0) (y) or to (x/2|y/2|z/2) (n)? [yes] ",true))
				g_bCenterZero = true;
					else g_bCenterZero = false;
			goto _norot;
		} else
		{
			switch(AskRangeInteger("    Put the center of the system to (0|0|0) (0), to (x/2|y/2|z/2) (1), or leave coords unchanged (2)? [1] ",0,2,1))
			{
				case 0:
					g_bCenterZero = true;
					g_bSaveCoordsUnchanged = false;
					break;
				case 1:
					g_bCenterZero = false;
					g_bSaveCoordsUnchanged = false;
					break;
				case 2:
					g_bSaveCoordsUnchanged = true;
					break;
			}

			if (!g_bSaveCoordsUnchanged)
				g_bSaveJustCenter = AskYesNo("\n    Put a specific atom into the box center (y/n)? [no] ",false);
					else g_bSaveJustCenter = false;
		}

		if (g_bSaveJustCenter)
		{
			if (g_oaMolecules.GetSize() > 1)
			{
//				sprintf(buf,"    Choose center atom from which of the molecules (");
				buf.sprintf("    Choose center atom from which of the molecules (");
				for (z=0;z<g_oaMolecules.GetSize();z++)
				{
//					sprintf(buf2,"%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
//					strcat(buf,buf2);
//					if (z < g_oaMolecules.GetSize()-1)
//						strcat(buf,", ");
					buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
					buf.strcat(buf2);
					if (z < g_oaMolecules.GetSize()-1)
						buf.strcat(", ");
				}
//				strcat(buf,")? ");
				buf.strcat(")? ");
				g_iSaveJustMol = AskRangeInteger_ND(buf,1,g_oaMolecules.GetSize()) - 1;
			} else
			{
				g_iSaveJustMol = 0;
				mprintf("    Choosing center atom from %s.\n",((CMolecule*)g_oaMolecules[0])->m_sName);
			}
			if (((CMolecule*)g_oaMolecules[g_iSaveJustMol])->m_laSingleMolIndex.GetSize() > 1)
				g_iSaveJustSM = AskRangeInteger("    Take which representant from %s (1-%d)? [1] ",1,((CMolecule*)g_oaMolecules[g_iSaveJustMol])->m_laSingleMolIndex.GetSize(),1,((CMolecule*)g_oaMolecules[g_iSaveJustMol])->m_sName,((CMolecule*)g_oaMolecules[g_iSaveJustMol])->m_laSingleMolIndex.GetSize()) - 1;
					else g_iSaveJustSM = 0;
_proccenter:
			AskString("    Which atom from %s %d to put into the box center? [#2] ",&buf,"#2",((CMolecule*)g_oaMolecules[g_iSaveJustMol])->m_sName,g_iSaveJustSM+1);
			if (!ParseAtom(buf,g_iSaveJustMol,g_iSaveJustAtomType,g_iSaveJustRealAtomType,g_iSaveJustAtom))
				goto _proccenter;
		}
		g_iSaveGesAtoms = 0;
		if (!AskYesNo("    Save all atoms in the system (y/n)? [yes] ",true))
		{
			mprintf("\n");
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				if (AskYesNo("    Save (some/all) atoms from molecule %s (y/n)? [yes] ",true,((CMolecule*)g_oaMolecules[z])->m_sName))
				{
					try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
					if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					g_oaSaveMolecules.Add(ag);
_savejust1:			mprintf("    Which atoms to save from %s (e.g. \"C1,C3-5,H\")? [all] ",((CMolecule*)g_oaMolecules[z])->m_sName);
					inpprintf("! Which atoms to save from %s (e.g. \"C1,C3-5,H\")? [all]\n",((CMolecule*)g_oaMolecules[z])->m_sName);
					myget(&buf);
					if (strlen(buf)==0)
						ag->AddAllAtoms((CMolecule*)g_oaMolecules[z],g_bSaveVirtAtoms);
					else if (!ag->ParseAtoms((CMolecule*)g_oaMolecules[z],buf))
						goto _savejust1;
					g_iSaveGesAtoms += ag->m_iAtomGes * ((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();
				}
			}
		} else
		{
_norot:
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
				if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				g_oaSaveMolecules.Add(ag);
				ag->AddAllAtoms((CMolecule*)g_oaMolecules[z],g_bSaveVirtAtoms);
				g_iSaveGesAtoms += ag->m_iAtomGes * ((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();
			}
		}

		if (g_bUnknownElements)
		{
			g_bUnwrapWannier = AskYesNo("    Write out Wannier centers together with trajectory (y/n)? [no] ",false);

			if (g_bUnwrapWannier)
				ParseDipole();

		} else g_bUnwrapWannier = false;

		if ((!g_bSaveCoordsUnchanged) && (!g_bSaveTrajNoRot))
			g_bUnwrap = AskYesNo("    Try to unwrap the trajectory (y/n)? [no] ",false);
				else g_bUnwrap = false;

		if (g_bBoxNonOrtho)
		{
			mprintf("\n");
			g_bWriteOrtho = AskYesNo("    Write \"orthonormalized\" coordinates (y/n)? [no] ",false);
			if (g_bWriteOrtho)
			{
				g_fWriteOrthoFac = AskFloat("      Enter scaling factor (x,y,z coordinates will range from 0 to this value) in pm: [%.4f] ",pow(g_fBoxVolume,1.0/3.0)*100.0,pow(g_fBoxVolume,1.0/3.0)*100.0);
				mprintf("\n");
			}
		}

		if (g_bUnwrap)
			mprintf("\n    Attention! This feature only works if the time step\n    in the trajectory is quite small. One particle may never\n    move further than half of the box length in one time step!\n");

		mprintf("\nSaving %d atoms per step:\n",g_iSaveGesAtoms);
		for (z2=0;z2<g_oaSaveMolecules.GetSize();z2++)
		{
			ag = (CAtomGroup*)g_oaSaveMolecules[z2];
			mprintf("  - In %s: ",ag->m_pMolecule->m_sName);
			for (z=0;z<ag->m_baAtomType.GetSize();z++)
			{
				for (z3=0;z3<((CxIntArray*)ag->m_oaAtoms[z])->GetSize();z3++)
				{
					mprintf("%s%d",((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z]])->m_sName,((CxIntArray*)ag->m_oaAtoms[z])->GetAt(z3)+1);
					if (z3 < ((CxIntArray*)ag->m_oaAtoms[z])->GetSize()-1)
						mprintf(", ");
				}
				if (z+1 < ag->m_baAtomType.GetSize())
					mprintf(", ");
			}
			mprintf("\n");
		}
//		g_bTrajAtomwise = (AskRangeInteger("\n    Save trajectory ordered molecule-wise (0) or atom-wise (1)? [0] ",0,1,0)!=0);

		mprintf("\n");

		if (g_bAdvanced2)
		{
			if ((int)g_iSaveGesAtoms == g_iGesAtomCount)
			{
				switch(AskRangeInteger("    Sort output coordinates by molecules (0), element types (1), or input file order (2)? [0] ",0,2,0))
				{
					case 0:
						g_bWriteAtomwise = false;
						g_bWriteInputOrder = false;
						break;
					case 1:
						g_bWriteAtomwise = true;
						g_bWriteInputOrder = false;
						break;
					case 2:
						g_bWriteAtomwise = false;
						g_bWriteInputOrder = true;
						break;
					default:
						eprintf("Should never happen ^^\n");
						abort();
				}
			} else 
			{
				g_bWriteAtomwise = (AskRangeInteger("    Sort output coordinates by molecules (0) or by element types (1)? [0] ",0,1,0)!=0);
				g_bWriteInputOrder = false;
			}
		} else
		{
			g_bWriteAtomwise = false;
			g_bWriteInputOrder = false;
		}

		mprintf(WHITE,"\n<<<< End of Process Trajectory <<<\n\n");
	} else g_bUnwrap = false;

	if (g_bUnwrap)
	{
		g_vaUnwrapArray.SetSize(g_oaSingleMolecules.GetSize());
		for (z=0;z<g_vaUnwrapArray.GetSize();z++)
			g_vaUnwrapArray[z] = 0;
	}


	g_iBinning = 2;

	if (g_bPeriodic && (g_bSDF || g_bPlProj || g_bSaveRefEnv || g_bCutCluster || g_bRevSDF) && (!g_bUnwrap))
	{
		g_bFold = true;
//		mprintf("\n    Wrapping molecules into the box.\n");
	} else if (g_bSaveJustTraj && !g_bSaveTrajNoRot && !g_bUnwrap)
	{
		if (g_bPeriodic && (!g_bUnwrap) && (!g_bSaveCoordsUnchanged))
		{
			g_bFold = AskYesNo("\n    Wrap the molecules into the box (y/n)? [yes] ",true);
		} else
		{
//			mprintf("\n    Not wrapping molecules into the box.\n");
			g_bFold = false;
		}
	} else
	{
		g_bFold = false;
//		mprintf("\n    Not wrapping molecules into the box.\n");
	}

	if (g_bAdvanced2 && g_bFold)
	{
		if (g_bSDF || g_bPlProj)
			g_bFoldAtomwise = (AskRangeInteger("    Wrap molecule-wise (0) or atom-wise (1)? [1] ",0,1,1)!=0);
		else
			g_bFoldAtomwise = (AskRangeInteger("    Wrap molecule-wise (0) or atom-wise (1)? [0] ",0,1,0)!=0);
	} else
	{
		if (g_bSDF || g_bPlProj)
			g_bFoldAtomwise = true;
		else
			g_bFoldAtomwise = false;
	}

	if (g_bFold)
	{
		if (g_bFoldAtomwise)
			mprintf("\n    Wrapping molecules into the box (atom-wise).\n\n");
				else mprintf("\n    Wrapping molecules into the box (molecule-wise).\n\n");
	} else mprintf("\n    Not wrapping molecules into the box.\n\n");

	if (g_bMSD || g_bPairMSD || g_bSaveJustTraj || g_bACF || (g_bIRSpec && g_bDipoleRefFixed) || g_bVHDF/* || g_bThermo*/)
	{
		if (g_bUnwrap || g_bBoxNonOrtho)
			g_bRemoveCOM = AskYesNo("    Remove center of mass movement of the box (y/n)? [no] ",false);
				else g_bRemoveCOM = AskYesNo("    Remove center of mass movement of the box (y/n)? [yes] ",true);
		mprintf("\n");
	} else
	{
		mprintf("    Not removing center of mass movement of the box.\n\n");
		g_bRemoveCOM = false;
	}

	if (g_bUseVelocities && (g_fTimestepLength > 1.0))
	{
		mprintf(RED,"    Warning: ");
		mprintf("Atom velocities are obtained from positions by numerical differentiation.\n");
		mprintf("             This is only reliable for trajectory timestep distances <= approx. 1 fs.\n");
		mprintf("             Your timestep distance is %.2f fs, which is larger than 1 fs.\n\n",g_fTimestepLength);

		if (!AskYesNo("    Proceed anyway (y/n)? [yes] ",true))
			return false;
	}

	if (g_TimeStep.GetCommentNumberCount() > 0)
	{
		g_bSkipDoubleSteps = AskYesNo("    Skip repeated time steps (y/n)? [no] ",false);
		if (g_bSkipDoubleSteps)
		{
			mprintf("\n");
			mprintf("    The comment line is \"%s\".\n",g_TimeStep.m_pComment);
			if (g_TimeStep.GetCommentNumberCount() > 1)
			{
				if (!AskYesNo("    Does this line contain the time step number (y/n)? [yes] ",true))
					g_bSkipDoubleSteps = false;
				if (g_bSkipDoubleSteps)
				{
//					sprintf(buf,"    Choose the time step number: ");
					buf.sprintf("    Choose the time step number: ");
					for (z=0;z<g_TimeStep.GetCommentNumberCount();z++)
					{
						if (z < g_TimeStep.GetCommentNumberCount()-1)
						{
//							sprintf(buf2,"%ld (%d), ",g_TimeStep.ExtractNumber(z),z+1);
//							strcat(buf,buf2);
							buf2.sprintf("%ld (%d), ",g_TimeStep.ExtractNumber(z),z+1);
							buf.strcat(buf2);
						} else
						{
//							sprintf(buf2,"%ld (%d)",g_TimeStep.ExtractNumber(z),z+1);
//							strcat(buf,buf2);
							buf2.sprintf("%ld (%d)",g_TimeStep.ExtractNumber(z),z+1);
							buf.strcat(buf2);
						}
					}
//					strcat(buf," [1] ");
					buf.strcat(" [1] ");
					g_iNumberPos = AskRangeInteger(buf,1,g_TimeStep.GetCommentNumberCount(),1) - 1;
				} else
				{
					if (!AskYesNo("    Is %d the time step number of the first time step (y/n)? [yes] ",true,g_TimeStep.ExtractNumber(0)))
						g_bSkipDoubleSteps = false;
				}
				if (!g_bSkipDoubleSteps)
					mprintf("    NOT skipping repeated time steps.\n");
			}
			mprintf("\n");
		}
	}

	if (g_bAdvanced2)
	{
		if (AskYesNo("    Perform a multi-interval analysis (y/n)? [no] ",false))
		{
			g_bMultiInterval = true;
			if (AskYesNo("    Use equidistant intervals (y) or enter data manually (n)? [yes] ",true))
			{
				g_iMultiIntervalBegin = AskUnsignedInteger("    In which time step to start processing the trajectory? [1] ",1) - 1;
				if (g_iMultiIntervalBegin == -1)
					g_iMultiIntervalBegin = 0;
				ti = AskUnsignedInteger("    How many intervals to create? [5] ",5);
				if (g_iTrajSteps != -1)
				{
					g_iMultiIntervalStride = AskUnsignedInteger("    Start a new interval each n-th step? [%d] ",g_iTrajSteps/ti,g_iTrajSteps/ti);
					g_iMultiIntervalLength = AskUnsignedInteger("    Enter the length in steps of one interval: [%d] ",g_iTrajSteps/ti,g_iTrajSteps/ti);
				} else
				{
					g_iMultiIntervalStride = AskUnsignedInteger("    Start a new interval each n-th step? [%d] ",10000/ti,10000/ti);
					g_iMultiIntervalLength = AskUnsignedInteger("    Enter the length in steps of one interval: [%d] ",10000/ti,10000/ti);
				}
				for (z=0;z<ti;z++)
				{
					g_laMultiIntervalStart.Add(g_iMultiIntervalBegin+z*g_iMultiIntervalStride);
					g_laMultiIntervalEnd.Add(g_iMultiIntervalBegin+z*g_iMultiIntervalStride+g_iMultiIntervalLength);
				}
			} else
			{
	_intervalagain:
				ti = AskInteger("    Enter starting time step of %d. interval: [quit] ",-1,g_laMultiIntervalStart.GetSize()+1);
				if (ti == 0)
				{
					eprintf("Please enter a number >= 1.\n");
					goto _intervalagain;
				}
				if (ti != -1)
				{
					g_laMultiIntervalStart.Add(ti-1);
					g_laMultiIntervalEnd.Add(AskUnsignedInteger_ND("    Enter ending time step of %d. interval: ",g_laMultiIntervalStart.GetSize())-1);
					goto _intervalagain;
				}
			}
			mprintf(WHITE,"\n* %d Intervals have been defined:\n",g_laMultiIntervalStart.GetSize());
			for (z=0;z<g_laMultiIntervalStart.GetSize();z++)
				mprintf("  - Interval %2d: Steps %6d - %6d.\n",z+1,g_laMultiIntervalStart[z]+1,g_laMultiIntervalEnd[z]+1);
			mprintf("\n");
		} else goto _nomulti;
	} else
	{
_nomulti:
		g_iBeginStep = AskUnsignedInteger("    In which time step to start processing the trajectory? [1] ",1) - 1;
		if (g_iBeginStep == -1)
			g_iBeginStep = 0;
		g_iMaxStep = AskUnsignedInteger("    How many time steps to use (from this position on)? [all] ",0);
		if (g_iMaxStep == 0)
			g_iMaxStep = -1; // Alle verwenden
	}
	g_iStride = AskUnsignedInteger("    Use every n-th time step from the trajectory? [1] ",1);

	if (g_iTrajSteps != -1)
	{
		if (g_iMaxStep == -1)
			mprintf(YELLOW,"\n    Using %d time steps: Every %d%s within range %d - approx. %d.\n",(g_iTrajSteps-g_iBeginStep)/g_iStride,g_iStride,(g_iStride==1)?"st":"th",g_iBeginStep+1,g_iTrajSteps);
		else mprintf(YELLOW,"\n    Using %d time steps: Every %d%s within range %d - %d.\n",g_iMaxStep/g_iStride,g_iStride,(g_iStride==1)?"st":"th",g_iBeginStep+1,g_iBeginStep+g_iMaxStep);
	} else if (g_iMaxStep != -1)
		mprintf(YELLOW,"\n    Using %d time steps: Every %d%s within range %d - %d.\n",g_iMaxStep/g_iStride,g_iStride,(g_iStride==1)?"st":"th",g_iBeginStep+1,g_iBeginStep+g_iMaxStep);
			else mprintf(YELLOW,"\n    Using every %d%s time step beginning from step %d.\n",g_iStride,(g_iStride==1)?"st":"th",g_iBeginStep+1);

	ti = -1;                                                                                                                                                                                                                                
	if (g_iTrajSteps != -1)                                                                                                                                                                                                                 
	{                                                                                                                                                                                                                                       
		if (g_iMaxStep == -1)                                                                                                                                                                                                           
			ti = (g_iTrajSteps-g_iBeginStep)/g_iStride;                                                                                                                                                                             
				else ti = g_iMaxStep/g_iStride;                                                                                                                                                                                 
	} else if (g_iMaxStep != -1)                                                                                                                                                                                                             
		ti = g_iMaxStep/g_iStride;                                                                                                                                                                                                      
                                                                                                                                                                                                                                        
	if (ti >= 0)                                                                                                                                                                                                                            
	{                                                                                                                                                                                                                                       
		if (g_bMSD)                                                                                                                                                                                                                     
		{                                                                                                                                                                                                                               
			for (z=0;z<g_oaObserv.GetSize();z++)                                                                                                                                                                                    
			{                                                                                                                                                                                                                       
				o = (CObservation*)g_oaObserv[z];                                                                                                                                                                               
				if (o->m_pMSD->m_iResolution > ti)                                                                                                                                                                              
				{                                                                                                                                                                                                               
					eprintf("\n    MSD resolution of %d requires at least %d time steps, but using only %d steps.\n    This will lead to erroneous results!\n\n",o->m_pMSD->m_iResolution,o->m_pMSD->m_iResolution,ti);     
					if (!AskYesNo("    Continue anyway (y/n)? [yes] ",true))                                                                                                                                                
						return false;                                                                                                                                                                                   
					goto _lessstepdone;                                                                                                                                                                                     
				}                                                                                                                                                                                                               
			}                                                                                                                                                                                                                       
		}                                                                                                                                                                                                                               
		if (g_bRDyn)                                                                                                                                                                                                                     
		{                                                                                                                                                                                                                               
			for (z=0;z<g_oaObserv.GetSize();z++)                                                                                                                                                                                    
			{                                                                                                                                                                                                                       
				o = (CObservation*)g_oaObserv[z];                                                                                                                                                                               
				if (o->m_pRDyn->m_iDepth*o->m_pRDyn->m_iStride > ti)                                                                                                                                                                              
				{                                                                                                                                                                                                               
					eprintf("\n    Reorientation dynamics resolution of %d requires at least %d time steps, but using only %d steps.\n    This will lead to erroneous results!\n\n",o->m_pRDyn->m_iDepth*o->m_pRDyn->m_iStride,o->m_pRDyn->m_iDepth*o->m_pRDyn->m_iStride,ti);     
					if (!AskYesNo("    Continue anyway (y/n)? [yes] ",true))                                                                                                                                                
						return false;                                                                                                                                                                                   
					goto _lessstepdone;                                                                                                                                                                                     
				}                                                                                                                                                                                                               
			}                                                                                                                                                                                                                       
		}                                                                                                                                                                                                                               
		if (g_bIRSpec)                                                                                                                                                                                                                     
		{                                                                                                                                                                                                                               
			for (z=0;z<g_oaObserv.GetSize();z++)                                                                                                                                                                                    
			{                                                                                                                                                                                                                       
				o = (CObservation*)g_oaObserv[z];                                                                                                                                                                               
				if (o->m_pIRSpec->m_iDepth*o->m_pIRSpec->m_iStride > ti)                                                                                                                                                                              
				{                                                                                                                                                                                                               
					eprintf("\n    IR spectrum resolution of %d requires at least %d time steps, but using only %d steps.\n    This will lead to erroneous results!\n\n",o->m_pIRSpec->m_iDepth*o->m_pIRSpec->m_iStride,o->m_pIRSpec->m_iDepth*o->m_pIRSpec->m_iStride,ti);     
					if (!AskYesNo("    Continue anyway (y/n)? [yes] ",true))                                                                                                                                                
						return false;                                                                                                                                                                                   
					goto _lessstepdone;                                                                                                                                                                                     
				}                                                                                                                                                                                                               
			}                                                                                                                                                                                                                       
		}                                                                                                                                                                                                                               
	}      
_lessstepdone:

	if (g_bScanVelocities)
	{
		mprintf(WHITE,"\n>>> Velocity Pre-Analysis >>>\n\n");
		g_iScanVelStart = AskUnsignedInteger("    Start in which time step for velocity pre-analysis? [0] ",0);
		g_iScanVelSteps = AskUnsignedInteger("    Use how many time steps from this point on for vel. pre-analysis? [all] ",0);
		g_iScanVelStride = AskUnsignedInteger("    Use every n-th step for velocity pre-analysis? [10] ",10);
		mprintf(WHITE,"\n<<< End of Velocity Pre-Analysis <<<\n");
	}

	if (g_bSaveRefEnv && (g_iNbhMode == 3))
	{
		mprintf(WHITE,"\n>>> Neighborhood Pre-Analysis >>>\n\n");
		g_iScanNbhStart = AskUnsignedInteger("    Start in which time step for neighborhood analysis? [0] ",0);
		g_iScanNbhSteps = AskUnsignedInteger("    Use how many time steps from this point on for nbh analysis? [all] ",0);
		g_iScanNbhStride = AskUnsignedInteger("    Use every n-th step for the neighborhood analysis? [10] ",10);
		mprintf(WHITE,"\n<<< End of Neighborhood Pre-Analysis <<<\n");
	}

	mprintf(WHITE,"\n########## All information collected ##########\n\n");
	BTOUT;
	return true;
}



/*	a = fopen("struct.xyz","wt");
	mfprintf(a,"%d\n\n",g_iGesAtomCount*20);

	FILE *b;
	for (z=0;z<20;z++)
	{
		g_pTempTimestep = new CTimeStep();
		g_pTempTimestep->CopyFrom(&g_TimeStep);
		g_pTempTimestep->CenterCOM();
		CxMatrix3 ma;
		CxVector3 vec1;
		vec1[0] = ((rand()%20001)-10000)/10000.0;
		vec1[1] = ((rand()%20001)-10000)/10000.0;
		vec1[2] = ((rand()%20001)-10000)/10000.0;
		vec1.Normalize();
		tf = (rand()%10000)/10000.0*2.0*Pi;
		ma.RotMat(vec1,tf);
		g_pTempTimestep->Transform(ma);

		vec1[0] = (rand()%20000)/10.0;
		vec1[1] = (rand()%20000)/10.0;
		vec1[2] = (rand()%20000)/10.0;

		for (z2=0;z2<g_iGesAtomCount;z2++)
			g_pTempTimestep->m_vaCoords[z2] += vec1;

		sprintf(buf,"coord_%02d.xyz",z+1);
		b = fopen(buf,"wt");
//		mfprintf(b,"%d\n\n",g_iGesAtomCount);
		for (z2=0;z2<g_iGesAtomCount;z2++)
		{
			mfprintf(a,"%s  %f  %f  %f\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z2]])->m_sName,g_pTempTimestep->m_vaCoords[z2][0]/100.0,g_pTempTimestep->m_vaCoords[z2][1]/100.0,g_pTempTimestep->m_vaCoords[z2][2]/100.0);
			mfprintf(b,"%s  %f  %f  %f\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z2]])->m_sName,g_pTempTimestep->m_vaCoords[z2][0]/100.0,g_pTempTimestep->m_vaCoords[z2][1]/100.0,g_pTempTimestep->m_vaCoords[z2][2]/100.0);
		}
		fclose(b);
	}
	fclose(a);*/
