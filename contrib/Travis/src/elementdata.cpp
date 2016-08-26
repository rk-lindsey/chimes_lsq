/*****************************************************************************
    TRAVIS - Trajectory Analyzer and Visualizer
    http://www.travis-analyzer.de/

    Copyright (c) 2009-2016 Martin Brehm
                  2012-2016 Martin Thomas

    This file written by Martin Brehm and Philipp di Dio.

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


/****************************************************************************

    Sources of Van-der-Waals radii:

    1) R. Scott Rowland, Robin Taylor: Intermolecular Nonbonded Contact Distances in Organic Crystal Structures: Comparison with Distances Expected from van der Waals Radii. In: J. Phys. Chem. 1996, 100, S. 7384-7391, doi:10.1021/jp953141+.
    2) A. Bondi: van der Waals Volumes and Radii. In: J. Phys. Chem. 1964, 68, S. 441-451, doi:10.1021/j100785a001.
    3) Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero, Christopher J. Cramer, Donald G. Truhlar: Consistent van der Waals Radii for the Whole Main Group. In: J. Phys. Chem. A. 2009, 113, S. 5806-5812, doi:10.1021/jp8111556. 


    Coherent Neutron Scattering Cross Sections:

	Published in "Neutron News, Vol. 3, No. 3, 1992, pp. 29-37".
    Taken from "http://www.ncnr.nist.gov/resources/n-lengths/".

*****************************************************************************/

#include "travis.h"

void AddElementData()
{
	// Element, Ord. Number, Mass, Covalent Radius [pm], VdW Radius [pm], ( Coherent Neutron Scattering Cross Section [barn] )

	// 1st Period
	AddElement("H",    1,   1.01f,  37.0f, 110.0f,  1.7568f);
	AddElement("D",    1,   2.01f,  37.0f, 110.0f);
	AddElement("He",   2,   4.00f,  32.0f, 140.0f);


	// 2nd Period
	AddElement("Li",   3,   6.94f, 134.0f, 182.0f);
	AddElement("Be",   4,   9.01f,  90.0f, 153.0f);
	AddElement("B",    5,  10.81f,  90.0f, 192.0f);
	AddElement("C",    6,  12.01f,  82.0f, 170.0f,  5.551f);
	AddElement("N",    7,  14.01f,  77.0f, 155.0f, 11.01f);
	AddElement("O",    8,  16.00f,  75.0f, 152.0f,  4.232f);
	AddElement("F",    9,  18.99f,  73.0f, 147.0f,  4.017f);
	AddElement("Ne",  10,  20.18f,  69.0f, 154.0f);


	// 3rd Period
	AddElement("Na",  11,  22.99f,  71.0f, 227.0f);
	AddElement("Mg",  12,  24.31f, 130.0f, 173.0f);
	AddElement("Al",  13,  26.98f, 154.0f, 184.0f);
	AddElement("Si",  14,  28.09f, 118.0f, 210.0f);
	AddElement("P",   15,  30.97f, 111.0f, 180.0f,  3.307f);
	AddElement("S",   16,  32.06f, 106.0f, 180.0f,  1.0186f);
	AddElement("Cl",  17,  35.45f, 102.0f, 175.0f, 11.5257f);
	AddElement("Ar",  18,  39.95f,  97.0f, 188.0f);


	// 4th Period
	AddElement("K",   19,  39.10f, 196.0f, 275.0f);
	AddElement("Ca",  20,  40.08f, 174.0f, 231.0f);

	AddElement("Sc",  21,  44.96f, 144.0f,   0.0f);
	AddElement("Ti",  22,  47.90f, 136.0f,   0.0f);
	AddElement("V",   23,  50.94f, 125.0f,   0.0f);
	AddElement("Cr",  24,  52.00f, 127.0f,   0.0f);
	AddElement("Mn",  25,  54.49f, 139.0f,   0.0f);
	AddElement("Fe",  26,  55.85f, 125.0f,   0.0f);
	AddElement("Co",  27,  58.93f, 126.0f,   0.0f);
	AddElement("Ni",  28,  58.71f, 121.0f, 163.0f);
	AddElement("Cu",  29,  63.54f, 138.0f, 140.0f);
	AddElement("Zn",  30,  65.37f, 131.0f, 139.0f);

	AddElement("Ga",  31,  69.72f, 126.0f, 187.0f);
	AddElement("Ge",  32,  72.59f, 122.0f, 211.0f);
	AddElement("As",  33,  74.92f, 121.0f, 185.0f);
	AddElement("Se",  34,  78.96f, 116.0f, 190.0f);
	AddElement("Br",  35,  79.91f, 114.0f, 185.0f);
	AddElement("Kr",  36,  83.80f, 110.0f, 202.0f);


	// 5th Period
	AddElement("Rb",  37,  85.47f, 211.0f, 303.0f);
	AddElement("Sr",  38,  87.62f, 192.0f, 249.0f);

	AddElement("Y",   39,  88.91f, 162.0f,   0.0f);
	AddElement("Zr",  40,  91.22f, 148.0f,   0.0f);
	AddElement("Nb",  41,  92.91f, 137.0f,   0.0f);
	AddElement("Mo",  42,  95.94f, 145.0f,   0.0f);
	AddElement("Tc",  43,  97.00f, 131.0f,   0.0f);
	AddElement("Ru",  44, 101.07f, 126.0f,   0.0f);
	AddElement("Rh",  45, 102.90f, 135.0f,   0.0f);
	AddElement("Pd",  46, 106.40f, 131.0f, 163.0f);
	AddElement("Ag",  47, 107.87f, 153.0f, 172.0f);
	AddElement("Cd",  48, 112.40f, 148.0f, 158.0f);

	AddElement("In",  49, 114.82f, 144.0f, 193.0f);
	AddElement("Sn",  50, 118.69f, 141.0f, 217.0f);
	AddElement("Sb",  51, 121.75f, 138.0f, 206.0f);
	AddElement("Te",  52, 127.60f, 135.0f, 206.0f);
	AddElement("I",   53, 126.90f, 133.0f, 198.0f);
	AddElement("Xe",  54, 131.30f, 130.0f, 216.0f);


	// 6th Period
	AddElement("Cs",  55, 132.91f, 225.0f, 343.0f);
	AddElement("Ba",  56, 137.34f, 198.0f, 268.0f);

	AddElement("La",  57, 138.91f, 169.0f,   0.0f);

	AddElement("Ce",  58, 140.12f, 204.0f,   0.0f);
	AddElement("Pr",  59, 140.91f, 203.0f,   0.0f);
	AddElement("Nd",  60, 144.24f, 201.0f,   0.0f);
	AddElement("Pm",  61, 146.90f, 199.0f,   0.0f);
	AddElement("Sm",  62, 150.36f, 198.0f,   0.0f);
	AddElement("Eu",  63, 151.96f, 198.0f,   0.0f);
	AddElement("Gd",  64, 157.25f, 196.0f,   0.0f);
	AddElement("Tb",  65, 158.93f, 194.0f,   0.0f);
	AddElement("Dy",  66, 162.50f, 192.0f,   0.0f);
	AddElement("Ho",  67, 164.93f, 192.0f,   0.0f);
	AddElement("Er",  68, 167.26f, 189.0f,   0.0f);
	AddElement("Tm",  69, 168.93f, 190.0f,   0.0f);
	AddElement("Yb",  70, 173.05f, 187.0f,   0.0f);
	AddElement("Lu",  71, 174.97f, 187.0f,   0.0f);

	AddElement("Hf",  72, 178.49f, 150.0f,   0.0f);
	AddElement("Ta",  73, 180.95f, 138.0f,   0.0f);
	AddElement("W",   74, 183.85f, 146.0f,   0.0f);
	AddElement("Re",  75, 186.20f, 159.0f,   0.0f);
	AddElement("Os",  76, 190.20f, 128.0f,   0.0f);
	AddElement("Ir",  77, 192.20f, 137.0f,   0.0f);
	AddElement("Pt",  78, 195.09f, 138.0f,   0.0f);
	AddElement("Au",  79, 196.97f, 144.0f,   0.0f);
	AddElement("Hg",  80, 200.59f, 149.0f,   0.0f);

	AddElement("Tl",  81, 204.37f, 148.0f,   0.0f);
	AddElement("Pb",  82, 207.19f, 146.0f,   0.0f);
	AddElement("Bi",  83, 208.98f, 146.0f,   0.0f);
	AddElement("Po",  84, 209.00f, 140.0f,   0.0f);
	AddElement("At",  85, 210.00f, 145.0f,   0.0f);
	AddElement("Rn",  86, 222.00f, 145.0f,   0.0f);


	/* 7th period */
	AddElement("Fr",  87, 223.00f, 260.0f,   0.0f);
	AddElement("Ra",  88, 226.03f, 221.0f,   0.0f);

	AddElement("Ac",  89, 227.00f, 215.0f,   0.0f);

	AddElement("Th",  90, 232.04f, 206.0f,   0.0f);
	AddElement("Pa",  91, 231.04f, 200.0f,   0.0f);
	AddElement("U",   92, 238.03f, 196.0f,   0.0f);
	AddElement("Np",  93, 237.05f, 190.0f,   0.0f);
	AddElement("Pu",  94, 244.10f, 187.0f,   0.0f);
	AddElement("Am",  95, 243.10f, 180.0f,   0.0f);
	AddElement("Cm",  96, 247.10f, 169.0f,   0.0f);
	AddElement("Bk",  97, 247.10f, 160.0f,   0.0f);
	AddElement("Cf",  98, 251.10f, 160.0f,   0.0f);
	AddElement("Es",  99, 254.10f, 160.0f,   0.0f);
	AddElement("Fm", 100, 257.10f, 160.0f,   0.0f);
	AddElement("Md", 101, 258.00f, 160.0f,   0.0f);
	AddElement("No", 102, 259.00f, 160.0f,   0.0f);
	AddElement("Lr", 103, 260.00f, 160.0f,   0.0f);


	// Virtual Atom
	AddElement("#",    0,   0.00f,   0.0f,   0.0f);

	// Colors for some common atoms. Other atoms have standard color.
	SetElementColor("H",  255, 255, 255, 150, 150, 150 );
	SetElementColor("C",  228, 113,   0, 180, 180, 180 );
	SetElementColor("N",    0,   0, 255, 140, 140, 140 );
	SetElementColor("O",  255,   0,   0, 160, 160, 160 );
	SetElementColor("S",  255, 255,   0, 200, 200, 200 );
}

