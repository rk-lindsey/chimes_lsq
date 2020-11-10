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

#ifndef RANDOM_H
#define RANDOM_H


#include "xobject.h"
#include "ziggurat.h"
#include <stdlib.h>


class CRandom : public CxObject
{
public:
	CRandom();
	~CRandom();

	float RandomUniform()
	{
//		return r4_uni(&seed);
		return ((float)(rand()%RAND_MAX)/RAND_MAX/RAND_MAX) + ((float)(rand()%RAND_MAX)/RAND_MAX);
	}

	float RandomNormal()
	{
		return r4_nor(&seed,kn,fn,wn);
	}

	float RandomExp()
	{
		return r4_exp(&seed,ke,fe,we);
	}

	float fn[128];
	int kn[128];
	float wn[128];
	float fe[256];
	int ke[256];
	float we[256];
	unsigned long int seed;
};

#endif
