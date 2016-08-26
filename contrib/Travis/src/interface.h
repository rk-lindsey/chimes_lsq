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

#ifndef INTERFACE_H
#define INTERFACE_H

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include "bintools.h"
#include "moltools.h"
#include "xobject.h"
#include "xobarray.h"
#include "asciiart.h"
#include "backtrace.h"
#include "grace.h"
#include "acf.h"
#include "fft.h"
#include "spectrum.h"
#include "analysis.h"
#include "analysisgroup.h"
#include "timestep.h"
#include "database.h"
#include "xintarray.h"
#include "globalvar.h"

void Interface_DefaultConf();

bool Interface_BeforeAnalysis();
bool Interface_BeforeAnalysis2();

bool Interface_Initialization();

void Interface_ProcessStep(CTimeStep *ts);

void Interface_AfterAnalysis();

void Interface_DecomposeModes(int n, CxObArray *ccr_matrix);

#endif

