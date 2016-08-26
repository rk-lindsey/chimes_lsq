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

#ifndef LMWRAPPER_H
#define LMWRAPPER_H

#include "xobject.h"
#include "lmmin.h"


class CLMWrapper : public CxObject
{
public:
	void Fit_PolyExp(int degree, int n, double *x, double *y, double *par, double &r, double &integral, double *fitcurve, int maxcall);
	void Fit_ExpSpectrum(int degree, double mi, double ma, int n, double *x, double *y, double *par, int fitcurvepoints, double *fitcurvex, double *fitcurve, const char *name, int maxcall, double zeroweight, bool evolve);
	void Fit_SD(int n, double *x, double *y, double *par, double *fitcurve, double *r, int maxcall);
	CLMWrapper();
	~CLMWrapper();
};

#endif
