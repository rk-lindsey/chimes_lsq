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


#ifndef CONFIG_H
#define CONFIG_H


#define SOURCE_VERSION "Feb 23 2016"

/* Please uncomment / comment out the flags you want to use / not to use. */

// You have to change this according to your target Operating System.
// For a generic platform-independent build, comment out both lines.
//#define TARGET_WINDOWS
#define TARGET_LINUX

// Use color for screen output?
#define USE_COLOR

// Maximum number of bonds any atom may form
#define MAX_BONDS 12

// Use the FFTW library? Otherwise, the built-in KISS-FFT routine is used (default)
//#define USE_FFTW

// Bound checking in dynamic arrays?
#define DEBUG_ARRAYS


// Handle volumetric data in single or double precision? Active=double, inactive=single
#define VORI_DOUBLE


#ifdef VORI_DOUBLE
  #define VORI_FLOAT double
#else
  #define VORI_FLOAT float
#endif


/***************************************************************************/
// Some "hardcore" debug flags
// Warning: They may drastically decrease performance
/***************************************************************************/


// Activate debug backtrace?
//#define DEBUG_BACKTRACE
//#define DEBUG_EXTENDED_BACKTRACE


//#define DEBUG_DATABASE
//#define DEBUG_CSTRING
//#define DEBUG_MATULTRA
//#define DEBUG_CVECTOR3
//#define DEBUG_CDVECTOR3


//#define DEBUG_COBARRAY
//#define DEBUG_CPTRARRAY
//#define DEBUG_CWORDARRAY
//#define DEBUG_CLONGARRAY
//#define DEBUG_CINTARRAY
//#define DEBUG_CFLOATARRAY
//#define DEBUG_CDOUBLEARRAY
//#define DEBUG_CVEC3ARRAY
//#define DEBUG_CDVEC3ARRAY


#endif

