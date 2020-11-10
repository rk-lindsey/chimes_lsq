/*****************************************************************************
    TRAVIS - Trajectory Analyzer and Visualizer
    http://www.travis-analyzer.de/

    Copyright (c) 2009-2016 Martin Brehm
                  2012-2016 Martin Thomas

    This file written by Martin Thomas.

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

#ifndef CONVERSION_H
#define CONVERSION_H

// Physical constants
#define CONST_ATOMIC_MASS_UNIT (1.660538921E-27) // CODATA 2010, Atomic mass unit in kg
#define CONST_BOHR_MAGNETON (927.400968E-26) // CODATA 2010, Bohr magneton in A*m^2
#define CONST_BOHR_RADIUS (0.52917721092E-10) // CODATA 2010, Bohr radius in m
#define CONST_ELEMENTARY_CHARGE (1.602176565E-19) // CODATA 2010, Elementary charge in C
#define CONST_SPEED_OF_LIGHT (299792458.0) // CODATA 2010, Vacuum speed of light in m/s

// Length conversion
#define LEN_AU2PM (52.917721092) // Atomic units to pm: CONST_BOHR_RADIUS*1E12

// Dipole moment conversion
#define DIP_DEBYE2CM (3.335640952E-30) // Debye to C*m: 1E-21/CONST_SPEED_OF_LIGHT
#define DIP_EPM2DEBYE (0.04803204506) // Elementary charge*pm to Debye: CONST_ELEMENTARY_CHARGE*CONST_SPEED_OF_LIGHT*1E9

// Total current conversion
#define CURR_AUFS2MBPM (9.142057808E-4) // Atomic unit of dipole moment/fs to Bohr magneton/pm: CONST_ELEMENTARY_CHARGE*CONST_BOHR_RADIUS*1E3/CONST_BOHR_MAGNETON
#define CURR_EMS2MBPM (1.727598547E-8) // Elementary charge*m/s to Bohr magneton/pm: CONST_ELEMENTARY_CHARGE*1E-12/CONST_BOHR_MAGNETON

// Magnetic moment conversion
#define MAG_AUPMFS2MB (9.142057808E-4) // Atomic unit of dipole moment*pm/fs to Bohr magneton: CONST_ELEMENTARY_CHARGE*CONST_BOHR_RADIUS*1E3/CONST_BOHR_MAGNETON
#define MAG_EPMMS2MB (1.727598547E-8) // Elementary charge*pm*m/s to Bohr magneton: CONST_ELEMENTARY_CHARGE*1E-12/CONST_BOHR_MAGNETON

#endif
