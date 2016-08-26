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

#include "bicgstab.h"

#include "3df.h"
#include "globalvar.h"

CSparseMatrix::CSparseMatrix() {
	m_size = 0;
	m_val = NULL;
	m_colInd = NULL;
	m_rowPtr = NULL;
}

CSparseMatrix::CSparseMatrix(const CSparseMatrix &matrix) {
	m_size = matrix.m_size;
	try { m_val = new double[matrix.m_rowPtr[matrix.m_size]]; } catch(...) { m_val = NULL; }
	if (m_val == NULL) NewException((double)matrix.m_rowPtr[matrix.m_size] * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	memcpy(m_val, matrix.m_val, matrix.m_rowPtr[matrix.m_size] * sizeof(double));
	try { m_colInd = new int[matrix.m_rowPtr[matrix.m_size]]; } catch(...) { m_colInd = NULL; }
	if (m_colInd == NULL) NewException((double)matrix.m_rowPtr[matrix.m_size] * sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	memcpy(m_colInd, matrix.m_colInd, matrix.m_rowPtr[matrix.m_size] * sizeof(int));
	try { m_rowPtr = new int[matrix.m_size + 1]; } catch(...) { m_rowPtr = NULL; }
	if (m_rowPtr == NULL) NewException((double)(matrix.m_size + 1) * sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	memcpy(m_rowPtr, matrix.m_rowPtr, (matrix.m_size + 1) * sizeof(int));
}

CSparseMatrix::~CSparseMatrix() {
	if (m_val != NULL)
		delete[] m_val;
	if (m_colInd != NULL)
		delete[] m_colInd;
	if (m_rowPtr != NULL)
		delete[] m_rowPtr;
}

void CSparseMatrix::setTest() {
	m_size = 9;
	m_val = new double[45];
	m_colInd = new int[45];
	m_rowPtr = new int[10];
	m_val[0] = 8.0;
	m_val[1] = -0.5;
	m_val[2] = -1.5;
	m_val[3] = -0.5;
	m_val[4] = -1.5;
	m_val[5] = -1.5;
	m_val[6] = 3.0;
	m_val[7] = -0.5;
	m_val[8] = -0.5;
	m_val[9] = -1.5;
	m_val[10] = -0.5;
	m_val[11] = -1.5;
	m_val[12] = 4.0;
	m_val[13] = -0.5;
	m_val[14] = -1.5;
	m_val[15] = -1.5;
	m_val[16] = 2.0;
	m_val[17] = -0.5;
	m_val[18] = -1.5;
	m_val[19] = -0.5;
	m_val[20] = -1.5;
	m_val[21] = -1.5;
	m_val[22] = 1.0;
	m_val[23] = -0.5;
	m_val[24] = -0.5;
	m_val[25] = -1.5;
	m_val[26] = -0.5;
	m_val[27] = -1.5;
	m_val[28] = 3.0;
	m_val[29] = -0.5;
	m_val[30] = -0.5;
	m_val[31] = -1.5;
	m_val[32] = 5.0;
	m_val[33] = -0.5;
	m_val[34] = -1.5;
	m_val[35] = -0.5;
	m_val[36] = -1.5;
	m_val[37] = -1.5;
	m_val[38] = 6.0;
	m_val[39] = -0.5;
	m_val[40] = -0.5;
	m_val[41] = -1.5;
	m_val[42] = -0.5;
	m_val[43] = -1.5;
	m_val[44] = 2.0;
	m_colInd[0] = 0;
	m_colInd[1] = 1;
	m_colInd[2] = 2;
	m_colInd[3] = 3;
	m_colInd[4] = 6;
	m_colInd[5] = 0;
	m_colInd[6] = 1;
	m_colInd[7] = 2;
	m_colInd[8] = 4;
	m_colInd[9] = 7;
	m_colInd[10] = 0;
	m_colInd[11] = 1;
	m_colInd[12] = 2;
	m_colInd[13] = 5;
	m_colInd[14] = 8;
	m_colInd[15] = 0;
	m_colInd[16] = 3;
	m_colInd[17] = 4;
	m_colInd[18] = 5;
	m_colInd[19] = 6;
	m_colInd[20] = 1;
	m_colInd[21] = 3;
	m_colInd[22] = 4;
	m_colInd[23] = 5;
	m_colInd[24] = 7;
	m_colInd[25] = 2;
	m_colInd[26] = 3;
	m_colInd[27] = 4;
	m_colInd[28] = 5;
	m_colInd[29] = 8;
	m_colInd[30] = 0;
	m_colInd[31] = 3;
	m_colInd[32] = 6;
	m_colInd[33] = 7;
	m_colInd[34] = 8;
	m_colInd[35] = 1;
	m_colInd[36] = 4;
	m_colInd[37] = 6;
	m_colInd[38] = 7;
	m_colInd[39] = 8;
	m_colInd[40] = 2;
	m_colInd[41] = 5;
	m_colInd[42] = 6;
	m_colInd[43] = 7;
	m_colInd[44] = 8;
	m_rowPtr[0] = 0;
	m_rowPtr[1] = 5;
	m_rowPtr[2] = 10;
	m_rowPtr[3] = 15;
	m_rowPtr[4] = 20;
	m_rowPtr[5] = 25;
	m_rowPtr[6] = 30;
	m_rowPtr[7] = 35;
	m_rowPtr[8] = 40;
	m_rowPtr[9] = 45;
}

void CSparseMatrix::matrixVectorProduct(const double *vect, double *result) const {
	int i;
	for (i = 0; i < m_size; i++) {
		result[i] = 0.0;
		int j;
		for (j = m_rowPtr[i]; j < m_rowPtr[i + 1]; j++) {
			result[i] += m_val[j] * vect[m_colInd[j]];
		}
	}
}

void CCurrentPDEDiscretizer::discretize(CSparseMatrix *pdeMatrix, const C3DF<VORI_FLOAT> &electronDensity, const C3DF<VORI_FLOAT> &electronDensityGradientX, const C3DF<VORI_FLOAT> &electronDensityGradientY, const C3DF<VORI_FLOAT> &electronDensityGradientZ) {
	if (pdeMatrix->m_val != NULL)
		delete[] pdeMatrix->m_val;
	if (pdeMatrix->m_colInd != NULL)
		delete[] pdeMatrix->m_colInd;
	if (pdeMatrix->m_rowPtr != NULL)
		delete[] pdeMatrix->m_rowPtr;
	
	pdeMatrix->m_size = electronDensity.m_iRes[0] * electronDensity.m_iRes[1] * electronDensity.m_iRes[2];
	int numEntries = 7 * pdeMatrix->m_size;
	try { pdeMatrix->m_val = new double[numEntries]; } catch(...) { pdeMatrix->m_val = NULL; }
	if (pdeMatrix->m_val == NULL) NewException((double)numEntries * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	try { pdeMatrix->m_colInd = new int[numEntries]; } catch(...) { pdeMatrix->m_colInd = NULL; }
	if (pdeMatrix->m_colInd == NULL) NewException((double)numEntries * sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	try { pdeMatrix->m_rowPtr = new int[pdeMatrix->m_size + 1]; } catch(...) { pdeMatrix->m_rowPtr = NULL; }
	if (pdeMatrix->m_rowPtr == NULL) NewException((double)(pdeMatrix->m_size + 1) * sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
// 	mprintf(GREEN, "%#.10g %#.10g %#.10g %#.10g %#.10g\n", electronDensity.m_pBin[0], electronDensityGradientX.m_pBin[0], electronDensityGradientY.m_pBin[0], electronDensityGradientZ.m_pBin[0], g_fCubeXStep);
	
	int j = 0;
	int i;
	for (i = 0; i < pdeMatrix->m_size; i++) {
		pdeMatrix->m_rowPtr[i] = j;
		if (i >= electronDensity.m_iRes[0] * electronDensity.m_iRes[1] * (electronDensity.m_iRes[2] - 1)) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeZStep * g_fCubeZStep) + 0.5 * electronDensityGradientZ.m_pBin[i] / g_fCubeZStep;
			pdeMatrix->m_colInd[j] = i - electronDensity.m_iRes[0] * electronDensity.m_iRes[1] * (electronDensity.m_iRes[2] - 1);
			j++;
		}
		if (i >= electronDensity.m_iRes[0] * electronDensity.m_iRes[1]) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeZStep * g_fCubeZStep) - 0.5 * electronDensityGradientZ.m_pBin[i] / g_fCubeZStep;
			pdeMatrix->m_colInd[j] = i - electronDensity.m_iRes[0] * electronDensity.m_iRes[1];
// 			if (j == 13600000) {
// 				mprintf(GREEN, "%#.10g %#.10g\n", electronDensity.m_pBin[i], electronDensityGradientZ.m_pBin[i]);
// 			}
			j++;
		}
		if (i % (electronDensity.m_iRes[0] * electronDensity.m_iRes[1]) >= electronDensity.m_iRes[0] * (electronDensity.m_iRes[1] - 1)) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeYStep * g_fCubeYStep) + 0.5 * electronDensityGradientY.m_pBin[i] / g_fCubeYStep;
			pdeMatrix->m_colInd[j] = i - electronDensity.m_iRes[0] * (electronDensity.m_iRes[1] - 1);
			j++;
		}
		if (i % (electronDensity.m_iRes[0] * electronDensity.m_iRes[1]) >= electronDensity.m_iRes[0]) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeYStep * g_fCubeYStep) - 0.5 * electronDensityGradientY.m_pBin[i] / g_fCubeYStep;
			pdeMatrix->m_colInd[j] = i - electronDensity.m_iRes[0];
			j++;
		}
		if (i % electronDensity.m_iRes[0] == electronDensity.m_iRes[0] - 1) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeXStep * g_fCubeXStep) + 0.5 * electronDensityGradientX.m_pBin[i] / g_fCubeXStep;
			pdeMatrix->m_colInd[j] = i - electronDensity.m_iRes[0] + 1;
			j++;
		}
		if (i % electronDensity.m_iRes[0] != 0) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeXStep * g_fCubeXStep) - 0.5 * electronDensityGradientX.m_pBin[i] / g_fCubeXStep;
			pdeMatrix->m_colInd[j] = i - 1;
			j++;
		}
		pdeMatrix->m_val[j] = -2.0 * electronDensity.m_pBin[i] / (g_fCubeXStep * g_fCubeXStep) - 2.0 * electronDensity.m_pBin[i] / (g_fCubeYStep * g_fCubeYStep) - 2.0 * electronDensity.m_pBin[i] / (g_fCubeZStep * g_fCubeZStep);
		pdeMatrix->m_colInd[j] = i;
		j++;
		if (i % electronDensity.m_iRes[0] != electronDensity.m_iRes[0] - 1) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeXStep * g_fCubeXStep) + 0.5 * electronDensityGradientX.m_pBin[i] / g_fCubeXStep;
			pdeMatrix->m_colInd[j] = i + 1;
			j++;
		}
		if (i % electronDensity.m_iRes[0] == 0) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeXStep * g_fCubeXStep) - 0.5 * electronDensityGradientX.m_pBin[i] / g_fCubeXStep;
			pdeMatrix->m_colInd[j] = i + electronDensity.m_iRes[0] - 1;
			j++;
		}
		if (i % (electronDensity.m_iRes[0] * electronDensity.m_iRes[1]) < electronDensity.m_iRes[0] * (electronDensity.m_iRes[1] - 1)) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeYStep * g_fCubeYStep) + 0.5 * electronDensityGradientY.m_pBin[i] / g_fCubeYStep;
			pdeMatrix->m_colInd[j] = i + electronDensity.m_iRes[0];
			j++;
		}
		if (i % (electronDensity.m_iRes[0] * electronDensity.m_iRes[1]) < electronDensity.m_iRes[0]) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeYStep * g_fCubeYStep) - 0.5 * electronDensityGradientY.m_pBin[i] / g_fCubeYStep;
			pdeMatrix->m_colInd[j] = i + electronDensity.m_iRes[0] * (electronDensity.m_iRes[1] - 1);
			j++;
		}
		if (i < electronDensity.m_iRes[0] * electronDensity.m_iRes[1] * (electronDensity.m_iRes[2] - 1)) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeZStep * g_fCubeZStep) + 0.5 * electronDensityGradientZ.m_pBin[i] / g_fCubeZStep;
			pdeMatrix->m_colInd[j] = i + electronDensity.m_iRes[0] * electronDensity.m_iRes[1];
			j++;
		}
		if (i < electronDensity.m_iRes[0] * electronDensity.m_iRes[1]) {
			pdeMatrix->m_val[j] = 1.0 * electronDensity.m_pBin[i] / (g_fCubeZStep * g_fCubeZStep) - 0.5 * electronDensityGradientZ.m_pBin[i] / g_fCubeZStep;
			pdeMatrix->m_colInd[j] = i + electronDensity.m_iRes[0] * electronDensity.m_iRes[1] * (electronDensity.m_iRes[2] - 1);
			j++;
		}
	}
	if (j != numEntries) {
		eprintf("Error in PDE discretizer\n");
		abort();
	}
	pdeMatrix->m_rowPtr[pdeMatrix->m_size] = numEntries;
// 	mprintf(GREEN, "%#.10g\n", pdeMatrix->m_val[13600000]);
}

// For this implementation, the Fortran code provided at http://www.staff.science.uu.nl/~vorst102/software.html (November 07, 2014) was adopted
/***********************************************************
 * This code is based on:                                  *
 *                                                         *
 * subroutine bicgstab2                                    *
 * Copyright (c) 1998 by M.A.Botchev                       *
 * Permission to copy all or part of this work is granted, *
 * provided that the copies are not made or distributed    *
 * for resale, and that the copyright notice and this      *
 * notice are retained.                                    *
 *                                                         *
 * subroutine bistbl                                       *
 * Copyright (c) 1995 by D.R. Fokkema.                     *
 * Permission to copy all or part of this work is granted, *
 * provided that the copies are not made or distributed    *
 * for resale, and that the copyright notice and this      *
 * notice are retained.                                    *
 ***********************************************************/
bool CCurrentPDESolver::bicgstabl(const int l, const CSparseMatrix *pdeMatrix, double *solution, const double *rightHandSide, const int maxIter, double *thresh, FILE *infoFile) {
	double **r;
	try { r = new double *[l + 1]; } catch(...) { r = NULL; }
	if (r == NULL) NewException((double)(l + 1) * sizeof(double *), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	double **u;
	try { u = new double *[l + 1]; } catch(...) { u = NULL; }
	if (u == NULL) NewException((double)(l + 1) * sizeof(double *), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	int i;
	for (i = 0; i < l + 1; i++) {
		try { r[i] = new double[pdeMatrix->m_size]; } catch(...) { r[i] = NULL; }
		if (r[i] == NULL) NewException((double)pdeMatrix->m_size * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		try { u[i] = new double[pdeMatrix->m_size]; } catch(...) { u[i] = NULL; }
		if (u[i] == NULL) NewException((double)pdeMatrix->m_size * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	double *r0;
	try { r0 = new double[pdeMatrix->m_size]; } catch(...) { r0 = NULL; }
	if (r0 == NULL) NewException((double)pdeMatrix->m_size * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	double *b;
	try { b = new double[pdeMatrix->m_size]; } catch(...) { b = NULL; }
	if (b == NULL) NewException((double)pdeMatrix->m_size * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	double *x;
	try { x = new double[pdeMatrix->m_size]; } catch(...) { x = NULL; }
	if (x == NULL) NewException((double)pdeMatrix->m_size * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	double alpha, beta, omega, sigma, rho0, rho1;
	double *rwork;
	try { rwork = new double[(l + 1) * (3 + l + 1)]; } catch(...) { rwork = NULL; }
	if (rwork == NULL) NewException((double)(l + 1) * (3 + l + 1) * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	double *rtemp;
	try { rtemp = new double[(l - 1) * (l - 1)]; } catch(...) { rtemp = NULL; }
	if (rtemp == NULL) NewException((double)(l - 1) * (l - 1) * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	int *rperm;
	try { rperm = new int[l - 1]; } catch(...) { rperm = NULL; }
	if (rperm == NULL) NewException((double)(l - 1) * sizeof(int), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	int z = 0;
	int y0 = z + (l + 1);
	int y1 = y0 + 1;
	int y = y1 + 1;
	double kappa0, kappal, varrho, hatgamma;
	double resNorm0, maxNormX, maxNormR;
	double retVal;
	
	if (infoFile != NULL) {
		fprintf(infoFile, "Starting BICGSTABL(%d) with %d maximum iterations and %g threshold\n", l, maxIter, *thresh);
	}
	
	CSparseMatrix precondMatrix(*pdeMatrix);
	CILUPreconditioner::precondition(&precondMatrix);
	
	pdeMatrix->matrixVectorProduct(solution, r[0]);
	for (i = 0; i < pdeMatrix->m_size; i++) {
		r[0][i] = rightHandSide[i] - r[0][i];
	}
	CILUPreconditioner::solve(&precondMatrix, r[0]);
// 	mprintf(GREEN, "%#.10g\n", r[0][10]);
	
	double resNorm = 0.0;
	for (i = 0; i < pdeMatrix->m_size; i++) {
		resNorm += r[0][i] * r[0][i];
		r0[i] = r[0][i];
		b[i] = r[0][i];
		x[i] = solution[i];
		solution[i] = 0.0;
	}
	resNorm = sqrt(resNorm);
	if (infoFile != NULL) {
		fprintf(infoFile, "%5d %20.14g\n", 0, resNorm);
		fflush(infoFile);
	}
	resNorm0 = resNorm;
	maxNormX = resNorm0;
	maxNormR = resNorm0;
	retVal = resNorm0;
	
	alpha = 0.0;
	omega = 1.0;
	sigma = 1.0;
	rho0 = 1.0;
	
	int numIter = 0;
	while (numIter < maxIter && resNorm > *thresh) {
		numIter++;
		rho0 *= -omega;
		for (i = 0; i < l; i++) {
			rho1 = 0.0;
			int j;
			for (j = 0; j < pdeMatrix->m_size; j++) {
				rho1 += r0[j] * r[i][j];
			}
			if (rho0 == 0.0) {
				eprintf("BiCGSTAB(l): rho0 is zero!\n");
				abort();
			}
			beta = alpha * (rho1 / rho0);
			rho0 = rho1;
			for (j = 0; j <= i; j++) {
				int k;
				for (k = 0; k < pdeMatrix->m_size; k++) {
					u[j][k] = r[j][k] - beta * u[j][k];
				}
			}
			pdeMatrix->matrixVectorProduct(u[i], u[i + 1]);
			CILUPreconditioner::solve(&precondMatrix, u[i + 1]);
			sigma = 0.0;
			for (j = 0; j < pdeMatrix->m_size; j++) {
				sigma += r0[j] * u[i + 1][j];
			}
			if (sigma == 0.0) {
				eprintf("BiCGSTAB(l): sigma is zero!\n");
				abort();
			}
			alpha = rho1 / sigma;
			for (j = 0; j < pdeMatrix->m_size; j++) {
				solution[j] += alpha * u[0][j];
			}
			for (j = 0; j <= i; j++) {
				int k;
				for (k = 0; k < pdeMatrix->m_size; k++) {
					r[j][k] -= alpha * u[j + 1][k];
				}
			}
			pdeMatrix->matrixVectorProduct(r[i], r[i + 1]);
			CILUPreconditioner::solve(&precondMatrix, r[i + 1]);
			resNorm = 0.0;
			for (j = 0; j < pdeMatrix->m_size; j++) {
				resNorm += r[0][j] * r[0][j];
			}
			resNorm = sqrt(resNorm);
// 			mprintf(GREEN, "%#.14g %#.10g\n", resNorm, r[0][1942857]);
			if (resNorm > maxNormX)
				maxNormX = resNorm;
			if (resNorm > maxNormR)
				maxNormR = resNorm;
		}
		
		for (i = 0; i < l + 1; i++) {
			int j;
			for (j = 0; j < l + 1 - i; j++) {
				rwork[(z + i) * (l + 1) + i + j] = 0.0;
				int k;
				for (k = 0; k < pdeMatrix->m_size; k++) {
					rwork[(z + i) * (l + 1) + i + j] += r[i + j][k] * r[i][k];
				}
			}
			for (j = 0; j < l - i; j++) {
				rwork[(z + i + 1 + j) * (l + 1) + i] = rwork[(z + i) * (l + 1) + i + 1 + j];
			}
		}
		for (i = 0; i < (l - 1); i++) {
			int j;
			for (j = 0; j < (l - 1); j++) {
				rtemp[i * (l - 1) + j] = rwork[(z + i + 1) * (l + 1) + j + 1];
			}
		}
		luFactorization(rtemp, rperm, l - 1);
		rwork[y0 * (l + 1)] = -1.0;
		for (i = 0; i < l - 1; i++) {
			rwork[y0 * (l + 1) + i + 1] = rwork[z * (l + 1) + i + 1];
		}
		forwardBackward(rtemp, rperm, l - 1, &rwork[y0 * (l + 1) + 1]);
		rwork[y0 * (l + 1) + l] = 0.0;
		rwork[y1 * (l + 1)] = 0.0;
		for (i = 0; i < l - 1; i++) {
			rwork[y1 * (l + 1) + i + 1] = rwork[(z + l) * (l + 1) + i + 1];
		}
		forwardBackward(rtemp, rperm, l - 1, &rwork[y1 * (l + 1) + 1]);
		rwork[y1 * (l + 1) + l] = -1.0;
		for (i = 0; i < l + 1; i++) {
			rwork[y * (l + 1) + i] = 0.0;
			int j;
			for (j = 0; j < l + 1; j++) {
				rwork[y * (l + 1) + i] += rwork[(z + i) * (l + 1) + j] * rwork[y0 * (l + 1) + j];
			}
		}
		kappa0 = 0.0;
		for (i = 0; i < l + 1; i++) {
			kappa0 += rwork[y0 * (l + 1) + i] * rwork[y * (l + 1) + i];
		}
		kappa0 = sqrt(kappa0);
		for (i = 0; i < l + 1; i++) {
			rwork[y * (l + 1) + i] = 0.0;
			int j;
			for (j = 0; j < l + 1; j++) {
				rwork[y * (l + 1) + i] += rwork[(z + i) * (l + 1) + j] * rwork[y1 * (l + 1) + j];
			}
		}
		kappal = 0.0;
		for (i = 0; i < l + 1; i++) {
			kappal += rwork[y1 * (l + 1) + i] * rwork[y * (l + 1) + i];
		}
		kappal = sqrt(kappal);
		for (i = 0; i < l + 1; i++) {
			rwork[y * (l + 1) + i] = 0.0;
			int j;
			for (j = 0; j < l + 1; j++) {
				rwork[y * (l + 1) + i] += rwork[(z + i) * (l + 1) + j] * rwork[y0 * (l + 1) + j];
			}
		}
		varrho = 0.0;
		for (i = 0; i < l + 1; i++) {
			varrho += rwork[y1 * (l + 1) + i] * rwork[y * (l + 1) + i];
		}
		varrho /= (kappa0 * kappal);

		if (fabs(varrho) < 0.7)
#ifdef TARGET_WINDOWS
			hatgamma = _copysign(1.0, varrho) * 0.7 * (kappa0 / kappal);
#else
			hatgamma = copysign(1.0, varrho) * 0.7 * (kappa0 / kappal);
#endif
		else
#ifdef TARGET_WINDOWS
			hatgamma = _copysign(1.0, varrho) * fabs(varrho) * (kappa0 / kappal);
#else
			hatgamma = copysign(1.0, varrho) * fabs(varrho) * (kappa0 / kappal);
#endif

		for (i = 0; i < l + 1; i++) {
			rwork[y0 * (l + 1) + i] -= hatgamma * rwork[y1 * (l + 1) + i];
		}
		omega = rwork[y0 * (l + 1) + l];
		for (i = 0; i < pdeMatrix->m_size; i++) {
			int j;
			for (j = 0; j < l; j++) {
				u[0][i] -= u[j + 1][i] * rwork[y0 * (l + 1) + j + 1];
			}
		}
		for (i = 0; i < pdeMatrix->m_size; i++) {
			int j;
			for (j = 0; j < l; j++) {
				solution[i] += r[j][i] * rwork[y0 * (l + 1) + j + 1];
			}
		}
		for (i = 0; i < pdeMatrix->m_size; i++) {
			int j;
			for (j = 0; j < l; j++) {
				r[0][i] -= r[j + 1][i] * rwork[y0 * (l + 1) + j + 1];
			}
		}
		for (i = 0; i < l + 1; i++) {
			rwork[y * (l + 1) + i] = 0.0;
			int j;
			for (j = 0; j < l + 1; j++) {
				rwork[y * (l + 1) + i] += rwork[(z + i) * (l + 1) + j] * rwork[y0 * (l + 1) + j];
			}
		}
		resNorm = 0.0;
		for (i = 0; i < l + 1; i++) {
			resNorm += rwork[y0 * (l + 1) + i] * rwork[y * (l + 1) + i];
		}
		resNorm = sqrt(resNorm);
		if (resNorm < retVal)
			retVal = resNorm;
		if (infoFile != NULL) {
			fprintf(infoFile, "%5d %20.14g\n", numIter, resNorm);
			fflush(infoFile);
		}
		if (resNorm > maxNormX)
			maxNormX = resNorm;
		if (resNorm > maxNormR)
			maxNormR = resNorm;
		
		if ((resNorm < 0.01 * maxNormR && resNorm0 < maxNormR) || (resNorm < 0.01 * resNorm0 && resNorm0 < maxNormX)) {
			if (infoFile != NULL) {
				fprintf(infoFile, "Updating residuals\n");
			}
			pdeMatrix->matrixVectorProduct(solution, r[0]);
			CILUPreconditioner::solve(&precondMatrix, r[0]);
			for (i = 0; i < pdeMatrix->m_size; i++) {
				r[0][i] = b[i] - r[0][i];
			}
			maxNormR = resNorm;
			if (resNorm < 0.01 * resNorm0 && resNorm0 < maxNormX) {
				if (infoFile != NULL) {
					fprintf(infoFile, "On-the-fly restart\n");
				}
				for (i = 0; i < pdeMatrix->m_size; i++) {
					x[i] += solution[i];
					solution[i] = 0.0;
					b[i] = r[0][i];
				}
				maxNormX = resNorm;
			}
		}
	}
	
	for (i = 0; i < pdeMatrix->m_size; i++) {
		solution[i] += x[i];
	}
	
	pdeMatrix->matrixVectorProduct(solution, r[0]);
	for (i = 0; i < pdeMatrix->m_size; i++) {
		r[0][i] = rightHandSide[i] - r[0][i];
	}
	CILUPreconditioner::solve(&precondMatrix, r[0]);
	resNorm = 0.0;
	for (i = 0; i < pdeMatrix->m_size; i++) {
		resNorm += r[0][i] * r[0][i];
	}
	resNorm = sqrt(resNorm);
	if (infoFile != NULL) {
		fprintf(infoFile, "Final %20.14g\n", resNorm);
		if (numIter == maxIter) {
			fprintf(infoFile, "Not converged\n");
		}
	}
	
	for (i = 0; i < l + 1; i++) {
		delete[] r[i];
		delete[] u[i];
	}
	delete[] r;
	delete[] u;
	delete[] r0;
	delete[] b;
	delete[] x;
	delete[] rwork;
	delete[] rtemp;
	delete[] rperm;
	
	*thresh = retVal;
	return !(numIter == maxIter);
}

double CCurrentPDESolver::calcResidual(const CSparseMatrix *pdeMatrix, const double *solution, const double *rightHandSide) {
	double *temp;
	try { temp = new double[pdeMatrix->m_size]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)pdeMatrix->m_size * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	CSparseMatrix precondMatrix(*pdeMatrix);
	CILUPreconditioner::precondition(&precondMatrix);
	
	pdeMatrix->matrixVectorProduct(solution, temp);
	int i;
	for (i = 0; i < pdeMatrix->m_size; i++) {
		temp[i] = rightHandSide[i] - temp[i];
	}
	CILUPreconditioner::solve(&precondMatrix, temp);
	double resNorm = 0.0;
	for (i = 0; i < pdeMatrix->m_size; i++) {
		resNorm += temp[i] * temp[i];
	}
	resNorm = sqrt(resNorm);
	
	return resNorm;
}

void CCurrentPDESolver::luFactorization(double *matrix, int *perm, int size) {
	int i;
	for (i = 0; i < size; i++) {
		double max = fabs(matrix[i * size + i]);
		int maxRow = i;
		int j;
		for (j = i; j < size; j++) {
			if (fabs(matrix[j * size + i]) > max) {
				max = fabs(matrix[j * size + i]);
				maxRow = j;
			}
		}
		if (maxRow != i) {
			int j;
			for (j = 0; j < size; j++) {
				double temp = matrix[i * size + j];
				matrix[i * size + j] = matrix[maxRow * size + j];
				matrix[maxRow * size + j] = temp;
			}
		}
		perm[i] = maxRow;
		int k;
		for (k = 0; k < i; k++) {
			matrix[i * size + k] /= matrix[k * size + k];
			int j;
			for (j = k + 1; j < size; j++) {
				matrix[i * size + j] -= matrix[i * size + k] * matrix[k * size + j];
			}
		}
	}
}

void CCurrentPDESolver::forwardBackward(double *luMatrix, int *perm, int size, double *sol) {
	int i;
	for (i = 0; i < size; i++) {
		if (perm[i] != i) {
			double temp = sol[i];
			sol[i] = sol[perm[i]];
			sol[perm[i]] = temp;
		}
	}
	for (i = 0; i < size; i++) {
		int j;
		for (j = 0; j < i; j++) {
			sol[i] -= luMatrix[i * size + j] * sol[j];
		}
	}
	for (i = size - 1; i >= 0; i--) {
		int j;
		for (j = size - 1; j > i; j--) {
			sol[i] -= luMatrix[i * size + j] * sol[j];
		}
		sol[i] /= luMatrix[i * size + i];
	}
}

void CILUPreconditioner::precondition(CSparseMatrix *matrix) {
// 	mprintf(GREEN, "%#.10g\n", matrix->m_val[13600000]);
	int *diagInd;
	try { diagInd = new int[matrix->m_size]; } catch(...) { diagInd = NULL; }
	if (diagInd == NULL) NewException((double)matrix->m_size * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	int i;
	for (i = 0; i < matrix->m_size; i++) {
		diagInd[i] = -1;
		int j;
		for (j = matrix->m_rowPtr[i]; j < matrix->m_rowPtr[i + 1]; j++) {
			if (i == matrix->m_colInd[j]) {
				diagInd[i] = j;
				break;
			}
		}
		if (diagInd[i] == -1) {
			printf("ILU preconditioner: Missing diagonal element in row %d!\n", i + 1);
			abort();
		}
	}
	for (i = 1; i < matrix->m_size; i++) {
		int k;
		for (k = matrix->m_rowPtr[i]; matrix->m_colInd[k] < i && k < matrix->m_rowPtr[i + 1]; k++) {
			matrix->m_val[k] /= matrix->m_val[diagInd[matrix->m_colInd[k]]];
			int j;
			for (j = k + 1; j < matrix->m_rowPtr[i + 1]; j++) {
				int l;
				for (l = matrix->m_rowPtr[matrix->m_colInd[k]]; l < matrix->m_rowPtr[matrix->m_colInd[k] + 1]; l++) {
					if (matrix->m_colInd[l] == matrix->m_colInd[j]) {
						matrix->m_val[j] -= matrix->m_val[k] * matrix->m_val[l];
						break;
					}
				}
			}
		}
	}
	delete[] diagInd;
}

void CILUPreconditioner::solve(CSparseMatrix *matrix, double *vect) {
	int i;
// 	mprintf(GREEN, "%#.10g\n", vector[1953124]);
	for (i = 0; i < matrix->m_size; i++) {
		int j;
		for (j = matrix->m_rowPtr[i]; matrix->m_colInd[j] < i && j < matrix->m_rowPtr[i + 1]; j++) {
			vect[i] -= matrix->m_val[j] * vect[matrix->m_colInd[j]];
// 			if (i == 1953124) {
// 				mprintf(GREEN, "%d %#.10g %#.10g\n", j, matrix->m_val[j], vector[matrix->m_colInd[j]]);
// 			}
		}
	}
// 	mprintf(GREEN, "%#.10g\n", vector[1953124]);
	for (i = matrix->m_size - 1; i >= 0; i--) {
		int j;
		for (j = matrix->m_rowPtr[i + 1] - 1; matrix->m_colInd[j] > i && j >= matrix->m_rowPtr[i]; j--) {
			vect[i] -= matrix->m_val[j] * vect[matrix->m_colInd[j]];
// 			if (i == 1953124) {
// 				mprintf(GREEN, "%d %#.10g %#.10g\n", j, matrix->m_val[j], vector[matrix->m_colInd[j]]);
// 			}
		}
		for (j = matrix->m_rowPtr[i]; j < matrix->m_rowPtr[i + 1]; j++) {
			if (matrix->m_colInd[j] == i) {
				vect[i] /= matrix->m_val[j];
				break;
			}
		}
	}
// 	mprintf(GREEN, "%#.10g\n", vector[1953124]);
}
