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

#include "normalmode.h"

#include "acf.h"
#include "globalvar.h"
#include "tools.h"
#include "xfloatarray.h"
#include "xobarray.h"

#include <stdio.h>

static float g_freqStep;
static CxFloatArray *g_transMatrix;

static double findRotationAngle(int n, CxObArray *cc_matrix, int i, int j) {
	CxFloatArray *aii = (CxFloatArray *)cc_matrix->GetAt(i * n + i);
	CxFloatArray *aij = (CxFloatArray *)cc_matrix->GetAt(i * n + j);
	CxFloatArray *ajj = (CxFloatArray *)cc_matrix->GetAt(j * n + j);
	double int1 = 0.0, int2 = 0.0, int3 = 0.0;
	int k;
	for(k = 0; k < aij->GetSize(); k++) {
// 		if(fabsf(aij->GetAt(k)) > 1.0e7f)
			int1 += (double)aij->GetAt(k) * (double)aij->GetAt(k);
// 		if((fabsf(aij->GetAt(k)) > 1.0e7f) && (fabsf(aii->GetAt(k)-ajj->GetAt(k)) > 1.0e7f))
			int2 += (double)aij->GetAt(k) * ((double)aii->GetAt(k) - (double)ajj->GetAt(k));
// 		if(fabsf(aii->GetAt(k)-ajj->GetAt(k)) > 1.0e7f)
			int3 += ((double)aii->GetAt(k) - (double)ajj->GetAt(k)) * ((double)aii->GetAt(k) - (double)ajj->GetAt(k));
	}
	int1 *= (double)g_freqStep;
	int2 *= (double)g_freqStep;
	int3 *= (double)g_freqStep;
// 	mprintf(RED, "---\n(%d, %d) int1 = %30g int2 = %30g int3 = %30g\n", i, j, int1, int2, int3);
	if(fabs(int2) < 1.0) {
// 		mprintf(RED, "int2 < 1.0\n");
		double f = 4.0 * int1 - int3;
// 		mprintf(RED, "f = %30g\n", f);
		if((fabs(int1) < 1.0) && (fabs(int2) < 1.0) && (fabs(int3) < 1.0))
			return 0.0;
		if(0.5 * f > 0.0)
			return 1.0;
		if(-2.0 * f > 0.0)
			return 0.0;
	} else {
// 		mprintf(RED, "int2 >= 1.0\n");
		double f = (4.0 * int1 - int3) / int2;
		double r[2][2];
		double s[2][2];
		for(k = 0; k < 2; k++) {
			int l;
			for(l = 0; l < 2; l++) {
// 	 			double r = (-f + (k == 0 ? -1.0 : 1.0) * sqrt(16.0 + f * f) + (l == 0 ? -1.0 : 1.0) * sqrt(2.0) * sqrt(16.0 + f * f - (k == 0 ? -1.0 : 1.0) * f * sqrt(16.0 + f * f))) / 4.0;
// 	 			double s = -2.0 * int2 / pow(r * r + 1.0, 4.0) * (2.0 * r * (pow(r, 4.0) - 14.0 * r * r + 9.0) + f * (3.0 * pow(r, 4.0) - 8.0 * r * r + 1.0));
// 	 			if((fabs(r) <= 1.0) && (s > 0.0))
// 	 				return (float)r;
				r[k][l] = (-f + (k == 0 ? -1.0 : 1.0) * sqrt(16.0 + f * f) + (l == 0 ? -1.0 : 1.0) * sqrt(2.0) * sqrt(16.0 + f * f - (k == 0 ? -1.0 : 1.0) * f * sqrt(16.0 + f * f))) / 4.0;
				s[k][l] = -2.0 * int2 / pow(r[k][l] * r[k][l] + 1.0, 4.0) * (2.0 * r[k][l] * (pow(r[k][l], 4.0) - 14.0 * r[k][l] * r[k][l] + 9.0) + f * (3.0 * pow(r[k][l], 4.0) - 8.0 * r[k][l] * r[k][l] + 1.0));
// 				mprintf(RED, "t = %10g, s = %30g\n", r[k][l], s[k][l]);
			}
		}
		for(k = 0; k < 2; k++) {
			int l;
			for(l = 0; l < 2; l++) {
				if((fabs(r[k][l]) <= 1.0) && (s[k][l] > 0.0))
					return r[k][l];
			}
		}
	}
// 	double f = (4.0 * int1 - int3) / int2;
// 	double r[2][2];
// 	double s[2][2];
// 	for(k = 0; k < 2; k++) {
// 		int l;
// 		for(l = 0; l < 2; l++) {
// 			double r = (-f + (k == 0 ? -1.0 : 1.0) * sqrt(16.0 + f * f) + (l == 0 ? -1.0 : 1.0) * sqrt(2.0) * sqrt(16.0 + f * f - (k == 0 ? -1.0 : 1.0) * f * sqrt(16.0 + f * f))) / 4.0;
// 			double s = -2.0 * int2 / pow(r * r + 1.0, 4.0) * (2.0 * r * (pow(r, 4.0) - 14.0 * r * r + 9.0) + f * (3.0 * pow(r, 4.0) - 8.0 * r * r + 1.0));
// 			if((fabs(r) <= 1.0) && (s >= 0.0))
// 				return (float)r;
// 			r[k][l] = (-f + (k == 0 ? -1.0 : 1.0) * sqrt(16.0 + f * f) + (l == 0 ? -1.0 : 1.0) * sqrt(2.0) * sqrt(16.0 + f * f - (k == 0 ? -1.0 : 1.0) * f * sqrt(16.0 + f * f))) / 4.0;
// 			s[k][l] = -2.0 * int2 / pow(r[k][l] * r[k][l] + 1.0, 4.0) * (2.0 * r[k][l] * (pow(r[k][l], 4.0) - 14.0 * r[k][l] * r[k][l] + 9.0) + f * (3.0 * pow(r[k][l], 4.0) - 8.0 * r[k][l] * r[k][l] + 1.0));
// 			mprintf(RED, "r = %10.6f, s = %30.2f\n", r[k][l], s[k][l]);
// 		}
// 	}
// 	for(k = 0; k < 2; k++) {
// 		int l;
// 		for(l = 0; l < 2; l++) {
// 			if((fabs(r[k][l]) <= 1.0f) && (s[k][l] >= 0.0f))
// 				return r[k][l];
// 		}
// 	}
	
	mprintf(RED, "No suitable rotation angle found. Aborting.\n");
	return 0.1;
// 	abort();
}

static double offDiagonalNorm(int n, CxObArray *cc_matrix) {
	double norm = 0.0;
	int i, j, k;
// 	mprintf("@");
	for(i = 0; i < n; i++) {
		for(j = i + 1; j < n; j++) {
			double integral = 0.0;
			CxFloatArray *fa = (CxFloatArray *)cc_matrix->GetAt(i * n + j);
			for(k = 0; k < fa->GetSize(); k++) {
// 				if(fabsf(fa->GetAt(k)) > 1.0e7f)
					integral += fa->GetAt(k) * fa->GetAt(k);
			}
			integral *= g_freqStep;
// 			mprintf("%14g ", integral);
			norm += integral;
		}
// 		mprintf("\n@");
	}
// 	mprintf("\n");
	return sqrt(norm);
}

void normalModeAnalysis(int n, CxObArray *cc_matrix) {
	mprintf("    Fourier transforming cross correlation matrix...\n");
	mprintf(WHITE, "      [");
	int i, j, k;

	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			CxFloatArray *temp;
			try { temp = new CxFloatArray("normalModeAnalysis():temp"); } catch(...) { temp = NULL; }
			if(temp == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			temp->CopyFrom((CxFloatArray *)cc_matrix->GetAt(i * n + j));
			if(g_pGlobalVACF->m_bWindowFunction) {
				for(k = 0; k < temp->GetSize(); k++) {
					temp->GetAt(k) *= powf(cosf((float)k / temp->GetSize() / 2.0f * Pi), 2.0f);
				}
			}
			if(g_pGlobalVACF->m_iZeroPadding > 0) {
				for(k = 0; k < g_pGlobalVACF->m_iZeroPadding; k++) {
					temp->Add(0.0f);
				}
			}
			if(g_pGlobalVACF->m_iMirror != 0) {
				int oldSize = temp->GetSize();
				temp->SetSize(2 * oldSize);
				for(k = 1; k < oldSize; k++) {
					temp->GetAt(oldSize + k) = temp->GetAt(oldSize - k);
				}
			}
			CFFT *fft;
			try { fft = new CFFT(); } catch(...) { fft = NULL; }
			if(fft == NULL) NewException((double)sizeof(fft), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			fft->PrepareFFT_C2C(temp->GetSize());
			for(k = 0; k < temp->GetSize(); k++) {
				fft->m_pInput[2 * k] = temp->GetAt(k);
				fft->m_pInput[2 * k + 1] = 0.0f;
			}
			fft->DoFFT();
			((CxFloatArray *)cc_matrix->GetAt(i * n + j))->SetSize(temp->GetSize());
			for(k = 0; k < temp->GetSize(); k++) {
				((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k) = fft->m_pOutput[2 * k];
			}
			delete fft;
			delete temp;
			
// 			CxFloatArray *temp;
// 			try { temp = new CxFloatArray(); } catch(...) { temp = NULL; }
// 			if(temp == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			temp->CopyFrom((CxFloatArray *)cc_matrix->GetAt(i * n + j));
// 			if(g_pGlobalVACF->m_bWindowFunction) {
// 				for(k = 0; k < temp->GetSize(); k++) {
// 					temp->GetAt(k) *= powf(cosf((float)k / temp->GetSize() / 2.0f * Pi), 2.0f);
// 				}
// 			}
// 			if(g_pGlobalVACF->m_iZeroPadding > 0) {
// 				for(k = 0; k < g_pGlobalVACF->m_iZeroPadding; k++) {
// 					temp->Add(0.0f);
// 				}
// 			}
// 			if(g_pGlobalVACF->m_iMirror != 0) {
// 				int oldSize = temp->GetSize();
// 				temp->SetSize(2 * oldSize);
// 				for(k = 1; k < oldSize; k++) {
// 					temp->GetAt(oldSize + k) = temp->GetAt(oldSize - k);
// 				}
// 			}
// 			CFFT *fft;
// 			try { fft = new CFFT(); } catch(...) { fft = NULL; }
// 			if(fft == NULL) NewException((double)sizeof(fft), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			fft->PrepareFFT_C2C(temp->GetSize());
// 			for(k = 0; k < temp->GetSize(); k++) {
// 				fft->m_pInput[2 * k] = temp->GetAt(k);
// 				fft->m_pInput[2 * k + 1] = 0.0f;
// 			}
// 			fft->DoFFT();
// 			((CxFloatArray *)cc_matrix->GetAt(i * n + j))->SetSize(temp->GetSize());
// 			for(k = 0; k < temp->GetSize(); k++) {
// 				((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k) = fft->m_pOutput[2 * k];
// 			}
// 			delete fft;
// 			delete temp;
// 			CACF *tempACF;
// 			try { tempACF = new CACF(); } catch(...) { tempACF = NULL; }
// 			if(tempACF == NULL) NewException((double)sizeof(CACF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			tempACF->m_iSize = ((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetSize();
// 			tempACF->m_bSpectrum = true;
// 			tempACF->m_bWindowFunction = g_pGlobalVACF->m_bWindowFunction;
// 			tempACF->m_fSpecWaveNumber = g_pGlobalVACF->m_fSpecWaveNumber;
// 			tempACF->m_iMirror = g_pGlobalVACF->m_iMirror;
// 			tempACF->m_iZeroPadding = g_pGlobalVACF->m_iZeroPadding;
// 			tempACF->Create();
// 			for(k = 0; k < tempACF->m_iSize; k++) {
// 				tempACF->m_pData[k] = ((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k);
// 			}
// 			if(tempACF->m_iMirror != 0)
// 				tempACF->Mirror(tempACF->m_iMirror);
// 			if(tempACF->m_bWindowFunction)
// 				tempACF->Window();
// 			CFFT *fft;
// 			try { fft = new CFFT(); } catch(...) { fft = NULL; }
// 			if(fft == NULL) NewException((double)sizeof(fft), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 			fft->PrepareFFT_C2C(tempACF->m_iSize + tempACF->m_iZeroPadding);
// 			tempACF->Transform(fft);
// 			delete fft;
// 			((CxFloatArray *)cc_matrix->GetAt(i * n + j))->SetSize(tempACF->m_pSpectrum->m_iSize);
// 			for(k = 0; k < tempACF->m_pSpectrum->m_iSize; k++) {
// 				((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k) = tempACF->m_pSpectrum->m_pComplex[2*k];
// 			}
// 			delete tempACF;
		}
		mprintf(WHITE, "#");
	}
	mprintf(WHITE, "]\n");
	g_freqStep = 1e15f / 299792458.0f / 100.0f / g_fTimestepLength / g_iStride / ((CxFloatArray *)cc_matrix->GetAt(0))->GetSize();
	
	mprintf("    Allocating transformation matrix...\n");
	try { g_transMatrix = new CxFloatArray("normalModeAnalysis():g_transMatrix"); } catch(...) { g_transMatrix = NULL; }
	if(g_transMatrix == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
// 	g_transMatrix->SetSize(n * n);
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			g_transMatrix->Add((i == j) ? 1.0f : 0.0f);
		}
	}
	
	mprintf("    Saving transformation matrix\n");
	FILE *matrixi_file = OpenFileWrite("transformation_matrix_i.dat", false);
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			fprintf(matrixi_file, " %14g", g_transMatrix->GetAt(i * n + j));
		}
		fprintf(matrixi_file, "\n");
	}
	
	CxFloatArray *speciSum;
	try { speciSum = new CxFloatArray("normalModeAnalysis():speciSum"); } catch(...) { speciSum = NULL; }
	if(speciSum == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	for(i = 0; i < ((CxFloatArray *)cc_matrix->GetAt(0))->GetSize(); i++)
		speciSum->Add(0.0f);
	for(i = 0; i < n; i++) {
		for(k = 0; k < ((CxFloatArray *)cc_matrix->GetAt(i * n + i))->GetSize(); k++) {
			speciSum->GetAt(k) += ((CxFloatArray *)cc_matrix->GetAt(i * n + i))->GetAt(k);
		}
	}
	FILE *speci_file = OpenFileWrite("normalmode_sum_initial.dat", false);
	for(k = 0;  k * g_freqStep < g_pGlobalVACF->m_fSpecWaveNumber; k++) {
		fprintf(speci_file, "%14G %20G\n", k * g_freqStep, speciSum->GetAt(k));
	}
	fclose(speci_file);
	delete(speciSum);
	
	for(i = 0; i < n; i++) {
		for(j = i; j < n; j++) {
			char name[128];

#ifdef TARGET_WINDOWS
			_snprintf(name, 128, "speci_%d_%d.dat", i, j);
#elif defined TARGET_LINUX
			snprintf(name, 128, "speci_%d_%d.dat", i, j);
#else
			sprintf(name, "speci_%d_%d.dat", i, j);
#endif

			FILE *speci_file = fopen(name, "w");
			for(k = 0; k * g_freqStep < g_pGlobalVACF->m_fSpecWaveNumber; k++) {
				fprintf(speci_file, "%14G %20G\n", k * g_freqStep, ((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k));
			}
			fclose(speci_file);
		}
	}
	
	mprintf("    Minimizing cross correlations...\n");
	double norm = offDiagonalNorm(n, cc_matrix);
	mprintf("      %4d %20g\n", 0, norm);
	double change = norm;
	int count = 0;
	CxFloatArray *tempMatrix;
	try { tempMatrix = new CxFloatArray("normalModeAnalysis():tempMatrix"); } catch(...) { tempMatrix = NULL; }
	if(tempMatrix == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	for(i = 0; i < n; i++) {
		tempMatrix->Add(0.0f);
	}
	while(fabs(change) > 1.0e-3 * norm) {
		count++;
		for(i = 0; i < n; i++) {
			for(j = i + 1; j < n; j++) {
				double t = findRotationAngle(n, cc_matrix, i, j);
				double c = 1.0 / sqrt(t * t + 1.0);
				double s = t * c;
// 				mprintf(RED, "(%d, %d) t = %10.6f c = %10.6f s = %10.6f\n", i, j, t, c, s);
				for(k = 0; k < ((CxFloatArray *)cc_matrix->GetAt(0))->GetSize(); k++) {
					float temp[3];
					int l;
					for(l = 0; l < j; l++) {
						if(l < i) {
							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) - s * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
							((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) = temp[0];
							((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k) = temp[1];
						} else if(l > i) {
							tempMatrix->GetAt(l) = ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
							((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k) = temp[1];
						}
					}
					for(l = i + 1; l < n; l++) {
						if(l < j) {
							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) - s * tempMatrix->GetAt(l);
							((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) = temp[0];
						} else if(l > j) {
							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) - s * ((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k);
							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k);
							((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) = temp[0];
							((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k) = temp[1];
						}
					}
					temp[0] = c * c * ((CxFloatArray *)cc_matrix->GetAt(i * n + i))->GetAt(k) + s * s * ((CxFloatArray *)cc_matrix->GetAt(j * n + j))->GetAt(k) - 2 * c * s * ((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k);
					temp[1] = (c * c - s * s) * ((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k) + c * s * (((CxFloatArray *)cc_matrix->GetAt(i * n + i))->GetAt(k) - ((CxFloatArray *)cc_matrix->GetAt(j * n + j))->GetAt(k));
					temp[2] = s * s * ((CxFloatArray *)cc_matrix->GetAt(i * n + i))->GetAt(k) + c * c * ((CxFloatArray *)cc_matrix->GetAt(j * n + j))->GetAt(k) + 2 * c * s * ((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k);
					((CxFloatArray *)cc_matrix->GetAt(i * n + i))->GetAt(k) = temp[0];
					((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k) = temp[1];
					((CxFloatArray *)cc_matrix->GetAt(j * n + j))->GetAt(k) = temp[2];
					
//					for(l = 0; l < n; l++) {
//						if(l < i) {
//							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) - s * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
//							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
//							((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) = temp[0];
//							((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k) = temp[1];
//						} else if(l < j) {
//							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) - s * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
//							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
//							((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) = temp[0];
//							((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k) = temp[1];
//						} else {
//// 					for(l = 0; l < n; l++) {
//							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) - s * ((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k);
//							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k);
//							((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) = temp[0];
//							((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k) = temp[1];
//						}
//					}
//					for(l = n - 1; l >= 0; l--) {
//						if(l > j) {
//							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) - s * ((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k);
//							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(j * n + l))->GetAt(k);
//							((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) = temp[0];
//							((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k) = temp[1];
//						} else if(l > i) {
//							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) - s * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
//							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(i * n + l))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
//							((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) = temp[0];
//							((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k) = temp[1];
//						} else {
//// 					for(l = 0; l < n; l++) {
//							temp[0] = c * ((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) - s * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
//							temp[1] = s * ((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) + c * ((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k);
////							((CxFloatArray *)cc_matrix->GetAt(l * n + i))->GetAt(k) = temp[0];
//							((CxFloatArray *)cc_matrix->GetAt(l * n + j))->GetAt(k) = temp[1];
//						}
//					}
				}
				float temp[2];
				int l;
				for(l = 0; l < n; l++) {
					temp[0] = c * g_transMatrix->GetAt(i * n + l) - s * g_transMatrix->GetAt(j * n + l);
					temp[1] = s * g_transMatrix->GetAt(i * n + l) + c * g_transMatrix->GetAt(j * n + l);
					g_transMatrix->GetAt(i * n + l) = temp[0];
					g_transMatrix->GetAt(j * n + l) = temp[1];
				}
// 				mprintf("@@     (%d, %d) %20g\n", i, j, offDiagonalNorm(n, cc_matrix));
			}
		}
		double newNorm = offDiagonalNorm(n, cc_matrix);
		change = newNorm - norm;
		norm = newNorm;
		mprintf("      %4d %20g\n", count, norm);
	}
	delete tempMatrix;
	
// 	float temp[9];
// 	for(k = 0; k < ((CxFloatArray *)cc_matrix->GetAt(0))->GetSize(); k++) {
// 		temp[0] = -0.0259f * ((CxFloatArray *)cc_matrix->GetAt(0))->GetAt(k) - 0.978f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.209f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k);
// 		temp[1] = -0.0259f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) - 0.978f * ((CxFloatArray *)cc_matrix->GetAt(30))->GetAt(k) + 0.209f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k);
// 		temp[2] = -0.0259f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k) - 0.978f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k) + 0.209f * ((CxFloatArray *)cc_matrix->GetAt(60))->GetAt(k);
// 		temp[3] = 0.986f * ((CxFloatArray *)cc_matrix->GetAt(0))->GetAt(k) + 0.009f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.164f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k);
// 		temp[4] = 0.986f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.009f * ((CxFloatArray *)cc_matrix->GetAt(30))->GetAt(k) + 0.164f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k);
// 		temp[5] = 0.986f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k) + 0.009f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k) + 0.164f * ((CxFloatArray *)cc_matrix->GetAt(60))->GetAt(k);
// 		temp[6] = -0.163f * ((CxFloatArray *)cc_matrix->GetAt(0))->GetAt(k) + 0.210f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.964f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k);
// 		temp[7] = -0.163f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.210f * ((CxFloatArray *)cc_matrix->GetAt(30))->GetAt(k) + 0.964f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k);
// 		temp[8] = -0.163f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k) + 0.210f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k) + 0.964f * ((CxFloatArray *)cc_matrix->GetAt(60))->GetAt(k);
// 		((CxFloatArray *)cc_matrix->GetAt(0))->GetAt(k) = temp[0];
// 		((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) = temp[1];
// 		((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k) = temp[2];
// 		((CxFloatArray *)cc_matrix->GetAt(12))->GetAt(k) = temp[4];
// 		((CxFloatArray *)cc_matrix->GetAt(15))->GetAt(k) = temp[5];
// 		((CxFloatArray *)cc_matrix->GetAt(24))->GetAt(k) = temp[8];
// 	}
// 	for(k = 0; k < ((CxFloatArray *)cc_matrix->GetAt(0))->GetSize(); k++) {
// 		temp[0] = -0.0259f * ((CxFloatArray *)cc_matrix->GetAt(0))->GetAt(k) - 0.978f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.209f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k);
// 		temp[1] = -0.0259f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) - 0.978f * ((CxFloatArray *)cc_matrix->GetAt(30))->GetAt(k) + 0.209f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k);
// 		temp[2] = -0.0259f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k) - 0.978f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k) + 0.209f * ((CxFloatArray *)cc_matrix->GetAt(60))->GetAt(k);
// 		temp[3] = 0.986f * ((CxFloatArray *)cc_matrix->GetAt(0))->GetAt(k) + 0.009f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.164f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k);
// 		temp[4] = 0.986f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.009f * ((CxFloatArray *)cc_matrix->GetAt(30))->GetAt(k) + 0.164f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k);
// 		temp[5] = 0.986f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k) + 0.009f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k) + 0.164f * ((CxFloatArray *)cc_matrix->GetAt(60))->GetAt(k);
// 		temp[6] = -0.163f * ((CxFloatArray *)cc_matrix->GetAt(0))->GetAt(k) + 0.210f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.964f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k);
// 		temp[7] = -0.163f * ((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) + 0.210f * ((CxFloatArray *)cc_matrix->GetAt(30))->GetAt(k) + 0.964f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k);
// 		temp[8] = -0.163f * ((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k) + 0.210f * ((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k) + 0.964f * ((CxFloatArray *)cc_matrix->GetAt(60))->GetAt(k);
// 		((CxFloatArray *)cc_matrix->GetAt(0))->GetAt(k) = temp[0];
// 		((CxFloatArray *)cc_matrix->GetAt(3))->GetAt(k) = temp[1];
// 		((CxFloatArray *)cc_matrix->GetAt(6))->GetAt(k) = temp[2];
// 		((CxFloatArray *)cc_matrix->GetAt(30))->GetAt(k) = temp[4];
// 		((CxFloatArray *)cc_matrix->GetAt(33))->GetAt(k) = temp[5];
// 		((CxFloatArray *)cc_matrix->GetAt(60))->GetAt(k) = temp[8];
// 	}
// 	for(k = 0; k < 4096; k++) {
// 		int l;
// 		for (l=0;l<81;l++)
// 		((CxFloatArray *)cc_matrix->GetAt(l))->GetAt(k) = 1000.0f;
// 	}
// 	mprintf(RED, "Old: %20g, New: %20g\n", norm, offDiagonalNorm(n, cc_matrix));
	
	CxFloatArray *specSum;
	try { specSum = new CxFloatArray("normalModeAnalysis():specSum"); } catch(...) { specSum = NULL; }
	if(specSum == NULL) NewException((double)sizeof(CxFloatArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	for(i = 0; i < ((CxFloatArray *)cc_matrix->GetAt(0))->GetSize(); i++)
		specSum->Add(0.0f);
	
	mprintf("    Saving normal mode spectra...\n");
	for(i = 0; i < n; i++) {
		char name[128];

#ifdef TARGET_WINDOWS
		_snprintf(name, 128, "normalmode_%d.dat", i);
#elif defined TARGET_LINUX
		snprintf(name, 128, "normalmode_%d.dat", i);
#else
		sprintf(name, "normalmode_%d.dat", i);
#endif

		FILE *spec_file = OpenFileWrite(name, false);
		for(k = 0; k * g_freqStep < g_pGlobalVACF->m_fSpecWaveNumber; k++) {
			fprintf(spec_file, "%14G %20G\n", k * g_freqStep, ((CxFloatArray *)cc_matrix->GetAt(i * n + i))->GetAt(k));
			specSum->GetAt(k) += ((CxFloatArray *)cc_matrix->GetAt(i * n + i))->GetAt(k);
		}
		fclose(spec_file);
	}
	
	mprintf("    Saving sum of normal mode spectra...\n");
	FILE *spec_file = OpenFileWrite("normalmode_sum.dat", false);
	for(k = 0; k * g_freqStep < g_pGlobalVACF->m_fSpecWaveNumber; k++) {
		fprintf(spec_file, "%14G %20G\n", k * g_freqStep, specSum->GetAt(k));
	}
	fclose(spec_file);
	
	delete specSum;
	
	for(i = 0; i < n; i++) {
		for(j = i + 1; j < n; j++) {
			char name[128];

#ifdef TARGET_WINDOWS
			_snprintf(name, 128, "spec_%d_%d.dat", i, j);
#elif defined TARGET_LINUX
			snprintf(name, 128, "spec_%d_%d.dat", i, j);
#else
			sprintf(name, "spec_%d_%d.dat", i, j);
#endif

			FILE *spec_file = fopen(name, "w");
			for(k = 0; k * g_freqStep < g_pGlobalVACF->m_fSpecWaveNumber; k++) {
				fprintf(spec_file, "%14G %20G\n", k * g_freqStep, ((CxFloatArray *)cc_matrix->GetAt(i * n + j))->GetAt(k));
			}
			fclose(spec_file);
		}
	}
	
	mprintf("    Saving transformation matrix\n");
	FILE *matrix_file = OpenFileWrite("transformation_matrix.dat", false);
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
			fprintf(matrix_file, " %14g", g_transMatrix->GetAt(i * n + j));
		}
		fprintf(matrix_file, "\n");
	}
	fclose(matrix_file);
	
	delete g_transMatrix;
}
