/* QR Least-Squares Solver with Column Pivoting
 * 
 * Note: This code will overwrite the input matrix! Create a copy if needed.
 *
 * Implemented in ANSI C. To compile with gcc: gcc -o qr_fact qr_fact.c -lm
 *
 * Author: Michael Mazack, <mazack @ yahoo . com>
 *
 * Date: April 27th, 2010
 *
 * License: Public Domain. Redistribution and modification without 
 * restriction is granted. If you find this code helpful, please let me know.
 */

#ifndef QR_FACT
#define QR_FACT


/* For updating the permuatation vector in (virtual) column swaps. */
void QR_swap_cols(int *p, int i, int j);

/* Backsolving of a trianglular system. */
void QR_back_solve(double **mat, double *rhs, int rows,
		int cols, double *sol, int *p);

/* Apply a Householder transform to the matrix at a given spot. */
void QR_householder(double **mat, int rows, int cols,
		 int row_pos, int col_pos, double *result);

/* Routine for applying the Householder transform to the matrix and the 
 * right hand side. */
void QR_apply_householder(double **mat, double *rhs, int rows, 
		       int cols, double *house, int row_pos,
		       int *p);

/* Get the column with the largest sub-norm starting from i = p[j] = row_pos. */
int QR_get_next_col(double **mat, int rows, int cols,
			  int row_pos, int *p);

/* Solve the least squares problem, sol = mat\rhs . */
void QR_least_squares(double **mat, double *rhs, double *sol, 
		      int rows, int cols);

/* A simple matrix-vector product routine. */
void QR_mat_vec(double **mat, int rows, int cols,
	     double *vec, double *rhs);

/* Routine for displaying a matrix. */
void QR_display_mat(double **mat, int rows, int cols,
		 int *p);

/* Routine for displaying a vector. */
void QR_display_vec(double *vec, int rows, int *p);

int QR_Test();


#endif

