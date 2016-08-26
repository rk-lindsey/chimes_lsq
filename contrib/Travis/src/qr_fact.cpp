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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qr_fact.h"
#include "tools.h"

int QR_Test()
{
  int i, j;

  /* Variables for the matrix, the right hand side, the solution, and copies. */
  double **A, **B;
  double *b, *c, *d;
  double *x;

  /* The residual. */
  double r;

  /* Dimensions for our matrix and vectors. */
  int M = 6;
  int N = 3;
  
  /* Allocate memory of the matrices and vectors. */
  A = (double**)malloc(sizeof(double *)*M);
  B = (double**)malloc(sizeof(double *)*M);
  b = (double*)malloc(sizeof(double)*M);
  c = (double*)malloc(sizeof(double)*M);
  d = (double*)malloc(sizeof(double)*M);
  x = (double*)malloc(sizeof(double)*N);

  /* Use a 2D array for the matrices. */
  for(i = 0; i < M; i++)
    {
      A[i] = (double*)malloc(sizeof(double)*N);
      B[i] = (double*)malloc(sizeof(double)*N);
    }

  /* Assign the matrix and the right hand side. Notice the format. */
  A[0][0] = 1;   A[0][1] = 2;   A[0][2] = 3;    b[0] = 3; 
  A[1][0] = 4;   A[1][1] = 5;   A[1][2] = 6;    b[1] = 9; 
  A[2][0] = 7;   A[2][1] = 8;   A[2][2] = 9;    b[2] = 15; 
  A[3][0] = 10;  A[3][1] = 11;  A[3][2] = 12;   b[3] = 22; 
  A[4][0] = 13;  A[4][1] = 14;  A[4][2] = 15;   b[4] = 27; 
  A[5][0] = 16;  A[5][1] = 17;  A[5][2] = -5;   b[5] = 33;

  /* Copy the matrix A into B since the QR routine will overwrite A. */
  for(i = 0; i < M; i++)
    for(j = 0; j < N; j++)
     B[i][j] = A[i][j];  

  /* Copy the vector b into d since the QR routine will overwrite b. */
  for(i = 0; i < M; i++)
    d[i] = b[i];

  mprintf("\n");
  mprintf("--------------------------------------------\n");
  mprintf("QR Least-Squares Solver with Column-Pivoting\n");
  mprintf("--------------------------------------------\n");

  mprintf("Matrix A: \n");
  QR_display_mat(A, M, N, NULL);

  mprintf("Right-hand side b: \n");
  QR_display_vec(b, M, NULL);

  /* Solve the least squares problem x = A\b. */
  QR_least_squares(A, b, x, M, N);

  mprintf("Solution x: \n");
  QR_display_vec(x, N, NULL);

  /* Use the copies of the initial matrix and vector to get the right hand side
   * which corresponds to the least squares solution. */
  QR_mat_vec(B, M, N, x, c);

  mprintf("Matrix-vector product A*x: \n");
  QR_display_vec(c, M, NULL);

  /* Compute the 2-norm of the difference between the original right hand side
   * and the right hand side computed from the least squares solution. */
  r = 0.0;
  for(i = 0; i < M; i++)
    r += (c[i] - d[i])*(c[i] - d[i]);
  r = sqrt(r);

  mprintf("Least squares residual: r = %lf\n", r);
  mprintf("\n");

  /* Collect garbage. */
  for(i = 0; i < M; i++)
    {
      free(A[i]);
      free(B[i]);
    }

  /* Collect more garbage. */
  free(A);
  free(B);
  free(b);
  free(c);
  free(d);
  free(x);

  return 0;
}

void QR_swap_cols(int *p, int i, int j)
{
  int temp;
  temp = p[i];
  p[i] = p[j];
  p[j] = temp;
}

void QR_back_solve(double **mat, double *rhs, int rows,
		int cols, double *sol, int *p)
{
  int i, j, bottom;
  double sum;

  (void)cols;
  bottom = 0;

  /* Fill the solution with zeros initially. */
  for(i = 0; i < cols; i++)
    sol[i] = 0.0;

  /* Find the first non-zero row from the bottom and start solving from here. */
  for(i = rows - 1; i >= 0; i--)
    if(fabs(mat[i][p[cols - 1]]) > 1e-7)
      {
	bottom = i;
	break;
      }

  /* Standard back solving routine starting at the first non-zero diagonal. */
  for(i = bottom; i >= 0; i--)
    {
      sum = 0.0;

      for(j = cols - 1; j >= 0; j--)
	if(j > i)
	  sum += sol[p[j]]*mat[i][p[j]];
      
      if(mat[i][p[i]] > 1e-7)
	sol[p[i]] = (rhs[i] - sum)/mat[i][p[i]];
      else
	sol[p[i]] = 0.0;
    }
}

void QR_householder(double **mat, int rows, int cols,
		 int row_pos, int col_pos, double *result)
{
  int i;
  double norm;

  (void)cols;

  norm = 0;
  for(i = row_pos; i < rows; i++)
    norm += mat[i][col_pos]*mat[i][col_pos];

  if(norm == 0)
    return;

  norm = sqrt(norm);

  result[0] = (mat[row_pos][col_pos] - norm);

  for(i = 1; i < (rows - row_pos); i++)
    result[i] = mat[i+row_pos][col_pos];

  norm = 0;
  for(i = 0; i < (rows - row_pos); i++)
    norm += result[i]*result[i];

  if(norm == 0)
    return;

  norm = sqrt(norm);

  for(i = 0; i < (rows - row_pos); i++)
    result[i] *= (1.0/norm);
}

void QR_apply_householder(double **mat, double *rhs, int rows, int cols, double *house, int row_pos, int *p)
{
  int i, j, k, n;
  double sum;
  double **hhmat;
  double **mat_cpy;
  double *rhs_cpy;

  // Get the dimensions for the Q matrix.
  n = rows - row_pos;

  // Allocate memory.
  hhmat = (double**)malloc(sizeof(double *)*n);
  for(i = 0; i < n; i++)
    hhmat[i] = (double*)malloc(sizeof(double)*n);

  mat_cpy = (double**)malloc(sizeof(double *)*rows);
  for(i = 0; i < rows; i++)
    mat_cpy[i] = (double*)malloc(sizeof(double)*cols);

  rhs_cpy = (double*)malloc(sizeof(double )*rows);
  
  // Copy the matrix.
  for(i = 0; i < rows; i++)
    for(j = 0; j < cols; j++)
      mat_cpy[i][j] = mat[i][j];

  // Copy the right hand side.
  for(i = 0; i < rows; i++)
    rhs_cpy[i] = rhs[i];

  // Build the Q matrix from the Householder transform.
  for(j = 0; j < n; j++)
    for(i = 0; i < n; i++)
      if(i != j)
	hhmat[i][j] = -2.0*house[j]*house[i];
      else
	hhmat[i][j] = 1.0 - 2.0*house[j]*house[i];

  // Multiply by the Q matrix.
  for(k = 0; k < cols; k++)
    for(j = 0; j < n; j++)
      {
	sum = 0.0;
	for(i = 0; i < n; i++)
	  sum += hhmat[j][i]*mat_cpy[i + row_pos][p[k]];
	  
	mat[j + row_pos][p[k]] = sum;
      }

  // Multiply the rhs by the Q matrix.
  for(j = 0; j < n; j++)
    {
      sum = 0.0;
      for(i = 0; i < n; i++)
	sum += hhmat[i][j]*rhs_cpy[i + row_pos];

      rhs[j + row_pos] = sum;
    }

  // Collect garbage.
  for(i = 0; i < (rows - row_pos); i++)
    free(hhmat[i]);

  for(i = 0; i < rows; i++)
    free(mat_cpy[i]);

  free(hhmat);
  free(mat_cpy);
  free(rhs_cpy);
}

int QR_get_next_col(double **mat, int rows, int cols,
			  int row_pos, int *p)
{
  int i, j, max_loc;
  double *col_norms;
  double max;

  max_loc = -1;
  col_norms = (double*)malloc(sizeof(double)*cols);

  // Compute the norms of the sub columns.
  for(j = 0; j < cols; j++)
    {
      col_norms[j] = 0;

      for(i = row_pos; i < rows; i++)
	col_norms[j] += mat[i][p[j]]*mat[i][p[j]];
    }

  // Find the maximum location.
  max = 1e-7;
  for(i = 0; i < cols; i++)
    if(col_norms[i] > max)
      {
	max = col_norms[i];
	max_loc = i;
      }

  // Collect garbge and return.
  free(col_norms);  
  return max_loc;
}

/* The star of the show. A QR least-squares solving routine for x = A\b.
 *
 * First argument : The row-major matrix (2D array), A.
 * Second argument: The right-hand side vector, b.
 * Third argument : The solution vector, x.
 * Fourth argument: The number of rows in A.
 * Fifth argument : The number of columns in A.
 *
 * WARNING: This routine will overwrite the matrix A and the right-hand side
 * vector b. In other words, A*x = b is solved using QR least-squares with, 
 * column pivoting, but neither the A nor b are what you started with. However,
 * the solution x corresponds to the solution of both the modified and original
 * systems. Please be aware of this.
 */
void QR_least_squares(double **mat, double *rhs, double *sol, 
		      int rows, int cols)
{
  int i, max_loc;
  int *p;
  double *v;

  /* Allocate memory for index vector and Householder transform vector. */
  p = (int*)malloc(sizeof(int)*cols);
  v = (double*)malloc(sizeof(double)*rows);

  /* Initial permutation vector. */
  for(i = 0; i < cols; i++)
    p[i] = i;
  
  /* Apply rotators to make R and Q'*b */
  for(i = 0; i < cols; i++)
    {
      max_loc = QR_get_next_col(mat, rows, cols, i, p);
      if(max_loc >= 0)
	QR_swap_cols(p, i, max_loc);

      QR_householder(mat, rows, cols, i, p[i], v);
      QR_apply_householder(mat, rhs, rows, cols, v, i, p);
    }

  /* Back solve Rx = Q'*b */
  QR_back_solve(mat, rhs, rows, cols, sol, p);

  /* Collect garbage. */
  free(p);
  free(v);
}

/* A very simple matrix vector product routine. */
void QR_mat_vec(double **mat, int rows, int cols,
	     double *vec, double *rhs)
{
  int i, j;
  double sum;

  for(i =  0; i < rows; i++)
    {
      sum = 0.0;
      for(j = 0; j < cols; j++)
	sum += mat[i][j]*vec[j];

      rhs[i] = sum;
    }
}

/* Simple routine for displaying a matrix. */
void QR_display_mat(double **mat, int rows, int cols,
		 int *p)
{
  int i, j;

  for(i = 0; i < rows; i++)
    {
      for(j = 0; j < cols; j++)
	if(p != NULL)
	  mprintf("\t%-3.5lf ", mat[i][p[j]]);
	else
	  mprintf("\t%-3.5lf ", mat[i][j]);
      
      mprintf("\n");
    }

  mprintf("\n");
}

/* Simple routine for displaying a vector. */
void QR_display_vec(double *vec, int rows, int *p)
{
  int i;

  for(i = 0; i < rows; i++)
    if(p != NULL)
      mprintf("\t%-3.5lf\n", vec[p[i]]);
    else
      mprintf("\t%-3.5lf\n", vec[i]);

  mprintf("\n");
}
