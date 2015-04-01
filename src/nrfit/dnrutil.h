double *vector(int,int);
double **matrix(int,int,int,int);
double **convert_matrix();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
double **submatrix();
void free_vector(double*,int,int);
void free_dvector();
void free_ivector();
void free_matrix(double**,int,int,int,int);
void free_dmatrix();
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror();