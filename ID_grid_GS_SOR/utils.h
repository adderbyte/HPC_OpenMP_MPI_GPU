#define abs_tolerance  0.000001
#define w  1.3
void init2d ( double ** A,double *b,double *p   ,int N );
void  Gauss_seidel(double **A, double ** local_a,double *t_upper,double *local_b,double* x, int row,int N,int rank);
int  convergence (double*x,int flag ,double *previous,int N);
double update (double** A, double* upper ,double* previous ,double* x, int N );

