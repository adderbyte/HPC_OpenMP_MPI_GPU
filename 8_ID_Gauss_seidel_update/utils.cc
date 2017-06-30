#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#include <time.h>


void init2d (double** A,double* b,double* p,int N ){

 // declare the  augmented array input here

 //
  // Linear memory allocation to make contiguous mem segment
 double * tempStore = new  double[N * N];
              	
   for(int i=0 ; i< N ;i++){
		A[i] = (tempStore + i * N);
		b[i] = ((double)rand() / (double)RAND_MAX) ;
		p[i] =0;}
  for (int y = 0; y < N; ++y){
	for (int x = 0; x < N; ++x){
				 A[y][x] =  ((double)rand() / (double)RAND_MAX)  + sin (y+x) ;
    
				}			
	
				} 

for (int i=0;i<N;i++)
	for ( int j=0;j<=N;j++)
	if(i ==j)
		{
		A[i][j] =  A[i][j] +( j * 0.5)  + (i*0.5) ;
		
		}

}



void  Gauss_seidel(double **A, double ** local_a,double *t_upper,double *local_b,double* x, int row,int N,int rank)
{
	double sum =0.0;
	if (rank > 0)
	{
	for (int i=0;i<row;i++){
			sum =0.0;
			for (int j=row+rank;j<N;j++) 
				{
				
				sum+= x[j] * local_a[i][j];
		
			
			
				}
			t_upper[i] = local_b[i] -sum; 

	
	

			    }
	
	}


	else {  int j=0;
		for (int i=0; i<row;i++)
					{
					sum= 0.0;
					//j=0;
					for ( j=i+1;j<N;j++)
						{ 
					     	sum+=A[i][j]*x[j];
						}
					t_upper[i] = local_b[i] -sum;
					}		






	     }
}





double update (double** A, double* upper , double* x, int N )
		{
		int j ;
		double sum =0.0;
		for(int i=0; i<N ; i++)
					{
				    	j=0;
					sum+= A[i][j] * x[j];
					while(j < i-1)					
					{  sum+= A[i][j] * x[j];
					 j++;    
					}
					x[i]= (upper[i] - sum )/A[i][j];
					}

		

		}

int  convergence (double*x,int flag ,double *previous,int N)
	{
	for (int i=0;i<=N;i++){
		if (abs(x[i]-previous[i])<abs_tolerance) 
		{
		flag++;


		}	
		


	}

	for (int j=0;j<=N;j++)
		{
		
		previous[j] = x[j];	

		}






	
	
	return flag;

	}


