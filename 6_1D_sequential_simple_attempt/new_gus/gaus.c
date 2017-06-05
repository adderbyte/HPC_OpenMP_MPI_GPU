#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include "utils.h"
#include <time.h>



int main  (int argc, char ** argv)
 {
// ----------decalre variables to monitor iterations and time-----------------------------//
 long int i,j,row;
 long int N=-999;
 int rank , np;
 float sum ; // store the sum of computations
 double num_t, start_t, end_t=0;

//------------------------Initiate Communication world-----------------------------------//
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&np); //communication world
MPI_Comm_rank(MPI_COMM_WORLD,&rank);  //rank in the communication world

//-------------get the size of array from user. Ensure valid problem size. -----------------------------------------// 
if (argc!=2) {
      if(rank == 0) {    
       fprintf(stderr,"Error: please supply  parameters ....");
       exit(-1);}
   }


else{
     N =  atoi(argv[1]);
    }

if (N<= 0)
    {
    if (rank == 0 ){
       fprintf(stderr,"Error: Invalid Array size ....");
       exit( -1);}
    }
 
// --------------------------------decalare array  to be used ----------------//
 double A[N][N]; // store array
 double b[N]  ; // store equation constant. The y's of the equation.
 double  upper[N]; //for getting input from neighbors
 double x[N]; // the parameters to be calculated
 double local_b[N];// beclare local b to be used for computation locally
 double local_a[N][N]; //declare array for computation locally in processoir
 double  t_upper[N]; //local neighbor communication

//------------diveide  computational load among load processors-----------------------//
 row = N/np; // for  load balancing;

// fill the matrice
//-------------Initialise  MASTER with the arrays ------------------------------------------------------//

if(rank==0){
// use [] notation to access array buckets 
// (THIS IS THE PREFERED WAY TO DO IT)
for(i=0; i < N; i++) {
  x[i] = 0.0;
  b[i] = ((double)rand() / (double)RAND_MAX) ;
  upper[i] =0; 
  for( j =0 ; j<N;j++ ){ A[i][j]=( ((double)rand() / (double)RAND_MAX)  + (i+10)*20 + j) ;   }
}

}
//----------Make array A  diagonally dominat before scattering ------------------------------------------------------//

if(rank==0){
    // Make the array diagonally dominant
         for (int i=0;i<N;i++){
                 for ( int j=0;j<=N;j++){
                              if(i ==j)
                              {    A[i][j] =  A[i][j] +( j+20000) *2  + (i+1)*10 ; } 

   }
   }

}

//-----------------------Scatter the array -----------------------------------------------------------//

MPI_Scatter(b, row, MPI_DOUBLE, local_b, row, MPI_DOUBLE,0,MPI_COMM_WORLD);
//MPI_Scatter(b, row, MPI_DOUBLE, local_b, row, MPI_DOUBLE,0,MPI_COMM_WORLD);
MPI_Scatter(A,row*N,MPI_DOUBLE,local_a,row*N, MPI_DOUBLE,0,MPI_COMM_WORLD);


//--------------------Begin Gauss Seide ------------------------------------------------------------//
//=--------------------------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------//


int flag;
for(flag=0; flag<4 ;flag++) {
//--------------Start timer and start computation -----------------------------------------------------//
   
	start_t = MPI_Wtime();
	MPI_Bcast(x,N,MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
	if (rank > 0)
	{
		for (i=0;i<row;i++)
		{
//------------------------compute sum per processor for rank>0-------------------------------------------------//
			sum=0.0;
				for (j=row+rank;j<N;j++) // ensure you compute sum for the other values except my  values
			{
					sum += x[j] * local_a[i][j];
			     
           }
//----------------------------------subtract sum from the constant y of equation just as in the equation ------------------------------//
	  t_upper[i] = local_b[i] -sum; //update this unto the variable t_upper
         }
        }
//------------------if rank is less than zero, typically this is for master processor" zero" computation. -----------------------------//
//-----The MAster also performs some computation . This fact is explained in the report ------------------------------------//
	else{for (i=0; i<row;i++)
		{sum =0.0;
			for (j=i+1;j<N;j++){
				sum+=A[i][j]*x[j];
	
			}
		t_upper[i] = local_b[i] -sum;
	}
	}
//-------------Gather results into new array called upper ---------------------------------------------------------------------------------//
	MPI_Gather(t_upper,row,MPI_DOUBLE,upper,row,MPI_DOUBLE,0,MPI_COMM_WORLD);

//------------------------------ Rank zero completes the computation after gathering ----------------------------------------------------------------------//
	if (rank == 0 ){
            //   double  temp [N];
                 
		for(i=0; i<N ; i++){

			j=0;
			sum=0.0;
			while (j<i-1){

				sum+= A[i][j] * x[j];
				j++;
			}

		x[i]= (upper[i] - sum )/A[i][j];
	
               
	}
            
         	} 
}     
                





		end_t = MPI_Wtime();
		num_t = end_t - start_t;


	
     // num_t = end_t - start_t;



if(rank==0){

 printf("\n Execution time: %lf \n ", num_t    );
 for(int i = 0 ; i< N ;i++)
	printf("%lf \n ",x[i] );

}




 MPI_Finalize();

 return 0;


}











