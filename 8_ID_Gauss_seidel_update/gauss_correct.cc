#include<iostream>
#include<iomanip>
#include<math.h>
#include<random>
#include <sstream>
#include <chrono>
#include "mpi.h"
#include "utils.h"

typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> second;


using namespace std;
int main( int argc , char **argv )
{
    //check argument input

     if (argc != 2)
    {
        std::cerr << "You have to specify the  array dimensions" << std::endl;
        return -1;
    }


    int N;

    N = std::stoi(argv[1]) ;

    // check input size array is not equals to zero

    if (N <= 0)
    {
        std::cerr << "Invalid dimension x" << std::endl;
        return -1;
    }
int np,rank;
MPI_Init(&argc,&argv);
MPI_Comm_size(MPI_COMM_WORLD,&np); //communication world
MPI_Comm_rank(MPI_COMM_WORLD,&rank);  //ran


double **A;
double *b ;
double *upper;
double *x;
double **local_a;
double *local_b;
double *t_upper;
double *previous; // previous value used in test for convergence



auto start = clk::now();
//set output precision
cout.precision(4);
cout.setf(ios::fixed);


// ---initialise array from rank 0 which serves as master ----//
if (rank==0){
 //initialise array to be used
	A = new double*[N];
	upper = new double[N];
	b= new double[N];  
	init2d(A, b, upper,N);
}



////------------diveide  computational load among load processors--//
//int row; // send count for b
int  row1; // for load balancing. send count for array A.

row1 = N/np;


//-----Scatter matrix -------//
local_b = new double[N];
local_a = new double*[N];
//allocate memory to local storage
double * tempStore = new  double[N * N];
for (int i = 0; i < N; ++i){
	local_a[i] = (tempStore + i * N);
}

 x      = new double[N];
t_upper = new double[N];

// allocate pointer
double * ptr =(double*) malloc(sizeof(double));

previous = new double[N];
if (rank == 0) {
	delete ptr;
//	delete ptr2;
	ptr = &A[0][0];
	

}

MPI_Scatter(&b[0], row1, MPI_DOUBLE, &local_b[0], row1, MPI_DOUBLE,0,MPI_COMM_WORLD);

MPI_Scatter(ptr,row1*(N),MPI_DOUBLE,&local_a[0][0],row1*(N), MPI_DOUBLE,0,MPI_COMM_WORLD);

/*uncomment to see output of local a
if (rank==1){
for(int i =0 ; i<N;i++)
	for(int j;j<N;j++)
	 	{cout<<local_a[i][j]<<endl;
		cout<<endl;}
}
*/

int flag =0; // for convergence;
//int * ptrflag = &flag;
for(int i ; i< N ;i++)
{
previous[i] = sin ( i + 0.04);
 
}


do
 {
	
//	MPI_Bcast(&x[0],N,MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&flag,1,MPI_INT, 0, MPI_COMM_WORLD);
	Gauss_seidel(A, local_a,t_upper ,local_b,x, row1, N,rank);
	MPI_Gather(&t_upper[0],row1,MPI_DOUBLE,&upper[0],row1,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if( rank ==0){  update ( A, upper , x,  N );	}
	MPI_Bcast(&x[0],N,MPI_DOUBLE, 0, MPI_COMM_WORLD);
 	flag=convergence(x,flag,previous,N);
 //  	flag++;
} while(flag<N);
/*
if(rank==0)
{
cout<<"\n--------------- Inside Gauss Siedel Function-------------------------------------------------------";
for(int i ; i<N;i++)
    cout<<b[i]<<endl;
//   fflush(stdout);
}
*/


//-----Get final time -----------//
if(rank ==0){
auto end = clk::now();
second time = end - start;
cout<<"Time taken: "<< time.count()<<endl;
cout <<"End."<< endl;
}
delete[] local_a[0];
delete[] local_a;
local_a=nullptr;
delete local_b;
delete previous;
delete x;
delete t_upper;
//delete ptrflag;
//-----clear memory ---------------//
if(rank==0)
{
	//free (ptr);
	delete[] A[0];
	delete[] A;
	A= nullptr;
	delete upper;
	delete b;}

MPI_Finalize();  

return 0;
}
