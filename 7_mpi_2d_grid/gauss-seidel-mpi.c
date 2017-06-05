#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include "utils.h"

void Gauss_Computation(double ** last_val, double ** present_val, int X_min, int X_max, int Y_min, int Y_max, double relaxation_parameter) {
	int i,j;
	// for more complete gauss-seidel here we can use version from the ID case-- kindly cross check
	for (i=X_min;i<X_max;i++)
		for (j=Y_min;j<Y_max;j++)
			present_val[i][j]=relaxation_parameter*present_val[i][j]+ (1-relaxation_parameter)*(last_val[i-1][j]) ;
}

int main(int argc, char ** argv) {
	//------------define  parametrs needed in the computation for loop and timing -------------//
    int rank,size;
    int global[2],local[2];
    int padded_all_com[2];
    int grid[2];
    int padded[2] = {0,0};
    int i,j,t;
    int convergence_check=0, converged=0;
    MPI_Datatype place_holder; //  place holder for MPI_data type

    double relaxation_parameter; // for convergence
    struct timeval tts,ttf,tcs,tcf;
    double ttotal=0,tcomp=0,total_time,execution_time;


 // -------define array for use in computation ------------------------//

    double ** U = malloc(1), ** present_val, ** last_val, ** swap; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //---- Get input from users----//

    if (argc!=5) {
        fprintf(stderr,"Supply 4 parameters______please");
        exit(-1);
    }
    else {
        global[0]=atoi(argv[1]);
        global[1]=atoi(argv[2]);
        grid[0]=atoi(argv[3]);
        grid[1]=atoi(argv[4]);
    }

    
    //---------------Communication World ---------------------------------------//
    MPI_Comm CART_COMM;
	int periods[2]={0,0};
    int rank_grid[2];
    MPI_Cart_create(MPI_COMM_WORLD,2,grid,periods,0,&CART_COMM);    //communicator creation
	MPI_Cart_coords(CART_COMM,rank,2,rank_grid);	                //

    // ---------------------initialise   array -------------------------//
    for (i=0;i<2;i++) {
        if (global[i]%grid[i]==0) {
            local[i]=global[i]/grid[i];
            padded_all_com[i]=global[i];
        }
        else {
            local[i]=(global[i]/grid[i])+1;
            padded_all_com[i]=local[i]*grid[i];
	    padded[i] = 1;
        }
    }


//------------parameter for convergence check--------------------------------------//


    relaxation_parameter=2.0/(1+sin(3.14/global[0]));

//---------allocate and iitialise U ----------------------------------------//
    if (rank==0) {
	free(U);
        U=allocate2d(padded_all_com[0],padded_all_com[1]);   
        init2d(U,global[0],global[1]);
    }

//-------------------------------------------initialise  blocks of array for local comptation-----------//
    present_val = allocate2d(local[0] + 2, local[1] + 2);
    last_val = allocate2d(local[0] + 2, local[1] + 2);
//---------------------------Initialise MPI data type ----------------------------//
//----------facilitates transpfer of data ---------------------------------------//
    MPI_Datatype mpi_block_all;
    MPI_Type_vector(local[0],local[1],padded_all_com[1],MPI_DOUBLE,&place_holder);
    MPI_Type_create_resized(place_holder,0,sizeof(double),&mpi_block_all);
    MPI_Type_commit(&mpi_block_all);
    MPI_Datatype local_block;
    MPI_Type_vector(local[0],local[1],local[1]+2,MPI_DOUBLE,&place_holder);
    MPI_Type_create_resized(place_holder,0,sizeof(double),&local_block);
    MPI_Type_commit(&local_block);

//----Master Rank zero scatters the matrix  -------------------------------------//

	int group_size = grid[0] * grid[1];
	int * sendcounts = malloc(group_size * sizeof(int));
	int * displs = malloc(group_size * sizeof(int));
 	double * ptr = malloc(1);

	if (rank == 0) {
		for (i = 0; i < grid[0]; i++) {
		    for (j = 0; j < grid[1]; j++) {
		       displs[grid[1]*i + j] = padded_all_com[1]*local[0]*i+ local[1]*j;
		       sendcounts[grid[1]*i + j] = 1;
		    }
		}
	}

	if (rank == 0) {
		free(ptr);
		ptr = &U[0][0];
	}

	MPI_Scatterv(ptr, sendcounts, displs, mpi_block_all, &present_val[1][1], 1, local_block, 0, MPI_COMM_WORLD);

//---------------------/initialize 2d array for local computatiomn  and free 2D-----------------------//
	init2d(last_val, local[0]+2, local[1]+2);
	init2d(present_val, local[0]+2, local[1]+2);

    if (rank==0) 
        free2d(U,padded_all_com[0],padded_all_com[1]);
//------------Define MPI data type for local computation-----------------------------//
    MPI_Datatype column;
    MPI_Type_vector(local[0] + 2, 1, local[1] + 2 , MPI_DOUBLE, &place_holder);
    MPI_Type_create_resized(place_holder, 0, sizeof(double), &column);
    MPI_Type_commit(&column);




//---- Ensure neighbor communication and take care of boundary conditions  ----------------//

    int north, south, east, west;

    MPI_Cart_shift(CART_COMM, 0, 1, &north, &south); // If -1 then neighbor doesn't exist
    MPI_Cart_shift(CART_COMM, 1, 1, &west, &east);	


    int i_min,i_max,j_min,j_max,req_length;

	i_min = 1;
	i_max = local[0]+1;
	j_min = 1;
	j_max = local[1]+1;
	req_length = 8;

	if (north < 0) {
		req_length -= 2;
		i_min++;
	}
	if (east < 0) {
		req_length -= 2;
		if (padded[1] == 1)
			j_max -= 2;
		else
			j_max --;
	}
	if (south < 0) {
		req_length -= 2;
		if (padded[0] == 1)
			i_max -= 2;
		else
			i_max--;
	}
	if (west < 0) {
		req_length -= 2;
		j_min++;
	}

//  ------------ allocate memory for convergence check
    int * converged_arr = malloc(sizeof(int));
    int * array_convergence_store = malloc(sizeof(int));

// -- MPI reuest  initialised for use with MPI_Send ------------------------//
	MPI_Request * requests = malloc(req_length * sizeof(MPI_Request));
	MPI_Status * statuses = malloc(req_length * sizeof(MPI_Status));

	gettimeofday(&tts,NULL);


	#ifdef TEST_CONV
    for (t=0;t<T && !convergence_check;t++) {
	#endif
	#ifndef TEST_CONV
	#undef T
	#define T 256
	for (t=0;t<T;t++) {
	#endif

	//  --- communicaion  and Gaus seidel computation -----------------//

		swap=last_val;
		last_val=present_val;
		present_val=swap;   

		req_length = 0;

		if (north >= 0) {
			MPI_Irecv(&present_val[0][0], local[1] + 2, MPI_DOUBLE, north, t, MPI_COMM_WORLD, &requests[req_length]);
			req_length++;
		}
		if (west >= 0) {
			MPI_Irecv(&present_val[0][0], 1, column, west, t, MPI_COMM_WORLD, &requests[req_length]);
			req_length++;
		}
		if (req_length) {
			MPI_Waitall(req_length, requests, statuses);
			req_length = 0;
		}
		gettimeofday(&tcs,NULL);

		Gauss_Computation(last_val,present_val,i_min,i_max,j_min,j_max, relaxation_parameter);

		gettimeofday(&tcf,NULL);
		tcomp+=(tcf.tv_sec-tcs.tv_sec)+(tcf.tv_usec-tcs.tv_usec)*0.000001;

		if (north >= 0) {
			MPI_Isend(&present_val[1][0], local[1] + 2, MPI_DOUBLE, north, t, MPI_COMM_WORLD, &requests[req_length]);
			req_length++;
		}
		if (south >= 0) {
			MPI_Isend(&present_val[local[0]][0], local[1] + 2, MPI_DOUBLE, south, t, MPI_COMM_WORLD, &requests[req_length]);
			req_length++;
			MPI_Irecv(&present_val[local[0]+1][0], local[1] + 2, MPI_DOUBLE, south, t, MPI_COMM_WORLD, &requests[req_length]);
			req_length++;
		}
		if (east >= 0) {
			MPI_Isend(&present_val[0][local[1]], 1, column, east, t, MPI_COMM_WORLD, &requests[req_length]);
			req_length++;
			MPI_Irecv(&present_val[0][local[1]+1], 1, column, east, t, MPI_COMM_WORLD, &requests[req_length]);
			req_length++;
		}
		if (west >= 0) {
			MPI_Isend(&present_val[0][1], 1, column, west, t, MPI_COMM_WORLD, &requests[req_length]);
			req_length++;
		}
		MPI_Waitall(req_length, requests, statuses);

		#ifdef TEST_CONV
        if (t%C==0) { 

			// ----- Test fo rconvergence and broadcast result ------------------------// 

			converged = converge(last_val, present_val, i_min, i_max, j_min, j_max);
			converged_arr[0] = converged;
			MPI_Reduce(&converged_arr[0], &array_convergence_store[0], 1, MPI_INT, MPI_PROD, 0, MPI_COMM_WORLD);
			MPI_Bcast(&array_convergence_store[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
			convergence_check = array_convergence_store[0];
		}
		#endif

    }
    gettimeofday(&ttf,NULL);

    // --- Total timetaken --------------------------//

    ttotal=(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001;

    MPI_Reduce(&ttotal,&total_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&tcomp,&execution_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

 // -----------------Rannk zero gathers results after computation ---------------------//
    if (rank==0) {
            U=allocate2d(padded_all_com[0],padded_all_com[1]);
    }

	if (rank == 0) {
		ptr = &U[0][0];
	}

	MPI_Gatherv(&present_val[1][1], 1, local_block, ptr, sendcounts, displs, mpi_block_all, 0, MPI_COMM_WORLD);


//---------print result out from rank zero only -----------------------------------------------------//
    if (rank==0) {
		printf(" Problem size X  %d by  Y %d and Grid Size  %d by  %d Iterations %d ComputationTime %lf ",global[0],global[1],grid[0],grid[1],t,execution_time);

		#ifdef PRINT_RESULTS
		 // print  results here
		#endif

    }

    MPI_Finalize();  

    return 0;
}
