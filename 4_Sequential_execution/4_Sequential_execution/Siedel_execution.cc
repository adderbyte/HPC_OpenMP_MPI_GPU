//Pivotisation(partial) to make the equations diagonally dominant//Simple Gauss Siedel Method sequential executiom
//
#include<iostream>
#include<iomanip>
#include<math.h>
#include<random>
#include <sstream>
#include <chrono>
// ---Tolerance---
typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> second;
double  *  Gauss_Siedel (float **input, int rowsize, double *ouput);

using namespace std;
int main( int argc , char *argv[] )
{
    //check argument input

     if (argc != 2)
    {
        std::cerr << "You have to specify the  array dimensions" << std::endl;
        return -1;
    }


    int sizeX, sizeY;
    int output_column = 1;// extend array by 1 column. This column indicate the ouput
    sizeX = std::stoi(argv[1])+ output_column;
    sizeY = std::stoi(argv[1]) ;

    // check input size array is not equals to zero

    if (sizeX <= 0)
    {
        std::cerr << "Invalid dimension x" << std::endl;
        return -1;
    }
    if (sizeY <= 0)
    {
        std::cerr << "Invalid dimension y" << std::endl;
        return -1;
    }
    int i,j,k; // initialise variables for use as counters.
    //Start timer 
    auto start = clk::now();
    //set output precision
    cout.precision(4);
    cout.setf(ios::fixed);
    // Initialise random timer for use in random number generation
    srand (time(0));
    // declare the  augmented array input here
    float ** augmented_Array_Input;

    // Linear memory allocation to make contiguous mem segment
     float * tempStore = new  float[sizeX * sizeY];

    // Allocate the pointers inside the array,
    // which will be used to index the linear memory
    augmented_Array_Input = new float*[sizeY];

    // Let the pointers inside the array point to the correct memory addresses
    for (int i = 0; i < sizeY; ++i)
    {
        augmented_Array_Input[i] = (tempStore + i * sizeX);
    }
    // Fill the augmented array
    for (int y = 0; y < sizeY; ++y)
    {
        for (int x = 0; x < sizeX; ++x)
        {
            augmented_Array_Input[y][x] =  ((double)rand() / (double)RAND_MAX)  + (y+10)*20 + x;
        }
    }

    // Make the array diagonally dominant
    for (int i=0;i<sizeY;i++)
        for ( int j=0;j<=sizeY-1;j++)
             if(i ==j)
                augmented_Array_Input[i][j] =  augmented_Array_Input[i][j] +( j+20000) *2  + (i+1)*10 ;

    //Pivotisation(partial) to make the equations diagonally dominant
    for (i=0;i<sizeY;i++)                    
        for (k=i+1;k<sizeY;k++)
            if (fabs(augmented_Array_Input[i][i])<fabs(augmented_Array_Input[k][i]))
                for (j=0;j<=sizeY;j++)
                {
                    double temp=augmented_Array_Input[i][j];
                    augmented_Array_Input[i][j]=augmented_Array_Input[k][j];
                    augmented_Array_Input[k][j]=temp;
                }
    ///instantiate and initialise array  to store solution for each iteration
     double *  solution_to_Equation = new double[sizeY]();

    /*Uncomment to store result and print later
     declare result array  to store results from  Gauss_Siedel Function
     double * result ;
     call the Gauss Sieded method
     result = Gauss_Siedel( augmented_Array_Input ,sizeY,solution_to_Equation);
     cout<<"\n The solution is as follows:\n";
     for (i=0;i<sizeY;i++)
        cout<<"solution_toEquation"<<i<<" = "<< *(result + i)<<endl;        //Print the contents of x[]
      */

     // Call Gauss siedel method
     Gauss_Siedel( augmented_Array_Input ,sizeY,solution_to_Equation);
    //End timer
     auto end = clk::now();
     // compute the time taken
     second time = end - start;
     cout<<endl;
     cout << "Total Implementation time:"<<endl;
     cout << time.count() << endl;
     
     /* Uncomment this block to view solution
    cout<<"\n The solution is as follows:\n";
    for (i=0;i<sizeY;i++)
        cout<<"solution_toEquation"<<i<<" = "<< *(result + i)<<endl;        //Print the contents of x[]*/


    delete[] augmented_Array_Input[0];
    delete[] augmented_Array_Input;
    augmented_Array_Input= nullptr;
    delete[] solution_to_Equation;
    return 0;
}
double * Gauss_Siedel (float **a , int rowsize, double *x){

    int i,j,flag=0,count=0; //counter  and control variables for the loops
    // y is a temporary storage and 
    double abs_tolerance,y;
    abs_tolerance =0.0000001 ;//
    cout<<"\n--------------- Inside Gauss Siedel Function-------------------------------------------------------";
    do                            //Perform iterations to calculate x1,x2,...xn
    {
       for (i=0;i<rowsize;i++)     //Loop that calculates x1,x2,...xn
        {
            y=x[i];
            x[i]=a[i][rowsize];
            for (j=0;j<rowsize;j++)
            {
                if (j!=i)
                x[i]=x[i]-a[i][j]*x[j];
            }
            x[i]=x[i]/a[i][i];
            if (abs(x[i]-y)<abs_tolerance)            //Compare the the value with the last value
                flag++;
           count++;
        }

       
    }while(flag<rowsize);
    cout<<endl;
    cout <<"Total iterations for Gauss methods convergence" << endl;
    cout<<count<<endl;
    return  x ;

      }
