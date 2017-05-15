/*Gauss_siedel seqquential implementation 
 ####    No relaxation parameter defined in this implementation. For fast convergence we would define relaxation parameter omega 
 ####    The implementation with relaxation parameter is called Side_sor. This improves convergence.The equation for this is in the PDF.
 ####    But an arbitrary convergence value was chosen. This is defined as abstol
*/
#include<iostream>
#include<iomanip>
#include<random>
#include<stdlib.h>
#include<math.h>
#include<sstream>
#include<chrono>
using namespace std;
typedef chrono::high_resolution_clock clk;
typedef chrono::duration<double> second;

int main(int argc , char *argv[] )
{
    cout.precision(4);
    cout.setf(ios::fixed);
    int n,i,j,k,flag=0,count=0;
   // cout<<"\nEnter the no. of equations\n";
   // if (argc != 2) usage(argv[0]);
    stringstream args(argv[1]);
    args>>n;                    //Input no. of equations
    double a[n][n+1];            //declare a 2d array for storing the elements of the augmented matrix
    double x[n];                //declare an array to store the values of variables
    double abstol,y;               // decalare tolerance and solution variable y

    abstol = 0.000000000001;
    srand (time(0));
   // Initialise random 2D array
    for (i=0;i<n;i++)
        for (j=0;j<=n;j++)
                a[i][j] = ((double)rand() / (double)RAND_MAX)  + (i+10)*20 + j  ;
  /*
   // Uncomment to view input variable
    for (i=0;i<n;i++)
        for (j=0;j<=n;j++)
                cout<< a[i][j] << endl;

   */
    cout << endl;
    // Increase the values at the diagonal to
    // ensure diagonally dominant matrix
    // Otherwise Nan is reported in solutions
    // And convergence will never be reached
    
    for (i=0;i<n;i++)
        for (j=0;j<=n-1;j++)
             if(i ==j)
                a[i][j] =  a[i][j] +( j+20000) *2  + (i+1)*10 ;

  /*Uncomment to see values after increasing diagonal element values
    for (i=0;i<n;i++)
        for (j=0;j<=n;j++)
                cout<< a[i][j] << endl;
  */
   
    cout<<"\n Initialising values of x variable to zero:\n";
    for (i=0;i<n;i++)
        x[i]=0;
   

    for (i=0;i<n;i++)                    //  ensure diagonal dominance as much as possible so we can have a solution
        for (k=i+1;k<n;k++)
            if (abs(a[i][i])<abs(a[k][i]))
                for (j=0;j<=n;j++)
                {
                    double temp=a[i][j];
                    a[i][j]=a[k][j];
                    a[k][j]=temp;
                }
    auto start = clk::now(); // start clock
   // cout<<"Iter"<<setw(10);// uncomment to verify output  for simple sample solutions
  /* Uncomment to view input variables 
 * for(i=0;i<n;i++)
        cout<<"x"<<i<<setw(18);
   */
    cout<<"\n--------------Performing Iterative Solutions   ----------------------------------";
    do                            //Perform iterations to calculate x1,x2,...xn
    {
       // cout<<"\n"<<count+1<<"."<<setw(16);
        for (i=0;i<n;i++)                //Loop that calculates x1,x2,...xn
        {
            y=x[i];
            x[i]=a[i][n];
            for (j=0;j<n;j++)
            {
                if (j!=i)
                x[i]=x[i]-a[i][j]*x[j];
            }
            x[i]=x[i]/a[i][i];
            if (abs(x[i]-y)<abstol)            // Arbitrary convergence criterion
                flag++;
           // cout<<x[i]<<setw(18);
        }
        count++;
    }while(flag<n); //If the values of all the variables don't differ from their previious values with error more than eps then flag must be n and hence stop the loop
    auto end = clk::now(); // end Timer
    second time = end - start;  // calculate duration
    cout << "Duration"<<endl;
    cout << time.count() << endl;
   /* cout<<"\n The solution is as follows:\n";
    for (i=0;i<n;i++)
    cout<<"x"<<i<<" = "<<x[i]<<endl;        //Print the contents of x[] // Umcomment to view solutions
    */
    cout <<"Total iterations" << endl;
    cout<<count<<endl;
    cout<< endl;
    cout<< "##################### finished ########" << endl;
    return 0;
}
