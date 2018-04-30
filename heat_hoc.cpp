#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "math.h"
using namespace std;

int main(){
 
/******************************************/
    double k;
    k=0.25;
    int n,m;
    n=20;//no. of grid points.
    double del_x = 1/(double(n)-1);
    m=1000;
    double del_t;
    del_t = 0.001;
    double R = k*(del_t/(del_x*del_x*2.0));   // Courant number
    cout << "\nR = " << R << endl << endl;  
    int m_counter;//time counter
    //double theta, alpha , beta;
    //theta =0.49;
     double onebyh_square= 1.0/(del_x*del_x);
    
    vector<double> x(n),u0(n),udiff(n), b(n);
    vector<double> u1(n),y(n),x1(n);
    
    double a,c,b1;
    c= 12.0/11.0;
    a= 3.0/44.0;

    double A[n][n];
    double B[n][n];
    double L[n][n];
    double U[n][n];
    double Cof1[n][n];
    double Cof2[n][n];
/***************************************/
    for(int j=0; j<n; j++){

        x1[j]=0.0+(j*del_x);

       }

    for (int i = 0; i < n; i++)
    {
        u0[i]=1.0;
        b[i]=0.0;
        u1[i]=0.0;
        udiff[i]=0.0;
    }   
/***************************************///very important to make sure all the elements of the declared 2d arrays to zero , otherwise it is taking garbage values instead of zero.
        for (int i=0; i<n; i++){

            for (int j = 0; j <n; j++)
            {
                U[i][j]=0.0;
                L[i][j]=0.0;
                A[i][j]=0.0;
                B[i][j]=0.0;
                Cof1[i][j]=0.0;
                Cof2[i][j]=0.0;
                
            } 

        }

/*****************************************
  The scheme used has been taken from sec 2.2.7 Lele.
  1.Compact Fdm schemes with spectral like resolution (1992)
  2.Compact Fdm schemes on non-uniform meshes.Application to direct numerical simulations of compressible flows.
/*********************************************/
    for(int j=2;j<n-2;j++)
    {
        A[j][j-1]=2.0/11.0;
        A[j][j]  =1.0;
        A[j][j+1]=2.0/11.0;
    }
    A[0][0]=1.0;
    A[0][1]= 11.0;
    A[1][0]=1.0/10.0;
    A[1][1]=1.0;
    A[1][2]=1.0/10.0;
    A[n-2][n-3]=1.0/10.0;
    A[n-2][n-2]=1.0;
    A[n-2][n-1]=1.0/10.0;
    A[n-1][n-2]=1.0;
    A[n-1][n-1]=11.0;


    /************************************/
    for (int i = 0; i < n; i++) {
 
        // Upper Triangular
        for (int k = i; k < n; k++) {
 
            // Summation of L(i, j) * U(j, k)
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum = sum + (L[i][j] * U[j][k]);
 
            // Evaluating U(i, k)
            U[i][k] = A[i][k] - sum;
        }
 
        // Lower Triangular
        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1; // Diagonal as 1
            else {
 
                // Summation of L(k, j) * U(j, i)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum = sum + (L[k][j] * U[j][i]);
 
                // Evaluating L(k, i)
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }

    
/****************************************************/

    //start of the m_counter loop 
    for(m_counter = 0; m_counter < m ; m_counter++)
    {
        /*************************/
             
            //u1[0]=1.0;
            //u1[n-1]=2.0;
        u0[0]=1.0;   u0[19]=2.0;

    /*****************************************************/
        for(int j=2; j<(n-2); j++){                  
 
        b[j] = onebyh_square*((c*(u0[j+1]-2.0*u0[j]+u0[j-1])) + (a*(u0[j+2]-2.0*u0[j]+u0[j-2])));
                
        }

        b[0]=   onebyh_square*((13.0*u0[0])-(27.0*u0[1])+(15.0*u0[2])-u0[3]);
        b[1]=   onebyh_square*(6.0/5.0)*(u0[0]-2.0*u0[1]+u0[2]);

      //  b[n-1]= onebyh_square*((13.0*u0[6])-(27.0*u0[7])+(15.0*u0[8])-u0[9]);
      //  b[n-2]= onebyh_square*(6.0/5.0)*(u0[7]-2.0*u0[8]+u0[9]);

        b[n-1]= onebyh_square*((13.0*u0[n-4])-(27.0*u0[n-3])+(15.0*u0[n-2])-u0[n-1]);
        b[n-2]= onebyh_square*(6.0/5.0)*(u0[n-3]-2.0*u0[n-2]+u0[n-1]);


//Inverting the Lower triangular Matrix

    for(int i=0;i<n;i++)
        {                   
        double sum=0;

        for(int j=0;j<i;j++)
            {sum=sum+L[i][j]*y[j];}

        y[i]= (b[i]- sum)/L[i][i];

        }

//  Inverting the upper triangular matix ********/


   // u1[n-1]=2.0;
    //u1[0]=1.0;

    for(int i=n-1;i>=0;i--)
    {
       double sum=0;
        for(int j=n-1;j>i;j--)
            {sum = sum + U[i][j]*udiff[j];
            }

        udiff[i]=(y[i]-sum)/U[i][i];
    }



/*********************************/
for (int i = 0; i < n; i++)
{
    u1[i] = u0[i]+ k*(del_t/2.0)*udiff[i];
}
/****************************************/   
//repeating the entire process again 

      u0 = u1;




//printing the necessary files
        if ( m_counter%10 == 0)
        {

            char filename5[150];
            sprintf(filename5,"hoc_results/temp_%d.txt",m_counter);
            std::ofstream fOuttemp;
            fOuttemp.open(filename5);
            for(int i=0;i<n;i++)
              { 
              //for(int j=0;j<n;j++)
                //{
                fOuttemp<<x1[i];
                fOuttemp<<"\t\t";
                fOuttemp<<u1[i];
                //}
                fOuttemp<<std::endl;
              }
        
            fOuttemp.close();
        }


    }//end m_counter

    cout << "Results from last time step u(x," << m*del_t << "): \n\n";
    cout << setprecision(3) << setw(5) << "x[j]" << "\t" << setprecision(3) << setw(5) << "u[j] \n\n";
    
    for(int j=0; j<n; j++){
        cout << setprecision(3) << setw(5) << x1[j] << "\t" << setprecision(3) << setw(5) << u1[j]<< endl;
        }


/******************************************/
return 0;
}