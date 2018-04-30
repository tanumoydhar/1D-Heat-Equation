#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
//#include <cmath>
#include "math.h"
using namespace std;

int main(){
 
/******************************************/
    double k;
    k=0.25;
    int n,m;
    n=10;
    double del_x = 1/(double(n)-1);
    m=10000;
    double del_t;
    del_t = 0.1;
    double R = k*(del_t/(del_x*del_x*2.0));    // Courant number
    cout << "\nR = " << R << endl << endl;  
    int m_counter;
    double theta;
    theta =0.5;//crank nicolson implementation

    vector<double> x(n),u0(n), b(n);
    vector<double> u1(n),y(n),x1(n);

    double d,a,c;

    double A[n][n];
    double L[n][n];
    double U[n][n];
/*******************************************/
    for(int j=0; j<n; j++){

        x1[j]=0.0+(j*del_x);

       }
    //Solving matrix Ax=b;
    d =(1.0+2.0*R*theta);
    a = -R*theta;
    c = -R*theta;
    cout<<"d="<<'\n'<<d<<endl;
    cout<<"a="<<'\n'<<a<<endl;
/***************************************/
        for (int i=0; i<n; i++){

            for (int j = 0; j <n; j++)
            {
                U[i][j]=0.0;
                L[i][j]=0.0;
                A[i][j]=0.0;
                
            }

        }
/********************************************/

    for(int i = 1; i< n-1; i++)
    {   

    		A[i][i] = d;
    		A[i][i-1] = a;
    		A[i][i+1] = c;
    }
    A[0][0] = 1.0;
    A[n-1][n-1] =1.0;
    /************************/

     for (int i = 0; i < n; i++)
    {
        u0[i]=1.0;
    }
    //LU decomposition starts 
    /******************     
    for(int i=0;i<n;i++)
    {
        L[i][i]=1;

        for(int j=0;j<i;j++)
        {   double sum1=0;
            for(int k=0;k<j;k++)
            {
            sum1=sum1+L[i][k]*U[k][j];
            }
            L[i][j]= (A[i][j]-sum1)/U[j][j];
        }

        for(int j=i;j<n;j++)
        {
            double sum2=0;
            for(int k=0;k<i;k++)
            {
            sum2=sum2+L[i][k]*U[k][j];
            }
            
            U[i][j]=A[i][j]-sum2;
        }
    
        
        
    }
    /**************************************/
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
    cout<<"A ="<<endl;
    for (int i = 0; i < n; i++)
    {
            for (int j = 0; j < n; j++)
            {
            cout<<A[i][j];
            cout<<'\t';

            }
            cout<<'\n';

    }
    cout<<"L ="<<endl;

    for (int i = 0; i < n; i++)
    {
            for (int j = 0; j < n; j++)
            {
            cout<<L[i][j];
            cout<<'\t';

            }
            cout<<'\n';
            
    }
    cout<<"U ="<<endl;

    for (int i = 0; i < n; i++)
    {
            for (int j = 0; j < n; j++)
            {
            cout<<U[i][j];
            cout<<'\t';

            }
            cout<<'\n';
            
    }
      

    /*****************************************************/
    //start of the m_counter loop 
    for(m_counter = 0; m_counter < m ; m_counter++)
    {
        
        for(int j=1; j<(n-1); j++){                  
 
        b[j] = (R*(1.0-theta))*u0[j+1]+ (1.0-2.0*R*(1.0-theta))*u0[j]+ (R*(1.0-theta))*u0[j-1];
                
        }
        b[0]=1.0;   b[n-1]=2.0;
        /*******************************************/  
    
//Inverting the Lower triangular Matrix

    for(int i=0;i<n;i++)
        {                   
        double sum=0;

        for(int j=0;j<i;j++)
            {sum=sum+L[i][j]*y[j];}

        y[i]= (b[i]- sum)/L[i][i];

        }
//  Inverting the upper triangular matix ********/


    u1[n-1]=2.0;
    u1[0]=1.0;

    for(int i=n-2;i>=1;i--)
    {
       double sum=0;
        for(int j=n-1;j>i;j--)
            {sum = sum + U[i][j]*u1[j];
            }
        u1[i]=(y[i]-sum)/U[i][i];
    }
     

//repeating the entire process again 

      u0 = u1;

//printing the necessary files
        if ( m_counter%10 == 0)
        {

            char filename5[150];
            sprintf(filename5,"implicit_results/temp_%d.txt",m_counter);
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

/***********************************************************************/    
        /****************************
        L[i][i]=1;
                
        for(int j=0;j<i;j++)
        {   double sum_1=0.0;
            for(int k=0;k<j;k++)
            {
            sum_1=sum_1+L[i][k]*U[k][j];
            }
            L[i][j]= (double)(A[i][j]-sum_1)/(double)U[j][j];
           // cout<<L[i][j]<<endl;
        }
       /******************************/ 

/**************************************************


        for (int i = 1; i < n-1; i++)
        {   
            
            A[i][i]= A[i][i]-((double)A[i][i-1]/(double)A[i-1][i-1])*A[i-1][i];

            b[i]= b[i]-((double)A[i][i-1]/(double)A[i-1][i-1])*b[i-1];
            
        }
/******************************************************/
        /**************************************************************
        
        u1[n-1] = 2.0; //b[n-1]/d[n-1];

        for(int i=n-2; i >= 1; i--){

             u1[i] = (double)(b[i]- A[i-1][i]*u1[i+1])/(double)A[i][i];
        }
        u1[0] = 1.0;
/********************************************************************/  