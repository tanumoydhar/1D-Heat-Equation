#include <iostream>
#include <fstream>
#include <vector>
//#include"vector.h"
#include <iomanip>
#include <cmath>
using namespace std;

int main(){
 
    // Setting PDE and discretization parameters
    // ---------------------------------------------------------------
    double k;
    k=0.25;
    //cout << "\nThermal conductivity (k) = "; cin >> k;
    int n,m;
    n=10;
    //cout << "\nNumber of grid points: "; cin >> n;
    double del_x = 1/(double(n)-1);
    //cout << "\nSpatial step size (del_x) = " << del_x << endl;
    m=10000;
    //cout << "\nNumber of time steps: "; cin >> m;
    double del_t;
    //cout << "\nTime step size (del_t) = "; cin >> del_t; cout;
    del_t = 0.01;
    double R = k*(del_t/pow(del_x,2));      // Courant number
    cout << "\nR = " << R << endl << endl;  // For FTCS this must be =< 1/2 for stable solution
    double m_counter;
    

 
    // Constructing initial condition vector and x-coordinate vector
    // ---------------------------------------------------------------
    vector<double> x0(n),u0(n);
    vector<double> u1(n), u2(n), u0temp(n), u1temp(n) ;
    for(int j=0; j<n; j++){
        x0[j]=0.0 + (j*del_x);
        
    }
    /**************************************************************/
    u1=u0;
   // cout<<"flag"<<endl;

    //for (int j=0; j<n; j++)
    //fout << x[j] << "," << 0 << "," << u[j] << "\n";    // print initial condition results in file 
    for (int m_counter = 0; m_counter < m; m_counter++)
    {
       
    //for(int i=0; i<m; i++){                              // time stepping 
        
        u2[0]=0.0;                                         // setting boundary condition at u(0)
        cout<<"flag"<<endl;
        for(int j=1; j<(n-1); j++){                  
 
            u2[j]= ((double)(1-(2.0*R))/(double)(1+(2.0*R)))*u0[j] + ((double)(2.0*R)/(double)((1+(2.0*R))))*(u1[j+1]+u1[j-1]);      // Forward time centred space difference scheme
    //        cout<<x0[j]<<'\t'<<u0[j]<<endl;
        }

        u2[n-1]=2.0;                                       // setting boundary condition at u(1)

        u0temp=u0;
        u1temp=u1;
        u1=u2;
        u0=u1temp;
 
        //for (int j=0; j<n; j++)
        //fout << x0[j] << "," << (i+1)*del_t << "," << u0[j] << "\n";
    //}
        //Time  filtering used by Aghor
    
    for (int j = 0; j < n; ++j)
    {
      u0[j]= 0.5*(u0temp[j]+u1temp[j]);
      u1[j]= 0.5*(u1temp[j]+u1[j]);
    //  cout << setprecision(3) << setw(5) << x0[j] << "\t" << setprecision(3) << setw(5) << u0[j]<< endl;
    }
    /***************************************************************/
    // Print last time step results on terminal
    // --------------------------------------------------------------
        if ( m_counter%10 == 0)
        {

            char filename5[150];
            sprintf(filename5,"dufortresults/temp_%d.txt",m_counter);
            std::ofstream fOuttemp;
            fOuttemp.open(filename5);
            for(int i=0;i<n;i++)
              { 
              //for(int j=0;j<n;j++)
                //{
                fOuttemp<<x0[i];
                fOuttemp<<"\t\t";
                fOuttemp<<u0[i];
                //}
                fOuttemp<<std::endl;
              }
        
            fOuttemp.close();
        }


    }
    cout << "Results from last time step u(x," << m*del_t << "): \n\n";
    //cout << setprecision(3) << setw(5) << "x[j]" << "\t" << setprecision(3) << setw(5) << "u[j] \n\n";
    for(int j=0; j<x0.size(); j++){
       cout << setprecision(3) << setw(5) << x0[j] << "\t" << setprecision(3) << setw(5) << u0[j]<< endl;
    }
        return 0;

}