#include <iostream>
#include <fstream>
#include <vector>
//#include"vector.h"
#include <iomanip>
#include <cmath>
using namespace std;
 

int main(){
 
    /******************************************/
    double k;
    k=0.25;
    //cout << "\nThermal conductivity (k) = "; cin >> k;
    int n,m;
    n=10;
    //cout << "\nNumber of grid points: "; cin >> n;
    double del_x = 1/(double(n)-1);
    //cout << "\nSpatial step size (del_x) = " << del_x << endl;
    m=1000;
    //cout << "\nNumber of time steps: "; cin >> m;
    double del_t;
    //cout << "\nTime step size (del_t) = "; cin >> del_t; cout;
    del_t = 0.001;
    double R = k*(del_t/pow(del_x,2));      // Courant number
    cout << "\nR = " << R << endl << endl;  // For FTCS this must be =< 1/2 for stable solution
    double m_counter;
 

    vector<double> x,u;
    for(int j=0; j<n; j++){
        x.push_back(j*del_x);
        u.push_back(1);
    }

    for (int m_counter = 0; m_counter < m; m_counter++)
    {
    //for (int j=0; j<n; j++)
    //fout << x[j] << "," << 0 << "," << u[j] << "\n";    // print initial condition results in file 
 
    //for(int i=0; i<m; i++){                              // time stepping 
 
        u[0]=1.0;                                         // setting boundary condition at u(0)
        for(int j=1; j<(n-1); j++){                  
 
            u[j]= u[j] + R*(u[j+1]-2*u[j]+u[j-1]);      // Forward time centred space difference scheme
        }
        u[n-1]=2.0;                                       // setting boundary condition at u(1)
 
    
        if ( m_counter%10 == 0)
        {

            char filename5[150];
            sprintf(filename5,"ftcsresults/temp_%d.txt",m_counter);
            std::ofstream fOuttemp;
            fOuttemp.open(filename5);
            for(int i=0;i<n;i++)
              { 
              //for(int j=0;j<n;j++)
                //{
                fOuttemp<<x[i];
                fOuttemp<<"\t\t";
                fOuttemp<<u[i];
                //}
                fOuttemp<<std::endl;
              }
        
            fOuttemp.close();
        }

    }
    cout << "Results from last time step u(x," << m*del_t << "): \n\n";
    cout << setprecision(3) << setw(5) << "x[j]" << "\t" << setprecision(3) << setw(5) << "u[j] \n\n";
    for(int j=0; j<x.size(); j++){
        cout << setprecision(3) << setw(5) << x[j] << "\t" << setprecision(3) << setw(5) << u[j]<< endl;
    }

    return 0;
}