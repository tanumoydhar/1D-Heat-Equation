#include <iostream>
#include <fstream>
#include <vector>
//#include"vector.h"
#include <iomanip>
#include <cmath>

/*#include "get_k1_T_rk4.cpp"
#include "get_k2_T_rk4.cpp"
#include "get_k3_T_rk4.cpp"
#include "get_k4_T_rk4.cpp"
#include "T_rk4.cpp"
/****************************************************/
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
    m=10000;
    //cout << "\nNumber of time steps: "; cin >> m;
    double del_t;
    //cout << "\nTime step size (del_t) = "; cin >> del_t; cout;
    del_t = 0.01;
    double R = k*(del_t/pow(del_x,2));      // Courant number
    cout << "\nR = " << R << endl << endl;  // For FTCS this must be =< 1/2 for stable solution
    double m_counter;
 

    vector<double> x(n),u0(n);
    vector<double> u1(n);
    vector<double> k1T(n),k2T(n),k3T(n),k4T(n),dummy_T(n);
    for(int j=0; j<n; j++){
        x[j]=0.0+(j*del_x);
        //u.push_back(1);
    }

    for (int m_counter = 0; m_counter < m; m_counter++)
    {
    

        u1[0]=1.0;                                         // setting boundary condition at u(0)
        for(int j=1; j<(n-1); j++){ 


            k1T[j]=(R*(u0[j+1]-2*u0[j]+u0[j-1]));
            //cout<<k1T[j]<<endl;
            //k2T[j]=(R*(k1T[j+1]-2*k1T[j]+k1T[j-1])); 
            k2T[j]=((u0[j+1]+0.5*k1T[j+1])-2*(u0[j]+0.5*k1T[j])+(u0[j-1]+0.5*k1T[j-1]))*R;
            k3T[j]=((u0[j+1]+0.5*k2T[j+1])-2*(u0[j]+0.5*k2T[j])+(u0[j-1]+0.5*k2T[j-1]))*R;  
            k4T[j]=((u0[j+1]+k3T[j+1])-2*(u0[j]+k3T[j])+(u0[j-1]+k3T[j-1]))*R;  
            /********************************************************/
            /*
            k2T[j]=(0.5*del_t)*(R*(u0[j+1]-2*u0[j]+u0[j-1])); 
            k3T[j]=(0.5*del_t)*k2T[j];  
            k4T[j]=del_t*k3T[j];
            /*********************************************************/
            //cout<<k4T[j]<<endl;

            dummy_T[j] = k4T[j];
            dummy_T[j] = dummy_T[j] + k3T[j];
            dummy_T[j] = dummy_T[j] + k3T[j];
            dummy_T[j] = dummy_T[j] + k2T[j];
            dummy_T[j] = dummy_T[j] + k2T[j];
            dummy_T[j] = dummy_T[j] + k1T[j];
            dummy_T[j] = ((double)1.0/(double)6.0)*del_t*dummy_T[j];

            u1[j] = (R*(u0[j+1]-2*u0[j]+u0[j-1])) + dummy_T[j];                  
 
            //u1[j]= u0[j] + ;      // Forward time centred space difference scheme
        }
        u1[n-1]=2.0;                                       // setting boundary condition at u(1)
 
    
        //T_rk4(u1, u0, R);


        u0=u1;
        //u1=u0;

        if ( m_counter%10 == 0)
        {

            char filename5[150];
            sprintf(filename5,"rk4results/temp_%d.txt",m_counter);
            std::ofstream fOuttemp;
            fOuttemp.open(filename5);
            for(int i=0;i<n;i++)
              { 
              //for(int j=0;j<n;j++)
                //{
                fOuttemp<<x[i];
                fOuttemp<<"\t\t";
                fOuttemp<<u1[i];
                //}
                fOuttemp<<std::endl;
              }
        
            fOuttemp.close();
        }

    }
    cout << "Results from last time step u(x," << m*del_t << "): \n\n";
    cout << setprecision(3) << setw(5) << "x[j]" << "\t" << setprecision(3) << setw(5) << "u[j] \n\n";
    for(int j=0; j<x.size(); j++){
        cout << setprecision(3) << setw(5) << x[j] << "\t" << setprecision(3) << setw(5) << u1[j]<< endl;
    }

    return 0;
}