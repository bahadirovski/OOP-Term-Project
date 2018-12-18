#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <cmath>
//#include "netcdf.h"

#define g_const 9.8065
#define epsilon 0.1

using namespace std;

template<class T = double,class M = int>
T **allocater(M row, M col, T inival = 0)
{
    T** A = new T*[col];
    T* B = new T[row*col];
    for (M i = 0; i < col; i++)
    {
        A[i] = &(B[i*row]);
    }

    for (M i = 0; i < col; i++)
    {
        for (M j = 0; j < row; j++)
        {
            A[i][j] = (T)inival;
        }
    }
    return A;
}


template<class T = double, class M = int>
class simulate
{
private:
    T x, y, dx, dy, delta_t, sim_t, h_west, h_east, h_north, h_south;
    T **h,**h_zero, **u, **v, **zeta, **zeta_star;
    M nx, ny;
public:
    simulate()
    {
        double inival = 10.;
        x = 500.;
        y = 500.;           
        dx = 10.;
        dy = 10.;
        delta_t = 0.01;
        sim_t = 100.;
        h_west = 0.; h_east = 0.; h_south = 0.; h_north = 0.;
        nx = M(x/dx) + 2;
        ny = M(y/dy) + 2;
        h = allocater<>(nx,ny,inival);
        h_zero = allocater<>(nx,ny,inival);
        u = allocater<>(nx,ny);
        v = allocater<>(nx,ny);
        zeta = allocater<>(nx,ny);
        zeta_star = allocater<>(nx,ny);        
    }
    
    simulate(T x, T y, T dx, T dy, T delta_t = 0.01, T sim_t = 100.)
    {
        double inival = 10.;
        (*this).x = x;
        (*this).y = y;
        (*this).dx = dx;
        (*this).dy = dy;
        (*this).delta_t = delta_t;
        (*this).sim_t = sim_t;
        h_west = 0.; h_east = 0.; h_south = 0.; h_north = 0.;
        nx = M(x/dx) + 2;
        ny = M(y/dy) + 2;
        h = allocater<>(nx,ny,inival);
        h_zero = allocater<>(nx,ny,inival);
        u = allocater<>(nx,ny);
        v = allocater<>(nx,ny);
        zeta = allocater<>(nx,ny);
        zeta_star = allocater<>(nx,ny);        
    }   
    
    void setDepth(T val, int sim_no = 0)
    {
        
        //(dx<dy)?dx:dy
        if (delta_t > (min(dx,dy))/pow(2*g_const*((T)(*(max_element(h[0], h[0] + nx*ny)))),1/2))
        {
            cout << "Stability Criteria is not satisfied" << endl;
        }
        else
        {
        	if(sim_no == 0)
        	{
        		for (int i=0; i<nx; i++)
        		{
        			for (int j=0; j<ny; j++)
        			{
        				h[i][j] = val;
        				h_zero[i][j] = val;
        			}
        		}
        	}
        	
        	else
        	{
        		for (int i=0; i<nx; i++)
        		{
        			double slope = (val - val/2)/(nx-1);
        			for (int j=0; j<ny; j++)
        			{
        				h[i][j] = val - slope*j;
        				h_zero[i][j] = val - slope*j;
        			}
        		}
        	}
        }
    }
    
    void setx(T x) 
    { 
        (*this).x = x; 
        nx = M(x/dx) + 2;
    }
    
    void sety(T y) 
    { 
        (*this).y = y; 
        ny = M(y/dy) + 2;
    }
    
    void setdx(T dx) 
    { 
        (*this).dx = dx; 
        nx = M(x/dx) + 2;
    }   
    
    void setdy(T dy) 
    { 
        (*this).dy = dy; 
        ny = M(y/dy) + 2;
    }
    
    T getx() { return x; }
    T gety() { return y; }
    T getdx() { return dx; }
    T getdy() { return dy; }
    
    void show()
    {
        cout << "Grid Resolution is = " << nx << "x" << ny << " = " << nx*ny << endl;
        cout << "Maximum Depth is = " << *(max_element(h[0], h[0] + nx*ny)) << endl;
        cout << "Minimum Depth is = " << *(min_element(h[0], h[0] + nx*ny)) << endl;
        cout << "Simulation Time is = " << sim_t << endl;
        cout << "Time Interval is = " << delta_t << endl;
    }

    void solve()
    {
        zeta[nx/2][ny/2] = 1.;
        for (T t=0; t < sim_t; t=t+delta_t)
        {
            for (int i = 1; i < nx-1; i++)
            {
                for (int j = 1; j < ny-1; j++)
                {
                    u[i][j] = u[i][j] - delta_t * g_const * ( zeta[i][j+1] - zeta[i][j] ) / dx;
                    v[i][j] = v[i][j] - delta_t * g_const * ( zeta[i+1][j] - zeta[i][j] ) / dy;
                }
            }

            for (int i = 1; i < nx-1; i++)
            {
                for (int j = 1; j < ny-1; j++)
                {
                    if (u[i][j] > 0) { h_east = h[i][j]; }
                    else { h_east = h[i][j+1]; }

                    if (u[i][j-1] > 0) { h_west = h[i][j-1]; }
                    else { h_west = h[i][j]; }

                    if (v[i][j] > 0) { h_north = h[i][j]; }
                    else { h_north = h[i+1][j]; }

                    if (v[i-1][j] > 0) { h_south = h[i-1][j]; }
                    else { h_south = h[i][j]; }                    

                    zeta_star[i][j] = zeta[i][j] - delta_t * ((u[i][j] * h_east - u[i][j-1] * h_west) / dx) 
                        - delta_t * ((v[i][j] * h_north - v[i-1][j] * h_south) / dy);
                }
            }

            for (int i = 1; i < nx-1; i++)
            {
                for (int j = 1; j < ny-1; j++)
                {
                    //zeta[i][j] = zeta_star[i][j];
                    zeta[i][j] = (1. - epsilon) * zeta_star[i][j] + 0.25 * epsilon * (zeta_star[i][j-1] 
                        + zeta_star[i][j+1] + zeta_star[i-1][j] + zeta_star[i+1][j]);
                }
            }

            for (int i = 1; i < nx-1; i++)
            {
                for (int j = 1; j < ny-1; j++)
                {
                    h[i][j] = h_zero[i][j] + zeta[i][j];
                }
            }            
        
            cout << t << "  " << zeta[26][26] << endl;
        	//if ((t-0.05)<pow(10.,-5))
        }
    }
};

int main()
{
    double x = 0., y = 0., dx = 0., dy = 0., delta_t = 0., sim_t = 0.;
    ifstream infile("params.txt");
    if (infile.is_open())
    {
    	while (infile >> x >> y >> dx >> dy >> delta_t >> sim_t) {}
    }

    simulate<>s1(x,y,dx,dy,delta_t,sim_t);
    // Maximum depth degismiyor, ya da minimum
    //double a = 15;
    //double b = 15; 
    //s1.setdx(a);
    //s1.setdy(b);
    //s1.setDepth(a,2);
    s1.show();
    s1.solve();

};