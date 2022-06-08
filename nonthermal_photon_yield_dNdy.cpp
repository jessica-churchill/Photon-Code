#include <iostream>
#include <cmath>
#include <math.h> 
#include <iomanip>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <ctime>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <gsl/gsl_integration.h>

using namespace std;

//**************************************************************************************************************
//Compile using: g++ -Wall -o nonthermal_photon_yield_dNdy.out nonthermal_photon_yield_dNdy.cpp -lgsl -lgslcblas
//**************************************************************************************************************

#define pi 3.14159265358979323846
const double hbarc = 0.1973269788;
#define COLS 4 // Number of columns in data
#define COLS2 2
double pz_const;
double pp_const;
double pz_tilde_const;
double eta;
double d_eta;
double tau;
double d_tau;
double Q_s;

double fq[10201];
double fg[10201];
int N = 101;

//double f0_array [] = {0.15, 0.30, 0.50, 0.70, 0.76, 0.85, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.66, 3.00};
//double f0_array [] = {0.15, 0.30, 0.76, 2.66};

double fd_q(int i, int j)
{
  return fq[(i*N)+j];
}

double f_q(double pz, double pp) //distribution function of quarks: fq
{
    double f;

    double pzmin = 1.000000e-02;
    double dpz = 3.000000e-02;
    double ppmin = 1.000000e-02;
    double dpp = 3.000000e-02;

    int i, j;
    i = (int )((fabs(pz)-pzmin)/dpz);
    j = (int )((fabs(pp)-ppmin)/dpp);

    if ( (fabs(pp) > 3.0) || (fabs(pz) > 3.0) )
    {
        f = 0;
    }

    else if (i > N || j > N || i+1 > N || j+1 > N) 
    {
        f = 0;
    }

    else if (fabs(pp) < ppmin)
    {
        f = fd_q(i,0);
    }

    else
    {
        f = (fd_q(i,j)+fd_q(i+1,j+1))/2;
    }

    return f;
}

double fd_g(int i, int j)
{
  return fg[(i*N)+j];
}

double f_g(double pz, double pp) //distribution function of gluons: fg
{
    double f;

    double pzmin = 1.000000e-02;
    double dpz = 3.000000e-02;
    double ppmin = 1.000000e-02;
    double dpp = 3.000000e-02;

    int i, j;
    i = (int )((fabs(pz)-pzmin)/dpz);
    j = (int )((fabs(pp)-ppmin)/dpp);

    if ( (fabs(pp) > 3.0) || (fabs(pz) > 3.0) )
    {
        f = 0;
    }

    else if (i > N || j > N|| i+1 > N || j+1 > N) 
    {
        f = 0;
    }

    else if (fabs(pp) < ppmin)
    {
        f = fd_g(i,0);
    }

    else
    {
        f = (fd_g(i,j)+fd_g(i+1,j+1))/2;
    }

    return f;
}

//**************************************************************************************************************


//************************************************************************************************************
//main function
//************************************************************************************************************

double yield(double pp_const, void * params)
{
    pz_tilde_const = pp_const*sinh(asinh(pz_const/pp_const)-eta);
    double R = 2.0*pp_const*f_q(pz_tilde_const,pp_const);
    return R;
}

int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    Q_s = 2.0;

    double f_0 = 1.00;
    // cout << "f_0? " << endl;
    // cin >> f_0;
    // cout << "f_0 = " << f_0 << endl;
    pz_const = 0.0;

        //double tau;

        double limit_array [] = {11, 31, 51, 72, 92};

        for (int f=0; f<5; f++)
        {
            //double f_0 = f0_array[f];
            //double f_0 = 0.250000;
            //cout << "f0 = "<< f_0 << endl;
            double sum = 0.0;
            double R = 0.0;

            double limit = limit_array[f];

            ofstream myfile1;
            myfile1.open("nonthermal_photon_yield_dNdy_f0=" + std::to_string(f_0) + ".dat", std::ios_base::app);

            fstream file2;
            vector < vector <double> > integral; // 2d array as a vector of vectors
            vector <double> rowVector2(COLS2); // vector to add into 'array' (represents a row)
            int row2 = 0; // Row counter
            file2.open("photon_yield_integral_f0=" + std::to_string(f_0) + ".dat");
            if (file2.is_open()) { 
                while (file2.good()) { 
                integral.push_back(rowVector2); // add a new row,
                for (int col2=0; col2<COLS2; col2++) {
                file2 >> std::scientific >> integral[row2][col2]; // fill the row with col elements
            }
            row2++; // Keep track of actual row 
            }
            }
            else cout << "Unable to open file" << endl;
            file2.close();

            //for (int k=0; k<491; k++) 10fm
            for (int k=0; k<limit; k++)
            {
                d_tau = 0.1;
                tau = 1.0+(d_tau*k);

                fstream file;
                vector < vector <double> > array; // 2d array as a vector of vectors
                vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
                int row = 0; // Row counter

                file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/f0=" + std::to_string(f_0) + "/" + std::to_string(tau) + "_tracking.dat", ios::in);
                if (file.is_open()) { 

                    while (file.good()) { 
                    array.push_back(rowVector); // add a new row,
                    for (int col=0; col<COLS; col++) {
                        file >> std::scientific >> array[row][col]; // fill the row with col elements
                    }
                    row++; // Keep track of actual row 
                    }
                }
                else cout << "Unable to open file" << endl;
                file.close();

                for (int l=0; l<10201; l++)
                {
                    fq[l]=array[l][3];
                    fg[l]=array[l][2];
                }

                for(int m=0; m<=100; m++)
                {   
                    d_eta = 0.1;
                    eta = d_eta*m;

                    double A = 73.90;
                    double L = 0.75;
                    double alphas = 0.4;
                    double alpha = 1.0/137.0;

                    gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

                    double result, error;

                    gsl_function F;
                    F.function = &yield;

                    double abs_error = 1.0e-5;  
                    double rel_error = 1.0e-5; 

                    gsl_integration_qag (&F, 0.0, 3.0, abs_error, rel_error, 1000, 4, w, &result, &error);

                    R = 2.0*A*(d_tau/Q_s)*(tau/Q_s)*d_eta*(40.0/(9.0*pi*pi))*alpha*alphas*L*integral[k][1]*(6.0/5.0)*pow(hbarc,-2.0)*pow(Q_s,4.0)*result;
                    sum = sum + R;
                    gsl_integration_workspace_free(w);
                }
                
            }
            myfile1 << std::scientific << tau*hbarc << "\t" << sum << "\n";
            cout << std::scientific << tau*hbarc << "\t" << sum << "\n";
            myfile1.close();
        }

duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
std::cout<<"Time: "<< duration <<'\n';
}
