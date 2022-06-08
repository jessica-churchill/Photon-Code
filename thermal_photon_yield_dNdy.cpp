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
//Compile using: g++ -Wall -o thermal_photon_yield_dNy.out thermal_photon_yield_dNdy.cpp -lgsl -lgslcblas
//**************************************************************************************************************

//**************************************************************************************************************
#define pi 3.14159265358979323846
const double hbarc = 0.1973269788;

double pz_const;
double pp_const;
double pz_tilde_const;
double eta;
double d_eta;
double tau;
double Q_s;
double T, mu;

double I(double pp, double pz) //integrand
{
    double alpha, alphas, L;

    alphas = 0.4;
    alpha = 1.0/137.0;
    L = 0.75;

    double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*(2.0*exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-2.0);
    //double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*exp(-sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))/T)
                //*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*exp(-sqrt((pp*pp)+(pz*pz))/T)*(1.0+exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-2.0); //annihilation
    //double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*exp(-sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))/T)
                //*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*exp(-sqrt((pp*pp)+(pz*pz))/T)*(1.0-exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-2.0); //compton
                

    return i;

}
//**************************************************************************************************************



//2d integration routine
//************************************************************************************************************
static double ppsav;
static double (*nrfunc)(double, double);

double quad2d(double (*func)(double, double), double pp1, double pp2)
{
    double f1(double pp, void *params);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    gsl_function integrand;
    integrand.function = &f1;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;

    nrfunc=func; 
    gsl_integration_qag (&integrand, pp1, pp2, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    return result;
}

double f1(double pp, void *params) 
{
    double f2(double pz, void *params);
    double pz_min(double);
    double pz_max(double);

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

    gsl_function integrand;
    integrand.function = &f2;

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
    double result;  
    double error;

    ppsav=pp;

    gsl_integration_qag (&integrand, -3.0, 3.0, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);

    return result;
}

double f2(double pz, void *params)  
{
    return (*nrfunc)(ppsav,pz); 
}

//************************************************************************************************************
//main function
//************************************************************************************************************

double yield(double pp_const, void * params)
{
    pz_tilde_const = pp_const*sinh(asinh(pz_const/pp_const)-eta);
    double R = 2.0*pp_const*exp(-(sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))-mu)/T);
    return R;
}

#define COLS 12 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 2.0;
    double f_0 = 1.00;

    fstream file;
    vector < vector <double> > energy_density; // 2d array as a vector of vectors
    vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
    int row = 0; // Row counter
    file.open("f0="+ std::to_string(f_0) + ".dat");
    if (file.is_open()) { 
        while (file.good()) { 
        energy_density.push_back(rowVector); // add a new row,
        for (int col=0; col<COLS; col++) {
            file >> std::scientific >> energy_density[row][col]; // fill the row with col elements
        }
        row++; // Keep track of actual row 
        }
    }
    else cout << "Unable to open file" << endl;
    file.close();

    ofstream myfile1;
    myfile1.open("thermal_photon_dNdy_2Qs_f0=" + std::to_string(f_0) + ".dat", std::ios_base::app);

    pz_const = 0.0;

    int limit_array [] = {11, 31, 52, 72, 92};

    //double f0_array [] = {0.15, 0.30, 0.76, 2.66};

    for (int f=0; f<5; f++)
    {
        int limit = limit_array[f];
        
        double tau;
        double sum = 0.0;
        double R = 0.0;

        //for (int k=0; k<489; k++)
        for (int k=0; k<limit; k++)
        {
            double d_tau = 0.1;
            tau = 1.0+(d_tau*k);
            //tau = energy_density[k*10][0];
            T = energy_density[k*10][1]*Q_s;
            mu = energy_density[k*10][2]*Q_s;

            for(int m=0; m<=100; m++)
            {
                d_eta = 0.1;
                eta = d_eta*m;

                //pz_tilde_const = pp_const*sinh(asinh(pz_const/pp_const)-eta);

                double A = 73.90;

                double pp1 = 0.01;
                double pp2 = 3.01;

                gsl_integration_workspace *w = gsl_integration_workspace_alloc (1000);

                double result, error;

                gsl_function F;
                F.function = &yield;

                double abs_error = 1.0e-5;  
                double rel_error = 1.0e-5; 

                gsl_integration_qag (&F, 0.0, 3.0, abs_error, rel_error, 1000, 4, w, &result, &error);

                R = 2.0*A*d_tau*tau*d_eta*result*quad2d(I, pp1, pp2)*(6.0/5.0);
                sum = sum + R;
                gsl_integration_workspace_free(w);
            }
            
        }
        myfile1 << std::scientific << tau*hbarc << "\t" << sum << "\n";
        cout << std::scientific << tau*hbarc << "\t" << sum << "\n";
        
    }
    myfile1.close();
    

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
