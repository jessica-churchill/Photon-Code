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
//Compile using: g++ -Wall -o thermal_photon_rateVStime.out thermal_photon_rateVStime.cpp -lgsl -lgslcblas
//**************************************************************************************************************

//**************************************************************************************************************
#define pi 3.14159265358979323846

const double pp_const = 0.5; //GeV
const double pz_const = 0.5; //GeV
const double hbarc = 0.1973269788; //GeV fm

double T, mu;

double I(double pp, double pz) //integrand, eq 69 in notes
{
    double alpha, alphas, L, g;

    alphas = 0.4;
    alpha = 1.0/137.0;
    g = sqrt(5.0);
    L = 0.75;

    double i = (40.0/(9.0*pi*pi))*alpha*alphas*exp(-sqrt((pp_const*pp_const)+(pz_const*pz_const))/T)*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*2.0*exp(-sqrt((pp*pp)+(pz*pz))/T)*exp((2.0*mu)/T)*pow(hbarc,-4.0);

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

#define COLS 4 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 1.0;

    // cout << "Enter number between 1-2 as a value for Q_s (GeV): " << endl;
    // cin >> Q_s;
    // cout << "Value of " << Q_s << "GeV chosen for Q_s." << endl;

    double f_0;
    cout << "f_0? " << endl;
    cin >> f_0;
    cout << "f_0 = " << f_0 << endl;

    ofstream myfile;
    myfile.open("thermal_photon_rateVStime_f0=" + std::to_string(f_0) + ".dat");

    fstream file;
    vector < vector <double> > energy_density; // 2d array as a vector of vectors
    vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
    int row = 0; // Row counter
    file.open("temperature_chempotential_f0="+ std::to_string(f_0) + ".dat");
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

    for (int i=0; i<489; i++)
    {
        double tau = energy_density[i*10][0];

        T = energy_density[i*10][1]*Q_s;
        mu = energy_density[i*10][2]*Q_s;

        double pp1 = 0.0130;
        double pp2 = 3.0130;
        double R = quad2d(I, pp1, pp2)*(6.0/5.0);

        myfile << std::scientific << tau/Q_s << "\t" << R*pow(Q_s, 4.0) << "\n";
        cout << std::scientific << tau/Q_s << "\t" << R*pow(Q_s, 4.0) << "\n";
    }

    myfile.close();

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
