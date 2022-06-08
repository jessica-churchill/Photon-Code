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
//Compile using: g++ -Wall -o nonthermal_photon_rate.out nonthermal_photon_rate.cpp -lgsl -lgslcblas
//**************************************************************************************************************


double fq[10201];
double fg[10201];
int N = 101;
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
#define pi 3.14159265358979323846

double pp_const;
const double pz_const = 0.5; //GeV

double I(double pp, double pz) //integrand
{
    double alpha, alphas, L, g, T, hc;

    alphas = 0.4;
    alpha = 1.0/137.0;
    g = sqrt(5.0);
    L = 0.75;
    T = 0.16; //Gev
    hc = 0.1973269788; //GeV fm

    double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*f_q(pz_const,pp_const)*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*(f_g(pz,pp)+f_q(pz,pp))*pow(hc,-4.0);

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

    double abs_error = 1.0e-7;  
    double rel_error = 1.0e-7;  
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

    double abs_error = 1.0e-7;  
    double rel_error = 1.0e-7;  
    double result;  
    double error;

    ppsav=pp;

    gsl_integration_qag (&integrand, -3.01, 3.01, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);

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

    fstream file;
    vector < vector <double> > array; // 2d array as a vector of vectors
    vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
    int row = 0; // Row counter

    double f_0;
    cout << "f_0? " << endl;
    cin >> f_0;
    cout << "f_0 = " << f_0 << endl;

    double time_step;
    cout << "Time step? " << endl;
    cin >> time_step;
    cout << "Time step = " << time_step << endl;

    //file.open("/Users/JessicaChurchill/Desktop/f0=" + std::to_string(f_0) + "/" + std::to_string(time_step) + "_tracking.dat", ios::in);
    //cout << "/Users/JessicaChurchill/Desktop/f0=" + std::to_string(f_0) + "/" + std::to_string(time_step) + "_tracking.dat" << endl;
    file.open("./f0=" + std::to_string(f_0) + "/" + std::to_string(time_step) + "_tracking.dat", ios::in);
    cout << "./f0=" + std::to_string(f_0) + "/" + std::to_string(time_step) + "_tracking.dat" << endl;
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

    for (int k=0; k<10201; k++)
    {
        fq[k]=array[k][3];
        fg[k]=array[k][2];
    }

    double Q_s = 1.0;
    
    // cout << "Enter number between 1-2 as a value for Q_s (GeV): " << endl;
    // cin >> Q_s;
    // cout << "Value of " << Q_s << "GeV chosen for Q_s." << endl;

    // cout << "p_perp:         " << "E dN/d^4xd^3p:     " << endl;

    // ofstream myfile;
    // myfile.open("nonthermal_photon_rate_" + std::to_string(Q_s) + "Qs.dat");

    ofstream myfile;
    myfile.open("nonthermal_photon_rate_f0=" + std::to_string(f_0) + "_t=" + std::to_string(time_step) +".dat");


    for (int i=0; i<=29; i++)
    {
        pp_const = 0.1000*i;

        double pp1 = 0.01;
        double pp2 = 3.01;
        double R = quad2d(I, pp1, pp2);

        
        cout << std::scientific << pp_const << "    " << R*(6.0/5.0)*pow(Q_s, 4.0) << endl;
        myfile << std::scientific << pp_const << "    " << R*(6.0/5.0)*pow(Q_s, 4.0) << endl;
    }
    myfile.close();

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
