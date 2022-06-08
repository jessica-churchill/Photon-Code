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
//Compile using: g++ -Wall -o photon_yield_integral.out photon_yield_integral.cpp -lgsl -lgslcblas
//**************************************************************************************************************

//const double T = 0.18001331534960535; 
//const double mu = -0.0032852430051302974;
const double L = 0.75;
const double hc = 0.1973269788; //GeV fm
//const double alphas = 0.4;
const double alphas = 0.234;
const double alpha = 1.0/137.0;
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
    double dpz = 6.000000e-02;
    double ppmin = 1.000000e-02;
    double dpp = 6.000000e-02;

    int i, j;
    i = (int )((fabs(pz)-pzmin)/dpz);
    j = (int )((fabs(pp)-ppmin)/dpp);

    if ( (fabs(pp) > 6.0) || (fabs(pz) > 6.0) )
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

    //double E = sqrt((pp*pp)+(pz*pz));
    //f = 1.0/((exp(E-mu)/T)+1.0);

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
    double dpz = 6.000000e-02;
    double ppmin = 1.000000e-02;
    double dpp = 6.000000e-02;

    int i, j;
    i = (int )((fabs(pz)-pzmin)/dpz);
    j = (int )((fabs(pp)-ppmin)/dpp);

    if ( (fabs(pp) > 6.0) || (fabs(pz) > 6.0) )
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

    //double E = sqrt((pp*pp)+(pz*pz));
    //f = 1.0/((exp(E-mu)/T)-1.0);

    return f;
}

//**************************************************************************************************************
#define pi 3.14159265358979323846

double pz_const;
double pp_const;

double I_tot(double pp, double pz) //integrand
{
    double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*(f_g(pz,pp)+f_q(pz,pp));
    //double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*f_g(pz,pp)*(1.0-f_q(pz,pp)); //compton
    //double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*f_q(pz,pp)*(1.0+f_g(pz,pp)); //annihilation

    return i;

}

double I_comp(double pp, double pz) //integrand
{
    //double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*(f_g(pz,pp)+f_q(pz,pp));
    double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*f_g(pz,pp)*(1.0-f_q(pz,pp)); //compton
    //double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*f_q(pz,pp)*(1.0+f_g(pz,pp)); //annihilation

    return i;

}

double I_annih(double pp, double pz) //integrand
{
    //double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*(f_g(pz,pp)+f_q(pz,pp));
    //double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*f_g(pz,pp)*(1.0-f_q(pz,pp)); //compton
    double i = (pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*f_q(pz,pp)*(1.0+f_g(pz,pp)); //annihilation

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

    double abs_error = 1.0e-4;  
    double rel_error = 1.0e-4;  
    double result;  
    double error;

    nrfunc=func; 
    gsl_integration_qag (&integrand, pp1, pp2, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);
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

    double abs_error = 1.0e-4;  
    double rel_error = 1.0e-4;  
    double result;  
    double error;

    ppsav=pp;

    gsl_integration_qag (&integrand, -6.0, 6.0, abs_error, rel_error, 1000, 4, work_ptr, &result, &error);
    gsl_integration_workspace_free(work_ptr);

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

    //double Q_s = 1.0;

    double f_0[] = {2.25, 3.81, 5.75, 6.0, 6.65, 7.0, 9.5, 10.25, 11.1, 11.75};
    double xi[] = {1.0, 1.5, 1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5};

    for (int f=9; f<10; f++)
    {
        ofstream myfile1;
        //myfile1.open("thermal_photon_rate_T=180_mu=-0.003.dat");
        myfile1.open("photon_yield_integral_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + ".dat");
        // ofstream myfile2;
        // myfile2.open("thermal_photon_rate_ratios_T=180_mu=-0.003.dat");
        //myfile2.open("thermal_photon_rate_ratios_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + ".dat");

        double pp1 = 0.01;
        double pp2 = 6.01;
        double tau;

        cout << f_0[f] << endl;

        for (int k=0; k<490; k++)
        {
        //for (int k=0; k<21; k++)
        //{
            //pz_const = 0.5;
            //pp_const = 0.2*k;

            double d_tau = 0.1;
            tau = 1.0+(d_tau*k);

            fstream file;
            vector < vector <double> > array; // 2d array as a vector of vectors
            vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
            int row = 0; // Row counter

            file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/f0_" + std::to_string(f_0[f]) + "_xi_" + std::to_string(xi[f]) + "_np_301_nk_64/" + std::to_string(tau) + "_tracking.dat", ios::in);
            //cout << "/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/f0_" + std::to_string(f_0[f]) + "_xi_" + std::to_string(xi[f]) + "_np_301_nk_64/" + std::to_string(tau) + "_tracking.dat" << endl;
            //file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/f0=" + std::to_string(f_0) + "/" + std::to_string(tau) + "_tracking.dat", ios::in);
            //file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/xi=1.5_f0=" + std::to_string(f_0) + "/" + std::to_string(tau) + "_tracking.dat", ios::in);
            //file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/xi=3.7_f0=" + std::to_string(f_0) + "/" + std::to_string(tau) + "_tracking.dat", ios::in);

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

            for (int l=1; l<10201; l++)
            {   
                fq[l]=array[l][3];
                fg[l]=array[l][2];
            }

            //double R_tot = ((40.0*alpha*alphas)/(9.0*pi*pi))*L*f_q(pz_const,pp_const)*quad2d(I_tot, pp1, pp2)*pow(hc,-4.0);
            //double R_comp = ((40.0*alpha*alphas)/(9.0*pi*pi))*L*f_q(pz_const,pp_const)*quad2d(I_comp, pp1, pp2)*pow(hc,-4.0);
            //double R_annih = ((40.0*alpha*alphas)/(9.0*pi*pi))*L*f_q(pz_const,pp_const)*quad2d(I_annih, pp1, pp2)*pow(hc,-4.0);
            double R_tot = quad2d(I_tot, pp1, pp2);
            double R_comp = quad2d(I_comp, pp1, pp2);
            double R_annih = quad2d(I_annih, pp1, pp2);

            //cout << quad2d(I_tot, pp1, pp2) << "\t" << quad2d(I_comp, pp1, pp2) << "\t" << quad2d(I_annih, pp1, pp2) << endl;

            myfile1 << std::scientific << tau << "\t" << R_tot << "\t" << R_comp << "\t" << R_annih << endl;
            //myfile2 << std::scientific << pp_const << "\t" << R_comp/R_tot << "\t" << R_annih/R_tot << endl;
            cout << std::scientific << tau << "\t" << R_tot << "\t" << R_comp << "\t" << R_annih << endl;
            //cout << std::scientific << pp_const << "\t" << R_comp/R_tot << "\t" << R_annih/R_tot <<  endl;
        }   
        myfile1.close();
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
