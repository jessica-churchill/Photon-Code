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
//Compile using: g++ -Wall -o thermal_photon_yield.out thermal_photon_yield.cpp -lgsl -lgslcblas
//**************************************************************************************************************

//**************************************************************************************************************
#define pi 3.14159265358979323846
const double L = 0.75;
const double hbarc = 0.1973269788; //GeV fm
const double alphas = 0.234;
const double alpha = 1.0/137.0;
double pp_const;
double pz_const;
double pz_tilde_const;
double eta;
double T, mu;

double f_q(double pz, double pp) //distribution function of quarks: fq
{
    double E = sqrt((pp*pp)+(pz*pz));
    //double f = 1.0/(exp(E-mu)/T);
    //double f = exp(-sqrt((pp*pp)+(pz*pz))/T)*exp(mu/T);
    //double f = exp(-E/T)*exp(mu/T);
    double f = 1.0/((exp(E-mu)/T)+1.0);

    return f;
}

double f_g(double pz, double pp) //distribution function of gluons: fg
{
    double E = sqrt((pp*pp)+(pz*pz));
    //double f = 1.0/(exp(E-mu)/T);
    //double f = exp(-E/T)*exp(mu/T);
    double f = 1.0/((exp(E-mu)/T)-1.0);

    return f;
}

double I(double pp, double pz) //integrand
{
    double i = f_q(pz_tilde_const,pp_const)*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*(f_g(pz,pp)+f_q(pz,pp));
    //double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*exp(-sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))/T)
                //*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*exp(-sqrt((pp*pp)+(pz*pz))/T)*(1.0+exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-4.0); //annihilation
    //double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*exp(-sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))/T)
                //*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*exp(-sqrt((pp*pp)+(pz*pz))/T)*(1.0-exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-4.0); //compton
                

    return i;

}

double I_comp(double pp, double pz) //integrand
{
    //double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*exp(-sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))/T)
    //            *(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*(2.0*exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-4.0);
    //double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*exp(-sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))/T)
                //*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*exp(-sqrt((pp*pp)+(pz*pz))/T)*(1.0+exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-4.0); //annihilation
    double i = f_q(pz_tilde_const,pp_const)*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*f_g(pz,pp)*(1.0-f_q(pz,pp)); //compton
                

    return i;

}

double I_annih(double pp, double pz) //integrand
{
    //double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*exp(-sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))/T)
      //          *(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*(2.0*exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-4.0);
    double i = f_q(pz_tilde_const,pp_const)*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*f_q(pz,pp)*(1.0+f_g(pz,pp)); //annihilation
    //double i = (40.0/(9.0*pi*pi))*alpha*alphas*L*exp(-sqrt((pp_const*pp_const)+(pz_tilde_const*pz_tilde_const))/T)
                //*(pp/pow(2.0*pi,2.0))*(1.0/sqrt((pp*pp)+(pz*pz)))*exp(-sqrt((pp*pp)+(pz*pz))/T)*(1.0-exp(-sqrt((pp*pp)+(pz*pz))/T))*exp((2.0*mu)/T)*pow(hbarc,-4.0); //compton
                

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

    double abs_error = 1.0e-8;  
    double rel_error = 1.0e-8;  
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

#define COLS 12 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start_clock;
    double duration;
    start_clock = std::clock();

    double f_0[] = {2.25, 3.81, 5.75, 6.65, 9.5, 11.1, 6.0, 7.0, 10.25, 11.75};

    double xi[] = {1.0, 1.5, 1.0, 1.0, 1.5, 1.5, 1.0, 1.0, 1.5, 1.5};

    double Q_s[] = {1.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0};

    int start;
    int stop;

    for (int f=0; f<10; f++)
    {
        fstream file;
        vector < vector <double> > energy_density; // 2d array as a vector of vectors
        vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
        int row = 0; // Row counter
        file.open("f0="+ std::to_string(f_0[f]) + ".dat");
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

        cout << f_0[f] << endl;

        ofstream myfile1;
        //myfile1.open("thermal_photon_yield_0.2-1.0_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + ".dat");
        myfile1.open("thermal_photon_yield_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + ".dat");

        if (Q_s[f]==1.0)
        {
            start = 0;
            stop = 11;
        }
        else 
        {
            start = 0;
            stop = 31;
        }

        for (int i=3; i<=20; i++)
        {
            pp_const = 0.2*i;
            pz_const = 0.5;

            double sum_tot = 0.0;
            double sum_comp = 0.0;
            double sum_annih = 0.0;
            double R_coeff = 0.0;
            double R_tot = 0.0;
            double R_comp = 0.0;
            double R_annih = 0.0;

            //for (int k=0; k<91; k++) //1.0fm/c for Q_s = 2 GeV
            for (int k=start; k<stop; k++) //1.0fm/c for Q_s = 1 GeV
            {
                double d_tau = 0.1;
                double tau = 1.0+(d_tau*k);

                for(int m=0; m<=100; m++)
                {
                    double d_eta = 0.1;
                    eta = d_eta*m;

                    pz_tilde_const = pp_const*sinh(asinh(pz_const/pp_const)-eta);

                    double A; 

                    if (f<3)
                    {
                        A = 100.58; //fm^2 for RHIC 0-20%
                        //A = 57.75; //fm^2 for RHIC 20-40%
                    }

                    else if (f>=3 || f<6)
                    {
                        A = 124.25; //fm^2 for LHC 0-20%
                    }
                    
                    else
                    {
                        A = 127.75; //fm^2 for LHC 0-20% @ 5.02 TeV
                    }

                    // if (f<3)
                    // {
                    //     A = 124.25; //fm^2 for LHC 0-20%
                    // }
                    
                    // else
                    // {
                    //     A = 127.75; //fm^2 for LHC 0-20% @ 5.02 TeV
                    // }

                    double pp1 = 0.010;
                    double pp2 = 6.010;

                    //tau = 0.5 + energy_density[k*10][0];
                    T = energy_density[1.0+(k*10)][1];
                    mu = energy_density[1.0+(k*10)][2];

                    R_coeff = 2.0*A*(40.0/(9.0*pi*pi))*alpha*alphas*L*(d_tau/Q_s[f])*(tau/Q_s[f])*d_eta*(6.0/5.0)*pow(Q_s[f],2.0)*pow(hbarc,-4.0);
                    R_tot = R_coeff*quad2d(I, pp1, pp2);
                    R_comp = R_coeff*quad2d(I_comp, pp1, pp2);
                    R_annih = R_coeff*quad2d(I_annih, pp1, pp2);
                    sum_tot = sum_tot + R_tot;
                    sum_comp = sum_comp + R_comp;
                    sum_annih = sum_annih + R_annih;
                    //cout << tau << "\t"<< eta << "\t"<< R_tot << "\t" << R_comp << "\t" << R_annih << endl;
                }
                //cout << tau << "\t" << T << "\t" << mu << endl;
            }
            myfile1 << std::scientific << pp_const << "\t" << sum_tot/pi << "\t" << sum_comp/pi << "\t" << sum_annih/pi << "\n";
            cout << std::scientific << pp_const << "\t" << sum_tot/pi << "\t" << sum_comp/pi << "\t" << sum_annih/pi << "\n";
        }
        myfile1.close();
    }

    duration = ( std::clock() - start_clock ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
