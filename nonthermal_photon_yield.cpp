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
//Compile using: g++ -Wall -o nonthermal_photon_yield.out nonthermal_photon_yield.cpp -lgsl -lgslcblas
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

    return f;
}

//**************************************************************************************************************


//************************************************************************************************************
//main function
//************************************************************************************************************

#define pi 3.14159265358979323846
const double hbarc = 0.1973269788;
#define COLS 4 // Number of columns in data
#define COLS2 4
double pz_const;
double pp_const;
double pz_tilde_const;
double eta;
int start;
int stop;
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    //double Q_s = 1.0;
    // double Q_s = 2.0;



    //New Runs
    double f_0[] = {2.25, 3.81, 5.75, 6.65, 9.5, 11.1, 6.0, 7.0, 10.25, 11.75};

    double xi[] = {1.0, 1.5, 1.0, 1.0, 1.5, 1.5, 1.0, 1.0, 1.5, 1.5};

    double Q_s[] = {1.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0};

    for (int f=0; f<10; f++)
    {
        fstream file2;
        vector < vector <double> > integral; // 2d array as a vector of vectors
        vector <double> rowVector2(COLS2); // vector to add into 'array' (represents a row)
        int row2 = 0; // Row counter

        file2.open("photon_yield_integral_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + ".dat");
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

        ofstream myfile1;

        //myfile1.open("nonthermal_photon_yield_0.1-0.6_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + ".dat");
        myfile1.open("nonthermal_photon_yield_20-40_0.2-0.4_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + ".dat");
        //myfile1.open("nonthermal_photon_yield_0.4_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + "_L=0.75.dat");
        //myfile1.open("nonthermal_photon_yield_0.1-1.0_xi=" + std::to_string(xi[f]) + "_f0=" + std::to_string(f_0[f]) + ".dat");

        cout << "f_0 = " << f_0[f] << endl; 
        cout << "xi = " << xi[f] << endl; 
        cout << "Q_s = " << Q_s[f] << endl; 

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
        cout << "start = " << start << endl;
        cout << "stop = " << stop << endl; 

        for (int i=3; i<=20; i++)
        {
            pp_const = 0.2000*i;
            pz_const = 0.5;

            double tau;
            double sum_tot = 0.0;
            double sum_comp = 0.0;
            double sum_annih = 0.0;
            double R_coeff = 0.0;
            double R_tot = 0.0;
            double R_comp = 0.0;
            double R_annih = 0.0;

            //for (int k=0; k<491; k++) 10fm
            //for (int k=0; k<31; k++) //0.4fm/c for Q_s = 1 GeV
            //for (int k=0; k<11; k++) //0.4fm/c for Q_s = 1 GeV

            //for (int k=0; k<51; k++) //0.6fm/c for Q_s = 2 GeV
            //for (int k=0; k<91; k++) //1.0fm/c for Q_s = 2 GeV


            for (int k=start; k<stop; k++) //1.0fm/c for Q_s = 1 GeV
            {
                double d_tau = 0.1;
                tau = 1.0+(d_tau*k);

                fstream file;
                vector < vector <double> > array; // 2d array as a vector of vectors
                vector <double> rowVector(COLS); // vector to add into 'array' (represents a row)
                int row = 0; // Row counter

                file.open("/Users/JessicaChurchill/Documents/Research/Final_Code/Paper_Code/f0_" + std::to_string(f_0[f]) + "_xi_" + std::to_string(xi[f]) + "_np_301_nk_64/" + std::to_string(tau) + "_tracking.dat", ios::in);

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

                for(int m=0; m<=100; m++)
                {
                    double d_eta = 0.1;
                    eta = d_eta*m;

                    pz_tilde_const = pp_const*sinh(asinh(pz_const/pp_const)-eta);

                    double A; 

                    if (f<2)
                    {
                        //A = 100.58; //fm^2 for RHIC 0-20%
                        A = 57.75; //fm^2 for RHIC 20-40%
                    }

                    else if (f>=2 || f<6)
                    {
                        //A = 124.25; //fm^2 for LHC 0-20%
                        A = 75.9; //fm^2 for LHC 20-40%
                    }
                    
                    else
                    {
                        //A = 127.75; //fm^2 for LHC 0-20% @ 5.02 TeV
                        A = 79.6; //fm^2 for LHC 20-40% @ 5.02 TeV
                    }

                    double L = 0.75;
                    //double L = 1.40;
                    //double alphas = 0.4;
                    const double alphas = 0.234;
                    double alpha = 1.0/137.0;

                    R_coeff = 2.0*A*(d_tau/Q_s[f])*(tau/Q_s[f])*d_eta*(40.0/(9.0*pi*pi))*alpha*alphas*L*f_q(pz_tilde_const,pp_const)*(6.0/5.0)*pow(hbarc,-4.0)*pow(Q_s[f],2.0);
                    R_tot = R_coeff*integral[k][1];
                    R_comp = R_coeff*integral[k][2];
                    R_annih = R_coeff*integral[k][3];
                    sum_tot = sum_tot + R_tot;
                    sum_comp = sum_comp + R_comp;
                    sum_annih = sum_annih + R_annih;

                }
            }
            myfile1 << std::scientific << pp_const << "\t" << sum_tot/pi << "\t" << sum_comp/pi << "\t" << sum_annih/pi << "\n";
            cout << std::scientific << pp_const << "\t" << sum_tot/pi << "\t" << sum_comp/pi << "\t" << sum_annih/pi << "\n";
        }
        myfile1.close();
    }
    

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}
