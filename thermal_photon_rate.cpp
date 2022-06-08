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

using namespace std;

//**************************************************************************************************************
//Compile using: g++ -Wall -o thermal_photon_rate.out thermal_photon_rate.cpp -lgsl -lgslcblas
//**************************************************************************************************************


//**************************************************************************************************************
#define pi 3.14159265358979323846

double pp;
const double pz = 0.5; //GeV
const double L = 0.75;
const double hc = 0.1973269788; //GeV fm
const double alphas = 0.4;
const double alpha = 1.0/137.0;
const double g = sqrt(5.0);
double T;
double mu;


double therm(double pp, double pz) //integrand
{
    double E = sqrt((pp*pp)+(pz*pz));

    double i = (6.0/5.0)*(5.0/9.0)*((alpha*alphas)/(2.0*pi*pi))*T*T*exp(-(E-(2.0*mu))/T)*log(((2.912*E)/(g*g*T))+1.0)*pow(hc,-4.0);

    return i;
}

double small_angle(double pp, double pz) //integrand
{
    double E = sqrt((pp*pp)+(pz*pz));

    double i = (6.0/5.0)*(5.0/(9.0*pi*pi))*alpha*alphas*T*T*L*exp(-(E-(2.0*mu))/T)*pow(hc,-4.0);

    return i;

}
//**************************************************************************************************************


//************************************************************************************************************
//main function
//************************************************************************************************************

//#define COLS 8 // Number of columns in data
int main (int argc, char **argv)
{
    std::clock_t start;
    double duration;
    start = std::clock();

    double Q_s = 1.0;
    
    // cout << "Enter number between 1-2 as a value for Q_s (GeV): " << endl;
    // cin >> Q_s;
    // cout << "Value of " << Q_s << "GeV chosen for Q_s." << endl;

    //Specify temperature and chemical potential values 
    // T = 0.218*Q_s; //Gev      
    // mu = -0.29*Q_s; //GeV

    // T = 0.24846604739923306; //Li
    // mu = 0.0;

    // T = 3.088201e-01; //Jess
    // mu = -1.149946e-01;

    // T = 0.21725690831702282; //Li
    // mu = -0.3514673634298637;

    //T = 2.176008e-01; //Jess
    //mu = -3.738523e-01;

// ***************************************************
    //T = 0.241905924775; //f_0=0.45, t=2.0
    //mu = 0.0; 

    //T = 0.20493822102; //f_0=0.45, t=5.0
    //mu = -0.29725;

    //T = 0.184111601759; //f_0=0.45, t=10.0
    //mu = -0.66875;

    //T = 0.16076603192; //f_0=0.45, t=25.0
    //mu = -1.18275;



    //T = 0.24077829111244176; //f_0=0.30, t=2.0
    //mu = -0.06482955488202494;

    //T = 0.21045784147819085; //f_0=0.30, t=5.0
    //mu = -0.1647358754170539;

    //T = 0.19357153422781972; //f_0=0.30, t=10.0
    //mu = -0.23785102268243352;

    // T = 0.1739638621922704; //f_0=0.30, t=25.0
    // mu = -0.3183103768463068;

    T = 0.2;
    mu = 0.0;





    cout << "p_perp:         " << "E dN/d^4xd^3p:     " << endl;

    ofstream myfile;
    myfile.open("thermal_photon_rate_T=200_mu=0.dat");

    for (int i=0; i<=29; i++)
    {
        pp = 0.1000*i;
        
        cout << std::scientific << pp << "\t" << therm(pp, pz)*pow(Q_s, 4.0) << "\t" << small_angle(pp, pz)*pow(Q_s, 4.0) << endl;
        myfile << std::scientific << pp << "\t" << therm(pp, pz)*pow(Q_s, 4.0) << "\t" << small_angle(pp, pz)*pow(Q_s, 4.0) << endl;
    }
    myfile.close();

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';
}