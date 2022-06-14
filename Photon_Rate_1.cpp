#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

#include <gsl/gsl_integration.h>

# define pi           3.14159265358979323846  /* pi */

//********************************************************************************************************
//Compile using: g++ -Wall -o Photon_Rate_1.out Photon_Rate_1.cpp -lgsl -lgslcblas
//Run: ./Photon_Rate_1.out
//********************************************************************************************************

double E; //GeV, photon energy 

double abs_error;	
double rel_error;

//Defining constants
const double T=0.2; //GeV
const double g = 5.0;
const double alphas = 0.4;
const double alpha = 1.0/137.0;
const double hc = 0.1973269788; //GeV fm

//*******************************************************************************************************
//Function being integrated
double f(double s,double t,double E1,double E2)
{
	double a, b, c, M, u, E2min, I;

    a = -1.0*pow((s+t),2.0);
    b = (s+t)*((E*s)-(E1*t));
    c = (s*t*(s+t))-pow((E*s)+(E1*t), 2.0);

    u = -s-t;

    //Matrix element
    M = ((320.0/3.0)*(16.0/3.0)*pi*pi*alpha*alphas*((((s*s)+(s*t))/(s*s))+(((s*s)+(s*t))/(u*u))))+(20.0*(128.0/9.0)*pi*pi*alpha*alphas*(((t*u)/(t*t))+((t*u)/(u*u))));

    double E2min_1 = (-t)/(4.0*E);
	double E2min_2 = (-b+sqrt(b*b - a*c))/a;
	double E2max = (-b-sqrt(b*b - a*c))/a;

    if( ((b*b) - (a*c)) >= 0.0) //check if points are within phase space 
	{
		if(E2min_1 < E2min_2)
		{
			E2min = E2min_2;
		}
		else  
		{
			E2min = E2min_1;
		}

		if(E2max < E2min) //if the lower bound for E2 is larger than the upper bound, skip
		{
			return 0.0;
		}
		else if ((E1+E2-E) < 0.0)
		{
			return 0.0;
		}
		else 
		{
			//I = ((1.0/pow(2.0*pi,7.0))*(1.0/(16.0*E)))*M*exp(-(E1+E2)/T)*(1.0/sqrt((a*pow(E2,2.0))+(2.0*b*E2)+c))*pow(hc,-4.0); //integrand
			I = ((1.0/pow(2.0*pi,7.0))*(1.0/(16.0*E)))*M*exp(-(E1+E2)/T)*(1.0+exp(-(E1+E2-E)/T))*(1.0/sqrt((a*pow(E2,2.0))+(2.0*b*E2)+c))*pow(hc,-4.0); //integrand
			return I;
		}
	}

	else
	{
		return 0.0; //out of phase space
	}

}
//*******************************************************************************************************


//*******************************************************************************************************
//Integration Limits
double tmin(double s)
{
	double kc = (g*T*T)/3.0; //infrared cutoff
	double t = -s + kc;
	return t;
}

double tmax(double s)
{
	double kc = (g*T*T)/3.0; //infrared cutoff
	double t = -kc;
	return t;
}

double E1min(double s,double t)
{
	double E1 = (s+t)/(4.0*E);
	return E1;
}


double E2min(double s,double t, double E1)
{
	double a, b, c, E2;

    a = -1.0*pow((s+t),2.0);
    b = (s+t)*((E*s)-(E1*t));
    c = (s*t*(s+t))-pow((E*s)+(E1*t), 2.0);

	double E2min_1 = (-t)/(4.0*E);
	double E2min_2 = (-b+sqrt(b*b - a*c))/a;

	//Choses the largest of two posssible lower limits
	if( ((b*b) - (a*c)) >= 0)
	{
		if(E2min_1 < E2min_2)
		{
			E2 = E2min_2;
		}
		else 
		{
			E2 = E2min_1;
		}
	}

	else
	{
		return 0.0;
	}

	return E2;

}

double E2max(double s,double t, double E1)
{
	double a, b, c, E2;

    a = -1.0*pow((s+t),2.0);
    b = (s+t)*((E*s)-(E1*t));
    c = (s*t*(s+t))-pow((E*s)+(E1*t), 2.0);

	E2 = (-b-sqrt(b*b - a*c))/a;

	return E2;
}
//*******************************************************************************************************

//*******************************************************************************************************
//Four dimensional integration routine
static double ssav,tsav,E1sav;
static double (*nrfunc)(double,double,double,double);

double quad4d(double (*func)(double, double, double, double), double smin, double smax)
{
	double f1(double s, void *params);

	gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

	gsl_function integrand;
	integrand.function = &f1;
		
    double result;	
    double error;	

	nrfunc=func;
	gsl_integration_qagiu (&integrand, smin, abs_error, rel_error, 1000, work_ptr, &result, &error);
	return result;
}

double f1(double s, void *params) 
{
	double f2(double t, void *params);
	double tmin(double),tmax(double);

	gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

	gsl_function integrand;
	integrand.function = &f2;

    double result;	
    double error;

	ssav=s;
	//cout << "s: "<< s << endl;
	gsl_integration_qag (&integrand, tmin(s), tmax(s), abs_error, rel_error, 1000, 3, work_ptr, 
	                       	&result, &error);
	return result;
}

double f2(double t, void *params) 
{
	double f3(double E1, void *params);
	double E1min(double,double),E1max(double,double);

	gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

	gsl_function integrand;
	integrand.function = &f3;

    double result;	
    double error;

	tsav=t;
	//cout << "t: " << t << endl;
	gsl_integration_qagiu (&integrand, E1min(ssav,t), abs_error, rel_error, 1000, work_ptr, &result, &error);
	return result;
}

double f3(double E1, void *params)  
{
	double f4(double E2, void *params);
	double E2min(double,double,double),E2max(double,double,double);

	gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

	gsl_function integrand;
	integrand.function = &f4;

    double result;	
    double error;

	E1sav=E1;
	//cout << "E1: " << E1 << endl;

	gsl_integration_qag (&integrand, E2min(ssav,tsav,E1), E2max(ssav,tsav,E1), abs_error, rel_error, 1000, 3, work_ptr, 
	                       	&result, &error);

    return result;

}

double f4(double E2, void *params)  
{
    double f = (*nrfunc)(ssav,tsav,E1sav,E2);
	//cout << "f: "<< f << endl;
	return f;
}

//*******************************************************************************************************


//*******************************************************************************************************
int main ()
{
	std::clock_t start;
    double duration;
    start = std::clock();

    ofstream myfile;
    myfile.open("photon_rate.dat");
		
	cout << "Energy (GeV):    " << "Result (fm^-4 GeV^-2):   " << "Analytic Solution:    " << "Ratio: " << endl;

	//for (int i=0; i<=58; i++)
	for (int i=1; i<=30; i++)
	{
		//E = 0.1000 + (0.0500*i);
		E = 0.1000*i;

        if (E < 2.1)
        {
        	abs_error = 5.0e-7;
            rel_error = 5.0e-7;
        }
        else if (1.95 <= E && E <= 2.05)
        {
        	abs_error = 1.0e-16;
            rel_error = 1.0e-16;
        }
        else if (2.1 <= E && E <= 2.25)
        {
        	abs_error = 5.0e-9;
            rel_error = 5.0e-9;
        }
        else if (2.25 <= E && E < 2.35)
        {
        	abs_error = 2.0e-8;
            rel_error = 2.0e-8;
        }
        else if (E == 2.45)
        {
        	abs_error = 1.0e-8;
            rel_error = 1.0e-8;
        }
        else if (E == 3.00)
        {
        	abs_error = 1.0e-9;
            rel_error = 1.0e-9;
        }
        else if (2.7 <= E && E <= 2.90)
        {
        	abs_error = 2.0e-10;
            rel_error = 2.0e-10;
        }
        else 
        {
        	abs_error = 1.6e-9;
            rel_error = 1.6e-9;
        }

	    double kc = (g*T*T)/3.0; //infrared cutoff
	    double smin = 2.0*kc; 
	    double smax = 3.0; //smax chosen randomly, taken to infinity with change of variables
	    double R;
        R = quad4d(f, smin, smax); 

        double y = (5.0/9.0)*((alpha*alphas)/(2.0*pi*pi))*T*T*exp(-E/T)*log(((2.912*E)/(g*T))+1.0)*pow(hc,-4.0);
        double ratio = y/R;

        cout << std::scientific << E << "     " << R << "             " << y << "          " << ratio <<endl;
        myfile << std::scientific << E << "     " << R << endl;
	}

	myfile.close();

	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"Time: "<< duration <<'\n';

}
//*******************************************************************************************************
