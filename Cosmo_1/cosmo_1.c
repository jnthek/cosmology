#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_integration.h>

#define K 1.38064852e-23
#define C 299792458.0
#define H 6.62607004e-34
#define PI 3.14
#define T_CMB 2.7

double number_density(double nu, void *params)
{
	double nn,T;
	double z = *(double*) params;
	T=T_CMB * (z + 1);
	nn = ((8*PI*nu*nu)/(C*C*C))*(1/(exp(H*nu/(K*T))-1));
	return nn;	
}

int main()
{
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000);
	double n,error,freq_low,freq_high;
	double z;
	z = 1000;
	
	gsl_function P_s;
	P_s.function = &number_density;	
	P_s.params = &z;
	
	freq_low = 0;
	freq_high = 1e20;
	
	//gsl_integration_qagiu (&P_s,0,0,1e-1,100000,w,&n,&error);	
	//gsl_integration_qags(&P_s,freq_low,freq_high,0,1e-13,100000,w,&n,&error);
	gsl_integration_qag(&P_s,freq_low,freq_high,0,1e-13,100000,GSL_INTEG_GAUSS61,w,&n,&error);
	
	printf("%.15e\n",n);
	return 0;
}
	
