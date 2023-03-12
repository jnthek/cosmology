#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_integration.h>

#define K 1.38064852e-16
#define C 2.997e10
#define H 6.62607004e-27
#define PI 3.14
#define T_CMB 2.7
#define MP 1.673e-24
#define G 6.674e-8

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
	double n, error, freq_low, freq_high ;
	double z, H0, rho_c, nH0, omega_b, h, nH;
	
	z = 1000;
	h = 0.6774;
	H0 = h * 3.24 * 1e-18;
	omega_b = 0.02207/(h*h);
	
	rho_c = (3*H0*H0)/(8*PI*G);
	nH0 = 0.75*(omega_b*rho_c/MP);
	nH = nH0*(1+z)*(1+z)*(1+z);
		
	gsl_function P_s;
	P_s.function = &number_density;	
	P_s.params = &z;
	
	freq_low = 0;
	freq_high = 1e15;
	
	//gsl_integration_qagiu (&P_s,0,0,1e-1,100000,w,&n,&error);	
	//gsl_integration_qags(&P_s,freq_low,freq_high,0,1e-13,100000,w,&n,&error);
	gsl_integration_qag(&P_s,freq_low,freq_high,0,1e-13,100000,GSL_INTEG_GAUSS61,w,&n,&error);
		
	printf("\nAt z = %5.2f, we have\n",z);
	printf("Radiation Density is %.15e cm^-3\n",n*(1+1000)*(1+1000)*(1+1000));
	printf("Hydrogen density is %e cm^-3\n",nH);
	printf("Radiation to matter density ratio is %e\n",n/nH);
	return 0;
}
	
