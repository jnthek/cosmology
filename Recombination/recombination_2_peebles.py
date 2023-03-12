# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 14:00:07 2017

@author: jishnu
"""

'''Looks like this code works'''

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

pi = np.pi
kB = 1.4e-16
c = 3e10
h = 6.6e-27
T_CMB_0 = 2.725
sigma = 5.7e-5
G = 6.7e-8
e = 4.8e-10
c = 3e10
me = 9.1e-28
mp = 1.7e-24
e0 = 13.6*1.6e-12 #ionization energy for Hydrogen

H0 = 2.3e-18
omega_m = 0.27
omega_b = 0.044
omega_r = 8.2e-5
omega_l = 1.0-omega_m-omega_r
rho_c = 3*(H0**2)/(8*pi*G)
rho_b_0 = omega_b*rho_c

lambda_2s1s = 8.227
alpha = 0.007297

def Tbf(a):
    return T_CMB_0/a

def H(a):
    return H0 * np.sqrt((omega_b * a**-3) + (omega_r * a**-4) + omega_l)

def rho_b(a):
    return rho_b_0/(a**3)
    
def nH(a):
    return rho_b(a)/mp

def n1s(Xe, a):
    return (1-Xe)*nH(a)

def lambda_alpha(Xe, a):
    #return H(a) * (3*e0)**3 / ((8*pi)**2 * n1s(Xe, a))
    return (H(a)/a) * (8*pi)/(n1s(Xe,a) * (1.216e-5)**3) 

def alpha2(Tb):
    return 9.78 * (e**4 / (me**2 * c**3)) * (e0/(kB*Tb))**0.5 * np.log(e0/(kB*Tb))
    #return 2.6e-13 * (Tb/1e4)**-0.8
    
def beta(Tb):
    return alpha2(Tb) * (2*pi*me*kB*Tb/h**2)**1.5 * np.exp(-e0/(kB*Tb))

def beta2(Tb):
    return alpha2(Tb) * (2*pi*me*kB*Tb/h**2)**1.5 * np.exp(((h*2.5e15)-e0)/(kB*Tb))
    #return beta(Tb) * np.exp(h*2.5e15/(kB*Tb))

def Cr(Tb, Xe, a):
    return (lambda_2s1s + lambda_alpha(Xe, a)) / (lambda_2s1s + lambda_alpha(Xe, a) + beta2(Tb))

def dy_dx(Xe, a):
    T = Tbf(a)
    value = (1/a) * (Cr(T, Xe, a)/H(a)) * (beta(T)*(1-Xe) - nH(a)*alpha2(T)*Xe**2)  
    return value
 
as1 = np.linspace((1/1587.4), 1/20.0, 10000)
y0 = 0.99 # the initial condition
ys = integrate.odeint(dy_dx, y0, as1)

z = (1/as1)-1

print ys[-1]

plt.semilogy(z, ys, label='Ionization Fraction')
plt.xlabel("Redshift")
plt.ylabel("Fraction Xe")
plt.legend()
plt.gca().invert_xaxis()
plt.show()