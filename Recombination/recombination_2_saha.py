# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 14:19:55 2017

@author: jishnu
"""

''' Tries Saha equation for Helium recombination He++ to He+
Cosmological parameters from WMAP first year results'''

import numpy as np
import matplotlib.pyplot as plt

pi = np.pi
kB = 1.4e-16
c = 3e10
h = 6.6e-27
T_CMB_0 = 2.725
sigma = 5.7e-5
G = 6.7e-8
me = 9.1e-28
mp = 1.7e-24
e0 = 54*1.6e-12 #ionization energy for Helium
alpha = 0.007297

H0 = 2.3e-18
omega_m = 0.27
omega_b = 0.044
omega_r = 8.2e-5
omega_l = 1.0-omega_m-omega_r
rho_c = 3*(H0**2)/(8*pi*G)

z = np.linspace(0,8000,1000)
a = 1/(1+z)

rho_b_0 = omega_b*rho_c

T_CMB = T_CMB_0*(1+z)
rho_b = rho_b_0*(a**-3)
H = H0*np.sqrt((omega_b*(a**-3)) + (omega_r*(a**-4)) + omega_l)
nb = rho_b*0.22/(2*mp) #assuming 22 percent helium
Tb = T_CMB #assuming equilibrium conditions during recombination

'''Saha recombination'''
y = (1/nb) * (2*pi*me*kB*Tb/(h**2))**1.5 * np.exp(-e0/(kB*Tb)) #just a handle for RHS of Saha equation
Xe = (-y + np.sqrt(y**2 + 4*y))/2

plt.figure()
#plt.semilogy(z,Xe)
plt.plot(z,Xe)
plt.xlabel('Redshift')
plt.ylabel('Ionisation Fraction')

plt.figure()
plt.semilogx(a,Xe)
plt.xlabel('Scale factor a')
plt.ylabel('Ionisation Fraction')

plt.figure()
plt.loglog(a,H)

plt.show()

#Tb = T_CMB_0*(1+1587.4)
#a = 1/(1+1587.4) 
#rho_b = rho_b_0*(a**-3)
#nb = rho_b/mp
#y = (1/nb) * (2*pi*me*kB*Tb/(h**2))**1.5 * np.exp(-e0/(kB*Tb)) #just a handle for RHS of Saha equation
#Xe = (-y + np.sqrt(y**2 + 4*y))/2
#print Xe



