# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 21:18:15 2017

@author: jishnu
"""

import scipy.constants as cst
import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt

c=cst.speed_of_light
k=cst.Boltzmann
h=cst.Planck
T_cmb=2.7
numb = np.array([])

for z in range(1000): 
    
    T = T_cmb*(1+z)
    u = lambda nu:((8*np.pi*nu*nu/(c**3))*1/(np.exp(h*nu/(k*T))))
    n=integ.quad(u,2.14e15,2.16e15)  
    numb = np.append(numb,n[0])
    
plt.plot(numb)
plt.show()    
print "Photon density was",n[0],"m^-3"