#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 23 15:34:00 2018

@author: menorahlam
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

G = 6.67e-11  # Gravitational Constant (m^3kg^-1s^-2)
c = 3.00e8  # Speed of Light (m/s)
h = 6.63e-34  # Planck's Constant (Js)
mp = 1.67e-27  # Proton mass (kg)
me = 9.11e-31  # Electron Mass (kg)
K_NR = ((h**2)/(5*me))*(3/(8*np.pi))**(2/3)  # EOS Non-Relativistic constant
K_R = ((1/mp)**(4/3)) * (h*c /4) * (3/(8*np.pi))**(1/3)  # EOS Relativistic constant
Ye = 0.5  # Z/A
K = (1.0036e13 * Ye**(5/3))
M_sun = 1.989e30  # Solar Mass (kg)
SN1A_mass = 1.44*M_sun

dr = 100
r0 = 0.
rf = 2455400  # Boundary Radius
h = dr  # Iteration Step


def dm (r, rho):

    return (4*np.pi*r**2)*rho


def dP (m, rho, r):
    return -G*m*rho / (r**2)

def EOS (P):
#    return (P/K_R)**(3/4)  # Relativistic
#    return mp*(P/K_NR)**(3/5)  # Non-Relativistic
    return (P * (1/K))**(3/5)


def P (rho):
#    return K_R*(rho**(4/3))  # Relativistic
#    return K_NR * (rho/mp)**(5/3)  # Non-Relativistic
    return K*rho**(5/3)

def k(h, dy):
    return h*dy

###############################
m_c = 0.  # Define center mass
rho_c = 1e9  # Define center density (This important initial condition varrys)
              # with different models
rho_carray = [1e9]
#rho_c = 1e13

###############################


#r = np.arange(r0, rf + dr, dr)  # Create the r array

# this dy has x=x_n + 1/2h and y=y_n + 1/2k1. Remember to correct this.
#def k2(h, dy):
#    return h*dy

# I introduce the 2nd order RK method for this session for simplicity
def yn (y, k1, k2, k3, k4):
    return y + (1/6)*(k1 + 2*k2 + 2*k3 + k4)

#i = 0
#while i < len(x)-1:
#    ri = r[i]
#    rhoi = rho[i]
#    mi = m[i]
#    k1_m = h*dm(ri, rhoi)
#    k1_P = h*dP(mi, rhoi, ri)
#    k2_m = h*dm(ri + 0.5*h, rhoi + 0.5*h)
#    k2_P = h*dP(mi + 0.5*k1_m, rhoi + 0.5*h, ri + 0.5*h)
#    k3_m = h*dm(ri + 0.5*h, rhoi + 0.5*h)
#    k3_P = h*dP(mi + 0.5*k2_m, rhoi + 0.5*h, ri + 0.5*h)
#    k4 = h*dy(yi + k3, xi + h)
#    y_new = 0.
#    y_new = yn(yi, k1, k2, k3, k4)
#    y_RK.append(y_new)
#    i += 1
#    yi = np.copy(y_new)

M = []
R = []

i = 0
j = 0
for j in range(len(rho_carray)):
    m = []  # Create the mass profile
    p = []  # Create the pressure profile
    rho = []  # Create the density profile
    r = []
    ri = 0
    rhoi = rho_carray[j]
    P_c = P(rhoi)
    mi = m_c
    Pi = P_c 
    while Pi >= 0.:
        #ri = r[i+1]
        #rhoi = rho[i]
        #mi = m[i]
        #Pi = P[i]
        r.append(ri)
        m.append(mi)
        p.append(Pi)
        rho.append(rhoi)
        ri += dr
        m_new = mi + dm(ri, rhoi)*dr
        P_new = Pi + dP(mi, rhoi, ri)*dr
        rho_new = EOS(P_new)
    
        mi = np.copy(m_new)
        Pi = np.copy(P_new)
        rhoi = np.copy(rho_new)
        #    mi = np.copy(m_new)
        #    Pi = np.copy(P_new)
        #    rhoi = np.copy(rho_new)
    M.append(m[-1])
    R.append(r[-1])
#P[:] > 0
#print (r)
#print (P)
#print (P_real)
#print (len(P_real))
plt.plot(R, M)
plt.show()
plt.clf()
    
    
    
plt.plot(r, m, label="Numerical Sol. (FD)")
plt.plot([r[0], r[-1]], [SN1A_mass, SN1A_mass], label="SN1A_WD_Mass")
plt.title("Mass Profile")
plt.xlabel("Radius/R")
plt.ylabel("Mass (kg)")
plt.legend()
plt.show()
plt.clf()

plt.plot(r, p, label="Numerical Sol. (FD)")
plt.title("Pressure Profile")
plt.xlabel("Radius/R")
plt.ylabel("Pressure (Pa)")
plt.legend()
plt.show()
plt.clf()

plt.plot(r, rho, label="Numerical Sol. (FD)")
plt.title("Density Profile")
plt.xlabel("Radius/R")
plt.ylabel("Density (kgm^-3)")
plt.legend()
plt.show()
plt.clf()


print (G*m[-1]/(c**2 * r[-1]))