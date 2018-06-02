import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

G = 6.67e-11  # Gravitational Constant (m^3kg^-1s^-2)
G_cgs = 6.67e-8  # Gravitational Constant cgs
c = 3.00e8  # Speed of Light (m/s)
h = 6.63e-34  # Planck's Constant (Js)
mp = 1.67e-27  # Proton mass (kg)
me = 9.11e-31  # Electron Mass (kg)
Ye = 0.5  # Z/A
K = (1.0036e13 * Ye**(5/3))
K_NR = ((h**2)/(5*me))*(3/(8*np.pi))**(2/3)  # EOS Non-Relativistic constant
K_R = ((1/mp)**(4/3)) * (h*c /4) * (3/(8*np.pi))**(1/3)  # EOS Relativistic constant
M_sun = 1.989e30  # Solar Mass (kg)
SN1A_mass = 1.44*M_sun

dr = 0.5  # K_NR and K_R
#dr = 100000  # K
h = dr  # Iteration Step

def dm (r, rho):
    return (4*np.pi*r**2)*rho


def dP (m, rho, r):
    return -G*m*rho / (r**2)
#    return -G_cgs*m*rho / (r**2)  # K


def EOS (P):
#    return (P * (1/K))**(3/5)
#    return (P/K_R)**(3/4)  # Relativistic
    return mp*(P/K_NR)**(3/5)  # Non-Relativistic


def P (rho):
#    return K*rho**(5/3)
#    return K_R*(rho**(4/3))  # Relativistic
    return K_NR * (rho/mp)**(5/3)  # Non-Relativistic


rho_c = 1e9
m_c = 1e-5

m = []  # Create the mass profile
p = []  # Create the pressure profile
rho = []  # Create the density profile
r = []

ri = 1e-5
rhoi = rho_c
P_c = P(rhoi)
mi = m_c
Pi = P_c 


while Pi >= 0.:
    r.append(ri)
    m.append(mi)
    p.append(Pi)
    rho.append(rhoi)
    m_new = mi + dm(ri, rhoi)*dr
    P_new = Pi + dP(mi, rhoi, ri)*dr
    rho_new = EOS(P_new)

    mi = np.copy(m_new)
    Pi = np.copy(P_new)
    rhoi = np.copy(rho_new)
    ri += dr

    
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