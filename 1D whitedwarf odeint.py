import scipy as sp
import numpy as np
import scipy.misc
from scipy.integrate import odeint
import matplotlib.pyplot as plt

k = 0.1 * 1.0036 * (1e13) * 0.5 ** (5 / 3)
def function(z,r):
    G = 6.674e-11
    k = 0.1*1.0036 * (1e13) * 0.5 ** (5 / 3)
    central_rho = 1e9


    rho = z[0]
    m = z[1]
    drhodr = (-G*m*rho)/(k*(5/3)*r*r*(rho)**(2/3))
    dmdr = 4*np.pi*r*r*rho
    return [drhodr, dmdr]

r = np.linspace(0.1,5e10,100000)
diu = odeint(function,[1e9, 0], r)


P = k*(diu[:,0])**(5/3)

plt.plot(r,diu[:,0])
plt.ylabel("density (m)")
plt.xlabel('distance from core (m)')
plt.show()

plt.plot(r,diu[:,1])
plt.ylabel("mass (kg)")
plt.xlabel('distance from core (m)')
plt.show()

plt.plot(r,P)
plt.ylabel("Pressure")
plt.xlabel('distance from core (m)')
plt.show()