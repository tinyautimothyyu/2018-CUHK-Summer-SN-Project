import numpy as np
import matplotlib.pyplot as plt

# Let theta be y, phi be z, and epsilon be x for the Lane–Emden equation
# Theta and phi are dimensionless variables, rho=rho_c*theta**n and r=a*epsilon
M_sun = 1.989e30  # Solar Mass (kg)
SN1A_mass = 1.44*M_sun
#rho_c = 1e9  # Central density (SI unit)
rho_c = 1e9 * 0.001  # Central density (cgs)
G = 6.67e-11  # Gravitational Constant (m^3kg^-1s^-2)
G_cgs = 6.67e-8  # Gravitational Constant cgs
c = 3.00e8  # Speed of Light (m/s)
h = 6.63e-34  # Planck's Constant (Js)
mp = 1.67e-27  # Proton mass (kg)
me = 9.11e-31  # Electron Mass (kg)
Ye = 0.5  # Z/A
K = 1.0036e13 * Ye**(5/3)  # Non-Relativistic EOS Constant
K_R = 1.2435e15 * Ye**(4/3)  # Relativistic EOS Constant
n = 3/2  # Non-relativistic constant for EOS
n_R = 3  # Relativistic Constant for EOS
square_a = (n+1)*K*rho_c**(1/n - 1)/(4*np.pi*G_cgs)
square_a_R = (n_R+1)*K_R*rho_c**(1/n_R - 1)/(4*np.pi*G_cgs)

dh = [(1e-5)/7, (1e-5)/6, 1e-5, 1e-3, 1e-1, 1.]  # Define iteration step size

# Define B.C for the ODEs
x_0 = 1e-10
y_0 = 1.
z_0 = 0.

# x, y, z arrays
x = [x_0]
y = [y_0]
z = [z_0]

# Define the equation that perform dimensionless variable theta
def rho (theta):
    return rho_c*theta**n
#    return rho_c*theta**n_R


# Define the equation that perform dimensionless variable epsilon
def r (epsilon):
    a = square_a**(1/2)
#    a = square_a_R**(1/2)
    return a*epsilon
    

def m (phi):
    return 4*np.pi*(square_a**(3/2)) * rho_c *phi
#    return 4*np.pi*(square_a_R**(3/2)) * rho_c *phi


def P (rho):
    return K*rho**(5/3)
#    return K_R*rho**(5/3)

# Define dtheta/depsilon = f(x, y, z)
def f (x, y, z):
    return -(z/(x**2))


# Define dphi/depsilon = g(x, y, z)
def g (x, y, z):
    return (y**n) * x**2
#    return (y**n_R) * x**2


# Main loop that get profiles
def main (x, y, z, dh):
    i = 0
    while y[i] >= 0:
        # RK-4 Method
        k1_f = dh*f(x[i], y[i], z[i])
        k1_g = dh*g(x[i], y[i], z[i])
        k2_f = dh*f(x[i]+(dh/2), y[i]+(k1_f/2), z[i]+(k1_g/2))
        k2_g = dh*g(x[i]+(dh/2), y[i]+(k1_f/2), z[i]+(k1_g/2))
        k3_f = dh*f(x[i]+(dh/2), y[i]+(k2_f/2), z[i]+(k2_g/2))
        k3_g = dh*g(x[i]+(dh/2), y[i]+(k2_f/2), z[i]+(k2_g/2))
        k4_f = dh*f(x[i]+(dh), y[i]+(k3_f), z[i]+(k3_g))
        k4_g = dh*g(x[i]+(dh), y[i]+(k3_f), z[i]+(k3_g)) 
        
        x.append(x[i]+dh)
        y.append(y[i]+(1/6)*(k1_f+k2_f*2+k3_f*2+k4_f))
        z.append(z[i]+(1/6)*(k1_g+k2_g*2+k3_g*2+k4_g))
        i += 1
    return x, y, z

# Convergence Test
#j=0
#theta_0 = []
#for j in range(len(dh)):
#    x = [x_0]
#    y = [y_0]
#    z = [z_0]
#    print(j)
#    print (dh[j])
#    xf, yf, zf = [], [], []
#    xf, yf, zf = main(x, y, z, dh[j])
#    xf_array = np.copy(np.array(xf))
#    print(xf_array.size)
#    yf_array = np.copy(np.array(yf))
#    print (yf_array.size)
#    print (yf_array[np.abs(xf_array[::] - 3.) <= 1e-8])
#    if j == 0:
#        thetai = np.copy(yf_array[np.abs(xf_array[::] - 3.) <= 1e-8][0])
#    else:
#        theta_0.append(yf_array[np.abs(xf_array[::] - 3.) <= 1e-8][0])

#theta_array = np.array(theta_0)
#dh_array = np.array(dh)
#print (theta_array)
#print (theta_array.size)


j = 1
xf, yf, zf = main(x, y, z, dh[j])

#Convergence Test Plot#
#plt.plot(np.log(dh_array[1::]), np.log(np.abs(theta_array - thetai)))
#plt.xlabel("ln(dr)")
#plt.ylabel("ln(error)")
#plt.title("Convergence Test")
#plt.savefig("ConvergenceTest")
#plt.clf()

plt.plot(xf, yf, label="theta-epsilon")
plt.xlabel("epsilon")
plt.ylabel("theta")
plt.legend()
plt.show()
plt.clf()

plt.plot(xf, zf, label="phi-epsilon")
plt.xlabel("epsilon")
plt.ylabel("phi")
plt.legend()
plt.show()

rho_array = rho(np.array(yf))*1000
m_array = m(np.array(zf))/1000
P_array = P(rho(np.array(yf)))*0.1

#Density Profile
plt.plot(np.array(xf)/xf[-1], rho_array)
plt.xlabel("radius/R*")
plt.ylabel("Density (kgm^-3)")
plt.title("Density Profile")
plt.savefig("Density Profile_R")
plt.clf()

# Mass Profile
plt.plot(np.array(xf)/xf[-1], m_array/M_sun, label="RK-4")
plt.xlabel("radius/R*")
#plt.ylabel("Mass (Solar Mass)")
plt.title("Mass Profile")
plt.savefig("Mass Profile_R")
plt.clf()

# Pressure Profile
plt.plot(np.array(xf)/xf[-1], P_array, label="RK-4")
plt.xlabel("radius/R*")
plt.ylabel("Pressure (Pa)")
plt.title("Pressure Profile")
plt.savefig("Pressure Profile_R")
plt.clf()

h=open("WD_R Parameters.txt", "w")
h.write("Epsilon = {}\n".format(xf[-1]))
h.write("Radius* = {} km\n".format(r(xf[-1])/100000))
h.write("Mass* = {} SolarMass\n".format(m_array[-1]/M_sun))
h.close()


print ("Epsilon = {}".format(xf[-1]))
print ("Radius* = {} km".format(r(xf[-1])/100000))
print (1.12e4*(rho_c/1e6)**(-1/6))


