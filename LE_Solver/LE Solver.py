import numpy as np
import matplotlib.pyplot as plt

n = [0, 1, 5]

# Define B.C for the ODEs #
x_0 = 1e-10
y_0 = 1.
z_0 = 0.

# Definition
num = 2000  # number of grids
tot_epsilon = 10  # total length

x = np.linspace (x_0, tot_epsilon, num)  # Radius array in code unit
y = np.zeros(num)  # Density array in code unit
y[0] = y_0
z = np.zeros(num)  # Mass array in code unit

dh = x[-1] / num


# Define dtheta/depsilon = f(x, y, z)
def f (x, y, z):
    return -(z/(x**2))


# Define dphi/depsilon = g(x, y, z)
def g (x, y, z, n):
    return (y**n) * x**2
    

# Analytic Soltion when n = 0
def theta0 (epsilon):
    return 1 - (1/6)*epsilon**2


# Analytic Solution when n = 1
def theta1 (epsilon):
    return np.sin(epsilon) / epsilon


# Analytic Solution when n = 5
def theta5 (epsilon):
    return (1 + epsilon**2 / 3)**(-1/2)


# LE Solver
def LE (n):
    # RK-4th order update x,y,z arrays #
    i = 0
    while i < len(x)-1:
        # RK-4 Method
        k1_f = dh*f(x[i], y[i], z[i])
        k1_g = dh*g(x[i], y[i], z[i], n)
        k2_f = dh*f(x[i]+(dh/2), y[i]+(k1_f/2), z[i]+(k1_g/2))
        k2_g = dh*g(x[i]+(dh/2), y[i]+(k1_f/2), z[i]+(k1_g/2), n)
        k3_f = dh*f(x[i]+(dh/2), y[i]+(k2_f/2), z[i]+(k2_g/2))
        k3_g = dh*g(x[i]+(dh/2), y[i]+(k2_f/2), z[i]+(k2_g/2), n)
        k4_f = dh*f(x[i]+(dh), y[i]+(k3_f), z[i]+(k3_g))
        k4_g = dh*g(x[i]+(dh), y[i]+(k3_f), z[i]+(k3_g), n)
        y[i+1] = (y[i]+(1/6)*(k1_f+k2_f*2+k3_f*2+k4_f))
        z[i+1] = (z[i]+(1/6)*(k1_g+k2_g*2+k3_g*2+k4_g))
        i += 1
    return y

# Main Loop to generate plots
for j in range(len(n)):
    y = LE(n[j])
    plt.plot(x, y, 'o', markersize=0.5, label="numerical sol., n={}".format(n[j]))
    if j == 0:
        y1 = y
    if j == 1:
        y2 = y
    else:
        y3 = y

# LE Sols. Plots
plt.xlabel("epsilon")
plt.ylabel("theta")
plt.title("LE eqn. sols. at n={}, Tot_epsilon={}".format(n, tot_epsilon))
plt.legend()
plt.show()
#plt.savefig("LE_tot{}".format(tot_epsilon))
plt.clf()

# Error Plots between Numerical and Analytical sols.
plt.plot(x, np.abs(y1 - theta0(x)), 'o', markersize=0.5, label="error, n=0")
plt.plot(x, np.abs(y2 - theta1(x)), 'o', markersize=0.5, label="error, n=1")
plt.plot(x, np.abs(y3 - theta5(x)), 'o', markersize=0.5, label="error, n=5")
plt.xlabel("epsilon")
plt.ylabel("error")
plt.legend()
plt.title("Abs. Error")
plt.show()
#plt.savefig("LE_error_tot{}".format(tot_epsilon))
