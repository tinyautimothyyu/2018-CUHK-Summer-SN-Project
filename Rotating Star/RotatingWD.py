import numpy as np
import matplotlib.pyplot as plt

piue = 2  # Chemical Coefficient 
a = 6.00e22  # dynes cm^-2
b = piue * 9.82e5  # gcm^-3
G = 6.67e-8  # Gravitational Constant cgs
Ye = 0.5  # Z/A
K = 1.0036e13 * Ye**(5/3)  # Non-Relativistic EOS Constant
K_R = 1.2435e15 * Ye**(4/3)  # Relativistic EOS Constant
n = 3/2  # Non-relativistic constant for EOS
n_R = 3  # Relativistic Constant for EOS

# Initiate Boundary Conditions
Hf_Polytropic = 0
Hf_WD = (8*a)/b
rho_max = 

# Factorial
def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

# Binomial Coefficient 
def BC (n, k):
    return factorial(n) / (factorial(k) * factorial(n-k))


# Function to compute Legendre Polynomials
def LP (x, n):
    P = 0
    if n == 0:
        P = 1
    elif n == 1:
        P = x
    else:
        k = 0
        while k <= n:
            P_k = BC(n, k)*BC((n+k), k)*((x-1)/2)**k
            #print (P)
            P += P_k
            k += 1
    return P

print (LP(2, 2))

# Polytropic Pressure equation
def EOS_Polytropic (K, rho, N):
    return K*rho**(1+1/N)

# WD Pressure Equation
def EOS_WD (x):  # x is a function of rho. a is constant
    return a*(x*(2*x**2 - 3)*np.sqrt(x**2 + 1) + 3*np.arcsinh(x))

# this x is used in EOS_WD
def x (rho):  # b is constant
    return (rho/b)**(1/3)

# Polytropic Enthalpy
def H_Polytropic (N, P, rho):  # P is pressure
    return ((1+N)*P)/rho


# WD Enthalpy
def H_WD (rho):  # a and b are the cosntant used in EOS_WD
    return (8*a/b)*np.sqrt(1+(rho/b)**(2/3))

# rho as a function of H_Polytropic (enthalpy)
def rho_Polytropic (H, K, N):
    return (H/(K*(1+N)))**N

# rho as a function of H_WD
def rho_WD (H):
    return b*((b*H/(8*a))**2 - 1)**(3/2)


# Polytropic equator radius
def Re_Polytropic (N, rho_max, Hhat_max):
    return np.sqrt((1+N)*K*G**(-1) * rho_max**(-1 + 1/N) * (1/Hhat_max))


# WD equator radius
def Re_WD (x_max, rho_max, Hhat_max):
    return np.sqrt(8*a*np.sqrt(1 + x_max**2)*(1/(b*G*rho_max*Hhat_max)))


# Dimensionless Treatment
def rho (rhohat):
    return rhohat*rho_max

def phi (phihat, Re):
    return phihat*G*rho_max*Re**2


def P (Phat, Re):
    return Phat*G*(Re*rho_max)**2

def H (Hhat):
    return Hhat**G*rho_max*Re**2

def omega (omegahat, Re):
    return omegahat*Re

def shaihat (omegahat):
    return -0.5*omegahat**2


def h0squarehat_Poly (A, B):
    return -(phihat(A)-phihat(B))/(shaihat(A)-shaihat(B))

def Chat_Poly (A, B):
    return phihat(A) + h0squarehat_Poly(A, B)*shaihat(A)


#Initiate Density Profile
KDIV = 129
rho_r = np.zeros(KDIV).reshape((KDIV,1))
rho_piu = np.zeros(KDIV).reshape((KDIV,1))
rho_array = np.column_stack((rho_piu, rho_r)) + 0.1