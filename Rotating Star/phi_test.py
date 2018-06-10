import numpy as np
import matplotlib.pyplot as plt
import time

KDIV = 129
NDIV = 129
LMAX = 16
r_max = 16/15


# Initiate Density Profile
rho = np.zeros(KDIV**2).reshape((KDIV,KDIV)) + 0.1



j = np.linspace(1, NDIV, NDIV)
i = np.linspace(1, KDIV, KDIV)

# Initiate r and piu arrays
r = r_max*(j - 1) / (NDIV - 1)
piu = (i - 1) / (KDIV - 1)

r[0] = 0.001
piu[0] = 0.001

def P2n_array ():
    P2n = np.zeros((LMAX+1)*(KDIV)).reshape((LMAX+1,KDIV))
    n = 0
    while n <= LMAX:
        i = 0
        while i <= KDIV-1:
            P2n[n][i] = LP(piu[i], 2*n)
            i += 1
        n += 1
    return P2n

# Initiate P2n_array
P2n = P2n_array()


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


def f2n (r1, r2, n):
    if r1 < r2:
        return (r2**(2*n+2)) / (r1**(2*n+1))
    elif r1 > r2:
        return (r1**(2*n)) / (r2**(2*n-1))
    else:
        return r1


def D1 (n, k):
    i = 0
    D = 0
    while i <= KDIV-2:
#        Di = (1/6)*(piu[i+2]-piu[i])*(LP(piu[i], 2*n)*rho[i][k]+4*LP(piu[i+1], 2*n)*rho[i+1][k]+LP(piu[i+2], 2*n)*rho[i+2][k])
        Di = (1/6)*(piu[i+2]-piu[i])*(P2n[n][i]*rho[i][k]+4*P2n[n][i+1]*rho[i+1][k]+P2n[n][i+2]*rho[i+2][k])
        D += Di
        i += 2
    return D

def D2 (n, j):
    k = 0
    D = 0
    while k <= NDIV-2:
        Di = (1/6)*(r[k+2]-r[k])*(f2n(r[k], r[j], n)*D1(n, k)+4*f2n(r[k+1], r[j], n)*D1(n, k+1)+f2n(r[k+2], r[j], n)*D1(n, k+2))
        D += Di
        k += 2
    return D

def PHI (i, j):
    start_time = time.time()
    n = 0
    PHI = 0
    while n <= LMAX:
        print (n)
#        PHIi = D2(n, j)*LP(piu[i], 2*n)
        PHIi = D2(n, j)*P2n[n][i]
#        print (LP(piu[i], 2*n))
        PHI += PHIi
        n += 1
    print("--- {}s seconds ---".format(time.time() - start_time))
    return -4*np.pi*PHI

        
print (P2n_array())



print (PHI(1,1))