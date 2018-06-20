#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 10:04:54 2018

@author: menorahlam
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# This exercise is a practice in using RK method to solve ODEs.
# Let's introduce a simple ODEs and the initial consition is y(0)=1 and an interval of x=[0, 5]

def dy (y, x):
    return y*(x**2)


# The analyetical solution is the following
def y (x):
    return np.exp((1/3)*x**3)

def k(h, dy):
    return h*dy

# this dy has x=x_n + 1/2h and y=y_n + 1/2k1. Remember to correct this.
#def k2(h, dy):
#    return h*dy

# I introduce the 2nd order RK method for this session for simplicity
def yn (y, k1, k2, k3, k4):
    return y + (1/6)*(k1 + 2*k2 + 2*k3 + k4)



dx = 0.001  # Define Resolution
x0 = 0  # Define x_0
xf = 5  # Define x_f
y0 = 1.  # Define y_0
h = dx  # Define step
y_RK = [y0]
x = np.arange(x0, xf+dx, dx)  # Create an uniform x array

yi = y0
i = 0
while i < len(x)-1:
    xi = x[i]
    k1 = h*dy(yi, xi)
    k2 = h*dy(yi + 0.5*k1, xi + h/2)
    k3 = h*dy(yi + 0.5*k2, xi + h/2)
    k4 = h*dy(yi + k3, xi + h)
    y_new = 0.
    y_new = yn(yi, k1, k2, k3, k4)
    y_RK.append(y_new)
    i += 1
    yi = np.copy(y_new)

plt.plot(x, y(x), label='Analyetical Sol.')
plt.plot(x, y_RK, label='RK (4-order) Sol.')
plt.xlim(0, 2.5)
plt.ylim(0, 175)
plt.legend()
plt.title('RK Method Analysis')
plt.plot(x, np.abs(y_RK - y(x)), label='Error')
plt.show()