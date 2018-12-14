# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 10:09:14 2018

@author: Kosar
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import scipy
import time
from scipy import linalg
from sympy.solvers import solve
def AnalyticalSolution(x):
    u=-np.exp(-10*x) - (1-np.exp(-10))*x + 1
    return u
n=10


A= np.zeros((n,n))
A[0,0]=2
A[0,1]=-1
A[n-1,n-1]=2
A[n-1,n-2]=-1
for i in range(1,n-1):
    A[i,i]=2
    A[i,i+1] = -1
    A[i,i-1]=-1
    
t_start = time.time()  
h = 1.0 / (n + 1)

f = np.zeros((n,1))
for i in range(1,n+1):
    f[i-1,0] = + 100*(np.exp(-10*h * i))

#t = np.linspace(0,1,n)
#Ainv = inv(A)
#v= h**2 * (np.dot(Ainv,f))
#h_i =[]

#u = np.zeros((n,1))
#e = np.zeros((n,1))
#t = np.linspace(0, 1, n)
start_time_1 = time.time()
LU = linalg.lu_factor(A)
x_LU = linalg.lu_solve(LU,f)
time_LU = time.time() - start_time_1
#for i in range(1,n+1):
#   u[i - 1, 0] = AnalyticalSolution(i*h)
#    
#for i in range(1,n+1):
#    e[i - 1, 0] = np.log10(np.abs((v[i-1]-u[i-1])/u[i-1]))
#for item in [10,10**2,10**3,10**4,10**5,10**6,10**7]:
#    h_i.append(np.log10(1/item))
#rel_e = [-1.18,-3.088,-5.08,-7.079,-9.004,-6.771,-5.96]

#plt.plot(t,u)

#plt.plot(t[1:-1],u,'r--')
#plt.plot(t[:,np.newaxis],v ,label = 'Numerical')
#plt.plot(t, AnalyticalSolution(t),'g--',label='Analytical')
#plt.grid(linestyle='--', linewidth=2)
#plt.xlabel('n = 100')
#plt.plot(h_i,rel_e)
#plt.xlabel('log(h)')
#plt.ylabel('log(e)')
##plt.legend()
#plt.show()

    


