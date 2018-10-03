# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 19:34:11 2018

@author: haninm_adm
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.sparse.linalg import spsolve_triangular


def AnalyticalSolution(x):
    u=-np.exp(-10*x) - (1-np.exp(-10))*x + 1
    return u


def Create_Matrix_A(a,b,c):
    n = len(b)
    A = np.zeros((n,n))
    A[0,0]=b[0]
    A[0,1] = c[0]
    A[n-1,n-1] = b[n-1]
    A[n-1,n-2]=a[n-2]

    for i in range(1,n-1):
        A[i,i]=b[i]
        A[i,i+1]=c[i]
        A[i,i-1]=a[i-1]
    return A;

def CreateSparseMatrix(a,b,c):
    n = len(b)
    row = np.arange(n) #main diag
    row = np.append(row, np.arange(n-1)) #diag + 1
    row = np.append(row, np.arange(1,n)) # diag - 1
    
    col = np.arange(n)
    col = np.append(col, np.arange(1,n))
    col = np.append(col, np.arange(0,n-1))
    
    data = np.array(([]))
    data = np.append(data, b)
    data = np.append(data, c)
    data = np.append(data, a)
    
    return csr_matrix((data, (row, col)))

def CreateSimplifiedMatrixA(a, b, n):
    diag_vector = b * np.ones((n,1))
    non_diag_vector = a * np.ones((n-1, 1))
    return CreateSparseMatrix(non_diag_vector, diag_vector, non_diag_vector)

def Construct_Solution_Vector(h):
    
    f = np.zeros((n,1))
    for i in range(1,n+1):
        f[i-1,0] = + 100*(np.exp(-10*h * i))
    return f

def OdeSparseSolver(a,b,n):
    h = 1 / (n + 1)
    A = CreateSimplifiedMatrixA(a,b,n)
    b = h**2 * Construct_Solution_Vector(h)
    v = spsolve(A, b)
    return v

def plot_relativ_error(a,b,n):
    v=[]
    h_i = []
    v = OdeSparseSolver(a,b,n)
    for i in range(1,n+1):
        u[i - 1, 0] = AnalyticalSolution(i*h)
    for i in range(1,n+1):
        e[i - 1, 0] = np.log10(np.abs((v[i-1]-u[i-1])/u[i-1]))
    for item in [10,10**2,10**3,10**4,10**5,10**6,10**7]:
        h_i.append(np.log10(1/item))
        
    rel_e = [-1.18,-3.088,-5.08,-7.079,-9.004,-6.771,-5.96]
    plt.plot(h_i,rel_e)

if __name__ == '__main__':
    start_time = time.time()
    n = 10000000
    u = np.zeros((n,1))
    e = np.zeros((n,1))
    t = np.linspace(0, 1, n+2)
    h = 1 / (n + 1)
    for i in range(1,n+1):
        u[i - 1, 0] = AnalyticalSolution(i*h)

    a = -1
    b = 2
    v1 = OdeSparseSolver(a, b, n)
    for i in range(1,n+1):
        e[i - 1, 0] = np.log10(np.abs((v1[i-1]-u[i-1])/u[i-1]))  
    print("--- %s seconds ---" % (time.time() - start_time))
#    print(v1)
#    print(v2)
    