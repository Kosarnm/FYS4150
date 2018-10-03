# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 09:30:45 2018

@author: Kosar
"""
import numpy as np
from scipy import linalg
 # Find largest off-diag. element a[k,l]
from scipy.sparse import csr_matrix
 
def maxoff(A,n):
    maximum = 0
    for i in range(0,n):
        for j in range(i+1,n):
            if(np.abs(A[i][j]) > maximum):
                maximum = np.abs(A[i][j])
                l = j
                k = i
    return maximum,l,k

def rotate(A,R,k,l,n):
    
    if (A[k][l] != 0):
        tau = (A[l][l] - A[k][k])/(2*A[k][l])
        t = - tau - np.sqrt(1 + tau**2)
        c = 1/np.sqrt(1+t**2)
        s = c*t
    else:
        c = 0.0
        s = 1.0


    for i in range(0,n):
        if i!=k and i!=l:
            a_ik = A[i][k]
            a_il = A[i][l]
            A[i][k] = c*a_ik - s*a_il
            A[k][i] = A[i][k]
            A[i][l] = c*a_il + s*a_ik
            A[l][i] = A[i][l]
        r_ik = R[i][k]
        r_il = R[i][l]
        R[i][k] = c*r_ik - s*r_il
        R[i][l] = c*r_il +s*r_ik
    
    a_kk = A[k][k]
    a_ll = A[l][l]
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll
    A[k][l] = 0.0
    A[l][k] = 0.0


def jacobi_method(A, R, n):
    for i in range(0,n):
        for j in range(0,n):
            if i==j:
                R[i][j] = 1.0
            else:
                R[i][j] = 0.0
    tol = 1.0e-10
    max_iteration = n**4
    iteration = 0
    max_offdiagonal,l,k = maxoff(A,n)
    while abs(max_offdiagonal) > tol and iteration < max_iteration :
        max_offdiagonal,l,k = maxoff(A,n)
        rotate(A,R,k,l,n)
        iteration +=1
    print('number of iteration: ' ,iteration , '\n')

    
    
def test_rotate():
    A = np.array([[3,2],[2,5]], dtype=float)
    print(np.linalg.eigvals(A))
    R = np.array([[0.,-1.],[1.,0.]])
    k = 0
    l = 1
    n = 2
    rotate(A,R,k,l,n)
    print(A)
    print('----')
    rotate(A,R,k,l,n)
    print(A)
    print('----')
    rotate(A,R,k,l,n)
#    print(R)
    print(A)

def test_jacobi():
    A = np.array([[3,2],[2,5]], dtype=float)
    print(np.linalg.eigvals(A))
    R = np.array([[0.,0.],[0.,0.]])
    n = 2
    jacobi_method(A, R, n)
def creatToeplitz_Matrix(ru_min,ru_max,n):
    h = (ru_max - ru_min)/n
    A= np.zeros((n,n), dtype=float)
    a = 2.0/h**2
    b = -1.0/h**2
    A[0,0]=a
    A[0,1]=b
    A[n-1,n-1]=a
    A[n-1,n-2]=b
    for i in range(1,n-1):
        A[i,i]=a
        A[i,i+1] = b
        A[i,i-1]=b
    return A
def creat_Matrix_oneElectron(ru_min,ru_max,n):
    h = (ru_max - ru_min)/n
    ru = []
    V =[]
    for i in range(0,n):
        ru.append(ru_min + i*h)
        V.append(ru[i]**2)
    A = np.zeros((n,n))
    A[0,0]=2/(h**2) + V[0]
    A[0,1]=-1/(h**2)
    A[n-1,n-1]=2/(h**2) + V[n-1]
    A[n-1,n-2]=-1/(h**2) 
    for i in range(1,n-1):
        A[i,i]=2/(h**2) + V[i]
        A[i,i+1] = -1/(h**2)
        A[i,i-1]=-1/(h**2)
    return A


def creat_matrix_TwoElectron(ru_min,ru_max,n,omega):
    h = (ru_max - ru_min)/n
    ru =[]
    V =[]
    for i in range(0,n):
        ru.append(ru_min + i*h)
        V.append(omega + 1/ru[i]**2)
    A = np.zeros((n,n))
    A[0,0]=2/(h**2) + V[0]
    A[0,1]=-1/(h**2)
    A[n-1,n-1]=2/(h**2) + V[n-1]
    A[n-1,n-2]=-1/(h**2) 
    for i in range(1,n-1):
        A[i,i]=2/(h**2) + V[i]
        A[i,i+1] = -1/(h**2)
        A[i,i-1]=-1/(h**2)
    return A
def Creat_rotate_matrix(n):
    R = np.zeros((n,n), dtype=float)
    return R
    
if __name__ == '__main__':
    n =100 
    R = Creat_rotate_matrix(n)
    A = creatToeplitz_Matrix(0,1,100) 
    B = linalg.eigvals(A)
    jacobi_method(A, R, n)

    
 
            
            
            
        
    