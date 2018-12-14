
"""
Solve linear system using LU decomposition and Gaussian elimination
"""

import numpy as np
from scipy.linalg import lu, inv
import matplotlib.pyplot as plt
import timeit

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

def Coeff_Matrix(a,b,n):
    A = np.zeros((n,n))
    A[0,0]=b
    A[0,1] = a
    A[n-1,n-1] = b
    A[n-1,n-2]=a

    for i in range(1,n-1):
        A[i,i]=b
        A[i,i+1]=a
        A[i,i-1]=a
    return A

def CreateSimplifiedMatrixA(a, b, n):
    diag_vector = b * np.ones((n,1))
    non_diag_vector = -a * np.ones((n, 1))
    return Create_Matrix_A(non_diag_vector, diag_vector, non_diag_vector)

def Construct_Solution_Vector(h):
    
    f = np.zeros((n,1))
    for i in range(1,n+1):
        f[i-1,0] = + 100*(np.exp(-10*h * i))
    return f
    
#a = []
#print('Enter a vector for a' )
#for x in input().split():
#    a.append(int(x))
#print('Enter a vector for b' )
#
#b = []
#for x in input().split():
#    b.append(int(x))
#    
#print('Enter a vector for c' )
#c = []
#for x in input().split():
#    c.append(int(x))   
#    
#n= int(input('enter an integer for n:'))
#
#A = np.zeros((n,n))
#A[0,0]=b[0]
#A[0,1] = c[0]
#A[n-1,n-1] = b[n-1]
#A[n-1,n-2]=a[n-2]
#
#for i in range(1,n-1):
#    A[i,i]=b[i]
#    A[i,i+1]=c[i]
#    A[i,i-1]=a[i-1]
    
#
#h = 1.0 / (n + 1)
#f = np.zeros((n,1))
#for i in range(1,n+1):
#    f[i-1,0] = + 100*(np.exp(-10*h * i))

def gausselim(A,B):
    """
    Solve Ax = B using Gaussian elimination and LU decomposition.
    A = LU   decompose A into lower and upper triangular matrices
    LUx = B  substitute into original equation for A
    Let y = Ux and solve:
    Ly = B --> y = (L^-1)B  solve for y using "forward" substitution
    Ux = y --> x = (U^-1)y  solve for x using "backward" substitution
    :param A: coefficients in Ax = B
    :type A: numpy.ndarray of size (m, n)
    :param B: dependent variable in Ax = B
    :type B: numpy.ndarray of size (m, 1)
    """
    # LU decomposition with pivot
    pl, u = lu(A, permute_l=True)
    # forward substitution to solve for Ly = B
    y = np.zeros(B.size)
    for m, b in enumerate(B.flatten()):
        y[m] = b
        # skip for loop if m == 0
        if m:
            for n in range(m):
                y[m] -= y[n] * pl[m,n]
        y[m] /= pl[m, m]

    # backward substitution to solve for y = Ux
    x = np.zeros(B.size)
    lastidx = B.size - 1  # last index
    for midx in range(B.size):
        m = B.size - 1 - midx  # backwards index
        x[m] = y[m]
        if midx:
            for nidx in range(midx):
                n = B.size - 1  - nidx
                x[m] -= x[n] * u[m,n]
        x[m] /= u[m, m]
    return x
def OdeSolve(a,b,n):
    h = 1 / (n + 1)
    A = CreateSimplifiedMatrixA(a,b,n)
    b = h**2*Construct_Solution_Vector(h)
    v = gausselim(A,b)
    return v

def OdeSolveSimplified(a,b,n):
    h = 1 / (n + 1)
    A = Coeff_Matrix(a,b,n)
    b = h**2*Construct_Solution_Vector(h)
    v = gausselim(A,b)
    return v
    
if __name__ == '__main__':
##    a = np.loadtxt('a-vector.txt')
##    b = np.loadtxt('b-vector.txt')
##    c = np.loadtxt ('c-vector.txt')
##    n = len(b)
#    e= np.zeros((n,1))
#    u = np.zeros((n,1))
##    h1 = np.zeros((n,1))
##    for i in range(0,n):
##        h1[i] = 1/(i+1)
#    a=-1*np.ones(n-1)
#    b=2*np.ones(n)
#    c=-1*np.ones(n-1)
#    A = Create_Matrix_A(a,b,c)
#    b = h**2*Construct_Solution_Vector()
#    v = gausselim(A,b)
#    
##    for i in range(0,n):
##        u[i] = -np.exp(-10*i) - (1-np.exp(-10))*i + 1
##        e[i]= np.log10(np.abs((v[i]-u[i])/u[i]))
##    print(e)
##    print (v)
#    t=np.linspace(0,1,n)
#    plt.plot(t, v ,'r', label = 'n=100')
#    plt.title('Numerical solution for different n')
##    plt.plot(np.log10(h1),e,'b-')
#    plt.legend()
#    plt.show()
    n = 100
    diag = 2
    non_diag = -1
    OdeSolve(non_diag,diag,n)
#    OdeSolveSimplified(non_diag,diag,n)
#    print(timeit.timeit('OdeSolve(-1,2,10)', globals=globals()))
    