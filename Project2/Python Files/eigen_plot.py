# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 21:06:48 2018

@author: Kosar

"""
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
#def plot_harmonic(l):
if __name__ == '__main__':
    ru_min = 0.
    ru_max = 10
    n=200
    ru= []
    V =[]
    U =[]
    ui_2=[]
    p_i=[]
    h = (ru_max - ru_min)/n
#    h = 1
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
        
    lamd, Ul = LA.eig(A)
    ui = Ul[:,35]
    for i in range(0,200):
        ui_2.append(ui[i]**2)
    for i in range(0,200):
        p_i.append(ru_min + h*i)
    plt.figure()
    plt.plot(p_i,ui_2)
#    for j in range(0,n):
#        for i in range(0,n):
#            U[]
    
#plot_harmonic(3)
 