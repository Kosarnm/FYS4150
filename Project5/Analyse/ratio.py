# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 20:44:22 2018

@author: Kosar
"""
import numpy as np
import matplotlib.pyplot as plt

first = np.loadtxt("first.txt")
sec = np.loadtxt("sec.txt")
third = np.loadtxt("third.txt")
fourth = np.loadtxt("forth.txt")
fifth = np.loadtxt("fifth.txt")
six = np.loadtxt("six.txt")

time = np.loadtxt("time.txt")
timeSI = time *1.00224e-13

csfont = {'fontname':'times new roman','size'   : 16}
hfont = {'fontname':'Helvetica'}

plt.plot(timeSI,first)
plt.plot(timeSI,sec)
plt.plot(timeSI,third)
plt.plot(timeSI,fourth)
plt.plot(timeSI,fifth)
plt.plot(timeSI,six)

#plt.title('title',**csfont)
plt.xlabel('Time', **csfont)
plt.ylabel(r'$\dfrac{T}{T_i}$',**csfont)
#plt.legend()
plt.show()


