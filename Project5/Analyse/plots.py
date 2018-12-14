# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 21:15:39 2018

@author: Kosar
"""

import numpy as np
import matplotlib.pyplot as plt
kinetic =np.loadtxt("kinetic.txt")
potential = np.loadtxt("potential.txt")
Temps = np.loadtxt("Temps.txt")

time = np.loadtxt("time.txt")

temp = np.linspace(50,1400,55)
#diffSI = diffusion* 0.998e-7 

timeSI = time *1.00224e-13
TempsSI = Temps *119.735
sumation = sum(kinetic)
mean = sumation/10000
T = 2*mean/(500*3)
ratio = T/TempsSI


E_total = kinetic + potential
#time = np.loadtxt("time.txt")
#timeSI = time *1.00224e-13
#plt.plot(timeSI,kinetic,label = r"$E_k$")
#plt.plot(timeSI,potential,label = r"$E_p$")
#plt.plot(timeSI,E_total,label = r"$E_total$")
plt.plot(timeSI,ratio)

csfont = {'fontname':'times new roman','size'   : 16}
hfont = {'fontname':'Helvetica'}

#plt.xlabel('Time', **csfont)
#plt.ylabel('Energy',**csfont)
#plt.legend()
#plt.show()
#plt.plot(timeSI,kinetic)
#csfont = {'fontname':'times new roman','size'   : 16}
#hfont = {'fontname':'Helvetica'}

#plt.title('title',**csfont)
#plt.xlabel('Time', **csfont)
#plt.ylabel(r'T/T_initial',**csfont)
##plt.legend()
#plt.show()

