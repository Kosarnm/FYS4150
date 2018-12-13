# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 21:15:39 2018

@author: Kosar
"""

import numpy as np
import matplotlib.pyplot as plt
kinetic =np.loadtxt("kinetic.txt")
potential = np.loadtxt("potential.txt")
Time = np.loadtxt("Times.txt")
Temps = np.loadtxt("Temps.txt")
temp50 = np.loadtxt("temp50.txt")
temp100 = np.loadtxt("temp100.txt")
temp150 = np.loadtxt("temp150.txt")
temp200 = np.loadtxt("temp200.txt")
temp250 = np.loadtxt("temp250.txt")
temp300 = np.loadtxt("temp300.txt")
temp350 = np.loadtxt("temp350.txt")
temp400 = np.loadtxt("temp400.txt")
temp500 = np.loadtxt("temp500.txt")
temp600= np.loadtxt("temp600.txt")
temp700 = np.loadtxt("temp700.txt")
temp800= np.loadtxt("temp800.txt")
temp900 = np.loadtxt("temp900.txt")
temp1000 = np.loadtxt("temp1000.txt")
temp1100 = np.loadtxt("temp1100.txt")
temp1200 = np.loadtxt("temp1200.txt")
temp1300= np.loadtxt("temp1300.txt")
temp1400 = np.loadtxt("temp1400.txt")

time = np.loadtxt("time.txt")

temp = np.linspace(50,1400,55)
#diffSI = diffusion* 0.998e-7 
temp50SI = temp50*119.735
temp100SI = temp100*119.735
temp150SI = temp150*119.735
temp200SI = temp200*119.735
temp250SI = temp250*119.735
temp300SI = temp300*119.735
temp350SI = temp350*119.735
temp400SI = temp400*119.735
temp500SI = temp500*119.735
temp600SI = temp600*119.735
temp700SI = temp700*119.735
temp800SI = temp800*119.735
temp900SI = temp900*119.735
temp1000SI = temp1000*119.735
temp1100SI = temp1100*119.735
temp1200SI = temp1200*119.735
temp1300SI = temp1300*119.735
temp1400SI = temp1400*119.735

timeSI = time *1.00224e-13
TempsSI = Temps *119.735
sumation = sum(kinetic)
mean = sumation/10000
T = 2*mean/(500*3)
ratio = T/TempsSI

#T = temp1200[0]/temp1200[9998]
R50SI = temp50[0]/temp50[9998]
R100SI = temp100[0]/temp100[9998]
R150SI = temp150[0]/temp150[9998]
R200SI = temp200[0]/temp200[9998]
R250SI = temp250[0]/temp250[9998]
R300SI = temp300[0]/temp300[9998]
R350SI = temp350[0]/temp350[9998]
R400SI = temp400[0]/temp400[9998]
R500SI = temp500[0]/temp500[9998]
R600SI = temp600[0]/temp600[9998]
R700SI = temp700[0]/temp700[9998]
R800SI = temp800[0]/temp800[9998]
R900SI = temp900[0]/temp900[9998]
R1000SI = temp1000[0]/temp1000[9998]
R1100SI = temp1100[0]/temp1100[9998]
R1200SI = temp1200[0]/temp1200[9998]
R1300SI = temp1300[0]/temp1300[9998]
R1400SI =temp1400[0]/temp1400[9998]
#ratio.append(R50SI)
#ratio.append(R100SI)
#ratio.append(R150SI)
#ratio.append(R200SI)
#ratio.append(R300SI)
#ratio.append(R350SI)
#ratio.append(R400SI)
#ratio.append(R500SI)
#ratio.append(R600SI)
#ratio.append(R700SI)
#ratio.append(R800SI)
#ratio.append(R900SI)
#ratio.append(R1000SI)
#ratio.append(R1100SI)
#ratio.append(R1200SI)
#ratio.append(R1300SI)
#ratio.append(R1400SI)

#
#print(len(temp1200))
#print(T)


#plt.plot(timeSI,temp50SI,label ="T = 50")
#plt.plot(timeSI,temp100SI,label ="T = 100")
#plt.plot(timeSI,temp150SI,label ="T = 150")
#plt.plot(timeSI,temp200SI,label ="T = 200")
#plt.plot(timeSI,temp250SI,label ="T = 250")
#plt.plot(timeSI,temp300SI,label ="T = 300")
#plt.plot(timeSI,temp350SI,label ="T = 350")
#plt.plot(timeSI,temp400SI,label ="T = 400")
#plt.plot(timeSI,temp500SI,label ="T = 500")
#plt.plot(timeSI,temp600SI,label ="T = 600")
#plt.plot(timeSI,temp700SI,label ="T = 700")
#plt.plot(timeSI,temp800SI,label ="T = 800")
#plt.plot(timeSI,temp900SI,label ="T = 900")
#plt.plot(timeSI,temp1000SI,label ="T = 1000")
#plt.plot(timeSI,temp1100SI,label ="T = 1100")
#plt.plot(timeSI,temp1200SI,label ="T = 1200")
#plt.plot(timeSI,temp1300SI,label ="T = 1300")
#plt.plot(timeSI,temp1400SI,label ="T = 1400")

#t = np.linspace(0,1,len(ratio))
#plt.plot(t,ratio)
E_total = kinetic + potential
#time = np.loadtxt("time.txt")
#timeSI = time *1.00224e-13
#plt.plot(timeSI,kinetic,label = r"$E_k$")
#plt.plot(timeSI,potential,label = r"$E_p$")
#plt.plot(timeSI,E_total,label = r"$E_total$")
plt.plot(timeSI,ratio)

csfont = {'fontname':'times new roman','size'   : 16}
hfont = {'fontname':'Helvetica'}

#plt.title('title',**csfont)
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

