# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 22:10:31 2018

@author: Kosar
"""
import numpy as np
import matplotlib.pyplot as plt

with open('statistics.txt') as infile:
    file = open("kinetic.txt", "w") 
    file2 = open("Temps.txt", "w")
    file3 = open("Times.txt","w")
    for line in infile:
#        coloumn2.append(line.split()[2])
        
        file.write(line.split()[3])
        file.write("\n")
        file2.write(line.split()[2])
        file2.write("\n")
        file3.write(line.split()[1])
        file3.write("\n")
file.close() 
#kinetic = np.loadtxt("kinetic.txt")
#potential = np.loadtxt("potential.txt")
#print(len(kinetic))
#print(len(potential))

#coloumn2 = []
#with open(r"50.txt", "r+") as f:
#    data = f.readlines()
##    print (data)
#    for line in data:
#        coloumn2.append(line.strip().split(" ")[2])
#
#E_total = kinetic + potential
#time = np.loadtxt("time.txt")
#timeSI = time *1.00224e-13
#plt.plot(timeSI,kinetic,label = r"$E_k$")
#plt.plot(timeSI,potential,label = r"$E_p$")
#plt.plot(timeSI,E_total,label = r"$E_total$")
#
#
#csfont = {'fontname':'times new roman','size'   : 16}
#hfont = {'fontname':'Helvetica'}
#
##plt.title('title',**csfont)
#plt.xlabel('Time', **csfont)
#plt.ylabel('Energy',**csfont)
#plt.legend()
#plt.show()