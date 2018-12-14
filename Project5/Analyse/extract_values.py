# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 22:10:31 2018

@author: Kosar
"""

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
