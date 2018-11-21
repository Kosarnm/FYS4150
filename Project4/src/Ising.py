# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 20:21:37 2018

@author: Kosar
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 10:01:17 2018

@author: Kosar
"""
#from __future__ import division
import numpy as np
from numpy.random import rand
import random
from random import randrange
import matplotlib.pyplot as plt
from numba import jit

class Ising:
    def __init__(self, lattice_dimention, temp):
        self._energy = 0
        self._magnetization = 0
        self._proboblity = 0
        self._heat_c = 0
        self._mean_2 = 0
        self._mean_mag_2 = 0
        self._susceptibility = 0
        self._dim = lattice_dimention
        self._grid = np.zeros((self._dim, self._dim), int)
        self._w = np.zeros((17))
        self._n_spins = self._dim ** 2
        self._mean = 0
        self._mean_mag =0
        self._no_of_accepted_states = 0
        self._temperature = temp
        self._expected_value = 0
        self._mean_mag_abs = 0
        self.initialize(temp)
                                       
    def periodic(self, x,dimension, translation):
            return (x + dimension + translation) % dimension
    
    def initialize(self, temp):
        random.seed(random.seed(100))
        
        for de in range(-8,9,4):
            self._w[de + 8] = np.exp(-de/temp) 
        
#        random initialization of spin matrix
        for i in range(self._dim):
            for j in range(self._dim):
                if (rand() > 0.5):
                    self._grid[i][j] = 1 
                else:
                    self._grid[i][j] = -1
                self._magnetization  += (self._grid[i][j] ) 
            
#        ordered initialization of spin matrix -all spins up
#        for i in range(self._dim):
#            for j in range(self._dim):
#                if (rand() > 0.5):
#                    self._grid[i][j] = 1 
#                else:
#                    self._grid[i][j] = 1
#                self._magnetization  += (self._grid[i][j] )
                

        for i in range(self._dim):
            for j in range(self._dim):
                a = self.periodic(i, self._dim, -1)
                b = self.periodic(j, self._dim, -1)
                self._energy -= self._grid[i][j] * (self._grid[a][j] + self._grid[i][b])
#                if rand() > 0.5 :
#                    self._magnetization  += self._grid[i][j] 
#                    
#                
#        for i in range(self._dim):
#            for j in range(self._dim):
#                a = self.periodic(i, self._dim, -1)
#                b = self.periodic(j, self._dim, -1)
#                
##                self._grid[i] [j] = self._grid[a][b]
#                  
#                self._magnetization  += self._grid[i][j] 
#                
                

#                
    def calcLocalEnergy(self, spin_matrix, a, b):
        E = -spin_matrix[a][b] * (spin_matrix[a][self.periodic(b,self._dim,-1)] +  
                      spin_matrix[self.periodic(a,self._dim,-1)][b] +
                      spin_matrix[a][self.periodic(b,self._dim,1)] +
                      spin_matrix[self.periodic(a,self._dim,1)][b])
        return E
       
    def expectedEnergy(self):
        num = 0.0
        den = 0.0
        b = 1.0 / self._temperature
        for i in range(self._dim):
            for j in range(self._dim):
               ei = -self.calcLocalEnergy(self._grid, i, j)
               pi = np.exp(-b * ei)
               den += pi
               num += ei * pi
        return num / den
   
    def Metropolis(self):
        for i in range(self._n_spins):
            b = randrange(self._dim)
            a = randrange(self._dim)
            deltaE = -2 * self.calcLocalEnergy(self._grid, a, b)

            if (rand() <= self._w[deltaE + 8]):
                self._grid[a][b] *= -1
                self._magnetization += 2*(self._grid[a][b])
#                print(self._magnetization)
                self._energy += deltaE
                self._no_of_accepted_states += 1
    @jit 
    def iterate_metropolis(self,iteration):
#        start = int(percentage * iteration)
        for i in range(iteration):
            self.Metropolis()
            if (-2 < self._energy <-1.6):
               self._proboblity += 1 
            self._mean += self._energy
            self._mean_mag += abs(self._magnetization)
#            self._mean_mag += self._magnetization
            self._mean_2 += self._energy**2
            self._mean_mag_2 += (abs(self._magnetization))**2
#            if i > start:
#                self._mean += self._energy
        self._mean /= iteration*self._n_spins
        self._proboblity /= iteration*self._n_spins*100
        self._mean_mag /= iteration*self._n_spins
        self._mean_mag_abs  /= iteration*self._n_spins
        self._mean_2 /= iteration*self._n_spins
        self._mean_mag_2 /= iteration*self._n_spins
        self._heat_c = (self._mean_2 - self._n_spins*self._mean**2)/self._temperature**2
        self._susceptibility = (self._mean_mag_2 -self._n_spins*self._mean_mag**2 )/(self._temperature)
        self._expected_value = self.expectedEnergy()
#
    def get_mean_energy(self):
        return self._mean
    
    def get_mean_mag(self):
        return self._mean_mag
    
    def get_expected_energy(self):
        return self._expected_value
    def get_heat_c(self):
        return self._heat_c
    def get_susceptibility(self):
        return self._susceptibility
    def get_accepted_spin(self):
        return self._no_of_accepted_states
    def get_probolity(self):
        return self._proboblity

if __name__ == "__main__":
    number_dims =20
    temp = 2.4
    iterations = np.linspace(1,1000, 100 )
    for i in iterations:
        ising_model = Ising(number_dims, temp)
        ising_model.iterate_metropolis(int(i))
        m_e = ising_model.get_mean_energy()
##        p = ising_model.get_probolity()
        heat = ising_model.get_heat_c()
        Susceptibility = ising_model.get_susceptibility()
        magnet = ising_model.get_mean_mag()
#        steady.append(ising_model.get_mean_energy())
#        steady2.append(ising_model.get_accepted_spin())
#        steady3.append(ising_model.get_heat_c())
#        steady4.append(ising_model.get_mean_mag())
    
#    plt.plot(steady) 
#    bins = np.array([-2,-1.98,-1.96,-1.94,-1.92,-1.9,-1.88,-1.86,-1.84,-1.82,-1.8,-1.78,-1.76,-1.74,-1.72,-1.7]) 
#    bins = np.array([-2,-1.95,-1.9,-1.85,-1.8,-1.75,-1.7]) 
#    plt.hist(steady,bins,normed=True)
#    csfont = {'fontname':'Times New Roman'}
#    plt.rcParams.update({'font.size': 14})
#    plt.title('Expected Energy as function of MC cycles' , **csfont)
#    plt.xlabel('MC cycles' , **csfont)
#    plt.ylabel('<E>' ,**csfont)
    steady = []
    steady_m =[]
    steady_s =[]
    steady_h =[]  
    iteration = 100   
    nt =20
    T =np.linspace(2, 2.3, nt)
    for i in range(nt):
        ising_model = Ising(40, T[i])
        ising_model.iterate_metropolis(int(iteration))
        heat = ising_model.get_heat_c()
        Susceptibility = ising_model.get_susceptibility()
        magnet = ising_model.get_mean_mag()
        steady.append(m_e)
        steady_m.append(magnet)
        steady_s.append(Susceptibility)
        steady_h.append(heat)