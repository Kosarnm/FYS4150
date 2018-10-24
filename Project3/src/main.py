# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 16:18:56 2018

@author: Kosar
"""
#from vec3 import TripleVector
from planet2 import Planet
from SolarSystem1 import SolarSystem
#from solv import Solver
import numpy as np
import matplotlib.pyplot as plt
def main():
    p1 = Planet(0,0,0.0,0,1,'sun')
    p2 = Planet(1,0,0,2*np.pi, 3e-6 ,'earth')
#    p3 = Planet(5.2,0,0,2*np.pi/np.sqrt(5.2),1,'jupiter')
#    p4 = Planet(1.52, 0,0, 2*np.pi/np.sqrt(1.52), 3.3e-7,'mars')
#    p5 = Planet(0.72, 0,0, 2*np.pi/np.sqrt(0.72), 2.4e-6,'Venus')
#    p6 = Planet (9.54, 0,0, 2*np.pi/np.sqrt(9.54), 2.75e-4,'Saturn')
#    p7 = Planet(0.3075 , 0,0, 12.44, 1.65e-7,'Mercury')
#    p8 = Planet(19.19, 0,0, 2*np.pi/np.sqrt(19.19) , 4.4e-5,'Uranus')
#    p9 = Planet(30.06, 0,0, 2*np.pi/np.sqrt(30.06), 5.15e-5,'Neptun')
#    p10 = Planet(39.53, 0,0, 2*np.pi/np.sqrt(39.53), 6.55e-9,'Pluto')
    S = SolarSystem()
    S.add_planet(p2)
#    S.add_planet(p3)
#    S.add_planet(p3)
#    S.add_planet(p4)
#    S.add_planet(p5)
#    S.add_planet(p6)
#    S.add_planet(p7)
#    S.add_planet(p8)
#    S.add_planet(p9)
#    S.add_planet(p10)
#    S.add_planet(p1)
    x_earth,y_earth = S.euler(p2,0.0001,1E5)
    x_sun,y_sun = S.euler(p1,0.0001,1E5)
#    x_jupiter,y_jupiter = S.Velocity_verlet(p3,0.001,1E4)
#    x_mars,y_mars = S.Velocity_verlet(p4,0.01,1E4)
#    x_Venus,y_Venus = S.Velocity_verlet(p5,0.01,1E4)
#    x_Saturn,y_Saturn = S.Velocity_verlet(p6,0.01,1E4)
#    x_Mercury,y_Mercury = S.Velocity_verlet(p7,0.01,1E5)
#    x_Uranus,y_Uranus = S.Velocity_verlet(p8,0.01,1E5)
#    x_Neptun,y_Neptun = S.Velocity_verlet(p9,0.01,1E5)
#    x_Pluto,y_Pluto = S.Velocity_verlet(p10,0.01,1E5)
    plt.plot(x_earth,y_earth,label = 'Earth')
#    plt.plot(x_sun,y_sun,label = 'Sun')
#    plt.plot(x_jupiter,y_jupiter,label = 'Jupiter')
#    plt.plot(x_mars,y_mars,label ="Mars")
#    plt.plot(x_Venus,y_Venus,label ="Venus")
#    plt.plot(x_Saturn,y_Saturn,label='Saturn')
#    plt.plot(x_Mercury,y_Mercury,label='Mercury')
#    plt.plot(x_Uranus,y_Uranus,label='Uranus')
#    plt.plot(x_Neptun,y_Neptun ,label='Neptun')
#    plt.plot(x_Pluto,y_Pluto ,label='Pluto')
   
#    plt.legend()
#    plt.show()
#    plt.plot(x_earth,y_earth)
    csfont = {'fontname':'Times New Roman'}
    plt.rcParams.update({'font.size': 18})
    plt.xlabel('x[AU]',**csfont)
    plt.ylabel('y[AU]',**csfont)
#    plt.title('Earth-Sun System,Euler-Cromer[dt=0.0001]',**csfont)
#    plt.plot(x_earth,y_earth,'g')
#    y = y_earth
#    plt.title(r'Three body System with Jupiter_mass = 1e- AU',**csfont)
#    print(y)
#    plt.plot(x_earth1,y_earth1,label='V = $2.4\pi$')
#    plt.plot(x_earth2,y_earth2,label='V = $2.6\pi$')
#    plt.plot(x_earth3,y_earth3,label='V = $2.8\pi$')
#    plt.plot(x_earth4,y_earth4,label='V = $2.83\pi$')
#    plt.legend()
    
#    plt.plot(x_sun,y_sun)
#    plt.plot(x_jupiter , y_jupiter)
#    plt.plot(x_mars,y_mars)
#    plt.show()
    
#    --total_energy
    
#    total = S.total_energy(p1,p2)
#    --angular momentum
#    angular_momentum = S.angular_momentum(p2,0.01,1E6)
#    --prehelion precsition
#    p =[]
#    xp_mercury,yp_mercury = S.prehelion_with_VV(p7,0.001,1E5)
#    for i in range(0,int(1E5)):
#         p.append(np.sqrt(xp_mercury[i]**2 + yp_mercury[i]**2))
#    for i in range(0,len(p)):
#         tan.append(6*np.pi * 4*np.pi**2/p[i]/365*174)
##     tan = min_value_y/min_value_x
#     time = np.linspace(0,100,len(p))
#     plt.plot(time,tan) 
    
main()
    
    
    