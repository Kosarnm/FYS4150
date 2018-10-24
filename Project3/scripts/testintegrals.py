# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 17:30:15 2018

@author: Kosar
"""

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import time

AU = 149597871000      #1 AU in Meters
sun_mass = 2E30        #Sun mass in kg

G = 4*np.pi**2



dt = 0.01
N = int(1E3)

m_sun = 1
 
x = np.zeros(N)
y = np.zeros(N)

vx = np.zeros(N)
vy = np.zeros(N)
ax = np.zeros(N)
ay = np.zeros(N)

x[0] = 1
y[0] = 0

vx[0] = 0
vy[0] = 2.0*np.pi

t1 = time.time()
def Grav_earth(x,y):
	r = np.sqrt(x**2+y**2)
	F_x = -G*x/r**3
	F_y = -G*y/r**3
	return F_x , F_y,r

#F_x0,F_y0 ,r = Grav_earth(x[0],y[0]) 
#ax[0] =  F_x0
#ay[0] =  F_y0

#Euler-cromer
for i in range(N-1):
	F_x, F_y, r = Grav_earth(x[i],y[i]) 
	vx[i+1] += vx[i] + F_x*dt  			#F_x = ax
	vy[i+1] += vy[i] + F_y*dt			#F_y =ay
	x[i+1] += x[i] + vx[i+1]*dt
	y[i+1] += y[i] + vy[i+1]*dt

#VelocityVerlet

#for i in range(N-1):
#    x[i+1] = x[i] + vx[i]*dt + 0.5*ax[i]*dt**2
#    y[i+1] = y[i] + vy[i]*dt + 0.5*ay[i]*dt**2
#    F_x, F_y,r= Grav_earth(x[i+1], y[i+1])
#    ax[i+1] = F_x
#    ay[i+1] =  F_y
#    vx[i+1] = vx[i] + 0.5*dt*ax[i] +  0.5*dt*ax[i+1]
#    vy[i+1] = vy[i] + 0.5*dt*ay[i] +  0.5*dt*ay[i+1]
#    x[i+1] = x[i] + vx[i]*dt 
#    y[i+1] = y[i] + vy[i]*dt 
#    F_x, F_y,r= Grav_earth(x[i+1], y[i+1])
#    ax[i+1] = F_x
#    ay[i+1] =  F_y
#    vx[i+1] = vx[i] + 0.5*dt*ax[i+1]
#    vy[i+1] = vy[i] + 0.5*dt*ay[i+1]

if __name__ == "__main__"
#csfont = {'fontname':'Times New Roman'}
##plt.plot(x,y)
#plt.rcParams.update({'font.size': 22})
#plt.xlabel('x[AU]',**csfont)
#plt.ylabel('y[AU]',**csfont)
#plt.title('Earth-Sun System Euler-Cromer[dt=0.0001]',**csfont)
t2 = time.time()
print(t2-t1)
#plt.hold
#plt.plot(x,y)
#plt.show()