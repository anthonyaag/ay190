#!/usr/bin/env python
import sys,math
import numpy as np
import matplotlib.pyplot as pl
import scipy as sp
from plot_defaults import *

# Global variables for the grid
x = np.arange(0,100,0.1)

dx = x[1] - x[0]
v = 0.1
n = len(x)
y = np.zeros(n)
cfl = 1.0
dt = .5
t = 0.0

# Initial Data
sigma = np.sqrt(1.0)
x0 = 30.0

def analytic(x,x0,sigma,t,v):
    mu = t*v+x0
    dem = 2*sigma **2
    return np.array(map(lambda i: math.exp( -1 * (i - mu) **2)/dem,x))

# Get initial analytic y data
y = analytic(x,x0,sigma,t,v)

def updateUpwind(y,dx,dt,v,n):
    temp = np.zeros(n)
    temp[1:n] = dt*v*(y[1:n] - y[0:n-1])/dx
    return y - temp

def updateFTCS(y,dx,dt,v,n):
    temp = np.zeros(n)
    temp[1:n-1] = dt*v*(y[2:n] - y[0:n-2])/(2*dx)
    return y - temp

def updateLFM(y,dx,dt,v,n):
    y[1:n-1] = (y[2:n] + y[0:n-2])/2 - dt*v*(y[2:n] - y[0:n-2])/(2*dx)
    return y
    
    

# evolve (and show evolution)
pl.ion()
pl.figure()
pl.plot(x,y,'x-') # numerical data
pl.plot(x,analytic(x,x0,sigma,t,v),'r-') # analytic data
pl.show()

error = []
time = []
# Implement a moving gaussian
for i in range(1200):
    t += dt
    
    y = updateFTCS(y,dx,dt,v,n)
    # Boundary Conditions
    y[0] = y[1] 
    y[-2] = y[-1]
    yana = analytic(x,x0,sigma,t,v)
    try:
        err = sum(abs(y-yana)) * dx
    except:
        break

    pl.clf()
    pl.plot(x,y,'b-')
    pl.plot(x,yana,'r-')
    pl.draw()
    
    # Save the time and errors
    time.append(t)
    error.append(err)


pl.show()

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.16)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.95)
myfig.subplots_adjust(right=0.96)

p1, = pl.semilogy(time,error,"b",linewidth=1.5)

# prepare x and y ranges

xmin = time[0]
xmax = time[-1]
ymin = 0
ymax = 1


# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()
'''
# set locators of tick marks
xminorLocator = pl.MultipleLocator(5)
xmajorLocator = pl.MultipleLocator(25)
yminorLocator = pl.MultipleLocator(0.05)
ymajorLocator = pl.MultipleLocator(.25)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
'''
# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("X",labelpad=2)
pl.ylabel("Y",labelpad=8)

# legend
pl.legend( [p1], ["Error"], loc=(0.4,0.6), frameon=False, labelspacing = 0.25 )

pl.savefig("4.pdf")
