#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as pl
from plot_defaults import *

# set up grid
xmin = 0.0
xmax = 1.0

# boundary values
A = 0
B = 0.1

def calc_rhs(u,xx):
    # rhs routine
    # rhs[0] is rhs for y
    # rhs[1] is rhs for u
    rhs = np.zeros(2)
    rhs[0] = u
    rhs[1] = 12 * xx - 4

    return rhs

def FE(z,x,npoints,dx):
    # forward-Euler integrator
    
    # make an array for all points
    # entry 0 contains y
    # entry 1 contains y'
    yy = np.zeros((npoints,2))

    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0

    for i in range(npoints-1):
        yy[i+1,:] = yy[i,:] + dx*calc_rhs(yy[i,1],x[i])

    return yy

def RK2(z,x,npoints,dx):
    # Runge-Kutta2 Integrator
    # make an array for all points
    # entry 0 contains y
    # entry 1 contains y'
    yy = np.zeros((npoints,2))

    yy[0,0] = A # boundary value A for y at x=0
    yy[0,1] = z # guessed boundary value for y' at x=0

    # RK2 integrator
    for i in range(npoints -1):
        k1 = dx * calc_rhs(yy[i,1],x[i])
        k2 = dx * calc_rhs(yy[i,1] + 0.5 * k1[0], x[i] + 0.5 * dx)
        yy[i+1,:] = yy[i,:] + k2

    return yy

def shooting(npoints,integrator = FE,z0 = -1100000.0,z1 = 10000000.0):

    # set up grid
    x = np.linspace(xmin,xmax,npoints)
    # dx based on x[1] and x[0] 
    dx = x[1] - x[0]

    # get initial guess for derivative
    yy0 = integrator(z0,x,npoints,dx)
    yy1 = integrator(z1,x,npoints,dx)
    phi0 = yy0[npoints-1,0] - B
    phi1 = yy1[npoints-1,0] - B
    dphidz = (phi1 - phi0)/(z1 - z0)

    i = 0
    itmax = 100
    err = float("inf")
    criterion = 1.0e-12

    z0 = z1
    phi0 = phi1
    while (err > 1.0e-12 and i < itmax):
        z1 = z0 - phi0 / dphidz # secant update
        yy = integrator(z1,x,npoints,dx)
        phi1 = yy[npoints-1,0] - B
        dphidz = (phi1 - phi0)/(z1 - z0)
        err = np.abs(phi1) # your error measure
        z0 = z1
        phi0 = phi1
        i += 1
        print i,z1,phi1
    return (x,yy)

# Generate Data
(x,yy_FE) = shooting(100)
(x,yy_RK2) = shooting(100,RK2)
print yy_FE
# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.16)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.95)
myfig.subplots_adjust(right=0.96)

p1, = pl.plot(x,yy_FE[:,0],"b",linewidth=1.5)
p2, = pl.plot(x,yy_RK2[:,0],"r",linewidth=1.5)
p3, = pl.plot(x,2.0*x**3 - 2*x**2 + 0.1*x,"k.",linewidth=1.5)

# prepare x and y ranges

xmin = 0
xmax = 1.
ymin = -.25
ymax = .25


# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set locators of tick marks
xminorLocator = pl.MultipleLocator(.05)
xmajorLocator = pl.MultipleLocator(.25)
yminorLocator = pl.MultipleLocator(0.05)
ymajorLocator = pl.MultipleLocator(.25)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("X",labelpad=2)
pl.ylabel("Y",labelpad=8)

# legend
pl.legend( [p1,p2,p3], ["Forward Euler","RK2","Data"], loc=(0.4,0.6), frameon=False, labelspacing = 0.25 )

pl.savefig("1.pdf")




