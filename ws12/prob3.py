#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, sqrt,pi
import matplotlib.pyplot as pl
from plot_defaults import *

data = np.transpose(np.loadtxt("presupernova.dat"))

count = data[0]
mass = data[1]
radius = data[2]
temperature = data[3]
density = data[4]
rvelocity = data[5]
electronFraction = data[6]
radialVelocity = data[7]
ggrav = 6.67e-8

def lin_interp(x,y,ran):
    ''' Takes input x and f(x) = y in lists, assumed to be sorted such that x
    is a sorted list in increasing order. Takes input, ran, a list of 
    points also assumed to be sorted in increasing order.'''
    ret = []
    x_iter = 1
    ran_iter = 0
    while ran_iter < len(ran): # while not done with all of the given points
        # Check to see if you're in the correct bin
        if ran[ran_iter] < x[x_iter] or x_iter == len(x) - 1:
            ret.append( y[x_iter - 1] + (ran[ran_iter] - x[x_iter - 1]) \
            *(y[x_iter] - y[x_iter -1]) / (x[x_iter] - x[x_iter -1]))
            ran_iter += 1

        # Otherwise go to the next bin and try this point again
        else:
            x_iter += 1
    return ret

## Use linear interpolation to map onto a evenly space grid.
r = np.arange(1e7,1e9,1e6)
rho = lin_interp(radius,density,r)

def calc_rhs(u,rho,r):
    # rhs routine
    # rhs[0] is rhs for phi'
    # rhs[1] is rhs for u'
    rhs = np.zeros(2)
    rhs[0] = u
    rhs[1] = 4*pi*ggrav*rho - 2*u/r

    return rhs

def FE(rho,r):
    # forward-Euler integrator
    dr = r[1] - r[0]
    # make an array for all points
    # entry 0 contains phi
    # entry 1 contains phi'
    npoints = len(r)
    ret = np.zeros((npoints,2))

    ret[0,0] = 0 # boundary value A for phi at r=0
    ret[0,1] = 0  # guessed boundary value for phi' at r=0

    for i in range(npoints-1):
        ret[i+1,:] = ret[i,:] + dr*calc_rhs(ret[i,1],rho[i],r[i])

    return ret

################################
## Test Numerical Convergence ##
################################
"""
#Run a low resolution
r = np.linspace(1e7,1e9,1001)
rho = array([1]*len(r)) # Assume constant density 1
phi = FE(rho,r)[:,0] # Calculate phi
phi -= (phi[-1] + ggrav * 4*pi*r[-1]**2/3) # Set outer boundary equal to -GM/R
analytic = 2*pi*ggrav* (r ** 2 - 3*r[-1]**2)/3

err1 = analytic[len(r)/2] - phi[len(r)/2]
h1 = r[1] - r[0]
#Run a high resolution
r = np.linspace(1e7,1e9,2001)
rho = array([1]*len(r)) # Assume constant density 1
phi = FE(rho,r)[:,0] # Calculate phi
phi -= (phi[-1] + ggrav * 4*pi*r[-1]**2/3) # Set outer boundary equal to -GM/R
analytic = 2*pi*ggrav* (r ** 2 - 3*r[-1]**2)/3

err2 = analytic[len(r)/2] - phi[len(r)/2]
h2  = r[1] - r[0]

err_ratio = float(err2)/float(err1)
test = float(h2)/float(h1)
for n in np.arange(0,3,0.01):
    if test **n < err_ratio:
        print n
        break

"""

################
## Make Plots ##
################

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# Estimate Pi using the btter random numbers


# make the plot
analytic = 2*pi*ggrav* (r ** 2 - 3*r[-1]**2)/3
p1, = pl.semilogy(r,-1*phi,"r",linewidth=1.5)
p2, = pl.semilogy(r,-1*analytic,"b",linewidth=1.5)

# prepare x and y ranges

xmin = r[0]
xmax = r[-1]
ymin = -1*analytic[-1]
ymax = -1* analytic[0]

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("Radius",labelpad=10)
pl.ylabel("Phi",labelpad=12)

# legend
pl.legend( [p1,p2], ["Phi Numeric"," Phi Analytic"], loc=(0.5,0.85), frameon=False, labelspacing = 0.25 )

pl.savefig("3a.pdf")


