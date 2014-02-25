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
x = np.arange(1e7,1e9,1e6)
interpRho = lin_interp(radius,density,x)


# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# Estimate Pi using the btter random numbers


# make the plot
p1, = pl.loglog(x,interpRho,"r",linewidth=1.5)
# prepare x and y ranges

xmin = 1e7
xmax = 1e9
ymin = 1e4
ymax = 2e10

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
pl.ylabel("Density",labelpad=12)

# legend
pl.legend( [p1], ["Density"], loc=(0.5,0.85), frameon=False,labelspacing = 0.25)

pl.savefig("2.pdf")

