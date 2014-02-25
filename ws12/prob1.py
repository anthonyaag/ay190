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


# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# Estimate Pi using the btter random numbers


# make the plot
p1, = pl.loglog(radius,density,"r",linewidth=1.5)
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
pl.legend( [p1], ["Density"], loc=(0.5,0.85), frameon=False, labelspacing = 0.25 )

pl.savefig("1.pdf")

