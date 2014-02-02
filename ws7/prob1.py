#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, sin,cos,sqrt,pi
import random
sr = random.SystemRandom()
import matplotlib.pyplot as pl
from plot_defaults import *


def mcPi(N):
    np.random.seed(1)
    n = 0.0
    for i in range(N):
        x = np.random.rand() - 0.5
        y = np.random.rand() - 0.5
        if x ** 2 + y ** 2 < 0.25:
            n += 1.0
    return 4 * n/N

def betterMcPi(N):
    n = 0.0
    for i in range(N):
        x = sr.random() - 0.5
        y = sr.random() - 0.5
        if x ** 2 + y ** 2 < 0.25:
            n += 1.0
    return 4 * n/N


### Using numpy random

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# Log the sigma base 10
n  = np.arange(10,10000,5)
est_pi = np.abs([mcPi(i)-pi for i in n])

# make the plot
p1, = pl.semilogy(n,est_pi,"r",linewidth=1.5)

# prepare x and y ranges
xmin = 10
xmax = 1000
ymin = 1e-4
ymax = 1.5

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("N",labelpad=10)
pl.ylabel("Error",labelpad=12)

# legend
pl.legend( [p1], ["Error"], loc=(0.5,0.80), frameon=False, labelspacing = 0.25 )

pl.savefig("1a.pdf")

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# Estimate Pi using the btter random numbers
est_pi = np.abs([betterMcPi(i)-pi for i in n])

# make the plot
pl.semilogy(n,est_pi,"r",linewidth=1.5)
# prepare x and y ranges

xmin = 10
xmax = 1000
ymin = 1e-4
ymax = 0.5

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("N",labelpad=10)
pl.ylabel("Error",labelpad=12)

# legend
pl.legend( [p1], ["Error"], loc=(0.5,0.85), frameon=False, labelspacing = 0.25 )

pl.savefig("1b.pdf")

