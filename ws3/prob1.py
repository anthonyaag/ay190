#!/usr/bin/python 

import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *

## Functions for making plots
def f1(x):
    return np.sin(x)

def f2(x):
    return x * np.sin(x)

# Interpolation rules
def midpoint(f,n,a,b):
    x = np.linspace(a, b, n)
    result = 0
    for i in range(len(x) - 1):
        result += (x[i + 1] - x[i]) * (f((x[i] + x[i + 1]) / 2))
    return result

def trapezoidal(f,n,a,b):
    x = np.linspace(a, b, n)
    result = 0
    for i in range(len(x) - 1):
        result += (x[i + 1] - x[i]) * (f(x[i]) + f(x[i+1])) / 2.0
    return result

def simpsons(f,n,a,b):
    x = np.linspace(a, b, n)
    result = 0
    for i in range(len(x) - 1):
        result += ((x[i + 1] - x[i])/6.0) * (f(x[i]) +\
            4*f((x[i] + x[i+1])/2.0) + f(x[i+1]))
    return result

#Generate data
ns = np.arange(2,900,1)
midp = []
trap = []
simp = []

#We subtract two when Integrating f1 from 0 to pi since the analytic answer is 2
for n in ns:
    midp.append( abs(2.0 - midpoint(f1,n,0,np.pi))    )
    trap.append( abs(2.0 - trapezoidal(f1,n,0,np.pi)) )
    simp.append( abs(2.0 - simpsons(f1,n,0,np.pi))    )

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# make the plot
p1, = pl.plot(1.0/ns,midp,"r",linewidth=1.5)
p2, = pl.plot(1.0/ns,trap,"b",linewidth=1.5)
p3, = pl.plot(1.0/ns,simp,"g",linewidth=1.5)

# prepare x and y ranges
xmin = 0
xmax = 0.5
ymin = 10** -12
ymax = 5

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()
ax.set_yscale('log')

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("Step Size",labelpad=10)
pl.ylabel("Absolute Error",labelpad=12)

# legend
pl.legend( (p1,p2,p3), ("Midpoint","Trapezoidal","Simpsons"),
           loc=(0.5,0.10), frameon=False, labelspacing = 0.25 )

pl.savefig("1a.pdf")

# b)

#Generate data
midp = []
trap = []
simp = []

#We subtract two when Integrating f1 from 0 to pi since the analytic answer is 2
for n in ns:
    midp.append( abs(np.pi - midpoint(f2,n,0,np.pi))    )
    trap.append( abs(np.pi - trapezoidal(f2,n,0,np.pi)) )
    simp.append( abs(np.pi - simpsons(f2,n,0,np.pi))    )

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# make the plot
p1, = pl.plot(1.0/ns,midp,"r",linewidth=1.5)
p2, = pl.plot(1.0/ns,trap,"b",linewidth=1.5)
p3, = pl.plot(1.0/ns,simp,"g",linewidth=1.5)

# prepare x and y ranges
xmin = 0
xmax = 0.5
ymin = 10** -12
ymax = 5

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()
ax.set_yscale('log')

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("Step Size",labelpad=10)
pl.ylabel("Absolute Error",labelpad=12)

# legend
pl.legend( (p1,p2,p3), ("Midpoint","Trapezoidal","Simpsons"),
           loc=(0.5,0.10), frameon=False, labelspacing = 0.25 )

pl.savefig("1b.pdf")
