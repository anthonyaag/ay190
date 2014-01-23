#!/usr/bin/python 

import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *
from scipy import special as sp
from numpy import exp 

#   [laguerre_roots, laguerre_weights] = sp.l_roots(n,0)
#   [legendre_roots, legendre_weights] = sp.p_roots(n,0)

# Problem a

#Integrand is x^2 / (exp(x) +1) so we exclude the weight function exp(-x) to get

def f(x):
    x = np.array(x) # make sure that the input is a numpy array
    return x** 2 / (1 + exp(-x)) 

def gauss_laguerre(f,n):
    [laguerre_roots, laguerre_weights] = sp.l_roots(n,0)
    return sum( f(laguerre_roots) * laguerre_weights )

# Try a few nodes

values = []
degree = range(2,50)
for i in degree:
    values.append(gauss_laguerre(f,i))

first_dif = values[0] - values[1]
last_dif = values[-2] - values[-1]

print first_dif, last_dif
#first_dif = 0.0247956061085
#last_dif =  4.19664303308e-14
print values[-1]
# last value 1.80308535474

# Get the error from most accurate term
error = abs(np.array(values) - values[-1])

### Now get a plot
# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# make the plot
p1, = pl.plot(degree,error,"r",linewidth=1.5)

# prepare x and y ranges
xmin = 0
xmax = 50
ymin = 10e-15
ymax = 1

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
pl.xlabel("Nodes",labelpad=10)
pl.ylabel("Absolute Error",labelpad=12)

# legend
pl.legend( [p1], ["Gauss-Laguerre Error"], loc=(0.3,0.8), frameon=False, labelspacing = 0.25 )

pl.savefig("2a.pdf")

# Problem b

def g(x):
    x = np.array(x) # make sure that the input is a numpy array
    return x **2 / ( exp(x) + 1 )

def gauss_legendre(f,n,a,b):
    [legendre_roots, legendre_weights] = sp.p_roots(n,0)
    trans = float(b-a)/2
    return trans * sum( \
                f(legendre_roots * trans + float(a+b)/2) * legendre_weights )

def step_gauss_leg(f,deg,A,B,step):
    steps = np.arange(A,B+step,step)
    ret = 0
    for i in range(len(steps) -1):
        ret += gauss_legendre(f,deg,steps[i], steps[i+1])
    return ret

values = []
degree = range(2,50)
for i in degree:
    values.append(step_gauss_leg(g,i,0,150,5))

first_dif = values[0] - values[1]
last_dif = values[-2] - values[-1]

print first_dif, last_dif
#first_dif = -0.0772043590512
#last_dif =  -3.99680288865e-15
print values[-1]
# last value 1.80308535474

# Again we see that it clearly convereges to the correct value. Now let us
# determine dn/dE

### Now get a plot
A = 0
B = 150
step = 0.1
degree = 8
bins = np.arange(A,B+step,step)
plot_bins = []


dn = []
for i in range(len(bins)-1):
    dn.append(gauss_legendre(g,degree,bins[i],bins[i+1])/step)
    plot_bins.append((bins[i] + bins[i+1])/2)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# make the plot
p1, = pl.plot(plot_bins,dn,"r",linewidth=1.5)

# prepare x and y ranges
xmin = -5
xmax = 150
ymin = -0.1
ymax = 0.5

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set locators of tick marks
xminorLocator = pl.MultipleLocator(10)
xmajorLocator = pl.MultipleLocator(50)
yminorLocator = pl.MultipleLocator(0.1)
ymajorLocator = pl.MultipleLocator(0.5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("Energy",labelpad=10)
pl.ylabel("dn/dE",labelpad=12)

# legend
pl.legend( [p1], ["dn/dE"], loc=(0.5,0.80), frameon=False, labelspacing = 0.25 )

pl.savefig("2b.pdf")

    
