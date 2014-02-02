#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, sin,cos,sqrt,pi
import matplotlib.pyplot as pl
from plot_defaults import *


def mcint(f,a,b,n,trials):
    ret = 0
    for i in range(trials):    
        # assumes that f is monotonically icreasing over (a,b) and b>a>0
        a1 = (b - a) * f(a)
        trialx = np.random.rand(n) * (b-a) + a
        trialy = np.random.rand(n) * (f(b) - f(a)) + f(a)
        count = 0
        for i in range(n):
            if  trialy[i] < f(trialx[i]):
                count += 1
        a2 = count * (f(b) - f(a)) * (b -a) / float(n)
        ret += abs((a1 + a2) - 22.0/3.0)
    return ret/trials

def f(x):
    return x**2 +1


### Using numpy random

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# Log the sigma base 10
n = np.arange(40,5000,40)
error = abs(np.array([mcint(f,2,3,i,200) for i in n]))

p1, = pl.semilogy(n,error,"r",linewidth=1.5)

# prepare x and y ranges
xmin = 10
xmax = 5000
ymin = 10**-2
ymax = 1

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

pl.savefig("3.pdf")
