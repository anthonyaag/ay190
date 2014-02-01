#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, exp, sin,cos,sqrt
from timeit import timeit
import matplotlib.pyplot as pl
from plot_defaults import *

# Problem b & c
mytimes = []
nptimes = []
trials = 100

for N in range(10,trials):
    mytimes.append( timeit("dft(x)", number = 20, \
        setup = "from proba import dft; import pylab; x = pylab.randn(%d)" %N))
    nptimes.append( timeit("np.fft.fft(x)", number = 20, \
        setup = "import numpy as np;    import pylab; x = pylab.randn(%d)" %N))

# Get x range
N = range(10,trials)
mytimes = array(mytimes)/20
nptimes = array(nptimes)/20
# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# make the plot
p1, = pl.semilogy(N,mytimes,"r",linewidth=1.5)
p2, = pl.semilogy(N,nptimes,"b",linewidth=1.5)

# prepare x and y ranges
xmin = 10
xmax = 100
ymin = 1e-6
ymax = 1e-1

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
pl.ylabel("Transform Time",labelpad=12)

# legend
pl.legend( (p1,p2) , ("My DFT", "FFT") , loc=(0.2,0.80), \
            frameon=False, labelspacing = 0.25 )

pl.savefig("1b.pdf")
pl.show()
