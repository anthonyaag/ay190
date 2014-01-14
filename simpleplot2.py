#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.14)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# read in the data
data = np.loadtxt("infile3.txt",comments="#")

# slice them
x = data[:,0]
y1 = data[:,1]
y2 = data[:,2]
y3 = np.sin(x) ** 2 * np.cos(x)

# make the plot
p1, = pl.plot(x,y1,"r",linewidth=2.5)
p2, = pl.plot(x,y2,"b",linewidth=2.5)
p3, = pl.plot(x,y3,"g",linewidth=2.5)

# prepare x and y ranges
xmin = 0.0
xmax = 19.0
ymin = -1.0
ymax = 1.0

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set locators of tick marks
xminorLocator = pl.MultipleLocator(1)
xmajorLocator = pl.MultipleLocator(5)
yminorLocator = pl.MultipleLocator(0.1)
ymajorLocator = pl.MultipleLocator(0.5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,2.0,0.75*2.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("X",labelpad=15)
pl.ylabel("Y",labelpad=-5)


# legend
# loc is the location of the legend in a
# coordinate system from (0,0) to (1,1)
# frameon=False turns of the stupid box around the legend
pl.legend( (p1,p2,p3), ("y1","y2","y3"), 
           loc=(0.4,0.025), frameon=False, labelspacing = 0.25 )

pl.text(0.14,0.88,"LaTeX: $\\Phi$",fontsize=30,
        horizontalalignment="left",rotation="horizontal",
        transform=ax.transAxes)

# uncomment the below line to get screen display
# pl.show()
# comment the below line to save as a pdf
pl.savefig("simpleplot2.pdf")
