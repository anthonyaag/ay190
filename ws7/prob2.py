#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, sin,cos,sqrt,pi
import matplotlib.pyplot as pl
from plot_defaults import *

def sameBirthday(N):
    birthdays = np.random.randint(365,size = N)
    if len(birthdays) != len(set(birthdays)):
        return 1 # At least two people have same birthday

trials = 5000
def birthdayParadox(num_people):
    counter = 0.0
    for i in range(trials):
        birthdays = np.random.randint(365,size = num_people)
        if len(birthdays) != len(set(birthdays)):
            counter += 1.0 # At least two people have same birthday
    return counter/trials

                    
# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# Log the sigma base 10
people  = np.arange(1,100)
prob    = [birthdayParadox(i) for i in people]
print prob

# make the plot
p1, = pl.plot(people,prob,"r",linewidth=1.5)

# prepare x and y ranges
'''
xmin = 2
xmax = 40
ymin = 0
ymax = 0.7

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])'''
# get axis object
ax = pl.gca()

# set locators of tick marks
xminorLocator = pl.MultipleLocator(10)
xmajorLocator = pl.MultipleLocator(50)
yminorLocator = pl.MultipleLocator(0.05)
ymajorLocator = pl.MultipleLocator(0.25)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("Number of People",labelpad=10)
pl.ylabel("Probability of Shared Birthday",labelpad=12)

# legend
pl.legend( [p1], ["Probability"],\
             loc=(0.5,0.80), frameon=False, labelspacing = 0.25 )

pl.savefig("2a.pdf")         
    
