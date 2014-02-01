#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, sin,cos,sqrt
import matplotlib.pyplot as pl
from plot_defaults import *
from math import log

# Part a  -- Read the file
f = open("m_sigma_table.dat",'r')

num_lines = 0
name = []
zshift = []
sigma = []
e_sigma = []
n_sigma = []
FWHM = []
e_FWHM = []
logL = []
e_logL = []
logM = []
E_logM = []
e_logM = []
for line in f.readlines():
    num_lines += 1
    temp = line.split()
    
    #Names
    name.append(temp[0] + ' ' + temp[1])
    #zshift
    zshift.append(float(temp[2]))
    #sigma and e_sigma
    sigma.append(float(temp[3]))
    e_sigma.append(float(temp[4]))

    #n_sigma
    if temp[5] == 'a':
        n_sigma.append('a')
    elif temp[5] == 'b':
        n_sigma.append('b')
    else:
        n_sigma.append(' ')
        temp = temp[0:5] + [' '] + temp[5:]

    #FWHM & e_FWHM
    FWHM.append(int(temp[6]))
    e_FWHM.append(int(temp[7]))

    #logL & elogL
    logL.append(float(temp[8]))
    e_logL.append(float(temp[9]))

    #logM, E_logM, & e_logM
    logM.append(float(temp[10]))
    E_logM.append(float(temp[11]))
    e_logM.append(float(temp[12]))
sigma = np.array(sigma)
e_sigma = np.array(e_sigma)

f.close()
    
# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)

# Log the sigma base 10
logSigma = [log(x,10) for x in sigma]
e_logSigma = np.array(e_sigma)/np.array(sigma)

# make the plot
p1, = pl.plot(logSigma,logM,"ko",linewidth=1.5)

# prepare x and y ranges
xmin = 1.4
xmax = 2.5
ymin = 4
ymax = 9

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set locators of tick marks
xminorLocator = pl.MultipleLocator(.1)
xmajorLocator = pl.MultipleLocator(.5)
yminorLocator = pl.MultipleLocator(0.5)
ymajorLocator = pl.MultipleLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("$log_{10}\\sigma_*$",labelpad=10)
pl.ylabel("$log_{10} M_{BH}$",labelpad=12)

# legend
pl.legend( [p1], [""], loc=(0.5,0.80), frameon=False, labelspacing = 0.25 )

pl.savefig("1a.pdf")

    
### Part b - Linear regression

def linreg(x,y):
    # Takes two equal length arrays 
    # Assumes y = a1 + a2*x and returns [a1,a2]
    x = np.array(x)
    y = np.array(y)
    a1 = (sum(y)*sum(x**2) - sum(x) * sum(x*y)) / ( sum(x**2) - sum(x) **2)
    a2 = (sum(y * x) - sum(x) * sum(y)) / ( sum(x**2) - sum(x) **2)
    return (a1, a2)

# prepare x and y ranges for best picture
'''
xmin = 1.4
xmax = 2.5
ymin = 4
ymax = 9
'''
# prepare x and y ranges for best comparison
xmin = 1.4
xmax = 2.75
ymin = 4.5
ymax = 10

# Do the regression 
intercept, slope = linreg(logSigma,logM)
regX = np.linspace(xmin,xmax,500)
regY = np.array([intercept + slope * x for x in regX])

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=.91)

# make the plot
p1, = pl.plot(logSigma,logM,"ko",linewidth=1.5)
p2, = pl.plot(regX,regY,"r",linewidth=1.5)


# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set locators of tick marks
xminorLocator = pl.MultipleLocator(.05)
xmajorLocator = pl.MultipleLocator(.25)
yminorLocator = pl.MultipleLocator(0.25)
ymajorLocator = pl.MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("$log_{10}\\sigma_*$",labelpad=10)
pl.ylabel("$log_{10} M_{BH}$",labelpad=12)

# legend
pl.legend( [p1,p2], ["data","linear fit"], \
        loc=(0.5,0.80), frameon=False, labelspacing = 0.25 )

pl.savefig("1b.pdf")


### Part c - Linear regression with errors

def linregErr(x,y,sx,sy):
    ''' Assumes inputs of lists of equal lengths. Returns (a1,a2,sa1,sa2) '''
    #First get reasonable approx of slope(dy/dx) using no error linear fit
    intercept, m = linreg(x,y)
    x = np.array(x)
    y = np.array(y)
    sx = np.array(sx)
    sy = np.array(sy)
    err = sqrt((m*sx) **2 + sy **2)
    s = 1/err**2
    S = sum(s)
    a1 = (sum(y * s) * sum(x **2 * s) - sum(x * s) * sum(x * y * s)) / \
            (S * sum(x **2 * s) - sum(x * s) ** 2)
    a2 = (S * sum(y * x *s ) - sum(x * s) * sum(y * s)) / \
            (S * sum(x **2 * s) - sum(x * s) ** 2) 
    sa1 = sqrt( sum(x **2 * s) / (S * sum(x **2 * s) - sum(x *s) **2) )
    sa2 = sqrt( S/ (S * sum(x **2 * s) - sum(x *s) **2) )
    return(a1,a2,sa1,sa2)


# prepare x and y ranges for best picture
'''
xmin = 1.4
xmax = 2.5
ymin = 4
ymax = 9
'''
# prepare x and y ranges for best comparison
xmin = 1.2
xmax = 2.75
ymin = 4
ymax = 10

# Do the regression 
intercept, slope, si, ss = linregErr(logSigma,logM,e_logSigma,e_logM)
regX = np.linspace(xmin,xmax,500)
regY = np.array([intercept + slope * x for x in regX])
regYlow = np.array([intercept - si + (slope - ss) * x for x in regX])
regYhigh = np.array([intercept + si + (slope + ss) * x for x in regX])

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=.91)

# make the plot
p1, = pl.plot(logSigma,logM,"ko",linewidth=1.5)
p2, = pl.plot(regX,regY,"r",linewidth=1.5)
p3, = pl.plot(regX,regYlow,"b",linewidth=1.5)
p4, = pl.plot(regX,regYhigh,"g",linewidth=1.5)

pl.errorbar(logSigma,logM, yerr = e_logM, xerr = e_logSigma, fmt = "ko", linewidth = 0.5)

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set locators of tick marks
xminorLocator = pl.MultipleLocator(.05)
xmajorLocator = pl.MultipleLocator(.25)
yminorLocator = pl.MultipleLocator(0.25)
ymajorLocator = pl.MultipleLocator(1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("$log_{10}\\sigma_*$",labelpad=10)
pl.ylabel("$log_{10} M_{BH}$",labelpad=12)

# legend
pl.legend( [p1,p2,p3,p4], ["Data","Center Fit","Low Fit","High Fit"], \
        loc=(0.07,0.68), frameon=False, labelspacing = 0.25 )

pl.savefig("1c.pdf")








