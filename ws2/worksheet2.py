#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import *
from numpy import float32, arange, array

# 1
oneThird         = float32(1.0/3)
thirteenThirds   = float32(13.0/3)
fourThirds       = float32(4.0/3)

def recursive_thirds(n):
    ret = [float32(1),oneThird]
    for i in range(2,n):
        ret.append(thirteenThirds * ret[-1] - fourThirds*ret[-2])
    return ret[:n]

def true_thirds(n):
    return [(1.0/3) ** x for x in range(n)] # Use max precision here

def absolute_error(n):
    error = []
    recursive = recursive_thirds(n)
    true = true_thirds(n)
    for i in range(n):
        error.append(true[i] - recursive[i])
    return error

def relative_error(n):
    error = []
    recursive = recursive_thirds(n)
    true = true_thirds(n)
    for i in range(n):
        error.append((true[i] - recursive[i])/true[i])
    return error

# At n = 15
# absolute error is -0.90082688410173606
# relative error is -4308627.0610251995

# print recursive_thirds(15)

# 2

"""

def f(x):
    return x ** 3 - 5 * x ** 2 + x

def fp(x):
    return 3 * x**2 - 10 *x + 1

def forward_diff(func,begin,end,h):
    ret = []
    for i in arange(begin,end,h):
        ret.append((func(i + h) - func(i))/h)
    return ret

def central_diff(func,begin,end,h):
    ret = []
    for i in arange(begin,end,h):
        ret.append((func(i + h) - func(i-h))/(2 *h))
    return ret

# Discritize [-2,6]
h1 = 0.1
h2 = h1/2

x1 = arange(-2,6, h1)
x2 = arange(-2,6, h2)
forward1 = array(forward_diff(f,-2,6,h1))
forward2 = array(forward_diff(f,-2,6,h2))
central1 = array(central_diff(f,-2,6,h1))
central2 = array(central_diff(f,-2,6,h2))
fp1 = array(map(fp,x1))
fp2 = array(map(fp,x2))

ferr1 = np.abs(forward1 - fp1)
ferr2 = np.abs(forward2 - fp2)
cerr1 = np.abs(central1 - fp1)
cerr2 = np.abs(central2 - fp2)

# 3 analytic exercise
"""
# 4 Interpolation: Cepheid Lightcurve
# a

# Need to import multiply function
from operator import mul

time = np.array([0.0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0])
mag  = np.array([0.302,0.185,0.106,0.093,0.24,0.579,0.561,0.468,0.302])

def lagrange(x,y,ran):
	'''Returns the values of the Lagrange interpolation using all of the
 	   points for given x and y lists over a given range'''
	# First generate a list of denominators
	denom = []
	for i in x:
		denom.append(reduce(mul,[i - j for j in x if j != i],1))
	# Multiply them by the values f(x_i) aka y
	denom = np.array(y) / np.array(denom)
	# get the return values
	ret = []
	for i in ran:
		numer = []
		for j in x:
			numer.append(reduce(mul,[i - k for k in x if j != k],1))
		ret.append(np.dot(np.array(numer),np.array(denom)))
	return ret	




# b

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

def quad_interp(x,y,ran):
    ''' Takes input x and f(x) = y in lists, assumed to be sorted such that x
    is a sorted list in increasing order. Takes input, ran, a list of 
    points also assumed to be sorted in increasing order.'''
    ret = []
    x_iter = 0
    ran_iter = 0
    while ran_iter < len(ran): # while not done with all of the given points
        # Check to see if you're in the correct bin
        if ran[ran_iter] < x[x_iter + 1] or x_iter == len(x) - 3:
            est = ran[ran_iter]
            (x_i,x_i1,x_i2) = x[x_iter:x_iter + 3]
            (y_i,y_i1,y_i2) = y[x_iter:x_iter + 3]
            temp =(est - x_i1)*(est - x_i2)* y_i /((x_i  - x_i1) *(x_i  - x_i2))
            temp+=(est - x_i) *(est - x_i2)* y_i1/((x_i1 - x_i) * (x_i1 - x_i2))
            temp+=(est - x_i) *(est - x_i1)* y_i2/((x_i2 - x_i) * (x_i2 - x_i1))
            ret.append(temp)
            ran_iter += 1
        else:
            x_iter += 1
    return ret



# 5

def hermite_interp(x,y,ran):
    ''' Takes input x and f(x) = y in lists, assumed to be sorted such that x
    is a sorted list in increasing order. Takes input, ran, a list of 
    points also assumed to be sorted in increasing order.'''
    ret = []
    x_iter = 0
    ran_iter = 0
    while ran_iter < len(ran): # while not done with all of the given points
        # Check to see if you're in the correct bin
        if ran[ran_iter] < x[x_iter + 1] or x_iter == len(x) - 3:
            est = ran[ran_iter]
            (x_i,x_i1,x_i2) = x[x_iter:x_iter + 3]
            (y_i,y_i1,y_i2) = y[x_iter:x_iter + 3]
            fp_i  = (y_i1 - y_i) /(x_i1 - x_i)
            fp_i1 = (y_i2 - y_i1)/(x_i2 - x_i1)
            z = (est - x_i) / (x_i1 - x_i)
            phi0z  = 2 * z**3 - 3* z**2 + 1
            phi0mz = 2 * (1-z)**3 - 3* (1-z)**2 + 1
            phi1z  = z**3 - 2 * z**2 + z
            phi1mz = (1-z)**3 - 2 * (1-z) **2 + (1- z)
            ret.append(y_i * phi0z + y_i1 * phi0mz + fp_i * (x_i1 - x_i) * \
                       phi1z - fp_i1 * (x_i1 - x_i) * phi1mz)
            ran_iter += 1
        else:
            x_iter += 1
    return ret


''' Natural Cubic Spline Interpolation can be acchieved using 
    scipy.interpolate.interp1d(x,y,kind = 'cubic')'''

from scipy.interpolate import interp1d

def nat_cubic_spline_interp(x,y,ran):
    f = interp1d(x,y,kind = 'cubic')
    return f(np.array(ran))

x = arange(0,1.01,0.01)
lin_inter = lin_interp(time,mag,x)
quad_inter = quad_interp(time,mag,x)
herm_inter = hermite_interp(time,mag,x)
nat_cub_inter = nat_cubic_spline_interp(time,mag,x)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.13)
myfig.subplots_adjust(bottom=0.11)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)


# make the plot
p1, = pl.plot(x,lin_inter,"r",linewidth=1.5)
p2, = pl.plot(x,quad_inter,"b",linewidth=1.5)
p3, = pl.plot(x,herm_inter,"g",linewidth=1.5)
p4, = pl.plot(x,nat_cub_inter,"c",linewidth=1.5)
p5, = pl.plot(time,mag,"ko",linewidth = 1.5)

# prepare x and y ranges
xmin = -0.2
xmax = 1.2
ymin = 0
ymax = 0.8

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set locators of tick marks
xminorLocator = pl.MultipleLocator(.1)
xmajorLocator = pl.MultipleLocator(.5)
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
pl.xlabel("X",labelpad=-20)
pl.ylabel("Y",labelpad=12)



# legend
# loc is the location of the legend in a
# coordinate system from (0,0) to (1,1)
# frameon=False turns of the stupid box around the legend
pl.legend( (p1,p2,p3,p4,p5), ("Linear","Quadratic","Cubic Hermite","Natrual Cubic", "Data"), 
           loc=(0.05,0.6), frameon=False, labelspacing = 0.25 )

#pl.text(0.14,0.88,"LaTeX: $\\Phi$",fontsize=30,
#        horizontalalignment="left",rotation="horizontal",
#        transform=ax.transAxes)

# uncomment the below line to get screen display
# pl.show()
# comment the below line to save as a pdf
pl.savefig("herm_natrual_cubic.pdf")
