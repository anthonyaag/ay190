#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, sin,cos,sqrt
import matplotlib.pyplot as pl


# Earth Parameters

T = 365.25635           # Period in Days
w = 0.01720212477395557 # Radial frequency in Days ^ -1
a = 1.496 * 10 ** 6     # semi-major axis in km
e = 0.0167              # ecentricity

def f(E,wt,e):
    return E - wt - e*sin(E)

def problema(f,t,E_guess,h):
    diff  = float('inf')
    E = E_guess
    counter = 0
    # Use secant method to find the root within relative error 10 ** -10
    while (abs(diff/E) > 10.0 ** -10):
        diff =  \
            float(f(E,w*t,e)) * h / (f(E + h/2,w*t,e) - f(E - h/2,w*t,e)) 
        E -= diff
        counter += 1
    return (a*cos(E), sqrt(a**2 - (e*a) ** 2) * sin(E),counter)

print problema(f,91,1,0.01)
# t = 91 days
# (x,y,steps) =  (-16898.400071799013, 1495695.9461840717,4)

print problema(f,182,1,0.01)
# t = 182 days
# (x,y,steps) =   (-1495915.5037150327, 15897.648882948801,4)

print problema(f,273,1,0.01)
# t = 273 days
# (x,y,steps) =  (-49209.341755800378, -1494981.9248908658,4)


## Something bad happend and now the Earth is highly ecentric 
e = 0.999999999

print problema(f,91,1,0.01)
# t = 91 days
# (x,y,steps) =  (-1004146.3196859331, 49.592461055262397, 7)

print problema(f,182,1,0.01)
# t = 182 days
# (x,y,steps) =  (-1495978.1642482972, 0.36147500183165632, 7)


print problema(f,273,1,0.01)
# t = 273 days
# (x,y,steps) =  (-1018362.0491944792, -49.009278030577299, 7)



