#!/usr/bin/env python

#####################
## Import Packages ##
#####################

import numpy as np
import scipy as sp
import matplotlib.pyplot as pl
from plot_defaults import *

###################
## Set Constants ##
###################

# global constants
ggrav = 6.67e-8
msun  = 1.99e33
radmax = 2.0e8 # 2000 km

# EOS parameters
# for white dwarfs:
polyG = 4.0/3.0
polyK = 1.244e15*0.5**polyG

#############################
## Define Helper Functions ##
#############################

### Integrators 
# tov_RHS(rad,h,p,rho,m):
# tov_integrate_FE(rad,h,p,rho,m):
# tov_integrate_RK2(rad,h,p,rho,m):
# tov_integrate_RK3(rad,h,p,rho,m):
# tov_integrate_RK4(rad,h,p,rho,m):

def get_rho(p):
    ''' Reverse EOS to get rho from p'''
    if p < 0:
        print p
    return (p/polyK) **(1/polyG)

def stellar_mass(npoints, radmax, integrator):
    ''' Returns the mass enclosed at the surface '''
    (radius, press, rho, mass,nsurf) = \
        stellar(npoints, radmax,integrator)
    return (mass[nsurf]/msun)

def stellar_radius(npoints, radmax, integrator):
    ''' Returns the mass enclosed at the surface '''
    (radius, press, rho, mass,nsurf) = \
        stellar(npoints, radmax,integrator)
    return (radius[nsurf]/1e5)

def self_convergance(radmax, integrator, n1,n2,n3):
    ''' Returns the self convergance factor of a specific integrator'''
    h3 = float(radmax)/n3
    h2 = float(radmax)/n2
    h1 = float(radmax)/n1
    yh3 = stellar_mass(n3,radmax,integrator)
    yh2 = stellar_mass(n2,radmax,integrator)
    yh1 = stellar_mass(n1,radmax,integrator)
    Q =  abs(yh3 - yh2)/ abs(yh2 - yh1)
    integ = str(integrator).split()[1][-3:]
    high = 8
    step = 0.01
    for n in np.arange(0,high,step):
        if (h3 ** n - h2**n) / (h2 **n - h1 **n) < Q:
            print "The self convergance rate for {} is {}".format(integ, n)
            return 
    print "The self convergance rate is higher than {}".format(high)
    return

######################
## Define Functions ##
######################

def tov_RHS(rad,p,rho,m):
    ''' Return the right hand side of the system of ODEs'''
    rhs = np.zeros(3)
    if p < 0:
        return rhs

    # RHS function
    rho = get_rho(p)

    if(rad > 1.0e-10):
        rhs[0] = -1*ggrav * m * rho/rad**2
        rhs[1] = 4* np.pi*rho*rad **2
        rhs[2] = rho

    else:
        rhs[0] = 0.0
        rhs[1] = 0.0
        rhs[2] = 1e10

    return rhs


### Integrators
def tov_integrate_FE(rad,h,p,rho,m):
    '''Forward-Euler Integrator'''

    new = np.zeros(3)
    old = np.zeros(3)
    old[0] = p
    old[1] = m
    old[2] = rho

    new = old + h * tov_RHS(rad,p,rho,m)
    
    # assign outputs
    pnew = new[0]
    mnew = new[1]
    
    return (pnew,mnew,rho)

def tov_integrate_RK2(rad,h,p,rho,m):
    '''Runge-Kutta2 Integrator'''
    new = np.zeros(3)
    k1  = np.zeros(3) #[p,m,rho]
    k2  = np.zeros(3)

    # Define k1, k2
    k1 = h * tov_RHS(rad,p,rho,m)
    k2 = h * tov_RHS(rad + h/2, p + 0.5 * k1[0],rho,m + 0.5*k1[1])
    

    # assign outputs
    pnew = p + k2[0]
    mnew = m + k2[1]
    rho = get_rho(pnew)
    
    return (pnew,mnew,rho)

def tov_integrate_RK3(rad,h,p,rho,m):
    '''Runge-Kutta3 Integrator'''

    new = np.zeros(2)
    k1  = np.zeros(2) #[p,m,rho]
    k2  = np.zeros(2) #[p,m,rho]
    k3  = np.zeros(2) #[p,m,rho]

    ## Define k1, k2, k3
    k1 = h * tov_RHS(rad,p,rho,m)
    k2 = h * tov_RHS(rad + h/2, p + 0.5* k1[0], rho, m + 0.5*k1[1])
    k3 = h * tov_RHS(rad+h, p - k1[0] + 2*k2[0],rho, m - k1[1] + 2*k2[1])

    # assign outputs
    pnew = p + (k1[0] + 4 * k2[0] + k3[0])/6
    mnew = m + (k1[1] + 4 * k2[1] + k3[1])/6
    rho = get_rho(pnew)
    
    return (pnew,mnew,rho)

def tov_integrate_RK4(rad,h,p,rho,m):
    '''Runge-Kutta4 Integrator'''
    new = np.zeros(3)
    k1  = np.zeros(3)
    k2  = np.zeros(3)
    k3  = np.zeros(3)
    k4  = np.zeros(3)

    # Define k1,k2,k3,k4
    k1 = h * tov_RHS(rad,p,rho,m)
    k2 = h * tov_RHS(rad + h/2, p + 0.5* k1[0],rho, m + 0.5*k1[1])
    k3 = h * tov_RHS(rad + h/2, p + 0.5* k2[0], rho, m + 0.5*k2[1])
    k4 = h * tov_RHS(rad + h, p + k3[0], rho, m + k3[1])

    # assign outputs
    pnew = p + (k1[0] + 2 * k2[0] + 2*k3[0] + k4[0])/6
    mnew = m + (k1[1] + 2 * k2[1] + 2*k3[1] + k4[1])/6
    rho  = get_rho(pnew)
    
    return (pnew,mnew,rho)


def stellar(npoints,radmax,integrator):
    ''' Takes a number of points, a maximum radius and a type of integrator
        and returns arrays of radius, pressure, density and mass enclosed as
        well as the index of the 'surface' '''
    radius = np.linspace(1e-10,radmax,npoints)
    dr = radius[1]-radius[0]

    # set up variables
    press = np.zeros(npoints)
    rho   = np.zeros(npoints)
    mass  = np.zeros(npoints)

    # set up central values
    rho[0]   = 1.0e10
    press[0] = polyK * rho[0]**polyG
    mass[0]  = 0.0

    # set up termination criterion
    press_min = 1.0e-8 * press[0]

    nsurf = 0
    for n in range(npoints-1):
        
        (press[n+1],mass[n+1], rho[n+1]) = integrator(radius[n],
                                                  dr,
                                                  press[n],
                                                  rho[n],mass[n])
        # check for termination criterion
        if(press[n+1] < press_min and nsurf==0):
            nsurf = n

        if(n+1 > nsurf and nsurf > 0):
            press[n+1] = press[nsurf]
            rho[n+1]   = rho[nsurf]
            mass[n+1]  = mass[nsurf]

    return (radius, press, rho, mass,nsurf)



##################################
## Perform Convergance Analysis ##
##################################

if(False):
    n = 1000
    print 'FE'
    print stellar_mass(n,radmax, tov_integrate_FE)
    print stellar_radius(n,radmax, tov_integrate_FE)
    self_convergance(radmax, tov_integrate_FE,101,201,401)
    self_convergance(radmax, tov_integrate_FE,1001,2001,4001)

    print 'RK2'
    print stellar_mass(n,radmax, tov_integrate_RK2)
    print stellar_radius(n,radmax, tov_integrate_RK2)
    self_convergance(radmax, tov_integrate_RK2,101,201,401)
    self_convergance(radmax, tov_integrate_RK2,1001,2001,4001)

    print 'RK3'
    print stellar_mass(n,radmax, tov_integrate_RK3)
    print stellar_radius(n,radmax, tov_integrate_RK3)
    self_convergance(radmax, tov_integrate_RK3,101,201,401)
    self_convergance(radmax, tov_integrate_RK3,1001,2001,4001)

    print 'RK4'
    print stellar_mass(n,radmax, tov_integrate_RK4)
    print stellar_radius(n,radmax, tov_integrate_RK4)
    self_convergance(radmax, tov_integrate_RK4,101,201,401)
    self_convergance(radmax, tov_integrate_RK4,1001,2001,4001)


###############
## Make Plot ##
###############

(radius, press, rho, mass, nsurf) = stellar(1000,radmax,tov_integrate_RK3)

radius = radius[0:nsurf] / 1e5
press  = press[0:nsurf] / press[nsurf]
rho    = rho[0:nsurf]   / rho[nsurf]
mass   = mass[0:nsurf]  / mass[nsurf]

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.16)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.95)
myfig.subplots_adjust(right=0.96)


p1, = pl.semilogy(radius,press,"b",linewidth=1.5)
p2, = pl.semilogy(radius,rho,"r",linewidth=1.5)
p3, = pl.semilogy(radius,mass,"k",linewidth=1.5)

# prepare x and y ranges

xmin = 1
xmax = 1.5e3
ymin = 1e-6
ymax = 1e10

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("Radius km",labelpad=10)
pl.ylabel("Ratio to Surface Value",labelpad=12)

# legend
pl.legend( (p1,p2,p3), ["Pressure","Density","Mass"], loc=(0.6,0.1), frameon=False, labelspacing = 0.25 )

pl.savefig("4.pdf")


