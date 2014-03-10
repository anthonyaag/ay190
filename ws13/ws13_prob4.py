#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as pl
import mpl_toolkits.mplot3d as mpl3d
from plot_defaults import *

def NbodyRHS(u,mass):
    rhs = np.zeros((len(mass),6))#6vector for each mass first 3 vel, next 3 acc
    rhs[:,0:3] = u[:,3:6] # Store velocities
    for i,mi in enumerate(mass): # For each object
        for j,mj in enumerate(mass): # For each object calculate acceleration
            if i != j: # Objects don't effect themselves
                xij = u[i,0:3] - u[j,0:3]
                rhs[i,3:6] += xij * mj/(np.dot(xij,xij) ** 1.5)#Effect of j on i
    rhs[:,3:6] *= -ggrav
    return rhs

def NbodyRK4(u,mass,dt):
    '''Runge-Kutta4 Integrator'''

    # Define k1,k2,k3 & k4
    k1 = dt * NbodyRHS(u,mass)
    k2 = dt * NbodyRHS(u + k1/2.0, mass)
    k3 = dt * NbodyRHS(u + k2/2.0, mass)
    k4 = dt * NbodyRHS(u + k3, mass)
    return u + (k1 + 2*k2 + 2*k3 + k4)/6.0

def TotalEnergy(u,mass):
    KE = 0
    PE = 0
    ## Add up kinetic energy of each particle then divide by two.
    for i,speed in enumerate(u[:,3:6]):
        KE += mass[i] * np.dot(speed,speed)
    KE /= 2.0
    for i,mi in enumerate(mass):
        for j,mj in enumerate(mass):
            if i < j: # account for double counting pairs
                xij = u[i,0:3] -u[j,0:3]
                PE -= mi * mj/np.sqrt(np.dot(xij, xij))
    PE *= ggrav 
    return PE + KE
           
   
# main program

# global constants
ggrav = 6.67e-8
msun  = 1.99e33
seconds_per_year = 24.*3600*365 # roughly
cm_per_pc = 3.1e18
distance_to_sgrAstar = 8e3 * cm_per_pc

final_data_file = "final_positions.asc"
plt.ion()
# system parameters
distance_unit_to_cm = 0.04 * cm_per_pc
time_unit_to_s = seconds_per_year
mass_unit_to_g = msun
Nsteps = 1e4
t0 = 0
t1 = 100*seconds_per_year
initial_data_file = "sgrAstar.asc"
(x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)


# convert from unitis in initial data file to cgs
x *= distance_unit_to_cm
y *= distance_unit_to_cm
z *= distance_unit_to_cm
vx *= distance_unit_to_cm / time_unit_to_s
vy *= distance_unit_to_cm / time_unit_to_s
vz *= distance_unit_to_cm / time_unit_to_s
mass *= mass_unit_to_g

xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
zmin = np.amin(z)
zmax = np.amax(z)
rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))

# use a single state vector to simplify the ODE code
# indices:
# u[:,0] = x
# u[:,1] = y
# u[:,2] = z
# u[:,3] = vx
# u[:,4] = vy
# u[:,5] = vz
u = np.array((x,y,z,vx,vy,vz)).transpose()

def gen_EnergyErr(u,t0,t1,Nsteps):
    Times = np.linspace(t0,t1,Nsteps+ 1)
    dt = Times[1] - Times[0]
    Energies = []
    ret_Times = []
    for it,time in enumerate(Times):
        u = NbodyRK4(u,mass,dt)
         
        if it % max(1,Nsteps/100) == 0:
          ret_Times.append(time)
          Energy = TotalEnergy(u,mass)
          Energies.append(Energy)

          radius = np.sqrt(np.dot(u[0,0:3] - u[1,0:3],u[0,0:3] - u[1,0:3]))
          print "it = %d, time = %g years, radius = %g, energy = %g" % \
                (it, time / seconds_per_year,radius,Energy)
          '''
                    plt.clf()
          fig = plt.gcf()
          ax = mpl3d.Axes3D(fig)
          ax.scatter(u[:,0],u[:,1],u[:,2])
          ax.set_xlim((-rmax,rmax))
          ax.set_ylim((-rmax,rmax))
          ax.set_zlim((-rmax,rmax))
          plt.draw()
          '''
    ret_Times = np.array(ret_Times)
    Energies = np.array(Energies)
    return (ret_Times/seconds_per_year, np.abs(Energies - Energies[0]), dt)

np.savetxt(final_data_file, u)    

(Times1, EnergiesErr1,dt1) = gen_EnergyErr(u,t0,t1,Nsteps)
'''
(Times2, EnergiesErr2,dt2) = gen_EnergyErr(u,t0,t1,Nsteps * 2)
(Times3, EnergiesErr3,dt3) = gen_EnergyErr(u,t0,t1,Nsteps * 4)

# set up the figure and control white space
myfig = pl.figure(figsize=(10,8))
myfig.subplots_adjust(left=0.17)
myfig.subplots_adjust(bottom=0.15)
myfig.subplots_adjust(top=0.97)
myfig.subplots_adjust(right=0.975)


# make the plot
p1, = pl.semilogy(Times1,EnergiesErr1,"r",linewidth=1.5)
p2, = pl.semilogy(Times2,EnergiesErr2,"b",linewidth=1.5)
p3, = pl.semilogy(Times3,EnergiesErr3,"g",linewidth=1.5)
# prepare x and y ranges

xmin = 1
xmax = 4
ymin = 1e34
ymax = 1e40

# set axis parameters
pl.axis([xmin,xmax,ymin,ymax])
# get axis object
ax = pl.gca()

# set the custom tick sizes we like
# these functions are defined in plot_defaults
set_ticklines(ax,1.0,0.75*1.0)
set_tick_sizes(ax, 13, 7)

# label the axes
pl.xlabel("Time (years)",labelpad=10)
pl.ylabel("Energy Error (ergs)",labelpad=12)

# legend
pl.legend( [p1,p2,p3], ["dt = " + str(dt1) + "s", "dt = " + str(dt2)  + "s", "dt = " + str(dt3)  + "s"], loc=(0.45,0.1), frameon=False,labelspacing = 0.25)

pl.savefig("4c.pdf")



# output result
file_header = "1:x 2:y 3:z 4:vx 5:vy 6:vz 7:mass"
'''

