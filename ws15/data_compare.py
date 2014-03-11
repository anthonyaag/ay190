#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as pl
from plot_defaults import * 

a = ['050','105','150','195']

for suffix in a:
  py = np.loadtxt('py_0' + suffix + '.dat').transpose()[1:3]
  f  = np.loadtxt('f_0' + suffix + '.dat').transpose()[1:3]

  # set up the figure and control white space
  myfig = pl.figure(figsize=(10,8))
  myfig.subplots_adjust(left=0.17)
  myfig.subplots_adjust(bottom=0.15)
  myfig.subplots_adjust(top=0.97)
  myfig.subplots_adjust(right=0.975)

  # make the plot
  p1, = pl.plot(py[0],py[1],"r",linewidth=1.5)
  p2, = pl.plot(f[0]-0.5,f[1],"b",linewidth=1.5)

  # prepare x and y ranges

  xmin = -0.5
  xmax = 0.5
  ymin = 0
  ymax = 1.2

  # set axis parameters
  pl.axis([xmin,xmax,ymin,ymax])
  # get axis object
  ax = pl.gca()

  # set the custom tick sizes we like
  # these functions are defined in plot_defaults
  set_ticklines(ax,1.0,0.75*1.0)
  set_tick_sizes(ax, 13, 7)

  # label the axes
  pl.xlabel("R",labelpad=10)
  pl.ylabel("Density",labelpad=12)

  # legend
  pl.legend( [p1,p2], ["SPH"," Exact Riemann"], loc=(0.5,0.75), frameon=False, labelspacing = 0.25 )

  pl.savefig("compare" + suffix + ".pdf")
