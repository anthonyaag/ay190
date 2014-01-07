#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as pl

# read in the data
data = np.loadtxt("infile3.txt",comments="#")

# slice them
x = data[:,0]
y1 = data[:,1]
y2 = data[:,2]

# make the plot
p1, = pl.plot(x,y1,"r",linewidth=2)
p2, = pl.plot(x,y2,"b",linewidth=2)

# set x and y ranges
pl.xlim(min(x),max(x))
pl.ylim(min(y1*1.05),max(y1*1.05))

# label the axes
pl.xlabel("X")
pl.ylabel("Y")

# legend
# loc is the location of the legend in a
# coordinate system from (0,0) to (1,1)
# frameon=False turns of the stupid box around the legend
pl.legend( (p1,p2), ("y1","y2"), loc=(0.7,0.85), frameon=False )

# comment the below line to get rid of screen display
pl.show()
# uncomment the below line to save as a pdf
#pl.savefig("simpleplot.pdf")
