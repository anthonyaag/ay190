#!/usr/bin/env python

import numpy as np
import timeit

# define a function that defines
# and fills a list, then squares each
# list entry by iterating throught the list
def listtest():
    L = range(1000)
    [i**2 for i in L]

# define a function that defines
# a numpy array by filling it and that
# then squares it
def arraytest():
    a = np.arange(1000)
    b = a**2

ntrials = 100000
# use the timeit module to execute 
# and time the listtest function
time = timeit.timeit('listtest()',\
                     setup="from __main__ import listtest",\
                     number=ntrials)
print "Time for Python list: %5.3f microseconds" % (time/ntrials*1.0e6)

# use the timeit module to execute 
# and time the arraytest function
time = timeit.timeit('arraytest()',\
                     setup="from __main__ import arraytest",
                     number=ntrials)
print "Time for NumPy array: %5.3f microseconds" % (time/ntrials*1.0e6)
