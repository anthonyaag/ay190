#!/usr/bin/env python

import numpy as np

# get a random number in [0.0,1.0)
a = np.random.random()

if a < 0.5:
    print "a = %5.6g is < 0.5!" % (a)
    print "This is fun!"
else:
    print "a = %5.6g is >= 0.5!" % (a)
    print "This is okay, I guess"

