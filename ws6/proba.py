#!/usr/bin/env python
import numpy as np
from numpy import arange, array
from numpy import abs, exp, sin,cos,sqrt
import matplotlib.pyplot as pl

# Define the imaginary number i

i = 1j

# Problem a
def dft(x):
    N = len(x)
    w = exp(-2 * np.pi * i / N)
    A = []
    for j in range(N):
        row = []
        for k in range(N):
            row.append(w ** (j * k))
        A.append(array(row))
    return np.dot(array(A),array(x))

test1 = np.random.randn(10)
test2 = np.random.randn(10)
test3 = np.random.randn(10)

#Check to see equivalence 
print abs(dft(test1) - np.fft.fft(test1)) < 0.00000001
print abs(dft(test2) - np.fft.fft(test2)) < 0.00000001
print abs(dft(test2) - np.fft.fft(test2)) < 0.00000001


    
