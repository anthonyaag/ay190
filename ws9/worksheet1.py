#!/usr/bin/env python
import numpy as np
from scipy.linalg import lu
import time

#### Read in matrix and vectors 
def read_matrix(f_string):
    # Reads a file 
    f = open(f_string,'r')
    matrix = []
    for line in f.readlines():
        row = np.array(map(float,line.split()))
        matrix.append(row)
    f.close()
    if not all(len(matrix) == len(row) for row in matrix):
        print 'NOT A SQUARE MATRIX!!!'
    elif np.linalg.cond(matrix) > 1e10: # condition -> infinity if singular
        print 'NOT AN INVERTABLE MATRIX!!!'
    else:
        print"An intvertable square matrix of dimension {}.".format(len(matrix))
    return np.array(matrix)

def read_vector(f_string):
    # Reads a file
    f = open(f_string,'r')
    vector = []
    for line in f.readlines():
        row = np.array(float(line))
        vector.append(row)
    f.close()
    print "This is a vector of length {}.".format(len(vector))
    return np.array(vector)

def gauss(A, b):
    ''' Performs Gauss elimination '''
    n = len(b)

    # Check for zero diagonal elements
    if any(np.diag(A)==0):
        raise ZeroDivisionError(('Division by zero will occur; '
                                  'pivoting currently not supported'))

    # Float Everything so we don't get integer division
    A = np.array([map(float,i) for i in A])
    b = map(float,b)

    # Forward Elimination
    for row in range(0, n-1):
        for i in range(row+1, n):
            factor = A[i,row] / A[row,row]
            A[i,row:n] -= factor * A[row,row:n]
            b[i] -= factor * b[row]

    # Back Elmination
    ret = np.zeros((n,1))
    ret[n-1] = b[n-1] / A[n-1, n-1]
    for row in range(n-2, -1, -1):
        sums = b[row]
        sums -= np.dot(A[row,row+1:n],ret[row+1:n])
        ret[row] = sums / A[row,row]

    return ret

## Test on a matrix we know the answer to (5,3,-2)
test_m = np.array([[1, 1, 1],[0,2,5],[0,1,-8]])
test_b = np.array([6,-4,19])
print gauss(test_m, test_b)

## Read in the matricies
m1 = read_matrix('LSE1_m.dat')
m2 = read_matrix('LSE2_m.dat')
m3 = read_matrix('LSE3_m.dat')
m4 = read_matrix('LSE4_m.dat')
m5 = read_matrix('LSE5_m.dat')

b1 = read_vector('LSE1_bvec.dat')
b2 = read_vector('LSE2_bvec.dat')
b3 = read_vector('LSE3_bvec.dat')
b4 = read_vector('LSE4_bvec.dat')
b5 = read_vector('LSE5_bvec.dat')

### Time the various sizes using Gauss elimination
repeats = 5
tic = time.clock()
for i in range(repeats*100):
    gauss(m1,b1)
toc = time.clock()
time1 = toc - tic

tic = time.clock()
for i in range(repeats):
    gauss(m2,b2)
toc = time.clock()
time2 = toc - tic

tic = time.clock()
for i in range(repeats):
    gauss(m3,b3)
toc = time.clock()
time3 = toc - tic

tic = time.clock()
for i in range(repeats):
    gauss(m4,b4)
toc = time.clock()
time4 = toc - tic

tic = time.clock()
for i in range(repeats):
    gauss(m5,b5)
toc = time.clock()
time5 = toc - tic

print np.array([time1/100, time2, time3, time4, time5])/repeats
# (0.00164, 0.138, 0.56, 18.148, 88.07)

## Compare to numpy linear solvers
repeats = 100
tic = time.clock()
for i in range(10*repeats):
    np.linalg.solve(m1,b1)
toc = time.clock()
time1 = toc - tic

tic = time.clock()
for i in range(10*repeats):
    np.linalg.solve(m2,b2)
toc = time.clock()
time2 = toc - tic

tic = time.clock()
for i in range(repeats*10):
    np.linalg.solve(m3,b3)
toc = time.clock()
time3 = toc - tic

tic = time.clock()
for i in range(repeats):
    np.linalg.solve(m4,b4)
toc = time.clock()
time4 = toc - tic

tic = time.clock()
for i in range(repeats):
    np.linalg.solve(m5,b5)
toc = time.clock()
time5 = toc - tic

print np.array([time1/10, time2/10, time3/10, time4, time5])/repeats
# (9e-5, 1.22e-3, 5.82e-3, 5.007e-1, 3.63)
