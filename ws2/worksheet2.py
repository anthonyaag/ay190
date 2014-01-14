#!/user/bin/env python

from numpy import float32, arange

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

print recursive_thirds(15)

# 2

def f(x):
    return x ** 3 - 5 * x ** 2 + x

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









