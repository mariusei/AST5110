#!/usr/bin/env python

# AST5110
# Mon 29/09/2014
# Marius Berge Eide

import matplotlib.pylab as plt
import matplotlib.animation as animation
from numpy import *


# Exercise 2

# Operators

# Derivative

def deriv(farray, dx, dim=0, direction=1):
    """Derivative of array FARRAY along dimension DIM in 
    direction DIRECTION, where:
     1: Shift 1/2 step forward
    -1: shift 1/2 step backward
    and DX is step length"""

    c =  3./640
    b = -1./24          - 5.*c
    a =  1.     - 3.*b  + 5.*c

    coeffs = array([a,b,c])  / dx

    p = +direction

    start= array([farray, \
            roll(farray, -p,    axis=dim), \
            roll(farray, -2*p,  axis=dim)  ])

    stop = array([roll(farray, p,    axis=dim), \
            roll(farray, 2*p,  axis=dim), \
            roll(farray, 3*p,  axis=dim)  ])

    return dot(coeffs, (start - stop))

def interp(farray, dim=0, direction=1):
    """Interpolates the array FARRAY along dimension DIM in
    direction DIRECTION returning the values at shifted 1/2 * DIRECTION
    location"""

    c =  3./256
    b = -1./16          - 3.*c
    a =  1./2  - b      - c

    coeffs = array([a,b,c])

    p = +direction

    start= array([farray, \
            roll(farray, -p,    axis=dim), \
            roll(farray, -2*p,  axis=dim)  ])

    stop = array([roll(farray, p,    axis=dim), \
            roll(farray, 2*p,  axis=dim), \
            roll(farray, 3*p,  axis=dim)  ])

    return dot(coeffs, (start + stop))
