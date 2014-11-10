
from numpy import *

def forwsubst(matrix, b, xinit):
    
    M       = max(shape(matrix))
    xout    = zeros(M)

    for i in xrange(0,M):
        xout[i] = ( b[i] - dot(matrix[i, 0:i], xinit[0:i]) ) / matrix[i, i]

    return xout


N = 1024

